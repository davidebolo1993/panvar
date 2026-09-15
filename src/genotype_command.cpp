#include "panvar/genotype_command.hpp"

#include "panvar/bubbles.hpp"
#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_locus.hpp"
#include "panvar/genotype_model.hpp"
#include "panvar/gfa.hpp"
#include "panvar/gfa_io.hpp"
#include "panvar/graph_sort.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/gz_reader.hpp"
#include "panvar/integrated_snarls.hpp"
#include "panvar/ls_hmm.hpp"
#include "panvar/marker_model.hpp"
#include "panvar/mosaic_beam.hpp"
#include "panvar/mosaic_spelling.hpp"
#include "panvar/output.hpp"

#include <algorithm>
#include <chrono>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

namespace panvar {
namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point from) {
    return std::chrono::duration<double>(Clock::now() - from).count();
}

// Named stage timings, because "it took 90 seconds" is not a performance result: the go/no-go gate
// asks which stage spent them.
struct StageTimings {
    double model = 0.0;
    double markers = 0.0;
    double proposal = 0.0;
    double spelling = 0.0;
    double rescoring = 0.0;
    double total = 0.0;
};

std::vector<SequencingRead> read_fastq(const std::string& path) {
    GzLineReader reader(path);
    if (!reader.ok()) throw std::runtime_error("cannot open reads file: " + path);
    std::vector<SequencingRead> reads;
    std::string line;
    while (reader.getline(line)) {
        if (line.empty()) continue;
        if (line[0] != '@') {
            throw std::runtime_error("expected a FASTQ record header in " + path + ", found: " +
                                     line.substr(0, 40));
        }
        SequencingRead record;
        record.name = line.substr(1);
        const std::size_t space = record.name.find_first_of(" \t");
        if (space != std::string::npos) record.name.resize(space);
        if (!reader.getline(record.bases)) {
            throw std::runtime_error("truncated FASTQ record in " + path);
        }
        std::string plus;
        if (!reader.getline(plus) || plus.empty() || plus[0] != '+') {
            throw std::runtime_error("malformed FASTQ separator in " + path);
        }
        if (!reader.getline(record.qualities)) {
            throw std::runtime_error("truncated FASTQ qualities in " + path);
        }
        for (char& base : record.bases) base = static_cast<char>(std::toupper(base));
        reads.push_back(std::move(record));
    }
    return reads;
}

std::vector<ReadPair> read_fragments(const std::string& first_path, const std::string& second_path) {
    const std::vector<SequencingRead> first = read_fastq(first_path);
    const std::vector<SequencingRead> second = read_fastq(second_path);
    if (first.size() != second.size()) {
        throw std::runtime_error(
            "paired read files hold different record counts (" + std::to_string(first.size()) +
            " and " + std::to_string(second.size()) + "); version one requires paired reads in "
            "matching order");
    }
    std::vector<ReadPair> fragments;
    fragments.reserve(first.size());
    for (std::size_t i = 0; i < first.size(); ++i) {
        fragments.push_back({first[i], second[i]});
    }
    return fragments;
}

// Declared library parameters, one `key value` per line. Estimation from the data is not
// implemented in version one, so the declaration is required and its provenance is recorded.
LibraryModel read_library(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open library model file: " + path);
    LibraryModel model;
    std::map<std::string, double> values;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream fields(line);
        std::string key;
        double value = 0.0;
        if (!(fields >> key >> value)) {
            throw std::runtime_error("malformed library model line: " + line);
        }
        values[key] = value;
    }
    const auto take = [&](const char* key) {
        const auto found = values.find(key);
        if (found == values.end()) {
            throw std::runtime_error(std::string("library model is missing '") + key + "'");
        }
        return found->second;
    };
    model.insert_mean_bp = take("insert_mean_bp");
    model.insert_sd_bp = take("insert_sd_bp");
    model.support_sigmas = take("support_sigmas");
    model.substitution_error = take("substitution_error");
    model.outlier_weight = take("outlier_weight");
    model.fragment_start_intensity = take("fragment_start_intensity");
    model.provenance = "declared: " + path;
    validate_library_model(model);
    return model;
}

std::ofstream open_output(const std::string& path) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write " + path);
    return out;
}

std::string join_alleles(const std::vector<std::size_t>& alleles) {
    std::string out;
    for (std::size_t i = 0; i < alleles.size(); ++i) {
        if (i != 0) out += ',';
        out += std::to_string(alleles[i]);
    }
    return out;
}

struct RankingSummary {
    std::size_t sequence_candidates = 0;
    std::size_t collapsed_beam_candidates = 0;
    std::string decided_by;
    double fragment_term_gap = 0.0;
    bool projection_ambiguous = false;
};

void write_status(
    const std::string& path,
    GenotypeRunStatus status,
    const std::string& fingerprint,
    const GenotypeCommandOptions& options,
    const LibraryModel& library,
    const StageTimings& timings,
    const LocusInputReport* model_report,
    const PreparedMarkerPanel* panel,
    const MosaicBeamResult* beam,
    const CertifiedRanking* ranking,
    const std::string& note,
    const RankingSummary* summary = nullptr) {
    std::ofstream out = open_output(path);
    out << "key\tvalue\n";
    out << "status\t" << genotype_run_status_name(status) << '\n';
    if (!note.empty()) out << "note\t" << note << '\n';
    out << "model_fingerprint\t" << fingerprint << '\n';
    out << "gfa\t" << options.gfa_path << '\n';
    out << "reference_path\t" << options.reference_path << '\n';
    out << "kmer_size\t" << options.kmer_size << '\n';
    out << "syncmer_s\t" << options.syncmer_s << '\n';
    out << "min_variant_bp\t" << options.min_variant_bp << '\n';
    out << "switch_probability\t" << options.switch_probability << '\n';
    out << "max_candidates\t" << options.max_candidates << '\n';
    out << "mismatch_band\t" << options.mismatch_band << '\n';
    // Version one is serial by contract: the scorer is not threaded before the serial algorithm
    // meets its performance gate. The requested value is recorded so a run is reproducible.
    out << "threads_requested\t" << options.threads << '\n';
    out << "threads_used\t1\n";
    out << "library_provenance\t" << library.provenance << '\n';
    out << "insert_mean_bp\t" << library.insert_mean_bp << '\n';
    out << "insert_sd_bp\t" << library.insert_sd_bp << '\n';
    out << "substitution_error\t" << library.substitution_error << '\n';
    out << "outlier_weight\t" << library.outlier_weight << '\n';
    out << "fragment_start_intensity\t" << library.fragment_start_intensity << '\n';
    if (model_report != nullptr) {
        out << "variant_blocks\t" << model_report->variant_blocks << '\n';
        out << "context_blocks\t" << model_report->context_blocks << '\n';
        out << "source_alleles\t" << model_report->total_source_alleles << '\n';
        out << "canonical_alleles\t" << model_report->total_canonical_alleles << '\n';
        out << "templates_forward\t" << model_report->forward_templates << '\n';
        out << "templates_reverse\t" << model_report->reverse_templates << '\n';
        out << "templates_refused\t" << model_report->refused_templates.size() << '\n';
        out << "boundary_bases_moved_to_context\t"
            << model_report->boundary_bases_moved_to_context << '\n';
    }
    if (panel != nullptr) {
        out << "marker_universe\t" << panel->codes.size() << '\n';
        out << "marker_anchors\t" << panel->anchor_slots.size() << '\n';
        out << "markers_dropped_multiblock\t" << panel->dropped_multiblock_markers << '\n';
        out << "alleles_scored\t" << panel->scored_alleles << '\n';
        out << "alleles_total\t" << panel->total_alleles << '\n';
        out << "every_allele_scored\t"
            << (panel->scored_alleles == panel->total_alleles ? "true" : "false") << '\n';
    }
    if (beam != nullptr) {
        out << "beam_width_used\t" << beam->beam_width_used << '\n';
        out << "raw_histories\t" << beam->final_raw_histories << '\n';
        out << "unique_candidates\t" << beam->candidates.size() << '\n';
        out << "beam_extensions\t" << beam->extensions_evaluated << '\n';
        out << "unique_target_reached\t" << (beam->unique_target_reached ? "true" : "false") << '\n';
    }
    if (ranking != nullptr) {
        out << "certified_gap\t" << ranking->certified_gap << '\n';
        out << "winner_set_size\t" << ranking->winner_set.size() << '\n';
        out << "score_gap_to_runner_up\t" << ranking->score_gap_to_runner_up << '\n';
        out << "candidate_set_log_weight\t" << ranking->candidate_set_log_weight << '\n';
        out << "winner_candidate_set_probability\t"
            << ranking->candidates.front().candidate_set_probability << '\n';
        out << "winner_numeric_error_bound\t" << ranking->candidates.front().numeric_error_bound
            << '\n';
        if (summary != nullptr) {
            out << "sequence_candidates\t" << summary->sequence_candidates << '\n';
            out << "beam_candidates_collapsed_by_sequence\t" << summary->collapsed_beam_candidates
                << '\n';
            // A winner separated only by the prior is a winner the reads did not choose. At C4 the
            // runner-up is often the same two haplotypes phased differently across a block the
            // inserts cannot span, and there the whole gap is the switch penalty.
            out << "winner_decided_by\t" << summary->decided_by << '\n';
            out << "fragment_term_gap_to_runner_up\t" << summary->fragment_term_gap << '\n';
            out << "winner_projection_ambiguous\t"
                << (summary->projection_ambiguous ? "true" : "false") << '\n';
        }
        const FragmentScoringCounters& counters = ranking->counters;
        out << "candidate_bases_scanned\t" << counters.candidate_bases_scanned << '\n';
        out << "piece_hits\t" << counters.piece_hits << '\n';
        out << "proposed_starts\t" << counters.proposed_starts << '\n';
        out << "unique_starts\t" << counters.unique_starts << '\n';
        out << "full_read_verifications\t" << counters.full_read_verifications << '\n';
        out << "accepted_mate_placements\t" << counters.accepted_mate_placements << '\n';
        out << "valid_fr_states\t" << counters.valid_fr_states << '\n';
        out << "physical_origins\t" << counters.physical_origins << '\n';
        out << "statistical_state_classes\t" << counters.statistical_state_classes << '\n';
    }
    out << std::fixed << std::setprecision(3);
    out << "seconds_model\t" << timings.model << '\n';
    out << "seconds_markers\t" << timings.markers << '\n';
    out << "seconds_proposal\t" << timings.proposal << '\n';
    out << "seconds_spelling\t" << timings.spelling << '\n';
    out << "seconds_rescoring\t" << timings.rescoring << '\n';
    out << "seconds_total\t" << timings.total << '\n';
}

} // namespace

const char* genotype_run_status_name(GenotypeRunStatus status) {
    switch (status) {
        case GenotypeRunStatus::Complete: return "COMPLETE";
        case GenotypeRunStatus::ModelIncomplete: return "MODEL_INCOMPLETE";
        case GenotypeRunStatus::CandidateIncomplete: return "CANDIDATE_INCOMPLETE";
        case GenotypeRunStatus::SpellingIncomplete: return "SPELLING_INCOMPLETE";
        case GenotypeRunStatus::RescoreIncomplete: return "RESCORE_INCOMPLETE";
    }
    return "MODEL_INCOMPLETE";
}

GenotypeRunStatus run_genotype(const GenotypeCommandOptions& options) {
    const Clock::time_point started = Clock::now();
    StageTimings timings;
    const std::string status_path = options.out_prefix + ".mosaic_status.tsv";

    if (options.gfa_path.empty() || options.reads1_path.empty() || options.reads2_path.empty() ||
        options.reference_path.empty() || options.library_path.empty()) {
        throw std::runtime_error(
            "genotype requires --gfa, --reads1, --reads2, --reference-path and --library");
    }

    const LibraryModel library = read_library(options.library_path);

    // ---- stage 1: structural model -------------------------------------------------------------
    Clock::time_point mark = Clock::now();
    GfaModel gfa = read_gfa_model(options.gfa_path);
    GraphSortOptions sort_options;
    sort_options.reference_path = options.reference_path;
    sort_options.flip = true;
    const std::string resolved = sort_graph_reference(gfa, sort_options).resolved_reference;
    const std::string staged = options.out_prefix + ".sorted.gfa";
    write_gfa_model(staged, gfa);
    ParseGfaOptions parse_options;
    const Graph graph = parse_gfa(staged, parse_options);
    validate_graph_paths(graph, "genotype", true, true);

    std::vector<Bubble> bubbles;
    if (options.bubbles_path.empty()) {
        BubbleCallOptions bubble_options;
        bubble_options.reference_path = resolved;
        bubble_options.quiet = true;
        bubble_options.min_variant_bp = options.min_variant_bp;
        bubble_options.snarl_pairs_override =
            find_top_level_snarls_cactus(snarl_input_from_model(gfa));
        bubble_options.snarl_source_supplied = true;
        bubbles = call_bubbles(graph, bubble_options);
    } else {
        bubbles = read_bubbles_csv(options.bubbles_path);
    }

    LocusInputOptions locus_options;
    // The locus identity is the graph's own reference frame, never the output path: a fingerprint
    // that moved when the user renamed their output would be useless for saying whether two runs
    // scored the same model.
    locus_options.locus_id = resolved;
    locus_options.reference_path_name = resolved;
    const LocusInputResult built = build_locus_input(graph, bubbles, locus_options);
    const PreparedGenotypeModel model = prepare_genotype_model(built.input);
    // The authoritative check, on real input, before any statistic is computed from these blocks.
    verify_complete_panel_round_trip(model, built.authoritative);
    const std::string fingerprint = model_fingerprint(model);
    timings.model = seconds_since(mark);

    if (model.blocks().size() < 2 || model.panel().template_count() < 2) {
        timings.total = seconds_since(started);
        write_status(status_path, GenotypeRunStatus::ModelIncomplete, fingerprint, options, library,
                     timings, &built.report, nullptr, nullptr, nullptr,
                     "the locus does not carry at least two blocks and two panel templates");
        return GenotypeRunStatus::ModelIncomplete;
    }

    // ---- stage 2: markers ----------------------------------------------------------------------
    mark = Clock::now();
    const std::vector<ReadPair> fragments =
        read_fragments(options.reads1_path, options.reads2_path);
    const PreparedMarkerPanel panel =
        build_marker_panel(model, options.kmer_size, options.syncmer_s);
    const MarkerObservations observations = count_marker_observations(panel, fragments);
    MarkerModelParameters marker_parameters;
    marker_parameters.kmer_size = options.kmer_size;
    marker_parameters.syncmer_s = options.syncmer_s;
    // Expected marker observations per haplotype copy, from the reads actually seen rather than a
    // declared depth: matched occurrences spread over the marker universe, halved for two copies.
    marker_parameters.expected_observations_per_copy =
        panel.codes.empty()
            ? 1.0
            : std::max(1e-6, static_cast<double>(observations.matched_occurrences) /
                                 (2.0 * static_cast<double>(panel.codes.size())));
    marker_parameters.background_observations = 0.5;
    marker_parameters.provenance = "estimated from matched marker occurrences";
    const PreparedProposalEmissions emissions =
        build_proposal_emissions(panel, observations, marker_parameters);
    timings.markers = seconds_since(mark);

    // ---- stage 3: proposal ---------------------------------------------------------------------
    mark = Clock::now();
    const std::vector<double> switch_probability(
        model.blocks().size() - 1, options.switch_probability);
    MosaicBeamParameters beam_parameters;
    beam_parameters.unique_candidate_target = options.max_candidates;
    beam_parameters.raw_history_limit = std::max<std::size_t>(options.max_candidates * 4, 64);
    const MosaicBeamResult beam = propose_mosaic_diplotypes(
        model.panel(), emissions, switch_probability, beam_parameters);
    timings.proposal = seconds_since(mark);

    if (beam.candidates.empty()) {
        timings.total = seconds_since(started);
        write_status(status_path, GenotypeRunStatus::CandidateIncomplete, fingerprint, options,
                     library, timings, &built.report, &panel, &beam, nullptr,
                     "the proposal produced no mosaic candidate");
        return GenotypeRunStatus::CandidateIncomplete;
    }

    // ---- stage 4: spelling and interning --------------------------------------------------------
    mark = Clock::now();
    std::vector<MosaicPath> paths;
    std::vector<std::string> sequences;
    std::vector<HaplotypeCandidate> pairs;
    std::map<std::vector<std::size_t>, std::size_t> slot_of_path;
    const auto slot_for = [&](const MosaicPath& path) {
        const auto found = slot_of_path.find(path.alleles);
        if (found != slot_of_path.end()) return found->second;
        const std::size_t slot = paths.size();
        slot_of_path.emplace(path.alleles, slot);
        paths.push_back(path);
        sequences.push_back(spell_mosaic(model, path).bases);
        return slot;
    };
    for (const MosaicCandidate& candidate : beam.candidates) {
        pairs.push_back({slot_for(candidate.diplotype.first),
                         slot_for(candidate.diplotype.second)});
    }
    const InternedHaplotypes interned = intern_haplotypes(paths, sequences);

    // Candidates are ranked as unordered pairs of SEQUENCES. Two beam candidates whose allele paths
    // differ but spell the same two haplotypes are one state to any read-based method; ranking them
    // separately would leave two candidates with identical fragment evidence and a gap made entirely
    // of how differently their paths happen to be represented.
    struct SequenceCandidate {
        std::size_t first = 0;     // interned sequence, first <= second
        std::size_t second = 0;
        std::vector<std::size_t> beam_candidates;
        std::size_t representative = 0;
        bool projection_ambiguous = false;
    };
    std::vector<SequenceCandidate> sequence_candidates;
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> slot_of_sequence_pair;
    std::vector<std::size_t> sequence_candidate_of_beam(beam.candidates.size(), 0);
    for (std::size_t i = 0; i < beam.candidates.size(); ++i) {
        std::size_t x = interned.index_of[pairs[i].first];
        std::size_t y = interned.index_of[pairs[i].second];
        if (y < x) std::swap(x, y);
        const auto inserted =
            slot_of_sequence_pair.emplace(std::make_pair(x, y), sequence_candidates.size());
        if (inserted.second) {
            SequenceCandidate fresh;
            fresh.first = x;
            fresh.second = y;
            sequence_candidates.push_back(fresh);
        }
        sequence_candidates[inserted.first->second].beam_candidates.push_back(i);
        sequence_candidate_of_beam[i] = inserted.first->second;
    }
    timings.spelling = seconds_since(mark);

    // ---- stage 5: exact prior, recomputed independently of the beam ------------------------------
    // Both priors are exact and neither reads a beam history. The path prior is kept as a manifest
    // diagnostic; the ranking uses the sequence-pair prior, which sums every allele path that spells
    // each haplotype and decides the log(2) multiplicity on bytes.
    std::vector<double> path_prior;
    path_prior.reserve(beam.candidates.size());
    for (const MosaicCandidate& candidate : beam.candidates) {
        path_prior.push_back(
            log_unordered_diplotype_prior(model.panel(), candidate.diplotype, switch_probability));
    }
    std::vector<double> haploid_sequence_prior(interned.sequences.size());
    for (std::size_t h = 0; h < interned.sequences.size(); ++h) {
        haploid_sequence_prior[h] =
            log_haploid_sequence_prior(model, interned.sequences[h], switch_probability);
    }
    std::vector<double> sequence_prior;
    std::vector<HaplotypeCandidate> ranked_pairs;
    for (SequenceCandidate& candidate : sequence_candidates) {
        sequence_prior.push_back(haploid_sequence_prior[candidate.first] +
                                 haploid_sequence_prior[candidate.second] +
                                 (candidate.first != candidate.second ? std::log(2.0) : 0.0));
        ranked_pairs.push_back({candidate.first, candidate.second});

        // Projection needs one allele path per haplotype. Take the member with the highest exact path
        // prior, ties to beam order, and flag the call if any member projects to different block
        // alleles: those bytes then do not determine the block genotype.
        std::size_t best = candidate.beam_candidates.front();
        for (const std::size_t member : candidate.beam_candidates) {
            if (path_prior[member] > path_prior[best]) best = member;
        }
        candidate.representative = best;
        const MosaicDiplotype& chosen = beam.candidates[best].diplotype;
        for (const std::size_t member : candidate.beam_candidates) {
            const MosaicDiplotype& other = beam.candidates[member].diplotype;
            for (std::size_t b = 0; b < model.blocks().size(); ++b) {
                auto lhs = std::minmax(chosen.first.alleles[b], chosen.second.alleles[b]);
                auto rhs = std::minmax(other.first.alleles[b], other.second.alleles[b]);
                if (lhs != rhs) candidate.projection_ambiguous = true;
            }
        }
    }

    // ---- stage 6: certified rescoring -----------------------------------------------------------
    mark = Clock::now();
    FragmentScoringOptions scoring_options;
    scoring_options.mismatch_band = options.mismatch_band;
    scoring_options.minimum_piece_bp = options.minimum_piece_bp;
    const CertifiedRanking ranking = certify_diplotype_ranking(
        interned, ranked_pairs, fragments, library, scoring_options, CertificationOptions{},
        sequence_prior);
    timings.rescoring = seconds_since(mark);

    RankingSummary summary;
    summary.sequence_candidates = sequence_candidates.size();
    summary.collapsed_beam_candidates = beam.candidates.size() - sequence_candidates.size();
    {
        const CertifiedCandidate& best = ranking.candidates.front();
        summary.decided_by = "fragments";
        summary.fragment_term_gap = std::numeric_limits<double>::infinity();
        for (std::size_t i = 1; i < ranking.candidates.size(); ++i) {
            const CertifiedCandidate& rival = ranking.candidates[i];
            const double margin = (best.fragment_term - best.numeric_error_bound) -
                                  (rival.fragment_term + rival.numeric_error_bound);
            if (i == 1) summary.fragment_term_gap = best.fragment_term - rival.fragment_term;
            if (margin <= 0.0) summary.decided_by = "prior";
        }
        if (ranking.candidates.size() == 1) summary.decided_by = "single_candidate";
        summary.projection_ambiguous = sequence_candidates[best.candidate_id].projection_ambiguous;
    }

    // ---- outputs --------------------------------------------------------------------------------
    {
        std::ofstream out = open_output(options.out_prefix + ".mosaic_candidates.tsv");
        out << "candidate_id\tunique_rank\thap1_alleles\thap2_alleles\tbest_proposal_score\t"
               "retained_logsum_proposal_score\tproposal_emission_component\t"
               "proposal_transition_component\tequivalent_retained_histories\tevery_allele_scored\t"
               "exact_path_prior\tsequence_candidate_id\n";
        for (std::size_t i = 0; i < beam.candidates.size(); ++i) {
            const MosaicCandidate& candidate = beam.candidates[i];
            out << candidate.candidate_id << '\t' << candidate.unique_rank << '\t'
                << join_alleles(candidate.diplotype.first.alleles) << '\t'
                << join_alleles(candidate.diplotype.second.alleles) << '\t'
                << candidate.best_proposal_score << '\t'
                << candidate.retained_logsum_proposal_score << '\t'
                << candidate.proposal_emission_component << '\t'
                << candidate.proposal_transition_component << '\t'
                << candidate.equivalent_retained_histories << '\t'
                << (candidate.every_allele_scored ? "true" : "false") << '\t'
                << path_prior[i] << '\t' << sequence_candidate_of_beam[i] << '\n';
        }
    }
    {
        // The proposal score is deliberately absent from this table's total: it comes from the same
        // reads as the fragment term, so adding both would count that evidence twice.
        std::ofstream out = open_output(options.out_prefix + ".mosaic_scores.tsv");
        out << "sequence_candidate_id\trepresentative_beam_candidate\tequivalent_beam_candidates\t"
               "score\tfragment_term\texposure_term\tsequence_pair_prior\tnumeric_error_bound\t"
               "fragments_with_mass\tfragments_total\tcandidate_set_probability\tmismatch_band\t"
               "projection_ambiguous\tin_winner_set\n";
        out << std::setprecision(12);
        for (const CertifiedCandidate& entry : ranking.candidates) {
            const SequenceCandidate& candidate = sequence_candidates[entry.candidate_id];
            const bool in_set =
                std::find(ranking.winner_set.begin(), ranking.winner_set.end(),
                          entry.candidate_id) != ranking.winner_set.end();
            out << entry.candidate_id << '\t' << candidate.representative << '\t'
                << candidate.beam_candidates.size() << '\t' << entry.score << '\t'
                << entry.fragment_term << '\t' << entry.exposure_term << '\t' << entry.prior << '\t'
                << entry.numeric_error_bound << '\t' << entry.fragments_with_mass << '\t'
                << fragments.size() << '\t' << entry.candidate_set_probability << '\t'
                << ranking.mismatch_band << '\t'
                << (candidate.projection_ambiguous ? "true" : "false") << '\t'
                << (in_set ? "true" : "false") << '\n';
        }
    }

    if (ranking.status != RescoreStatus::Complete) {
        timings.total = seconds_since(started);
        write_status(status_path, GenotypeRunStatus::RescoreIncomplete, fingerprint, options,
                     library, timings, &built.report, &panel, &beam, &ranking,
                     "the best candidate cannot be separated from every rival by more than "
                     "floating-point error; the equivalence set is reported instead of a call",
                     &summary);
        return GenotypeRunStatus::RescoreIncomplete;
    }

    const std::size_t winner = ranking.winner_set.front();
    const MosaicDiplotype& called =
        beam.candidates[sequence_candidates[winner].representative].diplotype;
    {
        // Only variant blocks are reporting units. Context blocks are modelled identically and are
        // what let the mosaic recombine, but they are not SV calls and are not presented as such.
        std::ofstream out = open_output(options.out_prefix + ".genotypes.tsv");
        out << "block_id\tkind\tallele1\tallele2\tallele1_bp\tallele2_bp\tallele1_sources\t"
               "allele2_sources\tzygosity\n";
        for (std::size_t b = 0; b < model.blocks().size(); ++b) {
            const PreparedGenotypeBlock& block = model.blocks()[b];
            if (block.kind != GenotypeBlockKind::Variant) continue;
            const std::size_t a1 = called.first.alleles[b];
            const std::size_t a2 = called.second.alleles[b];
            const auto sources = [&](std::size_t id) {
                std::string joined;
                for (const std::string& source : block.alleles[id].source_ids) {
                    if (!joined.empty()) joined += ',';
                    joined += source;
                }
                return joined;
            };
            out << block.block_id << "\tSV\t" << a1 << '\t' << a2 << '\t'
                << block.alleles[a1].sequence.size() << '\t'
                << block.alleles[a2].sequence.size() << '\t'
                << sources(a1) << '\t' << sources(a2) << '\t'
                << (a1 == a2 ? "HOM" : "HET") << '\n';
        }
    }
    {
        std::ofstream out = open_output(options.out_prefix + ".mosaic_haplotypes.fa");
        const HaplotypeSequence first = spell_mosaic(model, called.first);
        const HaplotypeSequence second = spell_mosaic(model, called.second);
        for (const auto& item : {std::make_pair("hap1", &first), std::make_pair("hap2", &second)}) {
            out << '>' << options.out_prefix << '_' << item.first
                << " fingerprint=" << fingerprint << " bp=" << item.second->bases.size() << '\n';
            for (std::size_t at = 0; at < item.second->bases.size(); at += 60) {
                out << item.second->bases.substr(at, 60) << '\n';
            }
        }
    }

    timings.total = seconds_since(started);
    write_status(status_path, GenotypeRunStatus::Complete, fingerprint, options, library, timings,
                 &built.report, &panel, &beam, &ranking, "", &summary);
    if (!options.quiet) {
        std::cerr << "[genotype] COMPLETE  fingerprint=" << fingerprint
                  << "  candidates=" << beam.candidates.size()
                  << "  gap=" << ranking.certified_gap
                  << "  total=" << timings.total << "s\n";
    }
    return GenotypeRunStatus::Complete;
}

int run_genotype_command(const std::vector<std::string>& args) {
    GenotypeCommandOptions options;
    const auto need = [&](std::size_t& i) {
        if (i + 1 >= args.size()) throw std::runtime_error("missing value for " + args[i]);
        return args[++i];
    };
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& flag = args[i];
        if (flag == "--gfa" || flag == "-i") options.gfa_path = need(i);
        else if (flag == "--bubbles") options.bubbles_path = need(i);
        else if (flag == "--reads1") options.reads1_path = need(i);
        else if (flag == "--reads2") options.reads2_path = need(i);
        else if (flag == "--reference-path" || flag == "-r") options.reference_path = need(i);
        else if (flag == "--out-prefix" || flag == "-o") options.out_prefix = need(i);
        else if (flag == "--library") options.library_path = need(i);
        else if (flag == "--threads") options.threads = std::stoul(need(i));
        else if (flag == "--max-candidates") options.max_candidates = std::stoul(need(i));
        else if (flag == "-q" || flag == "--quiet") options.quiet = true;
        else if (flag == "-h" || flag == "--help") {
            std::cout <<
                "Usage:\n"
                "  panvar genotype --gfa <graph.gfa> --reference-path <name> \\\n"
                "                  --reads1 <r1.fq> --reads2 <r2.fq> --library <params.txt> \\\n"
                "                  [--bubbles <bubbles.csv>] [--out-prefix <prefix>] \\\n"
                "                  [--threads <N>] [--max-candidates <K>]\n\n"
                "Calls a pair of whole-locus mosaic haplotypes and projects them to one allele\n"
                "pair per SV block. Paired reads only.\n\n"
                "Outputs <prefix>.genotypes.tsv, .mosaic_candidates.tsv, .mosaic_scores.tsv,\n"
                ".mosaic_status.tsv and .mosaic_haplotypes.fa. An incomplete run emits its\n"
                "diagnostics and says so in the status file; it never presents a partial ranking\n"
                "as a result.\n\n"
                "--library declares insert_mean_bp, insert_sd_bp, support_sigmas,\n"
                "substitution_error, outlier_weight and fragment_start_intensity, one per line.\n"
                "Version one does not estimate these from the data.\n";
            return 0;
        } else {
            throw std::runtime_error("unknown genotype option: " + flag);
        }
    }
    const GenotypeRunStatus status = run_genotype(options);
    return status == GenotypeRunStatus::Complete ? 0 : 1;
}

} // namespace panvar
