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
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
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
LibraryModel read_library(const std::string& path, std::optional<double>& marker_background) {
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
    // Optional. When absent the marker background is estimated jointly with the per-copy rate from
    // the invariant anchors, which is possible whenever those anchors carry two distinct copy numbers.
    const auto background = values.find("marker_background");
    if (background != values.end()) marker_background = background->second;
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

// The protected panel baseline and what the beam added beyond it.
struct CandidateSetReport {
    std::size_t panel_diplotypes = 0;
    bool panel_baseline_complete = false;
    std::size_t beam_diplotypes_found = 0;   // every diplotype the beam found, panel pairs included
    // Distinct haplotype sequences in the candidate set this status describes -- the panel's, plus at most two
    // per mosaic -- which is what a preflight predicts. It is 0 when no set was ranked. The shared scoring pass
    // may cover more (a prefix scores the union), so this counts the set, not the pass.
    std::size_t unique_haplotypes = 0;
    // True when no beam ran at all: K = 0, no phase audit and no positive development prefix, so the
    // candidate set is exactly the panel and nothing about it is conditional on a search. Reported as
    // candidate_scope PANEL_ONLY, which is neither COMPLETE nor CONDITIONAL -- both of those describe
    // the outcome of a search that did not happen.
    bool panel_only = false;
};

struct RankingSummary {
    std::size_t sequence_candidates = 0;
    std::size_t collapsed_beam_candidates = 0;
    // READS_WITHIN_CANDIDATES: no retained rival is a fragment contender -- a claim about this candidate set
    // only. READS_AUDITED: that, and the declared larger audit set agrees on winner, block genotypes and phase.
    // PRIOR_ONLY: some rival's fragment term is not
    // below the winner's beyond numerical error, so the prior chose among them. UNRESOLVED: the ranking
    // could not separate a winner at all.
    std::string phase_support;
    std::string phase_support_within_candidates;
    std::size_t phase_audit_candidates = 0;
    std::string phase_audit_status = "NA";
    std::string phase_audit_same_winner = "NA";
    std::string phase_audit_same_block_genotypes = "NA";
    std::string phase_audit_phase_support = "NA";
    std::size_t fragment_contenders = 0;
    double fragment_term_gap = 0.0;
    bool projection_ambiguous = false;
    // Per block: every fragment contender projects to the winner's unordered allele pair there.
    std::vector<bool> block_genotype_read_supported;
    // Per block, within the retained candidate set only: the posterior of the winner's unordered genotype, its
    // quality -10 log10(1 - posterior), and whether every candidate with another genotype there sits below the
    // winner beyond floating-point error. Conditional, never global, while the set omits candidates.
    std::vector<double> block_genotype_posterior;
    std::vector<double> block_genotype_quality;
    std::vector<bool> block_numerically_separated;
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
    const RankingSummary* summary = nullptr,
    const MarkerDepthEstimate* depth = nullptr,
    const CandidateSetReport* candidate_set = nullptr) {
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
    out << "max_mosaic_candidates\t" << options.max_mosaic_candidates << '\n';
    out << "mismatch_band\t" << options.mismatch_band << '\n';
    // Version one is serial by contract: the scorer is not threaded before the serial algorithm
    // meets its performance gate. The requested value is recorded so a run is reproducible.
    // Version one is serial, and a request for more is refused before any work; both are recorded.
    out << "threads_requested\t" << options.threads << '\n';
    out << "threads_used\t1\n";
    out << "max_beam_extensions\t" << options.max_beam_extensions << '\n';
    out << "max_piece_hits\t" << options.max_piece_hits << '\n';
    out << "max_proposed_starts_per_haplotype\t" << options.max_proposed_starts_per_haplotype << '\n';
    out << "max_physical_origins\t" << options.max_physical_origins << '\n';
    out << "max_state_classes\t" << options.max_state_classes << '\n';
    out << "max_mass_cells\t" << options.max_mass_cells << '\n';
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
        std::string dependent;
        for (const std::string& block : model_report->representation_dependent_blocks) {
            dependent += (dependent.empty() ? "" : ",") + block;
        }
        out << "representation_dependent_blocks\t" << (dependent.empty() ? "none" : dependent) << '\n';
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
    if (depth != nullptr) {
        out << "marker_depth_rate_per_copy\t" << depth->rate_per_copy << '\n';
        out << "marker_background\t" << std::setprecision(17) << depth->background << '\n';
        out << "marker_background_floor\t" << kMarkerBackgroundFloor << '\n';
        // proposal_marker_parameters passes the depth background to emission unchanged (gated).
        out << "marker_background_effective\t" << depth->background << std::setprecision(6) << '\n';
        out << "marker_background_estimated\t" << (depth->background_estimated ? "true" : "false")
            << '\n';
        out << "marker_depth_converged\t" << (depth->converged ? "true" : "false") << '\n';
        out << "marker_background_at_floor\t" << (depth->background_at_floor ? "true" : "false")
            << '\n';
        out << "marker_depth_iterations\t" << depth->iterations << '\n';
        out << "marker_depth_termination\t" << depth->termination << '\n';
        out << "marker_depth_gradient_background\t" << depth->gradient_background << '\n';
        out << "marker_depth_gradient_rate\t" << depth->gradient_rate << '\n';
        out << "marker_depth_log_likelihood\t" << std::setprecision(15) << depth->log_likelihood
            << std::setprecision(6) << '\n';
        out << "marker_depth_likelihood_resolution\t" << depth->likelihood_resolution << '\n';
        out << "marker_depth_newton_decrement\t" << depth->newton_decrement << '\n';
        out << "marker_depth_anchors_used\t" << depth->anchors_used << '\n';
        out << "marker_depth_distinct_multiplicities\t" << depth->distinct_multiplicities << '\n';
        out << "marker_depth_provenance\t" << depth->provenance << '\n';
    }
    if (beam != nullptr) {
        out << "proposal_state_limit\t" << beam->state_limit_used << '\n';
        out << "proposal_final_states\t" << beam->final_states << '\n';
        out << "mosaic_candidates\t" << beam->candidates.size() << '\n';
        out << "mosaic_candidates_found\t" << beam->unique_candidates_found << '\n';
        out << "proposal_states_pruned\t" << beam->states_pruned << '\n';
        out << "proposal_search_exhaustive\t" << (beam->search_exhaustive ? "true" : "false") << '\n';
        out << "mosaic_target_reached\t" << (beam->unique_target_reached ? "true" : "false") << '\n';
        // Any ranking over a set that omits candidates is conditional on that set.
        out << "candidate_set_conditional\t" << (beam->candidate_set_conditional ? "true" : "false")
            << '\n';
        out << "candidate_scope\t" << (beam->candidate_set_conditional ? "CONDITIONAL" : "COMPLETE") << '\n';
        out << "beam_extensions\t" << beam->extensions_evaluated << '\n';
    }
    if (candidate_set != nullptr) {
        // mosaic_candidates, mosaic_candidates_found and mosaic_target_reached above count off-panel mosaics;
        // the panel baseline is added to them whole.
        out << "panel_diplotypes\t" << candidate_set->panel_diplotypes << '\n';
        out << "panel_baseline_complete\t" << (candidate_set->panel_baseline_complete ? "true" : "false") << '\n';
        out << "beam_diplotypes_found\t" << candidate_set->beam_diplotypes_found << '\n';
        out << "unique_haplotypes\t" << candidate_set->unique_haplotypes << '\n';
        if (candidate_set->panel_only) {
            // No search ran, so neither COMPLETE nor CONDITIONAL applies: both describe what a beam
            // returned. The candidate set is the panel, entire.
            out << "candidate_scope\tPANEL_ONLY\n";
        }
        if (beam != nullptr && candidate_set->panel_baseline_complete) {
            out << "candidates_total\t" << candidate_set->panel_diplotypes + beam->candidates.size() << '\n';
        }
    }
    if (ranking != nullptr) {
        // Separation of the winner from every retained rival beyond floating-point error: numerical, within the
        // retained set. It is not a statistical confidence; the block posteriors and GQ are.
        out << "numeric_separation_gap\t" << ranking->certified_gap << '\n';
        // The gap is certified within the retained candidate set only. It is global only when the search
        // was exhaustive and nothing was truncated from its result.
        out << "global_certified\t"
            << (ranking->status == RescoreStatus::Complete && beam != nullptr && !beam->candidate_set_conditional
                    ? "true" : "false") << '\n';
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
            // Phase the reads did not choose is not reported as read-supported. Unordered block
            // genotypes can still be, where every fragment contender agrees on them.
            out << "phase_support\t" << summary->phase_support << '\n';
            out << "phase_support_within_candidates\t" << summary->phase_support_within_candidates << '\n';
            out << "phase_audit_candidates\t"
                << (summary->phase_audit_candidates > 0 ? std::to_string(summary->phase_audit_candidates) : "NA")
                << '\n';
            out << "phase_audit_status\t" << summary->phase_audit_status << '\n';
            out << "phase_audit_same_winner\t" << summary->phase_audit_same_winner << '\n';
            out << "phase_audit_same_block_genotypes\t" << summary->phase_audit_same_block_genotypes << '\n';
            out << "phase_audit_phase_support\t" << summary->phase_audit_phase_support << '\n';
            out << "fragment_contenders\t" << summary->fragment_contenders << '\n';
            std::size_t supported = 0;
            for (const bool block : summary->block_genotype_read_supported) supported += block ? 1 : 0;
            out << "blocks_genotype_read_supported\t" << supported << '\n';
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

std::string audited_phase_support(const std::string& within_candidates, bool audit_complete, bool same_winner,
                                  bool same_block_genotypes, const std::string& audit_phase_support) {
    if (within_candidates != "READS_WITHIN_CANDIDATES") return within_candidates;
    if (!audit_complete || !same_winner || !same_block_genotypes) return "UNRESOLVED";
    return audit_phase_support == "READS_WITHIN_CANDIDATES" ? "READS_AUDITED" : "PRIOR_ONLY";
}

std::vector<BlockGenotypeConfidence> block_genotype_confidence(
    const std::vector<double>& scores,
    const std::vector<double>& error_bounds,
    const std::vector<std::vector<std::pair<std::size_t, std::size_t>>>& genotypes) {
    if (scores.empty() || error_bounds.size() != scores.size() || genotypes.size() != scores.size()) {
        throw std::invalid_argument("block genotype confidence needs one score, bound and genotype row per candidate");
    }
    const std::size_t blocks = genotypes.front().size();
    double log_weight = -std::numeric_limits<double>::infinity();
    for (const double score : scores) log_weight = log_add_exp(log_weight, score);
    std::vector<double> log_other(blocks, -std::numeric_limits<double>::infinity());
    std::vector<BlockGenotypeConfidence> out(blocks);
    for (std::size_t c = 1; c < scores.size(); ++c) {
        if (genotypes[c].size() != blocks) throw std::invalid_argument("every candidate needs a genotype per block");
        for (std::size_t b = 0; b < blocks; ++b) {
            if (genotypes[c][b] == genotypes[0][b]) continue;
            log_other[b] = log_add_exp(log_other[b], scores[c]);
            if (scores[c] + error_bounds[c] >= scores[0] - error_bounds[0]) out[b].numerically_separated = false;
        }
    }
    for (std::size_t b = 0; b < blocks; ++b) {
        const double log_error = log_other[b] - log_weight;
        out[b].posterior = -std::expm1(log_error);
        out[b].quality = -10.0 * log_error / std::log(10.0);
    }
    return out;
}

namespace {

void remove_run_outputs(const std::string& out_prefix) {
    for (const char* suffix : {".genotypes.tsv", ".block_genotypes.tsv", ".mosaic_candidates.tsv",
                               ".mosaic_scores.tsv", ".mosaic_status.tsv", ".mosaic_haplotypes.fa"}) {
        std::remove((out_prefix + suffix).c_str());
    }
}

// The proposal a standalone run with target k returns. Neither the beam's pruning nor its ordering reads
// the target: it decides only the final truncation and the two flags below, recomputed here exactly as
// propose_mosaic_diplotypes computes them. Gated against standalone runs end to end.
MosaicBeamResult candidate_prefix(const MosaicBeamResult& beam, std::size_t k) {
    MosaicBeamResult prefix = beam;
    prefix.unique_target_reached = beam.unique_candidates_found >= k;
    prefix.candidate_set_conditional = !beam.search_exhaustive || beam.unique_candidates_found > k;
    if (prefix.candidates.size() > k) prefix.candidates.resize(k);
    return prefix;
}

} // namespace

GenotypeRunStatus run_genotype(const GenotypeCommandOptions& options) {
    const Clock::time_point started = Clock::now();
    StageTimings timings;
    const std::string status_path = options.out_prefix + ".mosaic_status.tsv";
    // A run owns its output prefix. Outputs an earlier run left under the same prefix are removed before
    // anything else, so a run that stops early can never sit beside a previous ranking and look complete.
    remove_run_outputs(options.out_prefix);
    if (!options.development_candidate_prefixes.empty()) {
        const std::vector<std::size_t>& prefixes = options.development_candidate_prefixes;
        for (std::size_t i = 0; i < prefixes.size(); ++i) {
            if (i > 0 && prefixes[i] <= prefixes[i - 1]) {
                throw std::runtime_error("development candidate prefixes must be strictly increasing");
            }
            remove_run_outputs(options.out_prefix + ".K" + std::to_string(prefixes[i]));
        }
        if (prefixes.back() != options.max_mosaic_candidates) {
            throw std::runtime_error("the largest development candidate prefix must equal max_mosaic_candidates");
        }
    }

    // Version one's scorer is serial. The option existed and was silently ignored; a request it cannot honour is
    // now refused rather than reported back as threads_used 1.
    if (options.threads != 1) {
        throw std::runtime_error("version one is serial: --threads must be 1");
    }

    if (options.phase_audit_candidates > 0) {
        const std::size_t largest = options.development_candidate_prefixes.empty()
                                        ? options.max_mosaic_candidates
                                        : options.development_candidate_prefixes.back();
        // At a primary set as large as the audit, the audit adds no rival and confirms only what the set shows.
        if (options.phase_audit_candidates < largest) {
            throw std::runtime_error("the phase audit set must be at least as large as every primary candidate set");
        }
    }

    if (options.gfa_path.empty() || options.reads1_path.empty() || options.reads2_path.empty() ||
        options.reference_path.empty() || options.library_path.empty()) {
        throw std::runtime_error(
            "genotype requires --gfa, --reads1, --reads2, --reference-path and --library");
    }

    std::optional<double> declared_marker_background;
    const LibraryModel library = read_library(options.library_path, declared_marker_background);

    // ---- stage 1: structural model -------------------------------------------------------------
    // Fail closed. Anything that stops the structural model being built, verified, or accepted as
    // representation-independent ends the run as MODEL_INCOMPLETE with the reason, rather than escaping
    // as an exception that leaves no status behind.
    Clock::time_point mark = Clock::now();
    Graph graph;
    LocusInputResult built;
    std::optional<PreparedGenotypeModel> prepared;
    std::string fingerprint = "NA";
    try {
        GfaModel gfa = read_gfa_model(options.gfa_path);
        GraphSortOptions sort_options;
        sort_options.reference_path = options.reference_path;
        sort_options.flip = true;
        const std::string resolved = sort_graph_reference(gfa, sort_options).resolved_reference;
        const std::string staged = options.out_prefix + ".sorted.gfa";
        write_gfa_model(staged, gfa);
        ParseGfaOptions parse_options;
        graph = parse_gfa(staged, parse_options);
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
        built = build_locus_input(graph, bubbles, locus_options);
        prepared.emplace(prepare_genotype_model(built.input));
        // The authoritative check, on real input, before any statistic is computed from these blocks.
        verify_complete_panel_round_trip(*prepared, built.authoritative);
        fingerprint = model_fingerprint(*prepared);
        // A model a node split would change is not scored. The round trip cannot see this case.
        require_representation_independent(built.report);
    } catch (const std::exception& error) {
        timings.model = seconds_since(mark);
        timings.total = seconds_since(started);
        write_status(status_path, GenotypeRunStatus::ModelIncomplete, fingerprint, options, library,
                     timings, built.input.blocks.empty() ? nullptr : &built.report, nullptr, nullptr,
                     nullptr, std::string("the structural model could not be built: ") + error.what());
        return GenotypeRunStatus::ModelIncomplete;
    }
    const PreparedGenotypeModel& model = *prepared;
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
    // Per-copy depth from the invariant anchors, whose diploid copy number is known whatever the
    // sample's genotype. Variable markers cannot be used: their expected count depends on the very
    // genotype being inferred.
    MarkerDepthEstimate depth;
    try {
        depth = estimate_marker_depth(panel, observations, declared_marker_background);
    } catch (const std::invalid_argument& error) {
        timings.markers = seconds_since(mark);
        timings.total = seconds_since(started);
        write_status(status_path, GenotypeRunStatus::CandidateIncomplete, fingerprint, options,
                     library, timings, &built.report, &panel, nullptr, nullptr,
                     std::string("marker depth could not be estimated: ") + error.what());
        return GenotypeRunStatus::CandidateIncomplete;
    }
    if (!depth.converged) {
        // A depth estimate that did not converge is not fed to the proposal: it changes which candidates
        // reach K. NA19240 used to stop at its iteration cap with a rate 3.4% high and still proceed.
        timings.markers = seconds_since(mark);
        timings.total = seconds_since(started);
        write_status(status_path, GenotypeRunStatus::CandidateIncomplete, fingerprint, options, library,
                     timings, &built.report, &panel, nullptr, nullptr,
                     "marker depth fit did not converge (" + depth.termination +
                         "); declare marker_background to use the one-dimensional fit",
                     nullptr, &depth);
        return GenotypeRunStatus::CandidateIncomplete;
    }
    // The fitted or declared background, unchanged: the fit and the emission share one model floor. This
    // used to raise a fitted 1e-9 to 1e-6 here, a log(1000) change per zero-copy marker observation.
    const MarkerModelParameters marker_parameters =
        proposal_marker_parameters(depth, options.kmer_size, options.syncmer_s);
    const PreparedProposalEmissions emissions =
        build_proposal_emissions(panel, observations, marker_parameters);
    timings.markers = seconds_since(mark);

    // ---- stage 3a: protected panel baseline --------------------------------------------------------------
    // Every unordered pair of panel templates is a candidate, whatever the beam retains. In leave-zero-out the
    // truth is such a pair, and the beam can prune it before later informative blocks: on four held-out C4
    // donors the true pair was dropped where template copies of a few locally better prefixes filled all 1024
    // merged states. The beam now supplies only recombinant mosaics beyond the panel, and K counts those.
    std::vector<MosaicCandidate> panel_candidates = panel_diplotypes(model.panel());
    for (MosaicCandidate& candidate : panel_candidates) candidate.every_allele_scored = emissions.every_allele_scored();
    CandidateSetReport candidate_report;
    candidate_report.panel_diplotypes = panel_candidates.size();
    // Predicted work is checked before any is done: one diploid combination per panel pair per fragment. A
    // baseline over budget is refused whole; a partial baseline is never scored or described as complete.
    if (!fragments.empty() &&
        panel_candidates.size() > options.max_panel_pair_fragment_evaluations / fragments.size()) {
        timings.total = seconds_since(started);
        write_status(status_path, GenotypeRunStatus::CandidateIncomplete, fingerprint, options, library,
                     timings, &built.report, &panel, nullptr, nullptr,
                     "the panel-diplotype baseline would exceed its declared scoring budget (" +
                         std::to_string(panel_candidates.size()) + " pairs x " +
                         std::to_string(fragments.size()) + " fragments); no partial baseline is scored",
                     nullptr, &depth, &candidate_report);
        return GenotypeRunStatus::CandidateIncomplete;
    }
    candidate_report.panel_baseline_complete = true;

    // A true panel-only run: no mosaics, no phase audit, and no positive development prefix. The beam
    // would propose nothing that is ever scored, so it is not built at all -- and its pre-allocation
    // budget then cannot refuse a call that never needed it. That refusal is what kept acot, gstm1 and
    // ankrd36c from running an arm whose own work fits every declared cap with room to spare.
    bool positive_prefix_requested = false;
    for (const std::size_t k : options.development_candidate_prefixes) {
        if (k > 0) positive_prefix_requested = true;
    }
    const bool panel_only = options.max_mosaic_candidates == 0 &&
                            options.phase_audit_candidates == 0 && !positive_prefix_requested;
    candidate_report.panel_only = panel_only;

    // ---- stage 3b: mosaic proposal --------------------------------------------------------------------------
    mark = Clock::now();
    const std::vector<double> switch_probability(
        model.blocks().size() - 1, options.switch_probability);
    MosaicBeamResult beam;
    if (!panel_only) {
        MosaicBeamParameters beam_parameters;
        // The beam returns every diplotype it found -- never more than its state limit -- so the number of
        // off-panel mosaics below is exact, not a count truncated by K.
        beam_parameters.unique_candidate_target =
            std::max(options.proposal_state_limit, options.max_mosaic_candidates + panel_candidates.size());
        beam_parameters.state_limit = options.proposal_state_limit;
        beam_parameters.max_extension_evaluations = options.max_beam_extensions;
        try {
            beam = propose_mosaic_diplotypes(model.panel(), emissions, switch_probability, beam_parameters);
        } catch (const std::length_error& error) {
            timings.proposal = seconds_since(mark);
            timings.total = seconds_since(started);
            write_status(status_path, GenotypeRunStatus::CandidateIncomplete, fingerprint, options, library,
                         timings, &built.report, &panel, nullptr, nullptr,
                         std::string("the proposal exceeded its declared resource budget: ") + error.what(),
                         nullptr, &depth, &candidate_report);
            return GenotypeRunStatus::CandidateIncomplete;
        }
    }
    timings.proposal = seconds_since(mark);
    candidate_report.beam_diplotypes_found = beam.unique_candidates_found;

    // Off-panel mosaics in proposal order; panel pairs the beam also found are already candidates.
    const MosaicBeamResult mosaics = off_panel_mosaics(beam, panel_candidates);

    // ---- stage 4: one spelling and scoring pass over every candidate a ranked prefix holds --------------
    // A diplotype's fragment score does not depend on which other candidates are scored with it, so the panel
    // and the largest ranked mosaic prefix are scored once and every prefix is certified from those scores.
    const std::vector<std::size_t> targets = options.development_candidate_prefixes.empty()
                                                 ? std::vector<std::size_t>{options.max_mosaic_candidates}
                                                 : options.development_candidate_prefixes;
    // Every target is ranked -- a mosaic-deficient one too, conditionally (review round 6) -- so the union holds
    // the panel and the largest set of mosaics any primary, prefix or audit set contains.
    std::size_t scored_mosaics = 0;
    const bool any_ranked = true;
    std::vector<std::size_t> score_targets = targets;
    if (options.phase_audit_candidates > 0) score_targets.push_back(options.phase_audit_candidates);
    for (const std::size_t k : score_targets) {
        scored_mosaics = std::max(scored_mosaics, candidate_prefix(mosaics, k).candidates.size());
    }
    mark = Clock::now();
    std::vector<MosaicPath> paths;
    std::vector<std::string> sequences;
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
    std::vector<const MosaicCandidate*> scored_union;
    if (any_ranked) {
        for (const MosaicCandidate& candidate : panel_candidates) scored_union.push_back(&candidate);
        for (std::size_t i = 0; i < scored_mosaics; ++i) scored_union.push_back(&mosaics.candidates[i]);
        for (const MosaicCandidate* candidate : scored_union) {
            slot_for(candidate->diplotype.first);
            slot_for(candidate->diplotype.second);
        }
    }
    const InternedHaplotypes interned = any_ranked ? intern_haplotypes(paths, sequences) : InternedHaplotypes{};
    // Candidates are ranked as unordered pairs of SEQUENCES. Two candidates whose allele paths differ but
    // spell the same two haplotypes are one state to any read-based method.
    const auto sequence_pair = [&](const MosaicCandidate& candidate) {
        std::size_t x = interned.index_of[slot_of_path.at(candidate.diplotype.first.alleles)];
        std::size_t y = interned.index_of[slot_of_path.at(candidate.diplotype.second.alleles)];
        if (y < x) std::swap(x, y);
        return std::make_pair(x, y);
    };
    std::map<std::pair<std::size_t, std::size_t>, std::size_t> scored_index;
    std::vector<HaplotypeCandidate> scored_pairs;
    for (const MosaicCandidate* candidate : scored_union) {
        const auto key = sequence_pair(*candidate);
        if (scored_index.emplace(key, scored_pairs.size()).second) scored_pairs.push_back({key.first, key.second});
    }
    // Exact haploid sequence priors, independent of the beam: each sums every allele path spelling the bytes.
    std::vector<double> haploid_sequence_prior(interned.sequences.size());
    for (std::size_t h = 0; h < interned.sequences.size(); ++h) {
        haploid_sequence_prior[h] = log_haploid_sequence_prior(model, interned.sequences[h], switch_probability);
    }
    timings.spelling = seconds_since(mark);

    mark = Clock::now();
    FragmentScoringOptions scoring_options;
    scoring_options.mismatch_band = options.mismatch_band;
    scoring_options.minimum_piece_bp = options.minimum_piece_bp;
    scoring_options.max_piece_hits = options.max_piece_hits;
    scoring_options.max_proposed_starts_per_haplotype = options.max_proposed_starts_per_haplotype;
    scoring_options.max_physical_origins = options.max_physical_origins;
    scoring_options.max_state_classes = options.max_state_classes;
    scoring_options.max_mass_cells = options.max_mass_cells;
    FragmentScoringResult scored;
    std::string rescore_failure;
    if (any_ranked) {
        try {
            scored = score_diplotypes(interned, scored_pairs, fragments, library, scoring_options);
        } catch (const ResourceExhausted& error) {
            rescore_failure = error.what();
        }
    }
    timings.rescoring = seconds_since(mark);

    // Certification of one candidate set -- the panel baseline plus a mosaic prefix -- from the shared scores.
    struct SequenceCandidate {
        std::size_t first = 0;     // interned sequence, first <= second
        std::size_t second = 0;
        std::vector<std::size_t> members;
        std::size_t representative = 0;
        bool projection_ambiguous = false;
    };
    struct SetRanking {
        std::vector<MosaicCandidate> candidates;
        std::vector<SequenceCandidate> sequence_candidates;
        std::vector<std::size_t> sequence_candidate_of;
        std::vector<double> path_prior;
        CertifiedRanking ranking;
        RankingSummary summary;
    };
    const auto rank_set = [&](const MosaicBeamResult& beam) {
        // The candidate set: the whole panel baseline first, then this target's mosaics.
        std::vector<MosaicCandidate> candidates = panel_candidates;
        candidates.insert(candidates.end(), beam.candidates.begin(), beam.candidates.end());

        std::vector<SequenceCandidate> sequence_candidates;
        std::map<std::pair<std::size_t, std::size_t>, std::size_t> slot_of_sequence_pair;
        std::vector<std::size_t> sequence_candidate_of(candidates.size(), 0);
        for (std::size_t i = 0; i < candidates.size(); ++i) {
            const auto key = sequence_pair(candidates[i]);
            const auto inserted = slot_of_sequence_pair.emplace(key, sequence_candidates.size());
            if (inserted.second) {
                SequenceCandidate fresh;
                fresh.first = key.first;
                fresh.second = key.second;
                sequence_candidates.push_back(fresh);
            }
            sequence_candidates[inserted.first->second].members.push_back(i);
            sequence_candidate_of[i] = inserted.first->second;
        }

        // ---- stage 5: exact prior, recomputed independently of the beam ------------------------------
        // The path prior is a manifest diagnostic; the ranking uses the sequence-pair prior, which decides
        // the log(2) multiplicity on bytes.
        std::vector<double> path_prior;
        path_prior.reserve(candidates.size());
        for (const MosaicCandidate& candidate : candidates) {
            path_prior.push_back(
                log_unordered_diplotype_prior(model.panel(), candidate.diplotype, switch_probability));
        }
        std::vector<double> sequence_prior;
        std::vector<DiplotypeScore> scores;
        for (SequenceCandidate& candidate : sequence_candidates) {
            sequence_prior.push_back(haploid_sequence_prior[candidate.first] +
                                     haploid_sequence_prior[candidate.second] +
                                     (candidate.first != candidate.second ? std::log(2.0) : 0.0));
            scores.push_back(scored.scores[scored_index.at({candidate.first, candidate.second})]);

            // Projection needs one allele path per haplotype. Take the member with the highest exact path
            // prior, ties to candidate order, and flag the call if any member projects to different block
            // alleles: those bytes then do not determine the block genotype.
            std::size_t best = candidate.members.front();
            for (const std::size_t member : candidate.members) {
                if (path_prior[member] > path_prior[best]) best = member;
            }
            candidate.representative = best;
            const MosaicDiplotype& chosen = candidates[best].diplotype;
            for (const std::size_t member : candidate.members) {
                const MosaicDiplotype& other = candidates[member].diplotype;
                for (std::size_t b = 0; b < model.blocks().size(); ++b) {
                    auto lhs = std::minmax(chosen.first.alleles[b], chosen.second.alleles[b]);
                    auto rhs = std::minmax(other.first.alleles[b], other.second.alleles[b]);
                    if (lhs != rhs) candidate.projection_ambiguous = true;
                }
            }
        }

        // ---- stage 6: certification of this set's scores ---------------------------------------------
        const CertifiedRanking ranking = certify_scored_ranking(
            scores, fragments.size(), scoring_options.mismatch_band, scored.counters, CertificationOptions{},
            sequence_prior);

        RankingSummary summary;
        summary.sequence_candidates = sequence_candidates.size();
        summary.collapsed_beam_candidates = candidates.size() - sequence_candidates.size();
        summary.block_genotype_read_supported.assign(model.blocks().size(), true);
        {
            const CertifiedCandidate& best = ranking.candidates.front();
            const MosaicDiplotype& winner_diplotype =
                candidates[sequence_candidates[best.candidate_id].representative].diplotype;
            summary.fragment_term_gap = ranking.candidates.size() > 1
                                            ? best.fragment_term - ranking.candidates[1].fragment_term
                                            : std::numeric_limits<double>::infinity();
            for (std::size_t i = 1; i < ranking.candidates.size(); ++i) {
                const CertifiedCandidate& rival = ranking.candidates[i];
                // A fragment contender: the reads alone do not put it below the winner.
                if (rival.fragment_term + rival.numeric_error_bound <
                    best.fragment_term - best.numeric_error_bound) {
                    continue;
                }
                ++summary.fragment_contenders;
                const MosaicDiplotype& other =
                    candidates[sequence_candidates[rival.candidate_id].representative].diplotype;
                for (std::size_t b = 0; b < model.blocks().size(); ++b) {
                    if (std::minmax(winner_diplotype.first.alleles[b], winner_diplotype.second.alleles[b]) !=
                        std::minmax(other.first.alleles[b], other.second.alleles[b])) {
                        summary.block_genotype_read_supported[b] = false;
                    }
                }
            }
            // A ranked set with no off-panel rival and no audit cannot speak to phase. Panel diplotypes
            // are unordered pairs of whole templates, so switching phase within a pair is not even
            // representable without a mosaic: having no fragment contender is then a property of the
            // candidate set, not evidence about phase. The condition is on what was RANKED, not on
            // whether a beam happened to be built -- a K = 0 view of a ladder ranks the same set as a
            // standalone K = 0 run and must make the same claim. A K = 0 arm carrying an audit is
            // different: the audit does compare against a larger set, so it keeps its claim.
            const bool no_rival_to_disagree =
                beam.candidates.empty() && options.phase_audit_candidates == 0;
            summary.phase_support = no_rival_to_disagree                        ? "UNRESOLVED"
                                    : ranking.status != RescoreStatus::Complete ? "UNRESOLVED"
                                    : summary.fragment_contenders == 0          ? "READS_WITHIN_CANDIDATES"
                                                                                : "PRIOR_ONLY";
            summary.projection_ambiguous = sequence_candidates[best.candidate_id].projection_ambiguous;
        }
        {
            // Per block, conditional on the retained candidate set: posterior and GQ of the winner's unordered
            // genotype, grouping candidates by their representative's genotype there.
            std::vector<double> ranked_scores;
            std::vector<double> ranked_bounds;
            std::vector<std::vector<std::pair<std::size_t, std::size_t>>> ranked_genotypes;
            for (const CertifiedCandidate& entry : ranking.candidates) {
                ranked_scores.push_back(entry.score);
                ranked_bounds.push_back(entry.numeric_error_bound);
                const MosaicDiplotype& d = candidates[sequence_candidates[entry.candidate_id].representative].diplotype;
                std::vector<std::pair<std::size_t, std::size_t>> row;
                for (std::size_t b = 0; b < model.blocks().size(); ++b) {
                    row.push_back(std::minmax(d.first.alleles[b], d.second.alleles[b]));
                }
                ranked_genotypes.push_back(std::move(row));
            }
            for (const BlockGenotypeConfidence& block :
                 block_genotype_confidence(ranked_scores, ranked_bounds, ranked_genotypes)) {
                summary.block_genotype_posterior.push_back(block.posterior);
                summary.block_genotype_quality.push_back(block.quality);
                summary.block_numerically_separated.push_back(block.numerically_separated);
            }
        }

        SetRanking out;
        out.candidates = std::move(candidates);
        out.sequence_candidates = std::move(sequence_candidates);
        out.sequence_candidate_of = std::move(sequence_candidate_of);
        out.path_prior = std::move(path_prior);
        out.ranking = ranking;
        out.summary = std::move(summary);
        return out;
    };

    // Everything after scoring, for one mosaic target. A standalone run calls it once; the development prefix
    // harness calls it once per prefix of the same proposal and scoring pass, so each prefix is refused,
    // certified and written by exactly the code a standalone run uses.
    const auto finish = [&](const MosaicBeamResult& beam, const GenotypeCommandOptions& options,
                            StageTimings timings) -> GenotypeRunStatus {
        const std::string status_path = options.out_prefix + ".mosaic_status.tsv";
        const std::string manifest_header =
            "candidate_id\tunique_rank\thap1_alleles\thap2_alleles\tbest_proposal_score\t"
            "retained_logsum_proposal_score\tproposal_emission_component\t"
            "proposal_transition_component\tequivalent_retained_histories\tevery_allele_scored\t"
            "exact_path_prior\tsequence_candidate_id\thistory_count_saturated\tsource\n";
        // A pruned search that found fewer off-panel mosaics than requested has omitted mosaics it never proposed,
        // and the status says so (mosaic_target_reached false, candidate_scope CONDITIONAL). Every panel diplotype
        // was still evaluated, so the call is made (review round 6): the panel-supported genotype is usable, and
        // the missing mosaics are explicit rather than a reason to withhold it.
        if (!rescore_failure.empty()) {
            timings.total = seconds_since(started);
            write_status(status_path, GenotypeRunStatus::RescoreIncomplete, fingerprint, options, library,
                         timings, &built.report, &panel, (panel_only ? nullptr : &beam), nullptr,
                         "rescoring stopped at a declared resource cap and returned no ranking: " +
                             rescore_failure,
                         nullptr, &depth, &candidate_report);
            return GenotypeRunStatus::RescoreIncomplete;
        }

        SetRanking primary = rank_set(beam);
        const std::vector<MosaicCandidate>& candidates = primary.candidates;
        const std::size_t panel_count = panel_candidates.size();
        const std::vector<SequenceCandidate>& sequence_candidates = primary.sequence_candidates;
        const std::vector<std::size_t>& sequence_candidate_of = primary.sequence_candidate_of;
        const std::vector<double>& path_prior = primary.path_prior;
        const CertifiedRanking& ranking = primary.ranking;
        RankingSummary& summary = primary.summary;
        CandidateSetReport set_report = candidate_report;
        {
            std::set<std::size_t> haplotypes;
            for (const SequenceCandidate& candidate : sequence_candidates) {
                haplotypes.insert(candidate.first);
                haplotypes.insert(candidate.second);
            }
            set_report.unique_haplotypes = haplotypes.size();
        }

        // Phase support from the declared audit set. A small candidate set can omit the rival that makes phase a
        // fragment tie: in the round-5 held-out cohort nine of ten READS claims at K = 8 vanished at K = 256. The
        // audit set is ranked from the same scores; READS_AUDITED needs the same winner, the same block genotypes
        // and read-supported phase there. Even so it is conditional on the audit set, not a global proof.
        summary.phase_support_within_candidates = summary.phase_support;
        if (options.phase_audit_candidates > 0) {
            summary.phase_audit_candidates = options.phase_audit_candidates;
            bool audit_complete = false;
            bool same_winner = false;
            bool same_blocks = false;
            const MosaicBeamResult audit_view = candidate_prefix(mosaics, options.phase_audit_candidates);
            if (!audit_view.search_exhaustive && !audit_view.unique_target_reached) {
                // An unreachable audit set leaves phase unresolved; it never withholds the block genotypes.
                summary.phase_audit_status = "MOSAIC_TARGET_UNREACHED";
            } else if (ranking.status != RescoreStatus::Complete) {
                summary.phase_audit_status = "NOT_RUN";
            } else {
                const SetRanking audit = rank_set(audit_view);
                audit_complete = audit.ranking.status == RescoreStatus::Complete;
                summary.phase_audit_status = audit_complete ? "COMPLETE" : "RESCORE_INCOMPLETE";
                summary.phase_audit_phase_support = audit.summary.phase_support;
                if (audit_complete) {
                    const SequenceCandidate& mine = sequence_candidates[ranking.winner_set.front()];
                    const SequenceCandidate& theirs = audit.sequence_candidates[audit.ranking.winner_set.front()];
                    same_winner = mine.first == theirs.first && mine.second == theirs.second;
                    const MosaicDiplotype& a = candidates[mine.representative].diplotype;
                    const MosaicDiplotype& b = audit.candidates[theirs.representative].diplotype;
                    same_blocks = true;
                    for (std::size_t block = 0; block < model.blocks().size(); ++block) {
                        if (std::minmax(a.first.alleles[block], a.second.alleles[block]) !=
                            std::minmax(b.first.alleles[block], b.second.alleles[block])) {
                            same_blocks = false;
                        }
                    }
                    summary.phase_audit_same_winner = same_winner ? "true" : "false";
                    summary.phase_audit_same_block_genotypes = same_blocks ? "true" : "false";
                }
            }
            summary.phase_support = audited_phase_support(summary.phase_support, audit_complete, same_winner,
                                                          same_blocks, summary.phase_audit_phase_support);
        }

        // ---- outputs --------------------------------------------------------------------------------
        {
            std::ofstream out = open_output(options.out_prefix + ".mosaic_candidates.tsv");
            out << manifest_header;
            for (std::size_t i = 0; i < candidates.size(); ++i) {
                const MosaicCandidate& candidate = candidates[i];
                out << candidate.candidate_id << '\t';
                if (i < panel_count) {
                    // A panel diplotype is not a proposal result: it has no rank or proposal score.
                    out << "NA\t" << join_alleles(candidate.diplotype.first.alleles) << '\t'
                        << join_alleles(candidate.diplotype.second.alleles) << "\tNA\tNA\tNA\tNA\tNA\t"
                        << (candidate.every_allele_scored ? "true" : "false") << '\t' << path_prior[i] << '\t'
                        << sequence_candidate_of[i] << "\tNA\tpanel\n";
                } else {
                    out << candidate.unique_rank << '\t'
                        << join_alleles(candidate.diplotype.first.alleles) << '\t'
                        << join_alleles(candidate.diplotype.second.alleles) << '\t'
                        << candidate.best_proposal_score << '\t'
                        << candidate.retained_logsum_proposal_score << '\t'
                        << candidate.proposal_emission_component << '\t'
                        << candidate.proposal_transition_component << '\t'
                        << candidate.equivalent_retained_histories << '\t'
                        << (candidate.every_allele_scored ? "true" : "false") << '\t'
                        << path_prior[i] << '\t' << sequence_candidate_of[i] << '\t'
                        << (candidate.history_count_saturated ? "true" : "false") << "\tmosaic\n";
                }
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
                    << candidate.members.size() << '\t' << entry.score << '\t'
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
                         library, timings, &built.report, &panel, (panel_only ? nullptr : &beam), &ranking,
                         "the best candidate cannot be separated from every rival by more than "
                         "floating-point error; the equivalence set is reported instead of a call",
                         &summary, &depth, &set_report);
            return GenotypeRunStatus::RescoreIncomplete;
        }

        const std::size_t winner = ranking.winner_set.front();
        const MosaicDiplotype& called = candidates[sequence_candidates[winner].representative].diplotype;
        {
            // Only variant blocks are reporting units. Context blocks are modelled identically and are
            // what let the mosaic recombine, but they are not SV calls and are not presented as such.
            std::ofstream out = open_output(options.out_prefix + ".genotypes.tsv");
            out << "block_id\tkind\tallele1\tallele2\tallele1_bp\tallele2_bp\tallele1_sources\t"
                   "allele2_sources\tzygosity\tconditional_genotype_posterior\tconditional_GQ\n";
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
                    << (a1 == a2 ? "HOM" : "HET") << '\t' << std::setprecision(12)
                    << summary.block_genotype_posterior[b] << '\t' << summary.block_genotype_quality[b]
                    << std::setprecision(6) << '\n';
            }
        }
        {
            // Every statistical block, for block-level evaluation. Context blocks are polymorphic and are
            // genotyped by the same mosaic; `reportable` marks the SV bubbles the SV-facing genotypes table
            // carries, so the interstitial blocks are distinguished rather than dropped.
            std::ofstream out = open_output(options.out_prefix + ".block_genotypes.tsv");
            out << "block_id\tkind\treportable\tallele1\tallele2\tallele1_bp\tallele2_bp\tzygosity\t"
                   "projection_ambiguous\tgenotype_support\tphase_support\tconditional_genotype_posterior\t"
                   "conditional_GQ\tnumerically_separated\tcandidate_scope\n";
            const bool ambiguous = sequence_candidates[winner].projection_ambiguous;
            for (std::size_t b = 0; b < model.blocks().size(); ++b) {
                const PreparedGenotypeBlock& block = model.blocks()[b];
                const bool variant = block.kind == GenotypeBlockKind::Variant;
                const std::size_t a1 = called.first.alleles[b];
                const std::size_t a2 = called.second.alleles[b];
                out << block.block_id << '\t' << (variant ? "SV" : "context") << '\t'
                    << (variant ? "true" : "false") << '\t' << a1 << '\t' << a2 << '\t'
                    << block.alleles[a1].sequence.size() << '\t' << block.alleles[a2].sequence.size() << '\t'
                    << (a1 == a2 ? "HOM" : "HET") << '\t' << (ambiguous ? "true" : "false") << '\t'
                    << (summary.block_genotype_read_supported[b] ? "READS" : "PRIOR_ONLY") << '\t'
                    << summary.phase_support << '\t' << std::setprecision(12)
                    << summary.block_genotype_posterior[b] << '\t' << summary.block_genotype_quality[b]
                    << std::setprecision(6) << '\t' << (summary.block_numerically_separated[b] ? "true" : "false")
                    << '\t' << (beam.candidate_set_conditional ? "CONDITIONAL" : "COMPLETE") << '\n';
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
                     &built.report, &panel, (panel_only ? nullptr : &beam), &ranking, "", &summary,
                     &depth, &set_report);
        if (!options.quiet) {
            std::cerr << "[genotype] COMPLETE  fingerprint=" << fingerprint
                      << "  candidates=" << candidates.size()
                      << "  separation=" << ranking.certified_gap
                      << "  total=" << timings.total << "s\n";
        }
        return GenotypeRunStatus::Complete;
    };

    if (options.development_candidate_prefixes.empty()) {
        return finish(candidate_prefix(mosaics, options.max_mosaic_candidates), options, timings);
    }
    GenotypeRunStatus last = GenotypeRunStatus::Complete;
    for (const std::size_t k : options.development_candidate_prefixes) {
        GenotypeCommandOptions prefix_options = options;
        prefix_options.max_mosaic_candidates = k;
        prefix_options.out_prefix = options.out_prefix + ".K" + std::to_string(k);
        prefix_options.development_candidate_prefixes.clear();
        remove_run_outputs(prefix_options.out_prefix);
        last = finish(candidate_prefix(mosaics, k), prefix_options, timings);
    }
    return last;
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
        else if (flag == "--threads") options.threads = std::stoul(need(i));   // refused unless 1, see run_genotype
        else if (flag == "--max-mosaic-candidates") options.max_mosaic_candidates = std::stoul(need(i));
        else if (flag == "--max-candidates") {
            throw std::runtime_error("--max-candidates is now --max-mosaic-candidates: every panel diplotype is "
                                     "always a candidate, and K counts only the off-panel mosaics added to it");
        }
        else if (flag == "-q" || flag == "--quiet") options.quiet = true;
        else if (flag == "-h" || flag == "--help") {
            std::cout <<
                "Usage:\n"
                "  panvar genotype --gfa <graph.gfa> --reference-path <name> \\\n"
                "                  --reads1 <r1.fq> --reads2 <r2.fq> --library <params.txt> \\\n"
                "                  [--bubbles <bubbles.csv>] [--out-prefix <prefix>] \\\n"
                "                  [--max-mosaic-candidates <K>]\n\n"
                "Calls a pair of whole-locus mosaic haplotypes and projects them to one allele\n"
                "pair per SV block. Paired reads only.\n\n"
                "Every unordered pair of panel haplotypes is always a candidate. --max-mosaic-candidates\n"
                "adds up to K recombinant mosaics from the proposal beam (0: the panel alone); when the\n"
                "beam finds fewer, the call is still made and marked CONDITIONAL.\n\n"
                "Outputs <prefix>.genotypes.tsv, .block_genotypes.tsv, .mosaic_candidates.tsv,\n"
                ".mosaic_scores.tsv,"
                ".mosaic_status.tsv and .mosaic_haplotypes.fa. An incomplete run emits its\n"
                "diagnostics and says so in the status file; it never presents a partial ranking\n"
                "as a result.\n\n"
                "--library declares insert_mean_bp, insert_sd_bp, support_sigmas,\n"
                "substitution_error, outlier_weight and fragment_start_intensity, one per line,\n"
                "and optionally marker_background (estimated from invariant anchors when absent).\n"
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
