#include "panvar/genotype_command.hpp"

#include "panvar/bubbles.hpp"
#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_locus.hpp"
#include "panvar/genotype_model.hpp"
#include "panvar/genotype_repair.hpp"
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
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>

namespace panvar {
namespace {

using Clock = std::chrono::steady_clock;

double seconds_since(Clock::time_point from) {
    return std::chrono::duration<double>(Clock::now() - from).count();
}

// Every output a run can leave under its prefix, final and partial. A run owns its prefix: these are
// removed before anything else, so a run that stops early can never sit beside an earlier result.
const char* const kFinalOutputs[] = {".genotypes.tsv", ".block_genotypes.tsv", ".mosaic_candidates.tsv",
                                     ".mosaic_scores.tsv", ".mosaic_haplotypes.fa", ".equivalence_set.tsv",
                                     ".repair_rounds.tsv", ".repair_memory.tsv"};
const char* const kPartialOutputs[] = {".repair_rounds.partial.tsv", ".repair_memory.partial.tsv"};

// The run's process-level guard. Status files are written whole (atomic rename), and the terminal
// status, or the start of the final outputs, CLAIMS the run under the mutex: from then on the
// watchdog cannot stop it, and before then a watchdog stop removes every final output first. So a
// stopped run can never leave a final-named table, and a completed one is never stamped as stopped.
struct RunWatch {
    std::mutex mutex;
    bool finished = false;
    std::string out_prefix;
    std::string fingerprint = "NA";
    std::string profile = "none";
    std::string profile_sha256 = "none";
    std::string phase = "starting";
    std::size_t repair_rounds_committed = 0;
    Clock::time_point started = Clock::now();
    std::uint64_t footprint_limit = 0;
    double time_limit = 0.0;
};
// The one run in flight in this process, or null. Runs are sequential within a process.
RunWatch* g_watch = nullptr;
// Provenance of the inputs, by content, for the run in flight.
std::vector<std::pair<std::string, std::string>> g_input_hashes;

std::string minimal_status(const RunWatch& watch, const std::string& state, const std::string& note) {
    const FootprintSample footprint = physical_footprint();
    std::ostringstream text;
    text << std::setprecision(12) << "key\tvalue\nstatus\t" << state << "\nnote\t" << note
         << "\ncall_class\t" << (state == "RUNNING" ? "pending" : "no_call") << "\nno_call_reason\t"
         << (state == "TIME_LIMIT" ? "time_limit" : state == "RESOURCE_REFUSAL" ? "footprint_limit" : "NA")
         << "\nmodel_fingerprint\t" << watch.fingerprint << "\nprofile\t" << watch.profile << "\nprofile_sha256\t"
         << watch.profile_sha256 << "\nphase\t" << watch.phase << "\nelapsed_seconds\t" << seconds_since(watch.started)
         << "\ntime_limit_seconds\t" << watch.time_limit << "\nfootprint_limit_bytes\t" << watch.footprint_limit
         << "\nphys_footprint_bytes\t" << footprint.current << "\nphys_footprint_peak_bytes\t" << footprint.peak
         << "\nrss_current_bytes\t" << current_rss_bytes() << "\nrss_peak_bytes\t"
         << static_cast<std::uint64_t>(peak_rss_gib() * 1073741824.0) << "\nrepair_rounds_committed\t"
         << watch.repair_rounds_committed << "\noutputs\t"
         << (state == "RUNNING" ? "IN PROGRESS: nothing under this prefix is a result yet"
                                : "PARTIAL: <prefix>.repair_rounds.partial.tsv holds committed rounds only; no "
                                  "genotype, not a call")
         << '\n';
    for (const auto& hash : g_input_hashes) text << "sha256_" << hash.first << '\t' << hash.second << '\n';
    return text.str();
}

// A RUNNING status, only while nothing has claimed the run.
void write_running_status(const std::string& phase) {
    if (g_watch == nullptr) return;
    std::lock_guard<std::mutex> lock(g_watch->mutex);
    if (g_watch->finished) return;
    g_watch->phase = phase;
    atomic_write(g_watch->out_prefix + ".mosaic_status.tsv", minimal_status(*g_watch, "RUNNING", phase));
}

// Claims the run: the watchdog can no longer stop it. Idempotent.
void claim_run() {
    if (g_watch == nullptr) return;
    std::lock_guard<std::mutex> lock(g_watch->mutex);
    g_watch->finished = true;
}

// Stops the process where it stands, if nothing has claimed the run. _Exit runs no destructor, so the
// staged graph and any final output are removed here; committed partial rounds stay as the record.
void stop_run(const std::string& state, const std::string& note, int code) {
    std::lock_guard<std::mutex> lock(g_watch->mutex);
    if (g_watch->finished) return;
    g_watch->finished = true;
    std::error_code ignored;
    for (const char* suffix : kFinalOutputs) std::filesystem::remove(g_watch->out_prefix + suffix, ignored);
    std::filesystem::remove(g_watch->out_prefix + ".sorted.gfa", ignored);
    atomic_write(g_watch->out_prefix + ".mosaic_status.tsv", minimal_status(*g_watch, state, note));
    std::cerr << "[genotype] " << state << ": " << note << '\n';
    std::_Exit(code);
}

// The status of a run is written whole, and a terminal status claims the run.
void commit_status(const std::string& path, const std::string& text) {
    if (g_watch != nullptr) {
        std::lock_guard<std::mutex> lock(g_watch->mutex);
        g_watch->finished = true;
        atomic_write(path, text);
        return;
    }
    atomic_write(path, text);
}

// One panel template name per line; '#' begins a comment and blank lines are ignored. Deliberately
// the dumbest format that can be reviewed in a diff and hashed: a curation decision is a fact about a
// panel, and it should be readable by whoever has to defend it later. A manifest that cannot be opened
// is an error rather than an empty curation -- the alternative silently scores the uncurated panel.
std::vector<std::string> read_curation_manifest(const std::string& path) {
    std::ifstream file(path);
    if (!file) {
        throw std::invalid_argument("curation manifest could not be opened: " + path);
    }
    std::vector<std::string> names;
    std::string line;
    while (std::getline(file, line)) {
        const std::size_t first = line.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) continue;
        // A '#' comments out a line only when it OPENS the line. Panel template names contain '#'
        // themselves -- sampleD#1#chr1:0-100, NA18989#1#haplotype1-0000031:36617797-36825011 -- so
        // stripping from the first '#' anywhere truncates every real name to its sample and the
        // manifest then names nothing that exists.
        if (line[first] == '#') continue;
        const std::size_t last = line.find_last_not_of(" \t\r\n");
        names.push_back(line.substr(first, last - first + 1));
    }
    if (names.empty()) {
        throw std::invalid_argument(
            "curation manifest names no templates: " + path +
            "; an empty manifest is indistinguishable from no curation and must not be passed");
    }
    return names;
}

// Named stage timings, because "it took 90 seconds" is not a performance result: the go/no-go gate
// asks which stage spent them.
struct StageTimings {
    double model = 0.0;
    double markers = 0.0;
    double proposal = 0.0;
    double spelling = 0.0;
    double rescoring = 0.0;
    // Accumulated, not assigned: each ranked set -- primary, every development prefix and the phase
    // audit -- passes through these stages, and a run's cost is their sum over all of them.
    double candidate_grouping = 0.0;   // grouping candidates into sequence-identical classes ONLY
    double prior = 0.0;        // exact path and sequence priors
    double ranking = 0.0;      // certification, block posteriors and the audit comparison
    // The three phases inside rescoring, which is 57-93% of a run. Reported separately because
    // "rescoring" alone cannot say whether the cost is finding placements on each haplotype or
    // combining those masses across every diploid candidate -- and that decides what to parallelise.
    double rescoring_fragment_preparation = 0.0;
    double rescoring_placement_and_mass = 0.0;
    double rescoring_diploid_combination = 0.0;
    double total = 0.0;
};

// The candidate set's scope, derived ONCE and shared by the status and every per-block row. PANEL_ONLY
// is not a kind of COMPLETE: both COMPLETE and CONDITIONAL describe what a search returned, and in a
// panel-only run no search happened. These disagreed until now -- the status said PANEL_ONLY while every
// block row said CONDITIONAL, taken from the empty beam -- which is harmless to the call and corrupting
// to anything that summarises the rows.
const char* candidate_scope_name(bool panel_only, bool candidate_set_conditional) {
    if (panel_only) return "PANEL_ONLY";
    return candidate_set_conditional ? "CONDITIONAL" : "COMPLETE";
}

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
    // Exact-path prior work, which R11 measured at 24-28 s per run on a 108,000-diplotype panel. The
    // haploid prior depends only on a complete allele path, so it is computed once per DISTINCT path and
    // every diplotype is formed from two cached values. `haploid_prior_evaluations` must equal
    // `distinct_haplotype_paths`; if it exceeds it, the cache is not being shared across ranked sets.
    std::size_t distinct_haplotype_paths = 0;
    std::size_t haploid_prior_evaluations = 0;
    std::size_t diploid_prior_lookups = 0;
    // Threading telemetry for the diploid combination phase: what was actually used, which is not what
    // was asked for when there are fewer candidates than workers, and the work done regardless.
    std::size_t threads_used_combination = 1;
    // Reported separately from combination: the two phases parallelise independently, and placement
    // falls back to one worker for full-origin recording even when more were asked for.
    std::size_t threads_used_placement = 1;
    std::uint64_t candidate_fragment_evaluations = 0;
    // Work actually performed when a cap refused the run: counters, the workers that ran, and the
    // phase timings up to the refusal. Without it a refusal reports nothing, and a breaker that
    // stopped at the cap is indistinguishable from one that stopped long past it.
    FragmentScoringObservation observed;
    bool observed_valid = false;
    // True when no beam ran at all: K = 0, no phase audit and no positive development prefix, so the
    // candidate set is exactly the panel and nothing about it is conditional on a search. Reported as
    // candidate_scope PANEL_ONLY, which is neither COMPLETE nor CONDITIONAL -- both of those describe
    // the outcome of a search that did not happen.
    bool panel_only = false;
    // The tie-triggered full-beam audit. Each of these is a separate fact and they are reported
    // separately: conflating "the search was exhaustive" with "the prefix was truncated" with "the full
    // beam was exactly scored" is what made a K-truncated prefix look like grounds for consensus.
    //
    // primary_prefix_k          -- K the primary ranking was truncated to
    // full_beam_candidates      -- off-panel mosaics the beam retained, before truncation
    // primary_prefix_truncated  -- the primary ranking saw fewer candidates than the beam retained
    // full_beam_audit_attempted -- the primary tied, the search was exhaustive, and candidates sat
    //                              outside K, so the complete retained beam was scored exactly
    // candidates_exactly_audited-- candidates in the audited set, panel baseline included
    // audit_completed           -- the audit ranked to completion (no resource refusal)
    // final_candidate_scope     -- the scope of the set that actually DETERMINED the call
    std::size_t primary_prefix_k = 0;
    std::size_t full_beam_candidates = 0;
    bool primary_prefix_truncated = false;
    bool full_beam_audit_attempted = false;
    std::size_t candidates_exactly_audited = 0;
    bool audit_completed = false;
    std::string final_candidate_scope = "NA";
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
    // Block consensus only, and zero on every other path. When no whole-locus winner can be certified,
    // these count the blocks on which every tied candidate agreed -- which were called -- and those on
    // which they disagreed, which were no-called individually. A reader must be able to tell a consensus
    // run that called ten of eleven blocks from one that called two, without parsing the rows.
    std::size_t consensus_blocks_called = 0;
    std::size_t consensus_blocks_no_called = 0;
};

// What the panel-anchored adaptive completion did, for the status. Every fact separately: the anchor,
// the search and why it stopped, the pool, and the decision made over it.
struct RepairReport {
    const GenotypeCommandOptions* options = nullptr;
    std::uint64_t neighbour_bound = 0;
    std::string union_status = "NA";           // the union ranking before repair
    std::size_t union_winner_set = 0;
    double union_certified_gap = 0.0;
    std::string anchor_candidate_id = "NA";
    std::string anchor_union_rank = "NA";
    double anchor_score = 0.0;
    std::string anchor_hap1 = "NA", anchor_hap2 = "NA";
    bool searched = false;
    std::size_t rounds_run = 0;
    std::size_t moves = 0;
    std::string stop_reason = "NA";
    bool converged = false;
    double final_incumbent_score = 0.0;
    std::string final_incumbent_hap1 = "NA", final_incumbent_hap2 = "NA";
    std::size_t pool_union_candidates = 0;
    std::uint64_t pool_repair_scored = 0;
    std::size_t pool_identities = 0;
    std::size_t pool_repair_only = 0;
    std::uint64_t pool_score_mismatches = 0;
    std::size_t pool_sequence_candidates = 0;
    std::string decision = "NA";
    std::size_t top_set_size = 0;
    std::string disagreeing_blocks = "none";
    std::string genotype_separation_gap = "NA";
    std::string genotype_separation_certified = "NA";
    std::string best_other_genotype_score = "NA";
    std::string same_genotype_score_gap = "NA";
    std::string same_genotype_fragment_gap = "NA";
    std::string same_genotype_prior_gap = "NA";
    std::string same_genotype_fragment_identical = "NA";
    std::string haplotype_representation = "none";
    // strict_call, genotype_equivalent_call or no_call, and why a no-call is one. Reporting only.
    std::string call_class = "no_call";
    std::string no_call_reason = "NA";
    // The union's exact ranking top, before any repair: a labelled diagnostic fallback for a run that
    // emits no genotype, never a call.
    std::string union_top_hap1 = "NA", union_top_hap2 = "NA";
    double union_top_score = 0.0;
    FragmentScoringCounters cumulative;
    std::uint64_t footprint_at_repair_start = 0;
    std::uint64_t origin_cache_bytes_end = 0;
    std::size_t enumeration_workers = 0;
    double seconds = 0.0;
    double seconds_scoring = 0.0;
};

void write_repair_report(std::ostream& out, const RepairReport& r) {
    const GenotypeCommandOptions& o = *r.options;
    const auto margin = [](std::uint64_t cap, std::uint64_t used) {
        return std::to_string(static_cast<std::int64_t>(cap) - static_cast<std::int64_t>(used));
    };
    out.unsetf(std::ios::floatfield);
    out << std::setprecision(17);
    out << "algorithm\tpanel_anchored_adaptive_completion\n";
    out << "repair_max_rounds\t" << o.repair_rounds << "\nrepair_local_pairs\t" << o.repair_local_pairs
        << "\nrepair_origin_cache_bytes\t" << o.repair_origin_cache_bytes << "\nrepair_origin_cache_policy\tplanned"
        << "\nrepair_max_neighbours\t" << o.repair_max_neighbours << "\nrepair_neighbour_bound\t" << r.neighbour_bound
        << "\nrepair_max_physical_origins\t" << o.repair_max_physical_origins << "\nrepair_max_state_classes\t"
        << o.repair_max_state_classes << "\nrepair_max_piece_hits\t" << o.repair_max_piece_hits
        << "\nrepair_max_latent_pairs\t" << o.repair_max_latent_pairs << '\n';
    out << "call_class\t" << r.call_class << "\nno_call_reason\t" << r.no_call_reason << '\n';
    out << "union_ranking_status\t" << r.union_status << "\nunion_winner_set_size\t" << r.union_winner_set
        << "\nunion_certified_gap\t" << r.union_certified_gap << "\nunion_top_score\t" << r.union_top_score
        << "\nunion_top_hap1\t" << r.union_top_hap1 << "\nunion_top_hap2\t" << r.union_top_hap2 << '\n';
    out << "repair_anchor_candidate_id\t" << r.anchor_candidate_id << "\nrepair_anchor_union_rank\t"
        << r.anchor_union_rank << "\nrepair_anchor_score\t" << r.anchor_score << "\nrepair_anchor_hap1\t"
        << r.anchor_hap1 << "\nrepair_anchor_hap2\t" << r.anchor_hap2 << '\n';
    out << "repair_searched\t" << (r.searched ? "true" : "false") << "\nrepair_rounds_run\t" << r.rounds_run
        << "\nrepair_moves\t" << r.moves << "\nrepair_stop_reason\t" << r.stop_reason << "\nrepair_converged\t"
        << (r.converged ? "true" : "false") << "\nrepair_final_incumbent_score\t" << r.final_incumbent_score
        << "\nrepair_final_incumbent_hap1\t" << r.final_incumbent_hap1 << "\nrepair_final_incumbent_hap2\t"
        << r.final_incumbent_hap2 << '\n';
    out << "pool_union_candidates\t" << r.pool_union_candidates << "\npool_repair_scored\t" << r.pool_repair_scored
        << "\npool_identities\t" << r.pool_identities << "\npool_repair_only\t" << r.pool_repair_only
        << "\npool_score_mismatches\t" << r.pool_score_mismatches << "\npool_sequence_candidates\t"
        << r.pool_sequence_candidates << '\n';
    out << "final_decision\t" << r.decision << "\ntop_set_size\t" << r.top_set_size
        << "\ntop_set_disagreeing_blocks\t" << r.disagreeing_blocks << "\ngenotype_separation_gap\t"
        << r.genotype_separation_gap << "\ngenotype_separation_certified\t" << r.genotype_separation_certified
        << "\nbest_other_genotype_score\t" << r.best_other_genotype_score
        << "\nsame_genotype_alternative_score_gap\t" << r.same_genotype_score_gap
        << "\nsame_genotype_alternative_fragment_gap\t" << r.same_genotype_fragment_gap
        << "\nsame_genotype_alternative_prior_gap\t" << r.same_genotype_prior_gap
        << "\nsame_genotype_alternative_fragment_identical\t" << r.same_genotype_fragment_identical << '\n';
    // What the haplotype FASTA is, when there is one. prior_selected: the reads give no preference
    // among phase alternatives and the prior chose this representative; it is not a read-phased result.
    out << "haplotype_representation\t" << r.haplotype_representation << '\n';
    // The pool is a local search from one anchor: never an exhaustive search of genotype space.
    out << "search_globally_exhaustive\tfalse\nconditional_gq_calibration\tuncalibrated\n";
    out << "resources_footprint_at_repair_start_bytes\t" << r.footprint_at_repair_start << '\n';
    out << "repair_cumulative_physical_origins\t" << r.cumulative.physical_origins
        << "\nrepair_margin_physical_origins\t" << margin(o.repair_max_physical_origins, r.cumulative.physical_origins)
        << "\nrepair_cumulative_state_classes\t" << r.cumulative.statistical_state_classes
        << "\nrepair_margin_state_classes\t" << margin(o.repair_max_state_classes, r.cumulative.statistical_state_classes)
        << "\nrepair_cumulative_piece_hits\t" << r.cumulative.piece_hits << "\nrepair_margin_piece_hits\t"
        << margin(o.repair_max_piece_hits, r.cumulative.piece_hits) << "\nrepair_cumulative_latent_pairs\t"
        << r.cumulative.latent_alignment_pairs << "\nrepair_margin_latent_pairs\t"
        << margin(o.repair_max_latent_pairs, r.cumulative.latent_alignment_pairs)
        << "\nrepair_origin_cache_bytes_end\t" << r.origin_cache_bytes_end << "\nthreads_repair_enumeration_workers\t"
        << r.enumeration_workers << '\n';
    out << std::fixed << std::setprecision(3) << "seconds_repair\t" << r.seconds << "\nseconds_repair_scoring\t"
        << r.seconds_scoring << '\n';
    out.unsetf(std::ios::floatfield);
    out << std::setprecision(6);
}

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
    const CandidateSetReport* candidate_set = nullptr,
    const RepairReport* repair = nullptr) {
    std::ostringstream out;
    out << "key\tvalue\n";
    out << "status\t" << genotype_run_status_name(status) << '\n';
    if (!note.empty()) out << "note\t" << note << '\n';
    out << "model_fingerprint\t" << fingerprint << '\n';
    // What configured the run, and what it ran on, by content.
    out << "profile\t" << (options.profile_path.empty() ? "none" : options.profile_path) << '\n';
    out << "profile_name\t" << (options.profile_name.empty() ? "none" : options.profile_name) << '\n';
    out << "profile_sha256\t" << (options.profile_sha256.empty() ? "none" : options.profile_sha256) << '\n';
    {
        std::string overridden;
        for (const std::string& key : options.overridden_settings) overridden += (overridden.empty() ? "" : ",") + key;
        out << "settings_overridden_after_profile\t" << (overridden.empty() ? "none" : overridden) << '\n';
    }
    for (const auto& hash : g_input_hashes) out << "sha256_" << hash.first << '\t' << hash.second << '\n';
    out << "footprint_limit_bytes\t" << options.footprint_limit_bytes << '\n';
    out << "time_limit_seconds\t" << options.time_limit_seconds << '\n';
    {
        // Worker and process memory at the moment the status is written.
        const FootprintSample footprint = physical_footprint();
        out << "threads_peak_active_workers\t" << peak_active_workers() << '\n';
        out << "resources_phys_footprint_peak_bytes\t" << footprint.peak << '\n';
        out << "resources_rss_peak_bytes\t" << static_cast<std::uint64_t>(peak_rss_gib() * 1073741824.0) << '\n';
    }
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
    // Named per phase rather than as one "threads_used": a run must never claim parallelism in a phase
    // that ran serially. Placement is still serial and says so.
    out << "threads_used_combination\t"
        << (candidate_set != nullptr ? candidate_set->threads_used_combination : 1) << '\n';
    out << "threads_used_placement\t"
        << (candidate_set != nullptr ? candidate_set->threads_used_placement : 1) << '\n';
    out << "max_beam_extensions\t" << options.max_beam_extensions << '\n';
    out << "max_piece_hits\t" << options.max_piece_hits << '\n';
    out << "max_proposed_starts_per_haplotype\t" << options.max_proposed_starts_per_haplotype << '\n';
    out << "max_physical_origins\t" << options.max_physical_origins << '\n';
    // Declared here so an overridden cap or an armed cache is part of what the run is, not a flag
    // someone has to remember was passed.
    out << "emission_cache_bytes\t" << options.emission_cache_bytes << '\n';
    if (options.emission_cache_bytes > 0) {
        // What the cache actually did, as distinct from the logical work counters below, which stay
        // identical whether or not it was armed. `payload_bytes` is what the declared budget is
        // spent against and is the figure to size a larger run from.
        const EmissionCacheCounters cache = last_emission_cache_counters();
        out << "emission_cache_dp_requests\t" << cache.dp_cache_requests << '\n';
        out << "emission_cache_dp_hits\t" << cache.dp_cache_hits << '\n';
        out << "emission_cache_dp_misses\t" << cache.dp_cache_misses << '\n';
        out << "emission_cache_entries\t" << cache.entries << '\n';
        out << "emission_cache_payload_bytes\t" << cache.payload_bytes << '\n';
        out << "emission_cache_executed_prefilters\t" << cache.executed_prefilters << '\n';
        out << "emission_cache_executed_dps\t" << cache.executed_dps << '\n';
        out << "emission_cache_entries_not_retained\t" << cache.entries_not_retained << '\n';
    }
    // Explicit memory accounts of the scoring pass. Counted from the containers themselves, so a
    // question about where RSS went is answered by an account rather than by reasoning backwards.
    // Peaks are per worker and folded by max: no two workers hold one haplotype's placements at
    // once, so summing them would describe a moment that never existed.
    {
        const ScorerMemoryAccounts memory = last_scorer_memory_accounts();
        out << "scorer_mass_matrix_bytes\t" << memory.mass_matrix_bytes << '\n';
        out << "scorer_haplotype_sequence_bytes\t" << memory.haplotype_sequence_bytes << '\n';
        out << "scorer_form_sequence_bytes\t" << memory.form_sequence_bytes << '\n';
        out << "scorer_piece_index_bytes\t" << memory.piece_index_bytes << '\n';
        out << "scorer_emission_cache_table_bytes\t" << memory.emission_cache_table_bytes << '\n';
        out << "scorer_emission_cache_arena_bytes\t" << memory.emission_cache_arena_bytes << '\n';
        out << "scorer_placement_rows_peak\t" << memory.placement_rows_peak << '\n';
        out << "scorer_placement_capacity_bytes_peak\t" << memory.placement_capacity_bytes_peak
            << '\n';
        out << "scorer_proposed_starts_peak\t" << memory.proposed_starts_peak << '\n';
        out << "scorer_proposed_starts_capacity_bytes_peak\t"
            << memory.proposed_starts_capacity_bytes_peak << '\n';
    }
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
        // A curated panel must announce itself. Without these a run on a curated panel and a run on the
        // whole panel differ only in the fingerprint, which says they differ but not how.
        out << "panel_curation_manifest\t"
            << (options.curation_manifest_path.empty() ? "none" : options.curation_manifest_path)
            << '\n';
        out << "templates_curated_out\t" << model_report->curated_exclusions.size() << '\n';
        out << "curation_digest\t"
            << (model_report->curation_digest.empty() ? "none" : model_report->curation_digest) << '\n';
        for (const LocusInputReport::CuratedExclusion& excluded : model_report->curated_exclusions) {
            out << "curated_exclusion\t" << excluded.template_name << '\t' << excluded.reason << '\n';
        }
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
        // The production resource contract, reported as predicted against observed so the prediction is
        // checkable from the run's own output rather than taken on trust.
        out << "proposal_payload_ceiling_bytes\t" << options.max_proposal_payload_bytes << '\n';
        out << "proposal_payload_predicted_bytes\t" << beam->allocation_predicted.total_bytes << '\n';
        out << "proposal_payload_observed_bytes\t" << beam->allocation_observed.total_bytes << '\n';
        out << "proposal_payload_predicted_rolling_bytes\t" << beam->allocation_predicted.rolling_vector_bytes
            << '\n';
        out << "proposal_payload_predicted_ancestry_bytes\t" << beam->allocation_predicted.ancestry_bytes
            << '\n';
        out << "proposal_payload_predicted_pair_bytes\t" << beam->allocation_predicted.pair_bytes << '\n';
        out << "proposal_payload_predicted_index_bytes\t" << beam->allocation_predicted.index_bytes << '\n';
        out << "proposal_payload_predicted_child_records\t" << beam->allocation_predicted.child_records
            << '\n';
        out << "proposal_payload_observed_child_records\t" << beam->allocation_observed.child_records
            << '\n';
        out << "mosaic_candidates\t" << beam->candidates.size() << '\n';
        out << "mosaic_candidates_found\t" << beam->unique_candidates_found << '\n';
        out << "proposal_states_pruned\t" << beam->states_pruned << '\n';
        out << "proposal_search_exhaustive\t" << (beam->search_exhaustive ? "true" : "false") << '\n';
        out << "mosaic_target_reached\t" << (beam->unique_target_reached ? "true" : "false") << '\n';
        // Any ranking over a set that omits candidates is conditional on that set.
        out << "candidate_set_conditional\t" << (beam->candidate_set_conditional ? "true" : "false")
            << '\n';
        out << "beam_extensions\t" << beam->extensions_evaluated << '\n';
    }
    if (candidate_set != nullptr) {
        // mosaic_candidates, mosaic_candidates_found and mosaic_target_reached above count off-panel mosaics;
        // the panel baseline is added to them whole.
        out << "panel_diplotypes\t" << candidate_set->panel_diplotypes << '\n';
        out << "panel_baseline_complete\t" << (candidate_set->panel_baseline_complete ? "true" : "false") << '\n';
        out << "beam_diplotypes_found\t" << candidate_set->beam_diplotypes_found << '\n';
        out << "unique_haplotypes\t" << candidate_set->unique_haplotypes << '\n';
        out << "distinct_haplotype_paths\t" << candidate_set->distinct_haplotype_paths << '\n';
        out << "haploid_prior_evaluations\t" << candidate_set->haploid_prior_evaluations << '\n';
        out << "diploid_prior_lookups\t" << candidate_set->diploid_prior_lookups << '\n';
        out << "candidate_fragment_evaluations\t" << candidate_set->candidate_fragment_evaluations << '\n';
        // The tie-triggered full-beam audit, reported as separate concepts on purpose. `mosaic_candidates`
        // and `proposal_search_exhaustive` above describe the PROPOSAL; these describe what was exactly
        // SCORED. candidate_set_conditional keeps its existing meaning and is not redefined here: a run
        // whose audit completed is no longer conditional on the prefix, and final_candidate_scope says so
        // without disturbing the key that describes the prefix.
        out << "primary_prefix_k\t" << candidate_set->primary_prefix_k << '\n';
        out << "full_beam_candidates\t" << candidate_set->full_beam_candidates << '\n';
        out << "primary_prefix_truncated\t"
            << (candidate_set->primary_prefix_truncated ? "true" : "false") << '\n';
        out << "full_beam_audit_attempted\t"
            << (candidate_set->full_beam_audit_attempted ? "true" : "false") << '\n';
        out << "candidates_exactly_audited\t" << candidate_set->candidates_exactly_audited << '\n';
        out << "audit_completed\t" << (candidate_set->audit_completed ? "true" : "false") << '\n';
        out << "final_candidate_scope\t" << candidate_set->final_candidate_scope << '\n';
        if (candidate_set->observed_valid) {
            // A refused run reports the work it performed, so the breaker can be held to a bound rather
            // than merely to a verdict -- and the workers and phase timings with it, so the refusal can
            // be diagnosed rather than only believed.
            out << "observed_piece_hits\t" << candidate_set->observed.counters.piece_hits << '\n';
            out << "observed_physical_origins\t" << candidate_set->observed.counters.physical_origins
                << '\n';
            out << "observed_threads_used_placement\t" << candidate_set->observed.threads_used_placement
                << '\n';
            out << "observed_threads_used_combination\t"
                << candidate_set->observed.threads_used_combination << '\n';
        }
        // Derived once, here, and used verbatim by every block row.
        out << "candidate_scope\t"
            << candidate_scope_name(candidate_set->panel_only,
                                    beam != nullptr && beam->candidate_set_conditional)
            << '\n';
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
            << (repair == nullptr && ranking->status == RescoreStatus::Complete && beam != nullptr &&
                        !beam->candidate_set_conditional
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
            out << "consensus_blocks_called\t" << summary->consensus_blocks_called << '\n';
            out << "consensus_blocks_no_called\t" << summary->consensus_blocks_no_called << '\n';
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
        // How the verifications split: the pair-HMM ran for one, the exact feasibility prefilter
        // proved the other impossible first. Reported because the rejected fraction is the number
        // that says whether the prefilter is earning its place on real data rather than on a
        // fixture built to contain unplaceable reads.
        out << "dp_evaluations\t" << counters.dp_evaluations << '\n';
        out << "prefilter_rejections\t" << counters.prefilter_rejections << '\n';
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
    out << "seconds_rescoring_fragment_preparation\t" << timings.rescoring_fragment_preparation << '\n';
    out << "seconds_rescoring_placement_and_mass\t" << timings.rescoring_placement_and_mass << '\n';
    out << "seconds_rescoring_diploid_combination\t" << timings.rescoring_diploid_combination << '\n';
    out << "seconds_candidate_grouping\t" << timings.candidate_grouping << '\n';
    out << "seconds_prior\t" << timings.prior << '\n';
    out << "seconds_ranking\t" << timings.ranking << '\n';
    out << "seconds_total\t" << timings.total << '\n';
    if (repair != nullptr) write_repair_report(out, *repair);
    commit_status(path, out.str());
}

} // namespace

const char* genotype_run_status_name(GenotypeRunStatus status) {
    switch (status) {
        case GenotypeRunStatus::Complete: return "COMPLETE";
        case GenotypeRunStatus::ModelIncomplete: return "MODEL_INCOMPLETE";
        case GenotypeRunStatus::CandidateIncomplete: return "CANDIDATE_INCOMPLETE";
        case GenotypeRunStatus::SpellingIncomplete: return "SPELLING_INCOMPLETE";
        case GenotypeRunStatus::RescoreIncomplete: return "RESCORE_INCOMPLETE";
        case GenotypeRunStatus::BlockConsensus: return "BLOCK_CONSENSUS";
        case GenotypeRunStatus::ConditionalGenotypeEquivalent: return "CONDITIONAL_GENOTYPE_EQUIVALENT";
        case GenotypeRunStatus::RepairIncomplete: return "REPAIR_INCOMPLETE";
        case GenotypeRunStatus::TimeLimit: return "TIME_LIMIT";
        case GenotypeRunStatus::ResourceRefusal: return "RESOURCE_REFUSAL";
    }
    return "MODEL_INCOMPLETE";
}

std::vector<BlockConsensusCall> block_consensus(
    const std::vector<std::vector<std::pair<std::size_t, std::size_t>>>& tied_genotypes) {
    if (tied_genotypes.empty()) return {};
    std::size_t blocks = tied_genotypes.front().size();
    for (const auto& row : tied_genotypes) blocks = std::min(blocks, row.size());
    std::vector<BlockConsensusCall> out(blocks);
    for (std::size_t b = 0; b < blocks; ++b) {
        const std::pair<std::size_t, std::size_t> here = tied_genotypes.front()[b];
        bool agreed = true;
        for (const auto& row : tied_genotypes) {
            if (row[b] != here) {
                agreed = false;
                break;
            }
        }
        out[b].agreed = agreed;
        // The pair is reported only when it is agreed; an unagreed block has no consensus genotype and
        // must not be given one that happens to belong to whichever candidate came first.
        if (agreed) {
            out[b].allele1 = here.first;
            out[b].allele2 = here.second;
        }
    }
    return out;
}

void apply_genotype_profile(const std::string& path, GenotypeCommandOptions& options) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("cannot open genotype profile: " + path);
    const auto whole_u64 = [&](const std::string& key, const std::string& text) {
        std::size_t used = 0;
        const unsigned long long value = std::stoull(text, &used);
        if (used != text.size() || text.empty() || text[0] == '-') {
            throw std::runtime_error("profile " + path + ": '" + key + "' is not a non-negative integer: " + text);
        }
        return static_cast<std::uint64_t>(value);
    };
    const auto whole_double = [&](const std::string& key, const std::string& text) {
        std::size_t used = 0;
        const double value = std::stod(text, &used);
        if (used != text.size() || !std::isfinite(value) || value < 0.0) {
            throw std::runtime_error("profile " + path + ": '" + key + "' is not a finite non-negative number: " + text);
        }
        return value;
    };
    std::set<std::string> seen;
    std::string line;
    while (std::getline(file, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        const std::size_t tab = line.find('\t');
        if (tab == std::string::npos) throw std::runtime_error("profile " + path + ": expected key<TAB>value: " + line);
        const std::string key = line.substr(0, tab);
        const std::string value = line.substr(tab + 1);
        if (!seen.insert(key).second) throw std::runtime_error("profile " + path + ": '" + key + "' is set twice");
        const auto u = [&]() { return whole_u64(key, value); };
        if (key == "profile") options.profile_name = value;
        else if (key == "proposal_state_limit") options.proposal_state_limit = u();
        else if (key == "max_mosaic_candidates") options.max_mosaic_candidates = u();
        else if (key == "phase_audit_candidates") options.phase_audit_candidates = u();
        else if (key == "threads") options.threads = u();
        else if (key == "kmer_size") options.kmer_size = u();
        else if (key == "syncmer_s") options.syncmer_s = u();
        else if (key == "min_variant_bp") options.min_variant_bp = u();
        else if (key == "switch_probability") options.switch_probability = whole_double(key, value);
        else if (key == "mismatch_band") options.mismatch_band = u();
        else if (key == "max_piece_hits") options.max_piece_hits = u();
        else if (key == "max_physical_origins") options.max_physical_origins = u();
        else if (key == "max_state_classes") options.max_state_classes = u();
        else if (key == "max_mass_cells") options.max_mass_cells = u();
        else if (key == "emission_cache_bytes") options.emission_cache_bytes = u();
        else if (key == "repair_rounds") options.repair_rounds = u();
        else if (key == "repair_local_pairs") options.repair_local_pairs = u();
        else if (key == "repair_origin_cache_bytes") options.repair_origin_cache_bytes = u();
        else if (key == "repair_max_neighbours") options.repair_max_neighbours = u();
        else if (key == "repair_max_physical_origins") options.repair_max_physical_origins = u();
        else if (key == "repair_max_state_classes") options.repair_max_state_classes = u();
        else if (key == "repair_max_piece_hits") options.repair_max_piece_hits = u();
        else if (key == "repair_max_latent_pairs") options.repair_max_latent_pairs = u();
        else if (key == "footprint_limit_bytes") options.footprint_limit_bytes = u();
        else if (key == "time_limit_seconds") options.time_limit_seconds = whole_double(key, value);
        else throw std::runtime_error("profile " + path + ": unknown setting '" + key + "'");
    }
    if (seen.count("profile") == 0) throw std::runtime_error("profile " + path + " does not name itself ('profile')");
    options.profile_path = path;
    options.profile_sha256 = sha256_file(path);
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
    for (const char* suffix : kFinalOutputs) std::remove((out_prefix + suffix).c_str());
    for (const char* suffix : kPartialOutputs) std::remove((out_prefix + suffix).c_str());
    std::remove((out_prefix + ".mosaic_status.tsv").c_str());
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
    if (options.threads < 1) {
        throw std::runtime_error("--threads must be at least 1");
    }

    // The adaptive completion is one algorithm over one candidate set: a development ladder of prefixes
    // or a phase audit beside it would each need their own search, and is refused rather than guessed.
    if (options.repair_rounds > 0) {
        if (!options.development_candidate_prefixes.empty()) {
            throw std::runtime_error("the repair search does not run over development candidate prefixes");
        }
        if (options.phase_audit_candidates != 0) {
            throw std::runtime_error("the repair search runs with phase audit 0; the pooled rerank replaces it");
        }
        if (options.repair_local_pairs < 1) throw std::runtime_error("repair_local_pairs must be at least 1");
    }

    // Provenance by content, recorded in every status this run writes. An input that cannot be read is
    // left to the stage that reads it, which fails with its own message.
    reset_peak_active_workers();
    g_input_hashes.clear();
    {
        const auto hash = [&](const std::string& name, const std::string& file) {
            if (file.empty()) return;
            try {
                g_input_hashes.emplace_back(name, sha256_file(file));
            } catch (const std::exception&) {
                g_input_hashes.emplace_back(name, "unreadable");
            }
        };
        hash("gfa", options.gfa_path);
        hash("bubbles", options.bubbles_path);
        hash("reads1", options.reads1_path);
        hash("reads2", options.reads2_path);
        hash("library", options.library_path);
        hash("panel_curation", options.curation_manifest_path);
    }
    // The staged sorted graph is an intermediate, removed however the run ends.
    struct StagedGraph {
        std::string path;
        ~StagedGraph() {
            std::error_code ignored;
            std::filesystem::remove(path, ignored);
        }
    } staged_graph{options.out_prefix + ".sorted.gfa"};

    // The run watch: RUNNING status and partial rounds for a repair run, and the watchdog over the
    // declared time and physical-footprint limits. Disarmed on every exit.
    std::unique_ptr<RunWatch> watch;
    std::thread watchdog;
    if (options.repair_rounds > 0 || options.footprint_limit_bytes != 0 || options.time_limit_seconds > 0.0) {
        watch = std::make_unique<RunWatch>();
        watch->out_prefix = options.out_prefix;
        watch->profile = options.profile_path.empty() ? "none" : options.profile_path;
        watch->profile_sha256 = options.profile_sha256.empty() ? "none" : options.profile_sha256;
        watch->started = started;
        watch->footprint_limit = options.footprint_limit_bytes;
        watch->time_limit = options.time_limit_seconds;
        g_watch = watch.get();
        write_running_status("model");
        if (options.footprint_limit_bytes != 0 || options.time_limit_seconds > 0.0) {
            watchdog = std::thread([limit = options.footprint_limit_bytes, seconds = options.time_limit_seconds,
                                    started]() {
                for (;;) {
                    {
                        std::lock_guard<std::mutex> lock(g_watch->mutex);
                        if (g_watch->finished) return;
                    }
                    if (seconds > 0.0 && seconds_since(started) > seconds) {
                        std::ostringstream why;
                        why << "wall " << seconds_since(started) << " s exceeded the declared " << seconds << " s";
                        stop_run("TIME_LIMIT", why.str(), 124);
                    }
                    if (limit != 0) {
                        const std::uint64_t now = physical_footprint().current;
                        if (now > limit) {
                            std::ostringstream why;
                            why << "physical footprint " << now << " bytes exceeded the declared " << limit;
                            stop_run("RESOURCE_REFUSAL", why.str(), 125);
                        }
                    }
                    std::this_thread::sleep_for(std::chrono::milliseconds(100));
                }
            });
        }
    }
    struct WatchDisarm {
        std::thread& thread;
        bool armed;
        ~WatchDisarm() {
            if (!armed) return;
            {
                std::lock_guard<std::mutex> lock(g_watch->mutex);
                // An exception that escapes the run before any terminal status would otherwise leave
                // RUNNING behind, which a reader could take for a run still in progress.
                if (!g_watch->finished && std::uncaught_exceptions() > 0) {
                    try {
                        atomic_write(g_watch->out_prefix + ".mosaic_status.tsv",
                                     minimal_status(*g_watch, "FAILED",
                                                    "an exception ended the run before a final status"));
                    } catch (...) {
                    }
                }
                g_watch->finished = true;
            }
            if (thread.joinable()) thread.join();
            g_watch = nullptr;
        }
    } watch_disarm{watchdog, static_cast<bool>(watch)};

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
        if (g_watch != nullptr) {
            // Under the watch lock: a watchdog stop either sees the finished file and removes it, or
            // stops the process before it is begun. It can never leave one half-written behind.
            std::lock_guard<std::mutex> lock(g_watch->mutex);
            write_gfa_model(staged, gfa);
        } else {
            write_gfa_model(staged, gfa);
        }
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
        if (!options.curation_manifest_path.empty()) {
            locus_options.excluded_templates = read_curation_manifest(options.curation_manifest_path);
        }
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
    if (g_watch != nullptr) {
        std::lock_guard<std::mutex> lock(g_watch->mutex);
        g_watch->fingerprint = fingerprint;
    }
    write_running_status("markers");

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
    // selection units. The beam now supplies only recombinant mosaics beyond the panel, and K counts
    // those; the grouped proposer's own unit is the allele-prefix pair.
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
    write_running_status("proposal");
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
        beam_parameters.max_proposal_payload_bytes = options.max_proposal_payload_bytes;
        try {
            // The grouped-prefix proposer. It carries every template attribution of a prefix exactly
            // rather than approximating a prefix by its best one, and it is cheaper doing so. The
            // state-based propose_mosaic_diplotypes is retained as the fixture oracle it is checked
            // against, and is no longer called here.
            //
            // The three storage accounts are left at their defaults on purpose. Capping them from
            // grouped_allocation_bounds with the same inputs the proposer uses would gate a prediction
            // against itself and always admit; the accounts have force in the fixture, where the
            // MEASURED use is compared against the prediction. Work remains governed by
            // max_beam_extensions, and the rolling storage is small -- about 13 MB at C4 width 1,024.
            beam = propose_mosaic_diplotypes_grouped(model.panel(), emissions, switch_probability,
                                                     beam_parameters);
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
    // The tie-triggered full-beam audit scores the complete retained beam, and it can only re-rank
    // candidates this shared pass already scored: rank_set resolves every candidate through slot_of_path
    // and scored_index, both built from scored_union below, so a candidate outside it would throw rather
    // than be audited. The set is therefore widened here, where the single scoring pass is defined,
    // rather than by scoring a second time later.
    //
    // Only when an audit could actually be triggered: the search must be exhaustive, or the beam omitted
    // candidates it never proposed and no audit is offered whatever the ranking does. On a locus whose
    // search is pruned -- which is every large locus today -- this widens nothing and costs nothing.
    if (mosaics.search_exhaustive) {
        scored_mosaics = std::max(scored_mosaics, mosaics.candidates.size());
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

    write_running_status("union_scoring");
    mark = Clock::now();
    FragmentScoringOptions scoring_options;
    // Candidates partition into fixed contiguous ranges, each worker owning complete score entries, so
    // the thread count changes the schedule and nothing else. Gated byte-identical at 1, 2 and 4.
    scoring_options.threads = options.threads;
    scoring_options.mismatch_band = options.mismatch_band;
    scoring_options.minimum_piece_bp = options.minimum_piece_bp;
    scoring_options.max_piece_hits = options.max_piece_hits;
    scoring_options.max_proposed_starts_per_haplotype = options.max_proposed_starts_per_haplotype;
    scoring_options.max_physical_origins = options.max_physical_origins;
    scoring_options.max_state_classes = options.max_state_classes;
    scoring_options.max_mass_cells = options.max_mass_cells;
    scoring_options.emission_cache_bytes = options.emission_cache_bytes;
    FragmentScoringResult scored;
    std::string rescore_failure;
    if (any_ranked) {
        try {
            scored = score_diplotypes(interned, scored_pairs, fragments, library, scoring_options,
                                      &candidate_report.observed);
            timings.rescoring_fragment_preparation = scored.seconds_fragment_preparation;
            timings.rescoring_placement_and_mass = scored.seconds_placement_and_mass;
            timings.rescoring_diploid_combination = scored.seconds_diploid_combination;
            candidate_report.threads_used_combination = scored.threads_used_combination;
            candidate_report.threads_used_placement = scored.threads_used_placement;
            candidate_report.candidate_fragment_evaluations = scored.candidate_fragment_evaluations;
        } catch (const ResourceExhausted& error) {
            rescore_failure = error.what();
            candidate_report.observed_valid = true;
            // The phases that did run, up to the refusal. Without these a refused run reports its
            // whole duration under seconds_rescoring and 0.000 for each of the three parts, which is
            // precisely the breakdown a refusal needs.
            timings.rescoring_fragment_preparation =
                candidate_report.observed.seconds_fragment_preparation;
            timings.rescoring_placement_and_mass = candidate_report.observed.seconds_placement_and_mass;
            timings.rescoring_diploid_combination =
                candidate_report.observed.seconds_diploid_combination;
            // And the worker counts, so one status file does not report 1 here and 2 under observed_.
            // Two authoritative-looking answers to the same question are worse than the understatement
            // they replaced: a reader can take either, and nothing in the file says which is right.
            candidate_report.threads_used_placement = candidate_report.observed.threads_used_placement;
            candidate_report.threads_used_combination =
                candidate_report.observed.threads_used_combination;
        }
    }
    timings.rescoring = seconds_since(mark);

    // Exact-path prior cache, shared by every ranked set: the primary one, each development prefix and
    // the phase audit. The haploid prior is a function of one complete allele path alone, so it is
    // computed once per distinct path here rather than twice per diplotype inside every ranking -- which
    // on a 108,000-diplotype panel was 24-28 s per run, a third to a half of the whole thing.
    //
    // Keyed on the allele path, NEVER on the spelled sequence: two different paths can spell identical
    // bytes and must keep distinct path priors, which is exactly the distinction the sequence priors
    // below are for.
    std::map<std::vector<std::size_t>, double> haploid_prior_cache;
    std::size_t haploid_prior_evaluations = 0;
    std::size_t diploid_prior_lookups = 0;
    const auto cached_haploid_prior = [&](const MosaicPath& path) {
        const auto found = haploid_prior_cache.find(path.alleles);
        if (found != haploid_prior_cache.end()) return found->second;
        ++haploid_prior_evaluations;
        const double value = log_haploid_mosaic_prior(model.panel(), path, switch_probability);
        haploid_prior_cache.emplace(path.alleles, value);
        return value;
    };

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

        const auto candidates_mark = Clock::now();
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
        timings.candidate_grouping += seconds_since(candidates_mark);

        // ---- stage 5: exact prior, recomputed independently of the beam ------------------------------
        // The path prior is a manifest diagnostic; the ranking uses the sequence-pair prior, which decides
        // the log(2) multiplicity on bytes.
        //
        // The haploid prior is a function of one complete allele path alone, so it is cached by that path
        // and each diplotype is formed from two cached values. This is the same arithmetic
        // log_unordered_diplotype_prior performs -- the -inf short circuit and the log(2) multiplicity on
        // PATH identity, never on spelled sequence -- and it changes no model. The cache lives outside
        // this lambda, so the phase audit and every development prefix reuse it instead of repeating a
        // full set of evaluations.
        const auto prior_mark = Clock::now();
        std::vector<double> path_prior;
        path_prior.reserve(candidates.size());
        for (const MosaicCandidate& candidate : candidates) {
            ++diploid_prior_lookups;
            const double first = cached_haploid_prior(candidate.diplotype.first);
            const double second = cached_haploid_prior(candidate.diplotype.second);
            const double negative_infinity = -std::numeric_limits<double>::infinity();
            if (first == negative_infinity || second == negative_infinity) {
                path_prior.push_back(negative_infinity);
            } else {
                path_prior.push_back(
                    first + second +
                    (candidate.diplotype.first != candidate.diplotype.second ? std::log(2.0) : 0.0));
            }
        }
        timings.prior += seconds_since(prior_mark);
        candidate_report.distinct_haplotype_paths = haploid_prior_cache.size();
        candidate_report.haploid_prior_evaluations = haploid_prior_evaluations;
        candidate_report.diploid_prior_lookups = diploid_prior_lookups;
        const auto ranking_mark = Clock::now();
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

        timings.ranking += seconds_since(ranking_mark);
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
    // `finish` takes its StageTimings BY VALUE, so each development prefix reports its own total without
    // the prefixes after it. That copy is made at the call site, before rank_set runs, so the stages
    // rank_set accumulates -- candidates, prior, ranking -- are not in it, and every run reported them as
    // 0.000. They are run-level accumulations, like the shared work counters, so each status refreshes
    // them from this alias, which names the outer accumulator the parameter below shadows.
    StageTimings& run_timings = timings;
    // ---- panel-anchored adaptive completion ----------------------------------------------------------
    // The union just ranked is the starting evidence, not the decision. Its best exactly scored PANEL
    // diplotype anchors the exact coordinate repair (a deterministic search anchor, not a claim about
    // the sample); the union and every candidate the repair scored are pooled and exactly reranked, and
    // the call is made from the pool's numerical top set. Nothing here reads truth.
    const auto adaptive_completion = [&](const MosaicBeamResult& beam, const SetRanking& chosen,
                                         CandidateSetReport set_report, StageTimings timings) -> GenotypeRunStatus {
        const Clock::time_point repair_mark = Clock::now();
        const std::vector<MosaicCandidate>& candidates = chosen.candidates;
        const std::size_t panel_count = panel_candidates.size();
        const CertifiedRanking& union_ranking = chosen.ranking;
        const std::size_t block_count = model.blocks().size();
        const MosaicBeamResult* beam_view = panel_only ? nullptr : &beam;
        const std::string& prefix = options.out_prefix;
        constexpr std::size_t npos = std::numeric_limits<std::size_t>::max();
        RepairReport report;
        report.options = &options;
        report.neighbour_bound = repair_neighbour_bound(model, options.repair_local_pairs);
        report.union_status = union_ranking.status == RescoreStatus::Complete ? "SEPARATED" : "TIED";
        report.union_winner_set = union_ranking.winner_set.size();
        report.union_certified_gap = union_ranking.certified_gap;
        report.pool_union_candidates = candidates.size();
        {
            const MosaicDiplotype& top =
                candidates[chosen.sequence_candidates[union_ranking.candidates.front().candidate_id].representative].diplotype;
            report.union_top_hap1 = join_alleles(top.first.alleles);
            report.union_top_hap2 = join_alleles(top.second.alleles);
            report.union_top_score = union_ranking.candidates.front().score;
        }
        // A local search from one anchor never makes the candidate set complete.
        set_report.final_candidate_scope = "CONDITIONAL";
        const auto stop_with = [&](GenotypeRunStatus status, const std::string& note, const CertifiedRanking* ranking,
                                   const RankingSummary* summary) {
            timings.candidate_grouping = run_timings.candidate_grouping;
            timings.prior = run_timings.prior;
            timings.ranking = run_timings.ranking;
            set_report.distinct_haplotype_paths = candidate_report.distinct_haplotype_paths;
            set_report.haploid_prior_evaluations = candidate_report.haploid_prior_evaluations;
            set_report.diploid_prior_lookups = candidate_report.diploid_prior_lookups;
            report.seconds = seconds_since(repair_mark);
            report.call_class = status == GenotypeRunStatus::Complete                        ? "strict_call"
                                : status == GenotypeRunStatus::ConditionalGenotypeEquivalent ? "genotype_equivalent_call"
                                                                                             : "no_call";
            if (report.call_class != "no_call") report.no_call_reason = "none";
            timings.total = seconds_since(started);
            write_status(status_path, status, fingerprint, options, library, timings, &built.report, &panel, beam_view,
                         ranking, note, summary, &depth, &set_report, &report);
            if (!options.quiet) {
                std::cerr << "[genotype] " << genotype_run_status_name(status) << "  fingerprint=" << fingerprint
                          << "  " << note.substr(0, 160) << "\n";
            }
            return status;
        };
        // A partial file is written only while nothing has claimed the run, and never torn.
        const auto commit_partial = [&](const std::string& path, const std::string& text) {
            if (g_watch == nullptr) {
                atomic_write(path, text);
                return;
            }
            std::lock_guard<std::mutex> lock(g_watch->mutex);
            if (g_watch->finished) return;
            atomic_write(path, text);
        };

        // ---- 1. the anchor: the best exactly scored panel diplotype, by the union's exact ranking
        std::size_t anchor = npos;
        std::size_t anchor_rank = 0;
        for (std::size_t p = 0; p < union_ranking.candidates.size() && anchor == npos; ++p) {
            for (const std::size_t member : chosen.sequence_candidates[union_ranking.candidates[p].candidate_id].members) {
                if (member >= panel_count) continue;
                if (anchor == npos || chosen.path_prior[member] > chosen.path_prior[anchor]) anchor = member;
            }
            if (anchor != npos) anchor_rank = p;
        }
        if (anchor == npos) {
            report.no_call_reason = "no_panel_anchor";
            return stop_with(GenotypeRunStatus::RepairIncomplete,
                             "the union holds no exactly scored panel diplotype to anchor the repair", nullptr, nullptr);
        }
        const MosaicDiplotype& anchor_diplotype = candidates[anchor].diplotype;
        const double anchor_score = union_ranking.candidates[anchor_rank].score;
        report.anchor_candidate_id = std::to_string(candidates[anchor].candidate_id);
        report.anchor_union_rank = std::to_string(anchor_rank + 1);
        report.anchor_score = anchor_score;
        report.anchor_hap1 = join_alleles(anchor_diplotype.first.alleles);
        report.anchor_hap2 = join_alleles(anchor_diplotype.second.alleles);
        if (options.repair_max_neighbours != 0 && report.neighbour_bound > options.repair_max_neighbours) {
            report.no_call_reason = "neighbourhood_bound";
            return stop_with(GenotypeRunStatus::RepairIncomplete,
                             "one repair round's neighbourhood (" + std::to_string(report.neighbour_bound) +
                                 ") exceeds its declared bound " + std::to_string(options.repair_max_neighbours) +
                                 "; refused before any allocation",
                             nullptr, nullptr);
        }

        // ---- 2. the pool: the exactly scored union, then every candidate the repair scores
        struct PoolEntry {
            PoolCandidate candidate;
            bool in_union = false;
            std::size_t union_index = 0;
            std::string source;
            std::string scored_in;
        };
        std::vector<PoolEntry> pool;
        std::map<DiplotypeIdentity, std::size_t> pool_at;
        std::uint64_t mismatches = 0;
        const auto bitwise_equal = [](double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; };
        const auto add = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b, SequenceDigest da,
                             SequenceDigest db, const DiplotypeScore& score, double sequence_prior, double path_prior_value,
                             const std::string& source, const std::string& scored_in, bool in_union,
                             std::size_t union_index) {
            const DiplotypeIdentity key = diplotype_identity(a, b);
            const auto found = pool_at.find(key);
            if (found != pool_at.end()) {
                // One candidate scored twice -- by the union pass and by the repair -- must score
                // identically, bit for bit. Anything else fails the run closed below.
                const PoolCandidate& held = pool[found->second].candidate;
                if (!bitwise_equal(held.score.log_likelihood, score.log_likelihood) ||
                    !bitwise_equal(held.score.fragment_term, score.fragment_term) ||
                    !bitwise_equal(held.score.exposure_term, score.exposure_term) ||
                    !bitwise_equal(held.score.absolute_term_sum, score.absolute_term_sum) ||
                    held.score.fragments_with_mass != score.fragments_with_mass ||
                    !bitwise_equal(held.sequence_prior, sequence_prior)) {
                    ++mismatches;
                }
                return;
            }
            PoolEntry entry;
            entry.candidate.first = a;
            entry.candidate.second = b;
            entry.candidate.first_sequence = da;
            entry.candidate.second_sequence = db;
            entry.candidate.score = score;
            entry.candidate.sequence_prior = sequence_prior;
            entry.candidate.path_prior = path_prior_value;
            entry.in_union = in_union;
            entry.union_index = union_index;
            entry.source = source;
            entry.scored_in = scored_in;
            pool_at.emplace(key, pool.size());
            pool.push_back(std::move(entry));
        };
        {
            std::vector<SequenceDigest> digest(interned.sequences.size());
            for (std::size_t h = 0; h < interned.sequences.size(); ++h) digest[h] = sequence_digest(interned.sequences[h]);
            for (std::size_t i = 0; i < candidates.size(); ++i) {
                const MosaicCandidate& c = candidates[i];
                const std::size_t x = interned.index_of[slot_of_path.at(c.diplotype.first.alleles)];
                const std::size_t y = interned.index_of[slot_of_path.at(c.diplotype.second.alleles)];
                const auto key = sequence_pair(c);
                const double prior = haploid_sequence_prior[key.first] + haploid_sequence_prior[key.second] +
                                     (key.first != key.second ? std::log(2.0) : 0.0);
                add(c.diplotype.first.alleles, c.diplotype.second.alleles, digest[x], digest[y],
                    scored.scores[scored_index.at(key)], prior, chosen.path_prior[i], i < panel_count ? "panel" : "mosaic",
                    "union", true, i);
            }
        }
        // The union pass's mass matrix is not needed again; release it before the repair context exists.
        std::vector<double>().swap(scored.log_mass);
        std::vector<OriginSummary>().swap(scored.origin_summary);
        // What the process carries into the repair. The union pass leaves allocator-retained pages behind,
        // and how many varies run to run: NA19240 30x entered the repair at 1.55 GB in one run and 1.15 GB
        // in another, with the same scores. Recorded so a footprint refusal can be attributed.
        report.footprint_at_repair_start = physical_footprint().current;

        // ---- 3. the exact coordinate repair from the anchor
        FragmentScoringOptions repair_scoring = scoring_options;
        repair_scoring.max_physical_origins = options.repair_max_physical_origins;
        repair_scoring.max_state_classes = options.repair_max_state_classes;
        repair_scoring.max_piece_hits = options.repair_max_piece_hits;
        repair_scoring.max_latent_alignment_pairs = options.repair_max_latent_pairs;
        RepairSearchOptions search_options;
        search_options.local_pairs = options.repair_local_pairs;
        search_options.max_rounds = options.repair_rounds;
        search_options.mode = RepairScoringMode::Delta;
        search_options.origin_cache = true;
        search_options.planned_cache = true;
        search_options.origin_cache_bytes = options.repair_origin_cache_bytes;
        search_options.switch_probability = options.switch_probability;
        search_options.declared_neighbour_cap = options.repair_max_neighbours;
        search_options.memory_log = true;
        write_running_status("repair");
        std::unique_ptr<CoordinateRepairSearch> search;
        std::string refusal;
        bool reproduced = true;
        // One scored candidate into the pool. The diagnostic also keeps the round's items, so a terminal
        // tied round can be reconstructed from the scores the search already computed.
        struct RoundItem {
            std::size_t round = 0;
            std::vector<std::size_t> first, second;   // as generated, relative to the round's incumbent
            SequenceDigest first_sequence, second_sequence;
            DiplotypeScore score;
            double prior = 0.0;
        };
        std::vector<RoundItem> production_items;
        // DIAGNOSTIC digest memo (branch hooks only): keyed by allele path, which fixes the spelled bytes
        // under this model, so a released and recycled slot or a freed buffer can never alias an entry.
        std::map<std::vector<std::size_t>, SequenceDigest> digest_memo;
        std::uint64_t digest_memo_hits = 0, digest_memo_misses = 0, digest_memo_mismatches = 0;
        const auto digest_of = [&](const std::vector<std::size_t>& path, const std::string& bases, bool memo) {
            if (!memo) return sequence_digest(bases);
            const auto found = digest_memo.find(path);
            if (found == digest_memo.end()) {
                ++digest_memo_misses;
                const SequenceDigest fresh = sequence_digest(bases);
                digest_memo.emplace(path, fresh);
                return fresh;
            }
            ++digest_memo_hits;
            if (options.diagnostic_digest_memo_check) {
                const SequenceDigest fresh = sequence_digest(bases);
                if (fresh < found->second || found->second < fresh) ++digest_memo_mismatches;
            }
            return found->second;
        };
        const auto pool_scored = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b,
                                     const std::string& sa, const std::string& sb, const DiplotypeScore& score,
                                     double prior, const std::string& source, const std::string& scored_in,
                                     std::vector<RoundItem>* items, std::size_t round, bool memo = false) {
            const double pa = cached_haploid_prior(MosaicPath{a});
            const double pb = cached_haploid_prior(MosaicPath{b});
            const double negative_infinity = -std::numeric_limits<double>::infinity();
            const double path = pa == negative_infinity || pb == negative_infinity
                                    ? negative_infinity
                                    : pa + pb + (a != b ? std::log(2.0) : 0.0);
            const SequenceDigest da = digest_of(a, sa, memo);
            const SequenceDigest db = digest_of(b, sb, memo);
            add(a, b, da, db, score, prior, path, source, scored_in, false, 0);
            if (items != nullptr) items->push_back({round, a, b, da, db, score, prior});
        };
        try {
            search = std::make_unique<CoordinateRepairSearch>(model, emissions, fragments, library, repair_scoring,
                                                              search_options);
            search->start(anchor_diplotype.first.alleles, anchor_diplotype.second.alleles, anchor_score);
            RepairHooks hooks;
            hooks.on_scored = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b,
                                  const std::string& sa, const std::string& sb, const DiplotypeScore& score,
                                  double prior) {
                ++report.pool_repair_scored;
                const std::size_t round = search->rounds().size() + 1;
                pool_scored(a, b, sa, sb, score, prior, "repair", "round" + std::to_string(round),
                            options.diagnostic_tie_branching ? &production_items : nullptr, round);
            };
            hooks.on_round = [&](const RepairRound& round) {
                commit_partial(prefix + ".repair_rounds.partial.tsv", search->render_rounds(search->rounds()));
                if (search->memory_on()) commit_partial(prefix + ".repair_memory.partial.tsv", search->memory_text());
                if (g_watch != nullptr) {
                    std::lock_guard<std::mutex> lock(g_watch->mutex);
                    g_watch->repair_rounds_committed = search->rounds().size();
                }
                write_running_status("repair round " + std::to_string(round.index) + " committed");
                if (options.test_interrupt_after_repair_rounds != 0 && round.moved &&
                    search->rounds().size() >= options.test_interrupt_after_repair_rounds) {
                    stop_run("TIME_LIMIT",
                             "test hook: interrupted after " + std::to_string(search->rounds().size()) +
                                 " committed repair rounds (the time-limit path, triggered deterministically)",
                             124);
                }
            };
            reproduced = search->run(hooks);
        } catch (const ResourceExhausted& error) {
            refusal = error.what();
        } catch (const std::exception& error) {
            report.no_call_reason = "repair_failed";
            return stop_with(GenotypeRunStatus::RepairIncomplete,
                             std::string("the repair search failed and returned no genotype: ") + error.what(),
                             nullptr, nullptr);
        }
        if (search) {
            report.searched = true;
            report.rounds_run = search->rounds().size();
            for (const RepairRound& round : search->rounds()) report.moves += round.moved ? 1 : 0;
            report.stop_reason = refusal.empty() ? search->stop_reason() : "resource_refusal";
            report.converged = refusal.empty() && reproduced && search->converged();
            report.final_incumbent_score = search->incumbent_total();
            report.final_incumbent_hap1 = join_alleles(search->incumbent_first());
            report.final_incumbent_hap2 = join_alleles(search->incumbent_second());
            if (search->context() != nullptr) {
                report.cumulative = search->context()->cumulative_counters();
                report.origin_cache_bytes_end = search->context()->origin_cache_bytes();
                report.enumeration_workers = search->context()->parallel_enumeration_stats().workers;
            }
            report.seconds_scoring = search->cumulative_seconds();
        }
        report.pool_identities = pool.size();
        for (const PoolEntry& entry : pool) report.pool_repair_only += entry.in_union ? 0 : 1;
        report.pool_score_mismatches = mismatches;
        if (!refusal.empty()) {
            report.no_call_reason = "cumulative_cap";
            // Committed rounds stay as <prefix>.repair_rounds.partial.tsv: a record, never a call.
            return stop_with(GenotypeRunStatus::RescoreIncomplete,
                             "the repair search stopped at a declared cumulative cap and returned no genotype: " +
                                 refusal,
                             nullptr, nullptr);
        }
        // Every stop that returned normally leaves its whole trajectory under the final name.
        const auto write_trajectory = [&]() {
            atomic_write(prefix + ".repair_rounds.tsv", search->render_rounds(search->rounds()));
            if (search->memory_on()) atomic_write(prefix + ".repair_memory.tsv", search->memory_text());
            for (const char* suffix : kPartialOutputs) std::remove((prefix + suffix).c_str());
        };
        claim_run();
        if (!reproduced) {
            report.no_call_reason = "anchor_not_reproduced";
            write_trajectory();
            std::ostringstream why;
            why << std::setprecision(17) << "the repair anchor did not reproduce its union score in round "
                << search->failed_round() << " (expected " << search->incumbent_expected() << ", scored "
                << search->incumbent_scored() << "); no genotype";
            return stop_with(GenotypeRunStatus::RepairIncomplete, why.str(), nullptr, nullptr);
        }
        if (!report.converged) {
            // A round tie that excludes the incumbent is a tie among improving moves: the search cannot
            // choose between them and fails closed. Every other non-converged stop keeps its own name.
            const std::string& stop = report.stop_reason;
            report.no_call_reason =
                stop == "rescore_incomplete" || stop == "tie_at_best" ? "improving_move_tie" : stop;
            write_trajectory();
            const GenotypeRunStatus status = stop_with(GenotypeRunStatus::RepairIncomplete,
                             "the repair search stopped without converging (" + report.stop_reason +
                                 "); the trajectory is kept and no genotype is emitted",
                             nullptr, nullptr);
            // ---- DIAGNOSTIC ONLY (v2/post-hoc): improving-move-tie branching. Every production output
            // above is already written and final; this writes <prefix>.diag.* files and nothing else.
            if (options.diagnostic_tie_branching && report.no_call_reason == "improving_move_tie") {
                const Clock::time_point diag_mark = Clock::now();
                std::uint64_t diag_footprint_max = physical_footprint().current;
                const std::vector<std::vector<std::pair<std::size_t, std::size_t>>> local_pairs = search->local_pairs();
                struct Tie {
                    std::size_t id = 0;
                    long branch = -1;                  // -1: the production search
                    RoundItem incumbent;
                    std::vector<RoundItem> tied;       // canonical identity order
                };
                // The caller's own numerical rule, over the terminal round exactly as the round certified
                // it: one entry per scored sequence pair, the incumbent included; the top set is the
                // winner set. An improving tie is a top set that excludes the incumbent.
                const auto terminal_tie = [&](const std::vector<RoundItem>& items, std::size_t round) {
                    Tie tie;
                    std::vector<const RoundItem*> here;
                    for (const RoundItem& item : items) {
                        if (item.round == round) here.push_back(&item);
                    }
                    tie.incumbent = *here.front();
                    std::map<std::pair<SequenceDigest, SequenceDigest>, std::size_t> slot;
                    std::vector<std::vector<const RoundItem*>> groups;
                    std::vector<DiplotypeScore> scores;
                    std::vector<double> priors;
                    for (const RoundItem* item : here) {
                        SequenceDigest x = item->first_sequence;
                        SequenceDigest y = item->second_sequence;
                        if (y < x) std::swap(x, y);
                        const auto found = slot.emplace(std::make_pair(x, y), groups.size());
                        if (found.second) {
                            groups.emplace_back();
                            scores.push_back(item->score);
                            priors.push_back(item->prior);
                        }
                        groups[found.first->second].push_back(item);
                    }
                    const CertifiedRanking ranking = certify_scored_ranking(
                        scores, fragments.size(), scoring_options.mismatch_band, FragmentScoringCounters{},
                        CertificationOptions{}, priors);
                    const DiplotypeIdentity incumbent_key = diplotype_identity(tie.incumbent.first, tie.incumbent.second);
                    for (const std::size_t g : ranking.winner_set) {
                        for (const RoundItem* item : groups[g]) {
                            if (diplotype_identity(item->first, item->second) != incumbent_key) tie.tied.push_back(*item);
                        }
                    }
                    std::sort(tie.tied.begin(), tie.tied.end(), [](const RoundItem& x, const RoundItem& y) {
                        return diplotype_identity(x.first, x.second) < diplotype_identity(y.first, y.second);
                    });
                    return tie;
                };
                const auto genotype_text = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
                    std::string text;
                    for (std::size_t block = 0; block < block_count; ++block) {
                        const auto pair = std::minmax(a[block], b[block]);
                        text += (block == 0 ? "" : ",") + std::to_string(pair.first) + "/" + std::to_string(pair.second);
                    }
                    return text;
                };
                std::vector<Tie> ties;
                std::ostringstream tie_table;
                tie_table << std::setprecision(17)
                          << "tie\tbranch\tneighbour\thap1_alleles\thap2_alleles\tunordered_genotype\trepaired_block\t"
                             "allele_x\tallele_y\tassignment\tlocal_pair_rank\ttotal\tlog_likelihood\tfragment_term\t"
                             "exposure_term\tprior\tabsolute_term_sum\tincumbent_total\timprovement\n";
                const auto record_tie = [&](Tie tie) {
                    tie.id = ties.size();
                    const RoundItem& inc = tie.incumbent;
                    const double inc_total = inc.score.log_likelihood + inc.prior;
                    for (std::size_t k = 0; k < tie.tied.size(); ++k) {
                        const RoundItem& n = tie.tied[k];
                        std::size_t block = 0;
                        for (std::size_t b = 0; b < block_count; ++b) {
                            if (n.first[b] != inc.first[b] || n.second[b] != inc.second[b]) block = b;
                        }
                        const auto pair = std::minmax(n.first[block], n.second[block]);
                        const auto& local = local_pairs[block];
                        const auto at = std::find(local.begin(), local.end(), pair);
                        const double total = n.score.log_likelihood + n.prior;
                        const DiplotypeIdentity key = diplotype_identity(n.first, n.second);
                        tie_table << tie.id << '\t' << tie.branch << '\t' << k << '\t' << join_alleles(key.first) << '\t'
                                  << join_alleles(key.second) << '\t' << genotype_text(n.first, n.second) << '\t'
                                  << model.blocks()[block].block_id << '\t' << pair.first << '\t' << pair.second << '\t'
                                  << (n.first[block] == pair.first ? 0 : 1) << '\t'
                                  << (at == local.end() ? std::string("NA") : std::to_string(at - local.begin() + 1))
                                  << '\t' << total << '\t' << n.score.log_likelihood << '\t' << n.score.fragment_term
                                  << '\t' << n.score.exposure_term << '\t' << n.prior << '\t'
                                  << n.score.absolute_term_sum << '\t' << inc_total << '\t' << total - inc_total << '\n';
                    }
                    ties.push_back(tie);
                    return ties.back();
                };
                const std::string diag = prefix + ".diag";
                // ---- branches: breadth-first over tied neighbours, canonical order, nested ties included
                struct Branch {
                    std::size_t id = 0;
                    long parent_tie = 0;
                    RoundItem start;
                    std::vector<DiplotypeIdentity> history;   // incumbents before it: never revisited
                };
                std::vector<DiplotypeIdentity> production_history;
                production_history.push_back(diplotype_identity(anchor_diplotype.first.alleles, anchor_diplotype.second.alleles));
                for (const RepairRound& round : search->rounds()) {
                    if (round.moved) production_history.push_back(diplotype_identity(round.first, round.second));
                }
                const std::size_t production_last_round = search->rounds().back().index;
                std::vector<Branch> worklist;
                std::set<DiplotypeIdentity> started;
                const auto enqueue = [&](const Tie& tie, const std::vector<DiplotypeIdentity>& history) {
                    std::vector<RoundItem> order = tie.tied;
                    if (options.diagnostic_reverse_branch_order) std::reverse(order.begin(), order.end());
                    for (const RoundItem& n : order) {
                        if (!started.insert(diplotype_identity(n.first, n.second)).second) continue;   // one branch per start
                        Branch branch;
                        branch.parent_tie = static_cast<long>(tie.id);
                        branch.start = n;
                        branch.history = history;
                        worklist.push_back(branch);
                    }
                };
                const Tie first_tie = record_tie(terminal_tie(production_items, production_last_round));
                enqueue(first_tie, production_history);
                // MEASUREMENT ONLY: a search's rounds and the windows its context enumerated or reused.
                const auto dump_reuse = [&](const std::string& tag, CoordinateRepairSearch& traced) {
                    if (!options.diagnostic_reuse_trace) return;
                    atomic_write(diag + ".reuse." + tag + ".rounds.tsv", traced.render_rounds(traced.rounds()));
                    std::ostringstream demand;
                    demand << std::setprecision(17)
                           << "order\tkey_hash\twindow_bytes\tat_start\tat_end\trel_block_begin\trel_new_end\trequests\t"
                              "hits\tmisses\tgroup_reuses\tadmitted\tdeclined\tentry_bytes\trecords\torigins_constructed\t"
                              "seconds_enumeration\tseconds_reenumeration\tevictions\n";
                    if (traced.context() != nullptr) {
                        for (const OriginDemandRow& row : traced.context()->origin_demand()) {
                            demand << row.order << '\t' << row.key_hash << '\t' << row.window_bytes << '\t'
                                   << row.at_start << '\t' << row.at_end << '\t' << row.rel_block_begin << '\t'
                                   << row.rel_new_end << '\t' << row.requests << '\t' << row.hits << '\t' << row.misses
                                   << '\t' << row.group_reuses << '\t' << row.admitted << '\t' << row.declined << '\t'
                                   << row.entry_bytes << '\t' << row.records << '\t' << row.origins_constructed << '\t'
                                   << row.seconds_enumeration << '\t' << row.seconds_reenumeration << '\t'
                                   << row.evictions << '\n';
                        }
                    }
                    atomic_write(diag + ".reuse." + tag + ".demand.tsv", demand.str());
                };
                if (search) dump_reuse("production", *search);
                // Warm branches: every tied neighbour is seeded from its parent incumbent's census by the exact
                // block delta of a promoted repair. Parents stay pinned (alive, never mutated) for the whole
                // diagnostic; a seed is pinned until its own branch has finished.
                const bool warm = options.diagnostic_warm_branches;
                PreparedScoringContext* shared = warm && search ? search->context() : nullptr;
                std::map<DiplotypeIdentity, std::pair<std::size_t, std::size_t>> seed_slots;
                std::set<std::size_t> parent_pins;
                std::map<std::size_t, std::size_t> seed_pins;   // slot -> pending branches seeded on it
                std::uint64_t seeding_origins = 0;
                const auto seed_tie = [&](CoordinateRepairSearch& parent, const Tie& tie) {
                    parent_pins.insert(parent.incumbent_slot_first());
                    parent_pins.insert(parent.incumbent_slot_second());
                    const std::uint64_t before = shared->cumulative_counters().physical_origins;
                    for (const RoundItem& n : tie.tied) {
                        const DiplotypeIdentity key = diplotype_identity(n.first, n.second);
                        if (seed_slots.count(key) != 0) continue;
                        const auto slots = parent.seed_neighbour_by_delta(n.first, n.second);
                        seed_slots.emplace(key, slots);
                        ++seed_pins[slots.first];
                        if (slots.second != slots.first) ++seed_pins[slots.second];
                    }
                    seeding_origins += shared->cumulative_counters().physical_origins - before;
                };
                const auto pins_now = [&]() {
                    std::vector<std::size_t> pins(parent_pins.begin(), parent_pins.end());
                    for (const auto& seed : seed_pins) pins.push_back(seed.first);
                    return pins;
                };
                std::string warm_failure;
                if (warm) {
                    try {
                        if (shared == nullptr) throw std::runtime_error("no production scoring context to share");
                        seed_tie(*search, first_tie);
                    } catch (const std::exception& error) {
                        warm_failure = std::string("seeding failed: ") + error.what();
                    }
                } else {
                    // The emission cache belongs to one scoring context at a time: the production search is
                    // finished and its outputs written, so it is released before any branch builds its own.
                    search.reset();
                }
                const std::uint64_t shared_origins_at_start =
                    shared != nullptr ? shared->cumulative_counters().physical_origins : 0;
                (void)shared_origins_at_start;
                std::ostringstream branch_table;
                branch_table << std::setprecision(17)
                             << "branch\tparent_tie\tstart_hap1\tstart_hap2\tstart_total\trounds\tmoves\tstop\tconverged\t"
                                "outcome\tterminal_hap1\tterminal_hap2\tterminal_genotype\tterminal_total\tseconds\t"
                                "physical_origins\tphys_footprint_at_end\n";
                std::string incomplete;                // non-empty: DIAGNOSTIC_INCOMPLETE, and why
                std::map<DiplotypeIdentity, std::size_t> branch_of_start;    // for the symmetry check
                std::vector<std::vector<RoundItem>> branch_items;
                std::uint64_t total_origins = 0;
                std::size_t branches_run = 0;
                std::set<std::string> terminal_genotypes;
                const auto write_diag_status = [&](const std::string& state, const std::string& extra) {
                    std::ostringstream text;
                    text << std::setprecision(17) << "key\tvalue\nstate\t" << state
                         << "\nlabel\tdiagnostic_tie_branching_v2_posthoc_not_a_call\nproduction_status\t"
                         << genotype_run_status_name(status) << "\nproduction_no_call_reason\t" << report.no_call_reason
                         << "\nmax_branches\t" << options.diagnostic_max_branches << "\ntime_limit_seconds\t"
                         << options.diagnostic_time_limit_seconds << "\nfootprint_limit_bytes\t"
                         << options.diagnostic_footprint_limit_bytes << "\nmax_total_physical_origins\t"
                         << options.diagnostic_max_total_physical_origins << "\nbranch_order\t"
                         << (options.diagnostic_reverse_branch_order ? "reversed" : "canonical")
                         << "\nties\t" << ties.size() << "\nbranches_started\t" << started.size() << "\nbranches_run\t"
                         << branches_run << "\ntotal_physical_origins\t" << total_origins << "\nseconds\t"
                         << seconds_since(diag_mark) << "\nphys_footprint_max_bytes\t" << diag_footprint_max << '\n'
                         << extra;
                    atomic_write(diag + ".status.tsv", text.str());
                    atomic_write(diag + ".terminal_ties.tsv", tie_table.str());
                    atomic_write(diag + ".branches.tsv", branch_table.str());
                };
                total_origins += seeding_origins;   // the first tie's seeds (zero when cold)
                write_diag_status("RUNNING", "");
                struct DiagnosticStop : std::runtime_error {
                    using std::runtime_error::runtime_error;
                };
                for (std::size_t next = 0; next < worklist.size() && incomplete.empty(); ++next) {
                    if (branches_run >= options.diagnostic_max_branches) {
                        incomplete = "branch cap " + std::to_string(options.diagnostic_max_branches) + " binds (" +
                                     std::to_string(worklist.size()) + " branches required)";
                        break;
                    }
                    Branch branch = worklist[next];
                    branch.id = branches_run++;
                    const Clock::time_point branch_mark = Clock::now();
                    std::vector<RoundItem> items;
                    std::unique_ptr<CoordinateRepairSearch> branch_search;
                    std::string branch_failure;
                    bool branch_reproduced = true;
                    const std::uint64_t shared_before =
                        shared != nullptr ? shared->cumulative_counters().physical_origins : 0;
                    const DiplotypeIdentity start_identity = diplotype_identity(branch.start.first, branch.start.second);
                    try {
                        if (!warm_failure.empty()) throw std::runtime_error(warm_failure);
                        RepairSearchOptions branch_options = search_options;
                        branch_options.request_order_seed = options.diagnostic_branch_request_order_seed;
                        branch_options.external_context = warm;
                        if (warm) branch_options.memory_log = false;
                        branch_search = std::make_unique<CoordinateRepairSearch>(model, emissions, fragments, library,
                                                                                 repair_scoring, branch_options);
                        const double start_total = branch.start.score.log_likelihood + branch.start.prior;
                        if (warm) {
                            const auto seeded = seed_slots.at(start_identity);
                            branch_search->attach_context(shared);
                            branch_search->set_pinned(pins_now());
                            branch_search->start_from_slots(branch.start.first, branch.start.second, start_total,
                                                            seeded.first, seeded.second);
                        } else {
                            branch_search->start(branch.start.first, branch.start.second, start_total);
                        }
                        for (const DiplotypeIdentity& seen : branch.history) branch_search->add_visited(seen.first, seen.second);
                        RepairHooks hooks;
                        hooks.on_scored = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b,
                                              const std::string& sa, const std::string& sb, const DiplotypeScore& score,
                                              double prior) {
                            const std::size_t round = branch_search->rounds().size() + 1;
                            pool_scored(a, b, sa, sb, score, prior, "branch",
                                        "branch" + std::to_string(branch.id) + "_round" + std::to_string(round), &items,
                                        round, options.diagnostic_digest_memo);
                        };
                        hooks.on_round = [&](const RepairRound&) {
                            diag_footprint_max = std::max(diag_footprint_max, physical_footprint().current);
                            if (seconds_since(diag_mark) > options.diagnostic_time_limit_seconds) {
                                throw DiagnosticStop("diagnostic wall limit reached");
                            }
                            if (diag_footprint_max > options.diagnostic_footprint_limit_bytes) {
                                throw DiagnosticStop("diagnostic footprint limit reached");
                            }
                        };
                        branch_reproduced = branch_search->run(hooks);
                    } catch (const ResourceExhausted& error) {
                        branch_failure = std::string("work cap: ") + error.what();
                    } catch (const DiagnosticStop& error) {
                        branch_failure = error.what();
                    } catch (const std::exception& error) {
                        branch_failure = std::string("branch failed: ") + error.what();
                    }
                    if (branch_search) dump_reuse("branch" + std::to_string(branch.id), *branch_search);
                    // Executed work only: in the shared context, what this branch added to the cumulative count.
                    const std::uint64_t origins =
                        warm ? (shared != nullptr ? shared->cumulative_counters().physical_origins - shared_before : 0)
                             : (branch_search && branch_search->context()
                                    ? branch_search->context()->cumulative_counters().physical_origins
                                    : 0);
                    branch_of_start[diplotype_identity(branch.start.first, branch.start.second)] = branch.id;
                    branch_items.push_back(items);
                    total_origins += origins;
                    std::string outcome;
                    if (!branch_failure.empty()) {
                        outcome = "incomplete";
                        incomplete = "branch " + std::to_string(branch.id) + ": " + branch_failure;
                    } else if (!branch_reproduced) {
                        outcome = "incomplete";
                        incomplete = "branch " + std::to_string(branch.id) + ": its start did not reproduce";
                    } else if (branch_search->converged()) {
                        outcome = "converged";
                        terminal_genotypes.insert(genotype_text(branch_search->incumbent_first(), branch_search->incumbent_second()));
                    } else if (branch_search->stop_reason() == "rescore_incomplete" || branch_search->stop_reason() == "tie_at_best") {
                        outcome = "split";   // a nested improving tie: its tied neighbours become branches
                        std::vector<DiplotypeIdentity> history = branch.history;
                        history.push_back(diplotype_identity(branch.start.first, branch.start.second));
                        for (const RepairRound& round : branch_search->rounds()) {
                            if (round.moved) history.push_back(diplotype_identity(round.first, round.second));
                        }
                        Tie nested = terminal_tie(items, branch_search->rounds().back().index);
                        nested.branch = static_cast<long>(branch.id);
                        const Tie& recorded = record_tie(nested);
                        enqueue(recorded, history);
                        if (warm) {
                            try {
                                const std::uint64_t before = seeding_origins;
                                seed_tie(*branch_search, recorded);
                                total_origins += seeding_origins - before;
                            } catch (const std::exception& error) {
                                incomplete = "branch " + std::to_string(branch.id) + ": seeding failed: " + error.what();
                            }
                        }
                    } else {
                        outcome = "incomplete";
                        incomplete = "branch " + std::to_string(branch.id) + " stopped at " + branch_search->stop_reason();
                    }
                    if (warm) {
                        const auto seeded = seed_slots.find(start_identity);
                        if (seeded != seed_slots.end()) {
                            for (const std::size_t slot : {seeded->second.first, seeded->second.second}) {
                                const auto pin = seed_pins.find(slot);
                                if (pin != seed_pins.end() && --pin->second == 0) seed_pins.erase(pin);
                                if (seeded->second.second == seeded->second.first) break;
                            }
                        }
                    }
                    if (incomplete.empty() && total_origins > options.diagnostic_max_total_physical_origins) {
                        incomplete = "total physical-origin cap binds across branches";
                    }
                    diag_footprint_max = std::max(diag_footprint_max, physical_footprint().current);
                    const DiplotypeIdentity start_key = diplotype_identity(branch.start.first, branch.start.second);
                    branch_table << branch.id << '\t' << branch.parent_tie << '\t' << join_alleles(start_key.first) << '\t'
                                 << join_alleles(start_key.second) << '\t'
                                 << branch.start.score.log_likelihood + branch.start.prior << '\t'
                                 << (branch_search ? branch_search->rounds().size() : 0) << '\t';
                    std::size_t moves = 0;
                    if (branch_search) {
                        for (const RepairRound& round : branch_search->rounds()) moves += round.moved ? 1 : 0;
                    }
                    const DiplotypeIdentity end_key = branch_search ? diplotype_identity(branch_search->incumbent_first(),
                                                                                         branch_search->incumbent_second())
                                                                    : start_key;
                    branch_table << moves << '\t' << (branch_search ? branch_search->stop_reason() : "none") << '\t'
                                 << (branch_search && branch_search->converged() ? "true" : "false") << '\t' << outcome
                                 << '\t' << join_alleles(end_key.first) << '\t' << join_alleles(end_key.second) << '\t'
                                 << genotype_text(end_key.first, end_key.second) << '\t'
                                 << (branch_search ? branch_search->incumbent_total() : 0.0) << '\t'
                                 << seconds_since(branch_mark) << '\t' << origins << '\t' << physical_footprint().current
                                 << '\n';
                    write_diag_status("RUNNING", "");   // atomic partial after every branch
                }
                // ---- phase-switch symmetry of every two-way tie whose two branches both ran
                {
                    std::ostringstream symmetry;
                    symmetry << std::setprecision(17)
                             << "tie\tswitch_block\tbranch_a\tbranch_b\ttied_pair_bit_identical\tequivalent\t"
                                "first_break_round\tfirst_break\tnumerically_equivalent\tnumeric_first_break_round\t"
                                "numeric_first_break\tmax_difference_over_bound\tdepths\n";
                    for (const Tie& tie : ties) {
                        if (tie.tied.size() != 2) continue;
                        const RoundItem& x = tie.tied[0];
                        const RoundItem& y = tie.tied[1];
                        const auto ax = branch_of_start.find(diplotype_identity(x.first, x.second));
                        const auto ay = branch_of_start.find(diplotype_identity(y.first, y.second));
                        if (ax == branch_of_start.end() || ay == branch_of_start.end()) continue;
                        // The block where the pair's allele assignment differs; a local switch there maps one onto
                        // the other exactly when it is the only difference.
                        std::size_t block = block_count;
                        bool local_switch = true;
                        for (std::size_t b = 0; b < block_count; ++b) {
                            if (x.first[b] == y.first[b] && x.second[b] == y.second[b]) continue;
                            if (x.first[b] == y.second[b] && x.second[b] == y.first[b] && block == block_count) block = b;
                            else local_switch = false;
                        }
                        if (!local_switch || block == block_count) {
                            symmetry << tie.id << "\tNA\t" << ax->second << '\t' << ay->second
                                     << "\tNA\tfalse\t0\tthe tied pair is not a local phase switch at one block\tfalse\t0\t"
                                        "the tied pair is not a local phase switch at one block\tNA\tNA\n";
                            continue;
                        }
                        const auto to_trace = [](const std::vector<RoundItem>& items) {
                            std::vector<SymmetryTraceItem> out;
                            for (const RoundItem& item : items) out.push_back({item.round, item.first, item.second, item.score, item.prior});
                            return out;
                        };
                        const SymmetryReport report = check_phase_switch_symmetry(
                            to_trace(branch_items[ax->second]), to_trace(branch_items[ay->second]), block,
                            fragments.size());
                        const bool pair_identical =
                            std::memcmp(&x.score.log_likelihood, &y.score.log_likelihood, sizeof(double)) == 0 &&
                            std::memcmp(&x.score.fragment_term, &y.score.fragment_term, sizeof(double)) == 0 &&
                            std::memcmp(&x.score.exposure_term, &y.score.exposure_term, sizeof(double)) == 0 &&
                            std::memcmp(&x.prior, &y.prior, sizeof(double)) == 0;
                        std::string depths;
                        double over_bound = 0.0;
                        for (const SymmetryDepth& d : report.depths) {
                            depths += (depths.empty() ? "" : ";") + std::to_string(d.round) + ":" +
                                      std::to_string(d.a_candidates) + "/" + std::to_string(d.b_candidates) + "/" +
                                      std::to_string(d.matched) + "/" + std::to_string(d.bit_identical) + "/" +
                                      std::to_string(d.indistinguishable);
                            over_bound = std::max(over_bound, d.max_difference_over_bound);
                        }
                        symmetry << tie.id << '\t' << model.blocks()[block].block_id << '\t' << ax->second << '\t'
                                 << ay->second << '\t' << (pair_identical ? "true" : "false") << '\t'
                                 << (report.equivalent ? "true" : "false") << '\t' << report.first_break_round << '\t'
                                 << (report.first_break.empty() ? "none" : report.first_break) << '\t'
                                 << (report.numerically_equivalent ? "true" : "false") << '\t'
                                 << report.numeric_first_break_round << '\t'
                                 << (report.numeric_first_break.empty() ? "none" : report.numeric_first_break) << '\t'
                                 << over_bound << '\t' << depths << '\n';
                    }
                    atomic_write(diag + ".symmetry.tsv", symmetry.str());
                }
                // ---- the decision: the production rule over the pooled numerical top set
                std::ostringstream extra;
                extra << std::setprecision(17) << "distinct_branch_terminal_genotypes\t" << terminal_genotypes.size() << '\n'
                      << "branch_seeding\t" << (warm ? "delta_promoted_shared_context" : "cold_full_placement") << '\n'
                      << "seeding_physical_origins\t" << seeding_origins << '\n'
                      << "digest_memo\t" << (options.diagnostic_digest_memo ? "allele_path" : "off") << '\n'
                      << "digest_memo_hits\t" << digest_memo_hits << "\ndigest_memo_misses\t" << digest_memo_misses << '\n'
                      << "digest_memo_check\t" << (options.diagnostic_digest_memo_check ? "on" : "off") << '\n'
                      << "digest_memo_mismatches\t" << digest_memo_mismatches << '\n'
                      << "branch_request_order_seed\t" << options.diagnostic_branch_request_order_seed << '\n';
                if (options.diagnostic_digest_memo_check && digest_memo_mismatches != 0 && incomplete.empty()) {
                    incomplete = "the digest memo disagreed with a recomputed digest";
                }
                if (options.diagnostic_pool_dump) {
                    std::vector<std::pair<DiplotypeIdentity, const PoolEntry*>> ordered;
                    for (const PoolEntry& entry : pool) {
                        ordered.push_back({diplotype_identity(entry.candidate.first, entry.candidate.second), &entry});
                    }
                    std::sort(ordered.begin(), ordered.end(),
                              [](const auto& x, const auto& y) { return x.first < y.first; });
                    std::ostringstream dump;
                    dump << std::hexfloat << "hap1\thap2\tlog_likelihood\tfragment_term\texposure_term\tabsolute_term_sum\t"
                                             "fragments_with_mass\tsequence_prior\tpath_prior\n";
                    for (const auto& item : ordered) {
                        const PoolCandidate& c = item.second->candidate;
                        dump << join_alleles(item.first.first) << '\t' << join_alleles(item.first.second) << '\t'
                             << c.score.log_likelihood << '\t' << c.score.fragment_term << '\t' << c.score.exposure_term
                             << '\t' << c.score.absolute_term_sum << '\t' << c.score.fragments_with_mass << '\t'
                             << c.sequence_prior << '\t' << c.path_prior << '\n';
                    }
                    atomic_write(diag + ".pool.tsv", dump.str());
                }
                if (!incomplete.empty()) {
                    extra << "decision\tDIAGNOSTIC_INCOMPLETE\nincomplete_reason\t" << incomplete << '\n';
                    write_diag_status("INCOMPLETE", extra.str());
                } else {
                    std::vector<PoolCandidate> pooled;
                    pooled.reserve(pool.size());
                    for (const PoolEntry& entry : pool) pooled.push_back(entry.candidate);
                    const PooledDecision decided = decide_pooled_genotype(std::move(pooled), block_count, fragments.size(),
                                                                          scoring_options.mismatch_band);
                    std::set<std::string> top_genotypes;
                    std::ostringstream top;
                    top << std::setprecision(17) << "rank\tsequence_candidate\thap1_alleles\thap2_alleles\t"
                                                    "unordered_genotype\tscored_in\tscore\tfragment_term\texposure_term\tprior\t"
                                                    "numeric_error_bound\n";
                    std::size_t rank = 0;
                    const std::set<std::size_t> in_top(decided.top_groups.begin(), decided.top_groups.end());
                    for (const CertifiedCandidate& c : decided.ranking.candidates) {
                        ++rank;
                        if (in_top.count(c.candidate_id) == 0) continue;
                        for (const std::size_t m : decided.groups[c.candidate_id]) {
                            const PoolCandidate& pc = decided.candidates[m];
                            top_genotypes.insert(genotype_text(pc.first, pc.second));
                            top << rank << '\t' << c.candidate_id << '\t' << join_alleles(pc.first) << '\t'
                                << join_alleles(pc.second) << '\t' << genotype_text(pc.first, pc.second) << '\t'
                                << pool[pool_at.at(diplotype_identity(pc.first, pc.second))].scored_in << '\t' << c.score
                                << '\t' << c.fragment_term << '\t' << c.exposure_term << '\t' << c.prior << '\t'
                                << c.numeric_error_bound << '\n';
                        }
                    }
                    atomic_write(diag + ".pooled_top.tsv", top.str());
                    const std::string outcome =
                        mismatches + decided.score_mismatches != 0          ? "DIAGNOSTIC_INCOMPLETE"
                        : decided.decision == TopSetDecision::Unique             ? "DIAGNOSTIC_UNIQUE"
                        : decided.decision == TopSetDecision::GenotypeEquivalent ? "DIAGNOSTIC_GENOTYPE_EQUIVALENT"
                                                                                 : "DIAGNOSTIC_DIVERGENT";
                    const PoolCandidate& rep = decided.candidates[decided.representative[decided.ranking.candidates.front().candidate_id]];
                    extra << "decision\t" << outcome << "\npool_identities\t" << pool.size() << "\npool_score_mismatches\t"
                          << mismatches + decided.score_mismatches << "\npooled_top_size\t" << decided.top_groups.size()
                          << "\npooled_top_genotypes\t" << top_genotypes.size() << "\npooled_top_score\t"
                          << decided.ranking.candidates.front().score << "\nrepresentative_hap1\t"
                          << join_alleles(rep.first) << "\nrepresentative_hap2\t" << join_alleles(rep.second)
                          << "\nrepresentative_genotype\t" << genotype_text(rep.first, rep.second)
                          << "\nphase_support\t" << decided.phase_support << "\ngenotype_separation_gap\t"
                          << decided.genotype_separation_gap << '\n';
                    write_diag_status("COMPLETE", extra.str());
                }
            }
            return status;
        }
        if (mismatches != 0) {
            report.no_call_reason = "pool_score_mismatch";
            write_trajectory();
            return stop_with(GenotypeRunStatus::RepairIncomplete,
                             std::to_string(mismatches) +
                                 " candidates scored by both the union pass and the repair did not score "
                                 "identically; the pool is not exact and no genotype is emitted",
                             nullptr, nullptr);
        }

        // ---- 4. the exact rerank of the pool, and the decision over its numerical top set
        const Clock::time_point decision_mark = Clock::now();
        std::vector<PoolCandidate> pooled;
        pooled.reserve(pool.size());
        for (const PoolEntry& entry : pool) pooled.push_back(entry.candidate);
        PooledDecision decided = decide_pooled_genotype(std::move(pooled), block_count, fragments.size(),
                                                        scoring_options.mismatch_band);
        run_timings.ranking += seconds_since(decision_mark);
        decided.ranking.counters = scored.counters;   // the union pass; the repair's are reported apart
        report.pool_sequence_candidates = decided.groups.size();
        if (decided.score_mismatches != 0) {
            report.no_call_reason = "sequence_score_mismatch";
            write_trajectory();
            return stop_with(GenotypeRunStatus::RepairIncomplete,
                             std::to_string(decided.score_mismatches) +
                                 " candidates spelling one sequence pair did not score identically; no genotype",
                             nullptr, nullptr);
        }
        const CertifiedCandidate& best = decided.ranking.candidates.front();
        std::ostringstream number;
        number << std::setprecision(17);
        const auto text = [&](double value) {
            number.str("");
            number << value;
            return number.str();
        };
        report.decision = top_set_decision_name(decided.decision);
        report.top_set_size = decided.top_groups.size();
        if (!decided.disagreeing_blocks.empty()) {
            report.disagreeing_blocks.clear();
            for (const std::size_t b : decided.disagreeing_blocks) {
                report.disagreeing_blocks += (report.disagreeing_blocks.empty() ? "" : ",") + model.blocks()[b].block_id;
            }
        }
        if (decided.has_other_genotype) {
            report.genotype_separation_gap = text(decided.genotype_separation_gap);
            report.genotype_separation_certified = decided.genotype_separation_certified ? "true" : "false";
            for (const CertifiedCandidate& c : decided.ranking.candidates) {
                if (c.candidate_id == decided.best_other_genotype_group) report.best_other_genotype_score = text(c.score);
            }
        } else {
            report.genotype_separation_gap = "inf";
            report.genotype_separation_certified = "true";
        }
        if (decided.has_same_genotype_alternative) {
            report.same_genotype_score_gap = text(decided.same_genotype_score_gap);
            report.same_genotype_fragment_gap = text(decided.same_genotype_fragment_gap);
            report.same_genotype_prior_gap = text(decided.same_genotype_prior_gap);
            report.same_genotype_fragment_identical = decided.same_genotype_fragment_identical ? "true" : "false";
        }

        RankingSummary summary;
        summary.sequence_candidates = decided.groups.size();
        summary.collapsed_beam_candidates = decided.candidates.size() - decided.groups.size();
        summary.phase_support = decided.phase_support;
        summary.phase_support_within_candidates = decided.phase_support;
        summary.fragment_contenders = decided.fragment_contenders;
        summary.fragment_term_gap = decided.fragment_term_gap;
        summary.projection_ambiguous = decided.projection_ambiguous[best.candidate_id];
        summary.block_genotype_read_supported = decided.block_genotype_read_supported;
        summary.block_genotype_posterior = decided.block_genotype_posterior;
        summary.block_genotype_quality = decided.block_genotype_quality;
        summary.block_numerically_separated = decided.block_numerically_separated;

        // Table identities: union candidates keep theirs; repair-only candidates follow, canonically.
        std::vector<std::size_t> table_id(decided.candidates.size());
        std::vector<const PoolEntry*> entry_of(decided.candidates.size());
        std::vector<std::size_t> group_of(decided.candidates.size());
        {
            std::size_t next = 0;
            for (const MosaicCandidate& c : candidates) next = std::max(next, c.candidate_id + 1);
            for (std::size_t i = 0; i < decided.candidates.size(); ++i) {
                const PoolCandidate& c = decided.candidates[i];
                entry_of[i] = &pool[pool_at.at(diplotype_identity(c.first, c.second))];
                table_id[i] = entry_of[i]->in_union ? candidates[entry_of[i]->union_index].candidate_id : next++;
            }
            for (std::size_t g = 0; g < decided.groups.size(); ++g) {
                for (const std::size_t m : decided.groups[g]) group_of[m] = g;
            }
        }
        const std::set<std::size_t> top(decided.top_groups.begin(), decided.top_groups.end());
        const auto carries_call = [&](std::size_t g) {
            for (const std::size_t m : decided.groups[g]) {
                for (std::size_t b = 0; b < block_count; ++b) {
                    if (std::minmax(decided.candidates[m].first[b], decided.candidates[m].second[b]) != decided.genotype[b]) {
                        return false;
                    }
                }
            }
            return true;
        };

        // ---- outputs. The run is claimed: the watchdog can no longer stop it half-written.
        {
            std::vector<std::size_t> order(decided.candidates.size());
            for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
            std::sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) { return table_id[x] < table_id[y]; });
            std::ofstream out = open_output(prefix + ".mosaic_candidates.tsv");
            out << std::setprecision(17)
                << "candidate_id\tunique_rank\thap1_alleles\thap2_alleles\tbest_proposal_score\t"
                   "retained_logsum_proposal_score\tproposal_emission_component\tproposal_transition_component\t"
                   "equivalent_retained_histories\tevery_allele_scored\texact_path_prior\tsequence_candidate_id\t"
                   "history_count_saturated\tsource\tscored_in\n";
            for (const std::size_t i : order) {
                const PoolCandidate& c = decided.candidates[i];
                const PoolEntry& entry = *entry_of[i];
                out << table_id[i] << '\t';
                if (entry.in_union && entry.union_index >= panel_count) {
                    const MosaicCandidate& m = candidates[entry.union_index];
                    out << m.unique_rank << '\t' << join_alleles(c.first) << '\t' << join_alleles(c.second) << '\t'
                        << m.best_proposal_score << '\t' << m.retained_logsum_proposal_score << '\t'
                        << m.proposal_emission_component << '\t' << m.proposal_transition_component << '\t'
                        << m.equivalent_retained_histories << '\t' << (m.every_allele_scored ? "true" : "false") << '\t'
                        << c.path_prior << '\t' << group_of[i] << '\t' << (m.history_count_saturated ? "true" : "false");
                } else {
                    out << "NA\t" << join_alleles(c.first) << '\t' << join_alleles(c.second) << "\tNA\tNA\tNA\tNA\tNA\t"
                        << (emissions.every_allele_scored() ? "true" : "false") << '\t' << c.path_prior << '\t'
                        << group_of[i] << "\tNA";
                }
                out << '\t' << entry.source << '\t' << entry.scored_in << '\n';
            }
        }
        {
            std::ofstream out = open_output(prefix + ".mosaic_scores.tsv");
            out << "sequence_candidate_id\trepresentative_beam_candidate\tequivalent_beam_candidates\t"
                   "score\tfragment_term\texposure_term\tsequence_pair_prior\tnumeric_error_bound\t"
                   "fragments_with_mass\tfragments_total\tcandidate_set_probability\tmismatch_band\t"
                   "projection_ambiguous\tin_winner_set\tcarries_called_genotype\n";
            out << std::setprecision(12);
            for (const CertifiedCandidate& entry : decided.ranking.candidates) {
                const std::size_t g = entry.candidate_id;
                out << g << '\t' << table_id[decided.representative[g]] << '\t' << decided.groups[g].size() << '\t'
                    << entry.score << '\t' << entry.fragment_term << '\t' << entry.exposure_term << '\t' << entry.prior
                    << '\t' << entry.numeric_error_bound << '\t' << entry.fragments_with_mass << '\t'
                    << fragments.size() << '\t' << entry.candidate_set_probability << '\t'
                    << decided.ranking.mismatch_band << '\t' << (decided.projection_ambiguous[g] ? "true" : "false")
                    << '\t' << (top.count(g) != 0 ? "true" : "false") << '\t' << (carries_call(g) ? "true" : "false")
                    << '\n';
            }
        }
        {
            // The complete numerical top set, every member of every tied sequence candidate.
            std::ofstream out = open_output(prefix + ".equivalence_set.tsv");
            out << "rank\tsequence_candidate_id\tcandidate_id\trepresentative\thap1_alleles\thap2_alleles\tsource\t"
                   "scored_in\tscore\tfragment_term\texposure_term\tsequence_pair_prior\tnumeric_error_bound\t"
                   "exact_path_prior\n";
            out << std::setprecision(17);
            std::size_t rank = 0;
            for (const CertifiedCandidate& entry : decided.ranking.candidates) {
                ++rank;
                if (top.count(entry.candidate_id) == 0) continue;
                for (const std::size_t m : decided.groups[entry.candidate_id]) {
                    const PoolCandidate& c = decided.candidates[m];
                    out << rank << '\t' << entry.candidate_id << '\t' << table_id[m] << '\t'
                        << (decided.representative[entry.candidate_id] == m ? "true" : "false") << '\t'
                        << join_alleles(c.first) << '\t' << join_alleles(c.second) << '\t' << entry_of[m]->source
                        << '\t' << entry_of[m]->scored_in << '\t' << entry.score << '\t' << entry.fragment_term << '\t'
                        << entry.exposure_term << '\t' << entry.prior << '\t' << entry.numeric_error_bound << '\t'
                        << c.path_prior << '\n';
                }
            }
        }
        write_trajectory();

        if (decided.decision == TopSetDecision::DistinctTie) {
            report.no_call_reason = "distinct_tie";
            return stop_with(GenotypeRunStatus::RescoreIncomplete,
                             "the pooled top set holds " + std::to_string(decided.top_groups.size()) +
                                 " candidates tied on total score that disagree at block(s) " + report.disagreeing_blocks +
                                 "; the equivalence set is reported and no genotype is called",
                             &decided.ranking, &summary);
        }
        const bool unique = decided.decision == TopSetDecision::Unique;
        const PoolCandidate& called = decided.candidates[decided.representative[best.candidate_id]];
        // Allele order: the representative's haplotypes for a unique call; the unordered pair otherwise.
        const auto alleles_at = [&](std::size_t b) {
            return unique ? std::make_pair(called.first[b], called.second[b]) : decided.genotype[b];
        };
        const auto sources_of = [&](const PreparedGenotypeBlock& block, std::size_t id) {
            std::string joined;
            for (const std::string& source : block.alleles[id].source_ids) joined += (joined.empty() ? "" : ",") + source;
            return joined;
        };
        {
            std::ofstream out = open_output(prefix + ".genotypes.tsv");
            out << "block_id\tkind\tallele1\tallele2\tallele1_bp\tallele2_bp\tallele1_sources\t"
                   "allele2_sources\tzygosity\tconditional_genotype_posterior\tconditional_GQ\n";
            for (std::size_t b = 0; b < block_count; ++b) {
                const PreparedGenotypeBlock& block = model.blocks()[b];
                if (block.kind != GenotypeBlockKind::Variant) continue;
                const auto [a1, a2] = alleles_at(b);
                out << block.block_id << "\tSV\t" << a1 << '\t' << a2 << '\t' << block.alleles[a1].sequence.size() << '\t'
                    << block.alleles[a2].sequence.size() << '\t' << sources_of(block, a1) << '\t' << sources_of(block, a2)
                    << '\t' << (a1 == a2 ? "HOM" : "HET") << '\t' << std::setprecision(12)
                    << summary.block_genotype_posterior[b] << '\t' << summary.block_genotype_quality[b]
                    << std::setprecision(6) << '\n';
            }
        }
        {
            std::ofstream out = open_output(prefix + ".block_genotypes.tsv");
            out << "block_id\tkind\treportable\tallele1\tallele2\tallele1_bp\tallele2_bp\tzygosity\t"
                   "projection_ambiguous\tgenotype_support\tphase_support\tconditional_genotype_posterior\t"
                   "conditional_GQ\tnumerically_separated\tcandidate_scope\n";
            for (std::size_t b = 0; b < block_count; ++b) {
                const PreparedGenotypeBlock& block = model.blocks()[b];
                const bool variant = block.kind == GenotypeBlockKind::Variant;
                const auto [a1, a2] = alleles_at(b);
                out << block.block_id << '\t' << (variant ? "SV" : "context") << '\t' << (variant ? "true" : "false")
                    << '\t' << a1 << '\t' << a2 << '\t' << block.alleles[a1].sequence.size() << '\t'
                    << block.alleles[a2].sequence.size() << '\t' << (a1 == a2 ? "HOM" : "HET") << '\t'
                    << (summary.projection_ambiguous ? "true" : "false") << '\t'
                    << (summary.block_genotype_read_supported[b] ? "READS" : "PRIOR_ONLY") << '\t'
                    << summary.phase_support << '\t' << std::setprecision(12) << summary.block_genotype_posterior[b]
                    << '\t' << summary.block_genotype_quality[b] << std::setprecision(6) << '\t'
                    << (summary.block_numerically_separated[b] ? "true" : "false") << '\t'
                    << set_report.final_candidate_scope << '\n';
            }
        }
        if (unique) {
            // A unique candidate's two haplotypes. Under PRIOR_ONLY the reads gave no preference among the
            // phase alternatives and the prior chose this one: a prior-selected haplotype representation,
            // not a read-phased result, and the header says so.
            report.haplotype_representation =
                summary.phase_support == "PRIOR_ONLY" ? "prior_selected" : "read_supported_within_candidates";
            std::ofstream out = open_output(prefix + ".mosaic_haplotypes.fa");
            const HaplotypeSequence first = spell_mosaic(model, MosaicPath{called.first});
            const HaplotypeSequence second = spell_mosaic(model, MosaicPath{called.second});
            for (const auto& item : {std::make_pair("hap1", &first), std::make_pair("hap2", &second)}) {
                out << '>' << prefix << '_' << item.first << " fingerprint=" << fingerprint
                    << " bp=" << item.second->bases.size() << " phase_support=" << summary.phase_support
                    << " representation=" << report.haplotype_representation << '\n';
                for (std::size_t at = 0; at < item.second->bases.size(); at += 60) {
                    out << item.second->bases.substr(at, 60) << '\n';
                }
            }
            return stop_with(GenotypeRunStatus::Complete,
                             "unique best candidate within the pooled candidate set; phase_support " +
                                 summary.phase_support,
                             &decided.ranking, &summary);
        }
        report.haplotype_representation = "none";
        return stop_with(GenotypeRunStatus::ConditionalGenotypeEquivalent,
                         std::to_string(decided.top_groups.size()) +
                             " candidates tied on total score carry one unordered genotype at every block; the "
                             "genotype is emitted unphased, phase is UNRESOLVED, the equivalence set is reported "
                             "and no haplotype FASTA is written",
                         &decided.ranking, &summary);
    };

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

        // ---- tie-triggered full-beam exact audit -----------------------------------------------------
        // A K-truncated prefix is not grounds for consensus. The candidates below K were constructed and
        // ordered by proposal mass, but never exactly scored, so their order is no proof that none of
        // them wins or disagrees blockwise. When the primary ranking ties, the search was exhaustive, and
        // the beam retained candidates outside K, the complete retained beam is scored exactly and the
        // output is derived from THAT ranking alone -- never from the prefix presented beside it.
        //
        // No second scoring pass and no new tuning parameter: the shared pass already scored the full
        // beam whenever the search was exhaustive, so this re-ranks candidates whose scores, sequence
        // interning and haploid-path priors are all already computed.
        SetRanking audited;
        const bool primary_tied = primary.ranking.status != RescoreStatus::Complete;
        const bool beam_has_more = mosaics.candidates.size() > beam.candidates.size();
        const bool audit_attempted = primary_tied && mosaics.search_exhaustive && beam_has_more;
        if (audit_attempted) audited = rank_set(mosaics);

        // A cap reached while scoring the audited set fails closed with no tables, and the gate for that
        // is the rescore_failure check ABOVE, not one here: score_diplotypes runs once, in the shared
        // pass, which already covers the full beam whenever the search was exhaustive. rank_set only
        // re-ranks scores that pass produced, so it cannot reach a cap and cannot set rescore_failure.
        // A check at this point would never fire -- an inert guard that reads as protection.

        // Everything below reads through these names, so binding them to the audited ranking is what
        // makes the candidate and score tables describe the set that determined the call.
        const SetRanking& chosen = audit_attempted ? audited : primary;
        const std::vector<MosaicCandidate>& candidates = chosen.candidates;
        const std::size_t panel_count = panel_candidates.size();
        const std::vector<SequenceCandidate>& sequence_candidates = chosen.sequence_candidates;
        const std::vector<std::size_t>& sequence_candidate_of = chosen.sequence_candidate_of;
        const std::vector<double>& path_prior = chosen.path_prior;
        const CertifiedRanking& ranking = chosen.ranking;
        RankingSummary& summary = audit_attempted ? audited.summary : primary.summary;
        CandidateSetReport set_report = candidate_report;
        // Eight separate facts, reported separately. candidate_set_conditional keeps its meaning -- it
        // describes the PREFIX -- and final_candidate_scope describes the set that decided the outcome.
        set_report.primary_prefix_k = beam.candidates.size();
        set_report.full_beam_candidates = mosaics.candidates.size();
        set_report.primary_prefix_truncated = beam_has_more;
        set_report.full_beam_audit_attempted = audit_attempted;
        set_report.candidates_exactly_audited = audit_attempted ? audited.candidates.size() : 0;
        set_report.audit_completed = audit_attempted;
        set_report.final_candidate_scope =
            audit_attempted ? candidate_scope_name(panel_only, false)
                            : candidate_scope_name(panel_only, beam.candidate_set_conditional);
        {
            std::set<std::size_t> haplotypes;
            for (const SequenceCandidate& candidate : sequence_candidates) {
                haplotypes.insert(candidate.first);
                haplotypes.insert(candidate.second);
            }
            set_report.unique_haplotypes = haplotypes.size();
        }
        if (options.repair_rounds > 0) return adaptive_completion(beam, chosen, set_report, timings);

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
            // Full round-trip precision for a double. At C4 these scores are around 1.7e6, so the
            // stream default of six significant digits rounded them to the nearest ten -- enough to
            // read, not enough to compare or reproduce. The standalone-versus-production agreement
            // gate gets its resolution from here, and a result artifact should not lose digits the
            // search actually computed.
            out << std::setprecision(17) << manifest_header;
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
            timings.candidate_grouping = run_timings.candidate_grouping;
            timings.prior = run_timings.prior;
            timings.ranking = run_timings.ranking;
            // The counters are run-level and the audit's ranking contributes to them, so they are read
            // here rather than from the copy taken before the audit ran -- which reported the primary
            // set's 130 paths and made "the audit evaluates only 130" look true when it is not.
            set_report.distinct_haplotype_paths = candidate_report.distinct_haplotype_paths;
            set_report.haploid_prior_evaluations = candidate_report.haploid_prior_evaluations;
            set_report.diploid_prior_lookups = candidate_report.diploid_prior_lookups;
            // BLOCK CONSENSUS. Scoring finished; what could not be done is separating a single
            // whole-locus diplotype. certify_scored_ranking sets RescoreIncomplete in exactly one place
            // -- more than one candidate inseparable from the best -- and every other incompleteness
            // (model, resource, candidate construction, spelling) has already returned above. Reaching
            // here therefore means a numerical tie and nothing else.
            //
            // Different phased histories can spell the same unordered genotype at every block, and short
            // reads may carry no information that separates them. Where every tied candidate agrees on a
            // block, that block's genotype is supported by the whole tied set and is emitted. Where they
            // disagree, only that block is no-called. Phase stays UNRESOLVED and no phased haplotype is
            // written: the reads did not choose one, and inventing it is the error being avoided.
            //
            // Refused when the candidate set itself is incomplete. A pruned search omitted candidates it
            // never proposed, so unanimity among those it happened to retain says nothing about those it
            // did not -- consensus would be agreement of an unrepresentative sample.
            // The audited full beam of an exhaustive search omits nothing, so consensus is permitted on
            // it. A K-truncated prefix is still refused: that is the whole point of auditing first.
            const bool candidates_complete =
                panel_only || audit_attempted || !beam.candidate_set_conditional;
            const std::size_t block_count = model.blocks().size();
            std::vector<std::vector<std::pair<std::size_t, std::size_t>>> tied_genotypes;
            if (candidates_complete) {
                for (const std::size_t id : ranking.winner_set) {
                    const MosaicDiplotype& tied =
                        candidates[sequence_candidates[id].representative].diplotype;
                    std::vector<std::pair<std::size_t, std::size_t>> row(block_count);
                    for (std::size_t b = 0; b < block_count; ++b) {
                        row[b] = std::minmax(tied.first.alleles[b], tied.second.alleles[b]);
                    }
                    tied_genotypes.push_back(std::move(row));
                }
            }
            // The rule itself is the gated pure function, so what ships is what the fixtures test.
            const std::vector<BlockConsensusCall> consensus_calls = block_consensus(tied_genotypes);
            std::vector<std::pair<std::size_t, std::size_t>> consensus(block_count, {0, 0});
            std::vector<bool> agreed(block_count, false);
            std::size_t consensus_called = 0;
            for (std::size_t b = 0; b < block_count && b < consensus_calls.size(); ++b) {
                agreed[b] = consensus_calls[b].agreed;
                consensus[b] = {consensus_calls[b].allele1, consensus_calls[b].allele2};
                consensus_called += agreed[b] ? 1 : 0;
            }

            if (candidates_complete && consensus_called > 0) {
                // The reads chose no whole-locus history, so no phase claim survives, whatever the
                // within-candidate analysis concluded before the tie was known.
                summary.phase_support = "UNRESOLVED";
                summary.consensus_blocks_called = consensus_called;
                summary.consensus_blocks_no_called = block_count - consensus_called;
                const auto sources_of = [&](const PreparedGenotypeBlock& block, std::size_t id) {
                    std::string joined;
                    for (const std::string& source : block.alleles[id].source_ids) {
                        if (!joined.empty()) joined += ',';
                        joined += source;
                    }
                    return joined;
                };
                {
                    std::ofstream out = open_output(options.out_prefix + ".genotypes.tsv");
                    out << "block_id\tkind\tallele1\tallele2\tallele1_bp\tallele2_bp\tallele1_sources\t"
                           "allele2_sources\tzygosity\tconditional_genotype_posterior\tconditional_GQ\n";
                    for (std::size_t b = 0; b < block_count; ++b) {
                        const PreparedGenotypeBlock& block = model.blocks()[b];
                        if (block.kind != GenotypeBlockKind::Variant) continue;
                        if (!agreed[b]) {
                            // Emitted as an explicit no-call rather than omitted: a missing row is
                            // indistinguishable from a block that was never modelled.
                            out << block.block_id << "\tSV\t.\t.\t.\t.\t.\t.\tNOCALL\tNA\tNA\n";
                            continue;
                        }
                        const std::size_t a1 = consensus[b].first;
                        const std::size_t a2 = consensus[b].second;
                        out << block.block_id << "\tSV\t" << a1 << '\t' << a2 << '\t'
                            << block.alleles[a1].sequence.size() << '\t'
                            << block.alleles[a2].sequence.size() << '\t' << sources_of(block, a1) << '\t'
                            << sources_of(block, a2) << '\t' << (a1 == a2 ? "HOM" : "HET") << '\t'
                            << std::setprecision(12) << summary.block_genotype_posterior[b] << '\t'
                            << summary.block_genotype_quality[b] << std::setprecision(6) << '\n';
                    }
                }
                {
                    std::ofstream out = open_output(options.out_prefix + ".block_genotypes.tsv");
                    out << "block_id\tkind\treportable\tallele1\tallele2\tallele1_bp\tallele2_bp\t"
                           "zygosity\tprojection_ambiguous\tgenotype_support\tphase_support\t"
                           "conditional_genotype_posterior\tconditional_GQ\tnumerically_separated\t"
                           "candidate_scope\n";
                    for (std::size_t b = 0; b < block_count; ++b) {
                        const PreparedGenotypeBlock& block = model.blocks()[b];
                        const bool variant = block.kind == GenotypeBlockKind::Variant;
                        out << block.block_id << '\t' << (variant ? "SV" : "context") << '\t'
                            << (variant ? "true" : "false") << '\t';
                        if (!agreed[b]) {
                            out << ".\t.\t.\t.\tNOCALL\tfalse\tNOCALL\tUNRESOLVED\tNA\tNA\tfalse\t"
                                << set_report.final_candidate_scope
                                << '\n';
                            continue;
                        }
                        const std::size_t a1 = consensus[b].first;
                        const std::size_t a2 = consensus[b].second;
                        out << a1 << '\t' << a2 << '\t' << block.alleles[a1].sequence.size() << '\t'
                            << block.alleles[a2].sequence.size() << '\t' << (a1 == a2 ? "HOM" : "HET")
                            << "\tfalse\t"
                            << (summary.block_genotype_read_supported[b] ? "READS" : "PRIOR_ONLY")
                            << "\tUNRESOLVED\t" << std::setprecision(12)
                            << summary.block_genotype_posterior[b] << '\t'
                            << summary.block_genotype_quality[b] << std::setprecision(6) << '\t'
                            << (summary.block_numerically_separated[b] ? "true" : "false") << '\t'
                            << set_report.final_candidate_scope << '\n';
                    }
                }
                timings.total = seconds_since(started);
                write_status(status_path, GenotypeRunStatus::BlockConsensus, fingerprint, options,
                             library, timings, &built.report, &panel, (panel_only ? nullptr : &beam),
                             &ranking,
                             "no whole-locus winner could be separated from its rivals beyond "
                             "floating-point error; blocks on which every tied candidate agrees are "
                             "called by consensus and the rest are no-called, with phase UNRESOLVED and "
                             "no phased haplotype emitted",
                             &summary, &depth, &set_report);
                if (!options.quiet) {
                    std::cerr << "[genotype] BLOCK_CONSENSUS  fingerprint=" << fingerprint
                              << "  tied=" << ranking.winner_set.size() << "  called=" << consensus_called
                              << "/" << block_count << "\n";
                }
                return GenotypeRunStatus::BlockConsensus;
            }

            timings.total = seconds_since(started);
            write_status(status_path, GenotypeRunStatus::RescoreIncomplete, fingerprint, options,
                         library, timings, &built.report, &panel, (panel_only ? nullptr : &beam), &ranking,
                         std::string("the best candidate cannot be separated from every rival by more "
                                     "than floating-point error; the equivalence set is reported instead "
                                     "of a call") +
                             (candidates_complete
                                  ? "; no block has the same genotype across every tied candidate, so "
                                    "there is no block consensus either"
                                  : "; the candidate set is conditional, so agreement among the "
                                    "candidates that were retained is not evidence about those that "
                                    "were not, and no block consensus is offered"),
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
                    << '\t' << set_report.final_candidate_scope << '\n';
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

        timings.candidate_grouping = run_timings.candidate_grouping;
        timings.prior = run_timings.prior;
        timings.ranking = run_timings.ranking;
        // Read after the audit, not from the copy taken before it: these count the whole run, and the
        // audit's ranking adds its own distinct mosaic paths to the shared cache.
        set_report.distinct_haplotype_paths = candidate_report.distinct_haplotype_paths;
        set_report.haploid_prior_evaluations = candidate_report.haploid_prior_evaluations;
        set_report.diploid_prior_lookups = candidate_report.diploid_prior_lookups;
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
    // The profile is applied first, wherever it appears; a flag given beside it overrides it and is
    // recorded as an override, so the status always says what the run actually ran under.
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (args[i] == "--profile") apply_genotype_profile(need(i), options);
    }
    const bool profiled = !options.profile_path.empty();
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& flag = args[i];
        if (flag == "--profile") { ++i; continue; }
        if (flag == "--gfa" || flag == "-i") options.gfa_path = need(i);
        else if (flag == "--bubbles") options.bubbles_path = need(i);
        else if (flag == "--reads1") options.reads1_path = need(i);
        else if (flag == "--reads2") options.reads2_path = need(i);
        else if (flag == "--reference-path" || flag == "-r") options.reference_path = need(i);
        else if (flag == "--out-prefix" || flag == "-o") options.out_prefix = need(i);
        else if (flag == "--library") options.library_path = need(i);
        else if (flag == "--threads") options.threads = std::stoul(need(i));   // workers for diploid combination
        else if (flag == "--panel-curation") options.curation_manifest_path = need(i);
        else if (flag == "--max-mosaic-candidates") options.max_mosaic_candidates = std::stoul(need(i));
        else if (flag == "--max-candidates") {
            throw std::runtime_error("--max-candidates is now --max-mosaic-candidates: every panel diplotype is "
                                     "always a candidate, and K counts only the off-panel mosaics added to it");
        }
        else if (flag == "-q" || flag == "--quiet") options.quiet = true;
        else if (flag == "-h" || flag == "--help") {
            std::cout <<
                "Usage:\n"
                "  panvar genotype --profile <profile.tsv> --gfa <graph.gfa> --reference-path <name> \\\n"
                "                  --reads1 <r1.fq> --reads2 <r2.fq> --library <params.txt> \\\n"
                "                  --out-prefix <prefix> [--panel-curation <manifest.tsv>] [--bubbles <bubbles.csv>]\n\n"
                "Calls a pair of whole-locus haplotypes and projects them to one allele pair per SV\n"
                "block. Paired reads only. The command never asks what kind of sample it is given.\n\n"
                "With a profile (the C4 production profile is profiles/genotype_c4.tsv) the call is\n"
                "panel-anchored adaptive completion:\n"
                "  1. every unordered pair of panel haplotypes and up to K grouped-prefix mosaics are\n"
                "     exactly scored;\n"
                "  2. the best exactly scored panel diplotype anchors an exact one-block coordinate\n"
                "     repair, which moves to the best improving neighbour until none improves;\n"
                "  3. the union and every candidate the repair scored are pooled and exactly reranked;\n"
                "  4. one candidate at the top calls (status COMPLETE); several tied on total score\n"
                "     with one unordered genotype at every block emit that genotype unphased\n"
                "     (CONDITIONAL_GENOTYPE_EQUIVALENT, no FASTA); tied candidates that disagree at\n"
                "     any block refuse. A round cap, resource refusal or interruption emits no genotype.\n"
                "The candidate set is always CONDITIONAL: the search is local, never exhaustive, and\n"
                "conditional_GQ is uncalibrated. phase_support is reported separately: PRIOR_ONLY\n"
                "means the reads give no preference between phase alternatives, and the FASTA of\n"
                "such a call is a prior-selected haplotype representation, not a read-phased result.\n\n"
                "Without a profile the unrepaired call from the proposal's candidate set is made\n"
                "(a development view). --max-mosaic-candidates sets K (0: the panel alone) and\n"
                "--threads the workers; given beside a profile, either is recorded as an override.\n\n"
                "Outputs <prefix>.genotypes.tsv, .block_genotypes.tsv, .mosaic_candidates.tsv,\n"
                ".mosaic_scores.tsv, .mosaic_status.tsv and .mosaic_haplotypes.fa; with a profile\n"
                "also .equivalence_set.tsv, .repair_rounds.tsv and .repair_memory.tsv. While running,\n"
                "the status says RUNNING and committed repair rounds are in .repair_rounds.partial.tsv;\n"
                "a run stopped by the profile's footprint or time limit keeps only those and says\n"
                "RESOURCE_REFUSAL or TIME_LIMIT. An incomplete run never presents a partial ranking\n"
                "as a result. Exit status 0 when a genotype was emitted, 1 otherwise.\n\n"
                "--panel-curation takes a manifest of panel templates to exclude: one template name\n"
                "per line, '#' opening a comment line. It is recorded with the reason, the manifest\n"
                "path and a digest, and a manifest naming an absent template is refused.\n\n"
                "--library declares insert_mean_bp, insert_sd_bp, support_sigmas,\n"
                "substitution_error, outlier_weight and fragment_start_intensity, one per line,\n"
                "and optionally marker_background (estimated from invariant anchors when absent).\n"
                "Version one does not estimate these from the data.\n";
            return 0;
        } else {
            throw std::runtime_error("unknown genotype option: " + flag);
        }
    }
    if (profiled) {
        // Every setting a flag gave beside the profile, with the value it took.
        for (const std::string& flag : args) {
            if (flag == "--threads") options.overridden_settings.push_back("threads=" + std::to_string(options.threads));
            if (flag == "--max-mosaic-candidates") {
                options.overridden_settings.push_back("max_mosaic_candidates=" +
                                                      std::to_string(options.max_mosaic_candidates));
            }
        }
    }
    const GenotypeRunStatus status = run_genotype(options);
    return status == GenotypeRunStatus::Complete || status == GenotypeRunStatus::ConditionalGenotypeEquivalent ? 0 : 1;
}

} // namespace panvar
