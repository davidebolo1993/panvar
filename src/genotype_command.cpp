#include "panvar/genotype_command.hpp"
#include "panvar/genotype_assessment.hpp"
#include "panvar/genotype_sequence_report.hpp"
#include "panvar/prepared_locus.hpp"
#include "panvar/prepared_locus_file.hpp"
#include "prepared_locus_parts.hpp"
#include "panvar/genotype_dosage_integrate.hpp"
#include "panvar/genotype_dosage_repair.hpp"
#include "panvar/gz_reader.hpp"
#include "panvar/library_estimate.hpp"

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
#include "panvar/genotype_block_support.hpp"
#include "panvar/genotype_joint_search.hpp"
#include "panvar/copying_transitions.hpp"
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
#include <numeric>
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
                                     ".repair_rounds.tsv", ".repair_memory.tsv", ".tie_branches.tsv",
                                     ".tie_terminal_ties.tsv", ".tie_branch_memory.tsv", ".library_sensitivity.tsv",
                                     ".library_sensitivity_blocks.tsv", ".settings.tsv"};
const char* const kPartialOutputs[] = {".repair_rounds.partial.tsv", ".repair_memory.partial.tsv",
                                       ".tie_branch_memory.partial.tsv"};
// The representation-conflict filter columns, appended to both genotype tables on every path; a path that
// does not evaluate the guard says so rather than claiming PASS.
const char* const kFilterHeader = "\tfilter\tfilter_reason\trepresentation_marker_margin\trepresentation_competitor_pair"
                                  "\trepresentation_fragment_diff\trepresentation_exposure_diff"
                                  "\trepresentation_likelihood_margin\trepresentation_conflict_statistic";
const char* const kFilterNotEvaluated = "\tNOT_EVALUATED\tguard_not_run_on_this_path\tNA\tNA\tNA\tNA\tNA\tNA";
// The combined block quality (block_quality), appended after the representation filter columns on every path.
const char* const kQualityHeader = "\tquality_filter\tquality_reason";
constexpr const char* kQualityRule = "block_quality_rule_v1";
constexpr const char* kRepresentationRule = "representation_conflict_rule_v1";

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
    // S2-A (early_incumbent experimental_v1): the last committed checkpoint's locus record (what a budget stop keeps,
    // rewritten with the stop), or a checkpoint publication that failed (an explicit failure, never a fallback)
    bool early_incumbent = false;
    bool checkpoint_committed = false;
    bool checkpoint_failed = false;
    std::string checkpoint_stage = "none";
    std::string checkpoint_locus;
    // M15-D.1: a sequence_first_v2 run reports its termination in .sequence_locus.tsv even without the early incumbent
    bool sequence_v2 = false;
};
// The one run in flight in this process, or null. Runs are sequential within a process.
RunWatch* g_watch = nullptr;
// Provenance of the inputs, by content, for the run in flight.
std::vector<std::pair<std::string, std::string>> g_input_hashes;
// I2-E: what the marker-proposal contribution reuse did in this run (reported when its key is set)
ProposalReuseReport g_proposal_reuse;
// S1-B (sequence_first_v1 only): this run's sequence report and its legacy classification (reset per run)
SequenceReport g_sequence_report;
// S1-D: the sequence-first schemas (v1 reports every prediction's FASTA; v2 adds panel references, the panel-membership
// FASTA rule and the diagnostic-only dosage contract)
bool sequence_first(const GenotypeCommandOptions& options) {
    return options.output_schema == "sequence_first_v1" || options.output_schema == "sequence_first_v2";
}
std::string dosage_mode_source(const GenotypeCommandOptions& options) { return options.dosage_mode_source; }
// S2-A: the early_incumbent contract (also checked before an API run touches its prefix).
std::optional<std::string> early_incumbent_refusal(const GenotypeCommandOptions& o) {
    if (o.early_incumbent != "experimental_v1") {
        if (o.union_progressive_schedule != "off") return std::string("union_progressive_schedule needs early_incumbent experimental_v1");
        if (o.test_interrupt_point != "none") return std::string("test_interrupt_point needs early_incumbent experimental_v1");
        if (o.test_perturb_early_incumbent_score) return std::string("test_perturb_early_incumbent_score needs early_incumbent experimental_v1");
        return std::nullopt;
    }
    const std::string what = "early_incumbent experimental_v1 ";
    if (o.output_schema != "sequence_first_v2") {
        return what + "needs output_schema sequence_first_v2 (its checkpoint follows the v2 panel-reference and "
                      "conditional-FASTA policy)";
    }
    if (o.repair_rounds == 0) return what + "needs the adaptive repair (repair rounds > 0)";
    if (o.diagnostic_union_only) return what + "is not supported with diagnostic_union_only";
    if (!o.development_candidate_prefixes.empty()) return what + "is not supported with development candidate prefixes";
    if (o.union_progressive_schedule == "experimental_v1") {   // S2-D
        const std::string p = "union_progressive_schedule experimental_v1 ";
        if (o.union_scoring_backend != "indexed_fine_experimental") return p + "needs union_scoring_backend indexed_fine_experimental";
        if (o.indexed_union_store_layout == "legacy_i2a") return p + "needs the exact indexed store layout";
        if (o.adaptive_emission_handoff) return p + "is not supported with adaptive_emission_handoff";
        if (o.diagnostic_union_evidence_dump || o.diagnostic_union_work_trace || o.diagnostic_union_seeding_profile > 0 ||
            !o.diagnostic_union_verification_census.empty()) {
            return p + "is not supported with union diagnostics that need the full mass matrix";
        }
    }
    return std::nullopt;
}
// S2-G: the blocked union route's contract (also checked before an API run touches its prefix).
std::optional<std::string> union_fragment_blocks_refusal(const GenotypeCommandOptions& o) {
    if (o.union_fragment_blocks == "off") {
        if (o.test_union_fragment_block_pairs != 0) {
            return std::string("test_union_fragment_block_pairs needs union_fragment_blocks fine_core_4096_experimental");
        }
        return std::nullopt;
    }
    const std::string p = "union_fragment_blocks fine_core_4096_experimental ";
    if (o.union_scoring_backend != "indexed_fine_experimental") return p + "needs union_scoring_backend indexed_fine_experimental";
    if (o.indexed_union_store_layout == "legacy_i2a") return p + "needs the exact indexed store layout";
    if (o.indexed_union_evidence_layout != "fine_views_experimental") {
        return p + "needs indexed_union_evidence_layout fine_views_experimental";
    }
    if (o.union_origin_cache_bytes == 0) return p + "needs union origin reuse (a positive union_origin_cache_bytes)";
    if (o.max_union_cached_origins_consumed == 0 || o.max_union_cached_states_consumed == 0) {
        return p + "needs both union origin consumption caps, positive";
    }
    if (o.union_progressive_schedule != "off") return p + "is not supported with union_progressive_schedule";
    if (o.adaptive_emission_handoff) return p + "is not supported with adaptive_emission_handoff";
    if (o.diagnostic_union_only) return p + "is not supported with diagnostic_union_only (its digest needs the full mass matrix)";
    if (o.diagnostic_union_evidence_dump || o.diagnostic_union_work_trace || o.diagnostic_union_seeding_profile > 0 ||
        !o.diagnostic_union_verification_census.empty()) {
        return p + "is not supported with union diagnostics that need the full mass matrix";
    }
    return std::nullopt;
}
// M3-D: the joint-search experiment's contract (checked before an API run touches its prefix).
// M4-B: the copying-transition specification's contract (checked before an API run touches its prefix). The distance mode
// is refused, never silently replaced by constant transitions, wherever it cannot be honoured.
std::optional<std::string> copying_transition_refusal(const GenotypeCommandOptions& o) {
    const bool length_set = o.setting_sources.count("copying_length_bp") != 0;
    if (o.copying_transition == "constant_v1") {
        if (length_set) return std::string("copying_length_bp is a distance_ref_core_v1_experimental parameter; it is refused with constant_v1");
        return std::nullopt;
    }
    const std::string p = "copying_transition distance_ref_core_v1_experimental ";
    if (o.copying_transition != "distance_ref_core_v1_experimental") return "copying_transition '" + o.copying_transition + "' is not supported";
    if (!length_set || !std::isfinite(o.copying_length_bp) || !(o.copying_length_bp > 0.0)) return p + "needs copying_length_bp > 0";
    const auto sw = o.setting_sources.find("switch_probability");
    if (sw != o.setting_sources.end() && sw->second.rfind("preset ", 0) != 0) {
        return p + "conflicts with switch_probability set by " + sw->second + " (a scalar switch setting and a distance specification are exclusive)";
    }
    if (o.repair_joint_search != "off") return p + "is not supported with repair_joint_search (joint search stays off)";
    if (!o.diagnostic_refinement_replay.empty()) return p + "is not supported with the refinement replay";
    if (o.diagnostic_tie_branching) return p + "is not supported with diagnostic tie branching";
    return std::nullopt;
}
CopyingTransitions g_transitions;   // M4-B: this run's one transition vector (set once the model exists)
LocusGeometry g_geometry;           // M4-B: the adapter's geometry of this run's model (when built here)
std::vector<double> run_switches(std::size_t blocks) {
    if (g_transitions.r.size() + 1 != blocks) {
        throw std::logic_error("the run's copying-transition vector does not match the block count (" + std::to_string(g_transitions.r.size()) +
                               " boundaries for " + std::to_string(blocks) + " blocks)");
    }
    return g_transitions.r;
}
std::optional<std::string> compound_annotation_refusal(const GenotypeCommandOptions& o) {   // M9-J
    if (o.compound_annotation == "off") {
        if (!o.compound_interval_declarations.empty()) return std::string("compound_interval_declarations needs compound_annotation contract_v2");
        return std::nullopt;
    }
    if (o.output_schema != "sequence_first_v2") return std::string("compound_annotation contract_v2 needs output_schema sequence_first_v2");
    return std::nullopt;
}
std::optional<CompoundDeclarations> g_compound_declarations;   // M9-J: this run's declarations (read before the prefix)
void compound_prepare(SequenceReport& r, const GenotypeCommandOptions& o) {
    if (o.compound_annotation != "contract_v2") return;
    r.compound.enabled = true;
    if (g_compound_declarations) {
        r.compound.declaration_path = g_compound_declarations->path;
        r.compound.declaration_sha256 = g_compound_declarations->sha256;
    }
}
// M9-J: the completed pool -- every sequence group, every recorded labelling, every group's full-precision mass.
void compound_pooled(SequenceReport& r, const PooledDecision& d, const PreparedGenotypeModel& model, const std::string& fingerprint,
                     const GenotypeCommandOptions& o) {
    if (o.compound_annotation != "contract_v2") return;
    const auto mark = std::chrono::steady_clock::now();
    compound_prepare(r, o);
    std::vector<std::string> ids;
    for (const auto& block : model.blocks()) ids.push_back(block.block_id);
    std::vector<std::pair<std::size_t, std::size_t>> representation;
    for (const RepresentationInterval& iv : r.intervals) representation.push_back({iv.first_block, iv.last_block});
    const std::vector<CompoundSpan> spans =
        resolve_compound_spans(g_compound_declarations ? &*g_compound_declarations : nullptr, fingerprint, ids, representation, &r.compound.declaration_status);
    std::vector<double> mass(d.groups.size(), 0.0);
    for (const CertifiedCandidate& c : d.ranking.candidates) mass[c.candidate_id] = c.candidate_set_probability;
    const std::set<std::size_t> top(d.top_groups.begin(), d.top_groups.end());
    std::vector<CompoundGroup> groups(d.groups.size());
    for (std::size_t g = 0; g < d.groups.size(); ++g) {
        groups[g].mass = mass[g];
        groups[g].decision_set = top.count(g) != 0;
        for (const std::size_t m : d.groups[g]) groups[g].members.push_back({&d.candidates[m].first, &d.candidates[m].second});
    }
    const CompoundAlleleBytes bytes = [&](std::size_t b, std::size_t a) -> const std::string& { return model.blocks()[b].alleles[a].sequence; };
    r.compound.rows.clear();
    for (const CompoundSpan& span : spans) {
        CompoundAnnotationRow row = compound_assess(bytes, span.first, span.last, groups, r.representative_first, r.representative_second);
        for (const std::string& s : span.sources) row.sources += (row.sources.empty() ? "" : ",") + s;
        r.compound.rows.push_back(std::move(row));
    }
    r.compound.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - mark).count();
}
// M9-J: a report without a completed pool (provisional): declared spans only, no bounds.
void compound_unassessed(SequenceReport& r, const PreparedGenotypeModel& model, const std::string& fingerprint, const GenotypeCommandOptions& o,
                         const std::string& status) {
    if (o.compound_annotation != "contract_v2") return;
    compound_prepare(r, o);
    std::vector<std::string> ids;
    for (const auto& block : model.blocks()) ids.push_back(block.block_id);
    const std::vector<CompoundSpan> spans =
        resolve_compound_spans(g_compound_declarations ? &*g_compound_declarations : nullptr, fingerprint, ids, {}, &r.compound.declaration_status);
    r.compound.rows.clear();
    for (const CompoundSpan& span : spans) {
        CompoundAnnotationRow row;
        row.first = span.first;
        row.last = span.last;
        row.status = status;
        for (const std::string& s : span.sources) row.sources += (row.sources.empty() ? "" : ",") + s;
        r.compound.rows.push_back(std::move(row));
    }
}
std::optional<std::string> junction_proposals_refusal(const GenotypeCommandOptions& o) {   // M9-G
    if (o.repair_junction_proposals == "off") return std::nullopt;
    const std::string p = "repair_junction_proposals experimental_v1 ";
    if (o.repair_rounds == 0) return p + "needs the adaptive repair (repair rounds > 0)";
    if (o.diagnostic_tie_branching) return p + "is not supported with diagnostic tie branching (it enumerates the index-ordered lists)";
    if (!o.diagnostic_refinement_replay.empty()) return p + "is not supported with the refinement replay";
    return std::nullopt;
}
std::optional<std::string> joint_search_refusal(const GenotypeCommandOptions& o) {
    if (o.repair_joint_search == "off") return std::nullopt;
    const std::string p = "repair_joint_search experimental_v1 ";
    if (o.repair_rounds == 0) return p + "needs the adaptive repair (repair rounds > 0)";
    if (o.library_path == "auto" && o.library_sensitivity_audit) {
        return p + "is not supported with the library sensitivity audit (set library_sensitivity_audit 0 or --library FILE)";
    }
    if (o.diagnostic_tie_branching) return p + "is not supported with diagnostic tie branching";
    if (!o.diagnostic_refinement_replay.empty()) return p + "is not supported with the refinement replay";
    if (o.diagnostic_union_only) return p + "is not supported with diagnostic_union_only";
    return std::nullopt;
}
// S2-H: the refinement diagnostics' contract (also checked before an API run touches its prefix).
std::optional<std::string> refinement_diagnostic_refusal(const GenotypeCommandOptions& o) {
    if (o.diagnostic_refinement_replay.empty()) return std::nullopt;
    const std::string p = "diagnostic_refinement_replay ";
    if (o.output_schema != "legacy") return p + "needs output_schema legacy (a diagnostic replay publishes no sequence report)";
    if (o.early_incumbent != "off") return p + "needs early_incumbent off (no early pass runs before a replayed round)";
    if (o.repair_rounds == 0) return p + "needs the adaptive repair (repair rounds > 0)";
    if (o.diagnostic_union_only || o.union_progressive_schedule != "off" || o.adaptive_emission_handoff) {
        return p + "is not supported with union diagnostics, the progressive schedule or the emission handoff";
    }
    if (!o.development_candidate_prefixes.empty()) return p + "is not supported with development candidate prefixes";
    return std::nullopt;
}
// M9-G: the run's junction-proposal account, summed over its searches (eligibility is a property of the model: the maximum).
void accumulate_junction(JunctionProposalReport& total, const JunctionProposalReport& one) {
    total.eligible_blocks = std::max(total.eligible_blocks, one.eligible_blocks);
    if (!one.block_eligible.empty()) {   // M11-C: a block is reported eligible when any search evaluated it as eligible
        if (total.block_eligible.empty()) total.block_eligible.assign(one.block_eligible.size(), 0);
        for (std::size_t b = 0; b < one.block_eligible.size() && b < total.block_eligible.size(); ++b) {
            total.block_eligible[b] = static_cast<unsigned char>(total.block_eligible[b] | one.block_eligible[b]);
        }
    }
    total.invocations += one.invocations;
    total.cache_hits += one.cache_hits;
    total.cutoff_tie_crossings += one.cutoff_tie_crossings;
    total.seconds += one.seconds;
}
// S2-H: the refinement's scoring options and search options, from the run's options and union scoring options -- ONE
// construction for the production search and the diagnostic replay.
struct RepairSetup {
    FragmentScoringOptions scoring;
    RepairSearchOptions search;
};
RepairSetup make_repair_setup(const GenotypeCommandOptions& options, const FragmentScoringOptions& scoring_options) {
    RepairSetup r;
    r.scoring = scoring_options;
    r.scoring.max_physical_origins = options.repair_max_physical_origins;
    r.scoring.max_state_classes = options.repair_max_state_classes;
    r.scoring.max_piece_hits = options.repair_max_piece_hits;
    r.scoring.max_latent_alignment_pairs = options.repair_max_latent_pairs;
    RepairSearchOptions& s = r.search;
    s.local_pairs = options.repair_local_pairs;
    s.max_rounds = options.repair_rounds;
    s.mode = RepairScoringMode::Delta;
    s.origin_cache = true;
    s.planned_cache = true;
    s.origin_cache_bytes = options.repair_origin_cache_bytes;
    if (options.adaptive_window_discovery == "fine_experimental") {
        s.fine_window_bytes = options.adaptive_fine_store_bytes;
        s.fine_window_served_cap = options.adaptive_fine_max_served_placements;
    }
    s.fused_seed = options.adaptive_fused_seed;
    s.switch_probability = options.switch_probability;
    if (g_transitions.distance_mode()) s.switch_vector = g_transitions.r;   // M4-B (constant mode: the scalar, as before)
    s.declared_neighbour_cap = options.repair_max_neighbours;
    s.memory_log = true;
    if (options.diagnostic_refinement_memory) s.memory_live_path = options.out_prefix + ".refinement_memory_live.tsv";
    s.round1_reuse_only_admission = options.repair_round1_origin_admission == "reuse_only_experimental";   // S2-J
    return r;
}
std::string hex_double(double value) {
    char buffer[64];
    std::snprintf(buffer, sizeof(buffer), "%a", value);
    return buffer;
}
// S2-H: the refinement setup as text (what a replay must match exactly; max_rounds as configured)
std::string repair_setup_text(const RepairSetup& r, std::size_t fragments) {
    const FragmentScoringOptions& o = r.scoring;
    const RepairSearchOptions& s = r.search;
    std::ostringstream t;
    t << "fragments=" << fragments << ";threads=" << o.threads << ";band=" << o.mismatch_band << ";piece=" << o.minimum_piece_bp
      << ";max_piece_hits=" << o.max_piece_hits << ";max_starts=" << o.max_proposed_starts_per_haplotype
      << ";max_origins=" << o.max_physical_origins << ";max_classes=" << o.max_state_classes << ";max_cells=" << o.max_mass_cells
      << ";max_latent=" << o.max_latent_alignment_pairs << ";emission_cache=" << o.emission_cache_bytes
      << ";prepared_terms=" << o.prepared_library_terms << ";gap_backend=" << static_cast<int>(o.gap_recursion_backend)
      << ";audits=" << o.audit_libraries.size() << ";union_placement_cache=" << o.union_placement_cache_bytes
      << ";union_origin_cache=" << o.union_origin_cache_bytes << ";union_prefilter=" << o.union_prefilter_cache_bytes
      << ";start_sort=" << o.start_sort_scratch_bytes << ";local_pairs=" << s.local_pairs << ";max_rounds=" << s.max_rounds
      << ";origin_cache_bytes=" << s.origin_cache_bytes << ";fine_window=" << s.fine_window_bytes << ";fine_served="
      << s.fine_window_served_cap << ";fused_seed=" << s.fused_seed << ";switch=" << hex_double(s.switch_probability)
      << ";neighbour_cap=" << s.declared_neighbour_cap;
    if (s.round1_reuse_only_admission) t << ";round1_admission=reuse_only_experimental";   // S2-J (absent by default)
    if (!s.switch_vector.empty()) t << ";transitions=" << g_transitions.vector_digest;   // M4-B (absent in constant mode)
    return t.str();
}
std::string emissions_digest(const PreparedProposalEmissions& emissions) {
    std::uint64_t h = 1469598103934665603ULL;
    const auto mix = [&h](const void* data, std::size_t n) {
        const unsigned char* b = static_cast<const unsigned char*>(data);
        for (std::size_t k = 0; k < n; ++k) h = (h ^ b[k]) * 1099511628211ULL;
    };
    for (const BlockProposalEmissions& e : emissions.blocks) {
        mix(&e.allele_count, sizeof e.allele_count);
        for (std::size_t x = 0; x < e.allele_count; ++x) {
            for (std::size_t y = x; y < e.allele_count; ++y) {
                const double v = e.at(x, y);
                mix(&v, sizeof v);
            }
        }
    }
    std::ostringstream t;
    t << std::hex << h;
    return t.str();
}
// S2-G: small physical-footprint checkpoints of a blocked-route run (stage, bytes, peak, seconds since the run began);
// recorded only when the route is on, reset at every run start, written to the status (final and stopped).
struct FootprintCheckpoint {
    std::string stage;
    std::uint64_t current = 0, peak = 0;
    double seconds = 0.0;
    std::optional<std::uint64_t> allocator;   // S2-K: the allocator's in-use bytes (read-only; unavailable elsewhere)
};
std::mutex g_checkpoint_mutex;
bool g_checkpoints_on = false;
Clock::time_point g_checkpoints_started;
std::vector<FootprintCheckpoint> g_footprint_checkpoints;
void note_footprint(const std::string& stage) {
    std::lock_guard<std::mutex> lock(g_checkpoint_mutex);
    if (!g_checkpoints_on) return;
    const FootprintSample fp = physical_footprint();
    g_footprint_checkpoints.push_back({stage, fp.current, fp.peak, seconds_since(g_checkpoints_started), allocator_bytes_in_use()});
}
std::string footprint_checkpoint_rows() {
    std::lock_guard<std::mutex> lock(g_checkpoint_mutex);
    if (!g_checkpoints_on) return "";
    std::ostringstream o;
    o << std::fixed << std::setprecision(3) << "footprint_checkpoints	" << g_footprint_checkpoints.size() << '\n';
    for (std::size_t i = 0; i < g_footprint_checkpoints.size(); ++i) {
        const FootprintCheckpoint& c = g_footprint_checkpoints[i];
        const std::string k = "footprint_checkpoint" + std::to_string(i + 1) + "_";
        o << k << "stage\t" << c.stage << '\n' << k << "phys_footprint_bytes\t" << c.current << '\n' << k
          << "phys_footprint_peak_bytes\t" << c.peak << '\n' << k << "seconds\t" << c.seconds << '\n' << k
          << "allocator_in_use_bytes\t" << (c.allocator ? std::to_string(*c.allocator) : std::string("unavailable")) << '\n';
    }
    return o.str();
}
// S2-A: the run's incumbent -- the best fully evaluated candidate so far -- and the early evaluation's account (reset
// at every run start). Scores and priors are the unchanged scorer's and the ranking's own.
struct EarlyIncumbent {
    bool enabled = false;
    bool have = false;
    std::vector<std::size_t> first, second;
    DiplotypeScore score;
    double sequence_prior = 0.0, total = 0.0;
    std::string stage = "none";
    std::string tied_groups = "NA", projection_ambiguous = "NA", explored_labellings = "NA";
    std::size_t evaluated = 0;   // fully evaluated labelled diplotypes so far
    bool dirty = false;          // held incumbent newer than the committed checkpoint
    std::size_t commits = 0;
    double seconds_first = std::numeric_limits<double>::quiet_NaN();
    double seconds_commit = std::numeric_limits<double>::quiet_NaN();
    bool failed = false;         // scoring inconsistency or checkpoint publication failure
    std::string failure = "none";
    // the early evaluation
    std::string early_status = "not_run";
    std::string selection = "NA";
    std::size_t pairs_examined = 0, prior_evaluations = 0;
    double seconds_preparation = 0.0, seconds_selection = 0.0, seconds_evaluation = 0.0;
    FragmentScoringCounters early_counters;
    std::uint64_t early_fragment_evaluations = 0;
    std::uint64_t footprint_after_early = 0;
    std::string union_check = "not_reached";   // reproduced | MISMATCH | not_reached
    std::vector<std::size_t> early_first, early_second;   // the early candidate itself (the incumbent moves on)
    DiplotypeScore early_score;
    double early_prior = 0.0;
    // the active panel spelled once per run, for v2 membership of every checkpoint
    bool templates_spelled = false;
    std::vector<std::string> template_sequences;
    std::vector<bool> template_reverse;
};
EarlyIncumbent g_incumbent;
// S2-D: the progressive union schedule of the run in flight (reset at every run start)
struct ProgressiveUnion {
    bool enabled = false, completed = false, prefix_committed = false;
    std::vector<ProgressiveStageReport> stages;
    double seconds_union_start = 0.0, seconds_prefix_checkpoint = std::numeric_limits<double>::quiet_NaN();
    std::size_t prefix_pairs = 0, prefix_tied_groups = 0;
    std::string refusal = "none";
};
ProgressiveUnion g_progressive;
// S1-D.1: true while execute_genotype_command runs the body (its dosage stage follows); a caller-only API run
// (run_genotype, run_genotype_prepared) has no dosage stage, and its v2 report says so rather than naming a mode.
bool g_dosage_stage_follows = false;
// S1-D: classifies the report's haplotypes against the active panel (exact sequence identity in the locus frame) and
// records the graph and model identity -- sequence_first_v2 only; reporting only.
void annotate_v2(SequenceReport& report, const PreparedGenotypeModel& model, const std::string& fingerprint) {
    report.schema = "sequence_first_v2";
    std::vector<std::string> sequences;
    std::vector<bool> reverse;
    const PanelAlleleMatrix& panel = model.panel();
    for (std::size_t t = 0; t < panel.template_count(); ++t) {
        sequences.push_back(spell_mosaic(model, MosaicPath{panel.template_alleles[t]}).bases);
        reverse.push_back(t < model.template_source_orientations().size() &&
                          model.template_source_orientations()[t] == ChainOrientation::Reverse);
    }
    annotate_panel_membership(report, panel.template_names, sequences, reverse);
    for (const auto& item : g_input_hashes) {
        if (item.first == "gfa") report.graph_sha256 = item.second;
    }
    report.model_fingerprint = fingerprint;
}
SequenceLegacy g_sequence_legacy;
bool g_sequence_published = false;   // S1-C: the report was published (or its failure recorded) for this run
std::vector<unsigned char> g_junction_block_flags;   // M11-C: this run's per-block JP1.1 eligibility (empty: not evaluated)
std::uint64_t g_union_end_footprint = 0;   // I2-K (measurement): physical footprint right after the union pass
// I3-A: what the status says about a prepared locus this run used (reset at every run start).
struct PreparedLocusUse {
    bool used = false;
    std::uint64_t id = 0;
    std::string key_sha256;
    std::string origin, file_checksum;   // I3-C: "prepared" or "loaded"; a loaded index's header SHA-256
    std::string panel_source;            // I3-D: "prepared", "loaded" or "rebuilt"
};
PreparedLocusUse g_prepared_use;

// M15-D.1: `exit_code` is the process exit status of a terminal state (stop or failure), written beside the sequence
// fields so the status, .sequence_locus.tsv and the process agree; -1 for the RUNNING record.
std::string minimal_status(const RunWatch& watch, const std::string& state, const std::string& note, int exit_code = -1) {
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
    text << footprint_checkpoint_rows();   // S2-G (blocked route only)
    if (exit_code >= 0) text << "exit_code\t" << exit_code << '\n';
    if (watch.early_incumbent) {   // S2-A
        text << "early_incumbent\texperimental_v1\ncheckpoint_committed\t" << (watch.checkpoint_committed ? "yes" : "no")
             << "\ncheckpoint_stage\t" << watch.checkpoint_stage << '\n';
    }
    if (watch.early_incumbent || watch.sequence_v2) {
        if (state == "FAILED") {   // an exception ended the run: an explicit failure, never a budget stop or fallback
            text << "sequence_prediction\tNO_PREDICTION\nsequence_completion\trun_failed\nsequence_exit_code\t2\n"
                    "checkpoint_retained\tno: a failed run keeps no checkpoint as a prediction\n";
        } else if (state != "RUNNING") {
            const bool kept = watch.checkpoint_committed && !watch.checkpoint_failed;
            text << "sequence_prediction\t" << (kept ? "PROVISIONAL" : "NO_PREDICTION") << "\nsequence_completion\t"
                 << (watch.checkpoint_failed ? "report_publication_failed"
                     : kept                  ? "budget_stop"
                                             : "budget_stop_before_first_evaluated_candidate")
                 << "\nsequence_exit_code\t" << (watch.checkpoint_failed ? 2 : kept ? 3 : exit_code >= 0 ? exit_code : 1)
                 << "\ncheckpoint_retained\t"
                 << (kept ? "yes: <prefix>.sequence_* hold the last committed checkpoint (PROVISIONAL; search incomplete)"
                          : "no") << '\n';
        }
    }
    return text.str();
}

// S2-A: the line `key<TAB>...` of a key/value record replaced (appended when absent).
std::string with_key_value(std::string text, const std::string& key, const std::string& value) {
    const std::string head = key + "\t";
    std::size_t at = text.rfind("\n" + head);
    if (at != std::string::npos) {
        ++at;
    } else if (text.rfind(head, 0) == 0) {
        at = 0;
    } else {
        return text + head + value + "\n";
    }
    const std::size_t end = text.find('\n', at);
    return text.replace(at, (end == std::string::npos ? text.size() : end) - at, head + value);
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
    if (g_watch == nullptr) return;
    std::lock_guard<std::mutex> lock(g_watch->mutex);
    if (g_watch->finished) return;
    g_watch->finished = true;
    std::error_code ignored;
    for (const char* suffix : kFinalOutputs) std::filesystem::remove(g_watch->out_prefix + suffix, ignored);
    if (g_watch->early_incumbent && g_watch->checkpoint_committed && !g_watch->checkpoint_failed) {
        // S2-A: the last committed checkpoint is kept as the PROVISIONAL prediction. Only its pre-rendered locus
        // record is rewritten with the stop -- no scoring and no reconstruction here (this may be a memory stop).
        std::string text = g_watch->checkpoint_locus;
        text = with_key_value(text, "completion", "budget_stop");
        text = with_key_value(text, "completion_note", state + ": " + note +
                                                           "; the last committed checkpoint is kept (search incomplete)");
        text = with_key_value(text, "legacy_run_status", state);
        text = with_key_value(text, "checkpoint_search_status", "budget_stop:" + state);
        code = 3;
        text = with_key_value(text, "exit_code", "3");   // M15-D.1: the record states the process exit
        atomic_write(g_watch->out_prefix + ".sequence_locus.tsv", text);
    } else if (!(g_watch->early_incumbent && g_watch->checkpoint_failed)) {
        // S1-C: the sequence-first reports are run-owned too; a stopped run never sits beside an earlier prediction
        for (const std::string& suffix : sequence_report_outputs()) std::filesystem::remove(g_watch->out_prefix + suffix, ignored);
        if (g_watch->early_incumbent || g_watch->sequence_v2) {   // S2-A: the budget expired before any candidate was fully evaluated
            // M15-D.1: exit_code is the process exit (124 time, 125 footprint), never a no-call's 1
            atomic_write(g_watch->out_prefix + ".sequence_locus.tsv",
                         "key\tvalue\nschema_version\tsequence_first_v2\nprediction\tNO_PREDICTION\nsequence_status\tNA\n"
                         "completion\tbudget_stop_before_first_evaluated_candidate\ncompletion_note\t" + state + ": " + note +
                         "\nlegacy_run_status\t" + state + "\ncheckpoint_stage\tnone\nexit_code\t" + std::to_string(code) + "\n");
        }
    } else {
        code = 2;   // S2-A: a checkpoint publication failed earlier; its failure record stays, never a fallback
    }
    std::filesystem::remove(g_watch->out_prefix + ".sorted.gfa", ignored);
    atomic_write(g_watch->out_prefix + ".mosaic_status.tsv", minimal_status(*g_watch, state, note, code));
    if (g_watch->sequence_v2) {   // M15-G: the stop record carries the primary view too (provisional or no prediction)
        try {
            assessment::write_assessment(g_watch->out_prefix);
        } catch (...) {
        }
    }
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
    // I2-C (measurement only): disjoint wall sub-stages, reported under the timing_detail key as seconds_detail_<name>
    std::vector<std::pair<std::string, double>> detail;
};

// I2-H: the profile key's value names.
const char* gap_recursion_backend_name(bool fixed_geometry) {
    return fixed_geometry ? "fixed_geometry_experimental" : "scaled_reference";
}

GapRecursionBackend gap_recursion_backend_of(const GenotypeCommandOptions& options) {
    return options.gap_recursion_fixed_geometry ? GapRecursionBackend::fixed_geometry_experimental
                                                : GapRecursionBackend::scaled_reference;
}

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

} // namespace

// A declared library (--library FILE), parsed strictly: the six required keys and the optional marker_background, each
// exactly once, one `key value` per line, a value optionally followed by a '#' comment. Blank lines and full-line '#'
// comments are skipped. Everything else is refused with the file, the line and the key, so a misspelt or repeated key
// can never silently fall back to another value. The gap model is not declarable.
LibraryModel read_library(const std::string& path, std::optional<double>& marker_background) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open library model file: " + path);
    static const std::vector<std::string> kRequired = {"insert_mean_bp", "insert_sd_bp", "support_sigmas",
                                                       "substitution_error", "outlier_weight", "fragment_start_intensity"};
    static const std::string kOptional = "marker_background";
    const auto accepted = [&] {
        std::string list;
        for (const std::string& k : kRequired) list += k + ", ";
        return list + kOptional + " (optional)";
    };
    std::map<std::string, std::pair<double, std::size_t>> values;   // key -> (value, line)
    std::string line;
    std::size_t number = 0;
    while (std::getline(input, line)) {
        ++number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const std::string where = path + ":" + std::to_string(number) + ": ";
        std::istringstream fields(line);
        std::string key, text, extra;
        if (!(fields >> key) || key[0] == '#') continue;   // blank line or full-line comment
        if (!(fields >> text)) throw std::runtime_error(where + "expected 'key value', found only '" + key + "'");
        if (fields >> extra && extra[0] != '#') {
            throw std::runtime_error(where + "unexpected content after the value of '" + key + "': '" + extra + "'");
        }
        const bool known = key == kOptional || std::find(kRequired.begin(), kRequired.end(), key) != kRequired.end();
        if (!known) throw std::runtime_error(where + "unknown library key '" + key + "' (accepted: " + accepted() + ")");
        const auto seen = values.find(key);
        if (seen != values.end()) {
            throw std::runtime_error(where + "library key '" + key + "' is repeated (first on line " +
                                     std::to_string(seen->second.second) + ")");
        }
        double value = 0.0;
        std::size_t used = 0;
        try {
            value = std::stod(text, &used);
        } catch (const std::exception&) {
            used = 0;
        }
        if (used == 0 || used != text.size() || !std::isfinite(value)) {
            throw std::runtime_error(where + "library key '" + key + "' is not a finite number: '" + text + "'");
        }
        values[key] = {value, number};
    }
    const auto take = [&](const std::string& key) {
        const auto found = values.find(key);
        if (found == values.end()) throw std::runtime_error(path + ": library model is missing '" + key + "'");
        return found->second.first;
    };
    LibraryModel model;
    model.insert_mean_bp = take("insert_mean_bp");
    model.insert_sd_bp = take("insert_sd_bp");
    model.support_sigmas = take("support_sigmas");
    model.substitution_error = take("substitution_error");
    model.outlier_weight = take("outlier_weight");
    model.fragment_start_intensity = take("fragment_start_intensity");
    model.provenance = "declared: " + path;
    try {
        validate_library_model(model);
    } catch (const std::exception& error) {
        throw std::runtime_error(path + ": invalid library model: " + error.what());
    }
    // Optional. When absent the marker background is estimated jointly with the per-copy rate from
    // the invariant anchors, which is possible whenever those anchors carry two distinct copy numbers.
    const auto background = values.find(kOptional);
    if (background != values.end()) marker_background = background->second.first;
    return model;
}

namespace {

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
    // Exact union placement reuse, as the union pass reported it (or as far as it got before a refusal).
    UnionPlacementCacheReport union_placement_cache;
    UnionOriginCacheReport union_origin_cache;
    UnionPrefilterCacheReport union_prefilter_cache;
    // I2-A: set when the initial union pass ran on the experimental indexed backend (whether it completed or refused)
    bool indexed_union = false;
    IndexedUnionReport indexed;
    // S2-G: the blocked union route (union_fragment_blocks), when it ran: its result without the scores (blocks,
    // completed and refused-block work, totals), the block size, the caps it ran under and its largest per-block footprint
    bool blocked_route = false;
    std::size_t blocked_block_pairs = 0;
    BlockedUnionResult blocked;
    FragmentScoringOptions blocked_options;
    std::uint64_t blocked_footprint_max = 0;
    StartSortReport start_sort;
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
    // M3-E (EXPERIMENTAL): full-data-likelihood block support (empty unless block_support_rule is experimental and the
    // summary path computed it); competitor paths for reporting, by group index.
    DataSupportResult data_support;
    DataSupportResultV2 data_support_v2;   // M3-E.1: block_support_rule full_data_likelihood_v2_experimental
    std::vector<std::string> data_support_competitor_paths;
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
    // The combined block quality of every block row written, in block order (empty when no table was written).
    std::vector<std::string> quality_block_ids;
    std::vector<BlockQuality> block_qualities;
};

// Records one block's combined quality in the summary and returns its two table columns.
std::string quality_columns(RankingSummary& summary, const std::string& block_id, const std::string& representation_filter,
                            const std::string& genotype_support, bool library_sensitive = false) {
    const BlockQuality quality = block_quality(representation_filter, genotype_support, library_sensitive);
    summary.quality_block_ids.push_back(block_id);
    summary.block_qualities.push_back(quality);
    return "\t" + quality.filter + "\t" + quality.reason;
}

// M3-E (EXPERIMENTAL): the versioned support columns appended to a block table (block_support_rule experimental only).
const char* const kDataSupportHeader =
    "\tgenotype_support_v1x\tquality_filter_v1x\tquality_reason_v1x\tdata_support_margin_v1x\tdata_support_certified_v1x\t"
    "competitor_group_v1x\tcompetitor_paths_v1x\tcompetitor_delta_fragment_v1x\tcompetitor_delta_exposure_v1x\t"
    "competitor_delta_likelihood_v1x\tcompetitor_delta_prior_v1x\tdisagreeing_groups_v1x\tmixed_groups_v1x";
const char* const kDataSupportHeaderV2 =
    "\tgenotype_support_v2x\tquality_filter_v2x\tquality_reason_v2x\tsupport_certified_v2x\tsupport_gap_v2x\tpoint_margin_v2x\t"
    "mixed_upper_minus_same_lower_v2x\tcompetitor_group_v2x\tcompetitor_paths_v2x\tcompetitor_delta_fragment_v2x\t"
    "competitor_delta_exposure_v2x\tcompetitor_delta_likelihood_v2x\tcompetitor_delta_prior_v2x\tdisagreeing_groups_v2x\tmixed_groups_v2x";
std::string data_support_header(const std::string& rule) {
    return rule == "full_data_likelihood_v1_experimental" ? kDataSupportHeader
           : rule == "full_data_likelihood_v2_experimental" ? kDataSupportHeaderV2
                                                             : "";
}
// v2 (M3-E.1): envelope certification; quality PASS / LOW_QUALITY / UNASSESSED (an absent comparison is never PASS).
std::string data_support_columns_v2(const RankingSummary& summary, std::size_t b, const std::string& representation_filter,
                                    bool library_sensitive, bool evaluated) {
    if (!evaluated || b >= summary.data_support_v2.blocks.size()) {
        return "\tNOT_EVALUATED\tNOT_EVALUATED\tsupport_not_evaluated\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
    }
    const BlockDataSupportV2& s = summary.data_support_v2.blocks[b];
    std::vector<std::string> failing;
    if (representation_filter == "LOWQUAL_REPRESENTATION_CONFLICT") failing.push_back("representation_conflict");
    else if (representation_filter != "PASS") failing.push_back("representation_not_evaluated");
    const std::string term = data_support_quality_reason_v2(s.label);
    const bool unassessed_term = term == "no_evaluated_alternative";
    if (!term.empty() && !unassessed_term) failing.push_back(term);
    if (library_sensitive) failing.push_back("library_sensitive");
    std::string reason;
    for (const std::string& f : failing) reason += (reason.empty() ? "" : ",") + f;
    if (unassessed_term) reason += (reason.empty() ? "" : ",") + term;
    const std::string filter = !failing.empty() ? "LOW_QUALITY" : unassessed_term ? "UNASSESSED" : "PASS";
    std::ostringstream out;
    out << std::setprecision(12) << '\t' << s.label << '\t' << filter << '\t' << (reason.empty() ? std::string("none") : reason) << '\t'
        << (s.label == "NO_EVALUATED_ALTERNATIVE" || s.label == "AMBIGUOUS_PROJECTION" ? "NA" : (s.certified ? "true" : "false"));
    const auto num = [&out](double v) {
        out << '\t';
        if (std::isnan(v)) out << "NA"; else out << v;
    };
    num(s.support_gap);
    num(s.point_margin);
    num(s.mixed_upper_minus_same_lower);
    out << '\t' << (s.competitor < 0 ? std::string("NA") : std::to_string(s.competitor)) << '\t'
        << (s.competitor < 0 ? std::string("NA") : summary.data_support_competitor_paths.at(static_cast<std::size_t>(s.competitor)));
    num(s.delta_fragment);
    num(s.delta_exposure);
    num(s.delta_likelihood);
    num(s.delta_prior);
    out << '\t' << s.disagreeing << '\t' << s.mixed;
    return out.str();
}
// The legacy quality combination with only the genotype-support term replaced by the experimental label's term.
std::string data_support_columns(const std::string& rule, const RankingSummary& summary, std::size_t b,
                                 const std::string& representation_filter, bool library_sensitive, bool evaluated) {
    if (rule == "full_data_likelihood_v2_experimental") {
        return data_support_columns_v2(summary, b, representation_filter, library_sensitive, evaluated);
    }
    if (!evaluated || b >= summary.data_support.blocks.size()) {
        return "\tNOT_EVALUATED\tNOT_EVALUATED\tsupport_not_evaluated\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA";
    }
    const BlockDataSupport& s = summary.data_support.blocks[b];
    std::string reason;
    const auto add = [&reason](const std::string& item) { reason += (reason.empty() ? "" : ",") + item; };
    if (representation_filter == "LOWQUAL_REPRESENTATION_CONFLICT") add("representation_conflict");
    else if (representation_filter != "PASS") add("representation_not_evaluated");
    if (!data_support_quality_reason(s.label).empty()) add(data_support_quality_reason(s.label));
    if (library_sensitive) add("library_sensitive");
    std::ostringstream out;
    out << std::setprecision(12) << '\t' << s.label << '\t' << (reason.empty() ? "PASS" : "LOW_QUALITY") << '\t'
        << (reason.empty() ? std::string("none") : reason) << '\t';
    const auto num = [&out](double v) {
        if (std::isnan(v)) out << "NA"; else out << v;
    };
    num(s.margin);
    out << '\t' << (s.competitor < 0 ? "NA" : (s.certified ? "true" : "false")) << '\t'
        << (s.competitor < 0 ? std::string("NA") : std::to_string(s.competitor)) << '\t'
        << (s.competitor < 0 ? std::string("NA") : summary.data_support_competitor_paths.at(static_cast<std::size_t>(s.competitor))) << '\t';
    num(s.delta_fragment);
    out << '\t';
    num(s.delta_exposure);
    out << '\t';
    num(s.delta_likelihood);
    out << '\t';
    num(s.delta_prior);
    out << '\t' << s.disagreeing << '\t' << s.mixed;
    return out.str();
}

// What the panel-anchored adaptive completion did, for the status. Every fact separately: the anchor,
// the search and why it stopped, the pool, and the decision made over it.
struct RepairReport {
    const GenotypeCommandOptions* options = nullptr;
    // S2-J: the production search's origin-admission account (reported when the experimental policy is on)
    bool admission_reported = false;
    OriginAdmissionReport admission;
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
    // v2 tie branching: what ran, and on what the call is conditional.
    std::string search_conditional = "none";
    bool tie_branched = false;
    std::size_t tie_ties = 0;
    std::size_t tie_branches_run = 0;
    std::size_t tie_distinct_terminal_genotypes = 0;
    std::string tie_incomplete_reason = "none";
    std::uint64_t tie_branch_physical_origins = 0;
    std::uint64_t pool_branch_scored = 0;
    std::uint64_t tie_digest_memo_hits = 0, tie_digest_memo_misses = 0;
    double seconds_tie_branching = 0.0;
    // M3-D (EXPERIMENTAL): the joint-block pass and its single continuation leg ("off" when the option is off).
    std::string joint_search_status = "off";
    JointSearchResult joint;
    // M9-G (EXPERIMENTAL): junction-aware proposals, summed over every search of the run (trunk, legs, tie branches).
    JunctionProposalReport junction;
    std::size_t junction_read_kmers = 0, junction_read_bytes = 0;
    double junction_rate = 0.0, junction_background = 0.0, junction_table_seconds = 0.0;
    bool joint_continued = false;
    std::size_t joint_leg_first_round = 0, joint_leg_rounds = 0, joint_leg_moves = 0;
    std::string joint_leg_stop = "none";
    bool joint_leg_converged = false;
    std::uint64_t pool_joint_scored = 0, pool_joint_leg_scored = 0;
    std::uint64_t joint_physical_origins = 0, joint_state_classes = 0, joint_leg_physical_origins = 0;
    std::uint64_t joint_footprint_peak = 0;
    double seconds_joint_search = 0.0, seconds_joint_leg = 0.0;
    // Branch seeds: materialised just in time (production) or eagerly at their tie (a gate hook). The work is
    // inside the cumulative counters and caps; these isolate it.
    std::uint64_t tie_seeds_materialised = 0;
    std::uint64_t tie_seed_haplotypes_placed = 0;
    std::uint64_t tie_seed_physical_origins = 0;
    std::uint64_t tie_seed_state_classes = 0;
    std::uint64_t tie_seed_piece_hits = 0;
    std::uint64_t tie_seed_latent_pairs = 0;
    double seconds_tie_seed_placement = 0.0;
    std::size_t tie_seed_peak_live_censuses = 0;
    std::size_t tie_peak_pending_descriptors = 0;
    std::uint64_t tie_seed_cache_built = 0, tie_seed_cache_admitted = 0, tie_seed_cache_rejected = 0;
    std::uint64_t tie_seed_cache_hits = 0, tie_seed_cache_physical_origins = 0, tie_seed_cache_peak_bytes = 0;
    // I4-C: a rejected seed's own new slots are released at once (their work stays charged); the bytes given back, the
    // most held at once counting a just-built seed before its admission decision, and the charged seed metadata
    std::uint64_t tie_seed_cache_rejected_released_bytes = 0, tie_seed_cache_transient_peak_bytes = 0;
    std::uint64_t tie_seed_cache_metadata_bytes = 0;
    double seconds_tie_seed_cache_build = 0.0;
    // The representation-conflict guard, over the called blocks (NA for a no-call: never relabelled LowQual).
    std::string locus_filter = "NA";
    std::string lowqual_blocks = "none";
    std::string not_evaluated_blocks = "none";
    std::string representation_guard_complete = "NA";   // NA for a no-call: the guard never runs on one
    std::size_t representation_competitors_scored = 0;
    std::uint64_t representation_guard_physical_origins = 0;
    double seconds_representation_guard = 0.0;
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
    FineWindowReport fine_window;   // I2-B: the adaptive search context's window enumeration (either path)
    // I2-C (measurement only): the adaptive search's wall, split without overlap
    bool timed = false;
    double seconds_pre_search = 0.0, seconds_search_construction = 0.0, seconds_search_run = 0.0, seconds_post_search = 0.0;
    RepairStageSeconds stages;
    std::uint64_t full_placements = 0;   // I2-D: full haplotype placements by the search's context
    std::size_t enumeration_workers = 0;
    double seconds = 0.0;
    double seconds_scoring = 0.0;
    // library_sensitivity_audit_v1: the rendered record (status lines), empty when the audit did not run.
    std::string library_sensitivity_text;
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
        << "\nrepair_max_neighbours\t" << (o.repair_max_neighbours_model ? std::string("model") : std::to_string(o.repair_max_neighbours))
        << "\nrepair_neighbour_bound\t" << r.neighbour_bound
        << "\nrepair_max_physical_origins\t" << o.repair_max_physical_origins << "\nrepair_max_state_classes\t"
        << o.repair_max_state_classes << "\nrepair_max_piece_hits\t" << o.repair_max_piece_hits
        << "\nrepair_max_latent_pairs\t" << o.repair_max_latent_pairs << '\n';
    if (o.setting_sources.count("adaptive_fused_seed") != 0) {
        out << "adaptive_fused_seed\t" << (o.adaptive_fused_seed ? "1" : "0") << "\nadaptive_full_placements\t"
            << r.full_placements << '\n';
    }
    if (o.timing_detail && r.timed) {
        // I2-C (measurement only): the adaptive search's wall without overlap (all wall seconds; worker-elapsed figures
        // carry the worker_seconds_ prefix elsewhere). pre_search: from the stage's start to the search's construction
        // (neighbourhood bound, the union pool); construction: the search and its scoring context (fragment preparation);
        // run: the search's rounds, split below; post_search: from the search's return to the status, less tie branching
        // and the representation guard (reported apart): the pooled decision, phase and outputs.
        const RepairStageSeconds& g = r.stages;
        out << std::fixed << std::setprecision(3)
            << "seconds_detail_adaptive_pre_search\t" << r.seconds_pre_search
            << "\nseconds_detail_adaptive_search_construction\t" << r.seconds_search_construction
            << "\nseconds_detail_adaptive_search_run\t" << r.seconds_search_run
            << "\nseconds_detail_adaptive_post_search\t" << r.seconds_post_search
            << "\nseconds_detail_adaptive_seed_placement_census\t" << g.seed_placement_census
            << "\nseconds_detail_adaptive_incumbent_spelling\t" << g.incumbent_spelling
            << "\nseconds_detail_adaptive_neighbour_generation\t" << g.neighbour_generation
            << "\nseconds_detail_adaptive_neighbour_spelling\t" << g.neighbour_spelling
            << "\nseconds_detail_adaptive_neighbour_priors\t" << g.neighbour_priors
            << "\nseconds_detail_adaptive_scoring\t" << g.scoring
            << "\nseconds_detail_adaptive_scored_pair_priors\t" << g.scored_pair_priors
            << "\nseconds_detail_adaptive_pool_hook\t" << g.pool_hook
            << "\nseconds_detail_adaptive_certification_and_move\t" << g.certification_and_move
            << "\nseconds_detail_adaptive_promotion_release\t" << g.promotion_release
            << "\nseconds_detail_adaptive_memory_snapshots\t" << g.memory_snapshots
            << "\nseconds_detail_adaptive_trajectory_hook\t" << g.trajectory_hook
            << "\nseconds_detail_adaptive_rounds_wall\t" << g.rounds_wall
            << "\nseconds_detail_adaptive_round_unattributed\t" << g.round_remainder()
            << "\nseconds_detail_adaptive_run_outside_rounds\t" << r.seconds_search_run - g.rounds_wall << '\n';
        out.unsetf(std::ios::floatfield);
        out << std::setprecision(17);
    }
    if (o.setting_sources.count("adaptive_window_discovery") != 0) {
        // I2-B (EXPERIMENTAL): the adaptive search's window enumeration, reported whenever the profile named its path
        // (either value), so a control run and a fine run report the same split. Provenance, work, memory and time only.
        const FineWindowReport& w = r.fine_window;
        out << "adaptive_window_discovery\t" << o.adaptive_window_discovery
            << "\nadaptive_window_discovery_scope\tadd_repaired_haplotypes_rounds"
            << "\nadaptive_window_enumerations\t" << w.enumerations
            << "\nworker_seconds_adaptive_window_rows\t" << w.seconds_rows
            << "\nworker_seconds_adaptive_window_pairing\t" << w.seconds_pairing
            << "\nadaptive_fine_store_budget_bytes\t" << w.budget_bytes
            << "\nadaptive_fine_max_served_placements\t" << w.served_cap
            << "\nadaptive_fine_rounds\t" << w.rounds
            << "\nadaptive_fine_worker_windows\t" << w.worker_windows
            << "\nadaptive_fine_fallback_windows\t" << w.fallback_windows
            << "\nadaptive_fine_units_discovered\t" << w.units_discovered
            << "\nadaptive_fine_units_served\t" << w.units_served
            << "\nadaptive_fine_units_not_retained\t" << w.units_not_retained
            << "\nadaptive_fine_placements_discovered\t" << w.placements_discovered
            << "\nadaptive_fine_placements_served\t" << w.placements_served
            << "\nadaptive_fine_hash_collisions\t" << w.hash_collisions
            << "\nadaptive_fine_round_peak_bytes\t" << w.round_peak_bytes
            << "\nadaptive_fine_round_peak_held_bytes\t" << w.round_peak_held_bytes
            << "\nadaptive_fine_round_peak_temporary_bytes\t" << w.round_peak_temporary_bytes
            << "\nadaptive_fine_ledger_scope\tcontainers_only_not_process_memory\n";
    }
    out << "call_class\t" << r.call_class << "\nno_call_reason\t" << r.no_call_reason << '\n';
    if (o.repair_junction_proposals != "off") {   // M9-G (EXPERIMENTAL): recorded only when the experiment is on
        out << "repair_junction_proposals\t" << o.repair_junction_proposals
            << "\njunction_proposal_rule\tJP1.1_all_canonical_31mers_incumbent_windows_poisson_content_tiebreak"
            << "\njunction_proposal_scope\tsearch_heuristic_not_likelihood_not_confidence_not_final_score"
            << "\njunction_eligible_blocks\t" << r.junction.eligible_blocks << "\njunction_invocations\t" << r.junction.invocations
            << "\njunction_cache_hits\t" << r.junction.cache_hits << "\njunction_cutoff_tie_crossings\t" << r.junction.cutoff_tie_crossings
            << "\njunction_seconds\t" << r.junction.seconds << "\njunction_read_kmers\t" << r.junction_read_kmers
            << "\njunction_read_table_bytes\t" << r.junction_read_bytes << "\njunction_read_table_seconds\t" << r.junction_table_seconds
            << "\njunction_rate_per_copy\t" << hex_double(r.junction_rate) << "\njunction_background\t" << hex_double(r.junction_background)
            << "\njunction_rate_source\tthe_run_marker_depth_fit_unchanged;no_runtime_all_kmer_rate_check\n";
    }
    if (o.repair_joint_search != "off") {   // M3-D (EXPERIMENTAL): recorded only when the experiment is on
        const JointSearchResult& j = r.joint;
        out << "repair_joint_search\t" << o.repair_joint_search << "\njoint_search_scope\tbounded_heuristic_neighbourhood_not_exhaustive"
            << "\njoint_search_status\t" << r.joint_search_status << "\njoint_search_decision\t" << j.decision
            << "\njoint_search_top\t4\njoint_search_single_limit\t2048\njoint_search_catalogue_limit\t2048"
            << "\njoint_search_single_alternatives\t" << j.bounds.single_alternatives << "\njoint_search_catalogue_raw\t"
            << j.bounds.catalogue_raw << "\njoint_search_bounds_within\t" << (j.bounds.within ? "true" : "false")
            << "\njoint_search_singles_scored\t" << j.singles_scored << "\njoint_search_promoted\t" << j.promoted
            << "\njoint_search_joint_requests\t" << j.joint_requests << "\njoint_search_catalogue_unique\t" << j.catalogue.size()
            << "\njoint_search_winner_set_size\t" << j.winner_set_size << "\njoint_search_incumbent_in_top_set\t"
            << (j.incumbent_in_top_set ? "true" : "false") << "\njoint_search_certified_gap\t" << j.certified_gap
            << "\njoint_search_moved\t" << (j.moved ? "true" : "false") << "\njoint_search_move\t"
            << (j.moved ? j.catalogue[j.winner].move : std::string("none")) << "\njoint_search_incumbent_total\t"
            << hex_double(j.incumbent_total) << "\njoint_search_winner_total\t"
            << (j.moved ? hex_double(j.catalogue[j.winner].score.log_likelihood + j.catalogue[j.winner].prior) : std::string("NA"))
            << "\njoint_search_origins_reused\t" << j.origins_reused << "\njoint_search_origins_recomputed\t" << j.origins_recomputed
            << "\njoint_search_dp_executions\t" << j.dp_executions << "\njoint_search_window_bases\t" << j.window_bases
            << "\njoint_search_full_placements\t" << j.full_placements << "\njoint_search_census_bytes_peak\t" << j.census_bytes_peak
            << "\njoint_search_transient_bytes_peak\t" << j.transient_bytes_peak << "\njoint_search_transient_bytes_after\t"
            << j.transient_bytes_after << "\njoint_search_slots_created\t" << j.slots_created << "\njoint_search_live_transient_after\t"
            << j.live_transient_after << "\njoint_search_physical_origins\t" << r.joint_physical_origins
            << "\njoint_search_state_classes\t" << r.joint_state_classes << "\njoint_search_footprint_peak\t" << r.joint_footprint_peak
            << "\npool_joint_scored\t" << r.pool_joint_scored << "\njoint_continued\t" << (r.joint_continued ? "true" : "false")
            << "\njoint_leg_first_round\t" << r.joint_leg_first_round << "\njoint_leg_rounds\t" << r.joint_leg_rounds
            << "\njoint_leg_moves\t" << r.joint_leg_moves << "\njoint_leg_stop\t" << r.joint_leg_stop << "\njoint_leg_converged\t"
            << (r.joint_leg_converged ? "true" : "false") << "\npool_joint_leg_scored\t" << r.pool_joint_leg_scored
            << "\njoint_leg_physical_origins\t" << r.joint_leg_physical_origins << "\nseconds_joint_search\t" << r.seconds_joint_search
            << "\nseconds_joint_search_singles\t" << j.seconds_singles << "\nseconds_joint_search_promotion\t" << j.seconds_promotion
            << "\nseconds_joint_search_joint\t" << j.seconds_joint << "\nseconds_joint_search_scoring\t" << j.seconds_scoring
            << "\nseconds_joint_search_priors\t" << j.seconds_priors << "\nseconds_joint_leg\t" << r.seconds_joint_leg << '\n';
    }
    out << "search_conditional\t" << r.search_conditional << "\nrepair_tie_branching\t"
        << (o.repair_tie_branching ? "true" : "false") << "\nrepair_tie_max_branches\t" << o.repair_tie_max_branches
        << "\nrepair_tie_seed_cache_bytes\t" << o.repair_tie_seed_cache_bytes
        << "\ntie_branch_order\t" << kTieBranchOrder << "\ntie_branch_seeding\t"
        << (o.test_eager_branch_seeds || o.test_retain_branch_parents ? "eager_delta_promoted_shared_context"
            : o.repair_tie_seed_cache_bytes > 0 ? "bounded_seed_cache_then_just_in_time_shared_context"
                                                : "census_free_just_in_time_shared_context")
        << "\ntie_branched\t" << (r.tie_branched ? "true" : "false") << "\ntie_ties\t" << r.tie_ties
        << "\ntie_branches_run\t" << r.tie_branches_run << "\ntie_distinct_terminal_genotypes\t"
        << r.tie_distinct_terminal_genotypes << "\ntie_incomplete_reason\t" << r.tie_incomplete_reason
        << "\ntie_branch_physical_origins\t" << r.tie_branch_physical_origins << "\npool_branch_scored\t"
        << r.pool_branch_scored << "\ntie_digest_memo_hits\t" << r.tie_digest_memo_hits << "\ntie_digest_memo_misses\t"
        << r.tie_digest_memo_misses << "\ntie_seeds_materialised\t" << r.tie_seeds_materialised
        << "\ntie_seed_haplotypes_placed\t" << r.tie_seed_haplotypes_placed << "\ntie_seed_physical_origins\t"
        << r.tie_seed_physical_origins << "\ntie_seed_state_classes\t" << r.tie_seed_state_classes
        << "\ntie_seed_piece_hits\t" << r.tie_seed_piece_hits << "\ntie_seed_latent_pairs\t" << r.tie_seed_latent_pairs
        << "\nseconds_tie_seed_placement\t" << r.seconds_tie_seed_placement << "\ntie_seed_peak_live_censuses\t"
        << r.tie_seed_peak_live_censuses << "\ntie_peak_pending_descriptors\t" << r.tie_peak_pending_descriptors
        << "\ntie_seed_cache_built\t" << r.tie_seed_cache_built << "\ntie_seed_cache_admitted\t" << r.tie_seed_cache_admitted
        << "\ntie_seed_cache_rejected\t" << r.tie_seed_cache_rejected << "\ntie_seed_cache_hits\t" << r.tie_seed_cache_hits
        << "\ntie_seed_cache_physical_origins\t" << r.tie_seed_cache_physical_origins << "\ntie_seed_cache_peak_bytes\t"
        << r.tie_seed_cache_peak_bytes << "\nseconds_tie_seed_cache_build\t" << r.seconds_tie_seed_cache_build
        << "\ntie_seed_cache_rejected_released_bytes\t" << r.tie_seed_cache_rejected_released_bytes
        << "\ntie_seed_cache_transient_peak_bytes\t" << r.tie_seed_cache_transient_peak_bytes
        << "\ntie_seed_cache_metadata_bytes\t" << r.tie_seed_cache_metadata_bytes << '\n';
    out << "representation_guard\t" << (o.representation_guard ? "true" : "false") << "\nrepresentation_guard_rule\t"
        << kRepresentationRule << "\nrepresentation_guard_bound_nats\t" << kRepresentationBoundNats << "\nlocus_filter\t"
        << r.locus_filter << "\nlowqual_blocks\t" << r.lowqual_blocks << "\nnot_evaluated_blocks\t"
        << r.not_evaluated_blocks << "\nrepresentation_guard_complete\t" << r.representation_guard_complete
        << "\nrepresentation_competitors_scored\t" << r.representation_competitors_scored
        << "\nrepresentation_guard_physical_origins\t" << r.representation_guard_physical_origins << '\n';
    out << "union_ranking_status\t" << r.union_status << "\nunion_winner_set_size\t" << r.union_winner_set
        << "\nunion_certified_gap\t" << r.union_certified_gap << "\nunion_top_score\t" << r.union_top_score
        << "\nunion_top_hap1\t" << r.union_top_hap1 << "\nunion_top_hap2\t" << r.union_top_hap2 << '\n';
    if (r.admission_reported) {   // S2-J
        const OriginAdmissionReport& a = r.admission;
        out << "repair_round1_origin_admission\treuse_only_experimental\nrepair_admission_restricted_plans\t" << a.restricted_plans
            << "\nrepair_admission_planned_windows\t" << a.planned_windows << "\nrepair_admission_singleton_windows\t"
            << a.singleton_windows << "\nrepair_admission_eligible_windows\t" << a.eligible_windows
            << "\nrepair_admission_misses_while_restricted\t" << a.misses_while_restricted
            << "\nrepair_admission_singleton_misses_not_retained\t" << a.singleton_misses_not_retained
            << "\nrepair_admission_retained_while_restricted\t" << a.retained_while_restricted
            << "\nrepair_admission_hits_while_restricted\t" << a.hits_while_restricted
            << "\nrepair_admission_group_reuses_while_restricted\t" << a.group_reuses_while_restricted << '\n';
    }
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
        << r.seconds_scoring << "\nseconds_tie_branching\t" << r.seconds_tie_branching
        << "\nseconds_representation_guard\t" << r.seconds_representation_guard << '\n';
    out.unsetf(std::ios::floatfield);
    out << std::setprecision(6);
}

// The effective value of an experimental profile key (the settings table's own mapping; M15-H: also the built-in
// defaults the advanced help prints, so the two can never disagree).
std::string experimental_setting_value(const GenotypeCommandOptions& options, const std::string& k) {
    return k == "union_scoring_backend"       ? options.union_scoring_backend
                                      : k == "indexed_union_store_bytes" ? std::to_string(options.indexed_union_store_bytes)
                                      : k == "adaptive_window_discovery" ? options.adaptive_window_discovery
                                      : k == "adaptive_fine_store_bytes" ? std::to_string(options.adaptive_fine_store_bytes)
                                      : k == "indexed_union_store_layout" ? options.indexed_union_store_layout
                                      : k == "indexed_union_evidence_layout" ? options.indexed_union_evidence_layout
                                      : k == "timing_detail"              ? std::string(options.timing_detail ? "1" : "0")
                                      : k == "adaptive_fused_seed"        ? std::string(options.adaptive_fused_seed ? "1" : "0")
                                      : k == "marker_proposal_reuse_table_bytes"
                                          ? std::to_string(options.marker_proposal_reuse_table_bytes)
                                      : k == "combination_prepared_library_terms"
                                          ? std::string(options.combination_prepared_library_terms ? "1" : "0")
                                      : k == "dosage_census_backend"
                                          ? std::string(options.dosage_census_interned ? "interned_experimental" : "reference")
                                      : k == "adaptive_emission_handoff"
                                          ? std::string(options.adaptive_emission_handoff ? "1" : "0")
                                      : k == "output_schema" ? options.output_schema
                                      : k == "early_incumbent" ? options.early_incumbent
                                      : k == "test_interrupt_point" ? options.test_interrupt_point
                                      : k == "diagnostic_union_evidence_dump" ? std::string(options.diagnostic_union_evidence_dump ? "1" : "0")
                                      : k == "diagnostic_union_kernel_dump" ? std::string(options.diagnostic_union_kernel_dump ? "1" : "0")
                                      : k == "union_progressive_schedule" ? options.union_progressive_schedule
                                      : k == "union_fragment_blocks" ? options.union_fragment_blocks
                                      : k == "diagnostic_refinement_memory" ? std::string(options.diagnostic_refinement_memory ? "1" : "0")
                                      : k == "diagnostic_refinement_replay" ? options.diagnostic_refinement_replay
                                      : k == "repair_round1_origin_admission" ? options.repair_round1_origin_admission
                                      : k == "repair_joint_search" ? options.repair_joint_search
                                      : k == "repair_junction_proposals" ? options.repair_junction_proposals
                                      : k == "compound_annotation" ? options.compound_annotation
                                      : k == "block_diagnostics" ? options.block_diagnostics
                                      : k == "compound_interval_declarations" ? options.compound_interval_declarations
                                      : k == "block_support_rule" ? options.block_support_rule
                                      : k == "copying_transition" ? options.copying_transition
                                      : k == "copying_length_bp" ? hex_double(options.copying_length_bp)
                                      : k == "gap_recursion_backend"
                                          ? std::string(gap_recursion_backend_name(options.gap_recursion_fixed_geometry))
                                                                         : std::to_string(options.adaptive_fine_max_served_placements);
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
    // Whether emission-cache lookups were timed; when not, no lookup time anywhere is a measurement.
    out << "emission_lookup_timing\t" << (options.emission_lookup_timing ? "measured" : "not_measured") << '\n';
    out << "time_limit_seconds\t" << options.time_limit_seconds << '\n';
    {
        // Worker and process memory at the moment the status is written.
        const FootprintSample footprint = physical_footprint();
        out << "threads_peak_active_workers\t" << peak_active_workers() << '\n';
        out << "resources_phys_footprint_peak_bytes\t" << footprint.peak << '\n';
        out << "resources_rss_peak_bytes\t" << static_cast<std::uint64_t>(peak_rss_gib() * 1073741824.0) << '\n';
    }
    note_footprint(std::string("status_written:") + genotype_run_status_name(status));   // S2-G: completion or refusal
    out << footprint_checkpoint_rows();
    out << "gfa\t" << options.gfa_path << '\n';
    out << "reference_path\t" << options.reference_path << '\n';
    if (!options.preset_name.empty()) {
        std::string keys;
        for (const std::string& k : options.preset_keys) keys += (keys.empty() ? "" : ",") + k;
        out << "preset\t" << options.preset_name << "\npreset_source\t" << (options.preset_explicit ? "flag" : "default")
            << "\npreset_keys\t" << (keys.empty() ? "none" : keys) << '\n';
    }
    // library_sensitivity_audit_v1: the record under --library auto; <prefix>.library_sensitivity.tsv on every run.
    std::string sensitivity = "library_sensitivity\tNOT_APPLICABLE\nlibrary_sensitivity_reason\tdeclared_library\n"
                              "library_sensitivity_status\tNOT_APPLICABLE\nlibrary_sensitivity_filter\tNOT_APPLICABLE\n";
    if (options.library_path == "auto") {
        out << "library_mode\tauto\nlibrary_insert_source\tplacement_estimate\nlibrary_intensity_source\tmarker_fit_estimate\n"
               "library_other_source\tshort_read_preset\n";
        out << "library_invocation\t" << (options.library_omitted ? "omitted_means_auto" : "explicit_auto") << '\n';
        if (repair != nullptr && !repair->library_sensitivity_text.empty()) {
            sensitivity = repair->library_sensitivity_text;
        } else {
            // No completed audit behind this run: never presented as evaluated.
            sensitivity = std::string("library_sensitivity\tNOT_EVALUATED\nlibrary_sensitivity_reason\t"
                                      "LIBRARY_SENSITIVITY_NOT_EVALUATED\nlibrary_sensitivity_status\tNOT_EVALUATED\n"
                                      "library_sensitivity_filter\tNOT_EVALUATED_LIBRARY_SENSITIVITY\n"
                                      "library_sensitivity_incomplete_endpoint\tnone\nlibrary_sensitivity_incomplete\t") +
                          (!options.library_sensitivity_audit ? "audit_off"
                           : status == GenotypeRunStatus::LibraryRefused ? "library_refused"
                           : options.repair_rounds == 0                  ? "not_adaptive_completion"
                           : repair != nullptr                           ? "primary_no_call"
                                                                         : "call_not_reached") +
                          "\n";
        }
        out << sensitivity;
    }
    {
        // Every public setting's effective value and where it came from (flag > profile > preset > default).
        std::string table = "key\tvalue\tsource\n";
        for (const std::string& key : genotype_public_settings()) {
            const auto source = options.setting_sources.find(key);
            table += key + "\t" + genotype_setting_value(options, key) + "\t" +
                     (source == options.setting_sources.end() ? std::string("default") : source->second) + "\n";
        }
        // I2-A: the experimental union backend, recorded whenever a profile set it (absent from a default run's table)
        for (const char* key : {"union_scoring_backend", "indexed_union_store_bytes", "adaptive_window_discovery",
                                "adaptive_fine_store_bytes", "adaptive_fine_max_served_placements",
                                "indexed_union_store_layout", "indexed_union_evidence_layout", "timing_detail", "adaptive_fused_seed",
                                "marker_proposal_reuse_table_bytes", "combination_prepared_library_terms",
                                "gap_recursion_backend",
                                "dosage_census_backend", "adaptive_emission_handoff", "output_schema", "early_incumbent",
                                "test_interrupt_point", "diagnostic_union_evidence_dump", "diagnostic_union_kernel_dump",
                                "union_progressive_schedule", "union_fragment_blocks", "diagnostic_refinement_memory",
                                "diagnostic_refinement_replay", "repair_round1_origin_admission", "repair_joint_search", "block_support_rule",
                                "repair_junction_proposals", "compound_annotation", "compound_interval_declarations",
                                "copying_transition", "copying_length_bp", "block_diagnostics"}) {
            const auto source = options.setting_sources.find(key);
            if (source == options.setting_sources.end()) continue;
            const std::string k = key;
            const std::string value = experimental_setting_value(options, k);
            table += k + "\t" + value + "\t" + source->second + "\n";
        }
        atomic_write(options.out_prefix + ".settings.tsv", table);
    }
    atomic_write(options.out_prefix + ".library_sensitivity.tsv",
                 "key\tvalue\nlibrary_mode\t" + std::string(options.library_path == "auto" ? "auto" : "file") + "\n" +
                     sensitivity);
    out << "kmer_size\t" << options.kmer_size << '\n';
    out << "syncmer_s\t" << options.syncmer_s << '\n';
    out << "min_variant_bp\t" << options.min_variant_bp << '\n';
    out << "switch_probability\t" << options.switch_probability << '\n';
    if (options.setting_sources.count("copying_transition") != 0) {   // M4-B: recorded only when a profile set the key
        const CopyingTransitions& tr = g_transitions;
        out << "copying_transition\t" << options.copying_transition << "\ncopying_redraw_measure\tuniform_template_rows_v1\n";
        if (options.copying_transition != "constant_v1") {
            std::ostringstream v;
            v << std::setprecision(17) << options.copying_length_bp;
            out << "copying_length_bp\t" << v.str() << "\ncopying_length_bp_hex\t" << hex_double(options.copying_length_bp) << '\n';
        }
        out << "locus_geometry_status\t" << g_geometry.status << "\nlocus_geometry_reference\t"
            << (g_geometry.reference_path_name.empty() ? "NA" : g_geometry.reference_path_name) << "\nlocus_geometry_digest\t"
            << (g_geometry.digest.empty() ? "NA" : g_geometry.digest) << "\nlocus_geometry_frame_digest\t"
            << (g_geometry.frame_digest.empty() ? "NA" : g_geometry.frame_digest) << '\n';
        if (!tr.r.empty()) {
            out << "transition_vector_sha256\t" << tr.vector_digest << "\ntransition_boundaries\t" << tr.r.size() << "\ntransition_r_min\t"
                << *std::min_element(tr.r.begin(), tr.r.end()) << "\ntransition_r_max\t" << *std::max_element(tr.r.begin(), tr.r.end())
                << "\ntransition_zero_distance_boundaries\t" << tr.zero_distance_boundaries
                << "\ntransition_zero_distance_informative_boundaries\t" << tr.zero_distance_informative_boundaries << '\n';
            if (tr.distance_mode()) {
                out << "locus_geometry_distances_bp\t";
                for (std::size_t b = 0; b < tr.distance_bp.size(); ++b) {   // exact halves of a bp
                    std::ostringstream v;
                    v << std::fixed << std::setprecision(1) << tr.distance_bp[b];
                    out << (b ? "," : "") << v.str();
                }
                out << '\n';
            }
        } else {
            out << "transition_vector_sha256\tNA\n";
        }
    }
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
        if (!model_report->curated_exclusions.empty()) {
            // The active panel: alleles only curated templates carried are not callable (not in the model).
            out << "curation_inactive_source_alleles\t" << model_report->total_curation_inactive_source_alleles << '\n';
            out << "curation_inactive_canonical_alleles\t" << model_report->total_curation_inactive_canonical_alleles << '\n';
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
        {
            // Exact union placement reuse: what it held, served and consumed. Executed seeding counters above
            // count only executed work; replayed placements are counted here and held to their own cap.
            const UnionPlacementCacheReport& u = candidate_set->union_placement_cache;
            out << "union_placement_cache\t" << (u.active ? "on" : "off") << "\nunion_placement_cache_budget_bytes\t"
                << u.budget_bytes << "\nmax_union_cached_placements_consumed\t" << u.max_placements_consumed
                << "\nunion_placement_cache_read_span_bp\t" << u.read_span_bp << "\nunion_placement_cache_max_gap_bp\t"
                << u.max_gap_bp << "\nunion_placement_cache_segment_occurrences\t" << u.segment_occurrences
                << "\nunion_placement_cache_distinct_segments\t" << u.distinct_segments
                << "\nunion_placement_cache_eligible_segments\t" << u.eligible_segments
                << "\nunion_placement_cache_admitted_segments\t" << u.admitted_segments
                << "\nunion_placement_cache_declined_segments\t" << u.declined_segments
                << "\nunion_placement_cache_evictions\t" << u.evictions << "\nunion_placement_cache_hits\t" << u.hits
                << "\nunion_placement_cache_misses\t" << u.misses << "\nunion_placement_cache_zone_starts_requested\t"
                << u.zone_starts_requested << "\nunion_placement_cache_zone_starts_served\t" << u.zone_starts_served
                << "\nunion_placement_cache_placements_enumerated\t" << u.placements_enumerated
                << "\nunion_placement_cache_placements_cached\t" << u.placements_cached
                << "\nunion_cached_placements_consumed\t" << u.placements_replayed
                << "\nunion_placement_cache_payload_bytes\t" << u.payload_bytes
                << "\nunion_placement_cache_index_bytes\t" << u.index_bytes
                << "\nunion_placement_cache_peak_bytes\t" << u.peak_bytes
                << "\nseconds_union_placement_cache_fill\t" << u.seconds_fill << '\n';
            // Exact union origin reuse: served origins and aggregates are not executed work, so they are counted
            // and capped here, apart from physical_origins and statistical_state_classes.
            const UnionOriginCacheReport& o = candidate_set->union_origin_cache;
            out << "union_origin_cache\t" << (o.active ? "on" : "off") << "\nunion_origin_cache_budget_bytes\t"
                << o.budget_bytes << "\nmax_union_cached_origins_consumed\t" << o.max_origins_consumed
                << "\nmax_union_cached_states_consumed\t" << o.max_states_consumed
                << "\nunion_origin_cache_segment_occurrences\t" << o.segment_occurrences
                << "\nunion_origin_cache_eligible_segments\t" << o.eligible_segments
                << "\nunion_origin_cache_admitted_segments\t" << o.admitted_segments
                << "\nunion_origin_cache_declined_segments\t" << o.declined_segments
                << "\nunion_origin_cache_evictions\t" << o.evictions << "\nunion_origin_cache_hits\t" << o.hits
                << "\nunion_origin_cache_misses\t" << o.misses
                << "\nunion_origin_cache_origins_constructed\t" << o.origins_constructed
                << "\nunion_origin_cache_state_aggregates_constructed\t" << o.state_aggregates_constructed
                << "\nunion_origin_cache_latent_pairs_constructed\t" << o.latent_pairs_constructed
                << "\nunion_origin_cache_valid_states_constructed\t" << o.valid_states_constructed
                << "\nunion_origin_cache_origin_records_held\t" << o.origin_records_held
                << "\nunion_origin_cache_origins_held\t" << o.origins_held
                << "\nunion_cached_origins_consumed\t" << o.origins_consumed
                << "\nunion_cached_states_consumed\t" << o.states_consumed
                << "\nunion_origin_cache_payload_bytes\t" << o.payload_bytes
                << "\nunion_origin_cache_index_bytes\t" << o.index_bytes
                << "\nunion_origin_cache_peak_bytes\t" << o.peak_bytes
                << "\nseconds_union_origin_cache_fill\t" << o.seconds_fill << '\n';
            // Exact union prefilter reuse: decisions served instead of a kernel are counted and capped here; the
            // logical prefilter counters above are unchanged, emission_cache_executed_prefilters counts kernels run.
            const UnionPrefilterCacheReport& f = candidate_set->union_prefilter_cache;
            out << "union_prefilter_cache\t" << (f.active ? "on" : "off") << "\nunion_prefilter_cache_budget_bytes\t"
                << f.budget_bytes << "\nmax_union_cached_prefilters_consumed\t" << f.max_consumed
                << "\nunion_prefilter_cache_read_length\t" << f.read_length << "\nunion_prefilter_cache_window_length\t"
                << f.window_length << "\nunion_prefilter_cache_window_contents\t" << f.window_contents
                << "\nunion_prefilter_cache_window_index_bytes\t" << f.window_index_bytes
                << "\nunion_prefilter_cache_positions_indexed\t" << f.positions_indexed
                << "\nunion_prefilter_cache_decision_table_bytes\t" << f.decision_table_bytes
                << "\nunion_prefilter_cache_decisions_held\t" << f.decisions_held
                << "\nunion_prefilter_cache_decisions_not_retained\t" << f.decisions_not_retained
                << "\nunion_cached_prefilters_consumed\t" << f.hits << "\nunion_prefilter_cache_misses\t" << f.misses
                << "\nworker_seconds_union_prefilter_index\t" << f.worker_seconds_index << '\n';
            const StartSortReport& q = candidate_set->start_sort;
            // The share and the peak a worker held depend on the worker count, so they sit with the threads_ keys.
            out << "start_sort_scratch_bytes\t" << q.budget_bytes << "\nthreads_start_sort_workers\t" << q.workers
                << "\nthreads_start_sort_share_bytes\t" << q.share_bytes << "\nstart_sort_out_of_place\t"
                << q.out_of_place << "\nstart_sort_in_place\t" << q.in_place << "\nstart_sort_fallbacks\t"
                << q.fallbacks << "\nthreads_start_sort_peak_worker_bytes\t" << q.peak_worker_bytes << '\n';
        }
        if (candidate_set->indexed_union && candidate_set->blocked_route) {
            // S2-G (EXPERIMENTAL): the blocked route's provenance, actual resource settings, per-block work and admission,
            // and on refusal the completed blocks, the refused block (direct/served breakdown unavailable) and the attempt.
            const IndexedUnionReport& x = candidate_set->indexed;
            const ContentIndexReport& ix = x.index;
            const BlockedUnionResult& r = candidate_set->blocked;
            const FragmentScoringOptions& o = candidate_set->blocked_options;
            const std::ios_base::fmtflags saved_flags = out.flags();
            const std::streamsize saved_precision = out.precision();
            out << "union_scoring_backend\tindexed_fine_experimental\n"
                << "union_scoring_backend_stages\tinitial_union_pass\n"
                << "union_scoring_backend_unchanged_stages\tadaptive_repair_prepared_context,dosage,phase_audit_rerank\n"
                << "indexed_union_refusal\t" << (x.refusal.empty() ? std::string("none") : x.refusal) << '\n'
                << "union_route\tblocked_fine_core_experimental\nunion_fragment_block_pairs\t" << candidate_set->blocked_block_pairs
                << "\nunion_origin_aggregate_policy\tfine_core_partition_experimental\nunion_evidence_layout\tfine_views"
                << "\nunion_store_layout\texact\nunion_origin_cache_bytes\t" << o.union_origin_cache_bytes
                << "\nunion_origin_whole_budget_bytes\t" << o.union_origin_cache_bytes / 2
                << "\nunion_origin_fine_budget_bytes\t" << o.union_origin_cache_bytes - o.union_origin_cache_bytes / 2
                << "\nindexed_union_store_budget_bytes\t" << x.store_budget_bytes
                << "\nunion_cap_physical_origins\t" << o.max_physical_origins << "\nunion_cap_state_classes\t" << o.max_state_classes
                << "\nunion_cap_piece_hits\t" << o.max_piece_hits << "\nunion_cap_latent_alignment_pairs\t" << o.max_latent_alignment_pairs
                << "\nunion_cap_mass_cells\t" << o.max_mass_cells
                << "\nunion_cap_cached_placements_consumed\t" << o.max_union_cached_placements_consumed
                << "\nunion_cap_cached_origins_consumed\t" << o.max_union_cached_origins_consumed
                << "\nunion_cap_cached_states_consumed\t" << o.max_union_cached_states_consumed
                << "\nunion_cap_cached_prefilters_consumed\t" << o.max_union_cached_prefilters_consumed
                << "\nunion_caps_scope\tcumulative_across_blocks"
                << "\nindexed_union_max_read_bp\t" << x.max_read_bp << "\nindexed_union_haplotypes\t" << ix.haplotypes
                << "\nindexed_union_index_metadata_bytes\t" << ix.unit_bytes + ix.occurrence_bytes
                << "\nindexed_union_interior_units\t" << ix.interior_units << "\nindexed_union_window_units\t" << ix.window_units
                << "\nindexed_union_fine_units\t" << ix.fine_units << "\nindexed_union_hash_collisions\t" << ix.hash_collisions
                << "\nunion_fragment_blocks_completed\t" << r.blocks.size()
                << "\nunion_fragment_blocks_run_completed\t" << (r.completed ? "yes" : "no")
                << "\nunion_fragment_block_footprint_max_bytes\t" << candidate_set->blocked_footprint_max
                << "\nunion_end_phys_footprint_bytes\t" << g_union_end_footprint << '\n';
            for (std::size_t i = 0; i < r.blocks.size(); ++i) {
                const BlockedBlockReport& b = r.blocks[i];
                const std::string k = "union_block" + std::to_string(i + 1) + "_";
                out << k << "fragments\t" << b.fragments << '\n' << k << "executed_physical_origins\t" << b.observed_physical_origins
                    << '\n' << k << "direct_physical_origins\t" << b.direct_physical_origins << '\n' << k << "served_physical_origins\t"
                    << b.served_physical_origins << '\n' << k << "whole_origins_constructed\t" << b.origins_constructed << '\n' << k
                    << "fine_origins_constructed\t" << b.fine_origins_constructed << '\n' << k << "whole_entities_admitted\t"
                    << b.origin_entities_admitted << '\n' << k << "whole_entities_declined\t" << b.origin_entities_declined << '\n'
                    << k << "whole_held_bytes\t" << b.origin_held_bytes << '\n' << k << "fine_windows_eligible\t"
                    << b.fine_windows_eligible << '\n' << k << "fine_windows_admitted\t" << b.fine_windows_admitted << '\n' << k
                    << "fine_windows_declined\t" << b.fine_windows_declined << '\n' << k << "fine_charged_bytes\t"
                    << b.fine_charged_bytes << '\n' << k << "fine_served_physical_origins\t" << b.fine_served_physical_origins
                    << '\n' << k << "executed_state_classes\t" << b.state_classes_executed << '\n' << k << "executed_dps\t"
                    << b.executed_dps << '\n' << k << "phys_footprint_after_bytes\t" << b.footprint_after << '\n';
            }
            const auto totals = [&](const char* name, const UnionCumulativeCharges& c) {
                out << name << "physical_origins\t" << c.physical_origins << '\n' << name << "piece_hits\t" << c.piece_hits << '\n'
                    << name << "cached_placements_consumed\t" << c.placements_consumed << '\n' << name
                    << "cached_origins_consumed\t" << c.origins_consumed << '\n' << name << "cached_states_consumed\t"
                    << c.states_consumed << '\n' << name << "cached_prefilters_consumed\t" << c.prefilters_consumed << '\n'
                    << name << "state_classes\t" << c.executed.statistical_state_classes << '\n';
            };
            totals("union_completed_blocks_", r.completed_totals);
            if (!r.completed) {
                out << "union_refused_block\t" << (r.has_refused_block ? std::to_string(r.blocks.size() + 1) : std::string("none")) << '\n';
                if (r.has_refused_block) {
                    const BlockedBlockReport& b = r.refused_block;
                    out << "union_refused_block_stage\t" << b.stage << "\nunion_refused_block_first_fragment\t" << b.first_fragment
                        << "\nunion_refused_block_fragments\t" << b.fragments
                        << "\nunion_refused_block_executed_physical_origins\t" << b.observed_physical_origins
                        << "\nunion_refused_block_piece_hits\t" << b.observed_piece_hits
                        << "\nunion_refused_block_direct_physical_origins\tunavailable"
                        << "\nunion_refused_block_served_physical_origins\tunavailable"
                        << "\nunion_refused_block_cached_origins_consumed\t" << b.origins_consumed
                        << "\nunion_refused_block_cached_states_consumed\t" << b.states_consumed
                        << "\nunion_refused_block_cached_placements_consumed\t" << b.placements_consumed
                        << "\nunion_refused_block_state_classes\t" << b.state_classes_executed
                        << "\nunion_refused_block_whole_origins_constructed\t" << b.origins_constructed
                        << "\nunion_refused_block_fine_origins_constructed\t" << b.fine_origins_constructed << '\n';
                }
                totals("union_attempt_", r.attempt_totals);
            }
            out << std::fixed << std::setprecision(3) << "seconds_indexed_union_index_build\t" << r.seconds_index_build
                << "\nseconds_union_block_combination\t" << r.seconds_combination << '\n';
            out.flags(saved_flags);
            out.precision(saved_precision);
        } else if (candidate_set->indexed_union) {
            // I2-A (EXPERIMENTAL): which stage used the indexed backend, its guards and what it did. Present only
            // when the profile selected it; every key here is backend provenance, executed work, memory or time.
            const IndexedUnionReport& x = candidate_set->indexed;
            const ContentIndexReport& ix = x.index;
            const ContentBatchReport& b = x.batch;
            const std::ios_base::fmtflags saved_flags = out.flags();
            const std::streamsize saved_precision = out.precision();
            out << "union_scoring_backend\tindexed_fine_experimental\n"
                << "union_scoring_backend_stages\tinitial_union_pass\n"
                << "union_scoring_backend_unchanged_stages\tadaptive_repair_prepared_context,dosage,phase_audit_rerank\n"
                << "indexed_union_refusal\t" << (x.refusal.empty() ? std::string("none") : x.refusal) << '\n'
                << "indexed_union_store_budget_bytes\t" << x.store_budget_bytes
                << "\nindexed_union_store_charged_bytes\t" << b.store_charged_bytes
                << "\nindexed_union_store_held_bytes_after_discovery\t" << b.placement_store_bytes
                << "\nindexed_union_store_layout\t" << b.store_layout
                << "\nindexed_union_evidence_layout\t" << b.evidence_layout
                << "\nindexed_union_view_coarse_units\t" << b.view_coarse_units
                << "\nindexed_union_view_rows_represented\t" << b.view_rows_represented
                << "\nindexed_union_view_table_bytes\t" << b.view_table_bytes
                << "\nindexed_union_view_zone_temporary_peak_bytes\t" << b.view_zone_temporary_peak_bytes
                << "\nindexed_union_view_zone_temporary_bytes_total\t" << b.view_zone_temporary_bytes_total
                << "\nseconds_indexed_union_view_zone_gather\t" << b.seconds_view_zone_gather
                << "\nindexed_union_store_budget_basis\t" << (b.store_layout == "exact" ? "allocated_capacity" : "occupied_rows")
                << "\nindexed_union_store_live_bytes_after_discovery\t" << b.store_live_bytes
                << "\nindexed_union_store_peak_bytes\t" << b.store_peak_bytes
                << "\nindexed_union_store_allocated_bytes_total\t" << b.store_allocated_bytes
                << "\nindexed_union_store_released_bytes_total\t" << b.store_released_bytes
                << "\nindexed_union_reconstruction_peak_bytes\t" << b.reconstruction_peak_bytes
                << "\nindexed_union_worker_scratch_bytes\t" << b.worker_scratch_bytes
                << "\nindexed_union_index_metadata_bytes\t" << ix.unit_bytes + ix.occurrence_bytes
                << "\nindexed_union_max_read_bp\t" << x.max_read_bp
                << "\nindexed_union_prefilter_status\t" << b.prefilter_status
                << "\nindexed_union_origin_reuse_status\t" << b.origin_reuse_status
                << "\nindexed_union_haplotypes\t" << ix.haplotypes
                << "\nindexed_union_segment_occurrences\t" << ix.segment_occurrences
                << "\nindexed_union_window_occurrences\t" << ix.window_occurrences
                << "\nindexed_union_interior_units\t" << ix.interior_units << "\nindexed_union_window_units\t" << ix.window_units
                << "\nindexed_union_fine_units\t" << ix.fine_units << "\nindexed_union_fine_occurrences\t" << ix.fine_occurrences
                << "\nindexed_union_unit_bytes\t" << ix.unit_bytes << "\nindexed_union_fine_unit_bytes\t" << ix.fine_unit_bytes
                << "\nindexed_union_coarse_interior_bytes\t" << ix.coarse_interior_bytes
                << "\nindexed_union_window_bytes\t" << ix.window_bytes << "\nindexed_union_hash_collisions\t" << ix.hash_collisions
                << "\nindexed_union_units_discovered_directly\t" << b.units_discovered_directly
                << "\nindexed_union_unit_bases_scanned\t" << b.unit_bases_scanned
                << "\nindexed_union_executed_piece_hits\t" << b.piece_hits
                << "\nindexed_union_executed_proposed_starts\t" << b.proposed_starts
                << "\nindexed_union_owned_starts_skipped_before_verification\t" << b.owned_starts_skipped_before_verification
                << "\nindexed_union_executed_unique_starts\t" << b.unique_starts
                << "\nindexed_union_executed_full_read_verifications\t" << b.full_read_verifications
                << "\nindexed_union_executed_prefilters\t" << b.executed_prefilters
                << "\nindexed_union_executed_dps\t" << b.executed_dps
                << "\nindexed_union_emission_cache_entries_held\t" << b.emission_cache_entries_held
                << "\nindexed_union_emission_cache_table_bytes\t" << b.emission_cache_table_bytes
                << "\nindexed_union_emission_cache_slab_bytes\t" << b.emission_cache_slab_bytes
                << "\nindexed_union_emission_cache_component_bytes\t" << b.emission_cache_component_bytes
                << "\nindexed_union_emission_cache_budgeted_bytes\t" << b.emission_cache_budgeted_bytes
                << "\nunion_end_phys_footprint_bytes\t" << g_union_end_footprint
                << "\nindexed_union_unit_placements_kept\t" << b.unit_placements_kept
                << "\nindexed_union_coarse_units_reconstructed\t" << b.coarse_units_reconstructed
                << "\nindexed_union_reconstructed_placements\t" << b.reconstructed_placements
                << "\nindexed_union_expanded_placements\t" << b.expanded_placements
                << "\nindexed_union_cached_placements_consumed\t" << b.placements_consumed_last
                << "\nindexed_union_executed_physical_origins\t" << b.counters_last.physical_origins
                << "\nindexed_union_direct_physical_origins\t" << b.direct_physical_origins
                << "\nindexed_union_served_physical_origins\t" << b.served_physical_origins
                << "\nindexed_union_origins_constructed\t" << b.origins_constructed
                << "\nindexed_union_origin_entities_eligible\t" << b.origin_entities_eligible
                << "\nindexed_union_origin_entities_admitted\t" << b.origin_entities_admitted
                << "\nindexed_union_origin_entities_declined\t" << b.origin_entities_declined
                << "\nindexed_union_origin_held_bytes\t" << b.origin_held_bytes
                << "\nindexed_union_columns_computed\t" << b.columns_computed_total
                << "\nindexed_union_pairs_computed\t" << b.pairs_computed_total
                << "\nthreads_indexed_union_discovery\t" << b.threads_used_discovery
                << "\nthreads_indexed_union_pairing\t" << b.threads_used_pairing
                << std::fixed << std::setprecision(3)
                << "\nseconds_indexed_union_index_build\t" << x.seconds_index_build
                << "\nseconds_indexed_union_preparation\t" << x.seconds_preparation
                << "\nseconds_indexed_union_discovery\t" << b.seconds_discovery
                << "\nseconds_indexed_union_reconstruction\t" << b.seconds_reconstruction
                << "\nseconds_indexed_union_origin_construction\t" << b.seconds_origin_construction
                << "\nseconds_indexed_union_expansion_and_pairing\t" << b.seconds_expansion_and_pairing
                << "\nseconds_indexed_union_combination\t" << b.seconds_candidate_scoring
                << "\nseconds_indexed_union_release\t" << x.seconds_release
                << "\nworker_seconds_indexed_union_verification\t" << b.worker_seconds_verification
                << "\nworker_seconds_indexed_union_pairing\t" << b.worker_seconds_pairing << '\n';
            out.flags(saved_flags);
            out.precision(saved_precision);
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
            // The combined block quality, counted and listed from the rows the tables carry.
            out << "quality_rule\t" << kQualityRule << '\n';
            if (summary->block_qualities.empty()) {
                out << "quality_pass_blocks\tNA\nquality_lowqual_blocks\tNA\nquality_pass_block_ids\tNA\n"
                       "quality_lowqual_block_ids\tNA\nquality_reason_representation_conflict\tNA\n"
                       "quality_reason_representation_not_evaluated\tNA\nquality_reason_no_read_support\tNA\n";
            } else {
                std::string pass_ids, low_ids;
                std::size_t pass = 0, low = 0, conflict = 0, unevaluated = 0, unsupported = 0, library_sensitive = 0;
                for (std::size_t i = 0; i < summary->block_qualities.size(); ++i) {
                    const BlockQuality& quality = summary->block_qualities[i];
                    std::string& ids = quality.filter == "PASS" ? pass_ids : low_ids;
                    ids += (ids.empty() ? "" : ",") + summary->quality_block_ids[i];
                    (quality.filter == "PASS" ? pass : low) += 1;
                    const std::string with_commas = "," + quality.reason + ",";
                    if (with_commas.find(",representation_conflict,") != std::string::npos) ++conflict;
                    if (with_commas.find(",representation_not_evaluated,") != std::string::npos) ++unevaluated;
                    if (with_commas.find(",no_read_support,") != std::string::npos) ++unsupported;
                    if (with_commas.find(",library_sensitive,") != std::string::npos) ++library_sensitive;
                }
                out << "quality_pass_blocks\t" << pass << "\nquality_lowqual_blocks\t" << low
                    << "\nquality_pass_block_ids\t" << (pass_ids.empty() ? "none" : pass_ids)
                    << "\nquality_lowqual_block_ids\t" << (low_ids.empty() ? "none" : low_ids)
                    << "\nquality_reason_representation_conflict\t" << conflict
                    << "\nquality_reason_representation_not_evaluated\t" << unevaluated
                    << "\nquality_reason_no_read_support\t" << unsupported << '\n';
                if (options.library_path == "auto") out << "quality_reason_library_sensitive\t" << library_sensitive << '\n';
            }
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
            if (options.block_support_rule != "legacy") {   // M3-E (EXPERIMENTAL): recorded only when the experiment is on
                std::map<std::string, std::size_t> counts;
                for (const BlockDataSupport& b : summary->data_support.blocks) ++counts[b.label];
                for (const BlockDataSupportV2& b : summary->data_support_v2.blocks) ++counts[b.label];
                const PhaseCompetitorClasses& phase_classes = options.block_support_rule == "full_data_likelihood_v2_experimental"
                                                                  ? summary->data_support_v2.phase : summary->data_support.phase;
                std::string labels;
                for (const auto& [k, v] : counts) labels += (labels.empty() ? "" : ",") + k + ":" + std::to_string(v);
                out << "block_support_rule\t" << options.block_support_rule << "\nblock_support_scope\t"
                    << "conditional_on_evaluated_candidates_not_calibrated\n" << (options.block_support_rule == "full_data_likelihood_v2_experimental" ? "block_support_labels_v2x" : "block_support_labels_v1x") << "\t"
                    << (labels.empty() ? std::string("NOT_EVALUATED") : labels) << "\nphase_competitors_genotype_changing\t"
                    << phase_classes.genotype_changing << "\nphase_competitors_genotype_preserving\t"
                    << phase_classes.genotype_preserving << "\nphase_competitors_mixed\t" << phase_classes.mixed << '\n';
            }
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
        // Served by union origin reuse, not executed (0 without it); merged_state_classes is the logical census.
        out << "reused_physical_origins\t" << counters.reused_physical_origins << '\n';
        out << "reused_state_aggregates\t" << counters.reused_state_aggregates << '\n';
        out << "reused_latent_alignment_pairs\t" << counters.reused_latent_alignment_pairs << '\n';
        out << "reused_valid_fr_states\t" << counters.reused_valid_fr_states << '\n';
        out << "merged_state_classes\t" << counters.merged_state_classes << '\n';
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
    if (options.setting_sources.count("output_schema") != 0) out << "output_schema\t" << options.output_schema << '\n';
    if (sequence_first(options)) out << sequence_report_status_keys(g_sequence_report);
    if (g_incumbent.enabled) {   // S2-A: the early evaluation's account and the checkpoints (work beside, never reset)
        const EarlyIncumbent& e = g_incumbent;
        const FragmentScoringCounters& c = e.early_counters;
        out << std::setprecision(6) << "early_incumbent\texperimental_v1\nearly_incumbent_status\t" << e.early_status
            << "\nearly_incumbent_selection\t" << e.selection << "\nearly_incumbent_pairs_examined\t" << e.pairs_examined
            << "\nearly_incumbent_path_prior_evaluations\t" << e.prior_evaluations
            << "\nearly_incumbent_seconds_preparation\t" << e.seconds_preparation
            << "\nearly_incumbent_seconds_selection\t" << e.seconds_selection << "\nearly_incumbent_seconds_evaluation\t"
            << e.seconds_evaluation << "\nearly_incumbent_candidate_fragment_evaluations\t" << e.early_fragment_evaluations
            << "\nearly_incumbent_piece_hits\t" << c.piece_hits << "\nearly_incumbent_full_read_verifications\t"
            << c.full_read_verifications << "\nearly_incumbent_dp_evaluations\t" << c.dp_evaluations
            << "\nearly_incumbent_accepted_mate_placements\t" << c.accepted_mate_placements
            << "\nearly_incumbent_physical_origins\t" << c.physical_origins << "\nearly_incumbent_statistical_state_classes\t"
            << c.statistical_state_classes << "\nearly_incumbent_footprint_after_bytes\t" << e.footprint_after_early
            << "\nearly_incumbent_union_check\t" << e.union_check;
        if (!e.early_first.empty()) {   // the early candidate and its components, exact (17 digits round-trip)
            out << std::defaultfloat << std::setprecision(17) << "\nearly_incumbent_hap1\t" << join_alleles(e.early_first)
                << "\nearly_incumbent_hap2\t" << join_alleles(e.early_second) << "\nearly_incumbent_log_likelihood\t"
                << e.early_score.log_likelihood << "\nearly_incumbent_fragment_term\t" << e.early_score.fragment_term
                << "\nearly_incumbent_exposure_term\t" << e.early_score.exposure_term
                << "\nearly_incumbent_absolute_term_sum\t" << e.early_score.absolute_term_sum
                << "\nearly_incumbent_fragments_with_mass\t" << e.early_score.fragments_with_mass
                << "\nearly_incumbent_sequence_prior\t" << e.early_prior << std::fixed << std::setprecision(6);
        }
        if (g_progressive.enabled) {   // S2-D: the schedule's stages, the prefix checkpoint, the repeated early evaluation
            const ProgressiveUnion& pu = g_progressive;
            out << "\nunion_progressive_schedule\texperimental_v1\nunion_progressive_completed\t" << (pu.completed ? "yes" : "no")
                << "\nunion_progressive_refusal\t" << pu.refusal << "\nunion_progressive_prefix_checkpoint\t"
                << (pu.prefix_committed ? "committed" : "none") << "\nunion_progressive_prefix_pairs\t" << pu.prefix_pairs
                << "\nunion_progressive_prefix_tied_groups\t" << pu.prefix_tied_groups
                << "\nunion_progressive_seconds_union_start\t" << pu.seconds_union_start
                << "\nunion_progressive_seconds_prefix_checkpoint\t" << pu.seconds_prefix_checkpoint;
            for (std::size_t i = 0; i < pu.stages.size(); ++i) {
                const ProgressiveStageReport& st = pu.stages[i];
                const std::string k = "\nunion_progressive_stage" + std::to_string(i + 1) + "_";
                out << k << "candidates\t" << st.candidates << k << "slots\t" << st.slots << k << "units_discovered\t"
                    << st.units_discovered << k << "units_reused\t" << st.units_reused << k << "columns_computed\t"
                    << st.columns_computed << k << "columns_reused\t" << st.columns_reused << k << "executed_dps\t" << st.executed_dps
                    << k << "direct_physical_origins\t" << st.direct_physical_origins << k << "served_physical_origins\t"
                    << st.served_physical_origins << k << "origins_constructed\t" << st.origins_constructed << k
                    << "store_live_bytes\t" << st.store_live_bytes << k << "origin_held_bytes\t" << st.origin_held_bytes << k
                    << "footprint_bytes\t" << st.footprint_bytes << k << "seconds_since_union_start\t" << st.seconds_since_start;
            }
            if (!pu.stages.empty()) {
                // stage 1 re-evaluates the S2-A early pair inside the union index: repeated work, reported, not hidden
                out << "\nunion_progressive_repeated_early_evaluation_columns\t" << pu.stages[0].columns_computed
                    << "\nunion_progressive_repeated_early_evaluation_executed_dps\t" << pu.stages[0].executed_dps;
            }
        }
        out << "\ncheckpoint_commits\t" << e.commits
            << "\ncheckpoint_last_stage\t" << e.stage << "\ncheckpoint_evaluated_candidates\t" << e.evaluated
            << "\ncheckpoint_seconds_first_evaluated_prediction\t" << e.seconds_first << "\ncheckpoint_failure\t"
            << e.failure << '\n';
        out << std::fixed << std::setprecision(3);
    }
    if (g_prepared_use.used) {
        // I3-A: this sample was analysed against a prepared locus (stage 1 and the marker panel were not rebuilt)
        out << "prepared_locus\t1\nprepared_locus_id\t" << g_prepared_use.id << "\nprepared_locus_key_sha256\t"
            << g_prepared_use.key_sha256 << '\n';
        // I3-C: runtime facts only (the index itself never records who loaded it)
        out << "prepared_locus_origin\t" << g_prepared_use.origin << '\n';
        if (!g_prepared_use.file_checksum.empty()) out << "prepared_locus_file_sha256\t" << g_prepared_use.file_checksum << '\n';
        out << "prepared_locus_panel_source\t" << g_prepared_use.panel_source << '\n';   // I3-D
        if (!options.panel_index_path.empty()) {   // M15-E: the public --panel-index route (runtime facts)
            out << "panel_index\t" << options.panel_index_path << "\npanel_index_file_bytes\t" << options.panel_index_file_bytes
                << "\npanel_index_load_seconds\t" << options.panel_index_load_seconds << "\npanel_index_panel_rebuild_reason\t"
                << (options.panel_index_panel_rebuild_reason.empty() ? std::string("none") : options.panel_index_panel_rebuild_reason)
                << '\n';
        }
    }
    if (options.setting_sources.count("adaptive_emission_handoff") != 0) {
        // I2-K (EXPERIMENTAL): what the union-to-adaptive emission handoff did (inherited entries are not stage work)
        const EmissionHandoffReport h = emission_handoff_report();
        out << "adaptive_emission_handoff\t" << (options.adaptive_emission_handoff ? "1" : "0")
            << "\nadaptive_emission_handoff_detached\t" << (h.detached ? "1" : "0")
            << "\nadaptive_emission_handoff_offered\t" << (h.offered ? "1" : "0")
            << "\nadaptive_emission_handoff_adopted\t" << (h.adopted ? "1" : "0")
            << "\nadaptive_emission_handoff_refusal\t" << h.refusal
            << "\nadaptive_emission_handoff_entries\t" << h.entries
            << "\nadaptive_emission_handoff_allocated_bytes\t" << h.allocated_bytes
            << "\nadaptive_emission_handoff_budgeted_bytes\t" << h.budgeted_bytes
            << "\nadaptive_emission_handoff_inherited_hits\t" << h.inherited_hits << '\n';
    }
    if (options.setting_sources.count("gap_recursion_backend") != 0) {
        out << "gap_recursion_backend\t" << gap_recursion_backend_name(options.gap_recursion_fixed_geometry)
            << "\ngap_recursion_backend_scope\tunion_adaptive_dosage\n";
    }
    if (options.setting_sources.count("combination_prepared_library_terms") != 0) {
        // I2-F (EXPERIMENTAL): what preparing the run's library terms removed (counted per combination call)
        const PreparedLibraryTermsTally t = prepared_library_terms_tally();
        out << "combination_prepared_library_terms\t" << (options.combination_prepared_library_terms ? "1" : "0")
            << "\ncombination_prepared_scope\tunion_and_adaptive_combination_not_dosage"
            << "\ncombination_prepared_calls\t" << t.combinations << "\ncombination_prepared_terms\t" << t.terms
            << "\ncombination_hypothetical_scalar_logs\t" << t.hypothetical_scalar_logs
            << "\ncombination_prepared_logs\t" << t.prepared_logs
            << "\ncombination_net_logs_avoided\t" << t.net_logs_avoided
            << "\ncombination_hypothetical_scalar_validations\t" << t.hypothetical_scalar_validations
            << "\ncombination_prepared_validations\t" << t.prepared_validations
            << "\ncombination_net_validations_avoided\t" << t.net_validations_avoided << '\n';
    }
    if (options.setting_sources.count("marker_proposal_reuse_table_bytes") != 0) {
        // I2-E (EXPERIMENTAL): the proposal emissions' contribution reuse -- work and storage only
        const ProposalReuseReport& q = g_proposal_reuse;
        out << "marker_proposal_reuse\t" << (q.enabled ? "on" : "off") << "\nmarker_proposal_reuse_table_budget_bytes\t"
            << q.table_budget_bytes << "\nmarker_proposal_contribution_visits\t" << q.visits
            << "\nmarker_proposal_contributions_executed\t" << q.executed << "\nmarker_proposal_contributions_reused\t"
            << q.reused << "\nmarker_proposal_contributions_fallback\t" << q.fallback
            << "\nmarker_proposal_reuse_markers_tabled\t" << q.markers_tabled << "\nmarker_proposal_reuse_markers_not_tabled\t"
            << q.markers_not_tabled << "\nmarker_proposal_reuse_peak_table_bytes\t" << q.peak_table_bytes
            << "\nmarker_proposal_reuse_peak_planning_bytes\t" << q.peak_planning_bytes
            << "\nmarker_proposal_reuse_budget_scope\tretained_table_and_index_only"
            << "\nseconds_marker_proposal_reuse_table_build\t" << q.seconds_table_build
            << "\nseconds_marker_proposal_emissions_total\t" << q.seconds_total << '\n';
    }
    if (options.timing_detail) {
        // I2-C (measurement only): disjoint wall sub-stages of the model and marker stages, and their unattributed rest
        double model_named = 0.0, markers_named = 0.0;
        for (const auto& [name, value] : timings.detail) {
            out << "seconds_detail_" << name << '\t' << value << '\n';
            (name.rfind("model_", 0) == 0 ? model_named : markers_named) += value;
        }
        out << "seconds_detail_model_unattributed\t" << timings.model - model_named << '\n';
        out << "seconds_detail_markers_unattributed\t" << timings.markers - markers_named << '\n';
    }
    if (repair != nullptr) g_junction_block_flags = repair->junction.block_eligible;   // M11-C: diagnostic output only
    if (repair != nullptr) write_repair_report(out, *repair);
    commit_status(path, out.str());
}

} // namespace

BlockQuality block_quality(const std::string& representation_filter, const std::string& genotype_support,
                           bool library_sensitive) {
    std::string reason;
    const auto add = [&reason](const char* item) { reason += (reason.empty() ? "" : ",") + std::string(item); };
    if (representation_filter == "LOWQUAL_REPRESENTATION_CONFLICT") {
        add("representation_conflict");
    } else if (representation_filter != "PASS") {
        add("representation_not_evaluated");   // NOT_EVALUATED, NA or anything else: fail closed
    }
    if (genotype_support != "READS") add("no_read_support");
    if (library_sensitive) add("library_sensitive");
    if (reason.empty()) return {"PASS", "none"};
    return {"LOW_QUALITY", reason};
}

std::string with_library_flag(const std::string& locus_filter, const std::string& flag) {
    if ((";" + locus_filter + ";").find(";" + flag + ";") != std::string::npos) return locus_filter;
    if (locus_filter.empty() || locus_filter == "PASS" || locus_filter == "NA") return flag;
    if (flag == "LOWQUAL_LIBRARY_SENSITIVE" && locus_filter == "NOT_EVALUATED") return flag;
    return locus_filter + ";" + flag;
}

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
        case GenotypeRunStatus::UnionDiagnosticOnly: return "UNION_DIAGNOSTIC_ONLY";
        case GenotypeRunStatus::RefinementReplayDiagnosticOnly: return "REFINEMENT_REPLAY_DIAGNOSTIC_ONLY";
        case GenotypeRunStatus::LibraryRefused: return "LIBRARY_REFUSED";
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

// One named setting, validated as a profile validates it; context names the source in every error.
void apply_genotype_setting(const std::string& context, const std::string& key, const std::string& value,
                            GenotypeCommandOptions& options) {
    const auto whole_u64 = [&](const std::string& key, const std::string& text) {
        std::size_t used = 0;
        unsigned long long value = 0;
        try {
            value = std::stoull(text, &used);
        } catch (const std::exception&) {
            throw std::runtime_error(context + ": '" + key + "' is not a non-negative integer: " + text);
        }
        if (used != text.size() || text.empty() || text[0] == '-') {
            throw std::runtime_error(context + ": '" + key + "' is not a non-negative integer: " + text);
        }
        return static_cast<std::uint64_t>(value);
    };
    const auto whole_double = [&](const std::string& key, const std::string& text) {
        std::size_t used = 0;
        double value = 0.0;
        try {
            value = std::stod(text, &used);
        } catch (const std::exception&) {
            throw std::runtime_error(context + ": '" + key + "' is not a finite non-negative number: " + text);
        }
        if (used != text.size() || !std::isfinite(value) || value < 0.0) {
            throw std::runtime_error(context + ": '" + key + "' is not a finite non-negative number: " + text);
        }
        return value;
    };
    const auto u = [&]() { return whole_u64(key, value); };
    if (key == "union_scoring_backend") {
        if (value != "reference" && value != "indexed_fine_experimental") {
            throw std::runtime_error(context + ": 'union_scoring_backend' is reference or indexed_fine_experimental: " + value);
        }
        options.union_scoring_backend = value;
    }
    else if (key == "indexed_union_store_bytes") options.indexed_union_store_bytes = u();
    else if (key == "indexed_union_store_layout") {
        if (value != "exact" && value != "legacy_i2a") {
            throw std::runtime_error(context + ": 'indexed_union_store_layout' is exact or legacy_i2a: " + value);
        }
        options.indexed_union_store_layout = value;
    }
    else if (key == "indexed_union_evidence_layout") {   // I4-B (EXPERIMENTAL)
        if (value != "reconstructed" && value != "fine_views_experimental") {
            throw std::runtime_error(context + ": 'indexed_union_evidence_layout' is reconstructed or fine_views_experimental: " + value);
        }
        options.indexed_union_evidence_layout = value;
    }
    else if (key == "marker_proposal_reuse_table_bytes") options.marker_proposal_reuse_table_bytes = u();
    else if (key == "output_schema") {   // S1-B (EXPERIMENTAL, reporting only)
        if (value != "legacy" && value != "sequence_first_v1" && value != "sequence_first_v2") {
            throw std::runtime_error(context + ": 'output_schema' is legacy, sequence_first_v1 or sequence_first_v2: " + value);
        }
        options.output_schema = value;
    }
    else if (key == "union_progressive_schedule") {   // S2-D (EXPERIMENTAL)
        if (value != "off" && value != "experimental_v1") {
            throw std::runtime_error(context + ": 'union_progressive_schedule' is off or experimental_v1: " + value);
        }
        options.union_progressive_schedule = value;
    }
    else if (key == "copying_transition") {   // M4-B (EXPERIMENTAL)
        if (value != "constant_v1" && value != "distance_ref_core_v1_experimental") {
            throw std::runtime_error(context + ": 'copying_transition' is constant_v1 or distance_ref_core_v1_experimental: " + value);
        }
        options.copying_transition = value;
    }
    else if (key == "copying_length_bp") {   // M4-B (EXPERIMENTAL)
        std::size_t used = 0;
        double v = 0.0;
        try { v = std::stod(value, &used); } catch (const std::exception&) { used = 0; }
        if (used != value.size() || !std::isfinite(v) || !(v > 0.0)) throw std::runtime_error(context + ": 'copying_length_bp' must be a finite number > 0: " + value);
        options.copying_length_bp = v;
    }
    else if (key == "block_support_rule") {   // M3-E (EXPERIMENTAL)
        if (value != "legacy" && value != "full_data_likelihood_v1_experimental" && value != "full_data_likelihood_v2_experimental") {
            throw std::runtime_error(context + ": 'block_support_rule' is legacy, full_data_likelihood_v1_experimental or "
                                               "full_data_likelihood_v2_experimental: " + value);
        }
        options.block_support_rule = value;
    }
    else if (key == "block_diagnostics") {   // M11-C (diagnostic output, default off)
        if (value != "off" && value != "diagnostic_v1") throw std::runtime_error(context + ": 'block_diagnostics' is off or diagnostic_v1: " + value);
        options.block_diagnostics = value;
    }
    else if (key == "compound_annotation") {   // M9-J (opt-in annotation)
        if (value != "off" && value != "contract_v2") throw std::runtime_error(context + ": 'compound_annotation' is off or contract_v2: " + value);
        options.compound_annotation = value;
    }
    else if (key == "compound_interval_declarations") options.compound_interval_declarations = value;   // M9-J
    else if (key == "repair_junction_proposals") {   // M9-G (EXPERIMENTAL)
        if (value != "off" && value != "experimental_v1") {
            throw std::runtime_error(context + ": 'repair_junction_proposals' is off or experimental_v1: " + value);
        }
        options.repair_junction_proposals = value;
    }
    else if (key == "repair_joint_search") {   // M3-D (EXPERIMENTAL)
        if (value != "off" && value != "experimental_v1") {
            throw std::runtime_error(context + ": 'repair_joint_search' is off or experimental_v1: " + value);
        }
        options.repair_joint_search = value;
    }
    else if (key == "repair_round1_origin_admission") {   // S2-J (EXPERIMENTAL)
        if (value != "existing" && value != "reuse_only_experimental") {
            throw std::runtime_error(context + ": 'repair_round1_origin_admission' is existing or reuse_only_experimental: " + value);
        }
        options.repair_round1_origin_admission = value;
    }
    else if (key == "diagnostic_refinement_memory") {   // S2-H (DIAGNOSTIC ONLY)
        if (value != "0" && value != "1") throw std::runtime_error(context + ": 'diagnostic_refinement_memory' is 0 or 1: " + value);
        options.diagnostic_refinement_memory = value == "1";
    }
    else if (key == "diagnostic_refinement_replay") {   // S2-H (DIAGNOSTIC ONLY)
        if (value.empty()) throw std::runtime_error(context + ": 'diagnostic_refinement_replay' needs an inputs file");
        options.diagnostic_refinement_replay = value;
    }
    else if (key == "union_fragment_blocks") {   // S2-G (EXPERIMENTAL)
        if (value != "off" && value != "fine_core_4096_experimental") {
            throw std::runtime_error(context + ": 'union_fragment_blocks' is off or fine_core_4096_experimental: " + value);
        }
        options.union_fragment_blocks = value;
    }
    else if (key == "diagnostic_union_kernel_dump") {   // S2-C (DIAGNOSTIC ONLY)
        if (value != "0" && value != "1") throw std::runtime_error(context + ": 'diagnostic_union_kernel_dump' is 0 or 1: " + value);
        options.diagnostic_union_kernel_dump = value == "1";
    }
    else if (key == "diagnostic_union_evidence_dump") {   // S2-B (DIAGNOSTIC ONLY)
        if (value != "0" && value != "1") throw std::runtime_error(context + ": 'diagnostic_union_evidence_dump' is 0 or 1: " + value);
        options.diagnostic_union_evidence_dump = value == "1";
    }
    else if (key == "early_incumbent") {   // S2-A (EXPERIMENTAL)
        if (value != "off" && value != "experimental_v1") {
            throw std::runtime_error(context + ": 'early_incumbent' is off or experimental_v1: " + value);
        }
        options.early_incumbent = value;
    }
    else if (key == "test_interrupt_point") {   // S2-A TEST HOOK
        if (value != "none" && value != "before_first_evaluation" && value != "after_checkpoint" &&
            value != "union_scoring" && value != "adaptive_refinement" && value != "after_prefix_checkpoint" &&
            value != "refinement_seed_census" && value != "post_refinement_guard" && value != "post_refinement_outputs" &&
            value != "footprint_before_first_evaluation") {
            throw std::runtime_error(context + ": 'test_interrupt_point' is none, before_first_evaluation, after_checkpoint, "
                                     "union_scoring or adaptive_refinement: " + value);
        }
        options.test_interrupt_point = value;
    }
    else if (key == "adaptive_emission_handoff") {
        if (value != "0" && value != "1") throw std::runtime_error(context + ": 'adaptive_emission_handoff' is 0 or 1: " + value);
        options.adaptive_emission_handoff = value == "1";
    }
    else if (key == "dosage_census_backend") {
        if (value == "reference") options.dosage_census_interned = false;
        else if (value == "interned_experimental") options.dosage_census_interned = true;
        else throw std::runtime_error(context + ": 'dosage_census_backend' is reference or interned_experimental: " + value);
    }
    else if (key == "gap_recursion_backend") {
        if (value == "scaled_reference") options.gap_recursion_fixed_geometry = false;
        else if (value == "fixed_geometry_experimental") options.gap_recursion_fixed_geometry = true;
        else throw std::runtime_error(context + ": 'gap_recursion_backend' is scaled_reference or fixed_geometry_experimental: " + value);
    }
    else if (key == "combination_prepared_library_terms") {
        if (value != "0" && value != "1") {
            throw std::runtime_error(context + ": 'combination_prepared_library_terms' is 0 or 1: " + value);
        }
        options.combination_prepared_library_terms = value == "1";
    }
    else if (key == "adaptive_fused_seed") {
        if (value != "0" && value != "1") throw std::runtime_error(context + ": 'adaptive_fused_seed' is 0 or 1: " + value);
        options.adaptive_fused_seed = value == "1";
    }
    else if (key == "timing_detail") {
        if (value != "0" && value != "1") throw std::runtime_error(context + ": 'timing_detail' is 0 or 1: " + value);
        options.timing_detail = value == "1";
    }
    else if (key == "adaptive_window_discovery") {
        if (value != "reference" && value != "fine_experimental") {
            throw std::runtime_error(context + ": 'adaptive_window_discovery' is reference or fine_experimental: " + value);
        }
        options.adaptive_window_discovery = value;
    }
    else if (key == "adaptive_fine_store_bytes") options.adaptive_fine_store_bytes = u();
    else if (key == "adaptive_fine_max_served_placements") options.adaptive_fine_max_served_placements = u();
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
    else if (key == "repair_max_neighbours") {
        if (value == "model") {   // M15-D: the model's structural neighbourhood maximum
            options.repair_max_neighbours_model = true;
            options.repair_max_neighbours = 0;
        } else {
            options.repair_max_neighbours_model = false;
            options.repair_max_neighbours = u();
        }
    }
    else if (key == "repair_max_physical_origins") options.repair_max_physical_origins = u();
    else if (key == "repair_max_state_classes") options.repair_max_state_classes = u();
    else if (key == "repair_max_piece_hits") options.repair_max_piece_hits = u();
    else if (key == "repair_max_latent_pairs") options.repair_max_latent_pairs = u();
    else if (key == "repair_tie_branching") {
        if (value != "0" && value != "1") throw std::runtime_error(context + ": 'repair_tie_branching' is 0 or 1: " + value);
        options.repair_tie_branching = value == "1";
    }
    else if (key == "repair_tie_max_branches") options.repair_tie_max_branches = u();
    else if (key == "repair_tie_seed_cache_bytes") options.repair_tie_seed_cache_bytes = u();
    else if (key == "union_placement_cache_bytes") options.union_placement_cache_bytes = u();
    else if (key == "max_union_cached_placements_consumed") options.max_union_cached_placements_consumed = u();
    else if (key == "union_origin_cache_bytes") options.union_origin_cache_bytes = u();
    else if (key == "max_union_cached_origins_consumed") options.max_union_cached_origins_consumed = u();
    else if (key == "max_union_cached_states_consumed") options.max_union_cached_states_consumed = u();
    else if (key == "union_prefilter_cache_bytes") options.union_prefilter_cache_bytes = u();
    else if (key == "start_sort_scratch_bytes") options.start_sort_scratch_bytes = u();
    else if (key == "emission_lookup_timing") {
        if (value != "0" && value != "1") throw std::runtime_error(context + ": 'emission_lookup_timing' is 0 or 1: " + value);
        options.emission_lookup_timing = value == "1";
    }
    else if (key == "max_union_cached_prefilters_consumed") options.max_union_cached_prefilters_consumed = u();
    else if (key == "test_library_sensitivity_fail_view") options.test_library_sensitivity_fail_view = u();
    else if (key == "library_sensitivity_audit") {
        if (value != "0" && value != "1") throw std::runtime_error(context + ": 'library_sensitivity_audit' is 0 or 1: " + value);
        options.library_sensitivity_audit = value == "1";
    }
    else if (key == "representation_guard") {
        if (value != "0" && value != "1") throw std::runtime_error(context + ": 'representation_guard' is 0 or 1: " + value);
        options.representation_guard = value == "1";
    }
    else if (key == "footprint_limit_bytes") options.footprint_limit_bytes = u();
    else if (key == "time_limit_seconds") options.time_limit_seconds = whole_double(key, value);
    else throw std::runtime_error(context + ": unknown setting '" + key + "'");
    // Ranges, the same for a profile key and its flag, checked before any work (syncmer_s <= kmer_size is checked once
    // every source has been applied: validate_genotype_settings).
    const auto in_range = [&](bool ok, const char* rule) {
        if (!ok) throw std::runtime_error(context + ": '" + key + "' " + rule + ": " + value);
    };
    if (key == "kmer_size") in_range(options.kmer_size >= 1 && options.kmer_size <= 31, "must be 1..31");
    if (key == "syncmer_s") in_range(options.syncmer_s >= 1 && options.syncmer_s <= 31, "must be 1..31");
    if (key == "min_variant_bp") in_range(options.min_variant_bp >= 1, "must be at least 1");
    if (key == "switch_probability") in_range(options.switch_probability > 0.0 && options.switch_probability <= 0.5, "must be in (0, 0.5]");
    if (key == "proposal_state_limit") in_range(options.proposal_state_limit >= 1, "must be at least 1");
    if (key == "repair_local_pairs") in_range(options.repair_local_pairs >= 1, "must be at least 1");
    if (key == "repair_tie_max_branches") in_range(options.repair_tie_max_branches >= 1, "must be at least 1");
    if (key == "threads") in_range(options.threads >= 1, "must be at least 1");
    options.setting_sources[key] = context;
}

const std::vector<std::string>& genotype_public_settings() {
    static const std::vector<std::string> keys = {
        // scientific, each also a named flag (panvar genotype --help-advanced)
        "kmer_size", "syncmer_s", "min_variant_bp", "switch_probability", "proposal_state_limit", "mismatch_band",
        "phase_audit_candidates", "max_mosaic_candidates", "repair_rounds", "repair_local_pairs", "repair_max_neighbours",
        "repair_tie_branching", "repair_tie_max_branches", "representation_guard", "threads",
        // resource caps, caches and limits: profile keys only
        "max_piece_hits", "max_physical_origins", "max_state_classes", "max_mass_cells", "emission_cache_bytes",
        "repair_origin_cache_bytes", "repair_max_physical_origins", "repair_max_state_classes", "repair_max_piece_hits",
        "repair_max_latent_pairs", "repair_tie_seed_cache_bytes", "union_placement_cache_bytes",
        "max_union_cached_placements_consumed", "union_origin_cache_bytes", "max_union_cached_origins_consumed",
        "max_union_cached_states_consumed", "union_prefilter_cache_bytes", "start_sort_scratch_bytes",
        "max_union_cached_prefilters_consumed", "emission_lookup_timing", "library_sensitivity_audit",
        "footprint_limit_bytes", "time_limit_seconds"};
    return keys;
}

std::string genotype_setting_value(const GenotypeCommandOptions& o, const std::string& key) {
    std::ostringstream out;
    out << std::setprecision(17);
    const auto b = [](bool x) { return std::string(x ? "1" : "0"); };
    if (key == "kmer_size") out << o.kmer_size;
    else if (key == "syncmer_s") out << o.syncmer_s;
    else if (key == "min_variant_bp") out << o.min_variant_bp;
    else if (key == "switch_probability") out << o.switch_probability;
    else if (key == "proposal_state_limit") out << o.proposal_state_limit;
    else if (key == "mismatch_band") out << o.mismatch_band;
    else if (key == "phase_audit_candidates") out << o.phase_audit_candidates;
    else if (key == "max_mosaic_candidates") out << o.max_mosaic_candidates;
    else if (key == "repair_rounds") out << o.repair_rounds;
    else if (key == "repair_local_pairs") out << o.repair_local_pairs;
    else if (key == "repair_max_neighbours") {
        if (o.repair_max_neighbours_model) return std::string("model");
        out << o.repair_max_neighbours;
    }
    else if (key == "repair_tie_branching") return b(o.repair_tie_branching);
    else if (key == "repair_tie_max_branches") out << o.repair_tie_max_branches;
    else if (key == "representation_guard") return b(o.representation_guard);
    else if (key == "threads") out << o.threads;
    else if (key == "max_piece_hits") out << o.max_piece_hits;
    else if (key == "max_physical_origins") out << o.max_physical_origins;
    else if (key == "max_state_classes") out << o.max_state_classes;
    else if (key == "max_mass_cells") out << o.max_mass_cells;
    else if (key == "emission_cache_bytes") out << o.emission_cache_bytes;
    else if (key == "repair_origin_cache_bytes") out << o.repair_origin_cache_bytes;
    else if (key == "repair_max_physical_origins") out << o.repair_max_physical_origins;
    else if (key == "repair_max_state_classes") out << o.repair_max_state_classes;
    else if (key == "repair_max_piece_hits") out << o.repair_max_piece_hits;
    else if (key == "repair_max_latent_pairs") out << o.repair_max_latent_pairs;
    else if (key == "repair_tie_seed_cache_bytes") out << o.repair_tie_seed_cache_bytes;
    else if (key == "union_placement_cache_bytes") out << o.union_placement_cache_bytes;
    else if (key == "max_union_cached_placements_consumed") out << o.max_union_cached_placements_consumed;
    else if (key == "union_origin_cache_bytes") out << o.union_origin_cache_bytes;
    else if (key == "max_union_cached_origins_consumed") out << o.max_union_cached_origins_consumed;
    else if (key == "max_union_cached_states_consumed") out << o.max_union_cached_states_consumed;
    else if (key == "union_prefilter_cache_bytes") out << o.union_prefilter_cache_bytes;
    else if (key == "start_sort_scratch_bytes") out << o.start_sort_scratch_bytes;
    else if (key == "max_union_cached_prefilters_consumed") out << o.max_union_cached_prefilters_consumed;
    else if (key == "emission_lookup_timing") return b(o.emission_lookup_timing);
    else if (key == "library_sensitivity_audit") return b(o.library_sensitivity_audit);
    else if (key == "footprint_limit_bytes") out << o.footprint_limit_bytes;
    else if (key == "time_limit_seconds") out << o.time_limit_seconds;
    else throw std::invalid_argument("not a public setting: " + key);
    return out.str();
}

void validate_genotype_settings(const GenotypeCommandOptions& o) {
    // I2-A: the experimental union backend's own requirements, before any work
    if (o.union_scoring_backend == "indexed_fine_experimental") {
        if (o.indexed_union_store_bytes == 0) {
            throw std::runtime_error("union_scoring_backend indexed_fine_experimental needs a positive indexed_union_store_bytes");
        }
        if (o.max_union_cached_placements_consumed == 0) {
            throw std::runtime_error("union_scoring_backend indexed_fine_experimental needs a positive "
                                     "max_union_cached_placements_consumed (every expanded placement is served from its store)");
        }
        if (o.diagnostic_union_work_trace || o.diagnostic_union_seeding_profile > 0 ||
            !o.diagnostic_union_verification_census.empty()) {
            throw std::runtime_error("the union diagnostics need union_scoring_backend reference");
        }
    } else if (o.indexed_union_store_bytes != 0) {
        throw std::runtime_error("indexed_union_store_bytes applies only to union_scoring_backend indexed_fine_experimental");
    }
    if (o.indexed_union_store_layout != "exact" && o.union_scoring_backend != "indexed_fine_experimental") {
        throw std::runtime_error("indexed_union_store_layout applies only to union_scoring_backend indexed_fine_experimental");
    }
    if (o.indexed_union_evidence_layout != "reconstructed") {   // I4-B
        if (o.union_scoring_backend != "indexed_fine_experimental") {
            throw std::runtime_error("indexed_union_evidence_layout applies only to union_scoring_backend indexed_fine_experimental");
        }
        if (o.indexed_union_store_layout != "exact") {
            throw std::runtime_error("indexed_union_evidence_layout fine_views_experimental needs indexed_union_store_layout exact");
        }
    }
    // I2-B: the experimental adaptive window discovery's own requirements
    if (o.adaptive_window_discovery == "fine_experimental") {
        if (o.adaptive_fine_store_bytes == 0 || o.adaptive_fine_max_served_placements == 0) {
            throw std::runtime_error("adaptive_window_discovery fine_experimental needs a positive adaptive_fine_store_bytes "
                                     "and a positive adaptive_fine_max_served_placements");
        }
    } else if (o.adaptive_fine_store_bytes != 0 || o.adaptive_fine_max_served_placements != 0) {
        throw std::runtime_error("adaptive_fine_store_bytes and adaptive_fine_max_served_placements apply only to "
                                 "adaptive_window_discovery fine_experimental");
    }
    if (o.syncmer_s > o.kmer_size) {
        throw std::runtime_error("syncmer_s (" + std::to_string(o.syncmer_s) + ", from " +
                                 (o.setting_sources.count("syncmer_s") ? o.setting_sources.at("syncmer_s") : std::string("default")) +
                                 ") must not exceed kmer_size (" + std::to_string(o.kmer_size) + ", from " +
                                 (o.setting_sources.count("kmer_size") ? o.setting_sources.at("kmer_size") : std::string("default")) + ")");
    }
}

std::vector<std::pair<std::string, std::string>> genotype_preset_settings(const std::string& name) {
    if (name == "long-read") throw std::runtime_error("--preset long-read: long-read scoring is not implemented");
    if (name != "short-read") throw std::runtime_error("unknown --preset '" + name + "' (available: short-read)");
    // The production algorithm (adaptive repair, exact tie branching, the representation guard), the short-read
    // technology constants and (M15-D) the retained release-candidate configuration, so an ordinary run needs no profile
    // (docs/guides/genotype_short_read.md, defaults table). Grouped by role:
    //   exact computational backends and routes (results gated identical on completed runs) with their cache budgets and
    //   consumption capacities; one candidate-generation setting (repair_junction_proposals, JP1.1); the early evaluated
    //   incumbent (a PROVISIONAL checkpoint on a budget stop; completed runs unchanged); the sequence-first v2 report with
    //   the compound annotation; the total-footprint guard (6 GiB) and the state-class work cap used by the release
    //   candidate; and a model-derived repair neighbourhood bound. No time limit, no dosage stage.
    // Preset-sourced features whose requirements a flag or profile removes are switched off and recorded as inactive
    // at argument resolution; the same features set explicitly keep their refusals.
    return {{"repair_rounds", "32"}, {"repair_tie_branching", "1"}, {"representation_guard", "1"},
            {"kmer_size", "31"}, {"syncmer_s", "15"}, {"min_variant_bp", "50"}, {"mismatch_band", "6"},
            {"switch_probability", "0.01"}, {"proposal_state_limit", "256"}, {"phase_audit_candidates", "0"},
            // exact computational backends and routes
            {"union_scoring_backend", "indexed_fine_experimental"}, {"adaptive_window_discovery", "fine_experimental"},
            {"indexed_union_store_layout", "exact"}, {"indexed_union_evidence_layout", "fine_views_experimental"},
            {"adaptive_fused_seed", "1"}, {"marker_proposal_reuse_table_bytes", "4194304"},
            {"combination_prepared_library_terms", "1"}, {"gap_recursion_backend", "fixed_geometry_experimental"},
            {"dosage_census_backend", "interned_experimental"}, {"union_fragment_blocks", "fine_core_4096_experimental"},
            {"repair_round1_origin_admission", "reuse_only_experimental"},
            // cache and store budgets (retained memory for reuse; they do not cap work)
            {"indexed_union_store_bytes", "2147483648"}, {"adaptive_fine_store_bytes", "536870912"},
            {"emission_cache_bytes", "536870912"}, {"repair_origin_cache_bytes", "805306368"},
            {"union_placement_cache_bytes", "536870912"}, {"union_origin_cache_bytes", "1073741824"},
            {"union_prefilter_cache_bytes", "134217728"}, {"start_sort_scratch_bytes", "335544320"},
            // consumption capacities of the reuse routes, and the state-class work cap of the release candidate
            {"adaptive_fine_max_served_placements", "2000000000"}, {"max_union_cached_placements_consumed", "2000000000"},
            {"max_union_cached_origins_consumed", "2000000000"}, {"max_union_cached_states_consumed", "2000000000"},
            {"max_union_cached_prefilters_consumed", "2000000000"}, {"max_state_classes", "600000000"},
            // candidate generation (affects predictions): junction proposals
            {"repair_junction_proposals", "experimental_v1"},
            // budget-stop behaviour: the early evaluated incumbent
            {"early_incumbent", "experimental_v1"},
            // reporting
            {"output_schema", "sequence_first_v2"}, {"compound_annotation", "contract_v2"}, {"timing_detail", "1"},
            // resources: total-footprint guard; model-derived repair neighbourhood bound
            {"footprint_limit_bytes", "6442450944"}, {"repair_max_neighbours", "model"}};
}

LibraryModel genotype_preset_library(const std::string& name) {
    genotype_preset_settings(name);   // refuses long-read and unknown names
    LibraryModel m;
    m.insert_mean_bp = std::numeric_limits<double>::quiet_NaN();
    m.insert_sd_bp = std::numeric_limits<double>::quiet_NaN();
    m.support_sigmas = 4.0;
    m.substitution_error = 0.002;
    m.outlier_weight = 0.001;
    // A declared technology convention, not an estimate: 30x diploid depth with 2 x 150 bp pairs, lambda = depth / (4 L).
    // Under --library auto the marker-fit intensity replaces it (docs/algorithms/genotype.md, stage 3).
    m.fragment_start_intensity = 0.05;
    m.provenance = "preset short-read";
    return m;
}

void apply_genotype_profile(const std::string& path, GenotypeCommandOptions& options) {
    std::ifstream file(path);
    if (!file) throw std::runtime_error("cannot open genotype profile: " + path);
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
        // M15-J: no profile reaches an older reporting contract or the dosage stage. output_schema may only restate the
        // fixed contract (the evidence profiles record it, and their bytes are cited by the sealed packages).
        if (key == "output_schema" && value != "sequence_first_v2") {
            throw std::runtime_error("profile " + path + ": 'output_schema' cannot select a reporting contract: genotyping "
                                     "reports sequence_first_v2 (a profile may only restate that value)");
        }
        if (key == "dosage_mode" || key == "dosage_audit" || key == "dosage_check") {
            throw std::runtime_error("profile " + path + ": '" + key + "' is not a profile setting: dosage is an audit only, "
                                     "requested with --dosage-check, and never replaces the prediction");
        }
        if (key == "profile") options.profile_name = value;
        else apply_genotype_setting("profile " + path, key, value, options);
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
    for (const std::string& suffix : sequence_report_outputs()) std::remove((out_prefix + suffix).c_str());   // S1-C
    for (const std::string& suffix : dosage::dosage_stage_outputs()) std::remove((out_prefix + suffix).c_str());
    for (const char* suffix : kPartialOutputs) std::remove((out_prefix + suffix).c_str());
    // S2-H: the refinement diagnostics are run-owned too (removed at the start only: they must survive a stop)
    std::remove((out_prefix + ".block_diagnostics.tsv").c_str());   // M11-C: run-owned diagnostic, removed at the start
    for (const char* suffix : {".refinement_memory_live.tsv", ".refinement_seed_report.tsv", ".refinement_replay_inputs.tsv",
                               ".refinement_replay.tsv",
                               ".refinement_replay_rounds.tsv", ".refinement_replay_scored.tsv"}) {
        std::remove((out_prefix + suffix).c_str());
    }
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

GenotypeRunStatus run_genotype(const GenotypeCommandOptions& options) { return run_genotype(options, nullptr); }


// ---- I3-A: the prepared locus -------------------------------------------------------------------------------------------
// Stage 1, shared by the ordinary run (built inline for one sample) and prepare_genotype_locus (built once). Fills `out`
// as it goes, so a caller that catches a failure can report how far it got.
struct LocusModelBuild {
    LocusInputResult built;
    std::optional<PreparedGenotypeModel> model;
    std::string fingerprint = "NA";
    std::string resolved;
    std::unordered_map<std::string, std::string> node_sequences;
};

void build_locus_model(const GenotypeCommandOptions& options, const std::string& staged, LocusModelBuild& out,
                       std::vector<std::pair<std::string, double>>& detail) {
    const Clock::time_point mark = Clock::now();
    GfaModel gfa = read_gfa_model(options.gfa_path);
    GraphSortOptions sort_options;
    sort_options.reference_path = options.reference_path;
    sort_options.flip = true;
    const std::string resolved = sort_graph_reference(gfa, sort_options).resolved_reference;
    if (g_watch != nullptr) {
        // Under the watch lock: a watchdog stop either sees the finished file and removes it, or
        // stops the process before it is begun. It can never leave one half-written behind.
        std::lock_guard<std::mutex> lock(g_watch->mutex);
        write_gfa_model(staged, gfa);
    } else {
        write_gfa_model(staged, gfa);
    }
    ParseGfaOptions parse_options;
    const Graph graph = parse_gfa(staged, parse_options);
    validate_graph_paths(graph, "genotype", true, true);
    detail.emplace_back("model_graph_load_sort_parse", seconds_since(mark));
    const Clock::time_point bubbles_mark = Clock::now();

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

    detail.emplace_back("model_bubbles", seconds_since(bubbles_mark));
    const Clock::time_point build_mark = Clock::now();
    LocusInputOptions locus_options;
    // The locus identity is the graph's own reference frame, never the output path: a fingerprint
    // that moved when the user renamed their output would be useless for saying whether two runs
    // scored the same model.
    locus_options.locus_id = resolved;
    locus_options.reference_path_name = resolved;
    if (!options.curation_manifest_path.empty()) {
        locus_options.excluded_templates = read_curation_manifest(options.curation_manifest_path);
    }
    out.built = build_locus_input(graph, bubbles, locus_options);
    out.model.emplace(prepare_genotype_model(out.built.input));
    // The authoritative check, on real input, before any statistic is computed from these blocks.
    verify_complete_panel_round_trip(*out.model, out.built.authoritative);
    out.fingerprint = model_fingerprint(*out.model);
    // A model a node split would change is not scored. The round trip cannot see this case.
    require_representation_independent(out.built.report);
    detail.emplace_back("model_locus_build_prepare_verify", seconds_since(build_mark));
    out.resolved = resolved;
    out.node_sequences = std::move(gfa.seq);
}

PreparedLocus::PreparedLocus(Parts&& parts) : parts_(std::make_unique<Parts>(std::move(parts))) {}
PreparedLocus::~PreparedLocus() = default;
const PreparedLocusKey& PreparedLocus::key() const { return parts_->key; }
const PreparedGenotypeModel& PreparedLocus::model() const { return *parts_->model; }
const std::string& PreparedLocus::fingerprint() const { return parts_->fingerprint; }
const LocusInputReport& PreparedLocus::report() const { return parts_->report; }
const std::string& PreparedLocus::resolved_reference() const { return parts_->resolved; }
const std::unordered_map<std::string, std::string>& PreparedLocus::node_sequences() const { return parts_->node_sequences; }
const PreparedMarkerPanel& PreparedLocus::panel() const { return *parts_->panel; }
const std::vector<std::pair<std::string, double>>& PreparedLocus::preparation_seconds() const { return parts_->seconds; }
std::uint64_t PreparedLocus::id() const { return parts_->id; }
const std::string& PreparedLocus::origin() const { return parts_->origin; }
const std::string& PreparedLocus::file_checksum() const { return parts_->file_checksum; }
const std::string& PreparedLocus::stored_provenance() const { return parts_->stored_provenance; }
const std::string& PreparedLocus::panel_source() const { return parts_->panel_source; }

std::string PreparedLocus::digest() const {
    std::ostringstream o;
    o << "model\t" << model_fingerprint(*parts_->model) << "\nfingerprint\t" << parts_->fingerprint << "\nresolved\t"
      << parts_->resolved << "\nkey\t" << parts_->key.describe() << '\n';
    const LocusInputReport& r = parts_->report;
    o << "report\t" << r.variant_blocks << ',' << r.context_blocks << ',' << r.total_source_alleles << ','
      << r.total_canonical_alleles << ',' << r.forward_templates << ',' << r.reverse_templates << ','
      << r.refused_templates.size() << ',' << r.blocks.size() << '\n';
    std::vector<const std::pair<const std::string, std::string>*> nodes;
    for (const auto& item : parts_->node_sequences) nodes.push_back(&item);
    std::sort(nodes.begin(), nodes.end(), [](const auto* a, const auto* b) { return a->first < b->first; });
    for (const auto* item : nodes) o << item->first << '\t' << item->second << '\n';
    const PreparedMarkerPanel& m = *parts_->panel;
    o << "panel\t" << m.kmer_size << ',' << m.syncmer_s << ',' << m.total_alleles << ',' << m.scored_alleles << ','
      << m.dropped_multiblock_markers << '\n';
    for (const std::uint64_t c : m.codes) o << c << ',';
    o << '\n';
    for (const std::uint32_t a : m.anchor_slots) o << a << ',';
    o << '\n';
    for (const std::uint32_t a : m.anchor_haploid_multiplicity) o << a << ',';
    o << '\n';
    for (const BlockMarkerInventory& b : m.blocks) {
        o << b.block_id << '\t' << b.allele_count << '\t';
        for (const std::uint32_t x : b.marker_slots) o << x << ',';
        o << '\t';
        for (const std::uint32_t x : b.fixed_haploid_multiplicity) o << x << ',';
        o << '\t';
        for (const auto& allele : b.allele_markers) {
            for (const auto& [slot, copies] : allele) o << slot << ':' << copies << ',';
            o << '|';
        }
        o << '\n';
    }
    return sha256_bytes(o.str());
}

std::string PreparedLocusKey::describe() const {
    std::ostringstream o;
    o << "gfa_sha256=" << gfa_sha256 << ";reference_path=" << reference_path << ";bubbles=" << bubbles
      << ";min_variant_bp=" << min_variant_bp << ";curation=" << curation << ";kmer_size=" << kmer_size
      << ";syncmer_s=" << syncmer_s;
    return o.str();
}

std::string PreparedLocusKey::mismatch(const PreparedLocusKey& other) const {
    if (gfa_sha256 != other.gfa_sha256) return "graph content (GFA SHA-256)";
    if (reference_path != other.reference_path) return "reference path";
    if (bubbles != other.bubbles) return "bubbles (supplied CSV content, or computed)";
    if (min_variant_bp != other.min_variant_bp) return "min_variant_bp";
    if (curation != other.curation) return "panel curation (manifest content)";
    if (kmer_size != other.kmer_size) return "marker kmer_size";
    if (syncmer_s != other.syncmer_s) return "marker syncmer_s";
    return "";
}

namespace {
std::atomic<std::uint64_t>& prepared_locus_counter() {
    static std::atomic<std::uint64_t> count{0};
    return count;
}
std::atomic<std::uint64_t>& prepared_locus_id_counter() {   // I3-C: ids are shared by prepared and loaded objects
    static std::atomic<std::uint64_t> count{0};
    return count;
}

// The key from a hasher of input files ("unreadable" never equals a real digest, so it refuses).
template <typename Hash>
PreparedLocusKey key_with(const GenotypeCommandOptions& options, Hash&& hash) {
    PreparedLocusKey key;
    key.gfa_sha256 = hash(options.gfa_path);
    key.reference_path = options.reference_path;
    key.bubbles = options.bubbles_path.empty() ? std::string("computed") : "sha256:" + hash(options.bubbles_path);
    key.min_variant_bp = options.min_variant_bp;
    key.curation = options.curation_manifest_path.empty() ? std::string("none") : "sha256:" + hash(options.curation_manifest_path);
    key.kmer_size = options.kmer_size;
    key.syncmer_s = options.syncmer_s;
    return key;
}

// The run's own provenance hashes (computed once per run for the status) give the same key without hashing again.
PreparedLocusKey key_from_run_hashes(const GenotypeCommandOptions& options) {
    return key_with(options, [&](const std::string& file) -> std::string {
        const std::string name = file == options.gfa_path ? "gfa" : file == options.bubbles_path ? "bubbles" : "panel_curation";
        for (const auto& item : g_input_hashes) {
            if (item.first == name) return item.second;
        }
        return std::string("unreadable");
    });
}

}  // namespace

PreparedLocusKey prepared_locus_key(const GenotypeCommandOptions& options) {
    return key_with(options, [](const std::string& file) {
        try {
            return sha256_file(file);
        } catch (const std::exception&) {
            return std::string("unreadable");
        }
    });
}

std::uint64_t prepared_locus_builds() { return prepared_locus_counter().load(); }
std::uint64_t next_prepared_locus_id() { return ++prepared_locus_id_counter(); }

BetweenSampleCleanup release_between_samples() {
    BetweenSampleCleanup out;
    const Clock::time_point mark = Clock::now();
    out.released = release_sample_retained_storage();
    out.free_pages_released = release_free_pages();
    out.seconds = seconds_since(mark);
    return out;
}

const char* batch_resource_verdict_name(BatchResourceVerdict verdict) {
    switch (verdict) {
        case BatchResourceVerdict::compliant: return "RESOURCE_COMPLIANT";
        case BatchResourceVerdict::ceiling_exceeded: return "RESOURCE_CEILING_EXCEEDED";
        case BatchResourceVerdict::unverified: return "RESOURCE_UNVERIFIED";
    }
    return "RESOURCE_UNVERIFIED";
}

BatchResourceLedger::BatchResourceLedger(std::uint64_t ceiling_bytes, std::string metric)
    : ceiling_(ceiling_bytes), metric_(std::move(metric)) {}

bool BatchResourceLedger::record(const std::string& phase, std::optional<std::uint64_t> current_bytes,
                                 std::optional<std::uint64_t> lifetime_peak_bytes) {
    Phase p;
    p.name = phase;
    p.current = current_bytes;
    p.peak = lifetime_peak_bytes;
    p.verified = current_bytes.has_value() && lifetime_peak_bytes.has_value();
    if (lifetime_peak_bytes) {
        const std::uint64_t before = max_peak_.value_or(0);
        p.peak_increase = *lifetime_peak_bytes > before ? *lifetime_peak_bytes - before : 0;
        max_peak_ = std::max(before, *lifetime_peak_bytes);
    }
    if (!p.verified) ++unverified_;
    // An observed breach, by either reading actually taken, is permanent; a missing reading is never zero.
    const bool breach = (lifetime_peak_bytes && *lifetime_peak_bytes > ceiling_) || (current_bytes && *current_bytes > ceiling_);
    if (!exceeded_ && breach) {
        exceeded_ = true;
        exceeded_during_ = phase;
    }
    phases_.push_back(p);
    return !exceeded_;
}

BatchResourceVerdict BatchResourceLedger::verdict() const {
    if (exceeded_) return BatchResourceVerdict::ceiling_exceeded;
    if (phases_.empty() || unverified_ > 0) return BatchResourceVerdict::unverified;
    return BatchResourceVerdict::compliant;
}

namespace {
// I3-C: the one preparation behind prepare_genotype_locus and prepare_genotype_locus_record. `record`, when given,
// receives the complete persistent content: the structural input MOVES into it (the locus never owns it), the report and
// node sequences are copied, and both complete content digests are computed from the prepared model and panel.
std::shared_ptr<const PreparedLocus> prepare_locus(const GenotypeCommandOptions& options, PreparedLocusRecord* record,
                                                   bool persist_panel) {
    PreparedLocus::Parts parts;
    parts.key = prepared_locus_key(options);
    const std::string staged = options.out_prefix + ".prepared_locus.sorted.gfa";
    struct Remove {
        std::string path;
        ~Remove() {
            std::error_code ignored;
            std::filesystem::remove(path, ignored);
        }
    } remove_staged{staged};
    LocusModelBuild built;
    build_locus_model(options, staged, built, parts.seconds);
    parts.model = std::move(built.model);
    parts.fingerprint = built.fingerprint;
    parts.resolved = built.resolved;
    if (record != nullptr) {
        record->key = parts.key;
        record->resolved_reference = built.resolved;
        record->preparation_algorithm_version = kPreparationAlgorithmVersion;
        // build_locus_model throws unless the round trip and the representation-independence check both passed
        record->round_trip_passed = true;
        record->representation_independent = true;
        record->provenance = current_prepared_locus_provenance();
        record->input = std::move(built.built.input);
        record->report = built.built.report;
        record->nodes.assign(built.node_sequences.begin(), built.node_sequences.end());
        std::sort(record->nodes.begin(), record->nodes.end(),
                  [](const auto& a, const auto& b) { return a.first < b.first; });
        record->model_fingerprint = built.fingerprint;
    }
    parts.report = std::move(built.built.report);   // the locus input and authoritative sequences are not kept
    parts.node_sequences = std::move(built.node_sequences);
    const Clock::time_point mark = Clock::now();
    parts.panel.emplace(build_marker_panel(*parts.model, options.kmer_size, options.syncmer_s));
    parts.seconds.emplace_back("markers_panel_build", seconds_since(mark));
    if (record != nullptr) {
        record->model_content_sha256 = prepared_locus_model_content_sha256(*parts.model);
        record->panel_content_sha256 = prepared_locus_panel_content_sha256(*parts.panel);
        if (persist_panel) record->panel = *parts.panel;   // I3-D: a transient copy for the writer only
    }
    parts.id = next_prepared_locus_id();
    ++prepared_locus_counter();
    return std::make_shared<const PreparedLocus>(std::move(parts));
}
}  // namespace

std::shared_ptr<const PreparedLocus> prepare_genotype_locus(const GenotypeCommandOptions& options) {
    return prepare_locus(options, nullptr, false);
}

PreparedLocusPreparation prepare_genotype_locus_record(const GenotypeCommandOptions& options, bool persist_panel) {
    PreparedLocusPreparation out;
    out.locus = prepare_locus(options, &out.record, persist_panel);
    return out;
}

GenotypeRunStatus run_genotype_body(const GenotypeCommandOptions& options, GenotypeDosageCapture* capture,
                                    const PreparedLocus* locus);

GenotypeRunStatus run_genotype(const GenotypeCommandOptions& options, GenotypeDosageCapture* capture) {
    return run_genotype_body(options, capture, nullptr);
}

GenotypeRunStatus run_genotype_prepared(const GenotypeCommandOptions& options, const PreparedLocus& locus,
                                        GenotypeDosageCapture* capture) {
    return run_genotype_body(options, capture, &locus);
}

// One body for both: `locus`, when given, replaces stage 1 and the marker-panel build (I3-A).
GenotypeRunStatus run_genotype_body(const GenotypeCommandOptions& options, GenotypeDosageCapture* capture,
                                    const PreparedLocus* locus) {
    if (locus != nullptr) {
        // I3-A: compatibility first, from the inputs' content and the preparation settings -- before any output is
        // removed or written, so a refused sample leaves its prefix exactly as it was.
        const std::string differs = locus->key().mismatch(prepared_locus_key(options));
        if (!differs.empty()) {
            throw PreparedLocusMismatch("the prepared locus does not match this sample's " + differs +
                                        "; prepare the locus again for these inputs and settings");
        }
    }
    if (const std::optional<std::string> refused = early_incumbent_refusal(options)) {   // S2-A: before the prefix
        throw std::invalid_argument(*refused);
    }
    if (const std::optional<std::string> refused = union_fragment_blocks_refusal(options)) {   // S2-G: before the prefix
        throw std::invalid_argument(*refused);
    }
    if (const std::optional<std::string> refused = refinement_diagnostic_refusal(options)) {   // S2-H: before the prefix
        throw std::invalid_argument(*refused);
    }
    if (const std::optional<std::string> refused = joint_search_refusal(options)) {   // M3-D: before the prefix
        throw std::invalid_argument(*refused);
    }
    if (const std::optional<std::string> refused = junction_proposals_refusal(options)) {   // M9-G: before the prefix
        throw std::invalid_argument(*refused);
    }
    if (const std::optional<std::string> refused = compound_annotation_refusal(options)) {   // M9-J: before the prefix
        throw std::invalid_argument(*refused);
    }
    g_compound_declarations.reset();   // M9-J: never a previous run's declarations; a malformed file refuses before the prefix
    if (options.compound_annotation == "contract_v2" && !options.compound_interval_declarations.empty()) {
        try {
            g_compound_declarations = read_compound_declarations(options.compound_interval_declarations);
        } catch (const std::runtime_error& error) {
            throw std::invalid_argument(error.what());
        }
    }
    if (const std::optional<std::string> refused = copying_transition_refusal(options)) {   // M4-B: before the prefix
        throw std::invalid_argument(*refused);
    }
    if (locus != nullptr && options.copying_transition != "constant_v1") {   // M4-B: before the prefix
        throw std::invalid_argument("copying_transition " + options.copying_transition + " needs the adapter's locus geometry; a prepared "
                                    "locus (built in-process or loaded from an index file) carries none -- run from the graph");
    }
    g_transitions = CopyingTransitions{};   // M4-B: never a previous run's vector
    g_geometry = LocusGeometry{};
    {
        std::lock_guard<std::mutex> lock(g_checkpoint_mutex);   // S2-G: this run's checkpoints only
        g_checkpoints_on = options.union_fragment_blocks != "off";
        g_checkpoints_started = Clock::now();
        g_footprint_checkpoints.clear();
    }
    g_incumbent = EarlyIncumbent{};   // S2-A
    g_progressive = ProgressiveUnion{};   // S2-D
    g_progressive.enabled = options.union_progressive_schedule == "experimental_v1";
    g_incumbent.enabled = options.early_incumbent == "experimental_v1";
    set_emission_lookup_timing(options.emission_lookup_timing);
    g_proposal_reuse = ProposalReuseReport{};   // I2-E: never a previous run's figures
    g_sequence_report = SequenceReport{};          // S1-B: likewise
    if (options.output_schema == "sequence_first_v2") g_sequence_report.schema = "sequence_first_v2";   // S1-D
    g_sequence_legacy = SequenceLegacy{};
    g_sequence_published = false;
    g_junction_block_flags.clear();   // M11-C: never a previous run's flags
    reset_prepared_library_terms_tally();         // I2-F: likewise
    g_union_end_footprint = 0;
    reset_emission_handoff_report();   // I2-K: nothing from an earlier run is reported or adopted
    g_prepared_use = PreparedLocusUse{};   // I3-A
    const Clock::time_point started = Clock::now();
    StageTimings timings;
    const std::string status_path = options.out_prefix + ".mosaic_status.tsv";
    // M15-D.1: a declared library (--library FILE) is read and validated in full BEFORE the prefix is touched, so a
    // missing, unreadable or invalid file refuses the run with every file at the prefix left as it was. Errors found
    // later (graph, reads, model, resources, writing) are not covered by this guarantee: by then the run owns the prefix.
    if (!options.library_path.empty() && options.library_path != "auto") {
        std::optional<double> preflight_background;
        (void)read_library(options.library_path, preflight_background);
    }
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
    // I2-A: the experimental union backend's requirements hold for library callers too, before any work
    if (options.union_scoring_backend != "reference" || options.indexed_union_store_bytes != 0 ||
        options.indexed_union_store_layout != "exact" || options.indexed_union_evidence_layout != "reconstructed" ||
        options.adaptive_window_discovery != "reference" || options.adaptive_fine_store_bytes != 0 ||
        options.adaptive_fine_max_served_placements != 0) {
        validate_genotype_settings(options);
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
        watch->early_incumbent = g_incumbent.enabled;   // S2-A
        watch->sequence_v2 = options.output_schema == "sequence_first_v2";   // M15-D.1
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
                        if ((g_watch->early_incumbent || g_watch->sequence_v2) && !g_watch->checkpoint_failed) {
                            // S2-A: a failure is never a fallback -- no committed checkpoint is left as a prediction
                            for (const std::string& suffix : sequence_report_outputs()) {
                                std::remove((g_watch->out_prefix + suffix).c_str());
                            }
                            atomic_write(g_watch->out_prefix + ".sequence_locus.tsv",
                                         "key\tvalue\nschema_version\tsequence_first_v2\nprediction\tNO_PREDICTION\n"
                                         "sequence_status\tNA\ncompletion\trun_failed\ncompletion_note\tan exception ended the "
                                         "run before a final status\nlegacy_run_status\tFAILED\ncheckpoint_stage\tnone\n"
                                         "exit_code\t2\n");
                            g_watch->checkpoint_committed = false;
                        }
                        atomic_write(g_watch->out_prefix + ".mosaic_status.tsv",
                                     minimal_status(*g_watch, "FAILED",
                                                    "an exception ended the run before a final status", 2));
                        if (g_watch->sequence_v2) assessment::write_assessment(g_watch->out_prefix);   // M15-G
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
            "genotype requires --gfa, --reads1, --reads2, --reference-path and --library (a file, or auto)");
    }

    // --library auto: the short-read preset's constants now, the insert estimated from the reads once the model and
    // the reads exist (stage 2), before anything is scored. --library FILE: exactly as declared.
    const bool auto_library = options.library_path == "auto";
    std::optional<double> declared_marker_background;
    LibraryModel library = auto_library ? genotype_preset_library("short-read") : read_library(options.library_path, declared_marker_background);
    if (auto_library) library.provenance = "auto: pending the insert estimate";

    // ---- stage 1: structural model -------------------------------------------------------------
    // Fail closed. Anything that stops the structural model being built, verified, or accepted as
    // representation-independent ends the run as MODEL_INCOMPLETE with the reason, rather than escaping
    // as an exception that leaves no status behind.
    Clock::time_point mark = Clock::now();
    // I3-A: stage 1 is build_locus_model, shared with prepare_genotype_locus; a compatible prepared locus replaces it.
    LocusModelBuild local;
    std::string fingerprint = "NA";
    if (locus != nullptr) {
        // Checked again against the run's own provenance hashes (an input changed since the first check refuses here).
        const std::string differs = locus->key().mismatch(key_from_run_hashes(options));
        if (!differs.empty()) {
            throw PreparedLocusMismatch("the prepared locus does not match this sample's " + differs +
                                        " (changed during the run); prepare the locus again");
        }
        fingerprint = locus->fingerprint();
        g_prepared_use.used = true;
        g_prepared_use.id = locus->id();
        g_prepared_use.key_sha256 = sha256_bytes(locus->key().describe());
        g_prepared_use.origin = locus->origin();
        g_prepared_use.file_checksum = locus->file_checksum();
        g_prepared_use.panel_source = locus->panel_source();
        timings.detail.emplace_back("model_prepared_locus_reused", seconds_since(mark));
        if (capture != nullptr) {
            capture->model = locus->model();
            capture->node_sequences = locus->node_sequences();
        }
    } else {
        try {
            build_locus_model(options, options.out_prefix + ".sorted.gfa", local, timings.detail);
            fingerprint = local.fingerprint;
            if (capture != nullptr) {
                capture->model = *local.model;
                capture->node_sequences = local.node_sequences;
            }
        } catch (const std::exception& error) {
            timings.model = seconds_since(mark);
            timings.total = seconds_since(started);
            write_status(status_path, GenotypeRunStatus::ModelIncomplete, local.fingerprint, options, library,
                         timings, local.built.input.blocks.empty() ? nullptr : &local.built.report, nullptr, nullptr,
                         nullptr, std::string("the structural model could not be built: ") + error.what());
            return GenotypeRunStatus::ModelIncomplete;
        }
    }
    const PreparedGenotypeModel& model = locus != nullptr ? locus->model() : *local.model;
    const LocusInputReport& locus_report = locus != nullptr ? locus->report() : local.built.report;
    timings.model = seconds_since(mark);
    if (g_watch != nullptr) {
        std::lock_guard<std::mutex> lock(g_watch->mutex);
        g_watch->fingerprint = fingerprint;
    }
    write_running_status("markers");

    if (model.blocks().size() < 2 || model.panel().template_count() < 2) {
        timings.total = seconds_since(started);
        write_status(status_path, GenotypeRunStatus::ModelIncomplete, fingerprint, options, library,
                     timings, &locus_report, nullptr, nullptr, nullptr,
                     "the locus does not carry at least two blocks and two panel templates");
        return GenotypeRunStatus::ModelIncomplete;
    }
    // M4-B: the run's ONE copying-transition vector, before any stage computes a prior.
    if (locus == nullptr) g_geometry = local.built.geometry;
    try {
        CopyingTransitionSpec spec;
        spec.mode = options.copying_transition;
        spec.constant_r = options.switch_probability;
        spec.copying_length_bp = options.copying_length_bp;
        g_transitions = build_copying_transitions(model.blocks().size(), spec.mode == "constant_v1" ? nullptr : &g_geometry, spec);
    } catch (const UnsupportedCopyingTransition& error) {
        timings.total = seconds_since(started);
        write_status(status_path, GenotypeRunStatus::ModelIncomplete, fingerprint, options, library, timings, &locus_report, nullptr,
                     nullptr, nullptr, std::string("unsupported copying transition: ") + error.what());
        return GenotypeRunStatus::ModelIncomplete;
    }

    // ---- stage 2: markers ----------------------------------------------------------------------
    mark = Clock::now();
    const std::vector<ReadPair> fragments =
        read_fragments(options.reads1_path, options.reads2_path);
    timings.detail.emplace_back("markers_read_fastq", seconds_since(mark));
    Clock::time_point marker_lap = Clock::now();
    if (auto_library) {
        LibraryEstimateOptions estimate_options;
        estimate_options.threads = options.threads;
        estimate_options.marker_kmer = options.kmer_size;
        estimate_options.marker_syncmer = options.syncmer_s;
        const LibraryEstimate estimate = estimate_library(model, fragments, estimate_options, nullptr);
        const AutoLibrary built_library = build_auto_library(estimate, genotype_preset_library("short-read"), estimate_options);
        atomic_write(options.out_prefix + ".library_estimate.tsv",
                     library_estimate_table(estimate, nullptr) + (built_library.accepted ? library_used_rows(built_library) : ""));
        if (!built_library.accepted) {
            timings.markers = seconds_since(mark);
            timings.total = seconds_since(started);
            library.provenance = "auto: refused";
            write_status(status_path, GenotypeRunStatus::LibraryRefused, fingerprint, options, library, timings, &locus_report, nullptr, nullptr,
                         nullptr, built_library.refusal);
            if (!options.quiet) std::cerr << "[genotype] LIBRARY_REFUSED  " << built_library.refusal << '\n';
            return GenotypeRunStatus::LibraryRefused;
        }
        library = built_library.library;
    }
    // library_sensitivity_audit_v1: the accepted auto library at 0.9x and 1.1x its estimated fragment-start intensity,
    // every other value identical. Scored beside it in every combination; never behind the reported genotype.
    std::vector<LibraryModel> audit_libraries;
    if (auto_library && options.library_sensitivity_audit) {
        for (const double factor : {0.9, 1.1}) {
            LibraryModel endpoint = library;
            endpoint.fragment_start_intensity = factor * library.fragment_start_intensity;
            audit_libraries.push_back(endpoint);
        }
    }
    if (capture != nullptr) capture->library = library;
    timings.detail.emplace_back("markers_library_estimate", seconds_since(marker_lap));
    marker_lap = Clock::now();
    std::optional<PreparedMarkerPanel> local_panel;   // I3-A: built here unless the prepared locus holds it
    if (locus == nullptr) local_panel.emplace(build_marker_panel(model, options.kmer_size, options.syncmer_s));
    const PreparedMarkerPanel& panel = locus != nullptr ? locus->panel() : *local_panel;
    timings.detail.emplace_back("markers_panel_build", seconds_since(marker_lap));
    marker_lap = Clock::now();
    const MarkerObservations observations = count_marker_observations(panel, fragments);
    timings.detail.emplace_back("markers_observation_counting", seconds_since(marker_lap));
    marker_lap = Clock::now();
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
                     library, timings, &locus_report, &panel, nullptr, nullptr,
                     std::string("marker depth could not be estimated: ") + error.what());
        return GenotypeRunStatus::CandidateIncomplete;
    }
    if (!depth.converged) {
        // A depth estimate that did not converge is not fed to the proposal: it changes which candidates
        // reach K. NA19240 used to stop at its iteration cap with a rate 3.4% high and still proceed.
        timings.markers = seconds_since(mark);
        timings.total = seconds_since(started);
        write_status(status_path, GenotypeRunStatus::CandidateIncomplete, fingerprint, options, library,
                     timings, &locus_report, &panel, nullptr, nullptr,
                     "marker depth fit did not converge (" + depth.termination +
                         "); declare marker_background to use the one-dimensional fit",
                     nullptr, &depth);
        return GenotypeRunStatus::CandidateIncomplete;
    }
    // The fitted or declared background, unchanged: the fit and the emission share one model floor. This
    // used to raise a fitted 1e-9 to 1e-6 here, a log(1000) change per zero-copy marker observation.
    timings.detail.emplace_back("markers_depth_fit", seconds_since(marker_lap));
    marker_lap = Clock::now();
    const MarkerModelParameters marker_parameters =
        proposal_marker_parameters(depth, options.kmer_size, options.syncmer_s);
    ProposalReuseOptions proposal_reuse;
    proposal_reuse.enabled = options.marker_proposal_reuse_table_bytes != 0;
    proposal_reuse.table_budget_bytes = options.marker_proposal_reuse_table_bytes;
    const PreparedProposalEmissions emissions =
        build_proposal_emissions(panel, observations, marker_parameters, proposal_reuse, &g_proposal_reuse);
    timings.detail.emplace_back("markers_proposal_emissions", seconds_since(marker_lap));
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
                     timings, &locus_report, &panel, nullptr, nullptr,
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

    // ---- S2-A: the early evaluated incumbent (early_incumbent experimental_v1) ------------------------------
    // A checkpoint is a PROVISIONAL sequence_first_v2 report of the incumbent: content files first, the locus commit
    // record last, under the watch mutex (the watchdog sees the previous checkpoint or this one, never a torn one).
    const auto checkpoint_report = [&](const std::string& completion, const std::string& note) {
        std::vector<std::string> block_ids;
        for (const auto& block : model.blocks()) block_ids.push_back(block.block_id);
        SequenceReport r = provisional_sequence_report(
            g_incumbent.first, g_incumbent.second, g_incumbent.total, completion, note, block_ids,
            [&](const std::vector<std::size_t>& labels) { return spell_mosaic(model, MosaicPath{labels}).bases; });
        if (!g_incumbent.templates_spelled) {   // the active panel, spelled once per run (as annotate_v2 spells it)
            const PanelAlleleMatrix& active = model.panel();
            for (std::size_t t = 0; t < active.template_count(); ++t) {
                g_incumbent.template_sequences.push_back(spell_mosaic(model, MosaicPath{active.template_alleles[t]}).bases);
                g_incumbent.template_reverse.push_back(t < model.template_source_orientations().size() &&
                                                       model.template_source_orientations()[t] == ChainOrientation::Reverse);
            }
            g_incumbent.templates_spelled = true;
        }
        r.schema = "sequence_first_v2";
        annotate_panel_membership(r, model.panel().template_names, g_incumbent.template_sequences, g_incumbent.template_reverse);
        for (const auto& item : g_input_hashes) {
            if (item.first == "gfa") r.graph_sha256 = item.second;
        }
        r.model_fingerprint = fingerprint;
        SequenceCheckpoint& c = r.checkpoint;
        c.present = true;
        c.stage = g_incumbent.stage;
        c.search_status = completion == "checkpoint_search_in_progress" ? "in_progress" : completion;
        c.selection = g_incumbent.selection;
        c.evaluated_candidates = g_incumbent.evaluated;
        c.commits = g_incumbent.commits;
        c.log_likelihood = g_incumbent.score.log_likelihood;
        c.fragment_term = g_incumbent.score.fragment_term;
        c.exposure_term = g_incumbent.score.exposure_term;
        c.sequence_prior = g_incumbent.sequence_prior;
        c.total_score = g_incumbent.total;
        c.fragments_with_mass = g_incumbent.score.fragments_with_mass;
        c.tied_groups = g_incumbent.tied_groups;
        c.projection_ambiguous = g_incumbent.projection_ambiguous;
        c.explored_labellings = g_incumbent.explored_labellings;
        c.seconds_first = g_incumbent.seconds_first;
        c.seconds_commit = g_incumbent.seconds_commit;
        return r;
    };
    // Commits the held incumbent; false when the run was already claimed or the publication failed (then the
    // failure is recorded, the failed report is the run's sequence report, and the caller ends the run).
    const auto commit_checkpoint = [&]() -> bool {
        const double now = seconds_since(started);
        ++g_incumbent.commits;
        g_incumbent.seconds_commit = now;
        if (std::isnan(g_incumbent.seconds_first)) g_incumbent.seconds_first = now;
        SequenceReport r = checkpoint_report("checkpoint_search_in_progress",
                                             "a fully evaluated candidate, the best among those evaluated so far; the "
                                             "search continues");
        SequenceLegacy legacy{"RUNNING", "pending", "NA", options.dosage_mode, dosage_mode_source(options)};
        if (!g_dosage_stage_follows) {
            legacy.dosage_mode = "not_run";
            legacy.dosage_mode_source = "caller_only_api_without_dosage_stage";
        }
        std::unique_lock<std::mutex> lock;
        if (g_watch != nullptr) {
            lock = std::unique_lock<std::mutex>(g_watch->mutex);
            if (g_watch->finished) {
                --g_incumbent.commits;
                return false;
            }
        }
        compound_unassessed(r, model, fingerprint, options, "NOT_ASSESSED_PROVISIONAL");   // M9-J (opt-in)
        if (!publish_sequence_report(r, legacy, options.out_prefix, "")) {
            g_incumbent.failed = true;
            g_incumbent.failure = "checkpoint_publication_failed: " + r.publication_error;
            g_sequence_report = r;   // the failed report (exit 2); its failure record is on disk
            g_sequence_published = true;
            if (g_watch != nullptr) {
                g_watch->checkpoint_committed = false;
                g_watch->checkpoint_failed = true;
            }
            return false;
        }
        std::ifstream committed(options.out_prefix + ".sequence_locus.tsv");
        std::stringstream text;
        text << committed.rdbuf();
        if (g_watch != nullptr) {
            g_watch->checkpoint_committed = true;
            g_watch->checkpoint_stage = g_incumbent.stage;
            g_watch->checkpoint_locus = text.str();
        }
        g_incumbent.dirty = false;
        return true;
    };
    if (g_incumbent.enabled) {
        write_running_status("early_incumbent");
        g_incumbent.seconds_preparation = seconds_since(started);
        // Selection: the panel pair maximising the existing proposal marker emissions plus the exact path prior (the
        // ranking's own prior arithmetic), ties to panel order. Table lookups only: nothing is fragment-scored here.
        const Clock::time_point select_mark = Clock::now();
        const std::vector<double> early_switch = run_switches(model.blocks().size());   // M4-B
        std::map<std::vector<std::size_t>, double> early_path_prior;
        const auto path_prior_of = [&](const MosaicPath& path) {
            const auto found = early_path_prior.find(path.alleles);
            if (found != early_path_prior.end()) return found->second;
            ++g_incumbent.prior_evaluations;
            const double value = log_haploid_mosaic_prior(model.panel(), path, early_switch);
            early_path_prior.emplace(path.alleles, value);
            return value;
        };
        std::size_t chosen = 0;
        double best = -std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < panel_candidates.size(); ++i) {
            const MosaicDiplotype& d = panel_candidates[i].diplotype;
            double value = 0.0;
            for (std::size_t b = 0; b < model.blocks().size(); ++b) value += emissions.blocks[b].at(d.first.alleles[b], d.second.alleles[b]);
            const double pa = path_prior_of(d.first), pb = path_prior_of(d.second);
            const double negative_infinity = -std::numeric_limits<double>::infinity();
            value = pa == negative_infinity || pb == negative_infinity
                        ? negative_infinity
                        : value + pa + pb + (d.first != d.second ? std::log(2.0) : 0.0);
            ++g_incumbent.pairs_examined;
            if (i == 0 || value > best) {
                best = value;
                chosen = i;
            }
        }
        g_incumbent.selection = "argmax_proposal_marker_emissions_plus_exact_path_prior;panel_pairs_examined=" +
                                std::to_string(g_incumbent.pairs_examined) + ";ties_to_panel_order;not_a_winner_claim";
        g_incumbent.seconds_selection = seconds_since(select_mark);
        if (options.test_interrupt_point == "before_first_evaluation") {
            stop_run("TIME_LIMIT", "test hook: before the first candidate evaluation completed", 124);
        }
        if (options.test_interrupt_point == "footprint_before_first_evaluation") {   // M15-D.1 TEST HOOK
            stop_run("RESOURCE_REFUSAL", "test hook: footprint limit before the first candidate evaluation completed", 125);
        }
        // Full evaluation of that one pair under the unchanged scorer: every fragment, the admitted placements and
        // paired geometry, exposure (through the union pass's own route and caps), and the ranking's sequence prior.
        const Clock::time_point evaluate_mark = Clock::now();
        const MosaicDiplotype& early = panel_candidates[chosen].diplotype;
        const HaplotypeSequence ea = spell_mosaic(model, early.first), eb = spell_mosaic(model, early.second);
        const InternedHaplotypes early_interned = intern_haplotypes({early.first, early.second}, {ea.bases, eb.bases});
        const std::vector<HaplotypeCandidate> early_pair{
            {std::min(early_interned.index_of[0], early_interned.index_of[1]),
             std::max(early_interned.index_of[0], early_interned.index_of[1])}};
        FragmentScoringOptions early_options;
        early_options.threads = options.threads;
        early_options.mismatch_band = options.mismatch_band;
        early_options.minimum_piece_bp = options.minimum_piece_bp;
        early_options.max_piece_hits = options.max_piece_hits;
        early_options.max_proposed_starts_per_haplotype = options.max_proposed_starts_per_haplotype;
        early_options.max_physical_origins = options.max_physical_origins;
        early_options.max_state_classes = options.max_state_classes;
        early_options.max_mass_cells = options.max_mass_cells;
        early_options.emission_cache_bytes = options.emission_cache_bytes;
        early_options.prepared_library_terms = options.combination_prepared_library_terms;
        early_options.gap_recursion_backend = gap_recursion_backend_of(options);
        early_options.union_prefilter_cache_bytes = options.union_prefilter_cache_bytes;
        early_options.start_sort_scratch_bytes = options.start_sort_scratch_bytes;
        early_options.max_union_cached_prefilters_consumed = options.max_union_cached_prefilters_consumed;
        try {
            FragmentScoringResult early_scored;
            if (options.union_scoring_backend == "indexed_fine_experimental") {
                std::vector<std::vector<std::uint32_t>> segments(early_interned.sequences.size());
                for (std::size_t slot = 0; slot < 2; ++slot) {
                    const HaplotypeSequence& spelled = slot == 0 ? ea : eb;
                    std::vector<std::uint32_t>& lengths = segments[early_interned.index_of[slot]];
                    if (!lengths.empty()) continue;
                    std::size_t cut = 0;
                    for (std::size_t b = 1; b < spelled.block_intervals.size(); ++b) {
                        const std::size_t next = std::max(cut, std::min(spelled.block_intervals[b].begin, spelled.bases.size()));
                        lengths.push_back(static_cast<std::uint32_t>(next - cut));
                        cut = next;
                    }
                    lengths.push_back(static_cast<std::uint32_t>(spelled.bases.size() - cut));
                }
                // the union call's own indexed settings (score identity is gated by the union check below)
                early_options.max_union_cached_placements_consumed = options.max_union_cached_placements_consumed;
                early_options.union_origin_cache_bytes = options.union_origin_cache_bytes;
                early_options.max_union_cached_origins_consumed = options.max_union_cached_origins_consumed;
                early_options.max_union_cached_states_consumed = options.max_union_cached_states_consumed;
                early_scored = score_diplotypes_indexed_union(early_interned, early_pair, fragments, library, early_options,
                                                              segments, options.indexed_union_store_bytes, nullptr, nullptr,
                                                              options.indexed_union_store_layout == "legacy_i2a",
                                                              options.indexed_union_evidence_layout == "fine_views_experimental");
            } else {
                early_scored = score_diplotypes(early_interned, early_pair, fragments, library, early_options);
            }
            g_incumbent.early_counters = early_scored.counters;
            g_incumbent.early_fragment_evaluations = early_scored.candidate_fragment_evaluations;
            const double ha = log_haploid_sequence_prior(model, ea.bases, early_switch);
            const double hb = log_haploid_sequence_prior(model, eb.bases, early_switch);
            g_incumbent.have = true;
            g_incumbent.first = early.first.alleles;
            g_incumbent.second = early.second.alleles;
            g_incumbent.score = early_scored.scores.front();
            g_incumbent.sequence_prior = ha + hb + (ea.bases != eb.bases ? std::log(2.0) : 0.0);
            g_incumbent.total = g_incumbent.score.log_likelihood + g_incumbent.sequence_prior;
            g_incumbent.stage = "early_panel_incumbent";
            g_incumbent.evaluated = 1;
            g_incumbent.dirty = true;
            g_incumbent.early_status = "evaluated";
            if (options.test_perturb_early_incumbent_score) {   // S2-B TEST HOOK: an inconsistency the union check must catch
                g_incumbent.score.log_likelihood =
                    std::nextafter(g_incumbent.score.log_likelihood, std::numeric_limits<double>::infinity());
                g_incumbent.total = g_incumbent.score.log_likelihood + g_incumbent.sequence_prior;
            }
            g_incumbent.early_first = g_incumbent.first;
            g_incumbent.early_second = g_incumbent.second;
            g_incumbent.early_score = g_incumbent.score;
            g_incumbent.early_prior = g_incumbent.sequence_prior;
        } catch (const ResourceExhausted& error) {
            g_incumbent.early_status = std::string("refused_at_declared_cap: ") + error.what();
        }
        g_incumbent.seconds_evaluation = seconds_since(evaluate_mark);
        g_incumbent.footprint_after_early = physical_footprint().current;
        if (g_incumbent.have && !commit_checkpoint() && g_incumbent.failed) {
            timings.total = seconds_since(started);
            write_status(status_path, GenotypeRunStatus::RescoreIncomplete, fingerprint, options, library, timings,
                         &locus_report, &panel, nullptr, nullptr,
                         "the early incumbent checkpoint could not be published (" + g_incumbent.failure +
                             "); an explicit failure, never a fallback",
                         nullptr, &depth, &candidate_report);
            return GenotypeRunStatus::RescoreIncomplete;
        }
        if (options.test_interrupt_point == "after_checkpoint") {
            stop_run("TIME_LIMIT", "test hook: immediately after the first checkpoint commit", 124);
        }
    }

    // ---- stage 3b: mosaic proposal --------------------------------------------------------------------------
    write_running_status("proposal");
    mark = Clock::now();
    const std::vector<double> switch_probability = run_switches(model.blocks().size());   // M4-B: the run's one vector
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
            if (g_incumbent.enabled && g_incumbent.have && !g_incumbent.failed) {   // S2-A: a declared cap stop
                g_sequence_report = checkpoint_report("proposal_cap", std::string("the mosaic proposal stopped at a declared "
                                                      "resource budget (") + error.what() + "); the best fully evaluated "
                                                      "candidate so far; search incomplete");
            }
            write_status(status_path, GenotypeRunStatus::CandidateIncomplete, fingerprint, options, library,
                         timings, &locus_report, &panel, nullptr, nullptr,
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
    std::vector<std::vector<std::uint32_t>> segment_lengths_of_slot;   // block segmentation: a reuse hint only
    std::map<std::vector<std::size_t>, std::size_t> slot_of_path;
    const auto slot_for = [&](const MosaicPath& path) {
        const auto found = slot_of_path.find(path.alleles);
        if (found != slot_of_path.end()) return found->second;
        const std::size_t slot = paths.size();
        slot_of_path.emplace(path.alleles, slot);
        paths.push_back(path);
        HaplotypeSequence spelled = spell_mosaic(model, path);
        // A partition from the block starts, so the lengths always sum to the haplotype (a hint either way).
        std::vector<std::uint32_t> lengths;
        std::size_t cut = 0;
        for (std::size_t b = 1; b < spelled.block_intervals.size(); ++b) {
            const std::size_t next = std::max(cut, std::min(spelled.block_intervals[b].begin, spelled.bases.size()));
            lengths.push_back(static_cast<std::uint32_t>(next - cut));
            cut = next;
        }
        lengths.push_back(static_cast<std::uint32_t>(spelled.bases.size() - cut));
        segment_lengths_of_slot.push_back(std::move(lengths));
        sequences.push_back(std::move(spelled.bases));
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
    if (options.test_interrupt_point == "union_scoring") {   // S2-A TEST HOOK: the union pass has not completed
        stop_run("TIME_LIMIT", "test hook: during union scoring (before the union pass completed)", 124);
    }
    mark = Clock::now();
    FragmentScoringOptions scoring_options;
    // I2-K (EXPERIMENTAL): owned here for this sample only; destroyed on every exit, adopted at most once
    EmissionCacheHandoff union_emissions;
    scoring_options.audit_libraries = audit_libraries;
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
    scoring_options.prepared_library_terms = options.combination_prepared_library_terms;   // I2-F
    scoring_options.gap_recursion_backend = gap_recursion_backend_of(options);   // I2-H
    // Exact union placement reuse: one segmentation per unique haplotype (its first path's blocks).
    std::vector<std::vector<std::uint32_t>> union_segments;
    if (options.union_placement_cache_bytes > 0 && any_ranked) {
        union_segments.assign(interned.sequences.size(), {});
        std::vector<char> filled(interned.sequences.size(), 0);
        for (std::size_t slot = 0; slot < paths.size(); ++slot) {
            const std::size_t h = interned.index_of[slot];
            if (filled[h] != 0) continue;
            filled[h] = 1;
            union_segments[h] = segment_lengths_of_slot[slot];
        }
        scoring_options.union_segment_lengths = &union_segments;
        scoring_options.union_placement_cache_bytes = options.union_placement_cache_bytes;
        scoring_options.max_union_cached_placements_consumed = options.max_union_cached_placements_consumed;
        scoring_options.union_origin_cache_bytes = options.union_origin_cache_bytes;
        scoring_options.max_union_cached_origins_consumed = options.max_union_cached_origins_consumed;
        scoring_options.max_union_cached_states_consumed = options.max_union_cached_states_consumed;
    }
    // I2-A: the indexed backend takes the same segmentation (first path's blocks per unique haplotype) as its own
    // argument; scoring_options -- which the repair stage copies -- is not touched.
    const bool indexed_union = options.union_scoring_backend == "indexed_fine_experimental";
    std::vector<std::vector<std::uint32_t>> indexed_segments;
    if (indexed_union && any_ranked) {
        indexed_segments.assign(interned.sequences.size(), {});
        std::vector<char> filled(interned.sequences.size(), 0);
        for (std::size_t slot = 0; slot < paths.size(); ++slot) {
            const std::size_t h = interned.index_of[slot];
            if (filled[h] != 0) continue;
            filled[h] = 1;
            indexed_segments[h] = segment_lengths_of_slot[slot];
        }
    }
    scoring_options.union_prefilter_cache_bytes = options.union_prefilter_cache_bytes;
    scoring_options.start_sort_scratch_bytes = options.start_sort_scratch_bytes;
    scoring_options.max_union_cached_prefilters_consumed = options.max_union_cached_prefilters_consumed;
    if (!options.diagnostic_refinement_replay.empty()) {
        // ---- S2-H DIAGNOSTIC ONLY: replay the first refinement round from a recorded anchor (no union, no early pass) ----
        const auto record = [&](const std::string& file) {
            std::map<std::string, std::string> kv;
            std::ifstream in(file);
            if (!in) throw std::runtime_error("diagnostic_refinement_replay: cannot read " + file);
            for (std::string line; std::getline(in, line);) {
                const std::size_t tab = line.find('\t');
                if (tab != std::string::npos) kv[line.substr(0, tab)] = line.substr(tab + 1);
            }
            return kv;
        };
        const std::map<std::string, std::string> in = record(options.diagnostic_refinement_replay);
        RepairSetup setup = make_repair_setup(options, scoring_options);
        // every identity first: model, inputs by content, proposal emissions, the refinement setup. S2-J: the ONE declared
        // override -- the round-1 origin-admission policy -- is compared at the recorded run's setting; every other field
        // of the setup must match exactly, and both setup texts are recorded.
        const std::string effective_setup = repair_setup_text(setup, fragments.size());
        RepairSetup comparable = setup;
        comparable.search.round1_reuse_only_admission = false;
        const std::string comparable_setup = repair_setup_text(comparable, fragments.size());
        const std::string admission_override = setup.search.round1_reuse_only_admission
                                                   ? std::string("repair_round1_origin_admission existing->reuse_only_experimental")
                                                   : std::string("none");
        std::vector<std::pair<std::string, std::string>> expected = {
            {"schema", "diagnostic_refinement_replay_inputs_v1"}, {"model_fingerprint", fingerprint},
            {"emissions_digest", emissions_digest(emissions)}, {"repair_setup", comparable_setup}};
        for (const auto& hash : g_input_hashes) expected.push_back({"sha256_" + hash.first, hash.second});
        std::string identities;
        for (const auto& [k, v] : expected) {
            const auto found = in.find(k);
            if (found == in.end() || found->second != v) {
                throw std::runtime_error("diagnostic_refinement_replay: the recorded " + k + " does not match this run (" +
                                         (found == in.end() ? std::string("absent") : found->second) + " vs " + v + ")");
            }
            identities += (identities.empty() ? "" : ",") + k;
        }
        const auto alleles = [&](const std::string& key) {
            std::vector<std::size_t> out;
            std::stringstream text(in.at(key));
            for (std::string item; std::getline(text, item, ',');) out.push_back(std::stoull(item));
            if (out.size() != model.blocks().size()) throw std::runtime_error("diagnostic_refinement_replay: " + key + " has the wrong block count");
            return out;
        };
        const std::vector<std::size_t> hap1 = alleles("anchor_hap1"), hap2 = alleles("anchor_hap2");
        const double anchor = std::strtod(in.at("anchor_score").c_str(), nullptr);
        if (hex_double(anchor) != in.at("anchor_score")) throw std::runtime_error("diagnostic_refinement_replay: the anchor score does not round-trip");
        setup.search.max_rounds = 1;   // declared: the first round only
        write_running_status("refinement replay");
        note_footprint("replay_refinement_context_start");
        CoordinateRepairSearch search(model, emissions, fragments, library, setup.scoring, setup.search);
        search.start(hap1, hap2, anchor);
        std::ostringstream scored_rows;
        scored_rows << "item\thap1\thap2\tlog_likelihood\tfragment_term\texposure_term\tabsolute_term_sum\tfragments_with_mass\tprior\n";
        std::size_t item = 0;
        RepairHooks hooks;
        hooks.on_checkpoint = [](const char* point) { note_footprint(point); };
        hooks.on_scored = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b, const std::string&,
                              const std::string&, const DiplotypeScore& score, double prior) {
            scored_rows << item++ << '\t' << join_alleles(a) << '\t' << join_alleles(b) << '\t' << hex_double(score.log_likelihood)
                        << '\t' << hex_double(score.fragment_term) << '\t' << hex_double(score.exposure_term) << '\t'
                        << hex_double(score.absolute_term_sum) << '\t' << score.fragments_with_mass << '\t' << hex_double(prior) << '\n';
        };
        const auto replay_mark = Clock::now();
        const bool reproduced = search.run(hooks);
        const double replay_seconds = seconds_since(replay_mark);
        note_footprint("replay_round_done");
        atomic_write(options.out_prefix + ".refinement_replay_rounds.tsv", search.render_rounds(search.rounds()));
        atomic_write(options.out_prefix + ".refinement_replay_scored.tsv", scored_rows.str());
        const FootprintSample fp = physical_footprint();
        std::ostringstream t;
        t << "key\tvalue\nschema\tdiagnostic_refinement_replay_v1\nnote\tDIAGNOSTIC ONLY: the first refinement round replayed "
             "from a recorded anchor in a fresh process; no union pass, no genotype\ninputs\t" << options.diagnostic_refinement_replay
          << "\nidentities_checked\t" << identities << "\ndeclared_setup_override\t" << admission_override
          << "\nrepair_setup_recorded\t" << in.at("repair_setup") << "\nrepair_setup_effective\t" << effective_setup
          << "\nanchor_declared\t" << in.at("anchor_score")
          << "\nanchor_reproduced\t" << (reproduced ? "yes" : "NO") << "\nanchor_scored\t"
          << (search.rounds().empty() ? std::string("NA") : hex_double(search.rounds().front().start_total))
          << "\nrounds_completed\t" << search.rounds().size() << "\nstop_reason\t" << search.stop_reason()
          << "\ncandidates_scored\t" << item << "\nseconds_replay\t" << replay_seconds << "\nphys_footprint_bytes\t" << fp.current
          << "\nphys_footprint_peak_bytes\t" << fp.peak << '\n';
        atomic_write(options.out_prefix + ".refinement_replay.tsv", t.str());
        timings.total = seconds_since(started);
        write_status(status_path, GenotypeRunStatus::RefinementReplayDiagnosticOnly, fingerprint, options, library, timings,
                     &locus_report, &panel, nullptr, nullptr,
                     "DIAGNOSTIC ONLY: the first refinement round was replayed from a recorded anchor; no union pass and no "
                     "genotype", nullptr, &depth, nullptr);
        return GenotypeRunStatus::RefinementReplayDiagnosticOnly;
    }
    FragmentScoringResult scored;
    std::string rescore_failure;
    std::vector<HaplotypeWork> union_work;
    if (options.diagnostic_union_work_trace) scoring_options.record_haplotype_work = &union_work;
    // DIAGNOSTIC ONLY: a seeding-profile run replays the seeding path instead of scoring the union at all.
    const bool seeding_profile_only =
        options.diagnostic_union_seeding_profile > 0 || !options.diagnostic_union_verification_census.empty();
    if (seeding_profile_only && !options.diagnostic_union_only) {
        throw std::invalid_argument("the union seeding profile and census are union-only diagnostics (set diagnostic_union_only)");
    }
    // S2-C/S2-D: a truth-free proposal score per sequence pair (the S2-A selection formula: proposal marker emissions plus
    // the exact path prior, maximised over the pair's labellings) -- the progressive schedule's order and the kernel dump's
    const auto pair_proposal_scores = [&]() {
        const std::vector<double> sw = run_switches(model.blocks().size());   // M4-B
        std::map<std::vector<std::size_t>, double> path_prior;
        const auto prior_of = [&](const MosaicPath& path) {
            const auto found = path_prior.find(path.alleles);
            if (found != path_prior.end()) return found->second;
            return path_prior.emplace(path.alleles, log_haploid_mosaic_prior(model.panel(), path, sw)).first->second;
        };
        std::vector<double> proposal(scored_pairs.size(), -std::numeric_limits<double>::infinity());
        for (const MosaicCandidate* candidate : scored_union) {
            const MosaicDiplotype& d = candidate->diplotype;
            double value = 0.0;
            for (std::size_t b = 0; b < model.blocks().size(); ++b) value += emissions.blocks[b].at(d.first.alleles[b], d.second.alleles[b]);
            const double pa = prior_of(d.first), pb = prior_of(d.second);
            const double ninf = -std::numeric_limits<double>::infinity();
            value = pa == ninf || pb == ninf ? ninf : value + pa + pb + (d.first != d.second ? std::log(2.0) : 0.0);
            double& slot = proposal[scored_index.at(sequence_pair(*candidate))];
            slot = std::max(slot, value);
        }
        return proposal;
    };
    if (options.diagnostic_union_kernel_dump && any_ranked && indexed_union) {
        // S2-C (DIAGNOSTIC ONLY): the union pass's exact inputs, for the progressive-kernel harness
        std::string bin;
        const auto put = [&](const auto& v) { bin.append(reinterpret_cast<const char*>(&v), sizeof v); };
        const auto put_library = [&](const LibraryModel& l) {
            put(l.insert_mean_bp); put(l.insert_sd_bp); put(l.support_sigmas); put(l.substitution_error);
            put(l.outlier_weight); put(l.fragment_start_intensity); put(l.gap_open); put(l.gap_extend);
            put(static_cast<std::uint64_t>(l.max_gap_bp));
        };
        bin.append("PVKI1\0\0\0", 8);
        put(static_cast<std::uint64_t>(interned.sequences.size()));
        for (std::size_t h = 0; h < interned.sequences.size(); ++h) {
            put(static_cast<std::uint64_t>(interned.sequences[h].size()));
            bin.append(interned.sequences[h]);
            put(static_cast<std::uint64_t>(indexed_segments[h].size()));
            for (const std::uint32_t x : indexed_segments[h]) put(x);
        }
        const std::vector<double> proposal = pair_proposal_scores();
        put(static_cast<std::uint64_t>(scored_pairs.size()));
        for (std::size_t i = 0; i < scored_pairs.size(); ++i) {
            const HaplotypeCandidate& c = scored_pairs[i];
            put(static_cast<std::uint64_t>(c.first));
            put(static_cast<std::uint64_t>(c.second));
            put(proposal[i]);
            put(haploid_sequence_prior[c.first] + haploid_sequence_prior[c.second] + (c.first != c.second ? std::log(2.0) : 0.0));
        }
        std::uint64_t early_first = std::numeric_limits<std::uint64_t>::max(), early_second = early_first;
        if (!g_incumbent.early_first.empty()) {
            std::size_t x = interned.index_of[slot_of_path.at(g_incumbent.early_first)];
            std::size_t y = interned.index_of[slot_of_path.at(g_incumbent.early_second)];
            if (y < x) std::swap(x, y);
            early_first = x;
            early_second = y;
        }
        put(early_first);
        put(early_second);
        put_library(library);
        const FragmentScoringOptions& o = scoring_options;
        for (const std::uint64_t x : {static_cast<std::uint64_t>(o.threads), static_cast<std::uint64_t>(o.mismatch_band),
                                      static_cast<std::uint64_t>(o.minimum_piece_bp), o.max_piece_hits, o.piece_hit_publish_batch,
                                      o.max_proposed_starts_per_haplotype, o.max_physical_origins, o.max_state_classes,
                                      o.max_mass_cells, o.max_latent_alignment_pairs, static_cast<std::uint64_t>(o.emission_cache_bytes),
                                      o.max_union_cached_placements_consumed, static_cast<std::uint64_t>(o.union_origin_cache_bytes),
                                      o.max_union_cached_origins_consumed, o.max_union_cached_states_consumed,
                                      static_cast<std::uint64_t>(o.union_prefilter_cache_bytes), o.max_union_cached_prefilters_consumed,
                                      static_cast<std::uint64_t>(o.start_sort_scratch_bytes),
                                      static_cast<std::uint64_t>(o.prepared_library_terms ? 1 : 0),
                                      static_cast<std::uint64_t>(o.gap_recursion_backend == GapRecursionBackend::scaled_reference ? 0 : 1),
                                      static_cast<std::uint64_t>(options.indexed_union_store_bytes),
                                      static_cast<std::uint64_t>(options.indexed_union_store_layout == "legacy_i2a" ? 1 : 0),
                                      static_cast<std::uint64_t>(options.indexed_union_evidence_layout == "fine_views_experimental" ? 1 : 0),
                                      static_cast<std::uint64_t>(o.audit_libraries.size())}) {
            put(x);
        }
        for (const LibraryModel& l : o.audit_libraries) put_library(l);
        atomic_write(options.out_prefix + ".union_kernel_inputs.bin", bin);
    }
    if (any_ranked && !seeding_profile_only) {
        try {
            candidate_report.indexed_union = indexed_union;
            if (indexed_union) {
                // I2-A (EXPERIMENTAL): the same pass through the indexed route; same result structure and refusals
                FragmentScoringOptions indexed_options = scoring_options;
                if (options.adaptive_emission_handoff) indexed_options.emission_handoff_out = &union_emissions;   // I2-K
                indexed_options.union_segment_lengths = nullptr;
                indexed_options.union_placement_cache_bytes = 0;
                indexed_options.max_union_cached_placements_consumed = options.max_union_cached_placements_consumed;
                indexed_options.union_origin_cache_bytes = options.union_origin_cache_bytes;
                indexed_options.max_union_cached_origins_consumed = options.max_union_cached_origins_consumed;
                indexed_options.max_union_cached_states_consumed = options.max_union_cached_states_consumed;
                if (options.union_fragment_blocks != "off") {
                    // ---- S2-G: the combined route -- the same indexed scorer over fragment blocks, each block owning its
                    // whole and fine-core aggregates; scores and audit views as the bulk route's (S2-F gates)
                    const std::size_t block = options.test_union_fragment_block_pairs != 0 ? options.test_union_fragment_block_pairs : 4096;
                    candidate_report.blocked_route = true;
                    candidate_report.blocked_block_pairs = block;
                    candidate_report.blocked_options = indexed_options;
                    BlockedUnionResult out = score_diplotypes_indexed_blocked(
                        interned, scored_pairs, fragments, library, indexed_options, indexed_segments, options.indexed_union_store_bytes,
                        true, block,
                        [&](BlockedBlockReport& rep, bool before) {
                            if (before) return;
                            rep.footprint_after = physical_footprint().current;
                            candidate_report.blocked_footprint_max = std::max(candidate_report.blocked_footprint_max, rep.footprint_after);
                        },
                        true);
                    candidate_report.indexed.index = out.index;
                    candidate_report.indexed.seconds_index_build = out.seconds_index_build;
                    candidate_report.indexed.store_budget_bytes = options.indexed_union_store_bytes;
                    candidate_report.indexed.refusal = out.refusal;
                    candidate_report.indexed.max_read_bp = 1;   // as the bulk adapter reports it
                    for (const ReadPair& f : fragments) {
                        candidate_report.indexed.max_read_bp =
                            std::max({candidate_report.indexed.max_read_bp, f.first.bases.size(), f.second.bases.size()});
                    }
                    double block_seconds = 0.0;
                    for (const BlockedBlockReport& b : out.blocks) block_seconds += b.seconds;
                    if (!out.completed) {
                        // no union score is published; the whole attempt's observed work is what the refusal reports
                        FragmentScoringObservation& o = candidate_report.observed;
                        o = FragmentScoringObservation{};
                        o.counters = out.attempt_totals.executed;
                        o.threads_used_placement = out.result.threads_used_placement;
                        o.threads_used_combination = out.result.threads_used_combination;
                        o.seconds_placement_and_mass = block_seconds - out.seconds_combination;
                        o.seconds_diploid_combination = out.seconds_combination;
                        o.union_origin_cache = out.result.union_origin_cache;
                        candidate_report.blocked = std::move(out);
                        note_footprint("union_refusal");
                        throw ResourceExhausted(candidate_report.blocked.refusal);   // the existing refusal path (the checkpoint stays)
                    }
                    FragmentScoringResult combined = std::move(out.result);
                    combined.seconds_placement_and_mass = block_seconds - out.seconds_combination;
                    out.result = FragmentScoringResult{};
                    candidate_report.blocked = std::move(out);
                    scored = std::move(combined);
                } else if (!g_progressive.enabled) {
                scored = score_diplotypes_indexed_union(interned, scored_pairs, fragments, library, indexed_options,
                                                        indexed_segments, options.indexed_union_store_bytes,
                                                        &candidate_report.indexed, &candidate_report.observed,
                                                        options.indexed_union_store_layout == "legacy_i2a",
                                                        options.indexed_union_evidence_layout == "fine_views_experimental");
                } else {
                    // ---- S2-D: the same union pass as a declared schedule on one selective batch (scheduling only)
                    g_progressive.seconds_union_start = seconds_since(started);
                    const std::vector<double> proposal = pair_proposal_scores();
                    const auto pair_of_labels = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
                        std::size_t x = interned.index_of[slot_of_path.at(a)], y = interned.index_of[slot_of_path.at(b)];
                        if (y < x) std::swap(x, y);
                        return scored_index.at({x, y});
                    };
                    const auto pair_prior = [&](std::size_t p) {
                        const HaplotypeCandidate& c = scored_pairs[p];
                        return haploid_sequence_prior[c.first] + haploid_sequence_prior[c.second] +
                               (c.first != c.second ? std::log(2.0) : 0.0);
                    };
                    std::vector<std::size_t> order(scored_pairs.size());
                    std::iota(order.begin(), order.end(), 0);
                    std::stable_sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return proposal[a] > proposal[b]; });
                    const std::size_t first_pair = !g_incumbent.early_first.empty()
                                                       ? pair_of_labels(g_incumbent.early_first, g_incumbent.early_second)
                                                       : order.front();
                    std::vector<std::vector<std::size_t>> stage_index(3);
                    stage_index[0].push_back(first_pair);
                    for (const std::size_t c : order) {
                        if (c != first_pair && stage_index[1].size() < 31) stage_index[1].push_back(c);
                    }
                    {
                        std::vector<char> placed(scored_pairs.size(), 0);
                        for (std::size_t st = 0; st < 2; ++st) for (const std::size_t c : stage_index[st]) placed[c] = 1;
                        for (std::size_t c = 0; c < scored_pairs.size(); ++c) if (!placed[c]) stage_index[2].push_back(c);
                    }
                    std::vector<std::vector<HaplotypeCandidate>> stages(3);
                    for (std::size_t st = 0; st < 3; ++st) for (const std::size_t c : stage_index[st]) stages[st].push_back(scored_pairs[c]);
                    std::vector<DiplotypeScore> pool_scores(scored_pairs.size());
                    std::vector<std::vector<DiplotypeScore>> pool_audits;
                    std::size_t stage_number = 0;
                    const auto on_stage = [&](ProgressiveStageReport& rep, const FragmentScoringResult& r) {
                        const FootprintSample fp = physical_footprint();
                        rep.footprint_bytes = fp.current;
                        rep.footprint_peak_bytes = fp.peak;
                        const std::vector<std::size_t>& idx = stage_index[stage_number];
                        if (pool_audits.size() != r.audit_scores.size()) {
                            pool_audits.assign(r.audit_scores.size(), std::vector<DiplotypeScore>(scored_pairs.size()));
                        }
                        for (std::size_t i = 0; i < idx.size(); ++i) {
                            pool_scores[idx[i]] = r.scores[i];
                            for (std::size_t k = 0; k < r.audit_scores.size(); ++k) pool_audits[k][idx[i]] = r.audit_scores[k][i];
                        }
                        g_progressive.stages.push_back(rep);
                        ++stage_number;
                        if (stage_number != 2) return;
                        // the 32-pair prefix checkpoint: its best fully evaluated sequence pair, with the prefix's numerical
                        // ties and the pair's projection ambiguity over its labellings in the union
                        std::vector<std::size_t> prefix = stage_index[0];
                        prefix.insert(prefix.end(), stage_index[1].begin(), stage_index[1].end());
                        std::vector<DiplotypeScore> ps;
                        std::vector<double> pp;
                        for (const std::size_t p : prefix) { ps.push_back(pool_scores[p]); pp.push_back(pair_prior(p)); }
                        const CertifiedRanking ranking = certify_scored_ranking(ps, fragments.size(), scoring_options.mismatch_band,
                                                                                FragmentScoringCounters{}, CertificationOptions{}, pp);
                        const std::size_t best = prefix[ranking.candidates.front().candidate_id];
                        g_progressive.prefix_pairs = prefix.size();
                        g_progressive.prefix_tied_groups = ranking.winner_set.size();
                        if (g_incumbent.enabled && !g_incumbent.failed) {
                            const std::vector<double> sw = run_switches(model.blocks().size());   // M4-B
                            const MosaicCandidate* representative = nullptr;
                            double representative_prior = -std::numeric_limits<double>::infinity();
                            std::size_t labellings = 0, prefix_labellings = 0;
                            std::vector<char> in_prefix(scored_pairs.size(), 0);
                            for (const std::size_t p : prefix) in_prefix[p] = 1;
                            std::vector<const MosaicCandidate*> members;
                            for (const MosaicCandidate* c : scored_union) {
                                const std::size_t p = scored_index.at(sequence_pair(*c));
                                prefix_labellings += in_prefix[p];
                                if (p != best) continue;
                                ++labellings;
                                members.push_back(c);
                                const double pa = log_haploid_mosaic_prior(model.panel(), c->diplotype.first, sw);
                                const double pb = log_haploid_mosaic_prior(model.panel(), c->diplotype.second, sw);
                                const double path = pa + pb + (c->diplotype.first != c->diplotype.second ? std::log(2.0) : 0.0);
                                if (representative == nullptr || path > representative_prior) { representative = c; representative_prior = path; }
                            }
                            bool ambiguous = false;
                            for (const MosaicCandidate* c : members) {
                                for (std::size_t b = 0; b < model.blocks().size(); ++b) {
                                    if (std::minmax(c->diplotype.first.alleles[b], c->diplotype.second.alleles[b]) !=
                                        std::minmax(representative->diplotype.first.alleles[b], representative->diplotype.second.alleles[b])) {
                                        ambiguous = true;
                                    }
                                }
                            }
                            const double total = pool_scores[best].log_likelihood + pair_prior(best);
                            g_incumbent.evaluated = std::max(g_incumbent.evaluated, prefix_labellings);
                            if (!g_incumbent.have || total >= g_incumbent.total) {
                                g_incumbent.have = true;
                                g_incumbent.first = representative->diplotype.first.alleles;
                                g_incumbent.second = representative->diplotype.second.alleles;
                                g_incumbent.score = pool_scores[best];
                                g_incumbent.sequence_prior = pair_prior(best);
                                g_incumbent.total = total;
                                g_incumbent.stage = "union_prefix_best";
                                g_incumbent.tied_groups = std::to_string(ranking.winner_set.size());
                                g_incumbent.projection_ambiguous = ambiguous ? "yes" : "no";
                                g_incumbent.explored_labellings = std::to_string(labellings);
                            }
                            if (commit_checkpoint()) {
                                g_progressive.prefix_committed = true;
                                g_progressive.seconds_prefix_checkpoint = seconds_since(started);
                            } else if (g_incumbent.failed) {
                                throw std::runtime_error("the progressive prefix checkpoint could not be published (" +
                                                         g_incumbent.failure + "); an explicit failure, never a fallback");
                            }
                        }
                        if (options.test_interrupt_point == "after_prefix_checkpoint") {   // S2-D TEST HOOK
                            stop_run("TIME_LIMIT", "test hook: right after the progressive 32-pair checkpoint", 124);
                        }
                    };
                    const ProgressiveResult out = score_diplotypes_indexed_progressive(
                        interned, stages, fragments, library, indexed_options, indexed_segments, options.indexed_union_store_bytes,
                        options.indexed_union_evidence_layout == "fine_views_experimental", on_stage, false);
                    candidate_report.indexed.index = out.index;
                    candidate_report.indexed.batch = out.batch;
                    candidate_report.indexed.seconds_index_build = out.seconds_index_build;
                    candidate_report.indexed.seconds_preparation = out.seconds_preparation;
                    candidate_report.indexed.store_budget_bytes = options.indexed_union_store_bytes;
                    candidate_report.indexed.max_read_bp = 1;   // as the bulk adapter reports it
                    for (const ReadPair& f : fragments) {
                        candidate_report.indexed.max_read_bp =
                            std::max({candidate_report.indexed.max_read_bp, f.first.bases.size(), f.second.bases.size()});
                    }
                    g_progressive.completed = out.completed;
                    if (!out.completed) {
                        g_progressive.refusal = out.refusal;
                        FragmentScoringObservation& o = candidate_report.observed;
                        o = FragmentScoringObservation{};
                        o.counters = out.summary.counters;
                        o.threads_used_placement = out.summary.threads_used_placement;
                        o.threads_used_combination = out.summary.threads_used_combination;
                        o.seconds_fragment_preparation = out.summary.seconds_fragment_preparation;
                        o.seconds_placement_and_mass = out.summary.seconds_placement_and_mass;
                        o.seconds_diploid_combination = out.summary.seconds_diploid_combination;
                        o.union_origin_cache = out.summary.union_origin_cache;
                        o.union_prefilter_cache = out.summary.union_prefilter_cache;
                        o.start_sort = out.summary.start_sort;
                        throw ResourceExhausted(out.refusal);   // the existing refusal path (the checkpoint stays)
                    }
                    FragmentScoringResult combined = out.summary;
                    combined.scores = std::move(pool_scores);
                    combined.audit_scores = std::move(pool_audits);
                    scored = std::move(combined);
                }
            } else {
                scored = score_diplotypes(interned, scored_pairs, fragments, library, scoring_options,
                                          &candidate_report.observed);
            }
            timings.rescoring_fragment_preparation = scored.seconds_fragment_preparation;
            timings.rescoring_placement_and_mass = scored.seconds_placement_and_mass;
            timings.rescoring_diploid_combination = scored.seconds_diploid_combination;
            candidate_report.threads_used_combination = scored.threads_used_combination;
            candidate_report.threads_used_placement = scored.threads_used_placement;
            candidate_report.candidate_fragment_evaluations = scored.candidate_fragment_evaluations;
            candidate_report.union_placement_cache = scored.union_placement_cache;
            candidate_report.union_origin_cache = scored.union_origin_cache;
            candidate_report.union_prefilter_cache = scored.union_prefilter_cache;
            candidate_report.start_sort = scored.start_sort;
            g_union_end_footprint = physical_footprint().current;
            note_footprint("union_completion");   // S2-G
        } catch (const ResourceExhausted& error) {
            rescore_failure = error.what();
            candidate_report.observed_valid = true;
            candidate_report.union_placement_cache = candidate_report.observed.union_placement_cache;
            candidate_report.union_origin_cache = candidate_report.observed.union_origin_cache;
            candidate_report.union_prefilter_cache = candidate_report.observed.union_prefilter_cache;
            candidate_report.start_sort = candidate_report.observed.start_sort;
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
    if (options.diagnostic_union_work_trace && any_ranked) {
        // Which unique haplotypes a panel template spells, and how many allele paths spell each one.
        std::set<std::vector<std::size_t>> panel_paths(model.panel().template_alleles.begin(),
                                                       model.panel().template_alleles.end());
        std::vector<std::size_t> paths_of(interned.sequences.size(), 0);
        std::vector<bool> is_panel(interned.sequences.size(), false);
        for (const auto& [alleles, slot] : slot_of_path) {
            const std::size_t h = interned.index_of[slot];
            ++paths_of[h];
            if (panel_paths.count(alleles) != 0) is_panel[h] = true;
        }
        std::ostringstream text;
        text << std::setprecision(9)
             << "haplotype\tclass\tallele_paths\trecorded\thaplotype_bp\tpiece_hits\tproposed_starts\tunique_starts\t"
                "full_read_verifications\tprefilter_rejections\tdp_evaluations\texecuted_dps\taccepted_mate_placements\t"
                "physical_origins\tstatistical_state_classes\tseconds_seeding_and_emission\tseconds_pairing_and_reduction\t"
                "reused_physical_origins\treused_state_aggregates\tmerged_state_classes\tseconds_cached_merge\t"
                "seconds_seed_scan\tseconds_start_sort\tseconds_verification\tseconds_replay\tseconds_row_sort\n";
        for (std::size_t h = 0; h < union_work.size(); ++h) {
            const HaplotypeWork& w = union_work[h];
            text << h << '\t' << (is_panel[h] ? "panel" : "mosaic") << '\t' << paths_of[h] << '\t'
                 << (w.recorded ? "true" : "false") << '\t' << w.haplotype_bp << '\t' << w.piece_hits << '\t'
                 << w.proposed_starts << '\t' << w.unique_starts << '\t' << w.full_read_verifications << '\t'
                 << w.prefilter_rejections << '\t' << w.dp_evaluations << '\t' << w.executed_dps << '\t'
                 << w.accepted_mate_placements << '\t' << w.physical_origins << '\t' << w.statistical_state_classes
                 << '\t' << w.seconds_seeding_and_emission << '\t' << w.seconds_pairing_and_reduction << '\t'
                 << w.reused_physical_origins << '\t' << w.reused_state_aggregates << '\t' << w.merged_state_classes
                 << '\t' << w.seconds_cached_merge << '\t' << w.seconds_seed_scan << '\t' << w.seconds_start_sort
                 << '\t' << w.seconds_verification << '\t' << w.seconds_replay << '\t' << w.seconds_row_sort << '\n';
        }
        atomic_write(options.out_prefix + ".union_haplotype_work.tsv", text.str());
    }
    if (!options.diagnostic_union_verification_census.empty() && any_ranked) {
        std::set<std::vector<std::size_t>> panel_paths(model.panel().template_alleles.begin(),
                                                       model.panel().template_alleles.end());
        std::vector<char> is_panel(interned.sequences.size(), 0);
        for (const auto& [alleles, slot] : slot_of_path) {
            if (panel_paths.count(alleles) != 0) is_panel[interned.index_of[slot]] = 1;
        }
        const UnionVerificationCensus census = census_union_verification(
            interned, fragments, library, scoring_options, is_panel, options.diagnostic_union_verification_census);
        std::ostringstream text;
        text << "# DIAGNOSTIC ONLY: exact census of the union's direct verification stream (phase-1 zones applied)\n"
             << "key\tpanel\tmosaic\n";
        const auto row = [&](const char* key, std::uint64_t UnionCensusClass::*field) {
            text << key << '\t' << census.panel.*field << '\t' << census.mosaic.*field << '\n';
        };
        row("haplotypes", &UnionCensusClass::haplotypes);
        row("requests", &UnionCensusClass::requests);
        row("ineligible", &UnionCensusClass::ineligible);
        row("accepted", &UnionCensusClass::accepted);
        row("rejected", &UnionCensusClass::rejected);
        row("first_accepted", &UnionCensusClass::first_accepted);
        row("first_rejected", &UnionCensusClass::first_rejected);
        row("repeated_accepted", &UnionCensusClass::repeated_accepted);
        row("repeated_rejected", &UnionCensusClass::repeated_rejected);
        row("columns_accepted", &UnionCensusClass::columns_accepted);
        row("columns_rejected", &UnionCensusClass::columns_rejected);
        row("columns_repeated_accepted", &UnionCensusClass::columns_repeated_accepted);
        row("columns_repeated_rejected", &UnionCensusClass::columns_repeated_rejected);
        row("prefilter_accepted_ns", &UnionCensusClass::prefilter_accepted_ns);
        row("prefilter_rejected_ns", &UnionCensusClass::prefilter_rejected_ns);
        row("identities_first_here", &UnionCensusClass::identities_first_here);
        text << "read_length\t" << census.read_length << "\t\nwindow_length\t" << census.window_length
             << "\t\nread_forms\t" << census.read_forms << "\t\nread_contents\t" << census.read_contents
             << "\t\nwindow_positions\t" << census.window_positions << "\t\nwindow_contents\t" << census.window_contents
             << "\t\nwindow_table_bytes\t" << census.window_table_bytes << "\t\nwindow_index_bytes\t"
             << census.window_index_bytes << "\t\nidentities\t" << census.identities << "\t\nrecord_bytes\t"
             << census.record_bytes << "\t\nfill_executed_prefilters\t" << census.fill_executed_prefilters
             << "\t\nfill_verifications\t" << census.fill_verifications << "\t\ndecisions_checked\t"
             << census.decisions_checked << "\t\nseconds_windows\t" << census.seconds_windows << "\t\nseconds_stream\t"
             << census.seconds_stream << "\t\nseconds_aggregate\t" << census.seconds_aggregate << "\t\n";
        for (const auto& [occurrences, identities] : census.multiplicity) {
            text << "multiplicity\t" << occurrences << '\t' << identities << '\n';
        }
        for (std::size_t i = 0; i < census.accumulation_panel.size(); ++i) {
            text << "accumulation_panel\t" << i + 1 << '\t' << census.accumulation_panel[i] << '\n';
        }
        for (std::size_t i = 0; i < census.accumulation_mosaic.size(); ++i) {
            text << "accumulation_mosaic\t" << i + 1 << '\t' << census.accumulation_mosaic[i] << '\n';
        }
        atomic_write(options.out_prefix + ".union_verification_census.tsv", text.str());
    }
    if (options.diagnostic_union_seeding_profile > 0 && any_ranked) {
        // DIAGNOSTIC ONLY: a deterministic replay of the seeding path over evenly spaced panel and mosaic haplotypes.
        std::set<std::vector<std::size_t>> panel_paths(model.panel().template_alleles.begin(),
                                                       model.panel().template_alleles.end());
        std::vector<char> is_panel(interned.sequences.size(), 0);
        for (const auto& [alleles, slot] : slot_of_path) {
            if (panel_paths.count(alleles) != 0) is_panel[interned.index_of[slot]] = 1;
        }
        std::vector<std::size_t> of_class[2];
        for (std::size_t h = 0; h < interned.sequences.size(); ++h) of_class[is_panel[h] != 0 ? 0 : 1].push_back(h);
        std::vector<std::size_t> sample;
        for (const std::vector<std::size_t>& members : of_class) {
            const std::size_t n = std::min(options.diagnostic_union_seeding_profile, members.size());
            for (std::size_t i = 0; i < n; ++i) sample.push_back(members[i * members.size() / n]);
        }
        const UnionSeedingProfile profile =
            profile_union_seeding(interned, fragments, library, scoring_options, sample, is_panel);
        std::ostringstream text;
        text << "# DIAGNOSTIC ONLY: single-thread replay of the union seeding path; nanoseconds over both strands\n"
             << "key\tpanel\tmosaic\n";
        const auto row = [&](const char* key, std::uint64_t UnionSeedingPhases::*field) {
            text << key << '\t' << profile.panel.*field << '\t' << profile.mosaic.*field << '\n';
        };
        row("haplotypes", &UnionSeedingPhases::haplotypes);
        row("production_ns", &UnionSeedingPhases::production_ns);
        row("scan_ns", &UnionSeedingPhases::scan_ns);
        row("scan_and_generate_ns", &UnionSeedingPhases::scan_and_generate_ns);
        row("start_sort_ns", &UnionSeedingPhases::start_sort_ns);
        row("start_sort_bucketed_ns", &UnionSeedingPhases::start_sort_bucketed_ns);
        row("window_check_ns", &UnionSeedingPhases::window_check_ns);
        row("prefilter_ns", &UnionSeedingPhases::prefilter_ns);
        row("prefilter_early_ns", &UnionSeedingPhases::prefilter_early_ns);
        row("emission_call_ns", &UnionSeedingPhases::emission_call_ns);
        row("dp_timed_ns", &UnionSeedingPhases::dp_timed_ns);
        row("storage_ns", &UnionSeedingPhases::storage_ns);
        row("replay_ns", &UnionSeedingPhases::replay_ns);
        row("row_sort_ns", &UnionSeedingPhases::row_sort_ns);
        row("piece_hits", &UnionSeedingPhases::piece_hits);
        row("proposals", &UnionSeedingPhases::proposals);
        row("unique_starts", &UnionSeedingPhases::unique_starts);
        row("prefilters", &UnionSeedingPhases::prefilters);
        row("prefilter_rejections", &UnionSeedingPhases::prefilter_rejections);
        row("dp_requests", &UnionSeedingPhases::dp_requests);
        row("executed_dps", &UnionSeedingPhases::executed_dps);
        row("dp_timed", &UnionSeedingPhases::dp_timed);
        row("placements", &UnionSeedingPhases::placements);
        row("replayed", &UnionSeedingPhases::replayed);
        for (const auto& [key, field] : std::vector<std::pair<const char*, std::uint64_t UnionSeedingPhases::*>>{
                 {"survivors", &UnionSeedingPhases::survivors}, {"key_full_ns", &UnionSeedingPhases::key_full_ns},
                 {"key_prefix_ns", &UnionSeedingPhases::key_prefix_ns}, {"lookup_ns", &UnionSeedingPhases::lookup_ns},
                 {"lookup_unclocked_ns", &UnionSeedingPhases::lookup_unclocked_ns},
                 {"lookup_nolock_ns", &UnionSeedingPhases::lookup_nolock_ns},
                 {"lookup_nocopy_ns", &UnionSeedingPhases::lookup_nocopy_ns},
                 {"probe_only_ns", &UnionSeedingPhases::probe_only_ns},
                 {"clock_pairs_ns", &UnionSeedingPhases::clock_pairs_ns},
                 {"batched_starts", &UnionSeedingPhases::batched_starts}, {"lanes2_ns", &UnionSeedingPhases::lanes2_ns},
                 {"lanes4_ns", &UnionSeedingPhases::lanes4_ns}, {"lanes8_ns", &UnionSeedingPhases::lanes8_ns}}) {
            row(key, field);
        }
        text << "bucketed_sort_lists_checked\t" << profile.bucketed_sort_checked << "\t\n"
             << "early_prefilter_decisions_checked\t" << profile.early_prefilter_checked << "\t\n"
             << "lane_decisions_checked\t" << profile.lane_decisions_checked << "\t\n"
             << "prefix_keys_checked\t" << profile.prefix_keys_checked << "\t\n"
             << "fill_workers\t" << profile.fill_workers << "\t\n"
             << "fill_seconds\t" << profile.fill_seconds << "\t\n";
        for (std::size_t i = 0; i < profile.fill_segments.size(); ++i) {
            text << "fill_segment\t" << profile.fill_segments[i].first << '\t' << profile.fill_segments[i].second << '\n';
        }
        for (std::size_t i = 0; i < profile.sampled.size(); ++i) {
            text << "sampled_haplotype\t" << profile.sampled[i] << '\t' << (is_panel[profile.sampled[i]] ? "panel" : "mosaic")
                 << ':' << std::hex << profile.sampled_hashes[i] << std::dec << '\n';
        }
        atomic_write(options.out_prefix + ".union_seeding_profile.tsv", text.str());
    }

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
    bool evidence_dump_pending = true;   // S2-B (DIAGNOSTIC ONLY): the union evidence dump is written once
    const auto rank_set = [&](const MosaicBeamResult& beam, std::size_t view = 0) {
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
            scores.push_back((view == 0 ? scored.scores : scored.audit_scores.at(view - 1))[scored_index.at({candidate.first, candidate.second})]);

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
        // S2-B (DIAGNOSTIC ONLY): the evaluator for the pre-placement bound prototype, once, from the primary ranking
        if (view == 0 && options.diagnostic_union_evidence_dump && evidence_dump_pending &&
            scored.log_mass.size() == fragments.size() * interned.sequences.size()) {
            evidence_dump_pending = false;
            std::string bin;
            const auto put = [&](const auto& v) { bin.append(reinterpret_cast<const char*>(&v), sizeof v); };
            bin.append("PVEV1\0\0\0", 8);
            put(static_cast<std::uint64_t>(interned.sequences.size()));
            put(static_cast<std::uint64_t>(fragments.size()));
            put(static_cast<std::uint64_t>(sequence_candidates.size()));
            put(library.insert_mean_bp);
            put(library.insert_sd_bp);
            put(library.support_sigmas);
            put(library.substitution_error);
            put(library.outlier_weight);
            put(library.fragment_start_intensity);
            put(library.gap_open);
            put(library.gap_extend);
            put(static_cast<std::uint64_t>(library.max_gap_bp));
            put(static_cast<std::uint64_t>(scoring_options.mismatch_band));
            for (const std::string& h : interned.sequences) {
                put(static_cast<std::uint64_t>(h.size()));
                bin.append(h);
            }
            for (const double x : scored.log_mass) put(x);
            for (const double x : scored.exposure) put(x);
            for (const ReadPair& f : fragments) put(log_background_probability(f));   // the scorer's own definition
            for (const double x : haploid_sequence_prior) put(x);
            for (std::size_t i = 0; i < sequence_candidates.size(); ++i) {
                put(static_cast<std::uint64_t>(sequence_candidates[i].first));
                put(static_cast<std::uint64_t>(sequence_candidates[i].second));
                put(sequence_prior[i]);
                put(scores[i].log_likelihood);
                put(scores[i].fragment_term);
                put(scores[i].exposure_term);
                put(scores[i].absolute_term_sum);
                put(static_cast<std::uint64_t>(scores[i].fragments_with_mass));
            }
            put(static_cast<std::uint64_t>(ranking.winner_set.size()));
            for (const std::size_t w : ranking.winner_set) put(static_cast<std::uint64_t>(w));
            put(static_cast<std::uint64_t>(ranking.candidates.front().candidate_id));
            std::uint64_t early_first = std::numeric_limits<std::uint64_t>::max(), early_second = early_first;
            if (!g_incumbent.early_first.empty()) {
                std::size_t x = interned.index_of[slot_of_path.at(g_incumbent.early_first)];
                std::size_t y = interned.index_of[slot_of_path.at(g_incumbent.early_second)];
                if (y < x) std::swap(x, y);
                early_first = x;
                early_second = y;
            }
            put(early_first);
            put(early_second);
            atomic_write(options.out_prefix + ".union_evidence.bin", bin);
        }

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
            if (options.block_support_rule != "legacy") {   // M3-E (EXPERIMENTAL): the same shared function, unpooled inputs
                const std::size_t blocks_here = model.blocks().size();
                std::vector<SupportGroup> groups(sequence_candidates.size());
                for (std::size_t g = 0; g < sequence_candidates.size(); ++g) {
                    std::set<BlockGenotypeRow> rows;
                    for (const std::size_t m : sequence_candidates[g].members) {
                        BlockGenotypeRow row;
                        for (std::size_t b = 0; b < blocks_here; ++b) {
                            row.push_back(std::minmax(candidates[m].diplotype.first.alleles[b], candidates[m].diplotype.second.alleles[b]));
                        }
                        rows.insert(row);
                    }
                    groups[g].rows.assign(rows.begin(), rows.end());
                }
                summary.data_support_competitor_paths.assign(sequence_candidates.size(), "NA");
                for (const CertifiedCandidate& c : ranking.candidates) {
                    SupportGroup& g = groups[c.candidate_id];
                    g.fragment = c.fragment_term;
                    g.exposure = c.exposure_term;
                    g.prior = c.prior;
                    g.bound = c.numeric_error_bound;
                    const MosaicDiplotype& rep = candidates[sequence_candidates[c.candidate_id].representative].diplotype;
                    summary.data_support_competitor_paths[c.candidate_id] = join_alleles(rep.first.alleles) + "|" + join_alleles(rep.second.alleles);
                }
                if (options.block_support_rule == "full_data_likelihood_v2_experimental") {
                    summary.data_support_v2 = block_data_support_v2(groups, best.candidate_id, blocks_here);
                } else {
                    summary.data_support = block_data_support_v1(groups, best.candidate_id, blocks_here);
                }
            }
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
        // I2-C: where the search began, was built and returned (measurement only)
        Clock::time_point search_begin = repair_mark, search_built = repair_mark, search_end = repair_mark;
        bool search_started = false, search_returned = false;
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
        // S2-K: the repair path's final outputs are prepared in memory while the resource watchdog is still armed, and
        // published together by flush_outputs: the run is claimed first (the watchdog can no longer stop it), then every
        // prepared file is written by atomic replacement, then the sequence report and the status (stop_with). A stop
        // before the claim removes every final output and keeps only what the existing contract keeps (the S2-A
        // checkpoint, as PROVISIONAL); nothing is written between a stop's cleanup and its exit, because the claim and
        // the stop exclude each other under the watch mutex.
        std::deque<std::pair<std::string, std::ostringstream>> pending_outputs;
        const auto defer_output = [&](const std::string& path) -> std::ostringstream& {
            pending_outputs.emplace_back(path, std::ostringstream{});
            return pending_outputs.back().second;
        };
        const auto flush_outputs = [&]() {
            note_footprint("outputs_prepared");
            if (options.test_interrupt_point == "post_refinement_outputs") {   // S2-K TEST HOOK: a memory stop before the claim
                stop_run("RESOURCE_REFUSAL", "test hook: post-refinement, outputs prepared and not yet published", 125);
            }
            claim_run();
            for (auto& [path, text] : pending_outputs) atomic_write(path, text.str());
            pending_outputs.clear();
            note_footprint("outputs_published");
        };
        const auto stop_with = [&](GenotypeRunStatus status, const std::string& note, const CertifiedRanking* ranking,
                                   const RankingSummary* summary) {
            timings.candidate_grouping = run_timings.candidate_grouping;
            timings.prior = run_timings.prior;
            timings.ranking = run_timings.ranking;
            set_report.distinct_haplotype_paths = candidate_report.distinct_haplotype_paths;
            set_report.haploid_prior_evaluations = candidate_report.haploid_prior_evaluations;
            set_report.diploid_prior_lookups = candidate_report.diploid_prior_lookups;
            report.seconds = seconds_since(repair_mark);
            if (search_started) {
                const auto between = [](Clock::time_point a, Clock::time_point b) {
                    return std::chrono::duration<double>(b - a).count();
                };
                const Clock::time_point stop = Clock::now();
                const Clock::time_point returned = search_returned ? search_end : stop;
                report.timed = true;
                report.seconds_pre_search = between(repair_mark, search_begin);
                report.seconds_search_construction = between(search_begin, search_built);
                report.seconds_search_run = between(search_built, returned);
                report.seconds_post_search = between(returned, stop) - report.seconds_tie_branching -
                                             report.seconds_representation_guard;
            }
            report.call_class = status == GenotypeRunStatus::Complete                        ? "strict_call"
                                : status == GenotypeRunStatus::ConditionalGenotypeEquivalent ? "genotype_equivalent_call"
                                                                                             : "no_call";
            if (report.call_class != "no_call") report.no_call_reason = "none";
            g_sequence_legacy = {genotype_run_status_name(status), report.call_class, report.no_call_reason,
                                 options.dosage_mode, dosage_mode_source(options)};   // S1-B, S1-D
            if (options.output_schema == "sequence_first_v2" && !g_dosage_stage_follows) {   // S1-D.1
                g_sequence_legacy.dosage_mode = "not_run";
                g_sequence_legacy.dosage_mode_source = "caller_only_api_without_dosage_stage";
            }
            // S1-C: the sequence report is published before the status commits it: its files first, its locus record
            // last; a failure leaves no partial report and the status says so.
            flush_outputs();   // S2-K: claim, then the prepared outputs (S2-A: the watchdog can no longer overwrite the report)
            if (sequence_first(options) && !g_sequence_published) {   // S2-A: a failed checkpoint already published
                compound_prepare(g_sequence_report, options);   // M9-J: enabled; rows only if a pool or a provisional set them
                publish_sequence_report(g_sequence_report, g_sequence_legacy, options.out_prefix,
                                        "no completed ranking and no safely available evaluated candidate (" +
                                            g_sequence_legacy.run_status + ", " + g_sequence_legacy.no_call_reason + ")");
                g_sequence_published = true;
                note_footprint("sequence_report_published");
            }
            timings.total = seconds_since(started);
            write_status(status_path, status, fingerprint, options, library, timings, &locus_report, &panel, beam_view,
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

        // ---- S2-A: the union pass completed: the early candidate must reproduce there bit for bit, and the union's
        // best fully evaluated candidate replaces the incumbent (a complete, exactly scored set)
        if (g_incumbent.enabled && !g_incumbent.failed) {
            const auto sequence_key = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
                std::size_t x = interned.index_of[slot_of_path.at(a)], y = interned.index_of[slot_of_path.at(b)];
                if (y < x) std::swap(x, y);
                return std::make_pair(x, y);
            };
            const auto union_prior = [&](const std::pair<std::size_t, std::size_t>& key) {
                return haploid_sequence_prior[key.first] + haploid_sequence_prior[key.second] +
                       (key.first != key.second ? std::log(2.0) : 0.0);
            };
            if (!g_incumbent.early_first.empty() && g_incumbent.union_check == "not_reached") {
                // S2-D: against the stored early candidate (the incumbent may already be a later checkpoint)
                const auto key = sequence_key(g_incumbent.early_first, g_incumbent.early_second);
                const DiplotypeScore& u = scored.scores[scored_index.at(key)];
                const DiplotypeScore& e = g_incumbent.early_score;
                const auto bits = [](double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; };
                const bool same = bits(u.log_likelihood, e.log_likelihood) && bits(u.fragment_term, e.fragment_term) &&
                                  bits(u.exposure_term, e.exposure_term) &&
                                  bits(u.absolute_term_sum, e.absolute_term_sum) &&
                                  u.fragments_with_mass == e.fragments_with_mass &&
                                  bits(union_prior(key), g_incumbent.early_prior);
                g_incumbent.union_check = same ? "reproduced" : "MISMATCH";
                if (!same) {
                    g_incumbent.failed = true;
                    g_incumbent.failure = "scoring_inconsistency: the union pass did not reproduce the early incumbent";
                    report.no_call_reason = "early_incumbent_not_reproduced";
                    std::ostringstream why;
                    why << std::setprecision(17) << "scoring inconsistency: the early incumbent scored " << e.log_likelihood
                        << " and the union pass " << u.log_likelihood << "; an explicit failure, never a fallback";
                    return stop_with(GenotypeRunStatus::RescoreIncomplete, why.str(), nullptr, nullptr);
                }
            }
            const CertifiedCandidate& top = union_ranking.candidates.front();
            const SequenceCandidate& group = chosen.sequence_candidates[top.candidate_id];
            const MosaicDiplotype& best = candidates[group.representative].diplotype;
            const std::pair<std::size_t, std::size_t> key{group.first, group.second};
            const DiplotypeScore& score = scored.scores[scored_index.at(key)];
            const double prior = union_prior(key);
            g_incumbent.evaluated = candidates.size();
            if (!g_incumbent.have || score.log_likelihood + prior >= g_incumbent.total) {
                g_incumbent.have = true;
                g_incumbent.first = best.first.alleles;
                g_incumbent.second = best.second.alleles;
                g_incumbent.score = score;
                g_incumbent.sequence_prior = prior;
                g_incumbent.total = score.log_likelihood + prior;
                g_incumbent.stage = "union_best";
                g_incumbent.tied_groups = std::to_string(union_ranking.winner_set.size());
                g_incumbent.projection_ambiguous = group.projection_ambiguous ? "yes" : "no";
                g_incumbent.explored_labellings = std::to_string(group.members.size());
            }
            if (!commit_checkpoint() && g_incumbent.failed) {
                report.no_call_reason = "checkpoint_publication_failed";
                return stop_with(GenotypeRunStatus::RescoreIncomplete,
                                 "the union checkpoint could not be published (" + g_incumbent.failure +
                                     "); an explicit failure, never a fallback",
                                 nullptr, nullptr);
            }
        }

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
            if (g_incumbent.enabled && g_incumbent.have && !g_incumbent.failed) {   // S2-A: a declared cap stop
                g_sequence_report = checkpoint_report("neighbourhood_bound", "the repair neighbourhood exceeds its declared "
                                                      "bound; the best fully evaluated candidate so far; search incomplete");
            }
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
        // library_sensitivity_audit_v1: every pooled identity's score at every view (the run's library first), from the
        // same combination passes that scored it. A second scoring of one identity must agree at every view, bit for bit.
        const std::size_t view_count = 1 + scoring_options.audit_libraries.size();
        const bool audit_on = view_count > 1;
        std::map<DiplotypeIdentity, std::vector<DiplotypeScore>> pool_views;
        std::uint64_t view_mismatches = 0;
        const auto same_score = [&](const DiplotypeScore& x, const DiplotypeScore& y) {
            return bitwise_equal(x.log_likelihood, y.log_likelihood) && bitwise_equal(x.fragment_term, y.fragment_term) &&
                   bitwise_equal(x.exposure_term, y.exposure_term) &&
                   bitwise_equal(x.absolute_term_sum, y.absolute_term_sum) &&
                   x.fragments_with_mass == y.fragments_with_mass;
        };
        const auto note_views = [&](const DiplotypeIdentity& key, const std::vector<DiplotypeScore>& views) {
            const auto placed = pool_views.emplace(key, views);
            if (placed.second) return;
            for (std::size_t v = 0; v < view_count; ++v) {
                if (placed.first->second.size() != view_count || views.size() != view_count ||
                    !same_score(placed.first->second[v], views[v])) {
                    ++view_mismatches;
                    break;
                }
            }
        };
        const auto union_views = [&](std::size_t index) {
            std::vector<DiplotypeScore> views{scored.scores[index]};
            for (const std::vector<DiplotypeScore>& audit : scored.audit_scores) views.push_back(audit[index]);
            return views;
        };
        // The candidate the last pool_scored call recorded, so the views hook that follows it can attach its views.
        DiplotypeIdentity last_scored;
        std::vector<TieItem>* last_items = nullptr;
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
                if (audit_on) {
                    note_views(diplotype_identity(c.diplotype.first.alleles, c.diplotype.second.alleles),
                               union_views(scored_index.at(key)));
                }
            }
        }
        note_footprint("candidate_pool_constructed");   // S2-G
        // The union pass's mass matrix is not needed again; release it before the repair context exists.
        std::vector<double>().swap(scored.log_mass);
        std::vector<OriginSummary>().swap(scored.origin_summary);
        // What the process carries into the repair. The union pass leaves allocator-retained pages behind,
        // and how many varies run to run: NA19240 30x entered the repair at 1.55 GB in one run and 1.15 GB
        // in another, with the same scores. Recorded so a footprint refusal can be attributed.
        report.footprint_at_repair_start = physical_footprint().current;

        // ---- 3. the exact coordinate repair from the anchor
        RepairSetup repair_setup = make_repair_setup(options, scoring_options);   // S2-H: shared with the replay
        FragmentScoringOptions& repair_scoring = repair_setup.scoring;
        if (options.adaptive_emission_handoff) repair_scoring.emission_handoff_in = &union_emissions;   // I2-K
        RepairSearchOptions& search_options = repair_setup.search;
        if (options.repair_junction_proposals == "experimental_v1") {   // M9-G: one sample-scoped read table per run
            const Clock::time_point table_mark = Clock::now();
            auto table = std::make_shared<const JunctionReadKmers>(JunctionReadKmers::from_fragments(fragments));
            report.junction_table_seconds = seconds_since(table_mark);
            report.junction_read_kmers = table->codes.size();
            report.junction_read_bytes = table->bytes();
            search_options.junction_reads = std::move(table);
            search_options.junction_rate = report.junction_rate = depth.rate_per_copy;
            search_options.junction_background = report.junction_background = depth.background;
        }
        if (options.diagnostic_refinement_memory) {   // S2-H: what a fresh-process replay of round 1 must match
            std::ostringstream t;
            t << "key\tvalue\nschema\tdiagnostic_refinement_replay_inputs_v1\nnote\tDIAGNOSTIC ONLY: replay inputs, not a result\n"
              << "model_fingerprint\t" << fingerprint << '\n';
            for (const auto& hash : g_input_hashes) t << "sha256_" << hash.first << '\t' << hash.second << '\n';
            t << "emissions_digest\t" << emissions_digest(emissions) << "\nrepair_setup\t" << repair_setup_text(repair_setup, fragments.size())
              << "\nanchor_hap1\t" << join_alleles(anchor_diplotype.first.alleles) << "\nanchor_hap2\t"
              << join_alleles(anchor_diplotype.second.alleles) << "\nanchor_score\t" << hex_double(anchor_score)
              << "\nanchor_union_rank\t" << anchor_rank + 1 << '\n';
            atomic_write(options.out_prefix + ".refinement_replay_inputs.tsv", t.str());
        }
        write_running_status("repair");
        std::unique_ptr<CoordinateRepairSearch> search;
        std::string refusal;
        bool reproduced = true;
        // One scored candidate into the pool. The diagnostic also keeps the round's items, so a terminal
        // tied round can be reconstructed from the scores the search already computed.
        using RoundItem = TieItem;
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
            if (items != nullptr) items->push_back({round, a, b, da, db, score, prior, {}});
            last_scored = diplotype_identity(a, b);
            last_items = items;
        };
        // The views hook: the candidate pool_scored just recorded, at every view.
        const auto record_views = [&](const std::vector<DiplotypeScore>& views) {
            note_views(last_scored, views);
            if (last_items != nullptr) last_items->back().views = views;
        };
        // Central tie branches, by start and history, as each ran: what an endpoint that decides identically reuses.
        struct CentralBranch {
            BranchResult result;
            std::vector<RepairRound> rounds;       // its rounds, each with its verdict at every view
            std::vector<double> terminal_totals;   // per view
        };
        std::map<std::pair<DiplotypeIdentity, std::vector<DiplotypeIdentity>>, CentralBranch> central_branches;
        try {
            search_begin = Clock::now();
            search_started = true;
            search = std::make_unique<CoordinateRepairSearch>(model, emissions, fragments, library, repair_scoring,
                                                              search_options);
            note_footprint("refinement_context_constructed");   // S2-G
            search->start(anchor_diplotype.first.alleles, anchor_diplotype.second.alleles, anchor_score);
            note_footprint("repair_search_started");   // S2-H: start() records the anchor; placement happens in run()
            search_built = Clock::now();
            RepairHooks hooks;
            hooks.on_checkpoint = [&](const char* point) {   // S2-H: the seed's actual placement and census
                note_footprint(point);
                if (options.test_interrupt_point == "refinement_seed_census" &&
                    std::string(point) == "initial_incumbent_census_retained") {   // S2-H TEST HOOK: a stop inside round 1
                    stop_run("TIME_LIMIT", "test hook: inside the first refinement round, after the seed census", 124);
                }
            };
            hooks.on_scored = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b,
                                  const std::string& sa, const std::string& sb, const DiplotypeScore& score,
                                  double prior) {
                ++report.pool_repair_scored;
                const std::size_t round = search->rounds().size() + 1;
                pool_scored(a, b, sa, sb, score, prior, "repair", "round" + std::to_string(round),
                            options.diagnostic_tie_branching || options.repair_tie_branching || audit_on
                                ? &production_items
                                : nullptr,
                            round);
                if (g_incumbent.enabled) {   // S2-A: a fully evaluated candidate; replaces the incumbent only if better
                    ++g_incumbent.evaluated;
                    if (score.log_likelihood + prior > g_incumbent.total) {
                        g_incumbent.have = true;
                        g_incumbent.first = a;
                        g_incumbent.second = b;
                        g_incumbent.score = score;
                        g_incumbent.sequence_prior = prior;
                        g_incumbent.total = score.log_likelihood + prior;
                        g_incumbent.stage = "repair_pool_best";
                        g_incumbent.tied_groups = g_incumbent.projection_ambiguous = g_incumbent.explored_labellings = "NA";
                        g_incumbent.dirty = true;
                    }
                }
            };
            if (audit_on) hooks.on_scored_views = record_views;
            hooks.on_round = [&](const RepairRound& round) {
                commit_partial(prefix + ".repair_rounds.partial.tsv", search->render_rounds(search->rounds()));
                if (search->memory_on()) commit_partial(prefix + ".repair_memory.partial.tsv", search->memory_text());
                if (g_watch != nullptr) {
                    std::lock_guard<std::mutex> lock(g_watch->mutex);
                    g_watch->repair_rounds_committed = search->rounds().size();
                }
                write_running_status("repair round " + std::to_string(round.index) + " committed");
                if (search->rounds().size() == 1) note_footprint("first_committed_refinement_round");   // S2-G
                if (g_incumbent.enabled && g_incumbent.dirty && !g_incumbent.failed) (void)commit_checkpoint();   // S2-A
                if (options.test_interrupt_point == "adaptive_refinement") {   // S2-A TEST HOOK
                    stop_run("TIME_LIMIT", "test hook: during adaptive refinement (after the first committed round)", 124);
                }
                if (options.test_interrupt_after_repair_rounds != 0 && round.moved &&
                    search->rounds().size() >= options.test_interrupt_after_repair_rounds) {
                    stop_run("TIME_LIMIT",
                             "test hook: interrupted after " + std::to_string(search->rounds().size()) +
                                 " committed repair rounds (the time-limit path, triggered deterministically)",
                             124);
                }
            };
            reproduced = search->run(hooks);
            accumulate_junction(report.junction, search->junction_report());   // M9-G
            search_end = Clock::now();
            search_returned = true;
            report.stages = search->stage_seconds();
            if (search_options.round1_reuse_only_admission && search->context() != nullptr) {   // S2-J
                report.admission_reported = true;
                report.admission = search->context()->origin_admission_report();
            }
        } catch (const ResourceExhausted& error) {
            search_end = Clock::now();
            search_returned = true;
            if (search) report.stages = search->stage_seconds();
            if (search && search_options.round1_reuse_only_admission && search->context() != nullptr) {   // S2-J
                report.admission_reported = true;
                report.admission = search->context()->origin_admission_report();
            }
            refusal = error.what();
        } catch (const std::exception& error) {
            report.no_call_reason = "repair_failed";
            return stop_with(GenotypeRunStatus::RepairIncomplete,
                             std::string("the repair search failed and returned no genotype: ") + error.what(),
                             nullptr, nullptr);
        }
        if (g_incumbent.enabled && g_incumbent.failed) {   // S2-A: a checkpoint publication failed during the search
            report.no_call_reason = "checkpoint_publication_failed";
            return stop_with(GenotypeRunStatus::RepairIncomplete,
                             "a repair checkpoint could not be published (" + g_incumbent.failure +
                                 "); an explicit failure, never a fallback",
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
                report.fine_window = search->context()->fine_window_report();
                report.full_placements = search->context()->full_placement_count();
                report.enumeration_workers = search->context()->parallel_enumeration_stats().workers;
            }
            report.seconds_scoring = search->cumulative_seconds();
        }
        report.pool_identities = pool.size();
        for (const PoolEntry& entry : pool) report.pool_repair_only += entry.in_union ? 0 : 1;
        report.pool_score_mismatches = mismatches;
        // S1-B (sequence_first_v1 only): a declared-cap stop keeps the last committed, exactly scored incumbent as
        // PROVISIONAL -- only after at least one committed round (the anchor reproduced and was rescored).
        const auto provisional = [&](const std::string& stop, const std::string& note) {
            if (g_incumbent.enabled) {   // S2-A: the best fully evaluated candidate so far (never the partial search state)
                if (g_incumbent.have && !g_incumbent.failed) g_sequence_report = checkpoint_report(stop, note);
                return;
            }
            if (!sequence_first(options) || !search || search->rounds().empty()) return;
            std::vector<std::string> block_ids;
            for (const auto& block : model.blocks()) block_ids.push_back(block.block_id);
            g_sequence_report = provisional_sequence_report(
                search->incumbent_first(), search->incumbent_second(), search->incumbent_total(), stop, note, block_ids,
                [&](const std::vector<std::size_t>& labels) { return spell_mosaic(model, MosaicPath{labels}).bases; });
            if (options.output_schema == "sequence_first_v2") annotate_v2(g_sequence_report, model, fingerprint);
            compound_unassessed(g_sequence_report, model, fingerprint, options, "NOT_ASSESSED_PROVISIONAL");   // M9-J
        };
        if (!refusal.empty()) {
            report.no_call_reason = "cumulative_cap";
            provisional("cumulative_cap", "the repair search stopped at a declared cumulative cap (" + refusal +
                                              "); the last committed incumbent, exactly scored; search incomplete");
            // Committed rounds stay as <prefix>.repair_rounds.partial.tsv: a record, never a call.
            return stop_with(GenotypeRunStatus::RescoreIncomplete,
                             "the repair search stopped at a declared cumulative cap and returned no genotype: " +
                                 refusal,
                             nullptr, nullptr);
        }
        // Every stop that returned normally leaves its whole trajectory under the final name. M3-D: a joint continuation
        // leg's rounds follow the trunk's, numbered on from it (its header line dropped).
        std::string joint_leg_rounds_text;
        const auto write_trajectory = [&]() {
            atomic_write(prefix + ".repair_rounds.tsv", search->render_rounds(search->rounds()) + joint_leg_rounds_text);
            if (search->memory_on()) atomic_write(prefix + ".repair_memory.tsv", search->memory_text());
            for (const char* suffix : kPartialOutputs) std::remove((prefix + suffix).c_str());
        };

        // ---- 3a. M3-D (EXPERIMENTAL, repair_joint_search experimental_v1): ONE bounded joint-block pass after a uniquely
        // resolved, converged trunk, in the refinement's own context (never a second live context; caps cumulative and
        // unchanged). Its complete catalogue joins the pool; an improving move resumes ordinary single-block refinement for
        // the remaining declared rounds (no second joint pass). Unresolved or tied trunks, a missing delta context and the
        // model-derived bounds skip it, recorded by name; a cap crossed inside it keeps the existing incomplete contract.
        if (options.repair_joint_search == "experimental_v1") {
            const Clock::time_point joint_mark = Clock::now();
            PreparedScoringContext* shared = search->context();
            JointSearchOptions joint_options;
            joint_options.switch_probability = options.switch_probability;
            joint_options.request_order_seed = search_options.request_order_seed;   // gate hook only (0 in production)
            if (options.test_joint_search_single_limit != 0) joint_options.single_limit = options.test_joint_search_single_limit;
            report.joint.bounds = joint_search_bounds(model, joint_options);
            const bool last_unique = !search->rounds().empty() && search->rounds().back().unique_best;
            const auto note_incumbent = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b,
                                            const DiplotypeScore& score, double prior, const char* stage) {
                if (!g_incumbent.enabled) return;   // S2-A: a fully evaluated candidate replaces the incumbent only if better
                ++g_incumbent.evaluated;
                if (score.log_likelihood + prior > g_incumbent.total) {
                    g_incumbent.have = true;
                    g_incumbent.first = a;
                    g_incumbent.second = b;
                    g_incumbent.score = score;
                    g_incumbent.sequence_prior = prior;
                    g_incumbent.total = score.log_likelihood + prior;
                    g_incumbent.stage = stage;
                    g_incumbent.tied_groups = g_incumbent.projection_ambiguous = g_incumbent.explored_labellings = "NA";
                    g_incumbent.dirty = true;
                }
            };
            if (!reproduced) {
                report.joint_search_status = "skipped_trunk_not_reproduced";
            } else if (!report.converged) {
                report.joint_search_status = "skipped_trunk_unresolved_" + report.stop_reason;
            } else if (!last_unique) {
                report.joint_search_status = "skipped_trunk_tied_at_incumbent";
            } else if (shared == nullptr || search_options.mode != RepairScoringMode::Delta) {
                report.joint_search_status = "skipped_no_delta_context";
            } else if (!report.joint.bounds.within) {
                report.joint_search_status = "skipped_bounds";
                report.joint.status = "skipped_bounds";
                report.joint.decision = "skipped";
            } else {
                write_running_status("joint search");
                const FragmentScoringCounters before = shared->cumulative_counters();
                try {
                    report.joint = run_joint_block_search(
                        *shared, model, search->incumbent_first(), search->incumbent_second(), search->incumbent_slot_first(),
                        search->incumbent_slot_second(), search->incumbent_total(), fragments.size(), repair_scoring.mismatch_band,
                        joint_options,
                        [&](const JointCandidate& c, const std::string& sa, const std::string& sb) {
                            ++report.pool_joint_scored;
                            pool_scored(c.first, c.second, sa, sb, c.score, c.prior, "joint", "joint_pass", nullptr, 0);
                            note_incumbent(c.first, c.second, c.score, c.prior, "joint_pool_best");
                        },
                        {});
                } catch (const ResourceExhausted& error) {
                    refusal = error.what();
                }
                report.joint_physical_origins = shared->cumulative_counters().physical_origins - before.physical_origins;
                report.joint_state_classes = shared->cumulative_counters().statistical_state_classes - before.statistical_state_classes;
                report.joint_footprint_peak = physical_footprint().peak;
                if (!refusal.empty()) {
                    report.joint_search_status = "refused_cumulative_cap";
                    report.no_call_reason = "cumulative_cap";
                    provisional("cumulative_cap", "the joint-block pass stopped at a declared cumulative cap (" + refusal +
                                                      "); the last committed incumbent, exactly scored; search incomplete");
                    return stop_with(GenotypeRunStatus::RescoreIncomplete,
                                     "the joint-block pass stopped at a declared cumulative cap and returned no genotype: " + refusal,
                                     nullptr, nullptr);
                }
                report.joint_search_status = report.joint.status + "_" + report.joint.decision;
                if (report.joint.moved) {
                    const JointCandidate& w = report.joint.catalogue[report.joint.winner];
                    const std::size_t first_round = search->rounds().back().index + 1;
                    report.joint_leg_first_round = first_round;
                    report.moves += 1;   // the joint move itself
                    report.final_incumbent_score = w.score.log_likelihood + w.prior;
                    report.final_incumbent_hap1 = join_alleles(w.first);
                    report.final_incumbent_hap2 = join_alleles(w.second);
                    if (first_round > options.repair_rounds) {
                        // No declared round remains to confirm the moved incumbent: it is not a converged optimum.
                        report.joint_leg_stop = "no_rounds_remaining";
                        report.stop_reason = "round_cap";
                        report.converged = false;
                    } else {
                        const Clock::time_point leg_mark = Clock::now();
                        RepairSearchOptions leg_options = search_options;
                        leg_options.external_context = true;
                        leg_options.memory_log = false;
                        leg_options.memory_live_path.clear();
                        leg_options.round1_reuse_only_admission = false;
                        leg_options.first_round = first_round;
                        leg_options.max_rounds = options.repair_rounds + 1 - first_round;
                        CoordinateRepairSearch leg(model, emissions, fragments, library, repair_scoring, leg_options);
                        leg.attach_context(shared);
                        leg.start_from_slots(w.first, w.second, w.score.log_likelihood + w.prior, w.slot_first, w.slot_second);
                        leg.add_visited(anchor_diplotype.first.alleles, anchor_diplotype.second.alleles);
                        for (const RepairRound& round : search->rounds()) {
                            if (round.moved) leg.add_visited(round.first, round.second);
                        }
                        RepairHooks leg_hooks;
                        leg_hooks.on_scored = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b,
                                                  const std::string& sa, const std::string& sb, const DiplotypeScore& score,
                                                  double prior) {
                            ++report.pool_joint_leg_scored;
                            pool_scored(a, b, sa, sb, score, prior, "joint_leg",
                                        "round" + std::to_string(first_round + leg.rounds().size()), nullptr, 0);
                            note_incumbent(a, b, score, prior, "joint_leg_pool_best");
                        };
                        leg_hooks.on_round = [&](const RepairRound& round) {
                            write_running_status("joint continuation round " + std::to_string(round.index) + " committed");
                        };
                        const FragmentScoringCounters leg_before = shared->cumulative_counters();
                        bool leg_reproduced = false;
                        try {
                            leg_reproduced = leg.run(leg_hooks);
                        } catch (const ResourceExhausted& error) {
                            refusal = error.what();
                        }
                        accumulate_junction(report.junction, leg.junction_report());   // M9-G
                        report.joint_continued = true;
                        report.joint_leg_physical_origins = shared->cumulative_counters().physical_origins - leg_before.physical_origins;
                        report.joint_leg_rounds = leg.rounds().size();
                        for (const RepairRound& round : leg.rounds()) report.joint_leg_moves += round.moved ? 1 : 0;
                        {
                            const std::string text = leg.render_rounds(leg.rounds());
                            const std::size_t newline = text.find('\n');
                            joint_leg_rounds_text = newline == std::string::npos ? std::string() : text.substr(newline + 1);
                        }
                        report.joint_footprint_peak = physical_footprint().peak;
                        report.seconds_joint_leg = seconds_since(leg_mark);
                        if (!refusal.empty()) {
                            report.joint_search_status += "_leg_refused_cumulative_cap";
                            report.no_call_reason = "cumulative_cap";
                            provisional("cumulative_cap", "the joint continuation stopped at a declared cumulative cap (" + refusal +
                                                              "); the last committed incumbent, exactly scored; search incomplete");
                            return stop_with(GenotypeRunStatus::RescoreIncomplete,
                                             "the joint continuation stopped at a declared cumulative cap and returned no genotype: " +
                                                 refusal,
                                             nullptr, nullptr);
                        }
                        if (!leg_reproduced) {
                            report.no_call_reason = "joint_leg_start_not_reproduced";
                            write_trajectory();
                            return stop_with(GenotypeRunStatus::RepairIncomplete,
                                             "the joint continuation's start did not reproduce its score; no genotype", nullptr,
                                             nullptr);
                        }
                        report.joint_leg_stop = leg.stop_reason();
                        report.joint_leg_converged = leg.converged();
                        report.stop_reason = leg.stop_reason();
                        report.converged = leg.converged();
                        report.rounds_run += leg.rounds().size();
                        report.moves += report.joint_leg_moves;
                        report.final_incumbent_score = leg.incumbent_total();
                        report.final_incumbent_hap1 = join_alleles(leg.incumbent_first());
                        report.final_incumbent_hap2 = join_alleles(leg.incumbent_second());
                        shared->release_transient({});
                    }
                }
            }
            report.seconds_joint_search = seconds_since(joint_mark);
            if (shared != nullptr) {   // the run's cumulative counters include the pass and its leg
                report.cumulative = shared->cumulative_counters();
                report.origin_cache_bytes_end = shared->origin_cache_bytes();
                report.full_placements = shared->full_placement_count();
            }
        }

        // ---- 3b. v2: an improving-move tie is branched and pooled, before the run is claimed (the watchdog
        // still bounds it). Branches share the production search's scoring context and its warm origin and
        // emission caches.
        //
        // Census-free frontier: a waiting branch is its compact descriptor only (allele paths, declared score
        // components and prior, sequence digests, parent tie, history identities) and holds nothing in the
        // context. Just before a branch runs, its two haplotypes are spelled, placed by the ordinary exact path and
        // censused, and its declared score is reproduced bit for bit before any neighbour is scored; a failure
        // fails the branching closed. At every settle point -- the root tie recorded, and each branch fully
        // processed -- every repair slot is released, so the context holds scoring state for the running branch
        // alone and the frontier grows only in descriptors. (The eager frontier, a census-holding seed per waiting
        // branch built by block delta at its tie, survives as a gate hook.)
        // Each branch has its own visited set, ranking state and candidates; everything any branch scores
        // joins the pool, and the pooled rule below decides -- never the order branches ran in.
        // M3-D: tie branching is rooted at the trunk search, so it never runs after a joint continuation; a continuation that
        // ends in an improving tie fails closed below (improving_move_tie), never choosing a tied sequence.
        const bool improving_tie = reproduced && !report.converged && !report.joint_continued &&
                                   (report.stop_reason == "rescore_incomplete" || report.stop_reason == "tie_at_best");
        std::string tie_ties_text, tie_branches_text;
        std::string tie_incomplete;
        if (improving_tie && options.repair_tie_branching) {
            const Clock::time_point tie_mark = Clock::now();
            write_running_status("tie branching");
            report.tie_branched = true;
            report.search_conditional = "tie_branched";
            PreparedScoringContext* shared = search->context();
            const std::uint64_t shared_before = shared != nullptr ? shared->cumulative_counters().physical_origins : 0;
            const std::uint64_t memo_hits_before = digest_memo_hits, memo_misses_before = digest_memo_misses;
            // Pending seeds. Production keeps each waiting branch as its compact descriptor only -- the BranchStart
            // the search holds: allele paths, declared score components and prior, sequence digests, parent tie and
            // history identities -- and nothing in the context. A branch's start is materialised just before it
            // runs (spelled, placed, censused, its declared score reproduced) and released once it is processed.
            // The eager frontier is a gate hook: a census-holding slot pair per seed, built by block delta at its tie
            // and pinned until its branch has run.
            const bool legacy_lifetime = options.test_retain_branch_parents;   // the pre-fix lifetime, test hook
            const bool eager_seeds = options.test_eager_branch_seeds || legacy_lifetime;
            // The bounded seed cache: seeds built by block delta when their tie is recorded, admitted in canonical tie
            // order while their slots' exact payload fits the budget, held until their branch runs. The eager frontier
            // is the same cache without a budget.
            const std::uint64_t seed_cache_budget =
                eager_seeds ? std::numeric_limits<std::uint64_t>::max() : options.repair_tie_seed_cache_bytes;
            std::map<DiplotypeIdentity, std::pair<std::size_t, std::size_t>> seed_slots;   // cached seeds
            std::set<std::size_t> legacy_parent_pins;
            std::map<std::size_t, std::size_t> seed_pins;        // slot -> cached seeds still waiting on it
            std::map<std::size_t, std::uint64_t> seed_slot_bytes;   // slot -> its payload when first pinned
            std::map<DiplotypeIdentity, std::uint64_t> seed_metadata;   // I4-C: each cached seed's charged bookkeeping
            std::uint64_t seed_cache_bytes = 0;
            std::set<DiplotypeIdentity> announced;          // every start a tie offered, as the worklist dedupes
            std::set<DiplotypeIdentity> pending;            // announced and not yet run: the frontier
            std::unique_ptr<CoordinateRepairSearch> live;   // the branch that ran last, alive for its own tie
            const std::vector<double> seed_switches = run_switches(block_count);   // M4-B
            const auto pins_now = [&]() {
                std::vector<std::size_t> pins(legacy_parent_pins.begin(), legacy_parent_pins.end());
                for (const auto& seed : seed_pins) pins.push_back(seed.first);
                return pins;
            };
            // The branch memory ledger: one row per settle point and per materialised seed, kept as a partial file
            // while branching runs (a refused run keeps it) and finalised with the run.
            std::ostringstream branch_memory;
            branch_memory << "event\tbranch\tlifetime\tslots_total\tslots_live\tslots_released\tsequence_bytes\t"
                             "segment_bytes\tmass_column_bytes\tcensus_count\tcensus_bytes\torigin_cache_payload\t"
                             "origin_cache_capacity\tmerge_working_bytes\tlogical_total\tcapacity_total\tphys_footprint\t"
                             "phys_footprint_peak\tpending_seed_slots\tpending_seeds\tseeding\tpending_descriptors\t"
                             "cumulative_physical_origins\tseed_cache_bytes\n";
            const auto ledger_row = [&](const char* event, long branch) {
                const ContextMemoryAccount account = shared->memory_account();
                const FootprintSample footprint = physical_footprint();
                std::string waiting;
                for (const auto& [slot, count] : seed_pins) {
                    waiting += (waiting.empty() ? "" : ",") + std::to_string(slot) + ":" + std::to_string(count);
                }
                branch_memory << event << '\t' << branch << '\t'
                              << (legacy_lifetime ? "retain_parents" : "branch_local") << '\t' << account.slots_total
                              << '\t' << account.slots_live << '\t' << account.slots_released << '\t'
                              << account.sequence_bytes << '\t' << account.segment_bytes << '\t'
                              << account.mass_column_bytes << '\t' << account.census_count << '\t'
                              << account.census_bytes << '\t' << account.origin_cache_payload << '\t'
                              << account.origin_cache_capacity << '\t' << account.merge_working_bytes << '\t'
                              << account.logical_total << '\t' << account.capacity_total << '\t' << footprint.current
                              << '\t' << footprint.peak << '\t' << seed_pins.size() << '\t'
                              << (waiting.empty() ? "none" : waiting) << '\t'
                              << (eager_seeds ? "eager" : seed_cache_budget > 0 ? "seed_cache" : "just_in_time") << '\t'
                              << pending.size() << '\t' << shared->cumulative_counters().physical_origins << '\t'
                              << seed_cache_bytes << '\n';
                atomic_write(prefix + ".tie_branch_memory.partial.tsv", branch_memory.str());
            };
            const auto settle = [&](long branch) {
                if (shared == nullptr) return;
                if (!legacy_lifetime) {
                    // Everything the pool needs was copied from the finished branch, and a split's children are
                    // recorded; the finished branch's search and every repair slot not seeding a waiting branch go
                    // (census-free: every repair slot).
                    live.reset();
                    shared->release_transient(pins_now());
                }
                ledger_row(branch < 0 ? "root_seeded" : "branch_settled", branch);
            };
            const auto on_tie = [&](const TieSet& tie) {
                if (shared == nullptr) throw std::runtime_error("no production scoring context to branch from");
                std::vector<const TieItem*> fresh;   // starts this tie offers first, in canonical tie order
                for (const TieItem& n : tie.tied) {
                    const DiplotypeIdentity key = diplotype_identity(n.first, n.second);
                    if (announced.insert(key).second) {
                        pending.insert(key);
                        fresh.push_back(&n);
                    }
                }
                report.tie_peak_pending_descriptors = std::max(report.tie_peak_pending_descriptors, pending.size());
                if (seed_cache_budget == 0) return;
                CoordinateRepairSearch& parent = tie.branch < 0 ? *search : *live;
                if (legacy_lifetime) {
                    legacy_parent_pins.insert(parent.incumbent_slot_first());
                    legacy_parent_pins.insert(parent.incumbent_slot_second());
                }
                for (const TieItem* n : fresh) {
                    const DiplotypeIdentity key = diplotype_identity(n->first, n->second);
                    const Clock::time_point seed_mark = Clock::now();
                    const std::uint64_t before = shared->cumulative_counters().physical_origins;
                    // I4-C: the live transient slots before the build, so the build's NEW slots are known exactly (a
                    // seed may resolve to an incumbent's slot or, through the sequence hash, to a slot already held)
                    const std::vector<std::size_t> live_before = eager_seeds ? std::vector<std::size_t>()
                                                                             : shared->live_transient_slots();
                    const auto slots = parent.seed_neighbour_by_delta(n->first, n->second);
                    ++report.tie_seed_cache_built;
                    report.tie_seed_cache_physical_origins += shared->cumulative_counters().physical_origins - before;
                    report.seconds_tie_seed_cache_build += seconds_since(seed_mark);
                    // What admitting it retains: every slot it needs that no admitted seed already holds, exactly.
                    std::uint64_t added = 0;
                    std::map<std::size_t, std::uint64_t> new_slots;
                    for (const std::size_t slot : {slots.first, slots.second}) {
                        if (seed_pins.count(slot) != 0 || new_slots.count(slot) != 0) continue;
                        new_slots[slot] = shared->slot_payload_bytes(slot);
                        added += new_slots[slot];
                    }
                    // I4-C: the cache's own bookkeeping for this seed is charged too -- its seed_slots entry (the
                    // identity's two allele paths) and one pin and one payload record per new slot (tree nodes)
                    const std::uint64_t node = 48;   // red-black node header plus allocator slack, per map entry
                    const std::uint64_t metadata =
                        node + sizeof(std::pair<const DiplotypeIdentity, std::pair<std::size_t, std::size_t>>) +
                        (key.first.capacity() + key.second.capacity()) * sizeof(std::size_t) +
                        new_slots.size() * 2 * (node + sizeof(std::pair<const std::size_t, std::uint64_t>));
                    added += metadata;
                    report.tie_seed_cache_transient_peak_bytes =
                        std::max(report.tie_seed_cache_transient_peak_bytes, seed_cache_bytes + added);
                    if (!eager_seeds && added > seed_cache_budget - seed_cache_bytes) {
                        // I4-C: rejected. Its own new slots go back now (not at the next settle point), except any an
                        // admitted seed or the parent's incumbent holds; the work it executed stays charged. The branch
                        // is materialised just in time as before.
                        ++report.tie_seed_cache_rejected;
                        std::vector<std::size_t> release;
                        for (const auto& [slot, bytes] : new_slots) {
                            if (std::binary_search(live_before.begin(), live_before.end(), slot)) continue;
                            if (slot == parent.incumbent_slot_first() || slot == parent.incumbent_slot_second()) continue;
                            if (legacy_parent_pins.count(slot) != 0) continue;
                            release.push_back(slot);
                            report.tie_seed_cache_rejected_released_bytes += bytes;
                        }
                        if (!release.empty()) shared->release_slots(release);
                        continue;
                    }
                    ++report.tie_seed_cache_admitted;
                    report.tie_seed_cache_metadata_bytes += metadata;
                    seed_cache_bytes += added;
                    seed_slot_bytes.insert(new_slots.begin(), new_slots.end());
                    report.tie_seed_cache_peak_bytes = std::max(report.tie_seed_cache_peak_bytes, seed_cache_bytes);
                    seed_slots.emplace(key, slots);
                    seed_metadata.emplace(key, metadata);
                    ++seed_pins[slots.first];
                    if (slots.second != slots.first) ++seed_pins[slots.second];
                }
                // Every cached seed holds its census until its branch runs.
                report.tie_seed_peak_live_censuses = std::max(report.tie_seed_peak_live_censuses, seed_pins.size());
            };
            // A branch start, just materialised or taken from the cache, must reproduce its declared complete score --
            // every component and the prior -- bit for bit before any neighbour is scored; otherwise the branching
            // fails closed as tie_branching_incomplete.
            const auto verify_start = [&](const BranchStart& branch, std::size_t slot_a, std::size_t slot_b,
                                          const std::string& a, const std::string& b, const std::string& which) {
                const FragmentScoringResult scored =
                    shared->score_candidates({{std::min(slot_a, slot_b), std::max(slot_a, slot_b)}});
                const DiplotypeScore& score = scored.scores.front();
                const double prior = log_unordered_sequence_pair_prior(model, a, b, seed_switches);
                DiplotypeScore declared = branch.start.score;
                if (options.test_perturb_branch_seed == static_cast<long>(branch.id)) {
                    // Gate hook: a declared score one ulp off, which the reproduction guard must refuse.
                    declared.log_likelihood = std::nextafter(declared.log_likelihood, 0.0);
                }
                if (!(score.log_likelihood == declared.log_likelihood && score.fragment_term == declared.fragment_term &&
                      score.exposure_term == declared.exposure_term &&
                      score.absolute_term_sum == declared.absolute_term_sum &&
                      score.fragments_with_mass == declared.fragments_with_mass && prior == branch.start.prior)) {
                    std::ostringstream why;
                    why << std::setprecision(17) << which << " did not reproduce its declared score: declared "
                        << declared.log_likelihood << " + " << branch.start.prior << ", scored " << score.log_likelihood
                        << " + " << prior;
                    throw std::runtime_error(why.str());
                }
            };
            const auto spell_start = [&](const BranchStart& branch, const std::string& which) {
                HaplotypeSequence a = spell_mosaic(model, MosaicPath{branch.start.first});
                HaplotypeSequence b = spell_mosaic(model, MosaicPath{branch.start.second});
                if (!(sequence_digest(a.bases) == branch.start.first_sequence) ||
                    !(sequence_digest(b.bases) == branch.start.second_sequence)) {
                    throw std::runtime_error(which + " does not spell its recorded sequences");
                }
                return std::make_pair(std::move(a.bases), std::move(b.bases));
            };
            // Just in time: the branch's two haplotypes spelled from its allele paths, placed in the shared context by
            // the ordinary exact path (so the work counts against the cumulative caps), their censuses retained, and
            // the declared complete score reproduced bit for bit before any neighbour is scored. Any failure throws,
            // which fails the branching closed as tie_branching_incomplete.
            const auto materialise = [&](const BranchStart& branch) {
                const std::string which = "branch " + std::to_string(branch.id) + " seed";
                const Clock::time_point seed_mark = Clock::now();
                const FragmentScoringCounters before = shared->cumulative_counters();
                const auto spelled = spell_start(branch, which);
                std::vector<std::size_t> slots;
                IncrementalBatchReport placed;
                try {
                    slots = shared->add_haplotypes({spelled.first, spelled.second}, &placed, true);
                    shared->retain_census(slots[0]);
                    if (slots[1] != slots[0]) shared->retain_census(slots[1]);
                } catch (const std::exception& error) {
                    throw std::runtime_error(which + " could not be materialised: " + error.what());
                }
                const FragmentScoringCounters& after = shared->cumulative_counters();
                ++report.tie_seeds_materialised;
                report.tie_seed_haplotypes_placed += placed.haplotypes_new;
                report.tie_seed_physical_origins += after.physical_origins - before.physical_origins;
                report.tie_seed_state_classes += after.statistical_state_classes - before.statistical_state_classes;
                report.tie_seed_piece_hits += after.piece_hits - before.piece_hits;
                report.tie_seed_latent_pairs += after.latent_alignment_pairs - before.latent_alignment_pairs;
                report.seconds_tie_seed_placement += seconds_since(seed_mark);
                verify_start(branch, slots[0], slots[1], spelled.first, spelled.second, which);
                // The running branch's start and the cached seeds are the only seed censuses in the context.
                const std::size_t live_seeds = (slots[1] != slots[0] ? 2 : 1) + seed_pins.size();
                report.tie_seed_peak_live_censuses = std::max(report.tie_seed_peak_live_censuses, live_seeds);
                ledger_row("seed_materialised", static_cast<long>(branch.id));
                return std::make_pair(slots[0], slots[1]);
            };
            RepairSearchOptions branch_options = search_options;
            branch_options.external_context = true;
            branch_options.memory_log = false;
            branch_options.memory_live_path.clear();               // S2-H/S2-J: production-search only
            branch_options.round1_reuse_only_admission = false;
            const auto run_branch = [&](const BranchStart& branch) {
                live.reset();
                BranchResult out;
                const Clock::time_point mark = Clock::now();
                const std::uint64_t before = shared->cumulative_counters().physical_origins;
                const DiplotypeIdentity key = diplotype_identity(branch.start.first, branch.start.second);
                pending.erase(key);
                const auto cached = seed_slots.find(key);
                const bool hit = cached != seed_slots.end();
                std::pair<std::size_t, std::size_t> seeded;
                if (hit) {
                    ++report.tie_seed_cache_hits;
                    seeded = cached->second;
                    const std::string which = "branch " + std::to_string(branch.id) + " cached seed";
                    const auto spelled = spell_start(branch, which);
                    verify_start(branch, seeded.first, seeded.second, spelled.first, spelled.second, which);
                    ledger_row("seed_cached", static_cast<long>(branch.id));
                } else {
                    seeded = materialise(branch);
                }
                live = std::make_unique<CoordinateRepairSearch>(model, emissions, fragments, library, repair_scoring,
                                                                branch_options);
                live->attach_context(shared);
                live->set_pinned(pins_now());
                live->start_from_slots(branch.start.first, branch.start.second,
                                       branch.start.score.log_likelihood + branch.start.prior, seeded.first, seeded.second);
                for (const DiplotypeIdentity& seen : branch.history) live->add_visited(seen.first, seen.second);
                RepairHooks hooks;
                hooks.on_scored = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b,
                                      const std::string& sa, const std::string& sb, const DiplotypeScore& score,
                                      double prior) {
                    ++report.pool_branch_scored;
                    const std::size_t round = live->rounds().size() + 1;
                    pool_scored(a, b, sa, sb, score, prior, "branch",
                                "branch" + std::to_string(branch.id) + "_round" + std::to_string(round), &out.items,
                                round, true);
                };
                if (audit_on) hooks.on_scored_views = record_views;
                hooks.on_round = [&](const RepairRound& round) {
                    write_running_status("tie branch " + std::to_string(branch.id) + " round " +
                                         std::to_string(round.index) + " committed");
                };
                const bool branch_reproduced = live->run(hooks);
                accumulate_junction(report.junction, live->junction_report());   // M9-G
                if (hit) {
                    // Consumed: the seed leaves the cache, and each slot no other cached seed holds with it.
                    seed_cache_bytes -= seed_metadata.at(key);   // I4-C: its charged bookkeeping
                    seed_metadata.erase(key);
                    seed_slots.erase(key);
                    for (const std::size_t slot : {seeded.first, seeded.second}) {
                        const auto pin = seed_pins.find(slot);
                        if (pin != seed_pins.end() && --pin->second == 0) {
                            seed_pins.erase(pin);
                            seed_cache_bytes -= seed_slot_bytes.at(slot);
                            seed_slot_bytes.erase(slot);
                        }
                        if (seeded.second == seeded.first) break;
                    }
                }
                out.rounds = live->rounds().size();
                for (const RepairRound& round : live->rounds()) {
                    if (!round.moved) continue;
                    ++out.moves;
                    out.moved_to.push_back(diplotype_identity(round.first, round.second));
                }
                out.stop = live->stop_reason();
                out.last_round = live->rounds().empty() ? 0 : live->rounds().back().index;
                out.terminal_first = live->incumbent_first();
                out.terminal_second = live->incumbent_second();
                out.terminal_total = live->incumbent_total();
                out.executed_physical_origins = shared->cumulative_counters().physical_origins - before;
                out.seconds = seconds_since(mark);
                if (!branch_reproduced) {
                    out.outcome = BranchOutcome::Incomplete;
                    out.why = "its start did not reproduce";
                } else if (live->converged()) {
                    out.outcome = BranchOutcome::Converged;
                } else if (out.stop == "rescore_incomplete" || out.stop == "tie_at_best") {
                    out.outcome = BranchOutcome::Split;
                } else {
                    out.outcome = BranchOutcome::Incomplete;
                    out.why = "stopped at " + out.stop;
                }
                if (audit_on) {
                    CentralBranch record;
                    record.result = out;
                    record.rounds = live->rounds();
                    record.terminal_totals.assign(view_count, std::numeric_limits<double>::quiet_NaN());
                    if (!live->rounds().empty() && !live->rounds().back().verdicts.empty()) {
                        for (std::size_t v = 0; v < view_count; ++v) {
                            const RoundVerdict& last = live->rounds().back().verdicts[v];
                            record.terminal_totals[v] = last.moved ? last.winner_total : last.incumbent_total;
                        }
                    }
                    central_branches[{key, branch.history}] = std::move(record);
                }
                return out;
            };
            TieBranchingReport branched;
            try {
                const TieSet root = terminal_tie(production_items, search->rounds().back().index, fragments.size(),
                                                 scoring_options.mismatch_band);
                std::vector<DiplotypeIdentity> root_history;
                root_history.push_back(diplotype_identity(anchor_diplotype.first.alleles, anchor_diplotype.second.alleles));
                for (const RepairRound& round : search->rounds()) {
                    if (round.moved) root_history.push_back(diplotype_identity(round.first, round.second));
                }
                branched = run_tie_branching(root, root_history, options.repair_tie_max_branches, fragments.size(),
                                             scoring_options.mismatch_band, run_branch, on_tie,
                                             options.test_reverse_tie_branch_order, settle);
            } catch (const std::exception& error) {
                branched.incomplete = std::string("tie branching failed: ") + error.what();
            }
            live.reset();
            // Every exit -- converged, capped, refused or failed -- closes with no seed slot or census retained;
            // the content-keyed origin and emission caches stay for the representation guard.
            if (shared != nullptr && !shared->refused()) {
                seed_pins.clear();
                seed_slots.clear();
                seed_slot_bytes.clear();
                seed_metadata.clear();
                seed_cache_bytes = 0;
                legacy_parent_pins.clear();
                shared->release_transient({});
                ledger_row("branching_closed", -1);
            }
            atomic_write(prefix + ".tie_branch_memory.tsv", branch_memory.str());
            std::remove((prefix + ".tie_branch_memory.partial.tsv").c_str());
            tie_incomplete = branched.incomplete;
            report.tie_ties = branched.ties.size();
            report.tie_branches_run = branched.branches.size();
            report.tie_distinct_terminal_genotypes = branched.terminal_genotypes.size();
            report.tie_incomplete_reason = tie_incomplete.empty() ? "none" : tie_incomplete;
            report.tie_branch_physical_origins =
                shared != nullptr ? shared->cumulative_counters().physical_origins - shared_before : 0;
            report.tie_digest_memo_hits = digest_memo_hits - memo_hits_before;
            report.tie_digest_memo_misses = digest_memo_misses - memo_misses_before;
            if (shared != nullptr) {
                report.cumulative = shared->cumulative_counters();   // the caps bind production and branches together
                report.origin_cache_bytes_end = shared->origin_cache_bytes();
                report.fine_window = shared->fine_window_report();
            }
            report.pool_identities = pool.size();
            report.pool_repair_only = 0;
            for (const PoolEntry& entry : pool) report.pool_repair_only += entry.in_union ? 0 : 1;
            report.pool_score_mismatches = mismatches;
            report.seconds_tie_branching = seconds_since(tie_mark);
            // The two tables, rendered now and written with the final outputs.
            const auto genotype_of = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
                std::string text;
                for (std::size_t block = 0; block < block_count; ++block) {
                    const auto pair = std::minmax(a[block], b[block]);
                    text += (block == 0 ? "" : ",") + std::to_string(pair.first) + "/" + std::to_string(pair.second);
                }
                return text;
            };
            std::ostringstream ties_out;
            ties_out << std::setprecision(17)
                     << "tie\tbranch\tneighbour\thap1_alleles\thap2_alleles\tunordered_genotype\trepaired_block\t"
                        "allele_x\tallele_y\tassignment\ttotal\tlog_likelihood\tfragment_term\texposure_term\tprior\t"
                        "incumbent_total\timprovement\n";
            for (const TieSet& tie : branched.ties) {
                const double inc_total = tie.incumbent.score.log_likelihood + tie.incumbent.prior;
                for (std::size_t k = 0; k < tie.tied.size(); ++k) {
                    const TieItem& n = tie.tied[k];
                    std::size_t block = 0;
                    for (std::size_t b = 0; b < block_count; ++b) {
                        if (n.first[b] != tie.incumbent.first[b] || n.second[b] != tie.incumbent.second[b]) block = b;
                    }
                    const auto pair = std::minmax(n.first[block], n.second[block]);
                    const DiplotypeIdentity key = diplotype_identity(n.first, n.second);
                    const double total = n.score.log_likelihood + n.prior;
                    ties_out << tie.id << '\t' << tie.branch << '\t' << k << '\t' << join_alleles(key.first) << '\t'
                             << join_alleles(key.second) << '\t' << genotype_of(n.first, n.second) << '\t'
                             << model.blocks()[block].block_id << '\t' << pair.first << '\t' << pair.second << '\t'
                             << (n.first[block] == pair.first ? 0 : 1) << '\t' << total << '\t'
                             << n.score.log_likelihood << '\t' << n.score.fragment_term << '\t' << n.score.exposure_term
                             << '\t' << n.prior << '\t' << inc_total << '\t' << total - inc_total << '\n';
                }
            }
            tie_ties_text = ties_out.str();
            std::ostringstream branches_out;
            branches_out << std::setprecision(17)
                         << "branch\tparent_tie\tstart_hap1\tstart_hap2\tstart_total\trounds\tmoves\tstop\toutcome\t"
                            "terminal_hap1\tterminal_hap2\tterminal_genotype\tterminal_total\texecuted_physical_origins\t"
                            "seconds\n";
            for (const BranchRecord& record : branched.branches) {
                const DiplotypeIdentity start = diplotype_identity(record.start.start.first, record.start.start.second);
                const BranchResult& r = record.result;
                const bool ended = !r.terminal_first.empty();
                const DiplotypeIdentity end = ended ? diplotype_identity(r.terminal_first, r.terminal_second) : start;
                branches_out << record.start.id << '\t' << record.start.parent_tie << '\t' << join_alleles(start.first)
                             << '\t' << join_alleles(start.second) << '\t'
                             << record.start.start.score.log_likelihood + record.start.start.prior << '\t' << r.rounds
                             << '\t' << r.moves << '\t' << r.stop << '\t'
                             << (r.outcome == BranchOutcome::Converged ? "converged"
                                 : r.outcome == BranchOutcome::Split   ? "split"
                                                                       : "incomplete")
                             << '\t' << join_alleles(end.first) << '\t' << join_alleles(end.second) << '\t'
                             << genotype_of(end.first, end.second) << '\t' << r.terminal_total << '\t'
                             << r.executed_physical_origins << '\t' << r.seconds << '\n';
            }
            tie_branches_text = branches_out.str();
        }
        // S2-K: no claim here -- the post-refinement scoring, audits and report preparation stay under the watchdog
        note_footprint("post_refinement_start");
        if (report.tie_branched) {
            defer_output(prefix + ".tie_terminal_ties.tsv") << tie_ties_text;
            defer_output(prefix + ".tie_branches.tsv") << tie_branches_text;
        }
        if (report.tie_branched && !tie_incomplete.empty()) {
            report.no_call_reason = "tie_branching_incomplete";
            write_trajectory();
            return stop_with(GenotypeRunStatus::RepairIncomplete,
                             "improving-move tie branching did not complete (" + tie_incomplete +
                                 "); the trajectory and the branch tables are kept and no genotype is emitted",
                             nullptr, nullptr);
        }
        if (!reproduced) {
            report.no_call_reason = "anchor_not_reproduced";
            write_trajectory();
            std::ostringstream why;
            why << std::setprecision(17) << "the repair anchor did not reproduce its union score in round "
                << search->failed_round() << " (expected " << search->incumbent_expected() << ", scored "
                << search->incumbent_scored() << "); no genotype";
            return stop_with(GenotypeRunStatus::RepairIncomplete, why.str(), nullptr, nullptr);
        }
        if (!report.converged && !report.tie_branched) {
            // A round tie that excludes the incumbent is a tie among improving moves: the search cannot
            // choose between them and fails closed (v1, or tie branching off). Every other non-converged stop
            // keeps its own name.
            const std::string& stop = report.stop_reason;
            report.no_call_reason =
                stop == "rescore_incomplete" || stop == "tie_at_best" ? "improving_move_tie" : stop;
            write_trajectory();
            if (stop == "round_cap") {
                provisional("round_cap", "the repair search reached its declared round cap without converging; the last "
                                         "committed incumbent, exactly scored; search incomplete");
            }
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
                        branch_options.memory_live_path.clear();               // S2-H/S2-J: production-search only
                        branch_options.round1_reuse_only_admission = false;
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
        if (options.block_support_rule != "legacy") {   // M3-E (EXPERIMENTAL): the shared full-data-likelihood support
            std::vector<SupportGroup> groups(decided.groups.size());
            for (std::size_t g = 0; g < decided.groups.size(); ++g) {
                std::set<BlockGenotypeRow> rows;
                for (const std::size_t m : decided.groups[g]) {
                    BlockGenotypeRow row;
                    for (std::size_t b = 0; b < block_count; ++b) {
                        row.push_back(std::minmax(decided.candidates[m].first[b], decided.candidates[m].second[b]));
                    }
                    rows.insert(row);
                }
                groups[g].rows.assign(rows.begin(), rows.end());
            }
            summary.data_support_competitor_paths.assign(decided.groups.size(), "NA");
            for (const CertifiedCandidate& c : decided.ranking.candidates) {
                SupportGroup& g = groups[c.candidate_id];
                g.fragment = c.fragment_term;
                g.exposure = c.exposure_term;
                g.prior = c.prior;
                g.bound = c.numeric_error_bound;
                const PoolCandidate& rep = decided.candidates[decided.representative[c.candidate_id]];
                summary.data_support_competitor_paths[c.candidate_id] = join_alleles(rep.first) + "|" + join_alleles(rep.second);
            }
            if (options.block_support_rule == "full_data_likelihood_v2_experimental") {
                summary.data_support_v2 = block_data_support_v2(groups, best.candidate_id, block_count);
            } else {
                summary.data_support = block_data_support_v1(groups, best.candidate_id, block_count);
            }
        }

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

        note_footprint("pool_decision_and_summary");   // S2-K
        // ---- outputs, prepared in memory (S2-K): published by flush_outputs once the run is claimed.
        {
            std::vector<std::size_t> order(decided.candidates.size());
            for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
            std::sort(order.begin(), order.end(), [&](std::size_t x, std::size_t y) { return table_id[x] < table_id[y]; });
            std::ostringstream& out = defer_output(prefix + ".mosaic_candidates.tsv");
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
            std::ostringstream& out = defer_output(prefix + ".mosaic_scores.tsv");
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
            std::ostringstream& out = defer_output(prefix + ".equivalence_set.tsv");
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
        // S1-B (sequence_first_v1 only): the sequence report of this completed ranking. Reporting only: it reads the
        // decision and spells the representative; nothing below depends on it.
        if (sequence_first(options)) {
            std::vector<std::string> block_ids;
            for (const auto& block : model.blocks()) block_ids.push_back(block.block_id);
            g_sequence_report = build_sequence_report(decided, block_ids, [&](const std::vector<std::size_t>& labels) {
                return spell_mosaic(model, MosaicPath{labels}).bases;
            });
            if (options.output_schema == "sequence_first_v2") annotate_v2(g_sequence_report, model, fingerprint);
            compound_pooled(g_sequence_report, decided, model, fingerprint, options);   // M9-J (opt-in)
        }

        if (decided.decision == TopSetDecision::DistinctTie) {
            // After tie branching the pooled top set spans the branches: disagreement there is a divergent tie.
            report.no_call_reason = report.tie_branched ? "divergent_tie" : "distinct_tie";
            return stop_with(GenotypeRunStatus::RescoreIncomplete,
                             "the pooled top set holds " + std::to_string(decided.top_groups.size()) +
                                 " candidates tied on total score that disagree at block(s) " + report.disagreeing_blocks +
                                 "; the equivalence set is reported and no genotype is called",
                             &decided.ranking, &summary);
        }
        const bool unique = decided.decision == TopSetDecision::Unique;
        const PoolCandidate& called = decided.candidates[decided.representative[best.candidate_id]];

        // ---- the representation-conflict guard: an annotation of the call just decided, never a change to it.
        // For each called block: the marker (depth) channel's best available pair at that block, its margin
        // over the called pair, and -- when that margin is decisive -- the called diplotype with that pair
        // substituted at the block (both assignments, the better by total), scored exactly (from the pool when
        // it is there, else in the live context, outside the pool). The block is LOWQUAL when the marker margin,
        // the fragment term's preference for the call and the exposure term's preference for the competitor
        // each reach ln(1000) nats: two depth readings against the call, the sequence term alone for it.
        struct BlockFilter {
            std::string filter = "PASS", reason = "none", competitor = "NA";
            double marker_margin = 0.0;
            double fragment_diff = std::numeric_limits<double>::quiet_NaN();
            double exposure_diff = std::numeric_limits<double>::quiet_NaN();
            double likelihood_margin = std::numeric_limits<double>::quiet_NaN();
            double statistic = std::numeric_limits<double>::quiet_NaN();
        };
        std::vector<BlockFilter> block_filter(block_count);
        {
            const Clock::time_point guard_mark = Clock::now();
            const double bound = std::isnan(options.test_representation_guard_bound_nats)
                                      ? kRepresentationBoundNats
                                      : options.test_representation_guard_bound_nats;
            PreparedScoringContext* guard_context = search ? search->context() : nullptr;
            const std::uint64_t guard_origins_before =
                guard_context != nullptr ? guard_context->cumulative_counters().physical_origins : 0;
            const std::vector<double> guard_switches = run_switches(block_count);   // M4-B
            for (std::size_t b = 0; b < block_count && options.representation_guard; ++b) {
                BlockFilter& bf = block_filter[b];
                const BlockProposalEmissions& emission = emissions.blocks[b];
                const auto pair = std::minmax(called.first[b], called.second[b]);
                std::pair<std::size_t, std::size_t> top = marker_best_pair(emission);
                if (options.test_representation_guard_runner_up && top == std::make_pair(pair.first, pair.second) &&
                    emission.allele_count > 1) {
                    top = marker_best_pair(emission, top);
                }
                bf.marker_margin = emission.at(top.first, top.second) - emission.at(pair.first, pair.second);
                if (top == std::make_pair(pair.first, pair.second) || !(bf.marker_margin >= bound)) continue;
                bf.competitor = std::to_string(top.first) + "/" + std::to_string(top.second);
                std::vector<RepresentationAssignment> assignments;
                std::string failure;
                for (int assignment = 0; assignment < 2; ++assignment) {
                    if (assignment == 1 && top.first == top.second) break;
                    std::vector<std::size_t> first = called.first;
                    std::vector<std::size_t> second = called.second;
                    first[b] = assignment == 0 ? top.first : top.second;
                    second[b] = assignment == 0 ? top.second : top.first;
                    DiplotypeScore score;
                    double prior = 0.0;
                    const auto held = options.test_representation_guard_force_context
                                          ? pool_at.end()
                                          : pool_at.find(diplotype_identity(first, second));
                    if (held != pool_at.end()) {
                        score = pool[held->second].candidate.score;
                        prior = pool[held->second].candidate.sequence_prior;
                    } else if (guard_context != nullptr) {
                        try {
                            const std::string sa = spell_mosaic(model, MosaicPath{first}).bases;
                            const std::string sb = spell_mosaic(model, MosaicPath{second}).bases;
                            const std::vector<std::size_t> slots = guard_context->add_haplotypes({sa, sb});
                            const FragmentScoringResult scored_pair =
                                guard_context->score_candidates({{std::min(slots[0], slots[1]), std::max(slots[0], slots[1])}});
                            score = scored_pair.scores.front();
                            prior = log_unordered_sequence_pair_prior(model, sa, sb, guard_switches);
                            ++report.representation_competitors_scored;
                        } catch (const std::exception& error) {
                            failure = error.what();
                            continue;
                        }
                    } else {
                        failure = "no scoring context";
                        continue;
                    }
                    assignments.push_back({score, prior});
                }
                if (assignments.empty()) {
                    bf.filter = "NOT_EVALUATED";
                    bf.reason = "competitor_not_scored: " + failure;
                    continue;
                }
                const RepresentationEvidence evidence =
                    representation_conflict(bf.marker_margin, called.score, called.sequence_prior, assignments, bound);
                bf.fragment_diff = evidence.fragment_diff;
                bf.exposure_diff = evidence.exposure_diff;
                bf.likelihood_margin = evidence.likelihood_margin;
                bf.statistic = evidence.statistic;
                if (evidence.conflict) {
                    bf.filter = "LOWQUAL_REPRESENTATION_CONFLICT";
                    bf.reason = "marker_and_exposure_prefer_competitor_fragment_term_prefers_call";
                }
            }
            if (!options.representation_guard) {
                for (BlockFilter& bf : block_filter) {
                    bf.filter = "NOT_EVALUATED";
                    bf.reason = "guard_off";
                }
            }
            // Fail closed: PASS only when every block completed the rule and passed. A flag wins the locus;
            // otherwise any unevaluated block makes it NOT_EVALUATED. Both facts are kept when they coexist.
            report.lowqual_blocks = "none";
            report.not_evaluated_blocks = "none";
            bool any_lowqual = false;
            bool any_unevaluated = false;
            for (std::size_t b = 0; b < block_count; ++b) {
                const std::string& id = model.blocks()[b].block_id;
                if (block_filter[b].filter == "LOWQUAL_REPRESENTATION_CONFLICT") {
                    any_lowqual = true;
                    report.lowqual_blocks = (report.lowqual_blocks == "none" ? std::string() : report.lowqual_blocks + ",") + id;
                } else if (block_filter[b].filter != "PASS") {
                    any_unevaluated = true;
                    report.not_evaluated_blocks =
                        (report.not_evaluated_blocks == "none" ? std::string() : report.not_evaluated_blocks + ",") + id;
                }
            }
            report.locus_filter = any_lowqual       ? "LOWQUAL_REPRESENTATION_CONFLICT"
                                  : any_unevaluated ? "NOT_EVALUATED"
                                                    : "PASS";
            report.representation_guard_complete = any_unevaluated ? "false" : "true";
            // The guard's own scoring, reported apart, and the cumulative repair counters refreshed so the
            // reported totals and margins include it (the caps bind it with everything else).
            if (guard_context != nullptr) {
                report.representation_guard_physical_origins =
                    guard_context->cumulative_counters().physical_origins - guard_origins_before;
                report.cumulative = guard_context->cumulative_counters();
                report.origin_cache_bytes_end = guard_context->origin_cache_bytes();
                report.fine_window = guard_context->fine_window_report();
            }
            report.seconds_representation_guard = seconds_since(guard_mark);
            note_footprint("representation_guard_done");   // S2-K
            if (options.test_interrupt_point == "post_refinement_guard") {   // S2-K TEST HOOK: a time stop after the guard
                stop_run("TIME_LIMIT", "test hook: post-refinement, after the representation guard", 124);
            }
        }
        // ---- library_sensitivity_audit_v2 (docs/algorithms/genotype.md, stage 9). The call is decided and
        // annotated; the same algorithm is now followed at 0.9x and 1.1x the estimated fragment-start intensity. Every
        // candidate the call scored already carries its score at both (one combination pass), so an endpoint that
        // decides as the call did -- its anchor, every round, every branch -- adds no scoring at all. Where it decides
        // differently it continues alone, in the shared context: what that holds is reused, only what is new is
        // placed. The common pool is then reranked at each view, and only unordered genotypes are compared.
        std::set<std::size_t> library_sensitive;
        bool library_not_evaluated = options.library_path == "auto" && !audit_on;   // auto without a completed audit
        if (audit_on) {
            const Clock::time_point audit_mark = Clock::now();
            PreparedScoringContext* shared = search->context();
            const FragmentScoringCounters work_before = shared->cumulative_counters();
            const std::size_t created_before = shared->haplotypes_created();
            const std::uint64_t dps_before = read_emission_cache_counters().executed_dps;
            const std::size_t band = scoring_options.mismatch_band;
            const std::vector<double> start_switches = run_switches(block_count);   // M4-B
            const char* const view_name[] = {"central", "low", "high"};
            std::vector<LibraryAuditView> views(view_count);
            std::vector<std::string> fork_of(view_count, "none"), trunk_stop(view_count, report.stop_reason);
            std::vector<std::size_t> reused_of(view_count, 0), run_of(view_count, 0);
            std::vector<TieItem> endpoint_items;                  // what the endpoints scored themselves
            // DIAGNOSTIC (recorded, never the verdict): each endpoint's own pool -- its union set, its trunk and every
            // branch it used -- which is what an independent run at that intensity pools and decides over.
            std::vector<std::vector<TieItem>> own_items(view_count);
            std::vector<std::unique_ptr<SetRanking>> endpoint_sets;   // the union set each endpoint ranked
            std::vector<const SetRanking*> set_of_view(view_count, nullptr);
            const auto identity_text = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
                const DiplotypeIdentity key = diplotype_identity(a, b);
                return join_alleles(key.first) + "|" + join_alleles(key.second);
            };
            views[0].anchor = identity_text(anchor_diplotype.first.alleles, anchor_diplotype.second.alleles);
            const DiplotypeIdentity central_anchor =
                diplotype_identity(anchor_diplotype.first.alleles, anchor_diplotype.second.alleles);
            // The union ranked at one view by the production rule, the full-beam tie audit included. Ranking is
            // bookkeeping only; the run's own accumulated counters and timings are restored.
            const auto rank_view = [&](std::size_t view) {
                const StageTimings saved = run_timings;
                const std::size_t saved_evaluations = haploid_prior_evaluations, saved_lookups = diploid_prior_lookups;
                const CandidateSetReport saved_report = candidate_report;
                auto ranked = std::make_unique<SetRanking>(rank_set(beam, view));
                if (ranked->ranking.status != RescoreStatus::Complete && mosaics.search_exhaustive &&
                    mosaics.candidates.size() > beam.candidates.size()) {
                    ranked = std::make_unique<SetRanking>(rank_set(mosaics, view));
                }
                run_timings = saved;
                haploid_prior_evaluations = saved_evaluations;
                diploid_prior_lookups = saved_lookups;
                candidate_report = saved_report;
                return ranked;
            };
            // An endpoint start placed in the shared context, its census retained, and its declared score at that
            // view reproduced bit for bit before any neighbour is scored.
            const auto materialise_at = [&](std::size_t view, const TieItem& start) {
                const std::string a = spell_mosaic(model, MosaicPath{start.first}).bases;
                const std::string b = spell_mosaic(model, MosaicPath{start.second}).bases;
                const std::vector<std::size_t> slots = shared->add_haplotypes({a, b}, nullptr, true);
                shared->retain_census(slots[0]);
                if (slots[1] != slots[0]) shared->retain_census(slots[1]);
                const FragmentScoringResult again =
                    shared->score_candidates({{std::min(slots[0], slots[1]), std::max(slots[0], slots[1])}});
                const DiplotypeScore& got = view == 0 ? again.scores.front() : again.audit_scores.at(view - 1).front();
                if (!same_score(got, start.views.at(view)) ||
                    log_unordered_sequence_pair_prior(model, a, b, start_switches) != start.prior) {
                    throw std::runtime_error("an endpoint start did not reproduce its declared score");
                }
                return std::make_pair(slots[0], slots[1]);
            };
            // One leg of an endpoint's own trajectory: the exact repair at that view, attached to the shared context.
            struct Leg {
                bool reproduced = false, converged = false;
                std::string stop = "none";
                std::size_t rounds = 0, last_round = 0;
                std::vector<DiplotypeIdentity> moves;
                std::vector<std::size_t> terminal_first, terminal_second;
                double terminal_total = 0.0;
            };
            const auto run_leg = [&](std::size_t view, const TieItem* start, const MosaicDiplotype* anchor_start,
                                     double anchor_total, std::size_t first_round,
                                     const std::vector<DiplotypeIdentity>& visited, std::vector<TieItem>& items) {
                RepairSearchOptions leg_options = search_options;
                leg_options.memory_live_path.clear();               // S2-H/S2-J: production-search only
                leg_options.round1_reuse_only_admission = false;
                leg_options.external_context = true;
                leg_options.memory_log = false;
                leg_options.score_view = view;
                leg_options.first_round = first_round;
                leg_options.max_rounds = options.repair_rounds + 1 - first_round;
                CoordinateRepairSearch leg_search(model, emissions, fragments, library, repair_scoring, leg_options);
                leg_search.attach_context(shared);
                if (start != nullptr) {
                    const auto slots = materialise_at(view, *start);
                    leg_search.start_from_slots(start->first, start->second,
                                                start->views.at(view).log_likelihood + start->prior, slots.first,
                                                slots.second);
                } else {
                    leg_search.start(anchor_start->first.alleles, anchor_start->second.alleles, anchor_total);
                }
                for (const DiplotypeIdentity& seen : visited) leg_search.add_visited(seen.first, seen.second);
                RepairHooks hooks;
                hooks.on_scored = [&](const std::vector<std::size_t>& a, const std::vector<std::size_t>& b,
                                      const std::string& sa, const std::string& sb, const DiplotypeScore& score,
                                      double prior) {
                    items.push_back({first_round + leg_search.rounds().size(), a, b, digest_of(a, sa, true),
                                     digest_of(b, sb, true), score, prior, {}});
                };
                hooks.on_scored_views = [&](const std::vector<DiplotypeScore>& all) { items.back().views = all; };
                Leg leg;
                leg.reproduced = leg_search.run(hooks);
                accumulate_junction(report.junction, leg_search.junction_report());   // M9-G
                leg.converged = leg_search.converged();
                leg.stop = leg_search.stop_reason();
                leg.rounds = leg_search.rounds().size();
                leg.last_round = leg_search.rounds().empty() ? first_round : leg_search.rounds().back().index;
                for (const RepairRound& round : leg_search.rounds()) {
                    if (round.moved) leg.moves.push_back(diplotype_identity(round.first, round.second));
                }
                leg.terminal_first = leg_search.incumbent_first();
                leg.terminal_second = leg_search.incumbent_second();
                leg.terminal_total = leg_search.incumbent_total();
                shared->release_transient({});
                return leg;
            };
            const auto at_view = [](TieItem item, std::size_t view) {
                item.score = item.views.at(view);
                return item;
            };
            // A recorded trajectory -- the call's trunk or one of its branches -- continued at view k: shared while
            // every round decides at k as it did, then, at the first round that does not, either stopped there (its
            // items are already scored) or forked from k's winner as a leg of its own. Nothing shared is scored again.
            struct Continued {
                std::vector<TieItem> items;             // at view k: the shared rounds', then the leg's
                std::vector<DiplotypeIdentity> moves;   // the incumbents moved to, in order
                std::string stop = "none";
                bool in_top = false, forked = false, leg = false;
                std::size_t rounds = 0, last_round = 0, fork_round = 0;
                std::vector<std::size_t> terminal_first, terminal_second;
                double terminal_total = 0.0;
            };
            const auto continue_at = [&](std::size_t k, const std::vector<RepairRound>& rounds,
                                         const std::vector<TieItem>& items, const std::string& recorded_stop,
                                         const std::vector<std::size_t>& start_first,
                                         const std::vector<std::size_t>& start_second,
                                         const std::vector<std::size_t>& recorded_first,
                                         const std::vector<std::size_t>& recorded_second, double recorded_total,
                                         std::vector<DiplotypeIdentity> visited) {
                Continued c;
                std::size_t r = npos;
                for (std::size_t i = 0; i < rounds.size() && r == npos; ++i) {
                    if (!rounds[i].verdicts.empty() && !same_decision(rounds[i].verdicts[k], rounds[i].verdicts[0])) r = i;
                }
                const std::size_t shared_rounds = r == npos ? rounds.size() : r + 1;
                const std::size_t shared_through = shared_rounds == 0 ? 0 : rounds[shared_rounds - 1].index;
                for (const TieItem& item : items) {
                    if (item.round <= shared_through) c.items.push_back(at_view(item, k));
                }
                std::vector<std::size_t> first = start_first, second = start_second;
                for (std::size_t i = 0; i < (r == npos ? rounds.size() : r); ++i) {
                    if (!rounds[i].moved) continue;
                    c.moves.push_back(diplotype_identity(rounds[i].first, rounds[i].second));
                    first = rounds[i].first;
                    second = rounds[i].second;
                }
                c.rounds = shared_rounds;
                if (r == npos) {
                    c.stop = recorded_stop;
                    c.in_top = !rounds.empty() && rounds.back().incumbent_in_top_set;
                    c.last_round = rounds.empty() ? 0 : rounds.back().index;
                    c.terminal_first = recorded_first;
                    c.terminal_second = recorded_second;
                    c.terminal_total = recorded_total;
                    return c;
                }
                const RoundVerdict& verdict = rounds[r].verdicts[k];
                c.forked = true;
                c.fork_round = rounds[r].index;
                c.last_round = rounds[r].index;
                c.terminal_first = first;
                c.terminal_second = second;
                c.terminal_total = verdict.incumbent_total;
                if (!verdict.moved) {
                    c.stop = verdict.stop;
                    c.in_top = verdict.incumbent_in_top_set;
                    return c;
                }
                const TieItem* winner = nullptr;
                for (const TieItem& item : items) {
                    if (item.round == rounds[r].index && diplotype_identity(item.first, item.second) == verdict.winner) {
                        winner = &item;
                        break;
                    }
                }
                if (winner == nullptr) throw std::logic_error("the endpoint's winning move was not recorded");
                c.moves.push_back(verdict.winner);
                c.terminal_first = winner->first;
                c.terminal_second = winner->second;
                c.terminal_total = verdict.winner_total;
                if (rounds[r].index >= options.repair_rounds) {
                    c.stop = "round_cap";
                    return c;
                }
                visited.insert(visited.end(), c.moves.begin(), c.moves.end());
                std::vector<TieItem> own;
                const Leg leg = run_leg(k, winner, nullptr, 0.0, rounds[r].index + 1, visited, own);
                c.leg = true;
                if (!leg.reproduced) throw std::runtime_error("fork_start_not_reproduced");
                endpoint_items.insert(endpoint_items.end(), own.begin(), own.end());
                c.items.insert(c.items.end(), own.begin(), own.end());
                c.moves.insert(c.moves.end(), leg.moves.begin(), leg.moves.end());
                c.rounds += leg.rounds;
                c.stop = leg.stop;
                c.in_top = leg.converged && leg.stop == "rescore_incomplete";
                c.last_round = leg.last_round;
                c.terminal_first = leg.terminal_first;
                c.terminal_second = leg.terminal_second;
                c.terminal_total = leg.terminal_total;
                return c;
            };
            std::vector<std::size_t> branch_forks_of(view_count, 0), legs_of(view_count, 0);
            for (std::size_t k = 1; k < view_count; ++k) {
                try {
                    if (options.test_library_sensitivity_fail_view == k) {
                        throw std::runtime_error("test hook: endpoint failed before any work");
                    }
                    // (1) the anchor at this view
                    const std::unique_ptr<SetRanking> ranked = rank_view(k);
                    std::size_t a = npos, a_rank = 0;
                    for (std::size_t p = 0; p < ranked->ranking.candidates.size() && a == npos; ++p) {
                        for (const std::size_t member :
                             ranked->sequence_candidates[ranked->ranking.candidates[p].candidate_id].members) {
                            if (member >= panel_count) continue;
                            if (a == npos || ranked->path_prior[member] > ranked->path_prior[a]) a = member;
                        }
                        if (a != npos) a_rank = p;
                    }
                    if (a == npos) throw std::runtime_error("no_panel_anchor");
                    const MosaicDiplotype anchor_k = ranked->candidates[a].diplotype;
                    const double anchor_total = ranked->ranking.candidates[a_rank].score;
                    views[k].anchor = identity_text(anchor_k.first.alleles, anchor_k.second.alleles);
                    const DiplotypeIdentity anchor_key =
                        diplotype_identity(anchor_k.first.alleles, anchor_k.second.alleles);
                    endpoint_sets.push_back(std::make_unique<SetRanking>(*ranked));
                    set_of_view[k] = endpoint_sets.back().get();

                    // (2) the trunk: shared with the call while every decision agrees, then its own
                    std::vector<TieItem> trunk;                        // at view k
                    std::vector<DiplotypeIdentity> history{anchor_key};   // its incumbents: the anchor, every move
                    std::string stop;
                    bool in_top = false;
                    std::size_t last_round = 0;
                    const auto take_leg = [&](const Leg& leg, std::vector<TieItem>& own) {
                        if (!leg.reproduced) throw std::runtime_error("start_not_reproduced");
                        for (const TieItem& item : own) trunk.push_back(item);
                        endpoint_items.insert(endpoint_items.end(), own.begin(), own.end());
                        history.insert(history.end(), leg.moves.begin(), leg.moves.end());
                        stop = leg.stop;
                        in_top = leg.converged && leg.stop == "rescore_incomplete";
                        last_round = leg.last_round;
                    };
                    const std::vector<RepairRound>& rounds = search->rounds();
                    if (anchor_key != central_anchor) {
                        fork_of[k] = "anchor";
                        ++legs_of[k];
                        std::vector<TieItem> own;
                        take_leg(run_leg(k, nullptr, &anchor_k, anchor_total, 1, {}, own), own);
                    } else {
                        const Continued c = continue_at(k, rounds, production_items, search->stop_reason(),
                                                        anchor_diplotype.first.alleles, anchor_diplotype.second.alleles,
                                                        search->incumbent_first(), search->incumbent_second(),
                                                        search->incumbent_total(), history);
                        if (c.forked) fork_of[k] = "round" + std::to_string(c.fork_round);
                        if (c.leg) ++legs_of[k];
                        trunk = c.items;
                        history.insert(history.end(), c.moves.begin(), c.moves.end());
                        stop = c.stop;
                        in_top = c.in_top;
                        last_round = c.last_round;
                    }
                    trunk_stop[k] = stop;
                    own_items[k] = trunk;
                    const bool converged = stop == "no_improvement" || stop == "no_unvisited_neighbour" ||
                                           (stop == "rescore_incomplete" && in_top);
                    const bool improving_tie = !converged && (stop == "rescore_incomplete" || stop == "tie_at_best");
                    if (!converged && !improving_tie) throw std::runtime_error("trunk_stopped_at_" + stop);
                    if (improving_tie && !options.repair_tie_branching) throw std::runtime_error("improving_move_tie");

                    // (3) exact tie branching at this view: a central branch is reused when it started identically
                    // and decided every round identically here; any other branch runs at this view
                    if (improving_tie) {
                        const TieSet root = terminal_tie(trunk, last_round, fragments.size(), band);
                        const auto run_branch_k = [&](const BranchStart& branch) {
                            const auto found =
                                central_branches.find({diplotype_identity(branch.start.first, branch.start.second), branch.history});
                            if (found != central_branches.end()) {
                                const CentralBranch& recorded = found->second;
                                std::vector<DiplotypeIdentity> visited = branch.history;
                                visited.push_back(diplotype_identity(branch.start.first, branch.start.second));
                                const Continued c = continue_at(
                                    k, recorded.rounds, recorded.result.items, recorded.result.stop, branch.start.first,
                                    branch.start.second, recorded.result.terminal_first, recorded.result.terminal_second,
                                    recorded.terminal_totals[k], visited);
                                ++(c.forked ? branch_forks_of[k] : reused_of[k]);
                                if (c.leg) ++legs_of[k];
                                BranchResult out = recorded.result;
                                out.items = c.items;
                                out.rounds = c.rounds;
                                out.moves = c.moves.size();
                                out.moved_to = c.moves;
                                out.stop = c.stop;
                                out.last_round = c.last_round;
                                out.terminal_first = c.terminal_first;
                                out.terminal_second = c.terminal_second;
                                out.terminal_total = c.terminal_total;
                                out.why.clear();
                                const bool converged = c.stop == "no_improvement" || c.stop == "no_unvisited_neighbour" ||
                                                       (c.stop == "rescore_incomplete" && c.in_top);
                                if (converged) {
                                    out.outcome = BranchOutcome::Converged;
                                } else if (c.stop == "rescore_incomplete" || c.stop == "tie_at_best") {
                                    out.outcome = BranchOutcome::Split;
                                } else {
                                    out.outcome = BranchOutcome::Incomplete;
                                    out.why = "stopped at " + c.stop;
                                }
                                own_items[k].insert(own_items[k].end(), out.items.begin(), out.items.end());
                                return out;
                            }
                            ++run_of[k];
                            ++legs_of[k];
                            BranchResult out;
                            std::vector<TieItem> own;
                            const Leg leg = run_leg(k, &branch.start, nullptr, 0.0, 1, branch.history, own);
                            endpoint_items.insert(endpoint_items.end(), own.begin(), own.end());
                            own_items[k].insert(own_items[k].end(), own.begin(), own.end());
                            out.items = own;
                            out.rounds = leg.rounds;
                            out.moves = leg.moves.size();
                            out.moved_to = leg.moves;
                            out.stop = leg.stop;
                            out.last_round = leg.last_round;
                            out.terminal_first = leg.terminal_first;
                            out.terminal_second = leg.terminal_second;
                            out.terminal_total = leg.terminal_total;
                            if (!leg.reproduced) {
                                out.outcome = BranchOutcome::Incomplete;
                                out.why = "its start did not reproduce";
                            } else if (leg.converged) {
                                out.outcome = BranchOutcome::Converged;
                            } else if (leg.stop == "rescore_incomplete" || leg.stop == "tie_at_best") {
                                out.outcome = BranchOutcome::Split;
                            } else {
                                out.outcome = BranchOutcome::Incomplete;
                                out.why = "stopped at " + leg.stop;
                            }
                            return out;
                        };
                        const TieBranchingReport branched_k = run_tie_branching(
                            root, history, options.repair_tie_max_branches, fragments.size(), band, run_branch_k,
                            [](const TieSet&) {}, false, [&](long) { shared->release_transient({}); });
                        if (!branched_k.incomplete.empty()) throw std::runtime_error("tie_branching_incomplete: " + branched_k.incomplete);
                    }
                } catch (const std::exception& error) {
                    views[k].incomplete = error.what();
                    if (views[k].incomplete.empty()) views[k].incomplete = "failed";
                }
            }

            // (4) the common pool: the call's pool, every union set an endpoint ranked, everything an endpoint scored
            struct AuditEntry {
                PoolCandidate candidate;
                std::vector<DiplotypeScore> views;
            };
            std::map<DiplotypeIdentity, AuditEntry> common;
            std::uint64_t audit_mismatches = view_mismatches;
            const auto admit = [&](const PoolCandidate& candidate, const std::vector<DiplotypeScore>& all) {
                if (all.size() != view_count) {
                    ++audit_mismatches;
                    return;
                }
                const DiplotypeIdentity key = diplotype_identity(candidate.first, candidate.second);
                const auto found = common.find(key);
                if (found == common.end()) {
                    common.emplace(key, AuditEntry{candidate, all});
                    return;
                }
                bool same = bitwise_equal(found->second.candidate.sequence_prior, candidate.sequence_prior);
                for (std::size_t v = 0; v < view_count; ++v) same = same && same_score(found->second.views[v], all[v]);
                if (!same) ++audit_mismatches;
            };
            for (const PoolEntry& entry : pool) {
                const auto found = pool_views.find(diplotype_identity(entry.candidate.first, entry.candidate.second));
                admit(entry.candidate, found == pool_views.end() ? std::vector<DiplotypeScore>{} : found->second);
            }
            const std::size_t central_pool = common.size();
            // A union haplotype's digest, computed at most once: a set holds thousands of candidates over far fewer
            // haplotypes of ~200 kb each.
            std::vector<std::optional<SequenceDigest>> union_digest(interned.sequences.size());
            const auto digest_at = [&](std::size_t h) {
                if (!union_digest[h]) union_digest[h] = sequence_digest(interned.sequences[h]);
                return *union_digest[h];
            };
            const auto union_candidate = [&](const SetRanking& set, std::size_t i) {
                const MosaicCandidate& c = set.candidates[i];
                const std::size_t x = interned.index_of[slot_of_path.at(c.diplotype.first.alleles)];
                const std::size_t y = interned.index_of[slot_of_path.at(c.diplotype.second.alleles)];
                const auto key = sequence_pair(c);
                PoolCandidate candidate;
                candidate.first = c.diplotype.first.alleles;
                candidate.second = c.diplotype.second.alleles;
                candidate.first_sequence = digest_at(x);
                candidate.second_sequence = digest_at(y);
                candidate.score = scored.scores[scored_index.at(key)];
                candidate.sequence_prior = haploid_sequence_prior[key.first] + haploid_sequence_prior[key.second] +
                                           (key.first != key.second ? std::log(2.0) : 0.0);
                candidate.path_prior = set.path_prior[i];
                return std::make_pair(candidate, union_views(scored_index.at(key)));
            };
            const auto item_candidate = [&](const TieItem& item) {
                PoolCandidate candidate;
                candidate.first = item.first;
                candidate.second = item.second;
                candidate.first_sequence = item.first_sequence;
                candidate.second_sequence = item.second_sequence;
                candidate.score = item.views.empty() ? item.score : item.views.front();
                candidate.sequence_prior = item.prior;
                const double pa = cached_haploid_prior(MosaicPath{item.first});
                const double pb = cached_haploid_prior(MosaicPath{item.second});
                const double negative_infinity = -std::numeric_limits<double>::infinity();
                candidate.path_prior = pa == negative_infinity || pb == negative_infinity
                                           ? negative_infinity
                                           : pa + pb + (item.first != item.second ? std::log(2.0) : 0.0);
                return candidate;
            };
            for (const std::unique_ptr<SetRanking>& set : endpoint_sets) {
                for (std::size_t i = 0; i < set->candidates.size(); ++i) {
                    const auto made = union_candidate(*set, i);
                    admit(made.first, made.second);
                }
            }
            // (4b) each endpoint's independent pipeline call: its own pool -- its union set, its trunk and every branch it
            // used or ran -- decided at its own view. This is what a run declaring that intensity reports.
            std::vector<PooledDecision> independents(view_count);
            std::vector<std::string> own_top_score(view_count, "NA"), common_top_score(view_count, "NA"),
                called_score(view_count, "NA");
            for (std::size_t k = 1; k < view_count; ++k) {
                if (!views[k].incomplete.empty() || set_of_view[k] == nullptr) continue;
                std::map<DiplotypeIdentity, PoolCandidate> own;
                for (std::size_t i = 0; i < set_of_view[k]->candidates.size(); ++i) {
                    const auto made = union_candidate(*set_of_view[k], i);
                    PoolCandidate candidate = made.first;
                    candidate.score = made.second.at(k);
                    own.emplace(diplotype_identity(candidate.first, candidate.second), candidate);
                }
                for (const TieItem& item : own_items[k]) {
                    PoolCandidate candidate = item_candidate(item);
                    candidate.score = item.views.at(k);
                    own.emplace(diplotype_identity(candidate.first, candidate.second), candidate);
                }
                std::vector<PoolCandidate> at;
                for (const auto& [key, candidate] : own) at.push_back(candidate);
                independents[k] = decide_pooled_genotype(std::move(at), block_count, fragments.size(), band);
                if (independents[k].score_mismatches != 0) {
                    views[k].incomplete = "independent_sequence_score_mismatch";
                    continue;
                }
                views[k].independent = &independents[k];
                std::ostringstream number;
                number << std::setprecision(17) << independents[k].ranking.candidates.front().score;
                own_top_score[k] = number.str();
            }
            for (const TieItem& item : endpoint_items) admit(item_candidate(item), item.views);

            // (5) the common pool reranked at each view, and the verdict. The reported call is never rescored from this
            // pool: endpoint candidates only describe the endpoints.
            std::vector<PooledDecision> decisions(view_count);
            const auto called_key = common.find(diplotype_identity(called.first, called.second));
            for (std::size_t v = 0; v < view_count; ++v) {
                if (audit_mismatches != 0) {
                    views[v].incomplete = "pool_score_mismatch";
                    continue;
                }
                if (!views[v].incomplete.empty()) continue;
                std::vector<PoolCandidate> at;
                at.reserve(common.size());
                for (const auto& [key, entry] : common) {
                    at.push_back(entry.candidate);
                    at.back().score = entry.views[v];
                }
                decisions[v] = decide_pooled_genotype(std::move(at), block_count, fragments.size(), band);
                if (decisions[v].score_mismatches != 0) {
                    views[v].incomplete = "sequence_score_mismatch";
                    continue;
                }
                views[v].decision = &decisions[v];
                std::ostringstream number;
                number << std::setprecision(17) << decisions[v].ranking.candidates.front().score;
                common_top_score[v] = number.str();
                if (called_key != common.end()) {
                    number.str("");
                    number << called_key->second.views[v].log_likelihood + called_key->second.candidate.sequence_prior;
                    called_score[v] = number.str();
                }
            }
            // views[0] is the call's own view; [1] and [2] are the endpoints, reported as low and high
            const LibrarySensitivity verdict = classify_library_sensitivity(decided.genotype, views);
            if (verdict.status == "SENSITIVE") library_sensitive.insert(verdict.sensitive_blocks.begin(), verdict.sensitive_blocks.end());

            const auto pairs_text = [](const std::set<std::pair<std::size_t, std::size_t>>& pairs) {
                std::string set;
                for (const auto& pair : pairs) set += (set.empty() ? "" : "+") + std::to_string(pair.first) + "/" + std::to_string(pair.second);
                return set.empty() ? std::string("NA") : set;
            };
            const auto ids_text = [&](const std::vector<std::size_t>& blocks_of) {
                std::string ids;
                for (const std::size_t b : blocks_of) ids += (ids.empty() ? "" : ",") + model.blocks()[b].block_id;
                return ids.empty() ? std::string("none") : ids;
            };
            const FragmentScoringCounters& work_after = shared->cumulative_counters();
            std::size_t diverged = 0, legs = 0;
            for (std::size_t v = 1; v < view_count; ++v) {
                diverged += (fork_of[v] != "none" ? 1 : 0) + branch_forks_of[v] + run_of[v];
                legs += legs_of[v];
            }
            std::ostringstream text;
            text << std::setprecision(17);
            std::string incomplete_endpoint = "none", incomplete_reason = verdict.incomplete;
            for (std::size_t v = 0; v < view_count && incomplete_endpoint == "none"; ++v) {
                if (!views[v].incomplete.empty()) {
                    incomplete_endpoint = view_name[v];
                    incomplete_reason = views[v].incomplete;
                }
            }
            library_not_evaluated = verdict.status == "NOT_EVALUATED";
            text << "library_sensitivity\t" << verdict.status << "\nlibrary_sensitivity_reason\t" << verdict.reason
                 << "\nlibrary_sensitivity_status\t" << verdict.status << "\nlibrary_sensitivity_filter\t"
                 << (verdict.status == "STABLE"      ? "PASS"
                     : verdict.status == "SENSITIVE" ? "LOWQUAL_LIBRARY_SENSITIVE"
                                                     : "NOT_EVALUATED_LIBRARY_SENSITIVITY")
                 << "\nlibrary_sensitivity_incomplete_endpoint\t" << incomplete_endpoint
                 << "\nlibrary_sensitivity_incomplete\t" << incomplete_reason << "\nlibrary_sensitivity_rule\t"
                 << kLibrarySensitivityRule << "\nlibrary_sensitivity_lambdas\t"
                 << scoring_options.audit_libraries[0].fragment_start_intensity << ','
                 << library.fragment_start_intensity << ',' << scoring_options.audit_libraries[1].fragment_start_intensity
                 << "\nlibrary_sensitive_blocks\t" << ids_text(verdict.sensitive_blocks)
                 << "\nlibrary_sensitivity_rank_blocks\t" << ids_text(verdict.rank_blocks)
                 << "\nlibrary_sensitivity_pipeline_blocks\t" << ids_text(verdict.pipeline_blocks);
            std::string reasons;
            for (const std::size_t b : verdict.sensitive_blocks) {
                reasons += (reasons.empty() ? "" : ",") + model.blocks()[b].block_id + ":" + verdict.block_reason[b];
            }
            text << "\nlibrary_sensitive_block_reasons\t" << (reasons.empty() ? "none" : reasons)
                 << "\nlibrary_sensitivity_phase_or_representative_only\t"
                 << (verdict.phase_or_representative_only ? "true" : "false")
                 << "\nlibrary_sensitivity_central_matches_primary\t" << (verdict.central_matches_primary ? "true" : "false")
                 << "\nlibrary_sensitivity_distinct_anchors\t" << verdict.distinct_anchors
                 << "\nlibrary_sensitivity_trajectory_forked\t" << (diverged > 0 ? "true" : "false")
                 << "\nlibrary_sensitivity_divergences\t" << diverged << "\nlibrary_sensitivity_fork_legs\t" << legs
                 << '\n';
            for (const std::size_t v : {std::size_t{1}, std::size_t{0}, std::size_t{2}}) {
                const std::string name = std::string("library_sensitivity_") + view_name[v];
                std::string rank_call, pipeline_call;
                for (std::size_t b = 0; b < block_count && !verdict.projections[v].empty(); ++b) {
                    rank_call += (rank_call.empty() ? "" : ",") + pairs_text(verdict.projections[v][b]);
                }
                for (std::size_t b = 0; b < block_count && !verdict.independent_projections[v].empty(); ++b) {
                    pipeline_call += (pipeline_call.empty() ? "" : ",") + pairs_text(verdict.independent_projections[v][b]);
                }
                text << name << "_lambda\t"
                     << (v == 0 ? library.fragment_start_intensity
                                : scoring_options.audit_libraries[v - 1].fragment_start_intensity)
                     << '\n'
                     << name << "_anchor\t" << views[v].anchor << '\n'
                     << name << "_fork\t" << fork_of[v] << '\n'
                     << name << "_trunk_stop\t" << trunk_stop[v] << '\n'
                     << name << "_incomplete\t" << (views[v].incomplete.empty() ? "none" : views[v].incomplete) << '\n'
                     << name << "_rank_call_class\t" << verdict.call_class[v] << '\n'
                     << name << "_rank_call\t" << (rank_call.empty() ? "NA" : rank_call) << '\n'
                     << name << "_pipeline_call_class\t" << verdict.independent_call_class[v] << '\n'
                     << name << "_pipeline_call\t" << (pipeline_call.empty() ? "NA" : pipeline_call) << '\n'
                     << name << "_top_set_size\t" << verdict.top_set_size[v] << '\n'
                     << name << "_top_set_digest\t" << verdict.top_set_digest[v] << '\n'
                     << name << "_representative_digest\t" << verdict.representative_digest[v] << '\n'
                     << name << "_pipeline_digest\t" << verdict.independent_digest[v] << '\n'
                     << name << "_reported_call_score\t" << called_score[v] << '\n'
                     << name << "_rank_top_score\t" << common_top_score[v] << '\n'
                     << name << "_pipeline_top_score\t" << own_top_score[v] << '\n'
                     << name << "_branches_reused\t" << reused_of[v] << '\n'
                     << name << "_branches_forked\t" << branch_forks_of[v] << '\n'
                     << name << "_branches_run\t" << run_of[v] << '\n'
                     << name << "_fork_legs\t" << legs_of[v] << '\n';
            }
            text << "library_sensitivity_pool_identities\t" << common.size() << "\nlibrary_sensitivity_pool_added\t"
                 << common.size() - central_pool << "\nlibrary_sensitivity_added_haplotypes_placed\t"
                 << shared->haplotypes_created() - created_before << "\nlibrary_sensitivity_added_physical_origins\t"
                 << work_after.physical_origins - work_before.physical_origins
                 << "\nlibrary_sensitivity_added_state_classes\t"
                 << work_after.statistical_state_classes - work_before.statistical_state_classes
                 << "\nlibrary_sensitivity_added_executed_dps\t" << read_emission_cache_counters().executed_dps - dps_before
                 << "\nseconds_library_sensitivity\t" << std::fixed << std::setprecision(3) << seconds_since(audit_mark)
                 << '\n';
            report.library_sensitivity_text = text.str();
            // The per-block table, whole or not at all, and only for a completed verdict.
            if (verdict.status != "NOT_EVALUATED") {
                std::ostringstream table;
                table << "block_id\treported\tlow_rank\tlow_pipeline\thigh_rank\thigh_pipeline\tlikelihood_rank_sensitive\t"
                         "pipeline_sensitive\treason\n";
                for (std::size_t b = 0; b < block_count; ++b) {
                    const auto pair = decided.genotype[b];
                    const bool rank = std::find(verdict.rank_blocks.begin(), verdict.rank_blocks.end(), b) != verdict.rank_blocks.end();
                    const bool pipeline =
                        std::find(verdict.pipeline_blocks.begin(), verdict.pipeline_blocks.end(), b) != verdict.pipeline_blocks.end();
                    table << model.blocks()[b].block_id << '\t' << pair.first << '/' << pair.second << '\t'
                          << pairs_text(verdict.projections[1][b]) << '\t' << pairs_text(verdict.independent_projections[1][b])
                          << '\t' << pairs_text(verdict.projections[2][b]) << '\t'
                          << pairs_text(verdict.independent_projections[2][b]) << '\t' << (rank ? "true" : "false") << '\t'
                          << (pipeline ? "true" : "false") << '\t' << verdict.block_reason[b] << '\n';
                }
                defer_output(prefix + ".library_sensitivity_blocks.tsv") << table.str();
            }
            if (!library_sensitive.empty()) report.locus_filter = with_library_flag(report.locus_filter, "LOWQUAL_LIBRARY_SENSITIVE");
            note_footprint("library_sensitivity_done");   // S2-K (only when the audit ran)
        }
        // An incomplete audit is a locus-level warning: no block is accused, block quality is unchanged.
        if (library_not_evaluated) report.locus_filter = with_library_flag(report.locus_filter, "NOT_EVALUATED_LIBRARY_SENSITIVITY");
        const auto filter_columns = [&](std::size_t b) {
            const BlockFilter& bf = block_filter[b];
            const auto number = [](double v) {
                if (std::isnan(v)) return std::string("NA");
                std::ostringstream text;
                text << std::setprecision(12) << v;
                return text.str();
            };
            return "\t" + bf.filter + "\t" + bf.reason + "\t" +
                   (bf.filter == "NOT_EVALUATED" && bf.reason == "guard_off" ? std::string("NA") : number(bf.marker_margin)) +
                   "\t" + bf.competitor + "\t" + number(bf.fragment_diff) + "\t" + number(bf.exposure_diff) + "\t" +
                   number(bf.likelihood_margin) + "\t" + number(bf.statistic);
        };
        // Allele order: the representative's haplotypes for a unique call; the unordered pair otherwise.
        const auto alleles_at = [&](std::size_t b) {
            return unique ? std::make_pair(called.first[b], called.second[b]) : decided.genotype[b];
        };
        const auto sources_of = [&](const PreparedGenotypeBlock& block, std::size_t id) {
            std::string joined;
            for (const std::string& source : block.alleles[id].source_ids) joined += (joined.empty() ? "" : ",") + source;
            return joined;
        };
        // The combined block quality of every block, recorded once in the summary; both tables carry it.
        std::vector<std::string> genotype_quality(block_count);
        for (std::size_t b = 0; b < block_count; ++b) {
            genotype_quality[b] = quality_columns(summary, model.blocks()[b].block_id, block_filter[b].filter,
                                                  summary.block_genotype_read_supported[b] ? "READS" : "PRIOR_ONLY",
                                                  library_sensitive.count(b) != 0);
        }
        {
            std::ostringstream& out = defer_output(prefix + ".genotypes.tsv");
            out << "block_id\tkind\tallele1\tallele2\tallele1_bp\tallele2_bp\tallele1_sources\t"
                   "allele2_sources\tzygosity\tconditional_genotype_posterior\tconditional_GQ" << kFilterHeader
                << kQualityHeader << '\n';
            for (std::size_t b = 0; b < block_count; ++b) {
                const PreparedGenotypeBlock& block = model.blocks()[b];
                if (block.kind != GenotypeBlockKind::Variant) continue;
                const auto [a1, a2] = alleles_at(b);
                out << block.block_id << "\tSV\t" << a1 << '\t' << a2 << '\t' << block.alleles[a1].sequence.size() << '\t'
                    << block.alleles[a2].sequence.size() << '\t' << sources_of(block, a1) << '\t' << sources_of(block, a2)
                    << '\t' << (a1 == a2 ? "HOM" : "HET") << '\t' << std::setprecision(12)
                    << summary.block_genotype_posterior[b] << '\t' << summary.block_genotype_quality[b]
                    << std::setprecision(6) << filter_columns(b) << genotype_quality[b] << '\n';
            }
        }
        {
            std::ostringstream& out = defer_output(prefix + ".block_genotypes.tsv");
            out << "block_id\tkind\treportable\tallele1\tallele2\tallele1_bp\tallele2_bp\tzygosity\t"
                   "projection_ambiguous\tgenotype_support\tphase_support\tconditional_genotype_posterior\t"
                   "conditional_GQ\tnumerically_separated\tcandidate_scope" << kFilterHeader << kQualityHeader
                << data_support_header(options.block_support_rule) << '\n';
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
                    << set_report.final_candidate_scope << filter_columns(b) << genotype_quality[b]
                    << (options.block_support_rule != "legacy"
                            ? data_support_columns(options.block_support_rule, summary, b, block_filter[b].filter, library_sensitive.count(b) != 0, true)
                            : std::string())
                    << '\n';
            }
        }
        if (unique) {
            // A unique candidate's two haplotypes. Under PRIOR_ONLY the reads gave no preference among the
            // phase alternatives and the prior chose this one: a prior-selected haplotype representation,
            // not a read-phased result, and the header says so.
            report.haplotype_representation =
                summary.phase_support == "PRIOR_ONLY" ? "prior_selected" : "read_supported_within_candidates";
            std::ostringstream& out = defer_output(prefix + ".mosaic_haplotypes.fa");
            const HaplotypeSequence first = spell_mosaic(model, MosaicPath{called.first});
            const HaplotypeSequence second = spell_mosaic(model, MosaicPath{called.second});
            for (const auto& item : {std::make_pair("hap1", &first), std::make_pair("hap2", &second)}) {
                out << '>' << prefix << '_' << item.first << " fingerprint=" << fingerprint
                    << " bp=" << item.second->bases.size() << " phase_support=" << summary.phase_support
                    << " representation=" << report.haplotype_representation << " locus_filter=" << report.locus_filter
                    << '\n';
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
            if (g_incumbent.enabled && g_incumbent.have && !g_incumbent.failed) {   // S2-A: a declared cap stop
                g_sequence_report = checkpoint_report("union_scoring_cap", "the union pass stopped at a declared "
                                                      "resource cap (" + rescore_failure + "); the best fully evaluated "
                                                      "candidate so far; search incomplete");
            }
            write_status(status_path, GenotypeRunStatus::RescoreIncomplete, fingerprint, options, library,
                         timings, &locus_report, &panel, (panel_only ? nullptr : &beam), nullptr,
                         "rescoring stopped at a declared resource cap and returned no ranking: " +
                             rescore_failure,
                         nullptr, &depth, &candidate_report);
            return GenotypeRunStatus::RescoreIncomplete;
        }

        if (options.diagnostic_union_only &&
            (options.diagnostic_union_seeding_profile > 0 || !options.diagnostic_union_verification_census.empty())) {
            timings.total = seconds_since(started);
            write_status(status_path, GenotypeRunStatus::UnionDiagnosticOnly, fingerprint, options, library, timings,
                         &locus_report, &panel, (panel_only ? nullptr : &beam), nullptr,
                         "DIAGNOSTIC ONLY: union seeding profile or census; the union pass was replayed, not scored",
                         nullptr, &depth, &candidate_report);
            return GenotypeRunStatus::UnionDiagnosticOnly;
        }
        if (options.diagnostic_union_only) {
            // DIAGNOSTIC ONLY: the union pass's certified ranking, exact to the bit, and a digest of every mass
            // cell and exposure; the run stops here -- no repair, no genotype, no haplotype, never a call.
            SetRanking diagnostic = rank_set(beam);
            std::uint64_t digest = 1469598103934665603ULL;
            const auto mix = [&digest](const double value) {
                unsigned char bytes[sizeof(double)];
                std::memcpy(bytes, &value, sizeof(double));
                for (const unsigned char b : bytes) digest = (digest ^ b) * 1099511628211ULL;
            };
            for (const double value : scored.log_mass) mix(value);
            for (const double value : scored.exposure) mix(value);
            std::ostringstream text;
            text << "# DIAGNOSTIC ONLY: union-pass ranking; not a genotype call\n"
                 << "# ranking_status\t" << (diagnostic.ranking.status == RescoreStatus::Complete ? "COMPLETE" : "TIED")
                 << "\n# winner_set_size\t" << diagnostic.ranking.winner_set.size() << "\n# mass_cells\t"
                 << scored.log_mass.size() << "\n# union_mass_digest\t" << std::hex << digest << std::dec
                 << "\nrank\thap1_sequence\thap2_sequence\tscore\tfragment_term\texposure_term\tprior\t"
                    "numeric_error_bound\tfragments_with_mass\n";
            const auto exact = [](const double value) {
                char buffer[64];
                std::snprintf(buffer, sizeof(buffer), "%a", value);
                return std::string(buffer);
            };
            for (std::size_t i = 0; i < diagnostic.ranking.candidates.size(); ++i) {
                const CertifiedCandidate& entry = diagnostic.ranking.candidates[i];
                const SequenceCandidate& pair = diagnostic.sequence_candidates[entry.candidate_id];
                text << i << '\t' << pair.first << '\t' << pair.second << '\t' << exact(entry.score) << '\t'
                     << exact(entry.fragment_term) << '\t' << exact(entry.exposure_term) << '\t' << exact(entry.prior)
                     << '\t' << exact(entry.numeric_error_bound) << '\t' << entry.fragments_with_mass << '\n';
            }
            atomic_write(options.out_prefix + ".union_diagnostic_ranking.tsv", text.str());
            timings.total = seconds_since(started);
            write_status(status_path, GenotypeRunStatus::UnionDiagnosticOnly, fingerprint, options, library, timings,
                         &locus_report, &panel, (panel_only ? nullptr : &beam), &diagnostic.ranking,
                         "DIAGNOSTIC ONLY: the union pass was scored and ranked; no repair and no genotype",
                         nullptr, &depth, &candidate_report);
            return GenotypeRunStatus::UnionDiagnosticOnly;
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
                // The combined block quality: this path does not run the representation guard (NOT_EVALUATED), and
                // a no-called block has no read support.
                std::vector<std::string> consensus_quality(block_count);
                for (std::size_t b = 0; b < block_count; ++b) {
                    consensus_quality[b] = quality_columns(
                        summary, model.blocks()[b].block_id, "NOT_EVALUATED",
                        !agreed[b] ? "NOCALL" : summary.block_genotype_read_supported[b] ? "READS" : "PRIOR_ONLY");
                }
                {
                    std::ofstream out = open_output(options.out_prefix + ".genotypes.tsv");
                    out << "block_id\tkind\tallele1\tallele2\tallele1_bp\tallele2_bp\tallele1_sources\t"
                           "allele2_sources\tzygosity\tconditional_genotype_posterior\tconditional_GQ" << kFilterHeader
                        << kQualityHeader << '\n';
                    for (std::size_t b = 0; b < block_count; ++b) {
                        const PreparedGenotypeBlock& block = model.blocks()[b];
                        if (block.kind != GenotypeBlockKind::Variant) continue;
                        if (!agreed[b]) {
                            // Emitted as an explicit no-call rather than omitted: a missing row is
                            // indistinguishable from a block that was never modelled.
                            out << block.block_id << "\tSV\t.\t.\t.\t.\t.\t.\tNOCALL\tNA\tNA" << kFilterNotEvaluated
                                << consensus_quality[b] << '\n';
                            continue;
                        }
                        const std::size_t a1 = consensus[b].first;
                        const std::size_t a2 = consensus[b].second;
                        out << block.block_id << "\tSV\t" << a1 << '\t' << a2 << '\t'
                            << block.alleles[a1].sequence.size() << '\t'
                            << block.alleles[a2].sequence.size() << '\t' << sources_of(block, a1) << '\t'
                            << sources_of(block, a2) << '\t' << (a1 == a2 ? "HOM" : "HET") << '\t'
                            << std::setprecision(12) << summary.block_genotype_posterior[b] << '\t'
                            << summary.block_genotype_quality[b] << std::setprecision(6) << kFilterNotEvaluated
                            << consensus_quality[b] << '\n';
                    }
                }
                {
                    std::ofstream out = open_output(options.out_prefix + ".block_genotypes.tsv");
                    out << "block_id\tkind\treportable\tallele1\tallele2\tallele1_bp\tallele2_bp\t"
                           "zygosity\tprojection_ambiguous\tgenotype_support\tphase_support\t"
                           "conditional_genotype_posterior\tconditional_GQ\tnumerically_separated\t"
                           "candidate_scope" << kFilterHeader << kQualityHeader << data_support_header(options.block_support_rule) << '\n';
                    for (std::size_t b = 0; b < block_count; ++b) {
                        const PreparedGenotypeBlock& block = model.blocks()[b];
                        const bool variant = block.kind == GenotypeBlockKind::Variant;
                        out << block.block_id << '\t' << (variant ? "SV" : "context") << '\t'
                            << (variant ? "true" : "false") << '\t';
                        if (!agreed[b]) {
                            out << ".\t.\t.\t.\tNOCALL\tfalse\tNOCALL\tUNRESOLVED\tNA\tNA\tfalse\t"
                                << set_report.final_candidate_scope << kFilterNotEvaluated << consensus_quality[b]
                                << (options.block_support_rule != "legacy" ? data_support_columns(options.block_support_rule, summary, b, "", false, false) : std::string())
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
                            << set_report.final_candidate_scope << kFilterNotEvaluated << consensus_quality[b]
                            << (options.block_support_rule != "legacy" ? data_support_columns(options.block_support_rule, summary, b, "", false, false) : std::string())
                            << '\n';
                    }
                }
                timings.total = seconds_since(started);
                write_status(status_path, GenotypeRunStatus::BlockConsensus, fingerprint, options,
                             library, timings, &locus_report, &panel, (panel_only ? nullptr : &beam),
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
                         library, timings, &locus_report, &panel, (panel_only ? nullptr : &beam), &ranking,
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
        // The combined block quality: this path does not run the representation guard (NOT_EVALUATED).
        std::vector<std::string> plain_quality(model.blocks().size());
        for (std::size_t b = 0; b < model.blocks().size(); ++b) {
            plain_quality[b] = quality_columns(summary, model.blocks()[b].block_id, "NOT_EVALUATED",
                                               summary.block_genotype_read_supported[b] ? "READS" : "PRIOR_ONLY");
        }
        {
            // Only variant blocks are reporting units. Context blocks are modelled identically and are
            // what let the mosaic recombine, but they are not SV calls and are not presented as such.
            std::ofstream out = open_output(options.out_prefix + ".genotypes.tsv");
            out << "block_id\tkind\tallele1\tallele2\tallele1_bp\tallele2_bp\tallele1_sources\t"
                   "allele2_sources\tzygosity\tconditional_genotype_posterior\tconditional_GQ" << kFilterHeader
                << kQualityHeader << '\n';
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
                    << std::setprecision(6) << kFilterNotEvaluated << plain_quality[b] << '\n';
            }
        }
        {
            // Every statistical block, for block-level evaluation. Context blocks are polymorphic and are
            // genotyped by the same mosaic; `reportable` marks the SV bubbles the SV-facing genotypes table
            // carries, so the interstitial blocks are distinguished rather than dropped.
            std::ofstream out = open_output(options.out_prefix + ".block_genotypes.tsv");
            out << "block_id\tkind\treportable\tallele1\tallele2\tallele1_bp\tallele2_bp\tzygosity\t"
                   "projection_ambiguous\tgenotype_support\tphase_support\tconditional_genotype_posterior\t"
                   "conditional_GQ\tnumerically_separated\tcandidate_scope" << kFilterHeader << kQualityHeader
                << data_support_header(options.block_support_rule) << '\n';
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
                    << '\t' << set_report.final_candidate_scope << kFilterNotEvaluated << plain_quality[b]
                    << (options.block_support_rule != "legacy" ? data_support_columns(options.block_support_rule, summary, b, "NOT_EVALUATED", false, true)
                                                               : std::string())
                    << '\n';
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
                     &locus_report, &panel, (panel_only ? nullptr : &beam), &ranking, "", &summary,
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

namespace {
// The advanced scientific settings exposed as named flags (panvar genotype --help-advanced). Each is exactly its profile
// key: same parser, same validation, and precedence flag > profile > preset > default whatever the argument order.
struct AdvancedFlag {
    const char* flag;
    const char* key;
    bool on_off;
    const char* help;
    bool development_only = false;   // refused under adaptive repair, the production default
};
const AdvancedFlag kAdvancedFlags[] = {
    {"--kmer-size", "kmer_size", false, "marker k-mer length, 1..31"},
    {"--syncmer-s", "syncmer_s", false, "marker syncmer s, 1..k"},
    {"--min-variant-bp", "min_variant_bp", false,
     "a candidate variant site is kept when some path spans at least N bp inside it"},
    {"--switch-probability", "switch_probability", false, "per-boundary template-switch prior, in (0, 0.5]"},
    {"--proposal-state-limit", "proposal_state_limit", false, "mosaic proposal beam width"},
    {"--mismatch-band", "mismatch_band", false,
     "per-mate edit budget of the banded gapped placement model (substitutions + gap bases)"},
    {"--repair-local-pairs", "repair_local_pairs", false, "allele pairs per block in a repair neighbourhood"},
    {"--repair-max-neighbours", "repair_max_neighbours", false, "declared bound on one repair neighbourhood (0: none)"},
    {"--repair-tie-branching", "repair_tie_branching", true, "exact branch-and-pool at improving-move ties"},
    {"--repair-tie-max-branches", "repair_tie_max_branches", false, "tie-branching cap; beyond it the run no-calls"},
    {"--representation-guard", "representation_guard", true, "the representation-conflict annotation"},
};
const AdvancedFlag* advanced_setting(const std::string& flag) {
    for (const AdvancedFlag& a : kAdvancedFlags) if (flag == a.flag) return &a;
    return nullptr;
}
} // namespace

// ---- M15-H: the public option inventory (one table; the advanced help, the guide's inventory and the tests derive from it)
namespace {
struct OptionEntry {
    const char* group;    // everyday | input | model | search | resources | reporting | legacy
    const char* option;   // the flag, or the profile key for profile-only controls
    const char* aliases;  // "" when none
    const char* key;      // the profile / settings key ("" when the option is not a setting)
    const char* values;
    const char* units;
    const char* effect;   // inference | computation | reporting | completion | input (combinations with +)
    const char* status;   // everyday | supported | expert profile only | legacy | development
    const char* text;
    const char* fixed_default;   // "" = computed from the preset and the built-in options
};
// Groups, in help order.
const std::vector<std::pair<std::string, std::string>> kOptionGroups = {
    {"everyday", "Everyday options"},
    {"input", "1. Input and panel overrides"},
    {"model", "2. Copying prior and marker / placement model"},
    {"search", "3. Candidate generation and bounded refinement"},
    {"resources", "4. Resource limits, work caps, cache budgets and exact backends"},
    {"reporting", "5. Reporting and diagnostics"},
    {"legacy", "6. Legacy and development compatibility"},
};
const OptionEntry kOptionEntries[] = {
    // everyday
    {"everyday", "--graph FILE", "--gfa, -i", "", "a GFA path", "-", "input", "everyday",
     "the locus graph; every P/W path is a panel haplotype unless curated out, and must cross the whole locus", "required"},
    {"everyday", "--reads1 FILE / --reads2 FILE", "", "", "FASTQ paths (.gz accepted)", "-", "input", "everyday",
     "paired reads, mates in the same order", "required"},
    {"everyday", "--out PREFIX", "--out-prefix, -o", "", "a path prefix", "-", "reporting", "everyday",
     "every report is written at the prefix; a run first removes the run-owned reports already there", "genotype"},
    {"everyday", "--reference-path NAME", "-r", "", "a path name", "-", "inference", "everyday",
     "the coordinate reference and orientation frame", "the single path of the graph's RS:Z sample (refused, candidates listed, when absent or ambiguous)"},
    {"everyday", "--threads N", "", "threads", "N >= 1", "worker threads", "computation", "everyday",
     "worker threads; results do not depend on them; memory is not divided per thread", ""},
    {"everyday", "--dosage-check", "", "", "flag", "-", "reporting", "everyday",
     "the optional copy-number audit of the reported alleles; never changes the prediction; adds run time; the only public dosage control (dosage is never repaired)", "off (dosage_check NOT_REQUESTED)"},
    {"everyday", "--build-panel-index FILE", "", "", "a new file path", "-", "computation", "everyday",
     "build the sample-independent panel index (structural model + marker panel) and stop; no reads, library or output prefix; never overwrites; refused beside --panel-index or a sample argument", "none"},
    {"everyday", "--panel-index FILE", "", "", "an index path", "-", "computation", "everyday",
     "genotype against a built index (saves locus preparation only); --graph is still required and read; integrity, versions and identity are verified before the output prefix is touched; never rebuilt silently", "none (prepared in memory)"},
    {"everyday", "-q, --quiet", "", "", "flag", "-", "reporting", "everyday", "no progress messages (warnings are still printed)", "off"},
    {"everyday", "-h, --help / --help-advanced", "", "", "flag", "-", "reporting", "everyday", "this help / every option by purpose", "-"},
    // 1. input and panel overrides
    {"input", "--library auto|FILE", "", "", "auto, or a library file", "-", "inference", "supported",
     "auto (also when omitted) estimates only the insert mean/sd (>= 500 decisive anchored pairs) and the fragment-start intensity (marker depth fit), and refuses (LIBRARY_REFUSED, never a fallback); support_sigmas 4, substitution_error 0.002, outlier_weight 0.001 and the gap model stay fixed. FILE declares exactly insert_mean_bp, insert_sd_bp, support_sigmas, substitution_error, outlier_weight, fragment_start_intensity (and optionally marker_background), each once; gap_open, gap_extend and the maximum displacement cannot be supplied", "auto"},
    {"input", "--panel-curation FILE", "", "", "a file of path names", "-", "inference", "supported",
     "paths removed from the active panel (they still define the blocks); part of a panel index's identity", "none"},
    {"input", "--bubbles FILE", "", "", "a bubbles CSV", "-", "inference", "supported",
     "precomputed variant sites instead of the graph's own decomposition; part of a panel index's identity", "computed from the graph"},
    {"input", "--min-variant-bp N", "", "min_variant_bp", "N >= 1", "bp", "inference", "supported",
     "SITE CONSTRUCTION, not call filtering: a candidate variant site is kept when some path spans at least N bp inside it; part of a panel index's identity", ""},
    {"input", "--preset short-read", "", "", "short-read", "-", "inference+computation", "supported",
     "the production algorithm, short-read constants and the accelerated configuration; always applied (long-read is refused, not implemented)", "short-read"},
    {"input", "--profile FILE", "", "", "key<TAB>value file naming itself", "-", "any", "supported",
     "optional expert overrides, never required: one file of literal key<TAB>value lines with a profile<TAB>name row; precedence default < preset < profile < flag; given once (a repeated --profile is refused); a key set twice is refused; output_schema may only restate sequence_first_v2 and dosage_mode is not a profile setting; not by itself a validated configuration", "none"},
    // 2. copying prior and marker/placement model
    {"model", "--switch-probability P", "", "switch_probability", "0 < P <= 0.5", "per block boundary", "inference", "supported",
     "the DONOR-REDRAW probability of the Li-Stephens copying prior: how often a haplotype is allowed to change panel donor at a block boundary; a prior regularisation, NOT a biological recombination probability", ""},
    {"model", "--kmer-size K", "", "kmer_size", "1 <= K <= 31", "bp", "inference", "supported", "marker k-mer length; part of a panel index's identity", ""},
    {"model", "--syncmer-s S", "", "syncmer_s", "1 <= S <= K", "bp", "inference", "supported", "marker syncmer s; part of a panel index's identity", ""},
    {"model", "--mismatch-band N", "", "mismatch_band", "N >= 0", "edits per mate", "inference", "supported",
     "per-mate edit budget of the banded gapped placement model (substitutions + gap bases)", ""},
    {"model", "copying_transition", "", "copying_transition", "constant_v1 | distance_ref_core_v1_experimental", "-", "inference", "expert profile only",
     "the transition law of the copying prior (distance-scaled is experimental, not in the tested configuration)", ""},
    {"model", "copying_length_bp", "", "copying_length_bp", "> 0 (hex double)", "bp", "inference", "expert profile only",
     "the length scale of the distance-scaled transition (used only by distance_ref_core_v1_experimental)", ""},
    // 3. candidate generation and bounded refinement
    {"search", "--mosaic-candidates K", "--max-mosaic-candidates", "max_mosaic_candidates", "K >= 0", "mosaics", "inference+computation", "supported",
     "the INITIAL off-panel mosaics added to every panel diplotype in the union ranking; it does not cap the mosaics the refinement generates later", ""},
    {"search", "--proposal-state-limit N", "", "proposal_state_limit", "N >= 1", "beam states", "inference", "supported", "the mosaic proposal beam width", ""},
    {"search", "--adaptive-repair on|off", "", "", "on | off", "-", "inference", "supported",
     "BOUNDED coordinate refinement: local moves scored exactly within declared caps; not a globally exhaustive search. off is the unrepaired development view (legacy report)", "on"},
    {"search", "--repair-rounds N", "", "repair_rounds", "N >= 0", "rounds", "inference+completion", "supported",
     "the refinement's round cap; 0 is the unrepaired view; resolved with --adaptive-repair and --repair-tie-branching, an explicit contradiction refused", ""},
    {"search", "--repair-local-pairs N", "", "repair_local_pairs", "N >= 1", "allele pairs per block", "inference", "supported", "allele pairs per block in one refinement neighbourhood", ""},
    {"search", "--repair-max-neighbours model|N", "", "repair_max_neighbours", "model, or N >= 0 (0: none)", "candidates per round", "completion", "supported",
     "model: the bound is the model's own structural maximum, sum over blocks of 2 * min(repair_local_pairs, a(a+1)/2), recorded as repair_neighbour_bound; N: a declared cap, a larger neighbourhood stops (neighbourhood_bound: PROVISIONAL with a checkpoint, else no prediction)", ""},
    {"search", "--repair-tie-branching on|off", "", "repair_tie_branching", "on | off", "-", "inference", "supported", "exact branch-and-pool at improving-move ties", ""},
    {"search", "--repair-tie-max-branches N", "", "repair_tie_max_branches", "N >= 1", "branches", "completion", "supported", "the tie-branching cap; beyond it the run does not complete (tie_branching_incomplete)", ""},
    {"search", "repair_junction_proposals", "", "repair_junction_proposals", "off | experimental_v1", "-", "inference", "expert profile only",
     "junction proposals during refinement (JP1.1, candidate generation): can change predictions; needs repair rounds > 0", ""},
    {"search", "early_incumbent", "", "early_incumbent", "off | experimental_v1", "-", "completion", "expert profile only",
     "keep the best fully evaluated candidate as a PROVISIONAL checkpoint on a budget stop; completed runs unchanged; needs sequence_first_v2 and adaptive repair", ""},
    {"search", "repair_joint_search", "", "repair_joint_search", "off | experimental_v1", "-", "inference", "expert profile only", "experimental joint two-block moves (not in the tested configuration)", ""},
    // 4. resources
    {"resources", "footprint_limit_bytes", "", "footprint_limit_bytes", "N >= 0 (0: none)", "bytes", "completion", "expert profile only",
     "total-footprint guard, SAMPLED every 100 ms (not an allocation ceiling); a stop exits 125, or 3 PROVISIONAL with a checkpoint", ""},
    {"resources", "time_limit_seconds", "", "time_limit_seconds", "S >= 0 (0: none)", "seconds", "completion", "expert profile only",
     "wall-clock guard; a stop exits 124, or 3 PROVISIONAL with a checkpoint", ""},
    {"resources", "max_physical_origins", "", "max_physical_origins", "N >= 0", "origins", "completion", "expert profile only", "union-pass work cap (typed refusal when reached)", ""},
    {"resources", "max_state_classes", "", "max_state_classes", "N >= 0", "state classes", "completion", "expert profile only", "union-pass work cap", ""},
    {"resources", "max_piece_hits", "", "max_piece_hits", "N >= 0", "piece hits", "completion", "expert profile only", "union-pass work cap", ""},
    {"resources", "max_mass_cells", "", "max_mass_cells", "N >= 0", "cells", "completion", "expert profile only", "union-pass work cap", ""},
    {"resources", "repair_max_physical_origins", "", "repair_max_physical_origins", "N >= 0", "origins", "completion", "expert profile only", "cumulative refinement work cap (cumulative_cap)", ""},
    {"resources", "repair_max_state_classes", "", "repair_max_state_classes", "N >= 0", "state classes", "completion", "expert profile only", "cumulative refinement work cap", ""},
    {"resources", "repair_max_piece_hits", "", "repair_max_piece_hits", "N >= 0", "piece hits", "completion", "expert profile only", "cumulative refinement work cap", ""},
    {"resources", "repair_max_latent_pairs", "", "repair_max_latent_pairs", "N >= 0", "pairs", "completion", "expert profile only", "cumulative refinement work cap", ""},
    {"resources", "max_union_cached_placements_consumed", "", "max_union_cached_placements_consumed", "N >= 0", "placements", "completion", "expert profile only", "served-work capacity of the placement reuse route", ""},
    {"resources", "max_union_cached_origins_consumed", "", "max_union_cached_origins_consumed", "N >= 0", "origins", "completion", "expert profile only", "served-work capacity of the origin reuse route", ""},
    {"resources", "max_union_cached_states_consumed", "", "max_union_cached_states_consumed", "N >= 0", "states", "completion", "expert profile only", "served-work capacity of the state reuse route", ""},
    {"resources", "max_union_cached_prefilters_consumed", "", "max_union_cached_prefilters_consumed", "N >= 0", "prefilters", "completion", "expert profile only", "served-work capacity of the prefilter reuse route", ""},
    {"resources", "adaptive_fine_max_served_placements", "", "adaptive_fine_max_served_placements", "N >= 0", "placements", "completion", "expert profile only", "served-work capacity of the fine adaptive discovery store", ""},
    {"resources", "emission_cache_bytes", "", "emission_cache_bytes", "N >= 0", "bytes", "computation", "expert profile only", "retained emission cache (reuse only; not a work cap)", ""},
    {"resources", "repair_origin_cache_bytes", "", "repair_origin_cache_bytes", "N >= 0", "bytes", "computation", "expert profile only", "retained refinement origin cache", ""},
    {"resources", "repair_tie_seed_cache_bytes", "", "repair_tie_seed_cache_bytes", "N >= 0", "bytes", "computation", "expert profile only", "retained tie-branch seed cache", ""},
    {"resources", "union_placement_cache_bytes", "", "union_placement_cache_bytes", "N >= 0", "bytes", "computation", "expert profile only", "retained union placement cache", ""},
    {"resources", "union_origin_cache_bytes", "", "union_origin_cache_bytes", "N >= 0", "bytes", "computation", "expert profile only", "retained union origin cache", ""},
    {"resources", "union_prefilter_cache_bytes", "", "union_prefilter_cache_bytes", "N >= 0", "bytes", "computation", "expert profile only", "retained union prefilter cache", ""},
    {"resources", "start_sort_scratch_bytes", "", "start_sort_scratch_bytes", "N >= 0", "bytes", "computation", "expert profile only", "start-sort scratch budget", ""},
    {"resources", "indexed_union_store_bytes", "", "indexed_union_store_bytes", "N >= 0", "bytes", "computation+completion", "expert profile only", "indexed union store; a store too small for the required data ends the run with a recorded resource refusal", ""},
    {"resources", "adaptive_fine_store_bytes", "", "adaptive_fine_store_bytes", "N >= 0", "bytes", "computation+completion", "expert profile only", "fine adaptive discovery store (refusal when too small)", ""},
    {"resources", "marker_proposal_reuse_table_bytes", "", "marker_proposal_reuse_table_bytes", "N >= 0", "bytes", "computation", "expert profile only", "marker-proposal reuse table", ""},
    {"resources", "union_scoring_backend", "", "union_scoring_backend", "reference | indexed_fine_experimental", "-", "computation", "expert profile only", "exact union scoring backend (results gated identical on completed runs)", ""},
    {"resources", "adaptive_window_discovery", "", "adaptive_window_discovery", "reference | fine_experimental", "-", "computation", "expert profile only", "exact adaptive window discovery", ""},
    {"resources", "indexed_union_store_layout", "", "indexed_union_store_layout", "exact | legacy_i2a", "-", "computation", "expert profile only", "indexed store layout", ""},
    {"resources", "indexed_union_evidence_layout", "", "indexed_union_evidence_layout", "reconstructed | fine_views_experimental", "-", "computation", "expert profile only", "indexed evidence layout", ""},
    {"resources", "adaptive_fused_seed", "", "adaptive_fused_seed", "0 | 1", "-", "computation", "expert profile only", "fused anchor seeding", ""},
    {"resources", "combination_prepared_library_terms", "", "combination_prepared_library_terms", "0 | 1", "-", "computation", "expert profile only", "library terms prepared once per combination", ""},
    {"resources", "gap_recursion_backend", "", "gap_recursion_backend", "scaled_reference | fixed_geometry_experimental", "-", "computation", "expert profile only", "gapped placement recursion", ""},
    {"resources", "dosage_census_backend", "", "dosage_census_backend", "reference | interned_experimental", "-", "computation", "expert profile only", "dosage census backend", ""},
    {"resources", "union_fragment_blocks", "", "union_fragment_blocks", "off | fine_core_4096_experimental", "-", "computation", "expert profile only", "blocked union route (needs the indexed backend, exact layout, fine views, origin cache and positive consumption caps)", ""},
    {"resources", "repair_round1_origin_admission", "", "repair_round1_origin_admission", "existing | reuse_only_experimental", "-", "computation", "expert profile only", "first-round origin admission", ""},
    {"resources", "union_progressive_schedule", "", "union_progressive_schedule", "off | experimental_v1", "-", "computation", "expert profile only", "progressive union schedule (not in the tested configuration)", ""},
    {"resources", "adaptive_emission_handoff", "", "adaptive_emission_handoff", "0 | 1", "-", "computation", "expert profile only", "union-to-adaptive emission handoff (not in the tested configuration)", ""},
    // 5. reporting and diagnostics
    {"reporting", "reporting contract (output_schema)", "", "output_schema", "sequence_first_v2", "-", "reporting", "fixed",
     "the reporting contract, recorded as output_schema in .settings.tsv; not a flag; a profile may only restate sequence_first_v2 (the evidence profiles do). The unrepaired development view (--adaptive-repair off) still reports legacy, recorded with source '(inactive: ...)'", ""},
    {"reporting", "--representation-guard on|off", "", "representation_guard", "on | off", "-", "reporting", "supported", "the representation-conflict check (feeds genotype_check CONFLICT)", ""},
    {"reporting", "compound_annotation", "", "compound_annotation", "off | contract_v2", "-", "reporting", "expert profile only", "the compound annotation (.sequence_compound.tsv); needs sequence_first_v2", ""},
    {"reporting", "compound_interval_declarations", "", "compound_interval_declarations", "a file", "-", "reporting", "expert profile only", "model-bound compound interval declarations (not in the tested configuration)", ""},
    {"reporting", "block_support_rule", "", "block_support_rule", "legacy | full_data_likelihood_v1_experimental | full_data_likelihood_v2_experimental", "-", "reporting", "expert profile only", "support semantics of the legacy block fields", ""},
    {"reporting", "library_sensitivity_audit", "", "library_sensitivity_audit", "0 | 1", "-", "reporting+computation", "expert profile only", "under --library auto, the 0.9x / 1.1x intensity audit behind LIBRARY_SENSITIVE", ""},
    {"reporting", "timing_detail", "", "timing_detail", "0 | 1", "-", "reporting", "expert profile only", "stage timers in the status", ""},
    {"reporting", "block_diagnostics", "", "block_diagnostics", "off | diagnostic_v1", "-", "reporting", "expert profile only", "per-block census route and JP tie diagnostics", ""},
    {"reporting", "emission_lookup_timing", "", "emission_lookup_timing", "0 | 1", "-", "reporting", "expert profile only", "measurement only (not for production calls)", ""},
    // 6. legacy and development compatibility
    {"legacy", "--max-candidates", "", "", "-", "-", "-", "legacy", "removed: refused with a pointer to --mosaic-candidates (every panel diplotype is always a candidate)", "refused"},
    {"legacy", "phase_audit_candidates", "", "phase_audit_candidates", "N >= 0", "candidates", "reporting", "expert profile only",
     "the phase audit set of the unrepaired development view; refused under adaptive repair", ""},
    {"legacy", "--output-schema", "", "", "-", "-", "-", "removed", "removed: genotyping reports the sequence-first v2 contract, recorded as output_schema in .settings.tsv (the unrepaired development view, --adaptive-repair off, still reports legacy, recorded the same way); historical contracts are not selectable", "refused"},
    {"legacy", "--dosage-mode, --dosage-audit", "", "", "-", "-", "-", "removed", "removed: dosage is an audit only (--dosage-check) and never replaces the prediction", "refused"},
    {"legacy", "--phase-audit-candidates", "", "", "-", "-", "-", "removed", "removed from the command line: the profile key phase_audit_candidates remains for the unrepaired development view", "refused"},
    {"legacy", "--repair-origin-cache-bytes, --repair-tie-seed-cache-bytes", "", "", "-", "-", "-", "removed", "removed from the command line: the profile keys repair_origin_cache_bytes and repair_tie_seed_cache_bytes remain", "refused"},
    {"legacy", "--dosage-allocation-test, --library-diagnostic", "", "", "-", "-", "-", "removed", "removed: development diagnostics, not available from the command line", "refused"},
    {"legacy", "diagnostic_* and test_* profile keys", "", "", "-", "-", "-", "development",
     "test hooks (diagnostic_union_evidence_dump, diagnostic_union_kernel_dump, diagnostic_refinement_memory, diagnostic_refinement_replay, test_interrupt_point, test_library_sensitivity_fail_view); not user settings", "off"},
};

std::string human_bytes(const std::string& raw) {
    if (raw.empty() || raw.find_first_not_of("0123456789") != std::string::npos) return raw;
    const unsigned long long n = std::stoull(raw);
    const unsigned long long gib = 1ull << 30, mib = 1ull << 20;
    if (n == 0) return "0";
    if (n % gib == 0) return std::to_string(n / gib) + " GiB";
    if (n % mib == 0) return std::to_string(n / mib) + " MiB";
    return raw;
}
std::string setting_display(const OptionEntry& e, const std::string& raw) {
    const std::string values = e.values;
    std::string v = raw;
    if (values == "on | off") v = raw == "1" ? "on" : raw == "0" ? "off" : raw;
    if (std::string(e.units) == "bytes") v = human_bytes(raw);
    return v;
}
std::string effective_default(const OptionEntry& e) {
    if (*e.fixed_default) return e.fixed_default;
    static const GenotypeCommandOptions built_in;
    const std::string key = e.key;
    const auto& pub = genotype_public_settings();
    const bool is_public = std::find(pub.begin(), pub.end(), key) != pub.end();
    const std::string b = setting_display(e, is_public ? genotype_setting_value(built_in, key) : experimental_setting_value(built_in, key));
    static const auto preset = genotype_preset_settings("short-read");
    const auto p = std::find_if(preset.begin(), preset.end(), [&](const auto& kv) { return kv.first == key; });
    if (p == preset.end()) return b;
    std::string shown = setting_display(e, p->second);
    if (key == "repair_max_neighbours" && p->second == "model") shown = "model (the model-derived bound)";
    return shown == b ? shown : shown + " (built-in " + b + ")";
}
std::string wrap(const std::string& text, std::size_t indent, std::size_t width) {
    std::string out, line;
    std::istringstream in(text);
    for (std::string word; in >> word;) {
        if (!line.empty() && indent + line.size() + 1 + word.size() > width) {
            out += std::string(indent, ' ') + line + "\n";
            line.clear();
        }
        line += (line.empty() ? "" : " ") + word;
    }
    if (!line.empty()) out += std::string(indent, ' ') + line + "\n";
    return out;
}
}  // namespace

std::vector<GenotypeOptionRow> genotype_option_inventory() {
    std::vector<GenotypeOptionRow> rows;
    for (const OptionEntry& e : kOptionEntries) {
        rows.push_back({e.group, e.option, e.aliases, e.key, e.values, effective_default(e), e.units, e.effect, e.status, e.text});
    }
    return rows;
}

std::string genotype_option_inventory_markdown() {
    std::ostringstream out;
    out << "| group | option | aliases | profile key | values | effective default (short-read) | units | changes | status |\n"
        << "|---|---|---|---|---|---|---|---|---|\n";
    const auto cell = [](std::string s) {
        std::string o;
        for (char c : s) o += c == '|' ? std::string("\\|") : std::string(1, c);
        return o.empty() ? std::string("—") : o;
    };
    for (const GenotypeOptionRow& r : genotype_option_inventory()) {
        out << "| " << r.group << " | `" << cell(r.option) << "` | " << cell(r.aliases) << " | " << cell(r.key) << " | " << cell(r.values)
            << " | " << cell(r.effective_default) << " | " << cell(r.units) << " | " << cell(r.effect) << " | " << r.status << " |\n";
    }
    return out.str();
}

std::string genotype_advanced_help() {
    std::ostringstream out;
    out << "panvar genotype: every option, by purpose\n\n"
           "Precedence: explicit flag > profile (--profile) > preset > built-in default, whatever the argument order;\n"
           "every effective setting and its source is written to <prefix>.settings.tsv. 'changes' says what an option\n"
           "can alter: inference (the prediction may change), computation (time or memory only; results gated\n"
           "identical), reporting (outputs or annotations only), completion (whether or how a run completes: caps,\n"
           "guards, budget stops). 'expert profile only' controls have no flag: set them in a --profile file.\n"
           "Defaults are the effective short-read values (the built-in value is shown where it differs).\n\n"
           "Profiles (--profile FILE): optional overrides, never required; an ordinary run needs none. One file per\n"
           "run (a repeated --profile is refused). Each line is a literal key<TAB>value; blank lines and lines\n"
           "starting with # are ignored; one row 'profile<TAB>NAME' is required and names the file in the records;\n"
           "a key set twice is refused; an unknown key is refused. A profile cannot select a reporting contract\n"
           "(output_schema may only restate sequence_first_v2) and cannot set the dosage stage (dosage_mode). A flag\n"
           "given beside the profile overrides it. Minimal working profile (TAB-separated):\n"
           "    profile<TAB>my_overrides\n"
           "    footprint_limit_bytes<TAB>4294967296\n";
    for (const auto& [group, title] : kOptionGroups) {
        if (group == "everyday") continue;
        out << "\n" << title << "\n";
        for (const OptionEntry& e : kOptionEntries) {
            if (group != e.group) continue;
            out << "  " << e.option << (*e.aliases ? std::string("   (alias ") + e.aliases + ")" : std::string()) << "\n";
            out << wrap(e.text, 6, 112);
            std::string facts = std::string("values: ") + e.values + "; default: " + effective_default(e);
            if (std::string(e.units) != "-") facts += std::string("; units: ") + e.units;
            facts += std::string("; changes: ") + e.effect + "; " + e.status;
            if (*e.key && std::string(e.option).rfind("--", 0) == 0) facts += std::string("; profile key ") + e.key;
            out << wrap(facts, 6, 112);
        }
    }
    out << "\nTerms (read before interpreting results)\n"
        << wrap("Donor redraw: --switch-probability is the copying prior's per-boundary donor-redraw probability, a "
                "regularisation of the panel model, not a biological recombination probability.", 2, 112)
        << wrap("Mosaic candidates: --mosaic-candidates sets the INITIAL off-panel mosaics; the bounded refinement may "
                "generate and score further mosaics beyond that number.", 2, 112)
        << wrap("Refinement is bounded: local moves are scored exactly within declared caps; it is not a globally "
                "exhaustive search, and a completed ranking is conditional on what was explored.", 2, 112)
        << wrap("min_variant_bp constructs the variant sites (blocks); it is not a filter on final calls.", 2, 112)
        << wrap("Library estimation learns the insert distribution and the fragment-start intensity only; the "
                "substitution, outlier and gap parameters stay fixed.", 2, 112)
        << wrap("Support is conditional on the explored candidates: READS means no explored rival projecting to a "
                "different genotype was within bounds of the winner's fragment term.", 2, 112)
        << wrap("Prior-dependent (read_support PRIOR_DEPENDENT): the reads alone did not separate the call from an "
                "explored rival. It does not establish that the prior alone caused the selection.", 2, 112)
        << wrap("Phase is assessed once per locus (the decision's phase support), not independently at each boundary.", 2, 112)
        << wrap("Locus summaries must be read with their coverage counts and not_assessed list: a positive value covers "
                "only the assessed blocks ([assessed/applicable] marks partial coverage).", 2, 112)
        << wrap("Nothing here is calibrated reliability; there is no quality threshold.", 2, 112);
    out << "\nFull inventory with defaults: docs/guides/genotype_short_read.md (option inventory).\n";
    return out.str();
}

namespace {
// M15-J: options removed from the public interface. Each is refused by name, with what replaces it, before anything is read
// or written; none is silently ignored. The behaviour behind them stays reachable through the internal API
// (GenotypeCommandOptions) where regression tests and diagnostics need it.
const std::map<std::string, std::string>& removed_genotype_options() {
    static const std::map<std::string, std::string> removed = {
        {"--output-schema", "--output-schema was removed: genotyping reports the sequence-first v2 contract, recorded as "
                            "output_schema in .settings.tsv; historical reporting contracts are not selectable"},
        {"--dosage-mode", "--dosage-mode was removed: dosage is an audit only and never replaces the prediction; use "
                          "--dosage-check to audit the reported alleles"},
        {"--dosage-audit", "--dosage-audit was removed: dosage is an audit only and never replaces the prediction; use "
                           "--dosage-check to audit the reported alleles"},
        {"--phase-audit-candidates", "--phase-audit-candidates was removed from the command line: set the profile key "
                                     "phase_audit_candidates (unrepaired development view only)"},
        {"--repair-origin-cache-bytes", "--repair-origin-cache-bytes was removed from the command line: set the profile key "
                                        "repair_origin_cache_bytes"},
        {"--repair-tie-seed-cache-bytes", "--repair-tie-seed-cache-bytes was removed from the command line: set the profile "
                                          "key repair_tie_seed_cache_bytes"},
        {"--dosage-allocation-test", "--dosage-allocation-test was removed: a development diagnostic of the dosage stage, "
                                     "not available from the command line"},
        {"--library-diagnostic", "--library-diagnostic was removed: a development diagnostic, not available from the "
                                 "command line"},
        {"--max-candidates", "--max-candidates is now --max-mosaic-candidates: every panel diplotype is always a candidate, "
                             "and K counts only the off-panel mosaics added to it"},
    };
    return removed;
}
}  // namespace

std::optional<int> parse_genotype_arguments(const std::vector<std::string>& args, GenotypeCommandOptions& options,
                                            std::ostream& out) {
    const auto need = [&](std::size_t& i) {
        if (i + 1 >= args.size()) throw std::runtime_error("missing value for " + args[i]);
        return args[++i];
    };
    // A technology preset is applied first, then the profile, wherever they appear; a flag given beside them overrides
    // both and is recorded as an override, so the status always says what the run actually ran under. An omitted
    // --preset means --preset short-read: the production algorithm is the default, and the unrepaired development view
    // needs an explicit --adaptive-repair off. The preset is applied once, the last --preset given or the default.
    options.preset_name = "short-read";
    options.preset_explicit = false;
    // Removed options are refused first, whatever else the command says (M15-J). Every argument is checked, so a value
    // spelled exactly like a removed option is refused too; no path or name needs that spelling.
    for (const std::string& a : args) {
        const auto it = removed_genotype_options().find(a);
        if (it != removed_genotype_options().end()) throw std::runtime_error(it->second);
    }
    // One profile per run: its file is what the status names (profile, profile_name, profile_sha256), so a second file
    // would apply settings the record does not identify.
    {
        std::size_t profiles = 0;
        for (const std::string& a : args) profiles += a == "--profile" ? 1 : 0;
        if (profiles > 1) {
            throw std::runtime_error("--profile was given " + std::to_string(profiles) + " times: give one profile file "
                                     "(merge the settings into one file; a key set twice is refused)");
        }
    }
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (args[i] != "--preset") continue;
        options.preset_name = need(i);
        options.preset_explicit = true;
    }
    {
        const std::string context = "preset " + options.preset_name + (options.preset_explicit ? "" : " (default)");
        for (const auto& [key, value] : genotype_preset_settings(options.preset_name)) {
            apply_genotype_setting(context, key, value, options);
            options.preset_keys.push_back(key);
        }
    }
    for (std::size_t i = 0; i < args.size(); ++i) {
        if (args[i] == "--profile") apply_genotype_profile(need(i), options);
    }
    const bool profiled = !options.profile_path.empty();
    std::optional<std::string> adaptive_repair_flag, repair_rounds_flag, tie_branching_flag;
    for (std::size_t i = 0; i < args.size(); ++i) {
        const std::string& flag = args[i];
        if (flag == "--profile" || flag == "--preset") { ++i; continue; }
        if (flag == "--gfa" || flag == "--graph" || flag == "-i") options.gfa_path = need(i);
        else if (flag == "--bubbles") options.bubbles_path = need(i);
        else if (flag == "--reads1") options.reads1_path = need(i);
        else if (flag == "--reads2") options.reads2_path = need(i);
        else if (flag == "--reference-path" || flag == "-r") options.reference_path = need(i);
        else if (flag == "--out-prefix" || flag == "--out" || flag == "-o") options.out_prefix = need(i);
        else if (flag == "--library") options.library_path = need(i);
        else if (flag == "--threads") apply_genotype_setting("flag --threads", "threads", need(i), options);   // workers
        else if (flag == "--panel-curation") options.curation_manifest_path = need(i);
        else if (flag == "--max-mosaic-candidates" || flag == "--mosaic-candidates") {
            apply_genotype_setting("flag " + flag, "max_mosaic_candidates", need(i), options);
        }
        else if (flag == "--repair-rounds") repair_rounds_flag = need(i);   // resolved with --adaptive-repair below
        else if (flag == "--adaptive-repair") {
            adaptive_repair_flag = need(i);
            if (*adaptive_repair_flag != "on" && *adaptive_repair_flag != "off") {
                throw std::runtime_error("--adaptive-repair takes on or off");
            }
        }
        else if (advanced_setting(flag) != nullptr) {
            // An advanced scientific setting: the profile key's own parser and validation, recorded as a flag.
            const AdvancedFlag& a = *advanced_setting(flag);
            std::string v = need(i);
            if (a.on_off) {
                if (v != "on" && v != "off") throw std::runtime_error(flag + " takes on or off: " + v);
                v = v == "on" ? "1" : "0";
            }
            apply_genotype_setting("flag " + flag, a.key, v, options);
            if (std::string(a.key) == "repair_tie_branching") tie_branching_flag = v;
        }
        else if (flag == "--dosage-check") options.dosage_check_flag = true;   // M15-D: resolved after every flag
        else if (flag == "--build-panel-index") options.build_panel_index_path = need(i);   // M15-E
        else if (flag == "--panel-index") options.panel_index_path = need(i);               // M15-E
        else if (flag == "-q" || flag == "--quiet") options.quiet = true;
        else if (flag == "-h" || flag == "--help") {
            out <<
                "Usage:\n"
                "  panvar genotype --graph <locus.gfa> --reads1 <R1.fq[.gz]> --reads2 <R2.fq[.gz]> --threads 4 --out <prefix>\n"
                "  panvar genotype --graph <locus.gfa> --build-panel-index <locus.pvpl>        (once per locus)\n"
                "  panvar genotype --graph <locus.gfa> --panel-index <locus.pvpl> --reads1 ... --reads2 ... --out <prefix>\n\n"
                "Genotypes one structural locus from paired short reads: a diplotype (two whole-locus haplotypes) over the\n"
                "graph's panel, projected to one allele pair per block. EXPERIMENTAL. No profile is needed: the short-read\n"
                "preset is applied, and the library is estimated from the reads (insert size and fragment-start intensity\n"
                "only; the error parameters stay fixed; an invalid estimate refuses the run, never a fallback).\n\n"
                "Everyday options:\n"
                "  --graph FILE               the locus graph (alias -i)\n"
                "  --reads1 FILE --reads2 FILE  paired FASTQ (.gz accepted), mates in the same order\n"
                "  --out PREFIX               output prefix (default genotype; alias -o)\n"
                "  --reference-path NAME      coordinate reference (default: the single path of the graph's RS:Z sample;\n"
                "                             refused with the candidates listed when absent or ambiguous)\n"
                "  --threads N                worker threads (default 1); results do not depend on them\n"
                "  --dosage-check             optional copy-number audit of the reported alleles. It never changes the\n"
                "                             prediction and adds run time (not necessarily cheaper than the call); without\n"
                "                             it dosage_check is NOT_REQUESTED\n"
                "  --build-panel-index FILE   build the sample-independent panel index for this graph, reference,\n"
                "                             curation, bubbles and marker settings, then stop; never overwritten\n"
                "  --panel-index FILE         genotype against a built index: it saves locus preparation, not the rest of\n"
                "                             the run. --graph is still required and still read and hashed; the index's\n"
                "                             integrity, versions and identity are verified before the output prefix is\n"
                "                             touched, and a missing, corrupt or incompatible index is refused (exit 2),\n"
                "                             never rebuilt. Rebuild it when the graph, reference, curation, bubbles,\n"
                "                             --min-variant-bp, --kmer-size or --syncmer-s change\n"
                "  -q, --quiet                no progress messages (warnings are still printed)\n"
                "  -h, --help                 this help\n"
                "  --help-advanced            every option by purpose: library file, curation and bubbles; the copying\n"
                "                             prior and marker model; candidates and bounded refinement; resource limits,\n"
                "                             caches and backends; reporting; legacy and development controls\n\n"
                "Outcomes: PREDICTION (a completed ranking; exit 0), PROVISIONAL (the best fully evaluated candidate after\n"
                "a budget stop or declared cap: the search is incomplete; exit 3), NO_PREDICTION (exit 1).\n"
                "Sequences: .sequence_haplotypes.tsv always gives each haplotype's SHA-256 and origin; when both haplotypes\n"
                "are on the panel, no FASTA is needed (extraction references point into the graph); when either is\n"
                "off-panel, .sequence_haplotypes.fa holds both.\n\n"
                "Primary assessment (assessment_v1), in .sequence_locus.tsv and .sequence_blocks.tsv (mirrored in the locus\n"
                "and block tables). Every check concerns the published prediction:\n"
                "  sequence_resolution  RESOLVED, GENOTYPE_RESOLVED_SEQUENCE_SET, SEQUENCE_TIED (sequence ambiguity)\n"
                "  genotype_check       NONE_DETECTED, CONFLICT or NOT_ASSESSED (+ coverage)\n"
                "  dosage_check         AGREES, SUSPECT, DISAGREES, UNAVAILABLE, NOT_REQUESTED or NOT_APPLICABLE; corrections\n"
                "                       are listed apart, as diagnostics (dosage_alternatives)\n"
                "  read_support         READS, PRIOR_DEPENDENT or NOT_ASSESSED: conditional on the explored candidates;\n"
                "                       prior-dependent does not mean the prior alone chose the call\n"
                "  phase                READ_SUPPORTED, PRIOR_ONLY, UNRESOLVED, NOT_ASSESSED or NOT_APPLICABLE, assessed once\n"
                "                       per locus (not independently at each boundary)\n"
                "  findings, not_assessed  every concrete finding and every check that could not run, with block ids\n"
                "Read a locus summary with its coverage and not_assessed list: a missing check is not a failure and never\n"
                "hides a finding. These assessments are not calibrated reliability, and there is no quality threshold. The\n"
                "older headline fields (genotype_filter, dosage_filter, phase_filter, call_summary) are kept as superseded\n"
                "compatibility diagnostics.\n\n"
                "Exit status: 0 PREDICTION, 3 PROVISIONAL, 1 NO_PREDICTION; 2 argument, input or output error; 124 time\n"
                "limit, 125 footprint limit (the 6 GiB footprint guard is sampled, not an allocation ceiling). The\n"
                "unrepaired development view (see --help-advanced) still reports the legacy contract: 0 a genotype was\n"
                "emitted, 1 no-call. Full guide: docs/guides/genotype_short_read.md.\n";
            return 0;
        }
        else if (flag == "--help-advanced") {
            out << genotype_advanced_help();
            return 0;
                } else {
            throw std::runtime_error("unknown genotype option: " + flag);
        }
    }
    // --adaptive-repair and --repair-rounds, resolved together so the result never depends on their order:
    //   neither: the preset's or profile's round cap; --adaptive-repair on: the adaptive repair, with the preset's 32
    //   rounds when nothing set a cap; --adaptive-repair off: no rounds and no tie branching (the unrepaired development
    //   view); --repair-rounds N alone: N > 0 the adaptive repair capped at N, N = 0 the unrepaired view.
    // An explicit contradiction between them is refused rather than resolved by whichever came last.
    if (repair_rounds_flag) apply_genotype_setting("flag --repair-rounds", "repair_rounds", *repair_rounds_flag, options);
    if (adaptive_repair_flag && *adaptive_repair_flag == "off") {
        if (repair_rounds_flag && options.repair_rounds > 0) {
            throw std::runtime_error("--adaptive-repair off contradicts --repair-rounds " + *repair_rounds_flag +
                                     ": the unrepaired view runs no repair rounds");
        }
        if (tie_branching_flag && *tie_branching_flag == "1") {
            throw std::runtime_error("--adaptive-repair off contradicts --repair-tie-branching on: tie branching is part "
                                     "of the adaptive repair");
        }
        apply_genotype_setting("flag --adaptive-repair", "repair_rounds", "0", options);
        apply_genotype_setting("flag --adaptive-repair", "repair_tie_branching", "0", options);
    } else if (adaptive_repair_flag) {
        if (repair_rounds_flag && options.repair_rounds == 0) {
            throw std::runtime_error("--adaptive-repair on contradicts --repair-rounds 0: the adaptive repair needs at "
                                     "least one round (use --adaptive-repair off for the unrepaired view)");
        }
        if (options.repair_rounds == 0) apply_genotype_setting("flag --adaptive-repair", "repair_rounds", "32", options);
    }
    // No repair rounds is the unrepaired view, whatever set it: tie branching belongs to the repair, so its effective value
    // is off too (in the settings record as in execution). An explicit --repair-tie-branching on beside it is refused.
    if (options.repair_rounds == 0 && options.repair_tie_branching) {
        if (tie_branching_flag && *tie_branching_flag == "1") {
            throw std::runtime_error("--repair-tie-branching on contradicts a repair round cap of 0 (" +
                                     options.setting_sources["repair_rounds"] +
                                     "): tie branching is part of the adaptive repair");
        }
        const auto source = options.setting_sources.find("repair_rounds");
        apply_genotype_setting(source == options.setting_sources.end() ? std::string("default") : source->second,
                               "repair_tie_branching", "0", options);
    }
    // M15-D: preset-sourced features whose requirements a flag or profile removed are switched off and recorded as
    // inactive (with the reason) instead of refusing an ordinary run; the same features set explicitly keep their refusals.
    {
        const auto from_preset = [&](const char* key) {
            const auto it = options.setting_sources.find(key);
            return it != options.setting_sources.end() && it->second.rfind("preset ", 0) == 0;
        };
        const auto deactivate = [&](const char* key, const char* value, const std::string& why) {
            const std::string source = options.setting_sources[key];
            apply_genotype_setting(source + " (inactive: " + why + ")", key, value, options);
        };
        if (options.repair_rounds == 0 && from_preset("output_schema") && options.output_schema != "legacy") {
            deactivate("output_schema", "legacy", "the unrepaired view has no pooled decision to report");
        }
        if (from_preset("early_incumbent") && options.early_incumbent != "off") {
            if (const auto why = early_incumbent_refusal(options)) deactivate("early_incumbent", "off", *why);
        }
        if (from_preset("compound_annotation") && options.compound_annotation != "off") {
            if (const auto why = compound_annotation_refusal(options)) deactivate("compound_annotation", "off", *why);
        }
        if (from_preset("repair_junction_proposals") && options.repair_junction_proposals != "off") {
            if (const auto why = junction_proposals_refusal(options)) deactivate("repair_junction_proposals", "off", *why);
        }
        if (from_preset("union_fragment_blocks") && options.union_fragment_blocks != "off") {
            if (const auto why = union_fragment_blocks_refusal(options)) deactivate("union_fragment_blocks", "off", *why);
        }
    }
    // M15-D: --dosage-check runs the existing dosage audit (never changing the prediction). It is consistent with an
    // explicit audit mode and contradicts any other explicit dosage setting; contradictions are refused before any output.
    if (options.dosage_check_flag) {
        if (options.dosage_mode_explicit && options.dosage_mode != "audit") {
            throw std::runtime_error("--dosage-check contradicts " + options.dosage_mode_source + " (dosage mode " + options.dosage_mode +
                                     "): --dosage-check runs the dosage audit, which never changes the prediction");
        }
        options.dosage_mode_source = options.dosage_mode_explicit ? options.dosage_mode_source + " and flag --dosage-check" : "flag --dosage-check";
        options.dosage_mode = "audit";
        options.dosage_mode_explicit = true;
    }
    // S1-B: sequence_first_v1 reports the adaptive repair's pooled decision; every other route is refused explicitly.
    if (sequence_first(options)) {
        if (options.repair_rounds == 0) {
            throw std::runtime_error("output_schema " + options.output_schema + " needs the adaptive repair (repair rounds > 0); the "
                                     "unrepaired view has no pooled decision to report");
        }
        if (options.diagnostic_union_only) {
            throw std::runtime_error("output_schema " + options.output_schema + " is not supported with diagnostic_union_only");
        }
    }
    // S1-D: under sequence_first_v2 dosage is diagnostic: audit unless a mode was declared; off is allowed (not
    // evaluated); repair is refused -- historical repair reproduction belongs to the legacy contract. S1-D.1: the one
    // resolver execute_genotype_command applies too.
    if (const std::optional<std::string> refused = resolve_dosage_contract(options)) throw std::runtime_error(*refused);
    if (const std::optional<std::string> refused = early_incumbent_refusal(options)) throw std::runtime_error(*refused);
    if (const std::optional<std::string> refused = union_fragment_blocks_refusal(options)) throw std::runtime_error(*refused);
    if (const std::optional<std::string> refused = refinement_diagnostic_refusal(options)) throw std::runtime_error(*refused);
    if (const std::optional<std::string> refused = joint_search_refusal(options)) throw std::runtime_error(*refused);
    if (const std::optional<std::string> refused = junction_proposals_refusal(options)) throw std::runtime_error(*refused);   // M9-G
    if (const std::optional<std::string> refused = compound_annotation_refusal(options)) throw std::runtime_error(*refused);   // M9-J
    if (profiled) {
        // Every setting a flag gave beside the profile, with the value it took.
        for (const std::string& flag : args) {
            if (flag == "--threads") options.overridden_settings.push_back("threads=" + std::to_string(options.threads));
            if (flag == "--repair-rounds" || flag == "--adaptive-repair") {
                options.overridden_settings.push_back("repair_rounds=" + std::to_string(options.repair_rounds));
            }
            if (flag == "--max-mosaic-candidates" || flag == "--mosaic-candidates") {
                options.overridden_settings.push_back("max_mosaic_candidates=" +
                                                      std::to_string(options.max_mosaic_candidates));
            }
            if (const AdvancedFlag* a = advanced_setting(flag)) {
                options.overridden_settings.push_back(std::string(a->key) + "=" + genotype_setting_value(options, a->key));
            }
        }
    }
    // M15-E: the panel-index modes. Building runs no sample, so a sample-only argument is refused rather than ignored;
    // building and loading in one invocation is refused (build once, then genotype each sample with --panel-index).
    if (!options.build_panel_index_path.empty()) {
        if (!options.panel_index_path.empty()) {
            throw std::runtime_error("--build-panel-index and --panel-index cannot be combined: build the index once, then "
                                     "genotype each sample with --panel-index");
        }
        static const std::set<std::string> sample_only = {"--reads1", "--reads2", "--library", "--out", "--out-prefix", "-o",
                                                          "--dosage-check"};
        for (const std::string& a : args) {
            if (sample_only.count(a) != 0) {
                throw std::runtime_error("--build-panel-index builds the sample-independent index only and runs no sample: " + a +
                                         " is a sample argument (give it when genotyping with --panel-index)");
            }
        }
    }
    validate_genotype_settings(options);   // every source applied: cross-field rules, before any work
    return std::nullopt;
}

namespace {
// The command's argument resolution: parsing, the graph's declared reference sample, the library default. Returns an
// exit code when parsing ends the command (help), otherwise fills `options`.
std::optional<int> resolve_genotype_command(const std::vector<std::string>& args, GenotypeCommandOptions& options) {
    if (const std::optional<int> help = parse_genotype_arguments(args, options, std::cout)) return *help;
    if (options.reference_path.empty() && !options.gfa_path.empty()) {
        // The graph's declared reference sample (a header RS:Z tag): its single path is the coordinate reference. No
        // declaration, or several paths for it, is refused with the candidates listed; nothing is picked silently.
        GzLineReader reader(options.gfa_path);
        if (!reader.ok()) throw std::runtime_error("cannot read " + options.gfa_path);
        std::string line, sample;
        std::vector<std::string> paths;
        while (reader.getline(line)) {
            if (line.rfind("H\t", 0) == 0) {
                const std::size_t at = line.find("RS:Z:");
                if (at != std::string::npos) sample = line.substr(at + 5, line.find('\t', at) == std::string::npos ? std::string::npos : line.find('\t', at) - at - 5);
            } else if (line.rfind("P\t", 0) == 0) {
                paths.push_back(line.substr(2, line.find('\t', 2) - 2));
            } else if (line.rfind("W\t", 0) == 0) {
                std::vector<std::string> f;
                std::stringstream row(line);
                for (std::string c; std::getline(row, c, '\t');) f.push_back(c);
                if (f.size() > 5) paths.push_back(f[1] + "#" + f[2] + "#" + f[3] + ":" + f[4] + "-" + f[5]);
            }
        }
        std::vector<std::string> of_sample;
        for (const std::string& p : paths) if (!sample.empty() && (p == sample || p.rfind(sample + "#", 0) == 0)) of_sample.push_back(p);
        if (of_sample.size() == 1) {
            options.reference_path = of_sample.front();
        } else {
            std::string listed;
            const std::vector<std::string>& show = of_sample.empty() ? paths : of_sample;
            for (std::size_t k = 0; k < show.size() && k < 8; ++k) listed += "\n  " + show[k];
            if (show.size() > 8) listed += "\n  ... (" + std::to_string(show.size()) + " paths)";
            throw std::runtime_error(std::string(sample.empty() ? "the graph declares no reference sample (no RS:Z header tag)"
                                                                : "the graph's reference sample '" + sample + "' has " + std::to_string(of_sample.size()) + " paths") +
                                     "; name the coordinate reference with --reference-path. Candidates:" + listed);
        }
    }
    // An omitted --library is --library auto: the estimate, its validation (LIBRARY_REFUSED, never a fallback) and the
    // sensitivity audit. --library FILE is the exact reproducibility route.
    if (options.library_path.empty()) {
        options.library_path = "auto";
        options.library_omitted = true;
    }
    return std::nullopt;
}
}  // namespace

std::optional<std::string> resolve_dosage_contract(GenotypeCommandOptions& options) {
    if (options.dosage_mode != "off" && options.dosage_mode != "audit" && options.dosage_mode != "repair") {
        return "dosage_mode takes off, audit or repair (got '" + options.dosage_mode + "')";
    }
    if (!options.dosage_mode_explicit && options.dosage_mode != "off") {
        options.dosage_mode_explicit = true;   // an API caller set a non-default mode: declared
        if (options.dosage_mode_source == "default") options.dosage_mode_source = "caller_options";
    }
    // M15-D: no schema implies a dosage stage. Without --dosage-check (or an explicit mode) dosage is off and reported as
    // not evaluated; the audit runs only when requested and never changes the prediction.
    if (options.output_schema != "sequence_first_v2") return std::nullopt;
    if (options.dosage_mode == "repair") {
        return "dosage_mode repair is not available under output_schema sequence_first_v2: dosage is diagnostic there "
               "(--dosage-check, or off) and never changes the prediction; historical repair results are reproduced only "
               "through the internal API with the legacy contract";
    }
    if (options.dosage_mode_explicit && options.dosage_mode_source == "default") {
        options.dosage_mode_source = "caller_options";   // S1-D.1: declared off through dosage_mode_explicit
    }
    return std::nullopt;
}

namespace {
// S1-D.1: the sequence_first_v2 FASTA policy over the whole prefix. The sequence report's FASTA is the only
// prediction FASTA: the legacy haplotype FASTA is suppressed (after the dosage stage and the locus table, its internal
// consumers, have read it), and so is every empty dosage FASTA. A non-empty dosage FASTA is an audit diagnostic
// (repair-candidate representatives) and stays, linked from the locus table.
void suppress_empty_dosage_fasta(const std::string& prefix) {
    for (const char* suffix : {".dosage_repaired.fa", ".dosage_original.fa", ".dosage_repaired_block_alleles.fa"}) {
        std::error_code ignored;
        const std::string file = prefix + suffix;
        if (std::filesystem::exists(file, ignored) && std::filesystem::file_size(file, ignored) == 0) {
            std::filesystem::remove(file, ignored);
        }
    }
}
void suppress_legacy_fasta(const std::string& prefix) {
    std::error_code ignored;
    std::filesystem::remove(prefix + ".mosaic_haplotypes.fa", ignored);
    std::filesystem::remove(prefix + ".original.mosaic_haplotypes.fa", ignored);
}
struct DosageStageScope {   // S1-D.1: marks the body run as followed by the dosage stage, and always resets it
    DosageStageScope() { g_dosage_stage_follows = true; }
    ~DosageStageScope() { g_dosage_stage_follows = false; }
};
}  // namespace

// The per-sample command after resolution -- genotype, dosage, locus table -- shared by the ordinary command and the
// prepared-locus route (I3-A).
int execute_genotype_command(const GenotypeCommandOptions& requested, const PreparedLocus* locus) {
    // The dosage stage and the locus table run after every existing output is committed (sections 20 and 22).
    // S1-D.1: the dosage contract is resolved here as at argument resolution, so a direct API caller and the
    // prepared-locus route get the CLI's defaults; a refused repair returns before the run touches its prefix.
    GenotypeCommandOptions options = requested;
    if (const std::optional<std::string> refused = resolve_dosage_contract(options)) {
        std::cerr << "[genotype] " << *refused << '\n';
        return 2;
    }
    const bool v2_fasta_policy = options.output_schema == "sequence_first_v2";
    const std::string& mode = options.dosage_mode;
    GenotypeDosageCapture capture;
    GenotypeRunStatus status;
    {
        const DosageStageScope scope;
        status = run_genotype_body(options, &capture, locus);
    }
    const bool called = status == GenotypeRunStatus::Complete || status == GenotypeRunStatus::ConditionalGenotypeEquivalent;
    const bool sequence_first = panvar::sequence_first(options);
    if (sequence_first && !g_sequence_published) {
        // S1-B/C: a run that ended without the adaptive completion's stop (no prediction, or an uncovered route)
        // publishes here; a watchdog stop never reaches here (no prediction, and its owned reports were removed).
        SequenceLegacy legacy = g_sequence_legacy;
        if (legacy.run_status == "NA") legacy.run_status = genotype_run_status_name(status);
        legacy.dosage_mode = options.dosage_mode;
        legacy.dosage_mode_source = dosage_mode_source(options);
        // A legacy call without a sequence report would be a route this schema does not cover: said, never hidden.
        const std::string reason =
            called && !g_sequence_report.built
                ? "UNSUPPORTED_ROUTE: the legacy caller emitted a genotype on a route " + options.output_schema + " does not report (" +
                      legacy.run_status + ")"
                : "no completed ranking and no safely available evaluated candidate (" + legacy.run_status +
                      (legacy.no_call_reason != "NA" ? ", " + legacy.no_call_reason : "") + ")";
        compound_prepare(g_sequence_report, options);   // M9-J: enabled; rows only if a pool or a provisional set them
        publish_sequence_report(g_sequence_report, legacy, options.out_prefix, reason);
        g_sequence_published = true;
    }
    if (sequence_first && g_sequence_report.publication_failed) {
        std::cerr << "[genotype] sequence report output failed: " << g_sequence_report.publication_error << '\n';
    }
    dosage::DosageShadowSummary summary;
    std::string audit_status = "NOT_RUN", audit_note = "none";
    std::vector<std::string> census_routes;   // M11-C: filled only by the census, only when block_diagnostics is on
    std::optional<double> declared_marker_background;
    std::optional<LibraryModel> library;
    std::vector<ReadPair> fragments;
    FragmentScoringOptions scoring;
    if (mode != "off") {
        // Audit writes only its own files and never changes the exit status, a call or an existing output. A failure
        // is recorded in its own status file.
        try {
            if (!called || !capture.model) {
                audit_note = std::string("the run emitted no genotype (") + genotype_run_status_name(status) + ")";
                atomic_write(options.out_prefix + ".dosage_status.tsv",
                             "key\tvalue\ndosage_audit\tshadow\ndosage_audit_status\tNOT_RUN\ndosage_audit_note\t" + audit_note + "\n");
            } else {
                // --library auto: the library the run estimated and scored with; --library FILE: re-read exactly as declared
                library = options.library_path == "auto" && capture.library ? *capture.library
                                                                           : read_library(options.library_path, declared_marker_background);
                fragments = read_fragments(options.reads1_path, options.reads2_path);
                scoring.mismatch_band = options.mismatch_band;
                scoring.minimum_piece_bp = options.minimum_piece_bp;
                scoring.max_piece_hits = options.max_piece_hits;
                scoring.max_proposed_starts_per_haplotype = options.max_proposed_starts_per_haplotype;
                scoring.max_physical_origins = options.max_physical_origins;
                scoring.max_state_classes = options.max_state_classes;
                scoring.max_mass_cells = options.max_mass_cells;
                scoring.emission_cache_bytes = options.emission_cache_bytes;
                scoring.gap_recursion_backend = gap_recursion_backend_of(options);   // I2-H: dosage too
                const std::vector<double> switches = run_switches(capture.model->blocks().size());   // M4-B
                dosage::DosageShadowOptions shadow_options;
                shadow_options.threads = options.threads;
                shadow_options.allocation_test = options.dosage_allocation_test;
                shadow_options.census_backend = options.dosage_census_interned ? dosage::CensusBackend::interned_experimental
                                                                               : dosage::CensusBackend::reference;   // I2-J
                shadow_options.census_backend_declared = options.setting_sources.count("dosage_census_backend") != 0;
                shadow_options.repair_max_physical_origins = options.repair_max_physical_origins;
                shadow_options.repair_max_state_classes = options.repair_max_state_classes;
                shadow_options.repair_max_piece_hits = options.repair_max_piece_hits;
                shadow_options.repair_max_latent_pairs = options.repair_max_latent_pairs;
                shadow_options.repair_origin_cache_bytes = 0;   // path on, nothing retained (measured: 2.77 vs 3.74 GiB)
                if (options.block_diagnostics == "diagnostic_v1") shadow_options.census_routes_out = &census_routes;   // M11-C
                summary = dosage::run_dosage_shadow(*capture.model, capture.node_sequences, fragments, *library, scoring, switches,
                                                    options.out_prefix, options.out_prefix, shadow_options);
                audit_status = summary.audit_status;
                audit_note = summary.note.empty() ? "none" : summary.note;
            }
        } catch (const std::exception& error) {
            audit_status = "FAILED";
            audit_note = error.what();
            atomic_write(options.out_prefix + ".dosage_status.tsv",
                         std::string("key\tvalue\ndosage_audit\tshadow\ndosage_audit_status\tFAILED\ndosage_audit_note\t") +
                             error.what() + "\n");
        }
        // M15-D.1: a requested audit that did not complete is said on stderr (also under --quiet: it is a warning, not
        // progress); the prediction is unchanged and DOSAGE_FILTER reports NOT_EVALUATED
        if (audit_status != "COMPLETE" && audit_status != "NOT_RUN") {
            std::cerr << "[genotype] WARNING: the requested dosage audit did not complete (" << audit_status << ": " << audit_note
                      << "); the prediction is unchanged and DOSAGE_FILTER is NOT_EVALUATED\n";
        }
    }
    if (options.block_diagnostics == "diagnostic_v1" && capture.model) {
        // M11-C (DIAGNOSTIC; nothing reads it): per-block census route and JP1.1 emission-tie condition. A value is written
        // only where it was computed; otherwise an explicit NOT_EVALUATED / UNAVAILABLE value with its reason.
        std::string census_absent;
        if (mode == "off") census_absent = "NOT_EVALUATED:dosage_audit_off";
        else if (!called) census_absent = "NOT_EVALUATED:no_genotype_emitted";
        else if (audit_status == "REFUSED_LIBRARY") census_absent = "NOT_EVALUATED:unsupported_library";
        else if (audit_status == "NOT_RUN") census_absent = "NOT_EVALUATED:no_rank1_call";
        else census_absent = "UNAVAILABLE:census_not_completed";
        const std::string jp_absent = options.repair_junction_proposals == "off" ? "NOT_EVALUATED:junction_proposals_off"
                                                                               : "NOT_EVALUATED:no_repair_search_report";
        std::ostringstream t;
        t << "# block diagnostics (M11-C, diagnostic_v1): census_route = the dosage census route of the block (model-derived; "
             "reported when the census completed, also if a later audit stage failed); jp_emission_tie = the JP1.1 eligibility "
             "condition (every proposal pair emission tied and more pairs than the local-pair limit) where a repair search "
             "evaluated it. Diagnostic only: not confidence, not correctness, not read by the caller.\n"
          << "# dosage_audit_status\t" << audit_status << "\n"
          << "block_index\tblock_id\tcensus_route\tjp_emission_tie\n";
        const auto& blocks = capture.model->blocks();
        for (std::size_t b = 0; b < blocks.size(); ++b) {
            const std::string route = census_routes.size() == blocks.size() ? census_routes[b] : census_absent;
            const std::string tie = g_junction_block_flags.size() == blocks.size() ? (g_junction_block_flags[b] ? "TIED_ELIGIBLE" : "NOT_ELIGIBLE")
                                                                                  : jp_absent;
            t << b << '\t' << blocks[b].block_id << '\t' << route << '\t' << tie << '\n';
        }
        atomic_write(options.out_prefix + ".block_diagnostics.tsv", t.str());
    }
    // Repair mode's rewrite builds every file before writing the first, so a refusal here leaves the caller's outputs
    // as they are; the locus table is the primary locus-level output in every mode.
    try {
        const PreparedGenotypeModel* model = capture.model ? &*capture.model : nullptr;
        dosage::ReportedFilterInputs filters;
        if (library) {
            filters.reads = &fragments;
            filters.library = &*library;
            filters.declared_marker_background = declared_marker_background;
            filters.kmer_size = options.kmer_size;
            filters.syncmer_s = options.syncmer_s;
            filters.scoring = scoring;
            filters.scoring.threads = options.threads;
            filters.repair_max_physical_origins = options.repair_max_physical_origins;
            filters.repair_max_state_classes = options.repair_max_state_classes;
            filters.repair_max_piece_hits = options.repair_max_piece_hits;
            filters.repair_max_latent_pairs = options.repair_max_latent_pairs;
            filters.switch_probability = options.switch_probability;
            if (g_transitions.distance_mode()) filters.switch_vector = g_transitions.r;   // M4-B
            filters.representation_guard = options.representation_guard;
        }
        dosage::DosageIntegration integration =
            dosage::apply_dosage_mode(model, mode, options.out_prefix, options.out_prefix,
                                      audit_status == "COMPLETE" ? &summary : nullptr, audit_status, audit_note,
                                      library ? &filters : nullptr);
        if (v2_fasta_policy) {   // S1-D.1: the audit has read the legacy FASTA; the locus table links what stays
            suppress_empty_dosage_fasta(options.out_prefix);
            integration.fasta_policy = "sequence_first_v2";
        }
        dosage::write_locus_genotype(model, options.out_prefix, integration);
        if (v2_fasta_policy) suppress_legacy_fasta(options.out_prefix);
        // M15-G: the primary assessment (assessment_v1) of the published v2 prediction, from this run's own outputs, once
        // every input it reads is final; a failure here is a reporting failure (exit 2, named below)
        if (v2_fasta_policy) assessment::write_assessment(options.out_prefix);
        if (options.library_diagnostic && model && options.library_path != "auto") {   // auto writes its own table
            std::optional<double> ignored;
            const LibraryModel declared = read_library(options.library_path, ignored);
            const std::vector<ReadPair> pairs = fragments.empty() ? read_fragments(options.reads1_path, options.reads2_path) : fragments;
            LibraryEstimateOptions estimate_options;
            estimate_options.threads = options.threads;
            estimate_options.marker_kmer = options.kmer_size;
            estimate_options.marker_syncmer = options.syncmer_s;
            atomic_write(options.out_prefix + ".library_estimate.tsv",
                         library_estimate_table(estimate_library(*model, pairs, estimate_options, &declared), &declared));
        }
    } catch (const std::exception& error) {
        std::cerr << "[genotype] locus output failed: " << error.what() << '\n';
        if (v2_fasta_policy) {   // S1-D.1: the policy holds on failure too
            suppress_empty_dosage_fasta(options.out_prefix);
            suppress_legacy_fasta(options.out_prefix);
        }
        // M15-D.1: the process exits 2 (an output error), so the records that state an exit say 2 as well. The sequence
        // report stays as published (its sequences are complete); the failure is named beside it.
        if (v2_fasta_policy) {   // M15-G: the primary view of the published prediction, best effort, from what was written
            try {
                assessment::write_assessment(options.out_prefix);
            } catch (...) {
            }
        }
        if (sequence_first) {
            const std::string why = std::string("FAILED: ") + error.what();
            for (const std::string file : {options.out_prefix + ".sequence_locus.tsv", options.out_prefix + ".mosaic_status.tsv"}) {
                try {
                    std::ifstream in(file);
                    if (!in) continue;
                    std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
                    const bool locus = file.size() > 19 && file.compare(file.size() - 19, 19, ".sequence_locus.tsv") == 0;
                    text = with_key_value(text, locus ? "exit_code" : "sequence_exit_code", "2");
                    if (!locus) text = with_key_value(text, "exit_code", "2");
                    text = with_key_value(text, "locus_output", why);
                    atomic_write(file, text);
                } catch (...) {
                }
            }
        }
        return 2;
    }
    if (sequence_first) return sequence_report_exit_code(g_sequence_report);   // S1-B: 0 / 3 / 1
    return called ? 0 : 1;
}

namespace {
// M15-E: --build-panel-index. The sample-independent preparation (structural model and marker panel) through the existing
// I3-C/I3-D writer: canonical bytes, atomic publication that never overwrites. No sample runs, no genotype prefix is
// cleared and no sample report is written; the only temporary (the staged sorted graph) sits next to the index and is
// removed. The build record goes to stdout.
int build_panel_index(const GenotypeCommandOptions& requested) {
    const std::string& path = requested.build_panel_index_path;
    try {
        std::error_code ec;
        if (std::filesystem::exists(path, ec)) {   // before any work; publish_prepared_locus refuses it too
            throw PreparedLocusFormatError(PreparedLocusFormatError::Kind::exists, path + " exists; an index is never overwritten");
        }
        GenotypeCommandOptions options = requested;
        options.out_prefix = path;   // the staged sorted graph: <index>.prepared_locus.sorted.gfa, removed after the build
        const Clock::time_point start = Clock::now();
        const PreparedLocusPreparation prepared = prepare_genotype_locus_record(options, /*persist_panel=*/true);
        const double prepare_seconds = seconds_since(start);
        const Clock::time_point write_start = Clock::now();
        const std::string bytes = encode_prepared_locus(prepared.record);
        publish_prepared_locus(bytes, path);
        const double write_seconds = seconds_since(write_start);
        const PreparedLocusKey& key = prepared.record.key;
        std::cout << std::setprecision(6) << "key\tvalue\npanel_index\t" << path << "\nschema_version\t" << kPreparedLocusSchemaVersion
                  << "\npreparation_algorithm_version\t" << kPreparationAlgorithmVersion << "\npanel_extension_version\t"
                  << kPanelExtensionVersion << "\nmarker_panel_algorithm_version\t" << kMarkerPanelAlgorithmVersion
                  << "\nidentity_gfa_sha256\t" << key.gfa_sha256 << "\nidentity_reference_path\t" << key.reference_path
                  << "\nresolved_reference\t" << prepared.record.resolved_reference << "\nidentity_bubbles\t" << key.bubbles
                  << "\nidentity_min_variant_bp\t" << key.min_variant_bp << "\nidentity_curation\t" << key.curation
                  << "\nidentity_kmer_size\t" << key.kmer_size << "\nidentity_syncmer_s\t" << key.syncmer_s
                  << "\nmodel_fingerprint\t" << prepared.record.model_fingerprint << "\nfile_bytes\t" << bytes.size()
                  << "\nfile_sha256\t" << sha256_bytes(bytes) << "\nprepare_seconds\t" << prepare_seconds << "\nwrite_seconds\t"
                  << write_seconds << "\ncontents\tsample-independent preparation only: structural model and marker panel\n";
        return 0;
    } catch (const PreparedLocusFormatError& error) {
        std::cerr << "[genotype] --build-panel-index " << path << " refused (" << prepared_locus_error_name(error.kind())
                  << "): " << error.what() << '\n';
        return 2;
    }
}

// M15-E: --panel-index. The index is loaded and verified in full (integrity, bounded decoding, schema and preparation
// versions, reconstruction digests) before anything else; the sample's own inputs are then checked against its identity
// (graph content, reference, bubbles, curation, min_variant_bp, k, s) by the prepared route, before the output prefix is
// touched. Any failure is a named refusal (exit 2); there is never a fallback to an in-memory preparation.
int run_with_panel_index(GenotypeCommandOptions& options) {
    const std::string& path = options.panel_index_path;
    std::shared_ptr<const PreparedLocus> locus;
    try {
        PreparedLocusLoadReport report;
        const Clock::time_point start = Clock::now();
        locus = load_prepared_locus(path, PreparedLocusLoadLimits{}, &report);
        options.panel_index_load_seconds = seconds_since(start);
        options.panel_index_file_bytes = report.file_bytes;
        options.panel_index_panel_source = report.panel_source;
        options.panel_index_panel_rebuild_reason = report.panel_rebuild_reason;
    } catch (const PreparedLocusFormatError& error) {
        std::cerr << "[genotype] --panel-index " << path << " refused (" << prepared_locus_error_name(error.kind())
                  << "): " << error.what() << "; nothing at the output prefix was touched\n";
        return 2;
    }
    if (locus->panel_source() != "loaded") {   // never described as a fully loaded panel
        std::cerr << "[genotype] --panel-index " << path << ": the marker panel was REBUILT from the stored model ("
                  << (options.panel_index_panel_rebuild_reason.empty() ? std::string("no panel section")
                                                                       : options.panel_index_panel_rebuild_reason)
                  << "); the structural model was loaded\n";
    }
    try {
        return execute_genotype_command(options, locus.get());
    } catch (const PreparedLocusMismatch& error) {
        std::cerr << "[genotype] --panel-index " << path << " is incompatible with this sample's inputs: " << error.what()
                  << "; nothing at the output prefix was touched (build an index for these inputs with --build-panel-index)\n";
        return 2;
    }
}
}  // namespace

int run_genotype_command(const std::vector<std::string>& args) {
    GenotypeCommandOptions options;
    if (const std::optional<int> done = resolve_genotype_command(args, options)) return *done;
    if (!options.build_panel_index_path.empty()) return build_panel_index(options);   // M15-E
    if (!options.panel_index_path.empty()) return run_with_panel_index(options);      // M15-E
    return execute_genotype_command(options, nullptr);
}

// I3-A: the same command against a prepared locus (no public CLI; development driver and tests).
int run_genotype_command_prepared(const std::vector<std::string>& args, const PreparedLocus& locus) {
    GenotypeCommandOptions options;
    if (const std::optional<int> done = resolve_genotype_command(args, options)) return *done;
    return execute_genotype_command(options, &locus);
}

std::shared_ptr<const PreparedLocus> prepare_genotype_locus(const std::vector<std::string>& args) {
    GenotypeCommandOptions options;
    if (resolve_genotype_command(args, options)) throw std::invalid_argument("prepare_genotype_locus: the command line ends in help");
    return prepare_genotype_locus(options);
}

PreparedLocusPreparation prepare_genotype_locus_record(const std::vector<std::string>& args, bool persist_panel) {
    GenotypeCommandOptions options;
    if (resolve_genotype_command(args, options)) throw std::invalid_argument("prepare_genotype_locus_record: the command line ends in help");
    return prepare_genotype_locus_record(options, persist_panel);
}

} // namespace panvar
