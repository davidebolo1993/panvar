#include "panvar/genotype_dosage_integrate.hpp"

#include "panvar/genotype_command.hpp"
#include "panvar/genotype_repair.hpp"
#include "panvar/ls_hmm.hpp"
#include "panvar/mosaic_spelling.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>

namespace panvar::dosage {

namespace {

bool exists(const std::string& path) { return static_cast<bool>(std::ifstream(path)); }

std::string slurp(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot read " + path);
    std::ostringstream s;
    s << in.rdbuf();
    return s.str();
}

std::vector<std::string> split(const std::string& line, char sep) {
    std::vector<std::string> f;
    std::string cell;
    std::stringstream row(line);
    while (std::getline(row, cell, sep)) f.push_back(cell);
    if (!line.empty() && line.back() == sep) f.push_back("");
    return f;
}

std::string join(const std::vector<std::string>& cells, const std::string& sep) {
    std::string out;
    for (std::size_t i = 0; i < cells.size(); ++i) out += (i ? sep : "") + cells[i];
    return out;
}

// A TSV kept as text: the header and every row's cells, so untouched cells are written back byte for byte.
struct Table {
    std::vector<std::string> header;
    std::vector<std::vector<std::string>> rows;
    std::size_t col(const std::string& name) const {
        const auto it = std::find(header.begin(), header.end(), name);
        if (it == header.end()) throw std::runtime_error("table has no column " + name);
        return static_cast<std::size_t>(it - header.begin());
    }
    bool has(const std::string& name) const { return std::find(header.begin(), header.end(), name) != header.end(); }
    std::string text() const {
        std::string out = join(header, "\t") + "\n";
        for (const auto& r : rows) out += join(r, "\t") + "\n";
        return out;
    }
};

Table read_table(const std::string& path) {
    Table t;
    std::istringstream in(slurp(path));
    std::string line;
    if (!std::getline(in, line)) return t;
    t.header = split(line, '\t');
    while (std::getline(in, line)) t.rows.push_back(split(line, '\t'));
    return t;
}

std::map<std::string, std::string> read_key_values(const std::string& path) {
    std::map<std::string, std::string> out;
    if (!exists(path)) return out;
    std::istringstream in(slurp(path));
    std::string line;
    while (std::getline(in, line)) {
        const std::size_t tab = line.find('\t');
        if (tab == std::string::npos) continue;
        out.emplace(line.substr(0, tab), line.substr(tab + 1));   // first occurrence wins, as the caller writes each once
    }
    return out;
}

std::string value(const std::map<std::string, std::string>& kv, const std::string& key) {
    const auto it = kv.find(key);
    return it == kv.end() ? "NA" : it->second;
}

// block_id -> named cells of <prefix>.dosage_blocks.tsv
std::map<std::string, std::map<std::string, std::string>> read_dosage_blocks(const std::string& path) {
    std::map<std::string, std::map<std::string, std::string>> out;
    if (!exists(path)) return out;
    const Table t = read_table(path);
    for (const auto& r : t.rows) {
        std::map<std::string, std::string> row;
        for (std::size_t i = 0; i < t.header.size() && i < r.size(); ++i) row[t.header[i]] = r[i];
        if (row.count("block_id")) out[row.at("block_id")] = row;
    }
    return out;
}

struct FastaRecord {
    std::string header;   // without '>'
    std::string sequence;
};

std::vector<FastaRecord> read_fasta(const std::string& path) {
    std::vector<FastaRecord> out;
    std::istringstream in(slurp(path));
    std::string line;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '>') out.push_back({line.substr(1), ""});
        else if (!out.empty()) out.back().sequence += line;
    }
    return out;
}

std::string format(double v) {
    std::ostringstream o;
    o << std::setprecision(17) << v;
    return o.str();
}

std::string basename_of(const std::string& path) {
    const std::size_t slash = path.find_last_of('/');
    return slash == std::string::npos ? path : path.substr(slash + 1);
}

std::string ids_or_none(const std::vector<std::string>& ids) { return ids.empty() ? "none" : join(ids, ","); }

std::string join_path(const MosaicPath& p) {
    std::string out;
    for (std::size_t i = 0; i < p.alleles.size(); ++i) out += (i ? "," : "") + std::to_string(p.alleles[i]);
    return out;
}

// the block-level dosage filter of section 20.4, from the stage's status and repair mode's application
std::pair<std::string, std::string> dosage_filter(const std::string& status, const std::string& reason, bool applied,
                                                  const std::string& not_applied_reason) {
    if (status == "CONSISTENT" || status == "NOT_SENSITIVE") return {"PASS", "none"};
    if (status == "REPAIRED") return applied ? std::make_pair(std::string("PASS"), std::string("none"))
                                             : std::make_pair(std::string("DOSAGE_REPAIR_NOT_APPLIED"), not_applied_reason);
    if (status == "CONTENT_REPAIRED") return {"DOSAGE_CONTENT_CORRECTED", reason};
    if (status == "TRIGGERED_UNREPAIRED") return {"DOSAGE_TRIGGERED_UNREPAIRED", reason};
    if (status == "AUTO_REPAIR") return {"DOSAGE_TRIGGERED_UNREPAIRED", "repair_not_run"};
    if (status == "SUSPECT") return {"DOSAGE_SUSPECT", reason};
    if (status == "UNDERPOWERED") return {"DOSAGE_UNDERPOWERED", reason};
    return {"DOSAGE_NOT_EVALUATED", reason};
}

} // namespace

const std::vector<std::string>& legacy_call_outputs() {
    static const std::vector<std::string> suffixes = {".genotypes.tsv",      ".block_genotypes.tsv", ".mosaic_haplotypes.fa",
                                                      ".mosaic_status.tsv",  ".mosaic_scores.tsv",   ".mosaic_candidates.tsv",
                                                      ".equivalence_set.tsv"};
    return suffixes;
}

const std::vector<std::string>& dosage_stage_outputs() {
    static const std::vector<std::string> suffixes = [] {
        std::vector<std::string> s = {".locus_genotype.tsv",
                                      ".dosage_census.tsv",
                                      ".dosage_features.tsv",
                                      ".dosage_feature_copies.tsv",
                                      ".dosage_blocks.tsv",
                                      ".dosage_repairs.tsv",
                                      ".dosage_repair_candidates.tsv",
                                      ".dosage_repair_scoring.tsv",
                                      ".dosage_allocation.tsv",
                                      ".dosage_allocation_equivalence.tsv",
                                      ".dosage_space.tsv",
                                      ".dosage_status.tsv",
                                      ".dosage_repaired.fa",
                                      ".dosage_original.fa",
                                      ".dosage_repaired_block_alleles.fa",
                                      ".dosage_provenance.tsv",
                                      ".dosage_provenance_panel.tsv"};
        for (const std::string& l : legacy_call_outputs()) s.push_back(".original" + l);
        return s;
    }();
    return suffixes;
}

const std::vector<std::string>& dosage_block_columns() {
    static const std::vector<std::string> columns = {
        "dosage_status",         "dosage_reason",      "dosage_filter",          "dosage_filter_reason", "repair_class",
        "dosage_repair_applied", "dosage_content_scope", "diplotype_scope",      "allocation_equivalence_scope",
        "serialization_representative", "original_allele1", "original_allele2", "original_content", "corrected_content",
        "admissible_content_vectors", "allocation_status", "original_GQ",       "gq_status",            "repair_separation",
        "conditional_repair_gap", "original_filter", "original_quality_filter", "original_genotype_support", "filter_basis",
        "quality_filter_basis"};
    return columns;
}

const std::vector<std::string>& locus_genotype_columns() {
    static const std::vector<std::string> columns = {
        "prefix", "run_status", "call_class", "no_call_reason", "locus_filter", "quality_rule", "quality_pass_blocks",
        "quality_lowqual_blocks", "filter_basis", "phase_support", "candidate_scope", "search_conditional", "reported_genotype_source",
        "hap1_id", "hap1_sha256", "hap1_bp", "hap1_block_path", "hap2_id", "hap2_sha256", "hap2_bp", "hap2_block_path",
        "unordered_genotype", "equivalence_members", "total_score", "fragment_term", "exposure_term", "sequence_pair_prior",
        "score_basis", "caller_total_score", "separation", "separation_basis", "separation_certified", "hap1_switches",
        "hap1_switch_blocks", "hap2_switches", "hap2_switch_blocks", "switch_basis", "dosage_mode", "dosage_audit_status",
        "dosage_complete_class_blocks", "dosage_complete_applied_blocks", "dosage_complete_applied_block_ids",
        "dosage_complete_not_applied_blocks", "dosage_complete_not_applied_detail", "dosage_content_corrected_blocks",
        "dosage_content_corrected_block_ids", "dosage_triggered_unrepaired_blocks", "dosage_triggered_unrepaired_block_ids",
        "dosage_suspect_blocks", "dosage_underpowered_blocks", "dosage_orientation_tied_block_ids", "locus_repair_scope",
        "block_table", "genotypes_table", "main_fasta", "original_fasta", "original_block_table", "dosage_blocks_table",
        "dosage_representative_fasta", "dosage_equivalence_table", "call_equivalence_set"};
    return columns;
}

std::optional<TemplateSwitches> template_switches(const PreparedGenotypeModel& model, const MosaicPath& path) {
    const PanelAlleleMatrix& panel = model.panel();
    TemplateSwitches out;
    std::vector<std::size_t> live;
    for (std::size_t b = 0; b < path.alleles.size(); ++b) {
        std::vector<std::size_t> carriers;
        for (std::size_t t = 0; t < panel.template_count(); ++t) {
            if (panel.template_alleles[t][b] == path.alleles[b]) carriers.push_back(t);
        }
        if (carriers.empty()) return std::nullopt;
        std::vector<std::size_t> kept;
        std::set_intersection(live.begin(), live.end(), carriers.begin(), carriers.end(), std::back_inserter(kept));
        if (b > 0 && kept.empty()) {
            ++out.count;
            out.blocks.push_back(b);
            live = carriers;
        } else {
            live = b == 0 ? carriers : kept;
        }
    }
    return out;
}

PreparedProposalEmissions rebuild_marker_emissions(const PreparedGenotypeModel& model, const ReportedFilterInputs& in) {
    const PreparedMarkerPanel panel = build_marker_panel(model, in.kmer_size, in.syncmer_s);
    const MarkerObservations observations = count_marker_observations(panel, *in.reads);
    const MarkerDepthEstimate depth = estimate_marker_depth(panel, observations, in.declared_marker_background);
    if (!depth.converged) throw std::runtime_error("marker depth fit did not converge");
    return build_proposal_emissions(panel, observations, proposal_marker_parameters(depth, in.kmer_size, in.syncmer_s));
}

std::vector<GuardRow> representation_guard_rows(const PreparedGenotypeModel& model, const PreparedProposalEmissions& emissions,
                                                const ReportedFilterInputs& in, const MosaicPath& first, const MosaicPath& second,
                                                const DiplotypeScore& called, double called_prior) {
    const std::size_t blocks = model.blocks().size();
    const std::vector<double> switches(blocks - 1, in.switch_probability);
    ScoringInputs si;
    si.model = &model; si.reads = in.reads; si.library = in.library; si.options = in.scoring;
    si.repair_max_physical_origins = in.repair_max_physical_origins; si.repair_max_state_classes = in.repair_max_state_classes;
    si.repair_max_piece_hits = in.repair_max_piece_hits; si.repair_max_latent_pairs = in.repair_max_latent_pairs;
    std::vector<GuardRow> out(blocks);
    for (std::size_t b = 0; b < blocks; ++b) {
        GuardRow& g = out[b];
        const BlockProposalEmissions& emission = emissions.blocks[b];
        const auto pair = std::minmax(first.alleles[b], second.alleles[b]);
        const std::pair<std::size_t, std::size_t> top = marker_best_pair(emission);
        g.marker_margin = emission.at(top.first, top.second) - emission.at(pair.first, pair.second);
        if (top == std::make_pair(pair.first, pair.second) || !(g.marker_margin >= kRepresentationBoundNats)) continue;
        g.competitor = std::to_string(top.first) + "/" + std::to_string(top.second);
        std::vector<std::pair<std::string, std::string>> diplotypes;
        for (int assignment = 0; assignment < 2; ++assignment) {
            if (assignment == 1 && top.first == top.second) break;
            MosaicPath a = first, c = second;
            a.alleles[b] = assignment == 0 ? top.first : top.second;
            c.alleles[b] = assignment == 0 ? top.second : top.first;
            diplotypes.push_back({spell_mosaic(model, a).bases, spell_mosaic(model, c).bases});
        }
        std::vector<RepresentationAssignment> assignments;
        try {
            const std::vector<DiplotypeScore> scores = score_diplotypes_exact(si, diplotypes);
            for (std::size_t i = 0; i < diplotypes.size(); ++i) {
                assignments.push_back({scores[i], log_unordered_sequence_pair_prior(model, diplotypes[i].first, diplotypes[i].second, switches)});
            }
        } catch (const std::exception& error) {
            g.filter = "NOT_EVALUATED";
            g.reason = std::string("competitor_not_scored: ") + error.what();
            continue;
        }
        const RepresentationEvidence e = representation_conflict(g.marker_margin, called, called_prior, assignments);
        g.fragment_diff = e.fragment_diff; g.exposure_diff = e.exposure_diff; g.likelihood_margin = e.likelihood_margin; g.statistic = e.statistic;
        if (e.conflict) {
            g.filter = "LOWQUAL_REPRESENTATION_CONFLICT";
            g.reason = "marker_and_exposure_prefer_competitor_fragment_term_prefers_call";
        }
    }
    return out;
}

std::vector<std::string> guard_cells(const GuardRow& g) {
    const auto number = [](double v) {
        if (std::isnan(v)) return std::string("NA");
        std::ostringstream text;
        text << std::setprecision(12) << v;
        return text.str();
    };
    return {g.filter, g.reason, g.reason == "guard_off" || g.filter == "NOT_EVALUATED_AFTER_DOSAGE_REPAIR" ? std::string("NA") : number(g.marker_margin),
            g.competitor, number(g.fragment_diff), number(g.exposure_diff), number(g.likelihood_margin), number(g.statistic)};
}

std::vector<bool> read_support(double fragment_term, double error_bound, const std::vector<Rival>& rivals, std::size_t blocks) {
    std::vector<bool> supported(blocks, true);
    for (const Rival& r : rivals) {
        if (r.fragment_term + r.error_bound < fragment_term - error_bound) continue;   // the reads alone put it below
        for (std::size_t b = 0; b < blocks; ++b) if (r.differs[b]) supported[b] = false;
    }
    return supported;
}

std::vector<RankedCandidate> read_caller_ranking(const std::string& prefix) {
    const Table scores = read_table(prefix + ".mosaic_scores.tsv"), candidates = read_table(prefix + ".mosaic_candidates.tsv");
    std::map<std::string, std::pair<std::string, std::string>> paths;
    for (const auto& r : candidates.rows) paths[r.at(candidates.col("candidate_id"))] = {r.at(candidates.col("hap1_alleles")), r.at(candidates.col("hap2_alleles"))};
    std::vector<RankedCandidate> out;
    for (const auto& r : scores.rows) {
        RankedCandidate c;
        c.fragment_term = std::stod(r.at(scores.col("fragment_term")));
        c.error_bound = std::stod(r.at(scores.col("numeric_error_bound")));
        c.winner = r.at(scores.col("in_winner_set")) == "true";
        const auto& p = paths.at(r.at(scores.col("representative_beam_candidate")));
        const std::vector<std::string> a = split(p.first, ','), b = split(p.second, ',');
        for (std::size_t k = 0; k < a.size(); ++k) c.genotype.push_back(std::minmax(std::stoul(a[k]), std::stoul(b[k])));
        out.push_back(c);
    }
    return out;
}

std::vector<std::pair<std::string, std::string>> locus_filter_summary(const std::vector<std::string>& block_ids,
                                                                      const std::vector<std::string>& filters,
                                                                      const std::vector<std::string>& supports,
                                                                      const std::vector<std::string>& qualities) {
    std::vector<std::string> lowqual, unevaluated, pass_ids, low_ids;
    std::size_t unsupported = 0, supported = 0;
    for (std::size_t b = 0; b < block_ids.size(); ++b) {
        if (filters[b] == "LOWQUAL_REPRESENTATION_CONFLICT") lowqual.push_back(block_ids[b]);
        else if (filters[b] != "PASS") unevaluated.push_back(block_ids[b]);
        (qualities[b] == "PASS" ? pass_ids : low_ids).push_back(block_ids[b]);
        if (supports[b] == "READS") ++supported; else ++unsupported;
    }
    return {{"locus_filter", !lowqual.empty() ? "LOWQUAL_REPRESENTATION_CONFLICT" : !unevaluated.empty() ? "NOT_EVALUATED" : "PASS"},
            {"lowqual_blocks", ids_or_none(lowqual)},
            {"not_evaluated_blocks", ids_or_none(unevaluated)},
            {"representation_guard_complete", unevaluated.empty() ? "true" : "false"},
            {"quality_pass_blocks", std::to_string(pass_ids.size())},
            {"quality_lowqual_blocks", std::to_string(low_ids.size())},
            {"quality_pass_block_ids", ids_or_none(pass_ids)},
            {"quality_lowqual_block_ids", ids_or_none(low_ids)},
            {"quality_reason_representation_conflict", std::to_string(lowqual.size())},
            {"quality_reason_representation_not_evaluated", std::to_string(unevaluated.size())},
            {"quality_reason_no_read_support", std::to_string(unsupported)},
            {"blocks_genotype_read_supported", std::to_string(supported)}};
}

DosageIntegration apply_dosage_mode(const PreparedGenotypeModel* model, const std::string& mode,
                                    const std::string& run_prefix, const std::string& out_prefix,
                                    const DosageShadowSummary* summary, const std::string& audit_status,
                                    const std::string& audit_note, const ReportedFilterInputs* filters) {
    if (mode != "off" && mode != "audit" && mode != "repair") throw std::runtime_error("dosage mode must be off, audit or repair");
    DosageIntegration d;
    d.mode = mode;
    d.summary = summary;
    d.audit_status = mode == "off" ? "NA" : audit_status;
    d.audit_note = audit_note;
    if (run_prefix != out_prefix) {
        for (const std::string& s : legacy_call_outputs()) {
            if (exists(run_prefix + s)) atomic_write(out_prefix + s, slurp(run_prefix + s));
        }
    }
    if (mode != "repair") return d;

    // Everything is decided and built in memory first; the originals are written, then the rewrites.
    const std::map<std::string, std::string> status = read_key_values(out_prefix + ".mosaic_status.tsv");
    const auto dosage = read_dosage_blocks(out_prefix + ".dosage_blocks.tsv");
    const bool stage_ran = summary != nullptr && audit_status == "COMPLETE";
    const std::optional<PersistedCall> call = read_persisted_call(out_prefix);

    // eligibility of each COMPLETE_DIPLOTYPE_REPAIR (section 22.2): fail closed on anything the frozen contract does
    // not cover exactly
    std::vector<const CompleteRepair*> eligible;
    if (stage_ran) {
        std::string refuse;
        if (value(status, "status") != "COMPLETE") refuse = "not_a_strict_call";
        else if (summary->complete_repairs.size() > 1) refuse = "several_complete_repairs_not_jointly_scored";
        else if (!call || !model) refuse = "reported_call_not_recovered";
        if (refuse.empty() && call && model) {
            const std::vector<FastaRecord> fa = exists(out_prefix + ".mosaic_haplotypes.fa") ? read_fasta(out_prefix + ".mosaic_haplotypes.fa")
                                                                                            : std::vector<FastaRecord>{};
            if (fa.size() != 2 || spell_mosaic(*model, call->first).bases != fa[0].sequence ||
                spell_mosaic(*model, call->second).bases != fa[1].sequence) {
                refuse = "reported_haplotypes_not_reproduced";
            }
        }
        for (const CompleteRepair& r : summary->complete_repairs) {
            if (!refuse.empty()) d.not_applied[r.block_id] = refuse;
            else if (!r.allele1 || !r.allele2) d.not_applied[r.block_id] = "repaired_allele_not_a_panel_allele";
            else eligible.push_back(&r);
        }
    }
    MosaicPath first, second;
    if (!eligible.empty()) {
        first = call->first;
        second = call->second;
        for (const CompleteRepair* r : eligible) {
            first.alleles.at(r->block) = *r->allele1;
            second.alleles.at(r->block) = *r->allele2;
            d.applied.push_back(r->block_id);
        }
        const CompleteRepair& r = *eligible.front();
        d.reported_first = first;
        d.reported_second = second;
        d.reported_fragment = r.fragment_term;
        d.reported_exposure = r.exposure_term;
        d.reported_prior = r.prior;
        d.reported_separation = r.repair_separation;
    }
    const std::set<std::string> applied(d.applied.begin(), d.applied.end());

    // The reported genotype's own filters (section 23). Only a run whose genotype a repair changed is re-annotated; the
    // guard and the read-support rule must first reproduce the caller's own annotation of its own call, or the reported
    // genotype's diagnostics are NOT_EVALUATED_AFTER_DOSAGE_REPAIR. An annotation only: nothing here feeds a decision.
    static const std::string kNotEvaluated = "NOT_EVALUATED_AFTER_DOSAGE_REPAIR";
    struct Reported { std::vector<std::string> guard; std::string support, quality, quality_reason; };
    std::map<std::string, Reported> reported;
    std::map<std::string, std::map<std::string, std::string>> caller;   // block_genotypes.tsv rows, as the caller wrote them
    if (exists(out_prefix + ".block_genotypes.tsv")) {
        const Table t = read_table(out_prefix + ".block_genotypes.tsv");
        for (const auto& r : t.rows) {
            std::map<std::string, std::string> row;
            for (std::size_t i = 0; i < t.header.size() && i < r.size(); ++i) row[t.header[i]] = r[i];
            caller[row.at("block_id")] = row;
        }
    }
    if (!applied.empty()) {
        const std::size_t nb = model->blocks().size();
        const CompleteRepair& r = *eligible.front();
        std::vector<GuardRow> rows(nb);
        std::string guard_note = "none", support_note = "none";
        const auto guard_not_evaluated = [&](const std::string& why) {
            for (GuardRow& g : rows) { g = GuardRow{}; g.filter = kNotEvaluated; g.reason = why; }
            guard_note = why;
        };
        if (!filters || !filters->reads || !filters->library) {
            guard_not_evaluated("inputs_unavailable");
        } else if (!filters->representation_guard) {
            for (GuardRow& g : rows) { g = GuardRow{}; g.filter = "NOT_EVALUATED"; g.reason = "guard_off"; }
        } else {
            try {
                const PreparedProposalEmissions emissions = rebuild_marker_emissions(*model, *filters);
                const DiplotypeScore original{call->score.fragment + call->score.exposure, call->score.fragment, call->score.exposure, 0.0, 0};
                const std::vector<GuardRow> again = representation_guard_rows(*model, emissions, *filters, call->first, call->second, original, call->prior);
                std::string mismatch;
                for (std::size_t b = 0; b < nb && mismatch.empty(); ++b) {
                    const auto& row = caller.at(model->blocks()[b].block_id);
                    const std::vector<std::string> cells = guard_cells(again[b]);
                    bool same = cells[0] == value(row, "filter") && cells[2] == value(row, "representation_marker_margin") &&
                                cells[3] == value(row, "representation_competitor_pair");
                    const double want = std::atof(value(row, "representation_conflict_statistic").c_str());
                    if (same && cells[7] != "NA") same = std::fabs(again[b].statistic - want) <= 1e-4 * std::max(1.0, std::fabs(want));
                    if (same && cells[7] == "NA") same = value(row, "representation_conflict_statistic") == "NA";
                    if (!same) mismatch = model->blocks()[b].block_id;
                }
                if (!mismatch.empty()) guard_not_evaluated("guard_not_reproduced_at_" + mismatch);
                else rows = representation_guard_rows(*model, emissions, *filters, first, second, DiplotypeScore{r.likelihood, r.fragment_term, r.exposure_term, 0.0, 0}, r.prior);
            } catch (const std::exception& error) {
                guard_not_evaluated(std::string("guard_failed: ") + error.what());
            }
        }
        std::vector<std::string> support(nb, kNotEvaluated);
        try {
            const std::vector<RankedCandidate> ranking = read_caller_ranking(out_prefix);
            if (ranking.empty() || !ranking.front().winner) throw std::runtime_error("no ranked winner");
            const auto genotype_of = [&](const MosaicPath& x, const MosaicPath& y) {
                std::vector<std::pair<std::size_t, std::size_t>> g;
                for (std::size_t b = 0; b < nb; ++b) g.push_back(std::minmax(x.alleles[b], y.alleles[b]));
                return g;
            };
            const auto rivals_of = [&](const std::vector<std::pair<std::size_t, std::size_t>>& gt, std::size_t from) {
                std::vector<Rival> out;
                for (std::size_t i = from; i < ranking.size(); ++i) {
                    Rival v{ranking[i].fragment_term, ranking[i].error_bound, std::vector<bool>(nb, false)};
                    for (std::size_t b = 0; b < nb; ++b) v.differs[b] = ranking[i].genotype[b] != gt[b];
                    out.push_back(v);
                }
                return out;
            };
            const std::vector<bool> again = read_support(ranking[0].fragment_term, ranking[0].error_bound,
                                                         rivals_of(genotype_of(call->first, call->second), 1), nb);
            std::string mismatch;
            for (std::size_t b = 0; b < nb && mismatch.empty(); ++b) {
                if ((again[b] ? "READS" : "PRIOR_ONLY") != value(caller.at(model->blocks()[b].block_id), "genotype_support")) mismatch = model->blocks()[b].block_id;
            }
            if (!mismatch.empty()) {
                support_note = "read_support_not_reproduced_at_" + mismatch;
            } else {
                // every caller candidate (the original winner included) and every other scored repair candidate is a
                // rival of the reported repaired diplotype; the repaired diplotype carries the caller winner's bound
                std::vector<Rival> rivals = rivals_of(genotype_of(first, second), 0);
                const auto unordered_seq = [](const std::string& x, const std::string& y) { return x <= y ? std::make_pair(x, y) : std::make_pair(y, x); };
                for (const CompleteRepair::Alternative& alt : r.alternatives) {
                    Rival v{alt.fragment_term, ranking[0].error_bound, std::vector<bool>(nb, false)};
                    v.differs[r.block] = unordered_seq(alt.seq1, alt.seq2) != unordered_seq(r.seq1, r.seq2);
                    rivals.push_back(v);
                }
                const std::vector<bool> s = read_support(r.fragment_term, ranking[0].error_bound, rivals, nb);
                for (std::size_t b = 0; b < nb; ++b) support[b] = s[b] ? "READS" : "PRIOR_ONLY";
            }
        } catch (const std::exception& error) {
            support_note = std::string("read_support_failed: ") + error.what();
        }
        for (std::size_t b = 0; b < nb; ++b) {
            Reported x;
            x.guard = guard_cells(rows[b]);
            x.support = support[b];
            if (rows[b].filter == kNotEvaluated || support[b] == kNotEvaluated) {
                x.quality = kNotEvaluated;
                std::string why;
                if (rows[b].filter == kNotEvaluated) why = "representation_not_evaluated_after_dosage_repair";
                if (support[b] == kNotEvaluated) why += (why.empty() ? "" : ",") + std::string("read_support_not_evaluated_after_dosage_repair");
                x.quality_reason = why;
            } else {
                const BlockQuality q = block_quality(rows[b].filter, support[b]);
                x.quality = q.filter; x.quality_reason = q.reason;
            }
            reported[model->blocks()[b].block_id] = x;
        }
        const bool complete = guard_note == "none" && support_note == "none";
        d.filter_basis = complete ? "recomputed_for_reported_genotype" : "not_evaluated_after_dosage_repair";
        d.filter_note = complete ? "none" : (guard_note != "none" ? guard_note : "") + (guard_note != "none" && support_note != "none" ? ";" : "") +
                                                (support_note != "none" ? support_note : "");
    }

    std::map<std::string, std::string> writes;   // path -> content, written only after every one is built
    const auto block_of = [&](const std::string& id) -> std::size_t {
        for (std::size_t b = 0; model && b < model->blocks().size(); ++b) {
            if (model->blocks()[b].block_id == id) return b;
        }
        throw std::runtime_error("unknown block " + id);
    };
    const auto rewrite_block_table = [&](const std::string& suffix, bool with_sources) {
        const std::string path = out_prefix + suffix;
        if (!exists(path)) return;
        Table t = read_table(path);
        const std::size_t c_id = t.col("block_id"), c_a1 = t.col("allele1"), c_a2 = t.col("allele2"), c_gq = t.col("conditional_GQ");
        for (auto& row : t.rows) {
            const std::string id = row.at(c_id);
            const std::string caller1 = row.at(c_a1), caller2 = row.at(c_a2), caller_gq = row.at(c_gq);
            const bool is_applied = applied.count(id) != 0;
            if (is_applied) {
                const std::size_t b = block_of(id);
                const PreparedGenotypeBlock& block = model->blocks()[b];
                const std::size_t a1 = first.alleles[b], a2 = second.alleles[b];
                row.at(c_a1) = std::to_string(a1);
                row.at(c_a2) = std::to_string(a2);
                row.at(t.col("allele1_bp")) = std::to_string(block.alleles[a1].sequence.size());
                row.at(t.col("allele2_bp")) = std::to_string(block.alleles[a2].sequence.size());
                row.at(t.col("zygosity")) = a1 == a2 ? "HOM" : "HET";
                row.at(t.col("conditional_genotype_posterior")) = "NA";
                row.at(c_gq) = "NA";
                if (with_sources) {
                    for (const auto& [col, a] : {std::make_pair("allele1_sources", a1), std::make_pair("allele2_sources", a2)}) {
                        row.at(t.col(col)) = join(block.alleles[a].source_ids, ",");
                    }
                }
            }
            const auto found = dosage.find(id);
            const std::map<std::string, std::string> none;
            const std::map<std::string, std::string>& g = found == dosage.end() ? none : found->second;
            const auto cell = [&](const std::string& k) { return value(g, k); };
            const std::string st = stage_ran ? cell("dosage_status") : "NA";
            const std::string reason = stage_ran ? cell("dosage_reason") : "NA";
            const auto na = d.not_applied.find(id);
            const auto [filter, filter_reason] =
                stage_ran ? dosage_filter(st, reason, is_applied, na == d.not_applied.end() ? "none" : na->second)
                          : std::make_pair(std::string("DOSAGE_NOT_EVALUATED"),
                                           std::string(summary == nullptr ? "dosage_stage_not_run" : "dosage_stage_" + audit_status));
            const std::string cls = stage_ran ? cell("repair_class") : "NA";
            std::vector<std::string> extra = {
                st, reason, filter, filter_reason, cls,
                cls == "COMPLETE_DIPLOTYPE_REPAIR" ? (is_applied ? "true" : "false") : "NA",
                cell("dosage_content_scope"), cell("diplotype_scope"), cell("allocation_equivalence_scope"),
                cell("serialization_representative"), caller1, caller2, cell("original_content"), cell("corrected_content"),
                cell("admissible_content_vectors"), cell("allocation_status"), caller_gq,
                is_applied ? "UNAVAILABLE_REPAIRED" : "UNCHANGED_CONDITIONAL", cell("repair_separation"),
                cell("conditional_repair_gap")};
            // the reported genotype's filters: the caller's where no repair changed the genotype, recomputed otherwise
            const std::string caller_filter = t.has("filter") ? row.at(t.col("filter")) : "NA";
            const std::string caller_quality = t.has("quality_filter") ? row.at(t.col("quality_filter")) : "NA";
            const std::string caller_support = caller.count(id) ? value(caller.at(id), "genotype_support") : "NA";
            std::string filter_basis = "caller", quality_basis = "caller";
            const auto rep = reported.find(id);
            if (rep != reported.end()) {
                static const char* guard_columns[] = {"filter", "filter_reason", "representation_marker_margin", "representation_competitor_pair",
                                                      "representation_fragment_diff", "representation_exposure_diff",
                                                      "representation_likelihood_margin", "representation_conflict_statistic"};
                for (std::size_t k = 0; k < 8; ++k) if (t.has(guard_columns[k])) row.at(t.col(guard_columns[k])) = rep->second.guard[k];
                if (t.has("genotype_support")) row.at(t.col("genotype_support")) = rep->second.support;
                if (t.has("quality_filter")) row.at(t.col("quality_filter")) = rep->second.quality;
                if (t.has("quality_reason")) row.at(t.col("quality_reason")) = rep->second.quality_reason;
                filter_basis = rep->second.guard[0] == kNotEvaluated ? "not_evaluated_after_dosage_repair"
                             : rep->second.guard[1] == "guard_off" ? "guard_off" : "recomputed_for_reported_genotype";
                quality_basis = rep->second.quality == kNotEvaluated ? "not_evaluated_after_dosage_repair" : "recomputed_for_reported_genotype";
            }
            for (const std::string& v : {caller_filter, caller_quality, caller_support, filter_basis, quality_basis}) row.push_back(v);
            row.insert(row.end() - 5, extra.begin(), extra.end());
        }
        t.header.insert(t.header.end(), dosage_block_columns().begin(), dosage_block_columns().end());
        writes[path] = t.text();
    };
    rewrite_block_table(".genotypes.tsv", true);
    rewrite_block_table(".block_genotypes.tsv", false);
    // locus-level summaries of the reported genotype, recomputed from the reported block rows exactly as the caller
    // derives them from its own (locus filter fails closed; quality counts over every block)
    std::vector<std::pair<std::string, std::string>> summaries;   // status key -> reported value, caller's order
    if (!reported.empty()) {
        std::vector<std::string> ids, fs, ss, qs;
        for (const PreparedGenotypeBlock& blk : model->blocks()) {
            const Reported& x = reported.at(blk.block_id);
            ids.push_back(blk.block_id); fs.push_back(x.guard[0]); ss.push_back(x.support); qs.push_back(x.quality);
        }
        summaries = locus_filter_summary(ids, fs, ss, qs);
    }
    const auto summary_of = [&](const std::string& key) {
        for (const auto& [k, v] : summaries) if (k == key) return v;
        return std::string();
    };
    if (!applied.empty()) {
        std::vector<FastaRecord> fa = read_fasta(out_prefix + ".mosaic_haplotypes.fa");
        std::ostringstream o;
        for (std::size_t h = 0; h < 2; ++h) {
            const std::string bases = spell_mosaic(*model, h == 0 ? first : second).bases;
            std::vector<std::string> tokens = split(fa[h].header, ' ');
            for (std::string& tok : tokens) {
                if (tok.rfind("bp=", 0) == 0) tok = "bp=" + std::to_string(bases.size());
                if (tok.rfind("locus_filter=", 0) == 0 && !summaries.empty()) tok = "locus_filter=" + summary_of("locus_filter");
            }
            o << '>' << join(tokens, " ") << " dosage_repaired_blocks=" << join(d.applied, ",") << " dosage_repair_role=RESOLVED\n";
            for (std::size_t at = 0; at < bases.size(); at += 60) o << bases.substr(at, 60) << '\n';
        }
        writes[out_prefix + ".mosaic_haplotypes.fa"] = o.str();
    }
    {
        std::vector<std::string> conditional, triggered;
        for (const auto& [id, g] : dosage) {
            if (!stage_ran) break;
            const std::string st = value(g, "dosage_status");
            if (st == "CONTENT_REPAIRED") conditional.push_back(id);
            if (st == "TRIGGERED_UNREPAIRED") triggered.push_back(id);
        }
        std::vector<std::string> refused;
        for (const auto& [id, why] : d.not_applied) refused.push_back(id + ":" + why);
        std::ostringstream o;
        o << "dosage_mode\trepair\ndosage_policy\tdosage_policy_v2\ndosage_audit_status\t" << audit_status << '\n';
        if (summary) {
            o << "dosage_molecule_span_bp\t" << summary->molecule_span_bp << "\ndosage_m\t" << summary->m << '\n';
            for (const auto& [k, v] : summary->repair_classes) o << "dosage_repair_class_" << k << '\t' << v << '\n';
        }
        o << "dosage_repaired_blocks\t" << ids_or_none(d.applied) << "\ndosage_repair_not_applied_blocks\t" << ids_or_none(refused)
          << "\ndosage_content_corrected_blocks\t" << ids_or_none(conditional) << "\ndosage_triggered_unrepaired_blocks\t"
          << ids_or_none(triggered) << "\ndosage_reported_genotype\t" << (applied.empty() ? "caller" : "caller_with_dosage_repair")
          << "\ndosage_search_winner_superseded\t" << (applied.empty() ? "false" : "true") << '\n'
          << "dosage_filters_basis\t" << d.filter_basis << "\ndosage_filters_note\t" << d.filter_note << '\n';
        const std::string path = out_prefix + ".mosaic_status.tsv";
        if (exists(path)) {
            // the caller's filter-derived lines describe the reported genotype; the caller's values follow as original_<key>
            std::string text = slurp(path), originals;
            if (!summaries.empty()) {
                std::vector<std::string> lines = split(text.substr(0, text.size() - (text.empty() || text.back() != '\n' ? 0 : 1)), '\n');
                for (std::string& line : lines) {
                    const std::size_t tab = line.find('\t');
                    if (tab == std::string::npos) continue;
                    const std::string key = line.substr(0, tab), replacement = summary_of(key);
                    if (replacement.empty()) continue;
                    originals += "original_" + line + "\n";
                    line = key + "\t" + replacement;
                }
                text = join(lines, "\n") + "\n";
            }
            writes[path] = text + o.str() + originals;
        }
    }
    for (const std::string& s : legacy_call_outputs()) {
        if (exists(out_prefix + s)) atomic_write(out_prefix + ".original" + s, slurp(out_prefix + s));
    }
    d.originals_written = true;
    for (const auto& [path, content] : writes) atomic_write(path, content);
    return d;
}

void write_locus_genotype(const PreparedGenotypeModel* model, const std::string& out_prefix, const DosageIntegration& d) {
    const std::map<std::string, std::string> status = read_key_values(out_prefix + ".mosaic_status.tsv");
    const std::string base = basename_of(out_prefix);
    const auto link = [&](const std::string& suffix) { return exists(out_prefix + suffix) ? base + suffix : std::string("NA"); };
    std::map<std::string, std::string> v;
    for (const std::string& c : locus_genotype_columns()) v[c] = "NA";
    v["prefix"] = base;
    v["run_status"] = value(status, "status");
    v["call_class"] = value(status, "call_class");
    v["no_call_reason"] = value(status, "no_call_reason");
    v["locus_filter"] = value(status, "locus_filter");
    v["quality_rule"] = value(status, "quality_rule");
    v["quality_pass_blocks"] = value(status, "quality_pass_blocks");
    v["quality_lowqual_blocks"] = value(status, "quality_lowqual_blocks");
    v["filter_basis"] = status.count("dosage_filters_basis") ? value(status, "dosage_filters_basis") : "caller";
    v["phase_support"] = value(status, "phase_support");
    v["candidate_scope"] = status.count("final_candidate_scope") ? value(status, "final_candidate_scope") : value(status, "candidate_scope");
    v["search_conditional"] = value(status, "search_conditional");
    v["dosage_mode"] = d.mode;
    v["dosage_audit_status"] = d.audit_status;
    const bool strict = v["call_class"] == "strict_call";
    const bool equivalent = v["call_class"] == "genotype_equivalent_call";
    const bool repaired = d.reported_first.has_value();
    v["reported_genotype_source"] = strict || equivalent ? (repaired ? "caller_with_dosage_repair" : "caller") : "none";

    if (strict || equivalent) {
        // the caller's call: the rank-1 representative of the equivalence set, its scores as the caller wrote them
        const Table eq = read_table(out_prefix + ".equivalence_set.tsv");
        std::size_t members = 0;
        for (const auto& r : eq.rows) {
            if (r.at(eq.col("rank")) != "1") continue;
            ++members;
            if (r.at(eq.col("representative")) != "true") continue;
            v["total_score"] = r.at(eq.col("score"));
            v["fragment_term"] = r.at(eq.col("fragment_term"));
            v["exposure_term"] = r.at(eq.col("exposure_term"));
            v["sequence_pair_prior"] = r.at(eq.col("sequence_pair_prior"));
            v["hap1_block_path"] = r.at(eq.col("hap1_alleles"));
            v["hap2_block_path"] = r.at(eq.col("hap2_alleles"));
        }
        v["equivalence_members"] = std::to_string(members);
        v["caller_total_score"] = v["total_score"];
        v["score_basis"] = "caller_winner";
        v["separation"] = value(status, "genotype_separation_gap");
        v["separation_basis"] = "caller_genotype_separation_gap";
        v["separation_certified"] = value(status, "genotype_separation_certified");
        if (repaired) {
            v["hap1_block_path"] = join_path(*d.reported_first);
            v["hap2_block_path"] = join_path(*d.reported_second);
            v["fragment_term"] = format(d.reported_fragment);
            v["exposure_term"] = format(d.reported_exposure);
            v["sequence_pair_prior"] = format(d.reported_prior);
            v["total_score"] = format(d.reported_fragment + d.reported_exposure + d.reported_prior);
            v["score_basis"] = "dosage_repaired_block_rescored";
            v["separation"] = format(d.reported_separation);
            v["separation_basis"] = "dosage_repair_separation";
            v["separation_certified"] = "NA";
        }
        // the unordered genotype at every block, from the reported block table
        const Table bg = read_table(out_prefix + ".block_genotypes.tsv");
        std::vector<std::string> gt;
        for (const auto& r : bg.rows) {
            std::size_t a = std::stoul(r.at(bg.col("allele1"))), b = std::stoul(r.at(bg.col("allele2")));
            if (a > b) std::swap(a, b);
            gt.push_back(r.at(bg.col("block_id")) + "=" + std::to_string(a) + "/" + std::to_string(b));
        }
        v["unordered_genotype"] = join(gt, ";");
    }
    if (strict && exists(out_prefix + ".mosaic_haplotypes.fa")) {
        const std::vector<FastaRecord> fa = read_fasta(out_prefix + ".mosaic_haplotypes.fa");
        for (std::size_t h = 0; h < 2 && h < fa.size(); ++h) {
            const std::string hap = h == 0 ? "hap1" : "hap2";
            v[hap + "_id"] = split(fa[h].header, ' ').front();
            v[hap + "_sha256"] = sha256_bytes(fa[h].sequence);
            v[hap + "_bp"] = std::to_string(fa[h].sequence.size());
        }
        if (model) {
            v["switch_basis"] = "panel_templates_min_switches_left_greedy";
            for (const std::string hap : {"hap1", "hap2"}) {
                MosaicPath p;
                for (const std::string& a : split(v[hap + "_block_path"], ',')) p.alleles.push_back(std::stoul(a));
                const std::optional<TemplateSwitches> s = template_switches(*model, p);
                if (!s) { v["switch_basis"] = "allele_not_carried_by_any_template"; continue; }
                v[hap + "_switches"] = std::to_string(s->count);
                std::vector<std::string> at;
                for (std::size_t b : s->blocks) at.push_back(model->blocks()[b].block_id);
                v[hap + "_switch_blocks"] = ids_or_none(at);
            }
        }
    } else if (equivalent) {
        // genotype-equivalent: the unordered genotype and its equivalence set only; no phase, no locus haplotypes
        v["hap1_block_path"] = v["hap2_block_path"] = "NA";
        v["switch_basis"] = "phase_not_claimed";
    }

    // dosage counts (NA when the stage did not run)
    const auto dosage = read_dosage_blocks(out_prefix + ".dosage_blocks.tsv");
    if (d.mode != "off" && d.audit_status == "COMPLETE") {
        std::size_t complete = 0, suspect = 0, underpowered = 0, triggered_any = 0;
        std::vector<std::string> conditional, triggered, tied;
        for (const auto& [id, g] : dosage) {
            const std::string st = value(g, "dosage_status"), cls = value(g, "repair_class");
            if (cls == "COMPLETE_DIPLOTYPE_REPAIR") ++complete;
            if (st == "CONTENT_REPAIRED") conditional.push_back(id);
            if (st == "TRIGGERED_UNREPAIRED" || st == "AUTO_REPAIR") triggered.push_back(id);
            if (st == "SUSPECT") ++suspect;
            if (st == "UNDERPOWERED") ++underpowered;
            if (st == "REPAIRED" || st == "CONTENT_REPAIRED" || st == "TRIGGERED_UNREPAIRED" || st == "AUTO_REPAIR") ++triggered_any;
            if (cls == "COMPLETE_DIPLOTYPE_REPAIR" && value(g, "orientation_tie") == "true") tied.push_back(id);
        }
        std::vector<std::string> refused;
        for (const auto& [id, why] : d.not_applied) refused.push_back(id + ":" + why);
        v["dosage_complete_class_blocks"] = std::to_string(complete);
        v["dosage_complete_applied_blocks"] = std::to_string(d.applied.size());
        v["dosage_complete_applied_block_ids"] = ids_or_none(d.applied);
        v["dosage_complete_not_applied_blocks"] = std::to_string(d.not_applied.size());
        v["dosage_complete_not_applied_detail"] = ids_or_none(refused);
        v["dosage_content_corrected_blocks"] = std::to_string(conditional.size());
        v["dosage_content_corrected_block_ids"] = ids_or_none(conditional);
        v["dosage_triggered_unrepaired_blocks"] = std::to_string(triggered.size());
        v["dosage_triggered_unrepaired_block_ids"] = ids_or_none(triggered);
        v["dosage_suspect_blocks"] = std::to_string(suspect);
        v["dosage_underpowered_blocks"] = std::to_string(underpowered);
        v["dosage_orientation_tied_block_ids"] = ids_or_none(tied);
        // whole-locus repair scope: COMPLETE only when every triggered block has a complete diplotype repair (applied,
        // in repair mode)
        const bool all_complete = triggered_any == complete && (d.mode != "repair" || d.applied.size() == complete);
        v["locus_repair_scope"] = triggered_any == 0 ? "NONE_TRIGGERED" : all_complete ? "COMPLETE" : "CONDITIONAL";
    } else {
        v["locus_repair_scope"] = "NOT_RUN";
    }
    v["block_table"] = link(".block_genotypes.tsv");
    v["genotypes_table"] = link(".genotypes.tsv");
    v["main_fasta"] = strict ? link(".mosaic_haplotypes.fa") : "NA";
    v["original_fasta"] = strict ? link(".original.mosaic_haplotypes.fa") : "NA";
    v["original_block_table"] = link(".original.block_genotypes.tsv");
    v["dosage_blocks_table"] = link(".dosage_blocks.tsv");
    v["dosage_representative_fasta"] = link(".dosage_repaired.fa");
    v["dosage_equivalence_table"] = link(".dosage_allocation_equivalence.tsv");
    v["call_equivalence_set"] = link(".equivalence_set.tsv");

    std::vector<std::string> row;
    for (const std::string& c : locus_genotype_columns()) row.push_back(v.at(c));
    atomic_write(out_prefix + ".locus_genotype.tsv", join(locus_genotype_columns(), "\t") + "\n" + join(row, "\t") + "\n");
}

} // namespace panvar::dosage
