#include "panvar/genotype_dosage_repair.hpp"

#include "panvar/genotype_repair.hpp"
#include "panvar/ls_hmm.hpp"
#include "panvar/mosaic_spelling.hpp"

#include <algorithm>
#include <chrono>
#include <set>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace panvar::dosage {

namespace {

std::vector<std::string> split_tabs(const std::string& line) {
    std::vector<std::string> f;
    std::stringstream row(line);
    for (std::string c; std::getline(row, c, '\t');) f.push_back(c);
    return f;
}

std::vector<std::size_t> parse_path(const std::string& joined) {
    std::vector<std::size_t> out;
    std::stringstream s(joined);
    for (std::string x; std::getline(s, x, ',');) out.push_back(std::stoul(x));
    return out;
}

// block_id -> the named columns of the run's block_genotypes.tsv
std::map<std::string, std::map<std::string, std::string>> read_block_genotypes(const std::string& path) {
    std::map<std::string, std::map<std::string, std::string>> out;
    std::ifstream in(path);
    if (!in) return out;
    std::string line;
    if (!std::getline(in, line)) return out;
    const std::vector<std::string> header = split_tabs(line);
    while (std::getline(in, line)) {
        const std::vector<std::string> f = split_tabs(line);
        std::map<std::string, std::string> row;
        for (std::size_t i = 0; i < header.size() && i < f.size(); ++i) row[header[i]] = f[i];
        if (row.count("block_id")) out[row["block_id"]] = row;
    }
    return out;
}

double seconds_since(std::chrono::steady_clock::time_point t) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
}

std::string field(const std::map<std::string, std::string>& row, const std::string& key) {
    const auto it = row.find(key);
    return it == row.end() ? "NA" : it->second;
}

} // namespace

std::optional<PersistedCall> read_persisted_call(const std::string& run_prefix) {
    const std::string path = run_prefix + ".equivalence_set.tsv";
    std::ifstream in(path);
    if (!in) return std::nullopt;
    std::string line;
    if (!std::getline(in, line)) return std::nullopt;
    const std::vector<std::string> header = split_tabs(line);
    const auto col = [&](const std::string& name) {
        const auto it = std::find(header.begin(), header.end(), name);
        if (it == header.end()) throw std::runtime_error(path + " has no column " + name);
        return static_cast<std::size_t>(it - header.begin());
    };
    const std::size_t c_rank = col("rank"), c_rep = col("representative"), c_h1 = col("hap1_alleles"), c_h2 = col("hap2_alleles"),
                      c_tot = col("score"), c_frag = col("fragment_term"), c_expo = col("exposure_term"),
                      c_prior = col("sequence_pair_prior");
    while (std::getline(in, line)) {
        const std::vector<std::string> f = split_tabs(line);
        if (f.size() <= std::max({c_rank, c_rep, c_h1, c_h2, c_tot, c_frag, c_expo, c_prior})) continue;
        if (f[c_rank] != "1" || f[c_rep] != "true") continue;
        PersistedCall call;
        call.hap1 = f[c_h1];
        call.hap2 = f[c_h2];
        call.first = MosaicPath{parse_path(call.hap1)};
        call.second = MosaicPath{parse_path(call.hap2)};
        call.score = {std::stod(f[c_tot]), std::stod(f[c_frag]), std::stod(f[c_expo])};
        call.prior = std::stod(f[c_prior]);
        return call;
    }
    return std::nullopt;
}

DosageShadowSummary run_dosage_shadow(const PreparedGenotypeModel& model, const NodeSequences& nodes,
                                      const std::vector<ReadPair>& reads, const LibraryModel& library,
                                      const FragmentScoringOptions& scoring, const std::vector<double>& switch_probability,
                                      const std::string& run_prefix, const std::string& out_prefix,
                                      const DosageShadowOptions& options, const std::optional<PersistedCall>& call_override) {
    const auto started = std::chrono::steady_clock::now();
    DosageShadowSummary summary;
    const DosageSettings settings;
    const auto write_status = [&]() {
        std::ostringstream o;
        o << std::setprecision(17) << "key\tvalue\n"
          << "dosage_audit\tshadow\n"
          << "dosage_audit_status\t" << summary.audit_status << '\n'
          << "dosage_audit_note\t" << (summary.note.empty() ? "none" : summary.note) << '\n'
          << "dosage_policy\tdosage_policy_v2\n"
          << "dosage_alpha\t" << settings.alpha << '\n'
          << "dosage_molecule_span_bp\t" << summary.molecule_span_bp << '\n'
          << "dosage_census_sensitive_blocks\t" << summary.sensitive_blocks << '\n'
          << "dosage_m\t" << summary.m << '\n'
          << "dosage_repair\t" << (options.repair ? "on" : "off") << '\n'
          << "dosage_repairs\t" << summary.repairs << '\n'
          << "dosage_repairs_conditional\t" << summary.repairs_conditional << '\n'
          << "dosage_allocation_test\t" << (options.allocation_test ? "on" : "off") << '\n'
          << "dosage_allocation_blocks_tested\t" << summary.allocation_blocks_tested << '\n';
        for (const auto& [k, v] : summary.status_counts) o << "dosage_status_" << k << '\t' << v << '\n';
        for (const auto& [k, v] : summary.reason_counts) o << "dosage_reason_" << k << '\t' << v << '\n';
        summary.seconds_total = seconds_since(started);
        o << "seconds_dosage_census\t" << summary.seconds_census << '\n'
          << "seconds_dosage_features\t" << summary.seconds_features << '\n'
          << "seconds_dosage_estimate\t" << summary.seconds_estimate << '\n'
          << "seconds_dosage_repair\t" << summary.seconds_repair << '\n'
          << "seconds_dosage_provenance\t" << summary.seconds_provenance << '\n'
          << "seconds_dosage_allocation\t" << summary.seconds_allocation << '\n'
          << "seconds_dosage_total\t" << summary.seconds_total << '\n';
        atomic_write(out_prefix + ".dosage_status.tsv", o.str());
    };
    const std::optional<std::size_t> span = molecule_span_bp(library.insert_mean_bp, library.insert_sd_bp, settings);
    if (!span) {
        summary.audit_status = "REFUSED_LIBRARY";
        summary.note = "the library declares no short-read insert distribution (unsupported_library)";
        std::ostringstream o;
        o << "block_id\tdosage_status\tdosage_reason\n";
        for (const PreparedGenotypeBlock& b : model.blocks()) {
            o << b.block_id << "\tNOT_EVALUATED\tunsupported_library\n";
            ++summary.status_counts["NOT_EVALUATED"];
            ++summary.reason_counts["unsupported_library"];
        }
        atomic_write(out_prefix + ".dosage_blocks.tsv", o.str());
        write_status();
        return summary;
    }
    summary.molecule_span_bp = *span;
    const std::optional<PersistedCall> persisted = call_override ? call_override : read_persisted_call(run_prefix);
    if (!persisted) {
        summary.audit_status = "NOT_RUN";
        summary.note = "the run reports no rank-1 representative (no call to audit)";
        write_status();
        return summary;
    }
    auto t = std::chrono::steady_clock::now();
    const std::vector<CensusBlock> census = run_census(model, nodes, *span, settings);
    summary.seconds_census = seconds_since(t);
    t = std::chrono::steady_clock::now();
    std::vector<BlockFeatureModel> features(model.blocks().size());
    for (const CensusBlock& c : census) {
        if (!c.sensitive()) continue;
        ++summary.sensitive_blocks;
        features[c.block] = pooled_features(model.blocks()[c.block], settings);
    }
    summary.seconds_features = seconds_since(t);
    const ReportedCall call = reported_call(model, persisted->first, persisted->second);
    t = std::chrono::steady_clock::now();
    DosageAudit audit = audit_dosage(model, census, features, call, reads, *span, settings, options.threads);
    summary.seconds_estimate = seconds_since(t);
    summary.m = audit.m;

    std::vector<std::pair<std::string, RepairResult>> repairs;
    std::map<std::size_t, std::size_t> repair_of;   // block -> index into repairs
    if (options.repair) {
        t = std::chrono::steady_clock::now();
        ScoringInputs in;
        in.model = &model; in.first = persisted->first; in.second = persisted->second; in.call = &call; in.reads = &reads;
        in.library = &library; in.options = scoring; in.options.threads = options.threads; in.switch_probability = switch_probability;
        in.persisted = persisted->score;
        RepairSettings rs;
        rs.chunk_haplotypes = options.chunk_haplotypes;
        for (BlockAudit& ba : audit.blocks) {
            if (ba.status != DosageStatus::AutoRepair) continue;
            repair_of[ba.block] = repairs.size();
            repairs.push_back({ba.block_id, repair_block(ba, census[ba.block], features[ba.block], nodes, in, settings, rs)});
            const RepairResult& r = repairs.back().second;
            ba.status = r.outcome == "REPAIRED" ? DosageStatus::Repaired : DosageStatus::TriggeredUnrepaired;
            ba.reason = r.reason;
            if (r.outcome == "REPAIRED") { ++summary.repairs; if (r.scope == "CONDITIONAL") ++summary.repairs_conditional; }
        }
        summary.seconds_repair = seconds_since(t);
    }

    // allocation (spec section 8.3 and Sol's 2026-09-28 semantics): every repeated-unit block's reported units, the
    // observability test when requested, and the final status -- UNRESOLVED when the test proves it or when several
    // allocations remain tied, NOT_TESTED when an eligible test was not run, NOT_APPLICABLE only where allocation does
    // not apply (pair route, or fewer than three units)
    std::map<std::size_t, AllocationTestResult> allocation;
    std::map<std::size_t, std::pair<std::string, std::string>> allocation_status;   // block -> (status, basis)
    {
        const auto ta = std::chrono::steady_clock::now();
        ScoringInputs in;
        in.model = &model; in.first = persisted->first; in.second = persisted->second; in.call = &call; in.reads = &reads;
        in.library = &library; in.options = scoring; in.options.threads = options.threads; in.switch_probability = switch_probability;
        in.persisted = persisted->score;
        RepairSettings rs;
        rs.chunk_haplotypes = options.chunk_haplotypes;
        for (const BlockAudit& ba : audit.blocks) {
            if (ba.route != CensusRoute::RepeatedUnit) { allocation_status[ba.block] = {"NOT_APPLICABLE", "not_repeated_unit"}; continue; }
            const UnitLibrary lib = build_unit_library(model.blocks()[ba.block], nodes, rs);
            std::vector<std::size_t> u1 = lib.panel.at(ba.call_alleles[0]), u2 = lib.panel.at(ba.call_alleles[1]);
            const RepairResult* r = repair_of.count(ba.block) ? &repairs[repair_of.at(ba.block)].second : nullptr;
            const bool repaired = r && r->outcome == "REPAIRED" && r->chosen >= 0;
            if (repaired) {
                u1 = r->candidates[static_cast<std::size_t>(r->chosen)].units_a;
                u2 = r->candidates[static_cast<std::size_t>(r->chosen)].units_b;
            }
            if (u1.size() + u2.size() < 3 || u1.empty() || u2.empty()) { allocation_status[ba.block] = {"NOT_APPLICABLE", "fewer_than_three_units"}; continue; }
            if (options.allocation_test) {
                try {
                    allocation[ba.block] = allocation_observability_test(lib, u1, u2, ba.block, in, rs);
                } catch (const std::exception& e) {
                    AllocationTestResult failed;
                    failed.status = "INCONCLUSIVE";
                    allocation[ba.block] = failed;
                }
                ++summary.allocation_blocks_tested;
            }
            if (repaired && r->allocation_tie > 1) allocation_status[ba.block] = {"UNRESOLVED", "tied_allocation_equivalence_set"};
            else if (allocation.count(ba.block)) allocation_status[ba.block] = {allocation.at(ba.block).status, "observability_test"};
            else allocation_status[ba.block] = {"NOT_TESTED", "test_not_run"};
        }
        summary.seconds_allocation = seconds_since(ta);
    }
    if (options.allocation_test) {
        std::ostringstream al;
        al << std::setprecision(17) << "block_id\tstatus\tcopies\tswitch_coordinate\treference_fnv64_1e9\tarrangement\tflank_order\tshape\t"
                                         "fragment_fnv64_1e9\tfragments_differing\texposure_delta\tlikelihood_delta\tidentical\n";
        for (const auto& [b, a] : allocation) {
            if (a.rows.empty()) {
                al << model.blocks()[b].block_id << '\t' << a.status << '\t' << a.copies << '\t' << a.switch_coordinate << "\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
            }
            for (const AllocationRow& r : a.rows) {
                al << model.blocks()[b].block_id << '\t' << a.status << '\t' << a.copies << '\t' << a.switch_coordinate << '\t' << std::hex
                   << a.reference_fnv64_1e9 << '\t' << std::dec << r.arrangement << '\t' << r.flank_order << '\t' << r.shape << '\t' << std::hex
                   << r.fragment_fnv64_1e9 << std::dec << '\t' << r.fragments_differing << '\t' << r.exposure_delta << '\t' << r.likelihood_delta
                   << '\t' << (r.identical ? "true" : "false") << '\n';
            }
        }
        atomic_write(out_prefix + ".dosage_allocation.tsv", al.str());
    }

    // the dosage block table: separate fields, never collapsed (spec section 11)
    const auto bg = read_block_genotypes(run_prefix + ".block_genotypes.tsv");
    std::ostringstream o;
    o << std::setprecision(17);
    o << "block_id\tcensus_route\tmolecule_span_bp\tdosage_status\tdosage_reason\tallocation_status\tread_support_status\t"
         "representation_status\talpha\tm\tz_critical\tfeatures_tested\tmax_robust_z\treference_kmers\treference_hits\tretained_pairs\t"
         "original_allele1\toriginal_allele2\trepaired_allele1\trepaired_allele2\trepair_route\tcandidate_source\trepair_scope\t"
         "repair_scope_reason\tcandidates_scored\ttop_set_size\torientation_tie\tallocation_tie\tcap_reached\t"
         "original_fragment_term\toriginal_exposure_term\toriginal_likelihood\toriginal_prior\toriginal_total\t"
         "repaired_fragment_term\trepaired_exposure_term\trepaired_likelihood\trepaired_prior\trepaired_prior_spellable\t"
         "repair_separation\tconditional_repair_gap\toriginal_GQ\toriginal_gq_status\trepaired_GQ\tgq_status\t"
         "allocation_basis\tallocation_equivalence_size\tallocation_representative\n";
    for (const BlockAudit& ba : audit.blocks) {
        const std::map<std::string, std::string> g = bg.count(ba.block_id) ? bg.at(ba.block_id) : std::map<std::string, std::string>{};
        const RepairResult* r = repair_of.count(ba.block) ? &repairs[repair_of.at(ba.block)].second : nullptr;
        const BlockCandidate* c = r && r->chosen >= 0 ? &r->candidates[static_cast<std::size_t>(r->chosen)] : nullptr;
        ++summary.status_counts[dosage_status_name(ba.status)];
        if (ba.reason != "none") ++summary.reason_counts[ba.reason];
        o << ba.block_id << '\t' << census_route_name(ba.route) << '\t' << *span << '\t' << dosage_status_name(ba.status) << '\t' << ba.reason
          << '\t' << allocation_status.at(ba.block).first << '\t' << field(g, "genotype_support") << '\t'
          << field(g, "filter") << '\t' << audit.alpha << '\t' << audit.m << '\t';
        if (audit.z_critical) o << *audit.z_critical; else o << "NA";
        o << '\t' << ba.features.size() << '\t' << ba.max_robust_z << '\t' << ba.reference_kmers << '\t' << ba.reference_hits << '\t'
          << ba.retained_pairs << '\t' << ba.call_alleles[0] << '\t' << ba.call_alleles[1] << '\t';
        if (c) o << c->ctx1_id << '\t' << c->ctx2_id << '\t';
        else o << "NA\tNA\t";
        if (r) {
            o << r->route << '\t' << (r->route == "pair" ? "panel_pair" : "unit_multiset") << '\t' << r->scope << '\t' << r->scope_reason
              << '\t' << r->candidates.size() << '\t' << r->top_set_size << '\t' << (r->orientation_tie ? "true" : "false") << '\t'
              << r->allocation_tie << '\t' << r->cap_reached << '\t' << r->original_fragment << '\t' << r->original_exposure << '\t'
              << r->original_likelihood << '\t';
        } else {
            o << "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t";
        }
        o << persisted->prior << '\t' << persisted->score.total << '\t';
        if (c) {
            const std::string h1 = call.haplotypes[0].substr(0, call.intervals[0][ba.block].begin) + c->ctx1_sequence +
                                   call.haplotypes[0].substr(call.intervals[0][ba.block].end);
            const std::string h2 = call.haplotypes[1].substr(0, call.intervals[1][ba.block].begin) + c->ctx2_sequence +
                                   call.haplotypes[1].substr(call.intervals[1][ba.block].end);
            const double prior = log_unordered_sequence_pair_prior(model, h1, h2, switch_probability);
            o << c->fragment_term << '\t' << c->exposure_term << '\t' << c->likelihood << '\t' << prior << '\t'
              << (std::isinf(prior) ? "false" : "true") << '\t' << r->repair_separation << '\t';
            if (r->gap_defined) o << r->conditional_repair_gap; else o << "NA";
        } else {
            o << "NA\tNA\tNA\tNA\tNA\tNA\tNA";
        }
        const bool repaired = ba.status == DosageStatus::Repaired;
        o << '\t' << field(g, "conditional_GQ") << '\t' << "conditional_uncalibrated" << '\t' << (repaired ? "NA" : field(g, "conditional_GQ"))
          << '\t' << (repaired ? "UNAVAILABLE_REPAIRED" : "UNCHANGED_CONDITIONAL") << '\t' << allocation_status.at(ba.block).second << '\t';
        if (r && r->route == "repeated_unit" && c) {
            o << r->allocation_tie << '\t' << (r->representative_only ? "REPRESENTATIVE_ONLY_NOT_RESOLVED" : "SINGLE_ARRANGEMENT_AT_PICK");
        } else {
            o << "NA\tNA";
        }
        o << '\n';
    }
    std::ostringstream cen;
    cen << "block_id\talleles\tmax_node_multiplicity\tmax_repeated_bp\tdosage_only_pairs\tcensus_route\tpooled_features\n";
    for (const CensusBlock& b : census) {
        cen << b.block_id << '\t' << b.alleles << '\t' << b.max_node_multiplicity << '\t' << b.max_repeated_bp << '\t' << b.dosage_only_pairs
            << '\t' << census_route_name(b.route) << '\t' << features[b.block].features.size() << '\n';
    }
    // provenance of every chosen repaired haplotype (spec section 9): only the panel supplied to the caller
    if (options.repair && options.provenance) {
        const auto tp = std::chrono::steady_clock::now();
        std::vector<PanelSequence> templates;
        for (std::size_t t = 0; t < model.panel().template_count(); ++t) {
            templates.push_back({model.panel().template_names[t], spell_panel_template(model, t).bases});
        }
        const std::string template_digest = panel_digest(templates);
        std::ostringstream searched;   // every panel sequence a nearest-sequence search may compare against
        searched << "level\tblock_id\tpanel_id\tsha256\tbp\n";
        for (const PanelSequence& t : templates) searched << "full_haplotype\tNA\t" << t.id << '\t' << sha256_hex(t.sequence) << '\t' << t.sequence.size() << '\n';
        std::ostringstream repaired_fa, original_fa, pv;
        pv << std::setprecision(17);
        pv << "block_id\trepair_scope\thomologue\trepaired_sha256\toriginal_sha256\tlevel\trank\tpanel_id\tpanel_sha256\traw_distance\t"
              "alignment_columns\tnormalized_distance\tcigar\tpanel_size\tpanel_digest\tdetail\n";
        std::map<std::string, std::string> display;   // display id -> full id, for collision detection
        for (const auto& entry : repairs) {
            const std::string& block_id = entry.first;
            const RepairResult& r = entry.second;
            if (r.outcome != "REPAIRED" || r.chosen < 0) continue;
            const BlockCandidate& c = r.candidates[static_cast<std::size_t>(r.chosen)];
            std::size_t b = 0;
            while (model.blocks()[b].block_id != block_id) ++b;
            const PreparedGenotypeBlock& block = model.blocks()[b];
            std::vector<PanelSequence> alleles;
            for (const CanonicalBlockAllele& a : block.alleles) alleles.push_back({std::to_string(a.id), a.sequence});
            const std::string allele_digest = panel_digest(alleles);
            for (const PanelSequence& a : alleles) searched << "block_allele\t" << block_id << '\t' << a.id << '\t' << sha256_hex(a.sequence) << '\t' << a.sequence.size() << '\n';
            std::optional<UnitLibrary> lib;
            if (r.route == "repeated_unit") lib = build_unit_library(block, nodes, RepairSettings{});
            for (int h = 0; h < 2; ++h) {
                const std::string& seq = h == 0 ? c.ctx1_sequence : c.ctx2_sequence;
                const SequenceInterval iv = call.intervals[h][b];
                const std::string full = call.haplotypes[h].substr(0, iv.begin) + seq + call.haplotypes[h].substr(iv.end);
                if (full == call.haplotypes[h]) continue;   // this homologue is unchanged: no new haplotype
                const std::string rid = sha256_hex(full), oid = sha256_hex(call.haplotypes[h]);
                const std::string shown = rid.substr(0, 12);
                if (!register_display_id(display, rid)) throw std::runtime_error("display_id_collision " + shown);
                const std::string hom = h == 0 ? "hap1" : "hap2";
                repaired_fa << '>' << rid << " display=" << shown << " block=" << block_id << " homologue=" << hom << " bp=" << full.size()
                            << " repair_scope=" << r.scope << " allocation_status=" << allocation_status.at(b).first;
                if (r.representative_only) {
                    repaired_fa << " allocation=REPRESENTATIVE_ONLY equivalence_set=" << out_prefix << ".dosage_allocation_equivalence.tsv#" << block_id
                                << " equivalence_size=" << r.allocation_tie;
                }
                repaired_fa << '\n';
                for (std::size_t at = 0; at < full.size(); at += 60) repaired_fa << full.substr(at, 60) << '\n';
                original_fa << '>' << oid << " block=" << block_id << " homologue=" << hom << " bp=" << call.haplotypes[h].size() << '\n';
                for (std::size_t at = 0; at < call.haplotypes[h].size(); at += 60) original_fa << call.haplotypes[h].substr(at, 60) << '\n';
                const auto emit = [&](const char* level, const std::vector<NearestPanel>& near, std::size_t size, const std::string& digest) {
                    for (std::size_t k = 0; k < near.size(); ++k) {
                        pv << block_id << '\t' << r.scope << '\t' << hom << '\t' << rid << '\t' << oid << '\t' << level << '\t' << k << '\t' << near[k].id << '\t'
                           << near[k].sha256 << '\t' << near[k].distance << '\t' << near[k].alignment_columns << '\t' << near[k].normalized << '\t'
                           << near[k].cigar << '\t' << size << '\t' << digest << "\tnearest\n";
                    }
                };
                emit("full_haplotype", nearest_panel_sequences(full, templates), templates.size(), template_digest);
                emit("block_allele", nearest_panel_sequences(seq, alleles), alleles.size(), allele_digest);
                if (lib) {
                    const std::vector<std::size_t>& units = h == 0 ? c.units_a : c.units_b;
                    for (std::size_t k = 0; k < units.size(); ++k) {
                        const Unit& u = lib->units.at(units[k]);
                        std::string carriers;
                        for (const auto& cr : u.carriers) {
                            carriers += (carriers.empty() ? "" : ";") + std::to_string(cr[0]) + ":" + std::to_string(cr[1]) + "/" + std::to_string(cr[2]);
                        }
                        pv << block_id << '\t' << r.scope << '\t' << hom << '\t' << rid << '\t' << oid << "\tunit\t" << k << '\t' << u.id << '\t' << sha256_hex(u.sequence)
                           << "\tNA\tNA\tNA\tNA\t" << lib->units.size() << "\tNA\tstate=" << u.state << ",carriers(allele:position/units)=" << carriers << '\n';
                    }
                } else {
                    const std::string& aid = h == 0 ? c.ctx1_id : c.ctx2_id;
                    pv << block_id << '\t' << r.scope << '\t' << hom << '\t' << rid << '\t' << oid << "\tsource_allele\t0\t" << aid << '\t'
                       << sha256_hex(block.alleles.at(std::stoul(aid)).sequence) << "\tNA\tNA\tNA\tNA\t" << alleles.size() << '\t'
                       << allele_digest << "\tpanel_pair\n";
                }
            }
        }
        atomic_write(out_prefix + ".dosage_repaired.fa", repaired_fa.str());
        atomic_write(out_prefix + ".dosage_original.fa", original_fa.str());
        atomic_write(out_prefix + ".dosage_provenance.tsv", pv.str());
        atomic_write(out_prefix + ".dosage_provenance_panel.tsv", searched.str());
        summary.seconds_provenance = seconds_since(tp);
    }
    atomic_write(out_prefix + ".dosage_census.tsv", cen.str());
    atomic_write(out_prefix + ".dosage_features.tsv", dosage_features_table(audit));
    {
        // per-allele copy vector of every tested feature, so any block genotype's admissibility can be checked
        std::ostringstream fc;
        fc << "block_id\tfeature\tpooled_index\tcopies_by_allele\tadmissible\n";
        for (const BlockAudit& ba : audit.blocks) {
            for (const FeatureAudit& fa : ba.features) {
                std::string copies, adm;
                for (std::uint16_t v : features[ba.block].features[fa.pooled_index].copies) copies += (copies.empty() ? "" : ",") + std::to_string(v);
                for (int a : fa.admissible) adm += (adm.empty() ? "" : ",") + std::to_string(a);
                fc << ba.block_id << '\t' << fa.feature << '\t' << fa.pooled_index << '\t' << copies << '\t' << adm << '\n';
            }
        }
        atomic_write(out_prefix + ".dosage_feature_copies.tsv", fc.str());
    }
    atomic_write(out_prefix + ".dosage_blocks.tsv", o.str());
    if (options.repair) {
        // every distinct unordered genotype of each repair's allocation-equivalence set; the representative only
        // serialises the set (smallest canonical SHA-256 pair) and is not a resolved allocation
        std::ostringstream eq;
        eq << std::setprecision(17) << "block_id\trepair_scope\tequivalence_size\tmember\tcanonical_sha256_a\tcanonical_sha256_b\t"
                                         "units_a\tunits_b\tstates\tpanel_log_prior\tlikelihood\tserialization_representative\n";
        for (const auto& entry : repairs) {
            const RepairResult& r = entry.second;
            if (r.route != "repeated_unit" || r.outcome != "REPAIRED") continue;
            std::vector<const BlockCandidate*> members;
            std::set<std::pair<std::string, std::string>> seen;
            for (const BlockCandidate& c : r.candidates) {
                if (!c.allocation_tied) continue;
                const auto key = c.ctx1_id <= c.ctx2_id ? std::make_pair(c.ctx1_id, c.ctx2_id) : std::make_pair(c.ctx2_id, c.ctx1_id);
                if (seen.insert(key).second) members.push_back(&c);
            }
            std::sort(members.begin(), members.end(), [](const BlockCandidate* a, const BlockCandidate* b) {
                const auto ka = std::minmax(a->ctx1_id, a->ctx2_id), kb = std::minmax(b->ctx1_id, b->ctx2_id);
                return ka < kb;
            });
            for (std::size_t i = 0; i < members.size(); ++i) {
                const BlockCandidate& c = *members[i];
                const auto key = std::minmax(c.ctx1_id, c.ctx2_id);
                std::string ua, ub;
                for (std::size_t u : c.units_a) ua += (ua.empty() ? "" : "-") + std::to_string(u);
                for (std::size_t u : c.units_b) ub += (ub.empty() ? "" : "-") + std::to_string(u);
                eq << entry.first << '\t' << r.scope << '\t' << members.size() << '\t' << i << '\t' << key.first << '\t' << key.second << '\t'
                   << ua << '\t' << ub << '\t' << c.states << '\t' << c.panel_log_prior << '\t' << c.likelihood << '\t'
                   << (c.chosen ? (r.representative_only ? "REPRESENTATIVE_ONLY_NOT_RESOLVED" : "SINGLE_ARRANGEMENT_AT_PICK") : "member") << '\n';
            }
        }
        atomic_write(out_prefix + ".dosage_allocation_equivalence.tsv", eq.str());
    }
    if (options.repair) {
        atomic_write(out_prefix + ".dosage_repairs.tsv", repair_blocks_table(repairs));
        atomic_write(out_prefix + ".dosage_repair_candidates.tsv", repair_candidates_table(repairs));
    }
    write_status();
    return summary;
}

} // namespace panvar::dosage
