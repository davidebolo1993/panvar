#include "panvar/genotype_sequence_report.hpp"

#include "panvar/fragment_likelihood.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

namespace panvar {

namespace {

using Genotype = std::pair<std::size_t, std::size_t>;
using Labels = std::vector<std::size_t>;

std::string join(const Labels& labels) {
    std::string out;
    for (std::size_t i = 0; i < labels.size(); ++i) out += (i ? "," : "") + std::to_string(labels[i]);
    return out;
}

std::string genotype_text(const Genotype& g) { return std::to_string(g.first) + "/" + std::to_string(g.second); }

std::string number(double value, int precision) {
    if (std::isnan(value)) return "NA";
    std::ostringstream out;
    out << std::setprecision(precision) << value;
    return out.str();
}

std::string digest_hex(const SequenceDigest& d) {
    char buffer[40];
    std::snprintf(buffer, sizeof buffer, "%016llx%016llx", static_cast<unsigned long long>(d.a),
                  static_cast<unsigned long long>(d.b));
    return buffer;
}

std::string basename_of(const std::string& prefix) {
    const std::size_t slash = prefix.find_last_of('/');
    return slash == std::string::npos ? prefix : prefix.substr(slash + 1);
}

// Orientation (S1-C): haplotype A is the one whose spelled sequence has the smaller SHA-256 -- the representative
// FASTA's hap1. Each group's two sequences are spelled once and their SHA-256 bound to the caller's digests, so a
// member's labellings are oriented by sequence identity, never by the internal digest order (which can disagree with
// SHA-256 order). For two identical sequences the lexicographically smaller labelling is A.
struct GroupOrientation {
    SequenceDigest a_digest;            // the digest of haplotype A's sequence
    bool identical = false;
    std::string a_sha256, b_sha256, a_bases, b_bases;
};

GroupOrientation orientation_of(const PoolCandidate& c, const HaplotypeSpeller& spell) {
    GroupOrientation o;
    std::string x = spell(c.first), y = spell(c.second);
    std::string hx = sha256_bytes(x), hy = sha256_bytes(y);
    o.identical = c.first_sequence == c.second_sequence;
    if (hy < hx) {
        o.a_digest = c.second_sequence;
        o.a_sha256 = std::move(hy);
        o.b_sha256 = std::move(hx);
        o.a_bases = std::move(y);
        o.b_bases = std::move(x);
    } else {
        o.a_digest = c.first_sequence;
        o.a_sha256 = std::move(hx);
        o.b_sha256 = std::move(hy);
        o.a_bases = std::move(x);
        o.b_bases = std::move(y);
    }
    return o;
}

std::pair<Labels, Labels> oriented(const PoolCandidate& c, const GroupOrientation& o) {
    Labels a = c.first, b = c.second;
    if (o.identical) {
        if (b < a) std::swap(a, b);
    } else if (!(c.first_sequence == o.a_digest)) {
        std::swap(a, b);
    }
    return {std::move(a), std::move(b)};
}

std::string fasta_record(const std::string& name, const std::string& header, const std::string& bases) {
    std::string out = ">" + name + " " + header + "\n";
    for (std::size_t i = 0; i < bases.size(); i += 60) out += bases.substr(i, 60) + "\n";
    if (bases.empty()) out += "\n";
    return out;
}

}  // namespace

SequenceReport build_sequence_report(const PooledDecision& decided, const std::vector<std::string>& block_ids,
                                     const HaplotypeSpeller& spell, std::size_t alternative_rows) {
    const std::vector<PoolCandidate>& candidates = decided.candidates;
    const std::size_t block_count = block_ids.size();
    if (decided.ranking.candidates.empty() || decided.top_groups.empty()) {
        throw std::invalid_argument("a sequence report needs a completed ranking with a decision set");
    }
    if (decided.block_genotype_quality.size() != block_count) {
        throw std::invalid_argument("the sequence report's block ids do not match the decision");
    }
    SequenceReport r;
    r.built = true;
    r.block_ids = block_ids;
    r.completion = "completed_ranking";
    r.completion_note = "a completed, valid exact ranking of the explored pool";
    r.prediction = "PREDICTION";
    r.legacy_phase_support = decided.phase_support;
    r.pool_groups = decided.groups.size();
    r.explored_labellings = candidates.size();
    r.decision_set_groups = decided.top_groups.size();

    // ranking entry per group
    std::vector<const CertifiedCandidate*> entry_of(decided.groups.size(), nullptr);
    std::vector<std::size_t> rank_of(decided.groups.size(), 0);
    for (std::size_t i = 0; i < decided.ranking.candidates.size(); ++i) {
        entry_of[decided.ranking.candidates[i].candidate_id] = &decided.ranking.candidates[i];
        rank_of[decided.ranking.candidates[i].candidate_id] = i + 1;
    }
    // explored genotype sets per group and block (every member labelling counted)
    std::vector<std::vector<std::set<Genotype>>> genos(decided.groups.size(), std::vector<std::set<Genotype>>(block_count));
    for (std::size_t g = 0; g < decided.groups.size(); ++g) {
        for (const std::size_t m : decided.groups[g]) {
            for (std::size_t b = 0; b < block_count; ++b) {
                genos[g][b].insert(std::minmax(candidates[m].first[b], candidates[m].second[b]));
            }
        }
    }

    // ---- the canonical representative: the highest-scoring decision-set group. Every decision-set group is tied with
    // the best under the caller's own certified numerical bounds (that is what puts it in the decision set), so the
    // tie is broken by the sorted SHA-256 pair of the spelled sequences -- never by discovery order, and never by
    // floating-point summation noise (sample order changes a prior in its last bits). Bitwise-equal scores are counted
    // for information only.
    const double best_score = decided.ranking.candidates.front().score;
    const std::vector<std::size_t> tied(decided.top_groups.begin(), decided.top_groups.end());
    r.exact_score_ties = 0;
    for (const std::size_t g : decided.top_groups) r.exact_score_ties += entry_of[g]->score == best_score ? 1 : 0;
    // every decision-set group is spelled once (its representative member); orientation is bound to SHA-256
    std::map<std::size_t, GroupOrientation> orient;
    for (const std::size_t g : decided.top_groups) {
        orient.emplace(g, orientation_of(candidates[decided.representative[g]], spell));
    }
    std::size_t rep = decided.ranking.candidates.front().candidate_id;
    if (tied.size() > 1) {
        r.representative_tie_break = "canonical_sha256";
        for (const std::size_t g : tied) {
            const GroupOrientation& o = orient.at(g), &best = orient.at(rep);
            if (std::tie(o.a_sha256, o.b_sha256) < std::tie(best.a_sha256, best.b_sha256)) rep = g;
        }
    }
    r.representative_group = rep;
    r.representative_member = decided.representative[rep];
    {
        GroupOrientation& o = orient.at(rep);
        const auto ab = oriented(candidates[r.representative_member], o);
        r.representative_first = ab.first;    // spells hap1 (A)
        r.representative_second = ab.second;  // spells hap2 (B)
        r.hap1_sha256 = o.a_sha256;
        r.hap2_sha256 = o.b_sha256;
        r.hap1_bases = o.a_bases;
        r.hap2_bases = o.b_bases;
    }
    r.representative_score = entry_of[rep]->score;
    r.representative_probability = entry_of[rep]->candidate_set_probability;
    r.decision_set_probability = 0.0;
    for (const std::size_t g : decided.top_groups) r.decision_set_probability += entry_of[g]->candidate_set_probability;

    // ---- representation intervals of every decision-set group: ONE conservative enclosing range per group, from a
    // linear scan -- the first and last block where any explored labelling differs from the first one (two labellings
    // differ at a block only if one of them differs from the first there). Unchanged blocks in between are inside the
    // range, so coupled alternatives keep their context; tuples are the whole observed labellings over the range.
    for (const std::size_t g : decided.top_groups) {
        const GroupOrientation& o = orient.at(g);
        std::vector<std::pair<Labels, Labels>> members;
        for (const std::size_t m : decided.groups[g]) members.push_back(oriented(candidates[m], o));
        std::size_t lo = block_count, hi = 0;
        for (std::size_t i = 1; i < members.size(); ++i) {
            for (std::size_t b = 0; b < block_count; ++b) {
                if (members[i].first[b] != members[0].first[b] || members[i].second[b] != members[0].second[b]) {
                    lo = std::min(lo, b);
                    hi = std::max(hi, b);
                }
            }
        }
        if (lo == block_count) continue;
        RepresentationInterval interval;
        interval.group = g;
        interval.first_block = lo;
        interval.last_block = hi;
        interval.hap_a_sha256 = o.a_sha256;
        interval.hap_b_sha256 = o.b_sha256;
        std::set<std::pair<Labels, Labels>> tuples;
        for (const auto& m : members) {
            Labels a(m.first.begin() + lo, m.first.begin() + hi + 1);
            Labels b(m.second.begin() + lo, m.second.begin() + hi + 1);
            if (o.identical && b < a) std::swap(a, b);
            tuples.insert({std::move(a), std::move(b)});
        }
        interval.tuples.assign(tuples.begin(), tuples.end());
        r.intervals.push_back(std::move(interval));
    }

    // ---- block rows
    std::set<std::size_t> top(decided.top_groups.begin(), decided.top_groups.end());
    bool all_identified = true;
    for (std::size_t b = 0; b < block_count; ++b) {
        SequenceBlockRow row;
        row.block_id = block_ids[b];
        row.legacy_representative_gq = decided.block_genotype_quality[b];
        row.legacy_numerically_separated = decided.block_numerically_separated[b];
        std::map<Genotype, std::vector<std::size_t>> observed;
        for (const std::size_t g : decided.top_groups) {
            for (const Genotype& x : genos[g][b]) observed[x].push_back(g);
        }
        std::string obs;
        for (const auto& [x, groups] : observed) {
            std::vector<std::size_t> sorted = groups;
            std::sort(sorted.begin(), sorted.end());
            std::string ids;
            for (const std::size_t g : sorted) ids += (ids.empty() ? "" : ",") + std::to_string(g);
            obs += (obs.empty() ? "" : ";") + genotype_text(x) + ":" + ids;
        }
        row.observed_genotypes = obs.empty() ? "NA" : obs;
        if (observed.size() == 1) {
            row.status = "IDENTIFIED";
            row.genotype = observed.begin()->first;
            double lower = 0.0, upper = 0.0;
            for (std::size_t g = 0; g < decided.groups.size(); ++g) {
                bool all_other = true, any_other = false;
                for (const Genotype& x : genos[g][b]) {
                    if (x == row.genotype) all_other = false;
                    else any_other = true;
                }
                const double p = entry_of[g]->candidate_set_probability;
                if (all_other) lower += p;
                if (any_other) upper += p;
            }
            row.error_mass_lower = lower;
            row.error_mass_upper = upper;
        } else {
            all_identified = false;
            row.status = decided.top_groups.size() == 1 ? "REPRESENTATION_AMBIGUOUS" : "GENOTYPE_UNRESOLVED";
        }
        r.blocks.push_back(std::move(row));
    }
    for (const RepresentationInterval& interval : r.intervals) {
        if (interval.group != rep) continue;
        const std::string id = block_ids[interval.first_block] + ".." + block_ids[interval.last_block];
        for (std::size_t b = interval.first_block; b <= interval.last_block; ++b) r.blocks[b].interval = id;
    }
    r.sequence_status = decided.top_groups.size() == 1 ? "RESOLVED"
                        : all_identified               ? "GENOTYPE_RESOLVED_SEQUENCE_SET"
                                                       : "SEQUENCE_TIED";

    // ---- alternatives: the decision set in rank order, then the next groups by rank
    std::vector<std::size_t> listed(decided.top_groups.begin(), decided.top_groups.end());
    std::size_t extra = 0;
    for (const CertifiedCandidate& c : decided.ranking.candidates) {
        if (extra >= alternative_rows) break;
        if (top.count(c.candidate_id) != 0) continue;
        listed.push_back(c.candidate_id);
        ++extra;
    }
    std::sort(listed.begin(), listed.end(), [&](std::size_t x, std::size_t y) { return rank_of[x] < rank_of[y]; });
    for (const std::size_t g : listed) {
        SequenceAlternative a;
        const CertifiedCandidate& e = *entry_of[g];
        a.rank = rank_of[g];
        a.group = g;
        a.in_decision_set = top.count(g) != 0;
        a.representative = g == rep;
        a.score = e.score;
        a.fragment_term = e.fragment_term;
        a.exposure_term = e.exposure_term;
        a.prior = e.prior;
        a.numeric_error_bound = e.numeric_error_bound;
        a.candidate_set_probability = e.candidate_set_probability;
        a.explored_labellings = decided.groups[g].size();
        a.projection_ambiguous = decided.projection_ambiguous[g];
        const PoolCandidate& c = candidates[decided.groups[g].front()];
        a.first = c.first_sequence < c.second_sequence ? c.first_sequence : c.second_sequence;
        a.second = c.first_sequence < c.second_sequence ? c.second_sequence : c.first_sequence;
        for (std::size_t b = 0; b < block_count; ++b) {
            if (genos[g][b] != genos[rep][b]) a.blocks_with_other_genotypes.push_back(b);
        }
        r.alternatives.push_back(std::move(a));
    }
    return r;
}

SequenceReport provisional_sequence_report(const std::vector<std::size_t>& first, const std::vector<std::size_t>& second,
                                           double total_score, const std::string& stop, const std::string& note,
                                           const std::vector<std::string>& block_ids, const HaplotypeSpeller& spell) {
    if (first.size() != block_ids.size() || second.size() != block_ids.size()) {
        throw std::invalid_argument("a provisional incumbent needs one allele per block on both haplotypes");
    }
    SequenceReport r;
    r.built = true;
    r.block_ids = block_ids;
    r.prediction = "PROVISIONAL";
    r.sequence_status = "NOT_DETERMINED";
    r.completion = stop;
    r.completion_note = note;
    std::string x = spell(first), y = spell(second);
    std::string hx = sha256_bytes(x), hy = sha256_bytes(y);
    r.representative_first = first;
    r.representative_second = second;
    if (hy < hx) {
        std::swap(hx, hy);
        std::swap(x, y);
        std::swap(r.representative_first, r.representative_second);
    }
    r.hap1_sha256 = hx;
    r.hap2_sha256 = hy;
    r.hap1_bases = std::move(x);
    r.hap2_bases = std::move(y);
    r.representative_score = total_score;
    r.representative_tie_break = "not_applicable";
    for (std::size_t b = 0; b < block_ids.size(); ++b) {
        SequenceBlockRow row;
        row.block_id = block_ids[b];
        row.status = "PROVISIONAL_LABEL";
        row.genotype = std::minmax(first[b], second[b]);
        r.blocks.push_back(std::move(row));
    }
    return r;
}

const std::vector<std::string>& sequence_report_outputs() {
    // exact owned suffixes: the five reports and atomic_write's staging name for each
    static const std::vector<std::string> suffixes = {
        ".sequence_locus.tsv", ".sequence_blocks.tsv", ".sequence_haplotypes.fa", ".representation_intervals.tsv",
        ".sequence_alternatives.tsv", ".sequence_locus.tsv.writing", ".sequence_blocks.tsv.writing",
        ".sequence_haplotypes.fa.writing", ".representation_intervals.tsv.writing", ".sequence_alternatives.tsv.writing",
        ".sequence_haplotypes.tsv", ".sequence_haplotypes.tsv.writing"};   // the last two: S1-D (sequence_first_v2)
    return suffixes;
}

void annotate_panel_membership(SequenceReport& r, const std::vector<std::string>& names,
                               const std::vector<std::string>& sequences, const std::vector<bool>& reverse) {
    if (names.size() != sequences.size() || names.size() != reverse.size()) {
        throw std::invalid_argument("panel membership needs one sequence and one orientation per template");
    }
    std::map<std::string, std::vector<std::size_t>> by_hash;
    for (std::size_t t = 0; t < sequences.size(); ++t) by_hash[sha256_bytes(sequences[t])].push_back(t);
    const std::string* bases[2] = {&r.hap1_bases, &r.hap2_bases};
    const std::string* hashes[2] = {&r.hap1_sha256, &r.hap2_sha256};
    for (int h = 0; h < 2; ++h) {
        HaplotypeOrigin o;
        o.evaluated = true;
        const auto found = by_hash.find(*hashes[h]);
        std::vector<std::pair<std::string, std::string>> matches;
        if (found != by_hash.end()) {
            for (const std::size_t t : found->second) {
                if (sequences[t] == *bases[h]) {   // exact verification, never the hash alone
                    matches.push_back({names[t], reverse[t] ? "reverse_complement" : "forward"});
                }
            }
        }
        std::sort(matches.begin(), matches.end());
        for (const auto& m : matches) {
            o.paths.push_back(m.first);
            o.orientations.push_back(m.second);
        }
        o.panel = !o.paths.empty();
        r.origin[h] = std::move(o);
    }
}

std::pair<bool, std::string> sequence_fasta_decision(const SequenceReport& r) {
    const bool predicted = r.built && r.prediction != "NO_PREDICTION" && !r.publication_failed;
    if (!predicted) return {false, "NO_PREDICTION"};
    if (r.schema != "sequence_first_v2") return {true, "SEQUENCE_FIRST_V1_ALWAYS"};
    const bool off = !r.origin[0].panel || !r.origin[1].panel;
    return off ? std::make_pair(true, std::string("OFF_PANEL_PAIR"))
               : std::make_pair(false, std::string("PANEL_REPRESENTED_PAIR"));
}

int sequence_report_exit_code(const SequenceReport& report) {
    if (report.publication_failed) return 2;
    if (report.prediction == "PREDICTION") return 0;
    if (report.prediction == "PROVISIONAL") return 3;
    return 1;
}

std::string sequence_report_status_keys(const SequenceReport& report) {
    std::ostringstream out;
    out << "sequence_prediction\t" << report.prediction << "\nsequence_status\t" << report.sequence_status
        << "\nsequence_completion\t" << report.completion << "\nsequence_exit_code\t" << sequence_report_exit_code(report)
        << "\nsequence_publication\t" << (report.publication_failed ? "FAILED" : "complete") << '\n';
    if (report.schema == "sequence_first_v2") {   // S1-D: intentional FASTA absence is never a publication failure
        const auto [emitted, reason] = sequence_fasta_decision(report);
        out << "sequence_fasta\t" << (emitted ? "emitted" : "not_emitted") << "\nsequence_fasta_reason\t" << reason << '\n';
    }
    return out.str();
}

void write_sequence_report(const SequenceReport& r, const SequenceLegacy& legacy, const std::string& prefix,
                           const std::string& no_prediction_reason) {
    const bool predicted = r.built && r.prediction != "NO_PREDICTION";
    std::size_t identified = 0, ambiguous = 0, unresolved = 0, provisional = 0;
    for (const SequenceBlockRow& row : r.blocks) {
        identified += row.status == "IDENTIFIED";
        ambiguous += row.status == "REPRESENTATION_AMBIGUOUS";
        unresolved += row.status == "GENOTYPE_UNRESOLVED";
        provisional += row.status == "PROVISIONAL_LABEL";
    }
    std::string intervals;
    for (const RepresentationInterval& i : r.intervals) {
        if (i.group != r.representative_group) continue;
        intervals += (intervals.empty() ? "" : ",") + r.block_ids[i.first_block] + ".." + r.block_ids[i.last_block];
    }
    std::ostringstream locus;
    locus << "key\tvalue\nschema_version\t" << r.schema << "\nprediction\t" << (predicted ? r.prediction : "NO_PREDICTION")
          << "\nsequence_status\t" << (predicted ? r.sequence_status : "NA") << "\nprovisional\t"
          << (r.prediction == "PROVISIONAL" ? "yes" : "no") << "\ncompletion\t"
          << (predicted ? r.completion : legacy.run_status) << "\ncompletion_note\t"
          << (predicted ? r.completion_note : no_prediction_reason) << "\ndecision_scope\t"
          << (r.prediction == "PREDICTION" ? "numerical_top_set_of_the_explored_exactly_scored_pool" : "NOT_DETERMINED")
          << "\nresolution_meaning\tRESOLVED_means_one_sequence_pair_within_the_declared_decision_set;"
             "not_a_confidence_or_correctness_claim\n"
             "representative_scope\tcanonical_member_of_the_numerically_inseparable_decision_set;"
             "not_necessarily_the_strictly_highest_floating_point_score\n";
    if (predicted) {
        locus << "representative_group\t" << (r.prediction == "PREDICTION" ? std::to_string(r.representative_group) : "NA")
              << "\nrepresentative_of_set\t" << (r.decision_set_groups > 1 ? "yes" : "no")
              << "\nrepresentative_tie_break\t" << r.representative_tie_break << "\nexact_score_ties\t"
              << (r.prediction == "PREDICTION" ? std::to_string(r.exact_score_ties) : "NA")
              << "\nrepresentative_hap1_sha256\t" << r.hap1_sha256 << "\nrepresentative_hap1_bp\t" << r.hap1_bases.size()
              << "\nrepresentative_hap2_sha256\t" << r.hap2_sha256 << "\nrepresentative_hap2_bp\t" << r.hap2_bases.size()
              << "\nrepresentative_labelling_hap1\t" << join(r.representative_first)
              << "\nrepresentative_labelling_hap2\t" << join(r.representative_second)
              << "\nrepresentative_labelling_scope\tone_explored_labelling;alternatives_in_representation_intervals"
              << "\nrepresentative_score\t" << number(r.representative_score, 12) << "\ndecision_set_groups\t"
              << r.decision_set_groups << "\npool_groups\t" << r.pool_groups << "\nexplored_labellings\t"
              << r.explored_labellings << "\nidentified_blocks\t" << identified << "\nrepresentation_ambiguous_blocks\t"
              << ambiguous << "\ngenotype_unresolved_blocks\t" << unresolved << "\nprovisional_label_blocks\t" << provisional
              << "\nrepresentation_intervals\t" << (intervals.empty() ? "none" : intervals)
              << "\nrepresentation_interval_scope\tconservative_enclosing_not_minimal"
              << "\nrepresentative_probability\t" << number(r.representative_probability, 6)
              << "\ndecision_set_probability\t" << number(r.decision_set_probability, 6) << '\n';
    }
    locus << "probability_scope\tcandidate_conditional_over_the_explored_pool;not_calibrated\nstructure_confidence\tNA\n"
             "structure_confidence_reason\tUNCALIBRATED\nprojection_coverage\tEXPLORED\n"
             "projection_coverage_scope\tuniqueness_and_error_mass_bounds_are_conditional_on_explored_labellings\n"
             "legacy_representative_gq_scope\tlegacy_diagnostic_one_labelling_per_group_not_projection_aware\n"
             "legacy_phase_support\t"
          << r.legacy_phase_support
          << "\nphase_assertion\tnone_beyond_legacy_phase_support\n"
             "dosage_scope\tcaller_decision_before_the_dosage_stage;the_dosage_stage_changes_legacy_outputs_only\n"
             "post_dosage_sequence\tnot_reported:when_dosage_repair_changes_the_legacy_reported_genotype_this_report_"
             "and_its_FASTA_are_NOT_the_post_dosage_result\n"
             "dosage_mode\t"
          << legacy.dosage_mode << "\n"
             "legacy_run_status\t"
          << legacy.run_status << "\nlegacy_call_class\t" << legacy.call_class << "\nlegacy_no_call_reason\t"
          << legacy.no_call_reason << "\nexit_code\t" << (predicted ? sequence_report_exit_code(r) : 1) << '\n';
    const auto fasta_decision = sequence_fasta_decision(r);
    if (r.schema == "sequence_first_v2") {
        // S1-D: panel references, the FASTA decision and the diagnostic dosage contract
        const auto origin = [&](int h) {
            return !predicted ? std::string("NA") : r.origin[h].panel ? std::string("PANEL") : std::string("OFF_PANEL");
        };
        const auto paths = [&](int h) {
            std::string out;
            for (const std::string& p : r.origin[h].paths) out += (out.empty() ? "" : ",") + p;
            return out.empty() ? std::string("NA") : out;
        };
        locus << "fasta_emitted\t" << (fasta_decision.first ? "yes" : "no") << "\nfasta_reason\t" << fasta_decision.second
              << "\nfasta_policy\tsequence_first_v2:FASTA_only_when_a_predicted_haplotype_is_not_spelled_by_an_active_panel_"
                 "template;panel_haplotypes_are_reported_by_extraction_reference\nhap1_origin\t"
              << origin(0) << "\nhap1_panel_paths\t" << paths(0) << "\nhap2_origin\t" << origin(1)
              << "\nhap2_panel_paths\t" << paths(1)
              << "\npanel_membership_rule\texact_sequence_identity_in_the_locus_frame_against_the_active_curated_panel;"
                 "SHA-256_lookup_with_byte_verification;no_truth_labels_or_similarity\ngraph_sha256\t"
              << r.graph_sha256 << "\nmodel_fingerprint\t" << r.model_fingerprint
              << "\nprepared_index_identity\tsee_status_prepared_locus_file_sha256_when_loaded\n"
                 "dosage_contract\tdiagnostic_audit_only;never_modifies_this_prediction;see_dosage_status.tsv\n"
                 "dosage_mode_source\t"
              << legacy.dosage_mode_source << '\n';
    }
    const std::string locus_text = locus.str();   // the commit record: published last (S1-C)

    std::ostringstream blocks;
    blocks << "block_id\tblock_status\tgenotype\tobserved_genotypes\terror_mass_lower\terror_mass_upper\t"
              "projection_coverage\trepresentation_interval\tlegacy_representative_GQ\tlegacy_numerically_separated\n";
    if (predicted) {
        for (const SequenceBlockRow& row : r.blocks) {
            const bool has_genotype = row.status == "IDENTIFIED" || row.status == "PROVISIONAL_LABEL";
            blocks << row.block_id << '\t' << row.status << '\t' << (has_genotype ? genotype_text(row.genotype) : "NA")
                   << '\t' << row.observed_genotypes << '\t' << number(row.error_mass_lower, 6) << '\t'
                   << number(row.error_mass_upper, 6) << '\t'
                   << (row.status == "PROVISIONAL_LABEL" ? "NOT_DETERMINED" : "EXPLORED") << '\t' << row.interval << '\t'
                   << number(row.legacy_representative_gq, 6) << '\t'
                   << (row.status == "PROVISIONAL_LABEL" ? "NA" : row.legacy_numerically_separated ? "true" : "false")
                   << '\n';
        }
    }
    atomic_write(prefix + ".sequence_blocks.tsv", blocks.str());

    std::ostringstream iv;
    // hapA/hapB: A is the group's sequence with the smaller SHA-256 (for the representative group: FASTA hap1, labelling
    // representative_labelling_hap1), B the other; their SHA-256 are written on every row.
    iv << "interval\tgroup\trepresentative_group\tfirst_block\tlast_block\tblocks\tinterval_scope\tobserved_tuples\t"
          "tuples_hapA_hapB\thapA_sha256\thapB_sha256\n";
    for (const RepresentationInterval& i : r.intervals) {
        std::string ids, tuples;
        for (std::size_t b = i.first_block; b <= i.last_block; ++b) ids += (ids.empty() ? "" : ",") + r.block_ids[b];
        for (const auto& t : i.tuples) tuples += (tuples.empty() ? "" : ";") + join(t.first) + "/" + join(t.second);
        iv << r.block_ids[i.first_block] << ".." << r.block_ids[i.last_block] << '\t' << i.group << '\t'
           << (i.group == r.representative_group ? "yes" : "no") << '\t' << r.block_ids[i.first_block] << '\t'
           << r.block_ids[i.last_block] << '\t' << ids << "\tconservative_enclosing_not_minimal\t" << i.tuples.size()
           << '\t' << tuples << '\t' << i.hap_a_sha256 << '\t' << i.hap_b_sha256 << '\n';
    }
    atomic_write(prefix + ".representation_intervals.tsv", iv.str());

    std::ostringstream alt;
    alt << "rank\tgroup\tin_decision_set\trepresentative\tscore\tfragment_term\texposure_term\tsequence_pair_prior\t"
           "numeric_error_bound\tcandidate_set_probability\texplored_labellings\tprojection_ambiguous_explored\t"
           "hapA_digest\thapA_bp\thapB_digest\thapB_bp\tblocks_with_other_genotypes\tother_genotype_blocks\n";
    for (const SequenceAlternative& a : r.alternatives) {
        std::string ids;
        for (const std::size_t b : a.blocks_with_other_genotypes) ids += (ids.empty() ? "" : ",") + r.block_ids[b];
        alt << a.rank << '\t' << a.group << '\t' << (a.in_decision_set ? "yes" : "no") << '\t'
            << (a.representative ? "yes" : "no") << '\t' << number(a.score, 12) << '\t' << number(a.fragment_term, 12)
            << '\t' << number(a.exposure_term, 12) << '\t' << number(a.prior, 12) << '\t'
            << number(a.numeric_error_bound, 6) << '\t' << number(a.candidate_set_probability, 6) << '\t'
            << a.explored_labellings << '\t' << (a.projection_ambiguous ? "true" : "false") << '\t'
            << digest_hex(a.first) << '\t' << a.first.length << '\t' << digest_hex(a.second) << '\t' << a.second.length
            << '\t' << a.blocks_with_other_genotypes.size() << '\t' << (ids.empty() ? "none" : ids) << '\n';
    }
    atomic_write(prefix + ".sequence_alternatives.tsv", alt.str());

    if (r.schema == "sequence_first_v2") {
        // S1-D: one row per predicted haplotype -- identity, origin, extraction reference, block-path provenance
        std::ostringstream hap;
        hap << "haplotype\tsha256\tlength_bp\torigin\tpanel_paths\tpath_orientations\textraction\tgraph_sha256\t"
               "model_fingerprint\tblock_path\trepresentation_intervals\tprovisional\n";
        if (predicted) {
            const std::string* sha[2] = {&r.hap1_sha256, &r.hap2_sha256};
            const std::string* bases[2] = {&r.hap1_bases, &r.hap2_bases};
            const std::vector<std::size_t>* labels[2] = {&r.representative_first, &r.representative_second};
            for (int h = 0; h < 2; ++h) {
                const HaplotypeOrigin& o = r.origin[h];
                std::string ps, os, ex;
                for (std::size_t k = 0; k < o.paths.size(); ++k) {
                    ps += (k ? "," : "") + o.paths[k];
                    os += (k ? "," : "") + o.orientations[k];
                }
                if (o.panel) {
                    ex = "gfa_path=" + o.paths.front() + ";interval=0-" + std::to_string(bases[h]->size()) +
                         ";orientation=" + o.orientations.front() + ";whole_locus_path";
                } else {
                    ex = "not_extractable:off_panel_spelled_candidate_sequence(not_an_assembly);see_FASTA";
                }
                hap << "hap" << (h + 1) << '\t' << *sha[h] << '\t' << bases[h]->size() << '\t'
                    << (o.panel ? "PANEL" : "OFF_PANEL") << '\t' << (ps.empty() ? "NA" : ps) << '\t'
                    << (os.empty() ? "NA" : os) << '\t' << ex << '\t' << r.graph_sha256 << '\t' << r.model_fingerprint
                    << '\t' << join(*labels[h]) << '\t' << (intervals.empty() ? "none" : intervals) << '\t'
                    << (r.prediction == "PROVISIONAL" ? "yes" : "no") << '\n';
            }
        }
        atomic_write(prefix + ".sequence_haplotypes.tsv", hap.str());
    }
    std::string fasta;
    if (fasta_decision.first && r.schema == "sequence_first_v2") {
        const std::string base = basename_of(prefix);
        std::ostringstream header;
        header << "prediction=" << r.prediction << " sequence_status=" << r.sequence_status
               << " representative_of_set=" << (r.decision_set_groups > 1 ? "yes" : "no")
               << " provisional=" << (r.prediction == "PROVISIONAL" ? "yes" : "no")
               << " representation_intervals=" << (intervals.empty() ? "none" : intervals)
               << " phase=not_asserted structure_confidence=NA scope=caller_decision_before_dosage_stage"
               << " fasta_reason=OFF_PANEL_PAIR spelled_candidate_not_assembly";
        const std::string* sha[2] = {&r.hap1_sha256, &r.hap2_sha256};
        const std::string* bases[2] = {&r.hap1_bases, &r.hap2_bases};
        for (int h = 0; h < 2; ++h) {
            const std::string origin = r.origin[h].panel ? "origin=PANEL panel_path=" + r.origin[h].paths.front()
                                                         : std::string("origin=OFF_PANEL");
            fasta += fasta_record(base + "_seq_hap" + std::to_string(h + 1), "sha256=" + *sha[h] + " bp=" +
                                  std::to_string(bases[h]->size()) + " " + origin + " " + header.str(), *bases[h]);
        }
        atomic_write(prefix + ".sequence_haplotypes.fa", fasta);
    } else if (r.schema == "sequence_first_v2") {
        std::remove((prefix + ".sequence_haplotypes.fa").c_str());   // intentional absence (fasta_reason), not a failure
    } else if (predicted) {
        const std::string base = basename_of(prefix);
        std::ostringstream header;
        header << "prediction=" << r.prediction << " sequence_status=" << r.sequence_status
               << " representative_of_set=" << (r.decision_set_groups > 1 ? "yes" : "no")
               << " provisional=" << (r.prediction == "PROVISIONAL" ? "yes" : "no")
               << " representation_intervals=" << (intervals.empty() ? "none" : intervals)
               << " phase=not_asserted structure_confidence=NA scope=caller_decision_before_dosage_stage"
               << " dosage_mode=" << legacy.dosage_mode;
        fasta = fasta_record(base + "_seq_hap1", "sha256=" + r.hap1_sha256 + " bp=" +
                                                     std::to_string(r.hap1_bases.size()) + " " + header.str(),
                             r.hap1_bases) +
                fasta_record(base + "_seq_hap2", "sha256=" + r.hap2_sha256 + " bp=" +
                                                     std::to_string(r.hap2_bases.size()) + " " + header.str(),
                             r.hap2_bases);
        atomic_write(prefix + ".sequence_haplotypes.fa", fasta);
    } else {
        std::remove((prefix + ".sequence_haplotypes.fa").c_str());
    }
    atomic_write(prefix + ".sequence_locus.tsv", locus_text);   // only once every required report is complete
}

bool publish_sequence_report(SequenceReport& r, const SequenceLegacy& legacy, const std::string& prefix,
                             const std::string& no_prediction_reason) {
    try {
        write_sequence_report(r, legacy, prefix, no_prediction_reason);
        return true;
    } catch (const std::exception& error) {
        // Never leave a partial report: remove every owned file of this run, then record the failure (best effort).
        for (const std::string& suffix : sequence_report_outputs()) std::remove((prefix + suffix).c_str());
        SequenceReport failed;
        failed.schema = r.schema;
        failed.publication_failed = true;
        failed.publication_error = error.what();
        failed.completion = "report_publication_failed";
        failed.legacy_phase_support = r.legacy_phase_support;
        r = failed;
        try {
            atomic_write(prefix + ".sequence_locus.tsv",
                         "key\tvalue\nschema_version\t" + failed.schema + "\nprediction\tNO_PREDICTION\nsequence_status\tNA\n"
                         "completion\treport_publication_failed\ncompletion_note\t" + std::string(error.what()) +
                         "\nlegacy_run_status\t" + legacy.run_status + "\nlegacy_call_class\t" + legacy.call_class +
                         "\nlegacy_no_call_reason\t" + legacy.no_call_reason + "\nexit_code\t2\n");
        } catch (const std::exception&) {
            std::remove((prefix + ".sequence_locus.tsv").c_str());
            std::remove((prefix + ".sequence_locus.tsv.writing").c_str());
        }
        return false;
    }
}

}  // namespace panvar
