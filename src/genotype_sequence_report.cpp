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

// A member's two labellings oriented so that A spells the group's first sequence (canonical digest order); for a
// pair of identical sequences, the lexicographically smaller labelling comes first.
std::pair<Labels, Labels> oriented(const PoolCandidate& c) {
    Labels a = c.first, b = c.second;
    if (c.second_sequence < c.first_sequence) {
        std::swap(a, b);
    } else if (c.first_sequence == c.second_sequence && b < a) {
        std::swap(a, b);
    }
    return {std::move(a), std::move(b)};
}

void write_or_throw(const std::string& path, const std::string& text) {
    const std::string temporary = path + ".tmp";
    {
        std::ofstream out(temporary, std::ios::binary);
        if (!out) throw std::runtime_error("cannot write " + temporary);
        out << text;
        if (!out) throw std::runtime_error("cannot write " + temporary);
    }
    if (std::rename(temporary.c_str(), path.c_str()) != 0) throw std::runtime_error("cannot publish " + path);
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
    std::size_t rep = decided.ranking.candidates.front().candidate_id;
    std::pair<std::string, std::string> rep_key;
    std::pair<std::string, std::string> rep_bases;
    const auto spelled = [&](std::size_t g) {
        const auto ab = oriented(candidates[decided.representative[g]]);
        std::string x = spell(ab.first), y = spell(ab.second);
        std::string hx = sha256_bytes(x), hy = sha256_bytes(y);
        if (hy < hx) {
            std::swap(hx, hy);
            std::swap(x, y);
        }
        return std::make_pair(std::make_pair(hx, hy), std::make_pair(std::move(x), std::move(y)));
    };
    if (tied.size() > 1) {
        r.representative_tie_break = "canonical_sha256";
        bool first = true;
        for (const std::size_t g : tied) {
            auto item = spelled(g);
            if (first || item.first < rep_key) {
                rep = g;
                rep_key = item.first;
                rep_bases = std::move(item.second);
                first = false;
            }
        }
    } else {
        auto item = spelled(rep);
        rep_key = item.first;
        rep_bases = std::move(item.second);
    }
    r.representative_group = rep;
    r.representative_member = decided.representative[rep];
    {
        const auto ab = oriented(candidates[r.representative_member]);
        r.representative_first = ab.first;
        r.representative_second = ab.second;
    }
    r.hap1_sha256 = rep_key.first;
    r.hap2_sha256 = rep_key.second;
    r.hap1_bases = std::move(rep_bases.first);
    r.hap2_bases = std::move(rep_bases.second);
    r.representative_score = entry_of[rep]->score;
    r.representative_probability = entry_of[rep]->candidate_set_probability;
    r.decision_set_probability = 0.0;
    for (const std::size_t g : decided.top_groups) r.decision_set_probability += entry_of[g]->candidate_set_probability;

    // ---- representation intervals of every decision-set group (conservative enclosing ranges)
    for (const std::size_t g : decided.top_groups) {
        std::vector<std::pair<Labels, Labels>> members;
        for (const std::size_t m : decided.groups[g]) members.push_back(oriented(candidates[m]));
        std::vector<std::pair<std::size_t, std::size_t>> spans;
        for (std::size_t i = 0; i < members.size(); ++i) {
            for (std::size_t j = i + 1; j < members.size(); ++j) {
                std::size_t lo = block_count, hi = 0;
                for (std::size_t b = 0; b < block_count; ++b) {
                    if (members[i].first[b] != members[j].first[b] || members[i].second[b] != members[j].second[b]) {
                        lo = std::min(lo, b);
                        hi = std::max(hi, b);
                    }
                }
                if (lo < block_count) spans.push_back({lo, hi});
            }
        }
        std::sort(spans.begin(), spans.end());
        std::vector<std::pair<std::size_t, std::size_t>> merged;
        for (const auto& s : spans) {
            if (!merged.empty() && s.first <= merged.back().second) {
                merged.back().second = std::max(merged.back().second, s.second);
            } else {
                merged.push_back(s);
            }
        }
        const bool same_sequences = candidates[decided.groups[g].front()].first_sequence ==
                                    candidates[decided.groups[g].front()].second_sequence;
        for (const auto& span : merged) {
            RepresentationInterval interval;
            interval.group = g;
            interval.first_block = span.first;
            interval.last_block = span.second;
            std::set<std::pair<Labels, Labels>> tuples;
            for (const auto& m : members) {
                Labels a(m.first.begin() + span.first, m.first.begin() + span.second + 1);
                Labels b(m.second.begin() + span.first, m.second.begin() + span.second + 1);
                if (same_sequences && b < a) std::swap(a, b);
                tuples.insert({std::move(a), std::move(b)});
            }
            interval.tuples.assign(tuples.begin(), tuples.end());
            r.intervals.push_back(std::move(interval));
        }
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

int sequence_report_exit_code(const SequenceReport& report) {
    if (report.prediction == "PREDICTION") return 0;
    if (report.prediction == "PROVISIONAL") return 3;
    return 1;
}

std::string sequence_report_status_keys(const SequenceReport& report) {
    std::ostringstream out;
    out << "sequence_prediction\t" << report.prediction << "\nsequence_status\t" << report.sequence_status
        << "\nsequence_completion\t" << report.completion << "\nsequence_exit_code\t" << sequence_report_exit_code(report)
        << '\n';
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
    locus << "key\tvalue\nschema_version\tsequence_first_v1\nprediction\t" << (predicted ? r.prediction : "NO_PREDICTION")
          << "\nsequence_status\t" << (predicted ? r.sequence_status : "NA") << "\nprovisional\t"
          << (r.prediction == "PROVISIONAL" ? "yes" : "no") << "\ncompletion\t"
          << (predicted ? r.completion : legacy.run_status) << "\ncompletion_note\t"
          << (predicted ? r.completion_note : no_prediction_reason) << "\ndecision_scope\t"
          << (r.prediction == "PREDICTION" ? "numerical_top_set_of_the_explored_exactly_scored_pool" : "NOT_DETERMINED")
          << "\nresolution_meaning\tRESOLVED_means_one_sequence_pair_within_the_declared_decision_set;"
             "not_a_confidence_or_correctness_claim\n";
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
             "dosage_scope\tcaller_decision_only;the_dosage_stage_changes_legacy_outputs_only\n"
             "legacy_run_status\t"
          << legacy.run_status << "\nlegacy_call_class\t" << legacy.call_class << "\nlegacy_no_call_reason\t"
          << legacy.no_call_reason << "\nexit_code\t" << (predicted ? sequence_report_exit_code(r) : 1) << '\n';
    write_or_throw(prefix + ".sequence_locus.tsv", locus.str());

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
    write_or_throw(prefix + ".sequence_blocks.tsv", blocks.str());

    std::ostringstream iv;
    iv << "interval\tgroup\trepresentative_group\tfirst_block\tlast_block\tblocks\tinterval_scope\tobserved_tuples\t"
          "tuples_hapA_hapB\n";
    for (const RepresentationInterval& i : r.intervals) {
        std::string ids, tuples;
        for (std::size_t b = i.first_block; b <= i.last_block; ++b) ids += (ids.empty() ? "" : ",") + r.block_ids[b];
        for (const auto& t : i.tuples) tuples += (tuples.empty() ? "" : ";") + join(t.first) + "/" + join(t.second);
        iv << r.block_ids[i.first_block] << ".." << r.block_ids[i.last_block] << '\t' << i.group << '\t'
           << (i.group == r.representative_group ? "yes" : "no") << '\t' << r.block_ids[i.first_block] << '\t'
           << r.block_ids[i.last_block] << '\t' << ids << "\tconservative_enclosing_not_minimal\t" << i.tuples.size()
           << '\t' << tuples << '\n';
    }
    write_or_throw(prefix + ".representation_intervals.tsv", iv.str());

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
    write_or_throw(prefix + ".sequence_alternatives.tsv", alt.str());

    std::string fasta;
    if (predicted) {
        const std::string base = basename_of(prefix);
        std::ostringstream header;
        header << "prediction=" << r.prediction << " sequence_status=" << r.sequence_status
               << " representative_of_set=" << (r.decision_set_groups > 1 ? "yes" : "no")
               << " provisional=" << (r.prediction == "PROVISIONAL" ? "yes" : "no")
               << " representation_intervals=" << (intervals.empty() ? "none" : intervals)
               << " phase=not_asserted structure_confidence=NA";
        fasta = fasta_record(base + "_seq_hap1", "sha256=" + r.hap1_sha256 + " bp=" +
                                                     std::to_string(r.hap1_bases.size()) + " " + header.str(),
                             r.hap1_bases) +
                fasta_record(base + "_seq_hap2", "sha256=" + r.hap2_sha256 + " bp=" +
                                                     std::to_string(r.hap2_bases.size()) + " " + header.str(),
                             r.hap2_bases);
        write_or_throw(prefix + ".sequence_haplotypes.fa", fasta);
    } else {
        std::remove((prefix + ".sequence_haplotypes.fa").c_str());
    }
}

}  // namespace panvar
