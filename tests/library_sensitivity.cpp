// Gates for the library-intensity sensitivity audit (library_sensitivity_audit_v1).
//
//   1  audit scores are exact: an audit library's score, computed beside the run's own in one combination pass, is
//      bit for bit the score a run declaring that library computes -- one-shot and prepared-context scoring, 1 and 4
//      threads -- and the run's own scores are unchanged by the audit.
//   2  the verdict compares unordered genotypes only:
//        a  the endpoints anchored differently but their top sets carry the reported genotype: STABLE;
//        b  an endpoint's top set carries another genotype at one block: SENSITIVE at exactly that block;
//        c  an endpoint's top set is a phase switch of the call (same unordered genotype): STABLE, and recorded as
//           phase_or_representative_only;
//        d  a tie of the call and its phase switch at one endpoint (genotype-equivalent): STABLE;
//        e  an endpoint tie that disagrees at a block: SENSITIVE there;
//        f  an endpoint whose exploration did not complete: NOT_EVALUATED, whatever the others say.

#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_model.hpp"
#include "panvar/genotype_repair.hpp"

#include <cstring>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

bool bitwise(double a, double b) { return std::memcmp(&a, &b, sizeof(double)) == 0; }

bool same_score(const panvar::DiplotypeScore& x, const panvar::DiplotypeScore& y) {
    return bitwise(x.log_likelihood, y.log_likelihood) && bitwise(x.fragment_term, y.fragment_term) &&
           bitwise(x.exposure_term, y.exposure_term) && bitwise(x.absolute_term_sum, y.absolute_term_sum) &&
           x.fragments_with_mass == y.fragments_with_mass;
}

std::string revcomp(const std::string& bases) {
    std::string out(bases.size(), 'N');
    for (std::size_t i = 0; i < bases.size(); ++i) {
        const char base = bases[bases.size() - 1 - i];
        out[i] = base == 'A' ? 'T' : base == 'C' ? 'G' : base == 'G' ? 'C' : 'A';
    }
    return out;
}

panvar::LibraryModel library(double intensity) {
    panvar::LibraryModel model;
    model.insert_mean_bp = 260.0;
    model.insert_sd_bp = 20.0;
    model.support_sigmas = 3.0;
    model.substitution_error = 0.02;
    model.outlier_weight = 0.005;
    model.fragment_start_intensity = intensity;
    model.provenance = "fixture: declared";
    return model;
}

panvar::FragmentScoringOptions options(std::size_t threads, const std::vector<panvar::LibraryModel>& audits = {}) {
    panvar::FragmentScoringOptions out;
    out.mismatch_band = 6;
    out.minimum_piece_bp = 10;
    out.threads = threads;
    out.audit_libraries = audits;
    return out;
}

std::string random_bases(std::mt19937_64& rng, std::size_t length) {
    static const char* alphabet = "ACGT";
    std::string out(length, 'A');
    for (char& base : out) base = alphabet[rng() % 4];
    return out;
}

void gate_exact_audit_scores() {
    std::mt19937_64 rng(20260929);
    const std::string backbone = random_bases(rng, 900);
    std::vector<std::string> haplotypes;
    for (int h = 0; h < 6; ++h) {
        std::string sequence = backbone;
        for (int edit = 0; edit < 2; ++edit) {
            const std::size_t at = (rng() % 850) + 10;
            sequence.replace(at, 12, random_bases(rng, 12));
        }
        haplotypes.push_back(sequence);
    }
    std::vector<panvar::ReadPair> fragments;
    for (int i = 0; i < 60; ++i) {
        const std::string& source = haplotypes[rng() % haplotypes.size()];
        const std::size_t start = rng() % (source.size() - 260);
        fragments.push_back({{"f" + std::to_string(i) + "/1", source.substr(start, 100), std::string(100, 'I')},
                             {"f" + std::to_string(i) + "/2", revcomp(source.substr(start + 160, 100)), std::string(100, 'I')}});
    }
    std::vector<panvar::MosaicPath> paths;
    for (std::size_t i = 0; i < haplotypes.size(); ++i) paths.push_back(panvar::MosaicPath{{i}});
    const panvar::InternedHaplotypes interned = panvar::intern_haplotypes(paths, haplotypes);
    std::vector<panvar::HaplotypeCandidate> pairs;
    for (std::size_t a = 0; a < interned.sequences.size(); ++a) {
        for (std::size_t b = a; b < interned.sequences.size(); ++b) pairs.push_back({a, b});
    }
    const double central = 0.0123456789;
    const std::vector<panvar::LibraryModel> audits{library(0.9 * central), library(1.1 * central)};
    for (const std::size_t threads : {std::size_t{1}, std::size_t{4}}) {
        const panvar::FragmentScoringResult plain =
            panvar::score_diplotypes(interned, pairs, fragments, library(central), options(threads));
        const panvar::FragmentScoringResult audited =
            panvar::score_diplotypes(interned, pairs, fragments, library(central), options(threads, audits));
        require(plain.audit_scores.empty() && audited.audit_scores.size() == 2, "one audit column per audit library");
        for (std::size_t c = 0; c < pairs.size(); ++c) {
            require(same_score(plain.scores[c], audited.scores[c]), "the audit changed a run's own score");
        }
        for (std::size_t k = 0; k < audits.size(); ++k) {
            const panvar::FragmentScoringResult declared =
                panvar::score_diplotypes(interned, pairs, fragments, audits[k], options(threads));
            for (std::size_t c = 0; c < pairs.size(); ++c) {
                require(same_score(declared.scores[c], audited.audit_scores[k][c]),
                        "a one-shot audit score is not the declared library's score");
            }
        }
        // The prepared context: the same, from retained columns.
        std::vector<panvar::DiplotypeScore> own;
        std::vector<std::vector<panvar::DiplotypeScore>> audit_columns;
        {
            panvar::PreparedScoringContext context(fragments, library(central), options(threads, audits));
            context.add_haplotypes(interned.sequences);
            const panvar::FragmentScoringResult scored = context.score_candidates(pairs);
            own = scored.scores;
            audit_columns = scored.audit_scores;
        }
        for (std::size_t c = 0; c < pairs.size(); ++c) {
            require(same_score(own[c], plain.scores[c]), "a context score changed under the audit");
        }
        for (std::size_t k = 0; k < audits.size(); ++k) {
            panvar::PreparedScoringContext context(fragments, audits[k], options(threads));
            context.add_haplotypes(interned.sequences);
            const panvar::FragmentScoringResult declared = context.score_candidates(pairs);
            for (std::size_t c = 0; c < pairs.size(); ++c) {
                require(same_score(declared.scores[c], audit_columns[k][c]),
                        "a context audit score is not the declared library's score");
            }
        }
    }
    std::cout << "PASS  audit scores are bit-identical to runs declaring each library (one-shot and context, 1 and 4 "
                 "threads); the run's own scores are unchanged\n";
}

// Two blocks. The call A = [0,0]|[1,1] carries 0/1 at both; its phase switch S = [0,1]|[1,0] carries the same
// unordered genotype as different sequences; B = [0,0]|[1,2] carries 0/2 at block 1.
panvar::PoolCandidate candidate(std::vector<std::size_t> first, std::vector<std::size_t> second, double total) {
    panvar::PoolCandidate c;
    std::string a = "seq", b = "seq";
    for (const std::size_t x : first) a += std::to_string(x);
    for (const std::size_t x : second) b += std::to_string(x);
    c.first = std::move(first);
    c.second = std::move(second);
    c.first_sequence = panvar::sequence_digest(a);
    c.second_sequence = panvar::sequence_digest(b);
    c.score.fragment_term = total;
    c.score.log_likelihood = total;
    c.score.absolute_term_sum = std::abs(total);
    c.score.fragments_with_mass = 10;
    return c;
}

panvar::PooledDecision decide(double a, double s, double b) {
    return panvar::decide_pooled_genotype({candidate({0, 0}, {1, 1}, a), candidate({0, 1}, {1, 0}, s),
                                           candidate({0, 0}, {1, 2}, b)},
                                          2, 10, 6);
}

void gate_classification() {
    const std::vector<std::pair<std::size_t, std::size_t>> call{{0, 1}, {0, 1}};
    const panvar::PooledDecision central = decide(-100.0, -103.0, -110.0);
    require(central.decision == panvar::TopSetDecision::Unique && central.genotype == call, "fixture: the call is A");
    const auto views = [&](const panvar::PooledDecision& low, const panvar::PooledDecision& high,
                           const std::string& low_anchor, const std::string& incomplete = "") {
        std::vector<panvar::LibraryAuditView> v(3);
        v[0].decision = &central;
        v[0].anchor = "anchor_A";
        v[1].decision = incomplete.empty() ? &low : nullptr;
        v[1].incomplete = incomplete;
        v[1].anchor = low_anchor;
        v[2].decision = &high;
        v[2].anchor = "anchor_A";
        return v;
    };
    // a: anchors differ, trajectories converge to the call
    {
        const panvar::PooledDecision low = decide(-100.0, -104.0, -106.0), high = decide(-100.0, -102.0, -112.0);
        const panvar::LibrarySensitivity r = panvar::classify_library_sensitivity(call, views(low, high, "anchor_B"));
        require(r.status == "STABLE" && r.sensitive_blocks.empty(), "a: converging trajectories are STABLE");
        require(r.distinct_anchors == 2 && !r.phase_or_representative_only && r.central_matches_primary,
                "a: two anchors recorded, nothing else changed");
    }
    // b: the low endpoint's top carries 0/2 at block 1
    {
        const panvar::PooledDecision low = decide(-100.0, -104.0, -95.0), high = decide(-100.0, -102.0, -112.0);
        const panvar::LibrarySensitivity r = panvar::classify_library_sensitivity(call, views(low, high, "anchor_A"));
        require(r.status == "SENSITIVE" && r.sensitive_blocks == std::vector<std::size_t>{1},
                "b: a differing endpoint genotype is SENSITIVE at exactly its block");
        require(r.call_class[1] == "strict_call" && r.projections[1][1].count({0, 2}) == 1, "b: the endpoint pair is recorded");
    }
    // c: the low endpoint's top is the phase switch
    {
        const panvar::PooledDecision low = decide(-100.0, -95.0, -110.0), high = decide(-100.0, -102.0, -112.0);
        const panvar::LibrarySensitivity r = panvar::classify_library_sensitivity(call, views(low, high, "anchor_A"));
        require(r.status == "STABLE" && r.phase_or_representative_only,
                "c: a phase-only change is STABLE and recorded as phase_or_representative_only");
        require(r.top_set_digest[1] != r.top_set_digest[0], "c: the top-set identities differ");
    }
    // d: the call and its phase switch tie at the high endpoint (genotype-equivalent)
    {
        const panvar::PooledDecision low = decide(-100.0, -104.0, -110.0), high = decide(-100.0, -100.0, -112.0);
        require(high.decision == panvar::TopSetDecision::GenotypeEquivalent, "fixture d: a genotype-equivalent tie");
        const panvar::LibrarySensitivity r = panvar::classify_library_sensitivity(call, views(low, high, "anchor_A"));
        require(r.status == "STABLE" && r.call_class[2] == "genotype_equivalent_call" && r.phase_or_representative_only,
                "d: a genotype-equivalent endpoint is STABLE");
    }
    // e: an endpoint tie between the call and B disagrees at block 1
    {
        const panvar::PooledDecision low = decide(-100.0, -104.0, -100.0), high = decide(-100.0, -102.0, -112.0);
        require(low.decision == panvar::TopSetDecision::DistinctTie, "fixture e: a distinct tie");
        const panvar::LibrarySensitivity r = panvar::classify_library_sensitivity(call, views(low, high, "anchor_A"));
        require(r.status == "SENSITIVE" && r.sensitive_blocks == std::vector<std::size_t>{1} && r.call_class[1] == "no_call",
                "e: an endpoint distinct tie is SENSITIVE at its disagreeing block");
    }
    // f: an incomplete endpoint makes the verdict NOT_EVALUATED even when the other endpoint differs
    {
        const panvar::PooledDecision low = decide(-100.0, -104.0, -110.0), high = decide(-100.0, -102.0, -90.0);
        const panvar::LibrarySensitivity r =
            panvar::classify_library_sensitivity(call, views(low, high, "anchor_A", "tie_branching_incomplete"));
        require(r.status == "NOT_EVALUATED" && r.sensitive_blocks.empty() && r.reason.find("tie_branching_incomplete") != std::string::npos,
                "f: an incomplete view is NOT_EVALUATED");
    }
    std::cout << "PASS  classification: converging anchors STABLE, a differing block SENSITIVE (exactly that block), "
                 "phase-only and genotype-equivalent STABLE, a distinct endpoint tie SENSITIVE, incomplete NOT_EVALUATED\n";
}

} // namespace

int main() {
    try {
        gate_exact_audit_scores();
        gate_classification();
    } catch (const std::exception& error) {
        std::cerr << "FAIL  " << error.what() << '\n';
        return 1;
    }
    return 0;
}
