#include "panvar/genotype_fragments.hpp"

#include "panvar/chain_kernel.hpp"

#include "panvar/md5.hpp"

#include "panvar/align.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/parallel.hpp"
#include "panvar/syncmer.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <mutex>
#include <stdexcept>
#include <array>
#include <cstring>
#include <unordered_map>
#include <unordered_set>

#include <zlib.h>
#include <kseq.h>

#include "edlib.h"

KSEQ_INIT(gzFile, gzread)

namespace panvar {

namespace {
// The per-mate alignment band, in edits. Defined once because THREE places need to agree on it: the
// block-local stage's partial-credit term, the haplotype scorer's placement band, and the
// band-boundary floor. The +1 matters -- two 150 bp mates at 5% give 8 + 8 = 16, not
// floor(0.05 * 300) = 15 -- and a floor computed from the fragment's total length instead of
// per mate is a different number that merely looks like the same one.
} // namespace

double InsertPrior::exposure(std::size_t hap_len) const {
    const long n = static_cast<long>(hap_len);
    double e = 0.0;
    for (long L = lo; L <= hi; ++L) {
        const double starts = static_cast<double>(std::max<long>(0, n - L + 1));
        if (starts > 0.0) e += std::exp(log_at(L)) * starts;
    }
    return e;
}

long fragment_insert_floor(const Fragment& f, bool allow_overlap) {
    return allow_overlap
        ? static_cast<long>(std::max(f.r1.size(), f.r2.size()))
        : static_cast<long>(f.r1.size() + f.r2.size());
}

long fragment_insert_floor(const std::vector<Fragment>& fragments, bool allow_overlap) {
    long m = 1;
    for (const Fragment& f : fragments)
        m = std::max<long>(m, fragment_insert_floor(f, allow_overlap));
    return m;
}

InsertPrior make_insert_prior(double mean, double sd, double discordant_rate,
                              int sigmas, long min_len) {
    InsertPrior ip;
    ip.lo = std::max<long>(min_len, static_cast<long>(mean - sigmas * sd));
    ip.hi = static_cast<long>(mean + sigmas * sd);
    if (ip.hi < ip.lo) ip.hi = ip.lo;
    const double span = static_cast<double>(ip.hi - ip.lo + 1);
    std::vector<double> w;
    double total = -std::numeric_limits<double>::infinity();
    const auto ladd = [](double a, double b) {
        if (a == -std::numeric_limits<double>::infinity()) return b;
        if (b == -std::numeric_limits<double>::infinity()) return a;
        const double hi2 = a > b ? a : b, lo2 = a > b ? b : a;
        return hi2 + std::log1p(std::exp(lo2 - hi2));
    };
    for (long L = ip.lo; L <= ip.hi; ++L) {
        const double z = (static_cast<double>(L) - mean) / sd;
        const double conc = std::log1p(-discordant_rate) - 0.5 * z * z
                          - std::log(sd) - 0.9189385332046727;
        const double disc = std::log(discordant_rate) - std::log(span);
        const double v = ladd(conc, disc);
        w.push_back(v);
        total = ladd(total, v);
    }
    for (const double v : w) ip.logp.push_back(v - total);
    // WHAT THE SUPPORT LEAVES OUT. The concordant component is Gaussian, so the residual beyond
    // each endpoint is a tail probability; the discordant component is uniform on the support and
    // has no tail by construction. Reported, not assumed away.
    const auto log_upper_tail = [&](double x) {
        const double z = (x - mean) / sd;
        if (z <= -8.0) return 0.0;
        const double q = 0.5 * std::erfc(z / std::sqrt(2.0));
        return q > 0.0 ? std::log(q) : -std::numeric_limits<double>::infinity();
    };
    const double lconc = std::log1p(-discordant_rate);
    ip.log_residual_above = lconc + log_upper_tail(static_cast<double>(ip.hi) + 0.5);
    const double zlo = (static_cast<double>(ip.lo) - 0.5 - mean) / sd;
    const double qlo = 0.5 * std::erfc(-zlo / std::sqrt(2.0));
    ip.log_residual_below =
        qlo > 0.0 ? lconc + std::log(qlo) : -std::numeric_limits<double>::infinity();
    return ip;
}

namespace {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();

// Placements of one mate that share a join coordinate. Two placements of the same mate with the
// same start are the same (start, L, strand) state; the Cartesian product enumerates both, so this
// sums their mass rather than deduplicating, and the two agree by construction instead of one of
// them quietly being the more correct model.
struct CoordAgg;

double log_add(double a, double b) {
    if (a == kNegInf) return b;
    if (b == kNegInf) return a;
    const double hi = a > b ? a : b;
    const double lo = a > b ? b : a;
    return hi + std::log1p(std::exp(lo - hi));
}

}  // namespace

// wgsim and `samtools fastq` both give the two mates one shared name with a /1 or /2 suffix, so the
// pairing can be recovered from the name alone -- which means interleaved input and split R1/R2
// files go through the same code and neither has to be declared.
//
// EXPORTED because the read counter has to agree with it exactly. Marker exclusion identifies a
// read by the FRAGMENT it belongs to, so a second copy of this rule in the counter would silently
// exclude the wrong reads the moment either drifted -- the same duplicated-derivation defect the
// shared linkage geometry exists to prevent.
std::string fragment_name(const char* raw, std::size_t len) {
    std::string name(raw, len);
    const std::size_t sp = name.find_first_of(" \t");
    if (sp != std::string::npos) name.resize(sp);
    if (name.size() > 2 && name[name.size() - 2] == '/' &&
        (name.back() == '1' || name.back() == '2')) {
        name.resize(name.size() - 2);
    }
    return name;
}

namespace {

// One read against one candidate context. `edits` is the infix (Hamming-window) edit distance, so
// the read may sit anywhere in the context with no end-gap penalty; `start` is where it landed,
// which is what the insert-size term needs.
struct ReadFit {
    bool ok = false;
    std::size_t edits = 0;
    std::size_t start = 0;
    std::size_t end = 0;
};

ReadFit infix_align(const std::string& read, const std::string& context, std::size_t max_edits) {
    ReadFit out;
    if (read.empty() || context.empty()) return out;
    const EdlibAlignConfig cfg = edlibNewAlignConfig(
        static_cast<int>(max_edits), EDLIB_MODE_HW, EDLIB_TASK_LOC, nullptr, 0);
    EdlibAlignResult res = edlibAlign(read.data(), static_cast<int>(read.size()),
                                      context.data(), static_cast<int>(context.size()), cfg);
    if (res.status == EDLIB_STATUS_OK && res.editDistance >= 0 && res.numLocations > 0) {
        out.ok = true;
        out.edits = static_cast<std::size_t>(res.editDistance);
        out.end = static_cast<std::size_t>(res.endLocations[0]);
        out.start = res.startLocations != nullptr
            ? static_cast<std::size_t>(res.startLocations[0])
            : (out.end >= read.size() ? out.end - read.size() + 1 : 0);
    }
    edlibFreeAlignResult(res);
    return out;
}

// Bounded: a Gaussian on the implied insert has unbounded influence and one mis-anchored pair can
// outvote thousands of ordinary ones. Mixed against a uniform, the worst a pair can say is
// log(discordant_rate / discordant_span).
// The SHARED normalised discrete prior, so the two scorers agree by construction. The old form was a
// continuous Gaussian plus a uniform term divided by a fixed span, which is neither normalised nor
// the same distribution the reference integrates.
double insert_ll(double implied, const InsertPrior& ip) {
    return ip.log_at(static_cast<long>(implied + 0.5));
}

struct CoordAgg {
    std::map<long, std::pair<double, double>> at;   // coordinate -> (summed log mass, best single)
    void add(long k, double v) {
        const auto it = at.find(k);
        if (it == at.end()) at.emplace(k, std::make_pair(v, v));
        else {
            it->second.first = log_add(it->second.first, v);
            if (v > it->second.second) it->second.second = v;
        }
    }
};

// The neighbouring-block sequence on each side of a block, taken from the allele the most panel
// haplotypes carry. It is identical for every candidate of the block, so it cannot shift any
// pairwise comparison -- its only job is to give a boundary-spanning fragment somewhere to land.
// That is precisely the evidence a marker built inside `allele_seq` cannot carry, and the reason an
// allele shorter than k is invisible to the marker model but visible here.
std::string majority_allele(const BlockAlleles& b) {
    if (b.allele_seq.empty()) return {};
    std::size_t best = 0;
    std::size_t best_n = 0;
    for (std::size_t i = 0; i < b.allele_seq.size(); ++i) {
        const std::size_t n = i < b.allele_haplotypes.size() ? b.allele_haplotypes[i] : 0;
        if (n > best_n) { best_n = n; best = i; }
    }
    return b.allele_seq[best];
}

// The locus spelled once through the majority allele of every block: a single sequence standing in
// for "everywhere else this fragment could have come from". Cheap, and it does not need to be the
// sample's own haplotype -- it only has to be close enough that a fragment from a paralogous copy
// finds a better home there than in the block under test.
std::string locus_consensus(const std::vector<BlockAlleles>& blocks,
                            std::vector<std::size_t>& offset,
                            std::vector<std::size_t>& length) {
    std::string acc;
    offset.assign(blocks.size(), 0);
    length.assign(blocks.size(), 0);
    for (std::size_t i = 0; i < blocks.size(); ++i) {
        const std::string m = majority_allele(blocks[i]);
        offset[i] = acc.size();
        length[i] = m.size();
        acc += m;
    }
    return acc;
}

std::string left_flank(const std::vector<BlockAlleles>& blocks, std::size_t bi, std::size_t want) {
    std::string acc;
    std::size_t i = bi;
    while (i > 0 && acc.size() < want) {
        --i;
        acc = majority_allele(blocks[i]) + acc;
    }
    if (acc.size() > want) acc.erase(0, acc.size() - want);
    return acc;
}

std::string right_flank(const std::vector<BlockAlleles>& blocks, std::size_t bi, std::size_t want) {
    std::string acc;
    for (std::size_t i = bi + 1; i < blocks.size() && acc.size() < want; ++i) {
        acc += majority_allele(blocks[i]);
    }
    if (acc.size() > want) acc.resize(want);
    return acc;
}

} // namespace

std::vector<Fragment> load_fragments(const std::vector<std::string>& paths,
                                     FragmentLoadStats* stats) {
    std::unordered_map<std::string, std::size_t> index;
    std::vector<Fragment> out;
    FragmentLoadStats st;

    for (const std::string& path : paths) {
        gzFile fp = gzopen(path.c_str(), "r");
        if (fp == nullptr) throw std::runtime_error("genotype-frag: cannot open reads " + path);
        kseq_t* seq = kseq_init(fp);
        while (kseq_read(seq) >= 0) {
            ++st.reads;
            std::string name = fragment_name(seq->name.s, seq->name.l);
            std::string bases(seq->seq.s, seq->seq.l);
            const auto it = index.find(name);
            if (it == index.end()) {
                index.emplace(name, out.size());
                Fragment f;
                f.name = std::move(name);
                f.r1 = std::move(bases);
                out.push_back(std::move(f));
            } else if (out[it->second].r2.empty()) {
                out[it->second].r2 = std::move(bases);
            } else {
                // A third read under one name. Not silently dropped: it means the input is not what
                // the pairing assumes, and every fragment count downstream would be wrong.
                ++st.over_paired;
            }
        }
        kseq_destroy(seq);
        gzclose(fp);
    }

    st.fragments = out.size();
    for (const Fragment& f : out) {
        if (f.r2.empty()) ++st.singleton; else ++st.paired;
    }
    if (stats != nullptr) *stats = st;
    return out;
}

std::vector<BlockFragmentResult> genotype_fragments(
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    const std::vector<Fragment>& fragments,
    const std::vector<std::size_t>& targets,
    const FragmentScoreOptions& options,
    const std::vector<int>* truth_allele1,
    const std::vector<int>* truth_allele2,
    std::size_t top_pairs_kept,
    double tie_eps) {

    const std::size_t k = options.kmer_size;
    const std::size_t s = options.syncmer_s != 0 ? options.syncmer_s : default_syncmer_s(k);

    // ---- contexts: flank + candidate allele + flank, one per surviving candidate -------------
    struct Target {
        std::size_t block = 0;
        std::string lf, rf;
        std::vector<std::size_t> candidates;     // panel allele indices, after coarse pruning
        std::vector<std::string> contexts;
        std::vector<std::size_t> allele_offset;  // where the allele starts inside the context
    };
    std::vector<Target> tg(targets.size());
    for (std::size_t t = 0; t < targets.size(); ++t) {
        tg[t].block = targets[t];
        tg[t].lf = left_flank(blocks, targets[t], options.flank_bp);
        tg[t].rf = right_flank(blocks, targets[t], options.flank_bp);
    }

    // ---- recruitment, and the factor each fragment belongs to ---------------------------------
    // ONE implementation, shared with the reference oracle and the incidence dump. A fragment is NOT
    // forced onto one block: at cyp2d6 the same syncmer occurs in blocks 3 and 5, and the fragment is
    // offered to both so its alignment can decide.
    const FragmentFactorIncidence incidence =
        classify_fragment_factors(blocks, targets, fragments, k, s, options.flank_bp,
                                  options.min_recruit_hits);
    const std::vector<std::vector<std::uint32_t>>& recruited = incidence.recruited;

    // ---- fragment -> factor incidence ---------------------------------------------------------
    // Recruitment offers a fragment to every block that could explain it, so the sets below are the
    // raw evidence relation, not yet an assignment. The factor a fragment BELONGS to follows from
    // the shape of that set:
    //
    //   0 blocks              unrecruited -- no factor, and it must not silently vanish
    //   1 block               a local block emission
    //   2 adjacent blocks     a transition factor across their boundary
    //   >2 consecutive        a higher-order path factor
    //   non-consecutive       ambiguous: shared evidence across repeated blocks, which must be
    //                         retained as shared rather than arbitrarily assigned or double counted
    //
    // Emitted so the "exactly once" invariant is checkable: one row per fragment, always.
    if (!options.incidence_path.empty()) {
        std::ofstream inc(options.incidence_path);
        if (!inc) throw std::runtime_error("genotype-frag: cannot write " + options.incidence_path);
        inc << "fragment\tn_blocks\tblocks\tfactor\tblock_a\tblock_b\n";
        const auto name_of = [](FragmentFactor f) {
            switch (f) {
                case FragmentFactor::Local:     return "local";
                case FragmentFactor::Boundary:  return "boundary";
                case FragmentFactor::Path:      return "path";
                case FragmentFactor::Ambiguous: return "ambiguous";
                default:                        return "unrecruited";
            }
        };
        for (std::size_t fi = 0; fi < fragments.size(); ++fi) {
            const std::vector<std::uint32_t>& tv = incidence.targets_of[fi];
            inc << fragments[fi].name << '\t' << tv.size() << '\t';
            for (std::size_t i = 0; i < tv.size(); ++i) inc << (i ? "," : "") << targets[tv[i]];
            if (tv.empty()) inc << '.';
            inc << '\t' << name_of(incidence.factor_of[fi]) << '\t'
                << (tv.empty() ? std::string(".") : std::to_string(targets[tv.front()])) << '\t'
                << (tv.size() < 2 ? std::string(".") : std::to_string(targets[tv.back()])) << '\n';
        }
        inc.flush();
        if (!inc) throw std::runtime_error("genotype-frag: write failed for " + options.incidence_path);
    }

    // ---- the rest of the locus, for the competitive background --------------------------------
    std::vector<std::size_t> cons_off, cons_len;
    const std::string consensus =
        options.compete ? locus_consensus(blocks, cons_off, cons_len) : std::string();

    // ---- coarse candidate selection, then fragment-level scoring ------------------------------
    std::vector<BlockFragmentResult> out(tg.size());

    run_parallel(tg.size(), options.threads, [&](std::size_t t) {
        Target& T = tg[t];
        const std::size_t bi = T.block;
        const BlockAlleles& B = blocks[bi];
        BlockFragmentResult& R = out[t];
        R.block_index = bi;
        R.kind = chain[bi].kind;
        R.bubble_id = chain[bi].bubble_id;
        R.n_alleles = B.allele_seq.size();
        R.n_fragments = recruited[t].size();

        if (B.allele_seq.empty()) return;

        // The codes the recruited fragments actually carry, used only to shortlist candidates. This
        // is the cheap vector-similarity stage -- a candidate generator, never the score.
        std::unordered_set<std::uint64_t> frag_codes;
        for (const std::uint32_t fi : recruited[t]) {
            for (const KmerOccurrence& o : collect_syncmers(fragments[fi].r1, k, s)) frag_codes.insert(o.code);
            for (const KmerOccurrence& o : collect_syncmers(fragments[fi].r2, k, s)) frag_codes.insert(o.code);
        }
        std::vector<std::pair<double, std::size_t>> ranked;
        ranked.reserve(B.allele_seq.size());
        for (std::size_t a = 0; a < B.allele_seq.size(); ++a) {
            const std::vector<KmerOccurrence> sy = collect_syncmers(B.allele_seq[a], k, s);
            double containment = 1.0;   // an allele with too few syncmers cannot be ranked coarsely
            if (sy.size() >= 3) {
                std::size_t hit = 0;
                for (const KmerOccurrence& o : sy) if (frag_codes.count(o.code)) ++hit;
                containment = static_cast<double>(hit) / static_cast<double>(sy.size());
            }
            ranked.emplace_back(containment, a);
        }
        std::stable_sort(ranked.begin(), ranked.end(),
                         [](const auto& x, const auto& y) { return x.first > y.first; });
        const std::size_t keep = std::min(options.max_alleles, ranked.size());
        for (std::size_t i = 0; i < keep; ++i) T.candidates.push_back(ranked[i].second);
        std::sort(T.candidates.begin(), T.candidates.end());
        R.n_candidates = T.candidates.size();

        T.contexts.resize(T.candidates.size());
        T.allele_offset.assign(T.candidates.size(), T.lf.size());
        for (std::size_t c = 0; c < T.candidates.size(); ++c) {
            T.contexts[c] = T.lf + B.allele_seq[T.candidates[c]] + T.rf;
        }

        const std::size_t nc = T.candidates.size();
        const std::size_t nf = recruited[t].size();
        if (nc == 0) return;

        // The locus with THIS block cut out. A fragment that belongs to the block has nowhere to go
        // in it; one that belongs to a paralogous copy elsewhere does, and is absorbed.
        std::string bg_ref;
        if (options.compete && !consensus.empty()) {
            bg_ref = consensus.substr(0, cons_off[bi]);
            bg_ref += consensus.substr(cons_off[bi] + cons_len[bi]);
        }

        // Same convention as the haplotype scorer and the reference: a specific mismatch is eps/3.
        const double log_eps = std::log(options.error_rate / 3.0);
        const double log_1meps = std::log1p(-options.error_rate);
        const double log_mix = std::log1p(-options.outlier_mix);
        const double log_out = std::log(options.outlier_mix);
        // THE SHARED RULE, over this haplotype's recruited subset. Computing the floor inline
        // here is how the bounded search kept the |r1| + |r2| floor after every other path moved
        // to max(|r1|, |r2|) -- the genotype_bounded_search fixture is what caught it, reporting
        // a support of [300, 650] where the rest of the model uses [150, 650].
        long min_frag_len_b = 1;
        for (const std::uint32_t fi : recruited[t]) {
            min_frag_len_b = std::max<long>(
                min_frag_len_b,
                fragment_insert_floor(fragments[fi], options.allow_overlapping_pairs));
        }
        const InsertPrior ins_prior = make_insert_prior(options.fragment_len, options.fragment_sd,
                                                        options.discordant_rate, options.insert_sigmas, min_frag_len_b);
        const auto read_ll = [&](std::size_t edits, std::size_t len) {
            return static_cast<double>(edits) * log_eps +
                   static_cast<double>(len - std::min(edits, len)) * log_1meps;
        };

        // [fragment][candidate] log P(fragment | candidate). One row per PHYSICAL fragment: this is
        // the whole point. A fragment enters the likelihood exactly once, so nothing downstream has
        // to discount for markers that shared a read.
        // RAW alignment log-likelihood, no outlier mixing: the mixture and the length normalisation
        // are properties of a PAIR, not of a candidate, so they are applied in the pair loop.
        std::vector<double> ll(nf * nc, kNegInf);
        std::vector<double> bg(nf, kNegInf);
        // Fragment start positions each context offers. The flank contributes to every candidate
        // equally, so what survives into a comparison is the allele-length difference alone.
        std::vector<double> log_len(nc, 0.0);
        for (std::size_t c = 0; c < nc; ++c) {
            const double span = static_cast<double>(T.contexts[c].size());
            const double occupied = options.use_insert_size ? options.fragment_len : 150.0;
            log_len[c] = std::log(std::max(1.0, span - occupied + 1.0));
        }
        const bool want_debug = options.debug_block >= 0 &&
                                static_cast<std::size_t>(options.debug_block) == bi &&
                                !options.debug_path.empty();
        std::vector<std::uint32_t> dbg_edits(want_debug ? nf * nc : 0, 0);

        for (std::size_t fi = 0; fi < nf; ++fi) {
            const Fragment& F = fragments[recruited[t][fi]];
            const std::size_t total_len = F.bases();
            if (total_len == 0) continue;
            const double floor_ll =
                read_ll(static_cast<std::size_t>(options.bg_divergence * static_cast<double>(total_len)),
                        total_len);

            // Orientation is a property of the fragment, not of the candidate, so it is resolved
            // once against the longest context and reused. Halves the alignment work.
            std::size_t probe = 0;
            for (std::size_t c = 1; c < nc; ++c) {
                if (T.contexts[c].size() > T.contexts[probe].size()) probe = c;
            }
            const std::string r1rc = reverse_complement(F.r1);
            const std::string r2rc = F.r2.empty() ? std::string() : reverse_complement(F.r2);
            const std::size_t band1 = mate_band_edits(options.max_divergence, F.r1.size());
            const std::size_t band2 = mate_band_edits(options.max_divergence, F.r2.size());

            const ReadFit p_fwd = infix_align(F.r1, T.contexts[probe], band1);
            const ReadFit p_rev = infix_align(r1rc, T.contexts[probe], band1);
            const bool r1_forward =
                (p_fwd.ok ? p_fwd.edits : band1 + 1) <= (p_rev.ok ? p_rev.edits : band1 + 1);
            const std::string& q1 = r1_forward ? F.r1 : r1rc;
            const std::string& q2 = r1_forward ? r2rc : F.r2;   // mates are on opposite strands

            // SCOPE OF THE "same model" CLAIM. The exhaustive reference and the final Hamming
            // haplotype scorer now score the same paired-fragment likelihood: only complete,
            // orientation-compatible mate placements. THIS STAGE DOES NOT. It accumulates each
            // mate's contribution independently and substitutes the band-boundary term for a mate
            // that fails to align, which is the very partial-credit state removed from the final
            // scorer -- a fragment with one aligning mate still contributes here.
            //
            // That is defensible as a PERMISSIVE SHORTLIST heuristic: this stage decides which
            // candidates are scored at all, and admitting too many costs time while admitting too
            // few loses the answer irrecoverably. It is not defensible as a likelihood, and nothing
            // downstream may treat its numbers as one. The consequence for evaluation is concrete:
            // a donor failure can arise HERE, by the floor pair never reaching the shortlist, and
            // that is a different defect from the final scorer ranking a shortlisted pair wrongly.
            // Every donor result must therefore report whether a sequence-floor pair survived the
            // shortlist, or the two failure modes are indistinguishable in the outcome.
            for (std::size_t c = 0; c < nc; ++c) {
                const std::string& ctx = T.contexts[c];
                const ReadFit f1 = infix_align(q1, ctx, band1);
                double lp = 0.0;
                bool any = false;
                if (f1.ok) { lp += read_ll(f1.edits, F.r1.size()); any = true; }
                else lp += read_ll(band1, F.r1.size());
                ReadFit f2;
                if (!F.r2.empty()) {
                    f2 = infix_align(q2, ctx, band2);
                    if (f2.ok) { lp += read_ll(f2.edits, F.r2.size()); any = true; }
                    else lp += read_ll(band2, F.r2.size());
                }
                if (options.use_insert_size && f1.ok && f2.ok) {
                    // The candidates differ in length, so the SAME physical fragment implies a
                    // different insert on each of them. That is how a 5-42 bp allele -- one that
                    // carries no k-mer of its own and is therefore invisible to the marker model --
                    // still gets voted on.
                    const std::size_t lo = std::min(f1.start, f2.start);
                    const std::size_t hi = std::max(f1.end, f2.end);
                    lp += insert_ll(static_cast<double>(hi - lo + 1), ins_prior);
                }
                if (want_debug) {
                    const std::size_t e = (f1.ok ? f1.edits : band1) +
                                          (F.r2.empty() ? 0 : (f2.ok ? f2.edits : band2));
                    dbg_edits[fi * nc + c] = static_cast<std::uint32_t>(e);
                }
                if (!any) lp = floor_ll;
                // Bounded loss at fragment level: no single fragment may express a preference
                // stronger than the likelihood ratio of a `bg_divergence` read.
                ll[fi * nc + c] = lp;
            }

            double bgl = floor_ll;
            if (!bg_ref.empty()) {
                const ReadFit b1 = infix_align(q1, bg_ref, band1);
                const ReadFit b2 = F.r2.empty() ? ReadFit{} : infix_align(q2, bg_ref, band2);
                double lb = 0.0;
                if (b1.ok) lb += read_ll(b1.edits, F.r1.size()); else lb += read_ll(band1, F.r1.size());
                if (!F.r2.empty()) {
                    if (b2.ok) lb += read_ll(b2.edits, F.r2.size());
                    else lb += read_ll(band2, F.r2.size());
                }
                if (b1.ok || b2.ok) bgl = std::max(bgl, lb);
            }
            bg[fi] = bgl;
        }

        if (want_debug) {
            std::ofstream df(options.debug_path);
            if (!df) throw std::runtime_error("genotype-frag: cannot write " + options.debug_path);
            df << "# block " << bi << " candidates";
            for (const std::size_t a : T.candidates) df << ' ' << a;
            df << "\n# allele_bp";
            for (const std::size_t a : T.candidates) df << ' ' << B.allele_seq[a].size();
            df << "\n# flank_bp " << T.lf.size() << ' ' << T.rf.size() << "\n";
            // The contexts themselves, so the alignments can be reproduced outside this binary.
            // A score that cannot be checked against the sequence it came from is not evidence.
            std::ofstream cf(options.debug_path + ".contexts.fa");
            if (!cf) throw std::runtime_error("genotype-frag: cannot write context FASTA");
            for (std::size_t c = 0; c < nc; ++c) {
                cf << ">block" << bi << "_cand" << c << "_allele" << T.candidates[c]
                   << " allele_bp=" << B.allele_seq[T.candidates[c]].size()
                   << " lf=" << T.lf.size() << " rf=" << T.rf.size() << '\n'
                   << T.contexts[c] << '\n';
            }
            df << "fragment\tname\tcandidate\tallele\tedits\tloglik\tbackground\n";
            for (std::size_t fi = 0; fi < nf; ++fi) {
                for (std::size_t c = 0; c < nc; ++c) {
                    df << fi << '\t' << fragments[recruited[t][fi]].name << '\t' << c << '\t'
                       << T.candidates[c] << '\t' << dbg_edits[fi * nc + c] << '\t'
                       << ll[fi * nc + c] << '\t' << bg[fi] << '\n';
                }
            }
        }

        // Informative = the fragment's likelihood is not flat across candidates. A flat fragment
        // contributes the same constant to every pair and cancels; counting it as evidence is how a
        // block with no discrimination comes to look well covered.
        for (std::size_t fi = 0; fi < nf; ++fi) {
            double lo = kNegInf, hi = kNegInf;
            for (std::size_t c = 0; c < nc; ++c) {
                const double v = ll[fi * nc + c];
                if (hi == kNegInf || v > hi) hi = v;
                if (lo == kNegInf || v < lo) lo = v;
            }
            if (hi - lo > 1e-9) ++R.n_informative;
        }

        // ---- diploid pair likelihood: each fragment once, mixed over the two homologues ---------
        // P(fragment | a, b) = (1-eta) * (P_a + P_b) / (L_a + L_b)  +  eta * P_background
        //
        // Three things are doing work here and each is separable by a flag. The numerator is the
        // mixture over the two homologues, so a fragment is explained if EITHER carries it. The
        // denominator is the length normalisation, so a pair that offers more sequence must explain
        // proportionally more fragments. The background absorbs fragments the rest of the locus
        // explains better, so they express no preference between candidates here.
        std::vector<PairScore> pairs;
        pairs.reserve(nc * (nc + 1) / 2);
        for (std::size_t ca = 0; ca < nc; ++ca) {
            for (std::size_t cb = ca; cb < nc; ++cb) {
                const double log_total_len = options.length_normalize
                    ? log_add(log_len[ca], log_len[cb]) : std::log(2.0);
                double total = 0.0;
                for (std::size_t fi = 0; fi < nf; ++fi) {
                    const double num = log_add(ll[fi * nc + ca], ll[fi * nc + cb]);
                    total += log_add(log_mix + num - log_total_len, log_out + bg[fi]);
                }
                pairs.push_back({T.candidates[ca], T.candidates[cb], total});
            }
        }
        if (pairs.empty()) return;
        std::stable_sort(pairs.begin(), pairs.end(),
                         [](const PairScore& x, const PairScore& y) { return x.score > y.score; });
        R.best_a = pairs[0].allele1;
        R.best_b = pairs[0].allele2;
        R.best_score = pairs[0].score;
        for (const PairScore& p : pairs) if (pairs[0].score - p.score <= tie_eps) ++R.top_class;
        for (std::size_t i = 0; i < std::min(top_pairs_kept, pairs.size()); ++i) {
            R.top_pairs.push_back(pairs[i]);
        }

        if (truth_allele1 != nullptr && truth_allele2 != nullptr &&
            bi < truth_allele1->size() && bi < truth_allele2->size()) {
            R.truth_a = (*truth_allele1)[bi];
            R.truth_b = (*truth_allele2)[bi];
            if (R.truth_a < 0 || R.truth_b < 0) {
                R.truth_rank = -1;    // not representable in the reduced panel
            } else {
                const std::size_t ta = std::min<std::size_t>(R.truth_a, R.truth_b);
                const std::size_t tb = std::max<std::size_t>(R.truth_a, R.truth_b);
                const bool in_a = std::find(T.candidates.begin(), T.candidates.end(), ta) != T.candidates.end();
                const bool in_b = std::find(T.candidates.begin(), T.candidates.end(), tb) != T.candidates.end();
                if (!in_a || !in_b) {
                    R.truth_rank = -2;   // representable, but coarse pruning dropped it
                } else {
                    for (const PairScore& p : pairs) {
                        if (p.allele1 != ta || p.allele2 != tb) continue;
                        int better = 0, ties = 0;
                        for (const PairScore& q : pairs) {
                            if (q.score > p.score + tie_eps) ++better;
                            else if (std::abs(q.score - p.score) <= tie_eps) ++ties;
                        }
                        R.truth_rank = better + 1;
                        R.truth_ties = ties - 1;
                        R.truth_delta = p.score - pairs[0].score;
                        break;
                    }
                }
            }
        }
    });

    return out;
}

FragmentFactorIncidence classify_fragment_factors(const std::vector<BlockAlleles>& blocks,
                                                  const std::vector<std::size_t>& targets,
                                                  const std::vector<Fragment>& fragments,
                                                  std::size_t kmer_size,
                                                  std::size_t syncmer_s,
                                                  std::size_t flank_bp,
                                                  std::size_t min_recruit_hits) {
    const std::size_t k = kmer_size;
    const std::size_t s = syncmer_s != 0 ? syncmer_s : default_syncmer_s(k);
    FragmentFactorIncidence out;
    out.recruited.assign(targets.size(), {});
    out.targets_of.assign(fragments.size(), {});
    out.factor_of.assign(fragments.size(), FragmentFactor::Unrecruited);

    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> recruit;
    recruit.reserve(1u << 20);
    {
        std::vector<std::uint64_t> codes;
        for (std::size_t t = 0; t < targets.size(); ++t) {
            codes.clear();
            const BlockAlleles& b = blocks[targets[t]];
            const auto add = [&](const std::string& seq) {
                for (const KmerOccurrence& o : collect_syncmers(seq, k, s)) codes.push_back(o.code);
            };
            for (const std::string& a : b.allele_seq) add(a);
            add(left_flank(blocks, targets[t], flank_bp));
            add(right_flank(blocks, targets[t], flank_bp));
            std::sort(codes.begin(), codes.end());
            codes.erase(std::unique(codes.begin(), codes.end()), codes.end());
            for (const std::uint64_t c : codes) recruit[c].push_back(static_cast<std::uint32_t>(t));
        }
    }
    std::unordered_map<std::uint32_t, std::uint32_t> hits;
    for (std::size_t fi = 0; fi < fragments.size(); ++fi) {
        hits.clear();
        const auto scan = [&](const std::string& r) {
            for (const KmerOccurrence& o : collect_syncmers(r, k, s)) {
                const auto it = recruit.find(o.code);
                if (it == recruit.end()) continue;
                for (const std::uint32_t t : it->second) ++hits[t];
            }
        };
        scan(fragments[fi].r1);
        scan(fragments[fi].r2);
        for (const auto& [t, n] : hits) {
            if (n >= min_recruit_hits) {
                out.recruited[t].push_back(static_cast<std::uint32_t>(fi));
                out.targets_of[fi].push_back(t);
            }
        }
        auto& tv = out.targets_of[fi];
        std::sort(tv.begin(), tv.end());
        if (tv.empty()) {
            out.factor_of[fi] = FragmentFactor::Unrecruited;
        } else if (tv.size() == 1) {
            out.factor_of[fi] = FragmentFactor::Local;
        } else {
            const bool consecutive = (tv.back() - tv.front() + 1) == tv.size();
            out.factor_of[fi] = !consecutive ? FragmentFactor::Ambiguous
                              : (tv.size() == 2 ? FragmentFactor::Boundary : FragmentFactor::Path);
        }
    }
    for (auto& v : out.recruited) std::sort(v.begin(), v.end());
    return out;
}

std::string chain_left_flank(const std::vector<BlockAlleles>& blocks, std::size_t bi,
                             std::size_t want) {
    return left_flank(blocks, bi, want);
}

std::string chain_right_flank(const std::vector<BlockAlleles>& blocks, std::size_t bi,
                              std::size_t want) {
    return right_flank(blocks, bi, want);
}

double reference_factor_loglik(const std::string& hap_a,
                               const std::string& hap_b,
                               const std::vector<Fragment>& fragments,
                               const ReferenceParams& params) {
    return reference_pair_loglik(hap_a, hap_b, fragments, params);
}

std::string chain_span_sequence(const std::vector<BlockAlleles>& blocks,
                                const std::vector<std::size_t>& targets,
                                std::size_t first_target,
                                std::size_t last_target,
                                const std::vector<int>& alleles,
                                std::size_t flank_bp) {
    if (first_target > last_target || last_target >= targets.size()) {
        throw std::runtime_error("genotype-frag: chain_span_sequence: target range out of order or "
                                 "past the end of the target list");
    }
    if (alleles.size() != last_target - first_target + 1) {
        throw std::runtime_error("genotype-frag: chain_span_sequence: one allele per target in the "
                                 "span is required (" + std::to_string(alleles.size()) + " given for " +
                                 std::to_string(last_target - first_target + 1) + " targets)");
    }
    const std::size_t bfirst = targets[first_target], blast = targets[last_target];
    // Context OUTSIDE the span, from the chain. This is the only place a majority allele is used,
    // and only for sequence the factor does not model; everything inside the span is explicit.
    std::string out = left_flank(blocks, bfirst, flank_bp);
    for (std::size_t bi = bfirst; bi <= blast; ++bi) {
        // A target's own allele is explicit; INTERVENING blocks (backbone between two bubble
        // targets) take the majority allele, because they are not part of the factor's state and the
        // incidence table does not treat them as targets either.
        bool is_target = false;
        std::size_t which = 0;
        for (std::size_t t = first_target; t <= last_target; ++t) {
            if (targets[t] == bi) { is_target = true; which = t - first_target; break; }
        }
        if (is_target) {
            const int a = alleles[which];
            if (a >= 0 && static_cast<std::size_t>(a) < blocks[bi].allele_seq.size()) {
                out += blocks[bi].allele_seq[static_cast<std::size_t>(a)];
            }
        } else {
            // An INTERVENING block, between two targets. Taking its majority allele would put a
            // GUESSED sequence inside a state the caller believes it stated explicitly -- the same
            // defect as majority-flank guessing, one level in. It is only safe when the block is
            // invariant, so that is asserted rather than assumed.
            if (blocks[bi].allele_seq.size() > 1) {
                throw std::runtime_error(
                    "genotype-frag: chain_span_sequence: block " + std::to_string(bi) +
                    " lies between two targets and carries " +
                    std::to_string(blocks[bi].allele_seq.size()) + " alleles, so its sequence cannot "
                    "be inferred. Make it an explicit target (--all-blocks) rather than letting the "
                    "span guess it.");
            }
            out += majority_allele(blocks[bi]);
        }
    }
    out += right_flank(blocks, blast, flank_bp);
    return out;
}

void write_fragment_results(const std::string& out_prefix,
                            const std::vector<BlockFragmentResult>& results,
                            bool have_truth) {
    const std::string blocks_path = out_prefix + ".frag_blocks.tsv";
    std::ofstream bf(blocks_path);
    if (!bf) throw std::runtime_error("genotype-frag: cannot write " + blocks_path);
    bf << "block\tkind\tbubble_id\tn_alleles\tn_candidates\tn_fragments\tn_informative\t"
          "best_a\tbest_b\tbest_score\ttop_class";
    if (have_truth) bf << "\ttruth_a\ttruth_b\ttruth_rank\ttruth_ties\ttruth_delta\texact";
    bf << '\n';
    for (const BlockFragmentResult& r : results) {
        const char* kind = r.kind == BlockKind::Bubble ? "bubble"
                         : r.kind == BlockKind::Backbone ? "backbone" : "flank";
        bf << r.block_index << '\t' << kind << '\t' << r.bubble_id << '\t' << r.n_alleles << '\t'
           << r.n_candidates << '\t' << r.n_fragments << '\t' << r.n_informative << '\t'
           << r.best_a << '\t' << r.best_b << '\t' << r.best_score << '\t' << r.top_class;
        if (have_truth) {
            bf << '\t' << r.truth_a << '\t' << r.truth_b << '\t' << r.truth_rank << '\t'
               << r.truth_ties << '\t' << r.truth_delta << '\t';
            if (r.truth_rank < 0) bf << "NA";
            else {
                const std::size_t ta = std::min<std::size_t>(r.truth_a, r.truth_b);
                const std::size_t tb = std::max<std::size_t>(r.truth_a, r.truth_b);
                bf << ((ta == std::min(r.best_a, r.best_b) && tb == std::max(r.best_a, r.best_b))
                       ? "1" : "0");
            }
        }
        bf << '\n';
    }
    bf.flush();
    if (!bf) throw std::runtime_error("genotype-frag: write failed for " + blocks_path);

    const std::string pairs_path = out_prefix + ".frag_pairs.tsv";
    std::ofstream pf(pairs_path);
    if (!pf) throw std::runtime_error("genotype-frag: cannot write " + pairs_path);
    pf << "block\trank\tallele1\tallele2\tscore\tdelta\tis_truth\n";
    for (const BlockFragmentResult& r : results) {
        for (std::size_t i = 0; i < r.top_pairs.size(); ++i) {
            const PairScore& p = r.top_pairs[i];
            const bool is_truth =
                r.truth_a >= 0 && r.truth_b >= 0 &&
                p.allele1 == std::min<std::size_t>(r.truth_a, r.truth_b) &&
                p.allele2 == std::max<std::size_t>(r.truth_a, r.truth_b);
            pf << r.block_index << '\t' << (i + 1) << '\t' << p.allele1 << '\t' << p.allele2 << '\t'
               << p.score << '\t' << (p.score - r.best_score) << '\t' << (is_truth ? 1 : 0) << '\n';
        }
    }
    pf.flush();
    if (!pf) throw std::runtime_error("genotype-frag: write failed for " + pairs_path);
}


// =================================================================================================
// WHOLE-HAPLOTYPE MODE
// =================================================================================================

namespace {

// One panel haplotype spelled end to end across the block chain, with the offset of every block so a
// pair score can be projected back onto blocks. A haplotype that bypasses a block contributes zero
// bases there, which is what a deletion spanning the site actually is.
struct HaplotypeSeq {
    std::string seq;
    std::vector<int> allele;        // per block, this haplotype's allele index (-1 = absent)
};
}   // close the anonymous namespace so the spelling check can use these
using HaplotypeSeqPub = HaplotypeSeq;
namespace {

HaplotypeSeq spell_haplotype(const std::vector<BlockAlleles>& blocks, const std::string& name) {
    HaplotypeSeq h;
    h.allele.assign(blocks.size(), -1);
    for (std::size_t bi = 0; bi < blocks.size(); ++bi) {
        const auto it = blocks[bi].allele_of.find(name);
        if (it == blocks[bi].allele_of.end()) continue;
        h.allele[bi] = static_cast<int>(it->second);
        if (it->second < blocks[bi].allele_seq.size()) h.seq += blocks[bi].allele_seq[it->second];
    }
    return h;
}

} // namespace

std::string spell_block_haplotype(const std::vector<BlockAlleles>& blocks,
                                  const std::string& name) {
    return spell_haplotype(blocks, name).seq;
}

std::string allele_catalogue_fingerprint(const std::vector<BlockAlleles>& blocks) {
    // Per-allele digests are folded in rather than the sequences themselves, so the cost does not
    // depend on holding the whole panel in one buffer. Block and allele indices are included
    // because it is the INDEX that the call table stores: two catalogues holding the same set of
    // sequences in a different order must not agree.
    std::string acc;
    acc.reserve(blocks.size() * 64);
    for (std::size_t bi = 0; bi < blocks.size(); ++bi) {
        acc += std::to_string(bi);
        acc += ':';
        acc += std::to_string(blocks[bi].allele_seq.size());
        acc += ':';
        acc += std::to_string(blocks[bi].bypass_allele);
        for (std::size_t ai = 0; ai < blocks[bi].allele_seq.size(); ++ai) {
            acc += ',';
            acc += std::to_string(ai);
            acc += '=';
            acc += md5_hex(blocks[bi].allele_seq[ai]);
        }
        acc += ';';
    }
    return md5_hex(acc);
}

void verify_block_spelling(const Graph& graph,
                           const std::vector<BlockAlleles>& blocks,
                           const std::vector<std::string>& haplotype_names) {
    const auto by_name = path_records_by_name(graph);
    std::size_t gapped = 0;
    for (const std::string& name : haplotype_names) {
        const auto it = by_name.find(name);
        if (it == by_name.end() || it->second == nullptr) {
            // Skipping this would exempt exactly the paths least likely to be trustworthy.
            throw std::runtime_error(
                "genotype-frag: path '" + name + "' is scored as a panel haplotype but has no record "
                "in the graph, so its block spelling cannot be verified against anything");
        }
        bool complete = false;
        const std::string raw = spell_path_steps_sequence(graph, it->second->steps, &complete);
        if (!complete) {
            throw std::runtime_error(
                "genotype-frag: the graph cannot fully spell path '" + name + "' (a step has no "
                "sequence), so whole-haplotype mode would score a haplotype it cannot verify. "
                "Previously this path was silently skipped by the round-trip check");
        }
        const HaplotypeSeq built = spell_haplotype(blocks, name);
        if (built.seq == raw) continue;
        // A path running antiparallel to the reference spells reference-oriented blocks, so its
        // concatenation is the exact reverse complement of its walk. That is a FRAME difference and
        // not a decomposition fault: the sequence is the same, and reads are aligned in both
        // orientations anyway. Refusing it kept c4, cyp2d6 and ankrd36c out of every experiment.
        if (built.seq == reverse_complement(raw)) continue;
        // A decomposition gap. This used to be fatal, and rightly so while the SCORED sequence was
        // the block concatenation -- the caller would have scored a sequence the panel does not
        // contain. Now the graph walk is authoritative and the decomposition is used only to project
        // the answer onto blocks, so a gap costs a projection, not a wrong call. Reported loudly and
        // counted; the caller decides.
        std::size_t at = 0;
        while (at < built.seq.size() && at < raw.size() && built.seq[at] == raw[at]) ++at;
        std::fprintf(stderr,
            "[genotype-frag] WARNING: the block decomposition does not reproduce path '%s': "
            "concatenated block alleles are %zu bp against a walk of %zu bp, first difference at "
            "offset %zu. The WALK is scored, so the call is unaffected; the per-block projection for "
            "this path is incomplete.\n",
            name.c_str(), built.seq.size(), raw.size(), at);
        ++gapped;
    }
    if (gapped > 0) {
        std::fprintf(stderr,
            "[genotype-frag] %zu of %zu panel paths have an incomplete block projection (see above). "
            "Scoring is unaffected.\n", gapped, haplotype_names.size());
    }
}

SpellFrame block_spelling_frame(const Graph& graph,
                                const std::vector<BlockAlleles>& blocks,
                                const std::string& name) {
    const auto by_name = path_records_by_name(graph);
    const auto it = by_name.find(name);
    if (it == by_name.end() || it->second == nullptr) return SpellFrame::Incomparable;
    bool complete = false;
    const std::string raw = spell_path_steps_sequence(graph, it->second->steps, &complete);
    if (!complete) return SpellFrame::Incomparable;
    const std::string built = spell_haplotype(blocks, name).seq;
    if (built == raw) return SpellFrame::Forward;
    if (built == reverse_complement(raw)) return SpellFrame::ReverseComplement;
    return SpellFrame::Incomparable;
}

HaplotypeResult genotype_haplotype_pairs(
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    const std::vector<std::string>& haplotype_names,
    const std::vector<Fragment>& fragments,
    const HaplotypeScoreOptions& options,
    const std::vector<int>* truth_allele1,
    const std::vector<int>* truth_allele2,
    std::size_t top_pairs_kept,
    const std::unordered_map<std::string, std::string>* walk_sequences) {

    HaplotypeResult out;
    out.n_fragments = fragments.size();
    const std::size_t k = options.kmer_size;
    const std::size_t s = options.syncmer_s != 0 ? options.syncmer_s : default_syncmer_s(k);

    // ---- the reads' own syncmer set, for the coarse shortlist ---------------------------------
    std::unordered_set<std::uint64_t> read_codes;
    read_codes.reserve(fragments.size() * 32);
    for (const Fragment& f : fragments) {
        for (const KmerOccurrence& o : collect_syncmers(f.r1, k, s)) read_codes.insert(o.code);
        for (const KmerOccurrence& o : collect_syncmers(f.r2, k, s)) read_codes.insert(o.code);
    }

    // Coarse selection: what fraction of a haplotype's own syncmers the reads contain. This is the
    // cheap vector-similarity stage and nothing more -- its score never reaches the output. A
    // candidate generator that is wrong here loses the answer outright, so it is deliberately
    // generous, and whether the truth survived it is reported rather than assumed.
    // ONE source of sequence for every candidate. A run that scored some candidates from walks and
    // others from block concatenation would be comparing two different reconstructions inside one
    // likelihood, which is precisely the class of defect this change exists to remove -- so a walk
    // map that does not cover every scored path is refused rather than silently filled in.
    std::vector<const std::string*> hap_seq(haplotype_names.size(), nullptr);
    std::vector<std::string> block_seq_store;
    if (walk_sequences != nullptr) {
        std::vector<std::string> missing;
        for (std::size_t i = 0; i < haplotype_names.size(); ++i) {
            const auto it = walk_sequences->find(haplotype_names[i]);
            if (it == walk_sequences->end() || it->second.empty()) {
                if (missing.size() < 5) missing.push_back(haplotype_names[i]);
                continue;
            }
            hap_seq[i] = &it->second;
        }
        if (!missing.empty()) {
            std::string m;
            for (const std::string& x : missing) { m += "\n    "; m += x; }
            throw std::runtime_error(
                "genotype-frag: the graph cannot completely spell every panel path, so scoring would "
                "mix walk-derived and block-derived sequence inside one likelihood. Refusing. "
                "Paths without a complete walk:" + m);
        }
    } else {
        block_seq_store.resize(haplotype_names.size());
        for (std::size_t i = 0; i < haplotype_names.size(); ++i) {
            block_seq_store[i] = spell_haplotype(blocks, haplotype_names[i]).seq;
            hap_seq[i] = &block_seq_store[i];
        }
    }
    const auto seq_of = [&](std::size_t i) { return hap_seq[i]; };

    std::vector<std::pair<double, std::size_t>> ranked;
    ranked.reserve(haplotype_names.size());
    {
        std::vector<double> containment(haplotype_names.size(), 0.0);
        run_parallel(haplotype_names.size(), options.threads, [&](std::size_t i) {
            const std::vector<KmerOccurrence> sy = collect_syncmers(*seq_of(i), k, s);
            if (sy.empty()) return;
            std::size_t hit = 0;
            for (const KmerOccurrence& o : sy) if (read_codes.count(o.code)) ++hit;
            containment[i] = static_cast<double>(hit) / static_cast<double>(sy.size());
        });
        for (std::size_t i = 0; i < haplotype_names.size(); ++i) ranked.emplace_back(containment[i], i);
    }
    std::stable_sort(ranked.begin(), ranked.end(),
                     [](const auto& a, const auto& b) { return a.first > b.first; });

    // Containment SATURATES. It asks what fraction of a haplotype's syncmers the reads contain, and
    // a panel haplotype differing from the sample only in WHICH combination of otherwise-present
    // alleles it carries scores 1.0 like every other -- the discriminating signal there is linkage,
    // which a presence score cannot see. Two consequences, both measured in the registered fixture,
    // where the correct pair was reachable only at full panel size:
    //
    //   * sequence-identical haplotypes each consume a slot, and a panel holds many;
    //   * once scores tie, the cut is decided by path order, which is arbitrary.
    //
    // So: collapse identical sequences to one candidate, and never cut through a tie -- extend until
    // the score strictly separates. This is the shortlist arm of the gstm1 donors whose answer never
    // reaches scoring.
    {
        std::unordered_set<std::string> seen;
        std::vector<std::pair<double, std::size_t>> unique_ranked;
        unique_ranked.reserve(ranked.size());
        for (const auto& r : ranked) {
            // STRAND-CANONICAL key: two paths whose sequences are reverse complements of one another
            // are the same haplotype seen from opposite strands, and the scorer is strand-symmetric
            // (every fragment is aligned in both orientations). Keying on the raw sequence would keep
            // both and let them consume two shortlist slots while carrying identical evidence.
            const std::string& sq = *seq_of(r.second);
            const std::string rc = reverse_complement(sq);
            if (seen.insert(sq < rc ? sq : rc).second) unique_ranked.push_back(r);
        }
        ranked.swap(unique_ranked);
    }
    // Forced haplotypes are moved to the front of the ranking, so they survive any cut. Their
    // containment score is left untouched -- only their position changes -- because the score is
    // reported and must stay the recruiter's own opinion.
    if (!options.force_haplotypes.empty()) {
        std::vector<std::pair<double, std::size_t>> forced, rest;
        for (const auto& r : ranked) {
            const bool want = std::find(options.force_haplotypes.begin(),
                                        options.force_haplotypes.end(),
                                        haplotype_names[r.second]) != options.force_haplotypes.end();
            (want ? forced : rest).push_back(r);
        }
        ranked.clear();
        ranked.insert(ranked.end(), forced.begin(), forced.end());
        ranked.insert(ranked.end(), rest.begin(), rest.end());
    }
    std::size_t nh = std::min(options.max_haplotypes, ranked.size());
    while (nh > 0 && nh < ranked.size() &&
           std::abs(ranked[nh].first - ranked[nh - 1].first) < 1e-12) {
        ++nh;   // the cut fell inside a tie; taking one side of it would be arbitrary
    }
    std::vector<HaplotypeSeq> haps(nh);
    // PER CANDIDATE, PER BLOCK. The old rule was all-or-nothing: any difference between the block
    // concatenation and the walk marked EVERY block of that haplotype unreliable. That was defensible
    // while a mismatch meant the decomposition was wrong, but a path may now legitimately END INSIDE
    // A BLOCK, and its concatenation is then a correct PREFIX of its walk. Under the old rule such a
    // path poisons all nineteen of its blocks because of one truncated terminal one.
    std::vector<std::vector<char>> proj_block(nh);
    out.haplotypes.resize(nh);
    // Reported so a run can say whether the shortlist was a real selection or a tie it could not
    // break. `nh` above the requested cap means it was extended through a tie.
    for (std::size_t i = 0; i < nh; ++i) {
        out.shortlist.push_back(haplotype_names[ranked[i].second]);
        // Alleles from the decomposition (needed to project the answer onto blocks); SEQUENCE from
        // the graph walk when available, because the walk is the haplotype and the concatenation is
        // only an approximation of it.
        // Alleles from the decomposition (only to PROJECT the answer onto blocks); SEQUENCE from
        // seq_of, which is the walk when walks were supplied and refuses to be partial.
        haps[i] = spell_haplotype(blocks, haplotype_names[ranked[i].second]);
        const std::string block_concat = haps[i].seq;   // before it is replaced by the walk
        haps[i].seq = *seq_of(ranked[i].second);
        // Does the decomposition reproduce this haplotype? If not, its per-block alleles cannot be
        // relied on anywhere, and every block of a projection involving it is unprojectable.
        {
            // Projectability comes from the shared frame: a block is reliable for this candidate
            // when the verified map covers it. Blocks outside a partial frame's window are not.
            const PathProjection pr =
                project_path_blocks(blocks, blocks, haplotype_names[ranked[i].second], haps[i].seq);
            proj_block[i].assign(blocks.size(), 0);
            if (pr.ok) {
                for (const PathBlockSlice& sl : pr.blocks) {
                    if (sl.block < proj_block[i].size()) proj_block[i][sl.block] = 1;
                }
            } else {
                // No verified map at all: nothing about this candidate's blocks can be trusted.
                const bool whole = (block_concat == haps[i].seq) ||
                                   (block_concat == reverse_complement(haps[i].seq));
                if (whole) proj_block[i].assign(blocks.size(), 1);
            }
        }
        out.haplotypes[i].name = haplotype_names[ranked[i].second];
        out.haplotypes[i].bp = haps[i].seq.size();
        out.haplotypes[i].containment = ranked[i].first;
    }
    if (nh == 0) return out;

    // ---- anchor index: syncmer -> (haplotype, position) --------------------------------------
    // Recruitment and placement in one structure. A read is not assigned to a haplotype; it is
    // offered every position any of its syncmers points at, on every shortlisted haplotype, and the
    // alignment decides. Syncmers commoner than `max_anchor_occ` anchor nothing, because inside a
    // tandem array they point everywhere and cost time without adding information.
    std::unordered_map<std::uint64_t, std::vector<std::pair<std::uint32_t, std::uint32_t>>> anchors;
    anchors.reserve(1u << 22);
    for (std::size_t hi = 0; hi < nh; ++hi) {
        std::unordered_map<std::uint64_t, std::uint32_t> per_hap;
        const std::vector<KmerOccurrence> sy = collect_syncmers(haps[hi].seq, k, s);
        for (const KmerOccurrence& o : sy) ++per_hap[o.code];
        for (const KmerOccurrence& o : sy) {
            ++out.completeness.anchor_occurrences_seen;
            if (per_hap[o.code] > options.max_anchor_occ) {
                ++out.completeness.anchor_occurrences_dropped;
                continue;
            }
            anchors[o.code].emplace_back(static_cast<std::uint32_t>(hi),
                                         static_cast<std::uint32_t>(o.start));
        }
    }

    const double log_1meps = std::log1p(-options.error_rate);
    // CONVENTION: `error_rate` is the total substitution probability at a base, so a SPECIFIC
    // observed mismatch has probability error_rate/3 -- which is what the simulator draws, choosing
    // uniformly among the other three bases. Using log(error_rate) per mismatch would be the
    // probability of "some substitution", counted once per observed base, and overstates every
    // mismatched read by log(3) per edit. It cancels between the two scorers, so it never showed in a
    // differential, but it is wrong for a calibrated likelihood and therefore for GQ.
    const double log_eps3 = std::log(options.error_rate / 3.0);
    const auto read_ll = [&](std::size_t edits, std::size_t len) {
        return static_cast<double>(edits) * log_eps3 +
               static_cast<double>(len - std::min(edits, len)) * log_1meps;
    };
    const long min_frag_len = fragment_insert_floor(fragments, options.allow_overlapping_pairs);
    const InsertPrior ins_prior = make_insert_prior(options.fragment_len, options.fragment_sd,
                                                    options.discordant_rate, options.insert_sigmas, min_frag_len);
    // Explicit 1/2 per strand, matching the reference. It is a per-fragment constant only while the
    // placement term is compared with itself; mixed against a background it changes the placement
    // scale and therefore the contrast between candidates.
    const double log_half_strand = std::log(0.5);

    // Actual omitted placement mass, so the tolerance is a measurement and not a promise. The help
    // previously said a failure to meet the bound would be reported and nothing reported it.
    double omitted_mass_max = 0.0, omitted_mass_sum = 0.0;
    std::size_t omitted_mass_n = 0;

    std::vector<double> ll(fragments.size() * nh, kNegInf);
    std::vector<double> floors(fragments.size(), kNegInf);
    std::vector<double> band_floors(fragments.size(), kNegInf);
    // Where each fragment landed on each haplotype, or -1 where it did not land at all. The depth
    // channel is built from this: a haplotype carrying sequence the sample does not have shows up as
    // a run of windows with no fragment in them.
    std::vector<std::int32_t> midpoint(fragments.size() * nh, -1);
    // How many of a fragment's two mates carried a usable anchor ON A GIVEN HAPLOTYPE. Per
    // (fragment, haplotype), not per fragment: a mate anchored only to some other shortlisted
    // haplotype is not seeded for the pair being scored, and counting it as seeded attributes its
    // deficit to the wrong stratum. The right unit is still the FRAGMENT -- one anchored mate can
    // rescue the other through the insert constraint -- but the haplotype has to match the dump.
    // THREE states, not one. "Seeded" only means a seed candidate survived for that mate on that
    // haplotype; it does not mean the candidate produced a placement, nor that the two mates formed a
    // valid FR fragment. The evidence can disappear at any of the three, and only splitting them says
    // which.
    //   bit 0/1 : mate 1 / mate 2 had a seed candidate
    //   bit 2/3 : mate 1 / mate 2 produced a successful placement
    //   bit 4   : a valid FR paired state was formed
    std::vector<std::uint8_t> mates_seeded(fragments.size() * nh, 0);

    // ---- zero-seed fallback index -----------------------------------------------------------
    // Pigeonhole: a placement with at most d mismatches cannot mismatch inside all of d+1 DISJOINT
    // pieces of the read, so at least one piece matches exactly and its occurrences propose the
    // start. Lossless only if the pieces are long enough to be worth looking up: the piece length is
    // floor(L / (d+1)) and d comes from max_divergence, so at the default 0.20 a 120 bp read gives
    // 26 pieces of 4 bp, which proposes essentially every position and is no filter at all. Below
    // kMinPiece the fallback scans every start instead. BOTH PATHS RETURN THE SAME PLACEMENTS -- the
    // pigeonhole is an acceleration of the exhaustive scan, never a different answer -- so the choice
    // between them cannot change a score, only the time taken to reach it.
    //
    // Occurrences are NOT capped. A capped index would make the fallback lossy in exactly the
    // repetitive places it exists to serve; repetition is carried as a longer occurrence list.
    constexpr std::size_t kMinPiece = 12;
    std::size_t zs_piece = 0;
    std::unordered_map<std::uint64_t, std::vector<std::pair<std::uint32_t, std::uint32_t>>> zs_index;
    const auto encode_piece = [](const std::string& t, std::size_t at, std::size_t P) -> std::uint64_t {
        std::uint64_t code = 0;
        for (std::size_t i = 0; i < P; ++i) {
            int b;
            switch (t[at + i]) {
                case 'A': case 'a': b = 0; break;
                case 'C': case 'c': b = 1; break;
                case 'G': case 'g': b = 2; break;
                case 'T': case 't': b = 3; break;
                default: return ~0ull;          // ambiguous base: this piece proposes nothing
            }
            code = (code << 2) | static_cast<std::uint64_t>(b);
        }
        return code;
    };
    if (options.zero_seed_fallback) {
        // One piece length for the whole run, the largest that stays lossless for EVERY read.
        std::size_t p_global = 64;
        for (const Fragment& F : fragments) {
            for (const std::string* r : {&F.r1, &F.r2}) {
                if (r->empty()) continue;
                const std::size_t d =
                    static_cast<std::size_t>(options.max_divergence * static_cast<double>(r->size())) + 1;
                p_global = std::min(p_global, r->size() / (d + 1));
            }
        }
        if (p_global >= kMinPiece && !options.zero_seed_exhaustive) {
            zs_piece = std::min<std::size_t>(p_global, 31);
            for (std::size_t hi = 0; hi < nh; ++hi) {
                const std::string& H = haps[hi].seq;
                if (H.size() < zs_piece) continue;
                for (std::size_t at = 0; at + zs_piece <= H.size(); ++at) {
                    const std::uint64_t c = encode_piece(H, at, zs_piece);
                    if (c == ~0ull) continue;
                    zs_index[c].emplace_back(static_cast<std::uint32_t>(hi),
                                             static_cast<std::uint32_t>(at));
                }
            }
        }
    }
    // EVERY distinct placement, not just the best one. A fragment compatible with several copies of a
    // repeat is evidence for a haplotype offering several, and pinning it to one arbitrary copy
    // leaves the others falsely empty -- which then charges the candidate for absence the placement
    // heuristic invented. Deduplicated by rounded midpoint so two syncmer anchors reaching the same
    // biological placement count once. CSR-packed: a dense array would be 48 haplotypes x every
    // fragment x every placement.
    // Orientation is part of the state. Two strands with identical coordinates are DIFFERENT states
    // under the contract and must not be merged when exact states are being retained.
    struct PlacementRec { std::int32_t mid, start, end; bool fwd; double ll; double log_mult = 0.0; };
    std::vector<std::vector<PlacementRec>> placements(
        options.joint_depth ? fragments.size() * nh : 0);

    run_parallel(fragments.size(), options.threads, [&](std::size_t fi) {
        const Fragment& F = fragments[fi];
        const std::size_t total_len = F.bases();
        if (total_len == 0) return;
        floors[fi] = read_ll(
            static_cast<std::size_t>(options.bg_divergence * static_cast<double>(total_len)), total_len);
        // The least penalty consistent with "did not place": exactly at the band edge, computed PER
        // MATE so it reproduces the block-local partial-credit term term-for-term. Never harsher than
        // the background floor, so this can only bound the influence, never raise it.
        {
            double b = 0.0;
            b += read_ll(mate_band_edits(options.max_divergence, F.r1.size()), F.r1.size());
            if (!F.r2.empty()) {
                b += read_ll(mate_band_edits(options.max_divergence, F.r2.size()), F.r2.size());
            }
            band_floors[fi] = std::max(floors[fi], b);
        }

        // Implied read start per (haplotype, orientation), gathered from the read's own syncmers.
        // Clustering on the implied START rather than on the match position is what lets one anchor
        // stand in for the whole read.
        struct Cand { std::uint32_t hap; bool fwd; long start; };
        std::uint64_t found_local = 0, kept_local = 0, anchor_hits_local = 0, combos_local = 0;
        std::uint64_t cart_local = 0, join_local = 0, rescue_pos_local = 0, rescue_placed_local = 0;
        std::uint64_t join_pick_local = 0, cart_pick_local = 0;
        std::uint64_t zs_inv_local = 0, zs_cand_local = 0, zs_ver_local = 0, zs_pl_local = 0;
        std::uint64_t zs_pigeon_local = 0, zs_exh_local = 0;
        const auto gather = [&](const std::string& r, bool fwd,
                                std::vector<Cand>& into) {
            for (const KmerOccurrence& o : collect_syncmers(r, k, s)) {
                const auto it = anchors.find(o.code);
                if (it == anchors.end()) continue;
                for (const auto& [hi, pos] : it->second) {
                    ++anchor_hits_local;
                    into.push_back({hi, fwd, static_cast<long>(pos) - static_cast<long>(o.start)});
                }
            }
        };
        const std::string r1rc = reverse_complement(F.r1);
        const std::string r2rc = F.r2.empty() ? std::string() : reverse_complement(F.r2);

        std::vector<Cand> c1, c2;
        gather(F.r1, true, c1);
        gather(r1rc, false, c1);
        if (!F.r2.empty()) { gather(F.r2, true, c2); gather(r2rc, false, c2); }

        bool trunc_local = false;
        const auto reduce = [&](std::vector<Cand>& c) {
            // Bucket implied starts to 64 bp and keep, per haplotype, the `placement_topk`
            // commonest. NOTE, so no downstream comment overstates it: this is NOT every plausible
            // placement. Anchors occurring more than `max_anchor_occ` times per haplotype were
            // already dropped, and of the clusters that remain only the top few survive here, so a
            // fragment inside a repeat with many copies is represented by a few of them, not all.
            // Whether that approximation matters is an open question, not a settled one.
            std::map<std::pair<std::uint32_t, long>, std::pair<int, long>> tally;
            for (const Cand& x : c) {
                auto& e = tally[{x.hap, (x.fwd ? 1 : -1) *
                                 (x.start / static_cast<long>(std::max<std::size_t>(1, options.placement_bin)) + 1)}];
                ++e.first;
                e.second = x.start;
            }
            std::map<std::uint32_t, std::vector<std::pair<int, std::pair<bool, long>>>> per_hap;
            for (const auto& [key, val] : tally) {
                per_hap[key.first].push_back({val.first, {key.second > 0, val.second}});
            }
            std::vector<Cand> keep;
            bool truncated = false;
            for (auto& [hi, v] : per_hap) {
                std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
                found_local += v.size();
                if (v.size() > options.placement_topk) truncated = true;
                for (std::size_t i = 0; i < std::min(options.placement_topk, v.size()); ++i) {
                    keep.push_back({hi, v[i].second.first, v[i].second.second});
                    ++kept_local;
                }
            }
            if (truncated) trunc_local = true;
            c.swap(keep);
        };
        reduce(c1);
        reduce(c2);
        for (const Cand& c : c1) mates_seeded[fi * nh + c.hap] |= 1u;
        for (const Cand& c : c2) mates_seeded[fi * nh + c.hap] |= 2u;
        {
            static std::mutex comp_mu;
            std::lock_guard<std::mutex> lk(comp_mu);
            out.completeness.clusters_found += found_local;
            out.completeness.clusters_kept += kept_local;
            out.completeness.anchor_hits += anchor_hits_local;
            if (trunc_local) ++out.completeness.fragments_truncated;
        }

        const std::size_t band1 =
            static_cast<std::size_t>(options.max_divergence * static_cast<double>(F.r1.size())) + 1;
        const std::size_t band2 = F.r2.empty() ? 1 :
            static_cast<std::size_t>(options.max_divergence * static_cast<double>(F.r2.size())) + 1;

        // `fwd` records which orientation the query was in, so only valid FR combinations are
        // paired. Without it, two mates in the SAME orientation could be combined into a "fragment".
        struct Placed { bool ok = false; std::size_t edits = 0; long start = 0; long end = 0; bool fwd = true; };
        const auto place_all = [&](const std::vector<Cand>& cands, std::uint32_t hi,
                                   const std::string& fwd_seq, const std::string& rev_seq,
                                   std::size_t band) {
            std::vector<Placed> found;
            for (const Cand& c : cands) {
                if (c.hap != hi) continue;
                const std::string& q = c.fwd ? fwd_seq : rev_seq;
                if (q.empty()) continue;
                const long lo = std::max<long>(0, c.start - static_cast<long>(options.anchor_slack));
                const long hi_end = std::min<long>(static_cast<long>(haps[hi].seq.size()),
                                                   c.start + static_cast<long>(q.size() + options.anchor_slack));
                if (hi_end - lo < static_cast<long>(q.size()) / 2) continue;
                Placed p;
                if (options.hamming_emission) {
                    // Fixed position, no sliding: exactly what the reference scores.
                    if (c.start < 0 || c.start + static_cast<long>(q.size()) >
                                       static_cast<long>(haps[hi].seq.size())) continue;
                    std::size_t mism = 0;
                    for (std::size_t bi2 = 0; bi2 < q.size(); ++bi2) {
                        if (q[bi2] != haps[hi].seq[static_cast<std::size_t>(c.start) + bi2]) ++mism;
                    }
                    if (mism > band) continue;
                    p.ok = true;
                    p.edits = mism;
                    p.start = c.start;
                    p.end = c.start + static_cast<long>(q.size()) - 1;
                    p.fwd = c.fwd;
                    found.push_back(p);
                    continue;
                }
                const std::string window = haps[hi].seq.substr(static_cast<std::size_t>(lo),
                                                               static_cast<std::size_t>(hi_end - lo));
                const ReadFit f = infix_align(q, window, band);
                if (!f.ok) continue;
                p.ok = true;
                p.edits = f.edits;
                p.start = lo + static_cast<long>(f.start);
                p.end = lo + static_cast<long>(f.end);
                p.fwd = c.fwd;
                found.push_back(p);
            }
            return found;
        };

        std::vector<std::uint32_t> touched;
        if (options.zero_seed_fallback) {
            // EVERY candidate haplotype, not only those an anchor pointed at. The loop below is the
            // only place the fallback can run, and a fragment whose mates seed nowhere has an empty
            // anchor set -- so restricting the loop to anchored haplotypes skips precisely the
            // fragments the fallback exists for. Measured before this fix: the fallback fired 1-11
            // times per cell against an unseeded stratum of 2-8 fragments over 6 haplotypes, found
            // at most 3 placements, and recovered exactly 0.00 nats of deficit.
            touched.resize(nh);
            for (std::size_t hx = 0; hx < nh; ++hx) touched[hx] = static_cast<std::uint32_t>(hx);
        } else {
            for (const Cand& c : c1) touched.push_back(c.hap);
            for (const Cand& c : c2) touched.push_back(c.hap);
            std::sort(touched.begin(), touched.end());
            touched.erase(std::unique(touched.begin(), touched.end()), touched.end());
        }

        for (const std::uint32_t hi : touched) {
            // Both mates are placed JOINTLY, over every combination of their candidate anchors,
            // rather than each taking its own best position. Independently, inside a duplication the
            // two mates settle on different copies and the pair then looks 13 kb long -- the mates
            // are one observation and have to be placed as one.
            std::vector<Placed> a1 = place_all(c1, hi, F.r1, r1rc, band1);
            std::vector<Placed> a2 = F.r2.empty() ? std::vector<Placed>{}
                                                  : place_all(c2, hi, F.r2, r2rc, band2);

            // ---- zero-seed fallback ------------------------------------------------------------
            // Runs ONLY when neither mate has a primary placement on this haplotype. A fragment that
            // recruitment placed is never touched by it, so it cannot perturb an existing
            // contribution -- that is a structural property of this guard, not a tuning choice.
            const bool zs_run = options.zero_seed_fallback
                && (options.zero_seed_complete || (a1.empty() && a2.empty()));
            if (zs_run) {
                ++zs_inv_local;
                const std::string& H = haps[hi].seq;
                const auto fallback_place = [&](const std::string& fwd, const std::string& rev,
                                                std::size_t band) {
                    std::vector<Placed> got;
                    for (int o = 0; o < 2; ++o) {
                        const std::string& q = (o == 0) ? fwd : rev;
                        if (q.empty() || H.size() < q.size()) continue;
                        const long last = static_cast<long>(H.size() - q.size());
                        std::vector<long> starts;
                        const std::size_t d =
                            static_cast<std::size_t>(options.max_divergence
                                                     * static_cast<double>(q.size())) + 1;
                        // Every one of the d+1 pieces must be encodable, or the pigeonhole
                        // guarantee is void. Skipping an unencodable piece leaves only d pieces, and
                        // d pieces cannot corner d mismatches. The case is real rather than
                        // hypothetical: where the read and the haplotype BOTH carry N at a position,
                        // this Hamming model scores it a match ('N' == 'N') while encode_piece
                        // rejects both sides, so the index would never propose the start. A
                        // haplotype N against a read base needs no special handling -- that already
                        // counts as a mismatch, so such a piece is not mismatch-free anyway.
                        bool all_encodable = zs_piece > 0 && q.size() >= (d + 1) * zs_piece;
                        if (all_encodable) {
                            for (std::size_t pc = 0; pc <= d; ++pc) {
                                if (encode_piece(q, pc * zs_piece, zs_piece) == ~0ull) {
                                    all_encodable = false; break;
                                }
                            }
                        }
                        if (all_encodable) {
                            ++zs_pigeon_local;
                            // d+1 disjoint pieces: at most d of them can carry a mismatch.
                            for (std::size_t pc = 0; pc <= d; ++pc) {
                                const std::size_t at = pc * zs_piece;
                                const std::uint64_t code = encode_piece(q, at, zs_piece);
                                const auto it = zs_index.find(code);
                                if (it == zs_index.end()) continue;
                                for (const auto& hp : it->second) {
                                    if (hp.first != hi) continue;
                                    const long st = static_cast<long>(hp.second)
                                                  - static_cast<long>(at);
                                    if (st < 0 || st > last) continue;
                                    starts.push_back(st);
                                    ++zs_cand_local;
                                }
                            }
                            std::sort(starts.begin(), starts.end());
                            starts.erase(std::unique(starts.begin(), starts.end()), starts.end());
                        } else {
                            ++zs_exh_local;
                            starts.reserve(static_cast<std::size_t>(last) + 1);
                            for (long st = 0; st <= last; ++st) starts.push_back(st);
                            zs_cand_local += static_cast<std::uint64_t>(last) + 1;
                        }
                        for (const long st : starts) {
                            ++zs_ver_local;
                            std::size_t mism = 0;
                            for (std::size_t bi2 = 0; bi2 < q.size() && mism <= band; ++bi2) {
                                if (q[bi2] != H[static_cast<std::size_t>(st) + bi2]) ++mism;
                            }
                            if (mism > band) continue;
                            Placed r;
                            r.ok = true; r.edits = mism; r.start = st;
                            r.end = st + static_cast<long>(q.size()) - 1;
                            r.fwd = (o == 0);
                            ++zs_pl_local;
                            got.push_back(r);
                        }
                    }
                    return got;
                };
                std::vector<Placed> f1 = fallback_place(F.r1, r1rc, band1);
                std::vector<Placed> f2 = F.r2.empty() ? std::vector<Placed>{}
                                                      : fallback_place(F.r2, r2rc, band2);
                if (options.zero_seed_complete) {
                    // Union, deduplicated on (start, strand). Under the fixed-position Hamming
                    // contract that pair identifies the state -- the end follows from the read
                    // length -- so a placement recruitment already found is not added twice.
                    const auto merge = [](std::vector<Placed>& into, const std::vector<Placed>& add) {
                        for (const Placed& p : add) {
                            bool seen = false;
                            for (const Placed& q : into) {
                                if (q.start == p.start && q.fwd == p.fwd) { seen = true; break; }
                            }
                            if (!seen) into.push_back(p);
                        }
                    };
                    merge(a1, f1);
                    merge(a2, f2);
                } else {
                    a1 = std::move(f1);
                    a2 = std::move(f2);
                }
                // Whatever the fallback produced now flows through the SAME machinery as any other
                // placement: one mate only goes to interval rescue below, both go to the adaptive
                // coordinate join. No second scoring path exists for it.
                if (!a1.empty() || !a2.empty()) mates_seeded[fi * nh + hi] |= 64u;
            }
            if (!a1.empty()) mates_seeded[fi * nh + hi] |= 4u;
            if (!a2.empty()) mates_seeded[fi * nh + hi] |= 8u;

            // ---- one-mate rescue ---------------------------------------------------------------
            // Exactly one mate placed: the other is searched in the interval FR orientation and the
            // insert prior allow, on this same haplotype. No truth is used -- the interval follows
            // from the anchored placement, the strand and the prior's support.
            std::vector<Placed> resc1, resc2;
            if (options.mate_rescue && !F.r2.empty() && (a1.empty() != a2.empty())) {
                const bool have1 = !a1.empty();
                const std::vector<Placed>& have = have1 ? a1 : a2;
                const std::string& want_fwd = have1 ? F.r2 : F.r1;
                const std::string& want_rev = have1 ? r2rc : r1rc;
                const std::size_t want_band = have1 ? band2 : band1;
                std::vector<Placed>& into = have1 ? resc2 : resc1;
                const long n = static_cast<long>(haps[hi].seq.size());
                // Grouped by the ANCHORED mate's strand, because that fixes both the partner's
                // orientation and which side of the anchor the interval lies on. Within a group the
                // intervals are merged, so each position is examined ONCE however many anchored
                // placements propose it.
                //
                // Merging also makes the rescued set duplicate-free BY CONSTRUCTION: merged
                // intervals within a group are disjoint, and the two groups differ in strand, so no
                // (start, strand) state can be proposed twice. Scanning per anchored placement could
                // push the same state once per proposing placement.
                //
                // That is stated as a property of the construction, NOT as a repair of an observed
                // defect. Removing those duplicates changed no score in any configuration measured
                // -- joint-marginal and plain, max and marginalised -- so whether they ever inflated
                // a likelihood is UNDEMONSTRATED, and the measured effect of merging here is a 9%
                // reduction in positions examined and nothing else. Rescue fires only when exactly
                // one mate placed, so the anchored list is short and its intervals largely do not
                // overlap; there is little to merge.
                for (int g = 0; g < 2; ++g) {
                    const bool anchored_fwd = (g == 0);
                    const std::string& q = anchored_fwd ? want_rev : want_fwd;
                    if (q.empty()) continue;
                    std::vector<std::pair<long, long>> iv;
                    for (const Placed& p : have) {
                        if (p.fwd != anchored_fwd) continue;
                        // FR: an anchored forward mate has its partner reverse and DOWNSTREAM; an
                        // anchored reverse mate has it forward and upstream.
                        long lo_s, hi_s;
                        if (anchored_fwd) {
                            lo_s = p.start + ins_prior.lo - static_cast<long>(q.size());
                            hi_s = p.start + ins_prior.hi - static_cast<long>(q.size());
                        } else {
                            lo_s = p.end - ins_prior.hi + 1;
                            hi_s = p.end - ins_prior.lo + 1;
                        }
                        lo_s = std::max<long>(0, lo_s);
                        hi_s = std::min<long>(hi_s, n - static_cast<long>(q.size()));
                        if (hi_s >= lo_s) iv.emplace_back(lo_s, hi_s);
                    }
                    if (iv.empty()) continue;
                    std::sort(iv.begin(), iv.end());
                    std::vector<std::pair<long, long>> merged;
                    for (const auto& e : iv) {
                        if (!merged.empty() && e.first <= merged.back().second + 1) {
                            merged.back().second = std::max(merged.back().second, e.second);
                        } else {
                            merged.push_back(e);
                        }
                    }
                    for (const auto& e : merged) {
                        rescue_pos_local += static_cast<std::uint64_t>(e.second - e.first + 1);
                        for (long st = e.first; st <= e.second; ++st) {
                            std::size_t mism = 0;
                            for (std::size_t bi2 = 0; bi2 < q.size() && mism <= want_band; ++bi2) {
                                if (q[bi2] != haps[hi].seq[static_cast<std::size_t>(st) + bi2]) ++mism;
                            }
                            if (mism > want_band) continue;
                            Placed r;
                            r.ok = true; r.edits = mism; r.start = st;
                            r.end = st + static_cast<long>(q.size()) - 1;
                            r.fwd = !anchored_fwd;
                            ++rescue_placed_local;
                            into.push_back(r);
                        }
                    }
                }
                if (!into.empty()) mates_seeded[fi * nh + hi] |= 32u;   // rescued
            }
            const std::vector<Placed>& b1 = resc1.empty() ? a1 : resc1;
            const std::vector<Placed>& b2 = resc2.empty() ? a2 : resc2;
            // The hypothetical product, counted whichever path runs, so the two costs are always
            // comparable rather than each mode reporting only its own.
            if (!b1.empty() && !b2.empty()) {
                cart_local += static_cast<std::uint64_t>(b1.size()) * b2.size();
            }

            // ---- coordinate join -----------------------------------------------------------
            // Forward placements keyed by start, reverse placements keyed by end; a fragment state
            // is (forward start s, reverse end e) with e - s + 1 inside the insert prior's support.
            //
            // COST, stated honestly: O(|F| + |R| + K), where K is the number of coordinate pairs
            // that fall within the insert support. It is NOT bounded by the library's insert width.
            // K is itself quadratic whenever many coordinates cluster inside one allowed interval --
            // a short tandem repeat whose whole array fits inside one insert is exactly that case.
            // On the array fixtures measured here K grew linearly in copy number, so the join saved
            // one factor of copies, but that is a measurement on those fixtures and not a bound.
            //
            // It is exact rather than a pruning: ins_prior.log_at returns -inf outside [lo, hi], so
            // every combination the product forms and the join skips carries exactly zero mass.
            //
            // ADAPTIVE DISPATCH. Both paths compute the same sum -- measured identical to 0.0 on
            // every diplotype -- so which one runs is a cost decision and not a model parameter, and
            // it can be made per fragment-haplotype from the actual coordinates rather than assumed
            // from the copy number. The product wins on small inputs, where building two ordered
            // maps costs more than the |b1| x |b2| combinations it avoids, and on clustered
            // coordinates, where K approaches the product anyway. K is counted by the same
            // two-pointer sweep run in counting mode, which is O(|F| + |R|) -- window sizes come
            // from pointer arithmetic, so counting the pairs does not require visiting them.
            const bool join_available = options.coordinate_join && !b1.empty() && !b2.empty();
            bool do_join = false;
            CoordAgg jf1, jr1, jf2, jr2;
            if (join_available) {
                for (const Placed& x : b1) {
                    (x.fwd ? jf1 : jr1).add(x.fwd ? x.start : x.end, read_ll(x.edits, F.r1.size()));
                }
                for (const Placed& y : b2) {
                    (y.fwd ? jf2 : jr2).add(y.fwd ? y.start : y.end, read_ll(y.edits, F.r2.size()));
                }
                // The valid-FR bit means "an FR pair exists with the reverse mate downstream",
                // which is what the product tests. It is deliberately BROADER than "within the
                // insert prior's support" -- keeping the same condition here keeps the stage
                // instrumentation comparable across the two paths.
                const auto any_fr = [](const CoordAgg& fwd, const CoordAgg& rev) {
                    return !fwd.at.empty() && !rev.at.empty() &&
                           rev.at.rbegin()->first >= fwd.at.begin()->first;
                };
                if (any_fr(jf1, jr2) || any_fr(jf2, jr1)) mates_seeded[fi * nh + hi] |= 16u;

                const auto count_pairs = [&](const CoordAgg& FA, const CoordAgg& RA) -> std::uint64_t {
                    if (FA.at.empty() || RA.at.empty()) return 0;
                    std::vector<long> RC;
                    RC.reserve(RA.at.size());
                    for (const auto& e : RA.at) RC.push_back(e.first);
                    std::size_t wlo = 0, whi = 0;
                    std::uint64_t k = 0;
                    for (const auto& fe : FA.at) {
                        const long elo = fe.first + ins_prior.lo - 1;
                        const long ehi = fe.first + ins_prior.hi - 1;
                        while (wlo < RC.size() && RC[wlo] < elo) ++wlo;
                        if (whi < wlo) whi = wlo;
                        while (whi < RC.size() && RC[whi] <= ehi) ++whi;
                        k += static_cast<std::uint64_t>(whi - wlo);
                    }
                    return k;
                };
                const std::uint64_t cart_cost =
                    static_cast<std::uint64_t>(b1.size()) * b2.size();
                const std::uint64_t k_pairs = count_pairs(jf1, jr2) + count_pairs(jf2, jr1);
                const std::uint64_t join_cost = k_pairs + jf1.at.size() + jr1.at.size()
                                              + jf2.at.size() + jr2.at.size();
                do_join = options.force_join || join_cost < cart_cost;
                if (do_join) ++join_pick_local; else ++cart_pick_local;
            }
            const auto for_each_join = [&](const std::function<void(double, double, long, long,
                                                                    bool)>& cb) {
                const auto sweep = [&](const CoordAgg& FA, const CoordAgg& RA, bool mate1_fwd) {
                    if (FA.at.empty() || RA.at.empty()) return;
                    // TWO-POINTER over sorted coordinates, NOT a probe of every allowed insert
                    // length. The probe form -- for each start, look up start + L - 1 for every L in
                    // the prior's support -- is what the join was first written as, and measured it
                    // was 12-29x MORE work than the product it replaced: the support is mean +/- 4sd,
                    // about 400 lengths, so it pays only above ~400 placements per mate, which is far
                    // beyond any copy number here. Advancing a window over the reverse coordinates
                    // that actually EXIST costs one visit per in-support combination instead. That is
                    // never more than the product and usually far less, but it is not a smaller
                    // COMPLEXITY class: if every reverse coordinate sits inside every forward
                    // window, K = |F| x |R| and the join degenerates to the product plus sorting.
                    std::vector<std::pair<long, std::pair<double, double>>> R(RA.at.begin(),
                                                                              RA.at.end());
                    std::size_t wlo = 0, whi = 0;   // [wlo, whi) = reverse ends inside the prior
                    // FA.at is a std::map, so starts arrive in increasing order and both window
                    // edges only ever move right: total pointer movement is O(|R|) for the sweep,
                    // not per start.
                    for (const auto& fe : FA.at) {
                        const long s2 = fe.first;
                        const long elo = s2 + ins_prior.lo - 1, ehi = s2 + ins_prior.hi - 1;
                        while (wlo < R.size() && R[wlo].first < elo) ++wlo;
                        if (whi < wlo) whi = wlo;
                        while (whi < R.size() && R[whi].first <= ehi) ++whi;
                        for (std::size_t t = wlo; t < whi; ++t) {
                            ++join_local;
                            const long e2 = R[t].first;
                            const double ins = options.use_insert_size
                                ? ins_prior.log_at(e2 - s2 + 1) : 0.0;
                            cb(fe.second.first + R[t].second.first + log_half_strand + ins,
                               fe.second.second + R[t].second.second + log_half_strand + ins,
                               s2, e2, mate1_fwd);
                        }
                    }
                };
                sweep(jf1, jr2, true);    // mate 1 forward, mate 2 reverse
                sweep(jf2, jr1, false);   // mate 2 forward, mate 1 reverse
            };
            double lp = kNegInf;      // the accumulated likelihood: max, or the sum when marginalising
            double best = kNegInf;    // always the best single placement, so the midpoint is a real one
            long mid = -1;
            const double miss1 = read_ll(band1, F.r1.size());
            const double miss2 = F.r2.empty() ? 0.0 : read_ll(band2, F.r2.size());
            const auto consider = [&](double v, long m) {
                if (v > best) { best = v; mid = m; }
                lp = options.marginalise_placements ? log_add(lp, v) : std::max(lp, v);
            };
            if (do_join) {
                for_each_join([&](double mass_v, double best_v, long s, long e, bool) {
                    // The summed mass and the best single placement are tracked separately: with
                    // several placements aggregated at one coordinate they are no longer the same
                    // number, and the midpoint must remain a real placement's.
                    if (best_v > best) { best = best_v; mid = (s + e) / 2; }
                    lp = options.marginalise_placements ? log_add(lp, mass_v)
                                                        : std::max(lp, best_v);
                });
            } else if (!b1.empty() && !b2.empty()) {
                // Both mates placed: the pair's placements are the COMBINATIONS, and the single-mate
                // terms below would double count them, so they are skipped.
                for (const Placed& x : b1) {
                    for (const Placed& y : b2) {
                        // FR only: the two mates of a fragment face each other. Combining same-strand
                        // placements invents fragments the library cannot produce.
                        if (x.fwd == y.fwd) continue;
                        const Placed& fw = x.fwd ? x : y;
                        const Placed& rv = x.fwd ? y : x;
                        // SHARED with the bounded search's rule, so the two cannot drift.
                        if (!fr_reverse_downstream(fw.start, rv.end)) continue;
                        mates_seeded[fi * nh + hi] |= 16u;
                        double v = read_ll(x.edits, F.r1.size()) + read_ll(y.edits, F.r2.size())
                                 + log_half_strand;
                        if (options.use_insert_size) {
                            v += insert_ll(static_cast<double>(rv.end - fw.start + 1), ins_prior);
                        }
                        consider(v, (fw.start + rv.end) / 2);
                    }
                }
            } else if (F.r2.empty()) {
            // SINGLE-END ONLY. For a PAIRED fragment, one placed mate plus a synthetic partner
            // pinned at the band boundary is a state the exhaustive reference does not have: the
            // reference scores complete (start, insert length, orientation) fragments, so a pair no
            // complete placement explains falls to the background mixture. Granting it partial
            // credit here made the accelerated scorer score HIGHER than the reference on exactly the
            // fragments with mates_placed=1 and valid_fr=0, uniformly, and only on diplotypes where
            // alternative repeat copies offer one-mate-only placements -- which is the non-truth
            // side of the panel. Measured: all 22 discrepant fragments carried that signature.
            // With no valid FR placement the mass stays -inf and falls through to floors[fi] below,
            // which IS the background mixture.
                for (const Placed& x : b1) consider(read_ll(x.edits, F.r1.size()) + miss2,
                                                    (x.start + x.end) / 2);
                for (const Placed& y : b2) consider(miss1 + read_ll(y.edits, F.r2.size()),
                                                    (y.start + y.end) / 2);
            }
            if (lp == kNegInf) {
                ll[fi * nh + hi] = options.band_floor ? band_floors[fi] : floors[fi];
                continue;
            }
            midpoint[fi * nh + hi] = static_cast<std::int32_t>(mid);
            ll[fi * nh + hi] = lp;

            if (options.joint_depth) {
                // Keep each distinct placement with its own likelihood. Deduplicate on the rounded
                // midpoint: two anchors that put the fragment in the same place are one placement,
                // and counting them twice would let anchor density masquerade as depth.
                auto& into = placements[fi * nh + hi];
                // With placement_dedup > 0 the old behaviour: collapse placements sharing a midpoint
                // within that radius, keeping the best. With 0, every distinct (start, end) survives,
                // which is what the exact reference enumerates -- two different mate pairings across
                // repeat copies can share a midpoint and are NOT the same state.
                const auto add = [&](long m, long st, long en, bool fwd, double v) {
                    const std::int32_t key = static_cast<std::int32_t>(m);
                    for (auto& e : into) {
                        if (options.placement_dedup > 0) {
                            if (std::abs(e.mid - key) <= static_cast<long>(options.placement_dedup)) {
                                e.ll = std::max(e.ll, v); return;
                            }
                        } else if (e.start == static_cast<std::int32_t>(st) &&
                                   e.end == static_cast<std::int32_t>(en) && e.fwd == fwd) {
                            e.ll = std::max(e.ll, v); return;
                        }
                    }
                    into.push_back({key, static_cast<std::int32_t>(st),
                                    static_cast<std::int32_t>(en), fwd, v});
                };
                if (do_join) {
                    // Placements keep the BEST value per state, exactly as the product's `add` does
                    // on collision, so the two enumerate the same states with the same values.
                    for_each_join([&](double, double best_v, long s, long e, bool m1f) {
                        add((s + e) / 2, s, e, m1f, best_v);
                    });
                } else if (!b1.empty() && !b2.empty()) {
                    for (const Placed& x : b1) {
                        for (const Placed& y : b2) {
                            ++combos_local;
                            if (x.fwd == y.fwd) continue;
                            const Placed& fw = x.fwd ? x : y;
                            const Placed& rv = x.fwd ? y : x;
                            if (!fr_reverse_downstream(fw.start, rv.end)) continue;
                            double v = read_ll(x.edits, F.r1.size()) + read_ll(y.edits, F.r2.size())
                                     + log_half_strand;
                            if (options.use_insert_size) {
                                v += insert_ll(static_cast<double>(rv.end - fw.start + 1), ins_prior);
                            }
                            add((fw.start + rv.end) / 2, fw.start, rv.end, x.fwd, v);
                        }
                    }
                } else if (F.r2.empty()) {
                    // Same restriction as the likelihood accumulation above: the partial-credit
                    // state exists for genuinely single-end fragments only.
                    for (const Placed& x : b1) add((x.start + x.end) / 2, x.start, x.end, x.fwd,
                                                   read_ll(x.edits, F.r1.size()) + miss2);
                    for (const Placed& y : b2) add((y.start + y.end) / 2, y.start, y.end, y.fwd,
                                                   miss1 + read_ll(y.edits, F.r2.size()));
                }
            }
        }
        {
            static std::mutex cm;
            std::lock_guard<std::mutex> lk(cm);
            out.completeness.mate_combinations += combos_local;
            out.completeness.cartesian_combinations += cart_local;
            out.completeness.join_operations += join_local;
            out.completeness.rescue_positions += rescue_pos_local;
            out.completeness.rescue_placements += rescue_placed_local;
            out.completeness.join_chosen += join_pick_local;
            out.completeness.cartesian_chosen += cart_pick_local;
            out.completeness.zs_invocations += zs_inv_local;
            out.completeness.zs_candidate_starts += zs_cand_local;
            out.completeness.zs_verified_starts += zs_ver_local;
            out.completeness.zs_placements += zs_pl_local;
            out.completeness.zs_pigeonhole += zs_pigeon_local;
            out.completeness.zs_exhaustive += zs_exh_local;
            if (options.joint_depth) {
                for (std::size_t hi = 0; hi < nh; ++hi) {
                    out.completeness.placements_before_grouping += placements[fi * nh + hi].size();
                }
            }
        }
        if (options.multiplicity_aware) {
            for (std::size_t hi = 0; hi < nh; ++hi) {
                auto& v = placements[fi * nh + hi];
                if (v.size() < 2) continue;
                // Group placements of EQUAL likelihood -- which is what every copy of a perfect
                // repeat produces -- into one carrying log P + log(multiplicity). Exact, and it turns
                // an N-copy array from N placements into one.
                std::sort(v.begin(), v.end(),
                          [](const PlacementRec& a, const PlacementRec& b) { return a.ll > b.ll; });
                std::vector<PlacementRec> grouped;
                for (const PlacementRec& p : v) {
                    if (!grouped.empty() && std::abs(grouped.back().ll - p.ll) < 1e-9) {
                        grouped.back().log_mult = log_add(grouped.back().log_mult, 0.0);
                        continue;
                    }
                    grouped.push_back(p);
                }
                // Prune by GROUP MASS, not by representative likelihood. A lower-likelihood group
                // with high multiplicity can carry more probability than a singleton above it, so
                // sorting on `ll` alone drops the wrong tail -- which is the same confusion between
                // placement COUNT and placement MASS this whole feature exists to remove.
                std::sort(grouped.begin(), grouped.end(),
                          [](const PlacementRec& a, const PlacementRec& b) {
                              return (a.ll + a.log_mult) > (b.ll + b.log_mult);
                          });
                // Prune by omitted MASS, not by count: drop the tail only while what it carries stays
                // under the tolerance.
                double total = kNegInf;
                for (const PlacementRec& p : grouped) total = log_add(total, p.ll + p.log_mult);
                if (total != kNegInf) {
                    const double keep_floor = total + std::log(options.mass_tolerance);
                    double dropped = kNegInf;
                    std::size_t keep = grouped.size();
                    while (keep > 1) {
                        const PlacementRec& last = grouped[keep - 1];
                        const double next = log_add(dropped, last.ll + last.log_mult);
                        if (next > keep_floor) break;
                        dropped = next;
                        --keep;
                    }
                    grouped.resize(keep);
                    // Every non-empty fragment-haplotype counts, contributing zero when nothing was
                    // pruned. Counting only the cases that dropped something makes the mean the mean
                    // OF THE DROPS, which is a different and much larger quantity than the mean
                    // omission per placement set -- and it is the latter the tolerance is about.
                    const double frac = (dropped == kNegInf) ? 0.0 : std::exp(dropped - total);
                    static std::mutex om;
                    std::lock_guard<std::mutex> lk(om);
                    omitted_mass_max = std::max(omitted_mass_max, frac);
                    omitted_mass_sum += frac;
                    ++omitted_mass_n;
                }
                v.swap(grouped);
            }
        }
        for (std::size_t hi = 0; hi < nh; ++hi) {
            if (ll[fi * nh + hi] == kNegInf) {
                ll[fi * nh + hi] = options.band_floor ? band_floors[fi] : floors[fi];
            }
        }
    });

    // State (3): drop fragments whose search never happened on some candidate. floors[fi] == kNegInf
    // is the existing "skip this fragment" convention in every scoring loop, so this reuses it rather
    // than adding a parallel flag that a loop could forget to check.
    std::size_t neutralised = 0;
    if (options.unplaced_neutral) {
        for (std::size_t fi = 0; fi < fragments.size(); ++fi) {
            if (floors[fi] == kNegInf) continue;
            for (std::size_t hi = 0; hi < nh; ++hi) {
                if (mates_seeded[fi * nh + hi] == 0) { floors[fi] = kNegInf; ++neutralised; break; }
            }
        }
    }
    out.n_neutralised = neutralised;

    // True-origin placement recall. wgsim names each fragment "<haplotype>_<start>_<end>_...", so when
    // the origin haplotype is in the shortlist the placement the recruiter SHOULD have found is known
    // exactly. Reported because raising placement_topk restores true placements and adds false ones at
    // the same time, and a retention percentage cannot distinguish those.
    if (options.joint_depth) {
        std::unordered_map<std::string, std::size_t> by_name;
        for (std::size_t h = 0; h < nh; ++h) by_name.emplace(out.shortlist[h], h);
        for (std::size_t fi = 0; fi < fragments.size(); ++fi) {
            const std::string& nm = fragments[fi].name;
            // trailing fields are _start_end_e1_e2_idx; the haplotype name carries no underscore
            const std::size_t p1 = nm.find('_');
            if (p1 == std::string::npos) continue;
            const std::size_t p2 = nm.find('_', p1 + 1);
            const std::size_t p3 = nm.find('_', p2 == std::string::npos ? p1 + 1 : p2 + 1);
            if (p2 == std::string::npos || p3 == std::string::npos) continue;
            const auto it = by_name.find(nm.substr(0, p1));
            if (it == by_name.end()) continue;
            long ts = 0, te = 0;
            try {
                ts = std::stol(nm.substr(p1 + 1, p2 - p1 - 1));
                te = std::stol(nm.substr(p2 + 1, p3 - p2 - 1));
            } catch (...) { continue; }
            const long truth_mid = (ts + te) / 2;
            ++out.completeness.truth_resolvable;
            // EITHER ORIENTATION. A haplotype's block-concatenated spelling can be the reverse
            // complement of the path spelling the reads were simulated from -- measured, one of
            // HG04036's two haplotypes is, and checking only the forward coordinate reported 48%
            // recall where the true figure is 100%. The alignment already tries both strands, so
            // this is a property of the coordinate frame and not of the placement.
            const long L = static_cast<long>(haps[it->second].seq.size());
            bool found = false;
            std::uint64_t extra = 0;
            for (const auto& pr : placements[fi * nh + it->second]) {
                const long m = static_cast<long>(pr.mid);
                if (std::abs(m - truth_mid) <= 150 || std::abs((L - m) - truth_mid) <= 150) found = true;
                else ++extra;
            }
            if (found) ++out.completeness.truth_recovered;
            out.completeness.spurious_placements += extra;
        }
        out.completeness.recall_measured = out.completeness.truth_resolvable > 0;
    }
    out.completeness.fragments_total = fragments.size();
    out.completeness.omitted_mass_max = omitted_mass_max;
    out.completeness.omitted_mass_mean =
        omitted_mass_n == 0 ? 0.0 : omitted_mass_sum / static_cast<double>(omitted_mass_n);
    if (options.joint_depth) {
        for (std::size_t i = 0; i < placements.size(); ++i) {
            out.completeness.groups_after_grouping += placements[i].size();
        }
    }

    for (std::size_t fi = 0; fi < fragments.size(); ++fi) {
        double lo = kNegInf, hi = kNegInf;
        for (std::size_t h = 0; h < nh; ++h) {
            const double v = ll[fi * nh + h];
            if (hi == kNegInf || v > hi) hi = v;
            if (lo == kNegInf || v < lo) lo = v;
        }
        if (hi - lo > 1e-9) ++out.n_informative;
    }

    for (std::size_t h = 0; h < nh; ++h) {
        double solo = 0.0;
        std::size_t placed = 0;
        for (std::size_t fi = 0; fi < fragments.size(); ++fi) {
            if (floors[fi] == kNegInf) continue;
            solo += ll[fi * nh + h];
            if (midpoint[fi * nh + h] >= 0) ++placed;
        }
        out.haplotypes[h].solo_ll = solo;
        out.haplotypes[h].placed = placed;
    }

    // ---- depth channel: per-haplotype window coverage ----------------------------------------
    std::vector<double> cov_ll(nh, 0.0);
    if (options.coverage_weight > 0.0 && options.coverage_window > 0) {
        for (std::size_t h = 0; h < nh; ++h) {
            const std::size_t nw = haps[h].seq.size() / options.coverage_window + 1;
            std::vector<double> count(nw, 0.0);
            for (std::size_t fi = 0; fi < fragments.size(); ++fi) {
                const std::int32_t m = midpoint[fi * nh + h];
                if (m < 0) continue;
                const std::size_t w = static_cast<std::size_t>(m) / options.coverage_window;
                if (w < nw) count[w] += 1.0;
            }
            // Rate fitted from the windows that ARE covered, so a homologue carrying half the depth
            // is not penalised for being one of two -- only for windows that should have reads and
            // have none.
            out.haplotypes[h].windows = nw;
            for (const double c : count) if (c == 0.0) ++out.haplotypes[h].zero_windows;
            std::vector<double> nz;
            for (const double c : count) if (c > 0.0) nz.push_back(c);
            if (nz.size() < 4) { cov_ll[h] = 0.0; continue; }
            std::nth_element(nz.begin(), nz.begin() + nz.size() / 2, nz.end());
            const double rate = std::max(1e-6, nz[nz.size() / 2]);
            double total = 0.0;
            for (const double c : count) {
                total += c * std::log(rate) - rate - std::lgamma(c + 1.0);
            }
            cov_ll[h] = options.coverage_weight * total;
            out.haplotypes[h].coverage_ll = cov_ll[h];
        }
    }

    // ---- diploid pair likelihood over whole haplotypes ---------------------------------------
    std::vector<double> log_len(nh, 0.0);
    for (std::size_t h = 0; h < nh; ++h) {
        const double occupied = options.use_insert_size ? options.fragment_len : 150.0;
        log_len[h] = std::log(std::max(1.0, static_cast<double>(haps[h].seq.size()) - occupied + 1.0));
    }
    const double log_mix = std::log1p(-options.outlier_mix);
    const double log_out = std::log(options.outlier_mix);

    // lambda, fitted once and outside every candidate: the observed fragment count over twice the
    // panel's median haplotype length. Using any candidate's own length here would make the term
    // self-fulfilling.
    double lambda = options.haploid_depth;
    if (options.total_depth && lambda <= 0.0) {
        std::vector<std::size_t> lens;
        for (std::size_t h = 0; h < nh; ++h) lens.push_back(haps[h].seq.size());
        std::sort(lens.begin(), lens.end());
        const double med = lens.empty() ? 0.0 : static_cast<double>(lens[lens.size() / 2]);
        if (med > 0.0) lambda = static_cast<double>(fragments.size()) / (2.0 * med);
    }

    std::vector<HaplotypePairScore> pairs;
    pairs.reserve(nh * (nh + 1) / 2);
    std::vector<std::pair<std::size_t, std::size_t>> pair_index;
    for (std::size_t a = 0; a < nh; ++a) {
        for (std::size_t b = a; b < nh; ++b) pair_index.emplace_back(a, b);
    }
    // PRODUCTION-PATH FRAGMENT DUMP. The existing dump lives inside the joint-depth rescoring block,
    // so the path that actually makes the call could not be inspected per fragment -- which
    // invalidated one diagnostic outright (a partition computed under --joint-depth was read as if it
    // described the production score; the two differ by an order of magnitude). This reproduces the
    // exact term the production loop sums, for one named pair.
    if (!options.dump_fragment_mass.empty() && !options.joint_depth &&
        !options.dump_mass_pair1.empty()) {
        long i1 = -1, i2 = -1;
        for (std::size_t i = 0; i < nh; ++i) {
            if (out.shortlist[i] == options.dump_mass_pair1) i1 = static_cast<long>(i);
            if (out.shortlist[i] == options.dump_mass_pair2) i2 = static_cast<long>(i);
        }
        if (i1 < 0 || i2 < 0) {
            throw std::runtime_error("genotype-frag: --dump-mass-pair named a haplotype that is not "
                                     "in the shortlist; the dump would describe a different pair");
        }
        const std::size_t a = static_cast<std::size_t>(i1), b = static_cast<std::size_t>(i2);
        const double log_total_len =
            options.length_normalize ? log_add(log_len[a], log_len[b]) : std::log(2.0);
        std::ofstream mf(options.dump_fragment_mass);
        if (!mf) throw std::runtime_error("genotype-frag: cannot write " + options.dump_fragment_mass);
        mf.precision(17);
        mf << "# pair\t" << out.shortlist[a] << '\t' << out.shortlist[b] << '\n';
        mf << "# production path (no --joint-depth); contrib is the exact term the pair score sums\n";
        mf << "fragment\tlog_mass\tll_a\tll_b\tfloor\tcontrib\n";
        double check = 0.0;
        for (std::size_t fi = 0; fi < fragments.size(); ++fi) {
            if (floors[fi] == kNegInf) continue;
            const double num = log_add(ll[fi * nh + a], ll[fi * nh + b]);
            const double contrib = log_add(log_mix + num - log_total_len, log_out + floors[fi]);
            check += contrib;
            mf << fragments[fi].name << '\t' << num << '\t' << ll[fi * nh + a] << '\t'
               << ll[fi * nh + b] << '\t' << floors[fi] << '\t' << contrib << '\n';
        }
        // FULL RECONCILIATION. A dump that accounts for only the fragment terms and waves at an
        // unexplained "per-pair constant" cannot be used as a gate: the constant is where the dosage
        // and coverage terms live, and those are exactly the terms under suspicion at a tandem array.
        // Every component is named and the residual is reported, so closure is checkable rather than
        // assumed.
        const double pair_bp_d = static_cast<double>(haps[a].seq.size() + haps[b].seq.size());
        double dosage_d = 0.0;
        if (options.truth_total_bp > 0.0) {
            dosage_d = -std::abs(pair_bp_d - options.truth_total_bp);
        } else if (options.total_depth && lambda > 0.0) {
            const double mean = lambda * pair_bp_d;
            const double n = static_cast<double>(fragments.size());
            dosage_d = n * std::log(mean) - mean - std::lgamma(n + 1.0);
        }
        const double total_d = check + cov_ll[a] + cov_ll[b] + dosage_d;
        mf << "# fragment_sum\t" << check << '\n';
        mf << "# coverage_a\t" << cov_ll[a] << '\n';
        mf << "# coverage_b\t" << cov_ll[b] << '\n';
        mf << "# dosage\t" << dosage_d << '\n';
        mf << "# total_score\t" << total_d << '\n';
        mf.flush();
        if (!mf) throw std::runtime_error("genotype-frag: write failed for " + options.dump_fragment_mass);
    }

    std::vector<double> pair_score(pair_index.size(), 0.0);
    run_parallel(pair_index.size(), options.threads, [&](std::size_t pi) {
        const auto [a, b] = pair_index[pi];
        const double log_total_len =
            options.length_normalize ? log_add(log_len[a], log_len[b]) : std::log(2.0);
        double total = 0.0;
        for (std::size_t fi = 0; fi < fragments.size(); ++fi) {
            if (floors[fi] == kNegInf) continue;
            const double num = log_add(ll[fi * nh + a], ll[fi * nh + b]);
            total += log_add(log_mix + num - log_total_len, log_out + floors[fi]);
        }
        const double pair_bp = static_cast<double>(haps[a].seq.size() + haps[b].seq.size());
        double dosage = 0.0;
        if (options.truth_total_bp > 0.0) {
            // Oracle arm: a sharp penalty on total-length error, standing in for perfect dosage
            // knowledge. Not a model, a bound.
            const double err = std::abs(pair_bp - options.truth_total_bp);
            dosage = -err;
        } else if (options.total_depth && lambda > 0.0) {
            const double mean = lambda * pair_bp;
            const double n = static_cast<double>(fragments.size());
            dosage = n * std::log(mean) - mean - std::lgamma(n + 1.0);
        }
        pair_score[pi] = total + cov_ll[a] + cov_ll[b] + dosage;
    });
    for (std::size_t pi = 0; pi < pair_index.size(); ++pi) {
        pairs.push_back({pair_index[pi].first, pair_index[pi].second, pair_score[pi], 0.0});
    }
    std::stable_sort(pairs.begin(), pairs.end(),
                     [](const HaplotypePairScore& x, const HaplotypePairScore& y) {
                         return x.score > y.score;
                     });

    // Per haplotype, how much of each window the observation process can actually see. A constant
    // lambda*window is wrong in three ways at once and each produces a FALSE empty window, which the
    // model then charges the candidate for:
    //
    //   * a fragment's midpoint cannot fall within half a fragment length of either end, so the
    //     first and last windows are partly unobservable;
    //   * the terminal window is usually shorter than the window size;
    //   * a window containing no anchorable syncmer can never receive a fragment at all, however
    //     much sequence the sample has there -- that is a property of recruitment, not of the sample.
    //
    // Windows of the third kind are marked unobservable and dropped from the likelihood rather than
    // scored as an observed zero.
    std::vector<std::vector<double>> exposure(options.joint_depth ? nh : 0);
    if (options.joint_depth) {
        const std::size_t W = options.joint_window;
        const double half = (options.use_insert_size ? options.fragment_len : 150.0) / 2.0;
        for (std::size_t h = 0; h < nh; ++h) {
            const double L = static_cast<double>(haps[h].seq.size());
            const std::size_t nw = haps[h].seq.size() / W + 1;
            exposure[h].assign(nw, 0.0);
            // which windows hold at least one syncmer that survived the anchor cap
            std::vector<char> anchorable(nw, 0);
            {
                std::unordered_map<std::uint64_t, std::uint32_t> per_hap;
                const std::vector<KmerOccurrence> sy = collect_syncmers(haps[h].seq, k, s);
                for (const KmerOccurrence& o : sy) ++per_hap[o.code];
                for (const KmerOccurrence& o : sy) {
                    if (per_hap[o.code] > options.max_anchor_occ) continue;
                    const std::size_t w = o.start / W;
                    if (w < nw) anchorable[w] = 1;
                }
            }
            for (std::size_t w = 0; w < nw; ++w) {
                if (!anchorable[w]) continue;                     // unobservable: excluded entirely
                const double lo = std::max(static_cast<double>(w * W), half);
                const double hi = std::min(static_cast<double>((w + 1) * W), std::max(0.0, L - half));
                exposure[h][w] = std::max(0.0, hi - lo);
            }
        }
    }

    // ---- joint fragment-assignment + window-depth rescoring ----------------------------------
    if (options.joint_depth && !pairs.empty()) {
        // lambda is SUPPLIED, never fitted. Fitting it from a candidate -- even the best one -- lets
        // the expectation follow the hypothesis it is supposed to test.
        const double lambda_joint = options.haploid_depth;

        // lgamma of small integer counts, precomputed: the inner loop evaluates it millions of times
        // and it is the whole cost of the model otherwise.
        std::vector<double> lgam(4096);
        for (std::size_t i = 0; i < lgam.size(); ++i) lgam[i] = std::lgamma(static_cast<double>(i) + 1.0);
        const auto lgam_of = [&](double n) {
            const std::size_t i = static_cast<std::size_t>(n + 0.5);
            return i < lgam.size() ? lgam[i] : std::lgamma(n + 1.0);
        };

        const std::size_t nrescore = options.joint_top_pairs == 0
            ? pairs.size() : std::min(options.joint_top_pairs, pairs.size());
        std::vector<double> joint(nrescore, kNegInf);
        std::vector<std::size_t> ambiguous(nrescore, 0), nulls(nrescore, 0);
        std::vector<std::size_t> conv_iters(nrescore, 0), conv_moves(nrescore, 0);

        run_parallel(nrescore, options.threads, [&](std::size_t pi) {
            const std::size_t a = pairs[pi].hap1, b = pairs[pi].hap2;
            const std::size_t W = options.joint_window;
            const std::size_t wa = haps[a].seq.size() / W + 1;
            const std::size_t wb = haps[b].seq.size() / W + 1;

            // A homozygous call is TWO copies of the same haplotype, so every window is present
            // twice and must expect twice the reads. Scoring it as one copy would make homozygosity
            // look half-covered everywhere.
            const bool homozygous = (a == b);

            // state: -1 null, otherwise an index into the fragment's placement list, with the
            // homologue it belongs to. Exactly one state per fragment, by construction.
            struct Opt { std::uint8_t side; std::int32_t w; double ll; };
            std::vector<std::vector<Opt>> opts(fragments.size());
            std::vector<int> z(fragments.size(), -1);
            std::vector<double> ca(wa, 0.0), cb(wb, 0.0);

            for (std::size_t fi = 0; fi < fragments.size(); ++fi) {
                if (floors[fi] == kNegInf) continue;
                auto& o = opts[fi];
                const auto gather = [&](std::size_t h, std::uint8_t side, std::size_t nw) {
                    for (const auto& pr : placements[fi * nh + h]) {
                        const std::int32_t mid = pr.mid; const double v = pr.ll + pr.log_mult;
                        if (mid < 0) continue;
                        const std::size_t w = static_cast<std::size_t>(mid) / W;
                        if (w >= nw) continue;
                        // Window observability is part of the window approximation, so rung zero must
                        // not apply it: under the contract every start is a state.
                        if (!options.rung_zero && exposure[h][w] <= 0.0) continue;
                        o.push_back({side, static_cast<std::int32_t>(w), v});
                    }
                };
                gather(a, 0, wa);
                if (!homozygous) gather(b, 1, wb);
                if (o.size() > 1) ++ambiguous[pi];
            }

            const auto wll = [&](double n, double expect) {
                return expect <= 0.0 ? 0.0 : n * std::log(expect) - lgam_of(n);
            };
            const auto expect_of = [&](std::uint8_t side, std::size_t w) {
                const double e = (side == 0 ? exposure[a][w] : exposure[b][w]) * lambda_joint;
                return homozygous ? 2.0 * e : e;
            };

            if (options.joint_marginal) {
                // Observed-data likelihood: sum over placements, no assignment at all.
                //
                // RUNG ZERO uses the CONTRACT's exposure -- SUM_L pi(L) * max(0, |H| - L + 1) -- in
                // place of the window approximation. The window form measures midpoints against
                // anchorable windows and an average fragment length, so its normaliser ranges over a
                // different state space from its event term, which is the invariant the contract
                // exists to enforce. Until they agree, a difference from the reference is a model
                // difference and no ladder of approximations below it means anything.
                //
                // A homozygous pair is TWO copies of one sequence. Its exposure is doubled below, and
                // its event INTENSITY must be doubled with it: a fragment could have come from either
                // copy, so the density it sees is 2 * lambda * sum_p, not lambda * sum_p. Omitting it
                // costs log(2) on every placed fragment -- about 14,000 nats over 20,000 fragments --
                // and silently penalises every homozygous call. Placements are gathered once because
                // the two copies are the same sequence and offer the same positions; the factor, not
                // a second gather, is what represents the second copy.
                double total = 0.0;
                double exposure_total = 0.0;
                if (options.rung_zero) {
                    const double ea = ins_prior.exposure(haps[a].seq.size());
                    const double eb = homozygous ? ea : ins_prior.exposure(haps[b].seq.size());
                    exposure_total = lambda_joint * (homozygous ? 2.0 * ea : ea + eb);
                } else {
                    for (std::size_t w = 0; w < wa; ++w) exposure_total += expect_of(0, w);
                    if (!homozygous) for (std::size_t w = 0; w < wb; ++w) exposure_total += expect_of(1, w);
                }
                total -= exposure_total;
                const double copies = homozygous ? 2.0 : 1.0;
                const double log_lam = std::log(std::max(lambda_joint, 1e-300));
                const double lmix = std::log1p(-options.outlier_mix);
                const double lout = std::log(options.outlier_mix);
                for (std::size_t fi = 0; fi < fragments.size(); ++fi) {
                    if (floors[fi] == kNegInf) continue;
                    double lse = kNegInf;
                    for (const Opt& o : opts[fi]) lse = log_add(lse, o.ll);
                    const double placed = (lse == kNegInf) ? kNegInf
                                        : lmix + log_lam + std::log(copies) + lse;
                    total += log_add(placed, lout + floors[fi]);
                }
                joint[pi] = total;
                conv_iters[pi] = 0;
                return;
            }

            // start from each fragment's best placement, then let depth redistribute it
            for (std::size_t fi = 0; fi < fragments.size(); ++fi) {
                if (opts[fi].empty()) continue;
                std::size_t best = 0;
                for (std::size_t j = 1; j < opts[fi].size(); ++j) {
                    if (opts[fi][j].ll > opts[fi][best].ll) best = j;
                }
                if (opts[fi][best].ll <= floors[fi]) continue;
                z[fi] = static_cast<int>(best);
                const Opt& o = opts[fi][best];
                (o.side == 0 ? ca : cb)[static_cast<std::size_t>(o.w)] += 1.0;
            }

            std::size_t iters_used = 0, last_moves = 0;
            for (std::size_t it = 0; it < options.joint_iterations; ++it) {
                std::size_t moves = 0;
                iters_used = it + 1;
                for (std::size_t oi = 0; oi < fragments.size(); ++oi) {
                    const std::size_t fi = options.joint_reverse_order
                        ? fragments.size() - 1 - oi : oi;
                    if (floors[fi] == kNegInf || opts[fi].empty()) continue;
                    if (z[fi] >= 0) {
                        const Opt& o = opts[fi][static_cast<std::size_t>(z[fi])];
                        (o.side == 0 ? ca : cb)[static_cast<std::size_t>(o.w)] -= 1.0;
                    }
                    int best = -1;
                    double bs = floors[fi];   // the null state
                    for (std::size_t j = 0; j < opts[fi].size(); ++j) {
                        const Opt& o = opts[fi][j];
                        auto& cnt = (o.side == 0 ? ca : cb);
                        const std::size_t w = static_cast<std::size_t>(o.w);
                        const double e = expect_of(o.side, w);
                        const double gain = wll(cnt[w] + 1.0, e) - wll(cnt[w], e);
                        if (o.ll + gain > bs) { bs = o.ll + gain; best = static_cast<int>(j); }
                    }
                    if (best != z[fi]) ++moves;
                    z[fi] = best;
                    if (best >= 0) {
                        const Opt& o = opts[fi][static_cast<std::size_t>(best)];
                        (o.side == 0 ? ca : cb)[static_cast<std::size_t>(o.w)] += 1.0;
                    }
                }
                last_moves = moves;
                if (moves == 0) break;
            }
            conv_iters[pi] = iters_used;
            conv_moves[pi] = last_moves;

            double total = 0.0;
            double assigned = 0.0;
            for (std::size_t fi = 0; fi < fragments.size(); ++fi) {
                if (floors[fi] == kNegInf) continue;
                if (z[fi] < 0) { total += floors[fi]; ++nulls[pi]; continue; }
                total += opts[fi][static_cast<std::size_t>(z[fi])].ll;
                assigned += 1.0;
            }
            // Every OBSERVABLE window of both homologues. Unobservable ones contribute nothing --
            // they are not evidence of absence.
            double summed = 0.0;
            for (std::size_t w = 0; w < wa; ++w) {
                const double e = expect_of(0, w);
                if (e <= 0.0) continue;
                total += ca[w] * std::log(e) - e - lgam_of(ca[w]);
                summed += ca[w];
            }
            if (!homozygous) {
                for (std::size_t w = 0; w < wb; ++w) {
                    const double e = expect_of(1, w);
                    if (e <= 0.0) continue;
                    total += cb[w] * std::log(e) - e - lgam_of(cb[w]);
                    summed += cb[w];
                }
            }
            // INVARIANT: window counts account for exactly the non-null fragments.
            if (std::abs(summed - assigned) > 1e-6) {
                std::fprintf(stderr, "genotype-frag: joint invariant violated -- windows hold %.0f "
                                     "fragments, %.0f were assigned\n", summed, assigned);
                std::abort();
            }
            joint[pi] = total;
        });

        out.convergence.pairs_rescored = nrescore;
        for (std::size_t pi = 0; pi < nrescore; ++pi) {
            if (conv_moves[pi] == 0) ++out.convergence.pairs_converged;
            out.convergence.max_iterations_used =
                std::max(out.convergence.max_iterations_used, conv_iters[pi]);
            out.convergence.total_moves_last_iteration += conv_moves[pi];
        }
        if (!options.dump_fragment_mass.empty() && nrescore > 0) {
            // The top pair's per-fragment placement mass, comparable term-for-term with the
            // reference's. Written before the re-sort so it belongs to a named pair.
            std::size_t a = pairs[0].hap1, b = pairs[0].hap2;
            if (!options.dump_mass_pair1.empty()) {
                // A NAMED pair, so the dump is comparable with a reference dump for the same pair.
                // Defaulting to whatever ranked first lets two files silently describe different
                // diplotypes and be differenced anyway.
                long i1 = -1, i2 = -1;
                for (std::size_t i = 0; i < nh; ++i) {
                    if (out.shortlist[i] == options.dump_mass_pair1) i1 = static_cast<long>(i);
                    if (out.shortlist[i] == options.dump_mass_pair2) i2 = static_cast<long>(i);
                }
                if (i1 < 0 || i2 < 0) {
                    throw std::runtime_error("genotype-frag: --dump-mass-pair named a haplotype that "
                                             "is not in the shortlist; the dump would describe a "
                                             "different pair from the one asked for");
                }
                a = static_cast<std::size_t>(i1);
                b = static_cast<std::size_t>(i2);
            }
            const bool homoz = (a == b);
            // ONE DEFINITION of a haplotype's placement mass, so the dump and the containment
            // certificate cannot describe different quantities. These are the records the final
            // score is built from: placements -> opts -> lse -> joint, and joint overwrites every
            // pair score at the end. ll[] is NOT the source of truth on this path.
            //
            // The records ALREADY carry log_half_strand -- both the coordinate-join and the
            // Cartesian path add it when the record is made -- so they are in the reference's
            // units as they stand. Adding it again here is exactly the -log 2 the first version of
            // the certificate reported as "two strand conventions".
            const auto hap_mass = [&](std::size_t fi2, std::size_t h) {
                double m = kNegInf;
                for (const auto& pr : placements[fi2 * nh + h]) m = log_add(m, pr.ll + pr.log_mult);
                return m;
            };
            std::ofstream mf(options.dump_fragment_mass);
            if (mf) {
                mf.precision(17);
                mf << "# pair\t" << out.shortlist[a] << '\t' << out.shortlist[b] << '\n';
                mf << "fragment\tlog_mass\tmates_seeded\tcontrib\tmates_placed\tvalid_fr\n";
                for (std::size_t fi = 0; fi < fragments.size(); ++fi) {
                    double lse = kNegInf;
                    for (const auto& pr : placements[fi * nh + a]) lse = log_add(lse, pr.ll + pr.log_mult);
                    if (!homoz) {
                        for (const auto& pr : placements[fi * nh + b]) lse = log_add(lse, pr.ll + pr.log_mult);
                    } else if (lse != kNegInf) {
                        lse += std::log(2.0);
                    }
                    // Seeding for THIS pair: a mate counts if it anchored on either homologue.
                    const std::uint8_t ma = mates_seeded[fi * nh + a];
                    const std::uint8_t mb = homoz ? ma : mates_seeded[fi * nh + b];
                    const std::uint8_t both = static_cast<std::uint8_t>(ma | mb);
                    const int nseed  = ((both & 1u) ? 1 : 0) + ((both & 2u) ? 1 : 0);
                    const int nplace = ((both & 4u) ? 1 : 0) + ((both & 8u) ? 1 : 0);
                    const int fr     = (both & 16u) ? 1 : 0;
                    // The full per-fragment CONTRIBUTION, so a reconciliation is possible: exposure
                    // is common to both arms at a fixed pair, so the per-fragment deltas must sum
                    // exactly to the whole-pair difference. A dump of placement mass alone cannot be
                    // reconciled, because a fragment with no placement has mass -inf and its real
                    // contribution is the finite background term.
                    const double contrib = log_add(
                        lse == kNegInf ? kNegInf
                                       : std::log1p(-options.outlier_mix) + std::log(lambda_joint) + lse,
                        std::log(options.outlier_mix) + floors[fi]);
                    mf << fragments[fi].name << '\t' << lse << '\t' << nseed << '\t'
                       << contrib << '\t' << nplace << '\t' << fr << '\n';
                }
            }
            // ---- THE CONTAINMENT CERTIFICATE -------------------------------------------------
            // Production is a LOWER bound on placement mass; the untruncated reference is its
            // UPPER side, summing every (start, L) including matches far outside the declared edit
            // band. A fixed-tolerance equality between them was never the contract, and it fails
            // wherever the band legitimately omits mass -- a fragment whose mates place nowhere IN
            // BAND still has finite Hamming mass at every position.
            //
            // PRODUCTION IS THE AUTHORITATIVE LOWER ENDPOINT. `lse` already carries the real origin
            // multiplicities and orientations, so it is used directly and never recomputed. It
            // does need converting into the REFERENCE's units: fragment_states_mass applies
            // log(0.5) per state because the reference AVERAGES the two library orientations,
            // while lse SUMS them. That is a derivation, not a fitted constant, and the identity
            // check below is what proves it -- with an exhaustive in-band search the converted
            // production mass must equal the enumerated in-band mass exactly.
            //
            // THE ASSERTIONS, in mass units and then in contribution units:
            //   L_h = lse_h + log(0.5)          production, in reference units
            //   T_h = certified omitted bound   same units
            //   U_h = logaddexp(L_h, T_h)
            //   (1) identity   L_h == in-band enumeration   (search completeness)
            //   (2) mass       L_h <= R_h <= U_h            per fragment and haplotype
            //   (3) contrib    C_lower <= C_ref <= C_upper  per fragment
            //   (4) sums       the same, summed like with like
            // The mass-level assertion is the one that matters: background dominance can hide a
            // log 2 mass error completely once the mixture is applied.
            if (!options.dump_containment.empty()) {
                std::ofstream cf(options.dump_containment);
                if (cf) {
                    cf.precision(17);
                    ReferenceParams rp;
                    rp.lambda = lambda_joint;
                    rp.eta = options.outlier_mix;
                    rp.error_rate = options.error_rate;
                    rp.fragment_len = options.fragment_len;
                    rp.fragment_sd = options.fragment_sd;
                    rp.bg_divergence = options.bg_divergence;
                    // The prior is passed in directly (ins_prior), so the parameters that would
                    // rebuild it are not copied here -- one prior, shared, rather than two
                    // constructions that must be kept in step.
                    cf << "# pair\t" << out.shortlist[a] << '\t' << out.shortlist[b] << '\n';
                    cf << "fragment\tL_a\tT_a\tU_a\tR_a\tL_b\tT_b\tU_b\tR_b"
                          "\tC_lower\tC_ref\tC_upper\tident_delta_a\tstratum\n";
                    std::size_t n_ident = 0, n_ident_bad = 0;
                    std::size_t n_both_inf = 0, n_prod_inf_only = 0, n_oracle_inf_only = 0,
                                n_finite_equal = 0, n_finite_unequal = 0;
                    std::size_t n_mass_bad = 0, n_contrib_bad = 0, n_strict = 0;
                    std::size_t s_paired = 0, s_one = 0, s_absent = 0;
                    double sum_lo = 0.0, sum_ref = 0.0, sum_hi = 0.0;
                    const double eps = 1e-9;
                    for (std::size_t fi = 0; fi < fragments.size(); ++fi) {
                        const Fragment& F = fragments[fi];
                        const std::string r1rc = reverse_complement(F.r1);
                        const std::string r2rc =
                            F.r2.empty() ? std::string() : reverse_complement(F.r2);
                        const std::size_t d1 = mate_band_edits(options.max_divergence, F.r1.size());
                        const std::size_t d2 = mate_band_edits(options.max_divergence, F.r2.size());
                        double Lm[2], Tm[2], Um[2], Rm[2];
                        double ident_delta = 0.0;
                        std::size_t hs[2] = {a, b};
                        for (int k = 0; k < 2; ++k) {
                            const std::size_t h = hs[k];
                            Lm[k] = hap_mass(fi, h);   // already in reference units
                            const TailLevel t =
                                tail_interval_level(F.r1, F.r2, r1rc, r2rc, haps[h].seq, d1, d2, 1,
                                                    ins_prior, log_eps3, log_1meps, nullptr);
                            Tm[k] = t.bound;
                            Um[k] = (Lm[k] == kNegInf) ? Tm[k]
                                   : (Tm[k] == kNegInf ? Lm[k] : log_add(Lm[k], Tm[k]));
                            Rm[k] = reference_fragment_on_haplotype(F, haps[h].seq, rp, ins_prior,
                                                                    r2rc, log_eps3, log_1meps);
                            // BOTH haplotypes, and -inf parity reported as its own category.
                            // Printing 0 for an infinite comparison is not a small delta, it is a
                            // different question answered with a number that looks like agreement.
                            ++n_ident;
                            const bool lp_inf = (Lm[k] == kNegInf), tp_inf = (t.lower == kNegInf);
                            if (lp_inf && tp_inf) ++n_both_inf;
                            else if (lp_inf) ++n_prod_inf_only;
                            else if (tp_inf) ++n_oracle_inf_only;
                            else if (std::abs(Lm[k] - t.lower) <= eps) ++n_finite_equal;
                            else ++n_finite_unequal;
                            if (lp_inf != tp_inf || (!lp_inf && std::abs(Lm[k] - t.lower) > eps))
                                ++n_ident_bad;
                            if (k == 0)
                                ident_delta = (lp_inf || tp_inf)
                                    ? std::numeric_limits<double>::quiet_NaN()
                                    : (Lm[0] - t.lower);
                            const bool lo_ok = (Lm[k] == kNegInf) || (Rm[k] >= Lm[k] - eps);
                            const bool hi_ok = (Rm[k] == kNegInf) || (Rm[k] <= Um[k] + eps);
                            if (!lo_ok || !hi_ok) ++n_mass_bad;
                            if (Lm[k] != kNegInf && Rm[k] > Lm[k] + eps && Rm[k] < Um[k] - eps)
                                ++n_strict;
                        }
                        const bool homoz2 = (a == b);
                        const MassInterval ia{Lm[0], Um[0]}, ib{Lm[1], Um[1]};
                        const MassInterval ci =
                            fragment_contribution(ia, ib, homoz2,
                                                  std::log1p(-options.outlier_mix),
                                                  std::log(lambda_joint),
                                                  std::log(options.outlier_mix), floors[fi]);
                        // TWO POINT INTERVALS, one per homologue. Collapsing {R_a, R_b} into a
                        // single interval and passing it for both is harmless on a homozygote and
                        // wrong everywhere else.
                        const MassInterval ra{Rm[0], Rm[0]}, rb2{Rm[1], Rm[1]};
                        const MassInterval cr =
                            fragment_contribution(ra, rb2, homoz2,
                                                  std::log1p(-options.outlier_mix),
                                                  std::log(lambda_joint),
                                                  std::log(options.outlier_mix), floors[fi]);
                        const double c_ref = cr.lower;
                        if (!(c_ref >= ci.lower - eps && c_ref <= ci.upper + eps)) ++n_contrib_bad;
                        sum_lo += ci.lower; sum_ref += c_ref; sum_hi += ci.upper;
                        const std::uint8_t both = static_cast<std::uint8_t>(
                            mates_seeded[fi * nh + a] | mates_seeded[fi * nh + b]);
                        const int nplace = ((both & 4u) ? 1 : 0) + ((both & 8u) ? 1 : 0);
                        const int fr = (both & 16u) ? 1 : 0;
                        const char* stratum = (nplace == 0) ? "absent"
                                            : (fr == 0 ? "one_mate_no_pair" : "paired");
                        if (nplace == 0) ++s_absent; else if (fr == 0) ++s_one; else ++s_paired;
                        cf << F.name << '\t' << Lm[0] << '\t' << Tm[0] << '\t' << Um[0] << '\t'
                           << Rm[0] << '\t' << Lm[1] << '\t' << Tm[1] << '\t' << Um[1] << '\t'
                           << Rm[1] << '\t' << ci.lower << '\t' << c_ref << '\t' << ci.upper
                           << '\t' << ident_delta << '\t' << stratum << '\n';
                    }
                    // MULTIPLICITY-AWARE PRUNES IN-BAND MASS, which would have to be added to
                    // the omitted bound before containment could hold. Refused rather than
                    // silently certified against a bound that does not cover it.
                    cf << "# multiplicity_aware\t" << (options.multiplicity_aware ? 1 : 0) << '\n';
                    cf << "# certificate_valid\t" << (options.multiplicity_aware ? 0 : 1) << '\n';
                    cf << "# identity_checked\t" << n_ident << '\n';
                    cf << "# both_neg_inf\t" << n_both_inf << '\n';
                    cf << "# production_only_neg_inf\t" << n_prod_inf_only << '\n';
                    cf << "# oracle_only_neg_inf\t" << n_oracle_inf_only << '\n';
                    cf << "# finite_equal\t" << n_finite_equal << '\n';
                    cf << "# finite_unequal\t" << n_finite_unequal << '\n';
                    cf << "# identity_failures\t" << n_ident_bad << '\n';
                    cf << "# mass_containment_failures\t" << n_mass_bad << '\n';
                    cf << "# contribution_containment_failures\t" << n_contrib_bad << '\n';
                    cf << "# strict_containments\t" << n_strict << '\n';
                    cf << "# stratum_paired\t" << s_paired << '\n';
                    cf << "# stratum_one_mate_no_pair\t" << s_one << '\n';
                    cf << "# stratum_absent\t" << s_absent << '\n';
                    cf << "# contrib_sum_lower\t" << sum_lo << '\n';
                    cf << "# contrib_sum_reference\t" << sum_ref << '\n';
                    cf << "# contrib_sum_upper\t" << sum_hi << '\n';
                    cf << "# sum_containment\t"
                       << ((sum_ref >= sum_lo - 1e-6 && sum_ref <= sum_hi + 1e-6) ? 1 : 0) << '\n';
                }
            }
        }
        for (std::size_t pi = 0; pi < nrescore; ++pi) pairs[pi].score = joint[pi];
        std::stable_sort(pairs.begin(), pairs.begin() + static_cast<long>(nrescore),
                         [](const HaplotypePairScore& x, const HaplotypePairScore& y) {
                             return x.score > y.score;
                         });
    }

    double norm = kNegInf;
    for (const HaplotypePairScore& p : pairs) norm = log_add(norm, p.score);
    for (HaplotypePairScore& p : pairs) p.posterior = std::exp(p.score - norm);

    // ---- what the evidence actually distinguishes ---------------------------------------------
    {
        EquivalenceSet& eq = out.equivalence;
        const double best = pairs.front().score;
        eq.margin = pairs.size() > 1 ? best - pairs[1].score : std::numeric_limits<double>::infinity();
        for (const HaplotypePairScore& p : pairs) {
            if (best - p.score > options.equivalence_tolerance) break;
            ++eq.size;
            eq.posterior_mass += p.posterior;
            if (eq.members.size() < options.equivalence_max_report) {
                eq.members.push_back(out.shortlist[p.hap1] + "/" + out.shortlist[p.hap2]);
            }
        }
        eq.blocks_total = chain.size();
    }

    // ---- project onto blocks ------------------------------------------------------------------
    out.blocks.resize(chain.size());
    for (std::size_t bi = 0; bi < chain.size(); ++bi) {
        BlockProjection& P = out.blocks[bi];
        P.block_index = bi;
        P.kind = chain[bi].kind;
        P.bubble_id = chain[bi].bubble_id;
        P.n_alleles = blocks[bi].allele_seq.size();
        // The posterior mass on each allele pair is computed either way, because it is what the
        // reported confidence means; only which one is CALLED depends on the projection.
        std::map<std::pair<int, int>, double> mass;
        for (const HaplotypePairScore& p : pairs) {
            if (p.posterior < 1e-12) continue;
            int a = haps[p.hap1].allele[bi];
            int b = haps[p.hap2].allele[bi];
            if (a > b) std::swap(a, b);
            mass[{a, b}] += p.posterior;
        }
        if (options.project_map && !pairs.empty()) {
            // No sort: allele1 stays with the pair's FIRST haplotype at every block, so the two
            // columns are two homologues rather than two sorted indices.
            P.allele1 = haps[pairs.front().hap1].allele[bi];
            P.allele2 = haps[pairs.front().hap2].allele[bi];
            const auto it = mass.find({std::min(P.allele1, P.allele2), std::max(P.allele1, P.allele2)});
            P.posterior = it == mass.end() ? 0.0 : it->second;
        } else {
            double best = -1.0;
            for (const auto& [key, m] : mass) {
                if (m > best) { best = m; P.allele1 = key.first; P.allele2 = key.second; }
            }
            P.posterior = best < 0.0 ? 0.0 : best;
        }
        // Determined = every member of the equivalence set agrees here. Computed even when the set
        // has one member, where it is trivially true, so the column means one thing throughout.
        {
            const double best = pairs.front().score;
            bool any_inexact = false;   // any member of the equivalence set the decomposition fails
            int fa = -2, fb = -2;
            // BLOCK EQUIVALENCE SIZE: the number of DISTINCT unordered allele pairs this block is
            // given by every locus pair within tolerance -- every one of them, not the truncated
            // reported-member list, which stops at equivalence_max_report and would understate the
            // ambiguity by however many members were not printed.
            std::set<std::pair<int, int>> block_pairs;
            for (const HaplotypePairScore& p : pairs) {
                if (best - p.score > options.equivalence_tolerance) break;
                // EVERY member contributing to the conclusion must be projectable AT THIS BLOCK.
                // Checking only the rank-one pair would declare a block determined on the strength
                // of members whose alleles here cannot be trusted; checking the whole haplotype
                // would discard blocks that are perfectly well mapped.
                const bool pa = bi < proj_block[p.hap1].size() && proj_block[p.hap1][bi];
                const bool pb = bi < proj_block[p.hap2].size() && proj_block[p.hap2][bi];
                if (!pa || !pb) any_inexact = true;
                int x = haps[p.hap1].allele[bi], y = haps[p.hap2].allele[bi];
                if (x > y) std::swap(x, y);
                if (x < 0 || y < 0) any_inexact = true;
                else if (pa && pb) block_pairs.insert({x, y});
                if (fa == -2) { fa = x; fb = y; }
            }
            const bool agree = block_pairs.size() == 1;
            // NA, not a count, when any contributing member is unprojectable here: a number would
            // be a claim about members whose alleles were never established.
            P.block_equivalence_size = any_inexact ? -1
                                                   : static_cast<long>(block_pairs.size());
            // A block where the selected pair has no allele is UNPROJECTABLE, not determined.
            // Agreement on -1 is agreement about nothing: the decomposition does not cover this path
            // here, so there is no block genotype to report. Emitting one would be a call the
            // evidence cannot support, which is exactly what the block-level product must not do.
            P.unprojectable = (fa < 0 || fb < 0) || any_inexact;
            // determined = projectable AND the block equivalence size is numeric and exactly 1.
            P.determined = !P.unprojectable && P.block_equivalence_size == 1;
            if (P.determined) ++out.equivalence.blocks_determined;
        }

        if (truth_allele1 != nullptr && truth_allele2 != nullptr &&
            bi < truth_allele1->size() && bi < truth_allele2->size()) {
            P.truth_a = (*truth_allele1)[bi];
            P.truth_b = (*truth_allele2)[bi];
            if (P.truth_a >= 0 && P.truth_b >= 0) {
                P.truth_representable = true;
                // Compared UNORDERED. The call is phased and the truth columns are not, so requiring
                // the same order would score a correct call wrong for the order it was written in.
                P.exact = (std::min(P.truth_a, P.truth_b) == std::min(P.allele1, P.allele2) &&
                           std::max(P.truth_a, P.truth_b) == std::max(P.allele1, P.allele2));
            }
        }
    }

    for (const auto& [n1, n2] : options.probe_pairs) {
        HaplotypeProbe pr;
        pr.name1 = n1;
        pr.name2 = n2;
        const auto find = [&](const std::string& n) {
            for (std::size_t i = 0; i < nh; ++i) if (out.shortlist[i] == n) return static_cast<long>(i);
            return -1L;
        };
        const long i1 = find(n1);
        const long i2 = find(n2);
        if (i1 < 0 || i2 < 0) { out.probes.push_back(pr); continue; }
        pr.in_shortlist = true;
        pr.placed1 = out.haplotypes[static_cast<std::size_t>(i1)].placed;
        pr.placed2 = out.haplotypes[static_cast<std::size_t>(i2)].placed;
        const std::size_t a = std::min<std::size_t>(i1, i2);
        const std::size_t b = std::max<std::size_t>(i1, i2);
        for (std::size_t r = 0; r < pairs.size(); ++r) {
            if (pairs[r].hap1 != a || pairs[r].hap2 != b) continue;
            pr.rank = static_cast<int>(r) + 1;
            pr.score = pairs[r].score;
            pr.delta = pairs[r].score - pairs.front().score;
            break;
        }
        out.probes.push_back(pr);
    }

    for (std::size_t i = 0; i < std::min(top_pairs_kept, pairs.size()); ++i) {
        out.top_pairs.push_back(pairs[i]);
    }
    return out;
}

MosaicFloors mosaic_floors(const std::vector<BlockAlleles>& blocks,
                           const std::vector<std::string>& haplotype_names,
                           const std::vector<std::string>& truth_block_seq,
                           const std::vector<double>& switch_penalties,
                           std::size_t threads) {
    MosaicFloors out;
    out.blocks = blocks.size();
    const std::size_t nb = blocks.size();
    const std::size_t nh = haplotype_names.size();

    // cost[block][haplotype] = edit distance from that haplotype's allele here to the truth's own
    // sequence here. Global alignment, not an alignment-derived approximation: these are block
    // alleles, short enough that the exact distance is affordable and the right thing to use.
    std::vector<std::vector<std::size_t>> cost(nb);
    std::vector<char> scored(nb, 0);
    run_parallel(nb, threads, [&](std::size_t bi) {
        cost[bi].assign(nh, 0);
        if (bi >= truth_block_seq.size()) return;
        const std::string& truth = truth_block_seq[bi];
        // Distinct alleles are aligned ONCE and reused across every haplotype carrying them; at a
        // locus with 466 haplotypes over a handful of alleles per block that is the whole cost.
        std::unordered_map<std::size_t, std::size_t> by_allele;
        bool any = false;
        for (std::size_t h = 0; h < nh; ++h) {
            const auto it = blocks[bi].allele_of.find(haplotype_names[h]);
            if (it == blocks[bi].allele_of.end() || it->second >= blocks[bi].allele_seq.size()) {
                // The haplotype does not traverse: it spells nothing, so it costs the truth's whole
                // length here. That is a real cost, not missing data.
                cost[bi][h] = truth.size();
                continue;
            }
            any = true;
            const std::size_t ai = it->second;
            const auto cached = by_allele.find(ai);
            if (cached != by_allele.end()) { cost[bi][h] = cached->second; continue; }
            const std::string& cand = blocks[bi].allele_seq[ai];
            std::size_t d;
            if (cand.empty() || truth.empty()) d = cand.size() + truth.size();
            else d = nw_edit_distance(cand, truth).edits;
            by_allele.emplace(ai, d);
            cost[bi][h] = d;
        }
        scored[bi] = (any && !truth.empty()) ? 1 : 0;
    });
    for (std::size_t bi = 0; bi < nb; ++bi) if (scored[bi]) ++out.blocks_scored;

    // complete: one haplotype for the whole locus
    std::size_t best_total = std::numeric_limits<std::size_t>::max();
    for (std::size_t h = 0; h < nh; ++h) {
        std::size_t t = 0;
        for (std::size_t bi = 0; bi < nb; ++bi) t += cost[bi][h];
        best_total = std::min(best_total, t);
    }
    out.complete = best_total == std::numeric_limits<std::size_t>::max() ? 0 : best_total;

    // free: the nearest allele at every block. LEXICOGRAPHIC -- minimise edits first, then minimise
    // SWITCHES among the paths achieving that minimum.
    //
    // Taking the first minimum-cost haplotype independently at each block gives the right edit total
    // and an arbitrary switch count: wherever several haplotypes carry the same optimal allele, the
    // choice among them is free, and picking by index invents switches that no path needs. The switch
    // count is the number that says whether a mosaic gain is buyable with evidence, so an inflated one
    // makes the mosaic layer look less attainable than it is.
    {
        std::size_t total = 0;
        for (std::size_t bi = 0; bi < nb; ++bi) {
            std::size_t best = std::numeric_limits<std::size_t>::max();
            for (std::size_t h = 0; h < nh; ++h) best = std::min(best, cost[bi][h]);
            if (best != std::numeric_limits<std::size_t>::max()) total += best;
        }
        out.free_mosaic = total;

        // Among optimal-cost paths only, a DP whose cost is the switch count. A haplotype is a legal
        // state at a block only if it attains that block's minimum; staying is free and moving costs
        // one, so this is the minimum number of source changes any optimal path requires.
        const std::size_t kInf = std::numeric_limits<std::size_t>::max() / 4;
        std::vector<std::size_t> dp(nh, kInf);
        bool started = false;
        for (std::size_t bi = 0; bi < nb; ++bi) {
            std::size_t best = std::numeric_limits<std::size_t>::max();
            for (std::size_t h = 0; h < nh; ++h) best = std::min(best, cost[bi][h]);
            if (best == std::numeric_limits<std::size_t>::max()) continue;
            if (!started) {
                for (std::size_t h = 0; h < nh; ++h) dp[h] = (cost[bi][h] == best) ? 0 : kInf;
                started = true;
                continue;
            }
            std::size_t global = kInf;
            for (std::size_t h = 0; h < nh; ++h) global = std::min(global, dp[h]);
            std::vector<std::size_t> next(nh, kInf);
            for (std::size_t h = 0; h < nh; ++h) {
                if (cost[bi][h] != best) continue;              // not an optimal state here
                const std::size_t stay = dp[h];
                const std::size_t move = global == kInf ? kInf : global + 1;
                next[h] = std::min(stay, move);
            }
            dp.swap(next);
        }
        std::size_t sw = kInf;
        for (std::size_t h = 0; h < nh; ++h) sw = std::min(sw, dp[h]);
        out.switches_free = (sw >= kInf || !started) ? 0 : sw;
    }

    // penalised: Viterbi over source haplotype, with a cost for switching. The min-plus structure
    // lets each block be done in O(n) rather than O(n^2): staying costs dp[h], switching costs
    // (best over all h) + penalty.
    for (const double pen : switch_penalties) {
        std::vector<double> dp(nh, 0.0);
        std::vector<std::size_t> sw_count(nh, 0);
        for (std::size_t bi = 0; bi < nb; ++bi) {
            double best = std::numeric_limits<double>::infinity();
            std::size_t best_h = 0;
            for (std::size_t h = 0; h < nh; ++h) {
                if (dp[h] < best) { best = dp[h]; best_h = h; }
            }
            std::vector<double> next(nh);
            std::vector<std::size_t> next_sw(nh);
            for (std::size_t h = 0; h < nh; ++h) {
                const double stay = dp[h];
                const double move = best + pen;
                if (bi > 0 && move < stay) { next[h] = move + cost[bi][h]; next_sw[h] = sw_count[best_h] + 1; }
                else { next[h] = stay + cost[bi][h]; next_sw[h] = sw_count[h]; }
            }
            dp.swap(next);
            sw_count.swap(next_sw);
        }
        double best = std::numeric_limits<double>::infinity();
        std::size_t best_h = 0;
        for (std::size_t h = 0; h < nh; ++h) if (dp[h] < best) { best = dp[h]; best_h = h; }
        out.penalised.emplace_back(pen, static_cast<std::size_t>(best + 0.5));
        out.penalised_switches.push_back(sw_count[best_h]);
    }
    return out;
}

void write_block_calls(const std::string& path, const std::vector<FragmentBlockCall>& calls) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("genotype: cannot write " + path);
    f << "block\tkind\tbubble_id\tstatus\tallele1\tallele2\tallele_set\tdiploid_dosage\t"
         "dosage_confidence\tallocation_determined\tlocal_gq\tlinkage_gq\tfit\t"
         "fragments_informing\n";
    const auto num = [](double v) { return v < 0.0 ? std::string(".") : std::to_string(v); };
    const auto idx = [](int v) { return v < 0 ? std::string(".") : std::to_string(v); };
    for (const FragmentBlockCall& c : calls) {
        const char* kind = c.kind == BlockKind::Bubble ? "bubble"
                         : c.kind == BlockKind::Backbone ? "backbone" : "flank";
        const char* st = c.status == FragmentBlockCallStatus::Called ? "CALLED"
                       : c.status == FragmentBlockCallStatus::Equivalent ? "EQUIVALENT"
                       : c.status == FragmentBlockCallStatus::DosageOnly ? "DOSAGE_ONLY" : "OFF_PANEL";
        f << c.block_index << '\t' << kind << '\t'
          << (c.bubble_id < 0 ? std::string(".") : std::to_string(c.bubble_id)) << '\t' << st
          << '\t' << idx(c.allele1) << '\t' << idx(c.allele2) << '\t';
        if (c.allele_set.empty()) {
            f << '.';
        } else {
            for (std::size_t i = 0; i < c.allele_set.size(); ++i) {
                f << (i ? "," : "") << c.allele_set[i].first << '/' << c.allele_set[i].second;
            }
        }
        f << '\t' << idx(c.diploid_dosage) << '\t' << num(c.dosage_confidence) << '\t'
          << (c.allocation_determined ? "yes" : "no") << '\t' << num(c.local_gq) << '\t'
          << num(c.linkage_gq) << '\t' << num(c.fit) << '\t' << c.fragments_informing << '\n';
    }
    f.flush();
    if (!f) throw std::runtime_error("genotype: write failed for " + path);
}

void spell_called_pair(const std::vector<BlockAlleles>& blocks,
                       const std::vector<int>& allele1,
                       const std::vector<int>& allele2,
                       std::string& seq1,
                       std::string& seq2) {
    seq1.clear();
    seq2.clear();
    for (std::size_t bi = 0; bi < blocks.size(); ++bi) {
        const auto take = [&](const std::vector<int>& a, std::string& into) {
            if (bi >= a.size() || a[bi] < 0) return;
            const std::size_t ai = static_cast<std::size_t>(a[bi]);
            if (ai < blocks[bi].allele_seq.size()) into += blocks[bi].allele_seq[ai];
        };
        take(allele1, seq1);
        take(allele2, seq2);
    }
}

void write_haplotype_results(const std::string& out_prefix,
                             const HaplotypeResult& result,
                             bool have_truth,
                             const std::string& catalogue_fingerprint) {
    const std::string bp = out_prefix + ".hap_blocks.tsv";
    std::ofstream bf(bp);
    if (!bf) throw std::runtime_error("genotype-frag: cannot write " + bp);
    // Provenance first, so --spell-calls can refuse a table whose allele indices index a different
    // catalogue than the one it is about to spell against.
    if (!catalogue_fingerprint.empty()) {
        bf << "# panvar-allele-catalogue\t" << catalogue_fingerprint << '\n';
    }
    // Appended, never inserted: adding a column mid-header once silently changed the meaning of
    // three assertions that indexed by position.
    bf << "block\tkind\tbubble_id\tn_alleles\tallele1\tallele2\tposterior\tdetermined"
          "\tprojectable\tblock_equivalence_size";
    if (have_truth) bf << "\ttruth_a\ttruth_b\trepresentable\texact";
    bf << '\n';
    for (const BlockProjection& p : result.blocks) {
        const char* kind = p.kind == BlockKind::Bubble ? "bubble"
                         : p.kind == BlockKind::Backbone ? "backbone" : "flank";
        bf << p.block_index << '\t' << kind << '\t' << p.bubble_id << '\t' << p.n_alleles << '\t'
           << p.allele1 << '\t' << p.allele2 << '\t' << p.posterior << '\t'
           << (p.determined ? 1 : 0) << '\t' << (p.unprojectable ? "NA" : "yes") << '\t';
        if (p.block_equivalence_size < 0) bf << "NA"; else bf << p.block_equivalence_size;
        if (have_truth) {
            bf << '\t' << p.truth_a << '\t' << p.truth_b << '\t' << (p.truth_representable ? 1 : 0)
               << '\t' << (p.truth_representable ? (p.exact ? "1" : "0") : "NA");
        }
        bf << '\n';
    }
    bf.flush();
    if (!bf) throw std::runtime_error("genotype-frag: write failed for " + bp);

    const std::string hp = out_prefix + ".hap_scores.tsv";
    std::ofstream hf(hp);
    if (!hf) throw std::runtime_error("genotype-frag: cannot write " + hp);
    hf << "haplotype\tbp\tcontainment\tplaced_fragments\twindows\tzero_windows\tcoverage_ll\tsolo_ll\n";
    for (const HaplotypeScore& h : result.haplotypes) {
        hf << h.name << '\t' << h.bp << '\t' << h.containment << '\t' << h.placed << '\t'
           << h.windows << '\t' << h.zero_windows << '\t' << h.coverage_ll << '\t' << h.solo_ll << '\n';
    }
    hf.flush();
    if (!hf) throw std::runtime_error("genotype-frag: write failed for " + hp);

    if (!result.probes.empty()) {
        const std::string qp = out_prefix + ".hap_probes.tsv";
        std::ofstream qf(qp);
        if (!qf) throw std::runtime_error("genotype-frag: cannot write " + qp);
        qf << "hap1\thap2\tin_shortlist\trank\tscore\tdelta\tplaced1\tplaced2\n";
        for (const HaplotypeProbe& pr : result.probes) {
            qf << pr.name1 << '\t' << pr.name2 << '\t' << (pr.in_shortlist ? 1 : 0) << '\t'
               << pr.rank << '\t' << pr.score << '\t' << pr.delta << '\t' << pr.placed1 << '\t'
               << pr.placed2 << '\n';
        }
        qf.flush();
        if (!qf) throw std::runtime_error("genotype-frag: write failed for " + qp);
    }

    {
        const std::string ep = out_prefix + ".equivalence.tsv";
        std::ofstream ef(ep);
        if (!ef) throw std::runtime_error("genotype-frag: cannot write " + ep);
        ef << "margin\tset_size\tposterior_mass\tblocks_determined\tblocks_total\n";
        ef << result.equivalence.margin << '\t' << result.equivalence.size << '\t'
           << result.equivalence.posterior_mass << '\t' << result.equivalence.blocks_determined
           << '\t' << result.equivalence.blocks_total << '\n';
        ef << "# members within tolerance of the best score\n";
        for (const std::string& m : result.equivalence.members) ef << m << '\n';
        ef.flush();
        if (!ef) throw std::runtime_error("genotype-frag: write failed for " + ep);
    }

    const std::string pp = out_prefix + ".hap_pairs.tsv";
    std::ofstream pf(pp);
    if (!pf) throw std::runtime_error("genotype-frag: cannot write " + pp);
    // Full precision: the reconciliation identity is exact arithmetic, so a 6-significant-digit score
    // makes a correct decomposition look like a mismatch at the 1e-6 gate.
    pf.precision(17);
    pf << "rank\thap1\thap2\tscore\tdelta\tposterior\n";
    const double best = result.top_pairs.empty() ? 0.0 : result.top_pairs.front().score;
    for (std::size_t i = 0; i < result.top_pairs.size(); ++i) {
        const HaplotypePairScore& p = result.top_pairs[i];
        pf << (i + 1) << '\t' << result.shortlist[p.hap1] << '\t' << result.shortlist[p.hap2] << '\t'
           << p.score << '\t' << (p.score - best) << '\t' << p.posterior << '\n';
    }
    pf.flush();
    if (!pf) throw std::runtime_error("genotype-frag: write failed for " + pp);
}

// =================================================================================================
// EXACT REFERENCE SCORER
// =================================================================================================

namespace {

// Hamming emission. The reference deliberately does NOT align: a placement is a start position, and
// the fragment either agrees with the sequence there or it does not. Indels inside a read would need
// alignment and would make "the set of starts" ambiguous, which is exactly the ambiguity the contract
// removes. Fixtures are built without them.
// `log_eps` here is log(error_rate/3): the probability of a SPECIFIC mismatching base, matching the
// simulator, which substitutes uniformly among the other three.
double reference_emission(const std::string& read, const std::string& target, std::size_t offset,
                          double log_eps, double log_1meps) {
    if (offset + read.size() > target.size()) return kNegInf;
    std::size_t mism = 0;
    for (std::size_t i = 0; i < read.size(); ++i) {
        if (read[i] != target[offset + i]) ++mism;
    }
    return static_cast<double>(mism) * log_eps +
           static_cast<double>(read.size() - mism) * log_1meps;
}

// The discrete insert prior pi(L), NORMALISED over exactly the range summed, and using the same
// concordant/discordant mixture the accelerated path uses. Both properties matter:
//
//   * unnormalised, the event term and the exposure disagree by a constant that depends on the range,
//     and the "likelihood" is then not a density over anything;
//   * a pure Gaussian here against a mixture in the fast path is a model difference, so a discrepancy
//     between them could not be attributed to acceleration.


// log sum over every start on ONE haplotype of P(f | start), with the insert prior integrated.
// One orientation of one fragment against one haplotype: `lead` maps forward at the start and
// `trail` maps forward at start + L - |trail|.
double reference_orientation(const std::string& lead, const std::string& trail,
                             const std::string& hap, const ReferenceParams& p,
                             const InsertPrior& ip,
                             double log_eps, double log_1meps) {
    if (hap.empty() || lead.empty()) return kNegInf;
    double acc = kNegInf;
    const long n = static_cast<long>(hap.size());
    // The state is (start, L). Exactly these states are what the exposure counts, so the event term
    // and its normaliser range over the same space -- the invariant the whole contract turns on.
    for (long s = 0; s + static_cast<long>(lead.size()) <= n; ++s) {
        const double e1 = reference_emission(lead, hap, static_cast<std::size_t>(s), log_eps, log_1meps);
        if (e1 == kNegInf) continue;
        if (trail.empty()) { acc = log_add(acc, e1); continue; }
        for (long L = ip.lo; L <= ip.hi; ++L) {
            if (s + L > n) break;                  // this start cannot host this insert length
            const long m2 = s + L - static_cast<long>(trail.size());
            if (m2 < 0) continue;
            const double e2 = reference_emission(trail, hap, static_cast<std::size_t>(m2), log_eps, log_1meps);
            if (e2 == kNegInf) continue;
            acc = log_add(acc, e1 + e2 + ip.log_at(L));
        }
    }
    return acc;
}

}  // namespace -- reference_fragment_on_haplotype is EXPORTED: the bounded search
   // needs the untruncated reference as the upper side of its in-band interval,
   // and a second copy of that integral in the command would be the duplicated-rule
   // failure this branch has paid for three times.
double reference_fragment_on_haplotype(const Fragment& f, const std::string& hap,
                                       const ReferenceParams& p, const InsertPrior& ip,
                                       const std::string& r2rc,
                                       double log_eps, double log_1meps) {
    if (hap.empty() || f.r1.empty()) return kNegInf;
    // BOTH STRANDS. A fragment arises from either strand at a start with probability 1/2 each, so a
    // haplotype and its reverse complement must score identically -- which they do not if only the
    // FR orientation is tried. Measured need: one of cyp2d6 HG04036's haplotypes is spelled
    // reverse-complemented by block concatenation, so a scorer that assumes one strand is wrong on
    // real panel sequence, not only on a contrived fixture.
    const std::string r1rc = reverse_complement(f.r1);
    const double fwd = reference_orientation(f.r1, r2rc, hap, p, ip, log_eps, log_1meps);
    const double rev = reference_orientation(f.r2, r1rc, hap, p, ip, log_eps, log_1meps);
    const double half = std::log(0.5);
    return log_add(fwd == kNegInf ? kNegInf : half + fwd,
                   rev == kNegInf ? kNegInf : half + rev);
}

double reference_pair_loglik(const std::string& hap_a, const std::string& hap_b,
                             const std::vector<Fragment>& fragments,
                             const ReferenceParams& params,
                             std::vector<double>* fragment_mass,
                             std::vector<double>* fragment_contrib) {
    // error_rate is the TOTAL substitution probability; a SPECIFIC observed mismatch has
    // probability error_rate/3, which is what the simulator draws. Must match the accelerated
    // scorer exactly or the two compute different likelihoods and no differential between them
    // is meaningful -- the per-fragment reconciliation would still close, because it sums each
    // model's own contributions.
    const double log_eps = std::log(params.error_rate / 3.0);
    const double log_1meps = std::log1p(-params.error_rate);

    // One prior, shared by the event term and the exposure. Its lower bound is the longest fragment
    // present, because an insert shorter than the two mates is not a state at all.
    const long min_len = fragment_insert_floor(fragments, params.allow_overlapping_pairs);
    const InsertPrior ip = make_insert_prior(params.fragment_len, params.fragment_sd,
                                            params.discordant_rate, params.insert_sigmas, min_len);

    // Exposure over the SAME (start, L) states the event term sums:
    //
    //     E_h = SUM over L of pi(L) * max(0, |H_h| - L + 1)
    //
    // Counting starts that admit only the SHORTEST insert credits a start near the end of a
    // haplotype with lengths it cannot host, so the normaliser would cover states the event term
    // never visits.
    const auto exposure_of = [&](const std::string& h) {
        const long n = static_cast<long>(h.size());
        double e = 0.0;
        for (long L = ip.lo; L <= ip.hi; ++L) {
            const double starts = static_cast<double>(std::max<long>(0, n - L + 1));
            if (starts > 0.0) e += std::exp(ip.log_at(L)) * starts;
        }
        return e;
    };
    const double exposure = exposure_of(hap_a) + exposure_of(hap_b);

    double total = -params.lambda * exposure;
    const double log_lam = std::log(params.lambda);
    const double log_mix = std::log1p(-params.eta);
    const double log_bg_w = std::log(params.eta);

    for (const Fragment& f : fragments) {
        const std::size_t len = f.bases();
        if (len == 0) continue;
        const std::string r2rc = f.r2.empty() ? std::string() : reverse_complement(f.r2);
        // The fragment may have come from either homologue: one sum over the union of their starts.
        double lse = reference_fragment_on_haplotype(f, hap_a, params, ip, r2rc, log_eps, log_1meps);
        lse = log_add(lse, reference_fragment_on_haplotype(f, hap_b, params, ip, r2rc, log_eps, log_1meps));

        const double bg = static_cast<double>(
                              static_cast<std::size_t>(params.bg_divergence * static_cast<double>(len))) * log_eps +
                          static_cast<double>(len - static_cast<std::size_t>(
                              params.bg_divergence * static_cast<double>(len))) * log_1meps;
        if (fragment_mass != nullptr) fragment_mass->push_back(lse);
        const double placed = (lse == kNegInf) ? kNegInf : log_mix + log_lam + lse;
        const double contrib = log_add(placed, log_bg_w + bg);
        if (fragment_contrib != nullptr) fragment_contrib->push_back(contrib);
        total += contrib;
    }
    return total;
}

namespace {
// Does the decomposition assign this path an allele at this block? Absence means the path does not
// traverse the block at all, which is different from traversing it and spelling nothing.
bool sl_block_has_no_allele(const std::vector<BlockAlleles>& blocks, std::uint32_t bi,
                            const std::string& name) {
    if (bi >= blocks.size()) return true;
    return blocks[bi].allele_of.find(name) == blocks[bi].allele_of.end();
}
}  // namespace

double omitted_mass_bound(std::size_t hap_len, std::size_t m1_len, std::size_t m2_len,
                          std::size_t d1, std::size_t d2, const InsertPrior& ip,
                          double log_eps, double log_1meps,
                          const std::vector<FragmentState>& found) {
    const long n = static_cast<long>(hap_len);
    // Found states per insert length, so the omitted count is exact at every L rather than lumped.
    std::unordered_map<long, std::size_t> found_at;
    for (const FragmentState& st : found) ++found_at[st.insert];

    const auto read_ll = [&](std::size_t edits, std::size_t len) {
        const std::size_t e = std::min(edits, len);
        return static_cast<double>(e) * log_eps + static_cast<double>(len - e) * log_1meps;
    };
    // ONE mate perfect, the other exactly one edit past its band. Both mates out is strictly worse.
    const double B = std::max(read_ll(d1 + 1, m1_len) + read_ll(0, m2_len),
                              read_ll(0, m1_len) + read_ll(d2 + 1, m2_len));

    // log SUM_L N_omitted(L) * pi(L), accumulated in log space.
    double acc = kNegInf;
    for (long L = ip.lo; L <= ip.hi; ++L) {
        const double starts = static_cast<double>(std::max<long>(0, n - L + 1));
        if (starts <= 0.0) continue;
        const double n_total_L = 2.0 * starts;               // the two library orientations
        const auto it = found_at.find(L);
        const double n_found_L = it == found_at.end() ? 0.0 : static_cast<double>(it->second);
        const double n_om = n_total_L - n_found_L;
        if (n_om <= 0.0) continue;
        acc = log_add(acc, std::log(n_om) + ip.log_at(L));
    }
    if (acc == kNegInf) return kNegInf;
    return std::log(0.5) + B + acc;
}

MassInterval fragment_contribution(const MassInterval& ma, const MassInterval& mb, bool homozygous,
                                   double log_mix, double log_lambda, double log_bg_weight,
                                   double log_p_bg) {
    const double log_two = std::log(2.0);
    const auto one = [&](double a, double b) {
        // TWO CHROMOSOME COPIES for a homozygote: 2*M_a, i.e. log M_a + log 2. Not M_a + M_a mixed
        // twice with the background -- that would count eta*P_bg once per copy.
        const double m = homozygous ? (a == kNegInf ? kNegInf : a + log_two) : log_add(a, b);
        const double sig = (m == kNegInf) ? kNegInf : log_mix + log_lambda + m;
        return log_add(sig, log_bg_weight + log_p_bg);
    };
    MassInterval out;
    out.lower = one(ma.lower, mb.lower);
    out.upper = one(ma.upper, mb.upper);
    return out;
}

MassInterval representative_total(const MassInterval& fragment_sum, double exposure_a,
                                  double exposure_b, bool homozygous, double lambda,
                                  double log_prior) {
    // ONCE, outside the fragment sum. Charging it per fragment would multiply it by the read count.
    const double e = homozygous ? 2.0 * exposure_a : exposure_a + exposure_b;
    const double term = -lambda * e + log_prior;
    MassInterval out;
    out.lower = fragment_sum.lower + term;
    out.upper = fragment_sum.upper + term;
    return out;
}

Certification certify(const std::vector<MassInterval>& classes, double tau) {
    Certification out;
    out.best_lower = kNegInf;
    for (const MassInterval& c : classes) out.best_lower = std::max(out.best_lower, c.lower);
    // tau applied in ONE place: a class stays plausible if its upper reaches within tau of the best
    // established lower bound. Mixing a strict test here with an approximate one elsewhere is how a
    // tolerance silently stops meaning anything.
    const double threshold = out.best_lower - tau;
    for (std::size_t i = 0; i < classes.size(); ++i) {
        if (classes[i].upper >= threshold) out.plausible.push_back(i);
    }
    out.certified = out.plausible.size() == 1;
    return out;
}

const char* verdict_name(Verdict v) {
    switch (v) {
        case Verdict::Certified:  return "CERTIFIED";
        case Verdict::Unresolved: return "UNRESOLVED";
        default:                  return "INCOMPLETE";
    }
}

Verdict verdict_of(const std::vector<MassInterval>& classes, const Certification& cert,
                   double tol) {
    if (cert.certified) return Verdict::Certified;
    // UNRESOLVED only if every plausible class actually reached the tolerance. Otherwise the
    // overlap may be an artefact of stopping early, and reporting it as ambiguity would attribute
    // a resource limit to the data.
    for (std::size_t i : cert.plausible) {
        if (classes[i].upper - classes[i].lower > tol) return Verdict::Incomplete;
    }
    return Verdict::Unresolved;
}

MassInterval aggregate_class(const std::vector<MassInterval>& reps,
                             const std::vector<double>& weights) {
    MassInterval out;
    out.lower = kNegInf;
    out.upper = kNegInf;
    for (std::size_t i = 0; i < reps.size(); ++i) {
        const double w = i < weights.size() ? weights[i] : 0.0;
        if (w <= 0.0) continue;
        const double lw = std::log(w);
        if (reps[i].lower != kNegInf) out.lower = log_add(out.lower, lw + reps[i].lower);
        if (reps[i].upper != kNegInf) out.upper = log_add(out.upper, lw + reps[i].upper);
    }
    return out;
}

EditClass edit_class_mass(const std::vector<FragmentState>& states,
                          std::uint32_t e1, std::uint32_t e2,
                          std::size_t m1_len, std::size_t m2_len,
                          const InsertPrior& ip, double log_eps, double log_1meps) {
    std::vector<FragmentState> sel;
    for (const FragmentState& st : states) {
        if (st.m1_edits == e1 && st.m2_edits == e2) sel.push_back(st);
    }
    EditClass out;
    out.count = sel.size();
    out.mass = fragment_states_mass(sel, m1_len, m2_len, ip, log_eps, log_1meps);
    return out;
}

TailLevel tail_interval_level(const std::string& r1, const std::string& r2,
                              const std::string& r1rc, const std::string& r2rc,
                              const std::string& hap, std::size_t d1, std::size_t d2,
                              std::size_t mult, const InsertPrior& ip,
                              double log_eps, double log_1meps, const PieceIndex* idx) {
    TailLevel out;
    const std::size_t D1 = d1 * mult, D2 = d2 * mult;
    const auto f1 = bounded_mate_placements(r1, hap, D1, nullptr, idx);
    const auto v1 = bounded_mate_placements(r1rc, hap, D1, nullptr, idx);
    const auto f2 = bounded_mate_placements(r2, hap, D2, nullptr, idx);
    const auto v2 = bounded_mate_placements(r2rc, hap, D2, nullptr, idx);
    const auto st = enumerate_fragment_states(0, f1, v1, f2, v2, r1.size(), r2.size(),
                                              ip.lo, ip.hi);
    const double m = fragment_states_mass(st, r1.size(), r2.size(), ip, log_eps, log_1meps);
    const double b = omitted_mass_bound(hap.size(), r1.size(), r2.size(), D1, D2, ip,
                                        log_eps, log_1meps, st);
    out.lower = m;
    out.upper = (m == kNegInf) ? b : (b == kNegInf ? m : log_add(m, b));
    out.bound = b;
    out.states = st.size();
    out.nonempty = !st.empty();
    // The PRODUCTION-band mass, recomputed at this depth from the same deeper state set: states
    // whose mates are both within the original d. It must not move as mult grows.
    std::vector<FragmentState> inband;
    for (const FragmentState& z : st) {
        if (z.m1_edits <= d1 && z.m2_edits <= d2) inband.push_back(z);
    }
    out.inband_at_d = fragment_states_mass(inband, r1.size(), r2.size(), ip, log_eps, log_1meps);
    return out;
}

TailInterval adaptive_tail_interval(const std::string& r1, const std::string& r2,
                                    const std::string& hap, std::size_t d1, std::size_t d2,
                                    const InsertPrior& ip, double log_eps, double log_1meps,
                                    double tolerance_nats, std::size_t max_depth_mult,
                                    const PieceIndex* idx) {
    TailInterval out;
    const std::string r1rc = reverse_complement(r1);
    const std::string r2rc = reverse_complement(r2);
    double prev_lower = kNegInf, prev_upper = std::numeric_limits<double>::infinity();
    double inband_ref = kNegInf;
    bool have_inband = false;
    for (std::size_t mult = 1; mult <= std::max<std::size_t>(1, max_depth_mult); ++mult) {
        const TailLevel lv = tail_interval_level(r1, r2, r1rc, r2rc, hap, d1, d2, mult, ip,
                                                 log_eps, log_1meps, idx);
        const double m = lv.lower, up = lv.upper;
        {
            const double mi = lv.inband_at_d;
            if (!have_inband) { inband_ref = mi; have_inband = true; }
            else if (!(mi == inband_ref ||
                       (mi != kNegInf && inband_ref != kNegInf &&
                        std::abs(mi - inband_ref) < 1e-9))) {
                out.inband_stable = false;
            }
            out.inband_at_d = inband_ref;
        }
        if (prev_lower != kNegInf && m != kNegInf && m < prev_lower - 1e-9) out.lower_monotone = false;
        if (prev_upper != std::numeric_limits<double>::infinity() && up != kNegInf &&
            up > prev_upper + 1e-9) out.upper_monotone = false;
        prev_lower = m; prev_upper = up;
        if (mult == 1) out.depth1_nonempty = lv.nonempty;
        out.lower = m; out.upper = up; out.bound = lv.bound;
        out.depth = mult; out.states = lv.states;
        const double width = (up == kNegInf || m == kNegInf) ? std::numeric_limits<double>::infinity()
                                                             : up - m;
        out.within_tolerance = width <= tolerance_nats;
        if (out.within_tolerance) break;
    }
    return out;
}

double fragment_states_mass(const std::vector<FragmentState>& states,
                            std::size_t m1_len, std::size_t m2_len,
                            const InsertPrior& ip, double log_eps, double log_1meps) {
    const auto read_ll = [&](std::uint32_t edits, std::size_t len) {
        return static_cast<double>(edits) * log_eps +
               static_cast<double>(len - edits) * log_1meps;
    };
    // log(0.5) PER STATE. Each state belongs to ONE of the two library orientations, and the
    // reference averages over them -- log_add(half + fwd, half + rev) -- rather than summing.
    // Omitting it makes this mass exactly log(2) = 0.6931 larger than the reference on every cell,
    // which is how the bound assertion caught it: an in-band sum cannot exceed the integral that
    // contains it, and it did, by that constant everywhere.
    const double log_half_strand = std::log(0.5);
    double acc = kNegInf;
    for (const FragmentState& st : states) {
        const double e1 = read_ll(st.m1_edits, m1_len);
        const double e2 = read_ll(st.m2_edits, m2_len);
        acc = log_add(acc, log_half_strand + e1 + e2 + ip.log_at(st.insert));
    }
    return acc;
}

std::vector<FragmentState> enumerate_fragment_states(
    std::uint32_t hap,
    const std::vector<MatePlacement>& m1_fwd, const std::vector<MatePlacement>& m1_rev,
    const std::vector<MatePlacement>& m2_fwd, const std::vector<MatePlacement>& m2_rev,
    std::size_t m1_len, std::size_t m2_len, long insert_lo, long insert_hi) {
    std::vector<FragmentState> out;
    // Orientation A: mate 1 forward, mate 2 reverse. Orientation B: mate 2 forward, mate 1 reverse.
    // Both are valid library orientations and BOTH must be formed; taking only one loses half the
    // states for every fragment whose mates happen to map the other way round.
    const auto build = [&](const std::vector<MatePlacement>& fwd, std::size_t fwd_len, bool fwd_is_m1,
                           const std::vector<MatePlacement>& rev, std::size_t rev_len) {
        for (const MatePlacement& f : fwd) {
            for (const MatePlacement& r : rev) {
                const long rev_end = r.start + static_cast<long>(rev_len) - 1;
                if (!valid_fr_coordinates(f.start, rev_end, insert_lo, insert_hi)) continue;
                FragmentState st;
                st.hap = hap;
                st.m1_start = fwd_is_m1 ? f.start : r.start;
                st.m1_fwd   = fwd_is_m1;
                st.m2_start = fwd_is_m1 ? r.start : f.start;
                st.m2_fwd   = !fwd_is_m1;
                st.m1_edits = fwd_is_m1 ? f.edits : r.edits;
                st.m2_edits = fwd_is_m1 ? r.edits : f.edits;
                st.frag_start = f.start;
                st.frag_end   = rev_end;
                st.insert     = rev_end - f.start + 1;
                out.push_back(st);
            }
        }
        (void)fwd_len;
    };
    build(m1_fwd, m1_len, true,  m2_rev, m2_len);
    build(m2_fwd, m2_len, false, m1_rev, m1_len);
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

std::vector<MatePlacement> exhaustive_mate_placements(const std::string& read,
                                                      const std::string& hap,
                                                      std::size_t max_edits) {
    std::vector<MatePlacement> out;
    if (read.empty() || hap.size() < read.size()) return out;
    for (std::size_t st = 0; st + read.size() <= hap.size(); ++st) {
        std::size_t mism = 0;
        for (std::size_t i = 0; i < read.size() && mism <= max_edits; ++i) {
            if (read[i] != hap[st + i]) ++mism;
        }
        if (mism <= max_edits) {
            out.push_back({static_cast<long>(st), static_cast<std::uint32_t>(mism)});
        }
    }
    return out;
}

namespace {
inline std::uint64_t encode_piece_at(const std::string& t, std::size_t at, std::size_t P) {
    std::uint64_t code = 0;
    for (std::size_t i = 0; i < P; ++i) {
        int b;
        switch (t[at + i]) {
            case 'A': case 'a': b = 0; break;
            case 'C': case 'c': b = 1; break;
            case 'G': case 'g': b = 2; break;
            case 'T': case 't': b = 3; break;
            default: return ~0ull;          // ambiguous: proposes nothing
        }
        code = (code << 2) | static_cast<std::uint64_t>(b);
    }
    return code;
}
}  // namespace

PieceIndex build_piece_index(const std::string& hap, std::size_t piece) {
    PieceIndex ix;
    if (piece == 0 || piece > 31 || hap.size() < piece) return ix;
    ix.piece = piece;
    ix.at.reserve(hap.size() * 2);
    for (std::size_t i = 0; i + piece <= hap.size(); ++i) {
        const std::uint64_t c = encode_piece_at(hap, i, piece);
        if (c == ~0ull) continue;
        ix.at[c].push_back(static_cast<std::uint32_t>(i));
    }
    return ix;
}

std::vector<MatePlacement> bounded_mate_placements(const std::string& read, const std::string& hap,
                                                   std::size_t max_edits, SearchWork* work,
                                                   const PieceIndex* idx) {
    SearchWork w;
    std::vector<MatePlacement> out;
    if (read.empty() || hap.size() < read.size()) { if (work) *work = w; return out; }
    if (max_edits >= read.size()) return bounded_mate_placements(read, hap, max_edits, work);
    const std::size_t npieces = max_edits + 1;
    const std::size_t P = read.size() / npieces;
    constexpr std::size_t kMinPiece = 12;
    const bool ambiguous = read.find_first_not_of("ACGTacgt") != std::string::npos;
    if (P < kMinPiece || ambiguous || idx == nullptr || !idx->usable() || idx->piece != P) {
        return bounded_mate_placements(read, hap, max_edits, work);   // same answer, slower path
    }
    w.pieces = npieces;
    std::vector<long> proposals;
    for (std::size_t pi = 0; pi < npieces; ++pi) {
        const std::size_t at = pi * P;
        const std::uint64_t c = encode_piece_at(read, at, P);
        if (c == ~0ull) continue;
        const auto it = idx->at.find(c);
        if (it == idx->at.end()) continue;
        for (const std::uint32_t pos : it->second) {
            const long st = static_cast<long>(pos) - static_cast<long>(at);
            if (st >= 0 && st + static_cast<long>(read.size()) <= static_cast<long>(hap.size())) {
                proposals.push_back(st);
                ++w.candidate_starts;
            }
        }
    }
    std::sort(proposals.begin(), proposals.end());
    proposals.erase(std::unique(proposals.begin(), proposals.end()), proposals.end());
    w.distinct_starts = proposals.size();
    for (const long st : proposals) {
        ++w.verified;
        std::size_t mism = 0;
        for (std::size_t i = 0; i < read.size() && mism <= max_edits; ++i) {
            if (read[i] != hap[static_cast<std::size_t>(st) + i]) ++mism;
        }
        if (mism <= max_edits) out.push_back({st, static_cast<std::uint32_t>(mism)});
    }
    std::sort(out.begin(), out.end());
    w.accepted = out.size();
    if (work) *work = w;
    return out;
}

std::vector<MatePlacement> bounded_mate_placements(const std::string& read, const std::string& hap,
                                                   std::size_t max_edits, SearchWork* work) {
    SearchWork w;
    std::vector<MatePlacement> out;
    if (read.empty() || hap.size() < read.size()) {
        if (work) *work = w;
        return out;
    }
    // TERMINAL DEPTH, PER MATE. Once max_edits reaches THIS read's length every position is
    // trivially in band and d+1 non-empty pieces cannot be constructed. Short-circuit to the
    // exhaustive scan. Applied per mate, so a short mate can be exhaustive while its longer partner
    // still uses pigeonhole recruitment.
    if (max_edits >= read.size()) {
        SearchWork w2;
        w2.pieces = 0;
        w2.exhaustive_fallback = true;
        auto all = exhaustive_mate_placements(read, hap, max_edits);
        w2.candidate_starts = w2.distinct_starts = w2.verified =
            hap.size() >= read.size() ? hap.size() - read.size() + 1 : 0;
        w2.accepted = all.size();
        if (work) *work = w2;
        return all;
    }
    const std::size_t npieces = max_edits + 1;
    const std::size_t P = read.size() / npieces;
    w.pieces = npieces;
    // Below this the pieces propose essentially every position and filter nothing; scanning is both
    // simpler and faster. The two paths must return the SAME placements -- this is an acceleration
    // of the exhaustive scan, never a different answer.
    constexpr std::size_t kMinPiece = 12;
    // AMBIGUITY FORCES THE EXHAUSTIVE PATH. A piece containing N proposes nothing, but an N in the
    // read that meets an N in the haplotype costs NO Hamming mismatch -- so a read can be well
    // inside the band while every one of its d+1 pieces is spoiled by ambiguity, and the pigeonhole
    // then proposes nothing at all and the placement is lost. Skipping such pieces without falling
    // back, which is what the first version did, is silently lossy exactly where reads are worst.
    const bool ambiguous = read.find_first_not_of("ACGTacgt") != std::string::npos;
    if (P < kMinPiece || ambiguous) {
        w.exhaustive_fallback = true;
        out = exhaustive_mate_placements(read, hap, max_edits);
        w.candidate_starts = w.distinct_starts = w.verified =
            hap.size() >= read.size() ? hap.size() - read.size() + 1 : 0;
        w.accepted = out.size();
        if (work) *work = w;
        return out;
    }
    std::vector<long> proposals;
    for (std::size_t pi = 0; pi < npieces; ++pi) {
        const std::size_t at = pi * P;
        // An ambiguous base makes this piece propose nothing; the pigeonhole still holds because
        // some OTHER piece must match exactly for an in-band placement -- unless the mismatches are
        // spread so that every piece is spoiled, which only ambiguity in the read can cause. That is
        // why the caller must treat N-containing reads as exhaustive.
        if (read.find_first_of("Nn", at) < at + P) continue;
        const std::string piece = read.substr(at, P);
        std::size_t pos = hap.find(piece, 0);
        while (pos != std::string::npos) {
            const long st = static_cast<long>(pos) - static_cast<long>(at);
            if (st >= 0 && st + static_cast<long>(read.size()) <= static_cast<long>(hap.size())) {
                proposals.push_back(st);
                ++w.candidate_starts;
            }
            pos = hap.find(piece, pos + 1);
        }
    }
    std::sort(proposals.begin(), proposals.end());
    proposals.erase(std::unique(proposals.begin(), proposals.end()), proposals.end());
    w.distinct_starts = proposals.size();
    for (const long st : proposals) {
        ++w.verified;
        std::size_t mism = 0;
        for (std::size_t i = 0; i < read.size() && mism <= max_edits; ++i) {
            if (read[i] != hap[static_cast<std::size_t>(st) + i]) ++mism;
        }
        if (mism <= max_edits) {
            out.push_back({st, static_cast<std::uint32_t>(mism)});
        }
    }
    std::sort(out.begin(), out.end());
    w.accepted = out.size();
    if (work) *work = w;
    return out;
}

PathProjection project_path_blocks(const std::vector<BlockAlleles>& projection_blocks,
                                   const std::vector<BlockAlleles>& catalogue_blocks,
                                   const std::string& name,
                                   const std::string& walk) {
    PathProjection out;
    const CandidateFrame f = build_candidate_frame(projection_blocks, name, walk);
    if (!f.ok) return out;                       // unprojectable: every block field is NA
    out.ok = true;
    out.partial = f.partial;
    out.reverse_frame = f.reverse_frame;
    if (f.mapped_lo > 0) out.unmapped.emplace_back(0, f.mapped_lo);
    if (f.mapped_hi < walk.size()) out.unmapped.emplace_back(f.mapped_hi, walk.size());

    for (std::size_t k = 0; k < f.offsets.size(); ++k) {
        const std::size_t b0 = f.offsets[k];
        const std::size_t b1 = (k + 1 < f.offsets.size()) ? f.offsets[k + 1] : f.mapped_hi;
        if (b1 < b0 || b0 < f.mapped_lo || b1 > f.mapped_hi) continue;
        // A block this path does not OCCUPY is not a block it spells emptily. Both look like a
        // zero-width slice: a bypassing allele (a deletion spanning the site) is a real, projectable
        // state with an entry in allele_of, while a block past the end of a truncated path has no
        // entry at all and lands at the boundary offset with zero width. Emitting the second marks
        // it projectable and lets a truncated path claim block coordinates it never reached --
        // measured on hapTRUNC, where the unprojectable flag then never fired anywhere.
        if (sl_block_has_no_allele(projection_blocks, k < f.block_at.size() ? f.block_at[k]
                                                                           : static_cast<std::uint32_t>(k),
                                   name)) continue;
        PathBlockSlice sl;
        sl.block = k < f.block_at.size() ? f.block_at[k] : static_cast<std::uint32_t>(k);
        sl.walk_begin = b0;
        sl.walk_end = b1;
        sl.reverse = f.reverse_frame;
        const std::string raw = walk.substr(b0, b1 - b0);
        // Emitted in REFERENCE/BLOCK orientation. The walk slice of an antiparallel path is the
        // reverse complement of the block's own sequence, and emitting walk bytes would make the
        // same allele hash differently depending on which strand the assembly happened to be on.
        sl.seq = f.reverse_frame ? reverse_complement(raw) : raw;
        sl.md5 = md5_hex(sl.seq);
        const std::string rc = reverse_complement(sl.seq);
        sl.canonical_md5 = md5_hex(sl.seq < rc ? sl.seq : rc);
        // Representability is decided against the CALLING catalogue by exact reference-oriented
        // sequence equality -- not by looking the path up in allele_of, which for a held-out path
        // finds it in its own held-out catalogue and answers a question nobody asked.
        if (sl.block < catalogue_blocks.size()) {
            const auto& cat = catalogue_blocks[sl.block].allele_seq;
            for (std::size_t ai = 0; ai < cat.size(); ++ai) {
                if (cat[ai] == sl.seq) {
                    sl.catalogue_allele = static_cast<long>(ai);
                    sl.catalogue_representable = true;
                    break;
                }
            }
        }
        out.blocks.push_back(std::move(sl));
    }
    std::sort(out.blocks.begin(), out.blocks.end(),
              [](const PathBlockSlice& a, const PathBlockSlice& b) { return a.block < b.block; });
    return out;
}

CandidateFrame build_candidate_frame(const std::vector<BlockAlleles>& blocks,
                                     const std::string& name,
                                     const std::string& walk) {
    CandidateFrame f;
    f.seq = walk;
    std::string concat;
    std::vector<std::size_t> off;
    for (std::size_t bi = 0; bi < blocks.size(); ++bi) {
        off.push_back(concat.size());
        const auto it = blocks[bi].allele_of.find(name);
        if (it != blocks[bi].allele_of.end() && it->second < blocks[bi].allele_seq.size()) {
            concat += blocks[bi].allele_seq[it->second];
        }
    }
    // Where `concat` sits inside `target`: 0 if it is a prefix, target.size()-concat.size() if a
    // suffix, npos otherwise. A truncated assembly loses sequence from an END, so those are the two
    // placements that can arise; anything else is a genuine disagreement and stays refused.
    const auto placement_in = [](const std::string& c, const std::string& target) -> std::size_t {
        if (c.empty() || c.size() > target.size()) return std::string::npos;
        if (target.compare(0, c.size(), c) == 0) return 0;
        if (target.compare(target.size() - c.size(), c.size(), c) == 0) {
            return target.size() - c.size();
        }
        return std::string::npos;
    };

    const std::size_t fwd_at = placement_in(concat, walk);
    if (fwd_at != std::string::npos) {
        f.offsets.resize(off.size());
        for (std::size_t k = 0; k < off.size(); ++k) f.offsets[k] = off[k] + fwd_at;
        f.block_at.resize(f.offsets.size());
        for (std::uint32_t k = 0; k < f.block_at.size(); ++k) f.block_at[k] = k;
        f.mapped_lo = fwd_at;
        f.mapped_hi = fwd_at + concat.size();
        f.partial = concat.size() != walk.size();
        f.ok = true;
        return f;
    }
    const std::string rcw = reverse_complement(walk);
    const std::size_t rev_at = placement_in(concat, rcw);
    if (rev_at != std::string::npos) {
        // Opposite frames. Block b occupies [off[b], off[b+1]) in the CONCAT frame, which maps to
        // [n - off[b+1], n - off[b]) in the walk frame; so the walk-frame start offsets are the
        // mirrored ends, in reverse block order.
        const std::size_t n = walk.size();
        std::vector<std::size_t> mirrored(off.size(), 0);
        std::vector<std::uint32_t> at(off.size(), 0);
        for (std::size_t bi = 0; bi < off.size(); ++bi) {
            const std::size_t end = (bi + 1 < off.size()) ? off[bi + 1] : concat.size();
            const std::size_t k = off.size() - 1 - bi;   // walk-order slot for block bi
            mirrored[k] = n - (rev_at + end);
            at[k] = static_cast<std::uint32_t>(bi);      // ...and it is still block bi
        }
        f.offsets = std::move(mirrored);
        f.block_at = std::move(at);
        f.reverse_frame = true;
        f.mapped_lo = n - (rev_at + concat.size());
        f.mapped_hi = n - rev_at;
        f.partial = concat.size() != walk.size();
        f.ok = true;
        return f;
    }
    // No trustworthy map. Reported rather than guessed: a candidate whose block boundaries are
    // unknown cannot contribute a scope, and pretending otherwise would put arbitrary blocks into
    // some fragment's dependency set.
    f.ok = false;
    return f;
}

std::string chain_oriented_sequence(const CandidateFrame& frame) {
    // An antiparallel candidate's walk bytes ARE the reverse complement of its chain-oriented
    // sequence; returning them raw makes a correct frame look like a different haplotype.
    return frame.reverse_frame ? reverse_complement(frame.seq) : frame.seq;
}

bool chain_oriented_block_span(const CandidateFrame& frame, std::uint32_t block,
                               std::size_t& lo, std::size_t& hi) {
    bool seen = false;
    std::size_t wlo = 0, whi = 0;
    for (std::size_t k = 0; k < frame.block_at.size(); ++k) {
        if (frame.block_at[k] != block) continue;
        const std::size_t s0 = frame.offsets[k];
        const std::size_t s1 = (k + 1 < frame.offsets.size()) ? frame.offsets[k + 1]
                                                              : frame.seq.size();
        if (!seen) { wlo = s0; whi = s1; seen = true; }
        else { wlo = std::min(wlo, s0); whi = std::max(whi, s1); }
    }
    if (!seen) return false;
    if (!frame.reverse_frame) { lo = wlo; hi = whi; return true; }
    // Mirror the interval: reversing the sequence maps [wlo, whi) to [n - whi, n - wlo).
    const std::size_t n = frame.seq.size();
    lo = n - whi;
    hi = n - wlo;
    return true;
}

FrameCoverage assess_frame_coverage(const Graph& graph,
                                    const std::vector<BlockAlleles>& blocks,
                                    const std::vector<std::string>& hmm_states) {
    FrameCoverage C;
    const auto by_name = path_records_by_name(graph);
    for (const std::string& nm : hmm_states) {
        const auto it = by_name.find(nm);
        if (it == by_name.end() || it->second == nullptr) {
            C.missing_names.push_back(nm);
            C.missing_reasons.push_back("absent-from-graph");
            continue;
        }
        bool spelled = false;
        const std::string walk = spell_path_steps_sequence(graph, it->second->steps, &spelled);
        if (!spelled) {
            C.missing_names.push_back(nm);
            C.missing_reasons.push_back("walk-not-spellable");
            continue;
        }
        CandidateFrame cf = build_candidate_frame(blocks, nm, walk);
        if (!cf.ok) {
            // The ONLY failure this function has: no trustworthy map at all. A partial frame is a
            // SUCCESS with an unmapped remainder, and never arrives here.
            C.missing_names.push_back(nm);
            C.missing_reasons.push_back("frame-unverified");
            continue;
        }
        (cf.partial ? C.partial_names : C.complete_names).push_back(nm);
        C.framed_names.push_back(nm);
        C.frames.push_back(std::move(cf));
    }
    // UNIQUENESS FIRST. Comparing sorted vectors alone would let a duplicated state name appear on
    // both sides and pass, so neither side may contain a repeat.
    const auto unique = [](std::vector<std::string> v) {
        std::sort(v.begin(), v.end());
        return std::adjacent_find(v.begin(), v.end()) == v.end();
    };
    C.names_unique = unique(hmm_states) && unique(C.framed_names);
    std::vector<std::string> a = hmm_states, b = C.framed_names;
    std::sort(a.begin(), a.end());
    std::sort(b.begin(), b.end());
    C.all_states_usable = C.names_unique && C.missing_names.empty() && a == b;
    return C;
}

double scope_restricted_pair_loglik(const CandidateFrame& frame_a, const CandidateFrame& frame_b,
                                    const std::vector<Fragment>& fragments,
                                    const std::vector<std::vector<std::uint32_t>>& scopes,
                                    const ReferenceParams& params,
                                    bool include_unmapped) {
    const std::string& hap_a = frame_a.seq;
    const std::string& hap_b = frame_b.seq;
    const double log_eps = std::log(params.error_rate / 3.0);
    const double log_1meps = std::log1p(-params.error_rate);
    const long min_len = fragment_insert_floor(fragments, params.allow_overlapping_pairs);
    const InsertPrior ip = make_insert_prior(params.fragment_len, params.fragment_sd,
                                             params.discordant_rate, params.insert_sigmas, min_len);
    // EXPOSURE ONCE, over the whole locus, exactly as the reference charges it. Charging it per
    // factor is what double-counted the latent start space in the recruitment version.
    const auto exposure_of = [&](const std::string& h) {
        const long n = static_cast<long>(h.size());
        double e = 0.0;
        for (long L = ip.lo; L <= ip.hi; ++L) {
            const double starts = static_cast<double>(std::max<long>(0, n - L + 1));
            if (starts > 0.0) e += std::exp(ip.log_at(L)) * starts;
        }
        return e;
    };
    double total = -params.lambda * (exposure_of(hap_a) + exposure_of(hap_b));
    const double log_lam = std::log(params.lambda);
    const double log_mix = std::log1p(-params.eta);
    const double log_bg_w = std::log(params.eta);

    for (std::size_t fi = 0; fi < fragments.size(); ++fi) {
        const Fragment& f = fragments[fi];
        const std::size_t len = f.bases();
        if (len == 0) continue;
        const std::vector<std::uint32_t>& sc =
            fi < scopes.size() ? scopes[fi] : std::vector<std::uint32_t>{};
        const auto in_scope = [&](std::uint32_t lo, std::uint32_t hi) {
            // Unattributable mass is ALWAYS included. It belongs to no block, so no scope can drop
            // it -- dropping it would make the restricted model lose likelihood that the whole-locus
            // model has, and the reconciliation residual would stop being zero.
            if (lo == kUnmappedBlock) return include_unmapped;
            for (std::uint32_t b = lo; b <= hi; ++b) {
                if (!std::binary_search(sc.begin(), sc.end(), b)) return false;
            }
            return true;
        };
        const std::string r2rc = f.r2.empty() ? std::string() : reverse_complement(f.r2);
        const std::string r1rc = reverse_complement(f.r1);
        const double half = std::log(0.5);
        double lse = kNegInf;
        for (int side = 0; side < 2; ++side) {
            const CandidateFrame& fr = side == 0 ? frame_a : frame_b;
            const std::string& hap = fr.seq;
            if (hap.empty()) continue;
            const long n = static_cast<long>(hap.size());
            for (int orient = 0; orient < 2; ++orient) {
                const std::string& lead = orient == 0 ? f.r1 : f.r2;
                const std::string& trail = orient == 0 ? r2rc : r1rc;
                if (lead.empty()) continue;
                for (long st = 0; st + static_cast<long>(lead.size()) <= n; ++st) {
                    const double e1 = reference_emission(lead, hap, static_cast<std::size_t>(st),
                                                         log_eps, log_1meps);
                    if (e1 == kNegInf) continue;
                    if (trail.empty()) {
                        const auto sp = ordered_block_span(fr, st,
                                            st + static_cast<long>(lead.size()) - 1);
                        if (in_scope(sp.first, sp.second)) {
                            lse = log_add(lse, half + e1);
                        }
                        continue;
                    }
                    for (long L = ip.lo; L <= ip.hi; ++L) {
                        if (st + L > n) break;
                        const long m2 = st + L - static_cast<long>(trail.size());
                        if (m2 < 0) continue;
                        const double e2 = reference_emission(trail, hap, static_cast<std::size_t>(m2),
                                                             log_eps, log_1meps);
                        if (e2 == kNegInf) continue;
                        const auto sq = ordered_block_span(fr, st, st + L - 1);
                        if (!in_scope(sq.first, sq.second)) continue;
                        lse = log_add(lse, half + e1 + e2 + ip.log_at(L));
                    }
                }
            }
        }
        const std::size_t bg_edits =
            static_cast<std::size_t>(params.bg_divergence * static_cast<double>(len));
        const double bg = static_cast<double>(bg_edits) * log_eps +
                          static_cast<double>(len - bg_edits) * log_1meps;
        total += log_add(lse == kNegInf ? kNegInf : log_mix + log_lam + lse, log_bg_w + bg);
    }
    return total;
}

std::pair<std::uint32_t, std::uint32_t> ordered_block_span(const CandidateFrame& frame,
                                                          long start, long end) {
    const auto at = [&](long pos) -> std::uint32_t {
        if (frame.offsets.empty()) return 0;
        // Outside the verified window there is no block, and saying so is the whole point: a
        // guessed block index here would put an arbitrary block into some fragment's dependency set.
        if (pos < static_cast<long>(frame.mapped_lo) ||
            pos >= static_cast<long>(frame.mapped_hi)) return kUnmappedBlock;
        std::size_t k = 0;
        while (k + 1 < frame.offsets.size() &&
               static_cast<long>(frame.offsets[k + 1]) <= pos) ++k;
        return k < frame.block_at.size() ? frame.block_at[k] : static_cast<std::uint32_t>(k);
    };
    const std::uint32_t a = at(start), b = at(end);
    // An origin touching unmapped sequence is unattributable, not partly attributable: reporting
    // the mapped endpoint alone would claim the whole origin depends on that block.
    if (a == kUnmappedBlock || b == kUnmappedBlock) return {kUnmappedBlock, kUnmappedBlock};
    return {std::min(a, b), std::max(a, b)};
}

OriginUniverse enumerate_fragment_origins(const Fragment& fragment,
                                          const std::vector<CandidateFrame>& frames,
                                          const ReferenceParams& params,
                                          std::size_t retain_topk,
                                          double scope_tol) {
    OriginUniverse out;
    if (fragment.r1.empty()) return out;
    const double log_eps = std::log(params.error_rate / 3.0);
    const double log_1meps = std::log1p(-params.error_rate);
    const long min_len = fragment_insert_floor(fragment, params.allow_overlapping_pairs);
    const InsertPrior ip = make_insert_prior(params.fragment_len, params.fragment_sd,
                                             params.discordant_rate, params.insert_sigmas,
                                             std::max<long>(1, min_len));
    const std::string r2rc = fragment.r2.empty() ? std::string()
                                                 : reverse_complement(fragment.r2);
    const std::string r1rc = reverse_complement(fragment.r1);
    const double half = std::log(0.5);

    // EXACTLY the states the exposure counts: (start, insert length, orientation). Enumerating a
    // different set here would make the oracle's sum incomparable with the reference's, which is the
    // one thing it exists to be compared against.
    for (std::uint32_t h = 0; h < frames.size(); ++h) {
        const std::string& hap = frames[h].seq;
        if (hap.empty() || !frames[h].ok) continue;
        const long n = static_cast<long>(hap.size());
        for (int orient = 0; orient < 2; ++orient) {
            const std::string& lead = orient == 0 ? fragment.r1 : fragment.r2;
            const std::string& trail = orient == 0 ? r2rc : r1rc;
            if (lead.empty()) continue;
            for (long st = 0; st + static_cast<long>(lead.size()) <= n; ++st) {
                const double e1 = reference_emission(lead, hap, static_cast<std::size_t>(st),
                                                     log_eps, log_1meps);
                if (e1 == kNegInf) continue;
                if (trail.empty()) {
                    FragmentOrigin o;
                    o.hap = h; o.start = st; o.insert = static_cast<std::int32_t>(lead.size());
                    o.fwd = orient == 0; o.emission_ll = half + e1; o.insert_ll = 0.0;
                    // ORDERED. For an antiparallel candidate block_of DECREASES along walk
                    // coordinates, so start/end come back reversed. A reversed interval makes the
                    // touched-block loop visit nothing and makes in_scope() accept the origin
                    // vacuously -- the origin would be silently exempt from every scope test.
                    {
                        const auto sp = ordered_block_span(frames[h], st,
                                            st + static_cast<long>(lead.size()) - 1);
                        o.block_lo = sp.first; o.block_hi = sp.second;
                    }
                    out.origins.push_back(o);
                    continue;
                }
                for (long L = ip.lo; L <= ip.hi; ++L) {
                    if (st + L > n) break;
                    const long m2 = st + L - static_cast<long>(trail.size());
                    if (m2 < 0) continue;
                    const double e2 = reference_emission(trail, hap, static_cast<std::size_t>(m2),
                                                         log_eps, log_1meps);
                    if (e2 == kNegInf) continue;
                    FragmentOrigin o;
                    o.hap = h; o.start = st; o.insert = static_cast<std::int32_t>(L);
                    o.fwd = orient == 0;
                    o.emission_ll = half + e1 + e2;
                    o.insert_ll = ip.log_at(L);
                    {
                        const auto sp = ordered_block_span(frames[h], st, st + L - 1);
                        o.block_lo = sp.first; o.block_hi = sp.second;
                    }
                    out.origins.push_back(o);
                }
            }
        }
    }

    double lse = kNegInf;
    for (const FragmentOrigin& o : out.origins) lse = log_add(lse, o.emission_ll + o.insert_ll);
    out.exact_lse = lse;

    // SCOPE, DECIDED PER CANDIDATE AND THEN UNIONED.
    //
    // The counterfactual must ask "can this block's allele change the likelihood", and that question
    // is only meaningful inside ONE candidate. Deleting a block's origins from the aggregate over
    // every panel candidate measures placement mass located there, which is a different thing, and
    // it fails in three ways: a block that matters to one low-mass candidate is drowned by unrelated
    // candidates carrying more mass; duplicating a panel path changes the answer; and the scope then
    // depends on which haplotypes the panel happens to contain.
    //
    // So: for each candidate separately, delete the origins touching block b and ask whether THAT
    // candidate's own logsumexp moves by more than scope_tol. Union over candidates. Mass is never
    // averaged across candidates before the decision.
    //
    // APPROXIMATION, AND ITS BOUND. The emission is finite at every position, so the exact structural
    // union is locus-wide for every fragment. This is therefore a THRESHOLDED scope, not the exact
    // minimal one, and the guarantee it carries is per candidate: for each candidate, the mass
    // dropped by restricting to the reported scope is at most scope_tol nats of that candidate's own
    // total. It is not a statement about the aggregate.
    // Recorded separately, before any scope decision, so the caller can tell a candidate whose
    // likelihood is fully block-attributable from one carrying mass no block can express.
    out.unmapped_lse = kNegInf;
    for (const FragmentOrigin& o : out.origins) {
        if (o.block_lo == kUnmappedBlock) {
            out.unmapped_lse = log_add(out.unmapped_lse, o.emission_ll + o.insert_ll);
        }
    }
    std::vector<std::uint32_t> touched;
    for (const FragmentOrigin& o : out.origins) {
        if (o.block_lo == kUnmappedBlock) continue;   // no block to touch
        for (std::uint32_t b = o.block_lo; b <= o.block_hi; ++b) touched.push_back(b);
    }
    std::sort(touched.begin(), touched.end());
    touched.erase(std::unique(touched.begin(), touched.end()), touched.end());

    std::vector<double> per_cand_lse(frames.size(), kNegInf);
    for (const FragmentOrigin& o : out.origins) {
        per_cand_lse[o.hap] = log_add(per_cand_lse[o.hap], o.emission_ll + o.insert_ll);
    }
    // The threshold applies to the FULL per-fragment contribution, not to the raw placement mass:
    //
    //     contrib(M) = log[ (1-eta) * lambda * M  +  eta * P_bg ]
    //
    // For a candidate the fragment does not belong to, the placement mass sits far below the
    // background and the term is background-dominated, so deleting ANY block moves the raw mass a
    // lot and the contribution not at all. Thresholding the raw mass therefore puts every block in
    // scope via candidates the fragment never came from -- measured: it put an untouched middle
    // block into a shared-sequence fragment's scope purely through a candidate that lacks the
    // shared unit entirely.
    const std::size_t flen = fragment.bases();
    const std::size_t bg_edits =
        static_cast<std::size_t>(params.bg_divergence * static_cast<double>(flen));
    const double bg_ll = static_cast<double>(bg_edits) * log_eps +
                         static_cast<double>(flen - bg_edits) * log_1meps;
    const double log_mix = std::log1p(-params.eta);
    const double log_bg_w = std::log(params.eta);
    const double log_lam = std::log(params.lambda);
    const auto contrib_of = [&](double mass) {
        return log_add(mass == kNegInf ? kNegInf : log_mix + log_lam + mass, log_bg_w + bg_ll);
    };

    out.panel_domain_scope.clear();
    for (const std::uint32_t b : touched) {
        bool matters = false;
        for (std::uint32_t h = 0; h < frames.size() && !matters; ++h) {
            if (per_cand_lse[h] == kNegInf) continue;
            double without = kNegInf;
            for (const FragmentOrigin& o : out.origins) {
                if (o.hap != h) continue;
                if (b >= o.block_lo && b <= o.block_hi) continue;
                without = log_add(without, o.emission_ll + o.insert_ll);
            }
            if (contrib_of(per_cand_lse[h]) - contrib_of(without) > scope_tol) matters = true;
        }
        if (matters) out.panel_domain_scope.push_back(b);
    }

    // JOINT BOUND. The per-block test above is necessary and not sufficient: it asks what removing
    // ONE block costs, while the restricted model removes them ALL. Ten blocks each under scope_tol
    // can sum to far more. So recompute each candidate's contribution using the whole excluded set
    // at once, and if the bound fails, add excluded blocks back -- greedily, most omitted mass first
    // -- until it holds. The achieved residual is reported rather than assumed.
    {
        const auto restricted_mass = [&](std::uint32_t h,
                                         const std::vector<std::uint32_t>& sc) {
            double m = kNegInf;
            for (const FragmentOrigin& o : out.origins) {
                if (o.hap != h) continue;
                bool inside = true;
                if (o.block_lo != kUnmappedBlock) {
                    for (std::uint32_t b = o.block_lo; b <= o.block_hi && inside; ++b) {
                        if (!std::binary_search(sc.begin(), sc.end(), b)) inside = false;
                    }
                }
                if (inside) m = log_add(m, o.emission_ll + o.insert_ll);
            }
            return m;
        };
        // DIPLOID residual. `contrib_of` is the haploid mixture, but no genotype is scored that
        // way: scope_restricted_pair_loglik sums BOTH haplotypes' placement mass into a single
        // log_add before mixing with the background. A per-candidate guarantee does not imply the
        // pair guarantee -- two candidates each losing a little mass lose it into the same sum,
        // and a pair whose partner contributes almost nothing amplifies the survivor's loss. So
        // the bound is taken over every unordered pair, INCLUDING a == b: the homozygote is a real
        // genotype and the scorer visits its frame twice, doubling both the full and the
        // restricted mass, which is not the same residual as the heterozygote's.
        const auto contrib_pair = [&](double x, double y) {
            const double s = log_add(x, y);
            return log_add(s == kNegInf ? kNegInf : log_mix + log_lam + s, log_bg_w + bg_ll);
        };
        const auto worst_residual = [&](const std::vector<std::uint32_t>& sc,
                                        std::uint32_t* wa, std::uint32_t* wb) {
            std::vector<double> restr(frames.size(), kNegInf);
            for (std::uint32_t h = 0; h < frames.size(); ++h) {
                if (per_cand_lse[h] != kNegInf) restr[h] = restricted_mass(h, sc);
            }
            double worst = 0.0;
            for (std::uint32_t a = 0; a < frames.size(); ++a) {
                for (std::uint32_t b = a; b < frames.size(); ++b) {
                    const double d = std::abs(contrib_pair(per_cand_lse[a], per_cand_lse[b]) -
                                              contrib_pair(restr[a], restr[b]));
                    if (d > worst) { worst = d; if (wa) *wa = a; if (wb) *wb = b; }
                }
            }
            return worst;
        };
        std::vector<std::uint32_t> excluded;
        for (const std::uint32_t b : touched) {
            if (!std::binary_search(out.panel_domain_scope.begin(), out.panel_domain_scope.end(), b)) excluded.push_back(b);
        }
        double worst = worst_residual(out.panel_domain_scope, &out.worst_pair_a, &out.worst_pair_b);
        out.initial_bound = worst;
        while (worst > scope_tol && !excluded.empty()) {
            // add back whichever excluded block carries the most omitted mass
            std::size_t best_i = 0;
            double best_mass = kNegInf;
            for (std::size_t i = 0; i < excluded.size(); ++i) {
                double m = kNegInf;
                for (const FragmentOrigin& o : out.origins) {
                    if (excluded[i] >= o.block_lo && excluded[i] <= o.block_hi) {
                        m = log_add(m, o.emission_ll + o.insert_ll);
                    }
                }
                if (m > best_mass) { best_mass = m; best_i = i; }
            }
            out.panel_domain_scope.push_back(excluded[best_i]);
            std::sort(out.panel_domain_scope.begin(), out.panel_domain_scope.end());
            excluded.erase(excluded.begin() + static_cast<long>(best_i));
            ++out.blocks_added;
            worst = worst_residual(out.panel_domain_scope, &out.worst_pair_a, &out.worst_pair_b);
        }
        out.achieved_bound = worst;
        out.bound_holds = worst <= scope_tol;
    }

    // What an accelerated representation keeping only the best `retain_topk` origins per candidate
    // would retain, and what it would drop. Reported as MASS, because a count says nothing about
    // whether the discarded origins mattered.
    if (retain_topk == 0) {
        out.retained_lse = out.exact_lse;
        out.omitted_lse = kNegInf;
    } else {
        std::vector<std::vector<double>> per_hap(frames.size());
        for (const FragmentOrigin& o : out.origins) {
            per_hap[o.hap].push_back(o.emission_ll + o.insert_ll);
        }
        double keep = kNegInf, drop = kNegInf;
        for (auto& v : per_hap) {
            std::sort(v.begin(), v.end(), std::greater<double>());
            for (std::size_t i = 0; i < v.size(); ++i) {
                if (i < retain_topk) keep = log_add(keep, v[i]);
                else drop = log_add(drop, v[i]);
            }
        }
        out.retained_lse = keep;
        out.omitted_lse = drop;
    }
    return out;
}


} // namespace panvar

namespace panvar {

ExposureCheck check_exposure(std::size_t window_len, const InsertPrior& ip) {
    ExposureCheck out;
    const double n = static_cast<double>(window_len);
    double mean_len = 0.0;
    for (long L = ip.lo; L <= ip.hi; ++L) {
        const double lp = ip.log_at(L);
        if (lp == kNegInf) continue;
        const double p = std::exp(lp);
        mean_len += p * static_cast<double>(L);
        out.exact += p * std::max(0.0, n - static_cast<double>(L) + 1.0);
    }
    out.affine = n + 1.0 - mean_len;
    // The clip max(0, n - L + 1) only bites once n < L - 1 for some L in support, so the two forms
    // coincide from hi - 1 upwards. Below that the affine surrogate is simply a different function
    // and the cis/trans cancellation it justifies does not hold.
    out.in_regime = static_cast<long>(window_len) >= ip.hi - 1;
    return out;
}

const char* owner_kind_name(OwnerKind k) {
    switch (k) {
        case OwnerKind::Unary:    return "unary";
        case OwnerKind::Linkage:  return "linkage";
        case OwnerKind::Wide:     return "wide";
        case OwnerKind::Invariant: return "invariant";
        default:                  return "unusable";
    }
}

FragmentOwner assign_fragment_owner_panel_domain(const Fragment& fragment,
                                                 const std::vector<CandidateFrame>& frames,
                                                 const std::vector<char>& block_variable,
                                                 const InsertPrior& ip, double max_divergence,
                                                 double log_eps, double log_1meps,
                                                 double scope_tol,
                                                 const std::vector<PieceIndex>* pidx) {
    FragmentOwner out;
    out.in_band = kNegInf;
    out.omitted_bound = kNegInf;
    out.unmapped = kNegInf;
    if (fragment.r1.empty() || fragment.r2.empty()) {
        out.why = UnusableReason::EmptyRead;   // set HERE: this return is the reachable one
        return out;
    }

    const std::size_t d1 = mate_band_edits(max_divergence, fragment.r1.size());
    const std::size_t d2 = mate_band_edits(max_divergence, fragment.r2.size());
    const std::string a1 = reverse_complement(fragment.r1);
    const std::string a2 = reverse_complement(fragment.r2);
    const double log_half_strand = std::log(0.5);

    // One origin per in-band state, carrying its own mass and the block span it touches. Collected
    // over EVERY candidate before anything is decided: the scope must not depend on which pair is
    // being scored.
    struct Origin { double mass; std::uint32_t lo, hi; };
    std::vector<Origin> origins;
    for (std::uint32_t h = 0; h < frames.size(); ++h) {
        if (!frames[h].ok || frames[h].seq.empty()) continue;
        const PieceIndex* ix = (pidx && pidx->size() == frames.size()) ? &(*pidx)[h] : nullptr;
        const std::string& hap = frames[h].seq;
        const auto f1 = bounded_mate_placements(fragment.r1, hap, d1, nullptr, ix);
        const auto v1 = bounded_mate_placements(a1, hap, d1, nullptr, ix);
        const auto f2 = bounded_mate_placements(fragment.r2, hap, d2, nullptr, ix);
        const auto v2 = bounded_mate_placements(a2, hap, d2, nullptr, ix);
        const auto st = enumerate_fragment_states(h, f1, v1, f2, v2, fragment.r1.size(),
                                                  fragment.r2.size(), ip.lo, ip.hi);
        // The band's own omitted-mass bound, per candidate, so what lies OUTSIDE the band is
        // bounded rather than assumed away. Summed across candidates: dropping it would make the
        // scope look certified when the evidence for it was never examined.
        const double b = omitted_mass_bound(hap.size(), fragment.r1.size(), fragment.r2.size(),
                                            d1, d2, ip, log_eps, log_1meps, st);
        if (b != kNegInf) out.omitted_bound = log_add(out.omitted_bound, b);
        for (const FragmentState& z : st) {
            const double e1 = static_cast<double>(z.m1_edits) * log_eps +
                              static_cast<double>(fragment.r1.size() - z.m1_edits) * log_1meps;
            const double e2 = static_cast<double>(z.m2_edits) * log_eps +
                              static_cast<double>(fragment.r2.size() - z.m2_edits) * log_1meps;
            const double m = log_half_strand + e1 + e2 + ip.log_at(z.insert);
            // frag_start/frag_end, NOT start+insert: the state already carries its ordered span,
            // and re-deriving it would get an antiparallel candidate's interval backwards -- the
            // defect that once turned a fragment's scope from {1} into {1,2}.
            const auto span = ordered_block_span(frames[h], z.frag_start, z.frag_end);
            out.in_band = log_add(out.in_band, m);
            if (span.first == kUnmappedBlock || span.second == kUnmappedBlock) {
                // Unattributable: it belongs to no block, so it can never be dropped by narrowing
                // the scope and can never justify widening it. Kept separate, never folded in.
                out.unmapped = log_add(out.unmapped, m);
                continue;
            }
            origins.push_back(Origin{m, span.first, span.second});
        }
    }
    out.origins = origins.size();
    if (origins.empty()) {
        out.why = UnusableReason::NoInBandOrigins;
        // NO RECLASSIFICATION HERE, DELIBERATELY.
        //
        // The obvious rule -- "omitted_bound is tiny, so the fragment is neutral" -- compares a
        // RAW PLACEMENT MASS against an absolute constant, and that is the same conflation this
        // codebase already had to unlearn once in interval scoring. What the caller actually
        // evaluates is
        //
        //     log( (1-eta) * lambda * (M_a + M_b) + eta * P_bg )
        //
        // so a tail of e^-55 is negligible only RELATIVE to eta * P_bg. If the background floor
        // sits at e^-200, that same tail dominates the mixture and moves the fragment's
        // contribution by ~145 nats. Nor does "no in-band origins" mean the tail is
        // candidate-INDEPENDENT: it means the production-band search found nothing, and the
        // tail-only fixture already showed states = 0 alongside finite candidate-specific
        // reference mass.
        //
        // A defensible disposition needs the contribution INTERVAL through the real mixture, per
        // candidate, summed over every such fragment and compared to a declared global tolerance
        // -- which is what the ownership audit computes. Until that audit passes, the fragment
        // stays Unusable and keeps the call INCOMPLETE.
        return out;
    }

    // THE ESSENTIAL ORIGIN SET. Origins are taken in decreasing mass until everything still
    // excluded -- the remaining origins, the unmapped mass and the out-of-band bound together --
    // costs at most scope_tol nats off the total. The scope is the union of the blocks THOSE
    // origins span. Anything outside it is provably worth less than the declared tolerance, which
    // is what makes this a certified scope rather than a recruitment heuristic.
    std::sort(origins.begin(), origins.end(),
              [](const Origin& x, const Origin& y) { return x.mass > y.mass; });
    double total = out.in_band;
    if (out.omitted_bound != kNegInf) total = log_add(total, out.omitted_bound);
    double kept = kNegInf;
    std::size_t n_keep = 0;
    for (; n_keep < origins.size(); ++n_keep) {
        if (total - kept <= scope_tol) break;
        kept = log_add(kept, origins[n_keep].mass);
    }
    out.dropped = (kept == kNegInf) ? std::numeric_limits<double>::infinity() : total - kept;
    out.panel_domain_certified = out.dropped <= scope_tol;

    std::vector<std::uint32_t> sc;
    for (std::size_t i = 0; i < n_keep; ++i) {
        for (std::uint32_t b = origins[i].lo; b <= origins[i].hi; ++b) sc.push_back(b);
    }
    std::sort(sc.begin(), sc.end());
    sc.erase(std::unique(sc.begin(), sc.end()), sc.end());
    out.panel_domain_scope = sc;

    if (!out.panel_domain_certified) {
        out.kind = OwnerKind::Unusable; out.why = UnusableReason::ScopeNotCertified; return out;
    }
    if (sc.empty()) {
        out.kind = OwnerKind::Unusable; out.why = UnusableReason::EmptyScope; return out;
    }
    // Unattributable mass above the tolerance means a block-factored model cannot express this
    // fragment at all. It is NOT quietly assigned to the nearest block.
    if (out.unmapped != kNegInf && total - out.unmapped < scope_tol) {
        out.kind = OwnerKind::Unusable;
        out.why = UnusableReason::UnmappedMassDominant;
        return out;
    }

    // ARITY COMES FROM THE VARIABLE BLOCKS ONLY. A block whose sequence is identical across every
    // candidate carries no genotype state, so spanning it costs no arity: it enters the factor as
    // sequence context. Counting physical blocks instead would make a bubble--backbone--bubble
    // fragment Wide and throw away most of the real linkage evidence at any locus whose bubbles
    // have reference sequence between them, which is most of them.
    for (std::uint32_t b : sc) {
        if (b < block_variable.size() && block_variable[b]) out.panel_domain_var_scope.push_back(b);
    }
    if (out.panel_domain_var_scope.empty()) {
        // No genotype dependence at all. Still owned, because it still carries depth and still
        // consumes normalisation -- dropping it here would silently unbalance both.
        out.kind = OwnerKind::Invariant;
        out.block_lo = sc.front();
        out.block_hi = sc.back();
        return out;
    }
    out.block_lo = out.panel_domain_var_scope.front();
    out.block_hi = out.panel_domain_var_scope.back();
    if (out.panel_domain_var_scope.size() == 1) out.kind = OwnerKind::Unary;
    else if (out.panel_domain_var_scope.size() == 2) out.kind = OwnerKind::Linkage;
    else out.kind = OwnerKind::Wide;   // three or more variables; never cropped to a pair
    return out;
}

LinkageGeometry build_linkage_geometry(const std::vector<CandidateFrame>& frames,
                                       const std::vector<std::vector<std::string>>& block_alleles,
                                       const std::vector<char>& block_variable,
                                       std::uint32_t block_a, std::uint32_t block_b,
                                       std::size_t flank_bp, const InsertPrior& ip) {
    LinkageGeometry g;
    g.block_a = block_a; g.block_b = block_b; g.flank_bp = flank_bp;
    if (frames.empty() || block_a >= block_alleles.size() || block_b >= block_alleles.size()) {
        g.refusal = "block index out of range";
        return g;
    }
    // CHAIN ORIENTATION THROUGHOUT. Spans and sequence both come from the chain-oriented view, so
    // an antiparallel candidate is compared on the same footing as a forward one. Working in raw
    // walk coordinates made `a3 > b2` fire on every correct antiparallel frame.
    std::vector<std::string> chain_seq(frames.size());
    for (std::size_t h = 0; h < frames.size(); ++h) {
        chain_seq[h] = chain_oriented_sequence(frames[h]);
    }
    const auto span = [&](std::size_t h, std::uint32_t b, std::size_t& lo, std::size_t& hi) {
        return chain_oriented_block_span(frames[h], b, lo, hi);
    };
    std::size_t alo = 0, ahi = 0, blo = 0, bhi = 0;
    if (!span(0, block_a, alo, ahi) || !span(0, block_b, blo, bhi)) {
        g.refusal = "block absent from the reference candidate's map";
        return g;
    }
    if (ahi > blo) { g.refusal = "blocks overlap or are out of order"; return g; }
    const std::string& w0 = chain_seq[0];
    g.context = w0.substr(ahi, blo - ahi);

    // ---- DERIVE THE FLANK LENGTHS -------------------------------------------------------------
    // Each side is bounded by the nearest EXTERNAL VARIABLE BLOCK (never borrow a third variable as
    // context) and by the COMMON VERIFIED MAPPED boundary (an accepted partial terminal frame
    // verifies less than the catalogue holds, and unverified bytes must not become context), then
    // capped by the requested flank_bp. Zero is a correct derived answer, not a failure.
    long prev_var = -1, next_var = -1;
    for (long b = static_cast<long>(block_a) - 1; b >= 0; --b) {
        if (static_cast<std::size_t>(b) < block_variable.size() && block_variable[b]) {
            prev_var = b; break;
        }
    }
    for (std::size_t b = block_b + 1; b < block_variable.size(); ++b) {
        if (block_variable[b]) { next_var = static_cast<long>(b); break; }
    }
    std::size_t lmax = flank_bp, rmax = flank_bp;
    for (std::size_t h = 0; h < frames.size(); ++h) {
        std::size_t a2 = 0, a3 = 0, b2 = 0, b3 = 0;
        if (!span(h, block_a, a2, a3) || !span(h, block_b, b2, b3)) {
            g.refusal = "block absent from candidate " + std::to_string(h);
            return g;
        }
        std::size_t mlo = frames[h].mapped_lo, mhi = frames[h].mapped_hi;
        if (frames[h].reverse_frame) {
            const std::size_t n = frames[h].seq.size();
            const std::size_t t = n - mhi;
            mhi = n - mlo;
            mlo = t;
        }
        std::size_t lo_limit = mlo, hi_limit = mhi;
        if (prev_var >= 0) {
            std::size_t p2 = 0, p3 = 0;
            if (span(h, static_cast<std::uint32_t>(prev_var), p2, p3)) {
                lo_limit = std::max(lo_limit, p3);
            }
        }
        if (next_var >= 0) {
            std::size_t n2 = 0, n3 = 0;
            if (span(h, static_cast<std::uint32_t>(next_var), n2, n3)) {
                hi_limit = std::min(hi_limit, n2);
            }
        }
        lmax = std::min(lmax, a2 > lo_limit ? a2 - lo_limit : 0);
        rmax = std::min(rmax, hi_limit > b3 ? hi_limit - b3 : 0);
    }
    g.lflank_bp = lmax;
    g.rflank_bp = rmax;
    g.lflank = w0.substr(alo - lmax, lmax);
    g.rflank = w0.substr(bhi, rmax);
    // THE CONTEXT AND FLANKS MUST AGREE ACROSS EVERY CANDIDATE. Where they do not, this edge is not
    // representable by a pairwise factor and is REFUSED. Guessing one candidate's context would let
    // the factor express a preference that belongs to a block outside it.
    for (std::size_t h = 1; h < frames.size(); ++h) {
        std::size_t a2 = 0, a3 = 0, b2 = 0, b3 = 0;
        if (!span(h, block_a, a2, a3) || !span(h, block_b, b2, b3)) {
            g.refusal = "block absent from candidate " + std::to_string(h); return g;
        }
        if (a3 > b2) { g.refusal = "blocks out of order on candidate " + std::to_string(h); return g; }
        const std::string& wh = chain_seq[h];
        if (wh.substr(a3, b2 - a3) != g.context) {
            g.refusal = "intervening context differs on candidate " + std::to_string(h); return g;
        }
        // THE BYTE COMPARISON REMAINS, though the flank is now derived structurally. "Invariant
        // per the block catalogue" and "identical in every authoritative walk" should agree; a
        // disagreement is a refusal, not something to reconcile silently.
        if (a2 < lmax || wh.substr(a2 - lmax, lmax) != g.lflank ||
            wh.substr(b3, rmax) != g.rflank) {
            g.refusal = "flank differs on candidate " + std::to_string(h); return g;
        }
    }
    g.alleles_a = block_alleles[block_a];
    g.alleles_b = block_alleles[block_b];
    const std::size_t na = g.alleles_a.size(), nb = g.alleles_b.size();
    if (na == 0 || nb == 0) { g.refusal = "a block offers no allele"; return g; }
    // EXPOSURE, ONCE PER CONFIGURATION, on the edge. Exact in both regimes; the affine surrogate is
    // only recorded so the cancellation claim stays auditable.
    g.window_len.assign(na * nb, 0);
    g.exposure.assign(na * nb, 0.0);
    bool affine = true;
    for (std::size_t al = 0; al < na; ++al) {
        for (std::size_t be = 0; be < nb; ++be) {
            const std::size_t idx = al * nb + be;
            g.window_len[idx] = g.lflank.size() + g.alleles_a[al].size() + g.context.size() +
                                g.alleles_b[be].size() + g.rflank.size();
            const ExposureCheck ec = check_exposure(g.window_len[idx], ip);
            g.exposure[idx] = ec.exact;
            if (!ec.in_regime) affine = false;
        }
    }
    g.exposure_affine = affine;
    g.ok = true;
    return g;
}

namespace {

// Every exact occurrence of `needle` in `hay`. Occurrences are NOT deduplicated: two identical
// repeat copies are two distinct origins, and collapsing them would lose multiplicity.
void find_all(const std::string& hay, const std::string& needle, std::vector<std::size_t>& out) {
    if (needle.empty() || hay.size() < needle.size()) return;
    for (std::size_t at = hay.find(needle); at != std::string::npos;
         at = hay.find(needle, at + 1)) {
        out.push_back(at);
    }
}

inline bool acgt_only(const std::string& s) {
    for (char c : s) {
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') return false;
    }
    return true;
}

}  // namespace

namespace {

// 2 bits per base; returns false on any non-ACGT, which is what sends the caller to the exhaustive
// fallback rather than to a quietly different search.
inline bool encode_piece(const std::string& s, std::size_t at, std::size_t len, std::uint64_t& out) {
    if (len == 0 || len > 32 || at + len > s.size()) return false;
    std::uint64_t v = 0;
    for (std::size_t i = 0; i < len; ++i) {
        int c;
        switch (s[at + i]) {
            case 'A': c = 0; break;
            case 'C': c = 1; break;
            case 'G': c = 2; break;
            case 'T': c = 3; break;
            default: return false;
        }
        v = (v << 2) | static_cast<std::uint64_t>(c);
    }
    out = v;
    return true;
}

void index_all_positions(const std::string& seq, std::size_t piece, std::uint32_t id,
                         std::unordered_map<std::uint64_t, std::vector<std::uint32_t>>& into) {
    if (seq.size() < piece) return;
    for (std::size_t at = 0; at + piece <= seq.size(); ++at) {
        std::uint64_t code = 0;
        if (!encode_piece(seq, at, piece, code)) continue;
        auto& v = into[code];
        if (v.empty() || v.back() != id) v.push_back(id);
    }
}

}  // namespace

void VirtualWindow::bind_pair(const LinkageGeometry& g, std::uint32_t a, std::uint32_t b) {
    geom = &g; alpha = a; beta = b;
    segments.clear();
    segments.push_back(&g.lflank);
    segments.push_back(&g.alleles_a[a]);
    segments.push_back(&g.context);
    segments.push_back(&g.alleles_b[b]);
    segments.push_back(&g.rflank);
}

void VirtualWindow::bind_chain(const std::string& lflank,
                               const std::vector<const std::string*>& alleles,
                               const std::vector<const std::string*>& contexts,
                               const std::string& rflank) {
    geom = nullptr;
    segments.clear();
    segments.push_back(&lflank);
    for (std::size_t i = 0; i < alleles.size(); ++i) {
        segments.push_back(alleles[i]);
        if (i + 1 < alleles.size()) segments.push_back(contexts[i]);
    }
    segments.push_back(&rflank);
}

std::size_t VirtualWindow::size() const {
    if (!segments.empty()) {
        std::size_t n = 0;
        for (const std::string* s : segments) n += s->size();
        return n;
    }
    if (geom == nullptr) return 0;
    return geom->lflank.size() + geom->alleles_a[alpha].size() + geom->context.size() +
           geom->alleles_b[beta].size() + geom->rflank.size();
}

char VirtualWindow::base_at(std::size_t pos) const {
    if (!segments.empty()) {
        for (const std::string* s : segments) {
            if (pos < s->size()) return (*s)[pos];
            pos -= s->size();
        }
        return 'N';   // out of range; callers bound-check before entering
    }
    const std::size_t lL = geom->lflank.size();
    if (pos < lL) return geom->lflank[pos];
    pos -= lL;
    const std::string& A = geom->alleles_a[alpha];
    if (pos < A.size()) return A[pos];
    pos -= A.size();
    if (pos < geom->context.size()) return geom->context[pos];
    pos -= geom->context.size();
    const std::string& B = geom->alleles_b[beta];
    if (pos < B.size()) return B[pos];
    pos -= B.size();
    return geom->rflank[pos];
}

std::size_t VirtualWindow::count_mismatches(const std::string& read, long start,
                                            std::size_t cap) const {
    const std::size_t n = size();
    if (start < 0 || static_cast<std::size_t>(start) + read.size() > n) return cap + 1;
    std::size_t mm = 0;
    for (std::size_t i = 0; i < read.size(); ++i) {
        if (base_at(static_cast<std::size_t>(start) + i) != read[i]) {
            if (++mm > cap) return mm;   // give up as soon as the band is exceeded
        }
    }
    return mm;
}

std::string VirtualWindow::materialize() const {
    if (!segments.empty()) {
        std::string out;
        out.reserve(size());
        for (const std::string* s : segments) out += *s;
        return out;
    }
    return geom->lflank + geom->alleles_a[alpha] + geom->context + geom->alleles_b[beta] +
           geom->rflank;
}

AlleleProductIndex build_allele_product_index(const LinkageGeometry& geom, std::size_t piece) {
    AlleleProductIndex ix;
    ix.piece = piece;
    if (!geom.ok || piece == 0 || piece > 32) return ix;
    const auto index_into = [&](const std::string& seq, std::uint32_t id,
                                std::unordered_map<std::uint64_t,
                                                   std::vector<AlleleSeedHit>>& into) {
        if (seq.size() < piece) return;
        for (std::size_t at = 0; at + piece <= seq.size(); ++at) {
            std::uint64_t code = 0;
            if (!encode_piece(seq, at, piece, code)) continue;
            into[code].push_back(AlleleSeedHit{id, static_cast<std::uint32_t>(at)});
        }
    };
    const auto index_plain = [&](const std::string& seq,
                                 std::unordered_map<std::uint64_t,
                                                    std::vector<std::uint32_t>>& into) {
        if (seq.size() < piece) return;
        for (std::size_t at = 0; at + piece <= seq.size(); ++at) {
            std::uint64_t code = 0;
            if (!encode_piece(seq, at, piece, code)) continue;
            into[code].push_back(static_cast<std::uint32_t>(at));
        }
    };
    index_plain(geom.lflank, ix.in_l);
    index_plain(geom.context, ix.in_c);
    index_plain(geom.rflank, ix.in_r);
    // ALLELE CONTEXTS. One string per allele, wide enough that every piece TOUCHING that allele
    // falls inside it: a piece touching A starts no earlier than |L|-(p-1) and, when the context is
    // at least p-1 long, ends no later than |L|+|A|+(p-1).
    const auto tail = [&](const std::string& x, std::size_t n) {
        return x.size() <= n ? x : x.substr(x.size() - n);
    };
    const auto head = [&](const std::string& x, std::size_t n) {
        return x.size() <= n ? x : x.substr(0, n);
    };
    const std::size_t back = piece - 1;
    for (std::uint32_t al = 0; al < geom.alleles_a.size(); ++al) {
        index_into(tail(geom.lflank, back) + geom.alleles_a[al] + head(geom.context, back), al,
                   ix.a_ctx);
    }
    for (std::uint32_t be = 0; be < geom.alleles_b.size(); ++be) {
        index_into(tail(geom.context, back) + geom.alleles_b[be] + head(geom.rflank, back), be,
                   ix.b_ctx);
    }
    // THE DIRECT A->B JUNCTION, split at every interior point. The two sides are looked up
    // separately and their PRODUCT is proposed -- which keeps the alpha-beta correlation that
    // flattening into independent flags destroys.
    ix.a_suffix.assign(piece, {});
    ix.b_prefix.assign(piece, {});
    for (std::size_t j = 1; j < piece; ++j) {
        for (std::uint32_t al = 0; al < geom.alleles_a.size(); ++al) {
            const std::string& A = geom.alleles_a[al];
            std::uint64_t code = 0;
            if (A.size() >= j && encode_piece(A, A.size() - j, j, code)) {
                ix.a_suffix[j][code].push_back(al);
            }
        }
        for (std::uint32_t be = 0; be < geom.alleles_b.size(); ++be) {
            std::uint64_t code = 0;
            if (encode_piece(geom.alleles_b[be], 0, piece - j, code)) {
                ix.b_prefix[j][code].push_back(be);
            }
        }
    }
    // THE COMPLETENESS PREDICATE, stated as the span shape it excludes. A piece that touches L and
    // reaches past C into B needs |A_alpha| + |C| < p - 1; symmetrically on the right. Everything
    // else is reachable: inside one component, inside an allele context, or across the A->B
    // junction above (which requires the piece to start in A and end in B, both guaranteed there).
    ix.complete = true;
    const std::size_t lC = geom.context.size();
    for (const std::string& A : geom.alleles_a)
        if (A.size() + lC < piece - 1) ix.complete = false;
    for (const std::string& B : geom.alleles_b)
        if (lC + B.size() < piece - 1) ix.complete = false;
    ix.ok = true;
    return ix;
}

AlleleProductSupport propose_allele_pairs(const Fragment& fragment, const LinkageGeometry& geom,
                                          const InsertPrior& ip, double max_divergence,
                                          const AlleleProductIndex* index) {
    AlleleProductSupport S;
    const std::size_t na = geom.alleles_a.size(), nb = geom.alleles_b.size();
    S.dense_pairs = na * nb;
    if (!geom.ok || na == 0 || nb == 0 || fragment.r1.empty() || fragment.r2.empty()) {
        S.exhaustive_fallback = true;
        return S;
    }
    if (index == nullptr || !index->ok || !index->complete) {
        // An incomplete index cannot be verified against: a piece it cannot reach is a placement
        // the direct verifier would silently never test.
        S.exhaustive_fallback = true; return S;
    }
    const std::size_t d1 = mate_band_edits(max_divergence, fragment.r1.size());
    const std::size_t d2 = mate_band_edits(max_divergence, fragment.r2.size());
    const std::string rc1 = reverse_complement(fragment.r1), rc2 = reverse_complement(fragment.r2);
    struct MateSpec { const std::string* seq; std::size_t d; };
    const MateSpec specs[4] = {{&fragment.r1, d1}, {&rc1, d1}, {&fragment.r2, d2}, {&rc2, d2}};
    for (const MateSpec& m : specs) {
        if (!acgt_only(*m.seq)) { S.exhaustive_fallback = true; return S; }
    }
    if (!acgt_only(geom.lflank) || !acgt_only(geom.context) || !acgt_only(geom.rflank)) {
        S.exhaustive_fallback = true; return S;
    }
    const std::size_t p = index->piece;
    for (const MateSpec& m : specs) {
        if (m.seq->size() / (m.d + 1) != p) { S.exhaustive_fallback = true; return S; }
    }
    const std::size_t lL = geom.lflank.size(), lC = geom.context.size();
    const std::size_t lR = geom.rflank.size();
    (void)lR;

    // POSITIONAL MATE STATES, per mate variant: (alpha, beta, start). A dimension a seed does not
    // constrain is EXPANDED into concrete tuples -- correctness first. Symbolic wildcards are a
    // later optimisation and only if profiling shows tuple construction dominating.
    std::vector<std::vector<std::pair<std::uint64_t, long>>> states(4);
    const auto key = [&](std::uint32_t al, std::uint32_t be) {
        return (static_cast<std::uint64_t>(al) << 32) | be;
    };
    for (int mi = 0; mi < 4; ++mi) {
        const std::string& r = *specs[mi].seq;
        const std::size_t d = specs[mi].d;
        auto& out = states[mi];
        for (std::size_t k = 0; k <= d; ++k) {
            if (k * p + p > r.size()) break;
            std::uint64_t code = 0;
            if (!encode_piece(r, k * p, p, code)) { S.exhaustive_fallback = true; return S; }
            const long back = static_cast<long>(k * p);
            // --- invariant components: no allele constrained, but the start IS fixed -----------
            {
                const auto it = index->in_l.find(code);
                if (it != index->in_l.end()) {
                    S.seed_occurrences += it->second.size();
                    for (std::uint32_t o : it->second) {
                        const long st = static_cast<long>(o) - back;
                        for (std::uint32_t x = 0; x < na; ++x)
                            for (std::uint32_t y = 0; y < nb; ++y) out.emplace_back(key(x, y), st);
                    }
                }
            }
            {
                const auto it = index->in_c.find(code);
                if (it != index->in_c.end()) {
                    S.seed_occurrences += it->second.size();
                    for (std::uint32_t o : it->second) {
                        // The context sits after A_alpha, so its window coordinate depends on alpha.
                        for (std::uint32_t x = 0; x < na; ++x) {
                            const long st = static_cast<long>(lL + geom.alleles_a[x].size() + o) -
                                            back;
                            for (std::uint32_t y = 0; y < nb; ++y) out.emplace_back(key(x, y), st);
                        }
                    }
                }
            }
            {
                const auto it = index->in_r.find(code);
                if (it != index->in_r.end()) {
                    S.seed_occurrences += it->second.size();
                    for (std::uint32_t o : it->second) {
                        for (std::uint32_t x = 0; x < na; ++x)
                            for (std::uint32_t y = 0; y < nb; ++y) {
                                const long st = static_cast<long>(
                                    lL + geom.alleles_a[x].size() + lC +
                                    geom.alleles_b[y].size() + o) - back;
                                out.emplace_back(key(x, y), st);
                            }
                    }
                }
            }
            // --- the A context: constrains alpha, beta free -----------------------------
            // The context string starts at window  |L| - min(|L|, p-1), which does NOT depend on
            // alpha: everything left of A is invariant. What depends on alpha is the B side below.
            {
                const auto it = index->a_ctx.find(code);
                if (it != index->a_ctx.end()) {
                    S.seed_occurrences += it->second.size();
                    const long base = static_cast<long>(lL - std::min(lL, p - 1));
                    for (const AlleleSeedHit& h : it->second) {
                        const long st = base + static_cast<long>(h.offset) - back;
                        for (std::uint32_t y = 0; y < nb; ++y) out.emplace_back(key(h.allele, y), st);
                    }
                }
            }
            // --- the B context: constrains beta; its coordinate SHIFTS by |A_alpha| ------------
            // This shift is the discrimination the positional search adds over a set of flags. It
            // cancels in the insert length when both mates sit inside B, so proposing every pair
            // there is a genuine answer, not a lost filter.
            {
                const auto it = index->b_ctx.find(code);
                if (it != index->b_ctx.end()) {
                    S.seed_occurrences += it->second.size();
                    const long base = static_cast<long>(lL + lC - std::min(lC, p - 1));
                    for (const AlleleSeedHit& h : it->second) {
                        for (std::uint32_t x = 0; x < na; ++x) {
                            const long st = base + static_cast<long>(geom.alleles_a[x].size()) +
                                            static_cast<long>(h.offset) - back;
                            out.emplace_back(key(x, h.allele), st);
                        }
                    }
                }
            }
            // --- the direct A->B junction, alpha and beta CORRELATED --------------------------
            if (lC < p) {
                for (std::size_t j = 1; j < p; ++j) {
                    std::uint64_t lc = 0, rc = 0;
                    if (!encode_piece(r, k * p, j, lc)) continue;
                    const std::size_t rlen = p - j;
                    if (lC != 0) {
                        if (rlen <= lC) continue;
                        if (r.compare(k * p + j, lC, geom.context) != 0) continue;
                        if (!encode_piece(r, k * p + j + lC, rlen - lC, rc)) continue;
                    } else {
                        if (!encode_piece(r, k * p + j, rlen, rc)) continue;
                    }
                    if (j >= index->a_suffix.size()) continue;
                    const auto ia = index->a_suffix[j].find(lc);
                    if (ia == index->a_suffix[j].end()) continue;
                    const std::size_t bj = lC != 0 ? j + lC : j;
                    if (bj >= index->b_prefix.size()) continue;
                    const auto ib = index->b_prefix[bj].find(rc);
                    if (ib == index->b_prefix[bj].end()) continue;
                    S.seed_occurrences += ia->second.size() * ib->second.size();
                    for (std::uint32_t x : ia->second) {
                        const long st = static_cast<long>(
                            lL + geom.alleles_a[x].size() - j) - back;
                        for (std::uint32_t y : ib->second) out.emplace_back(key(x, y), st);
                    }
                }
            }
        }
        S.positional_states_before_dedup += out.size();
        S.seed_start_proposals += out.size();
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
        S.positional_states_after_dedup += out.size();
        S.unique_seed_starts += out.size();
    }

    // THE VALID-FR JOIN, per orientation, applied BEFORE the pair set is formed. Intersecting mate
    // reachability first -- as the previous version did -- discards exactly the positional
    // relationship the insert prior needs, which is why it proposed every pair on real C4.
    std::vector<char> keep(na * nb, 0);
    const auto join = [&](const std::vector<std::pair<std::uint64_t, long>>& fwd,
                          const std::vector<std::pair<std::uint64_t, long>>& rev,
                          std::size_t rev_len) {
        std::size_t fi = 0, ri = 0;
        while (fi < fwd.size() && ri < rev.size()) {
            if (fwd[fi].first < rev[ri].first) { ++fi; continue; }
            if (rev[ri].first < fwd[fi].first) { ++ri; continue; }
            const std::uint64_t kk = fwd[fi].first;
            std::size_t fe = fi, re = ri;
            while (fe < fwd.size() && fwd[fe].first == kk) ++fe;
            while (re < rev.size() && rev[re].first == kk) ++re;
            const std::uint32_t al = static_cast<std::uint32_t>(kk >> 32);
            const std::uint32_t be = static_cast<std::uint32_t>(kk & 0xFFFFFFFFu);
            for (std::size_t a = fi; a < fe && !keep[al * nb + be]; ++a) {
                for (std::size_t b = ri; b < re; ++b) {
                    const long rev_end = rev[b].second + static_cast<long>(rev_len) - 1;
                    // THE SHARED RULE, not a second copy of it.
                    if (valid_fr_coordinates(fwd[a].second, rev_end, ip.lo, ip.hi)) {
                        ++S.seed_compatible_joins;
                        keep[al * nb + be] = 1;
                        break;
                    }
                }
            }
            fi = fe; ri = re;
        }
    };
    join(states[0], states[3], fragment.r2.size());   // r1 forward with r2 reverse-complemented
    join(states[2], states[1], fragment.r1.size());   // r2 forward with r1 reverse-complemented

    std::vector<std::pair<std::uint32_t, std::uint32_t>> pairs;
    for (std::uint32_t x = 0; x < na; ++x)
        for (std::uint32_t y = 0; y < nb; ++y)
            if (keep[x * nb + y]) pairs.emplace_back(x, y);
    S.seed_hits = S.seed_occurrences;
    S.proposed_states = pairs.size();
    S.proposals = std::move(pairs);
    // The positional states are RETAINED: the emission verifies at these starts rather than
    // searching each window again, which is the whole point of having derived them.
    S.mate_states = std::move(states);
    return S;
}

namespace {

// One (alpha, beta) window, scored. Shared by the dense oracle and the support-restricted form, so
// the two cannot disagree about how a cell is computed -- only about which cells are visited.
double score_window(const Fragment& fragment, const LinkageGeometry& geom, const InsertPrior& ip,
                    std::size_t al, std::size_t be, std::size_t d1, std::size_t d2,
                    const std::string& a1, const std::string& a2,
                    double log_eps, double log_1meps, std::uint32_t* n_states = nullptr) {
    const std::string win = geom.lflank + geom.alleles_a[al] + geom.context +
                            geom.alleles_b[be] + geom.rflank;
    // ONE PIECE INDEX PER WINDOW, used by all four mate searches. Reverted once on the arithmetic
    // that at 1.15M dense windows an O(|win|) build matches the scan it replaces -- true for the
    // dense path, wrong once the support search cuts the count. Profiling the running process is
    // what settled it.
    const std::size_t piece = fragment.r1.size() / (d1 + 1);
    PieceIndex widx;
    const PieceIndex* wp = nullptr;
    if (piece >= 12 && win.size() > 8 * piece) {
        widx = build_piece_index(win, piece);
        wp = &widx;
    }
    const auto f1 = bounded_mate_placements(fragment.r1, win, d1, nullptr, wp);
    const auto v1 = bounded_mate_placements(a1, win, d1, nullptr, wp);
    const auto f2 = bounded_mate_placements(fragment.r2, win, d2, nullptr, wp);
    const auto v2 = bounded_mate_placements(a2, win, d2, nullptr, wp);
    // THE SHARED valid-FR RULE, through enumerate_fragment_states.
    const auto st = enumerate_fragment_states(0, f1, v1, f2, v2, fragment.r1.size(),
                                              fragment.r2.size(), ip.lo, ip.hi);
    const double half = std::log(0.5);
    double m = kNegInf;
    if (n_states != nullptr) *n_states = static_cast<std::uint32_t>(st.size());
    // EVERY DISTINCT ORIGIN, with its multiplicity: identical repeat copies are separate states.
    for (const FragmentState& z : st) {
        const double e1 = static_cast<double>(z.m1_edits) * log_eps +
                          static_cast<double>(fragment.r1.size() - z.m1_edits) * log_1meps;
        const double e2 = static_cast<double>(z.m2_edits) * log_eps +
                          static_cast<double>(fragment.r2.size() - z.m2_edits) * log_1meps;
        m = log_add(m, half + e1 + e2 + ip.log_at(z.insert));
    }
    return m;
}

// Membership in the seed-compatible pair set, kept as a lookup rather than a second copy of the
// selection rule.
bool keep_pair(const AlleleProductSupport& sup, std::uint32_t al, std::uint32_t be,
               std::size_t n_b) {
    (void)n_b;
    const auto want = std::make_pair(al, be);
    return std::binary_search(sup.proposals.begin(), sup.proposals.end(), want);
}

void mark_informative(LinkageEmission& out) {
    out.informative = false;
    for (std::size_t al = 0; al < out.n_a && !out.informative; ++al) {
        for (std::size_t be = 1; be < out.n_b; ++be) {
            const double d = out.mass[al * out.n_b + be] - out.mass[al * out.n_b];
            if (!(std::abs(d) < 1e-12) &&
                !(out.mass[al * out.n_b + be] == kNegInf && out.mass[al * out.n_b] == kNegInf)) {
                out.informative = true; break;
            }
        }
    }
    if (!out.informative) return;
    bool varies_a = false;
    for (std::size_t be = 0; be < out.n_b && !varies_a; ++be) {
        for (std::size_t al = 1; al < out.n_a; ++al) {
            const double d = out.mass[al * out.n_b + be] - out.mass[be];
            if (!(std::abs(d) < 1e-12) &&
                !(out.mass[al * out.n_b + be] == kNegInf && out.mass[be] == kNegInf)) {
                varies_a = true; break;
            }
        }
    }
    out.informative = varies_a;
}

}  // namespace

std::size_t IntervalGeometry::cell_index(const std::vector<std::uint32_t>& choice) const {
    std::size_t idx = 0;
    for (std::size_t j = 0; j < alleles.size(); ++j) idx = idx * alleles[j].size() + choice[j];
    return idx;
}

void IntervalGeometry::cell_choice(std::size_t index, std::vector<std::uint32_t>& out) const {
    out.assign(alleles.size(), 0);
    for (std::size_t j = alleles.size(); j-- > 0;) {
        out[j] = static_cast<std::uint32_t>(index % alleles[j].size());
        index /= alleles[j].size();
    }
}

std::size_t IntervalGeometry::window_len(const std::vector<std::uint32_t>& choice) const {
    std::size_t n = lflank.size() + rflank.size();
    for (std::size_t j = 0; j < alleles.size(); ++j) {
        n += alleles[j][choice[j]].size();
        if (j + 1 < alleles.size()) n += contexts[j].size();
    }
    return n;
}

IntervalGeometry build_interval_geometry(const std::vector<CandidateFrame>& frames,
                                         const std::vector<std::vector<std::string>>& block_alleles,
                                         const std::vector<char>& block_variable,
                                         const std::vector<std::uint32_t>& blocks,
                                         std::size_t flank_bp, const InsertPrior& ip) {
    IntervalGeometry G;
    G.blocks = blocks;
    if (blocks.size() < 2) { G.refusal = "an interval factor needs at least two blocks"; return G; }
    for (std::size_t j = 0; j + 1 < blocks.size(); ++j) {
        if (blocks[j] >= blocks[j + 1]) {
            G.refusal = "blocks must be ascending"; return G;
        }
    }
    for (std::uint32_t b : blocks) {
        if (b >= block_variable.size() || !block_variable[b]) {
            G.refusal = "block " + std::to_string(b) + " is not variable"; return G;
        }
    }
    // THE CONSECUTIVE PAIRWISE GEOMETRIES supply the contexts and the outer flanks. Deriving them
    // independently here would be a second opinion about the same sequence, and the two would
    // eventually disagree -- which is how a non-adjacent geometry once went unnoticed.
    for (std::size_t j = 0; j + 1 < blocks.size(); ++j) {
        const LinkageGeometry g = build_linkage_geometry(frames, block_alleles, block_variable,
                                                         blocks[j], blocks[j + 1], flank_bp, ip);
        if (!g.ok) {
            G.refusal = "pair " + std::to_string(blocks[j]) + "-" +
                        std::to_string(blocks[j + 1]) + ": " + g.refusal;
            return G;
        }
        G.contexts.push_back(g.context);
        if (j == 0) G.lflank = g.lflank;
        if (j + 2 == blocks.size()) G.rflank = g.rflank;
    }
    for (std::uint32_t b : blocks) G.alleles.push_back(block_alleles[b]);
    // EVERY window, by complete enumeration. A sampled affine gate inverted a conclusion in this
    // work once; there is no sampled path here.
    std::vector<std::uint32_t> choice(G.alleles.size(), 0);
    G.min_window = SIZE_MAX; G.max_window = 0; G.windows_below_affine = 0;
    const long affine_from = std::max<long>(0, ip.hi - 1);
    bool done = false;
    while (!done) {
        const std::size_t w = G.window_len(choice);
        G.min_window = std::min(G.min_window, w);
        G.max_window = std::max(G.max_window, w);
        if (static_cast<long>(w) < affine_from) ++G.windows_below_affine;
        for (std::size_t j = 0; ; ++j) {
            if (j == choice.size()) { done = true; break; }
            if (++choice[j] < G.alleles[j].size()) break;
            choice[j] = 0;
        }
    }
    G.exposure_affine = G.windows_below_affine == 0;
    G.ok = true;
    return G;
}

namespace {

// THE PER-ORIGIN LOG CONTRIBUTION, in ONE place. Both paths sum the same terms, but floating-point
// addition is not associative: grouping them differently makes two mathematically identical
// contributions differ in the last bits, and then a per-origin comparison can only be made with a
// tolerance. Computing it here makes that comparison EXACT, which is what justifies the aggregate
// mass tolerance instead of the tolerance justifying itself.
inline double origin_contribution(std::uint32_t e1, std::size_t len1, std::uint32_t e2,
                                  std::size_t len2, long insert, double log_eps, double log_1meps,
                                  const InsertPrior& ip) {
    const double a = static_cast<double>(e1) * log_eps +
                     static_cast<double>(len1 - e1) * log_1meps;
    const double b = static_cast<double>(e2) * log_eps +
                     static_cast<double>(len2 - e2) * log_1meps;
    return std::log(0.5) + a + b + ip.log_at(insert);
}

// One verified placement of one mate: the alleles it pins, where it starts inside the FIRST block
// it touches, and its edit count. Blocks outside [first, last] are untouched and stay free.
struct IntervalPlacement {
    std::uint32_t first_block = 0, last_block = 0;
    std::vector<std::uint32_t> alleles;   // one per block in [first_block, last_block]
    long offset_in_first = 0;             // may be 0..|A_first|-1
    std::uint32_t edits = 0;
};

}  // namespace

IntervalEmission interval_emission(const Fragment& fragment, const IntervalGeometry& geom,
                                   const InsertPrior& ip, double max_divergence,
                                   double log_eps, double log_1meps, double log_p_bg,
                                   const IntervalSeedIndex* index, HybridWorkBudget* budget,
                                   bool want_signatures, bool want_origins) {
    IntervalEmission out;
    out.log_p_bg = log_p_bg;
    if (!geom.ok || fragment.r1.empty() || fragment.r2.empty()) return out;
    out.cells = geom.cells();
    out.mass.assign(out.cells, kNegInf);
    out.cell_states.assign(out.cells, 0);
    const std::size_t k = geom.alleles.size();
    const auto refuse = [&](const char* why) {
        out.mass.assign(out.cells, kNegInf);
        out.cell_states.assign(out.cells, 0);
        out.cell_signature.clear();
        out.work_refused = true;
        out.refusal = why;
        out.ok = false;
        return out;
    };
    // AN INCOMPLETE INDEX CANNOT BE VERIFIED AGAINST. A seed shape it cannot reach is a placement
    // this path would never test, and a silently missing placement is wrong mass -- so it refuses
    // and the caller falls back or declines under its own budget.
    if (index == nullptr || !index->ok) return refuse("interval-seed-index-missing");
    if (!index->complete) return refuse("interval-seed-index-incomplete");
    const std::size_t p = index->piece;
    const std::size_t d1 = mate_band_edits(max_divergence, fragment.r1.size());
    const std::size_t d2 = mate_band_edits(max_divergence, fragment.r2.size());
    const std::string rc1 = reverse_complement(fragment.r1), rc2 = reverse_complement(fragment.r2);
    const std::string* mv[4] = {&fragment.r1, &rc1, &fragment.r2, &rc2};
    const std::size_t bd[4] = {d1, d1, d2, d2};
    for (int m = 0; m < 4; ++m) {
        if (!acgt_only(*mv[m])) return refuse("non-acgt-read");
        if (mv[m]->size() / (bd[m] + 1) != p) return refuse("piece-length-mismatch");
    }
    for (const std::string& c : geom.contexts)
        if (!acgt_only(c)) return refuse("non-acgt-context");

    // Segment length helper: block j's allele plus the context that follows it.
    const auto seg = [&](std::size_t j, std::uint32_t a) {
        return geom.alleles[j][a].size() + (j + 1 < k ? geom.contexts[j].size() : 0);
    };
    // Verify one read over a concrete run of blocks, through the shared window walk.
    const auto verify_run = [&](const std::string& read, std::size_t first,
                                const std::vector<std::uint32_t>& al, long off,
                                std::size_t cap, std::size_t* edits) {
        // THE TRAILING CONTEXT BELONGS TO THE WINDOW. seg(j) counts block j's allele PLUS the
        // context that follows it, so a run's span includes that context -- and a read reaching
        // into it was being verified against a window that stopped at the allele, failing to fit
        // and vanishing. The span and the window must be built from the same definition.
        std::vector<const std::string*> alle, ctxp;
        for (std::size_t j = 0; j < al.size(); ++j) {
            alle.push_back(&geom.alleles[first + j][al[j]]);
            if (j + 1 < al.size()) ctxp.push_back(&geom.contexts[first + j]);
        }
        static const std::string kEmpty;
        const std::size_t last = first + al.size() - 1;
        const std::string& trail = last + 1 < k ? geom.contexts[last] : geom.rflank;
        VirtualWindow vw;
        vw.bind_chain(kEmpty, alle, ctxp, trail);
        const std::size_t mm = vw.count_mismatches(read, off, cap);
        if (mm > cap) return false;
        *edits = mm;
        return true;
    };
    // From a seed anchor, enumerate the concrete allele runs that cover the read, expanding left
    // while the start is negative and right while the read is not yet covered. Only the blocks the
    // read TOUCHES are enumerated; every other dimension stays free.
    std::vector<IntervalPlacement> place[4];
    for (int m = 0; m < 4; ++m) {
        const std::string& read = *mv[m];
        const std::size_t cap = bd[m];
        std::vector<std::pair<std::size_t, IntervalSeedHit>> hits;
        for (std::size_t q = 0; q <= cap; ++q) {
            if (q * p + p > read.size()) break;
            std::uint64_t code = 0;
            if (!encode_piece(read, q * p, p, code)) return refuse("non-acgt-seed");
            const auto ii = index->inside.find(code);
            if (ii != index->inside.end())
                for (const IntervalSeedHit& h : ii->second) hits.emplace_back(q, h);
            const auto ib = index->boundary.find(code);
            if (ib != index->boundary.end())
                for (const IntervalSeedHit& h : ib->second) hits.emplace_back(q, h);
        }
        out.seed_hits += hits.size();
        out.symbolic_states += hits.size();
        for (const auto& hq : hits) {
            const std::size_t q = hq.first;
            const IntervalSeedHit& h = hq.second;
            // Start of the read relative to the anchor block's own beginning.
            long s = static_cast<long>(h.offset) - static_cast<long>(q * p);
            std::size_t first = h.block;
            std::vector<std::uint32_t> al{h.allele};
            if (h.spans_boundary) al.push_back(h.next_allele);
            // Expand LEFT while the read starts before this block.
            std::function<void(std::size_t, std::vector<std::uint32_t>&, long)> grow_left =
                [&](std::size_t fb, std::vector<std::uint32_t>& acc, long soff) {
                    if (soff >= 0 || fb == 0) {
                        // Expand RIGHT until the read is covered, then verify.
                        std::function<void(std::vector<std::uint32_t>&)> grow_right =
                            [&](std::vector<std::uint32_t>& run) {
                                std::size_t span = 0;
                                for (std::size_t j = 0; j < run.size(); ++j)
                                    span += seg(fb + j, run[j]);
                                const long need = soff + static_cast<long>(read.size());
                                if (need > static_cast<long>(span) && fb + run.size() < k) {
                                    for (std::uint32_t a2 = 0;
                                         a2 < geom.alleles[fb + run.size()].size(); ++a2) {
                                        run.push_back(a2);
                                        ++out.tuple_expansions;
                                        grow_right(run);
                                        run.pop_back();
                                    }
                                    return;
                                }
                                if (soff < 0 || need > static_cast<long>(span)) return;
                                if (budget != nullptr &&
                                    !budget->charge_verification(read.size(),
                                                                 "interval-verification-limit")) {
                                    return;
                                }
                                ++out.full_read_verifications;
                                std::size_t e = 0;
                                if (!verify_run(read, fb, run, soff, cap, &e)) return;
                                ++out.accepted_placements;
                                IntervalPlacement pl;
                                pl.first_block = static_cast<std::uint32_t>(fb);
                                pl.last_block = static_cast<std::uint32_t>(fb + run.size() - 1);
                                pl.alleles = run;
                                pl.offset_in_first = soff;
                                pl.edits = static_cast<std::uint32_t>(e);
                                place[m].push_back(std::move(pl));
                            };
                        std::vector<std::uint32_t> run = acc;
                        grow_right(run);
                        return;
                    }
                    for (std::uint32_t a0 = 0; a0 < geom.alleles[fb - 1].size(); ++a0) {
                        std::vector<std::uint32_t> acc2;
                        acc2.push_back(a0);
                        acc2.insert(acc2.end(), acc.begin(), acc.end());
                        ++out.tuple_expansions;
                        grow_left(fb - 1, acc2, soff + static_cast<long>(seg(fb - 1, a0)));
                    }
                };
            grow_left(first, al, s);
            if (out.work_refused) return out;
        }
        // Identical placements arise from different seeds; a repeat origin is a DIFFERENT start and
        // must survive, so dedup is on the whole descriptor rather than on the cell.
        auto& v = place[m];
        out.placements_before_dedup += v.size();
        std::sort(v.begin(), v.end(), [](const IntervalPlacement& x, const IntervalPlacement& y) {
            return std::tie(x.first_block, x.last_block, x.alleles, x.offset_in_first, x.edits) <
                   std::tie(y.first_block, y.last_block, y.alleles, y.offset_in_first, y.edits);
        });
        v.erase(std::unique(v.begin(), v.end(),
                            [](const IntervalPlacement& x, const IntervalPlacement& y) {
                                return x.first_block == y.first_block &&
                                       x.last_block == y.last_block && x.alleles == y.alleles &&
                                       x.offset_in_first == y.offset_in_first &&
                                       x.edits == y.edits;
                            }),
                v.end());
        out.placements_after_dedup += v.size();
    }
    // ---- THE MATE JOIN ------------------------------------------------------------------------
    // Both FR predicates depend only on rev_end - fwd_start, so this works in RELATIVE coordinates
    // and never needs an absolute prefix -- which is what lets the blocks before the fragment stay
    // free. delta = rev_start - fwd_start = sum_{f0 <= j < g0} seg(j) + o2 - o1, signed by which
    // mate sits first, and every intervening allele and context is counted exactly once.
    std::vector<std::vector<std::array<std::uint32_t, 3>>> sigacc;
    if (want_signatures) sigacc.assign(out.cells, {});
    std::vector<std::vector<std::array<long, 5>>> orgacc;
    if (want_origins) { orgacc.assign(out.cells, {}); out.cell_contrib.assign(out.cells, {}); }
    // Absolute start of a block, for THIS cell. The symbolic path never needs this to decide a
    // state -- the FR predicate uses only the difference -- so it is computed solely to report the
    // origin, and the comparison against the oracle stays non-circular on the join arithmetic.
    const auto prefix_of = [&](const std::vector<std::uint32_t>& cc, std::size_t b) {
        long pfx = static_cast<long>(geom.lflank.size());
        for (std::size_t j = 0; j < b; ++j) pfx += static_cast<long>(seg(j, cc[j]));
        return pfx;
    };
    std::vector<int> pin(k, -1);
    std::vector<std::uint32_t> cell(k, 0);
    // Write one verified state into every cell its constraints cover: pinned dimensions are fixed,
    // free ones range over everything, because the state holds for all of them.
    const auto emit_state = [&](const std::vector<int>& pins, std::uint32_t e1, std::uint32_t e2,
                                long insert, std::size_t fb_blk, long fb_off, std::size_t rb_blk,
                                long rb_off, bool fwd_is_m1) {
        std::vector<std::size_t> freed;
        for (std::size_t j = 0; j < k; ++j) if (pins[j] < 0) freed.push_back(j);
        std::vector<std::uint32_t> idx(freed.size(), 0);
        const double e = origin_contribution(e1, fragment.r1.size(), e2, fragment.r2.size(),
                                             insert, log_eps, log_1meps, ip);
        bool done = false;
        while (!done) {
            for (std::size_t j = 0; j < k; ++j)
                cell[j] = pins[j] >= 0 ? static_cast<std::uint32_t>(pins[j]) : 0;
            for (std::size_t q = 0; q < freed.size(); ++q) cell[freed[q]] = idx[q];
            const std::size_t ci = geom.cell_index(cell);
            out.mass[ci] = log_add(out.mass[ci], e);
            if (want_origins) out.cell_contrib[ci].push_back(e);
            ++out.cell_states[ci];
            if (want_signatures)
                sigacc[ci].push_back({e1, e2, static_cast<std::uint32_t>(insert)});
            if (want_origins) {
                const long fs = prefix_of(cell, fb_blk) + fb_off;
                const long rs = prefix_of(cell, rb_blk) + rb_off;
                const long m1s = fwd_is_m1 ? fs : rs;
                const long m2s = fwd_is_m1 ? rs : fs;
                orgacc[ci].push_back({m1s, fwd_is_m1 ? 1L : 0L, m2s, fwd_is_m1 ? 0L : 1L, insert});
            }
            for (std::size_t q = 0; ; ++q) {
                if (q == freed.size()) { done = true; break; }
                if (++idx[q] < geom.alleles[freed[q]].size()) break;
                idx[q] = 0;
            }
            if (freed.empty()) done = true;
        }
    };
    const auto do_join = [&](int fi, int ri, bool fwd_is_m1) {
        const std::size_t rev_len = fwd_is_m1 ? fragment.r2.size() : fragment.r1.size();
        for (const IntervalPlacement& P : place[fi]) {
            for (const IntervalPlacement& Q : place[ri]) {
                ++out.joined_pairs;
                // 1. INTERSECT the pinned constraints; a conflict is not a state.
                std::fill(pin.begin(), pin.end(), -1);
                bool conflict = false;
                for (std::size_t j = 0; j < P.alleles.size(); ++j)
                    pin[P.first_block + j] = static_cast<int>(P.alleles[j]);
                for (std::size_t j = 0; j < Q.alleles.size() && !conflict; ++j) {
                    const std::size_t b = Q.first_block + j;
                    const int a = static_cast<int>(Q.alleles[j]);
                    if (pin[b] >= 0 && pin[b] != a) conflict = true;
                    else pin[b] = a;
                }
                if (conflict) continue;
                // 2. Only the blocks BETWEEN the two mates affect the insert; unpinned ones there
                //    are enumerated, and everything outside stays free.
                const std::size_t lo_b = std::min(P.first_block, Q.first_block);
                const std::size_t hi_b = std::max(P.first_block, Q.first_block);
                const int sign = Q.first_block >= P.first_block ? 1 : -1;
                std::vector<std::size_t> mid;
                for (std::size_t j = lo_b; j < hi_b; ++j) if (pin[j] < 0) mid.push_back(j);
                std::vector<std::uint32_t> mi(mid.size(), 0);
                bool mdone = false;
                while (!mdone) {
                    std::vector<int> pins = pin;
                    for (std::size_t q = 0; q < mid.size(); ++q)
                        pins[mid[q]] = static_cast<int>(mi[q]);
                    ++out.tuple_expansions;
                    long dsum = 0;
                    for (std::size_t j = lo_b; j < hi_b; ++j)
                        dsum += static_cast<long>(seg(j, static_cast<std::uint32_t>(pins[j])));
                    const long delta = sign * dsum + Q.offset_in_first - P.offset_in_first;
                    // 3. THE SHARED RULE, in relative coordinates.
                    if (valid_fr_coordinates(0, delta + static_cast<long>(rev_len) - 1,
                                             ip.lo, ip.hi)) {
                        ++out.verified_fr_states;
                        if (!mid.empty()) ++out.states_with_free_intermediate;
                        if (fwd_is_m1) ++out.states_orientation_a; else ++out.states_orientation_b;
                        emit_state(pins, fwd_is_m1 ? P.edits : Q.edits,
                                   fwd_is_m1 ? Q.edits : P.edits,
                                   delta + static_cast<long>(rev_len),
                                   P.first_block, P.offset_in_first,
                                   Q.first_block, Q.offset_in_first, fwd_is_m1);
                    }
                    for (std::size_t q = 0; ; ++q) {
                        if (q == mid.size()) { mdone = true; break; }
                        if (++mi[q] < geom.alleles[mid[q]].size()) break;
                        mi[q] = 0;
                    }
                    if (mid.empty()) mdone = true;
                }
            }
        }
    };
    do_join(0, 3, true);    // r1 forward with r2 reverse-complemented
    do_join(2, 1, false);   // r2 forward with r1 reverse-complemented
    for (double mm : out.mass) if (mm != kNegInf) ++out.finite_cells;
    if (want_origins) {
        out.cell_origin.assign(out.cells, std::string());
        for (std::size_t c = 0; c < out.cells; ++c) {
            std::sort(out.cell_contrib[c].begin(), out.cell_contrib[c].end());
            std::sort(orgacc[c].begin(), orgacc[c].end());
            std::string& b = out.cell_origin[c];
            b.resize(orgacc[c].size() * 40);
            for (std::size_t q = 0; q < orgacc[c].size(); ++q)
                std::memcpy(&b[q * 40], orgacc[c][q].data(), 40);
        }
    }
    if (want_signatures) {
        out.cell_signature.assign(out.cells, std::string());
        for (std::size_t c = 0; c < out.cells; ++c) {
            std::sort(sigacc[c].begin(), sigacc[c].end());
            std::string& b = out.cell_signature[c];
            b.resize(sigacc[c].size() * 12);
            for (std::size_t q = 0; q < sigacc[c].size(); ++q)
                std::memcpy(&b[q * 12], sigacc[c][q].data(), 12);
        }
    }
    out.ok = true;
    return out;
}

IntervalOracleCell interval_oracle_cell(const Fragment& fragment, const IntervalGeometry& geom,
                                        const InsertPrior& ip, double max_divergence,
                                        double log_eps, double log_1meps, std::size_t cell) {
    IntervalOracleCell out;
    out.mass = kNegInf;
    const std::size_t d1 = mate_band_edits(max_divergence, fragment.r1.size());
    const std::size_t d2 = mate_band_edits(max_divergence, fragment.r2.size());
    const std::string a1 = reverse_complement(fragment.r1);
    const std::string a2 = reverse_complement(fragment.r2);
    std::vector<std::uint32_t> choice;
    geom.cell_choice(cell, choice);
    std::vector<const std::string*> alle, ctxp;
    for (std::size_t j = 0; j < geom.alleles.size(); ++j)
        alle.push_back(&geom.alleles[j][choice[j]]);
    for (const std::string& cx : geom.contexts) ctxp.push_back(&cx);
    VirtualWindow vw;
    vw.bind_chain(geom.lflank, alle, ctxp, geom.rflank);
    const std::string win = vw.materialize();
    const auto f1 = bounded_mate_placements(fragment.r1, win, d1, nullptr, nullptr);
    const auto v1 = bounded_mate_placements(a1, win, d1, nullptr, nullptr);
    const auto f2 = bounded_mate_placements(fragment.r2, win, d2, nullptr, nullptr);
    const auto v2 = bounded_mate_placements(a2, win, d2, nullptr, nullptr);
    const auto st = enumerate_fragment_states(0, f1, v1, f2, v2, fragment.r1.size(),
                                              fragment.r2.size(), ip.lo, ip.hi);
    out.states = static_cast<std::uint32_t>(st.size());
    std::vector<std::array<long, 5>> og;
    std::vector<std::array<std::uint32_t, 3>> sg;
    for (const FragmentState& z : st) {
        const double c = origin_contribution(z.m1_edits, fragment.r1.size(), z.m2_edits,
                                             fragment.r2.size(), z.insert, log_eps, log_1meps, ip);
        out.mass = log_add(out.mass, c);
        out.contrib.push_back(c);
        og.push_back({z.m1_start, z.m1_fwd ? 1L : 0L, z.m2_start, z.m2_fwd ? 1L : 0L, z.insert});
        sg.push_back({z.m1_edits, z.m2_edits, static_cast<std::uint32_t>(z.insert)});
    }
    std::sort(out.contrib.begin(), out.contrib.end());
    std::sort(og.begin(), og.end());
    std::sort(sg.begin(), sg.end());
    out.origin.resize(og.size() * 40);
    for (std::size_t q = 0; q < og.size(); ++q) std::memcpy(&out.origin[q * 40], og[q].data(), 40);
    out.signature.resize(sg.size() * 12);
    for (std::size_t q = 0; q < sg.size(); ++q) std::memcpy(&out.signature[q * 12], sg[q].data(), 12);
    return out;
}

IntervalEmission interval_emission_oracle(const Fragment& fragment, const IntervalGeometry& geom,
                                          const InsertPrior& ip, double max_divergence,
                                          double log_eps, double log_1meps, double log_p_bg,
                                          bool want_signatures, bool want_origins) {
    IntervalEmission out;
    out.log_p_bg = log_p_bg;
    if (!geom.ok || fragment.r1.empty() || fragment.r2.empty()) return out;
    out.cells = geom.cells();
    out.mass.assign(out.cells, kNegInf);
    out.cell_states.assign(out.cells, 0);
    if (want_signatures) out.cell_signature.assign(out.cells, std::string());
    if (want_origins) {
        out.cell_origin.assign(out.cells, std::string());
        out.cell_contrib.assign(out.cells, {});
    }
    const std::size_t d1 = mate_band_edits(max_divergence, fragment.r1.size());
    const std::size_t d2 = mate_band_edits(max_divergence, fragment.r2.size());
    const std::string a1 = reverse_complement(fragment.r1);
    const std::string a2 = reverse_complement(fragment.r2);
    const double half = std::log(0.5);
    std::vector<std::uint32_t> choice;
    for (std::size_t c = 0; c < out.cells; ++c) {
        geom.cell_choice(c, choice);
        // THE SHARED WINDOW. Building the concatenation here by hand would be the second
        // implementation the segment walk exists to prevent.
        std::vector<const std::string*> alle, ctxp;
        for (std::size_t j = 0; j < geom.alleles.size(); ++j)
            alle.push_back(&geom.alleles[j][choice[j]]);
        for (const std::string& cx : geom.contexts) ctxp.push_back(&cx);
        VirtualWindow vw;
        vw.bind_chain(geom.lflank, alle, ctxp, geom.rflank);
        const std::string win = vw.materialize();
        const auto f1 = bounded_mate_placements(fragment.r1, win, d1, nullptr, nullptr);
        const auto v1 = bounded_mate_placements(a1, win, d1, nullptr, nullptr);
        const auto f2 = bounded_mate_placements(fragment.r2, win, d2, nullptr, nullptr);
        const auto v2 = bounded_mate_placements(a2, win, d2, nullptr, nullptr);
        const auto st = enumerate_fragment_states(0, f1, v1, f2, v2, fragment.r1.size(),
                                                  fragment.r2.size(), ip.lo, ip.hi);
        out.cell_states[c] = static_cast<std::uint32_t>(st.size());
        for (const FragmentState& z : st) {
            if (z.m1_fwd) ++out.states_orientation_a; else ++out.states_orientation_b;
        }
        double m = kNegInf;
        std::vector<std::array<std::uint32_t, 3>> sig;
        for (const FragmentState& z : st) {
            const double contrib = origin_contribution(z.m1_edits, fragment.r1.size(), z.m2_edits,
                                                       fragment.r2.size(), z.insert, log_eps,
                                                       log_1meps, ip);
            m = log_add(m, contrib);
            if (want_origins) out.cell_contrib[c].push_back(contrib);
            if (want_signatures) {
                sig.push_back({z.m1_edits, z.m2_edits, static_cast<std::uint32_t>(z.insert)});
            }
        }
        if (want_origins) std::sort(out.cell_contrib[c].begin(), out.cell_contrib[c].end());
        if (want_origins) {
            // The FULL origin, from the oracle's own absolute coordinates.
            std::vector<std::array<long, 5>> og;
            for (const FragmentState& z : st) {
                og.push_back({z.m1_start, z.m1_fwd ? 1L : 0L, z.m2_start, z.m2_fwd ? 1L : 0L,
                              z.insert});
            }
            std::sort(og.begin(), og.end());
            std::string& b = out.cell_origin[c];
            b.resize(og.size() * 40);
            for (std::size_t q = 0; q < og.size(); ++q) std::memcpy(&b[q * 40], og[q].data(), 40);
        }
        out.mass[c] = m;
        if (m != kNegInf) ++out.finite_cells;
        out.verified_fr_states += st.size();
        if (want_signatures) {
            std::sort(sig.begin(), sig.end());
            std::string& b = out.cell_signature[c];
            b.resize(sig.size() * 12);
            for (std::size_t q = 0; q < sig.size(); ++q) std::memcpy(&b[q * 12], sig[q].data(), 12);
        }
    }
    out.ok = true;
    return out;
}

IntervalSeedIndex build_interval_seed_index(const IntervalGeometry& geom, std::size_t piece) {
    IntervalSeedIndex ix;
    ix.piece = piece;
    if (!geom.ok || piece == 0 || piece > 32) return ix;
    const std::size_t k = geom.alleles.size();
    // THE COMPLETENESS PREDICATE FIRST, so a locus that cannot be indexed completely is known
    // before anything is built rather than after a placement has already been missed.
    ix.shortest_allele = SIZE_MAX;
    for (std::size_t j = 0; j < k; ++j) {
        for (const std::string& a : geom.alleles[j]) {
            if (a.size() < ix.shortest_allele) {
                ix.shortest_allele = a.size();
                ix.shortest_at_block = static_cast<std::uint32_t>(j);
            }
        }
    }
    if (ix.shortest_allele == SIZE_MAX) ix.shortest_allele = 0;
    // An allele shorter than piece - 1 could sit wholly inside one seed, making a three-allele
    // span reachable. Those shapes are not indexed, so the index is incomplete.
    ix.complete = ix.shortest_allele + 1 >= piece;
    // AND A CONTEXT AT LEAST AS LONG AS A PIECE can hold a seed WHOLLY INSIDE IT. Such a seed
    // constrains no allele and is indexed nowhere here, so a read whose every matching piece falls
    // in a context would be missed entirely -- pigeonhole guarantees a piece matches, not that the
    // matching piece is one this index carries. C4's contexts are empty, so this never bites there;
    // it is a predicate about the locus and the generic path must refuse rather than assume.
    // Every INVARIANT segment, flanks included: they are indexed nowhere either, so a long flank
    // hides the same gap as a long context.
    ix.longest_context = 0;
    for (const std::string& c : geom.contexts)
        ix.longest_context = std::max(ix.longest_context, c.size());
    ix.longest_context = std::max(ix.longest_context, geom.lflank.size());
    ix.longest_context = std::max(ix.longest_context, geom.rflank.size());
    if (ix.longest_context + 1 > piece) ix.complete = false;
    for (std::size_t j = 0; j < k; ++j) {
        for (std::uint32_t a = 0; a < geom.alleles[j].size(); ++a) {
            const std::string& A = geom.alleles[j][a];
            for (std::size_t at = 0; at + piece <= A.size(); ++at) {
                std::uint64_t code = 0;
                if (!encode_piece(A, at, piece, code)) continue;
                IntervalSeedHit h;
                h.block = static_cast<std::uint32_t>(j);
                h.allele = a;
                h.offset = static_cast<std::uint32_t>(at);
                ix.inside[code].push_back(h);
            }
        }
    }
    // BOUNDARIES, one per adjacent pair, keeping the two alleles CORRELATED in a single hit.
    for (std::size_t j = 0; j + 1 < k; ++j) {
        const std::string& C = geom.contexts[j];
        for (std::uint32_t a = 0; a < geom.alleles[j].size(); ++a) {
            const std::string& A = geom.alleles[j][a];
            const std::size_t ta = A.size() < piece - 1 ? A.size() : piece - 1;
            const std::string tail = A.substr(A.size() - ta);
            // A SEED THAT LEAVES THE ALLELE BUT STOPS INSIDE THE CONTEXT constrains block j ONLY.
            // Recording it as spanning the boundary would pin the NEXT allele it never touched --
            // once per candidate, so the same physical placement arrives with several different
            // run lengths and each emits its own state. That is a double count of exactly one
            // state, and it shows up as a log 2 mass difference.
            {
                const std::string ac = tail + C;
                for (std::size_t at = 0; at + piece <= ac.size(); ++at) {
                    if (at + piece <= ta) continue;          // wholly inside A: already in `inside`
                    std::uint64_t code = 0;
                    if (!encode_piece(ac, at, piece, code)) continue;
                    IntervalSeedHit h;
                    h.block = static_cast<std::uint32_t>(j);
                    h.allele = a;
                    h.offset = static_cast<std::uint32_t>(A.size() - ta + at);
                    h.spans_boundary = false;
                    ix.boundary[code].push_back(h);
                }
            }
            for (std::uint32_t b = 0; b < geom.alleles[j + 1].size(); ++b) {
                const std::string& B = geom.alleles[j + 1][b];
                const std::size_t hb = B.size() < piece - 1 ? B.size() : piece - 1;
                const std::string joined = tail + C + B.substr(0, hb);
                for (std::size_t at = 0; at + piece <= joined.size(); ++at) {
                    // Only seeds that genuinely REACH THE NEXT ALLELE pin it. Ones lying wholly in
                    // either allele are in `inside`; ones stopping in the context are above.
                    if (at + piece <= ta + C.size()) continue;
                    if (at >= ta + C.size()) continue;
                    std::uint64_t code = 0;
                    if (!encode_piece(joined, at, piece, code)) continue;
                    IntervalSeedHit h;
                    h.block = static_cast<std::uint32_t>(j);
                    h.allele = a;
                    // Offset back into block j's own coordinates: the tail began at |A| - ta.
                    h.offset = static_cast<std::uint32_t>(A.size() - ta + at);
                    h.next_allele = b;
                    h.spans_boundary = true;
                    ix.boundary[code].push_back(h);
                }
            }
        }
    }
    ix.ok = true;
    return ix;
}

IntervalGrouping build_interval_grouping(
    const IntervalGeometry& geom,
    const std::vector<std::vector<std::string>>& per_fragment_signatures) {
    IntervalGrouping G;
    const std::size_t k = geom.alleles.size();
    const std::size_t ncell = geom.cells();
    for (const auto& v : per_fragment_signatures) {
        if (v.size() != ncell) {
            G.refusal = "a fragment's signature vector does not span the cell product";
            return G;
        }
    }
    // The JOINT signature per cell, length-prefixed per fragment so two different splits cannot
    // alias, then interned to a dense id.
    std::unordered_map<std::string, std::uint32_t> intern;
    G.cell_signature_id.assign(ncell, 0);
    {
        std::string key;
        for (std::size_t c = 0; c < ncell; ++c) {
            key.clear();
            for (const auto& v : per_fragment_signatures) {
                const std::uint32_t n = static_cast<std::uint32_t>(v[c].size());
                key.append(reinterpret_cast<const char*>(&n), 4);
                key.append(v[c]);
            }
            G.cell_signature_id[c] = intern.emplace(key,
                static_cast<std::uint32_t>(intern.size())).first->second;
        }
    }
    G.distinct_cell_signatures = intern.size();
    // Per block: two alleles are equivalent when their whole SLICE of signature ids agrees. The
    // slice is taken over every combination of the other blocks, so this is equivalence across the
    // entire remaining product, not agreement at a sampled point.
    G.allele_class.resize(k);
    G.classes_per_block.assign(k, 0);
    std::vector<std::size_t> stride(k, 1);
    for (std::size_t j = k; j-- > 0;)
        stride[j] = (j + 1 < k) ? stride[j + 1] * geom.alleles[j + 1].size() : 1;
    for (std::size_t j = 0; j < k; ++j) {
        const std::size_t nj = geom.alleles[j].size();
        std::unordered_map<std::string, std::uint32_t> cls;
        G.allele_class[j].assign(nj, 0);
        for (std::uint32_t a = 0; a < nj; ++a) {
            std::string slice;
            slice.reserve((ncell / nj) * 4);
            for (std::size_t c = 0; c < ncell; ++c) {
                if ((c / stride[j]) % nj != a) continue;
                const std::uint32_t id = G.cell_signature_id[c];
                slice.append(reinterpret_cast<const char*>(&id), 4);
            }
            G.allele_class[j][a] = cls.emplace(slice,
                static_cast<std::uint32_t>(cls.size())).first->second;
        }
        G.classes_per_block[j] = cls.size();
    }
    // ---- THE JOINT EQUALITY CHECK -------------------------------------------------------------
    // Every original cell must carry the IDENTICAL signature to the representative cell of its
    // class tuple -- not merely agree marginally at each block. One representative per class tuple
    // is taken from the first allele of each class, and every cell is compared against it.
    {
        std::vector<std::vector<std::uint32_t>> first_of_class(k);
        for (std::size_t j = 0; j < k; ++j) {
            first_of_class[j].assign(G.classes_per_block[j],
                                     std::numeric_limits<std::uint32_t>::max());
            for (std::uint32_t a = 0; a < geom.alleles[j].size(); ++a) {
                const std::uint32_t cl = G.allele_class[j][a];
                if (first_of_class[j][cl] == std::numeric_limits<std::uint32_t>::max())
                    first_of_class[j][cl] = a;
            }
        }
        std::vector<std::uint32_t> ch, rep(k, 0);
        for (std::size_t c = 0; c < ncell; ++c) {
            geom.cell_choice(c, ch);
            for (std::size_t j = 0; j < k; ++j)
                rep[j] = first_of_class[j][G.allele_class[j][ch[j]]];
            ++G.cells_checked;
            if (G.cell_signature_id[c] != G.cell_signature_id[geom.cell_index(rep)])
                ++G.cells_disagreeing_with_representative;
        }
        G.joint_equality_verified = G.cells_disagreeing_with_representative == 0;
    }
    if (!G.joint_equality_verified) {
        G.refusal = "cells sharing a class tuple do not share a signature";
        return G;
    }
    const auto pairs = [](std::size_t n) { return n * (n + 1) / 2; };
    G.content_classes_raw = 1;
    G.content_classes_grouped = 1;
    G.factor_lookup_configurations = 1;
    for (std::size_t j = 0; j < k; ++j) {
        G.content_classes_raw *= pairs(geom.alleles[j].size());
        G.content_classes_grouped *= pairs(G.classes_per_block[j]);
        G.factor_lookup_configurations *= G.classes_per_block[j] * G.classes_per_block[j];
    }
    // Classes with at most one heterozygous block have a single biological phase and are EXACTLY
    // neutral, so they need never be stored. Counted here so the prediction below is honest about
    // how many classes actually carry a value.
    {
        std::size_t m0 = 1, m1 = 0;
        for (std::size_t j = 0; j < k; ++j) {
            const std::size_t R = G.classes_per_block[j];
            const std::size_t hom = R, het = R * (R - 1) / 2;
            m1 = m1 * hom + m0 * het;
            m0 *= hom;
        }
        G.classes_m_le_1 = m0 + m1;
        G.ordered_in_m_le_1 = m0 + 2 * m1;
    }
    // PREDICTED BYTES, before a single entry is allocated -- and it must count STORED PHASE
    // VALUES, not classes. A non-neutral class carries one value per biological phase, 2^(m-1) of
    // them, so "one double per class" understates it by the mean phase count. The totals follow
    // from sums already known: every class contributes 2^m ordered configurations and their sum is
    // exactly prod(R_j^2); the m <= 1 classes contribute m0 + 2*m1 of those and store nothing.
    G.stored_ordered_values = G.factor_lookup_configurations > G.ordered_in_m_le_1
                            ? G.factor_lookup_configurations - G.ordered_in_m_le_1 : 0;
    G.stored_canonical_values = G.stored_ordered_values / 2;   // two ordered per biological phase
    // Canonical storage: one key per class that carries a value, plus one double per phase.
    G.predicted_bytes =
        (G.content_classes_grouped - G.classes_m_le_1) *
            (sizeof(std::uint64_t) + 2 * sizeof(void*)) +
        G.stored_canonical_values * sizeof(double);
    G.ok = true;
    return G;
}

namespace {

// The index of the unordered pair (lo <= hi) among the C(R+1, 2) pairs of R classes, in a fixed
// enumeration: (0,0), (0,1), ... (0,R-1), (1,1), ...
inline std::size_t pair_ordinal(std::size_t lo, std::size_t hi, std::size_t R) {
    return lo * R - lo * (lo - 1) / 2 + (hi - lo);
}

}  // namespace

double IntervalFactorTable::log_psi(const std::vector<std::uint32_t>& hap1,
                                    const std::vector<std::uint32_t>& hap2) const {
    if (!ok || class_offset.empty()) return 0.0;
    const std::size_t k = allele_class.size();
    std::vector<std::uint32_t> c1(k), c2(k);
    for (std::size_t j = 0; j < k; ++j) {
        c1[j] = allele_class[j][hap1[j]];
        c2[j] = allele_class[j][hap2[j]];
    }
    return log_psi_classes(c1, c2);
}

// The two factor kinds behind one call. The pairwise arithmetic is the SparsePhaseClass
// convention from chain_kernel: straight pairs (amin,bmin) with (amax,bmax), crossed pairs
// (amin,bmax) with (amax,bmin). With a homozygous endpoint the two coincide and either is right.
double HybridHigherFactor::log_psi_classes(const std::vector<std::uint32_t>& c1,
                                           const std::vector<std::uint32_t>& c2) const {
    if (table != nullptr) return table->log_psi_classes(c1, c2);
    if (pairwise == nullptr || !pairwise->active || c1.size() != 2 || c2.size() != 2) return 0.0;
    const std::uint32_t a1 = c1[0], b1 = c1[1], a2 = c2[0], b2 = c2[1];
    const std::uint32_t amin = std::min(a1, a2), amax = std::max(a1, a2);
    const std::uint32_t bmin = std::min(b1, b2), bmax = std::max(b1, b2);
    const bool crossed = (a1 < a2 && b1 > b2) || (a1 > a2 && b1 < b2);
    for (const SparsePhaseClass& k : pairwise->classes) {
        if (k.amin != amin || k.amax != amax || k.bmin != bmin || k.bmax != bmax) continue;
        return std::log1p(crossed ? k.crossed_m1 : k.straight_m1);
    }
    return 0.0;   // absent means exactly neutral
}

double IntervalFactorTable::log_psi_classes(const std::vector<std::uint32_t>& cls1,
                                            const std::vector<std::uint32_t>& cls2) const {
    if (!ok || class_offset.empty()) return 0.0;
    const std::size_t k = allele_class.size();
    std::size_t ordinal = 0, m = 0;
    std::uint32_t bits = 0;
    for (std::size_t j = 0; j < k; ++j) {
        const std::uint32_t c1 = cls1[j];
        const std::uint32_t c2 = cls2[j];
        const std::size_t R = classes_per_block[j];
        const std::size_t lo = std::min(c1, c2), hi = std::max(c1, c2);
        ordinal = ordinal * (R * (R + 1) / 2) + pair_ordinal(lo, hi, R);
        // EFFECTIVE heterozygosity: a block is a phase dimension only when the two homologues
        // carry DIFFERENT CLASSES. Two distinct alleles inside one class are biologically
        // heterozygous and yet indistinguishable here, so swapping them stays neutral.
        if (c1 != c2) {
            if (c1 != lo) bits |= (1u << m);
            ++m;
        }
    }
    if (m <= 1) return 0.0;   // one biological phase: no preference is expressible
    // Canonicalise under global swap, which flips every bit.
    const std::uint32_t mask = (1u << m) - 1u;
    if (bits & 1u) bits = (~bits) & mask;
    const std::size_t phase = bits >> 1;
    const std::size_t off = class_offset[ordinal];
    if (off + phase >= class_offset[ordinal + 1]) return 0.0;
    return phase_value[off + phase];
}

IntervalFactorTable build_interval_factor(const IntervalGeometry& geom,
                                          const IntervalGrouping& grouping,
                                          const std::vector<IntervalEmission>& emissions,
                                          double lambda, double log_mix, double log_bg_weight,
                                          const std::vector<std::size_t>* context_positions,
                                          double* worst_context_delta,
                                          const std::vector<char>* cell_on_panel,
                                          double* observed_context_delta,
                                          double* observed_context_delta_panel) {
    IntervalFactorTable T;
    const std::size_t k = geom.alleles.size();
    if (!grouping.ok) { T.refusal = "grouping is not usable"; return T; }
    T.allele_class = grouping.allele_class;
    T.classes_per_block = grouping.classes_per_block;
    // One representative ALLELE per class, for reading the emission. Any member gives the same
    // signature, which is exactly what the joint-equality check established.
    std::vector<std::vector<std::uint32_t>> rep(k);
    for (std::size_t j = 0; j < k; ++j) {
        rep[j].assign(grouping.classes_per_block[j], std::numeric_limits<std::uint32_t>::max());
        for (std::uint32_t a = 0; a < geom.alleles[j].size(); ++a) {
            const std::uint32_t c = grouping.allele_class[j][a];
            if (rep[j][c] == std::numeric_limits<std::uint32_t>::max()) rep[j][c] = a;
        }
    }
    std::size_t n_classes = 1;
    std::vector<std::size_t> pairs_per_block(k);
    for (std::size_t j = 0; j < k; ++j) {
        const std::size_t R = grouping.classes_per_block[j];
        pairs_per_block[j] = R * (R + 1) / 2;
        n_classes *= pairs_per_block[j];
    }
    T.class_offset.assign(n_classes + 1, 0);
    // PREDICT BEFORE ALLOCATING, and predict the two accounts separately.
    T.predicted_payload_bytes = grouping.stored_canonical_values * sizeof(double);
    T.predicted_total_bytes = T.predicted_payload_bytes +
                              (n_classes + 1) * sizeof(std::size_t);
    for (std::size_t j = 0; j < k; ++j)
        T.predicted_total_bytes += geom.alleles[j].size() * sizeof(std::uint32_t);
    const double log_lambda = std::log(lambda);
    // THE SHARED COMBINER. Production uses it; the oracle must NOT, or a defect here would move
    // both sides of the comparison together.
    const auto mix = [&](double p, double q, double log_p_bg) {
        double sig = kNegInf;
        if (p != kNegInf) sig = p;
        if (q != kNegInf) sig = (sig == kNegInf) ? q : log_add(sig, q);
        if (sig != kNegInf) sig += log_mix + log_lambda;
        const double bg = log_bg_weight + log_p_bg;
        return (sig == kNegInf) ? bg : log_add(sig, bg);
    };

    // ---- THE CONTEXT CERTIFICATE ------------------------------------------------------------
    // How far can S move when ONLY the normalisation-context blocks change? Cells are grouped by
    // their choices at the evidence positions; within a group only context alleles differ. For
    // each fragment the extreme masses in its group bound that fragment's contribution swing
    // through the REAL mixture -- mix() with its own background floor -- and the sum over
    // fragments bounds the swing in S. This is a certificate, not an estimate: it is the worst
    // case over every group and every fragment.
    if (context_positions != nullptr && worst_context_delta != nullptr) {
        *worst_context_delta = 0.0;
        std::map<std::vector<std::uint32_t>, std::vector<std::size_t>> groups;
        std::vector<std::uint32_t> ch;
        for (std::size_t c = 0; c < geom.cells(); ++c) {
            geom.cell_choice(c, ch);
            std::vector<std::uint32_t> key;
            for (std::size_t j2 = 0; j2 < ch.size(); ++j2)
                if (std::find(context_positions->begin(), context_positions->end(), j2) ==
                    context_positions->end())
                    key.push_back(ch[j2]);
            groups[key].push_back(c);
        }
        // OBSERVED, not bounded: for each group, the actual S at every context choice, summed
        // over ALL fragments at the SAME configuration. The spread of those is a difference the
        // model can really exhibit. The per-fragment sum below is an upper bound instead: each
        // fragment is allowed its own worst group, which no single configuration need realise.
        if (observed_context_delta != nullptr) {
            double obs = 0.0, obs_panel = 0.0;
            for (const auto& kv : groups) {
                if (kv.second.size() < 2) continue;
                double lo3 = std::numeric_limits<double>::infinity(), hi3 = -lo3;
                double lo3p = lo3, hi3p = hi3;
                for (std::size_t c : kv.second) {
                    double acc2 = 0.0;
                    for (const IntervalEmission& E : emissions)
                        acc2 += mix(E.mass[c], E.mass[c], E.log_p_bg);
                    lo3 = std::min(lo3, acc2); hi3 = std::max(hi3, acc2);
                    if (cell_on_panel != nullptr && c < cell_on_panel->size() &&
                        (*cell_on_panel)[c]) {
                        lo3p = std::min(lo3p, acc2); hi3p = std::max(hi3p, acc2);
                    }
                }
                if (hi3 > lo3) obs = std::max(obs, hi3 - lo3);
                if (hi3p > lo3p) obs_panel = std::max(obs_panel, hi3p - lo3p);
            }
            *observed_context_delta = obs;
            if (observed_context_delta_panel != nullptr) *observed_context_delta_panel = obs_panel;
        }
        // THE PER-FRAGMENT-EXTREME SUM IS REMOVED. It summed each fragment's worst spread over
        // ITS OWN best group, which no single configuration need realise, and it was reported and
        // gated on as if it were an upper bound. On C4 it read 32.02 nats while the ACTUAL
        // group-wise difference was 333.61 -- a "bound" an order of magnitude below the thing it
        // claimed to bound. The observed group-wise maximum above is the max-norm difference over
        // the certification domain, and that is what the Lipschitz argument needs.
        if (worst_context_delta != nullptr && observed_context_delta != nullptr)
            *worst_context_delta = *observed_context_delta;
    }
    // Walk every content class in the same mixed-radix order log_psi uses.
    std::vector<std::size_t> pidx(k, 0);
    std::vector<std::uint32_t> lo(k, 0), hi(k, 0), h1(k, 0), h2(k, 0);
    std::vector<double> S;
    std::size_t cursor = 0;
    bool done = n_classes == 0;
    std::size_t ordinal = 0;
    while (!done) {
        // Decode this class's per-block unordered class pair from its pair ordinals.
        std::size_t m = 0;
        std::vector<std::size_t> hetj;
        for (std::size_t j = 0; j < k; ++j) {
            const std::size_t R = grouping.classes_per_block[j];
            std::size_t p = pidx[j], a = 0;
            while (p >= R - a) { p -= (R - a); ++a; }
            lo[j] = static_cast<std::uint32_t>(a);
            hi[j] = static_cast<std::uint32_t>(a + p);
            if (lo[j] != hi[j]) { hetj.push_back(j); ++m; }
        }
        T.class_offset[ordinal] = cursor;
        if (m >= 2) {
            const ContentClassNorm N = content_class_norm(m);
            const std::size_t nph = N.biological_phases;
            // ---- GLOBAL-SWAP COMPLETENESS, checked on the ORDERED configurations -------------
            // Every ordered configuration must have its partner under swapping both homologues,
            // and the two must carry the IDENTICAL raw score -- that equality is what licenses
            // folding 2^m ordered configurations onto 2^(m-1) biological phases at all. Averaging
            // over a class whose partners disagree would weight one phase more than the other, so
            // a disagreement is a REFUSAL and not a smaller factor.
            const std::size_t nord = N.ordered_configs;
            std::vector<double> So(nord, 0.0);
            for (std::size_t bits = 0; bits < nord; ++bits) {
                for (std::size_t j = 0; j < k; ++j) { h1[j] = rep[j][lo[j]]; h2[j] = rep[j][hi[j]]; }
                for (std::size_t b = 0; b < m; ++b) {
                    if (bits & (std::size_t(1) << b)) {
                        const std::size_t j = hetj[b];
                        h1[j] = rep[j][hi[j]]; h2[j] = rep[j][lo[j]];
                    }
                }
                const std::size_t d1 = geom.cell_index(h1), d2 = geom.cell_index(h2);
                double acc = 0.0;
                for (const IntervalEmission& E : emissions)
                    acc += mix(E.mass[d1], E.mass[d2], E.log_p_bg);
                So[bits] = acc;
            }
            const std::size_t omask = nord - 1;
            for (std::size_t bits = 0; bits < nord; ++bits) {
                const std::size_t partner = (~bits) & omask;
                if (std::abs(So[bits] - So[partner]) > 1e-9) ++T.swap_partner_failures;
            }
            if (T.swap_partner_failures != 0) {
                T.refusal = "global-swap partners disagree on the raw score";
                T.phase_value.clear(); T.class_offset.clear();
                return T;
            }
            S.assign(nph, 0.0);
            for (std::size_t ph = 0; ph < nph; ++ph) {
                // Canonical phase ph: bit 0 is always 0, the remaining bits are ph.
                const std::uint32_t bits = static_cast<std::uint32_t>(ph << 1);
                for (std::size_t j = 0; j < k; ++j) { h1[j] = rep[j][lo[j]]; h2[j] = rep[j][hi[j]]; }
                for (std::size_t b = 0; b < m; ++b) {
                    if (bits & (1u << b)) {
                        const std::size_t j = hetj[b];
                        h1[j] = rep[j][hi[j]]; h2[j] = rep[j][lo[j]];
                    }
                }
                const std::size_t c1 = geom.cell_index(h1), c2 = geom.cell_index(h2);
                // THE DIPLOID FORMULA, once per fragment: the two haplotype masses combined, the
                // background mixed in, and only THEN summed across fragments. Summing haploid log
                // masses, or centring per fragment, are different models.
                double acc = 0.0;
                for (const IntervalEmission& E : emissions)
                    acc += mix(E.mass[c1], E.mass[c2], E.log_p_bg);
                S[ph] = acc;
            }
            double mx = -std::numeric_limits<double>::infinity();
            for (double x : S) mx = std::max(mx, x);
            double sum = 0.0;
            for (double x : S) sum += std::exp(x - mx);
            const double Z = mx + std::log(sum);
            double top = -std::numeric_limits<double>::infinity();
            for (std::size_t ph = 0; ph < nph; ++ph) {
                const double v = S[ph] - Z + N.centring_canonical;
                T.phase_value.push_back(v);
                top = std::max(top, v);
            }
            T.max_log_psi = std::max(T.max_log_psi, top);
            const double slack = N.max_log_psi_bound - top;
            if (T.classes_stored == 0 || slack < T.worst_bound_slack) T.worst_bound_slack = slack;
            if (slack < -1e-9) {
                T.refusal = "a class exceeded its own log psi bound";
                T.phase_value.clear(); T.class_offset.clear();   // no indexable partial table
                return T;
            }
            // AND THE ARRANGEMENT COUNT ITSELF. A class that produced fewer phase values than
            // its heterozygosity demands is INCOMPLETE, and an incomplete factor must refuse
            // rather than present a smaller one -- "zero classes stored" is indistinguishable
            // from a legitimately neutral factor, which is exactly the ambiguity to avoid.
            if (T.phase_value.size() != cursor + nph) {
                T.refusal = "missing arrangement: a class produced " +
                            std::to_string(T.phase_value.size() - cursor) + " of " +
                            std::to_string(nph) + " biological phases";
                T.phase_value.clear(); T.class_offset.clear();
                return T;
            }
            cursor += nph;
            ++T.classes_stored;
            T.phase_values_stored += nph;
        } else {
            ++T.classes_neutral;
        }
        ++ordinal;
        for (std::size_t j = k; ; ) {
            if (j == 0) { done = true; break; }
            --j;
            if (++pidx[j] < pairs_per_block[j]) break;
            pidx[j] = 0;
            if (j == 0) { done = true; break; }
        }
    }
    T.class_offset[n_classes] = cursor;
    T.canonical_payload_bytes = T.phase_value.size() * sizeof(double);
    T.total_factor_bytes = T.canonical_payload_bytes +
                           T.class_offset.size() * sizeof(std::size_t);
    for (std::size_t j = 0; j < k; ++j)
        T.total_factor_bytes += T.allele_class[j].size() * sizeof(std::uint32_t);
    T.ok = true;
    return T;
}

ContentClassNorm content_class_norm(std::size_t m_het) {
    ContentClassNorm n;
    n.m_het = m_het;
    if (m_het > 30) { n.refusal = "too many heterozygous blocks to enumerate"; return n; }
    n.ordered_configs = static_cast<std::size_t>(1) << m_het;
    // m = 0 has ONE configuration which is its own global-swap partner, and one biological phase.
    n.biological_phases = m_het == 0 ? 1 : (static_cast<std::size_t>(1) << (m_het - 1));
    n.centring_ordered = std::log(static_cast<double>(n.ordered_configs));
    n.centring_canonical = std::log(static_cast<double>(n.biological_phases));
    // The maximum is log(number of biological phases) in BOTH forms: the ordered form's larger
    // constant is exactly cancelled by its two equal representatives inside the logsumexp.
    n.max_log_psi_bound = std::log(static_cast<double>(n.biological_phases));
    n.ok = true;
    return n;
}

SignatureMatrix build_signature_matrix(const std::vector<std::string>& cell_signatures,
                                       std::size_t n_a, std::size_t n_b, const char* model_tag) {
    SignatureMatrix M;
    M.n_a = n_a; M.n_b = n_b;
    M.model_tag = model_tag != nullptr ? model_tag : "";
    if (M.model_tag != kSignatureModelTag) {
        // The signature's exactness is a property of the emission model. A matrix built for a
        // different one is not a smaller answer, it is an unsound cache.
        M.refusal = "signature model tag '" + M.model_tag + "' is not '" +
                    std::string(kSignatureModelTag) + "'";
        return M;
    }
    if (cell_signatures.size() != n_a * n_b) {
        M.refusal = "cell signature count does not match the allele product";
        return M;
    }
    std::unordered_map<std::string, std::uint32_t> sid;
    M.cell_to_signature.resize(n_a * n_b);
    for (std::size_t k = 0; k < cell_signatures.size(); ++k) {
        M.cell_to_signature[k] = sid.emplace(cell_signatures[k],
            static_cast<std::uint32_t>(sid.size())).first->second;
    }
    M.n_signatures = sid.size();
    // ROW AND COLUMN CLASSES, over the COMPLETE vector of signature ids -- a row class is not
    // "these alleles agree somewhere", it is "these alleles agree at every B allele".
    const auto classify = [](std::size_t outer, std::size_t inner,
                             const std::vector<std::uint32_t>& cells, bool by_row,
                             std::size_t n_bb,
                             std::vector<std::uint32_t>& cls,
                             std::vector<std::vector<std::uint32_t>>& members,
                             std::vector<std::uint32_t>& lo,
                             std::vector<std::uint32_t>& hi) {
        std::unordered_map<std::string, std::uint32_t> id;
        cls.assign(outer, 0);
        for (std::size_t o = 0; o < outer; ++o) {
            std::string key(inner * 4, '\0');
            for (std::size_t i = 0; i < inner; ++i) {
                const std::uint32_t v = by_row ? cells[o * n_bb + i] : cells[i * n_bb + o];
                std::memcpy(&key[i * 4], &v, 4);
            }
            cls[o] = id.emplace(key, static_cast<std::uint32_t>(id.size())).first->second;
        }
        members.assign(id.size(), {});
        lo.assign(id.size(), std::numeric_limits<std::uint32_t>::max());
        hi.assign(id.size(), 0);
        for (std::size_t o = 0; o < outer; ++o) {
            const std::uint32_t c = cls[o];
            members[c].push_back(static_cast<std::uint32_t>(o));
            lo[c] = std::min(lo[c], static_cast<std::uint32_t>(o));
            hi[c] = std::max(hi[c], static_cast<std::uint32_t>(o));
        }
    };
    classify(n_a, n_b, M.cell_to_signature, true, n_b, M.row_class, M.row_members,
             M.row_min, M.row_max);
    classify(n_b, n_a, M.cell_to_signature, false, n_b, M.col_class, M.col_members,
             M.col_min, M.col_max);
    if (!M.validate()) {
        M.refusal = "a cell does not reconstruct its signature from its row and column classes";
        return M;
    }
    M.ok = true;
    return M;
}

bool SignatureMatrix::validate() const {
    if (cell_to_signature.size() != n_a * n_b) return false;
    if (row_class.size() != n_a || col_class.size() != n_b) return false;
    // MEMBERSHIP IS A PARTITION: every allele in exactly one class, every class non-empty.
    std::size_t seen = 0;
    for (const auto& m : row_members) { if (m.empty()) return false; seen += m.size(); }
    if (seen != n_a) return false;
    seen = 0;
    for (const auto& m : col_members) { if (m.empty()) return false; seen += m.size(); }
    if (seen != n_b) return false;
    // THE RECONSTRUCTION: members of one row class must agree at EVERY column, and members of one
    // column class at every row. If that fails the classes are not equivalence classes and every
    // count derived from them is meaningless.
    for (const auto& m : row_members) {
        for (std::size_t b = 0; b < n_b; ++b) {
            const std::uint32_t want = cell_to_signature[m[0] * n_b + b];
            for (std::uint32_t a : m)
                if (cell_to_signature[a * n_b + b] != want) return false;
        }
    }
    for (const auto& m : col_members) {
        for (std::size_t a = 0; a < n_a; ++a) {
            const std::uint32_t want = cell_to_signature[a * n_b + m[0]];
            for (std::uint32_t b : m)
                if (cell_to_signature[a * n_b + b] != want) return false;
        }
    }
    return true;
}

bool HybridWorkBudget::charge_cells(std::uint64_t n, const char* why) {
    if (exhausted) return false;
    proposed_cells += n;
    if (max_proposed_cells != 0 && proposed_cells > max_proposed_cells) {
        exhausted = true; reason = why;
        return false;
    }
    return true;
}

bool HybridWorkBudget::charge_verification(std::size_t read_len, const char* why) {
    if (exhausted) return false;
    ++full_read_verifications;
    bases_compared_upper_bound += read_len;
    if (max_full_read_verifications != 0 &&
        full_read_verifications > max_full_read_verifications) {
        exhausted = true; reason = why;
        return false;
    }
    return true;
}

LinkageEmission linkage_emission_supported(const Fragment& fragment, const LinkageGeometry& geom,
                                           const InsertPrior& ip, double max_divergence,
                                           double log_eps, double log_1meps, double log_p_bg,
                                           AlleleProductSupport* out_support,
                                           const AlleleProductIndex* index,
                                           HybridWorkBudget* budget,
                                           std::vector<std::string>* out_cell_signatures) {
    LinkageEmission out;
    out.log_p_bg = log_p_bg;
    if (!geom.ok || fragment.r1.empty() || fragment.r2.empty()) return out;
    out.n_a = geom.alleles_a.size();
    out.n_b = geom.alleles_b.size();
    // EVERY CELL EXISTS, including the empty ones: a pair with no placement is -inf, which is a
    // value the consumer needs, not a cell to omit.
    out.mass.assign(out.n_a * out.n_b, kNegInf);
    out.cell_states.assign(out.n_a * out.n_b, 0);
    const std::size_t d1 = mate_band_edits(max_divergence, fragment.r1.size());
    const std::size_t d2 = mate_band_edits(max_divergence, fragment.r2.size());
    const std::string a1 = reverse_complement(fragment.r1), a2 = reverse_complement(fragment.r2);
    AlleleProductSupport sup = propose_allele_pairs(fragment, geom, ip, max_divergence, index);
    // THE REFUSAL PATH, taken before any enumeration begins.
    const auto refuse = [&](const char* why) {
        out.mass.assign(out.n_a * out.n_b, kNegInf);
        out.cell_states.assign(out.n_a * out.n_b, 0);
        out.work_refused = true;
        out.work_refusal = why;
        out.ok = false;
        if (out_support != nullptr) *out_support = sup;
        return out;
    };
    if (budget != nullptr && budget->exhausted) return refuse(budget->reason.c_str());
    if (sup.exhaustive_fallback) {
        // THE FALLBACK MUST FIT BEFORE IT STARTS. An incomplete index or a non-ACGT read used to
        // drop straight into dense score_window() enumeration -- materialising one window per
        // allele pair -- with the budget consulted only afterwards, by which time every window had
        // been built. The dense cost is knowable in advance, so it is charged in advance.
        if (budget != nullptr) {
            const std::uint64_t cells = static_cast<std::uint64_t>(out.n_a) *
                                        static_cast<std::uint64_t>(out.n_b);
            // Upper bound on the whole-read comparisons the dense scan can perform: four mate
            // variants, each offered every start inside each window.
            std::uint64_t starts = 0;
            const std::size_t rlen = std::max(fragment.r1.size(), fragment.r2.size());
            for (std::size_t al = 0; al < out.n_a; ++al) {
                for (std::size_t be = 0; be < out.n_b; ++be) {
                    const std::size_t wlen = geom.lflank.size() + geom.alleles_a[al].size() +
                                             geom.context.size() + geom.alleles_b[be].size() +
                                             geom.rflank.size();
                    if (wlen >= rlen) starts += 4ull * (wlen - rlen + 1);
                }
            }
            if (!budget->cells_fit(cells) || !budget->verifications_fit(starts)) {
                budget->exhausted = true;
                budget->reason = "support-search-fallback-work-limit";
                return refuse("support-search-fallback-work-limit");
            }
            if (!budget->charge_cells(cells, "support-search-fallback-work-limit"))
                return refuse("support-search-fallback-work-limit");
        }
        // EXHAUSTIVE OVER THE VIRTUAL ALLELE PRODUCT -- every pair -- not over panel-carried pairs.
        for (std::size_t al = 0; al < out.n_a; ++al)
            for (std::size_t be = 0; be < out.n_b; ++be)
                out.mass[al * out.n_b + be] =
                    score_window(fragment, geom, ip, al, be, d1, d2, a1, a2, log_eps, log_1meps,
                                 &out.cell_states[al * out.n_b + be]);
    } else {
        if (budget != nullptr &&
            !budget->charge_cells(static_cast<std::uint64_t>(sup.proposals.size()),
                                  "support-search-proposed-cell-limit")) {
            return refuse("support-search-proposed-cell-limit");
        }
        // DIRECT POSITIONAL VERIFICATION. The positional search has already done the work the
        // window search would repeat: it produced, per mate variant, the starts at which that mate
        // COULD sit. So verify the whole mate at exactly those starts, through the shared virtual
        // window, and accumulate the fragment mass from the surviving placements. Nothing is
        // materialised, and score_window is not called at all on this path.
        struct Placed { long start; std::uint32_t edits; };
        std::vector<std::vector<std::pair<std::uint64_t, Placed>>> ver(4);
        const std::string* seqs[4] = {&fragment.r1, &a1, &fragment.r2, &a2};
        const std::size_t bands[4] = {d1, d1, d2, d2};
        for (int mi = 0; mi < 4; ++mi) {
            VirtualWindow vw;
            std::uint64_t cur = ~0ull;
            for (const auto& st : sup.mate_states[mi]) {
                const std::uint32_t al = static_cast<std::uint32_t>(st.first >> 32);
                const std::uint32_t be = static_cast<std::uint32_t>(st.first & 0xFFFFFFFFu);
                // Only pairs the seed-compatible join kept can carry a verified state: verified
                // starts are a SUBSET of seeded starts, so a pair with no seed-compatible join has
                // no verified one either. Skipping the rest is a sound prefilter, not a heuristic.
                if (!keep_pair(sup, al, be, out.n_b)) continue;
                if (st.first != cur) { vw.bind_pair(geom, al, be); cur = st.first; }
                // CHARGED BEFORE THE COMPARISON, so the limit bounds work done rather than work
                // already paid for.
                if (budget != nullptr &&
                    !budget->charge_verification(seqs[mi]->size(),
                                                 "support-search-verification-limit")) {
                    return refuse("support-search-verification-limit");
                }
                ++out.full_read_verifications;
                const std::size_t mm = vw.count_mismatches(*seqs[mi], st.second, bands[mi]);
                if (mm > bands[mi]) continue;
                ++out.accepted_mate_placements;
                ver[mi].push_back({st.first, Placed{st.second, static_cast<std::uint32_t>(mm)}});
            }
        }
        std::vector<std::vector<std::array<std::uint32_t, 3>>>* collect = nullptr;
        const double half = std::log(0.5);
        const auto read_ll = [&](std::uint32_t e, std::size_t len) {
            return static_cast<double>(e) * log_eps +
                   static_cast<double>(len - e) * log_1meps;
        };
        // THE SHARED valid-FR RULE again, now over VERIFIED placements: this is the point at which
        // a state stops being a proposal and becomes a fragment placement.
        const auto join = [&](int fi_v, int ri_v, bool fwd_is_m1) {
            const auto& F = ver[fi_v];
            const auto& R = ver[ri_v];
            const std::size_t rev_len = fwd_is_m1 ? fragment.r2.size() : fragment.r1.size();
            std::size_t fi = 0, ri = 0;
            while (fi < F.size() && ri < R.size()) {
                if (F[fi].first < R[ri].first) { ++fi; continue; }
                if (R[ri].first < F[fi].first) { ++ri; continue; }
                const std::uint64_t kk = F[fi].first;
                std::size_t fe = fi, re = ri;
                while (fe < F.size() && F[fe].first == kk) ++fe;
                while (re < R.size() && R[re].first == kk) ++re;
                const std::size_t cell = static_cast<std::size_t>(kk >> 32) * out.n_b +
                                         static_cast<std::size_t>(kk & 0xFFFFFFFFu);
                for (std::size_t a = fi; a < fe; ++a) {
                    for (std::size_t b = ri; b < re; ++b) {
                        const long rev_end = R[b].second.start + static_cast<long>(rev_len) - 1;
                        if (!valid_fr_coordinates(F[a].second.start, rev_end, ip.lo, ip.hi)) continue;
                        ++out.verified_fr_states;
                        ++out.cell_states[cell];
                        const std::uint32_t e1 = fwd_is_m1 ? F[a].second.edits : R[b].second.edits;
                        const std::uint32_t e2 = fwd_is_m1 ? R[b].second.edits : F[a].second.edits;
                        const long insert = rev_end - F[a].second.start + 1;
                        // THE SAME STATE that contributes the mass contributes the signature, so
                        // the two cannot come to describe different things.
                        if (collect != nullptr) {
                            (*collect)[cell].push_back(
                                {e1, e2, static_cast<std::uint32_t>(insert)});
                        }
                        out.mass[cell] = log_add(out.mass[cell],
                                                 half + read_ll(e1, fragment.r1.size()) +
                                                 read_ll(e2, fragment.r2.size()) +
                                                 ip.log_at(insert));
                    }
                }
                fi = fe; ri = re;
            }
        };
        // THE STRUCTURAL SIGNATURE, collected alongside the mass from the same states, so the two
        // cannot describe different things.
        std::vector<std::vector<std::array<std::uint32_t, 3>>> sig;
        if (out_cell_signatures != nullptr) sig.assign(out.n_a * out.n_b, {});
        collect = out_cell_signatures != nullptr ? &sig : nullptr;
        join(0, 3, true);    // r1 forward with r2 reverse-complemented
        join(2, 1, false);   // r2 forward with r1 reverse-complemented
        if (out_cell_signatures != nullptr) {
            out_cell_signatures->assign(out.n_a * out.n_b, std::string());
            for (std::size_t k = 0; k < sig.size(); ++k) {
                std::sort(sig[k].begin(), sig[k].end());
                std::string& b = (*out_cell_signatures)[k];
                b.resize(sig[k].size() * 12);
                for (std::size_t j = 0; j < sig[k].size(); ++j) {
                    std::memcpy(&b[j * 12], sig[k][j].data(), 12);
                }
            }
        }
    }
    for (const double m : out.mass)
        if (m != kNegInf) ++out.finite_emission_cells;
    if (out_support != nullptr) *out_support = sup;
    mark_informative(out);
    out.ok = true;
    return out;
}

LinkageEmission linkage_emission(const Fragment& fragment, const LinkageGeometry& geom,
                                 const InsertPrior& ip, double max_divergence,
                                 double log_eps, double log_1meps, double log_p_bg) {
    LinkageEmission out;
    out.log_p_bg = log_p_bg;
    if (!geom.ok || fragment.r1.empty() || fragment.r2.empty()) return out;
    out.n_a = geom.alleles_a.size();
    out.n_b = geom.alleles_b.size();
    out.mass.assign(out.n_a * out.n_b, kNegInf);
    out.cell_states.assign(out.n_a * out.n_b, 0);
    const std::size_t d1 = mate_band_edits(max_divergence, fragment.r1.size());
    const std::size_t d2 = mate_band_edits(max_divergence, fragment.r2.size());
    const std::string a1 = reverse_complement(fragment.r1);
    const std::string a2 = reverse_complement(fragment.r2);
    const double half = std::log(0.5);
    for (std::size_t al = 0; al < out.n_a; ++al) {
        for (std::size_t be = 0; be < out.n_b; ++be) {
            const std::string win = geom.lflank + geom.alleles_a[al] + geom.context +
                                    geom.alleles_b[be] + geom.rflank;
            // NOTE ON COST, measured rather than assumed. This materialises and searches EVERY
            // (alpha, beta) window: C4's ten edges need 1,154,642 of them and the construction did
            // not finish. Building a piece index per window was tried and does NOT fix it -- the
            // index costs O(|win|) to build, which at 1.15M windows is the same order as the
            // scanning it replaces. The real reduction has to come from proposing only the allele
            // pairs a fragment's seeds can reach, which is a separate bounded-complete search, so
            // this stays the straightforward form until that exists.
            const auto f1 = bounded_mate_placements(fragment.r1, win, d1, nullptr, nullptr);
            const auto v1 = bounded_mate_placements(a1, win, d1, nullptr, nullptr);
            const auto f2 = bounded_mate_placements(fragment.r2, win, d2, nullptr, nullptr);
            const auto v2 = bounded_mate_placements(a2, win, d2, nullptr, nullptr);
            const auto st = enumerate_fragment_states(0, f1, v1, f2, v2, fragment.r1.size(),
                                                      fragment.r2.size(), ip.lo, ip.hi);
            out.cell_states[al * out.n_b + be] = static_cast<std::uint32_t>(st.size());
            double m = kNegInf;
            for (const FragmentState& z : st) {
                const double e1 = static_cast<double>(z.m1_edits) * log_eps +
                                  static_cast<double>(fragment.r1.size() - z.m1_edits) * log_1meps;
                const double e2 = static_cast<double>(z.m2_edits) * log_eps +
                                  static_cast<double>(fragment.r2.size() - z.m2_edits) * log_1meps;
                m = log_add(m, half + e1 + e2 + ip.log_at(z.insert));
            }
            out.mass[al * out.n_b + be] = m;
        }
    }
    // INFORMATIVE means the mass varies with the COMBINATION, not merely with one endpoint. A
    // fragment reaching a distinguishing position in only one block spans the junction yet says
    // nothing about phase: measured at 16 of 116 on the phase fixture.
    for (std::size_t al = 0; al < out.n_a && !out.informative; ++al) {
        for (std::size_t be = 1; be < out.n_b; ++be) {
            const double d = out.mass[al * out.n_b + be] - out.mass[al * out.n_b];
            if (!(std::abs(d) < 1e-12) &&
                !(out.mass[al * out.n_b + be] == kNegInf && out.mass[al * out.n_b] == kNegInf)) {
                out.informative = true; break;
            }
        }
    }
    if (out.informative) {
        // ...and it must vary with A too, or it is a one-endpoint signal that the unary already owns.
        bool varies_a = false;
        for (std::size_t be = 0; be < out.n_b && !varies_a; ++be) {
            for (std::size_t al = 1; al < out.n_a; ++al) {
                const double d = out.mass[al * out.n_b + be] - out.mass[be];
                if (!(std::abs(d) < 1e-12) &&
                    !(out.mass[al * out.n_b + be] == kNegInf && out.mass[be] == kNegInf)) {
                    varies_a = true; break;
                }
            }
        }
        out.informative = varies_a;
    }
    out.ok = true;
    return out;
}

SparseEdgeLinkage make_sparse_kernel_edge(const SparseLinkageEdge& edge,
                                          const AlleleMapping& map_a, const AlleleMapping& map_b) {
    SparseEdgeLinkage out;   // inactive by default: the safe answer
    if (!edge.usable()) return out;
    if (map_a.status != MappingStatus::Ok || map_b.status != MappingStatus::Ok) return out;
    if (map_a.n_alleles != edge.n_a || map_b.n_alleles != edge.n_b) return out;
    if (map_a.allele.size() != map_b.allele.size()) return out;
    out.n_a = edge.n_a;
    out.n_b = edge.n_b;
    out.allele_a = map_a.allele;
    out.allele_b = map_b.allele;
    out.group_a.assign(edge.n_a, {});
    out.group_b.assign(edge.n_b, {});
    for (std::uint32_t h = 0; h < out.allele_a.size(); ++h) {
        if (out.allele_a[h] >= edge.n_a || out.allele_b[h] >= edge.n_b) return SparseEdgeLinkage{};
        out.group_a[out.allele_a[h]].push_back(h);
        out.group_b[out.allele_b[h]].push_back(h);
    }
    out.classes.reserve(edge.delta.size());
    for (const auto& kv : edge.delta) {
        const std::size_t amin = static_cast<std::size_t>((kv.first >> 48) & 0xFFFF);
        const std::size_t amax = static_cast<std::size_t>((kv.first >> 32) & 0xFFFF);
        const std::size_t bmin = static_cast<std::size_t>((kv.first >> 16) & 0xFFFF);
        const std::size_t bmax = static_cast<std::size_t>(kv.first & 0xFFFF);
        SparsePhaseClass c;
        c.amin = static_cast<std::uint32_t>(amin); c.amax = static_cast<std::uint32_t>(amax);
        c.bmin = static_cast<std::uint32_t>(bmin); c.bmax = static_cast<std::uint32_t>(bmax);
        // expm1, so a psi near one keeps its correction rather than losing it to cancellation.
        c.straight_m1 = std::expm1(edge.log_psi(amin, bmin, amax, bmax));
        c.crossed_m1 = std::expm1(edge.log_psi(amin, bmax, amax, bmin));
        out.classes.push_back(c);
    }
    out.active = true;
    return out;
}

HybridActivation plan_hybrid_activation_sparse(
    const std::vector<Fragment>& fragments,
    const std::vector<FragmentOwner>& owners,
    const std::vector<EdgeStatusEntry>& edge_status,
    const std::map<std::pair<std::uint32_t, std::uint32_t>, SparseLinkageEdge>& edges,
    const std::vector<AlleleMapping>& maps,
    std::size_t n_blocks,
    const std::vector<HigherFactorScope>& higher) {
    HybridActivation A;
    A.kernel_edges.assign(n_blocks, ChainEdgeLinkage{});
    A.sparse_kernel_edges.assign(n_blocks, SparseEdgeLinkage{});
    A.report = assess_hybrid_completeness(owners, edge_status, higher);
    A.ownership_complete = A.report.ownership_complete;
    if (!A.ownership_complete) {
        A.refusal = "hybrid model incomplete: " +
                    std::to_string(A.report.unconsumed_wide) + " wide, " +
                    std::to_string(A.report.unconsumed_refused_edge) + " on refused edges, " +
                    std::to_string(A.report.unconsumed_unusable) + " unusable";
        return A;
    }
    // BUILD EVERY EDGE FIRST. Nine successes and one failure is still a failure: the transaction is
    // over all edges, not each edge separately.
    std::vector<SparseEdgeLinkage> built(n_blocks);
    for (const auto& kv : edges) {
        const std::uint32_t a = kv.first.first, b = kv.first.second;
        // A SUPERSEDED EDGE IS NOT BUILT AT ALL. Building it and then declining to use it would
        // leave a second copy of the same evidence one wiring mistake away from being applied.
        bool gone = false;
        for (const HigherFactorScope& h : higher) if (h.supersedes(a, b)) { gone = true; break; }
        if (gone) continue;
        if (b >= n_blocks || a >= maps.size() || b >= maps.size()) {
            A.refusal = "edge " + std::to_string(a) + "-" + std::to_string(b) + " is out of range";
            return A;
        }
        const SparseEdgeLinkage k = make_sparse_kernel_edge(kv.second, maps[a], maps[b]);
        if (!k.active) {
            A.refusal = "edge " + std::to_string(a) + "-" + std::to_string(b) +
                        " passed completeness but could not be built for the kernel";
            return A;
        }
        built[b] = k;
    }
    for (std::size_t i = 0; i < owners.size() && i < fragments.size(); ++i) {
        // CONSUMED BY A HIGHER FACTOR: a Wide fragment inside a span, or the owner of an edge that
        // span supersedes. Excluded exactly like a pairwise consumer's, so that the excluded set
        // and the consumed set stay the same set.
        bool by_higher = false;
        for (const HigherFactorScope& h : higher) {
            if (owners[i].kind == OwnerKind::Wide && h.covers(owners[i].panel_domain_var_scope)) {
                by_higher = true; break;
            }
            if (owners[i].kind == OwnerKind::Linkage &&
                h.supersedes(owners[i].block_lo, owners[i].block_hi)) { by_higher = true; break; }
        }
        if (by_higher) {
            A.excluded_fragments.push_back(fragments[i].name);
            ++A.consumed_fragments;
            continue;
        }
        if (owners[i].kind != OwnerKind::Linkage) continue;
        const std::uint32_t b = owners[i].block_hi;
        // ACTIVE means built, whether or not it carries classes: a sparse-neutral edge is a
        // legitimate consumer that simply contributes no correction.
        if (b >= n_blocks || !built[b].active) {
            A.refusal = "fragment " + fragments[i].name + " is owned by an edge that is not active";
            return A;
        }
        A.excluded_fragments.push_back(fragments[i].name);
        ++A.consumed_fragments;
    }
    A.factors_buildable = true;
    A.sparse_kernel_edges = std::move(built);
    for (const SparseEdgeLinkage& k : A.sparse_kernel_edges) if (k.active) ++A.active_edges;
    A.hybrid_activated = true;
    A.call_status = HybridCallStatus::Complete;
    return A;
}

HybridActivation plan_hybrid_activation(const std::vector<Fragment>& fragments,
                                        const std::vector<FragmentOwner>& owners,
                                        const std::vector<EdgeStatusEntry>& edge_status,
                                        const std::map<std::pair<std::uint32_t, std::uint32_t>,
                                                       LinkageEdge>& edges,
                                        const std::vector<AlleleMapping>& maps,
                                        std::size_t n_blocks) {
    HybridActivation A;
    A.kernel_edges.assign(n_blocks, ChainEdgeLinkage{});   // inactive: the safe default
    A.report = assess_hybrid_completeness(owners, edge_status);
    A.ownership_complete = A.report.ownership_complete;
    if (!A.ownership_complete) {
        // TRANSACTIONAL REFUSAL. Nothing is subtracted and nothing is activated, so the marker
        // unaries are exactly what the legacy caller would see. Subtracting here and then failing
        // to activate would delete this evidence from BOTH models.
        A.refusal = "hybrid model incomplete: " +
                    std::to_string(A.report.unconsumed_wide) + " wide, " +
                    std::to_string(A.report.unconsumed_refused_edge) + " on refused edges, " +
                    std::to_string(A.report.unconsumed_unusable) + " unusable";
        return A;
    }
    // Build every kernel edge FIRST. If any required edge fails to build -- a mapping refused, an
    // allele count disagreeing -- the whole activation is abandoned rather than run with a hole.
    std::vector<ChainEdgeLinkage> built(n_blocks);
    for (const auto& kv : edges) {
        const std::uint32_t a = kv.first.first, b = kv.first.second;
        if (b >= n_blocks || a >= maps.size() || b >= maps.size()) {
            A.refusal = "edge " + std::to_string(a) + "-" + std::to_string(b) + " is out of range";
            return A;
        }
        const ChainEdgeLinkage k = make_kernel_edge(kv.second, maps[a], maps[b]);
        if (!k.active) {
            A.refusal = "edge " + std::to_string(a) + "-" + std::to_string(b) +
                        " passed completeness but could not be built for the kernel";
            return A;
        }
        built[b] = k;
    }
    // Only now, with every edge built, is the transaction allowed to commit. The exclusion list is
    // derived from fragments consumed by an ACTIVE edge -- not from "linkage-owned" -- so a refused
    // edge cannot satisfy the exclusion side.
    for (std::size_t i = 0; i < owners.size() && i < fragments.size(); ++i) {
        if (owners[i].kind != OwnerKind::Linkage) continue;
        const std::uint32_t b = owners[i].block_hi;
        if (b >= n_blocks || !built[b].active) {
            A.refusal = "fragment " + fragments[i].name + " is owned by an edge that is not active";
            return A;
        }
        A.excluded_fragments.push_back(fragments[i].name);
        ++A.consumed_fragments;
    }
    A.factors_buildable = true;
    A.kernel_edges = std::move(built);
    for (const ChainEdgeLinkage& k : A.kernel_edges) if (k.active) ++A.active_edges;
    A.hybrid_activated = true;
    A.call_status = HybridCallStatus::Complete;   // all four conditions hold
    return A;
}

ChainEdgeLinkage make_kernel_edge(const LinkageEdge& edge,
                                  const AlleleMapping& map_a, const AlleleMapping& map_b) {
    ChainEdgeLinkage out;   // inactive by default: the safe answer, not the convenient one
    if (!edge.usable()) return out;
    if (map_a.status != MappingStatus::Ok || map_b.status != MappingStatus::Ok) return out;
    if (map_a.n_alleles != edge.n_a || map_b.n_alleles != edge.n_b) return out;
    if (map_a.allele.size() != map_b.allele.size()) return out;
    if (edge.log_psi.size() != edge.n_a * edge.n_b * edge.n_a * edge.n_b) return out;
    out.active = true;
    out.n_a = edge.n_a;
    out.n_b = edge.n_b;
    out.allele_a = map_a.allele;
    out.allele_b = map_b.allele;
    out.log_psi = edge.log_psi;
    return out;
}

double estimate_fragment_lambda(std::size_t n_fragments,
                                const std::vector<std::size_t>& panel_lengths,
                                std::size_t* median_length_out) {
    if (median_length_out != nullptr) *median_length_out = 0;
    if (n_fragments == 0 || panel_lengths.empty()) return 0.0;
    std::vector<std::size_t> v = panel_lengths;
    std::sort(v.begin(), v.end());
    // MEDIAN, not mean: a panel with one truncated or one unusually long haplotype should not move
    // the depth scale, and at a length-variable locus the mean is exactly what does move.
    const std::size_t med = v.size() % 2 ? v[v.size() / 2]
                                         : (v[v.size() / 2 - 1] + v[v.size() / 2]) / 2;
    if (median_length_out != nullptr) *median_length_out = med;
    if (med == 0) return 0.0;
    return static_cast<double>(n_fragments) / (2.0 * static_cast<double>(med));
}

std::vector<HigherFactorScope> plan_higher_factors(
    const std::vector<FragmentOwner>& owners,
    const std::vector<EdgeStatusEntry>& edge_status,
    std::size_t n_blocks) {
    // ONE STATISTICAL FACTOR PER EXACT MINIMAL DEPENDENCY SCOPE.
    //
    // An earlier version dropped any scope contained in another, so fragments depending on
    // {3,4}, {4,5} and {3,4,5} were pooled into a single factor over {3,4,5} and normalised
    // ONCE. That silently changes the statistical model, and changes it PER DONOR: which scopes
    // happen to appear decides how the rest of the evidence is normalised. Mean-one
    // normalisation is defined over a factor's own content classes, so pooling distinct
    // dependency scopes before normalising is a different model, not an optimisation.
    //
    // So containment no longer merges anything. Each distinct scope keeps its own factor and its
    // own normalisation. Factors that overlap are evaluated together by the recurrence -- the
    // history already carries the common refinement of every factor touching a block -- and that
    // shared evaluation is an INFERENCE CLIQUE: it owns no evidence and performs no
    // normalisation, it only multiplies independently normalised potentials.
    std::set<std::vector<std::uint32_t>> scopes;
    for (const FragmentOwner& o : owners) {
        if (o.kind != OwnerKind::Wide || o.panel_domain_var_scope.empty()) continue;
        std::uint32_t lo = o.panel_domain_var_scope.front(), hi = o.panel_domain_var_scope.front();
        for (std::uint32_t b : o.panel_domain_var_scope) { lo = std::min(lo, b); hi = std::max(hi, b); }
        if (hi >= n_blocks) continue;
        std::vector<std::uint32_t> sp;              // contiguous closure, the only enlargement
        for (std::uint32_t b = lo; b <= hi; ++b) sp.push_back(b);
        scopes.insert(sp);
    }
    // A REFUSED pairwise edge leaves its owners with no consumer, and those owners depend on
    // exactly two blocks. That is its own minimal scope and therefore its own factor -- not
    // evidence to be folded into whichever wider factor happens to span it.
    std::set<std::vector<std::uint32_t>> from_edges;
    for (const EdgeStatusEntry& e : edge_status) {
        if (e.status == LinkageStatus::Ok) continue;
        if (e.block_b >= n_blocks || e.n_fragments == 0) continue;
        std::vector<std::uint32_t> sp;
        for (std::uint32_t b = e.block_a; b <= e.block_b; ++b) sp.push_back(b);
        scopes.insert(sp);
        from_edges.insert(sp);
    }
    std::vector<HigherFactorScope> out;
    for (const std::vector<std::uint32_t>& sp : scopes) {
        HigherFactorScope f;
        f.blocks = sp;
        // A factor built on a refused edge's scope supersedes THAT edge and nothing else.
        if (from_edges.count(sp) && sp.size() >= 2)
            f.superseded.push_back({sp.front(), sp.back()});
        out.push_back(std::move(f));
    }
    std::sort(out.begin(), out.end(),
              [](const HigherFactorScope& a, const HigherFactorScope& b) {
                  return a.blocks < b.blocks;
              });
    return out;
}

HybridCompletenessReport assess_hybrid_completeness(
    const std::vector<FragmentOwner>& owners,
    const std::vector<EdgeStatusEntry>& edges,
    const std::vector<HigherFactorScope>& higher) {
    HybridCompletenessReport R;
    R.owned_total = owners.size();
    std::map<std::pair<std::uint32_t, std::uint32_t>, LinkageStatus> st;
    for (const EdgeStatusEntry& e : edges) st[{e.block_a, e.block_b}] = e.status;
    for (const FragmentOwner& o : owners) {
        switch (o.kind) {
            case OwnerKind::Invariant:
                ++R.invariant;               // carries depth, not genotype evidence
                break;
            case OwnerKind::Unary:
                ++R.consumed_unary;          // the block's marker unary consumes it
                break;
            case OwnerKind::Linkage: {
                // A SUPERSEDED edge's owners belong to the factor that replaced it -- checked
                // before the edge's own status, because a superseded edge's refusal is no longer
                // anyone's problem.
                bool taken = false;
                for (const HigherFactorScope& h : higher)
                    if (h.supersedes(o.block_lo, o.block_hi)) { taken = true; break; }
                if (taken) { ++R.consumed_higher_superseded; break; }
                const auto it = st.find({o.block_lo, o.block_hi});
                if (it != st.end() && it->second == LinkageStatus::Ok) ++R.consumed_linkage;
                else ++R.unconsumed_refused_edge;
                break;
            }
            case OwnerKind::Wide: {
                // A higher factor consumes it when its WHOLE variable scope is inside the span.
                bool taken = false;
                for (const HigherFactorScope& h : higher)
                    if (h.covers(o.panel_domain_var_scope)) { taken = true; break; }
                if (taken) { ++R.consumed_higher_wide; break; }
                // Otherwise no consumer exists for three or more variables. Reported with its
                // scope rather than counted anonymously, and never cropped into a pair.
                ++R.unconsumed_wide;
                R.wide_scopes.push_back(o.panel_domain_var_scope);
                break;
            }
            default:
                ++R.unconsumed_unusable;     // explicitly reported missing evidence
                break;
        }
    }
    for (const EdgeStatusEntry& e : edges) {
        if (e.status == LinkageStatus::Ok) continue;
        bool superseded = false;
        for (const HigherFactorScope& h : higher)
            if (h.supersedes(e.block_a, e.block_b)) { superseded = true; break; }
        if (superseded) continue;   // replaced, so its refusal is not a defect in the model
        EdgeRefusal r;
        r.block_a = e.block_a; r.block_b = e.block_b;
        r.status = e.status; r.n_fragments = e.n_fragments; r.detail = e.detail;
        R.refusals.push_back(r);             // every one, in edge order
    }
    R.ownership_complete = (R.unconsumed_wide == 0 && R.unconsumed_refused_edge == 0 &&
                  R.unconsumed_unusable == 0);
    return R;
}

const char* hybrid_call_status_name(HybridCallStatus s) {
    return s == HybridCallStatus::Complete ? "COMPLETE" : "INCOMPLETE";
}

const char* mapping_status_name(MappingStatus s) {
    switch (s) {
        case MappingStatus::Ok:             return "ok";
        case MappingStatus::MissingMapping: return "missing-mapping";
        default:                            return "allele-out-of-range";
    }
}

AlleleMapping build_allele_mapping(const std::vector<int>& allele_of_block,
                                   std::size_t n_alleles, int bypass_allele) {
    AlleleMapping m;
    m.n_alleles = n_alleles;
    if (n_alleles == 0) {
        m.status = MappingStatus::OutOfRange;
        return m;
    }
    m.allele.resize(allele_of_block.size(), 0);
    for (std::size_t h = 0; h < allele_of_block.size(); ++h) {
        int a = allele_of_block[h];
        if (a < 0) {
            // A bypassing haplotype resolves to the block's bypass allele -- a real state, not
            // missing data. Without one there is nothing to resolve to, and inventing allele 0
            // would put the haplotype on somebody else's sequence.
            if (bypass_allele < 0) {
                m.status = MappingStatus::MissingMapping;
                m.first_bad_haplotype = h;
                m.first_bad_value = a;
                m.allele.clear();   // refused means NOTHING indexable, not a half-filled vector
                return m;
            }
            a = bypass_allele;
            ++m.n_bypass_resolved;
        }
        // VALIDATED BEFORE CONVERSION, both ends. a >= 0 is now established; the range check comes
        // before the cast, so no negative value can ever reach the unsigned index.
        if (a < 0 || static_cast<std::size_t>(a) >= n_alleles) {
            m.status = MappingStatus::OutOfRange;
            m.first_bad_haplotype = h;
            m.first_bad_value = a;
            m.allele.clear();   // same rule: a refused mapping exposes no partial answer
            return m;
        }
        m.allele[h] = static_cast<std::uint32_t>(a);
    }
    m.status = MappingStatus::Ok;
    return m;
}

const char* linkage_status_name(LinkageStatus s) {
    switch (s) {
        case LinkageStatus::Ok:                    return "ok";
        case LinkageStatus::ExposureDoesNotCancel: return "exposure-does-not-cancel";
        case LinkageStatus::InvalidEmissions:      return "invalid-emissions";
        case LinkageStatus::TooManyConfigurations: return "too-many-configurations";
        case LinkageStatus::ResourceExceeded:      return "resource-exceeded";
        case LinkageStatus::CountOverflow:         return "count-overflow";
        default:                                   return "not-computed";
    }
}

namespace {
inline std::uint64_t class_key(std::size_t amin, std::size_t amax,
                               std::size_t bmin, std::size_t bmax) {
    return (static_cast<std::uint64_t>(amin) << 48) | (static_cast<std::uint64_t>(amax) << 32) |
           (static_cast<std::uint64_t>(bmin) << 16) | static_cast<std::uint64_t>(bmax);
}
}  // namespace

double SparseLinkageEdge::log_psi(std::size_t a1, std::size_t b1,
                                  std::size_t a2, std::size_t b2) const {
    // A homozygous endpoint has no alternative phase: the class members are homologue swaps of one
    // another, S is equal across them, and mean-one centering gives exactly zero.
    if (a1 == a2 || b1 == b2) return 0.0;
    const std::size_t amin = std::min(a1, a2), amax = std::max(a1, a2);
    const std::size_t bmin = std::min(b1, b2), bmax = std::max(b1, b2);
    double d = 0.0;
    if (grouped) {
        // THE GROUPED LOOKUP. Alleles map to their signature classes and the class quadruple
        // carries the delta; the content class itself is never stored or visited.
        const auto itg = delta_class.find(class_key(row_class[amin], row_class[amax],
                                                    col_class[bmin], col_class[bmax]));
        if (itg == delta_class.end()) return 0.0;
        d = itg->second;
    } else {
        const auto it = delta.find(class_key(amin, amax, bmin, bmax));
        if (it == delta.end()) return 0.0;   // no corner carries in-band mass: exactly neutral
        d = it->second;
    }
    // STRAIGHT pairs amin with bmin. Homologue 1 carries (a1, b1), so the query is straight when
    // that pairing matches -- in either homologue order, which is the same biological phase.
    const bool straight = (a1 == amin && b1 == bmin) || (a1 == amax && b1 == bmax);
    const double x = straight ? -d : d;
    // log 2 - log(1 + e^x), guarded so a large |Delta| neither overflows nor loses the branch.
    const double kLog2 = std::log(2.0);
    if (x > 700.0) return kLog2 - x;
    if (x < -700.0) return kLog2;
    return kLog2 - std::log1p(std::exp(x));
}

SparseLinkageEdge build_sparse_linkage_edge_grouped(
    const std::vector<LinkageEmission>& emissions,
    const std::vector<std::vector<std::string>>& cell_signatures,
    const LinkageGeometry& geom, double lambda, double log_mix, double log_bg_weight,
    const SparseResourceLimits& limits, GroupedBuildStats* stats) {
    const auto t0 = std::chrono::steady_clock::now();
    SparseLinkageEdge E;
    GroupedBuildStats st;
    E.block_a = geom.block_a; E.block_b = geom.block_b;
    const std::size_t na = geom.alleles_a.size(), nb = geom.alleles_b.size();
    E.n_a = na; E.n_b = nb;
    E.theoretical_configs = na * nb * na * nb;
    if (na < 2 || nb < 2) { E.status = LinkageStatus::Ok; return E; }
    for (const LinkageEmission& m : emissions) {
        if (!m.ok) { E.status = LinkageStatus::InvalidEmissions; return E; }
    }
    if (cell_signatures.size() != emissions.size()) {
        E.status = LinkageStatus::InvalidEmissions; return E;
    }
    // Exposure must cancel, unchanged: the representation is grouped, the contract is not.
    // CANCELLATION IS ACCEPTED STRUCTURALLY, NOT NUMERICALLY.
    //
    // When every window is at least insert_hi - 1, exposure is AFFINE in the window length --
    // window_len = const + |A_alpha| + |B_beta| makes it additive, and the straight and crossed
    // sums coincide algebraically. `exposure_affine` is exactly that predicate, decided per
    // configuration when the geometry was built.
    //
    // The numerical difference is a DIAGNOSTIC and must not be the gate. A fixed absolute
    // tolerance rejected C4's edge 6-7 at 1.048e-9 with ZERO short windows -- rounding from
    // summing 650 prior terms instead of 401. Rescaling it would fix that case and introduce the
    // opposite failure: a genuine non-affine difference, small only because the exposure scale is
    // large, would pass. The algebraic precondition has neither failure mode, so it is the gate,
    // and the measured residual is carried alongside it for inspection.
    double asym = 0.0, scale = 1.0;
    for (const double e : geom.exposure) scale = std::max(scale, std::abs(e));
    for (std::size_t a1 = 0; a1 < na; ++a1)
    for (std::size_t b1 = 0; b1 < nb; ++b1)
    for (std::size_t a2 = 0; a2 < na; ++a2)
    for (std::size_t b2 = 0; b2 < nb; ++b2) {
        asym = std::max(asym, std::abs((geom.exposure[a1 * nb + b1] + geom.exposure[a2 * nb + b2]) -
                                       (geom.exposure[a1 * nb + b2] + geom.exposure[a2 * nb + b1])));
    }
    {
        char msg[128];
        std::snprintf(msg, sizeof msg,
                      "exposure residual %.6e absolute, %.6e relative to scale %.6e over "
                      "%zu x %zu cells", asym, asym / scale, scale, na, nb);
        E.refusal_detail = msg;   // a diagnostic on every edge, refused or not
    }
    if (!geom.exposure_affine) {
        E.status = LinkageStatus::ExposureDoesNotCancel;
        return E;
    }

    // ---- ESTIMATE BEFORE ALLOCATING ------------------------------------------------------------
    // Every major allocation is bounded first. An estimate computed after the fact is not a budget.
    std::size_t est_sig = 0;
    for (const auto& cs : cell_signatures) {
        if (cs.size() != na * nb) { E.status = LinkageStatus::InvalidEmissions; return E; }
        for (const std::string& x : cs) est_sig += x.size() + sizeof(std::string);
    }
    const std::size_t est_joint = na * nb * (sizeof(std::string) + 4 * cell_signatures.size());
    if (est_sig + est_joint > limits.max_bytes) {
        E.status = LinkageStatus::ResourceExceeded;
        E.bytes_grouped = est_sig + est_joint;
        return E;
    }
    st.estimates_checked = true;
    st.bytes_cell_signatures = est_sig;

    // The joint signature per cell: the concatenation over fragments, length-prefixed so two
    // different splits cannot alias.
    std::vector<std::string> joint(na * nb);
    for (const auto& cs : cell_signatures) {
        for (std::size_t k = 0; k < na * nb; ++k) {
            const std::uint32_t n = static_cast<std::uint32_t>(cs[k].size());
            joint[k].append(reinterpret_cast<const char*>(&n), 4);
            joint[k].append(cs[k]);
        }
    }
    const SignatureMatrix M = build_signature_matrix(joint, na, nb);
    if (!M.ok) { E.status = LinkageStatus::InvalidEmissions; return E; }
    const std::size_t RA = M.row_members.size(), RB = M.col_members.size();
    // THE KEY PACKS FOUR 16-BIT FIELDS. A locus with more classes than that would silently alias
    // two quadruples into one delta, so it is refused rather than approximated.
    if (RA > 65535 || RB > 65535) { E.status = LinkageStatus::CountOverflow; return E; }
    st.row_classes = RA; st.col_classes = RB;
    // Predicted grouped classes, checked BEFORE the delta map is built.
    std::size_t pred = 0;
    if (__builtin_mul_overflow(RA * RA, RB * RB, &pred)) {
        E.status = LinkageStatus::CountOverflow; return E;
    }
    E.predicted_classes = pred;
    if (pred > limits.max_classes) { E.status = LinkageStatus::ResourceExceeded; return E; }

    // Support cells, for the report; this is a per-cell scan, not a content-class one.
    for (std::size_t k = 0; k < na * nb; ++k) {
        for (const LinkageEmission& m : emissions)
            if (m.mass[k] != kNegInf) { ++E.support_cells; break; }
    }

    // ---- THE COLLAPSED TRAVERSAL ---------------------------------------------------------------
    // Ordered class pairs, kept only when some real member pair realises that order. Orientation
    // depends on it: for a class pair {i,j} which corner counts as straight turns on whether the
    // class-i member has the smaller allele index, and both orders occur because members are
    // scattered through the allele order.
    std::vector<std::pair<std::uint32_t, std::uint32_t>> apairs, bpairs;
    for (std::uint32_t i = 0; i < RA; ++i)
    for (std::uint32_t j = 0; j < RA; ++j) {
        const bool okp = (i == j) ? M.row_members[i].size() >= 2
                                  : M.row_min[i] < M.row_max[j];
        if (okp) apairs.emplace_back(i, j);
    }
    for (std::uint32_t u = 0; u < RB; ++u)
    for (std::uint32_t v = 0; v < RB; ++v) {
        const bool okp = (u == v) ? M.col_members[u].size() >= 2
                                  : M.col_min[u] < M.col_max[v];
        if (okp) bpairs.emplace_back(u, v);
    }
    const double log_lambda = std::log(lambda);
    const auto mix = [&](double p, double q, double log_p_bg) {
        double sig = kNegInf;
        if (p != kNegInf) sig = p;
        if (q != kNegInf) sig = (sig == kNegInf) ? q : log_add(sig, q);
        if (sig != kNegInf) sig += log_mix + log_lambda;
        const double bg = log_bg_weight + log_p_bg;
        return (sig == kNegInf) ? bg : log_add(sig, bg);
    };
    for (const auto& ap : apairs) {
        for (const auto& bp : bpairs) {
            ++st.representative_visits;
            // Representative alleles realising this ordered class quadruple. Their signature ids --
            // and therefore every corner mass -- are shared by every member of the same classes.
            const std::size_t a1 = M.row_members[ap.first].front();
            const std::size_t a2 = (ap.first == ap.second) ? M.row_members[ap.first][1]
                                                           : M.row_members[ap.second].front();
            const std::size_t b1 = M.col_members[bp.first].front();
            const std::size_t b2 = (bp.first == bp.second) ? M.col_members[bp.first][1]
                                                           : M.col_members[bp.second].front();
            double d = 0.0;
            for (const LinkageEmission& m : emissions) {
                d += mix(m.mass[a1 * nb + b1], m.mass[a2 * nb + b2], m.log_p_bg) -
                     mix(m.mass[a1 * nb + b2], m.mass[a2 * nb + b1], m.log_p_bg);
            }
            ++st.pattern_evaluations;
            if (d != 0.0) {
                E.delta_class.emplace(class_key(ap.first, ap.second, bp.first, bp.second), d);
            }
        }
    }
    E.grouped = true;
    E.row_class = M.row_class; E.col_class = M.col_class;
    E.n_row_classes = RA; E.n_col_classes = RB;
    E.stored_classes = E.delta_class.size();
    st.bytes_matrix = M.cell_to_signature.size() * sizeof(std::uint32_t);
    st.bytes_members = (M.row_class.size() + M.col_class.size()) * sizeof(std::uint32_t) +
                       (na + nb) * sizeof(std::uint32_t);
    st.bytes_delta_class = E.delta_class.size() * (sizeof(std::uint64_t) + sizeof(double) +
                                                   2 * sizeof(void*));
    E.bytes_grouped = st.bytes_matrix + st.bytes_members + st.bytes_delta_class +
                      (E.row_class.size() + E.col_class.size()) * sizeof(std::uint32_t);
    if (E.bytes_total() > limits.max_bytes) {
        E.status = LinkageStatus::ResourceExceeded; return E;
    }
    E.status = LinkageStatus::Ok;
    st.build_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    if (stats != nullptr) *stats = st;
    return E;
}

SparseLinkageEdge build_sparse_linkage_edge(const std::vector<LinkageEmission>& emissions,
                                            const LinkageGeometry& geom, double lambda,
                                            double log_mix, double log_bg_weight,
                                            const SparseResourceLimits& limits) {
    SparseLinkageEdge E;
    if (!geom.ok) return E;
    E.block_a = geom.block_a; E.block_b = geom.block_b;
    E.n_a = geom.alleles_a.size(); E.n_b = geom.alleles_b.size();
    const std::size_t na = E.n_a, nb = E.n_b;
    E.theoretical_configs = na * nb * na * nb;
    if (na < 2 || nb < 2) { E.status = LinkageStatus::Ok; return E; }   // no phase to carry
    for (const LinkageEmission& m : emissions) {
        if (!m.ok) { E.status = LinkageStatus::InvalidEmissions; return E; }
    }
    // Exposure must cancel, exactly as in the dense path: the contract is unchanged, only the
    // representation is.
    // CANCELLATION IS ACCEPTED STRUCTURALLY, NOT NUMERICALLY.
    //
    // When every window is at least insert_hi - 1, exposure is AFFINE in the window length --
    // window_len = const + |A_alpha| + |B_beta| makes it additive, and the straight and crossed
    // sums coincide algebraically. `exposure_affine` is exactly that predicate, decided per
    // configuration when the geometry was built.
    //
    // The numerical difference is a DIAGNOSTIC and must not be the gate. A fixed absolute
    // tolerance rejected C4's edge 6-7 at 1.048e-9 with ZERO short windows -- rounding from
    // summing 650 prior terms instead of 401. Rescaling it would fix that case and introduce the
    // opposite failure: a genuine non-affine difference, small only because the exposure scale is
    // large, would pass. The algebraic precondition has neither failure mode, so it is the gate,
    // and the measured residual is carried alongside it for inspection.
    double asym = 0.0, scale = 1.0;
    for (const double e : geom.exposure) scale = std::max(scale, std::abs(e));
    for (std::size_t a1 = 0; a1 < na; ++a1)
    for (std::size_t b1 = 0; b1 < nb; ++b1)
    for (std::size_t a2 = 0; a2 < na; ++a2)
    for (std::size_t b2 = 0; b2 < nb; ++b2) {
        asym = std::max(asym, std::abs((geom.exposure[a1 * nb + b1] + geom.exposure[a2 * nb + b2]) -
                                       (geom.exposure[a1 * nb + b2] + geom.exposure[a2 * nb + b1])));
    }
    {
        char msg[128];
        std::snprintf(msg, sizeof msg,
                      "exposure residual %.6e absolute, %.6e relative to scale %.6e over "
                      "%zu x %zu cells", asym, asym / scale, scale, na, nb);
        E.refusal_detail = msg;   // a diagnostic on every edge, refused or not
    }
    if (!geom.exposure_affine) {
        E.status = LinkageStatus::ExposureDoesNotCancel;
        return E;
    }

    // THE SUPPORT: haploid cells carrying in-band mass, and which fragments carry them. Everything
    // outside is all-background and cannot separate the phases.
    std::vector<std::vector<std::uint32_t>> at_cell(na * nb);
    for (std::uint32_t fi = 0; fi < emissions.size(); ++fi) {
        const LinkageEmission& m = emissions[fi];
        for (std::size_t k = 0; k < m.mass.size(); ++k) {
            if (m.mass[k] != kNegInf) at_cell[k].push_back(fi);
        }
    }
    for (const auto& v : at_cell) if (!v.empty()) ++E.support_cells;

    // PREDICT BEFORE BUILDING. Each supported cell (alpha, beta) induces at most
    // (n_a - 1) * (n_b - 1) classes, so the total is bounded before a single one is allocated. The
    // arithmetic is OVERFLOW-CHECKED: a wrapped estimate would silently authorise an unbounded
    // build, which is the opposite of a budget.
    {
        std::size_t per_cell = 0, bound = 0;
        if (__builtin_mul_overflow(na - 1, nb - 1, &per_cell) ||
            __builtin_mul_overflow(E.support_cells, per_cell, &bound)) {
            E.status = LinkageStatus::CountOverflow;
            return E;
        }
        E.predicted_classes = bound;
        if (bound > limits.max_classes) {
            E.status = LinkageStatus::ResourceExceeded;
            return E;
        }
    }

    // AFFECTED CLASSES are induced by that support: a class matters only if one of its four corners
    // carries mass somewhere. Enumerated from the support rather than over all allele quadruples.
    std::unordered_map<std::uint64_t, char> affected;
    for (std::size_t al = 0; al < na; ++al)
    for (std::size_t be = 0; be < nb; ++be) {
        if (at_cell[al * nb + be].empty()) continue;
        for (std::size_t a2 = 0; a2 < na; ++a2) {
            if (a2 == al) continue;
            for (std::size_t b2 = 0; b2 < nb; ++b2) {
                if (b2 == be) continue;
                affected[class_key(std::min(al, a2), std::max(al, a2),
                                   std::min(be, b2), std::max(be, b2))] = 1;
            }
        }
    }
    const double log_lambda = std::log(lambda);
    for (const auto& kv : affected) {
        const std::size_t amin = static_cast<std::size_t>((kv.first >> 48) & 0xFFFF);
        const std::size_t amax = static_cast<std::size_t>((kv.first >> 32) & 0xFFFF);
        const std::size_t bmin = static_cast<std::size_t>((kv.first >> 16) & 0xFFFF);
        const std::size_t bmax = static_cast<std::size_t>(kv.first & 0xFFFF);
        // Only fragments touching one of the four corners can differ between the phases; the rest
        // contribute the same all-background term to both and cancel exactly.
        std::vector<std::uint32_t> frags;
        for (const std::size_t c : {amin * nb + bmin, amin * nb + bmax,
                                    amax * nb + bmin, amax * nb + bmax}) {
            frags.insert(frags.end(), at_cell[c].begin(), at_cell[c].end());
        }
        std::sort(frags.begin(), frags.end());
        frags.erase(std::unique(frags.begin(), frags.end()), frags.end());
        double d = 0.0;
        for (std::uint32_t fi : frags) {
            const LinkageEmission& m = emissions[fi];
            const auto mix = [&](double p, double q) {
                double sig = kNegInf;
                if (p != kNegInf) sig = p;
                if (q != kNegInf) sig = (sig == kNegInf) ? q : log_add(sig, q);
                if (sig != kNegInf) sig += log_mix + log_lambda;
                const double bg = log_bg_weight + m.log_p_bg;
                return (sig == kNegInf) ? bg : log_add(sig, bg);
            };
            d += mix(m.mass[amin * nb + bmin], m.mass[amax * nb + bmax]) -
                 mix(m.mass[amin * nb + bmax], m.mass[amax * nb + bmin]);
        }
        if (d != 0.0) E.delta.emplace(kv.first, d);
    }
    E.stored_classes = E.delta.size();
    // MEASURED, including container overhead: bucket array, per-node next-pointer, key and value.
    // The Delta payload alone (8 bytes x classes) is not the allocation.
    E.bytes_delta = E.delta.bucket_count() * sizeof(void*) +
                    E.delta.size() * (sizeof(std::uint64_t) + sizeof(double) + sizeof(void*));
    E.bytes_support = at_cell.capacity() * sizeof(std::vector<std::uint32_t>);
    for (const auto& v : at_cell) E.bytes_support += v.capacity() * sizeof(std::uint32_t);
    // ...and the MEASURED footprint is checked too. The prediction is an upper bound on classes,
    // not on bytes, so a build that fits the class budget can still exceed the byte budget.
    if (E.bytes_total() > limits.max_bytes) {
        // NOTHING PARTIALLY BUILT REMAINS INDEXABLE: the tables are released, exactly as a refused
        // dense edge exposes none.
        E.delta.clear();
        E.stored_classes = 0;
        E.bytes_delta = 0;
        E.bytes_support = 0;
        E.status = LinkageStatus::ResourceExceeded;
        return E;
    }
    E.status = LinkageStatus::Ok;
    return E;
}

LinkageEdge aggregate_linkage_edge(const std::vector<LinkageEmission>& emissions,
                                   const LinkageGeometry& geom, double lambda,
                                   double log_mix, double log_bg_weight,
                                   std::size_t max_configs) {
    LinkageEdge E;
    if (!geom.ok) return E;
    E.block_a = geom.block_a; E.block_b = geom.block_b;
    E.n_a = geom.alleles_a.size(); E.n_b = geom.alleles_b.size();
    const std::size_t na = E.n_a, nb = E.n_b;
    // REFUSE, DO NOT TRUNCATE. At LPA scale (457 x 410) this table is 35.1 billion entries and
    // 561.7 GB; quietly restricting to the top few alleles would turn the marker shortlist into an
    // uncertified linkage cutoff, which is the defect this work exists to remove.
    const std::size_t cap = max_configs ? max_configs : 1000000u;
    const bool overflow = na != 0 && nb != 0 &&
                          (static_cast<double>(na) * nb * na * nb > static_cast<double>(cap));
    if (overflow) {
        E.status = LinkageStatus::TooManyConfigurations;
        return E;
    }
    const std::size_t ncfg = na * nb * na * nb;
    E.score.assign(ncfg, 0.0);
    E.log_psi.assign(ncfg, 0.0);
    E.n_fragments = emissions.size();
    for (const LinkageEmission& m : emissions) {
        if (!m.ok) ++E.n_invalid;
        else if (m.informative) ++E.n_informative;
    }
    const double log_lambda = std::log(lambda);

    // EXPOSURE MUST CANCEL WITHIN EVERY CONTENT CLASS, and is then not carried at all. Checking it
    // rather than assuming it is the whole point: the affine argument holds only above the insert
    // support, and a short window from a deletion or bypass breaks it. Where it fails the edge is
    // UNSUPPORTED -- keeping it would let an edge with no fragments move the model through an
    // exposure difference, which is exactly the guarantee being made here.
    for (std::size_t a1 = 0; a1 < na; ++a1)
    for (std::size_t b1 = 0; b1 < nb; ++b1)
    for (std::size_t a2 = 0; a2 < na; ++a2)
    for (std::size_t b2 = 0; b2 < nb; ++b2) {
        const double straight = geom.exposure[a1 * nb + b1] + geom.exposure[a2 * nb + b2];
        const double crossed  = geom.exposure[a1 * nb + b2] + geom.exposure[a2 * nb + b1];
        E.exposure_asymmetry = std::max(E.exposure_asymmetry, std::abs(straight - crossed));
    }
    E.exposure_cancels = E.exposure_asymmetry <= 1e-9;
    if (!E.exposure_cancels) {
        E.status = LinkageStatus::ExposureDoesNotCancel;
        return E;             // log_psi stays all-zero and must not be consumed
    }
    if (E.n_invalid > 0) {
        E.status = LinkageStatus::InvalidEmissions;
        return E;
    }
    // S_e: the two homologues combined ONCE, background mixed ONCE, summed over fragments, and
    // exposure charged ONCE for the edge -- outside the fragment loop, from the geometry.
    for (std::size_t a1 = 0; a1 < na; ++a1)
    for (std::size_t b1 = 0; b1 < nb; ++b1)
    for (std::size_t a2 = 0; a2 < na; ++a2)
    for (std::size_t b2 = 0; b2 < nb; ++b2) {
        const std::size_t c = ((a1 * nb + b1) * na + a2) * nb + b2;
        double acc = 0.0;
        for (const LinkageEmission& m : emissions) {
            // No !m.ok skip here: an invalid emission has already made the edge unusable above.
            // Skipping it silently would shrink the evidence set and report a confident answer from
            // fewer fragments than the ownership partition claims.
            const double m1 = m.mass[a1 * nb + b1];
            const double m2 = m.mass[a2 * nb + b2];
            double sig = kNegInf;
            if (m1 != kNegInf) sig = m1;
            if (m2 != kNegInf) sig = (sig == kNegInf) ? m2 : log_add(sig, m2);
            if (sig != kNegInf) sig += log_mix + log_lambda;
            const double bg = log_bg_weight + m.log_p_bg;
            acc += (sig == kNegInf) ? bg : log_add(sig, bg);
        }
        // NO EXPOSURE TERM. It is required to cancel (verified above) and is therefore not
        // carried, so an edge with zero fragments has S == 0 everywhere and log psi == 0 by
        // construction -- not by the accident of its alleles happening to be equal length.
        E.score[c] = acc;
    }
    // PHASE ONLY, as a MEAN-ONE LIKELIHOOD RATIO within each unordered-content class:
    //
    //     log psi(c) = S(c) - logmeanexp_{c' in C} S(c') = S(c) - logsumexp + log|C|
    //
    // NOT a sum-one conditional distribution. psi is MULTIPLIED INTO the existing Li-Stephens
    // transition, which already carries its own phase prior; a sum-one factor would count that
    // normalisation twice. Worse, it is not neutral when the edge says nothing: with S flat,
    // sum-one gives log psi = -log|C|, and the classes have DIFFERENT cardinalities -- 1 for
    // hom/hom, 2 for het/hom, 4 for het/het -- so a phase-uninformative edge would penalise
    // heterozygous content by up to log 4 = 1.3863 nats purely from class size. Mean-one centering
    // gives exactly 0 for every class size, so an uninformative edge is genuinely neutral.
    //
    // Done after aggregation, never per fragment -- per-fragment normalisation would let each
    // fragment choose its own configuration, the mosaic error one level down.
    //
    // WHAT THIS GUARANTEES, and no more: an edge with no phase information leaves the Li-Stephens
    // model untouched. It does NOT guarantee that an INFORMATIVE edge never changes content
    // ranking -- phase evidence can still move marginal content posteriors once it interacts with
    // nonuniform Li-Stephens weights, and claiming otherwise would overstate the centering.
    std::vector<char> done(ncfg, 0);
    for (std::size_t a1 = 0; a1 < na; ++a1)
    for (std::size_t b1 = 0; b1 < nb; ++b1)
    for (std::size_t a2 = 0; a2 < na; ++a2)
    for (std::size_t b2 = 0; b2 < nb; ++b2) {
        const std::size_t c = ((a1 * nb + b1) * na + a2) * nb + b2;
        if (done[c]) continue;
        // The class sharing this unordered endpoint content: the same {a1,a2} at A and {b1,b2} at B,
        // over both pairings and both homologue orders.
        std::vector<std::size_t> cls;
        const std::size_t A[2] = {a1, a2}, Bb[2] = {b1, b2};
        for (int p = 0; p < 2; ++p) for (int q = 0; q < 2; ++q) {
            const std::size_t x1 = A[p], y1 = Bb[q], x2 = A[1 - p], y2 = Bb[1 - q];
            cls.push_back(((x1 * nb + y1) * na + x2) * nb + y2);
        }
        std::sort(cls.begin(), cls.end());
        cls.erase(std::unique(cls.begin(), cls.end()), cls.end());
        double lse = kNegInf;
        for (std::size_t k : cls) lse = log_add(lse, E.score[k]);
        const double lmean = lse - std::log(static_cast<double>(cls.size()));
        for (std::size_t k : cls) { E.log_psi[k] = E.score[k] - lmean; done[k] = 1; }
    }
    E.status = LinkageStatus::Ok;
    return E;
}

OwnershipLedger ownership_ledger(const std::vector<FragmentOwner>& owners) {
    OwnershipLedger L;
    // Pooled in-band placement mass per class. This is a SIZE statistic: it says how much placement
    // mass a class holds, not how much any call depends on it. A fragment with a negligible share
    // can still carry a decisive likelihood RATIO between two candidates, and a ratio is what a
    // call turns on -- so this must never be quoted as the information a class contributes.
    double mass_all = kNegInf;
    double mass[5] = {kNegInf, kNegInf, kNegInf, kNegInf, kNegInf};
    const auto slot = [](OwnerKind k) -> int {
        switch (k) {
            case OwnerKind::Unary:     return 0;
            case OwnerKind::Linkage:   return 1;
            case OwnerKind::Wide:      return 2;
            case OwnerKind::Invariant: return 3;
            default:                   return 4;
        }
    };
    for (const FragmentOwner& o : owners) {
        switch (o.kind) {
            case OwnerKind::Unary:     ++L.unary; break;
            case OwnerKind::Linkage:   ++L.linkage; break;
            case OwnerKind::Wide:      ++L.wide; break;
            case OwnerKind::Invariant: ++L.invariant; break;
            default:                   ++L.unusable; break;
        }
        if (o.in_band == kNegInf) continue;
        mass_all = log_add(mass_all, o.in_band);
        const int i = slot(o.kind);
        mass[i] = log_add(mass[i], o.in_band);
    }
    L.total = owners.size();
    L.linkage_fragment_share = L.total ? static_cast<double>(L.linkage) /
                                         static_cast<double>(L.total) : 0.0;
    const auto share = [&](int i) {
        return (mass_all == kNegInf || mass[i] == kNegInf) ? 0.0 : std::exp(mass[i] - mass_all);
    };
    // EVERY class, so none of them can go unaccounted when the chain decides what to do with it.
    L.unary_in_band_mass_share     = share(0);
    L.linkage_in_band_mass_share   = share(1);
    L.wide_in_band_mass_share      = share(2);
    L.invariant_in_band_mass_share = share(3);
    L.unusable_in_band_mass_share  = share(4);
    return L;
}

void write_ownership_table(const std::string& path,
                           const std::vector<Fragment>& fragments,
                           const std::vector<FragmentOwner>& owners) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("genotype-frag: cannot write " + path);
    f.precision(10);
    // COLUMN NAMES CARRY THE DOMAIN. A reader who sees "certified" and does not know the
    // search only ever placed mates on panel haplotypes will draw a stronger conclusion
    // than the number supports; that is exactly the misreading that has to stop here.
    f << "fragment\towner\tblock_lo\tblock_hi\tpanel_domain_scope_size"
         "\tpanel_domain_var_scope_size\tspan_lo\tspan_hi"
         "\torigins\tin_band\tomitted_bound"
         "\tunmapped\tdropped_nats\tpanel_domain_certified\n";
    for (std::size_t i = 0; i < owners.size() && i < fragments.size(); ++i) {
        const FragmentOwner& o = owners[i];
        f << fragments[i].name << '\t' << owner_kind_name(o.kind) << '\t'
          << o.block_lo << '\t' << o.block_hi << '\t' << o.panel_domain_scope.size() << '\t'
          << o.panel_domain_var_scope.size() << '\t'
          << (o.panel_domain_scope.empty() ? 0u : o.panel_domain_scope.front()) << '\t'
          << (o.panel_domain_scope.empty() ? 0u : o.panel_domain_scope.back()) << '\t' << o.origins
          << '\t' << o.in_band << '\t' << o.omitted_bound << '\t' << o.unmapped << '\t'
          << o.dropped << '\t' << (o.panel_domain_certified ? 1 : 0) << '\n';
    }
    f.flush();
    if (!f) throw std::runtime_error("genotype-frag: write failed for " + path);
}


// ---- THE EXHAUSTIVE SCOPE ORACLE -------------------------------------------------------------
// Deliberately independent: its own log-sum-exp, its own placement walk, its own mixture. Sharing
// any of them with production would let one defect move both sides of every comparison together.
namespace {

double so_log_add(double a, double b) {
    if (a == -std::numeric_limits<double>::infinity()) return b;
    if (b == -std::numeric_limits<double>::infinity()) return a;
    const double hi = a > b ? a : b, lo = a > b ? b : a;
    return hi + std::log1p(std::exp(lo - hi));
}

char so_comp(char c) {
    switch (c) {
        case 'A': case 'a': return 'T';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        case 'T': case 't': return 'A';
        default: return 'N';
    }
}

std::string so_revcomp(const std::string& s) {
    std::string o(s.size(), 'N');
    for (std::size_t i = 0; i < s.size(); ++i) o[s.size() - 1 - i] = so_comp(s[i]);
    return o;
}

// Hamming at every start. No band, no seeding, no early exit: the point of an oracle is that it
// cannot miss a placement the production search would have skipped.
std::vector<std::size_t> so_edits_at_every_start(const std::string& read, const std::string& hap) {
    std::vector<std::size_t> out;
    if (read.empty() || hap.size() < read.size()) return out;
    out.resize(hap.size() - read.size() + 1, 0);
    for (std::size_t s = 0; s + read.size() <= hap.size(); ++s) {
        std::size_t e = 0;
        for (std::size_t i = 0; i < read.size(); ++i) if (hap[s + i] != read[i]) ++e;
        out[s] = e;
    }
    return out;
}

}  // namespace

double ScopeOracleParams::log_at(long L) const {
    if (L < insert_lo || L > insert_hi || insert_logp.empty())
        return -std::numeric_limits<double>::infinity();
    return insert_logp[static_cast<std::size_t>(L - insert_lo)];
}

std::vector<std::vector<std::uint32_t>> scope_oracle_domain(const ScopeOracleLocus& locus,
                                                            ScopeOracleDomain domain) {
    if (domain == ScopeOracleDomain::PanelOnly) return locus.panel_tuples;
    // THE FULL CARTESIAN PRODUCT ACROSS THE WHOLE LOCUS. Not a window around a previously inferred
    // scope: defining the domain from the old scope is the circularity that produced the bug.
    std::vector<std::vector<std::uint32_t>> out;
    const std::size_t n = locus.blocks();
    if (n == 0) return out;
    std::vector<std::uint32_t> t(n, 0);
    for (;;) {
        out.push_back(t);
        std::size_t j = n;
        while (j > 0) {
            --j;
            if (++t[j] < locus.block_alleles[j].size()) break;
            t[j] = 0;
            if (j == 0) return out;
        }
        if (j == 0 && t[0] == 0) break;
    }
    return out;
}

std::string scope_oracle_haplotype(const ScopeOracleLocus& locus,
                                   const std::vector<std::uint32_t>& tuple) {
    std::string h;
    for (std::size_t b = 0; b < locus.blocks() && b < tuple.size(); ++b) {
        const auto& al = locus.block_alleles[b];
        if (tuple[b] < al.size()) h += al[tuple[b]];
    }
    return h;
}

double scope_oracle_emission(const std::string& r1, const std::string& r2,
                             const std::string& hap, const ScopeOracleParams& prm,
                             std::size_t* n_origins, std::size_t* best_edits,
                             std::vector<std::size_t>* in_band_positions) {
    const double kNI = -std::numeric_limits<double>::infinity();
    if (n_origins) *n_origins = 0;
    if (best_edits) *best_edits = std::numeric_limits<std::size_t>::max();
    if (in_band_positions) in_band_positions->clear();
    if (r1.empty() || r2.empty() || hap.empty()) return kNI;

    // Four tables: each mate, forward and reverse-complemented, at every start.
    const std::vector<std::size_t> f1 = so_edits_at_every_start(r1, hap);
    const std::vector<std::size_t> v1 = so_edits_at_every_start(so_revcomp(r1), hap);
    const std::vector<std::size_t> f2 = so_edits_at_every_start(r2, hap);
    const std::vector<std::size_t> v2 = so_edits_at_every_start(so_revcomp(r2), hap);

    const auto read_ll = [&](std::size_t edits, std::size_t len) {
        return static_cast<double>(edits) * prm.log_eps +
               static_cast<double>(len - edits) * prm.log_1meps;
    };
    const double log_half_strand = std::log(0.5);
    double total = kNI;
    std::size_t norigins = 0, best = std::numeric_limits<std::size_t>::max();

    // BOTH LIBRARY ORIENTATIONS, enumerated the same way: one mate forward, the other reverse.
    // insert = rev_end - fwd_start + 1, the same coordinate rule production uses, so the fixture
    // stays comparable even though none of the code is shared.
    for (int which = 0; which < 2; ++which) {
        const std::vector<std::size_t>& fwd = which == 0 ? f1 : f2;
        const std::vector<std::size_t>& rev = which == 0 ? v2 : v1;
        const std::size_t flen = which == 0 ? r1.size() : r2.size();
        const std::size_t rlen = which == 0 ? r2.size() : r1.size();
        for (std::size_t fs = 0; fs < fwd.size(); ++fs) {
            for (long L = prm.insert_lo; L <= prm.insert_hi; ++L) {
                const double lp = prm.log_at(L);
                if (lp == kNI) continue;
                // rev_end = fwd_start + L - 1, so the reverse mate STARTS at rev_end - rlen + 1.
                const long rev_end = static_cast<long>(fs) + L - 1;
                const long rs = rev_end - static_cast<long>(rlen) + 1;
                if (rs < 0 || rs >= static_cast<long>(rev.size())) continue;
                if (rev_end < static_cast<long>(fs)) continue;          // FR: reverse downstream
                if (rev_end >= static_cast<long>(hap.size())) continue;
                const std::size_t e = fwd[fs] + rev[static_cast<std::size_t>(rs)];
                const double m = log_half_strand + read_ll(fwd[fs], flen) +
                                 read_ll(rev[static_cast<std::size_t>(rs)], rlen) + lp;
                total = so_log_add(total, m);
                ++norigins;
                if (e < best) best = e;
                if (e <= prm.band_edits && in_band_positions) {
                    in_band_positions->push_back(static_cast<std::size_t>(fs));
                    in_band_positions->push_back(static_cast<std::size_t>(rev_end));
                }
            }
        }
    }
    if (n_origins) *n_origins = norigins;
    if (best_edits) *best_edits = best;
    return total;
}

// How many DISTINCT physical origins attain the best edit count. One is an ordinary placement;
// more than one is a repeat, and collapsing them loses mass in proportion to the multiplicity.
std::size_t scope_oracle_best_multiplicity(const std::string& r1, const std::string& r2,
                                           const std::string& hap, const ScopeOracleParams& prm,
                                           std::size_t best_edits,
                                           std::vector<std::pair<std::size_t, std::size_t>>* out) {
    ScopeOracleParams tight = prm;
    tight.band_edits = best_edits;
    std::vector<std::size_t> pos;
    std::size_t no = 0, be = 0;
    scope_oracle_emission(r1, r2, hap, tight, &no, &be, &pos);
    std::set<std::pair<std::size_t, std::size_t>> ids;
    for (std::size_t k = 0; k + 1 < pos.size(); k += 2) ids.insert({pos[k], pos[k + 1]});
    if (out != nullptr) out->assign(ids.begin(), ids.end());
    return ids.size();
}

double scope_oracle_residual(const ScopeOracleFragmentResult& r,
                             const std::vector<std::uint32_t>& scope,
                             const ScopeOracleParams& prm) {
    // THE SUFFICIENCY TEST FOR A CLAIMED SCOPE. If the scope really captures everything the
    // fragment depends on, two tuples AGREEING on it must give the same diploid contribution.
    // Whatever they differ by is what the scope threw away -- measured through the same mixture
    // the demotion budget is declared in, so the two numbers are on one scale and comparable.
    const double kNI = -std::numeric_limits<double>::infinity();
    const double log_mix = std::log1p(-prm.eta), log_bgw = std::log(prm.eta);
    const double log_lambda = std::log(prm.lambda);
    const auto psi = [&](double ma, double mb) {
        double sig = kNI;
        if (ma != kNI) sig = ma;
        if (mb != kNI) sig = (sig == kNI) ? mb : so_log_add(sig, mb);
        if (sig != kNI) sig += log_mix + log_lambda;
        const double bg = log_bgw + prm.log_p_bg;
        return (sig == kNI) ? bg : so_log_add(sig, bg);
    };
    const auto agree = [&](const std::vector<std::uint32_t>& a,
                           const std::vector<std::uint32_t>& b) {
        for (std::uint32_t s : scope) if (s < a.size() && s < b.size() && a[s] != b[s]) return false;
        return true;
    };
    double worst = 0.0;
    // DIPLOID, because that is what the caller scores: hold one homologue fixed and vary the other
    // between two tuples the scope cannot tell apart.
    for (std::size_t i = 0; i < r.tuples.size(); ++i)
        for (std::size_t j = 0; j < r.tuples.size(); ++j) {
            if (!agree(r.tuples[i], r.tuples[j])) continue;
            for (std::size_t k = 0; k < r.tuples.size(); ++k)
                worst = std::max(worst, std::abs(psi(r.emission[i], r.emission[k]) -
                                                 psi(r.emission[j], r.emission[k])));
        }
    return worst;
}

ScopeOracleFragmentResult run_scope_oracle(const ScopeOracleLocus& locus,
                                           const std::string& name,
                                           const std::string& r1, const std::string& r2,
                                           const ScopeOracleParams& prm,
                                           ScopeOracleDomain domain,
                                           double demotion_tolerance,
                                           std::size_t max_tuples) {
    ScopeOracleFragmentResult R;
    R.name = name;
    R.tuples = scope_oracle_domain(locus, domain);
    const std::size_t nb = locus.blocks();
    R.structural.assign(nb, 0);
    R.block_delta.assign(nb, 0.0);
    R.block_witness.assign(nb, {0, 0});
    // THE RESOURCE BUDGET, CHECKED BEFORE ANY WORK. Exceeding it refuses outright: there is no
    // narrower domain to retreat to, because retreating to one is the defect.
    if (R.tuples.size() > max_tuples) {
        R.refusal = "domain of " + std::to_string(R.tuples.size()) +
                    " tuples exceeds the budget of " + std::to_string(max_tuples);
        return R;
    }

    // ---- EVERY TUPLE, MATERIALISED IN FULL ---------------------------------------------------
    std::vector<std::size_t> best_edits(R.tuples.size(), 0);
    std::set<std::uint32_t> span;
    for (std::size_t t = 0; t < R.tuples.size(); ++t) {
        const std::string hap = scope_oracle_haplotype(locus, R.tuples[t]);
        std::size_t no = 0, be = 0;
        std::vector<std::size_t> band_pos;
        R.emission.push_back(scope_oracle_emission(r1, r2, hap, prm, &no, &be, &band_pos));
        R.origins_enumerated += no;
        best_edits[t] = be;
        // IDENTITY, not statistics: (forward start, reverse end) separates two copies of a tandem
        // repeat that agree on edits and insert exactly.
        std::vector<std::pair<std::size_t, std::size_t>> ids;
        for (std::size_t k = 0; k + 1 < band_pos.size(); k += 2)
            ids.push_back({band_pos[k], band_pos[k + 1]});
        std::sort(ids.begin(), ids.end());
        R.in_band_origins.push_back(ids);
        R.best_origin_multiplicity.push_back(
            scope_oracle_best_multiplicity(r1, r2, hap, prm, be));
        if (be <= prm.band_edits) {
            R.any_in_band = true;
            // Which blocks an in-band origin physically touches -- the SPAN, which is not the
            // dependency set and is reported separately so the two can be compared.
            std::vector<std::size_t> bounds;
            std::size_t acc = 0;
            for (std::size_t b = 0; b < nb; ++b) {
                bounds.push_back(acc);
                acc += locus.block_alleles[b][R.tuples[t][b]].size();
            }
            bounds.push_back(acc);
            for (std::size_t k = 0; k + 1 < band_pos.size(); k += 2)
                for (std::size_t b = 0; b < nb; ++b)
                    if (band_pos[k] < bounds[b + 1] && band_pos[k + 1] >= bounds[b])
                        span.insert(static_cast<std::uint32_t>(b));
        }
    }
    R.apparent_span.assign(span.begin(), span.end());

    // ---- LAYER 1: STRUCTURAL DEPENDENCY, FROM HAPLOID EMISSION CHANGES ------------------------
    // Two tuples differing at EXACTLY ONE block. On the panel domain such a pair often does not
    // exist -- the panel ties blocks together -- and then the domain cannot even pose the question.
    // That is recorded as `unaskable` rather than silently read as "independent".
    for (std::size_t b = 0; b < nb; ++b) {
        bool askable = false, moves = false;
        for (std::size_t i = 0; i < R.tuples.size() && !moves; ++i)
            for (std::size_t j = i + 1; j < R.tuples.size(); ++j) {
                std::size_t diffs = 0;
                bool at_b = false;
                for (std::size_t k = 0; k < nb; ++k)
                    if (R.tuples[i][k] != R.tuples[j][k]) { ++diffs; if (k == b) at_b = true; }
                if (diffs != 1 || !at_b) continue;
                askable = true;
                if (R.emission[i] != R.emission[j]) { moves = true; break; }
            }
        R.structural[b] = moves ? 1 : 0;
        if (!askable) R.unaskable.push_back(static_cast<std::uint32_t>(b));
    }

    // ---- LAYER 2: CERTIFIED DEMOTION, THROUGH THE COMPLETE DIPLOID MIXTURE -------------------
    // psi = log( (1-eta) * lambda * (M_a + M_b) + eta * P_bg ), the contribution the caller
    // actually scores. NOT raw placement mass: a tail that is negligible against a 1e-6 floor is
    // decisive against e^-200, and this codebase has already paid for that conflation once.
    const double kNI = -std::numeric_limits<double>::infinity();
    const double log_mix = std::log1p(-prm.eta), log_bgw = std::log(prm.eta);
    const double log_lambda = std::log(prm.lambda);
    const auto psi = [&](double ma, double mb) {
        double sig = kNI;
        if (ma != kNI) sig = ma;
        if (mb != kNI) sig = (sig == kNI) ? mb : so_log_add(sig, mb);
        if (sig != kNI) sig += log_mix + log_lambda;
        const double bg = log_bgw + prm.log_p_bg;
        return (sig == kNI) ? bg : so_log_add(sig, bg);
    };
    // The worst swing block b can cause, over EVERY diploid pair in the domain, varying b on
    // either homologue while everything else is held fixed.
    for (std::size_t b = 0; b < nb; ++b) {
        double worst = 0.0;
        std::pair<std::size_t, std::size_t> wit{0, 0};
        for (std::size_t ia = 0; ia < R.tuples.size(); ++ia)
            for (std::size_t ib = 0; ib < R.tuples.size(); ++ib) {
                for (std::size_t ja = 0; ja < R.tuples.size(); ++ja) {
                    bool only_b = true;
                    for (std::size_t k = 0; k < nb; ++k)
                        if (k != b && R.tuples[ia][k] != R.tuples[ja][k]) { only_b = false; break; }
                    if (!only_b) continue;
                    const double d = std::abs(psi(R.emission[ia], R.emission[ib]) -
                                              psi(R.emission[ja], R.emission[ib]));
                    if (d > worst) { worst = d; wit = {ia, ja}; }
                }
            }
        R.block_delta[b] = worst;
        R.block_witness[b] = wit;      // the two tuples that force the dependency
    }

    // ---- THE SCOPE: STRUCTURAL AND NOT DEMOTABLE ---------------------------------------------
    // PER FRAGMENT, and therefore PROVISIONAL. Demotion is a locus-level claim; see
    // aggregate_scope_oracle, which is what a caller must actually act on.
    for (std::size_t b = 0; b < nb; ++b)
        if (R.structural[b] && R.block_delta[b] > demotion_tolerance)
            R.scope.push_back(static_cast<std::uint32_t>(b));
    return R;
}

ScopeOracleLocusResult aggregate_scope_oracle(
    const std::vector<ScopeOracleFragmentResult>& per_fragment, double locus_budget) {
    ScopeOracleLocusResult out;
    out.budget = locus_budget;
    std::size_t nb = 0;
    for (const ScopeOracleFragmentResult& r : per_fragment) nb = std::max(nb, r.structural.size());
    out.block_delta_sum.assign(nb, 0.0);

    // ---- THE STRUCTURAL SCOPE, RECORDED BEFORE ANY APPROXIMATION -----------------------------
    for (const ScopeOracleFragmentResult& r : per_fragment) {
        std::vector<std::uint32_t> sc;
        for (std::size_t b = 0; b < r.structural.size(); ++b)
            if (r.structural[b]) sc.push_back(static_cast<std::uint32_t>(b));
        out.structural_scope.push_back(std::move(sc));
    }
    for (const ScopeOracleFragmentResult& r : per_fragment)
        for (std::size_t b = 0; b < r.structural.size(); ++b)
            if (r.structural[b]) out.block_delta_sum[b] += r.block_delta[b];

    // ---- A REFUSED FRAGMENT INVALIDATES THE WHOLE CERTIFICATE --------------------------------
    // Not "contributes zero and the rest proceeds". Its scope is UNKNOWN, so it may depend on
    // blocks the allocation below just spent budget removing from everyone else. The entire
    // previous ownership ledger stands, unchanged, and the transaction does not activate.
    for (std::size_t i = 0; i < per_fragment.size(); ++i) {
        if (per_fragment[i].usable()) continue;
        out.refusal = "fragment " + per_fragment[i].name + " refused certification (" +
                      per_fragment[i].refusal + "); the aggregate certificate is void";
        out.effective_scope = out.structural_scope;   // nothing is demoted, nothing changes
        out.certified = false;
        return out;
    }

    // ---- ONE BUDGET, OVER (FRAGMENT, BLOCK) PAIRS -------------------------------------------
    // The atoms are per-fragment-per-block, because that is the granularity at which a dependency
    // is actually dropped. Aggregating to whole blocks first would force a fragment with a 1e-30
    // dependence on block 1 to carry it merely because some other fragment needs block 1 badly.
    //
    // ALLOCATION: cheapest first, deterministic, stopping at the budget. Ties break on
    // (fragment, block) so the result cannot depend on iteration order. This is exact rather than
    // optimal -- any subset whose bounds sum inside the budget carries the same guarantee -- and a
    // tighter allocator would be an optimization, not a change of contract.
    struct Cand { double bound; std::size_t frag; std::uint32_t block; };
    std::vector<Cand> cands;
    for (std::size_t i = 0; i < per_fragment.size(); ++i)
        for (std::size_t b = 0; b < per_fragment[i].structural.size(); ++b)
            if (per_fragment[i].structural[b])
                cands.push_back({per_fragment[i].block_delta[b], i,
                                 static_cast<std::uint32_t>(b)});
    std::sort(cands.begin(), cands.end(), [](const Cand& x, const Cand& y) {
        if (x.bound != y.bound) return x.bound < y.bound;
        if (x.frag != y.frag) return x.frag < y.frag;
        return x.block < y.block;
    });
    std::vector<std::set<std::uint32_t>> removed(per_fragment.size());
    for (const Cand& c : cands) {
        if (out.total_demoted_bound + c.bound > locus_budget) continue;
        out.total_demoted_bound += c.bound;
        out.removals.push_back({c.frag, c.block, c.bound});
        removed[c.frag].insert(c.block);
    }

    // ---- THE EFFECTIVE SCOPE: STRUCTURAL MINUS WHAT THE BUDGET PAID FOR ----------------------
    for (std::size_t i = 0; i < per_fragment.size(); ++i) {
        std::vector<std::uint32_t> sc;
        for (std::uint32_t b : out.structural_scope[i])
            if (!removed[i].count(b)) sc.push_back(b);
        out.effective_scope.push_back(std::move(sc));
    }
    out.certified = out.total_demoted_bound <= locus_budget;
    return out;
}



// ---- THE OPTIMIZED SCOPE CERTIFIER -----------------------------------------------------------
namespace {

struct CertPlacement {
    std::uint32_t block = 0; std::size_t off = 0;          // first base
    std::uint32_t end_block = 0; std::size_t end_off = 0;  // ONE PAST the last base
    std::vector<std::pair<std::uint32_t, std::uint32_t>> constraint;
    std::size_t edits = 0;
    bool operator<(const CertPlacement& o) const {
        return std::tie(block, off, end_block, end_off, constraint, edits) <
               std::tie(o.block, o.off, o.end_block, o.end_off, o.constraint, o.edits);
    }
};

// Build every template of `need` bases starting at (b, a, o), enumerating the alleles of the
// following blocks. THIS is the local region: it stops as soon as the read is covered, so an
// origin constrains a consecutive run and nothing else. Blocks past the run stay free.
void cert_expand(const ScopeOracleLocus& L, std::uint32_t b, std::uint32_t a, std::size_t o,
                 std::size_t need, std::string tmpl,
                 std::vector<std::pair<std::uint32_t, std::uint32_t>> cons,
                 std::vector<std::pair<std::string, CertPlacement>>& out,
                 std::size_t& work, const CertifierParams& cp, bool& overflow) {
    if (overflow || ++work > cp.max_contexts) { overflow = true; return; }
    if (b >= L.blocks() || a >= L.block_alleles[b].size()) return;
    if (cp.mut_no_multi_boundary && cons.size() >= 2) return;
    const std::string& seq = L.block_alleles[b][a];
    cons.push_back({b, a});
    const std::size_t take = std::min(seq.size() - std::min(o, seq.size()), need - tmpl.size());
    tmpl += seq.substr(std::min(o, seq.size()), take);
    if (tmpl.size() == need) {
        CertPlacement p;
        p.end_block = b; p.end_off = o + take;
        p.constraint = cons;
        out.push_back({tmpl, p});
        return;
    }
    // Not covered yet: continue into the next block, every allele of it. A short allele here is
    // what opens a multi-boundary case, and the work guard is what keeps it finite.
    for (std::uint32_t a2 = 0; a2 + 1 <= L.block_alleles[b + 1 < L.blocks() ? b + 1 : 0].size() &&
                               b + 1 < L.blocks(); ++a2)
        cert_expand(L, b + 1, a2, 0, need, tmpl, cons, out, work, cp, overflow);
}

// Every placement of `read` with at most `d` edits, over the WHOLE locus and the full allele
// product, found by scanning every start. Complete by construction -- no pigeonhole argument is
// needed because nothing is skipped -- and linear in the total allele length, not in the product.
std::vector<CertPlacement> cert_scan(const ScopeOracleLocus& L, const std::string& read,
                                     std::size_t d, const CertifierParams& cp, bool& overflow,
                                     std::size_t& contexts) {
    std::vector<CertPlacement> out;
    for (std::uint32_t b = 0; b < L.blocks(); ++b)
        for (std::uint32_t a = 0; a < L.block_alleles[b].size(); ++a)
            for (std::size_t o = 0; o < L.block_alleles[b][a].size(); ++o) {
                std::vector<std::pair<std::string, CertPlacement>> ctx;
                std::size_t work = 0;
                cert_expand(L, b, a, o, read.size(), std::string(), {}, ctx, work, cp, overflow);
                contexts += ctx.size();
                for (auto& cw : ctx) {
                    if (cw.first.size() != read.size()) continue;
                    std::size_t e = 0;
                    for (std::size_t i = 0; i < read.size() && e <= d; ++i)
                        if (cw.first[i] != read[i]) ++e;
                    if (e > d) continue;
                    cw.second.block = b; cw.second.off = o; cw.second.edits = e;
                    out.push_back(cw.second);
                }
            }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end(), [](const CertPlacement& x,
                                                     const CertPlacement& y) {
        return !(x < y) && !(y < x);
    }), out.end());
    return out;
}

// The same set, reached by PIGEONHOLE. A placement with at most d edits leaves at least one of
// its d+1 disjoint pieces exact, so a seed index over allele interiors and junctions must find
// it. Verified against cert_scan on the fixture: this is the path that makes C4 tractable, and
// an unverified fast path is how a search silently stops being complete.
std::vector<CertPlacement> cert_seeded(const ScopeOracleLocus& L, const LocusIndex& ix,
                                       const std::string& read, std::size_t d,
                                       const CertifierParams& cp, bool& overflow,
                                       std::size_t& contexts, std::size_t& seeds) {
    std::vector<CertPlacement> out;
    const std::size_t k = ix.k;
    if (k == 0 || read.size() < k * (d + 1)) { overflow = true; return out; }
    std::set<std::pair<std::uint32_t, std::pair<std::uint32_t, std::size_t>>> starts;
    for (std::size_t i = 0; i <= d; ++i) {
        const std::string piece = read.substr(i * k, k);
        const auto it = ix.at.find(piece);
        if (it == ix.at.end()) continue;
        for (const LocusIndex::Seed& sd : it->second) {
            ++seeds;
            // The read would start i*k bases before the seed. Inside this allele, or -- when that
            // runs off its front -- inside some allele of an earlier block, every one of which is
            // a candidate: the walk back must not assume which.
            const long want = static_cast<long>(sd.offset) - static_cast<long>(i * k);
            if (want >= 0) {
                starts.insert({sd.block, {sd.allele, static_cast<std::size_t>(want)}});
                continue;
            }
            long deficit = -want;
            for (long pb = static_cast<long>(sd.block) - 1; pb >= 0 && deficit > 0; --pb) {
                for (std::uint32_t pa = 0; pa < L.block_alleles[pb].size(); ++pa) {
                    const long plen = static_cast<long>(L.block_alleles[pb][pa].size());
                    if (plen >= deficit)
                        starts.insert({static_cast<std::uint32_t>(pb),
                                       {pa, static_cast<std::size_t>(plen - deficit)}});
                }
                long shortest = std::numeric_limits<long>::max();
                for (const std::string& s2 : L.block_alleles[pb])
                    shortest = std::min(shortest, static_cast<long>(s2.size()));
                deficit -= shortest;
            }
        }
    }
    for (const auto& st : starts) {
        std::vector<std::pair<std::string, CertPlacement>> ctx;
        std::size_t work = 0;
        cert_expand(L, st.first, st.second.first, st.second.second, read.size(), std::string(),
                    {}, ctx, work, cp, overflow);
        contexts += ctx.size();
        for (auto& cw : ctx) {
            if (cw.first.size() != read.size()) continue;
            std::size_t e = 0;
            for (std::size_t i = 0; i < read.size() && e <= d; ++i)
                if (cw.first[i] != read[i]) ++e;
            if (e > d) continue;
            cw.second.block = st.first; cw.second.off = st.second.second; cw.second.edits = e;
            out.push_back(cw.second);
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end(), [](const CertPlacement& x,
                                                     const CertPlacement& y) {
        return !(x < y) && !(y < x);
    }), out.end());
    return out;
}

}  // namespace

LocusIndex build_locus_index(const ScopeOracleLocus& locus, std::size_t k,
                             const CertifierParams& cp) {
    LocusIndex ix;
    ix.k = k;
    if (k == 0) { ix.refusal = "seed length zero"; return ix; }
    ix.allele_len.resize(locus.blocks());
    ix.len_classes.resize(locus.blocks());
    for (std::size_t b = 0; b < locus.blocks(); ++b) {
        std::set<std::size_t> ls;
        for (const std::string& s : locus.block_alleles[b]) {
            ix.allele_len[b].push_back(s.size());
            ls.insert(s.size());
        }
        ix.len_classes[b].assign(ls.begin(), ls.end());
    }
    // PANEL-INDEPENDENT BY CONSTRUCTION: a function of the alleles alone. Interiors, then every
    // adjacent junction, then the short runs a k-mer can cross. Nothing here has ever seen a path.
    std::size_t n = 0;
    for (std::uint32_t b = 0; b < locus.blocks(); ++b)
        for (std::uint32_t a = 0; a < locus.block_alleles[b].size(); ++a) {
            const std::string& s = locus.block_alleles[b][a];
            for (std::size_t o = 0; o + k <= s.size(); ++o) {
                ix.at[s.substr(o, k)].push_back({b, a, static_cast<std::uint32_t>(o)});
                if (++n > cp.max_seed_positions) {
                    ix.refusal = "seed index exceeds " + std::to_string(cp.max_seed_positions) +
                                 " positions";
                    return ix;
                }
            }
            // Junction and multi-boundary k-mers: start in this allele's tail and continue into
            // the following blocks, every allele, until k bases are covered.
            const std::size_t first = s.size() >= k ? s.size() - k + 1 : 0;
            for (std::size_t o = first; o < s.size() && !cp.mut_no_junction_seeds; ++o) {
                bool ovf = false;
                std::vector<std::pair<std::string, CertPlacement>> ctx;
                std::size_t work = 0;
                CertifierParams cp2 = cp;
                cp2.max_contexts = cp.max_short_run;
                cert_expand(locus, b, a, o, k, std::string(), {}, ctx, work, cp2, ovf);
                for (const auto& cw : ctx)
                    if (cw.first.size() == k) {
                        ix.at[cw.first].push_back({b, a, static_cast<std::uint32_t>(o)});
                        if (++n > cp.max_seed_positions) {
                            ix.refusal = "seed index exceeds " +
                                         std::to_string(cp.max_seed_positions) + " positions";
                            return ix;
                        }
                    }
            }
        }
    return ix;
}


CertifierResult certify_fragment_scope(const ScopeOracleLocus& locus, const LocusIndex& index,
                                       const std::string& r1, const std::string& r2,
                                       const ScopeOracleParams& prm, const CertifierParams& cp) {
    CertifierResult R;
    const double kNI = -std::numeric_limits<double>::infinity();
    const std::size_t nb = locus.blocks();
    R.block_delta.assign(nb, 0.0);
    if (!index.usable()) { R.refusal = "index: " + index.refusal; return R; }

    // ---- 1 + 2. SYMBOLIC ORIGIN DISCOVERY, LOCUS-WIDE ---------------------------------------
    // Both mates, both orientations, over the whole locus. Never seeded from a placement or a
    // prior scope: a fragment with origins in two distant repeats has BOTH found here, and a
    // search that started from one of them would call the other region independent.
    bool ovf = false;
    std::size_t ctxs = 0, seeds = 0;
    const std::string a1 = reverse_complement(r1), a2 = reverse_complement(r2);
    const auto find_all = [&](const std::string& rd) {
        return cp.use_seeds ? cert_seeded(locus, index, rd, cp.band_edits, cp, ovf, ctxs, seeds)
                            : cert_scan(locus, rd, cp.band_edits, cp, ovf, ctxs);
    };
    const std::vector<CertPlacement> f1 = find_all(r1);
    const std::vector<CertPlacement> v1 = find_all(a1);
    const std::vector<CertPlacement> f2 = find_all(r2);
    const std::vector<CertPlacement> v2 = find_all(a2);
    R.contexts_examined = ctxs; R.seeds_examined = seeds;
    // A REFUSAL MUST NAME EVERY BLOCK IT COULD NOT DECIDE. Returning an empty uncertified list
    // beside a refusal invites exactly the reading this whole exercise exists to stop: an empty
    // set of problems. Nothing was decided here, so nothing is certified.
    const auto refuse_all = [&](const std::string& why) {
        R.refusal = why;
        R.uncertified_blocks.clear();
        for (std::uint32_t b = 0; b < nb; ++b) R.uncertified_blocks.push_back(b);
        R.structural_scope.clear();
        R.witnesses.clear();
    };
    if (ovf) { refuse_all("context expansion exceeded the work budget"); return R; }

    // ---- THE JOIN, BY THE EXACT INSERT RULE -------------------------------------------------
    const auto read_ll = [&](std::size_t e, std::size_t len) {
        return static_cast<double>(e) * prm.log_eps +
               static_cast<double>(len - e) * prm.log_1meps;
    };
    const double log_half = std::log(0.5);
    const auto join = [&](const std::vector<CertPlacement>& fwd,
                          const std::vector<CertPlacement>& rev,
                          std::size_t flen, std::size_t rlen, bool m1_forward) {
        for (const CertPlacement& pf : fwd)
            for (const CertPlacement& pr : rev) {
                if (pr.end_block < pf.block) continue;          // reverse mate upstream: not FR
                // Consistency: where both pin a block they must agree, or this pair of placements
                // describes no haplotype at all.
                std::map<std::uint32_t, std::uint32_t> merged;
                bool clash = false;
                for (const auto& c : pf.constraint) merged[c.first] = c.second;
                for (const auto& c : pr.constraint) {
                    const auto it = merged.find(c.first);
                    if (it != merged.end() && it->second != c.second) { clash = true; break; }
                    merged[c.first] = c.second;
                }
                if (clash) continue;
                // insert = sum of block lengths from pf.block to pr.end_block-1, less the start
                // offset, plus the end offset. Blocks in that span that NEITHER mate reads still
                // move the insert whenever their alleles differ in length -- a real dependency
                // through pi(L) alone, recorded as a length constraint rather than ignored.
                long known = 0;
                std::vector<std::uint32_t> free_blocks;
                for (std::uint32_t c = pf.block; c < pr.end_block; ++c) {
                    const auto it = merged.find(c);
                    if (it != merged.end()) {
                        known += static_cast<long>(locus.block_alleles[c][it->second].size());
                    } else if (index.len_classes[c].size() == 1) {
                        known += static_cast<long>(index.len_classes[c][0]);
                    } else {
                        free_blocks.push_back(c);
                    }
                }
                std::size_t combos = 1;
                for (std::uint32_t c : free_blocks) combos *= index.len_classes[c].size();
                if (combos > cp.max_short_run * 16) {
                    R.uncertified_blocks.insert(R.uncertified_blocks.end(),
                                                free_blocks.begin(), free_blocks.end());
                    continue;
                }
                for (std::size_t ci = 0; ci < combos; ++ci) {
                    long extra = 0;
                    std::size_t rest = ci;
                    std::vector<std::pair<std::uint32_t, std::size_t>> lc;
                    for (std::uint32_t c : free_blocks) {
                        const std::size_t n = index.len_classes[c].size();
                        const std::size_t pick = rest % n; rest /= n;
                        extra += static_cast<long>(index.len_classes[c][pick]);
                        lc.push_back({c, index.len_classes[c][pick]});
                    }
                    const long insert = known + extra - static_cast<long>(pf.off) +
                                        static_cast<long>(pr.end_off);
                    if (insert < prm.insert_lo || insert > prm.insert_hi) continue;
                    SymbolicOrigin so;
                    for (const auto& kv : merged) so.allele_constraint.push_back(kv);
                    so.length_constraint = lc;
                    so.fwd_block = pf.block; so.fwd_off = pf.off;
                    so.rev_block = pr.end_block; so.rev_off = pr.end_off;
                    so.insert = insert; so.m1_forward = m1_forward;
                    so.edits = pf.edits + pr.edits;
                    // Mass is a function of the CONSTRAINT alone, which is what lets every tuple
                    // compatible with it share one value and the free blocks stay unenumerated.
                    (void)read_ll; (void)flen; (void)rlen; (void)log_half;
                    R.origins.push_back(std::move(so));
                }
            }
    };
    join(f1, v2, r1.size(), r2.size(), true);
    if (!cp.mut_one_orientation) join(f2, v1, r2.size(), r1.size(), false);
    // Deduplicate on PHYSICAL IDENTITY only. Two origins with equal edits and equal insert at
    // different places are two origins, and collapsing them loses exactly their multiplicity.
    std::vector<SymbolicOrigin> uniq;
    for (const SymbolicOrigin& o : R.origins) {
        bool dup = false;
        for (const SymbolicOrigin& u : uniq) if (u.same_place_as(o)) { dup = true; break; }
        if (!dup) uniq.push_back(o);
    }
    R.origins = std::move(uniq);
    if (cp.mut_drop_one_origin && !R.origins.empty()) {
        // The BEST origin, not an arbitrary one: with a wide band the tail holds thousands of
        // negligible placements and removing one of those is not a defect anyone could observe.
        // Losing one copy of a tandem repeat is, and this is what that looks like.
        std::size_t best = 0;
        for (std::size_t q = 1; q < R.origins.size(); ++q)
            if (R.origins[q].edits < R.origins[best].edits) best = q;
        R.origins.erase(R.origins.begin() + static_cast<long>(best));
    }
    if (cp.mut_nearest_region_only && !R.origins.empty()) {
        std::uint32_t lo = R.origins.front().fwd_block;
        for (const SymbolicOrigin& o : R.origins) lo = std::min(lo, o.fwd_block);
        std::vector<SymbolicOrigin> near;
        for (const SymbolicOrigin& o : R.origins) if (o.fwd_block == lo) near.push_back(o);
        R.origins = std::move(near);
    }

    // ---- 4. BOUND THE COMPLEMENT ------------------------------------------------------------
    // Everything past the verified band still has finite Hamming mass. Bound the NUMBER of
    // omitted (start, insert, orientation) states on one haplotype and their maximum emission,
    // and carry the product into the same global ledger every demotion is paid from. Calling the
    // unsearched region independent is the one thing that is not allowed.
    {
        std::size_t locus_len = 0;
        for (std::size_t b = 0; b < nb; ++b) {
            std::size_t longest = 0;
            for (const std::string& s : locus.block_alleles[b]) longest = std::max(longest, s.size());
            locus_len += longest;
        }
        const double n_states = static_cast<double>(locus_len) *
                                static_cast<double>(prm.insert_hi - prm.insert_lo + 1) * 2.0;
        const std::size_t total_len = r1.size() + r2.size();
        const std::size_t e = cp.band_edits + 1;                 // the cheapest omitted state
        const double max_emission = log_half + read_ll(e, total_len) +
                                    (prm.insert_logp.empty()
                                         ? 0.0
                                         : *std::max_element(prm.insert_logp.begin(),
                                                             prm.insert_logp.end()));
        R.omitted_bound = cp.mut_no_omitted_tail ? kNI : std::log(n_states) + max_emission;
    }

    // ---- 5. DEPENDENCY CLOSURE, TO A FIXED POINT --------------------------------------------
    // The ACTIVE SET is every block any origin constrains -- by allele or by length. A block in
    // no origin cannot change the in-band emission at all, so it is provably absent from the
    // scope and the product is taken over the active set alone. That is the whole optimization:
    // the domain enumerated is local to the origins, not to the locus.
    std::set<std::uint32_t> active;
    for (const SymbolicOrigin& o : R.origins) {
        for (const auto& c : o.allele_constraint) active.insert(c.first);
        for (const auto& c : o.length_constraint) active.insert(c.first);
    }
    if (cp.mut_truncate_closure) {
        // ONE REGION, NO EXPANSION: the active set taken from the best origin alone. This is what
        // a scope-guided search does, and on a fragment with origins in two distant repeats it
        // silently drops the second region.
        active.clear();
        std::size_t best = 0;
        for (std::size_t q = 1; q < R.origins.size(); ++q)
            if (R.origins[q].edits < R.origins[best].edits) best = q;
        if (!R.origins.empty()) {
            for (const auto& c : R.origins[best].allele_constraint) active.insert(c.first);
            for (const auto& c : R.origins[best].length_constraint) active.insert(c.first);
        }
    }
    for (int round = 0; !cp.mut_truncate_closure; ++round) {
        const std::size_t before = active.size();
        (void)round;
        for (const SymbolicOrigin& o : R.origins) {
            bool touches = false;
            for (const auto& c : o.allele_constraint) if (active.count(c.first)) touches = true;
            for (const auto& c : o.length_constraint) if (active.count(c.first)) touches = true;
            if (!touches) continue;
            for (const auto& c : o.allele_constraint) active.insert(c.first);
            for (const auto& c : o.length_constraint) active.insert(c.first);
        }
        if (active.size() == before) break;
    }
    std::vector<std::uint32_t> act(active.begin(), active.end());
    double combos = 1.0;
    for (std::uint32_t b : act) combos *= static_cast<double>(locus.block_alleles[b].size());
    if (combos > static_cast<double>(cp.max_contexts)) {
        // A RESOURCE LIMIT, NOT A FINDING. These blocks are reported uncertified; they are never
        // reported independent, and the caller must leave their ownership alone.
        refuse_all("active set of " + std::to_string(act.size()) +
                   " blocks exceeds the enumeration budget");
        return R;
    }
    // Emission for one assignment of the active set: the sum over origins compatible with it.
    const auto emission_of = [&](const std::vector<std::uint32_t>& pick) {
        double tot = kNI;
        for (const SymbolicOrigin& o : R.origins) {
            bool ok = true;
            for (const auto& c : o.allele_constraint) {
                const auto it = std::find(act.begin(), act.end(), c.first);
                if (it == act.end()) continue;
                if (pick[static_cast<std::size_t>(it - act.begin())] != c.second) {
                    ok = false; break;
                }
            }
            if (!ok) continue;
            for (const auto& c : o.length_constraint) {
                const auto it = std::find(act.begin(), act.end(), c.first);
                if (it == act.end()) continue;
                const std::uint32_t a = pick[static_cast<std::size_t>(it - act.begin())];
                if (locus.block_alleles[c.first][a].size() != c.second) { ok = false; break; }
            }
            if (!ok) continue;
            const std::size_t flen = o.m1_forward ? r1.size() : r2.size();
            const std::size_t rlen = o.m1_forward ? r2.size() : r1.size();
            const double m = log_half + read_ll(o.edits, flen + rlen) + prm.log_at(o.insert);
            tot = (tot == kNI) ? m : (m == kNI ? tot
                                              : (tot > m ? tot + std::log1p(std::exp(m - tot))
                                                         : m + std::log1p(std::exp(tot - m))));
        }
        return tot;
    };
    std::vector<std::vector<std::uint32_t>> picks;
    {
        std::vector<std::uint32_t> p(act.size(), 0);
        for (;;) {
            picks.push_back(p);
            std::size_t j = act.size();
            bool done = act.empty();
            while (j > 0) {
                --j;
                if (++p[j] < locus.block_alleles[act[j]].size()) break;
                p[j] = 0;
                if (j == 0) done = true;
            }
            if (done) break;
        }
    }
    std::vector<double> em;
    em.reserve(picks.size());
    for (const auto& p : picks) em.push_back(emission_of(p));
    R.active_blocks = act;
    R.active_picks = picks;
    R.active_emission = em;

    const double log_mix = std::log1p(-prm.eta), log_bgw = std::log(prm.eta);
    const double log_lambda = std::log(prm.lambda);
    const auto psi = [&](double ma, double mb) {
        double sig = kNI;
        if (ma != kNI) sig = ma;
        if (mb != kNI) sig = (sig == kNI) ? mb
                                          : (sig > mb ? sig + std::log1p(std::exp(mb - sig))
                                                      : mb + std::log1p(std::exp(sig - mb)));
        if (sig != kNI) sig += log_mix + log_lambda;
        const double bg = log_bgw + prm.log_p_bg;
        if (sig == kNI) return bg;
        return sig > bg ? sig + std::log1p(std::exp(bg - sig))
                        : bg + std::log1p(std::exp(sig - bg));
    };
    for (std::size_t jb = 0; jb < act.size(); ++jb) {
        const std::uint32_t b = act[jb];
        bool moves = false;
        double worst = 0.0;
        DependencyWitness wit;
        wit.block = b;
        for (std::size_t i = 0; i < picks.size(); ++i)
            for (std::size_t j = 0; j < picks.size(); ++j) {
                bool only_b = picks[i][jb] != picks[j][jb];
                for (std::size_t q = 0; q < act.size() && only_b; ++q)
                    if (q != jb && picks[i][q] != picks[j][q]) only_b = false;
                if (!only_b) continue;
                if (em[i] != em[j]) moves = true;
                for (std::size_t o = 0; o < picks.size(); ++o) {
                    const double d = std::abs(psi(em[i], em[o]) - psi(em[j], em[o]));
                    if (d > worst) {
                        worst = d;
                        wit.tuple_a = picks[i]; wit.tuple_b = picks[j];
                        wit.contribution_a = psi(em[i], em[o]);
                        wit.contribution_b = psi(em[j], em[o]);
                        wit.observed_difference = d;
                        std::size_t na = 0, nbb = 0;
                        for (const SymbolicOrigin& so : R.origins) {
                            bool ca = true, cb = true;
                            for (const auto& c : so.allele_constraint) {
                                const auto it = std::find(act.begin(), act.end(), c.first);
                                if (it == act.end()) continue;
                                const std::size_t q = static_cast<std::size_t>(it - act.begin());
                                if (picks[i][q] != c.second) ca = false;
                                if (picks[j][q] != c.second) cb = false;
                            }
                            na += ca ? 1u : 0u; nbb += cb ? 1u : 0u;
                        }
                        wit.origins_a = na; wit.origins_b = nbb;
                    }
                }
            }
        R.block_delta[b] = worst;
        if (moves) {
            R.structural_scope.push_back(b);
            R.witnesses.push_back(wit);
        }
    }
    std::sort(R.uncertified_blocks.begin(), R.uncertified_blocks.end());
    R.uncertified_blocks.erase(std::unique(R.uncertified_blocks.begin(),
                                           R.uncertified_blocks.end()),
                               R.uncertified_blocks.end());
    return R;
}

} // namespace panvar

namespace panvar {
namespace {

// log T_LS for the ordered diploid transition. The two homologues copy independently:
//   t(i -> i') = (1 - r) * [i == i'] + r / n_hap
// and T((i,j) -> (i2,j2)) = t(i->i2) * t(j->j2). Written out per pair here rather than using the
// factorised recursion, because a linkage edge is not separable and the oracle must see exactly the
// same potential the recursion uses.
inline double ls_log_t(std::size_t a, std::size_t b, std::size_t nh, double r) {
    const double v = (a == b ? (1.0 - r) : 0.0) + r / static_cast<double>(nh);
    return std::log(v);
}

inline double edge_log_phi(const HybridChain& c, std::size_t b,
                           std::size_t i, std::size_t j, std::size_t i2, std::size_t j2) {
    const std::size_t nh = c.n_hap;
    double lp = ls_log_t(i, i2, nh, c.recomb) + ls_log_t(j, j2, nh, c.recomb);
    const HybridEdge& e = c.edges[b];
    if (!e.has_linkage) return lp;
    // THE MAPPING. Homologue 1 takes (allele at A of i, allele at B of i2); homologue 2 takes
    // (allele at A of j, allele at B of j2). In that order, no implicit swap.
    const std::size_t a1 = e.allele_a[i],  b1 = e.allele_b[i2];
    const std::size_t a2 = e.allele_a[j],  b2 = e.allele_b[j2];
    const std::size_t cfg = ((a1 * e.n_b + b1) * e.n_a + a2) * e.n_b + b2;
    return lp + e.log_psi[cfg];
}

}  // namespace

HybridPosterior hybrid_forward_backward(const HybridChain& c) {
    HybridPosterior out;
    const std::size_t nh = c.n_hap, nb = c.n_blocks, ns = nh * nh;
    if (nh == 0 || nb == 0 || c.log_emission.size() != nb || c.edges.size() != nb) return out;

    // THROUGH THE SHARED KERNEL, not a second recursion. This path used to carry its own generic
    // log-space forward-backward; keeping it would have left two implementations of one recursion
    // with only the oracle to notice them drifting apart. The kernel works in scaled probabilities,
    // so each block's emissions are shifted by their own maximum and the shifts are added back into
    // the partition -- a per-block constant, which cancels in the normalisation.
    std::vector<double> shift(nb, 0.0);
    for (std::size_t b = 0; b < nb; ++b) {
        double m = kNegInf;
        for (double v : c.log_emission[b]) m = std::max(m, v);
        shift[b] = std::isfinite(m) ? m : 0.0;
    }
    std::vector<ChainEdgeLinkage> edges(nb);
    for (std::size_t b = 0; b < nb; ++b) {
        const HybridEdge& e = c.edges[b];
        if (!e.has_linkage) continue;
        edges[b].active = true;
        edges[b].n_a = e.n_a; edges[b].n_b = e.n_b;
        edges[b].allele_a = e.allele_a;
        edges[b].allele_b = e.allele_b;
        edges[b].log_psi = e.log_psi;
    }
    std::vector<std::vector<double>> fwd, bwd;
    ChainKernelStats st;
    chain_forward_backward(nh, nb, c.recomb,
                           [&](std::size_t b, std::vector<double>& ev) {
                               ev.assign(ns, 0.0);
                               for (std::size_t k = 0; k < ns; ++k) {
                                   const double v = c.log_emission[b][k];
                                   ev[k] = (v == kNegInf) ? 0.0 : std::exp(v - shift[b]);
                               }
                           },
                           &edges, fwd, bwd, &st);
    if (fwd.size() != nb) return out;
    out.log_partition_unnormalised = st.log_weight_sum;
    for (std::size_t b = 0; b < nb; ++b) out.log_partition_unnormalised += shift[b];
    out.factorised_edges = st.factorised_edges;
    out.linked_edges = st.linked_edges;
    out.log_marginal.assign(nb, std::vector<double>(ns, kNegInf));
    for (std::size_t b = 0; b < nb; ++b) {
        double z = 0.0;
        for (std::size_t k = 0; k < ns; ++k) z += fwd[b][k] * bwd[b][k];
        for (std::size_t k = 0; k < ns; ++k) {
            const double v = fwd[b][k] * bwd[b][k];
            out.log_marginal[b][k] = (v > 0.0 && z > 0.0) ? std::log(v / z) : kNegInf;
        }
    }
    out.ok = true;
    return out;
}


// ---------------------------------------------------------------------------------------------
// THE HIGHER-ORDER RECURRENCE
//
// The message carries, per homologue, the REFINED CLASS of every block whose run has already
// CLOSED, plus the current template. The open run needs no history -- its classes follow from the
// current template -- and the run start is implicit in how many blocks are closed.
//
// THE REFINEMENT IS NOT OPTIONAL. A block inside two factors gets two class partitions, and the
// history must hold their COMMON REFINEMENT: storing either factor's partition alone would discard
// information the other one needs, and a common refinement is never coarser than either input.
// ---------------------------------------------------------------------------------------------
namespace {

struct HoPlan {
    std::size_t nh = 0, nb = 0;
    std::vector<std::size_t> hist_block;                  // blocks needing history, ascending
    std::vector<std::size_t> hist_pos;                    // block -> index in hist_block, or npos
    std::vector<std::vector<std::uint32_t>> refined;      // [block][hap] -> refined class id
    std::vector<std::size_t> n_refined;                   // [block]
    // For factor f, block position j: refined class -> that factor's own class id.
    std::vector<std::vector<std::vector<std::uint32_t>>> to_factor;
    std::vector<std::vector<std::size_t>> ends_at;        // block -> factor indices ending there
    std::size_t last_factor_block = 0;
    std::vector<std::size_t> radix;                       // per hist position, n_refined
    // AFTER block b, the first history position any UNAPPLIED factor still needs. Everything
    // below it is DEAD: its factor has already consumed it, and carrying it multiplies the state
    // space by a product nothing will ever read. Positions below this are zeroed so the states
    // that differed only there MERGE.
    std::vector<std::size_t> base_after;
    std::vector<std::uint64_t> offset;                    // offset[k] = first index with k closed
    std::uint64_t per_hom = 0;                            // total per-homologue states
    std::uint64_t per_hom_hist = 0;                       // the same with the template removed
    bool ok = false;
    std::string refusal;
};

std::size_t ho_closed_before(const HoPlan& P, std::size_t b) {
    std::size_t k = 0;
    while (k < P.hist_block.size() && P.hist_block[k] < b) ++k;
    return k;
}

HoPlan ho_build_plan(const HybridChain& c) {
    HoPlan P;
    P.nh = c.n_hap; P.nb = c.n_blocks;
    if (c.hap_allele.size() != c.n_hap) { P.refusal = "hap_allele missing"; return P; }
    for (const auto& row : c.hap_allele)
        if (row.size() != c.n_blocks) { P.refusal = "hap_allele row wrong length"; return P; }
    P.ends_at.assign(c.n_blocks, {});
    P.n_refined.assign(c.n_blocks, 1);
    P.refined.assign(c.n_blocks, std::vector<std::uint32_t>(c.n_hap, 0));
    P.to_factor.resize(c.higher.size());

    for (std::size_t f = 0; f < c.higher.size(); ++f) {
        const auto& H = c.higher[f];
        if (H.blocks.empty() || !H.usable()) {
            P.refusal = "higher factor " + std::to_string(f) + " is unusable"; return P;
        }
        for (std::size_t j = 0; j + 1 < H.blocks.size(); ++j)
            if (H.blocks[j] + 1 != H.blocks[j + 1]) {
                P.refusal = "higher factor " + std::to_string(f) + " span is not contiguous";
                return P;
            }
        if (H.blocks.back() >= c.n_blocks) {
            P.refusal = "higher factor " + std::to_string(f) + " runs past the chain"; return P;
        }
        P.ends_at[H.blocks.back()].push_back(f);
        P.last_factor_block = std::max(P.last_factor_block, std::size_t(H.blocks.back()));
    }

    // The common refinement, per block, built from the tuple of every factor's class there.
    for (std::size_t b = 0; b < c.n_blocks; ++b) {
        std::map<std::vector<std::uint32_t>, std::uint32_t> seen;
        bool touched = false;
        for (std::size_t t = 0; t < c.n_hap; ++t) {
            std::vector<std::uint32_t> key;
            for (std::size_t f = 0; f < c.higher.size(); ++f) {
                const auto& H = c.higher[f];
                for (std::size_t j = 0; j < H.blocks.size(); ++j)
                    if (H.blocks[j] == b) {
                        touched = true;
                        key.push_back(H.class_of(j, c.hap_allele[t][b]));
                    }
            }
            auto it = seen.find(key);
            if (it == seen.end())
                it = seen.emplace(key, static_cast<std::uint32_t>(seen.size())).first;
            P.refined[b][t] = it->second;
        }
        P.n_refined[b] = touched ? seen.size() : 1;
        if (!touched) for (std::size_t t = 0; t < c.n_hap; ++t) P.refined[b][t] = 0;
    }

    // refined -> factor class. Well defined precisely BECAUSE the refinement refines each factor.
    for (std::size_t f = 0; f < c.higher.size(); ++f) {
        const auto& H = c.higher[f];
        P.to_factor[f].assign(H.blocks.size(), {});
        for (std::size_t j = 0; j < H.blocks.size(); ++j) {
            const std::size_t b = H.blocks[j];
            P.to_factor[f][j].assign(P.n_refined[b], 0);
            std::vector<char> set(P.n_refined[b], 0);
            for (std::size_t t = 0; t < c.n_hap; ++t) {
                const std::uint32_t rc = P.refined[b][t];
                const std::uint32_t fc = H.class_of(j, c.hap_allele[t][b]);
                if (set[rc] && P.to_factor[f][j][rc] != fc) {
                    P.refusal = "the refinement does not refine factor " + std::to_string(f);
                    return P;
                }
                P.to_factor[f][j][rc] = fc; set[rc] = 1;
            }
        }
    }

    P.hist_pos.assign(c.n_blocks, static_cast<std::size_t>(-1));
    for (std::size_t b = 0; b <= P.last_factor_block && b < c.n_blocks; ++b)
        if (P.n_refined[b] > 1) { P.hist_pos[b] = P.hist_block.size(); P.hist_block.push_back(b); }

    P.radix.resize(P.hist_block.size());
    for (std::size_t i = 0; i < P.hist_block.size(); ++i) P.radix[i] = P.n_refined[P.hist_block[i]];
    const std::size_t K = P.hist_block.size();
    P.base_after.assign(c.n_blocks, K);
    for (std::size_t b = 0; b < c.n_blocks; ++b) {
        std::size_t lo = K;
        for (std::size_t f = 0; f < c.higher.size(); ++f) {
            const auto& H = c.higher[f];
            if (H.blocks.back() <= b) continue;      // this factor is done with its span
            for (std::uint32_t q : H.blocks) {
                const std::size_t i = P.hist_pos[q];
                if (i != static_cast<std::size_t>(-1)) lo = std::min(lo, i);
            }
        }
        P.base_after[b] = lo;
    }
    P.offset.assign(P.hist_block.size() + 2, 0);
    long double prod = 1.0L, acc = 0.0L;
    for (std::size_t k = 0; k <= P.hist_block.size(); ++k) {
        P.offset[k] = static_cast<std::uint64_t>(acc);
        acc += prod * static_cast<long double>(P.nh);
        if (acc > 4.0e18L) { P.refusal = "per-homologue state space overflows"; return P; }
        if (k < P.radix.size()) prod *= static_cast<long double>(P.radix[k]);
    }
    P.offset[P.hist_block.size() + 1] = static_cast<std::uint64_t>(acc);
    P.per_hom = static_cast<std::uint64_t>(acc);
    P.per_hom_hist = P.per_hom / P.nh;
    if (P.per_hom == 0 || P.per_hom > 4.0e9) {
        P.refusal = "per-homologue state space too large (" + std::to_string(P.per_hom) + ")";
        return P;
    }
    P.ok = true;
    return P;
}

// One homologue's state, packed. `k` closed blocks, their digits, and the current template.
struct HoHom { std::size_t k = 0; std::uint64_t digits = 0; std::uint32_t t = 0; };

inline std::uint64_t ho_pack(const HoPlan& P, const HoHom& h) {
    return P.offset[h.k] + h.digits * P.nh + h.t;
}
inline HoHom ho_unpack(const HoPlan& P, std::uint64_t idx) {
    HoHom h;
    std::size_t k = 0;
    while (k + 1 <= P.hist_block.size() && idx >= P.offset[k + 1]) ++k;
    const std::uint64_t rel = idx - P.offset[k];
    h.k = k; h.t = static_cast<std::uint32_t>(rel % P.nh); h.digits = rel / P.nh;
    return h;
}
// Digit for hist position i, given `k` closed blocks and the packed digit word.
inline std::uint32_t ho_digit(const HoPlan& P, const HoHom& h, std::size_t i) {
    std::uint64_t d = h.digits;
    for (std::size_t q = 0; q < i; ++q) d /= P.radix[q];
    return static_cast<std::uint32_t>(d % P.radix[i]);
}
// A HISTORY-ONLY index: the same packing with the template dimension removed. The switch
// components aggregate over the old template, so their group key must not mention it.
inline std::uint64_t ho_hist_index(const HoPlan& P, const HoHom& h) {
    return P.offset[h.k] / P.nh + h.digits;      // offset[k] is a multiple of nh by construction
}
inline HoHom ho_from_hist(const HoPlan& P, std::uint64_t idx) {
    HoHom h; std::size_t k = 0;
    while (k + 1 <= P.hist_block.size() && idx >= P.offset[k + 1] / P.nh) ++k;
    h.k = k; h.digits = idx - P.offset[k] / P.nh; return h;
}
// DROP THE DEAD HISTORY. Positions below `base` belong to factors that have already applied, so
// nothing will read them again; zeroing them collapses every state that differed only there. The
// closed count moves up with them, or states that differ only in how much dead history they had
// closed would stay apart while meaning the same thing.
inline void ho_project(const HoPlan& P, HoHom& h, std::size_t b) {
    const std::size_t base = P.base_after[b];
    if (base == 0) return;
    std::uint64_t low = 1;
    for (std::size_t i = 0; i < base && i < P.radix.size(); ++i) low *= P.radix[i];
    h.digits = (h.digits / low) * low;          // zero every dead digit, keep the live ones in place
    if (h.k < base) h.k = std::min(base, P.hist_block.size());
}

// The refined class of block `b` for a homologue: history when closed, current template otherwise.
inline std::uint32_t ho_class_at(const HoPlan& P, const HoHom& h, std::size_t b) {
    const std::size_t i = P.hist_pos[b];
    if (i != static_cast<std::size_t>(-1) && i < h.k) return ho_digit(P, h, i);
    return P.refined[b][h.t];
}

}  // namespace


namespace {
// Overflow-checked arithmetic. A resource plan that silently wraps is worse than no plan: it
// reports a small number for something that will not fit.
inline bool ck_mul(std::uint64_t a, std::uint64_t b, std::uint64_t& out) {
    if (a != 0 && b > std::numeric_limits<std::uint64_t>::max() / a) return false;
    out = a * b; return true;
}
inline bool ck_add(std::uint64_t a, std::uint64_t b, std::uint64_t& out) {
    if (a > std::numeric_limits<std::uint64_t>::max() - b) return false;
    out = a + b; return true;
}
}  // namespace

HigherOrderPlan plan_higher_order(const HybridChain& c, std::size_t threads) {
    HigherOrderPlan pl;
    pl.threads = threads == 0 ? 1 : threads;
    const std::size_t nh = c.n_hap, nb = c.n_blocks;
    if (nh == 0 || nb == 0 || c.log_emission.size() != nb || c.edges.size() != nb) {
        pl.refusal = "malformed chain"; return pl;
    }
    const HoPlan P = ho_build_plan(c);
    if (!P.ok) { pl.refusal = P.refusal; return pl; }
    pl.n_blocks = nb; pl.n_hap = nh;
    pl.refinement_classes = P.n_refined;

    const double r = c.recomb;
    const bool has_stay = (1.0 - r) != 0.0, has_switch = (r / static_cast<double>(nh)) != 0.0;

    bool overflow = false;
    const auto mul = [&](std::uint64_t a, std::uint64_t b) {
        std::uint64_t o = 0; if (!ck_mul(a, b, o)) overflow = true; return o;
    };
    const auto add = [&](std::uint64_t& acc, std::uint64_t v) {
        std::uint64_t o = 0; if (!ck_add(acc, v, o)) overflow = true; else acc = o;
    };

    // THE HAPLOID REACHABILITY PASS. At most per_hom states, which is the product of the
    // refinement sizes times the panel -- enumerable even where the diploid message is not.
    std::vector<std::uint64_t> R;
    R.reserve(nh);
    for (std::uint32_t t = 0; t < nh; ++t) { HoHom h; h.t = t; R.push_back(ho_pack(P, h)); }
    std::sort(R.begin(), R.end());
    pl.haploid_states.assign(nb, 0);
    pl.message_entries.assign(nb, 0);
    pl.haploid_states[0] = R.size();
    pl.message_entries[0] = mul(R.size(), R.size());
    pl.peak_message_entries = pl.message_entries[0];
    pl.total_message_entries = pl.message_entries[0];
    if (!P.ends_at[0].empty())
        add(pl.factor_lookups, mul(pl.message_entries[0], P.ends_at[0].size()));

    std::uint64_t peak_temp_entries = 0;
    std::vector<std::uint64_t> nextR, Hset;
    for (std::size_t b = 1; b < nb; ++b) {
        const std::uint64_t Rp = R.size();
        // The extended histories, and the total number of digits the extension appends.
        Hset.clear();
        std::uint64_t D = 0;
        const std::size_t kb = ho_closed_before(P, b);
        for (std::uint64_t key : R) {
            const HoHom h = ho_unpack(P, key);
            HoHom e = h;
            for (std::size_t i = h.k; i < kb; ++i) {
                std::uint64_t place = 1;
                for (std::size_t q = 0; q < i; ++q) place *= P.radix[q];
                e.digits += static_cast<std::uint64_t>(P.refined[P.hist_block[i]][h.t]) * place;
                ++D;
            }
            e.k = std::max(h.k, kb);
            Hset.push_back(ho_hist_index(P, e));
        }
        std::sort(Hset.begin(), Hset.end());
        Hset.erase(std::unique(Hset.begin(), Hset.end()), Hset.end());
        const std::uint64_t H = has_switch ? Hset.size() : 0;

        // The destinations: stay keeps the state, switch frees the template over every panel entry.
        nextR.clear();
        if (has_stay) for (std::uint64_t key : R) nextR.push_back(key);
        if (has_switch)
            for (std::uint64_t hidx : Hset) {
                HoHom e = ho_from_hist(P, hidx);
                for (std::uint32_t n = 0; n < nh; ++n) { e.t = n; nextR.push_back(ho_pack(P, e)); }
            }
        for (std::uint64_t& key : nextR) {
            HoHom h = ho_unpack(P, key); ho_project(P, h, b); key = ho_pack(P, h);
        }
        std::sort(nextR.begin(), nextR.end());
        nextR.erase(std::unique(nextR.begin(), nextR.end()), nextR.end());

        // ---- THE COUNTS, in closed form from Rp, H and D ------------------------------------
        std::uint64_t upd = 0;
        if (has_stay) add(upd, mul(Rp, Rp));
        if (has_switch && has_stay) add(upd, mul(mul(2ull, mul(H, Rp)), nh));
        if (has_switch) add(upd, mul(mul(H, H), mul(nh, nh)));
        add(pl.forward_updates, upd);
        add(pl.adjoint_updates, upd);                       // the gather visits the same set
        add(pl.multiplier_reconstructions, mul(2ull, upd)); // one per sweep
        if (!P.ends_at[b].empty())
            add(pl.factor_lookups, mul(mul(2ull, upd), P.ends_at[b].size()));
        add(pl.dense_equivalent_updates, mul(mul(Rp, Rp), mul(nh + 1, nh + 1)));
        // Aggregate accumulations: three per source forward, three more in the adjoint's key pass.
        if (has_switch) add(pl.grouping_ops, mul(6ull, mul(Rp, Rp)));
        // History appends: twice per diploid source forward, four times in the adjoint (its key
        // pass and its broadcast each rebuild the keys rather than remembering them).
        if (has_switch) add(pl.history_ops, mul(6ull, mul(Rp, D)));

        std::uint64_t temp = 0;
        if (has_switch && has_stay) add(temp, mul(2ull, mul(H, Rp)));
        if (has_switch) add(temp, mul(H, H));
        peak_temp_entries = std::max(peak_temp_entries, temp);

        R.swap(nextR);
        pl.haploid_states[b] = R.size();
        pl.message_entries[b] = mul(R.size(), R.size());
        pl.peak_message_entries = std::max(pl.peak_message_entries, pl.message_entries[b]);
        add(pl.total_message_entries, pl.message_entries[b]);
        if (overflow) { pl.refusal = "the resource plan overflows 64 bits"; return pl; }
    }

    // ---- BYTES ------------------------------------------------------------------------------
    // PAYLOAD is the doubles. CONTAINER is what the hash actually costs around them: a node per
    // entry (key, value, next pointer, cached hash) plus a bucket array, sized here at one bucket
    // per entry. It is allocator-dependent and therefore an ESTIMATE, which is exactly why it is
    // reported apart from the payload rather than folded into it.
    pl.payload_bytes = mul(mul(2ull, pl.total_message_entries), sizeof(double));
    pl.container_bytes = mul(mul(2ull, pl.total_message_entries), 40ull);
    pl.temporary_bytes = mul(peak_temp_entries, 40ull);
    pl.per_thread_bytes = mul(pl.peak_message_entries, sizeof(double));
    std::uint64_t tot = 0;
    add(tot, pl.payload_bytes); add(tot, pl.container_bytes); add(tot, pl.temporary_bytes);
    add(tot, mul(pl.threads, pl.per_thread_bytes));
    pl.total_bytes = tot;
    if (overflow) { pl.refusal = "the resource plan overflows 64 bits"; return pl; }
    pl.ok = true;
    return pl;
}

// ONE CORE, TWO TRANSITIONS. The plan, the packing, the emissions, the factor application, the
// history extension, the marginals and the refusals are SHARED, so the only thing that differs
// between the production recurrence and its oracle is the transition enumeration -- which is the
// thing being compared. Anything else in common would be a second place for them to drift.
static HybridPosterior ho_run(const HybridChain& c, std::uint64_t max_message_entries,
                              HigherOrderStats* stats, bool dense, HigherOrderTrace* trace) {
    HybridPosterior out;
    HigherOrderStats st;
    const std::size_t nh = c.n_hap, nb = c.n_blocks, ns = nh * nh;
    if (nh == 0 || nb == 0 || c.log_emission.size() != nb || c.edges.size() != nb) {
        st.refused = true; st.refusal = "malformed chain";
        if (stats) *stats = st;
        return out;
    }
    const HoPlan P = ho_build_plan(c);
    if (!P.ok) {
        st.refused = true; st.refusal = P.refusal;
        if (stats) *stats = st;
        return out;
    }
    // REFUSAL BEFORE ALLOCATION, from the SAME plan the recurrence then follows. This is not a
    // guard that happens to agree with a separate estimate: the plan is computed here, the budget
    // is checked against it, and the fixtures require the plan's counts to equal the realised
    // ones exactly. Nothing has been allocated at this point beyond the plan itself.
    const HigherOrderPlan PL = plan_higher_order(c, 1);
    if (!PL.ok) {
        st.refused = true; st.refusal = "resource plan: " + PL.refusal;
        if (stats) *stats = st;
        return out;
    }
    st.planned_peak_message_entries = PL.peak_message_entries;
    st.planned_forward_updates = PL.forward_updates;
    st.planned_total_bytes = PL.total_bytes;
    if (max_message_entries != 0 && PL.peak_message_entries > max_message_entries) {
        st.refused = true;
        st.refusal = "planned peak message of " + std::to_string(PL.peak_message_entries) +
                     " entries exceeds max_message_entries (" +
                     std::to_string(max_message_entries) + ")";
        if (stats) *stats = st;
        return out;
    }

    std::vector<double> shift(nb, 0.0);
    for (std::size_t b = 0; b < nb; ++b) {
        double m = kNegInf;
        for (double v : c.log_emission[b]) m = std::max(m, v);
        shift[b] = std::isfinite(m) ? m : 0.0;
    }
    const double r = c.recomb;
    const double w_stay = 1.0 - r, w_switch = r / static_cast<double>(nh);

    // The emission of a diploid state, in scaled probability space.
    auto emit = [&](std::size_t b, std::uint32_t t1, std::uint32_t t2) {
        const double v = c.log_emission[b][static_cast<std::size_t>(t1) * nh + t2];
        return (v == kNegInf) ? 0.0 : std::exp(v - shift[b]);
    };

    // The factors that close at block b, applied to the joint class tuple of the span.
    std::vector<std::uint32_t> ca, cb;
    auto apply_factors = [&](std::size_t b, const HoHom& h1, const HoHom& h2) {
        double m = 1.0;
        for (std::size_t f : P.ends_at[b]) {
            const auto& H = c.higher[f];
            ca.resize(H.blocks.size()); cb.resize(H.blocks.size());
            for (std::size_t j = 0; j < H.blocks.size(); ++j) {
                const std::size_t q = H.blocks[j];
                ca[j] = P.to_factor[f][j][ho_class_at(P, h1, q)];
                cb[j] = P.to_factor[f][j][ho_class_at(P, h2, q)];
            }
            ++st.factor_lookups;
            const double lp = H.log_psi_classes(ca, cb);
            if (lp != 0.0) m *= std::exp(lp);
        }
        return m;
    };

    // Advance one homologue by one block. `sw` selects the switch component; `nt` is the new
    // template. On a switch every block whose run has just closed is appended to the history --
    // that is the only place history grows, and it is counted.
    auto advance = [&](const HoHom& cur, std::size_t b, bool sw, std::uint32_t nt) {
        HoHom nx = cur;
        if (sw) {
            const std::size_t kb = ho_closed_before(P, b);
            for (std::size_t i = cur.k; i < kb; ++i) {
                std::uint64_t place = 1;
                for (std::size_t q = 0; q < i; ++q) place *= P.radix[q];
                nx.digits += static_cast<std::uint64_t>(P.refined[P.hist_block[i]][cur.t]) * place;
                ++st.history_ops;
            }
            nx.k = std::max(cur.k, kb);
            nx.t = nt;
        }
        return nx;
    };

    // ---- FORWARD ------------------------------------------------------------------------------
    const auto t_fwd = std::chrono::steady_clock::now();
    std::vector<std::unordered_map<std::uint64_t, double>> msg(nb);
    for (std::uint32_t t1 = 0; t1 < nh; ++t1)
        for (std::uint32_t t2 = 0; t2 < nh; ++t2) {
            HoHom h1, h2; h1.t = t1; h2.t = t2;
            double w = emit(0, t1, t2);
            if (!P.ends_at[0].empty()) w *= apply_factors(0, h1, h2);
            if (w == 0.0) {
                // Dropped, and RECORDED. The plan counted this state; saying so is what turns
                // "the plan was exact" from a claim into a checkable condition.
                st.all_initial_emissions_finite = false;
                ++st.initial_states_dropped;
                continue;
            }
            msg[0][ho_pack(P, h1) * P.per_hom + ho_pack(P, h2)] += w;
        }
    st.peak_message_entries = msg[0].size();
    st.total_message_entries = msg[0].size();

    // The four transition components, kept apart because Li-Stephens FACTORISES:
    //   T = (1-r) I + (r/n) 11^T.
    // Enumerating (n1, n2) per source would cost O(entering x n_hap^2) -- the dense assumption that
    // makes a real locus four orders of magnitude too expensive. Instead each SWITCH component
    // AGGREGATES over the old template first (the closed run's classes are appended before the sum,
    // which is what lets the sum happen at all) and only then redistributes over the new template.
    const double c_ss = w_stay * w_stay, c_ws = w_switch * w_stay, c_ww = w_switch * w_switch;
    // STAY/STAY NEEDS NO AGGREGATE. Its group key IS the source key -- one source, one
    // destination -- so grouping it would allocate a second copy of the whole message to buy
    // nothing. Only the switch components actually merge sources.
    std::unordered_map<std::uint64_t, double> gWS, gSW, gWW;
    for (std::size_t b = 1; b < nb; ++b) {
        auto& nxt = msg[b];
        gWS.clear(); gSW.clear(); gWW.clear();
        // Counted for BOTH implementations: it is the yardstick, so it cannot live inside the
        // branch that is being measured against it.
        st.dense_equivalent_updates +=
            static_cast<std::uint64_t>(msg[b - 1].size()) *
            (static_cast<std::uint64_t>(nh) + 1) * (static_cast<std::uint64_t>(nh) + 1);
        if (dense) {
        for (const auto& kv : msg[b - 1]) {
            const HoHom c1 = ho_unpack(P, kv.first / P.per_hom);
            const HoHom c2 = ho_unpack(P, kv.first % P.per_hom);
            for (int s1 = 0; s1 < 2; ++s1) for (int s2 = 0; s2 < 2; ++s2) {
                const double w1 = s1 ? w_switch : w_stay, w2 = s2 ? w_switch : w_stay;
                if (w1 == 0.0 || w2 == 0.0) continue;
                const std::uint32_t lo1 = s1 ? 0 : c1.t, hi1 = s1 ? nh : c1.t + 1;
                const std::uint32_t lo2 = s2 ? 0 : c2.t, hi2 = s2 ? nh : c2.t + 1;
                for (std::uint32_t n1 = lo1; n1 < hi1; ++n1)
                for (std::uint32_t n2 = lo2; n2 < hi2; ++n2) {
                    const HoHom a = advance(c1, b, s1 != 0, n1);
                    const HoHom d = advance(c2, b, s2 != 0, n2);
                    double m = w1 * w2 * emit(b, n1, n2);
                    ++st.multiplier_reconstructions;
                    if (!P.ends_at[b].empty()) m *= apply_factors(b, a, d);
                    HoHom a2 = a, d2 = d;
                    ho_project(P, a2, b); ho_project(P, d2, b);
                    nxt[ho_pack(P, a2) * P.per_hom + ho_pack(P, d2)] += kv.second * m;
                    ++st.forward_updates;
                }
            }
        }
        } else {
        // ---- AGGREGATE ---------------------------------------------------------------------
        for (const auto& kv : msg[b - 1]) {
            const HoHom c1 = ho_unpack(P, kv.first / P.per_hom);
            const HoHom c2 = ho_unpack(P, kv.first % P.per_hom);
            // The closed runs are appended ONCE per source, before any aggregation -- that is what
            // makes the switch components summable at all. With r = 0 no run ever closes here and
            // the extension is not even computed.
            HoHom e1 = c1, e2 = c2;
            if (c_ws != 0.0 || c_ww != 0.0) {
                e1 = advance(c1, b, true, 0); e2 = advance(c2, b, true, 0);
            }
            const double w = kv.second;
            if (c_ws != 0.0) {
                gWS[ho_hist_index(P, e1) * P.per_hom + ho_pack(P, c2)] += w * c_ws;
                gSW[ho_pack(P, c1) * P.per_hom_hist + ho_hist_index(P, e2)] += w * c_ws;
                st.grouping_ops += 2;
            }
            if (c_ww != 0.0) {
                gWW[ho_hist_index(P, e1) * P.per_hom_hist + ho_hist_index(P, e2)] += w * c_ww;
                ++st.grouping_ops;
            }
        }
        // ---- REDISTRIBUTE ------------------------------------------------------------------
        const auto place = [&](HoHom a, HoHom d, double g) {
            double m = emit(b, a.t, d.t);
            ++st.multiplier_reconstructions;
            if (!P.ends_at[b].empty()) m *= apply_factors(b, a, d);
            ho_project(P, a, b); ho_project(P, d, b);
            nxt[ho_pack(P, a) * P.per_hom + ho_pack(P, d)] += g * m;
            ++st.forward_updates;
        };
        if (c_ss != 0.0)
            for (const auto& kv : msg[b - 1])
                place(ho_unpack(P, kv.first / P.per_hom), ho_unpack(P, kv.first % P.per_hom),
                      kv.second * c_ss);
        for (const auto& kv : gWS) {
            HoHom a = ho_from_hist(P, kv.first / P.per_hom);
            const HoHom d = ho_unpack(P, kv.first % P.per_hom);
            for (std::uint32_t n1 = 0; n1 < nh; ++n1) { a.t = n1; place(a, d, kv.second); }
        }
        for (const auto& kv : gSW) {
            const HoHom a = ho_unpack(P, kv.first / P.per_hom_hist);
            HoHom d = ho_from_hist(P, kv.first % P.per_hom_hist);
            for (std::uint32_t n2 = 0; n2 < nh; ++n2) { d.t = n2; place(a, d, kv.second); }
        }
        for (const auto& kv : gWW) {
            HoHom a = ho_from_hist(P, kv.first / P.per_hom_hist);
            HoHom d = ho_from_hist(P, kv.first % P.per_hom_hist);
            for (std::uint32_t n1 = 0; n1 < nh; ++n1) {
                a.t = n1;
                for (std::uint32_t n2 = 0; n2 < nh; ++n2) { d.t = n2; place(a, d, kv.second); }
            }
        }
        }
        // The plan said this would fit. The realised message may be SMALLER (dropped initial
        // states propagate), but never larger: the plan is an upper bound always, and exact when
        // every initial emission is finite. Exceeding it means the plan is wrong, and continuing
        // would allocate past a budget that was already approved -- so this refuses
        // transactionally, clearing every message rather than leaving partials.
        if (nxt.size() > PL.message_entries[b]) {
            for (auto& m : msg) m.clear();
            st.refused = true;
            st.refusal = "message at block " + std::to_string(b) + " has " +
                         std::to_string(nxt.size()) + " entries, more than the planned " +
                         std::to_string(PL.message_entries[b]);
            if (stats) *stats = st;
            return out;
        }
        st.peak_message_entries = std::max<std::uint64_t>(st.peak_message_entries, nxt.size());
        st.total_message_entries += nxt.size();
    }
    st.forward_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_fwd).count();

    double Z = 0.0;
    for (const auto& kv : msg[nb - 1]) Z += kv.second;
    if (!(Z > 0.0)) {
        st.refused = true; st.refusal = "zero partition weight";
        if (stats) *stats = st;
        return out;
    }

    // ---- THE TAPE-FREE ADJOINT ----------------------------------------------------------------
    // bar_b(s) = d Z / d msg_b(s). The update  dst += src * m  has adjoint  bar_src += bar_dst * m,
    // so this re-walks the SAME loop and RECOMPUTES m rather than reading it from a tape: one tape
    // entry per update would dominate every other cost at scale.
    //
    // PARALLEL SAFETY: the outer loop is over SOURCES and each iteration accumulates only into its
    // own source's adjoint, so threads over sources need no reduction and no atomic accumulation.
    const auto t_adj = std::chrono::steady_clock::now();
    std::vector<std::unordered_map<std::uint64_t, double>> bar(nb);
    for (const auto& kv : msg[nb - 1]) bar[nb - 1][kv.first] = 1.0;
    // THE ADJOINT OF A FACTORISED STEP. Forward was AGGREGATE then REDISTRIBUTE; reversed, that is
    // GATHER (each group's adjoint from its destinations) then BROADCAST (each source reads the
    // groups it fed). Both halves stay contention-free: the gather writes one entry per GROUP, the
    // broadcast writes one entry per SOURCE, and neither ever writes another's.
    if (dense) {
    for (std::size_t b = nb - 1; b-- > 0;) {
        auto& here = bar[b];
        const auto& next = bar[b + 1];
        for (const auto& kv : msg[b]) {
            const HoHom c1 = ho_unpack(P, kv.first / P.per_hom);
            const HoHom c2 = ho_unpack(P, kv.first % P.per_hom);
            double acc = 0.0;
            for (int s1 = 0; s1 < 2; ++s1) for (int s2 = 0; s2 < 2; ++s2) {
                const double w1 = s1 ? w_switch : w_stay, w2 = s2 ? w_switch : w_stay;
                if (w1 == 0.0 || w2 == 0.0) continue;
                const std::uint32_t lo1 = s1 ? 0 : c1.t, hi1 = s1 ? nh : c1.t + 1;
                const std::uint32_t lo2 = s2 ? 0 : c2.t, hi2 = s2 ? nh : c2.t + 1;
                for (std::uint32_t n1 = lo1; n1 < hi1; ++n1)
                for (std::uint32_t n2 = lo2; n2 < hi2; ++n2) {
                    const HoHom a = advance(c1, b + 1, s1 != 0, n1);
                    const HoHom d = advance(c2, b + 1, s2 != 0, n2);
                    double m = w1 * w2 * emit(b + 1, n1, n2);
                    ++st.multiplier_reconstructions;
                    if (!P.ends_at[b + 1].empty()) m *= apply_factors(b + 1, a, d);
                    HoHom a2 = a, d2 = d;
                    ho_project(P, a2, b + 1); ho_project(P, d2, b + 1);
                    const auto it = next.find(ho_pack(P, a2) * P.per_hom + ho_pack(P, d2));
                    if (it == next.end()) continue;
                    acc += it->second * m;
                    ++st.adjoint_updates;
                }
            }
            here[kv.first] = acc;
        }
    }
    } else {
    std::unordered_map<std::uint64_t, double> aWS, aSW, aWW;
    for (std::size_t b = nb - 1; b-- > 0;) {
        const std::size_t bb = b + 1;
        auto& here = bar[b];
        aWS.clear(); aSW.clear(); aWW.clear();
        // The group KEYS, recomputed from the sources -- a tape is exactly what we are refusing
        // to store, so they are derived again rather than remembered.
        const auto keys_of = [&](std::uint64_t key, std::uint64_t& kWS,
                                 std::uint64_t& kSW, std::uint64_t& kWW) {
            const HoHom c1 = ho_unpack(P, key / P.per_hom);
            const HoHom c2 = ho_unpack(P, key % P.per_hom);
            HoHom e1 = c1, e2 = c2;
            if (c_ws != 0.0 || c_ww != 0.0) {
                e1 = advance(c1, bb, true, 0); e2 = advance(c2, bb, true, 0);
            }
            kWS = ho_hist_index(P, e1) * P.per_hom + ho_pack(P, c2);
            kSW = ho_pack(P, c1) * P.per_hom_hist + ho_hist_index(P, e2);
            kWW = ho_hist_index(P, e1) * P.per_hom_hist + ho_hist_index(P, e2);
        };
        for (const auto& kv : msg[b]) {
            std::uint64_t kWS, kSW, kWW;
            keys_of(kv.first, kWS, kSW, kWW);
            if (c_ws != 0.0) {
                aWS.emplace(kWS, 0.0); aSW.emplace(kSW, 0.0); st.grouping_ops += 2;
            }
            if (c_ww != 0.0) { aWW.emplace(kWW, 0.0); ++st.grouping_ops; }
        }
        const auto& next = bar[bb];
        const auto gather = [&](HoHom a, HoHom d) {
            double m = emit(bb, a.t, d.t);
            ++st.multiplier_reconstructions;
            if (!P.ends_at[bb].empty()) m *= apply_factors(bb, a, d);
            ho_project(P, a, bb); ho_project(P, d, bb);
            const auto it = next.find(ho_pack(P, a) * P.per_hom + ho_pack(P, d));
            ++st.adjoint_updates;
            return it == next.end() ? 0.0 : it->second * m;
        };
        for (auto& kv : aWS) {
            HoHom a = ho_from_hist(P, kv.first / P.per_hom);
            const HoHom d = ho_unpack(P, kv.first % P.per_hom);
            double acc = 0.0;
            for (std::uint32_t n1 = 0; n1 < nh; ++n1) { a.t = n1; acc += gather(a, d); }
            kv.second = acc;
        }
        for (auto& kv : aSW) {
            const HoHom a = ho_unpack(P, kv.first / P.per_hom_hist);
            HoHom d = ho_from_hist(P, kv.first % P.per_hom_hist);
            double acc = 0.0;
            for (std::uint32_t n2 = 0; n2 < nh; ++n2) { d.t = n2; acc += gather(a, d); }
            kv.second = acc;
        }
        for (auto& kv : aWW) {
            HoHom a = ho_from_hist(P, kv.first / P.per_hom_hist);
            HoHom d = ho_from_hist(P, kv.first % P.per_hom_hist);
            double acc = 0.0;
            for (std::uint32_t n1 = 0; n1 < nh; ++n1) {
                a.t = n1;
                for (std::uint32_t n2 = 0; n2 < nh; ++n2) { d.t = n2; acc += gather(a, d); }
            }
            kv.second = acc;
        }
        // BROADCAST. Each source accumulates only into its own adjoint.
        for (const auto& kv : msg[b]) {
            std::uint64_t kWS, kSW, kWW;
            keys_of(kv.first, kWS, kSW, kWW);
            double acc = 0.0;
            if (c_ss != 0.0)
                acc += c_ss * gather(ho_unpack(P, kv.first / P.per_hom),
                                     ho_unpack(P, kv.first % P.per_hom));
            if (c_ws != 0.0) acc += c_ws * (aWS[kWS] + aSW[kSW]);
            if (c_ww != 0.0) acc += c_ww * aWW[kWW];
            here[kv.first] = acc;
        }
    }
    }
    st.adjoint_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_adj).count();

    // marginal(b, x) = SUM over states s with template pair x of msg_b(s) * bar_b(s).
    out.log_marginal.assign(nb, std::vector<double>(ns, kNegInf));
    for (std::size_t b = 0; b < nb; ++b) {
        std::vector<double> mg(ns, 0.0);
        for (const auto& kv : msg[b]) {
            const auto it = bar[b].find(kv.first);
            if (it == bar[b].end()) continue;
            const HoHom h1 = ho_unpack(P, kv.first / P.per_hom);
            const HoHom h2 = ho_unpack(P, kv.first % P.per_hom);
            mg[static_cast<std::size_t>(h1.t) * nh + h2.t] += kv.second * it->second;
        }
        double z = 0.0;
        for (double v : mg) z += v;
        for (std::size_t x = 0; x < ns; ++x)
            out.log_marginal[b][x] = (mg[x] > 0.0 && z > 0.0) ? std::log(mg[x] / z) : kNegInf;
    }
    out.log_partition_unnormalised = std::log(Z);
    for (std::size_t b = 0; b < nb; ++b) out.log_partition_unnormalised += shift[b];
    out.linked_edges = 0;
    out.factorised_edges = nb > 0 ? nb - 1 : 0;
    out.ok = true;
    // PAYLOAD ONLY -- the message and adjoint values. Not a memory total: it excludes the hash
    // index, allocator overhead, the factor tables, temporaries and any reduction buffers.
    st.predicted_payload_bytes = st.total_message_entries * 2ull * sizeof(double);
    // THE TRACE, for fixture comparison only: both implementations use the same HoPlan, so their
    // message keys are directly comparable and a disagreement can be localised to a state rather
    // than only showing up as a marginal that moved.
    if (trace) {
        trace->messages.assign(nb, {});
        trace->adjoints.assign(nb, {});
        for (std::size_t b = 0; b < nb; ++b) {
            for (const auto& kv : msg[b]) trace->messages[b].push_back(kv);
            for (const auto& kv : bar[b]) trace->adjoints[b].push_back(kv);
            std::sort(trace->messages[b].begin(), trace->messages[b].end());
            std::sort(trace->adjoints[b].begin(), trace->adjoints[b].end());
        }
    }
    if (stats) *stats = st;
    return out;
}

HybridPosterior hybrid_higher_order(const HybridChain& c, std::uint64_t max_message_entries,
                                    HigherOrderStats* stats, HigherOrderTrace* trace) {
    return ho_run(c, max_message_entries, stats, /*dense=*/false, trace);
}

HybridPosterior hybrid_higher_order_dense_oracle(const HybridChain& c,
                                                 std::uint64_t max_message_entries,
                                                 HigherOrderStats* stats,
                                                 HigherOrderTrace* trace) {
    return ho_run(c, max_message_entries, stats, /*dense=*/true, trace);
}

HybridPosterior hybrid_bruteforce(const HybridChain& c) {
    HybridPosterior out;
    const std::size_t nh = c.n_hap, nb = c.n_blocks, ns = nh * nh;
    if (nh == 0 || nb == 0 || c.log_emission.size() != nb || c.edges.size() != nb) return out;
    // ns^nb paths. Fixtures only; refuse rather than hang if someone points it at a real chain.
    double total_paths = 1.0;
    for (std::size_t b = 0; b < nb; ++b) total_paths *= static_cast<double>(ns);
    if (total_paths > 5.0e7) return out;

    out.log_marginal.assign(nb, std::vector<double>(ns, kNegInf));
    double z = kNegInf;
    std::vector<std::size_t> path(nb, 0);
    for (;;) {
        // Score this path from scratch: emissions at every block, plus one edge potential per edge.
        // Nothing is reused from the recursion, which is what makes this an independent check.
        double lp = 0.0;
        bool dead = false;
        for (std::size_t b = 0; b < nb && !dead; ++b) {
            const double e = c.log_emission[b][path[b]];
            if (e == kNegInf) { dead = true; break; }
            lp += e;
        }
        if (!dead) {
            for (std::size_t b = 1; b < nb; ++b) {
                const std::size_t i = path[b - 1] / nh, j = path[b - 1] % nh;
                const std::size_t i2 = path[b] / nh, j2 = path[b] % nh;
                // COMPUTED INDEPENDENTLY, not through edge_log_phi. Sharing that helper would make
                // a bug inside it -- the ordered-state mapping above all -- move both arms
                // identically and cancel out of the comparison, which is exactly the class of error
                // this oracle exists to catch.
                const double r = c.recomb;
                const double ti = (i == i2 ? (1.0 - r) : 0.0) + r / static_cast<double>(nh);
                const double tj = (j == j2 ? (1.0 - r) : 0.0) + r / static_cast<double>(nh);
                lp += std::log(ti) + std::log(tj);
                const HybridEdge& e = c.edges[b];
                if (e.has_linkage) {
                    // Homologue 1: allele of i at the LEFT block, allele of i2 at the RIGHT.
                    // Homologue 2: allele of j at the LEFT block, allele of j2 at the RIGHT.
                    const std::size_t h1a = e.allele_a[i],  h1b = e.allele_b[i2];
                    const std::size_t h2a = e.allele_a[j],  h2b = e.allele_b[j2];
                    lp += e.log_psi[((h1a * e.n_b + h1b) * e.n_a + h2a) * e.n_b + h2b];
                }
            }
            // THE HIGHER FACTORS, read off the whole path directly -- no message, no history, no
            // class refinement. Everything the recurrence carries state for is here a plain index
            // into the path, which is what makes the comparison independent.
            for (const HybridHigherFactor& H : c.higher) {
                if (!H.usable()) continue;
                std::vector<std::uint32_t> k1(H.blocks.size()), k2(H.blocks.size());
                for (std::size_t j = 0; j < H.blocks.size(); ++j) {
                    const std::size_t b = H.blocks[j];
                    k1[j] = H.class_of(j, c.hap_allele[path[b] / nh][b]);
                    k2[j] = H.class_of(j, c.hap_allele[path[b] % nh][b]);
                }
                lp += H.log_psi_classes(k1, k2);
            }
            z = log_add(z, lp);
            for (std::size_t b = 0; b < nb; ++b) {
                out.log_marginal[b][path[b]] = log_add(out.log_marginal[b][path[b]], lp);
            }
        }
        std::size_t k = 0;
        for (; k < nb; ++k) { if (++path[k] < ns) break; path[k] = 0; }
        if (k == nb) break;
    }
    out.log_partition_unnormalised = z;
    for (std::size_t b = 0; b < nb; ++b) {
        for (std::size_t s = 0; s < ns; ++s) {
            if (out.log_marginal[b][s] != kNegInf) out.log_marginal[b][s] -= z;
        }
    }
    out.ok = true;
    return out;
}

const char* unusable_reason_name(UnusableReason r) {
    switch (r) {
        case UnusableReason::EmptyRead:            return "empty_read";
        case UnusableReason::NoInBandOrigins:      return "no_in_band_origins";
        case UnusableReason::ScopeNotCertified:    return "scope_not_certified";
        case UnusableReason::EmptyScope:           return "empty_scope";
        case UnusableReason::UnmappedMassDominant: return "unmapped_mass_dominant";
        default:                                   return "none";
    }
}

}  // namespace panvar

