#include "panvar/genotype_fragments.hpp"

#include "panvar/align.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/parallel.hpp"
#include "panvar/syncmer.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include <zlib.h>
#include <kseq.h>

#include "edlib.h"

KSEQ_INIT(gzFile, gzread)

namespace panvar {

double InsertPrior::exposure(std::size_t hap_len) const {
    const long n = static_cast<long>(hap_len);
    double e = 0.0;
    for (long L = lo; L <= hi; ++L) {
        const double starts = static_cast<double>(std::max<long>(0, n - L + 1));
        if (starts > 0.0) e += std::exp(log_at(L)) * starts;
    }
    return e;
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
    return ip;
}

namespace {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();

double log_add(double a, double b) {
    if (a == kNegInf) return b;
    if (b == kNegInf) return a;
    const double hi = a > b ? a : b;
    const double lo = a > b ? b : a;
    return hi + std::log1p(std::exp(lo - hi));
}

// wgsim and `samtools fastq` both give the two mates one shared name with a /1 or /2 suffix, so the
// pairing can be recovered from the name alone -- which means interleaved input and split R1/R2
// files go through the same code and neither has to be declared.
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

    // ---- recruitment index: syncmer code -> which target blocks can explain it ----------------
    // A fragment is NOT forced onto one block. At cyp2d6 the same syncmer occurs in blocks 3 and 5,
    // and the marker model's answer to that is to delete the marker (over-expected) or the read's
    // contribution to it (confinement). Here the fragment is offered to both blocks and its
    // alignment decides -- because a fragment carrying a shared syncmer usually also carries a
    // block-specific one, and it is that co-occurrence, not a filter, that localises the evidence.
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> recruit;
    recruit.reserve(1u << 20);
    {
        std::vector<std::uint64_t> codes;
        for (std::size_t t = 0; t < tg.size(); ++t) {
            codes.clear();
            const BlockAlleles& b = blocks[tg[t].block];
            const auto add = [&](const std::string& seq) {
                for (const KmerOccurrence& o : collect_syncmers(seq, k, s)) codes.push_back(o.code);
            };
            for (const std::string& a : b.allele_seq) add(a);
            add(tg[t].lf);
            add(tg[t].rf);
            std::sort(codes.begin(), codes.end());
            codes.erase(std::unique(codes.begin(), codes.end()), codes.end());
            for (const std::uint64_t c : codes) recruit[c].push_back(static_cast<std::uint32_t>(t));
        }
    }

    std::vector<std::vector<std::uint32_t>> recruited(tg.size());
    {
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
                if (n >= options.min_recruit_hits) recruited[t].push_back(static_cast<std::uint32_t>(fi));
            }
        }
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

        const double log_eps = std::log(options.error_rate);
        const double log_1meps = std::log1p(-options.error_rate);
        const double log_mix = std::log1p(-options.outlier_mix);
        const double log_out = std::log(options.outlier_mix);
        long min_frag_len_b = 1;
        for (const std::uint32_t fi : recruited[t]) {
            min_frag_len_b = std::max<long>(min_frag_len_b,
                static_cast<long>(fragments[fi].r1.size() + fragments[fi].r2.size()));
        }
        const InsertPrior ins_prior = make_insert_prior(options.fragment_len, options.fragment_sd,
                                                        options.discordant_rate, 4, min_frag_len_b);
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
            const std::size_t band1 =
                static_cast<std::size_t>(options.max_divergence * static_cast<double>(F.r1.size())) + 1;
            const std::size_t band2 =
                static_cast<std::size_t>(options.max_divergence * static_cast<double>(F.r2.size())) + 1;

            const ReadFit p_fwd = infix_align(F.r1, T.contexts[probe], band1);
            const ReadFit p_rev = infix_align(r1rc, T.contexts[probe], band1);
            const bool r1_forward =
                (p_fwd.ok ? p_fwd.edits : band1 + 1) <= (p_rev.ok ? p_rev.edits : band1 + 1);
            const std::string& q1 = r1_forward ? F.r1 : r1rc;
            const std::string& q2 = r1_forward ? r2rc : F.r2;   // mates are on opposite strands

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

void verify_block_spelling(const Graph& graph,
                           const std::vector<BlockAlleles>& blocks,
                           const std::vector<std::string>& haplotype_names) {
    const auto by_name = path_records_by_name(graph);
    for (const std::string& name : haplotype_names) {
        const auto it = by_name.find(name);
        if (it == by_name.end() || it->second == nullptr) continue;
        bool complete = false;
        const std::string raw = spell_path_steps_sequence(graph, it->second->steps, &complete);
        if (!complete) continue;                    // a path the graph cannot spell is not our claim
        const HaplotypeSeq built = spell_haplotype(blocks, name);
        if (built.seq == raw) continue;
        std::size_t at = 0;
        while (at < built.seq.size() && at < raw.size() && built.seq[at] == raw[at]) ++at;
        throw std::runtime_error(
            "genotype-frag: block decomposition does not round-trip for path '" + name +
            "': concatenated block alleles are " + std::to_string(built.seq.size()) +
            " bp, the graph spells " + std::to_string(raw.size()) +
            " bp, first difference at offset " + std::to_string(at) +
            ". Whole-haplotype mode would score a sequence the panel does not contain.");
    }
}

HaplotypeResult genotype_haplotype_pairs(
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    const std::vector<std::string>& haplotype_names,
    const std::vector<Fragment>& fragments,
    const HaplotypeScoreOptions& options,
    const std::vector<int>* truth_allele1,
    const std::vector<int>* truth_allele2,
    std::size_t top_pairs_kept) {

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
    std::vector<std::pair<double, std::size_t>> ranked;
    ranked.reserve(haplotype_names.size());
    {
        std::vector<double> containment(haplotype_names.size(), 0.0);
        run_parallel(haplotype_names.size(), options.threads, [&](std::size_t i) {
            const HaplotypeSeq h = spell_haplotype(blocks, haplotype_names[i]);
            const std::vector<KmerOccurrence> sy = collect_syncmers(h.seq, k, s);
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
            const HaplotypeSeq h = spell_haplotype(blocks, haplotype_names[r.second]);
            if (seen.insert(h.seq).second) unique_ranked.push_back(r);
        }
        ranked.swap(unique_ranked);
    }
    std::size_t nh = std::min(options.max_haplotypes, ranked.size());
    while (nh > 0 && nh < ranked.size() &&
           std::abs(ranked[nh].first - ranked[nh - 1].first) < 1e-12) {
        ++nh;   // the cut fell inside a tie; taking one side of it would be arbitrary
    }
    std::vector<HaplotypeSeq> haps(nh);
    out.haplotypes.resize(nh);
    // Reported so a run can say whether the shortlist was a real selection or a tie it could not
    // break. `nh` above the requested cap means it was extended through a tie.
    for (std::size_t i = 0; i < nh; ++i) {
        out.shortlist.push_back(haplotype_names[ranked[i].second]);
        haps[i] = spell_haplotype(blocks, haplotype_names[ranked[i].second]);
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

    const double log_eps = std::log(options.error_rate);
    const double log_1meps = std::log1p(-options.error_rate);
    const auto read_ll = [&](std::size_t edits, std::size_t len) {
        return static_cast<double>(edits) * log_eps +
               static_cast<double>(len - std::min(edits, len)) * log_1meps;
    };
    long min_frag_len = 1;
    for (const Fragment& f : fragments) {
        min_frag_len = std::max<long>(min_frag_len, static_cast<long>(f.r1.size() + f.r2.size()));
    }
    const InsertPrior ins_prior = make_insert_prior(options.fragment_len, options.fragment_sd,
                                                    options.discordant_rate, 4, min_frag_len);
    // Explicit 1/2 per strand, matching the reference. It is a per-fragment constant only while the
    // placement term is compared with itself; mixed against a background it changes the placement
    // scale and therefore the contrast between candidates.
    const double log_half_strand = std::log(0.5);

    std::vector<double> ll(fragments.size() * nh, kNegInf);
    std::vector<double> floors(fragments.size(), kNegInf);
    // Where each fragment landed on each haplotype, or -1 where it did not land at all. The depth
    // channel is built from this: a haplotype carrying sequence the sample does not have shows up as
    // a run of windows with no fragment in them.
    std::vector<std::int32_t> midpoint(fragments.size() * nh, -1);
    // EVERY distinct placement, not just the best one. A fragment compatible with several copies of a
    // repeat is evidence for a haplotype offering several, and pinning it to one arbitrary copy
    // leaves the others falsely empty -- which then charges the candidate for absence the placement
    // heuristic invented. Deduplicated by rounded midpoint so two syncmer anchors reaching the same
    // biological placement count once. CSR-packed: a dense array would be 48 haplotypes x every
    // fragment x every placement.
    // Orientation is part of the state. Two strands with identical coordinates are DIFFERENT states
    // under the contract and must not be merged when exact states are being retained.
    struct PlacementRec { std::int32_t mid, start, end; bool fwd; double ll; };
    std::vector<std::vector<PlacementRec>> placements(
        options.joint_depth ? fragments.size() * nh : 0);

    run_parallel(fragments.size(), options.threads, [&](std::size_t fi) {
        const Fragment& F = fragments[fi];
        const std::size_t total_len = F.bases();
        if (total_len == 0) return;
        floors[fi] = read_ll(
            static_cast<std::size_t>(options.bg_divergence * static_cast<double>(total_len)), total_len);

        // Implied read start per (haplotype, orientation), gathered from the read's own syncmers.
        // Clustering on the implied START rather than on the match position is what lets one anchor
        // stand in for the whole read.
        struct Cand { std::uint32_t hap; bool fwd; long start; };
        const auto gather = [&](const std::string& r, bool fwd,
                                std::vector<Cand>& into) {
            for (const KmerOccurrence& o : collect_syncmers(r, k, s)) {
                const auto it = anchors.find(o.code);
                if (it == anchors.end()) continue;
                for (const auto& [hi, pos] : it->second) {
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

        std::uint64_t found_local = 0, kept_local = 0;
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
        {
            static std::mutex comp_mu;
            std::lock_guard<std::mutex> lk(comp_mu);
            out.completeness.clusters_found += found_local;
            out.completeness.clusters_kept += kept_local;
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
        for (const Cand& c : c1) touched.push_back(c.hap);
        for (const Cand& c : c2) touched.push_back(c.hap);
        std::sort(touched.begin(), touched.end());
        touched.erase(std::unique(touched.begin(), touched.end()), touched.end());

        for (const std::uint32_t hi : touched) {
            // Both mates are placed JOINTLY, over every combination of their candidate anchors,
            // rather than each taking its own best position. Independently, inside a duplication the
            // two mates settle on different copies and the pair then looks 13 kb long -- the mates
            // are one observation and have to be placed as one.
            const std::vector<Placed> a1 = place_all(c1, hi, F.r1, r1rc, band1);
            const std::vector<Placed> a2 = F.r2.empty() ? std::vector<Placed>{}
                                                        : place_all(c2, hi, F.r2, r2rc, band2);
            double lp = kNegInf;      // the accumulated likelihood: max, or the sum when marginalising
            double best = kNegInf;    // always the best single placement, so the midpoint is a real one
            long mid = -1;
            const double miss1 = read_ll(band1, F.r1.size());
            const double miss2 = F.r2.empty() ? 0.0 : read_ll(band2, F.r2.size());
            const auto consider = [&](double v, long m) {
                if (v > best) { best = v; mid = m; }
                lp = options.marginalise_placements ? log_add(lp, v) : std::max(lp, v);
            };
            if (!a1.empty() && !a2.empty()) {
                // Both mates placed: the pair's placements are the COMBINATIONS, and the single-mate
                // terms below would double count them, so they are skipped.
                for (const Placed& x : a1) {
                    for (const Placed& y : a2) {
                        // FR only: the two mates of a fragment face each other. Combining same-strand
                        // placements invents fragments the library cannot produce.
                        if (x.fwd == y.fwd) continue;
                        const Placed& fw = x.fwd ? x : y;
                        const Placed& rv = x.fwd ? y : x;
                        if (rv.end < fw.start) continue;          // reverse mate must lie downstream
                        double v = read_ll(x.edits, F.r1.size()) + read_ll(y.edits, F.r2.size())
                                 + log_half_strand;
                        if (options.use_insert_size) {
                            v += insert_ll(static_cast<double>(rv.end - fw.start + 1), ins_prior);
                        }
                        consider(v, (fw.start + rv.end) / 2);
                    }
                }
            } else {
                for (const Placed& x : a1) consider(read_ll(x.edits, F.r1.size()) + miss2,
                                                    (x.start + x.end) / 2);
                for (const Placed& y : a2) consider(miss1 + read_ll(y.edits, F.r2.size()),
                                                    (y.start + y.end) / 2);
            }
            if (lp == kNegInf) { ll[fi * nh + hi] = floors[fi]; continue; }
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
                if (!a1.empty() && !a2.empty()) {
                    for (const Placed& x : a1) {
                        for (const Placed& y : a2) {
                            if (x.fwd == y.fwd) continue;
                            const Placed& fw = x.fwd ? x : y;
                            const Placed& rv = x.fwd ? y : x;
                            if (rv.end < fw.start) continue;
                            double v = read_ll(x.edits, F.r1.size()) + read_ll(y.edits, F.r2.size())
                                     + log_half_strand;
                            if (options.use_insert_size) {
                                v += insert_ll(static_cast<double>(rv.end - fw.start + 1), ins_prior);
                            }
                            add((fw.start + rv.end) / 2, fw.start, rv.end, x.fwd, v);
                        }
                    }
                } else {
                    for (const Placed& x : a1) add((x.start + x.end) / 2, x.start, x.end, x.fwd,
                                                   read_ll(x.edits, F.r1.size()) + miss2);
                    for (const Placed& y : a2) add((y.start + y.end) / 2, y.start, y.end, y.fwd,
                                                   miss1 + read_ll(y.edits, F.r2.size()));
                }
            }
        }
        for (std::size_t hi = 0; hi < nh; ++hi) {
            if (ll[fi * nh + hi] == kNegInf) ll[fi * nh + hi] = floors[fi];
        }
    });

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
                        const std::int32_t mid = pr.mid; const double v = pr.ll;
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
            std::ofstream mf(options.dump_fragment_mass);
            if (mf) {
                mf << "# pair\t" << out.shortlist[a] << '\t' << out.shortlist[b] << '\n';
                mf << "fragment\tlog_mass\n";
                for (std::size_t fi = 0; fi < fragments.size(); ++fi) {
                    double lse = kNegInf;
                    for (const auto& pr : placements[fi * nh + a]) lse = log_add(lse, pr.ll);
                    if (!homoz) {
                        for (const auto& pr : placements[fi * nh + b]) lse = log_add(lse, pr.ll);
                    } else if (lse != kNegInf) {
                        lse += std::log(2.0);
                    }
                    mf << fragments[fi].name << '\t' << lse << '\n';
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
            bool agree = true;
            int fa = -2, fb = -2;
            for (const HaplotypePairScore& p : pairs) {
                if (best - p.score > options.equivalence_tolerance) break;
                int x = haps[p.hap1].allele[bi], y = haps[p.hap2].allele[bi];
                if (x > y) std::swap(x, y);
                if (fa == -2) { fa = x; fb = y; }
                else if (x != fa || y != fb) { agree = false; break; }
            }
            P.determined = agree;
            if (agree) ++out.equivalence.blocks_determined;
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
                             bool have_truth) {
    const std::string bp = out_prefix + ".hap_blocks.tsv";
    std::ofstream bf(bp);
    if (!bf) throw std::runtime_error("genotype-frag: cannot write " + bp);
    bf << "block\tkind\tbubble_id\tn_alleles\tallele1\tallele2\tposterior\tdetermined";
    if (have_truth) bf << "\ttruth_a\ttruth_b\trepresentable\texact";
    bf << '\n';
    for (const BlockProjection& p : result.blocks) {
        const char* kind = p.kind == BlockKind::Bubble ? "bubble"
                         : p.kind == BlockKind::Backbone ? "backbone" : "flank";
        bf << p.block_index << '\t' << kind << '\t' << p.bubble_id << '\t' << p.n_alleles << '\t'
           << p.allele1 << '\t' << p.allele2 << '\t' << p.posterior << '\t'
           << (p.determined ? 1 : 0);
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

} // namespace

double reference_pair_loglik(const std::string& hap_a, const std::string& hap_b,
                             const std::vector<Fragment>& fragments,
                             const ReferenceParams& params,
                             std::vector<double>* fragment_mass) {
    const double log_eps = std::log(params.error_rate);
    const double log_1meps = std::log1p(-params.error_rate);

    // One prior, shared by the event term and the exposure. Its lower bound is the longest fragment
    // present, because an insert shorter than the two mates is not a state at all.
    long min_len = 1;
    for (const Fragment& f : fragments) {
        min_len = std::max<long>(min_len, static_cast<long>(f.r1.size() + f.r2.size()));
    }
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
        total += log_add(placed, log_bg_w + bg);
    }
    return total;
}

} // namespace panvar
