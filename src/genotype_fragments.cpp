#include "panvar/genotype_fragments.hpp"

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
double insert_ll(double implied, const FragmentScoreOptions& o) {
    const double z = (implied - o.fragment_len) / o.fragment_sd;
    const double concordant = -0.5 * z * z - std::log(o.fragment_sd) - 0.9189385332046727;
    const double discordant = std::log(o.discordant_rate) - std::log(o.discordant_span);
    return log_add(std::log1p(-o.discordant_rate) + concordant, discordant);
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
                    lp += insert_ll(static_cast<double>(hi - lo + 1), options);
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
    const std::size_t nh = std::min(options.max_haplotypes, ranked.size());
    std::vector<HaplotypeSeq> haps(nh);
    out.haplotypes.resize(nh);
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
            if (per_hap[o.code] > options.max_anchor_occ) continue;
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

    std::vector<double> ll(fragments.size() * nh, kNegInf);
    std::vector<double> floors(fragments.size(), kNegInf);
    // Where each fragment landed on each haplotype, or -1 where it did not land at all. The depth
    // channel is built from this: a haplotype carrying sequence the sample does not have shows up as
    // a run of windows with no fragment in them.
    std::vector<std::int32_t> midpoint(fragments.size() * nh, -1);

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

        const auto reduce = [&](std::vector<Cand>& c) {
            // Bucket implied starts to 64 bp and keep, per haplotype, the two commonest.
            std::map<std::pair<std::uint32_t, long>, std::pair<int, long>> tally;
            for (const Cand& x : c) {
                auto& e = tally[{x.hap, (x.fwd ? 1 : -1) * (x.start / 64 + 1)}];
                ++e.first;
                e.second = x.start;
            }
            std::map<std::uint32_t, std::vector<std::pair<int, std::pair<bool, long>>>> per_hap;
            for (const auto& [key, val] : tally) {
                per_hap[key.first].push_back({val.first, {key.second > 0, val.second}});
            }
            std::vector<Cand> keep;
            for (auto& [hi, v] : per_hap) {
                std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
                for (std::size_t i = 0; i < std::min<std::size_t>(2, v.size()); ++i) {
                    keep.push_back({hi, v[i].second.first, v[i].second.second});
                }
            }
            c.swap(keep);
        };
        reduce(c1);
        reduce(c2);

        const std::size_t band1 =
            static_cast<std::size_t>(options.max_divergence * static_cast<double>(F.r1.size())) + 1;
        const std::size_t band2 = F.r2.empty() ? 1 :
            static_cast<std::size_t>(options.max_divergence * static_cast<double>(F.r2.size())) + 1;

        struct Placed { bool ok = false; std::size_t edits = 0; long start = 0; long end = 0; };
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
                const std::string window = haps[hi].seq.substr(static_cast<std::size_t>(lo),
                                                               static_cast<std::size_t>(hi_end - lo));
                const ReadFit f = infix_align(q, window, band);
                if (!f.ok) continue;
                Placed p;
                p.ok = true;
                p.edits = f.edits;
                p.start = lo + static_cast<long>(f.start);
                p.end = lo + static_cast<long>(f.end);
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
            double lp = kNegInf;
            long mid = -1;
            const double miss1 = read_ll(band1, F.r1.size());
            const double miss2 = F.r2.empty() ? 0.0 : read_ll(band2, F.r2.size());
            const auto consider = [&](double v, long m) { if (v > lp) { lp = v; mid = m; } };
            if (!a1.empty() && !a2.empty()) {
                for (const Placed& x : a1) {
                    for (const Placed& y : a2) {
                        double v = read_ll(x.edits, F.r1.size()) + read_ll(y.edits, F.r2.size());
                        if (options.use_insert_size) {
                            v += insert_ll(static_cast<double>(std::max(x.end, y.end) -
                                                               std::min(x.start, y.start) + 1), options);
                        }
                        consider(v, (std::min(x.start, y.start) + std::max(x.end, y.end)) / 2);
                    }
                }
            }
            for (const Placed& x : a1) consider(read_ll(x.edits, F.r1.size()) + miss2,
                                                (x.start + x.end) / 2);
            for (const Placed& y : a2) consider(miss1 + read_ll(y.edits, F.r2.size()),
                                                (y.start + y.end) / 2);
            if (lp == kNegInf) { ll[fi * nh + hi] = floors[fi]; continue; }
            midpoint[fi * nh + hi] = static_cast<std::int32_t>(mid);
            ll[fi * nh + hi] = lp;
        }
        for (std::size_t hi = 0; hi < nh; ++hi) {
            if (ll[fi * nh + hi] == kNegInf) ll[fi * nh + hi] = floors[fi];
        }
    });

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
        pair_score[pi] = total + cov_ll[a] + cov_ll[b];
    });
    for (std::size_t pi = 0; pi < pair_index.size(); ++pi) {
        pairs.push_back({pair_index[pi].first, pair_index[pi].second, pair_score[pi], 0.0});
    }
    std::stable_sort(pairs.begin(), pairs.end(),
                     [](const HaplotypePairScore& x, const HaplotypePairScore& y) {
                         return x.score > y.score;
                     });

    double norm = kNegInf;
    for (const HaplotypePairScore& p : pairs) norm = log_add(norm, p.score);
    for (HaplotypePairScore& p : pairs) p.posterior = std::exp(p.score - norm);

    // ---- project onto blocks ------------------------------------------------------------------
    out.blocks.resize(chain.size());
    for (std::size_t bi = 0; bi < chain.size(); ++bi) {
        BlockProjection& P = out.blocks[bi];
        P.block_index = bi;
        P.kind = chain[bi].kind;
        P.bubble_id = chain[bi].bubble_id;
        P.n_alleles = blocks[bi].allele_seq.size();
        // Marginal, not a read-out of the best pair: haplotype pairs that disagree elsewhere may
        // agree here, and that agreement is evidence.
        std::map<std::pair<int, int>, double> mass;
        for (const HaplotypePairScore& p : pairs) {
            if (p.posterior < 1e-12) continue;
            int a = haps[p.hap1].allele[bi];
            int b = haps[p.hap2].allele[bi];
            if (a > b) std::swap(a, b);
            mass[{a, b}] += p.posterior;
        }
        double best = -1.0;
        for (const auto& [key, m] : mass) {
            if (m > best) { best = m; P.allele1 = key.first; P.allele2 = key.second; }
        }
        P.posterior = best < 0.0 ? 0.0 : best;
        if (truth_allele1 != nullptr && truth_allele2 != nullptr &&
            bi < truth_allele1->size() && bi < truth_allele2->size()) {
            P.truth_a = (*truth_allele1)[bi];
            P.truth_b = (*truth_allele2)[bi];
            if (P.truth_a >= 0 && P.truth_b >= 0) {
                P.truth_representable = true;
                const int ta = std::min(P.truth_a, P.truth_b);
                const int tb = std::max(P.truth_a, P.truth_b);
                P.exact = (ta == P.allele1 && tb == P.allele2);
            }
        }
    }

    for (std::size_t i = 0; i < std::min(top_pairs_kept, pairs.size()); ++i) {
        out.top_pairs.push_back(pairs[i]);
    }
    return out;
}

void write_haplotype_results(const std::string& out_prefix,
                             const HaplotypeResult& result,
                             bool have_truth) {
    const std::string bp = out_prefix + ".hap_blocks.tsv";
    std::ofstream bf(bp);
    if (!bf) throw std::runtime_error("genotype-frag: cannot write " + bp);
    bf << "block\tkind\tbubble_id\tn_alleles\tallele1\tallele2\tposterior";
    if (have_truth) bf << "\ttruth_a\ttruth_b\trepresentable\texact";
    bf << '\n';
    for (const BlockProjection& p : result.blocks) {
        const char* kind = p.kind == BlockKind::Bubble ? "bubble"
                         : p.kind == BlockKind::Backbone ? "backbone" : "flank";
        bf << p.block_index << '\t' << kind << '\t' << p.bubble_id << '\t' << p.n_alleles << '\t'
           << p.allele1 << '\t' << p.allele2 << '\t' << p.posterior;
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

} // namespace panvar
