#include "panvar/gene_cn_kmer.hpp"

#include "panvar/align.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/minimap2_align.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <functional>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace panvar {
namespace {

inline int base_code(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1;   // N / other -> resets the rolling k-mer
    }
}

// Collect every canonical (min of forward / reverse-complement) k-mer of `seq` into `out`.
void collect_canonical_kmers(const std::string& seq, std::size_t k,
                             std::unordered_set<std::uint64_t>& out) {
    if (k == 0 || k > 31) return;
    const std::uint64_t mask = (1ULL << (2 * k)) - 1;
    const unsigned top_shift = static_cast<unsigned>(2 * (k - 1));
    std::uint64_t fwd = 0, rev = 0;
    std::size_t len = 0;
    for (char c : seq) {
        const int b = base_code(c);
        if (b < 0) { len = 0; fwd = 0; rev = 0; continue; }
        fwd = ((fwd << 2) | static_cast<std::uint64_t>(b)) & mask;
        rev = (rev >> 2) | (static_cast<std::uint64_t>(3 - b) << top_shift);
        if (++len >= k) out.insert(std::min(fwd, rev));
    }
}

// As above, but keeping how many times each canonical k-mer occurs: copy number is an occurrence
// ratio, so a set would lose the multiplicity the estimate is made of.
void count_canonical_kmers(const std::string& seq, std::size_t k,
                           std::unordered_map<std::uint64_t, long>& out) {
    if (k == 0 || k > 31 || seq.size() < k) return;
    const std::uint64_t mask = (1ULL << (2 * k)) - 1;
    const unsigned top_shift = static_cast<unsigned>(2 * (k - 1));
    std::uint64_t fwd = 0, rev = 0;
    std::size_t len = 0;
    for (char c : seq) {
        const int b = base_code(c);
        if (b < 0) { len = 0; fwd = 0; rev = 0; continue; }
        fwd = ((fwd << 2) | static_cast<std::uint64_t>(b)) & mask;
        rev = (rev >> 2) | (static_cast<std::uint64_t>(3 - b) << top_shift);
        if (++len >= k) ++out[std::min(fwd, rev)];
    }
}

// A run of hits is a candidate copy when it holds this share of the gene's private k-mers: a version of the gene
// can keep few of them, since one base differing from the reference removes up to k.
constexpr double kCopyMinShare = 0.2;
// A candidate is a copy when the gene aligns there at this identity or more, and the stretch it matched is
// closer to it than to any other gene of the locus. A related sequence that took a stretch of the gene's
// alleles lies further from it.
constexpr double kCopyMinIdentity = 0.98;

} // namespace

std::vector<std::vector<GeneCnEvidence>> resolve_gene_cn_kmer(
    const std::vector<std::string>& gene_markers,
    const std::vector<std::string>& hap_seqs,
    std::size_t k,
    const std::vector<std::string>* backgrounds,
    const std::vector<std::string>* spans,
    const std::vector<char>* read_runs) {

    const std::size_t n_genes = gene_markers.size();
    std::vector<std::vector<GeneCnEvidence>> out(
        hap_seqs.size(), std::vector<GeneCnEvidence>(n_genes));
    if (n_genes == 0 || k == 0 || k > 31) return out;

    // 1) canonical k-mer counts per gene. Counts, not a set: a k-mer occurring twice inside one copy
    //    of a gene gives two hits from that copy, so the denominator below must count it twice too.
    std::vector<std::unordered_map<std::uint64_t, long>> gene_counts(n_genes);
    for (std::size_t g = 0; g < n_genes; ++g)
        count_canonical_kmers(gene_markers[g], k, gene_counts[g]);

    // 2) private k-mers: absent from every sibling's marker and, when `backgrounds` is given, from the
    //    siblings' full spans, since a haplotype carrying a sibling carries its introns and UTRs and
    //    would count them as copies of this gene. Map private k-mer -> owning gene.
    //
    // The rest of the locus is not screened: a second, unannotated copy of the gene in the reference
    // would remove every marker the gene has.
    std::unordered_map<std::uint64_t, int> owner;
    std::vector<long> priv_size(n_genes, 0);
    for (std::size_t g = 0; g < n_genes; ++g) {
        std::unordered_set<std::uint64_t> background;
        if (backgrounds != nullptr && g < backgrounds->size())
            collect_canonical_kmers((*backgrounds)[g], k, background);
        for (const auto& [km, in_gene] : gene_counts[g]) {
            bool shared = background.count(km) != 0;
            for (std::size_t o = 0; o < n_genes && !shared; ++o)
                if (o != g && gene_counts[o].count(km)) shared = true;
            if (shared) continue;
            owner[km] = static_cast<int>(g);
            // The denominator counts occurrences per reference copy, so hits/denominator is the
            // number of copies; counting distinct k-mers would read a k-mer repeated within one copy
            // as extra copies of the gene.
            priv_size[g] += in_gene;
        }
    }

    // Each gene's sequence on the genome, both strands, for the identity of a run of hits.
    std::vector<std::string> span_fwd(n_genes), span_rev(n_genes);
    for (std::size_t g = 0; g < n_genes; ++g) {
        span_fwd[g] = spans != nullptr && g < spans->size() && !(*spans)[g].empty() ? (*spans)[g] : gene_markers[g];
        span_rev[g] = reverse_complement(span_fwd[g]);
    }
    const bool debug_runs = std::getenv("PANVAR_DEBUG_RUNS") != nullptr;

    // 3) per haplotype: count private-k-mer occurrences, derive dosage + CN.
    const std::uint64_t mask = (1ULL << (2 * k)) - 1;
    const unsigned top_shift = static_cast<unsigned>(2 * (k - 1));
    for (std::size_t h = 0; h < hap_seqs.size(); ++h) {
        std::vector<long> hits(n_genes, 0);
        std::vector<std::vector<std::size_t>> at(n_genes);   // where each hit lies on the haplotype
        std::uint64_t fwd = 0, rev = 0;
        std::size_t len = 0;
        const std::string& hs = hap_seqs[h];
        for (std::size_t i = 0; i < hs.size(); ++i) {
            const int b = base_code(hs[i]);
            if (b < 0) { len = 0; fwd = 0; rev = 0; continue; }
            fwd = ((fwd << 2) | static_cast<std::uint64_t>(b)) & mask;
            rev = (rev >> 2) | (static_cast<std::uint64_t>(3 - b) << top_shift);
            if (++len >= k) {
                const auto it = owner.find(std::min(fwd, rev));
                if (it != owner.end()) { ++hits[it->second]; at[it->second].push_back(i); }
            }
        }
        for (std::size_t g = 0; g < n_genes; ++g) {
            GeneCnEvidence& e = out[h][g];
            e.priv_kmers = priv_size[g];
            e.hits = hits[g];
            e.separable = priv_size[g] > 0;
            e.dosage = e.separable ? static_cast<double>(hits[g]) / static_cast<double>(priv_size[g]) : 0.0;
            e.cn = e.separable ? static_cast<int>(std::llround(e.dosage)) : 0;
            e.rounded = e.cn;
            if (!e.separable) continue;
            // Runs of hits: one copy's hits lie within the gene's length on the genome, so hits further apart
            // than that (and 10 kb) start a new run. A run holding a fifth of the private k-mers or more is a
            // candidate copy; the alignment decides (see kCopyMinIdentity).
            if (spans == nullptr || (read_runs != nullptr && !(*read_runs)[g])) continue;
            const std::size_t gene_len = span_fwd[g].size();
            const std::size_t gap = std::max<std::size_t>(10000, gene_len);
            std::vector<std::pair<std::size_t, std::size_t>> taken;   // the copies found, on the haplotype
            const auto check_run = [&](std::size_t j0, std::size_t j1) {   // hits [j0, j1) of gene g
                const double share = static_cast<double>(j1 - j0) / static_cast<double>(priv_size[g]);
                const std::size_t lo = at[g][j0], hi = at[g][j1 - 1];
                // The gene, on either strand, aligned anywhere a copy holding the run can lie: within the gene's
                // length (and a tenth, for indels) of both its ends. A run longer than that (hits of other
                // transcripts' exons, stray hits) is widened by the gene's length on both sides instead.
                const std::size_t reach = gene_len + gene_len / 10;
                std::size_t wlo = hi > reach ? hi - reach : 0, whi = std::min(hs.size() - 1, lo + reach);
                if (hi - lo + 1 > reach) { wlo = lo > reach ? lo - reach : 0; whi = std::min(hs.size() - 1, hi + reach); }
                const std::string window = hs.substr(wlo, whi - wlo + 1);
                // Edits allowed: a tenth of the gene, enough for a copy with long indels (each counts once in the
                // identity, but every base in the edit distance edlib bounds).
                const std::size_t bound = gene_len / 10;
                std::size_t s0 = 0, e0 = 0, s1 = 0, e1 = 0;
                const double idf = infix_identity(span_fwd[g], window, bound, s0, e0);
                const double idr = infix_identity(span_rev[g], window, bound, s1, e1);
                const bool fwd = idf >= idr;
                const double identity = std::max(idf, idr);
                // The stretch it matched, aligned to the other genes of the locus (their reference sequence with
                // this gene's own masked, since an annotation can hold another whole): one matching it as well
                // claims it -- a paralog's copy, or a longer gene holding a homologous piece.
                double other = -1.0;
                bool copy = false;
                // Where it matched, on the haplotype: it must hold at least half the run, and no other run's copy.
                const std::size_t p0 = wlo + (fwd ? s0 : s1), p1 = wlo + (fwd ? e0 : e1);
                const auto overlap = [](std::size_t a0, std::size_t a1, std::size_t b0, std::size_t b1) {
                    return std::min(a1, b1) >= std::max(a0, b0) ? std::min(a1, b1) - std::max(a0, b0) + 1 : std::size_t{0};
                };
                bool own = identity >= 0.0 && 2 * overlap(p0, p1, lo, hi) >= hi - lo + 1;
                for (const auto& [c0, c1] : taken)
                    if (own && 2 * overlap(p0, p1, c0, c1) > c1 - c0 + 1) own = false;
                if (identity >= 0.0) {
                    const std::string piece = fwd ? window.substr(s0, e0 - s0 + 1) : reverse_complement(window.substr(s1, e1 - s1 + 1));
                    if (backgrounds != nullptr && g < backgrounds->size() && !(*backgrounds)[g].empty()) {
                        std::size_t b0 = 0, b1 = 0;
                        other = std::max(infix_identity(piece, (*backgrounds)[g], piece.size() / 10, b0, b1),
                                         infix_identity(reverse_complement(piece), (*backgrounds)[g], piece.size() / 10, b0, b1));
                    }
                    copy = own && identity >= kCopyMinIdentity && other < identity;
                }
                if (copy) taken.emplace_back(p0, p1);
                if (debug_runs) {
                    const char* why = copy ? "copy" : identity < 0.0 ? "no_placement" : !own ? "not_this_run"
                                    : identity < kCopyMinIdentity ? "identity" : "other_gene";
                    if (FILE* fh = std::fopen(std::getenv("PANVAR_DEBUG_RUNS"), "a")) {
                        std::fprintf(fh, "%zu\t%zu\t%ld\t%zu\t%zu\t%zu\t%.4f\t%.4f\t%.4f\t%s\n", h, g, priv_size[g], lo, hi,
                                     j1 - j0, share, identity, other, why);
                        std::fclose(fh);
                    }
                }
                if (copy && e.copies >= 0) ++e.copies;
            };
            std::vector<std::pair<std::size_t, std::size_t>> candidates;
            std::size_t j0 = 0;
            for (std::size_t j = 1; j <= at[g].size(); ++j)
                if (j == at[g].size() || at[g][j] - at[g][j - 1] > gap) {
                    const double share = static_cast<double>(j - j0) / static_cast<double>(priv_size[g]);
                    if (share >= kCopyMinShare) candidates.emplace_back(j0, j);
                    else if (debug_runs)   // too small to be a copy: logged, not aligned
                        if (FILE* fh = std::fopen(std::getenv("PANVAR_DEBUG_RUNS"), "a")) {
                            std::fprintf(fh, "%zu\t%zu\t%ld\t%zu\t%zu\t%zu\t%.4f\tNA\tNA\tshare\n", h, g, priv_size[g], at[g][j0],
                                         at[g][j - 1], j - j0, share);
                            std::fclose(fh);
                        }
                    j0 = j;
                }
            // Aligned only where the candidates could raise the count; the debug log aligns them all, counting none
            // of the others.
            const bool could_raise = candidates.size() > static_cast<std::size_t>(e.cn);
            if (gene_len > 0 && (debug_runs || could_raise)) {
                if (could_raise) e.copies = 0;
                for (const auto& [a, b] : candidates) check_run(a, b);
            }
        }
    }
    return out;
}

namespace {

// Canonical 2-bit code of the length-k window at seq[i..i+k) (assumes ACGT; caller checks). Uses the same
// canonicalisation as collect_canonical_kmers so codes are directly comparable across the two.
inline std::uint64_t canon_code(const std::string& seq, std::size_t i, std::size_t k) {
    const std::uint64_t mask = (1ULL << (2 * k)) - 1;
    const unsigned top = static_cast<unsigned>(2 * (k - 1));
    std::uint64_t fwd = 0, rev = 0;
    for (std::size_t t = 0; t < k; ++t) {
        const int b = base_code(seq[i + t]);
        fwd = ((fwd << 2) | static_cast<std::uint64_t>(b)) & mask;
        rev = (rev >> 2) | (static_cast<std::uint64_t>(3 - b) << top);
    }
    return std::min(fwd, rev);
}

// Canonical k-mers of every window overlapping [lo, hi) in `seq` (a divergent site's span).
void collect_site_windows(const std::string& seq, std::size_t k, std::size_t lo, std::size_t hi,
                          std::unordered_set<std::uint64_t>& out) {
    if (seq.size() < k) return;
    const std::size_t last = seq.size() - k;
    const std::size_t start = (lo >= k - 1) ? lo - (k - 1) : 0;
    const std::size_t end = std::min(hi, last);
    for (std::size_t i = start; i <= end; ++i) {
        bool ok = true;
        for (std::size_t t = 0; t < k; ++t) if (base_code(seq[i + t]) < 0) { ok = false; break; }
        if (ok) out.insert(canon_code(seq, i, k));
    }
}

struct PairSiteEvidence {
    long sites = 0;
    int cn_a = 0, cn_b = 0;
    long hits_a = 0, hits_b = 0, priv_a = 0, priv_b = 0;
    // Whether this haplotype produced a usable split. `sites` and `priv_*` describe the marker sets,
    // which exist for every haplotype once the pair has divergent sites at all, so neither can answer
    // it; without this a haplotype with no usable site would publish cn_a = cn_b = 0 as a call.
    long usable_sites = 0;
    bool resolved = false;
};

// Per-site consensus split of a near-identical gene pair (a, b). Aligns b (query) to a (target) once to
// locate divergent sites, builds each site's A/B allele k-mers (private to the pair vs `other`), then per
// haplotype takes the median across sites of the site-local A fraction and splits total[hap] by it. Returns
// sites==0 when the pair has no alignment or no site with private k-mers on both sides, so the caller
// keeps the pooled dosage.
std::vector<PairSiteEvidence> resolve_pair_cn_persite(
    const std::string& marker_a, const std::string& marker_b,
    const std::unordered_set<std::uint64_t>& other,
    const std::vector<std::string>& hap_seqs,
    const std::vector<long>& total,
    std::size_t k) {

    std::vector<PairSiteEvidence> out(hap_seqs.size());
    if (marker_a.empty() || marker_b.empty()) return out;
    // asm5: a near-identical pair (that is why it is routed here) aligns tightly; the low-divergence preset
    // places the divergent columns most accurately, which is what the per-site split reads.
    const Minimap2Hit h = minimap2_best_hit("b", marker_b, "a", marker_a, "asm5");
    if (!h.ok || h.cigar.empty()) return out;

    std::unordered_set<std::uint64_t> aset, bset;
    collect_canonical_kmers(marker_a, k, aset);
    collect_canonical_kmers(marker_b, k, bset);

    // On a reverse-strand hit the CIGAR walks the target forward against the query's reverse
    // complement, from qlen - query_end; walking marker_b forward would read each site's B side at a
    // mirrored position and lose an inverted pair's sites. Canonical k-mers are the same on either
    // strand, so only the coordinates change.
    const std::string b_seq = h.reverse ? reverse_complement(marker_b) : marker_b;

    struct Site { std::vector<std::uint64_t> a_km, b_km; };
    std::vector<Site> sites;
    std::size_t tpos = h.target_start_bp;                                      // a (target)
    std::size_t qpos = h.reverse ? marker_b.size() - h.query_end_bp : h.query_start_bp;   // b (query)
    for (const auto& oplen : h.cigar) {
        const int op = oplen.first, len = oplen.second;
        if (op == 7) { tpos += len; qpos += len; continue; }         // '=' match
        std::size_t a_lo = tpos, a_hi = tpos, b_lo = qpos, b_hi = qpos;
        if (op == 8) { a_hi = tpos + len; b_hi = qpos + len; tpos += len; qpos += len; }   // 'X' sub
        else if (op == 1) { b_hi = qpos + len; qpos += len; }        // 'I' extra in b
        else if (op == 2) { a_hi = tpos + len; tpos += len; }        // 'D' extra in a
        else continue;                                               // clips etc.
        std::unordered_set<std::uint64_t> aw, bw;
        collect_site_windows(marker_a, k, a_lo, a_hi, aw);
        collect_site_windows(b_seq, k, b_lo, b_hi, bw);
        Site s;
        for (std::uint64_t km : aw) if (!bset.count(km) && !other.count(km)) s.a_km.push_back(km);
        for (std::uint64_t km : bw) if (!aset.count(km) && !other.count(km)) s.b_km.push_back(km);
        if (!s.a_km.empty() && !s.b_km.empty()) sites.push_back(std::move(s));
    }
    if (sites.empty()) return out;

    // k-mer -> every (site, 0=a / 1=b) it belongs to. Divergent columns closer than k bp share
    // windows, and every such site counts the k-mer in its denominator, so every one must get its
    // hits too. Close columns stay separate sites, each with its own vote, so a converted column is
    // not merged with an unconverted one into a site that outvotes the median.
    std::unordered_map<std::uint64_t, std::vector<std::pair<int, int>>> loc;
    for (std::size_t si = 0; si < sites.size(); ++si) {
        for (std::uint64_t km : sites[si].a_km) loc[km].emplace_back(static_cast<int>(si), 0);
        for (std::uint64_t km : sites[si].b_km) loc[km].emplace_back(static_cast<int>(si), 1);
    }

    const std::uint64_t mask = (1ULL << (2 * k)) - 1;
    const unsigned top = static_cast<unsigned>(2 * (k - 1));
    for (std::size_t hp = 0; hp < hap_seqs.size(); ++hp) {
        std::vector<long> ah(sites.size(), 0), bh(sites.size(), 0);
        std::uint64_t fwd = 0, rev = 0; std::size_t len = 0;
        for (char c : hap_seqs[hp]) {
            const int b = base_code(c);
            if (b < 0) { len = 0; fwd = 0; rev = 0; continue; }
            fwd = ((fwd << 2) | static_cast<std::uint64_t>(b)) & mask;
            rev = (rev >> 2) | (static_cast<std::uint64_t>(3 - b) << top);
            if (++len >= k) {
                const auto it = loc.find(std::min(fwd, rev));
                if (it != loc.end())
                    for (const auto& [si, side] : it->second) { if (side == 0) ++ah[si]; else ++bh[si]; }
            }
        }
        std::vector<double> fr;
        PairSiteEvidence& e = out[hp];
        e.sites = static_cast<long>(sites.size());
        for (std::size_t si = 0; si < sites.size(); ++si) {
            const long asz = static_cast<long>(sites[si].a_km.size());
            const long bsz = static_cast<long>(sites[si].b_km.size());
            e.hits_a += ah[si]; e.hits_b += bh[si]; e.priv_a += asz; e.priv_b += bsz;
            if (asz > 0 && bsz > 0) {
                const double na = static_cast<double>(ah[si]) / static_cast<double>(asz);
                const double nb = static_cast<double>(bh[si]) / static_cast<double>(bsz);
                if (na + nb >= 0.5) { fr.push_back(na / (na + nb)); ++e.usable_sites; }
            }
        }
        e.resolved = !fr.empty();
        if (!fr.empty()) {
            std::sort(fr.begin(), fr.end());
            const double med = (fr.size() % 2) ? fr[fr.size() / 2]
                                               : 0.5 * (fr[fr.size() / 2 - 1] + fr[fr.size() / 2]);
            const long tot = hp < total.size() ? total[hp] : 0;
            // No module total for this haplotype (it does not cross the site, or no CN was typed for
            // it): there is nothing to split, and reading it as 0 would publish cn_a = cn_b = 0 as a call.
            if (tot < 0) { e.resolved = false; continue; }
            e.cn_a = static_cast<int>(std::llround(med * static_cast<double>(tot)));
            e.cn_b = static_cast<int>(tot) - e.cn_a;
            if (e.cn_b < 0) { e.cn_b = 0; e.cn_a = static_cast<int>(tot); }
        }
    }
    return out;
}

} // namespace

std::vector<std::vector<GeneCnEvidence>> resolve_gene_cn(
    const std::vector<std::string>& markers,
    const std::vector<std::string>& hap_seqs,
    const std::vector<long>& total,
    double near_identical_jaccard,
    std::size_t k,
    const std::vector<std::string>* backgrounds,
    const std::vector<char>* in_copies,
    const std::vector<std::string>* spans) {

    std::vector<std::vector<GeneCnEvidence>> ev = resolve_gene_cn_kmer(markers, hap_seqs, k, backgrounds);
    const std::size_t ng = markers.size();
    std::vector<char> pooled(ng, 1);   // genes keeping the pooled count (not split per site as a pair)
    // Copies read as runs of hits raise a count that pooled hits round down, for the genes counted pooled. A gene in
    // the module's copies raised past the module total is not clipped here: the caller reports it as unreliable.
    const auto raise_by_runs = [&]() {
        if (spans == nullptr) return;
        const std::vector<std::vector<GeneCnEvidence>> runs =
            resolve_gene_cn_kmer(markers, hap_seqs, k, backgrounds, spans, &pooled);
        for (std::size_t hp = 0; hp < ev.size(); ++hp)
            for (std::size_t g = 0; g < ng; ++g) {
                GeneCnEvidence& e = ev[hp][g];
                if (!pooled[g] || !e.separable) continue;
                e.copies = runs[hp][g].copies;
                if (e.copies > e.cn) e.cn = e.copies;
            }
    };
    if (ng < 2 || k == 0 || k > 31) { raise_by_runs(); return ev; }

    std::vector<std::unordered_set<std::uint64_t>> gk(ng);
    for (std::size_t g = 0; g < ng; ++g) collect_canonical_kmers(markers[g], k, gk[g]);

    auto jaccard = [&](std::size_t i, std::size_t j) -> double {
        const auto& small = gk[i].size() <= gk[j].size() ? gk[i] : gk[j];
        std::size_t inter = 0;
        for (std::uint64_t km : small) if ((gk[i].size() <= gk[j].size() ? gk[j] : gk[i]).count(km)) ++inter;
        const std::size_t uni = gk[i].size() + gk[j].size() - inter;
        return uni ? static_cast<double>(inter) / static_cast<double>(uni) : 0.0;
    };

    // Near-identical adjacency -> connected components; route only components of size exactly 2.
    std::vector<std::vector<int>> adj(ng);
    for (std::size_t i = 0; i < ng; ++i)
        for (std::size_t j = i + 1; j < ng; ++j)
            if (!gk[i].empty() && !gk[j].empty() && jaccard(i, j) > near_identical_jaccard) {
                adj[i].push_back(static_cast<int>(j));
                adj[j].push_back(static_cast<int>(i));
            }
    std::vector<int> comp(ng, -1);
    int nc = 0;
    for (std::size_t s = 0; s < ng; ++s) {
        if (comp[s] != -1) continue;
        std::vector<int> stack{static_cast<int>(s)};
        comp[s] = nc;
        while (!stack.empty()) {
            const int u = stack.back(); stack.pop_back();
            for (int v : adj[u]) if (comp[v] == -1) { comp[v] = nc; stack.push_back(v); }
        }
        ++nc;
    }
    std::vector<std::vector<int>> members(nc);
    for (std::size_t g = 0; g < ng; ++g) members[comp[g]].push_back(static_cast<int>(g));

    for (int c = 0; c < nc; ++c) {
        if (members[c].size() != 2) continue;
        const int ia = members[c][0], ib = members[c][1];
        // The total counts the module's copies only, so a pair with a gene beside them is not split by it.
        if (in_copies != nullptr && (!(*in_copies)[static_cast<std::size_t>(ia)] || !(*in_copies)[static_cast<std::size_t>(ib)]))
            continue;
        std::unordered_set<std::uint64_t> other;
        for (std::size_t g = 0; g < ng; ++g)
            if (static_cast<int>(g) != ia && static_cast<int>(g) != ib) other.insert(gk[g].begin(), gk[g].end());
        const std::vector<PairSiteEvidence> pev =
            resolve_pair_cn_persite(markers[ia], markers[ib], other, hap_seqs, total, k);
        if (pev.empty() || pev[0].sites == 0) continue;   // per-site failed -> keep pooled dosage
        // The split hands each module copy to exactly one gene of the pair. When both genes sit in
        // every copy, the pooled dosages add up to twice the total and splitting would halve each gene.
        // The pooled sum is robust to gene conversion (a converted site moves hits between the two
        // genes, not out of the pair), so it decides: the split is skipped when the pair's median
        // dosage per module copy, over the haplotypes with a total, is 1.75 or more. The cut sits above
        // the midpoint because a module total short by a copy pushes a one-gene-per-copy pair toward
        // 1.5, and that pair must keep its split.
        {
            std::vector<double> per_copy;
            for (std::size_t hp = 0; hp < hap_seqs.size(); ++hp) {
                const long tot = hp < total.size() ? total[hp] : -1;
                if (tot <= 0 || !ev[hp][ia].separable || !ev[hp][ib].separable) continue;
                per_copy.push_back((ev[hp][ia].dosage + ev[hp][ib].dosage) / static_cast<double>(tot));
                        }
            if (!per_copy.empty()) {
                std::nth_element(per_copy.begin(), per_copy.begin() + per_copy.size() / 2, per_copy.end());
                            if (per_copy[per_copy.size() / 2] >= 1.75) continue;   // both genes per copy -> pooled
            }
        }
        pooled[static_cast<std::size_t>(ia)] = pooled[static_cast<std::size_t>(ib)] = 0;
        for (std::size_t hp = 0; hp < hap_seqs.size(); ++hp) {
            GeneCnEvidence& ea = ev[hp][ia];
            // Separability is per haplotype, not per pair: `sites` and `priv_*` exist for every
            // haplotype once the pair has sites anywhere, so only `resolved` says whether this
            // haplotype had evidence.
            if (!pev[hp].resolved) {
                // Leave the pooled estimate in place and mark it unreliable; the caller then reports
                // the module total for this haplotype rather than a split it has no basis for.
                ea.separable = false;
                ev[hp][ib].separable = false;
                continue;
            }
            ea.cn = ea.rounded = pev[hp].cn_a; ea.hits = pev[hp].hits_a; ea.priv_kmers = pev[hp].priv_a;
            ea.dosage = ea.priv_kmers > 0 ? static_cast<double>(ea.hits) / static_cast<double>(ea.priv_kmers) : 0.0;
            ea.separable = true;
            GeneCnEvidence& eb = ev[hp][ib];
            eb.cn = eb.rounded = pev[hp].cn_b; eb.hits = pev[hp].hits_b; eb.priv_kmers = pev[hp].priv_b;
            eb.dosage = eb.priv_kmers > 0 ? static_cast<double>(eb.hits) / static_cast<double>(eb.priv_kmers) : 0.0;
            eb.separable = true;
        }
    }
    raise_by_runs();
    return ev;
}

std::vector<std::vector<int>> find_hybrid_copies(
    const std::vector<std::vector<std::string>>& gene_exons,
    const std::vector<std::string>& gene_spans,
    const std::vector<std::string>& hap_seqs,
    std::size_t min_exons,
    std::size_t k) {

    const std::size_t ng = gene_exons.size();
    std::vector<std::vector<int>> out(hap_seqs.size(), std::vector<int>(ng, -1));
    if (ng < 2 || k == 0 || k > 31) return out;

    struct Pair {
        int a, b;
        std::size_t gap;                      // hits further apart are two copies: 10 kb, or the longer gene
        std::vector<int> slot_of_b;           // b's exon -> a's exon (-1: none)
        std::vector<std::size_t> inform;      // a's exons where either gene has private k-mers
        std::vector<long> priv_a, priv_b;     // private k-mers per a-exon slot
        std::vector<std::size_t> slot_bp;     // the longer of the two genes' exons at each slot
    };
    std::vector<Pair> pairs;
    // k-mer -> (pair, gene side 0/1, slot) entries
    std::unordered_map<std::uint64_t, std::vector<std::array<int, 3>>> where;
    std::vector<std::unordered_set<std::uint64_t>> span_k(ng);
    for (std::size_t g = 0; g < ng; ++g) collect_canonical_kmers(gene_spans[g], k, span_k[g]);
    for (std::size_t a = 0; a < ng; ++a)
        for (std::size_t b = a + 1; b < ng; ++b) {
            if (gene_exons[a].empty() || gene_exons[b].empty()) continue;
            std::vector<std::unordered_set<std::uint64_t>> a15(gene_exons[a].size());
            for (std::size_t e = 0; e < gene_exons[a].size(); ++e) collect_canonical_kmers(gene_exons[a][e], 15, a15[e]);
            Pair p{static_cast<int>(a), static_cast<int>(b),
                   std::max<std::size_t>({10000, gene_spans[a].size(), gene_spans[b].size()}),
                   std::vector<int>(gene_exons[b].size(), -1), {}, {}, {}, {}};
            std::size_t mapped = 0;
            for (std::size_t e = 0; e < gene_exons[b].size(); ++e) {
                std::unordered_set<std::uint64_t> b15;
                collect_canonical_kmers(gene_exons[b][e], 15, b15);
                std::size_t best = 0, best_n = 0;
                for (std::size_t x = 0; x < a15.size(); ++x) {
                    std::size_t n = 0;
                    for (std::uint64_t km : b15) n += a15[x].count(km);
                    if (n > best_n) { best_n = n; best = x; }
                }
                if (best_n >= 5) { p.slot_of_b[e] = static_cast<int>(best); ++mapped; }
            }
            if (2 * mapped < gene_exons[b].size()) continue;   // not paralogs exon for exon
            p.slot_bp.assign(gene_exons[a].size(), 0);
            for (std::size_t e = 0; e < gene_exons[a].size(); ++e) p.slot_bp[e] = gene_exons[a][e].size();
            for (std::size_t e = 0; e < gene_exons[b].size(); ++e)
                if (p.slot_of_b[e] >= 0) {
                    std::size_t& sb = p.slot_bp[static_cast<std::size_t>(p.slot_of_b[e])];
                    sb = std::max(sb, gene_exons[b][e].size());
                }
            const int pi = static_cast<int>(pairs.size());
            p.priv_a.assign(gene_exons[a].size(), 0);
            p.priv_b.assign(gene_exons[a].size(), 0);
            std::unordered_set<std::uint64_t> taken;
            for (int side = 0; side < 2; ++side) {
                const std::size_t g = side == 0 ? a : b;
                // private: in no other gene of the target, so a third paralog's copy does not read as this one
                const auto elsewhere = [&](std::uint64_t km) {
                    for (std::size_t x = 0; x < ng; ++x) if (x != g && span_k[x].count(km)) return true;
                    return false;
                };
                for (std::size_t e = 0; e < gene_exons[g].size(); ++e) {
                    const int slot = side == 0 ? static_cast<int>(e) : p.slot_of_b[e];
                    if (slot < 0) continue;
                    std::unordered_set<std::uint64_t> ks;
                    collect_canonical_kmers(gene_exons[g][e], k, ks);
                    for (std::uint64_t km : ks) {
                        if (elsewhere(km) || !taken.insert(km).second) continue;
                        where[km].push_back({pi, side, slot});
                        ++(side == 0 ? p.priv_a : p.priv_b)[static_cast<std::size_t>(slot)];
                    }
                }
            }
            for (std::size_t s = 0; s < p.priv_a.size(); ++s)
                if (p.priv_a[s] >= 5 || p.priv_b[s] >= 5) p.inform.push_back(s);
            pairs.push_back(std::move(p));
        }
    if (pairs.empty()) return out;

    // Per haplotype, per pair: each copy's word over the informative slots ('A', 'B' or '.'), and where the copy
    // lies on the haplotype (first and last hit).
    std::vector<std::vector<std::vector<std::string>>> words(hap_seqs.size(), std::vector<std::vector<std::string>>(pairs.size()));
    std::vector<std::vector<std::vector<std::pair<std::size_t, std::size_t>>>> spans(
        hap_seqs.size(), std::vector<std::vector<std::pair<std::size_t, std::size_t>>>(pairs.size()));
    const std::uint64_t mask = (1ULL << (2 * k)) - 1;
    const unsigned top = static_cast<unsigned>(2 * (k - 1));
    for (std::size_t h = 0; h < hap_seqs.size(); ++h) {
        std::vector<std::vector<std::pair<std::size_t, std::array<int, 3>>>> hits(pairs.size());
        std::vector<std::vector<std::uint64_t>> hit_km(pairs.size());
        std::uint64_t fwd = 0, rev = 0; std::size_t len = 0, pos = 0;
        for (char c : hap_seqs[h]) {
            ++pos;
            const int b = base_code(c);
            if (b < 0) { len = 0; fwd = 0; rev = 0; continue; }
            fwd = ((fwd << 2) | static_cast<std::uint64_t>(b)) & mask;
            rev = (rev >> 2) | (static_cast<std::uint64_t>(3 - b) << top);
            if (++len < k) continue;
            const std::uint64_t km = std::min(fwd, rev);
            const auto it = where.find(km);
            if (it == where.end()) continue;
            for (const auto& w : it->second) {
                hits[static_cast<std::size_t>(w[0])].emplace_back(pos, w);
                hit_km[static_cast<std::size_t>(w[0])].push_back(km);
            }
        }
        for (std::size_t pi = 0; pi < pairs.size(); ++pi) {
            const auto& hs = hits[pi];
            if (hs.empty()) continue;
            // Copies, each read on its own strand: the first change of exon slot inside a copy says which way it
            // runs. A copy ends at a long gap; where an exon it already holds comes again further on than the
            // exon is long (the next copy, as when two copies face each other or a partial copy follows); where
            // the slots run back more than two against its strand (the next copy in tandem); or where they run
            // back and keep going back (the next copy, on the other strand). Any other small step back is a
            // stray hit.
            std::vector<std::vector<std::size_t>> copies;
            int dir = 0, run_max = 0, cur = -1;
            std::size_t cur_from = 0;   // where the current slot's hits began
            std::unordered_map<int, std::size_t> slot_from;   // where each slot's hits began in this copy
            const Pair& pp = pairs[pi];
            for (std::size_t x = 0; x < hs.size(); ++x) {
                const int slot = hs[x].second[2];
                const std::size_t at = hs[x].first;
                bool same = !copies.empty() && at - hs[copies.back().back()].first < pp.gap;
                if (same && slot == cur) {
                    same = at - cur_from <= pp.slot_bp[static_cast<std::size_t>(slot)];
                } else if (same && dir == 0) {
                    dir = slot > cur ? 1 : -1;
                    run_max = std::max(cur * dir, slot * dir);
                } else if (same && (slot - cur) * dir > 0) {
                    run_max = std::max(run_max, slot * dir);
                } else if (same) {
                    const auto sf = slot_from.find(slot);
                    if (slot * dir <= run_max - 3) same = false;
                    else if (sf != slot_from.end() && at - sf->second > pp.slot_bp[static_cast<std::size_t>(slot)]) same = false;
                    else
                        for (std::size_t y = x + 1; y < hs.size() && hs[y].first - at < pp.gap; ++y)
                            if (hs[y].second[2] != slot) { same = (hs[y].second[2] - slot) * dir > 0; break; }
                }
                if (same) copies.back().push_back(x);
                else { copies.push_back({x}); dir = 0; slot_from.clear(); }
                if (!same || slot != cur) { cur = slot; cur_from = at; }
                slot_from.emplace(slot, at);
            }
            const Pair& p = pairs[pi];
            for (const auto& cp : copies) {
                std::vector<std::unordered_set<std::uint64_t>> sa(p.priv_a.size()), sb(p.priv_a.size());
                for (std::size_t x : cp) {
                    const auto& w = hs[x].second;
                    (w[1] == 0 ? sa : sb)[static_cast<std::size_t>(w[2])].insert(hit_km[pi][x]);
                }
                std::string word;
                for (std::size_t s : p.inform) {
                    const double fa = p.priv_a[s] ? static_cast<double>(sa[s].size()) / static_cast<double>(p.priv_a[s]) : 0.0;
                    const double fb = p.priv_b[s] ? static_cast<double>(sb[s].size()) / static_cast<double>(p.priv_b[s]) : 0.0;
                    word += fa >= 0.5 && fb < 0.5 ? 'A' : fb >= 0.5 && fa < 0.5 ? 'B' : '.';
                }
                words[h][pi].push_back(std::move(word));
                spans[h][pi].emplace_back(hs[cp.front()].first, hs[cp.back()].first);
            }
        }
    }

    // 'A', 'B', 'H' (hybrid) or 0 (too few exons read) for a copy, over the slots kept
    const auto label = [](const std::string& w, const std::vector<char>& keep, bool hybrids = true) -> char {
        std::string called;
        for (std::size_t i = 0; i < w.size(); ++i) if (keep[i] && w[i] != '.') called += w[i];
        if (called.size() < 3) return 0;
        const std::size_t na = static_cast<std::size_t>(std::count(called.begin(), called.end(), 'A'));
        const char maj = 2 * na >= called.size() ? 'A' : 'B';
        std::size_t first = called.size(), last = 0, n = 0;
        for (std::size_t i = 0; i < called.size(); ++i)
            if (called[i] != maj) { first = std::min(first, i); last = i; ++n; }
        const bool end_block = n > 0 && last - first + 1 == n && (first == 0 || last + 1 == called.size());
        return hybrids && end_block ? 'H' : maj;
    };
    // Per haplotype and gene: the hybrid copies found by any pair holding the gene, as spans on the haplotype.
    std::vector<std::vector<std::vector<std::pair<std::size_t, std::size_t>>>> hybrid_spans(
        hap_seqs.size(), std::vector<std::vector<std::pair<std::size_t, std::size_t>>>(ng));
    for (std::size_t pi = 0; pi < pairs.size(); ++pi) {
        const std::size_t ns = pairs[pi].inform.size();
        const std::vector<char> all(ns, 1);
        std::vector<char> keep(ns, 1);
        for (int pass = 0; pass < 2; ++pass) {
            std::vector<std::array<long, 2>> own(ns, {0, 0}), tot(ns, {0, 0});
            for (const auto& per_hap : words)
                for (const std::string& w : per_hap[pi]) {
                    // The second pass leaves out the copies that read as hybrids on every slot, so a hybrid
                    // common in the panel does not take its crossover slots out of the reading.
                    if (pass == 1 && label(w, all) == 'H') continue;
                    // by its majority: a hybrid's own gene still says what that gene's copies read
                    const char l = label(w, keep, false);
                    if (l != 'A' && l != 'B') continue;
                    const int side = l == 'A' ? 0 : 1;
                    // a slot a copy does not read (a partial copy) says nothing about the slot
                    for (std::size_t i = 0; i < ns; ++i)
                        if (w[i] != '.') { ++tot[i][side]; if (w[i] == l) ++own[i][side]; }
                }
            for (std::size_t i = 0; i < ns; ++i)
                keep[i] = tot[i][0] > 0 && tot[i][1] > 0 && 5 * own[i][0] >= 4 * tot[i][0] && 5 * own[i][1] >= 4 * tot[i][1];
        }
        if (static_cast<std::size_t>(std::count(keep.begin(), keep.end(), 1)) < min_exons) continue;
        const std::size_t a = static_cast<std::size_t>(pairs[pi].a), b = static_cast<std::size_t>(pairs[pi].b);
        for (std::size_t h = 0; h < hap_seqs.size(); ++h) {
            out[h][a] = std::max(out[h][a], 0);
            out[h][b] = std::max(out[h][b], 0);
            for (std::size_t c = 0; c < words[h][pi].size(); ++c)
                if (label(words[h][pi][c], keep) == 'H') {
                    hybrid_spans[h][a].push_back(spans[h][pi][c]);
                    hybrid_spans[h][b].push_back(spans[h][pi][c]);
                }
        }
    }
    // One hybrid copy can show in two pairs (a third paralog's pair reads it too): overlapping spans are
    // one copy, separate spans separate copies.
    for (std::size_t h = 0; h < hap_seqs.size(); ++h)
        for (std::size_t g = 0; g < ng; ++g) {
            auto& sp = hybrid_spans[h][g];
            if (sp.empty()) continue;
            std::sort(sp.begin(), sp.end());
            int n = 0;
            std::size_t end = 0;
            for (const auto& [from, to] : sp) {
                if (n == 0 || from > end) { ++n; end = to; }
                else end = std::max(end, to);
            }
            out[h][g] = n;
        }
    return out;
}

std::vector<int> paralog_groups(const std::vector<std::string>& seqs, std::size_t k, double min_shared) {
    const std::size_t n = seqs.size();
    std::vector<std::unordered_set<std::uint64_t>> ks(n);
    for (std::size_t g = 0; g < n; ++g) collect_canonical_kmers(seqs[g], k, ks[g]);
    std::vector<int> group(n);
    for (std::size_t g = 0; g < n; ++g) group[g] = static_cast<int>(g);
    const std::function<int(int)> root = [&](int g) { return group[g] == g ? g : group[g] = root(group[g]); };
    for (std::size_t a = 0; a < n; ++a)
        for (std::size_t b = a + 1; b < n; ++b) {
            const auto& small = ks[a].size() <= ks[b].size() ? ks[a] : ks[b];
            const auto& large = ks[a].size() <= ks[b].size() ? ks[b] : ks[a];
            if (small.empty()) continue;
            std::size_t shared = 0;
            for (std::uint64_t km : small) shared += large.count(km);
            if (static_cast<double>(shared) >= min_shared * static_cast<double>(small.size()))
                group[root(static_cast<int>(a))] = root(static_cast<int>(b));
        }
    for (std::size_t g = 0; g < n; ++g) group[g] = root(static_cast<int>(g));
    return group;
}

} // namespace panvar
