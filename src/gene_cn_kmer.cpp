#include "panvar/gene_cn_kmer.hpp"

#include "panvar/minimap2_align.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
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
    const std::uint64_t mask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
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

// As above, but keeping how MANY times each canonical k-mer occurs. Copy number is an occurrence
// ratio, so a set loses exactly the multiplicity the estimate is made of.
void count_canonical_kmers(const std::string& seq, std::size_t k,
                           std::unordered_map<std::uint64_t, long>& out) {
    if (k == 0 || k > 32 || seq.size() < k) return;
    const std::uint64_t mask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
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

} // namespace

std::vector<std::vector<GeneCnEvidence>> resolve_gene_cn_kmer(
    const std::vector<std::string>& gene_markers,
    const std::vector<std::string>& hap_seqs,
    std::size_t k) {

    const std::size_t n_genes = gene_markers.size();
    std::vector<std::vector<GeneCnEvidence>> out(
        hap_seqs.size(), std::vector<GeneCnEvidence>(n_genes));
    if (n_genes == 0 || k == 0 || k > 31) return out;

    // 1) canonical k-mer COUNTS per gene. Counts, not a set: a k-mer occurring twice inside one copy
    //    of a gene contributes two hits from that one copy, so the denominator below has to know.
    std::vector<std::unordered_map<std::uint64_t, long>> gene_counts(n_genes);
    for (std::size_t g = 0; g < n_genes; ++g)
        count_canonical_kmers(gene_markers[g], k, gene_counts[g]);

    // 2) private k-mers: unique to one gene against its siblings. Map private k-mer -> owning gene.
    //
    // Screening these against the WHOLE locus reference as well -- dropping any k-mer that also occurs
    // outside the gene's own marker -- was built and measured, and it is WORSE. Against pangene truth
    // it cost CYP2D6 3 of 127 exact calls (126 -> 123) and, once CDS-junction k-mers were kept rather
    // than dropped as absent, GSTM1 2 of 466 as well. It removes about a third of the markers, and the
    // dosage ratio loses more to the smaller denominator than it gains in specificity. Not retained,
    // switched off or otherwise; see docs/algorithms/call.md.
    std::unordered_map<std::uint64_t, int> owner;
    std::vector<long> priv_size(n_genes, 0);
    for (std::size_t g = 0; g < n_genes; ++g) {
        for (const auto& [km, in_gene] : gene_counts[g]) {
            bool shared = false;
            for (std::size_t o = 0; o < n_genes && !shared; ++o)
                if (o != g && gene_counts[o].count(km)) shared = true;
            if (shared) continue;
            owner[km] = static_cast<int>(g);
            // The denominator is OCCURRENCES per reference copy, so that hits/denominator is the
            // number of copies. Counting distinct k-mers instead made a k-mer repeated within one
            // copy read as extra copies of the gene.
            priv_size[g] += in_gene;
        }
    }

    // 3) per haplotype: count private-k-mer occurrences, derive dosage + CN.
    const std::uint64_t mask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
    const unsigned top_shift = static_cast<unsigned>(2 * (k - 1));
    for (std::size_t h = 0; h < hap_seqs.size(); ++h) {
        std::vector<long> hits(n_genes, 0);
        std::uint64_t fwd = 0, rev = 0;
        std::size_t len = 0;
        for (char c : hap_seqs[h]) {
            const int b = base_code(c);
            if (b < 0) { len = 0; fwd = 0; rev = 0; continue; }
            fwd = ((fwd << 2) | static_cast<std::uint64_t>(b)) & mask;
            rev = (rev >> 2) | (static_cast<std::uint64_t>(3 - b) << top_shift);
            if (++len >= k) {
                const auto it = owner.find(std::min(fwd, rev));
                if (it != owner.end()) ++hits[it->second];
            }
        }
        for (std::size_t g = 0; g < n_genes; ++g) {
            GeneCnEvidence& e = out[h][g];
            e.priv_kmers = priv_size[g];
            e.hits = hits[g];
            e.separable = priv_size[g] > 0;
            e.dosage = e.separable ? static_cast<double>(hits[g]) / static_cast<double>(priv_size[g]) : 0.0;
            e.cn = e.separable ? static_cast<int>(std::llround(e.dosage)) : 0;
        }
    }
    return out;
}

namespace {

// Canonical 2-bit code of the length-k window at seq[i..i+k) (assumes ACGT; caller checks). Uses the same
// canonicalisation as collect_canonical_kmers so codes are directly comparable across the two.
inline std::uint64_t canon_code(const std::string& seq, std::size_t i, std::size_t k) {
    const std::uint64_t mask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
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
    // Did THIS haplotype produce a usable split? `sites` and `priv_*` describe the marker sets, which
    // exist for every haplotype once the pair has divergent sites at all, so neither can answer it.
    // Without this the pair's global separability was read as per-haplotype evidence and a haplotype
    // with no usable site published cn_a=cn_b=0 as a confident call.
    long usable_sites = 0;
    bool resolved = false;
};

// Per-site consensus split of a near-identical gene pair (a, b). Aligns b (query) to a (target) ONCE to
// locate divergent sites, builds each site's A/B allele k-mers (private to the pair vs `other`), then per
// haplotype takes the median across sites of the site-local A fraction and splits total[hap] by it. Returns
// sites==0 when the pair could not be resolved (no alignment / no usable sites) so the caller keeps pooled.
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

    struct Site { std::vector<std::uint64_t> a_km, b_km; };
    std::vector<Site> sites;
    std::size_t tpos = h.target_start_bp, qpos = h.query_start_bp;   // a (target), b (query)
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
        collect_site_windows(marker_b, k, b_lo, b_hi, bw);
        Site s;
        for (std::uint64_t km : aw) if (!bset.count(km) && !other.count(km)) s.a_km.push_back(km);
        for (std::uint64_t km : bw) if (!aset.count(km) && !other.count(km)) s.b_km.push_back(km);
        if (!s.a_km.empty() && !s.b_km.empty()) sites.push_back(std::move(s));
    }
    if (sites.empty()) return out;

    std::unordered_map<std::uint64_t, std::pair<int, int>> loc;   // k-mer -> (site, 0=a / 1=b)
    for (std::size_t si = 0; si < sites.size(); ++si) {
        for (std::uint64_t km : sites[si].a_km) loc[km] = {static_cast<int>(si), 0};
        for (std::uint64_t km : sites[si].b_km) loc[km] = {static_cast<int>(si), 1};
    }

    const std::uint64_t mask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);
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
                if (it != loc.end()) { if (it->second.second == 0) ++ah[it->second.first]; else ++bh[it->second.first]; }
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
    std::size_t k) {

    std::vector<std::vector<GeneCnEvidence>> ev = resolve_gene_cn_kmer(markers, hap_seqs, k);
    const std::size_t ng = markers.size();
    if (ng < 2 || k == 0 || k > 31) return ev;

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
        std::unordered_set<std::uint64_t> other;
        for (std::size_t g = 0; g < ng; ++g)
            if (static_cast<int>(g) != ia && static_cast<int>(g) != ib) other.insert(gk[g].begin(), gk[g].end());
        const std::vector<PairSiteEvidence> pev =
            resolve_pair_cn_persite(markers[ia], markers[ib], other, hap_seqs, total, k);
        if (pev.empty() || pev[0].sites == 0) continue;   // per-site failed -> keep pooled dosage
        for (std::size_t hp = 0; hp < hap_seqs.size(); ++hp) {
            GeneCnEvidence& ea = ev[hp][ia];
            // Separability is per HAPLOTYPE, not per pair. `sites` says the pair could be separated
            // somewhere in the panel and `priv_*` are marker-set sizes that exist for every haplotype
            // regardless, so neither answers whether THIS haplotype had evidence. Only `resolved`
            // does. Testing priv_kmers instead was true whenever the pair had sites at all, so a
            // haplotype with no usable site still published cn=0 as a confident call -- and the check
            // could never fire, which is what made it look like the case did not arise.
            if (!pev[hp].resolved) {
                // Leave the pooled estimate in place and mark it unreliable; the caller then reports
                // the module total for this haplotype rather than a split it has no basis for.
                ea.separable = false;
                ev[hp][ib].separable = false;
                continue;
            }
            ea.cn = pev[hp].cn_a; ea.hits = pev[hp].hits_a; ea.priv_kmers = pev[hp].priv_a;
            ea.dosage = ea.priv_kmers > 0 ? static_cast<double>(ea.hits) / static_cast<double>(ea.priv_kmers) : 0.0;
            ea.separable = true;
            GeneCnEvidence& eb = ev[hp][ib];
            eb.cn = pev[hp].cn_b; eb.hits = pev[hp].hits_b; eb.priv_kmers = pev[hp].priv_b;
            eb.dosage = eb.priv_kmers > 0 ? static_cast<double>(eb.hits) / static_cast<double>(eb.priv_kmers) : 0.0;
            eb.separable = true;
        }
    }
    return ev;
}

} // namespace panvar
