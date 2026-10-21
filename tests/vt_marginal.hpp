// M13-A (DIAGNOSTIC): variation-tolerant structural scoring -- the truncated, group-factorised marginal.
//
// Model for one candidate arrangement G (two homologues, each = fixed context + a target block):
//   S      the block content of each homologue: its spelled block with substitutions at ADMITTED sites only; every other
//          position is fixed at the spelled base. One state per site, shared by every read (never redrawn per read).
//   prior  the C1 whole-component mixture (M6-F) over panel block alleles, a function of the sequence alone and identical
//          for every candidate:  P(X) = sum_k w_k prod_{x in block} e(X_x | b_k(x)),  e = 1-eps on match, eps/3 otherwise.
//          Substitution-only: only components of the block's length contribute (declared limitation). Homologues are
//          independent given G.
//   reads  the existing paired-fragment likelihood of the complete diplotype (score_diplotypes: placements summed,
//          sequencing error, exposure, outlier/background) -- the residual rate eps is NOT the sequencing error.
// Marginal (truncated to the admitted set; a lower bound on the full sequence marginal):
//   log m(G) = log sum_{k1,k2} w_k1 w_k2 sum_{X in S} prod_x e(X|b_k) L(R | X)
// Exactness argument: given (k1,k2) the prior factorises over positions, and a fragment's term depends only on sites its
// placements touch, so with sites grouped by shared fragments, L(X) = L(X0) prod_g exp(h_g(X_g)) and the sum over X is a
// product of per-group sums, each enumerated exactly. A fragment that gains or loses a placement only under a variant state
// can couple groups invisibly; the joint MAP-state rescoring (interaction_error) measures that.
#pragma once

#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_model.hpp"
#include "panvar/mosaic_spelling.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <vector>

namespace vt {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();

struct Settings {
    double eps_res = 1e-4;                 // M6-F's C1 residual divergence (its grid floor; a known limitation)
    std::size_t min_support = 2;           // fragments showing the same non-reference base at a block position
    std::size_t max_sites_per_homologue = 300;
    std::size_t max_group_sites = 12;      // 4,096 states
    std::size_t max_variant_haplotypes = 300;   // per candidate
    bool admit = true;                     // false: no admitted site (fixed content)
    bool dry = false;                      // admission and grouping only (preflight); no state scoring
    std::vector<int> forced_hom;           // gates only: when non-empty, these sites replace admission
    std::vector<std::size_t> forced_pos;
    std::vector<char> forced_alt;
};

struct Component {
    std::string block;   // a panel block allele
    double log_w = 0.0;  // log template weight
};

struct Candidate {
    std::string label;
    std::array<std::string, 2> hap;                     // complete homologues (spelled)
    std::array<panvar::SequenceInterval, 2> block;     // target block interval on each
    std::array<std::vector<Component>, 2> comps;       // panel components of the block's length, per homologue
};

struct Site {
    int hom = 0;            // homologue index 0/1
    std::size_t pos = 0;    // position on that homologue
    char ref = 'N', alt = 'N';
    std::size_t support = 0;
};

struct Result {
    std::string label;
    std::string status = "OK";   // OK | UNRESOLVED_<reason>
    double base_fragment = 0, base_exposure = 0, base_loglik = 0;
    double log_prior_panel = kNegInf;   // log P_C1(block1) + log P_C1(block2) at fixed content
    double log_m_panel = kNegInf;       // base loglik + panel prior (fixed content under the C1 prior)
    double log_m_res = kNegInf;         // truncated variation-tolerant marginal
    double eq_fit_gain = 0;             // E_q[log L] - log L(X0)
    double kl = 0;                      // KL(q || prior restricted), from log m = E_q[log L] - KL
    double expected_alt = 0;            // posterior expected substituted sites
    std::size_t sites = 0, groups = 0, max_group = 0, variant_haplotypes = 0, map_alt = 0;
    double interaction_error = 0;       // [log L(X*) - log L(X0)] - sum_g h_g(X*_g), at the joint MAP state
    std::vector<Site> site_list;
};

inline double lse(const std::vector<double>& v) {
    double m = kNegInf;
    for (double x : v) m = std::max(m, x);
    if (!std::isfinite(m)) return m;
    long double s = 0;
    for (double x : v) s += std::exp(static_cast<long double>(x - m));
    return m + static_cast<double>(std::log(s));
}

inline std::string revcomp(const std::string& x) {
    std::string r(x.rbegin(), x.rend());
    for (char& c : r) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
    return r;
}

// log e(x | b) for one position under the residual kernel
inline double log_e(char x, char b, double eps) { return x == b ? std::log1p(-eps) : std::log(eps / 3.0); }

// Sum over block positions of log e(X_x | b_x), skipping positions in `skip` (block coordinates).
inline double log_kernel(const std::string& x, const std::string& b, double eps, const std::set<std::size_t>& skip) {
    const double m = std::log1p(-eps), d = std::log(eps / 3.0);
    std::size_t same = 0, diff = 0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (skip.count(i)) continue;
        (x[i] == b[i] ? same : diff)++;
    }
    return static_cast<double>(same) * m + static_cast<double>(diff) * d;
}

inline double log_c1(const std::string& block, const std::vector<Component>& comps, double eps) {
    std::vector<double> t;
    for (const Component& c : comps) if (c.block.size() == block.size()) t.push_back(c.log_w + log_kernel(block, c.block, eps, {}));
    return lse(t);
}

// Site admission from the recorded origins on each base homologue (truth-free; the same rule for every candidate).
// origins[hom][f] = the fragment's accepted origins on that homologue's base sequence.
inline std::vector<Site> admit_sites(const Candidate& c, const std::vector<panvar::ReadPair>& fragments,
                                     const std::array<std::vector<std::vector<panvar::FragmentPlacement>>, 2>& origins,
                                     const Settings& s) {
    std::vector<Site> out;
    if (!s.forced_pos.empty()) {
        for (std::size_t i = 0; i < s.forced_pos.size(); ++i)
            out.push_back(Site{s.forced_hom[i], s.forced_pos[i], c.hap[s.forced_hom[i]][s.forced_pos[i]], s.forced_alt[i], 0});
        return out;
    }
    if (!s.admit) return out;
    for (int h = 0; h < 2; ++h) {
        const std::string& hap = c.hap[h];
        std::map<std::pair<std::size_t, char>, std::set<std::size_t>> ev;
        for (std::size_t f = 0; f < fragments.size(); ++f) {
            const auto& O = origins[h][f];
            if (O.empty()) continue;
            std::size_t best = std::numeric_limits<std::size_t>::max();
            for (const auto& o : O) best = std::min(best, o.forward_mismatches + o.reverse_mismatches);
            for (const auto& o : O) {
                if (o.forward_mismatches + o.reverse_mismatches != best) continue;
                const bool ff = o.orientation == panvar::MateOrientation::FirstForward;
                const std::string fwd = ff ? fragments[f].first.bases : fragments[f].second.bases;
                const std::string rev = revcomp(ff ? fragments[f].second.bases : fragments[f].first.bases);
                const auto scan = [&](const std::string& r, std::size_t start, std::size_t consumed, std::size_t mism) {
                    if (consumed != r.size() || start + r.size() > hap.size()) return;   // gapped mates are not admitted from
                    std::vector<std::size_t> pos;
                    for (std::size_t i = 0; i < r.size(); ++i) if (r[i] != hap[start + i]) pos.push_back(start + i);
                    if (pos.size() != mism) return;   // positional comparison must reproduce the scorer's mismatch count
                    for (std::size_t p : pos) {
                        if (p >= c.block[h].begin && p < c.block[h].end) ev[{p, r[p - start]}].insert(f);
                    }
                };
                scan(fwd, o.forward_start, o.forward_consumed_bp, o.forward_mismatches);
                scan(rev, o.reverse_start, o.reverse_consumed_bp, o.reverse_mismatches);
            }
        }
        std::map<std::size_t, std::pair<std::size_t, char>> best_alt;   // pos -> (support, alt)
        for (const auto& [key, fr] : ev) {
            if (fr.size() < s.min_support) continue;
            auto& b = best_alt[key.first];
            if (fr.size() > b.first || (fr.size() == b.first && key.second < b.second)) b = {fr.size(), key.second};
        }
        for (const auto& [p, sa] : best_alt) out.push_back(Site{h, p, hap[p], sa.second, sa.first});
    }
    return out;
}

// Groups of sites coupled by any fragment whose base-state origins touch both (mate spans only).
inline std::vector<std::vector<std::size_t>> group_sites(const std::vector<Site>& sites, std::size_t n_fragments,
                                                         const std::array<std::vector<std::vector<panvar::FragmentPlacement>>, 2>& origins) {
    std::vector<std::size_t> parent(sites.size());
    std::iota(parent.begin(), parent.end(), 0);
    const std::function<std::size_t(std::size_t)> find = [&](std::size_t x) { return parent[x] == x ? x : parent[x] = find(parent[x]); };
    std::array<std::map<std::size_t, std::size_t>, 2> at;   // pos -> site index per homologue
    for (std::size_t i = 0; i < sites.size(); ++i) at[sites[i].hom][sites[i].pos] = i;
    for (std::size_t f = 0; f < n_fragments; ++f) {
        std::vector<std::size_t> touched;
        for (int h = 0; h < 2; ++h) {
            for (const auto& o : origins[h][f]) {
                const std::array<std::pair<std::size_t, std::size_t>, 2> spans = {
                    std::make_pair(o.forward_start, o.forward_start + o.forward_consumed_bp),
                    std::make_pair(o.reverse_start, o.reverse_start + o.reverse_consumed_bp)};
                for (const auto& sp : spans) {
                    for (auto it = at[h].lower_bound(sp.first); it != at[h].end() && it->first < sp.second; ++it) touched.push_back(it->second);
                }
            }
        }
        for (std::size_t i = 1; i < touched.size(); ++i) parent[find(touched[i])] = find(touched[0]);
    }
    std::map<std::size_t, std::vector<std::size_t>> g;
    for (std::size_t i = 0; i < sites.size(); ++i) g[find(i)].push_back(i);
    std::vector<std::vector<std::size_t>> out;
    for (auto& [r, v] : g) out.push_back(v);
    return out;
}

inline std::string apply(const std::string& base, const std::vector<Site>& sites, const std::vector<std::size_t>& idx, unsigned mask, int hom) {
    std::string s = base;
    for (std::size_t j = 0; j < idx.size(); ++j) {
        const Site& t = sites[idx[j]];
        if (t.hom == hom && ((mask >> j) & 1U)) s[t.pos] = t.alt;
    }
    return s;
}

// Score a list of diplotypes given as byte pairs; returns log likelihood (fragment + exposure) per pair.
inline std::vector<std::pair<double, double>> score_pairs(const std::vector<std::pair<std::string, std::string>>& pairs,
                                                          const std::vector<panvar::ReadPair>& fragments, const panvar::LibraryModel& lib,
                                                          const panvar::FragmentScoringOptions& opt) {
    std::vector<panvar::MosaicPath> paths;
    std::vector<std::string> seqs;
    std::map<std::string, std::size_t> slot;
    const auto slot_of = [&](const std::string& s) {
        auto it = slot.find(s);
        if (it != slot.end()) return it->second;
        const std::size_t k = seqs.size();
        slot.emplace(s, k);
        paths.push_back(panvar::MosaicPath{{std::numeric_limits<std::size_t>::max(), k}});
        seqs.push_back(s);
        return k;
    };
    std::vector<std::pair<std::size_t, std::size_t>> ps;
    for (const auto& p : pairs) ps.push_back({slot_of(p.first), slot_of(p.second)});
    const panvar::InternedHaplotypes in = panvar::intern_haplotypes(paths, seqs);
    std::vector<panvar::HaplotypeCandidate> cands;
    for (const auto& p : ps) {
        std::size_t x = in.index_of[p.first], y = in.index_of[p.second];
        if (y < x) std::swap(x, y);
        cands.push_back({x, y});
    }
    const panvar::FragmentScoringResult r = panvar::score_diplotypes(in, cands, fragments, lib, opt);
    std::vector<std::pair<double, double>> out;
    for (const auto& sc : r.scores) out.push_back({sc.fragment_term, sc.exposure_term});
    return out;
}

// Base scoring with origin recording for a set of candidates (serial placement, as the scorer requires for recording).
struct BaseScan {
    std::map<std::string, std::vector<std::vector<panvar::FragmentPlacement>>> origins_by_bytes;   // haplotype bytes -> per fragment
};
inline BaseScan base_scan(const std::vector<Candidate>& cs, const std::vector<panvar::ReadPair>& fragments, const panvar::LibraryModel& lib,
                          panvar::FragmentScoringOptions opt) {
    std::vector<panvar::MosaicPath> paths;
    std::vector<std::string> seqs;
    std::map<std::string, std::size_t> slot;
    for (const Candidate& c : cs) for (const auto& h : c.hap) if (!slot.count(h)) {
        slot.emplace(h, seqs.size());
        paths.push_back(panvar::MosaicPath{{std::numeric_limits<std::size_t>::max(), seqs.size()}});
        seqs.push_back(h);
    }
    const panvar::InternedHaplotypes in = panvar::intern_haplotypes(paths, seqs);
    opt.record_origins_for_fragments.clear();
    for (std::size_t f = 0; f < fragments.size(); ++f) opt.record_origins_for_fragments.push_back(f);
    std::vector<panvar::HaplotypeCandidate> cands;
    for (std::size_t u = 0; u < in.sequences.size(); ++u) cands.push_back({u, u});
    const panvar::FragmentScoringResult r = panvar::score_diplotypes(in, cands, fragments, lib, opt);
    BaseScan b;
    for (const auto& [bytes, k] : slot) {
        const std::size_t u = in.index_of[k];
        auto& v = b.origins_by_bytes[bytes];
        v.assign(fragments.size(), {});
        for (std::size_t f = 0; f < fragments.size(); ++f) {
            const auto it = r.recorded_origins.find({f, u});
            if (it != r.recorded_origins.end()) v[f] = it->second;
        }
    }
    return b;
}

// Evaluate every candidate of one contrast. All state diplotypes are scored in one call.
inline std::vector<Result> evaluate(const std::vector<Candidate>& cs, const std::vector<panvar::ReadPair>& fragments, const panvar::LibraryModel& lib,
                                    const panvar::FragmentScoringOptions& opt, const Settings& s, const BaseScan& scan) {
    struct Work {
        std::vector<Site> sites;
        std::vector<std::vector<std::size_t>> groups;
        std::vector<std::vector<std::size_t>> pair_index;   // per group, per state mask: index into `pairs`
        std::size_t base_index = 0;
    };
    std::vector<Result> res(cs.size());
    std::vector<Work> work(cs.size());
    std::vector<std::pair<std::string, std::string>> pairs;
    for (std::size_t ci = 0; ci < cs.size(); ++ci) {
        const Candidate& c = cs[ci];
        Result& R = res[ci];
        R.label = c.label;
        Work& W = work[ci];
        W.base_index = pairs.size();
        pairs.push_back({c.hap[0], c.hap[1]});
        const std::array<std::vector<std::vector<panvar::FragmentPlacement>>, 2> o = {scan.origins_by_bytes.at(c.hap[0]), scan.origins_by_bytes.at(c.hap[1])};
        W.sites = admit_sites(c, fragments, o, s);
        std::size_t per[2] = {0, 0};
        for (const Site& t : W.sites) ++per[t.hom];
        R.sites = W.sites.size();
        R.site_list = W.sites;
        if (per[0] > s.max_sites_per_homologue || per[1] > s.max_sites_per_homologue) { R.status = "UNRESOLVED_TOO_MANY_SITES"; continue; }
        W.groups = group_sites(W.sites, fragments.size(), o);
        R.groups = W.groups.size();
        std::size_t variants = 0;
        for (const auto& g : W.groups) {
            R.max_group = std::max(R.max_group, g.size());
            variants += (std::size_t{1} << g.size()) - 1;
        }
        R.variant_haplotypes = variants;
        if (R.max_group > s.max_group_sites) { R.status = "UNRESOLVED_GROUP_TOO_LARGE"; continue; }
        if (variants > s.max_variant_haplotypes) { R.status = "UNRESOLVED_TOO_MANY_STATES"; continue; }
        if (s.dry) continue;
        for (const auto& g : W.groups) {
            std::vector<std::size_t> idx(std::size_t{1} << g.size());
            idx[0] = W.base_index;
            for (unsigned m = 1; m < (1U << g.size()); ++m) {
                idx[m] = pairs.size();
                pairs.push_back({apply(c.hap[0], W.sites, g, m, 0), apply(c.hap[1], W.sites, g, m, 1)});
            }
            W.pair_index.push_back(idx);
        }
    }
    if (s.dry) {
        for (Result& R : res) if (R.status == "OK") R.status = "DRY";
        return res;
    }
    const auto scored = score_pairs(pairs, fragments, lib, opt);
    // joint MAP check pairs, appended after the marginal is known
    std::vector<std::pair<std::string, std::string>> joint;
    std::vector<std::size_t> joint_of(cs.size(), std::numeric_limits<std::size_t>::max());
    std::vector<std::vector<unsigned>> map_masks(cs.size());
    for (std::size_t ci = 0; ci < cs.size(); ++ci) {
        const Candidate& c = cs[ci];
        Result& R = res[ci];
        const Work& W = work[ci];
        R.base_fragment = scored[W.base_index].first;
        R.base_exposure = scored[W.base_index].second;
        R.base_loglik = R.base_fragment + R.base_exposure;
        const std::string b1 = c.hap[0].substr(c.block[0].begin, c.block[0].size()), b2 = c.hap[1].substr(c.block[1].begin, c.block[1].size());
        R.log_prior_panel = log_c1(b1, c.comps[0], s.eps_res) + log_c1(b2, c.comps[1], s.eps_res);
        R.log_m_panel = R.base_loglik + R.log_prior_panel;
        if (R.status != "OK") continue;
        // block-coordinate site sets per homologue
        std::array<std::set<std::size_t>, 2> skip;
        for (const Site& t : W.sites) skip[t.hom].insert(t.pos - c.block[t.hom].begin);
        std::vector<double> terms;
        struct KTerm { double total; std::vector<std::vector<double>> post; };   // per group posterior over masks
        std::vector<KTerm> kt;
        std::vector<std::size_t> k1s, k2s;
        for (std::size_t k1 = 0; k1 < c.comps[0].size(); ++k1) {
            if (c.comps[0][k1].block.size() != b1.size()) continue;
            for (std::size_t k2 = 0; k2 < c.comps[1].size(); ++k2) {
                if (c.comps[1][k2].block.size() != b2.size()) continue;
                const std::array<const std::string*, 2> comp = {&c.comps[0][k1].block, &c.comps[1][k2].block};
                double t = c.comps[0][k1].log_w + c.comps[1][k2].log_w + log_kernel(b1, *comp[0], s.eps_res, skip[0]) +
                           log_kernel(b2, *comp[1], s.eps_res, skip[1]);
                KTerm K;
                for (std::size_t gi = 0; gi < W.groups.size(); ++gi) {
                    const auto& g = W.groups[gi];
                    std::vector<double> v(std::size_t{1} << g.size());
                    for (unsigned m = 0; m < v.size(); ++m) {
                        double pr = 0;
                        for (std::size_t j = 0; j < g.size(); ++j) {
                            const Site& st = W.sites[g[j]];
                            const char x = ((m >> j) & 1U) ? st.alt : st.ref;
                            pr += log_e(x, (*comp[st.hom])[st.pos - c.block[st.hom].begin], s.eps_res);
                        }
                        const auto& sp = scored[W.pair_index[gi][m]];
                        v[m] = pr + (sp.first + sp.second - R.base_loglik);
                    }
                    const double z = lse(v);
                    t += z;
                    for (double& x : v) x = std::exp(x - z);
                    K.post.push_back(v);
                }
                K.total = t;
                terms.push_back(t);
                kt.push_back(std::move(K));
                k1s.push_back(k1);
                k2s.push_back(k2);
            }
        }
        if (terms.empty()) { R.status = "UNRESOLVED_NO_EQUAL_LENGTH_COMPONENT"; continue; }
        const double Z = lse(terms);
        R.log_m_res = R.base_loglik + Z;
        // exact decomposition over the enumerated posterior
        double eh = 0, ealt = 0;
        std::size_t best_k = 0;
        for (std::size_t k = 0; k < kt.size(); ++k) {
            const double pk = std::exp(kt[k].total - Z);
            if (kt[k].total > kt[best_k].total) best_k = k;
            for (std::size_t gi = 0; gi < W.groups.size(); ++gi) {
                for (unsigned m = 0; m < kt[k].post[gi].size(); ++m) {
                    const auto& sp = scored[W.pair_index[gi][m]];
                    eh += pk * kt[k].post[gi][m] * (sp.first + sp.second - R.base_loglik);
                    ealt += pk * kt[k].post[gi][m] * static_cast<double>(__builtin_popcount(m));
                }
            }
        }
        R.eq_fit_gain = eh;
        R.expected_alt = ealt;
        R.kl = R.base_loglik + eh - R.log_m_res;
        // joint MAP state under the dominant component pair
        std::string h0 = c.hap[0], h1 = c.hap[1];
        for (std::size_t gi = 0; gi < W.groups.size(); ++gi) {
            const auto& p = kt[best_k].post[gi];
            const unsigned m = static_cast<unsigned>(std::max_element(p.begin(), p.end()) - p.begin());
            map_masks[ci].push_back(m);
            R.map_alt += __builtin_popcount(m);
            for (std::size_t j = 0; j < W.groups[gi].size(); ++j) {
                if (!((m >> j) & 1U)) continue;
                const Site& st = W.sites[W.groups[gi][j]];
                (st.hom == 0 ? h0 : h1)[st.pos] = st.alt;
            }
        }
        joint_of[ci] = joint.size();
        joint.push_back({h0, h1});
    }
    if (!joint.empty()) {
        const auto js = score_pairs(joint, fragments, lib, opt);
        for (std::size_t ci = 0; ci < cs.size(); ++ci) {
            if (joint_of[ci] == std::numeric_limits<std::size_t>::max()) continue;
            Result& R = res[ci];
            const Work& W = work[ci];
            double sum = 0;
            for (std::size_t gi = 0; gi < W.groups.size(); ++gi) {
                const auto& sp = scored[W.pair_index[gi][map_masks[ci][gi]]];
                sum += sp.first + sp.second - R.base_loglik;
            }
            R.interaction_error = (js[joint_of[ci]].first + js[joint_of[ci]].second - R.base_loglik) - sum;
        }
    }
    return res;
}

}  // namespace vt
