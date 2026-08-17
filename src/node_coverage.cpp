#include "panvar/node_coverage.hpp"

#include "panvar/graph_utils.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <thread>

#include <zlib.h>
#include <kseq.h>
#include <minimap.h>

namespace panvar {
namespace {

KSEQ_INIT(gzFile, gzread)

// Node names in these graphs are integers after graph_sort renumbers them, but nothing guarantees it
// for a graph produced elsewhere. Sort numerically when every name parses as a number and
// lexicographically otherwise, so the dense numbering is deterministic either way -- an unordered_map
// iteration order would make every coverage vector depend on the hash seed.
bool numeric_less(const std::string& a, const std::string& b) {
    const bool na = !a.empty() && a.find_first_not_of("0123456789") == std::string::npos;
    const bool nb = !b.empty() && b.find_first_not_of("0123456789") == std::string::npos;
    if (na && nb) {
        if (a.size() != b.size()) return a.size() < b.size();
        return a < b;
    }
    if (na != nb) return na;
    return a < b;
}

// Byte offset of each step along the spelled path, using the node's ACTUAL sequence length. Deliberately
// not `path_prefix_bp`, which substitutes 1 for a zero-length node: that helper feeds distance
// calculations where the substitution is harmless, but here the offsets must agree with the spelled
// string exactly or every projected interval is shifted.
std::vector<std::size_t> step_offsets(const Graph& graph, const PathRecord& path) {
    std::vector<std::size_t> pref(path.steps.size() + 1, 0);
    for (std::size_t i = 0; i < path.steps.size(); ++i) {
        const auto it = graph.nodes.find(path.steps[i].node_id);
        pref[i + 1] = pref[i] + (it == graph.nodes.end() ? 0 : it->second.sequence.size());
    }
    return pref;
}

} // namespace

NodeIndex build_node_index(const Graph& graph) {
    NodeIndex out;
    out.id.reserve(graph.nodes.size());
    for (const auto& [name, node] : graph.nodes) { (void)node; out.id.push_back(name); }
    std::sort(out.id.begin(), out.id.end(), numeric_less);
    out.length.reserve(out.id.size());
    out.of.reserve(out.id.size() * 2);
    for (std::uint32_t i = 0; i < out.id.size(); ++i) {
        out.of.emplace(out.id[i], i);
        out.length.push_back(static_cast<std::uint32_t>(graph.nodes.at(out.id[i]).sequence.size()));
    }
    return out;
}

PanelCoverage build_panel_coverage(const Graph& graph, const NodeIndex& index) {
    PanelCoverage out;
    out.path_names.reserve(graph.paths.size());
    out.by_path.reserve(graph.paths.size());
    out.path_seq.reserve(graph.paths.size());
    for (const PathRecord& p : graph.paths) {
        out.path_names.push_back(p.name);
        std::vector<std::uint32_t> v(index.size(), 0);
        for (const PathStep& s : p.steps) {
            const auto it = index.of.find(s.node_id);
            if (it != index.of.end()) ++v[it->second];
        }
        out.by_path.push_back(std::move(v));
        out.path_seq.push_back(spell_path_steps_sequence(graph, p.steps));
    }
    return out;
}

SampleCoverage inject_reads(
    const Graph& graph,
    const NodeIndex& index,
    const PanelCoverage& panel,
    const std::vector<std::string>& read_paths,
    const CoverageOptions& options) {

    SampleCoverage out;
    out.node.assign(index.size(), 0.0);
    if (panel.path_seq.empty()) return out;

    // Step offsets per path, once. Projection is a binary search into these plus a forward walk, so
    // this is the only per-path preprocessing needed.
    std::vector<std::vector<std::size_t>> offsets(graph.paths.size());
    for (std::size_t p = 0; p < graph.paths.size(); ++p) offsets[p] = step_offsets(graph, graph.paths[p]);

    mm_idxopt_t idx_opt;
    mm_mapopt_t map_opt;
    if (mm_set_opt(nullptr, &idx_opt, &map_opt) < 0) throw std::runtime_error("coverage: mm_set_opt failed");
    const std::string preset = options.preset.empty() ? std::string("sr") : options.preset;
    if (mm_set_opt(preset.c_str(), &idx_opt, &map_opt) < 0) {
        throw std::runtime_error("coverage: unknown minimap2 preset " + preset);
    }
    // Extended CIGAR, because coverage must be credited only where the read actually MATCHES.
    // Projecting the whole aligned interval smears a read across the private nodes of whatever
    // haplotype it aligned to: on the synthetic locus that put 35% of all coverage (2219 of 6420) on
    // the 99 nodes the sample does not traverse at all, because a read crosses a bubble on the wrong
    // allele with mismatches and still credits it. Crediting only '=' runs makes coverage
    // allele-specific, which is the whole point of a per-node vector.
    map_opt.flag |= (MM_F_CIGAR | MM_F_EQX);
    map_opt.flag &= ~MM_F_OUT_SAM;
    // A read from inside a tandem array matches every copy in every haplotype carrying it. Secondary
    // hits are the evidence here, not noise, so they must survive both the count cap and the
    // score-ratio filter that would otherwise discard near-equal placements.
    map_opt.best_n = std::max<int>(map_opt.best_n, static_cast<int>(options.best_n));
    map_opt.pri_ratio = 0.1f;

    std::vector<const char*> seqs(panel.path_seq.size());
    std::vector<const char*> names(panel.path_names.size());
    for (std::size_t i = 0; i < panel.path_seq.size(); ++i) {
        seqs[i] = panel.path_seq[i].c_str();
        names[i] = panel.path_names[i].c_str();
    }
    mm_idx_t* idx = mm_idx_str(idx_opt.w, idx_opt.k, (idx_opt.flag & MM_I_HPC) != 0 ? 1 : 0,
                               idx_opt.bucket_bits, static_cast<int>(seqs.size()), seqs.data(), names.data());
    if (idx == nullptr) throw std::runtime_error("coverage: could not index the panel paths");
    mm_mapopt_update(&map_opt, idx);
    // Repeat-seed filtering has to be switched off, and this is the single setting that decides
    // whether a tandem array is measurable at all. minimap2 discards minimizers that occur more often
    // than mid_occ, which mm_mapopt_update derives from the index -- sensible for mapping to a genome,
    // exactly wrong here. The panel holds hundreds of near-identical haplotypes and the array is a
    // cycle, so a seed inside it occurs (paths x copies) times: at lpa about 11,000. Left at the
    // default, 23% of reads failed to place and the array nodes came back at 2.59 times their
    // traversal count where every other multiplicity band sat at 12-15.
    map_opt.mid_occ = std::numeric_limits<int>::max();
    map_opt.max_occ = std::numeric_limits<int>::max();
    map_opt.mid_occ_frac = 0.0f;

    const std::size_t nthreads =
        std::max<std::size_t>(1, options.threads != 0 ? options.threads : std::thread::hardware_concurrency());
    std::mutex merge_mu;

    for (const std::string& path : read_paths) {
        gzFile fp = gzopen(path.c_str(), "r");
        if (fp == nullptr) { mm_idx_destroy(idx); throw std::runtime_error("coverage: cannot open " + path); }
        kseq_t* seq = kseq_init(fp);
        std::mutex read_mu;
        std::vector<std::thread> pool;
        for (std::size_t t = 0; t < nthreads; ++t) {
            pool.emplace_back([&] {
                mm_tbuf_t* tbuf = mm_tbuf_init();
                std::vector<double> local(index.size(), 0.0);
                std::unordered_map<std::uint32_t, double> per_read;
                std::vector<std::pair<std::size_t, std::size_t>> spans;
                SampleCoverage stat;
                std::vector<std::pair<std::string, std::string>> batch;
                for (;;) {
                    batch.clear();
                    {
                        std::lock_guard<std::mutex> lock(read_mu);
                        while (batch.size() < 1024 && kseq_read(seq) >= 0) {
                            batch.emplace_back(seq->name.s ? seq->name.s : "r",
                                               std::string(seq->seq.s, seq->seq.l));
                        }
                    }
                    if (batch.empty()) break;
                    for (const auto& [rname, rseq] : batch) {
                        ++stat.reads;
                        if (rseq.empty()) continue;
                        int n_regs = 0;
                        mm_reg1_t* regs = mm_map(idx, static_cast<int>(rseq.size()), rseq.c_str(),
                                                 &n_regs, tbuf, &map_opt, rname.c_str());
                        if (regs == nullptr || n_regs <= 0) { free(regs); continue; }
                        ++stat.aligned;
                        // Uniform split across placements. Each read then contributes a total of one
                        // read's worth of coverage however many copies it matched, which is what keeps
                        // the depth scale interpretable at an array: coverage on the repeat node grows
                        // because MORE READS come from a longer array, not because one read is counted
                        // many times.
                        const double w = options.placement == CoverageOptions::Placement::UniformSplit
                            ? 1.0 / static_cast<double>(n_regs) : 1.0;
                        // Under UnionMax the read's contributions are collected per node first and the
                        // largest kept, so a node reached by several of the read's alignments is
                        // charged once. per_read is cleared per read; it stays small because a read
                        // touches only the handful of nodes it actually spans.
                        per_read.clear();
                        for (int r = 0; r < n_regs; ++r) {
                            const mm_reg1_t& g = regs[r];
                            const std::size_t pi = static_cast<std::size_t>(g.rid);
                            if (pi >= offsets.size()) continue;
                            const std::size_t ts = static_cast<std::size_t>(g.rs);
                            const std::size_t te = static_cast<std::size_t>(g.re);
                            ++stat.placements;
                            stat.bases_placed += static_cast<double>(g.mlen);
                            const auto& pref = offsets[pi];
                            const auto& steps = graph.paths[pi].steps;
                            // Matched target intervals only, walked off the extended CIGAR.
                            spans.clear();
                            if (g.p != nullptr && g.p->n_cigar > 0) {
                                std::size_t tp = ts;
                                for (uint32_t ci = 0; ci < g.p->n_cigar; ++ci) {
                                    const uint32_t op = g.p->cigar[ci] & 0xf;
                                    const uint32_t ln = g.p->cigar[ci] >> 4;
                                    if (op == 7) { spans.emplace_back(tp, tp + ln); tp += ln; }
                                    else if (op == 8 || op == 2 || op == 3) { tp += ln; }
                                    // op 1 (I) consumes query only
                                }
                            } else {
                                spans.emplace_back(ts, te);
                            }
                            // First step whose end passes ts, then forward while the step starts
                            // before te. Both endpoints clipped, so a read covering part of a node
                            // contributes only the bases it actually covers.
                            for (const auto& [ss, se] : spans) {
                            std::size_t lo = static_cast<std::size_t>(
                                std::upper_bound(pref.begin(), pref.end(), ss) - pref.begin());
                            if (lo > 0) --lo;
                            for (std::size_t si = lo; si < steps.size() && pref[si] < se; ++si) {
                                const std::size_t a = std::max(ss, pref[si]);
                                const std::size_t b = std::min(se, pref[si + 1]);
                                if (b <= a) continue;
                                const auto it = index.of.find(steps[si].node_id);
                                if (it == index.of.end()) continue;
                                const double len = options.len_scale
                                    ? std::max<double>(1.0, index.length[it->second]) : 1.0;
                                const double c = w * static_cast<double>(b - a) / len;
                                if (options.placement == CoverageOptions::Placement::UnionMax) {
                                    auto ins = per_read.emplace(it->second, c);
                                    if (!ins.second) ins.first->second = std::max(ins.first->second, c);
                                } else {
                                    local[it->second] += c;
                                }
                            }
                            }
                        }
                        if (options.placement == CoverageOptions::Placement::UnionMax) {
                            for (const auto& [n, c] : per_read) local[n] += c;
                        }
                        free(regs);
                    }
                }
                mm_tbuf_destroy(tbuf);
                std::lock_guard<std::mutex> lock(merge_mu);
                for (std::size_t i = 0; i < local.size(); ++i) out.node[i] += local[i];
                out.reads += stat.reads;
                out.aligned += stat.aligned;
                out.placements += stat.placements;
                out.bases_placed += stat.bases_placed;
            });
        }
        for (std::thread& th : pool) th.join();
        kseq_destroy(seq);
        gzclose(fp);
    }
    mm_idx_destroy(idx);
    return out;
}

CoverageAudit audit_coverage(
    const Graph& graph,
    const NodeIndex& index,
    const PanelCoverage& panel,
    const SampleCoverage& sample,
    const std::vector<std::string>& sample_paths,
    const CoverageOptions& options) {

    CoverageAudit a;

    // Probe: take a synthetic read straight out of a path at a known offset and check it projects onto
    // exactly the nodes that path traverses there. No aligner involved -- this tests the projection
    // arithmetic alone, so a failure here is unambiguous.
    for (std::size_t pi = 0; pi < graph.paths.size() && pi < 8; ++pi) {
        const auto pref = step_offsets(graph, graph.paths[pi]);
        const std::size_t total = pref.back();
        if (total < 400) continue;
        for (int k = 1; k <= 5; ++k) {
            const std::size_t ts = (total / 6) * static_cast<std::size_t>(k);
            const std::size_t te = std::min(total, ts + 150);
            std::vector<std::string> want;
            for (std::size_t si = 0; si < graph.paths[pi].steps.size(); ++si) {
                if (pref[si] < te && pref[si + 1] > ts) want.push_back(graph.paths[pi].steps[si].node_id);
            }
            std::vector<std::string> got;
            std::size_t lo = static_cast<std::size_t>(
                std::upper_bound(pref.begin(), pref.end(), ts) - pref.begin());
            if (lo > 0) --lo;
            for (std::size_t si = lo; si < graph.paths[pi].steps.size() && pref[si] < te; ++si) {
                const std::size_t x = std::max(ts, pref[si]);
                const std::size_t y = std::min(te, pref[si + 1]);
                if (y > x) got.push_back(graph.paths[pi].steps[si].node_id);
            }
            ++a.probe_total;
            if (got == want) ++a.probe_exact;
        }
    }

    // Mass: with len-scale on, coverage * node length is the bases placed on that node, so the sum has
    // to equal the total placed bases. Off by a factor of two means double counting; short means the
    // projection is dropping intervals.
    if (options.len_scale && sample.bases_placed > 0.0) {
        double mass = 0.0;
        for (std::size_t i = 0; i < sample.node.size(); ++i) {
            mass += sample.node[i] * std::max<double>(1.0, index.length[i]);
        }
        // Under UniformSplit each read contributes one read's worth however many alignments it had, so
        // the expected mass is the mean placed bases per alignment times the number of reads. Under
        // UnionMax the same holds by construction. Under All every alignment counts in full.
        const double expect = options.placement == CoverageOptions::Placement::All
            ? sample.bases_placed
            : (sample.placements > 0 ? sample.bases_placed * static_cast<double>(sample.aligned) /
                                           static_cast<double>(sample.placements)
                                     : sample.bases_placed);
        a.mass_ratio = expect > 0.0 ? mass / expect : 0.0;
    }

    // Self-consistency: a sample's own reads against its own traversal vector. Correlation catches
    // gross errors, and the slope is the per-haplotype read depth, which catches the scale errors a
    // correlation cannot see -- a vector that is uniformly twice too large still correlates at 1.0.
    if (!sample_paths.empty()) {
        std::vector<double> expect(index.size(), 0.0);
        std::size_t found = 0;
        for (const std::string& name : sample_paths) {
            for (std::size_t p = 0; p < panel.path_names.size(); ++p) {
                if (panel.path_names[p] != name) continue;
                ++found;
                for (std::size_t i = 0; i < expect.size(); ++i) expect[i] += panel.by_path[p][i];
            }
        }
        if (found > 0) {
            double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
            std::size_t n = 0;
            for (std::size_t i = 0; i < expect.size(); ++i) {
                if (expect[i] == 0.0 && sample.node[i] == 0.0) continue;   // untraversed by either
                const double x = expect[i];
                const double y = sample.node[i];
                sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y; ++n;
            }
            if (n > 2) {
                const double nn = static_cast<double>(n);
                const double cov = sxy / nn - (sx / nn) * (sy / nn);
                const double vx = sxx / nn - (sx / nn) * (sx / nn);
                const double vy = syy / nn - (sy / nn) * (sy / nn);
                a.self_pearson = (vx > 0 && vy > 0) ? cov / std::sqrt(vx * vy) : 0.0;
                a.self_slope = vx > 0 ? cov / vx : 0.0;
            }
        }
    }
    return a;
}


std::vector<std::vector<std::uint32_t>> block_allele_node_vectors(
    const Graph& graph,
    const std::vector<BubblePathIndex>& path_indexes,
    const std::vector<Bubble>& bubbles,
    const Block& block,
    const BlockAlleles& alleles,
    const NodeIndex& index) {

    const Bubble* bubble = nullptr;
    if (block.kind == BlockKind::Bubble) {
        for (const Bubble& b : bubbles) if (b.id == block.bubble_id) { bubble = &b; break; }
    }
    std::vector<std::vector<std::uint32_t>> out(alleles.allele_haplotypes.size());
    std::vector<char> done(out.size(), 0);
    for (std::size_t pi = 0; pi < graph.paths.size(); ++pi) {
        const auto it = alleles.allele_of.find(graph.paths[pi].name);
        if (it == alleles.allele_of.end()) continue;
        const std::size_t ai = it->second;
        if (ai >= out.size() || done[ai]) continue;      // one representative per allele is enough:
        done[ai] = 1;                                    // an allele is a set of identical walks
        std::optional<std::vector<PathStep>> steps;
        if (bubble != nullptr) {
            steps = bubble_steps(graph.paths[pi], path_indexes[pi], *bubble);
            if (steps.has_value() && !block.trim_node.empty() && steps->size() > 1) {
                if (steps->front().node_id == block.trim_node) steps->erase(steps->begin());
                else if (steps->back().node_id == block.trim_node) steps->pop_back();
            }
        } else if (block.kind == BlockKind::Flank) {
            const bool leading = block.source.empty();
            steps = flank_steps(graph.paths[pi], path_indexes[pi], leading ? block.sink : block.source, leading);
        } else {
            steps = interval_interior_steps(graph.paths[pi], path_indexes[pi], block.source, block.sink);
        }
        if (!steps.has_value()) continue;
        std::vector<std::uint32_t> v(index.size(), 0);
        for (const PathStep& st : *steps) {
            const auto ni = index.of.find(st.node_id);
            if (ni != index.of.end()) ++v[ni->second];
        }
        out[ai] = std::move(v);
    }
    for (auto& v : out) if (v.empty()) v.assign(index.size(), 0);   // the bypass allele traverses nothing
    return out;
}

std::vector<CoverageScore> score_block_by_coverage(
    const std::vector<std::vector<std::uint32_t>>& allele_vec,
    const std::vector<std::size_t>& allele_bp,
    const SampleCoverage& sample,
    const NodeIndex& index,
    double lambda,
    double overdispersion) {

    // The block's allele vectors are already in node-index space, so the index itself is not needed
    // here. Kept in the signature because every other scorer takes it and callers pass it positionally.
    (void)index;

    // Only the nodes some allele of this block traverses. Everything else is another block's business,
    // and including it would add the same constant to every candidate while diluting the comparison.
    std::vector<std::uint32_t> nodes;
    for (const auto& v : allele_vec) {
        for (std::size_t i = 0; i < v.size(); ++i) if (v[i] > 0) nodes.push_back(static_cast<std::uint32_t>(i));
    }
    std::sort(nodes.begin(), nodes.end());
    nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());

    double obs_norm2 = 0.0, obs_sum = 0.0;
    for (const std::uint32_t n : nodes) { obs_norm2 += sample.node[n] * sample.node[n]; obs_sum += sample.node[n]; }
    const double obs_norm = std::sqrt(std::max(1e-12, obs_norm2));
    const double obs_mean = nodes.empty() ? 0.0 : obs_sum / static_cast<double>(nodes.size());
    const double mu = std::max(0.01, 0.02 * lambda);

    auto lognb = [&](double x, double mean) {
        if (mean <= 0.0) mean = 1e-9;
        if (overdispersion <= 0.0) return -mean + x * std::log(mean) - std::lgamma(x + 1.0);
        const double phi = overdispersion;
        return std::lgamma(x + phi) - std::lgamma(phi) - std::lgamma(x + 1.0)
             + phi * std::log(phi / (phi + mean)) + x * std::log(mean / (phi + mean));
    };

    std::vector<CoverageScore> out;
    const std::size_t na = allele_vec.size();
    out.reserve(na * (na + 1) / 2);
    for (std::size_t a = 0; a < na; ++a) {
        for (std::size_t b = a; b < na; ++b) {
            CoverageScore s;
            s.allele1 = a; s.allele2 = b;
            s.bp = (a < allele_bp.size() ? allele_bp[a] : 0) + (b < allele_bp.size() ? allele_bp[b] : 0);
            double ll = 0.0, dot = 0.0, pn2 = 0.0;
            double pmean = 0.0;
            for (const std::uint32_t n : nodes) pmean += allele_vec[a][n] + allele_vec[b][n];
            pmean = nodes.empty() ? 0.0 : pmean / static_cast<double>(nodes.size());
            double cdot = 0.0, cpx = 0.0, cpy = 0.0;
            for (const std::uint32_t n : nodes) {
                const double m = static_cast<double>(allele_vec[a][n] + allele_vec[b][n]);
                const double o = sample.node[n];
                ll += lognb(o, lambda * m + mu);
                dot += m * o; pn2 += m * m;
                // Pearson: the same score on centred vectors, which removes the component every
                // candidate carries and is what compresses cosine at a repeat-heavy block.
                cdot += (m - pmean) * (o - obs_mean);
                cpx += (m - pmean) * (m - pmean);
                cpy += (o - obs_mean) * (o - obs_mean);
            }
            s.loglik = ll;
            s.cosine = (pn2 > 0.0) ? dot / (std::sqrt(pn2) * obs_norm) : 0.0;
            s.pearson = (cpx > 0.0 && cpy > 0.0) ? cdot / std::sqrt(cpx * cpy) : 0.0;
            out.push_back(s);
        }
    }
    return out;
}

void write_node_coverage(
    const std::string& out_prefix,
    const NodeIndex& index,
    const PanelCoverage& panel,
    const SampleCoverage& sample) {

    const std::string p1 = out_prefix + ".nodecov.sample.tsv";
    std::ofstream f(p1);
    if (!f) throw std::runtime_error("coverage: cannot write " + p1);
    f << "node\tlength\tcoverage\n";
    for (std::size_t i = 0; i < index.size(); ++i) {
        f << index.id[i] << '\t' << index.length[i] << '\t' << sample.node[i] << '\n';
    }
    const std::string p2 = out_prefix + ".nodecov.panel.tsv";
    std::ofstream g(p2);
    if (!g) throw std::runtime_error("coverage: cannot write " + p2);
    g << "path";
    for (std::size_t i = 0; i < index.size(); ++i) g << '\t' << index.id[i];
    g << '\n';
    for (std::size_t p = 0; p < panel.path_names.size(); ++p) {
        g << panel.path_names[p];
        for (std::size_t i = 0; i < index.size(); ++i) g << '\t' << panel.by_path[p][i];
        g << '\n';
    }
}

} // namespace panvar
