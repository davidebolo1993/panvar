#include "panvar/genotype_dosage_repair.hpp"

#include "panvar/genotype_repair.hpp"
#include "panvar/ls_hmm.hpp"
#include "panvar/parallel.hpp"

#include <algorithm>
#include <cstring>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <tuple>
#include <unordered_set>

namespace panvar::dosage {

std::string sha256_hex(const std::string& bytes) { return sha256_bytes(bytes); }

namespace {

double seconds_since(std::chrono::steady_clock::time_point t) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
}

// numpy's float64 add.reduce over a contiguous run (pairwise summation, blocks of 8, PW_BLOCKSIZE 128). The prototype
// summed with it, so the port uses it wherever the prototype's order is part of a replayed value.
double numpy_sum(const double* a, std::size_t n) {
    if (n < 8) {
        double r = 0.0;
        for (std::size_t i = 0; i < n; ++i) r += a[i];
        return r;
    }
    if (n <= 128) {
        double r[8];
        for (int j = 0; j < 8; ++j) r[j] = a[j];
        std::size_t i = 8;
        for (; i < n - (n % 8); i += 8) for (int j = 0; j < 8; ++j) r[j] += a[i + j];
        double res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; i < n; ++i) res += a[i];
        return res;
    }
    std::size_t n2 = n / 2;
    n2 -= n2 % 8;
    return numpy_sum(a, n2) + numpy_sum(a + n2, n - n2);
}

std::string reverse_complement_bases(const std::string& s) {
    std::string r(s.rbegin(), s.rend());
    for (char& c : r) {
        switch (c) {
            case 'A': c = 'T'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            case 'T': c = 'A'; break;
            default: break;
        }
    }
    return r;
}

using Steps = std::vector<std::string>;   // "<node><+|->"

std::string spell(const Steps& steps, const NodeSequences& nodes) {
    std::string out;
    for (const std::string& s : steps) {
        const auto found = nodes.find(s.substr(0, s.size() - 1));
        if (found == nodes.end()) throw std::runtime_error("unit library: no sequence for node step " + s);
        out += s.back() == '+' ? found->second : reverse_complement_bases(found->second);
    }
    return out;
}

std::size_t node_bp(const std::string& step, const NodeSequences& nodes) {
    const auto found = nodes.find(step.substr(0, step.size() - 1));
    if (found == nodes.end()) throw std::runtime_error("unit library: no sequence for node step " + step);
    return found->second.size();
}

std::vector<Steps> step_walks(const CanonicalBlockAllele& allele) {
    std::vector<Steps> out;
    for (const std::vector<NodeStep>& w : all_walks(allele)) {
        Steps s;
        for (const NodeStep& n : w) s.push_back(n.node + (n.forward ? "+" : "-"));
        out.push_back(std::move(s));
    }
    return out;
}

template <typename T>
bool all_same(const std::vector<T>& v) {
    return std::all_of(v.begin(), v.end(), [&](const T& x) { return x == v.front(); });
}

// itertools.combinations_with_replacement(pool, k), in its order
void combinations_with_replacement(std::size_t n, std::size_t k, const std::function<void(const std::vector<std::size_t>&)>& fn) {
    std::vector<std::size_t> idx(k, 0);
    if (k == 0) { fn(idx); return; }
    if (n == 0) return;
    for (;;) {
        fn(idx);
        std::size_t i = k;
        while (i > 0 && idx[i - 1] == n - 1) --i;
        if (i == 0) return;
        const std::size_t v = idx[i - 1] + 1;
        for (std::size_t j = i - 1; j < k; ++j) idx[j] = v;
    }
}

} // namespace

// ---------------------------------------------------------------------------------------------- unit library

std::string UnitLibrary::block_sequence(const std::vector<std::size_t>& tuple) const {
    std::string s = flank_left;
    for (std::size_t u : tuple) s += units.at(u).sequence;
    s += flank_right;
    if (s.size() < trim_left + trim_right) throw std::runtime_error("unit tuple shorter than the anchor trims");
    return s.substr(trim_left, s.size() - trim_left - trim_right);
}

UnitLibrary build_unit_library(const PreparedGenotypeBlock& block, const NodeSequences& nodes, const RepairSettings& settings) {
    const std::size_t na = block.alleles.size();
    std::vector<std::vector<Steps>> walks(na);
    std::vector<Steps> primary(na);
    for (std::size_t a = 0; a < na; ++a) {
        walks[a] = step_walks(block.alleles[a]);
        if (walks[a].empty()) throw std::runtime_error("unit library: allele without a walk in " + block.block_id);
        primary[a] = walks[a].front();
    }
    // entry: the first revisited step, in walk order, of the allele that revisits most (first such allele)
    std::size_t top = 0, best = 0;
    for (std::size_t a = 0; a < na; ++a) {
        std::map<std::string, std::size_t> c;
        std::size_t mx = 0;
        for (const std::string& s : primary[a]) mx = std::max(mx, ++c[s]);
        if (a == 0 || mx > best) { best = mx; top = a; }
    }
    UnitLibrary lib;
    {
        std::map<std::string, std::size_t> c;
        for (const std::string& s : primary[top]) ++c[s];
        for (const std::string& s : primary[top]) if (c[s] > 1) { lib.entry_step = s; break; }
        if (lib.entry_step.empty()) throw std::runtime_error("unit library: no revisited node in " + block.block_id);
    }
    std::vector<const Steps*> allw;
    for (const auto& ws : walks) for (const Steps& w : ws) allw.push_back(&w);
    std::size_t lp = 0;
    for (;;) {
        bool ok = std::all_of(allw.begin(), allw.end(), [&](const Steps* w) { return w->size() > lp; });
        if (!ok) break;
        std::set<std::string> at;
        for (const Steps* w : allw) at.insert((*w)[lp]);
        if (at.size() != 1 || (*allw.front())[lp] == lib.entry_step) break;
        ++lp;
    }
    std::size_t ls = 0;
    for (;;) {
        bool ok = std::all_of(allw.begin(), allw.end(), [&](const Steps* w) { return w->size() > ls + lp; });
        if (!ok) break;
        std::set<std::string> at;
        for (const Steps* w : allw) at.insert((*w)[w->size() - 1 - ls]);
        if (at.size() != 1) break;
        ++ls;
    }
    struct Out { Steps left, right, pre; std::vector<Steps> modules; };
    std::vector<Out> out(na);
    for (std::size_t a = 0; a < na; ++a) {
        const Steps& w = primary[a];
        Out& o = out[a];
        o.left.assign(w.begin(), w.begin() + static_cast<long>(lp));
        o.right.assign(w.end() - static_cast<long>(ls), w.end());
        const Steps mid(w.begin() + static_cast<long>(lp), w.end() - static_cast<long>(ls));
        std::vector<std::size_t> cut;
        for (std::size_t i = 0; i < mid.size(); ++i) if (mid[i] == lib.entry_step) cut.push_back(i);
        o.pre.assign(mid.begin(), mid.begin() + static_cast<long>(cut.empty() ? mid.size() : cut.front()));
        for (std::size_t i = 0; i < cut.size(); ++i) {
            const std::size_t end = i + 1 < cut.size() ? cut[i + 1] : mid.size();
            o.modules.emplace_back(mid.begin() + static_cast<long>(cut[i]), mid.begin() + static_cast<long>(end));
        }
    }
    // optional steps: present in some module copies and absent from others
    std::size_t module_copies = 0;
    std::map<std::string, std::size_t> present;
    for (const Out& o : out) {
        for (const Steps& m : o.modules) {
            ++module_copies;
            for (const std::string& s : std::set<std::string>(m.begin(), m.end())) ++present[s];
        }
    }
    const auto optional = [&](const std::string& s) { return present[s] < module_copies; };
    const auto state_of = [&](const Steps& m) {
        std::size_t run = 0;
        bool large = false;
        for (const std::string& s : m) {
            if (optional(s)) { run += node_bp(s, nodes); }
            else { large |= run >= settings.large_optional_run_bp; run = 0; }
        }
        large |= run >= settings.large_optional_run_bp;
        return large ? 'I' : 'O';
    };
    // trims: the model block is the spelled walk trimmed by one common anchor overhang
    std::set<std::pair<std::size_t, std::size_t>> trims;
    std::set<std::string> left_trims, right_trims, flank_lefts, flank_rights;
    for (std::size_t a = 0; a < na; ++a) {
        const std::string sp = spell(primary[a], nodes);
        const std::string& seq = block.alleles[a].sequence;
        const std::size_t L = sp.find(seq);
        if (L == std::string::npos || sp.find(seq, L + 1) != std::string::npos) {
            throw std::runtime_error("unit library: allele " + std::to_string(a) + " does not spell its model sequence uniquely");
        }
        const std::size_t R = sp.size() - L - seq.size();
        trims.insert({L, R});
        left_trims.insert(sp.substr(0, L));
        right_trims.insert(sp.substr(sp.size() - R));
        Steps lpre = out[a].left;
        lpre.insert(lpre.end(), out[a].pre.begin(), out[a].pre.end());
        flank_lefts.insert(spell(lpre, nodes));
        flank_rights.insert(spell(out[a].right, nodes));
    }
    if (trims.size() != 1 || left_trims.size() != 1 || right_trims.size() != 1 || flank_lefts.size() != 1 || flank_rights.size() != 1) {
        throw std::runtime_error("unit library: alleles disagree on their flanks or anchor trims in " + block.block_id);
    }
    lib.trim_left = trims.begin()->first;
    lib.trim_right = trims.begin()->second;
    lib.flank_left = *flank_lefts.begin();
    lib.flank_right = *flank_rights.begin();
    const std::string left_trim = *left_trims.begin();
    if (out.front().right.empty()) throw std::runtime_error("unit library: empty right flank in " + block.block_id);
    const std::string right0 = out.front().right.front();
    std::set<std::pair<std::string, std::string>> edges;
    for (const Steps* w : allw) for (std::size_t i = 0; i + 1 < w->size(); ++i) edges.insert({(*w)[i], (*w)[i + 1]});
    std::map<std::string, std::size_t> canon;
    lib.panel.resize(na);
    for (std::size_t a = 0; a < na; ++a) {
        for (std::size_t pos = 0; pos < out[a].modules.size(); ++pos) {
            const Steps& m = out[a].modules[pos];
            const std::string s = spell(m, nodes);
            const char st = state_of(m);
            auto it = canon.find(s);
            if (it == canon.end()) {
                it = canon.emplace(s, lib.units.size()).first;
                Unit u;
                u.id = lib.units.size();
                u.sequence = s;
                u.state = st;
                lib.units.push_back(std::move(u));
            }
            Unit& u = lib.units[it->second];
            if (u.state != st) throw std::runtime_error("unit library: one unit sequence with two states");
            if (std::find(u.walks.begin(), u.walks.end(), m) == u.walks.end()) u.walks.push_back(m);
            u.carriers.push_back({a, pos + 1, out[a].modules.size()});
            lib.panel[a].push_back(it->second);
        }
        lib.kmax = std::max(lib.kmax, lib.panel[a].size());
    }
    for (Unit& u : lib.units) {
        for (const Steps& w : u.walks) {
            u.internal_capable |= edges.count({w.back(), lib.entry_step}) > 0;
            u.terminal_capable |= edges.count({w.back(), right0}) > 0;
        }
        u.anchor_compatible = (lib.flank_left + u.sequence).substr(0, lib.trim_left) == left_trim;
    }
    return lib;
}

// ---------------------------------------------------------------------------------------------- exact scoring

namespace {

struct GuardFailure : std::runtime_error { using std::runtime_error::runtime_error; };

struct ExactScore {
    double fragment = 0.0, exposure = 0.0, likelihood = 0.0;
};

std::string splice(const ReportedCall& call, int h, std::size_t block, const std::string& seq) {
    const SequenceInterval iv = call.intervals[h][block];
    return call.haplotypes[h].substr(0, iv.begin) + seq + call.haplotypes[h].substr(iv.end);
}

// Exact production scores of block candidates spliced into the reported homologues. The reported diplotype rides
// along in every chunk and must reproduce its persisted total, fragment and exposure terms within the guard tolerance
// before any candidate score is used (GuardFailure). A resource cap halves the chunk, down to 8 haplotypes.
std::vector<ExactScore> score_block_candidates(const ScoringInputs& in, std::size_t block,
                                               const std::vector<std::pair<std::string, std::string>>& ctx,
                                               const RepairSettings& rs, ExactScore& original,
                                               std::vector<std::vector<long long>>* grids = nullptr) {
    const ReportedCall& call = *in.call;
    const double winner_prior = log_unordered_sequence_pair_prior(*in.model, call.haplotypes[0], call.haplotypes[1], in.switch_probability);
    std::vector<ExactScore> out(ctx.size());
    std::vector<double> backgrounds;
    if (grids != nullptr) {
        grids->assign(ctx.size(), {});
        for (const ReadPair& f : *in.reads) backgrounds.push_back(log_background_probability(f));
    }
    std::size_t chunk_cap = rs.chunk_haplotypes;
    std::size_t next = 0;
    bool guarded = false;
    while (next < ctx.size() || !guarded) {
        std::map<std::pair<std::string, int>, std::size_t> slot;
        std::size_t end = next;
        for (; end < ctx.size(); ++end) {
            const std::size_t added = !slot.count({ctx[end].first, 1}) + !slot.count({ctx[end].second, 2});
            if (slot.size() + added > chunk_cap && end > next) break;
            for (const auto& key : {std::make_pair(ctx[end].first, 1), std::make_pair(ctx[end].second, 2)}) {
                if (!slot.count(key)) { const std::size_t n = slot.size(); slot[key] = n; }
            }
        }
        std::vector<MosaicPath> paths(slot.size());
        std::vector<std::string> sequences(slot.size());
        for (const auto& [key, i] : slot) {
            paths[i] = key.second == 1 ? in.first : in.second;
            sequences[i] = splice(call, key.second - 1, block, key.first);
        }
        paths.push_back(in.first); sequences.push_back(call.haplotypes[0]);
        paths.push_back(in.second); sequences.push_back(call.haplotypes[1]);
        const InternedHaplotypes interned = intern_haplotypes(paths, sequences);
        std::vector<HaplotypeCandidate> candidates;
        for (std::size_t k = next; k < end; ++k) {
            const std::size_t x = interned.index_of[slot.at({ctx[k].first, 1})], y = interned.index_of[slot.at({ctx[k].second, 2})];
            candidates.push_back({std::min(x, y), std::max(x, y)});
        }
        {
            const std::size_t x = interned.index_of[paths.size() - 2], y = interned.index_of[paths.size() - 1];
            candidates.push_back({std::min(x, y), std::max(x, y)});
        }
        FragmentScoringResult scored;
        try {
            scored = score_diplotypes(interned, candidates, *in.reads, *in.library, in.options);
        } catch (const ResourceExhausted&) {
            if (chunk_cap <= 8) throw;
            chunk_cap /= 2;
            continue;
        }
        const DiplotypeScore& w = scored.scores.back();
        const double dt = w.log_likelihood + winner_prior - in.persisted.total, df = w.fragment_term - in.persisted.fragment,
                     de = w.exposure_term - in.persisted.exposure;
        if (!(std::fabs(dt) <= rs.guard_tolerance && std::fabs(df) <= rs.guard_tolerance && std::fabs(de) <= rs.guard_tolerance)) {
            std::ostringstream m;
            m << std::setprecision(17) << "the reported call did not reproduce (total " << dt << ", fragment " << df << ", exposure " << de << ")";
            throw GuardFailure(m.str());
        }
        original = {w.fragment_term, w.exposure_term, w.log_likelihood};
        guarded = true;
        for (std::size_t k = next; k < end; ++k) {
            const DiplotypeScore& s = scored.scores[k - next];
            out[k] = {s.fragment_term, s.exposure_term, s.log_likelihood};
            if (grids != nullptr) {
                // per-fragment diploid contributions on a 1e-9-nat grid (the prototype's fragment_vector_fnv64_1e-9)
                const HaplotypeCandidate& c = candidates[k - next];
                std::vector<long long>& g = (*grids)[k];
                g.resize(in.reads->size());
                for (std::size_t f = 0; f < in.reads->size(); ++f) {
                    const double m = log_add_exp(scored.mass_of(f, c.first), scored.mass_of(f, c.second));
                    g[f] = std::llround(log_fragment_contribution(m, backgrounds[f], *in.library) * 1e9);
                }
            }
        }
        next = end;
    }
    return out;
}

// The context lifecycle (Sol, 2026-09-28, part C). Same scores as score_block_candidates, bit for bit (gated), from
// one PreparedScoringContext: the reported homologues are placed once with retained censuses and must reproduce their
// persisted scores first; each distinct repaired haplotype is placed once by the exact block delta from its homologue
// and reused by every diplotype that holds it; the emission and origin caches persist across batches; a column is
// released once no pending candidate needs it, and released slots are reused. Caps are cumulative over the repair
// (the profile's repair caps); a cap refuses the whole repair (ResourceExhausted), with no partial output.
std::vector<ExactScore> score_block_candidates_context(const ScoringInputs& in, std::size_t block,
                                                       const std::vector<std::pair<std::string, std::string>>& ctx,
                                                       const RepairSettings& rs, ExactScore& original,
                                                       std::vector<std::vector<long long>>* grids, RepairScoringStats* stats) {
    const auto t0 = std::chrono::steady_clock::now();
    RepairScoringStats local;
    RepairScoringStats& st = stats ? *stats : local;
    st.lifecycle = "context";
    st.candidates = ctx.size();
    const ReportedCall& call = *in.call;
    FragmentScoringOptions opt = in.options;
    if (in.repair_max_physical_origins) opt.max_physical_origins = in.repair_max_physical_origins;
    if (in.repair_max_state_classes) opt.max_state_classes = in.repair_max_state_classes;
    if (in.repair_max_piece_hits) opt.max_piece_hits = in.repair_max_piece_hits;
    if (in.repair_max_latent_pairs) opt.max_latent_alignment_pairs = in.repair_max_latent_pairs;
    PreparedScoringContext context(*in.reads, *in.library, opt);
    if (in.repair_origin_cache_path) {
        context.enable_origin_cache_path(true);
        context.set_origin_cache_limits(in.repair_origin_cache_bytes, 0);
    }
    context.begin_transient_round();
    const std::vector<std::size_t> base = context.add_haplotypes({call.haplotypes[0], call.haplotypes[1]});
    context.retain_census(base[0]);
    if (base[1] != base[0]) context.retain_census(base[1]);
    // the reproduction guard, first and exact, before any repaired haplotype exists
    {
        const FragmentScoringResult g = context.score_candidates({{std::min(base[0], base[1]), std::max(base[0], base[1])}});
        const double prior = log_unordered_sequence_pair_prior(*in.model, call.haplotypes[0], call.haplotypes[1], in.switch_probability);
        const DiplotypeScore& w = g.scores.front();
        const double dt = w.log_likelihood + prior - in.persisted.total, df = w.fragment_term - in.persisted.fragment,
                     de = w.exposure_term - in.persisted.exposure;
        if (!(std::fabs(dt) <= rs.guard_tolerance && std::fabs(df) <= rs.guard_tolerance && std::fabs(de) <= rs.guard_tolerance)) {
            std::ostringstream m;
            m << std::setprecision(17) << "the reported call did not reproduce (total " << dt << ", fragment " << df << ", exposure " << de << ")";
            throw GuardFailure(m.str());
        }
        original = {w.fragment_term, w.exposure_term, w.log_likelihood};
    }
    st.seconds_base = seconds_since(t0);
    std::array<std::string, 2> own;
    for (int h = 0; h < 2; ++h) {
        const SequenceInterval iv = call.intervals[h][block];
        own[h] = call.haplotypes[h].substr(iv.begin, iv.end - iv.begin);
    }
    using Key = std::pair<std::string, int>;
    std::map<Key, std::size_t> pending;           // references still to be scored, per repaired haplotype
    for (const auto& c : ctx) {
        for (int h = 0; h < 2; ++h) {
            const std::string& seq = h == 0 ? c.first : c.second;
            if (seq == own[h]) { ++st.base_references; continue; }
            ++pending[{seq, h}];
            ++st.column_references;
        }
    }
    st.unique_repaired_haplotypes = pending.size();
    std::map<Key, std::size_t> slot_of;
    std::vector<double> backgrounds;
    if (grids != nullptr) {
        grids->assign(ctx.size(), {});
        for (const ReadPair& f : *in.reads) backgrounds.push_back(log_background_probability(f));
    }
    std::vector<ExactScore> out(ctx.size());
    const std::size_t B = std::max<std::size_t>(1, rs.batch_candidates);
    for (std::size_t next = 0; next < ctx.size(); next += B) {
        const std::size_t end = std::min(ctx.size(), next + B);
        const auto tp = std::chrono::steady_clock::now();
        std::vector<HaplotypeCandidate> batch;
        for (std::size_t k = next; k < end; ++k) {
            std::array<std::size_t, 2> slots{};
            for (int h = 0; h < 2; ++h) {
                const std::string& seq = h == 0 ? ctx[k].first : ctx[k].second;
                if (seq == own[h]) { slots[h] = base[h]; continue; }
                auto it = slot_of.find({seq, h});
                if (it == slot_of.end()) {
                    const SequenceInterval iv = call.intervals[h][block];
                    BlockDeltaWork work;
                    const std::size_t slot = context.add_repaired_haplotype(base[h], iv.begin, iv.end, seq, &work, false);
                    st.origins_reused += work.origins_reused;
                    st.origins_recomputed += work.origins_recomputed;
                    st.origin_cache_hits += work.origin_cache_hits;
                    st.origin_cache_misses += work.origin_cache_misses;
                    ++st.columns_built;
                    it = slot_of.emplace(Key{seq, h}, slot).first;
                }
                slots[h] = it->second;
            }
            batch.push_back({std::min(slots[0], slots[1]), std::max(slots[0], slots[1])});
        }
        st.seconds_placement += seconds_since(tp);
        const auto tc = std::chrono::steady_clock::now();
        const FragmentScoringResult scored = context.score_candidates(batch);
        for (std::size_t k = next; k < end; ++k) {
            const DiplotypeScore& s = scored.scores[k - next];
            out[k] = {s.fragment_term, s.exposure_term, s.log_likelihood};
            if (grids != nullptr) {
                const HaplotypeCandidate& c = batch[k - next];
                std::vector<long long>& g = (*grids)[k];
                g.resize(in.reads->size());
                for (std::size_t f = 0; f < in.reads->size(); ++f) {
                    const double m = log_add_exp(context.mass_of(f, c.first), context.mass_of(f, c.second));
                    g[f] = std::llround(log_fragment_contribution(m, backgrounds[f], *in.library) * 1e9);
                }
            }
        }
        st.seconds_combination += seconds_since(tc);
        ++st.batches;
        st.peak_transient_bytes = std::max(st.peak_transient_bytes, context.transient_payload_bytes());
        // release every repaired column no pending candidate needs
        bool any_released = false;
        for (std::size_t k = next; k < end; ++k) {
            for (int h = 0; h < 2; ++h) {
                const std::string& seq = h == 0 ? ctx[k].first : ctx[k].second;
                if (seq == own[h]) continue;
                if (--pending[{seq, h}] == 0) any_released = true;
            }
        }
        if (any_released) {
            std::vector<std::size_t> keep = {base[0], base[1]};
            for (auto it = slot_of.begin(); it != slot_of.end();) {
                if (pending[it->first] == 0) it = slot_of.erase(it);
                else { keep.push_back(it->second); ++it; }
            }
            context.release_transient(keep);
            ++st.releases;
        }
    }
    const FragmentScoringCounters& cc = context.cumulative_counters();
    const EmissionCacheCounters& ec = context.cumulative_emission_cache();
    st.physical_origins = cc.physical_origins;
    st.piece_hits = cc.piece_hits;
    st.state_classes = cc.statistical_state_classes;
    st.dp_executions = ec.executed_dps;
    st.dp_requests = ec.dp_cache_requests;
    st.dp_hits = ec.dp_cache_hits;
    st.prefilters = ec.executed_prefilters;
    st.emission_cache_bytes = ec.payload_bytes;
    st.emission_cache_entries = ec.entries;
    st.origin_cache_entries = context.origin_cache_entries();
    st.origin_cache_bytes = context.origin_cache_bytes();
    st.slots = context.unique_haplotypes();
    st.haplotypes_created = context.haplotypes_created();
    st.footprint_peak_bytes = physical_footprint().peak;
    st.seconds_total = seconds_since(t0);
    if (context.refused()) throw ResourceExhausted("the repair context refused at a declared cap");
    return out;
}

// the tested features of a block for this call, recomputed exactly as audit_dosage builds them
std::vector<std::vector<Kmer>> tested_feature_kmers(const BlockAudit& audit, const BlockFeatureModel& fm, const ReportedCall& call,
                                                    const DosageSettings& settings, std::vector<std::size_t>& feature_index) {
    const std::size_t k = settings.kmer;
    std::unordered_set<Kmer> outside;
    for (int h = 0; h < 2; ++h) {
        const SequenceInterval iv = call.intervals[h][audit.block];
        for (const KmerAt& x : canonical_kmers(call.haplotypes[h], k)) {
            if (x.position + k <= iv.begin || x.position >= iv.end) outside.insert(x.kmer);
        }
    }
    std::vector<std::vector<Kmer>> tested;
    for (std::size_t j = 0; j < fm.features.size(); ++j) {
        std::vector<Kmer> kept;
        for (Kmer x : fm.features[j].kmers) if (!outside.count(x)) kept.push_back(x);
        if (kept.size() >= settings.min_feature_kmers) { tested.push_back(std::move(kept)); feature_index.push_back(j); }
    }
    if (tested.size() != audit.features.size()) throw std::runtime_error("repair: tested features do not match the audit");
    for (std::size_t t = 0; t < tested.size(); ++t) {
        if (tested[t].size() != audit.features[t].kmers) throw std::runtime_error("repair: tested feature sizes do not match the audit");
    }
    return tested;
}

bool admissible(const std::vector<int>& set, long x) { return std::find(set.begin(), set.end(), static_cast<int>(x)) != set.end(); }

using Genotype = std::pair<std::string, std::string>;
Genotype unordered(const std::string& a, const std::string& b) { return a <= b ? Genotype{a, b} : Genotype{b, a}; }

} // namespace

// ---------------------------------------------------------------------------------------------- repair

void classify_repair_scope(RepairResult& res, CensusRoute route) {
    // Sol, 2026-09-28: dosage content, exact diplotype and allocation completeness are separate facts
    const bool pair = route == CensusRoute::Pair;
    res.dosage_content_scope = pair || res.content_complete ? "COMPLETE" : "CONDITIONAL";
    res.diplotype_scope = res.scope == "COMPLETE" ? "COMPLETE" : "CONDITIONAL";
    res.allocation_equivalence_scope = pair ? "NOT_APPLICABLE" : res.diplotype_scope;
    if (res.outcome != "REPAIRED") {
        res.repair_class = "TRIGGERED_UNREPAIRED";
        res.serialization_representative = false;
        return;
    }
    if (res.dosage_content_scope != "COMPLETE") {
        // the dosage content itself is not established: keep the original genotype and claim no correction
        res.repair_class = "TRIGGERED_UNREPAIRED";
        res.outcome = "TRIGGERED_UNREPAIRED";
        res.reason = "dosage_content_conditional";
        res.chosen = -1;
        for (BlockCandidate& c : res.candidates) c.chosen = false;   // no candidate is claimed
        res.serialization_representative = false;
        return;
    }
    if (res.diplotype_scope == "COMPLETE") {
        res.repair_class = "COMPLETE_DIPLOTYPE_REPAIR";
        res.serialization_representative = res.representative_only;
    } else {
        res.repair_class = "CONDITIONAL_DOSAGE_CONTENT";
        res.serialization_representative = true;   // the exact allele sequences and allocation are not authoritative
    }
}

void choose_repair(RepairResult& res, CensusRoute route, double tie_tolerance_nats) {
    // merge by unordered spelled block genotype; the top set is every genotype within the tolerance of the maximum
    double best = -std::numeric_limits<double>::infinity();
    for (const BlockCandidate& c : res.candidates) best = std::max(best, c.likelihood);
    std::map<Genotype, std::vector<std::size_t>> top;
    for (std::size_t i = 0; i < res.candidates.size(); ++i) {
        if (res.candidates[i].likelihood >= best - tie_tolerance_nats) {
            res.candidates[i].in_top_set = true;
            top[unordered(res.candidates[i].ctx1_sequence, res.candidates[i].ctx2_sequence)].push_back(i);
        }
    }
    res.top_set_size = top.size();
    const auto refuse_tie = [&]() {
        res.outcome = "TRIGGERED_UNREPAIRED";
        res.reason = "top_genotype_tie";
        res.chosen = -1;
    };
    if (top.empty()) { res.outcome = "TRIGGERED_UNREPAIRED"; res.reason = "no_admissible_candidate"; return; }
    std::size_t pick = 0;
    if (route == CensusRoute::Pair) {
        if (top.size() > 1) return refuse_tie();
        const std::vector<std::size_t>& g = top.begin()->second;
        pick = g.front();
        for (std::size_t i : g) if (res.candidates[i].likelihood > res.candidates[pick].likelihood) pick = i;
        if (g.size() > 1) {   // both orientations of one genotype tied: orientation carries no block-genotype information
            res.orientation_tie = true;
            for (std::size_t i : g) if (res.candidates[i].orientation == "AB") pick = i;
        }
        res.candidates[pick].chosen_content = true;
    } else {
        // the content multiset of the top set; several multisets across distinct tied genotypes refuse
        std::set<std::vector<std::size_t>> multisets;
        for (const auto& [gt, idx] : top) {
            std::size_t b = idx.front();
            for (std::size_t i : idx) if (res.candidates[i].likelihood > res.candidates[b].likelihood) b = i;
            multisets.insert(res.candidates[b].multiset);
        }
        if (multisets.size() > 1) return refuse_tie();
        const std::vector<std::size_t> content = *multisets.begin();
        // allocation within the content multiset (Sol, 2026-09-28): the panel-structure prior first (exact), then every
        // arrangement within the tie tolerance of the best likelihood among them is one allocation-equivalence set;
        // the smallest canonical SHA-256 pair only serialises it and is not a resolved allocation
        double best_prior = -std::numeric_limits<double>::infinity();
        for (BlockCandidate& c : res.candidates) {
            if (c.multiset != content) continue;
            c.chosen_content = true;
            best_prior = std::max(best_prior, c.panel_log_prior);
        }
        double best_ll = -std::numeric_limits<double>::infinity();
        for (const BlockCandidate& c : res.candidates) {
            if (c.multiset == content && c.panel_log_prior == best_prior) best_ll = std::max(best_ll, c.likelihood);
        }
        std::set<Genotype> tied;
        std::pair<std::string, std::string> smallest;
        bool first = true;
        for (std::size_t i = 0; i < res.candidates.size(); ++i) {
            BlockCandidate& c = res.candidates[i];
            if (c.multiset != content || c.panel_log_prior != best_prior || c.likelihood < best_ll - tie_tolerance_nats) continue;
            c.allocation_tied = true;
            tied.insert(unordered(c.ctx1_sequence, c.ctx2_sequence));
            const auto key = c.ctx1_id <= c.ctx2_id ? std::make_pair(c.ctx1_id, c.ctx2_id) : std::make_pair(c.ctx2_id, c.ctx1_id);
            if (first || key < smallest) { smallest = key; pick = i; first = false; }
        }
        res.allocation_tie = tied.size();
        res.representative_only = tied.size() > 1;
    }
    res.chosen = static_cast<long>(pick);
    BlockCandidate& chosen = res.candidates[pick];
    chosen.chosen = true;
    res.repair_separation = chosen.likelihood - res.original_likelihood;
    const Genotype cg = unordered(chosen.ctx1_sequence, chosen.ctx2_sequence);
    double other = -std::numeric_limits<double>::infinity();
    for (const BlockCandidate& c : res.candidates) {
        if (unordered(c.ctx1_sequence, c.ctx2_sequence) != cg) other = std::max(other, c.likelihood);
    }
    res.gap_defined = std::isfinite(other);
    res.conditional_repair_gap = res.gap_defined ? chosen.likelihood - other : 0.0;
}


RepairResult repair_block(const BlockAudit& audit, const CensusBlock& census, const BlockFeatureModel& fm,
                          const NodeSequences& nodes, const ScoringInputs& in, const DosageSettings& settings,
                          const RepairSettings& rs) {
    RepairResult res;
    res.route = census_route_name(census.route);
    const PreparedGenotypeModel& model = *in.model;
    const ReportedCall& call = *in.call;
    const std::size_t b = audit.block;
    const PreparedGenotypeBlock& block = model.blocks()[b];
    const auto t_construct = std::chrono::steady_clock::now();
    std::vector<std::size_t> findex;
    const std::vector<std::vector<Kmer>> tested = tested_feature_kmers(audit, fm, call, settings, findex);
    const auto refuse = [&](const std::string& reason) {
        res.outcome = "TRIGGERED_UNREPAIRED";
        res.reason = reason;
        classify_repair_scope(res, census.route);
        return res;
    };

    if (census.route == CensusRoute::Pair) {
        std::size_t pairs = 0;
        for (std::size_t a = 0; a < block.alleles.size(); ++a) {
            for (std::size_t c = a; c < block.alleles.size(); ++c) {
                bool ok = true;
                for (std::size_t t = 0; t < tested.size() && ok; ++t) {
                    const PooledFeature& pf = fm.features[findex[t]];
                    ok = admissible(audit.features[t].admissible, static_cast<long>(pf.copies[a]) + static_cast<long>(pf.copies[c]));
                }
                if (!ok) continue;
                ++pairs;
                for (int o = 0; o < (a == c ? 1 : 2); ++o) {
                    BlockCandidate bc;
                    bc.label = "pair_" + std::to_string(a) + "_" + std::to_string(c);
                    bc.orientation = o == 0 ? "AB" : "BA";
                    bc.allele_a = a; bc.allele_b = c;
                    bc.rank = pairs - 1;
                    const std::size_t x = o == 0 ? a : c, y = o == 0 ? c : a;
                    bc.ctx1_sequence = block.alleles[x].sequence; bc.ctx2_sequence = block.alleles[y].sequence;
                    bc.ctx1_id = std::to_string(x); bc.ctx2_id = std::to_string(y);
                    res.candidates.push_back(std::move(bc));
                }
            }
        }
        res.diploid = pairs;
        if (pairs > rs.pair_cap) { res.candidates.clear(); res.cap_reached = "pairs"; return refuse("cap_pairs"); }
        if (pairs == 0) return refuse("no_admissible_candidate");
    } else if (census.route == CensusRoute::RepeatedUnit) {
        const UnitLibrary lib = build_unit_library(block, nodes, rs);
        const std::size_t K = settings.kmer;
        const std::size_t nm = lib.units.size();
        // module k-mer ids (first occurrence over units in id order), occurrences with positions
        std::unordered_map<Kmer, std::uint32_t> kid;
        std::vector<std::vector<std::pair<std::uint32_t, std::size_t>>> occ(nm);
        for (std::size_t u = 0; u < nm; ++u) {
            for (const KmerAt& x : canonical_kmers(lib.units[u].sequence, K)) {
                const auto [it, added] = kid.emplace(x.kmer, static_cast<std::uint32_t>(kid.size()));
                occ[u].push_back({it->second, x.position});
            }
        }
        const std::size_t nk = kid.size();
        // unit counts of each tested feature: the median over its k-mers of the k-mer's count in the unit
        std::vector<std::vector<double>> fmat(nm, std::vector<double>(tested.size()));
        for (std::size_t u = 0; u < nm; ++u) {
            std::unordered_map<std::uint32_t, std::uint32_t> cnt;
            for (const auto& [id, pos] : occ[u]) ++cnt[id];
            for (std::size_t t = 0; t < tested.size(); ++t) {
                std::vector<std::uint32_t> vals;
                vals.reserve(tested[t].size());
                for (Kmer x : tested[t]) {
                    const auto f = kid.find(x);
                    if (f == kid.end()) { vals.push_back(0); continue; }
                    const auto c = cnt.find(f->second);
                    vals.push_back(c == cnt.end() ? 0 : c->second);
                }
                std::sort(vals.begin(), vals.end());
                const std::size_t n = vals.size();
                fmat[u][t] = n % 2 ? static_cast<double>(vals[n / 2]) : (static_cast<double>(vals[n / 2 - 1]) + static_cast<double>(vals[n / 2])) / 2.0;
            }
        }
        std::vector<std::size_t> term, inte;
        for (const Unit& u : lib.units) { if (u.terminal_capable) term.push_back(u.id); if (u.internal_capable) inte.push_back(u.id); }
        res.internal_capable = inte.size(); res.terminal_capable = term.size();
        const std::size_t kmax = lib.kmax;
        for (std::size_t M = 2; M <= 2 * kmax; ++M) {
            bool ok = true;
            for (std::size_t t = 0; t < tested.size() && ok; ++t) {
                double tmin = std::numeric_limits<double>::infinity(), tmax = -tmin, imin = tmin, imax = -tmin;
                for (std::size_t u : term) { tmin = std::min(tmin, fmat[u][t]); tmax = std::max(tmax, fmat[u][t]); }
                for (std::size_t u : inte) { imin = std::min(imin, fmat[u][t]); imax = std::max(imax, fmat[u][t]); }
                const double lo = 2 * tmin + (M > 2 ? static_cast<double>(M - 2) * imin : 0.0);
                const double hi = 2 * tmax + (M > 2 ? static_cast<double>(M - 2) * imax : 0.0);
                bool any = false;
                for (int a : audit.features[t].admissible) any |= std::floor(lo) <= a && a <= std::ceil(hi);
                ok = any;
            }
            if (ok) res.m_range.push_back(static_cast<int>(M));
        }
        if (res.m_range.empty()) return refuse("no_admissible_candidate");

        // junction k-mers
        std::unordered_set<Kmer> junction;
        {
            std::set<std::string> starts;
            for (const Unit& u : lib.units) starts.insert(u.sequence.substr(0, K - 1));
            std::vector<std::string> joins;
            for (const std::string& st : starts) {
                joins.push_back(lib.flank_left.substr(lib.flank_left.size() - (K - 1)) + st);
                for (const Unit& u : lib.units) if (u.internal_capable) joins.push_back(u.sequence.substr(u.sequence.size() - (K - 1)) + st);
            }
            for (const Unit& u : lib.units) if (u.terminal_capable) joins.push_back(u.sequence.substr(u.sequence.size() - (K - 1)) + lib.flank_right.substr(0, K - 1));
            for (const std::string& j : joins) for (const KmerAt& x : canonical_kmers(j, K)) junction.insert(x.kmer);
        }
        // outside k-mers of this call
        std::unordered_set<Kmer> outside;
        for (int h = 0; h < 2; ++h) {
            const SequenceInterval iv = call.intervals[h][b];
            for (const KmerAt& x : canonical_kmers(call.haplotypes[h], K)) if (x.position + K <= iv.begin || x.position >= iv.end) outside.insert(x.kmer);
        }
        std::vector<char> excl(nk, 0);
        std::vector<Kmer> by_id(nk);
        for (const auto& [x, id] : kid) { by_id[id] = x; if (outside.count(x) || junction.count(x)) excl[id] = 1; }
        // profile classes: discriminating k-mers grouped by their count vector across units, classes in lexicographic order
        std::vector<std::vector<std::uint16_t>> col(nk, std::vector<std::uint16_t>(nm, 0));
        for (std::size_t u = 0; u < nm; ++u) for (const auto& [id, pos] : occ[u]) ++col[id][u];
        std::map<std::vector<std::uint16_t>, std::vector<std::uint32_t>> classes;
        for (std::size_t id = 0; id < nk; ++id) {
            if (excl[id] || all_same(col[id])) continue;
            classes[col[id]].push_back(static_cast<std::uint32_t>(id));
        }
        std::vector<long> kclass(nk, -1);
        std::size_t ncls = 0;
        for (const auto& [vec, ids] : classes) { for (std::uint32_t id : ids) kclass[id] = static_cast<long>(ncls); ++ncls; }
        std::size_t read_len = 0;
        for (const ReadPair& p : *in.reads) read_len = std::max(read_len, p.first.bases.size());
        const double w = static_cast<double>(read_len - K + 1);
        std::vector<std::vector<double>> U(ncls, std::vector<double>(nm, 0.0));
        for (std::size_t u = 0; u < nm; ++u) {
            std::vector<std::pair<long, std::size_t>> cp;
            for (const auto& [id, pos] : occ[u]) if (kclass[id] >= 0) cp.push_back({kclass[id], pos});
            std::sort(cp.begin(), cp.end());
            for (std::size_t i = 0; i < cp.size(); ++i) {
                const double gap = i == 0 || cp[i].first != cp[i - 1].first ? w : std::min(w, static_cast<double>(cp[i].second - cp[i - 1].second));
                U[static_cast<std::size_t>(cp[i].first)][u] += gap;
            }
        }
        // reference k-mers for beta: the context blocks adjacent to this block, on both homologues
        if (b == 0 || b + 1 >= model.blocks().size() || model.blocks()[b - 1].kind != GenotypeBlockKind::Context ||
            model.blocks()[b + 1].kind != GenotypeBlockKind::Context) {
            return refuse("no_reference");
        }
        std::array<std::unordered_map<Kmer, std::uint32_t>, 2> hc;
        std::array<std::unordered_set<Kmer>, 2> fl;
        for (int h = 0; h < 2; ++h) {
            for (const KmerAt& x : canonical_kmers(call.haplotypes[h], K)) {
                ++hc[h][x.kmer];
                for (std::size_t nb : {b - 1, b + 1}) {
                    const SequenceInterval iv = call.intervals[h][nb];
                    if (x.position >= iv.begin && x.position + K <= iv.end) fl[h].insert(x.kmer);
                }
            }
        }
        std::vector<Kmer> ref;
        for (Kmer x : fl[0]) if (fl[1].count(x) && hc[0][x] == 1 && hc[1][x] == 1 && !kid.count(x)) ref.push_back(x);
        std::sort(ref.begin(), ref.end());
        std::unordered_map<Kmer, std::size_t> ref_index;
        for (std::size_t i = 0; i < ref.size(); ++i) ref_index[ref[i]] = i;
        // observed fragments per class, and per reference k-mer
        std::vector<double> y(ncls, 0.0);
        std::vector<double> yref(ref.size(), 0.0);
        for (const ReadPair& p : *in.reads) {
            std::set<long> hit;
            std::set<std::size_t> rh;
            for (const std::string* s : {&p.first.bases, &p.second.bases}) {
                for (const KmerAt& x : canonical_kmers(*s, K)) {
                    const auto f = kid.find(x.kmer);
                    if (f != kid.end()) { if (kclass[f->second] >= 0) hit.insert(kclass[f->second]); }
                    else { const auto r = ref_index.find(x.kmer); if (r != ref_index.end()) rh.insert(r->second); }
                }
            }
            for (long c : hit) y[static_cast<std::size_t>(c)] += 1.0;
            for (std::size_t r : rh) yref[r] += 1.0;
        }
        if (yref.size() < 2) return refuse("no_reference");
        const double mu = numpy_sum(yref.data(), yref.size()) / static_cast<double>(yref.size());
        std::vector<double> dev(yref.size());
        for (std::size_t i = 0; i < yref.size(); ++i) dev[i] = (yref[i] - mu) * (yref[i] - mu);
        const double var = numpy_sum(dev.data(), dev.size()) / static_cast<double>(yref.size() - 1);
        if (!(mu > 0.0)) return refuse("no_reference_hits");
        const double beta = mu / (2.0 * w);
        const double phi = var > mu ? mu * mu / (var - mu) : std::numeric_limits<double>::infinity();
        const double e_k = static_cast<double>(K) * in.library->substitution_error / 3.0;
        res.beta = beta; res.phi = phi;
        std::vector<double> Umax(ncls, 0.0);
        for (std::size_t c = 0; c < ncls; ++c) Umax[c] = *std::max_element(U[c].begin(), U[c].end());
        // screening: Poisson EM over unit weights. Sparse over U's nonzero entries, visited in the dense loops' index
        // order: a skipped term is an exact zero, so every sum -- and so every weight -- is bit-identical to the dense form.
        const auto t_em = std::chrono::steady_clock::now();
        const double M0 = static_cast<double>(res.m_range.front());
        std::vector<double> wt(nm, M0 / static_cast<double>(nm)), wn(nm), denom(nm, 0.0), eps(ncls), lam(ncls), g(nm), ratio(ncls);
        std::vector<std::vector<std::pair<std::size_t, double>>> by_class(ncls), by_unit(nm);
        for (std::size_t c = 0; c < ncls; ++c) for (std::size_t u = 0; u < nm; ++u) if (U[c][u] != 0.0) by_class[c].push_back({u, U[c][u]});
        for (std::size_t u = 0; u < nm; ++u) for (std::size_t c = 0; c < ncls; ++c) if (U[c][u] != 0.0) by_unit[u].push_back({c, U[c][u]});
        for (std::size_t c = 0; c < ncls; ++c) eps[c] = e_k * beta * Umax[c] * M0;
        for (std::size_t u = 0; u < nm; ++u) { double s = 0.0; for (const auto& [c, x] : by_unit[u]) s += x; denom[u] = beta * s; }
        std::size_t it = 0;
        for (; it < rs.em_iterations; ++it) {
            for (std::size_t c = 0; c < ncls; ++c) { double s = 0.0; for (const auto& [u, x] : by_class[c]) s += x * wt[u]; lam[c] = beta * s + eps[c]; }
            for (std::size_t c = 0; c < ncls; ++c) ratio[c] = y[c] / lam[c];
            for (std::size_t u = 0; u < nm; ++u) { double s = 0.0; for (const auto& [c, x] : by_unit[u]) s += x * ratio[c]; g[u] = s; }
            double mx = 0.0;
            for (std::size_t u = 0; u < nm; ++u) { wn[u] = wt[u] * (beta * g[u]) / (denom[u] > 0 ? denom[u] : 1.0); mx = std::max(mx, std::fabs(wn[u] - wt[u])); }
            wt = wn;
            if (mx < rs.em_tolerance) break;
        }
        res.seconds_em = seconds_since(t_em);
        res.em_iterations = std::min(it + 1, rs.em_iterations);
        std::vector<std::size_t> rank_int = inte, rank_term = term;
        const auto by_weight = [&](std::size_t a, std::size_t c) { return wt[a] != wt[c] ? wt[a] > wt[c] : a < c; };
        std::sort(rank_int.begin(), rank_int.end(), by_weight);
        std::sort(rank_term.begin(), rank_term.end(), by_weight);
        const std::vector<std::size_t> pool_i(rank_int.begin(), rank_int.begin() + static_cast<long>(std::min(rs.unit_pool, rank_int.size())));
        const std::vector<std::size_t> pool_t(rank_term.begin(), rank_term.begin() + static_cast<long>(std::min(rs.unit_pool, rank_term.size())));
        // exhaustive admissible multisets over the pool, as compact fixed-width records (tens of millions of them at 30x:
        // three heap vectors each exceeded the footprint budget). Order and comparisons are the vector forms' exactly.
        const auto t_enum = std::chrono::steady_clock::now();
        constexpr std::size_t MAXM = 8;
        if (2 * kmax > MAXM || nm > 65535) throw std::runtime_error("repair: unit library exceeds the compact candidate record");
        struct Cand {
            std::array<std::uint16_t, 2> tp{};
            std::array<std::uint16_t, MAXM> ip{}, ms{};
            std::uint8_t ni = 0, nms = 0;
            std::vector<std::size_t> tp_vec() const { return {tp[0], tp[1]}; }
            std::vector<std::size_t> ip_vec() const { return std::vector<std::size_t>(ip.begin(), ip.begin() + ni); }
            std::vector<std::size_t> ms_vec() const { return std::vector<std::size_t>(ms.begin(), ms.begin() + nms); }
        };
        const auto ms_less = [](const Cand& a, const Cand& b) {   // std::vector<std::size_t>::operator<
            return std::lexicographical_compare(a.ms.begin(), a.ms.begin() + a.nms, b.ms.begin(), b.ms.begin() + b.nms);
        };
        std::vector<Cand> cand;
        const std::size_t nf = tested.size();
        for (int Mv : res.m_range) {
            if (Mv < 2 || static_cast<std::size_t>(Mv - 2) > 2 * (kmax - 1)) continue;
            std::vector<std::vector<std::size_t>> combos;
            std::vector<std::vector<double>> Fi;
            combinations_with_replacement(pool_i.size(), static_cast<std::size_t>(Mv - 2), [&](const std::vector<std::size_t>& idx) {
                std::vector<std::size_t> c;
                std::vector<double> f(nf, 0.0);
                for (std::size_t i : idx) { c.push_back(pool_i[i]); for (std::size_t t = 0; t < nf; ++t) f[t] += fmat[pool_i[i]][t]; }
                combos.push_back(std::move(c));
                Fi.push_back(std::move(f));
            });
            combinations_with_replacement(pool_t.size(), 2, [&](const std::vector<std::size_t>& tidx) {
                const std::size_t t0 = pool_t[tidx[0]], t1 = pool_t[tidx[1]];
                for (std::size_t r = 0; r < combos.size(); ++r) {
                    bool ok = true;
                    for (std::size_t t = 0; t < nf && ok; ++t) {
                        const double tot = std::nearbyint((Fi[r][t] + fmat[t0][t]) + fmat[t1][t]);
                        ok = admissible(audit.features[t].admissible, static_cast<long>(tot));
                    }
                    if (!ok) continue;
                    Cand c;
                    c.tp = {static_cast<std::uint16_t>(t0), static_cast<std::uint16_t>(t1)};
                    c.ni = static_cast<std::uint8_t>(combos[r].size());
                    for (std::size_t i = 0; i < combos[r].size(); ++i) c.ip[i] = static_cast<std::uint16_t>(combos[r][i]);
                    c.nms = static_cast<std::uint8_t>(2 + combos[r].size());
                    c.ms[0] = c.tp[0]; c.ms[1] = c.tp[1];
                    for (std::size_t i = 0; i < c.ni; ++i) c.ms[2 + i] = c.ip[i];
                    std::sort(c.ms.begin(), c.ms.begin() + c.nms);
                    cand.push_back(c);
                }
            });
        }
        res.multisets_evaluated = cand.size();
        res.seconds_enumerate = seconds_since(t_enum);
        if (cand.empty()) return refuse("no_admissible_candidate");
        // composite NB log-likelihood (prototype loglik_matrix), classes varying across the pool summed per candidate;
        // candidates are independent, so they are evaluated in parallel, each into its own slot
        const auto t_ll = std::chrono::steady_clock::now();
        std::vector<std::size_t> cols;
        { std::set<std::size_t> s(pool_i.begin(), pool_i.end()); s.insert(pool_t.begin(), pool_t.end()); cols.assign(s.begin(), s.end()); }
        std::vector<std::size_t> vary, fixed;
        for (std::size_t c = 0; c < ncls; ++c) {
            double mn = std::numeric_limits<double>::infinity(), mx = -mn;
            for (std::size_t u : cols) { mn = std::min(mn, U[c][u]); mx = std::max(mx, U[c][u]); }
            (mx != mn ? vary : fixed).push_back(c);
        }
        const auto loglik = [&](const std::vector<std::size_t>& cls, const std::vector<std::pair<std::size_t, double>>& n, std::vector<double>& term) {
            double Mtot = 0.0;
            for (const auto& [u, k] : n) Mtot += k;
            term.resize(cls.size());
            for (std::size_t i = 0; i < cls.size(); ++i) {
                const std::size_t c = cls[i];
                double dot = 0.0;
                for (const auto& [u, k] : n) dot += k * U[c][u];
                const double l = beta * dot + e_k * beta * Umax[c] * Mtot;
                term[i] = std::isinf(phi) ? y[c] * std::log(l) - l : y[c] * std::log(l) - (y[c] + phi) * std::log(l + phi);
            }
            double s = numpy_sum(term.data(), term.size());
            if (!std::isinf(phi)) s += (phi * std::log(phi)) * static_cast<double>(cls.size());
            return s;
        };
        std::map<int, double> const_M;
        {
            std::vector<double> term;
            for (int Mv : res.m_range) const_M[Mv] = fixed.empty() ? 0.0 : loglik(fixed, {{cols.front(), static_cast<double>(Mv)}}, term);
        }
        std::vector<double> ll(cand.size());
        const std::size_t ll_chunks = std::max<std::size_t>(1, std::min<std::size_t>(cand.size(), 256));
        run_parallel(ll_chunks, in.options.threads, [&](std::size_t ci) {
            std::vector<double> term;
            std::vector<std::pair<std::size_t, double>> n;
            for (std::size_t r = cand.size() * ci / ll_chunks; r < cand.size() * (ci + 1) / ll_chunks; ++r) {
                // unit counts in ascending unit order (cols and ms are both ascending), as the dense column scan gives
                n.clear();
                const Cand& cd = cand[r];
                for (std::size_t i = 0; i < cd.nms;) {
                    std::size_t j = i;
                    while (j < cd.nms && cd.ms[j] == cd.ms[i]) ++j;
                    n.push_back({cd.ms[i], static_cast<double>(j - i)});
                    i = j;
                }
                ll[r] = loglik(vary, n, term) + const_M.at(static_cast<int>(cd.nms));
            }
        });
        std::vector<std::uint32_t> order(cand.size());
        for (std::size_t i = 0; i < order.size(); ++i) order[i] = static_cast<std::uint32_t>(i);
        std::stable_sort(order.begin(), order.end(), [&](std::uint32_t a, std::uint32_t c) {
            return ll[a] != ll[c] ? ll[a] > ll[c] : ms_less(cand[a], cand[c]);
        });
        res.seconds_loglik = seconds_since(t_ll);
        // bounded structural construction, k-best order, stop at the first cap
        std::set<std::vector<std::size_t>> alleles;
        std::map<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>, std::pair<std::size_t, std::size_t>> diplos;  // key -> (multiset rank, cand)
        std::size_t retained = 0;
        bool refused = false;
        for (std::size_t pos = 0; pos < std::min(order.size(), rs.cap_multisets + 1); ++pos) {
            const std::size_t r = order[pos];
            if (pos == rs.cap_multisets) {
                res.cap_reached = "multisets";
                refused = ll[order[pos - 1]] - ll[r] == 0.0;
                break;
            }
            const Cand& cd = cand[r];
            const std::vector<std::size_t> cd_tp = cd.tp_vec(), cd_ip = cd.ip_vec();
            std::set<std::vector<std::size_t>> new_a;
            std::map<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>, std::pair<std::size_t, std::size_t>> new_d;
            std::set<std::pair<std::size_t, std::size_t>> term_orders = {{cd_tp[0], cd_tp[1]}, {cd_tp[1], cd_tp[0]}};
            for (const auto& [ta, tb] : term_orders) {
                const std::vector<std::size_t>& ipl = cd_ip;
                for (std::size_t mask = 0; mask < (std::size_t{1} << ipl.size()); ++mask) {
                    std::vector<std::size_t> A, Bm;
                    for (std::size_t i = 0; i < ipl.size(); ++i) ((mask >> i) & 1 ? A : Bm).push_back(ipl[i]);
                    if (A.size() + 1 > kmax || Bm.size() + 1 > kmax) continue;
                    std::sort(A.begin(), A.end()); std::sort(Bm.begin(), Bm.end());
                    do {
                        std::vector<std::size_t> Bp = Bm;
                        do {
                            std::vector<std::size_t> ha = A, hb = Bp;
                            ha.push_back(ta); hb.push_back(tb);
                            if (!(lib.units[ha[0]].anchor_compatible && lib.units[hb[0]].anchor_compatible)) continue;
                            const auto key = ha <= hb ? std::make_pair(ha, hb) : std::make_pair(hb, ha);
                            for (const auto* h : {&ha, &hb}) if (!alleles.count(*h)) new_a.insert(*h);
                            if (!diplos.count(key) && !new_d.count(key)) new_d[key] = {pos + 1, r};
                        } while (std::next_permutation(Bp.begin(), Bp.end()));
                    } while (std::next_permutation(A.begin(), A.end()));
                }
            }
            if (alleles.size() + new_a.size() > rs.cap_haploid || diplos.size() + new_d.size() > rs.cap_diploid) {
                res.cap_reached = alleles.size() + new_a.size() > rs.cap_haploid ? "haploid" : "diploid";
                refused = pos == 0 || ll[order[pos - 1]] - ll[r] == 0.0;
                break;
            }
            alleles.insert(new_a.begin(), new_a.end());
            diplos.insert(new_d.begin(), new_d.end());
            ++retained;
        }
        res.multisets_retained = retained;
        res.haploid = alleles.size();
        res.diploid = diplos.size();
        std::vector<std::string> why;
        if (inte.size() > rs.unit_pool || term.size() > rs.unit_pool) why.push_back("pool_truncated");
        if (res.cap_reached != "none") why.push_back("cap_" + res.cap_reached);
        if (retained < cand.size()) why.push_back("proposal_ranked");
        if (!why.empty()) {
            res.scope = "CONDITIONAL";
            res.scope_reason.clear();
            for (const std::string& x : why) res.scope_reason += (res.scope_reason.empty() ? "" : ",") + x;
        }
        {
            // exact sizes of the candidate space, without scoring (see CandidateSpace); always counted, because the
            // dosage-content scope of the repair rests on it
            const auto t_space = std::chrono::steady_clock::now();
            CandidateSpace& sp = res.space;
            sp.internal_capable = inte.size(); sp.terminal_capable = term.size();
            sp.pool_internal = pool_i.size(); sp.pool_terminal = pool_t.size();
            sp.pool_multisets = cand.size();
            sp.retained_multisets = retained; sp.constructed_haploid = alleles.size(); sp.constructed_diploid = diplos.size();
            sp.cap_reached = res.cap_reached;
            // feature sums in half-copies (unit counts are medians: whole or half)
            const auto halves = [&](std::size_t u, std::size_t t) { return static_cast<long>(std::llround(2.0 * fmat[u][t])); };
            const auto admissible_total = [&](const std::vector<long>& h2) {   // h2: summed half-copies per feature
                for (std::size_t t = 0; t < nf; ++t) {
                    if (!admissible(audit.features[t].admissible, static_cast<long>(std::nearbyint(static_cast<double>(h2[t]) / 2.0)))) return false;
                }
                return true;
            };
            std::set<int> mset(res.m_range.begin(), res.m_range.end());
            // (terminal pair, internal multiset) counts by DP over units: the prototype's candidate unit
            const auto count_multisets = [&](const std::vector<std::size_t>& I, const std::vector<std::size_t>& T) {
                const std::size_t kmax_int = 2 * (kmax - 1);
                std::vector<std::map<std::vector<long>, std::uint64_t>> dp(kmax_int + 1);
                dp[0][std::vector<long>(nf, 0)] = 1;
                for (std::size_t u : I) {   // multisets: add any number of copies of u, units in a fixed order
                    std::vector<std::map<std::vector<long>, std::uint64_t>> nx = dp;
                    for (std::size_t k = 0; k <= kmax_int; ++k) {
                        for (const auto& [sv, c] : dp[k]) {
                            std::vector<long> v = sv;
                            for (std::size_t add = 1; k + add <= kmax_int; ++add) {
                                for (std::size_t t = 0; t < nf; ++t) v[t] += halves(u, t);
                                nx[k + add][v] += c;
                            }
                        }
                    }
                    dp = std::move(nx);
                }
                std::map<std::vector<long>, std::uint64_t> tp;
                for (std::size_t a = 0; a < T.size(); ++a) for (std::size_t b = a; b < T.size(); ++b) {
                    std::vector<long> v(nf);
                    for (std::size_t t = 0; t < nf; ++t) v[t] = halves(T[a], t) + halves(T[b], t);
                    ++tp[v];
                }
                std::uint64_t total = 0;
                for (std::size_t k = 0; k <= kmax_int; ++k) {
                    if (!mset.count(static_cast<int>(k + 2))) continue;
                    for (const auto& [si, ci] : dp[k]) for (const auto& [st, ct] : tp) {
                        std::vector<long> v(nf);
                        for (std::size_t t = 0; t < nf; ++t) v[t] = si[t] + st[t];
                        if (admissible_total(v)) total += ci * ct;
                    }
                }
                return total;
            };
            sp.pool_multisets_dp = count_multisets(pool_i, pool_t);
            sp.full_multisets = count_multisets(inte, term);
            // every valid haploid tuple of the full library: (internal..., terminal), <= KMAX units, the first
            // anchor-compatible; histogrammed by (feature half-sums, units, I-state units)
            std::map<std::tuple<std::vector<long>, int, int>, std::uint64_t> hist;
            std::vector<std::size_t> first_int;
            for (std::size_t u : inte) if (lib.units[u].anchor_compatible) first_int.push_back(u);
            const auto add_haploid = [&](const std::vector<std::size_t>& tuple) {
                std::vector<long> v(nf, 0);
                int nI = 0;
                for (std::size_t u : tuple) { for (std::size_t t = 0; t < nf; ++t) v[t] += halves(u, t); nI += lib.units[u].state == 'I'; }
                ++hist[{v, static_cast<int>(tuple.size()), nI}];
                ++sp.full_haploid_valid;
            };
            for (std::size_t tt : term) {
                if (lib.units[tt].anchor_compatible) add_haploid({tt});
                if (kmax >= 2) for (std::size_t a : first_int) {
                    add_haploid({a, tt});
                    if (kmax >= 3) for (std::size_t b : inte) add_haploid({a, b, tt});
                    if (kmax >= 4) throw std::runtime_error("space audit: KMAX > 3 is not enumerated");
                }
            }
            std::vector<std::pair<std::tuple<std::vector<long>, int, int>, std::uint64_t>> H(hist.begin(), hist.end());
            std::vector<char> partnered(H.size(), 0);
            std::set<std::string> vectors, states;
            for (std::size_t a = 0; a < H.size(); ++a) {
                for (std::size_t b = a; b < H.size(); ++b) {
                    const auto& [ka, ca] = H[a];
                    const auto& [kb, cb] = H[b];
                    if (!mset.count(std::get<1>(ka) + std::get<1>(kb))) continue;
                    std::vector<long> v(nf);
                    for (std::size_t t = 0; t < nf; ++t) v[t] = std::get<0>(ka)[t] + std::get<0>(kb)[t];
                    if (!admissible_total(v)) continue;
                    sp.full_diploid += a == b ? ca * (ca + 1) / 2 : ca * cb;
                    partnered[a] = partnered[b] = 1;
                    std::string vs;
                    for (std::size_t t = 0; t < nf; ++t) vs += (t ? "," : "") + std::to_string(static_cast<long>(std::nearbyint(static_cast<double>(v[t]) / 2.0)));
                    vectors.insert(vs);
                    states.insert(std::to_string(std::get<1>(ka) + std::get<1>(kb)) + "u/" + std::to_string(std::get<2>(ka) + std::get<2>(kb)) + "I");
                }
            }
            for (std::size_t a = 0; a < H.size(); ++a) if (partnered[a]) sp.full_haploid_in_admissible += H[a].second;
            sp.dosage_vectors.assign(vectors.begin(), vectors.end());
            sp.module_state_counts.assign(states.begin(), states.end());
            sp.seconds = seconds_since(t_space);
            // the content every admissible candidate of the full library carries, when it is unique
            res.content_complete = sp.dosage_vectors.size() == 1 && sp.module_state_counts.size() == 1;
            if (rs.space_audit) {
                res.outcome = "SPACE_AUDIT";
                res.candidates.clear();
                return res;
            }
        }
        if (refused) return refuse("cap_refused");
        // structures in (multiset rank, key) order; one orientation (the prototype's: the smaller unit tuple first)
        std::vector<std::pair<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>, std::pair<std::size_t, std::size_t>>> st(diplos.begin(), diplos.end());
        std::stable_sort(st.begin(), st.end(), [](const auto& x, const auto& z) {
            return x.second.first != z.second.first ? x.second.first < z.second.first : x.first < z.first;
        });
        std::map<std::string, std::size_t> freq;
        for (const auto& t : lib.panel) { std::string s; for (std::size_t u : t) s += lib.units[u].state; ++freq[s]; }
        const double N = static_cast<double>(lib.panel.size());
        std::set<std::string> spelled;
        for (std::size_t n = 0; n < st.size(); ++n) {
            const auto& [key, meta] = st[n];
            BlockCandidate bc;
            bc.label = "s" + std::to_string(n);
            bc.rank = n;
            bc.units_a = key.first; bc.units_b = key.second;
            bc.multiset = cand[meta.second].ms_vec();
            bc.multiset_rank = meta.first;
            bc.proposal_loglik = ll[order[meta.first - 1]];
            std::string sa, sb;
            for (std::size_t u : key.first) sa += lib.units[u].state;
            for (std::size_t u : key.second) sb += lib.units[u].state;
            bc.states = sa + "|" + sb;
            bc.panel_log_prior = std::log((static_cast<double>(freq[sa]) + 0.5) / N) + std::log((static_cast<double>(freq[sb]) + 0.5) / N);
            bc.ctx1_sequence = lib.block_sequence(key.first);
            bc.ctx2_sequence = lib.block_sequence(key.second);
            for (const std::string* s : {&bc.ctx1_sequence, &bc.ctx2_sequence}) spelled.insert(*s);
            bc.ctx1_id = sha256_hex(bc.ctx1_sequence); bc.ctx2_id = sha256_hex(bc.ctx2_sequence);
            res.candidates.push_back(std::move(bc));
        }
        std::set<std::vector<std::size_t>> tuples;
        for (const auto& [key, meta] : st) { tuples.insert(key.first); tuples.insert(key.second); }
        if (spelled.size() != tuples.size()) throw std::runtime_error("repair: two unit tuples spell one block sequence");
    } else {
        throw std::runtime_error("repair_block on a block the census does not mark sensitive");
    }
    res.seconds_construct = seconds_since(t_construct);

    // exact scores
    const auto t_score = std::chrono::steady_clock::now();
    std::vector<std::pair<std::string, std::string>> ctx;
    for (const BlockCandidate& c : res.candidates) ctx.push_back({c.ctx1_sequence, c.ctx2_sequence});
    ExactScore original;
    std::vector<ExactScore> scores;
    try {
        scores = rs.lifecycle == "cold" ? score_block_candidates(in, b, ctx, rs, original)
                                          : score_block_candidates_context(in, b, ctx, rs, original, nullptr, &res.scoring);
        if (rs.lifecycle == "cold") res.scoring.lifecycle = "cold";
    } catch (const GuardFailure&) {
        res.candidates.clear();   // no partial output
        return refuse("winner_not_reproduced");
    } catch (const ResourceExhausted&) {
        res.candidates.clear();   // no partial output: a refused repair reports no candidate or score
        return refuse("resource_refusal");
    }
    res.seconds_score = seconds_since(t_score);
    res.original_fragment = original.fragment; res.original_exposure = original.exposure; res.original_likelihood = original.likelihood;
    for (std::size_t i = 0; i < scores.size(); ++i) {
        res.candidates[i].fragment_term = scores[i].fragment;
        res.candidates[i].exposure_term = scores[i].exposure;
        res.candidates[i].likelihood = scores[i].likelihood;
    }

    choose_repair(res, census.route, rs.tie_tolerance_nats);
    classify_repair_scope(res, census.route);
    return res;
}

// ---------------------------------------------------------------------------------------------- allocation test

AllocationTestResult allocation_observability_test(const UnitLibrary& lib, const std::vector<std::size_t>& hap1_units,
                                                   const std::vector<std::size_t>& hap2_units, std::size_t block,
                                                   const ScoringInputs& in, const RepairSettings& rs) {
    const auto started = std::chrono::steady_clock::now();
    AllocationTestResult res;
    struct Copy { std::string seq; bool terminal; };
    std::vector<Copy> copies;
    for (const auto* h : {&hap1_units, &hap2_units}) {
        for (std::size_t i = 0; i < h->size(); ++i) copies.push_back({lib.units.at((*h)[i]).sequence, i + 1 == h->size()});
    }
    res.copies = copies.size();
    if (copies.size() < 3 || hap1_units.empty() || hap2_units.empty()) { res.seconds = seconds_since(started); return res; }
    const std::size_t la = hap1_units.size(), M = copies.size();
    std::vector<std::size_t> A0, A1;
    for (std::size_t i = 0; i < M; ++i) (i < la ? A0 : A1).push_back(i);
    // the anchor: the unit bases inside the left anchor trim
    const long lock = static_cast<long>(lib.trim_left) - static_cast<long>(lib.flank_left.size());
    std::string anchor;
    for (const Unit& u : lib.units) if (u.anchor_compatible) { anchor = u.sequence.substr(0, static_cast<std::size_t>(std::max(0L, lock))); break; }
    // a common switch window: 2 x 600 bp of copy 0, centred at 250..20000 step 25, occurring exactly once in every copy
    const std::size_t FL = 600;
    std::vector<std::size_t> offs;
    {
        const std::string& ref = copies[0].seq;
        const std::size_t hi = std::min<std::size_t>(20000, ref.size() > FL ? ref.size() - FL : 0);
        for (std::size_t p = std::max<std::size_t>(250, FL); p < hi; p += 25) {
            const std::string w = ref.substr(p - FL, 2 * FL);
            std::vector<std::size_t> o;
            for (const Copy& c : copies) {
                const std::size_t i = c.seq.find(w);
                if (i == std::string::npos || c.seq.find(w, i + 1) != std::string::npos) break;
                o.push_back(i + FL);
            }
            if (o.size() == M) { offs = std::move(o); break; }
        }
    }
    if (offs.empty()) { res.status = "INCONCLUSIVE"; res.seconds = seconds_since(started); return res; }
    res.switch_coordinate = static_cast<long>(offs[0]);
    std::vector<std::string> start(M), body(M);
    for (std::size_t i = 0; i < M; ++i) { start[i] = copies[i].seq.substr(0, offs[i]); body[i] = copies[i].seq.substr(offs[i]); }
    std::map<std::size_t, std::string> nxt;
    const std::array<std::string, 2> flank_starts = {start[A0.front()], start[A1.front()]};
    for (const auto* h : {&A0, &A1}) for (std::size_t j = 0; j + 1 < h->size(); ++j) nxt[(*h)[j]] = start[(*h)[j + 1]];
    // alternative allocations: the construction rules (one terminal copy last per homologue, <= KMAX per homologue,
    // anchor-compatible first copies), deduplicated by sequence, in the prototype's enumeration order
    std::vector<std::size_t> terms, ints;
    for (std::size_t i = 0; i < M; ++i) (copies[i].terminal ? terms : ints).push_back(i);
    if (terms.size() != 2) throw std::runtime_error("allocation test: a genotype needs exactly two terminal copies");
    std::vector<std::pair<std::vector<std::size_t>, std::vector<std::size_t>>> arrangements;
    std::set<std::pair<std::vector<std::string>, std::vector<std::string>>> seen;
    for (const auto& [ta, tb] : {std::make_pair(terms[0], terms[1]), std::make_pair(terms[1], terms[0])}) {
        for (std::size_t mask = 0; mask < (std::size_t{1} << ints.size()); ++mask) {
            std::vector<std::size_t> A, Bm;
            for (std::size_t i = 0; i < ints.size(); ++i) ((mask >> i) & 1 ? A : Bm).push_back(ints[i]);
            if (A.size() + 1 > lib.kmax || Bm.size() + 1 > lib.kmax) continue;
            std::vector<std::size_t> oa = A;
            do {
                std::vector<std::size_t> ob = Bm;
                do {
                    std::vector<std::size_t> ha = oa, hb = ob;
                    ha.push_back(ta); hb.push_back(tb);
                    if (copies[ha[0]].seq.compare(0, anchor.size(), anchor) != 0 || copies[hb[0]].seq.compare(0, anchor.size(), anchor) != 0) continue;
                    std::pair<std::vector<std::string>, std::vector<std::string>> key;
                    for (std::size_t i : ha) key.first.push_back(copies[i].seq);
                    for (std::size_t i : hb) key.second.push_back(copies[i].seq);
                    if (!seen.insert(key).second) continue;
                    arrangements.push_back({ha, hb});
                } while (std::next_permutation(ob.begin(), ob.end()));
            } while (std::next_permutation(oa.begin(), oa.end()));
        }
    }
    const auto block_of = [&](const std::vector<std::string>& seqs) {
        std::string s = lib.flank_left;
        for (const std::string& x : seqs) s += x;
        s += lib.flank_right;
        return s.substr(lib.trim_left, s.size() - lib.trim_left - lib.trim_right);
    };
    std::vector<std::pair<std::string, std::string>> ctx;
    {
        std::vector<std::string> a, b;
        for (std::size_t i : A0) a.push_back(copies[i].seq);
        for (std::size_t i : A1) b.push_back(copies[i].seq);
        ctx.push_back({block_of(a), block_of(b)});
    }
    const auto keyA = A0 <= A1 ? std::make_pair(A0, A1) : std::make_pair(A1, A0);
    for (const auto& ap : arrangements) {
        const auto k = ap.first <= ap.second ? std::make_pair(ap.first, ap.second) : std::make_pair(ap.second, ap.first);
        if (k == keyA) continue;
        for (int fo = 0; fo < 2; ++fo) {
            std::array<std::vector<std::string>, 2> homs;
            const std::array<const std::vector<std::size_t>*, 2> hh = {&ap.first, &ap.second};
            for (int hi = 0; hi < 2; ++hi) {
                const std::vector<std::size_t>& h = *hh[static_cast<std::size_t>(hi)];
                for (std::size_t j = 0; j < h.size(); ++j) {
                    const std::string& st = j == 0 ? flank_starts[static_cast<std::size_t>((hi + fo) % 2)] : nxt.at(h[j - 1]);
                    homs[static_cast<std::size_t>(hi)].push_back(st + body[h[j]]);
                }
            }
            ctx.push_back({block_of(homs[0]), block_of(homs[1])});
            AllocationRow row;
            std::string arr;
            for (std::size_t i = 0; i < ap.first.size(); ++i) arr += (i ? "-" : "") + std::to_string(ap.first[i]);
            arr += "|";
            for (std::size_t i = 0; i < ap.second.size(); ++i) arr += (i ? "-" : "") + std::to_string(ap.second[i]);
            row.arrangement = arr;
            row.flank_order = fo;
            row.shape = std::to_string(ap.first.size()) + "+" + std::to_string(ap.second.size());
            res.rows.push_back(row);
        }
    }
    ExactScore original;
    std::vector<std::vector<long long>> grids;
    const std::vector<ExactScore> scores = rs.lifecycle == "cold" ? score_block_candidates(in, block, ctx, rs, original, &grids)
                                                                  : score_block_candidates_context(in, block, ctx, rs, original, &grids, nullptr);
    const auto fnv = [](const std::vector<long long>& g) {
        std::uint64_t h = 1469598103934665603ULL;
        for (long long q : g) {
            unsigned char b[sizeof(long long)];
            std::memcpy(b, &q, sizeof(long long));
            for (unsigned char x : b) { h ^= x; h *= 1099511628211ULL; }
        }
        return h;
    };
    res.reference_fnv64_1e9 = fnv(grids[0]);
    res.reference_likelihood = scores[0].likelihood;
    bool any = false;
    for (std::size_t r = 0; r < res.rows.size(); ++r) {
        AllocationRow& row = res.rows[r];
        const std::vector<long long>& g = grids[r + 1];
        row.fragment_fnv64_1e9 = fnv(g);
        for (std::size_t f = 0; f < g.size(); ++f) row.fragments_differing += g[f] != grids[0][f];
        row.exposure_delta = scores[r + 1].exposure - scores[0].exposure;
        row.likelihood_delta = scores[r + 1].likelihood - scores[0].likelihood;
        row.identical = row.fragments_differing == 0 && std::fabs(row.exposure_delta) <= 1e-9;
        any |= row.identical;
    }
    res.status = any ? "UNRESOLVED" : "RESOLVED";
    res.seconds = seconds_since(started);
    return res;
}

// ---------------------------------------------------------------------------------------------- tables

namespace {
std::string join_sizes(const std::vector<std::size_t>& v, const char* sep = "-") {
    std::string s;
    for (std::size_t x : v) s += (s.empty() ? "" : sep) + std::to_string(x);
    return s;
}
}

std::string repair_candidates_table(const std::vector<std::pair<std::string, RepairResult>>& repairs) {
    std::ostringstream o;
    o << std::setprecision(17);
    o << "block_id\troute\trepair_scope\tlabel\trank\torientation\tctx1_id\tctx2_id\tctx1_bp\tctx2_bp\tunits_a\tunits_b\tmultiset\tmultiset_rank\t"
         "proposal_loglik\tstates\tpanel_log_prior\tfragment_term\texposure_term\tlikelihood\tin_top_set\tchosen_content\t"
         "allocation_tied\tchosen\tserialization_representative\n";
    for (const auto& [block_id, r] : repairs) {
        for (const BlockCandidate& c : r.candidates) {
            o << block_id << '\t' << r.route << '\t' << r.scope << '\t' << c.label << '\t' << c.rank << '\t' << c.orientation << '\t' << c.ctx1_id << '\t'
              << c.ctx2_id << '\t' << c.ctx1_sequence.size() << '\t' << c.ctx2_sequence.size() << '\t'
              << (c.units_a.empty() ? "NA" : join_sizes(c.units_a)) << '\t' << (c.units_b.empty() ? "NA" : join_sizes(c.units_b)) << '\t'
              << (c.multiset.empty() ? "NA" : join_sizes(c.multiset, "+")) << '\t' << c.multiset_rank << '\t' << c.proposal_loglik << '\t'
              << (c.states.empty() ? "NA" : c.states) << '\t' << c.panel_log_prior << '\t' << c.fragment_term << '\t' << c.exposure_term
              << '\t' << c.likelihood << '\t' << (c.in_top_set ? "true" : "false") << '\t' << (c.chosen_content ? "true" : "false") << '\t'
              << (c.allocation_tied ? "true" : "false") << '\t' << (c.chosen ? "true" : "false") << '\t'
              << (c.chosen && r.representative_only ? "REPRESENTATIVE_ONLY_NOT_RESOLVED" : c.chosen ? "chosen" : "NA") << '\n';
        }
    }
    return o.str();
}

std::string repair_blocks_table(const std::vector<std::pair<std::string, RepairResult>>& repairs) {
    std::ostringstream o;
    o << std::setprecision(17);
    o << "block_id\troute\toutcome\treason\trepair_scope\trepair_scope_reason\tcap_reached\tcandidates\tmultisets_evaluated\t"
         "multisets_retained\thaploid\tdiploid\tm_range\tinternal_capable\tterminal_capable\tbeta\tphi\tem_iterations\ttop_set_size\t"
         "orientation_tie\tallocation_tie\tallocation_representative_only\trepair_class\tdosage_content_scope\tdiplotype_scope\t"
         "allocation_equivalence_scope\tserialization_representative\tchosen_label\trepaired_ctx1_id\trepaired_ctx2_id\toriginal_fragment_term\t"
         "original_exposure_term\toriginal_likelihood\trepaired_fragment_term\trepaired_exposure_term\trepaired_likelihood\t"
         "repair_separation\tconditional_repair_gap\n";
    for (const auto& [block_id, r] : repairs) {
        std::string mr;
        for (int m : r.m_range) mr += (mr.empty() ? "" : ",") + std::to_string(m);
        const BlockCandidate* c = r.chosen >= 0 ? &r.candidates[static_cast<std::size_t>(r.chosen)] : nullptr;
        o << block_id << '\t' << r.route << '\t' << r.outcome << '\t' << r.reason << '\t' << r.scope << '\t' << r.scope_reason << '\t'
          << r.cap_reached << '\t' << r.candidates.size() << '\t' << r.multisets_evaluated << '\t' << r.multisets_retained << '\t' << r.haploid
          << '\t' << r.diploid << '\t' << (mr.empty() ? "NA" : mr) << '\t' << r.internal_capable << '\t' << r.terminal_capable << '\t' << r.beta
          << '\t' << r.phi << '\t' << r.em_iterations << '\t' << r.top_set_size << '\t' << (r.orientation_tie ? "true" : "false") << '\t'
          << r.allocation_tie << '\t' << (r.representative_only ? "true" : "false") << '\t' << r.repair_class << '\t' << r.dosage_content_scope
          << '\t' << r.diplotype_scope << '\t' << r.allocation_equivalence_scope << '\t' << (r.serialization_representative ? "true" : "false")
          << '\t' << (c ? c->label + "/" + c->orientation : "NA") << '\t' << (c ? c->ctx1_id : "NA") << '\t'
          << (c ? c->ctx2_id : "NA") << '\t' << r.original_fragment << '\t' << r.original_exposure << '\t' << r.original_likelihood << '\t';
        if (c) o << c->fragment_term << '\t' << c->exposure_term << '\t' << c->likelihood << '\t' << r.repair_separation << '\t';
        else o << "NA\tNA\tNA\tNA\t";
        if (c && r.gap_defined) o << r.conditional_repair_gap;
        else o << "NA";
        o << '\n';
    }
    return o.str();
}

} // namespace panvar::dosage
