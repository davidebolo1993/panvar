#include "panvar/genotype_dosage.hpp"

#include "panvar/parallel.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace panvar::dosage {

// ------------------------------------------------------------------------------------------------ Python random

PythonRandom::PythonRandom(std::uint32_t seed) {
    // init_genrand(19650218) then init_by_array({seed}), exactly as CPython's random_seed for a small integer.
    std::array<std::uint32_t, 624>& mt = state_;
    mt[0] = 19650218u;
    for (std::uint32_t i = 1; i < 624; ++i) mt[i] = 1812433253u * (mt[i - 1] ^ (mt[i - 1] >> 30)) + i;
    const std::uint32_t key[1] = {seed};
    std::size_t i = 1, j = 0;
    for (std::size_t k = 624; k > 0; --k) {
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1664525u)) + key[j] + static_cast<std::uint32_t>(j);
        ++i; ++j;
        if (i >= 624) { mt[0] = mt[623]; i = 1; }
        if (j >= 1) j = 0;
    }
    for (std::size_t k = 623; k > 0; --k) {
        mt[i] = (mt[i] ^ ((mt[i - 1] ^ (mt[i - 1] >> 30)) * 1566083941u)) - static_cast<std::uint32_t>(i);
        ++i;
        if (i >= 624) { mt[0] = mt[623]; i = 1; }
    }
    mt[0] = 0x80000000u;
    index_ = 624;
}

std::uint32_t PythonRandom::next_uint32() {
    std::array<std::uint32_t, 624>& mt = state_;
    if (index_ >= 624) {
        for (std::size_t kk = 0; kk < 624; ++kk) {
            const std::uint32_t y = (mt[kk] & 0x80000000u) | (mt[(kk + 1) % 624] & 0x7fffffffu);
            mt[kk] = mt[(kk + 397) % 624] ^ (y >> 1) ^ ((y & 1u) ? 0x9908b0dfu : 0u);
        }
        index_ = 0;
    }
    std::uint32_t y = mt[index_++];
    y ^= y >> 11;
    y ^= (y << 7) & 0x9d2c5680u;
    y ^= (y << 15) & 0xefc60000u;
    y ^= y >> 18;
    return y;
}

std::uint32_t PythonRandom::randrange(std::uint32_t n) {
    if (n == 0) throw std::invalid_argument("randrange of an empty range");
    int bits = 0;
    for (std::uint32_t v = n; v != 0; v >>= 1) ++bits;
    for (;;) {
        const std::uint32_t r = next_uint32() >> (32 - bits);
        if (r < n) return r;
    }
}

// ------------------------------------------------------------------------------------------------ normal quantile

double normal_inverse_cdf(double p) {
    // CPython statistics._normal_dist_inv_cdf (Wichura AS241, PPND16), mu 0, sigma 1. CPython evaluates it in C built
    // with clang's default contraction on arm64, so this one function contracts (the rest of the file does not).
#pragma clang fp contract(on)
    if (!(p > 0.0 && p < 1.0)) throw std::invalid_argument("normal_inverse_cdf needs 0 < p < 1");
    const double q = p - 0.5;
    if (std::fabs(q) <= 0.425) {
        const double r = 0.180625 - q * q;
        const double num = (((((((2.5090809287301226727e+3 * r + 3.3430575583588128105e+4) * r + 6.7265770927008700853e+4) * r +
                                 4.5921953931549871457e+4) * r + 1.3731693765509461125e+4) * r + 1.9715909503065514427e+3) * r +
                              1.3314166789178437745e+2) * r + 3.3871328727963666080e+0) * q;
        const double den = (((((((5.2264952788528545610e+3 * r + 2.8729085735721942674e+4) * r + 3.9307895800092710610e+4) * r +
                                 2.1213794301586595867e+4) * r + 5.3941960214247511077e+3) * r + 6.8718700749205790830e+2) * r +
                              4.2313330701600911252e+1) * r + 1.0);
        return num / den;
    }
    double r = q <= 0.0 ? p : 1.0 - p;
    r = std::sqrt(-std::log(r));
    double x;
    if (r <= 5.0) {
        r = r - 1.6;
        const double num = (((((((7.74545014278341407640e-4 * r + 2.27238449892691845833e-2) * r + 2.41780725177450611770e-1) * r +
                                 1.27045825245236838258e+0) * r + 3.64784832476320460504e+0) * r + 5.76949722146069140550e+0) * r +
                              4.63033784615654529590e+0) * r + 1.42343711074968357734e+0);
        const double den = (((((((1.05075007164441684324e-9 * r + 5.47593808499534494600e-4) * r + 1.51986665636164571966e-2) * r +
                                 1.48103976427480074590e-1) * r + 6.89767334985100004550e-1) * r + 1.67638483018380384940e+0) * r +
                              2.05319162663775882187e+0) * r + 1.0);
        x = num / den;
    } else {
        r = r - 5.0;
        const double num = (((((((2.01033439929228813265e-7 * r + 2.71155556874348757815e-5) * r + 1.24266094738807843860e-3) * r +
                                 2.65321895265761230930e-2) * r + 2.96560571828504891230e-1) * r + 1.78482653991729133580e+0) * r +
                              5.46378491116411436990e+0) * r + 6.65790464350110377720e+0);
        const double den = (((((((2.04426310338993978564e-15 * r + 1.42151175831644588870e-7) * r + 1.84631831751005468180e-5) * r +
                                 7.86869131145613259100e-4) * r + 1.48753612908506148525e-2) * r + 1.36929880922735805310e-1) * r +
                              5.99832206555887937690e-1) * r + 1.0);
        x = num / den;
    }
    if (q < 0.0) x = -x;
    return x;
}

std::optional<std::size_t> molecule_span_bp(double insert_mean_bp, double insert_sd_bp, const DosageSettings& settings) {
    if (!(insert_mean_bp > 0.0) || !(insert_sd_bp >= 0.0) || !std::isfinite(insert_mean_bp) || !std::isfinite(insert_sd_bp)) {
        return std::nullopt;
    }
    return static_cast<std::size_t>(std::ceil(insert_mean_bp + settings.molecule_span_sigmas * insert_sd_bp));
}

// ------------------------------------------------------------------------------------------------ k-mers

namespace {

int base_code(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return -1;
    }
}

} // namespace

std::vector<KmerAt> canonical_kmers(const std::string& sequence, std::size_t k) {
    if (k == 0 || k > 32) throw std::invalid_argument("k-mer size must be in 1..32");
    std::vector<KmerAt> out;
    if (sequence.size() < k) return out;
    out.reserve(sequence.size() - k + 1);
    const Kmer mask = k == 32 ? ~Kmer{0} : ((Kmer{1} << (2 * k)) - 1);
    const unsigned shift = static_cast<unsigned>(2 * (k - 1));
    Kmer fwd = 0, rev = 0;
    std::size_t valid = 0;
    for (std::size_t i = 0; i < sequence.size(); ++i) {
        const int c = base_code(sequence[i]);
        if (c < 0) { valid = 0; fwd = rev = 0; continue; }
        fwd = ((fwd << 2) | static_cast<Kmer>(c)) & mask;
        rev = (rev >> 2) | (static_cast<Kmer>(3 - c) << shift);
        if (++valid >= k) out.push_back({std::min(fwd, rev), i + 1 - k});
    }
    return out;
}

// ------------------------------------------------------------------------------------------------ census

std::vector<std::vector<NodeStep>> all_walks(const CanonicalBlockAllele& allele) {
    std::vector<std::vector<NodeStep>> walks;
    for (const std::string& source : allele.source_ids) {
        const std::size_t at = source.find("walk:");
        if (at == std::string::npos) continue;
        std::vector<NodeStep> walk;
        std::stringstream steps(source.substr(at + 5));
        for (std::string step; std::getline(steps, step, ',');) {
            if (step.size() < 2 || (step.back() != '+' && step.back() != '-')) {
                throw std::runtime_error("malformed walk step '" + step + "' in allele source " + source);
            }
            walk.push_back({step.substr(0, step.size() - 1), step.back() == '+'});
        }
        walks.push_back(std::move(walk));
    }
    return walks;
}

std::vector<NodeStep> primary_walk(const CanonicalBlockAllele& allele) {
    std::vector<std::vector<NodeStep>> walks = all_walks(allele);
    if (walks.empty()) return {};
    return std::move(walks.front());
}

const char* census_route_name(CensusRoute route) {
    switch (route) {
        case CensusRoute::NotSensitive: return "not_sensitive";
        case CensusRoute::Pair: return "pair";
        case CensusRoute::RepeatedUnit: return "repeated_unit";
    }
    return "unknown";
}

std::vector<CensusBlock> run_census(const PreparedGenotypeModel& model, const NodeSequences& nodes,
                                    std::size_t molecule_span_bp, const DosageSettings& settings) {
    std::vector<CensusBlock> out;
    for (std::size_t b = 0; b < model.blocks().size(); ++b) {
        const PreparedGenotypeBlock& block = model.blocks()[b];
        CensusBlock row;
        row.block = b;
        row.block_id = block.block_id;
        std::vector<std::map<std::string, std::size_t>> counts;   // allele -> node -> visits
        for (const CanonicalBlockAllele& allele : block.alleles) {
            const std::vector<NodeStep> walk = primary_walk(allele);
            if (walk.empty()) continue;
            std::map<std::string, std::size_t> c;
            for (const NodeStep& s : walk) ++c[s.node];
            counts.push_back(std::move(c));
        }
        row.alleles = counts.size();
        const auto weight = [&](const std::string& node) -> std::size_t {
            const auto found = nodes.find(node);
            if (found == nodes.end()) throw std::runtime_error("census: walk node " + node + " has no sequence");
            return found->second.size();
        };
        std::vector<std::size_t> length(counts.size(), 0);
        for (std::size_t a = 0; a < counts.size(); ++a) {
            std::size_t repeated = 0;
            for (const auto& [node, visits] : counts[a]) {
                length[a] += visits * weight(node);
                row.max_node_multiplicity = std::max(row.max_node_multiplicity, visits);
                if (visits >= 2) repeated += weight(node);
            }
            row.max_repeated_bp = std::max(row.max_repeated_bp, repeated);
        }
        for (std::size_t i = 0; i < counts.size(); ++i) {
            for (std::size_t j = 0; j < counts.size(); ++j) {
                if (length[j] < length[i] + molecule_span_bp) continue;
                std::size_t novel = 0;   // bp of distinct nodes in i absent from j
                for (const auto& [node, visits] : counts[i]) {
                    if (!counts[j].count(node)) novel += weight(node);
                }
                if (novel <= settings.max_novel_bp) ++row.dosage_only_pairs;
            }
        }
        if (row.dosage_only_pairs > 0) {
            row.route = row.max_node_multiplicity >= 2 && row.max_repeated_bp > molecule_span_bp ? CensusRoute::RepeatedUnit
                                                                                                  : CensusRoute::Pair;
        }
        out.push_back(std::move(row));
    }
    return out;
}

// ------------------------------------------------------------------------------------------------ features

BlockFeatureModel pooled_features(const PreparedGenotypeBlock& block, const DosageSettings& settings) {
    const std::size_t na = block.alleles.size();
    for (std::size_t a = 0; a < na; ++a) {
        if (block.alleles[a].id != a) throw std::runtime_error("pooled_features: allele IDs are not dense in block " + block.block_id);
    }
    std::unordered_map<Kmer, std::uint32_t> kid;
    std::vector<Kmer> order;
    std::vector<std::vector<std::uint32_t>> per(na);
    for (std::size_t a = 0; a < na; ++a) {
        for (const KmerAt& x : canonical_kmers(block.alleles[a].sequence, settings.kmer)) {
            const auto [it, added] = kid.emplace(x.kmer, static_cast<std::uint32_t>(order.size()));
            if (added) order.push_back(x.kmer);
            per[a].push_back(it->second);
        }
    }
    const std::size_t nk = order.size();
    std::vector<std::uint16_t> counts(nk * na, 0);   // k-mer-major
    for (std::size_t a = 0; a < na; ++a) {
        for (std::uint32_t k : per[a]) {
            std::uint16_t& c = counts[static_cast<std::size_t>(k) * na + a];
            if (c == std::numeric_limits<std::uint16_t>::max()) throw std::runtime_error("k-mer count overflow");
            ++c;
        }
    }
    std::map<std::vector<std::uint16_t>, std::size_t> group_of;
    std::vector<std::vector<std::uint16_t>> columns;
    std::vector<std::vector<Kmer>> members;
    for (std::size_t k = 0; k < nk; ++k) {
        const std::uint16_t* col = &counts[k * na];
        if (std::all_of(col, col + na, [&](std::uint16_t v) { return v == col[0]; })) continue;
        std::vector<std::uint16_t> column(col, col + na);
        const auto [it, added] = group_of.emplace(column, columns.size());
        if (added) { columns.push_back(std::move(column)); members.emplace_back(); }
        members[it->second].push_back(order[k]);
    }
    BlockFeatureModel out;
    for (std::size_t g = 0; g < columns.size(); ++g) {
        if (members[g].size() < settings.min_feature_kmers) continue;
        PooledFeature f;
        f.copies = columns[g];
        f.kmers = std::move(members[g]);
        std::sort(f.kmers.begin(), f.kmers.end());
        out.features.push_back(std::move(f));
    }
    out.block_kmers = order;
    std::sort(out.block_kmers.begin(), out.block_kmers.end());
    return out;
}

ReportedCall reported_call(const PreparedGenotypeModel& model, const MosaicPath& first, const MosaicPath& second) {
    ReportedCall call;
    const HaplotypeSequence s1 = spell_mosaic(model, first), s2 = spell_mosaic(model, second);
    call.haplotypes = {s1.bases, s2.bases};
    call.intervals = {s1.block_intervals, s2.block_intervals};
    call.alleles = {first.alleles, second.alleles};
    return call;
}

const char* dosage_status_name(DosageStatus status) {
    switch (status) {
        case DosageStatus::Consistent: return "CONSISTENT";
        case DosageStatus::Repaired: return "REPAIRED";
        case DosageStatus::TriggeredUnrepaired: return "TRIGGERED_UNREPAIRED";
        case DosageStatus::AutoRepair: return "AUTO_REPAIR";
        case DosageStatus::Suspect: return "SUSPECT";
        case DosageStatus::Underpowered: return "UNDERPOWERED";
        case DosageStatus::NotSensitive: return "NOT_SENSITIVE";
        case DosageStatus::NotEvaluated: return "NOT_EVALUATED";
    }
    return "UNKNOWN";
}

// ------------------------------------------------------------------------------------------------ audit

namespace {

bool contains(const std::vector<Kmer>& sorted, Kmer x) { return std::binary_search(sorted.begin(), sorted.end(), x); }

struct PairHits {
    // dense hit counts per retained pair: index 0 the reference, j + 1 feature j
    std::size_t width = 0;
    std::vector<std::uint32_t> values;
    std::size_t pairs() const noexcept { return width == 0 ? 0 : values.size() / width; }
};

double estimate(std::uint64_t t_feature, std::size_t feature_size, double r0) {
    // 2 * (t / size) / r0, in the prototype's evaluation order
    return 2.0 * (static_cast<double>(t_feature) / static_cast<double>(feature_size)) / r0;
}

} // namespace

DosageAudit audit_dosage(const PreparedGenotypeModel& model, const std::vector<CensusBlock>& census,
                         const std::vector<BlockFeatureModel>& features, const ReportedCall& call,
                         const std::vector<ReadPair>& reads, std::size_t molecule_span_bp,
                         const DosageSettings& settings, std::size_t threads) {
    const std::size_t nb = model.blocks().size();
    if (census.size() != nb || features.size() != nb) throw std::runtime_error("audit_dosage: census/features do not match the model");
    DosageAudit out;
    out.molecule_span_bp = molecule_span_bp;
    out.alpha = settings.alpha;
    const std::size_t k = settings.kmer;

    // per homologue: every canonical k-mer with its position, and its full-haplotype counts
    std::array<std::vector<KmerAt>, 2> hk;
    std::array<std::unordered_map<Kmer, std::uint32_t>, 2> hcount;
    for (int h = 0; h < 2; ++h) {
        hk[h] = canonical_kmers(call.haplotypes[h], k);
        for (const KmerAt& x : hk[h]) ++hcount[h][x.kmer];
    }
    // reference candidates: k-mers inside census-invariant block intervals of each homologue
    std::array<std::unordered_set<Kmer>, 2> invariant;
    for (int h = 0; h < 2; ++h) {
        for (const CensusBlock& c : census) {
            if (c.sensitive()) continue;
            const SequenceInterval iv = call.intervals[h][c.block];
            for (const KmerAt& x : hk[h]) {
                if (x.position >= iv.begin && x.position + k <= iv.end) invariant[h].insert(x.kmer);
            }
        }
    }

    for (const CensusBlock& c : census) {
        BlockAudit ba;
        ba.block = c.block;
        ba.block_id = c.block_id;
        ba.route = c.route;
        ba.call_alleles = {call.alleles[0][c.block], call.alleles[1][c.block]};
        if (!c.sensitive()) { ba.status = DosageStatus::NotSensitive; out.blocks.push_back(std::move(ba)); continue; }
        const BlockFeatureModel& fm = features[c.block];

        // per-call exclusion: k-mers of either homologue entirely outside this block's interval
        std::unordered_set<Kmer> outside;
        for (int h = 0; h < 2; ++h) {
            const SequenceInterval iv = call.intervals[h][c.block];
            for (const KmerAt& x : hk[h]) {
                if (x.position + k <= iv.begin || x.position >= iv.end) outside.insert(x.kmer);
            }
        }
        struct Tested { std::size_t index; std::vector<Kmer> kmers; };
        std::vector<Tested> tested;
        for (std::size_t j = 0; j < fm.features.size(); ++j) {
            std::vector<Kmer> kept;
            for (Kmer x : fm.features[j].kmers) if (!outside.count(x)) kept.push_back(x);
            if (kept.size() >= settings.min_feature_kmers) tested.push_back({j, std::move(kept)});
        }
        if (tested.empty()) {
            ba.status = DosageStatus::NotEvaluated; ba.reason = "no_tested_feature"; out.blocks.push_back(std::move(ba)); continue;
        }
        // single-copy reference
        std::vector<Kmer> reference;
        for (Kmer x : invariant[0]) {
            if (!invariant[1].count(x)) continue;
            const auto c0 = hcount[0].find(x), c1 = hcount[1].find(x);
            if (c0 == hcount[0].end() || c1 == hcount[1].end() || c0->second != 1 || c1->second != 1) continue;
            if (contains(fm.block_kmers, x)) continue;
            reference.push_back(x);
        }
        std::sort(reference.begin(), reference.end());
        ba.reference_kmers = reference.size();
        if (reference.empty()) {
            ba.status = DosageStatus::NotEvaluated; ba.reason = "no_reference"; out.blocks.push_back(std::move(ba)); continue;
        }
        std::unordered_map<Kmer, std::uint32_t> lookup;
        lookup.reserve(reference.size() * 2);
        for (Kmer x : reference) lookup[x] = 0;
        for (std::size_t t = 0; t < tested.size(); ++t) for (Kmer x : tested[t].kmers) lookup[x] = static_cast<std::uint32_t>(t + 1);

        // hits per read pair, counted in contiguous ranges and concatenated in input order
        const std::size_t width = tested.size() + 1;
        const std::size_t chunks = std::max<std::size_t>(1, std::min<std::size_t>(reads.size(), 64));
        std::vector<std::vector<std::uint32_t>> part(chunks);
        run_parallel(chunks, threads, [&](std::size_t ci) {
            const std::size_t lo = reads.size() * ci / chunks, hi = reads.size() * (ci + 1) / chunks;
            std::vector<std::uint32_t> v(width);
            for (std::size_t r = lo; r < hi; ++r) {
                std::fill(v.begin(), v.end(), 0);
                bool any = false;
                for (const std::string* s : {&reads[r].first.bases, &reads[r].second.bases}) {
                    for (const KmerAt& x : canonical_kmers(*s, k)) {
                        const auto found = lookup.find(x.kmer);
                        if (found != lookup.end()) { ++v[found->second]; any = true; }
                    }
                }
                if (any) part[ci].insert(part[ci].end(), v.begin(), v.end());
            }
        });
        PairHits hits;
        hits.width = width;
        for (auto& p : part) hits.values.insert(hits.values.end(), p.begin(), p.end());
        const std::size_t n = hits.pairs();
        ba.retained_pairs = n;
        std::vector<std::uint64_t> total(width, 0);
        for (std::size_t i = 0; i < n; ++i) for (std::size_t j = 0; j < width; ++j) total[j] += hits.values[i * width + j];
        ba.reference_hits = total[0];
        if (total[0] == 0) {
            ba.status = DosageStatus::NotEvaluated; ba.reason = "no_reference_hits"; out.blocks.push_back(std::move(ba)); continue;
        }
        const double r0 = static_cast<double>(total[0]) / static_cast<double>(reference.size());
        std::vector<double> point(tested.size());
        for (std::size_t t = 0; t < tested.size(); ++t) point[t] = estimate(total[t + 1], tested[t].kmers.size(), r0);

        // bootstrap: indices drawn sequentially from one stream, re-seeded per block
        if (n > std::numeric_limits<std::uint32_t>::max()) throw std::runtime_error("bootstrap population too large");
        PythonRandom rng(settings.bootstrap_seed);
        const std::size_t R = settings.bootstrap_replicates;
        std::vector<std::vector<double>> boot(tested.size(), std::vector<double>(R));
        bool zero_reference = false;
        std::vector<std::uint64_t> t(width);
        for (std::size_t rep = 0; rep < R; ++rep) {
            std::fill(t.begin(), t.end(), 0);
            for (std::size_t i = 0; i < n; ++i) {
                const std::size_t pick = rng.randrange(static_cast<std::uint32_t>(n));
                for (std::size_t j = 0; j < width; ++j) t[j] += hits.values[pick * width + j];
            }
            if (t[0] == 0) { zero_reference = true; continue; }
            const double rr = static_cast<double>(t[0]) / static_cast<double>(reference.size());
            for (std::size_t f = 0; f < tested.size(); ++f) boot[f][rep] = estimate(t[f + 1], tested[f].kmers.size(), rr);
        }
        if (zero_reference) {
            ba.status = DosageStatus::NotEvaluated; ba.reason = "bootstrap_zero_reference"; out.blocks.push_back(std::move(ba)); continue;
        }
        static constexpr std::array<double, 5> P = {0.05, 0.25, 0.50, 0.75, 0.95};
        for (std::size_t f = 0; f < tested.size(); ++f) {
            const PooledFeature& pf = fm.features[tested[f].index];
            FeatureAudit fa;
            fa.feature = f;
            fa.kmers = tested[f].kmers.size();
            fa.copy_values = pf.copies;
            std::sort(fa.copy_values.begin(), fa.copy_values.end());
            fa.copy_values.erase(std::unique(fa.copy_values.begin(), fa.copy_values.end()), fa.copy_values.end());
            fa.feasible_min = 2 * static_cast<int>(fa.copy_values.front());
            fa.feasible_max = 2 * static_cast<int>(fa.copy_values.back());
            fa.call_cn = static_cast<int>(pf.copies.at(ba.call_alleles[0])) + static_cast<int>(pf.copies.at(ba.call_alleles[1]));
            fa.cn_estimate = point[f];
            std::vector<double> sorted = boot[f];
            double mean = 0.0;
            for (double x : sorted) mean += x;
            mean /= static_cast<double>(R);
            double ss = 0.0;
            for (double x : sorted) ss += (x - mean) * (x - mean);
            fa.bootstrap_sd = R > 1 ? std::sqrt(ss / static_cast<double>(R - 1)) : 0.0;
            std::sort(sorted.begin(), sorted.end());
            for (std::size_t q = 0; q < P.size(); ++q) fa.quantiles[q] = sorted[static_cast<std::size_t>(P[q] * static_cast<double>(R - 1))];
            const double lo = fa.quantiles[0], hi = fa.quantiles[4];
            for (long x = static_cast<long>(std::ceil(lo)); x <= static_cast<long>(std::floor(hi)); ++x) {
                if (x >= fa.feasible_min && x <= fa.feasible_max) fa.admissible.push_back(static_cast<int>(x));
            }
            if (fa.admissible.empty()) {
                int best = fa.feasible_min;
                double best_d = std::numeric_limits<double>::infinity();
                for (int x = fa.feasible_min; x <= fa.feasible_max; ++x) {
                    const double d = std::min(std::fabs(x - lo), std::fabs(x - hi));
                    if (d < best_d) { best_d = d; best = x; }
                }
                fa.admissible.push_back(best);
            }
            fa.call_admissible = std::find(fa.admissible.begin(), fa.admissible.end(), fa.call_cn) != fa.admissible.end();
            fa.robust_sd = (fa.quantiles[3] - fa.quantiles[1]) / settings.iqr_to_sd;
            const double dev = std::fabs(static_cast<double>(fa.call_cn) - fa.quantiles[2]);
            fa.robust_z = fa.robust_sd > 0.0 ? dev / fa.robust_sd : (dev == 0.0 ? 0.0 : std::numeric_limits<double>::infinity());
            ba.features.push_back(std::move(fa));
        }
        out.m += ba.features.size();
        out.blocks.push_back(std::move(ba));
    }

    // classification, once m is known
    if (out.m > 0) {
        const double zc = normal_inverse_cdf(1.0 - settings.alpha / (2.0 * static_cast<double>(out.m)));
        out.z_critical = zc;
        const double sd_max = 1.0 / (zc + normal_inverse_cdf(settings.power));
        for (BlockAudit& ba : out.blocks) {
            if (ba.features.empty()) continue;
            bool trig = false, susp = false, under = false;
            for (FeatureAudit& fa : ba.features) {
                fa.trigger = fa.robust_z >= zc && !fa.call_admissible;
                fa.suspect = fa.robust_z >= settings.suspect_z;
                fa.underpowered = fa.robust_sd > sd_max;
                trig |= fa.trigger; susp |= fa.suspect; under |= fa.underpowered;
                ba.max_robust_z = std::max(ba.max_robust_z, fa.robust_z);
            }
            ba.status = trig ? DosageStatus::AutoRepair : susp ? DosageStatus::Suspect : under ? DosageStatus::Underpowered
                                                                                               : DosageStatus::Consistent;
        }
    }
    return out;
}

// ------------------------------------------------------------------------------------------------ tables

namespace {

std::string join_ints(const std::vector<int>& v) {
    std::string s;
    for (int x : v) s += (s.empty() ? "" : ",") + std::to_string(x);
    return s;
}

std::string optional_double(const std::optional<double>& v) {
    if (!v) return "NA";
    std::ostringstream o;
    o << std::setprecision(17) << *v;
    return o.str();
}

} // namespace

std::string dosage_features_table(const DosageAudit& audit) {
    std::ostringstream o;
    o << std::setprecision(17);
    o << "block_id\tfeature\tkmers\tcount_values\tcall_cn\tcn_estimate\tbootstrap_sd\tq05\tq25\tq50\tq75\tq95\tadmissible\t"
         "call_admissible\trobust_sd\trobust_z\tm\talpha\tz_critical\ttrigger\tsuspect\tunderpowered\n";
    for (const BlockAudit& ba : audit.blocks) {
        for (const FeatureAudit& fa : ba.features) {
            std::string values;
            for (std::uint16_t v : fa.copy_values) values += (values.empty() ? "" : ",") + std::to_string(v);
            o << ba.block_id << '\t' << fa.feature << '\t' << fa.kmers << '\t' << values << '\t' << fa.call_cn << '\t' << fa.cn_estimate
              << '\t' << fa.bootstrap_sd;
            for (double q : fa.quantiles) o << '\t' << q;
            o << '\t' << join_ints(fa.admissible) << '\t' << (fa.call_admissible ? "true" : "false") << '\t' << fa.robust_sd << '\t'
              << fa.robust_z << '\t' << audit.m << '\t' << audit.alpha << '\t' << optional_double(audit.z_critical) << '\t'
              << (fa.trigger ? "true" : "false") << '\t' << (fa.suspect ? "true" : "false") << '\t'
              << (fa.underpowered ? "true" : "false") << '\n';
        }
    }
    return o.str();
}

std::string dosage_blocks_table(const DosageAudit& audit) {
    std::ostringstream o;
    o << std::setprecision(17);
    o << "block_id\tcensus_route\tmolecule_span_bp\tdosage_status\tdosage_reason\toriginal_allele1\toriginal_allele2\t"
         "features_tested\treference_kmers\treference_hits\tretained_pairs\tmax_robust_z\tm\talpha\tz_critical\n";
    for (const BlockAudit& ba : audit.blocks) {
        o << ba.block_id << '\t' << census_route_name(ba.route) << '\t' << audit.molecule_span_bp << '\t' << dosage_status_name(ba.status)
          << '\t' << ba.reason << '\t' << ba.call_alleles[0] << '\t' << ba.call_alleles[1] << '\t' << ba.features.size() << '\t'
          << ba.reference_kmers << '\t' << ba.reference_hits << '\t' << ba.retained_pairs << '\t' << ba.max_robust_z << '\t' << audit.m
          << '\t' << audit.alpha << '\t' << optional_double(audit.z_critical) << '\n';
    }
    return o.str();
}

} // namespace panvar::dosage
