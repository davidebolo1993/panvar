#include "panvar/genotype_reads.hpp"

#include "panvar/syncmer.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <fstream>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <unordered_map>

#include <zlib.h>
#include <kseq.h>

KSEQ_INIT(gzFile, gzread)

namespace panvar {
namespace {

std::uint64_t mix64(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

// Must match genotype_markers.cpp exactly, or reads and panel disagree about what an adjacency is.
std::uint64_t adjacency_key(std::uint64_t a, std::uint64_t b, std::size_t gap) {
    const std::uint64_t lo = std::min(a, b);
    const std::uint64_t hi = std::max(a, b);
    const std::uint64_t g = static_cast<std::uint64_t>(std::min<std::size_t>(gap, 0xffff));
    return mix64(mix64(lo) ^ (mix64(hi) * 0x100000001b3ULL) ^ (g * 0x9e3779b97f4a7c15ULL));
}

double median_of(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    return v[v.size() / 2];
}

double mean_of(const std::vector<double>& v) {
    if (v.empty()) return 0.0;
    return std::accumulate(v.begin(), v.end(), 0.0) / static_cast<double>(v.size());
}

// Mean of the central (1 - 2*trim) fraction. Anchor counts have a right tail from markers that are
// invariant across the panel but not unique in the genome, which is what a plain mean would follow and
// a median would ignore entirely; trimming keeps the resolution of a mean without the tail.
double trimmed_mean_of(std::vector<double> v, double trim) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const std::size_t drop = static_cast<std::size_t>(trim * static_cast<double>(v.size()));
    if (2 * drop >= v.size()) return median_of(v);
    return std::accumulate(v.begin() + static_cast<std::ptrdiff_t>(drop),
                           v.end() - static_cast<std::ptrdiff_t>(drop), 0.0) /
           static_cast<double>(v.size() - 2 * drop);
}

constexpr double kTrimFraction = 0.1;

double central_value(const std::vector<double>& v, DepthEstimator estimator) {
    switch (estimator) {
        case DepthEstimator::Mean:        return mean_of(v);
        case DepthEstimator::TrimmedMean: return trimmed_mean_of(v, kTrimFraction);
        case DepthEstimator::Median:      break;
    }
    return median_of(v);
}

} // namespace

ReadCounts count_reads(
    const std::vector<std::string>& read_paths,
    const ReadPanel& panel,
    std::size_t threads) {

    std::unordered_map<std::uint64_t, std::uint32_t> node_slot;
    std::unordered_map<std::uint64_t, std::uint32_t> edge_slot;
    node_slot.reserve(panel.node_codes.size() * 2);
    edge_slot.reserve(panel.edge_keys.size() * 2);
    for (std::uint32_t i = 0; i < panel.node_codes.size(); ++i) node_slot.emplace(panel.node_codes[i], i);
    for (std::uint32_t i = 0; i < panel.edge_keys.size(); ++i) edge_slot.emplace(panel.edge_keys[i], i);
    // Novelty is judged against every adjacency the panel carries, not the retained informative subset.
    const std::unordered_set<std::uint64_t> panel_adjacencies(panel.all_edge_keys.begin(),
                                                             panel.all_edge_keys.end());

    const std::size_t nthreads =
        std::max<std::size_t>(1, threads != 0 ? threads : std::thread::hardware_concurrency());

    ReadCounts total;
    total.node.assign(panel.node_codes.size(), 0);
    total.edge.assign(panel.edge_keys.size(), 0);

    const std::size_t k = panel.kmer_size;
    const std::size_t s = panel.syncmer_s != 0 ? panel.syncmer_s : default_syncmer_s(k);

    std::mutex merge_mu;
    for (const std::string& path : read_paths) {
        gzFile fp = gzopen(path.c_str(), "r");
        if (fp == nullptr) throw std::runtime_error("genotype: cannot open reads " + path);
        kseq_t* seq = kseq_init(fp);

        // One reader (kseq is not thread-safe) feeding batches to workers, each with its own
        // accumulator; `run_parallel` is not used here because it requires distinct output indices.
        std::mutex read_mu;
        std::vector<std::thread> pool;
        for (std::size_t t = 0; t < nthreads; ++t) {
            pool.emplace_back([&] {
                ReadCounts local;
                local.node.assign(panel.node_codes.size(), 0);
                local.edge.assign(panel.edge_keys.size(), 0);
                std::vector<std::string> batch;
                for (;;) {
                    batch.clear();
                    {
                        std::lock_guard<std::mutex> lock(read_mu);
                        while (batch.size() < 4096 && kseq_read(seq) >= 0) {
                            batch.emplace_back(seq->seq.s, seq->seq.l);
                        }
                    }
                    if (batch.empty()) break;
                    for (const std::string& r : batch) {
                        ++local.reads;
                        local.bases += r.size();
                        const std::vector<KmerOccurrence> sy = panel.all_kmers
                            ? collect_canonical_kmer_occurrences(r, k) : collect_syncmers(r, k, s);
                        local.syncmers += sy.size();
                        for (const KmerOccurrence& o : sy) {
                            const auto it = node_slot.find(o.code);
                            if (it != node_slot.end()) { ++local.node[it->second]; ++local.matched_syncmers; }
                        }
                        for (std::size_t i = 1; i < sy.size(); ++i) {
                            const std::uint64_t key =
                                adjacency_key(sy[i - 1].code, sy[i].code, sy[i].start - sy[i - 1].start);
                            const auto it = edge_slot.find(key);
                            if (it != edge_slot.end()) ++local.edge[it->second];
                            if (!panel_adjacencies.empty() && !panel_adjacencies.count(key))
                                ++local.novel_adjacencies;
                        }
                    }
                }
                std::lock_guard<std::mutex> lock(merge_mu);
                for (std::size_t i = 0; i < total.node.size(); ++i) total.node[i] += local.node[i];
                for (std::size_t i = 0; i < total.edge.size(); ++i) total.edge[i] += local.edge[i];
                total.reads += local.reads;
                total.bases += local.bases;
                total.syncmers += local.syncmers;
                total.matched_syncmers += local.matched_syncmers;
                total.novel_adjacencies += local.novel_adjacencies;
            });
        }
        for (std::thread& th : pool) th.join();
        kseq_destroy(seq);
        gzclose(fp);
    }
    return total;
}

std::vector<BlockDepth> estimate_depth(
    const ReadPanel& panel,
    const ReadCounts& counts,
    std::size_t min_anchors,
    double uneven_tolerance,
    DepthModel model,
    double depth_quantile,
    std::size_t region_bp,
    DepthEstimator estimator,
    DepthRegionStats* region_stats) {

    std::vector<BlockDepth> out(panel.anchor_slots.size());
    std::vector<double> all_anchor_counts;

    for (std::size_t bi = 0; bi < panel.anchor_slots.size(); ++bi) {
        BlockDepth& d = out[bi];
        d.block_index = bi;
        std::vector<double> v;
        v.reserve(panel.anchor_slots[bi].size());
        for (const std::uint32_t slot : panel.anchor_slots[bi]) {
            v.push_back(static_cast<double>(counts.node[slot]));
            all_anchor_counts.push_back(static_cast<double>(counts.node[slot]));
        }
        d.n_anchor = v.size();
        // Every block with ANY anchors gets a local estimate. The old `v.size() < min_anchors`
        // skip contradicted the shrinkage below: that pass exists precisely so a block with a handful
        // of anchors is pulled toward the region rather than falling off a cliff, but a block that
        // never computed a local value has nothing to pull. At lpa that discarded 61 anchors across
        // six bubble blocks. `min_anchors` now flags weak local evidence instead of gating it.
        if (v.empty()) continue;
        d.low_anchor = v.size() < min_anchors;
        // `anchor_median` stays the median whatever the estimator, so the raw audit column means the
        // same statistic across runs and remains comparable when the estimator changes. `median` is
        // the model's local centre and follows the estimator.
        const double med = median_of(v);
        d.anchor_median = med;
        d.local_center = central_value(v, estimator);
        d.median = d.local_center;
        d.local_available = true;
        std::vector<double> dev;
        dev.reserve(v.size());
        for (const double x : v) dev.push_back(std::fabs(x - med));
        d.mad = median_of(dev);
        d.lambda_hap = d.median / 2.0;      // diploid: an invariant marker is carried by both copies
        d.usable = d.lambda_hap > 0.0;
        d.uneven = d.median > 0.0 && (d.mad / d.median) > uneven_tolerance;
        // Provisional: the shrinkage pass below rewrites this for every block it touches. It survives
        // only when there is no region estimate to shrink toward at all.
        if (d.usable) { d.source = DepthSource::Local; d.region_shrink_weight = 0.0; }
    }

    // Region-wide fallback for blocks with too few anchors of their own, and shrinkage toward it for
    // the rest.
    //
    // A hard threshold is the wrong shape here: a block just over it gets its own noisy estimate at
    // full weight, a block just under it gets the region estimate at full weight, and the two answers
    // can differ by a third. That matters most exactly where it hurts most -- a tandem array has few
    // invariant markers (its unit syncmers vary with copy number, so they are not anchors), so its
    // depth rests on a handful of flanking anchors, and depth is the denominator that CONVERTS marker
    // multiplicity into copy number. Measured on a synthetic VNTR: 21 anchors gave lambda 15.5 against
    // a region-wide 11.5, and the 35% inflation came straight back out as a copy number 35% too low --
    // 12 true copies called as 9.
    //
    // Shrink instead, with a pseudo-count: a block with thousands of anchors keeps its own estimate,
    // one with a handful is pulled to the region. No cliff, and no tuning of where to put it.
    double global_median = central_value(all_anchor_counts, estimator);
    if (region_stats != nullptr) {
        region_stats->n_anchor = all_anchor_counts.size();
        region_stats->median = median_of(all_anchor_counts);
        region_stats->mean = mean_of(all_anchor_counts);
        region_stats->trimmed_mean = trimmed_mean_of(all_anchor_counts, kTrimFraction);
        region_stats->used = global_median;
    }

    if (model == DepthModel::Quantile || model == DepthModel::Bases) {
        double lambda = 0.0;
        if (model == DepthModel::Quantile) {
            std::vector<double> block_medians;
            for (const BlockDepth& d : out) {
                if (d.n_anchor >= min_anchors && d.median > 0.0) block_medians.push_back(d.median);
            }
            if (!block_medians.empty()) {
                std::sort(block_medians.begin(), block_medians.end());
                const std::size_t k = std::min(block_medians.size() - 1,
                    static_cast<std::size_t>(depth_quantile * static_cast<double>(block_medians.size())));
                lambda = block_medians[k] / 2.0;
            } else if (global_median > 0.0) {
                lambda = global_median / 2.0;
            }
        } else if (region_bp > 0 && counts.reads > 0) {
            // Expected count of a syncmer present once per haplotype: per-haplotype base depth times
            // the share of read positions a k-mer can start at.
            const double read_len = static_cast<double>(counts.bases) / static_cast<double>(counts.reads);
            const double per_hap_depth =
                static_cast<double>(counts.bases) / (2.0 * static_cast<double>(region_bp));
            const double k = static_cast<double>(panel.kmer_size);
            lambda = per_hap_depth * std::max(0.0, (read_len - k + 1.0)) / std::max(1.0, read_len);
        }
        if (lambda > 0.0) {
            const DepthSource src =
                model == DepthModel::Quantile ? DepthSource::Quantile : DepthSource::Bases;
            for (BlockDepth& d : out) {
                d.median = lambda * 2.0;
                d.lambda_hap = lambda;
                d.usable = true;
                d.source = src;
                d.region_shrink_weight = 1.0;   // one region-wide value, no block contributes its own
            }
            return out;
        }
    }

    constexpr double kPseudoAnchors = 200.0;
    for (BlockDepth& d : out) {
        if (global_median <= 0.0) continue;
        if (!d.local_available || !d.usable) {
            d.median = global_median;
            d.lambda_hap = global_median / 2.0;
            d.usable = true;
            d.source = DepthSource::RegionFallback;
            d.region_shrink_weight = 1.0;
            continue;
        }
        const double w = static_cast<double>(d.n_anchor);
        d.median = (w * d.median + kPseudoAnchors * global_median) / (w + kPseudoAnchors);
        d.lambda_hap = d.median / 2.0;
        d.source = DepthSource::Shrunk;
        d.region_shrink_weight = kPseudoAnchors / (w + kPseudoAnchors);
    }
    return out;
}

const char* depth_source_name(DepthSource source) {
    switch (source) {
        case DepthSource::Local:          return "LOCAL";
        case DepthSource::Shrunk:         return "SHRUNK";
        case DepthSource::RegionFallback: return "REGION_FALLBACK";
        case DepthSource::Quantile:       return "QUANTILE";
        case DepthSource::Bases:          return "BASES";
        case DepthSource::Joint:          return "JOINT";
        case DepthSource::None:           break;
    }
    return "NONE";
}

void write_marker_dump(
    const std::string& path,
    const std::vector<Block>& chain,
    const ReadPanel& panel,
    const ReadCounts& counts,
    const std::vector<BlockDepth>& depth,
    const std::vector<std::string>* truth_seq1,
    const std::vector<std::string>* truth_seq2) {

    std::ofstream f(path);
    if (!f) throw std::runtime_error("genotype: cannot write " + path);
    // `marker_class` separates the invariant depth anchors from the informative markers. Comparing the
    // two without it confounds efficiency with dosage: an array marker carried twenty times and a
    // single-copy anchor differ in expected count by twenty-fold before any efficiency effect.
    // `truth_mult` is this sample's own expected copy number when truth is supplied, and
    // `count_per_copy` divides it out, which is the quantity an efficiency model should see as flat.
    f << "block_index\tblock_kind\tbubble_id\tmarker_class\tslot\tcount\tgc\tfirst_pos\tclump"
         "\tn_alleles\tcarriers\tmult_min\tmult_max\tmult_mean\ttruth_mult\tcount_per_copy"
         "\ttruth_pos_mean\ttruth_pos_span"
         "\tvary_blocks\tocc_blocks\tactual\texpected\n";
    const bool have_dbg = panel.dbg_vary.size() == panel.node_codes.size();
    const double kFrag = panel.fragment_len > 0.0 ? panel.fragment_len : 350.0;
    std::unordered_map<std::uint64_t, std::uint32_t> slot_of;
    slot_of.reserve(panel.node_codes.size() * 2);
    for (std::uint32_t i = 0; i < panel.node_codes.size(); ++i) slot_of.emplace(panel.node_codes[i], i);

    // A truth haplotype's own copy number for every marker, read off the spelled sequence rather than
    // off a panel allele index. This is the whole point: at a held-out array no panel allele matches,
    // so an index-based dosage is missing at precisely the block being studied.
    // Positions as well as counts. `node_first_pos` records one occurrence from one allele, which for
    // a repeat marker carried twenty times is one of twenty places it sits -- useless as a positional
    // covariate at exactly the blocks the analysis is about. Taken from the truth sequence, every
    // occurrence is seen, so the mean and span are the sample's own.
    struct TruthPos { double n = 0, sum = 0, lo = 1e30, hi = -1e30; };
    auto dosage_from = [&](const std::string& seq, std::unordered_map<std::uint32_t, double>& out,
                           std::unordered_map<std::uint32_t, TruthPos>& pos_out) {
        if (seq.empty()) return;
        const std::vector<KmerOccurrence> occ =
            panel.all_kmers ? collect_canonical_kmer_occurrences(seq, panel.kmer_size)
                            : collect_syncmers(seq, panel.kmer_size, panel.syncmer_s);
        for (const KmerOccurrence& o : occ) {
            const auto it = slot_of.find(o.code);
            if (it == slot_of.end()) continue;
            out[it->second] += 1.0;
            TruthPos& tp = pos_out[it->second];
            // Normalized within its own sequence. The two haplotypes can differ in length by tens of
            // kilobases at an array, so pooling raw offsets would make position mean different things
            // for each and would correlate with copy number rather than with location.
            const double u = static_cast<double>(o.start) / static_cast<double>(std::max<std::size_t>(1, seq.size()));
            tp.n += 1.0; tp.sum += u;
            tp.lo = std::min(tp.lo, u); tp.hi = std::max(tp.hi, u);
        }
    };

    auto emit = [&](std::size_t bi, const char* kind, const char* cls, std::uint32_t slot,
                    std::size_t n_alleles, std::size_t carriers, std::uint32_t mmin,
                    std::uint32_t mmax, double mmean, double truth_mult, bool truth_known,
                    const TruthPos* tp) {
        if (slot >= counts.node.size() || slot >= panel.node_codes.size()) return;
        const std::string kmer = decode_kmer(panel.node_codes[slot], panel.kmer_size);
        std::size_t gc = 0;
        for (const char c : kmer) if (c == 'G' || c == 'C' || c == 'g' || c == 'c') ++gc;
        const std::uint32_t pos =
            slot < panel.node_first_pos.size() ? panel.node_first_pos[slot] : UINT32_MAX;
        f << chain[bi].index << '\t' << kind << '\t' << chain[bi].bubble_id << '\t' << cls << '\t'
          << slot << '\t' << counts.node[slot] << '\t'
          << (kmer.empty() ? 0.0 : static_cast<double>(gc) / static_cast<double>(kmer.size())) << '\t';
        if (pos == UINT32_MAX) f << "NA\tNA\t";
        else f << pos << '\t' << static_cast<long>(pos / kFrag) << '\t';
        f << n_alleles << '\t' << carriers << '\t' << mmin << '\t' << mmax << '\t' << mmean << '\t';
        // A marker the sample carries zero times has a defined dosage and an UNDEFINED per-copy rate;
        // reporting 0 there would put a background observation into the efficiency population.
        if (!truth_known) f << "NA\tNA\t";
        else if (truth_mult <= 0.0) f << truth_mult << "\tNA\t";
        else f << truth_mult << '\t' << static_cast<double>(counts.node[slot]) / truth_mult << '\t';
        if (tp == nullptr || tp->n <= 0.0) f << "NA\tNA\t";
        else f << (tp->sum / tp->n) << '\t' << (tp->hi - tp->lo) << '\t';
        if (have_dbg) {
            f << panel.dbg_vary[slot] << '\t' << panel.dbg_occ[slot] << '\t'
              << panel.dbg_actual[slot] << '\t' << panel.dbg_expected[slot] << '\n';
        } else {
            f << "NA\tNA\tNA\tNA\n";
        }
    };

    for (std::size_t bi = 0; bi < chain.size(); ++bi) {
        const char* kind = chain[bi].kind == BlockKind::Bubble ? "bubble"
                         : chain[bi].kind == BlockKind::Flank  ? "flank" : "backbone";
        const std::size_t na = bi < panel.by_block.size() ? panel.by_block[bi].size() : 0;

        // Per-slot multiplicity across this block's alleles, and this sample's own copy number.
        std::unordered_map<std::uint32_t, std::array<std::uint64_t, 3>> agg;  // carriers, sum, max
        std::unordered_map<std::uint32_t, std::uint32_t> mmin;
        std::unordered_map<std::uint32_t, double> tmult;
        std::unordered_map<std::uint32_t, TruthPos> tpos;
        // Truth is known for this block when sequences were supplied at all. A haplotype that does not
        // traverse contributes an empty string and therefore zero copies, which is a real dosage of 0
        // rather than an absence of information -- the distinction anchors depend on.
        // Both haplotypes bypassing is a dosage of ZERO, which is known, not unknown. Calling it
        // unknown discards the one case where a marker's absence is the informative observation.
        const bool block_has_truth =
            truth_seq1 != nullptr && truth_seq2 != nullptr &&
            bi < truth_seq1->size() && bi < truth_seq2->size();
        if (block_has_truth) {
            dosage_from((*truth_seq1)[bi], tmult, tpos);
            dosage_from((*truth_seq2)[bi], tmult, tpos);
        }
        if (bi < panel.by_block.size()) {
            for (std::size_t ai = 0; ai < panel.by_block[bi].size(); ++ai) {
                for (const auto& [slot, mult] : panel.by_block[bi][ai].nodes) {
                    auto& a = agg[slot];
                    ++a[0]; a[1] += mult; a[2] = std::max<std::uint64_t>(a[2], mult);
                    const auto it = mmin.find(slot);
                    if (it == mmin.end() || mult < it->second) mmin[slot] = mult;
                }
            }
        }

        auto truth_mult_of = [&](std::uint32_t slot) {
            const auto it = tmult.find(slot);
            return it == tmult.end() ? 0.0 : it->second;
        };
        auto truth_pos_of = [&](std::uint32_t slot) -> const TruthPos* {
            const auto it = tpos.find(slot);
            return it == tpos.end() ? nullptr : &it->second;
        };
        if (bi < panel.anchor_slots.size()) {
            for (const std::uint32_t slot : panel.anchor_slots[bi]) {
                // An anchor is carried once by every TRAVERSING allele, so a diploid expects two copies
                // only when both haplotypes traverse. Hard-coding 2 was wrong wherever a truth
                // haplotype takes a bypass or deletion allele, which is the case a depth anchor is
                // least able to absorb.
                emit(bi, kind, "anchor", slot, na, na, 1, 1, 1.0,
                     truth_mult_of(slot), block_has_truth, truth_pos_of(slot));
            }
        }
        for (const auto& [slot, a] : agg) {
            const double mean = static_cast<double>(a[1]) / static_cast<double>(std::max<std::uint64_t>(1, a[0]));
            emit(bi, kind, "informative", slot, na, a[0], mmin[slot],
                 static_cast<std::uint32_t>(a[2]), mean, truth_mult_of(slot), block_has_truth,
                 truth_pos_of(slot));
        }
    }
    if (!f) throw std::runtime_error("genotype: write failed for " + path);
}

void write_read_audit(
    const std::string& out_prefix,
    const std::vector<Block>& chain,
    const ReadPanel& panel,
    const ReadCounts& counts,
    const std::vector<BlockDepth>& depth) {

    const std::string path = out_prefix + ".reads.depth.tsv";
    std::ofstream f(path);
    if (!f) throw std::runtime_error("genotype: cannot write " + path);
    // A statistic that was never computed is NA, not 0. Zero is a legitimate anchor count, and writing
    // it for "no local estimate" is how an inherited depth came to look like a measured one.
    auto na = [](double v, bool available) {
        return available ? std::to_string(v) : std::string("NA");
    };
    // `raw_anchor_*` are this block's OWN observations and are 0 when it had too few anchors to make
    // any; `fitted_median` is what the model settled on. Writing the fitted value under the raw name
    // is how a block that inherited the region's depth came to look like a precise local measurement,
    // MAD 0 and all, at exactly the blocks where depth IS the answer.
    f << "block_index\tblock_kind\tbubble_id\tn_alleles\tn_anchor\traw_anchor_median\traw_anchor_mad"
         "\tlocal_center\tfitted_median\tlambda_hap\tdepth_source\tregion_shrink_weight"
         "\tlow_anchor\tusable\tuneven\n";
    for (std::size_t bi = 0; bi < chain.size() && bi < depth.size(); ++bi) {
        const BlockDepth& d = depth[bi];
        f << chain[bi].index << '\t' << (chain[bi].kind == BlockKind::Bubble ? "bubble" : chain[bi].kind == BlockKind::Flank ? "flank" : "backbone")
          << '\t' << chain[bi].bubble_id << '\t'
          << (bi < panel.by_block.size() ? panel.by_block[bi].size() : 0) << '\t' << d.n_anchor << '\t'
          << na(d.anchor_median, d.local_available) << '\t' << na(d.mad, d.local_available) << '\t'
          << na(d.local_center, d.local_available) << '\t' << d.median << '\t' << d.lambda_hap << '\t'
          << depth_source_name(d.source) << '\t' << d.region_shrink_weight << '\t'
          << (d.low_anchor ? 1 : 0) << '\t' << (d.usable ? 1 : 0) << '\t' << (d.uneven ? 1 : 0) << '\n';
    }
    (void)counts;
}

} // namespace panvar
