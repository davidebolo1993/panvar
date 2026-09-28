#include "panvar/library_estimate.hpp"

#include "panvar/marker_model.hpp"
#include "panvar/mosaic_spelling.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <sstream>
#include <thread>
#include <unordered_map>

namespace panvar {

namespace {

double seconds_since(std::chrono::steady_clock::time_point t) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count();
}

std::uint64_t fmix64(std::uint64_t h) {
    h ^= h >> 33; h *= 0xff51afd7ed558ccdULL; h ^= h >> 33; h *= 0xc4ceb9fe1a85ec53ULL; h ^= h >> 33;
    return h;
}

std::uint64_t fnv1a64(const std::string& s) {
    std::uint64_t h = 1469598103934665603ULL;
    for (unsigned char c : s) { h ^= c; h *= 1099511628211ULL; }
    return fmix64(h);
}

int code_of(char c) {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
    }
    return -1;
}

// Every k-mer of s with its start, the forward and canonical 2-bit codes; k-mers holding a non-ACGT base are skipped.
struct Kmer { std::size_t at; std::uint64_t forward, canonical; };
std::vector<Kmer> kmers_of(const std::string& s, std::size_t k) {
    std::vector<Kmer> out;
    if (s.size() < k) return out;
    const std::uint64_t mask = k >= 32 ? ~0ULL : ((1ULL << (2 * k)) - 1);
    std::uint64_t f = 0, r = 0;
    std::size_t valid = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        const int c = code_of(s[i]);
        if (c < 0) { valid = 0; f = r = 0; continue; }
        f = ((f << 2) | static_cast<std::uint64_t>(c)) & mask;
        r = (r >> 2) | (static_cast<std::uint64_t>(3 - c) << (2 * (k - 1)));
        if (++valid >= k) out.push_back({i + 1 - k, f, std::min(f, r)});
    }
    return out;
}

std::string reverse_complement(const std::string& s) {
    std::string out(s.rbegin(), s.rend());
    for (char& c : out) {
        switch (c) {
            case 'A': c = 'T'; break; case 'C': c = 'G'; break; case 'G': c = 'C'; break; case 'T': c = 'A'; break;
            default: c = 'N';
        }
    }
    return out;
}

std::string strip_mate_suffix(const std::string& name) {
    if (name.size() > 2 && name[name.size() - 2] == '/' && (name.back() == '1' || name.back() == '2')) return name.substr(0, name.size() - 2);
    return name;
}

double quantile_sorted(const std::vector<double>& v, double p) {
    if (v.empty()) return std::nan("");
    const double x = p * static_cast<double>(v.size() - 1);
    const std::size_t lo = static_cast<std::size_t>(std::floor(x)), hi = std::min(lo + 1, v.size() - 1);
    return v[lo] + (x - static_cast<double>(lo)) * (v[hi] - v[lo]);
}

std::string num(double v, int precision = 6) {
    if (std::isnan(v)) return "NA";
    std::ostringstream o;
    o << std::setprecision(precision) << v;
    return o.str();
}

// One mate's placement: orientation relative to the template frame and its implied start on every template.
struct MatePlacement { bool placed = false, consistent = true, forward = true; std::vector<long> start; };
// One pair's outcome, filled independently per pair (thread-safe), reduced in sample order.
struct PairOutcome {
    int state = 0;   // 0 unplaced, 1 mate inconsistent, 2 orientation discordant, 3 not decisive, 4 accepted
    long insert = 0;
    bool fwd_first = true;
    std::vector<long> fwd_start;   // the forward mate's start per template (for the mismatch count)
    std::vector<long> rev_start;
};

} // namespace

LibraryEstimate estimate_library(const PreparedGenotypeModel& model, const std::vector<ReadPair>& reads,
                                 const LibraryEstimateOptions& options, const LibraryModel* declared) {
    const auto started = std::chrono::steady_clock::now();
    LibraryEstimate est;
    const std::size_t k = options.k;
    // ---- anchor index over the panel templates (canonical block-chain frame)
    const std::size_t T = model.panel().template_count();
    est.templates = T;
    std::vector<std::string> templates(T);
    for (std::size_t t = 0; t < T; ++t) templates[t] = spell_panel_template(model, t).bases;
    struct Anchor { bool forward = true; std::vector<std::uint32_t> pos; };
    std::unordered_map<std::uint64_t, Anchor> anchors;
    {
        std::unordered_map<std::uint64_t, std::pair<std::uint32_t, int>> first;   // canonical -> (pos, count) in template 0
        std::unordered_map<std::uint64_t, bool> strand;
        for (const Kmer& x : kmers_of(templates.empty() ? std::string() : templates[0], k)) {
            if (fmix64(x.canonical) % options.anchor_sampling != 0) continue;
            auto& e = first[x.canonical];
            e.first = static_cast<std::uint32_t>(x.at);
            ++e.second;
            strand[x.canonical] = x.forward == x.canonical;
        }
        for (const auto& [code, e] : first) {
            if (e.second != 1) continue;
            Anchor a; a.forward = strand[code]; a.pos.assign(T, 0); a.pos[0] = e.first;
            anchors.emplace(code, std::move(a));
        }
        for (std::size_t t = 1; t < T && !anchors.empty(); ++t) {
            std::unordered_map<std::uint64_t, int> seen;
            for (const Kmer& x : kmers_of(templates[t], k)) {
                const auto it = anchors.find(x.canonical);
                if (it == anchors.end()) continue;
                int& n = seen[x.canonical];
                ++n;
                if ((x.forward == x.canonical) != it->second.forward) n += 100;   // other strand: not an anchor
                it->second.pos[t] = static_cast<std::uint32_t>(x.at);
            }
            for (auto it = anchors.begin(); it != anchors.end();) {
                const auto s = seen.find(it->first);
                if (s == seen.end() || s->second != 1) it = anchors.erase(it);
                else ++it;
            }
        }
    }
    est.anchors = anchors.size();
    est.seconds_index = seconds_since(started);

    // ---- read-level quantities over every pair
    est.pairs_total = reads.size();
    std::vector<double> len1, len2;
    std::uint64_t bases = 0, ambiguous = 0, q30 = 0;
    double qsum = 0.0, perr = 0.0;
    std::uint64_t qbases = 0;
    for (const ReadPair& p : reads) {
        len1.push_back(static_cast<double>(p.first.bases.size()));
        len2.push_back(static_cast<double>(p.second.bases.size()));
        for (const SequencingRead* r : {&p.first, &p.second}) {
            for (char c : r->bases) { ++bases; if (code_of(c) < 0) ++ambiguous; }
            for (char q : r->qualities) {
                const int phred = static_cast<int>(q) - 33;
                qsum += phred; perr += std::pow(10.0, -phred / 10.0); ++qbases; if (phred >= 30) ++q30;
            }
        }
    }
    std::sort(len1.begin(), len1.end());
    std::sort(len2.begin(), len2.end());
    est.read1_length_median = quantile_sorted(len1, 0.5);
    est.read2_length_median = quantile_sorted(len2, 0.5);
    est.phred_error_mean = qbases ? perr / static_cast<double>(qbases) : std::nan("");

    // ---- deterministic sample: the pairs with the smallest salted name hashes, in hash order
    std::vector<std::pair<std::uint64_t, std::size_t>> order;
    order.reserve(reads.size());
    for (std::size_t i = 0; i < reads.size(); ++i) order.push_back({fnv1a64(options.salt + "/" + strip_mate_suffix(reads[i].first.name)), i});
    std::sort(order.begin(), order.end());
    if (order.size() > options.max_pairs) order.resize(options.max_pairs);
    est.pairs_sampled = order.size();

    // ---- placement, per pair, in parallel over contiguous chunks of the sample
    const auto t_place = std::chrono::steady_clock::now();
    const auto place_mate = [&](const std::string& bases) {
        MatePlacement m;
        const std::size_t L = bases.size();
        for (const Kmer& x : kmers_of(bases, k)) {
            const auto it = anchors.find(x.canonical);
            if (it == anchors.end()) continue;
            const bool read_forward = x.forward == x.canonical;
            const bool forward = read_forward == it->second.forward;
            std::vector<long> start(T);
            for (std::size_t t = 0; t < T; ++t) {
                start[t] = forward ? static_cast<long>(it->second.pos[t]) - static_cast<long>(x.at)
                                   : static_cast<long>(it->second.pos[t]) - static_cast<long>(L - x.at - k);
            }
            if (!m.placed) { m.placed = true; m.forward = forward; m.start = start; }
            else if (m.forward != forward || m.start != start) m.consistent = false;
        }
        return m;
    };
    std::vector<PairOutcome> outcomes(order.size());
    // anchor k-mer depth over EVERY pair (the rate per copy), counted per thread and summed in thread order
    std::vector<std::uint64_t> anchor_ids;
    anchor_ids.reserve(anchors.size());
    for (const auto& [code, a] : anchors) anchor_ids.push_back(code);
    std::sort(anchor_ids.begin(), anchor_ids.end());
    std::unordered_map<std::uint64_t, std::size_t> anchor_index;
    for (std::size_t i = 0; i < anchor_ids.size(); ++i) anchor_index[anchor_ids[i]] = i;
    const std::size_t workers = std::max<std::size_t>(1, options.threads);
    std::vector<std::vector<std::uint32_t>> depth(workers, std::vector<std::uint32_t>(anchor_ids.size(), 0));
    {
        std::vector<std::thread> pool;
        for (std::size_t w = 0; w < workers; ++w) {
            pool.emplace_back([&, w]() {
                const std::size_t lo = order.size() * w / workers, hi = order.size() * (w + 1) / workers;
                for (std::size_t s = lo; s < hi; ++s) {
                    const ReadPair& p = reads[order[s].second];
                    PairOutcome& o = outcomes[s];
                    const MatePlacement a = place_mate(p.first.bases), b = place_mate(p.second.bases);
                    if (!a.placed || !b.placed) continue;
                    if (!a.consistent || !b.consistent) { o.state = 1; continue; }
                    if (a.forward == b.forward) { o.state = 2; continue; }
                    const MatePlacement& f = a.forward ? a : b;
                    const MatePlacement& r = a.forward ? b : a;
                    const long rl = static_cast<long>((a.forward ? p.second : p.first).bases.size());
                    std::vector<long> insert(T);
                    for (std::size_t t = 0; t < T; ++t) insert[t] = r.start[t] + rl - f.start[t];
                    if (insert[0] <= 0 || r.start[0] < f.start[0]) { o.state = 2; continue; }   // outward-facing
                    if (std::any_of(insert.begin(), insert.end(), [&](long v) { return v != insert[0]; })) { o.state = 3; continue; }
                    o.state = 4; o.insert = insert[0]; o.fwd_first = a.forward; o.fwd_start = f.start; o.rev_start = r.start;
                }
                const std::size_t dlo = reads.size() * w / workers, dhi = reads.size() * (w + 1) / workers;
                for (std::size_t i = dlo; i < dhi; ++i) {
                    for (const SequencingRead* m : {&reads[i].first, &reads[i].second}) {
                        for (const Kmer& x : kmers_of(m->bases, k)) {
                            const auto it = anchor_index.find(x.canonical);
                            if (it != anchor_index.end()) ++depth[w][it->second];
                        }
                    }
                }
            });
        }
        for (std::thread& t : pool) t.join();
    }
    std::vector<double> inserts;
    for (const PairOutcome& o : outcomes) {
        if (o.state >= 1) ++est.both_anchored;
        if (o.state == 1) ++est.mate_inconsistent;
        if (o.state == 2) ++est.orientation_discordant;
        if (o.state == 3) ++est.not_decisive;
        if (o.state >= 3) ++est.concordant;
        if (o.state == 4) { ++est.accepted; est.accepted_inserts.push_back(o.insert); inserts.push_back(static_cast<double>(o.insert)); }
    }
    // invariant-site mismatch rate over the first accepted pairs in sample order
    std::uint64_t sites = 0, mismatches = 0;
    std::size_t used = 0;
    for (std::size_t s = 0; s < outcomes.size() && used < options.error_pairs; ++s) {
        const PairOutcome& o = outcomes[s];
        if (o.state != 4) continue;
        ++used;
        const ReadPair& p = reads[order[s].second];
        const std::string fwd = o.fwd_first ? p.first.bases : p.second.bases;
        const std::string rev = reverse_complement(o.fwd_first ? p.second.bases : p.first.bases);
        for (const auto& [seq, start] : {std::make_pair(&fwd, &o.fwd_start), std::make_pair(&rev, &o.rev_start)}) {
            for (std::size_t j = 0; j < seq->size(); ++j) {
                if (code_of((*seq)[j]) < 0) continue;
                char base = 0;
                bool invariant = true;
                for (std::size_t t = 0; t < T && invariant; ++t) {
                    const long at = (*start)[t] + static_cast<long>(j);
                    if (at < 0 || at >= static_cast<long>(templates[t].size())) { invariant = false; break; }
                    const char c = templates[t][static_cast<std::size_t>(at)];
                    if (t == 0) base = c; else if (c != base) invariant = false;
                }
                if (!invariant) continue;
                ++sites;
                if ((*seq)[j] != base) ++mismatches;
            }
        }
    }
    est.invariant_mismatch_rate = sites ? static_cast<double>(mismatches) / static_cast<double>(sites) : std::nan("");
    est.seconds_placement = seconds_since(t_place);

    // ---- robust insert summary
    std::vector<double> sorted = inserts;
    std::sort(sorted.begin(), sorted.end());
    est.insert_median = quantile_sorted(sorted, 0.5);
    std::vector<double> dev;
    for (double x : sorted) dev.push_back(std::fabs(x - est.insert_median));
    std::sort(dev.begin(), dev.end());
    est.insert_robust_sd = 1.4826 * quantile_sorted(dev, 0.5);
    std::vector<double> core;
    for (double x : sorted) {
        if (est.insert_robust_sd > 0 && std::fabs(x - est.insert_median) > options.outlier_sds * est.insert_robust_sd) ++est.outliers;
        else core.push_back(x);
    }
    double m = 0.0, v = 0.0;
    for (double x : core) m += x;
    m = core.empty() ? std::nan("") : m / static_cast<double>(core.size());
    for (double x : core) v += (x - m) * (x - m);
    est.insert_mean = m;
    est.insert_sd = core.size() > 1 ? std::sqrt(v / static_cast<double>(core.size() - 1)) : std::nan("");
    // two-group split (Otsu over the core): separation in within-group sd
    double best_sep = 0.0, best_minor = 0.0;
    if (core.size() >= 20) {
        std::vector<double> prefix(core.size() + 1, 0.0), prefix2(core.size() + 1, 0.0);
        for (std::size_t i = 0; i < core.size(); ++i) { prefix[i + 1] = prefix[i] + core[i]; prefix2[i + 1] = prefix2[i] + core[i] * core[i]; }
        double best_between = -1.0;
        for (std::size_t i = 1; i < core.size(); ++i) {
            const double n1 = static_cast<double>(i), n2 = static_cast<double>(core.size() - i);
            const double m1 = prefix[i] / n1, m2 = (prefix.back() - prefix[i]) / n2;
            const double between = n1 * n2 * (m1 - m2) * (m1 - m2);
            if (between > best_between) {
                best_between = between;
                const double v1 = prefix2[i] / n1 - m1 * m1, v2 = (prefix2.back() - prefix2[i]) / n2 - m2 * m2;
                const double sd = std::sqrt(std::max(0.0, std::max(v1, v2)));
                best_sep = sd > 0 ? std::fabs(m2 - m1) / sd : 0.0;
                best_minor = std::min(n1, n2) / static_cast<double>(core.size());
            }
        }
    }
    est.multimodal = best_sep >= options.multimodal_separation_sds && best_minor >= options.multimodal_min_fraction;

    // ---- anchor depth and fragment-start intensity
    std::vector<double> per_anchor(anchor_ids.size(), 0.0);
    for (std::size_t i = 0; i < anchor_ids.size(); ++i) {
        std::uint64_t n = 0;
        for (std::size_t w = 0; w < workers; ++w) n += depth[w][i];
        per_anchor[i] = static_cast<double>(n);
    }
    std::sort(per_anchor.begin(), per_anchor.end());
    // an anchor occurs once per haplotype: two copies. The 5%-trimmed mean, not the median: counts are small integers at
    // low coverage, and the trim drops the locus-end and repeat-edge anchors
    const std::size_t trim = per_anchor.size() / 20;
    double depth_sum = 0.0;
    for (std::size_t i = trim; i + trim < per_anchor.size(); ++i) depth_sum += per_anchor[i];
    const std::size_t depth_n = per_anchor.size() > 2 * trim ? per_anchor.size() - 2 * trim : 0;
    est.anchor_depth_per_copy = depth_n ? depth_sum / static_cast<double>(depth_n) / 2.0 : std::nan("");
    const double depth_median = quantile_sorted(per_anchor, 0.5) / 2.0;
    const double read_len = (est.read1_length_median + est.read2_length_median) / 2.0;
    const double eps = std::isnan(est.invariant_mismatch_rate) ? 0.0 : est.invariant_mismatch_rate;
    const double window = read_len - static_cast<double>(k) + 1.0;
    // each fragment start on a haplotype yields two mates; a k-mer is read error-free with probability (1-eps)^k
    est.fragment_start_intensity = window > 0 ? est.anchor_depth_per_copy / (2.0 * window * std::pow(1.0 - eps, static_cast<double>(k))) : std::nan("");

    // the caller's own marker depth fit (invariant anchors of known copy number, Poisson MLE with background), with the
    // same read-window conversion: the second, independent route to the fragment-start intensity
    double fit_rate = std::nan(""), fit_intensity = std::nan("");
    std::string fit_note = "none";
    try {
        const PreparedMarkerPanel panel = build_marker_panel(model, options.marker_kmer, options.marker_syncmer);
        const MarkerDepthEstimate fit = estimate_marker_depth(panel, count_marker_observations(panel, reads), std::nullopt);
        if (fit.converged) {
            fit_rate = fit.rate_per_copy;
            const double w = read_len - static_cast<double>(options.marker_kmer) + 1.0;
            fit_intensity = w > 0 ? fit_rate / (2.0 * w * std::pow(1.0 - eps, static_cast<double>(options.marker_kmer))) : std::nan("");
        } else {
            fit_note = "marker_depth_not_converged:" + fit.termination;
        }
    } catch (const std::exception& error) {
        fit_note = std::string("marker_depth_failed:") + error.what();
    }
    est.marker_fit_intensity = fit_intensity;

    // ---- rows
    const std::size_t N = est.pairs_total, S = est.pairs_sampled, A = est.accepted;
    const auto add = [&](const std::string& p, const std::string& value, const std::string& method, std::size_t examined,
                         std::size_t accepted, const std::string& spread, const std::string& status, const std::string& reason,
                         const std::string& source) {
        est.rows.push_back({p, value, method, examined, accepted, spread, status, reason, source});
    };
    add("pairs_total", std::to_string(N), "count", N, N, "NA", "MEASURED", "none", "reads");
    add("read1_length_median", num(est.read1_length_median), "median", N, N,
        "min=" + num(len1.empty() ? std::nan("") : len1.front()) + ",max=" + num(len1.empty() ? std::nan("") : len1.back()), "MEASURED", "none", "reads");
    add("read2_length_median", num(est.read2_length_median), "median", N, N,
        "min=" + num(len2.empty() ? std::nan("") : len2.front()) + ",max=" + num(len2.empty() ? std::nan("") : len2.back()), "MEASURED", "none", "reads");
    add("base_quality_mean", num(qbases ? qsum / static_cast<double>(qbases) : std::nan("")), "mean_phred", N, N,
        "fraction_q30=" + num(qbases ? static_cast<double>(q30) / static_cast<double>(qbases) : std::nan("")), "MEASURED", "none", "reads");
    add("phred_error_mean", num(est.phred_error_mean), "mean_of_10^(-Q/10)", N, N, "NA", "MEASURED", "none", "reads");
    add("ambiguous_base_fraction", num(bases ? static_cast<double>(ambiguous) / static_cast<double>(bases) : std::nan("")), "fraction_non_ACGT",
        N, N, "NA", "MEASURED", "none", "reads");
    add("pairs_sampled", std::to_string(S), "smallest_salted_name_hash:" + options.salt, N, S, "NA", "MEASURED", "none", "reads");
    add("anchors", std::to_string(est.anchors), "kmers_once_per_template_one_strand,k=" + std::to_string(k) + ",hash_mod=" +
        std::to_string(options.anchor_sampling), 0, 0, "templates=" + std::to_string(T), "MEASURED", "none", "model");
    add("pairs_both_anchored", std::to_string(est.both_anchored), "anchor_placement", S, est.both_anchored, "NA", "MEASURED", "none", "placement_inference");
    add("pairs_mate_inconsistent", std::to_string(est.mate_inconsistent), "anchor_placement", S, 0, "NA", "MEASURED", "none", "placement_inference");
    add("pairs_orientation_discordant", std::to_string(est.orientation_discordant), "not_FR_inward", S, 0,
        "fraction_of_anchored=" + num(est.both_anchored ? static_cast<double>(est.orientation_discordant) / static_cast<double>(est.both_anchored) : std::nan("")),
        "MEASURED", "none", "placement_inference");
    add("pairs_not_decisive", std::to_string(est.not_decisive), "insert_differs_between_templates", S, 0, "NA", "MEASURED", "none", "placement_inference");
    add("pairs_accepted", std::to_string(A), "FR_concordant_decisive", S, A,
        "fraction_of_sampled=" + num(S ? static_cast<double>(A) / static_cast<double>(S) : std::nan("")), "MEASURED", "none", "placement_inference");
    const bool enough = A >= options.min_accepted_pairs;
    const double outlier_fraction = A ? static_cast<double>(est.outliers) / static_cast<double>(A) : std::nan("");
    std::string status = "ESTIMATED", reason = "none";
    if (!enough) { status = "FALLBACK"; reason = "insufficient_placements(" + std::to_string(A) + "<" + std::to_string(options.min_accepted_pairs) + ")"; }
    else if (est.multimodal) { status = "FALLBACK"; reason = "multimodal(separation=" + num(best_sep, 3) + "sd,minor=" + num(best_minor, 3) + ")"; }
    else if (outlier_fraction > options.max_outlier_fraction) { status = "WARN"; reason = "outlier_fraction=" + num(outlier_fraction, 3); }
    const bool fallback = status == "FALLBACK";
    const std::string fb_source = declared ? "explicit_or_preset_fallback" : "none_available";
    add("insert_mean_bp", fallback ? (declared ? num(declared->insert_mean_bp) : "NA") : num(est.insert_mean),
        "mean_within_" + num(options.outlier_sds, 2) + "_robust_sd_of_median", S, A, "median=" + num(est.insert_median) + ",sd=" + num(est.insert_sd),
        status, reason, fallback ? fb_source : "placement_inference");
    add("insert_sd_bp", fallback ? (declared ? num(declared->insert_sd_bp) : "NA") : num(est.insert_sd), "sd_within_outlier_bound", S, A,
        "mad_sd=" + num(est.insert_robust_sd) + ",q05=" + num(quantile_sorted(sorted, 0.05)) + ",q95=" + num(quantile_sorted(sorted, 0.95)),
        status, reason, fallback ? fb_source : "placement_inference");
    add("insert_median_bp", num(est.insert_median), "median", S, A, "robust_sd=" + num(est.insert_robust_sd), enough ? "MEASURED" : "FALLBACK",
        enough ? "none" : reason, "placement_inference");
    add("insert_outlier_fraction", num(outlier_fraction), "beyond_" + num(options.outlier_sds, 2) + "_robust_sd", S, A, "NA",
        outlier_fraction > options.max_outlier_fraction ? "WARN" : "MEASURED", "none", "placement_inference");
    add("insert_bimodality_separation_sd", num(best_sep, 4), "otsu_split_within_group_sd", S, A, "minor_fraction=" + num(best_minor, 4),
        est.multimodal ? "WARN" : "MEASURED", est.multimodal ? "multimodal" : "none", "placement_inference");
    add("substitution_error", num(est.invariant_mismatch_rate), "invariant_site_mismatch_rate", std::min(A, options.error_pairs),
        std::min(A, options.error_pairs), "sites=" + std::to_string(sites) + ",mismatches=" + std::to_string(mismatches) + ",phred_implied=" + num(est.phred_error_mean),
        sites ? "DIAGNOSTIC" : "FALLBACK", sites ? "preset_parameter_not_replaced" : "no_invariant_sites", "placement_inference");
    add("anchor_depth_per_copy", num(est.anchor_depth_per_copy), "trimmed_mean_anchor_kmer_count/2", N, N, "median=" + num(depth_median), "MEASURED",
        "none", "placement_inference");
    add("fragment_start_intensity", num(est.fragment_start_intensity), "anchor_depth/(2*(L-k+1)*(1-e)^k)", N, N, "NA",
        est.anchors ? "ESTIMATED" : "FALLBACK", est.anchors ? "none" : "no_anchors", "placement_inference");
    add("marker_depth_rate_per_copy", num(fit_rate), "caller_marker_depth_fit(k=" + std::to_string(options.marker_kmer) + ",s=" +
        std::to_string(options.marker_syncmer) + ")", N, N, "NA", std::isnan(fit_rate) ? "FALLBACK" : "MEASURED", fit_note, "placement_inference");
    add("fragment_start_intensity_marker_fit", num(fit_intensity), "marker_rate/(2*(L-k+1)*(1-e)^k)", N, N, "NA",
        std::isnan(fit_intensity) ? "FALLBACK" : "ESTIMATED", fit_note, "placement_inference");
    add("support_sigmas", declared ? num(declared->support_sigmas) : "NA", "not_estimated", 0, 0, "NA", "PRESET", "model_parameter", declared ? "explicit_or_preset" : "none_available");
    add("outlier_weight", declared ? num(declared->outlier_weight) : "NA", "not_estimated", 0, 0, "NA", "PRESET", "model_parameter", declared ? "explicit_or_preset" : "none_available");
    est.seconds_total = seconds_since(started);
    add("seconds_total", num(est.seconds_total), "wall", N, N, "index=" + num(est.seconds_index) + ",placement=" + num(est.seconds_placement), "MEASURED", "none", "timing");
    return est;
}

std::string library_estimate_table(const LibraryEstimate& estimate, const LibraryModel* declared) {
    const auto declared_value = [&](const std::string& p) -> std::string {
        if (!declared) return "NA";
        if (p == "insert_mean_bp") return num(declared->insert_mean_bp);
        if (p == "insert_sd_bp") return num(declared->insert_sd_bp);
        if (p == "insert_median_bp") return num(declared->insert_mean_bp);
        if (p == "substitution_error") return num(declared->substitution_error);
        if (p == "fragment_start_intensity" || p == "fragment_start_intensity_marker_fit") return num(declared->fragment_start_intensity);
        if (p == "support_sigmas") return num(declared->support_sigmas);
        if (p == "outlier_weight") return num(declared->outlier_weight);
        return "NA";
    };
    std::ostringstream o;
    o << "parameter\tvalue\tmethod\tpairs_examined\tpairs_accepted\tspread\tstatus\treason\tsource\tdeclared_value\n";
    for (const LibraryEstimateRow& r : estimate.rows) {
        o << r.parameter << '\t' << r.value << '\t' << r.method << '\t' << r.pairs_examined << '\t' << r.pairs_accepted << '\t' << r.spread
          << '\t' << r.status << '\t' << r.reason << '\t' << r.source << '\t' << declared_value(r.parameter) << '\n';
    }
    return o.str();
}

} // namespace panvar
