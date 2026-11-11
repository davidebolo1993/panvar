#include "panvar/recruit_reads.hpp"

#include "panvar/flat_map.hpp"
#include "panvar/kmer.hpp"
#include "panvar/read_input.hpp"
#include "panvar/recruit_calibration.hpp"
#include "panvar/recruit_scoring.hpp"

#include <zlib.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <thread>

namespace panvar::recruit {

using namespace scoring;


void Evidence::add(const Evidence& o) {
    sites_p += o.sites_p; sites_d += o.sites_d; decoy_out += o.decoy_out; runs_n += o.runs_n; seg_kmers += o.seg_kmers;
    n_p += o.n_p; n_s += o.n_s; n_r += o.n_r; n_d += o.n_d; n_n += o.n_n;
    kmers += o.kmers; windows += o.windows; bases += o.bases; multi += o.multi;
}

namespace {

// Recruited reads as gzip FASTQ, per locus and fate: <locus>.<fate>_1/_2.fq.gz for pairs, <locus>.<fate>.fq.gz
// for single reads. Files open when first written; past a cap they are closed and reopened to append (a
// gzip file of several members reads as one).
class ReadWriters {
public:
    explicit ReadWriters(std::string dir) : dir_(std::move(dir)) { std::filesystem::create_directories(dir_); }
    ~ReadWriters() { close_all(); }
    void write(const std::string& locus, const char* fate, const char* suffix, const std::string& header, const std::string& seq,
               const std::string& qual) {
        const std::string key = locus + "." + fate + suffix + ".fq.gz";
        auto it = open_.find(key);
        if (it == open_.end()) {
            if (open_.size() >= kMaxOpen) close_all();
            const bool again = started_.count(key) > 0;
            gzFile f = gzopen((std::filesystem::path(dir_) / key).string().c_str(), again ? "ab3" : "wb3");
            if (f == nullptr) throw std::runtime_error("recruit: cannot write " + dir_ + "/" + key);
            started_.insert(key);
            it = open_.emplace(key, f).first;
        }
        std::string rec = "@" + header + "\n" + seq + "\n+\n";
        rec += qual.size() == seq.size() ? qual : std::string(seq.size(), 'I');
        rec += '\n';
        if (gzwrite(it->second, rec.data(), static_cast<unsigned>(rec.size())) != static_cast<int>(rec.size()))
            throw std::runtime_error("recruit: cannot write " + dir_ + "/" + key);
    }
    void close_all() {
        for (auto& [key, f] : open_)
            if (gzclose(f) != Z_OK) throw std::runtime_error("recruit: cannot write " + dir_ + "/" + key);
        open_.clear();
    }

private:
    static constexpr std::size_t kMaxOpen = 200;
    std::string dir_;
    std::map<std::string, gzFile> open_;
    std::set<std::string> started_;
};

// The error rate from base qualities: the mean of 10^(-Q/10) over the bases that carry one, within [0.001,
// 0.08]; without qualities, or with qualities that say nothing (a mean above 30%, as from simulators that
// write Q0), 0.03 for long reads and 0.002 for short ones.
double estimated_error(const std::vector<Fragment>& batch, std::string& source) {
    double sum = 0;
    std::uint64_t n = 0, bases = 0, reads = 0;
    std::array<double, 94> p{};
    for (int q = 0; q < 94; ++q) p[q] = std::pow(10.0, -q / 10.0);
    for (const Fragment& f : batch)
        for (const std::string* q : {&f.qual1, &f.qual2}) {
            for (const char c : *q) {
                const int v = std::clamp(static_cast<int>(c) - 33, 0, 93);
                sum += p[v];
                ++n;
            }
        }
    for (const Fragment& f : batch) {
        bases += f.seq1.size() + f.seq2.size();
        reads += 1 + !f.seq2.empty();
    }
    const double raw = n ? sum / static_cast<double>(n) : 1.0;
    if (n >= 1000 && raw <= 0.3) {
        source = "estimated from " + std::to_string(n) + " base qualities";
        return std::clamp(raw, 0.001, 0.08);
    }
    const bool long_reads = reads > 0 && bases / reads > 1000;
    source = std::string(n >= 1000 ? "base qualities uninformative (all near Q0)" : "no base qualities") + ": the default for " +
             (long_reads ? "long" : "short") + " reads";
    return long_reads ? 0.03 : 0.002;
}

}  // namespace

ReadStats recruit_reads(const ReadOptions& options) {
    ReadOptions o = options;   // its error rate may be estimated below
    auto log = [&](const std::string& s) { if (o.log) o.log(s); };
    if (o.indexes.empty()) throw std::runtime_error("recruit: no locus index");
    std::vector<LocusTable> loci(o.indexes.size());
    std::vector<std::vector<Calibration>> calibration(o.indexes.size());
    std::vector<std::pair<std::uint64_t, std::uint32_t>> screen_items;
    for (std::size_t li = 0; li < o.indexes.size(); ++li) {
        const LocusIndex x = load(o.indexes[li]);
        loci[li] = table_of(x);
        calibration[li] = x.calibration;
        for (std::size_t lj = 0; lj < li; ++lj)
            if (loci[lj].name == x.name)
                throw std::runtime_error("recruit: two loci are named " + x.name + " (" + o.indexes[lj] + ", " + o.indexes[li] + "); rename one graph");
        for (const std::uint64_t s : x.screen) screen_items.emplace_back(s, static_cast<std::uint32_t>(li));
        log(x.name + ": " + std::to_string(x.codes.size()) + " k-mers, " + std::to_string(x.decoy.size()) + " decoys, " +
            std::to_string(x.screen.size()) + " screen syncmers");
    }
    // the screen: syncmer -> an entry of `screen_loci`, the loci it belongs to (any number of loci)
    std::sort(screen_items.begin(), screen_items.end());
    std::vector<std::vector<std::uint32_t>> screen_loci;
    std::map<std::vector<std::uint32_t>, std::uint32_t> set_id;   // few distinct sets: loci rarely share syncmers
    FlatMap<std::uint32_t> screen;
    screen.reserve(screen_items.size());
    for (std::size_t i = 0; i < screen_items.size();) {
        std::size_t j = i;
        std::vector<std::uint32_t> ls;
        while (j < screen_items.size() && screen_items[j].first == screen_items[i].first) ls.push_back(screen_items[j++].second);
        ls.erase(std::unique(ls.begin(), ls.end()), ls.end());
        const auto [it, fresh] = set_id.emplace(ls, static_cast<std::uint32_t>(screen_loci.size()));
        if (fresh) screen_loci.push_back(std::move(ls));
        bool existing = false;
        screen.slot(screen_items[i].first, existing) = it->second;
        i = j;
    }
    screen_items = {};

    // whether a read hits some locus's screen (the input layer asks, for coordinate-sorted files)
    auto screened = [&](const std::string& seq) {
        bool hit = false;
        kmer::for_each_closed_syncmer(seq.data(), seq.size(), kK, kS, [&](std::uint64_t c, std::size_t, bool) {
            if (!hit && screen.find(c) != nullptr) hit = true;
        });
        return hit;
    };
    if (o.inputs.empty()) throw std::runtime_error("recruit: no reads given");
    InputOptions io;
    io.cram_reference = o.cram_reference;
    io.interleaved = o.interleaved;
    io.threads = std::max<std::size_t>(1, o.threads / 4);   // decompression: a quarter of the scorers keeps up
    io.qualities = !o.reads_out.empty() || o.error_rate < 0;
    io.screened = screened;
    io.log = o.log;
    std::size_t next_input = 0;
    std::unique_ptr<FragmentSource> source;
    ReadStats st;
    std::vector<std::pair<std::string, InputStats>> input_stats;
    auto close_source = [&] {
        if (!source) return;
        const InputStats& s = source->stats();
        st.fragments += s.unscreened + s.dropped;   // passed over before scoring, unscreened
        input_stats.emplace_back(source->describe(), s);
        log(source->describe() + ": " + std::to_string(s.records) + " records, " + std::to_string(s.fragments + s.unscreened + s.dropped) +
            " fragments (" + std::to_string(s.pairs + s.unscreened) + " pairs)" +
            (s.secondary + s.supplementary + s.qc_fail > 0 ? "; skipped " + std::to_string(s.secondary) + " secondary, " +
                 std::to_string(s.supplementary) + " supplementary, " + std::to_string(s.qc_fail) + " QC-failed" : std::string()) +
            (s.fetched > 0 ? "; " + std::to_string(s.fetched) + " mates fetched" : std::string()) +
            (s.unmatched > 0 ? "; " + std::to_string(s.unmatched) + " reads without their mate" : std::string()));
        source.reset();
    };
    constexpr std::size_t kBatch = 20000;
    // The next batch is read while the workers score the current one.
    auto read_batch = [&](std::vector<Fragment>& batch) {   // the inputs one after another
        batch.clear();
        std::vector<Fragment> part;
        while (batch.size() < kBatch) {
            if (!source) {
                if (next_input == o.inputs.size()) break;
                const auto& [a, b] = o.inputs[next_input++];
                source = b.empty() ? open_input(a, io) : open_fastq_pair(a, b);
            }
            if (source->next(part, kBatch - batch.size()) == 0) {
                close_source();
                continue;
            }
            for (auto& f : part) batch.push_back(std::move(f));
        }
    };
    std::vector<Fragment> batch, next;
    read_batch(batch);
    if (o.error_rate < 0) {
        o.error_rate = estimated_error(batch, st.error_source);
    } else {
        st.error_source = "given";
    }
    st.error_rate = o.error_rate;
    {
        char b[64];
        std::snprintf(b, sizeof b, "%.4f", o.error_rate);
        log(std::string("error rate ") + b + " (" + st.error_source + ")");
    }
    // each locus's private-site cut: calibrated for the run's read class (from the first reads), or the given one
    std::vector<std::string> cut_note(loci.size());
    {
        std::uint64_t bases = 0, reads = 0;
        for (const Fragment& f : batch) {
            bases += f.seq1.size() + f.seq2.size();
            reads += !f.seq1.empty() + !f.seq2.empty();
        }
        const double mean = reads ? static_cast<double>(bases) / static_cast<double>(reads) : 0;
        for (std::size_t li = 0; li < loci.size(); ++li) {
            char b[256];
            if (!o.calibrated_cut || o.method != Method::Runs) {
                std::snprintf(b, sizeof b, "%u\tgiven", o.min_sites);
            } else {
                const Cut c = calibrated_cut(calibration[li], mean, o.error_rate, o.false_certain_rate);
                if (c.min_sites == 0) {
                    std::snprintf(b, sizeof b, "%u\tthe index carries no calibration", o.min_sites);
                } else {
                    loci[li].min_sites = c.min_sites;
                    const std::uint64_t scr = c.by->flank_screened + c.by->decoy_screened;
                    const std::uint64_t fc = c.by->flank_certain[c.min_sites] + c.by->decoy_certain[c.min_sites];
                    std::snprintf(b, sizeof b, "%u\tcalibrated (%s reads%s): held-out recall %.4f; false certain %llu of %llu screened off-target, 95%% bound %.3g",
                                  c.min_sites, c.by->read_class.c_str(), c.met ? "" : "; no cut reaches the rate, the largest is used",
                                  c.by->locus ? static_cast<double>(c.by->locus_certain[c.min_sites]) / static_cast<double>(c.by->locus) : 0.0,
                                  static_cast<unsigned long long>(fc), static_cast<unsigned long long>(scr), upper_bound(fc, scr));
                }
            }
            cut_note[li] = b;
            log(loci[li].name + ": private-site cut " + std::string(b).replace(std::string(b).find('\t'), 1, ", "));
        }
    }

    std::FILE* out = std::fopen(o.out.c_str(), "w");
    if (out == nullptr) throw std::runtime_error("recruit: cannot write " + o.out);
    std::fprintf(out, "fragment\tlocus\tfate\treason\tmates\tseg1\tseg2\tsites_p\tsites_d\tdecoy_out\truns_n\tseg_kmers\tn_p\tn_s\tn_r\tn_d\tn_n\tkmers\twindows\tbases\tmulti_p\tnodes\tother_loci\tcertain_loci\n");
    std::unique_ptr<ReadWriters> writers;
    if (!o.reads_out.empty()) writers = std::make_unique<ReadWriters>(o.reads_out);
    const Costs costs = costs_for(o.error_rate);
    const double gap_cost = costs.gap_cost, flank_share = costs.flank_share;

    // per locus: fragments by fate and reason, backbone k-mers in certain fragments, certain in several loci
    std::vector<std::map<std::string, std::uint64_t>> by_reason(loci.size());
    std::vector<std::uint64_t> backbone_seen(loci.size(), 0), multi_certain(loci.size(), 0);
    std::uint64_t read_bases = 0, reads_seen = 0;
    const std::size_t threads = std::max<std::size_t>(1, o.threads);
    std::vector<std::vector<Row>> rows;
    while (!batch.empty()) {
        rows.assign(batch.size(), {});
        std::atomic<std::size_t> cursor{0};
        std::vector<std::thread> pool;
        for (std::size_t t = 0; t < threads; ++t)
            pool.emplace_back([&] {
                std::vector<std::uint32_t> hit;
                for (std::size_t i; (i = cursor++) < batch.size();) {
                    const Fragment& fr = batch[i];
                    hit.clear();
                    for (const std::string* s : {&fr.seq1, &fr.seq2})
                        kmer::for_each_closed_syncmer(s->data(), s->size(), kK, kS, [&](std::uint64_t c, std::size_t, bool) {
                            if (const std::uint32_t* id = screen.find(c)) hit.insert(hit.end(), screen_loci[*id].begin(), screen_loci[*id].end());
                        });
                    std::sort(hit.begin(), hit.end());
                    hit.erase(std::unique(hit.begin(), hit.end()), hit.end());
                    for (const std::uint32_t li : hit) rows[i].push_back(score_fragment(fr, loci[li], li, o, gap_cost, flank_share));
                }
            });
        std::exception_ptr failed;
        try { read_batch(next); } catch (...) { failed = std::current_exception(); }
        for (auto& t : pool) t.join();
        if (failed) std::rethrow_exception(failed);
        for (std::size_t i = 0; i < batch.size(); ++i) {
            ++st.fragments;
            const Fragment& fr = batch[i];
            read_bases += fr.seq1.size() + fr.seq2.size();
            reads_seen += 1 + !fr.seq2.empty();
            if (rows[i].empty()) continue;
            ++st.screened;
            std::uint32_t certain_loci = 0;
            for (const Row& r : rows[i]) certain_loci += r.fate[0] == 'c';
            for (const Row& r : rows[i]) {
                const Evidence& v = r.ev;
                std::string others;
                for (const Row& x : rows[i])
                    if (x.locus != r.locus) others += (others.empty() ? "" : ",") + loci[x.locus].name;
                if (others.empty()) others = ".";
                std::fprintf(out, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%s\t%s\t%u\n", r.fragment.c_str(),
                             loci[r.locus].name.c_str(), r.fate, r.reason, r.mates, r.seg[0].c_str(), r.seg[1].c_str(), v.sites_p, v.sites_d,
                             v.decoy_out, v.runs_n, v.seg_kmers, v.n_p, v.n_s, v.n_r, v.n_d, v.n_n, v.kmers, v.windows, v.bases, v.multi,
                             r.nodes.c_str(), others.c_str(), certain_loci);
                ++by_reason[r.locus][std::string(r.fate) + (r.reason[0] == '.' ? "" : std::string("\t") + r.reason)];
                if (r.fate[0] == 'c') {
                    ++st.certain;
                    backbone_seen[r.locus] += r.backbone;
                    if (certain_loci > 1) ++multi_certain[r.locus];
                } else if (r.fate[0] == 'a') {
                    ++st.ambiguous;
                } else {
                    ++st.unrecruited;
                }
                if (writers && (r.fate[0] == 'c' || (r.fate[0] == 'a' && o.reads_ambiguous))) {
                    const std::string header = fr.name + " fate=" + r.fate + " reason=" + r.reason + " mates=" + r.mates +
                                               " sites_p=" + std::to_string(v.sites_p) + " other_loci=" + others;
                    const bool two = fr.paired && !fr.seq1.empty() && !fr.seq2.empty();
                    if (two) {
                        writers->write(loci[r.locus].name, r.fate, "_1", header, fr.seq1, fr.qual1);
                        writers->write(loci[r.locus].name, r.fate, "_2", header, fr.seq2, fr.qual2);
                    } else {
                        const bool first = !fr.seq1.empty();
                        writers->write(loci[r.locus].name, r.fate, "", header, first ? fr.seq1 : fr.seq2, first ? fr.qual1 : fr.qual2);
                    }
                }
            }
        }
        if (std::ferror(out) || std::fflush(out) != 0) throw std::runtime_error("recruit: cannot write " + o.out);
        batch.swap(next);
    }
    close_source();
    if (writers) writers->close_all();
    if (std::fclose(out) != 0) throw std::runtime_error("recruit: cannot write " + o.out);

    if (!o.summary.empty()) {
        std::FILE* sm = std::fopen(o.summary.c_str(), "w");
        if (sm == nullptr) throw std::runtime_error("recruit: cannot write " + o.summary);
        std::fprintf(sm, "section\tkey\tvalue\textra\n");
        std::fprintf(sm, "run\terror_rate\t%.5f\t%s\n", o.error_rate, st.error_source.c_str());
        std::fprintf(sm, "run\tfragments\t%llu\t\nrun\tscreened\t%llu\t\n", static_cast<unsigned long long>(st.fragments),
                     static_cast<unsigned long long>(st.screened));
        const double mean_read = reads_seen ? static_cast<double>(read_bases) / static_cast<double>(reads_seen) : 0;
        std::fprintf(sm, "run\tmean_read_length\t%.1f\t\n", mean_read);
        for (const auto& [name, s] : input_stats)
            std::fprintf(sm, "input\t%s\t%llu\trecords; skipped %llu secondary, %llu supplementary, %llu QC-failed; %llu duplicates, %llu unmapped kept; %llu fragments (%llu pairs), %llu unscreened pairs; %llu mates fetched, %llu without their mate\n",
                         name.c_str(), static_cast<unsigned long long>(s.records), static_cast<unsigned long long>(s.secondary),
                         static_cast<unsigned long long>(s.supplementary), static_cast<unsigned long long>(s.qc_fail),
                         static_cast<unsigned long long>(s.duplicates), static_cast<unsigned long long>(s.unmapped),
                         static_cast<unsigned long long>(s.fragments + s.unscreened + s.dropped), static_cast<unsigned long long>(s.pairs + s.unscreened),
                         static_cast<unsigned long long>(s.unscreened), static_cast<unsigned long long>(s.fetched), static_cast<unsigned long long>(s.unmatched));
        for (std::size_t li = 0; li < loci.size(); ++li) {
            for (const auto& [fate_reason, n] : by_reason[li]) {
                const auto tab = fate_reason.find('\t');
                const std::string fate = fate_reason.substr(0, tab), reason = tab == std::string::npos ? "." : fate_reason.substr(tab + 1);
                std::fprintf(sm, "locus\t%s\t%s\t%s %llu\n", loci[li].name.c_str(), fate.c_str(), reason.c_str(), static_cast<unsigned long long>(n));
            }
            // depth: backbone k-mers seen in certain fragments per backbone k-mer of the index, in k-mers
            // and, for reads of the mean length, in bases
            const double kdepth = loci[li].backbone ? static_cast<double>(backbone_seen[li]) / static_cast<double>(loci[li].backbone) : 0;
            const double bdepth = mean_read > kK ? kdepth * mean_read / (mean_read - kK + 1) : kdepth;
            std::fprintf(sm, "locus\t%s\tbackbone_depth\t%.2f k-mer depth, %.2f base depth over %llu backbone k-mers\n", loci[li].name.c_str(), kdepth, bdepth,
                         static_cast<unsigned long long>(loci[li].backbone));
            std::fprintf(sm, "locus\t%s\tcertain_in_other_loci\t%llu\n", loci[li].name.c_str(), static_cast<unsigned long long>(multi_certain[li]));
            std::fprintf(sm, "locus\t%s\tsite_cut\t%s\n", loci[li].name.c_str(), cut_note[li].c_str());
        }
        if (std::fclose(sm) != 0) throw std::runtime_error("recruit: cannot write " + o.summary);
    }
    return st;
}

}  // namespace panvar::recruit
