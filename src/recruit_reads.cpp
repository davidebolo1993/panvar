#include "panvar/recruit_reads.hpp"

#include "panvar/cli_utils.hpp"

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
#include <mutex>
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
// for single reads. Records are collected per file; each file's text is compressed in pieces of kPiece
// bytes (the last one shorter), each piece one gzip member, by several threads at each flush; a gzip file of
// several members reads as one. Pieces are cut by size alone, so the files depend neither on the thread count
// nor on the batches. Files open when first written; past a cap they are closed and reopened to append.
class ReadWriters {
public:
    ReadWriters(std::string dir, std::size_t threads) : dir_(std::move(dir)), threads_(std::max<std::size_t>(1, threads)) {
        std::filesystem::create_directories(dir_);
    }
    ~ReadWriters() {
        try { close_all(); } catch (...) {}
    }
    void write(const std::string& locus, const char* fate, const char* suffix, const std::string& header, const std::string& seq,
               const std::string& qual) {
        std::string& to = pending_[locus + "." + fate + suffix + ".fq.gz"];
        to += '@';
        to += header;
        to += '\n';
        to += seq;
        to += "\n+\n";
        if (qual.size() == seq.size()) to += qual;
        else to.append(seq.size(), 'I');
        to += '\n';
        pending_bytes_ += header.size() + 2 * seq.size() + 6;
        // bounds memory; the pieces stay the same. Text short of a whole piece waits (under kPiece per file), so
        // when a flush finds none whole, the next waits until kFlushAt more has come rather than scanning again
        if (pending_bytes_ >= flush_at_) {
            flush();
            flush_at_ = pending_bytes_ + kFlushAt;
        }
    }
    // Compresses each file's whole pieces written so far (with `all`, the rest too) and appends them.
    void flush(bool all = false) {
        struct Piece {
            const std::string* key;
            const char* data;
            std::size_t size;
            std::string gz;
        };
        std::vector<Piece> pieces;
        for (const auto& [key, text] : pending_)
            for (std::size_t a = 0; a < text.size() && (all || text.size() - a >= kPiece); a += kPiece)
                pieces.push_back({&key, text.data() + a, std::min(kPiece, text.size() - a), {}});
        std::atomic<std::size_t> cursor{0};
        std::exception_ptr failed;
        std::mutex failed_mutex;
        auto work = [&] {
            try {
                for (std::size_t i; (i = cursor++) < pieces.size();) pieces[i].gz = gzip_member(pieces[i].data, pieces[i].size);
            } catch (...) {
                std::lock_guard<std::mutex> lock(failed_mutex);
                if (!failed) failed = std::current_exception();
            }
        };
        std::vector<std::thread> pool;
        for (std::size_t t = 1; t < std::min(threads_, pieces.size()); ++t) pool.emplace_back(work);
        work();
        for (auto& t : pool) t.join();
        if (failed) std::rethrow_exception(failed);
        std::map<std::string, std::size_t> done;   // bytes compressed per file
        for (const Piece& p : pieces) {
            std::FILE* f = file(*p.key);
            if (std::fwrite(p.gz.data(), 1, p.gz.size(), f) != p.gz.size()) throw std::runtime_error("recruit: cannot write " + dir_ + "/" + *p.key);
            done[*p.key] += p.size;
        }
        for (const auto& [key, n] : done) {
            std::string& text = pending_[key];
            if (n == text.size()) {
                pending_.erase(key);
            } else {
                text.erase(0, n);
                if (text.capacity() > 2 * kPiece) text.shrink_to_fit();
            }
        }
        pending_bytes_ = 0;
        for (const auto& [key, text] : pending_) pending_bytes_ += text.size();
    }
    void close_all() {
        flush(true);
        close_files();
    }

private:
    static constexpr std::size_t kMaxOpen = 200;
    static constexpr std::size_t kPiece = std::size_t{4} << 20;
    static constexpr std::size_t kFlushAt = std::size_t{256} << 20;   // text waiting before a flush
    std::size_t flush_at_ = kFlushAt;

    static std::string gzip_member(const char* data, std::size_t size) {
        z_stream zs{};
        if (deflateInit2(&zs, 3, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY) != Z_OK) throw std::runtime_error("recruit: cannot start compression");
        std::string out(deflateBound(&zs, static_cast<uLong>(size)), '\0');
        zs.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(data));
        zs.avail_in = static_cast<uInt>(size);
        zs.next_out = reinterpret_cast<Bytef*>(out.data());
        zs.avail_out = static_cast<uInt>(out.size());
        const int ret = deflate(&zs, Z_FINISH);
        out.resize(out.size() - zs.avail_out);
        deflateEnd(&zs);
        if (ret != Z_STREAM_END) throw std::runtime_error("recruit: cannot compress the recruited reads");
        return out;
    }
    std::FILE* file(const std::string& key) {
        auto it = open_.find(key);
        if (it != open_.end()) return it->second;
        if (open_.size() >= kMaxOpen) close_files();
        const bool again = started_.count(key) > 0;
        std::FILE* f = std::fopen((std::filesystem::path(dir_) / key).string().c_str(), again ? "ab" : "wb");
        if (f == nullptr) throw std::runtime_error("recruit: cannot write " + dir_ + "/" + key);
        started_.insert(key);
        return open_.emplace(key, f).first->second;
    }
    void close_files() {
        std::string bad;
        for (auto& [key, f] : open_)
            if (std::fclose(f) != 0 && bad.empty()) bad = key;
        open_.clear();
        if (!bad.empty()) throw std::runtime_error("recruit: cannot write " + dir_ + "/" + bad);
    }

    std::string dir_;
    std::size_t threads_;
    std::map<std::string, std::string> pending_;
    std::size_t pending_bytes_ = 0;
    std::map<std::string, std::FILE*> open_;
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

// The read's seeds for the screen: closed syncmers, minimizers (the smallest hash of each window of kWindow
// consecutive k-mers, each reported once), or every canonical k-mer.
constexpr std::size_t kWindow = 20;
template <class F>
void for_each_seed(Screen mode, const std::string& s, F&& f) {
    if (mode == Screen::Syncmer) {
        kmer::for_each_closed_syncmer(s.data(), s.size(), kK, kS, [&](std::uint64_t c, std::size_t, bool) { f(c); });
        return;
    }
    if (mode == Screen::All) {
        kmer::for_each_canonical_kmer(s.data(), s.size(), kK, [&](std::uint64_t c, std::size_t, bool) { f(c); });
        return;
    }
    struct Seed { std::uint64_t hash, code; std::size_t at; };
    std::vector<Seed> w;   // the window's candidates, increasing hash (a monotone queue)
    std::size_t last = ~std::size_t{0}, prev = ~std::size_t{0}, emitted = ~std::size_t{0};
    kmer::for_each_canonical_kmer(s.data(), s.size(), kK, [&](std::uint64_t c, std::size_t at, bool) {
        if (prev != ~std::size_t{0} && at != prev + 1) w.clear();   // an N broke the run of k-mers
        if (w.empty()) last = at;
        prev = at;
        const std::uint64_t h = kmer::mix64(c);
        while (!w.empty() && w.back().hash >= h) w.pop_back();
        w.push_back({h, c, at});
        while (w.front().at + kWindow <= at) w.erase(w.begin());
        if (at + 1 >= last + kWindow && w.front().at != emitted) {   // a full window: its minimizer
            emitted = w.front().at;
            f(w.front().code);
        }
    });
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
        if (o.screen == Screen::Syncmer) {
            for (const std::uint64_t s : x.screen) screen_items.emplace_back(s, static_cast<std::uint32_t>(li));
        } else {
            for (std::size_t i = 0; i < x.codes.size(); ++i)
                if (x.cls[i] != kRepeat) screen_items.emplace_back(x.codes[i], static_cast<std::uint32_t>(li));
        }
        log("locus: " + x.name + "; k-mers: " + std::to_string(x.codes.size()) + "; decoys: " + std::to_string(x.decoy.size()) +
            "; screen syncmers: " + std::to_string(x.screen.size()));
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
        for_each_seed(o.screen, seq, [&](std::uint64_t c) {
            if (!hit && screen.find(c) != nullptr) hit = true;
        });
        return hit;
    };
    if (o.inputs.empty()) throw std::runtime_error("recruit: no reads given");
    InputOptions io;
    io.cram_reference = o.cram_reference;
    io.interleaved = o.interleaved;
    io.threads = std::max<std::size_t>(1, o.threads / 4);   // decompression: a quarter of the scorers keeps up
    io.fetch_threads = std::max<std::size_t>(1, o.threads);  // the mates fetched at the end: nothing else runs then
    io.parse_threads = std::max<std::size_t>(1, o.threads / 4);   // alignment records prepared beside the pairing
    io.qualities = true;   // a read's own error rate comes from its qualities
    io.screened = screened;
    io.log = o.log;
    std::size_t next_input = 0;
    std::unique_ptr<FragmentSource> source;
    std::unique_ptr<cli::ProgressBar> input_bar;   // how far through the current input
    ReadStats st;
    std::vector<std::pair<std::string, InputStats>> input_stats;
    auto close_source = [&] {
        if (!source) return;
        if (input_bar) input_bar->done();
        input_bar.reset();
        const InputStats& s = source->stats();
        st.fragments += s.unscreened + s.dropped;   // passed over before scoring, unscreened
        input_stats.emplace_back(source->describe(), s);
        if (s.fragments + s.unscreened + s.dropped == 0) log("warning: input: " + source->describe() + "; it holds no reads");
        log("input: " + source->describe() + "; records: " + std::to_string(s.records) + "; fragments: " +
            std::to_string(s.fragments + s.unscreened + s.dropped) + "; pairs: " + std::to_string(s.pairs + s.unscreened) +
            (s.secondary + s.supplementary + s.qc_fail > 0 ? "; skipped: secondary " + std::to_string(s.secondary) + ", supplementary " +
                 std::to_string(s.supplementary) + ", QC-failed " + std::to_string(s.qc_fail) : std::string()) +
            (s.fetched > 0 ? "; mates fetched: " + std::to_string(s.fetched) : std::string()) +
            (s.unmatched > 0 ? "; reads without their mate: " + std::to_string(s.unmatched) : std::string()));
        source.reset();
    };
    constexpr std::size_t kBatch = 20000;                        // fragments in a batch at most
    constexpr std::size_t kBatchBytes = std::size_t{1} << 30;    // and bases plus qualities: long reads fill this first
    constexpr std::size_t kBatchStep = 1000;                     // fragments asked for at a time, to stop near the bound
    // One sample's inputs must hold one read type: the error rate, read class and cut come from the first reads.
    // Each input's first reads are checked against the first input's: short (under 1 kb) or long.
    int first_long = -1;               // the first input's: 0 short, 1 long, -1 not yet seen
    std::string first_input;
    bool checked_input = false;        // the current input's first reads were checked
    auto check_type = [&](const std::vector<Fragment>& part) {
        std::uint64_t bases = 0, reads = 0;
        for (const Fragment& f : part) {
            bases += f.seq1.size() + f.seq2.size();
            reads += !f.seq1.empty() + !f.seq2.empty();
        }
        if (reads == 0) return;
        checked_input = true;
        const int is_long = static_cast<double>(bases) / static_cast<double>(reads) >= 1000 ? 1 : 0;
        if (first_long < 0) { first_long = is_long; first_input = source->describe(); return; }
        if (is_long != first_long)
            throw std::runtime_error("recruit: " + source->describe() + " holds " + (is_long ? "long" : "short") + " reads and " + first_input + " " +
                                     (is_long ? "short" : "long") + " ones: one sample's inputs must be one read type (give them as separate samples)");
    };
    // The next batch is read while the workers score the current one.
    auto read_batch = [&](std::vector<Fragment>& batch) {   // the inputs one after another
        batch.clear();
        std::vector<Fragment> part;
        std::size_t bytes = 0;
        while (batch.size() < kBatch && bytes < kBatchBytes) {
            if (!source) {
                if (next_input == o.inputs.size()) break;
                const auto& [a, b] = o.inputs[next_input++];
                source = b.empty() ? open_input(a, io) : open_fastq_pair(a, b);
                checked_input = false;
                if (o.log && source->progress_total() > 0)
                    input_bar = std::make_unique<cli::ProgressBar>(
                        cli::log_prefix("recruit") + "reading input" + (o.inputs.size() > 1 ? " " + std::to_string(next_input) + " of " + std::to_string(o.inputs.size()) : std::string()),
                        source->progress_total(), source->progress_in_bases() ? cli::ProgressBar::Unit::bases : cli::ProgressBar::Unit::bytes);
            }
            if (source->next(part, std::min(kBatchStep, kBatch - batch.size())) == 0) {
                close_source();
                continue;
            }
            if (input_bar) input_bar->set(source->progress_done());
            if (!checked_input) check_type(part);
            for (auto& f : part) {
                bytes += f.seq1.size() + f.qual1.size() + f.seq2.size() + f.qual2.size();
                batch.push_back(std::move(f));
            }
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
        log(std::string("error rate: ") + b + " (" + st.error_source + ")");
    }
    // each locus's decision for the run's read class (the preset's, or from the first reads): the score's tiers,
    // or by the first rules a private-site cut, calibrated or given. A locus without the score's calibration is
    // decided by the cut.
    std::vector<std::vector<std::string>> notes(loci.size());   // summary rows: key, value
    {
        std::uint64_t bases = 0, reads = 0;
        bool pairs = false;
        for (const Fragment& f : batch) {
            bases += f.seq1.size() + f.seq2.size();
            reads += !f.seq1.empty() + !f.seq2.empty();
            pairs = pairs || (f.paired && !f.seq1.empty() && !f.seq2.empty());
        }
        const double mean = reads ? static_cast<double>(bases) / static_cast<double>(reads) : 0;
        char b[768];
        // what a run unlike the calibration should know
        auto caveats = [&](const Calibration& c) {
            std::string note;
            if (o.error_rate > 2 * c.error_rate) {
                std::snprintf(b, sizeof b, "; the run's error rate %.4f is above the calibration's %.4f: its off-target reads may be certain more often",
                              o.error_rate, c.error_rate);
                note += b;
            }
            if (c.insert > 0 && !pairs) note += "; calibrated on pairs, the run reads single reads";
            return note;
        };
        for (std::size_t li = 0; li < loci.size(); ++li) {
            const Calibration* cal = calibration_for(calibration[li], o.read_class, mean, o.error_rate);
            const bool score = cal != nullptr && has_score(*cal);
            if (score) {
                loci[li].rates = Rates{cal->r_in, cal->r_out, cal->q_in, cal->q_out};
                const Tier ct = score_tier(*cal, o.false_certain_rate), lt = score_tier(*cal, o.likely_rate);
                loci[li].certain_bin = ct.bin;
                loci[li].likely_bin = std::min(lt.bin, ct.bin);
                loci[li].tiers = true;
                loci[li].by_score = !o.rules_v1 && o.calibrated_cut;
                if (loci[li].by_score) {
                    std::snprintf(b, sizeof b, "%s reads: r_in %.4g, r_out %.4g, q_in %.4g, q_out %.4g; off-target screened %llu, of them %llu from %zu held-out genomes",
                                  cal->read_class.c_str(), cal->r_in, cal->r_out, cal->q_in, cal->q_out,
                                  static_cast<unsigned long long>(cal->flank_screened + cal->decoy_screened + cal->genome_screened),
                                  static_cast<unsigned long long>(cal->genome_screened), cal->genomes.size());
                    notes[li].push_back("score");
                    notes[li].push_back(b + caveats(*cal));
                    for (const auto& [name, t] : {std::pair<const char*, const Tier*>{"certain", &ct}, {"likely", &lt}}) {
                        char th[32];
                        if (t->bin == kEveryBin) std::snprintf(th, sizeof th, "any");
                        else std::snprintf(th, sizeof th, "%.4f", t->bin * kScoreBin);
                        std::snprintf(b, sizeof b, "%s: held-out recall %.4f; off-target %llu of %llu screened (rate %.3g asked), 95%% bound %.3g, %s", th,
                                      t->locus_total ? static_cast<double>(t->locus) / static_cast<double>(t->locus_total) : 0.0,
                                      static_cast<unsigned long long>(t->off), static_cast<unsigned long long>(t->off_screened), t->rate, t->bound(),
                                      t->shown ? "within the rate" : t->off == 0 ? "too few to show the rate (none beyond it seen)"
                                                                     : "within the rate as seen, its bound above it");
                        notes[li].push_back(std::string(name) + "_threshold");
                        notes[li].push_back(b);
                    }
                    log("locus: " + loci[li].name + "; " + notes[li][1] + "; certain from " + notes[li][3] + "; likely from " + notes[li][5]);
                    continue;
                }
            }
            std::string note;
            if (!o.calibrated_cut) {
                note = std::to_string(o.min_sites) + " given";
            } else {
                const Cut c = calibrated_cut(calibration[li], o.read_class, mean, o.error_rate, o.false_certain_rate);
                const std::uint32_t fallback = std::max(o.min_sites, c.floor);
                if (c.by == nullptr) {
                    loci[li].min_sites = fallback;
                    note = std::to_string(fallback) + " (the index carries no calibration)";
                } else if (c.min_sites == 0) {
                    loci[li].min_sites = fallback;
                    note = std::to_string(fallback) + " (the calibration screened no off-target read)";
                } else {
                    loci[li].min_sites = c.min_sites;
                    const std::uint64_t scr = c.by->flank_screened + c.by->decoy_screened;
                    const std::uint64_t fc = c.by->flank_certain[c.min_sites] + c.by->decoy_certain[c.min_sites];
                    std::snprintf(b, sizeof b, "%u calibrated (%s reads): held-out recall %.4f; false certain %llu of %llu screened off-target, 95%% bound %.3g, %s",
                                  c.min_sites, c.by->read_class.c_str(),
                                  c.by->locus ? static_cast<double>(c.by->locus_certain[c.min_sites]) / static_cast<double>(c.by->locus) : 0.0,
                                  static_cast<unsigned long long>(fc), static_cast<unsigned long long>(scr), upper_bound(fc, scr),
                                  !c.met ? "no cut reaches the rate: the largest is used"
                                  : c.shown ? "within the rate" : "too few to show the rate (none beyond it seen)");
                    note = b;
                    if (c.floor > 2) note += "; " + c.by->read_class + " reads start at " + std::to_string(c.floor);
                    note += caveats(*c.by);
                    {   // the calibration scored by the default rules: a run scoring by others is not what it measured
                        const ReadOptions d;
                        if (o.decoy_share != d.decoy_share || o.decoy_share_long != d.decoy_share_long) {
                            std::snprintf(b, sizeof b, "; the run's decoy shares (%.3g, %.3g long) are not the calibration's (%.3g, %.3g)",
                                          o.decoy_share, o.decoy_share_long, d.decoy_share, d.decoy_share_long);
                            note += b;
                        }
                    }
                }
            }
            if (!o.rules_v1 && o.calibrated_cut) note += "; the index carries no score calibration for these reads: decided by the private-site cut";
            notes[li] = {"site_cut", note};
            log("locus: " + loci[li].name + "; site cut: " + note);
        }
    }

    // the table and summary are written under temporary names and renamed when complete; an earlier run's are
    // removed first, so a run that fails leaves none that looks finished
    const std::string out_tmp = o.out + ".tmp", summary_tmp = o.summary.empty() ? std::string() : o.summary + ".tmp";
    {
        std::error_code ec;
        std::filesystem::remove(o.out, ec);
        if (!o.summary.empty()) std::filesystem::remove(o.summary, ec);
    }
    std::FILE* out = std::fopen(out_tmp.c_str(), "w");
    if (out == nullptr) throw std::runtime_error("recruit: cannot write " + out_tmp);
    std::fprintf(out, "fragment\tlocus\tfate\treason\tmates\tseg1\tseg2\tsites_p\tsites_d\tdecoy_out\truns_n\tseg_kmers\tn_p\tn_s\tn_r\tn_d\tn_n\tkmers\twindows\tbases\tmulti_p\tnodes\tother_loci\tcertain_loci\tscore\ttier\n");
    std::unique_ptr<ReadWriters> writers;
    if (!o.reads_out.empty()) {
        // an earlier run's reads of these loci go first: a locus and fate with no reads now must not keep old ones
        std::error_code ec;
        std::size_t removed = 0;
        // only the names this run writes: another locus's name may start with this one's
        std::set<std::string> ours;
        for (const LocusTable& l : loci)
            for (const char* fate : {"certain", "likely", "ambiguous"})
                for (const char* suffix : {"", "_1", "_2"}) ours.insert(l.name + "." + fate + suffix + ".fq.gz");
        if (std::filesystem::is_directory(o.reads_out, ec))
            for (const auto& e : std::filesystem::directory_iterator(o.reads_out, ec))
                if (ours.count(e.path().filename().string())) { std::filesystem::remove(e.path(), ec); ++removed; }
        if (removed > 0) log("reads-out: " + o.reads_out + "; " + std::to_string(removed) + " files of an earlier run removed");
        writers = std::make_unique<ReadWriters>(o.reads_out, o.threads);
    }
    const Costs costs = costs_for(o.error_rate);
    const double gap_cost = costs.gap_cost, flank_share = costs.flank_share;

    // per locus: fragments by fate and reason, backbone k-mers in certain fragments, certain in several loci
    std::vector<std::map<std::string, std::uint64_t>> by_reason(loci.size());
    std::vector<std::uint64_t> backbone_seen(loci.size(), 0), multi_certain(loci.size(), 0);
    std::uint64_t read_bases = 0, reads_seen = 0;
    const std::size_t threads = std::max<std::size_t>(1, o.threads);
    std::vector<std::vector<Row>> rows;
    constexpr std::uint64_t kReportEvery = 10000000;
    std::uint64_t next_report = kReportEvery;
    while (!batch.empty()) {
        rows.assign(batch.size(), {});
        std::atomic<std::size_t> cursor{0};
        std::vector<std::thread> pool;
        std::mutex score_mu;
        std::exception_ptr score_failed;   // a scorer's error, rethrown once all have stopped
        for (std::size_t t = 0; t < threads; ++t)
            pool.emplace_back([&] {
              try {
                std::vector<std::uint32_t> hit;
                for (std::size_t i; (i = cursor++) < batch.size();) {
                    const Fragment& fr = batch[i];
                    hit.clear();
                    for (const std::string* s : {&fr.seq1, &fr.seq2})
                        for_each_seed(o.screen, *s, [&](std::uint64_t c) {
                            if (const std::uint32_t* id = screen.find(c)) hit.insert(hit.end(), screen_loci[*id].begin(), screen_loci[*id].end());
                        });
                    std::sort(hit.begin(), hit.end());
                    hit.erase(std::unique(hit.begin(), hit.end()), hit.end());
                    for (const std::uint32_t li : hit) rows[i].push_back(score_fragment(fr, loci[li], li, o, gap_cost, flank_share));
                }
              } catch (...) {
                std::lock_guard<std::mutex> lock(score_mu);
                if (!score_failed) score_failed = std::current_exception();
                cursor = batch.size();
              }
            });
        std::exception_ptr failed;
        try { read_batch(next); } catch (...) { failed = std::current_exception(); }
        for (auto& t : pool) t.join();
        if (score_failed) std::rethrow_exception(score_failed);
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
                char score[32] = ".";
                if (r.scored) std::snprintf(score, sizeof score, "%.2f", r.score);
                std::fprintf(out, "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%s\t%s\t%u\t%s\t%s\n", r.fragment.c_str(),
                             loci[r.locus].name.c_str(), r.fate, r.reason, r.mates, r.seg[0].c_str(), r.seg[1].c_str(), v.sites_p, v.sites_d,
                             v.decoy_out, v.runs_n, v.seg_kmers, v.n_p, v.n_s, v.n_r, v.n_d, v.n_n, v.kmers, v.windows, v.bases, v.multi,
                             r.nodes.c_str(), others.c_str(), certain_loci, score, r.tier);
                ++by_reason[r.locus][std::string(r.fate) + (r.reason[0] == '.' ? "" : std::string("\t") + r.reason)];
                if (r.fate[0] == 'c') {
                    ++st.certain;
                    backbone_seen[r.locus] += r.backbone;
                    if (certain_loci > 1) ++multi_certain[r.locus];
                } else if (r.fate[0] == 'l') {
                    ++st.likely;
                } else if (r.fate[0] == 'a') {
                    ++st.ambiguous;
                } else {
                    ++st.unrecruited;
                }
                if (writers && (r.fate[0] == 'c' || ((r.fate[0] == 'a' || r.fate[0] == 'l') && o.reads_ambiguous))) {
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
        if (writers) writers->flush();
        if (std::ferror(out) || std::fflush(out) != 0) throw std::runtime_error("recruit: cannot write " + o.out);
        if (st.fragments >= next_report) {   // a line every 10 million fragments: a stream has no total for a bar
            log("progress: fragments: " + std::to_string(st.fragments) + "; screened: " + std::to_string(st.screened) + "; certain: " +
                std::to_string(st.certain));
            while (next_report <= st.fragments) next_report += kReportEvery;
        }
        batch.swap(next);
    }
    close_source();
    if (writers) writers->close_all();
    if (std::fclose(out) != 0) throw std::runtime_error("recruit: cannot write " + out_tmp);
    std::filesystem::rename(out_tmp, o.out);

    if (!o.summary.empty()) {
        std::FILE* sm = std::fopen(summary_tmp.c_str(), "w");
        if (sm == nullptr) throw std::runtime_error("recruit: cannot write " + summary_tmp);
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
            for (std::size_t k = 0; k + 1 < notes[li].size(); k += 2)
                std::fprintf(sm, "locus\t%s\t%s\t%s\n", loci[li].name.c_str(), notes[li][k].c_str(), notes[li][k + 1].c_str());
        }
        if (std::fclose(sm) != 0) throw std::runtime_error("recruit: cannot write " + summary_tmp);
        std::filesystem::rename(summary_tmp, o.summary);
    }
    return st;
}

}  // namespace panvar::recruit
