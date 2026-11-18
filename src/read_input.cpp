#include "panvar/read_input.hpp"

#include "panvar/gz_reader.hpp"

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>
#include <pthread.h>

#include <algorithm>
#include <array>
#include <condition_variable>
#include <cstdlib>
#include <deque>
#include <exception>
#include <filesystem>
#include <map>
#include <mutex>
#include <functional>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <unordered_map>

namespace panvar::recruit {

namespace {

// A thread with an 8 MB stack, Linux's default: htslib's CRAM decoders keep large tables on the stack, more than
// the 512 KB macOS gives a std::thread.
class BigStackThread {
public:
    explicit BigStackThread(std::function<void()> f) : f_(std::move(f)) {
        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_attr_setstacksize(&attr, 8u << 20);
        const int r = pthread_create(&id_, &attr, &BigStackThread::run, this);
        pthread_attr_destroy(&attr);
        if (r != 0) throw std::runtime_error("recruit: cannot start a thread");
    }
    ~BigStackThread() { if (!joined_) pthread_join(id_, nullptr); }
    BigStackThread(const BigStackThread&) = delete;
    BigStackThread& operator=(const BigStackThread&) = delete;
    void join() { pthread_join(id_, nullptr); joined_ = true; }

private:
    static void* run(void* self) {
        static_cast<BigStackThread*>(self)->f_();
        return nullptr;
    }
    std::function<void()> f_;
    pthread_t id_{};
    bool joined_ = false;
};

bool read_fastq(GzLineReader& in, std::string& name, std::string& seq, std::string& qual) {
    std::string h, plus;
    while (in.getline(h)) {
        if (h.empty()) continue;
        if (h[0] != '@') throw std::runtime_error("recruit: not FASTQ (a record starts with '" + h.substr(0, 1) + "')");
        if (!in.getline(seq) || !in.getline(plus) || !in.getline(qual)) throw std::runtime_error("recruit: a truncated FASTQ record");
        if (!seq.empty() && seq.back() == '\r') seq.pop_back();
        if (!qual.empty() && qual.back() == '\r') qual.pop_back();
        if (plus.empty() || plus[0] != '+' || qual.size() != seq.size())
            throw std::runtime_error("recruit: a malformed FASTQ record at " + h.substr(0, 60) + " (multi-line FASTQ is not read)");
        const auto sp = h.find_first_of(" \t");
        name = h.substr(1, sp == std::string::npos ? std::string::npos : sp - 1);
        if (name.size() > 2 && name[name.size() - 2] == '/' && (name.back() == '1' || name.back() == '2')) name.resize(name.size() - 2);
        return true;
    }
    return false;
}

std::unique_ptr<GzLineReader> open_text(const std::string& path) {
    auto in = std::make_unique<GzLineReader>(path);
    if (!in->ok()) throw std::runtime_error("recruit: cannot read " + path);
    return in;
}

class FastqPair final : public FragmentSource {
public:
    FastqPair(std::string r1, std::string r2) : r1_(std::move(r1)), r2_(std::move(r2)), in1_(open_text(r1_)), in2_(open_text(r2_)) {}
    std::size_t next(std::vector<Fragment>& out, std::size_t n) override {
        out.resize(n);
        if (names2_.size() < n) names2_.resize(n);
        std::size_t n1 = 0, n2 = 0;
        std::exception_ptr failed1, failed2;
        std::thread mates([&] {
            try {
                while (n2 < n && read_fastq(*in2_, names2_[n2], out[n2].seq2, out[n2].qual2)) ++n2;
            } catch (...) {
                failed2 = std::current_exception();
            }
        });
        try {
            while (n1 < n && read_fastq(*in1_, out[n1].name, out[n1].seq1, out[n1].qual1)) out[n1++].paired = true;
        } catch (...) {
            failed1 = std::current_exception();
        }
        mates.join();
        if (failed1) std::rethrow_exception(failed1);
        if (failed2) std::rethrow_exception(failed2);
        if (n2 < n1) throw std::runtime_error("recruit: " + r2_ + " has fewer reads than " + r1_);
        if (n2 > n1) throw std::runtime_error("recruit: " + r2_ + " has more reads than " + r1_);
        if (n1 < n) {   // the first file ended: so must the second
            std::string name, seq, qual;
            if (read_fastq(*in2_, name, seq, qual)) throw std::runtime_error("recruit: " + r2_ + " has more reads than " + r1_);
        }
        for (std::size_t i = 0; i < n1; ++i)
            if (names2_[i] != out[i].name) throw std::runtime_error("recruit: mates out of step: " + out[i].name + " and " + names2_[i]);
        out.resize(n1);
        st_.fragments += n1;
        st_.pairs += n1;
        st_.records += 2 * n1;
        return n1;
    }
    const InputStats& stats() const override { return st_; }
    std::string describe() const override { return r1_ + " + " + r2_ + " (FASTQ pairs)"; }

private:
    std::string r1_, r2_;
    std::unique_ptr<GzLineReader> in1_, in2_;
    std::vector<std::string> names2_;
    InputStats st_;
};

class FastqSingle final : public FragmentSource {
public:
    FastqSingle(std::string path, bool interleaved) : path_(std::move(path)), interleaved_(interleaved), in_(open_text(path_)) {}
    std::size_t next(std::vector<Fragment>& out, std::size_t n) override {
        out.resize(n);
        std::size_t k = 0;
        std::string name2;
        while (k < n && read_fastq(*in_, out[k].name, out[k].seq1, out[k].qual1)) {
            ++st_.records;
            out[k].paired = interleaved_;
            out[k].seq2.clear();
            out[k].qual2.clear();
            if (interleaved_) {
                if (!read_fastq(*in_, name2, out[k].seq2, out[k].qual2)) throw std::runtime_error("recruit: " + path_ + ": an odd number of interleaved records");
                if (name2 != out[k].name) throw std::runtime_error("recruit: " + path_ + ": interleaved mates out of step: " + out[k].name + " and " + name2);
                ++st_.records;
                ++st_.pairs;
            }
            ++k;
        }
        out.resize(k);
        st_.fragments += k;
        return k;
    }
    const InputStats& stats() const override { return st_; }
    std::string describe() const override { return path_ + (interleaved_ ? " (interleaved FASTQ)" : " (FASTQ)"); }

private:
    std::string path_;
    bool interleaved_;
    std::unique_ptr<GzLineReader> in_;
    InputStats st_;
};

// One open alignment file with its header, a record buffer and, for CRAM, the reference and the fields read.
struct HtsFile {
    samFile* fp = nullptr;
    sam_hdr_t* hdr = nullptr;
    bam1_t* rec = nullptr;
    std::string path;
    HtsFile(const std::string& p, const InputOptions& o, htsThreadPool* pool) : path(p) {
        fp = sam_open(p.c_str(), "r");
        if (fp == nullptr) throw std::runtime_error("recruit: cannot open " + p);
        const htsFormat* f = hts_get_format(fp);
        if (f->format != sam && f->format != bam && f->format != cram) {
            sam_close(fp);
            throw std::runtime_error("recruit: " + p + " is not SAM, BAM or CRAM");
        }
        auto fail = [&](const std::string& why) {
            if (hdr) sam_hdr_destroy(hdr);
            sam_close(fp);
            throw std::runtime_error("recruit: " + p + ": " + why);
        };
        // a BGZF or CRAM file whose writer stopped between blocks lacks the end marker: truncated
        if (hts_check_EOF(fp) == 0) fail("truncated (no end-of-file marker)");
        if (pool != nullptr && pool->pool != nullptr) hts_set_thread_pool(fp, pool);
        if (f->format == cram) {
            if (!o.cram_reference.empty() && hts_set_fai_filename(fp, o.cram_reference.c_str()) != 0)
                fail("cannot use the CRAM reference " + o.cram_reference);
            hts_set_opt(fp, CRAM_OPT_REQUIRED_FIELDS,
                        SAM_QNAME | SAM_FLAG | SAM_RNAME | SAM_POS | SAM_RNEXT | SAM_PNEXT | SAM_SEQ | (o.qualities ? SAM_QUAL : 0));
        }
        hdr = sam_hdr_read(fp);
        if (hdr == nullptr) fail("cannot read the header");
        rec = bam_init1();
    }
    ~HtsFile() {
        if (rec) bam_destroy1(rec);
        if (hdr) sam_hdr_destroy(hdr);
        if (fp) sam_close(fp);
    }
    HtsFile(const HtsFile&) = delete;
    HtsFile& operator=(const HtsFile&) = delete;
    // The next record; false at the end; a damaged file (truncated BGZF, a CRAM whose reference cannot be
    // found) stops the run with the file named.
    bool read() { return read(rec); }
    bool read(bam1_t* into) {
        const int r = sam_read1(fp, hdr, into);
        if (r >= 0) return true;
        if (r == -1) return false;
        throw std::runtime_error("recruit: " + path + ": cannot read a record (truncated or damaged, or a CRAM whose reference cannot be found; see --cram-reference)");
    }
    // Sorted by coordinate: the header says so, or it says nothing and the file has an index.
    bool coordinate_sorted() const {
        kstring_t ks = KS_INITIALIZE;
        const bool tagged = sam_hdr_find_tag_hd(hdr, "SO", &ks) == 0 && ks.s != nullptr;
        const std::string so = tagged ? ks.s : "";
        ks_free(&ks);
        if (so == "coordinate") return true;
        if (tagged && so != "unknown") return false;
        hts_idx_t* idx = sam_index_load(fp, path.c_str());
        if (idx == nullptr) return false;
        hts_idx_destroy(idx);
        return true;
    }
    std::string sample() const {
        std::string sample;
        const int n = sam_hdr_count_lines(hdr, "RG");
        for (int i = 0; i < n; ++i) {
            kstring_t ks = KS_INITIALIZE;
            if (sam_hdr_find_tag_pos(hdr, "RG", i, "SM", &ks) == 0 && ks.s != nullptr) {
                if (sample.empty()) sample = ks.s;
                else if (sample != ks.s) sample = "several";
            }
            ks_free(&ks);
        }
        return sample;
    }
};

// Two bases per byte of BAM's 4-bit sequence, and the complement of a base (anything but ACGT is N).
const std::array<std::array<char, 2>, 256> kPairs = [] {
    std::array<std::array<char, 2>, 256> t{};
    for (int i = 0; i < 256; ++i) t[i] = {seq_nt16_str[i >> 4], seq_nt16_str[i & 15]};
    return t;
}();
const std::array<char, 256> kComplement = [] {
    std::array<char, 256> t{};
    t.fill('N');
    t['A'] = 'T'; t['C'] = 'G'; t['G'] = 'C'; t['T'] = 'A';
    return t;
}();

// The read as sequenced: a reverse-strand record's sequence reverse-complemented back, its qualities
// reversed (none when the record has none).
std::string sequence_of(const bam1_t* b, std::string* qual = nullptr) {
    const std::size_t n = static_cast<std::size_t>(b->core.l_qseq);
    std::string s(n, 'N');
    const std::uint8_t* q = bam_get_seq(b);
    for (std::size_t i = 0; i + 1 < n; i += 2) {
        const auto& p = kPairs[q[i / 2]];
        s[i] = p[0];
        s[i + 1] = p[1];
    }
    if (n % 2) s[n - 1] = kPairs[q[n / 2]][0];
    const bool reverse = b->core.flag & BAM_FREVERSE;
    if (reverse) {
        std::reverse(s.begin(), s.end());
        for (char& c : s) c = kComplement[static_cast<unsigned char>(c)];
    }
    if (qual != nullptr) {
        const std::uint8_t* bq = bam_get_qual(b);
        qual->clear();
        if (n > 0 && bq[0] != 0xff) {
            qual->resize(n);
            for (std::size_t i = 0; i < n; ++i) (*qual)[i] = static_cast<char>(std::min<int>(bq[i], 93) + 33);
            if (reverse) std::reverse(qual->begin(), qual->end());
        }
    }
    return s;
}

class Alignments final : public FragmentSource {
public:
    Alignments(const std::string& path, const InputOptions& o) : o_(o), path_(path) {
        if (o.threads > 1) {
            pool_.pool = hts_tpool_init(static_cast<int>(o.threads));
            if (pool_.pool == nullptr) throw std::runtime_error("recruit: cannot start htslib threads");
        }
        try {
            file_ = std::make_unique<HtsFile>(path, o, &pool_);
            coordinate_ = file_->coordinate_sorted();
            sample_ = file_->sample();
            if (coordinate_ && !o.screened) throw std::runtime_error("recruit: coordinate-sorted input needs the loci's screen");
        } catch (...) {
            file_.reset();
            if (pool_.pool) hts_tpool_destroy(pool_.pool);
            pool_.pool = nullptr;
            throw;
        }
    }
    ~Alignments() override {
        stop_pipeline();
        for (auto& c : chunks_) for (bam1_t* b : c->recs) bam_destroy1(b);
        file_.reset();
        if (pool_.pool) hts_tpool_destroy(pool_.pool);
    }
    std::size_t next(std::vector<Fragment>& out, std::size_t n) override {
        if (!started_) start_pipeline();
        while (ready_.size() < n && !eof_) {
            if (cur_ != nullptr && cur_pos_ < cur_->n) {
                take(cur_->out[cur_pos_++]);
                continue;
            }
            if (cur_ != nullptr) {
                const bool last = cur_->last;
                const std::exception_ptr error = cur_->error;
                release(cur_);
                cur_ = nullptr;
                if (error) {
                    stop_pipeline();
                    std::rethrow_exception(error);
                }
                if (last) {
                    stop_pipeline();
                    eof_ = true;
                    finish();
                    break;
                }
            }
            cur_ = wait_chunk(next_id_++);
            cur_pos_ = 0;
        }
        const std::size_t k = std::min(n, ready_.size());
        out.resize(k);
        for (std::size_t i = 0; i < k; ++i) {
            out[i] = std::move(ready_.front());
            ready_.pop_front();
        }
        return k;
    }
    const InputStats& stats() const override { return st_; }
    std::string describe() const override {
        const htsFormat* f = hts_get_format(file_->fp);
        const char* kind = f->format == cram ? "CRAM" : f->format == bam ? "BAM" : "SAM";
        return path_ + " (" + kind + (coordinate_ ? ", coordinate-sorted" : "") + (sample_.empty() ? "" : ", sample " + sample_) + ")";
    }

private:
    static constexpr std::int64_t kNear = 100000;   // a mate this far ahead is waited for
    static constexpr std::int64_t kFetchGap = 200000;   // CRAM mates this close are fetched by one query
    static constexpr std::size_t kHoldMax = 2000000;  // reads waiting for mates in an unsorted file before going lean
    struct Read {
        std::string seq, qual;
    };
    Read read_of(const bam1_t* b) const {
        Read r;
        r.seq = sequence_of(b, o_.qualities ? &r.qual : nullptr);
        return r;
    }
    struct Held {
        std::string seq, qual;
        bool read1 = true;
        bool candidate = false;   // kept because it hits a screen (else waiting for a near mate, unscreened)
        std::int32_t mtid = -1;
        std::int64_t mpos = -1;
    };
    struct Wanted {   // a kept read whose mate must be fetched
        std::string name, seq, qual;
        bool read1 = true;
        std::int32_t mtid = -1;
        std::int64_t mpos = -1;
    };

    // A record's fields that take() uses, prepared off the pairing thread.
    struct Rec {
        std::uint16_t flag = 0;
        std::int32_t tid = -1, mtid = -1;
        std::int64_t pos = -1, mpos = -1;
        std::string name, seq, qual;
    };
    void prepare(const bam1_t* b, Rec& r) const {
        r.flag = b->core.flag;
        r.tid = b->core.tid;
        r.pos = b->core.pos;
        r.mtid = b->core.mtid;
        r.mpos = b->core.mpos;
        r.name.assign(bam_get_qname(b));
        r.qual.clear();
        if (r.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL)) { r.seq.clear(); return; }   // skipped by take()
        r.seq = sequence_of(b, o_.qualities ? &r.qual : nullptr);
    }

    // The reading pipeline: one thread reads records into chunks (htslib decodes on its own threads), workers prepare
    // them (names, sequences, qualities) in parallel, and next() takes them in file order, so the pairing sees the
    // records exactly as one thread reading them would.
    struct Chunk {
        std::uint64_t id = 0;
        std::vector<bam1_t*> recs;
        std::size_t n = 0;
        bool last = false;               // the file ends with this chunk
        std::exception_ptr error;        // reading failed after its n records
        std::vector<Rec> out;
    };
    static constexpr std::size_t kChunk = 4096;
    void start_pipeline() {
        started_ = true;
        const std::size_t workers = std::max<std::size_t>(1, o_.parse_threads);
        for (std::size_t i = 0; i < 2 * workers + 4; ++i) {
            auto c = std::make_unique<Chunk>();
            c->recs.resize(kChunk);
            for (bam1_t*& b : c->recs) b = bam_init1();
            c->out.resize(kChunk);
            free_.push_back(c.get());
            chunks_.push_back(std::move(c));
        }
        producer_ = std::make_unique<BigStackThread>([this] { produce(); });
        for (std::size_t i = 0; i < workers; ++i) workers_.push_back(std::make_unique<BigStackThread>([this] { work(); }));
    }
    void stop_pipeline() {
        if (!started_ || stopped_) return;
        {
            std::lock_guard<std::mutex> lock(mu_);
            stop_ = true;
        }
        cv_.notify_all();
        if (producer_) producer_->join();
        for (auto& w : workers_) w->join();
        stopped_ = true;
    }
    void produce() {
        for (std::uint64_t id = 0;; ++id) {
            Chunk* c = nullptr;
            {
                std::unique_lock<std::mutex> lock(mu_);
                cv_.wait(lock, [&] { return stop_ || !free_.empty(); });
                if (stop_) break;
                c = free_.front();
                free_.pop_front();
            }
            c->id = id;
            c->n = 0;
            c->last = false;
            c->error = nullptr;
            try {
                while (c->n < kChunk) {
                    if (!file_->read(c->recs[c->n])) { c->last = true; break; }
                    ++c->n;
                }
            } catch (...) {
                c->error = std::current_exception();
                c->last = true;
            }
            {
                std::lock_guard<std::mutex> lock(mu_);
                todo_.push_back(c);
                if (c->last) produced_ = true;
            }
            cv_.notify_all();
            if (c->last) break;
        }
        std::lock_guard<std::mutex> lock(mu_);
        produced_ = true;
        cv_.notify_all();
    }
    void work() {
        for (;;) {
            Chunk* c = nullptr;
            {
                std::unique_lock<std::mutex> lock(mu_);
                cv_.wait(lock, [&] { return stop_ || !todo_.empty() || produced_; });
                if (todo_.empty()) {
                    if (stop_ || produced_) return;
                    continue;
                }
                c = todo_.front();
                todo_.pop_front();
            }
            for (std::size_t i = 0; i < c->n; ++i) prepare(c->recs[i], c->out[i]);
            {
                std::lock_guard<std::mutex> lock(mu_);
                done_[c->id] = c;
            }
            cv_.notify_all();
        }
    }
    Chunk* wait_chunk(std::uint64_t id) {
        std::unique_lock<std::mutex> lock(mu_);
        cv_.wait(lock, [&] { return done_.count(id) > 0; });
        Chunk* c = done_[id];
        done_.erase(id);
        return c;
    }
    void release(Chunk* c) {
        {
            std::lock_guard<std::mutex> lock(mu_);
            free_.push_back(c);
        }
        cv_.notify_all();
    }

    // Skip, count and turn a record into a fragment, or hold it for its mate.
    void take(Rec& r) {
        ++st_.records;
        const std::uint16_t flag = r.flag;
        if (flag & BAM_FSECONDARY) { ++st_.secondary; return; }
        if (flag & BAM_FSUPPLEMENTARY) { ++st_.supplementary; return; }
        if (flag & BAM_FQCFAIL) { ++st_.qc_fail; return; }
        if (flag & BAM_FDUP) ++st_.duplicates;
        if (flag & BAM_FUNMAP) ++st_.unmapped;
        std::string name = std::move(r.name);
        if (!(flag & BAM_FPAIRED)) {
            emit(std::move(name), Read{std::move(r.seq), std::move(r.qual)}, Read{}, false);
            return;
        }
        const bool read1 = !(flag & BAM_FREAD2);
        if (coordinate_) {
            check_order(r.tid, r.pos, name);
            if (r.tid != last_tid_) {
                sweep(r.tid, -1);
                last_tid_ = r.tid;
            } else if (st_.records % 65536 == 0) {
                sweep(r.tid, r.pos);
            }
        }
        const auto it = held_.find(name);
        if (it != held_.end()) {
            if (it->second.read1 != read1) {
                pair_up(std::move(name), std::move(it->second), Read{std::move(r.seq), std::move(r.qual)}, read1);
                held_.erase(it);
                return;
            }
            // a second read 1 (or read 2) of this name: not its mate; the held one goes on alone
            ++st_.unmatched;
            emit(it->first, Read{std::move(it->second.seq), std::move(it->second.qual)}, Read{}, false);
            held_.erase(it);
        }
        Held h;
        h.seq = std::move(r.seq);
        h.qual = std::move(r.qual);
        h.read1 = read1;
        h.mtid = r.mtid;
        h.mpos = r.mpos;
        if (!coordinate_) {
            if (lean_) {   // holding too much already: only a read that hits a screen is kept, its mate fetched later
                if (o_.screened(h.seq)) {
                    h.candidate = true;
                    h.mtid = -1;
                    h.mpos = -1;
                    held_.emplace(std::move(name), std::move(h));
                } else {
                    ++lean_dropped_;
                }
                return;
            }
            held_.emplace(std::move(name), std::move(h));
            if (held_.size() > hold_max_) go_lean();
            return;
        }
        const std::int32_t tid = r.tid;
        const std::int64_t pos = r.pos;
        // a mate without a known place (unplaced, or RNEXT/PNEXT unavailable) is waited for like a near one
        const bool mate_unplaced = h.mtid < 0 || h.mpos < 0 || tid < 0;
        const bool behind = !mate_unplaced && (h.mtid < tid || (h.mtid == tid && h.mpos < pos));
        if (behind) {   // its mate passed and was not kept: it did not hit; this one decides the pair
            if (o_.screened(h.seq)) wanted_.push_back({std::move(name), std::move(h.seq), std::move(h.qual), read1, h.mtid, h.mpos});
            else ++st_.unscreened;
            return;
        }
        const bool near = mate_unplaced || (h.mtid == tid && h.mpos - pos <= near_);
        if (near) {
            held_.emplace(std::move(name), std::move(h));
        } else if (o_.screened(h.seq)) {
            h.candidate = true;
            held_.emplace(std::move(name), std::move(h));
        }
    }

    // Pairing relies on the order the header claims: a record before the last one stops the run.
    void check_order(std::int32_t tid, std::int64_t pos, const std::string& qname) {
        if (tid < 0) { unplaced_seen_ = true; return; }
        if (unplaced_seen_ || tid < order_tid_ || (tid == order_tid_ && pos < order_pos_))
            throw std::runtime_error("recruit: " + path_ + " is not sorted by coordinate (" + qname +
                                     " comes after a later position), though it is read as sorted");
        order_tid_ = tid;
        order_pos_ = pos;
    }

    // Reads waiting for a mate whose place has passed lost it (filtered out): kept if they hit. Called at a
    // new contig (pos -1: any mate on another contig has passed) and now and then within one.
    void sweep(std::int32_t tid, std::int64_t pos) {
        for (auto it = held_.begin(); it != held_.end();) {
            Held& h = it->second;
            const bool known = h.mtid >= 0 && h.mpos >= 0;
            const bool passed = pos < 0 ? h.mtid != tid : (h.mtid < tid || (h.mtid == tid && h.mpos < pos));
            if (!h.candidate && known && passed) {
                if (o_.screened(h.seq)) wanted_.push_back({it->first, std::move(h.seq), std::move(h.qual), h.read1, h.mtid, h.mpos});
                else ++st_.dropped;
                it = held_.erase(it);
            } else {
                ++it;
            }
        }
    }

    // The end: every read still waiting either never met its mate (unsorted input: scored alone) or, in
    // coordinate-sorted input, has its mate fetched.
    // An unsorted file whose mates are far apart (all read 1 before all read 2, say) would hold most of the
    // sample: past kHoldMax waiting reads, only those hitting a screen are kept (as candidates, their mates
    // fetched by a second pass) and new reads are kept only if they hit. A pair is still scored when
    // either mate hits.
    void go_lean() {
        lean_ = true;
        if (o_.log) o_.log("input: " + path_ + "; reads waiting for their mates: over " + std::to_string(hold_max_) +
                           " (the file is not name-grouped): keeping only reads that hit a screen, their mates fetched by a second pass");
        for (auto it = held_.begin(); it != held_.end();) {
            if (o_.screened(it->second.seq)) {
                it->second.candidate = true;
                it->second.mtid = -1;
                it->second.mpos = -1;
                ++it;
            } else {
                ++lean_dropped_;
                it = held_.erase(it);
            }
        }
    }

    void finish() {
        if (!coordinate_ && !lean_) {
            for (auto& [name, h] : held_) {
                ++st_.unmatched;
                emit(name, Read{std::move(h.seq), std::move(h.qual)}, Read{}, false);
            }
            held_.clear();
            return;
        }
        for (auto& [name, h] : held_) {
            if (h.candidate || o_.screened(h.seq)) wanted_.push_back({name, std::move(h.seq), std::move(h.qual), h.read1, h.mtid, h.mpos});
            else ++st_.dropped;
        }
        held_.clear();
        // both mates waiting (stale mate fields): they are each other's mate
        std::unordered_map<std::string, std::size_t> first;
        std::vector<Wanted> rest;
        std::vector<char> paired(wanted_.size(), 0);
        for (std::size_t i = 0; i < wanted_.size(); ++i) {
            const auto [it, fresh] = first.emplace(wanted_[i].name, i);
            if (fresh) continue;
            Wanted& a = wanted_[it->second];
            Wanted& b = wanted_[i];
            if (paired[it->second] || a.read1 == b.read1) continue;
            Held h;
            h.seq = std::move(a.seq);
            h.qual = std::move(a.qual);
            h.read1 = a.read1;
            pair_up(b.name, std::move(h), Read{std::move(b.seq), std::move(b.qual)}, b.read1);
            paired[it->second] = paired[i] = 1;
        }
        for (std::size_t i = 0; i < wanted_.size(); ++i)
            if (!paired[i]) rest.push_back(std::move(wanted_[i]));
        wanted_.swap(rest);
        if (!wanted_.empty()) fetch();
        st_.unscreened += lean_dropped_ / 2;   // reads set aside whose mates were set aside too: unscreened pairs
        st_.dropped += lean_dropped_ % 2;
    }

    bool mate_of(const bam1_t* b, const Wanted& w) const {
        const std::uint16_t f = b->core.flag;
        if (f & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL)) return false;
        return static_cast<bool>(f & BAM_FREAD2) == w.read1 && w.name == bam_get_qname(b);
    }

    void fetch() {
        std::vector<char> found(wanted_.size(), 0);
        std::vector<Read> mates(wanted_.size());
        bool indexed = false;
        {
            HtsFile probe(path_, o_, nullptr);
            hts_idx_t* idx = sam_index_load(probe.fp, path_.c_str());
            indexed = idx != nullptr;
            if (idx) hts_idx_destroy(idx);
        }
        if (indexed) {
            // in the mates' order on the genome. A CRAM reads nearby mates by one query, so each reference contig
            // is loaded once (in the order the reads came, it reloads a contig at almost every mate); a BAM queries
            // each mate position alone, as a query decodes every record between its ends.
            const bool is_cram = hts_get_format(file_->fp)->format == cram;
            const std::int64_t gap = is_cram ? kFetchGap : 0;
            std::vector<std::size_t> order;
            for (std::size_t i = 0; i < wanted_.size(); ++i)
                if (wanted_[i].mtid >= 0 && wanted_[i].mpos >= 0) order.push_back(i);
            std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
                return std::tie(wanted_[a].mtid, wanted_[a].mpos, a) < std::tie(wanted_[b].mtid, wanted_[b].mpos, b);
            });
            if (o_.log && !order.empty()) o_.log("input: " + path_ + "; mates to fetch through the index: " + std::to_string(order.size()));
            std::vector<std::pair<std::size_t, std::size_t>> queries;   // ranges of `order`, one query each
            for (std::size_t a = 0; a < order.size();) {
                std::size_t b = a + 1;
                while (b < order.size() && wanted_[order[b]].mtid == wanted_[order[a]].mtid &&
                       wanted_[order[b]].mpos - wanted_[order[b - 1]].mpos <= gap)
                    ++b;
                queries.emplace_back(a, b);
                a = b;
            }
            // the queries are independent: each thread takes a run of them in genome order, with its own file and
            // index, and fills the slots of its mates; they are delivered in the reads' order once all are done
            const std::size_t n = std::max<std::size_t>(1, std::min(o_.fetch_threads, queries.size()));
            std::vector<std::size_t> cut{0};
            for (std::size_t t = 1; t < n; ++t) {   // about as many mates per thread
                std::size_t q = cut.back();
                while (q < queries.size() && queries[q].first < order.size() * t / n) ++q;
                cut.push_back(q);
            }
            cut.push_back(queries.size());
            std::vector<std::exception_ptr> errors(n);
            auto work = [&](std::size_t t) {
                try {
                    HtsFile f(path_, o_, nullptr);
                    hts_idx_t* idx = sam_index_load(f.fp, path_.c_str());
                    if (idx == nullptr) throw std::runtime_error("recruit: " + path_ + ": cannot load the index");
                    std::unordered_map<std::string, std::vector<std::size_t>> group;   // a query's wanted reads, by name
                    for (std::size_t q = cut[t]; q < cut[t + 1]; ++q) {
                        const auto [a, b] = queries[q];
                        group.clear();
                        for (std::size_t k = a; k < b; ++k) group[wanted_[order[k]].name].push_back(order[k]);
                        std::size_t left = b - a;
                        hts_itr_t* itr = sam_itr_queryi(idx, wanted_[order[a]].mtid, wanted_[order[a]].mpos, wanted_[order[b - 1]].mpos + 1);
                        if (itr == nullptr) continue;
                        int r = 0;
                        while (left > 0 && (r = sam_itr_next(f.fp, itr, f.rec)) >= 0) {
                            const auto it = group.find(bam_get_qname(f.rec));
                            if (it == group.end()) continue;
                            for (const std::size_t i : it->second)
                                if (!found[i] && f.rec->core.pos == wanted_[i].mpos && mate_of(f.rec, wanted_[i])) {
                                    found[i] = 1;
                                    mates[i] = read_of(f.rec);
                                    --left;
                                    break;
                                }
                        }
                        hts_itr_destroy(itr);
                        if (left > 0 && r < -1) {
                            hts_idx_destroy(idx);
                            throw std::runtime_error("recruit: " + path_ + ": cannot read a record while fetching mates");
                        }
                    }
                    hts_idx_destroy(idx);
                } catch (...) {
                    errors[t] = std::current_exception();
                }
            };
            std::vector<std::unique_ptr<BigStackThread>> pool;
            for (std::size_t t = 1; t < n; ++t) pool.push_back(std::make_unique<BigStackThread>([&work, t] { work(t); }));
            work(0);
            for (auto& th : pool) th->join();
            for (const auto& e : errors) if (e) std::rethrow_exception(e);
            for (std::size_t i = 0; i < wanted_.size(); ++i)   // in the order the reads came, as before
                if (found[i]) deliver(wanted_[i], std::move(mates[i]));
        }
        std::unordered_map<std::string, std::size_t> rest;   // a second pass for mates no index found
        for (std::size_t i = 0; i < wanted_.size(); ++i)
            if (!found[i]) rest.emplace(wanted_[i].name, i);   // no index, no mate place, or stale mate fields
        if (!rest.empty()) {
            if (o_.log) o_.log("input: " + path_ + "; second pass for mates: " + std::to_string(rest.size()) + (indexed ? " (not where the index put them)" : " (no index)"));
            HtsFile pass(path_, o_, &pool_);
            while (!rest.empty() && pass.read()) {
                const auto it = rest.find(bam_get_qname(pass.rec));
                if (it == rest.end() || !mate_of(pass.rec, wanted_[it->second])) continue;
                found[it->second] = 1;
                deliver(wanted_[it->second], read_of(pass.rec));
                rest.erase(it);
            }
        }
        for (std::size_t i = 0; i < wanted_.size(); ++i) {
            if (found[i]) continue;
            Wanted& w = wanted_[i];
            ++st_.unmatched;
            emit(std::move(w.name), Read{std::move(w.seq), std::move(w.qual)}, Read{}, false);
        }
        wanted_.clear();
    }

    void deliver(const Wanted& w, Read mate) {
        ++st_.fetched;
        if (lean_ && lean_dropped_ > 0) --lean_dropped_;   // the mate fetched was one set aside
        Held h;
        h.seq = w.seq;
        h.qual = w.qual;
        h.read1 = w.read1;
        pair_up(w.name, std::move(h), std::move(mate), !w.read1);
    }

    void pair_up(std::string name, Held first, Read second, bool read1) {
        // mate 1 first: by flags, else in arrival order
        Read held{std::move(first.seq), std::move(first.qual)};
        if (first.read1 || !read1) emit(std::move(name), std::move(held), std::move(second), true);
        else emit(std::move(name), std::move(second), std::move(held), true);
    }

    void emit(std::string name, Read r1, Read r2, bool paired) {
        Fragment f;
        f.name = std::move(name);
        f.seq1 = std::move(r1.seq);
        f.qual1 = std::move(r1.qual);
        f.seq2 = std::move(r2.seq);
        f.qual2 = std::move(r2.qual);
        f.paired = paired;
        ++st_.fragments;
        if (paired) ++st_.pairs;
        ready_.push_back(std::move(f));
    }

    InputOptions o_;
    std::string path_;
    htsThreadPool pool_{nullptr, 0};
    std::unique_ptr<HtsFile> file_;
    std::string sample_;
    bool coordinate_ = false, eof_ = false, lean_ = false;
    std::uint64_t lean_dropped_ = 0;
    std::size_t hold_max_ = std::getenv("PANVAR_TEST_HOLD_MAX") ? std::strtoull(std::getenv("PANVAR_TEST_HOLD_MAX"), nullptr, 10) : kHoldMax;
    std::int64_t near_ = std::getenv("PANVAR_TEST_NEAR") ? std::strtoll(std::getenv("PANVAR_TEST_NEAR"), nullptr, 10) : kNear;
    std::int32_t last_tid_ = -2;     // the contig of the last record (sweeps)
    std::int32_t order_tid_ = -1;    // the last placed record's position (the order check)
    std::int64_t order_pos_ = -1;
    bool unplaced_seen_ = false;
    std::unordered_map<std::string, Held> held_;
    std::vector<Wanted> wanted_;
    std::deque<Fragment> ready_;
    InputStats st_;
    // the pipeline
    std::mutex mu_;
    std::condition_variable cv_;
    std::vector<std::unique_ptr<Chunk>> chunks_;
    std::deque<Chunk*> free_, todo_;
    std::map<std::uint64_t, Chunk*> done_;
    std::unique_ptr<BigStackThread> producer_;
    std::vector<std::unique_ptr<BigStackThread>> workers_;
    bool started_ = false, stopped_ = false, stop_ = false, produced_ = false;
    std::uint64_t next_id_ = 0;
    Chunk* cur_ = nullptr;
    std::size_t cur_pos_ = 0;
};

}  // namespace

std::unique_ptr<FragmentSource> open_fastq_pair(const std::string& reads1, const std::string& reads2) {
    return std::make_unique<FastqPair>(reads1, reads2);
}

bool is_alignment_file(const std::string& path) {
    samFile* fp = hts_open(path.c_str(), "r");
    if (fp == nullptr) return false;
    const htsFormat* f = hts_get_format(fp);
    const bool yes = f->format == sam || f->format == bam || f->format == cram;
    hts_close(fp);
    return yes;
}

std::unique_ptr<FragmentSource> open_input(const std::string& path, const InputOptions& options) {
    if (!std::filesystem::is_regular_file(path)) throw std::runtime_error("recruit: no file at " + path);
    if (is_alignment_file(path)) return std::make_unique<Alignments>(path, options);
    return std::make_unique<FastqSingle>(path, options.interleaved);
}

InputStats count_records(const std::string& path, const InputOptions& options) {
    if (!is_alignment_file(path)) {
        FastqSingle f(path, false);
        std::vector<Fragment> batch;
        while (f.next(batch, 100000) > 0) {}
        InputStats st = f.stats();
        if (options.interleaved) st.pairs = st.records, st.fragments = 0;
        return st;
    }
    htsThreadPool pool{nullptr, 0};
    if (options.threads > 1) pool.pool = hts_tpool_init(static_cast<int>(options.threads));
    InputStats st;
    {
        HtsFile file(path, options, &pool);
        while (file.read()) {
            const std::uint16_t flag = file.rec->core.flag;
            ++st.records;
            if (flag & BAM_FSECONDARY) ++st.secondary;
            else if (flag & BAM_FSUPPLEMENTARY) ++st.supplementary;
            else if (flag & BAM_FQCFAIL) ++st.qc_fail;
            else {
                if (flag & BAM_FDUP) ++st.duplicates;
                if (flag & BAM_FUNMAP) ++st.unmapped;
                if (flag & BAM_FPAIRED) ++st.pairs;   // here: paired reads
                else ++st.fragments;                  // here: single reads
            }
        }
    }
    if (pool.pool) hts_tpool_destroy(pool.pool);
    return st;
}

}  // namespace panvar::recruit
