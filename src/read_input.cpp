#include "panvar/read_input.hpp"

#include "panvar/gz_reader.hpp"

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/thread_pool.h>

#include <algorithm>
#include <array>
#include <deque>
#include <exception>
#include <filesystem>
#include <stdexcept>
#include <thread>
#include <unordered_map>

namespace panvar::recruit {

namespace {

bool read_fastq(GzLineReader& in, std::string& name, std::string& seq) {
    std::string h, plus, qual;
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
                while (n2 < n && read_fastq(*in2_, names2_[n2], out[n2].seq2)) ++n2;
            } catch (...) {
                failed2 = std::current_exception();
            }
        });
        try {
            while (n1 < n && read_fastq(*in1_, out[n1].name, out[n1].seq1)) out[n1++].paired = true;
        } catch (...) {
            failed1 = std::current_exception();
        }
        mates.join();
        if (failed1) std::rethrow_exception(failed1);
        if (failed2) std::rethrow_exception(failed2);
        if (n2 < n1) throw std::runtime_error("recruit: " + r2_ + " has fewer reads than " + r1_);
        if (n2 > n1) throw std::runtime_error("recruit: " + r2_ + " has more reads than " + r1_);
        if (n1 < n) {   // the first file ended: so must the second
            std::string name, seq;
            if (read_fastq(*in2_, name, seq)) throw std::runtime_error("recruit: " + r2_ + " has more reads than " + r1_);
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
        while (k < n && read_fastq(*in_, out[k].name, out[k].seq1)) {
            ++st_.records;
            out[k].paired = interleaved_;
            out[k].seq2.clear();
            if (interleaved_) {
                if (!read_fastq(*in_, name2, out[k].seq2)) throw std::runtime_error("recruit: " + path_ + ": an odd number of interleaved records");
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
        if (pool != nullptr && pool->pool != nullptr) hts_set_thread_pool(fp, pool);
        if (f->format == cram) {
            if (!o.cram_reference.empty() && hts_set_fai_filename(fp, o.cram_reference.c_str()) != 0)
                throw std::runtime_error("recruit: " + p + ": cannot use the CRAM reference " + o.cram_reference);
            hts_set_opt(fp, CRAM_OPT_REQUIRED_FIELDS, SAM_QNAME | SAM_FLAG | SAM_RNAME | SAM_POS | SAM_RNEXT | SAM_PNEXT | SAM_SEQ);
        }
        hdr = sam_hdr_read(fp);
        if (hdr == nullptr) throw std::runtime_error("recruit: " + p + ": cannot read the header");
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
    bool read() {
        const int r = sam_read1(fp, hdr, rec);
        if (r >= 0) return true;
        if (r == -1) return false;
        throw std::runtime_error("recruit: " + path + ": cannot read a record (truncated or damaged, or a CRAM whose reference cannot be found; see --cram-reference)");
    }
    bool coordinate_sorted() const {
        kstring_t ks = KS_INITIALIZE;
        const bool yes = sam_hdr_find_tag_hd(hdr, "SO", &ks) == 0 && ks.s != nullptr && std::string(ks.s) == "coordinate";
        ks_free(&ks);
        return yes;
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

// The read as sequenced: a reverse-strand record's sequence reverse-complemented back.
std::string sequence_of(const bam1_t* b) {
    const std::size_t n = static_cast<std::size_t>(b->core.l_qseq);
    std::string s(n, 'N');
    const std::uint8_t* q = bam_get_seq(b);
    for (std::size_t i = 0; i + 1 < n; i += 2) {
        const auto& p = kPairs[q[i / 2]];
        s[i] = p[0];
        s[i + 1] = p[1];
    }
    if (n % 2) s[n - 1] = kPairs[q[n / 2]][0];
    if (b->core.flag & BAM_FREVERSE) {
        std::reverse(s.begin(), s.end());
        for (char& c : s) c = kComplement[static_cast<unsigned char>(c)];
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
        file_ = std::make_unique<HtsFile>(path, o, &pool_);
        coordinate_ = file_->coordinate_sorted();
        if (coordinate_ && !o.screened) throw std::runtime_error("recruit: coordinate-sorted input needs the loci's screen");
    }
    ~Alignments() override {
        file_.reset();
        if (pool_.pool) hts_tpool_destroy(pool_.pool);
    }
    std::size_t next(std::vector<Fragment>& out, std::size_t n) override {
        while (ready_.size() < n && !eof_) {
            if (!file_->read()) {
                eof_ = true;
                finish();
                break;
            }
            take(file_->rec);
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
        return path_ + " (" + kind + (coordinate_ ? ", coordinate-sorted" : "") + ")";
    }

private:
    static constexpr std::int64_t kNear = 100000;   // a mate this far ahead is waited for
    struct Held {
        std::string seq;
        bool read1 = true;
        bool candidate = false;   // kept because it hits a screen (else waiting for a near mate, unscreened)
        std::int32_t mtid = -1;
        std::int64_t mpos = -1;
    };
    struct Wanted {   // a kept read whose mate must be fetched
        std::string name, seq;
        bool read1 = true;
        std::int32_t mtid = -1;
        std::int64_t mpos = -1;
    };

    // Skip, count and turn a record into a fragment, or hold it for its mate.
    void take(const bam1_t* b) {
        ++st_.records;
        const std::uint16_t flag = b->core.flag;
        if (flag & BAM_FSECONDARY) { ++st_.secondary; return; }
        if (flag & BAM_FSUPPLEMENTARY) { ++st_.supplementary; return; }
        if (flag & BAM_FQCFAIL) { ++st_.qc_fail; return; }
        if (flag & BAM_FDUP) ++st_.duplicates;
        if (flag & BAM_FUNMAP) ++st_.unmapped;
        std::string name = bam_get_qname(b);
        if (!(flag & BAM_FPAIRED)) {
            emit(std::move(name), sequence_of(b), std::string(), false);
            return;
        }
        const bool read1 = !(flag & BAM_FREAD2);
        if (coordinate_ && b->core.tid != last_tid_) {
            sweep(b->core.tid);
            last_tid_ = b->core.tid;
        }
        const auto it = held_.find(name);
        if (it != held_.end()) {
            pair_up(std::move(name), std::move(it->second), sequence_of(b), read1);
            held_.erase(it);
            return;
        }
        Held h;
        h.seq = sequence_of(b);
        h.read1 = read1;
        h.mtid = b->core.mtid;
        h.mpos = b->core.mpos;
        if (!coordinate_) {
            held_.emplace(std::move(name), std::move(h));
            return;
        }
        const std::int32_t tid = b->core.tid;
        const std::int64_t pos = b->core.pos;
        const bool mate_unplaced = h.mtid < 0 || tid < 0;
        const bool behind = !mate_unplaced && (h.mtid < tid || (h.mtid == tid && h.mpos < pos));
        if (behind) {   // its mate passed and was not kept: it did not hit; this one decides the pair
            if (o_.screened(h.seq)) wanted_.push_back({std::move(name), std::move(h.seq), read1, h.mtid, h.mpos});
            else ++st_.unscreened;
            return;
        }
        const bool near = mate_unplaced || (h.mtid == tid && h.mpos - pos <= kNear);
        if (near) {
            held_.emplace(std::move(name), std::move(h));
        } else if (o_.screened(h.seq)) {
            h.candidate = true;
            held_.emplace(std::move(name), std::move(h));
        }
    }

    // A new contig: reads waiting for a mate on a passed contig lost it (filtered out): kept if they hit.
    void sweep(std::int32_t tid) {
        for (auto it = held_.begin(); it != held_.end();) {
            Held& h = it->second;
            if (!h.candidate && h.mtid >= 0 && h.mtid != tid) {
                if (o_.screened(h.seq)) wanted_.push_back({it->first, std::move(h.seq), h.read1, h.mtid, h.mpos});
                else ++st_.dropped;
                it = held_.erase(it);
            } else {
                ++it;
            }
        }
    }

    // The end: every read still waiting either never met its mate (unsorted input: scored alone) or, in
    // coordinate-sorted input, has its mate fetched.
    void finish() {
        if (!coordinate_) {
            for (auto& [name, h] : held_) {
                ++st_.unmatched;
                if (h.read1) emit(name, std::move(h.seq), std::string(), true);
                else emit(name, std::string(), std::move(h.seq), true);
            }
            held_.clear();
            return;
        }
        for (auto& [name, h] : held_) {
            if (h.candidate || o_.screened(h.seq)) wanted_.push_back({name, std::move(h.seq), h.read1, h.mtid, h.mpos});
            else ++st_.dropped;
        }
        held_.clear();
        if (wanted_.empty()) return;
        fetch();
    }

    bool mate_of(const bam1_t* b, const Wanted& w) const {
        const std::uint16_t f = b->core.flag;
        if (f & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY | BAM_FQCFAIL)) return false;
        return static_cast<bool>(f & BAM_FREAD2) == w.read1 && w.name == bam_get_qname(b);
    }

    void fetch() {
        std::vector<char> found(wanted_.size(), 0);
        HtsFile again(path_, o_, &pool_);
        hts_idx_t* idx = sam_index_load(again.fp, path_.c_str());
        if (idx != nullptr) {
            for (std::size_t i = 0; i < wanted_.size(); ++i) {
                const Wanted& w = wanted_[i];
                if (w.mtid < 0) continue;
                hts_itr_t* itr = sam_itr_queryi(idx, w.mtid, w.mpos, w.mpos + 1);
                if (itr == nullptr) continue;
                int r;
                while ((r = sam_itr_next(again.fp, itr, again.rec)) >= 0)
                    if (again.rec->core.pos == w.mpos && mate_of(again.rec, w)) { found[i] = 1; deliver(w, again.rec); break; }
                hts_itr_destroy(itr);
                if (r < -1) { hts_idx_destroy(idx); throw std::runtime_error("recruit: " + path_ + ": cannot read a record while fetching mates"); }
            }
            hts_idx_destroy(idx);
        }
        std::unordered_map<std::string, std::size_t> rest;   // a second pass for mates no index found
        for (std::size_t i = 0; i < wanted_.size(); ++i)
            if (!found[i] && (idx == nullptr || wanted_[i].mtid < 0)) rest.emplace(wanted_[i].name, i);
        if (!rest.empty()) {
            if (o_.log) o_.log(path_ + ": a second pass for " + std::to_string(rest.size()) + " mates" + (idx == nullptr ? " (no index)" : ""));
            HtsFile pass(path_, o_, &pool_);
            while (!rest.empty() && pass.read()) {
                const auto it = rest.find(bam_get_qname(pass.rec));
                if (it == rest.end() || !mate_of(pass.rec, wanted_[it->second])) continue;
                found[it->second] = 1;
                deliver(wanted_[it->second], pass.rec);
                rest.erase(it);
            }
        }
        for (std::size_t i = 0; i < wanted_.size(); ++i) {
            if (found[i]) continue;
            Wanted& w = wanted_[i];
            ++st_.unmatched;
            if (w.read1) emit(std::move(w.name), std::move(w.seq), std::string(), true);
            else emit(std::move(w.name), std::string(), std::move(w.seq), true);
        }
        wanted_.clear();
    }

    void deliver(const Wanted& w, const bam1_t* mate) {
        ++st_.fetched;
        Held h;
        h.seq = w.seq;
        h.read1 = w.read1;
        pair_up(w.name, std::move(h), sequence_of(mate), !w.read1);
    }

    void pair_up(std::string name, Held first, std::string seq, bool read1) {
        // mate 1 first: by flags, else in arrival order
        if (first.read1 || !read1) emit(std::move(name), std::move(first.seq), std::move(seq), true);
        else emit(std::move(name), std::move(seq), std::move(first.seq), true);
    }

    void emit(std::string name, std::string s1, std::string s2, bool paired) {
        Fragment f;
        f.name = std::move(name);
        f.seq1 = std::move(s1);
        f.seq2 = std::move(s2);
        f.paired = paired;
        ++st_.fragments;
        if (paired) ++st_.pairs;
        ready_.push_back(std::move(f));
    }

    InputOptions o_;
    std::string path_;
    htsThreadPool pool_{nullptr, 0};
    std::unique_ptr<HtsFile> file_;
    bool coordinate_ = false, eof_ = false;
    std::int32_t last_tid_ = -2;
    std::unordered_map<std::string, Held> held_;
    std::vector<Wanted> wanted_;
    std::deque<Fragment> ready_;
    InputStats st_;
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

std::string sample_of(const std::string& path) {
    samFile* fp = sam_open(path.c_str(), "r");
    if (fp == nullptr) return "";
    sam_hdr_t* hdr = sam_hdr_read(fp);
    std::string sample;
    bool agree = true;
    if (hdr != nullptr) {
        const int n = sam_hdr_count_lines(hdr, "RG");
        for (int i = 0; i < n && agree; ++i) {
            kstring_t ks = KS_INITIALIZE;
            if (sam_hdr_find_tag_pos(hdr, "RG", i, "SM", &ks) == 0 && ks.s != nullptr) {
                if (sample.empty()) sample = ks.s;
                else if (sample != ks.s) agree = false;
            }
            ks_free(&ks);
        }
        sam_hdr_destroy(hdr);
    }
    sam_close(fp);
    return agree ? sample : "";
}

InputStats count_records(const std::string& path, const InputOptions& options) {
    if (!is_alignment_file(path)) {
        FastqSingle f(path, false);
        std::vector<Fragment> batch;
        while (f.next(batch, 100000) > 0) {}
        return f.stats();
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
                if (flag & BAM_FPAIRED) ++st.pairs;   // mates, here: each pair counts twice
                else ++st.fragments;
            }
        }
    }
    if (pool.pool) hts_tpool_destroy(pool.pool);
    st.pairs /= 2;
    st.fragments += st.pairs;
    return st;
}

}  // namespace panvar::recruit
