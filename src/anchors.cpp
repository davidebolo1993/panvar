#include "panvar/anchors.hpp"

#include "panvar/background.hpp"
#include "panvar/gz_reader.hpp"
#include "panvar/kmer.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <cstdio>
#include <cstring>
#include <deque>
#include <filesystem>
#include <functional>
#include <map>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <utility>

#include <fcntl.h>
#include <unistd.h>

namespace panvar::anchors {

namespace {

namespace fs = std::filesystem;

constexpr std::uint64_t kSampleSeed = 0x9e3779b97f4a7c15ULL;
constexpr std::uint64_t kRefMagic = 0x3246455252434E41ULL;   // "ANCRREF2"
constexpr std::uint64_t kSegMagic = 0x3247455352434E41ULL;   // "ANCRSEG2"
constexpr std::size_t kPlaceBytes = 12;                      // genome u16, flags u16, contig u32, pos u32
constexpr std::uint16_t kForward = 1, kRepeated = 2;

// Streams a FASTA's contigs to `threads` workers: f(contig index, sequence). Contig names and lengths are
// collected in file order.
void for_each_contig(const std::string& fasta, std::size_t threads, std::vector<std::string>& names, std::vector<std::uint64_t>& lengths,
                     const std::function<void(std::uint32_t, const std::string&)>& f) {
    GzLineReader in(fasta);
    if (!in.ok()) throw std::runtime_error("anchors: cannot read " + fasta);
    std::deque<std::pair<std::uint32_t, std::string>> queue;
    std::mutex m;
    std::condition_variable not_empty, not_full;
    bool closed = false;
    std::exception_ptr failure;
    std::vector<std::thread> pool;
    for (std::size_t t = 0; t < std::max<std::size_t>(1, threads); ++t)
        pool.emplace_back([&] {
            while (true) {
                std::pair<std::uint32_t, std::string> job;
                {
                    std::unique_lock<std::mutex> lk(m);
                    not_empty.wait(lk, [&] { return !queue.empty() || closed; });
                    if (queue.empty()) return;
                    job = std::move(queue.front());
                    queue.pop_front();
                    not_full.notify_one();
                }
                try {
                    f(job.first, job.second);
                } catch (...) {
                    std::lock_guard<std::mutex> lk(m);
                    if (!failure) failure = std::current_exception();
                    closed = true;
                    queue.clear();
                    not_full.notify_all();
                }
            }
        });
    auto push = [&](std::uint32_t id, std::string seq) {
        std::unique_lock<std::mutex> lk(m);
        not_full.wait(lk, [&] { return queue.size() < std::max<std::size_t>(1, threads) || failure; });
        if (failure) return false;
        queue.emplace_back(id, std::move(seq));
        not_empty.notify_one();
        return true;
    };
    std::string line, seq;
    bool have = false;
    try {
        while (in.getline(line)) {
            if (line.empty()) continue;
            if (line[0] == '>') {
                if (have) { lengths.push_back(seq.size()); if (!push(static_cast<std::uint32_t>(names.size() - 1), std::move(seq))) break; }
                seq.clear();
                const auto sp = line.find_first_of(" \t");
                names.push_back(line.substr(1, sp == std::string::npos ? std::string::npos : sp - 1));
                have = true;
                continue;
            }
            seq += line;
        }
        if (have && lengths.size() < names.size()) { lengths.push_back(seq.size()); push(static_cast<std::uint32_t>(names.size() - 1), std::move(seq)); }
    } catch (...) {
        std::lock_guard<std::mutex> lk(m);
        if (!failure) failure = std::current_exception();
    }
    {
        std::lock_guard<std::mutex> lk(m);
        closed = true;
        not_empty.notify_all();
    }
    for (auto& t : pool) t.join();
    if (failure) std::rethrow_exception(failure);
}

struct RefHit {
    std::uint64_t code;
    std::uint32_t contig, pos;
    bool forward;
};

// Little-endian scalar writing and reading.
struct Writer {
    std::FILE* f;
    std::string path;
    void bytes(const void* p, std::size_t n) { if (n && std::fwrite(p, 1, n, f) != n) throw std::runtime_error("anchors: cannot write " + path); }
    void u64(std::uint64_t v) { bytes(&v, 8); }
    void str(const std::string& s) { u64(s.size()); bytes(s.data(), s.size()); }
};
struct Reader {
    const unsigned char* p;
    const unsigned char* end;
    std::string path;
    void need(std::size_t n) const { if (static_cast<std::size_t>(end - p) < n) throw std::runtime_error("anchors: " + path + " is truncated"); }
    std::uint64_t u64() { need(8); std::uint64_t v; std::memcpy(&v, p, 8); p += 8; return v; }
    std::string str() {
        const std::uint64_t n = u64();
        need(n);
        std::string s(reinterpret_cast<const char*>(p), n);
        p += n;
        return s;
    }
};

std::uint16_t get16(const unsigned char* p) { std::uint16_t v; std::memcpy(&v, p, 2); return v; }
std::uint32_t get32(const unsigned char* p) { std::uint32_t v; std::memcpy(&v, p, 4); return v; }

}  // namespace

bool sampled(std::uint64_t canonical_kmer) {
    return (kmer::mix64(canonical_kmer ^ kSampleSeed) >> 32) % kSampleOneIn == 0;
}

std::uint64_t ReferenceAnchors::fingerprint() const {
    std::uint64_t h = kmer::mix64(kK * 1000003 + kS * 1009 + kSampleOneIn);
    for (std::size_t i = 0; i < code.size(); ++i)
        h = kmer::mix64(h ^ code[i] ^ (static_cast<std::uint64_t>(at[i].contig) << 40) ^ (static_cast<std::uint64_t>(at[i].pos) << 1) ^ at[i].forward);
    for (const auto& c : contigs) for (const char ch : c) h = kmer::mix64(h ^ static_cast<unsigned char>(ch));
    return h;
}

ReferenceAnchors build_reference(const std::string& genome, const std::string& fasta, std::size_t threads,
                                 const std::function<void(const std::string&)>& log) {
    std::vector<RefHit> hits;
    std::mutex hm;
    ReferenceAnchors ref;
    ref.genome = genome;
    for_each_contig(fasta, threads, ref.contigs, ref.contig_lengths, [&](std::uint32_t contig, const std::string& seq) {
        if (seq.size() >= (1ULL << 32)) throw std::runtime_error("anchors: " + fasta + ": a contig of 4 Gb or more");
        std::vector<RefHit> local;
        kmer::for_each_closed_syncmer(seq.data(), seq.size(), kK, kS, [&](std::uint64_t c, std::size_t p, bool fwd) {
            if (sampled(c)) local.push_back({c, contig, static_cast<std::uint32_t>(p), fwd});
        });
        std::lock_guard<std::mutex> lk(hm);
        hits.insert(hits.end(), local.begin(), local.end());
    });
    std::sort(hits.begin(), hits.end(), [](const RefHit& a, const RefHit& b) { return a.code < b.code; });
    std::vector<RefHit> unique;
    for (std::size_t i = 0; i < hits.size();) {
        std::size_t e = i + 1;
        while (e < hits.size() && hits[e].code == hits[i].code) ++e;
        if (e - i == 1) unique.push_back(hits[i]);
        i = e;
    }
    std::sort(unique.begin(), unique.end(), [](const RefHit& a, const RefHit& b) {
        return a.contig != b.contig ? a.contig < b.contig : a.pos < b.pos;
    });
    for (const RefHit& h : unique) {
        ref.code.push_back(h.code);
        ref.at.push_back({h.contig, h.pos, h.forward});
    }
    if (log) log("anchors: unique: " + std::to_string(ref.code.size()) + "; sampled syncmers: " + std::to_string(hits.size()));
    return ref;
}

void save(const ReferenceAnchors& ref, const std::string& path) {
    const std::string tmp = path + ".tmp";
    Writer w{std::fopen(tmp.c_str(), "wb"), tmp};
    if (w.f == nullptr) throw std::runtime_error("anchors: cannot write " + tmp);
    try {
        w.u64(kRefMagic);
        w.u64(kK); w.u64(kS); w.u64(kSampleOneIn);
        w.str(ref.genome);
        w.u64(ref.contigs.size());
        for (std::size_t c = 0; c < ref.contigs.size(); ++c) { w.str(ref.contigs[c]); w.u64(ref.contig_lengths[c]); }
        w.u64(ref.code.size());
        for (std::size_t i = 0; i < ref.code.size(); ++i) {
            w.u64(ref.code[i]);
            w.u64((static_cast<std::uint64_t>(ref.at[i].contig) << 33) | (static_cast<std::uint64_t>(ref.at[i].pos) << 1) | ref.at[i].forward);
        }
        w.u64(ref.fingerprint());
    } catch (...) { std::fclose(w.f); throw; }
    if (std::fflush(w.f) != 0 || ::fsync(fileno(w.f)) != 0 || std::fclose(w.f) != 0) throw std::runtime_error("anchors: cannot write " + tmp);
    fs::rename(tmp, path);
}

ReferenceAnchors load_reference(const std::string& path) {
    const bg::Mapping m(path);
    Reader r{static_cast<const unsigned char*>(m.data()), static_cast<const unsigned char*>(m.data()) + m.size(), path};
    if (m.size() < 32 || r.u64() != kRefMagic || r.u64() != kK || r.u64() != kS || r.u64() != kSampleOneIn)
        throw std::runtime_error("anchors: " + path + " is not a reference anchor set of this version");
    ReferenceAnchors ref;
    ref.genome = r.str();
    const std::uint64_t nc = r.u64();
    for (std::uint64_t c = 0; c < nc; ++c) { ref.contigs.push_back(r.str()); ref.contig_lengths.push_back(r.u64()); }
    const std::uint64_t n = r.u64();
    r.need(n * 16);
    ref.code.resize(n);
    ref.at.resize(n);
    for (std::uint64_t i = 0; i < n; ++i) {
        ref.code[i] = r.u64();
        const std::uint64_t v = r.u64();
        ref.at[i] = {static_cast<std::uint32_t>(v >> 33), static_cast<std::uint32_t>((v >> 1) & 0xFFFFFFFFULL), (v & 1) != 0};
        if (ref.at[i].contig >= nc) throw std::runtime_error("anchors: " + path + " is corrupt");
    }
    if (r.u64() != ref.fingerprint()) throw std::runtime_error("anchors: " + path + " is corrupt (fingerprint)");
    return ref;
}

AnchorLookup::AnchorLookup(const ReferenceAnchors& ref) {
    std::uint64_t size = 1024;
    while (size < ref.code.size() * 2) size <<= 1;
    slots_.assign(size, 0);
    ids_.assign(size, 0);
    mask_ = size - 1;
    for (std::size_t i = 0; i < ref.code.size(); ++i) {
        std::uint64_t s = kmer::mix64(ref.code[i]) & mask_;
        while (slots_[s] != 0) s = (s + 1) & mask_;
        slots_[s] = (ref.code[i] << 1) | 1;
        ids_[s] = static_cast<std::uint32_t>(i);
    }
}

std::int64_t AnchorLookup::find(std::uint64_t code) const {
    const std::uint64_t key = (code << 1) | 1;
    for (std::uint64_t s = kmer::mix64(code) & mask_;; s = (s + 1) & mask_) {
        if (slots_[s] == 0) return -1;
        if (slots_[s] == key) return ids_[s];
    }
}

void scan_anchors(const AnchorLookup& lookup, const char* seq, std::size_t n, std::uint64_t offset,
                  std::uint32_t genome, std::uint32_t contig, std::vector<Hit>& out) {
    kmer::for_each_closed_syncmer(seq, n, kK, kS, [&](std::uint64_t c, std::size_t p, bool fwd) {
        if (!sampled(c)) return;
        const std::int64_t id = lookup.find(c);
        if (id >= 0) out.push_back({static_cast<std::uint32_t>(id), genome, contig, static_cast<std::uint32_t>(offset + p), fwd});
    });
}

std::vector<Hit> map_fasta(const AnchorLookup& lookup, const std::string& fasta, std::uint32_t genome, std::size_t threads,
                           std::vector<std::string>& contigs, std::vector<std::uint64_t>& lengths) {
    std::vector<Hit> hits;
    std::mutex hm;
    for_each_contig(fasta, threads, contigs, lengths, [&](std::uint32_t contig, const std::string& seq) {
        if (seq.size() >= (1ULL << 32)) throw std::runtime_error("anchors: " + fasta + ": a contig of 4 Gb or more");
        std::vector<Hit> local;
        scan_anchors(lookup, seq.data(), seq.size(), 0, genome, contig, local);
        std::lock_guard<std::mutex> lk(hm);
        hits.insert(hits.end(), local.begin(), local.end());
    });
    return hits;
}

// Segment anchor file: magic (written last), version, fingerprint, anchors, genomes, places, then per
// genome its name and contig names, then first[anchors + 1] (u64) and the places, 12 bytes each.
// The anchor file: header, genome and contig names, the per-anchor index `first` (places before each
// anchor, then the total), and the places, which `places` writes through its argument in anchor order.
void write_anchor_file(const std::string& path, const ReferenceAnchors& ref, const std::vector<std::string>& genome_names,
                       const std::vector<std::vector<std::string>>& contig_names,
                       const std::vector<std::vector<std::uint64_t>>& contig_lengths, const std::vector<std::uint64_t>& first,
                       const std::function<void(const std::function<void(const void*, std::size_t)>&)>& places) {
    const std::size_t n = ref.code.size();
    const std::string tmp = path + ".tmp";
    Writer w{std::fopen(tmp.c_str(), "wb"), tmp};
    if (w.f == nullptr) throw std::runtime_error("anchors: cannot write " + tmp);
    try {
        w.u64(0); w.u64(2); w.u64(ref.fingerprint()); w.u64(n); w.u64(genome_names.size()); w.u64(first[n]);
        std::string names;
        for (std::size_t g = 0; g < genome_names.size(); ++g) {
            auto put = [&](std::uint64_t v) { names.append(reinterpret_cast<const char*>(&v), 8); };
            put(genome_names[g].size()); names += genome_names[g];
            const auto& cn = g < contig_names.size() ? contig_names[g] : std::vector<std::string>{};
            const auto& cl = g < contig_lengths.size() ? contig_lengths[g] : std::vector<std::uint64_t>{};
            if (cl.size() != cn.size()) throw std::runtime_error("anchors: contig names and lengths differ in number");
            put(cn.size());
            for (std::size_t c = 0; c < cn.size(); ++c) { put(cn[c].size()); names += cn[c]; put(cl[c]); }
        }
        while (names.size() % 8) names += '\0';
        w.u64(names.size());
        w.bytes(names.data(), names.size());
        w.bytes(first.data(), first.size() * 8);
        places([&](const void* p, std::size_t len) { w.bytes(p, len); });
        const std::uint64_t magic = kSegMagic;
        if (std::fflush(w.f) != 0 || ::fsync(fileno(w.f)) != 0 || std::fseek(w.f, 0, SEEK_SET) != 0) throw std::runtime_error("anchors: cannot write " + tmp);
        w.u64(magic);
    } catch (...) { std::fclose(w.f); std::error_code ec; fs::remove(tmp, ec); throw; }
    if (std::fflush(w.f) != 0 || ::fsync(fileno(w.f)) != 0 || std::fclose(w.f) != 0) throw std::runtime_error("anchors: cannot write " + tmp);
    fs::rename(tmp, path);
}

void write_segment_anchors(const std::string& path, const ReferenceAnchors& ref, const std::vector<std::string>& genome_names,
                           const std::vector<std::vector<std::string>>& contig_names,
                           const std::vector<std::vector<std::uint64_t>>& contig_lengths, std::vector<Hit>& hits) {
    if (genome_names.size() > 65535) throw std::runtime_error("anchors: at most 65535 genomes in a segment");
    std::sort(hits.begin(), hits.end(), [](const Hit& a, const Hit& b) {
        if (a.anchor != b.anchor) return a.anchor < b.anchor;
        if (a.genome != b.genome) return a.genome < b.genome;
        return a.contig != b.contig ? a.contig < b.contig : a.pos < b.pos;
    });
    const std::size_t n = ref.code.size();
    std::vector<std::uint64_t> first(n + 1, 0);
    std::string places;
    std::size_t h = 0;
    for (std::size_t a = 0; a < n; ++a) {
        first[a] = places.size() / kPlaceBytes;
        while (h < hits.size() && hits[h].anchor == a) {
            std::size_t e = h;
            // capped per genome and contig: a graph genome's haplotypes are its contigs, each its own copy
            while (e < hits.size() && hits[e].anchor == a && hits[e].genome == hits[h].genome && hits[e].contig == hits[h].contig) ++e;
            const bool repeated = e - h > kMaxOccurrences;
            for (std::size_t i = h; i < (repeated ? h + 1 : e); ++i) {
                unsigned char rec[kPlaceBytes];
                const std::uint16_t g = static_cast<std::uint16_t>(hits[i].genome);
                const std::uint16_t flags = static_cast<std::uint16_t>((hits[i].forward ? kForward : 0) | (repeated ? kRepeated : 0));
                const std::uint32_t c = hits[i].contig, p = repeated ? 0 : hits[i].pos;
                std::memcpy(rec, &g, 2); std::memcpy(rec + 2, &flags, 2); std::memcpy(rec + 4, &c, 4); std::memcpy(rec + 8, &p, 4);
                places.append(reinterpret_cast<const char*>(rec), kPlaceBytes);
            }
            h = e;
        }
    }
    first[n] = places.size() / kPlaceBytes;
    write_anchor_file(path, ref, genome_names, contig_names, contig_lengths, first,
                      [&](const std::function<void(const void*, std::size_t)>& put) { put(places.data(), places.size()); });
}

SegmentAnchors::SegmentAnchors(const std::string& path) : map_(std::make_unique<bg::Mapping>(path)) {
    const auto* base = static_cast<const unsigned char*>(map_->data());
    Reader r{base, base + map_->size(), path};
    if (map_->size() < 56 || r.u64() != kSegMagic) throw std::runtime_error("anchors: " + path + " is incomplete or not an anchor file of this version");
    if (r.u64() != 2) throw std::runtime_error("anchors: " + path + ": unknown version");
    fingerprint_ = r.u64();
    anchors_ = r.u64();
    const std::uint64_t ng = r.u64(), np = r.u64(), names_bytes = r.u64();
    r.need(names_bytes);
    Reader nr{r.p, r.p + names_bytes, path};
    for (std::uint64_t g = 0; g < ng; ++g) {
        genomes_.push_back(nr.str());
        contigs_.emplace_back();
        lengths_.emplace_back();
        const std::uint64_t nc = nr.u64();
        for (std::uint64_t c = 0; c < nc; ++c) { contigs_.back().push_back(nr.str()); lengths_.back().push_back(nr.u64()); }
    }
    r.p += names_bytes;
    r.need((anchors_ + 1) * 8 + np * kPlaceBytes);
    first_ = reinterpret_cast<const std::uint64_t*>(r.p);
    places_ = r.p + (anchors_ + 1) * 8;
    if (first_[anchors_] != np) throw std::runtime_error("anchors: " + path + " is corrupt");
    for (std::size_t a = 0; a < anchors_; ++a)
        if (first_[a] > first_[a + 1]) throw std::runtime_error("anchors: " + path + " is corrupt");
}

SegmentAnchors::~SegmentAnchors() = default;

void SegmentAnchors::places(std::size_t a, std::vector<Place>& out) const {
    out.clear();
    for (std::uint64_t i = first_[a]; i < first_[a + 1]; ++i) {
        const unsigned char* rec = places_ + i * kPlaceBytes;
        const std::uint16_t flags = get16(rec + 2);
        Place p;
        p.genome = get16(rec);
        p.forward = (flags & kForward) != 0;
        p.repeated = (flags & kRepeated) != 0;
        p.contig = get32(rec + 4);
        p.pos = get32(rec + 8);
        if (p.genome >= genomes_.size() || p.contig >= contigs_[p.genome].size())
            throw std::runtime_error("anchors: a place outside its genome's contigs");
        out.push_back(p);
    }
}

std::uint64_t SegmentAnchors::bytes_touched(std::size_t a_begin, std::size_t a_end) const {
    return (a_end - a_begin + 1) * 8 + (first_[a_end] - first_[a_begin]) * kPlaceBytes;
}

void merge_segment_anchors(const std::vector<std::string>& inputs, const ReferenceAnchors& ref, const std::string& out) {
    std::vector<std::unique_ptr<SegmentAnchors>> in;
    std::vector<std::string> genomes;
    std::vector<std::vector<std::string>> contigs;
    std::vector<std::vector<std::uint64_t>> lengths;
    std::vector<std::uint32_t> genome_base;
    for (const auto& p : inputs) {
        in.push_back(std::make_unique<SegmentAnchors>(p));
        if (in.back()->fingerprint() != ref.fingerprint() || in.back()->anchors() != ref.code.size())
            throw std::runtime_error("anchors: " + p + " belongs to another anchor set");
        genome_base.push_back(static_cast<std::uint32_t>(genomes.size()));
        for (std::size_t g = 0; g < in.back()->genomes().size(); ++g) {
            genomes.push_back(in.back()->genomes()[g]);
            contigs.push_back(in.back()->contigs(g));
            lengths.emplace_back();
            for (std::size_t c = 0; c < contigs.back().size(); ++c) lengths.back().push_back(in.back()->contig_length(g, c));
        }
    }
    // Each input is anchor-major with its places capped per genome and contig, and inputs hold distinct
    // genomes: each anchor's places are copied straight through (genome ids shifted), to a temporary file,
    // so only the per-anchor index stays in memory however many genomes there are.
    const std::size_t n = ref.code.size();
    std::vector<std::uint64_t> first(n + 1, 0);
    const std::string spool = out + ".places.tmp";
    std::FILE* sp = std::fopen(spool.c_str(), "w+b");
    if (sp == nullptr) throw std::runtime_error("anchors: cannot write " + spool);
    try {
        Writer sw{sp, spool};
        std::vector<Place> ps;
        std::uint64_t total = 0;
        for (std::size_t a = 0; a < n; ++a) {
            first[a] = total;
            for (std::size_t s = 0; s < in.size(); ++s) {
                in[s]->places(a, ps);
                for (const Place& p : ps) {
                    unsigned char rec[kPlaceBytes];
                    const std::uint16_t g = static_cast<std::uint16_t>(genome_base[s] + p.genome);
                    const std::uint16_t flags = static_cast<std::uint16_t>((p.forward ? kForward : 0) | (p.repeated ? kRepeated : 0));
                    const std::uint32_t c = p.contig, pos = p.repeated ? 0 : p.pos;
                    std::memcpy(rec, &g, 2); std::memcpy(rec + 2, &flags, 2); std::memcpy(rec + 4, &c, 4); std::memcpy(rec + 8, &pos, 4);
                    sw.bytes(rec, kPlaceBytes);
                    ++total;
                }
            }
        }
        first[n] = total;
        if (genomes.size() > 65535) throw std::runtime_error("anchors: at most 65535 genomes in a segment");
        if (std::fflush(sp) != 0 || std::fseek(sp, 0, SEEK_SET) != 0) throw std::runtime_error("anchors: cannot read " + spool);
        write_anchor_file(out, ref, genomes, contigs, lengths, first, [&](const std::function<void(const void*, std::size_t)>& put) {
            std::vector<char> buf(1 << 20);
            for (std::size_t got; (got = std::fread(buf.data(), 1, buf.size(), sp)) > 0;) put(buf.data(), got);
            if (std::ferror(sp)) throw std::runtime_error("anchors: cannot read " + spool);
        });
    } catch (...) {
        std::fclose(sp);
        std::error_code ec;
        fs::remove(spool, ec);
        throw;
    }
    std::fclose(sp);
    std::error_code ec;
    fs::remove(spool, ec);
}

std::vector<Span> locate(const ReferenceAnchors& ref, const std::vector<const SegmentAnchors*>& segments, const std::string& contig,
                         std::uint64_t start, std::uint64_t end, const LocateOptions& o, std::uint64_t* bytes_read) {
    const auto cit = std::find(ref.contigs.begin(), ref.contigs.end(), contig);
    if (cit == ref.contigs.end()) throw std::runtime_error("anchors: the reference has no contig " + contig);
    const auto c = static_cast<std::uint32_t>(cit - ref.contigs.begin());
    if (start >= end) throw std::runtime_error("anchors: the interval's start must be below its end");
    if (end > ref.contig_lengths[c])
        throw std::runtime_error("anchors: " + contig + " is " + std::to_string(ref.contig_lengths[c]) + " bp; the interval ends at " + std::to_string(end));
    for (const auto* s : segments)
        if (s->fingerprint() != ref.fingerprint()) throw std::runtime_error("anchors: a segment's anchors belong to another anchor set");
    auto first_at_or_after = [&](std::uint64_t pos) {
        return static_cast<std::size_t>(std::lower_bound(ref.at.begin(), ref.at.end(), std::make_pair(c, pos),
            [](const RefPlace& a, const std::pair<std::uint32_t, std::uint64_t>& b) {
                return a.contig != b.first ? a.contig < b.first : a.pos < b.second;
            }) - ref.at.begin());
    };
    const std::size_t c_begin = first_at_or_after(0), c_end = first_at_or_after(~0ULL >> 1);
    const std::size_t in_begin = first_at_or_after(start);
    std::size_t in_end = in_begin;
    while (in_end < c_end && ref.at[in_end].pos + kK <= end) ++in_end;
    const std::size_t w_begin = in_begin >= c_begin + o.flank ? in_begin - o.flank : c_begin;
    const std::size_t w_end = std::min(c_end, in_end + o.flank);
    const std::size_t n_inside = in_end - in_begin;
    const std::size_t required = std::min(o.min_inside, n_inside);
    constexpr std::uint64_t kExtendSlack = 10000;   // a flank anchor bounds the copy only if it sits about where the reference puts it

    struct Point { std::uint64_t ref_pos, pos; bool inside, before; };
    std::vector<Span> out;
    std::vector<Place> ps;
    for (std::size_t s = 0; s < segments.size(); ++s) {
        const SegmentAnchors& sa = *segments[s];
        if (bytes_read != nullptr && w_end > w_begin) *bytes_read += sa.bytes_touched(w_begin, w_end);
        // per genome, per (contig, same strand): the points
        std::vector<std::map<std::pair<std::uint32_t, bool>, std::vector<Point>>> groups(sa.genomes().size());
        for (std::size_t a = w_begin; a < w_end; ++a) {
            sa.places(a, ps);
            for (const Place& p : ps) {
                if (p.repeated) continue;
                const bool same = p.forward == ref.at[a].forward;
                groups[p.genome][{p.contig, same}].push_back({ref.at[a].pos, p.pos, a >= in_begin && a < in_end, a < in_begin});
            }
        }
        for (std::size_t g = 0; g < groups.size(); ++g) {
            std::vector<Span> spans;
            for (auto& group : groups[g]) {
                const auto& key = group.first;
                auto& pts = group.second;
                const bool same = key.second;
                std::sort(pts.begin(), pts.end(), [](const Point& a, const Point& b) { return a.ref_pos != b.ref_pos ? a.ref_pos < b.ref_pos : a.pos < b.pos; });
                const std::size_t n = pts.size();
                std::vector<char> used(n, 0);
                std::size_t inside_left = 0;
                for (const Point& p : pts) inside_left += p.inside;
                while (true) {
                    if (required > 0 && inside_left < required) break;
                    // best-scoring chain over unused points: reference order, genome order by strand, gaps bounded
                    std::vector<double> score(n, 0);
                    std::vector<std::size_t> prev(n, SIZE_MAX);
                    std::size_t best = SIZE_MAX;
                    for (std::size_t j = 0; j < n; ++j) {
                        if (used[j]) continue;
                        score[j] = 1;
                        for (std::size_t i = 0; i < j; ++i) {
                            if (used[i] || pts[i].ref_pos >= pts[j].ref_pos) continue;
                            if (same ? pts[j].pos <= pts[i].pos : pts[j].pos >= pts[i].pos) continue;
                            const std::uint64_t dr = pts[j].ref_pos - pts[i].ref_pos;
                            const std::uint64_t dg = same ? pts[j].pos - pts[i].pos : pts[i].pos - pts[j].pos;
                            if (dg > dr + dr / 2 + o.slack || dr > dg + dg / 2 + o.slack) continue;
                            const double dd = static_cast<double>(dg > dr ? dg - dr : dr - dg);
                            const double sc = score[i] + 1 - o.gap_cost_per_kb * dd / 1000.0;
                            if (sc > score[j]) { score[j] = sc; prev[j] = i; }
                        }
                        if (best == SIZE_MAX || score[j] > score[best]) best = j;
                    }
                    if (best == SIZE_MAX) break;
                    std::vector<std::size_t> chain;
                    for (std::size_t j = best; j != SIZE_MAX; j = prev[j]) chain.push_back(j);
                    std::reverse(chain.begin(), chain.end());   // reference order
                    std::size_t inside = 0, flanks = 0;
                    std::size_t first_in = SIZE_MAX, last_in = SIZE_MAX, before = SIZE_MAX, after = SIZE_MAX;
                    for (const std::size_t j : chain) {
                        used[j] = 1;
                        if (pts[j].inside) { ++inside; --inside_left; if (first_in == SIZE_MAX) first_in = j; last_in = j; }
                        else { ++flanks; if (pts[j].before) before = j; else if (after == SIZE_MAX) after = j; }
                    }
                    if (inside < required || (required == 0 && (before == SIZE_MAX || after == SIZE_MAX))) continue;
                    // the copy on the genome, in reference terms: left = reference start side
                    std::uint64_t lo, hi, ilo, ihi;
                    auto gpos = [&](std::size_t j) { return pts[j].pos; };
                    auto extend = [&](std::size_t inner, std::size_t flank_j, std::uint64_t ref_boundary, bool left) {
                        // genome position bounding the copy on the reference's `left` side
                        const std::uint64_t ref_gap = left ? pts[inner].ref_pos - pts[flank_j == SIZE_MAX ? inner : flank_j].ref_pos
                                                           : pts[flank_j == SIZE_MAX ? inner : flank_j].ref_pos - pts[inner].ref_pos;
                        if (flank_j != SIZE_MAX) {
                            const std::uint64_t g_gap = gpos(inner) > gpos(flank_j) ? gpos(inner) - gpos(flank_j) : gpos(flank_j) - gpos(inner);
                            if (g_gap <= ref_gap + ref_gap / 2 + kExtendSlack) return gpos(flank_j);
                        }
                        // no trustworthy flank: extend by the reference's distance to the interval boundary
                        const std::uint64_t d = left ? (pts[inner].ref_pos > ref_boundary ? pts[inner].ref_pos - ref_boundary : 0)
                                                     : (ref_boundary > pts[inner].ref_pos + kK ? ref_boundary - pts[inner].ref_pos - kK : 0);
                        const bool towards_lower = left == same;   // on the genome, the reference's left is lower when same strand
                        if (towards_lower) return gpos(inner) > d ? gpos(inner) - d : 0;
                        return gpos(inner) + d;
                    };
                    if (first_in != SIZE_MAX) {
                        const std::uint64_t a = extend(first_in, before, start, true);
                        const std::uint64_t b = extend(last_in, after, end, false);
                        lo = std::min({a, b, gpos(first_in), gpos(last_in)});
                        hi = std::max({a, b, gpos(first_in), gpos(last_in)});
                        ilo = std::min(gpos(first_in), gpos(last_in));
                        ihi = std::max(gpos(first_in), gpos(last_in));
                    } else {   // an interval without inside anchors: between the nearest flanks
                        lo = std::min(gpos(before), gpos(after));
                        hi = std::max(gpos(before), gpos(after));
                        ilo = lo;
                        ihi = hi;
                    }
                    Span sp;
                    sp.genome = sa.genomes()[g];
                    sp.contig = sa.contigs(g)[key.first];
                    sp.start = lo;
                    sp.end = std::min(hi + kK, sa.contig_length(g, key.first));   // an estimated end stops at the contig's
                    sp.forward = same;
                    sp.anchors_inside = inside;
                    sp.anchors_flank = flanks;
                    sp.inner_start = ilo;
                    sp.inner_end = std::min(ihi + kK, sp.end);
                    spans.push_back(sp);
                }
            }
            // copies that overlap on a contig are one stretch of sequence: join them
            std::sort(spans.begin(), spans.end(), [](const Span& a, const Span& b) { return a.contig != b.contig ? a.contig < b.contig : a.start < b.start; });
            std::vector<Span> joined;
            for (const Span& sp : spans) {
                if (!joined.empty() && joined.back().contig == sp.contig && sp.start < joined.back().end) {
                    Span& j = joined.back();
                    if (sp.anchors_inside > j.anchors_inside) j.forward = sp.forward;
                    j.end = std::max(j.end, sp.end);
                    j.inner_start = std::min(j.inner_start, sp.inner_start);
                    j.inner_end = std::max(j.inner_end, sp.inner_end);
                    j.anchors_inside += sp.anchors_inside;
                    j.anchors_flank += sp.anchors_flank;
                    j.chains += sp.chains;
                } else {
                    joined.push_back(sp);
                }
            }
            std::stable_sort(joined.begin(), joined.end(), [](const Span& a, const Span& b) { return a.anchors_inside > b.anchors_inside; });
            out.insert(out.end(), joined.begin(), joined.end());
        }
    }
    return out;
}

}  // namespace panvar::anchors
