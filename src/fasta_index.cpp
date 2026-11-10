#include "panvar/fasta_index.hpp"

#include "panvar/gz_reader.hpp"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include <zlib.h>

namespace panvar::fasta {

namespace {

bool gzip_magic(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    unsigned char m[4] = {0, 0, 0, 0};
    in.read(reinterpret_cast<char*>(m), 4);
    return m[0] == 0x1f && m[1] == 0x8b;
}

bool bgzf_magic(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    unsigned char h[16] = {};
    in.read(reinterpret_cast<char*>(h), 16);
    return in.gcount() == 16 && h[0] == 0x1f && h[1] == 0x8b && (h[3] & 4) && h[12] == 'B' && h[13] == 'C';
}

bool exists(const std::string& p) { std::ifstream in(p); return static_cast<bool>(in); }

}  // namespace

bool IndexedFasta::indexed(const std::string& path) {
    if (!exists(path + ".fai")) return false;
    if (!gzip_magic(path)) return true;
    return bgzf_magic(path) && exists(path + ".gzi");
}

IndexedFasta::IndexedFasta(const std::string& path) : path_(path) {
    std::ifstream fai(path + ".fai");
    if (!fai) throw std::runtime_error("fasta: " + path + " has no .fai index");
    std::string line;
    while (std::getline(fai, line)) {
        std::istringstream ss(line);
        std::string name;
        Entry e{};
        if (!(ss >> name >> e.length >> e.offset >> e.line_bases >> e.line_width)) continue;
        fai_[name] = e;
    }
    if (gzip_magic(path)) {
        if (!bgzf_magic(path)) throw std::runtime_error("fasta: " + path + " is gzip but not bgzip; it cannot be read by interval");
        std::ifstream gzi(path + ".gzi", std::ios::binary);
        if (!gzi) throw std::runtime_error("fasta: " + path + " has no .gzi index");
        std::uint64_t n = 0;
        gzi.read(reinterpret_cast<char*>(&n), 8);
        gzi_.emplace_back(0, 0);
        for (std::uint64_t i = 0; i < n; ++i) {
            std::uint64_t c = 0, u = 0;
            gzi.read(reinterpret_cast<char*>(&c), 8);
            gzi.read(reinterpret_cast<char*>(&u), 8);
            if (!gzi) throw std::runtime_error("fasta: " + path + ".gzi is truncated");
            gzi_.emplace_back(c, u);
        }
        bgzf_ = true;
    }
    f_ = std::fopen(path.c_str(), "rb");
    if (f_ == nullptr) throw std::runtime_error("fasta: cannot read " + path);
}

IndexedFasta::~IndexedFasta() { if (f_ != nullptr) std::fclose(f_); }

std::uint64_t IndexedFasta::length(const std::string& contig) const {
    const auto it = fai_.find(contig);
    if (it == fai_.end()) throw std::runtime_error("fasta: " + path_ + " has no contig " + contig);
    return it->second.length;
}

void IndexedFasta::read_bytes(std::uint64_t offset, std::uint64_t n, std::string& out) {
    out.clear();
    if (n == 0) return;
    if (!bgzf_) {
        out.resize(n);
        if (std::fseek(f_, static_cast<long>(offset), SEEK_SET) != 0 || std::fread(&out[0], 1, n, f_) != n)
            throw std::runtime_error("fasta: cannot read " + path_);
        return;
    }
    // the BGZF block holding `offset`: the last .gzi entry at or before it, then blocks in turn
    auto it = std::upper_bound(gzi_.begin(), gzi_.end(), offset, [](std::uint64_t v, const auto& e) { return v < e.second; });
    --it;
    std::uint64_t coff = it->first, uoff = it->second;
    std::vector<unsigned char> comp;
    while (out.size() < n) {
        if (coff != cached_coff_) {
            unsigned char h[18];
            if (std::fseek(f_, static_cast<long>(coff), SEEK_SET) != 0 || std::fread(h, 1, 18, f_) != 18 || h[0] != 0x1f || h[1] != 0x8b ||
                h[12] != 'B' || h[13] != 'C')
                throw std::runtime_error("fasta: " + path_ + ": not a BGZF block at " + std::to_string(coff));
            const std::uint64_t xlen = h[10] | (h[11] << 8);
            const std::uint64_t bsize = (h[16] | (h[17] << 8)) + 1ULL;
            comp.resize(bsize - 18);
            if (std::fread(comp.data(), 1, comp.size(), f_) != comp.size()) throw std::runtime_error("fasta: " + path_ + " is truncated");
            const std::uint64_t cdata = bsize - xlen - 20;   // deflate data between the header and CRC32 + ISIZE
            const unsigned char* data = comp.data() + (xlen - 6);
            std::uint32_t isize = 0;
            std::memcpy(&isize, comp.data() + comp.size() - 4, 4);
            cached_.assign(isize, '\0');
            if (isize > 0) {
                z_stream zs{};
                if (inflateInit2(&zs, -15) != Z_OK) throw std::runtime_error("fasta: zlib");
                zs.next_in = const_cast<unsigned char*>(data);
                zs.avail_in = static_cast<uInt>(cdata);
                zs.next_out = reinterpret_cast<unsigned char*>(&cached_[0]);
                zs.avail_out = isize;
                const int rc = inflate(&zs, Z_FINISH);
                inflateEnd(&zs);
                if (rc != Z_STREAM_END) throw std::runtime_error("fasta: " + path_ + ": corrupt BGZF block");
            }
            cached_coff_ = coff;
            cached_uoff_ = uoff;
            cached_next_ = coff + bsize;
        }
        if (cached_.empty()) throw std::runtime_error("fasta: " + path_ + ": read past the end");
        const std::uint64_t block_end = cached_uoff_ + cached_.size();
        const std::uint64_t want = offset + out.size();
        if (want < block_end) {
            const std::uint64_t from = want - cached_uoff_;
            const std::uint64_t take = std::min<std::uint64_t>(n - out.size(), cached_.size() - from);
            out.append(cached_, from, take);
        }
        uoff = block_end;
        coff = cached_next_;
    }
}

std::string IndexedFasta::fetch(const std::string& contig, std::uint64_t start, std::uint64_t end) {
    const auto it = fai_.find(contig);
    if (it == fai_.end()) throw std::runtime_error("fasta: " + path_ + " has no contig " + contig);
    const Entry& e = it->second;
    end = std::min(end, e.length);
    if (start >= end) return {};
    auto file_pos = [&](std::uint64_t p) { return e.offset + (p / e.line_bases) * e.line_width + p % e.line_bases; };
    const std::uint64_t a = file_pos(start), b = file_pos(end - 1) + 1;
    std::string raw;
    read_bytes(a, b - a, raw);
    std::string out;
    out.reserve(end - start);
    for (const char c : raw) if (c != '\n' && c != '\r') out += c;
    if (out.size() != end - start) throw std::runtime_error("fasta: " + path_ + ": the .fai does not match the file");
    return out;
}

std::vector<std::string> fetch_intervals(const std::string& path, const std::vector<Interval>& iv) {
    std::vector<std::string> out(iv.size());
    if (iv.empty()) return out;
    if (IndexedFasta::indexed(path)) {
        IndexedFasta fa(path);
        for (std::size_t i = 0; i < iv.size(); ++i) out[i] = fa.fetch(iv[i].contig, iv[i].start, iv[i].end);
        return out;
    }
    // one pass: each needed contig held while it is read
    std::map<std::string, std::vector<std::size_t>> want;
    for (std::size_t i = 0; i < iv.size(); ++i) want[iv[i].contig].push_back(i);
    GzLineReader in(path);
    if (!in.ok()) throw std::runtime_error("fasta: cannot read " + path);
    std::string line, name, seq;
    bool keep = false;
    auto finish = [&] {
        if (!keep) return;
        for (const std::size_t i : want[name]) out[i] = iv[i].start < seq.size() ? seq.substr(iv[i].start, std::min<std::uint64_t>(iv[i].end, seq.size()) - iv[i].start) : std::string();
    };
    while (in.getline(line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            finish();
            const auto sp = line.find_first_of(" \t");
            name = line.substr(1, sp == std::string::npos ? std::string::npos : sp - 1);
            keep = want.count(name) != 0;
            seq.clear();
            continue;
        }
        if (keep) seq += line;
    }
    finish();
    return out;
}

}  // namespace panvar::fasta
