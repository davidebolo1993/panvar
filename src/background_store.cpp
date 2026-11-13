// The store's manifest, its lock, job claims, and reading a store as the sum of its listed segments.

#include "panvar/background.hpp"
#include "panvar/json_lite.hpp"

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

namespace panvar::bg {

namespace fs = std::filesystem;

namespace {

std::string hex64(std::uint64_t v) {
    char buf[17];
    std::snprintf(buf, sizeof buf, "%016llx", static_cast<unsigned long long>(v));
    return buf;
}

void fsync_path(const std::string& path) {
    const int fd = ::open(path.c_str(), O_RDONLY);
    if (fd >= 0) { ::fsync(fd); ::close(fd); }
}

int lock_fd(const std::string& path, bool wait) {
    const int fd = ::open(path.c_str(), O_RDWR | O_CREAT, 0644);
    if (fd < 0) throw std::runtime_error("background: cannot open " + path + ": " + std::strerror(errno));
    struct flock fl {};
    fl.l_type = F_WRLCK;
    fl.l_whence = SEEK_SET;
    while (::fcntl(fd, wait ? F_SETLKW : F_SETLK, &fl) != 0) {
        if (errno == EINTR) continue;
        const int e = errno;
        ::close(fd);
        if (!wait && (e == EACCES || e == EAGAIN)) return -1;
        throw std::runtime_error("background: cannot lock " + path + ": " + std::strerror(e));
    }
    return fd;
}

}  // namespace

StoreLock::StoreLock(const std::string& store_dir) : fd_(lock_fd((fs::path(store_dir) / ".lock").string(), true)) {}
StoreLock::~StoreLock() { if (fd_ >= 0) ::close(fd_); }

JobLock::JobLock(const std::string& dir) {
    fs::create_directories(dir);
    fd_ = lock_fd((fs::path(dir) / "job.lock").string(), false);
    if (fd_ < 0) throw std::runtime_error("background: " + dir + " is in use by another job");
}
JobLock::~JobLock() { if (fd_ >= 0) ::close(fd_); }

bool job_alive(const std::string& dir) {
    const std::string path = (fs::path(dir) / "job.lock").string();
    if (!fs::exists(path)) return false;
    const int fd = lock_fd(path, false);
    if (fd < 0) return true;
    ::close(fd);
    return false;
}

bool store_exists(const std::string& store_dir) { return fs::exists(fs::path(store_dir) / "manifest.json"); }

Manifest read_manifest(const std::string& store_dir) {
    const std::string path = (fs::path(store_dir) / "manifest.json").string();
    std::ifstream f(path);
    if (!f) throw std::runtime_error(store_dir + " is not a background store (no manifest.json)");
    std::stringstream ss;
    ss << f.rdbuf();
    json::Value v;
    try {
        v = json::parse(ss.str());
    } catch (const std::exception& e) {
        throw std::runtime_error("background: " + path + ": " + e.what());
    }
    if (!v.has("format") || v.at("format").str() != "panvar-background")
        throw std::runtime_error("background: " + path + " is not a background manifest");
    if (v.at("version").u64() != 2)
        throw std::runtime_error("background: " + path + " is format version " + v.at("version").text + "; this panvar reads version 2 (rebuild the store)");
    if (v.at("k").u64() != kK || v.at("minimizer").u64() != kM) throw std::runtime_error("background: " + path + ": k or minimizer differ");
    Manifest m;
    m.store_id = v.at("store_id").str();
    m.generation = v.at("generation").u64();
    m.partitions = static_cast<std::uint32_t>(v.at("partitions").u64());
    m.next_segment = v.at("next_segment").u64();
    for (const auto& r : v.at("coordinate_references").items)
        m.references.emplace_back(r.at("name").str(), std::stoull(r.at("fingerprint").str(), nullptr, 16));
    for (const auto& s : v.at("segments").items) {
        ManifestSegment ms;
        ms.dir = s.at("dir").str();
        ms.distinct = s.at("distinct").u64();
        ms.occurrences = s.at("occurrences").u64();
        for (const auto& g : s.at("genomes").items)
            ms.genomes.push_back({g.at("name").str(), g.at("path").str(), g.at("bytes").u64(), g.at("checksum").str()});
        m.segments.push_back(std::move(ms));
    }
    return m;
}

void write_manifest(const std::string& store_dir, const Manifest& m) {
    const std::string path = (fs::path(store_dir) / "manifest.json").string(), tmp = path + ".tmp";
    std::ostringstream o;
    o << "{\n  \"format\": \"panvar-background\",\n  \"version\": 2,\n  \"store_id\": " << json::quote(m.store_id)
      << ",\n  \"generation\": " << m.generation << ",\n  \"k\": " << kK << ",\n  \"minimizer\": " << kM
      << ",\n  \"partitions\": " << m.partitions << ",\n  \"next_segment\": " << m.next_segment << ",\n  \"coordinate_references\": [";
    for (std::size_t i = 0; i < m.references.size(); ++i)
        o << (i ? ", " : "") << "{\"name\": " << json::quote(m.references[i].first) << ", \"fingerprint\": \"" << hex64(m.references[i].second) << "\"}";
    o << "],\n  \"segments\": [";
    for (std::size_t s = 0; s < m.segments.size(); ++s) {
        const auto& ms = m.segments[s];
        o << (s ? "," : "") << "\n    {\"dir\": " << json::quote(ms.dir) << ", \"distinct\": " << ms.distinct
          << ", \"occurrences\": " << ms.occurrences << ", \"genomes\": [";
        for (std::size_t g = 0; g < ms.genomes.size(); ++g) {
            const auto& ge = ms.genomes[g];
            o << (g ? "," : "") << "\n      {\"name\": " << json::quote(ge.name) << ", \"path\": " << json::quote(ge.path)
              << ", \"bytes\": " << ge.bytes << ", \"checksum\": " << json::quote(ge.checksum) << "}";
        }
        o << "]}";
    }
    o << "\n  ]\n}\n";
    write_file_durably(tmp, o.str());
    fs::rename(tmp, path);
    fsync_path(store_dir);
}

void write_file_durably(const std::string& path, const std::string& content) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (f == nullptr) throw std::runtime_error("background: cannot write " + path + ": " + std::strerror(errno));
    const bool ok = std::fwrite(content.data(), 1, content.size(), f) == content.size() && std::fflush(f) == 0 && ::fsync(fileno(f)) == 0;
    const int e = errno;
    if (std::fclose(f) != 0 || !ok) {
        std::error_code ec;
        fs::remove(path, ec);
        throw std::runtime_error("background: cannot write " + path + ": " + std::strerror(ok ? errno : e));
    }
}

std::uint64_t content_id(const Manifest& m) {
    std::vector<std::string> g;
    for (const auto& s : m.segments) for (const auto& ge : s.genomes) g.push_back(ge.name + "\t" + ge.checksum);
    std::sort(g.begin(), g.end());
    std::uint64_t h = 0x84222325cbf29ce4ULL ^ m.partitions;
    for (const auto& x : g) for (const char c : x) h = kmer::mix64(h ^ static_cast<unsigned char>(c));
    return h;
}

Store::Store(const std::string& store_dir) : manifest_(read_manifest(store_dir)) {
    partitions_ = manifest_.partitions;
    for (const auto& ms : manifest_.segments) {
        segs_.push_back(std::make_unique<Segment>((fs::path(store_dir) / ms.dir).string()));
        if (segs_.back()->partitions() != partitions_)
            throw std::runtime_error("background: " + ms.dir + " has a partition count other than the store's");
    }
}

std::uint64_t Store::count(std::uint64_t canonical_kmer) const {
    const std::uint32_t p = partition_of(canonical_kmer, partitions_);
    std::uint64_t n = 0;
    for (const auto& s : segs_) n += s->count(canonical_kmer, p);
    return n;
}

std::vector<std::uint64_t> Store::count_batch(const std::vector<std::uint64_t>& kmers) const {
    std::vector<std::pair<std::uint32_t, std::size_t>> order(kmers.size());
    for (std::size_t i = 0; i < kmers.size(); ++i) order[i] = {partition_of(kmers[i], partitions_), i};
    std::sort(order.begin(), order.end(), [&](const auto& a, const auto& b) {
        return a.first != b.first ? a.first < b.first : kmers[a.second] < kmers[b.second];
    });
    std::vector<std::uint64_t> out(kmers.size(), 0);
    for (const auto& s : segs_)
        for (std::size_t i = 0; i < order.size();) {
            const std::uint32_t p = order[i].first;
            const Segment::Part pt = s->part(p);
            for (; i < order.size() && order[i].first == p; ++i) {
                const std::int64_t at = pt.codes.find(kmers[order[i].second]);
                if (at >= 0) out[order[i].second] += pt.counts.get(static_cast<std::uint64_t>(at));
            }
        }
    return out;
}

void Store::for_each(std::uint32_t partition, const std::function<void(std::uint64_t, std::uint64_t)>& f) const {
    std::vector<const Segment*> raw;
    for (const auto& s : segs_) raw.push_back(s.get());
    if (raw.empty()) return;
    merge_partition(raw, partition, f);
}

std::vector<std::uint64_t> Background::count_batch(const std::vector<std::uint64_t>& kmers) const {
    std::vector<std::uint64_t> out(kmers.size());
    for (std::size_t i = 0; i < kmers.size(); ++i) out[i] = count(kmers[i]);
    return out;
}

}  // namespace panvar::bg
