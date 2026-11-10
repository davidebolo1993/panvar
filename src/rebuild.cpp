#include "panvar/rebuild.hpp"
#include "panvar/module_fold.hpp"
#include "panvar/minimap2_align.hpp"

#include <algorithm>
#include <unordered_map>
#include <optional>
#include <iomanip>
#include <cstdlib>
#include <array>
#include <cctype>
#include <atomic>
#include <csignal>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <random>
#include <thread>
#include <unordered_set>
#include <utility>
#include <vector>

#include <zlib.h>
#if !defined(_WIN32)
#include <unistd.h>
#endif

#include "panvar/align.hpp"
#include "panvar/cli_utils.hpp"
#include "panvar/gfa.hpp"
#include "panvar/graph_utils.hpp"
#include "panvar/parallel.hpp"

#include <ksw2.h>

#include "gfa.h"
#include "minigraph.h"
#include "mgpriv.h" // mg_gchain_free: declared here rather than in the public header

namespace panvar {
namespace {

// Wall-clock HH:MM:SS, matching cli::RunLog's line prefix so rebuild's progress lines read the same as
// its RunLog summary lines and as every other module.
std::string hms() {
    const std::time_t t = std::time(nullptr);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    char buf[16];
    std::strftime(buf, sizeof(buf), "%H:%M:%S", &tm);
    return std::string(buf);
}

std::string spell(const Graph& g, const PathRecord& path) {
    std::string out;
    for (const PathStep& step : path.steps) {
        const auto it = g.nodes.find(step.node_id);
        // validate_graph_paths() has already refused a graph where this could miss; the guard means a
        // future caller cannot silently spell a short sequence.
        if (it == g.nodes.end())
            throw std::runtime_error("rebuild: path " + path.name + " references missing node " + step.node_id);
        out += step.reverse ? reverse_complement(it->second.sequence) : it->second.sequence;
    }
    return out;
}

// Degree per HANDLE, not per node. A bidirected graph node has two ends, and pooling them conflates two
// different shapes: a node with 25 neighbours on each side is a clean two-sided branch, while one with
// 50 on a single side is the tangle the gate exists to find -- yet pooled they both read 50. The node's
// degree is therefore the larger of its two handle degrees. Self-loops are counted separately: they are
// how a tandem array appears after folding, and they are not pathology.
void degree_stats(const Graph& g, std::size_t hub_degree, std::size_t& hubs, std::size_t& maxdeg,
                  std::size_t* selfloops) {
    hubs = 0;
    maxdeg = 0;
    std::size_t loops = 0;
    for (const auto& kv : g.nodes) {
        std::unordered_set<std::string> left, right;
        bool loop = false;
        for (const Neighbor& n : kv.second.start) {
            if (n.node_id == kv.first) { loop = true; continue; }
            left.insert(n.node_id);
        }
        for (const Neighbor& n : kv.second.end) {
            if (n.node_id == kv.first) { loop = true; continue; }
            right.insert(n.node_id);
        }
        if (loop) ++loops;
        const std::size_t deg = std::max(left.size(), right.size());
        if (deg > maxdeg) maxdeg = deg;
        if (deg >= hub_degree) ++hubs;
    }
    if (selfloops != nullptr) *selfloops = loops;
}

struct KmerStats {
    std::size_t distinct = 0; // how much DIFFERENT sequence the haplotype carries
    std::size_t total = 0;    // k-mer positions; total - distinct is repeat-copy redundancy
};

// Rolling 2-bit encoding: substr-per-position costs tens of millions of allocations at cohort scale.
// k <= 31 fits a uint64; above that fall back to substr.
//
// Three properties the richness metric has to have:
//   CANONICAL   -- a k-mer and its reverse complement count as one, or a haplotype's richness depends
//                  on which strand the GFA stores it on, and the seed changes when an input is flipped.
//   CONSISTENT  -- `total` counts only the windows `distinct` counts, so the redundancy figure
//                  (total - distinct) does not shift with how much ambiguity a haplotype carries.
//   UNIFORM     -- k > 31 behaves the same way, so crossing that boundary does not silently change
//                  what the ranking means.
KmerStats kmer_stats(const std::string& s, std::size_t k) {
    KmerStats st;
    if (k == 0 || s.size() < k) return st;
    if (k > 31) {
        std::unordered_set<std::string> set;
        set.reserve(s.size());
        for (std::size_t i = 0; i + k <= s.size(); ++i) {
            std::string up = s.substr(i, k);
            bool ok = true;
            for (char& c : up) {
                c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
                if (c != 'A' && c != 'C' && c != 'G' && c != 'T') { ok = false; break; }
            }
            if (!ok) continue;                        // ambiguity: not a countable window
            const std::string rc = reverse_complement(up);
            set.insert(up < rc ? up : rc);            // canonical
            ++st.total;
        }
        st.distinct = set.size();
        return st;
    }
    std::unordered_set<std::uint64_t> set;
    set.reserve(s.size());
    const std::uint64_t mask = (k < 32) ? ((1ULL << (2 * k)) - 1) : ~0ULL;
    const unsigned shift = static_cast<unsigned>(2 * (k - 1));
    std::uint64_t fwd = 0, rev = 0;
    std::size_t valid = 0;
    for (std::size_t i = 0; i < s.size(); ++i) {
        int b;
        switch (s[i]) {
            case 'A': case 'a': b = 0; break;
            case 'C': case 'c': b = 1; break;
            case 'G': case 'g': b = 2; break;
            case 'T': case 't': b = 3; break;
            default: valid = 0; fwd = 0; rev = 0; continue;  // ambiguity resets the roll
        }
        fwd = ((fwd << 2) | static_cast<std::uint64_t>(b)) & mask;
        rev = (rev >> 2) | (static_cast<std::uint64_t>(3 - b) << shift);
        if (++valid >= k) { set.insert(std::min(fwd, rev)); ++st.total; }
    }
    st.distinct = set.size();
    return st;
}

// A unique sibling of `final` -- outputs are staged here and renamed into place, so an interrupted or
// failed run never leaves a half-written graph where a complete one is expected.
std::filesystem::path staging_path(const std::string& final_path) {
    static std::atomic<unsigned> seq{0};
    std::random_device rd;
    const std::string tag = std::to_string(rd()) + "." + std::to_string(seq++);
    return std::filesystem::path(final_path + ".rebuild-tmp." + tag);
}

void commit_staged(const std::filesystem::path& staged, const std::string& final_path) {
    std::error_code ec;
    std::filesystem::rename(staged, final_path, ec);
    if (ec) {
        // Across filesystems rename fails; fall back to copy-then-remove, still not in place.
        std::filesystem::copy_file(staged, final_path,
                                   std::filesystem::copy_options::overwrite_existing, ec);
        std::filesystem::remove(staged);
        if (ec) throw std::runtime_error("rebuild: cannot move output into place: " + final_path);
    }
}

// Removes a staged file on scope exit unless it was committed, so an exception between staging and
// commit does not leave a complete rebuilt graph beside the output.
class StagedFile {
public:
    explicit StagedFile(std::filesystem::path p) : p_(std::move(p)) {}
    ~StagedFile() {
        if (p_.empty()) return;
        std::error_code ec;
        std::filesystem::remove(p_, ec);
    }
    StagedFile(const StagedFile&) = delete;
    StagedFile& operator=(const StagedFile&) = delete;
    void release() { p_.clear(); }

private:
    std::filesystem::path p_;
};

// The output's codec follows the name asked for, not the input's: gzip only when it ends in .gz.
bool wants_gzip(const std::string& path) {
    return path.size() >= 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
}

bool is_gzip_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    unsigned char m[2] = {0, 0};
    in.read(reinterpret_cast<char*>(m), 2);
    return in.gcount() == 2 && m[0] == 0x1f && m[1] == 0x8b;
}

void write_text(const std::filesystem::path& path, const std::string& text, bool gz) {
    if (gz) {
        gzFile f = gzopen(path.string().c_str(), "wb6");
        if (f == nullptr) throw std::runtime_error("rebuild: cannot write " + path.string());
        // gzwrite takes an unsigned length, so large outputs are written in chunks.
        bool ok = true;
        std::size_t chunk = std::size_t{1} << 30;
        const char* test_chunk = std::getenv("PANVAR_TEST_GZ_CHUNK");   // tests: a small chunk, and the count
        if (test_chunk != nullptr) chunk = std::max<std::size_t>(1, std::strtoull(test_chunk, nullptr, 10));
        std::size_t writes = 0;
        for (std::size_t at = 0; ok && at < text.size(); at += chunk, ++writes) {
            const std::size_t n = std::min(chunk, text.size() - at);
            ok = gzwrite(f, text.data() + at, static_cast<unsigned int>(n)) == static_cast<int>(n);
        }
        if (test_chunk != nullptr) std::cerr << "rebuild: gzip output written in " << writes << " chunk(s)\n";
        if (gzclose(f) != Z_OK || !ok) throw std::runtime_error("rebuild: failed writing " + path.string());
        return;
    }
    std::ofstream o(path, std::ios::binary);
    if (!o) throw std::runtime_error("rebuild: cannot write " + path.string());
    o << text;
    o.flush();
    if (!o) throw std::runtime_error("rebuild: failed writing " + path.string());
    o.close();
    if (!o) throw std::runtime_error("rebuild: failed closing " + path.string());
}

// `staged` may be given so the path is known in advance (the interrupt handler removes it).
void copy_file(const std::string& from, const std::string& to,
               const std::filesystem::path& staged_in = std::filesystem::path()) {
    // Staged then renamed, so a pass-through never truncates its own input, even when -i and -o name one
    // file.
    const std::filesystem::path staged = staged_in.empty() ? staging_path(to) : staged_in;
    StagedFile guard(staged);
    if (is_gzip_file(from) == wants_gzip(to)) {
        // Same codec on both sides: a byte copy, so a pass-through stays byte-identical.
        std::ifstream in(from, std::ios::binary);
        if (!in) throw std::runtime_error("rebuild: cannot open input: " + from);
        std::ofstream out(staged, std::ios::binary);
        if (!out) throw std::runtime_error("rebuild: cannot open output: " + staged.string());
        out << in.rdbuf();
        out.flush();
        if (!out) throw std::runtime_error("rebuild: failed writing " + staged.string());
    } else {
        // gzread passes a plain file through unchanged, so one reader serves both directions.
        gzFile in = gzopen(from.c_str(), "rb");
        if (in == nullptr) throw std::runtime_error("rebuild: cannot open input: " + from);
        std::string text;
        char buf[1 << 16];
        int n = 0;
        while ((n = gzread(in, buf, sizeof(buf))) > 0) text.append(buf, static_cast<std::size_t>(n));
        const bool read_ok = n == 0;
        gzclose(in);
        // A graph has at least a header line: reading nothing is a failed read, not an empty graph.
        if (!read_ok || text.empty()) throw std::runtime_error("rebuild: failed reading " + from);
        write_text(staged, text, wants_gzip(to));
    }
    commit_staged(staged, to);
    guard.release();
}

#if !defined(_WIN32)
// On SIGINT, SIGTERM or SIGHUP the run removes its scratch FASTAs and staged outputs. Every path it can
// create is listed before the handler is armed and never changes while it is, so the handler only reads it
// and calls unlink/rmdir, which are async-signal-safe; it then restores the default action and re-raises, so
// the exit status names the signal. A signal the caller set to be ignored (nohup's SIGHUP) stays ignored.
// SIGKILL cannot be caught.
class InterruptCleanup {
public:
    // `dirs` deepest first: the scratch directory, then any parents the run created for it. Each is removed
    // only if empty (rmdir), so a parent that gained other files is left alone. One armed while another is
    // (rebuild's second pass) takes over the other's paths too, and hands them back when it goes.
    InterruptCleanup(std::vector<std::string> files, std::vector<std::string> dirs)
        : files_(std::move(files)), dirs_(std::move(dirs)),
          prev_files_(g_files), prev_nfiles_(g_nfiles), prev_dirs_(g_dirs), prev_ndirs_(g_ndirs) {
        for (std::size_t i = 0; i < prev_nfiles_; ++i) files_.emplace_back(prev_files_[i]);
        for (std::size_t i = 0; i < prev_ndirs_; ++i) dirs_.emplace_back(prev_dirs_[i]);
        ptrs_.reserve(files_.size());
        for (const std::string& f : files_) ptrs_.push_back(f.c_str());
        dir_ptrs_.reserve(dirs_.size());
        for (const std::string& d : dirs_) dir_ptrs_.push_back(d.c_str());
        g_files = ptrs_.data();
        g_nfiles = ptrs_.size();
        g_dirs = dir_ptrs_.data();
        g_ndirs = dir_ptrs_.size();
        struct sigaction sa {};
        sa.sa_handler = &InterruptCleanup::on_signal;
        sigemptyset(&sa.sa_mask);
        for (std::size_t i = 0; i < kSignals.size(); ++i) {
            sigaction(kSignals[i], nullptr, &old_[i]);
            if (old_[i].sa_handler == SIG_IGN) continue;
            sigaction(kSignals[i], &sa, nullptr);
        }
    }
    ~InterruptCleanup() {
        for (std::size_t i = 0; i < kSignals.size(); ++i) sigaction(kSignals[i], &old_[i], nullptr);
        // Counts to zero first, so a signal arriving midway never pairs a list with a longer count.
        g_nfiles = 0;
        g_ndirs = 0;
        g_files = prev_files_;
        g_dirs = prev_dirs_;
        g_nfiles = prev_nfiles_;
        g_ndirs = prev_ndirs_;
    }
    InterruptCleanup(const InterruptCleanup&) = delete;
    InterruptCleanup& operator=(const InterruptCleanup&) = delete;

private:
    static void on_signal(int sig) {
        for (std::size_t i = 0; i < g_nfiles; ++i) unlink(g_files[i]);
        for (std::size_t i = 0; i < g_ndirs; ++i) rmdir(g_dirs[i]);
        std::signal(sig, SIG_DFL);
        std::raise(sig);
    }
    static constexpr std::array<int, 3> kSignals{SIGINT, SIGTERM, SIGHUP};
    static inline const char* const* volatile g_files = nullptr;
    static inline volatile std::size_t g_nfiles = 0;
    static inline const char* const* volatile g_dirs = nullptr;
    static inline volatile std::size_t g_ndirs = 0;
    std::vector<std::string> files_;
    std::vector<std::string> dirs_;
    std::vector<const char*> ptrs_;
    std::vector<const char*> dir_ptrs_;
    const char* const* prev_files_;
    std::size_t prev_nfiles_;
    const char* const* prev_dirs_;
    std::size_t prev_ndirs_;
    std::array<struct sigaction, 3> old_{};
};
#endif

// One-record FASTAs, deleted on scope exit: mg_ggen consumes files, so haplotypes round-trip through
// disk. Named by rank so the progressive order stays visible when debugging a run.
class FastaScratch {
public:
    explicit FastaScratch(const std::string& dir) : dir_(dir) {
        // create_directory (not create_directories) so an EXISTING directory is an error rather than
        // something this object will delete on scope exit. The name carries a random tag, so two
        // concurrent runs cannot collide and neither can adopt a directory it did not make.
        std::error_code ec;
        // Parents --tmp-dir names that do not exist yet are created, and removed again on exit if still
        // empty.
        for (std::filesystem::path p = std::filesystem::path(dir_).parent_path();
             !p.empty() && !std::filesystem::exists(p, ec); p = p.parent_path()) {
            made_parents_.push_back(p);
            if (p == p.parent_path()) break;
        }
        std::filesystem::create_directories(std::filesystem::path(dir_).parent_path(), ec);
        if (!std::filesystem::create_directory(dir_, ec) || ec)
            throw std::runtime_error("rebuild: cannot create scratch directory " + dir_ +
                                     (ec ? " (" + ec.message() + ")" : " (already exists)"));
    }
    ~FastaScratch() {
        std::error_code ec;
        std::filesystem::remove_all(dir_, ec);
        for (const std::filesystem::path& p : made_parents_) {   // deepest first; only if empty
            std::error_code pec;
            if (!std::filesystem::is_empty(p, pec) || pec) break;
            std::filesystem::remove(p, pec);
        }
        if (ec) {
            // No shell fallback: `rm -rf` on a path this class did not necessarily create is a worse
            // outcome than a leftover directory, and the path is interpolated into a command line.
            std::cerr << "[rebuild " << hms() << "] warning: could not remove scratch dir " << dir_
                      << " (" << ec.message() << ")\n";
        }
    }
    FastaScratch(const FastaScratch&) = delete;
    FastaScratch& operator=(const FastaScratch&) = delete;

    // `rank` is the position in richness order. Known before anything is written, so the interrupt
    // handler can be given every path in advance.
    std::string path_for(std::size_t rank) const {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "/h%04zu.fa", rank);
        return dir_ + buf;
    }

    std::string write(std::size_t rank, const std::string& name, const std::string& seq) {
        const std::string path = path_for(rank);
        std::ofstream o(path);
        if (!o) throw std::runtime_error("rebuild: cannot write " + path);
        o << '>' << name << '\n';
        for (std::size_t i = 0; i < seq.size(); i += 60) o << seq.substr(i, 60) << '\n';
        // A failed write (a full disk) stops the run rather than handing minigraph a truncated haplotype.
        o.flush();
        if (!o) throw std::runtime_error("rebuild: failed writing " + path);
        o.close();
        if (!o) throw std::runtime_error("rebuild: failed closing " + path);
        return path;
    }

    const std::vector<std::filesystem::path>& made_parents() const { return made_parents_; }

private:
    std::string dir_;
    std::vector<std::filesystem::path> made_parents_;   // created for dir_, deepest first
};

// Owns the gfa_t so every early return frees it.
class GfaHandle {
public:
    explicit GfaHandle(gfa_t* g) : g_(g) {}
    ~GfaHandle() { if (g_) gfa_destroy(g_); }
    GfaHandle(const GfaHandle&) = delete;
    GfaHandle& operator=(const GfaHandle&) = delete;
    gfa_t* get() const { return g_; }
    void reset(gfa_t* g) { if (g_) gfa_destroy(g_); g_ = g; }

private:
    gfa_t* g_;
};

// The ACGT bases of a haplotype's end that the walk's bases past the alignment chain do not hold. Both are read
// away from the chain. The walk part is aligned end to end against as much of the haplotype part as it could
// hold (its length, a tenth more and 50, in ACGT bases), with affine gaps (minimap2's scoring), so an insertion is one gap
// wherever it lies and haplotype bases past the walk's end are a gap at the end. Lost are the haplotype's bases
// past that stretch, those opposite a gap (bases the walk lacks), and those in foreign sequence (below).
// Scattered mismatches and short gaps are differences the rebuild collapses, and walk bases with no haplotype
// base opposite are never a loss.
long long unmatched_end_bases(const std::string& walk, const std::string& hap) {
    constexpr long long kForeign = 16;
    constexpr std::size_t kLongGap = 16;
    constexpr int kGap = 9;   // marks a haplotype-only column
    constexpr std::size_t kLongest = 50000;   // past this, the walk's part pairs base for base
    auto up = [](char c) { return static_cast<char>(std::toupper(static_cast<unsigned char>(c))); };
    auto acgt = [&](char c) { c = up(c); return c == 'A' || c == 'C' || c == 'G' || c == 'T'; };
    auto acgt_in = [&](std::size_t a, std::size_t b) {
        long long n = 0;
        for (std::size_t i = a; i < std::min(b, hap.size()); ++i) n += acgt(hap[i]);
        return n;
    };
    // the stretch holds that many ACGT bases: an N run there holds no variation and does not use up the room
    std::size_t span = 0;
    for (std::size_t room = walk.size() + walk.size() / 10 + 50; span < hap.size() && room > 0; ++span)
        room -= acgt(hap[span]);
    while (span < hap.size() && !acgt(hap[span])) ++span;
    if (walk.empty() || span == 0) return acgt_in(0, hap.size());
    if (walk.size() > kLongest) return acgt_in(walk.size(), hap.size());
    auto code = [&](char c) -> std::uint8_t {
        switch (up(c)) { case 'A': return 0; case 'C': return 1; case 'G': return 2; case 'T': return 3; default: return 4; }
    };
    std::vector<std::uint8_t> w(walk.size()), t(span);
    for (std::size_t i = 0; i < walk.size(); ++i) w[i] = code(walk[i]);
    for (std::size_t i = 0; i < span; ++i) t[i] = code(hap[i]);
    std::int8_t mat[25];
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < 5; ++j) mat[i * 5 + j] = static_cast<std::int8_t>((i == 4 || j == 4) ? -1 : i == j ? 2 : -4);
    const int band = static_cast<int>(std::max(span, walk.size()) - std::min(span, walk.size())) + 100;
    ksw_extz_t ez;
    std::memset(&ez, 0, sizeof(ez));
    ksw_extd2_sse(nullptr, static_cast<int>(w.size()), w.data(), static_cast<int>(t.size()), t.data(), 5, mat, 16, 4, 30,
                  1, band, -1, 0, 0, &ez);
    // per column: the haplotype base it holds (or none), and its weight toward a foreign stretch
    std::vector<std::pair<std::size_t, int>> cols;
    std::size_t i = 0, j = 0;
    for (int k = 0; k < ez.n_cigar; ++k) {
        const int op = static_cast<int>(ez.cigar[k] & 0xf);
        for (std::uint32_t n = ez.cigar[k] >> 4; n > 0; --n) {
            if (op == 0) {
                cols.emplace_back(j, !acgt(hap[j]) ? 0 : up(walk[i]) == up(hap[j]) ? -1 : 2);
                ++i; ++j;
            } else if (op == 2) {
                cols.emplace_back(j, kGap);
                ++j;
            } else {
                cols.emplace_back(SIZE_MAX, 0);
                ++i;
            }
        }
    }
    std::vector<char> lost(j, 0);
    for (int k = 0, at = 0; k < ez.n_cigar; ++k) {   // opposite a gap
        const int len = static_cast<int>(ez.cigar[k] >> 4);
        if ((ez.cigar[k] & 0xf) == 2)
            for (int m = 0; m < len; ++m) lost[cols[static_cast<std::size_t>(at + m)].first] = 1;
        at += len;
    }
    std::free(ez.cigar);
    // Foreign sequence: a stretch where a mismatch weighs 2, a match -1 and a gap column 1, whose weight climbs
    // kForeign or more. Unrelated sequences match at a third of their pairs or so, where the alignment can shift
    // to make them; SNPs, even dense, at a fifth or fewer. The stretch runs from where the weight starts to climb
    // to its top, so it takes in no matching flank. A gap of kLongGap or more weighs nothing: one the walk lacks
    // is counted above, one the haplotype lacks holds no haplotype base, and neither should join two stretches.
    std::vector<int> weight(cols.size(), 0);
    auto gap = [&](std::size_t k) { return cols[k].first == SIZE_MAX ? 1 : cols[k].second == kGap ? 2 : 0; };
    for (std::size_t k = 0; k < cols.size();) {
        if (!gap(k)) { weight[k] = cols[k].second; ++k; continue; }
        std::size_t e = k;
        while (e < cols.size() && gap(e) == gap(k)) ++e;
        if (e - k < kLongGap) std::fill(weight.begin() + static_cast<long>(k), weight.begin() + static_cast<long>(e), 1);
        k = e;
    }
    long long run = 0, top = 0;
    std::size_t from = 0, top_at = 0;
    auto close = [&]() {
        if (top >= kForeign)
            for (std::size_t k = from; k <= top_at; ++k)
                if (cols[k].first != SIZE_MAX) lost[cols[k].first] = 1;
    };
    for (std::size_t k = 0; k < cols.size(); ++k) {
        if (run == 0 && weight[k] <= 0) continue;
        if (run == 0) { from = k; top = 0; }
        run += weight[k];
        if (run > top) { top = run; top_at = k; }
        if (run <= 0) { close(); run = 0; }
    }
    if (run > 0) close();
    long long n = acgt_in(j, hap.size());
    for (std::size_t k = 0; k < j; ++k) n += lost[k] && acgt(hap[k]);
    return n;
}

} // namespace

namespace {

// A threshold as given: 0.975 stays 0.975, 0.95 stays 0.95.
std::string fmt_short(double v) {
    std::ostringstream o;
    o << v;
    return o.str();
}

std::string fmt2(double v) {
    char b[32];
    std::snprintf(b, sizeof(b), "%.4f", v);
    return b;
}

// One graph's shape, the same for the input and the rebuilt graph: hub nodes decide healthy or
// pathological; the rest is context.
std::string shape_line(const std::string& graph, bool pathological, std::size_t hubs, std::size_t hub_degree,
                       std::size_t min_hubs, std::size_t maxdeg, double density, double inflation, bool explain) {
    char b[64];
    std::snprintf(b, sizeof(b), "%.2f", inflation);
    return graph + " is " + (pathological ? "pathological" : "healthy") + ": hub nodes: " + std::to_string(hubs) +
           (explain ? " (" + std::to_string(hub_degree) + "+ neighbours on one side; " + std::to_string(min_hubs) +
                          " make a graph pathological)"
                    : std::string()) +
           "; max degree: " + std::to_string(maxdeg) + "; nodes/kb: " + std::to_string(static_cast<long>(density)) +
           "; graph bp / median haplotype bp: " + b;
}

// Enough digits to reproduce a threshold decision from the report.
std::string fmt_exact(double v) {
    char b[40];
    std::snprintf(b, sizeof(b), "%.10g", v);
    return b;
}

} // namespace

std::string rebuild_report_path(const std::string& out_path) {
    for (const char* ext : {".gfa.gz", ".gfa"}) {
        const std::string e = ext;
        if (out_path.size() > e.size() && out_path.compare(out_path.size() - e.size(), e.size(), e) == 0)
            return out_path.substr(0, out_path.size() - e.size()) + ".rebuild.report.tsv";
    }
    return out_path + ".rebuild.report.tsv";
}

RebuildSummary rebuild_graph(const RebuildOptions& options) {
    RebuildSummary sum;
    // The second pass repeats the first's lines up to the generated graph, which it reads back: those are not
    // printed again.
    bool echo = options.generated_gfa.empty();
    auto step = [&](const std::string& msg) {
        if (!options.quiet && echo) std::cerr << "[rebuild " << hms() << "] " << msg << std::endl;
    };
    // Warnings and a rejection are shown even with -q: they change what the output is.
    auto warn = [&](const std::string& msg) {
        if (echo) std::cerr << "[rebuild " << hms() << "] " << msg << std::endl;
    };

    step("reading: " + options.gfa_path);
    ParseGfaOptions po;
    po.include_paths = true;
    po.include_sequences = true;
    const Graph g = parse_gfa(options.gfa_path, po);
    validate_graph_paths(g, "rebuild", true, true);
    // Resolved before any work, with the resolver every module uses (match_reference_name), so a wrong name
    // fails at once.
    const std::string ref_name = resolve_reference_path_name(g, options.reference_path, "rebuild");

    // The report describes a rebuild: one left by an earlier run is removed when none is performed. A
    // directory where the graph or its report goes is refused up front.
    if (!options.out_path.empty()) {
        std::error_code dir_ec;
        for (const std::string& p : {options.out_path, rebuild_report_path(options.out_path)})
            if (std::filesystem::is_directory(p, dir_ec))
                throw std::runtime_error("rebuild: output '" + p + "' is an existing directory");
    }
    auto drop_stale_report = [&]() {
        if (options.out_path.empty()) return;
        const std::string stale = rebuild_report_path(options.out_path);
        std::error_code ec;
        if (std::filesystem::remove(stale, ec)) step("removed a report left by an earlier run: " + stale);
    };

    sum.raw_nodes = g.nodes.size();
    sum.haplotypes = g.paths.size();
    degree_stats(g, options.hub_degree, sum.raw_hubs, sum.raw_maxdeg, &sum.raw_selfloops);
    // Lengths from the node table: spelling every haplotype here is most of a pass-through's time and
    // memory, for a graph that is then copied unchanged.
    std::size_t max_len = 1;
    std::vector<std::size_t> path_len;
    for (const PathRecord& p : g.paths) {
        std::size_t len = 0;
        for (const PathStep& s : p.steps) len += g.nodes.at(s.node_id).sequence.size();
        max_len = std::max(max_len, len);
        path_len.push_back(len);
    }
    // Graph sequence over the median haplotype: about 1 for a well-merged locus, well above for one whose
    // recurring sequence was left unmerged. Logged to explain the gate's decision; not part of it.
    std::size_t graph_bp = 0;
    for (const auto& kv : g.nodes) graph_bp += kv.second.sequence.size();
    std::sort(path_len.begin(), path_len.end());
    const double inflation = path_len.empty() || path_len[path_len.size() / 2] == 0
                                 ? 0.0
                                 : static_cast<double>(graph_bp) / static_cast<double>(path_len[path_len.size() / 2]);
    sum.raw_density = static_cast<double>(sum.raw_nodes) / (static_cast<double>(max_len) / 1000.0);
    const std::size_t median_len = path_len.empty() ? 0 : path_len[path_len.size() / 2];
    std::size_t raw_edges = 0;
    for (const auto& kv : g.nodes) raw_edges += kv.second.start.size() + kv.second.end.size();
    raw_edges /= 2;

    // ---- Criteria A: the pathology gate ----
    sum.pathological = sum.raw_hubs >= options.min_hubs;
    step("input: haplotypes: " + std::to_string(sum.haplotypes) + "; nodes: " + std::to_string(sum.raw_nodes) +
         "; edges: " + std::to_string(raw_edges) + "; self-loops: " + std::to_string(sum.raw_selfloops));
    step(shape_line("input", sum.pathological, sum.raw_hubs, options.hub_degree, options.min_hubs, sum.raw_maxdeg,
                    sum.raw_density, inflation, true));
    if (!sum.pathological && !options.force) {
        step("passing the graph through unchanged (--force rebuilds a healthy graph)");
        if (!options.out_path.empty()) {
            // The staged copy is removed on an interrupt too, as a rebuild's outputs are.
            const std::filesystem::path staged = staging_path(options.out_path);
#if !defined(_WIN32)
            const InterruptCleanup interrupt_cleanup({staged.string()}, {});
#endif
            copy_file(options.gfa_path, options.out_path, staged);
        }
        drop_stale_report();
        sum.out_nodes = sum.raw_nodes;
        return sum;
    }
    if (g.paths.empty()) throw std::runtime_error("rebuild: no haplotype paths in " + options.gfa_path);
    if (!sum.pathological)
        warn("warning: --force on a healthy graph; variation shorter than --min-var (" +
             std::to_string(options.min_var) + " bp) is not kept");
    sum.ran = true;

    // Lexicographic, diversity first: (distinct k-mers, total k-mers) descending; abundance only separates
    // haplotypes carrying the same distinct sequence. Threads are capped at four per core, here and in the
    // per-haplotype phases: minigraph sets up per-thread state for every thread it is given.
    const unsigned hw = std::max(1u, std::thread::hardware_concurrency());
    unsigned nthreads = options.threads > 0
                            ? static_cast<unsigned>(std::min<std::size_t>(options.threads, 4u * hw))
                            : hw;
    if (options.threads > nthreads)
        step("--threads " + std::to_string(options.threads) + " capped at " + std::to_string(nthreads) +
             " (four per core)");
    step("ordering haplotypes by k-mer richness: " + std::to_string(sum.haplotypes));
    std::vector<std::string> hseq(g.paths.size());
    std::vector<std::pair<std::size_t, std::size_t>> score(g.paths.size());
    std::vector<std::size_t> idx(g.paths.size());
    run_parallel(g.paths.size(), nthreads, [&](std::size_t i) {
        hseq[i] = spell(g, g.paths[i]);
        // minigraph stores uppercase: masked (lowercase) input is uppercased so it does not read as
        // mismatches.
        for (char& c : hseq[i]) c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
        const KmerStats st = kmer_stats(hseq[i], options.kmer);
        score[i] = {st.distinct, st.total};
    });
    for (std::size_t i = 0; i < g.paths.size(); ++i) idx[i] = i;
    // Richness descending, then path name, then original index, so the seed and the build are deterministic:
    // ties are common among near-identical haplotypes.
    std::sort(idx.begin(), idx.end(), [&](std::size_t a, std::size_t b) {
        if (score[a] != score[b]) return score[a] > score[b];
        if (g.paths[a].name != g.paths[b].name) return g.paths[a].name < g.paths[b].name;
        return a < b;
    });
    // With a reference named, it seeds the graph and the rest follow in richness order: the seed is the one
    // haplotype the generator reproduces exactly, and downstream positions are read off the reference.
    if (!ref_name.empty()) {
        const auto it = std::find_if(idx.begin(), idx.end(),
                                     [&](std::size_t i) { return g.paths[i].name == ref_name; });
        std::rotate(idx.begin(), it, it + 1);
    }
    sum.seed = g.paths[idx.front()].name;
    step("seed: " + sum.seed + (ref_name.empty() ? "" : " (the reference, -r)"));

    const int n_threads = static_cast<int>(nthreads);

    // ---- progressive graph generation, in richness order ----
    // minigraph reads the first file as the seed graph and augments it with each subsequent one in turn,
    // so passing our richness order drives the progressive build directly.
    // Scratch FASTAs live beside --out by default, or under --tmp-dir in a dedicated subfolder (named
    // from the output basename) so cleanup removes only what we created, never the whole --tmp-dir.
    std::random_device scratch_rd;
    const std::string scratch_tag = ".rebuild.tmp." + std::to_string(scratch_rd());
    std::string scratch_dir = options.out_path + scratch_tag;
    if (!options.tmp_dir.empty()) {
        std::string base = options.out_path;
        const auto slash = base.find_last_of('/');
        if (slash != std::string::npos) base = base.substr(slash + 1);
        scratch_dir = options.tmp_dir + "/" + base + scratch_tag;
    }
    FastaScratch scratch(scratch_dir);
    // Staged names are fixed now so the interrupt handler knows every path this run can create.
    const std::string report_path = options.out_path.empty() ? std::string()
                                                            : rebuild_report_path(options.out_path);
    const std::filesystem::path out_staged =
        options.out_path.empty() ? std::filesystem::path() : staging_path(options.out_path);
    const std::filesystem::path report_staged =
        report_path.empty() ? std::filesystem::path() : staging_path(report_path);
#if !defined(_WIN32)
    std::vector<std::string> created;
    created.reserve(idx.size() + 2);
    for (std::size_t rank = 0; rank < idx.size(); ++rank) created.push_back(scratch.path_for(rank));
    if (!out_staged.empty()) created.push_back(out_staged.string());
    if (!report_staged.empty()) created.push_back(report_staged.string());
    created.push_back(scratch_dir + "/generated.gfa");   // saved for a second pass, below
    std::vector<std::string> dirs{scratch_dir};
    for (const std::filesystem::path& p : scratch.made_parents()) dirs.push_back(p.string());
    const InterruptCleanup interrupt_cleanup(std::move(created), std::move(dirs));
#endif
    std::vector<std::string> fasta_paths;
    fasta_paths.reserve(idx.size());
    for (std::size_t rank = 0; rank < idx.size(); ++rank) {
        const std::size_t h = idx[rank];
        fasta_paths.push_back(scratch.write(rank, g.paths[h].name, hseq[h]));
    }

    mg_idxopt_t ipt;
    mg_mapopt_t opt;
    mg_ggopt_t gpt;
    mg_opt_set(nullptr, &ipt, &opt, &gpt);
    if (mg_opt_set("ggs", &ipt, &opt, &gpt) < 0) throw std::runtime_error("rebuild: mg_opt_set(ggs) failed");
    gpt.min_var_len = static_cast<int>(options.min_var);
    // --base-align: minigraph's base-level graph generation (its -c); finds events the anchor-only mode
    // misses, at about twice the CPU.
    if (options.base_align) opt.flag |= MG_M_CIGAR;
    // minigraph's length gates assume chromosome-scale input: a chain under min_map_len (100 kb by default)
    // is dropped before min_var_len is consulted, so a locus would augment nothing. Both gates are scaled to
    // the locus, from the median haplotype length and capped by the shortest, so every haplotype can
    // contribute a chain covering half of itself. Loci above the defaults keep minigraph's values.
    std::vector<std::size_t> hlen;
    hlen.reserve(hseq.size());
    for (const std::string& h : hseq) if (!h.empty()) hlen.push_back(h.size());
    if (hlen.empty()) throw std::runtime_error("rebuild: every haplotype spelled an empty sequence");
    std::sort(hlen.begin(), hlen.end());
    const std::size_t min_len = hlen.front();
    const std::size_t med_len = hlen[hlen.size() / 2];
    const std::size_t seed_len = hseq[idx.front()].size();
    // The automatic value is never above half the shortest haplotype (nor under 500 bp); a requested one above
    // it is warned about, not changed.
    const int cap = static_cast<int>(std::max<std::size_t>(500, min_len / 2));
    if (options.min_align_len > 0) {
        gpt.min_map_len = static_cast<int>(options.min_align_len);
        if (gpt.min_map_len > cap) {
            warn("warning: --min-align-len " + std::to_string(options.min_align_len) +
                 " exceeds half the shortest haplotype (" + std::to_string(min_len) +
                 " bp); haplotypes shorter than twice it cannot clear the gate");
        }
        if (static_cast<std::size_t>(gpt.min_depth_len) > options.min_align_len) {
            gpt.min_depth_len = static_cast<int>(std::max<std::size_t>(500, options.min_align_len / 5));
        }
    } else {
        const int base = (med_len < 2 * static_cast<std::size_t>(gpt.min_map_len))
                             ? static_cast<int>(std::max<std::size_t>(1000, med_len / 2))
                             : gpt.min_map_len;
        const int chosen = std::min(base, cap);
        if (chosen != gpt.min_map_len) {
            gpt.min_map_len = chosen;
            if (med_len < 2 * static_cast<std::size_t>(gpt.min_depth_len))
                gpt.min_depth_len = static_cast<int>(std::max<std::size_t>(500, med_len / 10));
            gpt.min_depth_len = std::min(gpt.min_depth_len, gpt.min_map_len);
        }
    }
    step("haplotype lengths: " + std::to_string(min_len) + "-" + std::to_string(hlen.back()) + " bp (median: " +
         std::to_string(med_len) + "; seed: " + std::to_string(seed_len) + "); min alignment length: " +
         std::to_string(gpt.min_map_len) + " bp" + (options.min_align_len > 0 ? " (requested)" : ""));
    // Keep minigraph itself quiet (its per-sample logs are one block per haplotype, far too noisy at
    // cohort scale) and print our own throttled counter instead. mg_ggen augments one file per call
    // internally, so driving that loop here costs nothing extra and lets us report progress.
    mg_verbose = 1;

    GfaHandle handle(gfa_read(fasta_paths.front().c_str()));
    if (handle.get() == nullptr) {
        throw std::runtime_error("rebuild: minigraph could not seed from " + fasta_paths.front());
    }
    const std::size_t total = fasta_paths.size();
    step("generating graph: seed + " + std::to_string(total - 1) + " haplotypes (minigraph)");
    if (!options.generated_gfa.empty()) {
        handle.reset(gfa_read(options.generated_gfa.c_str()));
        if (handle.get() == nullptr) throw std::runtime_error("rebuild: cannot read back " + options.generated_gfa);
        echo = true;
        step("second pass, on the graph the first generated");
    } else {
        cli::ProgressBar bar(options.quiet ? std::string() : "[rebuild " + hms() + "] adding haplotypes", total - 1);
        for (std::size_t i = 1; i < fasta_paths.size(); ++i) {
            const char* f = fasta_paths[i].c_str();
            if (mg_ggen(handle.get(), 1, &f, &ipt, &opt, &gpt, n_threads) != 0) {
                throw std::runtime_error("rebuild: mg_ggen failed on " + fasta_paths[i]);
            }
            bar.tick();
        }
    }
    gfa_t* out = handle.get();
    // Saved for a second pass (below): generation is the slow step and its result does not change.
    const std::string generated_copy = scratch_dir + "/generated.gfa";
    if (options.generated_gfa.empty() && options.fold_modules) {
        std::FILE* fp = std::fopen(generated_copy.c_str(), "w");
        if (fp == nullptr) throw std::runtime_error("rebuild: cannot write " + generated_copy);
        gfa_print(out, fp, 0);
        if (std::fclose(fp) != 0) throw std::runtime_error("rebuild: cannot write " + generated_copy);
    }
    {
        // minigraph keeps each link as two arcs, one per direction, sharing a link id.
        std::unordered_set<std::uint64_t> link_ids;
        for (std::uint64_t k = 0; k < out->n_arc; ++k) link_ids.insert(out->arc[k].link_id);
        step("generated graph: nodes: " + std::to_string(out->n_seg) + "; edges: " + std::to_string(link_ids.size()));
    }

    // A graph with no arcs carries no variation of its own; whether the walks add any (a tip, a kept gap, an
    // indel the generator left out) is known once they are recovered and repaired, below.
    const bool generated_flat = out->n_arc == 0 || out->n_seg <= 1;

    // ---- path recovery: map every haplotype back to the graph, in memory (no GAF round-trip) ----
    // Each mg_llchain_t.v is an ORIENTED vertex (seg<<1|strand), which is what lets inversion
    // traversals survive into the P lines.
    mg_idxopt_t ipt_map;
    mg_mapopt_t opt_map;
    mg_ggopt_t gpt_map;
    mg_opt_set(nullptr, &ipt_map, &opt_map, &gpt_map);
    if (mg_opt_set("asm", &ipt_map, &opt_map, &gpt_map) < 0) {
        throw std::runtime_error("rebuild: mg_opt_set(asm) failed");
    }
    // Base-level alignment, so matched cover counts matching bases; a chain's mlen sums minimizer-anchor
    // spans, where each SNP costs about a k-mer of cover. The walks themselves do not change.
    opt_map.flag |= MG_M_CIGAR;
    mg_idx_t* gi = mg_index(out, &ipt_map, n_threads, &opt_map);
    if (gi == nullptr) throw std::runtime_error("rebuild: mg_index failed");
    step("aligning haplotypes back to the rebuilt graph: " + std::to_string(g.paths.size()));

    // Parallel like minigraph's own kt_for: index and options shared read-only, one mg_tbuf_t per
    // worker (what makes mg_map re-entrant). Results go to distinct indices, so output is
    // thread-order independent.
    // A recovered walk is a list of steps, each an oriented generator vertex and the part of it kept,
    // [a, b) in segment coordinates -- whole vertices except where a walk starts, ends, or (for a
    // joined inversion) turns.
    struct Step { std::uint32_t v, a, b; bool junction; };   // junction: a new link enters this step
    std::vector<std::uint32_t> seg_len(out->n_seg);
    for (std::uint32_t i = 0; i < out->n_seg; ++i) seg_len[i] = out->seg[i].len;
    // Sequence rebuild adds itself, beyond the generator's segments (ids continue after them): a
    // haplotype's own start or end where it runs past the graph (`tip`), and an assembly gap kept as
    // the haplotype's N bases (`gap`).
    std::vector<std::string> extra_seq, extra_kind;
    auto seg_seq = [&](std::uint32_t sg) -> const char* {
        return sg < out->n_seg ? out->seg[sg].seq : extra_seq[sg - out->n_seg].data();
    };
    auto add_extra = [&](std::string sq, const char* kind) {
        extra_seq.push_back(std::move(sq));
        extra_kind.push_back(kind);
        seg_len.push_back(static_cast<std::uint32_t>(extra_seq.back().size()));
        return static_cast<std::uint32_t>(seg_len.size() - 1);
    };
    // Assembly gaps: runs of 10 or more N (or other non-ACGT) bases.
    auto gap_runs = [](const std::string& q) {
        std::vector<std::pair<long long, long long>> runs;
        long long st = -1;
        for (long long i = 0; i <= static_cast<long long>(q.size()); ++i) {
            const bool n = i < static_cast<long long>(q.size()) && q[i] != 'A' && q[i] != 'C' && q[i] != 'G' && q[i] != 'T';
            if (n && st < 0) st = i;
            if (!n && st >= 0) { if (i - st >= 10) runs.push_back({st, i}); st = -1; }
        }
        return runs;
    };
    std::vector<std::vector<Step>> steps(g.paths.size());
    std::vector<int> inversion(g.paths.size(), 0);         // 1 candidate not verified, 2 joined
    std::vector<double> cover(g.paths.size(), 0.0);        // envelope (qe-qs)/len
    std::vector<double> matched(g.paths.size(), 0.0);      // matching bases / len
    std::vector<double> chain_id(g.paths.size(), 0.0);     // mlen / blen
    // Bases the walk spells beyond the haplotype's own start / end. A walk lists whole segments, and
    // the generator never cuts a segment where a haplotype ends, so a haplotype stopping partway into
    // its first or last segment is still given the rest of it. Read from the chain: ps/pe are where
    // the alignment starts and ends on the walk, qs/qe on the haplotype.
    std::vector<long long> lead_over(g.paths.size(), 0), trail_over(g.paths.size(), 0);
    // Where the alignment starts and ends on the walk (ps, pe) and on the haplotype (qs, qe).
    std::vector<std::array<long long, 4>> chain_span(g.paths.size(), {0, 0, 0, 0});
    // Largest |walk length - haplotype length| a haplotype may show: max(1 kb, 4 x --min-var, 0.5% of its
    // length). A lost SV of a few kb exceeds it, while the indels under --min-var collapsed by design, which
    // add up along a long haplotype, do not.
    auto length_bound_for = [&](std::size_t len) {
        return std::max<long long>({1000, 4 * static_cast<long long>(options.min_var),
                                    static_cast<long long>(len) / 200});
    };
    auto spell_steps = [&](const std::vector<Step>& w) {
        std::string sp;
        for (const Step& st : w) {
            const char* sq = seg_seq(st.v >> 1);
            const std::string part = sq ? std::string(sq + st.a, st.b - st.a) : std::string();
            sp += (st.v & 1) ? reverse_complement(part) : part;
        }
        return sp;
    };
    // Identity of `a` against `b`: 1 - edit distance / the longer length. Distance only: no traceback is
    // needed.
    auto identity_of = [](const std::string& a, const std::string& b, std::size_t* edits = nullptr) {
        const std::size_t ed = nw_distance(a, b);
        if (edits != nullptr) *edits = ed;
        const std::size_t len = std::max(a.size(), b.size());
        return len ? std::max(0.0, 1.0 - static_cast<double>(ed) / static_cast<double>(len)) : 0.0;
    };
    // A mapping chain, copied out of minigraph's result so it outlives the mapping call: haplotype
    // span (qs, qe), span on its own walk (ps, pe of plen), its vertices, and its base alignment.
    struct Chain {
        long long qs = 0, qe = 0, ps = 0, pe = 0, plen = 0;
        double mlen = 0, blen = 0;
        std::vector<std::uint32_t> v;
        std::vector<std::uint64_t> cigar;
        bool has_cigar = false;
    };
    auto make_chain = [](const mg_gchains_t& G, const mg_gchain_t& c) {
        Chain ch;
        ch.qs = c.qs; ch.qe = c.qe; ch.ps = c.ps; ch.pe = c.pe; ch.plen = c.plen;
        ch.mlen = c.p != nullptr ? c.p->mlen : c.mlen;   // the base alignment; the anchor estimate if none
        ch.blen = c.p != nullptr ? c.p->blen : c.blen;
        for (int32_t j = 0; j < c.cnt; ++j) ch.v.push_back(G.lc[c.off + j].v);
        if (c.p != nullptr) {
            ch.has_cigar = true;
            ch.cigar.assign(c.p->cigar, c.p->cigar + c.p->n_cigar);
        }
        return ch;
    };
    // Walk offset (from the chain's first vertex) of haplotype position x, read off the base alignment;
    // query_at is the inverse.
    auto path_at = [](const Chain& c, long long x) -> long long {
        if (!c.has_cigar || c.qe <= c.qs)
            return c.ps + (c.qe > c.qs ? (x - c.qs) * (c.pe - c.ps) / (c.qe - c.qs) : 0);
        long long q = c.qs, pp = c.ps;
        for (const std::uint64_t op_len : c.cigar) {
            const long long len = static_cast<long long>(op_len >> 4);
            const int op = static_cast<int>(op_len & 0xf);
            const bool cq = op == 0 || op == 1 || op == 7 || op == 8, cp = op == 0 || op == 2 || op == 3 || op == 7 || op == 8;
            if (cq && q + len >= x) return pp + (cp ? x - q : 0);
            if (cq) q += len;
            if (cp) pp += len;
        }
        return pp;
    };
    // As path_at, but past any deletion sitting exactly at x: an N run aligns as an insertion beside a
    // deletion of the graph sequence opposite it, and the gap's right edge lies after that deletion.
    auto path_at_right = [](const Chain& c, long long x) -> long long {
        if (!c.has_cigar || c.qe <= c.qs)
            return c.ps + (c.qe > c.qs ? (x - c.qs) * (c.pe - c.ps) / (c.qe - c.qs) : 0);
        long long q = c.qs, pp = c.ps;
        bool reached = false;
        for (const std::uint64_t op_len : c.cigar) {
            const long long len = static_cast<long long>(op_len >> 4);
            const int op = static_cast<int>(op_len & 0xf);
            const bool cq = op == 0 || op == 1 || op == 7 || op == 8, cp = op == 0 || op == 2 || op == 3 || op == 7 || op == 8;
            if (reached) {
                if (cq) return pp;            // the next base of the haplotype: stop here
                pp += len;                    // a deletion at x: step over it
                continue;
            }
            if (cq && q + len > x) return pp + (cp ? x - q : 0);
            if (cq) q += len;
            if (cp) pp += len;
            if (q == x) reached = true;
        }
        return pp;
    };
    auto query_at = [](const Chain& c, long long pp0) -> long long {
        if (!c.has_cigar || c.pe <= c.ps)
            return c.qs + (c.pe > c.ps ? (pp0 - c.ps) * (c.qe - c.qs) / (c.pe - c.ps) : 0);
        long long q = c.qs, pp = c.ps;
        for (const std::uint64_t op_len : c.cigar) {
            const long long len = static_cast<long long>(op_len >> 4);
            const int op = static_cast<int>(op_len & 0xf);
            const bool cq = op == 0 || op == 1 || op == 7 || op == 8, cp = op == 0 || op == 2 || op == 3 || op == 7 || op == 8;
            if (cp && pp + len >= pp0) return q + (cq ? pp0 - pp : 0);
            if (cq) q += len;
            if (cp) pp += len;
        }
        return q;
    };
    // Steps covering walk range [p0, p1) of chain c.
    auto slice = [&](const Chain& c, long long p0, long long p1) {
        std::vector<Step> w;
        long long at = 0;
        for (const std::uint32_t v : c.v) {
            const std::uint32_t len = seg_len[v >> 1];
            const long long s0 = at, s1 = at + len;
            at = s1;
            const long long o0 = std::max(p0, s0) - s0, o1 = std::min(p1, s1) - s0;
            if (o1 <= o0) continue;
            if (v & 1) w.push_back({v, static_cast<std::uint32_t>(len - o1), static_cast<std::uint32_t>(len - o0), false});
            else w.push_back({v, static_cast<std::uint32_t>(o0), static_cast<std::uint32_t>(o1), false});
        }
        return w;
    };
    // The part [from, to) of a walk, in walk offsets. A turn into a step survives only if the step's
    // start is kept; the caller marks the new turns it creates.
    auto trim_steps = [](const std::vector<Step>& w, long long from, long long to) {
        std::vector<Step> kept;
        long long at = 0;
        for (const Step& st : w) {
            const long long len = st.b - st.a, s0 = at, s1 = at + len;
            at = s1;
            const long long o0 = std::max(from, s0) - s0, o1 = std::min(to, s1) - s0;
            if (o1 <= o0) continue;
            Step t = st;
            if (st.v & 1) { t.a = static_cast<std::uint32_t>(st.b - o1); t.b = static_cast<std::uint32_t>(st.b - o0); }
            else { t.a = static_cast<std::uint32_t>(st.a + o0); t.b = static_cast<std::uint32_t>(st.a + o1); }
            t.junction = st.junction && o0 == 0;
            kept.push_back(t);
        }
        return kept;
    };
    auto chain_strand = [&](const Chain& c) {
        long long fwd = 0, rev = 0;
        for (const std::uint32_t v : c.v) (v & 1 ? rev : fwd) += seg_len[v >> 1];
        return rev > fwd ? 1 : 0;
    };
    // A point on the graph (segment, offset in segment coordinates) <-> a walk offset of chain c.
    using GraphPoint = std::pair<std::uint32_t, std::uint32_t>;
    auto to_point = [&](const Chain& c, long long pp) {
        long long at = 0;
        for (const std::uint32_t v : c.v) {
            const std::uint32_t len = seg_len[v >> 1];
            if (pp <= at + len) {
                const auto o = static_cast<std::uint32_t>(pp - at);
                return GraphPoint{v >> 1, (v & 1) ? len - o : o};
            }
            at += len;
        }
        return GraphPoint{0, 0};
    };
    auto from_point = [&](const Chain& c, const GraphPoint& gp, long long lo, long long hi) -> long long {
        long long at = 0;
        for (const std::uint32_t v : c.v) {
            const std::uint32_t len = seg_len[v >> 1];
            if ((v >> 1) == gp.first) {
                const long long pp = at + ((v & 1) ? len - gp.second : gp.second);
                if (pp >= lo && pp <= hi) return pp;
            }
            at += len;
        }
        return -1;
    };
    // ---- large inversions ----
    // The generator adds an inversion only when it lies inside one alignment chain, and chaining
    // bridges roughly 150 kb, so a larger inversion maps as three chains -- flank, reversed middle,
    // flank -- and nothing joins them. Joined here: three chains in haplotype order, the outer two on
    // one strand and the middle on the other, overlapping or abutting (inverted repeats at the
    // breakpoints make them overlap). Each turn is cut inside the overlap, located on the walks through
    // the base alignments, and becomes a new link. The joined walk is kept only if it spells the
    // haplotype within the identity and length bounds, so a rearrangement the sequence does not
    // support is never stitched in; the haplotype then fails as before, labelled `possible_inversion`.
    struct Candidate { Chain a, b, c; };
    std::vector<std::unique_ptr<Candidate>> candidate(g.paths.size());
    // Assembly gaps. A haplotype's N run aligns to whatever graph sequence lies opposite it, so its walk
    // spelled graph bases there -- reference bases downstream, where the assembly says "unknown". The
    // walk is re-routed through a piece holding the haplotype's own N bases instead. A gap inside one
    // chain is located through the base alignment; a gap long enough to split the alignment into two
    // chains on one strand, with the bases between them mostly N, joins the two across a gap piece.
    std::vector<std::unique_ptr<Chain>> gap_chain(g.paths.size());       // the single chain, if gapped
    std::vector<std::vector<std::pair<long long, long long>>> chain_q(g.paths.size());   // every chain's [qs, qe)
    std::vector<long long> clipped(g.paths.size(), 0);   // bases dropped by --clip-to-reference
    std::vector<long long> clipped_left(g.paths.size(), 0), clipped_right(g.paths.size(), 0);
    std::vector<std::size_t> clipped_gap(g.paths.size(), 0);   // N / ambiguity bases among them
    std::vector<std::size_t> input_len(g.paths.size());
    for (std::size_t h = 0; h < g.paths.size(); ++h) input_len[h] = hseq[h].size();
    struct GapPair { Chain a, c; };
    std::vector<std::unique_ptr<GapPair>> gap_pair(g.paths.size());
    auto find_gap_pair = [&](const mg_gchains_t& G, const std::string& q, long long best_span)
        -> std::unique_ptr<GapPair> {
        const int n = std::min<int>(G.n_gc, 20);
        long long best_kept = best_span;
        int bi = -1, bk = -1;
        for (int x = 0; x < n; ++x)
            for (int z = 0; z < n; ++z) {
                const mg_gchain_t &A = G.gc[x], &C = G.gc[z];
                if (x == z || C.qs <= A.qe || C.qe - A.qs <= best_kept) continue;
                long long nb = 0;
                for (long long i = A.qe; i < C.qs; ++i)
                    nb += q[i] != 'A' && q[i] != 'C' && q[i] != 'G' && q[i] != 'T';
                if (nb < 10 || nb * 10 < (C.qs - A.qe) * 8) continue;          // mostly N between them
                const Chain ca = make_chain(G, A), cc = make_chain(G, C);
                if (chain_strand(ca) != chain_strand(cc)) continue;
                // C must carry on along the graph where A stopped, read by the seed offsets its own segments
                // keep (rank 0), so two halves of a rotated or translocated haplotype are not joined across
                // the locus.
                long long xa = -1, xc = -1;
                for (auto it = ca.v.rbegin(); it != ca.v.rend() && xa < 0; ++it)
                    if (out->seg[*it >> 1].rank == 0) xa = out->seg[*it >> 1].soff;
                for (auto it = cc.v.begin(); it != cc.v.end() && xc < 0; ++it)
                    if (out->seg[*it >> 1].rank == 0) xc = out->seg[*it >> 1].soff;
                if (xa >= 0 && xc >= 0 && (chain_strand(ca) == 0 ? xc < xa : xc > xa)) continue;
                best_kept = C.qe - A.qs; bi = x; bk = z;
            }
        if (bi < 0) return nullptr;
        return std::make_unique<GapPair>(GapPair{make_chain(G, G.gc[bi]), make_chain(G, G.gc[bk])});
    };
    auto find_triple = [&](const mg_gchains_t& G, long long best_span) -> std::unique_ptr<Candidate> {
        const int n = std::min<int>(G.n_gc, 20);
        std::vector<int> ord(n);
        for (int i = 0; i < n; ++i) ord[i] = i;
        std::sort(ord.begin(), ord.end(), [&](int x, int y) { return G.gc[x].qs < G.gc[y].qs; });
        std::vector<Chain> ch;
        for (const int i : ord) ch.push_back(make_chain(G, G.gc[i]));
        const long long mv = static_cast<long long>(options.min_var);
        long long best_kept = best_span + mv;
        int bi = -1, bj = -1, bk = -1;
        for (int x = 0; x < n; ++x)
            for (int y = x + 1; y < n; ++y)
                for (int z = y + 1; z < n; ++z) {
                    const Chain &A = ch[x], &B = ch[y], &C = ch[z];
                    const int sa = chain_strand(A), sb = chain_strand(B), sc = chain_strand(C);
                    if (sa != sc || sa == sb) continue;
                    if (B.qs - A.qe > mv || C.qs - B.qe > mv) continue;           // gaps at the turns
                    if (B.qe <= A.qe || C.qe <= B.qe || B.qs <= A.qs || C.qs <= B.qs) continue;
                    if (C.qe - A.qs < best_kept) continue;
                    best_kept = C.qe - A.qs; bi = x; bj = y; bk = z;
                }
        if (bi < 0) return nullptr;
        return std::make_unique<Candidate>(Candidate{std::move(ch[bi]), std::move(ch[bj]), std::move(ch[bk])});
    };
    {
        std::atomic<std::size_t> next{0};
        // One bar for all workers; a haplotype counts once mapped, whatever came of it.
        cli::ProgressBar bar(options.quiet ? std::string() : "[rebuild " + hms() + "] aligning haplotypes",
                             g.paths.size());
        std::mutex bar_mutex;
        auto tick = [&] { const std::lock_guard<std::mutex> lock(bar_mutex); bar.tick(); };
        std::vector<std::thread> pool;
        const unsigned workers = static_cast<unsigned>(std::min<std::size_t>(nthreads, g.paths.size()));
        pool.reserve(workers);
        // A thread that cannot be started throws from emplace_back; the running ones are joined first,
        // since destroying a joinable std::thread calls std::terminate and skips every cleanup.
        try {
            for (unsigned t = 0; t < workers; ++t) {
                pool.emplace_back([&] {
                    mg_tbuf_t* tbuf = mg_tbuf_init();
                    for (std::size_t h = next.fetch_add(1); h < g.paths.size(); h = next.fetch_add(1)) {
                        const std::string& q = hseq[h];
                        if (q.empty()) { tick(); continue; }
                        mg_gchains_t* gc = mg_map(gi, static_cast<int>(q.size()), q.c_str(), tbuf, &opt_map,
                                                  g.paths[h].name.c_str());
                        if (gc != nullptr && gc->n_gc > 0) {
                            // One chain wins (concatenating chains would stitch an incoherent walk). Ranked
                            // by matching bases -- from the base alignment where there is one, else
                            // minigraph's anchor estimate -- then identity within the aligned block, then
                            // mapping quality, then a deterministic tie-break. The outer span would favour a
                            // chain reaching further through a large internal gap.
                            auto mlen_of = [](const mg_gchain_t& c) { return c.p != nullptr ? c.p->mlen : c.mlen; };
                            auto blen_of = [](const mg_gchain_t& c) { return c.p != nullptr ? c.p->blen : c.blen; };
                            auto better = [&](const mg_gchain_t& x, const mg_gchain_t& y) {
                                const auto xm = mlen_of(x), ym = mlen_of(y), xb = blen_of(x), yb = blen_of(y);
                                if (xm != ym) return xm > ym;
                                const double ix = xb > 0 ? static_cast<double>(xm) / xb : 0.0;
                                const double iy = yb > 0 ? static_cast<double>(ym) / yb : 0.0;
                                if (ix != iy) return ix > iy;
                                if (x.mapq != y.mapq) return x.mapq > y.mapq;
                                if (x.score != y.score) return x.score > y.score;
                                return x.qs < y.qs;
                            };
                            int best = 0;
                            for (int i = 1; i < gc->n_gc; ++i) if (better(gc->gc[i], gc->gc[best])) best = i;
                            const Chain c = make_chain(*gc, gc->gc[best]);
                            // Three different measures:
                            //   cover    -- where the alignment begins and ends (the outer envelope)
                            //   matched  -- how much of the query is genuinely aligned (envelope minus gaps)
                            //   chain_id -- how well it matches where it does align
                            // A chain spanning the query through a large internal gap scores high on the
                            // first and low on the second.
                            steps[h] = slice(c, 0, c.plen);
                            chain_span[h] = {c.ps, c.pe, c.qs, c.qe};
                            cover[h] = static_cast<double>(c.qe - c.qs) / static_cast<double>(q.size());
                            matched[h] = c.mlen / static_cast<double>(q.size());
                            chain_id[h] = c.blen > 0 ? c.mlen / c.blen : 0.0;
                            lead_over[h] = std::max<long long>(0, c.ps - c.qs);
                            trail_over[h] = std::max<long long>(
                                0, (c.plen - c.pe) - (static_cast<long long>(q.size()) - c.qe));
                            for (int i = 0; i < gc->n_gc; ++i) chain_q[h].push_back({gc->gc[i].qs, gc->gc[i].qe});
                            if (gc->n_gc >= 3) candidate[h] = find_triple(*gc, c.qe - c.qs);
                            if (!gap_runs(q).empty()) {
                                gap_chain[h] = std::make_unique<Chain>(c);
                                if (gc->n_gc >= 2) gap_pair[h] = find_gap_pair(*gc, q, c.qe - c.qs);
                            }
                        }
                        if (gc != nullptr) mg_gchain_free(gc);
                        tick();
                    }
                    mg_tbuf_destroy(tbuf);
                });
            }
        } catch (...) {
            for (std::thread& th : pool) th.join();
            throw;
        }
        for (std::thread& th : pool) th.join();
    }
    // The seed's walk is read off its own segments, not mapped back: it is the one haplotype the generator
    // reproduces exactly, and each segment records the seed offset it came from (rGFA rank 0). Mapped back,
    // it could route through a near-identical branch added later; with -r the seed is the reference, which
    // must be exact.
    {
        const std::size_t sh = idx.front();
        std::vector<std::pair<std::int32_t, std::uint32_t>> own;      // (offset, segment)
        for (std::uint32_t i = 0; i < out->n_seg; ++i)
            if (out->seg[i].rank == 0) own.push_back({out->seg[i].soff, i});
        std::sort(own.begin(), own.end());
        long long at = 0;
        std::vector<Step> w;
        bool tiles = !own.empty();
        for (const auto& [off, sg] : own) {
            if (off != at) { tiles = false; break; }
            w.push_back({sg << 1, 0, seg_len[sg], false});
            at += seg_len[sg];
        }
        if (tiles && at == static_cast<long long>(hseq[sh].size()) && spell_steps(w) == hseq[sh]) {
            steps[sh] = std::move(w);
            chain_span[sh] = {0, at, 0, at};
            cover[sh] = matched[sh] = chain_id[sh] = 1.0;
            lead_over[sh] = trail_over[sh] = 0;
            candidate[sh].reset();
            gap_chain[sh].reset();
            gap_pair[sh].reset();
        }
    }
    const long long mv = std::max<long long>(1, static_cast<long long>(options.min_var));
    auto wlen_of = [](const std::vector<Step>& w) {
        long long n = 0;
        for (const Step& st : w) n += st.b - st.a;
        return n;
    };
    // ---- assembly gaps: the block of the alignment an N run occupies ----
    // An N run aligns as mismatches, an insertion beside a deletion of the graph sequence opposite it, or
    // both, and an aligner may leave a base or two matching in between. The gap is the whole unmatched block
    // around the run, widened over matching islands shorter than 10 bp beside it; its haplotype bases [q0,
    // q1) replace the walk's [p0, p1), so the walk does not spell both the graph sequence and the N run.
    struct Op { long long q0, q1, p0, p1; int op; };
    auto block_of = [&](const Chain& c, long long g0, long long g1, long long& q0, long long& q1,
                        long long& p0, long long& p1) {
        if (!c.has_cigar) {
            q0 = g0; q1 = g1; p0 = path_at(c, g0); p1 = path_at_right(c, g1);
            return p1 >= p0;
        }
        std::vector<Op> ops;
        long long q = c.qs, pp = c.ps;
        for (const std::uint64_t op_len : c.cigar) {
            const long long len = static_cast<long long>(op_len >> 4);
            const int op = static_cast<int>(op_len & 0xf);
            const bool cq = op == 0 || op == 1 || op == 7 || op == 8, cp = op == 0 || op == 2 || op == 3 || op == 7 || op == 8;
            ops.push_back({q, q + (cq ? len : 0), pp, pp + (cp ? len : 0), op});
            if (cq) q += len;
            if (cp) pp += len;
        }
        auto matching = [](int op) { return op == 0 || op == 7; };
        long long i0 = -1, i1 = -1;
        for (long long i = 0; i < static_cast<long long>(ops.size()); ++i) {
            const Op& o = ops[i];
            const bool touches = o.q1 > o.q0 ? (o.q0 < g1 && o.q1 > g0) : (o.q0 >= g0 && o.q0 <= g1);
            if (!touches) continue;
            if (i0 < 0) i0 = i;
            i1 = i;
        }
        if (i0 < 0) return false;
        auto unmatched = [&](long long i) { return i >= 0 && i < static_cast<long long>(ops.size()) && !matching(ops[i].op); };
        auto island = [&](long long i) {
            return i >= 0 && i < static_cast<long long>(ops.size()) && matching(ops[i].op) && ops[i].q1 - ops[i].q0 < 10;
        };
        for (;;) {
            if (unmatched(i0 - 1)) --i0;
            else if (island(i0 - 1) && unmatched(i0 - 2)) i0 -= 2;
            else break;
        }
        for (;;) {
            if (unmatched(i1 + 1)) ++i1;
            else if (island(i1 + 1) && unmatched(i1 + 2)) i1 += 2;
            else break;
        }
        q0 = ops[i0].q0; p0 = ops[i0].p0; q1 = ops[i1].q1; p1 = ops[i1].p1;
        // A matching op at either end is in the block only for the part inside the run.
        if (matching(ops[i0].op) && q0 < g0) { p0 += g0 - q0; q0 = g0; }
        if (matching(ops[i1].op) && q1 > g1) { p1 -= q1 - g1; q1 = g1; }
        return q1 > q0 && p1 >= p0;
    };
    // The steps for walk range [P0, P1) of chain c, with every N run inside haplotype range [Q0, Q1)
    // replaced by a segment of the haplotype's own bases. `gap_bp` receives the bases so replaced.
    auto gapped_slice = [&](const Chain& c, long long P0, long long P1, long long Q0, long long Q1,
                            std::size_t h, long long& gap_bp) {
        std::vector<Step> w = slice(c, P0, P1);
        gap_bp = 0;
        auto runs = gap_runs(hseq[h]);
        std::reverse(runs.begin(), runs.end());          // right to left, so offsets stay valid
        for (const auto& [g0, g1] : runs) {
            if (g0 < std::max(Q0, static_cast<long long>(c.qs)) || g1 > std::min(Q1, static_cast<long long>(c.qe))) continue;
            long long q0, q1, p0, p1;
            if (!block_of(c, g0, g1, q0, q1, p0, p1) || q0 < Q0 || q1 > Q1 || p0 < P0 || p1 > P1) continue;
            // Settle each edge on the haplotype's own flanking bases: the first 48 bases next to the gap that
            // occur exactly once in the walk near where the alignment put the edge anchor it. An aligner can
            // scatter a long N run's deletion with matching stretches the block rule does not absorb, and the
            // bases right beside a gap are often repetitive.
            {
                constexpr long long K = 48;
                const std::string sp = spell_steps(w);
                const long long reach = std::max<long long>(1000, q1 - q0);
                const long long hl = static_cast<long long>(hseq[h].size());
                auto clean = [&](long long a) {
                    if (a < 0 || a + K > hl) return false;
                    for (long long i = a; i < a + K; ++i) {
                        const char ch = hseq[h][static_cast<std::size_t>(i)];
                        if (ch != 'A' && ch != 'C' && ch != 'G' && ch != 'T') return false;
                    }
                    return true;
                };
                // The single occurrence of hseq[a, a+K) in the walk within `want` +- reach, or -1.
                auto unique_at = [&](long long a, long long want) -> long long {
                    const std::string pat = hseq[h].substr(static_cast<std::size_t>(a), K);
                    const long long lo = std::max<long long>(0, want - reach);
                    long long hit = -1;
                    std::size_t pos = sp.find(pat, static_cast<std::size_t>(lo));
                    while (pos != std::string::npos && static_cast<long long>(pos) <= want + reach) {
                        if (hit >= 0) return -1;
                        hit = static_cast<long long>(pos);
                        pos = sp.find(pat, pos + 1);
                    }
                    return hit;
                };
                for (long long t = 0; t <= 2000; t += 8) {          // right edge: forward from q1
                    if (!clean(q1 + t)) continue;
                    const long long at = unique_at(q1 + t, p1 - P0 + t);
                    if (at >= 0) { p1 = P0 + at - t; break; }
                }
                for (long long t = 0; t <= 2000; t += 8) {          // left edge: back from q0
                    if (!clean(q0 - K - t)) continue;
                    const long long at = unique_at(q0 - K - t, p0 - P0 - K - t);
                    if (at >= 0) { p0 = P0 + at + K + t; break; }
                }
                if (p1 < p0) continue;
            }
            const long long wl = wlen_of(w);
            std::vector<Step> out_w = trim_steps(w, 0, p0 - P0);
            std::vector<Step> rest = trim_steps(w, p1 - P0, wl);
            const std::uint32_t id = add_extra(hseq[h].substr(static_cast<std::size_t>(q0), static_cast<std::size_t>(q1 - q0)), "gap");
            out_w.push_back({id << 1, 0, seg_len[id], true});
            if (!rest.empty()) rest.front().junction = true;
            out_w.insert(out_w.end(), rest.begin(), rest.end());
            w = std::move(out_w);
            gap_bp += q1 - q0;
        }
        return w;
    };
    // A joined walk is judged on what it will spell once its ends are settled: bases past the haplotype's
    // ends are cut later (end cuts), and the haplotype's bases past the walk's become tips.
    auto core_ok = [&](const std::vector<Step>& w, const std::array<long long, 4>& span, std::size_t h) {
        const std::string sp = spell_steps(w);
        const auto [ps, pe, qs, qe] = span;
        const long long qlen = static_cast<long long>(hseq[h].size()), wl = static_cast<long long>(sp.size());
        const long long m1 = std::min(ps, qs), m2 = std::min(wl - pe, qlen - qe);
        if (ps - m1 < 0 || pe + m2 > wl || qs - m1 < 0 || qe + m2 > qlen) return false;
        const std::string wc = sp.substr(static_cast<std::size_t>(ps - m1), static_cast<std::size_t>(pe + m2 - (ps - m1)));
        const std::string hc = hseq[h].substr(static_cast<std::size_t>(qs - m1), static_cast<std::size_t>(qe + m2 - (qs - m1)));
        const long long diff = static_cast<long long>(wc.size()) - static_cast<long long>(hc.size());
        const long long bound = length_bound_for(hseq[h].size());
        if (diff > bound || -diff > bound) return false;
        return identity_of(wc, hc) >= options.min_recovered_identity;
    };
    struct Plan { std::vector<Step> w; std::array<long long, 4> span{}; double cover = 0, matched = 0, cid = 0; int kind = 0; };
    std::vector<std::unique_ptr<Plan>> plan(g.paths.size());
    // Single-chain walks with gaps: replace them in place.
    for (std::size_t h = 0; h < g.paths.size(); ++h) {
        if (!gap_chain[h] || steps[h].empty()) continue;
        const Chain& c = *gap_chain[h];
        long long gbp = 0;
        std::vector<Step> w = gapped_slice(c, 0, c.plen, c.qs, c.qe, h, gbp);
        if (gbp == 0) continue;
        const long long wl = wlen_of(w);
        steps[h] = std::move(w);
        chain_span[h][1] = wl - (c.plen - c.pe);
        matched[h] = std::min(1.0, matched[h] + static_cast<double>(gbp) / static_cast<double>(hseq[h].size()));
    }
    // Turns, chosen in haplotype order so they are shared: a haplotype reuses a turn an earlier one took
    // whenever it lies inside its own overlap, else takes its overlap's midpoint and records it. Chosen per
    // haplotype, the carriers of one inversion would each get their own breakpoints.
    {
        std::vector<GraphPoint> turns_left, turns_right;   // on the first and on the last chain
        auto choose = [&](std::vector<GraphPoint>& seen, const Chain& flank, long long lo, long long hi,
                          long long mid) {
            for (const GraphPoint& gp : seen) {
                const long long pp = from_point(flank, gp, lo, hi);
                if (pp >= 0) return pp;
            }
            const long long pp = path_at(flank, mid);
            seen.push_back(to_point(flank, pp));
            return pp;
        };
        auto matched_part = [](const Chain& c, long long a, long long b, double& ml, double& bl) {
            const double f = c.qe > c.qs ? static_cast<double>(b - a) / (c.qe - c.qs) : 0.0;
            ml += f * c.mlen;
            bl += f * c.blen;
        };
        for (std::size_t h = 0; h < g.paths.size(); ++h) {
            const double qlen = static_cast<double>(hseq[h].size());
            if (candidate[h]) {
                inversion[h] = 1;
                const Chain &A = candidate[h]->a, &B = candidate[h]->b, &C = candidate[h]->c;
                long long xa = A.qe, xb = B.qs, xbe = B.qe, xc = C.qs;      // each side of a small gap
                long long pa = path_at(A, xa), pc = path_at(C, xc);
                if (B.qs < A.qe) {
                    pa = choose(turns_left, A, path_at(A, B.qs), path_at(A, A.qe), (B.qs + A.qe) / 2);
                    xa = xb = std::clamp<long long>(query_at(A, pa), B.qs, A.qe);
                }
                if (C.qs < B.qe) {
                    pc = choose(turns_right, C, path_at(C, C.qs), path_at(C, B.qe), (C.qs + B.qe) / 2);
                    xbe = xc = std::clamp<long long>(query_at(C, pc), C.qs, B.qe);
                }
                if (!(A.qs < xa && xb < xbe && xc < C.qe)) continue;
                // Cut the flanks at the recorded turn itself, so a shared turn is reused exactly.
                long long ga = 0, gb = 0, gcb = 0;
                std::vector<Step> w = gapped_slice(A, 0, pa, A.qs, xa, h, ga);
                std::vector<Step> wb = gapped_slice(B, path_at(B, xb), path_at(B, xbe), xb, xbe, h, gb);
                std::vector<Step> wc = gapped_slice(C, pc, C.plen, xc, C.qe, h, gcb);
                if (w.empty() || wb.empty() || wc.empty()) continue;
                wb.front().junction = true;
                wc.front().junction = true;
                w.insert(w.end(), wb.begin(), wb.end());
                w.insert(w.end(), wc.begin(), wc.end());
                auto pl = std::make_unique<Plan>();
                pl->span = {A.ps, wlen_of(w) - (C.plen - C.pe), A.qs, C.qe};
                double ml = 0, bl = 0;
                matched_part(A, A.qs, xa, ml, bl);
                matched_part(B, xb, xbe, ml, bl);
                matched_part(C, xc, C.qe, ml, bl);
                pl->cover = static_cast<double>(C.qe - A.qs) / qlen;
                pl->matched = std::min(1.0, (ml + ga + gb + gcb) / qlen);
                pl->cid = bl > 0 ? ml / bl : 0.0;
                pl->w = std::move(w);
                pl->kind = 2;
                plan[h] = std::move(pl);
            } else if (gap_pair[h]) {
                // Joined across a gap that split the alignment.
                const Chain &A = gap_pair[h]->a, &C = gap_pair[h]->c;
                long long ga = 0, gcb = 0;
                std::vector<Step> w = gapped_slice(A, 0, path_at(A, A.qe), A.qs, A.qe, h, ga);
                std::vector<Step> wc = gapped_slice(C, path_at(C, C.qs), C.plen, C.qs, C.qe, h, gcb);
                if (w.empty() || wc.empty()) continue;
                const std::string gap_sq = hseq[h].substr(static_cast<std::size_t>(A.qe), static_cast<std::size_t>(C.qs - A.qe));
                const std::uint32_t id = add_extra(gap_sq, "gap");
                w.push_back({id << 1, 0, seg_len[id], true});
                wc.front().junction = true;
                w.insert(w.end(), wc.begin(), wc.end());
                auto pl = std::make_unique<Plan>();
                pl->span = {A.ps, wlen_of(w) - (C.plen - C.pe), A.qs, C.qe};
                pl->cover = static_cast<double>(C.qe - A.qs) / qlen;
                pl->matched = std::min(1.0, (A.mlen + C.mlen + static_cast<double>(gap_sq.size()) + ga + gcb) / qlen);
                pl->cid = (A.blen + C.blen) > 0 ? (A.mlen + C.mlen) / (A.blen + C.blen) : 0.0;
                pl->w = std::move(w);
                pl->kind = 3;
                plan[h] = std::move(pl);
            }
        }
        // Verified in parallel: the joined walk must spell the haplotype within the bounds.
        std::vector<char> verified(g.paths.size(), 0);
        run_parallel(g.paths.size(), nthreads, [&](std::size_t h) {
            if (plan[h]) verified[h] = core_ok(plan[h]->w, plan[h]->span, h);
        });
        for (std::size_t h = 0; h < g.paths.size(); ++h) {
            if (!verified[h]) continue;
            Plan& pl = *plan[h];
            steps[h] = std::move(pl.w);
            chain_span[h] = pl.span;
            cover[h] = pl.cover;
            matched[h] = pl.matched;
            chain_id[h] = pl.cid;
            if (pl.kind == 2) inversion[h] = 2;
        }
    }
    // ---- tips ----
    // A haplotype running past the graph's end (or start) by --min-var or more: the generator never adds
    // sequence that close to a haplotype's end. The walk continues into a segment holding the haplotype's own
    // remaining bases, shared by haplotypes whose remainder is identical. Only at a genuine end: the walk
    // must stop at a segment end with no link onward, and no other chain may map the remainder, or the
    // remainder is graph sequence the alignment did not reach (a rotation, an unjoined inversion, a split
    // alignment). The opposite case, a walk running past the haplotype, is the end cut below.
    {
        auto leaves = [&](std::uint32_t v) {      // does any live arc leave oriented vertex v?
            if ((v >> 1) >= out->n_seg) return true;  // a rebuild segment: not a graph end
            const gfa_arc_t* av = gfa_arc_a(out, v);
            for (std::uint32_t i = 0; i < gfa_arc_n(out, v); ++i) if (!av[i].del) return true;
            return false;
        };
        auto mapped_elsewhere = [&](std::size_t h, long long a, long long b) {
            for (const auto& [cs, ce] : chain_q[h])
                if (std::min(ce, b) - std::max(cs, a) > mv) return true;
            return false;
        };
        // Bases of [a, b) that some chain of h aligns to the graph (the union of their query spans).
        auto mapped_bp = [&](std::size_t h, long long a, long long b) {
            std::vector<std::pair<long long, long long>> iv;
            for (const auto& [cs, ce] : chain_q[h]) {
                const long long lo = std::max(cs, a), hi = std::min(ce, b);
                if (hi > lo) iv.push_back({lo, hi});
            }
            std::sort(iv.begin(), iv.end());
            long long n = 0, end = a;
            for (const auto& [lo, hi] : iv) {
                if (hi <= end) continue;
                n += hi - std::max(lo, end);
                end = hi;
            }
            return n;
        };
        auto ambiguous_in = [](const std::string& sq, std::size_t a, std::size_t b) {
            std::size_t n = 0;
            for (std::size_t i = a; i < b && i < sq.size(); ++i) {
                const char c = sq[i];
                if (c != 'A' && c != 'C' && c != 'G' && c != 'T') ++n;
            }
            return n;
        };
        std::map<std::string, std::uint32_t> tip_left, tip_right;
        auto tip_id = [&](std::map<std::string, std::uint32_t>& seen, const std::string& sq) {
            const auto it = seen.find(sq);
            if (it != seen.end()) return it->second;
            const std::uint32_t id = add_extra(sq, "tip");
            seen.emplace(sq, id);
            return id;
        };
        for (std::size_t h = 0; h < g.paths.size(); ++h) {
            if (steps[h].empty()) continue;
            auto& [ps, pe, qs, qe] = chain_span[h];
            long long qlen = static_cast<long long>(hseq[h].size());
            long long wl = wlen_of(steps[h]);
            if ((qlen - qe) - (wl - pe) >= mv || qs - ps >= mv) {
                // The alignment stops a few bases short of where walk and haplotype stop agreeing: extend it
                // over them first, so a tip holds only what the graph lacks.
                const std::string sp = spell_steps(steps[h]);
                const std::string& hq = hseq[h];
                if ((qlen - qe) - (wl - pe) >= mv)
                    while (pe < wl && qe < qlen && sp[static_cast<std::size_t>(pe)] == hq[static_cast<std::size_t>(qe)]) { ++pe; ++qe; }
                if (qs - ps >= mv)
                    while (ps > 0 && qs > 0 && sp[static_cast<std::size_t>(ps - 1)] == hq[static_cast<std::size_t>(qs - 1)]) { --ps; --qs; }
            }
            // --clip-to-reference: the reference seeds the graph, so its ends are the graph's, and a
            // haplotype running past one by --min-var or more extends outside the locus -- typically
            // extraction carrying on into flanking sequence. The excess is dropped instead of being kept as a
            // tip; the report says how much. Only when it is mostly sequence the graph does not hold: excess
            // that other chains align back into the locus is a rearrangement, not flank.
            if (options.clip_to_reference) {
                const double mbases = matched[h] * static_cast<double>(qlen);
                if ((qlen - qe) - (wl - pe) >= mv && wl - pe < mv && !leaves(steps[h].back().v) &&
                    2 * mapped_bp(h, qe, qlen) < qlen - qe) {
                    steps[h] = trim_steps(steps[h], 0, pe);
                    clipped[h] += qlen - qe;
                    clipped_right[h] += qlen - qe;
                    clipped_gap[h] += ambiguous_in(hseq[h], static_cast<std::size_t>(qe), hseq[h].size());
                    hseq[h].resize(static_cast<std::size_t>(qe));
                }
                if (qs - ps >= mv && ps < mv && !leaves(steps[h].front().v ^ 1) &&
                    2 * mapped_bp(h, 0, qs) < qs) {
                    steps[h] = trim_steps(steps[h], ps, wlen_of(steps[h]));
                    clipped[h] += qs;
                    clipped_left[h] += qs;
                    clipped_gap[h] += ambiguous_in(hseq[h], 0, static_cast<std::size_t>(qs));
                    hseq[h].erase(0, static_cast<std::size_t>(qs));
                    pe -= ps; qe -= qs; ps = 0; qs = 0;
                }
                if (clipped[h] > 0) {
                    qlen = static_cast<long long>(hseq[h].size());
                    wl = wlen_of(steps[h]);
                    matched[h] = qlen > 0 ? std::min(1.0, mbases / static_cast<double>(qlen)) : 0.0;
                    cover[h] = qlen > 0 ? static_cast<double>(qe - qs) / static_cast<double>(qlen) : 0.0;
                    ++sum.clipped_haplotypes;
                    sum.clipped_bp += static_cast<std::size_t>(clipped[h]);
                }
            }
            bool right_tip = (qlen - qe) - (wl - pe) >= mv, left_tip = qs - ps >= mv;
            right_tip = right_tip && (qlen - qe) - (wl - pe) >= mv && wl - pe < mv &&
                        !leaves(steps[h].back().v) && !mapped_elsewhere(h, qe, qlen);
            left_tip = left_tip && qs - ps >= mv && ps < mv && !leaves(steps[h].front().v ^ 1) &&
                       !mapped_elsewhere(h, 0, qs);
            if (right_tip) {
                std::vector<Step> w = trim_steps(steps[h], 0, pe);
                const std::string tail = hseq[h].substr(static_cast<std::size_t>(qe));
                const std::uint32_t id = tip_id(tip_right, tail);
                w.push_back({id << 1, 0, seg_len[id], true});
                steps[h] = std::move(w);
                matched[h] += static_cast<double>(tail.size()) / static_cast<double>(qlen);
                cover[h] += static_cast<double>(tail.size()) / static_cast<double>(qlen);
            }
            if (left_tip) {
                std::vector<Step> w = trim_steps(steps[h], ps, wlen_of(steps[h]));
                if (!w.empty()) w.front().junction = true;
                const std::string head = hseq[h].substr(0, static_cast<std::size_t>(qs));
                const std::uint32_t id = tip_id(tip_left, head);
                w.insert(w.begin(), Step{id << 1, 0, seg_len[id], false});
                steps[h] = std::move(w);
                pe += qs - ps;
                ps = qs;
                matched[h] += static_cast<double>(head.size()) / static_cast<double>(qlen);
                cover[h] += static_cast<double>(head.size()) / static_cast<double>(qlen);
            }
            matched[h] = std::min(1.0, matched[h]);
            cover[h] = std::min(1.0, cover[h]);
        }
        // Overruns from the walks as they now stand.
        for (std::size_t h = 0; h < g.paths.size(); ++h) {
            if (steps[h].empty()) continue;
            const auto [ps, pe, qs, qe] = chain_span[h];
            const long long qlen = static_cast<long long>(hseq[h].size());
            lead_over[h] = std::max<long long>(0, ps - qs);
            trail_over[h] = std::max<long long>(0, (wlen_of(steps[h]) - pe) - (qlen - qe));
        }
    }
    // A walk is cut where its haplotype stops, but only where the haplotype genuinely stops there: its bases
    // outside the alignment must match the walk up to the cut. When they do not -- a haplotype lacking an
    // insertion the seed carries near the end -- cutting at the right length would end the walk in the seed's
    // sequence; those are left uncut and the overrun bound refuses them. Unless the haplotype's own unaligned end
    // is shorter than --min-var: then the walk is cut where the alignment ends, a difference under --min-var like
    // any other the rebuild collapses. The overrun may run over several nodes (a short haplotype ending early in
    // a long seed's flank): the walk is trimmed over as many as it takes, cutting inside the node where the
    // haplotype stops.
    auto tail_matches = [&](const std::string& walk_part, const std::string& hap_part) {
        if (hap_part.empty()) return true;
        if (walk_part.size() != hap_part.size()) return false;
        std::size_t ed = 0;
        const double id = identity_of(walk_part, hap_part, &ed);
        return ed <= 2 || id >= options.min_recovered_identity;
    };
    auto cut_ends = [&](std::size_t h) {
        if (steps[h].empty() || (lead_over[h] <= 0 && trail_over[h] <= 0)) return;
        const std::string walk_sp = spell_steps(steps[h]);
        auto& [ps, pe, qs, qe] = chain_span[h];
        const long long qlen = static_cast<long long>(hseq[h].size());
        const long long mvc = std::max<long long>(1, static_cast<long long>(options.min_var));
        const long long wl = static_cast<long long>(walk_sp.size());
        long long lead_cut = 0, trail_cut = 0;
        if (lead_over[h] > 0 && lead_over[h] < wl &&
            tail_matches(walk_sp.substr(static_cast<std::size_t>(ps - qs), static_cast<std::size_t>(qs)),
                         hseq[h].substr(0, static_cast<std::size_t>(qs)))) {
            lead_cut = lead_over[h];
        } else if (lead_over[h] > 0 && qs < mvc && ps < wl) {
            lead_cut = ps;
        }
        if (trail_over[h] > 0 && trail_over[h] < wl - lead_cut &&
            tail_matches(walk_sp.substr(static_cast<std::size_t>(pe), static_cast<std::size_t>(qlen - qe)),
                         hseq[h].substr(static_cast<std::size_t>(qe)))) {
            trail_cut = trail_over[h];
        } else if (trail_over[h] > 0 && qlen - qe < mvc && wl - pe < wl - lead_cut) {
            trail_cut = wl - pe;
        }
        if (lead_cut == 0 && trail_cut == 0) return;
        std::vector<Step> kept = trim_steps(steps[h], lead_cut, wl - trail_cut);
        if (kept.empty()) return;
        steps[h] = std::move(kept);
        ps -= lead_cut;
        pe -= lead_cut;
        if (lead_cut > 0) { lead_over[h] = 0; ++sum.end_cuts; }
        if (trail_cut > 0) { trail_over[h] = 0; ++sum.end_cuts; }
    };
    // Cut before the repair below, so a walk running on past its haplotype's end is not taken for a deletion.
    for (std::size_t h = 0; h < g.paths.size(); ++h) cut_ends(h);
    // Haplotypes with a copy gain or loss left for the fold; the second pass writes those the fold did not take.
    std::vector<char> held_back(g.paths.size(), 0);
    // Modules found from the haplotypes' own sequences when a repair needs them, then handed to the fold.
    std::optional<std::vector<ModuleUnit>> known_modules;
    // ---- indels the generator left out ----
    // The generator drops a difference near either end of an alignment chain or of a haplotype, so a haplotype
    // carrying a deletion there comes back through the other allele, its walk too long by the deletion (and one
    // carrying an insertion, too short). A walk whose length differs from its haplotype's by more than the bound
    // (whole, or over its aligned part), or whose alignment leaves the walk and the haplotype --min-var or more
    // apart in length at an end, is aligned to the haplotype with gap-aware scoring (minimap2), so a deletion is
    // one gap. Gaps no more than --min-var matching bases apart are one block, since differences shorter than
    // that collapse anyway; a block whose net length change reaches --min-var, less the bases both sides share at
    // its ends, replaces the walk's part with the haplotype's: a deletion is an edge past the skipped bases, an
    // insertion a node of the haplotype's own bases, shared by haplotypes inserting the same bases at the same
    // place. Past the alignment's ends, an indel leaves too few bases beyond it to align: there the shorter
    // side, --min-var or more, must match the longer side's end, and the rest of the longer side is the indel.
    // The repaired walk is kept only if it then passes the bounds. With --mode folded, a block holding half a
    // module's unit or more is a copy gained or lost, left for the fold.
    {
        // Walk [p0, p1) becomes haplotype [q0, q1); `end` 1 or 2 for a block before or after the alignment.
        struct Indel { long long p0, p1, q0, q1, gain; int end = 0; };
        std::vector<std::vector<Indel>> found(g.paths.size());
        std::vector<std::array<long long, 4>> span_of(g.paths.size(), {0, 0, 0, 0});   // the alignment's extent
        std::vector<long long> aligned_eq(g.paths.size(), 0);   // its matching bases
        const long long mvl = std::max<long long>(1, static_cast<long long>(options.min_var));
        run_parallel(g.paths.size(), nthreads, [&](std::size_t h) {
            if (steps[h].empty() || inversion[h] == 2) return;
            const auto [ps, pe, qs, qe] = chain_span[h];
            const std::string sp = spell_steps(steps[h]);
            const long long bound = length_bound_for(hseq[h].size());
            const bool core_off = pe > ps && qe > qs && std::llabs((pe - ps) - (qe - qs)) > bound;
            const long long wl = static_cast<long long>(sp.size()), ql = static_cast<long long>(hseq[h].size());
            const bool whole_off = std::llabs(wl - ql) > bound;
            const bool end_off = std::llabs(ps - qs) >= mvl || std::llabs((wl - pe) - (ql - qe)) >= mvl;
            if (!core_off && !whole_off && !end_off) return;
            // The haplotype's alignment to its walk, as the forward pieces that follow each other in both: an
            // assembly gap's N bases match nothing, so a long one splits the alignment, and what lies between two
            // pieces is compared like a gap inside one.
            std::vector<Minimap2Hit> all = minimap2_hits("haplotype", hseq[h], "walk", sp, "asm20", 16);
            std::sort(all.begin(), all.end(), [](const Minimap2Hit& x, const Minimap2Hit& y) {
                return x.query_end_bp - x.query_start_bp > y.query_end_bp - y.query_start_bp;
            });
            std::vector<const Minimap2Hit*> chain;
            for (const Minimap2Hit& c : all) {
                if (c.reverse || c.cigar.empty()) continue;
                const bool fits = std::all_of(chain.begin(), chain.end(), [&](const Minimap2Hit* x) {
                    const bool before = c.query_end_bp <= x->query_start_bp && c.target_end_bp <= x->target_start_bp;
                    const bool after = c.query_start_bp >= x->query_end_bp && c.target_start_bp >= x->target_end_bp;
                    return before || after;
                });
                if (fits) chain.push_back(&c);
            }
            if (chain.empty()) return;
            std::sort(chain.begin(), chain.end(), [](const Minimap2Hit* x, const Minimap2Hit* y) {
                return x->query_start_bp < y->query_start_bp;
            });
            span_of[h] = {static_cast<long long>(chain.front()->target_start_bp), static_cast<long long>(chain.back()->target_end_bp),
                          static_cast<long long>(chain.front()->query_start_bp), static_cast<long long>(chain.back()->query_end_bp)};
            // Outside the alignment: walk and haplotype bases before its start (lead) or after its end.
            auto end_block = [&](bool lead) {
                const long long t = lead ? span_of[h][0] : span_of[h][1], q = lead ? span_of[h][2] : span_of[h][3];
                const long long nw = lead ? t : wl - t, nq = lead ? q : ql - q, outer = std::min(nw, nq);
                if (outer < mvl || std::llabs(nw - nq) < mvl) return;
                const std::string w_out = lead ? sp.substr(0, static_cast<std::size_t>(outer)) : sp.substr(static_cast<std::size_t>(wl - outer));
                const std::string q_out = lead ? hseq[h].substr(0, static_cast<std::size_t>(outer))
                                               : hseq[h].substr(static_cast<std::size_t>(ql - outer));
                if (!tail_matches(w_out, q_out)) return;
                // Gain: the inserted bases and the outer ones, none of which the alignment counted.
                if (lead) found[h].push_back({outer, t, outer, q, q, 1});
                else found[h].push_back({t, wl - outer, q, ql - outer, ql - q, 2});
            };
            end_block(true);
            long long t0 = -1, q0 = 0, t1 = 0, q1 = 0, eq_in = 0, eq_pending = 0;
            auto close_block = [&]() {
                if (t0 < 0) return;
                long long a = t0, b = t1, c = q0, d = q1, trimmed = 0;
                while (a < b && c < d && sp[static_cast<std::size_t>(a)] == hseq[h][static_cast<std::size_t>(c)]) { ++a; ++c; ++trimmed; }
                while (a < b && c < d && sp[static_cast<std::size_t>(b - 1)] == hseq[h][static_cast<std::size_t>(d - 1)]) { --b; --d; ++trimmed; }
                if (std::llabs((d - c) - (b - a)) >= mvl)
                    found[h].push_back({a, b, c, d, std::max<long long>(0, (d - c) - std::max<long long>(0, eq_in - trimmed))});
                t0 = -1;
                eq_in = eq_pending = 0;
            };
            for (std::size_t k = 0; k < chain.size(); ++k) {
                const Minimap2Hit& hit = *chain[k];
                if (k > 0) {   // between two pieces
                    t0 = static_cast<long long>(chain[k - 1]->target_end_bp); q0 = static_cast<long long>(chain[k - 1]->query_end_bp);
                    t1 = static_cast<long long>(hit.target_start_bp); q1 = static_cast<long long>(hit.query_start_bp);
                    close_block();
                }
                long long t = static_cast<long long>(hit.target_start_bp), q = static_cast<long long>(hit.query_start_bp);
                for (const auto& [op, len] : hit.cigar) if (op == 7) aligned_eq[h] += len;
                for (const auto& [op, len] : hit.cigar) {
                    const long long n = len;
                    if (op == 7) {                          // '=': matching bases
                        if (t0 >= 0) { if (n > mvl) close_block(); else eq_pending += n; }
                        t += n; q += n;
                        continue;
                    }
                    if (op == 8 && t0 < 0) { t += n; q += n; continue; }   // a mismatch alone is no indel
                    if (t0 < 0) { t0 = t; q0 = q; }
                    eq_in += eq_pending; eq_pending = 0;
                    if (op != 1) t += n;                    // 'X' and 'D' (walk only) step the walk
                    if (op != 2) q += n;                    // 'X' and 'I' (haplotype only) step the haplotype
                    t1 = t; q1 = q;
                }
                close_block();
            }
            end_block(false);
        });
        // Copy changes are left for the fold: a block's longer side holding half a module's unit or more at the
        // copy identity floor. The modules are found once, from the haplotypes' own sequences.
        const std::unordered_set<std::string> copies_written(options.copies_for.begin(), options.copies_for.end());
        bool any = false;
        for (const auto& f : found) any = any || !f.empty();
        if (any && options.fold_modules) {
            std::vector<std::string> names(g.paths.size());
            for (std::size_t h = 0; h < g.paths.size(); ++h) names[h] = g.paths[h].name;
            ModuleFoldOptions mfo;
            mfo.threads = nthreads;
            mfo.min_identity = options.fold_identity;
            mfo.min_copy_identity = std::min(0.90, options.fold_identity);
            mfo.min_copy_bp = options.fold_min_copy_bp;
            known_modules = discover_modules(names, hseq, idx.front(), mfo);
            run_parallel(g.paths.size(), nthreads, [&](std::size_t h) {
                if (found[h].empty() || copies_written.count(g.paths[h].name)) return;
                const std::string sp = spell_steps(steps[h]);
                std::vector<Indel> kept;
                for (const Indel& in : found[h]) {
                    const bool walk_side = in.p1 - in.p0 >= in.q1 - in.q0;
                    const std::string side = walk_side ? sp.substr(static_cast<std::size_t>(in.p0), static_cast<std::size_t>(in.p1 - in.p0))
                                                       : hseq[h].substr(static_cast<std::size_t>(in.q0), static_cast<std::size_t>(in.q1 - in.q0));
                    bool copy = false;
                    for (const ModuleUnit& m : *known_modules) {
                        if (2 * side.size() < m.unit.size()) continue;
                        for (const Minimap2Hit& u : minimap2_hits("block", side, "unit", m.unit, "asm20", 4))
                            if (!u.reverse && 2 * (u.target_end_bp - u.target_start_bp) >= m.unit.size() &&
                                u.gc_identity >= mfo.min_copy_identity) { copy = true; break; }
                        if (copy) break;
                    }
                    if (copy) held_back[h] = 1;
                    else kept.push_back(in);
                }
                found[h] = std::move(kept);
            });
        }
        // Insertion nodes shared by haplotypes inserting the same bases at the same point of the graph.
        std::map<std::string, std::uint32_t> insertion_ids;
        auto point_of = [&](const std::vector<Step>& w, long long off) {
            long long at = 0;
            for (const Step& st : w) {
                const long long len = st.b - st.a;
                if (off <= at + len) {
                    const long long in = off - at;
                    return std::to_string(st.v >> 1) + ":" + std::to_string((st.v & 1) ? st.b - in : st.a + in) +
                           ((st.v & 1) ? "-" : "+");
                }
                at += len;
            }
            return std::string("end");
        };
        for (std::size_t h = 0; h < g.paths.size(); ++h) {
            if (found[h].empty()) continue;
            std::vector<Step> w = steps[h];
            long long gain = 0, delta = 0;
            bool lead_fixed = false, trail_fixed = false;
            for (auto it = found[h].rbegin(); it != found[h].rend(); ++it) {   // right to left: offsets stay valid
                const long long wl = wlen_of(w);
                const std::string at = point_of(w, it->p0) + "|" + point_of(w, it->p1);
                std::vector<Step> left = trim_steps(w, 0, it->p0), right = trim_steps(w, it->p1, wl);
                if (it->q1 > it->q0) {
                    std::string ins = hseq[h].substr(static_cast<std::size_t>(it->q0), static_cast<std::size_t>(it->q1 - it->q0));
                    const std::string key = at + "|" + ins;
                    auto found_id = insertion_ids.find(key);
                    if (found_id == insertion_ids.end())
                        found_id = insertion_ids.emplace(key, add_extra(std::move(ins), "insertion")).first;
                    left.push_back({found_id->second << 1, 0, seg_len[found_id->second], true});
                }
                if (!right.empty()) right.front().junction = true;
                left.insert(left.end(), right.begin(), right.end());
                w = std::move(left);
                gain += it->gain;
                lead_fixed = lead_fixed || it->end == 1;
                trail_fixed = trail_fixed || it->end == 2;
                if (it->end != 2) delta += (it->q1 - it->q0) - (it->p1 - it->p0);   // a trailing block lies past span[1]
            }
            auto span = span_of[h];
            span[1] += delta;
            // A repaired end block makes the walk agree with the haplotype out to that end.
            if (lead_fixed) span[0] = span[2] = 0;
            if (trail_fixed) { span[1] = wlen_of(w); span[3] = static_cast<long long>(hseq[h].size()); }
            if (!core_ok(w, span, h)) continue;
            steps[h] = std::move(w);
            chain_span[h] = span;
            // Matched cover from this alignment, base by base: its matching bases and those the repair makes match.
            matched[h] = std::max(matched[h], std::min(1.0, static_cast<double>(aligned_eq[h] + gain) /
                                                                static_cast<double>(hseq[h].size())));
            const long long qlen = static_cast<long long>(hseq[h].size());
            lead_over[h] = std::max<long long>(0, span[0] - span[2]);
            trail_over[h] = std::max<long long>(0, (wlen_of(steps[h]) - span[1]) - (qlen - span[3]));
            sum.indels_added += found[h].size();
        }
        step("indels the generator left out, written into walks: " + std::to_string(sum.indels_added));
    }
    mg_idx_destroy(gi);
    // No variation at all, nor a copy change left for the fold: every haplotype would spell the seed. Keep the
    // original.
    if (generated_flat && extra_seq.empty() && sum.indels_added == 0 &&
        std::find(held_back.begin(), held_back.end(), 1) == held_back.end()) {
        step("the rebuilt graph holds no variation: passing the input through unchanged");
        if (!options.out_path.empty()) copy_file(options.gfa_path, options.out_path, out_staged);
        drop_stale_report();
        sum.out_nodes = sum.raw_nodes;
        sum.ran = false;
        sum.degenerate = true;
        return sum;
    }
    double cover_sum = 0.0, matched_sum = 0.0, chainid_sum = 0.0;
    double matched_min = 1.0;
    std::size_t n_recovered = 0;
    for (std::size_t h = 0; h < steps.size(); ++h) {
        if (!steps[h].empty()) {
            ++sum.paths_recovered;
            ++n_recovered;
            chainid_sum += chain_id[h];
        }
        cover_sum += cover[h];
        matched_sum += matched[h];
        matched_min = std::min(matched_min, matched[h]);
    }
    sum.mean_query_cover =
        g.paths.empty() ? 0.0 : cover_sum / static_cast<double>(g.paths.size());
    sum.mean_matched_cover =
        g.paths.empty() ? 0.0 : matched_sum / static_cast<double>(g.paths.size());
    sum.min_matched_cover = g.paths.empty() ? 0.0 : matched_min;
    sum.mean_chain_identity = n_recovered ? chainid_sum / static_cast<double>(n_recovered) : 0.0;
    step("recovered walks: " + std::to_string(sum.paths_recovered) + "/" + std::to_string(g.paths.size()) +
         "; matched cover mean: " + fmt2(sum.mean_matched_cover) + "; min: " + fmt2(sum.min_matched_cover) +
         " (threshold: " + fmt_short(options.min_matched_cover) + ")");

    // ---- the emitted graph: minigraph's segments, cut where haplotypes end, minus what no walk uses ----
    // A walk lists whole segments, and the generator never cuts one where a haplotype ends, so a haplotype
    // stopping partway into its first or last segment would be given the rest of it. The chain says where the
    // haplotype starts and ends on its walk, so the end segment is cut there. Pieces and links no walk uses
    // are dropped: a haplotype is mapped back against the finished graph and can take another route than the
    // branch its generation-time alignment added.
    struct Piece { std::uint32_t seg, beg, end; };   // [beg, end) of minigraph segment `seg`
    const auto n_all = static_cast<std::uint32_t>(seg_len.size());   // generator's segments + rebuild's
    std::vector<std::vector<std::uint32_t>> cuts(n_all);
    // Walks the repair changed are cut again (see cut_ends above).
    for (std::size_t h = 0; h < g.paths.size(); ++h) {
        if (steps[h].empty()) continue;
        cut_ends(h);
        for (const Step& st : steps[h]) {
            cuts[st.v >> 1].push_back(st.a);
            cuts[st.v >> 1].push_back(st.b);
        }
        if (inversion[h] == 2) ++sum.inversions_joined;
    }
    {
        std::size_t unconfirmed = 0;
        for (const int iv : inversion) unconfirmed += iv == 1;
        step("large inversions joined: " + std::to_string(sum.inversions_joined) +
                 "; candidates their sequence does not confirm: " + std::to_string(unconfirmed));
    }
    std::vector<Piece> pieces;
    std::vector<std::uint32_t> first_piece(n_all + 1, 0);
    for (std::uint32_t i = 0; i < n_all; ++i) {
        first_piece[i] = static_cast<std::uint32_t>(pieces.size());
        std::vector<std::uint32_t>& c = cuts[i];
        c.push_back(0);
        c.push_back(seg_len[i]);
        std::sort(c.begin(), c.end());
        c.erase(std::unique(c.begin(), c.end()), c.end());
        for (std::size_t k = 0; k + 1 < c.size(); ++k) pieces.push_back({i, c[k], c[k + 1]});
        if (c.size() == 1) pieces.push_back({i, 0, 0});            // empty segment: keep it addressable
    }
    first_piece[n_all] = static_cast<std::uint32_t>(pieces.size());
    auto first_of = [&](std::uint32_t seg) { return first_piece[seg]; };
    auto last_of = [&](std::uint32_t seg) { return first_piece[seg + 1] - 1; };
    // Walks re-expressed over pieces: oriented ids (piece << 1 | reverse). The turns of a joined
    // inversion are links the generator never made; they are collected to be added.
    std::vector<std::vector<std::uint32_t>> walks(g.paths.size());
    std::vector<std::pair<std::uint32_t, std::uint32_t>> junctions;
    for (std::size_t h = 0; h < g.paths.size(); ++h) {
        if (steps[h].empty()) continue;
        std::vector<std::uint32_t> pw;
        bool pending_junction = false;
        for (const Step& st : steps[h]) {
            const std::uint32_t sg = st.v >> 1, rev = st.v & 1;
            std::vector<std::uint32_t> ids;
            for (std::uint32_t pid = first_of(sg); pid <= last_of(sg); ++pid)
                if (pieces[pid].beg >= st.a && pieces[pid].end <= st.b && pieces[pid].end > pieces[pid].beg)
                    ids.push_back(pid);
            if (rev) std::reverse(ids.begin(), ids.end());
            pending_junction = pending_junction || st.junction;
            if (ids.empty()) continue;
            if (pending_junction && !pw.empty()) junctions.push_back({pw.back(), ids.front() << 1 | rev});
            pending_junction = false;
            for (const std::uint32_t pid : ids) pw.push_back(pid << 1 | rev);
        }
        walks[h] = std::move(pw);
    }
    // Every link, both orientations, over pieces: minigraph's arcs re-attached to the end pieces, plus
    // the joins between consecutive pieces of a cut segment. For v -> w the dual is (w^1) -> (v^1).
    auto arc_key = [](std::uint32_t x, std::uint32_t y) { return (static_cast<std::uint64_t>(x) << 32) | y; };
    std::unordered_set<std::uint64_t> arcs;
    arcs.reserve(static_cast<std::size_t>(out->n_arc) * 4 + pieces.size() * 2);
    auto add_arc = [&](std::uint32_t x, std::uint32_t y) {
        arcs.insert(arc_key(x, y));
        arcs.insert(arc_key(y ^ 1, x ^ 1));
    };
    for (std::uint64_t k = 0; k < out->n_arc; ++k) {
        const gfa_arc_t& a = out->arc[k];
        if (a.del) continue;
        const std::uint32_t v = static_cast<std::uint32_t>(a.v_lv >> 32), w = a.w;
        const std::uint32_t from = (v & 1) ? (first_of(v >> 1) << 1 | 1) : (last_of(v >> 1) << 1);
        const std::uint32_t to = (w & 1) ? (last_of(w >> 1) << 1 | 1) : (first_of(w >> 1) << 1);
        add_arc(from, to);
    }
    for (std::uint32_t i = 0; i < n_all; ++i)
        for (std::uint32_t pid = first_of(i); pid < last_of(i); ++pid) add_arc(pid << 1, (pid + 1) << 1);
    for (const auto& [x, y] : junctions) add_arc(x, y);
    // What the walks use, canonical per link so each is emitted once.
    std::vector<char> piece_used(pieces.size(), 0);
    std::unordered_set<std::uint64_t> used_arcs;
    for (const auto& w : walks) {
        for (std::size_t j = 0; j < w.size(); ++j) {
            piece_used[w[j] >> 1] = 1;
            if (j == 0) continue;
            const std::uint64_t k1 = arc_key(w[j - 1], w[j]), k2 = arc_key(w[j] ^ 1, w[j - 1] ^ 1);
            if (arcs.count(k1)) used_arcs.insert(std::min(k1, k2));
        }
    }
    std::size_t all_links = 0;
    for (const std::uint64_t k : arcs) {
        const std::uint64_t dual = arc_key(static_cast<std::uint32_t>(k) ^ 1, static_cast<std::uint32_t>(k >> 32) ^ 1);
        if (k <= dual) ++all_links;
    }
    std::vector<std::uint32_t> out_id(pieces.size(), 0);   // 1-based GFA id of each kept piece
    std::uint32_t n_kept = 0;
    for (std::size_t pid = 0; pid < pieces.size(); ++pid)
        if (piece_used[pid]) out_id[pid] = ++n_kept;
    for (std::size_t pid = 0; pid < pieces.size(); ++pid) {
        const std::uint32_t sg = pieces[pid].seg;
        if (sg < out->n_seg) { if (!piece_used[pid]) ++sum.pruned_segments; continue; }
        if (!piece_used[pid] || pieces[pid].beg != 0) continue;   // count each used rebuild segment once
        if (extra_kind[sg - out->n_seg] == "gap") ++sum.gaps_kept;
        else if (extra_kind[sg - out->n_seg] == "tip") ++sum.tips_added;
    }
    sum.pruned_links = all_links - used_arcs.size();
    // A clipped haplotype named with a PanSN interval (name:start-end) is emitted with that interval narrowed
    // by what was clipped from each end. Rewritten only when the interval's length is the haplotype's own, so
    // a name using another convention is left alone.
    auto emitted_name = [&](std::size_t h) -> std::string {
        const std::string& name = g.paths[h].name;
        if (clipped[h] == 0) return name;
        const std::size_t colon = name.rfind(':');
        const std::size_t dash = name.rfind('-');
        if (colon == std::string::npos || dash == std::string::npos || dash < colon + 2 || dash + 1 >= name.size())
            return name;
        const std::string a = name.substr(colon + 1, dash - colon - 1), b = name.substr(dash + 1);
        const auto digits = [](const std::string& t) {
            return !t.empty() && t.size() < 19 &&
                   std::all_of(t.begin(), t.end(), [](char c) { return c >= '0' && c <= '9'; });
        };
        if (!digits(a) || !digits(b)) return name;
        const long long start = std::stoll(a), end = std::stoll(b);
        const long long given = static_cast<long long>(hseq[h].size()) + clipped[h];
        // Either convention: 1-based inclusive (end - start + 1 bases) or half-open (end - start bases); both
        // narrow the same way.
        if (!((start >= 1 && end - start + 1 == given) || (start >= 0 && end - start == given))) return name;
        return name.substr(0, colon + 1) + std::to_string(start + clipped_left[h]) + "-" +
               std::to_string(end - clipped_right[h]);
    };
    std::vector<std::uint64_t> emit_arcs(used_arcs.begin(), used_arcs.end());
    std::sort(emit_arcs.begin(), emit_arcs.end());
    step("clipped past the reference's ends (--clip-to-reference): haplotypes: " +
             std::to_string(sum.clipped_haplotypes) + "; bp: " + std::to_string(sum.clipped_bp));
    step("assembly gaps kept as N nodes: " + std::to_string(sum.gaps_kept) +
             "; tip nodes added for haplotypes running past the graph: " + std::to_string(sum.tips_added));
    step("node ends cut where a haplotype stops: " + std::to_string(sum.end_cuts) +
             "; unused nodes dropped: " + std::to_string(sum.pruned_segments) +
             "; unused edges dropped: " + std::to_string(sum.pruned_links));
    auto piece_seq = [&](std::uint32_t pid) {
        const Piece& pc = pieces[pid];
        const char* sq = seg_seq(pc.seg);
        return sq ? std::string(sq + pc.beg, pc.end - pc.beg) : std::string();
    };

    // ---- emit plain GFA: S + oriented L + one P per haplotype ----
    // Orientation is carried through: flattening links to +/+ would destroy every inversion downstream.
    // Segments are renamed to consecutive integers, since tools that store node ids as integers reject
    // minigraph's "s1", "s2", ...; the original name is kept in a user-defined `mn` tag, with the piece's
    // range when the segment was cut (rGFA's SN names a stable source sequence, which these are not).
    StagedFile staged_guard(out_staged);   // removed on any exit that does not commit it
    std::filesystem::path staged_out;
    std::string emitted_gfa;   // kept for the module fold (--mode folded)
    std::optional<std::array<std::size_t, 5>> fold_summary;   // nodes, links, maxdeg, hubs, selfloops after folding
    std::size_t folded_bp = 0;                                // the folded graph's sequence
    if (!options.out_path.empty()) {
        const std::filesystem::path& staged = out_staged;
        std::ostringstream o;
        o << "H\tVN:Z:1.0\n";
        for (std::uint32_t pid = 0; pid < pieces.size(); ++pid) {
            if (!piece_used[pid]) continue;
            const Piece& pc = pieces[pid];
            const std::string sq = piece_seq(pid);
            o << "S\t" << out_id[pid] << '\t' << (sq.empty() ? "*" : sq);
            const char* name = pc.seg < out->n_seg ? out->seg[pc.seg].name : extra_kind[pc.seg - out->n_seg].c_str();
            if (name != nullptr) {
                o << "\tmn:Z:" << name;
                if (pc.beg != 0 || pc.end != seg_len[pc.seg]) o << ':' << pc.beg << '-' << pc.end;
            }
            o << '\n';
        }
        for (const std::uint64_t k : emit_arcs) {
            const std::uint32_t x = static_cast<std::uint32_t>(k >> 32), y = static_cast<std::uint32_t>(k);
            o << "L\t" << out_id[x >> 1] << '\t' << "+-"[x & 1] << '\t' << out_id[y >> 1] << '\t'
              << "+-"[y & 1] << "\t0M\n";
        }
        for (std::size_t h = 0; h < g.paths.size(); ++h) {
            if (walks[h].empty()) continue;
            o << "P\t" << emitted_name(h) << '\t';
            for (std::size_t j = 0; j < walks[h].size(); ++j) {
                if (j) o << ',';
                o << out_id[walks[h][j] >> 1] << ((walks[h][j] & 1) ? '-' : '+');
            }
            o << "\t*\n";
        }
        if (options.fold_modules) emitted_gfa = o.str();
        write_text(staged, options.fold_modules ? emitted_gfa : o.str(), wants_gzip(options.out_path));
        staged_out = staged;   // committed only if the acceptance contract below is satisfied
    }

    std::vector<double> walk_id;
    std::vector<char> walk_exact(g.paths.size(), 0);      // walk spells the haplotype byte for byte
    std::vector<std::size_t> walk_edits(g.paths.size(), 0), walk_aln(g.paths.size(), 0);
    std::vector<long long> len_diff(g.paths.size(), 0);   // spelled walk length - haplotype length
    std::vector<std::size_t> walk_ambiguous(g.paths.size(), 0);   // N (non-ACGT) bases the walk spells
    // ---- recovered-walk identity: the end-to-end check ----
    // Re-spell each recovered walk from the rebuilt graph and compare it with the haplotype the input
    // spelled. Chain cover and identity describe the alignment; this describes what a caller gets back, so it
    // catches an error anywhere in mapping, chain choice or emission.
    {
        std::vector<std::string> seg_seq(pieces.size());
        for (std::uint32_t pid = 0; pid < pieces.size(); ++pid) seg_seq[pid] = piece_seq(pid);
        walk_id.assign(g.paths.size(), -1.0);
        std::vector<double>& wid = walk_id;
        run_parallel(g.paths.size(), nthreads, [&](std::size_t h) {
            if (walks[h].empty() || hseq[h].empty()) return;
            std::string spelled;
            for (const std::uint32_t v : walks[h]) {
                const std::uint32_t seg = v >> 1;
                if (seg >= seg_seq.size()) return;
                const std::string& sq = seg_seq[seg];
                spelled += (v & 1) ? reverse_complement(sq) : sq;
            }
            if (spelled.empty()) return;
            walk_exact[h] = spelled == hseq[h];
            for (const char c : spelled)
                if (c != 'A' && c != 'C' && c != 'G' && c != 'T') ++walk_ambiguous[h];
            len_diff[h] = static_cast<long long>(spelled.size()) - static_cast<long long>(hseq[h].size());
            wid[h] = identity_of(spelled, hseq[h], &walk_edits[h]);
            walk_aln[h] = std::max(spelled.size(), hseq[h].size());   // the identity's denominator
        });
        double wsum = 0.0, wmin = 1.0;
        std::size_t wn = 0, wpoor = 0;
        for (const double v : wid) {
            if (v < 0.0) continue;
            wsum += v; wmin = std::min(wmin, v); ++wn;
            if (v < 0.99) ++wpoor;   // the count an acceptance threshold would act on
        }
        sum.walks_below_99 = wpoor;
        sum.walk_identity_checked = wn;
        sum.mean_walk_identity = wn ? wsum / static_cast<double>(wn) : 0.0;
        sum.min_walk_identity = wn ? wmin : 0.0;
        step("walk identity mean: " + fmt2(sum.mean_walk_identity) + "; min: " + fmt2(sum.min_walk_identity) +
             " (threshold: " + fmt_short(options.min_recovered_identity) + "); walks below 0.99: " +
             std::to_string(wpoor));
    }

    // ---- acceptance contract, and the report that explains its verdict ----
    //
    // minigraph augments variation above --min-var, so differences below it are collapsed and a walk need not
    // be byte-identical. The contract is structural plus a threshold: every path comes back, its steps are
    // joined by links that exist, and it spells what it spelled before within the bounds. Otherwise the
    // rebuilt graph is discarded and the original written, since a graph that silently drops or truncates a
    // haplotype is worse than no rebuild.
    {
        // Every consecutive pair of steps on an emitted walk needs a link that actually exists, or the
        // P line describes a traversal the graph does not permit. `arcs` holds both orientations of
        // every link (a walk may traverse one the other way round), built with the pieces above.
        // Assembly gaps: N (or any non-ACGT) bases count against both measures where the walk spells graph
        // bases opposite them. A failing haplotype whose shortfall the gap accounts for is labelled
        // `assembly_gap`; it still fails, since its walk would carry the backbone across the gap, which reads
        // downstream as the reference allele rather than as unknown. Counted over the haplotype as given,
        // clipped N bases included.
        std::vector<std::size_t> gap_bp(g.paths.size(), 0);
        for (std::size_t h = 0; h < g.paths.size(); ++h) {
            for (const char c : hseq[h])
                if (c != 'A' && c != 'C' && c != 'G' && c != 'T') ++gap_bp[h];
            gap_bp[h] += clipped_gap[h];
        }
        // Two absolute bounds beside the relative ones, since identity is a fraction of the whole haplotype
        // and 2% of a long locus is an SV: an end overrun of --min-var or more, or a length difference over
        // the bound, is a structural error whatever the identity says.
        const long long overrun_bound = std::max<long long>(1, static_cast<long long>(options.min_var));
        auto overruns = [&](std::size_t h) {
            return lead_over[h] >= overrun_bound || trail_over[h] >= overrun_bound;
        };
        // The opposite: haplotype bases past the alignment chain that the walk does not match, an indel or a
        // divergent end the repair could not write in. --min-var or more of them is variation the graph lost, not
        // a difference it collapses by design. Only ACGT bases count: an assembly gap at an end holds no
        // variation to lose.
        std::vector<long long> lead_short(g.paths.size(), 0), trail_short(g.paths.size(), 0);
        for (std::size_t h = 0; h < g.paths.size(); ++h) {
            if (steps[h].empty()) continue;
            const auto [ps, pe, qs, qe] = chain_span[h];
            const std::string& q = hseq[h];
            const std::string w = spell_steps(steps[h]);
            const auto sz = [](long long x) { return static_cast<std::size_t>(x); };
            lead_short[h] = unmatched_end_bases(std::string(w.rbegin() + static_cast<long>(w.size() - sz(ps)), w.rend()),
                                                std::string(q.rbegin() + static_cast<long>(q.size() - sz(qs)), q.rend()));
            trail_short[h] = unmatched_end_bases(w.substr(std::min(w.size(), sz(pe))), q.substr(std::min(q.size(), sz(qe))));
        }
        auto falls_short = [&](std::size_t h) {
            return lead_short[h] >= overrun_bound || trail_short[h] >= overrun_bound;
        };
        auto length_off = [&](std::size_t h) {
            const long long bound = length_bound_for(hseq[h].size());
            return len_diff[h] > bound || -len_diff[h] > bound;
        };
        auto fails = [&](std::size_t h) {
            const double id = h < walk_id.size() ? walk_id[h] : -1.0;
            return matched[h] < options.min_matched_cover || id < 0.0 || id < options.min_recovered_identity ||
                   overruns(h) || falls_short(h) || length_off(h);
        };
        // What the N bases can explain. Cover: every N of the sequence compared (none aligns). Identity: only
        // the N the walk does not spell; those kept as gap or tip pieces cost no edits, and clipped ones are
        // not compared.
        auto gap_explains = [&](std::size_t h) {
            std::size_t gap = 0;
            for (const char c : hseq[h])
                if (c != 'A' && c != 'C' && c != 'G' && c != 'T') ++gap;
            const std::size_t unspelled = gap > walk_ambiguous[h] ? gap - walk_ambiguous[h] : 0;
            if (gap == 0 || walk_aln[h] <= unspelled || hseq[h].size() <= gap) return false;
            const double id = 1.0 - static_cast<double>(walk_edits[h] > unspelled ? walk_edits[h] - unspelled : 0) /
                                         static_cast<double>(walk_aln[h] - unspelled);
            const double cov = matched[h] * static_cast<double>(hseq[h].size()) /
                               static_cast<double>(hseq[h].size() - gap);
            return id >= options.min_recovered_identity && cov >= options.min_matched_cover &&
                   !overruns(h) && !falls_short(h) && !length_off(h);
        };
        // Why a haplotype falls below the thresholds, for the verdict's first example.
        auto why = [&](std::size_t h) -> std::string {
            const double id = h < walk_id.size() ? walk_id[h] : -1.0;
            return (overruns(h) ? " (walk overruns the haplotype's " +
                                      std::string(lead_over[h] >= overrun_bound ? "start" : "end") + " by " +
                                      std::to_string(std::max(lead_over[h], trail_over[h])) + " bp)"
                    : length_off(h) ? " (walk is " + std::to_string(len_diff[h] < 0 ? -len_diff[h] : len_diff[h]) +
                                          " bp " + (len_diff[h] < 0 ? "shorter" : "longer") + " than the haplotype)"
                    : falls_short(h) && id >= options.min_recovered_identity && matched[h] >= options.min_matched_cover
                        ? " (walk stops " + std::to_string(std::max(lead_short[h], trail_short[h])) +
                              " bp short of the haplotype's " +
                              std::string(lead_short[h] >= overrun_bound ? "start" : "end") + ")"
                    : id < 0.0 ? " (walk identity could not be computed)"
                               : " (walk identity " + fmt2(id) + ", matched cover " + fmt2(matched[h]) +
                                     (gap_explains(h) ? "; " + std::to_string(gap_bp[h]) + " bp of N, an assembly gap"
                                                      : std::string()) +
                                     ")") +
                   (inversion[h] == 1 ? ", a possible inversion its sequence does not confirm" : "");
        };
        std::size_t missing_paths = 0, dangling = 0, failing = 0;
        std::vector<char> dangling_walk(g.paths.size(), 0);
        std::string first_bad;
        for (std::size_t h = 0; h < g.paths.size(); ++h) {
            if (walks[h].empty()) {
                ++missing_paths;
                if (first_bad.empty()) first_bad = g.paths[h].name + " (no walk recovered)";
                continue;
            }
            for (std::size_t j = 1; j < walks[h].size(); ++j) {
                const std::uint64_t key = (static_cast<std::uint64_t>(walks[h][j - 1]) << 32) | walks[h][j];
                if (!arcs.count(key)) {
                    ++dangling;
                    dangling_walk[h] = 1;
                    if (first_bad.empty()) first_bad = g.paths[h].name + " (a step pair has no edge)";
                    break;
                }
            }
            // An identity that could not be computed is no evidence of a good recovery: it fails.
            if (fails(h)) {
                ++failing;
                if (first_bad.empty()) first_bad = g.paths[h].name + why(h);
            }
        }
        sum.paths_failing = failing;
        sum.dangling_steps = dangling;

        // The reference, if named, must come back exactly: every position and REF base downstream is read off
        // it. It seeds the graph, which makes that possible. --allow-loss does not override this: it accepts
        // lost haplotype variation, not wrong coordinates.
        bool ref_ok = true;
        if (!ref_name.empty()) {
            // Resolved at the start, so it is present.
            std::size_t ri = 0;
            while (g.paths[ri].name != ref_name) ++ri;
            if (walks[ri].empty()) {
                ref_ok = false;
                sum.reject_reason = "reference path was not recovered: " + g.paths[ri].name;
            } else if (!walk_exact[ri]) {
                const double rid = ri < walk_id.size() ? walk_id[ri] : -1.0;
                ref_ok = false;
                sum.reject_reason = "reference path is not reproduced exactly: " + g.paths[ri].name +
                                    " (identity " + fmt_exact(rid < 0 ? 0.0 : rid) + ", " +
                                    std::to_string(walk_edits[ri]) + " edits)";
            }
        }
        sum.reference_exact = ref_name.empty() || ref_ok;
        // Narrowing a clipped name can give it another path's name; two P lines of one name make a graph no
        // reader takes, so such a rebuild is refused whatever --allow-loss says.
        bool names_unique = true;
        {
            std::unordered_set<std::string> names;
            for (std::size_t h = 0; h < g.paths.size() && names_unique; ++h) {
                const std::string n = emitted_name(h);
                if (!names.insert(n).second) {
                    names_unique = false;
                    if (sum.reject_reason.empty())
                        sum.reject_reason = "two haplotypes would be emitted as " + n +
                                            " (a clipped name narrowed onto another path's)";
                }
            }
        }

        if (sum.reject_reason.empty()) {
            if (missing_paths > 0)
                sum.reject_reason = "haplotypes not recovered: " + std::to_string(missing_paths) + "; first: " + first_bad;
            else if (dangling > 0)
                sum.reject_reason = "walks with a step pair that has no edge: " + std::to_string(dangling) + "; first: " + first_bad;
            else if (failing > 0)
                sum.reject_reason = "haplotypes below the thresholds: " + std::to_string(failing) + "; first: " + first_bad;
        }
        sum.accepted = ref_ok && names_unique && sum.reject_reason.empty();
        const bool may_override = options.allow_loss && sum.reference_exact && names_unique;
        sum.overridden = !sum.accepted && may_override;

        // ---- copy-number modules (--mode folded) ----
        // Folded whenever the rest of the contract holds. A haplotype failing the per-path bounds whose folded
        // path spells it exactly then passes: the folded path, not the rebuilt walk, is what is written.
        std::vector<FoldedModule> folded_modules;
        // 0 not folded; 1 spells its haplotype exactly; 2 does not; 3 does not, but within the bounds walks are held to
        std::vector<char> after_fold(g.paths.size(), 0);
        auto rescued = [&](std::size_t h) { return after_fold[h] == 1 || after_fold[h] == 3; };
        // Whether what is written for a haplotype fails: its walk, unless a folded path replaces it.
        auto fails_written = [&](std::size_t h) { return fails(h) ? !rescued(h) : after_fold[h] == 2; };
        bool fold_logged = false;
        if (options.fold_modules && !staged_out.empty() && !emitted_gfa.empty() && ref_ok && names_unique &&
            missing_paths == 0 && dangling == 0) {
            std::vector<std::pair<std::string, std::string>> own;   // as emitted: clipped where clipped
            for (std::size_t h = 0; h < g.paths.size(); ++h)
                if (!walks[h].empty()) own.emplace_back(emitted_name(h), hseq[h]);
            ModuleFoldOptions fo;
            fo.threads = nthreads;
            fo.min_identity = options.fold_identity;
            fo.min_copy_identity = std::min(0.90, options.fold_identity);
            fo.min_copy_bp = options.fold_min_copy_bp;
            fo.gate_hub_degree = options.hub_degree;
            fo.hub_degree = options.hub_degree;   // one hub definition, the gate's
            const ModuleFoldResult fr = fold_modules(emitted_gfa, own, emitted_name(idx.front()), fo,
                                                     known_modules ? &*known_modules : nullptr);
            folded_modules = fr.modules;
            sum.modules_found = fr.modules.size();
            for (const FoldedModule& m : fr.modules) {
                sum.modules_folded += m.folded ? 1 : 0;
                if (!m.folded)
                    step("copy-number module at " + std::to_string(m.ref_start) + "-" + std::to_string(m.ref_end) +
                         (ref_name.empty() ? " of the seed not folded: " : " of the reference not folded: ") + m.reason);
            }
            if (sum.modules_folded > 0) {
                // Only a path the fold rewrote has an after-fold verdict; the others are written as walked.
                std::unordered_map<std::string, bool> spells(fr.spells_own.begin(), fr.spells_own.end());
                std::unordered_set<std::string> folded_names;
                for (const FoldedModule& m : fr.modules)
                    if (m.folded) folded_names.insert(m.folded_paths.begin(), m.folded_paths.end());
                for (std::size_t h = 0; h < g.paths.size(); ++h) {
                    if (!folded_names.count(emitted_name(h))) continue;
                    const auto it = spells.find(emitted_name(h));
                    if (it != spells.end()) after_fold[h] = it->second ? 1 : 2;
                }
                // Every folded path that does not spell its haplotype exactly is judged on what is written, by the
                // bounds walks are held to: outside the module it still spells differences under --min-var as the
                // backbone's, so exactness would rescue almost none. A haplotype whose walk fails passes on it; one
                // whose walk passes fails if it does not, since the walk is not what is written.
                {
                    std::vector<std::size_t> recheck;
                    // The folded path keeps the walk's ends, so an overrun or shortfall there is not rescued.
                    for (std::size_t h = 0; h < g.paths.size(); ++h)
                        if (after_fold[h] == 2 && !(fails(h) && (overruns(h) || falls_short(h)))) recheck.push_back(h);
                    if (!recheck.empty()) {
                        std::unordered_map<std::string, std::string> seg;
                        std::unordered_map<std::string, std::string> path_of;
                        std::unordered_set<std::string> wanted;
                        for (const std::size_t h : recheck) wanted.insert(emitted_name(h));
                        std::istringstream in(fr.gfa);
                        for (std::string line; std::getline(in, line);) {
                            if (line.rfind("S\t", 0) == 0) {
                                const std::size_t t1 = line.find('\t', 2), t2 = line.find('\t', t1 + 1);
                                seg[line.substr(2, t1 - 2)] = line.substr(t1 + 1, t2 == std::string::npos ? std::string::npos : t2 - t1 - 1);
                            } else if (line.rfind("P\t", 0) == 0) {
                                const std::size_t t1 = line.find('\t', 2), t2 = line.find('\t', t1 + 1);
                                const std::string nm = line.substr(2, t1 - 2);
                                if (wanted.count(nm)) path_of[nm] = line.substr(t1 + 1, t2 == std::string::npos ? std::string::npos : t2 - t1 - 1);
                            }
                        }
                        std::vector<char> within(recheck.size(), 0);
                        run_parallel(recheck.size(), nthreads, [&](std::size_t k) {
                            const std::size_t h = recheck[k];
                            const auto pit = path_of.find(emitted_name(h));
                            if (pit == path_of.end()) return;
                            std::string spelled;
                            std::stringstream steps_in(pit->second);
                            for (std::string st; std::getline(steps_in, st, ',');) {
                                if (st.size() < 2) return;
                                const auto sit = seg.find(st.substr(0, st.size() - 1));
                                if (sit == seg.end()) return;
                                spelled += st.back() == '-' ? reverse_complement(sit->second) : sit->second;
                            }
                            const long long diff = static_cast<long long>(spelled.size()) - static_cast<long long>(hseq[h].size());
                            within[k] = std::llabs(diff) <= length_bound_for(hseq[h].size()) &&
                                        identity_of(spelled, hseq[h]) >= options.min_recovered_identity;
                        });
                        for (std::size_t k = 0; k < recheck.size(); ++k) if (within[k]) after_fold[recheck[k]] = 3;
                    }
                }
                // The verdict on what is written. The reference must come back exactly after folding too, which
                // --allow-loss does not override.
                std::size_t ri = g.paths.size();
                for (std::size_t h = 0; h < g.paths.size() && !ref_name.empty(); ++h)
                    if (g.paths[h].name == ref_name) { ri = h; break; }
                if (ri < g.paths.size() && after_fold[ri] >= 2) {
                    sum.accepted = false;
                    sum.overridden = false;
                    sum.reference_exact = false;
                    sum.reject_reason = "reference path is not reproduced exactly after folding: " + g.paths[ri].name;
                } else {
                    std::size_t still = 0;
                    std::string first;
                    for (std::size_t h = 0; h < g.paths.size(); ++h)
                        if (fails_written(h)) {
                            ++still;
                            if (first.empty())
                                first = g.paths[h].name + (fails(h) ? why(h) : std::string(" (its folded path is outside the bounds)"));
                        }
                    sum.paths_failing = still;
                    if (still == 0) {
                        if (!sum.accepted) step("the folded paths of every haplotype below the thresholds pass them");
                        sum.accepted = true;
                        sum.overridden = false;
                        sum.reject_reason.clear();
                    } else {
                        sum.accepted = false;
                        sum.overridden = may_override;
                        sum.reject_reason = "haplotypes below the thresholds after folding: " + std::to_string(still) + "; first: " + first;
                    }
                }
                write_text(staged_out, fr.gfa, wants_gzip(options.out_path));
                fold_summary = {fr.nodes, fr.links, fr.maxdeg, fr.hubs, fr.selfloops};
                {
                    std::string copies;
                    for (const FoldedModule& m : fr.modules)
                        if (m.folded)
                            copies += std::string(copies.empty() ? "" : "; ") + "copies: " + std::to_string(m.min_copies) +
                                      "-" + std::to_string(m.max_copies) + (ref_name.empty() ? ", seed: " : ", reference: ") +
                                      std::to_string(m.reference_copies);
                    step("folded copy-number modules: " + std::to_string(sum.modules_folded) + " (" + copies + ")");
                    fold_logged = true;
                }
                for (std::size_t i = 0; i < fr.gfa.size();) {   // the folded graph's sequence, for its shape
                    const std::size_t e = std::min(fr.gfa.find('\n', i), fr.gfa.size());
                    if (fr.gfa.compare(i, 2, "S\t") == 0) {
                        const std::size_t t1 = fr.gfa.find('\t', i + 2);
                        const std::size_t t2 = std::min(fr.gfa.find('\t', t1 + 1), e);
                        if (t1 < e) folded_bp += t2 - t1 - 1;
                    }
                    i = e + 1;
                }
            }
        }

        if (!fold_logged)
            step(std::string("folded copy-number modules: 0") +
                 (!options.fold_modules ? " (--mode minigraph)"
                  : sum.modules_found > 0 ? ""
                  : missing_paths > 0 || dangling > 0 || !ref_ok || !names_unique ? " (not tried: the walks fail first)"
                  : ""));

        // A haplotype still failing whose copy gain or loss was left for the fold, which did not take it: a second
        // pass, on the graph already generated, writes those copy changes into the walks as indels.
        if (!sum.accepted && options.generated_gfa.empty() && options.fold_modules) {
            std::vector<std::string> retry;
            for (std::size_t h = 0; h < g.paths.size(); ++h)
                if (held_back[h] && fails_written(h)) retry.push_back(g.paths[h].name);
            // Only when it can lead to acceptance: every haplotype still failing is one it rewrites.
            std::size_t failing_now = 0;
            for (std::size_t h = 0; h < g.paths.size(); ++h) failing_now += fails_written(h);
            const bool fixable = ref_ok && names_unique && missing_paths == 0 && dangling == 0;
            if (!retry.empty() && ((fixable && retry.size() == failing_now) || options.allow_loss)) {
                step("copy changes the fold did not take: haplotypes: " + std::to_string(retry.size()) +
                     "; written in on a second pass");
                std::error_code ec;
                std::filesystem::remove(staged_out, ec);
                RebuildOptions second = options;
                second.generated_gfa = generated_copy;
                second.copies_for = std::move(retry);
                return rebuild_graph(second);
            }
        }

        // ---- report sidecar: one row per path, so the verdict can be read rather than trusted ----
        // Always <prefix>.rebuild.report.tsv, so it can never name the input or output graph.
        const std::string& report = report_path;
        StagedFile report_guard(report_staged);   // installed after the graph, below
        if (!report.empty()) {
            std::ofstream a(report_staged);
            if (!a) throw std::runtime_error("rebuild: cannot write report " + report_staged.string());
            {
                const char* const columns[] = {
                    "path", "original_bp", "recovered_steps", "envelope_cover", "matched_cover", "chain_identity",
                    "walk_identity", "gap_bp", "start_overrun_bp", "end_overrun_bp", "start_short_bp", "end_short_bp",
                    "length_diff_bp", "clipped_bp", "inversion_joined", "after_fold", "status"};
                for (const char* c : columns) a << (c == columns[0] ? "" : "\t") << c;
                a << '\n';
                for (std::size_t h = 0; h < g.paths.size(); ++h) {
                    const double id = h < walk_id.size() ? walk_id[h] : -1.0;
                    std::string status = "ok";
                    if (walks[h].empty()) status = "not_recovered";
                    else if (dangling_walk[h]) status = "dangling_walk";
                    else if (overruns(h)) status = "end_overrun";
                    else if (length_off(h)) status = "length_mismatch";
                    else if (matched[h] < options.min_matched_cover) status = "low_cover";
                    else if (id < 0.0) status = "identity_unavailable";
                    else if (id < options.min_recovered_identity) status = "low_identity";
                    else if (falls_short(h)) status = "end_short";   // the only reason: the other measures pass
                    if ((status == "low_cover" || status == "low_identity") && gap_explains(h))
                        status = "assembly_gap";
                    // Only where the inversion would explain the failure (the walk falls short of or differs
                    // from the haplotype): a dangling walk, an overrun or an assembly gap is its own reason.
                    if ((status == "low_cover" || status == "low_identity" || status == "length_mismatch") &&
                        inversion[h] == 1)
                        status = "possible_inversion";
                    if (h < g.paths.size() && g.paths[h].name == ref_name && walks[h].size() > 0 &&
                        (!walk_exact[h] || after_fold[h] >= 2))
                        status = "reference_inexact";
                    else if (status != "ok" && rescued(h))
                        status = "ok_after_fold";   // the walk fails, but its folded path, which is written, passes
                    else if (status == "ok" && after_fold[h] == 2)
                        status = "fails_after_fold";   // the walk passes, but its folded path, which is written, does not
                    // fmt_exact, as for the thresholds: two values both printed as 0.9900 can fall on
                    // opposite sides of a threshold, so the report could not reproduce its verdict.
                    a << g.paths[h].name << '\t' << input_len[h] << '\t' << walks[h].size() << '\t'
                      << fmt_exact(cover[h]) << '\t' << fmt_exact(matched[h]) << '\t'
                      << fmt_exact(chain_id[h]) << '\t'
                      << (id < 0.0 ? std::string("NA") : fmt_exact(id)) << '\t' << gap_bp[h] << '\t'
                      << lead_over[h] << '\t' << trail_over[h] << '\t' << lead_short[h] << '\t' << trail_short[h]
                      << '\t' << len_diff[h] << '\t'
                      << clipped[h] << '\t' << (inversion[h] == 2 ? 1 : 0) << '\t'
                      << (after_fold[h] == 0 ? "." : after_fold[h] == 1 ? "exact" : after_fold[h] == 3 ? "within" : "differs")
                      << '\t' << status << '\n';
                }
                // The verdict beside the per-path rows, recording what was done with the output: with
                // --allow-loss the rebuilt graph is written although the contract failed.
                a << "#verdict\t"
                  << (sum.accepted ? "accepted"
                                   : (may_override ? "accepted_with_override"
                                                                                : "rejected"))
                  << '\n';
                a << "#reason\t" << (sum.reject_reason.empty() ? std::string("-") : sum.reject_reason)
                  << '\n';
                a << "#min_recovered_identity\t" << fmt_exact(options.min_recovered_identity) << '\n';
                a << "#min_matched_cover\t" << fmt_exact(options.min_matched_cover) << '\n';
                a << "#max_end_overrun_bp\t" << (overrun_bound - 1) << '\n';
                a << "#max_end_short_bp\t" << (overrun_bound - 1) << '\n';
                a << "#max_length_diff_bp\t" << length_bound_for(0) << " (or 0.5% of a haplotype's length, if larger)\n";
                a << "#dangling_walks\t" << sum.dangling_steps << '\n';
                a << "#end_cuts\t" << sum.end_cuts << '\n';
                a << "#inversions_joined\t" << sum.inversions_joined << '\n';
                a << "#gaps_kept\t" << sum.gaps_kept << '\n';
                a << "#tips_added\t" << sum.tips_added << '\n';
                a << "#indels_added\t" << sum.indels_added << '\n';
                a << "#clipped_haplotypes\t" << sum.clipped_haplotypes << '\n';
                a << "#clipped_bp\t" << sum.clipped_bp << '\n';
                a << "#pruned_segments\t" << sum.pruned_segments << '\n';
                a << "#pruned_links\t" << sum.pruned_links << '\n';
                // One line per copy-number module found (--mode folded): where it lies on the reference,
                // its unit and period, and whether it was folded -- with the copies it holds, or why not.
                a << "#mode\t" << (options.fold_modules ? "folded" : "minigraph") << '\n';
                for (const FoldedModule& m : folded_modules) {
                    a << "#module\t" << m.ref_start << '-' << m.ref_end << "\tunit_bp=" << m.unit_bp
                      << "\tperiod_bp=" << m.period << "\tunit_from=" << m.unit_source;
                    if (m.folded)
                        a << "\tfolded\tcopies=" << m.min_copies << '-' << m.max_copies
                          << "\treference_copies=" << m.reference_copies << "\tcopy_identity=" << fmt_exact(m.copy_identity)
                          << "\tpaths_folded=" << m.paths_folded << '/' << m.paths
                          << "\texact_stretches=" << m.paths_exact << '/' << m.paths_folded
                          << "\tstretch_templates=" << m.stretch_templates
                          << "\tdifference_templates=" << m.difference_templates << "\thubs=" << m.hubs;
                    else
                        a << "\tnot_folded\t" << m.reason;
                    a << '\n';
                }
                a.flush();
                if (!a) throw std::runtime_error("rebuild: failed writing report " + report_staged.string());
                a.close();
                if (!a) throw std::runtime_error("rebuild: failed closing report " + report_staged.string());
                // Not installed here: the graph is committed first and the report after it, so a report never
                // describes an output that was not written.
            }
        }

        if (!staged_out.empty()) {
            // A report left by an earlier run goes before the graph is replaced, so an interruption cannot
            // leave the old report beside the new graph.
            drop_stale_report();
            if (sum.accepted || may_override) {
                if (sum.accepted) step("accepted: every haplotype within the thresholds");
                else warn("warning: accepted with --allow-loss: " + sum.reject_reason);
                commit_staged(staged_out, options.out_path);
                staged_guard.release();
            } else {
                std::error_code ec;
                std::filesystem::remove(staged_out, ec);
                warn("warning: rejected: " + sum.reject_reason + "; the input graph is written unchanged");
                // The same staged name, so the interrupt handler covers this copy too.
                copy_file(options.gfa_path, options.out_path, out_staged);
            }
        }

        // The graph is in place, so the report now describes something that exists. Failing to install it is
        // an error: it is the one record of what the rebuild accepted or lost.
        if (!report_staged.empty() && std::filesystem::exists(report_staged)) {
            try {
                commit_staged(report_staged, report);
                report_guard.release();
                sum.report_written = true;
            } catch (const std::exception& e) {
                throw std::runtime_error(std::string("rebuild: the graph was written to ") + options.out_path +
                                         " but its report could not be: " + e.what());
            }
        }
    }

    // ---- structural summary of what we emitted ----
    sum.out_nodes = n_kept;
    {
        // Per handle, as the input gate measures it, so the before and after numbers compare like with like.
        std::vector<std::array<std::unordered_set<std::uint32_t>, 2>> hs(pieces.size());
        std::vector<char> has_loop(pieces.size(), 0);   // NODES with a self-loop, as the input side counts
        for (const std::uint64_t k : emit_arcs) {
            const std::uint32_t v = static_cast<std::uint32_t>(k >> 32);   // oriented piece
            const std::uint32_t w = static_cast<std::uint32_t>(k);
            const std::uint32_t vs = v >> 1, ws = w >> 1;
            if (vs == ws) { has_loop[vs] = 1; continue; }
            // v is left by its end when forward, by its start when reverse; w is entered at its start
            // when forward, at its end when reverse.
            hs[vs][(v & 1) ? 0 : 1].insert(ws);
            hs[ws][(w & 1) ? 1 : 0].insert(vs);
        }
        sum.out_edges = emit_arcs.size();
        sum.out_selfloops = static_cast<std::size_t>(std::count(has_loop.begin(), has_loop.end(), 1));
        for (const auto& h : hs) {
            const std::size_t deg = std::max(h[0].size(), h[1].size());
            sum.out_maxdeg = std::max(sum.out_maxdeg, deg);
            if (deg >= options.hub_degree) ++sum.out_hubs;
        }
    }
    if (fold_summary) {   // what was written is the folded graph
        sum.out_nodes = (*fold_summary)[0];
        sum.out_edges = (*fold_summary)[1];
        sum.out_maxdeg = (*fold_summary)[2];
        sum.out_hubs = (*fold_summary)[3];
        sum.out_selfloops = (*fold_summary)[4];
    }
    if (sum.paths_recovered < sum.haplotypes) {
        std::cerr << "[rebuild " << hms() << "] warning: " << (sum.haplotypes - sum.paths_recovered)
                  << " of " << sum.haplotypes << " haplotypes could not be mapped back and have no P line\n";
    }
    {
        // The rebuilt graph's shape, as the input's was reported; when rejected, it is the graph not written.
        std::size_t out_bp = folded_bp;
        if (!fold_summary)
            for (std::size_t pid = 0; pid < pieces.size(); ++pid)
                if (piece_used[pid]) out_bp += pieces[pid].end - pieces[pid].beg;
        const double out_density = static_cast<double>(sum.out_nodes) / (static_cast<double>(max_len) / 1000.0);
        const double out_inflation = median_len == 0 ? 0.0 : static_cast<double>(out_bp) / static_cast<double>(median_len);
        step(std::string(sum.accepted || sum.overridden ? "rebuilt graph: " : "rebuilt graph (not written): ") +
             "haplotypes: " + std::to_string(sum.paths_recovered) + "; nodes: " + std::to_string(sum.out_nodes) +
             "; edges: " + std::to_string(sum.out_edges) + "; self-loops: " + std::to_string(sum.out_selfloops));
        step(shape_line(sum.accepted || sum.overridden ? "rebuilt graph" : "rebuilt graph (not written)",
                        sum.out_hubs >= options.min_hubs, sum.out_hubs, options.hub_degree, options.min_hubs,
                        sum.out_maxdeg, out_density, out_inflation, false));
    }
    return sum;
}

} // namespace panvar
