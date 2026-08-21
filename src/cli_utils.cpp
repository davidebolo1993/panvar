#include "panvar/cli_utils.hpp"

#include <algorithm>
#include <atomic>
#include <random>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <filesystem>
#include <map>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <utility>

#if defined(_WIN32)
#include <io.h>
#else
#include <unistd.h>
#endif

namespace panvar::cli {

double elapsed_seconds(const std::chrono::steady_clock::time_point& start) {
    const auto now = std::chrono::steady_clock::now();
    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
    return static_cast<double>(ms) / 1000.0;
}

namespace {
// Current wall-clock time as HH:MM:SS for log-line prefixes.
std::string clock_hms() {
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
} // namespace

RunLog::RunLog(std::string module, bool quiet)
    : module_(std::move(module)), quiet_(quiet), start_(std::chrono::steady_clock::now()) {}

void RunLog::line(const std::string& message) const {
    if (quiet_) {
        return;
    }
    std::cerr << '[' << module_ << ' ' << clock_hms() << "] " << message << '\n';
    std::cerr.flush();
}

void RunLog::info(const std::string& message) const { line(message); }

void RunLog::wrote(const std::vector<std::string>& paths) const {
    if (quiet_ || paths.empty()) {
        return;
    }
    // Compact form: if every path shares one parent directory, name it once and list
    // the basenames; otherwise fall back to the full paths.
    const std::filesystem::path dir0 = std::filesystem::path(paths.front()).parent_path();
    bool same_dir = true;
    for (const std::string& p : paths) {
        if (std::filesystem::path(p).parent_path() != dir0) {
            same_dir = false;
            break;
        }
    }
    std::string joined;
    for (std::size_t i = 0; i < paths.size(); ++i) {
        if (i != 0) {
            joined += ", ";
        }
        joined += same_dir ? std::filesystem::path(paths[i]).filename().string() : paths[i];
    }
    std::string head = "wrote " + std::to_string(paths.size()) + " file" +
                       (paths.size() == 1 ? "" : "s");
    if (same_dir && !dir0.empty()) {
        head += " in " + dir0.string() + "/";
    }
    line(head + ": " + joined);
}

void RunLog::done() const {
    char buf[32];
    std::snprintf(buf, sizeof(buf), "done in %.1fs", elapsed_seconds(start_));
    line(buf);
}

bool stderr_is_tty() {
#if defined(_WIN32)
    return _isatty(_fileno(stderr)) != 0;
#else
    return isatty(fileno(stderr)) != 0;
#endif
}

ProgressBar::ProgressBar(std::string label, std::size_t total)
    : label_(std::move(label)), total_(total) {
    // Only animate on an interactive stderr; skip entirely for empty work or when
    // suppressed by an empty label (the mechanism used to honor --quiet).
    active_ = total_ > 0 && !label_.empty() && stderr_is_tty();
    if (active_) {
        render();
    }
}

ProgressBar::~ProgressBar() {
    done();
}

void ProgressBar::tick() {
    if (current_ < total_) {
        ++current_;
    }
    if (active_ && !finished_) {
        render();
    }
}

void ProgressBar::done() {
    if (finished_) {
        return;
    }
    finished_ = true;
    if (active_) {
        current_ = total_;
        render();
        std::cerr << '\n';
        std::cerr.flush();
    }
}

void ProgressBar::render() {
    constexpr std::size_t kWidth = 30;
    const std::size_t filled =
        total_ == 0 ? kWidth : (current_ * kWidth) / total_;
    std::cerr << '\r' << label_ << " [";
    for (std::size_t i = 0; i < kWidth; ++i) {
        std::cerr << (i < filled ? '#' : '-');
    }
    std::cerr << "] " << current_ << '/' << total_;
    std::cerr.flush();
}

void print_general_help() {
    std::cout
        << "panvar - modular pangenome graph toolkit\n\n"
        << "Usage:\n"
        << "  panvar <subcommand> [options]\n\n"
        << "Subcommands:\n"
        << "  rebuild    Module 0: re-induce a pathological locus graph before decomposition\n"
        << "  bubble     Module 1: extract/refine bubble sites from a GFA (no external tools)\n"
        << "  panphorte  Module 2: normalize tandem-repeat bubbles into a compact GFA\n"
        << "  refine     Module 2b: POA-realign bubble interiors to remove pggb alignment artifacts\n"
        << "  call       Module 3: graph-native structural variant calling\n"
#ifdef PANVAR_ENABLE_EXPERIMENTAL_GENOTYPE
        << "  genotype   Module 3b: per-bubble genotyping of a short-read sample"
           " (EXPERIMENTAL, not release-reviewed)\n"
#endif
        << "  describe   Module 4: per-bubble k-mer feature description\n"
        << "  associate  Module 5: GWAS on describe genotypes vs a phenotype/covariate table\n"
        << "  inspect    Utility: inspect path walks through one or all called bubbles\n"
        << "  benchmark  Utility: round-trip QV of the caller's own output (cosigt-style bands)\n\n"
        << "Run 'panvar <subcommand> --help' for options.\n";
}

std::string trim_ascii(const std::string& text) {
    std::size_t lo = 0;
    while (lo < text.size() && std::isspace(static_cast<unsigned char>(text[lo]))) {
        ++lo;
    }
    std::size_t hi = text.size();
    while (hi > lo && std::isspace(static_cast<unsigned char>(text[hi - 1]))) {
        --hi;
    }
    return text.substr(lo, hi - lo);
}

std::size_t parse_size_arg(const std::string& name, const std::string& value) {
    try {
        return static_cast<std::size_t>(std::stoull(value));
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid value for " + name + ": " + value);
    }
}

double parse_similarity_arg(const std::string& name, const std::string& value) {
    double parsed = 0.0;
    try {
        parsed = std::stod(value);
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid value for " + name + ": " + value);
    }
    if (parsed > 1.0 && parsed <= 100.0) {
        parsed /= 100.0;
    }
    if (!(parsed > 0.0 && parsed <= 1.0)) {
        throw std::runtime_error(name + " must be in (0,1] or (0,100]");
    }
    return parsed;
}

double parse_unit_fraction_arg(const std::string& name, const std::string& value) {
    double parsed = 0.0;
    try {
        parsed = std::stod(value);
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid value for " + name + ": " + value);
    }
    if (!(parsed >= 0.0 && parsed <= 1.0)) {
        throw std::runtime_error(name + " must be in [0,1]");
    }
    return parsed;
}

std::string join_with_comma(const std::vector<std::string>& values) {
    if (values.empty()) {
        return {};
    }
    std::string out;
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i > 0) {
            out += ", ";
        }
        out += values[i];
    }
    return out;
}

void ensure_parent_dir_for_file(const std::string& path_text) {
    const std::filesystem::path p(path_text);
    const auto parent = p.parent_path();
    if (!parent.empty()) {
        std::filesystem::create_directories(parent);
    }
}

} // namespace panvar::cli

namespace panvar {
namespace cli {

std::string staging_path(const std::string& final_path, const std::string& module) {
    static std::atomic<unsigned> seq{0};
    std::random_device rd;
    const std::string tag = std::to_string(rd()) + "." + std::to_string(seq++);
    return final_path + "." + module + "-tmp." + tag;
}

void commit_staged(const std::string& staged, const std::string& final_path) {
    std::error_code ec;
    std::filesystem::rename(staged, final_path, ec);
    if (!ec) return;

    // The old fallback copied STRAIGHT ONTO the destination and then removed the staged file whether
    // or not the copy reported an error. That is the one thing this function exists to prevent: a
    // copy interrupted half-way leaves a truncated file where a complete one used to be, and the
    // caller's rollback has nothing to restore from because the reserve is a different file entirely.
    //
    // A staged file is always a sibling of its destination, so a failed rename here is not the
    // ordinary cross-filesystem case -- it means something unexpected. Copy to a SECOND sibling and
    // rename that complete copy into place, so the destination is only ever replaced atomically. If
    // any step fails the destination is untouched and the staged file is left for the caller's
    // cleanup.
    const std::string bridge = staging_path(final_path, "commit");
    std::error_code copy_ec;
    std::filesystem::copy_file(staged, bridge, std::filesystem::copy_options::overwrite_existing,
                               copy_ec);
    if (copy_ec) {
        std::error_code rm;
        std::filesystem::remove(bridge, rm);
        throw std::runtime_error("cannot move output into place: " + final_path + " (" + ec.message() +
                                 "; the copy fallback also failed: " + copy_ec.message() + ")");
    }
    std::error_code ren_ec;
    std::filesystem::rename(bridge, final_path, ren_ec);
    if (ren_ec) {
        std::error_code rm;
        std::filesystem::remove(bridge, rm);
        throw std::runtime_error("cannot move output into place: " + final_path + " (" +
                                 ren_ec.message() + ")");
    }
    std::error_code rm;
    std::filesystem::remove(staged, rm);
}

void reject_output_collisions(const std::string& module,
                              const std::vector<std::string>& outputs,
                              const std::vector<std::string>& inputs) {
    const auto canon = [](const std::string& p) {
        std::error_code ec;
        const auto c = std::filesystem::weakly_canonical(p, ec);
        return (ec || c.empty()) ? std::filesystem::path(p).string() : c.string();
    };
    std::map<std::string, std::string> seen;   // canonical -> the spelling the user gave
    for (const std::string& out : outputs) {
        if (out.empty()) continue;
        const auto [it, fresh] = seen.emplace(canon(out), out);
        if (!fresh) {
            throw std::runtime_error(module + ": two outputs would be written to the same file: " +
                                     (it->second == out ? out : it->second + " and " + out));
        }
    }
    for (const std::string& in : inputs) {
        if (in.empty()) continue;
        const auto it = seen.find(canon(in));
        if (it != seen.end()) {
            throw std::runtime_error(module + ": output '" + it->second +
                                     "' is the same file as input '" + in +
                                     "'; refusing to overwrite it");
        }
    }
}

StagedOutputs::StagedOutputs(std::string module) : module_(std::move(module)) {}

StagedOutputs::~StagedOutputs() {
    std::error_code ec;
    for (const auto& [staged, _final] : pending_) std::filesystem::remove(staged, ec);
}

std::string StagedOutputs::stage(const std::string& final_path) {
    ensure_parent_dir_for_file(final_path);
    const std::string staged = staging_path(final_path, module_);
    pending_.emplace_back(staged, final_path);
    return staged;
}

void StagedOutputs::commit() {
    // Installing sequentially made every module's output family non-atomic: a rename that failed on
    // the fourth of six files left three new outputs beside three stale ones, with a non-zero exit
    // and nothing saying which was which. Each module page promises the opposite -- a refusal leaves
    // the previous complete output intact -- so the guarantee belongs here, once, rather than in
    // eight partial re-implementations.
    //
    // Any destination that already exists is moved aside first, so a failure can put it back. On
    // success the reserves are dropped; on failure this run's installs are removed and the reserves
    // restored, in reverse order.
    // The entry is recorded as soon as the SET-ASIDE succeeds, not after the install. Recording it
    // afterwards left a window with no owner: if installation failed between moving the old file to
    // its reserve and putting the new one in place, that destination was absent from `installed`, so
    // rollback skipped it -- the old output stayed in a `*-prev` file under a name nobody would look
    // for, and the destination was simply gone. `new_installed` distinguishes "this reserve needs
    // restoring" from "a file of ours is sitting there and must be removed first".
    struct Installed {
        std::string final_path;
        std::string reserve;            // empty when the destination did not exist before this run
        bool new_installed = false;     // our file actually reached final_path
    };
    std::vector<Installed> installed;

    // Fault injection, test-only, at two distinct points. The rollback branch is otherwise
    // unreachable from a test -- every natural way to make a rename fail here either fails earlier,
    // at staging, or is absorbed by the set-aside -- and an untested rollback is worth about as much
    // as no rollback. AT fails BEFORE the Nth set-aside; AFTER_SETASIDE fails between the Nth
    // set-aside and its install, which is the window described above and the one the previous
    // fixture could not reach.
    std::size_t fail_at = 0;
    std::size_t fail_after_setaside = 0;
    if (const char* env = std::getenv("PANVAR_TEST_FAIL_COMMIT_AT")) {
        fail_at = static_cast<std::size_t>(std::strtoull(env, nullptr, 10));
    }
    if (const char* env = std::getenv("PANVAR_TEST_FAIL_COMMIT_AFTER_SETASIDE")) {
        fail_after_setaside = static_cast<std::size_t>(std::strtoull(env, nullptr, 10));
    }
    std::size_t installs = 0;

    try {
        for (const auto& [staged, final_path] : pending_) {
            if (!std::filesystem::exists(staged)) continue;   // an optional output nobody wrote
            ++installs;
            if (fail_at != 0 && installs == fail_at) {
                throw std::runtime_error("PANVAR_TEST_FAIL_COMMIT_AT: injected failure installing "
                                         + final_path);
            }
            std::string reserve;
            if (std::filesystem::exists(final_path)) {
                reserve = staging_path(final_path, module_ + "-prev");
                std::error_code ec;
                std::filesystem::rename(final_path, reserve, ec);
                if (ec) {
                    throw std::runtime_error("cannot set aside the previous output before replacing it: "
                                             + final_path);
                }
            }
            installed.push_back({final_path, reserve, false});
            if (fail_after_setaside != 0 && installs == fail_after_setaside) {
                throw std::runtime_error("PANVAR_TEST_FAIL_COMMIT_AFTER_SETASIDE: injected failure "
                                         "between set-aside and install for " + final_path);
            }
            commit_staged(staged, final_path);
            installed.back().new_installed = true;
        }
    } catch (...) {
        // A restore that itself fails is the one outcome the caller cannot recover from, so it is
        // named rather than swallowed: the previous output is gone and its reserve copy is the only
        // remaining version.
        std::vector<std::string> unrestored;
        for (auto it = installed.rbegin(); it != installed.rend(); ++it) {
            std::error_code ec;
            // Remove OUR file first, whether or not there was a predecessor: an output this run
            // created where none existed has to go too, or a failed run leaves half a new family
            // behind. Only then can a reserve be renamed back over the name.
            if (it->new_installed) std::filesystem::remove(it->final_path, ec);
            if (it->reserve.empty()) continue;
            std::error_code ren;
            std::filesystem::rename(it->reserve, it->final_path, ren);
            if (ren) unrestored.push_back(it->final_path + " (kept as " + it->reserve + ")");
        }
        if (!unrestored.empty()) {
            std::string msg = "output commit failed AND the previous output could not be restored for:";
            for (const auto& u : unrestored) msg += "\n  " + u;
            std::throw_with_nested(std::runtime_error(msg));
        }
        throw;   // pending_ is left intact so the destructor still clears the staged files
    }

    std::error_code ec;
    for (const auto& entry : installed) {
        if (!entry.reserve.empty()) std::filesystem::remove(entry.reserve, ec);
    }
    pending_.clear();
}

} // namespace cli
} // namespace panvar
