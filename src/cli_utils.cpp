#include "panvar/cli_utils.hpp"

#include <algorithm>
#include <atomic>
#include <random>
#include <cctype>
#include <csignal>
#include <cstring>
#include <functional>
#include <mutex>
#include <unordered_map>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <filesystem>
#include <fstream>
#include <map>
#include <regex>
#include <set>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <utility>

#if defined(_WIN32)
#include <io.h>
#else
#include <cerrno>
#include <fcntl.h>
#include <signal.h>
#include <sys/file.h>
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

std::string log_prefix(const std::string& module) { return "[" + module + " " + clock_hms() + "] "; }

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
    // Unwinding from an error: drawing the finished bar would show N/N for a run that failed.
    if (std::uncaught_exceptions() > 0) {
        if (active_ && !finished_) std::cerr << '\n';
        return;
    }
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
        << "panvar - variant sites, structural variants and copy number from a pangenome graph\n\n"
        << "Usage:\n"
        << "  panvar <command> [options]\n\n"
        << "Commands, in pipeline order:\n"
        << "  rebuild    Re-induce a fragmented locus graph before finding sites\n"
        << "  bubble     Find variant sites (bubbles), ordered along the reference\n"
        << "  panphorte  Fold tandem-repeat sites onto their repeat unit\n"
        << "  refine     Realign the haplotypes inside each site\n"
        << "  call       Call structural variants and copy number\n"
        << "  describe   Turn the sites into association features\n"
        << "  associate  Test the features for association with a phenotype\n\n"
        << "Read recruitment (in progress):\n"
        << "  background Count every k-mer of many genomes once, with anchors on a reference\n"
        << "  recruit    Build locus indexes against a background (--index-only)\n\n"
        << "Checks:\n"
        << "  inspect    Show every haplotype's allele and walk at each site\n"
        << "  benchmark  Rebuild each haplotype from call's output and list what is left\n\n"
        << "Run 'panvar <command> --help' for its options; 'panvar --version' prints the version.\n";
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
    // Digits only: stoull wraps a negative number ("-1" -> 18446744073709551615), stops at the first
    // non-digit ("1e6" -> 1, "10k" -> 10), and skips leading whitespace (newlines included) and a sign,
    // so each would be read as a different number than the one typed.
    if (value.empty() || !std::all_of(value.begin(), value.end(),
                                      [](char c) { return std::isdigit(static_cast<unsigned char>(c)); }))
        throw std::runtime_error("Invalid value for " + name + ": " + value + " (expected a non-negative integer)");
    try {
        std::size_t used = 0;
        const unsigned long long v = std::stoull(value, &used);
        if (used != value.size())
            throw std::runtime_error("trailing characters");
        return static_cast<std::size_t>(v);
    } catch (const std::exception&) {
        throw std::runtime_error("Invalid value for " + name + ": " + value + " (expected a non-negative integer)");
    }
}

namespace {
// The whole value must be the number: stod stops at the first character it cannot use, so "0.5%" and
// "0.8abc" read as 0.5 and 0.8, and it accepts hex floats ("0x0.8p0").
double parse_whole_double(const std::string& name, const std::string& value, const char* range) {
    const bool plain = !value.empty() &&
                       value.find_first_not_of("0123456789.eE+-") == std::string::npos;
    try {
        std::size_t used = 0;
        const double parsed = std::stod(value, &used);
        if (plain && used == value.size()) return parsed;
    } catch (const std::exception&) {
    }
    throw std::runtime_error("Invalid value for " + name + ": " + value + " (" + name + " " + range + ")");
}
}  // namespace

double parse_similarity_arg(const std::string& name, const std::string& value) {
    double parsed = parse_whole_double(name, value, "must be in (0,1] or (0,100]");
    if (parsed > 1.0 && parsed <= 100.0) {
        parsed /= 100.0;
    }
    if (!(parsed > 0.0 && parsed <= 1.0)) {
        throw std::runtime_error(name + " must be in (0,1] or (0,100]");
    }
    return parsed;
}

double parse_unit_fraction_arg(const std::string& name, const std::string& value) {
    const double parsed = parse_whole_double(name, value, "must be in [0,1]");
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

namespace {
// Staged path -> destination, so a write error names the file the user asked for.
std::mutex g_final_mutex;
std::unordered_map<std::string, std::string> g_final_of;
}  // namespace

void close_output(std::ofstream& out, const std::string& what) {
    out.flush();
    out.close();
    if (out) return;
    // Name the destination, not the staging sibling (see staging_path) the writer was handed.
    {
        std::lock_guard<std::mutex> lock(g_final_mutex);
        const auto it = g_final_of.find(what);
        if (it != g_final_of.end())
            throw std::runtime_error("Failed writing " + it->second + " (disk full or I/O error?)");
    }
    std::string name = what;
    const std::size_t slash = name.find_last_of('/');
    const std::size_t tmp = name.find("-tmp.", slash == std::string::npos ? 0 : slash);
    if (tmp != std::string::npos) {
        const std::size_t dot = name.rfind('.', tmp);
        if (dot != std::string::npos && (slash == std::string::npos || dot > slash)) name.erase(dot);
    }
    throw std::runtime_error("Failed writing " + name + " (disk full or I/O error?)");
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
    const std::string suffix = "." + module + "-tmp." + std::to_string(rd()) + "." + std::to_string(seq++);
    // A sibling of the destination, so installing it is a rename. A name too long for the suffix (most
    // file systems allow 255 bytes) is replaced by its hash.
    const std::filesystem::path p(final_path);
    const std::string name = p.filename().string();
    if (name.size() + suffix.size() <= 255) return final_path + suffix;
    const std::string short_name = "." + std::to_string(std::hash<std::string>{}(name)) + suffix;
    return p.has_parent_path() ? (p.parent_path() / short_name).string() : short_name;
}

namespace {
// Staged files to delete if the run is interrupted (Ctrl-C, kill), which skips the destructors that
// otherwise remove them. The handler may only call async-signal-safe functions, so the list is a
// push-only chain of C strings, never freed, and the handler only unlinks: a name already installed or
// removed no longer exists and is skipped by the system.
struct StagedNode {
    StagedNode* next;
    char path[1];
};
std::atomic<StagedNode*> g_staged{nullptr};
// While outputs are being installed, a signal is only recorded: the install finishes (or rolls back) and
// then re-raises it, so an interrupted run leaves one complete output family, never a mix.
volatile std::sig_atomic_t g_committing = 0, g_pending = 0;

#if !defined(_WIN32)
extern "C" void remove_staged_and_reraise(int sig) {
    if (g_committing) { g_pending = sig; return; }
    for (StagedNode* n = g_staged.load(); n != nullptr; n = n->next) ::unlink(n->path);
    std::signal(sig, SIG_DFL);
    std::raise(sig);
}

#endif

void register_staged(const std::string& path) {
#if !defined(_WIN32)
    static std::once_flag installed;
    std::call_once(installed, [] {
        // Only a signal left at its default is taken over: one the caller ignores (nohup) stays ignored.
        for (const int sig : {SIGINT, SIGTERM, SIGHUP}) {
            struct sigaction old {};
            if (::sigaction(sig, nullptr, &old) == 0 && old.sa_handler == SIG_DFL) {
                struct sigaction act {};
                act.sa_handler = remove_staged_and_reraise;
                sigemptyset(&act.sa_mask);
                ::sigaction(sig, &act, nullptr);
            }
        }
    });
#endif
    auto* n = static_cast<StagedNode*>(std::malloc(sizeof(StagedNode) + path.size()));
    if (n == nullptr) return;
    std::memcpy(n->path, path.c_str(), path.size() + 1);
    n->next = g_staged.load();
    while (!g_staged.compare_exchange_weak(n->next, n)) {}
}
}  // namespace

void commit_staged(const std::string& staged, const std::string& final_path) {
    std::error_code ec;
    std::filesystem::rename(staged, final_path, ec);
    if (!ec) return;

    // Never copy straight onto the destination: a copy interrupted half-way would leave a truncated
    // file where a complete one was, which is the one thing this function exists to prevent.
    //
    // A staged file is always a sibling of its destination, so a failed rename here is not the ordinary
    // cross-filesystem case -- it means something unexpected. Copy to a SECOND sibling and rename that
    // complete copy into place, so the destination is only ever replaced atomically. If any step fails
    // the destination is untouched and the staged file is left for the caller's cleanup.
    const std::string bridge = staging_path(final_path, "commit");
    register_staged(bridge);
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

namespace {
// The key two paths name the same file by. weakly_canonical resolves links, `.`, `..` and the case of
// components that exist, but not the case of the file names an output will create, so on macOS's
// case-insensitive volumes X.BUBBLES.CSV and x.bubbles.csv would pass as two files and one output would
// be lost. Case is folded there; a case-sensitive volume on macOS is rare enough that refusing two
// outputs differing only in case costs nothing.
std::string same_file_key(const std::string& p) {
    std::error_code ec;
    const auto c = std::filesystem::weakly_canonical(p, ec);
    std::string key = (ec || c.empty()) ? std::filesystem::path(p).string() : c.string();
#ifdef __APPLE__
    for (char& ch : key) ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
#endif
    return key;
}
}  // namespace

void require_regular_input(const std::string& module, const std::string& flag, const std::string& path) {
    std::error_code ec;
    if (!std::filesystem::exists(path, ec)) return;   // the reader reports a missing file in its own words
    // A path that still resolves to /dev/stdin or /dev/fd/N is an open descriptor, even when it is
    // redirected from a regular file: reopening it can share the descriptor's offset (macOS), so the
    // second read starts at the end. Where the system resolves it to the file itself, it reopens normally.
    const std::string resolved = std::filesystem::weakly_canonical(path, ec).string();
    const bool descriptor = !ec && (resolved == "/dev/stdin" || resolved.rfind("/dev/fd/", 0) == 0);
    if (descriptor || !std::filesystem::is_regular_file(path, ec))
        throw std::runtime_error(module + ": " + flag + " " + path + " is not a regular file that can be "
                                 "reopened; " + module + " reads it more than once, so it cannot come from "
                                 "a pipe, a process substitution or standard input -- write it to a file first");
}

std::string renamed_input_flag(const std::string& arg, const std::string& module) {
    static const std::map<std::string, std::string> renamed = {
        {"--bubble-prefix-in", "--bubble-prefix"}, {"--bubbles-csv-in", "--bubbles-csv"},
        {"--panphorte-prefix-in", "--panphorte-prefix"}, {"--snarls-in", "--snarls"}};
    static std::set<std::string> warned;
    const auto it = renamed.find(arg);
    if (it == renamed.end()) return arg;
    if (warned.insert(arg).second)
        std::cerr << log_prefix(module) << "warning: " << arg << " is now " << it->second << "; the old name will stop working\n";
    return it->second;
}

void reject_output_collisions(const std::string& module,
                              const std::vector<std::string>& outputs,
                              const std::vector<std::string>& inputs) {
    const auto canon = [](const std::string& p) { return same_file_key(p); };
    std::map<std::string, std::string> seen;   // canonical -> the spelling the user gave
    for (const std::string& out : outputs) {
        if (out.empty()) continue;
        std::error_code dir_ec;
        if (std::filesystem::is_directory(out, dir_ec)) {
            throw std::runtime_error(module + ": output '" + out + "' is an existing directory");
        }
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

namespace {
#if !defined(_WIN32)
// The process a hold file names (<output>.<module>-hold.<pid>-tmp.<rand>.<n>, or -holding.<pid>- while it is
// being made), or 0 for any other name. `making` tells which.
pid_t hold_owner(const std::string& name, bool& making) {
    static const std::regex hold_re(R"(-hold(ing)?\.([0-9]{1,9})-tmp\.[0-9]+\.[0-9]+$)");
    std::smatch m;
    if (!std::regex_search(name, m, hold_re)) return 0;
    making = m[1].matched;
    return static_cast<pid_t>(std::stol(m[2].str()));
}

// Remove the hold files left in `dir` by runs that ended without their cleanup (SIGKILL, abort). One
// goes only when its process no longer exists and its lock is free: a run on another host or in another
// pid namespace sharing the directory keeps its lock, though its pid means nothing here.
void remove_stale_holds(const std::filesystem::path& dir) {
    std::error_code ec;
    for (std::filesystem::directory_iterator it(dir, ec), end; !ec && it != end; it.increment(ec)) {
        bool making = false;
        const pid_t pid = hold_owner(it->path().filename().string(), making);
        if (pid <= 0 || ::kill(pid, 0) == 0 || errno != ESRCH) continue;
        // A hold still being made is unlocked for a moment; only one a minute old is a leftover.
        if (making) {
            const auto age = std::filesystem::file_time_type::clock::now() -
                             std::filesystem::last_write_time(it->path(), ec);
            if (ec || age < std::chrono::minutes(1)) { ec.clear(); continue; }
        }
        const int fd = ::open(it->path().c_str(), O_RDWR | O_CLOEXEC);
        if (fd < 0) continue;
        if (::flock(fd, LOCK_EX | LOCK_NB) == 0) ::unlink(it->path().c_str());
        ::close(fd);
    }
}
#endif

// Remove a hold file, then release its lock.
void drop_hold(const std::pair<std::string, int>& hold) {
    std::error_code ec;
    std::filesystem::remove(hold.first, ec);
#if !defined(_WIN32)
    if (hold.second >= 0) ::close(hold.second);
#endif
}
}  // namespace

StagedOutputs::StagedOutputs(std::string module) : module_(std::move(module)) {}

StagedOutputs::~StagedOutputs() {
    std::error_code ec;
    for (const auto& [staged, _final] : pending_) std::filesystem::remove(staged, ec);
    for (const auto& hold : holds_) drop_hold(hold);
    // A run that did not commit takes back the directories it made (commit() forgets them); only empty
    // ones go, deepest first (a child's spelling is longer than its parent's). Another run still using
    // one has its hold file there, so the directory is not empty and stays.
    std::sort(made_dirs_.begin(), made_dirs_.end(),
              [](const std::string& a, const std::string& b) { return a.size() > b.size(); });
    for (const std::string& dir : made_dirs_) {
        std::error_code is_empty_ec;
        if (std::filesystem::is_empty(dir, is_empty_ec) && !is_empty_ec) std::filesystem::remove(dir, ec);
    }
}

std::string StagedOutputs::stage(const std::string& final_path) {
    // A directory at the destination would otherwise be set aside at commit like a previous output and
    // replaced by the file, with no message.
    std::error_code dir_ec;
    if (std::filesystem::is_directory(final_path, dir_ec)) {
        throw std::runtime_error(module_ + ": output '" + final_path + "' is an existing directory");
    }
    const std::string name = std::filesystem::path(final_path).filename().string();
    if (name.size() > 255)
        throw std::runtime_error(module_ + ": output name too long (" + std::to_string(name.size()) +
                                 " bytes; file systems allow 255): " + final_path);
    // Note the directories this creates, so a run that fails before committing can remove them.
    std::error_code exists_ec;
    for (std::filesystem::path dir = std::filesystem::path(final_path).parent_path();
         !dir.empty() && !std::filesystem::exists(dir, exists_ec); dir = dir.parent_path()) {
        made_dirs_.push_back(dir.string());
        if (dir == dir.parent_path()) break;
    }
    ensure_parent_dir_for_file(final_path);
    // An empty file in each output directory until this run commits or ends. A concurrent run that made
    // the directory and then fails removes it only when it is empty, so the hold keeps it for this run
    // even before any output is written there. Its name carries this process's id and it stays locked
    // while the run lives, so a later run can tell and remove the hold of a run that was killed.
    const std::filesystem::path parent = std::filesystem::path(final_path).parent_path();
    if (!parent.empty() && held_dirs_.insert(same_file_key(parent.string())).second) {
#if !defined(_WIN32)
        remove_stale_holds(parent);
        // Made and locked under a name no cleaner reads, then renamed: a hold is never seen unlocked.
        std::string hold = staging_path(final_path, module_ + "-hold." + std::to_string(::getpid()));
        const std::string making = staging_path(final_path, module_ + "-holding." + std::to_string(::getpid()));
        register_staged(making);
        register_staged(hold);
        int fd = ::open(making.c_str(), O_RDWR | O_CREAT | O_CLOEXEC, 0644);
        if (fd < 0) {   // the directory went between the check above and now: make it again
            ensure_parent_dir_for_file(final_path);
            fd = ::open(making.c_str(), O_RDWR | O_CREAT | O_CLOEXEC, 0644);
        }
        if (fd >= 0) ::flock(fd, LOCK_EX | LOCK_NB);
        if (fd < 0 || ::rename(making.c_str(), hold.c_str()) != 0) hold = making;   // still holds the directory
        holds_.emplace_back(hold, fd);
#else
        const std::string hold = staging_path(final_path, module_ + "-hold");
        register_staged(hold);
        std::ofstream made(hold);
        if (!made) {
            ensure_parent_dir_for_file(final_path);
            made.open(hold);
        }
        holds_.emplace_back(hold, -1);
#endif
    }
    const std::string staged = staging_path(final_path, module_);
    register_staged(staged);
    {
        std::lock_guard<std::mutex> lock(g_final_mutex);
        g_final_of[staged] = final_path;
    }
    pending_.emplace_back(staged, final_path);
    return staged;
}

void StagedOutputs::commit() {
    // A module's whole output family installs together or not at all: every module page promises that a
    // refusal leaves the previous complete output intact, so the guarantee lives here once rather than
    // in a partial re-implementation per module.
    //
    // Any destination that already exists is moved aside first, so a failure can put it back. On
    // success the reserves are dropped; on failure this run's installs are removed and the reserves
    // restored, in reverse order.
    //
    // The entry is recorded as soon as the set-aside succeeds, not after the install: otherwise a
    // failure between the two would leave that destination unowned, rollback would skip it, and the old
    // output would stay parked under a name nobody would look for. `new_installed` distinguishes "this
    // reserve needs restoring" from "a file of ours is sitting there and must be removed first".
    struct Installed {
        std::string final_path;
        std::string reserve;            // empty when the destination did not exist before this run
        bool new_installed = false;     // our file actually reached final_path
    };
    std::vector<Installed> installed;

    // Fault injection, test-only, at two distinct points. The rollback branch is otherwise
    // unreachable from a test: every natural way to make a rename fail here either fails earlier, at
    // staging, or is absorbed by the set-aside. AT fails before the Nth set-aside; AFTER_SETASIDE fails
    // between the Nth set-aside and its install, the window described above.
    std::size_t fail_at = 0;
    std::size_t fail_after_setaside = 0;
    if (const char* env = std::getenv("PANVAR_TEST_FAIL_COMMIT_AT")) {
        fail_at = static_cast<std::size_t>(std::strtoull(env, nullptr, 10));
    }
    if (const char* env = std::getenv("PANVAR_TEST_FAIL_COMMIT_AFTER_SETASIDE")) {
        fail_after_setaside = static_cast<std::size_t>(std::strtoull(env, nullptr, 10));
    }
    std::size_t installs = 0;
    // PANVAR_TEST_RAISE_AT_COMMIT: the run is killed with every staged file written.
    // PANVAR_TEST_RAISE_IN_COMMIT=N: killed after the Nth install, inside the window a signal must wait.
    if (std::getenv("PANVAR_TEST_RAISE_AT_COMMIT") != nullptr) std::raise(SIGTERM);
    std::size_t raise_in = 0;
    if (const char* env = std::getenv("PANVAR_TEST_RAISE_IN_COMMIT"))
        raise_in = static_cast<std::size_t>(std::strtoull(env, nullptr, 10));
    // A signal arriving from here on waits until the install or its rollback is complete.
    struct Deferred {
        Deferred() { g_committing = 1; }
        ~Deferred() {
            g_committing = 0;
            if (const int sig = g_pending) { g_pending = 0; std::raise(sig); }
        }
    } deferred;

    std::set<std::string> installed_keys;
    try {
        for (const auto& [staged, final_path] : pending_) {
            if (!std::filesystem::exists(staged)) continue;   // an optional output nobody wrote
            // Two outputs of this run naming one file (on macOS, possibly differing only in case):
            // the second would set the first aside as the "previous output" and then delete it.
            if (!installed_keys.insert(same_file_key(final_path)).second) {
                throw std::runtime_error(module_ + ": two outputs would be written to the same file: " +
                                         final_path);
            }
            ++installs;
            if (fail_at != 0 && installs == fail_at) {
                throw std::runtime_error("PANVAR_TEST_FAIL_COMMIT_AT: injected failure installing "
                                         + final_path);
            }
            std::string reserve;
            // symlink_status, not exists(): a dangling symlink at the destination does not "exist",
            // so it would not be set aside, and a rollback could not put it back.
            // A directory made at the destination after staging is not a previous output to replace.
            if (std::filesystem::is_directory(std::filesystem::symlink_status(final_path)))
                throw std::runtime_error(module_ + ": output '" + final_path + "' is an existing directory");
            if (std::filesystem::exists(std::filesystem::symlink_status(final_path))) {
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
            if (raise_in != 0 && installs == raise_in) std::raise(SIGTERM);
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
    made_dirs_.clear();   // they hold this run's outputs now
    for (const auto& hold : holds_) drop_hold(hold);
    holds_.clear();
}

} // namespace cli
} // namespace panvar
