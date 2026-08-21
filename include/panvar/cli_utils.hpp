#pragma once

#include <chrono>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace panvar::cli {

void print_general_help();

// True when stderr is an interactive terminal (gate live progress so piped logs stay clean).
bool stderr_is_tty();

// Wall-clock seconds elapsed since `start`, used for progress logging.
double elapsed_seconds(const std::chrono::steady_clock::time_point& start);

// Minimal single-threaded progress bar printed to stderr (keeps stdout clean for
// the machine-readable summary). Renders as `\rlabel [####----] 42/120` and is
// rewritten in place on each tick; call done() (or let the destructor run) to emit
// the trailing newline. Rendering is skipped when stderr is not a TTY so piped logs
// stay clean.
class ProgressBar {
public:
    ProgressBar(std::string label, std::size_t total);
    ~ProgressBar();

    ProgressBar(const ProgressBar&) = delete;
    ProgressBar& operator=(const ProgressBar&) = delete;

    // Advance by one unit and redraw.
    void tick();
    // Finish: redraw at 100% and emit a newline (idempotent).
    void done();

private:
    void render();

    std::string label_;
    std::size_t total_ = 0;
    std::size_t current_ = 0;
    bool active_ = false;
    bool finished_ = false;
};

// Consistent, timestamped run logging for the command layer. Every line goes to
// stderr (so stdout stays free for any machine-readable output) as
// `[module HH:MM:SS] message`. Honors --quiet (emits nothing when quiet). Construct
// once per command; call info()/wrote() for events and done() for the closing line.
class RunLog {
public:
    RunLog(std::string module, bool quiet);
    // One timestamped event line.
    void info(const std::string& message) const;
    // "wrote: a, b, c" — the full set of output files, reported together.
    void wrote(const std::vector<std::string>& paths) const;
    // Closing line: "done in X.Ys" (wall time since construction).
    void done() const;

private:
    void line(const std::string& message) const;
    std::string module_;
    bool quiet_ = false;
    std::chrono::steady_clock::time_point start_;
};

std::string trim_ascii(const std::string& text);
std::size_t parse_size_arg(const std::string& name, const std::string& value);
double parse_similarity_arg(const std::string& name, const std::string& value);
double parse_unit_fraction_arg(const std::string& name, const std::string& value);
std::string join_with_comma(const std::vector<std::string>& values);
void ensure_parent_dir_for_file(const std::string& path_text);

// A unique sibling of `final_path`: outputs are written here and renamed into place, so an
// interrupted or failed run never leaves a half-written file where a complete one is expected.
std::string staging_path(const std::string& final_path, const std::string& module);

// Rename a staged file into place, falling back to copy-then-remove across filesystems.
void commit_staged(const std::string& staged, const std::string& final_path);

// Every output of one command, staged together and committed only once the command has succeeded.
// Anything still staged when the set is destroyed is removed, so a throw leaves no partial output.
// Refuse, before anything is opened, two outputs that name one file or any output that names an
// input. Empty entries are ignored, so a caller can pass optional destinations unconditionally.
// Comparison is by weakly_canonical path, so a symlink or a second mount point is caught too.
void reject_output_collisions(const std::string& module,
                              const std::vector<std::string>& outputs,
                              const std::vector<std::string>& inputs);

class StagedOutputs {
public:
    explicit StagedOutputs(std::string module);
    ~StagedOutputs();
    StagedOutputs(const StagedOutputs&) = delete;
    StagedOutputs& operator=(const StagedOutputs&) = delete;

    // Register a destination and return the path to write to instead.
    std::string stage(const std::string& final_path);
    // Move every staged file to its destination, in registration order, as ONE transaction: a
    // failure part-way through removes what this run installed and restores whatever it displaced,
    // so a caller never leaves a half-updated output family behind. Throws on failure.
    void commit();

private:
    std::string module_;
    std::vector<std::pair<std::string, std::string>> pending_;   // (staged, final)
};

} // namespace panvar::cli
