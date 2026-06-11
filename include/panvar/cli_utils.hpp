#pragma once

#include <chrono>
#include <cstddef>
#include <string>
#include <vector>

namespace panvar::cli {

void print_general_help();

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

std::string trim_ascii(const std::string& text);
std::size_t parse_size_arg(const std::string& name, const std::string& value);
double parse_similarity_arg(const std::string& name, const std::string& value);
double parse_unit_fraction_arg(const std::string& name, const std::string& value);
std::string join_with_comma(const std::vector<std::string>& values);
void ensure_parent_dir_for_file(const std::string& path_text);

} // namespace panvar::cli
