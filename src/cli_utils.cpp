#include "panvar/cli_utils.hpp"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <exception>
#include <filesystem>
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
        << "  bubble     Module 1: extract/refine bubble sites from a GFA (no external tools)\n"
        << "  panphorte  Module 2: normalize tandem-repeat bubbles into a compact GFA\n"
        << "  call       Module 3: graph-native structural variant calling\n"
        << "  describe   Module 4: per-bubble k-mer feature description\n"
        << "  inspect    Utility: inspect path walks through one or all called bubbles\n\n"
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
