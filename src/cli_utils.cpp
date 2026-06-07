#include "panvar/cli_utils.hpp"

#include <algorithm>
#include <cctype>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace panvar::cli {

void print_general_help() {
    std::cout
        << "panvar - modular pangenome graph toolkit\n\n"
        << "Usage:\n"
        << "  panvar <subcommand> [options]\n\n"
        << "Subcommands:\n"
        << "  bubble    Module 1: refine/import sites from 'vg snarls'\n"
        << "  inspect   Utility: inspect path walks through one called bubble\n"
        << "  allele    Module 2: allele extraction and clustering from module-1 sites\n"
        << "  call      Module 3: variant calling on module-2 clustered alleles\n"
        << "  describe  Module 4: per-bubble k-mer feature description\n\n"
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
