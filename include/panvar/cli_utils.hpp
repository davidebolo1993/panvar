#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace panvar::cli {

void print_general_help();

std::string trim_ascii(const std::string& text);
std::size_t parse_size_arg(const std::string& name, const std::string& value);
double parse_similarity_arg(const std::string& name, const std::string& value);
double parse_unit_fraction_arg(const std::string& name, const std::string& value);
std::string join_with_comma(const std::vector<std::string>& values);
void ensure_parent_dir_for_file(const std::string& path_text);

} // namespace panvar::cli
