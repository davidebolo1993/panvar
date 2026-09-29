#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace panvar {

// Parse a cosigt sample->haplotype table: "sample<TAB>hap1,hap2<TAB>..." (header row tolerated).
// Returns (haplotype-path -> samples carrying it, sample column order). A sample assigned the same
// haplotype twice appears twice under that haplotype. Shared by every module taking --samples so the
// diploid summation is defined in one place. `module` prefixes messages; when `known_haplotypes` is
// given, a haplotype outside it is an error naming `known_label`.
std::pair<std::unordered_map<std::string, std::vector<std::string>>, std::vector<std::string>>
read_cosigt_table(const std::string& path,
                  const std::string& module,
                  const std::unordered_set<std::string>* known_haplotypes = nullptr,
                  const char* known_label = "");

} // namespace panvar
