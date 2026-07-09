#pragma once

#include <string>
#include <vector>

namespace panvar {

// `panvar benchmark`: round-trip QV of the caller's own output. For each called bubble and
// each haplotype that traverses it, reconstruct the haplotype on the passed graph (its own
// bubble walk if we called a variant for it, else the reference walk), align to the true walk
// (edlib NW), and score QV = -10*log10(max(0.5, delta)/S). Reports the percentage of
// haplotypes in the cosigt QV bands (<17, 17-23, 23-33, >33).
int run_benchmark_command(const std::vector<std::string>& args);

} // namespace panvar
