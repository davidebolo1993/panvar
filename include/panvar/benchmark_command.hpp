#pragma once

#include <string>
#include <vector>

namespace panvar {

// `panvar benchmark`: rebuilds each haplotype from the reference and call's events for it (then from
// call's VCF), aligns each rebuild to the haplotype, and reports the large differences left, each with
// its likely cause, plus a check of every carrier against its record.
int run_benchmark_command(const std::vector<std::string>& args);

} // namespace panvar
