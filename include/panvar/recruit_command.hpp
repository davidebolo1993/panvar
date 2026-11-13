#pragma once

#include <string>
#include <vector>

namespace panvar {

// `panvar recruit`: builds (or reuses) each locus's index against a background, then recruits reads.
int run_recruit_command(const std::vector<std::string>& args);

} // namespace panvar
