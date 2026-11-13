#pragma once

#include <string>
#include <vector>

namespace panvar {

// `panvar background`: builds the background count store (every canonical 31-mer of a set of genomes with
// its count) that read recruitment classifies a locus's k-mers against.
int run_background_command(const std::vector<std::string>& args);

} // namespace panvar
