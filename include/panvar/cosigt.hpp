#pragma once

#include <functional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace panvar {

// A cosigt sample -> haplotype table: "sample<TAB>hap1<TAB>hap2..." (a header row is tolerated), a
// haplotype assembled in several contigs written "contig1,contig2" in its column. Shared by every module
// taking --samples, so a sample's haplotypes mean one thing everywhere.
struct CosigtTable {
    std::unordered_map<std::string, std::vector<std::string>> path_to_samples;   // contig -> sample, per slot naming it
    std::vector<std::string> sample_order;
    std::unordered_map<std::string, std::vector<std::vector<std::string>>> slots;   // sample -> its slots' contigs
};

// Read the table. `module` prefixes messages; when `known_haplotypes` is given, a contig outside it is an
// error naming `known_label`. `quiet` drops the mixed-ploidy note.
CosigtTable read_cosigt_table(const std::string& path,
                              const std::string& module,
                              const std::unordered_set<std::string>* known_haplotypes = nullptr,
                              const char* known_label = "",
                              bool quiet = false);

// Which contig of a slot stands for it where `observed` says which contigs were seen: exactly one. None
// seen, or two (overlapping contigs, which would count the same sequence twice), leaves the slot -- and so
// the sample -- unknown there, even where the two would agree (call --diploid keeps agreeing ones).
const std::string* slot_contig(const std::vector<std::string>& slot,
                               const std::function<bool(const std::string&)>& observed);

} // namespace panvar
