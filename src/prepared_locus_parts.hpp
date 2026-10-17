#pragma once

// Internal to panvarlib: the storage of a PreparedLocus, shared by its preparation (genotype_command.cpp) and its
// loading from a persistent index (prepared_locus_file.cpp, I3-C). Not installed, not part of any API.

#include "panvar/prepared_locus.hpp"

#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace panvar {

struct PreparedLocus::Parts {
    PreparedLocusKey key;
    std::optional<PreparedGenotypeModel> model;
    std::string fingerprint, resolved;
    LocusInputReport report;
    std::unordered_map<std::string, std::string> node_sequences;
    std::optional<PreparedMarkerPanel> panel;
    std::vector<std::pair<std::string, double>> seconds;
    std::uint64_t id = 0;
    // I3-C, runtime only (never written to an index): how this object came to exist.
    std::string origin = "prepared";   // "prepared" or "loaded"
    std::string file_checksum;         // loaded: the index's header SHA-256
    std::string stored_provenance;     // loaded: the index's provenance, described
    std::string panel_source = "prepared";   // I3-D: "prepared", "loaded" or "rebuilt"
};

// A process-unique id for a new PreparedLocus (prepared or loaded).
std::uint64_t next_prepared_locus_id();

}  // namespace panvar
