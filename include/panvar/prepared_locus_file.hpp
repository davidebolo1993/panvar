#pragma once

// I3-C (EXPERIMENTAL, internal; no CLI): the persistent prepared-locus index, v1. Design and its corrections:
// tests/analysis/i3c/DESIGN.md.
//
// One policy for the structural source: a PreparedLocus never owns the StructuralLocusInput, prepared or loaded. The
// complete persistent content lives in a PreparedLocusRecord, a transient value produced beside a preparation (to be
// written) or by decoding a file (to be loaded or rewritten). Loading releases the read buffer after decoding and the
// structural input once the reconstructed model's digest verifies, before the panel is built.
//
// The file is canonical: the same inputs and the same build give byte-identical files. No runtime id, timing, path or
// timestamp is written. Provenance (a deterministic build description) is recorded, never used for compatibility, and
// preserved byte for byte when an unchanged index is rewritten. Compatibility is the schema version and the
// preparation_algorithm_version; identity is the input-content key, checked against a sample's supplied inputs by the
// existing sample API (which still hashes the original inputs: this is not yet an index-only invocation).

#include "panvar/genotype_locus.hpp"
#include "panvar/genotype_model.hpp"
#include "panvar/marker_model.hpp"
#include "panvar/prepared_locus.hpp"

#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace panvar {

inline constexpr std::uint32_t kPreparedLocusSchemaVersion = 1;
// Raised whenever the semantics of graph sort/flip, parsing, bubble calling, build_locus_input, prepare_genotype_model,
// model_fingerprint, build_marker_panel, the content digests or the persisted field set change.
inline constexpr std::uint32_t kPreparationAlgorithmVersion = 1;
// I3-D: the optional PANEL_MARKER_TABLE section (id 16; tests/analysis/i3d/DESIGN.md). Its own layout version and the
// version of the build_marker_panel semantics that produced the stored panel; separate from the base versions above.
inline constexpr std::uint32_t kPanelExtensionVersion = 1;
inline constexpr std::uint32_t kMarkerPanelAlgorithmVersion = 1;

// Producer provenance (M15-K): who wrote an index -- the writer, the source revision actually built (source_commit: a git
// commit or "unavailable"), whether its tracked files matched that commit (source_tree: clean | dirty | unknown), the
// compiler and the build type. Reported, never compared when an index is loaded: compatibility is the identity and the
// versions, checked separately.
struct PreparedLocusProvenance {
    std::string writer, source_commit, source_tree, compiler, build_type;
    bool operator==(const PreparedLocusProvenance& o) const {
        return writer == o.writer && source_commit == o.source_commit && source_tree == o.source_tree &&
               compiler == o.compiler && build_type == o.build_type;
    }
};
// This binary's provenance (deterministic for a given build).
PreparedLocusProvenance current_prepared_locus_provenance();

struct PreparedLocusRecord {
    PreparedLocusKey key;
    std::string resolved_reference;
    std::uint32_t preparation_algorithm_version = kPreparationAlgorithmVersion;
    bool round_trip_passed = false;
    bool representation_independent = false;
    PreparedLocusProvenance provenance;
    StructuralLocusInput input;
    LocusInputReport report;
    std::vector<std::pair<std::string, std::string>> nodes;   // every GFA segment, ids strictly increasing
    std::string model_fingerprint, model_content_sha256, panel_content_sha256;
    // I3-D: the persisted panel (section 16), when present and understood; decoded straight into the value a loaded
    // PreparedLocus takes ownership of.
    std::optional<PreparedMarkerPanel> panel;
    // Optional sections this binary does not understand (an unknown id, or section 16 at an unsupported extension or
    // panel-algorithm version), kept verbatim so a rewrite preserves them; ids strictly increasing.
    std::vector<std::pair<std::uint32_t, std::string>> unknown_optional_sections;
};

struct PreparedLocusPreparation {
    std::shared_ptr<const PreparedLocus> locus;
    PreparedLocusRecord record;
};
// Prepares as prepare_genotype_locus does, and also returns the complete record (provenance = this build). I3-D:
// `persist_panel` also copies the prepared panel into the record, so the written index carries section 16.
PreparedLocusPreparation prepare_genotype_locus_record(const GenotypeCommandOptions& options, bool persist_panel = false);
PreparedLocusPreparation prepare_genotype_locus_record(const std::vector<std::string>& args, bool persist_panel = false);

// Complete content digests (streaming SHA-256 over a canonical field-by-field encoding of every field).
std::string prepared_locus_model_content_sha256(const PreparedGenotypeModel& model);
std::string prepared_locus_panel_content_sha256(const PreparedMarkerPanel& panel);

class PreparedLocusFormatError : public std::runtime_error {
public:
    enum class Kind {
        io, truncated, bad_magic, unsupported_schema, invalid_header, checksum_mismatch, invalid_length, invalid_count,
        overflow, invalid_enum, invalid_section, unknown_required_section, invalid_structure, incompatible_preparation,
        limit_exceeded, reconstruction_mismatch, exists,
        // I3-D: the target WAS published (complete, linked after fsync) and is left in place; a later step (removing the
        // temporary, or opening/fsyncing the directory) failed, so a temporary may remain or the entry's durability is
        // unconfirmed.
        post_publication
    };
    PreparedLocusFormatError(Kind kind, const std::string& message);
    Kind kind() const { return kind_; }
private:
    Kind kind_;
};
const char* prepared_locus_error_name(PreparedLocusFormatError::Kind kind);

// The canonical file bytes (header + payload) of a record.
std::string encode_prepared_locus(const PreparedLocusRecord& record);
// Atomic publication: a temporary opened O_CREAT|O_EXCL, written, fsynced, linked to `path` (never overwrites: an
// existing `path` refuses with `exists` and is untouched), the temporary unlinked and the directory fsynced. A failure
// BEFORE the link removes the temporary and leaves `path` uncreated. A failure AFTER the link (unlinking the temporary,
// or opening/fsyncing the directory) throws `post_publication`: `path` then already holds the complete, valid index and
// is left in place (never deleted), and a temporary may remain or the entry's durability is unconfirmed.
void publish_prepared_locus(const std::string& file_bytes, const std::string& path);
void write_prepared_locus(const PreparedLocusRecord& record, const std::string& path);

struct PreparedLocusLoadLimits {
    std::uint64_t max_file_bytes = 1ull << 30;
    std::uint64_t max_elements = 1ull << 28;            // any single decoded count
    std::uint64_t max_accounted_bytes = 7ull << 29;     // 3.5 GiB, cumulative over every load phase (DESIGN D)
};

// Runtime only (never written): what a read or load did.
struct PreparedLocusLoadReport {
    std::vector<std::pair<std::string, double>> seconds;     // per phase
    std::uint64_t file_bytes = 0;
    std::string file_checksum;                               // the header SHA-256, hex
    std::uint64_t accounted_peak_bytes = 0, accounted_final_bytes = 0, accounting_limit = 0;
    std::uint64_t model_projection_bytes = 0, panel_projection_bytes = 0;
    std::uint64_t panel_syncmer_occurrences = 0;
    std::uint64_t structural_input_accounted_bytes = 0;     // what the release of the structural input returned
    std::string peak_phase;
    // I3-D: where the panel came from ("loaded" from section 16, or "rebuilt" through the v1 path) and why it was rebuilt
    std::string panel_source, panel_rebuild_reason;
    std::uint64_t panel_decoded_accounted_bytes = 0;
    // consecutive charges under one label form one entry: (label, the highest accounted total reached under it)
    std::vector<std::pair<std::string, std::uint64_t>> phase_peaks;
};

// Container, checksum, bounded decode and the structural invariants (DESIGN A.1, C, D phases 1-3). No reconstruction.
PreparedLocusRecord read_prepared_locus_record(const std::string& path, const PreparedLocusLoadLimits& limits,
                                               PreparedLocusLoadReport* report = nullptr);
// The full load: read_prepared_locus_record, then reconstruction through the existing constructors, verification of
// the fingerprint and both content digests, and the release of the structural source (DESIGN B, D phases 4-8). Refuses
// with PreparedLocusFormatError; never rebuilds.
std::shared_ptr<const PreparedLocus> load_prepared_locus(const std::string& path, const PreparedLocusLoadLimits& limits,
                                                         PreparedLocusLoadReport* report = nullptr);

// I3-D: executions of the panel-reconstruction pre-count (the v1 route's collect_syncmers pass) in this process.
std::uint64_t prepared_locus_panel_precounts();
// I3-D: "" when every PreparedMarkerPanel field is equal, else the first differing field (tests, write verification).
std::string prepared_marker_panel_difference(const PreparedMarkerPanel& a, const PreparedMarkerPanel& b);

// The decode stage alone over in-memory file bytes (tests: corruption with a recomputed checksum, accounting controls).
PreparedLocusRecord decode_prepared_locus(const std::string& file_bytes, const PreparedLocusLoadLimits& limits,
                                          PreparedLocusLoadReport* report = nullptr);
// Recomputes the header checksum of `file_bytes` in place (tests only: semantic corruption behind a valid checksum).
void reseal_prepared_locus(std::string& file_bytes);

}  // namespace panvar
