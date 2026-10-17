#pragma once

// I3-A: a prepared locus -- the sample-independent preparation of `panvar genotype` (the structural model and the
// prepared marker panel), built once and analysed against by successive samples, sequentially.
//
// What is prepared, and why exactly this (tests/results/i3a_prepared_locus_20261016/DEPENDENCIES.md):
//   stage 1  graph read, reference sort and resolution, parse and path validation, bubbles (computed from the graph,
//            min_variant_bp and the Cactus snarls, or read from a supplied CSV), locus build with the curation
//            manifest, model preparation, the panel round trip, the representation-independence check, the model
//            fingerprint and the graph's node sequences;
//   stage 2  the marker panel: the model with kmer_size and syncmer_s only.
// Nothing that reads the reads or the library is prepared: library estimation, marker observations, depth, proposals,
// the union index (its haplotypes are the sample's candidates and their junctions), placements, the adaptive search,
// and the dosage census (it depends on the molecule span, i.e. on the library).
//
// The object is immutable after construction and owns all of its storage; it holds no pointer into a graph, a staged
// file or a sample. Its compatibility key is computed from input CONTENT (SHA-256) and the preparation settings; a
// sample whose options give a different key is refused (PreparedLocusMismatch), never silently analysed against it.
// Sequential reuse only: one sample at a time.

#include "panvar/genotype_command.hpp"
#include "panvar/fragment_likelihood.hpp"
#include "panvar/genotype_locus.hpp"
#include "panvar/genotype_model.hpp"
#include "panvar/marker_model.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace panvar {

struct PreparedLocusKey {
    std::string gfa_sha256;
    std::string reference_path;        // as requested; the prepared object records what it resolved to
    std::string bubbles;               // "computed" or the supplied CSV's SHA-256
    std::size_t min_variant_bp = 0;    // read only when bubbles are computed (recorded either way)
    std::string curation;              // "none" or the manifest's SHA-256
    std::size_t kmer_size = 0, syncmer_s = 0;
    // The first differing field ("" when equal), named, for the refusal message.
    std::string mismatch(const PreparedLocusKey& other) const;
    std::string describe() const;
};

class PreparedLocusMismatch : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

class PreparedLocus {
public:
    struct Parts;
    explicit PreparedLocus(Parts&& parts);
    ~PreparedLocus();
    PreparedLocus(const PreparedLocus&) = delete;
    PreparedLocus& operator=(const PreparedLocus&) = delete;

    const PreparedLocusKey& key() const;
    const PreparedGenotypeModel& model() const;
    const std::string& fingerprint() const;
    const LocusInputReport& report() const;
    const std::string& resolved_reference() const;
    const std::unordered_map<std::string, std::string>& node_sequences() const;
    const PreparedMarkerPanel& panel() const;
    // One-time preparation timings (seconds_detail-style names) and this object's process-unique id.
    const std::vector<std::pair<std::string, double>>& preparation_seconds() const;
    std::uint64_t id() const;
    // I3-C, runtime only: "prepared" (built in this process) or "loaded" (from a persistent index); for a loaded one, the
    // index's header SHA-256 and its stored provenance (described), else empty.
    const std::string& origin() const;
    const std::string& file_checksum() const;
    const std::string& stored_provenance() const;
    // I3-D, runtime only: "prepared" (built here), "loaded" (from an index's panel section) or "rebuilt" (from an index
    // without one, or with one this binary does not support).
    const std::string& panel_source() const;
    // A SHA-256 over SELECTED content -- the model fingerprint (recomputed), a few report counts, the node sequences
    // and every marker-panel field -- for immutability checks within a process. LIMITATION (I3-A.1): it is NOT a
    // complete field-by-field integrity digest. It covers the model only through model_fingerprint (which does not hash,
    // for example, every allele's source walk ids), and only some LocusInputReport fields (counts, not the per-block,
    // coverage, partial-template or curation records). A persistent format must not rely on it: it needs a canonical,
    // complete payload-integrity scheme and a schema/build compatibility policy.
    std::string digest() const;

private:
    std::unique_ptr<Parts> parts_;
};

// The compatibility key a command's options give (hashes the GFA, the bubbles CSV and the curation manifest).
PreparedLocusKey prepared_locus_key(const GenotypeCommandOptions& options);
// Completed preparations in this process (tests prove that a batch prepares once). Loads are not preparations.
std::uint64_t prepared_locus_builds();
// I3-C: completed loads from a persistent index in this process.
std::uint64_t prepared_locus_loads();

// Prepares from a `panvar genotype` command line, resolved exactly as the command resolves it (reference sample,
// library default). Throws when the structural model cannot be built (the ordinary command would write
// MODEL_INCOMPLETE for that sample instead).
std::shared_ptr<const PreparedLocus> prepare_genotype_locus(const std::vector<std::string>& args);
std::shared_ptr<const PreparedLocus> prepare_genotype_locus(const GenotypeCommandOptions& options);

// ---- I3-A.1: successive samples in one process ---------------------------------------------------------------------------
// Between two FINISHED samples (their workers and watchdog joined): release the process-global storage the sample left
// (release_sample_retained_storage) and ask the allocator to return free pages (release_free_pages, macOS only; 0
// elsewhere). The prepared locus is untouched.
struct BetweenSampleCleanup {
    double seconds = 0.0;
    RetainedProcessStorage released;      // what the process-global containers held before the release
    std::uint64_t free_pages_released = 0;  // as reported by the allocator (0 where unsupported)
};
BetweenSampleCleanup release_between_samples();

// The batch's resource verdict, from the process-lifetime peak read at every phase boundary (preparation, each sample,
// each cleanup, teardown). The peak cannot be reset, so attribution is by phase: the phase during which the lifetime peak
// first rose above the ceiling is the one between the last boundary at or below it and the first above it. This is
// DETECTION after an excursion, not prevention: no mechanism here stops an allocation that crosses the ceiling, and the
// per-run watchdog is sampled (100 ms) and can miss brief excursions. Once exceeded the verdict cannot recover, and the
// caller must not start another sample. I3-B: availability is explicit at EVERY boundary: a reading that is missing at
// any boundary (std::nullopt -- never read as zero) makes the batch UNVERIFIED unless a breach is observed; an observed
// breach is CEILING_EXCEEDED permanently, whatever other boundaries were unverified.
enum class BatchResourceVerdict { compliant, ceiling_exceeded, unverified };
const char* batch_resource_verdict_name(BatchResourceVerdict verdict);

class BatchResourceLedger {
public:
    BatchResourceLedger(std::uint64_t ceiling_bytes, std::string metric);
    // Records a phase boundary; either reading may be unavailable (nullopt). Returns whether another sample may start
    // (false only after an observed breach; an unverified boundary does not stop the batch, it only forbids COMPLIANT).
    bool record(const std::string& phase, std::optional<std::uint64_t> current_bytes,
                std::optional<std::uint64_t> lifetime_peak_bytes);
    BatchResourceVerdict verdict() const;
    const std::string& exceeded_during() const { return exceeded_during_; }   // "" unless exceeded
    std::uint64_t ceiling() const { return ceiling_; }
    // The largest lifetime peak actually READ (unavailable boundaries contribute nothing); nullopt if none was read.
    std::optional<std::uint64_t> max_peak() const { return max_peak_; }
    const std::string& metric() const { return metric_; }
    std::size_t unverified_boundaries() const { return unverified_; }
    struct Phase {
        std::string name;
        std::optional<std::uint64_t> current, peak;
        // The rise of the lifetime peak over the largest peak READ at any earlier boundary; nullopt when this
        // boundary's peak is unavailable (a rise across unverified boundaries is attributed to the next read one).
        std::optional<std::uint64_t> peak_increase;
        bool verified = false;   // both readings available
    };
    const std::vector<Phase>& phases() const { return phases_; }

private:
    std::uint64_t ceiling_ = 0;
    std::optional<std::uint64_t> max_peak_;
    std::string metric_;
    bool exceeded_ = false;
    std::size_t unverified_ = 0;
    std::string exceeded_during_;
    std::vector<Phase> phases_;
};

// The ordinary per-sample command (genotype, dosage, locus table) against a prepared locus: the same implementation
// as run_genotype_command, which builds the model and panel inline. Throws PreparedLocusMismatch, before writing any
// output, when the sample's options are not compatible with the object.
int run_genotype_command_prepared(const std::vector<std::string>& args, const PreparedLocus& locus);
GenotypeRunStatus run_genotype_prepared(const GenotypeCommandOptions& options, const PreparedLocus& locus,
                                        GenotypeDosageCapture* capture = nullptr);

}  // namespace panvar
