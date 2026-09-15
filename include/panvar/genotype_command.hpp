#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace panvar {

// The only public surface of the mosaic caller. Debug probes live in unit tests or a developer
// executable, not here: `genotype-frag`, `--hybrid-call` and `--mosaic-rescore` are deliberately
// gone and must not come back.
struct GenotypeCommandOptions {
    std::string gfa_path;
    // Optional. A bubbles CSV from `panvar bubble`; when absent the bubbles are called internally
    // with the same SV threshold, so the block alphabet is the same either way.
    std::string bubbles_path;
    std::string reads1_path;
    std::string reads2_path;
    std::string reference_path;
    std::string out_prefix = "genotype";
    // Declared library parameters. Version one does not estimate the insert distribution from the
    // data, so this is required and its provenance is recorded in the status file.
    std::string library_path;
    std::size_t threads = 1;
    std::size_t max_candidates = 64;

    std::size_t kmer_size = 31;
    std::size_t syncmer_s = 15;
    std::size_t min_variant_bp = 50;
    double switch_probability = 0.01;
    // The declared mismatch band of the truncated fragment likelihood. It is part of the model, not a
    // search depth, so it is fixed for the whole run and recorded in the status file.
    std::size_t mismatch_band = 6;
    std::size_t minimum_piece_bp = 12;
    bool quiet = false;
};

// A run that cannot complete emits its diagnostics and says so. It never presents a partial ranking
// as a result.
enum class GenotypeRunStatus {
    Complete,
    ModelIncomplete,
    CandidateIncomplete,
    SpellingIncomplete,
    RescoreIncomplete,
};

const char* genotype_run_status_name(GenotypeRunStatus status);

GenotypeRunStatus run_genotype(const GenotypeCommandOptions& options);

int run_genotype_command(const std::vector<std::string>& args);

} // namespace panvar
