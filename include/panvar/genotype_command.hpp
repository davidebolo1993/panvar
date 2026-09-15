#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
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
    // Declared resource budgets. Reaching one ends the run incomplete with no ranking and no call; the
    // values are recorded in the status file. Not command-line options: the public surface is fixed.
    std::uint64_t max_beam_extensions = 2000000000ULL;
    // Merged proposal states retained per block: a resource cap on the grouped-state search. A search
    // that prunes and still misses max_candidates ends CANDIDATE_INCOMPLETE, whatever this is set to.
    std::size_t proposal_state_limit = 1024;
    // Panel pairs times fragments: the predicted diploid-combination work of the protected panel baseline.
    // A baseline over this is refused whole, never scored partially.
    std::uint64_t max_panel_pair_fragment_evaluations = 20000000000ULL;
    std::uint64_t max_piece_hits = 2000000000ULL;
    std::uint64_t max_proposed_starts_per_haplotype = 50000000ULL;
    std::uint64_t max_physical_origins = 500000000ULL;
    std::uint64_t max_state_classes = 200000000ULL;
    std::uint64_t max_mass_cells = 200000000ULL;
    bool quiet = false;
    // Development harness only; not a command-line option. When non-empty (strictly increasing, the
    // last equal to max_candidates), the model, markers and proposal are computed once, and each prefix
    // k is then refused, ranked and written under <out_prefix>.K<k> by the same code a standalone run
    // with max_candidates = k uses: K only truncates the final proposal ordering. Every candidate a
    // prefix contains is scored; the production refusal rule still applies to each k.
    std::vector<std::size_t> development_candidate_prefixes;
    // Not a command-line option. When positive (and larger than every primary candidate set), the same proposal
    // and scores are also ranked with this many mosaics, and read-supported phase is reported READS_AUDITED only
    // if that larger set agrees on the winner, the block genotypes and read-supported phase. Without an audit a
    // read-supported phase is READS_WITHIN_CANDIDATES.
    std::size_t phase_audit_candidates = 0;
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

// The phase claim after an audit. A read-supported phase within the primary candidates becomes READS_AUDITED
// only when the audit set was ranked to completion, agrees on the winner and on every block genotype, and is
// itself read-supported; with the same call but a phase tie there it is PRIOR_ONLY; otherwise UNRESOLVED. Any
// other primary claim is unchanged. Exposed for its gate.
std::string audited_phase_support(const std::string& within_candidates, bool audit_complete, bool same_winner,
                                  bool same_block_genotypes, const std::string& audit_phase_support);

// Per-block genotype confidence within a retained candidate set. Candidate 0 is the winner; scores are total
// log scores and error_bounds their floating-point bounds; genotypes[c][b] is candidate c's unordered allele
// pair at block b. The posterior is the set mass of candidates sharing the winner's genotype there; the quality
// is -10 log10 of the rest, summed in log space so a posterior indistinguishable from one keeps a finite GQ;
// numerically_separated says every candidate with another genotype is below the winner beyond error.
// Exposed for its gate.
struct BlockGenotypeConfidence {
    double posterior = 1.0;
    double quality = 0.0;
    bool numerically_separated = true;
};
std::vector<BlockGenotypeConfidence> block_genotype_confidence(
    const std::vector<double>& scores,
    const std::vector<double>& error_bounds,
    const std::vector<std::vector<std::pair<std::size_t, std::size_t>>>& genotypes);

GenotypeRunStatus run_genotype(const GenotypeCommandOptions& options);

int run_genotype_command(const std::vector<std::string>& args);

} // namespace panvar
