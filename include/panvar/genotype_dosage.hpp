#pragma once

// Dosage audit (docs/algorithms/genotype_dosage.md, sections 2-7): a graph census of copy-number-sensitive blocks,
// pooled copy-count features, a single-copy reference, a bootstrap copy-number estimate and the family-wise
// classification (policy v2). Truth-free: the inputs are the prepared model, the graph's node sequences, the call the
// caller reports, the reads and the library. Nothing here changes a call.

#include "panvar/genotype_model.hpp"
#include "panvar/mosaic_spelling.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace panvar::dosage {

// Frozen constants of the specification. They are named here, once, so a changed value is a visible design change.
struct DosageSettings {
    std::size_t kmer = 31;
    std::size_t min_feature_kmers = 500;
    std::size_t max_novel_bp = 50;             // census: structural threshold, declared apart from the molecule span
    double molecule_span_sigmas = 4.0;         // census: span = ceil(insert mean + 4 insert sd), short reads only
    std::size_t bootstrap_replicates = 200;
    std::uint32_t bootstrap_seed = 20260925;
    double alpha = 0.001;
    double suspect_z = 2.5;
    double power = 0.90;
    double iqr_to_sd = 1.349;
};

// The Python `random.Random(seed)` stream (MT19937, init_by_array seeding, randrange by getrandbits rejection). The
// bootstrap is frozen to this stream so the port reproduces the sealed prototype bit for bit.
class PythonRandom {
public:
    explicit PythonRandom(std::uint32_t seed);
    std::uint32_t next_uint32();
    // random.randrange(n) for 0 < n < 2^32
    std::uint32_t randrange(std::uint32_t n);

private:
    std::array<std::uint32_t, 624> state_{};
    std::size_t index_ = 624;
};

// statistics.NormalDist().inv_cdf (Wichura AS241), the function the prototype used.
double normal_inverse_cdf(double p);

// Molecule span of a short-read library: ceil(insert_mean + 4 * insert_sd). Refuses a library that declares no
// insert distribution (returns nullopt), so a long-read model never silently inherits this formula.
std::optional<std::size_t> molecule_span_bp(double insert_mean_bp, double insert_sd_bp, const DosageSettings& settings);

// 2-bit packed canonical k-mers (A<C<G<T, so numeric order is the prototype's lexicographic order). A k-mer with a
// base outside ACGT is skipped: it can never match. Positions are start offsets.
using Kmer = std::uint64_t;
struct KmerAt {
    Kmer kmer;
    std::size_t position;
};
std::vector<KmerAt> canonical_kmers(const std::string& sequence, std::size_t k);

// One allele's primary walk, as the locus builder recorded it in the allele's first `walk:` source.
struct NodeStep {
    std::string node;
    bool forward = true;
};
using NodeSequences = std::unordered_map<std::string, std::string>;   // node id -> forward sequence
std::vector<NodeStep> primary_walk(const CanonicalBlockAllele& allele);
std::vector<std::vector<NodeStep>> all_walks(const CanonicalBlockAllele& allele);

enum class CensusRoute { NotSensitive, Pair, RepeatedUnit };
const char* census_route_name(CensusRoute route);

struct CensusBlock {
    std::size_t block = 0;
    std::string block_id;
    std::size_t alleles = 0;
    std::size_t max_node_multiplicity = 0;
    std::size_t dosage_only_pairs = 0;
    std::size_t max_repeated_bp = 0;
    CensusRoute route = CensusRoute::NotSensitive;
    bool sensitive() const noexcept { return route != CensusRoute::NotSensitive; }
};

std::vector<CensusBlock> run_census(const PreparedGenotypeModel& model, const NodeSequences& nodes,
                                    std::size_t molecule_span_bp, const DosageSettings& settings);

// Pooled features of one block (graph only): groups of varying k-mers sharing one copy-count column, >= min k-mers,
// ordered by first occurrence (alleles by ascending ID, positions left to right).
struct PooledFeature {
    std::vector<std::uint16_t> copies;   // per allele ID
    std::vector<Kmer> kmers;             // sorted
};
struct BlockFeatureModel {
    std::vector<PooledFeature> features;
    std::vector<Kmer> block_kmers;       // every k-mer of every panel allele of the block, sorted
};
BlockFeatureModel pooled_features(const PreparedGenotypeBlock& block, const DosageSettings& settings);

// The call the caller reports: its two homologues, their block intervals and block alleles.
struct ReportedCall {
    std::array<std::string, 2> haplotypes;
    std::array<std::vector<SequenceInterval>, 2> intervals;
    std::array<std::vector<std::size_t>, 2> alleles;
};
ReportedCall reported_call(const PreparedGenotypeModel& model, const MosaicPath& first, const MosaicPath& second);

enum class DosageStatus { Consistent, Repaired, TriggeredUnrepaired, AutoRepair, Suspect, Underpowered, NotSensitive,
                          NotEvaluated, ContentRepaired };
const char* dosage_status_name(DosageStatus status);

struct FeatureAudit {
    std::size_t feature = 0;
    std::size_t pooled_index = 0;             // index into the block's pooled features (graph)
    std::size_t kmers = 0;
    std::vector<std::uint16_t> copy_values;   // distinct values of the copy vector, ascending
    int call_cn = 0;
    int feasible_min = 0, feasible_max = 0;
    double cn_estimate = 0.0;
    double bootstrap_sd = 0.0;
    std::array<double, 5> quantiles{};        // 0.05 0.25 0.50 0.75 0.95
    std::vector<int> admissible;
    bool call_admissible = false;
    double robust_sd = 0.0;
    double robust_z = 0.0;
    bool trigger = false, suspect = false, underpowered = false;
};

struct BlockAudit {
    std::size_t block = 0;
    std::string block_id;
    CensusRoute route = CensusRoute::NotSensitive;
    DosageStatus status = DosageStatus::NotSensitive;
    std::string reason = "none";
    std::size_t reference_kmers = 0;
    std::uint64_t reference_hits = 0;
    std::size_t retained_pairs = 0;
    std::array<std::size_t, 2> call_alleles{};
    std::vector<FeatureAudit> features;
    double max_robust_z = 0.0;
};

struct DosageAudit {
    std::size_t molecule_span_bp = 0;
    std::size_t m = 0;
    std::optional<double> z_critical;
    double alpha = 0.0;
    std::vector<BlockAudit> blocks;
};

DosageAudit audit_dosage(const PreparedGenotypeModel& model, const std::vector<CensusBlock>& census,
                         const std::vector<BlockFeatureModel>& features, const ReportedCall& call,
                         const std::vector<ReadPair>& reads, std::size_t molecule_span_bp,
                         const DosageSettings& settings, std::size_t threads);

std::string dosage_features_table(const DosageAudit& audit);
std::string dosage_blocks_table(const DosageAudit& audit);

} // namespace panvar::dosage
