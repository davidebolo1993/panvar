#pragma once

// The calibration of a locus index (section Calibration in locus_index.hpp) and what it gives a run: the score's
// thresholds, the lowest at which the off-target fragments reach a tier at no more than the asked rate; and, for
// the first rules, the lowest private-site cut that does.

#include "panvar/flat_map.hpp"
#include "panvar/gfa.hpp"
#include "panvar/locus_index.hpp"
#include "panvar/recruit_scoring.hpp"

#include <climits>
#include <cstdint>
#include <string>
#include <vector>

namespace panvar::recruit {

// An off-target stretch: a decoy window with context around it (reads overlap [core_start, core_end)),
// or a flank of the coordinate reference (all of it).
struct OffTarget {
    bool decoy = false;
    std::string seq;
    std::uint64_t core_start = 0, core_end = 0;
};

// The bases of context kept around a decoy window (long reads reach this far beyond it).
constexpr std::uint64_t kCalibrationContext = 15000;
// The coordinate reference's bases simulated on each side of the locus.
constexpr std::uint64_t kCalibrationFlank = 5000000;

// The bases of a calibration genome read on each side of each of its copies of the locus.
constexpr std::uint64_t kCalibrationGenomeFlank = 1000000;

// A stretch of a calibration genome is a copy of the locus, not off-target sequence, when over a window of this many
// k-mer starts this share are the locus's own (private once the genome is held out, or shared, not repeats): the
// density the index's pieces take. The index cannot place such a stretch (no anchors in order, not next to a
// copy), so its k-mers stay shared; a sample's reads from it are the locus's all the same.
constexpr std::uint64_t kCopyWindow = 5000;
constexpr double kCopyDensity = 0.95;

// A stretch of a calibration genome read as off-target sequence: a homologous window with context around it
// (fragments overlap [core_start, core_end)) or a flank of one of its copies of the locus (all of it). Offsets
// into one of the genome's blocks, which are disjoint, so no base is counted twice.
struct OffStretch {
    std::size_t block = 0;
    std::uint64_t start = 0, end = 0, core_start = 0, core_end = 0;
    bool window = false;
};

// Which codes may be k-mers or decoys of a locus: one bit per hashed code, so that most k-mers of a genome, which
// are neither, cost one lookup.
class DecoyFilter {
public:
    DecoyFilter(const std::vector<std::uint64_t>& codes, const std::vector<std::uint64_t>& decoys);
    bool maybe(std::uint64_t code) const;
private:
    std::vector<std::uint64_t> bits_;
    std::uint64_t mask_ = 0;
};

// What a calibration genome's off-target sequence gives one read class while the genome is held out of the
// background: the locus k-mers only it holds outside the locus read as private, the decoy k-mers only it holds as
// unknown. Fragments are laid over its stretches at a fixed stride; one holding no private k-mer cannot be called
// and is only counted, the others are scored from the labels with the class's errors, each error leaving the
// k-mers over it unknown.
struct GenomeClassScan {
    std::uint64_t fragments = 0, screened = 0, with_private = 0;
    std::uint64_t copies = 0;                      // fragments left out: over a copy of the locus
    double sites_p = 0, sites_d = 0, kmers = 0;    // over the screened fragments
    std::vector<scoring::Judged> judged;           // the screened fragments that can be called
};
struct GenomeScan {
    CalibrationGenome stats;
    std::vector<std::uint64_t> homolog;            // a bit per locus k-mer, then per decoy: the ones it holds
    std::vector<GenomeClassScan> classes;          // by read class, in the calibration's order
};
GenomeScan scan_genome(const LocusIndex& x, const FlatMap<std::uint32_t>& code_index, const DecoyFilter& filter,
                       const std::vector<std::uint64_t>& decoy_counts, const std::string& genome,
                       const std::vector<std::string>& blocks, const std::vector<OffStretch>& stretches);

// The sample a path belongs to (the PanSN name's first field); `owner` below numbers samples by their first path.
std::string calibration_sample(const std::string& path_name);

// Calibrates a locus index: simulated reads of held-out samples (`owner`: per k-mer, the one sample holding it, or
// -1) and of the off-target stretches, and the scans of the calibration genomes chosen for it.
void calibrate(LocusIndex& x, const Graph& g, const std::vector<OffTarget>& off, const std::vector<GenomeScan>& genomes,
               const std::vector<std::int32_t>& owner, std::size_t threads);

// The calibration of the run's read class (`read_class`, or when empty the nearest: short below 1 kb mean read
// length; long: hifi at an error rate up to 1%, else ont); null when the index has none.
const Calibration* calibration_for(const std::vector<Calibration>& c, const std::string& read_class, double mean_read_length,
                                   double error_rate);

// One tier's threshold: the lowest score bin at which the off-target fragments reaching it are at most `rate` of
// those screened. Two kinds of off-target fragment add up: simulated ones of the indexed genomes, which gain
// private sites by sequencing errors, and scanned ones of the held-out genomes, which hold them as a new sample
// would; each kind's share is of its own screened fragments. When the samples are too small to show the rate (no
// fragment beyond the threshold would still leave its 95% bound above the rate), the threshold is extrapolated:
// read off an exponential tail fitted to the off-target scores (see score_tier), never below the one observed.
struct Tier {
    std::int32_t bin = 0;                     // fragments whose score bin is at least this reach the tier
    double rate = 0;                          // asked
    std::uint64_t off = 0, off_screened = 0;  // off-target fragments reaching it, of those screened (both kinds)
    double share = 0;                         // the off-target share reaching it (the kinds' shares added)
    double share_bound = 1;                   // its 95% upper bound
    std::uint64_t locus = 0, locus_total = 0; // held-out fragments reaching it, of those simulated
    bool shown = false;                       // the bound is within the rate
    bool extrapolated = false;                // the threshold is the fitted tail's
    std::int32_t observed_bin = 0;            // the threshold the fragments seen give
    double tail_scale = 0;                    // the fitted tail: the score over which it falls e-fold,
    double tail_site_share = 0;               // the share of off-target fragments with a private site it comes from,
    std::int32_t tail_from = 0;               // and the bin of the score it is hung from
    double bound() const { return share_bound; }
};
constexpr std::int32_t kEveryBin = INT32_MIN;   // a threshold every callable fragment reaches
Tier score_tier(const Calibration& c, double rate);
// Whether a calibration carries the score (an index of this version with screened fragments of both kinds).
bool has_score(const Calibration& c);

struct Cut {
    std::uint32_t min_sites = 0;          // 0: no calibrated cut (no calibration of this class, or no screened negative)
    const Calibration* by = nullptr;      // the read class it came from
    bool met = true;                      // false: no cut reaches the rate (the largest is used)
    bool shown = false;                   // the false-certain share's 95% bound is within the rate: the
                                          // simulation was large enough to show it, not only to see none
    std::uint32_t floor = 2;              // the lowest cut the read class allows (lowest_cut)
};

// The lowest cut a read class takes: 2, and 3 for ont reads. Ultra-long noisy reads cross enough sequence to
// carry two chance private sites far more often than the simulated off-target reads show.
std::uint32_t lowest_cut(const std::string& read_class);

// The run's read class (`read_class`, or when empty the nearest: short below 1 kb mean read length; long:
// hifi at an error rate up to 1%, else ont) and its lowest cut, never below lowest_cut, whose observed off-target
// false-certain share is at most `rate`. With fewer screened negatives than 1 / rate this is the lowest cut
// with none observed.
Cut calibrated_cut(const std::vector<Calibration>& c, const std::string& read_class, double mean_read_length, double error_rate,
                   double rate);

// The share's one-sided 95% upper bound with `events` of `n` (Poisson).
double upper_bound(std::uint64_t events, std::uint64_t n);

}  // namespace panvar::recruit
