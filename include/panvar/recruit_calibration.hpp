#pragma once

// The calibration of a locus index (section Calibration in locus_index.hpp) and what it gives a run: the score's
// thresholds, the lowest at which the off-target fragments reach a tier at no more than the asked rate; and, for
// the first rules, the lowest private-site cut that does.

#include "panvar/gfa.hpp"
#include "panvar/locus_index.hpp"

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

// The calibration genomes chosen per locus by default, at most.
constexpr std::size_t kCalibrationGenomesMax = 8;
// The bases of a calibration genome read on each side of each of its copies of the locus.
constexpr std::uint64_t kCalibrationGenomeFlank = 1000000;

// A calibration genome's off-target sequence for one locus, and the labels that change while it is held out.
struct GenomeOffTarget {
    std::size_t genome = 0;                   // index into the locus's calibration_genomes
    std::vector<OffTarget> off;               // its homologous windows with context, and its flanks
    std::vector<std::uint32_t> relabel;       // locus k-mers (indices) it alone holds outside the locus: private
    std::vector<std::uint64_t> unlabel;       // decoy k-mers it alone holds: unknown
};

void calibrate(LocusIndex& x, const Graph& g, const std::vector<OffTarget>& off, const std::vector<GenomeOffTarget>& genomes,
               std::size_t threads);

// The calibration of the run's read class (`read_class`, or when empty the nearest: short below 1 kb mean read
// length; long: hifi at an error rate up to 1%, else ont); null when the index has none.
const Calibration* calibration_for(const std::vector<Calibration>& c, const std::string& read_class, double mean_read_length,
                                   double error_rate);

// One tier's threshold: the lowest score bin at which the off-target fragments reaching it are at most `rate` of
// those screened.
struct Tier {
    std::int32_t bin = 0;                     // fragments whose score bin is at least this reach the tier
    double rate = 0;                          // asked
    std::uint64_t off = 0, off_screened = 0;  // off-target fragments reaching it, of those screened
    std::uint64_t locus = 0, locus_total = 0; // held-out fragments reaching it, of those simulated
    bool shown = false;                       // the off-target share's 95% bound is within the rate
    double bound() const;
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
