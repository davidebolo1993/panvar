#pragma once

// The calibration of a locus index (section Calibration in locus_index.hpp) and the private-site cut it
// gives a run: the lowest cut at which the off-target reads reach certain at no more than the asked rate.

#include "panvar/gfa.hpp"
#include "panvar/locus_index.hpp"

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

void calibrate(LocusIndex& x, const Graph& g, const std::vector<OffTarget>& off, std::size_t threads);

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

// The run's read class: `read_class`, or when empty the nearest (short below 1 kb mean read length; long: hifi at
// an error rate up to 1%, else ont).
std::string run_read_class(const std::string& read_class, double mean_read_length, double error_rate);

// The run's read class (run_read_class) and its lowest cut, never below lowest_cut, whose observed off-target
// false-certain share is at most `rate`. With fewer screened negatives than 1 / rate this is the lowest cut
// with none observed.
Cut calibrated_cut(const std::vector<Calibration>& c, const std::string& read_class, double mean_read_length, double error_rate,
                   double rate);

// The share's one-sided 95% upper bound with `events` of `n` (Poisson).
double upper_bound(std::uint64_t events, std::uint64_t n);

}  // namespace panvar::recruit
