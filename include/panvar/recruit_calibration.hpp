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
    std::uint32_t min_sites = 0;          // 0: the index carries no calibration
    const Calibration* by = nullptr;      // the read class it came from
    bool met = true;                      // false: no cut reaches the rate (the largest is used)
};

// The read class nearest a run (short below 1 kb mean read length; long: hifi at an error rate up to 1%,
// else ont) and its lowest cut, never below 2, whose off-target false-certain share is at most `rate`.
Cut calibrated_cut(const std::vector<Calibration>& c, double mean_read_length, double error_rate, double rate);

// The share's one-sided 95% upper bound with `events` of `n` (Poisson).
double upper_bound(std::uint64_t events, std::uint64_t n);

}  // namespace panvar::recruit
