#include "panvar/copying_transitions.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>

#include "panvar/genotype_repair.hpp"   // sha256_bytes

namespace panvar {

namespace {
std::string digest_of(const std::string& mode, const std::vector<double>& r) {
    std::string bytes = mode;
    for (const double v : r) {
        char buffer[32];
        std::snprintf(buffer, sizeof buffer, "|%a", v);
        bytes += buffer;
    }
    return sha256_bytes(bytes);
}
}  // namespace

CopyingTransitions build_copying_transitions(std::size_t blocks, const LocusGeometry* geometry, const CopyingTransitionSpec& spec) {
    if (blocks < 1) throw UnsupportedCopyingTransition("copying transitions need at least one block");
    CopyingTransitions t;
    t.mode = spec.mode;
    if (spec.mode == "constant_v1") {
        if (!std::isfinite(spec.constant_r) || !(spec.constant_r > 0.0) || spec.constant_r > 0.5) {
            throw UnsupportedCopyingTransition("constant_v1: switch_probability must be in (0, 0.5]");
        }
        t.constant_r = spec.constant_r;
        t.r.assign(blocks - 1, spec.constant_r);
    } else if (spec.mode == "distance_ref_core_v1_experimental") {
        if (!std::isfinite(spec.copying_length_bp) || !(spec.copying_length_bp > 0.0)) {
            throw UnsupportedCopyingTransition("distance_ref_core_v1_experimental: copying_length_bp must be finite and > 0");
        }
        if (geometry == nullptr) throw UnsupportedCopyingTransition("distance_ref_core_v1_experimental: no locus geometry was supplied");
        if (!geometry->available()) throw UnsupportedCopyingTransition("distance_ref_core_v1_experimental: locus geometry is " + geometry->status);
        if (geometry->block_ids.size() != blocks) {
            throw UnsupportedCopyingTransition("distance_ref_core_v1_experimental: the geometry describes " + std::to_string(geometry->block_ids.size()) +
                                               " blocks, the model " + std::to_string(blocks));
        }
        t.copying_length_bp = spec.copying_length_bp;
        t.geometry_digest = geometry->digest;
        t.distance_bp = geometry->distances();
        for (std::size_t b = 0; b < t.distance_bp.size(); ++b) {
            const double d = t.distance_bp[b];
            if (!std::isfinite(d) || d < 0.0) throw UnsupportedCopyingTransition("distance_ref_core_v1_experimental: a negative or non-finite distance");
            const double r = d > 0.0 ? -std::expm1(-d / spec.copying_length_bp) : 0.0;
            t.r.push_back(std::min(1.0, std::max(0.0, r)));
            if (d == 0.0) {
                ++t.zero_distance_boundaries;
                if (geometry->frame_distinct_alleles[b] > 1 && geometry->frame_distinct_alleles[b + 1] > 1) ++t.zero_distance_informative_boundaries;
            }
        }
    } else {
        throw UnsupportedCopyingTransition("copying_transition '" + spec.mode + "' is not constant_v1 or distance_ref_core_v1_experimental");
    }
    t.vector_digest = digest_of(t.mode, t.r);
    return t;
}

double hazard_matched_copying_length(const LocusGeometry& geometry, double constant_r) {
    if (!(constant_r > 0.0) || constant_r >= 1.0) throw UnsupportedCopyingTransition("hazard-matched length needs 0 < r < 1");
    const std::vector<double> d = geometry.distances();
    if (d.empty()) throw UnsupportedCopyingTransition("hazard-matched length needs at least one boundary");
    double total = 0.0;
    for (const double x : d) total += x;
    if (!(total > 0.0)) throw UnsupportedCopyingTransition("hazard-matched length needs a positive total distance");
    return total / (static_cast<double>(d.size()) * -std::log1p(-constant_r));
}

}  // namespace panvar
