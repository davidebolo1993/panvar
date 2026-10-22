// The combined per-block quality (block_quality_rule_v1): PASS only for representation PASS with READS support;
// otherwise LOW_QUALITY with every applying reason, in canonical order. A pure function of two call-time facts.
#include "panvar/genotype_command.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

namespace {
void expect(const std::string& filter, const std::string& support, const std::string& quality, const std::string& reason) {
    const panvar::BlockQuality got = panvar::block_quality(filter, support);
    if (got.filter != quality || got.reason != reason) {
        throw std::runtime_error("block_quality(" + filter + ", " + support + ") = " + got.filter + "/" + got.reason +
                                 ", expected " + quality + "/" + reason);
    }
    std::cout << "  ok: " << filter << " + " << support << " -> " << quality << " / " << reason << '\n';
}
}  // namespace

int main() {
    try {
        expect("PASS", "READS", "PASS", "none");
        expect("PASS", "PRIOR_ONLY", "LOW_QUALITY", "no_read_support");
        expect("LOWQUAL_REPRESENTATION_CONFLICT", "READS", "LOW_QUALITY", "representation_conflict");
        expect("NOT_EVALUATED", "READS", "LOW_QUALITY", "representation_not_evaluated");
        // Combined reasons, canonical order: representation first, then support.
        expect("LOWQUAL_REPRESENTATION_CONFLICT", "PRIOR_ONLY", "LOW_QUALITY", "representation_conflict,no_read_support");
        expect("NOT_EVALUATED", "PRIOR_ONLY", "LOW_QUALITY", "representation_not_evaluated,no_read_support");
        expect("NOT_EVALUATED", "NOCALL", "LOW_QUALITY", "representation_not_evaluated,no_read_support");
        // Fail closed: an absent or unknown filter is not PASS, and nor is any support but READS.
        expect("NA", "READS", "LOW_QUALITY", "representation_not_evaluated");
        expect("PASS", "NA", "LOW_QUALITY", "no_read_support");
        std::cout << "block quality gates passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "block quality gate FAILED: " << error.what() << '\n';
        return 1;
    }
}
