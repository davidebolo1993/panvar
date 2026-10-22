// I3-A.1: the batch resource contract and the between-sample release, deterministically (no real memory pressure).
#include "panvar/fragment_likelihood.hpp"
#include "panvar/prepared_locus.hpp"

#include <iostream>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error(message);
}

using panvar::BatchResourceLedger;
using panvar::BatchResourceVerdict;

void test_ledger() {
    const std::optional<std::uint64_t> na;   // an unavailable reading: never zero
    {   // compliant: every boundary verified and at or below the ceiling, the ceiling itself included
        BatchResourceLedger l(100, "test_peak");
        require(l.record("start", 1, 1) && l.record("prepare", 40, 60) && l.record("sample", 70, 100) &&
                    l.record("cleanup", 20, 100) && l.record("teardown", 10, 100),
                "a batch whose lifetime peak reaches but never exceeds the ceiling may continue");
        require(l.verdict() == BatchResourceVerdict::compliant && l.exceeded_during().empty() && l.max_peak() == 100u &&
                    l.unverified_boundaries() == 0,
                "compliant at exactly the ceiling");
        require(l.phases()[1].peak_increase == 59u && l.phases()[2].peak_increase == 40u && l.phases()[3].peak_increase == 0u,
                "the lifetime peak's rise is attributed to the phase during which it rose");
    }
    {   // exceeded during preparation: no sample may start; later boundaries are still recorded; no recovery
        BatchResourceLedger l(100, "test_peak");
        require(l.record("start", 1, 1), "start");
        require(!l.record("prepare", 50, 101), "a lifetime peak above the ceiling stops the batch");
        require(!l.record("teardown", 5, 101), "nothing recovers the verdict");
        require(l.verdict() == BatchResourceVerdict::ceiling_exceeded && l.exceeded_during() == "prepare",
                "exceeded during preparation (the phase between the last boundary at or below and the first above)");
    }
    {   // exceeded during the second sample: the first sample's result stands, the third is not started
        BatchResourceLedger l(100, "test_peak");
        require(l.record("start", 1, 1) && l.record("sample", 80, 95) && l.record("cleanup", 30, 95), "first sample");
        require(!l.record("sample", 60, 130), "the second sample's excursion stops the batch");
        require(l.exceeded_during() == "sample" && l.phases().size() == 4 && l.phases()[3].peak_increase == 35u,
                "attributed to the second sample, by the peak's rise across it");
    }
    {   // exceeded at teardown: every sample completed, the batch is still not compliant
        BatchResourceLedger l(100, "test_peak");
        require(l.record("sample", 90, 99) && !l.record("teardown", 10, 105), "a teardown excursion is detected too");
        require(l.verdict() == BatchResourceVerdict::ceiling_exceeded && l.exceeded_during() == "teardown", "teardown");
    }
    {   // a current reading above the ceiling counts even if the peak reading lags
        BatchResourceLedger l(100, "test_peak");
        require(!l.record("sample", 120, 99), "a current footprint above the ceiling is an excursion");
    }
    {   // I3-B: available -> unavailable -> available: the gap forbids COMPLIANT; the batch is not stopped by it
        BatchResourceLedger l(100, "test_peak");
        require(l.record("start", 1, 1) && l.record("prepare", 40, 60), "verified boundaries");
        require(l.record("sample", na, na), "an unverified boundary does not stop the batch");
        require(l.record("cleanup", 30, 90) && l.record("teardown", 10, 90), "verified again");
        require(l.verdict() == BatchResourceVerdict::unverified && l.unverified_boundaries() == 1,
                "one unverified boundary makes the batch UNVERIFIED, never COMPLIANT, though every reading taken was under");
        require(!l.phases()[2].peak && !l.phases()[2].current && !l.phases()[2].peak_increase && !l.phases()[2].verified,
                "the missing readings stay missing (not zero)");
        require(l.phases()[3].peak_increase == 30u && l.max_peak() == 90u,
                "the rise across the gap is attributed to the next boundary actually read; the maximum counts only reads");
    }
    {   // I3-B: one reading missing at a boundary (the other present) is still unverified
        BatchResourceLedger l(100, "test_peak");
        require(l.record("sample", 50, na), "a boundary with only its current reading");
        require(l.verdict() == BatchResourceVerdict::unverified, "a boundary missing either reading is unverified");
    }
    {   // I3-B: unavailable -> observed breach: permanent CEILING_EXCEEDED despite the earlier unverified boundary
        BatchResourceLedger l(100, "test_peak");
        require(l.record("start", na, na) && l.record("prepare", na, na), "unverified boundaries");
        require(!l.record("sample", 80, 140), "an observed breach stops the batch");
        require(!l.record("teardown", na, na), "and stays stopped, readings or not");
        require(l.verdict() == BatchResourceVerdict::ceiling_exceeded && l.exceeded_during() == "sample" &&
                    l.unverified_boundaries() == 3,
                "an observed breach is CEILING_EXCEEDED whatever other boundaries were unverified");
        require(l.phases()[2].peak_increase == 140u, "with nothing read before, the whole peak is this boundary's rise");
    }
    {   // nothing measured anywhere, and an empty ledger: never compliant
        BatchResourceLedger l(100, "unavailable");
        require(l.record("start", na, na) && l.record("sample", na, na), "an unmeasured batch is not stopped");
        require(l.verdict() == BatchResourceVerdict::unverified && !l.max_peak(), "but it is UNVERIFIED, with no peak read");
        require(std::string(panvar::batch_resource_verdict_name(l.verdict())) == "RESOURCE_UNVERIFIED", "name");
        BatchResourceLedger empty(100, "test_peak");
        require(empty.verdict() == BatchResourceVerdict::unverified, "a ledger with no boundary is not compliant");
    }
    std::cout << "  ok: batch resource ledger (compliant at the ceiling; excursions in preparation, a sample and teardown; "
                 "no recovery; attribution by peak rise; per-boundary availability: a gap forbids COMPLIANT, missing is "
                 "never zero, a later breach is permanent CEILING_EXCEEDED)\n";
}

std::string random_bases(std::mt19937_64& rng, std::size_t n) {
    std::string s(n, 'A');
    for (char& c : s) c = "ACGT"[rng() % 4];
    return s;
}

std::string revcomp(const std::string& s) {
    std::string out(s.rbegin(), s.rend());
    for (char& c : out) c = c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : 'A';
    return out;
}

void test_release() {
    std::mt19937_64 rng(31);
    const std::string hap = random_bases(rng, 3000);
    std::vector<panvar::ReadPair> reads;
    for (std::size_t start = 0; start + 300 <= hap.size(); start += 37) {
        panvar::ReadPair p;
        p.first = {"r" + std::to_string(start) + "/1", hap.substr(start, 100), std::string(100, 'I')};
        p.second = {"r" + std::to_string(start) + "/2", revcomp(hap.substr(start + 200, 100)), std::string(100, 'I')};
        reads.push_back(p);
    }
    panvar::LibraryModel lib;
    lib.insert_mean_bp = 300.0;
    lib.insert_sd_bp = 20.0;
    lib.support_sigmas = 4.0;
    lib.substitution_error = 0.002;
    lib.outlier_weight = 0.001;
    lib.fragment_start_intensity = 0.05;
    lib.provenance = "fixture: declared";
    panvar::FragmentScoringOptions o;
    o.mismatch_band = 6;
    o.minimum_piece_bp = 10;
    o.threads = 1;
    o.emission_cache_bytes = 1u << 24;
    {
        panvar::PreparedScoringContext ctx(reads, lib, o);
        (void)ctx.add_haplotypes({hap});
        bool refused = false;
        try {
            (void)panvar::release_sample_retained_storage();
        } catch (const std::logic_error&) {
            refused = true;
        }
        require(refused, "release must refuse while a scoring context pins the emission cache");
    }
    const panvar::RetainedProcessStorage held = panvar::retained_process_storage();
    require(!held.emission_cache_armed && held.emission_cache_table_bytes > 0 && held.calling_thread_scratch_bytes > 0,
            "after the context ends, its cache's capacity and this thread's scratch are still held (the measured residue)");
    const panvar::RetainedProcessStorage released = panvar::release_sample_retained_storage();
    require(released.total() == held.total(), "the release reports what was held");
    const panvar::RetainedProcessStorage after = panvar::retained_process_storage();
    require(after.total() == 0, "nothing is left in the process-global containers after the release");
    const panvar::BetweenSampleCleanup again = panvar::release_between_samples();
    require(again.released.total() == 0 && again.seconds >= 0.0, "a second release finds nothing");
    // the cache still works after a release
    panvar::PreparedScoringContext ctx(reads, lib, o);
    (void)ctx.add_haplotypes({hap});
    std::cout << "  ok: between-sample release refused while pinned; releases " << held.total()
              << " bytes of retained capacity to zero; the cache works afterwards\n";
}

}  // namespace

int main() {
    try {
        test_ledger();
        test_release();
        std::cout << "batch resource: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "batch resource: FAIL: " << error.what() << '\n';
        return 1;
    }
}
