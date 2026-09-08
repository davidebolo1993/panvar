#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

namespace panvar {

// ---------------------------------------------------------------------------------------------
// THE ONE CHAIN INFERENCE KERNEL.
//
// There were two forward-backward implementations -- production's factorised, scaled-probability
// recurrence and a generic log-space one written for the hybrid fixture -- and adding a hybrid
// variant inside the production caller would have made three. Three implementations of one
// recursion is how a factorised branch and a linked branch end up disagreeing without either being
// obviously wrong. This is the single kernel; the brute-force path oracle stays deliberately
// independent, because an oracle that shares code with the thing it checks proves nothing.
//
// TWO EDGE PATHS, and a fixture must exercise BOTH IN ONE CHAIN or a branch can be present but
// never run:
//
//   linkage-free edge  the existing O(n_h^2) factorised arithmetic, preserved exactly. Each
//                      homologue switches independently, so a step is row sums then column sums.
//                      Bit-identical to the legacy recurrence: that is what makes "hybrid disabled
//                      reproduces the legacy caller" a structural claim rather than a tolerance.
//
//   linkage edge       the generic UNNORMALISED psi transition, O(n_h^4). It is never
//                      row-normalised: forward-backward normalises its messages globally, and
//                      row-normalising would change the model and make one edge's linkage depend on
//                      unrelated outgoing states.
//
// NUMERICAL NOTE. An earlier version of this comment claimed exp(log psi) would OVERFLOW because
// the phase spread reaches 16043 nats on the two-bubble fixture. That was wrong. psi is mean-one
// within each content class, so over a class of size |C| the values average to 1 and therefore
//
//     max log psi <= log|C| <= log 4 = 1.3863
//
// -- measured: the informative self-test case peaks at exactly log 2 = 0.693147. The large
// magnitudes are all NEGATIVE (its minimum is -39.3), losing configurations underflowing toward
// zero, which is the arithmetically correct answer for something e^-39 less likely. exp() of a
// mean-one psi cannot exceed |C|.
//
// The per-edge shift below is therefore NOT rescuing an overflow. It is kept because it costs
// nothing, cancels exactly in the per-block normalisation, and makes this kernel safe if it is ever
// handed an arbitrary unnormalised log-potential rather than a mean-one one -- at which point the
// bound above no longer holds. It is accumulated into the weight sum so the partition stays
// recoverable. The mean-one invariant itself is asserted by --linkage-selftest, not assumed here.
struct ChainEdgeLinkage {
    bool active = false;
    std::size_t n_a = 0, n_b = 0;
    std::vector<std::uint32_t> allele_a;   // per haplotype -> allele index at the LEFT block
    std::vector<std::uint32_t> allele_b;   // per haplotype -> allele index at the RIGHT block
    std::vector<double> log_psi;           // [((a1*n_b + b1)*n_a + a2)*n_b + b2]
};

// ---------------------------------------------------------------------------------------------
// THE GROUPED SPARSE EDGE, the third path. Instead of an n_h^4 loop over a dense table:
//
//     F'(y) = F'_LS(y) + SUM_x F(x) T(x,y) [psi(x,y) - 1]
//
// The first term is the EXISTING factorised O(n_h^2) recurrence, untouched. Only non-neutral phase
// classes contribute a correction, and psi - 1 is formed with expm1 so a psi near one does not lose
// its correction to cancellation.
//
// T(x,y) factorises as (1-r)I + (r/n) 11^T per homologue, so the correction over a source group
// expands into four terms -- identity/identity, identity/switch, switch/identity, switch/switch --
// each computed from group row sums, column sums and totals rather than by revisiting equivalent
// haplotype transitions. A haplotype is indexed by its LEFT and RIGHT allele, so a class visits
// only the compatible state groups.
struct SparsePhaseClass {
    std::uint32_t amin = 0, amax = 0, bmin = 0, bmax = 0;
    double straight_m1 = 0.0;   // expm1(log psi_straight)
    double crossed_m1 = 0.0;    // expm1(log psi_crossed)
};

struct SparseEdgeLinkage {
    bool active = false;
    std::size_t n_a = 0, n_b = 0;
    std::vector<std::uint32_t> allele_a, allele_b;              // per haplotype
    std::vector<std::vector<std::uint32_t>> group_a, group_b;   // allele -> haplotypes
    std::vector<SparsePhaseClass> classes;
    // An edge with no classes is EXACTLY the linkage-free case and must take the factorised path
    // unchanged -- not a correction loop that happens to add zero.
    bool has_corrections() const { return active && !classes.empty(); }
};

// Counters, so a fixture can REQUIRE that both paths ran. A branch that exists but is never taken
// is not covered, however many assertions surround it.
struct ChainKernelStats {
    std::size_t factorised_edges = 0;
    std::size_t linked_edges = 0;
    std::size_t grouped_edges = 0;        // edges taking the sparse grouped-correction path
    std::size_t corrections_applied = 0;  // ordered source-group -> target-group corrections
    // Roundoff can drive a corrected weight slightly negative when psi is near zero and the
    // correction nearly cancels the baseline. A TINY negative is clamped and counted; a materially
    // negative one is an internal error, not something to absorb.
    std::size_t clamped_negatives = 0;
    double worst_negative = 0.0;
    // THE LOG OF AN UNNORMALISED WEIGHT SUM, not a probability-model partition function. The initial
    // ordered states are given weight 1 each, NOT a uniform prior 1/n_hap^2, so this exceeds the
    // partition of a model carrying that prior by exactly 2*log(n_hap). Block posteriors are
    // unaffected -- the difference is a constant that divides out -- and the brute-force oracle uses
    // the same convention, so their agreement is real. It matters only if this value is ever used
    // for absolute fit, calibration, or comparison across panels of different sizes, which is why it
    // is named for what it is rather than called a partition.
    double log_weight_sum = 0.0;
};

// `block_emissions(bi, e)` must fill `e` with n_hap*n_hap unnormalised emission WEIGHTS (not logs)
// for block bi, indexed i * n_hap + j. `edges` is either null/empty (every edge factorised) or of
// size n_blocks, where edges[b] joins block b-1 to b; edges[0] is unused.
// `sparse`, when given, is consulted first: an edge with corrections takes the grouped path, an
// active edge without them takes the factorised path, and anything else falls back to `edges`.
void chain_forward_backward(
    std::size_t n_hap, std::size_t n_blocks, double recomb,
    const std::function<void(std::size_t, std::vector<double>&)>& block_emissions,
    const std::vector<ChainEdgeLinkage>* edges,
    std::vector<std::vector<double>>& fwd,
    std::vector<std::vector<double>>& bwd,
    ChainKernelStats* stats = nullptr,
    const std::vector<SparseEdgeLinkage>* sparse = nullptr);

}  // namespace panvar
