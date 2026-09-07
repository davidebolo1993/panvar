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
// NUMERICAL NOTE, and it is not optional. The production recurrence works in scaled PROBABILITIES,
// while psi is a log ratio whose spread reaches thousands of nats on real fixtures (16043 measured
// on the two-bubble phase fixture). exp() of that overflows outright. So each linkage edge is
// shifted by its own maximum log psi before exponentiating; the shift is a constant per edge, it
// cancels in the per-block normalisation, and it is accumulated into `log_scale` so the partition
// function is still recoverable. Configurations far below the maximum underflow to zero, which is
// the arithmetically correct answer for something e^-16000 times less likely.
struct ChainEdgeLinkage {
    bool active = false;
    std::size_t n_a = 0, n_b = 0;
    std::vector<std::uint32_t> allele_a;   // per haplotype -> allele index at the LEFT block
    std::vector<std::uint32_t> allele_b;   // per haplotype -> allele index at the RIGHT block
    std::vector<double> log_psi;           // [((a1*n_b + b1)*n_a + a2)*n_b + b2]
};

// Counters, so a fixture can REQUIRE that both paths ran. A branch that exists but is never taken
// is not covered, however many assertions surround it.
struct ChainKernelStats {
    std::size_t factorised_edges = 0;
    std::size_t linked_edges = 0;
    double log_scale = 0.0;    // sum of log normalisers plus every linkage shift: log Z
};

// `block_emissions(bi, e)` must fill `e` with n_hap*n_hap unnormalised emission WEIGHTS (not logs)
// for block bi, indexed i * n_hap + j. `edges` is either null/empty (every edge factorised) or of
// size n_blocks, where edges[b] joins block b-1 to b; edges[0] is unused.
void chain_forward_backward(
    std::size_t n_hap, std::size_t n_blocks, double recomb,
    const std::function<void(std::size_t, std::vector<double>&)>& block_emissions,
    const std::vector<ChainEdgeLinkage>* edges,
    std::vector<std::vector<double>>& fwd,
    std::vector<std::vector<double>>& bwd,
    ChainKernelStats* stats = nullptr);

}  // namespace panvar
