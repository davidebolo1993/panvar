#pragma once

#include "panvar/genotype_blocks.hpp"
#include "panvar/genotype_reads.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace panvar {

struct GenotypeOptions {
    double recomb_rate = 1.0;          // Li-Stephens switch scaling; 1.0 = one expected switch per locus
    double error_background = 0.0;     // 0 = estimate from the observed marker counts
    double overdispersion = 0.0;       // negative-binomial phi; 0 = estimate from the depth anchors
    // Adjacent syncmers share reads, so treating markers as independent observations overstates
    // evidence and produces overconfident GQ. Discount the log-likelihood by the effective number of
    // independent observations: roughly the block span divided by the fragment length, not the marker
    // count. 0 disables the correction.
    double fragment_len = 350.0;
    double min_gq = 20.0;
    double min_explained = 0.5;        // below this, the panel does not account for what was observed
    double min_detected = 0.5;         // below this, the reads do not account for what the call predicts
    std::size_t max_alleles_per_block = 64;   // prune by detected-marker fraction before pairing
    // Down-weight markers by how many of the block's alleles carry them. At a block with hundreds of
    // alleles the marker set is swamped by markers shared across many of them: one carried by 200 of
    // 450 discriminates almost nothing, yet thousands of such terms outvote the few allele-specific
    // ones. Weight = (n_alleles / carriers)^beta, renormalized to mean 1 so the likelihood scale (and
    // with it the ESS discount and GQ) stays comparable. 0 disables.
    double carrier_weight = 0.0;
    // Attribute each call to the blocks that determined it. Costs one extra forward-backward pass per
    // block plus a cache of every block's emissions (n_blocks * n_haplotypes^2 doubles), so it is a
    // diagnostic rather than something to leave on for a cohort.
    bool provenance = false;
    std::size_t threads = 0;
};

struct BlockCall {
    // Diagnostic: where the truth's allele pair ranks by emission alone, ignoring linkage. Separates
    // "the emission is wrong" from "linkage overrode a correct emission".
    int truth_allele1 = -1;
    int truth_allele2 = -1;
    int truth_emission_rank = -1;
    double truth_emission_delta = 0.0;   // log-likelihood gap to the best-emission pair
    std::size_t block_index = 0;
    bool is_bubble = true;
    std::size_t bubble_id = 0;
    std::size_t allele1 = 0;
    std::size_t allele2 = 0;
    double gq = 0.0;
    std::size_t n_markers = 0;   // distinct panel markers surviving in this block
    double explained = 0.0;      // share of observed marker mass the called pair accounts for
    // The complement, and the one that catches a coverage dropout. `explained` asks how much of what
    // we SAW the call accounts for, so when almost nothing was seen every call explains all of it and
    // the metric reads 1.0. `detected` asks how much of what the call PREDICTS actually turned up. A
    // block whose reads are missing predicts counts that are not there, and the likelihood then picks
    // whichever allele predicts least -- collapsing a copy-number call to its minimum, confidently,
    // because every alternative fits even worse.
    double detected = 0.0;
    std::size_t hap1 = 0;        // most probable panel haplotype pair at this block
    std::size_t hap2 = 0;
    double hap_posterior = 0.0;
    std::string filter = "PASS";     // quality only: PASS / LOWGQ / OFFPANEL / NOCALL
    // Where the call came from, kept separate from quality. `local` = this block's own markers decided
    // it; `linked` = it has none and the chain carried it. With --provenance this refines to
    // self / neighbours / distant / none.
    std::string evidence = "local";
    // Which blocks this call actually depended on, found by neutralizing each block in turn and
    // seeing whose call changes. "self" means the block's own markers decide it; "neighbours" means
    // the adjacent blocks carry it; "distant" means the evidence came from elsewhere in the chain.
    // Empty and "none" unless GenotypeOptions::provenance is set.
    std::vector<std::size_t> influencers;
    std::string provenance = "none";
};

struct GenotypeSummary {
    std::size_t blocks = 0;
    std::size_t called = 0;
    std::size_t no_calls = 0;
    std::size_t off_panel = 0;
    double mean_gq = 0.0;
    double lambda_hap = 0.0;
    double overdispersion = 0.0;
    double error_background = 0.0;
};

// Diploid Li-Stephens over the block chain: hidden state is a pair of panel haplotypes, emission is
// the negative-binomial likelihood of the observed marker counts given the allele pair those two
// haplotypes induce, transitions allow each haplotype to switch independently. Forward-backward, so
// every block gets a posterior rather than only a best path.
std::vector<BlockCall> genotype_sample(
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    const ReadPanel& panel,
    const ReadCounts& counts,
    const std::vector<BlockDepth>& depth,
    const std::vector<std::string>& haplotype_names,
    const GenotypeOptions& options,
    GenotypeSummary* summary = nullptr,
    const std::vector<int>* truth_allele1 = nullptr,
    const std::vector<int>* truth_allele2 = nullptr);

void write_genotypes(
    const std::string& out_prefix,
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    const std::vector<BlockCall>& calls,
    const std::vector<std::string>& haplotype_names);

} // namespace panvar
