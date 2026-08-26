#pragma once

#include "panvar/genotype_blocks.hpp"
#include "panvar/genotype_reads.hpp"

#include <cstddef>
#include <limits>
#include <cstdint>
#include <string>
#include <vector>

namespace panvar {

struct GenotypeOptions {
    double recomb_rate = 1.0;          // Li-Stephens switch scaling; 1.0 = one expected switch per locus
    // How much block-local emission a state may give up and still be reachable by the chain. Any
    // diploid state losing more than this to the block's best is excluded BEFORE forward-backward,
    // so the posterior and GQ are computed under the constraint rather than patched after it.
    //
    // Measured motivation: over 6 loci x 20 donors, linkage moved off a UNIQUE local optimum 93
    // times -- 20 of those rescued the call and 73 broke it, and the two separate cleanly by how far
    // the chain moved (rescues median 0.15, overrides median 1.57).
    //
    // Infinity (the default) means unrestricted and reproduces the unconstrained chain exactly. Note
    // 0 is NOT "off": it admits only states tied with the block optimum, which still lets linkage
    // resolve a tie -- the case where it demonstrably earns its place.
    double max_linkage_emission_loss = std::numeric_limits<double>::infinity();
    double error_background = 0.0;     // 0 = estimate from the observed marker counts
    double overdispersion = 0.0;       // negative-binomial phi; 0 = estimate from the depth anchors
    // Adjacent syncmers share reads, so treating markers as independent observations overstates
    // evidence and produces overconfident GQ. Discount the log-likelihood by the effective number of
    // independent observations: roughly the block span divided by the fragment length, not the marker
    // count. 0 disables the correction.
    double fragment_len = 350.0;
    // Measured, not guessed: sweeping the threshold, precision plateaus from about 10 upward, so
    // every point above it costs recall and buys nothing.
    double min_gq = 10.0;
    double min_explained = 0.5;        // below this, the panel does not account for what was observed
    // Two-sided bound on obs/pred. Below it the reads do not account for what the call predicts (a
    // dropout); above its reciprocal the call does not account for what the reads carry (multiplicity
    // under-stated). 0 disables.
    double min_detected = 0.5;
    std::size_t max_alleles_per_block = 64;   // prune by detected-marker fraction before pairing
    // Score the block's TOTAL count as its own term, at its own effective sample size (fragments over
    // the block, not markers), with each pair's overall scale profiled out of the per-marker product so
    // the two do not double count.
    //
    // OFF by default: measured, and it does not help. At its honest weight the term is far too small to
    // overturn a composition preference, and forced hard enough to win it makes identity to the truth
    // worse -- it chases each haplotype's nearest LENGTH independently of sequence, which is not
    // jointly achievable. Kept so the negative result stays reproducible.
    double mass_weight = 0.0;
    // Probability that a marker's count comes from sequence the candidate allele does not model rather
    // than from the allele itself, mixed against a flat component (unmodelled sequence sits at an
    // unknown copy number). Without it a marker the candidate lacks is predicted at the error
    // background, so a real count there costs tens of log units and becomes an unbounded veto -- and
    // under leave-one-out every candidate lacks some of the sample's sequence.
    //
    // OFF by default: the diagnosis is right and the remedy is not enough. It moves the most identical
    // available pair up the emission ranking, monotonically in eps, so the mechanism is real -- but it
    // never flips a call, and pushing eps far enough to beat the residual preference for the longer
    // allele discards real evidence. Kept so the measurement is reproducible.
    double marker_outlier = 0.0;
    // Hard window, in relative deviations, on a block's total marker multiplicity, applied only at
    // tandem arrays that have no bypass allele. Candidates outside it are not scored. 0 disables.
    //
    // Why a constraint is wanted at all: a likelihood maximising read explanation always prefers a
    // SUPERSET allele at an array, since an allele containing the sample's array plus extra units
    // explains every read and the extra copies cost only a modest over-prediction. The emission is not
    // wrong to prefer it -- total length is information the emission does not have.
    //
    // OFF by default, for a decisive reason: it breaks leave-zero-out, refusing haplotypes the panel
    // demonstrably contains. What it does buy is the array's TOTAL LENGTH, to within a fraction of a
    // repeat unit -- so constrain the number you report, in the copy-number output, not the allele you
    // call. A soft penalty instead of a window cannot work: at its honest precision it is worth a
    // couple of log units against composition preferences two orders of magnitude larger.
    double mass_window = 0.0;
    // Huber threshold on the standardised residual, in standard deviations. 0 disables.
    //
    // A Poisson or negative-binomial term has UNBOUNDED influence: a marker the candidate does not
    // carry is predicted at the error background, so real counts there cost tens of log units, while a
    // few percent error on a well-predicted marker costs almost nothing. A handful of absent markers
    // can therefore outvote a thousand carrying nearly all the observed mass. Under leave-one-out every
    // candidate lacks some of the sample's sequence and a longer array lacks less, so that asymmetry
    // has a direction -- toward calling too long.
    //
    // A bounded loss is the standard answer to contaminated data: quadratic while the residual is
    // small, linear beyond. The multinomial does NOT solve this -- o * log(p) with p near zero is just
    // as unbounded.
    double robust_c = 0.0;
    // Score a block by the SHAPE of its count vector rather than its magnitude: model the observed
    // counts as a multinomial draw whose proportions come from the candidate pair, instead of as
    // independent Poissons around absolute predicted means.
    //
    // This is what "allele balance" means formally, and it is how SNV callers decide zygosity, because
    // a ratio survives what absolute depth does not -- a contaminated marker set, a mis-estimated
    // lambda, unusual coverage. It is also the only form in which a marker set a paralogue has inflated
    // remains usable at all.
    bool compositional = false;
    // Weight on a syncmer-ADJACENCY term alongside the per-node one. Adjacencies are computed for
    // every allele and, until this existed, never read by the emission -- a whole evidence class was
    // carried and discarded.
    //
    // It needs its own weight rather than entering at full strength: nodes and adjacencies come from
    // the SAME reads, so summing both at weight 1 double counts, which was measured to be strictly
    // worse at low depth. 0 disables, and that is the default until a stratified sweep says otherwise.
    double edge_weight = 0.0;
    // Weight on the SCALE half of the compositional emission, relative to the shape half. The total's
    // nominal precision is the number of fragments crossing the block, but that assumes counting noise
    // is the only error -- and it is not: lambda carries about 7% uncertainty of its own, and markers
    // a paralogue has inflated break the count scale while leaving the ratio usable. Swept, not assumed.
    double scale_weight = 1.0;
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

// Node coverage as an alternative emission for selected blocks. The marker path stays the default and
// is untouched; this replaces only what the emission is computed FROM, so blocks, pruning, chaining,
// filters and scoring are shared and any difference is attributable to the evidence.
struct CoverageEvidence {
    // [block][allele][node] traversal multiplicity, empty for blocks that keep the marker emission.
    std::vector<std::vector<std::vector<std::uint32_t>>> block_allele_nodes;
    std::vector<double> node;            // the sample's per-node coverage
    double lambda = 0.0;                 // depth per haplotype copy, from invariant nodes
    std::vector<char> use_block;         // per block: 1 = score from coverage
    // Total sequence across both haplotypes that the coverage implies at each block, and its standard
    // error. Used by the combined mode, where the marker emission still CHOOSES the alleles and this
    // only says how much sequence there should be between them.
    std::vector<double> target_bp;
    std::vector<double> target_sd;
};

struct BlockCall {
    // Diagnostic: where the truth's allele pair ranks by emission alone, ignoring linkage. Separates
    // "the emission is wrong" from "linkage overrode a correct emission".
    int truth_allele1 = -1;
    int truth_allele2 = -1;
    int truth_emission_rank = -1;
    double truth_emission_delta = 0.0;   // log-likelihood gap to the best-emission pair
    // How many pairs share the top emission value. Rank counts STRICTLY better pairs, so a block
    // whose markers cannot separate anything gives every pair the same score and reports rank 1 --
    // which reads as "the emission got it right" when the emission said nothing at all. Without this
    // count, an identifiability failure and a linkage override are indistinguishable in the output.
    int truth_emission_ties = -1;
    // How many candidates the block actually SCORED, after pruning to --max-alleles. The tie count
    // is over these, not over every allele in the block, so dividing it by the full allele count
    // understates how tied a large block is -- badly, above the 64-allele default. Emitted so no
    // consumer has to infer min(n_alleles, budget), which is only correct at the default.
    int n_scored_alleles = -1;
    // Where the CALLED pair sits by emission alone. Together with the truth's rank this measures the
    // chain in both directions, which a one-sided count cannot: a delta of 0 means linkage kept the
    // block-local optimum, a large negative delta means it overrode strong local evidence, and a
    // block where the emission's own optimum is wrong while the call is right is linkage EARNING its
    // place. Counting only the overrides would price linkage without its benefit.
    int called_emission_rank = -1;
    double called_emission_delta = 0.0;   // called pair's emission minus the best pair's
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
    double detected = 0.0;   // obs/pred over the called pair's markers; 1.0 is a perfect fit
    // Total sequence this block's reads account for, across both haplotypes, measured WITHOUT going
    // through the panel's allele list: the observed marker mass divided by the depth gives the
    // multiplicity the sample carries, and the called pair's own multiplicity-per-base converts it to
    // bases. It is the answer to "how much sequence is there", which at a tandem array is the question,
    // and it is continuous where an allele call is quantised to whatever lengths the panel happens to
    // hold. `mass_bp_sd` is its standard error, from the number of independent fragments over the span.
    double mass_bp = 0.0;
    double mass_bp_sd = 0.0;
    double called_bp = 0.0;  // the called pair's own total, for comparison
    // Highest multiplicity any of this block's markers reaches in any allele. 1-2 is ordinary
    // variation; a large value means a tandem array, and there the two length figures above answer
    // different questions -- see `is_array`.
    // Total sequence the ALIGNMENT coverage implies at this block, across both haplotypes. Independent
    // of the markers and of which allele pair was called, so it is the arbiter when the two evidence
    // paths disagree about how much sequence is present. 0 where coverage was not computed.
    // How many of this block's alleles carry no informative marker at all. Such an allele is invisible
    // to the emission: it predicts nothing, so a pair containing it is separated from a pair containing
    // its neighbour only by the ABSOLUTE count level of the other allele's markers, never by
    // composition. That makes the call rest entirely on the depth scale being right.
    //
    // Measured on the paralogous fixture, where one allele of a two-allele deletion bubble has zero
    // markers: a true heterozygote reads 1.78 x lambda on the other allele's five markers and a true
    // homozygote reads 2.33 x lambda, against predictions of 1.0 and 2.0. The two classes are half as
    // far apart as the model believes and the heterozygote sits nearer the homozygous prediction, so it
    // is called homozygous at GQ 20 and PASS.
    std::size_t alleles_without_markers = 0;
    // True when the call is homozygous at a block that also holds an allele with no markers. The
    // competing heterozygous hypothesis -- one copy of the called allele and one of the invisible one --
    // predicts exactly HALF the counts on exactly the same markers, so the two differ only by the
    // absolute depth scale and never by which markers are present. Every other genotype decision this
    // model makes rests on composition, which is scale-free; this one rests on lambda being right.
    bool scale_only = false;
    double cov_bp = 0.0;
    std::size_t max_copies = 1;
    // True when the block is a tandem array. It changes how the row should be read, which is why it is
    // a column rather than something the reader has to infer:
    //   allele1/allele2 and called_bp -- the closest panel alleles BY CONTENT. At an array the panel
    //     rarely holds the sample's own arrangement, and the likelihood favours whichever allele
    //     carries most of the sample's repeat-unit variants, which is not the one with the right
    //     number of copies, so at a large array called_bp can miss by several repeat units.
    //   mass_bp +- mass_bp_sd -- how much sequence the READS say is there, which is the copy number,
    //     and which lands within a couple of percent where called_bp does not.
    // Taking called_bp for the copy number at an array is the mistake this column exists to prevent.
    bool is_array = false;
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
    const std::vector<int>* truth_allele2 = nullptr,
    const CoverageEvidence* coverage = nullptr);

// `probe_target` records whether the truth_* columns describe this sample's truth or a pair supplied
// by --probe-pair, so a diagnostic run can never be mistaken for a scored one downstream.
void write_genotypes(
    const std::string& out_prefix,
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    const std::vector<BlockCall>& calls,
    const std::vector<std::string>& haplotype_names,
    bool probe_target = false);

} // namespace panvar
