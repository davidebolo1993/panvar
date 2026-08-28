#pragma once

// FRAGMENT-LEVEL GENOTYPE PROTOTYPE (arm F4 of experiments/FRAGMENT_EVIDENCE_PREREGISTRATION.md).
//
// The production emission reduces the reads to a vector of marker counts in `count_reads`, which
// keeps the sequence of each read and discards its name, its mate and its offsets. Every mechanism
// added since -- confinement, over-expected, marker clumps, the rho discount, the adjacency channel
// -- reconstructs, indirectly, the linkage that reduction destroyed. This path keeps it instead:
// the unit of observation is the physical fragment, and each fragment contributes exactly one term
// to a candidate's likelihood.
//
// It is a SEPARATE path on purpose. Nothing here is reachable from `panvar genotype`, no production
// default moves, and the two can be run over the same reads and the same panel so that a difference
// is attributable to the observation unit alone.
//
// What it deliberately does NOT do yet: no chain, no HMM, no linkage between blocks, no depth or
// dosage channel. A block is scored on its own fragments. Those are the next layers and they are
// only worth building if this one passes its gate.

#include "panvar/genotype_blocks.hpp"
#include "panvar/gfa.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace panvar {

// One physical fragment: the two mates of a pair, or a single read. Mates are joined by name, so
// this works for interleaved input and for split R1/R2 files without being told which it is.
struct Fragment {
    std::string name;
    std::string r1;
    std::string r2;   // empty for a single-end fragment
    std::size_t bases() const { return r1.size() + r2.size(); }
};

struct FragmentLoadStats {
    std::size_t reads = 0;
    std::size_t fragments = 0;
    std::size_t paired = 0;
    std::size_t singleton = 0;
    // Names seen more than twice. wgsim and samtools fastq both give the two mates one shared name,
    // so a third occurrence means the input is not what this assumes and the pairing is unreliable.
    std::size_t over_paired = 0;
};

// Read every FASTA/FASTQ in `paths`, strip a trailing /1 or /2 and anything after the first space,
// and group by the remaining name.
std::vector<Fragment> load_fragments(
    const std::vector<std::string>& paths,
    FragmentLoadStats* stats = nullptr);

// ONE insert prior, used by the exact reference and by the accelerated scorer alike. They previously
// had different ones -- a normalised discrete distribution in the reference against a continuous
// Gaussian plus a uniform term divided by a fixed span in the fast path -- which is a MODEL
// difference, so any discrepancy between them could not be attributed to acceleration.
struct InsertPrior {
    long lo = 0, hi = 0;
    std::vector<double> logp;          // indexed by L - lo, normalised to sum to 1
    double log_at(long L) const {
        if (L < lo || L > hi || logp.empty()) return -std::numeric_limits<double>::infinity();
        return logp[static_cast<std::size_t>(L - lo)];
    }
    // SUM over L of pi(L) * max(0, |H| - L + 1): the number of (start, L) states a haplotype offers,
    // which is what the contract normalises by.
    double exposure(std::size_t hap_len) const;
};

InsertPrior make_insert_prior(double mean, double sd, double discordant_rate,
                              int sigmas, long min_len);

struct FragmentScoreOptions {
    std::size_t kmer_size = 31;
    std::size_t syncmer_s = 0;          // 0 = default_syncmer_s(k)
    // Bases of neighbouring-block sequence glued to each side of a candidate allele. A fragment
    // overlapping the block boundary only has somewhere to land if the candidate carries its
    // flank; a 150 bp read needs at least a read length, and an insert needs a fragment length.
    // The flank is the SAME sequence for every candidate of a block, so it cancels out of every
    // pairwise comparison -- it buys reachability, it cannot buy preference.
    std::size_t flank_bp = 500;
    std::size_t min_recruit_hits = 2;   // syncmer hits against a block's contexts to recruit
    std::size_t max_alleles = 64;       // candidates kept per block, by coarse containment
    double error_rate = 0.01;           // per-base probability of an edit
    // The divergence at which a fragment stops being able to express a preference. Read alignment
    // likelihood is unbounded below, so one fragment from sequence no candidate models can outvote
    // hundreds of ordinary ones -- and under leave-one-out EVERY candidate lacks some of the
    // sample's sequence, so that asymmetry has a direction. Flooring each fragment's contribution
    // at the likelihood of a `bg_divergence` read is the bounded-loss answer, at fragment level
    // rather than at marker level where `--marker-outlier` measured it and found it too late.
    double bg_divergence = 0.10;
    double outlier_mix = 0.05;
    // Divide each candidate's likelihood by the number of fragment start positions its context
    // offers. Without it P(fragment|allele) is not a generative likelihood and a SUPERSET allele is
    // never penalised: every read from a short allele also fits an allele that contains it plus
    // extra, and the extra costs nothing because nothing says reads should have come from it. That
    // is the same "prefers the longer allele" pathology the marker emission has, and it is why
    // --mass-window exists there. Here it falls out of the model instead of being bolted on.
    bool length_normalize = true;
    // Score each fragment against the REST OF THE LOCUS as well, and let that be its background.
    // A block scored in isolation has to explain every fragment handed to it, so a read from a
    // paralogous copy elsewhere in the locus votes here -- and at CYP2D6 it votes for whichever
    // candidate is most paralog-like. Giving the fragment somewhere else to belong is what stops it,
    // and it needs no filter and no truth input: the rest of the locus is panel sequence.
    bool compete = true;          // prior weight on that background component
    bool use_insert_size = true;
    double fragment_len = 350.0;
    double fragment_sd = 50.0;
    // Probability that a pair is discordant -- wrong orientation, wrong copy of a repeat, chimeric.
    // Without this the insert term is a Gaussian and therefore UNBOUNDED, and the damage is not
    // hypothetical: measured at cyp2d6, NA18939's own haplotype 1 carries a 13.6 kb duplication, its
    // two mates anchored to different copies, the implied insert came out at 13.6 kb, and z^2/2 alone
    // cost 35,000 nats -- so the sample's own haplotype scored 15x worse than an unrelated one. A
    // bounded loss is the standard answer and it is the same one this project already reached for
    // marker counts.
    double discordant_rate = 0.01;
    double discordant_span = 10000.0;   // the uniform the discordant component spreads over
    // Alignment band, as a fraction of read length. Beyond it edlib bails and the read is scored at
    // the background floor for that candidate.
    double max_divergence = 0.20;
    std::size_t threads = 0;
    // Write the per-fragment, per-candidate log-likelihood matrix for this chain index, plus each
    // fragment's best edit distance. A block-level score is a sum over hundreds of fragments and
    // nothing about it can be checked from the sum alone -- whether the contexts are right, whether
    // the reads reach them, whether a preference comes from ten fragments or four hundred. -1 = off.
    long debug_block = -1;
    std::string debug_path;
};

// What one block's fragments say about one candidate pair.
struct PairScore {
    std::size_t allele1 = 0;
    std::size_t allele2 = 0;
    double score = 0.0;   // sum over fragments of log(0.5 P(f|a1) + 0.5 P(f|a2))
};

struct BlockFragmentResult {
    std::size_t block_index = 0;
    BlockKind kind = BlockKind::Bubble;
    std::size_t bubble_id = 0;
    std::size_t n_alleles = 0;          // in the panel, before pruning
    std::size_t n_candidates = 0;       // scored
    std::size_t n_fragments = 0;        // recruited
    // Recruited fragments whose log-likelihood is not flat across the candidates. A block can recruit
    // hundreds of fragments and still be told nothing by any of them; the two numbers separate "no
    // reads" from "no discrimination", which a single count cannot.
    std::size_t n_informative = 0;
    std::size_t best_a = 0;
    std::size_t best_b = 0;
    double best_score = 0.0;
    std::size_t top_class = 0;          // pairs within `tie_eps` of the best
    // Diagnostics, filled only with --truth-haplotypes. -1 = the truth allele is not representable
    // in the reduced panel; -2 = it is, but coarse pruning dropped it before scoring. Those are
    // different failures and must never collapse into one another.
    int truth_a = -1;
    int truth_b = -1;
    int truth_rank = -2;
    int truth_ties = 0;
    double truth_delta = 0.0;
    std::vector<PairScore> top_pairs;   // best `top_pairs_kept`, descending
};

// Score every block in `targets` (chain indices) from the fragments, independently of one another.
std::vector<BlockFragmentResult> genotype_fragments(
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    const std::vector<Fragment>& fragments,
    const std::vector<std::size_t>& targets,
    const FragmentScoreOptions& options,
    const std::vector<int>* truth_allele1 = nullptr,
    const std::vector<int>* truth_allele2 = nullptr,
    std::size_t top_pairs_kept = 20,
    double tie_eps = 1e-6);

// ---------------------------------------------------------------------------------------------
// WHOLE-HAPLOTYPE MODE
//
// Measured, at cyp2d6 block 5, HG00096, leave-zero-out: the block-local model above cannot be right
// there, and the reason is not its scoring. The panel's allele 7 (556 bp) contains the truth label's
// allele 0 (261 bp) and both start at the same position in the sample's own haplotype -- the extra
// 295 bp of the sample's real sequence lives in the NEIGHBOURING block. Any block-local context has
// to glue some flank there, the sample's flank is unknown (that is the genotyping problem), and a
// fragment spanning the junction then prefers whichever candidate supplies more real sequence before
// the guessed flank starts. The evidence that decides the block sits across the block boundary.
//
// So the candidate is a whole locus haplotype, scored end to end, and per-block calls are a
// projection of the answer rather than the unit of inference. No flank is guessed, because each
// haplotype carries its own.
//
// Still absent, deliberately: mosaics. A pair of complete panel haplotypes cannot represent a sample
// that is haplotype A across one block and haplotype B across the next, so under leave-one-out this
// is bounded by the panel's mosaic ceiling and will lose blocks to it. That bound is measurable
// here, which is the point of building this before any factor graph.

struct HaplotypeScoreOptions : FragmentScoreOptions {
    std::size_t max_haplotypes = 48;   // shortlist kept after coarse syncmer containment
    // Occurrences PER HAPLOTYPE above which a syncmer anchors nothing, because inside a tandem
    // array it points everywhere. Per haplotype, not across the shortlist: a global cap makes
    // anchoring depend on how many haplotypes were shortlisted, and it did -- at cyp2d6 NA18939 the
    // same reads gave the truth rank 2 at --max-haplotypes 48 and rank 1 at 96, purely because the
    // larger shortlist pushed more codes past a shared cap.
    // LEFT AT 8. Raising it to 64 was measured and NOT adopted. Inside a tandem array a cap of 8 discards most of the
    // array's own syncmers, so its fragments never anchor and the array contributes nothing.
    //
    // Measured over all 10 lpa donors, not the sweep sample: total excess 793,921 -> 542,595, -32%.
    // The effect is strongly heterogeneous and the two donors in the 8/16/32/64/128 sweep were the
    // two best cases, so the -91% and -81% they showed are NOT representative:
    //
    //   -90% NA18508   -91% HG01106   -81% HG02572   -38% HG01975   -18% HG00146
    //    -2% HG04184     0% HG00735     0% HG02391     0% NA19240   +45% HG03239
    //
    // A 32% aggregate gain with one donor 45% worse is evidence for an occurrence-aware anchoring
    // redesign, not for moving a global default -- the same conclusion, on the same shape of
    // evidence, that this project reached for --edge-weight, where a locus-dependent optimum meant
    // the default stayed put. A single cap cannot be right for a locus with and without arrays at
    // once. The replacement keeps repetitive anchors and weights them by inverse occurrence, seeding
    // each fragment from its rarest syncmers, and caps candidate start bins rather than marker
    // multiplicity. Until that exists the default does not move.
    //
    // Consequence for the joint model's baseline: lpa's target excess is the cap-8 figure, 793,921,
    // and 542,595 is what a different anchoring policy would already recover. Both are recorded so
    // the model is not credited with anchoring's share.
    //
    // lpa remains 542,595 above a panel floor of 7,835. The cap was a contributor, not the cause.
    std::size_t max_anchor_occ = 8;
    std::size_t anchor_slack = 40;     // bases of window either side of an anchored read start
    // Haplotype pairs to report the likelihood's opinion of, by name. Probing costs nothing: the
    // pair scores already exist.
    std::vector<std::pair<std::string, std::string>> probe_pairs;
    // Sum a fragment's likelihood over EVERY placement it has on a haplotype instead of taking its
    // best one:
    //
    //   P(f|h) = (1/N_h) * SUM over valid placements p of P(f|h,p)
    //
    // A fragment compatible with ten positions is evidence for a haplotype that offers ten, and the
    // maximum represents that identically to a haplotype offering one. The measured reason to try it:
    // at cyp2d6 the ceiling pair is shortlisted in 10 of 10 donors with 100% of fragments placing on
    // it, and still loses by a median of 2,323 nats -- so neither candidate generation nor placement
    // availability is the failure, and the geometry of the likelihood is what is left.
    //
    // The 1/N_h is not optional with this on. Summing placements without it rewards a repetitive
    // haplotype for offering more places to land, which is the superset pathology in a new costume --
    // so this implies length normalisation unless it is explicitly overridden.
    bool marginalise_placements = false;
    // Implied-start bins kept per mate per haplotype. 2 is what best-placement scoring used; a sum
    // over placements is only meaningful if the placements are actually enumerated.
    std::size_t placement_topk = 2;
    // Width of the implied-start bucket used to collapse anchors into one placement. 1 disables the
    // binning entirely, which rung zero needs.
    std::size_t placement_bin = 64;
    // RUNG ZERO: the accelerated event enumeration with none of the approximations -- exact contract
    // exposure instead of anchor-dependent windows, no anchor cap, no start binning, no topk. If this
    // does not reproduce the exact reference then the difference is in the model, not in the
    // acceleration, and no ladder of approximations below it means anything.
    // NOT "no approximations": the anchor cap, start binning and top-k are removed, but placements
    // are still reached by syncmer recruitment and banded local alignment rather than by enumerating
    // every (start, L, strand) state. Those two remain until this matches the reference exactly.
    bool rung_zero = false;
    // Emission by Hamming distance at the anchor's implied start, with no alignment. The reference
    // scores a placement at a FIXED position; banded local alignment can slide a read to a better
    // offset and so create placements the reference rejects, which is a model difference and not an
    // acceleration. This exists to remove it from the comparison.
    bool hamming_emission = false;
    // Write each fragment's total placement log-mass for the top pair. The ladder needs retained
    // placement PROBABILITY MASS, not a count of placements: dropping ten negligible placements and
    // dropping one dominant one are the same number and completely different facts.
    std::string dump_fragment_mass;
    // Which pair the mass dump describes. Empty = whatever ranked first, which is NOT comparable with
    // a reference dump for a named pair -- the two files can silently describe different diplotypes.
    std::string dump_mass_pair1, dump_mass_pair2;
    // Merge radius, in bp, for collapsing placements that share a midpoint. 0 keeps every distinct
    // (start, end) state.
    //
    // This is an APPROXIMATION IN ITS OWN RIGHT and was previously hard-coded at 16 bp inside what was
    // called "unrestricted recruitment". Inside a tandem array, mate pairings across non-adjacent
    // copies can share a midpoint and were being collapsed by max, so a mass loss attributed to
    // syncmer recruitment may have been this instead. It gets its own rung.
    std::size_t placement_dedup = 16;
    // MULTIPLICITY-PRESERVING COMPRESSION, replacing the two fixed pruning heuristics.
    //
    // Measured: on a 10-copy array `--max-anchor-occ 8` retains 4% of the placement mass and
    // `--placement-topk 2` retains 7%, and each flips the call. Both discard evidence to bound work.
    // Raising them instead trades correctness for unbounded runtime, which is not a fix.
    //
    // The compression is exact where it matters: inside a perfect repeat every copy yields the SAME
    // likelihood, so N identical placements are one group carrying log P + log(N). An array costs one
    // group rather than N placements, and nothing is thrown away.
    //
    // Pruning is then by RETAINED MASS rather than by count: groups are dropped only while the
    // omitted probability mass stays under `mass_tolerance`. If that cannot be met the run says so
    // instead of silently keeping two placements.
    // STATUS on the 10-copy fixture: the grouping is EXACT. Scored under --rung-zero, so that the
    // exact contract exposure is used on both sides, it gives 0.0022 nats against the reference --
    // identical to the same recruitment without grouping. The cap retains 4% of the placement mass
    // and top-k 7%, and both flip the call; this retains 100% and keeps the correct winner.
    //
    // An earlier 0.19-nat residual was attributed to the grouping and that was wrong: the arm had
    // been run WITHOUT --rung-zero, so it changed the grouping and the exposure model at once, and
    // the difference was entirely the window exposure.
    //
    // Restricted to the global marginal model, where a placement's position does not enter the score.
    // Under the windowed model equal-likelihood placements in different windows are not
    // interchangeable and collapsing them would move coverage between windows.
    bool multiplicity_aware = false;
    double mass_tolerance = 1e-3;
    // How a haplotype-pair posterior becomes a per-block allele pair.
    //
    //   map      take the best pair's alleles. The answer is then a real pair some haplotype pair
    //            actually realises.
    //   marginal take, per block, the allele pair with the most posterior mass summed over pairs.
    //
    // `marginal` sounds strictly better -- pairs that disagree elsewhere may agree here, and that
    // agreement is evidence. Measured at cyp2d6 under leave-one-out it is much worse: the MAP pair
    // reconstructs the donor's own sequence to a median of 130 edits and the marginal projection of
    // the SAME posterior to 385. Taking a per-block argmax of marginals assembles a combination no
    // single pair realises, and a chimera of good pairs is not a good sequence. Verified not to be a
    // spelling artefact: concatenating a haplotype's own block alleles reproduces its sequence at 0
    // edits over 213 kb.
    bool project_map = true;
    // TOTAL-DOSAGE CHANNEL, and deliberately not another coverage weight.
    //
    // Alignment identity cannot see copy number: reads from a sample's single copy align perfectly
    // to a candidate carrying two near-identical copies. Measured, cyp2d6 HG04036 carries a
    // haplotype 12,121 bp shorter than the panel median and the caller misses by 12,136 -- the
    // absent sequence almost exactly.
    //
    // The per-haplotype window channel could not fix it, for two reasons, and only one was a tuning
    // problem. The rate was fitted per haplotype, which normalises absolute depth away; and every
    // fragment was counted on EVERY haplotype it placed on, so one real fragment made two candidate
    // copies both look covered. No shared rate repairs the second -- the reads have already been
    // duplicated.
    //
    // This asks a coarser question that needs neither: does the pair's TOTAL length explain how many
    // fragments were seen at all?
    //
    //     N_fragments ~ Poisson( lambda * (L(h1) + L(h2)) )
    //
    // with lambda estimated ONCE, outside any candidate, from the observed fragment count over twice
    // the panel's median haplotype length. Each fragment is counted once by construction, because
    // the observation is a single scalar.
    bool total_depth = false;
    double haploid_depth = 0.0;      // fragments per bp per haplotype copy; 0 = estimate externally
    // Diagnostic ONLY: score pairs on how close their total length is to the sample's true total,
    // which no caller can know. It bounds what perfect dosage knowledge could buy, so a failure here
    // means the tail is not a dosage problem at all.
    double truth_total_bp = 0.0;
    // JOINT FRAGMENT-ASSIGNMENT + WINDOW-DEPTH MODEL.
    //
    // The hypothesis, pre-registered in docs/reports/genotype-fragment-joint-model.md: localised
    // depth with each fragment assigned EXACTLY ONCE recovers error that neither alignment nor total
    // dosage can. Both cheap answers are already measured out -- a correctly calibrated scalar total
    // is worth 8.6 nats on a 12 kb error against alignment differences of 10^3-10^4, and the
    // per-haplotype window channel counted every fragment on every haplotype it placed on, so one
    // real fragment made every candidate copy look covered.
    //
    // What changes here is the constraint, not the weight. For a candidate PAIR, each fragment is
    // assigned to one homologue or to a null state, window counts are built from those assignments
    // alone, and the objective is a single likelihood:
    //
    //     sum_f log P(fragment | its assignment)  +  sum_w log Poisson(n_w | lambda * window)
    //
    // There is deliberately no weight between the two terms. Coordinate ascent: assign, recount,
    // repeat. A fragment that both homologues explain equally contributes to whichever window needs
    // it, which is the mechanism a copy-number difference is visible through and the one an
    // independent per-haplotype count destroys.
    bool joint_depth = false;
    // Sum over placements instead of assigning each fragment to one of them:
    //
    //   log L(a,b) = -lambda * E(a,b)
    //              + sum_f log[ eta * P_bg(f) + (1-eta) * lambda * sum_p P(f | p) ]
    //
    // where p ranges over the fragment's placements on EITHER homologue and E is the pair's total
    // exposure. This is the observed-data likelihood of a Poisson process: each physical fragment
    // enters once, and there is no freedom to pick a window assignment that flatters the depth
    // profile.
    //
    // In the IDEAL form, placement multiplicity would carry copy number on its own -- a fragment with
    // two placements on a two-copy haplotype contributes twice the density. **This implementation is
    // truncated and does not have that property.** `max_anchor_occ` drops common anchors and
    // `placement_topk` keeps only the top few implied-start clusters per mate, so inside a repeat the
    // placement count saturates well below the true copy number. Multiplicity carries copy number
    // only up to that truncation, and the cost of the truncation is unmeasured.
    //
    // Parameter-free relative to the hard assignment: it removes a degree of freedom rather than
    // adding a regulariser. It exists because the hard-assignment arm gained about +575 nats where a
    // real length difference exists and lost about -437 nats where none does, and assignment freedom
    // is the suspect -- untested until this arm runs.
    bool joint_marginal = false;
    std::size_t joint_top_pairs = 32;   // candidate pairs re-scored jointly, by alignment rank
    // Cap, not a schedule: the optimiser stops early when no fragment moves. Whether it actually
    // reached that state is reported rather than assumed -- a stable fragment-count invariant says
    // nothing about having reached an optimum.
    std::size_t joint_iterations = 50;
    // Reverse the order fragments are visited in. Coordinate ascent is order-dependent, so a result
    // that changes under this is optimiser instability rather than evidence about the model.
    bool joint_reverse_order = false;
    // Score tolerance, in nats, defining the equivalence set. Pairs within this of the best are
    // reported as not distinguished by the evidence rather than ranked against each other. 0 means
    // exact ties only.
    double equivalence_tolerance = 0.0;
    std::size_t equivalence_max_report = 16;
    std::size_t joint_window = 500;
    // Sequence compatibility and copy number are different signals and a read alignment cannot carry
    // both. Measured, at cyp2d6 leave-ZERO-out: NA18939's haplotype 1 is 13.6 kb longer than the
    // panel's typical haplotype -- a duplication -- and the reads from the extra copy align perfectly
    // to the single copy in a shorter haplotype, so alignment alone sees nothing and the shorter one
    // wins. What separates them is that the longer haplotype's extra 13.6 kb would be UNCOVERED if
    // the sample did not carry it. That is a depth observation, scored here in non-overlapping
    // windows of a haplotype's own length, with the rate fitted per haplotype so a homologue at half
    // depth is not penalised for being one of two.
    //
    // It reuses the same fragments as the sequence term, so the two are not independent and summing
    // them at full weight double counts -- the defect --edge-weight was measured to have. The weight
    // is exposed rather than fixed for that reason.
    // OFF by default, on the same reasoning this project applied to --edge-weight. The channel is
    // real and the mechanism is right, but it reuses the fragments the sequence term already used, so
    // at full weight the two double count -- and measured at cyp2d6 leave-ZERO-out it is what broke
    // NA18939, moving its own true pair from rank 1 to rank 2. A term that fails when the answer is
    // provably in the panel does not get to be a default. Kept switchable so the sweep is possible.
    double coverage_weight = 0.0;
    std::size_t coverage_window = 500;
    // A global 1/L factor is the WRONG way to price length here and is off by default in this mode:
    // with the depth rate profiled out it is unidentifiable, and any fixed proxy for it penalises
    // exactly the true haplotype in the duplication case above. Left switchable so that is
    // reproducible.
};

struct HaplotypePairScore {
    std::size_t hap1 = 0;
    std::size_t hap2 = 0;
    double score = 0.0;
    double posterior = 0.0;
};

// Per block, what the haplotype-pair posterior implies about this block's allele pair. The call is
// marginal, not a read-out of the single best pair: two haplotype pairs disagreeing everywhere else
// may agree here, and that agreement is evidence the best-pair read-out throws away.
struct BlockProjection {
    std::size_t block_index = 0;
    BlockKind kind = BlockKind::Bubble;
    std::size_t bubble_id = 0;
    std::size_t n_alleles = 0;
    // PHASED: allele1 is always the first haplotype of the called pair and allele2 the second,
    // across every block. Sorting them per block by index -- which both this and production's table
    // used to do -- is invisible to an unordered allele-pair comparison and fatal to a sequence one:
    // concatenating allele1 down the chain then switches homologue at every block and produces a
    // chimera of the two. Measured at cyp2d6: the same called pair reconstructs the donor at 130
    // edits phased and 385 unphased.
    int allele1 = -1;
    int allele2 = -1;
    double posterior = 0.0;
    // True when every member of the equivalence set carries this same allele pair here. A block can
    // be determined even where the haplotype pair is not, and that is the part of the answer worth
    // reporting when allocation is ambiguous.
    bool determined = false;
    int truth_a = -1;
    int truth_b = -1;
    bool exact = false;
    bool truth_representable = false;
};

// Per shortlisted haplotype, what each channel says about it on its own. A pair score is a sum of
// tens of thousands of terms and nothing about it is checkable from the sum: whether a haplotype
// lost on sequence or on depth, whether its reads placed at all, whether the coarse stage even gave
// it a chance. These are the columns that answer that.
// How much of the placement evidence the recruiter actually kept. Reported so that a read-length
// ladder cannot mistake a RECRUITMENT limit for an INFORMATION limit: if completeness is low, a
// result saying "350 bp cannot resolve this" is a statement about max_anchor_occ and placement_topk,
// not about the reads. Nothing may be called an information limit while these are far from complete.
struct PlacementCompleteness {
    // POST-ANCHOR retention. These clusters are what survives max_anchor_occ, so a high figure here
    // says placement_topk was not binding -- it does NOT say the recruiter saw everything. Naming it
    // "completeness" invites exactly that misreading, which is why the two stages are separate.
    std::uint64_t clusters_found = 0;               // implied-start clusters offered, per mate, post-anchor
    std::uint64_t clusters_kept = 0;                // ...surviving placement_topk
    std::uint64_t fragments_truncated = 0;          // fragments where placement_topk actually bound
    std::uint64_t fragments_total = 0;
    // The stage BEFORE that, which the above cannot see.
    std::uint64_t anchor_occurrences_seen = 0;
    std::uint64_t anchor_occurrences_dropped = 0;   // excluded by max_anchor_occ
    // Ground truth, available only when the reads were simulated from haplotypes that are IN the
    // panel: wgsim encodes each fragment's origin in its name, so the placement the recruiter should
    // have found is known. Raising placement_topk restores true placements AND adds false ones, so a
    // retention figure alone cannot say whether the evidence improved.
    bool recall_measured = false;
    std::uint64_t truth_resolvable = 0;   // fragments whose origin haplotype is in the shortlist
    std::uint64_t truth_recovered = 0;    // ...with a retained placement at the true position
    // Retained placements elsewhere on the ORIGIN haplotype. Not necessarily errors: inside a repeat
    // the other copies are legitimate alternatives, and they are exactly what a copy-number model
    // needs. Counted separately from recall because raising placement_topk trades one against the
    // other, and a single "completeness" number hides which is moving.
    std::uint64_t spurious_placements = 0;
    // What --mass-tolerance omitted from the placements RECRUITMENT FOUND. It is not placement
    // completeness: states recruitment never generated are not in the denominator, and measured on
    // stochastic reads that omission is far larger. Naming this "omitted mass" without the
    // qualification invites reading "bound met" as "nothing was lost".
    double omitted_mass_max = 0.0;
    double omitted_mass_mean = 0.0;
    // COST, so that "correct" and "tractable" are separate claims. The compression is DOWNSTREAM:
    // anchors are still expanded, every recruited placement is still aligned, and mate combinations
    // are still formed -- only then are equal likelihoods grouped. So it restores the evidence the
    // cap and top-k destroyed, but it may retain the work they were introduced to avoid. An array
    // costs one group rather than N placements only AFTER enumeration.
    std::uint64_t anchor_hits = 0;                  // index lookups that yielded a position
    std::uint64_t mate_combinations = 0;            // (mate1, mate2) pairs actually scored
    std::uint64_t placements_before_grouping = 0;
    std::uint64_t groups_after_grouping = 0;
};

struct HaplotypeScore {
    std::string name;
    std::size_t bp = 0;
    std::size_t placed = 0;         // fragments that landed on it
    std::size_t zero_windows = 0;   // depth windows with no fragment
    std::size_t windows = 0;
    double coverage_ll = 0.0;
    double solo_ll = 0.0;           // sum over fragments of log P(fragment | this haplotype)
    double containment = 0.0;       // the coarse shortlist score, reported not trusted
};

// What the likelihood thinks of a NAMED haplotype pair. Production has `--probe-pair` for the same
// reason and the reasoning carries over: a score carries a baseline that shifts whenever the reads or
// the shortlist do, so scores from different runs are not on a common scale, while two pairs probed
// in the SAME run are. It also separates three failures a truncated top-N list cannot: the pair was
// never shortlisted, the pair was shortlisted and scored badly, or the pair was scored and is close.
struct HaplotypeProbe {
    std::string name1, name2;
    bool in_shortlist = false;
    int rank = -2;              // -2 = not shortlisted, so never scored
    double score = 0.0;
    double delta = 0.0;         // score minus the best pair's
    std::size_t placed1 = 0;    // fragments that landed on each, for telling a placement failure
    std::size_t placed2 = 0;    // apart from a likelihood failure
};

// Whether the coordinate ascent actually converged, per rescored pair. Reported because "it ran five
// iterations" and "it reached a fixed point" are different claims and only one of them licenses
// reading the result as the model's opinion.
struct JointConvergence {
    std::size_t pairs_rescored = 0;
    std::size_t pairs_converged = 0;
    std::size_t max_iterations_used = 0;
    std::size_t total_moves_last_iteration = 0;
};

// What the evidence actually distinguishes, as opposed to which pair happened to rank first.
//
// Measured reason this is not a nicety: at cyp2d6 HG04036 the sequence-correct pair sits 100-262 nats
// behind the winner out of a score near 200,000, and the two differ by 12,161 edits of reconstruction.
// Reporting one pair there asserts a distinction the data does not support. An equivalence set states
// the ambiguity instead of resolving it arbitrarily, and it is the honest short-read output at an
// array whether or not any depth model is ever adopted.
struct EquivalenceSet {
    double margin = 0.0;            // best score minus the runner-up: how much the call is actually won by
    std::size_t size = 0;           // pairs within `tolerance` of the best
    double posterior_mass = 0.0;    // their combined posterior
    // Alleles the set AGREES on, per block: where every member carries the same pair, the call is
    // determined even though the haplotype pair is not. This is the part of the answer that survives
    // the ambiguity, and at a locus where allocation is unidentifiable it may be most of it.
    std::size_t blocks_determined = 0;
    std::size_t blocks_total = 0;
    std::vector<std::string> members;   // names, capped for reporting
};

struct HaplotypeResult {
    std::vector<std::string> shortlist;         // haplotype names actually scored
    std::vector<HaplotypeScore> haplotypes;
    std::vector<HaplotypePairScore> top_pairs;
    std::size_t n_fragments = 0;
    std::size_t n_informative = 0;
    std::vector<BlockProjection> blocks;
    std::vector<HaplotypeProbe> probes;
    JointConvergence convergence;
    EquivalenceSet equivalence;
    PlacementCompleteness completeness;
};

// Whole-haplotype mode rests on an assumption it never checked:
//
//     concatenated block alleles for path X  ==  the raw GFA spelling of path X
//
// When it fails the caller scores a sequence the panel does not contain, silently. Measured on a
// tandem array built as a chain of identical nodes: one haplotype spelled 1100 bp against its true
// 1500, another spelled 0 bp, and a full table of plausible numbers came out anyway. A test-only
// check is not enough -- any graph whose decomposition does not round-trip produces this.
//
// Throws naming the path, both lengths and the first mismatching offset.
void verify_block_spelling(
    const Graph& graph,
    const std::vector<BlockAlleles>& blocks,
    const std::vector<std::string>& haplotype_names);

HaplotypeResult genotype_haplotype_pairs(
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    const std::vector<std::string>& haplotype_names,
    const std::vector<Fragment>& fragments,
    const HaplotypeScoreOptions& options,
    const std::vector<int>* truth_allele1 = nullptr,
    const std::vector<int>* truth_allele2 = nullptr,
    std::size_t top_pairs_kept = 20);

void write_haplotype_results(
    const std::string& out_prefix,
    const HaplotypeResult& result,
    bool have_truth);

// Spell a per-block call table into the two sequences it claims the sample carries.
//
// Needed because the two callers report in different shapes and cannot otherwise be compared on the
// objective that matters. Production emits a per-block allele pair -- a mosaic, not a haplotype pair
// -- so it has no "called haplotype" to align. Concatenating its own called alleles along the chain
// gives exactly the sequence its output asserts, and that can be aligned to the donor's truth like
// anything else. `allele1`/`allele2` of -1 (a block with no call) contribute nothing, which is the
// same convention a bypass allele already has.
void spell_called_pair(
    const std::vector<BlockAlleles>& blocks,
    const std::vector<int>& allele1,
    const std::vector<int>& allele2,
    std::string& seq1,
    std::string& seq2);

// ---------------------------------------------------------------------------------------------
// FLOORS
//
// Three different questions about how well the PANEL could possibly do, which the single
// "best complete pair" number conflates:
//
//   complete   one panel haplotype per homologue across the whole locus. What a non-mosaic caller
//              can reach at best.
//   free       the nearest panel allele at every block, chosen independently. A mathematical
//              optimum: it may switch source haplotype at every boundary with no evidence for any
//              of the switches, so it is a bound and not a target.
//   penalised  the same, with a cost for changing source haplotype between blocks. Between the two,
//              and the only one of the three whose switches are constrained at all.
//
// H_mosaic = complete - free is what a mosaic model could recover AT MOST. Small, and a factor graph
// is not worth building; large, and the next question is whether the switches have spanning-fragment
// evidence, which is a separate measurement again.
struct MosaicFloors {
    std::size_t blocks = 0;
    std::size_t blocks_scored = 0;      // where the truth traverses and the panel offers an allele
    std::size_t complete = 0;
    std::size_t free_mosaic = 0;
    std::size_t switches_free = 0;      // source-haplotype changes the free optimum uses
    std::vector<std::pair<double, std::size_t>> penalised;   // (penalty per switch, total)
    std::vector<std::size_t> penalised_switches;
};

// Per homologue, so the two can differ and be reported separately: one may be represented well and
// the other not, and a summed figure hides that entirely.
MosaicFloors mosaic_floors(
    const std::vector<BlockAlleles>& blocks,
    const std::vector<std::string>& haplotype_names,
    const std::vector<std::string>& truth_block_seq,
    const std::vector<double>& switch_penalties,
    std::size_t threads);

// ---------------------------------------------------------------------------------------------
// EXACT REFERENCE SCORER
//
// The model of docs/reports/genotype-fragment-model-contract.md, implemented as directly and as
// slowly as possible: every fragment start on both haplotypes is enumerated, the insert prior is
// summed over, and there is no syncmer index, no anchor cap and no placement_topk anywhere in it.
//
// It exists to be an ORACLE, not a caller. It is O(fragments x haplotype length x insert range) and
// is meant for small synthetic haplotypes where the answer can also be worked out by hand. The fast
// path is accepted only if it reproduces this ranking and these likelihood differences within a
// stated bound -- and the fast path is never adjusted to make them agree, because this defines the
// model and the fast path only approximates it.
struct ReferenceParams {
    // kept in step with FragmentScoreOptions so both scorers build the same prior
    double lambda = 0.05;        // fragments per start position
    double eta = 0.05;           // background weight
    double error_rate = 0.01;
    double fragment_len = 350.0;
    double fragment_sd = 50.0;
    double bg_divergence = 0.10; // the background's implied per-base disagreement
    int insert_sigmas = 4;       // how far into the insert prior's tails to sum
    // Same concordant/discordant mixture the accelerated path uses. A pure Gaussian here would be a
    // MODEL difference, so any discrepancy between the two scorers could not be blamed on
    // acceleration -- which is the only thing the differential test is for.
    double discordant_rate = 0.01;
};

// log L(a,b) under the contract. `hap_a` and `hap_b` are the two homologues; pass the same sequence
// twice for a homozygous pair and the exposure and placement set both double, as they must.
double reference_pair_loglik(
    const std::string& hap_a,
    const std::string& hap_b,
    const std::vector<Fragment>& fragments,
    const ReferenceParams& params,
    // Optional: per-fragment log SUM_p P(f|p) over both homologues, in fragment order. This is the
    // exact placement mass the accelerated recruiter is trying to retain.
    std::vector<double>* fragment_mass = nullptr,
    // The full per-fragment likelihood term, log[(1-eta)*lambda*e^m + eta*P_bg]. Needed so the
    // decomposition reconciles: at a fixed pair the exposure cancels, so per-fragment deltas must sum
    // exactly to the whole-pair difference. Placement mass alone cannot do that -- a fragment with no
    // placement has mass -inf while its real contribution is the finite background term.
    std::vector<double>* fragment_contrib = nullptr);

void write_fragment_results(
    const std::string& out_prefix,
    const std::vector<BlockFragmentResult>& results,
    bool have_truth);

} // namespace panvar
