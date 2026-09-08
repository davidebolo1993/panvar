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

#include "panvar/chain_kernel.hpp"
#include "panvar/candidate_frame.hpp"
#include "panvar/genotype_blocks.hpp"
#include "panvar/gfa.hpp"

#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <tuple>
#include <vector>

namespace panvar {

// One physical fragment: the two mates of a pair, or a single read. Mates are joined by name, so
// this works for interleaved input and for split R1/R2 files without being told which it is.
// THE ONE FRAGMENT-NAMING RULE. Strips a trailing /1 or /2 and anything after whitespace, so both
// mates map to one fragment. Exported because marker exclusion identifies reads by fragment, and a
// second copy of this rule would exclude the wrong reads the moment the two drifted apart.
std::string fragment_name(const char* raw, std::size_t len);

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
    // A fragment that cannot be placed within the band is charged at bg_divergence, which ASSERTS a
    // divergence the data never showed -- all that is known is that it exceeds max_divergence. The
    // gap between those two is a cliff at the band edge, and a haplotype carrying more copies of a
    // tandem unit collects the fragments that fall off it. Measured at LPA: 267 net fragments, 0.9%
    // of the data, carrying 84% of a 99,148-edit error at ~227 nats each.
    //
    // With this on, an unplaceable fragment is charged at the BAND BOUNDARY instead: the least
    // penalty consistent with having failed to place. That is what the block-local candidate stage
    // has always done (it adds read_ll(band, len) for a mate that fails to align); the final scorer
    // dropped it, and the two have disagreed ever since.
    bool band_floor = false;
    // THREE fragment-haplotype states, not two. "Unplaceable" currently conflates:
    //
    //   (1) a valid placement                      -> the normal likelihood;
    //   (2) seeds existed, alignment was tried and failed -> genuine negative evidence, charged at
    //       the band boundary (--band-floor) rather than at an asserted bg_divergence;
    //   (3) NO seed was ever offered on this haplotype -> the search did not happen. That is missing
    //       evidence manufactured by the recruiter, and charging it penalises a haplotype for the
    //       approximation rather than for the data.
    //
    // With this on, a fragment unseeded on ANY shortlisted haplotype is dropped from the likelihood
    // entirely. Dropping globally rather than per pair is deliberate: pair totals must sum over the
    // SAME fragment set or they are not comparable, and a per-pair set would silently reintroduce
    // the asymmetry this exists to remove.
    bool unplaced_neutral = false;
    // Haplotypes forced into the shortlist regardless of their containment rank. This is how a
    // known-good pair gets SCORED rather than merely reported as absent: --probe-haplotypes can only
    // describe a pair the shortlist already holds, so it cannot separate "candidate generation lost
    // it" from "the likelihood ranks it below the wrong answer". Forcing removes the first
    // explanation and leaves the second exposed.
    std::vector<std::string> force_haplotypes;
    // Where to write the fragment -> factor incidence table. Step 3 of the block-evidence plan:
    // before any factor is scored, it must be possible to say which fragments inform which factor,
    // and to check that each fragment informs EXACTLY ONE. "A fragment enters the likelihood exactly
    // once" is already the standing invariant; this makes it inspectable rather than assumed.
    std::string incidence_path;
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
    // ONE-MATE RESCUE. When only one mate anchors, the other is currently searched nowhere and the
    // fragment falls to the single-mate path, losing the paired evidence the exhaustive reference
    // finds. Measured: that stratum -- one mate seeded, that mate placed, no valid FR pair -- carries
    // 87-96% of the whole recruitment deficit.
    //
    // The rescue is truth-independent by construction: given an anchored mate's placement on a
    // candidate haplotype, FR orientation and the insert prior's support fix an interval on that same
    // haplotype where the other mate must lie, and it is searched there and nowhere else.
    bool mate_rescue = false;
    // Join mate placements by coordinate instead of by Cartesian product. Exact, not an
    // approximation: the insert prior has strictly bounded support [lo, hi] and log_at returns -inf
    // outside it, so a combination whose implied insert falls outside contributes exactly zero mass
    // and enumerating it is wasted work, not omitted evidence.
    //
    // ON BY DEFAULT, with adaptive dispatch deciding per fragment-haplotype. Justified by
    // measurement rather than by the asymptotics: on a 2-copy fixture where dispatch sends all
    // 20550 decisions to the product -- so the maps and the K count are pure overhead -- end-to-end
    // time is unchanged (0.053s either way), and on a 32-copy array it is 0.293s -> 0.159s. Turn it
    // off with --no-coordinate-join, which is also how the Cartesian diagnostic arm is selected.
    bool coordinate_join = true;
    // DIAGNOSTIC: take the join whenever it is available, bypassing the cost dispatch. Only for
    // tests that must prove the join path itself ran -- with adaptive dispatch a small fixture can
    // send every decision to the product, and an equality assertion then compares the product with
    // itself and passes while testing nothing.
    bool force_join = false;
    // NAMING, because the flag is misleading and the distinction matters: this is a
    // NO-PRIMARY-PLACEMENT fallback, not strictly a zero-seed one. The guard tests
    // a1.empty() && a2.empty() AFTER placement, so it fires when recruitment produced no surviving
    // placement -- which includes fragments that seeded fine but whose anchors all failed to align
    // within the band, not only fragments with no seed. The flag name is kept because the tests and
    // the recorded measurements use it.
    //
    // For a candidate haplotype on which neither mate obtained a primary placement, find placements
    // directly. Exhaustive under the fixed-position Hamming
    // contract -- it considers every start and keeps those within the band, which is exactly what
    // the primary path would keep had recruitment proposed them. Truth-independent, and it
    // introduces no new seed length to tune.
    bool zero_seed_fallback = false;
    // DIAGNOSTIC: disable the pigeonhole filter so the fallback scans every start. The two must
    // return identical placements at the same band -- that is what makes the filter an acceleration
    // rather than a heuristic -- and without this flag the only way to compare them is to change
    // max_divergence, which changes the band and therefore compares two different models.
    bool zero_seed_exhaustive = false;
    // DIAGNOSTIC, not a production proposal: run the fallback for EVERY fragment-haplotype pair and
    // union its placements with the primary ones. The shipped guard only fires when recruitment
    // found nothing at all, so a fragment that recruited one placement but missed its other
    // repeat-copy placements keeps an incomplete list. This measures whether that incompleteness is
    // what the growing non-truth residual is made of.
    bool zero_seed_complete = false;
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
    // Distinct unordered allele pairs this block is given by EVERY locus pair within
    // equivalence_tolerance -- not the truncated reported-member list, which stops at
    // equivalence_max_report and would understate the ambiguity by the members it never printed.
    // -1 means NA: some contributing member could not be projected HERE, so a count would be a
    // claim about alleles that were never established. This draws the distinction the block table
    // needs: a caller can select the right allele and still not have DETERMINED it, when equally
    // supported locus pairs disagree at that block.
    long block_equivalence_size = -1;
    double posterior = 0.0;
    // True when every member of the equivalence set carries this same allele pair here. A block can
    // be determined even where the haplotype pair is not, and that is the part of the answer worth
    // reporting when allocation is ambiguous.
    bool determined = false;
    // The selected pair has NO allele at this block -- the block decomposition does not cover that
    // path there. Distinct from "determined": every top pair agreeing on -1 is agreement about
    // nothing, and reporting it as a determined call would emit a block genotype the decomposition
    // cannot support. Measured at cyp2d6, where one path's chain stops 1978 bp short of its walk.
    // No trustworthy block call here. TWO causes, both fatal to a block-level product:
    //   * the selected pair has no allele at this block (-1); or
    //   * the decomposition does not reproduce one of the selected haplotypes at all -- its block
    //     concatenation differs from its graph walk, so NO per-block allele of that haplotype can be
    //     relied on, even where an allele index exists.
    // The second is the common case and the narrower -1 test misses it: measured on a fixture whose
    // block spelling is 1300 bp against a 1400 bp walk, every block still carried an allele index.
    bool unprojectable = false;
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
    // Separated so that "correct" and "tractable" cannot be argued from one number. The Cartesian
    // count is what the product WOULD cost and grows as (copies)^2 per fragment; the join count is
    // what the coordinate join actually performs, which is O(|F| + |R| + K) for K coordinate pairs
    // inside the insert support. K is NOT bounded by the library's insert width: it is quadratic
    // again whenever coordinates cluster inside one allowed interval, which is what a short tandem
    // array fitting inside a single insert looks like. Measured linear in copy number on the array
    // fixtures, which is a measurement and not a bound. Rescue work is counted apart from both,
    // because the interval scan is linear in the prior's width per anchored placement and would
    // otherwise hide inside whichever total it was added to. The three counters also have DIFFERENT
    // UNITS -- a rescue position is a Hamming comparison with early exit, a join probe is a
    // log-sum-exp -- so their ratio is not a runtime ratio.
    std::uint64_t cartesian_combinations = 0;       // |b1| x |b2|, formed or not
    std::uint64_t join_operations = 0;              // (forward start, insert length) probes
    std::uint64_t rescue_positions = 0;             // interval positions examined by mate rescue
    std::uint64_t rescue_placements = 0;            // of those, ones that passed the band
    std::uint64_t join_chosen = 0;                  // fragment-haplotype decisions dispatched to
    std::uint64_t cartesian_chosen = 0;             // each path by the adaptive rule
    // Zero-seed fallback, reported apart from every other stage so that what it recovers and what it
    // costs are both attributable to it alone.
    std::uint64_t zs_invocations = 0;               // (fragment, haplotype) pairs it ran for
    std::uint64_t zs_candidate_starts = 0;          // starts proposed, before deduplication
    std::uint64_t zs_verified_starts = 0;           // distinct starts actually Hamming-verified
    std::uint64_t zs_placements = 0;                // of those, within the band
    std::uint64_t zs_pigeonhole = 0;                // invocations accelerated by the q-gram filter
    std::uint64_t zs_exhaustive = 0;                // invocations that scanned every start
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
    // Fragments dropped by --unplaced-neutral: unseeded on at least one shortlisted haplotype, so
    // the search never happened there and the fragment cannot discriminate honestly.
    std::size_t n_neutralised = 0;
    std::vector<BlockProjection> blocks;
    std::vector<HaplotypeProbe> probes;
    JointConvergence convergence;
    EquivalenceSet equivalence;
    PlacementCompleteness completeness;
};

// Compare each path's BLOCK concatenation against its GFA walk, and report the disagreements.
//
// This no longer gates scoring. Whole-haplotype mode scores the WALK -- the haplotype is the walk,
// and the decomposition is only used to project an answer onto blocks -- so a decomposition gap
// costs a projection, not a wrong scored sequence. It therefore WARNS on a gap and counts them,
// rather than throwing.
//
// It still THROWS when the graph cannot spell a path at all, because then there is no authoritative
// sequence to fall back on.
//
// History: while the block concatenation WAS the scored sequence this had to be fatal. Measured on
// a tandem array built as a chain of identical nodes -- one haplotype spelled 1100 bp against its
// true 1500, another spelled 0 bp, and a full table of plausible numbers came out anyway.
void verify_block_spelling(
    const Graph& graph,
    const std::vector<BlockAlleles>& blocks,
    const std::vector<std::string>& haplotype_names);

// Which FRAME a path's block spelling is in, relative to its own GFA walk.
//
// The block chain is oriented along the reference, so a path that runs ANTIPARALLEL to the reference
// spells reference-oriented blocks: its concatenated alleles are the exact reverse complement of its
// walk. Both are the same sequence; they differ only in frame. Measured: 60 of 127 paths at cyp2d6
// and at c4, 19 of 465 at ankrd36c -- and the round-trip invariant refused all three loci outright,
// which is what kept the two pharmacogene loci out of every experiment.
enum class SpellFrame { Forward, ReverseComplement, Incomparable };

// Compare a path's block spelling to its walk spelling and report the frame. Incomparable means the
// graph cannot fully spell the walk, or the two differ by something other than orientation -- which
// IS a decomposition fault and must still be refused.
SpellFrame block_spelling_frame(
    const Graph& graph,
    const std::vector<BlockAlleles>& blocks,
    const std::string& name);

// This path's allele in every block of the chain, concatenated, with a block the path bypasses
// contributing nothing.
//
// NOT the scored sequence. Whole-haplotype mode scores the graph WALK; this is the DECOMPOSITION's
// reconstruction of the same path, which is not guaranteed to reproduce it -- measured at cyp2d6,
// where NA18989#1 concatenates to 205236 bp against a walk of 207214, an exact prefix. Kept because
// the two disagreeing is a finding worth reporting, and because block alleles are what the answer is
// projected onto.
std::string spell_block_haplotype(
    const std::vector<BlockAlleles>& blocks,
    const std::string& name);

// A fingerprint of the ALLELE CATALOGUE that integer allele indices are indices into.
//
// An allele index means nothing on its own. It is meaningful only for one exact combination of
// graph, bubble decomposition, and excluded panel paths -- and the catalogue is the product of all
// three, so hashing it captures every input at once without having to enumerate them.
//
// This exists because spelling a call table against the wrong catalogue is silent. Measured: the
// LPA pilot's completion arm reported 56042 edits when the true value was 21, purely because the
// spelling run omitted the --exclude-haplotypes the scoring run had used. The numbers were wrong
// and plausible, which is the same failure mode as the representation drift one layer above.
std::string allele_catalogue_fingerprint(const std::vector<BlockAlleles>& blocks);

// WHOLE-HAPLOTYPE MODE IS A DIAGNOSTIC ORACLE, NOT THE PRODUCT.
//
// Frozen 2026-08-31. It picks one complete panel pair by short-read likelihood, and the LPA
// investigation established that this objective diverges from sequence reconstruction off-panel:
// across 29 development donors the panel represents every one within 940 edits while the caller
// lands a median ~11000 edits above the floor, and a mosaic oracle showed that added
// representational flexibility recovers only ~134 of that. The divergence is not a bug to tune out
// -- argmax P(reads|pair) and argmin d(pair,truth) are simply different pairs at a tandem array.
//
// It remains valuable and is kept: it generates certified floors, the injection experiments, the
// equivalence machinery and an internal baseline. What it must NOT become is the module's output
// contract. That is FragmentBlockCall above. Do not add scoring coefficients here in the hope of closing the
// gap; that was tried and is recorded in the ledger (coverage weighting, refuted across 16 donors).
HaplotypeResult genotype_haplotype_pairs(
    const std::vector<Block>& chain,
    const std::vector<BlockAlleles>& blocks,
    const std::vector<std::string>& haplotype_names,
    const std::vector<Fragment>& fragments,
    const HaplotypeScoreOptions& options,
    const std::vector<int>* truth_allele1 = nullptr,
    const std::vector<int>* truth_allele2 = nullptr,
    std::size_t top_pairs_kept = 20,
    // AUTHORITATIVE SEQUENCE. The graph's own walk for each path, name -> sequence, already
    // canonicalised into the block frame.
    //
    // The haplotype IS the walk; the graph is only a convention for storing it, and spelling a P
    // line end to end recovers it exactly (this is what `odgi paths -f` does). Concatenating
    // per-block alleles is a DECOMPOSITION artifact and is not guaranteed to reproduce it: measured
    // at cyp2d6, NA18989#1 spells 205236 bp against a walk of 207214, the block spelling being an
    // exact PREFIX -- the chain simply stops 1978 bp early and the tail is dropped silently.
    //
    // So scoring uses the walk when it is supplied, and the block decomposition is kept only for
    // PROJECTING the answer back onto blocks. A decomposition gap then costs a projection, never a
    // wrong scored sequence.
    const std::unordered_map<std::string, std::string>* walk_sequences = nullptr);

void write_haplotype_results(
    const std::string& out_prefix,
    const HaplotypeResult& result,
    bool have_truth,
    const std::string& catalogue_fingerprint = std::string());

// ---------------------------------------------------------------------------------------------
// THE BLOCK CALL CONTRACT
//
// What `genotype` reports per block once the fragment evidence engine replaces the marker emission.
// Written as a schema first, deliberately, because the shape of the output is the architectural
// decision: the module reports what the reads determine about each block, and says so when they
// determine only part of it. It does NOT primarily report two named whole-locus haplotypes -- that
// output is retained as a diagnostic oracle.
//
// The LPA investigation is the reason for every field here:
//
//   * dosage separate from composition -- at a KIV-2 array the reads routinely determine the diploid
//     TOTAL while leaving the split between homologues unidentifiable. Reporting a confident pair
//     there was the central failure: the caller returned 29+11 against a truth of 18+18 with the
//     total nearly right.
//   * allele_set, not one pair -- 457 alleles at LPA block 13, many observationally identical. The
//     honest unit is the equivalence class; refining inside it needs linkage that often is not there.
//   * local_gq and linkage_gq SEPARATE -- they fail independently, and a single margin conflates
//     them. Measured: margin is anti-correlated with correctness across the development cohort.
//   * status -- absolute fit, not score separation. A confidently wrong call reported margin 8181
//     while a correct one reported 3807, so a wide margin is not evidence of a good call.
//
// Every field is always emitted. Unknown is "." and never an absent column: a reader must be able
// to tell "not determined" from "not reported", which is exactly the distinction a vacuous audit
// column destroys.
enum class FragmentBlockCallStatus {
    Called,        // an allele pair the fragments support
    Equivalent,    // a set of pairs the fragments cannot separate
    DosageOnly,    // total copy number determined, composition or allocation not
    OffPanel,      // no allele pair explains the fragments adequately
};

// NAMED FragmentBlockCall, not BlockCall. genotype.hpp declares a DIFFERENT BlockCall in the same
// namespace -- the marker caller's -- and the two headers could not be included together, which
// blocked the genotype command from using anything in this file. Two unrelated types sharing a name
// in one namespace was the defect; renaming the fragment-side one is the fix.
struct FragmentBlockCall {
    std::size_t block_index = 0;
    BlockKind kind = BlockKind::Bubble;
    long bubble_id = -1;
    FragmentBlockCallStatus status = FragmentBlockCallStatus::OffPanel;

    // Composition. allele1/allele2 are set only when status == Called; otherwise allele_set carries
    // the members the evidence admits, and -1 means "not determined".
    int allele1 = -1;
    int allele2 = -1;
    std::vector<std::pair<int, int>> allele_set;

    // Dosage, which is often determined when composition is not.
    int diploid_dosage = -1;          // total copies across both homologues, -1 = not determined
    double dosage_confidence = -1.0;

    // Allocation between homologues, reported ONLY where linkage supports it.
    bool allocation_determined = false;

    // Two confidences, because they fail independently: local_gq is "is this the right content",
    // linkage_gq is "is this the right assignment to homologues".
    double local_gq = -1.0;
    double linkage_gq = -1.0;

    // Absolute fit, the basis for OffPanel. Not a margin.
    double fit = -1.0;                // e.g. posterior-predictive p-value; -1 = not computed
    std::size_t fragments_informing = 0;
};

// Write the contract. One row per block, every column always present, "." for not determined.
void write_block_calls(const std::string& path, const std::vector<FragmentBlockCall>& calls);

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

// ---------------------------------------------------------------------------------------------
// EXACT LOCAL REFERENCE SCORERS  (plan step B, prerequisite)
//
// The same contract as reference_pair_loglik, restricted to ONE block and to ONE ADJACENT PAIR of
// blocks. These exist so that a fragment-derived block emission can be checked against an oracle
// BEFORE it is wired into the chain -- the whole-haplotype reference cannot do that, because it
// scores a complete haplotype and says nothing about what any single block is worth.
//
// They are deliberately separate functions because the incidence table distinguishes the evidence:
//
//   local fragment     -> reference_block_loglik            (a unary block emission)
//   boundary fragment  -> reference_block_pair_loglik       (a transition factor)
//
// Putting boundary evidence into a unary emission would double count it or assign it arbitrarily,
// so the oracle keeps the two apart from the start rather than discovering the distinction later.
//
// CONTEXT. A block's sequence alone is not scoreable: a fragment overlapping its edge needs the
// neighbouring sequence to align against. Both take `flank_bp` and build
// left_flank + allele(s) + right_flank from the block chain's majority alleles, which is exactly
// what the block-local candidate stage already does.
//
// Cost is that of the whole-haplotype reference on a context-sized sequence, so these are oracles
// for small blocks and small fragment sets, not callers.

// EXPLICIT STATE, EXPLICIT EVIDENCE. Both take the two homologue sequences already constructed by
// the caller and the fragment subset the caller selected. Nothing is guessed:
//
//   * no majority-allele flanks. A majority flank is a GUESSED haplotype context, and guessing the
//     flank is the defect already recorded at CYP2D6 block 5 -- an oracle that guesses is not an
//     oracle. The caller states the context it means.
//   * no implicit fragment set. The unary oracle must be given only the fragments the incidence
//     table classified `local` to that block; the transition oracle only those classified
//     `boundary` for that adjacent pair. Handing either the whole read set makes the unary score
//     include boundary evidence and the transition score include everything, so their difference
//     measures nothing.
//
// With those supplied, both are exactly reference_pair_loglik over the stated sequences and the
// stated fragments -- which is the point: the local scorer must BE the reference model restricted to
// a factor, not a second model that resembles it.
double reference_factor_loglik(
    const std::string& hap_a,           // homologue 1's sequence for this factor's span
    const std::string& hap_b,           // homologue 2's
    const std::vector<Fragment>& fragments,   // the incidence-selected subset for this factor
    const ReferenceParams& params);

// ---------------------------------------------------------------------------------------------
// FRAGMENT -> FACTOR INCIDENCE, computed ONCE
//
// Recruitment offers a fragment to every target whose alleles or flanks share enough syncmers with
// it. The SHAPE of that target set is what decides which factor the fragment belongs to. Three
// places need this answer -- the incidence table, the reference oracle, and eventually the block
// emission -- and a second implementation of the rule classifies fragments differently: an
// approximation that counted alleles rather than syncmer hits selected ZERO boundary fragments where
// the real rule selects 177, which would have made a phase gate pass while measuring nothing.
enum class FragmentFactor {
    Unrecruited,   // no target: outside the modelled evidence, counted explicitly
    Local,         // exactly one target -> a unary block emission
    Boundary,      // two ADJACENT targets -> a transition factor
    Path,          // more than two consecutive targets -> a higher-order factor
    Ambiguous,     // non-consecutive targets: shared evidence across repeated blocks
};

struct FragmentFactorIncidence {
    std::vector<std::vector<std::uint32_t>> targets_of;  // per fragment, sorted target ranks
    std::vector<FragmentFactor> factor_of;               // per fragment
    std::vector<std::vector<std::uint32_t>> recruited;   // per target, the fragment indices
};

// `targets` are BLOCK indices, in chain order; a fragment's targets_of holds RANKS into that list,
// so adjacency means adjacent TARGETS (consecutive bubble targets are raw blocks 1 and 3).
FragmentFactorIncidence classify_fragment_factors(
    const std::vector<BlockAlleles>& blocks,
    const std::vector<std::size_t>& targets,
    const std::vector<Fragment>& fragments,
    std::size_t kmer_size,
    std::size_t syncmer_s,
    std::size_t flank_bp,
    std::size_t min_recruit_hits);

// The chain's context sequence either side of a block, from majority alleles. Exposed so a caller
// classifying fragments can build the SAME recruitment index the block-local stage builds -- an
// approximation of it would classify fragments differently from the incidence table.
std::string chain_left_flank(const std::vector<BlockAlleles>& blocks, std::size_t bi,
                             std::size_t want);
std::string chain_right_flank(const std::vector<BlockAlleles>& blocks, std::size_t bi,
                              std::size_t want);

// Build one homologue's sequence over a span of the TARGET CHAIN, from explicit alleles.
//
// `targets` is the same ordered target list the incidence table uses, so "adjacent" means adjacent
// TARGETS -- consecutive bubble blocks are targets 1 and 3 by raw index, and any backbone sequence
// between them is included here. An oracle indexed by raw block index would disagree with the
// incidence table about what a boundary is.
std::string chain_span_sequence(
    const std::vector<BlockAlleles>& blocks,
    const std::vector<std::size_t>& targets,
    std::size_t first_target,           // inclusive
    std::size_t last_target,            // inclusive; == first_target for a unary factor
    const std::vector<int>& alleles,    // one per target in [first_target, last_target]
    std::size_t flank_bp);              // context taken from the chain OUTSIDE the span

// One candidate, ready for the oracle: the AUTHORITATIVE walk bytes plus a verified map from that
// walk's coordinates to block indices.
//
// The oracle must not rebuild candidates by concatenating block alleles. The walk is the haplotype;
// the concatenation is the decomposition's reconstruction of it and is not guaranteed to reproduce
// it -- reverse-complemented for antiparallel paths, and short for unprojectable ones. A locally
// rebuilt candidate would reintroduce exactly the defect the walk rule was established to remove.
//
// `ok` is false when the two cannot be reconciled, i.e. there is no trustworthy coordinate map. Such
// a candidate must not silently contribute a scope: its block boundaries are unknown.
// THE PRODUCTION BAND, in one place. floor(divergence * len) + 1, and the +1 matters: at 5% a
// 120 bp read gives 7 here and 6 if the expression is re-derived without it, so a search certified
// against the re-derived value silently loses exactly the boundary placements. It was file-local,
// and the first --bounded-search gate duplicated the formula and got it wrong.
inline std::size_t mate_band_edits(double max_divergence, std::size_t len) {
    if (len == 0) return 0;
    return static_cast<std::size_t>(max_divergence * static_cast<double>(len)) + 1;
}

// ---------------------------------------------------------------------------------------------
// BOUNDED-COMPLETE SINGLE-MATE PLACEMENT (stage 1, Hamming).
//
// Every start where `read` matches `hap` within `max_edits` MISMATCHES -- fixed position, no gaps,
// the same emission reference_emission implements. Pigeonhole: a placement with at most d mismatches
// cannot mismatch inside all of d+1 DISJOINT pieces, so at least one piece matches exactly and its
// occurrences propose the start. Every proposal is then VERIFIED by counting mismatches, so the
// pigeonhole only ever restricts which starts are examined -- it cannot admit a placement the
// exhaustive scan would reject, and cannot reject one it would accept.
//
// Occurrences are NOT capped and there is no top-k. A cap makes completeness a function of a tuning
// parameter, which is the property this search exists to remove.
struct MatePlacement {
    long start = 0;             // offset into hap
    std::uint32_t edits = 0;    // Hamming distance at that offset
    bool operator<(const MatePlacement& o) const { return start < o.start; }
    bool operator==(const MatePlacement& o) const { return start == o.start; }
};

// Work done, so a search that skips can be caught saying so.
struct SearchWork {
    std::uint64_t pieces = 0;            // d+1 per read
    std::uint64_t candidate_starts = 0;  // proposals before dedup
    std::uint64_t distinct_starts = 0;   // after dedup
    std::uint64_t verified = 0;          // Hamming actually computed
    std::uint64_t accepted = 0;          // within band
    bool exhaustive_fallback = false;    // piece too short to filter; scanned every start
};

// REUSABLE UNCAPPED OCCURRENCE INDEX over one haplotype at one piece length.
//
// Without it, bounded_mate_placements locates each piece with hap.find(piece, pos) in a loop, which
// RESCANS the whole haplotype per piece: at 9 pieces x 2 mates x 2 strands that is ~36 full scans of
// a 226 kb sequence per fragment-candidate cell. Measured: ~72 ms per cell, which is 20 hours for
// 23953 fragments x 131 candidates and is why the full-panel run did not finish in 12 hours. The
// index is built ONCE per (haplotype, piece length) and reused across every fragment.
//
// EVERY occurrence is retained -- no cap, no top-k. A capped index would destroy exactly the
// copy-number information the multiplicity gates exist to protect.
struct PieceIndex {
    std::size_t piece = 0;
    std::unordered_map<std::uint64_t, std::vector<std::uint32_t>> at;
    bool usable() const { return piece > 0 && !at.empty(); }
};
PieceIndex build_piece_index(const std::string& hap, std::size_t piece);

// Same contract as the un-indexed overload; `idx` must have been built over `hap` at the piece
// length this call would use, or it is ignored.
std::vector<MatePlacement> bounded_mate_placements(
    const std::string& read, const std::string& hap, std::size_t max_edits, SearchWork* work,
    const PieceIndex* idx);

std::vector<MatePlacement> bounded_mate_placements(
    const std::string& read, const std::string& hap, std::size_t max_edits, SearchWork* work);

// The reference this is certified against: every start, Hamming counted, nothing skipped.
std::vector<MatePlacement> exhaustive_mate_placements(
    const std::string& read, const std::string& hap, std::size_t max_edits);

// THE COORDINATE HALF of the valid-FR rule. Named for what it actually checks: the previous name,
// valid_fr_state, claimed three conditions while taking no strand arguments at all.
//   * checked HERE: the reverse mate DOWNSTREAM -- its end at or after the forward mate's start --
//     and the implied insert inside the prior's support [lo, hi];
//   * enforced by the CALLER, through which vectors it is given: opposite strands. A caller passing
//     two forward vectors gets same-strand pairs and this function cannot tell.
// So the synthetic case table tests COORDINATES; orientation is covered separately by asserting
// that both library orientations actually occur among the enumerated states.
// The DOWNSTREAM half on its own, because the accelerated Cartesian path needs exactly this and not
// the support test: it applies the support through insert_ll returning -inf, which carries the same
// mass but keeps the `mates_seeded` valid-FR bit deliberately BROADER than the prior's support.
// Substituting the full predicate there would silently narrow that bit, so the shared thing is the
// rule both actually agree on.
inline bool fr_reverse_downstream(long fwd_start, long rev_end) { return rev_end >= fwd_start; }

inline bool valid_fr_coordinates(long fwd_start, long rev_end, long insert_lo, long insert_hi) {
    const long insert = rev_end - fwd_start + 1;
    return fr_reverse_downstream(fwd_start, rev_end) && insert >= insert_lo && insert <= insert_hi;
}

// The allowed reverse-END interval for a forward start, shared so the coordinate join's window
// arithmetic and this predicate cannot drift apart.
inline std::pair<long, long> fr_reverse_end_window(long fwd_start, long insert_lo, long insert_hi) {
    return {fwd_start + insert_lo - 1, fwd_start + insert_hi - 1};
}

// One fragment state. The KEY is every field: two placements differing in any of them are distinct
// states, and collapsing any component changes the likelihood. Notably `hap` is part of it -- the
// same coordinates on two haplotypes are two states, not one.
struct FragmentState {
    std::uint32_t hap = 0;
    long m1_start = 0; bool m1_fwd = true;
    long m2_start = 0; bool m2_fwd = false;
    long frag_start = 0, frag_end = 0;
    long insert = 0;
    std::uint32_t m1_edits = 0, m2_edits = 0;
    bool operator<(const FragmentState& o) const {
        return std::tie(hap, m1_start, m1_fwd, m2_start, m2_fwd, frag_start, frag_end, insert) <
               std::tie(o.hap, o.m1_start, o.m1_fwd, o.m2_start, o.m2_fwd, o.frag_start, o.frag_end,
                        o.insert);
    }
    bool operator==(const FragmentState& o) const {
        return !(*this < o) && !(o < *this);
    }
};

// Every valid-FR state for one fragment on one haplotype, from mate placements found by either
// search. Both library orientations are formed: (m1 forward, m2 reverse) and (m2 forward, m1
// reverse). Same-strand, reverse-upstream and out-of-support combinations are rejected by
// valid_fr_state, never by a second rule written here.
std::vector<FragmentState> enumerate_fragment_states(
    std::uint32_t hap,
    const std::vector<MatePlacement>& m1_fwd, const std::vector<MatePlacement>& m1_rev,
    const std::vector<MatePlacement>& m2_fwd, const std::vector<MatePlacement>& m2_rev,
    std::size_t m1_len, std::size_t m2_len, long insert_lo, long insert_hi);

// AN UPPER BOUND ON THE MASS THIS SEARCH DOES NOT REACH.
//
// An omitted state is one where some mate exceeds its edit band. Its per-state likelihood is at
// most the best such state: one mate at exactly d+1 edits and the other perfect, whichever pairing
// is larger, times the largest insert prior in support. The COUNT of omitted states is bounded by
// the total (start, insert, orientation) combinations minus those found -- so the bound is
//
//     log M_omitted <= log(0.5) + B + log SUM_L N_omitted(L) * pi(L)
//
// SUMMED PER INSERT LENGTH, not N_omitted * max_L pi(L). The latter charges the largest insert
// probability to every omitted state, and pi varies by orders of magnitude across the support, so
// it inflates the bound by roughly that ratio. Measured on the tail-only fixture: the flat form
// gave -39.6 against a true reference of -133.0, loose by 93 nats.
//
// N_omitted(L) is taken from the SAME universe exposure integrates over,
//     N_total(L) = 2 * max(0, |H| - L + 1),
// minus the distinct found states at that L, so the mass and exposure invariant is preserved.
//
// B is the best a state OUTSIDE the band can achieve: ONE mate perfect and the other exactly one
// edit past its band -- not both mates out, which is strictly worse and would understate nothing
// but is not the maximum.
//
// It counts the AGGREGATE of omitted origins, not one representative: a repeat with many omitted
// copies contributes all of them through N_omitted. That is the property that makes it safe at an
// array, where bounding a single origin would understate the tail by the copy number.
// STATES AT EXACTLY ONE MISMATCH LEVEL, and their summed mass.
//
// For the multiplicity property, "reference minus in-band" is the WRONG quantity: extra repeat
// copies also add unrelated starts, junction states, insert lengths and orientations, all with
// finite emission, so that difference measures everything the extra sequence brought. Subtraction
// in log space is fragile near cancellation besides, and shows nothing about WHICH states supplied
// the mass. So the class is isolated directly: exactly the states whose two mates carry (e1, e2)
// mismatches, counted and summed.
struct EditClass {
    std::size_t count = 0;
    double mass = 0.0;
};
EditClass edit_class_mass(const std::vector<FragmentState>& states,
                          std::uint32_t e1, std::uint32_t e2,
                          std::size_t m1_len, std::size_t m2_len,
                          const InsertPrior& ip, double log_eps, double log_1meps);

// ADAPTIVE TAIL. The bound above charges B -- one mate at d+1, the other perfect -- to EVERY
// omitted state, when most omitted states are far worse (a random position mismatches ~3/4 of the
// read). That is where the looseness lives: per-insert-length summation removed only ~0.8 nats of a
// 93-nat slack, because pi sums to 1 and the flat form over-counted by roughly a factor of two.
//
// So: enumerate EXACTLY through a deeper band D >= d, and bound only what lies beyond D, using
// D+1 mismatches. Each increase in D moves B down by about one log_eps per edit while the exact
// part absorbs the states that used to be bounded. D grows until the interval width meets the
// declared tolerance or the cap is reached, at which point the fragment is reported UNCERTIFIABLE
// rather than certified on a width nobody checked.
struct TailInterval {
    double lower = 0.0;        // exact mass through D
    double upper = 0.0;        // logadd(lower, bound beyond D)
    double bound = 0.0;        // the beyond-D bound alone
    std::size_t depth = 0;     // the D actually used
    std::size_t states = 0;    // exact states through D
    bool within_tolerance = false;
    // REFINEMENT INVARIANTS, recorded per fragment so deepening cannot go wrong silently.
    // The last is the one that matters most: deepening D must only REALLOCATE tail accounting, and
    // must leave the production-band mass untouched. If it changes, D is altering the in-band model
    // rather than resolving its tail, and every in-band number computed at a different D is a
    // different quantity.
    bool lower_monotone = true;      // exact mass never decreases as D grows
    bool upper_monotone = true;      // upper never increases as D grows
    bool inband_stable = true;       // mass at the PRODUCTION band d is unchanged at every D
    double inband_at_d = 0.0;        // that production-band mass, for reporting
    // Whether the FIRST iteration -- D = d, the production band -- found any state. Reported here so
    // callers do not recompute the production band separately, which doubled the dominant loop.
    bool depth1_nonempty = false;
};

TailInterval adaptive_tail_interval(const std::string& r1, const std::string& r2,
                                    const std::string& hap, std::size_t d1, std::size_t d2,
                                    const InsertPrior& ip, double log_eps, double log_1meps,
                                    double tolerance_nats, std::size_t max_depth_mult,
                                    const PieceIndex* idx = nullptr);

// ONE LEVEL of that ladder, exposed so a LAZY scorer can deepen a cell by exactly one step instead
// of recomputing every shallower level from scratch. adaptive_tail_interval is literally this
// function in a loop plus the monotonicity and in-band-stability checks, so the eager and lazy arms
// cannot drift apart: they call the same code for the same D.
struct TailLevel {
    double lower = 0.0;        // exact mass through D = d * mult
    double upper = 0.0;        // logadd(lower, bound beyond D)
    double bound = 0.0;
    double inband_at_d = 0.0;  // mass restricted to the PRODUCTION band d, recomputed at this D
    std::size_t states = 0;
    bool nonempty = false;
};
TailLevel tail_interval_level(const std::string& r1, const std::string& r2,
                              const std::string& r1rc, const std::string& r2rc,
                              const std::string& hap, std::size_t d1, std::size_t d2,
                              std::size_t mult, const InsertPrior& ip,
                              double log_eps, double log_1meps,
                              const PieceIndex* idx = nullptr);

double omitted_mass_bound(std::size_t hap_len, std::size_t m1_len, std::size_t m2_len,
                          std::size_t d1, std::size_t d2, const InsertPrior& ip,
                          double log_eps, double log_1meps,
                          const std::vector<FragmentState>& found);

// EXACT IN-BAND MASS over a state set, in the reference's own terms: for each state,
// e1 + e2 + log_pi(insert), summed in log space. reference_pair_loglik accumulates exactly
// e1 + e2 + ip.log_at(L) per (start, insert) and this mirrors it, so agreement is a statement about
// the same quantity rather than about two similar ones.
//
// NOTE ON THE COMPARISON THIS SUPPORTS. The reference integrates EVERY start and insert and does
// NOT truncate at max_divergence, so the in-band mass computed here can never equal the reference's
// total -- out-of-band states carry small but nonzero probability. The correct relation is
//     M_in_band <= M_reference <= M_in_band + M_omitted_bound
// and demanding equality with the untruncated reference would either be impossible or would force
// the bounded search to evaluate the tail exactly, which defeats its purpose.
// The untruncated reference for one fragment on one haplotype: every start, every insert in the
// prior's support, no divergence band. This is the upper side of the in-band interval.
double reference_fragment_on_haplotype(const Fragment& f, const std::string& hap,
                                       const ReferenceParams& p, const InsertPrior& ip,
                                       const std::string& r2rc,
                                       double log_eps, double log_1meps);

double fragment_states_mass(const std::vector<FragmentState>& states,
                            std::size_t m1_len, std::size_t m2_len,
                            const InsertPrior& ip, double log_eps, double log_1meps);

// ---------------------------------------------------------------------------------------------
// DIPLOID PROPAGATION AND CLASS AGGREGATION.
//
// ORDER MATTERS, and getting it wrong builds a different model:
//   1. sum every fragment's contribution for ONE FIXED representative diplotype g_i;
//   2. add that representative's exact exposure and any declared prior;
//   3. only THEN aggregate representatives of a biological class C.
// Aggregating per FRAGMENT instead would let a different full-path representative explain each
// fragment -- an implicit mosaic model, not a diplotype likelihood.
//
// The homozygous factor 2*M_a represents TWO PHYSICAL CHROMOSOME COPIES. It is unrelated to any
// log 2 from enumerating (a,b) and (b,a); the pair is unordered and is counted once.
struct MassInterval {
    double lower = 0.0;
    double upper = 0.0;
};

// One fragment's diploid contribution: the two candidate masses are combined ONCE, then mixed with
// the background ONCE. A homozygote contributes 2*M_a -- not two independently mixed background
// terms, which would double-count eta*P_bg.
MassInterval fragment_contribution(const MassInterval& ma, const MassInterval& mb, bool homozygous,
                                   double log_mix, double log_lambda, double log_bg_weight,
                                   double log_p_bg);

// One representative's TOTAL: the summed fragment contributions plus its exposure, charged ONCE
// outside the fragment loop. Exposure is -lambda * (E_a + E_b), or -2*lambda*E_a for a homozygote
// -- two chromosome copies, matching the 2*M_a factor. It belongs to the REPRESENTATIVE and must be
// added before class aggregation, because equivalent output representatives can differ in length
// and therefore in exposure.
MassInterval representative_total(const MassInterval& fragment_sum, double exposure_a,
                                  double exposure_b, bool homozygous, double lambda,
                                  double log_prior);

// CONDITION G, stated once so no caller re-derives it:
//
//     B      = max over classes C of L(C)
//     P_tau  = { C : U(C) >= B - tau }
//     certified iff |P_tau| == 1
//
// tau only ever ENLARGES the plausible set -- raising it lowers the threshold, so a class already
// in P can never leave it. It is applied in exactly ONE place; a tolerance meaning one thing in the
// plausible test and another in the tie test silently stops meaning anything.
//
// PRUNING USES THE SAME THRESHOLD. A class is safely prunable only when
//
//     U(C) < B - tau
//
// and NOT merely when U(C) < B. Pruning on the tighter test would discard classes that legitimately
// belong to the tolerance-expanded equivalence set, so the lazy arm would report a SMALLER plausible
// set than the exhaustive one -- an efficiency shortcut silently changing the answer. Since B only
// rises during refinement, a safely pruned class can never become plausible again.
inline bool safely_prunable(const MassInterval& c, double best_lower, double tau) {
    return c.upper < best_lower - tau;
}
struct Certification {
    std::vector<std::size_t> plausible;   // indices of classes that could still win
    bool certified = false;               // exactly one plausible class
    double best_lower = 0.0;
};
Certification certify(const std::vector<MassInterval>& classes, double tau);

// THREE OUTCOMES, and the third must never be read as biological ambiguity:
//   CERTIFIED  -- exactly one plausible class;
//   UNRESOLVED -- several plausible, but every one reached the tolerance or an exact endpoint;
//   INCOMPLETE -- refinement stopped early, so overlap may be budget rather than data.
enum class Verdict { Certified, Unresolved, Incomplete };
const char* verdict_name(Verdict v);
Verdict verdict_of(const std::vector<MassInterval>& classes, const Certification& cert,
                   double tol);

// Aggregate representatives of one biological class: L(C) = log SUM_i w_i exp(L(g_i)).
// Weights must sum to 1 for sequence-identical aliases, so adding a duplicate alias cannot change
// the class score -- otherwise a genotype gains confidence purely from catalogue multiplicity.
MassInterval aggregate_class(const std::vector<MassInterval>& reps,
                             const std::vector<double>& weights);

// ---------------------------------------------------------------------------------------------
// AUTHORITATIVE PATH -> BLOCK PROJECTION. The one place block coordinates are derived.
//
// Bytes are SLICED from the graph walk through the verified CandidateFrame, never rebuilt by
// concatenating alleles: concatenation is reverse-complemented for an antiparallel path and short
// for a truncated one.
//
// The COORDINATE AUTHORITY that is shared is build_candidate_frame(): this function and the
// origin/scope diagnostics both derive from it, so there is one rule for where a block sits. They do
// not all call this function -- today only --dump-scored-sequences does -- and claiming otherwise
// would overstate it.
struct PathBlockSlice {
    std::uint32_t block = 0;
    // kind/bubble_id are NOT here. They belong to the chain, and this function's one job is
    // authoritative path-to-block coordinates and bytes; the writer joins the rest.
    std::size_t walk_begin = 0;        // half-open, in AUTHORITATIVE WALK coordinates
    std::size_t walk_end = 0;
    bool reverse = false;              // the block runs antiparallel to the walk
    std::string seq;                   // in REFERENCE/BLOCK orientation, not walk orientation
    std::string md5;                   // of seq as emitted
    std::string canonical_md5;         // md5(min(seq, revcomp(seq))) -- orientation-independent
    long catalogue_allele = -1;        // this arm's numeric index, or -1 for NA. PROVENANCE, not
                                       // identity: exclusion renumbers it.
    bool catalogue_representable = false;
};

struct PathProjection {
    bool ok = false;                   // false => unprojectable; every block field is NA
    bool partial = false;
    bool reverse_frame = false;
    std::vector<PathBlockSlice> blocks;
    // Half-open walk intervals belonging to NO block, for a partial terminal frame. Emitted
    // explicitly; nothing is silently assigned to the nearest block.
    std::vector<std::pair<std::size_t, std::size_t>> unmapped;
};

// GEOMETRY AND CATALOGUE MEMBERSHIP ARE SEPARATE INPUTS, and conflating them was a real defect: a
// held-out path projected with the held-out block set had its alleles tested against that same
// held-out catalogue, where they are representable almost by construction, instead of against the
// reduced calling panel -- which is the only set that answers "could the caller have produced this".
//   projection_blocks -- contains the named path; derives the verified frame and the coordinates.
//   catalogue_blocks  -- the reduced calling panel; decides representability and supplies the index.
// For a retained path the two are the same set.
//
// Matching is by EXACT reference-oriented sequence equality, not canonical md5. Canonicalisation is
// right for identity across whole-path orientation, but inside a block catalogue that is already
// oriented, collapsing a sequence with its reverse complement would call two different alleles one.
PathProjection project_path_blocks(const std::vector<BlockAlleles>& projection_blocks,
                                   const std::vector<BlockAlleles>& catalogue_blocks,
                                   const std::string& name,
                                   const std::string& walk);

// Build the frame for one path. `walk` is the authoritative sequence supplied by the caller (from
// the shared accessor); this function only verifies it against the decomposition and derives the
// coordinate map, mirroring the offsets when the frames are opposite.

// THE RECONCILIATION GATE.
//
// Given each fragment's candidate-independent scope, score a candidate pair by summing, per
// fragment, only the origins whose blocks lie inside that fragment's scope -- and charge exposure
// ONCE for the whole locus rather than per factor.
//
// If scope is correct this must equal reference_pair_loglik for EVERY candidate pair, not merely
// differ from it by a constant: a scope that contains all the mass loses nothing, so there is no
// constant left to absorb. Recruitment cropping failed exactly here, dropping origins whose
// availability depended on phase (measured: 2 x 63.4905 nats on two fragments).
//
// `scopes[i]` is fragment i's scope, from enumerate_fragment_origins over the whole candidate set.
// `offsets_a` / `offsets_b` give the block start offsets along each homologue.

// The block span an origin touches, ALWAYS ordered low..high.
//
// For an antiparallel candidate the block index decreases along walk coordinates, so a raw
// (block_of(start), block_of(end)) pair comes back reversed -- and a reversed interval makes every
// `for (b = lo; b <= hi; ++b)` loop visit nothing, so the origin is silently exempt from the scope
// test rather than failing it. Every caller goes through this one helper: the same rule
// reimplemented in four places is how the two previous frame bugs happened.

// ---------------------------------------------------------------------------------------------
// EVIDENCE OWNERSHIP -- the rule the hybrid block caller rests on.
//
// Marker unary factors and fragment linkage factors must NOT count the same read evidence twice, so
// every fragment is OWNED by exactly one factor:
//
//   Unary(A)        depends on ONE variable block A -> feeds A's marker unary;
//   Linkage(A,B)    depends on TWO variable blocks -> feeds their linkage factor ONCE, sequence
//                   evidence included;
//   Wide            depends on THREE OR MORE variable blocks -> needs a wider factor. It is NEVER
//                   cropped into narrower factors: tests/genotype_frag_factorisation.sh already
//                   established that cropping deletes real mass;
//   Invariant       depends on NO variable block. It still carries depth and still consumes
//                   normalisation, so it is owned rather than discarded, but it supplies no
//                   genotype evidence to any factor;
//   Unusable        no certified scope at all, or mass that belongs to no block.
//
// ARITY IS COUNTED IN VARIABLE BLOCKS, NOT PHYSICAL BLOCKS CROSSED, and the distinction is not
// cosmetic -- getting it wrong throws away most of the real linkage evidence, because real bubbles
// commonly have reference sequence between them. A block whose sequence is IDENTICAL across every
// candidate carries no genotype state and no phase decision: a fragment spanning
//
//     variable A -- fixed backbone -- variable B
//
// still defines the pairwise factor psi(A, B). The backbone enters as sequence CONTEXT and nothing
// more. This is the distinction target-chain ranks and chain_span_sequence() already embodied:
// invariant intervening sequence is carried through, variable intervening sequence must become an
// explicit target or a wider factor.
//
// A fragment is genuinely Wide only when three or more VARIABLE blocks are involved. An intervening
// block whose sequence varies is itself variable, so it lands in the variable scope and pushes the
// arity to three by construction -- there is no separate rule for it, and no way for it to be
// silently carried as context.
//
// Ownership is a PARTITION, which is why "counted twice" and "silently dropped" are one assertion.
//
// CANDIDATE INDEPENDENCE is load-bearing and is why the scope is computed over the WHOLE candidate
// set at once, never over the pair being scored. A candidate may reweight an origin; it must never
// decide which variables a factor depends on, or the factor topology becomes another
// genotype-dependent approximation -- the defect this whole line of work exists to remove.
//
// WHY THE BAND AND NOT THE EXACT ORACLE. enumerate_fragment_origins is O(|hap|) per candidate and
// exists to certify, not to run: the exact emission is finite at every position, so every fragment
// has an origin in every block and a union-of-spans scope is the whole locus. Inside the production
// band only real placements survive, and what lies outside it is bounded rather than ignored --
// so the scope below is certified, with `dropped` reporting exactly what restricting to it costs.
enum class OwnerKind { Unary, Linkage, Wide, Invariant, Unusable };
const char* owner_kind_name(OwnerKind k);

struct FragmentOwner {
    OwnerKind kind = OwnerKind::Unusable;
    // The VARIABLE blocks the factor is over. Unary: lo == hi. Linkage: lo < hi, and every block
    // strictly between them is fixed -- they need NOT be physically adjacent.
    std::uint32_t block_lo = 0, block_hi = 0;
    std::vector<std::uint32_t> scope;           // physical blocks touched, certified, ascending
    std::vector<std::uint32_t> var_scope;       // the variable subset; arity comes from THIS
    double in_band = 0.0;         // logsumexp over in-band origins, all candidates
    double omitted_bound = 0.0;   // certified bound on everything outside the band
    double unmapped = 0.0;        // in-band mass belonging to NO block (partial frames)
    double dropped = 0.0;         // NATS the scope restriction costs; certified <= scope_tol
    std::size_t origins = 0;
    bool certified = false;
};

// One fragment's owning factor. `pidx`, when given, is the per-candidate piece index the bounded
// search uses; it changes speed only.
// ---------------------------------------------------------------------------------------------
// NORMALISATION, DECIDED EXPLICITLY. Partitioning the observed fragments does NOT by itself
// partition the likelihood normalisation: exposure and background are sums over states, including
// states no fragment occupies, and they can be double-charged or dropped independently of who owns
// which read.
//
// THE CHOICE MADE HERE: linkage is a CONDITIONAL PHASE SCORE, and the exposure term cancels.
//
// It cancels EXACTLY, not approximately, and the reason is worth stating because it is what makes
// the choice safe. exposure(n) = SUM_L pi(L) * max(0, n - L + 1), and for any candidate longer than
// the insert support every term is positive, so
//
//     exposure(n) = n + 1 - E[L]
//
// -- affine in length. A junction's competing configurations are phase assignments of the SAME
// allele multiset: cis = (a1b1, a2b2) against trans = (a1b2, a2b1). Both carry
// len(a1)+len(a2)+len(b1)+len(b2) in total, so E_a + E_b is identical and the difference is zero by
// construction rather than by tolerance.
//
// TWO PRECONDITIONS, both of which an earlier version of this comment got wrong.
//
// (1) THE BACKGROUND DOES NOT CANCEL, and must stay INSIDE each configuration's mixture. It is
//     genotype-independent, but it sits inside the log, so it is not a common additive term:
//
//         log(A_x + eta*P_bg) - log(A_y + eta*P_bg)  !=  log A_x - log A_y
//
//     It is a FLOOR, and it decides how much a weakly-placing fragment is allowed to say about
//     phase. Measured: for a contrast well above the floor the two forms agree to 0.0000 nats; for
//     one near the floor, dropping the background turns a true 0.0064-nat contrast into 10.0000 --
//     overstating it by 9.9936, essentially manufacturing the entire signal. This is the same
//     failure mode the band floor already exists to prevent, arriving by a different route.
//
// (2) EXPOSURE IS AFFINE ONLY WHILE EVERY CONSTRUCTED WINDOW EXCEEDS THE INSERT SUPPORT. The
//     max(0, n - L + 1) clips below that, and the cancellation fails with it. Measured against the
//     fixture's [300, 550] prior: exact and affine agree to 2e-13 at n = 549 (= hi - 1) and above,
//     and differ by 5.0e-05 at n = 548, 63.9 at n = 300 and 263.9 at n = 100. A deletion or bypass
//     allele can easily produce a short window,
//     so this is a live case and not a corner: check_exposure() below reports both values and
//     whether the regime holds, and every phase configuration must be checked. Outside the regime,
//     either enlarge the fixed flanks or retain the exact exposure difference -- never assume it.
//
// WHAT THIS DELIBERATELY AVOIDS: scoring a local junction window with the whole-locus contribution
// formula. A shorter window has a different exposure, so an absolute score computed that way would
// carry a length preference that has nothing to do with phase -- the same class of defect as the
// exposure/length cancellation still open at C4. A ratio over configurations on one fixed window
// never forms that normaliser at all.
//
// AN INTENTIONAL INFORMATION TRADEOFF, recorded rather than glossed. Excluding linkage-owned
// fragments from the marker unaries AND normalising their linkage factor conditional on endpoint
// content discards their CONTENT evidence: only phase survives. That is what makes the partition
// safe against double counting, but it is NOT a lossless likelihood factorisation.
//
// WHAT THE MASS SHARE IS NOT. ownership_ledger() pools in-band placement mass in log space across
// fragments and candidates. That is a MASS share, not information and not a likelihood
// contribution: a fragment holding a tiny share of pooled mass can still carry a decisive
// likelihood RATIO between two candidates, which is the quantity a call actually turns on. So the
// share is a size statistic and nothing more. Only the block-content regression -- C4 and the exact
// leave-zero-out controls -- can establish what information the exclusion costs, and the mass share
// must never be quoted as evidence that the loss is harmless.
//
// EVERY OWNERSHIP CLASS NEEDS A DISPOSITION, or the partition leaks somewhere unexamined:
//
//   Unary      marker CONTENT evidence, into that block's unary factor;
//   Linkage    conditional PHASE evidence, into one edge factor, content deliberately given up;
//   Wide       three or more variables. NOT representable by a pairwise transition, so it must be
//              reported UNSUPPORTED/INCOMPLETE -- never silently deleted, and never cropped into a
//              pair, which is the defect tests/genotype_frag_factorisation.sh refuted. NOT
//              "unresolved": unresolved means the model evaluated the evidence and could not
//              separate the states, a statement about the DATA, whereas Wide means the pairwise
//              model could not consume the evidence at all, a statement about the MODEL's reach.
//              Reporting the second as the first blames the data for a modelling limit -- the same
//              error as calling a depth-capped run "ambiguous";
//   Invariant  no genotype dependence, so ignorable for RANKING, but it still carries depth and
//              still belongs to absolute-fit calibration;
//   Unusable   explicitly reported MISSING evidence, not silence.

// The exact exposure and its affine surrogate side by side, so precondition (2) is CHECKED rather
// than assumed. They coincide once the window reaches hi - 1; below that the clipping bites.
struct ExposureCheck {
    double exact = 0.0;
    double affine = 0.0;      // n + 1 - E[L]
    bool in_regime = false;   // window >= hi - 1, where the two coincide and cancellation is exact
};
ExposureCheck check_exposure(std::size_t window_len, const InsertPrior& ip);
//
// `block_variable[b]` is nonzero when block b's sequence differs between candidates. A fixed block
// is context, never a factor variable.
FragmentOwner assign_fragment_owner(const Fragment& fragment,
                                    const std::vector<CandidateFrame>& frames,
                                    const std::vector<char>& block_variable,
                                    const InsertPrior& ip, double max_divergence,
                                    double log_eps, double log_1meps, double scope_tol,
                                    const std::vector<PieceIndex>* pidx = nullptr);

// The partition's headline counts, plus each class's share of pooled in-band placement mass.
//
// THE MASS SHARES ARE SIZE STATISTICS, NOT INFORMATION. They say how much in-band placement mass a
// class pools, not how much a call depends on it -- a low-mass fragment can still carry a decisive
// likelihood ratio. Named for what they measure so they cannot be quoted as a loss figure.
// EVERY class is reported, so no class can go unaccounted.
struct OwnershipLedger {
    std::size_t total = 0, unary = 0, linkage = 0, wide = 0, invariant = 0, unusable = 0;
    double linkage_fragment_share = 0.0;   // linkage-owned fragments / all fragments
    double unary_in_band_mass_share = 0.0;
    double linkage_in_band_mass_share = 0.0;
    double wide_in_band_mass_share = 0.0;
    double invariant_in_band_mass_share = 0.0;
    double unusable_in_band_mass_share = 0.0;
};
OwnershipLedger ownership_ledger(const std::vector<FragmentOwner>& owners);

// ---------------------------------------------------------------------------------------------
// THE LINKAGE POTENTIAL, psi_f, for one fragment owned by one variable-target edge (A, B).
//
// A linkage-owned fragment's mass depends on the two blocks' alleles and on nothing else, so its
// whole contribution is a table over (allele at A, allele at B). The window scored for each entry is
//
//     flank + allele_alpha(A) + FIXED CONTEXT + allele_beta(B) + flank
//
// where the context is the concatenation of the intervening blocks, every one of which is fixed --
// that is what made the fragment pairwise rather than Wide. The flanks are taken from a candidate's
// actual walk and VERIFIED identical across all candidates over the bases used; where they differ
// the fragment is not representable by this factor and is refused, never guessed.
//
// THE MIXTURE KEEPS ITS BACKGROUND. Each configuration's contribution is
//
//     log[ (1-eta) * lambda * (m(alpha1,beta1) + m(alpha2,beta2)) + eta * P_bg ]
//
// with P_bg INSIDE the log. It does not cancel between configurations -- it is a floor that decides
// how much a weakly-placing fragment may say about phase -- and dropping it turns a 0.0064-nat
// contrast into 10.0000 nats of manufactured signal.
//
// EXPOSURE IS NOT IN HERE. It belongs to the edge and is charged ONCE, outside the fragment sum;
// folding it into psi_f would multiply it by the fragment count. `exposure` below is the per
// (alpha, beta) window exposure for the edge to use once, computed EXACTLY rather than assumed to
// cancel: the affine cancellation holds only above the insert support, and a deletion or bypass
// allele can drop a window below it. Using the exact value is correct in both regimes and needs no
// precondition; `exposure_affine` records whether the regime happened to hold, so the claim stays
// auditable instead of assumed.
// ONE FRAGMENT'S HAPLOID EMISSION TABLE. This is NOT psi: it is the raw m_f(alpha, beta) a single
// homologue would produce. The edge potential still has to combine the two homologues once, mix the
// background once, sum over every fragment owned by the edge, charge exposure once, and only then
// remove the content baseline. LinkageEdge below does that; nothing here may be used as a factor.
struct LinkageEmission {
    std::size_t n_a = 0, n_b = 0;              // allele counts at A and B
    std::vector<double> mass;                  // [alpha * n_b + beta] -> log placement mass m
    // This fragment's own background, SUPPLIED BY THE CALLER on the same definition
    // fragment_contribution uses (bg_divergence over the fragment's length). It is not derivable
    // here -- bg_divergence is not a parameter of this function -- and inventing a floor locally
    // would silently give phase contrasts a different floor from genotype contrasts.
    double log_p_bg = 0.0;
    bool informative = false;   // mass actually varies with the (alpha, beta) COMBINATION
    bool ok = false;
};

// THE SHARED GEOMETRY. Window, context and flank construction lives HERE and nowhere else. It was
// briefly inlined in the genotype-frag diagnostic, which would have forced the genotype command to
// re-derive the same rule -- duplicated geometry is precisely the pattern that produced the earlier
// block-coordinate defects, and the reason build_candidate_frame became the single authority.
struct LinkageGeometry {
    std::uint32_t block_a = 0, block_b = 0;
    std::vector<std::string> alleles_a, alleles_b;
    std::string context, lflank, rflank;
    std::size_t flank_bp = 0;
    // Every window's length, [alpha * n_b + beta]. EXPOSURE LIVES ON THE EDGE, not on a fragment:
    // it is charged once per configuration, and a per-fragment copy makes accidental multiplication
    // by the fragment count easy however loudly the comments forbid it.
    std::vector<std::size_t> window_len;
    std::vector<double> exposure;              // EXACT per configuration, computed once
    bool exposure_affine = false;
    bool ok = false;   // false: flanks or context differ across candidates -- refuse, never guess
    std::string refusal;
};

LinkageGeometry build_linkage_geometry(const std::vector<CandidateFrame>& frames,
                                       const std::vector<std::vector<std::string>>& block_alleles,
                                       std::uint32_t block_a, std::uint32_t block_b,
                                       std::size_t flank_bp, const InsertPrior& ip);

LinkageEmission linkage_emission(const Fragment& fragment, const LinkageGeometry& geom,
                                 const InsertPrior& ip, double max_divergence,
                                 double log_eps, double log_1meps, double log_p_bg);

// THE EDGE POTENTIAL, aggregated in the one order that is a diploid likelihood:
//
//   S_e(config) = SUM_f log[(1-eta)*lambda*(m_f(a1,b1) + m_f(a2,b2)) + eta*P_bg,f]
//                 - lambda*(E(a1,b1) + E(a2,b2))
//
// combine the two homologues ONCE, mix the background ONCE, sum over fragments, charge exposure
// ONCE for the edge. Then the phase factor is a MEAN-ONE LIKELIHOOD RATIO within each
// unordered-content class:
//
//   log psi_e(c) = S_e(c) - logmeanexp over configurations of the SAME UNORDERED content
//                = S_e(c) - logsumexp + log|C|
//
// A PHASE LIKELIHOOD RATIO, NOT A CONDITIONAL DISTRIBUTION. psi multiplies INTO the Li-Stephens
// transition, which already carries a phase prior; a sum-one factor would count that normalisation
// a second time. It would also not be neutral on an uninformative edge: with S flat, sum-one gives
// log psi = -log|C|, and content classes have different cardinalities -- 1 for hom/hom, 2 for
// het/hom, 4 for het/het -- so an edge carrying NO phase information would penalise heterozygous
// content by up to log 4 = 1.3863 nats from class size alone. Mean-one centering gives exactly 0
// for every class size.
//
// GUARANTEED: an edge with no phase information leaves the Li-Stephens model unchanged.
// NOT GUARANTEED: that an informative edge never changes content ranking. Phase evidence can still
// move marginal content posteriors through its interaction with nonuniform Li-Stephens weights, so
// "neutral when phase-uninformative" is the claim, and "never changes content ranking" is not.
//
// NORMALISING PER FRAGMENT WOULD BE A DIFFERENT MODEL: it would let each fragment pick its own
// phase configuration independently, the mosaic error excluded at the diplotype level reappearing
// one level down.
// THE CONTRACT, chosen explicitly: LINKAGE IS OBSERVED-FRAGMENT PHASE EVIDENCE. Exposure is
// therefore REQUIRED to cancel within every content class, and is not carried in S at all -- so an
// edge with no fragments gives log psi == 0 by construction, not by luck of allele lengths.
//
// The alternative -- treating the ABSENCE of edge fragments as phase evidence -- would need the
// edge exposures to form a non-overlapping partition across the whole chain, or overlapping windows
// double-charge the zero-event exposure. That is not established here, so it is not assumed.
//
// WHERE CANCELLATION FAILS the edge is UNSUPPORTED, and the result must say so. Exposure is affine
// only above the insert support; a deletion or bypass can drop a window below it, and then
// E(a1,b1) + E(a2,b2) genuinely differs between phases. Silently keeping such an edge would let a
// zero-fragment edge move the model, which is the property this contract exists to guarantee.
//
// AN OWNED FRAGMENT WITH NO EMISSION MAKES THE EDGE INCOMPLETE. It must never be skipped: dropping
// it would quietly shrink the evidence set and report a confident answer from less data than the
// ownership partition claims.
// ONE STATUS, and usability DERIVED from it. There were two booleans -- `ok` meaning "computed"
// and `usable` meaning "safe to consume" -- and an unsupported edge deliberately set ok=true,
// usable=false. Any consumer testing the wrong one silently consumes a refused edge. A single enum
// removes the confusion rather than testing for it.
enum class LinkageStatus {
    NotComputed,
    Ok,
    ExposureDoesNotCancel,   // a short window broke the cancellation the contract requires
    InvalidEmissions,        // an owned fragment produced no emission -> INCOMPLETE
    TooManyConfigurations,   // the dense table would not fit; REFUSED, never truncated
};
const char* linkage_status_name(LinkageStatus s);

// THE DENSE-STATE LIMIT, and why it is a refusal and not a cap.
//
// A LinkageEdge is dense over n_A * n_B * n_A * n_B configurations. That is fine at fixture scale
// and impossible at locus scale:
//
//     2 x 2      ->             16 configs      0.000 GB
//     32 x 32    ->      1,048,576 configs      0.017 GB
//     64 x 64    ->     16,777,216 configs      0.268 GB
//     457 x 410  -> 35,107,516,900 configs    561.720 GB      (LPA-sized)
//
// So the dense representation CANNOT claim to handle all loci, and a sparse or on-demand
// active-allele representation is required before real loci. Until it exists this REFUSES with
// TooManyConfigurations rather than truncating, because a silent truncation here would let the
// marker top-K become an uncertified linkage cutoff -- reintroducing, at the factor level, exactly
// the uncertified-shortlist defect this whole line of work exists to remove.
//
// An arbitrary diploid edge potential also destroys the Li-Stephens O(n_h^2) transition
// factorisation, making a linkage edge O(n_h^4). That is a separate cost from the table size and is
// equally a locus-scale blocker.
struct LinkageEdge {
    std::uint32_t block_a = 0, block_b = 0;
    std::size_t n_a = 0, n_b = 0;
    std::size_t n_fragments = 0;     // owned by this edge
    std::size_t n_informative = 0;   // ...whose emission varies with the combination: the real
                                     // evidence entering psi. The rest consume normalisation only.
    std::size_t n_invalid = 0;       // ...owned but with no formable emission -> INCOMPLETE
    // Indexed [((a1 * n_b + b1) * n_a + a2) * n_b + b2].
    std::vector<double> score;       // S_e, the summed fragment contributions (no exposure term)
    std::vector<double> log_psi;     // MEAN-ONE phase ratio within each content class
    double exposure_asymmetry = 0.0; // worst |E(a1,b1)+E(a2,b2) - E(a1,b2)+E(a2,b1)| in a class
    bool exposure_cancels = false;
    LinkageStatus status = LinkageStatus::NotComputed;
    // THE ONLY TEST A CONSUMER SHOULD MAKE, and it must be made BEFORE touching score or log_psi.
    // On a refusal those vectors are EMPTY, not zero-filled -- TooManyConfigurations refuses
    // precisely so the allocation never happens -- while n_a and n_b are still set. Indexing them
    // by n_a/n_b without checking this reads out of bounds; that is not hypothetical, it segfaulted
    // the self-test the moment the refusal case was added.
    bool usable() const { return status == LinkageStatus::Ok; }
};

// ---------------------------------------------------------------------------------------------
// FIXTURE-SCALE HYBRID CHAIN. Ordered diploid states over a haplotype panel, with
//
//     log phi_b(s, s') = log T_LS(s, s') + log psi_b(s, s')
//
// as an UNNORMALISED edge potential. It is never row-normalised: forward-backward normalises its
// messages globally, and row-normalising T_hybrid would change the model and make the linkage at
// one edge depend on unrelated outgoing states.
//
// THE HIDDEN STATE STAYS ORDERED, (i, j), so phase survives inference. Unordered genotypes are
// formed only when reporting.
//
// THE MAPPING GATED HERE is the one place catalogue indices, path identities and ordered homologues
// can be confused:
//
//     ordered state (i, j) at block A, (i2, j2) at block B
//       -> allele_a[i], allele_a[j]      (allele index at A, via BlockAlleles::allele_of)
//       -> allele_b[i2], allele_b[j2]    (allele index at B)
//       -> config ((a1*n_b + b1)*n_a + a2)*n_b + b2
//
// with homologue 1 taking (allele_a[i], allele_b[i2]) and homologue 2 taking (allele_a[j],
// allele_b[j2]) -- in that order, with no implicit swap. A haplotype that BYPASSES a block must map
// to that block's bypass_allele, never to -1.
// ---------------------------------------------------------------------------------------------
// HAPLOTYPE -> ALLELE MAPPING, built through one validated path.
//
// THE DANGEROUS BOUNDARY IS int -> unsigned. BlockAlleles reports a haplotype's allele as an int
// with -1 for "no allele here", while ChainEdgeLinkage indexes log_psi with unsigned values. A bare
// cast turns -1 into 4294967295 and the subsequent log_psi[...] reads far out of bounds -- silently,
// because nothing downstream range-checks it. So every value is validated BEFORE conversion:
//
//   a < 0  and the block HAS a bypass allele  -> resolve to bypass_allele. A haplotype that does not
//          traverse the block is not missing data; it bypasses, which is a well-defined state;
//   a < 0  and the block has NO bypass allele -> MissingMapping. The edge is REFUSED and the result
//          INCOMPLETE. Never cast, never guessed, never defaulted to allele 0;
//   a >= n_alleles                            -> OutOfRange, refused for the same reason;
//   otherwise                                 -> converted, and it is now known to be in range.
//
// A HAPLOTYPE INDEX IS NOT AN ALLELE INDEX. The panel can hold hundreds of haplotypes over a
// handful of alleles, and the two coincide only in a fixture that was not designed to catch this.
enum class MappingStatus { Ok, MissingMapping, OutOfRange };
const char* mapping_status_name(MappingStatus s);

struct AlleleMapping {
    // Per panel haplotype. EMPTY unless status == Ok -- a refused mapping exposes no partial answer,
    // for the same reason a refused LinkageEdge exposes no table: a half-filled vector invites a
    // consumer that checked the wrong thing to index it.
    std::vector<std::uint32_t> allele;
    std::size_t n_alleles = 0;
    std::size_t n_bypass_resolved = 0;   // haplotypes that bypass and were resolved, not dropped
    MappingStatus status = MappingStatus::MissingMapping;
    std::size_t first_bad_haplotype = 0;
    int first_bad_value = 0;
};

// `allele_of_block[h]` is BlockAlleles' int for haplotype h, -1 where it found none.
AlleleMapping build_allele_mapping(const std::vector<int>& allele_of_block,
                                   std::size_t n_alleles, int bypass_allele);

// ---------------------------------------------------------------------------------------------
// HYBRID COMPLETENESS.
//
// THREE SITUATIONS THAT MUST STAY DISTINCT, and the third is the one that gets lost:
//
//   1. no linkage-owned evidence at this edge  -> a valid factorised (legacy) edge;
//   2. linkage evidence successfully represented -> a valid linked edge;
//   3. linkage or Wide evidence EXISTS but cannot be represented -> INCOMPLETE / UNSUPPORTED.
//
// An inactive ChainEdgeLinkage cannot express the third: the kernel reads it exactly like case 1,
// so a refused edge would silently become a neutral one and the posterior would be presented as
// complete while a reduced evidence model actually ran. Refusal is therefore carried HERE, outside
// anything the kernel sees, and a refused edge is never handed to it at all.
//
// NOT "UNRESOLVED". Unresolved means the model evaluated the evidence and could not separate the
// states. Here the model did not evaluate it -- the same distinction Wide already carries.
//
// THE INVARIANT:
//
//     hybrid COMPLETE  <=>  every non-invariant owned fragment has a SUPPORTED CONSUMER
//
// Unary fragments are consumed by their block's marker unary; Linkage fragments by their edge, but
// only while that edge is usable; Invariant fragments need no consumer, since they carry no
// genotype evidence. Wide and Unusable fragments have no consumer in a pairwise model at all.
// Block-content calls may still be emitted -- the marker unaries are unaffected -- but the PHASE
// across an affected edge, and the hybrid call as a whole, are INCOMPLETE.
struct EdgeStatusEntry {
    std::uint32_t block_a = 0, block_b = 0;
    LinkageStatus status = LinkageStatus::NotComputed;
    std::size_t n_fragments = 0;
};

struct EdgeRefusal {
    std::uint32_t block_a = 0, block_b = 0;
    LinkageStatus status = LinkageStatus::NotComputed;
    std::size_t n_fragments = 0;   // owned fragments left without a consumer by this refusal
};

struct HybridCompletenessReport {
    // OWNERSHIP completeness only: every fragment has a valid disposition. This is NOT the call
    // status -- an ownership-complete run can still fail to build its factors, and reporting this
    // field as "complete" is how an INCOMPLETE call comes to look finished.
    bool ownership_complete = false;
    std::size_t owned_total = 0;
    std::size_t invariant = 0;              // need no consumer
    std::size_t consumed_unary = 0;
    std::size_t consumed_linkage = 0;
    std::size_t unconsumed_wide = 0;        // three or more variables: no pairwise consumer
    std::size_t unconsumed_refused_edge = 0;
    std::size_t unconsumed_unusable = 0;
    // EVERY refusal, not the last one. A single "reason" field loses the others exactly when
    // several edges fail for different reasons, which is when the report matters most.
    std::vector<EdgeRefusal> refusals;
    // The variable scope of Wide fragments, so an unsupported fragment is reported with what it
    // actually depended on rather than as an anonymous count.
    std::vector<std::vector<std::uint32_t>> wide_scopes;
};

HybridCompletenessReport assess_hybrid_completeness(const std::vector<FragmentOwner>& owners,
                                                    const std::vector<EdgeStatusEntry>& edges);

// THE ONE PLACE A KERNEL EDGE IS BUILT. A refused LinkageEdge, or either endpoint's mapping having
// been refused, yields an INACTIVE entry carrying no table -- so a refused edge can never reach the
// kernel, which has no way to tell one from a legitimately linkage-free edge. The refusal is
// carried by HybridCompletenessReport instead, where it makes the run INCOMPLETE.
ChainEdgeLinkage make_kernel_edge(const LinkageEdge& edge,
                                  const AlleleMapping& map_a, const AlleleMapping& map_b);

// ---------------------------------------------------------------------------------------------
// EVERY PARAMETER THE LINKAGE FACTORS USE, in one object, reported with the result.
//
// These came from literals scattered through the caller -- 0.10, 0.05, 0.05 -- and a real-data
// result obtained that way is not interpretable, because nothing records what model produced it.
//
// THEY DO NOT ALL COME FROM THE SAME PLACE, and grouping them must not imply they do:
//
//   lambda           FRAGMENT-START INTENSITY. NOT the marker model's lambda_hap, which is a
//                    per-marker depth in different units; copying that across would be a unit error
//                    wearing the appearance of plumbing. Supplied explicitly, or left at the
//                    documented default -- and never fitted to the winning candidate, which would
//                    make the parameter depend on the answer;
//   bg_divergence    emission model: the background's implied per-base disagreement;
//   outlier_mix      emission model: the weight eta on that background;
//   error_rate       alignment/emission: per-base edit probability;
//   max_divergence   alignment: the band a placement must fall within;
//   fragment_len/sd, discordant_rate, insert_sigmas
//                    library geometry, inferred or supplied, from which the insert prior is built.
struct HybridLinkageParameters {
    double lambda = 0.05;
    bool lambda_supplied = false;      // false: the documented default, recorded as such
    double bg_divergence = 0.10;
    double outlier_mix = 0.05;
    double error_rate = 0.001;
    double max_divergence = 0.05;
    double fragment_len = 350.0;
    double fragment_sd = 50.0;
    double discordant_rate = 0.01;
    int insert_sigmas = 4;
};

// ---------------------------------------------------------------------------------------------
// TRANSACTIONAL ACTIVATION.
//
// Subtracting linkage-owned fragments from the marker unaries and activating their edges are ONE
// transaction. Doing the first without the second makes those fragments disappear from BOTH models:
// removed from the marker counts, and not consumed by any edge because the edge was refused. The
// run would then be quietly weaker than the legacy caller it is meant to extend.
//
// So: build ownership and every required edge, validate mappings and completeness, and only if the
// hybrid model is COMPLETE subtract and activate. If it is not, nothing is subtracted, no edge is
// activated, and the untouched legacy block-content call is what gets reported -- with the hybrid
// result labelled INCOMPLETE and never presented as the primary call.
//
// THE EQUALITY, phrased against ACTIVE edges rather than merely linkage-owned fragments, so a
// refused edge cannot satisfy the exclusion side by accident:
//
//     {excluded fragments} == {fragments consumed by ACTIVE linkage edges}
//
// and each excluded fragment is consumed exactly once.
// FOUR DISTINCT CONDITIONS, kept apart because collapsing them is how an INCOMPLETE call comes to
// look finished. The mapping_refused case is ownership-complete yet factor-incomplete, and its final
// status must be INCOMPLETE:
//
//   ownership_complete   every fragment has a valid disposition;
//   factors_buildable    every required edge and mapping was constructed;
//   hybrid_activated     evidence subtraction and linkage activation actually committed;
//   call_status          COMPLETE only when all of the above hold.
enum class HybridCallStatus { Complete, Incomplete };
const char* hybrid_call_status_name(HybridCallStatus s);

struct HybridActivation {
    bool ownership_complete = false;
    bool factors_buildable = false;
    bool hybrid_activated = false;
    HybridCallStatus call_status = HybridCallStatus::Incomplete;
    HybridCompletenessReport report;
    std::vector<std::string> excluded_fragments;   // empty unless hybrid_activated
    std::vector<ChainEdgeLinkage> kernel_edges;    // all inactive unless hybrid_activated
    std::size_t active_edges = 0;
    std::size_t consumed_fragments = 0;            // by active edges
    std::string refusal;                           // why the hybrid did not activate
};

// `edges` and `maps` are indexed by the edge's (block_lo, block_hi); `maps[b]` is block b's
// haplotype -> allele mapping. `n_blocks` sizes the returned kernel edge vector.
HybridActivation plan_hybrid_activation(const std::vector<Fragment>& fragments,
                                        const std::vector<FragmentOwner>& owners,
                                        const std::vector<EdgeStatusEntry>& edge_status,
                                        const std::map<std::pair<std::uint32_t, std::uint32_t>,
                                                       LinkageEdge>& edges,
                                        const std::vector<AlleleMapping>& maps,
                                        std::size_t n_blocks);

struct HybridEdge {
    bool has_linkage = false;
    std::size_t n_a = 0, n_b = 0;
    std::vector<std::uint32_t> allele_a;   // per panel haplotype: allele index at the LEFT block
    std::vector<std::uint32_t> allele_b;   // per panel haplotype: allele index at the RIGHT block
    std::vector<double> log_psi;           // from LinkageEdge; empty when has_linkage is false
};

struct HybridChain {
    std::size_t n_hap = 0;
    std::size_t n_blocks = 0;
    double recomb = 0.0;                          // Li-Stephens switch probability r
    std::vector<std::vector<double>> log_emission; // [block][i * n_hap + j], ordered states
    std::vector<HybridEdge> edges;                 // edges[b] joins block b-1 to b; edges[0] unused
};

struct HybridPosterior {
    // UNNORMALISED. Initial ordered states carry weight 1, not a uniform prior 1/n_hap^2, so this
    // exceeds a prior-carrying model's partition by exactly 2*log(n_hap). Block marginals are
    // unaffected and the brute-force oracle uses the same convention, so the two agreeing is a real
    // check. Named for what it is, because an absolute-fit or cross-panel use would be wrong.
    double log_partition_unnormalised = 0.0;
    std::vector<std::vector<double>> log_marginal;  // [block][i * n_hap + j], normalised
    // Which kernel edge paths actually RAN. A fixture must exercise both in one chain, or the
    // factorised branch or the linked branch can be present and never taken.
    std::size_t factorised_edges = 0;
    std::size_t linked_edges = 0;
    bool ok = false;
};

// Forward-backward over the ordered diploid states.
HybridPosterior hybrid_forward_backward(const HybridChain& chain);

// THE ORACLE. Enumerates EVERY ordered diploid state path and sums
//   SUM_b log E_b(s_b) + SUM_b log T_LS(s_{b-1}, s_b) + SUM_b log psi_b(s_{b-1}, s_b)
// independently of the recursion. Exponential in the chain length, so it is for fixtures only --
// but it is the only check that catches linkage applied at the wrong edge, either factor applied
// twice, an accidental row-normalisation, an ordered-state mapping error, or a correct best call
// produced from wrong posterior mass. A best-path comparison alone catches none of those.
HybridPosterior hybrid_bruteforce(const HybridChain& chain);

// `max_configs` bounds the dense table; 0 means the default. Exceeding it is a REFUSAL.
LinkageEdge aggregate_linkage_edge(const std::vector<LinkageEmission>& emissions,
                                   const LinkageGeometry& geom, double lambda,
                                   double log_mix, double log_bg_weight,
                                   std::size_t max_configs = 0);

void write_ownership_table(const std::string& path,
                           const std::vector<Fragment>& fragments,
                           const std::vector<FragmentOwner>& owners);

// ---------------------------------------------------------------------------------------------
// PLACEMENT-DEPENDENCY ORACLE
//
// The blocker (tests/genotype_frag_factorisation.sh) established that recruitment scoping is not a
// valid factorisation: a fragment recruited to one block can hold likelihood mass at another, and
// cropping deletes it. What replaces recruitment is the fragment's ORIGIN UNIVERSE.
//
// A correct generative likelihood SUMS over alternative origins -- P(f|G) = SUM_z P(f,z|G) -- so an
// N-copy array legitimately offering ~N origins is evidence, not a bug. The defects are narrower:
// top-k discarding origins without preserving their mass; cropping removing origins whose
// availability depends on another block or on phase; and summed origin mass not balanced by the
// matching exposure. So this is built to PRESERVE CALIBRATED MASS, not to suppress origin count.
//
// CANDIDATE INDEPENDENCE is the load-bearing property. The universe and the scope derived from it
// are enumerated over the whole candidate set, never over the pair being scored. A candidate may
// enable, disable or reweight a PREDEFINED origin; it must never decide which variables a factor is
// allowed to depend on, or the factor topology becomes another genotype-dependent approximation.
struct FragmentOrigin {
    std::uint32_t hap = 0;          // index into the candidate list this origin lies on
    std::int64_t start = 0;         // leading mate's start
    std::int32_t insert = 0;        // insert length L; the state is (start, L, orientation)
    bool fwd = true;
    double emission_ll = 0.0;       // both mates
    double insert_ll = 0.0;         // log pi(L)
    std::uint32_t block_lo = 0;     // blocks this origin's span touches, in chain order
    std::uint32_t block_hi = 0;
};

struct OriginUniverse {
    std::vector<FragmentOrigin> origins;
    double exact_lse = 0.0;               // logsumexp over every origin, every candidate
    // Block variables ABLE TO CHANGE the sum, which is not the same as "touched by some origin".
    // The exact emission is finite at every position, so every fragment has an origin in every block
    // and a union-of-spans scope is always the whole locus -- measured, and it made the first
    // version of this oracle useless. A block is in scope when deleting every origin that touches it
    // moves the exact logsumexp by more than `scope_tol` nats.
    std::vector<std::uint32_t> scope;
    double retained_lse = 0.0;            // mass the accelerated representation keeps
    double omitted_lse = 0.0;             // mass it drops
    // Mass from origins that belong to NO block, because they fall outside a partial frame's
    // verified window. Reported SEPARATELY rather than folded into exact_lse silently: a
    // reconciliation that simply always includes it is zero by construction and proves nothing
    // about whether a block-level model can represent it. Two candidates with identical block
    // alleles can differ here, and then their block genotype is NOT determined by block evidence.
    double unmapped_lse = 0.0;
    // The bound ACHIEVED, verified jointly: the largest |contribution(full) - contribution(scope)|
    // over candidates, after restricting to `scope`. Testing blocks one at a time does not bound
    // their combined removal -- ten blocks can each move the contribution by under scope_tol while
    // together they move it far more -- so the scope is grown until this actually holds.
    // The residual is DIPLOID: the scorer sums both haplotypes' mass into one mixture before
    // taking the log, so a guarantee proved one candidate at a time does not imply it. Both
    // figures below are the worst over all candidate PAIRS, homozygotes included.
    double initial_bound = 0.0;   // before any greedy add-back
    double achieved_bound = 0.0;  // after
    std::uint32_t blocks_added = 0;
    std::uint32_t worst_pair_a = 0, worst_pair_b = 0;
    bool bound_holds = false;
};

// One fragment's complete origin universe over `candidates`. `block_offsets[h]` gives, for candidate
// h, the cumulative start offset of every block along that candidate, so an origin's span maps to
// block variables. `retain_topk` mirrors --placement-topk so the omitted mass can be reported; 0
// keeps everything.
OriginUniverse enumerate_fragment_origins(
    const Fragment& fragment,
    const std::vector<CandidateFrame>& frames,   // sequence AND its verified block map, together
    const ReferenceParams& params,
    std::size_t retain_topk,
    double scope_tol = 1e-6);

// `include_unmapped` selects WHICH question is being asked, and the two answers differ:
//   true  -- the EXACT likelihood. Unattributable mass is retained, so this reconciles to the
//            whole-locus reference. That equality is necessary and, on its own, uninformative:
//            always including the mass makes the residual zero by construction.
//   false -- the BLOCK-LEVEL likelihood, the part a block-factored model could actually express.
//            Two candidates projecting to the same block alleles MUST score equally here. If their
//            exact scores differ while these agree, the difference is a global/unprojectable
//            component and the block genotype is ambiguous, not confidently determined.
double scope_restricted_pair_loglik(
    const CandidateFrame& frame_a,
    const CandidateFrame& frame_b,
    const std::vector<Fragment>& fragments,
    const std::vector<std::vector<std::uint32_t>>& scopes,
    const ReferenceParams& params,
    bool include_unmapped = true);

void write_fragment_results(
    const std::string& out_prefix,
    const std::vector<BlockFragmentResult>& results,
    bool have_truth);

} // namespace panvar
