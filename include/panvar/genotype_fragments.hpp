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
    int allele1 = -1;
    int allele2 = -1;
    double posterior = 0.0;
    int truth_a = -1;
    int truth_b = -1;
    bool exact = false;
    bool truth_representable = false;
};

// Per shortlisted haplotype, what each channel says about it on its own. A pair score is a sum of
// tens of thousands of terms and nothing about it is checkable from the sum: whether a haplotype
// lost on sequence or on depth, whether its reads placed at all, whether the coarse stage even gave
// it a chance. These are the columns that answer that.
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

struct HaplotypeResult {
    std::vector<std::string> shortlist;         // haplotype names actually scored
    std::vector<HaplotypeScore> haplotypes;
    std::vector<HaplotypePairScore> top_pairs;
    std::size_t n_fragments = 0;
    std::size_t n_informative = 0;
    std::vector<BlockProjection> blocks;
    std::vector<HaplotypeProbe> probes;
};

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

void write_fragment_results(
    const std::string& out_prefix,
    const std::vector<BlockFragmentResult>& results,
    bool have_truth);

} // namespace panvar
