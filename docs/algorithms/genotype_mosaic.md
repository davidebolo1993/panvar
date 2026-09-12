# Mosaic genotyper: design, implementation plan, and review gates

This is the authoritative development contract for the replacement short-read genotyper. It combines
the implementation plan, the audit of reusable code, the required correctness gates, and the LZO/LOO
evaluation plan. A phase is not complete because its code compiles or a fixture is green: it is
complete only when its named acceptance and mutation gates pass and the result has been reviewed.

`panvar genotype` is not exposed until the end-to-end caller reaches the command-integration phase.
The current production `main` branch therefore remains a Panvar build without a genotype command.

## Current state

The clean starting point is commit `1d87162`. The first model foundation is commit `fefd4f9` on
`genotype-mosaic`:

- `genotype_model` defines paired reads, a provenanced library model, fragment-specific insert
  support, and a precomputed discrete insert distribution;
- overlapping mates use `max(|r1|, |r2|)` as their physical insert floor;
- read lengths are observations, not one run-wide constant;
- `ls_hmm` computes the exact Li-Stephens prior of a fixed allele mosaic in
  `O(blocks * templates)`;
- an independent exhaustive oracle covers `r=0`, intermediate `r`, `r=1`, recombinants, unreachable
  paths, unordered homologue multiplicity, variable read lengths, and insert normalization;
- the fixture explicitly distinguishes exact marginalization from keeping only the best template
  history.

Phase 2 now provides the immutable panel/block model and its refusal/round-trip gates. The next
implementation phase is **Phase 3: strict mosaic spelling**.

Historical code is available from Git, not from the current source tree:

```text
archive/genotype-factor-2026-09-11  final higher-order-factor research state
archive/cpp-audit-2026-09-11        non-genotype audit refinements for later review
```

Inspect an old file with `git show <tag>:<path>`. Never restore or cherry-pick an entire legacy
genotype source file. Extract one primitive at a time, put it behind a new narrow interface, and add a
focused independent gate before using it.

## Objective and evidence flow

The caller is a bounded two-stage model:

```text
paired reads
   |
   +-- canonical syncmer counts + block emissions + Li-Stephens transitions
   |       |
   |       +-- deterministic top-K whole-locus mosaic diplotypes
   |
   +-- exact/certified paired-read likelihood on only those K candidates
           |
           +-- final candidate ranking -> one allele pair per block
```

The design combines two useful ideas:

- PanGenie-like proposal: graph markers and a Li-Stephens HMM can construct mosaics that are not
  complete panel haplotypes.
- Locityper-like ranking: paired reads are evaluated against complete candidate locus sequences, so
  phase, insert size, repeated placements, and sequence context enter one coherent likelihood.

The caller can recombine block alleles already represented in the graph. It does not invent a new
block allele sequence.

## What is actually genotyped

The latent call is a **pair of whole-locus mosaic haplotypes**, not a collection of independent block
genotypes. A haploid mosaic chooses one block allele at every ordered block. Its hidden panel
template may remain the same or switch at a block boundary under the Li-Stephens transition. The
workflow is therefore:

1. infer a bounded set of phased, whole-locus mosaic pairs;
2. spell both complete locus sequences, including invariant contexts;
3. score paired reads against each complete diploid locus candidate;
4. choose the winning locus candidate or a certified equivalence set;
5. project that result back to one unordered allele pair per block for output and evaluation.

Blocks are the state alphabet, legal switch coordinates, and reporting units. They are not called
independently. A fragment crossing a boundary can therefore reject a locally plausible but
whole-locus-incoherent combination.

This distinction also states the graph's representation limit precisely. Every original complete
panel path must respell byte-for-byte, so the block decomposition may not corrupt sequences already
in the panel. A new mosaic, however, is only a mechanically valid concatenation of represented block
alleles; it is not automatically a biologically observed haplotype. The switch prior penalizes
unnecessary recombination and the whole-locus fragment likelihood tests its junctions. The caller
still cannot recover sequence absent from every block allele, and it can switch only at declared
block boundaries. LOO availability ceilings must report that representation limit separately from
caller error.

### Relationship to PanGenie and Locityper

[PanGenie](https://www.nature.com/articles/s41588-022-01043-w) supplies the proposal-side idea: its
HMM follows pairs of panel haplotype paths along graph bubbles and permits recombination between
bubbles, thereby modelling the sample as a mosaic. It then obtains local genotype likelihoods from
that HMM. The replacement Panvar caller retains mosaic generation but makes the retained
**whole-locus diplotype** the unit of final comparison.

[Locityper](https://pmc.ncbi.nlm.nih.gov/articles/PMC11844405/) supplies the ranking-side idea: align
and assign paired reads to complete locus haplotypes while accounting for edit likelihood, insert
size, depth, and competing placements. Locityper selects among supplied reference haplotypes; the
Panvar candidate set may also contain off-panel mosaics assembled from represented block alleles.

The statement that PanGenie “does not use multiplicity” needs qualification. PanGenie uses observed
k-mer counts, but deliberately selects k-mers that occur at most once within an allele and nowhere
outside the bubble. That avoids ambiguous repeated origins; it does not model the repeat-copy
placement multiplicity that is essential at loci such as LPA. Panvar keeps four different quantities
separate:

| Quantity | Where it belongs | Required treatment |
|---|---|---|
| Sequence-identical graph routes | structural model | collapse to one callable allele, retain every source alias |
| Panel templates carrying an allele | structural model / prior | retain carrier count; do not mistake it for read evidence |
| Repeated marker occurrences inside an allele | proposal model | retain exact occurrence count when informative; never reduce blindly to presence/absence |
| Distinct physical placements of a fragment | final likelihood | sum all valid origins; compress only with exact log-multiplicity |

Thus the intended model is not merely “PanGenie plus Locityper.” It is a PanGenie-like, multiplicity-
aware mosaic proposal followed by a Locityper-like whole-locus comparison whose fragment likelihood
preserves repeated physical origins exactly.

## Non-negotiable statistical contract

The marker/HMM score selects candidates only. It is not added to the final likelihood: marker counts
and fragment likelihood come from the same reads, so adding both would count the evidence twice.

For a candidate with haplotypes `a,b`, the final ranking score is:

```text
exact paired-read whole-locus log likelihood(a,b)
+ exact Li-Stephens prior(a,b), recomputed independently of the beam
```

A depth or dosage term may be added only after it has its own non-duplication proof and is named
separately in the output. Version one has no such term.

The exact fragment contract is:

```text
-lambda * (E_a + E_b)
+ sum_f log((1-eta) * lambda * (M_f,a + M_f,b) + eta * P_background(f))
+ exact_LS_log_prior(a,b)
```

where:

- `M_f,h` sums every accepted physical placement of fragment `f` on haplotype `h`;
- `E_h` is the haplotype exposure under the same normalized insert distribution;
- `lambda` is fragment-start intensity, not marker depth;
- `eta` and `P_background` form the declared robust background mixture;
- a homozygote counts two chromosome copies in both mass and exposure.

The proposal can be approximate. Final scoring must be exact relative to the declared model, or emit
a certified lower/upper interval. A winner exists only when its lower bound exceeds every competing
upper bound.

Candidate selection uses the same reads as rescoring. Until omitted candidate-set mass is bounded,
report:

- `candidate_set_weight`;
- `candidate_set_probability`;
- score gap to the next retained candidate;
- candidate-set stability as K increases.

Do not call these a global posterior or standard calibrated GQ.

## Code boundaries

New work is split by contract:

| Component | Responsibility |
|---|---|
| `genotype_model` | Immutable panel, block, read, library, candidate and status types |
| `marker_model` | Canonical syncmer inventory, observations and per-block proposal emissions |
| `ls_hmm` | Li-Stephens transition definition and exact fixed-mosaic prior |
| `mosaic_beam` | Deterministic top-K proposal with compact backpointers |
| `mosaic_spelling` | Strict block-allele path spelling and coordinate map |
| `fragment_likelihood` | Candidate-neutral placement, exposure, background and diploid score |
| `genotype_command` | Minimal option parsing and orchestration only |

Shared definitions of read emission, insert probability, valid-FR geometry, exposure, fragment
contribution and log-addition live in `fragment_likelihood`; callers may not carry private copies.
Truth names, truth sequences, evaluation distances, and LZO/LOO labels may not enter product data
structures.

## Reuse audit and constraints

### Canonical syncmers: reuse

The retained `syncmer` implementation is technically coherent:

- forward and reverse-complement k-mers have one canonical code;
- non-ACGT bases break a rolling window;
- positions remain positions in the original sequence;
- `k <= 31` prevents overflow of the two-bit code.

Before accepting user parameters, add one shared validation rule:

```text
1 <= syncmer_s <= kmer_size <= 31
```

An invalid `s` must refuse, not silently produce an empty marker set.

### Marker model: proposal only

Reuse the useful properties of the old marker path, but extract and test them:

- collapse sequence-identical graph walks;
- retain informative multiplicity differences;
- score candidates against one common marker universe;
- make unordered allele-pair emissions symmetric;
- keep invariant markers separate when estimating coverage.

Syncmers from one fragment are correlated and the model is approximate. That is acceptable for
proposal ranking, not for a second final likelihood term.

Do not inherit the old default cap of 64 alleles. C4 contains blocks with 118--126 alleles. Version one
scores every block allele; if a future resource policy prunes, each candidate must report whether its
allele was scored or received a fallback, and a fallback may not masquerade as evidence.

Adjacency evidence remains disabled. The old adjacency key sorted canonical codes, omitted their
orientations, hashed into 64 bits, and was duplicated between panel and read code. Before enabling it,
create one shared collision-checkable key containing both oriented codes and their gap.

### Li-Stephens transition: reuse the definition

At boundary `b`, for `n` templates:

```text
T_b(i -> j) = (1-r_b) I(i=j) + r_b/n.
```

The two homologues transition independently. Beam scores are proposal diagnostics only. For a fixed
haploid allele path `a_0,...,a_B`, compute its prior independently:

```text
v_0(h) = I[allele(h,0)=a_0] / n

v_b(h) = I[allele(h,b)=a_b] *
         ((1-r_b) v_(b-1)(h) + (r_b/n) sum_t v_(b-1)(t))
```

Sum the final vector. For an unordered pair, add `log(2)` exactly when the two canonical allele
mosaics differ. Sequence-equivalent representations must be canonicalized before this pair-level
multiplicity is assigned.

This recurrence is implemented and exhaustively tested in `ls_hmm`. Do not replace it with the best
beam history or a constant switch penalty.

### Fragment machinery: extract primitives, not the old scorer

The archived exhaustive substitution-only fragment model is the oracle. The old accelerated
whole-haplotype API is not a production base. It deliberately contained:

- repetitive-anchor occurrence caps;
- placement top-K truncation;
- start binning and midpoint deduplication;
- disabled fallback and mate rescue;
- a best-placement rather than summed-placement mode;
- only the first optimal edlib location;
- all-pairs `O(H^2)` scoring;
- experimental depth modes.

These approximations are particularly unsafe at LPA, where physical repeat-origin multiplicity is
evidence.

Extract and independently gate the ideas behind bounded mate placement, fragment-state enumeration,
state mass, omitted-mass bounds, fragment mixture contribution, exposure, and score certification.
Do not copy the old types or command surface.

The archived code also records traps that require permanent regression fixtures:

- encoding orientation in the sign of a binned start can flip strand for negative implied starts;
- the old `FragmentState` ordering omitted edit counts despite documentation calling every field a
  key; include them or prove geometrically equal states have equal edits;
- a run-wide maximum read length incorrectly raises the insert floor for shorter fragments;
- the exact single-end reverse-orientation arm was ambiguous. Version one requires paired reads;
- several seeds rediscovering one physical placement must deduplicate, while distinct repeat copies
  with identical statistics must retain multiplicity.

## Performance design fixed before implementation

### Invert the placement index

Do not build a large position hash for every candidate and query it once per read. Instead:

1. Split both orientations of every read into `d+1` pigeonhole pieces.
2. Build `piece -> [(fragment, mate, orientation, piece_offset)]`.
3. Scan each unique candidate haplotype once with rolling packed codes.
4. Turn each hit into an implied read start.
5. Sort/unique starts per read and candidate.
6. Verify full reads with the bounded Hamming rule.

Bucket by piece length when reads have mixed lengths. Count candidate bases scanned, rather than
calling every hash lookup “work.”

### Score unique haplotypes once

K diploid candidates contain at most 2K haploid mosaics. Intern haplotypes using hash, byte equality,
and allele-path equality. Compute `M[fragment][unique_haplotype]` once, then combine only the two
columns requested by each diploid candidate. Never enumerate every pair of unique haplotypes.

### Compress statistical states exactly

Physical origins sharing `(insert length, mate-1 edits, mate-2 edits, orientation)` may be one stored
class with exact log-multiplicity. Keep separate counters for:

- physical origin count;
- stored statistical class count;
- total placement mass.

The exhaustive oracle retains full origins. Compression is accepted only after both give identical
mass.

### Use flat storage

Use flat vectors, sorting, linear reduction, and two-pointer mate joins. Do not allocate `std::map` or
one `unordered_map<key, vector<...>>` inside fragment/candidate loops. A packed sorted occurrence
vector with `equal_range` is preferred when its cost is measured lower.

### Refine only contenders

At the initial band, produce lower and upper scores for all candidates:

1. find the highest established lower bound;
2. discard candidates whose upper bound is lower;
3. deepen only candidates with overlapping intervals;
4. stop when one lower bound exceeds every other upper bound;
5. otherwise emit the equivalence set or `RESCORE_INCOMPLETE`.

Do not perform exact tail work for a candidate already proven unable to win.

### Keep the beam allocation-efficient

For C4, 131 squared ordered template states, K=256, and ten transitions imply roughly 44 million
extensions, which should be practical. Use:

- backpointers instead of copied path vectors;
- one flat score buffer per block;
- `nth_element` or partial selection followed by deterministic final sorting;
- precomputed stay/switch combinations;
- no hash map per extension.

## Implementation phases

Every phase ends in its own commit. Do not combine a refactor, new model, and real-data experiment in
one change.

### Phase 1 — model foundation: complete

Implemented at `fefd4f9`:

- fragment-specific overlapping-pair insert support;
- provenanced library parameters;
- a discrete normalized insert table built once per fragment/read-length bucket;
- exact fixed-mosaic Li-Stephens prior;
- independent exhaustive tests and input refusal.

Reviewer checkpoint: verify the exhaustive oracle is independent, the zero-probability arm is not
skipped, and the fixture distinguishes summing histories from retaining only the best.

### Phase 2 — immutable panel/block model: complete

`genotype_model` now builds one prepared model containing:

- ordered blocks and their invariant contexts;
- every block allele sequence and stable allele ID;
- panel template names;
- `template -> block -> allele` mapping;
- chain orientation and complete/partial-frame status;
- a deliberately structural-only interface; marker universe and observations remain absent until
  the structural model passes.

Acceptance gates implemented in `genotype_mosaic_core`:

1. Every complete panel template has exactly one allele at every block.
2. Duplicate names, missing mappings, inconsistent block counts, or partial frames refuse explicitly.
3. Allele IDs are stable under input record reordering.
4. A complete index round-trip recovers `(template, block, allele)` exactly.
5. No truth or experiment label is present in the prepared model.

Additional multiplicity contract: sequence-identical source walks collapse by exact sequence, but
their sorted source aliases and panel-carrier counts remain attached to the canonical allele. Empty
bypass alleles remain explicit states. Block order stays semantic; allele records and template rows
are canonicalized independently of their input order.

The graph-facing adapter and authoritative byte-for-byte panel-path check belong to Phase 3, where
the speller exists. No marker emission code may be added before that round-trip passes.

### Phase 3 — strict mosaic spelling

Create one `MosaicPath`/`HaplotypeSequence` speller. It must handle:

- chain orientation and reverse frames;
- bypass and empty alleles;
- invariant sequence between blocks exactly once;
- block-to-sequence coordinate intervals;
- recombinant paths with no graph P-line;
- explicit partial-frame refusal.

Acceptance gates:

1. Derive and respell every complete panel path; bytes must equal its authoritative sequence.
2. Spell an off-panel recombinant with a hand-computed expected sequence.
3. Cover reverse orientation, boundaries, empty alleles, and variable allele lengths.
4. Mutation: duplicate or omit an invariant context; round-trip must fail.
5. Mutation: use raw walk coordinates on a reverse frame; mirrored interval gate must fail.

### Phase 4 — marker proposal model

Extract canonical syncmer inventory and read counting into `marker_model`; do not restore the old
command. Prepare block emissions once and pass the immutable object to the candidate generator.

Acceptance gates:

1. Validate `1 <= s <= k <= 31`.
2. Panel and read paths use the same marker-key implementation.
3. Score every C4 allele; report `scored_alleles == total_alleles` for every block.
4. Emission matrices are symmetric under homologue swap.
5. Sequence-identical alleles have identical proposal emissions.
6. Candidate-independent constants cannot change ranking.
7. An independent hand-counted fixture checks marker presence and multiplicity.
8. Adjacency evidence is absent from version one.

The marker likelihood is labelled a proposal score, not a final calibrated likelihood.

### Phase 5 — deterministic top-K mosaics

At block `b`, a beam state is an ordered panel-template pair `(h1,h2)`. Extension score is exactly:

```text
block proposal emission + diploid Li-Stephens transition
```

Requirements:

- configurable unique candidate target K and explicit raw-history limit;
- deterministic tie-breaking independent of thread scheduling;
- preserve homologue orientation during extension;
- canonicalize a global homologue swap only in final diploid candidates;
- collapse candidates spelling the same two allele paths;
- log-sum retained equivalent-history proposal mass for diagnostics;
- retain best-history score separately;
- continue raw expansion until K unique mosaics or the raw limit is reached.

Candidate manifest fields:

```text
candidate_id, unique_rank, hap1_alleles, hap2_alleles,
best_proposal_score, retained_logsum_proposal_score,
proposal_emission_component, proposal_transition_component,
equivalent_retained_histories, every_allele_scored
```

Acceptance gate: exhaustive enumeration on a small locus across `r=0`, intermediate `r`, `r=1`,
ties, bypass alleles, homologue swaps, and multiple templates sharing alleles. Mutations of a
transition, tie order, swap canonicalization, and premature allele-level merging must each fail a
named assertion.

### Phase 6 — candidate-neutral fragment likelihood

Extend `LibraryModel` with explicit error, insert, background, outlier and fragment-start parameters
plus provenance. Read lengths remain fragment observations.

The API receives:

```text
unique spelled haplotypes
exactly K requested diploid pairs
paired fragments
one validated library model
```

It scores only those K pairs. Initial C4 model:

- paired reads only;
- substitutions only and fixed-position Hamming emission;
- observed read lengths and overlapping pairs;
- normalized discrete insert prior;
- both library orientations;
- every in-band placement;
- no occurrence cap, placement top-K, start bin, or midpoint deduplication;
- physical-origin multiplicity preserved;
- analytic exposure;
- candidate-independent background and lambda;
- no depth, dosage, truth-length, or LPA-specific term.

Required gates:

1. Exact in-band mass per fragment/haplotype against exhaustive scoring.
2. Exact exposure.
3. Both library orientations populated.
4. Repeat-copy multiplicity preserved.
5. One placement found through several pieces counted once.
6. Distinct origins with identical statistics counted with multiplicity.
7. A homozygote receives two chromosome copies.
8. Candidate reordering leaves raw scores unchanged.
9. One candidate alone has its identical raw score from a larger list.
10. Changing proposal scores cannot change fragment scores.
11. Score intervals contain exhaustive results.
12. A winner is emitted only when its lower bound beats all competing upper bounds.
13. Unsupported indel/non-ACGT behavior refuses with a named reason.
14. Resource exhaustion returns no partial ranking.

Permanent prior mutations must fail:

- retain only the best template history;
- omit the unordered-pair `log(2)`;
- add the proposal score to the final score;
- substitute one constant recombination penalty for the exact transition.

Runtime counters:

```text
candidate_bases_scanned, piece_hits, proposed_starts, unique_starts,
full_read_verifications, accepted_mate_placements, valid_fr_states,
physical_origins, statistical_state_classes, candidates_deepened
```

### Phase 7 — command integration

The new caller becomes the only implementation of:

```text
panvar genotype
```

Do not reintroduce `genotype-frag`, `--hybrid-call`, or `--mosaic-rescore`. Keep the public surface
small:

```text
--gfa, --bubbles, --reads1, --reads2, --reference-path,
--out-prefix, --threads, --max-candidates
```

Library parameters are estimated and recorded; a single optional model/config input may provide
declared simulation parameters. Debug probes live in unit tests or a developer executable.

Execution order:

1. load and validate the structural model;
2. count markers and construct proposal emissions;
3. generate unique mosaic candidates;
4. spell and intern candidate haplotypes;
5. recompute exact Li-Stephens priors independently;
6. compute/certify fragment likelihoods for requested pairs;
7. choose a certified winner;
8. project its two allele paths to standard per-block output.

Outputs:

```text
*.genotypes.tsv
*.mosaic_candidates.tsv
*.mosaic_scores.tsv
*.mosaic_status.tsv
*.mosaic_haplotypes.fa
```

Status is one of:

- `COMPLETE`;
- `MODEL_INCOMPLETE`;
- `CANDIDATE_INCOMPLETE`;
- `SPELLING_INCOMPLETE`;
- `RESCORE_INCOMPLETE`.

An incomplete run emits diagnostics but no ranking presented as complete.

### Phase 8 — one-donor C4 go/no-go

Use one frozen LZO donor. First require at K=64:

```text
proposal + spelling + rescoring <= 90 seconds
peak RSS <= 4 GB
```

Report wall time separately for model preparation, proposal, spelling, read-piece indexing, scanning,
verification/joining, candidate combination, refinement, and total. If the gate fails, profile the
real recurrence/scorer before changing the model or adding threads.

Then run K = 8,16,32,64,128,256, recording:

- raw template histories and unique mosaics;
- whether a sequence-equivalent truth is present;
- truth proposal and fragment-score ranks;
- winning per-block true edlib distance;
- phase-aware haplotype-path result;
- score gap and candidate-set probability;
- runtime, RSS, and winner stability between K and 2K.

Stop conditions:

- truth absent at K=256: candidate generator inadequate;
- truth present but loses after rescoring: fragment likelihood problem;
- winning mosaic correct but block calls wrong: spelling/projection problem;
- K=64 exceeds time or memory gate: profile before any donor cohort.

This—not another large instrument suite—is the first biological go/no-go result.

### Phase 9 — development pilot and held-out LZO

The eight previously inspected donors are development data, not untouched evaluation data.

Pilot rules:

- two workers initially;
- clean checkpoint signatures including binary, manifest, parameters, scripts, simulator, and seed;
- COMPLETE runs only enter accuracy tables;
- candidate manifests frozen before truth comparison;
- one K and one model for all donors;
- 30x ceiling control plus at least one nested lower coverage;
- stop if any sequence-equivalent truth is absent at K=256.

After the pilot, freeze K and all model parameters. Evaluate the other 56 donors using nested subsets
from one reproducible 30x simulation, preserving mate pairs, for example 30x, 15x, 7.5x, and 3.75x.

Compare:

- marker proposal alone;
- marker proposal plus fragment rescoring;
- if desired, an archived old-caller binary as an external research baseline, never as shipped code.

Primary result: per-block unordered diploid true edlib distance, displayed as one violin per block.
Companion results: phased path distance/switch error, truth-in-candidate rate, COMPLETE rate,
LOWCOV/no-call rate, runtime/RSS, and stability with K. Expect saturation at 30x; lower coverages are
the main accuracy experiment.

### Phase 10 — LOO, then other loci

Only after LZO is frozen, remove both donor haplotypes from every proposal source before marker
catalogue construction, template HMM construction, candidate generation, and rescoring.

Evaluate against:

1. actual truth distance;
2. availability floor: closest retained diploid block alleles after removal.

Report excess error above the availability floor. Do not penalize a caller for failing to emit an
allele removed from its model.

Once C4 LZO/LOO is accepted, proceed to the other loci except LPA. Treat LPA last because repeat-copy
multiplicity and dosage may require an explicitly new, independently tested model. Long-read
integration comes only after the short-read caller is general across loci.

## Critical-review protocol

At every phase Claude reports, in one attributable checkpoint:

1. exact commit, binary hash, dirty/clean state, and parameters;
2. what production behavior changed and what did not;
3. expected-versus-observed gate table with units and domains;
4. non-vacuity evidence for every important fixture arm;
5. each mutation, the binary/source change, and the assertion it breaks;
6. refusal behavior and proof that no partial output is presented as complete;
7. measured runtime/memory counters rather than extrapolated viability;
8. claims still untested.

The reviewer checks source as well as reports. In particular:

- a green test over zero rows is a failure of the test;
- identical paths through shared code are not independent validation;
- a threshold may not be selected after seeing the result it gates;
- comparison quantities must use the same units, normalization, candidate domain, and read universe;
- oracle manifests are frozen independently of output;
- a mutation that changes no observable is a fixture gap, not a passing mutation;
- a resource estimate is not a performance result;
- an accuracy result from an INCOMPLETE arm is excluded, not coerced into a distance;
- no cohort starts before the one-donor runtime gate passes.

I will act as this critical reviewer between phases. A request to continue should include the report
above; the next phase begins only after discrepancies are resolved or explicitly accepted as limits.

## Deliberate non-goals for version one

- higher-order fragment factors, ownership scopes, and factor planning;
- full-domain factor certification or higher-order inference;
- the `genotype-frag` command or historical hybrid/audit/probe flags;
- new syncmer definitions or adjacency evidence;
- learned parameter tuning;
- experimental depth/copy-number likelihoods;
- indel-aware emission without an independent exact placement oracle;
- single-end support;
- global posterior/GQ claims from a truncated candidate set;
- LPA-specific corrections;
- long-read integration;
- scorer multithreading before the serial algorithm meets its performance gate.

The immediate target is narrow: **one frozen C4 donor, K=64, truth-in-candidate reported, certified
fragment ranking and per-block edit distances emitted in under 90 seconds and under 4 GB RSS**.
