# Higher-order inference: elimination scopes, written before allocation

## The exact model

    P  =  prod_b U_b(X_b)  *  prod_b T(X_{b-1}, X_b)
          *  F1( c2(X_2), c3(X_3), c4(X_4), c5(X_5) )
          *  F2( c4(X_4), c5(X_5), c6(X_6) )

`X_b` is the ORDERED PAIR of panel-template identities in force at block b. `c_b` maps a template
pair to the signature-class pair the factor is indexed by. `U_b` is the marker unary, `T` the
Li-Stephens transition, and the eight retained pairwise linkage factors are absorbed into their own
`T`-adjacent terms exactly as today.

## Two facts that constrain every plan

**F1 is genuinely fourth-order.** Li-Stephens permits a chromosome to copy a different template at
each block, so `X_2 .. X_5` are four distinct variables. The allele tuple stays well defined under
switching, but the factor still couples four time steps. Grouping does NOT reduce its order.

**Classes are sufficient to look the factor up, not to carry a message.** `T` decomposes as

    T = (1 - r) I  +  (r / n) 1 1^T

The stay term is diagonal in TEMPLATE IDENTITY. Two templates inside one signature class may have
different continuation probabilities, so a message collapsed to classes loses the information the
stay term needs. Any collapse must come with a contraction proof, not an assumption.

## Sizes, named honestly

  factor_lookup_configurations   F1  (2*16*10*8)^2 = 6,553,600      F2  144^2 = 20,736
  factor-label separator {4,5}   10*8 = 80 haploid, 6,400 ordered diploid ALLELE assignments

F1 distinguishes all 10 block-4 and all 8 block-5 alleles while F2 uses 6 and 8, so F1 refines F2
there and the common refinement is F1's own -- a common refinement can never be coarser than either
input. These are FACTOR-LOOKUP sizes. The exact inference separator may still require structured
template information and is not yet derived.

## Elimination scopes, symbolically

Order the blocks 0..10. Eliminating left to right, the running scope after block b is the set of
variables still needed by an un-applied factor:

  after X_2   {X_2}                     F1 pending, needs c2
  after X_3   {X_2, X_3}                F1 pending
  after X_4   {X_2, X_3, X_4}           F1 pending; F2 begins, needs c4
  after X_5   {X_2..X_5} then F1 APPLIES and X_2, X_3 drop; {X_4, X_5} remain for F2
  after X_6   {X_4, X_5, X_6} then F2 APPLIES and all drop

So the maximal elimination scope is four consecutive block variables, each an ordered template pair.
Over raw template identities that is (131^2)^4 -- not constructible. Over signature classes it is
6,553,600, but that is exactly the collapse the stay term forbids without proof.

## The structured contraction to derive

Expanding `T = (1-r)I + (r/n)11^T` across the span turns the sum into run-partition terms, each
partitioning the span into RUNS of constant template. The count is `4^(span-1)`, NOT `2^(span-1)`:
the expansion is per homologue, and the ORDERED DIPLOID state has four stay/switch combinations per
edge. F1's four-block span therefore has 64 joint transition-component patterns, not eight. Within a run the template is fixed, so that
run's class values are determined by ONE template identity, and the sum over templates inside a run
can be grouped by the class tuple that template induces over the run's blocks. Between runs the
switch term is rank one and contributes a grouped sum.

The message therefore wants to be a mixture: exact template identity where a run is still open, and
grouped class history where it has closed.

THE TWO HOMOLOGUES CANNOT BE CONTRACTED INDEPENDENTLY AND SQUARED. Marker unaries couple them at
every block, F1 and F2 couple their class histories, and the stay components retain exact template
identity on each side. So `(131 * 2,560)^2` is a WARNING BOUND on how large a naive joint
representation could get -- not a justified message layout. The term-by-term stay/switch expansion
may avoid ever materialising that product, because a run that has CLOSED can be summed immediately
rather than carried. The quantity to predict before allocating is the largest
augmented message -- its exact-template entries, its class-history entries, the operation count and
the bytes -- and to refuse if it exceeds the limit.

That derivation is NOT done. Until it is, the only sound implementation is the unoptimised one on a
small fixture, checked against complete brute force.

## Gate order

0. THE FIXTURE ITSELF must contain, or the oracle proves little:
   a four-block factor and an overlapping three-block factor SIMULTANEOUSLY;
   at least three templates;
   several templates mapping to ONE signature class while carrying UNEQUAL unary weights;
   non-identity class mappings;
   asymmetric marker emissions, plus a separate SYMMETRIC homologue-swap arm;
   r = 0, an intermediate value, and r = 1;
   non-neutral F1 and F2, plus zero-factor controls;
   the unnormalised partition weight AND every block marginal compared.
   The brute-force oracle enumerates complete ordered-template-pair histories independently. It may
   share the already-validated factor tables; it may NOT share the proposed contraction.
1. unoptimised higher-order inference on a small fixture, against complete brute force
2. the structured contraction, with its own proof and mutations
3. real-C4 prediction: largest augmented message, exact-template entries, class-history entries,
   operations, bytes -- refusing before allocation if over limit
4. two elimination plans agreeing with each other and with the same independent oracle
5. the factor-assignment ledger: every unary, every LS transition, the eight retained pairwise
   factors, F1 and F2 -- each exactly once; edges 3-4 and 4-5 absent from the kernel entirely
6. the 4-5 LS transition caught separately when omitted and when doubled
7. feature-disabled and zero-higher-order-factor cases reproducing the existing kernel byte for byte

## Predicted resources for the structured contraction (derived, not guessed)

Message state inside a span, per homologue: (classes at blocks whose run has CLOSED) x (current
template identity). The OPEN run needs no history -- its classes follow from the current template --
and the run start is implicit in the history length, so it is not a separate dimension.

                              F1 {2,3,4,5}          F2 {4,5,6}
  signature classes           2 x 16 x 10 x 8       6 x 8 x 3
  peak STORED message         17,658,669 = 141 MB   634,957 = 5.1 MB
  last position, STREAMED     1,774,945,069         40,173,901
  operations across the span  1.79e9                4.08e7
  outgoing after the factor   17,161                17,161
  warning bound (131*2560)^2  900 GB                --

TWO ARITHMETIC CORRECTIONS MADE WHILE DERIVING THIS, both worth keeping:

  * A first estimate put F1 at 3e13 operations -- four orders too high -- because it assumed dense
    O(S^2) transitions. Li-Stephens factorises as (1-r)I + (r/n)11^T and the existing kernel already
    exploits it. With the factorisation the STAY component leaves (history, X) untouched while the
    SWITCH component sums over X, closes the run, appends its classes and redistributes uniformly:
    O(entering) plus O(new histories x DIP), never O(entering x DIP).

  * The peak is NOT at the last position. Its 1.77e9 states would be 14 GB, but the factor applies
    there and collapses everything to DIP = 17,161, so they are generated and consumed streaming
    and never stored. Peak STORAGE is one position earlier, at 141 MB.

What makes this fit at all is that a CLOSED run can be summed immediately rather than carried. That
is the mechanism which keeps the message off the 900 GB warning bound, and any implementation that
carries full class history alongside full template identity will hit that bound instead.

These are PRE-ALLOCATION estimates. The implementation must compute them first and REFUSE above its
limit rather than allocate and discover.

## Both elimination orders, predicted (and they differ per factor)

  factor          order    peak entries      peak MB          streamed          updates
  F1 {2,3,4,5}    fwd        17,658,669        141.3     1,774,945,069    1,792,689,543
  F1 {2,3,4,5}    bwd       110,945,865        887.6    28,227,528,265   28,339,589,595
  F2 {4,5,6}      fwd           634,957          5.1        40,173,901       40,808,858
  F2 {4,5,6}      bwd           171,610          1.4        10,056,346       10,227,956

F1 is 16x cheaper FORWARD and F2 is 4x cheaper BACKWARD **IN ISOLATION** -- and those two
preferences DO NOT COMPOSE. F1 and F2 overlap on blocks 4 and 5, so "F1 forward, F2 backward" is not
a schedule but two incompatible eliminations. A single global order must be chosen and costed whole:

  schedule       peak stored       MB   peak at        streamed          updates
  forward         17,572,864    140.6   block 4   1,757,286,400    3,769,756,870
  backward       109,830,400    878.6   block 3  28,116,582,400   56,473,144,224

The global schedule is FORWARD: 6x cheaper in memory and 15x in updates. Its 3.77e9 updates are
about twice F1's 1.79e9 alone, because F2 contributes its own span rather than coming free.

AND THE OVERLAP FORCES A REFINEMENT. Block 4 is needed by both factors with DIFFERENT class
groupings -- F1 gives it 10 classes, F2 gives 6 -- so the carried history must hold the COMMON
REFINEMENT, which is F1's 10 because F1 refines F2 there. Storing F2's 6 would discard information
F1 requires. A common refinement can never be coarser than either input.

The per-factor figures below are retained only to show WHY the orders differ. The mechanism is mechanical: history cost is
the product of class counts for the blocks whose runs have CLOSED, so the cheaper order is the one
that puts the SMALLEST class counts first. F1's classes are (2,16,10,8) and forward begins at 2;
F2's are (6,8,3) and backward begins at 3. Assuming a single order for both would have left F2 four
times more expensive than necessary, which is why both are predicted rather than one chosen.

## The adjoint sweep is the production marginal algorithm

Z is an arithmetic circuit in the forward messages and is LINEAR in every message entry, so the
reverse-mode adjoint gives all block marginals at forward cost:

    marginal(b, x) = SUM over s with X(s) = x of  msg_b(s) * bar_b(s)

with bar seeded at 1 on the last block and swept backwards through the SAME updates the forward pass
made. This replaces the independent backward contraction, which stays only as an oracle: the
backward direction builds a right-to-left history the forward direction never needs.

PRODUCTION IS TAPE-FREE. Recording one tape entry per update would be ~1.9e9 entries at C4 scale and
would dominate every other cost, so the reverse sweep RE-WALKS the forward loop over the retained
messages and RECOMPUTES each multiplier. That preserves the message-entry update count exactly; it
does NOT preserve total CPU work, which is why the auxiliary categories are counted separately.

On the six-block fixture (four arms: r=0, r=0.25, r=1, and symmetric emissions), at r=0.25:

    message-entry updates      forward   22,608
                               adjoint   22,608
                               total     45,216
    multiplier reconstructions           22,608   per sweep
    factor lookups                       19,584   per sweep
    class/history append ops             74,088   per sweep   <-- 3.3x the update count

The last line is the point of separating them: history mapping, not the multiply-add, is the largest
per-sweep operation category, and an update count alone would have hidden that by a factor of three.

The fixture's forward count is also PREDICTED before either sweep runs, as
sum over b of |msg_b| * fanout^2 with fanout = (stay ? 1 : 0) + (switch ? NH : 0), and the adjoint
must match it exactly -- an inequality means the two enumerations have drifted apart, which a
marginal check would not necessarily catch if the drift were mass-preserving.

## Parallelization: the adjoint partitions by source, the forward pass does not

The reverse sweep's outer loop is over SOURCES, and each iteration accumulates only into its own
source's adjoint while reading the next block's adjoints read-only. So threads over sources need no
reduction buffer, no locks, and NO ATOMIC FLOATING-POINT ACCUMULATION -- which would be
non-deterministic in ordering and is excluded outright.

This is asserted as a GATE, not a design intent: visiting sources in the opposite order must give
BITWISE identical adjoints, and does in all four arms. Mutating the accumulation so two sources share
one destination breaks it in three arms.

The contended direction is the FORWARD pass, where many sources reach one destination. That is where
a deterministic thread-local reduction is required, and its buffers must be counted in the memory
prediction before allocation.

## What the operation counter measures

MESSAGE-ENTRY UPDATES: one accumulation into a state of the next message -- a multiply and an add in
linear space. NOT floating-point instructions, and not comparable to a FLOP count. It is also not
total CPU work: multiplier reconstructions, factor lookups and class/history mapping operations are
counted and reported SEPARATELY, because the tape-free adjoint trades storage for recomputation and
because history mapping outnumbers the updates themselves on the fixture.

PREDICTED INFERENCE PAYLOAD: the message and adjoint arrays only. It EXCLUDES allocator overhead,
the factor lookup structures, temporary buffers, and the forward pass's thread-local reduction
arrays. It is therefore a lower bound on inference memory, not a total, and must not be quoted as
one; the benchmark reports RSS separately.

## STATUS OF THE C4 FIGURES IN THIS DOCUMENT

The per-factor and per-schedule tables above are PROSE, not program output: no code in the tree
derives them, and they cannot be re-derived from the class counts stated alongside them (17,658,669
is not a diploid square, so it does not follow from the (2,16,10,8) x 131 state model as written).
They are retained for the QUALITATIVE conclusions they support, all of which are independently
argued: forward beats backward globally, the peak is not at the last position, a closed run can be
summed immediately, and block 4 must carry the common refinement. The absolute counts and byte
figures are NOT to be quoted until the benchmark on a representative C4 slice replaces them with
measured forward time, adjoint time, actual updates in each of the four categories, peak entries,
predicted payload, and RSS.

## THE CERTIFICATION DOMAIN IS WRONG, AND EVERY OWNERSHIP CLASS INHERITS IT

`assign_fragment_owner` enumerates a fragment's origins by placing its two mates on the PANEL
haplotype frames, one at a time. The scope it certifies is therefore the set of blocks that matter
GIVEN THAT THE TRUE SEQUENCE IS A PANEL HAPLOTYPE. That is not the model the caller runs.
Li-Stephens permits a switch at every block, so the reachable set is the full Cartesian product of
per-block alleles, and it contains recombinant tuples no panel path carries.

The two domains are not close to each other. Measured on C4, enlarging factor {3,4} to {2,3,4}:

    factor      observed max, full product     observed max, panel tuples only
    3_4                    333.61 nats                        0 nats
    3_4_5                  340.10 nats                        0 nats
    4_5                         0                             0
    4_5_6                       0                             0

Block 2 is EXACTLY flat across all 105 panel-carried tuples of `{2,3,4}` -- and moves the diploid
score by 333.61 nats somewhere in the other 15,895 cells the factor scores. A scope certified on
the left column is silent about the right one. 105 of 16,000 cells for `{2,3,4}` and 105 of 128,000
for `{2,3,4,5}` are panel-carried, so the certified region is under one percent of the domain.

### What this invalidates

Every ownership class -- Unary, Linkage, Wide, and each Unusable reason -- was derived by this
machinery, so any of them may be too small: Unary may be Linkage, Linkage may be Wide, and a Wide
fragment's scope may be wider than recorded. Factor topology, exclusion sets and resource plans all
follow from the scopes, so none of them is settled either. The C4 COMPLETE run is not withdrawn --
its arithmetic is unchanged and its explicit F1 span happens to contain block 2 -- but its status is
PROVISIONAL until every fragment has been re-audited against the full model domain.

Unaffected: the legacy caller, which does not use this machinery at all, and the higher-order
recurrence's own tests, which are statements about a factorisation given its tables.

### The correction, in two layers

Conflating these two is what produced the contradiction, so they are named separately:

  STRUCTURAL DEPENDENCY   some accepted, model-supported origin depends on the block. Decided over
                          the full local allele product, never over panel tuples.
  CERTIFIED DEMOTION      the dependency exists but its complete mixture effect -- through mix(),
                          with the real background floor, summed over every owned fragment -- is
                          bounded below a declared locus tolerance, so it may be dropped.

Exact collapse to one signature class remains SUFFICIENT for demotion and is not required. What is
required is a bound that holds over the whole domain: a quantity that is only an average, or only a
per-fragment extreme, is not one. This document previously reported such a quantity as a bound; see
the note below.

### The false bound

The earlier context statistic summed, over fragments, each fragment's worst spread within ITS OWN
best group. No single configuration need realise that combination, so it was not an upper bound on
anything. On C4 it read 32.02 nats beside an actual group-wise difference of 333.61. It has been
removed rather than renamed, and the gate now uses the observed group-wise maximum over the
certification domain, doubled by the Lipschitz constant of the mean-one normalisation.
