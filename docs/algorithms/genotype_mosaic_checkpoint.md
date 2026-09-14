# Mosaic genotyper: state audit and plan for phases 4-8

Companion to `genotype_mosaic.md`, which remains the authoritative contract. This document records
what is actually built, what its own gates still owe, and the work order to reach the Phase 8
go/no-go.

## 1. Verified state

Commit `d842deb` plus an uncommitted working tree. Build: Release, `-DPANVAR_BUILD_TESTING=ON`.
Both genotype test binaries pass.

| Phase | Component | Code | Gates |
|---|---|---|---|
| 1 model foundation | `genotype_model`, `ls_hmm` | committed `fefd4f9` | pass |
| 2 immutable panel/block | `genotype_model` | committed `715e8ed` | pass |
| 3 strict spelling | `mosaic_spelling` | committed `d842deb` | pass; graph checkpoint met at `51b6a49` |
| 4 marker proposal | `marker_model` | committed `c3b8c38` | pass; gate 3 met on real C4 at `51b6a49` |
| 5 deterministic top-K | `mosaic_beam` | committed `c3b8c38` | pass against exhaustive oracle |
| — graph adapter | `genotype_locus` | committed `51b6a49` | pass on real C4 |
| 6 fragment likelihood | `fragment_likelihood` | committed `98324c9`..`6bb382c` | pass; oracle-certified |
| 7 command integration | `genotype_command` | committed `7c5a42b` | pass end-to-end on C4 |

`panvar genotype` is now registered in `main.cpp`, which the contract permits from Phase 7.

The k-way merge over already-sorted allele lists is present in `build_marker_panel`
(`src/marker_model.cpp:198`), replacing the second N log N sort.

`tests/genotype_mosaic_proposal.cpp` is substantive, not decorative: it carries an independent
exhaustive enumerator over ordered template-pair histories and compares the beam against it at
r=0, r=0.25/0.6, and r=1, asserting best score, retained log-sum, emission and transition
components, and history counts separately. Non-vacuity is asserted explicitly
(`summed_histories_are_live`), and chunk-size variation stands in for scheduling nondeterminism.

## 2. The C4 block structure, measured

Panvar is an SV caller, so `panvar bubble`'s default `--min-variant-bp 50` is the correct block
alphabet: it keeps bubbles whose interior span is SV-sized. On `tests/real_data/c4.gfa.gz` that is
**5 SV bubbles** out of 1005 snarls in the complete decomposition.

The other 1000 snarls are all polymorphic, and every one of them lies in one of the six spans
*between* the SV bubbles — none is nested inside one:

```text
gap  0: 101 snarls    gap  6: 313 snarls
gap  2: 132 snarls    gap  8: 100 snarls
gap  4:  10 snarls    gap 10: 344 snarls
```

So at C4 no inter-block span is invariant, and `invariant_segments` — one shared string per gap —
cannot respell 131 templates byte-for-byte. Phase 3's round-trip gate forces the adapter to **tile**
the locus instead: the six context spans become blocks in their own right, alternating with the five
SV blocks. That yields **11 blocks**, which is the figure the contract's performance paragraph
assumes.

Measured per-block census after canonicalising orientation:

| Block | Kind | Alleles | Min bp | Max bp |
|---|---|---|---|---|
| 0 | context | 68 | 12291 | 12313 |
| 1 | SV | 68 | 277 | 7947 |
| 2 | context | 100 | 20784 | 20808 |
| 3 | SV | 16 | 203 | 506 |
| 4 | context | 10 | 274 | 324 |
| 5 | SV | 8 | 88 | 205 |
| 6 | context | 118 | 47756 | 47808 |
| 7 | SV | 119 | 26667 | 98508 |
| 8 | context | 62 | 13857 | 13862 |
| 9 | SV | 5 | 387 | 456 |
| 10 | context | 126 | 63407 | 63478 |

700 alleles over 11 blocks. The contract's "C4 contains blocks with 118--126 alleles" matches blocks
6, 7 and 10 exactly, which confirms the tiled decomposition is the one it was written against.

Two further structural facts the adapter must handle, both already anticipated by the contract:

- **0 partial frames.** All 131 panel paths traverse every one of the ten SV boundary nodes exactly
  once. No template is refused.
- **60 of 131 templates are antiparallel** to the reference chain. `mosaic_spelling`'s reverse-frame
  handling is load-bearing on real data, not a defensive branch.

130 of the 131 templates spell a distinct whole-locus sequence; the locus is 194,089--265,891 bp.

## 3. Defect found and fixed

`run_beam` stored each retained layer by moving `best`, which after `flush()`'s `best.swap(pool)`
owns the expansion pool's buffer — capacity `expansion_chunk_size`, default 262144 nodes. Layers are
retained for every block, so one chunk-sized allocation stayed alive per block. Storing exactly the
retained nodes instead, measured at C4's shape (131 templates, 11 blocks, r=0.01):

| K | Extensions | RSS before | RSS after | Wall |
|---|---|---|---|---|
| 64 | 11,000,201 | 121 MB | 26 MB | 0.08 s |
| 256 | 43,949,321 | 121 MB | 26 MB | 0.20 s |

Both test binaries remain green after the fix.

Two things worth recording. First, the contract's estimate — "131 squared ordered template states,
K=256, and ten transitions imply roughly 44 million extensions" — is exact: the measured figure is
43,949,321. The proposal stage is not a performance risk, and the stay/switch factorization the
contract permits is not needed for C4.

Second, the waste is proportional to block count, so the fix matters more elsewhere than here. On a
locus with a thousand blocks the same code reached 4.58 GB at K=64, over the Phase 8 memory gate on
its own. At C4's eleven blocks it is 121 MB, which no gate would have caught. It is fixed either way.

## 4. Steps 1 and 2: complete

Both landed at `51b6a49`, gated by `tests/genotype_locus_adapter.cpp` against the real C4 panel.

### What the adapter does

It is the only place where graph structure becomes statistical input: every later stage reads
`StructuralLocusInput` and never touches nodes, edges or walks. It tiles the locus — SV bubbles at
the production `--min-variant-bp`, interleaved with context blocks covering the spans between and
beyond them — and emits one source allele per distinct walk, leaving `prepare_genotype_model` to
collapse sequence-identical walks and retain their aliases.

Four decisions were carried over from the archived block chain
(`archive/genotype-factor-2026-09-11:src/genotype_blocks.cpp`), each of which the archive records as
a defect it had already paid for:

- bubbles contained in another bubble are dropped, because a chain is a partition and a contained
  bubble puts the same sequence in two blocks, after which marker confinement strips both bare;
- a bubble's `source`/`sink` are ordered by node position rather than taken at face value, because a
  reverse-oriented bubble stores them swapped and the neighbouring context spans would then run
  straight across its interior;
- abutting bubbles sharing one boundary node trim it from the later one, or the chain counts those
  bases twice;
- a template that does not walk a block **bypasses** it, with an explicit empty allele. Treating
  that as missing data makes a template whose deletion removed a block's own boundaries
  unrepresentable at that block, so it could never be called and never be scored.

### The canonicalisation defect, found by the gate

The node-split gate failed on the first run. Splitting every node in half left the concatenated
locus byte-identical and the census identical — 11 blocks, 700 alleles, 71 forward / 60 reverse —
but every block's sequence shifted, some by hundreds of bases, and the fingerprint with it.

The cause is that a bubble's source and sink are traversed by every path through it, so their bases
sit inside the variant block for every allele. Where the graph cuts those boundary nodes is
representation: a split moves bases between a variant block and its context neighbour.

The fix is to trim each variant block down to the span over which its alleles actually differ,
pushing the shared prefix and suffix outward into the context neighbours. The maximal trim is the
same whichever way the boundary node was cut, so both representations converge; the shared flanks
are identical for every allele, so no template's spelled locus changes; and the result is also the
honest definition of an SV allele.

The two representations move different amounts and still converge, which is the evidence that the
trim is doing the work: the original C4 graph moves 1,159 bases outward, the node-split graph moves
660 — because its context blocks already absorbed half of each cut boundary node — and both arrive
at the same 11 blocks, the same 700 allele sequences and the same fingerprint. Blocks carrying a
bypass allele have no common flank and are left alone; that is a remaining limitation, not exercised
at C4.

### Measured C4 census, after canonicalisation

| Block | Kind | Alleles | Min bp | Max bp |
|---|---|---|---|---|
| ctx:0 | context | 68 | 12326 | 12348 |
| sv:1 | SV | 68 | 0 | 7670 |
| ctx:2 | context | 100 | 20811 | 20835 |
| sv:3 | SV | 16 | 117 | 420 |
| ctx:4 | context | 10 | 325 | 375 |
| sv:5 | SV | 8 | 0 | 117 |
| ctx:6 | context | 118 | 47951 | 48003 |
| sv:7 | SV | 119 | 26346 | 98187 |
| ctx:8 | context | 62 | 14267 | 14272 |
| sv:9 | SV | 5 | 0 | 69 |
| ctx:10 | context | 126 | 63396 | 63467 |

736 source walks collapse to 700 canonical alleles. 71 forward and 60 antiparallel templates, zero
partial frames, zero bypasses, zero contained bubbles, zero abutting boundaries. The three SV blocks
with a zero-length allele are deletions: after trimming, the shortest allele of a deletion bubble is
empty, which is what it should be.

### Gate results

1. **Round trip.** All 131 templates respell byte-for-byte in their own frame, with the 60
   antiparallel ones canonicalised. Orientation is derived structurally, from the order in which a
   template visits the chain's boundary nodes — never by comparing bytes against the authoritative
   sequence, which would make this gate tautological.
2. **Fingerprint invariance.** Splitting every graph node of length >= 2 in half — changing every node
   identifier, every walk signature and every edge — leaves the fingerprint at `a58e2656e6209e0c`.
   The split model also passes its own round trip. Synthetic arms cover walk renaming, record
   reordering and sequence-identical aliases as invariants, and allele-sequence and block-kind
   changes as non-invariants.
3. **Block intervals tile** the spelled haplotype with no gap or overlap.
4. **Partial frames refuse**, proven synthetically since C4 has none.
5. **Phase 4 gate 3.** The marker panel covers exactly the 700 canonical alleles.

### Marker results on real C4

At k=31, s=15 over the prepared model: marker universe 47,960, of which 25,167 are invariant
anchors; **114 markers dropped for varying in more than one block** (0.24% of the universe); no block
is left without a discriminating marker; 694 of 700 alleles carry at least one.

That answers the `dropped_multiblock_markers` question with a measurement rather than a guess: the
C4A/C4B duplication does not cost the proposal model anything material, so the current semantics
stand.

The six alleles without a marker are three empty deletion alleles and three non-empty alleles
shorter than k, which carry no syncmer at all. Those three are invisible to the marker proposal and
can only be separated by the fragment likelihood. This is the same mechanism that makes GSTM1 hard,
now quantified at C4 and small.

One honesty note for the reviewer: `scored_alleles == total_alleles` holds by construction, because
version one has no allele cap anywhere in the code. On its own that equality proves nothing. The
load-bearing assertions are that the panel covers exactly the model's canonical allele set and the
per-allele marker coverage census above.

## 5. Phase 6 and Phase 7: complete

### The fragment likelihood

Four commits, in the order the contract names them.

`fragment_likelihood` owns the single definition of read emission, insert probability, valid-FR
geometry, exposure, fragment contribution and log-addition. `log_add_exp` had three copies in the
tree before this and now has one. The exhaustive `enumerate_placements_exhaustive` is certified
against hand arithmetic, not against a second implementation of the same idea, and then certifies
the accelerated scorer: every `(fragment, haplotype)` mass agrees to 1e-9 on planted substitutions,
tandem repeats, mixed read lengths and twelve random seeds.

The accelerated path inverts the index — pieces of the reads are indexed and each haplotype is
scanned once, rather than a position hash built per candidate. Piece keys are prefixes of the
pigeonhole chunks, which preserves the band guarantee because an exact chunk has an exact prefix.
Reads are bucketed by piece length so mixed read lengths do not lose placements. Unique haplotypes
are interned by bytes. Statistical classes are compressed by `(insert, edits, edits, orientation)`
with an exact log-multiplicity, which is the sum rather than an approximation of it; physical
origins and stored classes are counted separately.

There is no occurrence cap, placement top-K, start bin, midpoint deduplication or best-placement
mode anywhere in this path.

### Two findings the gates forced

**Exactness is unreachable through the band.** The scorer is exact only when the band is at least
as wide as the longest mate, and at that band the pigeonhole index has no pieces left to key on. So
the upper bound is what certifies, and the reference it is certified against is the oracle on an
independent code path — not the scorer in a different mode, which would have been circular.

**One unplaceable fragment can swallow every real candidate separation.** A fragment that places on
no candidate at the working band contributes the outlier arm to the lower bound but the entire
omitted-tail bound to the upper one. Measured on C4 at 30x: exactly **one fragment in 34,415**
placed nowhere, and that single fragment produced a **244-nat interval against an 18-nat candidate
separation** — identical to four decimals across all five candidates, so unable to reorder anything,
but wide enough that no band the index can support would ever close it.

The resolution is a declared model statement, not a relaxed proof: a fragment that places on no
candidate is an outlier, so its contribution is the outlier arm exactly and its omitted tail is
zero. That closes the interval from 244 nats to 0.52. It is named, counted in the status file as
`fragments_declared_outlier`, and can be switched off for a fully unconditional bound. The
containment gate tests the bound with the policy off; a separate gate proves the policy touches only
fragments that place nowhere, and never the score itself — only the upper bound.

### The command

`panvar genotype` with the eight declared options and no others. It runs the declared execution
order, including the byte-for-byte panel round trip on real input before any statistic is computed
from those blocks. The exact Li-Stephens prior is passed *into* certification rather than added
after it, because pruning on the fragment term alone would discard candidates the prior would have
rescued. The proposal score is not passed and may never be. Context blocks are modelled but never
reported as calls. The fingerprint is keyed on the graph's reference frame, not the output prefix.

### Dry run on real C4

One frozen donor simulated at 30x from two panel templates (`HG00096#1` and `HG00268#2`, the latter
one of the 60 antiparallel ones), 34,415 read pairs, 100 bp mates, insert 350±40, 0.2% substitution
error, declared library parameters.

| K | Status | Certified gap | Total | Peak RSS |
|---|---|---|---|---|
| 16 | COMPLETE | 17.02 nats | 5.0 s | — |
| 64 | COMPLETE | 17.02 nats | **6.3 s** | **544 MB** |

Against the Phase 8 gate of 90 s and 4 GB: 14x under time, 7.5x under memory. Stage split at K=64:
model 0.96 s, markers 0.95 s, proposal 0.37 s, spelling 0.00 s, rescoring 4.05 s. Rescoring
dominates, as expected.

**Both called haplotypes are byte-exact matches to the two simulated truth haplotypes** — 226,137
and 232,743 bp, the second in reverse complement, which is correct because that template is
antiparallel and the caller emits in canonical chain orientation.

Three things this result is *not*. The truth haplotypes are panel templates, so this is the easy
end of the problem — Phase 8 is LZO and Phase 10 is LOO, and neither is this. The reads were
simulated under exactly the declared library model. And the beam returned only 8 unique candidates
at K=64, so the candidate set is much narrower than K; whether that is the switch prior, the marker
emissions, or a real property of C4 is a question for the first LZO run.

## 6. Remaining work order

### Step 1 — Phase 8 one-donor C4 go/no-go

The runtime half of the gate already passes on the dry run above. What remains is the accuracy half:
the K = 8..256 ladder with truth-in-candidate rank, per-block true edlib distance, score gap,
candidate-set probability, and winner stability between K and 2K, on a frozen LZO donor rather than
a panel member.

### Step 2 — Phases 9 and 10

The development pilot, the held-out LZO cohort, then LOO. Unchanged from the contract.

### Known limitations, carried forward

- **Library parameters are declared, not estimated.** `--library` is required. The contract asks for
  estimation with declaration as an option; version one has only the declaration. For a simulated
  donor this is the right input anyway, but it is a real gap before any real sequencing run.
- **`--threads` is recorded and not used.** The contract forbids scorer multithreading before the
  serial algorithm meets its performance gate. It now does, so threading is unblocked but unbuilt.
- **Blocks carrying a bypass allele are not divergence-trimmed**, so node-split invariance is not
  proven for them. C4 has no bypasses.
- **The beam's bypass-allele arm is exercised only through the adapter**, not in the beam's own
  exhaustive-oracle fixture.

## 7. Reviewer questions

1. **Tiling.** Confirm that promoting the six polymorphic context spans to blocks is the intended
   reading of "ordered blocks and their invariant contexts". At C4 no gap is invariant, so the
   alternative — one shared context string per gap — cannot pass Phase 3's byte-for-byte gate, and
   would inject systematic false bases into every spelled candidate that the fragment likelihood
   would then have to absorb. Related: confirm context blocks are modelled but excluded from
   per-block genotype output, since they are not SV calls.
2. **Long blocks and switch resolution.** SV block 7 reaches 98.5 kb and context block 10 reaches
   63.5 kb. Inside such a block the mosaic cannot recombine at all, so the caller is effectively
   choosing a whole template over those spans. Is that the intended resolution for C4, or should
   long blocks be sub-divided at a declared coordinate? This changes the legal switch set and so is
   a model change needing its own gate.
3. **`dropped_multiblock_markers` — answered, confirmation requested.** 114 of 47,960 markers are
   dropped for varying in more than one block, and no block loses its discriminating markers. The
   C4A/C4B duplication costs the proposal model nothing material, so the current semantics should
   stand. Confirm rather than decide.
4. **Divergence trim.** The adapter canonicalises block boundaries by trimming each variant block to
   the span over which its alleles differ. This is a block-boundary definition, so by the contract's
   own rule it is a model choice and not merely an implementation detail. It is what makes the
   node-split fingerprint gate pass. Confirm it, and confirm that leaving blocks with a bypass allele
   untrimmed is acceptable for version one.
5. **Oracle provenance.** Confirm that re-deriving the Phase 6 oracle independently, rather than
   porting the archived exhaustive scorer, is the intended reading of "the archived exhaustive
   substitution-only fragment model is the oracle". It is now certified against hand arithmetic and
   in turn certifies the accelerated path on twelve random seeds.
6. **The declared-outlier policy.** This is the one place where a model statement was added to make
   certification reachable, and it is the question most worth a reviewer's attention. A fragment
   that places on no candidate at the working band is declared an outlier, closing a 244-nat
   interval to 0.52. The alternative is that every real run returns `RESCORE_INCOMPLETE`, because
   the unconditional bound cannot be closed by any band the pigeonhole index supports. Confirm the
   policy, its default, and that reporting `fragments_declared_outlier` is sufficient auditing.
7. **Candidate-set narrowness.** The beam returns 8 unique candidates at K=64 on the C4 dry run.
   Whether that is the switch prior at 0.01, the marker emissions, or a real property of C4 is not
   yet established, and it bears directly on the truth-in-candidate rate the go/no-go will measure.
