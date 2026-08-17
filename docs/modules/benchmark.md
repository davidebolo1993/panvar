# Module `benchmark`

CLI: `panvar benchmark`

## What it does

Answers two different questions about the caller's own output, and keeps them apart.

**What was there to be found, and what became of it** — the truth event ledger. Independently of any call, the haplotype's walk over a bubble is node-aligned to the reference walk, and each maximal divergent block between shared anchors is one truth event, sized `max(ref_bp, hap_bp)` from the walks. Every event is then classified:

| class | definition |
|-------|------------|
| `called` | size `≥ --min-sv-bp` and a specific emitted record shares at least one node with the block |
| `missed` | size `≥ --min-sv-bp` and no emitted record covers it |
| `below_threshold` | size `< --min-sv-bp` — variation the caller never set out to emit |

Event size is a property of the two walks. It does not move with the alignment, with the calls, or with anything else the run does, which is what makes it usable as truth.

`called` is **attribution, not coverage**: a record shares at least one node with the block. It does not establish that the record spans the block, nor that it represents it correctly. So `missed` is a firm negative — nothing the caller emitted touches this event at all — while `called` is an upper bound on discovery. It is also silent on whether *this haplotype* was genotyped as carrying it; that is the carrier table below.

**How much of the sequence comes back** — four reconstructions over the same bubbles, against the same truth, on the same QV scale, differing only in what each is allowed to use.

| level | reconstruction built from | answers |
|-------|---------------------------|---------|
| `graph` | reference walk + the haplotype's own steps at every block sharing a node with **any** call at that bubble | optimistic upper bound: can the graph hold this haplotype, and did the caller flag the divergent blocks |
| `called` | the same substitution restricted to blocks a **specific** record is attributed to and that reach `--min-sv-bp` | the ceiling the retained records would reach if each reproduced its block exactly |
| `carrier` | `called` **plus** the requirement that this haplotype's `GT` names one of the records overlapping the block (`--vcf`) | the same ceiling, after carrier assignment is taken into account |
| `genotype` | reference sequence + only the edits this haplotype's `GT` names, from the VCF alone (`--vcf`) | what a consumer reconstructing this sample from the VCF actually gets |

`graph` is deliberately optimistic and is not a genotyping score: sharing a node with some call is not the same as matching one, so a called block can authorise copying a neighbouring uncalled allele, and no genotype is read at all. `called` narrows that to per-record attribution and the size threshold, but it still splices in the haplotype's **true** block rather than the record's `REF`/`ALT` effect — so it too is a ceiling, just a much tighter one. Only `genotype` reconstructs what the records actually say.

`graph ≥ called ≥ carrier` are **nested ceilings**: each substitutes a subset of the previous level's blocks, and all three implant the haplotype's own true sequence, so the ordering is guaranteed. **`genotype` is not bounded by them.** It applies every record the haplotype carries rather than only attributed eligible blocks, so it can beat `carrier` in places — which is exactly why the last two loss terms are signed. Treat it as a fourth measurement, not the bottom of a chain.

What the steps separate: `called` and `carrier` implant the **same true blocks** and differ only in whether the haplotype was genotyped as a carrier — that step is a clean isolation. The `carrier`-to-`genotype` step is not symmetric: `carrier` substitutes only blocks that are attributed *and* eligible, while `genotype` applies **every** record the haplotype's `GT` names, including ones covering no eligible truth block. So that step mixes two things, and they are reported separately — `representation` where the haplotype has an eligible truth event, `false_positive_damage` where it does not.

That is emitted as a `loss_bp` partition of the genotype residual — five consecutive steps from truth down to the reconstruction, which **sum exactly** to it (`sum_check`, asserted, because an attribution built on arithmetic that has stopped closing is worthless):

| term | bases lost to |
|------|---------------|
| `out_of_scope` | truth events **below** `--min-sv-bp`. `call` was never asked to emit these, so they are not a caller failure |
| `discovery_or_attribution` | eligible truth events no retained record covers |
| `carrier_missed` | covered by a record, but this haplotype is not genotyped as carrying it — a false negative |
| `representation` | right record, right carrier: what the VCF's encoding of it fails to reproduce |
| `false_positive_damage` | the VCF edits a haplotype that has no eligible truth event there |

**Every comparative figure is computed over the common set** — the `(haplotype, bubble)` observations *all* levels could score, reported as `common_set/observations`. `graph` and `called` cover every graph haplotype, `carrier` only those joined to a VCF column, and `genotype` only those whose bubble could also be placed; totalling each over its own population and subtracting compares different denominators. With one haplotype left out of the VCF that produced `carrier_walk` 100% sitting *above* `called` 83.3%, and a carrier loss of **minus** 20 bp.

The last two terms are **signed**. `carrier` only substitutes blocks that are eligible *and* attributed, while `genotype` applies every record the haplotype carries, so a record can improve a region `carrier` left alone. A negative term means the VCF beat the walk-based ceiling there; clamping it to zero broke the partition on four of six loci.

Measured on the six reference loci, this is lopsided: out-of-scope 13–78 kb, discovery **0 bp everywhere**, missed carriers 0 at four loci and at most 15.6 kb, false-positive damage small (and negative at C4), and representation 43 kb to 8.57 Mb (the maximum is LPA; GSTM1 is next at 2.7 Mb). The caller finds the eligible variation and puts it on the right haplotypes; the VCF's encoding of it is where the sequence is lost.

A carrier at a block is judged against **every** record overlapping it, not only the one the ledger attributes it to. One divergent region can be described by several records and different carriers take different ones, so asking about the single attributed record would really be asking whether the haplotype carries whichever record sorted first.

A do-nothing **baseline** — plain reference, no edits — is always computed alongside, and is the denominator of the headline:

```
gap_closed          = (baseline_delta - genotype_delta) / (baseline_delta - graph_delta)
variation_recovered = (baseline_delta - <level>_delta)  /  baseline_delta
```

`gap_closed` measures against the *achievable* bound; `variation_recovered` against all of the variation there was. Reported for every level, the pair says which of representation or graph is the limit. **When `baseline_delta == graph_delta` there is no distance to close and `gap_closed` is `NA`, not `1.0`** — the two are equal precisely when a bubble was missed entirely, so reporting `1.0` scored a total miss as a perfect result.

The baseline also makes the metric self-checking: a haplotype genotyped as carrying nothing must score byte-identically to it, so a coordinate or orientation error in the projection shows up immediately instead of hiding inside a plausible QV.

### Carrier calls

Per (haplotype, bubble): does the VCF genotype this haplotype as carrying an alt, against whether its walk really diverges by at least `--min-sv-bp`. Truth comes from the ledger, so it is a property of the walks. Reported for `ALL` only — a bubble-level judgement has no single SV type, and fanning it out over every type present in the bubble reported one observation several times. The per-type partition is the ledger, where each event maps to at most one record.

`FN` here means the record exists but this haplotype's `GT` is `0` at every record of the bubble: a **genotyping** miss, distinct from a `missed` event, which is a **discovery** miss.

### On the residual split

`resid_run_lt_bp` / `resid_run_ge_bp` divide the residual by contiguous alignment-run length. **This is a property of the edit path edlib returns, not of any event.** A clean 60 bp deletion of ordinary non-repetitive sequence comes back as fourteen runs of 1–10 bases, because co-optimal paths distribute the gap over chance matches; every real locus therefore reported `over_threshold_bp = 0`. The columns are kept as a description of the residual's shape. The called/missed question belongs to the ledger.

### On `--min-sv-bp`

It sets the event-size threshold for the ledger and for the `called` reconstruction. It does **not** re-run discovery, so lowering it here does not find anything new — it reclassifies. A real threshold experiment runs `call` at each threshold and `benchmark` on that run's output at the same threshold, over a fixed bubble set: `scripts/validate_benchmark_threshold.py` does exactly that and checks the invariants that must hold (the truth events themselves cannot move, the classes must partition, `below_threshold` cannot rise as the threshold falls, the allele VCF must stay exact).

Pass it the same `--call-extra` the pipeline uses, or it measures a different caller. Measured at C4 and CYP2D6 with `--cn`: every newly eligible event is attributed to a record at every step, but genotype reconstruction is **not** monotone — at C4 it falls 67.4% → 41.7% between thresholds 50 and 10, with the count of haplotypes made worse than plain reference rising 7 → 43. More records mean more merged and overlapping records describing one walk, and independent records do not compose. The script reports those separately from the hard invariants rather than failing on them, because they are call-side interactions and are the point of running the sweep.

Algorithm and worked trace: [algorithms/benchmark.md](../algorithms/benchmark.md).

## Required inputs

- `-i, --gfa <graph.gfa>` — the same passed graph the calls were made on. Validated: every path step must name a node the graph has, consecutive steps must be joined by a link, path names must be unique and link overlaps must be zero, because every walk is spelled by concatenating whole segments.
- one of `-b, --bubble-prefix-in <prefix>` or `-c, --bubbles-csv-in <path>`. Bubble ids must be unique and every node named must exist in the graph.
- `--variant-nodes <path>` — `call`'s `<prefix>.variant_nodes.tsv`. Strictly validated, because a call's node set is the whole basis on which a truth event is attributed to it: a stale node produces exactly the output a genuine caller miss produces. The header must be `call`'s, every row must have the header's field count, `variant_id` must be non-empty and unique, and **every node must belong to the bubble the row names**. Any bubble id absent from the CSV is refused, not just the case where all of them are.
- `-r, --reference-path <name>` — the diff baseline, resolved through the shared resolver (exact name, else a unique case-insensitive substring; ambiguity is an error, never settled by file order).
- `-o, --out-prefix <prefix>`.

Optional: `--vcf <path>` — either VCF `call` produces. **The two are different contracts and the mode is detected explicitly**, reported as `vcf_mode`. A *region* VCF has one row per call and its IDs must match `variant_nodes.tsv` completely, agreeing on `BUBBLE_ID`; anything less is two different runs and is refused. An *allele* VCF has one row per bubble (`bubbleN_ALLELES`) and shares no ID with a per-call row, so the `carrier` level and the per-call `loss_bp` terms **do not exist** there and are reported `NA` / `not_applicable` rather than as zeros — reporting them anyway turned a perfect 0 bp reconstruction into a 1.4 Mb missed-carrier loss with negative representation. In allele mode the meaningful outputs are the genotype reconstruction and `loss_bp/serialization_residual`.

 `<prefix>.region.vcf` scores the merged, interpreted output; `<prefix>.alleles.vcf` (from `call --allele-vcf`) scores the lossless one, which is a serialization ceiling rather than a call-sensitivity score. Sample columns are joined to graph paths by exact name; each column is one haplotype, so genotypes are haploid and nothing has to be phased. Also validated: unique record IDs, a `BUBBLE_ID` on every record (without it a record cannot be placed and would silently pile up against bubble 0), a field count matching the `#CHROM` header, and haploid `GT`. Duplicate sample columns are refused. A diploid `GT` is refused rather than parsed — one column is one haplotype, and `0/1` read as an integer becomes `0`, silently scoring a heterozygous carrier as reference. A partial join is accepted and **prominently reported**, because a QV over a subset is not the QV of the run.

All outputs are staged and committed only on success, and no output may name an input.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--vcf <path>` | add the genotype level and the baseline it is measured against | off |
| `--called-bubbles-only` | score only bubbles carrying `≥ 1` call. Ascertainment-biased by construction: the bubble the caller said nothing about is exactly where a miss lives | off — every reference-traversed bubble is scored |
| `--min-sv-bp <N>` | event-size threshold for the ledger and the `called` reconstruction; set it to what `call` ran with | `50` |
| `--dup-model <cn\|cnbp>` | lay down `CN` copies of `RU_LEN` in place of `REF_CN`, or apply the per-sample `CNBP` delta. Both tile an inferred reference span, so a DUP is reconstructed at the right *length* out of approximately right sequence and is counted `heuristic` | `cnbp` |
| `--no-truth-events` | skip the per-event ledger table | off |
| `--vcf` (again) | also enables the `carrier` level and the `loss_bp` partition, both of which need genotypes | — |
| `--threads <N>` | worker threads over haplotypes | `0` (auto) |
| `-q, --quiet` | disable progress logs | — |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.truth_events.tsv` | one row per (haplotype, bubble, truth event): `ref_bp`, `hap_bp`, `size_bp`, `class`, and the `variant_id` / `svtype` it maps to. The ledger everything else aggregates |
| `<prefix>.qv.tsv` | per (haplotype, bubble): `svtypes`, `is_carrier`, `δ`, `aln_len`, `qv`, the run split `resid_run_lt_bp` / `resid_run_ge_bp`, the ledger counts `truth_events` / `truth_called` / `truth_missed` / `truth_below` / `truth_missed_bp`, `called_delta` / `called_qv`, and `carrier_delta` / `carrier_qv`. With `--vcf`, also `gt_delta` / `gt_aln_len` / `gt_qv`, the baseline `ref_delta`, and `gt_called_carrier` / `gt_true_carrier` |
| `<prefix>.qv_by_haplotype.tsv` | the same per haplotype: `identity`, the run split, the ledger, `called_sum_delta` / `called_qv` / `called_quintile`, `carrier_sum_delta` / `carrier_qv`, plus `qv`, cosigt `band`, `qv_max`, `qv_ratio`, `quintile`. With `--vcf`, the `gt_*` mirror plus `ref_sum_delta`, `ref_qv` and `gap_closed` (literal `NA` where undefined) |
| `<prefix>.qv_summary.tsv` | `truth_event` / `truth_bp` (the ledger, `ALL` and per svtype), `quintile` / `haplotype` / `called_quintile` distributions, `residual_run`, `excluded` (bubbles with no reference walk, haplotype-bubble pairs not traversed, pairs decomposed coarsely, haplotypes with no VCF column, VCF samples with no path), `called_recon`. With `--vcf`: `gt_quintile` / `gt_haplotype`, `gt_gap` (pooled deltas, `gap_closed_pooled`, `gap_closed_mean`, `gap_closed_undefined`, `worse_than_baseline`), `variation_recovered` for every level, `carrier_quintile` / `carrier_recon` and the `loss_bp` partition, `gt_carrier` (TP/FP/FN/TN with precision and recall) and `gt_records` (`applied` / `unplaceable` / `clamped` / `unhandled` / `ref_mismatch` / `heuristic`) |

Nothing is dropped silently: every haplotype-bubble pair excluded from a denominator is counted in the `excluded` scope, so no rate is read as though it covered everything.

Results can be plotted via `scripts/plot_benchmark.R`.
