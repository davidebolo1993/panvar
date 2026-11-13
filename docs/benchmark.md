# Module `benchmark`

CLI: `panvar benchmark`

## What it does

Checks how much of each haplotype's variation `call`'s output accounts for, using the haplotypes the graph was built from as the truth:

1. The variation `call` reported for a haplotype is put into the reference, and the result is aligned to the haplotype. Where `call` reported the variation, the two match; where it did not, the difference is reported.
2. To say how much that is, the plain reference, a rebuild with nothing put in, is aligned the same way: its differences are all there was to find. The share recovered is 1 − the bases left / the bases present.

Each haplotype is rebuilt twice: from `call`'s own events for it (`carrier_events.tsv`, stage `call`), where what is left was missed by `call`; and from the VCF records it carries (stage `vcf`), where what is left is one of those misses again or was lost in writing the VCF. Differences of `--min-sv-bp` or more are listed; smaller ones (SNVs, small indels) are only totalled.

## How it works

Alignments use [edlib](https://github.com/Martinsos/edlib).

1. A haplotype is scored only where it overlaps the reference, so a partial assembly is not charged for what it lacks: between the first and last anchors, nodes both walk once and in the same order. A walk starting or ending on the reference's first or last node is scored to the region's edge, so a deletion right there counts.
2. Each rebuild is the reference edited with the haplotype's events, in reference order:
   - `call`'s events: each replaces a reference interval with the haplotype's own bases. A copy-number carrier's event is its copy stretch, so the events inside it are not applied again: where `call` aligned the allele, its first copy to its last on either walk, cut where the alignment cuts the two walks, so no base is both in the stretch and in another event of the carrier; where it compared the allele node by node, up to the nodes its walk shares with the reference. An insertion at the stretch's edge is inside when counting it brings the stretch's length nearer what `CN` implies.
   - VCF records: every record the haplotype carries, as written: `REF`/`ALT`, or a symbolic record's `INSSEQ`; an `INV` is turned round. A record that cannot be spelled is skipped and counted, and the warning names the cause: an insertion or replacement over `call`'s `--max-seq-bp` (20,000 bp), whose sequence `call` leaves out of the VCF (`call --max-seq-bp 0` writes them all); a copy-number record with no reference copy to paste and no copies named, which no setting writes; or any other.
   - Of two overlapping events the longer is applied and the other counted as a conflict; events past the scored stretch are left out and counted.
3. A `MODULE` record that names its carriers' copies (`CN_FORMS`, `CNFORM`) has each carrier's copies pasted as written: the record and the carrier's other records then rebuild its walk, with nothing to infer. A form `call` wrote without its sequence (over `--max-seq-bp`) is skipped and counted. Otherwise a copy-number record holds a count, not a sequence, so the VCF rebuild pastes or removes `CN − REF_CN` copies of a reference copy, beside the reference's copies, which run from `POS`+1 to `END`. Where copies come in forms (diverged paralogs, long and short copies), a carrier gets the form closest to its own change per copy (`CNBP`). A gained copy shorter than every form is explained by a deletion record inside the copy when the sizes match; when two deletions could explain it, the VCF does not say which, and none is applied.
4. Each rebuild is aligned to the haplotype end to end, by edit distance. To keep long regions fast, the alignment is pinned at exact blocks, runs of 500 bases or more identical in both and found once in each, and only the stretches between them are aligned. A repeat holds no such block, so it is aligned in one piece and its copies go where they fit best. The plain reference is aligned the same way, once.
5. Gaps `--merge-gap` columns apart or fewer form one difference. A stretch dense in mismatches (one column in five, over `--min-sv-bp` columns) is one too, since an inversion or an equal-length replacement aligns almost without gaps. A difference is sized as `call` sizes it, by the larger side.
6. Each difference left is looked up against the haplotype's own events, the records (its own and the others'), and the sites, to say why it is left. An insertion or deletion inside a repeat can be written at any copy, and the alignment picks one end, so it is looked up anywhere it can slide to.

   Left by `call`'s own events (`call`): the first cause that applies.

   | cause | the difference |
   |-------|----------------|
   | `incomplete` | it lies at a record of this haplotype whose event did not rebuild all of it, or could not be placed on the reference |
   | `not_carrier` | a record is there, but this haplotype is not among its carriers |
   | `not_reported` | `call` reported nothing there; `bubbles` says whether it lies in a site or outside every site, where `call` never looked |

   Left by the VCF (`vcf`), each difference whole: one at `call`'s misses, in the same direction, is those misses again, with the cause of the one it overlaps most. One that goes beyond them by `--min-sv-bp` or more was lost in writing the VCF: the first reason that applies.

   | reason | the difference |
   |--------|----------------|
   | `no_sequence` | the record has no sequence to put in |
   | `copy_count` | the copies pasted for a copy-number record differ from the carrier's own |
   | `representative` | the record holds another carrier's version of the event |
   | `other` | none of these, as when the record's breakpoints are approximate (`IMPRECISE`, tried before `representative`) or the alignment places the difference further along a repeat |

   The `vcf` rows with a VCF reason are what only the VCF leaves.
7. Each carrier's own change is also compared with the record it is filed under. The first rebuild uses each carrier's own event, so a carrier filed under a record of the wrong size, or given the wrong copy number, still rebuilds perfectly there; the VCF rebuild shows the loss, but not which carrier and record it comes from. The check names them: the carrier's own events against the record's `SVLEN` (`--min-sv-bp` or more apart is a mismatch), and at a copy-number record its `CN` against the copies its walk holds, and its copy stretch's length against what `CN` implies (`cn_stretch`). Where the record names the carrier's copies, the stretch is instead made from the reference's with those copies and the carrier's own events inside it put in, and aligned to the stretch its walk holds (`cn_exact`): a difference of `--min-sv-bp` or more, or copies and events that overlap, is a mismatch. Differences under `--min-sv-bp` are not, since `call` does not report them; two insertions at one point are tried in both orders, as the VCF does not say which comes first.

## Required inputs

- `-i, --gfa <path>`: the graph `call` ran on.
- `-b, --bubble-prefix <prefix>` or `-c, --bubbles-csv <path>`: its site table, for attributing what is left.
- `-r, --reference-path <name>`: the reference path, as given to `call`.
- `--calls <prefix>`: `call`'s output from one run without `--diploid`, `<prefix>.region.vcf` and `<prefix>.carrier_events.tsv`: each haplotype is scored against its own VCF column, so a VCF with one column per sample is refused.
- `-o, --out-prefix <prefix>`.

## Key options

| flag | what it does | default |
|------|--------------|---------|
| `--min-sv-bp <N>` | smallest difference listed; set below `call`'s to also see what its threshold leaves out | `50` |
| `--merge-gap <N>` | differences N or fewer aligned bases apart count as one | `10` |
| `-t, --threads <N>` | worker threads over haplotypes | `0` (all cores) |
| `-q, --quiet` | no progress logs | off |

## Outputs

| file | contents |
|------|----------|
| `<prefix>.summary.tsv` | totals over the haplotypes scored, as `key`, `value` |
| `<prefix>.haplotypes.tsv` | the same totals, one row per haplotype |
| `<prefix>.events.tsv` | every large difference left, with the rebuild that left it and why |
| `<prefix>.present.tsv` | every large difference between the plain reference and each haplotype: what there was to find |
| `<prefix>.carrier_checks.tsv` | carriers whose own change disagrees with their record |

How each rebuild was made (events applied, those left out for overlapping another or lying past the scored stretch, records that could not be spelled) is in the log, as warnings where it can mean a loss.

`summary.tsv`:

| key | meaning |
|-----|---------|
| `min_sv_bp`, `merge_gap` | the run's settings |
| `haplotypes_scored`, `haplotypes_skipped_<reason>` | haplotypes scored, and those skipped, by reason |
| `large_present`, `large_present_bp` | large differences there were to find: count and bases |
| `large_left_call`, `large_left_call_bp` | large differences `call`'s own events leave |
| `large_left_call_<cause>` | the same, by cause (step 6) |
| `large_left_vcf`, `large_left_vcf_bp` | large differences the VCF's records leave, `call`'s misses again included |
| `large_left_vcf_only`, `large_left_vcf_only_bp` | those with a VCF reason: what only the VCF leaves |
| `large_left_vcf_<reason>` | the same, by reason (step 6) |
| `recovered_by_call_fraction`, `recovered_by_vcf_fraction` | share of the bases present each rebuild accounts for: 1 − left / present; negative when the rebuild is further from the haplotype than the reference |
| `haplotypes_clean_call`, `haplotypes_clean_vcf` | haplotypes with nothing large left |
| `small_bp_present`, `small_bp_left_call`, `small_bp_left_vcf` | bases in smaller indels: present, and left by each rebuild |
| `snv_present`, `snv_left_call`, `snv_left_vcf` | SNVs: present, and left by each rebuild |
| `carriers_checked`, `carrier_mismatches` | carriers checked against their record (step 7), and those that disagree |
| `cn_exact_checked` | carriers whose copies their record names, checked by `cn_exact` (step 7) |

`haplotypes.tsv`, one row per haplotype (`.` for one not scored):

| column | meaning |
|--------|---------|
| `haplotype`, `status` | the haplotype, and `scored` or why not |
| `ref_start`, `ref_end` | the reference stretch scored (step 1) |
| `large_present`, `large_present_bp` | large differences there were to find |
| `large_left_call`, `large_left_call_bp` | those `call`'s own events leave |
| `large_left_vcf`, `large_left_vcf_bp` | those the VCF's records leave |
| `large_left_vcf_only`, `large_left_vcf_only_bp` | those only the VCF leaves |
| `small_bp_present`, `small_bp_left_call`, `small_bp_left_vcf` | bases in smaller indels: present, and left by each rebuild |
| `snv_present`, `snv_left_call`, `snv_left_vcf` | SNVs: present, and left by each rebuild |

`events.tsv`, one row per large difference left, and `present.tsv`, one row per large difference there was to find:

| column | meaning |
|--------|---------|
| `haplotype` | the haplotype |
| `stage` | `events.tsv` only: the rebuild that leaves it, `call` or `vcf` |
| `ref_start`, `ref_end` | where it lies on the reference, from the base before it to the base after |
| `type`, `size_bp`, `del_bp`, `ins_bp` | `DEL`, `INS` or `SUB`; its size, as `call` sizes it, and the bases removed and added |
| `reason` | `events.tsv` only: why it is left, a cause of `call`'s miss or, at `vcf`, a VCF reason (step 6) |
| `bubbles` | the sites it lies in, or can slide to inside a repeat; `.` for none |
| `records` | `events.tsv` only: the records it lies at, `.` for none |

`carrier_checks.tsv`, one row per carrier that disagrees with its record:

| column | meaning |
|--------|---------|
| `variant_id`, `haplotype` | the record and the carrier |
| `kind` | what was compared: the event's size (its kind), the copy number (`cn`), the copy stretch's length (`cn_stretch`), or the stretch made from the carrier's named copies and events (`cn_exact`) |
| `reported`, `own_event` | what the record says, and what the carrier's own change is |
| `difference_bp` | own minus reported (in copies, for `cn`) |

### Visualization

`scripts/plot_benchmark.R` plots one benchmark run or several, each with the `--min-sv-bp` `call` ran with.

- A bar is a share: of the differences present (haplotype against reference), those still left (rebuild against haplotype), stacked by reason.
- Each run has two groups: `call`'s own events (`call`), then the VCF (`vcf`). Above each bar, what is left of what was present. Benchmarked below `call`'s threshold, each group has two bars, all the differences and those of `call`'s threshold or more, named under them; they differ by what `call` was not asked for.
- Four panels: share of the differences and of their bases, over the whole haplotype and inside the sites only (`bubbles`). Options keep one row or column of panels and wrap many runs over several rows; the numbers go to a table beside the figure.
- Over 100%, the rebuild differs from the haplotype more than the reference does, as where pasted reference copies differ from the carrier's own.
- Above 10 runs (`--max-grouped`), `plot_benchmark.R` draws one panel with one bar per run, for one measure (differences, or bases with `--only-bases`), one scope (the whole haplotype, or the sites with `--only-sites`), one rebuild (`--stage call` or `vcf`, default `vcf`) and one size class (`--size all` or `above`, default `above`), 25 runs to a row (`--per-row`).

```
Rscript scripts/plot_benchmark.R --out benchmark.png <bench prefix>[:<call min bp>] [...]
```

`scripts/benchmark_thresholds.sh` shows what a lower `call` threshold recovers. Each example graph's sites are found once at 10 bp, then `call` runs alone on them at 50, 40, 30, 20 and 10 bp, each run benchmarked at 10 bp, so what there is to find is the same for all. Each threshold gets one plot of every graph, laid out as `benchmark.png` (`thresholds.<bp>.png`).

Each script lists its options under `--help`.

## Limitations

- The truth is the panel the graph was built from: how `call` does on other haplotypes is out of scope.
- It measures how well the output rebuilds each haplotype, not whether each call is right: a spurious event leaves a difference that reads as missed.
- A copy-number record that does not name its carriers' copies (a `UNIT` record, or a VCF without `CN_FORMS`) is rebuilt as copies of a reference copy, so how real copies differ shows as `copy_count`. A short copy is inferred from `CNBP` and a deletion record, and the VCF does not say which copy is the short one.
- An imprecise inversion is applied between its `POS` and `END`, so a few kb can be left at its ends.
- A haplotype that starts or ends inside the region is scored between its first and last shared anchors, so a deletion at its very end looks like an assembly that stops there.
- Inside a repeat, a difference with changes besides one insertion or deletion is looked up where the alignment places it, which can be another copy than `call`'s; it then shows as small differences or as `other`.
- Sizes follow `call`'s, except that very unequal replacements, inversions in AT-rich sequence, and deletions of a few units of a diverged repeat can be sized differently.
- The carrier check compares sizes and counts, not positions or sequence: a carrier filed under another record of the same size passes, and differences inside a copy stretch that cancel in length pass too.

## Example

```bash
# what call's output leaves of each haplotype (call run without --diploid)
panvar benchmark -i locus_r.normalized.sorted.gfa -b locus_r -r GRCh38 --calls locus_call -o locus_bench
# differences under call's 50 bp listed too
panvar benchmark -i locus_r.normalized.sorted.gfa -b locus_r -r GRCh38 --calls locus_call -o locus_bench_10 --min-sv-bp 10
```

`scripts/benchmark_examples.sh` benchmarks `call`'s output on every graph of `tests/examples/call` and plots what it leaves, per graph, in `tests/examples/benchmark/benchmark.png`; benchmarked at 10 bp, so the differences under `call`'s 50 bp count too.
