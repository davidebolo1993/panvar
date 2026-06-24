# SMBE 2026 — Talk speech draft
**Genotyping, SV Calling and Association Testing using Pangenome Graphs**
Davide Bolognini — 12 min talk + 3 min Q&A

> Pacing note: ~1,550 words ≈ 12 min at a calm technical rate. This is content-dense
> for the slot — trim markers `[CUT IF LONG]` show where to drop ~60–90 s if you run over.
> Numbers in brackets are slide indices in `SMBE2026_Bolognini.pdf`.

---

## 0. Opening + the merge — slides [1]–[4] · ~0:45

[1] Good morning, I'm Davide Bolognini, from Human Technopole and GenoGra.

[2] Before I start — a quick plug. If you've not yet seen Chiara Paleni's poster on
haplotype-based structural-variant detection, you still have time: half of what I'll show
in the second part of this talk is her work.

[3] Because, originally, I was asked to talk only about *genotyping*. Then we were asked
to merge the two stories — so this is now a bit of a fusion.

[4] Which is why the title changed: **genotyping, SV calling, and association testing,
all on pangenome graphs.** So: I'll start with our genotyper, COSIGT, show a couple of
benchmarks — including one on ancient DNA that's relevant for this audience — a population
application on HLA, then our amylase work, and finally a new tool, panvar, that takes us
from a graph all the way to an association test.

## 1. COSIGT — the method — slides [5]–[6] · ~1:30

[5] The problem COSIGT solves: at complex, repetitive loci, linear-reference genotyping
breaks down. We instead work *locally* and *graph-natively*.

Walk through it left to right. For a target locus we take the high-quality assembled
haplotypes — say from the HPRC — and build a **local pangenome graph** (1). Each haplotype
is a path, and its node-coverage vector is just how many times it traverses each node.
Then we take a short-read sample, align it to *all* haplotypes at once (2), and get the
read coverage over the same nodes (3b). Now the trick (4): for every possible **pair** of
haplotypes — every candidate diploid genotype — we add their two coverage vectors to get
the coverage that genotype *would* produce. We score each candidate by **cosine similarity**
against the observed read coverage (5–6), and the best-matching pair is the call. So
genotyping becomes: which combination of known haplotypes best explains what I sequenced.

[6] [CUT IF LONG] The payoff is that this is reference-free *within* the locus — it
handles paralogs, copy-number and structural diversity that a single linear reference
simply can't represent. This is joint work with Andrea Guarracino and the team.

## 2. A worked locus + benchmark design — slides [7]–[9] · ~1:30

[7]–[8] Here's a real region — the CYP2D6/CYP2D7 pharmacogene cluster. Each row is a
haplotype, and you can see the coverage structure tracks the gene annotation: the duplicated,
recombining segments that make this locus notoriously hard are exactly where the haplotypes
diverge. This is the kind of place COSIGT is built for.

[9] To benchmark properly we used two designs. In **leave-zero-out** the sample's true
haplotypes *are* present in the graph — this measures raw genotyping accuracy when the right
answer is available, scored by an alignment-based quality value, QV. In **leave-all-out** we
do the opposite: we **exclude** the sample's true haplotypes entirely, so COSIGT must fall
back on the most similar *available* haplotype — this is the honest test of what happens with
novel or absent variation, and we score it as the fraction of the maximum achievable quality.
We ran this across three region sets — challenging medically-relevant genes, pharmaco- and
immunogenes, and a set of structural variants — under full coverage, downsampling to 1, 2
and 5×, and a simulated ancient-DNA condition.

## 3. Benchmark results — slides [10]–[11] · ~1:15

[10] Headline result, leave-zero-out: across the medically-relevant genes and the SV set,
COSIGT (green) matches or beats locityper (the leading competitor), and crucially it **holds
up as coverage drops** — at 1× over 93% of calls reach mid-or-higher quality versus ~85% for
locityper, and the edge holds at 2×. That low-coverage robustness is the whole point: it's
what makes this *population-scalable*.

[11] [CUT IF LONG] And leave-all-out — the harder case, where the sample's true haplotype is
*removed* — COSIGT still recovers roughly 88–91% of the best achievable quality across genes
and SVs. So even for variation that isn't in the panel, it lands on a close neighbour rather
than failing.

## 4. Ancient DNA — the SMBE-relevant benchmark — slides [12]–[15] · ~1:30

[12] Now the part most relevant for this room. Ancient DNA is the hardest case: very low
coverage, short fragments, deamination damage, and modern-human contamination. We simulated
exactly that with a purpose-built simulator, ancestralsim, wrapping Gargammel — single-end
reads at 1× and 2×, ~70 bp fragments, single-stranded deamination, and 0% or 10%
contamination — on 48 genes enriched for pharmacogenetic and immunogenetic relevance: the
HLA, CYP and MUC families.

[13] Per gene — across CYP, HLA and MUC families — COSIGT stays mostly in the high-QV bins
even at 1× with contamination. This is where the gap really opens: at 1× COSIGT keeps about
94% of calls at mid-or-higher quality while locityper collapses to roughly 46%; at 2× it's
~95% versus ~60%. And it's barely worse than matched *modern* DNA at the same depth — so
cosine similarity's insensitivity to coverage magnitude is exactly what carries it through
degraded, low-depth ancient samples.

[14]–[15] [CUT IF LONG] One more worry for this audience specifically: an ancient sample
isn't just *damaged*, it's *evolutionarily diverged* from the modern haplotypes in our
pangenome — so is divergence itself a confounder? We tested that directly with a
split-demography simulation — using msprime to add ancestral substitutions onto the
haplotypes under a modern/ancient population split — and COSIGT holds up even then. So the
accuracy isn't an artefact of simulating only damage; genuine evolutionary distance doesn't
break it. The message: pangenome genotyping is viable on ancient samples, which opens these
complex loci to evolutionary questions through time.

## 5. Population application — HLA in MoliSANI — slide [16] · ~1:00

[16] To show it works at population scale on real data: we genotyped classical HLA genes in
~1,085 whole-genome samples from the Moli-sani cohort — a large Italian population study —
and compared COSIGT (green) against t1k, a dedicated HLA tool (purple). The bars are the
percentage of calls matching imputed reference types — so the 100% line is perfect
concordance. Mean haplotype-level accuracy is 89.5% for COSIGT versus 90.8% for t1k —
essentially on par across HLA-A, B, C, DP and DQB1, with a *general-purpose* graph method,
not an HLA-specific one. DQA1 is hard for
both tools, so that's a shared limitation, not ours. The point: one framework, competitive
with specialist callers, at cohort scale.

## 6. Amylase — evolution over time — slides [17]–[18] · ~1:15

[17] The amylase locus is a beautiful test case: highly copy-number-variable, structurally
diverse, and tied to diet. From long-read assemblies we resolved **28 distinct structural
haplotypes** — differing not just in copy number but in the arrangement and orientation of
the duplicated units. And the striking result from the phylogeny is **recurrence**: nearly
identical structures have arisen *independently*, again and again, on different haplotype
backgrounds — driven by non-allelic homologous recombination between the long paralogous
segments, at a structural turnover rate orders of magnitude above the SNP rate. The ancestral
haplotype that seeded all of this traces back to roughly **280,000 years ago**, before the
out-of-Africa expansion. So this is an old, restless locus that keeps reinventing the same
architectures.

[18] And because the genotyping works on ancient samples, we could put this on a **time
axis** across West-Eurasian ancient genomes. Two things stand out. First, the haplotype
frequencies turn over: the high-copy, duplication-containing haplotypes rise more than
**sevenfold — from about 12% to about 86%** over the last twelve thousand years — while the
ancestral haplotypes collapse from roughly 88% to 14%. That's the H1ᵃ-down, H3ʳ-up pattern on
the plot. Second, that shows up directly as copy number: AMY1, AMY2A and AMY2B all increase
significantly through time — the p-values on the slide, around 10⁻⁶ for AMY1 and AMY2A and
0.003 for AMY2B — with AMY1 gaining on average about three copies. The onset, around nine
thousand years ago, and a selection coefficient near 0.02 point to **positive, directional
selection following the spread of agriculture**. So we read diet-driven selection straight
off ancient genomes — that's what this whole approach unlocks: not just "what's the genotype,"
but "how did it change through time." [Bolognini et al., 2024]

## 7. panvar — from graph to variants — slides [19]–[24] · ~1:45

[19] So far: genotyping. Now the second half — calling and testing variants directly on the
graph. Our driving example is **LPA**, which carries the **KIV-2** tandem repeat — its copy
number is one of the strongest genetic determinants of Lp(a) and cardiovascular risk, and
it's essentially invisible to standard short-read calling.

[20] In the graph the KIV-2 repeat shows up as this characteristic loop — a self-cycling
structure. That topology *is* the copy number.

[21]–[24] panvar reads the graph the way COSIGT does — per-path node coverage — and
decomposes each variable site into a **bubble**. A bubble can be a simple substitution, or a
tandem expansion: here you can literally read the alleles — CG, CGCG, CGCGCG — as growing
repeat copies. panphorte normalizes the repeat into an explicit copy-number-bearing
structure, and in the coverage matrix the KIV-2 bubble lights up as a single column whose
value *is* the number of repeat copies each haplotype carries. [Coggi et al., 2026]

## 8. SV typology + CN validation — slides [25]–[30] · ~1:30

[25]–[27] panvar classifies the full SV vocabulary straight from graph topology:
deletions, insertions — merged when they're close or similar enough — inversions, and
copy-number. The clever part is copy number, which we derive three ways depending on how the
locus is represented: from an explicit panphorte-normalized tandem repeat; from **node
multiplicity** — the peak of a self-loop — when it isn't normalized; or from **coverage**,
when paralogs are folded together, as path-length over unit-length. Three topologies, one
copy-number answer.

[28]–[29] And it's accurate. Against ground-truth copy number across four hard loci — C4,
CYP2D6, GSTM1, and LPA — we get correlations of essentially 1.0 with exact-match rates at or
near 100% for C4, GSTM1 and LPA, including the full KIV-2 gradient up to ~30 copies. CYP2D6
is the hardest — the recombinant locus from earlier — and even there we're at r≈0.84, ~92%.
And where we disagree, panvar tends to call slightly *more* copies than the truth set: those
excess calls are most likely CYP2D6–CYP2D7 hybrid alleles that the ground-truth annotation
didn't resolve — so it's arguably panvar seeing real structure, not making errors.

[30] Calling everything together, each haplotype gets a genotype at every bubble — colored
by event type, with copy number as a shade. The KIV-2 DUP bubble (highlighted) is the one
carrying that continuous copy-number signal, sample by sample. And critically, panvar exports
this as a standard dosage matrix — so it feeds straight into association.

## 9. Association — slide [31] · ~0:45
> **NOTE: figure pending** — region-scan / Manhattan + KIV-2 effect. Cluster is down;
> regenerate `assoc_graph_quant.manhattan.png` (+ the KIV-2 dosage-vs-phenotype panel) and
> drop it in. Speech below describes what the figure shows.

[31] The final step closes the loop: we test those dosages against a phenotype. On our LPA
example, scanning the region, the signal concentrates exactly on the KIV-2 copy-number
bubble — more copies, lower Lp(a), the known biological direction — and it stays the lead hit
after we control for population structure with principal components, which collapses the
inflation back to a well-behaved null. So we go end-to-end — graph, to copy number, to a
calibrated, structure-corrected association — on a variant that conventional pipelines miss
entirely. And because the output is a standard dosage matrix, we reproduce the same hit with
an external tool like GEMMA.

## 10. Wrap-up — slides [32]–[36] · ~0:45

[32] This is very much a team effort — across Human Technopole and PopMed, TGEN, UTHSC, UCB,
UNIPV and Neuromed. Thank you all.

[33] Everything is open: COSIGT, the ancient-DNA pipeline, panphorte and panvar are all on
GitHub.

[34]–[35] [CUT IF LONG] It's all part of making pangenome analyses accessible through
GenoGra — and yes, we're hiring, so come talk to me if pangenomics is your thing.

[36] To sum up: pangenome graphs let us genotype complex loci from low coverage — even
ancient DNA — at population scale, call the structural variation conventional methods miss,
and take it all the way to association. Thank you — I'm happy to take questions.

---

## Anticipated Q&A (3 min)
- **Runtime / scalability?** Local per-locus graphs → embarrassingly parallel; cohort-scale
  on MoliSANI shown.
- **Why cosine similarity vs a likelihood model?** Robust to coverage scale, no error-model
  assumptions; works down to 1×.
- **aDNA damage handling?** Damage/contamination simulated in the benchmark; split strategy
  for short fragments.
- **panvar vs graph SV callers (vg/PanGenie)?** We target copy-number/tandem topology
  explicitly and emit association-ready dosages.
- **Is the LPA association real or simulated?** Be ready to say clearly: the example
  phenotype is **simulated** on a real ~10k cosigt-genotyped cohort to demo the pipeline —
  not a real Lp(a) GWAS finding (yet).
