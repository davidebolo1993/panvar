# References

Tools `panvar` builds on (GitHub) and the papers behind the ideas in each module.  The module pages ([modules/](modules/)) and the algorithm notes ([algorithms/](algorithms/)) link here.

## General on pangenome graphs

- [pangenome/pggb](https://github.com/pangenome/pggb) — the pangenome graph builder whose `*.smooth.final.gfa` is the usual input to `bubble`.
- [pangenome/odgi](https://github.com/pangenome/odgi) — graph sort / visualization; `panvar`'s internal reference sort mirrors the odgi-style ordering used downstream.
- [vgteam/vg](https://github.com/vgteam/vg) — `vg snarls` is the reference behaviour `bubble`'s internal snarl finder reproduces.
- [pangenome/PanSN-spec](https://github.com/pangenome/PanSN-spec) — the `sample#hap#contig:start-end` path naming `--gtf` annotation relies on.

## `bubble`

GitHub: [vgteam/vg](https://github.com/vgteam/vg) (snarls), [BubbleGun](https://github.com/fawaz-dabbaghieh/bubble_gun) (superbubble).

- Paten B, Eizenga JM, Rosen YM, Novak AM, Garrison E, Hickey G. **Superbubbles, Ultrabubbles, and Cacti.** *J. Comput. Biol.* 25(7):649–663, 2018. <https://doi.org/10.1089/cmb.2017.0251> — the snarl / cactus decomposition `bubble` uses.
- Onodera T, Sadakane K, Shibuya T. **Detecting Superbubbles in Assembly Graphs.** WABI 2013. — the superbubble definition .
- Dabbaghie F, Ebler J, Marschall T. **BubbleGun: enumerating bubbles and superbubbles in genome graphs.** *Bioinformatics* 38(17):4217–4219, 2022. <https://doi.org/10.1093/bioinformatics/btac448>

## `panphorte`

GitHub: [GenoGra/Panphorte](https://github.com/GenoGra/Panphorte) — the original tool whose tandem-repeat normalization idea `panphorte` re-implements.

- Coggi M, Basile L, Branchini B, Amodeo G, Di Donato GW, Santambrogio MD. **On the optimization of copy number variations representation in pangenome graphs.** *Front. Bioinform.* 2026;6:1811916. <https://doi.org/10.3389/fbinf.2026.1811916> — the original panphorte tandem-repeat normalization.

## `call`

GitHub: [lh3/minimap2](https://github.com/lh3/minimap2) — used for INS subtype realignment and, in the per-gene copy-number resolver, the one-time reference alignment of a near-identical paralog pair's coding sequences (to locate their divergent columns for the per-site split). Per-haplotype copy number is read by k-mer dosage, without alignment.

- Li H. **Minimap2: pairwise alignment for nucleotide sequences.** *Bioinformatics* 34(18):3094–3100, 2018. <https://doi.org/10.1093/bioinformatics/bty191>.

## `describe`

The closed-syncmer sampling and dosage-based (multiplicity-aware) encoding behind `describe`'s BIMBAM substrates — carrying copy number through, so a fixed-copy repeat's dosage stays testable where a presence/absence k-mer encoding would miss it.

- Edgar R. **Syncmers are more sensitive than minimizers for selecting conserved k-mers in biological sequences.** *PeerJ* 9:e10805, 2021. <https://doi.org/10.7717/peerj.10805> — the closed-syncmer sampling `describe` uses by default.
- Rahman A, Hallgrímsdóttir I, Eisen M, Pachter L. **Association mapping from sequencing reads using k-mers.** *eLife* 7:e32920, 2018. <https://doi.org/10.7554/eLife.32920> — testing marker multiplicity (copy number).

## `inspect`

GitHub: [marbl/Mash](https://github.com/marbl/Mash) — the bottom-k MinHash sketch used to estimate walk Jaccard.

- Broder AZ. **On the resemblance and containment of documents.** SEQUENCES 1997. — MinHash, the sketch behind clustering.
- Ondov BD, Treangen TJ, Melsted P, et al. **Mash: fast genome and metagenome distance estimation using MinHash.** *Genome Biology* 17:132, 2016. <https://doi.org/10.1186/s13059-016-0997-x> — bottom-k sketch Jaccard estimation.

## `associate`

GitHub: [genetics-statistics/GEMMA](https://github.com/genetics-statistics/GEMMA) (reference LMM GWAS; `describe`'s BIMBAM is its native format), [libeigen/eigen](https://gitlab.com/libeigen/eigen) (the linear algebra behind the LMM eigendecomposition).

- Zhou X, Stephens M. **Genome-wide efficient mixed-model analysis for association studies.** *Nature Genetics* 44:821–824, 2012. <https://doi.org/10.1038/ng.2310> — GEMMA tool.
- Kang HM, Sul JH, Service SK, et al. **Variance component model to account for sample structure in genome-wide association studies (EMMAX).** *Nature Genetics* 42:348–354, 2010. <https://doi.org/10.1038/ng.548> — the eigendecompose-once, test-each-marker fast-LMM approximation.
- Benjamini Y, Hochberg Y. **Controlling the false discovery rate.** *J. R. Stat. Soc. B* 57(1):289–300, 1995. — the BH FDR.

## Biology (LPA)

- Schmidt K, Noureen A, Kronenberg F, Utermann G. **Structure, function, and genetics of lipoprotein(a).** *J. Lipid Res.* 57(8):1339–1359, 2016. <https://doi.org/10.1194/jlr.R067314>.
- Coassin S, Kronenberg F. **Lipoprotein(a) beyond the kringle IV repeat polymorphism: the complexity of genetic variation in the LPA gene.** *Atherosclerosis* 349:17–35, 2022. — the inverse KIV-2 / Lp(a) relationship
- **Moli-sani** study (Italian / Southern-European cohort): Lipoprotein(a) as an early marker of ardiovascular events. *Frontiers in Cardiovascular Medicine* 12:1571395, 2025. — Southern-European cohort context for the simulated covariate/Lp(a) ranges
