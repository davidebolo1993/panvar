# References

Background literature for the concepts panvar builds on. panvar is self-contained at runtime (no
external graph or alignment tools); these are cited for the *ideas*, not as runtime dependencies.

## Graph / bubble structure

- Paten B, Eizenga JM, Rosen YM, Novak AM, Garrison E, Hickey G. **Superbubbles, Ultrabubbles, and
  Cacti.** *Journal of Computational Biology* 25(7):649–663, 2018. doi:10.1089/cmb.2017.0251.
  — the snarl / ultrabubble decomposition `bubble` reproduces internally.
- Garrison E, Sirén J, Novak AM, et al. **Variation graph toolkit improves read mapping by representing
  genetic variation in the reference.** *Nature Biotechnology* 36:875–879, 2018. doi:10.1038/nbt.4227.
  — variation graphs and the `vg snarls` JSONL format `bubble` can ingest.
- Guarracino A, Heumos S, Nahnsen S, Prins P, Garrison E. **ODGI: understanding pangenome graphs.**
  *Bioinformatics* 38(13):3319–3326, 2022. doi:10.1093/bioinformatics/btac308.
  — the sort/flip conventions panphorte mirrors when emitting a call-ready graph.
- Garrison E, Guarracino A, Heumos S, et al. **Building pangenome graphs (PGGB).** *Nature Methods*
  2024. — how the test graphs (`tests/real_data/*.gfa.gz`) were built.

## k-mer features

- Edgar RC. **Syncmers are more sensitive than minimizers for selecting conserved k-mers in biological
  sequences.** *PeerJ* 9:e10805, 2021. doi:10.7717/peerj.10805. — closed-syncmer sampling used by
  `describe` (`--feature-mode syncmer`).

## Copy-number vs presence/absence association

- Rahman A, Hallgrímsdóttir I, Eisen M, Pachter L. **Association mapping from sequencing reads using
  k-mers.** *eLife* 7:e32920, 2018. — the count/abundance idea behind testing marker *multiplicity*
  (copy number), which presence/absence GWAS cannot see at a fixed-copy repeat. panvar's
  `associate` module tests this dosage directly (see [gwas_example.md](gwas_example.md)).
- Lees JA, Galardini M, Bentley SD, Weiser JN, Corander J. **pyseer: a comprehensive tool for microbial
  pangenome-wide association studies.** *Bioinformatics* 34(24):4310–4312, 2018. — background on
  pangenome-wide association; its k-mer mode is presence/absence, which is *why* panvar carries dosage
  through to BIMBAM and tests multiplicity directly rather than emitting a presence/absence k-mer table.
- Zhou X, Stephens M. **Genome-wide efficient mixed-model analysis for association studies (GEMMA).**
  *Nature Genetics* 44:821–824, 2012. — the BIMBAM mean-genotype dosage format `describe` exports and the
  LMM that `associate --model lmm` mirrors.

## Lp(a) / KIV-2 biology (the worked example)

- Schmidt K, Noureen A, Kronenberg F, Utermann G. **Structure, function, and genetics of
  lipoprotein(a).** *Journal of Lipid Research* 57(8):1339–1359, 2016. doi:10.1194/jlr.R067314.
- Coassin S, Kronenberg F. **Lipoprotein(a) beyond the kringle IV repeat polymorphism: The complexity
  of genetic variation in the LPA gene.** *Atherosclerosis* 349:17–35, 2022.
  — the inverse KIV-2 ↔ Lp(a) relationship, the log-normal ~0.3–300 mg/dL (median ~10–12) distribution, and
  the small-isoform (≤22 KIV) ~5× effect used to ground `make_lpa_phenotype.py`.
- **Moli-sani** study (Italian / Southern-European cohort): Lipoprotein(a) as an early marker of
  cardiovascular events. *Frontiers in Cardiovascular Medicine* 12:1571395, 2025. — a Southern-European
  cohort context for the simulated covariate/Lp(a) ranges; >50 mg/dL marks clinical high risk.

## Population structure / mixed models (the `associate` LMM)

- Kang HM, Sul JH, Service SK, et al. **Variance component model to account for sample structure in
  genome-wide association studies (EMMAX).** *Nature Genetics* 42:348–354, 2010. — the eigendecompose-once,
  test-each-marker fast-LMM approximation `panvar associate --model lmm` implements.
- Price AL, Patterson NJ, Plenge RM, et al. **Principal components analysis corrects for stratification in
  genome-wide association studies.** *Nature Genetics* 38:904–909, 2006. — the PC-covariate approach behind
  `--pca`.
