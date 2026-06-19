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

- Lees JA, Galardini M, Bentley SD, Weiser JN, Corander J. **pyseer: a comprehensive tool for microbial
  pangenome-wide association studies.** *Bioinformatics* 34(24):4310–4312, 2018. — the fsm-lite
  `<feature> | strain:count` table `describe` emits is the format such presence/absence-based tools read.
- Rahman A, Hallgrímsdóttir I, Eisen M, Pachter L. **Association mapping from sequencing reads using
  k-mers.** *eLife* 7:e32920, 2018. — the count/abundance idea behind testing marker *multiplicity*
  (copy number), which presence/absence GWAS cannot see at a fixed-copy repeat. panvar's
  `scripts/gwas_demo.py` is a transparent stand-in that runs both tests on the same file.

## Lp(a) / KIV-2 biology (the worked example)

- Schmidt K, Noureen A, Kronenberg F, Utermann G. **Structure, function, and genetics of
  lipoprotein(a).** *Journal of Lipid Research* 57(8):1339–1359, 2016. doi:10.1194/jlr.R067314.
- Coassin S, Kronenberg F. **Lipoprotein(a) beyond the kringle IV repeat polymorphism: The complexity
  of genetic variation in the LPA gene.** *Atherosclerosis* 349:17–35, 2022.
  — the inverse relationship between KIV-2 copy number and plasma Lp(a) used in `gwas_example.md`.
