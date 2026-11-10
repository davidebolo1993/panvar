# Real loci

pggb graphs of human loci, one folder per kind of locus. Each folder's loci are checked the same way.

| folder | loci | checked by |
|---|---|---|
| `paralogs/` | c4, gstm, cyp2d, acot | copy number per gene against the pangene tables (`<locus>.bed`; `gstm1.bed` holds GSTM1 copies alone) |
| `tandem_repeat/` | lpa | KIV-2 copies against `lpa.repeats.tsv`; calls and benchmark |
| `inversion/` | ankrd36c | the 335 kb inversion; calls and benchmark |
| `rebuild/` | myom2 | rebuild: the raw graph is too tangled for panphorte |
| `hla/` | 20 HLA regions | calls and benchmark |
| `populations/` | — | `samples.tsv`: 1000 Genomes population and sex of the graphs' samples; `population_pcs.tsv`: mean and spread of genotype PC1-PC5 per population, from MAGE's joint PC table ([Taylor et al. 2024](https://doi.org/10.1038/s41586-024-07708-2)); both feed the simulated cohorts of `dev/associate_inputs.sh` |

The reference is the GRCh38 path. The `.snarls.jsonl` files are `vg snarls` output for `bubble --snarls`. `Homo_sapiens.GRCh38.116.loci.gtf.gz` holds the Ensembl 116 GTF's genes that overlap these loci (within 50 kb), every line of each; it gives `call --gtf` the same results as the whole file, which is too large to commit.

HLA graphs are named after the genes they cover on the GRCh38 path (chr6, 1-based):

| graph | region | genes |
|---|---|---|
| `hla_f` | 29,711,814-29,738,528 | HLA-F (partly) |
| `hla_g` | 29,762,429-29,843,545 | HLA-G; pseudogenes HLA-P, HLA-V |
| `hla_a` | 29,863,213-29,974,447 | HLA-A; pseudogenes HLA-H, -K, -T, -U, -W |
| `hla_j` | 29,985,061-30,020,220 | pseudogene HLA-J |
| `hla_l` | 30,259,163-30,293,014 | pseudogene HLA-L |
| `hla_n` | 30,346,483-30,357,061 | pseudogene HLA-N |
| `hla_e` | 30,480,903-30,496,852 | HLA-E |
| `hla_c` | 31,250,837-31,286,044 | HLA-C |
| `hla_b` | 31,335,469-31,369,490 | HLA-B |
| `hla_s` | 31,377,181-31,387,181 | pseudogene HLA-S |
| `hla_dra` | 32,437,462-32,459,814 | HLA-DRA |
| `hla_drb1_drb5` | 32,484,347-32,603,538 | HLA-DRB1, HLA-DRB5; pseudogene HLA-DRB6 |
| `hla_dqa1` | 32,627,874-32,662,782 | HLA-DQA1; HLA-DQB1 (partly) |
| `hla_dqa2` | 32,738,592-32,750,685 | HLA-DQA2 |
| `hla_dqb2` | 32,754,398-32,764,815 | HLA-DQB2 |
| `hla_dmb` | 32,932,075-32,943,563 | HLA-DMB (partly) |
| `hla_dma` | 32,948,433-32,969,322 | HLA-DMA |
| `hla_doa` | 33,001,667-33,011,912 | HLA-DOA |
| `hla_dpa1_dpb1` | 33,061,247-33,105,579 | HLA-DPA1, HLA-DPB1; pseudogene HLA-DPA2 |
| `hla_dpb2` | 33,111,577-33,129,591 | pseudogene HLA-DPB2 |
