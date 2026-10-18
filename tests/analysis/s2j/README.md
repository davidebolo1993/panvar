# S2-J: reuse-aware origin-cache admission in refinement round 1 (experimental, default off)

- **Setting:** profile key `repair_round1_origin_admission reuse_only_experimental`. In round 1 of the production search, an origin-window entry is retained only if the round's complete plan requests it for two or more distinct new haplotypes. Later rounds, tie branches and legs keep the existing policy.
- **Design:** `DESIGN.md`. **Gates:** `panvar_genotype_round1_admission`, `test_round1_admission_caller`. **Tools:** runner `s2j_run.py`, comparator `s2j_compare.py`.

**Outcome:**
- **ACOT 30× replay:** round 1 commits at a 3.59 GB peak, against the existing policy's > 4.29 GB refusal. The origin cache holds 155 vs about 795 MB, with no extra round-1 enumeration.
- **Full ACOT 30× caller:** COMPLETE in 103 s. Round 1 commits at 4.10 GB, but the process peak is 4.36 GB, reached after the run claimed itself (post-refinement output).
- **C4 and ACOT 5×:** science identical and zero whole-run recomputation (both converge in round 1). The compare units failed only on 4 status keys missing from the frozen allowances.

Results: `tests/results/s2j_origin_admission_20261018/README.md`.
