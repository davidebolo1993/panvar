# S2-L: bounded cross-locus validation and performance-baseline closure

Branch `s2-l-cross-locus` from `fe9025d`. No C++ change: validation and consolidation only.

## 1. Comparator v3 (`s2l_compare.py`)
Every earlier comparator and its failed results are unchanged.

**What v3 protects:**
- every scientific output file, byte for byte;
- the rounds' scientific columns;
- every effective setting **value**, except explicitly declared settings; the key sets must be equal;
- every status key not named by an enabled category, including configured caps and budgets and every prefilter configuration row.

**Categories:**
- **`operational`:** timing; process resources; scheduling-variable emission and DP counters; and only the five audited prefilter-cache counters (consumed hits, completed-block consumed hits, misses, decisions held, decisions not retained). S2-K showed their split and admission depend on scheduling, while decisions are pure functions of content.
- **`profile`:** profile path, name and SHA-256, accepted **only if** the effective-settings check passed.
- **`route`:** blocked versus ordinary union rows and union execution work.
- **`admission`:** the round-1 policy rows, refinement work and margins, and the origin-cache end bytes.
- **`checkpoints`:** footprint checkpoint rows.
- **`envelope`:** `footprint_limit_bytes`.

**Controls** (`compare_controls_v3.py`, on relocated copies of the saved S2-K 6 GiB ACOT 30× run):
- negative controls must fail for their own reason: score, candidate, trajectory, filter, prediction, effective scientific setting (`mismatch_band`), resource cap in settings (`max_physical_origins`), resource cap in status (`max_union_cached_prefilters_consumed`), prefilter configuration (cache budget), undeclared setting, undeclared status value, missing field, new field;
- positive controls must pass: an audited prefilter counter, a timing row;
- it also recompares the saved S2-K 6 GiB run with S2-J's completed run (the S2-K failure), declaring `footprint_limit_bytes`.

Fixes to S2-K tools: the controls' label column is no longer overwritten; trailing whitespace removed from `DESIGN.md`.

## 2–3. Frozen matrix (`s2l_run.py`)
**Cells:** C4, CYP2D6 and ACOT at 5× and 30×, using the existing development inputs (P1 cells). No new samples, truth or tuning.

**Every run:**
- the cell's P1 'opt' profile: its scientific settings and every work and consumption cap;
- `sequence_first_v2`;
- S2-A `early_incumbent experimental_v1`;
- `footprint_limit_bytes 6442450944`;
- time limit = declared timeout − 20 s;
- 4 threads, one process at a time;
- parked options off.

**The two arms:**
- **control:** ordinary (unblocked) indexed union; existing origin admission;
- **combined:** adds `union_fragment_blocks fine_core_4096_experimental` and `repair_round1_origin_admission reuse_only_experimental`. These two settings are the only declared differences.

**Pairs:**
- compared by v3 with categories `operational,profile,route,admission,checkpoints` only when both arms completed (COMPLETE, CONDITIONAL_GENOTYPE_EQUIVALENT or BLOCK_CONSENSUS);
- a refusal is reported as a resource limitation; a genuine tie is a completed outcome;
- ACOT 30× runs the combined arm first. If the budget permits, the combined ACOT 30× is repeated once and compared with itself (categories `operational,profile,checkpoints`).

**Summary:** `s2l_summary.py` reports, per run:
- outcome and first binding limit;
- prediction, completion, phase support;
- total, wall, union, refinement and post-search time;
- process peak;
- rounds, moves, branching, library sensitivity;
- executed and served work, and margins;
- admission counts;
- identities: model fingerprint, input digest, effective-settings digest, prepared locus.
