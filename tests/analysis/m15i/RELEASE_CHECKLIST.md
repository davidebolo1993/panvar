# Release checklist: genotype short-read experimental release 1 (closed in M15-I)

## Decisions closed as recommended

| item | decision |
|---|---|
| superseded `headline_filters_v2` columns | **retained, unchanged**, for this experimental release (documented as superseded, with their known misreadings) |
| block checks for legacy-tied v2 predictions (NOT_ASSESSED) | **deferred** (new scoring and dosage work) |
| assessment view for the caller-only API | **deferred** (documented) |
| run-guard flags and other aliases (guide §10) | **deferred** (not implemented) |
| compatibility removals (guide §10) | **deferred** (after this release) |
| known limitations (CYP2D6 audit refusal, ACOT 30× equivalence gap, sampled footprint guard, no default time limit, non-transactional run-time errors, platform) | **retained explicitly** in the release notes |

## Ready / blocked

| item | state |
|---|---|
| functionality frozen at M15-H (`cd51d37`); suite 64/64; outputs preserved | ready |
| release notes with evidence and limitations | ready (`docs/releases/genotype-short-read-experimental-1.md`) |
| integration plan; dry merge into `main` clean; panphorte fix retained | ready (`INTEGRATION_PLAN.md`) |
| integration merge and the checks of plan §5 on the merged tree | **blocked:** needs authorization (A1) |
| version bump to 0.2.0 in `CMakeLists.txt` | **blocked:** needs a decision (A2) |
| fast-forward `main` and `genotype-mosaic`; push | **blocked:** needs authorization (A3), after A1 passes |
| tag `v0.2.0-exp.1`; publication | **blocked:** needs authorization (A4) |
| reconcile the primary checkout's WIP (`tests/genotype_probe.cpp`) with the candidate's 24 changes to that file | **blocked:** Davide's decision (not needed for the merge, which uses a separate worktree) |

## Authorizations needed (exact)

- **A1.** "Create worktree `integrate-m15` on a new branch `integrate/genotype-sr-exp-1` from `main`; merge `m15-i-release`
  with `--no-ff`; run the integration plan §5 checks." This is local only: no push, tag or cleanup.
- **A2.** Choose between:
  - **(a)** a release commit bumping `project(panvar VERSION 0.1.0)` to `0.2.0` on the integration branch, with a rebuild
    and the suite; or
  - **(b)** no bump.
- **A3.** After A1 (and A2) pass, push the integration branch to `origin/main` (a fast-forward from `e91b73c`), and
  optionally push `m15-i-release` to `origin/genotype-mosaic` (a fast-forward from `cee58de`).
- **A4.** Create and push the annotated tag `v0.2.0-exp.1` on the integrated commit, and publish the release notes.
