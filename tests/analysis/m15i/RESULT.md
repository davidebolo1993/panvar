# M15-I result: release documentation and a non-mutating integration plan; ready for authorization

**Logistics:**
- **Branch:** `m15-i-release`. Candidate head `5eaca6a`: `cd51d37` (the M15-H functional freeze) plus documentation only.
  There is no code, test or CMake difference.
- **Not done:** no merge, cherry-pick, worktree switch, tag, push, publication or cleanup. WIP (`311a41ad…`), stash
  `dbb170c5` and the sealed results are preserved.

## Deliverables

1. **Release notes:** `docs/releases/genotype-short-read-experimental-1.md`. They cover:
   - the accelerated ordinary command;
   - automatic library estimation;
   - `--dosage-check`;
   - the public panel index;
   - the conditional FASTA and extraction references;
   - `assessment_v1`;
   - help and inventory;
   - the correctness fixes and the panphorte fix;
   - a breaking-change table against `genotype-mosaic` `cee58de`;
   - the retained legacy compatibility.
2. **Evidence and limitations** (same file):
   - the reserved ACOT figures (18 prior-exposed families, simulated 5×; median / q90 / q99 / max excess 0 / 8 / 98 /
     478 bp; 0 catastrophes, with the rate bound stated);
   - the two unflagged errors of about 478 bp (sv:9, HG02572 and NA19468);
   - C4 and CYP2D6 as development-only evidence;
   - the automatic library checked in development only;
   - no calibrated confidence and no reliable-call category;
   - the audit, assessment, resource, API, interface and error limitations;
   - **tested platform: macOS 26.5.2 arm64, Apple clang 21.0.0, CMake 4.3.3, Python 3.14.7 only**.
3. **Checklist closed** (`RELEASE_CHECKLIST.md`):
   - superseded headlines retained;
   - legacy-tie assessments, API expansion, aliases and removals deferred;
   - known limitations retained explicitly.
4. **Integration plan** (`INTEGRATION_PLAN.md`):
   - **Candidate** `5eaca6a` into **target** `main` (= `origin/main` `e91b73c`): 571 ahead, 1 behind, merge base `1d87162`.
   - **A dry `git merge-tree` is clean.** The merged tree keeps `main`'s panphorte fix (blob `ca07e8c`), its fixture, and
     both test registrations. Retention is through a `--no-ff` merge commit.
   - **Expected textual conflicts:** none. The candidate line never touched the fix's three files.
   - **What lands:** the whole mosaic-caller line (768 files; 532 under `tests/analysis`). `tests/results` stays
     git-ignored.
   - **Hazard:** the primary checkout's WIP `tests/genotype_probe.cpp` overlaps 24 candidate commits. So the integration
     uses a new worktree, and `genotype-mosaic` is not fast-forwarded locally.
   - **Required checks:** a clean build; the full suite; genotype outputs identical to the candidate's; the chr3 panphorte
     reproduction; interface smoke checks. About 1,200 s.
   - **Proposed release:** *panvar genotype: short-read experimental release 1*, tag `v0.2.0-exp.1`, with an optional
     CMake `VERSION` bump from 0.1.0 to 0.2.0.

## Ready / blocked

| item | state |
|---|---|
| functionality frozen; suite 64/64 at `cd51d37`; outputs preserved | **ready** |
| release notes, evidence, limitations, platform | **ready** |
| integration plan; dry merge clean; panphorte fix retained | **ready** |
| integration merge plus checks on the merged tree | **blocked: authorization A1** |
| CMake version 0.2.0 | **blocked: decision A2** |
| push `main` (and optionally `genotype-mosaic`) | **blocked: authorization A3**, after A1 |
| tag `v0.2.0-exp.1`; publication | **blocked: authorization A4** |
| reconcile the primary checkout's WIP with the candidate's `genotype_probe.cpp` | **Davide's decision** (does not block A1) |

## Exact authorizations

- **A1.** Create worktree `panvar-worktrees/integrate-m15` on a new branch `integrate/genotype-sr-exp-1` from `main`;
  merge `m15-i-release` with `--no-ff`; run the checks of plan §5. Local only.
- **A2.** Choose between (a) a release commit bumping `project(panvar VERSION 0.1.0)` to `0.2.0` on the integration branch
  (rebuild and suite), and (b) no bump.
- **A3.** After A1 and A2 pass, push the integration branch to `origin/main` (a fast-forward from `e91b73c`); optionally
  push `m15-i-release` to `origin/genotype-mosaic` (a fast-forward from `cee58de`).
- **A4.** Create and push the annotated tag `v0.2.0-exp.1` on the integrated commit, and publish the release notes.

**Note:** a check repeated during this batch exposed a shell trap. In zsh, `"$T:src/…"` is parsed as a parameter modifier
(`:s`), so a `git cat-file` of a variable tree path read nothing. The literal tree ID confirms that the fix is present.

**Stopped for review.**
