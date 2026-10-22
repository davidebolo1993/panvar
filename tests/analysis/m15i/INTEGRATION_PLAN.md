# M15-I integration plan (non-mutating; nothing below has been executed)

## 1. Revisions

| role | ref | commit | relation |
|---|---|---|---|
| **functional freeze** | `m15-h-help` | `cd51d37` | the last code change. `panvar` SHA-256 `3f7e03a41c9d507d8e4600df819140202aa069b8feb717cad67abc96b17b462b`; full suite 64/64 |
| **release candidate** | `m15-i-release` | the M15-I head (documentation only; see §6) | `cd51d37` plus release notes and this plan. **No code, test or CMake change**, verifiable with `git diff --stat cd51d37 m15-i-release -- src include tests/*.cpp tests/*.sh tests/*.py CMakeLists.txt` (empty) |
| **intended target** | `main` = `origin/main` | `e91b73c` | `1d87162` plus **one** commit, the panphorte fix `e91b73c`. The candidate is 570 commits ahead and 1 behind (merge base `1d87162`) |
| genotype release branch | `genotype-mosaic` = `origin/genotype-mosaic` | `cee58de` | the candidate descends from it (330 ahead, 0 behind): a fast-forward |

## 2. How the panphorte fix is retained

The fix `e91b73c` changes `src/panphorte.cpp`, `tests/panphorte_stats.sh` and `docs/algorithms/panphorte.md`. **The
candidate line never touches these three files** (0 commits in `1d87162..cd51d37`).

**A dry merge confirms it is clean:** `git merge-tree --write-tree main cd51d37` exited 0 with no conflicts. That command
writes unreferenced objects only, and changes no branch, index, worktree or ref. Its tree `1e1547af…`:
- takes `main`'s `src/panphorte.cpp` (blob `ca07e8c`, with the fragment-sharing code);
- carries the new `panphorte_stats` fixture;
- registers both `panvar_panphorte_stats` and `panvar_genotype_assessment`.

**Recommended method:** a **merge commit** of the candidate into `main` (`--no-ff`). The fix stays as its own commit in
history; no cherry-pick or rebase is needed.

## 3. Expected conflicts and affected components

**Textual conflicts:** none expected (the dry merge is clean).

**Semantic interaction:** none expected.
- panphorte (the fix) and genotype (the candidate) share no source file.
- The genotype command does not call panphorte.
- The suite runs `panphorte_stats` with the new fixture.

**What lands on `main`:**
- **The full mosaic-caller line since `1d87162`:** 768 files, about 171 k insertions.
  - **Code:** `src/genotype_*`, `src/prepared_locus*`, `src/library_estimate*`, `include/panvar` (24 headers).
  - **Tests:** C++ and script tests, and the CMake registrations.
  - **Docs:** `docs/modules/genotype.md`, `docs/algorithms/genotype.md`, `docs/guides/genotype_short_read.md`,
    `docs/releases/`.
  - **`scripts/genotype_report/`.**
  - **`tests/analysis/`:** 532 files of batch designs, harnesses and results notes. This is the largest share.
- **Not included:** `tests/results/` is git-ignored, so the sealed evidence packages stay local and unpublished.

**Hazard: the primary checkout.** `/Users/davide.bolognini/panvar` is on `genotype-mosaic` (`cee58de`) with an uncommitted
WIP edit of `tests/genotype_probe.cpp` (SHA-256 prefix `311a41adc1eddc0c`).
- **The candidate line changes that same file** in 24 commits (+3,295 / −163).
- **So `genotype-mosaic` must not be fast-forwarded in the primary checkout** until the WIP is reconciled. Stash `dbb170c5` (29 September) holds a related but **different** version of this WIP (SHA-256 prefix
  `31119144`, not the working copy's `311a41ad`); both are kept.
- **The integration below uses a separate worktree** and never touches the primary checkout.

**Occupied worktrees:**
- `panvar-worktrees/main-1d87162` has `main` checked out (clean apart from three untracked build logs).
- `panvar-worktrees/s1-d-panel-references` has the candidate line.
- **Neither is switched.** A new worktree is used.

## 4. Steps (each needs Davide's authorization; none executed)

```
# 1. integration worktree on a new branch from main (no occupied worktree is switched)
git worktree add ../panvar-worktrees/integrate-m15 -b integrate/genotype-sr-exp-1 main
cd ../panvar-worktrees/integrate-m15
git merge --no-ff m15-i-release -m "Merge genotype short-read experimental release 1 (M15-I) into main"
# 2. build and test (section 5)
# 3. after approval: fast-forward main to the merge, push, tag (section 6)
git push origin integrate/genotype-sr-exp-1:main            # main: e91b73c -> merge commit (fast-forward)
git push origin m15-i-release:genotype-mosaic               # optional: the genotype branch, fast-forward from cee58de
git tag -a v0.2.0-exp.1 <merge commit> -m "panvar 0.2.0-exp.1: genotype short-read experimental release 1"
git push origin v0.2.0-exp.1
```

**The local `main` branch is checked out in `main-1d87162`.** It is updated afterwards in that worktree, by Davide or with
his explicit authorization (`git -C ../panvar-worktrees/main-1d87162 merge --ff-only integrate/genotype-sr-exp-1`).

## 5. Checks required on the integrated tree (before any push)

1. **A clean configure and build of every target** (`cmake -DCMAKE_BUILD_TYPE=Release`, testing on), with
   `git status --porcelain` empty. Record the `panvar` SHA-256.
2. **The full enabled suite, serially, 900 s per test.** Expect every enabled test to pass, including:
   - `panvar_panphorte_stats`, with the shared-fragment fixture;
   - `panvar_refine_stats`;
   - `panvar_genotype_assessment`;
   - `panvar_genotype_cli_contract`, with the guide-inventory check;
   - `panvar_genotype_command_e2e`.

   Report the disabled tests separately.
3. **Genotype outputs are unchanged by the merge.** The no-profile C4, CYP2D6 and ACOT development commands are identical
   to the candidate's runs (M15-H `a2.ex.preserved_c4`, M15-G `a3.ex.complete_*`), apart from names, timing and the exact
   caches' scheduling-dependent counters.
4. **Panphorte behaves as the fix intends.** The chr3:46654531–46940288 pipeline (rebuild → bubble → panphorte 0.97 →
   refine) gives:
   - 6 fragment nodes at bubble 5;
   - identical walks for identical haplotypes;
   - 466 of 466 paths spelling the same bases with unchanged REP copy counts.

   This is as recorded when the fix was made (input: the colleague's GFA).
5. **Smoke checks of the public interface:**
   - `--help` and `--help-advanced` match the M15-H captures;
   - a `--build-panel-index` / `--panel-index` round trip on C4 (cold = indexed);
   - one `--dosage-check` run;
   - one refused contradiction and one refused missing library.
6. **The platform recorded:** macOS arm64 only (the release notes say so).

**Estimated execution:** about 1,200 s (the suite takes about 1,000 s).

## 6. Proposed release name and version (for Davide's approval)

- **Name:** *panvar genotype: short-read experimental release 1*.
- **Version and tag:** **`v0.2.0-exp.1`**, an annotated tag on the integration merge commit. There are no earlier release
  tags.
- **The CMake version:** `CMakeLists.txt` declares `project(panvar VERSION 0.1.0)`.
  - **Recommended:** a one-line version bump to `0.2.0` in a release commit on the integration branch, before tagging. That
    needs approval, a rebuild and the suite.
  - **Alternative:** tag without the bump and record the mismatch in the notes.
- **Release notes:** `docs/releases/genotype-short-read-experimental-1.md`.
