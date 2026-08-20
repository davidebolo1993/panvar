#!/usr/bin/env bash
# Does --depth-estimator mean help, hurt, or do nothing once the ladder stops being saturated?
#
# At 30x the ladder sits at 36/36 in nearly every case under both estimators, so it can detect a
# regression and cannot detect a gain. Lower depth is where the estimator should matter: anchor counts
# fall, so the integer lattice a median is confined to becomes relatively coarser. At 30x the region
# median is 23 and the lattice step is 1 count, about 4 percent of lambda; at 5x the median is about 4
# and the same step is roughly 25 percent.
#
# CRITERION, fixed before running, because a default change hangs on it:
#
#   flip      the mean never scores worse than the median outside noise at any depth, AND scores
#             better somewhere at 5x or 10x
#   hold      the two are indistinguishable everywhere -- the case for the mean then rests on theory
#             alone (Fisher consistency, the lattice argument, agreement with simulation theory to
#             0.07 percent), which is defensible but is not what this sweep was run to show
#   revert    the mean is worse anywhere it is not clearly noise
#
# "Noise" here is not hand-waved: the ladder is 4 pairs and both regimes, so a one-block difference on
# one case at one depth is within draw-to-draw variation and is not evidence either way. A difference
# that reproduces across depths or across cases is.
#
# The paralogous case matters most. With one allele stripped of markers it can separate heterozygous
# from homozygous ONLY by the absolute count level, which is exactly what the depth estimator sets.
#
# EVERY case carries --twin-divergence, and that is not optional. The generator emits each design twice
# as an EXACT twin by default (defect V15), so holding out a sample's two haplotypes leaves an
# identical one in the panel and leave-one-out has a perfectly representable answer. Measured: with
# exact twins the ladder scores 36/36 and 16/16 at every case and at 30x, 10x AND 5x -- 2.5x per
# haplotype -- because picking an identical panel haplotype barely needs coverage. The ladder was never
# saturated by generous depth; it was saturated because the task is trivial, and lowering depth cannot
# make a trivial task discriminate. Diverging the twin is what turns leave-one-out into an off-panel
# test at all.
#
# Carried caveat: synthetic_bench.sh still shares one wgsim seed across both homologues (defect V14,
# unfixed). Both arms inherit it identically, so the comparison is sound while neither arm is a clean
# off-panel measurement in absolute terms.
#
#   SEED_DIR=<scratch> genotype_depth_ladder.sh [depths...]
set -uo pipefail

R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SP="${SEED_DIR:?set SEED_DIR to a scratch directory}"
export PYTHON="${PYTHON:-$HOME/miniconda3/bin/python}"
DEPTHS=("$@"); [[ ${#DEPTHS[@]} -eq 0 ]] && DEPTHS=(30 10 5)
# Integer SNP count, not a fraction. Passing 0.002 here silently yields an empty result.
TWIN="${TWIN:-40}"

run_case() {
  local name="$1" genflags="$2" depth="$3"
  for EST in median mean; do
    local out="$SP/ladder/${name}_${depth}x_${EST}"
    mkdir -p "$out"
    local summary
    summary=$(GEN_EXTRA="$genflags --twin-divergence $TWIN" GENOTYPE_EXTRA="--depth-estimator $EST" SEED=11 \
              bash "$R/tests/synthetic_bench.sh" "$out" 4 "$depth" 0.001 2>&1 \
              | grep -E "^TOTAL leave" | sed -E 's/TOTAL leave-(zero|one)-out *: */\1 /' | tr '\n' '|')
    printf '%-16s %3sx %-7s %s\n' "$name" "$depth" "$EST" "${summary:-<no TOTAL line -- run failed>}"
  done
}

for d in "${DEPTHS[@]}"; do
  echo "======================================== depth ${d}x"
  run_case clean            ""                          "$d"
  run_case paralogous       "--paralog-del --paralog-ins" "$d"
  run_case segdup           "--segdups 2"                "$d"
  run_case folded           "--dup-fold-divergence 0.02" "$d"
  run_case per-design-vntr  "--dup-per-design-divergence 0.02" "$d"
done
echo "ALL DONE"
