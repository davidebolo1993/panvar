#!/usr/bin/env bash
# Data collection for docs/reports/genotype-seed-preregistration.md.
#
# Does the GC slope reproduce across read draws, or was it one seed?
#
# The association was measured once. r-squared is about 0.001, the per-block slopes disagree in sign
# for 2 of 11 blocks, and the reads come from wgsim, which has no GC-selection mechanism at all. Before
# any efficiency model is fitted, the slope has to be shown to exist across draws.
#
# Fitted on ANCHORS only. An anchor is carried at multiplicity 1 by every traversing allele, so dosage
# confounding is structurally absent and no truth alleles are needed. Anchors are also the population
# that actually sets lambda, so a GC effect on them is precisely the one that biases the depth.
#
# SEED and SEED+466k select the same path indices at p=0 while giving different wgsim seeds, so the
# genotype is held fixed and only the reads change. Uses genotype_sim.sh unmodified rather than
# reimplementing the read simulation.
set -uo pipefail

R="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SP="${SEED_DIR:?set SEED_DIR to a scratch directory}"
N="${1:-55}"   # k=0..4 exploratory, k=5..54 the confirmatory 50
mkdir -p "$SP/seeds"

for k in $(seq 0 $((N - 1))); do
  S=$((42 + 466 * k))
  dump="$SP/seeds/m_$S.tsv"
  [[ -s "$dump" ]] && { echo "seed $S: cached"; continue; }
  d="$SP/seeds/run_$S"
  mkdir -p "$d"
  OUT="$d" LOO=1 SEED="$S" PYTHON="$HOME/miniconda3/bin/python" \
    GENOTYPE_EXTRA="--dump-markers $dump" \
    bash "$R/tests/genotype_sim.sh" lpa 1 30 0.001 >"$d/log" 2>&1
  if [[ -s "$dump" ]]; then
    # Keep the columns the pre-registered model needs, for markers this sample actually carries.
    # The full dump is ~5 MB a seed and none of the rest is used.
    awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i
                       print "block\tclass\tgc\tcount\tmult\tupos"; next}
                $(h["truth_mult"])!="NA" && $(h["truth_mult"])>0 && $(h["truth_pos_mean"])!="NA" {
                  print $(h["block_index"])"\t"$(h["marker_class"])"\t"$(h["gc"])"\t"$(h["count"])"\t"$(h["truth_mult"])"\t"$(h["truth_pos_mean"])}' \
        "$dump" > "$dump.keep" && mv "$dump.keep" "$dump"
    echo "seed $S: $(( $(wc -l < "$dump") - 1 )) markers"
  else
    echo "seed $S: FAILED (see $d/log)"
  fi
  rm -rf "$d/fa" "$d"/lpa/r_*.fq.gz 2>/dev/null
done
echo "ALL DONE"
