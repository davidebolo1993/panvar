#!/usr/bin/env bash
# Make the inputs of scripts/associate_examples.sh, in tests/synthetic_data/associate/<locus>/: a cohort of 1,000
# samples drawn from the graph's haplotypes within 1000 Genomes populations (simulate_cohort.py); call's records
# described for it; SNVs (paftools) and SVs (svim-asm, merged by truvari) from the haplotypes aligned to the
# reference (assembly_variants.py). Each at the sample and haplotype level, plus the haplotype-level graph
# matrix validate_gemma.sh reads.
#   dev/associate_inputs.sh [panvar binary] [data dir] [refine_examples.sh output dir]
#                           [call_examples.sh output dir] [dir with minimap2, k8, paftools.js, samtools, svim-asm]
#                           [dir with bcftools, tabix, truvari] [work dir]
set -euo pipefail
here=$PWD
cd "$(dirname "$0")/.."
arg() { if [ -z "$1" ]; then echo "$2"; elif [ "${1#/}" != "$1" ]; then echo "$1"; else echo "$here/$1"; fi; }
run() { local log=$1; shift; "$@" 2> "$log" || { echo "failed: $* (log: $log)"; tail -n 3 "$log"; exit 1; }; }
command -v python3 > /dev/null || { echo "python3 not found on PATH"; exit 1; }

PANVAR=$(arg "${1:-}" build/panvar)
DATA=$(arg "${2:-}" tests/synthetic_data/associate)
REFINED=$(arg "${3:-}" tests/examples/refine)
CALLED=$(arg "${4:-}" tests/examples/call)
svim=$(command -v svim-asm || true)
TOOLS=$(arg "${5:-}" "${svim%/svim-asm}")
[ -x "$TOOLS/svim-asm" ] || { echo "no svim-asm: give the directory holding minimap2, k8, paftools.js, samtools and svim-asm"; exit 1; }
truvari=$(command -v truvari || true)
MERGE=$(arg "${6:-}" "${truvari%/truvari}")
[ -x "$MERGE/truvari" ] || { echo "no truvari: give the directory holding bcftools, tabix and truvari"; exit 1; }
WORK=$(arg "${7:-}" build/associate_inputs)
POP=tests/real_data/populations

# <name> <simulate_cohort.py arguments...>
inputs() {
  local name=$1 w=$WORK/$1 d=$DATA/$1; shift
  local gfa=$REFINED/$name/$name.normalized.sorted.gfa call=$CALLED/$name/$name
  [ -f "$call.region.vcf" ] || { echo "no $call.region.vcf: run scripts/call_examples.sh first"; exit 1; }
  mkdir -p "$w" "$d"
  python3 dev/simulate_cohort.py --gfa "$gfa" --populations $POP/samples.tsv --population-pcs $POP/population_pcs.tsv \
    -o "$w/cohort" "$@" > "$w/cohort.log"
  # sex only (no PCs), and the PCs with the simulated genetic part and the population offsets (a null run)
  cut -f1-3 "$w/cohort.phenotype.tsv" > "$w/cohort.phenotype.no_pcs.tsv"
  paste "$w/cohort.phenotype.tsv" <(cut -f5,7 "$w/cohort.truth.tsv") > "$w/cohort.phenotype.null.tsv"
  run "$w/describe.log" "$PANVAR" describe -i "$gfa" -b "$REFINED/$name/$name" -r GRCh38 -o "$w/describe" \
    --variant-vcf "$call.region.vcf" --variant-nodes "$call.variant_nodes.tsv" --samples "$w/cohort.samples.tsv" --threads 2
  run "$w/assembly.log" python3 dev/assembly_variants.py --gfa "$gfa" --reference GRCh38 --tools "$TOOLS" \
    --merge-tools "$MERGE" --out-dir "$w/assembly" --samples "$w/cohort.samples.tsv" --threads 2
  cp "$w"/cohort.samples.tsv "$w"/cohort.phenotype.tsv "$w"/cohort.phenotype.no_pcs.tsv "$w"/cohort.phenotype.null.tsv \
    "$w"/cohort.truth.tsv "$d/"
  local lv
  for lv in sample haplotype; do
    mkdir -p "$d/$lv"
    cp -R "$w/describe/$lv/variant" "$w/assembly/$lv/snv" "$w/assembly/$lv/svim" "$d/$lv/"
  done
  cp -R "$w/describe/haplotype/graph" "$d/haplotype/"
  echo "$name: inputs in $d"
}

# GSTM1 expression follows its copies; Lp(a) falls with the shorter allele's KIV-2 copies, higher in African
# populations at the same copies.
inputs gstm --copies tests/real_data/paralogs/gstm.bed --gene GSTM1
inputs lpa --copies tests/real_data/tandem_repeat/lpa.repeats.tsv --model isoform --rate 0.1 \
  --continent-offset AFR=0.9 --genetic-r2 0.5
