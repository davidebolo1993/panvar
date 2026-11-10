#!/usr/bin/env bash
# Association on two simulated cohorts (tests/synthetic_data/associate, made by dev/associate_inputs.sh): panvar's
# records beside SNVs and SVs called from the same assemblies, with a permutation threshold. Four runs per locus:
#   pcs      covariates sex and PC1-PC5; each source's rare records also tested together
#   no_pcs   sex only, so ancestry confounding is left in
#   null     pcs plus the simulated genetic part and population offsets as covariates: a false hit is expected
#            in about 5% of such runs
#   alleles  as pcs, on the haplotype-level matrices, each row also as the min and max of a sample's haplotypes
#   scripts/associate_examples.sh [panvar binary] [output dir] [data dir]
set -euo pipefail
here=$PWD
cd "$(dirname "$0")/.."
# Relative paths given as arguments are read from the calling directory.
arg() { if [ -z "$1" ]; then echo "$2"; elif [ "${1#/}" != "$1" ]; then echo "$1"; else echo "$here/$1"; fi; }
# Run a step with stderr to a log; on failure show the log's end and stop.
run() { local log=$1; shift; "$@" 2> "$log" || { echo "failed: $* (log: $log)"; tail -n 3 "$log"; exit 1; }; }

command -v Rscript > /dev/null || { echo "Rscript not found on PATH"; exit 1; }
PANVAR=$(arg "${1:-}" build/panvar)
OUT=$(arg "${2:-}" tests/examples/associate)
DATA=$(arg "${3:-}" tests/synthetic_data/associate)

example() {  # <name>
  local name=$1 d=$DATA/$1 dir=$OUT/$1
  [ -f "$d/cohort.samples.tsv" ] || { echo "no $d/cohort.samples.tsv: run dev/associate_inputs.sh first"; exit 1; }
  mkdir -p "$dir/associate"
  # The three matrices, at the sample and at the haplotype level.
  local ms=() mh=() layer m
  for layer in variant snv svim; do
    m=$d/sample/$layer
    ms+=(--genotypes "$m/bimbam_$layer.bimbam.gz" --samples "$m/samples.txt.gz" --feature-annot "$m/feature_annot.$layer.tsv.gz")
    m=$d/haplotype/$layer
    mh+=(--genotypes "$m/bimbam_$layer.bimbam.gz" --samples "$m/samples.txt.gz" --feature-annot "$m/feature_annot.$layer.tsv.gz")
  done
  local tests=(pcs:cohort.phenotype.tsv no_pcs:cohort.phenotype.no_pcs.tsv null:cohort.phenotype.null.tsv alleles:cohort.phenotype.tsv) t
  for t in "${tests[@]}"; do
    local args=("${ms[@]}")
    case ${t%%:*} in
      pcs) args+=(--group layer) ;;
      alleles) args=("${mh[@]}" --haplotypes "$d/cohort.samples.tsv" --encoding sum,min,max --group layer) ;;
    esac
    run "$dir/associate/${t%%:*}.log" "$PANVAR" associate "${args[@]}" --phenotype "$d/${t#*:}" --permutations 1000 \
      -o "$dir/associate/${t%%:*}"
    [ -s "$OUT/sources.tsv" ] || printf 'locus\trun\t%s\n' "$(head -n 1 "$dir/associate/${t%%:*}.layers.tsv")" > "$OUT/sources.tsv"
    awk -F'\t' -v n="$name" -v r="${t%%:*}" 'NR > 1 { print n "\t" r "\t" $0 }' "$dir/associate/${t%%:*}.layers.tsv" >> "$OUT/sources.tsv"
  done

  # Per run: each source's best row (partial r2, role after conditioning) and the signals selected.
  for t in "${tests[@]}"; do
    local r=${t%%:*}
    echo "$name $r: $(awk -F'\t' 'NR > 1 { printf "%s%s %.2f %s", (NR > 2 ? "; " : ""), $1, $10, $13 }' \
      "$dir/associate/$r.layers.tsv"); signals $(awk -F'\t' '$1 == "cojo_independent_signals" { print $2 }' \
      "$dir/associate/$r.summary.tsv"), rows past the permutation threshold $(awk -F'\t' \
      '$1 == "significant_permutation" { print $2 }' "$dir/associate/$r.summary.tsv")"
  done
}

mkdir -p "$OUT"
: > "$OUT/sources.tsv"
example gstm
example lpa
Rscript scripts/plot_associate.R --out "$OUT/associate.region.png" gstm="$OUT/gstm/associate/pcs" lpa="$OUT/lpa/associate/pcs" > /dev/null
echo "per-source summary: $OUT/sources.tsv; plot: $OUT/associate.region.png"
