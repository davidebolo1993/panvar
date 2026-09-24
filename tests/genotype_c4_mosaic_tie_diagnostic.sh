#!/usr/bin/env bash
# DIAGNOSTIC (v2/post-hoc): the 16 synthetic-mosaic no-call cells of C4 production-v1, one at a time.
#
#   tests/genotype_c4_mosaic_tie_diagnostic.sh <build dir> <campaign dir> <out dir>
#
# The 15 improving-move ties run the production call with the tie-branching diagnostic on (production
# outputs unchanged, the diagnostic's in <prefix>.diag.*). The one union-cap refusal (dense_switching_1 30x)
# is rerun once with an EXPERIMENT-ONLY doubled union physical-origin cap, recorded as an override; no
# profile or default moves. Truth enters only afterwards: stage 2 evaluates each diagnostic representative
# (and the union-cap rerun's call) against the plan's staged truth, after the truth-free tables are hashed.
# Inputs are the campaign's own bound inputs (inputs.tsv), re-verified before use. Resumable: a cell with a
# finished diagnostic status is not rerun.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${1:?usage: genotype_c4_mosaic_tie_diagnostic.sh <build> <campaign> <out>}"
CAMPAIGN="${2:?}"
OUT="${3:?}"
PROBE="$BUILD/genotype_probe"
PROFILE="$REPO/profiles/genotype_c4_campaign_local.tsv"
REF='HG00096#1#haplotype1-0000024:31848049-32074186'
UNION_CAP_EXPERIMENT=1400000000
mkdir -p "$OUT"
sha() { shasum -a 256 "$1" | cut -d' ' -f1; }
bound() {   # the campaign bound this input by content; refuse if it moved
  local want; want="$(awk -F'\t' -v p="$1" '$1==p{print $2}' "$CAMPAIGN/inputs.tsv")"
  [ -n "$want" ] && [ "$(sha "$1")" = "$want" ] || { echo "REFUSED: $1 is not the bound input" >&2; exit 1; }
}
cells="$(awk -F'\t' 'NR>1 && $2=="mosaic" && $5=="no_call"{print $1"\t"$6}' "$CAMPAIGN/completeness.tsv" | sort)"
[ "$(echo "$cells" | wc -l | tr -d ' ')" = "16" ] || { echo "expected 16 mosaic no-call cells" >&2; exit 1; }
if [ "$(uname)" = "Darwin" ]; then TIME=(/usr/bin/time -l); else TIME=(/usr/bin/time -v); fi

# ---- stage 1: truth-free
while IFS=$'\t' read -r job reason; do
  row="$(awk -F'\t' -v j="$job" '$1==j' "$CAMPAIGN/jobs.tsv")"
  IFS=$'\t' read -r _ _ _ _ _ _ model r1 r2 lib cur _ <<< "$row"
  d="$OUT/$job"; mkdir -p "$d"
  if [ -f "$d/call.diag.status.tsv" ] && ! grep -q $'^state\tRUNNING' "$d/call.diag.status.tsv"; then continue; fi
  if { [ "$reason" = union_scoring_cap ] || [ "$reason" = rescore_incomplete ]; } && \
     grep -q '^settings_overridden_after_profile.*max_physical_origins=' "$d/call.mosaic_status.tsv" 2>/dev/null; then continue; fi
  for f in "$model" "$r1" "$r2" "$lib" "$cur"; do bound "$f"; done
  extra=()
  # The campaign ledger names this pre-repair refusal rescore_incomplete (its consolidation, union_scoring_cap).
  [ "$reason" = union_scoring_cap ] || [ "$reason" = rescore_incomplete ] && extra=(--union-max-physical-origins "$UNION_CAP_EXPERIMENT")
  echo "[$(date -u +%H:%M:%S)] $job ($reason)"
  "${TIME[@]}" "$PROBE" tie-branch --profile "$PROFILE" --threads 8 --gfa "$model" --reference "$REF" --reads1 "$r1" \
    --reads2 "$r2" --library "$lib" --panel-curation "$cur" --out-prefix "$d/call" ${extra[@]+"${extra[@]}"} > "$d/log" 2> "$d/time"
  echo "  $(tail -1 "$d/log")"
done <<< "$cells"

# ---- truth-free tables frozen before truth is loaded
( cd "$OUT" && find . -name 'call.diag.*.tsv' | sort | xargs shasum -a 256 ) > "$OUT/truth_free_tables.sha256"

# ---- stage 2: evaluation, truth only now
while IFS=$'\t' read -r job reason; do
  row="$(awk -F'\t' -v j="$job" '$1==j' "$CAMPAIGN/jobs.tsv")"
  IFS=$'\t' read -r _ _ _ _ scheme label model r1 r2 lib cur truth _ <<< "$row"
  d="$OUT/$job"
  if [ "$reason" = union_scoring_cap ] || [ "$reason" = rescore_incomplete ]; then
    # The experiment-only rerun's own outcome, whatever it is: evaluate it as its status says.
    "$PROBE" evaluate --scheme "$scheme" --donor "$label" --reference "$REF" --run-prefix "$d/call" --model-gfa "$model" \
      --truth-gfa "$truth" --panel-curation "$cur" --out "$d/call.eval" > "$d/eval.log" 2>&1
  fi
  if [ -f "$d/call.diag.status.tsv" ]; then
    decision="$(awk -F'\t' '$1=="decision"{print $2}' "$d/call.diag.status.tsv")"
    case "$decision" in DIAGNOSTIC_UNIQUE|DIAGNOSTIC_GENOTYPE_EQUIVALENT) ;; *) rm -f "$d"/*.sorted.gfa; continue ;; esac
    rep="$(awk -F'\t' '$1=="representative_hap1"{a=$2} $1=="representative_hap2"{b=$2} END{print a"|"b}' "$d/call.diag.status.tsv")"
    "$PROBE" evaluate --scheme "$scheme" --donor "$label" --reference "$REF" --run-prefix "$d/call" --model-gfa "$model" \
      --truth-gfa "$truth" --panel-curation "$cur" --fallback-alleles "$rep" \
      --evaluated-label diagnostic_tie_branching_v2_posthoc_not_a_call --out "$d/call.diag_eval" > "$d/eval.log" 2>&1
  fi
  rm -f "$d"/*.sorted.gfa
done <<< "$cells"
( cd "$OUT" && shasum -a 256 -c truth_free_tables.sha256 > /dev/null ) && echo "truth-free tables unchanged by evaluation" \
  || { echo "DEFECT: evaluation changed a truth-free table" >&2; exit 1; }
python3 "$REPO/scripts/c4_mosaic_tie_summary.py" --campaign "$CAMPAIGN" --out "$OUT"
