#!/usr/bin/env bash
# DIAGNOSTIC (v2/post-hoc): the tie-branching diagnostic on named C4 production-v1 improving-move-tie cells,
# one at a time, then truth evaluation of each diagnostic representative after the truth-free tables are hashed.
#
#   tests/genotype_c4_tie_diagnostic.sh <build dir> <campaign dir> <out dir> <job> [job ...]
#
# Production outputs are written unchanged; the diagnostic's go to <prefix>.diag.*. Inputs are the campaign's
# bound inputs, re-verified. Resumable: a cell whose diagnostic finished is not rerun. Extra tie-branch
# arguments (e.g. --trace-branches 1) may be passed in TIE_BRANCH_ARGS.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${1:?usage: genotype_c4_tie_diagnostic.sh <build> <campaign> <out> <job>...}"
CAMPAIGN="${2:?}"; OUT="${3:?}"; shift 3
[ "$#" -gt 0 ] || { echo "name at least one job" >&2; exit 1; }
PROBE="$BUILD/genotype_probe"
PROFILE="$REPO/profiles/genotype_c4_campaign_local.tsv"
REF='HG00096#1#haplotype1-0000024:31848049-32074186'
read -r -a EXTRA <<< "${TIE_BRANCH_ARGS:-}"
mkdir -p "$OUT"
sha() { shasum -a 256 "$1" | cut -d' ' -f1; }
bound() {
  local want; want="$(awk -F'\t' -v p="$1" '$1==p{print $2}' "$CAMPAIGN/inputs.tsv")"
  [ -n "$want" ] && [ "$(sha "$1")" = "$want" ] || { echo "REFUSED: $1 is not the bound input" >&2; exit 1; }
}
for job in "$@"; do
  row="$(awk -F'\t' -v j="$job" '$1==j' "$CAMPAIGN/jobs.tsv")"
  [ -n "$row" ] || { echo "no job $job" >&2; exit 1; }
  IFS=$'\t' read -r _ _ _ _ _ _ model r1 r2 lib cur _ <<< "$row"
  d="$OUT/$job"; mkdir -p "$d"
  if [ -f "$d/call.diag.status.tsv" ] && ! grep -q $'^state\tRUNNING' "$d/call.diag.status.tsv"; then continue; fi
  CUR=(); [ "$cur" != none ] && CUR=(--panel-curation "$cur")
  for f in "$model" "$r1" "$r2" "$lib"; do bound "$f"; done
  [ "$cur" != none ] && bound "$cur"
  echo "[$(date -u +%H:%M:%S)] $job"
  /usr/bin/time -l "$PROBE" tie-branch --profile "$PROFILE" --threads 8 --gfa "$model" --reference "$REF" --reads1 "$r1" \
    --reads2 "$r2" --library "$lib" ${CUR[@]+"${CUR[@]}"} --out-prefix "$d/call" ${EXTRA[@]+"${EXTRA[@]}"} > "$d/log" 2> "$d/time"
  echo "  $(tail -1 "$d/log")"
done
( cd "$OUT" && find . -name 'call.diag.*.tsv' | sort | xargs shasum -a 256 ) > "$OUT/truth_free_tables.sha256"
for job in "$@"; do
  row="$(awk -F'\t' -v j="$job" '$1==j' "$CAMPAIGN/jobs.tsv")"
  IFS=$'\t' read -r _ _ _ _ scheme label model _ _ _ cur truth _ <<< "$row"
  d="$OUT/$job"
  CUR=(); [ "$cur" != none ] && CUR=(--panel-curation "$cur")
  decision="$(awk -F'\t' '$1=="decision"{print $2}' "$d/call.diag.status.tsv" 2>/dev/null)"
  case "$decision" in DIAGNOSTIC_UNIQUE|DIAGNOSTIC_GENOTYPE_EQUIVALENT) ;; *) continue ;; esac
  [ -f "$d/call.diag_eval.summary.tsv" ] && continue
  rep="$(awk -F'\t' '$1=="representative_hap1"{a=$2} $1=="representative_hap2"{b=$2} END{print a"|"b}' "$d/call.diag.status.tsv")"
  "$PROBE" evaluate --scheme "$scheme" --donor "$label" --reference "$REF" --run-prefix "$d/call" --model-gfa "$model" \
    --truth-gfa "$truth" ${CUR[@]+"${CUR[@]}"} --fallback-alleles "$rep" \
    --evaluated-label diagnostic_tie_branching_v2_posthoc_not_a_call --out "$d/call.diag_eval" > "$d/eval.log" 2>&1
  rm -f "$d"/*.sorted.gfa
done
( cd "$OUT" && shasum -a 256 -c truth_free_tables.sha256 > /dev/null ) && echo "truth-free tables unchanged by evaluation" \
  || { echo "DEFECT: evaluation changed a truth-free table" >&2; exit 1; }
