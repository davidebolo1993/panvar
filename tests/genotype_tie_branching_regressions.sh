#!/usr/bin/env bash
# v2 production tie branching: five frozen C4 regressions, compared with their frozen v2/post-hoc diagnostics.
#
#   tests/genotype_tie_branching_regressions.sh <build dir> <campaign dir> <diagnostics_v2 dir> <out dir> [extra panvar args]
#
# Each cell was an improving-move-tie no-call of C4 production-v1. `panvar genotype` with the v2 profile now
# branches and pools it itself; the result must equal what the frozen diagnostic reported for that cell:
#   outcome   DIAGNOSTIC_UNIQUE -> COMPLETE strict_call (FASTA written); DIAGNOSTIC_GENOTYPE_EQUIVALENT ->
#             CONDITIONAL_GENOTYPE_EQUIVALENT (no FASTA); search_conditional tie_branched
#   pool      the same identities, no score mismatch; the same pooled top set, member for member and score for
#             score at 17 digits, and the same representative
#   branches  the same branch starts and endpoints, and the same terminal ties
#   repair    the production trajectory equal to the frozen campaign job's
# Truth is not used. Inputs are the campaign's bound inputs, re-verified.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${1:?usage: genotype_tie_branching_regressions.sh <build> <campaign> <diagnostics_v2> <out> [extra args]}"
CAMPAIGN="${2:?}"; DIAG="${3:?}"; OUT="${4:?}"; shift 4
PANVAR="$BUILD/panvar"
PROFILE="$REPO/profiles/genotype_c4_v2_campaign_local.tsv"
REF='HG00096#1#haplotype1-0000024:31848049-32074186'
CELLS="mosaic.both_homologues_2.2:mosaic_ties_c5d341a mosaic.dense_switching_1.5:mosaic_ties_c5d341a
loo.HG01505.5:loo_ties_648476c loo.HG01596.2:loo_ties_648476c loo.HG01573.5:loo_ties_648476c"
mkdir -p "$OUT"
fail=0
same() { if [ "$2" = "$3" ]; then echo "  ok: $1"; else echo "  REGRESSION: $1 ($2 vs $3)"; fail=1; fi; }
key() { awk -F'\t' -v k="$1" '$1==k{print $2; exit}' "$2"; }
sha() { if [ "$#" -gt 0 ]; then shasum -a 256 "$1"; else shasum -a 256; fi | cut -d' ' -f1; }
if [ "$(uname)" = "Darwin" ]; then TIME=(/usr/bin/time -l); else TIME=(/usr/bin/time -v); fi
for spec in $CELLS; do
  job="${spec%%:*}"; frozen="$DIAG/${spec#*:}/$job/call"
  row="$(awk -F'\t' -v j="$job" '$1==j' "$CAMPAIGN/jobs.tsv")"
  IFS=$'\t' read -r _ _ _ _ _ _ model r1 r2 lib cur _ _ _ _ _ _ prefix <<< "$row"
  for f in "$model" "$r1" "$r2" "$lib"; do
    [ "$(sha "$f")" = "$(awk -F'\t' -v p="$f" '$1==p{print $2}' "$CAMPAIGN/inputs.tsv")" ] || { echo "REFUSED: $f is not the bound input"; exit 1; }
  done
  CUR=(); [ "$cur" != none ] && CUR=(--panel-curation "$cur")
  d="$OUT/$job"; rm -rf "$d"; mkdir -p "$d"
  echo "== $job"
  "${TIME[@]}" "$PANVAR" genotype --profile "$PROFILE" --threads 8 --gfa "$model" --reference-path "$REF" \
    --reads1 "$r1" --reads2 "$r2" --library "$lib" ${CUR[@]+"${CUR[@]}"} --out-prefix "$d/call" --quiet "$@" \
    > "$d/log" 2> "$d/time"
  s="$d/call.mosaic_status.tsv"
  echo "  wall $(awk '/ real/{print $1}' "$d/time") s, tie branching $(key seconds_tie_branching "$s") s, peak footprint $(awk '/peak memory footprint/{printf "%.2f GiB", $1/2^30}' "$d/time")"
  decision="$(key decision "$frozen.diag.status.tsv")"
  case "$decision" in
    DIAGNOSTIC_UNIQUE) want_status=COMPLETE; want_class=strict_call; want_fasta=1 ;;
    DIAGNOSTIC_GENOTYPE_EQUIVALENT) want_status=CONDITIONAL_GENOTYPE_EQUIVALENT; want_class=genotype_equivalent_call; want_fasta=0 ;;
    *) want_status=UNEXPECTED; want_class=UNEXPECTED; want_fasta=0 ;;
  esac
  same "status" "$(key status "$s")" "$want_status"
  same "call class" "$(key call_class "$s")" "$want_class"
  same "search_conditional" "$(key search_conditional "$s")" "tie_branched"
  same "tie branch order" "$(key tie_branch_order "$s")" "bfs_canonical_identity_v1"
  same "FASTA written" "$([ -f "$d/call.mosaic_haplotypes.fa" ] && echo 1 || echo 0)" "$want_fasta"
  same "branches run" "$(key tie_branches_run "$s")" "$(key branches_run "$frozen.diag.status.tsv")"
  same "ties" "$(key tie_ties "$s")" "$(key ties "$frozen.diag.status.tsv")"
  same "pool identities" "$(key pool_identities "$s")" "$(key pool_identities "$frozen.diag.status.tsv")"
  same "pool score mismatches" "$(key pool_score_mismatches "$s")" "0"
  same "pooled top set, member and score" \
    "$(awk -F'\t' 'NR>1{print $5"\t"$6"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' "$d/call.equivalence_set.tsv" | sort | sha)" \
    "$(awk -F'\t' 'NR>1{print $3"\t"$4"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' "$frozen.diag.pooled_top.tsv" | sort | sha)"
  same "representative" "$(awk -F'\t' '$4=="true"{print $5"|"$6; exit}' "$d/call.equivalence_set.tsv")" \
    "$(key representative_hap1 "$frozen.diag.status.tsv")|$(key representative_hap2 "$frozen.diag.status.tsv")"
  same "branch starts and endpoints" \
    "$(awk -F'\t' 'NR>1{print $3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' "$d/call.tie_branches.tsv" | sort | sha)" \
    "$(awk -F'\t' 'NR>1{print $3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14}' "$frozen.diag.branches.tsv" | sort | sha)"
  same "terminal ties" \
    "$(awk -F'\t' 'NR>1{print $1"\t"$2"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$15"\t"$16"\t"$17}' "$d/call.tie_terminal_ties.tsv" | sort | sha)" \
    "$(awk -F'\t' 'NR>1{print $1"\t"$2"\t"$4"\t"$5"\t"$7"\t"$8"\t"$9"\t"$10"\t"$12"\t"$16"\t"$18"\t"$19}' "$frozen.diag.terminal_ties.tsv" | sort | sha)"
  same "production trajectory" "$(cut -f1-28 "$d/call.repair_rounds.tsv" | sha)" "$(cut -f1-28 "$prefix.repair_rounds.tsv" | sha)"
  rm -f "$d"/*.sorted.gfa
done
[ "$fail" -eq 0 ] && echo "tie-branching regressions: PASS" || { echo "tie-branching regressions: FAILED"; exit 1; }
