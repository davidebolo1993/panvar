#!/usr/bin/env bash
# The five C4 production regressions: the real `panvar genotype --profile profiles/genotype_c4.tsv`
# command, one invocation per case, run serially (each uses the profile's eight threads).
#
#   tests/genotype_production_c4_regressions.sh <build dir> <work dir> [case ...]
#
# Cases (the expected values come from the multi-start evidence of commit 9a0f4c0, restricted to what
# the production rule pools: the exactly scored union and the panel-anchor branch; nothing here reruns
# multi-start):
#
#   na20509_10      NA20509 LOO 10x   score -151090.979646, distance 114, floor 91
#   na20509_30      NA20509 LOO 30x   score -414918.317265, distance 113, floor 91, excess 22; a unique
#                                     winner within the pool (NOT CONDITIONAL_GENOTYPE_EQUIVALENT), 4.40
#                                     nats to the best different genotype, phase PRIOR_ONLY: the best
#                                     same-genotype alternative has an identical fragment term and is
#                                     7.66 nats lower through the prior alone
#   na19240_30      NA19240 LOO 30x   distance 4 = floor 4, the ctx:6 repair in the trajectory
#   na20509_lzo_30  NA20509 LZO 30x   exact (distance 0), no repair movement
#   bridgeable_30   synthetic mosaic  11/11 blocks exact; the union's mosaic winner retained or reached
#                                     from the panel anchor
#
# Every case must also: reproduce the recorded panel-anchor trajectory in every search column
# (tests/genotype_production_c4_expected/<case>.rounds.tsv), stop converged, pool with zero score
# mismatches, report candidate_scope and final_candidate_scope CONDITIONAL, and never claim global
# certification. Truth enters only the evaluation, after the call; the call's inputs are the model
# graph, the reads, the library and (for the synthetic mosaic) the curation manifest.
#
# Inputs are the frozen, gitignored packages under tests/results, bound here by the SHA-256 of their
# uncompressed content; a moved input refuses. Writes <work dir>/regressions.tsv.
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${1:?usage: genotype_production_c4_regressions.sh <build dir> <work dir> [case ...]}"
WORK="${2:?usage: genotype_production_c4_regressions.sh <build dir> <work dir> [case ...]}"
shift 2
CASES=("$@")
[ "${#CASES[@]}" -gt 0 ] || CASES=(na20509_10 na20509_30 na19240_30 na20509_lzo_30 bridgeable_30)
PANVAR="$BUILD/panvar"
PROBE="$BUILD/genotype_probe"
PROFILE="$REPO/profiles/genotype_c4.tsv"
EXPECTED="$REPO/tests/genotype_production_c4_expected"
PKG="$REPO/tests/results/genotype_loo_c4_development_7e37bb4/inputs"
FRZ="$REPO/tests/results/genotype_production_c4_inputs"
REF='HG00096#1#haplotype1-0000024:31848049-32074186'
[ -x "$PANVAR" ] && [ -x "$PROBE" ] || { echo "panvar and genotype_probe must be built in $BUILD" >&2; exit 1; }
if [ ! -f "$PKG/c4.gfa.gz" ] || [ ! -f "$FRZ/NA19240/c4.loo.gfa.gz" ]; then
  echo "skipped: the frozen inputs under tests/results are absent (gitignored)."; exit 0
fi
mkdir -p "$WORK"
sha() { if command -v sha256sum > /dev/null; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi; }
# stage <compressed or plain source> <destination> <expected sha256 of the content>
stage() {
  if [ ! -f "$2" ]; then
    case "$1" in *.gz) gzip -dc "$1" > "$2" ;; *) cp "$1" "$2" ;; esac
  fi
  local got; got="$(sha "$2")"
  [ "$got" = "$3" ] || { echo "REFUSED: $2 has sha256 $got, expected $3" >&2; exit 1; }
}
key() { awk -F'\t' -v k="$2" '$1==k{print $2; exit}' "$1"; }

# ---- the shared graph, and the cases' inputs by content
stage "$PKG/c4.gfa.gz" "$WORK/c4.gfa" b13af390c28dd8b43a913c43f123bb1fcf639e0569f99783fa12a778820c8733
case_inputs() {   # sets MODEL R1 R2 LIB CUR TRUTH DONOR SCHEME for case $1
  local d="$WORK/inputs/$1"; mkdir -p "$d"; CUR=""; TRUTH="$WORK/c4.gfa"
  case "$1" in
    na20509_10)
      DONOR=NA20509; SCHEME=loo
      stage "$PKG/NA20509/c4.loo.gfa.gz" "$d/model.gfa" 559c63f67ff34e1bcfbbacafbde8568e097382ee50a37ee41a9721b7355595b1
      stage "$PKG/NA20509/cov10/reads_1.fq.gz" "$d/r1.fq" c31ad6bd56569a0d8dc980e639f3d6ecca34452b466027f4b9a8de7b56beec6a
      stage "$PKG/NA20509/cov10/reads_2.fq.gz" "$d/r2.fq" 079711b3f64cd31a5b31bee15639a86e75dac1097bf7dea6eb20a73c25b46e75
      stage "$PKG/NA20509/cov10/library.txt" "$d/library.txt" 0c8938ee6521939f7e3fa1cbe05c431c4ef3232836edf4f730bc00c4f8f0a029 ;;
    na20509_30|na20509_lzo_30)
      DONOR=NA20509; SCHEME=loo
      if [ "$1" = na20509_lzo_30 ]; then
        SCHEME=lzo; cp "$WORK/c4.gfa" "$d/model.gfa"
      else
        stage "$PKG/NA20509/c4.loo.gfa.gz" "$d/model.gfa" 559c63f67ff34e1bcfbbacafbde8568e097382ee50a37ee41a9721b7355595b1
      fi
      stage "$PKG/NA20509/cov30/reads_1.fq.gz" "$d/r1.fq" 2e55d2910b3655bb0c22a6f285038b4005195459025faf364c58475897fb99de
      stage "$PKG/NA20509/cov30/reads_2.fq.gz" "$d/r2.fq" 81c6a2bcb0da92898769b4290c18353d0eae4361ac697272bc0cad9045d3daa9
      stage "$PKG/NA20509/cov30/library.txt" "$d/library.txt" 5fe66a66e99388a15f47c01907cb204710b809d1a3c305d41a4f4c9cf09f2a32 ;;
    na19240_30)
      DONOR=NA19240; SCHEME=loo
      stage "$FRZ/NA19240/c4.loo.gfa.gz" "$d/model.gfa" 05f9eb4b632369e7145d4fe741efbb0bae7c28df5727aa569bde5fa018e455b0
      stage "$FRZ/NA19240/reads_1.fq.gz" "$d/r1.fq" 64b908d4d0e4c401e1d6a78b85c146037e9bd21326fbab5fc796a1d3b1398e9f
      stage "$FRZ/NA19240/reads_2.fq.gz" "$d/r2.fq" d2cff2441f02fc40743f35c8a3ece8ea5fead05fd025d1effdd8626b23adb937
      stage "$FRZ/NA19240/library.txt" "$d/library.txt" 5fe66a66e99388a15f47c01907cb204710b809d1a3c305d41a4f4c9cf09f2a32 ;;
    bridgeable_30)
      DONOR=MOSAIC; SCHEME=mosaic
      stage "$FRZ/bridgeable/staged.gfa.gz" "$d/model.gfa" 5914a6a6b828f97297d24924e1e24c0fab33635edba4a8a77b3cde03d4f277e2
      stage "$FRZ/bridgeable/reads_1.fq.gz" "$d/r1.fq" 7070f18e8cbfb741a8978a86f6b0818500e2d7a026c15573f051a1bbe064506f
      stage "$FRZ/bridgeable/reads_2.fq.gz" "$d/r2.fq" 42b698ab27fe1947a34650ad1b39cf6ab0b1c33bf7aa3ebf0cb18217d4b444c8
      stage "$FRZ/bridgeable/library.txt" "$d/library.txt" 43a1b4f1dd2c2f870b48cc26347fafae93287889d8a6f54d79759d8db35f2859
      stage "$FRZ/bridgeable/mosaic_exclude.tsv" "$d/mosaic_exclude.tsv" 81eed2d6643caaaf6cf9d735c1cb215c371ff42ae3887d8142690bc34e3cb175
      CUR="$d/mosaic_exclude.tsv"; TRUTH="$d/model.gfa" ;;
    *) echo "unknown case $1" >&2; exit 1 ;;
  esac
  MODEL="$d/model.gfa"; R1="$d/r1.fq"; R2="$d/r2.fq"; LIB="$d/library.txt"
}
# expected <case>: score distance floor status phase anchor
expected() {
  case "$1" in
    na20509_10)     echo "-151090.979646 114 91 COMPLETE PRIOR_ONLY 329" ;;
    na20509_30)     echo "-414918.317265 113 91 COMPLETE PRIOR_ONLY 280" ;;
    na19240_30)     echo "-426151.472546 4 4 COMPLETE PRIOR_ONLY 7793" ;;
    na20509_lzo_30) echo "-409264.331248 0 0 COMPLETE PRIOR_ONLY 8461" ;;
    bridgeable_30)  echo "-405132.712785 0 0 COMPLETE PRIOR_ONLY 5356" ;;
  esac
}
SEARCH_COLS="round start_total moved repaired_block local_pair_rank assignment allele_x allele_y hap1_alleles hap2_alleles fragment_term exposure_term prior total delta_total status unique_best separation neighbours_generated neighbours_scored duplicate suppressed_visited self_skipped rejected_infinite_prior scored_pairs order_invariant swap_invariant origins_reused origins_recomputed fragments_reduced fragments_copied"
pick() {
  awk -F'\t' -v want="$SEARCH_COLS" 'BEGIN{n=split(want,w," ")} NR==1{for(i=1;i<=NF;i++)h[$i]=i}
    {line=""; for(j=1;j<=n;j++){line=line (j>1?"\t":"") $h[w[j]]} print line}' "$1"
}
if [ "$(uname)" = "Darwin" ]; then TIME=(/usr/bin/time -l); else TIME=(/usr/bin/time -v); fi

fail=0
check() { if [ "$2" = "$3" ]; then echo "  ok: $1 = $2"; else echo "  REGRESSION FAILED: $1 is '$2', expected '$3'"; fail=1; fi; }
TABLE="$WORK/regressions.tsv"
printf "case\tstatus\tanchor\trounds\tmoves\tstop\trepaired_blocks\tscore\tfinal_decision\ttop_set\tphase_support\tcandidate_scope\tfinal_candidate_scope\tgenotype_separation_gap\tsame_genotype_fragment_identical\tsame_genotype_prior_gap\twinner_source\twinner_candidate\tdistance\tfloor\texcess\tblocks_exact\tblocks\tpool_identities\tpool_repair_only\tpool_score_mismatches\ttrajectory\twall_seconds\tseconds_union_scoring\tseconds_repair\tpeak_rss_bytes\tpeak_phys_footprint_bytes\n" > "$TABLE"
for c in "${CASES[@]}"; do
  echo "== $c =="
  case_inputs "$c"
  o="$WORK/$c/call"; mkdir -p "$WORK/$c"; rm -f "$o".*
  cmd=("$PANVAR" genotype --profile "$PROFILE" --gfa "$MODEL" --reference-path "$REF" --reads1 "$R1" --reads2 "$R2"
       --library "$LIB" --out-prefix "$o" --quiet)
  [ -n "$CUR" ] && cmd+=(--panel-curation "$CUR")
  # Truth is never an input to the call: the command line names no truth file.
  case " ${cmd[*]} " in *" $TRUTH "*) [ "$TRUTH" = "$MODEL" ] || { echo "REFUSED: truth on the call's command line" >&2; exit 1; } ;; esac
  printf '%q ' "${cmd[@]}" > "$WORK/$c/command.txt"
  "${TIME[@]}" "${cmd[@]}" > "$WORK/$c/log" 2> "$WORK/$c/time"
  code=$?
  st="$o.mosaic_status.tsv"
  # ---- evaluation: truth enters here, after the call, and only here.
  ev=("$PROBE" evaluate --scheme "$SCHEME" --donor "$DONOR" --reference "$REF" --run-prefix "$o" --model-gfa "$MODEL"
      --truth-gfa "$TRUTH" --out "$WORK/$c/eval")
  [ -n "$CUR" ] && ev+=(--panel-curation "$CUR")
  "${ev[@]}" > "$WORK/$c/eval.log" 2>&1 || echo "  evaluation failed (see $WORK/$c/eval.log)"
  es="$WORK/$c/eval.summary.tsv"
  read -r want_score want_distance want_floor want_status want_phase want_anchor <<< "$(expected "$c")"
  score="$(awk -F'\t' 'NR==2{print $4}' "$o.mosaic_scores.tsv" 2>/dev/null)"
  distance="$(key "$es" unordered_distance_total)"; floor="$(key "$es" availability_floor_total)"
  winner_rep="$(awk -F'\t' 'NR==2{print $2}' "$o.mosaic_scores.tsv" 2>/dev/null)"
  winner_source="$(awk -F'\t' -v id="$winner_rep" 'NR>1&&$1==id{print $14}' "$o.mosaic_candidates.tsv" 2>/dev/null)"
  blocks="$(awk 'END{print NR-1}' "$o.block_genotypes.tsv" 2>/dev/null)"
  repaired="$(awk -F'\t' 'NR>1&&$3=="true"{printf "%s%s", (n++?",":""), $4}' "$o.repair_rounds.tsv" 2>/dev/null)"
  if diff -q <(pick "$EXPECTED/$c.rounds.tsv") <(pick "$o.repair_rounds.tsv") > /dev/null 2>&1; then trajectory=identical; else trajectory=DIFFERS; fi
  wall="$(awk '/ real /{print $1}' "$WORK/$c/time")"
  rss="$(awk '/maximum resident set size/{print $1}' "$WORK/$c/time")"
  footprint="$(awk '/peak memory footprint/{print $1}' "$WORK/$c/time")"
  check "exit" "$code" "0"
  check "status" "$(key "$st" status)" "$want_status"
  check "anchor (best exactly scored panel diplotype)" "$(key "$st" repair_anchor_candidate_id)" "$want_anchor"
  check "trajectory vs the recorded panel-anchor branch" "$trajectory" "identical"
  check "converged" "$(key "$st" repair_converged)" "true"
  check "final decision" "$(key "$st" final_decision)" "UNIQUE"
  check "score" "$score" "$want_score"
  check "distance" "$distance" "$want_distance"
  check "floor" "$floor" "$want_floor"
  check "phase support" "$(key "$st" phase_support)" "$want_phase"
  check "candidate_scope" "$(key "$st" candidate_scope)" "CONDITIONAL"
  check "final_candidate_scope" "$(key "$st" final_candidate_scope)" "CONDITIONAL"
  check "global_certified" "$(key "$st" global_certified)" "false"
  check "pool score mismatches" "$(key "$st" pool_score_mismatches)" "0"
  check "haplotype representation" "$(key "$st" haplotype_representation)" "prior_selected"
  check "FASTA header carries PRIOR_ONLY" "$(grep -c 'phase_support=PRIOR_ONLY representation=prior_selected' "$o.mosaic_haplotypes.fa" 2>/dev/null)" "2"
  case "$c" in
    na20509_30)
      check "excess over floor" "$((distance - floor))" "22"
      check "gap to the best different genotype (2 dp)" "$(printf '%.2f' "$(key "$st" genotype_separation_gap)")" "4.40"
      check "best same-genotype alternative: identical fragment term" "$(key "$st" same_genotype_alternative_fragment_identical)" "true"
      check "best same-genotype alternative: prior gap (2 dp)" "$(printf '%.2f' "$(key "$st" same_genotype_alternative_prior_gap)")" "7.66"
      check "not genotype-equivalent" "$(key "$st" status)" "COMPLETE" ;;
    na19240_30)
      check "ctx:6 repaired" "$(echo ",$repaired," | grep -c ',ctx:6,')" "1" ;;
    na20509_lzo_30)
      check "repair moves" "$(key "$st" repair_moves)" "0" ;;
    bridgeable_30)
      check "blocks exact" "$(key "$es" blocks_exact)/$blocks" "11/11"
      check "winner is the union's mosaic candidate" "$winner_source" "mosaic" ;;
  esac
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$c" "$(key "$st" status)" "$(key "$st" repair_anchor_candidate_id)" "$(key "$st" repair_rounds_run)" \
    "$(key "$st" repair_moves)" "$(key "$st" repair_stop_reason)" "${repaired:-none}" "$score" "$(key "$st" final_decision)" \
    "$(key "$st" top_set_size)" "$(key "$st" phase_support)" "$(key "$st" candidate_scope)" "$(key "$st" final_candidate_scope)" \
    "$(key "$st" genotype_separation_gap)" "$(key "$st" same_genotype_alternative_fragment_identical)" \
    "$(key "$st" same_genotype_alternative_prior_gap)" "$winner_source" "$winner_rep" "$distance" "$floor" \
    "$((distance - floor))" "$(key "$es" blocks_exact)" "$blocks" "$(key "$st" pool_identities)" "$(key "$st" pool_repair_only)" \
    "$(key "$st" pool_score_mismatches)" "$trajectory" "$wall" "$(key "$st" seconds_rescoring)" "$(key "$st" seconds_repair)" \
    "$rss" "$footprint" >> "$TABLE"
done
column -t -s$'\t' "$TABLE"
[ "$fail" -eq 0 ] && echo "production regressions: PASS" || { echo "production regressions: FAILED"; exit 1; }
