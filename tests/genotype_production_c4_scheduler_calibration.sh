#!/usr/bin/env bash
# The scheduler calibration for the final C4 campaign. It ran once, on 2026-09-23 under the v1 campaign profile, and
# selected JOBS=1 THREADS=8; it now calls the campaign's production profile (v6 campaign-local), so a new calibration
# measures what the campaign runs. It is not rerun by the campaign; a plan with more than one concurrent job requires
# its calibration.tsv (tests/genotype_production_c4_campaign.sh plan, SCHEDULER_CALIBRATION). Two representative 10x jobs
# on frozen inputs -- NA20509 LOO and NA20509 ordinary LZO -- first one after the other at 8 threads each
# (JOBS=1, THREADS=8), then both at once at 4 threads each (JOBS=2, THREADS=4). No other combination is run.
#
#   tests/genotype_production_c4_scheduler_calibration.sh <build dir> <work dir> [reference LOO 10x prefix]
#
# Selects JOBS=2 THREADS=4 only if all hold: two-job throughput at least 30% higher; the combined physical
# footprint below 75% of physical memory; neither process reaching the profile's footprint limit. Otherwise JOBS=1
# THREADS=8. Scientific outputs must be byte-identical between the configurations, and the LOO job to the
# given reference (an NA20509 LOO 10x call made under the same profile; the v1 production regression's call is not
# one). Writes <work>/calibration.tsv.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${1:?usage: genotype_production_c4_scheduler_calibration.sh <build dir> <work dir> [reference prefix]}"
WORK="${2:?usage: genotype_production_c4_scheduler_calibration.sh <build dir> <work dir> [reference prefix]}"
REFERENCE_LOO="${3:-}"
PROFILE="$REPO/profiles/genotype_c4_v6_campaign_local.tsv"
PKG="$REPO/tests/results/genotype_loo_c4_development_7e37bb4/inputs"
REF='HG00096#1#haplotype1-0000024:31848049-32074186'
mkdir -p "$WORK"
sha() { shasum -a 256 "$1" | cut -d' ' -f1; }
stage() { [ -f "$2" ] || gzip -dc "$1" > "$2"; [ "$(sha "$2")" = "$3" ] || { echo "REFUSED: $2 moved" >&2; exit 1; }; }
key() { awk -F'\t' -v k="$2" '$1==k{print $2; exit}' "$1"; }
stage "$PKG/c4.gfa.gz" "$WORK/c4.gfa" b13af390c28dd8b43a913c43f123bb1fcf639e0569f99783fa12a778820c8733
stage "$PKG/NA20509/c4.loo.gfa.gz" "$WORK/c4.loo.gfa" 559c63f67ff34e1bcfbbacafbde8568e097382ee50a37ee41a9721b7355595b1
stage "$PKG/NA20509/cov10/reads_1.fq.gz" "$WORK/r1.fq" c31ad6bd56569a0d8dc980e639f3d6ecca34452b466027f4b9a8de7b56beec6a
stage "$PKG/NA20509/cov10/reads_2.fq.gz" "$WORK/r2.fq" 079711b3f64cd31a5b31bee15639a86e75dac1097bf7dea6eb20a73c25b46e75
LIB="$PKG/NA20509/cov10/library.txt"
call() {   # call <name> <model> <threads>
  "/usr/bin/time" -l "$BUILD/panvar" genotype --profile "$PROFILE" --threads "$3" --gfa "$2" --reference-path "$REF" \
    --reads1 "$WORK/r1.fq" --reads2 "$WORK/r2.fq" --library "$LIB" --out-prefix "$WORK/$1" --quiet \
    > "$WORK/$1.log" 2> "$WORK/$1.time"
  echo $? > "$WORK/$1.exit"
}
science() {   # the scientific outputs of one call, one stream
  local p="$1"
  cat "$p.genotypes.tsv" "$p.block_genotypes.tsv" "$p.mosaic_scores.tsv" "$p.mosaic_candidates.tsv" \
      "$p.equivalence_set.tsv" 2>/dev/null
  grep -v '^>' "$p.mosaic_haplotypes.fa" 2>/dev/null
  cut -f1-28 "$p.repair_rounds.tsv" 2>/dev/null
}
now() { python3 -c 'import time; print(f"{time.time():.3f}")'; }

echo "== A: JOBS=1 THREADS=8, one job after the other =="
t0="$(now)"; call loo_t8 "$WORK/c4.loo.gfa" 8; call lzo_t8 "$WORK/c4.gfa" 8; t1="$(now)"
echo "== B: JOBS=2 THREADS=4, both at once =="
t2="$(now)"; call loo_t4 "$WORK/c4.loo.gfa" 4 & a=$!; call lzo_t4 "$WORK/c4.gfa" 4 & b=$!
peak_rss=0
while kill -0 "$a" 2>/dev/null || kill -0 "$b" 2>/dev/null; do
  now_rss="$(ps -A -o rss=,command= | awk '/panvar genotype --profile/ && !/awk/ {s+=$1} END{print s*1024+0}')"
  [ "$now_rss" -gt "$peak_rss" ] && peak_rss="$now_rss"
  sleep 1
done
wait "$a"; wait "$b"; t3="$(now)"

memory="$(sysctl -n hw.memsize)"; cores="$(sysctl -n hw.logicalcpu)"
row() {   # row <name> <config>
  local n="$1" t="$WORK/$1.time" s="$WORK/$1.mosaic_status.tsv"
  local wall user sys fp rss
  wall="$(awk '/ real /{print $1}' "$t")"; user="$(awk '/ real /{print $3}' "$t")"; sys="$(awk '/ real /{print $5}' "$t")"
  fp="$(awk '/peak memory footprint/{print $1}' "$t")"; rss="$(awk '/maximum resident set size/{print $1}' "$t")"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%.2f\t%s\t%s\t%s\t%s\t%s\n" "$n" "$2" "$(cat "$WORK/$n.exit")" "$(key "$s" status)" "$wall" \
    "$(awk -v u="$user" -v s="$sys" 'BEGIN{print u+s}')" "$(awk -v u="$user" -v s="$sys" -v w="$wall" 'BEGIN{print (u+s)/w}')" \
    "$(key "$s" threads_peak_active_workers)" "$fp" "$rss" "$(key "$s" settings_overridden_after_profile)" "$(key "$s" call_class)"
}
{
  printf "job\tconfig\texit\tstatus\twall_seconds\tcpu_seconds\tcpu_per_wall\tpeak_active_workers\tpeak_footprint_bytes\tpeak_rss_bytes\toverrides\tcall_class\n"
  row loo_t8 "JOBS=1 THREADS=8"; row lzo_t8 "JOBS=1 THREADS=8"; row loo_t4 "JOBS=2 THREADS=4"; row lzo_t4 "JOBS=2 THREADS=4"
} > "$WORK/calibration_jobs.tsv"
column -t -s$'\t' "$WORK/calibration_jobs.tsv"

elapsed_a="$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f", b-a}')"
elapsed_b="$(awk -v a="$t2" -v b="$t3" 'BEGIN{printf "%.2f", b-a}')"
gain="$(awk -v a="$elapsed_a" -v b="$elapsed_b" 'BEGIN{printf "%.4f", a/b - 1}')"
fp_loo4="$(awk '/peak memory footprint/{print $1}' "$WORK/loo_t4.time")"; fp_lzo4="$(awk '/peak memory footprint/{print $1}' "$WORK/lzo_t4.time")"
combined_fp=$((fp_loo4 + fp_lzo4))   # an upper bound: the two peaks need not coincide
limit_hit=false
for n in loo_t8 lzo_t8 loo_t4 lzo_t4; do [ "$(key "$WORK/$n.mosaic_status.tsv" status)" = RESOURCE_REFUSAL ] && limit_hit=true; done
same() { diff -q <(science "$1") <(science "$2") > /dev/null && echo identical || echo DIFFERS; }
id_loo="$(same "$WORK/loo_t8" "$WORK/loo_t4")"; id_lzo="$(same "$WORK/lzo_t8" "$WORK/lzo_t4")"
id_ref="not_given"; [ -n "$REFERENCE_LOO" ] && id_ref="$(same "$REFERENCE_LOO" "$WORK/loo_t8")"
cpu_b="$(awk -F'\t' 'NR>1 && $2=="JOBS=2 THREADS=4"{s+=$6} END{print s}' "$WORK/calibration_jobs.tsv")"
cpu_a="$(awk -F'\t' 'NR>1 && $2=="JOBS=1 THREADS=8"{s+=$6} END{print s}' "$WORK/calibration_jobs.tsv")"
choice="JOBS=1 THREADS=8"
if awk -v g="$gain" 'BEGIN{exit !(g >= 0.30)}' && [ "$combined_fp" -lt "$((memory * 3 / 4))" ] && [ "$limit_hit" = false ]; then
  choice="JOBS=2 THREADS=4"
fi
[ "$id_loo" = identical ] && [ "$id_lzo" = identical ] && [ "$id_ref" != DIFFERS ] || choice="REFUSED: outputs differ"
{
  printf "key\tvalue\n"
  printf "elapsed_two_jobs_jobs1_threads8_seconds\t%s\nelapsed_two_jobs_jobs2_threads4_seconds\t%s\n" "$elapsed_a" "$elapsed_b"
  printf "throughput_gain\t%s\nthroughput_gain_required\t0.30\n" "$gain"
  printf "cpu_utilization_jobs1_threads8\t%s\n" "$(awk -v c="$cpu_a" -v e="$elapsed_a" -v n="$cores" 'BEGIN{printf "%.3f", c/(e*n)}')"
  printf "cpu_utilization_jobs2_threads4\t%s\n" "$(awk -v c="$cpu_b" -v e="$elapsed_b" -v n="$cores" 'BEGIN{printf "%.3f", c/(e*n)}')"
  printf "logical_cores\t%s\nphysical_memory_bytes\t%s\n" "$cores" "$memory"
  printf "combined_peak_footprint_jobs2_bytes\t%s\ncombined_footprint_ceiling_bytes\t%s\n" "$combined_fp" "$((memory * 3 / 4))"
  printf "combined_sampled_rss_peak_jobs2_bytes\t%s\nper_process_limit_reached\t%s\n" "$peak_rss" "$limit_hit"
  printf "total_active_workers_jobs2\t%s\n" "$(awk -F'\t' 'NR>1 && $2=="JOBS=2 THREADS=4"{s+=$8} END{print s}' "$WORK/calibration_jobs.tsv")"
  printf "identity_loo_threads4_vs_threads8\t%s\nidentity_lzo_threads4_vs_threads8\t%s\nidentity_loo_vs_regression\t%s\n" \
    "$id_loo" "$id_lzo" "$id_ref"
  printf "selected\t%s\n" "$choice"
} > "$WORK/calibration.tsv"
column -t -s$'\t' "$WORK/calibration.tsv"
case "$choice" in REFUSED*) exit 1 ;; esac
