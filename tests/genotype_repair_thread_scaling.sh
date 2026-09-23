#!/usr/bin/env bash
# Thread-scaling benchmark: the frozen NA20509 30x repair search at several thread counts.
#
#   tests/genotype_repair_thread_scaling.sh <build dir> <work dir> [thread counts...]
#   tests/genotype_repair_thread_scaling.sh --dry-run <build dir> <work dir> [thread counts...]
#   tests/genotype_repair_thread_scaling.sh --summarize <out prefix> <threads> <time file>
#
# Default thread counts: 1 2 4 8 16. 32 runs only if requested AND the node has at least 32 physical
# cores. Any count above the node's physical cores is skipped, not run: a benchmark on oversubscribed
# cores measures the scheduler, not the scorer.
#
# Inputs come only from the frozen, signature-bound package tests/results/genotype_loo_c4_development_7e37bb4
# (gitignored: copy it to the node). Nothing is resimulated. Each run is the configuration of the
# converged local benchmark: 384 MiB planned origin cache, up to 11 repair rounds, the declared caps.
#
# Reported per thread count: physical and logical cores, workers requested, workers actually active at
# peak (the process-wide gauge), enumeration time (sum of the per-round window-enumeration seconds),
# fixed time (wall minus enumeration), total wall, peak RSS, peak physical footprint (macOS only), and
# whether the trajectory is identical to the one-thread run and to the recorded converged result
# (-414918.317265, distance 113, floor 91).
set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PKG="$REPO/tests/results/genotype_loo_c4_development_7e37bb4"
REF='HG00096#1#haplotype1-0000024:31848049-32074186'
EXPECTED_TOTAL="-414918.317265"
EXPECTED_DISTANCE="113"
EXPECTED_FLOOR="91"

physical_cores() {
  if [ "$(uname)" = "Darwin" ]; then sysctl -n hw.physicalcpu
  elif command -v lscpu > /dev/null; then lscpu -p=Core,Socket | grep -v '^#' | sort -u | wc -l | tr -d ' '
  else getconf _NPROCESSORS_ONLN; fi
}
logical_cores() {
  if [ "$(uname)" = "Darwin" ]; then sysctl -n hw.logicalcpu; else getconf _NPROCESSORS_ONLN; fi
}
key() { awk -F'\t' -v k="$2" '$1==k{print $2}' "$1"; }

# One result row from a finished (or stopped) run.
summarize() {   # summarize <out prefix> <threads> <time file>
  local o="$1" t="$2" time_file="$3"
  local rounds="$o.rounds.tsv"; [ -f "$rounds" ] || rounds="$o.rounds.partial.tsv"
  local enumeration wall rss footprint status workers ewk total distance floor
  enumeration="$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}{s+=$h["seconds_enumeration_all"]}END{printf "%.3f",s}' "$rounds")"
  wall="$(key "$o.run_status.tsv" elapsed_seconds)"
  status="$(key "$o.run_status.tsv" status)"
  workers="$(key "$o.summary.tsv" peak_active_workers 2>/dev/null)"
  ewk="$(key "$o.summary.tsv" enumeration_workers 2>/dev/null)"
  total="$(key "$o.summary.tsv" final_total 2>/dev/null)"
  distance="$(key "$o.summary.tsv" final_distance 2>/dev/null)"
  floor="$(key "$o.summary.tsv" availability_floor_total 2>/dev/null)"
  if [ "$(uname)" = "Darwin" ]; then
    rss="$(awk '/maximum resident set size/{print $1}' "$time_file")"
    footprint="$(awk '/peak memory footprint/{print $1}' "$time_file")"
  else
    rss="$(awk -F': ' '/Maximum resident set size/{print $2 * 1024}' "$time_file")"
    footprint="NA"
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.3f\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$t" "$(physical_cores)" "$(logical_cores)" "${workers:-NA}" "${ewk:-NA}" "$status" "$wall" "$enumeration" \
    "$(awk -v w="$wall" -v e="$enumeration" 'BEGIN{print w-e}')" "${rss:-NA}" "${footprint:-NA}" \
    "${total:-NA}" "${distance:-NA}" "${floor:-NA}" \
    "$([ "${total:-}" = "$EXPECTED_TOTAL" ] && [ "${distance:-}" = "$EXPECTED_DISTANCE" ] && [ "${floor:-}" = "$EXPECTED_FLOOR" ] && echo expected || echo DIFFERS)"
}
HEADER="threads\tphysical_cores\tlogical_cores\tpeak_active_workers\tenumeration_workers\tstatus\twall_seconds\tenumeration_seconds\tfixed_seconds\tpeak_rss_bytes\tpeak_phys_footprint_bytes\tfinal_total\tfinal_distance\tfloor\tresult_vs_recorded"

if [ "${1:-}" = "--summarize" ]; then
  printf "%b\n" "$HEADER"
  summarize "$2" "$3" "$4"
  exit 0
fi
DRY=0
if [ "${1:-}" = "--dry-run" ]; then DRY=1; shift; fi
BUILD="${1:?usage: genotype_repair_thread_scaling.sh [--dry-run] <build dir> <work dir> [threads...]}"
WORK="${2:?usage: genotype_repair_thread_scaling.sh [--dry-run] <build dir> <work dir> [threads...]}"
shift 2
THREADS=("$@")
[ "${#THREADS[@]}" -gt 0 ] || THREADS=(1 2 4 8 16)
PROBE="$BUILD/genotype_probe"
[ -x "$PROBE" ] || { echo "genotype_probe not built at $PROBE" >&2; exit 1; }
[ -f "$PKG/NA20509/call.30.K128.mosaic_status.tsv" ] || { echo "the frozen package is absent: $PKG" >&2; exit 1; }
PHYSICAL="$(physical_cores)"
echo "node: $(uname -s) $(uname -m), physical cores $PHYSICAL, logical cores $(logical_cores)"

mkdir -p "$WORK"
R1="$WORK/r1.fq"; R2="$WORK/r2.fq"; L="$WORK/c4.loo.gfa"; G="$WORK/c4.gfa"
if [ "$DRY" -eq 0 ]; then
  [ -f "$R1" ] || gzip -dc "$PKG/inputs/NA20509/cov30/reads_1.fq.gz" > "$R1"
  [ -f "$R2" ] || gzip -dc "$PKG/inputs/NA20509/cov30/reads_2.fq.gz" > "$R2"
  [ -f "$L" ] || gzip -dc "$PKG/inputs/NA20509/c4.loo.gfa.gz" > "$L"
  [ -f "$G" ] || gzip -dc "$PKG/inputs/c4.gfa.gz" > "$G"
  # The reads are bound by content in the campaign settings; a changed file refuses.
  sha() { if command -v sha256sum > /dev/null; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi; }
  want1="$(awk -F'\t' '$1=="source"&&$2=="NA20509"&&$3=="reads1_sha256"{print $4}' "$PKG/settings.tsv")"
  want2="$(awk -F'\t' '$1=="source"&&$2=="NA20509"&&$3=="reads2_sha256"{print $4}' "$PKG/settings.tsv")"
  [ "$(sha "$R1")" = "$want1" ] && [ "$(sha "$R2")" = "$want2" ] || { echo "frozen reads do not match settings.tsv" >&2; exit 1; }
fi
if [ "$(uname)" = "Darwin" ]; then TIME=(/usr/bin/time -l); else TIME=(/usr/bin/time -v); fi
[ -x /usr/bin/time ] || { echo "/usr/bin/time is required for peak memory (GNU time on Linux)" >&2; exit 1; }

TABLE="$WORK/thread_scaling.tsv"
[ "$DRY" -eq 1 ] || printf "%b\n" "$HEADER" > "$TABLE"
for t in "${THREADS[@]}"; do
  if [ "$t" -gt "$PHYSICAL" ]; then
    echo "skip: $t threads exceeds this node's $PHYSICAL physical cores"; continue
  fi
  o="$WORK/t$t/NA20509.30.K128"
  cmd=("${TIME[@]}" "$PROBE" repair-search --reference "$REF" --donor NA20509
       --model-gfa "$L" --run "$PKG/NA20509/call.30.K128" --reads1 "$R1" --reads2 "$R2"
       --library "$PKG/inputs/NA20509/cov30/library.txt"
       --threads "$t" --delta 1 --origin-cache 1 --origin-cache-bytes 402653184 --origin-cache-policy planned
       --max-physical-origins 4900000000 --max-state-classes 4200000000
       --max-piece-hits 14000000000 --max-latent-pairs 14000000000
       --truth-gfa "$G" --local-pairs 32 --rounds 11 --max-neighbors 670
       --time-limit-seconds 1200 --out "$o")
  if [ "$DRY" -eq 1 ]; then echo "would run: ${cmd[*]}"; continue; fi
  mkdir -p "$WORK/t$t"
  echo "running $t thread(s)..."
  "${cmd[@]}" > "$WORK/t$t/out.log" 2> "$WORK/t$t/time"
  summarize "$o" "$t" "$WORK/t$t/time" >> "$TABLE"
done
[ "$DRY" -eq 1 ] && exit 0

# The trajectory must not depend on the thread count: every run's scientific columns against the first.
first="${THREADS[0]}"
for t in "${THREADS[@]}"; do
  [ -f "$WORK/t$t/NA20509.30.K128.rounds.tsv" ] || continue
  if diff -q <(cut -f1-28 "$WORK/t$first/NA20509.30.K128.rounds.tsv") <(cut -f1-28 "$WORK/t$t/NA20509.30.K128.rounds.tsv") > /dev/null; then
    echo "trajectory at $t threads: identical to $first thread(s)"
  else
    echo "trajectory at $t threads: DIFFERS from $first thread(s)"
  fi
done
column -t -s$'\t' "$TABLE"
