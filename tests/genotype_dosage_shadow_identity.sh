#!/usr/bin/env bash
# Shadow byte-identity gate (docs/algorithms/genotype.md, stage 11): one C4 campaign cell re-run with the
# dosage flag absent and with --dosage-audit shadow. Every pre-existing output must equal the frozen 11fe681 campaign
# output byte for byte, except the status file's wall-time, memory and binary-provenance lines; shadow may add only
# <prefix>.dosage_*.tsv. Usage: genotype_dosage_shadow_identity.sh <panvar binary> <campaign root> <cell> <scratch dir>
set -uo pipefail
PANVAR=$1; C=$2; CELL=$3; OUT=$4
mkdir -p "$OUT"
IFS=$'\t' read -r job camp donor cov scheme lab model r1 r2 lib cur truth x1 x2 x3 x4 salt prefix <<< "$(awk -F'\t' -v j="$CELL" '$1==j' "$C/jobs.tsv")"
[ -n "${prefix:-}" ] || { echo "no job $CELL"; exit 2; }
PROFILE=$(grep -o -- '--profile [^ ]*' "$prefix.command.txt" | cut -d' ' -f2)
args=(--profile "$PROFILE" --threads 8 --gfa "$model" --reference-path 'HG00096#1#haplotype1-0000024:31848049-32074186'
      --reads1 "$r1" --reads2 "$r2" --library "$lib" --quiet)
[ "$cur" != none ] && args+=(--panel-curation "$cur")
"$PANVAR" genotype "${args[@]}" --out-prefix "$OUT/off" > "$OUT/off.log" 2>&1; rc_off=$?
"$PANVAR" genotype "${args[@]}" --dosage-audit shadow --out-prefix "$OUT/shadow" > "$OUT/shadow.log" 2>&1; rc_shadow=$?
fail=0; n=0
[ $rc_off -eq $rc_shadow ] || { echo "FAIL exit status off=$rc_off shadow=$rc_shadow"; fail=1; }
# Normalisation, as the tie-branch gate defines it: the status apart from wall time, process memory, worker gauges,
# the scheduling-dependent emission-cache split and path-bearing provenance lines; repair rounds apart from their
# wall and memory columns (1-28 kept); the FASTA with the output path in its header replaced; the memory ledger
# (repair_memory.tsv, tie_branch_memory.tsv) is process memory only and is not compared; tie branches are compared
# apart from their wall-time column (1-14 kept).
status_science() {
  grep -vE '^(seconds_|resources_|threads_|phys_footprint|rss_|emission_cache_dp_hits|emission_cache_dp_misses|emission_cache_executed_dps|emission_cache_executed_prefilters|scorer_emission_cache_arena_bytes|union_cached_prefilters_consumed|union_prefilter_cache_misses|worker_seconds_|sha256_|profile|preset|settings_overridden)' "$1" | grep -v "$2"
}
norm() {   # norm <file> <prefix>
  case "$1" in
    *.mosaic_status.tsv) status_science "$1" "$2" ;;
    *.repair_rounds.tsv) cut -f1-28 "$1" ;;
    *.tie_branches.tsv) cut -f1-14 "$1" ;;
    *.mosaic_haplotypes.fa) sed "s#^>$2#>PREFIX#" "$1" ;;
    *) cat "$1" ;;
  esac
}
for f in "$prefix".*; do
  suf=${f#"$prefix"}
  case "$suf" in .eval.*|.evaluate.log|.command.txt|.log|.time|.done|.repair_memory.tsv|.tie_branch_memory.tsv) continue;; esac
  a=$(norm "$f" "$prefix" | shasum -a 256 | cut -c1-16)
  for arm in off shadow; do
    b=$(norm "$OUT/$arm$suf" "$OUT/$arm" 2>/dev/null | shasum -a 256 | cut -c1-16)
    [ "$a" = "$b" ] || { echo "FAIL $arm$suf differs from the campaign"; fail=1; }
  done
  n=$(( n + 1 ))
done
echo "compared $n campaign outputs against both arms"
extra=$(cd "$OUT" && ls shadow.* | grep -v '^shadow\.log$' | while read -r x; do [ -e "off${x#shadow}" ] || echo "$x"; done | grep -v '\.dosage_' || true)
[ -z "$extra" ] || { echo "FAIL shadow wrote non-dosage files: $extra"; fail=1; }
missing=$(cd "$OUT" && ls off.* | while read -r x; do [ -e "shadow${x#off}" ] || echo "$x"; done || true)
[ -z "$missing" ] || { echo "FAIL shadow lacks: $missing"; fail=1; }
echo "dosage files: $(cd "$OUT" && ls shadow.dosage_* | tr '\n' ' ')"
grep -E "dosage_audit_status|dosage_m|dosage_status_" "$OUT/shadow.dosage_status.tsv" | tr '\n' ' '; echo
[ $fail -eq 0 ] && echo "PASS $CELL" || echo "FAILED $CELL"
exit $fail
