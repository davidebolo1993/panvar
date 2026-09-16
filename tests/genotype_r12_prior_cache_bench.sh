#!/usr/bin/env bash
# Developer harness, not a regression test: does caching the exact-path prior remove the fixed cost R11
# measured, without changing a single statistical output?
#
#   BUILD=<build dir> tests/genotype_r12_prior_cache_bench.sh <out-dir>
#
# Reads are NOT simulated. Every case reuses reads already on disk -- the frozen 2x pilot reads and the
# compressed c4 source -- so the only difference from the R11 measurements is the code. Nothing is
# written into an existing results directory; the R11 pilot outputs are the "before" side of the identity
# comparison and must survive untouched.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${BUILD:-$REPO/build}"
OUT="${1:?usage: BUILD=<build> tests/genotype_r12_prior_cache_bench.sh <out-dir>}"
PROBE="$BUILD/genotype_probe"
PILOT="$REPO/tests/results"
C4_SOURCE="$PILOT/genotype_coverage_c4_held_out_abb35a0/HG00096/source"
# Workers for the diploid combination phase. One invocation uses one setting, so a 1,2,2,1 ordering is
# four invocations into four directories -- which exposes machine drift instead of hiding it inside a
# single run, and the drift is real: this machine has already moved ~1.95x between sessions.
THREADS="${THREADS:-1}"

mkdir -p "$OUT"
[ -z "$(ls -A "$OUT" 2>/dev/null)" ] || { echo "$OUT is not empty" >&2; exit 1; }

{
  printf 'key\tvalue\n'
  printf 'git_commit\t%s\n' "$(git -C "$REPO" rev-parse HEAD)"
  printf 'git_tree_dirty\t%s\n' "$(git -C "$REPO" status --porcelain | grep -q . && echo true || echo false)"
  printf 'probe_sha256\t%s\n' "$(shasum -a 256 "$PROBE" | cut -d' ' -f1)"
  printf 'script_sha256\t%s\n' "$(shasum -a 256 "${BASH_SOURCE[0]}" | cut -d' ' -f1)"
  printf 'threads\t%s\n' "$THREADS"
} > "$OUT/signature.tsv"

# c4 reads live compressed after the R11 disk cleanup; decompress once, into our own directory.
mkdir -p "$OUT/c4_reads"
gzip -dc "$C4_SOURCE/reads_1.fq.gz" > "$OUT/c4_reads/reads_1.fq"
gzip -dc "$C4_SOURCE/reads_2.fq.gz" > "$OUT/c4_reads/reads_2.fq"
cp "$C4_SOURCE/library.txt" "$OUT/c4_reads/library.txt"
gzip -dc "$REPO/tests/real_data/c4.gfa.gz" > "$OUT/c4.gfa"

SCRATCH="${SCRATCH:-/private/tmp/claude-502/-Users-davide-bolognini-panvar/bf7c1521-d4d1-4ddc-ac70-2206578e62d9/scratchpad/preflight}"
C4_REF="HG00096#1#haplotype1-0000024:31848049-32074186"

# case | gfa | reference | reads dir | prefixes | audit
CASES=(
  "c4_audit0|$OUT/c4.gfa|$C4_REF|$OUT/c4_reads|0|0"
  "c4_audit256|$OUT/c4.gfa|$C4_REF|$OUT/c4_reads|0,8|256"
  "acot_2x|$SCRATCH/acot.gfa|GRCh38#0#chr14:73416151-73668669|$PILOT/genotype_pilot_acot_2x_d24e471/source|0|0"
  "gstm1_2x|$SCRATCH/gstm1.gfa|GRCh38#0#chr1:109561616-109802331|$PILOT/genotype_pilot_gstm1_2x_d24e471/source|0|0"
  "ankrd36c_2x|$SCRATCH/ankrd36c.gfa|GRCh38#0#chr2:95666662-96152714|$PILOT/genotype_pilot_ankrd36c_2x_d24e471/source|0|0"
  # Same inputs and same command as the first case, run last. Any difference between the two is this
  # machine, not this code -- which is what a total moving by 2x on work the cache barely touches demands
  # before any timing change is attributed to the change.
  "c4_audit0_repeat|$OUT/c4.gfa|$C4_REF|$OUT/c4_reads|0|0"
)

# ONLY="name name" restricts the run to named cases, for a bounded review batch that must not spend the
# machine on cases the reviewer did not ask for.
if [ -n "${ONLY:-}" ]; then
  SELECTED=()
  for spec in "${CASES[@]}"; do
    IFS='|' read -r name _ <<< "$spec"
    case " $ONLY " in *" $name "*) SELECTED+=("$spec") ;; esac
  done
  [ "${#SELECTED[@]}" -gt 0 ] || { echo "ONLY matched no case" >&2; exit 1; }
  CASES=("${SELECTED[@]}")
  echo "restricted to: $ONLY"
fi

echo "=== runs ==="
for spec in "${CASES[@]}"; do
  IFS='|' read -r name gfa ref reads prefixes audit <<< "$spec"
  mkdir -p "$OUT/$name"
  /usr/bin/time -l "$PROBE" genotype --gfa "$gfa" --reference "$ref" \
    --reads1 "$reads/reads_1.fq" --reads2 "$reads/reads_2.fq" --library "$reads/library.txt" \
    --out-prefix "$OUT/$name/run" --state-limit 1024 --prefixes "$prefixes" --phase-audit "$audit" \
    --threads "$THREADS" \
    > "$OUT/$name/run.log" 2> "$OUT/$name/time.txt"
  echo "  $name done"
done

field() { awk -F'\t' -v k="$2" '$1==k{print $2}' "$1" 2>/dev/null; }
rss_mb() { awk '/maximum resident set size/{printf "%.0f", $1/1048576}' "$1" 2>/dev/null; }

echo
echo "=== benchmark table ==="
printf '%-12s %8s %8s %10s %10s %9s %9s %9s %9s %8s\n' \
  case pairs paths hap_evals dip_lookups s_prior s_group s_rank s_total rss_MB
for spec in "${CASES[@]}"; do
  IFS='|' read -r name _ _ _ prefixes _ <<< "$spec"
  k="${prefixes%%,*}"
  s="$OUT/$name/run.K$k.mosaic_status.tsv"
  printf '%-12s %8s %8s %10s %10s %9s %9s %9s %9s %8s\n' \
    "$name" "$(field "$s" panel_diplotypes)" "$(field "$s" distinct_haplotype_paths)" \
    "$(field "$s" haploid_prior_evaluations)" "$(field "$s" diploid_prior_lookups)" \
    "$(field "$s" seconds_prior)" "$(field "$s" seconds_candidate_grouping)" "$(field "$s" seconds_ranking)" \
    "$(field "$s" seconds_total)" "$(rss_mb "$OUT/$name/time.txt")"
done

echo
echo "=== prior-evaluation invariant: haploid evaluations must equal distinct paths ==="
bad=0
for spec in "${CASES[@]}"; do
  IFS='|' read -r name _ _ _ prefixes _ <<< "$spec"
  k="${prefixes%%,*}"
  s="$OUT/$name/run.K$k.mosaic_status.tsv"
  p="$(field "$s" distinct_haplotype_paths)"; e="$(field "$s" haploid_prior_evaluations)"
  if [ "$p" = "$e" ]; then echo "  OK   $name: $e evaluations for $p distinct paths"
  else echo "  FAIL $name: $e evaluations for $p distinct paths"; bad=$((bad+1)); fi
done

# ---- identity against the R11 measurements, which used the same reads and the unoptimised code -------
# Expected to differ: timing fields, the three new counters, and the per-block candidate_scope column,
# which R11 reported wrongly. Everything statistical must be byte-identical.
drop_status() { grep -vE '^(seconds_|distinct_haplotype_paths|haploid_prior_evaluations|diploid_prior_lookups|threads_|candidate_fragment_evaluations)' "$1" 2>/dev/null; }
drop_last()   { awk -F'\t' '{NF--; print}' OFS='\t' "$1" 2>/dev/null; }
seqs_only()   { grep -v '^>' "$1" 2>/dev/null; }

echo
echo "=== output identity vs the R11 pilots (same reads, unoptimised code) ==="
for L in acot gstm1 ankrd36c; do
  old="$PILOT/genotype_pilot_${L}_2x_d24e471/panel.K0"
  new="$OUT/${L}_2x/run.K0"
  # A restricted run compares only what it actually ran; silence here would look like a pass.
  [ -f "$new.genotypes.tsv" ] || { echo "  skipped    $L (not run in this batch)"; continue; }
  for what in genotypes scores status blocks haplotypes; do
    case "$what" in
      genotypes)  a="$old.genotypes.tsv";       b="$new.genotypes.tsv";       cmp_out=$(diff <(cat "$a") <(cat "$b") 2>&1 || true) ;;
      scores)     a="$old.mosaic_scores.tsv";   b="$new.mosaic_scores.tsv";   cmp_out=$(diff <(cat "$a") <(cat "$b") 2>&1 || true) ;;
      status)     a="$old.mosaic_status.tsv";   b="$new.mosaic_status.tsv";   cmp_out=$(diff <(drop_status "$a") <(drop_status "$b") 2>&1 || true) ;;
      blocks)     a="$old.block_genotypes.tsv"; b="$new.block_genotypes.tsv"; cmp_out=$(diff <(drop_last "$a") <(drop_last "$b") 2>&1 || true) ;;
      haplotypes) a="$old.mosaic_haplotypes.fa"; b="$new.mosaic_haplotypes.fa"; cmp_out=$(diff <(seqs_only "$a") <(seqs_only "$b") 2>&1 || true) ;;
    esac
    if [ -z "$cmp_out" ]; then echo "  identical  $L $what"
    else echo "  DIFFERS    $L $what"; echo "$cmp_out" | head -6 | sed 's/^/      /'; bad=$((bad+1)); fi
  done
done

echo
echo "=== scope consistency: every block row must equal its status ==="
for spec in "${CASES[@]}"; do
  IFS='|' read -r name _ _ _ prefixes _ <<< "$spec"
  k="${prefixes%%,*}"
  s="$(field "$OUT/$name/run.K$k.mosaic_status.tsv" candidate_scope)"
  rows="$(awk -F'\t' 'NR>1{print $NF}' "$OUT/$name/run.K$k.block_genotypes.tsv" 2>/dev/null | sort -u | tr '\n' ' ')"
  if [ "$(echo $rows)" = "$s" ]; then echo "  OK   $name: status $s, rows $rows"
  else echo "  FAIL $name: status $s, rows $rows"; bad=$((bad+1)); fi
done

echo
echo "=== rescoring split: which phase actually dominates ==="
printf '%-14s %10s %10s %10s %10s %9s\n' case rescoring prep place combine 'combine%'
for spec in "${CASES[@]}"; do
  IFS='|' read -r name _ _ _ prefixes _ <<< "$spec"
  k="${prefixes%%,*}"
  s="$OUT/$name/run.K$k.mosaic_status.tsv"
  awk -F'\t' -v c="$name" '
    $1=="seconds_rescoring"{r=$2}
    $1=="seconds_rescoring_fragment_preparation"{p=$2}
    $1=="seconds_rescoring_placement_and_mass"{m=$2}
    $1=="seconds_rescoring_diploid_combination"{d=$2}
    END{printf "%-14s %10.2f %10.2f %10.2f %10.2f %8.1f%%\n", c, r, p, m, d, (r>0 ? 100*d/r : 0)}' "$s"
done

echo
echo "failures: $bad"
[ "$bad" -eq 0 ] || exit 1
