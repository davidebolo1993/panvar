#!/usr/bin/env bash
# genotype_interval_checkpoint.sh - stop/resume equivalence for --interval-score.
#
#   genotype_interval_checkpoint.sh <panvar> <out-dir>
#
# WHY THIS EXISTS, concretely: a full-panel interval run was killed after 13.5 hours having produced
# NOTHING, because output was written only at the end. A multi-hour computation that loses
# everything on interruption is not usable, and "it was still running" is not a result.
#
# THE CONTRACT:
#   * checkpoints are written atomically (temp + rename), so a kill leaves either the previous
#     complete checkpoint or the new one, never a torn file;
#   * resume is REFUSED unless binary, reads, graph, candidate manifest and every parameter match --
#     splicing two different computations together is worse than starting over;
#   * an interrupted run reports INCOMPLETE and must never be interpreted as ambiguity;
#   * an uninterrupted run and a stop/resume run produce IDENTICAL intervals and verdicts.
set -uo pipefail
BIN="${1:?usage: genotype_interval_checkpoint.sh <panvar> <outdir>}"
OUT="${2:?}"; mkdir -p "$OUT"; OUT="$OUT/run.$$"; rm -rf "$OUT"; mkdir -p "$OUT"
fails=0
ok()  { printf '  ok   %s\n' "$*"; }
bad() { printf '  FAIL %s\n' "$*"; fails=$((fails+1)); }
seq_of() { awk -v n="$1" -v s="$2" 'BEGIN{x=s; b="ACGT";
           for(i=0;i<n;i++){x=(1103515245*x+12345)%2147483648; printf "%s", substr(b,(int(x/65536)%4)+1,1)} }'; }

U=$(seq_of 3000 21); V=$(seq_of 3000 22)
HA="${U}${V}"; HB="${U}$(seq_of 3000 23)"
{ printf 'H\tVN:Z:1.0\n'; printf 'S\t1\t%s\n' "$HA"; printf 'S\t2\t%s\n' "$HB"
  printf 'L\t1\t+\t2\t+\t0M\n'
  printf 'P\thapA\t1+\t*\n'; printf 'P\thapB\t2+\t*\n'; } > "$OUT/g.gfa"
"$BIN" bubble -i "$OUT/g.gfa" -r hapA -o "$OUT/b" --min-variant-bp 0 -q >/dev/null 2>&1
: > "$OUT/r1.fq"; : > "$OUT/r2.fq"
python3 - "$HA" "$OUT" <<'PY'
import sys
h,out=sys.argv[1],sys.argv[2]
comp={'A':'T','C':'G','G':'C','T':'A'}
with open(out+"/r1.fq","w") as f1, open(out+"/r2.fq","w") as f2:
    for i in range(0, 5800, 2):
        a=h[i:i+150]; b=h[i+200:i+350]
        if len(a)<150 or len(b)<150: break
        rc="".join(comp[c] for c in reversed(b))
        f1.write("@f%d/1\n%s\n+\n%s\n"%(i,a,"I"*150))
        f2.write("@f%d/2\n%s\n+\n%s\n"%(i,rc,"I"*150))
PY
N=$(( $(wc -l < "$OUT/r1.fq") / 4 ))
[ "$N" -gt 100 ] && ok "fixture has $N fragments" || bad "only $N fragments"

run() { "$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/b" -o "$OUT/$1" \
  -R "$OUT/r1.fq" -R "$OUT/r2.fq" -t 2 --max-divergence 0.05 --fragment-len 350 \
  --fragment-sd 50 --error-rate 0.001 --interval-score "$OUT/$1.tsv" --interval-tol 1.0 \
  --interval-candidates "hapA,hapB" "${@:2}" -q; }

run full >/dev/null 2>&1
[ -s "$OUT/full.tsv" ] && ok "uninterrupted run produced a result" || bad "no uninterrupted result"

# Interrupt mid-run, then resume from the checkpoint.
run part --interval-checkpoint "$OUT/ck.txt" --interval-batch 5 >/dev/null 2>&1 &
BP=$!; sleep 0.35; kill -TERM $BP 2>/dev/null; wait $BP 2>/dev/null
CKN=$(grep -vc '^#' "$OUT/ck.txt" 2>/dev/null || echo 0)
[ "${CKN:-0}" -gt 0 ] \
  && ok "an interrupted run left a checkpoint with $CKN fragments" \
  || bad "no checkpoint survived the interrupt"
run resume --interval-checkpoint "$OUT/ck.txt" --interval-batch 5 >/dev/null 2>&1
RS=$(grep -E "^# resumed" "$OUT/resume.tsv" | cut -f2)
UF=$(grep -E "^# unfinished" "$OUT/resume.tsv" | cut -f2)
# NON-VACUITY: the interrupt must land MID-run. Resuming 400 of 400 means the first run finished
# before the kill and nothing was actually interrupted -- measured, that is what the first version
# of this fixture did.
{ [ "${RS:-0}" -gt 0 ] && [ "${RS:-0}" -lt "$N" ]; } \
  && ok "the resumed run reused $RS of $N fragments -- a genuine partial" \
  || bad "resume reused ${RS:-0} of $N; the interrupt did not land mid-run"
[ "${UF:-1}" = 0 ] && ok "and finished every fragment" || bad "${UF:-?} fragments unfinished"

# THE LOAD-BEARING ASSERTION: identical results either way.
if diff <(grep -v '^#' "$OUT/full.tsv") <(grep -v '^#' "$OUT/resume.tsv") >/dev/null 2>&1; then
  ok "stop/resume gives IDENTICAL class intervals to an uninterrupted run"
else
  bad "stop/resume intervals differ from the uninterrupted run"
fi
for k in verdict robust_gap n_plausible; do
  a=$(grep "^# $k" "$OUT/full.tsv" | cut -f2); b=$(grep "^# $k" "$OUT/resume.tsv" | cut -f2)
  [ "$a" = "$b" ] && ok "$k identical ($a)" || bad "$k differs: $a vs $b"
done

# RESUME MUST REFUSE a mismatched signature rather than splice two computations.
M=$("$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/b" -o "$OUT/mis" -R "$OUT/r1.fq" -R "$OUT/r2.fq" \
    -t 2 --max-divergence 0.05 --fragment-len 350 --fragment-sd 50 --error-rate 0.001 \
    --interval-score "$OUT/mis.tsv" --interval-tol 0.25 --interval-candidates "hapA,hapB" \
    --interval-checkpoint "$OUT/ck.txt" 2>&1 | grep -ci "signature does not match")
[ "${M:-0}" -ge 1 ] \
  && ok "a checkpoint from different parameters is REFUSED, not spliced" \
  || bad "a mismatched checkpoint was accepted"

echo
if [ "$fails" -eq 0 ]; then echo "interval checkpoint: all assertions passed"; else
  echo "interval checkpoint: $fails assertion(s) failed"; fi
exit "$fails"
