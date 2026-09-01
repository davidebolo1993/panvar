#!/usr/bin/env bash
# genotype_frag_factorisation.sh - KNOWN FAILURE. BLOCKS PLAN STEP B.
#
#   genotype_frag_factorisation.sh <panvar-binary> <out-dir>
#
# NOT registered in CMakeLists, deliberately. It does not pass, and it must not be dressed up as an
# expected failure: WILL_FAIL would hide a crash and would have to be unpicked the moment the model
# is fixed. Run it by hand; B stays disabled until it becomes a normal passing gate.
#
# WHAT IT ASKS. The exact oracle is exact per factor -- reference_factor_loglik equals
# reference_pair_loglik on the same sequences and fragments, asserted in the registered suite. That
# is necessary and NOT sufficient. For the factors to be a valid factorisation of the locus
# likelihood, this must also hold:
#
#     whole_locus(candidate) - [ sum(local factors) + sum(boundary factors) ]  ==  constant
#
# constant across CANDIDATES and across PHASES. A constant cancels in every comparison, so ranking is
# preserved. A difference that moves with the candidate means the factors implement a different
# model, and any emission built on them inherits that.
#
# WHAT IT MEASURES (2026-09-01, compact fixture, 3 candidates):
#
#     flank  100   diffs  107.07  -83.40   43.58    spread 190.47
#     flank  200   diffs   63.89  -63.09   63.89    spread 126.98
#     flank  400   diffs   79.89  -47.09   79.89    spread 126.98
#     flank  800   diffs   79.89  -47.09   79.89    spread 126.98
#
# It CONVERGES as the flank grows, so it is not a context-truncation artifact, and it settles at
# 126.98 rather than 0. Note which candidate differs: cAD/cCB and cAB/cAB agree at 79.89 while
# cAB/cCD sits at -47.09 -- the opposite PHASE of the same allele content. Identical alleles,
# identical lengths, therefore identical exposure. So this cannot be repaired by subtracting
# overlapping-flank exposure.
#
# THE DIAGNOSIS THAT FOLLOWS. The moving part is placement mass, not exposure. The whole model sums a
# fragment's placements across the complete diploid locus; the factor model assigns the fragment by
# RECRUITMENT and then sums placements only inside that factor's cropped context. Recruitment
# assigns every observed fragment exactly once -- which the registered suite asserts -- but it does
# not prove that all of that fragment's likelihood mass lies in the assigned context. A fragment
# labelled `local` can still have alternative whole-locus placements whose availability depends on
# phase, and cropping them changes the model.
#
# CONFIRMED PER FRAGMENT (2026-09-01). The whole 126.98 is two fragments classified `local`:
#
#     cB_161 / cB_181   whole-locus contrib, phase AD/CB   -9.0361
#                       whole-locus contrib, phase AB/CD  -72.5266
#                       cropped local-factor contrib       -9.0361  (identical in both phases)
#                       phase-dependent term lost         +63.4905  each, 126.9809 together
#
# The local factor at target 0 holds alleles (0,1) in BOTH phases -- the two candidates carry the
# same alleles at each block and differ only in pairing -- so its score cannot depend on phase by
# construction. The fragments' WHOLE-LOCUS mass does depend on it (-5.989 against -85.698), because
# they have alternative placements elsewhere in the diploid whose availability changes with phase.
# Cropping to the local context deletes exactly that.
#
# So: recruitment says where a fragment was DISCOVERED; placement mass says everywhere it could have
# arisen, and only the latter determines which genotype variables the likelihood depends on.
# Recruitment incidence is sound for diagnostics and acceleration and is NOT by itself a valid
# likelihood factorisation.
#
# NOT PROVEN, and worth stating because the same fixture invites the wrong inference: the 40
# boundary-class fragments reconcile to 0.0000 here, and that does NOT establish that boundary-class
# fragments are generally safe. ANY class may hold placements outside its cropped scope. The rule is
# the candidate-independent union of plausible placement dependencies, never the recruitment label.
#
# Exposure/dosage must likewise be charged once over a disjoint start-state space, and a fragment
# must enter as a single log-sum over its plausible origins: log SUM_z P(f,z|G) is not SUM_z log
# P(f,z|G), so a fragment cannot be split across factors at all.
#
# WHY AN EARLIER RUN LOOKED FINE. On the main fixture -- 350 bp fragments, 400 bp bubbles -- the
# evidence is 112 local against 78 boundary and the spread came out 0.000000. That fixture barely
# exercises the transition factor. This one is 2 local against 40 boundary and stresses exactly the
# factor whose composition is in question. A fixture that cannot see the failure is not evidence.
# EXIT CODES, and they must be distinguishable:
#   0  the factorisation reconciles within tolerance -- this script has become a passing gate
#   1  a valid measurement, and a candidate-dependent mismatch remains (the current state)
#   2  the instrument failed -- a scorer errored or produced nothing, so nothing was measured
#
# NOT `set -e` alone: a scorer exiting 1 would then terminate this script with status 1, which is
# indistinguishable from "mismatch reproduced". Every invocation is wrapped instead.
set -uo pipefail
BIN="${1:?usage: genotype_frag_factorisation.sh <panvar> <outdir>}"
OUT="${2:?}"; mkdir -p "$OUT"; OUT="$OUT/run.$$"; rm -rf "$OUT"; mkdir -p "$OUT"
PY_BIN="${PYTHON:-python3}"
seq_of() { awk -v n="$1" -v s="$2" 'BEGIN{x=s; b="ACGT";
           for(i=0;i<n;i++){x=(1103515245*x+12345)%2147483648; printf "%s", substr(b,(int(x/65536)%4)+1,1)} }'; }

L=$(seq_of 200 71); M=$(seq_of 40 72); N=$(seq_of 200 73)
X1=$(seq_of 60 74); X2=$(seq_of 60 75); Y1=$(seq_of 60 76); Y2=$(seq_of 60 77)
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\n' "$L"; printf 'S\t2\t%s\n' "$X1"; printf 'S\t3\t%s\n' "$X2"
  printf 'S\t4\t%s\n' "$M"; printf 'S\t5\t%s\n' "$Y1"; printf 'S\t6\t%s\n' "$Y2"
  printf 'S\t7\t%s\n' "$N"
  for a in 2 3; do printf 'L\t1\t+\t%s\t+\t0M\nL\t%s\t+\t4\t+\t0M\n' "$a" "$a"; done
  for a in 5 6; do printf 'L\t4\t+\t%s\t+\t0M\nL\t%s\t+\t7\t+\t0M\n' "$a" "$a"; done
  printf 'P\tref\t1+,2+,4+,5+,7+\t*\n'
  printf 'P\tcAB\t1+,2+,4+,5+,7+\t*\n'; printf 'P\tcCD\t1+,3+,4+,6+,7+\t*\n'
  printf 'P\tcAD\t1+,2+,4+,6+,7+\t*\n'; printf 'P\tcCB\t1+,3+,4+,5+,7+\t*\n'
} > "$OUT/g.gfa"
"$BIN" bubble -i "$OUT/g.gfa" -r ref -o "$OUT/b" --min-variant-bp 0 -q >/dev/null 2>&1

H1="${L}${X1}${M}${Y2}${N}"; H2="${L}${X2}${M}${Y1}${N}"
emit() { awk -v s="$1" -v tag="$2" 'BEGIN{ rl=60; ins=160; n=length(s);
    for(i=1;i+ins-1<=n;i+=20){ r1=substr(s,i,rl); r2=substr(s,i+ins-rl,rl);
      rc=""; for(j=length(r2);j>0;j--){c=substr(r2,j,1);
        rc=rc (c=="A"?"T":c=="C"?"G":c=="G"?"C":"A")}
      printf ">%s_%d/1\n%s\n>%s_%d/2\n%s\n", tag,i,r1,tag,i,rc } }'; }
emit "$H1" cA >  "$OUT/reads.fa"
emit "$H2" cB >> "$OUT/reads.fa"

"$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/b" -o "$OUT/d" \
  --dump-scored-sequences "$OUT/d" -q >/dev/null 2>&1
"$PY_BIN" - "$OUT/d.scored_sequences.fa" "$OUT" <<'PYEOF'
import sys
n=None;buf=[];seq={}
for l in open(sys.argv[1]):
    if l[0]=='>':
        if n: seq[n]="".join(buf)
        n=l[1:].split()[0]; buf=[]
    else: buf.append(l.strip())
if n: seq[n]="".join(buf)
for k in ("cAB","cCD","cAD","cCB"):
    if k in seq: open(f"{sys.argv[2]}/{k}.fa","w").write(f">{k}\n{seq[k]}\n")
PYEOF

echo "fragment composition of this fixture:"
"$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/b" -o "$OUT/i" -R "$OUT/reads.fa" \
  --incidence "$OUT/inc.tsv" --fragment-len 160 --fragment-sd 25 -q >/dev/null 2>&1
awk -F'\t' 'NR>1{c[$4]++} END{for(k in c) printf "  %-12s %s\n", k, c[k]}' "$OUT/inc.tsv"
echo

CF=(--error-rate 0.01 --fragment-len 160 --fragment-sd 25 --frag-outlier 0.05)

# Any failure here is INSTRUMENT failure (exit 2), never the model mismatch this script reports.
run_score() {
  local label="$1"; shift
  local out status
  out=$("$@" 2>/dev/null | tail -1); status=$?
  if [ "$status" -ne 0 ]; then
    echo "INSTRUMENT FAILURE: $label exited $status" >&2
    printf '  command:'; printf ' %q' "$@"; printf '\n' >&2
    exit 2
  fi
  if [ -z "$out" ]; then
    echo "INSTRUMENT FAILURE: $label produced no output" >&2
    exit 2
  fi
  printf '%s' "$out"
}
MAXSPREAD=0
TOL=0.001
printf "%8s  %-36s %s\n" flank "whole - sum(factors), per candidate" spread
for FB in 100 200 400 800; do
  DL=""
  while read -r P Q a0 a1 b0 b1; do
    [ -z "$P" ] && continue
    W=$(run_score "whole $P/$Q" "$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/b" -o "$OUT/w" \
         -R "$OUT/reads.fa" --reference-score "$OUT/$P.fa" "$OUT/$Q.fa" "${CF[@]}" -q)
    A=$(run_score "local0 $P/$Q" "$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/b" -o "$OUT/f0" \
         -R "$OUT/reads.fa" --reference-block 0 "$a0" "$b0" --flank-bp "$FB" "${CF[@]}" -q)
    B=$(run_score "local1 $P/$Q" "$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/b" -o "$OUT/f1" \
         -R "$OUT/reads.fa" --reference-block 1 "$a1" "$b1" --flank-bp "$FB" "${CF[@]}" -q)
    C=$(run_score "boundary $P/$Q" "$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/b" -o "$OUT/fb" \
         -R "$OUT/reads.fa" --reference-block-pair 0 "$a0" "$a1" "$b0" "$b1" --flank-bp "$FB" \
         "${CF[@]}" -q)
    D=$("$PY_BIN" -c "print('%.4f' % (float('$W')-(float('$A')+float('$B')+float('$C'))))") || {
      echo "INSTRUMENT FAILURE: could not difference the scores" >&2; exit 2; }
    DL="$DL $D"
  done <<'EOF'
cAD cCB 0 1 1 0
cAB cCD 0 0 1 1
cAB cAB 0 0 0 0
EOF
  SP=$("$PY_BIN" -c "
v=[float(x) for x in '$DL'.split()]
print('%.6f' % (max(v)-min(v)))") || { echo "INSTRUMENT FAILURE: spread" >&2; exit 2; }
  printf '%8s  %-36s %s\n' "$FB" "$DL" "$SP"
  MAXSPREAD=$("$PY_BIN" -c "print('%.6f' % max(float('$MAXSPREAD'), float('$SP')))")
done
echo
if "$PY_BIN" -c "import sys; sys.exit(0 if float('$MAXSPREAD') < float('$TOL') else 1)"; then
  echo "PASS: whole-locus minus sum(factors) is candidate-independent (max spread $MAXSPREAD)."
  echo "This script has become a normal gate. Plan step B is unblocked with respect to it."
  exit 0
fi
echo "EXPECTED: spread 0 (a candidate-independent constant). OBSERVED: max spread $MAXSPREAD."
echo "cAD/cCB and cAB/cCD hold identical alleles and lengths, so exposure is identical;"
echo "the phase-dependent term is placement mass cropped by the factor's context."
echo
echo "KNOWN FAILURE. Plan step B stays disabled until this is a passing gate."
exit 1
