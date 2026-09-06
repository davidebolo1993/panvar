#!/usr/bin/env bash
# c4_na19240_certified_interval.sh - the six-diplotype C4 experiment, frozen.
#
#   c4_na19240_certified_interval.sh <out-dir> <repo> <panvar-binary>
#
# WHAT IT PINS. Under the accelerated production scorer, c4/NA19240 leave-ZERO-out returns the
# COMPETITOR pair with the truth-equivalent pair 461.55 nats behind (seed 7). Under bounded-complete
# recruitment with certified placement mass and the SAME Hamming observation model, the truth-
# equivalent pair wins and is CERTIFIED.
#
# THE CONCLUSION THIS SUPPORTS, and no more: in this three-haplotype experiment, incomplete
# recruitment is SUFFICIENT to explain the production miscall. It shows the emission is adequate for
# THIS contrast, not that it is generally defect-free, and it is not full-caller certification --
# that needs the entire shortlist, where other candidates may enter the plausible set.
#
# ATTRIBUTION, decomposed rather than assumed. 12697 fragment-candidate cells widened under complete
# recruitment, including fragments outside the audited 76, so "the 20 reclassified fragments caused
# the reversal" does not follow from the reversal alone:
#
#   group                            n   accelerated     complete       swing
#   A original 20 reclassified      20      -3384.87         0.05    +3384.93
#   B remaining 54 audited          54      +2150.70      +649.86    -1500.84
#   C other 23877 (outside audit) 23877       +735.00     +2052.52    +1317.52
#     discordant 2                   2         -0.00        -0.00       +0.00
#   TOTAL                        23953       -499.17     +2702.43    +3201.60
#
# Group A dominates, but B moves AGAINST truth and C is a third of the swing. Only C lies outside
# the audited set.
#
# AN OPEN HYPOTHESIS, deliberately not asserted: C's +1317.52 against an exposure delta of -1317.65
# agree to 0.13 nats in ~1300. That suggests complete recruitment restores the diffuse event mass
# compensating the LONGER truth haplotype's exposure penalty (252695 bp against 226342), with the
# audited fragments supplying the concentrated effect. A contrast with the OPPOSITE length
# difference would test it -- the cancellation should reverse sign. Until then it is a coincidence
# worth chasing, not a mechanism.
#
# TOLERANCES, not floating-point strings: reads, haplotypes, parameters, fragment count and class
# membership are pinned; scores are asserted by inequality. The binary hash is recorded as
# PROVENANCE only -- making it a precondition would break this on every unrelated rebuild.
set -uo pipefail
OUT="${1:?out}"; REPO="${2:?repo}"; BIN="${3:?bin}"
mkdir -p "$OUT"; PY="${PYTHON:-python3}"
md5of() { if command -v md5 >/dev/null 2>&1; then md5 -q "$1"; else md5sum "$1" | cut -d' ' -f1; fi; }
fails=0
ok()  { printf '  ok   %s\n' "$*"; }
bad() { printf '  FAIL %s\n' "$*"; fails=$((fails+1)); }

G="$REPO/results/real_data/c4/bubble/bubble.sorted.gfa"
P="$REPO/results/real_data/c4/bubble/bubble"
T1="NA19240#1#haplotype1-0000009:31981961-32208331"
TE="NA19238#1#haplotype1-0000026:31825469-32078164"
WC="HG00171#1#haplotype1-0000033:140888062-141114404"
SEED=7; COV=30; DIV=0.05; FLEN=350; FSD=50; ERR=0.001; TOL=1.0
EXP_FRAGS=23953
EXP_TRUTH_MD5=cee31ff8f87da6caa4ed278d1da2e7ed
EXP_R1_MD5=d4241737877d77c816ff2021ea9c411e

"$BIN" genotype-frag -i "$G" -b "$P" -o "$OUT/d" --dump-scored-sequences "$OUT/d" -q >/dev/null 2>&1
"$PY" - "$OUT/d.scored_sequences.fa" "$OUT/d.scored_sequences.tsv" "$OUT" "$T1" "$TE" "$WC" \
  >/dev/null <<'PYEOF'
import sys, hashlib
fa,tsv,out,t1,te,wc=sys.argv[1:7]
frame={};h=None
for l in open(tsv):
    f=l.rstrip('\n').split('\t')
    if h is None: h=f; continue
    d=dict(zip(h,f)); frame[d['name']]=d.get('frame','fwd')
seqs={};n=None;b=[]
def flush():
    if n is not None: seqs[n]=''.join(b)
for l in open(fa):
    if l[0]=='>': flush(); n=l[1:].split()[0]; b=[]
    else: b.append(l.strip())
flush()
def rc(s): return s[::-1].translate(str.maketrans('ACGTN','TGCAN'))
for role,nm in (("truth_h1",t1),("truth_equiv_h2",te),("competitor",wc)):
    s=seqs[nm].upper()
    o=rc(s) if frame.get(nm)=='rc' else s
    open(f"{out}/{role}.fa","w").write(">%s\n%s\n"%(role,o))
PYEOF
# The record NAMES must match the audit script exactly: they enter the FASTA and
# therefore the wgsim read names, so different labels give different reads and the
# pinned md5s stop linking the two experiments.
cat "$OUT/truth_h1.fa" "$OUT/truth_equiv_h2.fa" > "$OUT/truth.fa"
[ "$(md5of "$OUT/truth.fa")" = "$EXP_TRUTH_MD5" ] \
  && ok "truth FASTA reproduces its pinned md5" \
  || bad "truth FASTA md5 is $(md5of "$OUT/truth.fa"), expected $EXP_TRUTH_MD5"
L=$(awk '/^>/{next}{n+=length($0)}END{print n}' "$OUT/truth.fa"); N=$(( COV * L / 600 ))
wgsim -N "$N" -1 150 -2 150 -d "$FLEN" -s "$FSD" -e "$ERR" -r 0 -R 0 -X 0 -S "$SEED" \
  "$OUT/truth.fa" "$OUT/r1.fq" "$OUT/r2.fq" >/dev/null 2>&1
[ "$(md5of "$OUT/r1.fq")" = "$EXP_R1_MD5" ] \
  && ok "reads reproduce their pinned md5 ($N pairs)" \
  || bad "read md5 is $(md5of "$OUT/r1.fq"), expected $EXP_R1_MD5"

"$BIN" genotype-frag -i "$G" -b "$P" -o "$OUT/is" -R "$OUT/r1.fq" -R "$OUT/r2.fq" \
  --max-divergence "$DIV" --fragment-len "$FLEN" --fragment-sd "$FSD" --error-rate "$ERR" \
  --interval-score "$OUT/is.tsv" --interval-tol "$TOL" --interval-tau 0.0 \
  --interval-candidates "$T1,$TE,$WC" -q >/dev/null 2>&1
if [ ! -s "$OUT/is.tsv" ]; then
  bad "no interval score produced"
else
  NF=$(grep -E "^# fragments" "$OUT/is.tsv" | cut -f2)
  VD=$(grep -E "^# verdict" "$OUT/is.tsv" | cut -f2)
  GAP=$(grep -E "^# robust_gap" "$OUT/is.tsv" | cut -f2)
  NP=$(grep -E "^# n_plausible" "$OUT/is.tsv" | cut -f2)
  AG=$(grep -E "^# leaders_agree" "$OUT/is.tsv" | cut -f2)
  LEAD=$(grep -E "^# robust_leader" "$OUT/is.tsv" | cut -f2)
  NCLS=$(awk -F'\t' 'NR>1 && $1!~/^#/{n++} END{print n+0}' "$OUT/is.tsv")
  [ "${NF:-0}" = "$EXP_FRAGS" ] && ok "all $EXP_FRAGS fragments scored" \
                               || bad "scored ${NF:-0} fragments, expected $EXP_FRAGS"
  [ "${NCLS:-0}" = 6 ] && ok "all six unordered diplotypes scored" \
                       || bad "${NCLS:-0} diplotypes scored, expected 6"
  [ "${VD:-}" = CERTIFIED ] && ok "verdict is CERTIFIED" || bad "verdict is ${VD:-none}"
  [ "${NP:-0}" = 1 ] && ok "exactly one plausible class" || bad "${NP:-?} plausible classes"
  [ "${AG:-0}" = 1 ] && ok "robust leader and nominal winner agree" \
                     || bad "robust leader and nominal winner DIFFER"
  case "$LEAD" in
    *NA19238*) ok "the leader is the TRUTH-EQUIVALENT pair" ;;
    *) bad "the leader is not the truth-equivalent pair: $LEAD" ;;
  esac
  # INEQUALITIES, not exact strings: the gap must be comfortably positive and widths small.
  "$PY" -c "import sys; sys.exit(0 if float('${GAP:-0}') > 500.0 else 1)" 2>/dev/null \
    && ok "robust gap is comfortably positive ($GAP nats)" \
    || bad "robust gap is ${GAP:-?}, expected > 500"
  MW=$(awk -F'\t' 'NR>1 && $1!~/^#/ && $7>m{m=$7} END{printf "%.4f",m}' "$OUT/is.tsv")
  "$PY" -c "import sys; sys.exit(0 if float('$MW') < 1.0 else 1)" 2>/dev/null \
    && ok "every class interval is below tolerance (widest $MW)" \
    || bad "widest class interval is $MW, expected < 1.0"
  # The production-band classification must still reproduce the original audit partition.
  PB=$(awk -F'\t' -v te="$TE" 'NR>1 && $2==te && $3==1{n++} END{print n+0}' "$OUT/is.tsv.placements.tsv")
  # Measured 19930 of 23953 (~83%). The floor is a sanity bound on the classification still
  # working, not a pinned count -- a pinned count would break on any read-simulation change.
  [ "${PB:-0}" -gt 19000 ] \
    && ok "production-band placements on the truth-equivalent haplotype: $PB of $EXP_FRAGS" \
    || bad "only ${PB:-0} production-band placements, expected > 19000"
fi
echo "  provenance: binary $(md5of "$BIN"), commit $(cd "$REPO" && git rev-parse --short HEAD 2>/dev/null)"
echo
if [ "$fails" -eq 0 ]; then echo "c4 certified interval: all assertions passed"; else
  echo "c4 certified interval: $fails assertion(s) failed"; fi
exit "$fails"
