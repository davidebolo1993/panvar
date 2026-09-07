#!/usr/bin/env bash
# genotype_lazy_refinement.sh - candidate-level LAZY refinement against the EAGER scorer.
#
#   genotype_lazy_refinement.sh <panvar> <out-dir>
#
# WHAT LAZY REFINEMENT IS. Every fragment x candidate cell gets a cheap production-band bound. Class
# intervals are formed from those bounds, and only cells participating in a class that is still
# plausible are deepened. A CANDIDATE is discarded only when EVERY genotype class containing it
# satisfies U(C) < B - tau -- never on a per-candidate score. A candidate can be hopeless with one
# partner and decisive with another, so any per-haplotype cutoff would reopen the recruitment defect
# this line of work exists to characterise.
#
# WHAT IS ASSERTED HERE, and the asymmetry is deliberate:
#   * the two arms agree on the ROBUST LEADER, the PLAUSIBLE SET and the VERDICT;
#   * for classes lazy ELIMINATED, its interval need only CONTAIN the eager interval -- demanding
#     numerical agreement there would forbid the pruning that is the entire point;
#   * for the classes that decide the answer -- the leader and the rival that sets the robust gap --
#     the intervals must agree numerically, and the gap with them;
#   * elimination is NON-VACUOUS: some round leaves fewer candidates alive than the panel holds;
#   * lazy never deepens a cell the eager arm would not have deepened, and does strictly fewer
#     depth increments;
#   * raising tau only ever ENLARGES the plausible set. A lazy arm that returns a SMALLER set than
#     the eager one at the same tau has silently changed the answer, not accelerated it.
#
# WHAT THIS FIXTURE CANNOT DISCRIMINATE, recorded because it was measured rather than assumed.
# Replacing the screen's safely_prunable(C, B, tau) with safely_prunable(C, B, 0.0) -- pruning at
# U < B, the rule the bounded-search contract shows is unsafe -- leaves every assertion above
# passing. The reason is a property of the fixture, not a hole in the assertions: the lazy arm
# deepens ZERO cells here, so the screen has nothing to schedule and its tolerance cannot matter.
# All 1502 cells the eager arm deepens turn out to contribute exactly zero width to all 21 classes
# -- wide in raw mass, but entirely below the background floor once mixed -- which is also why both
# arms report the gap identically to ten digits.
#
# That mutation is also not a SOUNDNESS defect in this implementation, and calling it one would
# misplace where the safety lives. The screen's tau schedules refinement; it never removes a class
# from the output. Every class is re-tested at the end by certify() at the declared tau, from
# whatever intervals it ended with, so pruning too eagerly can only leave a class coarse and widen
# the final plausible set -- an UNRESOLVED where eager said CERTIFIED, which the n_plausible
# assertion does catch. It cannot manufacture a certification. Using tau here is a COMPLETENESS
# safeguard against that spurious UNRESOLVED, not the soundness rule.
#
# NOT ASSERTED HERE: anything about the real C4 panel. This is an analytic fixture; the production
# comparison against the frozen 48-candidate oracle is a separate, much longer run.
set -uo pipefail
BIN="${1:?usage: genotype_lazy_refinement.sh <panvar> <outdir>}"
OUT="${2:?}"; mkdir -p "$OUT"; OUT="$OUT/run.$$"; rm -rf "$OUT"; mkdir -p "$OUT"
PY="${PYTHON:-python3}"
fails=0
ok()  { printf '  ok   %s\n' "$*"; }
bad() { printf '  FAIL %s\n' "$*"; fails=$((fails+1)); }

# SIX haplotypes on a shared backbone. hapA and hapB differ by a handful of substitutions, so their
# classes contest; hapC..hapF carry an unrelated tail each, so their classes should fall out on the
# coarse bound alone. Reads come from hapA only, making A|A the truth.
"$PY" - "$OUT" <<'PYEOF'
import sys, random
out = sys.argv[1]
rnd = random.Random(20260907)
B = "ACGT"
def seq(n, r): return "".join(r.choice(B) for _ in range(n))
back = seq(3000, random.Random(11))
tailA = seq(3000, random.Random(12))
def mutate(s, k, r):
    s = list(s)
    for p in r.sample(range(len(s)), k):
        s[p] = r.choice([c for c in B if c != s[p]])
    return "".join(s)
haps = {
    "hapA": back + tailA,
    "hapB": back + mutate(tailA, 40, random.Random(13)),
    "hapC": back + seq(3000, random.Random(14)),
    "hapD": back + seq(3000, random.Random(15)),
    "hapE": back + seq(3000, random.Random(16)),
    "hapF": back + seq(3000, random.Random(17)),
}
names = list(haps)
with open(out + "/g.gfa", "w") as g:
    g.write("H\tVN:Z:1.0\n")
    for i, n in enumerate(names, 1):
        g.write("S\t%d\t%s\n" % (i, haps[n]))
    for i in range(1, len(names)):
        g.write("L\t%d\t+\t%d\t+\t0M\n" % (i, i + 1))
    for i, n in enumerate(names, 1):
        g.write("P\t%s\t%d+\t*\n" % (n, i))
comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
h = haps["hapA"]
with open(out + "/r1.fq", "w") as f1, open(out + "/r2.fq", "w") as f2:
    k = 0
    for i in range(0, len(h) - 360, 8):
        a = h[i:i + 150]; b = h[i + 200:i + 350]
        if len(a) < 150 or len(b) < 150: break
        rc = "".join(comp[c] for c in reversed(b))
        f1.write("@f%d/1\n%s\n+\n%s\n" % (k, a, "I" * 150))
        f2.write("@f%d/2\n%s\n+\n%s\n" % (k, rc, "I" * 150))
        k += 1
open(out + "/cands.txt", "w").write(",".join(names))
PYEOF
N=$(( $(wc -l < "$OUT/r1.fq") / 4 ))
CAND=$(cat "$OUT/cands.txt")
[ "$N" -gt 100 ] && ok "fixture has $N fragments over 6 candidates" || bad "only $N fragments"

run() { # run <tag> <tau> [extra...]
  "$BIN" genotype-frag -i "$OUT/g.gfa" -b "$OUT/b" -o "$OUT/$1" \
    -R "$OUT/r1.fq" -R "$OUT/r2.fq" -t 2 --max-divergence 0.05 --fragment-len 350 \
    --fragment-sd 50 --error-rate 0.001 --interval-score "$OUT/$1.tsv" --interval-tol 1.0 \
    --interval-tau "$2" --interval-candidates "$CAND" "${@:3}" -q; }

"$BIN" bubble -i "$OUT/g.gfa" -r hapA -o "$OUT/b" --min-variant-bp 0 -q >/dev/null 2>&1
run eager 0.0 >/dev/null 2>&1
run lazy  0.0 --interval-lazy --interval-rounds "$OUT/rounds.tsv" >/dev/null 2>&1
[ -s "$OUT/eager.tsv" ] && [ -s "$OUT/lazy.tsv" ] \
  && ok "both arms produced a result" || { bad "an arm produced nothing"; echo; exit "$fails"; }

hdr() { grep -E "^# $2	" "$OUT/$1.tsv" | cut -f2; }
for f in verdict robust_leader n_plausible; do
  E=$(hdr eager "$f"); L=$(hdr lazy "$f")
  [ "$E" = "$L" ] && ok "$f agrees: $E" || bad "$f differs: eager=$E lazy=$L"
done
[ "$(hdr eager verdict)" = CERTIFIED ] \
  && ok "the fixture certifies (otherwise this compares two non-answers)" \
  || bad "the fixture did not certify -- verdict $(hdr eager verdict)"

# The ledger. deepened_cells must be a SUBSET of what eager would deepen, and strictly smaller.
CN=$(hdr lazy coarse_nontolerant); LD=$(hdr lazy deepened_cells)
ES=$(hdr eager level_steps);        LS=$(hdr lazy level_steps)
EC=$(hdr eager deepened_cells)
[ "${CN:-0}" = "${EC:-1}" ] \
  && ok "the two arms agree on the eligible set: lazy coarse_nontolerant $CN = eager deepened $EC" \
  || bad "lazy says $CN cells miss the tolerance at level 1, eager deepened $EC -- one is wrong"
[ "${LD:-0}" -le "${CN:-0}" ] \
  && ok "lazy deepened $LD cells, within the $CN the eager rule would deepen" \
  || bad "lazy deepened $LD cells but only $CN are eligible -- it went deeper than eager"
[ "${LD:-1}" -lt "${EC:-0}" ] \
  && ok "lazy deepened strictly fewer cells than eager ($LD < $EC)" \
  || bad "lazy deepened $LD cells, eager $EC -- no saving"
[ "${LS:-1}" -lt "${ES:-0}" ] \
  && ok "lazy did strictly fewer depth increments than eager ($LS < $ES)" \
  || bad "lazy did $LS depth increments, eager $ES -- no saving"

# Elimination must actually happen, and must never reverse.
if [ -s "$OUT/rounds.tsv" ]; then
  MINC=$(awk -F'\t' 'NR>1 && $1=="screen" && $4>0 {if(m==""||$4<m)m=$4} END{print m+0}' "$OUT/rounds.tsv")
  [ "${MINC:-6}" -lt 6 ] \
    && ok "elimination is non-vacuous: $MINC of 6 candidates alive at the tightest round" \
    || bad "no candidate was ever eliminated -- the coarse bounds screened nothing"
  NONINC=$(awk -F'\t' 'NR>1 && $1=="screen"{if(p!=""&&$3>p)b=1; p=$3} END{print b+0}' "$OUT/rounds.tsv")
  [ "${NONINC:-1}" = 0 ] \
    && ok "the plausible set never grows across screening rounds" \
    || bad "a class re-entered the plausible set -- pruning is not monotone"
else
  bad "no round log written"
fi

# Interval agreement, asymmetric by design.
"$PY" - "$OUT/eager.tsv" "$OUT/lazy.tsv" "$(hdr eager robust_leader)" <<'PYEOF'
import sys
def load(p):
    rows, hdr = {}, {}
    for l in open(p):
        if l.startswith("#"):
            f = l[2:].rstrip("\n").split("\t")
            if len(f) == 2: hdr[f[0]] = f[1]
            continue
        f = l.rstrip("\n").split("\t")
        if f[0] == "class": continue
        rows[f[0]] = (float(f[4]), float(f[5]))
    return rows, hdr
E, EH = load(sys.argv[1]); L, LH = load(sys.argv[2]); lead = sys.argv[3]
bad = 0
def ok(m): print("  ok   " + m)
def no(m):
    global bad; bad += 1; print("  FAIL " + m)
if set(E) != set(L): no("the two arms scored different class sets")
else: ok("both arms scored the same %d classes" % len(E))
slop = 1e-6
viol = [c for c in E if not (L[c][0] <= E[c][0] + slop and L[c][1] >= E[c][1] - slop)]
if viol: no("%d lazy intervals do not contain the eager interval, e.g. %s" % (len(viol), viol[0]))
else: ok("every lazy interval safely contains the eager interval")
# The leader and the gap-setting rival must agree numerically, not merely contain.
rival = max((c for c in E if c != lead), key=lambda c: E[c][1])
tol = float(EH.get("interval_tol", "1.0"))
for role, c in (("leader", lead), ("gap-setting rival", rival)):
    d = max(abs(L[c][0] - E[c][0]), abs(L[c][1] - E[c][1]))
    if d <= tol: ok("%s interval agrees within %.4g nats (tolerance %g)" % (role, d, tol))
    else: no("%s interval differs by %.4g nats, tolerance %g" % (role, d, tol))
ge, gl = float(EH["robust_gap"]), float(LH["robust_gap"])
if abs(ge - gl) <= 2 * tol:
    ok("robust gap agrees: eager %.4f, lazy %.4f (within 2x tolerance)" % (ge, gl))
else:
    no("robust gap differs: eager %.4f, lazy %.4f" % (ge, gl))
sys.exit(bad)
PYEOF
fails=$(( fails + $? ))

# TAU MONOTONICITY. tau only enlarges the plausible set; the lazy arm must never return fewer
# classes than the eager one at the same tau, which is how an "optimisation" hides a wrong answer.
# tau=300 exceeds this fixture's 252-nat gap ON PURPOSE: it puts the rival back into the
# plausible set, so the two arms have a MULTI-CLASS set to disagree about. At a tau below the gap
# both arms trivially report one class and the comparison asserts nothing.
run eager_t 300.0 >/dev/null 2>&1
run lazy_t  300.0 --interval-lazy >/dev/null 2>&1
PE=$(grep -E "^# n_plausible	" "$OUT/eager_t.tsv" | cut -f2)
PL=$(grep -E "^# n_plausible	" "$OUT/lazy_t.tsv" | cut -f2)
[ "${PE:-0}" = "${PL:-1}" ] \
  && ok "at tau=300 both arms report the same plausible set size ($PE)" \
  || bad "at tau=300 eager reports $PE plausible classes, lazy $PL"
P0=$(hdr eager n_plausible)
[ "${PE:-0}" -gt "${P0:-0}" ] \
  && ok "raising tau from 0 to 300 did not shrink the plausible set ($P0 -> $PE)" \
  || bad "tau=300 did not enlarge the plausible set ($P0 -> $PE) -- the tau assertions are vacuous"

echo
if [ "$fails" -eq 0 ]; then echo "lazy refinement: all assertions passed"; else
  echo "lazy refinement: $fails assertion(s) failed"; fi
exit "$fails"
