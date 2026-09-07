#!/usr/bin/env bash
# genotype_hybrid_chain.sh - the hybrid block caller: marker unaries + fragment linkage.
#
#   genotype_hybrid_chain.sh <panvar> <out-dir>
#
# THE ESSENTIAL RULE, pinned before the chain is built: marker unary factors and fragment linkage
# factors MUST NOT COUNT THE SAME READ EVIDENCE TWICE. The design is EVIDENCE OWNERSHIP --
#
#   1. a fragment depending on ONE variable block feeds that block's marker unary;
#   2. a fragment depending on TWO variable blocks feeds their linkage factor exactly once,
#      sequence evidence included;
#   3. a fragment depending on THREE OR MORE variable blocks gets a WIDER FACTOR, or an explicit
#      refusal. It is never CROPPED into narrower factors -- cropping is what
#      tests/genotype_frag_factorisation.sh already refuted;
#   4. boundary fragments are EXCLUDED from the marker counts feeding the unaries.
#
# ARITY IS COUNTED IN VARIABLE BLOCKS, NOT PHYSICAL BLOCKS CROSSED. A block whose sequence is
# identical across every candidate carries no genotype state and no phase decision, so a fragment
# spanning `variable A -- fixed backbone -- variable B` still defines the PAIRWISE factor psi(A,B):
# the backbone enters as sequence context and nothing more. Counting physical blocks would make that
# fragment Wide and throw away most of the real linkage evidence, because real bubbles commonly have
# reference sequence between them. An intervening block that DOES vary is itself a variable, so it
# raises the arity to three by construction -- there is no way for it to be carried as context.
#
# Ownership is a PARTITION, so "counted twice" and "silently dropped" are the same assertion from
# two sides, and gate 3 checks both directions at once.
#
# WHAT GATE 3 DOES NOT YET COVER, stated so it is not read as more than it is. That a fragment has
# one owner is partly structural -- FragmentOwner holds a single kind -- so the assertion mainly
# catches fragments DROPPED or MIS-CLASSIFIED, not a consumer that double counts. The rule that
# actually needs guarding is 4: the marker counts feeding the unaries must EXCLUDE the linkage
# fragments. That is only testable once the chain exists, and it is the first thing gate 1 must
# pin. Measured meanwhile: collapsing an origin's block span to its start block turns all 2960
# fragments unary and fires two assertions here, so the linkage detection itself is not vacuous.
#
# THE GATES. 3 and 4 run now, against --ownership-table. 1, 2 and 7 need the chain itself and are
# held behind --hybrid-call; this file FAILS rather than skips once that entry point exists.
#   1. with NO linkage fragments the hybrid output equals the existing marker caller EXACTLY;
#   2. with identical per-block allele CONTENT but different PHASE, boundary fragments select the
#      correct phase -- the one thing a per-block caller cannot do by construction;
#   3. every fragment is owned by exactly one factor;                                    [ACTIVE]
#   4. no factor loses candidate-dependent placement mass;                               [ACTIVE]
#   5. two variable bubbles separated by a FIXED backbone still form a pairwise factor;   [ACTIVE]
#   6. an intervening VARIABLE block gives a three-variable factor, never a cropped pair; [ACTIVE]
#   7. the exposure precondition holds, and is measured rather than assumed;             [ACTIVE]
#   8. every ownership class is accounted, and its mass share reported;                  [ACTIVE]
#   9. ambiguous evidence yields an equivalence set or UNRESOLVED, never a confident guess.
#
# NORMALISATION is settled in genotype_fragments.hpp. Linkage is a CONDITIONAL PHASE SCORE, and TWO
# preconditions come with it -- both of which an earlier version of this contract got wrong:
#
#   * THE BACKGROUND DOES NOT CANCEL. log(A_x + eta*P_bg) - log(A_y + eta*P_bg) != log A_x - log A_y.
#     It is genotype-independent but sits inside the log, so it is a FLOOR that decides how much a
#     weakly-placing fragment may say about phase. Measured: a contrast near the floor is 0.0064
#     nats with the background retained and 10.0000 without -- dropping it manufactures the entire
#     signal. Each configuration must keep its own background inside its mixture. Held behind
#     --hybrid-call, since only the chain forms the mixture.
#   * EXPOSURE IS AFFINE ONLY ABOVE THE INSERT SUPPORT, so the cis/trans cancellation is conditional
#     on every constructed window clearing it. A deletion or bypass allele can produce a short
#     window; then the cancellation is simply false. Gate 7 measures the boundary directly.
#
# THE TRADEOFF, recorded not glossed: excluding linkage-owned fragments from the unaries while
# normalising the linkage factor conditional on endpoint content discards their CONTENT evidence and
# keeps only phase. That is not a lossless factorisation.
#
# GATE 8 MEASURES SIZE, NOT INFORMATION. Pooled in-band placement mass says how much mass a class
# holds, not how much a call depends on it -- a fragment with a negligible share can still carry a
# decisive likelihood RATIO, and a ratio is what a call turns on. So no mass figure here may be
# quoted as evidence that the content loss is harmless. Only the block-content regression, on C4 and
# the exact leave-zero-out controls, can establish that.
#
# EVERY CLASS NEEDS A DISPOSITION or the partition leaks: Unary is marker content; Linkage is
# conditional phase, content given up; Wide is NOT representable by a pairwise transition and must
# be reported UNSUPPORTED/INCOMPLETE rather than deleted or cropped; Invariant is ignorable for
# ranking but belongs to absolute-fit calibration; Unusable is explicitly reported missing evidence.
#
# WIDE IS NOT "UNRESOLVED". Unresolved means the implemented model evaluated the evidence and could
# not separate the states -- a statement about the data. Wide means the pairwise model could not
# CONSUME the evidence at all -- a statement about the model's reach. Reporting the second as the
# first would blame the data for a modelling limit, which is the same error as calling a run stopped
# by a depth cap "ambiguous".
#
# Gates 5 (C4 block 7 keeps the marker caller's 48/48) and 6 (exact leave-zero-out controls do not
# regress) are real-panel regressions and live in tests/regressions/: they need the C4 graph and
# simulated cohorts, and this file must stay runnable in ctest.
set -uo pipefail
BIN="${1:?usage: genotype_hybrid_chain.sh <panvar> <outdir>}"
OUT="${2:?}"; mkdir -p "$OUT"; OUT="$OUT/run.$$"; rm -rf "$OUT"; mkdir -p "$OUT"
PY="${PYTHON:-python3}"
fails=0
ok()  { printf '  ok   %s\n' "$*"; }
bad() { printf '  FAIL %s\n' "$*"; fails=$((fails+1)); }

# ---------------------------------------------------------------------------------------------
# THE PHASE FIXTURE. Two bubbles that ABUT through a shared boundary node, so they are ADJACENT
# blocks in the chain and a fragment crossing the junction owns exactly two of them:
#
#     X --> ( A1 | B1 ) --> J --> ( A2 | B2 ) --> Z
#
#     hapAA = A1,A2   hapAB = A1,B2   hapBA = B1,A2   hapBB = B1,B2
#
# The sample is hapAB / hapBA. Its per-block allele CONTENT is {A1,B1} at block 1 and {A2,B2} at
# block 2 -- identical to what hapAA/hapBB would give. Marker counts alone therefore cannot separate
# the two diplotypes at either block: only a fragment crossing J can say whether A1 sits in cis with
# A2 or with B2.
#
# THE BUBBLES MUST ABUT. Separated by a backbone stretch, a fragment spanning variant to variant
# touches THREE blocks and is Wide, not Linkage -- and chaining phase through the backbone carries
# nothing, because a constant block has no allele to phase against. This is a real constraint on the
# minimal adjacent-pair chain, not a fixture convenience, and it is why J is a shared node here.
"$PY" - "$OUT" <<'PYEOF'
import sys, random
out = sys.argv[1]; B = "ACGT"
def seq(n, s):
    r = random.Random(s); return "".join(r.choice(B) for _ in range(n))
def mut(s, k, sd):
    r = random.Random(sd); s = list(s)
    for p in r.sample(range(len(s)), k): s[p] = r.choice([c for c in B if c != s[p]])
    return "".join(s)
X = seq(4000, 1); A1 = seq(600, 2); B1 = mut(A1, 90, 3); J = seq(40, 4)
A2 = seq(600, 5); B2 = mut(A2, 90, 6); Z = seq(4000, 7)
segs = [("1", X), ("2", A1), ("3", B1), ("4", J), ("5", A2), ("6", B2), ("7", Z)]
links = [("1","2"),("1","3"),("2","4"),("3","4"),("4","5"),("4","6"),("5","7"),("6","7")]
paths = {"hapAA": ["1","2","4","5","7"], "hapAB": ["1","2","4","6","7"],
         "hapBA": ["1","3","4","5","7"], "hapBB": ["1","3","4","6","7"]}
with open(out + "/g.gfa", "w") as g:
    g.write("H\tVN:Z:1.0\n")
    for n, s in segs: g.write("S\t%s\t%s\n" % (n, s))
    for a, b in links: g.write("L\t%s\t+\t%s\t+\t0M\n" % (a, b))
    for n, st in paths.items(): g.write("P\t%s\t%s\t*\n" % (n, ",".join(x + "+" for x in st)))
sq = {n: "".join(dict(segs)[x] for x in st) for n, st in paths.items()}
comp = {"A":"T","C":"G","G":"C","T":"A"}
n_span = 0
with open(out + "/r1.fq","w") as f1, open(out + "/r2.fq","w") as f2:
    k = 0
    for nm in ("hapAB", "hapBA"):
        h = sq[nm]
        for i in range(0, len(h) - 360, 6):
            a = h[i:i+150]; b = h[i+200:i+350]
            if len(b) < 150: break
            rc = "".join(comp[c] for c in reversed(b))
            f1.write("@%s_%d/1\n%s\n+\n%s\n" % (nm, k, a, "I"*150))
            f2.write("@%s_%d/2\n%s\n+\n%s\n" % (nm, k, rc, "I"*150))
            if i < len(X)+len(A1) and i + 350 > len(X)+len(A1)+len(J): n_span += 1
            k += 1
open(out + "/junction.txt","w").write("%d %d" % (len(X)+len(A1), n_span))
PYEOF
"$BIN" bubble -i "$OUT/g.gfa" -r hapAA -o "$OUT/b" --min-variant-bp 0 -q >/dev/null 2>&1
NB=$(( $(wc -l < "$OUT/b.bubbles.csv") - 1 ))
[ "$NB" = 2 ] && ok "fixture decomposes into 2 abutting bubbles" || bad "$NB bubbles, expected 2"

"$BIN" genotype-frag -i "$OUT/b.sorted.gfa" -b "$OUT/b" -o "$OUT/o" -R "$OUT/r1.fq" -R "$OUT/r2.fq" \
  -t 2 --max-divergence 0.05 --fragment-len 350 --fragment-sd 50 --error-rate 0.001 \
  --ownership-table "$OUT/own.tsv" -q >/dev/null 2>&1
[ -s "$OUT/own.tsv" ] && ok "ownership table written" || { bad "no ownership table"; exit "$fails"; }

NFRAG=$(( $(wc -l < "$OUT/r1.fq") / 4 ))
# GATE 3: ownership is a PARTITION over every loaded fragment.
"$PY" - "$OUT/own.tsv" "$NFRAG" <<'PYEOF'
import sys, collections
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])][1:]
want = int(sys.argv[2]); bad = 0
def ok(m): print("  ok   " + m)
def no(m):
    global bad; bad += 1; print("  FAIL " + m)
seen = collections.Counter(r[0] for r in rows)
dup = [f for f, n in seen.items() if n != 1]
if dup: no("%d fragments own more than one factor, e.g. %s" % (len(dup), dup[0]))
elif len(seen) != want: no("%d fragments owned, %d loaded -- some were dropped" % (len(seen), want))
else: ok("all %d fragments owned by exactly one factor (a partition, no double counting)" % want)
kinds = collections.Counter(r[1] for r in rows)
print("       ownership: " + ", ".join("%s=%d" % kv for kv in sorted(kinds.items())))
for r in rows:
    # Linkage names two DISTINCT VARIABLE blocks. They need not be physically adjacent: a fixed
    # backbone between them is context. What must hold is the arity, not the block distance.
    if r[1] == "linkage" and (int(r[5]) != 2 or int(r[3]) <= int(r[2])):
        no("linkage factor with var arity %s over blocks %s..%s" % (r[5], r[2], r[3])); break
    if r[1] == "unary" and (int(r[5]) != 1 or r[2] != r[3]):
        no("unary factor with var arity %s over blocks %s..%s" % (r[5], r[2], r[3])); break
    if r[1] == "wide" and int(r[5]) < 3:
        no("wide factor with var arity %s -- wide must mean three or more variables" % r[5]); break
else:
    ok("unary owners name one variable, linkage two, wide three or more")
if kinds.get("linkage", 0) == 0:
    no("no fragment crosses the junction -- the phase gate would be vacuous")
else:
    ok("%d fragments own the junction linkage factor" % kinds["linkage"])
# Both homologues must supply boundary evidence; a one-sided set could fix phase by depth alone.
side = collections.Counter(r[0].split("_")[0] for r in rows if r[1] == "linkage")
if len(side) == 2 and min(side.values()) > 0:
    ok("both homologues supply junction fragments (%s)" % ", ".join("%s=%d" % kv for kv in sorted(side.items())))
else:
    no("junction fragments come from one homologue only: %s" % dict(side))
# GATE 4: no factor loses candidate-dependent placement mass. `dropped_nats` is what restricting a
# fragment to its owning factor's scope actually costs, measured over the WHOLE candidate set.
TOL = 1e-6
worst = max(float(r[12]) for r in rows)
if worst <= TOL: ok("no factor loses placement mass (worst restriction costs %.3e nats, tol %g)" % (worst, TOL))
else: no("a factor loses %.3e nats of placement mass, tolerance %g" % (worst, TOL))
uncert = [r for r in rows if r[13] != "1"]
if uncert: no("%d fragments have an UNCERTIFIED scope yet were still assigned" % len(uncert))
else: ok("every assigned scope is certified against the out-of-band bound")
sys.exit(bad)
PYEOF
fails=$(( fails + $? ))

# ---------------------------------------------------------------------------------------------
# GATES 5 AND 6. Two more fixtures, differing ONLY in whether the block between the two variable
# bubbles varies:
#
#   bb   variable A -- FIXED backbone -- variable B     must still give a PAIRWISE factor
#   vm   variable A -- VARIABLE middle -- variable B    must give a THREE-variable factor
#
# The pair is the point. Counting physical blocks crossed makes both of them Wide and discards the
# bb linkage evidence, which at a real locus is most of it; counting variable dependencies keeps bb
# pairwise and still refuses to crop vm.
"$PY" - "$OUT" <<'PYEOF2'
import sys, random
out = sys.argv[1]; B = "ACGT"
def seq(n, s):
    r = random.Random(s); return "".join(r.choice(B) for _ in range(n))
def mut(s, k, sd):
    r = random.Random(sd); s = list(s)
    for p in r.sample(range(len(s)), k): s[p] = r.choice([c for c in B if c != s[p]])
    return "".join(s)
comp = {"A":"T","C":"G","G":"C","T":"A"}
def emit(name, segs, links, paths, sample):
    d = dict(segs)
    with open("%s/%s.gfa" % (out, name), "w") as g:
        g.write("H\tVN:Z:1.0\n")
        for n, x in segs: g.write("S\t%s\t%s\n" % (n, x))
        for a, b in links: g.write("L\t%s\t+\t%s\t+\t0M\n" % (a, b))
        for n, st in paths.items(): g.write("P\t%s\t%s\t*\n" % (n, ",".join(x + "+" for x in st)))
    sq = {n: "".join(d[x] for x in st) for n, st in paths.items()}
    with open("%s/%s.r1.fq" % (out, name), "w") as f1, open("%s/%s.r2.fq" % (out, name), "w") as f2:
        k = 0
        for nm in sample:
            h = sq[nm]
            for i in range(0, len(h) - 360, 6):
                a = h[i:i+150]; b = h[i+200:i+350]
                if len(b) < 150: break
                rc = "".join(comp[c] for c in reversed(b))
                f1.write("@%s_%d/1\n%s\n+\n%s\n" % (nm, k, a, "I"*150))
                f2.write("@%s_%d/2\n%s\n+\n%s\n" % (nm, k, rc, "I"*150))
                k += 1
X = seq(4000,1); A1 = seq(600,2); B1 = mut(A1,90,3)
A2 = seq(600,5); B2 = mut(A2,90,6); Z = seq(4000,7)
S1 = seq(20,10); MID = seq(80,11); S2 = seq(20,12)
emit("bb",
     [("1",X),("2",A1),("3",B1),("4",S1),("5",MID),("6",S2),("7",A2),("8",B2),("9",Z)],
     [("1","2"),("1","3"),("2","4"),("3","4"),("4","5"),("5","6"),("6","7"),("6","8"),("7","9"),("8","9")],
     {"hapAA":["1","2","4","5","6","7","9"],"hapAB":["1","2","4","5","6","8","9"],
      "hapBA":["1","3","4","5","6","7","9"],"hapBB":["1","3","4","5","6","8","9"]},
     ("hapAB","hapBA"))
M1 = seq(80,20); M2 = mut(M1,20,21)
emit("vm",
     [("1",X),("2",A1),("3",B1),("4",S1),("5",M1),("6",M2),("7",S2),("8",A2),("9",B2),("10",Z)],
     [("1","2"),("1","3"),("2","4"),("3","4"),("4","5"),("4","6"),("5","7"),("6","7"),
      ("7","8"),("7","9"),("8","10"),("9","10")],
     {"hAMA":["1","2","4","5","7","8","10"],"hAMB":["1","2","4","5","7","9","10"],
      "hBNA":["1","3","4","6","7","8","10"],"hBNB":["1","3","4","6","7","9","10"]},
     ("hAMB","hBNA"))
PYEOF2
for FX in bb vm; do
  REF=$(grep "^P" "$OUT/$FX.gfa" | head -1 | cut -f2)
  "$BIN" bubble -i "$OUT/$FX.gfa" -r "$REF" -o "$OUT/$FX.b" --min-variant-bp 0 -q >/dev/null 2>&1
  "$BIN" genotype-frag -i "$OUT/$FX.b.sorted.gfa" -b "$OUT/$FX.b" -o "$OUT/$FX.o" \
    -R "$OUT/$FX.r1.fq" -R "$OUT/$FX.r2.fq" -t 2 --max-divergence 0.05 --fragment-len 350 \
    --fragment-sd 50 --error-rate 0.001 --ownership-table "$OUT/$FX.own.tsv" -q >/dev/null 2>&1
  [ -s "$OUT/$FX.own.tsv" ] || { bad "$FX: no ownership table"; continue; }
  "$PY" - "$OUT/$FX.own.tsv" "$OUT/$FX.own.tsv.blocks.tsv" "$FX" <<'PYEOF3'
import sys, collections
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])][1:]
blk  = [l.rstrip("\n").split("\t") for l in open(sys.argv[2])][1:]
fx = sys.argv[3]; bad = 0
def ok(m): print("  ok   " + m)
def no(m):
    global bad; bad += 1; print("  FAIL " + m)
var = {int(b[0]): b[4] == "1" for b in blk}
nvar = sum(1 for v in var.values() if v)
kinds = collections.Counter(r[1] for r in rows)
print("       %s: %d of %d blocks variable; %s" %
      (fx, nvar, len(var), ", ".join("%s=%d" % kv for kv in sorted(kinds.items()))))
# A fragment physically spanning the two outer bubbles.
span = [r for r in rows if int(r[6]) <= 1 and int(r[7]) >= 3]
if not span: no("%s: no fragment spans the outer variable blocks" % fx)
elif fx == "bb":
    # GATE 5: the intervening backbone is FIXED, so these are PAIRWISE, not wide.
    if not all(v is False for b, v in var.items() if b == 2):
        no("bb: block 2 should be a FIXED backbone but is marked variable")
    elif any(r[1] != "linkage" for r in span):
        no("bb: %d/%d backbone-spanning fragments are not pairwise (%s) -- linkage evidence lost"
           % (sum(1 for r in span if r[1] != "linkage"), len(span),
              collections.Counter(r[1] for r in span)))
    else:
        ok("bb: %d fragments across a FIXED backbone form a PAIRWISE factor psi(%s,%s)"
           % (len(span), span[0][2], span[0][3]))
        ok("bb: the fixed backbone is carried as context, costing no arity")
else:
    # GATE 6: the intervening block VARIES, so these must be three-variable, never cropped.
    if any(r[1] != "wide" or int(r[5]) < 3 for r in span):
        no("vm: a fragment over a VARIABLE middle was cropped to %s"
           % collections.Counter(r[1] for r in span))
    else:
        ok("vm: %d fragments over a VARIABLE middle give three-variable factors, not cropped pairs"
           % len(span))
    if kinds.get("linkage", 0) == 0:
        no("vm: adjacent variable pairs produced no pairwise factor at all")
    else:
        ok("vm: adjacent variable pairs still form pairwise factors (%d)" % kinds["linkage"])
seen = collections.Counter(r[0] for r in rows)
if any(n != 1 for n in seen.values()): no("%s: ownership is not a partition" % fx)
else: ok("%s: ownership is a partition over %d fragments" % (fx, len(seen)))
sys.exit(bad)
PYEOF3
  fails=$(( fails + $? ))
done

# ---------------------------------------------------------------------------------------------
# GATE 7: THE EXPOSURE PRECONDITION. exposure(n) = SUM_L pi(L) max(0, n-L+1) equals the affine
# surrogate n+1-E[L] only once the window clears the insert support; the clip bites below that and
# the cis/trans cancellation goes with it. Asserted at the boundary rather than trusted.
"$BIN" genotype-frag -i "$OUT/b.sorted.gfa" -b "$OUT/b" -o "$OUT/e" -R "$OUT/r1.fq" -R "$OUT/r2.fq" \
  --fragment-len 350 --fragment-sd 50 --error-rate 0.001 \
  --exposure-probe 100,300,548,549,550,4000 -q > "$OUT/expo.tsv" 2>/dev/null
if [ -s "$OUT/expo.tsv" ]; then
  "$PY" - "$OUT/expo.tsv" <<'PYEOF4'
import sys
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])][1:]
d = {int(r[0]): (float(r[1]), float(r[2]), float(r[3]), r[4] == "1", int(r[6])) for r in rows}
bad = 0
def ok(m): print("  ok   " + m)
def no(m):
    global bad; bad += 1; print("  FAIL " + m)
hi = next(iter(d.values()))[4]
inside  = [n for n, v in d.items() if v[3]]
outside = [n for n, v in d.items() if not v[3]]
if not inside or not outside:
    no("the probe does not straddle the regime boundary (insert_hi=%d)" % hi)
else:
    ok("probe straddles the regime boundary at insert_hi=%d: in %s, out %s"
       % (hi, sorted(inside), sorted(outside)))
# Inside the regime the two forms must agree to floating point.
worst_in = max(abs(d[n][2]) for n in inside)
if worst_in < 1e-6: ok("inside the regime exact and affine agree (worst %.2e nats)" % worst_in)
else: no("inside the regime they differ by %.4g nats -- cancellation is not exact" % worst_in)
# Outside it they must NOT: a gate that passes because the difference is negligible everywhere
# would be asserting nothing.
worst_out = max(abs(d[n][2]) for n in outside)
if worst_out > 1.0:
    ok("outside the regime they diverge (worst %.4g nats) -- the precondition is real" % worst_out)
else:
    no("outside the regime they differ by only %.4g nats -- gate is vacuous" % worst_out)
# The boundary must sit exactly at hi-1, not near it.
if d.get(hi - 1, (0,0,0,False,0))[3] and not d.get(hi - 2, (0,0,0,True,0))[3]:
    ok("the regime begins exactly at hi-1 = %d" % (hi - 1))
else:
    no("the regime boundary is not at hi-1 = %d" % (hi - 1))
sys.exit(bad)
PYEOF4
  fails=$(( fails + $? ))
else
  bad "exposure probe produced nothing"
fi

# GATE 8: THE DISCARDED-CONTENT TRADEOFF, quantified. Linkage-owned fragments leave the marker
# unaries; under a factor conditional on endpoint content only their phase information survives.
# This is not lossless, so its size is recorded here and block-content calls must not regress.
if [ -s "$OUT/own.tsv.ledger.tsv" ]; then
  "$PY" - "$OUT/own.tsv.ledger.tsv" <<'PYEOF5'
import sys
d = dict(l.rstrip("\n").split("\t") for l in open(sys.argv[1]))
bad = 0
def ok(m): print("  ok   " + m)
def no(m):
    global bad; bad += 1; print("  FAIL " + m)
n = int(d["fragments"]); link = int(d["linkage"])
fr = float(d["linkage_fragment_share"])
CLASSES = ("unary", "linkage", "wide", "invariant", "unusable")
tot = sum(int(d[k]) for k in CLASSES)
if tot == n: ok("the ledger's classes sum to every fragment (%d)" % n)
else: no("ledger classes sum to %d, %d fragments loaded" % (tot, n))
# EVERY class must be reported, so none can go unaccounted when the chain decides its disposition.
missing = [c for c in CLASSES if "%s_in_band_mass_share" % c not in d]
if missing: no("no in-band mass share reported for: %s" % ", ".join(missing))
else:
    ok("every ownership class reports its share: " +
       ", ".join("%s %.2f%%" % (c, 100 * float(d["%s_in_band_mass_share" % c])) for c in CLASSES))
shares = sum(float(d["%s_in_band_mass_share" % c]) for c in CLASSES)
if abs(shares - 1.0) < 1e-6: ok("the class mass shares account for all in-band mass (sum %.6f)" % shares)
else: no("class mass shares sum to %.6f, not 1 -- mass is unaccounted" % shares)
# NON-VACUITY, not an acceptance cutoff. An arbitrary "< 25%" threshold would be a tuned pass
# condition on a fixture; what matters is that the gate measures something real and that the
# BLOCK-CONTENT regression -- C4 and the exact leave-zero-out controls -- decides whether the
# deliberate content loss is acceptable. The mass share is a size statistic and cannot decide it:
# a fragment holding a negligible share can still carry a decisive likelihood ratio.
if link > 0 and fr > 0:
    ok("linkage-owned: %d of %d fragments (%.2f%%), holding %.2f%% of pooled in-band mass"
       % (link, n, 100 * fr, 100 * float(d["linkage_in_band_mass_share"])))
    ok("      (a SIZE statistic -- only the block-content regression measures information lost)")
else:
    no("no linkage fragments -- the tradeoff gate measures nothing")
sys.exit(bad)
PYEOF5
  fails=$(( fails + $? ))
else
  bad "no ownership ledger written"
fi

if ! "$BIN" genotype-frag --help 2>&1 | grep -q -- "--hybrid-call"; then
  echo
  echo "  PENDING, held behind --hybrid-call (the chain that forms the mixtures does not exist):"
  echo "    * hybrid DISABLED -> exactly the legacy marker / Li-Stephens result;"
  echo "    * hybrid ENABLED but NO linkage-owned fragments -> exactly the legacy result;"
  echo "    * hybrid ENABLED with linkage fragments -> content calls pass the block-content"
  echo "      regression. Reproducing the legacy result here is NOT the gate and cannot be:"
  echo "      those fragments' content evidence was deliberately removed from the unaries.;"
  echo "    * with linkage, correct phase while preserving unordered block content;"
  echo "    * the background RETAINED inside each configuration's mixture;"
  echo "    * Li-Stephens prior and fragment linkage each applied exactly once;"
  echo "    * a global swap of the two homologues leaves the output unchanged;"
  echo "    * Wide-owned fragments are reported UNSUPPORTED/INCOMPLETE, never silently dropped;"
  echo "    * ambiguous evidence gives an equivalence set or UNRESOLVED, not a confident guess."
  echo "  This file must FAIL, not skip, once --hybrid-call appears without satisfying them."
else
  # THE SENTINEL. Without this branch the paragraph above is a promise nothing enforces: the moment
  # --hybrid-call exists the `if` is simply skipped and the file exits 0 having tested none of the
  # chain. That is the vacuous-gate pattern this suite exists to prevent, so the arrival of the
  # entry point FAILS here until the assertions above replace this block.
  bad "--hybrid-call exists but none of the chain gates are implemented."
  echo "       Replace this sentinel with the real assertions. Until then its presence is the"
  echo "       failure: an entry point that no gate constrains is worse than no entry point."
fi
echo
if [ "$fails" -eq 0 ]; then echo "hybrid chain: all active assertions passed"; else
  echo "hybrid chain: $fails assertion(s) failed"; fi
exit "$fails"
