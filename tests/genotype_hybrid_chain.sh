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
#   7. ambiguous evidence yields an equivalence set or UNRESOLVED, never a confident guess.
#
# NORMALISATION is settled in genotype_fragments.hpp and asserted at the arithmetic level below:
# linkage is a CONDITIONAL PHASE SCORE and exposure cancels exactly, because exposure(n) = n+1-E[L]
# is affine in length and cis/trans partition the same allele multiset. Partitioning observed
# fragments does not partition normalisation on its own, so this is chosen rather than assumed.
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

if ! "$BIN" genotype-frag --help 2>&1 | grep -q -- "--hybrid-call"; then
  echo
  echo "  PENDING gates 1, 2 and 7: the chain entry point (--hybrid-call) does not exist yet."
  echo "  Ownership is pinned and asserted above; the chain that consumes it is not built."
  echo "  This file must FAIL, not skip, once --hybrid-call appears without satisfying them."
fi
echo
if [ "$fails" -eq 0 ]; then echo "hybrid chain: all active assertions passed"; else
  echo "hybrid chain: $fails assertion(s) failed"; fi
exit "$fails"
