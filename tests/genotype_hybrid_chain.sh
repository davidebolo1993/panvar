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
#   9. the linkage EMISSION keeps its background inside the mixture, exposure exact;      [ACTIVE]
#  10. psi is formed by edge-level diploid aggregation, sums to 1 per content class, and is
#      invariant under a global homologue swap;                                           [ACTIVE]
#  11. a three-allele fixture pins the indexing a 0/1 fixture cannot reach;                [ACTIVE]
#  12. an edge carrying NO phase information is exactly neutral: zero fragments, flat emissions
#      or all-unplaced emissions each give an identically zero log factor;                  [ACTIVE]
#  13. exposure that does not cancel, an unformable emission, and an oversized dense table are each
#      REFUSED with a status, never truncated or silently skipped;                          [ACTIVE]
#  14. forward-backward agrees with a BRUTE-FORCE path oracle on the log partition and on every
#      block marginal, with a linkage-free control and a non-vacuity check;                 [ACTIVE]
#  15. ONE chain exercises BOTH kernel edge paths -- factorised and linked -- so neither branch can
#      be present but unused;                                                               [ACTIVE]
#  17. marker occurrence exclusion subtracts a linkage fragment's OCCURRENCES and never deletes a
#      marker shared with unary-owned fragments;                                            [ACTIVE]
#  24. an EXACT reverse-complement duplicate is indistinguishable in CHAIN orientation -- same
#      sequence, same block spans -- while still flagged antiparallel, and a mixed forward/reverse
#      panel builds geometry and recovers phase;                                            [ACTIVE]
#  22. every state in the marker HMM's DECLARED universe has a verified candidate frame, with the
#      universe enumerated by name and checked UNIQUE before comparison; a verified partial terminal
#      frame is ACCEPTED and reported separately, never counted as missing; a shortfall is an
#      INCOMPLETE hybrid-SUBSTRATE result labelled reason=candidate-frame-coverage, never biological
#      ambiguity. Measured on real panels: C4 131 complete, CYP2D6 126 complete + 1 accepted partial
#      (NA18989#1#haplotype1), LPA 466 complete; 0 missing everywhere;                      [ACTIVE]
#  21. the four status conditions stay distinct -- ownership_complete, factors_buildable,
#      hybrid_activated -- and call_status is COMPLETE only when all three hold. A run can be
#      ownership-complete and factor-INcomplete, and must then report INCOMPLETE;           [ACTIVE]
#
# EXIT-CODE POLICY, so a pipeline can tell the two apart: a model-level INCOMPLETE is a VALID program
# result and exits zero with an explicit status. Malformed input or internal inconsistency exits
# nonzero. "The model cannot represent this evidence" is not "the command failed".
#
#  20. activation is TRANSACTIONAL: nothing is subtracted from the marker unaries unless every
#      required edge was built and the model is complete, and
#          {excluded fragments} == {fragments consumed by ACTIVE linkage edges}, each once;  [ACTIVE]
#  19. hybrid COMPLETE <=> every non-invariant owned fragment has a supported consumer; a refused
#      edge or Wide fragment makes the run INCOMPLETE (not "unresolved"), all reasons are retained,
#      and no refused edge ever reaches the kernel;                                         [ACTIVE]
#  18. the haplotype -> allele mapping validates every value BEFORE the int -> unsigned conversion:
#      permuted and normal mappings resolve, a bypassing haplotype resolves to bypass_allele, and a
#      missing or out-of-range value is refused with nothing left indexable;                [ACTIVE]
#
# STILL OWED AT STEP 4, and easy to lose: under --hybrid-call the exclusion list must be GENERATED
# from the same ownership ledger that builds the linkage factors, and the identity
#     excluded fragment IDs == linkage-owned fragment IDs
# asserted end to end. Ownership and counting can each be correct while disagreeing about which
# fragments they cover, and --exclude-fragments today only validates the counting mechanism.
#  29. the allele-product SUPPORT SEARCH equals the dense oracle cell for cell -- finite cells, log
#      mass and informative classification -- including an off-panel (alpha,beta) combination no
#      panel haplotype carries, free-dimension expansion, repeat multiplicity, the exhaustive
#      non-ACGT fallback and zero-state cells;                                              [ACTIVE]
#  28. the DENSE-EMISSION window guard (dense_emission_window_alignments = fragments x n_A x n_B
#      summed over edges) is operational, not statistical, and refuses transactionally with its own
#      reason. Named for the construction it bounds, so the forthcoming support search reports its
#      own counters rather than reusing this budget;                                        [ACTIVE]
#  27. the GROUPED sparse contraction equals the dense kernel on weight sums and marginals, an edge
#      with no classes takes the factorised path unchanged, and an adversarial case placing the
#      dominant Li-Stephens mass in the phase psi drives to zero still agrees;             [ACTIVE]
#  26. sparse psi equals dense psi on EVERY configuration of every densely-testable edge; a class
#      with no in-band corner is proved exactly neutral and never stored; and a C4-scale edge that
#      the dense table refuses (197,177,764 configurations) completes.                      [ACTIVE]
#  16. every generated linkage edge is finite, mean-one per content class, and bounded by
#      max log psi <= log|C| <= log 4, so exp(log psi) cannot overflow. An earlier version of this
#      contract claimed exp would OVERFLOW at the 16043-nat spread; that spread is entirely on the
#      NEGATIVE side (losing configurations underflowing), and the claim was wrong.         [ACTIVE]
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
# GATE 9: THE LINKAGE POTENTIAL psi_f, gated before any chain consumes it. Two properties:
# the background stays INSIDE each configuration's mixture, and exposure is computed EXACTLY per
# window rather than assumed to cancel.
"$BIN" genotype-frag -i "$OUT/b.sorted.gfa" -b "$OUT/b" -o "$OUT/lp" -R "$OUT/r1.fq" -R "$OUT/r2.fq" \
  -t 2 --max-divergence 0.05 --fragment-len 350 --fragment-sd 50 --error-rate 0.001 \
  --ownership-table "$OUT/own.tsv" --linkage-potential "$OUT/lpot.tsv" -q >/dev/null 2>&1
if [ -s "$OUT/lpot.tsv" ]; then
  "$PY" - "$OUT/lpot.tsv" <<'PYEOF6'
import sys, math, collections
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])][1:]
bad = 0
def ok(m): print("  ok   " + m)
def no(m):
    global bad; bad += 1; print("  FAIL " + m)
def num(x):
    try: return float(x)
    except ValueError: return float("nan")
if not rows: no("no linkage potentials emitted"); sys.exit(1)
ok("%d linkage potentials emitted, 0 refused" % len(rows))
with_bg = [num(r[7]) for r in rows]
no_bg   = [num(r[8]) for r in rows]
# THE BACKGROUND IS LOAD-BEARING. With it inside the mixture every contrast is finite; without it a
# fragment that fails to place on a configuration gives an INFINITE log-ratio and can veto that
# configuration single-handedly. This is the property, not a rounding detail.
fin_w = [x for x in with_bg if math.isfinite(x)]
if len(fin_w) == len(rows): ok("every contrast is finite with the background inside the mixture")
else: no("%d contrasts are non-finite WITH the background -- the floor is not bounding them"
         % (len(rows) - len(fin_w)))
nonfin_n = [x for x in no_bg if not math.isfinite(x)]
if nonfin_n:
    ok("dropping the background makes %d of %d contrasts INFINITE (%.0f%%) -- each would veto a "
       "configuration outright" % (len(nonfin_n), len(rows), 100.0*len(nonfin_n)/len(rows)))
else:
    no("dropping the background changes no contrast to infinite -- the gate is vacuous here, so it "
       "is not demonstrating that the background must be retained")
# Phase must actually be decided: all junction fragments should agree on one configuration, and by
# a margin, or there is no linkage signal to carry into the chain.
# PHASE-INFORMATIVE means a NON-ZERO contrast. A fragment can span the junction and still carry no
# phase information -- if it reaches a distinguishing position in only one of the two blocks, its
# mass is identical under both configurations and the contrast is exactly 0. Those fragments are
# legitimately edge-owned (they still consume normalisation) but must not be counted as evidence,
# and requiring them to "agree" would be asserting a signal that cannot exist.
sgn = collections.Counter("neg" if x < 0 else ("pos" if x > 0 else "zero") for x in fin_w)
inform = [x for x in fin_w if x != 0.0]
if not inform:
    no("no finite contrast is non-zero -- the junction carries no phase signal at all")
elif len(set(x < 0 for x in inform)) == 1:
    ok("all %d PHASE-INFORMATIVE contrasts agree on one configuration, mean |contrast| %.1f nats"
       % (len(inform), sum(abs(x) for x in inform)/len(inform)))
else:
    no("phase-informative fragments disagree on the configuration: %s" % dict(sgn))
# OVER-ASSIGNMENT, measured. Ownership scopes by the blocks an origin SPANS, which is a superset of
# the blocks its mass actually depends on. Fragments in that gap are pulled out of the marker
# unaries -- losing their content evidence -- and give no phase evidence in return. Recorded so the
# cost of the span-union rule is visible rather than assumed to be zero.
z = sgn.get("zero", 0)
print("       %d of %d edge-owned fragments are phase-UNINFORMATIVE (%.0f%%): edge-owned by span, "
      "but their mass does not depend on the combination" % (z, len(rows), 100.0*z/len(rows)))
# Exposure: recorded per window, exact by construction. On THIS fixture every window clears the
# insert support, so the short-window branch is not exercised here -- gate 7 covers the boundary.
aff = collections.Counter(r[6] for r in rows)
wmin = min(int(r[5]) for r in rows)
ok("windows are %d bp and up; exposure regime affine=%s (short-window branch is gate 7's job)"
   % (wmin, dict(aff)))
sys.exit(bad)
PYEOF6
  fails=$(( fails + $? ))
else
  bad "no linkage potentials written"
fi

# GATE 10: EDGE-LEVEL AGGREGATION. psi is formed only AFTER the two homologues are combined once,
# the background mixed once, every owned fragment summed, and exposure charged once for the edge.
# Normalising per fragment would let each fragment pick its own phase configuration -- the mosaic
# error excluded at the diplotype level, reappearing one level down.
if [ -s "$OUT/lpot.tsv.edges.tsv" ]; then
  "$PY" - "$OUT/lpot.tsv.edges.tsv" "$OUT/lpot.tsv" <<'PYEOF7'
import sys
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])][1:]
frs  = [l.rstrip("\n").split("\t") for l in open(sys.argv[2])][1:]
bad = 0
def ok(m): print("  ok   " + m)
def no(m):
    global bad; bad += 1; print("  FAIL " + m)
if not rows: no("no linkage edge aggregated"); sys.exit(1)
for r in rows:
    na, nb = int(r[2]), int(r[3])
    nfrag, ninf, ninv, ncfg = int(r[4]), int(r[5]), int(r[6]), int(r[7])
    psisum, spread, swap = float(r[9]), float(r[11]), float(r[12])
    usable, status = r[13] == "1", r[14]
    if not usable:
        no("edge %s-%s is %s -- it must not be consumed" % (r[0], r[1], status)); continue
    if ninv: no("edge %s-%s has %d owned fragments with no emission yet is usable" % (r[0],r[1],ninv))
    tag = "edge %s-%s" % (r[0], r[1])
    if ncfg == na*nb*na*nb: ok("%s: %d configurations over %dx%d alleles" % (tag, ncfg, na, nb))
    else: no("%s: %d configurations, expected %d" % (tag, ncfg, na*nb*na*nb))
    # psi is a MEAN-ONE likelihood ratio within each unordered-content class -- NOT a sum-one
    # conditional distribution. It multiplies into the Li-Stephens transition, which already carries
    # a phase prior, and sum-one would both double that normalisation and leave a -log|C| content
    # penalty on an edge that says nothing (classes have sizes 1, 2 and 4).
    if psisum < 1e-9: ok("%s: psi is mean-one within every content class (worst dev %.2e)" % (tag, psisum))
    else: no("%s: psi mean deviates from 1 by %.2e within a content class" % (tag, psisum))
    # GLOBAL HOMOLOGUE SWAP: (a1,b1,a2,b2) and (a2,b2,a1,b1) are one diploid state written twice.
    if swap < 1e-9: ok("%s: global homologue swap leaves psi unchanged (%.2e)" % (tag, swap))
    else: no("%s: swapping the homologues changes psi by %.4g -- ordered state leaked in" % (tag, swap))
    if spread > 0.0: ok("%s: phase is decided, spread %.1f nats within a content class" % (tag, spread))
    else: no("%s: psi is flat -- the edge decides no phase at all" % tag)
    # BOTH COUNTS. Edge-owned is what leaves the unaries; phase-informative is what actually enters
    # psi. Reporting only the first overstates the evidence the factor carries.
    ok("%s: %d fragments edge-owned, %d phase-informative (%.0f%%)"
       % (tag, nfrag, ninf, 100.0*ninf/nfrag if nfrag else 0.0))
# The per-fragment informative flag must agree with the edge's count -- two paths, one answer.
finf = sum(1 for r in frs if r[10] == "1")
einf = sum(int(r[5]) for r in rows if r[13] == "1")
if finf == einf: ok("per-fragment and edge-level informative counts agree (%d)" % finf)
else: no("per-fragment informative %d, edge-level %d" % (finf, einf))
sys.exit(bad)
PYEOF7
  fails=$(( fails + $? ))
else
  bad "no edge aggregate written"
fi

# ---------------------------------------------------------------------------------------------
# GATE 25: THE DERIVED INVARIANT FLANK. A pairwise A-B factor may include only sequence that is
# invariant OUTSIDE A and B. A fixed reach borrows the neighbouring variable block, and the geometry
# then correctly refuses -- which is what refused every C4 edge once orientation was fixed. Each side
# is bounded by BOTH the nearest external VARIABLE boundary and the COMMON VERIFIED MAPPED boundary
# (a partial terminal frame verifies less than the catalogue holds), then capped by the request.
#
# THE FIXTURE has three bubbles separated by GENUINE backbone blocks, so invariant sequence exists
# outside the scored pair -- every earlier fixture had its variable blocks spanning to the chain
# ends, where the derived flank is zero for a structural reason and proves nothing about the cap.
"$PY" - "$OUT" <<'PYEOFL'
import sys, random
out = sys.argv[1]; B = "ACGT"
def seq(n, s):
    r = random.Random(s); return "".join(r.choice(B) for _ in range(n))
def mut(s, k, sd):
    r = random.Random(sd); s = list(s)
    for p in r.sample(range(len(s)), k): s[p] = r.choice([c for c in B if c != s[p]])
    return "".join(s)
X = seq(1500,1); A = seq(400,2); Bb = mut(A,60,3)
S1 = seq(20,11); M1 = seq(150,4); S2 = seq(20,12)
C = seq(400,5); D = mut(C,60,6)
S3 = seq(20,13); M2 = seq(300,7); S4 = seq(20,14)
E = seq(400,8); F = mut(E,60,9); Z = seq(1500,10)
segs = [("1",X),("2",A),("3",Bb),("4",S1),("5",M1),("6",S2),("7",C),("8",D),
        ("9",S3),("10",M2),("11",S4),("12",E),("13",F),("14",Z)]
links = [("1","2"),("1","3"),("2","4"),("3","4"),("4","5"),("5","6"),("6","7"),("6","8"),
         ("7","9"),("8","9"),("9","10"),("10","11"),("11","12"),("11","13"),("12","14"),("13","14")]
paths = {"h1":["1","2","4","5","6","7","9","10","11","12","14"],
         "h2":["1","3","4","5","6","8","9","10","11","13","14"],
         "h3":["1","2","4","5","6","8","9","10","11","12","14"],
         "h4":["1","3","4","5","6","7","9","10","11","13","14"]}
d = dict(segs)
for tag, extra in (("fl", None), ("fx", "13")):
    with open("%s/%s.gfa" % (out, tag), "w") as g:
        g.write("H\tVN:Z:1.0\n")
        for n, x in segs:
            if extra is not None and n == extra:
                x = "".join(("A" if c != "A" else "T") if i % 7 == 0 else c
                            for i, c in enumerate(x))
            g.write("S\t%s\t%s\n" % (n, x))
        for a, b in links: g.write("L\t%s\t+\t%s\t+\t0M\n" % (a, b))
        for n, st in paths.items(): g.write("P\t%s\t%s\t*\n" % (n, ",".join(y+"+" for y in st)))
comp = {"A":"T","C":"G","G":"C","T":"A"}
h = "".join(d[x] for x in paths["h1"])
k = 0
with open(out+"/fl.r1.fq","w") as f1, open(out+"/fl.r2.fq","w") as f2:
    for i in range(0, len(h)-360, 8):
        a = h[i:i+150]; b = h[i+200:i+350]
        if len(b) < 150: break
        rc = "".join(comp[c] for c in reversed(b))
        f1.write("@f%d/1\n%s\n+\n%s\n" % (k, a, "I"*150))
        f2.write("@f%d/2\n%s\n+\n%s\n" % (k, rc, "I"*150)); k += 1
open(out+"/fl.sizes.txt","w").write("%d %d" % (len(M1), len(M2)))
PYEOFL
for TAG in fl fx; do
  "$BIN" bubble -i "$OUT/$TAG.gfa" -r h1 -o "$OUT/$TAG.b" --min-variant-bp 0 -q >/dev/null 2>&1
  "$BIN" genotype -i "$OUT/$TAG.b.sorted.gfa" -b "$OUT/$TAG.b" -r h1 -o "$OUT/$TAG.g" \
    -R "$OUT/fl.r1.fq" -R "$OUT/fl.r2.fq" --fragment-len 350 \
    --hybrid-geometry-probe "$OUT/$TAG.geom.tsv" -q >/dev/null 2>&1
done
if [ -s "$OUT/fl.geom.tsv" ] && [ -s "$OUT/fx.geom.tsv" ]; then
  "$PY" - "$OUT/fl.geom.tsv" "$OUT/fx.geom.tsv" "$OUT/fl.sizes.txt" <<'PYEOFM'
import sys
def load(p):
    return [l.rstrip("\n").split("\t") for l in open(p)][1:]
A, Bx = load(sys.argv[1]), load(sys.argv[2])
m1, m2 = (int(x) for x in open(sys.argv[3]).read().split())
bad = 0
def ok(m): print("  ok   " + m)
def no(m):
    global bad; bad += 1; print("  FAIL " + m)
def row(rows, a, b, fl):
    for r in rows:
        if int(r[0]) == a and int(r[1]) == b and int(r[6]) == fl: return r
    return None
if any(r[9] != "1" for r in A):
    no("some edge refuses geometry: %s" % [(r[0], r[1], r[6], r[10]) for r in A if r[9] != "1"])
else:
    ok("every adjacent-variable edge builds geometry at all requested widths")
# EXTERNAL invariant block contributes real context, bounded by the next variable.
r = row(A, 1, 3, 550)
if r and int(r[8]) == m2: ok("external invariant block contributes exactly %d bp of right flank" % m2)
else: no("edge 1-3 right flank is %s, expected %d" % (r[8] if r else "?", m2))
r = row(A, 3, 5, 550)
if r and int(r[7]) == m1: ok("external invariant block contributes exactly %d bp of left flank" % m1)
else: no("edge 3-5 left flank is %s, expected %d" % (r[7] if r else "?", m1))
# ADJACENT to the chain end -> derived zero, for a structural reason.
r = row(A, 1, 3, 550)
if r and int(r[7]) == 0: ok("no invariant sequence available on the outer side -> derived flank 0")
else: no("edge 1-3 left flank is %s, expected 0" % (r[7] if r else "?"))
# REQUEST SMALLER THAN AVAILABLE -> exactly the request.
r = row(A, 1, 3, 100)
if r and int(r[8]) == 100: ok("a request of 100 below the available %d yields exactly 100" % m2)
else: no("capped right flank is %s, expected 100" % (r[8] if r else "?"))
r = row(A, 3, 5, 100)
if r and int(r[7]) == 100: ok("a request of 100 below the available %d yields exactly 100" % m1)
else: no("capped left flank is %s, expected 100" % (r[7] if r else "?"))
# CHANGING AN EXTERNAL VARIABLE ALLELE must not touch A-B geometry.
p, q = row(A, 1, 3, 550), row(Bx, 1, 3, 550)
if p and q and p[7:10] == q[7:10]:
    ok("mutating an EXTERNAL variable allele leaves edge 1-3 geometry unchanged (%s/%s)"
       % (p[7], p[8]))
else:
    no("external allele change altered edge 1-3 geometry: %s vs %s"
       % (p[7:10] if p else "?", q[7:10] if q else "?"))
sys.exit(bad)
PYEOFM
  fails=$(( fails + $? ))
else
  bad "the derived-flank fixture produced no geometry probe"
fi
# ---------------------------------------------------------------------------------------------
# GATE 29c: EXACT SIGNATURE GROUPING -- member-weighted HMM equivalence, not just factor values.
# Two allele pairs may share a linkage contribution only when the whole vector of per-fragment
# structural signatures agrees. This checks that the grouped representation reproduces the
# UNCOMPRESSED one through psi, the partition weight and every marginal, under adversarial weights:
# unequal marker unaries among alleles sharing a class, unequal class sizes, r at both extremes.
# Three mutations must FAIL -- dropping mismatch counts, merging straight with crossed, and
# ignoring membership -- or the assertions above are measuring nothing.
if "$BIN" genotype-frag -i /dev/null -b none -o "$OUT/gr" --grouping-selftest \
     > "$OUT/group.txt" 2>/dev/null; then
  while IFS=$'\t' read -r v m; do
    [ "$v" = ok ] && ok "$m" || { [ -n "${m:-}" ] && bad "$m"; }
  done < <(grep -E '^(ok|FAIL)\t' "$OUT/group.txt")
else
  bad "grouping selftest reported failures"
  sed -n 's/^FAIL\t/  /p' "$OUT/group.txt"
fi

# GATE 29b: THE OPERATIONAL WORK BUDGET, at unit level. The production gate above proves the
# refusal is transactional; this proves it is charged BEFORE the work, that an unlimited budget is
# bit-inert, and that an unaffordable dense fallback never enumerates a window.
if "$BIN" genotype-frag -i /dev/null -b none -o "$OUT/bg" --budget-selftest \
     > "$OUT/budget.txt" 2>/dev/null; then
  while IFS=$'\t' read -r v m; do
    [ "$v" = ok ] && ok "$m" || { [ -n "${m:-}" ] && bad "$m"; }
  done < <(grep -E '^(ok|FAIL)\t' "$OUT/budget.txt")
else
  bad "budget selftest reported failures"
  sed -n 's/^FAIL\t/  /p' "$OUT/budget.txt"
fi

# GATE 29a: COORDINATE VALIDATION. Direct verification cannot recover from a wrong derived start:
# it does not fail loudly, it verifies the wrong bases and silently drops the placement. So the
# coordinate is checked on its own -- the view against the materialised window, the |A2|-|A1| shift
# of a B-side seed, and an EMPTY allele, which the earlier six-map index could not seed at all.
if "$BIN" genotype-frag -i /dev/null -b none -o "$OUT/co" --coordinate-selftest \
     > "$OUT/coord.txt" 2>/dev/null; then
  while IFS=$'\t' read -r v m; do
    [ "$v" = ok ] && ok "$m" || { [ -n "${m:-}" ] && bad "$m"; }
  done < <(grep -E '^(ok|FAIL)\t' "$OUT/coord.txt")
else
  bad "coordinate selftest reported failures"
  sed -n 's/^FAIL\t/  /p' "$OUT/coord.txt"
fi

# GATE 29: THE ALLELE-PRODUCT SUPPORT SEARCH must equal the DENSE ORACLE cell for cell -- same
# finite cells, same log mass, same informative classification. "The same phase call" would pass
# while multiplicity or an off-panel combination went missing.
#
# THE DECISIVE CASE is an (alpha, beta) combination NO panel haplotype carries. Ownership placements
# are placements on complete panel candidates, so deciding support from them would collapse the
# hybrid back toward complete-panel haplotypes -- and would pass every other gate here.
"$BIN" genotype-frag -i /dev/null -b none -o "$OUT/ss" --support-selftest \
  > "$OUT/support.tsv" 2>/dev/null
if [ -s "$OUT/support.tsv" ]; then
  "$PY" - "$OUT/support.tsv" <<'PYEOFS'
import sys
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])]
h, rows = rows[0], rows[1:]
d = {r[0]: dict(zip(h, r)) for r in rows if len(r) == len(h)}
bad = 0
def ok(m): print("  ok   " + m)
def no(m):
    global bad; bad += 1; print("  FAIL " + m)
need = ("offpanel_A2_B1", "junction_crossing_seed", "seed_spans_A_ctx_B", "A_only_constraint",
        "B_only_constraint", "invariant_only_seed", "reverse_strand", "nine_identical_origins",
        "non_acgt_fallback", "zero_state_cells", "c4_scale_118x119", "c4_scale_dup_alleles")
miss = [k for k in need if k not in d]
if miss: no("support cases missing: %s" % ", ".join(miss)); sys.exit(1)
ok("support search covered over %d cases" % len(need))
# STAGE COUNTERS must be populated, or a stage is silently not running.
st = d["c4_scale_118x119"]
stages = ("seed_start_proposals", "unique_seed_starts", "seed_compatible_joins",
          "full_read_verifications", "accepted_mate_placements", "verified_fr_states",
          "finite_emission_cells")
empty = [k for k in stages if int(st[k]) == 0]
if not empty:
    ok("every positional stage is live: " + ", ".join("%s=%s" % (k, st[k]) for k in stages))
else:
    no("positional stages are empty on c4_scale_118x119: %s" % ", ".join(empty))
# DIRECT VERIFICATION, not a window rescan: the whole read must be checked at FEWER starts than
# were proposed, and every verified FR state must come from an accepted placement.
if int(st["full_read_verifications"]) <= int(st["unique_seed_starts"]):
    ok("full-read verification runs at %s of %s seeded starts -- the window is never rebuilt"
       % (st["full_read_verifications"], st["unique_seed_starts"]))
else:
    no("more verifications (%s) than seeded starts (%s)"
       % (st["full_read_verifications"], st["unique_seed_starts"]))
if int(st["accepted_mate_placements"]) <= int(st["full_read_verifications"]):
    ok("accepted placements are a subset of the reads actually verified")
else:
    no("accepted (%s) exceeds verified (%s)"
       % (st["accepted_mate_placements"], st["full_read_verifications"]))
# nine_identical_origins keeps EVERY repeat origin: 9 cells but many more FR states.
ni = d["nine_identical_origins"]
if int(ni["verified_fr_states"]) > int(ni["finite_emission_cells"]):
    ok("identical repeat origins stay distinct states (%s states over %s cells)"
       % (ni["verified_fr_states"], ni["finite_emission_cells"]))
else:
    no("repeat origins collapsed: %s states over %s cells"
       % (ni["verified_fr_states"], ni["finite_emission_cells"]))
# EXACTNESS, everywhere.
diff = [k for k, v in d.items() if v["cells_differ"] != "0"]
if diff: no("finite-support cells differ from the dense oracle in: %s" % ", ".join(diff))
else: ok("identical finite-support cells in every case")
worst = max(float(v["worst_mass_diff"]) for v in d.values())
if worst == 0.0: ok("identical log mass in every case (exactly 0 difference)")
else: no("log mass differs by %.4g" % worst)
mism = [k for k, v in d.items() if v["informative_match"] != "1"]
if mism: no("informative classification differs in: %s" % ", ".join(mism))
else: ok("identical informative/uninformative classification in every case")
# The off-panel combination must actually be FOUND, not merely agreed upon as absent.
op = d["offpanel_A2_B1"]
if int(op["finite_supported"]) > 0:
    ok("the off-panel (A2,B1) combination is found (%s finite cell) -- no panel path carries it"
       % op["finite_supported"])
else:
    no("the off-panel combination has no finite cell; the fixture proves nothing")
# A free dimension is EXPANDED, not discarded.
for k, lab in (("A_only_constraint", "beta"), ("B_only_constraint", "alpha")):
    v = d[k]
    if int(v["finite_supported"]) > 1 and int(v["verified"]) > 1:
        ok("%s: the free %s dimension is expanded (%s verified, %s finite)"
           % (k, lab, v["verified"], v["finite_supported"]))
    else:
        no("%s: the free %s dimension was discarded" % (k, lab))
# An invariant-only seed constrains neither and must force the full product.
iv = d["invariant_only_seed"]
if int(iv["verified"]) == int(iv["dense_pairs"]):
    ok("an invariant-only seed constrains neither allele and expands to all %s pairs"
       % iv["dense_pairs"])
else:
    no("an invariant-only seed verified only %s of %s pairs" % (iv["verified"], iv["dense_pairs"]))
# Non-ACGT must take the EXHAUSTIVE fallback and still match.
na = d["non_acgt_fallback"]
if na["fallback"] == "1" and na["cells_differ"] == "0":
    ok("a non-ACGT read falls back to the exhaustive allele product and still matches exactly")
else:
    no("non-ACGT: fallback=%s cells_differ=%s" % (na["fallback"], na["cells_differ"]))
# Zero-state cells are still emitted.
z = d["zero_state_cells"]
if z["finite_dense"] == "0" and z["finite_supported"] == "0" and z["cells_differ"] == "0":
    ok("a fragment with no placement still emits every cell, as -inf")
else:
    no("zero-state handling differs: %s" % z)
# Multiplicity: nine identical origins must be preserved, not collapsed.
ni = d["nine_identical_origins"]
if int(ni["seed_occurrences"]) > 100 and ni["worst_mass_diff"] == "0":
    ok("nine identical repeat origins are preserved with identical mass (%s seed occurrences)"
       % ni["seed_occurrences"])
else:
    no("repeat multiplicity: seed_occurrences=%s mass_diff=%s"
       % (ni["seed_occurrences"], ni["worst_mass_diff"]))
# THE REDUCTION, at C4 scale.
c4 = d["c4_scale_118x119"]
if float(c4["reduction"]) > 5.0 and c4["cells_differ"] == "0":
    ok("at C4 scale %s dense pairs collapse to %s verified (%.1fx) with no loss"
       % (c4["dense_pairs"], c4["verified"], float(c4["reduction"])))
else:
    no("C4-scale reduction is only %sx" % c4["reduction"])
sys.exit(bad)
PYEOFS
  fails=$(( fails + $? ))
else
  bad "the support self-test produced no output"
fi

# ---------------------------------------------------------------------------------------------
# GATE 24: CHAIN ORIENTATION. A candidate's walk may run ANTIPARALLEL to the chain; its bytes are
# then the reverse complement of the reference-oriented sequence and its block spans run backwards.
# Anything comparing sequence or ordering intervals across candidates must work in CHAIN
# coordinates. Testing raw walk offsets refused every C4 edge with "blocks out of order", on frames
# that were perfectly correct.
#
# THE DECISIVE CASE is an EXACT reverse-complement duplicate: same steps reversed, every sign
# flipped. In chain orientation it must be indistinguishable from its original -- same span, same
# bytes -- while still being FLAGGED reverse_frame, and a mixed forward/reverse panel must build
# geometry and call normally.
"$PY" - "$OUT" <<'PYEOFJ'
import sys
d = sys.argv[1]
lines = open(d + "/g.gfa").read().rstrip("\n").split("\n")
out = []
for l in lines:
    out.append(l)
    f = l.split("\t")
    if f[0] == "P" and f[1] == "hapAA":
        steps = f[2].split(",")
        rc = ",".join(s[:-1] + ("-" if s[-1] == "+" else "+") for s in reversed(steps))
        out.append("P\thapAA_rc\t%s\t*" % rc)
open(d + "/g_rc.gfa", "w").write("\n".join(out) + "\n")
PYEOFJ
"$BIN" bubble -i "$OUT/g_rc.gfa" -r hapAA -o "$OUT/brc" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" genotype -i "$OUT/brc.sorted.gfa" -b "$OUT/brc" -r hapAA -o "$OUT/orc" \
  -R "$OUT/r1.fq" -R "$OUT/r2.fq" --fragment-len 350 \
  --hybrid-orientation-probe "$OUT/orient.tsv" -q >/dev/null 2>&1
"$BIN" genotype -i "$OUT/brc.sorted.gfa" -b "$OUT/brc" -r hapAA -o "$OUT/hrc" \
  -R "$OUT/r1.fq" -R "$OUT/r2.fq" --fragment-len 350 --hybrid-call -q >/dev/null 2>&1
if [ -s "$OUT/orient.tsv" ] && [ -s "$OUT/hrc.genotypes.tsv" ]; then
  "$PY" - "$OUT/orient.tsv" "$OUT/hrc.genotypes.tsv" <<'PYEOFK'
import sys
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])][1:]
d = {r[0]: r for r in rows}
bad = 0
def ok(m): print("  ok   " + m)
def no(m):
    global bad; bad += 1; print("  FAIL " + m)
if "hapAA" not in d or "hapAA_rc" not in d:
    no("the reverse-complement duplicate is not in the panel"); sys.exit(1)
a, b = d["hapAA"], d["hapAA_rc"]
if a[1] == "0" and b[1] == "1":
    ok("the duplicate is FLAGGED antiparallel (reverse_frame 0 vs 1), not silently normalised")
else:
    no("reverse_frame flags are %s and %s" % (a[1], b[1]))
if a[2] == b[2] and a[3] == b[3]:
    ok("chain-oriented SEQUENCE is identical (md5 %s, %s bp)" % (a[2][:12], a[3]))
else:
    no("chain-oriented sequence differs: %s/%s vs %s/%s" % (a[2][:12], a[3], b[2][:12], b[3]))
if a[4] == b[4]:
    ok("chain-oriented BLOCK SPANS are identical (%s)" % a[4])
else:
    no("chain-oriented block spans differ: %s vs %s" % (a[4], b[4]))
g = [l.rstrip("\n").split("\t") for l in open(sys.argv[2])][1:]
pairs = {(r[9], r[10]) for r in g}
if pairs and pairs <= {("hapAB", "hapBA"), ("hapBA", "hapAB")}:
    ok("a MIXED forward/reverse panel still recovers the correct phase (%s)" % sorted(pairs))
else:
    no("mixed panel phase is %s, expected hapAB/hapBA" % sorted(pairs))
sys.exit(bad)
PYEOFK
  fails=$(( fails + $? ))
  "$BIN" genotype -i "$OUT/brc.sorted.gfa" -b "$OUT/brc" -r hapAA -o "$OUT/grc" \
    -R "$OUT/r1.fq" -R "$OUT/r2.fq" --fragment-len 350 \
    --hybrid-geometry-probe "$OUT/geom_rc.tsv" -q >/dev/null 2>&1
  # COLUMN LOOKED UP BY NAME. This read column 8 by position and silently started reading
  # lflank_derived when the derived-flank columns were added -- the same staleness that made the
  # psi mean/sum check report a phantom failure.
  GOK=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="ok") c=i; next}
                     $1==1 && $2==2 && $7==0 {print $c}' "$OUT/geom_rc.tsv" 2>/dev/null)
  [ "${GOK:-0}" = "1" ] \
    && ok "linkage geometry builds across the variable pair on a mixed-orientation panel" \
    || bad "geometry still refuses the variable pair on a mixed panel (ok=${GOK:-none})"
else
  bad "the orientation probe or the mixed-panel call produced nothing"
fi
# GATE 22: CANDIDATE-FRAME COVERAGE OVER THE HMM STATE UNIVERSE. The requirement is over the
# states the marker HMM actually declares, NOT the raw panel paths: a recorded state reduction is
# legitimate, silently dropping candidates whose frame construction failed is not -- that would make
# the linkage topology and the posterior depend on an undocumented change of state space.
#
# A panel that fails is an INCOMPLETE hybrid-SUBSTRATE result, not evidence that the biological
# genotype is ambiguous, and the report must say so:
#     hybrid_status INCOMPLETE / reason candidate-frame-coverage / legacy_call_status AVAILABLE
"$BIN" genotype -i "$OUT/b.sorted.gfa" -b "$OUT/b" -r hapAA -o "$OUT/pf" \
  -R "$OUT/r1.fq" -R "$OUT/r2.fq" --hybrid-preflight "$OUT/preflight.tsv" -q >/dev/null 2>&1
if [ -s "$OUT/preflight.tsv" ]; then
  "$PY" - "$OUT/preflight.tsv" <<'PYEOFH'
import sys
kv = {}; states = []; missing = []
for l in open(sys.argv[1]):
    f = l.rstrip("\n").split("\t")
    if f[0] == "hmm_state": states.append(f[1])
    elif f[0] == "missing": missing.append((f[1], f[2] if len(f) > 2 else "?"))
    elif len(f) >= 2: kv[f[0]] = f[1]
bad = 0
def ok(m): print("  ok   " + m)
def no(m):
    global bad; bad += 1; print("  FAIL " + m)
for k in ("raw_panel_paths","hmm_states","framed_states","complete_frames",
          "accepted_partial_frames","missing_states","names_unique","all_states_usable"):
    if k not in kv: no("preflight missing metric %s" % k); sys.exit(1)
# THE UNIVERSE IS NAMED, not counted, so it is reproducible and a later disagreement is attributable.
if len(states) == int(kv["hmm_states"]):
    ok("the HMM state universe is enumerated by name (%d states), not just counted" % len(states))
else:
    no("%d state names listed but hmm_states says %s" % (len(states), kv["hmm_states"]))
if kv["framed_states"] == kv["hmm_states"] and kv["missing_states"] == "0":
    ok("every HMM state has a verified candidate frame (%s/%s)"
       % (kv["framed_states"], kv["hmm_states"]))
else:
    no("only %s of %s HMM states framed; missing: %s"
       % (kv["framed_states"], kv["hmm_states"], missing[:3]))
# UNIQUENESS BEFORE COMPARISON. Sorted-vector equality alone would let a duplicated state name
# appear on both sides and cancel out, so neither side may contain a repeat.
if kv["names_unique"] == "1":
    ok("HMM state names and framed names are each unique -- duplicates cannot cancel out")
else:
    no("a duplicate name exists; a sorted comparison would pass on both sides regardless")
if len(set(states)) == len(states):
    ok("the enumerated state universe contains no repeated name (%d distinct)" % len(set(states)))
else:
    no("the enumerated state universe repeats a name")
# Counts agreeing is not enough: the NAMES must be the same set, or a state could be replaced.
if kv["all_states_usable"] == "1":
    ok("all_states_usable: unique names, nothing missing, framed set equals the state set")
else:
    no("coverage incomplete -- this is a hybrid-SUBSTRATE result, not biological ambiguity")
# A VERIFIED PARTIAL TERMINAL FRAME IS ACCEPTED, not missing. Reported separately so the case stays
# visible; on CYP2D6 exactly one path (NA18989#1#haplotype1) ends inside a block.
if int(kv["complete_frames"]) + int(kv["accepted_partial_frames"]) == int(kv["framed_states"]):
    ok("framed = complete (%s) + accepted partial (%s), reported separately"
       % (kv["complete_frames"], kv["accepted_partial_frames"]))
else:
    no("complete %s + partial %s != framed %s"
       % (kv["complete_frames"], kv["accepted_partial_frames"], kv["framed_states"]))
# "frame-partial" must never appear as a REFUSAL reason: it was unreachable code, and it implied
# partial frames are rejected when the contract accepts them.
if any(r == "frame-partial" for _, r in missing):
    no("a state was refused with reason frame-partial; verified partial frames are ACCEPTED")
else:
    ok("no state is refused as frame-partial (that reason was unreachable and is gone)")
if len(missing) != int(kv["missing_states"]):
    no("%d missing rows but missing_states says %s" % (len(missing), kv["missing_states"]))
else:
    ok("every missing state would be listed with a reason (%d here)" % len(missing))
sys.exit(bad)
PYEOFH
  fails=$(( fails + $? ))
else
  bad "hybrid preflight produced no output"
fi

# GATE 20: TRANSACTIONAL ACTIVATION. Subtracting linkage-owned fragments from the marker unaries
# and activating their edges is ONE transaction. Half of it makes those fragments vanish from BOTH
# models -- removed from the marker counts, consumed by no edge because the edge was refused --
# leaving a run quietly WEAKER than the legacy caller it extends. So nothing is subtracted unless
# every required edge was built and the model is complete.
#
# THE EQUALITY is phrased against ACTIVE edges, not merely linkage-owned fragments, so a refused
# edge cannot satisfy the exclusion side by accident:
#     {excluded fragments} == {fragments consumed by ACTIVE linkage edges},  each exactly once.
"$BIN" genotype-frag -i /dev/null -b none -o "$OUT/ac" --activation-selftest \
  > "$OUT/act.tsv" 2>/dev/null
if [ -s "$OUT/act.tsv" ]; then
  "$PY" - "$OUT/act.tsv" <<'PYEOFG'
import sys
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])]
d = {r[0]: r for r in rows[1:]}
bad = 0
def ok(m): print("  ok   " + m)
def no(m):
    global bad; bad += 1; print("  FAIL " + m)
need = ("no_linkage","complete","refused_edge","wide_present","mapping_refused","unusable_present")
miss = [k for k in need if k not in d]
if miss: no("activation self-test missing: %s" % ", ".join(miss)); sys.exit(1)
# THE EQUALITY must hold in every case, activated or not.
viol = [k for k in need if d[k][9] != "1"]
if viol: no("excluded != consumed-by-active-edges in: %s" % ", ".join(viol))
else: ok("excluded IDs equal active-edge-consumed IDs in all %d cases" % len(need))
# Each excluded fragment consumed exactly once: the excluded list has no duplicates.
dup = [k for k in need if d[k][6] != d[k][8]]
if dup: no("duplicate exclusions in: %s" % ", ".join(dup))
else: ok("each excluded fragment is consumed exactly once (no duplicates)")
# No linkage evidence -> activated, nothing excluded: the inference is legacy by construction.
n0 = d["no_linkage"]
if n0[3] == "1" and n0[6] == "0" and n0[5] == "0":
    ok("no linkage evidence: activates with 0 active edges and 0 exclusions -- legacy inference")
else: no("no_linkage gave activated=%s edges=%s excluded=%s" % (n0[3], n0[5], n0[6]))
# A complete hybrid excludes exactly the fragments its active edges consume.
c = d["complete"]
if c[3] == "1" and c[5] == "1" and c[6] == "2" and c[7] == "2":
    ok("complete hybrid: 1 active edge consuming 2 fragments, both excluded")
else: no("complete gave activated=%s edges=%s excluded=%s consumed=%s" % (c[3],c[5],c[6],c[7]))
# THE TRANSACTIONAL PROPERTY: every failure path excludes NOTHING. Otherwise the fragments would be
# gone from the markers and consumed by nobody.
for case in ("refused_edge","wide_present","mapping_refused","unusable_present"):
    r = d[case]
    if r[3] == "0" and r[6] == "0" and r[5] == "0" and r[4] == "INCOMPLETE":
        ok("%s: no activation, NOTHING subtracted, status INCOMPLETE, reason recorded" % case)
    else:
        no("%s: activated=%s excluded=%s active_edges=%s status=%s -- a partial transaction committed"
           % (case, r[3], r[6], r[5], r[4]))
# Completeness passing is NOT sufficient: an edge that cannot be built must still abandon the whole
# transaction rather than run with a hole.
# THE FOUR CONDITIONS MUST STAY DISTINCT. mapping_refused is OWNERSHIP-complete and
# FACTOR-incomplete, and its final call_status must be INCOMPLETE -- reporting a single "complete"
# here is how an incomplete call comes to look finished.
m = d["mapping_refused"]
if m[1] == "1" and m[2] == "0" and m[3] == "0" and m[4] == "INCOMPLETE" and m[6] == "0":
    ok("mapping_refused: ownership-complete but factor-INcomplete -> call_status INCOMPLETE, "
       "nothing subtracted")
else:
    no("mapping_refused: ownership=%s factors=%s activated=%s status=%s excluded=%s"
       % (m[1], m[2], m[3], m[4], m[6]))
# call_status is COMPLETE only when all three conditions hold, never on a subset.
for case in need:
    r = d[case]
    allthree = (r[1] == "1" and r[2] == "1" and r[3] == "1")
    if (r[4] == "COMPLETE") != allthree:
        no("%s: call_status %s but conditions were %s/%s/%s" % (case, r[4], r[1], r[2], r[3]))
        break
else:
    ok("call_status is COMPLETE exactly when ownership, factors and activation all hold")
for case in need:
    if d[case][3] == "0" and d[case][10] == "-":
        no("%s: refused activation without recording a reason" % case)
sys.exit(bad)
PYEOFG
  fails=$(( fails + $? ))
else
  bad "activation self-test produced no output"
fi

# GATE 19: HYBRID COMPLETENESS. Three situations must stay distinct: no linkage-owned evidence (a
# valid legacy edge), linkage successfully represented (a valid linked edge), and linkage or Wide
# evidence that EXISTS but cannot be represented (INCOMPLETE/UNSUPPORTED). An inactive
# ChainEdgeLinkage cannot express the third -- the kernel reads it exactly like the first -- so a
# refused edge would silently become neutral and the posterior would be presented as complete while
# a reduced evidence model ran.
#
#     hybrid COMPLETE  <=>  every non-invariant owned fragment has a SUPPORTED CONSUMER
#
# NOT "unresolved": the model did not evaluate the evidence, it could not consume it.
"$BIN" genotype-frag -i /dev/null -b none -o "$OUT/cp" --completeness-selftest \
  > "$OUT/comp.tsv" 2>/dev/null
if [ -s "$OUT/comp.tsv" ]; then
  "$PY" - "$OUT/comp.tsv" <<'PYEOFF'
import sys
lines = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])]
comp = {}; kern = {}
mode = None
for f in lines:
    if f[0] == "case": mode = "comp"; continue
    if f[0] == "kernel_edge_case": mode = "kern"; continue
    (comp if mode == "comp" else kern)[f[0]] = f
bad = 0
def ok(m): print("  ok   " + m)
def no(m):
    global bad; bad += 1; print("  FAIL " + m)
def want(case, complete):
    r = comp.get(case)
    if r is None: no("completeness case %s missing" % case); return None
    got = r[1] == "1"
    if got == complete:
        ok("%s -> %s" % (case, "COMPLETE" if complete else "INCOMPLETE"))
    else:
        no("%s reported %s, expected %s"
           % (case, "COMPLETE" if got else "INCOMPLETE", "COMPLETE" if complete else "INCOMPLETE"))
    return r
want("no_linkage", True); want("linkage_ok", True); want("invariant_only", True)
for c in ("too_many_configs","exposure_refused","invalid_emission","mapping_refused",
          "wide_fragment","unusable_fragment"):
    want(c, False)
# EVERY refusal survives, not just the last -- which is when the report matters most.
m = comp.get("multi_refusal")
if m and m[9] == "3" and m[10].count(";") == 2 and len({x.split(":")[1] for x in m[10].split(";")}) == 3:
    ok("multi_refusal retains all 3 refusals with 3 DISTINCT reasons (%s)" % m[10])
else:
    no("multi_refusal lost refusals or reasons: n=%s reasons=%s" % (m[9] if m else "?", m[10] if m else "?"))
# A Wide fragment is reported with the variable scope it actually depended on.
w = comp.get("wide_fragment")
if w and w[6] == "1" and w[11] != "-":
    ok("a Wide fragment is unsupported and reports its variable scope (%s)" % w[11])
else:
    no("Wide fragment scope not reported: wide=%s scopes=%s" % (w[6] if w else "?", w[11] if w else "?"))
# An invariant fragment needs no consumer and must not make the run incomplete.
iv = comp.get("invariant_only")
if iv and iv[3] == "2" and iv[1] == "1":
    ok("invariant fragments need no consumer and keep the run complete")
else:
    no("invariant handling wrong: %s" % (iv if iv else "missing"))
# NO REFUSED EDGE TABLE REACHES THE KERNEL.
if kern.get("usable_edge", ["","0"])[1] == "1":
    ok("a usable edge builds an ACTIVE kernel entry (psi %s)" % kern["usable_edge"][4])
else:
    no("a usable edge did not build an active kernel entry")
leaked = [k for k, v in kern.items() if k != "usable_edge" and (v[1] != "0" or v[4] != "0")]
if leaked:
    no("refused edges reached the kernel with a table: %s" % ", ".join(leaked))
else:
    ok("every refused edge and refused mapping yields an INACTIVE entry with no table (%d cases)"
       % (len(kern) - 1))
sys.exit(bad)
PYEOFF
  fails=$(( fails + $? ))
else
  bad "completeness self-test produced no output"
fi

# GATE 18: HAPLOTYPE -> ALLELE MAPPING. The int -> unsigned boundary is where this breaks silently:
# BlockAlleles reports -1 for "no allele here", and a bare cast makes that 4294967295, which then
# indexes log_psi far out of bounds without any downstream check noticing. Every value is validated
# BEFORE conversion, and a refused mapping exposes an EMPTY vector so a consumer that checked the
# wrong thing has nothing to index.
"$BIN" genotype-frag -i /dev/null -b none -o "$OUT/mp" --mapping-selftest > "$OUT/map.tsv" 2>/dev/null
if [ -s "$OUT/map.tsv" ]; then
  "$PY" - "$OUT/map.tsv" <<'PYEOFE'
import sys
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])]
hdr, rows = rows[0], rows[1:]
d = {r[0]: r for r in rows}
bad = 0
def ok(m): print("  ok   " + m)
def no(m):
    global bad; bad += 1; print("  FAIL " + m)
need = ("normal","permuted","bypass","missing","missing_large","out_of_range","bad_bypass")
miss = [k for k in need if k not in d]
if miss: no("mapping self-test missing cases: %s" % ", ".join(miss)); sys.exit(1)
# A HAPLOTYPE INDEX IS NOT AN ALLELE INDEX. The permuted case maps 3 haplotypes onto 2 alleles with
# no fixed point at haplotype 0, so anything treating one as the other is wrong here.
if d["normal"][1] == "ok" and d["normal"][6] == "0,1,0":
    ok("normal mapping resolves (%s)" % d["normal"][6])
else: no("normal mapping failed: %s %s" % (d["normal"][1], d["normal"][6]))
if d["permuted"][1] == "ok" and d["permuted"][6] == "1,0,1":
    ok("permuted mapping resolves (%s) -- haplotype index is never used as an allele index"
       % d["permuted"][6])
else: no("permuted mapping failed: %s %s" % (d["permuted"][1], d["permuted"][6]))
# A BYPASSING haplotype resolves to bypass_allele -- a real state -- not to -1 and not to allele 0.
b = d["bypass"]
if b[1] == "ok" and b[3] == "1" and b[6] == "0,2,1":
    ok("a bypassing haplotype resolves to bypass_allele (%s), not -1 and not allele 0" % b[6])
else: no("bypass mapping gave status=%s resolved=%s mapping=%s" % (b[1], b[3], b[6]))
# REFUSALS: never cast, never wrapped, never partially exposed.
for case, want in (("missing","missing-mapping"), ("missing_large","missing-mapping"),
                   ("out_of_range","allele-out-of-range"), ("bad_bypass","allele-out-of-range")):
    r = d[case]
    if r[1] != want:
        no("%s: status %s, expected %s" % (case, r[1], want)); continue
    if r[7] != "0":
        no("%s: refused but still exposes a %s-element vector to index" % (case, r[7])); continue
    ok("%s: refused as %s at haplotype %s (value %s), nothing indexable" % (case, r[1], r[4], r[5]))
# The large negative must be REFUSED, not wrapped into a huge unsigned index.
if d["missing_large"][5] == "-999" and d["missing_large"][1] == "missing-mapping":
    ok("a large negative (-999) is refused at the boundary, not cast to 4294966297")
else:
    no("a large negative was not caught at the int -> unsigned boundary")
sys.exit(bad)
PYEOFE
  fails=$(( fails + $? ))
else
  bad "mapping self-test produced no output"
fi

# GATE 17: MARKER OCCURRENCE EXCLUSION. A fragment owned by a linkage edge contributes its sequence
# there and must leave the marker counts, or the same read is counted twice. What must be subtracted
# is its OCCURRENCES, not the markers themselves: a marker it shares with unary-owned fragments has
# to keep their counts and stay in the panel. Deleting the marker instead would remove evidence that
# was never double counted, and at a marker-poor block that is the difference between a call and a
# no-call.
#
# The claim is exact and is asserted as such: counts(excluded) == counts(all) - counts(that fragment
# alone), slot by slot. Overlapping fragments guarantee shared markers, which is the case that
# matters; a fixture whose fragments shared nothing would assert nothing.
"$PY" - "$OUT" <<'PYEOFC'
import sys, random
out = sys.argv[1]; B = "ACGT"
def seq(n, s):
    r = random.Random(s); return "".join(r.choice(B) for _ in range(n))
def mut(s, k, sd):
    r = random.Random(sd); s = list(s)
    for p in r.sample(range(len(s)), k): s[p] = r.choice([c for c in B if c != s[p]])
    return "".join(s)
X = seq(3000,1); A1 = seq(500,2); B1 = mut(A1,70,3); J = seq(40,4)
A2 = seq(500,5); B2 = mut(A2,70,6); Z = seq(3000,7)
segs = [("1",X),("2",A1),("3",B1),("4",J),("5",A2),("6",B2),("7",Z)]
links = [("1","2"),("1","3"),("2","4"),("3","4"),("4","5"),("4","6"),("5","7"),("6","7")]
paths = {"hapAA":["1","2","4","5","7"],"hapAB":["1","2","4","6","7"],
         "hapBA":["1","3","4","5","7"],"hapBB":["1","3","4","6","7"]}
d = dict(segs)
with open(out+"/mx.gfa","w") as g:
    g.write("H\tVN:Z:1.0\n")
    for n,x in segs: g.write("S\t%s\t%s\n"%(n,x))
    for a,b in links: g.write("L\t%s\t+\t%s\t+\t0M\n"%(a,b))
    for n,st in paths.items(): g.write("P\t%s\t%s\t*\n"%(n,",".join(x+"+" for x in st)))
h = "".join(d[x] for x in paths["hapAB"])
comp = {"A":"T","C":"G","G":"C","T":"A"}
k = 0
with open(out+"/mx.r1.fq","w") as f1, open(out+"/mx.r2.fq","w") as f2:
    for i in range(0, len(h)-360, 4):
        a = h[i:i+150]; b = h[i+200:i+350]
        if len(b) < 150: break
        rc = "".join(comp[c] for c in reversed(b))
        f1.write("@f%d/1\n%s\n+\n%s\n"%(k,a,"I"*150))
        f2.write("@f%d/2\n%s\n+\n%s\n"%(k,rc,"I"*150))
        k += 1
PYEOFC
"$BIN" bubble -i "$OUT/mx.gfa" -r hapAA -o "$OUT/mxb" --min-variant-bp 0 -q >/dev/null 2>&1
echo "f400" > "$OUT/mx.excl.txt"
awk '/^@f400\/1$/{p=1} p{print; if(++n==4) exit}' "$OUT/mx.r1.fq" > "$OUT/mx.one.r1.fq"
awk '/^@f400\/2$/{p=1} p{print; if(++n==4) exit}' "$OUT/mx.r2.fq" > "$OUT/mx.one.r2.fq"
"$BIN" genotype -i "$OUT/mxb.sorted.gfa" -b "$OUT/mxb" -r hapAA -o "$OUT/mx.all" \
  -R "$OUT/mx.r1.fq" -R "$OUT/mx.r2.fq" --dump-markers "$OUT/mx.all.tsv" -q >/dev/null 2>&1
"$BIN" genotype -i "$OUT/mxb.sorted.gfa" -b "$OUT/mxb" -r hapAA -o "$OUT/mx.one" \
  -R "$OUT/mx.one.r1.fq" -R "$OUT/mx.one.r2.fq" --dump-markers "$OUT/mx.one.tsv" -q >/dev/null 2>&1
"$BIN" genotype -i "$OUT/mxb.sorted.gfa" -b "$OUT/mxb" -r hapAA -o "$OUT/mx.ex" \
  -R "$OUT/mx.r1.fq" -R "$OUT/mx.r2.fq" --exclude-fragments "$OUT/mx.excl.txt" \
  --dump-markers "$OUT/mx.ex.tsv" -q >/dev/null 2>&1
if [ -s "$OUT/mx.all.tsv" ] && [ -s "$OUT/mx.one.tsv" ] && [ -s "$OUT/mx.ex.tsv" ]; then
  "$PY" - "$OUT/mx.all.tsv" "$OUT/mx.one.tsv" "$OUT/mx.ex.tsv" <<'PYEOFD'
import sys
def load(p):
    rows = {}
    for i, l in enumerate(open(p)):
        if i == 0: continue
        f = l.rstrip("\n").split("\t")
        rows[(f[0], f[3], f[4])] = int(f[5])
    return rows
A, O, E = (load(x) for x in sys.argv[1:4])
bad = 0
def ok(m): print("  ok   " + m)
def no(m):
    global bad; bad += 1; print("  FAIL " + m)
if set(A) == set(O) == set(E):
    ok("no marker was deleted: all %d slots present in every arm" % len(A))
else:
    no("the marker panel differs between arms -- exclusion removed markers, not occurrences")
mism = [k for k in A if E[k] != A[k] - O[k]]
if not mism:
    ok("exclusion subtracts exactly the fragment's occurrences, slot by slot (%d slots)" % len(A))
else:
    no("%d slots where excluded != all - fragment_alone, e.g. %s" % (len(mism), mism[0]))
touched = [k for k in A if O[k] > 0]
shared = [k for k in touched if A[k] > O[k]]
zeroed = [k for k in shared if E[k] == 0]
if not shared:
    no("the excluded fragment shares no marker with any other -- the gate asserts nothing")
elif zeroed:
    no("%d SHARED markers were zeroed by exclusion; they must keep the other fragments' counts"
       % len(zeroed))
else:
    kk = shared[0]
    ok("%d of %d touched markers are shared, none zeroed (e.g. %d -> %d, minus %d)"
       % (len(shared), len(touched), A[kk], E[kk], O[kk]))
ta, to, te = sum(A.values()), sum(O.values()), sum(E.values())
if te == ta - to: ok("total occurrences: %d - %d = %d" % (ta, to, te))
else: no("total occurrences %d != %d - %d" % (te, ta, to))
sys.exit(bad)
PYEOFD
  fails=$(( fails + $? ))
else
  bad "marker exclusion fixture produced no dumps"
fi

# GATE 14: THE BRUTE-FORCE PATH ORACLE. Forward-backward is compared against an INDEPENDENT
# enumeration of every ordered diploid state path, scoring
#     SUM_b log E_b(s_b) + SUM_b log T_LS(s_{b-1},s_b) + SUM_b log psi_b(s_{b-1},s_b)
# from scratch. Both the total log partition and every block marginal must agree. A best-path
# comparison would catch none of: linkage at the wrong edge, either factor applied twice, an
# accidental transition row-normalisation, an ordered-state mapping error, or a correct best call
# resting on wrong posterior mass.
#
# The oracle does NOT share the recursion's potential helper. Sharing it would make a bug inside
# that helper move both arms identically and cancel out of the comparison -- the mapping error above
# all, which is the single place catalogue indices, path identities and ordered homologues can be
# confused. Measured with the shared helper deliberately restored, all four mutations below are
# caught: linkage at the wrong edge 0.056, linkage twice 0.106, mapping swapped 0.056, Li-Stephens
# twice 0.446 nats of log-partition disagreement.
"$BIN" genotype-frag -i /dev/null -b none -o "$OUT/ho" --hybrid-oracle > "$OUT/oracle.tsv" 2>/dev/null
if [ -s "$OUT/oracle.tsv" ]; then
  "$PY" - "$OUT/oracle.tsv" <<'PYEOFB'
import sys
# ONLY the key/value lines: the probe also emits wide grouped-sparse rows, and a bare
# dict(split) over every line fails on them.
d = {}
for _l in open(sys.argv[1]):
    _f = _l.rstrip("\n").split("\t")
    if len(_f) == 2: d[_f[0]] = _f[1]
bad = 0
def ok(m): print("  ok   " + m)
def no(m):
    global bad; bad += 1; print("  FAIL " + m)
TOL = 1e-9
for key, label in (("log_weight_sum_abs_diff", "log weight sum"),
                   ("worst_marginal_abs_diff", "every block marginal")):
    v = float(d[key])
    if v < TOL: ok("%s agrees with the brute-force path oracle (%.2e)" % (label, v))
    else: no("%s disagrees with the oracle by %.6g" % (label, v))
# A LINKAGE-FREE control separates "the recursion is right" from "the linkage table is right".
for key, label in (("log_weight_sum_abs_diff_no_linkage", "log weight sum"),
                   ("worst_marginal_abs_diff_no_linkage", "block marginals")):
    v = float(d[key])
    if v < TOL: ok("linkage-free control: %s agrees (%.2e)" % (label, v))
    else: no("linkage-free control: %s disagrees by %.6g" % (label, v))
sd = float(d["worst_marginal_sum_dev"])
if sd < 1e-9: ok("block marginals are distributions (worst sum deviation %.2e)" % sd)
else: no("a block marginal sums to 1 +/- %.4g" % sd)
# NON-VACUITY: linkage must actually move the posterior, or the agreement above proves only that
# two implementations of plain Li-Stephens match.
le = float(d["linkage_marginal_effect"])
if le > 1e-6: ok("linkage moves the posterior by %.4f -- the comparison is not vacuous" % le)
else: no("linkage changes the posterior by only %.2e; the oracle asserts nothing about psi" % le)
# THE GROUPED SPARSE CONTRACTION must equal the dense kernel exactly:
#     F'(y) = F'_LS(y) + SUM_x F(x) T(x,y) [psi(x,y) - 1]
# across the regimes a single recombination rate or a single class would conceal: r = 0 (identity
# term only), r = 1 (uniform switch only), an intermediate r with all four expanded terms live,
# several classes on one edge, several haplotypes per allele, TWO linked edges in one chain, and a
# near-total cancellation.
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])]
hi = [k for k, r in enumerate(rows) if r and r[0] == "sparse_case"]
sp = {}
if hi:
    h = rows[hi[0]]
    for r in rows[hi[0]+1:]:
        if len(r) == len(h): sp[r[0]] = dict(zip(h, r))
if not sp:
    no("no grouped-sparse comparison was reported")
else:
    need = ("r0_identity_only", "r1_uniform_switch", "r_intermediate", "two_linked_edges",
            "no_linkage_control", "multi_class_3x3", "multi_class_r0", "multi_class_r1",
            "adversarial_mass_in_zero_phase", "swapsym_no_linkage", "swapsym_linked",
            "swapsym_multi_class")
    miss = [k for k in need if k not in sp]
    if miss: no("grouped-kernel cases missing: %s" % ", ".join(miss))
    else: ok("grouped kernel covered over %d cases: r=0, r=1, intermediate, multi-class, two edges"
             % len(need))
    worst_w = max(float(v["logw_absdiff"]) for v in sp.values())
    worst_m = max(float(v["marg_absdiff"]) for v in sp.values())
    if worst_w < 1e-9 and worst_m < 1e-9:
        ok("grouped == dense everywhere: weight sums %.1e, marginals %.1e" % (worst_w, worst_m))
    else:
        bad_cases = [k for k, v in sp.items()
                     if float(v["logw_absdiff"]) >= 1e-9 or float(v["marg_absdiff"]) >= 1e-9]
        no("grouped differs from dense in: %s" % ", ".join(bad_cases))
    # Each regime must actually exercise its terms: r=0 and r=1 must give DIFFERENT posteriors, or
    # the transition terms are not being distinguished at all.
    if abs(float(sp["r0_identity_only"]["block1_max_marginal"]) -
           float(sp["r1_uniform_switch"]["block1_max_marginal"])) > 1e-6:
        ok("r=0 and r=1 give distinct posteriors (%.4f vs %.4f) -- both terms are live"
           % (float(sp["r0_identity_only"]["block1_max_marginal"]),
              float(sp["r1_uniform_switch"]["block1_max_marginal"])))
    else:
        no("r=0 and r=1 give the same posterior; the transition terms are not distinguished")
    # An edge with no classes takes the factorised path unchanged.
    c = sp["no_linkage_control"]
    if c["classes"] == "0" and c["corrections"] == "0":
        ok("an edge with zero classes applies NO corrections -- factorised path unchanged")
    else:
        no("the no-linkage control applied %s corrections over %s classes"
           % (c["corrections"], c["classes"]))
    # Several classes, and several corrections per class.
    mc = sp["multi_class_3x3"]
    if int(mc["classes"]) > 1 and int(mc["corrections"]) > int(sp["r_intermediate"]["corrections"]):
        ok("a 3x3 edge carries %s classes and %s corrections, against %s for the 2x2 single class"
           % (mc["classes"], mc["corrections"], sp["r_intermediate"]["corrections"]))
    else:
        no("the multi-class case is not richer than the single-class one")
    # TWO LINKED EDGES in one chain, not one linked plus one factorised.
    if int(sp["two_linked_edges"]["linked_edges"]) == 2 and \
       int(sp["two_linked_edges"]["corrections"]) == 2 * int(sp["r_intermediate"]["corrections"]):
        ok("two linked edges in one chain apply exactly twice the corrections")
    else:
        no("two linked edges: %s edges, %s corrections"
           % (sp["two_linked_edges"]["linked_edges"], sp["two_linked_edges"]["corrections"]))
    # HOMOLOGUE-SWAP INVARIANCE is a property of the KERNEL and can only be tested against
    # emissions that are themselves symmetric. Measured against the deliberately asymmetric default
    # emissions it fails even with NO linkage (0.108), which says nothing about the kernel.
    asym = float(sp["no_linkage_control"]["swap_absdiff"])
    sym = max(float(sp[k]["swap_absdiff"]) for k in
              ("swapsym_no_linkage", "swapsym_linked", "swapsym_multi_class"))
    if sym < 1e-12 and asym > 1e-3:
        ok("with symmetric emissions the posterior is swap-invariant (%.1e); the asymmetric fixture "
           "differs by %.3f even unlinked, so that arm tests the fixture, not the kernel"
           % (sym, asym))
    elif sym >= 1e-12:
        no("the posterior is not swap-invariant under symmetric emissions: %.4g" % sym)
    else:
        no("the asymmetric control shows no swap difference; the comparison proves nothing")
    # NEAR-TOTAL CANCELLATION must remain a distinct regime and produce no material negative weight.
    adv = sp["adversarial_mass_in_zero_phase"]
    if abs(float(adv["block1_max_marginal"]) -
           float(sp["r_intermediate"]["block1_max_marginal"])) < 1e-6:
        no("the adversarial case is indistinguishable from the ordinary one")
    elif int(adv["clamped"]) == 0:
        ok("near-total cancellation stays a distinct regime (%.4f) with NO negative weights"
           % float(adv["block1_max_marginal"]))
    else:
        ok("near-total cancellation clamped %s tiny negatives, worst %s -- recorded, not absorbed"
           % (adv["clamped"], adv["worst_negative"]))

# BOTH KERNEL EDGE PATHS MUST RUN IN ONE CHAIN.# BOTH KERNEL EDGE PATHS MUST RUN IN ONE CHAIN. There is a single inference kernel with a
# factorised O(n_h^2) path and a linked O(n_h^4) path; a branch that is present but never taken is
# not covered. The linkage-free control must take the factorised path only.
fac, lnk = int(d["factorised_edges"]), int(d["linked_edges"])
if fac > 0 and lnk > 0:
    ok("one chain exercised BOTH kernel paths: %d factorised edge(s), %d linked" % (fac, lnk))
else:
    no("the fixture took only one kernel path: %d factorised, %d linked" % (fac, lnk))
fac0, lnk0 = int(d["factorised_edges_no_linkage"]), int(d["linked_edges_no_linkage"])
if lnk0 == 0 and fac0 > 0:
    ok("the linkage-free control takes the factorised path only (%d edges)" % fac0)
else:
    no("the linkage-free control took %d linked edge(s)" % lnk0)
sys.exit(bad)
PYEOFB
  fails=$(( fails + $? ))
else
  bad "hybrid oracle produced no output"
fi

# GATE 12: NEUTRALITY OF AN UNINFORMATIVE EDGE. This is the assertion that matters most, and it
# cannot be produced from a read fixture: it needs an edge with NO fragments, and one whose emissions
# are identical across every configuration. --linkage-selftest constructs both directly, over content
# classes of all three cardinalities (1, 2 and 4) -- exactly where a sum-one normalisation leaks a
# -log|C| penalty against heterozygous content.
"$BIN" genotype-frag -i /dev/null -b none -o "$OUT/st" --linkage-selftest > "$OUT/self.tsv" 2>/dev/null
if [ -s "$OUT/self.tsv" ]; then
  "$PY" - "$OUT/self.tsv" <<'PYEOFA'
import sys
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])][1:]
d = {r[0]: r for r in rows}
bad = 0
def ok(m): print("  ok   " + m)
def no(m):
    global bad; bad += 1; print("  FAIL " + m)
need = ("zero_fragments", "flat_emissions", "all_unplaced", "informative",
        "unequal_exposure", "unequal_exposure_with_frags", "invalid_emission",
        "too_many_configs", "lambda_crossover_lo", "lambda_crossover_hi")
miss = [k for k in need if k not in d]
if miss: no("self-test missing cases: %s" % ", ".join(miss)); sys.exit(1)
# THE DENSE TABLE MUST REFUSE, NOT TRUNCATE. At LPA scale (457 x 410) it is 35.1 billion
# configurations and 561.7 GB, so the dense form cannot claim to handle every locus. A silent
# truncation would let the marker shortlist become an uncertified linkage cutoff.
tm = d["too_many_configs"]
if tm[6] == "0" and tm[8] == "too-many-configurations":
    ok("too_many_configs: a 64x64 edge (16.8M configs) is REFUSED, not truncated")
else:
    no("too_many_configs: usable=%s status=%s -- an oversized edge was accepted" % (tm[6], tm[8]))
sizes = d["zero_fragments"][5]
if "1x" in sizes and "2x" in sizes and "4x" in sizes:
    ok("self-test covers content classes of all three cardinalities (%s)" % sizes)
else:
    no("self-test does not cover class sizes 1, 2 and 4: %s" % sizes)
# THE NEUTRALITY PROPERTY. No fragments, flat emissions, or emissions that are all unplaced (pure
# background) must ALL give an identically zero log factor -- otherwise the edge moves the
# Li-Stephens model while carrying no phase information.
for case in ("zero_fragments", "flat_emissions", "all_unplaced"):
    v = float(d[case][2])
    if v == 0.0: ok("%s: log psi is identically zero -- the edge is exactly neutral" % case)
    else: no("%s: max |log psi| = %.10g, so an edge with no phase information still moves the "
             "model (sum-one leaks log 4 = 1.3863 here)" % (case, v))
# Mean-one, and swap symmetry, on the informative case.
inf = d["informative"]
if float(inf[2]) > 0.0: ok("informative: log psi is non-zero (max %.2f nats) -- the gate is not vacuous"
                           % float(inf[2]))
else: no("informative: log psi is zero -- the self-test asserts nothing")
usable_cases = [k for k in need if d[k][6] == "1"]
for case in usable_cases:
    if float(d[case][3]) > 1e-9:
        no("%s: mean(exp(log psi)) deviates from 1 by %.2e" % (case, float(d[case][3])))
        break
else:
    ok("every USABLE case is mean-one within each content class (worst dev %.2e)"
       % max(float(d[k][3]) for k in usable_cases))
# A C4-SCALE EDGE MUST COMPLETE. 118 x 119 is 197,177,764 ordered configurations -- the dense table
# refuses it outright -- and the sparse form must build it without truncation or any score cutoff.
lc = {}
for l in open(sys.argv[1]):
    f = l.rstrip("\n").split("\t")
    if len(f) >= 7 and "x" in f[0] and f[0].replace("x", "").isdigit():
        lc[f[0]] = f
if lc:
    for k, f in sorted(lc.items()):
        if f[6] == "ok" and int(f[4]) > 0:
            ok("a %s edge completes: %s configurations -> %s stored classes (%.0fx fewer), %s bytes"
               % (k, f[3], f[4], float(f[3]) / max(1.0, float(f[4])), f[8]))
        else:
            no("a %s edge did not complete: status %s, stored %s" % (k, f[6], f[4]))
        # THE DENSE CAP MUST NOT BE CONSULTED. The same edge built through the dense path with a
        # deliberately tiny configuration cap refuses; the sparse path must be indifferent to it, or
        # production is still gated on the obsolete limit.
        if len(f) > 10 and f[10] == "too-many-configurations" and f[6] == "ok":
            ok("      the same edge refuses at a tiny DENSE cap while sparse succeeds -- production "
               "does not consult it")
        elif len(f) > 10:
            no("dense-cap control for %s: dense %s, sparse %s" % (k, f[10], f[6]))
else:
    no("no large-edge case was reported; the scale claim is untested")

# THE MEAN-ONE BOUND, asserted rather than reasoned about.
import math as _m
# SPARSE psi IS EXACT, NOT AN APPROXIMATION. Its structure removes the dense ordered table:
# a homozygous endpoint has no alternative phase (log psi = 0 by swap invariance), each het x het
# class has only TWO biological phases, so ONE contrast Delta = S_straight - S_crossed reconstructs
# both mean-one ratios; and a class whose four haploid corners all lack in-band mass has both phases
# equal to the same all-background sum, so it is exactly neutral and never stored. None of that is
# threshold pruning -- it follows from the emission support.
for case in usable_cases:
    sd = float(d[case][14]); nchk = int(d[case][15])
    if nchk == 0:
        no("%s: sparse psi was compared on ZERO configurations" % case)
    elif sd > 1e-9:
        no("%s: sparse psi differs from dense by %.4g" % (case, sd))
    else:
        ok("%s: sparse == dense on all %d configurations (%.1e), %s stored vs %s theoretical"
           % (case, nchk, sd, d[case][16], d[case][17]))
# A flat edge has SUPPORT but stores nothing: the corners exist, the phases are equal, Delta is 0.
fl = d.get("flat_emissions")
if fl and int(fl[18]) > 0 and int(fl[16]) == 0:
    ok("flat_emissions has %s support cells yet stores 0 classes -- proven neutral, not pruned"
       % fl[18])
elif fl:
    no("flat_emissions: support %s, stored %s -- a provably neutral class was stored" % (fl[18], fl[16]))

# LAMBDA MUST REACH THE FACTORS, asserted where it is observable: near the background crossover,
# where the mixture is in transition and lambda is not a removable constant.
if "lambda_crossover_lo" in d and "lambda_crossover_hi" in d:
    lo_mx, lo_mn = float(d["lambda_crossover_lo"][10]), float(d["lambda_crossover_lo"][11])
    hi_mx, hi_mn = float(d["lambda_crossover_hi"][10]), float(d["lambda_crossover_hi"][11])
    if abs(hi_mx - lo_mx) > 1e-6 or abs(hi_mn - lo_mn) > 1e-6:
        ok("lambda reaches the factors: at the crossover psi moves [%.4f,%.4f] -> [%.4f,%.4f]"
           % (lo_mn, lo_mx, hi_mn, hi_mx))
    else:
        no("lambda does not change psi even at the background crossover -- it is decoration")
else:
    no("the lambda crossover cases are missing; propagation is not asserted anywhere")
for case in usable_cases:
    mx, mn = float(d[case][10]), float(d[case][11])
    nf, ob = int(d[case][12]), int(d[case][13])
    if nf: no("%s: %d non-finite log psi values" % (case, nf))
    elif ob: no("%s: %d configurations exceed log|C| -- psi is not mean-one" % (case, ob))
    elif mx > _m.log(4.0) + 1e-9:
        no("%s: max log psi %.6f exceeds log 4 = %.6f" % (case, mx, _m.log(4.0)))
    else:
        ok("%s: max log psi %.6f <= log|C| (min %.4g, all finite) -- exp cannot overflow"
           % (case, mx, mn))
for case in usable_cases:
    if float(d[case][4]) > 1e-12:
        no("%s: global homologue swap changes log psi by %.2e" % (case, float(d[case][4]))); break
else:
    ok("global homologue swap leaves log psi unchanged in every case")
# EXPOSURE MUST CANCEL, or the edge is UNSUPPORTED. Equal-length alleles make exposure
# phase-invariant and so cannot detect this; the unequal_exposure cases set it directly. Without the
# requirement, a ZERO-FRAGMENT edge moves the model by 0.105 nats on pure exposure asymmetry.
for case in ("unequal_exposure", "unequal_exposure_with_frags"):
    if d[case][6] == "0" and d[case][8] == "exposure-does-not-cancel":
        ok("%s: edge refused as UNSUPPORTED (exposure asymmetry %s)" % (case, d[case][9]))
    else:
        no("%s: usable=%s status=%s -- a non-cancelling exposure was accepted"
           % (case, d[case][6], d[case][8]))
    if float(d[case][2]) != 0.0:
        no("%s: log psi is non-zero (%.4g) on an unusable edge" % (case, float(d[case][2])))
# A REFUSED EDGE MUST NOT BE INDEXED. score/log_psi are EMPTY on refusal while n_a/n_b remain set,
# so a consumer that checks the wrong thing reads out of bounds -- this segfaulted the self-test the
# moment the refusal case was added, which is why usable() is now the single derived test.
for case in ("unequal_exposure", "invalid_emission", "too_many_configs"):
    if d[case][5] != "-":
        no("%s: class sizes were computed on a refused edge (%s)" % (case, d[case][5]))
        break
else:
    ok("refused edges expose no class structure -- consumers cannot index them by accident")
# AN OWNED FRAGMENT WITH NO EMISSION makes the edge INCOMPLETE; it must never be silently skipped.
iv = d["invalid_emission"]
if iv[6] == "0" and iv[7] == "1" and iv[8] == "invalid-emissions":
    ok("invalid_emission: 1 unformable emission makes the edge INCOMPLETE, not smaller")
else:
    no("invalid_emission: usable=%s n_invalid=%s status=%s -- the fragment was skipped silently"
       % (iv[6], iv[7], iv[8]))
for case in ("zero_fragments", "flat_emissions", "all_unplaced", "informative"):
    if d[case][6] != "1":
        no("%s: a well-formed edge was marked unusable (%s)" % (case, d[case][8])); break
else:
    ok("well-formed edges stay usable; only the refused cases are unusable")
sys.exit(bad)
PYEOFA
  fails=$(( fails + $? ))
else
  bad "linkage self-test produced no output"
fi

# GATE 11: THREE ALLELES PER BLOCK. The 2-allele fixture exercises only the 0/1 cis-trans pair, so
# an indexing defect touching allele 2 or above would pass unseen -- and real blocks have hundreds.
# The sample's true haplotypes are DELIBERATELY absent from the panel: the factor must express a
# combination no single panel haplotype carries, which is the point of factorising at all.
"$PY" - "$OUT" <<'PYEOF8'
import sys, random
out = sys.argv[1]; B = "ACGT"
def seq(n, s):
    r = random.Random(s); return "".join(r.choice(B) for _ in range(n))
def mut(s, k, sd):
    r = random.Random(sd); s = list(s)
    for p in r.sample(range(len(s)), k): s[p] = r.choice([c for c in B if c != s[p]])
    return "".join(s)
X = seq(4000,1); J = seq(40,4); Z = seq(4000,7)
A1 = seq(600,2); A2 = mut(A1,90,3); A3 = mut(A1,90,31)
B1 = seq(600,5); B2 = mut(B1,90,6); B3 = mut(B1,90,61)
segs = [("1",X),("2",A1),("3",A2),("4",A3),("5",J),("6",B1),("7",B2),("8",B3),("9",Z)]
links = [("1","2"),("1","3"),("1","4"),("2","5"),("3","5"),("4","5"),
         ("5","6"),("5","7"),("5","8"),("6","9"),("7","9"),("8","9")]
paths = {"hap1":["1","2","5","6","9"],"hap2":["1","3","5","7","9"],"hap3":["1","4","5","8","9"]}
d = dict(segs)
with open(out+"/a3.gfa","w") as g:
    g.write("H\tVN:Z:1.0\n")
    for n,x in segs: g.write("S\t%s\t%s\n"%(n,x))
    for a,b in links: g.write("L\t%s\t+\t%s\t+\t0M\n"%(a,b))
    for n,st in paths.items(): g.write("P\t%s\t%s\t*\n"%(n,",".join(x+"+" for x in st)))
comp = {"A":"T","C":"G","G":"C","T":"A"}
truth = {"tA3B2": X+A3+J+B2+Z, "tA2B3": X+A2+J+B3+Z}
with open(out+"/a3.r1.fq","w") as f1, open(out+"/a3.r2.fq","w") as f2:
    k = 0
    for nm,h in truth.items():
        for i in range(0,len(h)-360,6):
            a = h[i:i+150]; b = h[i+200:i+350]
            if len(b) < 150: break
            rc = "".join(comp[c] for c in reversed(b))
            f1.write("@%s_%d/1\n%s\n+\n%s\n"%(nm,k,a,"I"*150))
            f2.write("@%s_%d/2\n%s\n+\n%s\n"%(nm,k,rc,"I"*150))
            k += 1
PYEOF8
"$BIN" bubble -i "$OUT/a3.gfa" -r hap1 -o "$OUT/a3b" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" genotype-frag -i "$OUT/a3b.sorted.gfa" -b "$OUT/a3b" -o "$OUT/a3o" \
  -R "$OUT/a3.r1.fq" -R "$OUT/a3.r2.fq" -t 2 --max-divergence 0.05 --fragment-len 350 \
  --fragment-sd 50 --error-rate 0.001 --ownership-table "$OUT/a3.own.tsv" \
  --linkage-potential "$OUT/a3.lpot.tsv" -q >/dev/null 2>&1
if [ -s "$OUT/a3.lpot.tsv.edges.tsv" ]; then
  "$PY" - "$OUT/a3.lpot.tsv.edges.tsv" <<'PYEOF9'
import sys
rows = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])][1:]
bad = 0
def ok(m): print("  ok   " + m)
def no(m):
    global bad; bad += 1; print("  FAIL " + m)
if not rows: no("3-allele fixture produced no edge"); sys.exit(1)
r = rows[0]; na, nb = int(r[2]), int(r[3])
if na >= 3 and nb >= 3: ok("3-allele fixture: %dx%d alleles, %s configurations" % (na, nb, r[6]))
else: no("3-allele fixture collapsed to %dx%d alleles" % (na, nb))
# The truth is (A3,B2)+(A2,B3) = allele indices (2,1) and (1,2), so the winning configuration must
# be a1=2,b1=1,a2=1,b2=2 (or its global swap). Decoding it pins the INDEXING, which a 2-allele
# fixture cannot: index 2 is used at both blocks.
best = int(r[10])
b2 = best % nb; t = best // nb; a2 = t % na; t //= na; b1 = t % nb; a1 = t // nb
got = ((a1,b1),(a2,b2)); want = {((2,1),(1,2)), ((1,2),(2,1))}
if got in want: ok("3-allele fixture recovers the true phase (a1=%d,b1=%d | a2=%d,b2=%d), using "
                   "allele index 2 at both blocks" % (a1,b1,a2,b2))
else: no("3-allele fixture chose (a1=%d,b1=%d | a2=%d,b2=%d); truth is (2,1)|(1,2)" % (a1,b1,a2,b2))
if float(r[12]) < 1e-9: ok("3-allele fixture: swap symmetry holds (%.2e)" % float(r[12]))
else: no("3-allele fixture: swap asymmetry %.4g" % float(r[12]))
sys.exit(bad)
PYEOF9
  fails=$(( fails + $? ))
else
  bad "3-allele fixture produced no edge aggregate"
fi

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

# ---------------------------------------------------------------------------------------------
# GATE 23: THE END-TO-END HYBRID CALL. The sentinel that used to stand here checked
# `genotype-frag --help`, but the entry point landed on `genotype`, where the architecture put it --
# so it guarded a flag that could never appear and had gone vacuous. It is replaced by the real
# assertions, and the flag is looked for on the command that actually has it.
if ! "$BIN" genotype --help 2>&1 | grep -q -- "--hybrid-call"; then
  bad "--hybrid-call is missing from the genotype command; the end-to-end gates cannot run"
else
  # 1. NO FLAG: byte-identical legacy output. Structural, not a tolerance -- with no linkage edges
  #    the kernel takes the factorised path at every edge.
  "$BIN" genotype -i "$OUT/b.sorted.gfa" -b "$OUT/b" -r hapAA -o "$OUT/e2e.leg" \
    -R "$OUT/r1.fq" -R "$OUT/r2.fq" --fragment-len 350 -q >/dev/null 2>&1
  "$BIN" genotype -i "$OUT/b.sorted.gfa" -b "$OUT/b" -r hapAA -o "$OUT/e2e.hyb" \
    -R "$OUT/r1.fq" -R "$OUT/r2.fq" --fragment-len 350 --hybrid-call -q >/dev/null 2>&1
  if [ -s "$OUT/e2e.leg.genotypes.tsv" ] && [ -s "$OUT/e2e.hyb.genotypes.tsv" ]; then
    ok "both arms produced a call table"
    # 2. CONTENT PRESERVED, PHASE CORRECTED. The fixture is built so per-block allele content is
    #    identical under either phase; only a junction-spanning fragment can decide between them.
    "$PY" - "$OUT/e2e.leg.genotypes.tsv" "$OUT/e2e.hyb.genotypes.tsv" <<'PYEOFI'
import sys
def load(p):
    rows = []
    for i, l in enumerate(open(p)):
        f = l.rstrip("\n").split("\t")
        if i == 0: continue
        rows.append(f)
    return rows
L, H = load(sys.argv[1]), load(sys.argv[2])
bad = 0
def ok(m): print("  ok   " + m)
def no(m):
    global bad; bad += 1; print("  FAIL " + m)
if len(L) != len(H): no("block counts differ: %d vs %d" % (len(L), len(H))); sys.exit(1)
hdr = open(sys.argv[1]).readline().rstrip("\n").split("\t")
if len(hdr) > 11 and hdr[11] == "best_ordered_pair_posterior":
    ok("column 12 is best_ordered_pair_posterior, and is described as such")
else:
    no("column 12 is %s -- the test's description no longer matches the schema"
       % (hdr[11] if len(hdr) > 11 else "absent"))
content_same = all(l[7] == h[7] and l[8] == h[8] for l, h in zip(L, H))
if content_same:
    ok("block CONTENT is unchanged by linkage (allele1/allele2 identical at every block)")
else:
    no("linkage changed block content: %s" % [(l[0], l[7], l[8], h[7], h[8])
                                              for l, h in zip(L, H) if l[7] != h[7] or l[8] != h[8]])
lp = {(l[9], l[10]) for l in L}
hp = {(h[9], h[10]) for h in H}
truth = {("hapAB", "hapBA"), ("hapBA", "hapAB")}
if lp & truth:
    no("the legacy caller already recovered the phase; the gate would prove nothing")
else:
    ok("the marker caller alone gets the phase WRONG (%s) -- as the fixture intends" % sorted(lp))
if hp <= truth and hp:
    ok("the hybrid caller recovers the CORRECT phase (%s)" % sorted(hp))
else:
    no("hybrid phase is %s, expected hapAB/hapBA" % sorted(hp))
# Column 12 is best_ordered_pair_posterior: the largest ORDERED state, NOT the unordered diplotype
# posterior, which would aggregate both homologue orders and any sequence-equivalent
# representatives. Describing it as "the pair posterior" conflated the two. The phase-recovery
# result above does not depend on this number -- the NAMED pair is what changed -- but the quantity
# must be labelled for what it is.
try:
    lm = max(float(l[11]) for l in L); hm = max(float(h[11]) for h in H)
    if hm > lm:
        ok("best_ORDERED_pair posterior rises %.3f -> %.3f (not the unordered diplotype posterior, "
           "which is not emitted)" % (lm, hm))
    else:
        no("best_ordered_pair_posterior did not rise: %.3f -> %.3f" % (lm, hm))
except (ValueError, IndexError):
    no("could not read the best_ordered_pair_posterior column")
sys.exit(bad)
PYEOFI
    fails=$(( fails + $? ))
  else
    bad "an end-to-end arm produced no call table"
  fi
  # 2b. THE WORK GUARDS are operational, not statistical: counts of what the support search does,
  #     no likelihood anywhere in them. There are now TWO, because after direct positional
  #     verification one number cannot name both -- proposed cells are allele pairs retained,
  #     full-read verifications are positional starts actually compared. Exceeding either must
  #     refuse TRANSACTIONALLY: nothing subtracted, legacy call intact, rather than appearing to
  #     hang. Both are charged BEFORE the work; tests/genotype_hybrid_chain's unit companion
  #     (--budget-selftest) proves the ordering, this proves the transaction.
  for guard in "--hybrid-max-proposed-cells 1" "--hybrid-max-full-read-verifications 1"; do
  "$BIN" genotype -i "$OUT/b.sorted.gfa" -b "$OUT/b" -r hapAA -o "$OUT/e2e.wg" \
    -R "$OUT/r1.fq" -R "$OUT/r2.fq" --fragment-len 350 --hybrid-call \
    $guard --hybrid-status "$OUT/wg.tsv" -q >/dev/null 2>&1
  if [ -s "$OUT/wg.tsv" ] && [ -s "$OUT/e2e.wg.genotypes.tsv" ]; then
    WST=$(awk -F'\t' '$1=="hybrid_status"{print $2}' "$OUT/wg.tsv")
    WEX=$(awk -F'\t' '$1=="fragments_excluded"{print $2}' "$OUT/wg.tsv")
    WAC=$(awk -F'\t' '$1=="active_edges"{print $2}' "$OUT/wg.tsv")
    WRE=$(awk -F'\t' '$1=="reason"{print $2}' "$OUT/wg.tsv")
    if [ "$WST" = "INCOMPLETE" ] && [ "${WEX:-1}" = "0" ] && [ "${WAC:-1}" = "0" ]; then
      ok "$guard refuses transactionally: INCOMPLETE, 0 active, 0 excluded"
    else
      bad "work guard $guard: status=$WST active=$WAC excluded=$WEX"
    fi
    case "$WRE" in
      support-search-*-limit*) ok "the refusal names WHICH limit ($(echo "$WRE" | cut -c1-46)...)" ;;
      *) bad "the work refusal lost its reason: '$WRE'" ;;
    esac
    # And the legacy call must be untouched by a refusal.
    if cmp -s "$OUT/e2e.leg.genotypes.tsv" "$OUT/e2e.wg.genotypes.tsv"; then
      ok "a $guard refusal leaves the legacy call byte-identical"
    else
      bad "a $guard refusal changed the legacy call"
    fi
  else
    bad "the work-guard arm produced nothing for $guard"
  fi
  done
  # THE DEPRECATED ALIAS still refuses, and says it is deprecated rather than silently meaning
  # something new.
  "$BIN" genotype -i "$OUT/b.sorted.gfa" -b "$OUT/b" -r hapAA -o "$OUT/e2e.al" \
    -R "$OUT/r1.fq" -R "$OUT/r2.fq" --fragment-len 350 --hybrid-call \
    --hybrid-max-verified-windows 1 --hybrid-status "$OUT/al.tsv" > "$OUT/al.log" 2>&1
  AST=$(awk -F'\t' '$1=="hybrid_status"{print $2}' "$OUT/al.tsv" 2>/dev/null)
  if [ "$AST" = "INCOMPLETE" ] && grep -qi "deprecated" "$OUT/al.log"; then
    ok "--hybrid-max-verified-windows still refuses, and is reported as DEPRECATED"
  else
    bad "the deprecated alias: status=$AST, deprecation notice $(grep -ci deprecated "$OUT/al.log")"
  fi
  # THE WORK ACTUALLY DONE is reported beside the limits it was held to -- on an UNREFUSED run,
  # since a refused one has by construction stopped counting.
  "$BIN" genotype -i "$OUT/b.sorted.gfa" -b "$OUT/b" -r hapAA -o "$OUT/e2e.wk" \
    -R "$OUT/r1.fq" -R "$OUT/r2.fq" --fragment-len 350 --hybrid-call \
    --hybrid-status "$OUT/wk.tsv" -q >/dev/null 2>&1
  WP=$(awk -F'\t' '$1=="work_proposed_cells"{print $2}' "$OUT/wk.tsv" 2>/dev/null)
  WV=$(awk -F'\t' '$1=="work_full_read_verifications"{print $2}' "$OUT/wk.tsv" 2>/dev/null)
  WB=$(awk -F'\t' '$1=="work_bases_compared_upper_bound"{print $2}' "$OUT/wk.tsv" 2>/dev/null)
  if [ "${WP:-0}" -gt 0 ] && [ "${WV:-0}" -gt 0 ] && [ "${WB:-0}" -ge "${WV:-1}" ]; then
    ok "an unrefused run reports its work: $WP cells, $WV verifications, $WB bases (upper bound)"
  else
    bad "work counters not reported: cells=$WP verifications=$WV bases=$WB"
  fi

  # 3. THE LINKAGE PARAMETER CONTRACT. Every parameter reported with the result, and actually USED
  #    rather than echoed -- changing lambda must change the answer, or the report is decoration.
  "$BIN" genotype -i "$OUT/b.sorted.gfa" -b "$OUT/b" -r hapAA -o "$OUT/e2e.p1" \
    -R "$OUT/r1.fq" -R "$OUT/r2.fq" --fragment-len 350 --hybrid-call \
    --hybrid-status "$OUT/p1.tsv" -q >/dev/null 2>&1
  "$BIN" genotype -i "$OUT/b.sorted.gfa" -b "$OUT/b" -r hapAA -o "$OUT/e2e.p2" \
    -R "$OUT/r1.fq" -R "$OUT/r2.fq" --fragment-len 350 --hybrid-call --hybrid-lambda 0.4 \
    --hybrid-status "$OUT/p2.tsv" -q >/dev/null 2>&1
  if [ -s "$OUT/p1.tsv" ] && [ -s "$OUT/p2.tsv" ]; then
    MISSING=""
    for k in param_lambda param_lambda_source param_bg_divergence param_outlier_mix \
             param_error_rate param_max_divergence param_fragment_len param_fragment_sd \
             param_discordant_rate param_insert_sigmas; do
      grep -q "^$k	" "$OUT/p1.tsv" || MISSING="$MISSING $k"
    done
    [ -z "$MISSING" ] && ok "every linkage parameter is reported with the result" \
                      || bad "linkage parameters not reported:$MISSING"
    L1=$(awk -F'\t' '$1=="param_lambda"{print $2}' "$OUT/p1.tsv")
    S1=$(awk -F'\t' '$1=="param_lambda_source"{print $2}' "$OUT/p1.tsv")
    L2=$(awk -F'\t' '$1=="param_lambda"{print $2}' "$OUT/p2.tsv")
    S2=$(awk -F'\t' '$1=="param_lambda_source"{print $2}' "$OUT/p2.tsv")
    [ "$S1" = "default" ] && [ "$S2" = "supplied" ] \
      && ok "lambda records its provenance (default $L1, supplied $L2)" \
      || bad "lambda source not recorded: '$S1' then '$S2'"
    # USED, NOT MERELY ECHOED -- but asserted where it is OBSERVABLE. lambda cancels wherever every
    # configuration sits far above or far below the background: there it is a constant per fragment
    # and the mean-one centering removes it exactly. This fixture is saturated (contrasts ~160 nats
    # per fragment), and measured: lambda over 0.0005..5.0, a 10,000x range, leaves the call table
    # byte-identical. That is the model behaving correctly, not the parameter being ignored, so the
    # propagation gate lives at the crossover instead -- see the lambda_crossover cases below.
    if cmp -s "$OUT/e2e.p1.genotypes.tsv" "$OUT/e2e.p2.genotypes.tsv"; then
      ok "on this SATURATED fixture lambda cancels (mean-one centering removes a per-fragment"
      ok "      constant); propagation is asserted at the background crossover instead"
    else
      ok "changing --hybrid-lambda changes the call table"
    fi
  else
    bad "hybrid status report not written"
  fi

  # 4. THE EXCLUSION IDENTITY, end to end: reads removed == 2 x fragments consumed by active edges.
  HL=$("$BIN" genotype -i "$OUT/b.sorted.gfa" -b "$OUT/b" -r hapAA -o "$OUT/e2e.h2" \
       -R "$OUT/r1.fq" -R "$OUT/r2.fq" --fragment-len 350 --hybrid-call 2>&1)
  NEX=$(echo "$HL" | sed -nE 's/.*hybrid: .*, ([0-9]+) excluded from markers.*/\1/p')
  NRD=$(echo "$HL" | sed -nE 's/.*marker exclusion: ([0-9]+) fragment\(s\) named, ([0-9]+) reads.*/\2/p')
  NAC=$(echo "$HL" | sed -nE 's/.*hybrid: [A-Z]+, [0-9]+ fragments, [0-9]+ candidate edge\(s\), ([0-9]+) active.*/\1/p')
  STAT=$(echo "$HL" | sed -nE 's/.*hybrid: ([A-Z]+),.*/\1/p')
  [ "${STAT:-}" = "COMPLETE" ] && ok "hybrid status is COMPLETE on the fixture" \
                               || bad "hybrid status is ${STAT:-none}"
  [ "${NAC:-0}" -ge 1 ] && ok "$NAC linkage edge(s) active" || bad "no linkage edge activated"
  if [ -n "${NEX:-}" ] && [ -n "${NRD:-}" ] && [ "$NRD" = "$(( NEX * 2 ))" ]; then
    ok "exclusion is exact end to end: $NEX fragments -> $NRD reads (both mates, no more)"
  else
    bad "excluded $NEX fragments but $NRD reads left the counts (expected $(( ${NEX:-0} * 2 )))"
  fi
fi

echo
if [ "$fails" -eq 0 ]; then echo "hybrid chain: all active assertions passed"; else
  echo "hybrid chain: $fails assertion(s) failed"; fi
exit "$fails"
