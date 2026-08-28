#!/usr/bin/env bash
# genotype_frag_reference.sh - the exact reference scorer against analytically checkable fixtures.
#
#   genotype_frag_reference.sh <panvar-binary> <out-dir>
#
# The reference implements docs/reports/genotype-fragment-model-contract.md directly: every fragment
# start on both haplotypes, the insert prior summed over, both strands, and no syncmer index, anchor
# cap or placement_topk anywhere. These fixtures pin it to values that are derivable rather than
# recorded, so that it can then serve as the oracle the fast path is tested against.
#
# Each case states WHAT IT WOULD CATCH, because a fixture whose failure mode is unclear is a fixture
# nobody can act on.
set -uo pipefail
BIN="${1:?usage: genotype_frag_reference.sh <panvar> <outdir>}"
OUT="${2:?}"; mkdir -p "$OUT"; OUT="$OUT/run.$$"; rm -rf "$OUT"; mkdir -p "$OUT"
fails=0
ok()  { printf "  ok   %s\n" "$1"; }
bad() { printf "  FAIL %s\n" "$1"; fails=$((fails+1)); }
PY="${PYTHON:-python3}"
P="--haploid-depth 0.05 --fragment-len 350 --fragment-sd 50 --error-rate 0.01"

if ! "$BIN" genotype-frag --help >/dev/null 2>&1; then
  echo "  FAIL binary has no genotype-frag; every assertion below would be vacuous"; exit 1
fi

score() { "$BIN" genotype-frag --reference-score "$1" "$2" -R "$3" $P 2>/dev/null; }
# a > b, by more than tol
gt() { "$PY" -c "import sys; a,b,t=float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]); sys.exit(0 if a-b>t else 1)" "$1" "$2" "${3:-0.5}"; }
near() { "$PY" -c "import sys; a,b,t=float(sys.argv[1]),float(sys.argv[2]),float(sys.argv[3]); sys.exit(0 if abs(a-b)<=t else 1)" "$1" "$2" "$3"; }

"$PY" - "$OUT" <<'PYEOF'
import random, sys, os
out = sys.argv[1]
random.seed(11)
def rnd(n): return ''.join(random.choice("ACGT") for _ in range(n))
def rc(s):  return s.translate(str.maketrans("ACGT","TGCA"))[::-1]
def wr(name, seq): open(os.path.join(out, name), "w").write(f">{name}\n{seq}\n")

CORE = rnd(1400)
wr("A.fa", CORE)
wr("Arc.fa", rc(CORE))

# one substitution every 70 bp -> a clearly different haplotype of the SAME length
B = list(CORE)
for i in range(0, len(B), 70): B[i] = "ACGT"[("ACGT".index(B[i]) + 1) % 4]
wr("B.fa", "".join(B))

# a 300 bp segment duplicated in tandem: same content, longer, higher copy number
SEG, PRE, POST = CORE[400:700], CORE[:400], CORE[700:]
wr("dup.fa", PRE + SEG + SEG + POST)

# a 300 bp deletion: shorter, strictly less content
wr("del.fa", PRE + POST)

# same TOTAL length as CORE but wrong composition: swap two 200 bp blocks' contents for noise
SW = CORE[:500] + rnd(200) + CORE[700:]
wr("swap.fa", SW)

def frags(seq, tag, step=25, rl=120, ins=350, n=None):
    out_r = []
    for i in range(0, len(seq) - ins, step):
        out_r.append(f">{tag}{i}/1\n{seq[i:i+rl]}\n>{tag}{i}/2\n{rc(seq[i+ins-rl:i+ins])}\n")
        if n and len(out_r) >= n: break
    return out_r

open(os.path.join(out, "reads_core.fa"), "w").write("".join(frags(CORE, "c")))
open(os.path.join(out, "reads_dup.fa"),  "w").write("".join(frags(PRE + SEG + SEG + POST, "d")))
# reads from CORE plus a handful of fragments from unrelated sequence (background)
open(os.path.join(out, "reads_bg.fa"), "w").write(
    "".join(frags(CORE, "c")) + "".join(frags(rnd(1400), "x", n=5)))
print(sum(1 for l in open(os.path.join(out,"reads_core.fa")) if l.startswith(">")) // 2)
PYEOF

R="$OUT/reads_core.fa"

# ---------------------------------------------------------------- the truth outranks alternatives
AA=$(score "$OUT/A.fa" "$OUT/A.fa" "$R")
AB=$(score "$OUT/A.fa" "$OUT/B.fa" "$R")
BB=$(score "$OUT/B.fa" "$OUT/B.fa" "$R")
if gt "$AA" "$AB" && gt "$AB" "$BB"; then
  ok "homozygous truth > one wrong homologue > both wrong ($AA > $AB > $BB)"
else
  bad "ordering wrong: A/A=$AA A/B=$AB B/B=$BB"
fi

# ------------------------------------------------------------------ strand symmetry of the model
# Would catch: a scorer that only tries the FR orientation. Real panel sequence needs this -- one of
# cyp2d6 HG04036's haplotypes is spelled reverse-complemented by block concatenation.
ARC=$(score "$OUT/Arc.fa" "$OUT/Arc.fa" "$R")
near "$AA" "$ARC" 0.001 \
  && ok "a haplotype and its reverse complement score identically ($AA vs $ARC)" \
  || bad "reverse-complemented haplotype scores differently: $AA vs $ARC"

# ------------------------------------------------------------- homozygous doubling, analytically
# Contract: a homozygous pair doubles BOTH exposure and the placement set. Scored against an empty
# second homologue the difference must be exactly -lambda*starts(A) + n_placed*log(2).
# Would catch: doubling exposure without doubling intensity -- the log(2)-per-fragment bug.
printf '>empty\n\n' > "$OUT/empty.fa"
A0=$(score "$OUT/A.fa" "$OUT/empty.fa" "$R")
NF=$(grep -c '/1$' "$R")
"$PY" - "$AA" "$A0" "$NF" <<'PYEOF'
import sys, math
aa, a0, nf = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
# Exposure is over (start, L), so E = SUM_L pi(L) * (N - L + 1), with pi the NORMALISED
# concordant/discordant mixture over the range actually summed. Computed here rather than
# hard-coded, so that the prediction tracks the contract instead of a previous implementation.
N, mu, sd, disc = 1400, 350.0, 50.0, 0.01
lo, hi = int(mu - 4*sd), int(mu + 4*sd)
span = hi - lo + 1
w = []
for L in range(lo, hi + 1):
    z = (L - mu) / sd
    conc = (1 - disc) * math.exp(-0.5*z*z) / (sd * math.sqrt(2*math.pi))
    w.append(conc + disc / span)
tot = sum(w)
E = sum((wi / tot) * max(0, N - L + 1) for wi, L in zip(w, range(lo, hi + 1)))
expected = -0.05 * E + nf * math.log(2)
got = aa - a0
print(f"  .... exposure {E:.1f} start-length states; predicted {expected:.3f}, observed {got:.3f}")
sys.exit(0 if abs(expected - got) < 0.6 else 1)
PYEOF
[ $? -eq 0 ] && ok "homozygous doubling matches the analytic prediction (-lambda*starts + n*log2)" \
             || bad "homozygous doubling does not match the analytic prediction"

# --------------------------------------------------------------- copy number from placements alone
# Reads from a 2-copy haplotype: the 2-copy candidate must beat the 1-copy one. Nothing in the model
# is a depth term -- this must fall out of placement multiplicity and exposure by itself.
# Would catch: a model that cannot express dosage without a bolted-on coverage channel.
RD="$OUT/reads_dup.fa"
D2=$(score "$OUT/dup.fa" "$OUT/dup.fa" "$RD")
D1=$(score "$OUT/A.fa"   "$OUT/A.fa"   "$RD")
gt "$D2" "$D1" \
  && ok "reads from a duplicated segment prefer the 2-copy haplotype ($D2 > $D1), with no depth term" \
  || bad "the 2-copy haplotype does not win on its own reads: $D2 vs $D1"

# --------------------------------------------------------------------- deletion is not free either
DEL=$(score "$OUT/del.fa" "$OUT/del.fa" "$R")
gt "$AA" "$DEL" \
  && ok "a 300 bp deletion is rejected on reads that span it ($AA > $DEL)" \
  || bad "the deletion haplotype was not penalised: $AA vs $DEL"

# ------------------------------------------------- right total length, wrong composition, rejected
# Would catch: a likelihood that has drifted into scoring total length rather than sequence.
SWP=$(score "$OUT/swap.fa" "$OUT/swap.fa" "$R")
gt "$AA" "$SWP" \
  && ok "equal total length with wrong composition is rejected ($AA > $SWP)" \
  || bad "a same-length wrong-composition haplotype was not penalised: $AA vs $SWP"

# ----------------------------------------------------------------- background fragments are bounded
# Five fragments of unrelated sequence must cost the BOUND, not the emission they would otherwise
# earn. The bound is derivable rather than a taste threshold: a fragment the pair cannot explain
# falls to eta * B_f, so it costs (placed - floor) each, about 107 nats for a 240 bp fragment at
# eta=0.05, bg_divergence=0.10, error_rate=0.01. Unbounded it would be roughly 800 nats each.
# Would catch: an unbounded emission -- the defect where one foreign fragment outvoted hundreds of
# ordinary ones -- and equally a background so generous that unexplained fragments become free.
BG=$(score "$OUT/A.fa" "$OUT/A.fa" "$OUT/reads_bg.fa")
"$PY" - "$AA" "$BG" <<'PYEOF'
import sys, math
aa, bg = float(sys.argv[1]), float(sys.argv[2])
# error_rate is the TOTAL substitution probability, so a specific mismatching base has eps/3 --
# the same convention the scorer and the simulator use.
eta, eps, bgdiv, n = 0.05, 0.01, 0.10, 240
e_bg = int(bgdiv * n)
floor = math.log(eta) + e_bg*math.log(eps/3) + (n-e_bg)*math.log1p(-eps)
placed = math.log(1-eta) + math.log(0.05) + math.log(0.5) - 4.83
predicted = 5 * (placed - floor)
e_rand = int(0.75 * n)
unbounded = 5 * (placed - (e_rand*math.log(eps/3) + (n-e_rand)*math.log1p(-eps)))
drop = aa - bg
print(f"  .... five unexplained fragments cost {drop:.1f} nats; bound predicts ~{predicted:.0f}, "
      f"unbounded would be ~{unbounded:.0f}")
sys.exit(0 if abs(drop - predicted) < 0.25 * predicted else 1)
PYEOF
[ $? -eq 0 ] && ok "unexplained fragments cost the derived bound, not the emission they would earn" \
             || bad "background cost does not match the derived bound"

# ------------------------------------------------- the emission convention, with REAL mismatches
# Every other fixture here is error-free, and an error-free fixture cannot detect a wrong per-mismatch
# constant: with zero edits, log(eps) and log(eps/3) give the same score. Measured consequence of that
# blind spot: the reference and the accelerated scorer disagreed on this constant for a whole round of
# experiments, the reconciliation still closed -- it sums each model's own contributions -- and a
# rescue result was reported from a comparison of two different likelihoods.
#
# One read placed at a known position with a known number of mismatches. Its contribution is derivable:
#   log[(1-eta)*lambda*(1/2)*pi(L)*eps3^m*(1-eps)^(n-m) + eta*Pbg]   summed over the placements it has.
# Asserting the SIGN of the per-mismatch slope is enough to pin the constant: adding one mismatch must
# cost log((1-eps)/eps3) = log(0.99/0.003333) = 5.69 nats, not log(0.99/0.01) = 4.60.
"$PY" - "$OUT" <<'PYEOF'
import sys, os, random
out = sys.argv[1]
random.seed(3)
core = "".join(random.choice("ACGT") for _ in range(1200))
open(os.path.join(out, "mm_hap.fa"), "w").write(">h\n" + core + "\n")
def rc(s): return s.translate(str.maketrans("ACGT","TGCA"))[::-1]
def sub(c): return random.choice([b for b in "ACGT" if b != c])
# the SAME fragment, placed identically, differing only in how many bases mismatch
for m in (0, 1):
    st, ins, rl = 300, 350, 120
    r1 = list(core[st:st+rl]); r2 = rc(core[st+ins-rl:st+ins])
    for i in range(m):
        r1[10+i] = sub(r1[10+i])
    open(os.path.join(out, "mm%d.fa" % m), "w").write(
        ">f/1\n" + "".join(r1) + "\n>f/2\n" + r2 + "\n")
PYEOF
S0=$(score "$OUT/mm_hap.fa" "$OUT/mm_hap.fa" "$OUT/mm0.fa")
S1=$(score "$OUT/mm_hap.fa" "$OUT/mm_hap.fa" "$OUT/mm1.fa")
"$PY" - "$S0" "$S1" <<'PYEOF'
import sys, math
s0, s1 = float(sys.argv[1]), float(sys.argv[2])
eps = 0.01
cost = s0 - s1
want_third = math.log((1-eps) / (eps/3))   # 5.69 -- a specific mismatching base
want_plain = math.log((1-eps) / eps)       # 4.60 -- "some substitution", the wrong convention
print(f"  .... one mismatch costs {cost:.3f} nats; eps/3 predicts {want_third:.3f}, eps predicts {want_plain:.3f}")
sys.exit(0 if abs(cost - want_third) < 0.15 else 1)
PYEOF
[ $? -eq 0 ] && ok "a mismatch costs log((1-eps)/(eps/3)), the specific-base convention" \
             || bad "the per-mismatch constant does not match the eps/3 convention"

echo
if [ "$fails" -eq 0 ]; then echo "reference scorer: all assertions passed"; else
  echo "reference scorer: $fails assertion(s) failed"; fi
exit "$fails"
