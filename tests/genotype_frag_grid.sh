#!/usr/bin/env bash
# genotype_frag_grid.sh - does multiplicity-aware compression scale, as well as being correct?
#
#   genotype_frag_grid.sh <panvar-binary> <out-dir> [copy-counts] [repeat-lengths] [seeds]
#
# Correctness is settled on one fixture: grouping adds no measurable error and retains all placement
# mass where the anchor cap keeps 4% and top-k 7%. Tractability is NOT settled. The compression is
# DOWNSTREAM -- anchors are expanded, every recruited placement is aligned, and mate combinations are
# formed, and only then are equal likelihoods grouped. Measured on a 10-copy array it evaluates 171x
# more mate combinations than cap 8 + top-k 2 while storing 27x fewer placements. Mate combinations
# grow as (copies)^2 per fragment, so this grid exists to find where that stops being affordable.
#
# Reads are STOCHASTIC and drawn from the same model the exact scorer assumes: uniform start, insert
# from the library prior, either strand, per-base errors at the declared rate. Deterministic evenly
# spaced reads make a fixture's margin an artefact of the spacing.
#
# For near-ties the criterion is NOT an identical winner: the compressed optimum must lie within a
# predeclared equivalence tolerance of the reference best, and the score error must meet its bound.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/lib/experiment.sh"
BIN="${1:?usage: genotype_frag_grid.sh <panvar> <outdir> [copies] [replens] [seeds]}"
OUT="${2:?}"; COPIES="${3:-2,4,8}"; REPLENS="${4:-100}"; SEEDS="${5:-2}"
mkdir -p "$OUT"; OUT="$OUT/run.$$"; rm -rf "$OUT"; mkdir -p "$OUT"
PY="${PYTHON:-python3}"
exp_init "frag-grid" "$BIN" "$OUT/provenance.txt"
EQ_TOL=1.0        # predeclared: the compressed optimum must be within this of the reference best
MASS_BOUND=1e-3   # predeclared: omitted placement mass

printf "  equivalence tolerance %s nats, mass bound %s; reads stochastic from the scorer's own model\n\n" "$EQ_TOL" "$MASS_BOUND"
printf "  %-5s %-7s %-5s | %-10s %-9s %-8s %-7s | %-9s %-9s %-6s\n" \
       copies replen seed "ref_best" "err" "in_class" "mass%" "combos" "groups" "sec"

for CN in ${COPIES//,/ }; do
for RL in ${REPLENS//,/ }; do
for SD in $(seq 1 "$SEEDS"); do
  W="$OUT/c${CN}_r${RL}_s${SD}"; mkdir -p "$W"
  "$PY" - "$W" "$CN" "$RL" "$SD" <<'PYEOF'
import random, sys, os
out, cn, rl, seed = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
random.seed(1000 + seed)
def rnd(n): return "".join(random.choice("ACGT") for _ in range(n))
def rc(s): return s.translate(str.maketrans("ACGT","TGCA"))[::-1]
L, N, SPC, SPCALT = rnd(200), rnd(200), rnd(100), rnd(100)
SEG = rnd(rl)
A_many, A_few, A_ref = SEG*cn + SPC, SEG*(cn-1) + SPC, SEG + SPCALT
many, few, ref = L+A_many+N, L+A_few+N, L+A_ref+N
for nm, sq in (("many",many),("fewer",few),("ref",ref)):
    open(os.path.join(out,f"{nm}.fa"),"w").write(f">{nm}\n{sq}\n")
with open(os.path.join(out,"g.gfa"),"w") as g:
    g.write("H\tVN:Z:1.0\n")
    for i,sq in enumerate((L, A_ref, A_few, A_many, N), start=1): g.write(f"S\t{i}\t{sq}\n")
    for a in (2,3,4): g.write(f"L\t1\t+\t{a}\t+\t0M\nL\t{a}\t+\t5\t+\t0M\n")
    g.write("P\tref\t1+,2+,5+\t*\nP\tfewer\t1+,3+,5+\t*\nP\tmany\t1+,4+,5+\t*\n")
# STOCHASTIC reads from the scorer's own model: uniform start, insert ~ N(150,20), either strand,
# per-base error 0.01. Homozygous for `many`, so two independent homologue streams.
mu, sd_i, rlen, eps, lam = 150.0, 20.0, 60, 0.01, 0.05
recs = []
for hom in (0, 1):
    n_frag = int(lam * max(0, len(many) - int(mu) + 1))
    for i in range(n_frag):
        ins = max(2*rlen, int(random.gauss(mu, sd_i)))
        if len(many) - ins < 1: continue
        st = random.randrange(0, len(many) - ins + 1)
        frag = many[st:st+ins]
        r1, r2 = frag[:rlen], rc(frag[-rlen:])
        if random.random() < 0.5: r1, r2 = rc(frag[-rlen:]), frag[:rlen]
        def err(x): return "".join(random.choice("ACGT") if random.random() < eps else c for c in x)
        recs.append(f">h{hom}_{i}/1\n{err(r1)}\n>h{hom}_{i}/2\n{err(r2)}\n")
open(os.path.join(out,"reads.fa"),"w").write("".join(recs))
PYEOF
  "$BIN" bubble -i "$W/g.gfa" -r ref -o "$W/bub" --min-variant-bp 0 -q >/dev/null 2>&1
  P="--haploid-depth 0.05 --fragment-len 150 --fragment-sd 20 --error-rate 0.01"
  : > "$W/ref.tsv"
  for i in ref fewer many; do for j in ref fewer many; do
    [[ "$i" > "$j" ]] && continue
    v=$("$BIN" genotype-frag --reference-score "$W/$i.fa" "$W/$j.fa" -R "$W/reads.fa" $P 2>/dev/null)
    printf '%s/%s\t%s\n' "$i" "$j" "$v" >> "$W/ref.tsv"
  done; done
  t0=$(date +%s)
  /usr/bin/time -l "$BIN" genotype-frag -i "$W/g.gfa" -b "$W/bub" -o "$W/f" -R "$W/reads.fa" \
    --haplotype-mode --joint-marginal --hamming-emission --joint-top-pairs 0 --top-pairs 20 \
    --rung-zero --multiplicity-aware --mass-tolerance "$MASS_BOUND" $P > "$W/log.txt" 2> "$W/time.txt" || true
  t1=$(date +%s)
  "$PY" - "$W" "$CN" "$RL" "$SD" "$EQ_TOL" "$((t1-t0))" <<'PYEOF'
import sys, os, re
W, cn, rl, sd, tol, secs = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], float(sys.argv[5]), sys.argv[6]
ref = {}
for l in open(os.path.join(W,"ref.tsv")):
    k,v = l.rstrip("\n").split("\t")
    try: ref[k] = float(v)
    except ValueError: pass
if not ref: print(f"  {cn:<5} {rl:<7} {sd:<5} | reference produced nothing"); sys.exit()
best = max(ref.values())
top_class = {k for k,v in ref.items() if best - v <= tol}
pf = os.path.join(W,"f.hap_pairs.tsv")
if not os.path.exists(pf): print(f"  {cn:<5} {rl:<7} {sd:<5} | accelerated run produced nothing"); sys.exit()
rows = [l.rstrip("\n").split("\t") for l in open(pf)][1:]
if not rows: print(f"  {cn:<5} {rl:<7} {sd:<5} | no pairs"); sys.exit()
a,b = rows[0][1], rows[0][2]
key = "/".join(sorted((a,b)))
got = float(rows[0][3])
refv = ref.get(key, float("nan"))
err = abs(got - refv) if refv == refv else float("nan")
log = open(os.path.join(W,"log.txt")).read() + open(os.path.join(W,"time.txt")).read()
def grab(pat, d="?"):
    m = re.search(pat, log); return m.group(1) if m else d
combos = grab(r"(\d+) mate combinations")
groups = grab(r"(\d+) groups after")
massmax = grab(r"max ([0-9.e+-]+),")
print(f"  {cn:<5} {rl:<7} {sd:<5} | {best:10.2f} {err:9.3f} {('yes' if key in top_class else 'NO'):<8} "
      f"{massmax:<7} | {combos:<9} {groups:<9} {secs:<6}")
PYEOF
done; done; done
echo
echo "  in_class = the accelerated optimum lies within the equivalence tolerance of the reference best."
echo "  combos is the tractability number: it grows as (copies)^2 per fragment and is the work the"
echo "  anchor cap and top-k were introduced to avoid."
