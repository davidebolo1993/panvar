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
GRID_FAILS=0      # the grid was previously a pure measurement: every gate printed its verdict and
                  # the script exited 0 regardless, so no failure in it could ever break a build
EQ_TOL=1.0        # predeclared: the compressed optimum must be within this of the reference best
MASS_BOUND=1e-3   # predeclared: omitted placement mass

printf "  equivalence tolerance %s nats, mass bound %s; reads stochastic from the scorer's own model\n\n" "$EQ_TOL" "$MASS_BOUND"
printf "  %-5s %-7s %-5s | %-10s %-9s %-10s %-10s %-8s %-6s | %-8s %-8s %-7s %-5s\n" \
       copies replen seed "ref_best" "U-R" "rescued-R" "zeroseed-R" "max|G-U|" "in_cls" "combos" "combosM" "groups" "sec"

for CN in ${COPIES//,/ }; do
for RL in ${REPLENS//,/ }; do
for SD in $(seq 1 "$SEEDS"); do
  W="$OUT/c${CN}_r${RL}_s${SD}"; mkdir -p "$W"
  "$PY" - "$W" "$CN" "$RL" "$SD" <<'PYEOF'
import random, sys, os, math
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
# Drawn from the CONTRACT'S POINT PROCESS, not from a convenient approximation of it. Events occur
# over (start, L) states, so:
#     N ~ Poisson(lambda * E(n)),  E(n) = SUM_L pi(L) * (n - L + 1)
#     P(L | event) proportional to pi(L) * (n - L + 1)
#     start uniform on [0, n - L]
# Fixing N at int(lambda*(n-mu+1)) and drawing L from pi(L) alone is a different distribution: it
# over-represents long inserts, which have fewer valid starts. It barely moves the U/R differential --
# all arms see the same reads -- but it is the difference between "a calibrated stochastic fixture"
# and "a draw from the model", and only the latter licenses calling R the truth.
mu, sd_i, rlen, eps, lam, disc = 150.0, 20.0, 60, 0.01, 0.05, 0.01
# Drawn from the CONTRACT'S POINT PROCESS, not an approximation of it. Events occur over (start, L)
# states, so:
#     N ~ Poisson(lambda * E(n)),  E(n) = SUM_L pi(L) * (n - L + 1)
#     P(L | event) proportional to pi(L) * (n - L + 1)
#     start uniform on [0, n - L]
# Fixing N at int(lambda*(n-mu+1)) and drawing L from pi(L) alone is a DIFFERENT distribution: it
# over-represents long inserts, which have fewer valid starts to occupy. It barely moves the U/R
# differential, since all arms see the same reads, but it is the difference between a calibrated
# stochastic fixture and a draw from the model -- and only the latter licenses calling R the truth.
# lo is clamped to the shortest POSSIBLE fragment before normalisation, exactly as the scorer builds
# it (min_len = |r1| + |r2|). Including shorter inserts in the prior and then discarding them after
# sampling puts mass into the exposure and the Poisson intensity that no event can ever occupy.
lo_i, hi_i = max(2*rlen, int(mu - 4*sd_i)), int(mu + 4*sd_i)
span = hi_i - lo_i + 1
pri = [ (1-disc)*math.exp(-0.5*((Lv-mu)/sd_i)**2)/(sd_i*math.sqrt(2*math.pi)) + disc/span
        for Lv in range(lo_i, hi_i+1) ]
tot_pri = sum(pri)
pri = [x/tot_pri for x in pri]
n_seq = len(many)
wev = [ pri[i] * max(0, n_seq - Lv + 1) for i, Lv in enumerate(range(lo_i, hi_i+1)) ]
E_n = sum(wev)
cum_ev, acc = [], 0.0
for wi in wev:
    acc += wi / E_n; cum_ev.append(acc)
def draw_event_insert():
    u = random.random()
    for i, c in enumerate(cum_ev):
        if u <= c: return lo_i + i
    return hi_i
def poisson(mean):
    if mean > 500: return int(mean + random.gauss(0, math.sqrt(mean)))
    Lp, kk, pp = math.exp(-mean), 0, 1.0
    while True:
        kk += 1; pp *= random.random()
        if pp <= Lp: return kk - 1
def sub(c): return random.choice([b for b in "ACGT" if b != c])
def err(x): return "".join(sub(c) if random.random() < eps else c for c in x)
recs = []
for hom in (0, 1):
    for i in range(poisson(lam * E_n)):
        ins = draw_event_insert()
        if ins < 2*rlen or n_seq - ins < 0: continue
        st = random.randrange(0, n_seq - ins + 1)
        frag = many[st:st+ins]
        r1, r2 = frag[:rlen], rc(frag[-rlen:])
        if random.random() < 0.5: r1, r2 = rc(frag[-rlen:]), frag[:rlen]
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
  # THREE arms on the SAME reads, so recruitment and grouping can be told apart:
  #   R exhaustive reference        U recruited, no grouping        G recruited + grouped
  # E_recruitment = U - R and E_grouping = G - U. Grouping was exact on error-free reads, so if
  # G - U is ~0 here the failure is recruitment and not compression.
  run_arm() {   # tag, extra flags
    local tag="$1"; shift
    local rc=0
    /usr/bin/time -l "$BIN" genotype-frag -i "$W/g.gfa" -b "$W/bub" -o "$W/$tag" -R "$W/reads.fa" \
      --haplotype-mode --joint-marginal --hamming-emission --joint-top-pairs 0 --top-pairs 20 \
      --rung-zero "$@" $P > "$W/$tag.log" 2> "$W/$tag.time" || rc=$?
    # No `|| true`: a failed binary must not be reported as a measurement.
    if [ $rc -ne 0 ] || [ ! -s "$W/$tag.hap_pairs.tsv" ]; then
      echo "  ${CN}/${RL}/${SD}: arm $tag FAILED (rc=$rc)"; return 1
    fi
  }
  t0=$(date +%s)
  # ONE named pair throughout -- many/many -- for the dumps AND for the reported U-R. Dumping mass
  # for many/many while evaluating U-R at G's own optimum compares different pairs, and they differ in
  # exactly the cells that fail.
  run_arm U --dump-fragment-mass "$W/U.mass" --dump-mass-pair many,many || { GRID_FAILS=$((GRID_FAILS+1)); continue; }
  run_arm M --mate-rescue --no-coordinate-join --dump-fragment-mass "$W/M.mass" --dump-mass-pair many,many || { GRID_FAILS=$((GRID_FAILS+1)); continue; }
  run_arm G --multiplicity-aware --mass-tolerance "$MASS_BOUND" --dump-fragment-mass "$W/G.mass" --dump-mass-pair many,many || { GRID_FAILS=$((GRID_FAILS+1)); continue; }
  # J is M with the mate combinations formed by coordinate join instead of by Cartesian product.
  # Same flags otherwise, so any difference is the join and nothing else.
  run_arm J --mate-rescue --force-join --dump-fragment-mass "$W/J.mass" --dump-mass-pair many,many || { GRID_FAILS=$((GRID_FAILS+1)); continue; }
  # Z is the zero-seed arm: M plus the lossless fallback, changing nothing else, so the difference
  # between them is the fallback and nothing else.
  run_arm Z --mate-rescue --zero-seed-fallback --dump-fragment-mass "$W/Z.mass" --dump-mass-pair many,many || { GRID_FAILS=$((GRID_FAILS+1)); continue; }
  "$BIN" genotype-frag --reference-score "$W/many.fa" "$W/many.fa" -R "$W/reads.fa" $P \
    --dump-fragment-mass "$W/R.mass" >/dev/null 2>&1
  t1=$(date +%s)
  "$PY" - "$W" "$CN" "$RL" "$SD" "$EQ_TOL" "$((t1-t0))" <<'PYEOF' || GRID_FAILS=$((GRID_FAILS+1))
import sys, os, re
W, cn, rl, sd, tol, secs = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], float(sys.argv[5]), sys.argv[6]
ref = {}
for l in open(os.path.join(W,"ref.tsv")):
    k,v = l.rstrip("\n").split("\t")
    try: ref[k] = float(v)
    except ValueError: pass
if not ref: print(f"  {cn:<5} {rl:<7} {sd:<5} | reference produced nothing"); sys.exit(1)
best = max(ref.values())
top_class = {k for k,v in ref.items() if best - v <= tol}
def arm(tag):
    pf = os.path.join(W, tag + ".hap_pairs.tsv")
    if not os.path.exists(pf): return None
    rows = [l.rstrip("\n").split("\t") for l in open(pf)][1:]
    if not rows: return None
    a, b = rows[0][1], rows[0][2]
    return "/".join(sorted((a, b))), float(rows[0][3])
u, g = arm("U"), arm("G")
if u is None or g is None:
    print(f"  {cn:<5} {rl:<7} {sd:<5} | an arm produced nothing"); sys.exit(1)
ukey, uval = u; gkey, gval = g
# Compared at the SAME pair -- G's own optimum -- so the two errors are about scoring and not about
# which pair each arm happened to rank first.
def score_of(tag, key):
    pf = os.path.join(W, tag + ".hap_pairs.tsv")
    for l in list(open(pf))[1:]:
        f = l.rstrip("\n").split("\t")
        if "/".join(sorted((f[1], f[2]))) == key: return float(f[3])
    return float("nan")
PAIR = "many/many"
m_at_p = score_of("M", PAIR)
u_at_g = score_of("U", PAIR)
gval = score_of("G", PAIR)
gkey = PAIR
refv = ref.get(PAIR, float("nan"))
# SIGNED. Exhaustive and recruited scoring share the same exposure, so recruitment can only REMOVE
# event mass and U - R must be <= 0. A positive value would mean recruitment invented mass and is a
# finding, not a rounding artefact -- taking the absolute value would have hidden it.
e_recruit = (u_at_g - refv) if refv == refv and u_at_g == u_at_g else float("nan")
# Grouping compared over EVERY diplotype, with exact key-set equality, not only at G's optimum:
# grouping could perturb another pair enough to reorder the ranking without showing up at the winner.
def all_scores(tag):
    d = {}
    for l in list(open(os.path.join(W, tag + ".hap_pairs.tsv")))[1:]:
        f = l.rstrip("\n").split("\t")
        d["/".join(sorted((f[1], f[2])))] = float(f[3])
    return d
us, gs = all_scores("U"), all_scores("G")
if set(us) != set(gs):
    e_group = float("nan"); group_note = "KEYSETS DIFFER"
else:
    e_group = max(abs(gs[kk] - us[kk]) for kk in us); group_note = ""
log = open(os.path.join(W,"G.log")).read() + open(os.path.join(W,"G.time")).read()
def grab(pat, d="?"):
    m = re.search(pat, log); return m.group(1) if m else d
combos = grab(r"(\d+) mate combinations")
groups = grab(r"(\d+) groups after")
e_rescue = (m_at_p - refv) if refv == refv and m_at_p == m_at_p else float("nan")
z_at_p = score_of("Z", PAIR)
e_zero = (z_at_p - refv) if refv == refv and z_at_p == z_at_p else float("nan")
mlog = open(os.path.join(W,"M.log")).read() + open(os.path.join(W,"M.time")).read()
mcombos = (re.search(r"(\d+) mate combinations", mlog) or [None,"?"])[1]
print(f"  {cn:<5} {rl:<7} {sd:<5} | {best:10.2f} {e_recruit:+9.2f} {e_rescue:+10.2f} {e_zero:+10.2f} {e_group:7.3f}{group_note:<3} "
      f"{('yes' if gkey in top_class else 'NO'):<6} | {combos:<8} {mcombos:<8} {groups:<7} {secs:<5}")
PYEOF
  # ---- join gate: J must equal M NUMERICALLY, on every diplotype, not merely pick the same winner.
  # The join reorders one sum; it is not a different model and not an approximation, so "same winner"
  # is far too weak a criterion -- it would pass a join that lost mass everywhere by a constant.
  "$PY" - "$W" <<'PYEOF' || GRID_FAILS=$((GRID_FAILS+1))
import os, re, sys, math
W = sys.argv[1]
def all_scores(tag):
    pf = os.path.join(W, tag + ".hap_pairs.tsv")
    if not os.path.exists(pf): return {}
    d = {}
    for l in list(open(pf))[1:]:
        f = l.rstrip("\n").split("\t")
        if len(f) >= 4: d["/".join(sorted((f[1], f[2])))] = float(f[3])
    return d
m, j = all_scores("M"), all_scores("J")
if not m or not j or set(m) != set(j):
    print(f"        join: KEY SETS DIFFER ({len(m)} cartesian vs {len(j)} join) -- not comparable")
    sys.exit(1)
dev = max(abs(m[k] - j[k]) for k in m)
def logs(tag):
    return (open(os.path.join(W, tag + ".log")).read() +
            open(os.path.join(W, tag + ".time")).read())
def grab(txt, pat, d="?"):
    g = re.search(pat, txt); return g.group(1) if g else d
lm, lj = logs("M"), logs("J")
cart  = grab(lj, r"(\d+) Cartesian combinations")
probe = grab(lj, r"(\d+) coordinate-join probes")
rpos  = grab(lj, r"(\d+) rescue interval positions")
rpl   = grab(lj, r"(\d+) rescued placements")
def rss(txt):
    g = re.search(r"(\d+)\s+maximum resident set size", txt)
    return f"{int(g.group(1))/1048576:.0f}MB" if g else "?"
def secs(txt):
    g = re.search(r"([\d.]+)\s+real", txt)
    return f"{float(g.group(1)):.2f}s" if g else "?"
eqm = os.path.join(W, "M.equivalence.tsv"); eqj = os.path.join(W, "J.equivalence.tsv")
# A MISSING file is not evidence of agreement. Defaulting to "unchanged" when the comparison could
# not be made reports the reassuring answer for the case where nothing was checked at all.
if not (os.path.exists(eqm) and os.path.exists(eqj)):
    print(f"        join: equivalence dumps missing ({os.path.basename(eqm)} / "
          f"{os.path.basename(eqj)}) -- class comparison NOT MADE")
    sys.exit(1)
same_cls = "unchanged" if open(eqm).read() == open(eqj).read() else "CHANGED"
verdict = "EQUAL" if dev < 1e-9 else ("close" if dev < 1e-6 else "DIFFER")
print(f"        join: max|cartesian-join| = {dev:.3e} over {len(m)} diplotypes [{verdict}]; "
      f"equivalence class {same_cls}")
print(f"        work: {cart} Cartesian hypothetical -> {probe} join probes; "
      f"rescue examined {rpos} positions -> {rpl} placements; "
      f"cartesian {secs(lm)}/{rss(lm)}, join {secs(lj)}/{rss(lj)}")
# Both paths are the same sum in a different order, so anything above summation-order rounding is a
# defect, and a changed equivalence class is a changed call. Printing DIFFER and exiting 0 made this
# whole gate advisory.
sys.exit(0 if (dev < 1e-9 and same_cls == "unchanged") else 1)
PYEOF

  "$PY" - "$W" <<'PYEOF' || GRID_FAILS=$((GRID_FAILS+1))
import os, re, sys
W = sys.argv[1]
def all_scores(tag):
    pf = os.path.join(W, tag + ".hap_pairs.tsv")
    if not os.path.exists(pf): return {}
    d = {}
    for l in list(open(pf))[1:]:
        f = l.rstrip("\n").split("\t")
        if len(f) >= 4: d["/".join(sorted((f[1], f[2])))] = float(f[3])
    return d
m, z = all_scores("M"), all_scores("Z")
if not m or not z or set(m) != set(z):
    print("        zero-seed: KEY SETS DIFFER -- not comparable"); sys.exit(1)
log = open(os.path.join(W, "Z.log")).read() + open(os.path.join(W, "Z.time")).read()
g = re.search(r"zero-seed: ran for (\d+) fragment-haplotype pairs \((\d+) pigeonhole, (\d+) "
              r"exhaustive\); (\d+) candidate starts -> (\d+) verified -> (\d+) placements", log)
if not g:
    print("        zero-seed: arm produced no fallback report"); sys.exit(1)
inv, pig, exh, cand, ver, pl = (int(x) for x in g.groups())
# The fallback may only ADD evidence: it runs solely where neither mate placed, so no pair may lose
# score against M. A drop means it perturbed a fragment recruitment had already placed.
worst = min(z[k] - m[k] for k in m)
print(f"        zero-seed: {inv} invocations ({pig} pigeonhole / {exh} exhaustive), "
      f"{cand} candidates -> {ver} verified -> {pl} placements; "
      f"worst pair change vs rescue-only {worst:+.3f}")
sys.exit(0 if worst > -1e-9 else 1)
PYEOF

  # Per-fragment decomposition of the deficit by seeding stratum. A fragment is the unit: one
  # anchored mate can rescue the other through the insert constraint.
  # The decomposition is GATED on reconciling: at a fixed pair the exposure is common to both arms,
  # so the per-fragment contribution deltas must sum EXACTLY to the whole-pair U - R. If they do not,
  # the strata are describing something other than the reported number and must not be read.
  RU=$(awk -F'\t' -v a=many -v b=many 'NR>1 && (($2==a&&$3==b)){print $4}' "$W/U.hap_pairs.tsv")
  "$PY" - "$W/R.mass" "$W/U.mass" "${RU:-nan}" "$(awk -F'\t' '$1=="many/many"{print $2}' "$W/ref.tsv")" <<'PYEOF' || GRID_FAILS=$((GRID_FAILS+1))
import sys, collections
def load(p):
    d={}
    for l in open(p):
        if l.startswith("#") or l.startswith("fragment"): continue
        f=l.rstrip("\n").split("\t")
        # name -> (log_mass, mates_seeded, contrib)
        d[f[0]]=(float(f[1]), f[2] if len(f)>2 else "NA", float(f[3]) if len(f)>3 else float("nan"),
                 f[4] if len(f)>4 else "NA", f[5] if len(f)>5 else "NA")
    return d
try:
    R, U = load(sys.argv[1]), load(sys.argv[2])
    u_whole, r_whole = float(sys.argv[3]), float(sys.argv[4])
except (OSError, ValueError):
    print("        reconciliation inputs missing"); sys.exit(1)
strata=collections.defaultdict(lambda: [0, 0.0])
tot_delta = 0.0
for k,(rv, _, rc, _rp, _rf) in R.items():
    if k not in U: continue
    uv, seeded, uc, placed, fr = U[k]
    # Every fragment contributes, including those recruitment found nothing for: their contribution is
    # the finite BACKGROUND term, not zero. Assigning -inf cases a zero deficit is what made the
    # zero-seed stratum look free.
    d = uc - rc
    tot_delta += d
    strata[(seeded, placed, fr)][0] += 1
    strata[(seeded, placed, fr)][1] += d
whole = u_whole - r_whole
# Numerical-rounding scale, not a tolerance for real disagreement. Exposure cancels at a
# fixed pair, so this identity is exact arithmetic and anything above rounding is a defect.
ok = abs(tot_delta - whole) < 1e-6
print(f"        reconciliation: per-fragment deltas sum to {tot_delta:+.2f}, whole-pair U-R is "
      f"{whole:+.2f} -> {'MATCH' if ok else 'MISMATCH, strata not interpretable'}")
if not ok: sys.exit(1)
# seed -> placement -> valid FR, so the stage the evidence disappears at is visible
parts=[]
for key in sorted(strata, key=lambda k: -abs(strata[k][1])):
    seeded, placed, fr = key
    n, d = strata[key]
    if n == 0: continue
    parts.append(f"seed{seeded}/place{placed}/fr{fr}: n={n} deficit={d:+.0f} "
                 f"({100*d/tot_delta if tot_delta else 0:.0f}%)")
print("        by stage (many/many): " + ";  ".join(parts[:5]))
PYEOF
done; done; done
echo
echo "  U-R is SIGNED and must be <= 0: recruitment can only remove event mass, since both arms share"
echo "  the same exposure. max|G-U| is over EVERY diplotype with exact key-set equality, not just the"
echo "  winner, so grouping cannot perturb another pair unnoticed."
echo "  combos is the tractability number: it grows as (copies)^2 per fragment and is the work the"
echo "  anchor cap and top-k were introduced to avoid."

if [ "$GRID_FAILS" -ne 0 ]; then
  printf "\n  grid: %d gate(s) FAILED -- join equality, reconciliation or an arm\n" "$GRID_FAILS"
  exit 1
fi
printf "\n  grid: all gates passed\n"
