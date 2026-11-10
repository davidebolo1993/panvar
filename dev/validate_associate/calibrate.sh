#!/usr/bin/env bash
# Null calibration of what no tool computes, on validate.sh's simulated genotypes (run validate.sh first):
#   permutations  200 noise phenotypes x 200 shuffles: runs with any row past the permutation threshold (5%)
#   logistic      200 noise case/control phenotypes at 30% and 3% cases: rows below 0.05 and 0.01; exact
#                 counted tails apart, normal and saddlepoint pooled (a row's own |z| picks between them)
#   groups        300 noise case/control phenotypes at 2% cases: burden, SKAT, SKAT-O below 0.05 and 0.01
#   delta_r2      100 noise phenotypes, 300 noise records against one: intervals excluding 0
# Each rate has its Wilson 95% interval and is flagged when all of it lies above the target; a flag the docs list
# as a limitation is "documented" and does not fail. Refused replicates are skipped and counted. The summary is
# also <out dir>/calibration.txt; CALIBRATE_SUMMARY_ONLY=1 reprints it.
#   dev/validate_associate/calibrate.sh <panvar> <validate.sh out dir>
set -euo pipefail
PANVAR=$1 OUT=$2
U=$OUT/unrelated C=$OUT/calibration
[ -s "$U/g.bimbam" ] || { echo "no $U/g.bimbam: run validate.sh first"; exit 1; }
mkdir -p "$C"
python3 - "$U" "$C" <<'PY'
import random, sys
U, C = sys.argv[1], sys.argv[2]
rnd = random.Random(20261107)
rows = [l.rstrip("\n").split("\t") for l in open(f"{U}/pheno.q.tsv")]
head, rows = rows[0], rows[1:]
def write(path, ys):
    with open(path, "w") as f:
        f.write("\t".join(head) + "\n")
        for r, y in zip(rows, ys):
            f.write("\t".join([r[0], str(y)] + r[2:]) + "\n")
for k in range(200):
    write(f"{C}/q{k}.tsv", [round(rnd.gauss(0, 1), 6) for _ in rows])
    for frac, tag in ((0.30, "b30"), (0.03, "b03")):
        write(f"{C}/{tag}_{k}.tsv", [int(rnd.random() < frac) for _ in rows])
for k in range(300):
    write(f"{C}/b02_{k}.tsv", [int(rnd.random() < 0.02) for _ in rows])
# delta_r2: 300 noise records in layer A, one in layer B, 40 samples
for k in range(100):
    n = 40
    with open(f"{C}/d{k}.bimbam", "w") as g, open(f"{C}/d{k}.annot", "w") as a:
        a.write("feature_id\tlayer\tencoding\tbubbles\tnodes\tsvtype\tgene\tAF\tAN\n")
        for j in range(301):
            v = [sum(rnd.random() < 0.3 for _ in range(2)) for _ in range(n)]
            v[0] = 1 if len(set(v)) < 2 else v[0]
            fid = f"A{j}" if j < 300 else "B0"
            g.write(fid + ", X, Y, " + ", ".join(map(str, v)) + "\n")
            a.write(f"{fid}\t{'A' if j < 300 else 'B'}\tdosage\t.\t.\tSNV\t.\t.\t.\n")
    with open(f"{C}/d{k}.pheno", "w") as f:
        f.write("sample\tphenotype\n" + "".join(f"s{i}\t{rnd.gauss(0, 1):.6f}\n" for i in range(n)))
with open(f"{C}/d.samples", "w") as f:
    f.write("".join(f"s{i}\n" for i in range(40)))
PY
if [ "${CALIBRATE_SUMMARY_ONLY:-0}" != 1 ]; then
base=(--genotypes "$U/g.bimbam" --samples "$U/samples.txt" --feature-annot "$U/annot.tsv" --unit variant --min-maf 0 --finemap-l 0 -q)
: > "$C/perm.txt"; : > "$C/logistic.txt"; : > "$C/groups.txt"; : > "$C/delta.txt"; : > "$C/refused.txt"
pv() { "$PANVAR" associate "$@" 2>> "$C/refused.log" || { echo "$*" >> "$C/refused.txt"; return 1; }; }
for k in $(seq 0 199); do
  pv "${base[@]}" --phenotype "$C/q$k.tsv" --permutations 200 --seed $((k + 1)) -o "$C/r" &&
  awk -F'\t' '$1 == "significant_permutation" { print ($2 > 0) }' "$C/r.summary.tsv" >> "$C/perm.txt"
  for t in b30 b03; do
    pv "${base[@]}" --phenotype "$C/${t}_$k.tsv" -o "$C/r" || continue
    awk -F'\t' -v t=$t 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next } $(c["p"]) != "NA" { print t "\t" $(c["p_method"]) "\t" $(c["p"]) }' "$C/r.assoc.tsv" >> "$C/logistic.txt"
  done
done
for k in $(seq 0 299); do
  pv "${base[@]}" --phenotype "$C/b02_$k.tsv" --group gene --rare-maf 0.025 -o "$C/r" || continue
  awk -F'\t' 'NR > 1 { print $4 "\t" $5 "\t" $6 }' "$C/r.groups.tsv" >> "$C/groups.txt"
done
for k in $(seq 0 99); do
  pv --genotypes "$C/d$k.bimbam" --samples "$C/d.samples" --feature-annot "$C/d$k.annot" --phenotype "$C/d$k.pheno" \
    --unit variant --min-maf 0 --finemap-l 0 -q -o "$C/r" || continue
  awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) c[$i] = i; next } NR == 3 { print ($(c["delta_r2_low"]) > 0 || $(c["delta_r2_high"]) < 0) }' "$C/r.layers.tsv" >> "$C/delta.txt"
done
fi
python3 - "$C" <<'PY' | tee "$OUT/calibration.txt"
import math, sys
C = sys.argv[1]
def ci(k, n):   # Wilson 95%
    if n == 0: return (float("nan"),) * 2
    p, z = k / n, 1.96
    c = (p + z * z / (2 * n)) / (1 + z * z / n); h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h
flags = 0
def line(label, k, n, target, documented=False):
    global flags
    lo, hi = ci(k, n)
    bad = n > 0 and lo > target
    flags += bad and not documented
    mark = ("   << above (documented limitation)" if documented else "   << above") if bad else ""
    print(f"  {label:56s} {k:6d}/{n:<7d} {100 * k / max(n, 1):6.2f}%   95% CI {100 * lo:5.2f}-{100 * hi:5.2f}%   target {100 * target:g}%{mark}")
print("null calibration (a rate is flagged when its whole 95% interval lies above the target)")
print(f"  replicates refused by associate (covariates separating the cases): {sum(1 for _ in open(f'{C}/refused.txt'))}")
v = [int(x) for x in open(f"{C}/perm.txt")]
line("permutation threshold: runs with a false hit", sum(v), len(v), 0.05)
rows = [l.split("\t") for l in open(f"{C}/logistic.txt")]
for t, lab in (("b30", "30% cases"), ("b03", "3% cases")):
    for kind, ms in (("normal or saddlepoint tail", ("score", "score_spa")), ("counted tail", ("score_exact_carriers",))):
        ps = [float(r[2]) for r in rows if r[0] == t and r[1] in ms]
        if not ps: continue
        for a in (0.05, 0.01):
            line(f"logistic, {lab}, {kind}: p < {a}", sum(p < a for p in ps), len(ps), a)
g = [[float(x) for x in l.split("\t")] for l in open(f"{C}/groups.txt")]
for j, name in enumerate(("burden", "SKAT", "SKAT-O")):
    ps = [r[j] for r in g if math.isfinite(r[j])]
    for a in (0.05, 0.01):
        line(f"groups, 2% cases, {name}: p < {a}", sum(p < a for p in ps), len(ps), a, documented=name != "burden")
v = [int(x) for x in open(f"{C}/delta.txt")]
line("delta_r2 on noise: intervals excluding 0", sum(v), len(v), 0.05)
sys.exit(1 if flags else 0)
PY
exit "${PIPESTATUS[0]}"
