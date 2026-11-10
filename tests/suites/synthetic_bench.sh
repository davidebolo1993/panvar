#!/usr/bin/env bash
# synthetic_bench.sh - end-to-end check on a generated locus whose answers are known by construction.
#
# tests/tools/make_synthetic_locus.py plants four SVs (DEL, INS, a DUP with 1-3 copies of a unit, INV)
# in a backbone and writes the graph, the designs and where each SV sits. The pipeline bubble ->
# panphorte -> call --cn must find exactly one site per SV, none nested, and call each SV once, at its
# planted position, carried by exactly the haplotypes whose design has it, with every haplotype's DUP
# copy number equal to its design's. Real loci cannot check this: their truth is whatever the graph
# encodes.
#
#   synthetic_bench.sh <panvar_bin> <out_dir>
#
# Env: PYTHON (default python3), SEED (default 7), GEN_EXTRA (extra flags for make_synthetic_locus.py).
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <panvar_bin> <out_dir>" >&2
  exit 2
fi
BIN="$1"; OUT="$2"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PY="${PYTHON:-python3}"
SEED="${SEED:-7}"

mkdir -p "$OUT"
# shellcheck disable=SC2086  # GEN_EXTRA is a flag list
"$PY" "$HERE/tools/make_synthetic_locus.py" -o "$OUT" --seed "$SEED" ${GEN_EXTRA:-} >/dev/null

"$BIN" bubble -i "$OUT/graph.gfa" -r ref -o "$OUT/bub" --quiet
"$BIN" panphorte -i "$OUT/bub.sorted.gfa" -b "$OUT/bub" -r ref -o "$OUT/pp" --threads 2 --quiet
# call's warnings are kept out of the test log, but shown when it fails.
"$BIN" call -i "$OUT/pp.normalized.sorted.gfa" -b "$OUT/pp" -r ref -o "$OUT/call" --cn --threads 2 \
  --quiet 2> "$OUT/call.stderr" || { cat "$OUT/call.stderr" >&2; echo "FAIL [call] exited non-zero" >&2; exit 1; }

"$PY" - "$OUT" <<'PY'
import csv, sys
out = sys.argv[1]
fails = []
def check(cond, msg):
    print(("ok  " if cond else "FAIL") + " " + msg)
    if not cond: fails.append(msg)

sites = {r["site"]: r for r in csv.DictReader(open(f"{out}/truth.bubbles.tsv"), delimiter="\t")}
designs = {r["design"]: r for r in csv.DictReader(open(f"{out}/truth.variants.tsv"), delimiter="\t")}
haps = {r["haplotype"]: r["design"] for r in csv.DictReader(open(f"{out}/truth.haplotypes.tsv"), delimiter="\t")}

# bubble: one site per planted SV, none inside another
got = list(csv.DictReader(open(f"{out}/bub.bubbles.csv")))
check(len(got) == len(sites), f"[bubble] {len(got)} sites for {len(sites)} planted SVs")
inside = {g["bubble_id"]: set(g["inside_nodes"].split(";")) for g in got}
nested = [a for a in inside for b in inside if a != b and inside[a] and inside[a] <= inside[b]]
check(not nested, f"[bubble] no site nested in another (nested: {nested})")

# call: each SV once, at its planted position, with the right carriers and copy numbers
hdr, recs = None, []
for line in open(f"{out}/call.region.vcf"):
    if line.startswith("##"): continue
    f = line.rstrip("\n").split("\t")
    if line.startswith("#"): hdr = f; continue
    recs.append(f)
def info(f): return dict(k.split("=", 1) if "=" in k else (k, "") for k in f[7].split(";"))
def fmt(f, i):
    return dict(zip(f[8].split(":"), f[9 + i].split(":")))
samples = hdr[9:]

def carriers(f):
    return {s for i, s in enumerate(samples) if fmt(f, i).get("GT") == "1"}
def expected(field, test):
    return {h for h, d in haps.items() if test(designs[d][field])}

# (site, VCF SVTYPE, planted POS, planted END, field in truth.variants.tsv, carrier test)
pos = {k: int(v["ref_pos"]) for k, v in sites.items()}
size = {k: int(v["size"]) for k, v in sites.items()}
want = [
    ("DEL", "DEL", pos["DEL"], pos["DEL"] + size["DEL"], "DEL", lambda v: v == "1"),
    # the payload is inserted after a reference node as long as itself
    ("INS", "INS", pos["INS"] + size["INS"], pos["INS"] + size["INS"], "INS", lambda v: v == "1"),
    ("DUP", "CNV", pos["DUP"], pos["DUP"] + size["DUP"], "DUP_copies", lambda v: v != "1"),
    ("INV", "INV", pos["INV"], pos["INV"] + size["INV"], "INV", lambda v: v == "1"),
]
for site, svtype, p, e, field, test in want:
    rs = [f for f in recs if info(f).get("SVTYPE") == svtype]
    check(len(rs) == 1, f"[call] {site}: one {svtype} record (got {len(rs)})")
    if len(rs) != 1: continue
    f = rs[0]; i = info(f)
    # A left-aligned event may sit up to CIPOS bases left of where it was planted.
    slack = int(i.get("CIPOS", "0,0").split(",")[1])
    shift = p - int(f[1])
    check(0 <= shift <= slack and int(i["END"]) == e - shift,
          f"[call] {site}: POS {f[1]} END {i['END']} (planted {p}-{e}, slack {slack})")
    check(carriers(f) == expected(field, test), f"[call] {site}: carriers are the designs that have it")
    if site == "DUP":
        cn = {s: fmt(f, k).get("CN") for k, s in enumerate(samples) if s in haps}
        bad = [s for s, c in cn.items() if c != designs[haps[s]]["DUP_copies"]]
        check(not bad, f"[call] DUP: every haplotype's CN is its design's copy count (wrong: {bad})")

print()
if fails:
    print(f"SYNTHETIC BENCH FAILED ({len(fails)} assertion(s))")
    sys.exit(1)
print("synthetic bench: OK")
PY
