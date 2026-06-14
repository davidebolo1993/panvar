#!/usr/bin/env bash
# Synthetic smoke: exercises every call event type + confounders on tiny hand-built
# graphs (tests/synthetic_data/), asserting the EXACT expected records. No vg required
# (snarls are checked in). Verifies we catch what we claim, with no double-counting and
# no copy-number event also reported as an INS.
#   syn.gfa   : P-lines; DEL / INS(novel) / INV / self-loop DUP / peak DUP + confounders.
#   syn_w.gfa : W-lines; INV + a multi-copy locus in MIXED orientation (+,-,+).
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <panvar_bin> <synthetic_data_dir> <out_dir>" >&2
  exit 2
fi
PANVAR_BIN="$1"; DATA_DIR="$2"; OUT_DIR="$3"
mkdir -p "$OUT_DIR"

run_case() {  # <gfa> <snarls> <out_prefix>
  local gfa="$1" snarls="$2" pref="$3"
  [[ -f "$gfa" && -f "$snarls" ]] || { echo "error: missing $gfa / $snarls" >&2; exit 1; }
  "$PANVAR_BIN" bubble -i "$gfa" -o "${pref}.bub" --snarls-in "$snarls" --quiet
  # 'synref' is a unique substring of the reference path (also tests --reference-path matching).
  "$PANVAR_BIN" call -i "$gfa" --bubble-prefix-in "${pref}.bub" \
    --reference-path synref -o "${pref}.call" --cn-from-multiplicity --classify-ins --quiet
}

run_case "$DATA_DIR/syn.gfa"   "$DATA_DIR/syn.snarls.jsonl"   "$OUT_DIR/p"
run_case "$DATA_DIR/syn_w.gfa" "$DATA_DIR/syn_w.snarls.jsonl" "$OUT_DIR/w"

python3 - "$OUT_DIR/p.call.region.vcf" "$OUT_DIR/w.call.region.vcf" <<'PY'
import sys

def load(path):
    hdr, recs = None, []
    for line in open(path):
        if line.startswith("#CHROM"): hdr = line.rstrip().split("\t")
        if line.startswith("#"): continue
        recs.append(line.rstrip().split("\t"))
    return hdr, recs

def info(f): return {k.split('=')[0]: (k.split('=')[1] if '=' in k else '') for k in f[7].split(';')}
def carriers(hdr, f):
    return sorted(s.split('#')[0] for s, v in zip(hdr[9:], f[9:])
                  if v.split(":")[0] not in ("0/0", "./.", "0", ".", "0|0"))

fails = []
def check(cond, msg):
    print(("ok  " if cond else "FAIL") + " " + msg)
    if not cond: fails.append(msg)

# ---- P-line graph: the seven isolated events + confounders --------------------
hdr, recs = load(sys.argv[1])
by = {}
for f in recs:
    by.setdefault((info(f).get("SVTYPE"), frozenset(carriers(hdr, f))), []).append(f)
def find(sv, c): return by.get((sv, frozenset(c)), [])

check(len(find("DEL", ["s_del"])) == 1, "[P] DEL of DEL1 carried by s_del")
ins = find("INS", ["s_ins"])
check(len(ins) == 1 and info(ins[0]).get("INS_SUBTYPE") == "NOVEL",
      "[P] INS of INS1 carried by s_ins, INS_SUBTYPE=NOVEL")
check(len(find("INV", ["s_inv"])) == 1, "[P] INV of INV1 carried by s_inv")
sd = find("DUP", ["s_selfdup"])
check(len(sd) == 1 and info(sd[0]).get("REF_CN") == "1",
      "[P] self-loop DUP carried by s_selfdup, REF_CN=1")
pk = find("DUP", ["s_peakdup1", "s_peakdup2"])
check(len(pk) == 1 and info(pk[0]).get("REF_CN") == "2",
      "[P] peak-multiplicity DUP carried by s_peakdup1+s_peakdup2, REF_CN=2")
check(len(find("DEL", ["s_null"])) == 1, "[P] G x1 in s_null is a DEL")
check(len([f for f in recs if info(f).get("SVTYPE") == "DUP"]) == 2,
      "[P] exactly two DUP records; invariant fold B yields no DUP")
check(not any(info(f).get("SVTYPE") == "DUP" and "s_null" in carriers(hdr, f) for f in recs),
      "[P] s_null is never reported as a DUP")
for hap in ["s_selfdup", "s_peakdup1", "s_peakdup2"]:
    n = sum(1 for f in recs if hap in carriers(hdr, f))
    check(n == 1, f"[P] {hap} appears in exactly one record (no double-count), got {n}")

# ---- W-line graph: parsing + mixed-orientation multi-copy ---------------------
hdr, recs = load(sys.argv[2])
inv = [f for f in recs if info(f).get("SVTYPE") == "INV" and "s_winv" in carriers(hdr, f)]
check(len(inv) == 1, "[W] INV called from W-line walk (s_winv)")
mix = [f for f in recs if info(f).get("SVTYPE") == "DUP" and "s_wmixdup" in carriers(hdr, f)]
check(len(mix) == 1 and info(mix[0]).get("REF_CN") == "2",
      "[W] mixed-orientation (+,-,+) multi-copy is a single DUP, REF_CN=2")
n = sum(1 for f in recs if "s_wmixdup" in carriers(hdr, f))
check(n == 1, f"[W] s_wmixdup in exactly one record (mixed orientation not double-counted / not an INV), got {n}")

print()
if fails:
    print(f"SYNTHETIC SMOKE FAILED ({len(fails)} assertion(s))")
    sys.exit(1)
print("synthetic smoke: OK")
PY
