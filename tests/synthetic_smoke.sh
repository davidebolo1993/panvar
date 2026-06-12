#!/usr/bin/env bash
# Synthetic smoke: exercises every call event type + confounders on a tiny hand-built
# graph (tests/synthetic_data/), asserting the EXACT expected records. No vg required
# (snarls are checked in). Verifies we catch what we claim, with no double-counting and
# no copy-number event also reported as an INS.
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <panvar_bin> <synthetic_data_dir> <out_dir>" >&2
  exit 2
fi
PANVAR_BIN="$1"; DATA_DIR="$2"; OUT_DIR="$3"
GFA="$DATA_DIR/syn.gfa"
SNARLS="$DATA_DIR/syn.snarls.jsonl"

for f in "$GFA" "$SNARLS"; do
  [[ -f "$f" ]] || { echo "error: missing $f" >&2; exit 1; }
done
mkdir -p "$OUT_DIR"

"$PANVAR_BIN" bubble -i "$GFA" -o "$OUT_DIR/bub" --snarls-in "$SNARLS" --quiet
# 'synref' is a unique substring of the reference path name (tests --reference-path matching).
"$PANVAR_BIN" call -i "$GFA" --bubble-prefix-in "$OUT_DIR/bub" \
  --reference-path synref -o "$OUT_DIR/call" --cn-from-multiplicity --classify-ins --quiet

python3 - "$OUT_DIR/call.region.vcf" <<'PY'
import sys
vcf = sys.argv[1]
hdr = None
recs = []
for line in open(vcf):
    if line.startswith("#CHROM"): hdr = line.rstrip().split("\t")
    if line.startswith("#"): continue
    recs.append(line.rstrip().split("\t"))

def info(f): return {k.split('=')[0]: (k.split('=')[1] if '=' in k else '') for k in f[7].split(';')}
def carriers(f):
    return sorted(s.split('#')[0] for s, v in zip(hdr[9:], f[9:])
                  if v.split(":")[0] not in ("0/0", "./.", "0", ".", "0|0"))

fails = []
def check(cond, msg):
    print(("ok  " if cond else "FAIL") + " " + msg)
    if not cond: fails.append(msg)

# index records by (svtype, frozenset(carriers))
by = {}
for f in recs:
    by.setdefault((info(f).get("SVTYPE"), frozenset(carriers(f))), []).append(f)

def find(svtype, carrs):
    return by.get((svtype, frozenset(carrs)), [])

# --- the seven expected events, each by type + exact carrier set ---
check(len(find("DEL", ["s_del"])) == 1, "DEL of DEL1 carried by s_del")
ins = find("INS", ["s_ins"])
check(len(ins) == 1 and info(ins[0]).get("INS_SUBTYPE") == "NOVEL",
      "INS of INS1 carried by s_ins, INS_SUBTYPE=NOVEL")
check(len(find("INV", ["s_inv"])) == 1, "INV of INV1 carried by s_inv")
selfdup = find("DUP", ["s_selfdup"])
check(len(selfdup) == 1 and info(selfdup[0]).get("REF_CN") == "1",
      "self-loop DUP carried by s_selfdup with REF_CN=1")
peak = find("DUP", ["s_peakdup1", "s_peakdup2"])
check(len(peak) == 1 and info(peak[0]).get("REF_CN") == "2",
      "peak-multiplicity DUP carried by s_peakdup1+s_peakdup2 with REF_CN=2")
check(len(find("DEL", ["s_null"])) == 1, "G x1 in s_null is a DEL")

# --- confounders / no double-counting ---
check(len([f for f in recs if info(f).get("SVTYPE") == "DUP"]) == 2,
      "exactly two DUP records (self-loop + peak); invariant fold B yields no DUP")
# s_null must NOT appear in any DUP (low multiplicity is not a duplication)
check(not any(info(f).get("SVTYPE") == "DUP" and "s_null" in carriers(f) for f in recs),
      "s_null is never reported as a DUP")
# each duplication carrier appears in EXACTLY ONE record (no DUP-also-as-INS double count)
for hap in ["s_selfdup", "s_peakdup1", "s_peakdup2"]:
    n = sum(1 for f in recs if hap in carriers(f))
    check(n == 1, f"{hap} appears in exactly one record (no double-count), got {n}")

print()
if fails:
    print(f"SYNTHETIC SMOKE FAILED ({len(fails)} assertion(s))")
    sys.exit(1)
print("synthetic smoke: OK (%d records)" % len(recs))
PY
