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

# ---- flip check: a graph whose reference traverses a node in reverse -----------
# The bundled real graphs are already reference-forward, so flip is a no-op there.
# This tiny graph forces the reverse-complement path: the reference walks node 2 as
# '2-', so the internal sort+flip must reverse-complement node 2, rewrite the L-lines
# and path strands, and leave every spelled sequence unchanged.
FLIP_IN="$OUT_DIR/flip_in.gfa"
FLIP_OUT="$OUT_DIR/flip_out.sorted.gfa"
{
  printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tACGT\n'
  printf 'S\t2\tGGGG\n'
  printf 'S\t3\tTTAACC\n'
  printf 'L\t1\t+\t2\t-\t0M\n'
  printf 'L\t2\t-\t3\t+\t0M\n'
  printf 'P\tsynref\t1+,2-,3+\t*,*\n'
  printf 'P\talt\t1+,2-,3+\t*,*\n'
} > "$FLIP_IN"
"$PANVAR_BIN" bubble -i "$FLIP_IN" --reference-path synref \
  -o "$OUT_DIR/flip" --sorted-gfa-out "$FLIP_OUT" --quiet >/dev/null

python3 - "$FLIP_OUT" <<'PY'
import sys
def rc(s): return s.translate(str.maketrans("ACGT","TGCA"))[::-1]
nodes={}; ref=None
for line in open(sys.argv[1]):
    f=line.rstrip("\n").split("\t")
    if f[0]=="S": nodes[f[1]]=f[2]
    elif f[0]=="P" and f[1]=="synref":
        ref=[(s[:-1],s[-1]) for s in f[2].split(",")]
fails=[]
def check(c,m):
    print(("ok  " if c else "FAIL")+" "+m)
    if not c: fails.append(m)
check(ref is not None, "[flip] reference path 'synref' present in output")
check(all(o=="+" for _,o in ref), "[flip] reference is all-forward after flip")
spelled="".join(nodes[n] if o=="+" else rc(nodes[n]) for n,o in ref)
check(spelled=="ACGTCCCCTTAACC", f"[flip] reference spelling preserved (got {spelled})")
check("CCCC" in nodes.values() and "GGGG" not in nodes.values(),
      "[flip] node walked as '2-' was reverse-complemented (GGGG -> CCCC)")
if fails:
    print(f"\nFLIP CHECK FAILED ({len(fails)} assertion(s))"); sys.exit(1)
print("flip check: OK")
PY

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
