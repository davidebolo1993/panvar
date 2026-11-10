#!/usr/bin/env bash
# Synthetic smoke: exercises every call event type + confounders on tiny hand-built
# graphs (tests/synthetic_data/), asserting the exact expected records and their positions.
# Verifies we catch what we claim, with no double-counting and no copy-number event also
# reported as an INS.
#   syn.gfa   : P-lines; DEL / INS(novel) / INV / self-loop DUP / peak DUP + confounders.
#   syn_w.gfa : W-lines; INV + a multi-copy locus in MIXED orientation (+,-,+).
# The main runs use bubble's own site finder along the reference; the checked-in vg snarls are
# read once and must give the same records. A run without --cn checks that the folded module
# is then reported as plain events. The last checks cover the shared input/output rules.
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <panvar_bin> <synthetic_data_dir> <out_dir>" >&2
  exit 2
fi
PANVAR_BIN="$1"; DATA_DIR="$2"; OUT_DIR="$3"
mkdir -p "$OUT_DIR"

# The events sit a few dozen bp apart, under bubble's default merge distance, so --merge-nearby-bp 0
# keeps each at its own site.
run_case() {  # <gfa> <out_prefix>
  local gfa="$1" pref="$2"
  [[ -f "$gfa" ]] || { echo "error: missing $gfa" >&2; exit 1; }
  # 'synref' is a unique substring of the reference path (also tests --reference-path matching).
  "$PANVAR_BIN" bubble -i "$gfa" -r synref --merge-nearby-bp 0 -o "${pref}.bub" --quiet
  "$PANVAR_BIN" call -i "${pref}.bub.sorted.gfa" --bubble-prefix "${pref}.bub" \
    --reference-path synref -o "${pref}.call" --cn --classify-ins --quiet
}

run_case "$DATA_DIR/syn.gfa"   "$OUT_DIR/p"
run_case "$DATA_DIR/syn_w.gfa" "$OUT_DIR/w"

# The same P-line graph with its sites from the checked-in vg snarls (the graph is not re-sorted then).
"$PANVAR_BIN" bubble -i "$DATA_DIR/syn.gfa" --snarls "$DATA_DIR/syn.snarls.jsonl" -r synref \
  --merge-nearby-bp 0 -o "$OUT_DIR/ps.bub" --quiet
"$PANVAR_BIN" call -i "$DATA_DIR/syn.gfa" --bubble-prefix "$OUT_DIR/ps.bub" \
  --reference-path synref -o "$OUT_DIR/ps.call" --cn --classify-ins --quiet

# Without --cn: the folded G module is no copy-number record, its gains and loss are plain events.
"$PANVAR_BIN" call -i "$OUT_DIR/p.bub.sorted.gfa" --bubble-prefix "$OUT_DIR/p.bub" \
  --reference-path synref -o "$OUT_DIR/p.nocn" --classify-ins --quiet

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

# ---- shared input/output rules -----------------------------------------------------
CLI_FAILS=0
cli_check() {  # <condition status> <message>
  if [[ "$1" == 0 ]]; then echo "ok   [cli] $2"; else echo "FAIL [cli] $2"; CLI_FAILS=$((CLI_FAILS + 1)); fi
}
# A run that fails before writing anything leaves no directory it made behind.
rm -rf "$OUT_DIR/cli"
rc=0; "$PANVAR_BIN" bubble -i "$OUT_DIR/missing.gfa" -r synref -o "$OUT_DIR/cli/a/b" --quiet 2>/dev/null || rc=$?
ok=0; [[ $rc -ne 0 && ! -e "$OUT_DIR/cli" ]] || ok=1
cli_check $ok "a failed bubble run leaves no output directory (rc=$rc)"
# The same for a refine run that fails before staging, on an unreadable sites table.
printf 'not a sites table\n' > "$OUT_DIR/bad.bubbles.csv"
rm -rf "$OUT_DIR/cli"
rc=0; "$PANVAR_BIN" refine -i "$OUT_DIR/p.bub.sorted.gfa" --bubbles-csv "$OUT_DIR/bad.bubbles.csv" \
  --reference-path synref -o "$OUT_DIR/cli/r/b" --quiet 2>/dev/null || rc=$?
ok=0; [[ $rc -ne 0 && ! -e "$OUT_DIR/cli" ]] || ok=1
cli_check $ok "a refine run refused before staging leaves no output directory (rc=$rc)"
# A refine run that fails after its outputs were written (at install) removes the files and the
# directories it made.
rm -rf "$OUT_DIR/cli"
rc=0; PANVAR_TEST_FAIL_COMMIT_AT=1 "$PANVAR_BIN" refine -i "$OUT_DIR/p.bub.sorted.gfa" \
  --bubble-prefix "$OUT_DIR/p.bub" --reference-path synref -o "$OUT_DIR/cli/r/b" --quiet 2>/dev/null || rc=$?
ok=0; [[ $rc -ne 0 && ! -e "$OUT_DIR/cli" ]] || ok=1
cli_check $ok "a refine run that fails after writing its outputs leaves no output directory (rc=$rc)"
# A directory that existed before the run is the user's: a failed run keeps it, even when empty.
rm -rf "$OUT_DIR/cli_pre"; mkdir -p "$OUT_DIR/cli_pre/empty"
rc=0; "$PANVAR_BIN" bubble -i "$OUT_DIR/missing.gfa" -r synref -o "$OUT_DIR/cli_pre/empty/b" --quiet 2>/dev/null || rc=$?
rc2=0; PANVAR_TEST_FAIL_COMMIT_AT=1 "$PANVAR_BIN" refine -i "$OUT_DIR/p.bub.sorted.gfa" \
  --bubble-prefix "$OUT_DIR/p.bub" --reference-path synref -o "$OUT_DIR/cli_pre/empty/r" --quiet 2>/dev/null || rc2=$?
ok=0; [[ $rc -ne 0 && $rc2 -ne 0 && -d "$OUT_DIR/cli_pre/empty" ]] || ok=1
[[ $ok -ne 0 ]] || [[ -z "$(ls -A "$OUT_DIR/cli_pre/empty")" ]] || ok=1
cli_check $ok "failed runs keep a pre-existing empty output directory, and leave nothing in it (rc=$rc,$rc2)"
# Two runs share a new directory: A makes it, B stages into it, then A fails before B writes anything.
# B must still succeed. Each run reads its graph from a FIFO, which it opens only after staging, so a
# writer's open completing tells the order.
rm -rf "$OUT_DIR/cli_race"
ok=0; python3 - "$PANVAR_BIN" "$DATA_DIR/syn.gfa" "$OUT_DIR/cli_race" <<'PY' || ok=1
import errno, os, subprocess, sys, time
bin_, gfa, d = sys.argv[1:4]
os.makedirs(d)
procs = {}
def start(name):
    fifo = os.path.join(d, name + ".fifo")
    os.mkfifo(fifo)
    procs[name] = subprocess.Popen([bin_, "bubble", "-i", fifo, "-r", "synref", "-o",
                                    os.path.join(d, "shared", name), "--quiet"],
                                   stdout=subprocess.DEVNULL, stderr=open(os.path.join(d, name + ".err"), "w"))
    deadline = time.time() + 60
    while True:   # a non-blocking open for writing succeeds once the reader has opened the FIFO
        try:
            fd = os.open(fifo, os.O_WRONLY | os.O_NONBLOCK)
            os.set_blocking(fd, True)
            return fd
        except OSError as e:
            if e.errno != errno.ENXIO or procs[name].poll() is not None or time.time() > deadline:
                for p in procs.values(): p.kill()
                sys.exit(f"run {name} did not open its input (rc={procs[name].poll()})")
            time.sleep(0.05)
fa = start("a")
if not os.path.isdir(os.path.join(d, "shared")):
    for p in procs.values(): p.kill()
    sys.exit("bubble did not make its output directory before reading its input; this check needs that order")
fb = start("b")
try:
    os.write(fa, b"S\t1\n"); os.close(fa)
    ra = procs["a"].wait(timeout=120)
    with open(gfa, "rb") as f: data = f.read()
    os.write(fb, data); os.close(fb)
    rb = procs["b"].wait(timeout=120)
except (OSError, subprocess.TimeoutExpired) as e:
    # Say which run was still going and where it waited, before killing it.
    for name, p in procs.items():
        state = "running" if p.poll() is None else f"rc={p.returncode}"
        print(f"run {name}: {state}; stderr: {open(os.path.join(d, name + '.err')).read().strip()!r}")
        if p.poll() is None and sys.platform == "darwin":
            subprocess.run(["sample", str(p.pid), "1"], stdout=sys.stdout, stderr=subprocess.DEVNULL)
    for p in procs.values(): p.kill()
    sys.exit(f"race runs did not finish: {e}")
ok = ra != 0 and rb == 0 and os.path.isfile(os.path.join(d, "shared", "b.bubbles.csv"))
print(f"run a rc={ra}, run b rc={rb}")
sys.exit(0 if ok else 1)
PY
cli_check $ok "a failed run keeps a new output directory a concurrent run is using"
# A run killed with SIGKILL cannot remove its hold file; the next run in that directory removes it, and
# any hold whose process is gone and whose lock is free, but keeps one of a live process or a held lock
# (a run on another host sharing the directory). A hold left while it was being made (-holding) goes only
# once it is a minute old. The lock check below shows a live run holds its lock; it cannot catch the
# moment between making the file and locking it.
rm -rf "$OUT_DIR/cli_stale"
ok=0; python3 - "$PANVAR_BIN" "$DATA_DIR/syn.gfa" "$OUT_DIR/cli_stale" <<'PY' || ok=1
import errno, fcntl, os, re, signal, subprocess, sys, time
bin_, gfa, d = sys.argv[1:4]
k = os.path.join(d, "k")
os.makedirs(d)
fifo = os.path.join(d, "in.fifo")
os.mkfifo(fifo)
p = subprocess.Popen([bin_, "bubble", "-i", fifo, "-r", "synref", "-o", os.path.join(k, "x"), "--quiet"],
                     stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
deadline = time.time() + 60
while True:   # the run has staged once it opens its input; it then waits on it until killed
    try:
        writer = os.open(fifo, os.O_WRONLY | os.O_NONBLOCK); break
    except OSError as e:
        if e.errno != errno.ENXIO or p.poll() is not None or time.time() > deadline:
            p.kill(); sys.exit(f"bubble did not open its input (rc={p.poll()})")
        time.sleep(0.05)
hold_re = re.compile(r"-hold\.([0-9]+)-tmp\.")
killed = [f for f in os.listdir(k) if hold_re.search(f)]
def locked(path):   # another process holds the lock
    fd = os.open(path, os.O_RDWR)
    try:
        fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB); return False
    except BlockingIOError:
        return True
    finally:
        os.close(fd)
live_locked = bool(killed) and all(locked(os.path.join(k, f)) for f in killed)
os.kill(p.pid, signal.SIGKILL); p.wait(); os.close(writer)
if not live_locked:
    sys.exit(f"the live run's hold is not locked: {killed}")
if not killed or any(hold_re.search(f).group(1) != str(p.pid) for f in killed):
    sys.exit(f"the run's hold does not name its process {p.pid}: {sorted(os.listdir(k))}")
dead = subprocess.Popen(["true"]); dead.wait()
fakes = {"live": f"y.bubble-hold.{os.getpid()}-tmp.1.0", "locked": f"z.bubble-hold.{dead.pid}-tmp.1.0",
         "stale": f"w.bubble-hold.{dead.pid}-tmp.1.0", "other": "v.txt",
         "making_old": f"u.bubble-holding.{dead.pid}-tmp.1.0", "making_new": f"t.bubble-holding.{dead.pid}-tmp.1.0"}
for f in fakes.values(): open(os.path.join(k, f), "w").close()
old = time.time() - 120
os.utime(os.path.join(k, fakes["making_old"]), (old, old))
lock = os.open(os.path.join(k, fakes["locked"]), os.O_RDWR)
fcntl.flock(lock, fcntl.LOCK_EX)
rc = subprocess.run([bin_, "bubble", "-i", gfa, "-r", "synref", "-o", os.path.join(k, "x2"), "--quiet"],
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode
left = set(os.listdir(k))
os.close(lock)
gone = [f for f in killed + [fakes["stale"], fakes["making_old"]] if f in left]
kept = [f for f in (fakes["live"], fakes["locked"], fakes["other"], fakes["making_new"]) if f in left]
print(f"rc={rc}, stale holds left: {gone}, others kept: {len(kept)}/4")
sys.exit(0 if rc == 0 and not gone and len(kept) == 4 and "x2.bubbles.csv" in left else 1)
PY
cli_check $ok "a later run removes the hold of a killed run and keeps those of live runs"
# A command that reads its graph twice cannot take it on standard input: refused by name, or, where the
# system reopens /dev/stdin as the file itself, accepted and run normally.
for cmd in panphorte refine; do
  rc=0
  "$PANVAR_BIN" "$cmd" -i /dev/stdin --bubble-prefix "$OUT_DIR/p.bub" --reference-path synref \
    -o "$OUT_DIR/stdin_$cmd" --quiet < "$OUT_DIR/p.bub.sorted.gfa" > /dev/null 2> "$OUT_DIR/stdin_$cmd.err" || rc=$?
  ok=0; [[ $rc -eq 0 ]] || grep -q "reads it more than once" "$OUT_DIR/stdin_$cmd.err" || ok=1
  cli_check $ok "$cmd -i /dev/stdin is refused as a re-read input, or runs (rc=$rc)"
done
[[ $CLI_FAILS -eq 0 ]] || echo "CLI CHECKS FAILED ($CLI_FAILS)"

# The record checks run whatever the CLI checks gave, so one run reports every failure.
REC_RC=0
python3 - "$OUT_DIR/p.call.region.vcf" "$OUT_DIR/w.call.region.vcf" "$OUT_DIR/ps.call.region.vcf" \
  "$OUT_DIR/p.nocn.region.vcf" "$DATA_DIR/syn.gfa" "$DATA_DIR/syn_w.gfa" <<'PY' || REC_RC=1
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

# Where each node lies on the reference walk (original node ids), for the expected positions.
def ref_layout(gfa):
    lens, walk = {}, None
    for line in open(gfa):
        f = line.rstrip("\n").split("\t")
        if f[0] == "S": lens[f[1]] = len(f[2])
        elif f[0] == "P" and f[1].startswith("synref"): walk = [s[:-1] for s in f[2].split(",")]
        elif f[0] == "W" and f[1] == "synref":
            walk = [s for s in f[6].replace(">", " ").replace("<", " ").split()]
    spans, at = {}, 0
    for n in walk:
        spans.setdefault(n, []).append((at, at + lens[n])); at += lens[n]
    return spans
def pos_end(f): return int(f[1]), int(info(f)["END"])

# ---- P-line graph: the seven isolated events + confounders --------------------
hdr, recs = load(sys.argv[1])
by = {}
for f in recs:
    by.setdefault((info(f).get("SVTYPE"), frozenset(carriers(hdr, f))), []).append(f)
def find(sv, c): return by.get((sv, frozenset(c)), [])

# syn.gfa's reference is named chrS:1-2000, read as a 0-based start of 1. A record's POS is the base
# before the event (its 0-based offset, plus the start) and END its last base; an insertion sits at
# the last base of the node it follows.
sp = ref_layout(sys.argv[5]); S = 1
def at(n, k=0): return S + sp[n][k][0]
def to(n, k=-1): return S + sp[n][k][1]
want_p = {("DEL", "s_del"): (at("10"), to("10")), ("INS", "s_ins"): (to("9"), to("9")),
          ("INV", "s_inv"): (at("12"), to("12")), ("CNV", "s_selfdup"): (at("13"), to("13")),
          ("CNV", "s_peakdup1"): (at("15"), to("15")), ("DEL", "s_confdel"): (at("18"), to("18"))}
for (sv, hap), pe in want_p.items():
    got = [pos_end(f) for f in recs if info(f).get("SVTYPE") == sv and hap in carriers(hdr, f)]
    check(got == [pe], f"[P] {sv} carried by {hap} at POS,END {pe} (got {got})")
p_records = sorted((f[1], info(f)["END"], info(f)["SVTYPE"], tuple(carriers(hdr, f))) for f in recs)

check(len(find("DEL", ["s_del"])) == 1, "[P] DEL of DEL1 carried by s_del")
ins = find("INS", ["s_ins"])
check(len(ins) == 1 and info(ins[0]).get("INS_SUBTYPE") == "NOVEL",
      "[P] INS of INS1 carried by s_ins, INS_SUBTYPE=NOVEL")
check(len(find("INV", ["s_inv"])) == 1, "[P] INV of INV1 carried by s_inv")
sd = find("CNV", ["s_selfdup"])
check(len(sd) == 1 and info(sd[0]).get("REF_CN") == "1",
      "[P] self-loop DUP carried by s_selfdup, REF_CN=1")
# Under --cn the folded G module (G-C-G, reference G x2, no self-loop) is one total-module DUP REF_CN=2
# carrying the gains AND the loss s_null: the loss folds into the per-sample CN, so it is not also a
# DEL, and the duplicated content folds into the DUP, so it is not also an INS.
gpk = [f for f in recs if info(f).get("SVTYPE") == "CNV" and info(f).get("REF_CN") == "2"
       and {"s_peakdup1", "s_peakdup2", "s_null"} <= set(carriers(hdr, f))]
check(len(gpk) == 1, "[P] folded G module is one DUP REF_CN=2 carrying the gains and the loss s_null")
check(len(gpk) == 1 and info(gpk[0]).get("CN_METHOD") == "MODULE", "[P] the G module DUP is CN_METHOD=MODULE")
check(not any(info(f).get("SVTYPE") == "DEL" and "s_null" in carriers(hdr, f) for f in recs),
      "[P] s_null copy-number loss is in the DUP, not also a DEL")
gbub = info(gpk[0]).get("BUBBLE_ID") if gpk else None
check(not any(info(f).get("SVTYPE") == "INS" and info(f).get("BUBBLE_ID") == gbub
             and {"s_peakdup1", "s_peakdup2"} & set(carriers(hdr, f)) for f in recs),
      "[P] the duplicated content folds into the DUP, not also an INS")
check(len([f for f in recs if info(f).get("SVTYPE") == "CNV"]) == 2,
      "[P] exactly two DUP records (self-loop + folded G module); invariant fold B yields no DUP")
for hap in ["s_selfdup", "s_peakdup1", "s_peakdup2"]:
    n = sum(1 for f in recs if hap in carriers(hdr, f))
    check(n == 1, f"[P] {hap} appears in exactly one record (no double-count), got {n}")

# ---- W-line graph: parsing + mixed-orientation multi-copy ---------------------
hdr, recs = load(sys.argv[2])
inv = [f for f in recs if info(f).get("SVTYPE") == "INV" and "s_winv" in carriers(hdr, f)]
check(len(inv) == 1, "[W] INV called from W-line walk (s_winv)")
spw = ref_layout(sys.argv[6])   # the W line starts at 0
check([pos_end(f) for f in inv] == [(spw["4"][0][0], spw["4"][0][1])], "[W] INV at the inverted node's span")
mix = [f for f in recs if info(f).get("SVTYPE") == "CNV" and "s_wmixdup" in carriers(hdr, f)]
check(len(mix) == 1 and info(mix[0]).get("REF_CN") == "2",
      "[W] mixed-orientation (+,-,+) multi-copy is a single DUP, REF_CN=2")
check([pos_end(f) for f in mix] == [(spw["5"][0][0], spw["5"][-1][1])],
      "[W] the DUP spans the reference's copies of the unit")
n = sum(1 for f in recs if "s_wmixdup" in carriers(hdr, f))
check(n == 1, f"[W] s_wmixdup in exactly one record (mixed orientation not double-counted / not an INV), got {n}")

# ---- sites from vg snarls: the same records ---------------------------------------
hdr, recs = load(sys.argv[3])
ps_records = sorted((f[1], info(f)["END"], info(f)["SVTYPE"], tuple(carriers(hdr, f))) for f in recs)
check(ps_records == p_records, "[snarls] the same records (POS, END, type, carriers) as the main P run")

# ---- without --cn: no module record; the gains are an INS and the loss a DEL --------
hdr, recs = load(sys.argv[4])
check(not any(info(f).get("CN_METHOD") == "MODULE" for f in recs), "[no --cn] no CN_METHOD=MODULE record")
check(any(info(f).get("SVTYPE") == "INS" and {"s_peakdup1", "s_peakdup2"} <= set(carriers(hdr, f)) for f in recs),
      "[no --cn] the G module gains are an INS carried by s_peakdup1 and s_peakdup2")
check(any(info(f).get("SVTYPE") == "DEL" and carriers(hdr, f) == ["s_null"] for f in recs),
      "[no --cn] the G module loss is a DEL carried by s_null")

print()
if fails:
    print(f"SYNTHETIC SMOKE FAILED ({len(fails)} assertion(s))")
    sys.exit(1)
print("synthetic smoke: records OK")
PY
[[ $CLI_FAILS -eq 0 && $REC_RC -eq 0 ]] || exit 1
echo "synthetic smoke: OK"
