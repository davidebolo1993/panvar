#!/usr/bin/env bash
# rebuild_stats.sh - contract assertions for `panvar rebuild`.
#
# These assert the properties a caller relies on: that the input is never destroyed, that a healthy graph
# comes back byte-identical, that a malformed graph is refused rather than silently repaired, that each
# repair keeps what the haplotypes hold, and that nothing is left behind.
#
#   rebuild_stats.sh <panvar-binary> <out-dir> [input.gfa]
set -uo pipefail

BIN="${1:?usage: rebuild_stats.sh <panvar> <outdir> [input.gfa]}"
PY="${PYTHON:-python3}"
OUT="${2:?}"
SRC="${3:-}"
# A UNIQUE directory per invocation, under the requested output path. Two runs sharing one directory is
# not merely untidy: a concurrent invocation removing it while minigraph is reading scratch FASTAs
# segfaults the reader, which looks exactly like a rebuild race and is not one. The parent is left
# alone so a caller's directory is never removed.
mkdir -p "$OUT"
OUT="$OUT/run.$$.$(date +%s)"
rm -rf "$OUT"
mkdir -p "$OUT"
fails=0
ok()  { printf "  ok   %s\n" "$1"; }
bad() { printf "  FAIL %s\n" "$1"; fails=$((fails + 1)); }

# A small but real graph: a handful of haplotypes over shared nodes.
GFA="$OUT/in.gfa"
if [ -n "$SRC" ] && [ -s "$SRC" ]; then
  case "$SRC" in *.gz) gzip -dc "$SRC" > "$GFA" ;; *) cp "$SRC" "$GFA" ;; esac
else
  {
    printf 'H\tVN:Z:1.0\n'
    printf 'S\t1\tACGTACGTACGTACGTACGT\nS\t2\tTTTTGGGGCCCCAAAATTTT\nS\t3\tGGGGCCCCAAAATTTTGGGG\n'
    printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
    printf 'P\thapA\t1+,2+,3+\t*\nP\thapB\t1+,3+\t*\n'
  } > "$GFA"
fi

# ---------------------------------------------------------------- the input must survive
before=$(wc -c < "$GFA")
"$BIN" rebuild -i "$GFA" -o "$GFA" >/dev/null 2>&1
after=$(wc -c < "$GFA")
[ "$before" = "$after" ] && [ "$after" -gt 0 ] \
  && ok "same --gfa and --out leaves the input intact ($after bytes)" \
  || bad "same --gfa and --out changed the input: $before -> $after bytes"
e=$("$BIN" rebuild -i "$GFA" -o "$GFA" 2>&1 | grep -c "same file")
[ "$e" -gt 0 ] && ok "same --gfa and --out is refused with a clear message" \
               || bad "same --gfa and --out was not refused"

ln -sf "$GFA" "$OUT/link.gfa"
before=$(wc -c < "$GFA")
"$BIN" rebuild -i "$GFA" -o "$OUT/link.gfa" >/dev/null 2>&1
[ "$(wc -c < "$GFA")" = "$before" ] && ok "a symlink to the input is refused too" \
                                    || bad "writing through a symlink destroyed the input"

# ---------------------------------------------------------------- healthy graph passes through
"$BIN" rebuild -i "$GFA" -o "$OUT/pass.gfa" >/dev/null 2>&1
if cmp -s "$GFA" "$OUT/pass.gfa"; then
  ok "a healthy graph passes through byte-identically"
else
  # The gate may legitimately rebuild a synthetic graph; then require every path to survive instead.
  pin=$(grep -c '^P' "$GFA"); pout=$(grep -c '^P' "$OUT/pass.gfa" 2>/dev/null || echo 0)
  [ "$pin" = "$pout" ] && ok "graph was rebuilt, and every path survived ($pout/$pin)" \
                       || bad "path count changed: $pin -> $pout"
fi
left=$(find "$OUT" -maxdepth 1 -name '*rebuild-tmp*' -o -maxdepth 1 -name '*rebuild.tmp*' | wc -l | tr -d ' ')
[ "$left" = "0" ] && ok "no staging files or scratch directories left behind" \
                  || bad "$left temporary artefacts left in the output directory"

# ---------------------------------------------------------------- malformed input is refused
printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTAC\nS\t2\tTTTTGGGGCC\nP\tp1\t1+,2+,9+\t*\n' > "$OUT/e_missing.gfa"
printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTAC\nS\t2\tTTTTGGGGCC\nP\tp1\t1+,2+\t*\nP\tp1\t2+,1+\t*\n' > "$OUT/e_dup.gfa"
printf 'H\tVN:Z:1.0\nS\t1\t*\nS\t2\tTTTTGGGGCC\nP\tp1\t1+,2+\t*\n' > "$OUT/e_nostar.gfa"
printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTAC\nS\t2\tTTTTGGGGCC\nL\t1\t+\t2\t+\t0M\n' > "$OUT/e_nopath.gfa"
for c in "e_missing:missing node" "e_dup:duplicate path name" "e_nostar:node with no sequence" \
         "e_nopath:no paths"; do
  f=${c%%:*}; what=${c##*:}
  e=$("$BIN" rebuild -i "$OUT/$f.gfa" -o "$OUT/$f.out.gfa" 2>&1 | grep -c "Error")
  [ "$e" -gt 0 ] && ok "refused: $what" || bad "$what was accepted"
  [ ! -s "$OUT/$f.out.gfa" ] && ok "  ... and wrote no output" || bad "$what left an output file behind"
done

# ---------------------------------------------------------------- the acceptance contract
# minigraph augments variation above --min-var, so sub-threshold differences are collapsed by
# construction and a recovered walk is never byte-identical. The contract is therefore structural plus
# a threshold: every path back, every step pair backed by an edge, identity/cover above the bounds --
# otherwise the rebuilt graph is discarded and the ORIGINAL is written.
if [ -n "$SRC" ] && [ -s "$SRC" ]; then
  in_nodes=$(grep -c '^S' "$GFA"); in_paths=$(grep -c '^P' "$GFA")

  # An identity bound nothing can satisfy must reject and roll back.
  "$BIN" rebuild -i "$GFA" -o "$OUT/strict.gfa" --force --min-recovered-identity 0.999999 \
    > "$OUT/strict.log" 2>&1
  grep -q "warning: rejected:" "$OUT/strict.log" && ok "an unsatisfiable contract is rejected" \
                                       || bad "an unsatisfiable contract was accepted"
  [ "$(grep -c '^S' "$OUT/strict.gfa")" = "$in_nodes" ] \
    && ok "  ... and the ORIGINAL graph was written ($in_nodes nodes)" \
    || bad "rollback did not restore the original graph"
  [ "$(grep -c '^P' "$OUT/strict.gfa")" = "$in_paths" ] \
    && ok "  ... with every path intact ($in_paths)" || bad "rollback lost paths"

  # A deletion of a few kb the generator leaves out near a haplotype's start is written into its walk, so
  # the forced rebuild is accepted with the default bounds.
  "$BIN" rebuild -i "$GFA" -o "$OUT/indel.gfa" --force -t 2 > "$OUT/indel.log" 2>&1
  r=$(awk -F'\t' '$1=="#verdict"{v=$2} $1=="#indels_added"{n=$2} END{print v, (n>0 ? "added" : "none")}' "$OUT/indel.rebuild.report.tsv")
  [ "$r" = "accepted added" ] && ok "a deletion the generator left out is written into the walk; the rebuild is accepted" \
                              || bad "the forced rebuild was not accepted with the indel written in: $r"

  # --allow-loss accepts the same run and says so.
  "$BIN" rebuild -i "$GFA" -o "$OUT/loss.gfa" --force --min-recovered-identity 0.999999 \
    --allow-loss > "$OUT/loss.log" 2>&1
  grep -q "warning: accepted with --allow-loss" "$OUT/loss.log" && ok "--allow-loss accepts but warns" \
                                    || bad "--allow-loss did not warn"

  # A reference that is not in the graph must reject, whatever else passes.
  "$BIN" rebuild -i "$GFA" -o "$OUT/badref.gfa" --force -r "NOT_A_REAL_PATH_XYZ" \
    > "$OUT/badref.log" 2>&1
  grep -q "reference path not found" "$OUT/badref.log" \
    && ok "a missing reference path is rejected" || bad "a missing reference path was accepted"

  # The report sidecar must have one row per path and name a status for each.
  aud="$OUT/strict.rebuild.report.tsv"
  if [ -s "$aud" ]; then
    # data rows only: the trailing #verdict/#reason lines are metadata, not paths
    rows=$(grep -vc '^#' "$aud"); rows=$(( rows - 1 ))
    [ "$rows" = "$in_paths" ] && ok "report has one row per path ($rows)" \
                              || bad "report has $rows rows for $in_paths paths"
    bad_status=$(awk -F'\t' 'NR>1 && $0 !~ /^#/ && $NF!="ok" && $NF!="not_recovered" && $NF!="low_cover" && $NF!="low_identity" && $NF!="identity_unavailable" && $NF!="assembly_gap" && $NF!="reference_inexact" && $NF!="end_overrun" && $NF!="end_short" && $NF!="length_mismatch" && $NF!="possible_inversion" && $NF!="dangling_walk" && $NF!="ok_after_fold" && $NF!="fails_after_fold"' "$aud" | wc -l | tr -d ' ')
    [ "$bad_status" = "0" ] && ok "every report row carries a known status" \
                            || bad "$bad_status report rows have an unrecognised status"
  else
    bad "no report sidecar was written"
  fi
fi


# ---------------------------------------------------------------- reproducibility
# A rebuild that depends on thread count or on which strand the input happens to be stored on is not
# reproducible, and neither property is visible in a single run's output.
if [ -n "$SRC" ] && [ -s "$SRC" ]; then
  seed_of() { "$BIN" rebuild -i "$1" -o "$2" --force --min-recovered-identity 0.97 2>&1 \
                | grep -oE 'seed: [^ ]*' | head -1; }
  seed_of_k() { "$BIN" rebuild -i "$1" -o "$2" --force --min-recovered-identity 0.97 --kmer "$3" 2>&1 \
                  | grep -oE 'seed: [^ ]*' | head -1; }

  s1=$(seed_of "$GFA" "$OUT/rep1.gfa"); s2=$(seed_of "$GFA" "$OUT/rep2.gfa")
  [ -n "$s1" ] && [ "$s1" = "$s2" ] && ok "repeated runs choose the same seed" \
                                    || bad "seed differs between runs: '$s1' vs '$s2'"

  "$BIN" rebuild -i "$GFA" -o "$OUT/th1.gfa" --force --min-recovered-identity 0.97 -t 1 -q >/dev/null 2>&1
  "$BIN" rebuild -i "$GFA" -o "$OUT/th8.gfa" --force --min-recovered-identity 0.97 -t 8 -q >/dev/null 2>&1
  cmp -s "$OUT/th1.gfa" "$OUT/th8.gfa" && ok "output is identical with 1 and 8 threads" \
                                       || bad "thread count changes the output"
  cmp -s "$OUT/th1.rebuild.report.tsv" "$OUT/th8.rebuild.report.tsv" \
    && ok "  ... and so is the report" || bad "thread count changes the report"

  # Reverse-complement every segment and flip every step: the same graph, stored on the other strand.
  # k-mer richness must be canonical, or the seed (and so the whole rebuild) follows the storage.
  "$PY" - "$GFA" "$OUT/rc.gfa" <<'PYEOF'
import sys
src, dst = sys.argv[1], sys.argv[2]
comp = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')
def rc(s): return s.translate(comp)[::-1]
flip = {'+': '-', '-': '+'}
out = []
for line in open(src):
    f = line.rstrip('\n').split('\t')
    if f[0] == 'S' and len(f) > 2:
        f[2] = rc(f[2])
    elif f[0] == 'P' and len(f) > 2:
        f[2] = ','.join(t[:-1] + flip.get(t[-1], t[-1]) for t in f[2].split(','))
    elif f[0] == 'L' and len(f) > 4:
        f[2], f[4] = flip.get(f[2], f[2]), flip.get(f[4], f[4])
    out.append('\t'.join(f))
open(dst, 'w').write('\n'.join(out) + '\n')
PYEOF
  s3=$(seed_of "$OUT/rc.gfa" "$OUT/rc.out.gfa")
  [ "$s1" = "$s3" ] && ok "seed is the same on the reverse-complemented graph (canonical k-mers)" \
                    || bad "seed follows the strand the graph is stored on: '$s1' vs '$s3'"
fi

# ---------------------------------------------------------------- ambiguous bases
# N must neither crash the k-mer roll nor be counted as sequence.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tACGTNNNNACGTACGTACGT\nS\t2\tTTTTGGGGCCCCAAAANNNN\nS\t3\tGGGGCCCCAAAATTTTGGGG\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\thapA\t1+,2+,3+\t*\nP\thapB\t1+,3+\t*\n'; } > "$OUT/amb.gfa"
"$BIN" rebuild -i "$OUT/amb.gfa" -o "$OUT/amb.out.gfa" >/dev/null 2>&1
rc=$?
[ "$rc" -eq 0 ] && [ -s "$OUT/amb.out.gfa" ] && ok "ambiguous bases (N) are handled without failure" \
                                             || bad "a graph containing N failed (exit $rc)"


# ---------------------------------------------------------------- traversal and overlap validation
# A path that steps between nodes with no link spells a sequence no walk could produce; a non-zero
# overlap would double-count the overlapping bases in every length and identity figure downstream.
printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTAC\nS\t2\tTTTTGGGGCC\nS\t3\tGGGGCCCCAA\nL\t1\t+\t2\t+\t0M\nP\tp1\t1+,2+,3+\t*\n' > "$OUT/e_nolink.gfa"
printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTAC\nS\t2\tTTTTGGGGCC\nL\t1\t+\t2\t+\t5M\nP\tp1\t1+,2+\t*\n' > "$OUT/e_ovl.gfa"
for c in "e_nolink:a step pair with no link" "e_ovl:a non-zero overlap"; do
  f=${c%%:*}; what=${c##*:}
  e=$("$BIN" rebuild -i "$OUT/$f.gfa" -o "$OUT/$f.out.gfa" 2>&1 | grep -c "Error")
  [ "$e" -gt 0 ] && ok "refused: $what" || bad "$what was accepted"
done
# ...but a VALID reverse traversal must still pass.
printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTAC\nS\t2\tTTTTGGGGCC\nL\t1\t+\t2\t-\t0M\nP\tp1\t1+,2-\t*\nP\tp2\t1+,2-\t*\n' > "$OUT/rev.gfa"
"$BIN" rebuild -i "$OUT/rev.gfa" -o "$OUT/rev.out.gfa" >/dev/null 2>&1
[ -s "$OUT/rev.out.gfa" ] && ok "a valid reverse-orientation traversal is accepted" \
                          || bad "a valid reverse traversal was rejected"

# ---------------------------------------------------------------- k > 31 uses the same rules as k <= 31
# Both branches must be canonical and ambiguity-aware, or crossing k=31 silently changes what richness
# means. Same graph, two k values either side of the boundary: the seed must not change.
if [ -n "$SRC" ] && [ -s "$SRC" ]; then
  # The k>31 branch must be canonical and ambiguity-aware like the k<=31 one. Two seeds agreeing at
  # different k is incidental -- the property is that at a GIVEN k above 31 the ranking does not follow
  # the strand the graph is stored on, and that ambiguity does not break it. Test that directly.
  for kk in 21 41; do
    a1=$(seed_of_k "$GFA" "$OUT/kf.$kk.gfa" "$kk")
    a2=$(seed_of_k "$OUT/rc.gfa" "$OUT/kr.$kk.gfa" "$kk")
    [ -n "$a1" ] && [ "$a1" = "$a2" ] \
      && ok "k=$kk: seed is strand-independent (canonical on both branches)" \
      || bad "k=$kk: seed follows the stored strand: '$a1' vs '$a2'"
  done

  # An ambiguous --reference-path must be refused, not resolved by file order.
  e=$("$BIN" rebuild -i "$GFA" -o "$OUT/amb.out.gfa" --force -r "haplotype" 2>&1 | grep -c "ambiguous")
  [ "$e" -gt 0 ] && ok "an ambiguous --reference-path is refused" \
                 || bad "an ambiguous --reference-path was resolved silently"

  # The report must carry the global verdict beside the per-path rows.
  "$BIN" rebuild -i "$GFA" -o "$OUT/verd.gfa" --force --min-recovered-identity 0.999999 >/dev/null 2>&1
  grep -q "^#verdict" "$OUT/verd.rebuild.report.tsv" 2>/dev/null \
    && ok "report records the global verdict and reason" || bad "report has no verdict row"
fi

# ---------------------------------------------------------------- equal richness is ordered, not arbitrary
# Two byte-identical haplotypes tie on every richness measure. The seed must be decided by name, not by
# whatever std::sort happens to do.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tACGTACGTACGTACGTACGT\nS\t2\tTTTTGGGGCCCCAAAATTTT\nS\t3\tGGGGCCCCAAAATTTTGGGG\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tzzz\t1+,2+,3+\t*\nP\taaa\t1+,2+,3+\t*\nP\tmmm\t1+,3+\t*\n'; } > "$OUT/tie.gfa"
t1=$("$BIN" rebuild -i "$OUT/tie.gfa" -o "$OUT/tie1.gfa" --force 2>&1 | grep -oE 'seed: [^ ]*' | head -1)
t2=$("$BIN" rebuild -i "$OUT/tie.gfa" -o "$OUT/tie2.gfa" --force 2>&1 | grep -oE 'seed: [^ ]*' | head -1)
[ -n "$t1" ] && [ "$t1" = "$t2" ] && ok "equal richness is broken deterministically ($t1)" \
                                  || bad "equal-richness seed is unstable: '$t1' vs '$t2'"
case "$t1" in *aaa*) ok "  ... and by path name, not input order" ;;
              *) bad "tie broken by something other than name: $t1" ;; esac


# ---------------------------------------------------------------- disposition, precision, units
if [ -n "$SRC" ] && [ -s "$SRC" ]; then
  # --allow-loss WRITES the rebuilt graph, so the report must not call that "rejected".
  "$BIN" rebuild -i "$GFA" -o "$OUT/ovr.gfa" --force --min-recovered-identity 0.999999 --allow-loss \
    > /dev/null 2>&1
  v=$(awk -F'\t' '$1=="#verdict"{print $2}' "$OUT/ovr.rebuild.report.tsv" 2>/dev/null)
  [ "$v" = "accepted_with_override" ] && ok "--allow-loss records disposition '$v'" \
                                      || bad "--allow-loss recorded verdict '$v', expected accepted_with_override"
  # The recorded threshold must be precise enough to reproduce the decision: 0.999999 must not print
  # as 1.0000.
  t=$(awk -F'\t' '$1=="#min_recovered_identity"{print $2}' "$OUT/ovr.rebuild.report.tsv" 2>/dev/null)
  case "$t" in 0.999999*) ok "report records the threshold at full precision ($t)" ;;
               *) bad "report records threshold as '$t', losing the digits that decided it" ;; esac
fi

# A forced rebuild of a graph with no variation.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tACGTACGTACGTACGTACGT\nS\t2\tTTTTGGGGCCCCAAAATTTT\n'
  printf 'L\t1\t+\t2\t+\t0M\n'
  printf 'P\tp1\t1+,2+\t*\nP\tp2\t1+,2+\t*\n'; } > "$OUT/flat.gfa"
"$BIN" rebuild -i "$OUT/flat.gfa" -o "$OUT/flat.out.gfa" --force > "$OUT/flat.log" 2>&1

# An unknown ('*') overlap means unverified, not zero.
printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTAC\nS\t2\tTTTTGGGGCC\nL\t1\t+\t2\t+\t*\nP\tp1\t1+,2+\t*\n' > "$OUT/e_star.gfa"
e=$("$BIN" rebuild -i "$OUT/e_star.gfa" -o "$OUT/e_star.out.gfa" 2>&1 | grep -c "UNKNOWN overlap")
[ "$e" -gt 0 ] && ok "an unknown ('*') overlap is refused, not read as 0M" \
               || bad "an unknown ('*') overlap was accepted as zero"


# ---------------------------------------------------------------- a walk must not overrun its haplotype
# A walk lists whole segments and the generator never cuts one where a haplotype ends. Here the seed (the
# richest haplotype) carries a 1.5 kb insertion 1 kb from its end, which nobody else has, and the others'
# walks spell it. Where a haplotype's last 1 kb matches the walk's, the insertion is a deletion written into
# its walk (endins); where it does not (endins2: the others end in other sequence), the walk still runs past
# the haplotype, and the overrun bound must refuse that.
"$PY" - "$OUT/endins.gfa" "$OUT/endins2.gfa" <<'PYEOF'
import random, sys
r = random.Random(7)
base = ''.join(r.choice('ACGT') for _ in range(100000))
ins = ''.join(r.choice('ACGT') for _ in range(1500))
other = ''.join(r.choice('ACGT') for _ in range(1000))
def write(path, seqs):
    with open(path, 'w') as o:
        o.write('H\tVN:Z:1.0\n')
        for i, (n, s) in enumerate(seqs.items(), 1):
            o.write(f'S\t{i}\t{s}\n')
        for i, n in enumerate(seqs, 1):
            o.write(f'P\t{n}\t{i}+\t*\n')
# h1's 2 kb deletion mid-locus gives the generator a real event, so the graph is not degenerate.
write(sys.argv[1], {'ref': base, 'h1': base[:50000] + base[52000:], 'h2': base, 'h3': base[:99000] + ins + base[99000:]})
b2 = base[:99000] + other
write(sys.argv[2], {'ref': b2, 'h1': b2[:50000] + b2[52000:], 'h2': b2, 'h3': base[:99000] + ins + base[99000:]})
PYEOF
"$BIN" rebuild -i "$OUT/endins.gfa" -o "$OUT/endins.out.gfa" --force -t 2 > "$OUT/endins.log" 2>&1
a=$(awk -F'\t' '$1=="#verdict"{v=$2} NR>1 && $0 !~ /^#/ && $7==1 {n++} END{print v, n+0}' "$OUT/endins.out.rebuild.report.tsv" 2>/dev/null)
[ "$a" = "accepted 4" ] && ok "a seed-private insertion near the end is skipped by the others' walks, all 4 exact" \
                        || bad "seed-private insertion near the end: '$a', expected 'accepted 4'"
"$BIN" rebuild -i "$OUT/endins2.gfa" -o "$OUT/endins2.out.gfa" --force -t 2 > "$OUT/endins2.log" 2>&1
v=$(awk -F'\t' '$1=="#verdict"{print $2}' "$OUT/endins2.out.rebuild.report.tsv" 2>/dev/null)
n=$(awk -F'\t' '$NF=="end_overrun" || $NF=="length_mismatch"' "$OUT/endins2.out.rebuild.report.tsv" 2>/dev/null | wc -l | tr -d ' ')
[ "$v" = "rejected" ] && [ "$n" -gt 0 ] \
  && ok "seed-private end sequence copied into other walks is refused ($n walk(s) overrun)" \
  || bad "seed-private end sequence was accepted (verdict '$v', $n overrunning walks)"

# ...whereas a haplotype that genuinely stops early (a truncated assembly) gets its end segment cut, so
# its walk ends where it does, and the rebuild is accepted with that walk exact.
"$PY" - "$OUT/trunc.gfa" <<'PYEOF'
import random, sys
r = random.Random(11)
base = ''.join(r.choice('ACGT') for _ in range(100000))
seqs = {'ref': base, 'h1': base[:50000] + base[52000:], 'h2': base[:97000]}
with open(sys.argv[1], 'w') as o:
    o.write('H\tVN:Z:1.0\n')
    for i, (n, s) in enumerate(seqs.items(), 1):
        o.write(f'S\t{i}\t{s}\n')
    for i, n in enumerate(seqs, 1):
        o.write(f'P\t{n}\t{i}+\t*\n')
PYEOF
"$BIN" rebuild -i "$OUT/trunc.gfa" -o "$OUT/trunc.out.gfa" --force -t 2 > "$OUT/trunc.log" 2>&1
v=$(awk -F'\t' '$1=="#verdict"{print $2}' "$OUT/trunc.out.rebuild.report.tsv" 2>/dev/null)
cuts=$(awk -F'\t' '$1=="#end_cuts"{print $2}' "$OUT/trunc.out.rebuild.report.tsv" 2>/dev/null)
wid=$(awk -F'\t' '$1=="h2"{print $7}' "$OUT/trunc.out.rebuild.report.tsv" 2>/dev/null)
[ "$v" = "accepted" ] && [ "${cuts:-0}" -gt 0 ] && [ "$wid" = "1" ] \
  && ok "a truncated haplotype's end segment is cut and its walk is exact ($cuts cut)" \
  || bad "truncated haplotype: verdict '$v', cuts '${cuts:-}', walk identity '$wid'"

# ---------------------------------------------------------------- an indel near a haplotype's end
# The generator drops a difference near a haplotype's end, and one of --min-var or more is under the length
# bound there. A 400 bp insertion 100 bp from the start and a 400 bp deletion 300 bp before the end are written
# into their walks, so both come back exact; one whose outer bases differ from the walk cannot be placed, and
# its haplotype is flagged rather than passed.
"$PY" - "$OUT/nearend.gfa" "$OUT/divend.gfa" <<'PYEOF'
import random, sys
r = random.Random(7)
rnd = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
def mutate(s, rate): return ''.join(r.choice([c for c in 'ACGT' if c != x]) if r.random() < rate else x for x in s)
def write(path, haps):
    with open(path, 'w') as o:
        o.write('H\tVN:Z:1.0\n')
        for i, (n, s) in enumerate(haps.items(), 1): o.write(f'S\t{i}\t{s}\n')
        for i, n in enumerate(haps, 1): o.write(f'P\t{n}\t{i}+\t*\n')
b = rnd(100000)
write(sys.argv[1], {'ref': b, 'h1': b[:50000] + b[52000:], 'ins': b[:100] + rnd(400) + b[100:],
                    'del': b[:99300] + b[99700:]})
write(sys.argv[2], {'ref': b, 'h1': b[:50000] + b[52000:], 'div': mutate(b[:200], 0.3) + rnd(400) + b[200:]})
PYEOF
"$BIN" rebuild -i "$OUT/nearend.gfa" -o "$OUT/nearend.out.gfa" --force -r ref -t 2 > /dev/null 2>&1
A="$OUT/nearend.out.rebuild.report.tsv"
a=$(awk -F'\t' '$1=="#verdict"{v=$2} $1=="#indels_added"{n=$2} ($1=="ins"||$1=="del"){w=w" "$1":"$7":"$NF} END{print v, n w}' "$A" 2>/dev/null)
[ "$a" = "accepted 2 ins:1:ok del:1:ok" ] \
  && ok "a 400 bp insertion and deletion near the ends are written into their walks, both exact" \
  || bad "near-end indels: '$a', expected 'accepted 2 ins:1:ok del:1:ok'"
"$BIN" rebuild -i "$OUT/divend.gfa" -o "$OUT/divend.out.gfa" --force -r ref -t 2 > /dev/null 2>&1
a=$(awk -F'\t' '$1=="#verdict"{v=$2} $1=="div"{s=$NF} END{print v, s}' "$OUT/divend.out.rebuild.report.tsv" 2>/dev/null)
[ "$a" = "rejected end_short" ] && ok "  ... one that cannot be placed is flagged end_short, not passed as ok" \
                                || bad "  ... an unplaced near-end insertion gave '$a', expected 'rejected end_short'"

# ---------------------------------------------------------------- an N run at a haplotype's end
# A haplotype starting or ending in an assembly gap inside the locus: its walk stops where the N run begins,
# or runs part way into it, and those bases are not variation the graph lost. 300 bases the graph lacks on
# either side of the N run still are, also where the walk spells other bases opposite them, and the report
# gives the shortfall; so is an end that differs from the walk's with no N run at all. An end dense in SNPs,
# or a short insertion the walk lacks, is a difference the rebuild collapses and is not short. Each layout at
# both ends, and with the haplotype stored on the other strand.
"$PY" - "$OUT" <<'PYEOF'
import random, sys
r = random.Random(7)
rnd = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
rc = lambda s: s[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))
def write(path, haps, rev=()):
    with open(path, 'w') as o:
        o.write('H\tVN:Z:1.0\n')
        for i, (n, s) in enumerate(haps.items(), 1): o.write(f'S\t{i}\t{rc(s) if n in rev else s}\n')
        for i, n in enumerate(haps, 1): o.write(f'P\t{n}\t{i}{"-" if n in rev else "+"}\t*\n')
b = rnd(100000)
x = rnd(300)
n5 = 'N' * 500
# h1's 2 kb deletion puts node boundaries at 60 and 62 kb, where the gapped haplotypes stop or start; an N
# run 400 bp short of a boundary leaves the walk spelling those 400 bases opposite it
cases = {
    'nend': b[:60000] + n5, 'nstart': n5 + b[62000:],                     # the N run at the boundary
    'ninto': b[:59600] + n5, 'ninto_s': n5 + b[62400:],                   # the walk runs into the N run
    'nmixed': b[:60000] + x + n5, 'nmixed_s': n5 + x + b[62000:],         # lost bases, then the N run
    'nover': b[:59600] + n5 + x, 'nover_s': x + n5 + b[62400:],           # the N run, then lost bases
    'xin': b[:59600] + x + n5, 'xin_s': n5 + x + b[62400:],               # lost bases opposite the walk
    'mm': b[:59600] + rnd(400), 'mm_s': rnd(400) + b[62400:],             # a divergent end, no N run
}
def snp(s, p, seed):
    rr = random.Random(seed); o = list(s)
    for i in range(len(o)):
        if rr.random() < p: o[i] = rr.choice([c for c in 'ACGT' if c != o[i]])
    return ''.join(o)
cases.update({   # SNPs at 10% over the last or first 200 bases, and at 8% over 1 kb
    'snp10': b[:99800] + snp(b[99800:], 0.10, 1), 'snp10_s': snp(b[:200], 0.10, 2) + b[200:],
    'snp8': b[:99000] + snp(b[99000:], 0.08, 3), 'snp8_s': snp(b[:1000], 0.08, 4) + b[1000:],
})
# Near the boundary: 30 inserted bases 30 bp before it, after SNPs at 4%; 45 bases the haplotype lacks 60 bp
# before it, among SNPs at 10%; 60 inserted bases 40 bp before it
e, s_ = snp(b[59400:60000], 0.10, 5), snp(b[62000:62600], 0.10, 6)
cases.update({
    'insnear': b[:59400] + snp(b[59400:59970], 0.04, 7) + rnd(30) + b[59970:60000],
    'insnear_s': b[62000:62030] + rnd(30) + snp(b[62030:62600], 0.04, 8) + b[62600:],
    'delsnp': b[:59400] + e[:495] + e[540:], 'delsnp_s': s_[:60] + s_[105:] + b[62600:],
    'ins60': b[:59960] + rnd(60) + b[59960:60000], 'ins60_s': b[62000:62040] + rnd(60) + b[62040:],
    # an N run inserted 200 bp before the boundary
    'nrep': b[:59800] + 'N' * 150 + b[59800:60000], 'nrep_s': b[62000:62200] + 'N' * 150 + b[62200:],
})
for c, gap in cases.items():
    for strand, rev in (('', ()), ('_rev', ('gap',))):
        write(f'{sys.argv[1]}/{c}{strand}.gfa', {'ref': b, 'h1': b[:60000] + b[62000:], 'gap': gap}, rev)
# a 30 bp insertion just before the end, where the walk's last bases align equally well to the inserted ones
c40, a30 = 'C' * 40, 'A' * 30
for c, haps in {'tie': {'ref': b + c40, 'h1': b[:60000] + b[62000:] + c40, 'gap': b + a30 + c40},
                'tie_s': {'ref': c40 + b, 'h1': c40 + b[:60000] + b[62000:], 'gap': c40 + a30 + b}}.items():
    for strand, rev in (('', ()), ('_rev', ('gap',))):
        write(f'{sys.argv[1]}/{c}{strand}.gfa', haps, rev)
PYEOF
NCASES="nend nstart ninto ninto_s nmixed nmixed_s nover nover_s xin xin_s mm mm_s snp10 snp10_s snp8 snp8_s tie tie_s insnear insnear_s delsnp delsnp_s ins60 ins60_s nrep nrep_s"
for c in $NCASES; do for s in "" _rev; do
  "$BIN" rebuild -i "$OUT/$c$s.gfa" -o "$OUT/$c$s.out.gfa" --force -r ref -t 2 -q >/dev/null 2>&1
done; done
# verdict, then the gapped haplotype's start_short_bp, end_short_bp and status, read by column name
short() { awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i} $1=="#verdict"{v=$2}
  $1=="gap"{s=$c["start_short_bp"]" "$c["end_short_bp"]" "$NF} END{print v, s}' "$OUT/$1.out.rebuild.report.tsv" 2>/dev/null; }
nshort() {   # nshort <expected> <case>...: each case on both strands, compared over the expected fields
  local e="$1" got="" c s f; shift
  f=$(echo "$e" | wc -w | tr -d ' ')
  for c in "$@"; do for s in "" _rev; do
    [ "$(short "$c$s" | cut -d' ' -f1-"$f")" = "$e" ] || got="$got $c$s='$(short "$c$s")'"
  done; done
  printf '%s' "$got"
}
g=$(nshort "accepted 0 0 ok" nend nstart ninto ninto_s nrep nrep_s)
[ -z "$g" ] && ok "an N run at a haplotype's start or end is not counted as a walk stopping short, either strand" \
            || bad "an N run at an end:$g, expected 'accepted 0 0 ok'"
# verdict and shortfalls only: the status names one failing check, and at the start the shorter haplotype also
# falls under the identity bound. A few of the lost bases may match the walk's by chance.
lost() {   # lost <start|end> <least> <case>...: rejected, with at least <least> bp short at that end only
  local side=$1 least=$2 got="" c s v a z; shift 2
  for c in "$@"; do for s in "" _rev; do
    read -r v a z _ <<< "$(short "$c$s")"
    [ "$side" = end ] && { local t=$a; a=$z; z=$t; }
    [ "$v" = rejected ] && [ "${a:-0}" -ge "$least" ] 2>/dev/null && [ "$z" = 0 ] || got="$got $c$s='$(short "$c$s")'"
  done; done
  printf '%s' "$got"
}
g="$(lost end 280 nmixed nover xin)$(lost start 280 nmixed_s nover_s xin_s)"
[ -z "$g" ] && ok "  ... while 300 bases the graph lacks before or past it are, with the shortfall in the report" \
            || bad "  ... 300 lost bases beside an end N run:$g, expected 'rejected' with a 300 bp shortfall"
g="$(lost end 390 mm)$(lost start 390 mm_s)"
[ -z "$g" ] && ok "  ... and so are 400 end bases that differ from the bases the walk spells opposite them" \
            || bad "  ... a 400 bp divergent end:$g, expected 'rejected' with a 400 bp shortfall"
g=$(nshort "accepted 0 0 ok" snp10 snp10_s snp8 snp8_s delsnp delsnp_s)
[ -z "$g" ] && ok "  ... but an end with SNPs at 10% over 200 bp or 8% over 1 kb, or bases the haplotype lacks, is not short" \
            || bad "  ... an end dense in SNPs, or one the haplotype lacks bases of:$g, expected 'accepted 0 0 ok'"
g="$(lost end 54 ins60)$(lost start 54 ins60_s)"
[ -z "$g" ] && ok "  ... while 60 bases inserted 40 bp from the boundary are" \
            || bad "  ... 60 inserted bases near the boundary:$g, expected 'rejected' with a 60 bp shortfall"
g=""
for c in insnear insnear_s insnear_rev insnear_s_rev; do
  read -r v a z st <<< "$(short "$c")"
  [ "$v" = accepted ] && [ "$st" = ok ] && [ $((a + z)) -le 36 ] 2>/dev/null || g="$g $c='$(short "$c")'"
done
[ -z "$g" ] && ok "  ... and 30 bases inserted 30 bp from it count as 30 at most, not with the bases after them" \
            || bad "  ... 30 inserted bases near the boundary:$g, expected accepted, at most about 30 short"
g=""
for c in tie tie_s tie_rev tie_s_rev; do
  read -r v a z st <<< "$(short "$c")"
  [ "$v" = accepted ] && [ "$st" = ok ] && [ $((a + z)) -le 30 ] 2>/dev/null || g="$g $c='$(short "$c")'"
done
[ -z "$g" ] && ok "  ... nor one with a 30 bp insertion the walk's last bases align to either side of" \
            || bad "  ... a 30 bp insertion at an end:$g, expected accepted, at most 30 short"

# ---------------------------------------------------------------- a large inversion is kept
# The generator adds an inversion only inside one alignment chain, which bridges at most ~150 kb, so a
# 200 kb inversion between inverted repeats maps as three chains. They are joined into one walk with
# two new reversing links, and kept only because the haplotype's own sequence confirms them; a fake
# (the reversed middle 10% mutated) is not joined, and is labelled.
"$PY" - "$OUT/inv.gfa" <<'PYEOF'
import random, sys
r = random.Random(5)
rnd = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
rc = lambda s: s.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
def mutate(s, rate):
    return ''.join(r.choice([c for c in 'ACGT' if c != b]) if r.random() < rate else b for b in s)
L, IR, M, R = rnd(40000), rnd(3000), rnd(200000), rnd(40000)
ref = L + IR + M + rc(IR) + R
inv = L + IR + rc(M) + rc(IR) + R
haps = {'ref': ref, 'h1': ref[:20000] + ref[22000:], 'h2': inv, 'h3': mutate(inv, 0.002),
        'h4': L + IR + mutate(rc(M), 0.10) + rc(IR) + R}
with open(sys.argv[1], 'w') as o:
    o.write('H\tVN:Z:1.0\n')
    for i, (n, s) in enumerate(haps.items(), 1): o.write(f'S\t{i}\t{s}\n')
    for i, n in enumerate(haps, 1): o.write(f'P\t{n}\t{i}+\t*\n')
PYEOF
"$BIN" rebuild -i "$OUT/inv.gfa" -o "$OUT/inv.out.gfa" --force -r ref -t 2 > "$OUT/inv.log" 2>&1
A="$OUT/inv.out.rebuild.report.tsv"
j=$(awk -F'\t' '($1=="h2"||$1=="h3") && $(NF-2)==1 && $NF=="ok"' "$A" 2>/dev/null | wc -l | tr -d ' ')
w2=$(awk -F'\t' '$1=="h2"{print $7}' "$A" 2>/dev/null)
[ "$j" = "2" ] && [ "$w2" = "1" ] && ok "a 200 kb inversion is joined into one walk (h2 exact)" \
                                  || bad "200 kb inversion not joined: $j joined, h2 identity '$w2'"
f4=$(awk -F'\t' '$1=="h4"{print $(NF-2)":"$NF}' "$A" 2>/dev/null)
[ "$f4" = "0:possible_inversion" ] && ok "  ... while a fake (10% mutated middle) is not, and is labelled" \
                                   || bad "  ... the fake inversion gave '$f4'"

# ---------------------------------------------------------------- tips and assembly gaps
# A haplotype running past the seed's start or end keeps its own extra sequence as a tip piece (the
# generator never adds sequence that close to an end), and an assembly gap is kept as a piece of the
# haplotype's own N bases rather than the graph's sequence opposite it. Every walk then spells its
# haplotype exactly.
"$PY" - "$OUT/tg.gfa" <<'PYEOF'
import random, sys
r = random.Random(3)
rnd = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
pre, base, post = rnd(2000), rnd(100000), rnd(3000)
haps = {'ref': base, 'h1': base[:50000] + base[52000:], 'h2': base + post, 'h3': pre + base,
        'h4': base[:40000] + 'N' * 7400 + base[47400:], 'h5': base[:60000] + 'N' * 800 + base[60800:]}
with open(sys.argv[1], 'w') as o:
    o.write('H\tVN:Z:1.0\n')
    for i, (n, s) in enumerate(haps.items(), 1): o.write(f'S\t{i}\t{s}\n')
    for i, n in enumerate(haps, 1): o.write(f'P\t{n}\t{i}+\t*\n')
PYEOF
"$BIN" rebuild -i "$OUT/tg.gfa" -o "$OUT/tg.out.gfa" --force -r ref -t 2 > "$OUT/tg.log" 2>&1
A="$OUT/tg.out.rebuild.report.tsv"
v=$(awk -F'\t' '$1=="#verdict"{print $2}' "$A" 2>/dev/null)
exact=$(awk -F'\t' 'NR>1 && $0!~/^#/ && $7==1' "$A" 2>/dev/null | wc -l | tr -d ' ')
tips=$(awk -F'\t' '$1=="#tips_added"{print $2}' "$A" 2>/dev/null); gaps=$(awk -F'\t' '$1=="#gaps_kept"{print $2}' "$A" 2>/dev/null)
[ "$v" = "accepted" ] && [ "$exact" = "6" ] && [ "$tips" = "2" ] && [ "$gaps" = "2" ] \
  && ok "tips and assembly gaps: accepted, 6/6 walks exact, 2 tips, 2 gap pieces" \
  || bad "tips/gaps: verdict '$v', $exact exact walks, tips '$tips', gaps '$gaps'"
# --clip-to-reference: the same two extensions are dropped instead (3000 + 2000 bp), every walk exact.
"$BIN" rebuild -i "$OUT/tg.gfa" -o "$OUT/tgc.gfa" --force -r ref --clip-to-reference -t 2 > /dev/null 2>&1
cb=$(awk -F'\t' '$1=="#clipped_bp"{print $2}' "$OUT/tgc.rebuild.report.tsv" 2>/dev/null)
ex=$(awk -F'\t' 'NR>1 && $0!~/^#/ && $7==1' "$OUT/tgc.rebuild.report.tsv" 2>/dev/null | wc -l | tr -d ' ')
[ "$cb" = "5000" ] && [ "$ex" = "6" ] && ok "--clip-to-reference drops exactly the 5000 bp past the reference's ends" \
                                  || bad "--clip-to-reference: clipped '$cb' bp, $ex exact walks"
nn=$(awk -F'\t' '$1=="S" && $4=="mn:Z:gap" && $3 !~ /[^N]/' "$OUT/tg.out.gfa" | wc -l | tr -d ' ')
[ "$nn" = "2" ] && ok "  ... and each gap is a piece of N bases, not graph sequence" \
                || bad "  ... gap pieces holding only N: $nn, expected 2"

# A gap whose deletion the aligner places on the far side of a matching base must still replace the
# graph sequence opposite it (the walk spelled both the graph sequence and the N run), and a rotated
# haplotype -- its "extension" past the graph is graph sequence the alignment did not reach -- must not
# become a private tip that passes it.
"$PY" - "$OUT/gisl.gfa" "$OUT/rot.gfa" <<'PYEOF'
import random, sys
def write(path, haps):
    with open(path, 'w') as o:
        o.write('H\tVN:Z:1.0\n')
        for i, (n, s) in enumerate(haps.items(), 1): o.write(f'S\t{i}\t{s}\n')
        for i, n in enumerate(haps, 1): o.write(f'P\t{n}\t{i}+\t*\n')
r = random.Random(1)
rnd = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
b = list(rnd(100000)); b[41999] = b[39999]; b = ''.join(b)       # equal-cost deletion placement
write(sys.argv[1], {'ref': b, 'h1': b[:60000] + b[62000:], 'h4': b[:40000] + 'N' * 2000 + b[42000:]})
r = random.Random(21)
b = rnd(100000)
write(sys.argv[2], {'ref': b, 'h1': b[:30000] + b[32000:], 'rot': b[60000:] + b[:60000]})
PYEOF
"$BIN" rebuild -i "$OUT/gisl.gfa" -o "$OUT/gisl.out.gfa" --force -r ref -t 2 > /dev/null 2>&1
a=$(awk -F'\t' '$1=="h4"{print $7":"$NF}' "$OUT/gisl.out.rebuild.report.tsv" 2>/dev/null)
[ "$a" = "1:ok" ] && ok "a gap with a matching base inside its alignment block is replaced exactly" \
                 || bad "gap with a matching island: h4 gave '$a'"
"$BIN" rebuild -i "$OUT/rot.gfa" -o "$OUT/rot.out.gfa" --force -r ref -t 2 > /dev/null 2>&1
a=$(awk -F'\t' '$1=="rot"{print $NF}' "$OUT/rot.out.rebuild.report.tsv" 2>/dev/null)
tp=$(awk -F'\t' '$1=="#tips_added"{print $2}' "$OUT/rot.out.rebuild.report.tsv" 2>/dev/null)
[ "$a" != "ok" ] && [ "${tp:-0}" = "0" ] && ok "a rotated haplotype is not passed through a duplicate tip ($a)" \
                                     || bad "rotated haplotype: status '$a', tips '$tp'"

# ---------------------------------------------------------------- options checked before any work
# --audit is refused: a user-chosen report path could name the input or the output graph, and the report's
# final rename would replace that graph with a table.
e=$("$BIN" rebuild -i "$GFA" -o "$OUT/au.gfa" --audit "$GFA" 2>&1 | grep -c "audit was removed")
[ "$e" -gt 0 ] && ok "--audit is refused with a pointer to <prefix>.rebuild.report.tsv" \
               || bad "--audit was accepted"
[ "$(wc -c < "$GFA")" -gt 0 ] && head -1 "$GFA" | grep -q '^H' \
  && ok "  ... and the input graph is untouched" || bad "--audit damaged the input graph"

# A reference name that matches nothing is an error at startup, on the pass-through path too, not a
# pass-through with exit 0 after the whole build. Case does not matter, as elsewhere.
"$BIN" rebuild -i "$GFA" -o "$OUT/nr.gfa" -r "NOT_A_REAL_PATH_XYZ" > "$OUT/nr.log" 2>&1
rc=$?
[ "$rc" -ne 0 ] && grep -q "reference path not found" "$OUT/nr.log" \
  && ok "an unknown --reference-path fails at startup (exit $rc)" \
  || bad "an unknown --reference-path did not fail (exit $rc)"
[ ! -e "$OUT/nr.gfa" ] && ok "  ... before writing any output" || bad "an unknown reference still wrote output"
p1=$(awk -F'\t' '$1=="P"{print $2; exit}' "$GFA")
up=$(printf '%s' "$p1" | tr '[:lower:]' '[:upper:]')
if [ "$up" != "$p1" ]; then
  "$BIN" rebuild -i "$GFA" -o "$OUT/ci.gfa" -r "$up" >/dev/null 2>&1
  [ "$?" -eq 0 ] && ok "--reference-path is matched case-insensitively ($up -> $p1)" \
                 || bad "--reference-path $up did not resolve to $p1"
fi

# Sizes are non-negative integers that fit what they are passed to; anything else is refused rather
# than wrapped or truncated into a different number.
for c in "--min-var:-1" "--min-var:1e3" "--min-align-len:4294968296" "--kmer:0" "--kmer:5" "--threads:4294967296"; do
  f=${c%%:*}; v=${c#*:}
  "$BIN" rebuild -i "$GFA" -o "$OUT/sz.gfa" "$f" "$v" >/dev/null 2>&1
  [ "$?" -ne 0 ] && ok "refused: $f $v" || bad "$f $v was accepted"
done

# ---------------------------------------------------------------- output codec follows -o
is_gz() { [ "$(head -c 2 "$1" | od -An -tx1 | tr -d ' \n')" = "1f8b" ]; }
gzip -c "$GFA" > "$OUT/in.gfa.gz"
"$BIN" rebuild -i "$OUT/in.gfa.gz" -o "$OUT/plain.gfa" >/dev/null 2>&1
if is_gz "$OUT/plain.gfa"; then
  bad "a gzipped input passed through into a .gfa output is still gzip"
elif cmp -s "$GFA" "$OUT/plain.gfa"; then
  ok "a gzipped input passes through into a .gfa output as plain text, same content"
else
  pin=$(grep -c '^P' "$GFA"); pout=$(grep -c '^P' "$OUT/plain.gfa" 2>/dev/null || echo 0)
  [ "$pin" = "$pout" ] && ok "a gzipped input is written as plain text (graph was rebuilt)" \
                       || bad "plain output from a gzipped input lost paths"
fi
"$BIN" rebuild -i "$GFA" -o "$OUT/z.gfa.gz" >/dev/null 2>&1
is_gz "$OUT/z.gfa.gz" \
  && ok "a .gz output is gzip-compressed" || bad "a .gz output is not gzip"
gzip -dc "$OUT/z.gfa.gz" 2>/dev/null | cmp -s - "$OUT/pass.gfa" \
  && ok "  ... and decompresses to the same graph as the plain run" || bad ".gz output content differs"

# ---------------------------------------------------------------- no stale report beside a pass-through
printf '#verdict\taccepted\n' > "$OUT/st.rebuild.report.tsv"
"$BIN" rebuild -i "$GFA" -o "$OUT/st.gfa" >/dev/null 2>&1
if cmp -s "$GFA" "$OUT/st.gfa"; then
  [ ! -e "$OUT/st.rebuild.report.tsv" ] && ok "a pass-through removes a report left by an earlier run" \
                                          || bad "a stale report survived beside a passed-through graph"
fi

# A degenerate rebuild is reported as such, not as a healthy graph passed through.
if grep -q "holds no variation" "$OUT/flat.log" 2>/dev/null; then
  grep -q "rebuilds a healthy graph" "$OUT/flat.log" \
    && bad "a degenerate rebuild is reported as a healthy graph passed through" \
    || ok "a degenerate rebuild is reported as such, not as a healthy graph"
fi

if [ -n "$SRC" ] && [ -s "$SRC" ]; then
  # Soft-masked input: minigraph stores uppercase, so lowercase bases could read as mismatches and refuse
  # the rebuild. The same graph lowercased must rebuild the same way: same per-path recovery and
  # verdict. (The graphs themselves are compared only if accepted -- a rejected run writes its input,
  # which differs in case by construction.)
  awk -F'\t' 'BEGIN{OFS="\t"} $1=="S"{$3=tolower($3)} {print}' "$GFA" > "$OUT/lc.gfa"
  "$BIN" rebuild -i "$OUT/lc.gfa" -o "$OUT/lc.out.gfa" --force --min-recovered-identity 0.97 -t 1 -q \
    >/dev/null 2>&1
  aud_core() { awk -F'\t' '$1 ~ /^#verdict/ || ($0 !~ /^#/ && NR>1){print $1, $3, $5, $7, $NF}' "$1"; }
  if cmp -s <(aud_core "$OUT/lc.out.rebuild.report.tsv") <(aud_core "$OUT/th1.rebuild.report.tsv") \
     && { ! grep -q '^#verdict.accepted$' "$OUT/th1.rebuild.report.tsv" || cmp -s "$OUT/lc.out.gfa" "$OUT/th1.gfa"; }; then
    ok "a soft-masked (lowercase) graph rebuilds exactly as the uppercase one"
  else
    bad "lowercase input changes the rebuild"
  fi

  # With -r, the reference seeds the graph and must come back byte for byte. --base-align, so that every
  # indel of the test graph is found and the rebuild is accepted; the checks below need an accepted one.
  rp=$(awk -F'\t' '$1=="P"{n++; if(n==3){print $2; exit}}' "$GFA")
  "$BIN" rebuild -i "$GFA" -o "$OUT/rs.gfa" --force --base-align -t 4 -r "$rp" \
    > "$OUT/rs.log" 2>&1
  grep -q '^#verdict.accepted$' "$OUT/rs.rebuild.report.tsv" \
    && ok "a forced --base-align rebuild with -r is accepted" \
    || bad "the forced --base-align rebuild with -r was not accepted: $(awk -F'\t' '$1=="#reason"{print $2}' "$OUT/rs.rebuild.report.tsv")"

  # Segment names from minigraph go in a user tag, not rGFA's SN.
  ! grep -q $'\tSN:Z:' "$OUT/rs.gfa" && grep -q $'\tmn:Z:' "$OUT/rs.gfa" \
    && ok "rebuilt segments keep minigraph's name as mn:Z, not rGFA SN:Z" \
    || bad "rebuilt segments still carry SN:Z (or no name tag)"
  grep -q "seed: $rp" "$OUT/rs.log" && ok "-r seeds the graph from the reference ($rp)" \
                                   || bad "-r did not seed from $rp: $(grep -o 'seed: [^ ]*' "$OUT/rs.log" | head -1)"
  rid=$(awk -F'\t' -v p="$rp" '$1==p{print $7}' "$OUT/rs.rebuild.report.tsv" 2>/dev/null)
  [ "$rid" = "1" ] && ok "  ... and its walk spells it exactly (walk_identity $rid)" \
                   || bad "  ... but its walk identity is '$rid', not 1"

  # An interrupted rebuild removes its scratch FASTAs and staged files.
  mkdir -p "$OUT/sig"
  "$BIN" rebuild -i "$GFA" -o "$OUT/sig/s.gfa" --force -t 1 -q >/dev/null 2>&1 &
  pid=$!
  for _ in $(seq 1 600); do
    [ -n "$(find "$OUT/sig" -name '*.fa' 2>/dev/null | head -1)" ] && break
    kill -0 "$pid" 2>/dev/null || break
    sleep 0.1
  done
  kill -TERM "$pid" 2>/dev/null
  wait "$pid"; rc=$?
  left=$(find "$OUT/sig" -name '*rebuild-tmp*' -o -name '*rebuild.tmp*' | wc -l | tr -d ' ')
  [ "$rc" -eq 143 ] && [ "$left" = "0" ] && ok "SIGTERM mid-rebuild leaves no scratch or staged files" \
                    || bad "SIGTERM mid-rebuild: exit $rc, $left temporary artefact(s) left"

  # A healthy graph is passed through: its staged copy is removed on an interrupt too. The window is
  # the copy itself, so the test waits for the staged file rather than sleeping; if the copy finishes
  # first, it is skipped, not passed.
  mkdir -p "$OUT/sigp"
  "$BIN" rebuild -i "$GFA" -o "$OUT/sigp/p.gfa.gz" -q >/dev/null 2>&1 &
  pid=$!
  seen=0
  for _ in $(seq 1 3000); do
    [ -n "$(find "$OUT/sigp" -name '*rebuild-tmp*' 2>/dev/null | head -1)" ] && { seen=1; break; }
    kill -0 "$pid" 2>/dev/null || break
    sleep 0.01
  done
  kill -TERM "$pid" 2>/dev/null
  wait "$pid"; rc=$?
  left=$(find "$OUT/sigp" -name '*rebuild-tmp*' | wc -l | tr -d ' ')
  if [ "$seen" = "1" ] && [ "$rc" -eq 143 ]; then
    [ "$left" = "0" ] && ok "SIGTERM during a pass-through copy leaves no staged file" \
                      || bad "SIGTERM during a pass-through copy left $left staged file(s)"
  else
    printf "  skip SIGTERM during a pass-through copy: the copy finished first (exit %s)\n" "$rc"
  fi

  # A gzip output is written in chunks (gzwrite takes an unsigned length); with a 4 kb chunk the
  # pass-through copy is many chunks and decompresses to the input byte for byte.
  chunks=$(PANVAR_TEST_GZ_CHUNK=4096 "$BIN" rebuild -i "$GFA" -o "$OUT/chunk.gfa.gz" -q 2>&1 >/dev/null |
           sed -n 's/^rebuild: gzip output written in \([0-9]*\) chunk(s)$/\1/p')
  bytes=$( (gzip -dc < "$GFA" 2>/dev/null || cat "$GFA") | wc -c | tr -d ' ')
  if [ "$chunks" = "$(( (bytes + 4095) / 4096 ))" ] &&
     cmp -s <(gzip -dc < "$GFA" 2>/dev/null || cat "$GFA") <(gzip -dc < "$OUT/chunk.gfa.gz"); then
    ok "a gzip output written in $chunks 4 kb chunks is the input byte for byte"
  else
    bad "a gzip output in 4 kb chunks: ${chunks:-no} chunk(s) for $bytes bytes, or it differs from the input"
  fi
fi

# ---------------------------------------------------------------- clipping, names, short junk ends, IUPAC, options
"$PY" - "$OUT/r2" <<'PYEOF'
import random, sys
out = sys.argv[1]
r = random.Random(23)
b = ''.join(r.choice('ACGT') for _ in range(100000))
flank = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
def write(name, seqs, pnames=None, extra=''):
    with open(f'{out}.{name}.gfa', 'w') as o:
        o.write('H\tVN:Z:1.0\n')
        for i, (n, s) in enumerate(seqs.items(), 1):
            o.write(f'S\t{i}\t{s}\n')
        for i, n in enumerate(seqs, 1):
            o.write(f'P\t{(pnames or {}).get(n, n)}\t{i}+\t*\n')
        o.write(extra)
# A haplotype rotated by 60 kb: its last 40 kb are locus sequence, mapped by a second chain.
write('rot', {'ref': b, 'h1': b[:50000] + b[52000:], 'rot': b[60000:] + b[:60000]})
# Haplotypes running past the reference's ends into flank, named with PanSN intervals.
write('name', {'ref': b, 'h1': b[:50000] + b[52000:], 'S2': b + flank(3000), 'S3': flank(2000) + b},
      {'S2': 'S2#1#ctg2:1001-104000', 'S3': 'S3#1#ctg3:5001-107000'})
# A truncated haplotype ending in 30 bp that match nothing: under --min-var, so not a reason to refuse.
write('junk', {'ref': b, 'h1': b[:50000] + b[52000:], 'h2': b[:97000] + flank(30)})
# The reference walks a segment carrying IUPAC codes in reverse; its complement must come back.
amb = 'R' * 40
mid = b[40000:40200]
with open(f'{out}.iupac.gfa', 'w') as o:
    rc = lambda s: s[::-1].translate(str.maketrans('ACGTRY', 'TGCAYR'))
    o.write('H\tVN:Z:1.0\n')
    o.write(f'S\t1\t{b[:40000]}\nS\t2\t{rc(mid + amb)}\nS\t3\t{b[40000:]}\nS\t4\t{b[:40000]}\nS\t5\t{b[40200:]}\n')
    o.write('L\t1\t+\t2\t-\t0M\nL\t2\t-\t3\t+\t0M\nL\t4\t+\t5\t+\t0M\n')
    o.write('P\tref\t1+,2-,3+\t*\nP\tdel\t4+,5+\t*\n')
PYEOF
aud() { awk -F'\t' -v k="$2" '$1==k{print $'"${3:-2}"'}' "${1%.gfa}.rebuild.report.tsv" 2>/dev/null; }

"$BIN" rebuild -i "$OUT/r2.rot.gfa" -o "$OUT/r2.rot.out.gfa" --force -r ref --clip-to-reference -t 2 >/dev/null 2>&1
[ "$(aud "$OUT/r2.rot.out.gfa" '#verdict')" = "rejected" ] \
  && ok "--clip-to-reference does not clip locus sequence another chain maps (rotated haplotype refused)" \
  || bad "--clip-to-reference clipped a rotated haplotype's locus sequence (verdict $(aud "$OUT/r2.rot.out.gfa" '#verdict'))"

"$BIN" rebuild -i "$OUT/r2.name.gfa" -o "$OUT/r2.name.out.gfa" --force -r ref --clip-to-reference -t 2 >/dev/null 2>&1
names=$(awk -F'\t' '$1=="P" && $2 ~ /ctg/ {print $2}' "$OUT/r2.name.out.gfa" | sort | paste -sd' ' -)
[ "$names" = "S2#1#ctg2:1001-101000 S3#1#ctg3:7001-107000" ] \
  && ok "a clipped haplotype's PanSN interval is narrowed by what was clipped" \
  || bad "clipped haplotype names: '$names', expected S2#1#ctg2:1001-101000 S3#1#ctg3:7001-107000"

"$BIN" rebuild -i "$OUT/r2.junk.gfa" -o "$OUT/r2.junk.out.gfa" --force -t 2 >/dev/null 2>&1
[ "$(aud "$OUT/r2.junk.out.gfa" '#verdict')" = "accepted" ] && [ "$(aud "$OUT/r2.junk.out.gfa" h2 NF)" = "ok" ] \
  && ok "a truncated haplotype ending in under --min-var junk bases is cut, not a reason to refuse" \
  || bad "junk tail: verdict '$(aud "$OUT/r2.junk.out.gfa" '#verdict')', h2 '$(aud "$OUT/r2.junk.out.gfa" h2 NF)'"

"$BIN" rebuild -i "$OUT/r2.iupac.gfa" -o "$OUT/r2.iupac.out.gfa" --force -r ref -t 2 >/dev/null 2>&1
spell() { "$PY" - "$1" "$2" <<'PYEOF'
import sys
seq, path = {}, None
for l in open(sys.argv[1]):
    f = l.rstrip('\n').split('\t')
    if f[0] == 'S': seq[f[1]] = f[2]
    if f[0] == 'P' and f[1] == sys.argv[2]: path = f[2]
rc = lambda s: s[::-1].translate(str.maketrans('ACGTRYKMBVDHSWNacgtrykmbvdhswn', 'TGCAYRMKVBHDSWNtgcayrmkvbhdswn'))
print(''.join(seq[s[:-1]] if s[-1] == '+' else rc(seq[s[:-1]]) for s in path.split(',')))
PYEOF
}
[ "$(spell "$OUT/r2.iupac.gfa" ref)" = "$(spell "$OUT/r2.iupac.out.gfa" ref)" ] \
  && [ "$(aud "$OUT/r2.iupac.out.gfa" '#verdict')" = "accepted" ] \
  && ok "IUPAC codes on a reverse step come back complemented (Y), not as N" \
  || bad "IUPAC codes on a reverse step: the reference does not spell back exactly"

mkdir -p "$OUT/r2.dir.rebuild.report.tsv"
"$BIN" rebuild -i "$OUT/r2.junk.gfa" -o "$OUT/r2.dir.gfa" --force -q -t 2 >/dev/null 2>&1
rc=$?
[ "$rc" -ne 0 ] && [ ! -e "$OUT/r2.dir.gfa" ] \
  && ok "a report path that cannot be written is refused before anything is installed" \
  || bad "an unwritable report path: exit $rc, graph $( [ -e "$OUT/r2.dir.gfa" ] && echo installed || echo absent)"

"$BIN" rebuild -i "$OUT/r2.junk.gfa" -o "$OUT/r2.h0.gfa" --min-hubs 0 >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "--min-hubs 0 is refused (it forced a rebuild without --force's warning)" \
               || bad "--min-hubs 0 was accepted"

"$BIN" rebuild -i "$OUT/r2.junk.gfa" -o "$OUT/r2.td.gfa" --force -t 2 --tmp-dir "$OUT/td_new/deep" >/dev/null 2>&1
[ ! -e "$OUT/td_new" ] && ok "--tmp-dir parents the run created are removed with the scratch" \
                       || bad "--tmp-dir left $(find "$OUT/td_new" | head -3 | paste -sd' ' -)"

# ---------------------------------------------------------------- piped input, name narrowing, gaps beside divergence
"$PY" - "$OUT/r3" <<'PYEOF'
import random, sys
d = sys.argv[1]
def write(path, haps):
    with open(path, 'w') as o:
        o.write('H\tVN:Z:1.0\n')
        for i, (n, s) in enumerate(haps, 1): o.write(f'S\t{i}\t{s}\n')
        for i, (n, s) in enumerate(haps, 1): o.write(f'P\t{n}\t{i}+\t*\n')
r = random.Random(3)
rnd = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
def mutate(s, rate): return ''.join(r.choice([c for c in 'ACGT' if c != b]) if r.random() < rate else b for b in s)
pre, base, post = rnd(2000), rnd(100000), rnd(3000)
# half-open PanSN names, as the validated loci use them (end - start == length)
write(d + '.half.gfa', [('grch38#1#chr1:500000-600000', base), ('s1#1#ctg1:0-98000', base[:50000] + base[52000:]),
                        ('s2#1#ctg2:1000-106000', pre + base + post)])
# two names that narrow onto one another after clipping
left, right = rnd(5000), rnd(5000)
write(d + '.coll.gfa', [('ref#1#chr1:1-100000', base), ('sample#1#chr1:1-105000', left + base),
                        ('sample#1#chr1:5001-110000', base + right), ('variant#1#chr1:1-102000', base[:50000] + rnd(2000) + base[50000:])])
# 3% divergent beside a 4 kb gap the walk keeps as N
m = mutate(base, 0.03)
write(d + '.div.gfa', [('ref', base), ('h1', base[:50000] + base[52000:]), ('div', m[:30000] + 'N' * 4000 + m[34000:])])
# a rotation whose breakpoint is an N run
write(d + '.rotn.gfa', [('ref', base), ('h1', base[:30000] + base[32000:]), ('rot', base[60000:] + 'N' * 3000 + base[3000:60000])])
PYEOF
# A pipe gives its bytes once; the original is copied from the input when passed through, so a streamed
# input would write an empty graph.
cat "$OUT/r3.half.gfa" | "$BIN" rebuild -i /dev/stdin -o "$OUT/r3.stdin.gfa.gz" -t 2 -q >/dev/null 2>&1; rc=$?
[ "$rc" -ne 0 ] && [ ! -e "$OUT/r3.stdin.gfa.gz" ] && ok "a streamed input is refused (it cannot be read twice)" \
  || bad "streamed input: exit $rc, output $(ls -la "$OUT/r3.stdin.gfa.gz" 2>&1 | awk '{print $5}')"
# --clip-to-reference narrows half-open names as it does 1-based ones.
"$BIN" rebuild -i "$OUT/r3.half.gfa" -o "$OUT/r3.half.out.gfa" --force -r grch38 --clip-to-reference -t 2 -q >/dev/null 2>&1
[ "$(awk -F'\t' '$1=="P" && $2 ~ /^s2#/{print $2}' "$OUT/r3.half.out.gfa")" = "s2#1#ctg2:3000-103000" ] \
  && ok "a clipped half-open PanSN name is narrowed (ctg2:1000-106000 -> 3000-103000)" \
  || bad "clipped half-open name: $(awk -F'\t' '$1=="P"{print $2}' "$OUT/r3.half.out.gfa" | paste -sd' ' -)"
# Two names narrowed onto one: refused, --allow-loss or not (a graph with two P lines of one name is unreadable).
"$BIN" rebuild -i "$OUT/r3.coll.gfa" -o "$OUT/r3.coll.out.gfa" --force -r 'ref#1#chr1:1-100000' \
  --clip-to-reference --allow-loss -t 2 -q >/dev/null 2>&1
[ "$(aud "$OUT/r3.coll.out.gfa" '#verdict')" = "rejected" ] && grep -q "two haplotypes would be emitted as" "$OUT/r3.coll.out.rebuild.report.tsv" \
  && [ -z "$(awk -F'\t' '$1=="P"{print $2}' "$OUT/r3.coll.out.gfa" | sort | uniq -d)" ] \
  && ok "clipped names that collide are refused, even with --allow-loss" \
  || bad "name collision: verdict $(aud "$OUT/r3.coll.out.gfa" '#verdict')"
# N bases the walk spells as N cost no edits, so they do not excuse real differences.
"$BIN" rebuild -i "$OUT/r3.div.gfa" -o "$OUT/r3.div.out.gfa" --force -r ref -t 2 -q >/dev/null 2>&1
[ "$(awk -F'\t' '$1=="div"{print $NF}' "$OUT/r3.div.out.rebuild.report.tsv")" = "low_identity" ] \
  && ok "a divergent haplotype beside a kept gap reads low_identity, not assembly_gap" \
  || bad "divergent beside a gap: status $(awk -F'\t' '$1=="div"{print $NF}' "$OUT/r3.div.out.rebuild.report.tsv")"
# Two chains joined across a gap must follow each other along the graph: a rotation with an N run at its
# breakpoint must not be stitched back to back across a link spanning the locus.
"$BIN" rebuild -i "$OUT/r3.rotn.gfa" -o "$OUT/r3.rotn.out.gfa" --force -r ref -t 2 -q >/dev/null 2>&1
[ "$(awk -F'\t' '$1=="rot"{print $NF}' "$OUT/r3.rotn.out.rebuild.report.tsv")" != "ok" ] \
  && [ "$(aud "$OUT/r3.rotn.out.gfa" '#verdict')" = "rejected" ] \
  && ok "a rotated haplotype with an N run at its breakpoint is not joined" \
  || bad "rotation across a gap: rot $(awk -F'\t' '$1=="rot"{print $NF}' "$OUT/r3.rotn.out.rebuild.report.tsv"), verdict $(aud "$OUT/r3.rotn.out.gfa" '#verdict')"

echo
if [ "$fails" -eq 0 ]; then echo "rebuild_stats: all assertions passed"; exit 0; fi
echo "rebuild_stats: $fails assertion(s) FAILED"; exit 1
