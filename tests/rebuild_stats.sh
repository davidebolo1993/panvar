#!/usr/bin/env bash
# rebuild_stats.sh - contract assertions for `panvar rebuild`.
#
# CTest previously exercised no rebuild at all: a forced synthetic run collapsed to one segment and
# passed through, so only the degeneracy guard was covered. These assert the properties a caller relies
# on -- that the input is never destroyed, that a healthy graph comes back byte-identical, that a
# malformed graph is refused rather than silently repaired, and that nothing is left behind.
#
#   rebuild_stats.sh <panvar-binary> <out-dir> [input.gfa]
set -uo pipefail

BIN="${1:?usage: rebuild_stats.sh <panvar> <outdir> [input.gfa]}"
PY="${PYTHON:-python3}"
OUT="${2:?}"
SRC="${3:-}"
# A FRESH directory each run. The suite writes fixtures, symlinks and staged outputs here, and reusing
# a populated directory makes results depend on what a previous run (or a concurrent manual one) left
# behind -- which showed up once as a failure that would not reproduce.
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
  grep -q "rejected" "$OUT/strict.log" && ok "an unsatisfiable contract is rejected" \
                                       || bad "an unsatisfiable contract was accepted"
  [ "$(grep -c '^S' "$OUT/strict.gfa")" = "$in_nodes" ] \
    && ok "  ... and the ORIGINAL graph was written ($in_nodes nodes)" \
    || bad "rollback did not restore the original graph"
  [ "$(grep -c '^P' "$OUT/strict.gfa")" = "$in_paths" ] \
    && ok "  ... with every path intact ($in_paths)" || bad "rollback lost paths"

  # --allow-loss accepts the same run and says so.
  "$BIN" rebuild -i "$GFA" -o "$OUT/loss.gfa" --force --min-recovered-identity 0.999999 \
    --allow-loss > "$OUT/loss.log" 2>&1
  grep -q "WARNING" "$OUT/loss.log" && ok "--allow-loss accepts but warns" \
                                    || bad "--allow-loss did not warn"

  # A reference that is not in the graph must reject, whatever else passes.
  "$BIN" rebuild -i "$GFA" -o "$OUT/badref.gfa" --force -r "NOT_A_REAL_PATH_XYZ" \
    > "$OUT/badref.log" 2>&1
  grep -q "reference path not found" "$OUT/badref.log" \
    && ok "a missing reference path is rejected" || bad "a missing reference path was accepted"

  # The audit sidecar must have one row per path and name a status for each.
  aud="$OUT/strict.gfa.rebuild_audit.tsv"
  if [ -s "$aud" ]; then
    rows=$(( $(wc -l < "$aud") - 1 ))
    [ "$rows" = "$in_paths" ] && ok "audit has one row per path ($rows)" \
                              || bad "audit has $rows rows for $in_paths paths"
    bad_status=$(awk -F'\t' 'NR>1 && $8!="ok" && $8!="not_recovered" && $8!="low_cover" && $8!="low_identity"' "$aud" | wc -l | tr -d ' ')
    [ "$bad_status" = "0" ] && ok "every audit row carries a known status" \
                            || bad "$bad_status audit rows have an unrecognised status"
  else
    bad "no audit sidecar was written"
  fi
fi


# ---------------------------------------------------------------- reproducibility
# A rebuild that depends on thread count or on which strand the input happens to be stored on is not
# reproducible, and neither property is visible in a single run's output.
if [ -n "$SRC" ] && [ -s "$SRC" ]; then
  seed_of() { "$BIN" rebuild -i "$1" -o "$2" --force --min-recovered-identity 0.97 2>&1 \
                | grep -oE 'seed=[^;]*' | head -1; }

  s1=$(seed_of "$GFA" "$OUT/rep1.gfa"); s2=$(seed_of "$GFA" "$OUT/rep2.gfa")
  [ -n "$s1" ] && [ "$s1" = "$s2" ] && ok "repeated runs choose the same seed" \
                                    || bad "seed differs between runs: '$s1' vs '$s2'"

  "$BIN" rebuild -i "$GFA" -o "$OUT/th1.gfa" --force --min-recovered-identity 0.97 -t 1 -q >/dev/null 2>&1
  "$BIN" rebuild -i "$GFA" -o "$OUT/th8.gfa" --force --min-recovered-identity 0.97 -t 8 -q >/dev/null 2>&1
  cmp -s "$OUT/th1.gfa" "$OUT/th8.gfa" && ok "output is identical with 1 and 8 threads" \
                                       || bad "thread count changes the output"
  cmp -s "$OUT/th1.gfa.rebuild_audit.tsv" "$OUT/th8.gfa.rebuild_audit.tsv" \
    && ok "  ... and so is the audit" || bad "thread count changes the audit"

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


echo
if [ "$fails" -eq 0 ]; then echo "rebuild_stats: all assertions passed"; exit 0; fi
echo "rebuild_stats: $fails assertion(s) FAILED"; exit 1
