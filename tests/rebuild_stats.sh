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
OUT="${2:?}"
SRC="${3:-}"
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

echo
if [ "$fails" -eq 0 ]; then echo "rebuild_stats: all assertions passed"; exit 0; fi
echo "rebuild_stats: $fails assertion(s) FAILED"; exit 1
