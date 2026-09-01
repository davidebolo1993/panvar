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
  # Match the override's OWN wording, not the bare word WARNING: low matched cover also warns, so
  # a bare grep would pass on that alone and this assertion would survive --allow-loss going silent.
  grep -q "accepting a rebuild that fails the contract" "$OUT/loss.log" \
    && ok "--allow-loss accepts but warns" || bad "--allow-loss did not warn"

  # ---- matched cover is ADVISORY, and must not be able to reject on its own ----
  # It scores a recovery from landmarks -- seeds sampled along the haplotype, with any seed too
  # common to locate anything discarded before chaining. A repeat-dense locus therefore leaves much
  # of itself unlandmarked and the fraction falls with no sequence lost, so it cannot separate "a gap
  # was found here" from "nothing here could be checked". 1.0 is unsatisfiable by construction (a
  # rebuild is never byte-identical), so if cover could still reject, this would reject.
  "$BIN" rebuild -i "$GFA" -o "$OUT/cov.gfa" --force --min-matched-cover 1.0 \
    --min-recovered-identity 0.97 > "$OUT/cov.log" 2>&1
  grep -q "rejected" "$OUT/cov.log" \
    && bad "--min-matched-cover 1.0 rejected the rebuild; cover is meant to be advisory" \
    || ok "an unsatisfiable --min-matched-cover does not reject"
  grep -qi "min-matched-cover" "$OUT/cov.log" \
    && ok "  ... but it is still reported, so the signal is not silently dropped" \
    || bad "low matched cover was neither rejected nor reported -- the signal vanished"

  # The complement, so the pair above cannot both pass for a run that checks nothing: identity on the
  # same fixture still rejects. Without this, a build that ignored every bound would look correct.
  "$BIN" rebuild -i "$GFA" -o "$OUT/idb.gfa" --force --min-matched-cover 1.0 \
    --min-recovered-identity 0.999999 > "$OUT/idb.log" 2>&1
  grep -q "identity bound" "$OUT/idb.log" \
    && ok "identity still rejects, and the message names the bound that failed" \
    || bad "identity did not reject, or the message does not name the bound"

  # A path failing identity must SAY low_identity even when its cover is also low. The status test
  # used to run cover first, so every real failure was reported as the advisory reason instead.
  aud_i="$OUT/idb.gfa.rebuild_audit.tsv"
  if [ -s "$aud_i" ]; then
    n_li=$(awk -F'\t' 'NR>1 && $0 !~ /^#/ && $8=="low_identity"' "$aud_i" | wc -l | tr -d ' ')
    [ "$n_li" -gt 0 ] \
      && ok "a path failing identity reports low_identity, not the advisory low_cover" \
      || bad "no path reported low_identity under an unsatisfiable identity bound"
  fi

  # A reference that is not in the graph must reject, whatever else passes.
  "$BIN" rebuild -i "$GFA" -o "$OUT/badref.gfa" --force -r "NOT_A_REAL_PATH_XYZ" \
    > "$OUT/badref.log" 2>&1
  grep -q "reference path not found" "$OUT/badref.log" \
    && ok "a missing reference path is rejected" || bad "a missing reference path was accepted"

  # The audit sidecar must have one row per path and name a status for each.
  aud="$OUT/strict.gfa.rebuild_audit.tsv"
  if [ -s "$aud" ]; then
    # data rows only: the trailing #verdict/#reason lines are metadata, not paths
    rows=$(grep -vc '^#' "$aud"); rows=$(( rows - 1 ))
    [ "$rows" = "$in_paths" ] && ok "audit has one row per path ($rows)" \
                              || bad "audit has $rows rows for $in_paths paths"
    bad_status=$(awk -F'\t' 'NR>1 && $0 !~ /^#/ && $8!="ok" && $8!="not_recovered" && $8!="low_cover" && $8!="low_identity" && $8!="identity_unavailable" && $8!="dangling_walk"' "$aud" | wc -l | tr -d ' ')
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
  seed_of_k() { "$BIN" rebuild -i "$1" -o "$2" --force --min-recovered-identity 0.97 --kmer "$3" 2>&1 \
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

  # The audit must carry the global verdict beside the per-path rows.
  "$BIN" rebuild -i "$GFA" -o "$OUT/verd.gfa" --force --min-recovered-identity 0.999999 >/dev/null 2>&1
  grep -q "^#verdict" "$OUT/verd.gfa.rebuild_audit.tsv" 2>/dev/null \
    && ok "audit records the global verdict and reason" || bad "audit has no verdict row"
fi

# ---------------------------------------------------------------- equal richness is ordered, not arbitrary
# Two byte-identical haplotypes tie on every richness measure. The seed must be decided by name, not by
# whatever std::sort happens to do.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tACGTACGTACGTACGTACGT\nS\t2\tTTTTGGGGCCCCAAAATTTT\nS\t3\tGGGGCCCCAAAATTTTGGGG\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tzzz\t1+,2+,3+\t*\nP\taaa\t1+,2+,3+\t*\nP\tmmm\t1+,3+\t*\n'; } > "$OUT/tie.gfa"
t1=$("$BIN" rebuild -i "$OUT/tie.gfa" -o "$OUT/tie1.gfa" --force 2>&1 | grep -oE 'seed=[^;]*' | head -1)
t2=$("$BIN" rebuild -i "$OUT/tie.gfa" -o "$OUT/tie2.gfa" --force 2>&1 | grep -oE 'seed=[^;]*' | head -1)
[ -n "$t1" ] && [ "$t1" = "$t2" ] && ok "equal richness is broken deterministically ($t1)" \
                                  || bad "equal-richness seed is unstable: '$t1' vs '$t2'"
case "$t1" in *aaa*) ok "  ... and by path name, not input order" ;;
              *) bad "tie broken by something other than name: $t1" ;; esac


# ---------------------------------------------------------------- disposition, precision, units
if [ -n "$SRC" ] && [ -s "$SRC" ]; then
  # --allow-loss WRITES the rebuilt graph, so the audit must not call that "rejected".
  "$BIN" rebuild -i "$GFA" -o "$OUT/ovr.gfa" --force --min-recovered-identity 0.999999 --allow-loss \
    > /dev/null 2>&1
  v=$(awk -F'\t' '$1=="#verdict"{print $2}' "$OUT/ovr.gfa.rebuild_audit.tsv" 2>/dev/null)
  [ "$v" = "accepted_with_override" ] && ok "--allow-loss records disposition '$v'" \
                                      || bad "--allow-loss recorded verdict '$v', expected accepted_with_override"
  # The recorded threshold must be precise enough to reproduce the decision: 0.999999 must not print
  # as 1.0000.
  t=$(awk -F'\t' '$1=="#min_recovered_identity"{print $2}' "$OUT/ovr.gfa.rebuild_audit.tsv" 2>/dev/null)
  case "$t" in 0.999999*) ok "audit records the threshold at full precision ($t)" ;;
               *) bad "audit records threshold as '$t', losing the digits that decided it" ;; esac
fi

# An unchanged graph must not be called "untangled".
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tACGTACGTACGTACGTACGT\nS\t2\tTTTTGGGGCCCCAAAATTTT\n'
  printf 'L\t1\t+\t2\t+\t0M\n'
  printf 'P\tp1\t1+,2+\t*\nP\tp2\t1+,2+\t*\n'; } > "$OUT/flat.gfa"
"$BIN" rebuild -i "$OUT/flat.gfa" -o "$OUT/flat.out.gfa" --force > "$OUT/flat.log" 2>&1
if grep -q "structure:" "$OUT/flat.log"; then
  grep -q "untangled)" "$OUT/flat.log" && ! grep -q "unchanged\|NOT untangled" "$OUT/flat.log" \
    && bad "a graph with unchanged degree metrics was reported as untangled" \
    || ok "an unimproved graph is not reported as untangled"
else
  ok "an unimproved graph is not reported as untangled (no rebuild attempted)"
fi

# An unknown ('*') overlap means unverified, not zero.
printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTAC\nS\t2\tTTTTGGGGCC\nL\t1\t+\t2\t+\t*\nP\tp1\t1+,2+\t*\n' > "$OUT/e_star.gfa"
e=$("$BIN" rebuild -i "$OUT/e_star.gfa" -o "$OUT/e_star.out.gfa" 2>&1 | grep -c "UNKNOWN overlap")
[ "$e" -gt 0 ] && ok "an unknown ('*') overlap is refused, not read as 0M" \
               || bad "an unknown ('*') overlap was accepted as zero"


echo
if [ "$fails" -eq 0 ]; then echo "rebuild_stats: all assertions passed"; exit 0; fi
echo "rebuild_stats: $fails assertion(s) FAILED"; exit 1
