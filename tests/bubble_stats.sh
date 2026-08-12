#!/usr/bin/env bash
# bubble_stats.sh - contract assertions for `panvar bubble`.
#
# CTest previously exercised the default internal cactus finder nowhere: synthetic_smoke.sh supplies
# checked-in external snarls, and real_smoke.sh makes its bubble assertions conditional on finding a
# first bubble -- and a header-only CSV is non-empty, so a total discovery regression would escape.
#
# Every graph here is small enough to reason about by hand, so each expected number is derivable rather
# than recorded from a previous run.
#
#   bubble_stats.sh <panvar-binary> <out-dir>
set -uo pipefail

BIN="${1:?usage: bubble_stats.sh <panvar> <outdir>}"
OUT="${2:?}"
mkdir -p "$OUT"
OUT="$OUT/run.$$.$(date +%s)"
rm -rf "$OUT"; mkdir -p "$OUT"
fails=0
ok()  { printf "  ok   %s\n" "$1"; }
bad() { printf "  FAIL %s\n" "$1"; fails=$((fails + 1)); }
# value of <column> in the first data row
cel() { awk -F, -v want="$2" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} NR==2{print $(c[want])}' "$1"; }
rows() { tail -n +2 "$1" 2>/dev/null | wc -l | tr -d ' '; }

# ---------------------------------------------------------------- a pure deletion is an allele
# Paths 1,2,3 and 1,3 over a 3-node graph: two haplotypes traverse the site, one carrying interior node
# 2 (4 bp) and one deleting it. So path_support is 2 and the SHORTEST interior span is 0, not 4.
{ printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTAC\nS\t2\tTTTT\nS\t3\tGGGGCCCCAA\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tfull\t1+,2+,3+\t*\nP\tdel\t1+,3+\t*\n'; } > "$OUT/del.gfa"
"$BIN" bubble -i "$OUT/del.gfa" -r full -o "$OUT/del" --min-variant-bp 0 -q >/dev/null 2>&1
[ "$(rows "$OUT/del.bubbles.csv")" = "1" ] && ok "one bubble found on the deletion graph" \
                                           || bad "expected 1 bubble, got $(rows "$OUT/del.bubbles.csv")"
[ "$(cel "$OUT/del.bubbles.csv" path_support)" = "2" ] \
  && ok "a pure-deletion allele counts as path support (2)" \
  || bad "path_support=$(cel "$OUT/del.bubbles.csv" path_support), expected 2"
[ "$(cel "$OUT/del.bubbles.csv" min_inside_bp)" = "0" ] \
  && ok "min_inside_bp is 0 -- the deletion, not the shortest insertion" \
  || bad "min_inside_bp=$(cel "$OUT/del.bubbles.csv" min_inside_bp), expected 0"
"$BIN" bubble -i "$OUT/del.gfa" -r full -o "$OUT/del2" --min-variant-bp 0 --min-path-support 2 -q >/dev/null 2>&1
[ "$(rows "$OUT/del2.bubbles.csv")" = "1" ] && ok "--min-path-support 2 keeps it" \
                                            || bad "--min-path-support 2 dropped a 2-path bubble"

# allele columns: two distinct walks, the reference's and the deletion's, one path each
[ "$(cel "$OUT/del.bubbles.csv" distinct_alleles)" = "2" ] && ok "distinct_alleles = 2" \
  || bad "distinct_alleles=$(cel "$OUT/del.bubbles.csv" distinct_alleles), expected 2"
[ "$(cel "$OUT/del.bubbles.csv" ref_allele_support)" = "1" ] && ok "ref_allele_support = 1" \
  || bad "ref_allele_support=$(cel "$OUT/del.bubbles.csv" ref_allele_support), expected 1"
[ "$(cel "$OUT/del.bubbles.csv" alt_allele_support_max)" = "1" ] && ok "alt_allele_support_max = 1" \
  || bad "alt_allele_support_max=$(cel "$OUT/del.bubbles.csv" alt_allele_support_max), expected 1"

# ---------------------------------------------------------------- boundaries carry reference order
# The cactus pair is unordered, but every consumer reads source/sink as an interval in reference order:
# `call` anchors coordinates on the source and merging joins one sink to the next source. On the
# deletion graph the reference is 1,2,3, so the bubble must be reported source=1, sink=3.
[ "$(cel "$OUT/del.bubbles.csv" source)" = "1" ] && [ "$(cel "$OUT/del.bubbles.csv" sink)" = "3" ] \
  && ok "boundaries are in reference order (source=1, sink=3)" \
  || bad "source=$(cel "$OUT/del.bubbles.csv" source) sink=$(cel "$OUT/del.bubbles.csv" sink), expected 1 and 3"

# ---------------------------------------------------------------- --superbubbles tests the GRAPH
# Interior nodes 2 and 3 are joined both ways, so the interior contains a directed cycle -- but no
# stored path walks it, so every path-derived signal (self-loop, revisit, both-orientation use) says
# acyclic. A superbubble's interior must be a DAG regardless of what anyone walked.
{ printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTAC\nS\t2\tTTTTTTTT\nS\t3\tCCCCCCCC\nS\t4\tGGGGCCCCAA\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'L\t2\t+\t3\t+\t0M\nL\t3\t+\t2\t+\t0M\n'
  printf 'P\tpA\t1+,2+,4+\t*\nP\tpB\t1+,3+,4+\t*\n'; } > "$OUT/cyc.gfa"
"$BIN" bubble -i "$OUT/cyc.gfa" -r pA -o "$OUT/cyc_all" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" bubble -i "$OUT/cyc.gfa" -r pA -o "$OUT/cyc_sb"  --min-variant-bp 0 --superbubbles -q >/dev/null 2>&1
[ "$(rows "$OUT/cyc_all.bubbles.csv")" = "1" ] && ok "the cyclic site is found without --superbubbles" \
                                               || bad "cyclic site not found at all"
[ "$(rows "$OUT/cyc_sb.bubbles.csv")" = "0" ] \
  && ok "an interior cycle no path walks is excluded by --superbubbles" \
  || bad "--superbubbles kept a site whose interior contains a directed cycle"
# ...and an acyclic site must survive, or the test is only measuring over-rejection.
"$BIN" bubble -i "$OUT/del.gfa" -r full -o "$OUT/del_sb" --min-variant-bp 0 --superbubbles -q >/dev/null 2>&1
[ "$(rows "$OUT/del_sb.bubbles.csv")" = "1" ] && ok "an acyclic site survives --superbubbles" \
                                              || bad "--superbubbles dropped an acyclic site"

# ---------------------------------------------------------------- filters apply after merging
# Two 1 bp bubbles 10 bp apart. Each satisfies --max-variant-bp 1; the fusion spans 12 bp and does not.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tC\nS\t3\tG\nS\t4\tTTTTTTTTTT\nS\t5\tA\nS\t6\tT\nS\t7\tCCCCCCCCCC\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'L\t4\t+\t5\t+\t0M\nL\t4\t+\t6\t+\t0M\nL\t5\t+\t7\t+\t0M\nL\t6\t+\t7\t+\t0M\n'
  printf 'P\tpA\t1+,2+,4+,5+,7+\t*\nP\tpB\t1+,3+,4+,6+,7+\t*\n'; } > "$OUT/mg.gfa"
"$BIN" bubble -i "$OUT/mg.gfa" -r pA -o "$OUT/mg0" --min-variant-bp 0 --max-variant-bp 1 -q >/dev/null 2>&1
[ "$(rows "$OUT/mg0.bubbles.csv")" = "2" ] && ok "two 1 bp bubbles pass --max-variant-bp 1 unmerged" \
                                           || bad "expected 2 unmerged bubbles, got $(rows "$OUT/mg0.bubbles.csv")"
"$BIN" bubble -i "$OUT/mg.gfa" -r pA -o "$OUT/mg1" --min-variant-bp 0 --max-variant-bp 1 \
  --merge-nearby-bp 50 -q >/dev/null 2>&1
big=$(awk -F, 'NR>1{for(i=1;i<=NF;i++);} NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} $(c["max_inside_bp"])+0 > 1' "$OUT/mg1.bubbles.csv" | wc -l | tr -d ' ')
[ "$big" = "0" ] && ok "no merged bubble exceeds --max-variant-bp (filters re-applied)" \
                 || bad "$big merged bubble(s) exceed --max-variant-bp 1"

# ---------------------------------------------------------------- an empty result is a result
printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTACGTACGT\nP\tp1\t1+\t*\n' > "$OUT/single.gfa"
"$BIN" bubble -i "$OUT/single.gfa" -r p1 -o "$OUT/single" -q >/dev/null 2>&1
rc=$?
[ "$rc" -eq 0 ] && [ -s "$OUT/single.bubbles.csv" ] && [ "$(rows "$OUT/single.bubbles.csv")" = "0" ] \
  && ok "a graph with no snarl emits an empty table and exits 0" \
  || bad "single-node graph: exit $rc, $(rows "$OUT/single.bubbles.csv") rows"

# ---------------------------------------------------------------- alt support discriminates, traversal does not
# Three haplotypes: two share one alternate walk, one takes another. Traversal support is 3 for the
# bubble whatever the alternates look like; alt support separates them.
{ printf 'H\tVN:Z:1.0\nS\t1\tAAAAAAAAAA\nS\t2\tCCCC\nS\t3\tGGGG\nS\t4\tTTTTTTTTTT\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t3\t+\t4\t+\t0M\nL\t1\t+\t4\t+\t0M\n'
  printf 'P\tref\t1+,4+\t*\nP\tha\t1+,2+,4+\t*\nP\thb\t1+,2+,4+\t*\nP\thc\t1+,3+,4+\t*\n'; } > "$OUT/alt.gfa"
"$BIN" bubble -i "$OUT/alt.gfa" -r ref -o "$OUT/alt" --min-variant-bp 0 -q >/dev/null 2>&1
[ "$(cel "$OUT/alt.bubbles.csv" path_support)" = "4" ] && ok "traversal support counts all 4 paths" \
  || bad "path_support=$(cel "$OUT/alt.bubbles.csv" path_support), expected 4"
[ "$(cel "$OUT/alt.bubbles.csv" alt_allele_support_max)" = "2" ] \
  && ok "alt_allele_support_max = 2 (the shared alternate)" \
  || bad "alt_allele_support_max=$(cel "$OUT/alt.bubbles.csv" alt_allele_support_max), expected 2"
"$BIN" bubble -i "$OUT/alt.gfa" -r ref -o "$OUT/alt3" --min-variant-bp 0 --min-alt-support 3 -q >/dev/null 2>&1
[ "$(rows "$OUT/alt3.bubbles.csv")" = "0" ] \
  && ok "--min-alt-support 3 drops it, where --min-path-support 3 would not" \
  || bad "--min-alt-support 3 kept a bubble whose best alternate has 2"
"$BIN" bubble -i "$OUT/alt.gfa" -r ref -o "$OUT/altp3" --min-variant-bp 0 --min-path-support 3 -q >/dev/null 2>&1
[ "$(rows "$OUT/altp3.bubbles.csv")" = "1" ] \
  && ok "  ... confirming --min-path-support 3 keeps it (traversal, not allele)" \
  || bad "--min-path-support 3 unexpectedly dropped it"

echo
if [ "$fails" -eq 0 ]; then echo "bubble_stats: all assertions passed"; exit 0; fi
echo "bubble_stats: $fails assertion(s) FAILED"; exit 1
