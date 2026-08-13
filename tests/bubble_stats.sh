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

# ---------------------------------------------------------------- merging is defined on the reference
# Three bubbles along a reference: 1..4, then 4..7 SHARING boundary node 4 (gap 0), then 9..12 behind a
# 30 bp spacer. --merge-nearby-bp 0 means merging off, so the smallest meaningful threshold is 1: it must
# fuse the shared-boundary pair and nothing else. That is also the B11 assertion -- the old distance was
# seeded with the whole length of the starting boundary node (10 bp here), so two bubbles sharing a
# boundary had a nominal gap of 10 and a threshold of 1 could never join them.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tC\nS\t3\tG\nS\t4\tTTTTTTTTTT\nS\t5\tA\nS\t6\tT\nS\t7\tCCCCCCCCCC\n'
  printf 'S\t8\tGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\nS\t9\tAAAAAAAAAA\nS\t10\tA\nS\t11\tC\nS\t12\tTTTTTTTTTT\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'L\t4\t+\t5\t+\t0M\nL\t4\t+\t6\t+\t0M\nL\t5\t+\t7\t+\t0M\nL\t6\t+\t7\t+\t0M\n'
  printf 'L\t7\t+\t8\t+\t0M\nL\t8\t+\t9\t+\t0M\n'
  printf 'L\t9\t+\t10\t+\t0M\nL\t9\t+\t11\t+\t0M\nL\t10\t+\t12\t+\t0M\nL\t11\t+\t12\t+\t0M\n'
  printf 'P\tpA\t1+,2+,4+,5+,7+,8+,9+,10+,12+\t*\n'
  printf 'P\tpB\t1+,3+,4+,6+,7+,8+,9+,11+,12+\t*\n'; } > "$OUT/mr.gfa"
"$BIN" bubble -i "$OUT/mr.gfa" -r pA -o "$OUT/mr0" --min-variant-bp 0 -q >/dev/null 2>&1
[ "$(rows "$OUT/mr0.bubbles.csv")" = "3" ] && ok "three bubbles found along the reference" \
                                           || bad "expected 3 bubbles, got $(rows "$OUT/mr0.bubbles.csv")"
"$BIN" bubble -i "$OUT/mr.gfa" -r pA -o "$OUT/mrA" --min-variant-bp 0 --merge-nearby-bp 1 -q >/dev/null 2>&1
[ "$(rows "$OUT/mrA.bubbles.csv")" = "2" ] \
  && ok "a threshold of 1 fuses the pair sharing a boundary (gap 0, not the boundary's length)" \
  || bad "shared-boundary merge gave $(rows "$OUT/mrA.bubbles.csv") bubbles, expected 2"
"$BIN" bubble -i "$OUT/mr.gfa" -r pA -o "$OUT/mrB" --min-variant-bp 0 --merge-nearby-bp 50 -q >/dev/null 2>&1
[ "$(rows "$OUT/mrB.bubbles.csv")" = "1" ] \
  && ok "--merge-nearby-bp 50 spans the 30 bp connector and fuses all three" \
  || bad "50 bp merge gave $(rows "$OUT/mrB.bubbles.csv") bubbles, expected 1"
# A 20 bp threshold must NOT cross a 30 bp connector: the distance is the sequence strictly between the
# facing boundaries, so it cannot be shrunk by the length of a boundary node.
"$BIN" bubble -i "$OUT/mr.gfa" -r pA -o "$OUT/mrC" --min-variant-bp 0 --merge-nearby-bp 20 -q >/dev/null 2>&1
[ "$(rows "$OUT/mrC.bubbles.csv")" = "2" ] \
  && ok "a 20 bp threshold does not cross a 30 bp connector" \
  || bad "20 bp merge gave $(rows "$OUT/mrC.bubbles.csv") bubbles, expected 2"

# Every fused bubble must contain only nodes inside its own reference span.
"$BIN" bubble -i "$OUT/mr.gfa" -r pA -o "$OUT/mrD" --min-variant-bp 0 --merge-nearby-bp 50 -q >/dev/null 2>&1
src=$(cel "$OUT/mrD.bubbles.csv" source); snk=$(cel "$OUT/mrD.bubbles.csv" sink)
[ "$src" = "1" ] && [ "$snk" = "12" ] \
  && ok "the fusion spans source=1 to sink=12, the outermost boundaries" \
  || bad "fused bubble is $src..$snk, expected 1..12"


# ---- B5: the interior is the graph's, not the panel's -------------------------------------------
# Node 5 is a third allele between the same two boundaries that NEITHER stored path takes. Collecting
# the interior from path intervals alone made it invisible: an allele nobody in the cohort carries is
# still part of the site, and everything computed from `inside` inherited the omission.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGG\nS\t4\tTTTTTTTTTT\nS\t5\tACGTACGTAC\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'L\t1\t+\t5\t+\t0M\nL\t5\t+\t4\t+\t0M\n'
  printf 'P\tpA\t1+,2+,4+\t*\nP\tpB\t1+,3+,4+\t*\n'; } > "$OUT/uw.gfa"
"$BIN" bubble -i "$OUT/uw.gfa" -r pA -o "$OUT/uw" --min-variant-bp 0 -q >/dev/null 2>&1
[ "$(cel "$OUT/uw.bubbles.csv" inside_node_count)" = "3" ] \
  && ok "an allele no path walks is still part of the interior (3 inside, not 2)" \
  || bad "unwalked branch: inside_node_count $(cel "$OUT/uw.bubbles.csv" inside_node_count), expected 3"

# And the case B4 alone could not reach: the only directed cycle sits on that unwalked branch, so the
# acyclicity search had nothing to find it in. --superbubbles must exclude the site.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGG\nS\t4\tTTTTTTTTTT\nS\t5\tACGTACGTAC\nS\t6\tTGCATGCATG\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'L\t1\t+\t5\t+\t0M\nL\t5\t+\t6\t+\t0M\nL\t6\t+\t5\t+\t0M\nL\t6\t+\t4\t+\t0M\n'
  printf 'P\tpA\t1+,2+,4+\t*\nP\tpB\t1+,3+,4+\t*\n'; } > "$OUT/uwc.gfa"
"$BIN" bubble -i "$OUT/uwc.gfa" -r pA -o "$OUT/uwc0" --min-variant-bp 0 -q >/dev/null 2>&1
[ "$(rows "$OUT/uwc0.bubbles.csv")" = "1" ] && ok "the cyclic site is reported without --superbubbles" \
                                            || bad "expected 1 bubble, got $(rows "$OUT/uwc0.bubbles.csv")"
"$BIN" bubble -i "$OUT/uwc.gfa" -r pA -o "$OUT/uwc1" --min-variant-bp 0 --superbubbles -q >/dev/null 2>&1
[ "$(rows "$OUT/uwc1.bubbles.csv")" = "0" ] \
  && ok "--superbubbles excludes a cycle carried on a branch no path walks" \
  || bad "cyclic unwalked branch survived --superbubbles: $(rows "$OUT/uwc1.bubbles.csv") bubbles"

# ---- the reference name must be the RESOLVED one -------------------------------------------------
# sort_graph_reference accepts an exact name, a unique case-insensitive match or a unique substring,
# but everything downstream compares path names exactly. Keeping the user's spelling meant `-r FULL`
# for a path named `full` sorted correctly and then matched nothing: boundaries came back unoriented
# and reference allele support read 0, with no error anywhere.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGG\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tfull\t1+,2+,3+\t*\nP\tdel\t1+,3+\t*\n'; } > "$OUT/case.gfa"
"$BIN" bubble -i "$OUT/case.gfa" -r FULL -o "$OUT/case" --min-variant-bp 0 -q >/dev/null 2>&1
[ "$(cel "$OUT/case.bubbles.csv" source)" = "1" ] && [ "$(cel "$OUT/case.bubbles.csv" sink)" = "3" ] \
  && ok "a case-insensitive reference alias still orients the boundaries" \
  || bad "-r FULL gave $(cel "$OUT/case.bubbles.csv" source)..$(cel "$OUT/case.bubbles.csv" sink), expected 1..3"
[ "$(cel "$OUT/case.bubbles.csv" ref_allele_support)" = "1" ] \
  && ok "a case-insensitive reference alias still finds the reference allele" \
  || bad "-r FULL gave ref_allele_support $(cel "$OUT/case.bubbles.csv" ref_allele_support), expected 1"

# ---- the graph contract applies through BOTH doors ------------------------------------------------
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGG\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tfull\t1+,2+,3+\t*\nP\tdel\t1+,3+\t*\n'; } > "$OUT/bad.gfa"
printf '{"start": {"node_id": "1"}, "end": {"node_id": "3"}}\n' > "$OUT/snarls.jsonl"
rm -f "$OUT/ext".*
"$BIN" bubble -i "$OUT/bad.gfa" --snarls-in "$OUT/snarls.jsonl" -o "$OUT/ext" --min-variant-bp 0 -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "--snarls-in refuses a path step with no link behind it" \
               || bad "--snarls-in accepted a malformed graph"
[ "$(ls "$OUT/ext".* 2>/dev/null | wc -l | tr -d ' ')" = "0" ] \
  && ok "a refused --snarls-in run writes no output" \
  || bad "--snarls-in left $(ls "$OUT/ext".* 2>/dev/null | wc -l | tr -d ' ') file(s) behind"

# ---- nothing lands until the run succeeds ---------------------------------------------------------
# The sorted GFA is written before discovery finishes, so a later failure used to leave a
# complete-looking graph for the next command in a pipeline to consume.
rm -f "$OUT/stage".*
"$BIN" bubble -i "$OUT/bad.gfa" -r full -o "$OUT/stage" --min-variant-bp 0 -q >/dev/null 2>&1
[ "$(ls "$OUT/stage".* 2>/dev/null | wc -l | tr -d ' ')" = "0" ] \
  && ok "a failed run leaves no partial output, not even the sorted GFA" \
  || bad "failed run left: $(ls "$OUT/stage".* 2>/dev/null | tr '\n' ' ')"

# ---- an alternate connector belongs to a fused bubble ---------------------------------------------
# Node 5 is a connector branch no path walks, so it is in neither part's interior and not on the
# reference. Assembling a fusion from its parts plus the reference connector dropped it.
{ printf 'H\tVN:Z:1.0\n'
  for n in 1 2 3 4 5 6 7 8 9 10; do printf "S\t$n\tACGTACGTAC\n"; done
  printf 'L\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'L\t4\t+\t6\t+\t0M\nL\t6\t+\t7\t+\t0M\nL\t4\t+\t5\t+\t0M\nL\t5\t+\t7\t+\t0M\n'
  printf 'L\t7\t+\t8\t+\t0M\nL\t7\t+\t9\t+\t0M\nL\t8\t+\t10\t+\t0M\nL\t9\t+\t10\t+\t0M\n'
  printf 'P\tpA\t1+,2+,4+,6+,7+,8+,10+\t*\nP\tpB\t1+,3+,4+,6+,7+,9+,10+\t*\n'; } > "$OUT/conn.gfa"
"$BIN" bubble -i "$OUT/conn.gfa" -r pA -o "$OUT/conn" --min-variant-bp 0 --min-alt-support 1 \
       --merge-nearby-bp 50 -q >/dev/null 2>&1
case "$(cel "$OUT/conn.bubbles.csv" inside_nodes)" in
  *5*) ok "a fused bubble contains the alternate connector branch, not just the reference route" ;;
  *)   bad "fused interior $(cel "$OUT/conn.bubbles.csv" inside_nodes) omits the alternate connector 5" ;;
esac

# ---- a cycle through a boundary handle ------------------------------------------------------------
# The acyclicity search is induced over the strict interior, so a cycle running through a BOUNDARY is
# not in the graph it searches. It is caught anyway, because a boundary carrying a cycle stops being a
# cut vertex and the decomposition promotes it into the interior of the enclosing snarl. Pinned here
# so a change to the decomposition cannot quietly reopen the hole.
{ printf 'H\tVN:Z:1.0\n'
  for n in 1 2 3 4 5 6 7; do printf "S\t$n\tACGTACGTAC\n"; done
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\nL\t4\t+\t5\t+\t0M\nL\t5\t+\t6\t+\t0M\n'
  printf 'L\t2\t+\t7\t+\t0M\nL\t7\t+\t4\t+\t0M\nL\t3\t+\t2\t+\t0M\n'
  printf 'P\tpA\t1+,2+,3+,4+,5+,6+\t*\nP\tpB\t1+,2+,7+,4+,5+,6+\t*\n'; } > "$OUT/bcyc.gfa"
"$BIN" bubble -i "$OUT/bcyc.gfa" -r pA -o "$OUT/bcyc" --min-variant-bp 0 --superbubbles -q >/dev/null 2>&1
[ "$(rows "$OUT/bcyc.bubbles.csv")" = "0" ] \
  && ok "--superbubbles excludes a cycle that runs through a boundary handle" \
  || bad "boundary-handle cycle survived --superbubbles: $(rows "$OUT/bcyc.bubbles.csv") bubbles"

echo
if [ "$fails" -eq 0 ]; then echo "bubble_stats: all assertions passed"; exit 0; fi
echo "bubble_stats: $fails assertion(s) FAILED"; exit 1
