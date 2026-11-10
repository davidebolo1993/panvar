#!/usr/bin/env bash
# bubble_stats.sh - contract assertions for `panvar bubble`.
#
# These run the default internal cactus finder unless a check gives --snarls. A header-only CSV is
# non-empty, so only exact counts catch a discovery regression, and a check expecting no site also needs
# the run to have succeeded.
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
# every row's value of <column>, sorted numerically, space-joined
cel_all() { awk -F, -v want="$2" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next} {print $(c[want])}' "$1" | sort -n | paste -sd' ' -; }

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

# ---------------------------------------------------------------- a loop inside a site
# Interior nodes 2 and 3 are joined both ways, so the interior contains a directed cycle no stored path
# walks. A site with a loop inside is a site like any other.
{ printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTAC\nS\t2\tTTTTTTTT\nS\t3\tCCCCCCCC\nS\t4\tGGGGCCCCAA\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'L\t2\t+\t3\t+\t0M\nL\t3\t+\t2\t+\t0M\n'
  printf 'P\tpA\t1+,2+,4+\t*\nP\tpB\t1+,3+,4+\t*\n'; } > "$OUT/cyc.gfa"
"$BIN" bubble -i "$OUT/cyc.gfa" -r pA -o "$OUT/cyc_all" --min-variant-bp 0 -q >/dev/null 2>&1
[ "$(rows "$OUT/cyc_all.bubbles.csv")" = "1" ] && ok "a site with a loop inside is found" \
                                               || bad "cyclic site not found at all"

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
[ $? -eq 0 ] && [ -s "$OUT/alt3.bubbles.csv" ] && [ "$(rows "$OUT/alt3.bubbles.csv")" = "0" ] \
  && ok "--min-alt-support 3 drops it, where --min-path-support 3 would not" \
  || bad "--min-alt-support 3 failed or kept a bubble whose best alternate has 2"
"$BIN" bubble -i "$OUT/alt.gfa" -r ref -o "$OUT/altp3" --min-variant-bp 0 --min-path-support 3 -q >/dev/null 2>&1
[ "$(rows "$OUT/altp3.bubbles.csv")" = "1" ] \
  && ok "  ... confirming --min-path-support 3 keeps it (traversal, not allele)" \
  || bad "--min-path-support 3 unexpectedly dropped it"

# ---------------------------------------------------------------- merging is defined on the reference
# Three bubbles along a reference: 1..4, then 4..7 SHARING boundary node 4 (10 bp), then 9..12 behind a
# 30 bp spacer. The distance is the reference sequence between the two interiors, facing boundaries
# included: 10 bp for the pair sharing node 4 (a shared boundary counts its own length), 10+30+10 = 50
# for the pair across the spacer. Measured strictly between the facing boundaries instead, a shared
# boundary of any length would read as 0, and sites hundreds of bp apart across long boundaries would fuse.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tC\nS\t3\tG\nS\t4\tTTTTTTTTTT\nS\t5\tA\nS\t6\tT\nS\t7\tCCCCCCCCCC\n'
  printf 'S\t8\tGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\nS\t9\tAAAAAAAAAA\nS\t10\tA\nS\t11\tC\nS\t12\tTTTTTTTTTT\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'L\t4\t+\t5\t+\t0M\nL\t4\t+\t6\t+\t0M\nL\t5\t+\t7\t+\t0M\nL\t6\t+\t7\t+\t0M\n'
  printf 'L\t7\t+\t8\t+\t0M\nL\t8\t+\t9\t+\t0M\n'
  printf 'L\t9\t+\t10\t+\t0M\nL\t9\t+\t11\t+\t0M\nL\t10\t+\t12\t+\t0M\nL\t11\t+\t12\t+\t0M\n'
  printf 'P\tpA\t1+,2+,4+,5+,7+,8+,9+,10+,12+\t*\n'
  printf 'P\tpB\t1+,3+,4+,6+,7+,8+,9+,11+,12+\t*\n'; } > "$OUT/mr.gfa"
"$BIN" bubble -i "$OUT/mr.gfa" -r pA -o "$OUT/mr0" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
[ "$(rows "$OUT/mr0.bubbles.csv")" = "3" ] && ok "three bubbles found along the reference (merging off)" \
                                           || bad "expected 3 bubbles, got $(rows "$OUT/mr0.bubbles.csv")"
# The default is 100 bp, so with no option the 30 bp connector is bridged too.
"$BIN" bubble -i "$OUT/mr.gfa" -r pA -o "$OUT/mrD" --min-variant-bp 0 -q >/dev/null 2>&1
[ "$(rows "$OUT/mrD.bubbles.csv")" = "1" ] && ok "by default (100 bp) all three fuse into one site" \
                                           || bad "default merging gave $(rows "$OUT/mrD.bubbles.csv") bubbles, expected 1"
"$BIN" bubble -i "$OUT/mr.gfa" -r pA -o "$OUT/mrA" --min-variant-bp 0 --merge-nearby-bp 1 -q >/dev/null 2>&1
[ "$(rows "$OUT/mrA.bubbles.csv")" = "3" ] \
  && ok "a threshold of 1 does not bridge a 10 bp shared boundary" \
  || bad "1 bp merge gave $(rows "$OUT/mrA.bubbles.csv") bubbles, expected 3"
"$BIN" bubble -i "$OUT/mr.gfa" -r pA -o "$OUT/mrA10" --min-variant-bp 0 --merge-nearby-bp 10 -q >/dev/null 2>&1
[ "$(rows "$OUT/mrA10.bubbles.csv")" = "2" ] \
  && ok "a threshold of 10 fuses the pair sharing the 10 bp boundary, and nothing else" \
  || bad "10 bp merge gave $(rows "$OUT/mrA10.bubbles.csv") bubbles, expected 2"
"$BIN" bubble -i "$OUT/mr.gfa" -r pA -o "$OUT/mrB" --min-variant-bp 0 --merge-nearby-bp 50 -q >/dev/null 2>&1
[ "$(rows "$OUT/mrB.bubbles.csv")" = "1" ] \
  && ok "--merge-nearby-bp 50 spans the 30 bp connector and fuses all three" \
  || bad "50 bp merge gave $(rows "$OUT/mrB.bubbles.csv") bubbles, expected 1"
# A 49 bp threshold must NOT cross the 30 bp connector: with the two 10 bp boundaries the interiors are
# 50 bp apart.
"$BIN" bubble -i "$OUT/mr.gfa" -r pA -o "$OUT/mrC" --min-variant-bp 0 --merge-nearby-bp 49 -q >/dev/null 2>&1
[ "$(rows "$OUT/mrC.bubbles.csv")" = "2" ] \
  && ok "a 49 bp threshold does not cross a 30 bp connector between 10 bp boundaries" \
  || bad "49 bp merge gave $(rows "$OUT/mrC.bubbles.csv") bubbles, expected 2"
# Long boundaries: two sites sharing nothing but 1 bp between a 300 bp and a 400 bp boundary are 701 bp
# apart, and the default (100) leaves them alone.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\nS\t2\tC\nS\t3\tG\nS\t4\t%s\nS\t5\tA\nS\t6\t%s\nS\t7\tA\nS\t8\tC\nS\t9\t%s\n' \
    "$(printf 'ACGT%.0s' $(seq 1 10))" "$(printf 'ACGGT%.0s' $(seq 1 60))" "$(printf 'TTGCA%.0s' $(seq 1 80))" "$(printf 'GATC%.0s' $(seq 1 10))"
  printf 'L\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'L\t4\t+\t5\t+\t0M\nL\t5\t+\t6\t+\t0M\n'
  printf 'L\t6\t+\t7\t+\t0M\nL\t6\t+\t8\t+\t0M\nL\t7\t+\t9\t+\t0M\nL\t8\t+\t9\t+\t0M\n'
  printf 'P\tpA\t1+,2+,4+,5+,6+,7+,9+\t*\n'
  printf 'P\tpB\t1+,3+,4+,5+,6+,8+,9+\t*\n'; } > "$OUT/mrL.gfa"
"$BIN" bubble -i "$OUT/mrL.gfa" -r pA -o "$OUT/mrL" --min-variant-bp 0 -q >/dev/null 2>&1
[ "$(rows "$OUT/mrL.bubbles.csv")" = "2" ] \
  && ok "the default does not fuse sites 1 bp apart between a 300 bp and a 400 bp boundary" \
  || bad "long-boundary pair gave $(rows "$OUT/mrL.bubbles.csv") bubbles by default, expected 2"

# Every fused bubble must contain only nodes inside its own reference span.
"$BIN" bubble -i "$OUT/mr.gfa" -r pA -o "$OUT/mrD" --min-variant-bp 0 --merge-nearby-bp 50 -q >/dev/null 2>&1
src=$(cel "$OUT/mrD.bubbles.csv" source); snk=$(cel "$OUT/mrD.bubbles.csv" sink)
[ "$src" = "1" ] && [ "$snk" = "12" ] \
  && ok "the fusion spans source=1 to sink=12, the outermost boundaries" \
  || bad "fused bubble is $src..$snk, expected 1..12"


# ---- B5: the interior is the graph's, not the panel's -------------------------------------------
# Node 5 is a third allele between the same two boundaries that NEITHER stored path takes. Collected
# from path intervals alone, the interior would miss it: an allele nobody in the cohort carries is
# still part of the site, and everything computed from `inside` depends on it.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGG\nS\t4\tTTTTTTTTTT\nS\t5\tACGTACGTAC\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'L\t1\t+\t5\t+\t0M\nL\t5\t+\t4\t+\t0M\n'
  printf 'P\tpA\t1+,2+,4+\t*\nP\tpB\t1+,3+,4+\t*\n'; } > "$OUT/uw.gfa"
"$BIN" bubble -i "$OUT/uw.gfa" -r pA -o "$OUT/uw" --min-variant-bp 0 -q >/dev/null 2>&1
[ "$(cel "$OUT/uw.bubbles.csv" inside_node_count)" = "3" ] \
  && ok "an allele no path walks is still part of the interior (3 inside, not 2)" \
  || bad "unwalked branch: inside_node_count $(cel "$OUT/uw.bubbles.csv" inside_node_count), expected 3"

# The only directed cycle sits on that unwalked branch; the site is still reported.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGG\nS\t4\tTTTTTTTTTT\nS\t5\tACGTACGTAC\nS\t6\tTGCATGCATG\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'L\t1\t+\t5\t+\t0M\nL\t5\t+\t6\t+\t0M\nL\t6\t+\t5\t+\t0M\nL\t6\t+\t4\t+\t0M\n'
  printf 'P\tpA\t1+,2+,4+\t*\nP\tpB\t1+,3+,4+\t*\n'; } > "$OUT/uwc.gfa"
"$BIN" bubble -i "$OUT/uwc.gfa" -r pA -o "$OUT/uwc0" --min-variant-bp 0 -q >/dev/null 2>&1
[ "$(rows "$OUT/uwc0.bubbles.csv")" = "1" ] && ok "a site with a loop on an unwalked branch is reported" \
                                            || bad "expected 1 bubble, got $(rows "$OUT/uwc0.bubbles.csv")"

# ---- the reference name must be the RESOLVED one -------------------------------------------------
# sort_graph_reference accepts an exact name, a unique case-insensitive match or a unique substring,
# but everything downstream compares path names exactly. Kept as the user spelled it, `-r FULL` for a
# path named `full` would sort correctly and then match nothing: boundaries unoriented and reference
# allele support 0, with no error anywhere.
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
"$BIN" bubble -i "$OUT/bad.gfa" --snarls "$OUT/snarls.jsonl" -o "$OUT/ext" --min-variant-bp 0 -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "--snarls refuses a path step with no link behind it" \
               || bad "--snarls accepted a malformed graph"
[ "$(ls "$OUT/ext".* 2>/dev/null | wc -l | tr -d ' ')" = "0" ] \
  && ok "a refused --snarls run writes no output" \
  || bad "--snarls left $(ls "$OUT/ext".* 2>/dev/null | wc -l | tr -d ' ') file(s) behind"

# ---- nothing lands until the run succeeds ---------------------------------------------------------
# The sorted GFA is written before discovery finishes, so a later failure must not leave a
# complete-looking graph for the next command in a pipeline to consume.
rm -f "$OUT/stage".*
"$BIN" bubble -i "$OUT/bad.gfa" -r full -o "$OUT/stage" --min-variant-bp 0 -q >/dev/null 2>&1
[ "$(ls "$OUT/stage".* 2>/dev/null | wc -l | tr -d ' ')" = "0" ] \
  && ok "a failed run leaves no partial output, not even the sorted GFA" \
  || bad "failed run left: $(ls "$OUT/stage".* 2>/dev/null | tr '\n' ' ')"

# ---- an alternate connector belongs to a fused bubble ---------------------------------------------
# Node 5 is a connector branch no path walks, so it is in neither part's interior and not on the
# reference. A fusion assembled from its parts plus the reference connector alone would miss it.
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

# ---- a graph-only allele when every path takes the deletion -------------------------------------
# Node 2 is an allele the graph carries and no stored path walks, and every path crosses the site
# directly source->sink. The endpoint finder must accept adjacent boundaries: that interval is what tells
# the graph-derived search which side of each boundary faces inward, and without it the site would give
# no bubble at all rather than a bubble with one interior node.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGG\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tpA\t1+,3+\t*\nP\tpB\t1+,3+\t*\n'; } > "$OUT/gonly.gfa"
"$BIN" bubble -i "$OUT/gonly.gfa" -r pA -o "$OUT/gonly" --min-variant-bp 0 -q >/dev/null 2>&1
[ "$(rows "$OUT/gonly.bubbles.csv")" = "1" ] \
  && ok "a graph-only allele is found when every path takes the direct route" \
  || bad "expected 1 bubble for the graph-only allele, got $(rows "$OUT/gonly.bubbles.csv")"
[ "$(cel "$OUT/gonly.bubbles.csv" inside_node_count)" = "1" ] \
  && ok "the graph-only allele is the interior of that bubble" \
  || bad "graph-only bubble has inside_node_count $(cel "$OUT/gonly.bubbles.csv" inside_node_count), expected 1"
[ "$(cel "$OUT/gonly.bubbles.csv" source)" = "1" ] && [ "$(cel "$OUT/gonly.bubbles.csv" sink)" = "3" ] \
  && ok "the adjacent interval still establishes reference order" \
  || bad "graph-only bubble is $(cel "$OUT/gonly.bubbles.csv" source)..$(cel "$OUT/gonly.bubbles.csv" sink), expected 1..3"

# ---- --snarls resolves the reference the same way the internal door does ----------------------
# Kept as the user spelled it, `-r FULL` for a path named `full` would match nothing: an imported
# reversed pair would stay reversed and reference allele support read 0, with no error.
printf '{"start": {"node_id": "3"}, "end": {"node_id": "1"}}\n' > "$OUT/revsnarl.jsonl"
"$BIN" bubble -i "$OUT/case.gfa" --snarls "$OUT/revsnarl.jsonl" -r FULL -o "$OUT/snx" \
       --min-variant-bp 0 -q >/dev/null 2>&1
[ "$(cel "$OUT/snx.bubbles.csv" source)" = "1" ] && [ "$(cel "$OUT/snx.bubbles.csv" sink)" = "3" ] \
  && ok "--snarls orients an imported reversed pair from a case-insensitive alias" \
  || bad "--snarls -r FULL gave $(cel "$OUT/snx.bubbles.csv" source)..$(cel "$OUT/snx.bubbles.csv" sink), expected 1..3"
[ "$(cel "$OUT/snx.bubbles.csv" ref_allele_support)" = "1" ] \
  && ok "--snarls finds the reference allele through an alias" \
  || bad "--snarls -r FULL gave ref_allele_support $(cel "$OUT/snx.bubbles.csv" ref_allele_support), expected 1"

# An ambiguous reference must be an error rather than resolved by file order.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGG\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tsampleA\t1+,2+,3+\t*\nP\tsampleB\t1+,3+\t*\n'; } > "$OUT/amb.gfa"
"$BIN" bubble -i "$OUT/amb.gfa" --snarls "$OUT/revsnarl.jsonl" -r sample -o "$OUT/amb" \
       --min-variant-bp 0 -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "--snarls refuses an ambiguous reference" \
               || bad "--snarls resolved an ambiguous reference by file order"
"$BIN" bubble -i "$OUT/amb.gfa" --snarls "$OUT/revsnarl.jsonl" -r nope -o "$OUT/amb2" \
       --min-variant-bp 0 -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "--snarls refuses a reference that is not there" \
               || bad "--snarls accepted a reference that does not exist"

# ---------------------------------------------------------------------------------------------
# Emitted sites must be pairwise disjoint. Two retained sites claiming the same interior node
# describe the same sequence twice, and panphorte refuses to fold both -- so an overlap does not
# degrade a locus, it stops it running at all.
#
# A nested pair: an outer site 1..6 and an inner site 2..5 whose interiors both contain node 3.
# Only one may survive, and it must be the SMALLER one, which keeps the variants inside apart.
# ---------------------------------------------------------------------------------------------
nb() { local n=$1 c=$2 out=""; while [ ${#out} -lt "$n" ]; do out="$out$c"; done; printf '%s' "${out:0:$n}"; }
{ printf 'H\tVN:Z:1.0\n'
  printf "S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t3b\t%s\nS\t4\t%s\nS\t5\t%s\nS\t6\t%s\n" \
         "$(nb 200 A)" "$(nb 100 C)" "$(nb 120 G)" "$(nb 120 T)" "$(nb 100 A)" "$(nb 100 G)" "$(nb 200 T)"
  for e in "1 2" "2 3" "2 3b" "3 4" "3b 4" "4 5" "5 6" "2 4" "1 5"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"
  done
  printf 'P\tnref#0#chrN:1-940\t1+,2+,3+,4+,5+,6+\t*\n'
  printf 'P\tnhapA#1#chrN\t1+,2+,3b+,4+,5+,6+\t*\n'
  printf 'P\tnhapB#1#chrN\t1+,2+,4+,5+,6+\t*\n'
  printf 'P\tnhapC#1#chrN\t1+,5+,6+\t*\n'; } > "$OUT/nest.gfa"
# The internal finder emits TOP-LEVEL snarls only, so it cannot produce a nested overlapping pair --
# a fixture fed through that door yields one site and would assert nothing. --snarls can: it takes
# the pairs it is given. Supplying an outer 1..6 and an inner 2..5, whose interiors both contain 3,
# is the reproducer, and it is also the door a user can actually hit.
printf '{"start": {"node_id": "1"}, "end": {"node_id": "6"}}\n{"start": {"node_id": "2"}, "end": {"node_id": "5"}}\n' \
  > "$OUT/overlap.jsonl"
"$BIN" bubble -i "$OUT/nest.gfa" --snarls "$OUT/overlap.jsonl" -r 'nref#0#chrN:1-940' \
       -o "$OUT/nest" --min-variant-bp 0 -q >/dev/null 2>&1
if [ ! -s "$OUT/nest.bubbles.csv" ]; then
  bad "the overlapping-snarl fixture produced no bubbles CSV"
else
  nsites=$(tail -n +2 "$OUT/nest.bubbles.csv" | wc -l | tr -d ' ')
  clashes=$(tail -n +2 "$OUT/nest.bubbles.csv" | awk -F'"' '{print $2}' |
            tr ';' '\n' | grep -v '^$' | sort | uniq -d | wc -l | tr -d ' ')
  [ "$clashes" = "0" ] \
    && ok "two overlapping imported snarls are resolved to a disjoint set ($nsites site(s))" \
    || bad "$clashes interior node(s) claimed by more than one emitted site"
  # the SMALLER site is the one kept
  kept_inside=$(tail -n +2 "$OUT/nest.bubbles.csv" | head -1 | awk -F'"' '{print $2}')
  case "$kept_inside" in
    *3*) [ "${#kept_inside}" -lt 7 ] \
           && ok "the smaller of two overlapping sites is the one retained ($kept_inside)" \
           || bad "the enclosing site was retained ($kept_inside), not the smaller one" ;;
    *) bad "unexpected retained interior '$kept_inside'" ;;
  esac
fi

# ------------------------------------------- a repeated boundary is not clipped at 64 occurrences
# A path revisiting a boundary many times must have its interval run to the last visit, not close after
# a fixed number of endpoints. The graph-derived interior would hide a clipped interval in
# `inside_node_count` -- every unit variant is still reachable between the boundaries -- so the
# observable is the walked span, `max_inside_bp`.
#
# 70 copies of a 20 bp unit variant between two boundaries that each recur 70 times. The boundaries are
# inside the loop, so the site is widened to the nodes around it (1..4), and a copy is the unit variant
# with both boundaries: 70 x 60 = 4200 bp; a search capped at 64 occurrences would report 64 x 60 = 3840.
{ printf 'H\tVN:Z:1.0\nS\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\n' \
    "$(printf 'A%.0s' $(seq 60))" "$(printf 'C%.0s' $(seq 20))" \
    "$(printf 'G%.0s' $(seq 20))" "$(printf 'T%.0s' $(seq 60))"
  for i in $(seq 1 70); do printf 'S\t%d\tACGTACGTACGTACGTACGT\n' $((100 + i)); done
  printf 'L\t1\t+\t2\t+\t0M\nL\t3\t+\t4\t+\t0M\nL\t3\t+\t2\t+\t0M\n'
  for i in $(seq 1 70); do printf 'L\t2\t+\t%d\t+\t0M\nL\t%d\t+\t3\t+\t0M\n' $((100 + i)) $((100 + i)); done
  printf 'P\tarray\t1+'; for i in $(seq 1 70); do printf ',2+,%d+,3+' $((100 + i)); done; printf ',4+\t*\n'
  printf 'P\tsingle\t1+,2+,101+,3+,4+\t*\n'; } > "$OUT/repeat70.gfa"
printf '{"start": {"node_id": "2"}, "end": {"node_id": "3"}}\n' > "$OUT/repeat70.jsonl"
"$BIN" bubble -i "$OUT/repeat70.gfa" -r array --snarls "$OUT/repeat70.jsonl" \
       -o "$OUT/repeat70" --min-variant-bp 0 -q >/dev/null 2>&1
r70_max=$(cel "$OUT/repeat70.bubbles.csv" max_inside_bp)
[ "$r70_max" = "4200" ] \
  && ok "a 70-fold repeated boundary keeps its full interior (max_inside_bp 4200)" \
  || bad "max_inside_bp=$r70_max, expected 4200 (64-endpoint clipping reports 3840)"
r70_in=$(cel "$OUT/repeat70.bubbles.csv" inside_node_count)
[ "$r70_in" = "72" ] && ok "all 70 unit variants and the two old boundaries are interior nodes" \
                     || bad "inside_node_count=$r70_in, expected 72"
# the single-copy haplotype still supplies the shortest interior: one copy, 60 bp
r70_min=$(cel "$OUT/repeat70.bubbles.csv" min_inside_bp)
[ "$r70_min" = "60" ] && ok "the one-copy haplotype still sets min_inside_bp (60)" \
                      || bad "min_inside_bp=$r70_min, expected 60"

# ------------------------------------------- two outputs may not name one file
# Two outputs naming one file are refused, as an output naming an input is: otherwise the run succeeds
# and leaves whichever file happened to be installed last. Both rules live in cli_utils, shared with
# `inspect`.
"$BIN" bubble -i "$OUT/del.gfa" -r full -o "$OUT/coll" --min-variant-bp 0 \
       --bubbles-csv "$OUT/both.csv" --bandage-csv "$OUT/both.csv" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && [ ! -f "$OUT/both.csv" ] \
  && ok "two outputs naming one file are refused, and nothing is written" \
  || bad "colliding --bubbles-csv/--bandage-csv were accepted"
"$BIN" bubble -i "$OUT/del.gfa" -r full -o "$OUT/coll2" --min-variant-bp 0 \
       --bubbles-csv "$OUT/del.gfa" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "an output naming the input graph is refused" \
              || bad "an output aliasing the input was accepted"

# ------------------------------------------- a failed commit leaves the previous family untouched
# StagedOutputs is shared by every module, so this asserts a cross-module contract through the module
# that happens to write the smallest multi-file family. Outputs are staged and installed one after
# another, so a failure on the fourth of six must not leave three new files beside three stale ones,
# with a non-zero exit and nothing saying which is which.
#
# The failure is injected, because it cannot be reached otherwise: the natural ways to make a rename
# fail here either fail earlier, at staging, or are absorbed by the set-aside that protects the
# previous file. Run one writes a family from the deletion graph; run two writes a DIFFERENT graph to
# the same prefix and fails installing the second file, so a partial install would be plainly visible.
"$BIN" bubble -i "$OUT/del.gfa" -r full -o "$OUT/txn" --min-variant-bp 0 -q >/dev/null 2>&1
txn_csv_before=$(cksum < "$OUT/txn.bubbles.csv"); txn_gfa_before=$(cksum < "$OUT/txn.sorted.gfa")
"$BIN" bubble -i "$OUT/cyc.gfa" -r pA -o "$OUT/other" --min-variant-bp 0 -q >/dev/null 2>&1
if [ "$(cksum < "$OUT/other.bubbles.csv")" = "$txn_csv_before" ]; then
  bad "the rollback fixture is not discriminating: both graphs write the same bubbles.csv"
else
  PANVAR_TEST_FAIL_COMMIT_AT=2 "$BIN" bubble -i "$OUT/cyc.gfa" -r pA -o "$OUT/txn" \
      --min-variant-bp 0 -q >/dev/null 2>&1
  [ "$?" -ne 0 ] && ok "an injected commit failure exits non-zero" \
                 || bad "an injected commit failure exited 0"
  [ "$(cksum < "$OUT/txn.bubbles.csv")" = "$txn_csv_before" ] \
    && ok "the previous bubbles.csv survives a failed commit unchanged" \
    || bad "bubbles.csv was replaced by a run that did not complete"
  [ "$(cksum < "$OUT/txn.sorted.gfa")" = "$txn_gfa_before" ] \
    && ok "the previous sorted.gfa survives a failed commit unchanged" \
    || bad "sorted.gfa was replaced by a run that did not complete"
  residue=$(ls "$OUT" | grep -c -- '-tmp\.\|-prev\.')
  [ "$residue" = "0" ] && ok "a failed commit leaves no staging or set-aside residue" \
                       || bad "$residue staging/set-aside file(s) left behind"
fi

# The window between moving the old file aside and installing the new one. Rollback must own it: the
# previous output must come back from its `*-prev` file rather than stay under a name nobody would look
# for, with the destination gone.
TX2="$OUT/tx2"; mkdir -p "$TX2"
"$BIN" bubble -i "$OUT/del.gfa" -r full -o "$TX2/b" --min-variant-bp 0 -q >/dev/null 2>&1
tx2_csv=$(cksum < "$TX2/b.bubbles.csv"); tx2_gfa=$(cksum < "$TX2/b.sorted.gfa")
tx2_n=$(ls "$TX2" | wc -l | tr -d ' ')
PANVAR_TEST_FAIL_COMMIT_AFTER_SETASIDE=2 "$BIN" bubble -i "$OUT/cyc.gfa" -r pA -o "$TX2/b" \
    --min-variant-bp 0 -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a failure between set-aside and install exits non-zero" \
               || bad "a failure between set-aside and install exited 0"
[ "$(ls "$TX2" | wc -l | tr -d ' ')" = "$tx2_n" ] \
  && ok "every previous output is still present after that failure" \
  || bad "the family lost a file: $(ls "$TX2" | tr '\n' ' ')"
[ "$(cksum < "$TX2/b.bubbles.csv" 2>/dev/null)" = "$tx2_csv" ] \
  && [ "$(cksum < "$TX2/b.sorted.gfa" 2>/dev/null)" = "$tx2_gfa" ] \
  && ok "and each is byte-unchanged, restored from its reserve" \
  || bad "a restored file does not match what was there before"
tx2_res=$(ls "$TX2" | grep -c -- '-tmp\.\|-prev\.')
[ "$tx2_res" = "0" ] && ok "no set-aside reserve is orphaned" \
                     || bad "$tx2_res orphaned reserve(s); the previous ordering left exactly this"

# An output this run CREATED, where none existed, must also be removed on rollback -- otherwise a
# failed run leaves half a new family behind.
TX3="$OUT/tx3"; mkdir -p "$TX3"
PANVAR_TEST_FAIL_COMMIT_AT=3 "$BIN" bubble -i "$OUT/del.gfa" -r full -o "$TX3/fresh" \
    --min-variant-bp 0 -q >/dev/null 2>&1
tx3_left=$(ls "$TX3" 2>/dev/null | wc -l | tr -d ' ')
[ "$tx3_left" = "0" ] && ok "a failed FIRST run leaves no partial new family" \
                      || bad "$tx3_left file(s) left by a run that never completed: $(ls "$TX3" | tr '\n' ' ')"

# A run killed with its outputs staged (Ctrl-C, kill) removes them on the way out, then dies of the
# signal as usual: exit status 128 + 15.
TX4="$OUT/tx4"; mkdir -p "$TX4"
PANVAR_TEST_RAISE_AT_COMMIT=1 "$BIN" bubble -i "$OUT/del.gfa" -r full -o "$TX4/k" --min-variant-bp 0 -q \
    >/dev/null 2>&1
tx4_rc=$?
tx4_left=$(ls "$TX4" 2>/dev/null | wc -l | tr -d ' ')
[ "$tx4_rc:$tx4_left" = "143:0" ] && ok "a killed run leaves no staged file and exits by the signal" \
                                  || bad "killed run: exit $tx4_rc, $tx4_left file(s) left: $(ls "$TX4" | tr '\n' ' ')"
# A staged file is named after its destination plus a suffix; a destination near the 255-byte limit
# gets a shorter staged name, so an output name that fits is written.
TX5="$OUT/tx5"; mkdir -p "$TX5"
long=$(printf 'x%.0s' $(seq 1 230))
"$BIN" bubble -i "$OUT/del.gfa" -r full -o "$TX5/$long" --min-variant-bp 0 -q >/dev/null 2>&1
[ -s "$TX5/$long.bubbles.csv" ] && ok "a 230-character output prefix is written" \
                                || bad "a 230-character output prefix failed"
# Killed while outputs are being installed: the signal waits for the install, so the new family is whole and
# no previous output is left parked under a temporary name.
TX6="$OUT/tx6"; mkdir -p "$TX6"
"$BIN" bubble -i "$OUT/del.gfa" -r full -o "$TX6/k" --min-variant-bp 0 -q >/dev/null 2>&1
cp "$TX6/k.bubbles.csv" "$TX6/prev.csv"; printf 'old\n' > "$TX6/k.bubbles.csv"
PANVAR_TEST_RAISE_IN_COMMIT=1 "$BIN" bubble -i "$OUT/del.gfa" -r full -o "$TX6/k" --min-variant-bp 0 -q \
    >/dev/null 2>&1
tx6="$?:$(cmp -s "$TX6/k.bubbles.csv" "$TX6/prev.csv" && echo new || echo old):$(ls "$TX6" | grep -c -- '-tmp\.')"
[ "$tx6" = "143:new:0" ] && ok "a run killed while installing finishes the install first" \
                         || bad "killed while installing: exit:bubbles.csv:leftovers = $tx6"
# A P line ending in a comma has an empty step, refused here as by every other module.
{ printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTAC\nS\t2\tTTTTTTTTTT\nL\t1\t+\t2\t+\t0M\n'
  printf 'P\tref\t1+,2+,\t*\n'; } > "$OUT/trail.gfa"
e=$("$BIN" bubble -i "$OUT/trail.gfa" -r ref -o "$OUT/trail" -q 2>&1 >/dev/null)
case "$e" in *"ends with a comma"*) ok "a P line ending in a comma is refused" ;; *) bad "trailing comma: ${e:-accepted}" ;; esac
# A link to a segment the file does not define is refused, not dropped.
{ printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTAC\nS\t2\tTTTTTTTTTT\nL\t1\t+\t2\t+\t0M\nL\t2\t+\t9\t+\t0M\n'
  printf 'P\tref\t1+,2+\t*\n'; } > "$OUT/dang.gfa"
printf 'bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,1,+,2,+,1,0,0,\n' \
  > "$OUT/dang.csv"
e=$("$BIN" inspect -i "$OUT/dang.gfa" -c "$OUT/dang.csv" -o "$OUT/dang" -q 2>&1 >/dev/null)
case "$e" in *"segment 9 is not in the graph"*) ok "a link to an undefined segment is refused" ;; *) bad "dangling link: ${e:-accepted}" ;; esac
# An empty file is named as such, before any reference lookup.
: > "$OUT/empty.gfa"
e=$("$BIN" bubble -i "$OUT/empty.gfa" -r full -o "$OUT/emp" -q 2>&1 >/dev/null)
case "$e" in *"no segments and no paths"*) ok "an empty GFA is reported as empty" ;;
             *) bad "empty GFA: $e" ;; esac

# ---------------------------------------------------------------- damaged compressed input
# Every GFA and GTF goes through one reader. zlib reports a stream that stops early only when asked,
# and a BGZF write interrupted between blocks is a valid stream, so a damaged graph could parse as a
# shorter valid one -- GFA writers put the path lines last, so what would go missing is haplotypes.
gz_err() { "$BIN" bubble -i "$1" -r ref --min-variant-bp 0 -o "$OUT/gzchk" -q 2>&1 >/dev/null; }
gzip -c "$OUT/alt.gfa" > "$OUT/gzok.gfa.gz"
e=$(gz_err "$OUT/gzok.gfa.gz")
[ -z "$e" ] && ok "an intact gzipped GFA is read" || bad "an intact gzipped GFA was rejected: $e"
gzsz=$(wc -c < "$OUT/gzok.gfa.gz" | tr -d ' ')
head -c $((gzsz - 6)) "$OUT/gzok.gfa.gz" > "$OUT/gzcut.gfa.gz"
e=$(gz_err "$OUT/gzcut.gfa.gz")
case "$e" in *"the file is truncated"*) ok "a truncated gzipped GFA is rejected as truncated" ;;
             *) bad "a truncated gzipped GFA was accepted: $e" ;; esac
cp "$OUT/gzok.gfa.gz" "$OUT/gzflip.gfa.gz"
python3 -I -c "import sys;p=sys.argv[1];b=bytearray(open(p,'rb').read());b[len(b)//2]^=0xFF;open(p,'wb').write(b)" \
  "$OUT/gzflip.gfa.gz"
e=$(gz_err "$OUT/gzflip.gfa.gz")
case "$e" in *"cannot decompress"*) ok "a corrupted gzipped GFA is rejected" ;;
             *) bad "a corrupted gzipped GFA was accepted: $e" ;; esac
if command -v bgzip >/dev/null; then
  bgzip -c "$OUT/alt.gfa" > "$OUT/bgok.gfa.gz"
  e=$(gz_err "$OUT/bgok.gfa.gz")
  [ -z "$e" ] && ok "an intact BGZF GFA is read" || bad "an intact BGZF GFA was rejected: $e"
  bgsz=$(wc -c < "$OUT/bgok.gfa.gz" | tr -d ' ')
  head -c $((bgsz - 28)) "$OUT/bgok.gfa.gz" > "$OUT/bgnoeof.gfa.gz"
  e=$(gz_err "$OUT/bgnoeof.gfa.gz")
  case "$e" in *"no end-of-file marker"*) ok "a BGZF GFA cut on a block boundary is rejected" ;;
               *) bad "a BGZF GFA without its EOF marker was accepted: $e" ;; esac
else
  printf "  --   BGZF EOF-marker checks (skipped: bgzip not installed)\n"
fi


# ---------------------------------------------------------------- a failed fusion keeps its parts
# Two 600 bp sites pass --max-variant-bp 1000 on their own; fused they span 1210 bp and fail. The
# fusion is undone and both sites kept, not erased with it.
"${PYTHON:-python3}" - "$OUT/fuse.gfa" <<'PYEOF'
import random, sys
r = random.Random(13)
rnd = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
S = {'1': rnd(3500), '2': rnd(600), '3': rnd(600), '4': rnd(10), '5': rnd(600), '6': rnd(600), '7': rnd(3000)}
L = [('1','2'),('1','3'),('2','4'),('3','4'),('4','5'),('4','6'),('5','7'),('6','7')]
P = {'pA': '1+,2+,4+,5+,7+', 'pB': '1+,3+,4+,6+,7+'}
with open(sys.argv[1], 'w') as o:
    o.write('H\tVN:Z:1.0\n')
    for k, v in S.items(): o.write(f'S\t{k}\t{v}\n')
    for a, b in L: o.write(f'L\t{a}\t+\t{b}\t+\t0M\n')
    for k, v in P.items(): o.write(f'P\t{k}\t{v}\t*\n')
PYEOF
"$BIN" bubble -i "$OUT/fuse.gfa" -r pA -o "$OUT/fuse" -q --max-variant-bp 1000 --merge-nearby-bp 50 \
  > "$OUT/fuse.log" 2>&1
n=$(( $(wc -l < "$OUT/fuse.bubbles.csv") - 1 ))
[ "$n" = "2" ] && ok "a fusion that fails the filters is undone, keeping both parts" \
               || bad "a failed fusion left $n site(s), expected 2"
grep -q "merged only the runs of their parts" "$OUT/fuse.log" && ok "  ... and says so even with -q" \
                                                 || bad "  ... without saying so"

# ---------------------------------------------------------------- sites whose variants can move into each other
# A 72 bp (AC)n microsatellite holds two sites: a 2 bp deletion (nodes 2..5, under the size floor) and
# an expansion (5..8, kept). Either change can slide along the whole repeat and stay the same
# haplotype, so the two sites are one locus: the expansion's site is fused with the small one and runs
# from node 2; --no-fuse-ambiguous keeps it at 5..8 (4..7 once sorting has renumbered the nodes along
# the reference). The control has the same shape built from homopolymer runs (10 A, a C, 9 T, ...),
# where neither change can move, so the sites stay apart.
trgfa() {   # trgfa <out> <n2> <n3> <n5> <n6> <n7> <n8>
  { printf 'H\tVN:Z:1.0\n'
    "${PYTHON:-python3}" -c 'import random; r=random.Random(29); print("S\t1\t"+"".join(r.choice("ACGT") for _ in range(300))); print("S\t9\t"+"".join(r.choice("ACGT") for _ in range(300)))'
    printf 'S\t2\t%s\nS\t3\t%s\nS\t5\t%s\nS\t6\t%s\nS\t7\t%s\nS\t8\t%s\n' "$2" "$3" "$4" "$5" "$6" "$7"
    for e in "1 2" "2 3" "3 5" "2 5" "5 6" "5 7" "6 8" "7 8" "8 9"; do
      set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
    printf 'P\tref\t1+,2+,3+,5+,6+,8+,9+\t*\nP\thdel\t1+,2+,5+,6+,8+,9+\t*\nP\thexp\t1+,2+,3+,5+,7+,8+,9+\t*\n'; } > "$1"
}
ac() { local o=""; for _ in $(seq "$1"); do o="${o}AC"; done; printf '%s' "$o"; }
trgfa "$OUT/tr.gfa" "$(ac 10)" "AC" "$(ac 10)" "$(ac 5)" "$(ac 40)" "$(ac 10)"
"$BIN" bubble -i "$OUT/tr.gfa" -r ref -o "$OUT/tr" --merge-nearby-bp 0 -q >/dev/null 2>&1
"$BIN" bubble -i "$OUT/tr.gfa" -r ref -o "$OUT/trn" --merge-nearby-bp 0 --no-fuse-ambiguous -q >/dev/null 2>&1
[ "$(rows "$OUT/tr.bubbles.csv"):$(cel "$OUT/tr.bubbles.csv" source):$(cel "$OUT/tr.bubbles.csv" sink)" = "1:2:7" ] \
  && ok "two sites whose changes slide into each other in a repeat are fused (rows:source:sink 1:2:7)" \
  || bad "placement fusion gave $(rows "$OUT/tr.bubbles.csv"):$(cel "$OUT/tr.bubbles.csv" source):$(cel "$OUT/tr.bubbles.csv" sink), expected 1:2:7"
[ "$(cel "$OUT/trn.bubbles.csv" source):$(cel "$OUT/trn.bubbles.csv" sink)" = "4:7" ] \
  && ok "--no-fuse-ambiguous leaves it at its own boundaries (4:7)" \
  || bad "--no-fuse-ambiguous gave $(cel "$OUT/trn.bubbles.csv" source):$(cel "$OUT/trn.bubbles.csv" sink), expected 4:7"
# The same microsatellite with one base changed between the two sites (AC -> GC): exact sliding stops
# there, but a change placed across it costs one mismatch, less than another gap, so the sites are fused.
trgfa "$OUT/tri.gfa" "$(ac 10)" "AC" "$(ac 4)GC$(ac 5)" "$(ac 5)" "$(ac 40)" "$(ac 10)"
"$BIN" bubble -i "$OUT/tri.gfa" -r ref -o "$OUT/tri" --merge-nearby-bp 0 -q >/dev/null 2>&1
[ "$(rows "$OUT/tri.bubbles.csv"):$(cel "$OUT/tri.bubbles.csv" source):$(cel "$OUT/tri.bubbles.csv" sink)" = "1:2:7" ] \
  && ok "an interrupted repeat: sites on either side of a changed base are fused (rows:source:sink 1:2:7)" \
  || bad "interrupted repeat gave $(rows "$OUT/tri.bubbles.csv"):$(cel "$OUT/tri.bubbles.csv" source):$(cel "$OUT/tri.bubbles.csv" sink), expected 1:2:7"
trgfa "$OUT/hp.gfa" "AAAAAAAAAACTTTTTTTTT" "GG" "CCCCCCCCCCATTTTTTTTT" "GGGGGGGGGG" \
      "GGGGGGGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCCGGGGGGGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCC" "AAAAAAAAAACGGGGGGGGG"
"$BIN" bubble -i "$OUT/hp.gfa" -r ref -o "$OUT/hp" --merge-nearby-bp 0 -q >/dev/null 2>&1
[ "$(cel "$OUT/hp.bubbles.csv" source):$(cel "$OUT/hp.bubbles.csv" sink)" = "4:7" ] \
  && ok "sites whose changes cannot move into each other stay apart (4:7)" \
  || bad "unmovable changes were fused: $(cel "$OUT/hp.bubbles.csv" source):$(cel "$OUT/hp.bubbles.csv" sink), expected 4:7"
# A repeat the reference does not have: it carries ONE copy (node 3) of a 30 bp unit, and h1/h2 three
# more, which the graph put around it -- one before it (node 6, a 30 bp site under the size floor) and
# two after (node 7, a 60 bp site, kept). Each insertion can slide across the reference copy, so the
# sites are one: S..K, nodes 2..6 once sorted (4..6 without the fusion).
"${PYTHON:-python3}" - "$OUT/nr.gfa" <<'PYEOF'
import random, sys
r = random.Random(31)
rnd = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
U = rnd(30)
S = {'1': rnd(300), '2': rnd(40), '3': U, '4': rnd(40), '5': rnd(300), '6': U, '7': U + U}
L = [('1','2'), ('2','3'), ('3','4'), ('4','5'), ('2','6'), ('6','3'), ('3','7'), ('7','4')]
P = {'ref': '1+,2+,3+,4+,5+', 'h1': '1+,2+,6+,3+,7+,4+,5+', 'h2': '1+,2+,6+,3+,7+,4+,5+'}
with open(sys.argv[1], 'w') as o:
    o.write('H\tVN:Z:1.0\n')
    for k, v in S.items(): o.write(f'S\t{k}\t{v}\n')
    for a, b in L: o.write(f'L\t{a}\t+\t{b}\t+\t0M\n')
    for k, v in P.items(): o.write(f'P\t{k}\t{v}\t*\n')
PYEOF
"$BIN" bubble -i "$OUT/nr.gfa" -r ref -o "$OUT/nr" --merge-nearby-bp 0 -q >/dev/null 2>&1
"$BIN" bubble -i "$OUT/nr.gfa" -r ref -o "$OUT/nrn" --merge-nearby-bp 0 --no-fuse-ambiguous -q >/dev/null 2>&1
[ "$(rows "$OUT/nr.bubbles.csv"):$(cel "$OUT/nr.bubbles.csv" source):$(cel "$OUT/nr.bubbles.csv" sink):$(cel "$OUT/nrn.bubbles.csv" source):$(cel "$OUT/nrn.bubbles.csv" sink)" = "1:2:6:4:6" ] \
  && ok "an expansion the graph split around the reference's single copy is one site (fused 2:6, unfused 4:6)" \
  || bad "single-copy expansion: $(rows "$OUT/nr.bubbles.csv"):$(cel "$OUT/nr.bubbles.csv" source):$(cel "$OUT/nr.bubbles.csv" sink):$(cel "$OUT/nrn.bubbles.csv" source):$(cel "$OUT/nrn.bubbles.csv" sink), expected 1:2:6:4:6"
"$BIN" call -i "$OUT/nr.sorted.gfa" -c "$OUT/nr.bubbles.csv" -r ref -o "$OUT/nrc" -q >/dev/null 2>&1
nr_rec=$(grep -v '^#' "$OUT/nrc.region.vcf" | awk -F'\t' '{n=split($8,a,";"); for(i=1;i<=n;i++){split(a[i],kv,"="); if(kv[1]=="SVLEN") l=kv[2]; if(kv[1]=="AC") c=kv[2]} print $2":"l":"c}' | paste -sd, -)
[ "$nr_rec" = "340:90:2" ] \
  && ok "  ... and call reports it as one 90 bp INS on both haplotypes (POS:SVLEN:AC 340:90:2)" \
  || bad "single-copy expansion call: $nr_rec, expected 340:90:2"

# ---------------------------------------------------------------- placement fusion: what it must not do
# Unique sequence: a retained 60 bp deletion and a 1 bp deletion 6 bp from it, neither able to move even
# one base, so they stay apart (the site stays 1..3).
cat > "$OUT/uniq.gfa" <<'GFAEOF'
H	VN:Z:1.0
S	1	TTTCCTCATGCAATTCAAAACCATGTCCGT
S	2	AATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGAC
S	3	CTAACC
S	4	G
S	5	TGAGGTAAACCAGGTCTCTCCGCCCCCTTA
L	1	+	2	+	0M
L	1	+	3	+	0M
L	2	+	3	+	0M
L	3	+	4	+	0M
L	3	+	5	+	0M
L	4	+	5	+	0M
P	REF	1+,2+,3+,4+,5+	*
P	h1	1+,3+,4+,5+	*
P	h2	1+,2+,3+,5+	*
P	h3	1+,3+,5+	*
P	h4	1+,2+,3+,4+,5+	*
GFAEOF
"$BIN" bubble -i "$OUT/uniq.gfa" -r REF -o "$OUT/uniq" -q >/dev/null 2>&1
[ "$(tail -n +2 "$OUT/uniq.bubbles.csv" | cut -d, -f2,4 | paste -sd' ' -)" = "1,3" ] \
  && ok "changes in unique sequence a few bases apart are not fused (1,3)" \
  || bad "unique-sequence changes were fused: $(tail -n +2 "$OUT/uniq.bubbles.csv" | cut -d, -f2,4 | paste -sd' ' -), expected 1,3"
# Two short deletions in unique sequence, each pinning one of the two placement rules. In ug1 each
# could reach the other across three differing bases, but only by both moving: the tolerance is spent
# once, by the change that moves. In ug2 a placement past a mismatch that ends on another mismatch is
# not a placement: the stretch ends only on one that matches.
ugfa() {   # ugfa <out> <n1> <n2> <n3> <n4> <n5>
  { printf 'H\tVN:Z:1.0\n'
    printf 'S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\nS\t5\t%s\n' "$2" "$3" "$4" "$5" "$6"
    for e in "1 2" "1 3" "2 3" "3 4" "3 5" "4 5"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
    printf 'P\tREF\t1+,2+,3+,4+,5+\t*\nP\th1\t1+,3+,4+,5+\t*\nP\th2\t1+,2+,3+,5+\t*\n'
    printf 'P\th3\t1+,2+,3+,4+,5+\t*\nP\th4\t1+,3+,4+,5+\t*\nP\th5\t1+,2+,3+,5+\t*\n'; } > "$1"
}
ugfa "$OUT/ug1.gfa" CTGATAACACAGCTTTAAGC GA GGATACAT CGA TCATCTGTAATCGGCGGGCTGCGGTAC
ugfa "$OUT/ug2.gfa" GGCTTAATGTGCTTAAATAC AA TGG AC ACTTTAGCGTTGGTGATTTAGAAACTTTGGCAT
for g in ug1 ug2; do
  "$BIN" bubble -i "$OUT/$g.gfa" -r REF -o "$OUT/$g" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
done
[ "$(tail -n +2 "$OUT/ug1.bubbles.csv" | cut -d, -f2,4 | paste -sd' ' -)/$(tail -n +2 "$OUT/ug2.bubbles.csv" | cut -d, -f2,4 | paste -sd' ' -)" = "1,3 3,5/1,3 3,5" ] \
  && ok "deletions that meet only if both move, or past a mismatch, are not fused (1,3 3,5 twice)" \
  || bad "placement rules: ug1 $(tail -n +2 "$OUT/ug1.bubbles.csv" | cut -d, -f2,4 | paste -sd' ' -), ug2 $(tail -n +2 "$OUT/ug2.bubbles.csv" | cut -d, -f2,4 | paste -sd' ' -), expected 1,3 3,5 each"
# In a CA microsatellite: a 2 bp deletion site under the floor (1..3), a retained site whose only
# alternative spells the reference's sequence (3..6: no change, so no stretch), 120 bp of CA, then a
# retained site deleting one CA (8..11). The first and last can be placed in each other across the
# middle site; fused, they would enclose it and the fusion would be dropped in its favour, losing 8..11.
# The fusion is not made: both retained sites stay.
cat > "$OUT/loss.gfa" <<'GFAEOF'
H	VN:Z:1.0
S	1	GATTCGGTACCTGAAGTCTG
S	2	CA
S	3	CACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA
S	4	CACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA
S	5	CACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA
S	6	CA
S	7	CACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA
S	8	CA
S	9	CACACACACACACACACACACACACACACACACACACACACACACACACACACACACACA
S	10	CACACACACACACACACACACACACACACACACACACACACACACACACACACACACA
S	11	TTGGACCTAGCATGCAGTAG
L	1	+	2	+	0M
L	1	+	3	+	0M
L	10	+	11	+	0M
L	2	+	3	+	0M
L	3	+	4	+	0M
L	3	+	5	+	0M
L	4	+	6	+	0M
L	5	+	6	+	0M
L	6	+	7	+	0M
L	7	+	8	+	0M
L	8	+	10	+	0M
L	8	+	9	+	0M
L	9	+	11	+	0M
P	REF	1+,2+,3+,4+,6+,7+,8+,9+,11+	*
P	h1	1+,3+,4+,6+,7+,8+,9+,11+	*
P	h2	1+,2+,3+,5+,6+,7+,8+,10+,11+	*
P	h3	1+,2+,3+,4+,6+,7+,8+,9+,11+	*
P	h4	1+,2+,3+,4+,6+,7+,8+,10+,11+	*
GFAEOF
"$BIN" bubble -i "$OUT/loss.gfa" -r REF -o "$OUT/loss" -q >/dev/null 2>&1
[ "$(tail -n +2 "$OUT/loss.bubbles.csv" | cut -d, -f2,4 | paste -sd' ' -)" = "3,6 8,11" ] \
  && ok "a fusion that would enclose another retained site is not made (3,6 8,11)" \
  || bad "enclosed site: got $(tail -n +2 "$OUT/loss.bubbles.csv" | cut -d, -f2,4 | paste -sd' ' -), expected 3,6 8,11"
# A site enclosing the repeat: hbyp replaces it whole, and a cycle (hdup) with truncated haplotypes
# makes it a locus wrapper, retained with its children. The overlap pass drops it for the smaller sites
# inside, so it neither joins nor blocks their fusion: the 2 bp deletion and the expansion are one site
# (5..10), and hexp2 (-40 -10 +80 against the reference, +30 overall) is not called for the +70 INS.
"${PYTHON:-python3}" - "$OUT/wrap.gfa" <<'PYEOF'
import random, sys
R = random.Random(1); r = lambda n: ''.join(R.choice('ACGT') for _ in range(n))
ac = lambda n: 'AC' * n
S = {'1': r(300), '2': r(200), '4': r(300), '7': r(300), '6': r(10),
     'm2': r(40) + ac(10), 'm3': ac(20), 'm5': ac(10), 'm6': ac(5), 'm7': ac(40), 'm8': ac(10) + r(60)}
L = [('1','2'), ('2','4'), ('4','m2'), ('m2','m3'), ('m3','m5'), ('m2','m5'), ('m5','m6'), ('m5','m7'),
     ('m6','m8'), ('m7','m8'), ('m8','7'), ('4','6'), ('6','7'), ('7','1')]
core = 'm2 m3 m5 m6 m8'
P = {'REF': f'1 2 4 {core} 7', 'hdel': '1 2 4 m2 m5 m6 m8 7', 'hexp': '1 2 4 m2 m3 m5 m7 m8 7',
     'hexp2': '1 2 4 m2 m5 m7 m8 7', 'hbyp': '1 2 4 6 7', 'hdup': f'1 2 4 {core} 7 1 2 4 6 7'}
for i in range(8):
    P[f't{i}'] = ['1 2 4 m2 m3 m5', '1 2 4 m2 m3 m5 m6 m8', '1 2 4 m2 m5 m6', '1 2 4 m2 m3 m5 m7'][i % 4]
with open(sys.argv[1], 'w') as o:
    o.write('H\tVN:Z:1.0\n')
    for k, v in S.items(): o.write(f'S\t{k}\t{v}\n')
    for a, b in L: o.write(f'L\t{a}\t+\t{b}\t+\t0M\n')
    for k, v in P.items(): o.write(f'P\t{k}\t' + ','.join(x + '+' for x in v.split()) + '\t*\n')
PYEOF
"$BIN" bubble -i "$OUT/wrap.gfa" -r REF -o "$OUT/wrap" -q >/dev/null 2>&1
"$BIN" call -i "$OUT/wrap.sorted.gfa" -c "$OUT/wrap.bubbles.csv" -r REF -o "$OUT/wrapc" -q >/dev/null 2>&1
wrap_gt=$(awk -F'\t' '/^#CHROM/{for(i=10;i<=NF;i++) c[$i]=i} !/^#/ && /SVLEN=70;/{split($c["hexp"],a,":"); split($c["hexp2"],b,":"); print a[1]":"b[1]}' "$OUT/wrapc.region.vcf")
[ "$(tail -n +2 "$OUT/wrap.bubbles.csv" | cut -d, -f2,4 | paste -sd' ' -):$wrap_gt" = "5,10:1:0" ] \
  && ok "a site enclosing a fusable repeat neither joins nor blocks the fusion (5,10; +70 INS GT hexp:hexp2 1:0)" \
  || bad "wrapper: got $(tail -n +2 "$OUT/wrap.bubbles.csv" | cut -d, -f2,4 | paste -sd' ' -):$wrap_gt, expected 5,10:1:0"
# The same wrapper, its sink W now inside the AC run, 120 bp past the expansion, and a retained site after
# it (hk3 expands W..z1 by 90 bp). The expansion and that site can be placed in each other across W. The
# overlap pass keeps, of two sites sharing nodes, the one with fewer: with a one-node bypass the wrapper
# has fewer nodes than the fused site would, so it would drop the fused site and both retained sites in
# it, and the fusion is not made (7,10 12,15). With a six-node bypass the wrapper is the larger, dropped
# whatever is fused, and the fusion is made (7,20).
straddle_gfa() {   # straddle_gfa <out> <bypass nodes>
  "${PYTHON:-python3}" - "$1" "$2" <<'PYEOF'
import random, sys
R = random.Random(1); r = lambda n: ''.join(R.choice('ACGT') for _ in range(n))
ac = lambda n: 'AC' * n
nb = int(sys.argv[2])   # bypass chain length
S = {'1': r(300), '2': r(200), '4': r(300), '7': r(300),
     'm2': r(40) + ac(10), 'm3': ac(20), 'm5': ac(10), 'm6': ac(5), 'm7': ac(40), 'm8': ac(10),
     'g2': ac(60), 'W': ac(5), 'k1': ac(5), 'ka': ac(50), 'z1': r(60)}
byp = [f'b{i}' for i in range(nb)]
for b in byp: S[b] = r(10)
L = [('1','2'), ('2','4'), ('4','m2'), ('m2','m3'), ('m3','m5'), ('m2','m5'), ('m5','m6'), ('m5','m7'),
     ('m6','m8'), ('m7','m8'), ('m8','g2'), ('g2','W'), ('W','k1'), ('W','ka'),
     ('k1','z1'), ('ka','z1'), ('z1','7'), ('7','1')]
chain = ['4'] + byp + ['W']
L += list(zip(chain, chain[1:]))
core = 'm2 m3 m5 m6 m8 g2 W'
P = {'REF': f'1 2 4 {core} k1 z1 7', 'hdel': '1 2 4 m2 m5 m6 m8 g2 W k1 z1 7',
     'hexp': '1 2 4 m2 m3 m5 m7 m8 g2 W ka z1 7', 'hexp2': '1 2 4 m2 m3 m5 m7 m8 g2 W k1 z1 7',
     'hbyp': '1 2 4 ' + ' '.join(byp) + ' W k1 z1 7', 'hk3': f'1 2 4 {core} ka z1 7',
     'hdup': f'1 2 4 {core} k1 z1 7 1 2 4 ' + ' '.join(byp) + ' W k1 z1 7'}
for i in range(8):
    P[f't{i}'] = ['1 2 4 m2 m3 m5', '1 2 4 m2 m3 m5 m6 m8', '1 2 4 m2 m5 m6', '1 2 4 m2 m3 m5 m7'][i % 4]
with open(sys.argv[1], 'w') as o:
    o.write('H\tVN:Z:1.0\n')
    for k, v in S.items(): o.write(f'S\t{k}\t{v}\n')
    for a, b in L: o.write(f'L\t{a}\t+\t{b}\t+\t0M\n')
    for k, v in P.items(): o.write(f'P\t{k}\t' + ','.join(x + '+' for x in v.split()) + '\t*\n')
PYEOF
}
for d in "1 7,10 12,15" "6 7,20"; do
  n=${d%% *}; want=${d#* }
  straddle_gfa "$OUT/straddle$n.gfa" "$n"
  "$BIN" bubble -i "$OUT/straddle$n.gfa" -r REF -o "$OUT/straddle$n" -q >/dev/null 2>&1
  got=$(tail -n +2 "$OUT/straddle$n.bubbles.csv" | cut -d, -f2,4 | paste -sd' ' -)
  [ "$got" = "$want" ] && ok "a fusion across a wrapper's edge with a $n-node bypass: $got" \
                       || bad "straddle $n: got $got, expected $want"
done
# A retained site off the reference, nested in an expansion's allele (hnest inserts 70 bp inside it),
# has no reference span; it still shares nodes with the fused site the AC expansion would make, and,
# smaller, would drop it with the 80 bp expansion beside it. The fusion is not made: 5,8 and 9,12 stay.
"${PYTHON:-python3}" - "$OUT/nest.gfa" 1 <<'PYEOF'
import random, sys
R = random.Random(7); r = lambda n: ''.join(R.choice('ACGT') for _ in range(n))
ac = lambda n: 'AC' * n
nest = sys.argv[2] == '1'
S = {'f1': r(300), 'n2': r(40)+ac(10), 'n3': 'AC', 'n5': ac(10), 'n6': ac(5),
     'n7a': ac(20), 'n7b': ac(20), 'n8': ac(10), 'n9': ac(5), 'n10': ac(40), 'n11': ac(10)+r(60), 'f2': r(300)}
if nest: S['n7x'] = r(70)
L = [('f1','n2'),('n2','n3'),('n3','n5'),('n2','n5'),('n5','n6'),('n5','n7a'),('n7a','n7b'),('n6','n8'),('n7b','n8'),
     ('n8','n9'),('n8','n10'),('n9','n11'),('n10','n11'),('n11','f2')]
if nest: L += [('n7a','n7x'),('n7x','n7b')]
P = {'REF': 'f1 n2 n3 n5 n6 n8 n9 n11 f2', 'hdel': 'f1 n2 n5 n6 n8 n9 n11 f2',
     'hexp': 'f1 n2 n3 n5 n7a n7b n8 n9 n11 f2', 'hk': 'f1 n2 n3 n5 n6 n8 n10 n11 f2',
     'hexp2': 'f1 n2 n3 n5 n7a n7b n8 n10 n11 f2'}
if nest: P['hnest'] = 'f1 n2 n3 n5 n7a n7x n7b n8 n9 n11 f2'
with open(sys.argv[1], 'w') as o:
    o.write('H\tVN:Z:1.0\n')
    for k, v in S.items(): o.write(f'S\t{k}\t{v}\n')
    for a, b in L: o.write(f'L\t{a}\t+\t{b}\t+\t0M\n')
    for k, v in P.items(): o.write(f'P\t{k}\t' + ','.join(x + '+' for x in v.split()) + '\t*\n')
PYEOF
"$BIN" bubble -i "$OUT/nest.gfa" -r REF -o "$OUT/nest0" --merge-nearby-bp 0 --emit-snarls-jsonl "$OUT/nest0.jsonl" -q >/dev/null 2>&1
{ cat "$OUT/nest0.jsonl"; echo '{"start": {"node_id": "5"}, "end": {"node_id": "8"}}'; } > "$OUT/nest.jsonl"
"$BIN" bubble -i "$OUT/nest0.sorted.gfa" -r REF --snarls "$OUT/nest.jsonl" -o "$OUT/nest" --merge-nearby-bp 0 -q >/dev/null 2>&1
[ "$(tail -n +2 "$OUT/nest.bubbles.csv" | cut -d, -f2,4 | paste -sd' ' -)" = "5,8 9,12" ] \
  && ok "a nested site off the reference blocks a fusion that would drop its neighbour (5,8 9,12)" \
  || bad "nested site: got $(tail -n +2 "$OUT/nest.bubbles.csv" | cut -d, -f2,4 | paste -sd' ' -), expected 5,8 9,12"
# Two CA microsatellites, each a 2 bp deletion site under the floor beside an expansion site. With
# --max-variant-bp 100, locus A fuses into an 82 bp site; locus B's fused site (112 bp) fails the cap and
# is undone on its own, its 90 bp expansion site kept. Undoing every fusion at the locus lost A's too.
"${PYTHON:-python3}" - "$OUT/undo.gfa" <<'PYEOF2'
import random, sys
r = random.Random(53)
rnd = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
ac = lambda k: 'AC' * k
S, L, P = {}, [], {'ref': [], 'hdel': [], 'hexp': []}
def node(seq):
    k = str(len(S) + 1); S[k] = seq; return k
def locus(alt_ac):
    n2, n3, n5, n6, n7, n8 = node(ac(10)), node('AC'), node(ac(10)), node(ac(5)), node(ac(alt_ac)), node(ac(10))
    L.extend([(n2, n3), (n3, n5), (n2, n5), (n5, n6), (n5, n7), (n6, n8), (n7, n8)])
    P['ref'] += [n2, n3, n5, n6, n8]; P['hdel'] += [n2, n5, n6, n8]; P['hexp'] += [n2, n3, n5, n7, n8]
    return n2, n8
f1 = node(rnd(300)); a2, a8 = locus(30); m = node(rnd(300)); b2, b8 = locus(45); f2 = node(rnd(300))
L.extend([(f1, a2), (a8, m), (m, b2), (b8, f2)])
for k in P: P[k] = [f1] + P[k][:5 if k == 'ref' else len(P[k]) // 2] + [m] + P[k][5 if k == 'ref' else len(P[k]) // 2:] + [f2]
with open(sys.argv[1], 'w') as o:
    o.write('H\tVN:Z:1.0\n')
    for k, v in S.items(): o.write(f'S\t{k}\t{v}\n')
    for a, b in L: o.write(f'L\t{a}\t+\t{b}\t+\t0M\n')
    for k, v in P.items(): o.write(f'P\t{k}\t' + ','.join(x + '+' for x in v) + '\t*\n')
PYEOF2
"$BIN" bubble -i "$OUT/undo.gfa" -r ref -o "$OUT/undo" --merge-nearby-bp 0 --max-variant-bp 100 -q >/dev/null 2>&1
[ "$(cel_all "$OUT/undo.bubbles.csv" max_inside_bp)" = "82 90" ] \
  && ok "a fused site failing a filter is undone on its own; the other fusion stands (max_inside 82 90)" \
  || bad "undo: max_inside_bp $(cel_all "$OUT/undo.bubbles.csv" max_inside_bp), expected 82 90"

# ---------------------------------------------------------------- a site in a cycle-rooted component
# A cycle (hdup walks 7 -> 1) roots the decomposition in a cycle. Most haplotypes end inside site 4..7,
# which REF, h4 and hdup cross. Holding most path ends does not make it a locus wrapper: it is reported,
# as in the same graph without the cycle (both 2..4 and 4..7).
cat > "$OUT/g3cyc.gfa" <<'GFAEOF'
H	VN:Z:1.0
S	1	CAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGTGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACCGCGGTGTTAAGTGTCGAGCTACATCACTTCTCATGTAGCCAGAAGGCTGCAACTCATCGACTCTATGTAGTGACCGCGTCGATGTCAAACCCCGGGGGGAGCTCAGATATCCGATACAGGGATGAAGAAATAACCTCATCCCATTGGTGACGAAAGGTTGTAAGTAGCTGGCCGCCGAGA
S	2	TAGCTGAGCGGCGAACCACTAGAAAAGGTTCAGACCCCGGAGCCCAGCCGTCACGATTGTTATGCGTATAAGCCCGGTTCACTACGTCCGTTCTGGCAAGCCGGGGCTAATCCGTCATTGTCAAGAGACATCTTTCGTCTCATTAGGCTACTAACGCCGCCGGGTCGTTACTCGAAAAGCAGGTGGAATTGGTGTATTCA
S	3	GCTTGCTCGATTTGATCGATCTGCAAGGTGCTGTCTAGATAGATACCATGGCCCGGAAGTACGGGCTTCTGGCGCATGTCGCACTCGTCCCTGGTCACGAACTGTACAAACATTGGACAC
S	4	TCTTTCCCGTTCTGGTACAAAATGTGCTCCAATCATGCATGAAACAGATACATCGCTTGGGCCACGTAGTCTAGAGCACACTAAATGAGACATCTTAGAGGAGATAGGCGTAGATCCGGTTACTAGCCGTGATGCAAGGTGGGGGAACGGGATGTTGTAACATGCGGGTGTGCACGCCACTAAGACGAAACCTAGTGCCTCTTGCTAGTCATTATTAGTACGAAGGGTTGTGCTCCGATAGTTGAAAATGTGGTGTTATGCTCACGGCGTGGTGTGTCTTTAACCCCAAGCTATCAATAC
S	5	TGAATAGGCTACATATGTTATACTCCGTGTCGTAAGGATGACGGCTCCGCTACTGGTGGTCTGTCGCCTCAGCCGTTGACCGCAACACCGTGAAGCACGGGTAAGGCAGCAGAAAGGCGAGAACTGCAGGAGAGCGTATTTGCGCAACCC
S	6	TGAGGGTCTA
S	7	GAGAGTCCACCTGGGCCTTTACGGAACTATATTGGTTTAATAAAACGGGTCCAGCAAGTGGATTTGGGTCCAGACTGAATCTCTCACGGCTTGTCTTTATGCCATTAAACTTGCCAGATTCTACTCCGCACCTACTCACACTTAATAATACAAGTGTCCGTTCTTCTGGCGGCAGGCGGGGTGTACCGCCACTCCTTCAACAATTTCCACTCGCTGCCGCGTGAGCTAGAGTGAAGCCAATCCTACTCGAACTTCGACCTGTTGTACCATATCTGCAAATTCCCTGCCGAGATACCGTAA
L	1	+	2	+	0M
L	2	+	3	+	0M
L	3	+	4	+	0M
L	2	+	4	+	0M
L	4	+	5	+	0M
L	4	+	6	+	0M
L	5	+	7	+	0M
L	6	+	7	+	0M
L	7	+	1	+	0M
P	REF	1+,2+,3+,4+,5+,7+	*
P	h4	1+,2+,4+,6+,7+	*
P	h1	1+,2+,4+,6+	*
P	h2	1+,2+,3+,4+,6+	*
P	h3	1+,2+,4+,5+	*
P	h5	1+,2+,3+,4+,6+	*
P	hdup	1+,2+,3+,4+,5+,7+,1+,2+,4+,6+,7+	*
GFAEOF
"$BIN" bubble -i "$OUT/g3cyc.gfa" -r REF -o "$OUT/g3cyc" --min-variant-bp 0 --merge-nearby-bp 0 --no-fuse-ambiguous -q >/dev/null 2>&1
[ "$(tail -n +2 "$OUT/g3cyc.bubbles.csv" | cut -d, -f2,4 | paste -sd' ' -)" = "2,4 4,7" ] \
  && ok "a crossed site where most haplotypes end is kept in a cycle-rooted component (2,4 4,7)" \
  || bad "cycle-rooted site: got $(tail -n +2 "$OUT/g3cyc.bubbles.csv" | cut -d, -f2,4 | paste -sd' ' -), expected 2,4 4,7"

# ---------------------------------------------------------------- shared input/output contract
# These go through `bubble` because it reads and writes with the shared code every module uses.
io_gfa() {   # io_gfa <out> [extra lines...]: a three-segment deletion site, plus any extra lines
  { printf 'H\tVN:Z:1.0\nS\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGG\n'
    printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
    printf 'P\tfull\t1+,2+,3+\t*\nP\tdel\t1+,3+\t*\n'
    shift; for l in "$@"; do printf '%s\n' "$l"; done; } > "$1"
}
io_refused() {   # io_refused <label> <gfa>: bubble must exit non-zero on it
  "$BIN" bubble -i "$2" -r full -o "$OUT/io_r" --min-variant-bp 0 -q >/dev/null 2>&1
  [ "$?" -ne 0 ] && ok "$1" || bad "$1 (accepted)"
}
io_gfa "$OUT/io_dupS.gfa" "$(printf 'S\t2\tTTTTTTTT')"
io_refused "a segment defined twice is refused (not resolved by file order)" "$OUT/io_dupS.gfa"
io_gfa "$OUT/io_ov.gfa"; sed -i.bak 's/0M$/0M2I/' "$OUT/io_ov.gfa"
io_refused "an overlap other than '*' or <n>M is refused (0M2I read as 0M)" "$OUT/io_ov.gfa"
io_gfa "$OUT/io_lx.gfa"; sed -i.bak '6s/\t+\t3/\tx\t3/' "$OUT/io_lx.gfa"
io_refused "a link orientation other than '+'/'-' is refused" "$OUT/io_lx.gfa"
{ printf 'H\tVN:Z:1.0\nS\t1\tAAAAAAAAAA\nS\t3\tGGGGGGGGGG\nL\t1\t+\t3\t+\t0M\n'
  printf 'W\tfull\t0\tchr1\t0\t20\t9>1>3\n'; } > "$OUT/io_w.gfa"
io_refused "a W walk with text before its first step is refused" "$OUT/io_w.gfa"

# A damaged member in a multi-member gzip ended the file there with a clean EOF: a shorter graph, exit 0.
io_gfa "$OUT/io_m.gfa"
# Every line but the last path in the first member, so a reader that stops there still has a graph.
sed '$d' "$OUT/io_m.gfa" | gzip -c > "$OUT/io_p1.gz"; tail -1 "$OUT/io_m.gfa" | gzip -c > "$OUT/io_p2.gz"
cat "$OUT/io_p1.gz" "$OUT/io_p2.gz" > "$OUT/io_multi.gfa.gz"
"$BIN" bubble -i "$OUT/io_multi.gfa.gz" -r full -o "$OUT/io_mo" --min-variant-bp 0 -q >/dev/null 2>&1
[ "$(grep -c '^P' "$OUT/io_mo.sorted.gfa" 2>/dev/null)" = "2" ] \
  && ok "a multi-member gzip is read through every member" || bad "a multi-member gzip lost paths"
python3 -c "import sys; b=bytearray(open(sys.argv[1],'rb').read()); b[0]=0x00; open(sys.argv[2],'wb').write(bytes(b))" \
  "$OUT/io_p2.gz" "$OUT/io_p2bad.gz"
cat "$OUT/io_p1.gz" "$OUT/io_p2bad.gz" > "$OUT/io_multibad.gfa.gz"
io_refused "a damaged gzip member is refused, not read as the end of the file" "$OUT/io_multibad.gfa.gz"

# A write that fails (a full disk; a file-size limit stands in) left a truncated CSV installed, exit 0.
io_gfa "$OUT/io_big.gfa"
awk 'BEGIN{for(i=0;i<4000;i++) printf "P\tp%d\t1+,2+,3+\t*\n", i}' >> "$OUT/io_big.gfa"
( trap '' XFSZ; ulimit -f 20; exec "$BIN" bubble -i "$OUT/io_big.gfa" -r full -o "$OUT/io_fz" \
    --min-variant-bp 0 -q ) >/dev/null 2>&1
rc=$?
[ "$rc" -ne 0 ] && [ ! -e "$OUT/io_fz.sorted.gfa" ] \
  && ok "a failed write is an error and installs nothing" \
  || bad "a failed write: exit $rc, sorted GFA $( [ -e "$OUT/io_fz.sorted.gfa" ] && echo installed || echo absent)"

# Outputs that differ only in case are one file on macOS; a directory at an output path is not moved aside.
io_gfa "$OUT/io_c.gfa"
mkdir -p "$OUT/io_case" "$OUT/io_dir/results"; echo keep > "$OUT/io_dir/results/keep.txt"
"$BIN" bubble -i "$OUT/io_c.gfa" -r full -o "$OUT/io_case/x" --min-variant-bp 0 \
  --sorted-gfa-out "$OUT/io_case/X.BUBBLES.CSV" -q >/dev/null 2>&1
rc=$?
if [ "$(uname)" = Darwin ]; then
  [ "$rc" -ne 0 ] && ok "two outputs differing only in case are refused on macOS" \
                 || bad "two outputs differing only in case were accepted"
fi
"$BIN" bubble -i "$OUT/io_c.gfa" -r full -o "$OUT/io_dir/y" --min-variant-bp 0 \
  --sorted-gfa-out "$OUT/io_dir/results" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && [ -f "$OUT/io_dir/results/keep.txt" ] \
  && ok "an output path naming an existing directory is refused and the directory kept" \
  || bad "a directory at an output path was not refused or not kept"

# ---- fixtures in synthetic_data/bubble_cases --------------------------------------------------------
CASES="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/synthetic_data/bubble_cases"
# allsites <csv>: every row's source..sink, space-joined
allsites() { tail -n +2 "$1" 2>/dev/null | awk -F, '{printf "%s%s..%s", (NR>1?" ":""), $2, $4}'; }

# The wrapper rule reads "most paths" over the snarl's own component. Counted over the file, 30 paths
# of an unrelated two-node component would outvote the locus's 14 and make its two sites one wrapper.
"$BIN" bubble -i "$CASES/wrapper_two_components.gfa" -r REF -o "$OUT/wc" --min-variant-bp 0 \
  --merge-nearby-bp 0 --no-fuse-ambiguous -q >/dev/null 2>&1
[ "$(allsites "$OUT/wc.bubbles.csv")" = "5..7 7..10" ] \
  && ok "an unrelated component's paths do not change a cycle-rooted locus's sites" \
  || bad "with an unrelated component the sites are '$(allsites "$OUT/wc.bubbles.csv")', expected '5..7 7..10'"

# A --merge-nearby-bp fusion that the overlap pass would drop takes a site with it: here the fused site
# meets the wrapper kept beside the sites inside it, and the site on the wrapper's far side must survive.
"$BIN" bubble -i "$CASES/merge_wrapper.gfa" -r REF -o "$OUT/mw" -q >/dev/null 2>"$OUT/mw.err"
[ "$(allsites "$OUT/mw.bubbles.csv")" = "3..6 8..12" ] \
  && ok "a merge that would lose a site to the overlap pass is undone, both sites kept" \
  || bad "merge next to a wrapper gave '$(allsites "$OUT/mw.bubbles.csv")', expected '3..6 8..12'"
grep -q "would lose a site to the overlap" "$OUT/mw.err" \
  && ok "the undone merge is reported" || bad "the undone merge was not reported"

# Only a crossing read through both inward sides is one: `3+ 1+` past the 1+..3+ site walks the outer
# sides of both boundaries and counted as a deletion allele.
"$BIN" bubble -i "$CASES/outward_walk.gfa" -r ref --snarls "$CASES/outward_walk.jsonl" -o "$OUT/ow" \
  --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
[ "$(cel "$OUT/ow.bubbles.csv" path_support)" = "1" ] && [ "$(cel "$OUT/ow.bubbles.csv" distinct_alleles)" = "1" ] \
  && ok "a walk through the boundaries' outer sides is not a crossing" \
  || bad "outward walk: path_support $(cel "$OUT/ow.bubbles.csv" path_support), distinct_alleles $(cel "$OUT/ow.bubbles.csv" distinct_alleles), expected 1 and 1"

# A path crossing an empty-interior site twice (1 3 1 3) carries two copies of the deletion, not the
# deletion: three distinct walks, best alternate support 1, so --min-alt-support 2 removes the site.
"$BIN" bubble -i "$CASES/repeated_empty.gfa" -r ref --snarls "$CASES/repeated_empty.jsonl" -o "$OUT/re" \
  --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
[ "$(cel "$OUT/re.bubbles.csv" distinct_alleles)" = "3" ] && [ "$(cel "$OUT/re.bubbles.csv" alt_allele_support_max)" = "1" ] \
  && ok "two empty crossings are one allele spanning both" \
  || bad "repeated empty crossings: distinct $(cel "$OUT/re.bubbles.csv" distinct_alleles), alt max $(cel "$OUT/re.bubbles.csv" alt_allele_support_max), expected 3 and 1"
# Numeric node ids are the same ids as quoted ones.
printf '{"start": {"node_id": 1}, "end": {"node_id": 3}}\n' > "$OUT/num.jsonl"
"$BIN" bubble -i "$CASES/repeated_empty.gfa" -r ref --snarls "$OUT/num.jsonl" -o "$OUT/num" \
  --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
[ "$(rows "$OUT/num.bubbles.csv")" = "1" ] \
  && ok "--snarls reads unquoted numeric node ids" || bad "numeric node ids gave $(rows "$OUT/num.bubbles.csv") bubbles, expected 1"

# Emitted snarls carry their boundary handles' orientation, so re-importing them gives the same sites.
# Voting the side from the walks instead would, on this 3-node cycle, bring back a site the internal run
# does not have (1..3 holding 2).
"$BIN" bubble -i "$CASES/outward_walk.gfa" -r ref -o "$OUT/rt" --min-variant-bp 0 --merge-nearby-bp 0 \
  --no-fuse-ambiguous --emit-snarls-jsonl "$OUT/rt.jsonl" -q >/dev/null 2>&1
"$BIN" bubble -i "$OUT/rt.sorted.gfa" -r ref --snarls "$OUT/rt.jsonl" -o "$OUT/rti" --min-variant-bp 0 \
  --merge-nearby-bp 0 --no-fuse-ambiguous -q >/dev/null 2>&1
grep -q '"backward"' "$OUT/rt.jsonl" && [ "$(allsites "$OUT/rt.bubbles.csv")" = "$(allsites "$OUT/rti.bubbles.csv")" ] \
  && ok "emitted snarls re-import to the same sites" \
  || bad "round trip: internal '$(allsites "$OUT/rt.bubbles.csv")', re-imported '$(allsites "$OUT/rti.bubbles.csv")'"

# An internal run that finds no snarl still writes the JSONL it was asked for, replacing an old one.
echo stale > "$OUT/ns.jsonl"
"$BIN" bubble -i "$CASES/no_snarl.gfa" -r ref -o "$OUT/ns" --emit-snarls-jsonl "$OUT/ns.jsonl" -q \
  >/dev/null 2>"$OUT/ns.err"
[ -f "$OUT/ns.jsonl" ] && [ ! -s "$OUT/ns.jsonl" ] && ! grep -q "no-op" "$OUT/ns.err" \
  && ok "--emit-snarls-jsonl with no snarl writes an empty file over the old one" \
  || bad "--emit-snarls-jsonl with no snarl left '$(head -c 20 "$OUT/ns.jsonl")' / $(head -c 80 "$OUT/ns.err")"

# Two snarls on the same two nodes, on different sides: only one can be a site (a site is named by its
# boundary nodes), and the run says the other was not reported.
"$BIN" bubble -i "$CASES/shared_pair.gfa" -r REF -o "$OUT/sp" --merge-nearby-bp 0 --no-fuse-ambiguous \
  -q >/dev/null 2>"$OUT/sp.err"
grep -q "sharing both boundary nodes" "$OUT/sp.err" \
  && ok "a snarl collapsed onto another's boundary nodes is reported" \
  || bad "the collapsed snarl was not reported"
# The two snarls lie inside a loop h2 and h3 walk twice (X a Y b X a Y b), so the site is widened to the
# nodes around the loop, W..Z, and holds both snarls' interiors; `call` places its alleles.
[ "$(awk -F, 'NR==2{print $2 $3 $4 $5, $15, $10}' "$OUT/sp.bubbles.csv")" = '1+8+ "2;3;4;5;6;7" 1' ] \
  && "$BIN" call -i "$OUT/sp.sorted.gfa" -b "$OUT/sp" -r REF -o "$OUT/spc" -q >/dev/null 2>&1 \
  && [ "$(grep -vc '^#' "$OUT/spc.region.vcf")" -gt 0 ] \
  && ok "two snarls on one node pair inside a loop become one site around it, and call places it" \
  || bad "shared pair row: $(awk -F, 'NR==2{print $2 $3 $4 $5, $15, $10}' "$OUT/sp.bubbles.csv"), records $(grep -vc '^#' "$OUT/spc.region.vcf" 2>/dev/null)"

# An undefined node id is refused before the sort renumbers it: h1 10+,2+,12+ has no segment 2, and once
# renumbered its 2 would be the reference's middle node, so h1 would read as the reference.
"$BIN" bubble -i "$CASES/undefined_node.gfa" -r ref -o "$OUT/und" -q >/dev/null 2>"$OUT/und.err"; rc=$?
[ "$rc" -ne 0 ] && grep -q "not in the graph: 2" "$OUT/und.err" \
  && ok "a walk through an undefined node is refused before renumbering" \
  || bad "undefined node: exit $rc, $(head -c 120 "$OUT/und.err")"

# The stale-file sweep never removes a file this run wrote: a site table written as <prefix>.bandage_genes.csv.
"$BIN" bubble -i "$CASES/shared_pair.gfa" -r REF -o "$OUT/own" --bubbles-csv "$OUT/own.bandage_genes.csv" -q >/dev/null 2>&1
[ -s "$OUT/own.bandage_genes.csv" ] && head -1 "$OUT/own.bandage_genes.csv" | grep -q '^bubble_id,' \
  && ok "the stale sweep keeps this run's own output" \
  || bad "own output named like a stale file was removed"

# --no-flip still numbers along the reference: it reads some nodes reversed, and taking their stored
# orientation as the walking direction numbered them against it (1+,4-,3-,5+).
"$BIN" bubble -i "$CASES/noflip.gfa" -r ref -o "$OUT/nf" --min-variant-bp 0 --no-flip -q >/dev/null 2>&1
[ "$(awk -F'\t' '$1=="P" && $2=="ref"{print $3}' "$OUT/nf.sorted.gfa")" = "1+,2-,4-,5+" ] \
  && ok "--no-flip keeps node ids increasing along the reference" \
  || bad "--no-flip reference is $(awk -F'\t' '$1=="P" && $2=="ref"{print $3}' "$OUT/nf.sorted.gfa"), expected 1+,2-,4-,5+"

# With --snarls no sorted graph is written, so one an earlier run left at this prefix is removed;
# likewise a gene table this run did not write.
printf 'stale\n' > "$OUT/sti.sorted.gfa"; printf 'stale\n' > "$OUT/sti.bandage_genes.csv"
"$BIN" bubble -i "$CASES/repeated_empty.gfa" -r ref --snarls "$CASES/repeated_empty.jsonl" -o "$OUT/sti" \
  --min-variant-bp 0 -q >/dev/null 2>&1
[ ! -e "$OUT/sti.sorted.gfa" ] && [ ! -e "$OUT/sti.bandage_genes.csv" ] \
  && ok "default-named outputs this run did not write are removed" \
  || bad "a stale sorted.gfa or bandage_genes.csv survived"
# A graph read through /dev/stdin may be that sorted graph itself, which no identity test can see: it is kept.
"$BIN" bubble -i "$CASES/repeated_empty.gfa" -r ref -o "$OUT/sin" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" bubble -i /dev/stdin -r ref --snarls "$CASES/repeated_empty.jsonl" -o "$OUT/sin" --min-variant-bp 0 -q \
  < "$OUT/sin.sorted.gfa" >/dev/null 2>&1
[ "$?" = 0 ] && [ -s "$OUT/sin.sorted.gfa" ] && ok "  ... but not a sorted graph read through /dev/stdin" \
  || bad "the sorted graph read through /dev/stdin was removed"
# A folded graph re-run at its own panphorte prefix keeps the files panphorte wrote there: its
# unit-graph table describes that graph by node id, which re-snarling does not change.
FD="$OUT/fold"; mkdir -p "$FD"
PCASE="$(dirname "$CASES")/panphorte_cases/case_period.gfa"
"$BIN" bubble -i "$PCASE" -r arr -o "$FD/b" -q >/dev/null 2>&1
"$BIN" panphorte -i "$FD/b.sorted.gfa" -b "$FD/b" -r arr -o "$FD/P" --min-similarity 0.95 -q >/dev/null 2>&1
fold_kept() { for x in .normalized.sorted.gfa .panphorte.unit_graph.tsv .panphorte.report.tsv; do
                [ -s "$FD/P$x" ] || printf '%s ' "$x"; done; }
if [ -s "$FD/P.panphorte.unit_graph.tsv" ]; then
  "$BIN" bubble -i "$FD/P.normalized.sorted.gfa" -r arr -o "$FD/P" -q >/dev/null 2>&1
  [ "$?" = 0 ] && [ -z "$(fold_kept)" ] && ok "a re-run at a folded graph's prefix keeps its unit-graph table" \
    || bad "re-run at the folded graph's prefix lost: $(fold_kept)"
  "$BIN" bubble -i /dev/stdin -r arr -o "$FD/P" -q < "$FD/P.normalized.sorted.gfa" >/dev/null 2>&1
  [ "$?" = 0 ] && [ -z "$(fold_kept)" ] && ok "  ... and the folded graph read through /dev/stdin" \
    || bad "re-run from /dev/stdin at the folded graph's prefix lost: $(fold_kept)"
else
  bad "panphorte wrote no unit-graph table for the folded-prefix test"
fi

# An output directory that cannot be written gives the system's reason, not a disk-full hint.
mkdir -p "$OUT/rodir"; chmod 555 "$OUT/rodir"
if [ ! -w "$OUT/rodir" ]; then
  ro_err=$("$BIN" bubble -i "$PCASE" -r arr -o "$OUT/rodir/x" -q 2>&1 >/dev/null)
  case "$ro_err" in *"rodir/x.sorted.gfa: Permission denied"*) ok "an unwritable output names the file and the reason";;
    *) bad "unwritable output message: $ro_err";; esac
fi
chmod 755 "$OUT/rodir"

# ---- shared GFA / gzip input --------------------------------------------------------------------------
IO="$OUT/io"; mkdir -p "$IO"
check() { [ "$2" = "$3" ] && ok "$1 = $3" || bad "$1: expected $3, got ${2:-<missing>}"; }
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tCAGATTTTCATATTATGCAGAAAATCTACTTCGCCTGATACGAGTCGGTTATCTTCGGATACTGTATAGTCCCACCTGGT\n'
  printf 'S\t2\tA\nS\t3\tG\n'
  printf 'S\t4\tGATCCTATGCTTGTGAGTACCCAGAAAATAGCGACGGACCGCGGTGTTAAGTGTCGAGCT\n'
  printf 'S\t5\tACATCACTTCTCATGTAGCCAGAAGGCTGCAACTCATCGACTCTATGTAGTGACCGCGTCGATGTCAAACCCCGGGGGGA\n'
  printf 'S\t6\tGCTCAGATATCCGATACAGGGATGAAGAAATAACCTCATCCCATTGGTGACGAAAGGTTGTAAGTAGCTG\n'
  for l in "1 2" "1 3" "2 4" "3 4" "4 5" "4 6" "6 5" "2 5"; do set -- $l; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
} > "$IO/head.gfa"
wgfa() {   # wgfa <ref start> <ref end>: the graph with W lines, the reference walking 221 bp
  cat "$IO/head.gfa"
  printf 'W\tGRCh38\t0\tchr1\t%s\t%s\t>1>2>4>5\n' "$1" "$2"
  printf 'W\tHG1\t1\tchr1\t*\t*\t>1>3>4>6>5\nW\tHG2\t1\tchr1\t0\t141\t>1>2>5\nW\tHG3\t1\tchr1\t*\t*\t>1>3>4>5\n'
}
wpos() {   # wpos <gfa> <label>: POS of the first record call writes on it
  "$BIN" bubble -i "$1" -r GRCh38 -o "$IO/$2" -q >/dev/null 2>&1
  "$BIN" call -i "$IO/$2.sorted.gfa" -c "$IO/$2.bubbles.csv" -r GRCh38 -o "$IO/$2_c" --no-per-bubble-vcf -q >/dev/null 2>&1
  grep -v '^#' "$IO/$2_c.region.vcf" 2>/dev/null | head -1 | cut -f1,2
}
wgfa 1000 1221 > "$IO/w.gfa"; wgfa 1000 1220 > "$IO/wm.gfa"
# A W line's start is 0-based, so its walk starts at 1001 and the deletion of node 2 anchors at 1081, also
# when the interval is one base short of the walk (the shifted name must not be read as 0-based again).
check "a W-line reference: POS from its 0-based start, with or without a length mismatch" \
      "$(wpos "$IO/w.gfa" w) $(wpos "$IO/wm.gfa" wm)" "$(printf 'chr1\t1081 chr1\t1081')"
"$BIN" bubble -i "$IO/w.gfa" -r 'GRCh38#0#chr1:1000-1221' -o "$IO/wl" -q >/dev/null 2>&1
check "-r may name a W line by its own coordinates" "$?" "0"
# A gzip graph through /dev/stdin from a file: the BGZF check must not reopen the path and move the
# shared offset, or the reader decompresses from the wrong place.
# (over the reader's 64 KiB buffer, so a second read happens at the moved offset: an unconnected 300 kb
# segment of incompressible sequence)
{ cat "$IO/w.gfa"; printf 'S\t99\t'; awk 'BEGIN{srand(7); for(i=0;i<300000;i++) printf "%s", substr("ACGT", int(rand()*4)+1, 1); print ""}'; } \
  | gzip -c > "$IO/w.gfa.gz"
"$BIN" bubble -i /dev/stdin -r GRCh38 -o "$IO/stdin" -q < "$IO/w.gfa.gz" >/dev/null 2>&1
check "a gzip graph read through /dev/stdin redirected from a file" "$?" "0"
# CRLF line ends are read as LF: no sequence, walk step or overlap field keeps the carriage return.
sed 's/$/\r/' "$IO/w.gfa" > "$IO/wcr.gfa"
"$BIN" bubble -i "$IO/w.gfa" -r GRCh38 -o "$IO/wlf" -q >/dev/null 2>&1
"$BIN" bubble -i "$IO/wcr.gfa" -r GRCh38 -o "$IO/wcr" -q >/dev/null 2>&1
check "a CRLF graph with W lines gives the same sites and sorted graph" \
  "$(cmp -s "$IO/wlf.bubbles.csv" "$IO/wcr.bubbles.csv" && cmp -s "$IO/wlf.sorted.gfa" "$IO/wcr.sorted.gfa" && echo same)" "same"
sed 's/$/\r/' "$OUT/del.gfa" > "$IO/delcr.gfa"
"$BIN" bubble -i "$IO/delcr.gfa" -r full -o "$IO/delcr" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" bubble -i "$OUT/del.gfa" -r full -o "$IO/dellf" --min-variant-bp 0 -q >/dev/null 2>&1
check "a CRLF graph with P lines gives the same sites and sorted graph" \
  "$(cmp -s "$IO/dellf.bubbles.csv" "$IO/delcr.bubbles.csv" && cmp -s "$IO/dellf.sorted.gfa" "$IO/delcr.sorted.gfa" && echo same)" "same"
# With --snarls the input itself is the graph read, not a sorted copy.
printf '{"start": {"node_id": "1"}, "end": {"node_id": "3"}}\n' > "$IO/del.jsonl"
"$BIN" bubble -i "$OUT/del.gfa" -r full --snarls "$IO/del.jsonl" -o "$IO/dellfs" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" bubble -i "$IO/delcr.gfa" -r full --snarls "$IO/del.jsonl" -o "$IO/delcrs" --min-variant-bp 0 -q >/dev/null 2>&1
check "a CRLF graph read with --snarls gives the same sites" \
  "$(cmp -s "$IO/dellfs.bubbles.csv" "$IO/delcrs.bubbles.csv" && echo same)" "same"
# Malformed input is refused, not crashed on or half-read.
{ cat "$IO/w.gfa"; printf 'S\t9\t*\n'; } > "$IO/star.gfa"
"$BIN" bubble -i "$IO/star.gfa" -r GRCh38 -o "$IO/star" -q >/dev/null 2>&1
check "a segment with no sequence is refused (the snarl finder has no bases to work on)" "$?" "2"
{ cat "$IO/head.gfa"; printf 'P\tGRCh38#0#chr1:1000-1221\t1+,2+,4+,5+\t0M,0M,5M\nP\th\t1+,3+,4+,5+\t*\n'; } > "$IO/povl.gfa"
"$BIN" bubble -i "$IO/povl.gfa" -r GRCh38 -o "$IO/povl" -q >/dev/null 2>&1
check "a P line with a non-zero overlap is refused" "$?" "2"
{ cat "$IO/head.gfa"; printf 'P\tGRCh38#0#chr1:1000-1221\t1+,,2+,4+,5+\t*\nP\th\t1+,3+,4+,5+\t*\n'; } > "$IO/pempty.gfa"
"$BIN" bubble -i "$IO/pempty.gfa" -r GRCh38 -o "$IO/pempty" -q >/dev/null 2>&1
check "a P line with an empty step is refused" "$?" "2"
# An insertion split across adjacent snarls, none reaching --min-variant-bp alone: the pieces are fused,
# and kept when an allele differs in length from the reference's by the floor. Three 20 bp insertions
# 10 bp apart = one 60 bp site.
# The same chain with equal-length alleles (substitutions) stays below the floor: no site.
frag() {   # frag <label> <alt piece length>
  local lb=$1 alt=$2
  { printf 'H\tVN:Z:1.0\n'
    printf 'S\tS\t%s\nS\tT\t%s\n' "$(awk 'BEGIN{srand(11); for(i=0;i<300;i++) printf "%s", substr("ACGT",int(rand()*4)+1,1)}')" \
                                     "$(awk 'BEGIN{srand(12); for(i=0;i<300;i++) printf "%s", substr("ACGT",int(rand()*4)+1,1)}')"
    for n in a b c; do
      printf 'S\t%s1\t%s\nS\t%s2\t%s\n' $n "$(awk -v s=$n 'BEGIN{srand(length(s)+20); for(i=0;i<5;i++) printf "%s", substr("ACGT",int(rand()*4)+1,1)}')" \
             $n "$(awk -v s=$n -v l=$alt 'BEGIN{srand(length(s)+30); for(i=0;i<l;i++) printf "%s", substr("ACGT",int(rand()*4)+1,1)}')"
    done
    printf 'S\tx\tGATTACAGAT\nS\ty\tCCGGTTAACC\n'
    for e in "S a1" "S a2" "a1 x" "a2 x" "x b1" "x b2" "b1 y" "b2 y" "y c1" "y c2" "c1 T" "c2 T"; do
      set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
    printf 'P\tGRCh38#0#chr1:0-640\tS+,a1+,x+,b1+,y+,c1+,T+\t*\n'
    printf 'P\th1#1#chr1\tS+,a1+,x+,b1+,y+,c1+,T+\t*\n'
    printf 'P\th2#1#chr1\tS+,a2+,x+,b2+,y+,c2+,T+\t*\n'; } > "$IO/$lb.gfa"
  "$BIN" bubble -i "$IO/$lb.gfa" -r GRCh38 -o "$IO/$lb" -q >/dev/null 2>&1 || { echo "bubble failed: exit $?"; return; }
  tail -n +2 "$IO/$lb.bubbles.csv" | wc -l | tr -d ' '
}
check "an insertion split across three snarls under the floor is one site" "$(frag fragins 25)" "1"
check "a chain of equal-length substitutions under the floor stays no site" "$(frag fragsub 5)" "0"
# The floor is measured from the reference's allele: a short and a long alternate 60 bp apart, each 30 bp
# from the reference, are no site (call would report nothing there).
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tS\t%s\nS\tT\t%s\n' "$(awk 'BEGIN{srand(11); for(i=0;i<300;i++) printf "%s", substr("ACGT",int(rand()*4)+1,1)}')" \
                                   "$(awk 'BEGIN{srand(12); for(i=0;i<300;i++) printf "%s", substr("ACGT",int(rand()*4)+1,1)}')"
  for n in a b c; do
    for k in 1:5 2:25 3:15; do
      printf 'S\t%s%s\t%s\n' $n ${k%:*} "$(awk -v s=$n$k -v l=${k#*:} 'BEGIN{srand(length(s)*7+l); for(i=0;i<l;i++) printf "%s", substr("ACGT",int(rand()*4)+1,1)}')"
    done
  done
  printf 'S\tx\tGATTACAGAT\nS\ty\tCCGGTTAACC\n'
  for k in 1 2 3; do
    for e in "S a$k" "a$k x" "x b$k" "b$k y" "y c$k" "c$k T"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  done
  printf 'P\tGRCh38#0#chr1:0-665\tS+,a3+,x+,b3+,y+,c3+,T+\t*\n'
  printf 'P\th1#1#chr1\tS+,a1+,x+,b1+,y+,c1+,T+\t*\n'
  printf 'P\th2#1#chr1\tS+,a2+,x+,b2+,y+,c2+,T+\t*\n'; } > "$IO/fragmid.gfa"
"$BIN" bubble -i "$IO/fragmid.gfa" -r GRCh38 -o "$IO/fragmid" -q >/dev/null 2>&1
check "alternates differing by the floor only from each other are no site" \
  "$?:$(tail -n +2 "$IO/fragmid.bubbles.csv" 2>/dev/null | wc -l | tr -d ' ')" "0:0"
# ---- inversion sites and the rescue on spec fixtures (tests/tools/spec_gfa.py) ----
SPEC_GFA="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/tools/spec_gfa.py"
spec_run() {   # spec_run <label> [bubble args] < spec: prints "id source..sink kind parent" per site
  local lb=$1; shift
  cat > "$IO/$lb.spec"
  python3 -I "$SPEC_GFA" "$IO/$lb.spec" "$IO/$lb" >/dev/null
  # A failed run prints that, so a check expecting no site cannot pass on it.
  "$BIN" bubble -i "$IO/$lb/g.gfa" -r ref -o "$IO/$lb/b" -q "$@" >/dev/null 2>"$IO/$lb/b.log" \
    || { echo "bubble failed: exit $?"; return; }
  tail -n +2 "$IO/$lb/b.bubbles.csv" | awk -F, '{printf "%s %s..%s %s %s;", $1, $2, $4, $16, $17}'
}
INV_NODES='N S 2000 1
N I 300 2
N M1 3000 3
N Ma 1500 4
N X 60 8
N Mb 1500 9
N M3 3000 5
N T 2000 6'
INV_PATHS='P ref#0#chr1 S+,I+,M1+,Ma+,X+,Mb+,M3+,I-,T+
P h1#1#chr1 S+,I+,M1+,Ma+,Mb+,M3+,I-,T+
P hI#1#chr1 S+,I+,M3-,Mb-,Ma-,M1-,I-,T+
P hJ#1#chr1 S+,I+,M3-,Mb-,X-,Ma-,M1-,I-,T+
P h2#1#chr1 S+,I+,M1+,Ma+,X+,Mb+,M3+,I-,T+'
# A 9.6 kb inversion between inverted repeats (I) has no snarl of its own: it gets an inversion site after the
# ordinary ones, and the deletion site inside names it as parent.
check "an inversion between collapsed inverted repeats gets a site enclosing the site inside it" \
  "$(printf '%s\n%s\n' "$INV_NODES" "$INV_PATHS" | spec_run inv)" "1 4..6 site 2;2 1..8 inversion ;"
check "inversion sites are made with the size floor off" \
  "$(printf '%s\n%s\n' "$INV_NODES" "$INV_PATHS" | spec_run inv0 --min-variant-bp 0 | grep -o 'inversion' | head -1)" "inversion"
# A partial haplotype does not stop a backbone node from being a flank, and a branch only it walks still
# belongs to the site inside.
check "a partial haplotype does not prevent the inversion site" \
  "$(printf '%s\n%s\nP hP#1#chr1 Ma+,X+,Mb+,M3+\n' "$INV_NODES" "$INV_PATHS" | spec_run invpart)" "1 4..6 site 2;2 1..8 inversion ;"
check "a site whose branch only a partial haplotype walks is still inside the inversion site" \
  "$(printf '%s\nN Y 70 21\n%s\nP hP#1#chr1 Ma+,Y+,Mb+\n' "$INV_NODES" "$INV_PATHS" | spec_run invbranch)" "1 4..7 site 2;2 1..9 inversion ;"
# A site whose boundaries are the inversion's first and last nodes does not hold it: carriers cross it
# reversed end to end, which reads as its reference allele.
check "an inversion reaching a site's boundaries still gets an inversion site" \
  "$(spec_run invbound <<'EOF'
N S 2000 1
N I 300 2
N M1 3000 3
N J 250 10
N Ma 1500 4
N X 60 8
N Mb 1500 9
N M3 3000 5
N T 2000 6
P ref#0#chr1 S+,I+,M1+,J+,Ma+,X+,Mb+,J-,M3+,I-,T+
P h1#1#chr1 S+,I+,M1+,J+,Ma+,Mb+,J-,M3+,I-,T+
P hI#1#chr1 S+,I+,M3-,J+,Mb-,Ma-,J-,M1-,I-,T+
P hJ#1#chr1 S+,I+,M3-,J+,Mb-,X-,Ma-,J-,M1-,I-,T+
P hN#1#chr1 S+,I+,M1+,J+,Mb-,X-,Ma-,J-,M3+,I-,T+
P h2#1#chr1 S+,I+,M1+,J+,Ma+,X+,Mb+,J-,M3+,I-,T+
EOF
)" "1 3..8 site 2;2 1..9 inversion ;"
# An inversion site passes the support filters and --max-variant-bp like any site (the size floor does not
# apply). Two haplotypes carry the 9.6 kb inversion, five cross it.
check "--max-variant-bp drops an inversion site over the cap" \
  "$(printf '%s\n%s\n' "$INV_NODES" "$INV_PATHS" | spec_run invcap --max-variant-bp 5000)" "1 4..6 site ;"
check "--min-alt-support drops an inversion site fewer haplotypes carry" \
  "$(printf '%s\n%s\n' "$INV_NODES" "$INV_PATHS" | spec_run invalt --min-alt-support 3)" ""
check "--min-path-support drops an inversion site fewer haplotypes cross" \
  "$(printf '%s\n%s\n' "$INV_NODES" "$INV_PATHS" | spec_run invps --min-path-support 6)" ""
check "an inversion site under the cap and the support floors is kept" \
  "$(printf '%s\n%s\n' "$INV_NODES" "$INV_PATHS" | spec_run invok --max-variant-bp 20000 --min-alt-support 2 --min-path-support 5)" \
  "1 4..6 site 2;2 1..8 inversion ;"
# A one-node inversion is an ordinary site; dropped by the cap, it does not come back as an inversion site.
check "a site dropped by --max-variant-bp is not made an inversion site" \
  "$(spec_run invnode --max-variant-bp 100 <<'EOF'
N A 100 1
N B 300 2
N C 100 3
P ref A+,B+,C+
P h A+,B-,C+
EOF
)" ""
# An inversion near the locus's ends keeps its short flanks: the 1 kb inversion between 30 bp end nodes has
# no snarl and is read from the walks.
check "an inversion with flanks shorter than the size floor gets a site" \
  "$(spec_run invshort <<'EOF'
N S 30 1
N Ma 500 4
N Mb 500 9
N T 30 6
P ref#0#chr1 S+,Ma+,Mb+,T+
P h1#1#chr1 S+,Ma+,Mb+,T+
P hI#1#chr1 S+,Mb-,Ma-,T+
EOF
)" "1 1..4 inversion ;"
# A carrier stored on the other strand whose first node is short and walked as the reference walks it (an
# inverted end node) is still read on the reference's strand: hR carries hI's inversion, stored reversed.
{ printf '%s\nN E 20 7\n' "$INV_NODES"; sed 's/,T+$/,T+,E+/' <<< "$INV_PATHS"
  printf 'P hR#1#chr1 E+,T-,I+,M1+,Ma+,Mb+,M3+,I-,S-\n'; } | spec_run invrev >/dev/null
check "a reverse-strand carrier with a short inverted end node supports the inversion site" \
  "$(awk -F, '$16=="inversion"{print $11 " of " $8}' "$IO/invrev/b.bubbles.csv")" "3 of 6"
# A carrier's strand comes from where its ends lie on the reference, so an end flank under the size floor
# reads the same on both strands: hR is hI stored reversed, with the short flank first.
INV_SHORT_T=$(sed 's/^N T 2000 6$/N T 30 6/' <<< "$INV_NODES")
INV_SHORT_S=$(sed 's/^N S 2000 1$/N S 30 1/' <<< "$INV_NODES")
INV_REF3='P ref#0#chr1 S+,I+,M1+,Ma+,X+,Mb+,M3+,I-,T+
P h1#1#chr1 S+,I+,M1+,Ma+,Mb+,M3+,I-,T+
P h2#1#chr1 S+,I+,M1+,Ma+,X+,Mb+,M3+,I-,T+'
INV_FWD='P hI#1#chr1 S+,I+,M3-,Mb-,X-,Ma-,M1-,I-,T+'
INV_REV='P hR#1#chr1 T-,I+,M1+,Ma+,X+,Mb+,M3+,I-,S-'
check "a reverse-stored carrier behind a short right flank gets the inversion site" \
  "$(printf '%s\n%s\n%s\n' "$INV_SHORT_T" "$INV_REF3" "$INV_REV" | spec_run invshr)" "1 4..6 site 2;2 1..8 inversion ;"
printf '%s\n%s\n%s\n%s\n' "$INV_SHORT_T" "$INV_REF3" "$INV_FWD" "$INV_REV" | spec_run invshb >/dev/null
check "  ... and counts as a carrier beside the same haplotype stored forward" \
  "$(awk -F, '$16=="inversion"{print $11 " of " $8}' "$IO/invshb/b.bubbles.csv")" "2 of 5"
check "a reverse-stored carrier behind a short left flank gets the inversion site" \
  "$(printf '%s\n%s\n%s\n' "$INV_SHORT_S" "$INV_REF3" "$INV_REV" | spec_run invshl)" "1 4..6 site 2;2 1..8 inversion ;"
# Two inversions on one haplotype, one beside each short end flank: the haplotype is still read forward,
# and stored reversed it carries both too. The two sites share only their boundary node M.
INV_TWO='N S 30 1
N I1 300 2
N A 2000 3
N M 5000 4
N I2 300 5
N B 2000 6
N T 30 7
P ref#0#chr1 S+,I1+,A+,I1-,M+,I2+,B+,I2-,T+
P h1#1#chr1 S+,I1+,A+,I1-,M+,I2+,B+,I2-,T+'
printf '%s\nP hX#1#chr1 S+,I1+,A-,I1-,M+,I2+,B-,I2-,T+\nP hR#1#chr1 T-,I2+,B+,I2-,M-,I1+,A+,I1-,S-\n' "$INV_TWO" \
  | spec_run invtwo >/dev/null
check "two inversions beside short end flanks on one haplotype, stored either way, get a site each" \
  "$(awk -F, '$16=="inversion"{printf "%s..%s %s;", $2, $4, $11}' "$IO/invtwo/b.bubbles.csv")" "1..4 2;4..7 2;"
check "inversion sites sharing only a boundary node are both kept" \
  "$(printf '%s\nP hX#1#chr1 S+,I1+,A-,I1-,M+,I2+,B+,I2-,T+\nP hY#1#chr1 S+,I1+,A+,I1-,M+,I2+,B-,I2-,T+\n' "$INV_TWO" \
     | spec_run invadj)" "1 1..4 inversion ;2 4..7 inversion ;"
# A short node inside the inversion walked as the reference walks it does not break the falling order: every
# step pair across it falls too. hI carries hJ's inversion with s left forward, hR is hI stored reversed.
spec_run invsn >/dev/null <<'EOF'
N S 2000 1
N I 300 2
N M1 3000 3
N Ma 1500 4
N s 20 11
N X 60 8
N Mb 1500 9
N M3 3000 5
N T 2000 6
P ref#0#chr1 S+,I+,M1+,Ma+,s+,X+,Mb+,M3+,I-,T+
P h1#1#chr1 S+,I+,M1+,Ma+,s+,Mb+,M3+,I-,T+
P hJ#1#chr1 S+,I+,M3-,Mb-,X-,s-,Ma-,M1-,I-,T+
P h2#1#chr1 S+,I+,M1+,Ma+,s+,X+,Mb+,M3+,I-,T+
P hI#1#chr1 S+,I+,M3-,Mb-,X-,s+,Ma-,M1-,I-,T+
P hR#1#chr1 T-,I+,M1+,Ma+,s-,X+,Mb+,M3+,I-,S-
EOF
check "a short forward node inside an inversion counts in its falling order on both strands" \
  "$(awk -F, '$16=="inversion"{print $11 " of " $8}' "$IO/invsn/b.bubbles.csv")" "3 of 6"
# Short nodes out of place after a carrier's last long run, absorbed into it, do not decide its strand: y and
# e0 sit in the wrong place at hF's end, and hR is hF stored reversed. Both carry the inversion.
printf '%s\nN e0 30 12\nN y 30 13\nN z 30 14\n%s\n' "$INV_NODES" 'P ref#0#chr1 e0+,S+,I+,M1+,Ma+,X+,Mb+,M3+,I-,T+,y+,z+
P h1#1#chr1 e0+,S+,I+,M1+,Ma+,Mb+,M3+,I-,T+,y+,z+
P h2#1#chr1 e0+,S+,I+,M1+,Ma+,X+,Mb+,M3+,I-,T+,y+,z+
P hF#1#chr1 S+,I+,M3-,Mb-,X-,Ma-,M1-,I-,T+,y-,e0+,z-
P hR#1#chr1 z+,e0-,y+,T-,I+,M1+,Ma+,X+,Mb+,M3+,I-,S-' | spec_run invabs >/dev/null
check "short nodes absorbed after a carrier's last long run do not decide its strand" \
  "$(awk -F, '$16=="inversion"{print $11 " of " $8}' "$IO/invabs/b.bubbles.csv")" "2 of 5"
# An inversion whose flank lies in a site's interior with no backbone node beyond is lost, and says so
# under -q: hB and hB2 invert everything between the end nodes, which the site 1..3 then holds.
spec_run invnf >/dev/null <<'EOF'
N E1 30 1
N I1 300 2
N A 2000 3
N M1 5000 4
N I3 300 5
N C 2000 6
N M2 5000 7
N I2 300 8
N B 2000 9
N E2 30 10
P ref#0#chr1 E1+,I1+,A+,I1-,M1+,I3+,C+,I3-,M2+,I2+,B+,I2-,E2+
P h1#1#chr1 E1+,I1+,A+,I1-,M1+,I3+,C+,I3-,M2+,I2+,B+,I2-,E2+
P hX#1#chr1 E1+,I1+,A-,I1-,M1+,I3+,C+,I3-,M2+,I2+,B+,I2-,E2+
P hB#1#chr1 E1+,I1+,B-,I2-,M2-,I3+,C-,I3-,M1-,I1+,A-,I2-,E2+
P hB2#1#chr1 E1+,I1+,B-,I2-,M2-,I3+,C-,I3-,M1-,I1+,A-,I2-,E2+
EOF
check "an inversion that cannot be made a site is a warning under -q" \
  "$(grep -c 'warning: inversion not made a site' "$IO/invnf/b.log")" "1"
# The rescue never chains across a retained site: a SNP before a 100 bp insertion site does not stop the two
# 30 bp pieces after it from being fused. The rescued site then merges with that site, 20 bp away.
check "a split insertion is rescued beside a retained site and merged with it" \
  "$(spec_run rsbar <<'EOF'
N b0 500 1
N b1 20 2
N s =A
N sx =C
N b2 20 3
N r 100 4
N b3 20 5
N i2 30 6
N b4 20 7
N i3 30 8
N b5 20 9
N b6 500 10
P ref b0+,b1+,s+,b2+,b3+,b4+,b5+,b6+
P hA b0+,b1+,sx+,b2+,r+,b3+,i2+,b4+,i3+,b5+,b6+
P hB b0+,b1+,s+,b2+,r+,b3+,b4+,b5+,b6+
P hC b0+,b1+,s+,b2+,b3+,i2+,b4+,i3+,b5+,b6+
EOF
)" "1 5..11 site ;"
# SNP snarls do not link small insertions far apart into one site.
check "small insertions far apart are not chained through SNPs" \
  "$(spec_run rschain <<'EOF'
N F 500 1
N b0 60 2
N b1 60 3
N b2 60 4
N b3 60 5
N b4 60 6
N G 500 7
N i0 30 8
N s1 =A
N t1 =C
N s2 =A
N t2 =C
N i3 30 9
P ref F+,b0+,b1+,s1+,b2+,s2+,b3+,b4+,G+
P hA F+,b0+,i0+,b1+,t1+,b2+,t2+,b3+,i3+,b4+,G+
P hB F+,b0+,b1+,t1+,b2+,s2+,b3+,b4+,G+
EOF
)" ""
# panphorte reads its graph twice: a pipe is refused rather than read as empty the second time.
"$BIN" panphorte -i <(cat "$IO/w.sorted.gfa") --bubble-prefix "$IO/w" -o "$IO/pp" --reference-path GRCh38 -q >/dev/null 2>&1
check "panphorte refuses a graph from a pipe (it reads it twice)" "$?" "2"
check "--version prints the name and a version" "$("$BIN" --version | grep -c '^panvar [0-9][0-9.]*$')" "1"

# ---- loops, merge chains and --snarls order ----
# A tandem copy stored as a loop of two nodes (2 3, back edge 3->2): the snarls inside the loop have no
# interior, so the site is the one around it, 1..4, measured over every copy (60 bp a copy).
{ printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTAC\nS\t2\t%s\nS\t3\t%s\nS\t4\tGGGGCCCCAA\n' \
    "$(printf 'T%.0s' $(seq 30))" "$(printf 'GA%.0s' $(seq 15))"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t3\t+\t2\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'P\tref\t1+,2+,3+,4+\t*\nP\th1\t1+,2+,3+,2+,3+,4+\t*\nP\th2\t1+,2+,3+,2+,3+,2+,3+,4+\t*\n'; } > "$IO/loop2.gfa"
"$BIN" bubble -i "$IO/loop2.gfa" -r ref -o "$IO/loop2" -q >/dev/null 2>&1
check "a two-node loop is one site around it, measured over every copy" \
  "$(tail -n +2 "$IO/loop2.bubbles.csv" | cut -d, -f2,4,13,14 | paste -sd' ' -)" "1,4,60,180"
# A copy whose unit holds a site (2 (3|4) 5, back edge 5->2): the snarl 2..5 inside the loop is widened to
# 1..6, so the expansion is measured whole, not as its 1 bp interior.
{ u=$(printf 'ACGTTGCA%.0s' $(seq 15))
  printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTAC\nS\t2\t%s\nS\t3\tA\nS\t4\tC\nS\t5\t%s\nS\t6\tGGGGCCCCAA\n' "$u" "$u"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t3\t+\t5\t+\t0M\nL\t4\t+\t5\t+\t0M\nL\t5\t+\t2\t+\t0M\nL\t5\t+\t6\t+\t0M\n'
  printf 'P\tref\t1+,2+,3+,5+,6+\t*\nP\th1\t1+,2+,3+,5+,2+,4+,5+,6+\t*\nP\th2\t1+,2+,4+,5+,2+,3+,5+,2+,4+,5+,6+\t*\n'; } > "$IO/loopsite.gfa"
"$BIN" bubble -i "$IO/loopsite.gfa" -r ref -o "$IO/loopsite" -q >/dev/null 2>&1
check "a site inside a loop is widened and its expansion measured whole" \
  "$(tail -n +2 "$IO/loopsite.bubbles.csv" | cut -d, -f2,4,13,14 | paste -sd' ' -)" "1,6,241,723"
# A merge chain that fails a cap merges the longest runs of its parts that pass: three 1 bp sites 10 bp
# apart fused span 23 bp, over --max-variant-bp 15; the first two (12 bp) still merge.
{ printf 'H\tVN:Z:1.0\nS\t1\tAAAAAAAAAA\nS\t2\tC\nS\t3\tG\nS\t4\tTTTTTTTTTT\nS\t5\tA\nS\t6\tT\n'
  printf 'S\t7\tCCCCCCCCCC\nS\t8\tA\nS\t9\tG\nS\t10\tGTGTGTGTGT\n'
  for e in "1 2" "1 3" "2 4" "3 4" "4 5" "4 6" "5 7" "6 7" "7 8" "7 9" "8 10" "9 10"; do
    set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tpA\t1+,2+,4+,5+,7+,8+,10+\t*\nP\tpB\t1+,3+,4+,6+,7+,9+,10+\t*\n'; } > "$IO/chain3.gfa"
"$BIN" bubble -i "$IO/chain3.gfa" -r pA -o "$IO/chain3" --min-variant-bp 0 --max-variant-bp 15 -q >/dev/null 2>&1
check "a failed merge chain keeps the run of its parts that passes" \
  "$(tail -n +2 "$IO/chain3.bubbles.csv" | cut -d, -f2,4 | paste -sd' ' -)" "1,7 7,10"
# With --snarls the graph is not sorted: ids still follow the reference (9 2 1 7 4), not the node ids.
{ printf 'H\tVN:Z:1.0\nS\t9\tACGTACGTAC\nS\t2\tTTTT\nS\t3\tGG\nS\t1\tCCCCCCCCCCCC\nS\t7\tAAAA\nS\t8\tC\nS\t4\tGATTACAGAT\n'
  for e in "9 2" "9 3" "2 1" "3 1" "1 7" "1 8" "7 4" "8 4"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref\t9+,2+,1+,7+,4+\t*\nP\th1\t9+,3+,1+,8+,4+\t*\n'; } > "$IO/order.gfa"
printf '{"start": {"node_id": "9"}, "end": {"node_id": "1"}}\n{"start": {"node_id": "1"}, "end": {"node_id": "4"}}\n' \
  > "$IO/order.jsonl"
"$BIN" bubble -i "$IO/order.gfa" -r ref --snarls "$IO/order.jsonl" -o "$IO/order" --min-variant-bp 0 \
  --merge-nearby-bp 0 -q >/dev/null 2>&1
check "with --snarls, site ids follow the reference" \
  "$(tail -n +2 "$IO/order.bubbles.csv" | cut -d, -f1,2,4 | paste -sd' ' -)" "1,9,1 2,1,4"
# A haplotype that enters a site, turns back out through its source and reaches the sink from outside does
# not cross it: h walks 1+ 2+ 2- 1- 4+ 3-, so the site 1..3 has the reference's allele only.
{ printf 'H\tVN:Z:1.0\nS\t1\tAAAAAAAAAA\nS\t2\t%s\nS\t3\tGGGGGGGGGG\nS\t4\tTTTTT\n' "$(printf 'C%.0s' $(seq 60))"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t2\t+\t2\t-\t0M\nL\t2\t-\t1\t-\t0M\nL\t1\t-\t4\t+\t0M\nL\t4\t+\t3\t-\t0M\n'
  printf 'P\tref\t1+,2+,3+\t*\nP\th\t1+,2+,2-,1-,4+,3-\t*\n'; } > "$IO/hairpin.gfa"
printf '{"start": {"node_id": "1"}, "end": {"node_id": "3"}}\n' > "$IO/hairpin.jsonl"
"$BIN" bubble -i "$IO/hairpin.gfa" -r ref --snarls "$IO/hairpin.jsonl" -o "$IO/hairpin" -q >/dev/null 2>&1
check "a walk reaching a site's sink from outside is no crossing" \
  "$?:$(tail -n +2 "$IO/hairpin.bubbles.csv" 2>/dev/null | cut -d, -f2,4,8,9,14)" "0:1,3,1,1,60"
# A walk with a real crossing and such a turn elsewhere keeps the crossing: h crosses through 5 (70 bp), then
# turns at 2 and reaches 3 from outside. The same walk stored on the other strand gives the same.
hairpin2() {
  local name="$1" walk="$2"
  { printf 'H\tVN:Z:1.0\nS\t1\tAAAAAAAAAA\nS\t2\t%s\nS\t3\tGGGGGGGGGG\nS\t4\tTTTTT\n' "$(printf 'C%.0s' $(seq 60))"
    printf 'S\t5\t%s\nS\t6\tTTTTTCCCCC\n' "$(printf 'A%.0s' $(seq 70))"
    for e in "1 + 2 +" "2 + 3 +" "2 + 2 -" "2 - 1 -" "1 - 4 +" "4 + 3 -" "1 + 5 +" "5 + 3 +" "3 + 6 +" "6 + 1 +" "1 + 3 +"; do
      set -- $e; printf 'L\t%s\t%s\t%s\t%s\t0M\n' "$1" "$2" "$3" "$4"; done
    printf 'P\tref\t1+,2+,3+\t*\nP\th\t%s\t*\n' "$walk"; } > "$IO/$name.gfa"
  "$BIN" bubble -i "$IO/$name.gfa" -r ref --snarls "$IO/hairpin.jsonl" -o "$IO/$name" --min-variant-bp 0 -q >/dev/null 2>&1
  echo "$?:$(tail -n +2 "$IO/$name.bubbles.csv" 2>/dev/null | cut -d, -f8,9,11,14)"
}
check "  ... but a real crossing elsewhere on that walk is kept" \
  "$(hairpin2 hairpin2f 1+,5+,3+,6+,1+,2+,2-,1-,4+,3-)" "0:2,2,1,70"
check "  ... on either strand" "$(hairpin2 hairpin2r 3+,4-,1+,2+,2-,1-,6-,3-,5-,1-)" "0:2,2,1,70"
# Narrowed to its crossings, such a walk reads as it would without the turn: a deletion crossing beside
# one through the interior is left out of the window either way, which spans the interior crossings only.
check "  ... and reads as the same walk without the turn" \
  "$(hairpin2 hairpin2d 1+,3+,6+,1+,2+,3+,6+,1+,2+,2-,1-,4+,3-) $(hairpin2 hairpin2e 1+,3+,6+,1+,2+,3+)" \
  "0:2,1,0,60 0:2,1,0,60"
# --snarls on a W-line graph: the checked-in vg snarls give the two sites with their reference order.
SYN="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/synthetic_data"
"$BIN" bubble -i "$SYN/syn_w.gfa" --snarls "$SYN/syn_w.snarls.jsonl" -r synref --merge-nearby-bp 0 \
  -o "$IO/synw" -q >/dev/null 2>&1
check "--snarls on a W-line graph" \
  "$?:$(tail -n +2 "$IO/synw.bubbles.csv" 2>/dev/null | cut -d, -f1,2,4,8,9,13,14 | paste -sd' ' -)" \
  "0:1,1,2,4,2,80,80 2,2,3,4,2,135,210"
# Fusing by placement links only sites with an insertion or deletion. A 60 bp (AC)n insertion (kept) and a
# 2 bp AC deletion (under the floor) can each be placed onto a SNP between them, two bases of TT away, but not
# into each other. They must not be linked through the SNP into one site; the insertion's site stays as it is.
"${PYTHON:-python3}" - "$IO/bridge.gfa" <<'PYEOF'
import random, sys
r = random.Random(3); seq = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
S = {'f1': seq(300), 'insA': 'AC' * 30, 'r1': 'AC' * 10 + 'TT', 'cR': 'A', 'cA': 'G', 'r2a': 'TT' + 'AC' * 4,
     'dB': 'AC', 'r2b': 'AC' * 5, 'f2': seq(300)}
L = [('f1','insA'),('insA','r1'),('f1','r1'),('r1','cR'),('r1','cA'),('cR','r2a'),('cA','r2a'),('r2a','dB'),
     ('dB','r2b'),('r2a','r2b'),('r2b','f2')]
P = {'ref': 'f1 r1 cR r2a dB r2b f2', 'hA': 'f1 insA r1 cR r2a dB r2b f2', 'hC': 'f1 r1 cA r2a dB r2b f2',
     'hB': 'f1 r1 cR r2a r2b f2', 'hA2': 'f1 insA r1 cR r2a dB r2b f2'}
with open(sys.argv[1], 'w') as o:
    o.write('H\tVN:Z:1.0\n')
    for k, v in S.items(): o.write(f'S\t{k}\t{v}\n')
    for a, b in L: o.write(f'L\t{a}\t+\t{b}\t+\t0M\n')
    for k, v in P.items(): o.write(f'P\t{k}\t' + ','.join(x + '+' for x in v.split()) + '\t*\n')
PYEOF
"$BIN" bubble -i "$IO/bridge.gfa" -r ref -o "$IO/bridge" --merge-nearby-bp 0 -q >/dev/null 2>&1
check "two indels are not fused through a SNP both could be placed onto" \
  "$(tail -n +2 "$IO/bridge.bubbles.csv" | cut -d, -f2,4,13,14 | paste -sd' ' -)" "1,3,0,60"
# An allele longer than the locus's flanks: its cycle with the other allele (5 kb) is longer than the path
# between the locus's ends (2 kb). Rooting the decomposition on that cycle would lose the site; the path is the root.
"${PYTHON:-python3}" - "$IO/longallele.gfa" <<'PYEOF'
import random, sys
r = random.Random(4); seq = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
S = {'1': seq(1000), '2': seq(20), '3': seq(5000), '4': seq(10), '5': seq(20), '6': seq(1000)}
L = [('1','2'),('2','3'),('2','4'),('3','5'),('4','5'),('5','6')]
P = {'ref': '1 2 4 5 6', 'h1': '1 2 3 5 6', 'h2': '1 2 4 5 6', 'h3': '1 2 3 5 6'}
with open(sys.argv[1], 'w') as o:
    o.write('H\tVN:Z:1.0\n')
    for k, v in S.items(): o.write(f'S\t{k}\t{v}\n')
    for a, b in L: o.write(f'L\t{a}\t+\t{b}\t+\t0M\n')
    for k, v in P.items(): o.write(f'P\t{k}\t' + ','.join(x + '+' for x in v.split()) + '\t*\n')
PYEOF
"$BIN" bubble -i "$IO/longallele.gfa" -r ref -o "$IO/longallele" -q >/dev/null 2>&1
check "an allele longer than the locus's flanks keeps its site" \
  "$(tail -n +2 "$IO/longallele.bubbles.csv" | cut -d, -f2,4,13,14 | paste -sd' ' -)" "2,5,10,5000"
# The same, with one haplotype stored on the other strand whose first node, a private 10 bp tip, is walked
# forward: its strand comes from where its ends lie, so it does not mark the backbone as walked both ways
# and the 10 kb cycle of the two alleles still does not root the decomposition.
spec_run rootend >/dev/null <<'EOF'
N S 100 1
N A 5000 2
N B 5000 3
N T 100 4
N e 10 5
N f 10 6
P ref S+,A+,T+
P h1 S+,B+,T+
P h2 S+,A+,T+
P hX f+,T-,A-,S-,e-
EOF
check "a reverse-stored haplotype with a forward end tip does not make a long allele root the decomposition" \
  "$(tail -n +2 "$IO/rootend/b.bubbles.csv" | cut -d, -f2,4,8,14 | paste -sd' ' -)" "2,5,4,5000"
# One haplotype walking the locus twice (an edge from its last node back to its first) crosses every site
# twice. Each crossing is measured alone, so the SNP sites stay SNPs and the 75 bp indel sites stay apart.
locus2x() {   # locus2x <out> <twice>: 150 sites, an 80 vs 5 bp one at 25, 75 and 125, SNPs elsewhere
  "${PYTHON:-python3}" - "$1" "$2" <<'PYEOF'
import random, sys
r = random.Random(1); seq = lambda n: ''.join(r.choice('ACGT') for _ in range(n))
out, twice, N = sys.argv[1], sys.argv[2] == '1', 150
L = ['H\tVN:Z:1.0']
bb = lambda i: 3 * i + 1
for i in range(N):
    L.append(f'S\t{bb(i)}\t{seq(30)}')
    a, b = (seq(80), seq(5)) if i % 50 == 25 else ('A', 'G')
    L += [f'S\t{3*i+2}\t{a}', f'S\t{3*i+3}\t{b}']
    for x in (2, 3): L += [f'L\t{bb(i)}\t+\t{3*i+x}\t+\t0M', f'L\t{3*i+x}\t+\t{bb(i+1)}\t+\t0M']
L.append(f'S\t{bb(N)}\t{seq(30)}')
ref = [s for i in range(N) for s in (f'{bb(i)}+', f'{3*i+2}+')] + [f'{bb(N)}+']
alt = [s for i in range(N) for s in (f'{bb(i)}+', f'{3*i+3}+')] + [f'{bb(N)}+']
L.append('P\tref\t' + ','.join(ref) + '\t*')
for k in range(6): L.append(f'P\th{k}\t' + ','.join(ref if k % 2 else alt) + '\t*')
if twice: L += [f'L\t{bb(N)}\t+\t{bb(0)}\t+\t0M', 'P\tdup\t' + ','.join(ref + alt) + '\t*']
open(out, 'w').write('\n'.join(L) + '\n')
PYEOF
}
locus2x "$IO/once.gfa" 0; locus2x "$IO/twice.gfa" 1
"$BIN" bubble -i "$IO/once.gfa" -r ref -o "$IO/once" -q >/dev/null 2>&1
"$BIN" bubble -i "$IO/twice.gfa" -r ref -o "$IO/twice" -q >/dev/null 2>&1
check "a haplotype walking the locus twice does not merge its sites" \
  "$(tail -n +2 "$IO/twice.bubbles.csv" | cut -d, -f2,4 | paste -sd' ' -)" \
  "$(tail -n +2 "$IO/once.bubbles.csv" | cut -d, -f2,4 | paste -sd' ' -)"
check "  ... which are the three indel sites" "$(rows "$IO/twice.bubbles.csv")" "3"

echo
if [ "$fails" -eq 0 ]; then echo "bubble_stats: all assertions passed"; exit 0; fi
echo "bubble_stats: $fails assertion(s) FAILED"; exit 1
