#!/usr/bin/env bash
# inspect_stats.sh - contract assertions for `panvar inspect`.
#
# real_smoke.sh only checks that inspect's files exist and that the gzip is readable, which passes
# whether or not the contents are right; these check the contents.
#
# Every fixture here is small enough that the expected output is derivable by hand.
#
#   inspect_stats.sh <panvar-binary> <out-dir>
set -uo pipefail

BIN="${1:?usage: inspect_stats.sh <panvar> <outdir>}"
OUT="${2:?}"
# Absolute, so a test that runs from another folder still finds both.
case "$BIN" in /*) ;; */*) BIN="$PWD/$BIN" ;; esac
mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)/run.$$.$(date +%s)"
rm -rf "$OUT"; mkdir -p "$OUT"
fails=0
ok()  { printf "  ok   %s\n" "$1"; }
bad() { printf "  FAIL %s\n" "$1"; fails=$((fails + 1)); }
heads() { gzip -dc "$1" 2>/dev/null | grep -c '^>' | tr -d ' '; }
hdr()   { gzip -dc "$1" 2>/dev/null | grep '^>' | sed 's/^>//'; }

# ---------------------------------------------------------------- a pure deletion is an allele
# `full` carries interior node 2; `del` goes straight from source to sink. Both are real alleles of the
# same site and `bubble` reports distinct_alleles=2: a walk with no interior node is a crossing too.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGG\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tfull\t1+,2+,3+\t*\nP\tdel\t1+,3+\t*\n'; } > "$OUT/del.gfa"
"$BIN" bubble -i "$OUT/del.gfa" -r full -o "$OUT/db" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" inspect -i "$OUT/db.sorted.gfa" -c "$OUT/db.bubbles.csv" --bubble-id 1 -o "$OUT/d" -q >/dev/null 2>&1

[ "$(heads "$OUT/d/bubble_1/d.bubble_1.paths.fa.gz")" = "2" ] \
  && ok "both alleles are emitted, including the one with no interior node" \
  || bad "expected 2 sequences, got $(heads "$OUT/d/bubble_1/d.bubble_1.paths.fa.gz")"
hdr "$OUT/d/bubble_1/d.bubble_1.paths.fa.gz" | grep -q '^del ' \
  && ok "the pure-deletion haplotype appears by name" \
  || bad "the deletion haplotype is missing from the FASTA"
# 1 + 3 = 20 bp with node 2 (10 bp) deleted; the full walk is 30.
hdr "$OUT/d/bubble_1/d.bubble_1.paths.fa.gz" | grep -q '^del .*length_bp=20' \
  && ok "the deletion allele spells 20 bp (source + sink, interior dropped)" \
  || bad "deletion length wrong: $(hdr "$OUT/d/bubble_1/d.bubble_1.paths.fa.gz" | grep '^del ')"
# The fallback still has to report interval metadata, not leave it at whatever was last set.
hdr "$OUT/d/bubble_1/d.bubble_1.paths.fa.gz" | grep -q '^del .*interval=0-1' \
  && ok "the empty-interior fallback reports its own interval" \
  || bad "deletion interval metadata wrong: $(hdr "$OUT/d/bubble_1/d.bubble_1.paths.fa.gz" | grep '^del ')"
# The count matrix must gain a row for it, reading zero for the deleted node.
[ "$(awk 'NR>1' "$OUT/d/bubble_1/d.bubble_1.node_counts.tsv" | wc -l | tr -d ' ')" = "2" ] \
  && ok "the count matrix has a row per allele" \
  || bad "count matrix has $(awk 'NR>1' "$OUT/d/bubble_1/d.bubble_1.node_counts.tsv" | wc -l | tr -d ' ') rows, expected 2"
# Columns: path_name, path_length_bp, crossing, crossings, then one per interior node.
awk '$1=="del"{print $5}' "$OUT/d/bubble_1/d.bubble_1.node_counts.tsv" | grep -q '^0:0:0$' \
  && ok "the deleted node reads 0:0:0 for the carrier of the deletion" \
  || bad "deleted node count is $(awk '$1=="del"{print $5}' "$OUT/d/bubble_1/d.bubble_1.node_counts.tsv"), expected 0:0:0"

# ---------------------------------------------------------------- reverse traversal
# `rev` crosses the same site on the opposite strand. The walk is canonicalized source->sink, so it must
# spell the same 30 bp as the forward haplotype and be reported as sink->source.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGG\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tfull\t1+,2+,3+\t*\nP\tdel\t1+,3+\t*\nP\trev\t3-,2-,1-\t*\n'; } > "$OUT/rev.gfa"
"$BIN" bubble -i "$OUT/rev.gfa" -r full -o "$OUT/rb" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" inspect -i "$OUT/rb.sorted.gfa" -c "$OUT/rb.bubbles.csv" --bubble-id 1 -o "$OUT/r" -q >/dev/null 2>&1
[ "$(heads "$OUT/r/bubble_1/r.bubble_1.paths.fa.gz")" = "3" ] \
  && ok "a reverse traversal is emitted alongside the forward ones" \
  || bad "expected 3 sequences, got $(heads "$OUT/r/bubble_1/r.bubble_1.paths.fa.gz")"
hdr "$OUT/r/bubble_1/r.bubble_1.paths.fa.gz" | grep -q '^rev .*source_to_sink=0' \
  && ok "a reverse traversal is reported as sink->source" \
  || bad "reverse traversal metadata wrong: $(hdr "$OUT/r/bubble_1/r.bubble_1.paths.fa.gz" | grep '^rev ')"
hdr "$OUT/r/bubble_1/r.bubble_1.paths.fa.gz" | grep -q '^rev .*length_bp=30' \
  && ok "a reverse traversal spells the same length as the forward walk" \
  || bad "reverse traversal length wrong: $(hdr "$OUT/r/bubble_1/r.bubble_1.paths.fa.gz" | grep '^rev ')"

# ---------------------------------------------------------------- a repeated node inside the site
# A boundary cannot repeat within a discovered snarl -- a boundary carrying a cycle stops being a cut
# vertex, and the decomposition then emits no snarl there at all. The reachable case is a repeated INTERIOR node,
# which is what a tandem array is: `rep` walks 1,2,3,2,4 and visits node 2 twice.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGG\nS\t4\tTTTTTTTTTT\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t1\t+\t4\t+\t0M\n'
  printf 'L\t2\t+\t3\t+\t0M\nL\t3\t+\t2\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'P\tfull\t1+,2+,4+\t*\nP\tdel\t1+,4+\t*\nP\trep\t1+,2+,3+,2+,4+\t*\n'; } > "$OUT/rep.gfa"
"$BIN" bubble -i "$OUT/rep.gfa" -r full -o "$OUT/pb" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" inspect -i "$OUT/pb.sorted.gfa" -c "$OUT/pb.bubbles.csv" -o "$OUT/p" -q >/dev/null 2>&1
[ "$(heads "$OUT/p/bubble_1/p.bubble_1.paths.fa.gz")" = "3" ] \
  && ok "all three alleles are emitted when one revisits an interior node" \
  || bad "expected 3 sequences at the repeated-node site, got $(heads "$OUT/p/bubble_1/p.bubble_1.paths.fa.gz")"
[ "$(hdr "$OUT/p/bubble_1/p.bubble_1.paths.fa.gz" | grep -c 'length_bp=0')" = "0" ] \
  && ok "no allele spells an empty sequence at a repeated interior node" \
  || bad "a zero-length walk was emitted at a repeated interior node"
# The revisited node must be COUNTED twice, not deduplicated: multiplicity is the copy number.
awk '$1=="rep"{print $5}' "$OUT/p/bubble_1/p.bubble_1.node_counts.tsv" | grep -q '^2:2:0$' \
  && ok "a revisited interior node is counted twice (2:2:0), not collapsed" \
  || bad "revisited node count is $(awk '$1=="rep"{print $5}' "$OUT/p/bubble_1/p.bubble_1.node_counts.tsv"), expected 2:2:0"

# ---------------------------------------------------------------- the graph contract
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGG\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tfull\t1+,2+,3+\t*\nP\tfull\t1+,2+,3+\t*\n'; } > "$OUT/dup.gfa"
"$BIN" inspect -i "$OUT/dup.gfa" -c "$OUT/db.bubbles.csv" --bubble-id 1 -o "$OUT/dupo" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a duplicate path name is refused" \
               || bad "a duplicate path name was accepted (two rows would carry the same label)"

# A non-zero overlap means adjacent nodes share bases, so spelling by concatenation overcounts every
# length inspect reports.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGG\n'
  printf 'L\t1\t+\t2\t+\t1M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tfull\t1+,2+,3+\t*\nP\tdel\t1+,3+\t*\n'; } > "$OUT/ovl.gfa"
"$BIN" inspect -i "$OUT/ovl.gfa" -c "$OUT/db.bubbles.csv" --bubble-id 1 -o "$OUT/ovlo" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a non-zero link overlap is refused" \
               || bad "a 1M overlap was accepted; every path_length_bp would be an overcount"

# ---------------------------------------------------------------- mismatched inputs
# Nothing ties the CSV to the GFA, so a CSV from another graph must be refused, not run to a header-only
# matrix with node lengths of 0.
printf 'bubble_id,source,sink,path_support,min_inside_bp,max_inside_bp,inside_nodes\n99,900,902,2,10,10,"901"\n' \
  > "$OUT/ghost.csv"
"$BIN" inspect -i "$OUT/del.gfa" -c "$OUT/ghost.csv" --bubble-id 99 -o "$OUT/gh" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a bubbles CSV naming nodes absent from the GFA is refused" \
               || bad "a mismatched CSV ran to completion"
# Outputs go into the folder <prefix>/, and site files can sit beside it as <prefix>.*: neither may exist.
left=$(find "$OUT" -maxdepth 1 \( -name gh -o -name 'gh.*' \) | wc -l | tr -d ' ')
[ "$left" = "0" ] \
  && ok "a refused run leaves no partial per-bubble outputs" \
  || bad "mismatched run left $left output(s) at the prefix"

# ---------------------------------------------------------------- output aliasing
"$BIN" inspect -i "$OUT/db.sorted.gfa" -c "$OUT/db.bubbles.csv" --bubble-id 1 \
       --fasta-out "$OUT/same.out" --table-out "$OUT/same.out" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "two outputs naming the same file are refused" \
               || bad "two outputs were allowed to name the same file"
before=$(wc -c < "$OUT/db.sorted.gfa" | tr -d ' ')
"$BIN" inspect -i "$OUT/db.sorted.gfa" -c "$OUT/db.bubbles.csv" --bubble-id 1 \
       --fasta-out "$OUT/db.sorted.gfa" -q >/dev/null 2>&1
[ "$(wc -c < "$OUT/db.sorted.gfa" | tr -d ' ')" = "$before" ] \
  && ok "an output aliasing the input GFA leaves the input intact" \
  || bad "the input GFA was overwritten by an output"

# ---------------------------------------------------------------- long-walk clustering
# Walks thousands of steps long are compared exactly. A cyclic array: node 2 (10 bp) carries a self-loop,
# walked 3000 times by `long` and 500 by `short`, so the shorter walk's visits are all in the longer one.
# The length-weighted Dice similarity is 2 x 5000 / (30000 + 5000) = 0.286; the assertions below bracket
# it in [0.22, 0.30].
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tACGTACGTAC\nS\t2\tTTTTGGGGCC\nS\t3\tGGGGGGGGGG\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\n'
  units=""; for i in $(seq 1 3000); do units="$units,2+"; done
  printf 'P\tlong\t1+%s,3+\t*\n' "$units"
  units=""; for i in $(seq 1 500); do units="$units,2+"; done
  printf 'P\tshort\t1+%s,3+\t*\n' "$units"; } > "$OUT/arr.gfa"
"$BIN" bubble -i "$OUT/arr.gfa" -r long -o "$OUT/ab" --min-variant-bp 0 -q >/dev/null 2>&1
nclust() { awk 'NR>1{print $1}' "$1" 2>/dev/null | sort -u | wc -l | tr -d ' '; }

"$BIN" inspect -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" --cluster --cluster-similarity 0.22 \
       -o "$OUT/a22" -q >/dev/null 2>&1
[ "$(nclust "$OUT/a22/bubble_1/a22.bubble_1.clusters.tsv")" = "1" ] \
  && ok "long walks at a 6x length ratio (similarity 0.286) cluster at 0.22" \
  || bad "expected 1 cluster at 0.22, got $(nclust "$OUT/a22/bubble_1/a22.bubble_1.clusters.tsv")"

"$BIN" inspect -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" --cluster --cluster-similarity 0.30 \
       -o "$OUT/a30" -q >/dev/null 2>&1
[ "$(nclust "$OUT/a30/bubble_1/a30.bubble_1.clusters.tsv")" = "2" ] \
  && ok "and separate at 0.30, bracketing the similarity 0.286" \
  || bad "expected 2 clusters at 0.30, got $(nclust "$OUT/a30/bubble_1/a30.bubble_1.clusters.tsv")"

# ------------------------------------------- copies are compared by their bases
# Two and three copies of a 30 bp unit between the same boundaries: the second walk holds every visit of
# the first and one more copy, so they share 60 of 60 and 90 bases, a similarity of 2 x 60 / 150 = 0.80.
# 0.79 and 0.85 bracket it: the pair clusters at the first and separates at the second.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\tAAAAAAAAAAAAAAAAAAAA\nS\tA\tCCCCCCCCCC\nS\tB\tGGGGGGGGGG\nS\tC\tTTTTTTTTTT\n'
  printf 'S\tG\tAAAAAAAAAAAAAAAAAAAA\n'
  printf 'L\tF\t+\tA\t+\t0M\nL\tA\t+\tB\t+\t0M\nL\tB\t+\tC\t+\t0M\nL\tC\t+\tA\t+\t0M\nL\tC\t+\tG\t+\t0M\n'
  printf 'P\th1\tF+,A+,B+,C+,A+,B+,C+,G+\t*\n'
  printf 'P\th3\tF+,A+,B+,C+,A+,B+,C+,A+,B+,C+,G+\t*\n'; } > "$OUT/exact.gfa"
"$BIN" bubble -i "$OUT/exact.gfa" -r h1 -o "$OUT/ex" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" inspect -i "$OUT/ex.sorted.gfa" -c "$OUT/ex.bubbles.csv" --cluster --cluster-similarity 0.79 \
       -o "$OUT/ex79" -q >/dev/null 2>&1
[ "$(nclust "$OUT/ex79/bubble_1/ex79.bubble_1.clusters.tsv")" = "1" ] \
  && ok "two and three copies (similarity 0.80) cluster at 0.79" \
  || bad "expected 1 cluster at 0.79, got $(nclust "$OUT/ex79/bubble_1/ex79.bubble_1.clusters.tsv")"
"$BIN" inspect -i "$OUT/ex.sorted.gfa" -c "$OUT/ex.bubbles.csv" --cluster --cluster-similarity 0.85 \
       -o "$OUT/ex85" -q >/dev/null 2>&1
[ "$(nclust "$OUT/ex85/bubble_1/ex85.bubble_1.clusters.tsv")" = "2" ] \
  && ok "and separate at 0.85" \
  || bad "expected 2 clusters at 0.85, got $(nclust "$OUT/ex85/bubble_1/ex85.bubble_1.clusters.tsv")"

# Similarity is weighted by node length: two alleles 1 kb + SNP + 1 kb apart share 2000 of 2001 bases
# (0.9995) and cluster at 0.90, while one more 1 kb copy (0.80) separates. Counted per step, the SNP pair
# would share 2 of 3 elements (0.67) and separate too.
kb=$(printf 'ACGTTGCA%.0s' $(seq 1 125))
{ printf 'H\tVN:Z:1.0\nS\tS\tGATTACAGATTACAGATTAC\nS\tX\t%s\nS\tA\tA\nS\tC\tC\nS\tY\t%s\nS\tT\tTTGGCCAATTGGCCAATTGG\n' "$kb" "$kb"
  for e in "S X" "X A" "X C" "A Y" "C Y" "Y Y" "Y T"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tsnpA\tS+,X+,A+,Y+,T+\t*\nP\tsnpC\tS+,X+,C+,Y+,T+\t*\nP\tdup\tS+,X+,A+,Y+,Y+,T+\t*\n'; } > "$OUT/lw.gfa"
printf 'bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,S,+,T,+,3,2001,3001,"X;A;C;Y"\n' \
  > "$OUT/lw.csv"
"$BIN" inspect -i "$OUT/lw.gfa" -c "$OUT/lw.csv" --cluster --cluster-similarity 0.90 -o "$OUT/lw" -q >/dev/null 2>&1
lw=$(tail -n +2 "$OUT/lw/bubble_1/lw.bubble_1.clusters.tsv" 2>/dev/null | cut -f4 | sort | paste -sd' ' -)
[ "$lw" = "dup snpA;snpC" ] \
  && ok "walks are compared by bases: a SNP between 1 kb nodes clusters, an extra 1 kb copy does not" \
  || bad "length-weighted clusters: '$lw', expected 'dup snpA;snpC'"

# ------------------------------------------- the representative does not follow GFA record order
# Identical walks are pooled; taking the first member met as `representative_path` would make the
# representative, and through it the sort key of the whole TSV, depend on the order of the P records. The
# representative is the lexicographically smallest member of the most carried walk that founded the group.
{ printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTAC\nS\t2\tTTTT\nS\t3\tGGGGCCCCAA\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tzzz\t1+,2+,3+\t*\nP\taaa\t1+,2+,3+\t*\nP\tref\t1+,3+\t*\n'; } > "$OUT/repord.gfa"
{ printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTAC\nS\t2\tTTTT\nS\t3\tGGGGCCCCAA\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\taaa\t1+,2+,3+\t*\nP\tzzz\t1+,2+,3+\t*\nP\tref\t1+,3+\t*\n'; } > "$OUT/repord2.gfa"
for g in repord repord2; do
  "$BIN" bubble -i "$OUT/$g.gfa" -r ref -o "$OUT/$g" --min-variant-bp 0 -q >/dev/null 2>&1
  "$BIN" inspect -i "$OUT/$g.sorted.gfa" -c "$OUT/$g.bubbles.csv" --cluster --cluster-similarity 0.90 \
         -o "$OUT/${g}_i" -q >/dev/null 2>&1
done
rep1=$(awk -F'\t' 'NR>1 && $2==2{print $3}' "$OUT/repord_i/bubble_1/repord_i.bubble_1.clusters.tsv")
rep2=$(awk -F'\t' 'NR>1 && $2==2{print $3}' "$OUT/repord2_i/bubble_1/repord2_i.bubble_1.clusters.tsv")
[ -n "$rep1" ] && [ "$rep1" = "$rep2" ] && [ "$rep1" = "aaa" ] \
  && ok "the representative is order-independent and lexicographically smallest (aaa)" \
  || bad "representative depends on P-record order: '$rep1' vs '$rep2' (expected aaa both ways)"

# ------------------------------------------- two output-name rules, pinned directly

# 1. A DERIVED output collides with an explicit one. The preflight has to run AFTER bubble selection,
#    because the derived names contain the bubble ids; an explicit --table-out naming
#    <prefix>.bubble_N.node_lengths.tsv would otherwise clobber one of them.
"$BIN" inspect -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" --bubble-id 1 \
       --table-out "$OUT/coll/bubble_1/coll.bubble_1.node_lengths.tsv" -o "$OUT/coll" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "an explicit output colliding with a DERIVED one is refused" \
               || bad "an explicit --table-out naming a derived output was accepted"

# 2. The bubbles CSV names nodes that DO exist in the graph, but no path crosses the pair. A header-only
#    count table would read as "this site has no variation" rather than "this CSV does not describe
#    this graph".
hdr=$(head -1 "$OUT/ab.bubbles.csv")
{ echo "$hdr"
  awk -F, -v h="$hdr" 'NR==2{n=split(h,c,","); for(i=1;i<=n;i++){ if(c[i]=="source") s=i; if(c[i]=="sink") k=i }
        $s=$k; print }' OFS=, "$OUT/ab.bubbles.csv"; } > "$OUT/nocross.csv"
"$BIN" inspect -i "$OUT/ab.sorted.gfa" -c "$OUT/nocross.csv" -o "$OUT/nocross" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a selected bubble no path crosses is refused, not reported as empty" \
               || bad "a bubble with present nodes and no crossing path exited 0"

# ------------------------------------------- a stale panel is reported, not silently rescored
# Requiring at least one crossing catches a CSV from a DIFFERENT graph. It does not catch a stale one:
# same node ids, some paths missing, so crossings still happen -- just fewer. The CSV records how many
# `bubble` saw, so the two are compared and a mismatch is named.
{ printf 'H\tVN:Z:1.0\nS\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGG\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tfull\t1+,2+,3+\t*\nP\tdel\t1+,3+\t*\nP\tthird\t1+,2+,3+\t*\n'; } > "$OUT/three.gfa"
"$BIN" bubble -i "$OUT/three.gfa" -r full -o "$OUT/b3" --min-variant-bp 0 -q >/dev/null 2>&1
{ printf 'H\tVN:Z:1.0\nS\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGG\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tfull\t1+,2+,3+\t*\nP\tdel\t1+,3+\t*\n'; } > "$OUT/two.gfa"
warn=$("$BIN" inspect -i "$OUT/two.gfa" -c "$OUT/b3.bubbles.csv" -o "$OUT/stale" 2>&1 | grep -c "differ from the table's path_support")
[ "$warn" -ge 1 ] && ok "a stale panel (2 crossings against path_support=3) is reported" \
                  || bad "no warning when the graph carries fewer paths than the CSV records"
# ...and a matching panel must stay quiet, or the warning is noise.
quiet=$("$BIN" inspect -i "$OUT/b3.sorted.gfa" -c "$OUT/b3.bubbles.csv" -o "$OUT/fresh" 2>&1 | grep -c "differ from the table's path_support")
[ "$quiet" = "0" ] && ok "and a matching panel produces no such warning" \
                   || bad "the stale-panel warning fires on a matching panel"

# ------------------------------------------- the shared boundaries do not make alleles alike
# Two SNP alleles and a pure deletion share nothing but the 200 bp boundaries. Counted, the boundaries
# would make every pair ~0.99 alike; between the boundaries they share nothing: three clusters.
long_a=$(printf 'A%.0s' $(seq 1 200)); long_g=$(printf 'G%.0s' $(seq 1 200))
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\nS\t2\tC\nS\t4\tT\nS\t3\t%s\n' "$long_a" "$long_g"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t4\t+\t0M\nL\t4\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tsnpC\t1+,2+,3+\t*\nP\tsnpT\t1+,4+,3+\t*\nP\tdel\t1+,3+\t*\n'; } > "$OUT/hub.gfa"
printf 'bubble_id,source,sink,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,1,3,3,0,1,"2;4"\n' \
  > "$OUT/hub.csv"
"$BIN" inspect -i "$OUT/hub.gfa" -c "$OUT/hub.csv" --cluster -o "$OUT/hub" -q >/dev/null 2>&1
[ "$(nclust "$OUT/hub/bubble_1/hub.bubble_1.clusters.tsv")" = "3" ] \
  && ok "a pure deletion does not join two unrelated alleles into one cluster" \
  || bad "expected 3 clusters with a deletion beside two SNP alleles, got $(nclust "$OUT/hub/bubble_1/hub.bubble_1.clusters.tsv")"

# ------------------------------------------- every crossing of a path is its own row
# `twice` crosses bubble 1..3 twice, through spacer node 4 between; `delthenins` deletes node 2 on the
# first crossing and keeps it on the second. One window from the first source to the last sink would spell
# the spacer as part of the allele and lose the deletion crossing, so each crossing is a row labelled <path>@k.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tCCCCC\nS\t3\tGGGGGGGGGG\nS\t4\tTTTTTTTTTTTTTTTTTTTT\nS\t5\tACACACACAC\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\nL\t4\t+\t1\t+\t0M\nL\t3\t+\t5\t+\t0M\n'
  printf 'P\tref\t1+,2+,3+,5+\t*\nP\ttwice\t1+,2+,3+,4+,1+,2+,3+,5+\t*\n'
  printf 'P\tdelthenins\t1+,3+,4+,1+,2+,3+,5+\t*\nP\tdel\t1+,3+,5+\t*\n'; } > "$OUT/mc.gfa"
printf 'bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,1,+,3,+,4,0,5,"2"\n' \
  > "$OUT/mc.csv"
"$BIN" inspect -i "$OUT/mc.gfa" -c "$OUT/mc.csv" --cluster -o "$OUT/mc" -q >/dev/null 2>"$OUT/mc.err"
rows=$(awk -F'\t' 'NR>1{print $1":"$2":"$3"/"$4":"$5}' "$OUT/mc/bubble_1/mc.bubble_1.node_counts.tsv" | tr '\n' ' ')
[ "$rows" = "ref:25:1/1:1:1:0 twice@1:25:1/2:1:1:0 twice@2:25:2/2:1:1:0 delthenins@1:20:1/2:0:0:0 delthenins@2:25:2/2:1:1:0 del:20:1/1:0:0:0 " ] \
  && ok "each crossing is a row: name@k, its own length and counts, crossing k of n" \
  || bad "per-crossing rows wrong: $rows"
# The longest allele is 25 bp (1+2+3); a window over both crossings of `twice` spells 70.
[ "$(hdr "$OUT/mc/bubble_1/mc.bubble_1.paths.fa.gz" | sed 's/.*length_bp=\([0-9]*\).*/\1/' | sort -n | tail -1)" = "25" ] \
  && ok "no allele spells the spacer between two crossings" \
  || bad "an allele spelled across the spacer between two crossings"
hdr "$OUT/mc/bubble_1/mc.bubble_1.paths.fa.gz" | grep -q '^delthenins@1 .*crossing=1/2' \
  && ok "the FASTA header gives the crossing as k/n" \
  || bad "FASTA crossing metadata wrong: $(hdr "$OUT/mc/bubble_1/mc.bubble_1.paths.fa.gz" | grep '^delthenins@1 ')"
awk -F'\t' 'NR>1 && $3=="del"{print $4}' "$OUT/mc/bubble_1/mc.bubble_1.clusters.tsv" | grep -qx 'del;delthenins@1' \
  && ok "a deletion crossing clusters with the deletion, the path's other crossing elsewhere" \
  || bad "cluster of the deletion crossing wrong: $(cat "$OUT/mc/bubble_1/mc.bubble_1.clusters.tsv")"
# path_support counts paths, not rows: 4 paths cross, 6 rows, no stale-panel warning.
! grep -q "differ from the table's path_support" "$OUT/mc.err" \
  && ok "the panel check counts paths, not crossings" \
  || bad "a path crossing twice tripped the stale-panel warning"

# A crossing walked backwards, sink- ... source-, is the same bubble crossed on the other strand:
# `inv` crosses forward, then again inverted through an inversion edge.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tCCCCC\nS\t3\tGGGGGGGGGG\nS\t4\tTTTTTTTTTTTTTTTTTTTT\nS\t5\tACACACACAC\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\nL\t4\t+\t3\t-\t0M\nL\t1\t-\t5\t+\t0M\n'
  printf 'P\tref\t1+,2+,3+,4+\t*\nP\tinv\t1+,2+,3+,4+,3-,1-,5+\t*\n'; } > "$OUT/inv.gfa"
printf 'bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,1,+,3,+,2,0,5,"2"\n' \
  > "$OUT/inv.csv"
"$BIN" inspect -i "$OUT/inv.gfa" -c "$OUT/inv.csv" -o "$OUT/inv" -q >/dev/null 2>&1
hdr "$OUT/inv/bubble_1/inv.bubble_1.paths.fa.gz" | grep -q '^inv@2 .*length_bp=20 source_to_sink=0 .*crossing=2/2' \
  && ok "an inverted second crossing is its own row, read source->sink" \
  || bad "inverted crossing wrong: $(hdr "$OUT/inv/bubble_1/inv.bubble_1.paths.fa.gz" | tr '\n' '|')"

# ------------------------------------------- argument refusals and a stale clusters table
"$BIN" inspect -i "$OUT/mc.gfa" -c "$OUT/mc.csv" --bubble-id 0 -o "$OUT/b0" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "--bubble-id 0 is refused" \
               || bad "--bubble-id 0 was accepted"
"$BIN" inspect -i "$OUT/mc.gfa" -c "$OUT/mc.csv" --cluster-similarity 0.5 -o "$OUT/cs" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "--cluster-similarity without --cluster is refused" \
               || bad "--cluster-similarity without --cluster was silently ignored"
"$BIN" inspect -i "$OUT/mc.gfa" -c "$OUT/mc.csv" -o "$OUT/mc" -q >/dev/null 2>&1
[ ! -e "$OUT/mc/bubble_1/mc.bubble_1.clusters.tsv" ] \
  && ok "a rerun without --cluster removes the earlier clusters table" \
  || bad "a stale clusters.tsv survived a rerun without --cluster"


# ---------------------------------------------------------------- a path touching a site from outside
# `outer` walks 1- 2+ 3- against a 1+..3+ site and `junc` passes 3+ 1+ outside both boundaries; neither
# crosses the site (`ref` does, so the CSV's orientations hold). Neither gets a row: the first would read
# as an allele, the second as a reverse-complemented "deletion".
{ printf 'H\tVN:Z:1.0\nS\t1\tAAAAAAAAAC\nS\t2\tCCCCC\nS\t3\tGGGGGGGGGT\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t-\t2\t+\t0M\nL\t2\t+\t3\t-\t0M\nL\t3\t+\t1\t+\t0M\n'
  printf 'P\tref\t1+,2+,3+\t*\nP\touter\t1-,2+,3-\t*\nP\tjunc\t3+,1+\t*\n'; } > "$OUT/touch.gfa"
printf 'bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,1,+,3,+,1,5,5,"2"\n' > "$OUT/touch.csv"
"$BIN" inspect -i "$OUT/touch.gfa" -c "$OUT/touch.csv" --bubble-id 1 -o "$OUT/touch" -q >/dev/null 2>&1
[ "$(hdr "$OUT/touch/bubble_1/touch.bubble_1.paths.fa.gz" | cut -d' ' -f1 | paste -sd' ' -)" = "ref" ] \
  && ok "paths touching the boundaries from outside are not crossings" \
  || bad "touching paths: rows '$(hdr "$OUT/touch/bubble_1/touch.bubble_1.paths.fa.gz" | cut -d' ' -f1 | paste -sd' ' -)', expected 'ref'"
"$BIN" inspect -i "$OUT/touch.gfa" -c "$OUT/touch.csv" --bubble-id 1 -o "$OUT/touchn" >/dev/null 2>"$OUT/touchn.err"
grep -q "without crossing it (no row): 2" "$OUT/touchn.err" \
  && ok "paths touching both boundaries without a row are counted in a note" \
  || bad "no note for the paths touching the site from outside"
# A CSV without orientation columns reads every boundary '+'. `inv` crosses 1+..3+ while ref, del and ref2
# cross 1+..3-: read strictly, that one path would drop the other three, so without the columns the
# orientations are not used.
{ printf 'H\tVN:Z:1.0\nS\t1\tAAAAAAAAAC\nS\t2\tCCCCC\nS\t3\tGGGGGGGGGT\nS\t4\tTTTTTTTTTA\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t-\t0M\nL\t3\t-\t4\t+\t0M\nL\t1\t+\t3\t-\t0M\nL\t2\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'P\tref\t1+,2+,3-,4+\t*\nP\tdel\t1+,3-,4+\t*\nP\tinv\t1+,2+,3+,4+\t*\nP\tref2\t1+,2+,3-,4+\t*\n'; } > "$OUT/legacy.gfa"
printf 'bubble_id,source,sink,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,1,3,0,0,0,2\n' > "$OUT/legacy.csv"
"$BIN" inspect -i "$OUT/legacy.gfa" -c "$OUT/legacy.csv" --bubble-id 1 -o "$OUT/legacy" -q >/dev/null 2>&1
[ "$(hdr "$OUT/legacy/bubble_1/legacy.bubble_1.paths.fa.gz" | cut -d' ' -f1 | paste -sd' ' -)" = "ref del inv ref2" ] \
  && ok "a CSV without orientation columns keeps every path's allele" \
  || bad "orientation-less CSV: rows '$(hdr "$OUT/legacy/bubble_1/legacy.bubble_1.paths.fa.gz" | cut -d' ' -f1 | paste -sd' ' -)'"
# The sweep takes ids as inspect writes them: bubble_007 is not one.
mkdir -p "$OUT/touch/bubble_007"; printf 'mine\n' > "$OUT/touch/bubble_007/touch.bubble_007.node_counts.tsv"
"$BIN" inspect -i "$OUT/touch.gfa" -c "$OUT/touch.csv" -o "$OUT/touch" -q >/dev/null 2>&1
[ -e "$OUT/touch/bubble_007/touch.bubble_007.node_counts.tsv" ] && ok "a zero-padded bubble file is not swept" \
  || bad "the sweep removed a zero-padded name inspect never writes"

# ---------------------------------------------------------------- a copy-number ladder does not chain
# 1..32 copies of one looped unit. Walks one copy apart are ~97% alike, so connected components would chain
# them into one cluster at 0.90; each member must be within the threshold of its representative.
{ printf 'H\tVN:Z:1.0\nS\t1\tGATTACA\nS\t2\tACGTTGCAAGG\nS\t3\tTTGGCC\nL\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\n'
  for n in $(seq 1 32); do w="1+"; for ((k = 0; k < n; k++)); do w="$w,2+"; done; printf 'P\tcn%02d\t%s,3+\t*\n' "$n" "$w"; done; } > "$OUT/kiv.gfa"
printf 'bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,1,+,3,+,32,11,352,"2"\n' > "$OUT/kiv.csv"
"$BIN" inspect -i "$OUT/kiv.gfa" -c "$OUT/kiv.csv" --cluster -o "$OUT/kiv" -q >/dev/null 2>&1
spread=$(tail -n +2 "$OUT/kiv/bubble_1/kiv.bubble_1.clusters.tsv" 2>/dev/null | awk -F'\t' '{n=split($4,a,";"); d=substr(a[n],3)-substr(a[1],3); if (d > m) m=d} END{print (NR>0 && m<=8) ? "bounded" : "chained " m}')
[ "$spread" = "bounded" ] \
  && ok "a copy-number ladder splits: no cluster spans more than 8 copies at 0.90" \
  || bad "copy-number ladder clusters: $spread ($(tail -n +2 "$OUT/kiv/bubble_1/kiv.bubble_1.clusters.tsv" 2>/dev/null | wc -l | tr -d ' ') clusters)"

# ---------------------------------------------------------------- labels and outputs
# A path named "x@1" crossing once beside a path "x" crossing twice: the labels must not collide.
{ printf 'H\tVN:Z:1.0\nS\t1\tAAAAAAAAAC\nS\t2\tCCCCC\nS\t3\tGGGGGGGGGT\nS\t4\tTTTTT\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\nL\t4\t+\t1\t+\t0M\n'
  printf 'P\tx\t1+,2+,3+,4+,1+,2+,3+\t*\nP\tx@1\t1+,2+,3+\t*\n'; } > "$OUT/lab.gfa"
printf 'bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,1,+,3,+,2,5,5,"2"\n' > "$OUT/lab.csv"
"$BIN" inspect -i "$OUT/lab.gfa" -c "$OUT/lab.csv" --bubble-id 1 -o "$OUT/lab" -q >/dev/null 2>&1
[ "$(tail -n +2 "$OUT/lab/bubble_1/lab.bubble_1.node_counts.tsv" 2>/dev/null | cut -f1 | paste -sd' ' -)" = "x@1 x@2 x@1@1" ] \
  && ok "a path named like a crossing label keeps its own suffix (x@1 x@2 x@1@1)" \
  || bad "crossing labels: '$(tail -n +2 "$OUT/lab/bubble_1/lab.bubble_1.node_counts.tsv" 2>/dev/null | cut -f1 | paste -sd' ' -)'"
# An explicit output named like a derived clusters table is this run's output, not a stale one to sweep.
"$BIN" inspect -i "$OUT/lab.gfa" -c "$OUT/lab.csv" --bubble-id 1 -o "$OUT/ex" --table-out "$OUT/ex/bubble_1/ex.bubble_1.clusters.tsv" -q >/dev/null 2>&1
[ -s "$OUT/ex/bubble_1/ex.bubble_1.clusters.tsv" ] && ok "an output this run wrote is never removed as stale" \
  || bad "the run removed its own --table-out"
# Files compress by name.
"$BIN" inspect -i "$OUT/lab.gfa" -c "$OUT/lab.csv" --bubble-id 1 -o "$OUT/gz" --fasta-out "$OUT/gz.fa" \
       --table-out "$OUT/gz.tsv.gz" -q >/dev/null 2>&1
[ "$(head -c 1 "$OUT/gz.fa" 2>/dev/null)" = ">" ] && gzip -t "$OUT/gz.tsv.gz" 2>/dev/null \
  && ok "outputs are gzip exactly when named .gz" || bad "compression does not follow the output names"
# A rerun over a smaller site table removes the files of the bubble it no longer has.
cp "$OUT/lab.csv" "$OUT/two.csv"; printf '2,3,+,1,+,1,5,5,"4"\n' >> "$OUT/two.csv"
"$BIN" inspect -i "$OUT/lab.gfa" -c "$OUT/two.csv" -o "$OUT/fam" --cluster -q >/dev/null 2>&1
had=$(ls "$OUT"/fam/bubble_2/ 2>/dev/null | wc -l | tr -d ' ')
printf 'old\n' > "$OUT/fam.bubble_1.node_counts.tsv"   # a site file beside the prefix
"$BIN" inspect -i "$OUT/lab.gfa" -c "$OUT/lab.csv" -o "$OUT/fam" -q >/dev/null 2>&1
[ "$had" = "6" ] && [ ! -e "$OUT/fam/bubble_2" ] && [ ! -e "$OUT/fam/bubble_1/fam.bubble_1.clusters.tsv" ] \
  && [ ! -e "$OUT/fam.bubble_1.node_counts.tsv" ] && [ -s "$OUT/fam/bubble_1/fam.bubble_1.node_counts.tsv" ] \
  && ok "a rerun at the prefix leaves only its own family (bubble 2's folder, clusters and flat files removed)" \
  || bad "stale family: had $had bubble-2 files, left: $(ls -d "$OUT"/fam/bubble_2 "$OUT"/fam/bubble_1/fam.bubble_1.clusters.tsv "$OUT"/fam.bubble_1.node_counts.tsv 2>/dev/null | paste -sd' ' -)"
# Repeatable --bubble-id; a CSV with id 0 or a repeated id is refused.
"$BIN" inspect -i "$OUT/lab.gfa" -c "$OUT/two.csv" --bubble-id 1 --bubble-id 2 -o "$OUT/rep" -q >/dev/null 2>&1
[ -s "$OUT/rep/bubble_1/rep.bubble_1.node_counts.tsv" ] && [ -s "$OUT/rep/bubble_2/rep.bubble_2.node_counts.tsv" ] \
  && ok "--bubble-id is repeatable" || bad "--bubble-id 1 --bubble-id 2 did not inspect both"
sed 's/^2,/0,/' "$OUT/two.csv" > "$OUT/zero.csv"; sed 's/^2,/1,/' "$OUT/two.csv" > "$OUT/dup.csv"
z=0; for c in zero dup; do "$BIN" inspect -i "$OUT/lab.gfa" -c "$OUT/$c.csv" -o "$OUT/$c" -q >/dev/null 2>&1 || z=$((z + 1)); done
[ "$z" = 2 ] && ok "a bubbles CSV with id 0 or a repeated id is refused" || bad "bad CSV ids accepted ($z of 2 refused)"
# clusters.tsv separates members with ';', so such a name is refused when clustering.
sed 's/^P\tx@1/P\tx;y/' "$OUT/lab.gfa" > "$OUT/semi.gfa"
"$BIN" inspect -i "$OUT/semi.gfa" -c "$OUT/lab.csv" -o "$OUT/semi" --cluster -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a ';' in a path name is refused with --cluster" || bad "a ';' path name was listed in clusters.tsv"
# A prefix ending in '/' or given as '.' still names its files after the folder: none starts with a dot,
# which the plot script's file search would skip.
"$BIN" inspect -i "$OUT/lab.gfa" -c "$OUT/lab.csv" -o "$OUT/slash/" --genotype-table -q >/dev/null 2>&1
mkdir -p "$OUT/here"
(cd "$OUT/here" && "$BIN" inspect -i "$OUT/lab.gfa" -c "$OUT/lab.csv" -o . -q >/dev/null 2>&1)
[ -s "$OUT/slash/slash.sites.tsv" ] && [ -s "$OUT/slash/bubble_1/slash.bubble_1.node_counts.tsv" ] \
  && [ -s "$OUT/here/bubble_1/here.bubble_1.node_counts.tsv" ] \
  && [ -z "$(find "$OUT/slash" "$OUT/here" -name '.*' -type f)" ] \
  && ok "a prefix 'dir/' or '.' names its files after the folder" \
  || bad "prefix 'dir/' or '.': $(find "$OUT/slash" "$OUT/here" -type f | sed "s|$OUT/||" | paste -sd' ' -)"
# ---- --genotype-table --------------------------------------------------------------------------------
# Site S..T: interior A (reference) or B (alternative); hD deletes the interior; hX does not cross the site.
# h2 crosses twice -- its allele is both crossings. The reference walks the site sink -> source (P line on
# the minus strand), so under -r the alleles are reverse-complemented and the nodes placed on it.
gt_check() { [ "$2" = "$3" ] && ok "$1 = $3" || bad "$1: expected $3, got ${2:-<missing>}"; }
{ printf 'H\tVN:Z:1.0\nS\tS\tACGTACGTAA\nS\tA\tGGGGCCCC\nS\tB\tTTTTAAAAT\nS\tT\tCATGCATGCC\nS\tZ\tAAAAAAAAAA\n'
  for e in "S + A +" "A + T +" "S + B +" "B + T +" "S + T +" "T + S +" "Z + T +"; do set -- $e; printf 'L\t%s\t%s\t%s\t%s\t0M\n' "$1" "$2" "$3" "$4"; done
  printf 'P\tref#0#chr1:100-128\tT-,A-,S-\t*\n'
  printf 'P\th1#1#chr1\tS+,B+,T+\t*\n'
  printf 'P\th1#2#chr1\tS+,A+,T+\t*\n'
  printf 'P\thD#1#chr1\tS+,T+\t*\n'
  printf 'P\th2#1#chr1\tS+,A+,T+,S+,B+,T+\t*\n'
  printf 'P\thX#1#chr1\tZ+,T+\t*\n'; } > "$OUT/gt.gfa"
printf 'bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,S,+,T,+,5,0,9,"A;B"\n' > "$OUT/gt.csv"
printf 's1\th1#1#chr1\th1#2#chr1\ns2\thD#1#chr1\thX#1#chr1\n' > "$OUT/gt.samples.tsv"
"$BIN" inspect -i "$OUT/gt.gfa" -c "$OUT/gt.csv" -o "$OUT/gt" --genotype-table -r ref -q >/dev/null 2>&1
GT="$OUT/gt/bubble_1/gt.bubble_1"
gt_check "haplotype rows: one allele each, NA for a path not crossing (h1#1 h2 hX)" \
  "$(gunzip -c "$GT.genotypes.tsv.gz" | awk -F'\t' '$1=="h1#1#chr1"||$1=="h2#1#chr1"||$1=="hX#1#chr1"{s=0; na=0; for(i=2;i<=NF;i++){if($i=="NA")na=1; else s+=$i} printf "%s ", (na?"NA":s)}')" "1 1 NA "
gt_check "four alleles: A, B, the deletion, and h2's two crossings" "$(gunzip -c "$GT.alleles.fa.gz" | grep '^>' | cut -d' ' -f1 | sort -u | wc -l | tr -d ' ')" "4"
gt_check "a two-crossing allele is written a record per crossing" \
  "$(gunzip -c "$GT.alleles.fa.gz" | grep -c 'h2#1#chr1 crossing=')" "2"
gt_check "under -r (reference sink -> source) alleles are on the reference strand: A reads as its reverse complement" \
  "$(gunzip -c "$GT.alleles.fa.gz" | grep -A1 'ref#0#chr1' | tail -1)" "GGCATGCATGGGGGCCCCTTACGTACGT"
gt_check "the node dictionary places the reference's nodes on its interval" \
  "$(gunzip -c "$GT.nodes.json.gz" | grep -o '"reference_strand": "[a-z_]*"\|"A": {[^}]*"ref_start": [0-9]*' | sed 's/.*: //' | paste -sd' ' -)" \
  "\"sink_to_source\" 111"
"$BIN" inspect -i "$OUT/gt.gfa" -c "$OUT/gt.csv" -o "$OUT/gts" --samples "$OUT/gt.samples.tsv" -q >/dev/null 2>&1
gt_check "sample rows count alleles over the sample's haplotypes; a sample with one haplotype not crossing is NA" \
  "$(gunzip -c "$OUT/gts/bubble_1/gts.bubble_1.genotypes.tsv.gz" | awk -F'\t' 'NR>1{s=0; na=0; for(i=2;i<=NF;i++){if($i=="NA")na=1; else s+=$i} printf "%s=%s ", $1, (na?"NA":s)}')" "s1=2 s2=NA "
gt_check "sites.tsv: 6 haplotypes, 5 crossing, 1 NA, none touching only, the reference crosses, 4 alleles" \
  "$(awk -F'\t' 'NR==2{print $4, $5, $6, $7, $8, $9, $10, $11}' "$OUT/gt/gt.sites.tsv")" "6 5 1 0 1 4 NA NA"
gt_check "sites.tsv with --samples: 2 samples, 1 NA; no -r, so the reference column is NA" \
  "$(awk -F'\t' 'NR==2{print $8, $10, $11}' "$OUT/gts/gts.sites.tsv")" "NA 2 1"
"$BIN" inspect -i "$OUT/touch.gfa" -c "$OUT/touch.csv" -o "$OUT/touchg" --genotype-table -q >/dev/null 2>&1
gt_check "sites.tsv counts paths touching both boundaries without crossing among the NA rows" \
  "$(awk -F'\t' 'NR==2{print $5, $6, $7}' "$OUT/touchg/touchg.sites.tsv")" "1 2 2"
"$BIN" inspect -i "$OUT/gt.gfa" -c "$OUT/gt.csv" -o "$OUT/gts" -q >/dev/null 2>&1
gt_check "a run without --genotype-table removes the sites.tsv an earlier run left" \
  "$([ -e "$OUT/gts/gts.sites.tsv" ] && echo kept || echo removed)" "removed"
gt_check "an allele crossing twice has a walk record per crossing, as in the allele FASTA" \
  "$(gunzip -c "$GT.alleles.walks.fa.gz" | grep -c 'h2#1#chr1 crossing=')" "2"
"$BIN" inspect -i "$OUT/gt.gfa" -c "$OUT/gt.csv" -o "$OUT/gtx" --genotype-table -r '' -q >/dev/null 2>&1
gt_check "an empty -r is refused" "$?" "2"
# A samples table, or a file of the user's, under a name inspect sweeps is not removed: the sweep takes
# sites.tsv only when it starts with inspect's header.
mkdir -p "$OUT/gtk/bubble_1"; cp "$OUT/gt.samples.tsv" "$OUT/gtk/bubble_1/gtk.bubble_1.clusters.tsv"
"$BIN" inspect -i "$OUT/gt.gfa" -c "$OUT/gt.csv" -o "$OUT/gtk" --samples "$OUT/gtk/bubble_1/gtk.bubble_1.clusters.tsv" -q >/dev/null 2>&1
gt_check "the sweep keeps the samples table" "$([ -e "$OUT/gtk/bubble_1/gtk.bubble_1.clusters.tsv" ] && echo kept || echo removed)" "kept"
printf 'chrom\tpos\n' > "$OUT/gtk/gtk.sites.tsv"
"$BIN" inspect -i "$OUT/gt.gfa" -c "$OUT/gt.csv" -o "$OUT/gtk" -q >/dev/null 2>&1
gt_check "the sweep keeps a sites.tsv inspect did not write" "$(cat "$OUT/gtk/gtk.sites.tsv" 2>/dev/null)" "$(printf 'chrom\tpos')"
# A site table without orientation columns is read from each path's first boundary to its last: h2 crosses
# S+..T+ twice and gets one window, not a row per crossing found by guessing '+' orientations.
printf 'bubble_id,source,sink,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,S,T,5,0,9,"A;B"\n' > "$OUT/gtn.csv"
"$BIN" inspect -i "$OUT/gt.gfa" -c "$OUT/gtn.csv" -o "$OUT/gtn" -q >/dev/null 2>&1
gt_check "without orientation columns a path crossing twice is one first-to-last row" \
  "$(gunzip -c "$OUT/gtn/bubble_1/gtn.bubble_1.paths.fa.gz" | grep '^>h2' | cut -d' ' -f1 | paste -sd' ' -)" ">h2#1#chr1"
# Under -r, a multi-copy allele lists its copies in the order the reference meets them, however the site
# table orients the site: here written against the paths (4-..1-), the reference still reads 2 before 3.
{ printf 'H\tVN:Z:1.0\nS\t1\tAAAAAAAAAA\nS\t2\tCCCC\nS\t3\tGGGGGG\nS\t4\tTTTTTTTTTT\nS\t5\tACGTACGT\n'
  for e in "1 2" "1 3" "2 4" "3 4" "4 5" "5 1"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref\t1+,2+,4+,5+,1+,3+,4+\t*\nP\talt\t1+,3+,4+,5+,1+,2+,4+\t*\nP\trefrev\t4-,3-,1-,5-,4-,2-,1-\t*\n'; } > "$OUT/ro.gfa"
printf 'bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,4,-,1,-,3,4,6,"2;3"\n' > "$OUT/ro.csv"
"$BIN" inspect -i "$OUT/ro.gfa" -c "$OUT/ro.csv" -o "$OUT/ro" --genotype-table -r ref -q >/dev/null 2>&1
"$BIN" inspect -i "$OUT/ro.gfa" -c "$OUT/ro.csv" -o "$OUT/ror" --genotype-table -r refrev -q >/dev/null 2>&1
gt_check "-r keeps the reference's own copy order (ref; refrev, stored reversed)" \
  "$(gunzip -c "$OUT/ro/bubble_1/ro.bubble_1.alleles.walks.fa.gz" | grep -A1 '^>allele1' | grep -v '^>' | paste -sd' ' -); $(gunzip -c "$OUT/ror/bubble_1/ror.bubble_1.alleles.walks.fa.gz" | grep -A1 '^>allele1' | grep -v '^>' | paste -sd' ' -)" \
  "1+,2+,4+ 1+,3+,4+; 4-,3-,1- 4-,2-,1-"
# -r matches as written before it matches in any case: "ref" names ref#0#..., not Ref#9#... as well.
{ cat "$OUT/gt.gfa"; printf 'P\tRef#9#chr1\tS+,A+,T+\t*\n'; } > "$OUT/gtcase.gfa"
"$BIN" inspect -i "$OUT/gtcase.gfa" -c "$OUT/gt.csv" -o "$OUT/gtcase" --genotype-table -r ref -q >/dev/null 2>&1
gt_check "-r ref beside Ref#9: the path matching as written" \
  "$(gunzip -c "$OUT/gtcase/bubble_1/gtcase.bubble_1.nodes.json.gz" 2>/dev/null | grep -o '"reference_path": "[^"]*"' | head -1)" \
  '"reference_path": "ref#0#chr1:100-128"'
# A flag's earlier name still works, with a warning naming the new one.
cp "$OUT/gt.csv" "$OUT/gtold.bubbles.csv"
e=$("$BIN" inspect -i "$OUT/gt.gfa" --bubble-prefix-in "$OUT/gtold" -o "$OUT/gtold" -q 2>&1 >/dev/null)
gt_check "--bubble-prefix-in is read as --bubble-prefix, with a warning" \
  "$([ -s "$OUT/gtold/bubble_1/gtold.bubble_1.paths.fa.gz" ] && echo read) $(printf '%s' "$e" | grep -q 'is now --bubble-prefix' && echo warned)" \
  "read warned"
gt_check "the walks FASTA spells what the paths FASTA spells (one walk per crossing record)" \
  "$(gunzip -c "$OUT/gt/bubble_1/gt.bubble_1.walks.fa.gz" | grep -c '^>') $(gunzip -c "$OUT/gt/bubble_1/gt.bubble_1.paths.fa.gz" | grep -c '^>')" "6 6"
"$BIN" inspect -i "$OUT/gt.gfa" -c "$OUT/gt.csv" -o "$OUT/gtr" -r ref -q >/dev/null 2>&1
gt_check "-r without --genotype-table is refused" "$?" "2"

# ---- inversion sites. A 9.7 kb inversion between collapsed inverted repeats (hI, hJ) holds a site
# where h1 and hI lose X. bubble makes the inversion an inversion site enclosing that site (see call_stats).
INVNEST="$OUT/invnest"; mkdir -p "$INVNEST"
cat > "$INVNEST.spec" <<'SPEC'
N S 2000 1
N I 300 2
N M1 3000 3
N Ma 1500 4
N X 60 8
N Mb 1500 9
N M3 3000 5
N T 2000 6
P ref#0#chr1:0-13660 S+,I+,M1+,Ma+,X+,Mb+,M3+,I-,T+
P h1#1#chr1 S+,I+,M1+,Ma+,Mb+,M3+,I-,T+
P hI#1#chr1 S+,I+,M3-,Mb-,Ma-,M1-,I-,T+
P hJ#1#chr1 S+,I+,M3-,Mb-,X-,Ma-,M1-,I-,T+
P h2#1#chr1 S+,I+,M1+,Ma+,X+,Mb+,M3+,I-,T+
SPEC
"${PYTHON:-python3}" -I "$(cd "$(dirname "$0")/.." && pwd)/tools/spec_gfa.py" "$INVNEST.spec" "$INVNEST" >/dev/null
"$BIN" bubble -i "$INVNEST/g.gfa" -r 'ref#0#chr1:0-13660' -o "$INVNEST/b" -q >/dev/null 2>&1
"$BIN" inspect -i "$INVNEST/b.sorted.gfa" -b "$INVNEST/b" -o "$INVNEST/i" --genotype-table -q >/dev/null 2>&1
# Shown on its own nodes, the inversion site's alleles are the two orientations of the stretch it turns.
gt_check "an inversion site's alleles are its two orientations, carriers hI and hJ" \
      "$(gunzip -c "$INVNEST/i/bubble_2/i.bubble_2.genotypes.tsv.gz" | awk -F'\t' 'NR > 1 && $3 == 1 { print $1 }' | paste -sd, -)" \
      "hI#1#chr1,hJ#1#chr1"
# The inner site's deletion junction (its source straight to its sink) is that site's, not the inversion site's.
gt_check "the inner site's deletion edge is not counted at the inversion site" \
      "$(head -1 "$INVNEST/i/bubble_2/i.bubble_2.edge_counts.tsv" | tr '\t' '\n' | grep -c '^edge\.4+>6+$')" "0"

# Columns follow the graph's order: a sorted graph numbers nodes along the reference, so node 9 comes
# before node 10 (as text, "10" would sort first and a heatmap would read the site out of order).
ORD="$OUT/colorder"; mkdir -p "$ORD"
cat > "$ORD/spec" <<'EOF2'
N 8 60 1
N 9 70 2
N 10 55 3
N 11 60 4
P ref 8+,9+,11+
P alt 8+,10+,11+
P both 8+,9+,10+,11+
B 8 11 9;10
EOF2
python3 "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/tools/spec_gfa.py" "$ORD/spec" "$ORD"
"$BIN" inspect -i "$ORD/g.gfa" -b "$ORD/g" -o "$ORD/i" -q >/dev/null 2>&1
cols="$(head -1 "$ORD/i/bubble_1/i.bubble_1.node_counts.tsv" | cut -f5- | tr '\t' ' ')"
lens="$(tail -n +2 "$ORD/i/bubble_1/i.bubble_1.node_lengths.tsv" | cut -f1 | paste -sd' ' -)"
edges="$(head -1 "$ORD/i/bubble_1/i.bubble_1.edge_counts.tsv" | cut -f5- | tr '\t' ' ')"
[ "$cols" = "node.9 node.10" ] && [ "$lens" = "9 10" ] && [ "$edges" = "edge.8+>9+ edge.8+>10+ edge.9+>10+ edge.9+>11+ edge.10+>11+" ] \
  && ok "node and edge columns follow the graph's numeric order" \
  || bad "column order: nodes '$cols', lengths '$lens', edges '$edges'"

# A haplotype stored on the other strand meets its two crossings in reverse order; its allele must be the
# same as the forward copy's. h4 carries the two copies in the other order on the forward strand: a
# different allele.
STR="$OUT/strand"; mkdir -p "$STR"
cat > "$STR/spec" <<'EOF2'
N 1 =AAAAAAAAAC
N 2 =CCCCCCCCCG
N 3 =GGGGGGGGGT
N 4 =TTTTTTTTTA
N 5 =ACGTACGTAG
P h1 1+,2+,3+,4+,1+,5+,3+
P h2 3-,5-,1-,4-,3-,2-,1-
P h3 1+,2+,3+,4+,1+,5+,3+
P h4 1+,5+,3+,4+,1+,2+,3+
B 1 3 2;5
EOF2
python3 "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/tools/spec_gfa.py" "$STR/spec" "$STR"
"$BIN" inspect -i "$STR/g.gfa" -c "$STR/g.bubbles.csv" --bubble-id 1 --genotype-table -o "$STR/i" -q >/dev/null 2>&1
gt="$(gzip -dc "$STR/i/bubble_1/i.bubble_1.genotypes.tsv.gz" 2>/dev/null | tail -n +2 | tr '\t' ' ' | paste -sd';' -)"
[ "$gt" = "h1 1 0;h2 1 0;h3 1 0;h4 0 1" ] \
  && ok "a reverse-stored haplotype crossing twice has its forward copy's allele; copies swapped do not" \
  || bad "allele by strand: '$gt' (expected h1 1 0;h2 1 0;h3 1 0;h4 0 1)"
# A path name with a space would cut every FASTA id at the space: refused.
sed 's/^P\th3\t/P\th 3\t/' "$STR/g.gfa" > "$STR/sp.gfa"
"$BIN" inspect -i "$STR/sp.gfa" -c "$STR/g.bubbles.csv" -o "$STR/sp" -q >/dev/null 2>"$STR/sp.err"; rc=$?
[ "$rc" -ne 0 ] && grep -q "whitespace" "$STR/sp.err" \
  && ok "a path name with a space is refused" \
  || bad "path name with a space: exit $rc, $(head -1 "$STR/sp.err")"

echo
if [ "$fails" -eq 0 ]; then echo "inspect_stats: all assertions passed"; exit 0; fi
echo "inspect_stats: $fails assertion(s) FAILED"; exit 1
