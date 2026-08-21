#!/usr/bin/env bash
# inspect_stats.sh - contract assertions for `panvar inspect`.
#
# real_smoke.sh only checks that inspect's files exist and that the gzip is readable, which passes
# whether or not the contents are right -- it would not have caught any of the defects below: a whole
# allele missing from the FASTA, a graph with duplicate path names accepted, or a bubbles CSV from a
# different graph running to completion with a header-only matrix.
#
# Every fixture here is small enough that the expected output is derivable by hand.
#
#   inspect_stats.sh <panvar-binary> <out-dir>
set -uo pipefail

BIN="${1:?usage: inspect_stats.sh <panvar> <outdir>}"
OUT="${2:?}"
mkdir -p "$OUT"
OUT="$OUT/run.$$.$(date +%s)"
rm -rf "$OUT"; mkdir -p "$OUT"
fails=0
ok()  { printf "  ok   %s\n" "$1"; }
bad() { printf "  FAIL %s\n" "$1"; fails=$((fails + 1)); }
heads() { gzip -dc "$1" 2>/dev/null | grep -c '^>' | tr -d ' '; }
hdr()   { gzip -dc "$1" 2>/dev/null | grep '^>' | sed 's/^>//'; }

# ---------------------------------------------------------------- a pure deletion is an allele
# `full` carries interior node 2; `del` goes straight from source to sink. Both are real alleles of the
# same site and `bubble` reports distinct_alleles=2, but inspect used the inside-node-only interval
# finder, which has no interval for a walk with no interior -- so it emitted one sequence of two.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGG\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tfull\t1+,2+,3+\t*\nP\tdel\t1+,3+\t*\n'; } > "$OUT/del.gfa"
"$BIN" bubble -i "$OUT/del.gfa" -r full -o "$OUT/db" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" inspect -i "$OUT/db.sorted.gfa" -c "$OUT/db.bubbles.csv" --bubble-id 1 -o "$OUT/d" -q >/dev/null 2>&1

[ "$(heads "$OUT/d.bubble_1.paths.fa.gz")" = "2" ] \
  && ok "both alleles are emitted, including the one with no interior node" \
  || bad "expected 2 sequences, got $(heads "$OUT/d.bubble_1.paths.fa.gz")"
hdr "$OUT/d.bubble_1.paths.fa.gz" | grep -q '^del ' \
  && ok "the pure-deletion haplotype appears by name" \
  || bad "the deletion haplotype is missing from the FASTA"
# 1 + 3 = 20 bp with node 2 (10 bp) deleted; the full walk is 30.
hdr "$OUT/d.bubble_1.paths.fa.gz" | grep -q '^del .*length_bp=20' \
  && ok "the deletion allele spells 20 bp (source + sink, interior dropped)" \
  || bad "deletion length wrong: $(hdr "$OUT/d.bubble_1.paths.fa.gz" | grep '^del ')"
# The fallback still has to report interval metadata, not leave it at whatever was last set.
hdr "$OUT/d.bubble_1.paths.fa.gz" | grep -q '^del .*interval=0-1' \
  && ok "the empty-interior fallback reports its own interval" \
  || bad "deletion interval metadata wrong: $(hdr "$OUT/d.bubble_1.paths.fa.gz" | grep '^del ')"
# The count matrix must gain a row for it, reading zero for the deleted node.
[ "$(awk 'NR>1' "$OUT/d.bubble_1.node_counts.tsv" | wc -l | tr -d ' ')" = "2" ] \
  && ok "the count matrix has a row per allele" \
  || bad "count matrix has $(awk 'NR>1' "$OUT/d.bubble_1.node_counts.tsv" | wc -l | tr -d ' ') rows, expected 2"
awk '$1=="del"{print $3}' "$OUT/d.bubble_1.node_counts.tsv" | grep -q '^0:0:0$' \
  && ok "the deleted node reads 0:0:0 for the carrier of the deletion" \
  || bad "deleted node count is $(awk '$1=="del"{print $3}' "$OUT/d.bubble_1.node_counts.tsv"), expected 0:0:0"

# ---------------------------------------------------------------- reverse traversal
# `rev` crosses the same site on the opposite strand. The walk is canonicalized source->sink, so it must
# spell the same 30 bp as the forward haplotype and be reported as sink->source.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGG\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tfull\t1+,2+,3+\t*\nP\tdel\t1+,3+\t*\nP\trev\t3-,2-,1-\t*\n'; } > "$OUT/rev.gfa"
"$BIN" bubble -i "$OUT/rev.gfa" -r full -o "$OUT/rb" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" inspect -i "$OUT/rb.sorted.gfa" -c "$OUT/rb.bubbles.csv" --bubble-id 1 -o "$OUT/r" -q >/dev/null 2>&1
[ "$(heads "$OUT/r.bubble_1.paths.fa.gz")" = "3" ] \
  && ok "a reverse traversal is emitted alongside the forward ones" \
  || bad "expected 3 sequences, got $(heads "$OUT/r.bubble_1.paths.fa.gz")"
hdr "$OUT/r.bubble_1.paths.fa.gz" | grep -q '^rev .*source_to_sink=0' \
  && ok "a reverse traversal is reported as sink->source" \
  || bad "reverse traversal metadata wrong: $(hdr "$OUT/r.bubble_1.paths.fa.gz" | grep '^rev ')"
hdr "$OUT/r.bubble_1.paths.fa.gz" | grep -q '^rev .*length_bp=30' \
  && ok "a reverse traversal spells the same length as the forward walk" \
  || bad "reverse traversal length wrong: $(hdr "$OUT/r.bubble_1.paths.fa.gz" | grep '^rev ')"

# ---------------------------------------------------------------- a repeated node inside the site
# A boundary cannot repeat within a discovered snarl -- a boundary carrying a cycle stops being a cut
# vertex, and the decomposition then emits no snarl there at all (checked: the fixture below with the
# cycle moved onto node 3..4 yields zero bubbles). The reachable case is a repeated INTERIOR node,
# which is what a tandem array is: `rep` walks 1,2,3,2,4 and visits node 2 twice.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGG\nS\t4\tTTTTTTTTTT\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t1\t+\t4\t+\t0M\n'
  printf 'L\t2\t+\t3\t+\t0M\nL\t3\t+\t2\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'P\tfull\t1+,2+,4+\t*\nP\tdel\t1+,4+\t*\nP\trep\t1+,2+,3+,2+,4+\t*\n'; } > "$OUT/rep.gfa"
"$BIN" bubble -i "$OUT/rep.gfa" -r full -o "$OUT/pb" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" inspect -i "$OUT/pb.sorted.gfa" -c "$OUT/pb.bubbles.csv" -o "$OUT/p" -q >/dev/null 2>&1
[ "$(heads "$OUT/p.bubble_1.paths.fa.gz")" = "3" ] \
  && ok "all three alleles are emitted when one revisits an interior node" \
  || bad "expected 3 sequences at the repeated-node site, got $(heads "$OUT/p.bubble_1.paths.fa.gz")"
[ "$(hdr "$OUT/p.bubble_1.paths.fa.gz" | grep -c 'length_bp=0')" = "0" ] \
  && ok "no allele spells an empty sequence at a repeated interior node" \
  || bad "a zero-length walk was emitted at a repeated interior node"
# The revisited node must be COUNTED twice, not deduplicated: multiplicity is the copy number.
awk '$1=="rep"{print $3}' "$OUT/p.bubble_1.node_counts.tsv" | grep -q '^2:2:0$' \
  && ok "a revisited interior node is counted twice (2:2:0), not collapsed" \
  || bad "revisited node count is $(awk '$1=="rep"{print $3}' "$OUT/p.bubble_1.node_counts.tsv"), expected 2:2:0"

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
# Nothing ties the CSV to the GFA, so a CSV from another graph used to produce a header-only matrix,
# node lengths of 0, and exit 0.
printf 'bubble_id,source,sink,path_support,min_inside_bp,max_inside_bp,inside_nodes\n99,900,902,2,10,10,"901"\n' \
  > "$OUT/ghost.csv"
rm -f "$OUT/gh".*
"$BIN" inspect -i "$OUT/del.gfa" -c "$OUT/ghost.csv" --bubble-id 99 -o "$OUT/gh" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a bubbles CSV naming nodes absent from the GFA is refused" \
               || bad "a mismatched CSV ran to completion"
[ "$(ls "$OUT/gh".* 2>/dev/null | wc -l | tr -d ' ')" = "0" ] \
  && ok "a refused run leaves no partial per-bubble outputs" \
  || bad "mismatched run left $(ls "$OUT/gh".* 2>/dev/null | wc -l | tr -d ' ') file(s)"

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
# The sketch holds 512 shingles, so only a walk with more than that is truncated -- and truncation is
# the whole reason the estimator matters. An earlier version of this fixture had 60 shingles and
# therefore tested nothing: both estimators are exact when nothing is dropped.
#
# A cyclic array: node 2 carries a self-loop, walked 3000 times by `long` and 500 by `short`, so the
# shingle sets are 3000 and 500 with the smaller nested in the larger. True Jaccard is 500/3000 =
# 0.1667, i.e. an identity of 2J/(1+J) = 0.286.
#
#   old estimator (intersection over union of two truncated sketches): J 0.0847, identity 0.156
#   new estimator (bottom-k of the union):                             J 0.1580, identity 0.273
#
# The assertions below bracket the identity in [0.22, 0.30], which the old estimator's 0.156 cannot
# satisfy -- so this fails if the biased form is ever restored.
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
[ "$(nclust "$OUT/a22.bubble_1.clusters.tsv")" = "1" ] \
  && ok "truncated sketches at a 6x length ratio still cluster at 0.22 (needs the unbiased estimator)" \
  || bad "expected 1 cluster at 0.22, got $(nclust "$OUT/a22.bubble_1.clusters.tsv"); the biased estimator reads 0.156 here"

"$BIN" inspect -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" --cluster --cluster-similarity 0.30 \
       -o "$OUT/a30" -q >/dev/null 2>&1
[ "$(nclust "$OUT/a30.bubble_1.clusters.tsv")" = "2" ] \
  && ok "and separate at 0.30, bracketing the estimate at the true 0.286" \
  || bad "expected 2 clusters at 0.30, got $(nclust "$OUT/a30.bubble_1.clusters.tsv")"

# ------------------------------------------- COMPLETE sketches are compared exactly, not estimated
# The bottom-k estimator subsamples even when both sketches already hold every shingle the walks have,
# so it answers with sampling error where the exact value is sitting in memory. This is the worked
# example from docs/algorithms/inspect.md: two and three copies of a 3-node unit between the same
# boundaries. W1 has 6 shingles, W3 has 9, they share 6, so the multiset Jaccard is 6/9 = 0.667 and
# the identity is 2(0.667)/1.667 = 0.80. The estimator reads 0.5 and therefore 0.667.
#
# 0.79 and 0.85 bracket the exact identity, so the pair must cluster at the first and separate at the
# second. With the estimator both readings are "separate", which is what makes this discriminating.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\tAAAAAAAAAAAAAAAAAAAA\nS\tA\tCCCCCCCCCC\nS\tB\tGGGGGGGGGG\nS\tC\tTTTTTTTTTT\n'
  printf 'S\tG\tAAAAAAAAAAAAAAAAAAAA\n'
  printf 'L\tF\t+\tA\t+\t0M\nL\tA\t+\tB\t+\t0M\nL\tB\t+\tC\t+\t0M\nL\tC\t+\tA\t+\t0M\nL\tC\t+\tG\t+\t0M\n'
  printf 'P\th1\tF+,A+,B+,C+,A+,B+,C+,G+\t*\n'
  printf 'P\th3\tF+,A+,B+,C+,A+,B+,C+,A+,B+,C+,G+\t*\n'; } > "$OUT/exact.gfa"
"$BIN" bubble -i "$OUT/exact.gfa" -r h1 -o "$OUT/ex" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" inspect -i "$OUT/ex.sorted.gfa" -c "$OUT/ex.bubbles.csv" --cluster --cluster-similarity 0.79 \
       -o "$OUT/ex79" -q >/dev/null 2>&1
[ "$(nclust "$OUT/ex79.bubble_1.clusters.tsv")" = "1" ] \
  && ok "two complete sketches cluster at 0.79, so the exact 0.80 identity was used" \
  || bad "expected 1 cluster at 0.79, got $(nclust "$OUT/ex79.bubble_1.clusters.tsv"); the estimator reads 0.667 here"
"$BIN" inspect -i "$OUT/ex.sorted.gfa" -c "$OUT/ex.bubbles.csv" --cluster --cluster-similarity 0.85 \
       -o "$OUT/ex85" -q >/dev/null 2>&1
[ "$(nclust "$OUT/ex85.bubble_1.clusters.tsv")" = "2" ] \
  && ok "and separate at 0.85, bracketing the exact identity rather than an estimate of it" \
  || bad "expected 2 clusters at 0.85, got $(nclust "$OUT/ex85.bubble_1.clusters.tsv")"

# ------------------------------------------- the representative does not follow GFA record order
# Identical walks are pooled in encounter order and the first member became `representative_path`, so
# reordering otherwise identical P records changed the reported representative and, through it, the
# sort key of the whole TSV. The rule is now the lexicographically smallest member of the medoid walk.
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
rep1=$(awk -F'\t' 'NR>1 && $2==2{print $3}' "$OUT/repord_i.bubble_1.clusters.tsv")
rep2=$(awk -F'\t' 'NR>1 && $2==2{print $3}' "$OUT/repord2_i.bubble_1.clusters.tsv")
[ -n "$rep1" ] && [ "$rep1" = "$rep2" ] && [ "$rep1" = "aaa" ] \
  && ok "the representative is order-independent and lexicographically smallest (aaa)" \
  || bad "representative depends on P-record order: '$rep1' vs '$rep2' (expected aaa both ways)"

# ------------------------------------------- the two second-pass branches, pinned directly
# Both were implemented in the second review pass and neither had a focused assertion, so a
# regression in either would have shown up only as a wrong result somewhere downstream.

# 1. A DERIVED output collides with an explicit one. The preflight has to run AFTER bubble selection,
#    because the derived names contain the bubble ids; before that it saw only the explicit flags and
#    an explicit --table-out naming <prefix>.bubble_N.node_lengths.tsv silently clobbered one of them.
"$BIN" inspect -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" --bubble-id 1 \
       --table-out "$OUT/coll.bubble_1.node_lengths.tsv" -o "$OUT/coll" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "an explicit output colliding with a DERIVED one is refused" \
               || bad "an explicit --table-out naming a derived output was accepted"

# 2. The bubbles CSV names nodes that DO exist in the graph, but no path crosses the pair. That used
#    to produce a header-only count table and exit 0, which reads as "this site has no variation"
#    rather than "this CSV does not describe this graph".
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
warn=$("$BIN" inspect -i "$OUT/two.gfa" -c "$OUT/b3.bubbles.csv" -o "$OUT/stale" 2>&1 | grep -c "path_support=3")
[ "$warn" -ge 1 ] && ok "a stale panel (2 crossings against path_support=3) is reported" \
                  || bad "no warning when the graph carries fewer paths than the CSV records"
# ...and a matching panel must stay quiet, or the warning is noise.
quiet=$("$BIN" inspect -i "$OUT/b3.sorted.gfa" -c "$OUT/b3.bubbles.csv" -o "$OUT/fresh" 2>&1 | grep -c "path_support=")
[ "$quiet" = "0" ] && ok "and a matching panel produces no such warning" \
                   || bad "the stale-panel warning fires on a matching panel"

echo
if [ "$fails" -eq 0 ]; then echo "inspect_stats: all assertions passed"; exit 0; fi
echo "inspect_stats: $fails assertion(s) FAILED"; exit 1
