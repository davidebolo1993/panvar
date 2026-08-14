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

echo
if [ "$fails" -eq 0 ]; then echo "inspect_stats: all assertions passed"; exit 0; fi
echo "inspect_stats: $fails assertion(s) FAILED"; exit 1
