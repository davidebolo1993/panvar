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

# The paths-walk companion is also derived from an explicit FASTA name, and the genotype mode adds
# two more derived companions. They belong in the same preflight: StagedOutputs is transactional, but
# registering the same destination twice cannot define which logical output should win.
"$BIN" inspect -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" --bubble-id 1 \
       --fasta-out "$OUT/collwalk.fa.gz" --table-out "$OUT/collwalk.walks.fa.gz" \
       -o "$OUT/collwalk" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "an explicit output colliding with the derived paths walk is refused" \
               || bad "an explicit output naming the derived paths walk was accepted"
"$BIN" inspect -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" --bubble-id 1 \
       --table-out "$OUT/collgt.bubble_1.nodes.json.gz" --genotype-table \
       -o "$OUT/collgt" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "an explicit output colliding with a derived genotype companion is refused" \
               || bad "an explicit output naming the derived node dictionary was accepted"

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

# ------------------------------------------------------------- --genotype-table / --samples
# full and full2 spell the same allele; rev crosses on the opposite strand, so once canonicalized it
# spells the same 30 bp and must share that allele too; del is the pure deletion; other never crosses.
# Node 3 ends in T so a reverse complement is distinguishable from the forward spelling.
{ printf 'H\tVN:Z:1.0\nS\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGT\n'
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\n'
  printf 'P\tfull\t1+,2+,3+\t*\nP\tfull2\t1+,2+,3+\t*\nP\tdel\t1+,3+\t*\nP\trev\t3-,2-,1-\t*\nP\tother\t2+\t*\n'; } > "$OUT/gt.gfa"
"$BIN" bubble -i "$OUT/gt.gfa" -r full -o "$OUT/gb" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" inspect -i "$OUT/gb.sorted.gfa" -c "$OUT/gb.bubbles.csv" --bubble-id 1 -o "$OUT/g0" -q >/dev/null 2>&1
"$BIN" inspect -i "$OUT/gb.sorted.gfa" -c "$OUT/gb.bubbles.csv" --bubble-id 1 -o "$OUT/g" --genotype-table -q >/dev/null 2>&1
GT="$OUT/g.bubble_1.genotypes.tsv.gz"; AL="$OUT/g.bubble_1.alleles.fa.gz"
[ "$(heads "$AL")" = "2" ] && ok "identical sequences share an allele (5 paths, 2 alleles)" \
                          || bad "expected 2 alleles, got $(heads "$AL")"
hdr "$AL" | grep -qx 'allele1 full,full2,rev' \
  && ok "allele1 is the 3-haplotype allele and lists its haplotypes, the reverse walk included" \
  || bad "allele1 header wrong: $(hdr "$AL" | head -1)"
hdr "$AL" | grep -qx 'allele2 del' && ok "allele2 is the deletion" || bad "allele2 header wrong: $(hdr "$AL" | sed -n 2p)"
[ "$(gzip -dc "$GT" | head -1)" = "$(printf 'haplotype\tallele1\tallele2')" ] \
  && ok "table columns are the FASTA ids" || bad "table header wrong: $(gzip -dc "$GT" | head -1)"
[ "$(gzip -dc "$GT" | awk -F'\t' '$1=="del"{print $2 $3} $1=="rev"{print $2 $3}' | tr '\n' ' ')" = "01 10 " ] \
  && ok "haplotype rows are 0/1 indicators of the allele they spell" \
  || bad "haplotype rows wrong: $(gzip -dc "$GT")"
[ "$(gzip -dc "$GT" | awk -F'\t' '$1=="other"{print $2 $3}')" = "NANA" ] \
  && ok "a haplotype that does not cross the bubble reads NA" || bad "non-crossing row wrong: $(gzip -dc "$GT")"
[ "$(gzip -dc "$AL" | awk 'NR==2')" = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGT" ] \
  && ok "without --reference-path alleles are spelled source->sink" \
  || bad "allele1 sequence wrong: $(gzip -dc "$AL" | awk 'NR==2')"
cmp -s <(gzip -dc "$OUT/g0.bubble_1.paths.fa.gz") <(gzip -dc "$OUT/g.bubble_1.paths.fa.gz") \
  && cmp -s "$OUT/g0.bubble_1.node_counts.tsv" "$OUT/g.bubble_1.node_counts.tsv" \
  && ok "--genotype-table leaves the existing outputs byte-identical" \
  || bad "--genotype-table changed paths.fa.gz or node_counts.tsv"

"$BIN" inspect -i "$OUT/gb.sorted.gfa" -c "$OUT/gb.bubbles.csv" --bubble-id 1 -o "$OUT/grev" --genotype-table -r rev -q >/dev/null 2>&1
[ "$(gzip -dc "$OUT/grev.bubble_1.alleles.fa.gz" | awk 'NR==2')" = "ACCCCCCCCCGGGGGGGGGGTTTTTTTTTT" ] \
  && ok "a reference crossing sink->source reverse-complements every allele" \
  || bad "reference-oriented allele1 wrong: $(gzip -dc "$OUT/grev.bubble_1.alleles.fa.gz" | awk 'NR==2')"
"$BIN" inspect -i "$OUT/gb.sorted.gfa" -c "$OUT/gb.bubbles.csv" --bubble-id 1 -o "$OUT/gfwd" --genotype-table -r full2 -q >/dev/null 2>&1
[ "$(gzip -dc "$OUT/gfwd.bubble_1.alleles.fa.gz" | awk 'NR==2')" = "AAAAAAAAAACCCCCCCCCCGGGGGGGGGT" ] \
  && ok "a reference crossing source->sink keeps the spelling" \
  || bad "forward reference flipped allele1: $(gzip -dc "$OUT/gfwd.bubble_1.alleles.fa.gz" | awk 'NR==2')"
# Walks: every record of the paths FASTA and of the allele FASTA has a companion walk under the same
# header, and the walk spells exactly that record's sequence -- in whatever orientation -r gave it.
walks_spell() {   # walks_spell <gfa> <sequence fasta.gz> <walk fasta.gz>
  python3 - "$1" "$2" "$3" <<'PY'
import gzip, sys
nodes = {}
for line in open(sys.argv[1]):
    f = line.rstrip("\n").split("\t")
    if f[0] == "S": nodes[f[1]] = f[2]
def records(path):
    out, name = [], None
    for line in gzip.open(path, "rt"):
        line = line.rstrip("\n")
        if line.startswith(">"): name = line[1:]; out.append([name, ""])
        else: out[-1][1] += line
    return out
comp = str.maketrans("ACGT", "TGCA")
seqs, walks = records(sys.argv[2]), records(sys.argv[3])
headers_seq = [n for n, _ in seqs]; headers_walk = [n for n, _ in walks]
if headers_seq != headers_walk: print(f"headers differ: {headers_seq} vs {headers_walk}"); sys.exit(1)
for (name, seq), (_, walk) in zip(seqs, walks):
    spelled = "".join(nodes[t[:-1]] if t[-1] == "+" else nodes[t[:-1]][::-1].translate(comp) for t in walk.split(","))
    if spelled != seq: print(f"{name}: walk {walk} spells {spelled}, not {seq}"); sys.exit(1)
PY
}
msg=$(walks_spell "$OUT/gb.sorted.gfa" "$OUT/g.bubble_1.paths.fa.gz" "$OUT/g.bubble_1.walks.fa.gz") \
  && ok "each path's walk (1+,2+,3+ form) sits under its FASTA header and spells its sequence" || bad "path walks: $msg"
msg=$(walks_spell "$OUT/gb.sorted.gfa" "$OUT/g0.bubble_1.paths.fa.gz" "$OUT/g0.bubble_1.walks.fa.gz") \
  && ok "the path walks are written without --genotype-table too" || bad "path walks without --genotype-table: $msg"
msg=$(walks_spell "$OUT/gb.sorted.gfa" "$OUT/g.bubble_1.alleles.fa.gz" "$OUT/g.bubble_1.alleles.walks.fa.gz") \
  && ok "each allele's walk sits under the allele FASTA header and spells the allele" || bad "allele walks: $msg"
msg=$(walks_spell "$OUT/gb.sorted.gfa" "$OUT/grev.bubble_1.alleles.fa.gz" "$OUT/grev.bubble_1.alleles.walks.fa.gz") \
  && ok "a sink->source reference reverses the allele walks with the alleles" || bad "reversed allele walks: $msg"
[ ! -e "$OUT/g0.bubble_1.nodes.json.gz" ] && [ ! -e "$OUT/g0.bubble_1.alleles.walks.fa.gz" ] \
  && ok "the node dictionary and allele walks are written only with --genotype-table" \
  || bad "node dictionary or allele walks written without --genotype-table"
# Node dictionary: every walk node, its GFA sequence, visiting haplotypes, and reference position.
"$BIN" inspect -i "$OUT/gb.sorted.gfa" -c "$OUT/gb.bubbles.csv" --bubble-id 1 -o "$OUT/gdel" --genotype-table -r del -q >/dev/null 2>&1
msg=$(python3 - "$OUT/gb.sorted.gfa" "$OUT/gfwd.bubble_1.nodes.json.gz" "$OUT/gdel.bubble_1.nodes.json.gz" "$OUT/g.bubble_1.nodes.json.gz" <<'PY'
import gzip, json, sys
seq = {f[1]: f[2] for f in (l.rstrip("\n").split("\t") for l in open(sys.argv[1])) if f[0] == "S"}
by_seq = {v: k for k, v in seq.items()}
A, C, G = by_seq["AAAAAAAAAA"], by_seq["CCCCCCCCCC"], by_seq["GGGGGGGGGT"]
full, dele, bare = (json.load(gzip.open(p, "rt")) for p in sys.argv[2:5])
def fail(m): print(m); sys.exit(1)
for d in (full, dele, bare):
    if set(d["nodes"]) != {A, C, G}: fail(f"nodes {sorted(d['nodes'])}")
    for k, v in d["nodes"].items():
        if v["seq"] != seq[k] or v["length"] != len(seq[k]): fail(f"node {k} sequence wrong")
    if (d["nodes"][A]["haplotypes"], d["nodes"][C]["haplotypes"], d["nodes"][G]["haplotypes"]) != (4, 3, 4): fail("haplotype counts")
n = full["nodes"]
if not all(n[x]["on_reference"] for x in (A, C, G)): fail("full visits every node")
if [n[x]["ref_offset"] for x in (A, C, G)] != [0, 10, 20] or [n[x]["order"] for x in (A, C, G)] != [0, 1, 2]: fail("offsets/order along full")
if full["reference_strand"] != "source_to_sink" or "ref_start" in n[A]: fail("strand or a coordinate without an interval")
m = dele["nodes"]
if m[C]["on_reference"] or m[C]["between"] != [A, G] or m[C]["between_ref_offsets"] != [9, 10]: fail(f"off-reference node {m[C]}")
if [m[x]["order"] for x in (A, C, G)] != [0, 1, 2]: fail("order with an off-reference node")
if bare["reference_path"] is not None or bare["nodes"][A]["on_reference"] is not None: fail("no reference must mean null positions")
PY
) && ok "the node dictionary carries every walk node, its sequence, haplotype counts and reference positions" \
  || bad "node dictionary: $msg"
# A path named chr:start-end gives 1-based positions on that interval.
{ printf 'H\tVN:Z:1.0\nS\t1\tAAAAAAAAAA\nS\t2\tCCCCCCCCCC\nS\t3\tGGGGGGGGGT\nS\t4\tTTTTTTTTTT\nS\t0\tACGTACGTAC\n'
  printf 'L\t0\t+\t1\t+\t0M\nL\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\n'
  printf 'P\tref#1#chr1:101-150\t0+,1+,2+,3+,4+\t*\nP\talt\t0+,1+,3+,4+\t*\n'; } > "$OUT/pos.gfa"
"$BIN" bubble -i "$OUT/pos.gfa" -r 'ref#1#chr1:101-150' -o "$OUT/pb" --min-variant-bp 0 -q >/dev/null 2>&1
"$BIN" inspect -i "$OUT/pb.sorted.gfa" -c "$OUT/pb.bubbles.csv" -o "$OUT/pi" --genotype-table -r 'ref#1#chr1:101-150' -q >/dev/null 2>&1
pos=$(python3 - "$OUT/pb.sorted.gfa" "$OUT"/pi.bubble_*.nodes.json.gz <<'PY'
import gzip, json, sys
seq = {f[1]: f[2] for f in (l.rstrip("\n").split("\t") for l in open(sys.argv[1])) if f[0] == "S"}
d = json.load(gzip.open(sys.argv[2], "rt"))
print(" ".join(f"{seq[k][:1]}{v.get('ref_start')}-{v.get('ref_end')}" for k, v in d["nodes"].items()))
PY
)
[ "$pos" = "A111-120 C121-130 G131-140" ] && ok "positions are 1-based on the reference path's chr:start-end interval" \
  || bad "interval positions wrong: $pos"

# grep -c, not grep -q: -q exits at the first match and pipefail then reports inspect's SIGPIPE.
miss=$("$BIN" inspect -i "$OUT/gb.sorted.gfa" -c "$OUT/gb.bubbles.csv" --bubble-id 1 -o "$OUT/goth" --genotype-table -r other 2>&1 \
  | grep -c "does not cross bubble 1")
[ "$miss" -ge 1 ] && ok "a reference that does not cross the bubble is named" \
                  || bad "no warning for a reference that misses the bubble"
"$BIN" inspect -i "$OUT/gb.sorted.gfa" -c "$OUT/gb.bubbles.csv" --bubble-id 1 -o "$OUT/gnr" -r full -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "--reference-path without --genotype-table is refused, not ignored" \
               || bad "--reference-path without --genotype-table exited 0"

# Samples: het, the same haplotype twice, two haplotypes sharing an allele, and one missing the bubble.
printf 'sample\thaplotype_1\thaplotype_2\nhet\tfull\tdel\nsame\tfull\tfull\nshared\tfull\trev\nmiss\tfull\tother\n' > "$OUT/cosigt.tsv"
"$BIN" inspect -i "$OUT/gb.sorted.gfa" -c "$OUT/gb.bubbles.csv" --bubble-id 1 -o "$OUT/gs" --samples "$OUT/cosigt.tsv" -q >/dev/null 2>&1
[ "$(gzip -dc "$OUT/gs.bubble_1.genotypes.tsv.gz" | awk -F'\t' 'NR>1{printf "%s=%s%s ", $1, $2, $3}')" = "het=11 same=20 shared=20 miss=NANA " ] \
  && ok "sample rows count allele copies (sum 2) and read NA when a haplotype misses the bubble" \
  || bad "sample rows wrong: $(gzip -dc "$OUT/gs.bubble_1.genotypes.tsv.gz")"
[ "$(gzip -dc "$OUT/gs.bubble_1.genotypes.tsv.gz" | head -1 | cut -f1)" = "sample" ] \
  && ok "--samples implies --genotype-table and labels rows as samples" || bad "sample table header wrong"
printf 'sample\thaplotype_1\thaplotype_2\nbad\tfull\tnope\n' > "$OUT/cosigt_bad.tsv"
"$BIN" inspect -i "$OUT/gb.sorted.gfa" -c "$OUT/gb.bubbles.csv" --bubble-id 1 -o "$OUT/gsb" --samples "$OUT/cosigt_bad.tsv" -q >/dev/null 2>&1
[ "$?" -ne 0 ] && ok "a cosigt haplotype absent from the graph is refused" \
               || bad "an unknown cosigt haplotype exited 0"

echo
if [ "$fails" -eq 0 ]; then echo "inspect_stats: all assertions passed"; exit 0; fi
echo "inspect_stats: $fails assertion(s) FAILED"; exit 1
