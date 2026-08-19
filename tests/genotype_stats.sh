#!/usr/bin/env bash
# genotype_stats.sh - contract assertions for `panvar genotype`'s depth estimation and reporting.
#
# CTest exercised `genotype` nowhere at all: every other module has a stats test and this one, the
# largest in the tree, had none. The depth code now branches over local against shrunk against region
# fallback, three estimators, NA against zero, and the direct against the indexed route, and every one
# of those branches was reached only by hand.
#
# The fixture makes every expected number derivable rather than recorded. Reads are exact copies of the
# sample's two haplotypes, R of each, so an anchor -- a marker every allele carries once -- is seen
# exactly 2R times and lambda is exactly R by construction, not by observation.
#
#   genotype_stats.sh <panvar-binary> <out-dir>
set -uo pipefail

BIN="${1:?usage: genotype_stats.sh <panvar> <outdir>}"
OUT="${2:?}"
mkdir -p "$OUT"
OUT="$OUT/run.$$.$(date +%s)"
rm -rf "$OUT"; mkdir -p "$OUT"
fails=0
ok()  { printf "  ok   %s\n" "$1"; }
bad() { printf "  FAIL %s\n" "$1"; fails=$((fails + 1)); }
# value of <column> for the row whose block_index is <blk>, in a depth audit
dep() { awk -F'\t' -v b="$2" -v w="$3" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i;next}
        $(c["block_index"])==b{print $(c[w]); exit}' "$1"; }

# Deterministic pseudo-random ACGT without python: a seeded linear congruential generator. The
# sequences must be fixed across runs or the syncmer set, and therefore every count, would move.
seq_of() { awk -v n="$1" -v s="$2" 'BEGIN{x=s; b="ACGT";
           for(i=0;i<n;i++){x=(1103515245*x+12345)%2147483648; printf "%s", substr(b,(int(x/65536)%4)+1,1)} }'; }

L=$(seq_of 400 1); B=$(seq_of 400 2)
A1=$(seq_of 120 3); A2=$(seq_of 120 4); A3=$(seq_of 120 5)

# One bubble with three disjoint alleles between two shared backbones. The alleles share no sequence,
# so the bubble's own interior contributes no invariant marker; the shared backbones do, and they sit
# inside every allele's walk, which is where this block's anchors come from. The flanks outside the
# bubble are empty, so they have zero anchors -- the case that produced the original defect.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\n' "$L"; printf 'S\t2\t%s\n' "$A1"; printf 'S\t3\t%s\n' "$A2"
  printf 'S\t4\t%s\n' "$A3"; printf 'S\t5\t%s\n' "$B"
  for a in 2 3 4; do printf 'L\t1\t+\t%s\t+\t0M\nL\t%s\t+\t5\t+\t0M\n' "$a" "$a"; done
  printf 'P\tref\t1+,2+,5+\t*\n'
  for i in 1 2 3 4 5 6; do printf 'P\thapA%d\t1+,2+,5+\t*\n' "$i"; done
  for i in 1 2 3 4 5 6; do printf 'P\thapB%d\t1+,3+,5+\t*\n' "$i"; done
  for i in 1 2 3 4 5 6; do printf 'P\thapC%d\t1+,4+,5+\t*\n' "$i"; done
} > "$OUT/g.gfa"
"$BIN" bubble -i "$OUT/g.gfa" -r ref -o "$OUT/bub" --min-variant-bp 0 -q >/dev/null 2>&1

R=8
: > "$OUT/reads.fa"
for i in $(seq 1 $R); do
  printf '>a%d\n%s%s%s\n>b%d\n%s%s%s\n' "$i" "$L" "$A1" "$B" "$i" "$L" "$A2" "$B"
done >> "$OUT/reads.fa"

# ------------------------------------------------------------------ depth is exact by construction
"$BIN" genotype -i "$OUT/g.gfa" -b "$OUT/bub" -r ref -o "$OUT/exact" \
  -R "$OUT/reads.fa" --depth-model median -q >/dev/null 2>&1
[ "$(dep "$OUT/exact.reads.depth.tsv" 1 lambda_hap)" = "$R" ] \
  && ok "lambda is exactly R=$R when the reads are R copies of each haplotype" \
  || bad "lambda_hap=$(dep "$OUT/exact.reads.depth.tsv" 1 lambda_hap), expected $R"
[ "$(dep "$OUT/exact.reads.depth.tsv" 1 raw_anchor_median)" = "16.000000" ] \
  && ok "raw_anchor_median is exactly 2R=16" \
  || bad "raw_anchor_median=$(dep "$OUT/exact.reads.depth.tsv" 1 raw_anchor_median), expected 16.000000"

# --------------------------------------------- a block with no anchors of its own says so, with NA
# The defect this test exists for: the audit wrote the model's fitted value under the raw column, so a
# block that inherited the region's depth reported a precise-looking local measurement with MAD 0.
for blk in 0 2; do
  [ "$(dep "$OUT/exact.reads.depth.tsv" $blk raw_anchor_median)" = "NA" ] \
    && ok "block $blk has no anchors: raw_anchor_median is NA, not 0" \
    || bad "block $blk raw_anchor_median=$(dep "$OUT/exact.reads.depth.tsv" $blk raw_anchor_median), expected NA"
  [ "$(dep "$OUT/exact.reads.depth.tsv" $blk raw_anchor_mad)" = "NA" ] \
    && ok "block $blk raw_anchor_mad is NA, not a spurious 0" \
    || bad "block $blk raw_anchor_mad=$(dep "$OUT/exact.reads.depth.tsv" $blk raw_anchor_mad), expected NA"
  [ "$(dep "$OUT/exact.reads.depth.tsv" $blk depth_source)" = "REGION_FALLBACK" ] \
    && ok "block $blk depth_source is REGION_FALLBACK" \
    || bad "block $blk depth_source=$(dep "$OUT/exact.reads.depth.tsv" $blk depth_source)"
  [ "$(dep "$OUT/exact.reads.depth.tsv" $blk region_shrink_weight)" = "1" ] \
    && ok "block $blk takes the region's depth whole (weight 1)" \
    || bad "block $blk region_shrink_weight=$(dep "$OUT/exact.reads.depth.tsv" $blk region_shrink_weight)"
  [ "$(dep "$OUT/exact.reads.depth.tsv" $blk n_anchor)" = "0" ] \
    && ok "block $blk really has 0 anchors, so NA is the honest value" \
    || bad "block $blk n_anchor=$(dep "$OUT/exact.reads.depth.tsv" $blk n_anchor), fixture assumption broken"
done

# ------------------------------------------------- a block WITH anchors is shrunk by its own weight
n=$(dep "$OUT/exact.reads.depth.tsv" 1 n_anchor)
[ "$(dep "$OUT/exact.reads.depth.tsv" 1 depth_source)" = "SHRUNK" ] \
  && ok "the anchored block is SHRUNK, not inherited" \
  || bad "anchored block depth_source=$(dep "$OUT/exact.reads.depth.tsv" 1 depth_source)"
want=$(awk -v n="$n" 'BEGIN{printf "%.6g", 200.0/(n+200.0)}')
got=$(awk -v g="$(dep "$OUT/exact.reads.depth.tsv" 1 region_shrink_weight)" 'BEGIN{printf "%.6g", g}')
[ "$want" = "$got" ] \
  && ok "shrinkage coefficient is 200/(n+200) = $want with n=$n anchors" \
  || bad "region_shrink_weight=$got, expected 200/($n+200)=$want"

# ------------------------------------------------- min_anchors flags thin evidence, it does not gate
# Before the cliff was removed, a block below --min-anchors computed no local value at all, so the
# shrinkage pass had nothing to shrink and the block took the region's depth whole.
"$BIN" genotype -i "$OUT/g.gfa" -b "$OUT/bub" -r ref -o "$OUT/ma1" \
  -R "$OUT/reads.fa" --depth-model median --min-anchors 1 -q >/dev/null 2>&1
"$BIN" genotype -i "$OUT/g.gfa" -b "$OUT/bub" -r ref -o "$OUT/ma9" \
  -R "$OUT/reads.fa" --depth-model median --min-anchors 100000 -q >/dev/null 2>&1
[ "$(dep "$OUT/ma1.reads.depth.tsv" 1 fitted_median)" = "$(dep "$OUT/ma9.reads.depth.tsv" 1 fitted_median)" ] \
  && ok "--min-anchors does not change the fitted depth (no cliff)" \
  || bad "fitted_median moved with --min-anchors: $(dep "$OUT/ma1.reads.depth.tsv" 1 fitted_median) vs $(dep "$OUT/ma9.reads.depth.tsv" 1 fitted_median)"
[ "$(dep "$OUT/ma9.reads.depth.tsv" 1 low_anchor)" = "1" ] && [ "$(dep "$OUT/ma1.reads.depth.tsv" 1 low_anchor)" = "0" ] \
  && ok "--min-anchors flips low_anchor without gating the estimate" \
  || bad "low_anchor did not track --min-anchors (got $(dep "$OUT/ma1.reads.depth.tsv" 1 low_anchor) and $(dep "$OUT/ma9.reads.depth.tsv" 1 low_anchor))"
[ "$(dep "$OUT/ma9.reads.depth.tsv" 1 depth_source)" = "SHRUNK" ] \
  && ok "a block below --min-anchors still contributes, still SHRUNK" \
  || bad "below --min-anchors the block became $(dep "$OUT/ma9.reads.depth.tsv" 1 depth_source)"

# ---------------------------------------------------- the estimator reaches the model, not just the log
# Extra copies of the first half of the left backbone, so those anchors read higher than the rest and
# the count distribution is skewed. With a flat distribution median and mean coincide and the
# comparison below would pass while proving nothing.
E=6
cp "$OUT/reads.fa" "$OUT/skew.fa"
for i in $(seq 1 $E); do printf '>x%d\n%s\n' "$i" "${L:0:200}"; done >> "$OUT/skew.fa"
for m in median mean; do
  # Deliberately NOT -q: the region summary is a log line, and quieting it would make the two
  # assertions below pass on an empty file.
  "$BIN" genotype -i "$OUT/g.gfa" -b "$OUT/bub" -r ref -o "$OUT/est_$m" \
    -R "$OUT/skew.fa" --depth-estimator "$m" >"$OUT/est_$m.log" 2>&1
done
lm=$(dep "$OUT/est_median.reads.depth.tsv" 1 lambda_hap)
la=$(dep "$OUT/est_mean.reads.depth.tsv" 1 lambda_hap)
awk -v a="$lm" -v b="$la" 'BEGIN{exit !(b > a)}' \
  && ok "--depth-estimator reaches the DEFAULT joint model (lambda $lm -> $la)" \
  || bad "lambda unchanged by --depth-estimator under the default model: $lm vs $la"
# The raw column must mean one statistic whatever the estimator, or it stops being comparable.
[ "$(dep "$OUT/est_median.reads.depth.tsv" 1 raw_anchor_median)" = "$(dep "$OUT/est_mean.reads.depth.tsv" 1 raw_anchor_median)" ] \
  && ok "raw_anchor_median is the median under every estimator" \
  || bad "raw_anchor_median moved with the estimator"
[ "$(dep "$OUT/est_median.reads.depth.tsv" 1 local_center)" != "$(dep "$OUT/est_mean.reads.depth.tsv" 1 local_center)" ] \
  && ok "local_center follows the estimator, unlike raw_anchor_median" \
  || bad "local_center did not follow the estimator"
grep -q "median .*, mean .*, trimmed " "$OUT/est_median.log" \
  && ok "all three estimator candidates are reported side by side" \
  || bad "the region summary does not report median, mean and trimmed together"
grep -q "selected_anchor_center" "$OUT/est_median.log" \
  && ok "the anchor centre is not reported as lambda" \
  || bad "the region summary still names its anchor centre lambda"

# ------------------------------------------------ the audit describes the depth the emission USED
# Joint refinement replaces every block's fitted depth. Written before that pass, the audit described
# a state the model never used.
[ "$(dep "$OUT/est_mean.reads.depth.tsv" 1 depth_source)" = "JOINT" ] \
  && ok "the audit carries the final joint depth, not the first pass" \
  || bad "depth_source=$(dep "$OUT/est_mean.reads.depth.tsv" 1 depth_source), expected JOINT"

# ---------------------------------------------------------- direct and indexed must not drift apart
# --build-index exits after writing the index but still requires -o; supplied so the mode is reachable.
"$BIN" genotype -i "$OUT/g.gfa" -b "$OUT/bub" -r ref -o "$OUT/idxpfx" \
  --build-index "$OUT/idx.bin" -q >/dev/null 2>&1
if [ -s "$OUT/idx.bin" ]; then
  "$BIN" genotype --index "$OUT/idx.bin" -o "$OUT/viaidx" -R "$OUT/skew.fa" \
    --depth-model median --depth-estimator mean -q >/dev/null 2>&1
  "$BIN" genotype -i "$OUT/g.gfa" -b "$OUT/bub" -r ref -o "$OUT/direct" -R "$OUT/skew.fa" \
    --depth-model median --depth-estimator mean -q >/dev/null 2>&1
  if [ -s "$OUT/viaidx.reads.depth.tsv" ] && [ -s "$OUT/direct.reads.depth.tsv" ]; then
    [ "$(dep "$OUT/viaidx.reads.depth.tsv" 1 lambda_hap)" = "$(dep "$OUT/direct.reads.depth.tsv" 1 lambda_hap)" ] \
      && ok "the indexed route honours --depth-estimator identically to the direct one" \
      || bad "indexed lambda $(dep "$OUT/viaidx.reads.depth.tsv" 1 lambda_hap) != direct $(dep "$OUT/direct.reads.depth.tsv" 1 lambda_hap)"
    # marker_clumps was not serialized, so an indexed run fell back to span/fragment_len for its
    # effective sample size and produced a different GQ from the same reads and the same calls. That is
    # the symptom the index defect actually presented as, so GQ is what has to match. Compare the WHOLE
    # file: a first version of this cut columns 1-12 and claimed to include GQ, which is column 13.
    if diff "$OUT/direct.genotypes.tsv" "$OUT/viaidx.genotypes.tsv" >/dev/null 2>&1; then
      ok "indexed and direct genotype files are byte-identical, GQ included"
    else
      bad "indexed and direct genotype files differ (a panel field is not serialized)"
    fi
    # An index carries the fragment length its clumps were built at. Running at a different one would
    # mix the two: clumps from the index, every other emission term from the caller.
    "$BIN" genotype --index "$OUT/idx.bin" -o "$OUT/fl_bad" -R "$OUT/skew.fa" \
      --depth-model median --fragment-len 1400 -q >"$OUT/fl_bad.log" 2>&1
    grep -qi "index was built at" "$OUT/fl_bad.log" \
      && ok "an index refuses a --fragment-len it was not built at" \
      || bad "an index accepted a conflicting --fragment-len, mixing two clump definitions"
    "$BIN" genotype --index "$OUT/idx.bin" -o "$OUT/fl_ok" -R "$OUT/skew.fa" \
      --depth-model median -q >/dev/null 2>&1
    [ -s "$OUT/fl_ok.genotypes.tsv" ] \
      && ok "an index run without --fragment-len inherits the length it was built at" \
      || bad "an index run without --fragment-len failed"
  else
    bad "indexed or direct run produced no depth audit"
  fi
else
  bad "--build-index wrote no index"
fi

# ------------------------------------------------ --fragment-len reaches the clumping, not just the GQ
# The clump count IS the effective sample size the emission divides by, and it was computed at a
# hard-coded 350 while the flag reached only the emission. A longer fragment must merge clumps.
for FL in 350 1400; do
  "$BIN" genotype -i "$OUT/g.gfa" -b "$OUT/bub" -r ref -o "$OUT/fl$FL" -R "$OUT/reads.fa" \
    --depth-model median --fragment-len "$FL" --dump-markers "$OUT/fl$FL.tsv" -q >/dev/null 2>&1
done
c350=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}{c[$(h["clump"])]=1} END{print length(c)}' "$OUT/fl350.tsv")
c1400=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}{c[$(h["clump"])]=1} END{print length(c)}' "$OUT/fl1400.tsv")
[ "$c350" -gt "$c1400" ] \
  && ok "--fragment-len reaches marker clumping (350 -> $c350 clumps, 1400 -> $c1400)" \
  || bad "clump count did not fall with a longer fragment ($c350 vs $c1400): the flag is not reaching it"

# ---------------------------------------------------- the marker dump separates dosage from efficiency
"$BIN" genotype -i "$OUT/g.gfa" -b "$OUT/bub" -r ref -o "$OUT/dump" -R "$OUT/reads.fa" \
  --depth-model median --truth-haplotypes 'hapA1,hapB1' --dump-markers "$OUT/markers.tsv" -q >/dev/null 2>&1
if [ -s "$OUT/markers.tsv" ]; then
  na=$(awk -F'\t' 'NR>1 && $4=="anchor"' "$OUT/markers.tsv" | wc -l | tr -d ' ')
  [ "$na" = "$n" ] && ok "the dump has one anchor row per anchor ($n)" \
                   || bad "dump has $na anchor rows, expected $n"
  ni=$(awk -F'\t' 'NR>1 && $4=="informative"' "$OUT/markers.tsv" | wc -l | tr -d ' ')
  [ "$ni" -gt 0 ] && ok "the dump also covers informative markers ($ni), not anchors alone" \
                  || bad "no informative marker rows: dosage cannot be separated from efficiency"
  odd=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $4=="anchor" && $(h["count"])!=16{k++} END{print k+0}' "$OUT/markers.tsv")
  [ "$odd" = "0" ] && ok "every anchor count is exactly 2R=16, as the construction requires" \
                   || bad "$odd anchors deviate from 16 on exact-copy reads"
  # An anchor's expected bound was reported as 0 on every row, because it is absent from by_block and
  # the informative-marker expectation map has no entry for it. It is once per traversing haplotype.
  ae=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $4=="anchor"{print $(h["expected"]); exit}' "$OUT/markers.tsv")
  aa=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $4=="anchor"{print $(h["actual"]); exit}' "$OUT/markers.tsv")
  { [ "$ae" != "0" ] && [ "$ae" = "$aa" ]; } \
    && ok "an anchor's expected bound is its own ($ae), not 0, and matches actual" \
    || bad "anchor expected=$ae actual=$aa; expected should be the anchor bound and equal actual here"
  # The point of the whole dump: count divided by this sample's own copy number is the efficiency
  # signal with dosage removed. On exact-copy reads it must be exactly lambda for EVERY marker,
  # anchor and informative alike, or the normalization is wrong.
  bad_cpc=$(awk -F'\t' -v r="$R" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
            $(h["truth_mult"])!="NA" && $(h["truth_mult"])>0 && $(h["count_per_copy"])!=r{k++} END{print k+0}' "$OUT/markers.tsv")
  [ "$bad_cpc" = "0" ] \
    && ok "count_per_copy is exactly lambda=$R for every marker: dosage is fully divided out" \
    || bad "$bad_cpc markers have count_per_copy != $R, so the dosage normalization is wrong"
  # Dosage is summed PER TRUTH HAPLOTYPE from the spelled sequences, so a marker only one of them
  # carries must read 1, not 2 and not NA. hapA1 carries allele A1 and hapB1 carries A2, so both
  # single-haplotype cases exist in this fixture by construction.
  n1=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $(h["truth_mult"])=="1"{k++} END{print k+0}' "$OUT/markers.tsv")
  [ "$n1" -gt 0 ] \
    && ok "markers carried by ONE truth haplotype report dosage 1 ($n1 of them)" \
    || bad "no marker has truth_mult=1, so per-haplotype dosage is not being summed"
  # A marker the sample carries zero times has a defined dosage and an UNDEFINED per-copy rate.
  # Reporting 0 would put a background observation into the efficiency population.
  z=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $(h["truth_mult"])=="0"{k++} END{print k+0}' "$OUT/markers.tsv")
  zbad=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
         $(h["truth_mult"])=="0" && $(h["count_per_copy"])!="NA"{k++} END{print k+0}' "$OUT/markers.tsv")
  { [ "$z" -gt 0 ] && [ "$zbad" = "0" ]; } \
    && ok "zero-dosage markers ($z) report count_per_copy NA, not 0" \
    || bad "$zbad of $z zero-dosage markers report a numeric count_per_copy"
  # Position drives clump membership, and without it a GC trend cannot be told from a position trend.
  npos=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $(h["first_pos"])!="NA"{k++} END{print k+0}' "$OUT/markers.tsv")
  [ "$npos" -gt 0 ] && ok "markers carry a position and a clump ($npos rows)" \
                    || bad "no marker has a position: GC and position stay confounded"
else
  bad "--dump-markers wrote nothing"
fi

echo
if [ "$fails" -eq 0 ]; then echo "genotype_stats: all assertions passed"; else echo "genotype_stats: $fails FAILED"; fi
exit "$fails"
