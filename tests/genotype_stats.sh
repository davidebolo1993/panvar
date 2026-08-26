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
    # Inheritance has to be checked at a NON-DEFAULT length, and at 0. The previous version of this
    # assertion used an index built at the default 350 and only checked that a file appeared, so it
    # could not detect either failure. 0 is the sharp case: it MEANS "disable clumping", and guarding
    # inheritance on `> 0` sent such an index back to 350 -- GQ 33.03 against the direct route's 99.
    for FL in 1400 0; do
      "$BIN" genotype -i "$OUT/g.gfa" -b "$OUT/bub" -r ref -o "$OUT/ixb$FL" \
        --build-index "$OUT/ix$FL.bin" --fragment-len "$FL" -q >/dev/null 2>&1
      "$BIN" genotype -i "$OUT/g.gfa" -b "$OUT/bub" -r ref -o "$OUT/d$FL" -R "$OUT/skew.fa" \
        --depth-model median --fragment-len "$FL" -q >/dev/null 2>&1
      "$BIN" genotype --index "$OUT/ix$FL.bin" -o "$OUT/v$FL" -R "$OUT/skew.fa" \
        --depth-model median -q >/dev/null 2>&1
      if diff "$OUT/d$FL.genotypes.tsv" "$OUT/v$FL.genotypes.tsv" >/dev/null 2>&1; then
        ok "an index built at --fragment-len $FL is inherited exactly (files identical)"
      else
        bad "index built at $FL, run without the flag, differs from direct $FL"
      fi
    done
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

# --------------------------------------- --max-multiplicity 0 means NO cap, in the audit as well
# The bubble audit wrote `if (mult > max_multiplicity) continue`, so at the default 0 it discarded
# every candidate before region filtering ever ran. The region-loss column then read 0 whatever the
# panel contained, and that vacuous zero was quoted as evidence that the filter dropped nothing.
#
# A fixture where nothing is dropped cannot catch this -- the whole defect was a vacuous zero -- so
# the left flank here CONTAINS a copy of allele A1. A1's syncmers are then seen across the panel far
# more often than the bubble's alleles account for, and the filter must remove them.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s%s%s\n' "$(seq_of 200 1)" "$A1" "$(seq_of 200 8)"
  printf 'S\t2\t%s\n' "$A1"; printf 'S\t3\t%s\n' "$A2"; printf 'S\t4\t%s\n' "$A3"
  printf 'S\t5\t%s\n' "$B"
  for a in 2 3 4; do printf 'L\t1\t+\t%s\t+\t0M\nL\t%s\t+\t5\t+\t0M\n' "$a" "$a"; done
  printf 'P\tref\t1+,2+,5+\t*\n'
  for i in 1 2 3 4 5 6; do printf 'P\thapA%d\t1+,2+,5+\t*\n' "$i"; done
  for i in 1 2 3 4 5 6; do printf 'P\thapB%d\t1+,3+,5+\t*\n' "$i"; done
  for i in 1 2 3 4 5 6; do printf 'P\thapC%d\t1+,4+,5+\t*\n' "$i"; done
} > "$OUT/rep.gfa"
"$BIN" bubble -i "$OUT/rep.gfa" -r ref -o "$OUT/repbub" --min-variant-bp 0 -q >/dev/null 2>&1
col() { awk -F'\t' -v w="$2" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
        {s+=$(h[w])} END{print s+0}' "$1"; }
"$BIN" genotype -i "$OUT/rep.gfa" -b "$OUT/repbub" -r ref -o "$OUT/au0" --audit -q >/dev/null 2>&1
"$BIN" genotype -i "$OUT/rep.gfa" -b "$OUT/repbub" -r ref -o "$OUT/auN" --audit \
  --max-multiplicity 1000000 -q >/dev/null 2>&1
"$BIN" genotype -i "$OUT/rep.gfa" -b "$OUT/repbub" -r ref -o "$OUT/au1" --audit \
  --max-multiplicity 1 -q >/dev/null 2>&1
if [ -s "$OUT/au0.audit.alleles.tsv" ] && [ -s "$OUT/auN.audit.alleles.tsv" ]; then
  r0=$(col "$OUT/au0.audit.alleles.tsv" n_retained_nodes)
  [ "$r0" -gt 0 ] \
    && ok "the audit retains candidates at the default cap ($r0 markers), so its filters see input" \
    || bad "the audit retained 0 markers at --max-multiplicity 0: the cap discarded every candidate"
  if diff "$OUT/au0.audit.alleles.tsv" "$OUT/auN.audit.alleles.tsv" >/dev/null 2>&1; then
    ok "--max-multiplicity 0 and an unreachable cap produce identical audits, as 0 = no cap requires"
  else
    bad "the default cap and an unreachable cap disagree, so 0 is being treated as a cap of zero"
  fi
  # And the flag is not inert in the other direction: a cap of 1 must be able to remove something,
  # or the assertion above would pass on a build where the option never reached the audit at all.
  r1=$(col "$OUT/au1.audit.alleles.tsv" n_retained_nodes)
  [ "$r1" -lt "$r0" ] \
    && ok "a cap of 1 retains strictly fewer than no cap ($r1 < $r0): the option reaches the audit" \
    || bad "a cap of 1 retained $r1 against no cap's $r0; the flag is not reaching the audit"
else
  bad "--audit wrote no allele table"
fi

# ------------------------------------------- the ledger accounts for every candidate, kept or not
# Inferring why a block is marker-poor from what survived is how one filter gets blamed for another
# filter's losses. The ledger is taken before the erase, so retained + dropped must be the whole
# candidate set.
#
# TWO fixtures, because one cannot exercise both outcomes. `rep.gfa` produces a block where every
# candidate survives; the first version of this test used only that, asserted a partition over 82
# retained and 0 dropped, and therefore never executed the drop branches at all -- the same vacuity
# it was written to catch. `dup.gfa` repeats a bubble allele in a second bubble, so its candidates
# are seen across the panel more often than the blocks account for and every one is dropped.
A4=$(seq_of 120 6); T=$(seq_of 400 7)
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\n' "$L"; printf 'S\t2\t%s\n' "$A1"; printf 'S\t3\t%s\n' "$A2"
  printf 'S\t5\t%s\n' "$B"; printf 'S\t6\t%s\n' "$A1"; printf 'S\t7\t%s\n' "$A4"
  printf 'S\t8\t%s\n' "$T"
  for a in 2 3; do printf 'L\t1\t+\t%s\t+\t0M\nL\t%s\t+\t5\t+\t0M\n' "$a" "$a"; done
  for a in 6 7; do printf 'L\t5\t+\t%s\t+\t0M\nL\t%s\t+\t8\t+\t0M\n' "$a" "$a"; done
  printf 'P\tref\t1+,2+,5+,6+,8+\t*\n'
  for i in 1 2 3 4; do printf 'P\thapW%d\t1+,2+,5+,6+,8+\t*\n' "$i"; done
  for i in 1 2 3 4; do printf 'P\thapX%d\t1+,3+,5+,7+,8+\t*\n' "$i"; done
  for i in 1 2 3 4; do printf 'P\thapY%d\t1+,2+,5+,7+,8+\t*\n' "$i"; done
} > "$OUT/dup.gfa"
"$BIN" bubble -i "$OUT/dup.gfa" -r ref -o "$OUT/dupbub" --min-variant-bp 0 -q >/dev/null 2>&1

# Which chain index carries a block's candidates is the block builder's business, not this test's.
find_ledger() {   # <gfa> <bubprefix> <tag> -> echoes the path of the first non-empty ledger
  for b in 0 1 2 3 4; do
    "$BIN" genotype -i "$1" -b "$2" -r ref -o "$OUT/$3$b" --ledger-block "$b" -q >/dev/null 2>&1
    if [ -s "$OUT/$3$b.block$b.ledger.tsv" ]; then echo "$OUT/$3$b.block$b.ledger.tsv"; return; fi
  done
}
fates() { awk -F'\t' -v w="$2" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
          $(h["fate"])~w{k++} END{print k+0}' "$1"; }

LED=$(find_ledger "$OUT/rep.gfa" "$OUT/repbub" led)
DROPLED=$(find_ledger "$OUT/dup.gfa" "$OUT/dupbub" dled)
if [ -n "$LED" ] && [ -n "$DROPLED" ]; then
  for f in "$LED" "$DROPLED"; do
    tot=$(awk 'NR>1' "$f" | wc -l | tr -d ' ')
    kept=$(fates "$f" '^retained$'); drop=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
                                            $(h["fate"])!="retained"{k++} END{print k+0}' "$f")
    [ "$tot" = "$((kept + drop))" ] \
      && ok "ledger partitions its rows: $kept retained + $drop dropped = $tot ($(basename "$f"))" \
      || bad "ledger fates do not partition its $tot rows ($(basename "$f"))"
    bf=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
         $(h["fate"])!="retained" && $(h["fate"])!="multi_block" &&
         $(h["fate"])!="over_expected" && $(h["fate"])!="both"{k++} END{print k+0}' "$f")
    [ "$bf" = "0" ] || bad "$bf rows carry an unrecognized fate ($(basename "$f"))"
    # A row's fate must follow from the numbers printed beside it, or the ledger cannot be audited.
    incons=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
             { m=($(h["vary_blocks"])>1); o=($(h["actual"])>$(h["expected"]));
               want=(m&&o)?"both":(m?"multi_block":(o?"over_expected":"retained"));
               if ($(h["fate"])!=want) k++} END{print k+0}' "$f")
    [ "$incons" = "0" ] \
      && ok "each fate follows from its own vary_blocks and actual/expected ($(basename "$f"))" \
      || bad "$incons rows have a fate their columns do not imply ($(basename "$f"))"
  done

  # COVERAGE GAP, stated rather than papered over: these fixtures exercise `retained` and
  # `over_expected` only. `multi_block` needs one variable sequence living in two SEPARATE chain
  # blocks, and the block builder collapses these hand-written bubbles into a single block, so the
  # multi_block and both branches are not reached here. They do occur on real panels (gstm1 block 3
  # is 7 multi-block of 11; cyp2d6 block 5 is 16 both of 36). Anyone extending this file should aim
  # at that gap first.
  # The point of the second fixture: the drop branches must actually execute somewhere.
  nd=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $(h["fate"])!="retained"{k++} END{print k+0}' "$DROPLED")
  [ "$nd" -gt 0 ] \
    && ok "a dropped-marker fixture records $nd non-retained rows, so the drop branches run" \
    || bad "no ledger row records a drop: the fate logic is never exercised"
  nk=$(fates "$LED" '^retained$')
  [ "$nk" -gt 0 ] && ok "and the other fixture records $nk retained rows, so both branches run" \
                  || bad "no ledger row records a retained marker"

  # Adjacencies pass through the same two rules, and a node-only ledger cannot say whether longer
  # context escapes the filters that remove single syncmers.
  ne=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $(h["unit"])=="edge"{k++} END{print k+0}' "$LED")
  nn=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $(h["unit"])=="node"{k++} END{print k+0}' "$LED")
  { [ "$ne" -gt 0 ] && [ "$nn" -gt 0 ]; } \
    && ok "the ledger covers both marker units ($nn node, $ne edge)" \
    || bad "ledger has $nn node and $ne edge rows; it must cover both"

  # The ledger's fates must match what actually reached the panel the model scores against.
  # "No dropped marker leaked into the panel" is NOT enough: at a block where everything is dropped
  # the panel is empty and the check passes on nothing. Assert the EQUALITY instead -- the retained
  # NODE rows must be exactly the markers the dump shows surviving -- which is non-vacuous whichever
  # way the block went. (`conf.tsv` carries nodes only, so edges are excluded from the count.)
  for pair in "$LED:rep:repbub" "$DROPLED:dup:dupbub"; do
    f="${pair%%:*}"; rest="${pair#*:}"; gfa="${rest%%:*}"; bub="${rest#*:}"
    BLK=$(basename "$f" | sed -E 's/.*\.block([0-9]+)\.ledger\.tsv/\1/')
    "$BIN" genotype -i "$OUT/$gfa.gfa" -b "$OUT/$bub" -r ref -o "$OUT/c$gfa" -R "$OUT/reads.fa" \
      --dump-block "$BLK" -q >/dev/null 2>&1
    C="$OUT/c$gfa.block$BLK.conf.tsv"
    if [ -s "$C" ]; then
      # Compare the SETS, not their sizes. Equal counts with different members would pass a count
      # check while meaning the ledger describes a different panel than the one being scored.
      awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
           $(h["fate"])=="retained" && $(h["unit"])=="node"{print $(h["slot"])}' "$f" \
        | sort -u > "$OUT/$gfa.want.slots"
      awk -F'\t' 'NR>1{print $2}' "$C" | sort -u > "$OUT/$gfa.got.slots"
      if diff -q "$OUT/$gfa.want.slots" "$OUT/$gfa.got.slots" >/dev/null 2>&1; then
        ok "the ledger's retained node SET is exactly the scored panel's ($(wc -l < "$OUT/$gfa.want.slots" | tr -d ' '), $gfa)"
      else
        bad "ledger's retained node set differs from the scored panel's ($gfa): $(diff "$OUT/$gfa.want.slots" "$OUT/$gfa.got.slots" | head -3 | tr '\n' ' ')"
      fi
    else
      bad "--dump-block wrote no conf table for $gfa"
    fi
  done

else
  bad "--ledger-block wrote nothing for one of the two fixtures"
fi

# ------------------------------------------------------ --probe-pair, and the sentinel it can forge
# The emission-rank diagnostics report on the truth pair, but at an unrepresentable block the truth
# is not in the panel and the pair worth asking about is the certified-optimal one. --probe-pair
# points them at any pair. Its danger is that rank -2 means "pruned before scoring", so an
# out-of-range allele index would read as a pruning finding -- the exact conclusion the flag exists
# to establish. It must be rejected, not silently reported.
"$BIN" genotype -i "$OUT/g.gfa" -b "$OUT/bub" -r ref -o "$OUT/pt" -R "$OUT/reads.fa" \
  --truth-haplotypes 'hapA1,hapB1' -q >/dev/null 2>&1
"$BIN" genotype -i "$OUT/g.gfa" -b "$OUT/bub" -r ref -o "$OUT/pp" -R "$OUT/reads.fa" \
  --truth-haplotypes 'hapA1,hapB1' --probe-pair 1:0,2 -q >/dev/null 2>&1
tgt() { awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}{print $(h["rank_target"]); exit}' "$1"; }
{ [ "$(tgt "$OUT/pt.genotypes.tsv")" = "truth" ] && [ "$(tgt "$OUT/pp.genotypes.tsv")" = "probe" ]; } \
  && ok "rank_target distinguishes a probe run from a scored one" \
  || bad "rank_target does not mark probe rows: a diagnostic run could be read as a scored one"
# A probe changes only the reported ranks, never the call itself.
c1=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}{print $(h["allele1"])","$(h["allele2"])}' "$OUT/pt.genotypes.tsv" | tr '\n' ' ')
c2=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}{print $(h["allele1"])","$(h["allele2"])}' "$OUT/pp.genotypes.tsv" | tr '\n' ' ')
[ "$c1" = "$c2" ] && ok "a probe run reports the same calls: it observes, it does not steer" \
                  || bad "the calls changed under --probe-pair; the probe is affecting inference"
probe_refuses() {
  "$BIN" genotype -i "$OUT/g.gfa" -b "$OUT/bub" -r ref -o "$OUT/px" -R "$OUT/reads.fa" \
    --probe-pair "$1" -q >"$OUT/px.log" 2>&1
  if [ $? -ne 0 ] && grep -qi "$2" "$OUT/px.log"; then ok "--probe-pair rejects $3"
  else bad "--probe-pair accepted $3 (exit $?), which rank -2 would report as pruning"; fi
}
probe_refuses "1:0,99"    "allele out of range" "an allele index past the end"
probe_refuses "99:0,1"    "out of range"        "a block index past the end"
probe_refuses "1:0,1,2"   "exactly BLOCK"       "a third field"
probe_refuses "1:x,1"     "whole number"        "a non-numeric allele"
probe_refuses "1:0,1junk" "whole number"        "trailing junk"

# ----------------------------------- the certified oracle's edit columns, and which row owns them
# excess_total_edits is (called - best) EDITS, and only the edit_distance row's `best` minimises
# edits. On this fixture the length_error row's "best" pair sits 69 edits from truth while the
# edit-optimal one sits at 0 -- so a number there would read as a certified optimum and be wrong by
# the whole distance. The other rows must say NA.
"$BIN" genotype -i "$OUT/g.gfa" -b "$OUT/bub" -r ref -o "$OUT/orc" -R "$OUT/reads.fa" \
  --truth-haplotypes 'hapA1,hapB1' --certified-oracle -q >/dev/null 2>&1
ORC="$OUT/orc.oracle.tsv"
if [ -s "$ORC" ]; then
  ocol() { awk -F'\t' -v c="$2" -v w="$3" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
           $(h["criterion"])==c{print $(h[w]); exit}' "$1"; }
  # One edit_distance row per scored block, so an aggregate cannot triple-count.
  nb=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}{b[$(h["block_index"])]=1} END{print length(b)}' "$ORC")
  ne=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next} $(h["criterion"])=="edit_distance"{k++} END{print k+0}' "$ORC")
  [ "$nb" = "$ne" ] && ok "the oracle writes exactly one edit_distance row per scored block ($nb)" \
                    || bad "$nb blocks but $ne edit_distance rows: an aggregate would miscount"
  # The edit columns belong to that row alone.
  na=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
       $(h["criterion"])!="edit_distance" && ($(h["called_total_edits"])!="NA" || $(h["excess_total_edits"])!="NA"){k++}
       END{print k+0}' "$ORC")
  [ "$na" = "0" ] \
    && ok "only the edit_distance row carries edit columns; the others say NA" \
    || bad "$na non-edit_distance rows carry an edit number that is not a certified optimum"
  # Excess is a minimum over all pairs, so it can never be negative.
  neg=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
        $(h["excess_total_edits"])!="NA" && $(h["excess_total_edits"])+0<0{k++} END{print k+0}' "$ORC")
  [ "$neg" = "0" ] && ok "no certified excess is negative (it is a minimum over every pair)" \
                   || bad "$neg rows report a negative certified excess"
  # This fixture's truth pair IS in the panel and IS what was called, which pins all three at 0.
  # Every block here has a representable truth that was called exactly, so all three must be 0
  # everywhere -- checked across blocks, not on whichever row happens to come first.
  nz=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
       $(h["criterion"])=="edit_distance" &&
       ($(h["best_total_edits"])+0!=0 || $(h["called_total_edits"])+0!=0 ||
        $(h["excess_total_edits"])+0!=0){k++} END{print k+0}' "$ORC")
  [ "$nz" = "0" ] \
    && ok "a representable truth called exactly gives best=called=excess=0 edits, at every block" \
    || bad "$nz blocks have a nonzero edit count though truth is in the panel and was called"
  # And the criteria really do disagree SOMEWHERE, so the NA rule is not guarding a non-issue.
  # Asked per block rather than off the first row: which block comes first is the chain builder's
  # business, and an earlier version of this test silently started reading a block where the two
  # criteria happen to agree.
  dis=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
        { b=$(h["block_index"]); c=$(h["criterion"]); e=$(h["best_total_edits"])
          if (c=="edit_distance") ed[b]=e; if (c=="length_error") le[b]=e }
        END{ for (b in ed) if (le[b] != "" && le[b]+0 > ed[b]+0) { print b" "le[b]" "ed[b]; exit } }' "$ORC")
  [ -n "$dis" ] \
    && ok "a length-optimal pair is farther in edits than the edit-optimal one (block $dis): criteria disagree" \
    || bad "no block where length_error and edit_distance disagree: the fixture cannot test the NA rule"
  # An EMPTY truth sequence is a genotype (the deletion/bypass allele), not a missing value. The
  # evaluation used the empty string as a missing-value sentinel and silently skipped such blocks,
  # which is precisely the case a deletion caller has to get right. This fixture's flanks spell
  # empty, so they must appear.
  nblk=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}{b[$(h["block_index"])]=1} END{print length(b)}' "$ORC")
  nchain=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}{b[$(h["block_index"])]=1} END{print length(b)}' "$OUT/orc.genotypes.tsv")
  [ "$nblk" = "$nchain" ] \
    && ok "the oracle scores every block in the chain ($nblk), empty-truth ones included" \
    || bad "the oracle scored $nblk of $nchain blocks: an empty truth is being read as no truth"
else
  bad "--certified-oracle wrote no table"
fi

# ------------------------------------- the linkage constraint: what the chain may and may not do
# Measured motivation: linkage moved off a UNIQUE block-local optimum 93 times over the cohort, 20
# rescues against 73 overrides, separating cleanly by how far it moved. The constraint excludes
# states losing more than tau BEFORE forward-backward, so posterior and GQ are recomputed under it.
lk() {  # <tag> [flags...] -> md5 of the genotypes table
  "$BIN" genotype -i "$OUT/g.gfa" -b "$OUT/bub" -r ref -o "$OUT/$1" -R "$OUT/reads.fa" \
    --truth-haplotypes 'hapA1,hapB1' "${@:2}" -q >/dev/null 2>&1
  md5 -q "$OUT/$1.genotypes.tsv" 2>/dev/null || md5sum "$OUT/$1.genotypes.tsv" | cut -d" " -f1
}
LK_ABSENT=$(lk lk_a)
LK_INF=$(lk lk_i --max-linkage-emission-loss inf)
[ "$LK_ABSENT" = "$LK_INF" ] \
  && ok "the constraint is absent by default: no flag and explicit inf are byte-identical" \
  || bad "default output moved: --max-linkage-emission-loss inf differs from omitting it"

# The defining property, asserted exactly rather than by eyeball: with tau, no reported call may sit
# more than tau below its block's emission optimum. called_delta measures precisely that.
worst_delta() { awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
                { d=$(h["called_delta"]); if (d != "" && d < m) m = d } END{printf "%.6f", m+0}' "$1"; }
for TAU in 0 0.25 2; do
  lk "lk_t$TAU" --max-linkage-emission-loss "$TAU" >/dev/null
  wd=$(worst_delta "$OUT/lk_t$TAU.genotypes.tsv")
  # awk has no abs on a string compare; test -wd <= tau via awk to keep float semantics.
  okk=$(awk -v w="$wd" -v t="$TAU" 'BEGIN{print (-w <= t + 1e-6) ? 1 : 0}')
  [ "$okk" = "1" ] \
    && ok "at tau=$TAU no call falls further than tau below its block optimum (worst $wd)" \
    || bad "at tau=$TAU a call sits $wd below its block optimum, outside the constraint"
done

# tau=0 is NOT "off": it admits every state TIED with the optimum, so linkage can still resolve a
# tie. At a block with no markers every pair ties, so tau=0 must change nothing there.
"$BIN" genotype -i "$OUT/dup.gfa" -b "$OUT/dupbub" -r ref -o "$OUT/lk_tie0" -R "$OUT/reads.fa" \
  --truth-haplotypes 'hapW1,hapX1' -q >/dev/null 2>&1
"$BIN" genotype -i "$OUT/dup.gfa" -b "$OUT/dupbub" -r ref -o "$OUT/lk_tie1" -R "$OUT/reads.fa" \
  --truth-haplotypes 'hapW1,hapX1' --max-linkage-emission-loss 0 -q >/dev/null 2>&1
if diff -q "$OUT/lk_tie0.genotypes.tsv" "$OUT/lk_tie1.genotypes.tsv" >/dev/null 2>&1; then
  ok "tau=0 leaves a fully-tied block untouched: tied states stay admissible, linkage still decides"
else
  bad "tau=0 changed a block whose states all tie, so it is excluding the optimum itself"
fi

# The constraint must reach the posterior, not just the argmax: GQ is recomputed under it.
gqs() { awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}{s+=$(h["gq"])} END{printf "%.4f", s}' "$1"; }
g_un=$(gqs "$OUT/lk_a.genotypes.tsv"); g_c=$(gqs "$OUT/lk_t0.genotypes.tsv")
if [ "$LK_ABSENT" = "$(lk lk_t0b --max-linkage-emission-loss 0)" ]; then
  ok "tau=0 reproduces the unconstrained run on this fixture (its optima are already unique)"
else
  [ "$g_un" != "$g_c" ] \
    && ok "the constraint reaches the posterior: total GQ moves ($g_un -> $g_c), not just the call" \
    || bad "the call changed under tau=0 but total GQ did not: GQ is inherited, not recomputed"
fi

# And it must not depend on how the work is divided.
t1=$(lk lk_j1 --max-linkage-emission-loss 0.25 --threads 1)
t8=$(lk lk_j8 --max-linkage-emission-loss 0.25 --threads 8)
[ "$t1" = "$t8" ] && ok "the constrained run is thread-independent (1 vs 8)" \
                  || bad "constrained output differs between 1 and 8 threads"

# --------------------------------------- the emission's tie count, without which rank 1 is a lie
# truth_rank counts only STRICTLY better pairs, so a block whose markers separate nothing reports
# rank 1 for every pair -- and read as "the emission ranked truth first" that turns an evidence
# failure into an apparent linkage failure. On the cohort it would have said 94% of leave-one-out
# errors were the chain overriding a correct emission. truth_ties is what tells them apart.
ties_of() { awk -F'\t' -v b="$2" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
            $(h["block_index"])==b{print $(h["truth_ties"]); exit}' "$1"; }
col_of()  { awk -F'\t' -v b="$2" -v w="$3" 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
            $(h["block_index"])==b{print $(h[w]); exit}' "$1"; }

# (a) a block with NO retained markers: every pair scores alike, so all A(A+1)/2 of them tie.
"$BIN" genotype -i "$OUT/dup.gfa" -b "$OUT/dupbub" -r ref -o "$OUT/tie0" -R "$OUT/reads.fa" \
  --truth-haplotypes 'hapW1,hapX1' -q >/dev/null 2>&1
if [ -s "$OUT/tie0.genotypes.tsv" ]; then
  found=0
  for b in 0 1 2 3 4; do
    mk=$(col_of "$OUT/tie0.genotypes.tsv" "$b" n_markers)
    na=$(col_of "$OUT/tie0.genotypes.tsv" "$b" n_alleles)
    ti=$(ties_of "$OUT/tie0.genotypes.tsv" "$b")
    [ -n "$mk" ] && [ "$mk" = "0" ] && [ -n "$na" ] && [ "$na" -gt 1 ] || continue
    found=1
    want=$(( na * (na + 1) / 2 ))
    [ "$ti" = "$want" ] \
      && ok "a block with no markers ties every pair ($ti = ${na}x$((na+1))/2), so rank 1 means nothing" \
      || bad "block $b has 0 markers and $na alleles: expected $want ties, got $ti"
    break
  done
  [ "$found" = "1" ] || bad "no zero-marker multi-allele block in the fixture to test ties against"
else
  bad "the tie fixture produced no genotypes"
fi

# (b) a block whose markers DO separate: exactly one pair at the top. Without this the assertion
# above would pass on a build that reported A(A+1)/2 unconditionally.
"$BIN" genotype -i "$OUT/g.gfa" -b "$OUT/bub" -r ref -o "$OUT/tie1" -R "$OUT/reads.fa" \
  --truth-haplotypes 'hapA1,hapB1' --depth-model median -q >/dev/null 2>&1
if [ -s "$OUT/tie1.genotypes.tsv" ]; then
  uniq_seen=0
  for b in 0 1 2 3; do
    mk=$(col_of "$OUT/tie1.genotypes.tsv" "$b" n_markers)
    ti=$(ties_of "$OUT/tie1.genotypes.tsv" "$b")
    [ -n "$mk" ] && [ "$mk" -gt 0 ] && [ "$ti" = "1" ] && uniq_seen=1
  done
  [ "$uniq_seen" = "1" ] \
    && ok "a block whose markers separate reports exactly one pair at the top" \
    || bad "no well-supplied block reports a unique optimum: the tie count is not discriminating"
  # Bounds hold everywhere: at least one pair is at the top, never more than there are pairs.
  bad_b=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i;next}
          { t=$(h["truth_ties"]); a=$(h["n_alleles"]); if (t=="-1") next;
            if (t < 1 || t > a*(a+1)/2) k++ } END{print k+0}' "$OUT/tie1.genotypes.tsv")
  [ "$bad_b" = "0" ] && ok "every tie count lies in [1, A(A+1)/2]" \
                     || bad "$bad_b blocks report a tie count outside [1, A(A+1)/2]"
else
  bad "the unique-optimum fixture produced no genotypes"
fi

# ------------------------------------------------ an EMPTY allele is an allele, and must be dumped
# The bypass a deletion takes has no sequence. --dump-block skipped such alleles, so any consumer
# counting FASTA records silently dropped the one allele a deletion is about -- which is how an
# identifiability experiment came to score a block as 13 alleles when the panel has 14.
#
# The invariant that broke is not "an empty record exists" but "the dump and the panel agree on how
# many alleles the block has". That holds for every block, needs no special fixture, and fails
# against the old code on any block containing a bypass.
"$BIN" genotype -i "$OUT/rep.gfa" -b "$OUT/repbub" -r ref -o "$OUT/hm" \
  --dump-haplotype-alleles "$OUT/hm.tsv" -q >/dev/null 2>&1
if [ -s "$OUT/hm.tsv" ]; then
  # How many blocks the PANEL says exist. Counting only blocks whose dump is non-empty is what let
  # the old bug hide: it emptied two dumps entirely, and a loop that skips empty files then checks
  # fewer blocks and still reports success.
  nblocks=$(awk -F'\t' 'NR>1{b[$2]=1} END{print length(b)}' "$OUT/hm.tsv")
  agreed=0; checked=0
  for b in 0 1 2 3; do
    "$BIN" genotype -i "$OUT/rep.gfa" -b "$OUT/repbub" -r ref -o "$OUT/fa$b" -R "$OUT/reads.fa" \
      --dump-block "$b" -q >/dev/null 2>&1
    F="$OUT/fa$b.block$b.fa"
    npanel=$(awk -F'\t' -v bb="$b" 'NR>1 && $2==bb{a[$4]=1} END{print length(a)}' "$OUT/hm.tsv")
    [ "$npanel" -gt 0 ] || continue
    nrec=$([ -s "$F" ] && grep -c '^>' "$F" || echo 0)
    checked=$((checked + 1))
    [ "$nrec" = "$npanel" ] && agreed=$((agreed + 1)) \
      || bad "block $b: --dump-block wrote $nrec alleles, the panel matrix has $npanel"
    # Every record must also carry a sequence line, empty or not, or it cannot be re-read.
    nhdr="$nrec"; nline=$([ -s "$F" ] && awk 'END{print NR}' "$F" || echo 0)
    [ "$nline" = "$((nhdr * 2))" ] \
      || bad "block $b: $nhdr headers but $nline lines; a record is missing its sequence line"
  done
  { [ "$checked" = "$nblocks" ] && [ "$agreed" = "$checked" ]; } \
    && ok "--dump-block and the panel matrix agree on the allele count (all $nblocks blocks)" \
    || bad "checked $checked of the panel's $nblocks blocks, $agreed agreed: a dump is short or missing"
else
  bad "--dump-haplotype-alleles wrote nothing"
fi

# ----------------------------------------------------- the ledger refuses questions it cannot answer
# A refusal must FAIL, not merely print. A message on stderr with exit 0 is a success as far as any
# caller is concerned, and every one of these contracts exists to stop a silent no-op.
refuses() {   # <logfile> <expected text> <description> -- run the command before calling
  if [ "$1" -ne 0 ] && grep -qi "$3" "$2"; then ok "$4"
  elif [ "$1" -eq 0 ]; then bad "$4 -- but it exited 0, so a caller sees success"
  else bad "$4 -- exited nonzero with an unexpected message"
  fi
}
"$BIN" genotype -i "$OUT/rep.gfa" -b "$OUT/repbub" -r ref -o "$OUT/lb" \
  --ledger-block 99 -q >"$OUT/lb_range.log" 2>&1
refuses $? "$OUT/lb_range.log" "out of range" \
  "an out-of-range --ledger-block is rejected, not silently ignored"
"$BIN" genotype -i "$OUT/rep.gfa" -b "$OUT/repbub" -r ref -o "$OUT/lb2" \
  --ledger-block 1 --no-region-unique -q >"$OUT/lb_nru.log" 2>&1
refuses $? "$OUT/lb_nru.log" "cannot be combined" \
  "--ledger-block refuses --no-region-unique, whose fates would all read retained"
# The ledger is a property of the panel, so it must not require reads to produce one.
"$BIN" genotype -i "$OUT/rep.gfa" -b "$OUT/repbub" -r ref -o "$OUT/lb3" \
  --ledger-block 1 -q >/dev/null 2>&1
ls "$OUT"/lb3.block*.ledger.tsv >/dev/null 2>&1 \
  && ok "--ledger-block needs no reads: the ledger describes the panel, not a sample" \
  || bad "--ledger-block produced nothing without reads, though it describes only the panel"

echo
if [ "$fails" -eq 0 ]; then echo "genotype_stats: all assertions passed"; else echo "genotype_stats: $fails FAILED"; fi
exit "$fails"
