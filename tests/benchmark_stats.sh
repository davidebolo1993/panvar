#!/usr/bin/env bash
# benchmark_stats.sh - contract assertions for `panvar benchmark`.
#
# real_smoke.sh checks that benchmark RUNS. It asserts nothing about what the numbers mean, and the
# two headline claims -- "this much variation was missed" and "this much of the gap was closed" --
# were both wrong in ways no end-to-end run could show:
#
#   * the residual split classified by contiguous ALIGNMENT RUN length. A clean 60 bp deletion of
#     non-repetitive sequence comes back from edlib as fourteen runs of 1-10 bases, because the
#     co-optimal edit path distributes the gap over chance matches. Every real locus therefore
#     reported over_threshold_bp = 0: the bucket meant to hold missed callable-size events was
#     structurally empty, and the carrier truth flag built on it was empty with it.
#   * gap_closed divided by (baseline - graph). When a bubble is missed entirely those are equal, so
#     the ratio is undefined -- and it was reported as 1.0, i.e. a total miss scored 100%.
#
# Both fixtures below are small enough that every expected number is derivable by hand.
#
#   benchmark_stats.sh <panvar-binary> <out-dir>
set -uo pipefail

BIN="${1:?usage: benchmark_stats.sh <panvar> <outdir>}"
OUT="${2:?}"
mkdir -p "$OUT"
OUT="$OUT/run.$$.$(date +%s)"
rm -rf "$OUT"; mkdir -p "$OUT"
fails=0
ok()  { printf "  ok   %s\n" "$1"; }
bad() { printf "  FAIL %s\n" "$1"; fails=$((fails + 1)); }
check() { [ "$2" = "$3" ] && ok "$1 = $3" || bad "$1: expected $3, got ${2:-<missing>}"; }

# Deterministic pseudo-random sequence (LCG). NOT the repeating k-mer the other harnesses use: the
# point of several assertions here is that a clean event in ordinary, non-repetitive sequence is
# classified correctly, and periodic sequence would make that untestable.
gblk() { local n=$1 s=${2:-1} out="" i=0 b
  while [ $i -lt "$n" ]; do s=$(( (s * 1103515245 + 12345) % 2147483648 ))
    b=$(( (s / 65536) % 4 ))
    case $b in 0) out="${out}A";; 1) out="${out}C";; 2) out="${out}G";; *) out="${out}T";; esac
    i=$((i+1)); done; printf '%s' "$out"; }
sumcol() { awk -F'\t' -v s="$1" -v k="$3" '$1==s && $3==k {print $4; exit}' "$2"; }
pctcol() { awk -F'\t' -v s="$1" -v k="$3" '$1==s && $3==k {print $5; exit}' "$2"; }

# ---------------------------------------------------------------------------------------------
# Fixture A: three independent bubbles on one reference chain, so one call run can present a called
# event, an uncalled above-threshold event and a below-threshold event at once.
#
#   f1 100    1..100
#   dA  60  101..160   deleted by hapdel and hapall            (60 bp, callable)
#   f2 100  161..260
#   dB  60  261..320   deleted by hapmiss and hapall           (60 bp, callable -- but see below)
#   f3 100  321..420
#   dC  20  421..440   deleted by hapsmall and hapall          (20 bp, BELOW a 50 bp threshold)
#   f4 100  441..540
# ---------------------------------------------------------------------------------------------
F1=$(gblk 100 11); DA=$(gblk 60 22); F2=$(gblk 100 33); DB=$(gblk 60 44)
F3=$(gblk 100 55); DC=$(gblk 20 66); F4=$(gblk 100 77)
REF='ref#0#chr1:1-540'
{ printf 'H\tVN:Z:1.0\n'
  printf "S\tf1\t%s\nS\tdA\t%s\nS\tf2\t%s\nS\tdB\t%s\n" "$F1" "$DA" "$F2" "$DB"
  printf "S\tf3\t%s\nS\tdC\t%s\nS\tf4\t%s\n" "$F3" "$DC" "$F4"
  for e in "f1 dA" "dA f2" "f1 f2" "f2 dB" "dB f3" "f2 f3" "f3 dC" "dC f4" "f3 f4"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"; done
  printf 'P\t%s\tf1+,dA+,f2+,dB+,f3+,dC+,f4+\t*\n' "$REF"
  printf 'P\thapref#1#chr1\tf1+,dA+,f2+,dB+,f3+,dC+,f4+\t*\n'
  printf 'P\thapdel#1#chr1\tf1+,f2+,dB+,f3+,dC+,f4+\t*\n'
  printf 'P\thapmiss#1#chr1\tf1+,dA+,f2+,f3+,dC+,f4+\t*\n'
  printf 'P\thapsmall#1#chr1\tf1+,dA+,f2+,dB+,f3+,f4+\t*\n'; } > "$OUT/a.gfa"

"$BIN" bubble -i "$OUT/a.gfa" -r "$REF" -o "$OUT/ab" --min-variant-bp 1 -q >/dev/null 2>&1
nb=$(tail -n +2 "$OUT/ab.bubbles.csv" 2>/dev/null | wc -l | tr -d ' ')
check "fixture A presents three bubbles" "$nb" "3"

"$BIN" call -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" -r "$REF" -o "$OUT/ac" -q >/dev/null 2>&1
"$BIN" benchmark -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" -r "$REF" \
   --variant-nodes "$OUT/ac.variant_nodes.tsv" --vcf "$OUT/ac.region.vcf" -o "$OUT/A" -q >/dev/null 2>&1
AS="$OUT/A.qv_summary.tsv"
AQ="$OUT/A.qv.tsv"
AE="$OUT/A.truth_events.tsv"
[ -s "$AS" ] || { bad "benchmark wrote no summary for fixture A"; printf "%d assertion(s) failed\n" "$((fails+1))"; exit 1; }

# --- The truth ledger sizes events from the WALKS, so a clean deletion is one event of its own size.
ev_hapdel=$(awk -F'\t' '$1 ~ /^hapdel/ {print $5"/"$6}' "$AE" | tr '\n' ' ' | sed 's/ $//')
check "a clean 60 bp deletion is ONE 60 bp truth event" "$ev_hapdel" "60/called"
ev_small=$(awk -F'\t' '$1 ~ /^hapsmall/ {print $5"/"$6}' "$AE" | tr '\n' ' ' | sed 's/ $//')
check "a 20 bp deletion is one below-threshold event" "$ev_small" "20/below_threshold"

tm=$(sumcol truth_event "$AS" missed)
tc=$(sumcol truth_event "$AS" called)
tb=$(sumcol truth_event "$AS" below_threshold)
check "fixture A: no above-threshold event is missed" "$tm" "0"
check "fixture A: two above-threshold events are called" "$tc" "2"
check "fixture A: one below-threshold event" "$tb" "1"

# --- Carrier truth comes from the walks: every deletion carrier is a true carrier of its own bubble.
tcol=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="gt_true_carrier") print i}' "$AQ")
tcar=$(awk -F'\t' -v c="$tcol" 'NR>1 && $c==1' "$AQ" | wc -l | tr -d ' ')
check "true carriers = the two callable-size deletion haplotypes only" "$tcar" "2"

# --- The four reconstructions are a chain of upper bounds. `carrier` substitutes a SUBSET of the
# blocks `called` does (the same true blocks, minus those this haplotype is not genotyped as
# carrying), so it can never beat it, and `graph` can never be beaten by either.
vg=$(pctcol variation_recovered "$AS" graph)
vc=$(pctcol variation_recovered "$AS" called)
vw=$(pctcol variation_recovered "$AS" carrier_walk)
vt=$(pctcol variation_recovered "$AS" genotype)
ordered=$(awk -v g="$vg" -v c="$vc" -v w="$vw" -v t="$vt" \
  'BEGIN{print (g>=c-1e-9 && c>=w-1e-9 && w>=t-1e-9) ? "yes" : "no"}')
check "graph >= called >= carrier >= genotype ($vg/$vc/$vw/$vt)" "$ordered" "yes"

# On this fixture every carrier is genotyped correctly, so carrier_walk must EQUAL called and the
# assignment term of the loss decomposition must be exactly zero.
check "with correct genotypes, carrier_walk equals called" "$(awk -v a="$vc" -v b="$vw" 'BEGIN{print (a==b)?"eq":"ne"}')" "eq"
check "and carrier-assignment loss is zero" "$(sumcol loss_bp "$AS" carrier_assignment)" "0"

# --- and the negative control, without which the two levels could be identical for a trivial reason.
# Strip every carrier: rewrite each sample's GT subfield to 0, leaving records, nodes and coordinates
# untouched. Discovery is unchanged, so `called` must not move; every carrier is now missing, so
# `carrier_walk` must fall to the do-nothing baseline and the assignment term must absorb the whole
# difference. This is what proves carrier_walk measures ASSIGNMENT and nothing else.
awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print;next}
  {for(i=10;i<=NF;i++){n=split($i,p,":"); p[1]="0"; s=p[1]; for(j=2;j<=n;j++) s=s":"p[j]; $i=s} print}' \
  "$OUT/ac.region.vcf" > "$OUT/nocarrier.vcf"
"$BIN" benchmark -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" -r "$REF" \
   --variant-nodes "$OUT/ac.variant_nodes.tsv" --vcf "$OUT/nocarrier.vcf" -o "$OUT/NC" \
   --no-truth-events -q >/dev/null 2>&1
NS="$OUT/NC.qv_summary.tsv"
check "stripping carriers leaves the called level untouched (discovery is unchanged)" \
      "$(pctcol variation_recovered "$NS" called)" "$vc"
ncw=$(pctcol variation_recovered "$NS" carrier_walk)
check "but carrier_walk collapses to the baseline" \
      "$(awk -v w="$ncw" 'BEGIN{print (w+0==0)?"zero":"nonzero: " w}')" "zero"
check "and the whole loss moves into the carrier-assignment term" \
      "$(awk -v a="$(sumcol loss_bp "$NS" carrier_assignment)" 'BEGIN{print (a+0>0)?"positive":"zero"}')" "positive"
check "with every true carrier now a false negative" \
      "$(awk -v tp="$(sumcol gt_carrier "$NS" TP)" -v fn="$(sumcol gt_carrier "$NS" FN)" \
          'BEGIN{print tp"/"fn}')" "0/2"

# ---------------------------------------------------------------------------------------------
# Fixture B: ONE bubble, one 60 bp deletion, and `call` run at a threshold that refuses to emit it.
# The default must still score the bubble and must report the event as MISSED -- and gap_closed must
# be undefined rather than 100%, since with nothing called the graph bound equals the baseline.
# ---------------------------------------------------------------------------------------------
BF1=$(gblk 100 101); BA=$(gblk 60 202); BF2=$(gblk 100 303)
BREF='ref#0#chr2:1-260'
{ printf 'H\tVN:Z:1.0\n'
  printf "S\tg1\t%s\nS\tga\t%s\nS\tg2\t%s\n" "$BF1" "$BA" "$BF2"
  for e in "g1 ga" "ga g2" "g1 g2"; do set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"; done
  printf 'P\t%s\tg1+,ga+,g2+\t*\n' "$BREF"
  printf 'P\thapref#1#chr2\tg1+,ga+,g2+\t*\n'
  printf 'P\thapdel#1#chr2\tg1+,g2+\t*\n'; } > "$OUT/b.gfa"
"$BIN" bubble -i "$OUT/b.gfa" -r "$BREF" -o "$OUT/bb" -q >/dev/null 2>&1
"$BIN" call -i "$OUT/bb.sorted.gfa" -c "$OUT/bb.bubbles.csv" -r "$BREF" -o "$OUT/bc" \
      --min-sv-bp 200 -q >/dev/null 2>&1
nrec=$(grep -v '^#' "$OUT/bc.region.vcf" 2>/dev/null | wc -l | tr -d ' ')
check "call at --min-sv-bp 200 emits no record for a 60 bp deletion" "$nrec" "0"

"$BIN" benchmark -i "$OUT/bb.sorted.gfa" -c "$OUT/bb.bubbles.csv" -r "$BREF" \
   --variant-nodes "$OUT/bc.variant_nodes.tsv" --vcf "$OUT/bc.region.vcf" -o "$OUT/B" -q >/dev/null 2>&1
BS="$OUT/B.qv_summary.tsv"
nscored=$(awk -F'\t' '$1=="excluded" && $3=="bubbles_no_reference_walk"{print $4}' "$BS")
check "the uncalled bubble is still reference-traversed" "$nscored" "0"
check "the uncalled 60 bp deletion is reported MISSED" "$(sumcol truth_event "$BS" missed)" "1"
check "and its bases are counted" "$(sumcol truth_bp "$BS" missed_bp)" "60"
check "gap_closed is UNDEFINED, not 100%" "$(pctcol gt_gap "$BS" gap_closed_pooled)" "NA"
check "and the undefined haplotypes are counted" "$(sumcol gt_gap "$BS" gap_closed_undefined)" "2"

# The reproducer for why the ledger had to replace the run split. Nothing is called here, so the
# reconstruction IS the reference and the whole 60 bp deletion is residual. edlib's co-optimal edit
# path scatters it over chance matches, so EVERY base of a single clean 60 bp event lands in the
# "shorter than the threshold" bucket -- which is what the old classification called
# "variation that could not have been called".
lt=$(sumcol residual_run "$BS" short_runs)
ge=$(sumcol residual_run "$BS" long_runs)
check "the alignment-run split puts the whole 60 bp deletion in short runs" "$lt/$ge" "60/0"
check "while the ledger, sized from the walks, calls it one 60 bp event" \
      "$(sumcol truth_bp "$BS" missed_bp)/$(sumcol truth_event "$BS" missed)" "60/1"

# --called-bubbles-only reproduces the old blind spot, and must say so by scoring nothing.
"$BIN" benchmark -i "$OUT/bb.sorted.gfa" -c "$OUT/bb.bubbles.csv" -r "$BREF" \
   --variant-nodes "$OUT/bc.variant_nodes.tsv" --vcf "$OUT/bc.region.vcf" -o "$OUT/Bc" \
   --called-bubbles-only -q >/dev/null 2>&1
cm=$(sumcol truth_event "$OUT/Bc.qv_summary.tsv" missed)
check "--called-bubbles-only cannot see the miss (that is what it is for)" "${cm:-0}" "0"

# ---------------------------------------------------------------------------------------------
# Fixture C: a symbolic INV must be reconstructed, not counted unhandled.
# ---------------------------------------------------------------------------------------------
CF1=$(gblk 100 7); CI=$(gblk 80 8); CF2=$(gblk 100 9)
CREF='ref#0#chr3:1-280'
{ printf 'H\tVN:Z:1.0\n'
  printf "S\th1\t%s\nS\thi\t%s\nS\th2\t%s\n" "$CF1" "$CI" "$CF2"
  for e in "h1 + hi +" "hi + h2 +" "h1 + hi -" "hi - h2 +"; do
    set -- $e; printf "L\t%s\t%s\t%s\t%s\t0M\n" "$1" "$2" "$3" "$4"; done
  printf 'P\t%s\th1+,hi+,h2+\t*\n' "$CREF"
  printf 'P\thapref#1#chr3\th1+,hi+,h2+\t*\n'
  printf 'P\thapinv#1#chr3\th1+,hi-,h2+\t*\n'; } > "$OUT/c.gfa"
"$BIN" bubble -i "$OUT/c.gfa" -r "$CREF" -o "$OUT/cb" -q >/dev/null 2>&1
"$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$CREF" -o "$OUT/cc" -q >/dev/null 2>&1
ninv=$(awk -F'\t' '$0!~/^#/ && $8~/SVTYPE=INV/' "$OUT/cc.region.vcf" | wc -l | tr -d ' ')
check "fixture C yields a symbolic INV record" "$ninv" "1"
"$BIN" benchmark -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$CREF" \
   --variant-nodes "$OUT/cc.variant_nodes.tsv" --vcf "$OUT/cc.region.vcf" -o "$OUT/C" -q >/dev/null 2>&1
CS="$OUT/C.qv_summary.tsv"
check "the INV is applied, not left unhandled" "$(sumcol gt_records "$CS" unhandled)" "0"
check "and it is applied exactly (no heuristic tiling)" "$(sumcol gt_records "$CS" heuristic)" "0"
# Reconstructing an inversion from the reference is exact, so the VCF must reach the graph bound.
check "an inverted haplotype reconstructs perfectly from the VCF" "$(sumcol gt_gap "$CS" genotype_delta)" "0"

# ---------------------------------------------------------------------------------------------
# Input contracts. Each must be REFUSED, and must leave no output behind.
# ---------------------------------------------------------------------------------------------
refused() { # refused <label> <outprefix> <command...>
  local label="$1" pfx="$2"; shift 2
  "$@" >/dev/null 2>&1
  local rc=$?
  local left=0
  for s in .qv.tsv .qv_by_haplotype.tsv .qv_summary.tsv .truth_events.tsv; do
    [ -e "$pfx$s" ] && left=1
  done
  if [ $rc -ne 0 ] && [ $left -eq 0 ]; then ok "$label"
  elif [ $rc -eq 0 ]; then bad "$label: accepted (exit 0)"
  else bad "$label: refused but left output behind"; fi
}

# a bubbles CSV describing a different graph
awk -F, 'BEGIN{OFS=","} NR==2{$2="NOSUCHNODE"} {print}' "$OUT/ab.bubbles.csv" > "$OUT/wrong.bubbles.csv"
refused "a bubbles CSV naming a node the graph does not have is refused" "$OUT/E1" \
  "$BIN" benchmark -i "$OUT/ab.sorted.gfa" -c "$OUT/wrong.bubbles.csv" -r "$REF" \
     --variant-nodes "$OUT/ac.variant_nodes.tsv" -o "$OUT/E1" -q

# an output that is also an input
# An output that is also an input. Here the file to protect IS the collision, so "left no output
# behind" is the wrong check -- what must hold is that the input survives untouched.
cp "$OUT/ac.variant_nodes.tsv" "$OUT/coll.qv.tsv"
"$BIN" benchmark -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" -r "$REF" \
   --variant-nodes "$OUT/coll.qv.tsv" -o "$OUT/coll" -q >/dev/null 2>&1
rc=$?
if [ $rc -ne 0 ] && cmp -s "$OUT/coll.qv.tsv" "$OUT/ac.variant_nodes.tsv"; then
  ok "an output path that is also an input is refused, and the input is untouched"
elif [ $rc -eq 0 ]; then bad "an output that is also an input was accepted"
else bad "an output that is also an input was refused AFTER overwriting the input"; fi

# a VCF with a duplicated sample column
awk -F'\t' 'BEGIN{OFS="\t"} /^#CHROM/{$NF=$(NF-1)} {print}' "$OUT/ac.region.vcf" > "$OUT/dup.vcf"
refused "a VCF with a duplicate sample column is refused" "$OUT/E3" \
  "$BIN" benchmark -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" -r "$REF" \
     --variant-nodes "$OUT/ac.variant_nodes.tsv" --vcf "$OUT/dup.vcf" -o "$OUT/E3" -q

# --- variant_nodes.tsv contracts. Each of these previously ran to completion with exit 0, and the
# first two SILENTLY RECLASSIFIED a called event as missed -- the same output a genuine caller miss
# produces, which is what makes them worth pinning individually.
awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print;next} {$2=$2+900000; print}' "$OUT/ac.variant_nodes.tsv" > "$OUT/other.tsv"
refused "variant nodes whose bubble ids are absent from the CSV are refused" "$OUT/E4" \
  "$BIN" benchmark -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" -r "$REF" \
     --variant-nodes "$OUT/other.tsv" -o "$OUT/E4" -q

# ONE stale node, real bubble id kept: turned called=1 into missed=1 with exit 0.
awk -F'\t' 'BEGIN{OFS="\t"} NR>1{$4="NOSUCHNODE"} {print}' "$OUT/ac.variant_nodes.tsv" > "$OUT/stalenode.tsv"
refused "a variant_nodes node the bubble does not contain is refused" "$OUT/E4a" \
  "$BIN" benchmark -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" -r "$REF" \
     --variant-nodes "$OUT/stalenode.tsv" -o "$OUT/E4a" -q

# A node that EXISTS in the graph but belongs to another bubble: graph membership alone would pass it.
foreign=$(awk -F'\t' -v B="$(awk -F'\t' 'NR==2{print $2}' "$OUT/ac.variant_nodes.tsv")" \
  'NR>1 && $2!=B {print $4; exit}' "$OUT/ac.variant_nodes.tsv")
if [ -n "$foreign" ]; then
  awk -F'\t' -v N="$foreign" 'BEGIN{OFS="\t"} NR==2{$4=N} {print}' "$OUT/ac.variant_nodes.tsv" > "$OUT/foreign.tsv"
  refused "a node from a DIFFERENT bubble is refused, not just one absent from the graph" "$OUT/E4b" \
    "$BIN" benchmark -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" -r "$REF" \
       --variant-nodes "$OUT/foreign.tsv" -o "$OUT/E4b" -q
else
  bad "fixture A did not yield two bubbles with distinct called nodes"
fi

# No header: the first data row is eaten as one, dropping a real call.
tail -n +2 "$OUT/ac.variant_nodes.tsv" > "$OUT/nohdr.tsv"
refused "a variant_nodes file with no header is refused" "$OUT/E4c" \
  "$BIN" benchmark -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" -r "$REF" \
     --variant-nodes "$OUT/nohdr.tsv" -o "$OUT/E4c" -q

# A short row used to be skipped silently.
{ head -1 "$OUT/ac.variant_nodes.tsv"; printf 'truncated\t1\n'; tail -n +2 "$OUT/ac.variant_nodes.tsv"; } > "$OUT/short.tsv"
refused "a malformed variant_nodes row is refused, not skipped" "$OUT/E4d" \
  "$BIN" benchmark -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" -r "$REF" \
     --variant-nodes "$OUT/short.tsv" -o "$OUT/E4d" -q

{ cat "$OUT/ac.variant_nodes.tsv"; tail -1 "$OUT/ac.variant_nodes.tsv"; } > "$OUT/dupid.tsv"
refused "a duplicate variant_id is refused" "$OUT/E4e" \
  "$BIN" benchmark -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" -r "$REF" \
     --variant-nodes "$OUT/dupid.tsv" -o "$OUT/E4e" -q

# --- VCF contracts. A diploid GT is the dangerous one: stoi("0/1") is 0, so a heterozygous carrier
# scored as reference with nothing said.
sed 's/\t0:/\t0\/1:/g' "$OUT/ac.region.vcf" > "$OUT/dip.vcf"
refused "a diploid GT is refused, not read as its first allele" "$OUT/E4f" \
  "$BIN" benchmark -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" -r "$REF" \
     --variant-nodes "$OUT/ac.variant_nodes.tsv" --vcf "$OUT/dip.vcf" -o "$OUT/E4f" -q

sed 's/BUBBLE_ID=[0-9]*;//' "$OUT/ac.region.vcf" > "$OUT/nobid.vcf"
refused "a VCF record with no BUBBLE_ID is refused" "$OUT/E4g" \
  "$BIN" benchmark -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" -r "$REF" \
     --variant-nodes "$OUT/ac.variant_nodes.tsv" --vcf "$OUT/nobid.vcf" -o "$OUT/E4g" -q

awk -F'\t' 'BEGIN{OFS="\t"} !/^#/ && NR>0 {NF=NF-1} {print}' "$OUT/ac.region.vcf" > "$OUT/shortvcf.vcf"
refused "a VCF data line with the wrong field count is refused" "$OUT/E4h" \
  "$BIN" benchmark -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" -r "$REF" \
     --variant-nodes "$OUT/ac.variant_nodes.tsv" --vcf "$OUT/shortvcf.vcf" -o "$OUT/E4h" -q

# an ambiguous reference query, resolved through the SHARED resolver rather than a private copy
refused "an ambiguous --reference-path is refused, not resolved by file order" "$OUT/E5" \
  "$BIN" benchmark -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" -r "chr1" \
     --variant-nodes "$OUT/ac.variant_nodes.tsv" -o "$OUT/E5" -q

# ---------------------------------------------------------------------------------------------
# Determinism: the same run at 1 and 8 threads must be byte-identical.
# ---------------------------------------------------------------------------------------------
"$BIN" benchmark -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" -r "$REF" \
   --variant-nodes "$OUT/ac.variant_nodes.tsv" --vcf "$OUT/ac.region.vcf" -o "$OUT/T1" --threads 1 -q >/dev/null 2>&1
"$BIN" benchmark -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" -r "$REF" \
   --variant-nodes "$OUT/ac.variant_nodes.tsv" --vcf "$OUT/ac.region.vcf" -o "$OUT/T8" --threads 8 -q >/dev/null 2>&1
same=1
for s in .qv.tsv .qv_by_haplotype.tsv .qv_summary.tsv .truth_events.tsv; do
  cmp -s "$OUT/T1$s" "$OUT/T8$s" || same=0
done
check "output is byte-identical at 1 and 8 threads" "$same" "1"

# ---------------------------------------------------------------------------------------------
# The allele VCF is a serialization ceiling, not a call-sensitivity score: it must reach 0 residual
# even where the region VCF does not, because it spells every allele out.
# ---------------------------------------------------------------------------------------------
"$BIN" call -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" -r "$REF" -o "$OUT/av" \
      --allele-vcf -q >/dev/null 2>&1
if [ -s "$OUT/av.alleles.vcf" ]; then
  "$BIN" benchmark -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" -r "$REF" \
     --variant-nodes "$OUT/av.variant_nodes.tsv" --vcf "$OUT/av.alleles.vcf" -o "$OUT/AV" -q >/dev/null 2>&1
  check "the allele VCF reconstructs every haplotype exactly" \
        "$(sumcol gt_gap "$OUT/AV.qv_summary.tsv" genotype_delta)" "0"
else
  bad "call --allele-vcf wrote no alleles VCF"
fi

printf "%d assertion(s) failed\n" "$fails"
[ "$fails" -eq 0 ]
