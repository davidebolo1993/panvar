#!/usr/bin/env bash
# describe_stats.sh - contract assertions for `panvar describe`.
#
# Small fixtures with known answers: which paths cross a site (a deletion's direct route included), which
# features are kept, the dosages and NA in every matrix, the sample level, the variant matrix, and the
# refusals.
#
# Every expected number below is derivable by hand from the fixture.
#
#   describe_stats.sh <panvar-binary> <out-dir>
set -uo pipefail

BIN="${1:?usage: describe_stats.sh <panvar> <outdir>}"
OUT="${2:?}"
mkdir -p "$OUT"
OUT="$OUT/run.$$.$(date +%s)"
rm -rf "$OUT"; mkdir -p "$OUT"
fails=0
ok()  { printf "  ok   %s\n" "$1"; }
bad() { printf "  FAIL %s\n" "$1"; fails=$((fails + 1)); }
check() { [ "$2" = "$3" ] && ok "$1 = $3" || bad "$1: expected $3, got ${2:-<missing>}"; }

# Deterministic pseudo-random sequence (LCG): distinct, non-repetitive node content.
gblk() { local n=$1 s=${2:-1} out="" i=0 b
  while [ $i -lt "$n" ]; do s=$(( (s * 1103515245 + 12345) % 2147483648 ))
    b=$(( (s / 65536) % 4 ))
    case $b in 0) out="${out}A";; 1) out="${out}C";; 2) out="${out}G";; *) out="${out}T";; esac
    i=$((i+1)); done; printf '%s' "$out"; }

# dosage <bimbam.gz> <feature-id> <column-index-from-1>
dosage() { gunzip -c "$1" | awk -F', ' -v f="$2" -v c="$3" '$1==f {print $(c+3); exit}'; }
colof()  { gunzip -c "$1" | grep -n "^$2$" | cut -d: -f1; }

refused() { # refused <label> <outdir> <command...>
  local label="$1" dir="$2"; shift 2
  "$@" >/dev/null 2>&1
  local rc=$?
  if [ $rc -ne 0 ]; then ok "$label"; else bad "$label: accepted (exit 0)"; fi
}

# ---------------------------------------------------------------------------------------------
# Fixture A: two bubbles. Bubble 1 has an interior node some haplotypes DELETE outright (the direct
# source->sink edge); bubble 2 is the same shape further along. `trunc` stops before bubble 2, so it
# is genuinely unobservable there -- the control that keeps NA meaningful.
#
#   f1 200 | a 120 | f2 200 | f3 200 | b 120 | f4 200
# ---------------------------------------------------------------------------------------------
F1=$(gblk 200 1); A=$(gblk 120 2); F2=$(gblk 200 3)
F3=$(gblk 200 4); B=$(gblk 120 5); F4=$(gblk 200 6)
REF='ref#0#chr1:1-1040'
{ printf 'H\tVN:Z:1.0\n'
  printf "S\tf1\t%s\nS\ta\t%s\nS\tf2\t%s\nS\tf3\t%s\nS\tb\t%s\nS\tf4\t%s\n" "$F1" "$A" "$F2" "$F3" "$B" "$F4"
  for e in "f1 a" "a f2" "f1 f2" "f2 f3" "f3 b" "b f4" "f3 f4"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"; done
  printf 'P\t%s\tf1+,a+,f2+,f3+,b+,f4+\t*\n' "$REF"
  for i in 1 2 3; do printf 'P\tfull%s#1#chr1\tf1+,a+,f2+,f3+,b+,f4+\t*\n' "$i"; done
  for i in 1 2 3; do printf 'P\tdel%s#1#chr1\tf1+,f2+,f3+,f4+\t*\n' "$i"; done
  printf 'P\ttrunc#1#chr1\tf1+,a+,f2+\t*\n'; } > "$OUT/a.gfa"

"$BIN" bubble -i "$OUT/a.gfa" -r "$REF" -o "$OUT/ab" --merge-nearby-bp 0 -q >/dev/null 2>&1
check "fixture A presents two bubbles" "$(tail -n +2 "$OUT/ab.bubbles.csv" | wc -l | tr -d ' ')" "2"

"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/A" \
   --no-wide-matrix -q >/dev/null 2>&1
AI="$OUT/A/describe.index.tsv"
GB="$OUT/A/haplotype/graph/bimbam_graph.bimbam.gz"
GS="$OUT/A/haplotype/graph/samples.txt.gz"

# A deletion allele crosses the site: ref + 3 full + 3 del + trunc = 8, every one crossing bubble 1
# (trunc ends just after it).
check "bubble 1 counts all 8 paths, not just the 5 that keep the interior" \
      "$(awk -F'\t' 'NR==2{print $3}' "$AI")" "8"
check "and it yields features instead of nothing" \
      "$(awk -F'\t' 'NR==2{print ($13>0 && $16>0) ? "yes" : "no"}' "$AI")" "yes"

# ...and the features are the right ones: a matrix of flank k-mers every haplotype carries would not be
# empty either. The signal a deletion contributes is the junction it creates: syncmers spanning the
# direct source->sink join, carried by every bypass haplotype and by no path that keeps the interior.
KB="$OUT/A/haplotype/kmers/bimbam_kmers.bimbam.gz"
KS="$OUT/A/haplotype/kmers/samples.txt.gz"
alt_exclusive=$(gunzip -c "$KB" | awk -F', *' -v samples="$(gunzip -c "$KS" | tr '\n' ' ')" '
  BEGIN { n = split(samples, s, " ") }
  {
    del_all = 1; other_any = 0
    for (i = 1; i <= n; i++) {
      v = $(i + 3)
      carries = (v != "NA" && v != "" && v + 0 > 0)
      if (s[i] ~ /^del/) { if (!carries) del_all = 0 }
      else if (carries)  { other_any = 1 }
    }
    if (del_all && !other_any) c++
  }
  END { print c + 0 }')
[ "${alt_exclusive:-0}" -ge 2 ]   && ok "the deletion contributes ALT-exclusive junction markers ($alt_exclusive of them)"   || bad "expected >=2 ALT-exclusive junction syncmers, found ${alt_exclusive:-0}"


# --- zero vs missing, in both directions, on hand-checkable cells.
cdel=$(colof "$GS" 'del1#1#chr1'); cfull=$(colof "$GS" 'full1#1#chr1'); ctr=$(colof "$GS" 'trunc#1#chr1')
b1node=$(awk -F',' 'NR==2{gsub(/"/,"",$15); split($15,x,";"); print x[1]}' "$OUT/ab.bubbles.csv")
check "a deletion haplotype reads 0 at the deleted node (traverses, does not carry)" \
      "$(dosage "$GB" "$b1node" "$cdel")" "0"
check "a carrier reads 1 there" "$(dosage "$GB" "$b1node" "$cfull")" "1"
b2node=$(awk -F',' 'NR==3{gsub(/"/,"",$15); split($15,x,";"); print x[1]}' "$OUT/ab.bubbles.csv")
check "a path that never reaches bubble 2 reads NA there, not 0" \
      "$(dosage "$GB" "$b2node" "$ctr")" "NA"
check "and it reads a real dosage at the bubble it does reach" \
      "$(dosage "$GB" "$b1node" "$ctr")" "1"

# --- the sidecar names the real subtype, which is what separates a node from an edge feature.
ENC="$OUT/A/haplotype/graph/feature_annot.graph.tsv.gz"
check "graph features are labelled node/edge, not 'count'" \
      "$(gunzip -c "$ENC" | awk -F'\t' 'NR>1{print $3}' | sort -u | tr '\n' ',')" "edge,node,"
check "k-mer features are labelled by their sampling mode" \
      "$(gunzip -c "$OUT/A/haplotype/kmers/feature_annot.kmers.tsv.gz" | awk -F'\t' 'NR>1{print $3}' | sort -u)" "syncmer"

# --- determinism.
"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/T1" \
   --no-wide-matrix --threads 1 -q >/dev/null 2>&1
"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/T8" \
   --no-wide-matrix --threads 8 -q >/dev/null 2>&1
same=1
for f in haplotype/graph/bimbam_graph.bimbam.gz haplotype/kmers/bimbam_kmers.bimbam.gz \
         haplotype/graph/feature_annot.graph.tsv.gz haplotype/kmers/feature_annot.kmers.tsv.gz; do
  cmp -s <(gunzip -c "$OUT/T1/$f") <(gunzip -c "$OUT/T8/$f") || same=0
done
check "output is byte-identical at 1 and 8 threads" "$same" "1"

# ---------------------------------------------------------------------------------------------
# Fixture B: a short variant beside very long neighbours. --variant-flank-bp admits bases, not whole
# nodes: at 30 it takes 30 bases of each 4000 bp neighbour, not all their k-mers.
# kmer_candidates is pre-filter, so it tracks kept bases; feature counts do not, because flank k-mers
# are identical across paths and the discriminative filter removes them.
# ---------------------------------------------------------------------------------------------
LF=$(gblk 4000 11); V=$(gblk 60 12); VA=$(gblk 60 13); RF=$(gblk 4000 14)
BREF='ref#0#chr2:1-8060'
{ printf 'H\tVN:Z:1.0\n'
  printf "S\tg1\t%s\nS\tv\t%s\nS\tvalt\t%s\nS\tg2\t%s\n" "$LF" "$V" "$VA" "$RF"
  for e in "g1 v" "v g2" "g1 valt" "valt g2"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"; done
  printf 'P\t%s\tg1+,v+,g2+\t*\n' "$BREF"
  for i in 1 2 3; do printf 'P\tr%s#1#chr2\tg1+,v+,g2+\t*\n' "$i"; done
  for i in 1 2 3; do printf 'P\tq%s#1#chr2\tg1+,valt+,g2+\t*\n' "$i"; done; } > "$OUT/b.gfa"
"$BIN" bubble -i "$OUT/b.gfa" -r "$BREF" -o "$OUT/bb" --merge-nearby-bp 0 -q >/dev/null 2>&1
bid=$(tail -n +2 "$OUT/bb.bubbles.csv" | cut -d, -f1 | head -1)
ins=$(tail -n +2 "$OUT/bb.bubbles.csv" | head -1 | sed 's/.*"\(.*\)".*/\1/' | tr ';' ',')
printf 'variant_id\tbubble_id\tsvtype\tnode_ids\nv1\t%s\tINS\t%s\n' "$bid" "$ins" > "$OUT/vn.tsv"
cand() { awk -F'\t' 'NR==2{print $4}' "$1/describe.index.tsv"; }
for fl in 0 30 200; do
  "$BIN" describe -i "$OUT/bb.sorted.gfa" --bubble-prefix "$OUT/bb" --out-dir "$OUT/F$fl" \
     --variant-nodes "$OUT/vn.tsv" --variant-flank-bp $fl --no-wide-matrix -q >/dev/null 2>&1
done
c0=$(cand "$OUT/F0"); c30=$(cand "$OUT/F30"); c200=$(cand "$OUT/F200")
check "the flank admits BASES: candidates grow with it and stay far below the 4000 bp neighbours" \
      "$(awk -v a="$c0" -v b="$c30" -v c="$c200" 'BEGIN{print (a<b && b<c && c<500) ? "yes" : "no ("a","b","c")"}')" "yes"
check "the graph substrate stays node-granular (same node count at 30 and 200)" \
      "$(awk -F'\t' 'NR==2{print $12}' "$OUT/F30/describe.index.tsv")" \
      "$(awk -F'\t' 'NR==2{print $12}' "$OUT/F200/describe.index.tsv")"

# ---------------------------------------------------------------------------------------------
# Input contracts. Each must be refused; none may be silently skipped or repaired.
# ---------------------------------------------------------------------------------------------
refused "a --bubble-id that does not exist is refused" "$OUT/E1" \
  "$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/E1" --bubble-id 999 -q

awk -F, 'BEGIN{OFS=","} NR==2{$2="NOSUCHNODE"} {print}' "$OUT/ab.bubbles.csv" > "$OUT/badcsv.csv"
refused "a bubbles CSV naming a node the graph lacks is refused" "$OUT/E2" \
  "$BIN" describe -i "$OUT/ab.sorted.gfa" -c "$OUT/badcsv.csv" --out-dir "$OUT/E2" -q

# variant_nodes goes through the SHARED strict reader; a stale node silently removes real features.
printf 'variant_id\tbubble_id\tsvtype\tnode_ids\nv1\t%s\tINS\tNOSUCHNODE\n' "$bid" > "$OUT/vnbad.tsv"
refused "a variant_nodes node the bubble does not contain is refused" "$OUT/E3" \
  "$BIN" describe -i "$OUT/bb.sorted.gfa" --bubble-prefix "$OUT/bb" --out-dir "$OUT/E3" \
     --variant-nodes "$OUT/vnbad.tsv" -q
tail -n +2 "$OUT/vn.tsv" > "$OUT/vnnohdr.tsv"
refused "a variant_nodes file with no header is refused" "$OUT/E4" \
  "$BIN" describe -i "$OUT/bb.sorted.gfa" --bubble-prefix "$OUT/bb" --out-dir "$OUT/E4" \
     --variant-nodes "$OUT/vnnohdr.tsv" -q

# --samples: a duplicate sample id makes the haplotype assignment depend on file order.
printf 'sample\thap1\thap2\ns1\tfull1#1#chr1\tdel1#1#chr1\ns1\tfull2#1#chr1\tdel2#1#chr1\n' > "$OUT/dupsam.tsv"
refused "a duplicate sample id in --samples is refused" "$OUT/E5" \
  "$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/E5" \
     --samples "$OUT/dupsam.tsv" -q

"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/E6" \
     --only-kmers --only-graph --no-wide-matrix -q >/dev/null 2>&1
check "two --only-* flags combine (k-mers and graph)" "$(ls "$OUT/E6/haplotype" 2>/dev/null | sort | tr '\n' ' ')" "graph kmers "
refused "--variant-flank-bp without --variant-nodes is refused as inert" "$OUT/E7" \
  "$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/E7" \
     --variant-flank-bp 30 -q

# ---------------------------------------------------------------------------------------------
# Transactionality: a rerun with a different substrate selection must not leave the old family
# looking current, a failing run must not disturb what is there, and neither may touch user files.
# ---------------------------------------------------------------------------------------------
TX="$OUT/tx"
"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$TX" --no-wide-matrix -q >/dev/null 2>&1
check "a full run writes both graph substrates" "$(ls "$TX/haplotype" | sort | tr '\n' ' ')" "graph kmers "
touch "$TX/USER_FILE.txt"
"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$TX" --only-graph --no-wide-matrix -q >/dev/null 2>&1
check "a --only-graph rerun removes the now-stale k-mer family" "$(ls "$TX/haplotype" | sort | tr '\n' ' ')" "graph "
check "and leaves an unrelated user file alone" "$([ -f "$TX/USER_FILE.txt" ] && echo kept || echo REMOVED)" "kept"
before=$(ls "$TX" | sort | tr '\n' ' ')
"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$TX" --bubble-id 999 -q >/dev/null 2>&1
check "a failing run leaves the directory untouched" "$(ls "$TX" | sort | tr '\n' ' ')" "$before"
check "and leaves no staging directory behind" "$(ls -d "$TX".describe-staging.* 2>/dev/null | wc -l | tr -d ' ')" "0"

# ---------------------------------------------------------------------------------------------
# --variant-nodes keeps the deletion. A DEL's EVENT_NODES are the nodes its carriers lack, so on a carrier
# only the site's boundaries give it a positive signal: the bypass edge source->sink and the junction
# k-mers.
# ---------------------------------------------------------------------------------------------
dbid=$(tail -n +2 "$OUT/ab.bubbles.csv" | head -1 | cut -d, -f1)
dins=$(tail -n +2 "$OUT/ab.bubbles.csv" | head -1 | sed 's/.*"\(.*\)".*/\1/' | tr ';' ',')
dsrc=$(tail -n +2 "$OUT/ab.bubbles.csv" | head -1 | cut -d, -f2)
dsnk=$(tail -n +2 "$OUT/ab.bubbles.csv" | head -1 | cut -d, -f4)
printf 'variant_id\tbubble_id\tsvtype\tnode_ids\nd1\t%s\tDEL\t%s\n' "$dbid" "$dins" > "$OUT/vndel.tsv"
"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/VD" \
   --variant-nodes "$OUT/vndel.tsv" --variant-flank-bp 30 --no-wide-matrix -q >/dev/null 2>&1
VDB="$OUT/VD/haplotype/graph/bimbam_graph.bimbam.gz"
check "the ALT-specific bypass edge survives --variant-nodes" \
      "$(gunzip -c "$VDB" | cut -d, -f1 | grep -c "^${dsrc}+>${dsnk}+$")" "1"
check "and the deletion still yields junction k-mers" \
      "$(awk -v a="$(gunzip -c "$OUT/VD/haplotype/kmers/bimbam_kmers.bimbam.gz" | wc -l)" 'BEGIN{print (a+0>0)?"yes":"no"}')" "yes"

# ---------------------------------------------------------------------------------------------
# Sample level: a diploid genotype is complete or it is NA. A sample with one haplotype unobserved is NA,
# not half its dosage, and stays out of AC and AN.
# ---------------------------------------------------------------------------------------------
"$BIN" call -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" -r "$REF" -o "$OUT/ac" -q >/dev/null 2>&1
printf 'sample\thap1\thap2\nS1\tfull1#1#chr1\tdel1#1#chr1\nS2\tfull2#1#chr1\ttrunc#1#chr1\n' > "$OUT/sam.tsv"
"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/SV" \
   --variant-vcf "$OUT/ac.region.vcf" --samples "$OUT/sam.tsv" --no-wide-matrix -q >/dev/null 2>&1
SVB="$OUT/SV/sample/variant/bimbam_variant.bimbam.gz"
if [ -s "$SVB" ]; then
  s2col=$(gunzip -c "$OUT/SV/sample/variant/samples.txt.gz" | grep -n '^S2$' | cut -d: -f1)
  # S2's second haplotype does not traverse bubble 2, so at a bubble-2 record S2 must be NA.
  na2=$(gunzip -c "$SVB" | awk -F', ' -v c="$s2col" 'NR>0 && $(c+3)=="NA"{n++} END{print n+0}')
  check "a half-observed diploid sample is NA somewhere, not a halved dosage" \
        "$(awk -v n="$na2" 'BEGIN{print (n+0>0)?"yes":"no"}')" "yes"
else
  bad "no sample-level variant BIMBAM was written"
fi

# --samples must not accept a haplotype that matches nothing: it contributes silently nothing, which
# looks exactly like a sample with no signal.
printf 'sample\thap1\thap2\nS1\tfull1#1#chr1\tGHOST#9#chr1\n' > "$OUT/ghost.tsv"
refused "an unknown haplotype name in --samples is refused" "$OUT/E8" \
  "$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/E8" \
     --samples "$OUT/ghost.tsv" -q

# --no-bimbam must mean no BIMBAM anywhere, not just the pooled graph/k-mer ones.
refused "--no-bimbam with --variant-vcf is refused" "$OUT/E9" \
  "$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/E9" \
     --no-bimbam --variant-vcf "$OUT/ac.region.vcf" -q
refused "--no-bimbam with --samples is refused" "$OUT/E10" \
  "$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/E10" \
     --no-bimbam --samples "$OUT/sam.tsv" -q

# --variant-vcf parser: each of these is refused rather than read some way.
mkvcf() { grep -v '^##' "$OUT/ac.region.vcf" > /dev/null 2>&1; }
awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print;next} {$9="CN"; print}' "$OUT/ac.region.vcf" > "$OUT/nogt.vcf"
refused "a FORMAT with no GT is refused, not read from the first field" "$OUT/E11" \
  "$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/E11" \
     --variant-vcf "$OUT/nogt.vcf" --only-variant -q
awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print;next} {for(i=10;i<=NF;i++){n=split($i,p,":"); p[1]="7"; s=p[1]; for(j=2;j<=n;j++) s=s":"p[j]; $i=s} print}' \
  "$OUT/ac.region.vcf" > "$OUT/badidx.vcf"
refused "an allele index the record does not have is refused" "$OUT/E12" \
  "$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/E12" \
     --variant-vcf "$OUT/badidx.vcf" --only-variant -q
awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print;next} {$8=$8";NALLELES=abc"; print}' "$OUT/ac.region.vcf" > "$OUT/badnall.vcf"
refused "a malformed NALLELES is refused, not atoi'd to a number" "$OUT/E13" \
  "$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/E13" \
     --variant-vcf "$OUT/badnall.vcf" --only-variant -q
: > "$OUT/empty.vcf"
refused "a headerless/empty variant VCF is refused, not a successful empty run" "$OUT/E14" \
  "$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/E14" \
     --variant-vcf "$OUT/empty.vcf" --only-variant -q

# Output ownership is exact: only bubble_<digits> is generated, so a user's bubble_notes must survive.
TXO="$OUT/own"
"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$TXO" --no-wide-matrix -q >/dev/null 2>&1
mkdir -p "$TXO/bubble_notes"; touch "$TXO/bubble_notes/mine.txt" "$TXO/bubble_readme.md"
"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$TXO" --only-graph --no-wide-matrix -q >/dev/null 2>&1
check "a user's bubble_notes/ is not treated as generated output" \
      "$([ -f "$TXO/bubble_notes/mine.txt" ] && echo kept || echo DELETED)" "kept"
check "nor is bubble_readme.md" "$([ -f "$TXO/bubble_readme.md" ] && echo kept || echo DELETED)" "kept"
check "while the generated bubble_1 is still installed" "$([ -d "$TXO/bubble_1" ] && echo yes || echo no)" "yes"

# ---------------------------------------------------------------- variant-substrate input contracts
# A minimal --variant-vcf: two haplotype columns, one biallelic record. `describe` needs BUBBLE_ID in
# INFO and GT first in FORMAT; everything else here exists to be made wrong one field at a time.
vcf_head() {
  printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\thapA\thapB\n'
}
# panvar exits 2 on a refusal, so assert the CLASS of outcome rather than a specific code
verdict() { [ "$1" -eq 0 ] && echo accepted || echo refused; }
VV="$OUT/variant"; mkdir -p "$VV"
{ vcf_head; printf 'chr1\t100\trec1\tA\tG\t.\tPASS\tSVTYPE=SNV;BUBBLE_ID=1\tGT\t1\t0\n'; } > "$VV/ok.vcf"
"$BIN" describe --only-variant --variant-vcf "$VV/ok.vcf" --out-dir "$VV/ok_out" -q >/dev/null 2>&1
check "a well-formed --variant-vcf is accepted" "$(verdict $?)" "accepted"

# --samples names a haplotype the VCF does not have. The graph substrate validates every assignment
# against graph path names; this mode called the same reader without the VCF's column set, so an
# unknown name survived and produced a diploid missing one homologue, reported covered at half dosage.
printf 'sample\thap1\thap2\nS1\thapA\thapZZZ\n' > "$VV/bad_samples.tsv"
"$BIN" describe --only-variant --variant-vcf "$VV/ok.vcf" --samples "$VV/bad_samples.tsv" \
       --out-dir "$VV/bad_out" -q >/dev/null 2>&1
check "--samples naming a haplotype absent from the VCF is refused" "$(verdict $?)" "refused"
printf 'sample\thap1\thap2\nS1\thapA\thapB\n' > "$VV/good_samples.tsv"
"$BIN" describe --only-variant --variant-vcf "$VV/ok.vcf" --samples "$VV/good_samples.tsv" \
       --out-dir "$VV/good_out" -q >/dev/null 2>&1
check "--samples naming only real haplotype columns is accepted" "$(verdict $?)" "accepted"
printf 'sample\thap1\thap2\nS1\thapA\thapB\nS2\thapA\n' > "$VV/mixed_samples.tsv"
check "the mixed-ploidy note carries the log prefix" "$("$BIN" describe --only-variant --variant-vcf "$VV/ok.vcf" --samples "$VV/mixed_samples.tsv" --out-dir "$VV/mixed_out" 2>&1 >/dev/null | grep -c '^\[describe [0-9:]*\] note: --samples has mixed ploidy')" "1"

# NALLELES counts REF plus the ALTs and decides how many dosage rows are emitted and which GT indices
# they test, so a value disagreeing with the ALT column invents or omits alleles.
{ vcf_head; printf 'chr1\t100\trec1\tA\tG,T\t.\tPASS\tBUBBLE_ID=1;NALLELES=4\tGT\t1\t2\n'; } > "$VV/nall.vcf"
"$BIN" describe --only-variant --variant-vcf "$VV/nall.vcf" --out-dir "$VV/nall_out" -q >/dev/null 2>&1
check "NALLELES disagreeing with the ALT cardinality is refused" "$(verdict $?)" "refused"
{ vcf_head; printf 'chr1\t100\trec1\tA\tG,T\t.\tPASS\tBUBBLE_ID=1;NALLELES=3\tGT\t1\t2\n'; } > "$VV/nall_ok.vcf"
"$BIN" describe --only-variant --variant-vcf "$VV/nall_ok.vcf" --out-dir "$VV/nall_ok_out" -q >/dev/null 2>&1
check "NALLELES agreeing with the ALT cardinality is accepted" "$(verdict $?)" "accepted"
# Without NALLELES, the ALT column alone gives one row per ALT, and AF one value per ALT.
{ vcf_head; printf 'chr1\t100\trec1\tA\tG,T\t.\tPASS\tBUBBLE_ID=1;AF=0.5,0.5\tGT\t1\t2\n'; } > "$VV/two_alt.vcf"
"$BIN" describe --only-variant --variant-vcf "$VV/two_alt.vcf" --out-dir "$VV/two_alt_out" -q >/dev/null 2>&1
check "a two-ALT record without NALLELES gives one row per ALT (rows)" \
      "$(gunzip -c "$VV/two_alt_out/haplotype/variant/bimbam_variant.bimbam.gz" 2>/dev/null | tr -d ' ' | tr '\n' ' ')" \
      "rec1_a1,A,B,1,0 rec1_a2,A,B,0,1 "

# A generated per-ALT id and a record ID share one namespace in both directions: a record later in the
# file named rec1_a1 collides with the rec1_a1 generated above it.
{ vcf_head
  printf 'chr1\t100\trec1\tA\tG,T\t.\tPASS\tBUBBLE_ID=1;NALLELES=3\tGT\t1\t2\n'
  printf 'chr1\t200\trec1_a1\tC\tA\t.\tPASS\tBUBBLE_ID=2\tGT\t1\t0\n'; } > "$VV/collide.vcf"
"$BIN" describe --only-variant --variant-vcf "$VV/collide.vcf" --out-dir "$VV/coll_out" -q >/dev/null 2>&1
check "a later record colliding with an earlier generated per-ALT id is refused" "$(verdict $?)" "refused"

# ---------------------------------------------------------------- output transaction
# describe installs its files from a staging directory, file by file. A failed run restores the files it
# replaced and removes the ones it installed.
TXD="$OUT/txd"
"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$TXD/out" \
   --no-wide-matrix -q >/dev/null 2>&1
txd_before=$(ls "$TXD/out" | sort | tr '\n' ',')
PANVAR_TEST_FAIL_COMMIT_AT=2 "$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" \
   --out-dir "$TXD/out" --no-wide-matrix -q >/dev/null 2>&1
check "an injected commit failure exits non-zero" "$(verdict $?)" "refused"
check "the previous describe directory is restored entry for entry" \
      "$(ls "$TXD/out" | sort | tr '\n' ',')" "$txd_before"
check "no describe backup directory is left behind" \
      "$(ls "$TXD" | grep -c describe-backup)" "0"

# A run that fails with NO previous output must leave nothing installed.
PANVAR_TEST_FAIL_COMMIT_AT=2 "$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" \
   --out-dir "$TXD/fresh" --no-wide-matrix -q >/dev/null 2>&1
check "a failed FIRST run installs no partial describe family" \
      "$(ls "$TXD/fresh" 2>/dev/null | wc -l | tr -d ' ')" "0"

# An input that is one of describe's own files would be replaced or removed by the commit; any other file
# in its folders is not describe's, so it is read and left alone.
"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/alias_dir" --no-wide-matrix -q >/dev/null 2>&1
cp "$OUT/ab.bubbles.csv" "$OUT/alias_dir/describe.index.tsv"
for o in "$OUT/alias_dir" "$OUT/alias_dir/"; do
  check "an input that is one of describe's files is refused (-o ${o#"$OUT"/})" \
        "$("$BIN" describe -i "$OUT/ab.sorted.gfa" --bubbles-csv "$OUT/alias_dir/describe.index.tsv" \
           --out-dir "$o" --no-wide-matrix -q 2>&1 | grep -c "is one of describe's own files")" "1"
done
check "and that input still exists" \
      "$(head -c 9 "$OUT/alias_dir/describe.index.tsv")" "$(head -c 9 "$OUT/ab.bubbles.csv")"
cp "$OUT/ab.bubbles.csv" "$OUT/alias_dir/bubble_1/b.csv"
"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubbles-csv "$OUT/alias_dir/bubble_1/b.csv" --out-dir "$OUT/alias_dir" \
   --no-wide-matrix -q >/dev/null 2>&1
check "an input in a site folder under another name is read and kept" \
      "$(verdict $?):$([ -s "$OUT/alias_dir/bubble_1/b.csv" ] && echo kept)" "accepted:kept"

# ---------------------------------------------------------------------------------------------
# Fixture M: two bubbles a..b and d..c with a spacer s between them. h4 and h5 cross both bubbles twice
# (f0 a v1 b s d w1 c a v1 b s d w2 c f9). Their walk across bubble 1 runs from the first a to the last
# b and so would hold s, d, w1 and c too: those belong to bubble 2 and would be counted twice. Fixture P: one 80 bp sequence X is the alternative allele
# of bubble 1 for p0-p3 and of bubble 2 for q0-q3, so each X k-mer is discriminative in each bubble but
# its pooled row reads 1 for everyone.
# ---------------------------------------------------------------------------------------------
mkfx() { # mkfx <prefix> <M|P|G>
  "${PYTHON:-python3}" - "$1" "$2" <<'PYEOF'
import sys, random
pre, kind = sys.argv[1], sys.argv[2]
random.seed(7)
rnd = lambda n: ''.join(random.choice('ACGT') for _ in range(n))
if kind == 'P':
    random.seed(11); X = rnd(80); random.seed(7)
    nodes = {'f0': 60, 'a': 40, 'v1': X, 'v2': 70, 'b': 40, 's': 50, 'd': 40, 'w1': X, 'w2': 65, 'c': 40, 'f9': 60}
    paths = {f'p{i}#1#c': 'f0+,a+,v1+,b+,s+,d+,w2+,c+,f9+' for i in range(4)}
    paths.update({f'q{i}#1#c': 'f0+,a+,v2+,b+,s+,d+,w1+,c+,f9+' for i in range(4)})
else:
    nodes = {'f0': 60, 'a': 40, 'v1': 50, 'v2': 55, 'b': 40, 's': 30, 'd': 40, 'w1': 50, 'w2': 45, 'c': 40, 'f9': 60}
    twice = 'f0+,a+,v1+,b+,s+,d+,w1+,c+,a+,v1+,b+,s+,d+,w2+,c+,f9+'
    paths = {'ref#0#c': 'f0+,a+,v1+,b+,s+,d+,w1+,c+,f9+', 'h1#1#c': 'f0+,a+,v1+,b+,s+,d+,w1+,c+,f9+',
             'h2#1#c': 'f0+,a+,v2+,b+,s+,d+,w2+,c+,f9+', 'h3#1#c': 'f0+,a+,v2+,b+,s+,d+,w2+,c+,f9+',
             'h4#1#c': twice, 'h5#1#c': twice}
    if kind == 'G':   # a node id holding the oriented-edge separator
        nodes['v>1'] = nodes.pop('v1')
        paths = {k: v.replace('v1+', 'v>1+') for k, v in paths.items()}
inside = {1: [n for n in nodes if n[0] == 'v'], 2: [n for n in nodes if n[0] == 'w']}
seqs = {n: v if isinstance(v, str) else rnd(v) for n, v in nodes.items()}
edges = set()
for p in paths.values():
    st = p.split(',')
    for x, y in zip(st, st[1:]): edges.add((x[:-1], x[-1], y[:-1], y[-1]))
with open(pre + '.gfa', 'w') as f:
    f.write('H\tVN:Z:1.0\n')
    for n in nodes: f.write(f'S\t{n}\t{seqs[n]}\n')
    for e in sorted(edges): f.write('L\t%s\t%s\t%s\t%s\t0M\n' % e)
    for name, p in paths.items(): f.write(f'P\t{name}\t{p}\t*\n')
with open(pre + '.bubbles.csv', 'w') as f:
    f.write('bubble_id,source,source_orient,sink,sink_orient,inside_node_count,total_node_count,path_support,'
            'distinct_alleles,ref_allele_support,alt_allele_support_max,alt_allele_support_min,min_inside_bp,'
            'max_inside_bp,inside_nodes\n')
    for bid, s, t in ((1, 'a', 'b'), (2, 'd', 'c')):
        ins = inside[bid]
        f.write(f'{bid},{s},+,{t},+,{len(ins)},{len(ins)+2},{len(paths)},2,1,1,1,0,0,"{";".join(ins)}"\n')
PYEOF
}
mkfx "$OUT/mx" M; mkfx "$OUT/px" P; mkfx "$OUT/gx" G
"$BIN" describe -i "$OUT/mx.gfa" -c "$OUT/mx.bubbles.csv" -o "$OUT/mx_d" --no-wide-matrix --threads 2 -q >/dev/null 2>&1
MG="$OUT/mx_d/haplotype/graph/bimbam_graph.bimbam.gz"
H4=$(colof "$OUT/mx_d/haplotype/graph/samples.txt.gz" 'h4#1#c')
check "a path crossing bubble 1 twice carries v1 twice, a twice, and w1 once (v1:a:w1)" \
      "$(dosage "$MG" v1 "$H4"):$(dosage "$MG" a "$H4"):$(dosage "$MG" w1 "$H4")" "2:2:1"
check "the spacer between the two bubbles is no feature of either" "$(gunzip -c "$MG" | grep -c '^s,')" "0"
check "bubble 1's graph features are its own nodes and edges" \
      "$(gunzip -c "$OUT/mx_d/bubble_1/graph_features.tsv.gz" | tail -n +2 | cut -f4 | grep -cvE '^(a|b|v1|v2)$|^(a|v1|v2)[+-]>(v1|v2|b)[+-]$')" "0"
check "describe.index.tsv names its files relative to the output directory" \
      "$(awk -F'\t' 'NR==2{print $7" "$18}' "$OUT/mx_d/describe.index.tsv")" \
      "bubble_1/kmer_features.tsv.gz bubble_1/graph_features.tsv.gz"

"$BIN" describe -i "$OUT/px.gfa" -c "$OUT/px.bubbles.csv" -o "$OUT/px_d" --feature-mode all --no-wide-matrix \
       --threads 1 -q >/dev/null 2>&1
check "a k-mer discriminative in two bubbles whose pooled row is constant is not written" \
      "$(gunzip -c "$OUT/px_d/haplotype/kmers/bimbam_kmers.bimbam.gz" |
         awk -F', ' '{d=0; for (i=5;i<=NF;i++) if ($i!=$4) d=1; if (!d) n++} END{print n+0}')" "0"
check "and the row is left out of the feature_annot sidecar too (rows match)" \
      "$(gunzip -c "$OUT/px_d/haplotype/kmers/feature_annot.kmers.tsv.gz" | tail -n +2 | wc -l | tr -d ' ')" \
      "$(gunzip -c "$OUT/px_d/haplotype/kmers/bimbam_kmers.bimbam.gz" | wc -l | tr -d ' ')"

# An error inside a worker is reported as an error, with no staging left behind.
"$BIN" describe -i "$OUT/gx.gfa" -c "$OUT/gx.bubbles.csv" -o "$OUT/gx_d" --threads 4 -q >/dev/null 2>&1
check "an error in a worker thread is a clean refusal (exit 2)" "$?" "2"
check "and leaves no output or staging directory" "$(ls -d "$OUT"/gx_d* 2>/dev/null | wc -l | tr -d ' ')" "0"

# The three substrates by default (the variant one from --variant-vcf); --only-* flags combine in any order.
"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/all3" \
     --variant-vcf "$OUT/ac.region.vcf" --no-wide-matrix -q >/dev/null 2>&1
check "by default all three substrates are written" "$(ls "$OUT/all3/haplotype" 2>/dev/null | sort | tr '\n' ' ')" "graph kmers variant "
"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/kv" --only-variant \
     --variant-vcf "$OUT/ac.region.vcf" --only-kmers --no-wide-matrix -q >/dev/null 2>&1
check "--only-variant --only-kmers writes those two" "$(ls "$OUT/kv/haplotype" 2>/dev/null | sort | tr '\n' ' ')" "kmers variant "
check "and the same variant matrix as the full run" \
      "$(gunzip -c "$OUT/kv/haplotype/variant/bimbam_variant.bimbam.gz" | cksum)" \
      "$(gunzip -c "$OUT/all3/haplotype/variant/bimbam_variant.bimbam.gz" | cksum)"
refused "--only-variant with the graph substrates but no --variant-vcf is refused" "$OUT/kv2" \
        "$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/kv2" \
        --only-kmers --only-variant -q
"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix "$OUT/ab" --out-dir "$OUT/nov" --no-wide-matrix \
     > /dev/null 2> "$OUT/nov.err"
check "without --variant-vcf the default run says the variant matrix was skipped" \
      "$(grep -c 'variant matrix: skipped' "$OUT/nov.err")" "1"
refused "--variant-vcf with --only-kmers is refused (no variant substrate is written)" "$OUT/inert1" \
        "$BIN" describe -i "$OUT/mx.gfa" -c "$OUT/mx.bubbles.csv" -o "$OUT/inert1" --only-kmers \
        --variant-vcf "$OUT/ac.region.vcf" -q
printf 'variant_id\tbubble_id\tsvtype\tnode_ids\nx\t1\tDEL\tv1\n' > "$OUT/one.vn"
refused "a node table from another call run is refused under --only-variant too" "$OUT/inert2" \
        "$BIN" describe --only-variant --variant-vcf "$OUT/ac.region.vcf" --variant-nodes "$OUT/one.vn" -o "$OUT/inert2" -q
printf 'variant_id\tbubble_id\tsvtype\tnode_ids\n' > "$OUT/none.vn"
refused "a --variant-nodes file listing no variant is refused" "$OUT/inert3" \
        "$BIN" describe -i "$OUT/mx.gfa" -c "$OUT/mx.bubbles.csv" --variant-nodes "$OUT/none.vn" -o "$OUT/inert3" -q

# A cohort simulated from a panel draws each sample's haplotypes with replacement, so one
# haplotype genotypes several samples: each counts it.
printf 'sample\thap1\thap2\nS1\th1#1#c\th2#1#c\nS2\th2#1#c\th3#1#c\n' > "$OUT/two.tsv"
"$BIN" describe -i "$OUT/mx.gfa" -c "$OUT/mx.bubbles.csv" --samples "$OUT/two.tsv" -o "$OUT/two_d" --no-wide-matrix -q >/dev/null 2>&1
check "one haplotype in two samples counts for each (exit, S1 v1 = h1+h2, S2 v1 = h2+h3)" \
      "$?:$(dosage "$OUT/two_d/sample/graph/bimbam_graph.bimbam.gz" v1 1):$(dosage "$OUT/two_d/sample/graph/bimbam_graph.bimbam.gz" v1 2)" \
      "0:$(( $(dosage "$MG" v1 "$(colof "$OUT/mx_d/haplotype/graph/samples.txt.gz" 'h1#1#c')") + $(dosage "$MG" v1 "$(colof "$OUT/mx_d/haplotype/graph/samples.txt.gz" 'h2#1#c')") )):$(( $(dosage "$MG" v1 "$(colof "$OUT/mx_d/haplotype/graph/samples.txt.gz" 'h2#1#c')") + $(dosage "$MG" v1 "$(colof "$OUT/mx_d/haplotype/graph/samples.txt.gz" 'h3#1#c')") ))"
printf 'sample\thap1\thap2\nS1\th1#1#c\th1#1#c\nS2\th2#1#c\th3#1#c\n' > "$OUT/hom.tsv"
"$BIN" describe -i "$OUT/mx.gfa" -c "$OUT/mx.bubbles.csv" --samples "$OUT/hom.tsv" -o "$OUT/hom_d" \
       --no-wide-matrix -q >/dev/null 2>&1
check "one haplotype listed twice in one sample still counts twice (homozygous; S1 v1)" \
      "$(dosage "$OUT/hom_d/sample/graph/bimbam_graph.bimbam.gz" v1 1)" "2"

{ printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\thA\thB\thC\n'
  printf 'c\t10\tcnv1\tA\t<CNV>\t.\t.\tSVTYPE=CNV;BUBBLE_ID=1\tGT:CN\t1:3\t1:3\t.:.\n'
  printf 'c\t20\tcnv2\tA\t<CNV>\t.\t.\tSVTYPE=CNV;BUBBLE_ID=2\tGT:CN\t0:2\t1:5\t1:8\n'; } > "$OUT/sc.vcf"
"$BIN" describe --only-variant --variant-vcf "$OUT/sc.vcf" -o "$OUT/sc_d" --scale-dosage -q >/dev/null 2>&1
check "--scale-dosage puts a feature with one dosage throughout at 0, and a varying one on 0..2" \
      "$(gunzip -c "$OUT/sc_d/haplotype/variant/bimbam_variant.bimbam.gz" | tr '\n' ';')" \
      "cnv1, A, B, 0, 0, NA;cnv2, A, B, 0, 1, 2;"

# ---- crossings: every crossing of a bubble, each read in the bubble's orientation ----------------------
# Site 1+..4+ with interior 2 and 3. hF crosses forward through 2, then backward through 3 (an inverted
# second copy): two crossings, each with its own node. hJ leaves through the sink straight into the source
# again: the 4+>1+ edge between its two crossings is no allele's. hH deletes the interior twice; under
# --variant-nodes every crossing keeps its boundaries, so both 1+>4+ edges count. hX only passes the boundaries' outer sides, which is no crossing.
"${PYTHON:-python3}" - "$OUT/cr" <<'PYEOF'
import sys, random
pre = sys.argv[1]; random.seed(5)
seqs = {n: ''.join(random.choice('ACGT') for _ in range(40)) for n in '012345'}
paths = {'hA#1#c': '0+,1+,2+,4+,5+', 'hB#1#c': '0+,1+,3+,4+,5+', 'hD#1#c': '0+,1+,4+,5+',
         'hF#1#c': '0+,1+,2+,4+,5+,4-,3-,1-,0-', 'hH#1#c': '0+,1+,4+,5+,1+,4+,5+',
         'hJ#1#c': '0+,1+,2+,4+,1+,3+,4+,5+', 'hX#1#c': '2+,4+,1+,3+'}
edges = set()
for p in paths.values():
    st = p.split(',')
    for x, y in zip(st, st[1:]): edges.add((x[:-1], x[-1], y[:-1], y[-1]))
with open(pre + '.gfa', 'w') as f:
    f.write('H\tVN:Z:1.0\n')
    for n, q in seqs.items(): f.write(f'S\t{n}\t{q}\n')
    for e in sorted(edges): f.write('L\t%s\t%s\t%s\t%s\t0M\n' % e)
    for name, p in paths.items(): f.write(f'P\t{name}\t{p}\t*\n')
with open(pre + '.bubbles.csv', 'w') as f:
    f.write('bubble_id,source,source_orient,sink,sink_orient,inside_node_count,total_node_count,path_support,'
            'distinct_alleles,ref_allele_support,alt_allele_support_max,alt_allele_support_min,min_inside_bp,'
            'max_inside_bp,inside_nodes\n1,1,+,4,+,2,4,7,3,1,1,1,0,40,"2;3"\n')
PYEOF
printf 'variant_id\tbubble_id\tsvtype\tnode_ids\nv1\t1\tDEL\t2\nv2\t1\tINS\t3\n' > "$OUT/cr.vn"
gm() {   # gm <describe dir> <path> <feature label>: that path's per-bubble graph count, '-' if no such feature
  "${PYTHON:-python3}" - "$1/bubble_1" "$2" "$3" <<'PYEOF'
import csv, gzip, sys
d, path, lab = sys.argv[1:]
feat = {r['label']: r['feature_name'] for r in csv.DictReader(gzip.open(d + '/graph_features.tsv.gz', 'rt'), delimiter='\t')}
rows = {r['path_name']: r for r in csv.DictReader(gzip.open(d + '/graph_matrix.tsv.gz', 'rt'), delimiter='\t')}
print(rows[path][feat[lab]] if lab in feat and path in rows else '-')
PYEOF
}
"$BIN" describe -i "$OUT/cr.gfa" -c "$OUT/cr.bubbles.csv" -o "$OUT/cr_d" --min-paths 0 -k 15 -q >/dev/null 2>&1
"$BIN" describe -i "$OUT/cr.gfa" -c "$OUT/cr.bubbles.csv" -o "$OUT/cr_v" --min-paths 0 -k 15 --variant-nodes "$OUT/cr.vn" -q >/dev/null 2>&1
check "a forward and a backward crossing are both counted (hF node 2:node 3)" \
      "$(gm "$OUT/cr_d" 'hF#1#c' 2):$(gm "$OUT/cr_d" 'hF#1#c' 3)" "1:1"
check "no edge joins one crossing to the next (4+>1+ is no feature)" "$(gm "$OUT/cr_d" 'hJ#1#c' '4+>1+')" "-"
check "a double deletion keeps both junction edges, with and without --variant-nodes (hH 1+>4+)" \
      "$(gm "$OUT/cr_d" 'hH#1#c' '1+>4+'):$(gm "$OUT/cr_v" 'hH#1#c' '1+>4+')" "2:2"
check "a path past the boundaries' outer sides is no traverser (hX)" \
      "$(gunzip -c "$OUT/cr_d/bubble_1/graph_matrix.tsv.gz" | cut -f4 | grep -c '^hX#1#c$' | tr -d ' ')" "0"

# Nested sites: the inner one's nodes would be counted once per site.
{ head -1 "$OUT/cr.bubbles.csv"; tail -n +2 "$OUT/cr.bubbles.csv"; printf '2,0,+,5,+,4,6,7,3,1,1,1,0,120,"1;2;3;4"\n'; } > "$OUT/nest.csv"
refused "nested sites are refused" "$OUT/nest_d" \
        "$BIN" describe -i "$OUT/cr.gfa" -c "$OUT/nest.csv" -o "$OUT/nest_d" -q
"$BIN" describe -i "$OUT/cr.gfa" -c "$OUT/nest.csv" --bubble-id 1 -o "$OUT/nest_1" -q >/dev/null 2>&1
check "selecting one of them still works" "$?" "0"

# A haplotype assembled in two contigs, one per bubble, is one haplotype: at each bubble it is the contig
# that crosses there, so the sample is observed at both.
printf 'sample\thap1\thap2\nS\th1#1#c\th2#1#c\nT\th3#1#c,h5#1#c\th2#1#c\n' > "$OUT/ctg.tsv"
"${PYTHON:-python3}" - "$OUT/mx.gfa" "$OUT/ctg.gfa" <<'PYEOF'
import sys
lines = open(sys.argv[1]).read().splitlines()
out = [l for l in lines if not (l.startswith('P\th5#1#c') or l.startswith('P\th3#1#c'))]
out.append('P\th3#1#c\tf0+,a+,v2+,b+,s+\t*')   # bubble 1 only
out.append('P\th5#1#c\ts+,d+,w2+,c+,f9+\t*')   # bubble 2 only
open(sys.argv[2], 'w').write('\n'.join(out) + '\n')
PYEOF
"$BIN" describe -i "$OUT/ctg.gfa" -c "$OUT/mx.bubbles.csv" --samples "$OUT/ctg.tsv" -o "$OUT/ctg_d" --no-wide-matrix -q >/dev/null 2>&1
check "a haplotype in two contigs genotypes its sample (T: v2, w2, no NA)" \
      "$(dosage "$OUT/ctg_d/sample/graph/bimbam_graph.bimbam.gz" v2 2):$(dosage "$OUT/ctg_d/sample/graph/bimbam_graph.bimbam.gz" w2 2):$(gunzip -c "$OUT/ctg_d/sample/graph/bimbam_graph.bimbam.gz" | awk -F', ' '$5=="NA"' | wc -l | tr -d ' ')" "2:2:0"

# Options whose output the run does not produce, and the variant substrate under --bubble-id.
refused "--syncmer-s with --feature-mode all is refused" "$OUT/opt1" \
        "$BIN" describe -i "$OUT/mx.gfa" -c "$OUT/mx.bubbles.csv" --feature-mode all --syncmer-s 5 -o "$OUT/opt1" -q
refused "--scale-dosage with --no-bimbam is refused" "$OUT/opt2" \
        "$BIN" describe -i "$OUT/mx.gfa" -c "$OUT/mx.bubbles.csv" --scale-dosage --no-bimbam -o "$OUT/opt2" -q
refused "--bubble-id with no row in --variant-nodes is refused" "$OUT/opt3" \
        "$BIN" describe -i "$OUT/cr.gfa" -c "$OUT/nest.csv" --bubble-id 2 --variant-nodes "$OUT/cr.vn" -o "$OUT/opt3" -q
"$BIN" describe --only-variant --variant-vcf "$OUT/sc.vcf" --bubble-id 2 -o "$OUT/bv_d" -q >/dev/null 2>&1
check "--bubble-id restricts the variant substrate" \
      "$(gunzip -c "$OUT/bv_d/haplotype/variant/bimbam_variant.bimbam.gz" | cut -d, -f1 | tr '\n' ' ')" "cnv2 "

# A record's nodes in the variant annotation come from variant_nodes.tsv when it is given: EVENT_NODES of a
# copy-number record is its boundary alone.
{ printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\thA#1#c\thB#1#c\thD#1#c\n'
  printf 'c\t10\tv1\tA\t<DEL>\t.\t.\tSVTYPE=DEL;BUBBLE_ID=1;EVENT_NODES=1\tGT\t0\t0\t1\n'
  printf 'c\t10\tv2\tA\t<INS>\t.\t.\tSVTYPE=INS;BUBBLE_ID=1;EVENT_NODES=1\tGT\t0\t1\t0\n'; } > "$OUT/cr.vcf"
"$BIN" describe -i "$OUT/cr.gfa" -c "$OUT/cr.bubbles.csv" --variant-vcf "$OUT/cr.vcf" --variant-nodes "$OUT/cr.vn" \
  -o "$OUT/cn_d" --min-paths 0 -k 15 -q >/dev/null 2>&1
check "the variant annotation's nodes are the record's variant_nodes.tsv row" \
      "$(gunzip -c "$OUT/cn_d/haplotype/variant/feature_annot.variant.tsv.gz" | awk -F'\t' 'NR>1 && $1=="v1" {print $5}')" "2"

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
"$BIN" describe -i "$INVNEST/b.sorted.gfa" -b "$INVNEST/b" -o "$INVNEST/d" -q >/dev/null 2>&1
# Described on its own nodes by oriented edges; canonical k-mers do not see an inversion between inverted repeats.
check "an inversion site has no k-mer substrate" \
      "$(awk -F'\t' '$1 == 2 { print $11 }' "$INVNEST/d/describe.index.tsv")" "skipped:inversion-site"
check "and its carriers' oriented edges are graph features (2+>7-, carried by hI and hJ)" \
      "$(gunzip -c "$INVNEST/d/bubble_2/graph_features.tsv.gz" | awk -F'\t' '$4 == "2+>7-" { print $5 }')" "2"

# ---- pooling, refusals and the variant matrix.
SPECGFA="$(cd "$(dirname "$0")/.." && pwd)/tools/spec_gfa.py"
# A k-mer two sites share, at sites no haplotype crosses together: its pooled row would be all NA.
FX1="$OUT/fx1"; mkdir -p "$FX1"
cat > "$FX1.spec" <<'SPEC'
N f1 200 1
N a 120 5
N f2 200 2
N f3 200 3
N b @a:0:1
N f4 200 4
P p1 f1+,a+,f2+
P p2 f1+,a+,f2+
P p3 f1+,a+,f2+
P p4 f1+,f2+
P p5 f1+,f2+
P p6 f1+,f2+
P q1 f3+,b+,f4+
P q2 f3+,b+,f4+
P q3 f3+,b+,f4+
P q4 f3+,f4+
P q5 f3+,f4+
P q6 f3+,f4+
B f1 f2 a
B f3 f4 b
SPEC
"${PYTHON:-python3}" -I "$SPECGFA" "$FX1.spec" "$FX1" >/dev/null
"$BIN" describe -i "$FX1/g.gfa" -c "$FX1/g.bubbles.csv" -o "$FX1/d" -q >/dev/null 2>&1
check "a pooled k-mer no haplotype observes at all its sites is not written (all-NA rows)" \
      "$(gunzip -c "$FX1/d/haplotype/kmers/bimbam_kmers.bimbam.gz" | awk -F', ' '{na=0; for (i=4;i<=NF;i++) na+=($i=="NA"); if (na==NF-3) n++} END{print n+0}')" "0"

# Two adjacent sites share a boundary (s2): its node and the k-mers inside it are counted at one site.
FX2="$OUT/fx2"; mkdir -p "$FX2"
cat > "$FX2.spec" <<'SPEC'
N s1 200 1
N a 120 5
N s2 200 2
N b 120 6
N s3 200 3
P n1 s1+,a+,s2+,b+,s3+
P n2 s1+,a+,s2+,b+,s3+
P n3 s1+,a+,s2+,b+,s3+
P d1 s1+,a+,s2+,b+,s3+,s1+,a+,s2+,b+,s3+
P d2 s1+,a+,s2+,b+,s3+,s1+,a+,s2+,b+,s3+
P d3 s1+,a+,s2+,b+,s3+,s1+,a+,s2+,b+,s3+
B s1 s2 a
B s2 s3 b
SPEC
"${PYTHON:-python3}" -I "$SPECGFA" "$FX2.spec" "$FX2" >/dev/null
"$BIN" describe -i "$FX2/g.gfa" -c "$FX2/g.bubbles.csv" -o "$FX2/d" -q >/dev/null 2>&1
check "a boundary two sites share reads its copy count (s2 on n1:d1)" \
      "$(dosage "$FX2/d/haplotype/graph/bimbam_graph.bimbam.gz" s2 1):$(dosage "$FX2/d/haplotype/graph/bimbam_graph.bimbam.gz" s2 4)" "1:2"
check "and so do the k-mers inside it (distinct k-mer rows)" \
      "$(gunzip -c "$FX2/d/haplotype/kmers/bimbam_kmers.bimbam.gz" | cut -d, -f4- | sort -u | tr -d ' ')" "1,1,1,2,2,2"

# -r places each row on the reference: node a lies at 201-320 of n1 (after s1's 200 bases); without -r, '.'.
"$BIN" describe -i "$FX2/g.gfa" -c "$FX2/g.bubbles.csv" -r n1 -o "$FX2/dr" -q >/dev/null 2>&1
check "with -r a node row gets its place on the reference (pos-end), without it '.'" \
      "$(gunzip -c "$FX2/dr/haplotype/graph/feature_annot.graph.tsv.gz" | awk -F'\t' '$1=="a"{print $7"-"$8}'):$(gunzip -c "$FX2/d/haplotype/graph/feature_annot.graph.tsv.gz" | awk -F'\t' '$1=="a"{print $7}')" "201-320:."

# Only one of the two sites repeated: L repeats site 1 (the owner of s2), R site 2. Every L and R path walks
# s2 twice, so s2 and the k-mers inside it read 2 for both, whichever site's crossings make the visit.
FX3="$OUT/fx3"; mkdir -p "$FX3"
cat > "$FX3.spec" <<'SPEC'
N s1 200 1
N a 120 5
N a2 @a:3:9
N s2 200 2
N b 120 6
N b2 @b:3:10
N s3 200 3
P n1 s1+,a+,s2+,b+,s3+
P n2 s1+,a+,s2+,b+,s3+
P L1 s1+,a+,s2+,s1+,a2+,s2+,b+,s3+
P L2 s1+,a+,s2+,s1+,a2+,s2+,b+,s3+
P R1 s1+,a+,s2+,b+,s3+,s2+,b2+,s3+
P R2 s3-,b2-,s2-,s3-,b-,s2-,a-,s1-
B s1 s2 a;a2
B s2 s3 b;b2
SPEC
"${PYTHON:-python3}" -I "$SPECGFA" "$FX3.spec" "$FX3" >/dev/null
"$BIN" describe -i "$FX3/g.gfa" -c "$FX3/g.bubbles.csv" -o "$FX3/d" -k 15 --feature-mode all -q >/dev/null 2>&1
check "a shared boundary reads every visit when only the other site repeats (s2 on n1 L1 L2 R1 R2)" \
      "$(gunzip -c "$FX3/d/haplotype/graph/bimbam_graph.bimbam.gz" | awk -F', ' '$1=="s2"{print $4,$5,$6,$7,$8,$9}')" "1 1 2 2 2 2"
# The k-mers on s2 alone in a run: how many are annotated, how many of those have a matrix row, and how many
# of the rows do not read `want` on every path.
s2_kmers() {
    gunzip -c "$1/haplotype/kmers/feature_annot.kmers.tsv.gz" | awk -F'\t' 'NR>1 && $5=="s2"{print $1}' > "$1/s2.ids"
    printf '%s:' "$(wc -l < "$1/s2.ids" | tr -d ' ')"
    gunzip -c "$1/haplotype/kmers/bimbam_kmers.bimbam.gz" |
        awk -F', ' -v want="$2" 'NR==FNR{k[$1]=1; next} ($1 in k){f++; v=$4; for(i=5;i<=NF;i++) v=v" "$i; if(v!=want) n++} END{print f+0 ":" n+0}' "$1/s2.ids" -
}
check "and so do the k-mers inside it (k-mers on s2 alone: annotated, in the matrix, not reading 1 1 2 2 2 2)" \
      "$(s2_kmers "$FX3/d" "1 1 2 2 2 2")" "186:186:0"
# With --variant-nodes and a flank of at least k, the two sites keep different ends of s2 (each its own
# interior's side), yet every k-mer kept inside it still reads each visit, both ends on every walk.
printf 'variant_id\tbubble_id\tsvtype\tnode_ids\nbubble1_INS_a2\t1\tINS\ta2\nbubble2_INS_b2\t2\tINS\tb2\n' > "$FX3/vn.tsv"
"$BIN" describe -i "$FX3/g.gfa" -c "$FX3/g.bubbles.csv" -o "$FX3/v" -k 15 --feature-mode all -q \
    --variant-nodes "$FX3/vn.tsv" --variant-flank-bp 40 >/dev/null 2>&1
check "--variant-nodes: the k-mers kept on a shared boundary read every visit (annotated, in the matrix, not reading 1 1 2 2 2 2)" \
      "$(s2_kmers "$FX3/v" "1 1 2 2 2 2")" "52:52:0"
# Without L, the owner (site 1) walks s2 once on every path, a constant it does not keep at its own site;
# R's extra visits are counted at site 2. The pooled rows still sum both, so s2 reads 1 1 2 2, not 0 0 1 1.
mkdir -p "$FX3/noL"
grep -v '^P L' "$FX3.spec" > "$FX3/noL.spec"
"${PYTHON:-python3}" -I "$SPECGFA" "$FX3/noL.spec" "$FX3/noL" >/dev/null
"$BIN" describe -i "$FX3/noL/g.gfa" -c "$FX3/noL/g.bubbles.csv" -o "$FX3/noL/d" -k 15 --feature-mode all -q >/dev/null 2>&1
check "a shared boundary constant at its owner still reads every visit (s2 on n1 n2 R1 R2)" \
      "$(gunzip -c "$FX3/noL/d/haplotype/graph/bimbam_graph.bimbam.gz" | awk -F', ' '$1=="s2"{print $4,$5,$6,$7}')" "1 1 2 2"
check "and so do the k-mers inside it (annotated, in the matrix, not reading 1 1 2 2)" \
      "$(s2_kmers "$FX3/noL/d" "1 1 2 2")" "186:186:0"

# --variant-nodes on a deletion of part of the interior: c is walked by everyone, v only by r1-r3. The
# carriers' junction c+>t+ and its k-mers are kept, and c, beside the deleted node, is read on every walk.
FX4="$OUT/fx4"; mkdir -p "$FX4"
cat > "$FX4.spec" <<'SPEC'
N s 200 1
N c 60 2
N v 150 3
N t 200 4
P r1 s+,c+,v+,t+
P r2 s+,c+,v+,t+
P r3 t-,v-,c-,s-
P d1 s+,c+,t+
P d2 s+,c+,t+
P d3 t-,c-,s-
B s t c;v
SPEC
"${PYTHON:-python3}" -I "$SPECGFA" "$FX4.spec" "$FX4" >/dev/null
printf 'variant_id\tbubble_id\tsvtype\tnode_ids\nbubble1_DEL_v\t1\tDEL\tv\n' > "$FX4/vn.tsv"
"$BIN" describe -i "$FX4/g.gfa" -c "$FX4/g.bubbles.csv" -o "$FX4/d" -k 15 --variant-nodes "$FX4/vn.tsv" -q >/dev/null 2>&1
check "--variant-nodes keeps an inner deletion's junction (c+>t+ on r1 d1)" \
      "$(dosage "$FX4/d/haplotype/graph/bimbam_graph.bimbam.gz" c+\>t+ 1):$(dosage "$FX4/d/haplotype/graph/bimbam_graph.bimbam.gz" c+\>t+ 4)" "0:1"
check "and no row for c or s+>c+, which every path walks once" \
      "$(gunzip -c "$FX4/d/haplotype/graph/bimbam_graph.bimbam.gz" | cut -d, -f1 | grep -c '^c$\|^s+>c+$')" "0"
check "and k-mers only the carriers have" \
      "$(gunzip -c "$FX4/d/haplotype/kmers/bimbam_kmers.bimbam.gz" | cut -d, -f4- | tr -d ' ' | grep -c '^0,0,0,1,1,1$' | awk '{print ($1>0)?"yes":"no"}')" "yes"
# The node before the junction (c, 5 bp) is shorter than the flank, and the one before it (a2) only the
# carriers walk, so no walk has a2 beside v: the flank is kept around the carriers' junction itself.
cat > "$FX4.b.spec" <<'SPEC'
N s 200 1
N a1 40 5
N a2 @a1:4:9
N c 5 2
N v 150 3
N t 200 4
P r1 s+,a1+,c+,v+,t+
P r2 s+,a1+,c+,v+,t+
P r3 t-,v-,c-,a1-,s-
P d1 s+,a2+,c+,t+
P d2 s+,a2+,c+,t+
P d3 t-,c-,a2-,s-
B s t a1;a2;c;v
SPEC
"${PYTHON:-python3}" -I "$SPECGFA" "$FX4.b.spec" "$FX4/b" >/dev/null
"$BIN" describe -i "$FX4/b/g.gfa" -c "$FX4/b/g.bubbles.csv" -o "$FX4/bd" -k 15 --feature-mode all \
   --variant-nodes "$FX4/vn.tsv" -q >/dev/null 2>&1
"$BIN" describe -i "$FX4/b/g.gfa" -c "$FX4/b/g.bubbles.csv" -o "$FX4/bf" -k 15 --feature-mode all -q >/dev/null 2>&1
jk() { gunzip -c "$1/haplotype/kmers/feature_annot.kmers.tsv.gz" | awk -F'\t' '$5=="a2;c;t"' | wc -l | tr -d ' '; }
check "the junction k-mers run back into a node only the carriers walk (a2;c;t k-mers, as without the option)" \
      "$(jk "$FX4/bd")" "$(jk "$FX4/bf")"

# The output directory: only describe's own files are replaced or removed; another file in a folder named
# like describe's stays, and so does its folder.
mkdir -p "$OUT/foreign/sample" "$OUT/foreign/bubble_2"
echo keep > "$OUT/foreign/sample/notes.txt"; echo keep > "$OUT/foreign/bubble_2/notes.txt"
"$BIN" describe -i "$FX2/g.gfa" -c "$FX2/g.bubbles.csv" -o "$OUT/foreign" --bubble-id 1 -q >/dev/null 2>&1
check "a run beside other files in describe-named folders succeeds" "$(verdict $?)" "accepted"
check "and leaves those files in place" \
      "$(cat "$OUT/foreign/sample/notes.txt" "$OUT/foreign/bubble_2/notes.txt" 2>/dev/null | tr '\n' ' ')" "keep keep "
"$BIN" describe -i "$FX2/g.gfa" -c "$FX2/g.bubbles.csv" -o "$FX2/d" --bubble-id 1 -q >/dev/null 2>&1
check "a previous run's folders are replaced (rerun, one site)" "$(ls -d "$FX2"/d/bubble_* | wc -l | tr -d ' ')" "1"

# Inputs refused: a site table without boundary orientation, and a P line ending in an empty step.
cut -d, -f1,2,4,6- "$FX2/g.bubbles.csv" > "$FX2/noorient.csv"
refused "a site table without orientation columns is refused" "$OUT/no1" \
        "$BIN" describe -i "$FX2/g.gfa" -c "$FX2/noorient.csv" -o "$OUT/no1" -q
awk -F'\t' -v OFS='\t' '$1 == "P" && $2 == "n1" { $3 = $3 "," } 1' "$FX2/g.gfa" > "$FX2/trail.gfa"
refused "a P line ending in an empty step is refused" "$OUT/no2" \
        "$BIN" describe -i "$FX2/trail.gfa" -c "$FX2/g.bubbles.csv" -o "$OUT/no2" -q

# Options: k=1 needs --feature-mode all; a flank of 0 is refused without --variant-nodes too.
"$BIN" describe -i "$FX2/g.gfa" -c "$FX2/g.bubbles.csv" -o "$OUT/k1" -k 1 -q > "$OUT/k1.err" 2>&1
check "-k 1 in syncmer mode names the mode" "$(grep -c 'feature-mode syncmer needs' "$OUT/k1.err")" "1"
"$BIN" describe -i "$FX2/g.gfa" -c "$FX2/g.bubbles.csv" -o "$OUT/k1a" -k 1 --feature-mode all -q >/dev/null 2>&1
check "-k 1 with --feature-mode all runs" "$?" "0"
refused "--variant-flank-bp 0 without --variant-nodes is refused" "$OUT/fl0" \
        "$BIN" describe -i "$FX2/g.gfa" -c "$FX2/g.bubbles.csv" -o "$OUT/fl0" --variant-flank-bp 0 -q

# The variant substrate on call's current FORMAT: CN on CNV records, '.' elsewhere, CNBP/CNRESID, an INV with
# non-crossers '.', and GENES percent-encoded.
vhead() { printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\th1\th2\th3\n'; }
{ vhead
  printf 'c\t10\tc1\tA\t<CNV>\t.\t.\tSVTYPE=CNV;BUBBLE_ID=1;GENES=C4A,HLA%%3AX\tGT:CN:CNBP:CNRESID\t0:2:0:0\t1:3:900:0\t1:1:-900:0\n'
  printf 'c\t10\td1\tA\t<DEL>\t.\t.\tSVTYPE=DEL;BUBBLE_ID=1\tGT:CN\t0:.\t1:.\t0:.\n'
  printf 'c\t10\ti1\tA\t<INV>\t.\t.\tSVTYPE=INV;BUBBLE_ID=2;PATH_INVERSION\tGT:CN\t0:.\t1:.\t.:.\n'; } > "$OUT/cur.vcf"
"$BIN" describe --only-variant --variant-vcf "$OUT/cur.vcf" -o "$OUT/cur_d" -q >/dev/null 2>&1
check "call's current FORMAT reads as CN, GT and NA (c1 d1 i1)" \
      "$(gunzip -c "$OUT/cur_d/haplotype/variant/bimbam_variant.bimbam.gz" | cut -d, -f4- | tr -d ' ' | tr '\n' ' ')" "2,3,1 0,1,0 0,1,NA "
check "GENES is decoded and joined by ';'" \
      "$(gunzip -c "$OUT/cur_d/haplotype/variant/feature_annot.variant.tsv.gz" | awk -F'\t' '$1=="c1" {print $7}')" "C4A;HLA:X"
{ vhead; printf 'c\t10\tm1\tA\t<CNV>,<CNV>\t.\t.\tSVTYPE=CNV;BUBBLE_ID=1;NALLELES=3\tGT:CN\t0:2\t1:3\t2:4\n'; } > "$OUT/mcnv.vcf"
"$BIN" describe --only-variant --variant-vcf "$OUT/mcnv.vcf" -o "$OUT/mcnv_d" -q >/dev/null 2>&1
check "a copy-number record with two ALTs is one CN row" \
      "$(gunzip -c "$OUT/mcnv_d/haplotype/variant/bimbam_variant.bimbam.gz" | tr -d ' ' | tr '\n' ' ')" "m1,A,B,2,3,4 "
{ vhead; printf 'c\t10\t.\tA\t<DEL>\t.\t.\tSVTYPE=DEL;BUBBLE_ID=1\tGT\t0\t1\t0\n'; } > "$OUT/dot.vcf"
refused "a record with no ID is refused" "$OUT/dot_d" "$BIN" describe --only-variant --variant-vcf "$OUT/dot.vcf" -o "$OUT/dot_d" -q
{ vhead; printf 'c\t10\tm2\tA\t<DEL>,<INS>\t.\t.\tSVTYPE=DEL;BUBBLE_ID=1;NALLELES=3;AF=0.4\tGT\t0\t1\t2\n'; } > "$OUT/af.vcf"
refused "one AF value on a record with two ALTs is refused" "$OUT/af_d" \
        "$BIN" describe --only-variant --variant-vcf "$OUT/af.vcf" -o "$OUT/af_d" -q

# The VCF and the node table must come from one call run.
printf 'variant_id\tbubble_id\tsvtype\tnode_ids\nv1\t1\tDEL\t2\n' > "$OUT/cr1.vn"
refused "a VCF record with no node-table row is refused" "$OUT/vn1" \
        "$BIN" describe -i "$OUT/cr.gfa" -c "$OUT/cr.bubbles.csv" --variant-vcf "$OUT/cr.vcf" --variant-nodes "$OUT/cr1.vn" -o "$OUT/vn1" -k 15 -q
printf 'variant_id\tbubble_id\tsvtype\tnode_ids\nv1\t1\tINS\t2\nv2\t1\tINS\t3\n' > "$OUT/cr2.vn"
refused "a record of another type in the node table is refused" "$OUT/vn2" \
        "$BIN" describe -i "$OUT/cr.gfa" -c "$OUT/cr.bubbles.csv" --variant-vcf "$OUT/cr.vcf" --variant-nodes "$OUT/cr2.vn" -o "$OUT/vn2" -k 15 -q

# -q silences the mixed-ploidy note; the log counts only the files written (no k-mer files at an inversion site).
printf 'sA\tn1\tn2\nsB\tn3\n' > "$FX2/mixed.tsv"
check "-q silences the --samples ploidy note" \
      "$("$BIN" describe -i "$FX2/g.gfa" -c "$FX2/g.bubbles.csv" -o "$OUT/mx_q" --samples "$FX2/mixed.tsv" -q 2>&1 | wc -l | tr -d ' ')" "0"
"$BIN" describe -i "$INVNEST/b.sorted.gfa" -b "$INVNEST/b" -o "$INVNEST/d2" 2> "$INVNEST/d2.err" >/dev/null
check "the log's file count is the files written" \
      "$(sed -n 's/.*wrote \([0-9]*\) files.*/\1/p' "$INVNEST/d2.err")" "$(find "$INVNEST/d2" -type f | wc -l | tr -d ' ')"

# --only-variant reads the node table for the variant rows' nodes, as a full run does; a flank there widens
# nothing and is refused.
"$BIN" describe --only-variant --variant-vcf "$OUT/cr.vcf" --variant-nodes "$OUT/cr.vn" -o "$OUT/ov_d" -q >/dev/null 2>&1
check "--only-variant with the node table gives the full run's variant annotation" \
      "$(gunzip -c "$OUT/ov_d/haplotype/variant/feature_annot.variant.tsv.gz" | cksum)" \
      "$(gunzip -c "$OUT/cn_d/haplotype/variant/feature_annot.variant.tsv.gz" | cksum)"
refused "--variant-flank-bp with --only-variant is refused" "$OUT/ov_f" \
        "$BIN" describe --only-variant --variant-vcf "$OUT/cr.vcf" --variant-nodes "$OUT/cr.vn" --variant-flank-bp 5 -o "$OUT/ov_f" -q
# A samples table with no sample would give a sample level with no columns.
printf 'sample\thap1\thap2\n' > "$OUT/nosam.tsv"
refused "a --samples table with no sample is refused" "$OUT/nosam_d" \
        "$BIN" describe -i "$FX2/g.gfa" -c "$FX2/g.bubbles.csv" -o "$OUT/nosam_d" --samples "$OUT/nosam.tsv" -q
# A diploid VCF (call --diploid) is refused, pointing to the per-haplotype VCF and --samples.
awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print;next} {$10="0|1"; print}' "$OUT/cr.vcf" > "$OUT/dip.vcf"
check "a diploid VCF is refused with a pointer to --samples" \
      "$("$BIN" describe --only-variant --variant-vcf "$OUT/dip.vcf" -o "$OUT/dip_d" -q 2>&1 | grep -c -- '--samples')" "1"

# Two sites share a boundary b shorter than k-1: the k-mers from x1 across b into x2 are counted once, at
# site 1, which counts b; site 2 spells only b's bases. Without the bases past b, neither site has them.
"${PYTHON:-python3}" - "$OUT/sb" <<'PYEOF'
import sys, random
pre = sys.argv[1]; random.seed(11)
seq = {n: ''.join(random.choice('ACGT') for _ in range(L)) for n, L in (('a', 40), ('x1', 20), ('b', 6), ('x2', 20), ('c', 40))}
paths = {'h1': 'a x1 b x2 c', 'h2': 'a x1 b x2 c', 'h3': 'a x1 b x2 c', 'h4': 'a b x2 c', 'h5': 'a x1 b c', 'h6': 'a b c'}
edges = set()
for p in paths.values():
    st = p.split(); edges.update(zip(st, st[1:]))
with open(pre + '.gfa', 'w') as f:
    f.write('H\tVN:Z:1.0\n')
    for n, s in seq.items(): f.write('S\t%s\t%s\n' % (n, s))
    for x, y in sorted(edges): f.write('L\t%s\t+\t%s\t+\t0M\n' % (x, y))
    for n, p in paths.items(): f.write('P\t%s#1#c\t%s\t*\n' % (n, ','.join(s + '+' for s in p.split())))
with open(pre + '.bubbles.csv', 'w') as f:
    f.write('bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n')
    f.write('1,a,+,b,+,6,0,20,"x1"\n2,b,+,c,+,6,0,20,"x2"\n')
k = seq['x1'][-5:] + seq['b'] + seq['x2'][:4]
rc = k[::-1].translate(str.maketrans('ACGT', 'TGCA'))
open(pre + '.kmer', 'w').write(min(k, rc) + '\n')
PYEOF
"$BIN" describe -i "$OUT/sb.gfa" -c "$OUT/sb.bubbles.csv" -o "$OUT/sb_d" -k 15 --feature-mode all -q >/dev/null 2>&1
sbk=$(cat "$OUT/sb.kmer")
check "a k-mer across a short shared boundary is counted at one site, in its carriers (site 1:site 2)" \
      "$(gunzip -c "$OUT/sb_d/bubble_1/kmer_features.tsv.gz" | awk -F'\t' -v k="$sbk" '$4 == k { print $5 }'):$(gunzip -c "$OUT/sb_d/bubble_2/kmer_features.tsv.gz" | awk -F'\t' -v k="$sbk" '$4 == k { n++ } END { print n+0 }')" \
      "3:0"

# ---- the sample level, --scale-dosage on pooled rows, the wide-matrix cap, option and input refusals.
# Every sample below carries the boundary a twice and so reads 2 there, while h4 and h5 (in no sample)
# make it vary across the haplotypes: the haplotype level keeps it, the sample level drops it.
printf 'sample\thap1\thap2\nS1\th1#1#c\th2#1#c\nS2\th1#1#c\th3#1#c\nS3\th2#1#c\th3#1#c\n' > "$OUT/cst.tsv"
"$BIN" describe -i "$OUT/mx.gfa" -c "$OUT/mx.bubbles.csv" --samples "$OUT/cst.tsv" -o "$OUT/cst_d" --no-wide-matrix -q >/dev/null 2>&1
check "a row varying across haplotypes but constant across samples: kept at the haplotype level, dropped at the sample level" \
      "$(gunzip -c "$OUT/cst_d/haplotype/graph/bimbam_graph.bimbam.gz" | grep -c '^a,'):$(gunzip -c "$OUT/cst_d/sample/graph/bimbam_graph.bimbam.gz" | grep -c '^a,')" "1:0"
const_rows() { gunzip -c "$1" | awk -F', ' '{v=""; d=0; for (i=4;i<=NF;i++) if ($i!="NA") { if (v=="") v=$i; else if ($i!=v) d=1 } if (!d) n++} END{print n+0}'; }
check "no sample-level k-mer or graph row is constant (kmers:graph)" \
      "$(const_rows "$OUT/cst_d/sample/kmers/bimbam_kmers.bimbam.gz"):$(const_rows "$OUT/cst_d/sample/graph/bimbam_graph.bimbam.gz")" "0:0"
check "and the sample-level sidecar keeps one row per matrix row" \
      "$(gunzip -c "$OUT/cst_d/sample/kmers/feature_annot.kmers.tsv.gz" | tail -n +2 | wc -l | tr -d ' ')" \
      "$(gunzip -c "$OUT/cst_d/sample/kmers/bimbam_kmers.bimbam.gz" | wc -l | tr -d ' ')"

# --scale-dosage on the pooled k-mer and graph rows, at both levels: each row is 2 (v - lo) / (hi - lo) of
# the unscaled row, lo and hi over its observed columns only (h5 does not cross site 1 in this graph, so
# reads NA at a; counted as 0 it would move a's other paths off 0).
"$BIN" describe -i "$OUT/ctg.gfa" -c "$OUT/mx.bubbles.csv" --samples "$OUT/ctg.tsv" -o "$OUT/raw_d" --no-wide-matrix -q >/dev/null 2>&1
"$BIN" describe -i "$OUT/ctg.gfa" -c "$OUT/mx.bubbles.csv" --samples "$OUT/ctg.tsv" -o "$OUT/scl_d" --no-wide-matrix --scale-dosage -q >/dev/null 2>&1
scaled_ok() { # scaled_ok <raw.gz> <scaled.gz>: rows that differ from the min-max rescaling of the raw row
  paste -d'|' <(gunzip -c "$1") <(gunzip -c "$2") | awk -F'|' '{
    nr = split($1, r, ", "); ns = split($2, s, ", "); bad = (nr != ns || r[1] != s[1]); lo = ""; hi = ""
    for (i = 4; i <= nr; i++) if (r[i] != "NA") { if (lo == "" || r[i] + 0 < lo) lo = r[i] + 0; if (hi == "" || r[i] + 0 > hi) hi = r[i] + 0 }
    for (i = 4; i <= nr && !bad; i++) {
      if (r[i] == "NA") { bad = (s[i] != "NA"); continue }
      e = (hi > lo) ? 2 * (r[i] - lo) / (hi - lo) : 0; d = s[i] - e; if (d < 0) d = -d; if (d > 1e-4) bad = 1 }
    n += bad; rows++ } END { print (rows > 0 ? n : "empty") }'; }
check "--scale-dosage on pooled rows is the min-max of the observed columns (hap kmers:graph, sample kmers:graph)" \
      "$(scaled_ok "$OUT/raw_d/haplotype/kmers/bimbam_kmers.bimbam.gz" "$OUT/scl_d/haplotype/kmers/bimbam_kmers.bimbam.gz"):$(scaled_ok "$OUT/raw_d/haplotype/graph/bimbam_graph.bimbam.gz" "$OUT/scl_d/haplotype/graph/bimbam_graph.bimbam.gz"),$(scaled_ok "$OUT/raw_d/sample/kmers/bimbam_kmers.bimbam.gz" "$OUT/scl_d/sample/kmers/bimbam_kmers.bimbam.gz"):$(scaled_ok "$OUT/raw_d/sample/graph/bimbam_graph.bimbam.gz" "$OUT/scl_d/sample/graph/bimbam_graph.bimbam.gz")" "0:0,0:0"
H4S=$(colof "$OUT/scl_d/haplotype/graph/samples.txt.gz" 'h4#1#c'); H1S=$(colof "$OUT/scl_d/haplotype/graph/samples.txt.gz" 'h1#1#c')
H5S=$(colof "$OUT/scl_d/haplotype/graph/samples.txt.gz" 'h5#1#c')
check "a scaled: h4 (2 copies) at 2, h1 (1 copy) at 0, h5 (not crossing) NA" \
      "$(dosage "$OUT/scl_d/haplotype/graph/bimbam_graph.bimbam.gz" a "$H4S"):$(dosage "$OUT/scl_d/haplotype/graph/bimbam_graph.bimbam.gz" a "$H1S"):$(dosage "$OUT/scl_d/haplotype/graph/bimbam_graph.bimbam.gz" a "$H5S")" "2:0:NA"

# --max-wide-features: a site with more k-mers than the cap gets no dense matrix, unless --force-wide.
"$BIN" describe -i "$OUT/mx.gfa" -c "$OUT/mx.bubbles.csv" -o "$OUT/cap_d" --max-wide-features 1 -q >/dev/null 2>&1
check "above --max-wide-features: matrix_reason skipped:feature-cap, no matrix file (reason:path:file)" \
      "$(awk -F'\t' 'NR==2{print $11":"$8}' "$OUT/cap_d/describe.index.tsv"):$([ -e "$OUT/cap_d/bubble_1/kmer_matrix.tsv.gz" ] && echo present || echo absent)" \
      "skipped:feature-cap:.:absent"
"$BIN" describe -i "$OUT/mx.gfa" -c "$OUT/mx.bubbles.csv" -o "$OUT/capf_d" --max-wide-features 1 --force-wide -q >/dev/null 2>&1
check "with --force-wide it is written anyway" \
      "$(awk -F'\t' 'NR==2{print $11}' "$OUT/capf_d/describe.index.tsv"):$([ -s "$OUT/capf_d/bubble_1/kmer_matrix.tsv.gz" ] && echo present || echo absent)" \
      "written:present"

# An explicit --syncmer-s 0 is refused, not read as the default s.
refused "--syncmer-s 0 is refused" "$OUT/s0_d" \
        "$BIN" describe -i "$OUT/mx.gfa" -c "$OUT/mx.bubbles.csv" -o "$OUT/s0_d" --syncmer-s 0 -k 15 -q

# Under --only-variant, an unreadable or headerless --variant-nodes is named as such, not as a mismatch.
check "--only-variant with a missing --variant-nodes says it cannot open it" \
      "$("$BIN" describe --only-variant --variant-vcf "$OUT/ac.region.vcf" --variant-nodes "$OUT/no_such.vn" -o "$OUT/vnx_d" -q 2>&1 | grep -c 'failed to open variant nodes')" "1"
tail -n +2 "$OUT/vndel.tsv" > "$OUT/nohead.vn"
check "and a table without call's header is refused for its header" \
      "$("$BIN" describe --only-variant --variant-vcf "$OUT/ac.region.vcf" --variant-nodes "$OUT/nohead.vn" -o "$OUT/vnh_d" -q 2>&1 | grep -c 'header is not call')" "1"

# A site no path crosses is skipped:no-paths and gets no folder.
FX5="$OUT/fx5"; mkdir -p "$FX5"
{ cat "$FX2.spec"; printf 'N x1 50 7\nN y 40 8\nN x2 50 9\nB x1 x2 y\n'; } > "$FX5.spec"
"${PYTHON:-python3}" -I "$SPECGFA" "$FX5.spec" "$FX5" >/dev/null
"$BIN" describe -i "$FX5/g.gfa" -c "$FX5/g.bubbles.csv" -o "$FX5/d" -q >/dev/null 2>&1
check "a skipped:no-paths site gets no bubble folder (status:folder)" \
      "$(awk -F'\t' '$1==3{print $2}' "$FX5/d/describe.index.tsv"):$([ -e "$FX5/d/bubble_3" ] && echo present || echo absent)" \
      "skipped:no-paths:absent"

# A failed run removes the directories it made for its output, parents included.
mkdir -p "$OUT/fresh_parent"
"$BIN" describe --only-variant --variant-vcf "$OUT/no_such.vcf" -o "$OUT/fresh_parent/new/a/b" -q >/dev/null 2>&1
check "a failed run leaves none of the directories it made" "$(ls "$OUT/fresh_parent" | wc -l | tr -d ' ')" "0"
"$BIN" describe -i "$FX2/g.gfa" -c "$OUT/no_such.csv" -o "$OUT/fresh_parent/new/c" -q >/dev/null 2>&1
check "nor does one failing in the graph pass" "$(ls "$OUT/fresh_parent" | wc -l | tr -d ' ')" "0"

# --samples and --variant-nodes are read by both passes, so a pipe is refused for what it is.
check "--samples from a pipe, read by both passes, is refused as a pipe" \
      "$("$BIN" describe -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" --variant-vcf "$OUT/ac.region.vcf" \
         --samples <(cat "$OUT/sam.tsv") -o "$OUT/pipe_d" -q 2>&1 | grep -c 'not a regular file')" "1"
check "and so is --variant-nodes" \
      "$("$BIN" describe -i "$OUT/ab.sorted.gfa" -c "$OUT/ab.bubbles.csv" --variant-vcf "$OUT/ac.region.vcf" \
         --variant-nodes <(cat "$OUT/ac.variant_nodes.tsv") -o "$OUT/pipe_d" -q 2>&1 | grep -c 'not a regular file')" "1"

# The closing line names the directory once and its entries; -t is --threads.
"$BIN" describe -i "$FX2/g.gfa" -c "$FX2/g.bubbles.csv" -o "$OUT/log_d/" -t 2 > /dev/null 2> "$OUT/log_d.err"
check "the closing line is 'wrote N files in <dir>/: ...'" \
      "$(grep -c "\] wrote [0-9]* files in $OUT/log_d/: describe.index.tsv, haplotype/, 2 bubble_<id>/$" "$OUT/log_d.err")" "1"

printf "%d assertion(s) failed\n" "$fails"
[ "$fails" -eq 0 ]
