#!/usr/bin/env bash
# describe_stats.sh - contract assertions for `panvar describe`.
#
# The registered C4 smoke test checks that describe produces output. It asserts nothing about the
# dosages, and the module's central defect was invisible to it: all three step-building sites required
# an interior bubble node, so a path taking the direct source->sink edge -- a pure deletion, the most
# common SV class -- was not counted as a traverser at all. On a 7-path fixture where 3 haplotypes
# delete the interior the bubble reported 4 paths, and the node that discriminates them was then
# discarded as non-discriminative BECAUSE ONLY ITS CARRIERS HAD BEEN OBSERVED. The bubble emitted zero
# features, and every file still existed.
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

"$BIN" bubble -i "$OUT/a.gfa" -r "$REF" -o "$OUT/ab" -q >/dev/null 2>&1
check "fixture A presents two bubbles" "$(tail -n +2 "$OUT/ab.bubbles.csv" | wc -l | tr -d ' ')" "2"

"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix-in "$OUT/ab" --out-dir "$OUT/A" \
   --no-wide-matrix -q >/dev/null 2>&1
AI="$OUT/A/describe.index.tsv"
GB="$OUT/A/haplotype/graph/bimbam_graph.bimbam.gz"
GS="$OUT/A/haplotype/graph/samples.txt.gz"

# --- THE reproducer: a deletion allele is a traversal.
# ref + 3 full + 3 del + trunc = 8, and every one of them crosses bubble 1 (trunc ends just after it).
# Before the fix only the 5 that keep the interior node were counted.
check "bubble 1 counts all 8 paths, not just the 5 that keep the interior" \
      "$(awk -F'\t' 'NR==2{print $3}' "$AI")" "8"
check "and it yields features instead of nothing" \
      "$(awk -F'\t' 'NR==2{print ($13>0 && $16>0) ? "yes" : "no"}' "$AI")" "yes"

# --- zero vs missing, in both directions, on hand-checkable cells.
cdel=$(colof "$GS" 'del1#1#chr1'); cfull=$(colof "$GS" 'full1#1#chr1'); ctr=$(colof "$GS" 'trunc#1#chr1')
b1node=$(awk -F',' 'NR==2{gsub(/"/,"",$NF); split($NF,x,";"); print x[1]}' "$OUT/ab.bubbles.csv")
check "a deletion haplotype reads 0 at the deleted node (traverses, does not carry)" \
      "$(dosage "$GB" "$b1node" "$cdel")" "0"
check "a carrier reads 1 there" "$(dosage "$GB" "$b1node" "$cfull")" "1"
b2node=$(awk -F',' 'NR==3{gsub(/"/,"",$NF); split($NF,x,";"); print x[1]}' "$OUT/ab.bubbles.csv")
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
"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix-in "$OUT/ab" --out-dir "$OUT/T1" \
   --no-wide-matrix --threads 1 -q >/dev/null 2>&1
"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix-in "$OUT/ab" --out-dir "$OUT/T8" \
   --no-wide-matrix --threads 8 -q >/dev/null 2>&1
same=1
for f in haplotype/graph/bimbam_graph.bimbam.gz haplotype/kmers/bimbam_kmers.bimbam.gz \
         haplotype/graph/feature_annot.graph.tsv.gz haplotype/kmers/feature_annot.kmers.tsv.gz; do
  cmp -s <(gunzip -c "$OUT/T1/$f") <(gunzip -c "$OUT/T8/$f") || same=0
done
check "output is byte-identical at 1 and 8 threads" "$same" "1"

# ---------------------------------------------------------------------------------------------
# Fixture B: a SHORT variant beside VERY LONG neighbours. --variant-flank-bp must admit bases, not
# whole nodes: at 30 it used to pull in both 4000 bp neighbours entire, with every k-mer in them.
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
"$BIN" bubble -i "$OUT/b.gfa" -r "$BREF" -o "$OUT/bb" -q >/dev/null 2>&1
bid=$(tail -n +2 "$OUT/bb.bubbles.csv" | cut -d, -f1 | head -1)
ins=$(tail -n +2 "$OUT/bb.bubbles.csv" | head -1 | sed 's/.*"\(.*\)".*/\1/' | tr ';' ',')
printf 'variant_id\tbubble_id\tsvtype\tnode_ids\nv1\t%s\tINS\t%s\n' "$bid" "$ins" > "$OUT/vn.tsv"
cand() { awk -F'\t' 'NR==2{print $4}' "$1/describe.index.tsv"; }
for fl in 0 30 200; do
  "$BIN" describe -i "$OUT/bb.sorted.gfa" --bubble-prefix-in "$OUT/bb" --out-dir "$OUT/F$fl" \
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
  "$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix-in "$OUT/ab" --out-dir "$OUT/E1" --bubble-id 999 -q

awk -F, 'BEGIN{OFS=","} NR==2{$2="NOSUCHNODE"} {print}' "$OUT/ab.bubbles.csv" > "$OUT/badcsv.csv"
refused "a bubbles CSV naming a node the graph lacks is refused" "$OUT/E2" \
  "$BIN" describe -i "$OUT/ab.sorted.gfa" -c "$OUT/badcsv.csv" --out-dir "$OUT/E2" -q

# variant_nodes goes through the SHARED strict reader; a stale node silently removes real features.
printf 'variant_id\tbubble_id\tsvtype\tnode_ids\nv1\t%s\tINS\tNOSUCHNODE\n' "$bid" > "$OUT/vnbad.tsv"
refused "a variant_nodes node the bubble does not contain is refused" "$OUT/E3" \
  "$BIN" describe -i "$OUT/bb.sorted.gfa" --bubble-prefix-in "$OUT/bb" --out-dir "$OUT/E3" \
     --variant-nodes "$OUT/vnbad.tsv" -q
tail -n +2 "$OUT/vn.tsv" > "$OUT/vnnohdr.tsv"
refused "a variant_nodes file with no header is refused" "$OUT/E4" \
  "$BIN" describe -i "$OUT/bb.sorted.gfa" --bubble-prefix-in "$OUT/bb" --out-dir "$OUT/E4" \
     --variant-nodes "$OUT/vnnohdr.tsv" -q

# --samples: a duplicate sample id makes the haplotype assignment depend on file order.
printf 'sample\thap1\thap2\ns1\tfull1#1#chr1\tdel1#1#chr1\ns1\tfull2#1#chr1\tdel2#1#chr1\n' > "$OUT/dupsam.tsv"
refused "a duplicate sample id in --samples is refused" "$OUT/E5" \
  "$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix-in "$OUT/ab" --out-dir "$OUT/E5" \
     --samples "$OUT/dupsam.tsv" -q

refused "two --only-* flags are refused, not silently intersected to nothing" "$OUT/E6" \
  "$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix-in "$OUT/ab" --out-dir "$OUT/E6" \
     --only-kmers --only-graph -q
refused "--variant-flank-bp without --variant-nodes is refused as inert" "$OUT/E7" \
  "$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix-in "$OUT/ab" --out-dir "$OUT/E7" \
     --variant-flank-bp 30 -q

# ---------------------------------------------------------------------------------------------
# Transactionality: a rerun with a different substrate selection must not leave the old family
# looking current, a failing run must not disturb what is there, and neither may touch user files.
# ---------------------------------------------------------------------------------------------
TX="$OUT/tx"
"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix-in "$OUT/ab" --out-dir "$TX" --no-wide-matrix -q >/dev/null 2>&1
check "a full run writes both graph substrates" "$(ls "$TX/haplotype" | sort | tr '\n' ' ')" "graph kmers "
touch "$TX/USER_FILE.txt"
"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix-in "$OUT/ab" --out-dir "$TX" --only-graph --no-wide-matrix -q >/dev/null 2>&1
check "a --only-graph rerun removes the now-stale k-mer family" "$(ls "$TX/haplotype" | sort | tr '\n' ' ')" "graph "
check "and leaves an unrelated user file alone" "$([ -f "$TX/USER_FILE.txt" ] && echo kept || echo REMOVED)" "kept"
before=$(ls "$TX" | sort | tr '\n' ' ')
"$BIN" describe -i "$OUT/ab.sorted.gfa" --bubble-prefix-in "$OUT/ab" --out-dir "$TX" --bubble-id 999 -q >/dev/null 2>&1
check "a failing run leaves the directory untouched" "$(ls "$TX" | sort | tr '\n' ' ')" "$before"
check "and leaves no staging directory behind" "$(ls -d "$TX".describe-staging.* 2>/dev/null | wc -l | tr -d ' ')" "0"

# ---------------------------------------------------------------------------------------------
# Provenance: two runs that produce different matrices must have distinguishable params.
# ---------------------------------------------------------------------------------------------
check "params.json records which substrates were emitted" \
      "$(grep -c '"emit_kmers"\|"emit_graph"\|"emit_variant"' "$TX/describe.params.json")" "3"
check "and the variant restriction, flank and thread count" \
      "$(~/miniconda3/bin/python -c "import json;d=json.load(open('$TX/describe.params.json'));print(sum(k in d for k in ('variant_nodes','variant_flank_bp','threads','samples','scale_dosage','input_sizes_bytes')))")" "6"

printf "%d assertion(s) failed\n" "$fails"
[ "$fails" -eq 0 ]
