#!/usr/bin/env bash
# call_stats.sh - coordinate contract assertions for `panvar call`.
#
# real_smoke.sh checks that call RUNS and that its files exist; it asserts nothing about where a
# record lands. The coordinates are the part of a VCF every downstream consumer trusts silently, so
# they need a fixture whose every expected number is derivable by hand rather than read back from
# the caller's own output.
#
# The reference chain below starts at genomic 1, so a genomic coordinate IS an offset along the path:
#
#   n1 100 bp    1..100     flank
#   n2  70 bp  101..170     deleted by hapdel
#   n3 100 bp  171..270     flank
#   n4  30 bp  271..300     flank
#   n5  80 bp  301..380     inverted by hapinv
#   n6 100 bp  381..480     flank
#   n7  60 bp  481..540     flank
#   nI  60 bp     --        inserted by hapins, between n6 and n7
#
# A symbolic DEL/INV is anchored on the base PRECEDING the event, so REF carries one real base and
# the event occupies POS+1 onwards. That makes the expected records:
#
#   DEL  POS 100  END 170  SVLEN -70    REF = last base of n1
#   INV  POS 300  END 380  SVLEN  80    REF = last base of n4
#   INS  POS 480  END 480  SVLEN  60    REF = last base of n6
#
#   call_stats.sh <panvar-binary> <out-dir>
set -uo pipefail

BIN="${1:?usage: call_stats.sh <panvar> <outdir>}"
OUT="${2:?}"
mkdir -p "$OUT"
OUT="$OUT/run.$$.$(date +%s)"
rm -rf "$OUT"; mkdir -p "$OUT"
fails=0
ok()  { printf "  ok   %s\n" "$1"; }
failed=()
bad() { printf "  FAIL %s\n" "$1"; fails=$((fails + 1)); failed+=("${1:0:300}"); }

# Deterministic sequence, distinct per node ($2 seeds it): identical flanks would give the snarl finder a
# different graph than intended, and a repeated motif is a tandem repeat, which bubble extends sites to
# and call reports whole alleles at. A Park-Miller generator, exact in any awk's doubles.
blk() { awk -v n="$1" -v seed="${2:-0}" 'BEGIN{ b="ACGT"; x=(seed*7919+1)%2147483647;
          for(i=0;i<n;i++){ x=(x*16807)%2147483647; printf "%s", substr(b, int(x*4/2147483647)+1, 1) } }'; }

N1=$(blk 100 0); N2=$(blk 70 3);  N3=$(blk 100 7); N4=$(blk 30 11)
N5=$(blk 80 17); N6=$(blk 100 23); N7=$(blk 60 29); NI=$(blk 60 31)
REF_NAME='ref#0#chr1:1-540'

{ printf 'H\tVN:Z:1.0\n'
  printf "S\tn1\t%s\nS\tn2\t%s\nS\tn3\t%s\nS\tn4\t%s\n" "$N1" "$N2" "$N3" "$N4"
  printf "S\tn5\t%s\nS\tn6\t%s\nS\tn7\t%s\nS\tnI\t%s\n" "$N5" "$N6" "$N7" "$NI"
  # forward links along the reference chain, plus the deletion bypass and the insertion detour
  for e in "n1 + n2 +" "n2 + n3 +" "n1 + n3 +" "n3 + n4 +" "n4 + n5 +" "n5 + n6 +" \
           "n6 + n7 +" "n6 + nI +" "nI + n7 +" "n4 + n5 -" "n5 - n6 +"; do
    set -- $e; printf "L\t%s\t%s\t%s\t%s\t0M\n" "$1" "$2" "$3" "$4"
  done
  printf 'P\t%s\tn1+,n2+,n3+,n4+,n5+,n6+,n7+\t*\n' "$REF_NAME"
  printf 'P\thapref#1#chr1\tn1+,n2+,n3+,n4+,n5+,n6+,n7+\t*\n'
  printf 'P\thapdel#1#chr1\tn1+,n3+,n4+,n5+,n6+,n7+\t*\n'
  printf 'P\thapinv#1#chr1\tn1+,n2+,n3+,n4+,n5-,n6+,n7+\t*\n'
  printf 'P\thapins#1#chr1\tn1+,n2+,n3+,n4+,n5+,n6+,nI+,n7+\t*\n'; } > "$OUT/coord.gfa"

"$BIN" bubble -i "$OUT/coord.gfa" -r "$REF_NAME" -o "$OUT/cb" --merge-nearby-bp 0 -q >/dev/null 2>&1
nb=$(tail -n +2 "$OUT/cb.bubbles.csv" 2>/dev/null | wc -l | tr -d ' ')
[ "$nb" = "3" ] && ok "the fixture presents three bubbles (deletion, inversion, insertion)" \
                || bad "expected 3 bubbles in the fixture, got $nb"

"$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" \
      -o "$OUT/cc" -q >/dev/null 2>&1
VCF="$OUT/cc.region.vcf"
[ -s "$VCF" ] || { bad "call wrote no region VCF"; printf "%d assertion(s) failed\n" "$((fails + 1))"; exit 1; }

# field <svtype> <column>  -- POS=2 REF=4; END/SVLEN come out of INFO
field() { awk -F'\t' -v t="$1" -v c="$2" '$0!~/^#/ && $8~("SVTYPE="t";") {print $c; exit}' "$VCF"; }
info()  { awk -F'\t' -v t="$1" -v k="$2" '$0!~/^#/ && $8~("SVTYPE="t";"){
            n=split($8,a,";"); for(i=1;i<=n;i++){split(a[i],kv,"="); if(kv[1]==k){print kv[2]; exit}}}' "$VCF"; }

check() { # check <label> <got> <want>
  [ "$2" = "$3" ] && ok "$1 = $3" || bad "$1: expected $3, got ${2:-<missing>}"; }

# The event coordinates, each derivable from the node lengths at the top of this file.
check "DEL POS (base before the deleted span)"  "$(field DEL 2)"      "100"
check "DEL END (last deleted base)"             "$(info  DEL END)"    "170"
check "DEL SVLEN"                               "$(info  DEL SVLEN)"  "-70"
check "INV POS (base before the inverted span)" "$(field INV 2)"      "300"
check "INV END (last inverted base)"            "$(info  INV END)"    "380"
check "INV SVLEN"                               "$(info  INV SVLEN)"  "80"
check "INS POS (base before the insertion)"     "$(field INS 2)"      "480"
check "INS END (INS does not span reference)"   "$(info  INS END)"    "480"
check "INS SVLEN"                               "$(info  INS SVLEN)"  "60"

# REF must be the reference base AT POS, which is the last base of the preceding flank. Anchoring on
# the first affected base instead would put the first base of n2/n5 here and still look well formed,
# so this is what separates the two conventions.
check "DEL REF is the last base of n1" "$(field DEL 4)" "${N1: -1}"
check "INV REF is the last base of n4" "$(field INV 4)" "${N4: -1}"
check "INS REF is the last base of n6" "$(field INS 4)" "${N6: -1}"

check "every record is FILTER PASS, which the header defines" \
  "$(awk -F'\t' '!/^#/{print $7}' "$VCF" | sort -u | paste -sd, -):$(grep -c '^##FILTER=<ID=PASS,' "$VCF")" "PASS:1"
check "the IMPRECISE header gives its three causes (cluster, N run, PATH_INVERSION)" \
  "$(grep '^##INFO=<ID=IMPRECISE,' "$VCF" | grep -c 'CLUSTER.*an N in the reference.*PATH_INVERSION')" "1"
# Fields that repeat others are not written: START_NODE and END_NODE (EVENT_NODES), NMERGED (AC), CN_SCOPE and
# REF_CN_SOURCE (CN_METHOD), CN_REF_FOLD (REF_CN), CN_SHARED_BP (CN_UNIT_BP x REF_CN), MERGE_SIZE_RATIO
# (SVLEN_RANGE).
check "no redundant fields in the header" \
  "$(grep -cE '^##INFO=<ID=(START_NODE|END_NODE|NMERGED|CN_SCOPE|REF_CN_SOURCE|CN_REF_FOLD|CN_SHARED_BP|MERGE_SIZE_RATIO),' "$VCF")" "0"
# The header states the --max-seq-bp the sequences were written with: 60 leaves out the 70 bp DELSEQ.
"$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" -o "$OUT/cmax" --max-seq-bp 60 -q >/dev/null 2>&1
check "--max-seq-bp 60: the header's limit and the DEL's DELSEQ" \
  "$(grep -c '^##INFO=<ID=DELSEQ.*left out above 60 bp' "$OUT/cmax.region.vcf"):$(awk -F'\t' '!/^#/ && /SVTYPE=DEL;/ && /DELSEQ=/' "$OUT/cmax.region.vcf" | wc -l | tr -d ' ')" "1:0"

# The interval in the reference's name: pggb, odgi and impg write a 0-based start with a half-open end
# (end - start = the path's length), panvar's own W-line names are 1-based inclusive (end - start + 1).
# The path's length tells them apart; one matching neither is read as half-open, with a warning. Read
# as 1-based, the six real loci came out 1 bp low: their REF bases match GRCh38 one base further on.
for nm in 'ref#0#chr1:1000-1540' 'ref#0#chr1:1001-1540' 'ref#0#chr1:1000-1600'; do
  sed "s/$REF_NAME/$nm/" "$OUT/cb.sorted.gfa" > "$OUT/cname.gfa"
  "$BIN" call -i "$OUT/cname.gfa" -c "$OUT/cb.bubbles.csv" -r "$nm" -o "$OUT/cname" -q >/dev/null 2>"$OUT/cname.err"
  printf '%s:%s ' "$(awk -F'\t' '!/^#/ && /SVTYPE=DEL;/{print $2; exit}' "$OUT/cname.region.vcf")" \
    "$(grep -c 'does not match the path' "$OUT/cname.err" | tr -d ' ')"
done > "$OUT/cname.txt"
check "reference name interval: half-open, 1-based, neither (DEL POS:warnings)" "$(cat "$OUT/cname.txt")" "1100:0 1100:0 1100:1 "

# END-POS must equal the deleted length, i.e. the record's own interval closes on the event.
dp=$(field DEL 2); de=$(info DEL END)
[ -n "$dp" ] && [ -n "$de" ] && [ "$((de - dp))" = "70" ] \
  && ok "DEL END-POS equals the deleted length (70)" \
  || bad "DEL END-POS: expected 70, got $((${de:-0} - ${dp:-0}))"
ip=$(field INV 2); ie=$(info INV END)
[ -n "$ip" ] && [ -n "$ie" ] && [ "$((ie - ip))" = "80" ] \
  && ok "INV END-POS equals the inverted length (80)" \
  || bad "INV END-POS: expected 80, got $((${ie:-0} - ${ip:-0}))"

# DELSEQ/INVSEQ are the sequences the record says it changes; they must be what the reference
# actually holds at POS+1..END. Independent of merging, and it pins the offset directly.
[ "$(info DEL DELSEQ)" = "$N2" ] && ok "DELSEQ is exactly n2, the deleted node" \
                                 || bad "DELSEQ does not equal n2"
[ "$(info INS INSSEQ)" = "$NI" ] && ok "INSSEQ is exactly nI, the inserted node" \
                                 || bad "INSSEQ does not equal nI"

# ---------------------------------------------------------------------------------------------
# Input contracts. Each of these previously exited 0 and wrote a header-only VCF, which is
# indistinguishable from a locus with no variation -- the failure mode that makes a silent contract
# breach worse than a crash.
# ---------------------------------------------------------------------------------------------

# refuses <label> <expected-substring> <args...> -- must exit non-zero, say why, and write nothing
refuses() {
  local label="$1" want="$2"; shift 2
  local pfx="$OUT/refuse.$((++refuse_n))"
  local out rc
  out=$("$@" -o "$pfx" 2>&1); rc=$?
  if [ "$rc" -eq 0 ]; then bad "$label: expected a non-zero exit, got 0"; return; fi
  case "$out" in *"$want"*) ;; *) bad "$label: message did not mention '$want' (got: $(printf '%s' "$out" | head -1))"; return;; esac
  if [ -e "$pfx.region.vcf" ]; then bad "$label: refused but still wrote a region VCF"; return; fi
  ok "$label"
}
refuse_n=0

refuses "--bubble-id naming no such bubble is refused" "not in" \
  "$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" --bubble-id 999 -q

# a bubbles CSV whose nodes are not in this graph
sed 's/^\([0-9]*\),\([0-9]*\),/\1,77\2,/' "$OUT/cb.bubbles.csv" > "$OUT/alien.csv"
refuses "a bubbles CSV from another graph is refused" "does not belong to this GFA" \
  "$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/alien.csv" -r "$REF_NAME" -q

{ head -1 "$OUT/cb.bubbles.csv"; sed -n '2p' "$OUT/cb.bubbles.csv"; sed -n '2p' "$OUT/cb.bubbles.csv"; } > "$OUT/dup.csv"
refuses "a duplicate bubble id is refused" "duplicate bubble id" \
  "$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/dup.csv" -r "$REF_NAME" -q

# a graph whose path steps name a node with no S line: spell() would skip it and silently shorten
# every sequence derived from that walk, so the run must not start
awk -F'\t' '!($1=="S" && $2=="3")' "$OUT/cb.sorted.gfa" > "$OUT/missing_node.gfa"
refuses "a graph missing a referenced node is refused" "not in the graph" \
  "$BIN" call -i "$OUT/missing_node.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" -q

# an output that would overwrite an input
cp "$OUT/cb.sorted.gfa" "$OUT/alias.region.vcf"
out=$("$BIN" call -i "$OUT/alias.region.vcf" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" -o "$OUT/alias" -q 2>&1)
case "$out" in
  *"is also an input"*) [ -s "$OUT/alias.region.vcf" ] \
      && ok "an output that aliases an input is refused, input intact" \
      || bad "refused the alias but the input was destroyed anyway" ;;
  *) bad "an output aliasing an input was not refused (got: $(printf '%s' "$out" | head -1))" ;;
esac

# Stale per-bubble files. A --bubble-id run must not leave the previous full run's per-bubble VCFs
# lying beside it, where nothing distinguishes them from output this run produced.
"$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" -o "$OUT/stale" -q >/dev/null 2>&1
before=$(ls "$OUT"/stale.bubble_*.vcf 2>/dev/null | wc -l | tr -d ' ')
"$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" -o "$OUT/stale" --bubble-id 1 -q >/dev/null 2>&1
after=$(ls "$OUT"/stale.bubble_*.vcf 2>/dev/null | wc -l | tr -d ' ')
[ "$before" = "3" ] && [ "$after" = "1" ] \
  && ok "a narrowed rerun clears the previous run's per-bubble VCFs ($before -> $after)" \
  || bad "stale per-bubble VCFs: expected 3 -> 1, got $before -> $after"

# Number=A: AC and SVLEN are emitted as one value per ALT on a multiallelic record, so Number=1 was
# a lie any spec-compliant reader would trip over.
for k in AC AF SVLEN; do
  got=$(grep -o "##INFO=<ID=$k,Number=[^,]*" "$VCF" | head -1 | sed 's/.*Number=//')
  [ "$got" = "A" ] && ok "INFO/$k is declared Number=A" || bad "INFO/$k: expected Number=A, got ${got:-<missing>}"
done

# Option contracts. A fraction that is not a fraction used to be accepted: nan then compared false
# against everything and silently disabled the filter it was meant to set.
for bad in nan 2 -0.5; do
  refuses "--min-alt-af $bad is refused" "must be in [0,1]" \
    "$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" --min-alt-af "$bad" -q
done
refuses "--minimap-preset must be one of the advertised three" "must be asm5" \
  "$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" --minimap-preset bogus -q

# --min-alt-af is tested on the record's AF, not a minor-allele frequency; the old name is gone.
refuses "--min-maf is no longer an option" "Unknown option" \
  "$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" --min-maf 0.1 -q

# A rerun with --no-per-bubble-vcf must not leave the previous run's per-bubble files looking current.
"$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" -o "$OUT/npb" -q >/dev/null 2>&1
npb_before=$(ls "$OUT"/npb.bubble_*.vcf 2>/dev/null | wc -l | tr -d ' ')
"$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" -o "$OUT/npb" --no-per-bubble-vcf -q >/dev/null 2>&1
npb_after=$(ls "$OUT"/npb.bubble_*.vcf 2>/dev/null | wc -l | tr -d ' ')
[ "$npb_before" = "3" ] && [ "$npb_after" = "0" ] \
  && ok "--no-per-bubble-vcf clears the previous run's per-bubble VCFs ($npb_before -> $npb_after)" \
  || bad "--no-per-bubble-vcf stale files: expected 3 -> 0, got $npb_before -> $npb_after"

# Output must not depend on thread count. The per-bubble loop is parallel and writes into per-bubble
# slots, so a shared-state regression would show up as a reordering rather than a crash.
"$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" -o "$OUT/t1" --threads 1 -q >/dev/null 2>&1
"$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" -o "$OUT/t8" --threads 8 -q >/dev/null 2>&1
if cmp -s "$OUT/t1.region.vcf" "$OUT/t8.region.vcf"; then
  ok "1 thread and 8 threads give byte-identical region VCFs"
else bad "output differs between 1 and 8 threads"; fi


# ---------------------------------------------------------------------------------------------
# Per-gene copy number at a self-loop UNIT node carrying BOTH genes of a near-identical pair in every
# copy (GREP = GA + GLINK + GB). The per-site split hands each module copy to one gene, so here it
# halved both: ref (2 + 2) read 1/1 and ghapB (3 + 3) read 2/1, while the rows' own dosage said 2 and
# 3. The pooled dosages add up to twice the total, so the pair now stays pooled, and ghapD, which
# deletes the module, reads a MEASURED zero: none of either gene's private k-mers.
# ---------------------------------------------------------------------------------------------
# Deterministic pseudo-random sequence. The repeating-motif helper above will not do here: a gene
# built from one 40-mer has only ~40 distinct k-mers, so the private sets are degenerate and the
# per-site split cannot resolve for anyone. Same seed always gives the same sequence.
gblk() { awk -v n="$1" -v seed="${2:-1}" 'BEGIN{
           b="ACGT"; x=seed*1103515245+12345;
           for(i=0;i<n;i++){ x=(x*1103515245+12345)%2147483648; if(x<0)x=-x;
                             printf "%s", substr(b, (int(x/65536)%4)+1, 1) } }'; }

GFL=$(gblk 400 5); GLINK=$(gblk 200 13); GFR=$(gblk 400 19)
GA=$(gblk 1200 2)
# gene B: gene A with five substituted bases, so the pair is near-identical (routed to per-site)
# but still has divergent sites to split on
GB="$GA"
for p in 100 300 500 700 900; do
  c="${GB:$p:1}"; case "$c" in A) n=C;; C) n=G;; G) n=T;; *) n=A;; esac
  GB="${GB:0:$p}$n${GB:$((p+1))}"
done
GREP="$GA$GLINK$GB"
gstart=$((${#GFL} + 1))
ga1=$((gstart + ${#GA} - 1))
gb0=$((gstart + ${#GA} + ${#GLINK}))
gb1=$((gb0 + ${#GB} - 1))
gtotal=$((${#GFL} + 2 * ${#GREP} + ${#GFR}))

{ printf 'H\tVN:Z:1.0\n'
  printf "S\tf1\t%s\nS\trep\t%s\nS\tf2\t%s\n" "$GFL" "$GREP" "$GFR"
  printf 'L\tf1\t+\trep\t+\t0M\nL\trep\t+\trep\t+\t0M\nL\trep\t+\tf2\t+\t0M\nL\tf1\t+\tf2\t+\t0M\n'
  printf 'P\tgref#0#chr9:1-%d\tf1+,rep+,rep+,f2+\t*\n' "$gtotal"
  printf 'P\tghapA#1#chr9\tf1+,rep+,rep+,f2+\t*\n'
  printf 'P\tghapB#1#chr9\tf1+,rep+,rep+,rep+,f2+\t*\n'
  printf 'P\tghapD#1#chr9\tf1+,f2+\t*\n'; } > "$OUT/gene.gfa"
{ printf 'chr9\ttest\tgene\t%d\t%d\t.\t+\t.\tgene_id "GENEA"; gene_name "GENEA";\n' "$gstart" "$ga1"
  printf 'chr9\ttest\tCDS\t%d\t%d\t.\t+\t0\tgene_id "GENEA"; gene_name "GENEA";\n' "$gstart" "$ga1"
  printf 'chr9\ttest\tgene\t%d\t%d\t.\t+\t.\tgene_id "GENEB"; gene_name "GENEB";\n' "$gb0" "$gb1"
  printf 'chr9\ttest\tCDS\t%d\t%d\t.\t+\t0\tgene_id "GENEB"; gene_name "GENEB";\n' "$gb0" "$gb1"; } > "$OUT/gene.gtf"

"$BIN" bubble -i "$OUT/gene.gfa" -r 'gref#0#chr9' -o "$OUT/gb" --merge-nearby-bp 0 -q >/dev/null 2>&1
"$BIN" call -i "$OUT/gb.sorted.gfa" -c "$OUT/gb.bubbles.csv" -r 'gref#0#chr9' -o "$OUT/gc" \
      --cn --gtf "$OUT/gene.gtf" -q >/dev/null 2>&1
GTSV="$OUT/gc.dup_gene_cn.tsv"
if [ ! -s "$GTSV" ]; then
  bad "the two-gene UNIT fixture produced no dup_gene_cn table"
else
  # a haplotype that DOES carry the module resolves and is reliable
  carrier_rel=$(awk -F'\t' '$3 ~ /ghapB/ && $4=="GENEA" {print $6; exit}' "$GTSV")
  [ "$carrier_rel" = "1" ] && ok "a haplotype carrying the module reports a reliable per-gene split" \
    || bad "carrier per-gene split: expected reliable=1, got ${carrier_rel:-<missing>}"
  gcn() { awk -F'\t' -v h="$2" -v g="$3" '$3 ~ h && $4==g {print $5"/"$6"/"$7; exit}' "$1"; }
  got="$(gcn "$GTSV" gref GENEA) $(gcn "$GTSV" gref GENEB) $(gcn "$GTSV" ghapB GENEA) $(gcn "$GTSV" ghapB GENEB)"
  [ "$got" = "2/1/2.00 2/1/2.00 3/1/3.00 3/1/3.00" ] \
    && ok "both genes in every copy: each gene counts every copy (2/2, 3/3), not half of them" \
    || bad "two genes per copy: expected 2/1/2.00 2/1/2.00 3/1/3.00 3/1/3.00 (cn/reliable/dosage), got $got"
  [ "$(gcn "$GTSV" ghapD GENEA)" = "0/1/0.00" ] \
    && ok "a haplotype without the module reads a measured zero (no private k-mer hit)" \
    || bad "module-loss haplotype: expected 0/1/0.00, got $(gcn "$GTSV" ghapD GENEA)"

  # Exons: the same genes annotated as three CDS pieces each. Joined end to end, every exon boundary
  # made 30 k-mers that exist nowhere in the genome and counted as markers no haplotype could hit, so
  # ghapB read about 2.8 copies of each gene instead of 3.
  : > "$OUT/gene3.gtf"
  for gn in GENEA:$gstart:$ga1 GENEB:$gb0:$gb1; do
    IFS=: read -r name lo hi <<< "$gn"
    printf 'chr9\ttest\tgene\t%d\t%d\t.\t+\t.\tgene_id "%s"; gene_name "%s";\n' "$lo" "$hi" "$name" "$name" >> "$OUT/gene3.gtf"
    for piece in "$lo $((lo + 349))" "$((lo + 450)) $((lo + 799))" "$((lo + 900)) $hi"; do
      set -- $piece
      printf 'chr9\ttest\tCDS\t%d\t%d\t.\t+\t0\tgene_id "%s"; gene_name "%s";\n' "$1" "$2" "$name" "$name" >> "$OUT/gene3.gtf"
    done
  done
  "$BIN" call -i "$OUT/gb.sorted.gfa" -c "$OUT/gb.bubbles.csv" -r 'gref#0#chr9' -o "$OUT/gc3" \
        --cn --gtf "$OUT/gene3.gtf" -q >/dev/null 2>&1
  [ "$(gcn "$OUT/gc3.dup_gene_cn.tsv" ghapB GENEA)" = "3/1/3.00" ] \
    && ok "a multi-exon gene's dosage is exact (no exon-junction k-mers)" \
    || bad "multi-exon gene: expected 3/1/3.00 for ghapB, got $(gcn "$OUT/gc3.dup_gene_cn.tsv" ghapB GENEA)"

  # A gene lying only in the left flank (the bubble's source) is not in the module: it was given the
  # module's copy number, reliable, though no copy of it varies.
  { cat "$OUT/gene.gtf"
    printf 'chr9\ttest\tgene\t50\t300\t.\t+\t.\tgene_id "FLANKX"; gene_name "FLANKX";\n'; } > "$OUT/geneF.gtf"
  "$BIN" call -i "$OUT/gb.sorted.gfa" -c "$OUT/gb.bubbles.csv" -r 'gref#0#chr9' -o "$OUT/gcF" \
        --cn --gtf "$OUT/geneF.gtf" -q >/dev/null 2>&1
  [ -s "$OUT/gcF.dup_gene_cn.tsv" ] && ! grep -q FLANKX "$OUT/gcF.dup_gene_cn.tsv" \
    && ! grep -v '^#' "$OUT/gcF.region.vcf" | grep -q 'FLANKX' \
    && ok "a gene only on the site's flank gets no per-gene row and no GENES tag" \
    || bad "a flank-only gene was reported: $(grep FLANKX "$OUT/gcF.dup_gene_cn.tsv" | head -2 | tr '\n' '|')"

  # A sibling's intron carrying part of GENEA, as GSTM2/4/5 carry 340 of GSTM1's 754 CDS k-mers. GENEB
  # here is a single-copy gene whose node X every haplotype has; 300 bp of GENEA sit in its intron.
  # Private against GENEB's CDS only, those k-mers read as a quarter of a copy of GENEA on ghapD, which
  # has none; screened against GENEB's whole span, they are not markers.
  GX=$(gblk 1000 31)
  GXI="${GX:0:400}${GA:0:300}${GX:700}"
  xstart=$((${#GFL} + 1)); xcds1=$((xstart + 399)); xend=$((xstart + ${#GXI} - 1))
  astart=$((xend + 1)); aend=$((astart + ${#GA} - 1))
  bgtotal=$((${#GFL} + ${#GXI} + 2 * (${#GA} + ${#GLINK}) + ${#GFR}))
  { printf 'H\tVN:Z:1.0\n'
    printf "S\tf1\t%s\nS\tx\t%s\nS\trep\t%s\nS\tf2\t%s\n" "$GFL" "$GXI" "$GA$GLINK" "$GFR"
    printf 'L\tf1\t+\tx\t+\t0M\nL\tx\t+\trep\t+\t0M\nL\trep\t+\trep\t+\t0M\nL\trep\t+\tf2\t+\t0M\nL\tx\t+\tf2\t+\t0M\n'
    printf 'P\tgref#0#chr9:1-%d\tf1+,x+,rep+,rep+,f2+\t*\n' "$bgtotal"
    printf 'P\tghapB#1#chr9\tf1+,x+,rep+,rep+,rep+,f2+\t*\n'
    printf 'P\tghapD#1#chr9\tf1+,x+,f2+\t*\n'; } > "$OUT/geneBG.gfa"
  { printf 'chr9\ttest\tgene\t%d\t%d\t.\t+\t.\tgene_id "GENEB"; gene_name "GENEB";\n' "$xstart" "$xend"
    printf 'chr9\ttest\tCDS\t%d\t%d\t.\t+\t0\tgene_id "GENEB"; gene_name "GENEB";\n' "$xstart" "$xcds1"
    printf 'chr9\ttest\tgene\t%d\t%d\t.\t+\t.\tgene_id "GENEA"; gene_name "GENEA";\n' "$astart" "$aend"
    printf 'chr9\ttest\tCDS\t%d\t%d\t.\t+\t0\tgene_id "GENEA"; gene_name "GENEA";\n' "$astart" "$aend"; } > "$OUT/geneBG.gtf"
  # The site is written by hand, flank to flank, so that GENEB's node is inside it as GSTM2/4/5 are
  # inside the GSTM module; `bubble` would make x the site's boundary, and boundary genes are not
  # targets.
  printf 'bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,f1,+,f2,+,3,0,%d,"rep;x"\n' \
    "$bgtotal" > "$OUT/gbBG.bubbles.csv"
  "$BIN" call -i "$OUT/geneBG.gfa" -c "$OUT/gbBG.bubbles.csv" -r 'gref#0#chr9' -o "$OUT/gcBG" \
        --cn --gtf "$OUT/geneBG.gtf" -q >/dev/null 2>&1
  bg="$(gcn "$OUT/gcBG.dup_gene_cn.tsv" ghapD GENEA) $(gcn "$OUT/gcBG.dup_gene_cn.tsv" ghapB GENEA)"
  [ "$bg" = "0/1/0.00 3/1/3.00" ] \
    && ok "k-mers a sibling's intron carries are not counted as copies of the gene (0.00 and 3.00)" \
    || bad "sibling-intron k-mers: expected '0/1/0.00 3/1/3.00' for ghapD, ghapB GENEA, got '$bg'"
fi

# ---------------------------------------------------------------------------------------------
# Per-site split of a near-identical pair with ONE gene per module copy: shared chunks, 1 bp allele
# nodes at nine divergent sites, a link node closing each copy. h3 carries A and two copies of B in
# which sites 0-3 were converted to A's allele, so its pooled dosage reads 1.89 A / 1.11 B; the
# per-site median across sites out-votes the converted ones and gives the truth, 1 A / 2 B. With B
# inverted in the reference (the pair on opposite strands), the alignment of B to A is reverse-strand;
# reading its CIGAR forward lost every site and fell back to the pooled 2/1 without notice.
# `noflank` carries the copies but not the flanks, so it does not cross the site and has no module
# total: it was split as a total of 0, a confident 0/0.
# `own` puts every copy on nodes of its own (same sequence), apart from its first chunk and the link: no
# path revisits the genes' nodes, as where a graph kept near-identical copies apart. The genes still lie
# in the copies; measured on the revisited nodes they did not, and the pair was left to its pooled 2/1.
# The graph's only cycle is the module, which `bubble` does not report as a site on its own (see the
# open points), so the bubbles CSV is written by hand: the site between the two flanks.
# ---------------------------------------------------------------------------------------------
cat > "$OUT/pair_gen.py" <<'PYEOF'
import random, sys, json
out = sys.argv[1]; haps = json.loads(sys.argv[2]); D = int(sys.argv[3]); rnd = random.Random(int(sys.argv[4]))
CONV = set(int(x) for x in sys.argv[5].split(',') if x)
OWN = len(sys.argv) > 6 and sys.argv[6] == 'own'   # each copy on its own nodes, apart from c0 and lk
def rs(n): return ''.join(rnd.choice('ACGT') for _ in range(n))
def rc(s): return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
L = 1500
pos = [100 + i * (L - 200) // (D - 1) for i in range(D)]
G = rs(L); segs = {}; prev = 0
for i, p in enumerate(pos):
    segs['c%d' % i] = G[prev:p]; a = G[p]
    segs['a%d' % i] = a; segs['b%d' % i] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[a]; prev = p + 1
segs['c%d' % D] = G[prev:]; segs['f1'] = rs(500); segs['f2'] = rs(500); segs['lk'] = rs(300)
def copy_steps(c):
    st = []
    for i in range(D):
        al = 'a' if c[0] == 'A' or (c[0] == 'C' and i in CONV) else 'b'
        st += ['c%d' % i, '%s%d' % (al, i)]
    st.append('c%d' % D)
    return [s + '+' for s in st] if c[1] == '+' else [s + '-' for s in reversed(st)]
def own_nodes(st, tag):
    out = []
    for s in st:
        n = s[:-1]
        if n != 'c0':
            segs[n + tag] = segs[n]; n += tag
        out.append(n + s[-1])
    return out
def path_steps(name, copies, flanks):
    st = ['f1+'] if flanks else []
    for i, c in enumerate(copies):
        cs = copy_steps(c)
        st += (own_nodes(cs, '_%s%d' % (name, i)) if OWN else cs) + ['lk+']
    return st + (['f2+'] if flanks else [])
seq = lambda st: ''.join(segs[s[:-1]] if s[-1] == '+' else rc(segs[s[:-1]]) for s in st)
paths = {n: path_steps(n, c, not n.startswith('noflank')) for n, c in haps.items()}
edges = set((x, y) for st in paths.values() for x, y in zip(st, st[1:]))
reflen = len(seq(paths['ref']))
used = set(s[:-1] for st in paths.values() for s in st)
segs = {k: v for k, v in segs.items() if k in used}
with open(out + '.gfa', 'w') as f:
    f.write('H\tVN:Z:1.0\n')
    for k, v in segs.items(): f.write('S\t%s\t%s\n' % (k, v))
    for x, y in sorted(edges): f.write('L\t%s\t%s\t%s\t%s\t0M\n' % (x[:-1], x[-1], y[:-1], y[-1]))
    for n, st in paths.items():
        f.write('P\t%s\t%s\t*\n' % ('ref#0#chr9:1-%d' % reflen if n == 'ref' else n + '#1#chr9', ','.join(st)))
inside = sorted(k for k in segs if k not in ('f1', 'f2'))
with open(out + '.csv', 'w') as f:
    f.write('bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n')
    f.write('1,f1,+,f2,+,%d,0,%d,"%s"\n' % (len(paths), reflen, ';'.join(inside)))
off = len(segs['f1']) + 1
with open(out + '.gtf', 'w') as f:
    for gi, c in enumerate(haps['ref'][:2]):
        name = 'GENEA' if gi == 0 else 'GENEB'; s = off; e = off + L - 1
        for kind in ('gene', 'CDS'):
            f.write('chr9\tt\t%s\t%d\t%d\t.\t%s\t.\tgene_id "%s"; gene_name "%s";\n' % (kind, s, e, c[1], name, name))
        off = e + 1 + len(segs['lk'])
PYEOF
for t in fwd inv own; do
  o=$([ "$t" = inv ] && echo - || echo +)
  python3 "$OUT/pair_gen.py" "$OUT/pair_$t" \
    "{\"ref\":[\"A+\",\"B$o\"],\"h3\":[\"A+\",\"C$o\",\"C$o\"],\"h4\":[\"A+\",\"A+\",\"B$o\"],\"noflankE\":[\"A+\",\"B$o\"]}" \
    9 7 0,1,2,3 $t
  "$BIN" call -i "$OUT/pair_$t.gfa" -c "$OUT/pair_$t.csv" -r 'ref#0#chr9' -o "$OUT/pc_$t" \
        --gtf "$OUT/pair_$t.gtf" --cn -q >/dev/null 2>&1
  PT="$OUT/pc_$t.dup_gene_cn.tsv"
  pcn() { awk -F'\t' -v h="$1" -v g="$2" '$3 ~ h && $4==g {print $5"/"$6; exit}' "$PT"; }
  got="$(pcn h3 GENEA) $(pcn h3 GENEB) $(pcn h4 GENEA) $(pcn h4 GENEB)"
  [ "$got" = "1/1 2/1 2/1 1/1" ] \
    && ok "per-site split, pair $t: converted copies are out-voted (h3 1 A / 2 B, h4 2 A / 1 B)" \
    || bad "per-site split, pair $t: expected '1/1 2/1 2/1 1/1' (cn/reliable), got '$got'"
  [ "$(pcn noflankE GENEA)" = "./0" ] \
    && ok "per-site split, pair $t: a haplotype with no module total is unresolved (./0), not 0" \
    || bad "per-site split, pair $t: no-total haplotype expected ./0, got $(pcn noflankE GENEA)"
done

# ---------------------------------------------------------------------------------------------
# Hybrid copies, read exon by exon. Two paralogs of six 150 bp coding exons; their introns are nodes every
# copy shares, their exons differ at every 24th base. hhyb carries a third copy, GENEA with GENEB's last
# exon: one crossover, a hybrid, marked for both genes while the counts stay as the dosage gives them.
# hconv carries GENEA with GENEB's third exon: conversion inside the copy, no hybrid.
# ---------------------------------------------------------------------------------------------
cat > "$OUT/hyb_gen.py" <<'PYEOF'
import json, random, sys
out = sys.argv[1]; rnd = random.Random(11)
step = int(sys.argv[2]) if len(sys.argv) > 2 else 24   # paralog exons differ at every step-th base
rs = lambda n: ''.join(rnd.choice('ACGT') for _ in range(n))
mut = lambda s: ''.join({'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[c] if i % step == 5 else c for i, c in enumerate(s))
segs = {'f1': rs(800), 'f2': rs(800), 'lk': rs(400)}
for i in range(6):
    segs['ea%d' % i] = rs(150); segs['eb%d' % i] = mut(segs['ea%d' % i])
    if i < 5: segs['in%d' % i] = rs(250)
def copy(word):   # word: per exon 'a', 'b' or 'c'; upper case drops the intron after that exon; '~' a 12 kb spacer;
    if word == '~': return ['sp+']   # a leading 'r' writes the copy reverse-complemented
    if word.startswith('r'):
        return [s[:-1] + ('-' if s[-1] == '+' else '+') for s in reversed(copy(word[1:]))]
    st = []
    for i, c in enumerate(word):
        st.append('e%s%d+' % (c.lower(), i))
        if i < 5 and c.islower(): st.append('in%d+' % i)
    return st + ['lk+']
haps = {'ref': ['aaaaaa', 'bbbbbb'], 'h1': ['aaaaaa', 'bbbbbb'], 'h2': ['aaaaaa', 'aaaaaa', 'bbbbbb'],
        'h3': ['aaaaaa', 'bbbbbb', 'bbbbbb'], 'hhyb': ['aaaaaa', 'aaaaab', 'bbbbbb'],
        'hconv': ['aaaaaa', 'aabaaa', 'bbbbbb']}
haps.update(json.loads(sys.argv[3]) if len(sys.argv) > 3 else {})
if any(w == '~' for ws in haps.values() for w in ws): segs['sp'] = rs(12000)
# 'c': a version of GENEA differing from it at two bases per exon, so it keeps about half its private k-mers
if any('c' in w.lower() for ws in haps.values() for w in ws):
    for i in range(6):
        segs['ec%d' % i] = ''.join({'A': 'G', 'C': 'T', 'G': 'A', 'T': 'C'}[c] if j % 60 == 33 else c
                                   for j, c in enumerate(segs['ea%d' % i]))
paths = {n: ['f1+'] + [s for w in ws for s in copy(w)] + ['f2+'] for n, ws in haps.items()}
edges = set((x, y) for st in paths.values() for x, y in zip(st, st[1:]))
rc = lambda q: q[::-1].translate(str.maketrans('ACGT', 'TGCA'))
seq = lambda st: ''.join(segs[s[:-1]] if s[-1] == '+' else rc(segs[s[:-1]]) for s in st)
reflen = len(seq(paths['ref']))
with open(out + '.gfa', 'w') as f:
    f.write('H\tVN:Z:1.0\n')
    for k, v in segs.items(): f.write('S\t%s\t%s\n' % (k, v))
    for x, y in sorted(edges): f.write('L\t%s\t%s\t%s\t%s\t0M\n' % (x[:-1], x[-1], y[:-1], y[-1]))
    for n, st in paths.items():
        f.write('P\t%s\t%s\t*\n' % ('ref#0#chr9:1-%d' % reflen if n == 'ref' else n + '#1#chr9', ','.join(st)))
with open(out + '.csv', 'w') as f:
    f.write('bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n')
    f.write('1,f1,+,f2,+,%d,0,%d,"%s"\n' % (len(paths), reflen, ';'.join(sorted(k for k in segs if k not in ('f1', 'f2')))))
with open(out + '.gtf', 'w') as f:
    at = len(segs['f1']) + 1
    for name, c in (('GENEA', 'a'), ('GENEB', 'b')):
        lo = at
        for i in range(6):
            hi = at + 149
            f.write('chr9\tt\tCDS\t%d\t%d\t.\t+\t.\tgene_id "%s"; gene_name "%s";\n' % (at, hi, name, name))
            at = hi + 1 + (250 if i < 5 else 0)
        f.write('chr9\tt\tgene\t%d\t%d\t.\t+\t.\tgene_id "%s"; gene_name "%s";\n' % (lo, at - 1, name, name))
        at += len(segs['lk'])
PYEOF
python3 "$OUT/hyb_gen.py" "$OUT/hyb"
"$BIN" call -i "$OUT/hyb.gfa" -c "$OUT/hyb.csv" -r 'ref#0#chr9' -o "$OUT/hyb" --gtf "$OUT/hyb.gtf" --cn -q >/dev/null 2>&1
hyb_of() { awk -F'\t' -v h="$1" -v g="$2" 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "hybrid_copies") c = i }
  $3 ~ ("^" h "#") && $4 == g { print $5 "/" $c; exit }' "$OUT/hyb.dup_gene_cn.tsv"; }
check "hybrid copies: GENEA and GENEB cn/hybrids for h2, hhyb, hconv" \
  "$(hyb_of h2 GENEA) $(hyb_of h2 GENEB) $(hyb_of hhyb GENEA) $(hyb_of hhyb GENEB) $(hyb_of hconv GENEA) $(hyb_of hconv GENEB)" \
  "2/0 1/0 2/1 1/1 2/0 1/0"

# GENES on a module of paralogs (exons differing at every 60th base): an event inside a gene of the copies
# names the paralogs its carriers hold, not the reference copy's gene. Each carrier drops one intron: hda in
# its GENEA copy, hdb in its GENEB copy, hbb (two GENEB copies) in its first, which the walk places on GENEA.
python3 "$OUT/hyb_gen.py" "$OUT/par" 60 \
  '{"hda": ["aaAaaa", "bbbbbb"], "hdb": ["aaaaaa", "bbBbbb"], "hbb": ["bbBbbb", "bbbbbb"]}'
"$BIN" call -i "$OUT/par.gfa" -c "$OUT/par.csv" -r 'ref#0#chr9' -o "$OUT/par" --gtf "$OUT/par.gtf" --cn -q >/dev/null 2>&1
par_genes() { awk -F'\t' -v h="$1" '/^#CHROM/ { for (i = 10; i <= NF; i++) if ($i ~ ("^" h "#")) c = i; next }
  !/^#/ && $8 !~ /SVTYPE=CNV/ && $c ~ /^1/ { match($8, /GENES=[^;]*/); print substr($8, RSTART + 6, RLENGTH - 6) }' \
  "$OUT/par.region.vcf"; }
check "GENES on paralog copies: the carriers' genes (hda, hdb, hbb)" \
  "$(par_genes hda) $(par_genes hdb) $(par_genes hbb)" "GENEA,GENEB GENEA,GENEB GENEB"

# Copies of a GENEA version keeping about half its private k-mers ('c'): two of them pool to about one copy's
# worth of hits. Read as runs, copies 12 kb apart count each (hcc, and hc1 beside a reference-like copy);
# adjacent copies are one run and keep the pooled count (hcn). A GENEB copy holding two GENEA exons is a run of
# GENEA's hits closer to GENEB, not a copy (hcv). A hybrid copy keeps the pooled count, though its run and the
# two others read as copies (hch).
python3 "$OUT/hyb_gen.py" "$OUT/low" 24 \
  '{"hcc": ["cccccc", "~", "cccccc", "bbbbbb"], "hc1": ["aaaaaa", "~", "cccccc"], "hcn": ["cccccc", "cccccc", "bbbbbb"],
    "hcv": ["aaaaaa", "~", "bbaabb"], "hch": ["cccccc", "~", "cccccc", "~", "aaaaab", "bbbbbb"]}'
"$BIN" call -i "$OUT/low.gfa" -c "$OUT/low.csv" -r 'ref#0#chr9' -o "$OUT/low" --gtf "$OUT/low.gtf" --cn -q >/dev/null 2>&1
low_cn() { awk -F'\t' -v h="$1" 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "copies") c = i }
  $3 ~ ("^" h "#") && $4 == "GENEA" { print $5 "/" $c; exit }' "$OUT/low.dup_gene_cn.tsv"; }
check "low-share copies of GENEA, cn/copies: hcc, hc1, hcn" "$(low_cn hcc) $(low_cn hc1) $(low_cn hcn)" "2/2 2/2 1/."
check "runs of GENEA's hits that are not copies of it, cn/copies: hcv (closer to GENEB), hch (a hybrid)" \
  "$(low_cn hcv) $(low_cn hch)" "1/1 2/3"
# The same low-share copies on the other strand: one of the two reversed (hrc), or both (hrr).
python3 "$OUT/hyb_gen.py" "$OUT/lowr" 24 \
  '{"hrc": ["cccccc", "~", "rcccccc", "bbbbbb"], "hrr": ["rcccccc", "~", "rcccccc", "bbbbbb"]}'
"$BIN" call -i "$OUT/lowr.gfa" -c "$OUT/lowr.csv" -r 'ref#0#chr9' -o "$OUT/lowr" --gtf "$OUT/lowr.gtf" --cn -q >/dev/null 2>&1
lowr_cn() { awk -F'\t' -v h="$1" 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "copies") c = i }
  $3 ~ ("^" h "#") && $4 == "GENEA" { print $5 "/" $c; exit }' "$OUT/lowr.dup_gene_cn.tsv"; }
check "low-share copies of GENEA on the reverse strand, cn/copies: hrc, hrr" "$(lowr_cn hrc) $(lowr_cn hrr)" "2/2 2/2"

# ---------------------------------------------------------------------------------------------
# A reference that REVISITS an event node. Every node->position map in `call` records a node's FIRST
# occurrence, so an event sitting at a later copy was anchored at the earlier one -- in sequence the
# haplotype still carries.
#
#   f1 100 bp    1..100
#   R   50 bp  101..150     first occurrence of R, nowhere near the event
#   f2 100 bp  151..250
#   A   70 bp  251..320  \  the deletion covers 251..370: A and the SECOND R
#   R   50 bp  321..370  /
#   f3 100 bp  371..470
#
# Anchored on the preceding base: POS 250, END 370, SVLEN -120. Taking R's first occurrence instead
# put the record at POS 100 / END 220, 150 bp upstream.
# ---------------------------------------------------------------------------------------------
# Seeds chosen so neither event can left-align (the base before each differs from its last base).
RF1=$(blk 100 0); RR=$(blk 50 8); RF2=$(blk 100 13); RA=$(blk 70 21); RF3=$(blk 100 29)
RI=$(gblk 90 4)
RREF_NAME='rref#0#chr1:1-470'
{ printf 'H\tVN:Z:1.0\n'
  printf "S\tf1\t%s\nS\tR\t%s\nS\tf2\t%s\nS\tA\t%s\nS\tf3\t%s\nS\tI\t%s\n" \
         "$RF1" "$RR" "$RF2" "$RA" "$RF3" "$RI"
  for e in "f1 R" "R f2" "f2 A" "A R" "R f3" "f2 f3" "R I" "I f3"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"
  done
  printf 'P\t%s\tf1+,R+,f2+,A+,R+,f3+\t*\n' "$RREF_NAME"
  printf 'P\trhapref#1#chr1\tf1+,R+,f2+,A+,R+,f3+\t*\n'
  printf 'P\trhapdel#1#chr1\tf1+,R+,f2+,f3+\t*\n'
  printf 'P\trhapdel2#1#chr1\tf1+,R+,f2+,f3+\t*\n'
  # an insertion after the SECOND visit to R, i.e. anchored at 370, not at R's first visit (150)
  printf 'P\trhapins#1#chr1\tf1+,R+,f2+,A+,R+,I+,f3+\t*\n'
  printf 'P\trhapins2#1#chr1\tf1+,R+,f2+,A+,R+,I+,f3+\t*\n'; } > "$OUT/rep.gfa"

"$BIN" bubble -i "$OUT/rep.gfa" -r "$RREF_NAME" -o "$OUT/rb" --merge-nearby-bp 0 -q >/dev/null 2>&1
"$BIN" call -i "$OUT/rb.sorted.gfa" -c "$OUT/rb.bubbles.csv" -r "$RREF_NAME" -o "$OUT/rc" -q >/dev/null 2>&1
RVCF="$OUT/rc.region.vcf"
if [ ! -s "$RVCF" ]; then
  bad "the repeated-anchor fixture produced no region VCF"
else
  rpos=$(awk -F'\t' '$0!~/^#/ && $8~/SVTYPE=DEL/{print $2; exit}' "$RVCF")
  rend=$(awk -F'\t' '$0!~/^#/ && $8~/SVTYPE=DEL/{n=split($8,a,";");
           for(i=1;i<=n;i++){split(a[i],kv,"="); if(kv[1]=="END"){print kv[2]; exit}}}' "$RVCF")
  [ "$rpos" = "250" ] \
    && ok "a deletion at the SECOND visit to a node anchors there, not at the first (POS 250)" \
    || bad "repeated-anchor POS: expected 250, got ${rpos:-<missing>} (100 = anchored at the first visit)"
  [ "$rend" = "370" ] \
    && ok "its END follows the same occurrence (370)" \
    || bad "repeated-anchor END: expected 370, got ${rend:-<missing>}"
  # The insertion sits after the SECOND visit to R, so it anchors on that visit's last base (370).
  # Anchoring on R's first visit would put it at 150, inside sequence the haplotype also carries.
  ipos=$(awk -F'\t' '$0!~/^#/ && $8~/SVTYPE=INS/{print $2; exit}' "$RVCF")
  [ "$ipos" = "370" ] \
    && ok "an insertion after the SECOND visit anchors there too (POS 370)" \
    || bad "repeated-anchor INS POS: expected 370, got ${ipos:-<missing>} (150 = R's first visit)"
  # The DEL and the INS sit at different occurrences and must not have been merged into one record.
  nrec=$(awk -F'\t' '$0!~/^#/' "$RVCF" | wc -l | tr -d ' ')
  [ "$nrec" = "2" ] \
    && ok "events at different occurrences stay separate records ($nrec)" \
    || bad "expected 2 records at distinct occurrences, got $nrec"
fi

# The same repeated-node graph on the REVERSE strand, discovered with --no-flip so the bubble really
# is reverse-oriented (source_orient/sink_orient '-') rather than normalised to forward. The walk is
# then traversed in decreasing coordinate, so an index derived as `left + j` would be wrong; the
# coordinates must come out identical to the forward case.
{ printf 'H\tVN:Z:1.0\n'
  printf "S\tf1\t%s\nS\tR\t%s\nS\tf2\t%s\nS\tA\t%s\nS\tf3\t%s\n" "$RF1" "$RR" "$RF2" "$RA" "$RF3"
  for e in "f3 R" "R A" "A f2" "f2 R" "R f1" "f3 f2"; do
    set -- $e; printf "L\t%s\t-\t%s\t-\t0M\n" "$1" "$2"
  done
  printf 'P\trvref#0#chr1:1-470\tf3-,R-,A-,f2-,R-,f1-\t*\n'
  printf 'P\trvhapref#1#chr1\tf3-,R-,A-,f2-,R-,f1-\t*\n'
  printf 'P\trvhapdel#1#chr1\tf3-,f2-,R-,f1-\t*\n'
  printf 'P\trvhapdel2#1#chr1\tf3-,f2-,R-,f1-\t*\n'; } > "$OUT/rev.gfa"
"$BIN" bubble -i "$OUT/rev.gfa" -r 'rvref#0#chr1:1-470' -o "$OUT/rvb" --no-flip --merge-nearby-bp 0 -q >/dev/null 2>&1
rvorient=$(tail -1 "$OUT/rvb.bubbles.csv" 2>/dev/null | cut -d, -f3)
[ "$rvorient" = "-" ] \
  && ok "--no-flip yields a genuinely reverse-oriented bubble to test against" \
  || bad "expected a reverse-oriented bubble (source_orient '-'), got '${rvorient:-<none>}'"
"$BIN" call -i "$OUT/rvb.sorted.gfa" -c "$OUT/rvb.bubbles.csv" -r 'rvref#0#chr1:1-470' \
      -o "$OUT/rvc" -q >/dev/null 2>&1
rvpos=$(awk -F'\t' '$0!~/^#/ && $8~/SVTYPE=DEL/{print $2; exit}' "$OUT/rvc.region.vcf")
rvend=$(awk -F'\t' '$0!~/^#/ && $8~/SVTYPE=DEL/{n=split($8,a,";");
          for(i=1;i<=n;i++){split(a[i],kv,"="); if(kv[1]=="END"){print kv[2]; exit}}}' "$OUT/rvc.region.vcf")
{ [ "$rvpos" = "100" ] && [ "$rvend" = "220" ]; } \
  && ok "a reverse-oriented bubble gives the same coordinates as the forward one ($rvpos/$rvend)" \
  || bad "reverse-oriented coordinates: expected 100/220, got ${rvpos:-?}/${rvend:-?}"

# ---------------------------------------------------------------------------------------------
# A module the reference holds once. The reference walks u x (260 bp) once; pkdup walks u x u x u, NON-
# consecutively so no self-loop claims it: two whole copies and a third u. Counted at the module's anchor
# it has 2 copies, and the extra u (60 bp) is an insertion. The busiest-node count (the PEAK route) read 3.
# ---------------------------------------------------------------------------------------------
PF1=$(gblk 300 41); PU=$(gblk 60 43); PX=$(gblk 200 47); PF2=$(gblk 300 53)
{ printf 'H\tVN:Z:1.0\n'
  printf "S\tpf1\t%s\nS\tpu\t%s\nS\tpx\t%s\nS\tpf2\t%s\n" "$PF1" "$PU" "$PX" "$PF2"
  for e in "pf1 pu" "pu px" "px pu" "px pf2" "pu pf2"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"
  done
  printf 'P\tpkref#0#chr8:1-860\tpf1+,pu+,px+,pf2+\t*\n'
  printf 'P\tpkflat#1#chr8\tpf1+,pu+,px+,pf2+\t*\n'
  printf 'P\tpkdup#1#chr8\tpf1+,pu+,px+,pu+,px+,pu+,pf2+\t*\n'
  printf 'P\tpkdup2#1#chr8\tpf1+,pu+,px+,pu+,px+,pu+,pf2+\t*\n'; } > "$OUT/peak.gfa"
"$BIN" bubble -i "$OUT/peak.gfa" -r 'pkref#0#chr8:1-860' -o "$OUT/pkb" --merge-nearby-bp 0 -q >/dev/null 2>&1
"$BIN" call -i "$OUT/pkb.sorted.gfa" -c "$OUT/pkb.bubbles.csv" -r 'pkref#0#chr8:1-860' \
      -o "$OUT/pkc" --cn -q >/dev/null 2>&1
PKV="$OUT/pkc.region.vcf"
pk_method=$(awk -F'\t' '$0!~/^#/ && $8~/SVTYPE=CNV/{n=split($8,a,";");
              for(i=1;i<=n;i++){split(a[i],kv,"="); if(kv[1]=="CN_METHOD"){print kv[2]; exit}}}' "$PKV")
[ "$pk_method" = "MODULE" ] \
  && ok "a non-consecutive module duplication is a MODULE record" \
  || bad "expected CN_METHOD=MODULE, got '${pk_method:-<no DUP record>}'"
# the duplicating haplotype must read 2 copies, the flat one 1 -- a positive AND a negative
pk_cn() { awk -F'\t' -v s="$1" '/^#CHROM/{for(i=10;i<=NF;i++) if($i ~ s) c=i}
          $0!~/^#/ && $8~/SVTYPE=CNV/{split($9,k,":"); for(j in k) if(k[j]=="CN") ci=j
            split($c,v,":"); print v[ci]; exit}' "$PKV"; }
[ "$(pk_cn pkdup)" = "2" ] \
  && ok "MODULE reports 2 copies for the haplotype with two whole copies and a third u" \
  || bad "MODULE CN for the duplicating haplotype: expected 2, got '$(pk_cn pkdup)'"
[ "$(pk_cn pkflat)" = "1" ] \
  && ok "MODULE reports 1 copy for the haplotype that does not duplicate" \
  || bad "MODULE CN for the non-duplicating haplotype: expected 1, got '$(pk_cn pkflat)'"
# The third u lies between the duplicating haplotype's copies: it goes with them, in the module record's forms.
check "the third u is in the duplicating haplotype's copies, no insertion of its own (INS records, pkdup's CNFORM)" \
  "$(grep -v '^#' "$PKV" | grep -c 'SVTYPE=INS;') $(awk -F'\t' '/^#CHROM/{for(i=10;i<=NF;i++) if($i ~ /^pkdup#/) c=i}
     !/^#/ && /SVTYPE=CNV;/{n=split($9,k,":"); split($c,v,":"); for(j=1;j<=n;j++) if(k[j]=="CNFORM") print v[j]}' "$PKV")" "0 1,2"

# ---------------------------------------------------------------------------------------------
# Within-haplotype coalescing must also judge nearness at the OCCURRENCE, not at a node's first
# visit. Reference f1 R f2 R f3 visits R twice; one haplotype inserts after the FIRST R and again
# after the SECOND. By first occurrence both anchors resolve to the same base, so the gap reads 0
# and they coalesce into a single insertion carrying both sequences at the earlier position.
#
# Also the DEL-cluster contract, which had no focused test: two deletions with reference retained
# between them coalesce into one record that is IMPRECISE, whose END spans the whole affected
# interval while SVLEN stays the net deleted bases.
# ---------------------------------------------------------------------------------------------
CF1=$(blk 100 2); CR=$(blk 50 6); CF2=$(blk 60 10); CF3=$(blk 100 14)
CI1=$(gblk 70 61); CI2=$(gblk 80 67)
CD1=$(blk 60 18); CK=$(blk 40 22); CD2=$(blk 60 26)
{ printf 'H\tVN:Z:1.0\n'
  printf "S\tc1\t%s\nS\tcR\t%s\nS\tc2\t%s\nS\tc3\t%s\nS\ti1\t%s\nS\ti2\t%s\n" \
         "$CF1" "$CR" "$CF2" "$CF3" "$CI1" "$CI2"
  printf "S\td1\t%s\nS\tck\t%s\nS\td2\t%s\n" "$CD1" "$CK" "$CD2"
  for e in "c1 cR" "cR c2" "c2 cR" "cR c3" "cR i1" "i1 c2" "cR i2" "i2 c3" \
           "c3 d1" "d1 ck" "ck d2" "d2 c1x" "c3 ck" "ck c1x" "c3 c1x"; do
    set -- $e; printf "L\t%s\t+\t%s\t+\t0M\n" "$1" "$2"
  done
  printf "S\tc1x\t%s\n" "$(blk 100 30)"
  printf 'P\tcref#0#chr9:1-620\tc1+,cR+,c2+,cR+,c3+,d1+,ck+,d2+,c1x+\t*\n'
  printf 'P\tchapref#1#chr9\tc1+,cR+,c2+,cR+,c3+,d1+,ck+,d2+,c1x+\t*\n'
  printf 'P\tchapins#1#chr9\tc1+,cR+,i1+,c2+,cR+,i2+,c3+,d1+,ck+,d2+,c1x+\t*\n'
  printf 'P\tchapdel#1#chr9\tc1+,cR+,c2+,cR+,c3+,ck+,c1x+\t*\n'
  # a haplotype that removes the whole span, so c3..c1x is ONE snarl containing d1, ck and d2 --
  # otherwise ck is a cut vertex, the two deletions land in different bubbles, and within-haplotype
  # coalescing (which is per bubble) never gets to see them
  printf 'P\tchapwipe#1#chr9\tc1+,cR+,c2+,cR+,c3+,c1x+\t*\n'; } > "$OUT/coal.gfa"
"$BIN" bubble -i "$OUT/coal.gfa" -r 'cref#0#chr9:1-620' -o "$OUT/cob" --merge-nearby-bp 0 -q >/dev/null 2>&1
# Coalescing joins pieces the node comparison found apart, so it is tested there (--no-align-alleles);
# aligned at base level, the two deletions are two (checked below).
"$BIN" call -i "$OUT/cob.sorted.gfa" -c "$OUT/cob.bubbles.csv" -r 'cref#0#chr9:1-620' \
      -o "$OUT/coc" --no-align-alleles -q >/dev/null 2>&1
"$BIN" call -i "$OUT/cob.sorted.gfa" -c "$OUT/cob.bubbles.csv" -r 'cref#0#chr9:1-620' \
      -o "$OUT/coca" -q >/dev/null 2>&1
check "aligned at base level, chapdel's two 60 bp deletions with retained sequence between are two DEL records" \
  "$(grep -v '^#' "$OUT/coca.region.vcf" | grep -c 'SVLEN=-60;')" "2"
CVCF="$OUT/coc.region.vcf"
if [ ! -s "$CVCF" ]; then
  bad "the coalescing fixture produced no region VCF"
else
  n_ins=$(awk -F'\t' '$0!~/^#/ && $8~/SVTYPE=INS/' "$CVCF" | wc -l | tr -d ' ')
  [ "$n_ins" = "2" ] \
    && ok "two insertions at distinct occurrences stay two records (not coalesced into one)" \
    || bad "expected 2 INS records at distinct anchors, got $n_ins"
  # the DEL cluster: coalesced across retained sequence, so IMPRECISE with END wider than |SVLEN|
  # the cluster is the DEL carried by the haplotype that kept `ck`: net 120 bp over a 160 bp span
  dsv=$(awk -F'\t' '$0!~/^#/ && $8~/SVLEN=-120/{n=split($8,a,";");
          for(i=1;i<=n;i++){split(a[i],kv,"="); if(kv[1]=="SVLEN"){print -kv[2]; exit}}}' "$CVCF")
  dspan=$(awk -F'\t' '$0!~/^#/ && $8~/SVLEN=-120/{n=split($8,a,";");
            for(i=1;i<=n;i++){split(a[i],kv,"="); if(kv[1]=="END"){print kv[2]-$2; exit}}}' "$CVCF")
  dimp=$(awk -F'\t' '$0!~/^#/ && $8~/SVLEN=-120/ && $8~/IMPRECISE/{print "yes"; exit}' "$CVCF")
  if [ -z "$dsv" ]; then
    bad "the fixture produced no DEL record to test the cluster contract on"
  else
    [ "$dsv" = "120" ] \
      && ok "a coalesced DEL cluster reports the NET deleted bases as SVLEN (120)" \
      || bad "DEL cluster SVLEN: expected 120 (two 60 bp deletions), got $dsv"
    [ "$dspan" -gt "$dsv" ] \
      && ok "its END spans the whole affected interval ($dspan > $dsv), including retained sequence" \
      || bad "DEL cluster END-POS=$dspan should exceed |SVLEN|=$dsv"
    [ "$dimp" = "yes" ] \
      && ok "and it is flagged IMPRECISE, so END-POS == |SVLEN| still holds everywhere else" \
      || bad "a DEL cluster whose span exceeds |SVLEN| was not flagged IMPRECISE"
    # the 120 bases it removes lie somewhere in the span: each breakpoint up to the slack inward
    dci=$(awk -F'\t' '$0!~/^#/ && $8~/SVLEN=-120/{n=split($8,a,";"); for(i=1;i<=n;i++){split(a[i],kv,"=");
            if(kv[1]=="CIPOS") p=kv[2]; if(kv[1]=="CIEND") e=kv[2]} print p" "e; exit}' "$CVCF")
    check "its CIPOS/CIEND give the slack, END-POS less |SVLEN|" "$dci" "0,$((dspan - dsv)) -$((dspan - dsv)),0"
  fi
fi


# ---------------------------------------------------------------- repeated module boundary
# CN, CNBP, CN_MODULE_REF_BP, POS and END share one occurrence-aware resolver, which deliberately
# selects the WIDEST span (first source to last sink) when a boundary recurs. No reference locus
# produces a CN record with a repeated boundary, so nothing exercised that choice.
#
# The reference here visits S and K TWICE, with a full module between each pair. Nothing is
# self-adjacent, so the UNIT route cannot fire and the collapsed-module route is used.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%s\nS\tP\t%s\nS\tQ\t%s\nS\tK\t%s\nS\tG\t%s\n' \
    "$(gblk 150 21)" "$(gblk 50 22)" "$(gblk 150 23)" "$(gblk 150 24)" "$(gblk 50 25)" "$(gblk 150 26)"
  for e in "F S" "S P" "P Q" "Q K" "K S" "K G" "Q P" "P K"; do
    set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tgrch38#0#chr1:1-1100\tF+,S+,P+,Q+,K+,S+,P+,Q+,K+,G+\t*\n'
  printf 'P\th1#1#chr1\tF+,S+,P+,Q+,P+,Q+,K+,S+,P+,Q+,K+,G+\t*\n'
  printf 'P\th2#1#chr1\tF+,S+,P+,K+,S+,P+,Q+,K+,G+\t*\n'; } > "$OUT/modrep.gfa"
printf 'bubble_id,source,source_orient,sink,sink_orient,inside_node_count,total_node_count,path_support,distinct_alleles,ref_allele_support,alt_allele_support_max,alt_allele_support_min,min_inside_bp,max_inside_bp,inside_nodes\n1,S,+,K,+,2,4,3,3,1,1,1,150,450,"P;Q"\n' > "$OUT/modrep.bubbles.csv"
"$BIN" call -i "$OUT/modrep.gfa" -c "$OUT/modrep.bubbles.csv" \
   --reference-path "grch38#0#chr1:1-1100" -o "$OUT/modrep" --cn --min-sv-bp 20 -q >/dev/null 2>&1
MRV="$OUT/modrep.region.vcf"
inf() { awk -F'\t' -v k="$1" '$0!~/^#/ && $8~/SVTYPE=CNV/{n=split($8,a,";");
          for(i=1;i<=n;i++){split(a[i],kv,"="); if(kv[1]==k){print kv[2]; exit}}}' "$MRV"; }
fmtv() { awk -F'\t' -v s="$1" -v k="$2" '$0~/^#CHROM/{for(i=10;i<=NF;i++) if($i==s) col=i}
          $0!~/^#/ && $8~/SVTYPE=CNV/{split($9,f,":"); for(j in f) if(f[j]==k) ki=j
            split($col,v,":"); print v[ki]; exit}' "$MRV"; }
[ "$(inf CN_METHOD)" = "MODULE" ] \
  && ok "a repeated module boundary still routes to MODULE" \
  || bad "CN_METHOD=$(inf CN_METHOD), expected MODULE"
[ "$(inf CN_SPAN_AMBIGUOUS)" = "3" ] \
  && ok "CN_SPAN_AMBIGUOUS counts every traverser with a repeated boundary (3)" \
  || bad "CN_SPAN_AMBIGUOUS=$(inf CN_SPAN_AMBIGUOUS), expected 3"
[ "$(inf CN_MODULE_REF_BP)" = "600" ] \
  && ok "CN_MODULE_REF_BP sums the module nodes over both visits (600)" \
  || bad "CN_MODULE_REF_BP=$(inf CN_MODULE_REF_BP), expected 600"
[ "$(fmtv h1#1#chr1 CNBP)" = "300" ] && [ "$(fmtv h2#1#chr1 CNBP)" = "-150" ] \
  && ok "CNBP is measured over the same node set (+300 / -150)" \
  || bad "CNBP h1=$(fmtv h1#1#chr1 CNBP) h2=$(fmtv h2#1#chr1 CNBP), expected 300 / -150"
# h1 gains a whole copy; h2 keeps the anchor P but loses Q in one copy: a partial copy CN does not show.
[ "$(fmtv h1#1#chr1 CNRESID)" = "0" ] && [ "$(fmtv h2#1#chr1 CNRESID)" = "-150" ] \
  && ok "MODULE CNRESID: 0 for a whole copy, -150 for a partial one" \
  || bad "CNRESID h1=$(fmtv h1#1#chr1 CNRESID) h2=$(fmtv h2#1#chr1 CNRESID), expected 0 / -150"
# POS/END are the widest SPAN, which with a repeated boundary encloses the intervening K and S. So
# END-POS EXCEEDS CN_MODULE_REF_BP by exactly that intervening sequence. The identity
# END-POS == CN_MODULE_REF_BP holds only when every boundary is visited once, and CN_SPAN_AMBIGUOUS
# is precisely the flag that says it does not.
mr_pos=$(awk -F'\t' '$0!~/^#/ && $8~/SVTYPE=CNV/{print $2; exit}' "$MRV")
mr_span=$(( $(inf END) - mr_pos ))
[ "$mr_span" = "700" ] \
  && ok "END-POS is the widest span (700), not the node sum" \
  || bad "END-POS=$mr_span, expected 700"
[ "$(( mr_span - $(inf CN_MODULE_REF_BP) ))" = "100" ] \
  && ok "and the excess is exactly the boundary visits between the outermost ones (K+S=100)" \
  || bad "span excess $(( mr_span - $(inf CN_MODULE_REF_BP) )), expected 100"

# A haplotype crossing the module sink to source opens its span on the sink. hrev passes S once upstream
# (F S X), then walks G- K- (Q- P-)x2 S- F-: its first S comes before its K, and reading the direction from
# that took the span S X G K, which holds no module node, so hrev read as CN 0 (-600 bp) instead of 2.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%s\nS\tP\t%s\nS\tQ\t%s\nS\tK\t%s\nS\tG\t%s\nS\tX\t%s\n' \
    "$(gblk 150 81)" "$(gblk 50 82)" "$(gblk 150 83)" "$(gblk 150 84)" "$(gblk 50 85)" "$(gblk 150 86)" "$(gblk 60 87)"
  for e in "F + S +" "S + P +" "P + Q +" "Q + P +" "Q + K +" "K + G +" "S + X +" "X + G -"; do
    set -- $e; printf 'L\t%s\t%s\t%s\t%s\t0M\n' "$1" "$2" "$3" "$4"; done
  printf 'P\tref#0#chr1:1-1000\tF+,S+,P+,Q+,P+,Q+,K+,G+\t*\n'
  printf 'P\threv#1#chr1\tF+,S+,X+,G-,K-,Q-,P-,Q-,P-,S-,F-\t*\n'
  printf 'P\thgain#1#chr1\tF+,S+,P+,Q+,P+,Q+,P+,Q+,K+,G+\t*\n'; } > "$OUT/modrev.gfa"
printf '%s\n1,S,+,K,+,2,4,3,2,1,1,1,300,900,"P;Q"\n' "$(head -1 "$OUT/modrep.bubbles.csv")" > "$OUT/modrev.bubbles.csv"
"$BIN" call -i "$OUT/modrev.gfa" -c "$OUT/modrev.bubbles.csv" -r 'ref#0#chr1:1-1000' -o "$OUT/modrev" \
  --cn --no-per-bubble-vcf -q >/dev/null 2>&1
MRV="$OUT/modrev.region.vcf"
check "a reverse crossing opens the module span on the sink (hrev CN:CNBP, AC)" \
  "$(fmtv hrev#1#chr1 CN):$(fmtv hrev#1#chr1 CNBP),$(inf AC)" "2:0,1"

# ---------------------------------------------------------------- short UNIT unit, large copy change
# A self-loop is removed from the DEL/INS/INV walks, so the UNIT route is the only one that can see
# its copies. It used to skip any unit shorter than --min-sv-bp, which silently erased a 30 bp unit
# gaining 8 copies (+240 bp); with --cn the inferred routes then answered with a partial second
# record. The record's SVLEN must also not depend on the order the haplotypes are listed in.
rep_gfa() {   # rep_gfa <unit_bp> <ref_copies> <out.gfa> <hap>:<copies> ...
  local u=$1 rc=$2 out=$3; shift 3
  local R1; R1=$(blk 100 0); local RU; RU=$(blk "$u" 5); local R3; R3=$(blk 100 13)
  reps() { local k=$1 s=""; for ((i=0;i<k;i++)); do s="$s,r+"; done; printf '%s' "$s"; }
  { printf 'H\tVN:Z:1.0\nS\tn1\t%s\nS\tr\t%s\nS\tn3\t%s\n' "$R1" "$RU" "$R3"
    printf 'L\tn1\t+\tr\t+\t0M\nL\tr\t+\tr\t+\t0M\nL\tr\t+\tn3\t+\t0M\nL\tn1\t+\tn3\t+\t0M\n'
    printf 'P\tref#0#chr1:1-%d\tn1+%s,n3+\t*\n' $((200 + rc * u)) "$(reps "$rc")"
    for h in "$@"; do printf 'P\t%s#1#chr1\tn1+%s,n3+\t*\n' "${h%%:*}" "$(reps "${h#*:}")"; done
  } > "$out"
}
rep_call() {   # rep_call <gfa> <prefix> [call options...]
  local g=$1 p=$2; shift 2
  "$BIN" bubble -i "$g" -r ref --min-variant-bp 0 -o "$p.b" --merge-nearby-bp 0 --quiet >/dev/null 2>&1
  "$BIN" call -i "$p.b.sorted.gfa" --bubble-prefix "$p.b" -r ref -o "$p" "$@" --quiet >/dev/null 2>&1
}
nrec() { grep -vc '^#' "$1"; }
rinfo() { grep -v '^#' "$1" | head -1 | tr ';\t' '\n\n' | awk -F= -v k="$2" '$1==k{print $2}'; }
rcn() { awk -F'\t' -v s="$2" '/^#CHROM/{for(i=10;i<=NF;i++) if($i==s) c=i} !/^#/{split($9,f,":");
          for(j in f) if(f[j]=="CN") k=j; split($c,v,":"); print v[k]; exit}' "$1"; }
rep_gfa 30 4 "$OUT/rep30.gfa" hexp:12 hcon:1 hzero:0
for mode in "" "--cn"; do
  rep_call "$OUT/rep30.gfa" "$OUT/rep30$mode" $mode
  V="$OUT/rep30$mode.region.vcf"
  [ "$(nrec "$V")" = "1" ] && [ "$(rinfo "$V" CN_METHOD)" = "UNIT" ] \
    && ok "30 bp unit with an SV-sized copy change: one UNIT record ${mode:-(default)}" \
    || bad "30 bp unit ${mode:-(default)}: $(nrec "$V") record(s), method $(rinfo "$V" CN_METHOD), expected 1 UNIT"
  [ "$(rcn "$V" 'hexp#1#chr1'):$(rcn "$V" 'hcon#1#chr1'):$(rcn "$V" 'hzero#1#chr1')" = "12:1:0" ] \
    && ok "and every haplotype's copy number is measured ${mode:-(default)}" \
    || bad "copy numbers $(rcn "$V" 'hexp#1#chr1'):$(rcn "$V" 'hcon#1#chr1'):$(rcn "$V" 'hzero#1#chr1'), expected 12:1:0"
done
# hexp gains 8 copies, hcon and hzero lose 3 and 4: carriers change in both directions, so no one
# size describes the record. SVLEN is left out and SVLEN_RANGE spans the carriers (-4..+8 x 30 bp).
[ "$(rinfo "$OUT/rep30.region.vcf" SVLEN):$(rinfo "$OUT/rep30.region.vcf" SVLEN_RANGE)" = ":-120,240" ] \
  && ok "mixed gains and losses: no SVLEN, SVLEN_RANGE -120,240" \
  || bad "mixed UNIT record: SVLEN '$(rinfo "$OUT/rep30.region.vcf" SVLEN)', SVLEN_RANGE '$(rinfo "$OUT/rep30.region.vcf" SVLEN_RANGE)', expected none and -120,240"
# Carriers that all gain keep SVLEN, the largest change; the range shows the smaller one too.
rep_gfa 30 4 "$OUT/repg.gfa" hexp:12 hmid:6 hsame:4
rep_call "$OUT/repg.gfa" "$OUT/repg"
[ "$(rinfo "$OUT/repg.region.vcf" SVLEN):$(rinfo "$OUT/repg.region.vcf" SVLEN_RANGE)" = "240:60,240" ] \
  && ok "gains only: SVLEN is the largest change (240), SVLEN_RANGE 60,240" \
  || bad "gain-only UNIT record: SVLEN '$(rinfo "$OUT/repg.region.vcf" SVLEN)', SVLEN_RANGE '$(rinfo "$OUT/repg.region.vcf" SVLEN_RANGE)', expected 240 and 60,240"
rep_gfa 60 4 "$OUT/repA.gfa" hlo:0 hhi:8
rep_gfa 60 4 "$OUT/repB.gfa" hhi:8 hlo:0
rep_call "$OUT/repA.gfa" "$OUT/repA"; rep_call "$OUT/repB.gfa" "$OUT/repB"
a=$(rinfo "$OUT/repA.region.vcf" SVLEN):$(rinfo "$OUT/repA.region.vcf" SVLEN_RANGE)
b=$(rinfo "$OUT/repB.region.vcf" SVLEN):$(rinfo "$OUT/repB.region.vcf" SVLEN_RANGE)
[ "$a" = ":-240,240" ] && [ "$b" = ":-240,240" ] \
  && ok "size fields are independent of haplotype order (a -4/+4 record: no SVLEN, range -240,240)" \
  || bad "size fields depend on haplotype order or rule: '$a' vs '$b', expected ':-240,240' for both"

# ---------------------------------------------------------------- module losses are carriers
# A module the reference holds once: h0 loses it (CN 0, a carrier, as on UNIT), h2 and h3 gain copies.
# Sequences are random so that distinct nodes do not look alike to the sequence-identity merge.
# a seeded random sequence, the same on every platform (awk's rand differs between awk implementations)
RSEQ_PY="$(cd "$(dirname "$0")/.." && pwd)/tools/rseq.py"
rseq() { "${PYTHON:-python3}" -I "$RSEQ_PY" "$1" "$2"; }
cn_call() {   # cn_call <gfa> <prefix>
  "$BIN" bubble -i "$1" -r ref --min-variant-bp 0 -o "$2.b" --merge-nearby-bp 0 --quiet >/dev/null 2>&1
  "$BIN" call -i "$2.b.sorted.gfa" --bubble-prefix "$2.b" -r ref -o "$2" --cn --quiet >/dev/null 2>&1
}
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tf0\t%s\nS\tf1\t%s\nS\tma\t%s\nS\tmb\t%s\n' "$(rseq 1000 21)" "$(rseq 10000 22)" "$(rseq 1500 23)" "$(rseq 500 24)"
  for e in "f0 ma" "ma mb" "mb f1" "mb ma" "f0 f1"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  M="ma+,mb+"
  printf 'P\tref#0#chr1:1-13000\tf0+,%s,f1+\t*\n' "$M"
  printf 'P\th0#1#chr1\tf0+,f1+\t*\n'
  printf 'P\th1#1#chr1\tf0+,%s,f1+\t*\n' "$M"
  printf 'P\th2#1#chr1\tf0+,%s,%s,f1+\t*\n' "$M" "$M"
  printf 'P\th3#1#chr1\tf0+,%s,%s,%s,f1+\t*\n' "$M" "$M" "$M"
} > "$OUT/peak.gfa"
cn_call "$OUT/peak.gfa" "$OUT/peak"; V="$OUT/peak.region.vcf"
[ "$(grep -v '^#' "$V" | grep -c 'CN_METHOD=MODULE;')" = "1" ] \
  && ok "--cn: a single-copy reference module is one MODULE DUP" || bad "--cn: no single MODULE DUP in $V"
h0gtcn=$(awk -F'\t' '/^#CHROM/{for(i=10;i<=NF;i++) if($i=="h0#1#chr1") c=i}
  !/^#/ && /CN_METHOD=MODULE;/{split($9,k,":"); split($c,v,":"); for(j in k) f[k[j]]=v[j]; print f["GT"]":"f["CN"]}' "$V")
[ "$h0gtcn" = "1:0" ] && ok "MODULE: the total loss is a carrier (GT 1) with CN 0" \
                      || bad "MODULE: the total loss reads GT:CN $h0gtcn, expected 1:0"
# The reference runs a module of two anchor nodes (a3 20 bp, a2 40 bp) twice. htest deletes a2 of the second
# copy and the 60 bp node x after it: one copy fewer at the anchor (CN 1), but its module changes by 40 bp, under
# --min-sv-bp, so it is no carrier of the module record. Its 100 bp deletion is then reported as a deletion, and
# its CN on the module record is the reference's, so a dosage read off CN does not count the change twice.
{ printf 'N L 2000 81\nN A3 20 82\nN A2 40 83\nN S 40 84\nN X 60 85\nN R 2000 86\n'
  printf 'P ref#0#chr1:0-4220 L+,A3+,A2+,S+,A3+,A2+,X+,R+\nP hfull#1#chr1 L+,A3+,A2+,X+,R+\n'
  printf 'P hgain#1#chr1 L+,A3+,A2+,S+,A3+,A2+,S+,A3+,A2+,X+,R+\nP htest#1#chr1 L+,A3+,A2+,S+,A3+,R+\n'
  for i in 1 2 3 4 5 6 7 8 9 10; do printf 'P hr%s#1#chr1 L+,A3+,A2+,S+,A3+,A2+,X+,R+\n' "$i"; done; } > "$OUT/modlow.spec"
python3 -I "$(cd "$(dirname "$0")/.." && pwd)/tools/spec_gfa.py" "$OUT/modlow.spec" "$OUT/modlow"
cn_call "$OUT/modlow/g.gfa" "$OUT/modlow"; V="$OUT/modlow.region.vcf"
modlow=$(awk -F'\t' '/^#CHROM/{for(i=10;i<=NF;i++) if($i=="htest#1#chr1") c=i}
  !/^#/ {split($c,g,":"); split($8,f,";"); printf "%s:%s:%s ", f[2], g[1], g[2]}' "$V")
[ "$modlow" = "SVTYPE=CNV:0:2 SVTYPE=DEL:1:. " ] \
  && ok "MODULE: a copy change under --min-sv-bp is no carrier, reads REF_CN and keeps its own 100 bp deletion" \
  || bad "MODULE non-carrier: records 'SVTYPE:GT:CN' for htest are '$modlow', expected 'SVTYPE=CNV:0:2 SVTYPE=DEL:1:. '"
# MODBP, the module's own change, is what GT is decided on: htest's -40 bp is under --min-sv-bp, hfull's -100 is not.
check "MODULE: MODBP decides GT (htest, hfull GT:CN:MODBP)" \
  "$(awk -F'\t' '/^#CHROM/{for(i=10;i<=NF;i++) c[$i]=i} !/^#/ && /SVTYPE=CNV;/{n=split($9,k,":"); for(j=1;j<=n;j++) x[k[j]]=j
       for (h in c) if (h == "htest#1#chr1" || h == "hfull#1#chr1") { split($c[h],v,":"); r[h]=v[x["GT"]]":"v[x["CN"]]":"v[x["MODBP"]] }
       print r["htest#1#chr1"], r["hfull#1#chr1"]}' "$V")" "0:2:-40 1:1:-100"
# gt_cn <vcf> <sample> <svtype>: that sample's GT:CN on the first record of the type
gt_cn() { awk -F'\t' -v s="$2" -v t="$3" '/^#CHROM/{for(i=10;i<=NF;i++) if($i==s) c=i}
  !/^#/ && $8 ~ ("SVTYPE="t";") {split($c,g,":"); print g[1]":"g[2]; exit}' "$1"; }
# htest's second a2 is replaced by q, a copy of a2 with 2 bases changed: one copy fewer at the anchor (CN 1) and
# a module change under --min-sv-bp, but no record of its own holds that change, so CN stays as measured.
sed -e 's/^P htest#1#chr1 .*/P htest#1#chr1 L+,A3+,A2+,S+,A3+,Q+,X+,R+/' -e 's/^N A2 40 83$/N A2 40 83\nN Q @A2:2:7/' \
  "$OUT/modlow.spec" > "$OUT/modpara.spec"
python3 -I "$(cd "$(dirname "$0")/.." && pwd)/tools/spec_gfa.py" "$OUT/modpara.spec" "$OUT/modpara"
cn_call "$OUT/modpara/g.gfa" "$OUT/modpara"
check "MODULE non-carrier whose change no record holds keeps its measured CN (GT:CN, records)" \
  "$(gt_cn "$OUT/modpara.region.vcf" 'htest#1#chr1' CNV) $(grep -vc '^#' "$OUT/modpara.region.vcf")" "0:1 1"
# hbig deletes the second copy and x after it, 120 bp: as much unique sequence as copy, so still a carrier of one
# copy fewer, and x (60 bp) a deletion of its own.
sed 's/^P htest#1#chr1 .*/P hbig#1#chr1 L+,A3+,A2+,S+,R+/' "$OUT/modlow.spec" > "$OUT/modbig.spec"
python3 -I "$(cd "$(dirname "$0")/.." && pwd)/tools/spec_gfa.py" "$OUT/modbig.spec" "$OUT/modbig"
cn_call "$OUT/modbig/g.gfa" "$OUT/modbig"
check "MODULE carrier whose deletion holds as much unique sequence as copy: a carrier and the unique part (CNV GT:CN, DEL GT:SVLEN)" \
  "$(gt_cn "$OUT/modbig.region.vcf" 'hbig#1#chr1' CNV) $(gt_cn "$OUT/modbig.region.vcf" 'hbig#1#chr1' DEL | cut -d: -f1):$(awk -F'\t' '!/^#/ && /SVTYPE=DEL;/{match($8,/SVLEN=-?[0-9]+/); print substr($8,RSTART+6,RLENGTH-6)}' "$OUT/modbig.region.vcf")" "1:1 1:-60"
# --diploid: sample s1 holds a carrier (one copy lost) and a non-carrier whose own deletion holds its change. The
# non-carrier's CN reads REF_CN, but its CNRESID comes from its measured copies, as when no record holds it.
sed -e 's/^P hfull#1#chr1 /P s1#1#chr1 /' -e 's/^P htest#1#chr1 /P s1#2#chr1 /' "$OUT/modlow.spec" > "$OUT/moddip.spec"
python3 -I "$(cd "$(dirname "$0")/.." && pwd)/tools/spec_gfa.py" "$OUT/moddip.spec" "$OUT/moddip"
"$BIN" bubble -i "$OUT/moddip/g.gfa" -r ref --min-variant-bp 0 -o "$OUT/moddip.b" --merge-nearby-bp 0 --quiet >/dev/null 2>&1
"$BIN" call -i "$OUT/moddip.b.sorted.gfa" --bubble-prefix "$OUT/moddip.b" -r ref -o "$OUT/moddip" --cn --diploid --quiet >/dev/null 2>&1
check "--diploid: one carrier and one non-carrier holding its change (s1 GT:CN:CNBP:CNRESID:MODBP:CNFORM:HCN)" \
  "$(awk -F'\t' '/^#CHROM/{for(i=10;i<=NF;i++) if($i=="s1") c=i} !/^#/ && /SVTYPE=CNV;/{print $c; exit}' "$OUT/moddip.region.vcf")" "1|0:3:-140:20:-140:2|.:1|2"
"$BIN" call -i "$OUT/moddip.b.sorted.gfa" --bubble-prefix "$OUT/moddip.b" -r ref -o "$OUT/moddip_af" --cn --diploid \
  --min-alt-af 0.1 --quiet >/dev/null 2>&1
check "--diploid, the holding deletion filtered out: the non-carrier's measured CN (s1 GT:CN:...:HCN)" \
  "$(awk -F'\t' '/^#CHROM/{for(i=10;i<=NF;i++) if($i=="s1") c=i} !/^#/ && /SVTYPE=CNV;/{print $c; exit}' "$OUT/moddip_af.region.vcf")" "1|0:2:-140:20:-140:2|.:1|1"
# htest's deletion holds its change only while it is written: --min-alt-af 0.1 leaves it out (AF 1/14), and the
# module record then reads its measured CN; so does --diploid, where the sample's own column then reads 1|1.
"$BIN" call -i "$OUT/modlow.b.sorted.gfa" --bubble-prefix "$OUT/modlow.b" -r ref -o "$OUT/modlow_af" --cn --min-alt-af 0.1 \
  --quiet >/dev/null 2>&1
check "MODULE non-carrier whose holding record is filtered out keeps its measured CN (records, GT:CN)" \
  "$(grep -vc '^#' "$OUT/modlow_af.region.vcf") $(gt_cn "$OUT/modlow_af.region.vcf" 'htest#1#chr1' CNV)" "1 0:1"
# htest walks T for the reference's S, an extra S after x, and q, a copy of a2 with 2 bases changed, for a2. Its
# deletion of S and insertion of S are on module nodes but not on the anchor: no record holds its copy change.
{ printf 'N L 2000 81\nN A3 20 82\nN A2 40 83\nN Q @A2:2:7\nN S 80 84\nN T 20 89\nN X 60 85\nN R 2000 86\n'
  printf 'P ref#0#chr1:0-4300 L+,A3+,A2+,S+,A3+,A2+,X+,R+\nP hfull#1#chr1 L+,A3+,A2+,X+,R+\n'
  printf 'P hgain#1#chr1 L+,A3+,A2+,S+,A3+,A2+,S+,A3+,A2+,X+,R+\nP htest#1#chr1 L+,A3+,A2+,T+,A3+,Q+,X+,S+,R+\n'
  for i in 1 2 3 4 5 6 7 8 9 10; do printf 'P hr%s#1#chr1 L+,A3+,A2+,S+,A3+,A2+,X+,R+\n' "$i"; done; } > "$OUT/modoff.spec"
python3 -I "$(cd "$(dirname "$0")/.." && pwd)/tools/spec_gfa.py" "$OUT/modoff.spec" "$OUT/modoff"
cn_call "$OUT/modoff/g.gfa" "$OUT/modoff"
check "MODULE non-carrier whose records are off the anchor keeps its measured CN (GT:CN)" \
  "$(gt_cn "$OUT/modoff.region.vcf" 'htest#1#chr1' CNV)" "0:1"
# The carrier hfull's own deletion of one copy is its copy change: it stays off the deletion record.
modfull=$(awk -F'\t' '/^#CHROM/{for(i=10;i<=NF;i++) if($i=="hfull#1#chr1") c=i}
  !/^#/ && $8 ~ /SVTYPE=DEL;/ {split($c,g,":"); n=split($8,f,";"); for(i=1;i<=n;i++) if(f[i] ~ /^AC=/) printf "%s ", f[i]; printf "hfull:%s", g[1]}' "$V")
[ "$modfull" = "AC=1 hfull:0" ] \
  && ok "MODULE carrier: its copy loss is not also a deletion record" \
  || bad "MODULE carrier: the deletion record reads '$modfull', expected 'AC=1 hfull:0'"
# The reference runs a 20 bp module (ma, mb: 10 bp each) three times -- one copy below --min-sv-bp, so
# no module record: its copy changes are sequence. hgain's three extra copies are a 60 bp insertion,
# hloss's lost ones a 60 bp deletion, and hnet's, replaced by a 20 bp node x of its own, a substitution
# of net -40 (hloss, inside whose deletion it lies, has no allele there).
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tf0\t%s\nS\tf1\t%s\nS\tma\t%s\nS\tmb\t%s\nS\tx\t%s\n' "$(rseq 1000 31)" "$(rseq 2000 32)" "$(rseq 10 33)" \
    "$(rseq 10 34)" "$(rseq 20 35)"
  for e in "f0 ma" "ma mb" "mb ma" "mb f1" "f0 f1" "f0 x" "x f1"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  M="ma+,mb+"
  printf 'P\tref#0#chr1:1-3060\tf0+,%s,%s,%s,f1+\t*\n' "$M" "$M" "$M"
  printf 'P\thgain#1#chr1\tf0+,%s,%s,%s,%s,%s,%s,f1+\t*\n' "$M" "$M" "$M" "$M" "$M" "$M"
  printf 'P\thnet#1#chr1\tf0+,x+,f1+\t*\n'
  printf 'P\thloss#1#chr1\tf0+,f1+\t*\n'
} > "$OUT/peaknet.gfa"
cn_call "$OUT/peaknet.gfa" "$OUT/peaknet"; V="$OUT/peaknet.region.vcf"
check "a module of 20 bp copies has no record; its copy changes are INS, SUB, DEL (SVTYPE SVLEN GT per record)" \
  "$(grep -v '^#' "$V" | awk -F'\t' '{ match($8, /SVTYPE=[A-Z]+;SVLEN=-?[0-9]+/); g = ""
      for (i = 10; i <= NF; i++) g = g substr($i, 1, 1); printf "%s:%s ", substr($8, RSTART, RLENGTH), g }')" \
  "SVTYPE=INS;SVLEN=60:0100 SVTYPE=SUB;SVLEN=-40:001. SVTYPE=DEL;SVLEN=-60:0001 "

# ---------------------------------------------------------------- copies and the other events of a haplotype
# The reference runs a module of two anchor nodes (A3 60 bp, A2 100 bp) twice, S (40 bp) between the copies and
# X1 X2 X3 (80, 80, 240 bp) after them. Each haplotype's copies go to the module record and its other changes to
# records of their own:
#   hgain, hloss   one copy more, one fewer: the module record alone
#   hhalf          the second copy and X1 deleted: one copy fewer, and X1 (half a copy) a deletion of its own
#   heq            the second copy and X1 X2: as much unique sequence as copy, still one copy fewer and a deletion
#   hmore          the second copy and X1 X2 X3: more unique sequence than copy, so the deletion's carrier
#   hdir           one copy more where X1 X2 X3 are deleted: the deletion goes the other way, so a copy gain
#   hinv           one copy more, turned round
#   hyy            Y inserted in each copy: two insertions, one record each
{ printf 'N L 2000 81\nN A3 60 82\nN A2 100 83\nN S 40 84\nN X1 80 85\nN X2 80 86\nN X3 240 87\nN Y 70 89\nN R 2000 88\n'
  printf 'P ref#0#chr1:0-4680 L+,A3+,A2+,S+,A3+,A2+,X1+,X2+,X3+,R+\n'
  printf 'P hgain#1#chr1 L+,A3+,A2+,S+,A3+,A2+,S+,A3+,A2+,X1+,X2+,X3+,R+\n'
  printf 'P hloss#1#chr1 L+,A3+,A2+,X1+,X2+,X3+,R+\n'
  printf 'P hhalf#1#chr1 L+,A3+,A2+,S+,X2+,X3+,R+\nP heq#1#chr1 L+,A3+,A2+,S+,X3+,R+\nP hmore#1#chr1 L+,A3+,A2+,S+,R+\n'
  printf 'P hdir#1#chr1 L+,A3+,A2+,S+,A3+,A2+,S+,A3+,A2+,R+\n'
  printf 'P hinv#1#chr1 L+,A3+,A2+,S+,A3+,A2+,S-,A2-,A3-,X1+,X2+,X3+,R+\n'
  printf 'P hyy#1#chr1 L+,A3+,Y+,A2+,S+,A3+,Y+,A2+,X1+,X2+,X3+,R+\n'
  for i in 1 2 3 4 5 6 7 8; do printf 'P hr%s#1#chr1 L+,A3+,A2+,S+,A3+,A2+,X1+,X2+,X3+,R+\n' "$i"; done; } > "$OUT/sep.spec"
python3 -I "$(cd "$(dirname "$0")/.." && pwd)/tools/spec_gfa.py" "$OUT/sep.spec" "$OUT/sep"
cn_call "$OUT/sep/g.gfa" "$OUT/sep"
# hrec <vcf> <haplotype>: the module record's GT:CN, then SVTYPE:SVLEN of every other record it carries
hrec() { awk -F'\t' -v s="$2" '/^#CHROM/{for(i=10;i<=NF;i++) if($i==s) c=i}
  !/^#/ {split($c,g,":"); if ($8 ~ /SVTYPE=CNV;/) cnv=g[1]":"g[2]
         else if (g[1]=="1") { match($8,/SVTYPE=[A-Z]+/); t=substr($8,RSTART+7,RLENGTH-7); match($8,/SVLEN=-?[0-9]+/)
                               o=o" "t":"substr($8,RSTART+6,RLENGTH-6) } }
  END{print cnv o}' "$1"; }
V="$OUT/sep.region.vcf"
check "copies: a gain and a loss are the module record alone (hgain, hloss)" "$(hrec "$V" 'hgain#1#chr1')|$(hrec "$V" 'hloss#1#chr1')" "1:3|1:1"
check "copies: a loss with half a copy of unique sequence is a carrier and a deletion of that sequence (hhalf)" \
  "$(hrec "$V" 'hhalf#1#chr1')" "1:1 DEL:-80"
check "copies: a loss with as much unique sequence as copy is still a carrier (heq)" "$(hrec "$V" 'heq#1#chr1')" "1:1 DEL:-160"
check "copies: a loss inside more unique sequence than copy is that deletion's carrier (hmore)" \
  "$(hrec "$V" 'hmore#1#chr1')" "0:2 DEL:-560"
check "copies: a gain beside a deletion of the other direction is a carrier and that deletion (hdir)" \
  "$(hrec "$V" 'hdir#1#chr1')" "1:3 DEL:-400"
check "copies: an inverted gained copy is the module record alone (hinv)" "$(hrec "$V" 'hinv#1#chr1')" "1:3"
check "copies: one insertion in each copy is two records, one event each (hyy; its carrier events)" \
  "$(hrec "$V" 'hyy#1#chr1') $(awk -F'\t' '$NF ~ /hyy/ && $3 != "cn"' "$OUT/sep.carrier_events.tsv" | cut -f1 | sort | uniq -c | awk '{printf "%s", $1}')" "0:2 INS:70 INS:70 11"
check "copies: the module record names each carrier's copies; hhalf's is the second copy's deletion (its CN_FORMS entry)" \
  "$(awk -F'\t' '/^#CHROM/{for(i=10;i<=NF;i++) if($i=="hhalf#1#chr1") c=i} !/^#/ && /SVTYPE=CNV;/{n=split($9,k,":"); split($c,v,":")
       for(j=1;j<=n;j++) if(k[j]=="CNFORM") f=v[j]; match($8,/CN_FORMS=[^;]*/); split(substr($8,RSTART+9,RLENGTH-9),fs,","); print fs[f]}' "$V")" "2201-2360:"
# The same graph, the reference walked backwards: the records hold the same changes.
awk '$1=="P"{n=split($3,s,","); o=""; for(i=n;i>=1;i--){x=s[i]; c=substr(x,length(x)); x=substr(x,1,length(x)-1) (c=="+"?"-":"+")
            o=o (o==""?"":",") x}; print $1, $2, o; next} {print}' "$OUT/sep.spec" > "$OUT/seprev.spec"
python3 -I "$(cd "$(dirname "$0")/.." && pwd)/tools/spec_gfa.py" "$OUT/seprev.spec" "$OUT/seprev"
cn_call "$OUT/seprev/g.gfa" "$OUT/seprev"
check "copies, the reference walked backwards: hhalf, heq, hmore, hdir as forward" \
  "$(for h in hhalf heq hmore hdir; do hrec "$OUT/seprev.region.vcf" "$h#1#chr1"; done | paste -sd'|' -)" \
  "1:1 DEL:-80|1:1 DEL:-160|0:2 DEL:-560|1:3 DEL:-400"
# Two pieces of one change inside a copy, each under --min-sv-bp: hsm deletes P1 (30 bp) and, 8 bp on, replaces
# P2 (30 bp) by Q2 (10 bp). Joined they are one 50 bp change, a record; apart each would be dropped.
{ printf 'N L 2000 81\nN P1 30 90\nN K 8 91\nN P2 30 92\nN Q2 10 93\nN A 100 83\nN S 40 84\nN X 200 87\nN R 2000 88\n'
  printf 'P ref#0#chr1:0-4600 L+,P1+,K+,P2+,A+,S+,P1+,K+,P2+,A+,X+,R+\nP hsm#1#chr1 L+,P1+,K+,P2+,A+,S+,K+,Q2+,A+,X+,R+\n'
  printf 'P hgain#1#chr1 L+,P1+,K+,P2+,A+,S+,P1+,K+,P2+,A+,S+,P1+,K+,P2+,A+,X+,R+\n'
  for i in 1 2 3 4 5 6; do printf 'P hr%s#1#chr1 L+,P1+,K+,P2+,A+,S+,P1+,K+,P2+,A+,X+,R+\n' "$i"; done; } > "$OUT/sepsm.spec"
python3 -I "$(cd "$(dirname "$0")/.." && pwd)/tools/spec_gfa.py" "$OUT/sepsm.spec" "$OUT/sepsm"
cn_call "$OUT/sepsm/g.gfa" "$OUT/sepsm"
check "copies: two pieces under --min-sv-bp 8 bp apart inside a copy are one record (hsm)" \
  "$(hrec "$OUT/sepsm.region.vcf" 'hsm#1#chr1')" "0:2 DEL:-50"
# A gained copy repeats the reference's copies, so its form names those stretches instead of spelling them, and
# --max-seq-bp under the copy leaves it spelled: hgain's copy forward, hinv's turned round.
"$BIN" call -i "$OUT/sep.b.sorted.gfa" --bubble-prefix "$OUT/sep.b" -r ref -o "$OUT/sep_short" --cn --max-seq-bp 100 --quiet >/dev/null 2>&1
check "copies: forms of gained copies name the reference stretches they repeat (hgain, hinv CN_FORMS entries)" \
  "$(awk -F'\t' '/^#CHROM/{for(i=10;i<=NF;i++) c[$i]=i} !/^#/ && /SVTYPE=CNV;/{n=split($9,k,":"); match($8,/CN_FORMS=[^;]*/); split(substr($8,RSTART+9,RLENGTH-9),fs,",")
       for (h in c) if (h=="hgain#1#chr1" || h=="hinv#1#chr1") { split($c[h],v,":"); for(j=1;j<=n;j++) if(k[j]=="CNFORM") r[h]=fs[v[j]] }
       print r["hgain#1#chr1"], r["hinv#1#chr1"]}' "$OUT/sep_short.region.vcf")" "2001-2000:=2001-2200 2361-2360:~2001-2200"
# --diploid: hhalf and heq as one sample; each haplotype's copies are listed, joined as GT is.
sed -e 's/^P hhalf#1#chr1 /P s9#1#chr1 /' -e 's/^P heq#1#chr1 /P s9#2#chr1 /' "$OUT/sep.spec" > "$OUT/sepdip.spec"
python3 -I "$(cd "$(dirname "$0")/.." && pwd)/tools/spec_gfa.py" "$OUT/sepdip.spec" "$OUT/sepdip"
"$BIN" bubble -i "$OUT/sepdip/g.gfa" -r ref --min-variant-bp 0 -o "$OUT/sepdip.b" --merge-nearby-bp 0 --quiet >/dev/null 2>&1
"$BIN" call -i "$OUT/sepdip.b.sorted.gfa" --bubble-prefix "$OUT/sepdip.b" -r ref -o "$OUT/sepdip" --cn --diploid --quiet >/dev/null 2>&1
check "copies, --diploid: a sample's CNFORM joins its haplotypes' as GT does (s9 GT:CNFORM)" \
  "$(awk -F'\t' '/^#CHROM/{for(i=10;i<=NF;i++) if($i=="s9") c=i} !/^#/ && /SVTYPE=CNV;/{n=split($9,k,":"); split($c,v,":")
       for(j=1;j<=n;j++) { if(k[j]=="GT") g=v[j]; if(k[j]=="CNFORM") f=v[j] }; print g":"f}' "$OUT/sepdip.region.vcf")" "1|1:3|3"
# No reference copy: two copies of A B (300 bp) gained by some haplotypes, one by others; the module record
# carries their sequence.
{ printf 'N L 2000 91\nN A 150 92\nN B 150 93\nN C 30 94\nN R 2000 95\nP ref#0#chr1:0-4030 L+,C+,R+\n'
  for i in 1 2 3 4; do printf 'P hone%s#1#chr1 L+,A+,B+,C+,R+\n' "$i"; done
  for i in 1 2; do printf 'P htwo%s#1#chr1 L+,A+,B+,A+,B+,C+,R+\n' "$i"; done
  for i in 1 2 3 4 5 6; do printf 'P hr%s#1#chr1 L+,C+,R+\n' "$i"; done; } > "$OUT/sep0.spec"
python3 -I "$(cd "$(dirname "$0")/.." && pwd)/tools/spec_gfa.py" "$OUT/sep0.spec" "$OUT/sep0"
cn_call "$OUT/sep0/g.gfa" "$OUT/sep0"
check "copies, no reference copy: the gained copies are spelled (CN_FORMS without '*')" \
  "$(awk -F'\t' '!/^#/ && /SVTYPE=CNV;/{print (($8 ~ /CN_FORMS=[^;]*[ACGT]/ && $8 !~ /CN_FORMS=[^;]*\*/) ? "spelled" : "not")}' "$OUT/sep0.region.vcf")" "spelled"
"$BIN" call -i "$OUT/sep0.b.sorted.gfa" --bubble-prefix "$OUT/sep0.b" -r ref -o "$OUT/sep0_short" --cn --max-seq-bp 100 --quiet >/dev/null 2>&1
check "copies, --max-seq-bp under a copy: the form is written without its sequence ('*')" \
  "$(awk -F'\t' '!/^#/ && /SVTYPE=CNV;/{print (($8 !~ /CN_FORMS=/) ? "none" : ($8 ~ /CN_FORMS=[^;]*\*/) ? "unspelled" : "spelled")}' "$OUT/sep0_short.region.vcf")" "unspelled"

# ---------------------------------------------------------------- inverted duplications
# A fold-back walk (A+ A-) adds a reverse-complemented copy. The UNIT route counts it as CN 2 exactly
# like a tandem copy, which is the right dosage but loses the orientation: FORMAT/CNINV and INFO/INVDUP
# record it. An inverted copy held in a separate node is an INS, which --classify-ins flags the same way.
rc_seq() { printf '%s' "$1" | rev | tr ACGT TGCA; }
ID_C=$(rseq 300 33)
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tf0\t%s\nS\tf1\t%s\nS\tf2\t%s\nS\tf3\t%s\n' "$(rseq 1000 31)" "$(rseq 1000 32)" "$(rseq 1000 34)" "$(rseq 1000 35)"
  printf 'S\ta\t%s\nS\tb\t%s\nS\tc\t%s\nS\tcrc\t%s\n' "$(rseq 300 36)" "$(rseq 300 37)" "$ID_C" "$(rc_seq "$ID_C")"
  printf 'L\tf0\t+\ta\t+\t0M\nL\ta\t+\tf1\t+\t0M\nL\ta\t+\ta\t-\t0M\nL\ta\t-\tf1\t+\t0M\n'
  printf 'L\tf1\t+\tb\t+\t0M\nL\tb\t+\tf2\t+\t0M\nL\tb\t+\tb\t+\t0M\n'
  printf 'L\tf2\t+\tc\t+\t0M\nL\tc\t+\tf3\t+\t0M\nL\tc\t+\tcrc\t+\t0M\nL\tcrc\t+\tf3\t+\t0M\n'
  printf 'P\tref#0#chr1:1-4900\tf0+,a+,f1+,b+,f2+,c+,f3+\t*\n'
  printf 'P\thfold#1#chr1\tf0+,a+,a-,f1+,b+,f2+,c+,f3+\t*\n'
  printf 'P\thtan#1#chr1\tf0+,a+,f1+,b+,b+,f2+,c+,f3+\t*\n'
  printf 'P\thinv#1#chr1\tf0+,a+,f1+,b+,f2+,c+,crc+,f3+\t*\n'
} > "$OUT/invdup.gfa"
"$BIN" bubble -i "$OUT/invdup.gfa" -r ref --min-variant-bp 0 -o "$OUT/invdup.b" --merge-nearby-bp 0 --quiet >/dev/null 2>&1
"$BIN" call -i "$OUT/invdup.b.sorted.gfa" --bubble-prefix "$OUT/invdup.b" -r ref -o "$OUT/invdup" \
  --classify-ins --quiet >/dev/null 2>&1
V="$OUT/invdup.region.vcf"
fval() {   # fval <vcf> <record-id-substring> <sample> <FORMAT key>
  awk -F'\t' -v r="$2" -v s="$3" -v k="$4" '/^#CHROM/{for(i=10;i<=NF;i++) if($i==s) c=i}
    !/^#/ && index($3, r) {split($9,kk,":"); split($c,v,":"); for(j in kk) if(kk[j]==k) print v[j]}' "$1"; }
has_flag() { awk -F'\t' -v r="$2" '!/^#/ && index($3, r) { if ($8 ~ /(^|;)INVDUP(;|$)/) print "yes"; else print "no" }' "$1"; }
fold_id=$(awk -F'\t' '!/^#/ && $3 ~ /CNV/ {print $3}' "$V" | while read -r i; do
  [ "$(fval "$V" "$i" 'hfold#1#chr1' GT)" = "1" ] && echo "$i"; done | head -1)
tan_id=$(awk -F'\t' '!/^#/ && $3 ~ /CNV/ {print $3}' "$V" | while read -r i; do
  [ "$(fval "$V" "$i" 'htan#1#chr1' GT)" = "1" ] && echo "$i"; done | head -1)
ins_id=$(awk -F'\t' '!/^#/ && $8 ~ /SVTYPE=INS;/ {print $3}' "$V" | head -1)
[ -n "$fold_id" ] && [ "$(fval "$V" "$fold_id" 'hfold#1#chr1' CN)" = "2" ] \
  && [ "$(fval "$V" "$fold_id" 'hfold#1#chr1' CNINV)" = "1" ] && [ "$(has_flag "$V" "$fold_id")" = "yes" ] \
  && ok "a fold-back copy stays a UNIT DUP (CN 2) and is flagged INVDUP with CNINV 1" \
  || bad "fold-back: record '$fold_id' CN=$(fval "$V" "$fold_id" 'hfold#1#chr1' CN) CNINV=$(fval "$V" "$fold_id" 'hfold#1#chr1' CNINV) INVDUP=$(has_flag "$V" "$fold_id")"
[ -n "$tan_id" ] && [ "$(fval "$V" "$tan_id" 'htan#1#chr1' CNINV)" = "0" ] && [ "$(has_flag "$V" "$tan_id")" = "no" ] \
  && ok "a tandem copy reads CNINV 0 and carries no INVDUP" \
  || bad "tandem: record '$tan_id' CNINV=$(fval "$V" "$tan_id" 'htan#1#chr1' CNINV) INVDUP=$(has_flag "$V" "$tan_id")"
[ -n "$ins_id" ] && grep -v '^#' "$V" | awk -F'\t' -v r="$ins_id" '$3==r' | grep -q 'INS_SUBTYPE=DUP;' \
  && [ "$(has_flag "$V" "$ins_id")" = "yes" ] \
  && ok "--classify-ins flags an inverted copy held in its own node as INS_SUBTYPE=DUP + INVDUP" \
  || bad "--classify-ins did not flag the inverted INS '$ins_id' as an inverted duplication"

# --classify-ins without --cn: hdup inserts a copy of the reference's 400 bp node a with every 25th base
# changed (96% identity), hnew 400 novel bases. The copy is a DUP at the default --ins-dup-min-identity 0.90
# and NOVEL at 0.99; the novel insertion is NOVEL either way.
CI_A=$(rseq 400 41)
CI_D=$(printf '%s' "$CI_A" | awk '{ s = ""; for (i = 1; i <= length($0); i++) { c = substr($0, i, 1)
  if (i % 25 == 0) c = (c == "A" ? "C" : c == "C" ? "G" : c == "G" ? "T" : "A"); s = s c }; print s }')
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tf0\t%s\nS\ta\t%s\nS\tf1\t%s\nS\td\t%s\nS\tn\t%s\n' "$(rseq 1000 42)" "$CI_A" "$(rseq 1000 43)" "$CI_D" "$(rseq 400 44)"
  printf 'L\tf0\t+\ta\t+\t0M\nL\ta\t+\tf1\t+\t0M\nL\ta\t+\td\t+\t0M\nL\td\t+\tf1\t+\t0M\nL\ta\t+\tn\t+\t0M\nL\tn\t+\tf1\t+\t0M\n'
  printf 'P\tref#0#chr1:1-2400\tf0+,a+,f1+\t*\nP\thdup#1#chr1\tf0+,a+,d+,f1+\t*\nP\thnew#1#chr1\tf0+,a+,n+,f1+\t*\n'
} > "$OUT/cins.gfa"
"$BIN" bubble -i "$OUT/cins.gfa" -r ref --min-variant-bp 0 -o "$OUT/cins.b" --merge-nearby-bp 0 --quiet >/dev/null 2>&1
ci_sub() {   # ci_sub <vcf> <sample>: INS_SUBTYPE of the INS that sample carries
  awk -F'\t' -v s="$2" '/^#CHROM/{for(i=10;i<=NF;i++) if($i==s) c=i}
    !/^#/ && $8 ~ /SVTYPE=INS;/ && substr($c, 1, 1) == "1" { match($8, /INS_SUBTYPE=[A-Z]+/); print substr($8, RSTART + 12, RLENGTH - 12) }' "$1"; }
for id in 0.90 0.99; do
  "$BIN" call -i "$OUT/cins.b.sorted.gfa" --bubble-prefix "$OUT/cins.b" -r ref -o "$OUT/cins_$id" \
    --classify-ins --ins-dup-min-identity "$id" --minimap-preset asm20 --minimap-best-n 4 --quiet >/dev/null 2>&1
done
check "--classify-ins without --cn: the 96% copy at identity 0.90 and 0.99, the novel insertion at both" \
  "$(ci_sub "$OUT/cins_0.90.region.vcf" 'hdup#1#chr1') $(ci_sub "$OUT/cins_0.99.region.vcf" 'hdup#1#chr1') $(ci_sub "$OUT/cins_0.90.region.vcf" 'hnew#1#chr1') $(ci_sub "$OUT/cins_0.99.region.vcf" 'hnew#1#chr1')" \
  "DUP NOVEL NOVEL NOVEL"
refuses "--minimap-best-n 0 is refused" "must be >= 1" \
  "$BIN" call -i "$OUT/cins.b.sorted.gfa" --bubble-prefix "$OUT/cins.b" -r ref -o "$OUT/cins_bad" --classify-ins --minimap-best-n 0 -q

# ---------------------------------------------------------------- UNIT DUP placement
# A UNIT record describes one node, so it is placed on that node. It used to take the bubble's source
# as its flank, which put it upstream of the array whenever other sequence sat in between -- on the
# real loci always, the LPA KIV-2 array by 4.5 kb. With no reference copy the record spans no
# reference sequence, so END = POS.
rep_place_gfa() {   # rep_place_gfa <out.gfa> <ref copies> <hap>:<copies> ...
  local out=$1 rc=$2; shift 2
  reps() { local k=$1 s=""; for ((i=0;i<k;i++)); do s="$s,r+"; done; printf '%s' "$s"; }
  { printf 'H\tVN:Z:1.0\n'
    printf 'S\tn1\t%s\nS\tp\t%s\nS\tr\t%s\nS\tn3\t%s\n' "$(rseq 1000 41)" "$(rseq 200 42)" "$(rseq 60 43)" "$(rseq 1000 44)"
    printf 'L\tn1\t+\tp\t+\t0M\nL\tp\t+\tr\t+\t0M\nL\tr\t+\tr\t+\t0M\nL\tr\t+\tn3\t+\t0M\nL\tp\t+\tn3\t+\t0M\nL\tn1\t+\tr\t+\t0M\n'
    printf 'P\tref#0#chr1:1-%d\tn1+,p+%s,n3+\t*\n' $((2200 + rc * 60)) "$(reps "$rc")"
    for h in "$@"; do printf 'P\t%s#1#chr1\tn1+,p+%s,n3+\t*\n' "${h%%:*}" "$(reps "${h#*:}")"; done
    # deletes p, so the bubble opens at n1 and p sits between its source and the array
    printf 'P\thdp#1#chr1\tn1+%s,n3+\t*\n' "$(reps 3)"
  } > "$out"
}
rep_pos_end() { grep -v '^#' "$1" | awk -F'\t' '/CN_METHOD=UNIT/{match($8,/END=[0-9]+/); print $2":"substr($8,RSTART+4,RLENGTH-4); exit}'; }
rep_place_gfa "$OUT/rpl.gfa" 4 h2:2 h7:7
rep_call "$OUT/rpl.gfa" "$OUT/rpl"
[ "$(rep_pos_end "$OUT/rpl.region.vcf")" = "1200:1440" ] \
  && ok "a UNIT DUP covers its array (POS 1200, END 1440: four 60 bp copies), not the bubble's source" \
  || bad "UNIT DUP placed at POS:END $(rep_pos_end "$OUT/rpl.region.vcf"), expected 1200:1440"
rep_place_gfa "$OUT/rpl0.gfa" 0 h3:3 h5:5
rep_call "$OUT/rpl0.gfa" "$OUT/rpl0"
[ "$(rep_pos_end "$OUT/rpl0.region.vcf")" = "1200:1200" ] \
  && ok "with no reference copy it sits at the insertion point and spans no reference (END = POS)" \
  || bad "no-reference-copy UNIT DUP at POS:END $(rep_pos_end "$OUT/rpl0.region.vcf"), expected 1200:1200"

# ---------------------------------------------------------------- collapsed-module summary and routing
# Reference S P Q P Q K: the module {P,Q} is visited twice, 300 bp a copy. htri carries three copies;
# hhalf carries P Q P, two passes through the anchor P but a copy short of Q: CN 2, no carrier, and the
# missing Q a 150 bp deletion beside the module record.
MOD_HDR='bubble_id,source,source_orient,sink,sink_orient,inside_node_count,total_node_count,path_support,distinct_alleles,ref_allele_support,alt_allele_support_max,alt_allele_support_min,min_inside_bp,max_inside_bp,inside_nodes'
mod_gfa() {   # mod_gfa <out.gfa> <path line>...  -- the module graph above with the given haplotypes
  local out=$1; shift
  local -a haps=("$@")   # the link loop below reuses the positional parameters
  { printf 'H\tVN:Z:1.0\n'
    printf 'S\tF\t%s\nS\tS\t%s\nS\tP\t%s\nS\tQ\t%s\nS\tX\t%s\nS\tK\t%s\nS\tG\t%s\n' "$(rseq 500 51)" \
      "$(rseq 40 52)" "$(rseq 150 53)" "$(rseq 150 54)" "$(rseq 30 55)" "$(rseq 40 56)" "$(rseq 500 57)"
    for e in "F S" "S P" "P Q" "Q P" "Q K" "P K" "P X" "X P" "K G"; do
      set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
    printf 'P\tref#0#chr1:1-1680\tF+,S+,P+,Q+,P+,Q+,K+,G+\t*\n'
    for h in "${haps[@]}"; do printf 'P\t%s\t*\n' "$h"; done; } > "$out"
}
mod_gfa "$OUT/mdec.gfa" $'htri#1#chr1\tF+,S+,P+,Q+,P+,Q+,P+,Q+,K+,G+' $'hhalf#1#chr1\tF+,S+,P+,Q+,P+,K+,G+'
printf '%s\n1,S,+,K,+,2,4,3,3,1,1,1,0,900,"P;Q"\n' "$MOD_HDR" > "$OUT/mdec.bubbles.csv"
"$BIN" call -i "$OUT/mdec.gfa" -c "$OUT/mdec.bubbles.csv" -r ref -o "$OUT/mdec" --cn -q >/dev/null 2>&1
check "a partial copy: CN counts whole copies at the anchor, the missing part is a DEL (records: SVTYPE SVLEN GT:CN)" \
  "$(grep -v '^#' "$OUT/mdec.region.vcf" | awk -F'\t' '{ match($8, /SVTYPE=[A-Z]+(;SVLEN=-?[0-9]+)?/); g = ""
      for (i = 10; i <= NF; i++) { split($i, f, ":"); g = g " " f[1] (f[2] == "." || f[2] == "" ? "" : ":" f[2]) }
      printf "%s%s | ", substr($8, RSTART, RLENGTH), g }')" \
  "SVTYPE=CNV 0:2 1:3 0:2 | SVTYPE=DEL;SVLEN=-150 0 0 1 | "

# hpeak walks P X P Q P: three passes through P but the module's 600 bp, no change in size, so no carrier
# and no record. The busiest-node count (the PEAK route) read it as a gain.
mod_gfa "$OUT/minv.gfa" $'hpeak#1#chr1\tF+,S+,P+,X+,P+,Q+,P+,K+,G+'
printf '%s\n1,S,+,K,+,3,5,2,2,1,1,1,0,900,"P;Q;X"\n' "$MOD_HDR" > "$OUT/minv.bubbles.csv"
"$BIN" call -i "$OUT/minv.gfa" -c "$OUT/minv.bubbles.csv" -r ref -o "$OUT/minv" --cn -q >/dev/null 2>&1
! grep -v '^#' "$OUT/minv.region.vcf" | grep -q 'SVTYPE=CNV' \
  && ok "a module whose bp does not change has no carrier and no DUP" \
  || bad "CN-invariant module: expected no DUP record, got $(grep -v '^#' "$OUT/minv.region.vcf" | grep -o 'CN_METHOD=[A-Z_]*' | head -1)"

# ---------------------------------------------------------------- a path crossing a bubble twice
# htwice runs S A K, back to S, then S B K: three traversals (S..K, the direct K->S step, S..K). Its
# walk is the widest source..sink interval, so it includes the first crossing and the K->S return:
# the INS it yields is S,K,B rather than B. Not corrected, but every record of the bubble must say
# so, in both VCFs; hdel traverses once.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%s\nS\tA\t%s\nS\tB\t%s\nS\tK\t%s\nS\tG\t%s\n' "$(rseq 500 61)" "$(rseq 40 62)" \
    "$(rseq 100 63)" "$(rseq 120 64)" "$(rseq 40 65)" "$(rseq 500 66)"
  for e in "F S" "S A" "A K" "S B" "B K" "K S" "S K" "K G"; do
    set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref#0#chr1:1-1180\tF+,S+,A+,K+,G+\t*\n'
  printf 'P\thdel#1#chr1\tF+,S+,K+,G+\t*\n'
  printf 'P\thtwice#1#chr1\tF+,S+,A+,K+,S+,B+,K+,G+\t*\n'; } > "$OUT/mc.gfa"
printf '%s\n1,S,+,K,+,2,4,3,3,1,1,1,0,120,"A;B"\n' "$MOD_HDR" > "$OUT/mc.bubbles.csv"
mc_log=$("$BIN" call -i "$OUT/mc.gfa" -c "$OUT/mc.bubbles.csv" -r ref -o "$OUT/mc" 2>&1 >/dev/null)
mc_n=$(grep -v '^#' "$OUT/mc.region.vcf" | wc -l | tr -d ' ')
mc_f=$(grep -v '^#' "$OUT/mc.region.vcf" | grep -c 'MULTI_CROSSING=1[^0-9]')
[ "$mc_n" -ge 1 ] && [ "$mc_f" = "$mc_n" ] \
  && grep -q '^##INFO=<ID=MULTI_CROSSING,' "$OUT/mc.region.vcf" \
  && printf '%s\n' "$mc_log" | grep -q 'paths crossing a site more than once: 1; sites: 1' \
  && ok "a path crossing a bubble twice flags all $mc_n record(s) MULTI_CROSSING=1" \
  || bad "double crossing: $mc_f of $mc_n record(s) carry MULTI_CROSSING=1"
! grep -v '^#' "$VCF" | grep -q 'MULTI_CROSSING' \
  && ok "MULTI_CROSSING is absent where every path crosses once" \
  || bad "MULTI_CROSSING emitted on the single-crossing coordinate fixture"

# ---------------------------------------------------------------- bubble reversed against the reference
# The same graph called with its bubble row in reference order (fwd) and with source and sink swapped
# (rev): a CSV whose orientation disagrees with call's -r path, as when -r names a path other than the
# one the bubbles were oriented to. Walks are spelled source->sink, so rev wrote DELSEQ/INVSEQ/INSSEQ
# reverse-complemented and put the INS one base late (the first base of the downstream flank). Every
# one of those records must now be the same line in both runs. The substitution's INS arm is excluded: it anchors
# on the walk-order flank, the far end of its DEL arm in a reverse bubble.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\nS\t5\t%s\nS\t6\t%s\n' "$(rseq 40 71)" "$(rseq 60 72)" \
    "$(rseq 70 73)" "$(rseq 40 74)" "$(rseq 25 75)" "$(rseq 55 76)"
  for e in "1 + 2 +" "2 + 4 +" "4 + 5 +" "2 + 3 +" "3 + 4 +" "1 + 4 +" "1 + 2 -" "2 - 4 +" "1 + 6 +" "6 + 4 +"; do
    set -- $e; printf 'L\t%s\t%s\t%s\t%s\t0M\n' "$1" "$2" "$3" "$4"; done
  printf 'P\tref#0#chr1:1001-1165\t1+,2+,4+,5+\t*\n'
  printf 'P\thINS#1#c1\t1+,2+,3+,4+,5+\t*\n'       # INS of node 3 after node 2
  printf 'P\thDEL#1#c1\t1+,4+,5+\t*\n'             # DEL of node 2
  printf 'P\thINV#1#c1\t1+,2-,4+,5+\t*\n'          # INV of node 2
  printf 'P\thSUB#1#c1\t1+,6+,4+,5+\t*\n'; } > "$OUT/rv.gfa"   # substitution 2 -> 6
printf '%s\n1,1,+,4,+,3,5,5,5,1,1,1,0,130,"2;3;6"\n' "$MOD_HDR" > "$OUT/rv_fwd.csv"
printf '%s\n1,4,+,1,+,3,5,5,5,1,1,1,0,130,"2;3;6"\n' "$MOD_HDR" > "$OUT/rv_rev.csv"
for o in fwd rev; do
  "$BIN" call -i "$OUT/rv.gfa" -c "$OUT/rv_$o.csv" -r ref -o "$OUT/rv_$o" --min-sv-bp 20 -q >/dev/null 2>&1
done
rv_line() { grep -v '^#' "$1" | awk -F'\t' -v id="$2" '$3==id'; }
for id in bubble1_DEL_2 bubble1_INV_2 bubble1_INS_3; do
  f=$(rv_line "$OUT/rv_fwd.region.vcf" "$id"); r=$(rv_line "$OUT/rv_rev.region.vcf" "$id")
  [ -n "$f" ] && [ "$f" = "$r" ] \
    && ok "reversed bubble: $id is the same record (POS, END, sequence) as in reference order" \
    || bad "reversed bubble: $id differs -- fwd '$(printf '%s' "$f" | cut -f2-5)' rev '$(printf '%s' "$r" | cut -f2-5)'"
done
[ "$(rv_line "$OUT/rv_fwd.region.vcf" bubble1_INS_3 | cut -f2)" = "1100" ] \
  && ok "the INS sits after the last base of node 2 (POS 1100)" \
  || bad "INS POS $(rv_line "$OUT/rv_fwd.region.vcf" bubble1_INS_3 | cut -f2), expected 1100"

# ---------------------------------------------------------------- traversals vs boundary visits
# MULTI_CROSSING counts traversals -- boundary to the OTHER boundary through interior nodes only --
# not boundary visits. htouchS deletes A, then reaches S again from outside the bubble (X -> S) and
# leaves it for Y; htouchK meets K before ever entering the bubble. Each traverses S..K once.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%s\nS\tA\t%s\nS\tK\t%s\nS\tX\t%s\nS\tY\t%s\nS\tG\t%s\n' "$(rseq 500 61)" "$(rseq 40 62)" \
    "$(rseq 100 63)" "$(rseq 40 65)" "$(rseq 80 67)" "$(rseq 90 68)" "$(rseq 500 66)"
  for e in "F S" "S A" "A K" "S K" "K G" "K X" "X S" "S Y" "Y G" "F K"; do
    set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref#0#chr1:1-1180\tF+,S+,A+,K+,G+\t*\n'
  printf 'P\thdel#1#chr1\tF+,S+,K+,G+\t*\n'
  printf 'P\thtouchS#1#chr1\tF+,S+,K+,X+,S+,Y+,G+\t*\n'
  printf 'P\thtouchK#1#chr1\tF+,K+,X+,S+,K+,G+\t*\n'; } > "$OUT/mct.gfa"
printf '%s\n1,S,+,K,+,1,3,4,2,1,3,3,0,100,"A"\n' "$MOD_HDR" > "$OUT/mct.bubbles.csv"
mct_log=$("$BIN" call -i "$OUT/mct.gfa" -c "$OUT/mct.bubbles.csv" -r ref -o "$OUT/mct" 2>&1 >/dev/null)
[ "$(grep -v '^#' "$OUT/mct.region.vcf" | awk -F'\t' '$3=="bubble1_DEL_A"' | grep -o 'AC=[0-9]*')" = "AC=3" ] \
  && ok "boundary visits outside the walk: the DEL of A is still called on all three carriers" \
  || bad "touch-only fixture: DEL of A not called with AC=3"
! grep -v '^#' "$OUT/mct.region.vcf" | grep -q 'MULTI_CROSSING' \
  && ! printf '%s\n' "$mct_log" | grep -q 'more than once' \
  && ok "a path touching a boundary without traversing the bubble is not flagged MULTI_CROSSING (VCF and summary)" \
  || bad "touch-only paths flagged: $(grep -v '^#' "$OUT/mct.region.vcf" | grep -o 'MULTI_CROSSING=[0-9]*' | head -1)"

# Traversals are counted over the whole path, not only its walk. hback crosses S..K with A deleted,
# leaves for X, and comes back K- A- S-: its walk is one crossing, the other lies outside it. hempty
# crosses S..K twice with nothing between, leaving for X in between. Both are flagged; hdel is not.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%s\nS\tA\t%s\nS\tK\t%s\nS\tX\t%s\nS\tG\t%s\n' "$(rseq 500 61)" "$(rseq 40 62)" \
    "$(rseq 100 63)" "$(rseq 40 65)" "$(rseq 80 67)" "$(rseq 500 66)"
  for e in "F + S +" "S + A +" "A + K +" "S + K +" "K + G +" "K + X +" "X + S +" "X + K -" "S - G +"; do
    set -- $e; printf 'L\t%s\t%s\t%s\t%s\t0M\n' "$1" "$2" "$3" "$4"; done
  printf 'P\tref#0#chr1:1-1180\tF+,S+,A+,K+,G+\t*\n'
  printf 'P\thdel#1#chr1\tF+,S+,K+,G+\t*\n'
  printf 'P\thback#1#chr1\tF+,S+,K+,X+,K-,A-,S-,G+\t*\n'
  printf 'P\themptyx#1#chr1\tF+,S+,K+,X+,S+,K+,G+\t*\n'; } > "$OUT/mcw.gfa"
printf '%s\n1,S,+,K,+,1,3,4,2,1,3,3,0,100,"A"\n' "$MOD_HDR" > "$OUT/mcw.bubbles.csv"
mcw_log=$("$BIN" call -i "$OUT/mcw.gfa" -c "$OUT/mcw.bubbles.csv" -r ref -o "$OUT/mcw" 2>&1 >/dev/null)
mcw_n=$(grep -v '^#' "$OUT/mcw.region.vcf" | wc -l | tr -d ' ')
mcw_f=$(grep -v '^#' "$OUT/mcw.region.vcf" | grep -c 'MULTI_CROSSING=2[^0-9]')
[ "$mcw_n" -ge 1 ] && [ "$mcw_f" = "$mcw_n" ] \
  && printf '%s\n' "$mcw_log" | grep -q 'paths crossing a site more than once: 2; sites: 1' \
  && ok "a second traversal outside the walk, and two empty traversals, flag MULTI_CROSSING=2 on all $mcw_n record(s)" \
  || bad "traversals outside the walk: $mcw_f of $mcw_n record(s) carry MULTI_CROSSING=2"

# ---------------------------------------------------------------- orientation from the reference walk
# The reference meets the sink K once upstream (a repeat copy), then crosses the bubble S..K forward:
#
#   F 500  1..500   K 40  501..540   X 80  541..620   S 40  621..660   A 100  661..760   K 40  761..800
#
# hdel deletes A; hins inserts B (70 bp) between A and K. Orientation used to be read from the
# boundaries' FIRST occurrences (K at 501 < S at 621), so this forward bubble was treated as reversed and
# the INS sat at 661. Every record must be what the forward bubble says.
RK_S=$(rseq 40 62); RK_A=$(rseq 100 63); RK_B=$(rseq 70 64)
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%s\nS\tA\t%s\nS\tB\t%s\nS\tK\t%s\nS\tX\t%s\nS\tG\t%s\n' "$(rseq 500 61)" "$RK_S" \
    "$RK_A" "$RK_B" "$(rseq 40 65)" "$(rseq 80 67)" "$(rseq 500 66)"
  for e in "F K" "K X" "X S" "S A" "A K" "S K" "K G" "A B" "B K"; do
    set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref#0#chr1:1-1300\tF+,K+,X+,S+,A+,K+,G+\t*\n'
  printf 'P\thdel#1#chr1\tF+,K+,X+,S+,K+,G+\t*\n'
  printf 'P\thins#1#chr1\tF+,K+,X+,S+,A+,B+,K+,G+\t*\n'; } > "$OUT/rk.gfa"
printf '%s\n1,S,+,K,+,2,4,3,3,1,1,1,0,100,"A;B"\n' "$MOD_HDR" > "$OUT/rk.bubbles.csv"
"$BIN" call -i "$OUT/rk.gfa" -c "$OUT/rk.bubbles.csv" -r ref -o "$OUT/rk" -q >/dev/null 2>&1
rk_get() {   # rk_get <vcf> <id> <column> [INFO key]
  awk -F'\t' -v id="$2" -v c="$3" -v k="${4:-}" '!/^#/ && $3==id { if (k=="") { print $c; exit }
    n=split($8,a,";"); for(i=1;i<=n;i++){split(a[i],kv,"="); if(kv[1]==k){print kv[2]; exit}} }' "$1"; }
seqform() { [ "$1" = "$2" ] && { echo forward; return; }   # seqform <got> <want>
            [ -n "$1" ] && [ "$1" = "$(rc_seq "$2")" ] && { echo reverse-complemented; return; }
            echo "other:${1:0:24}"; }
V="$OUT/rk.region.vcf"
check "forward bubble, sink repeated upstream: DEL POS:END"  "$(rk_get "$V" bubble1_DEL_A 2):$(rk_get "$V" bubble1_DEL_A 8 END)" "660:760"
check "forward bubble, sink repeated upstream: DELSEQ (node A)" "$(seqform "$(rk_get "$V" bubble1_DEL_A 8 DELSEQ)" "$RK_A")" forward
check "forward bubble, sink repeated upstream: INS POS:END:REF" \
      "$(rk_get "$V" bubble1_INS_B 2):$(rk_get "$V" bubble1_INS_B 8 END):$(rk_get "$V" bubble1_INS_B 4)" "760:760:${RK_A: -1}"
check "forward bubble, sink repeated upstream: INSSEQ (node B)" "$(seqform "$(rk_get "$V" bubble1_INS_B 8 INSSEQ)" "$RK_B")" forward
! grep -v '^#' "$OUT/rk.region.vcf" | grep -q 'MULTI_CROSSING' \
  && ok "the reference's own upstream visit to K is no traversal: no MULTI_CROSSING" \
  || bad "the orientation fixture is flagged MULTI_CROSSING"

# The module route's POS/END under the same layout: K also sits upstream of S, so the module span is
# S(641..680) .. K(1281..1340) with P Q P Q between. Which boundary the span opens on was decided by
# first occurrences (K first), so POS/END were shifted by the length difference of K and S (20 bp).
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tK\t%s\nS\tX\t%s\nS\tS\t%s\nS\tP\t%s\nS\tQ\t%s\nS\tG\t%s\n' "$(rseq 500 51)" "$(rseq 60 56)" \
    "$(rseq 80 58)" "$(rseq 40 52)" "$(rseq 150 53)" "$(rseq 150 54)" "$(rseq 500 57)"
  for e in "F K" "K X" "X S" "S P" "P Q" "Q P" "Q K" "K G"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref#0#chr1:1-1840\tF+,K+,X+,S+,P+,Q+,P+,Q+,K+,G+\t*\n'
  printf 'P\thtri#1#chr1\tF+,K+,X+,S+,P+,Q+,P+,Q+,P+,Q+,K+,G+\t*\n'; } > "$OUT/rkmod.gfa"
printf '%s\n1,S,+,K,+,2,4,2,2,1,1,1,0,900,"P;Q"\n' "$MOD_HDR" > "$OUT/rkmod.bubbles.csv"
"$BIN" call -i "$OUT/rkmod.gfa" -c "$OUT/rkmod.bubbles.csv" -r ref -o "$OUT/rkmod" --cn -q >/dev/null 2>&1
V="$OUT/rkmod.region.vcf"
check "module DUP, sink repeated upstream: POS:END:CN_MODULE_REF_BP" \
      "$(rk_get "$V" bubble1_CNV_S 2):$(rk_get "$V" bubble1_CNV_S 8 END):$(rk_get "$V" bubble1_CNV_S 8 CN_MODULE_REF_BP)" "680:1280:600"

# ---------------------------------------------------------------- a reference on the minus strand
# The reference walks the bubble S..K as G- K- A- S- F-, i.e. sink to source on '-' steps:
#
#   G 500  1..500   K- 40  501..540   A- 100  541..640   S- 40  641..680   F 500  681..1180
#
# The DEL of A is placed on the reference walk: POS 540, the last reference base of K-, and DELSEQ the
# reference's own bases. hdel deletes A; hins inserts B between A and K.
RM_K="A$(rseq 38 65)A"; RM_A=$(rseq 100 63); RM_B=$(rseq 70 64)
RM_F=$(rseq 500 61); RM_S=$(rseq 40 62); RM_G=$(rseq 500 66)
RM_REF=$(rc_seq "$RM_F$RM_S$RM_A$RM_K$RM_G")
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%s\nS\tA\t%s\nS\tB\t%s\nS\tK\t%s\nS\tG\t%s\n' "$RM_F" "$RM_S" "$RM_A" "$RM_B" "$RM_K" "$RM_G"
  for e in "F S" "S A" "A K" "S K" "K G" "A B" "B K"; do
    set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref#0#chr1:1-1180\tG-,K-,A-,S-,F-\t*\n'
  printf 'P\thdel#1#chr1\tG-,K-,S-,F-\t*\n'
  printf 'P\thins#1#chr1\tG-,K-,B-,A-,S-,F-\t*\n'; } > "$OUT/rm.gfa"
printf '%s\n1,S,+,K,+,2,4,3,3,1,1,1,0,100,"A;B"\n' "$MOD_HDR" > "$OUT/rm.bubbles.csv"
"$BIN" call -i "$OUT/rm.gfa" -c "$OUT/rm.bubbles.csv" -r ref -o "$OUT/rm" -q >/dev/null 2>&1
V="$OUT/rm.region.vcf"
check "minus-strand reference: DEL POS:END:DELSEQ is the reference span" \
      "$(rk_get "$V" bubble1_DEL_A 2):$(rk_get "$V" bubble1_DEL_A 8 END):$(rk_get "$V" bubble1_DEL_A 8 DELSEQ)" \
      "540:640:${RM_REF:540:100}"

# ---------------------------------------------------------------- UNIT placement from the reference walk
# A UNIT DUP is placed on the bubble's own reference occurrences, not on a node's first occurrence on
# the whole reference, which can be a copy outside the bubble. Node lengths: F 500, U 80, S 40, A 100,
# R 60 (the self-loop unit), K 40, D 90, G 500.
#   rep_up    ref F R U S A R R K G: R also sits upstream at 501..560; the bubble's copies are at
#             781..840 and 841..900, so POS:END is 780:840 (it was 500:560)
#   rep0_src  ref F S U S A K G, carriers insert R right after the bubble's source S (621..660): no
#             reference copy, so POS = END = 660 (the upstream S at 501..540 gave 540)
#   rep0_rev  ref G- S- D- K- A- S- F-, minus strand, source also downstream: carriers insert R between
#             A- (671..770) and the bubble's S- (771..810), so POS = END = 770 (it was 540, then 500)
#   rep0_skip ref F S A D A K G, carriers skip the first A and insert R after the second (731..830):
#             POS = END = 830, which counting visits to A got wrong (640)
#   rep0_rskip the same on the minus strand, ref G- K- A- D- A- S- F-: carriers insert R between K-
#             (501..540) and the first A- while skipping the second, so POS = END = 540
#   rep0_loop ref F S A K G, carriers run A twice (A is looped too) and then R: the copies follow
#             A (541..640), not S, although A is not aligned, so POS = END = 640
rep_gfa() {   # rep_gfa <name> <ref walk> <hap walk>...
  local name=$1 i=0 n p; shift
  { printf 'H\tVN:Z:1.0\n'
    printf 'S\tF\t%s\nS\tU\t%s\nS\tS\t%s\nS\tA\t%s\nS\tR\t%s\nS\tK\t%s\nS\tD\t%s\nS\tG\t%s\n' "$(rseq 500 41)" \
      "$(rseq 80 42)" "$(rseq 40 43)" "$(rseq 100 44)" "$(rseq 60 45)" "$(rseq 40 46)" "$(rseq 90 47)" "$(rseq 500 48)"
    for p in "$@"; do printf '%s\n' "$p" | tr ',' '\n' | awk 'NR > 1 { print prev "\t" $0 } { prev = $0 }'; done \
      | sort -u | awk -F'\t' '{ printf "L\t%s\t%s\t%s\t%s\t0M\n", substr($1, 1, length($1) - 1), substr($1, length($1)),
                                 substr($2, 1, length($2) - 1), substr($2, length($2)) }'
    for p in "$@"; do
      if [ $i = 0 ]; then n='ref#0#chr1'; else n="h$i#1#chr1"; fi
      printf 'P\t%s\t%s\t*\n' "$n" "$p"; i=$((i + 1))
    done; } > "$OUT/$name.gfa"
  printf '%s\n1,S,+,K,+,2,4,3,3,1,1,1,0,900,"A;R"\n' "$MOD_HDR" > "$OUT/$name.bubbles.csv"
  "$BIN" call -i "$OUT/$name.gfa" -c "$OUT/$name.bubbles.csv" -r ref -o "$OUT/$name" --cn -q >/dev/null 2>&1
}
rep_gfa rep_up 'F+,R+,U+,S+,A+,R+,R+,K+,G+' 'F+,R+,U+,S+,A+,R+,R+,R+,K+,G+' 'F+,R+,U+,S+,A+,R+,K+,G+' \
  'F+,R+,U+,S+,A+,R+,R+,R+,R+,K+,G+'
rep_gfa rep0_src 'F+,S+,U+,S+,A+,K+,G+' 'F+,S+,U+,S+,R+,R+,A+,K+,G+' 'F+,S+,U+,S+,R+,R+,R+,A+,K+,G+' \
  'F+,S+,U+,S+,A+,K+,G+'
rep_gfa rep0_rev 'G-,S-,D-,K-,A-,S-,F-' 'G-,S-,D-,K-,A-,R-,R-,S-,F-' 'G-,S-,D-,K-,A-,R-,R-,R-,S-,F-' \
  'G-,S-,D-,K-,A-,S-,F-'
rep_gfa rep0_skip 'F+,S+,A+,D+,A+,K+,G+' 'F+,S+,D+,A+,R+,R+,K+,G+' 'F+,S+,D+,A+,R+,R+,R+,K+,G+' \
  'F+,S+,A+,D+,A+,K+,G+'
rep_gfa rep0_rskip 'G-,K-,A-,D-,A-,S-,F-' 'G-,K-,R-,R-,A-,D-,S-,F-' 'G-,K-,R-,R-,R-,A-,D-,S-,F-' \
  'G-,K-,A-,D-,A-,S-,F-'
rep_gfa rep0_loop 'F+,S+,A+,K+,G+' 'F+,S+,A+,A+,R+,R+,K+,G+' 'F+,S+,A+,A+,R+,R+,R+,K+,G+' 'F+,S+,A+,K+,G+'
for c in "rep_up 780:900:2" "rep0_src 660:660:0" "rep0_rev 770:770:0" "rep0_skip 830:830:0" \
         "rep0_rskip 540:540:0" "rep0_loop 640:640:0"; do
  set -- $c; V="$OUT/$1.region.vcf"
  check "UNIT DUP from the reference walk ($1): POS:END:REF_CN" \
        "$(rk_get "$V" bubble1_CNV_R 2):$(rk_get "$V" bubble1_CNV_R 8 END):$(rk_get "$V" bubble1_CNV_R 8 REF_CN)" "$2"
done

# ---------------------------------------------------------------- a copy change counted once (de-dup)
# A 300 bp module (ma, mb) the reference runs twice, then a unique anchor k, a 100 bp node w and a 200 bp
# node y. The
# MODULE DUP counts copies; the walk diff can also type an extra or missing copy as INS/DEL. A carrier
# leaves an INS/DEL only when its own copy change accounts for it:
#   hgain   three copies: the copy INS restates its CN gain -> no INS record carries it
#   hgainx  three copies plus a novel 200 bp node x: the copy INS goes, the novel INS stays
#   hnovel  two copies plus x: CN-neutral, so its INS of x stays
#   hloss   one copy and w deleted: the copy DEL goes (CN 1 says it); the 100 bp DEL of w stays
#   hgain3  five copies: a 900 bp INS, three units, still restates the gain -> gone. The old rule kept
#           anything beyond twice the unit, which double-counted multi-copy gains
#   hneut   two copies, y deleted: a CN-neutral 200 bp DEL stays. The old rule dropped every unlinked
#           DEL of half a unit or more in a CN bubble, copy change or not
dd_walk() { printf 'f0+,%s,k+,%s,f1+' "$1" "$2"; }
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tf0\t%s\nS\tf1\t%s\nS\tma\t%s\nS\tmb\t%s\nS\tk\t%s\nS\tw\t%s\nS\tx\t%s\nS\ty\t%s\n' "$(rseq 1000 71)" \
    "$(rseq 2000 72)" "$(rseq 150 73)" "$(rseq 150 74)" "$(rseq 60 75)" "$(rseq 100 76)" "$(rseq 200 77)" "$(rseq 200 78)"
  M="ma+,mb+"
  W_ref=$(dd_walk "$M,$M" "w+,y+"); W_gain=$(dd_walk "$M,$M,$M" "w+,y+"); W_gainx=$(dd_walk "$M,$M,$M,x+" "w+,y+")
  W_novel=$(dd_walk "$M,$M,x+" "w+,y+"); W_loss=$(dd_walk "$M" "y+"); W_gain3=$(dd_walk "$M,$M,$M,$M,$M" "w+,y+")
  W_neut=$(dd_walk "$M,$M" "w+")
  for p in "$W_ref" "$W_gain" "$W_gainx" "$W_novel" "$W_loss" "$W_gain3" "$W_neut"; do printf '%s\n' "$p" | tr ',' '\n' | awk 'NR > 1 { print prev "\t" $0 } { prev = $0 }'; done \
    | sort -u | awk -F'\t' '{ printf "L\t%s\t%s\t%s\t%s\t0M\n", substr($1, 1, length($1) - 1), substr($1, length($1)),
                               substr($2, 1, length($2) - 1), substr($2, length($2)) }'
  printf 'P\tref#0#chr1:1-3960\t%s\t*\n' "$W_ref"
  printf 'P\thgain#1#chr1\t%s\t*\nP\thgainx#1#chr1\t%s\t*\nP\thnovel#1#chr1\t%s\t*\nP\thloss#1#chr1\t%s\t*\n' \
    "$W_gain" "$W_gainx" "$W_novel" "$W_loss"
  printf 'P\thgain3#1#chr1\t%s\t*\nP\thneut#1#chr1\t%s\t*\n' "$W_gain3" "$W_neut"
} > "$OUT/dedup.gfa"
# one bubble f0..f1 over all of it, so the 100 bp DEL sits in the same CN bubble as the copy loss
printf '%s\n1,f0,+,f1,+,6,8,7,7,1,1,1,0,2000,"ma;mb;k;w;x;y"\n' "$MOD_HDR" > "$OUT/dedup.bubbles.csv"
"$BIN" call -i "$OUT/dedup.gfa" -c "$OUT/dedup.bubbles.csv" -r ref -o "$OUT/dedup" --cn -q >/dev/null 2>&1
V="$OUT/dedup.region.vcf"
dd_carriers() {   # dd_carriers <SVTYPE> <SVLEN>: carriers of that record, sorted, comma-joined
  awk -F'\t' -v t="$1" -v l="$2" '/^#CHROM/{for(i=10;i<=NF;i++) n[i]=$i}
    !/^#/ && $8 ~ ("SVTYPE=" t ";") && $8 ~ ("SVLEN=" l ";") {
      for(i=10;i<=NF;i++){split($i,g,":"); if(g[1]!="0" && g[1]!=".") {sub(/#.*/,"",n[i]); print n[i]}}}' "$V" | sort | paste -sd, -; }
check "de-dup: the module DUP stands (MODULE records)" "$(grep -v '^#' "$V" | grep -c 'CN_METHOD=MODULE;' | tr -d ' ')" "1"
check "de-dup: no carrier keeps a copy-sized INS (300 bp) whose CN gain explains it" "$(dd_carriers INS 300)" ""
check "de-dup: the novel 200 bp INS keeps both its carriers, CN-neutral and gainer" "$(dd_carriers INS 200)" "hgainx,hnovel"
check "de-dup: the copy-loss DEL (300 bp) is gone" "$(dd_carriers DEL -300)" ""
check "de-dup: the 100 bp DEL beside the copy loss stays" "$(dd_carriers DEL -100)" "hloss"
check "de-dup: a multi-copy restatement (900 bp INS, +3 copies) is gone" "$(dd_carriers INS 900)" ""
check "de-dup: a CN-neutral DEL larger than half a unit (200 bp) stays" "$(dd_carriers DEL -200)" "hneut"
# A copy gain (h1 repeats X B, 55 bp) and a 45 bp novel insertion (h2-h4 take N, which spells the start of
# X B). Before de-dup the INS took h1's 55 bp; once h1 is taken off it, its size is 45 bp, under --min-sv-bp,
# and its size filter must be recomputed: the INS goes, the CNV stays.
FBX=$(rseq 30 31); FBB=$(rseq 25 32); FBN=$(printf '%s%s' "$FBX" "$FBB" | cut -c1-45)
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%s\nS\tA\t%s\nS\tX\t%s\nS\tB\t%s\nS\tN\t%s\nS\tK\t%s\nS\tG\t%s\n' "$(rseq 300 1)" \
    "$(rseq 40 2)" "$(rseq 40 3)" "$FBX" "$FBB" "$FBN" "$(rseq 40 4)" "$(rseq 300 5)"
  for e in "F S" "S A" "A X" "X B" "B K" "B X" "B N" "N K" "K G"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref#0#chr1:1-775\tF+,S+,A+,X+,B+,K+,G+\t*\nP\th1#1#chr1\tF+,S+,A+,X+,B+,X+,B+,K+,G+\t*\n'
  for h in h2 h3 h4; do printf 'P\t%s#1#chr1\tF+,S+,A+,X+,B+,N+,K+,G+\t*\n' $h; done
  printf 'P\th5#1#chr1\tF+,S+,A+,X+,B+,K+,G+\t*\n'; } > "$OUT/fb.gfa"
"$BIN" bubble -i "$OUT/fb.gfa" -r ref --min-variant-bp 0 -o "$OUT/fb.b" --merge-nearby-bp 0 -q >/dev/null 2>&1
"$BIN" call -i "$OUT/fb.b.sorted.gfa" --bubble-prefix "$OUT/fb.b" -r ref -o "$OUT/fb" --cn -q >/dev/null 2>&1
check "de-dup: an INS left under --min-sv-bp once a copy gain is taken off it is dropped (records)" \
  "$(grep -v '^#' "$OUT/fb.region.vcf" | grep -o 'SVTYPE=[A-Z]*' | cut -c8- | paste -sd, -)" "CNV"
# A copy lost whose remaining copy ends in another version of B (Z, walked twice by htwo, so a module node).
# The reference runs X a1 B X a2 B p; hflav runs X a1 p Z. The walk comparison cuts it at p: a 470 bp DEL
# (B X a2 B), then Z inserted after p. Joined, they are the copy lost (-330); apart, the DEL alone was taken
# for it and Z was left as a 140 bp INS beside the copy loss.
sf_walks="ref:f0+,X+,a1+,B+,X+,a2+,B+,p+,k+,f1+ hflav:f0+,X+,a1+,p+,Z+,k+,f1+
  htwo:f0+,X+,a1+,Z+,X+,a2+,Z+,p+,k+,f1+ hloss:f0+,X+,a1+,B+,p+,k+,f1+"
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tf0\t%s\nS\tf1\t%s\nS\tX\t%s\nS\ta1\t%s\nS\ta2\t%s\nS\tB\t%s\nS\tp\t%s\nS\tZ\t%s\nS\tk\t%s\n' "$(rseq 1000 61)" \
    "$(rseq 1000 62)" "$(rseq 150 63)" "$(rseq 40 64)" "$(rseq 40 65)" "$(rseq 140 66)" "$(rseq 10 67)" "$(rseq 140 68)" "$(rseq 60 69)"
  for w in $sf_walks; do printf '%s\n' "${w#*:}" | tr ',' '\n' | awk 'NR > 1 { print prev "\t" $0 } { prev = $0 }'; done \
    | sort -u | awk -F'\t' '{ printf "L\t%s\t%s\t%s\t%s\t0M\n", substr($1, 1, length($1) - 1), substr($1, length($1)),
                               substr($2, 1, length($2) - 1), substr($2, length($2)) }'
  for w in $sf_walks; do
    n=${w%%:*}; [ "$n" = ref ] && n='ref#0#chr1:1-2730' || n="$n#1#chr1"
    printf 'P\t%s\t%s\t*\n' "$n" "${w#*:}"
  done
} > "$OUT/subflav.gfa"
printf '%s\n1,f0,+,f1,+,7,9,4,4,1,1,1,0,2000,"X;a1;a2;B;p;Z;k"\n' "$MOD_HDR" > "$OUT/subflav.bubbles.csv"
"$BIN" call -i "$OUT/subflav.gfa" -c "$OUT/subflav.bubbles.csv" -r ref -o "$OUT/subflav" --cn -q >/dev/null 2>&1
V="$OUT/subflav.region.vcf"
check "de-dup: a copy loss cut in two by a short match takes both pieces (hflav CN; its INS records)" \
  "$(awk -F'\t' '/^#CHROM/{for(i=10;i<=NF;i++) if($i ~ /^hflav/) c=i} !/^#/ && /CN_METHOD=MODULE;/{split($c,g,":"); print g[2]}' "$V"):$(dd_carriers INS 140)" "1:"
# The same copy loss with the pieces kept apart by a 60 bp match (p), and Z an 86% version of B: hvers's
# remaining copy holds another version of what it lost, so its INS of Z is part of the copy change. hnew
# inserts a novel N there instead, which stays. hfrag inserts a 60 bp version of part of B: a version, but
# far short of what its deletion took beyond the copy lost, so it stays too.
vers() { printf '%s' "$1" | awk '{n=split($0,c,""); m["A"]="C"; m["C"]="G"; m["G"]="T"; m["T"]="A"
  for(i=1;i<=n;i++) printf "%s", (i % 7 == 3 ? m[c[i]] : c[i])}'; }
VB=$(rseq 140 66)
sv_walks="ref:f0+,X+,a1+,B+,X+,a2+,B+,p+,k+,f1+ hvers:f0+,X+,a1+,p+,Z+,k+,f1+
  htwo:f0+,X+,a1+,Z+,X+,a2+,Z+,p+,k+,f1+ hnew:f0+,X+,a1+,p+,N+,k+,f1+ hloss:f0+,X+,a1+,B+,p+,k+,f1+
  hfrag:f0+,X+,a1+,p+,Q+,k+,f1+"
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tf0\t%s\nS\tf1\t%s\nS\tX\t%s\nS\ta1\t%s\nS\ta2\t%s\nS\tB\t%s\nS\tp\t%s\nS\tZ\t%s\nS\tN\t%s\nS\tQ\t%s\nS\tk\t%s\n' \
    "$(rseq 1000 61)" "$(rseq 1000 62)" "$(rseq 150 63)" "$(rseq 40 64)" "$(rseq 40 65)" "$VB" "$(rseq 60 67)" "$(vers "$VB")" \
    "$(rseq 140 68)" "$(vers "${VB:0:60}")" "$(rseq 60 69)"
  for w in $sv_walks; do printf '%s\n' "${w#*:}" | tr ',' '\n' | awk 'NR > 1 { print prev "\t" $0 } { prev = $0 }'; done \
    | sort -u | awk -F'\t' '{ printf "L\t%s\t%s\t%s\t%s\t0M\n", substr($1, 1, length($1) - 1), substr($1, length($1)),
                               substr($2, 1, length($2) - 1), substr($2, length($2)) }'
  for w in $sv_walks; do
    n=${w%%:*}; [ "$n" = ref ] && n='ref#0#chr1:1-2780' || n="$n#1#chr1"
    printf 'P\t%s\t%s\t*\n' "$n" "${w#*:}"
  done
} > "$OUT/subvers.gfa"
printf '%s\n1,f0,+,f1,+,9,11,6,6,1,1,1,0,2000,"X;a1;a2;B;p;Z;N;Q;k"\n' "$MOD_HDR" > "$OUT/subvers.bubbles.csv"
"$BIN" call -i "$OUT/subvers.gfa" -c "$OUT/subvers.bubbles.csv" -r ref -o "$OUT/subvers" --cn -q >/dev/null 2>&1
V="$OUT/subvers.region.vcf"
check "de-dup: an INS beside a copy loss that is a version of the lost copy goes; a novel one stays (INS 140 carriers)" \
  "$(dd_carriers INS 140)" "hnew"
check "de-dup: a short version beside a copy loss, far under its excess, stays (INS 60 carriers)" "$(dd_carriers INS 60)" "hfrag"

# ---------------------------------------------------------------- module copy scale and partial copies
# The reference visits a 1000 bp node P twice and a 50 bp node Q three times (Q P Q P Q), then a unique
# anchor k and a 600 bp node R once. REF_CN is 2 (P, the anchor, revisited by every path holding two
# copies and the most bp), not the 3 one short node reaches.
#   hgain  one more P Q: CN 3
#   hloss  one P Q fewer: CN 1
#   hpart  one P Q fewer AND R deleted: still CN 1 (whole copies); R's 600 bp DEL stays, since no copy
#          change explains it
ms_walk() { printf 'f0+,%s,k+,%s,f1+' "$1" "$2"; }
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tf0\t%s\nS\tf1\t%s\nS\tP\t%s\nS\tQ\t%s\nS\tk\t%s\nS\tR\t%s\n' "$(rseq 1000 81)" "$(rseq 2000 82)" \
    "$(rseq 1000 83)" "$(rseq 50 84)" "$(rseq 60 85)" "$(rseq 600 86)"
  W_ref=$(ms_walk "Q+,P+,Q+,P+,Q+" "R+"); W_gain=$(ms_walk "Q+,P+,Q+,P+,Q+,P+,Q+" "R+")
  W_loss=$(ms_walk "Q+,P+,Q+" "R+"); W_part="f0+,Q+,P+,Q+,k+,f1+"
  for p in "$W_ref" "$W_gain" "$W_loss" "$W_part"; do printf '%s\n' "$p" | tr ',' '\n' | awk 'NR > 1 { print prev "\t" $0 } { prev = $0 }'; done \
    | sort -u | awk -F'\t' '{ printf "L\t%s\t%s\t%s\t%s\t0M\n", substr($1, 1, length($1) - 1), substr($1, length($1)),
                               substr($2, 1, length($2) - 1), substr($2, length($2)) }'
  printf 'P\tref#0#chr1:1-5810\t%s\t*\nP\thgain#1#chr1\t%s\t*\nP\thloss#1#chr1\t%s\t*\nP\thpart#1#chr1\t%s\t*\n' \
    "$W_ref" "$W_gain" "$W_loss" "$W_part"
} > "$OUT/mscale.gfa"
printf '%s\n1,f0,+,f1,+,4,6,4,4,1,1,1,0,9000,"P;Q;k;R"\n' "$MOD_HDR" > "$OUT/mscale.bubbles.csv"
"$BIN" call -i "$OUT/mscale.gfa" -c "$OUT/mscale.bubbles.csv" -r ref -o "$OUT/mscale" --cn -q >/dev/null 2>&1
V="$OUT/mscale.region.vcf"
ms_fmt() { awk -F'\t' -v s="$1" -v key="$2" '/^#CHROM/{for(i=10;i<=NF;i++) if($i==s) c=i}
  !/^#/ && /CN_METHOD=MODULE;/{n=split($9,k,":"); split($c,v,":"); for(j=1;j<=n;j++) if(k[j]==key) print v[j]}' "$V"; }
check "module scale: REF_CN from passes through the anchor (REF_CN:CN_METHOD)" \
      "$(rk_get "$V" "$(grep -v '^#' "$V" | awk -F'\t' '/CN_METHOD=MODULE;/{print $3; exit}')" 8 REF_CN):$(grep -v '^#' "$V" | grep -o 'CN_METHOD=[A-Z]*' | head -1 | cut -d= -f2)" \
      "2:MODULE"
check "module scale: CN of the gain, the loss and the partial loss" \
      "$(ms_fmt 'hgain#1#chr1' CN):$(ms_fmt 'hloss#1#chr1' CN):$(ms_fmt 'hpart#1#chr1' CN)" "3:1:1"
check "partial copy: the 600 bp DEL no copy change explains stays" \
      "$(grep -v '^#' "$V" | grep -c 'SVTYPE=DEL;SVLEN=-600;' | tr -d ' ')" "1"

# ---------------------------------------------------------------- twin-aware UNIT counting
# One element on two nodes: L carries a self-loop, T spells the same 30 bp, and both sit between the
# same neighbours -- the pair refine's twin merge folds into one node. The reference has the element
# twice, both on T (F S T A T K G), so counting L alone reads every path routed through L as copies
# the reference lacks -- the C4 22 bp record. Counted together, the reference has 2 copies: h1 (L
# twice) is CN 2 and no carrier, h2 (L four times) the one gain, h3 (T once) a sub-threshold loss,
# and the record sits on T's first reference copy (541..570).
#   tw_move   the element moved past A (ref F S T A K G, haps F S A L K G): T and L share no
#             neighbours, so they are two places, not twins, and the move stays a DEL plus a REF_CN 0 DUP
#   tw_peak   --cn, ref T A T and haps L A L L: +30 bp, below threshold. UNIT stands down, and the module, which
#             counts the twin with its loop too, must not bring back a record
#   tw_ord1/2 no reference copy; one carrier walks only T, the other L: the record's POS must not
#             depend on which carrier comes first in the GFA
#   tw_loops  both nodes looped (h2 walks T T, h1 L L L): refine merges them, so call counts them as one
#             unit too -- one record on L, REF_CN 2, h1 the only carrier, and none on T
twin_gfa() {   # twin_gfa <name> <element bp> <call options> <ref walk> <hap walk>...
  local name=$1 el opts=$3 i=0 n p; el=$(rseq "$2" 61); shift 3
  { printf 'H\tVN:Z:1.0\n'
    printf 'S\tF\t%s\nS\tS\t%s\nS\tA\t%s\nS\tL\t%s\nS\tT\t%s\nS\tK\t%s\nS\tG\t%s\n' "$(rseq 500 62)" \
      "$(rseq 40 63)" "$(rseq 100 64)" "$el" "$el" "$(rseq 40 65)" "$(rseq 500 66)"
    for p in "$@"; do printf '%s\n' "$p" | tr ',' '\n' | awk 'NR > 1 { print prev "\t" $0 } { prev = $0 }'; done \
      | sort -u | awk -F'\t' '{ printf "L\t%s\t%s\t%s\t%s\t0M\n", substr($1, 1, length($1) - 1), substr($1, length($1)),
                                 substr($2, 1, length($2) - 1), substr($2, length($2)) }'
    for p in "$@"; do
      if [ $i = 0 ]; then n='ref#0#chr1'; else n="h$i#1#chr1"; fi
      printf 'P\t%s\t%s\t*\n' "$n" "$p"; i=$((i + 1))
    done; } > "$OUT/$name.gfa"
  printf '%s\n1,S,+,K,+,3,5,4,4,1,1,1,0,900,"A;L;T"\n' "$MOD_HDR" > "$OUT/$name.bubbles.csv"
  "$BIN" call -i "$OUT/$name.gfa" -c "$OUT/$name.bubbles.csv" -r ref -o "$OUT/$name" -q $opts >/dev/null 2>&1
}
twin_gfa tw30 30 '' 'F+,S+,T+,A+,T+,K+,G+' 'F+,S+,L+,A+,L+,K+,G+' 'F+,S+,L+,A+,L+,L+,L+,K+,G+' 'F+,S+,T+,A+,K+,G+'
twin_gfa tw_move 60 '' 'F+,S+,T+,A+,K+,G+' 'F+,S+,A+,L+,K+,G+' 'F+,S+,A+,L+,L+,K+,G+' 'F+,S+,T+,A+,K+,G+'
twin_gfa tw_peak 30 --cn 'F+,S+,T+,A+,T+,K+,G+' 'F+,S+,L+,A+,L+,L+,K+,G+' 'F+,S+,L+,A+,L+,L+,K+,G+' 'F+,S+,T+,A+,T+,K+,G+'
twin_gfa tw_ord1 30 '' 'F+,S+,A+,K+,G+' 'F+,S+,T+,A+,T+,K+,G+' 'F+,S+,L+,A+,L+,L+,K+,G+'
twin_gfa tw_ord2 30 '' 'F+,S+,A+,K+,G+' 'F+,S+,L+,A+,L+,L+,K+,G+' 'F+,S+,T+,A+,T+,K+,G+'
twin_gfa tw_loops 30 '' 'F+,S+,T+,A+,T+,K+,G+' 'F+,S+,L+,A+,L+,L+,L+,K+,G+' 'F+,S+,T+,T+,A+,T+,K+,G+'
V="$OUT/tw30.region.vcf"
check "twin UNIT: L and T counted as one unit (REF_CN:AC:UNIT_TWIN_NODES)" \
      "$(rk_get "$V" bubble1_CNV_L 8 REF_CN):$(rk_get "$V" bubble1_CNV_L 8 AC):$(rk_get "$V" bubble1_CNV_L 8 UNIT_TWIN_NODES)" "2:1:T"
check "twin UNIT: covers the reference's copies, the first via the twin, and A between them (POS:END)" \
      "$(rk_get "$V" bubble1_CNV_L 2):$(rk_get "$V" bubble1_CNV_L 8 END)" "540:700"
check "twin UNIT: a moved element is not a twin (REF_CN, twin field, DEL_T kept)" \
      "$(rk_get "$OUT/tw_move.region.vcf" bubble1_CNV_L 8 REF_CN):$(grep -c UNIT_TWIN_NODES= "$OUT/tw_move.region.vcf" | tr -d ' '):$(grep -c 'bubble1_DEL_T' "$OUT/tw_move.region.vcf" | tr -d ' ')" "0:0:1"
check "twin UNIT: under --cn, the module counts the twin too and writes no record" \
      "$(grep -v '^#' "$OUT/tw_peak.region.vcf" | grep -c '_CNV_' | tr -d ' ')" "0"
check "twin UNIT: two looped twins are one unit (REF_CN:AC:twins, DUP records)" \
      "$(rk_get "$OUT/tw_loops.region.vcf" bubble1_CNV_L 8 REF_CN):$(rk_get "$OUT/tw_loops.region.vcf" bubble1_CNV_L 8 AC):$(rk_get "$OUT/tw_loops.region.vcf" bubble1_CNV_L 8 UNIT_TWIN_NODES):$(grep -v '^#' "$OUT/tw_loops.region.vcf" | grep -c '_CNV_' | tr -d ' ')" "2:1:T:1"
check "twin UNIT: a REF_CN 0 record's POS does not depend on carrier order" \
      "$(rk_get "$OUT/tw_ord1.region.vcf" bubble1_CNV_L 2):$(rk_get "$OUT/tw_ord1.region.vcf" bubble1_CNV_L 8 AC)" \
      "$(rk_get "$OUT/tw_ord2.region.vcf" bubble1_CNV_L 2):$(rk_get "$OUT/tw_ord2.region.vcf" bubble1_CNV_L 8 AC)"

# ---------------------------------------------------------------- loop nodes in the DEL/INS walks
# A loop node leaves the DEL/INS/INV walks only where a surviving record counts its copies.
#   w1_sub   C4 after refine: a 22 bp looped element sits once in each copy of a module (A L), the
#            reference has two copies and no haplotype changes it by 50 bp, so no UNIT record. hdel
#            loses a whole copy: the DEL must be the full 122 bp and precise; dropping the loop node
#            made it 100 bp
#   w1_mod   --cn, T twin of L in each copy of a module (T A): the module counts every copy, so nothing
#            restates the copies as INS/DEL; the loop is inside the copies and stays in the walks, so h3's
#            two extra L in one copy are a 60 bp INS
#   w1_route T twin of L, no record (a +30 bp change at most): h1 runs L where the reference runs T,
#            the same sequence, so it is no variant -- not DEL T plus INS L
twin_gfa w1_sub 22 '' 'F+,S+,A+,L+,A+,L+,K+,G+' 'F+,S+,A+,L+,L+,A+,L+,K+,G+' 'F+,S+,A+,L+,K+,G+' 'F+,S+,A+,L+,A+,L+,K+,G+'
twin_gfa w1_mod 30 --cn 'F+,S+,T+,A+,T+,A+,K+,G+' 'F+,S+,L+,A+,L+,A+,K+,G+' 'F+,S+,T+,A+,T+,A+,T+,A+,K+,G+' \
  'F+,S+,L+,L+,L+,A+,L+,A+,K+,G+' 'F+,S+,T+,A+,K+,G+'
twin_gfa w1_route 30 '' 'F+,S+,T+,A+,T+,K+,G+' 'F+,S+,L+,A+,L+,K+,G+' 'F+,S+,L+,L+,A+,L+,K+,G+' 'F+,S+,T+,A+,T+,K+,G+'
V="$OUT/w1_sub.region.vcf"
check "kept loop: no UNIT record, and the lost copy is one precise 122 bp DEL (DUPs:SVLEN:IMPRECISE)" \
      "$(grep -v '^#' "$V" | grep -c '_CNV_' | tr -d ' '):$(grep -v '^#' "$V" | grep -o 'SVTYPE=DEL;SVLEN=-[0-9]*' | sed 's/.*=//' | paste -sd, -):$(grep -v '^#' "$V" | grep -c IMPRECISE | tr -d ' ')" "0:-122:0"
V="$OUT/w1_mod.region.vcf"
check "--cn module: no INS/DEL restates its copies; one copy's +60 bp of L is an INS (MODULE:other records)" \
      "$(grep -v '^#' "$V" | grep -c 'CN_METHOD=MODULE;' | tr -d ' '):$(grep -v '^#' "$V" | grep -v '_CNV_' | grep -o 'SVTYPE=[A-Z]*;SVLEN=-\?[0-9]*' | paste -sd, -)" \
      "1:SVTYPE=INS;SVLEN=60"
# h1 runs both module copies through L, the twin of the reference's T: its copies are whole copies.
check "--cn module: a copy routed through the twin is a whole copy (h1 GT:CN:CNBP:CNRESID:MODBP:CNFORM)" \
      "$(awk -F'\t' '/^#CHROM/{for(i=10;i<=NF;i++) if($i=="h1#1#chr1") c=i} !/^#/ && /CN_METHOD=MODULE;/{print $c}' "$V")" "0:2:0:0:0:."
check "a path routed through the twin of a kept loop is no variant" \
      "$(grep -vc '^#' "$OUT/w1_route.region.vcf" | tr -d ' ')" "0"

# ---------------------------------------------------------------- a UNIT record beside a CN-invariant module
# Reference S P Q T P Q K: {P,Q} folds twice (the module applies) and T is a 12 bp self-loop. hrep
# runs T seven times: every haplotype has the module's copy number, but the short-unit UNIT route
# still writes a DUP. The summary must not then claim the bubble has no CN record.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%s\nS\tP\t%s\nS\tQ\t%s\nS\tT\t%s\nS\tK\t%s\nS\tG\t%s\n' "$(rseq 500 51)" \
    "$(rseq 40 52)" "$(rseq 150 53)" "$(rseq 150 54)" "$(rseq 12 55)" "$(rseq 40 56)" "$(rseq 500 57)"
  for e in "F S" "S P" "P Q" "Q T" "T T" "T P" "Q K" "K G"; do
    set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref#0#chr1:1-1692\tF+,S+,P+,Q+,T+,P+,Q+,K+,G+\t*\n'
  printf 'P\threp#1#chr1\tF+,S+,P+,Q+,T+,T+,T+,T+,T+,T+,T+,P+,Q+,K+,G+\t*\n'; } > "$OUT/mrep.gfa"
printf '%s\n1,S,+,K,+,3,5,2,2,1,1,1,0,900,"P;Q;T"\n' "$MOD_HDR" > "$OUT/mrep.bubbles.csv"
mrep_log=$("$BIN" call -i "$OUT/mrep.gfa" -c "$OUT/mrep.bubbles.csv" -r ref -o "$OUT/mrep" --cn 2>&1 >/dev/null)
grep -v '^#' "$OUT/mrep.region.vcf" | grep -q 'CN_METHOD=UNIT;' \
  && ! printf '%s\n' "$mrep_log" | grep -q 'every haplotype has the reference' \
  && ok "a UNIT DUP beside a CN-invariant module: the bubble is not counted as having no CN record" \
  || bad "CN-invariant module with a UNIT record: UNIT=$(grep -vc '^#' "$OUT/mrep.region.vcf") summary='$(printf '%s\n' "$mrep_log" | grep 'every haplotype')'"

# ---------------------------------------------------------------- a deletion at the region's first base
# Reference S A S B K, bubble S..K, S being the region's first node. hdel walks S B K, so the deletion
# is S A: the first 160 bases, with no preceding base to anchor on. POS is then the first affected
# base and END the last, i.e. END = POS + |SVLEN| - 1; it used to be one base past the deletion.
C8_S=$(rseq 60 80)
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tS\t%s\nS\tA\t%s\nS\tB\t%s\nS\tK\t%s\nS\tG\t%s\n' "$C8_S" "$(rseq 100 81)" "$(rseq 70 84)" \
    "$(rseq 40 82)" "$(rseq 500 83)"
  for e in "S A" "A S" "S B" "B K" "K G"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref#0#chr1:1-830\tS+,A+,S+,B+,K+,G+\t*\n'
  printf 'P\thdel#1#chr1\tS+,B+,K+,G+\t*\n'; } > "$OUT/first.gfa"
printf '%s\n1,S,+,K,+,2,4,2,2,1,1,1,0,170,"A;B"\n' "$MOD_HDR" > "$OUT/first.bubbles.csv"
# The node comparison's placement (the reference revisits S, so the aligned deletion is the same
# haplotype one base later, anchored on base 1 as VCF normalisation leaves it).
"$BIN" call -i "$OUT/first.gfa" -c "$OUT/first.bubbles.csv" -r ref -o "$OUT/first" --min-sv-bp 20 --no-align-alleles -q >/dev/null 2>&1
VCF="$OUT/first.region.vcf"
check "DEL at the region's first base: POS:END:SVLEN" "$(field DEL 2):$(info DEL END):$(info DEL SVLEN)" "1:160:-160"
check "DEL at the region's first base: REF is that first base" "$(field DEL 4)" "${C8_S:0:1}"
grep -v '^#' "$VCF" | grep -q 'IMPRECISE' && bad "DEL at the region's first base is marked IMPRECISE" \
  || ok "DEL at the region's first base is PRECISE (END - POS = |SVLEN| - 1, by the convention)"
VCF="$OUT/cc.region.vcf"

# ---------------------------------------------------------------- at most one module-route DUP per bubble
# A bubble's module copy number is one record; two report its CN twice. Checked over every region VCF the
# --cn fixtures above wrote.
mod_n=0; mod_dup=""
for v in "$OUT"/*.region.vcf; do
  r=$(grep -v '^#' "$v" | awk -F'\t' '$8 ~ /CN_METHOD=MODULE;/ { match($8, /BUBBLE_ID=[0-9]+/);
        n[substr($8, RSTART, RLENGTH)]++; t++ } END { m=0; for (b in n) if (n[b] > 1) m++; print t+0, m }')
  mod_n=$((mod_n + ${r% *})); [ "${r#* }" = "0" ] || mod_dup="$mod_dup $(basename "$v")"
done
[ "$mod_n" -ge 2 ] && [ -z "$mod_dup" ] \
  && ok "no bubble carries more than one MODULE DUP ($mod_n such records checked)" \
  || bad "module-route DUPs: $mod_n record(s); bubbles with more than one in:${mod_dup:- (none, but too few records to check)}"

# ---------------------------------------------------------------- assembly gaps and letter case
# Reference F S A K G; hsub replaces A (100 bp) by B (160 bp), a substitution written as an INS of its
# net 60 bp; hlow by b, which is B written in lower case on a node of its own, and hgap by a 60 bp run
# of N (an assembly gap). The gap haplotype's allele is unknown: GT '.', out of AN, and no record of N
# bases. hlow's allele is hsub's: one INS of AC 2, spelled in upper case.
GB=$(rseq 160 74)
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%s\nS\tA\t%s\nS\tB\t%s\nS\tb\t%s\nS\tN\t%s\nS\tK\t%s\nS\tG\t%s\n' \
    "$(rseq 300 71)" "$(rseq 40 72)" "$(rseq 100 73)" "$GB" "$(printf '%s' "$GB" | tr 'ACGT' 'acgt')" \
    "$(printf 'N%.0s' $(seq 60))" "$(rseq 40 75)" "$(rseq 300 76)"
  for e in "F S" "S A" "A K" "S B" "B K" "S b" "b K" "S N" "N K" "K G"; do
    set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref#0#chr1:1-780\tF+,S+,A+,K+,G+\t*\n'
  printf 'P\thsub#1#chr1\tF+,S+,B+,K+,G+\t*\n'
  printf 'P\thlow#1#chr1\tF+,S+,b+,K+,G+\t*\n'
  printf 'P\thgap#1#chr1\tF+,S+,N+,K+,G+\t*\n'; } > "$OUT/gap.gfa"
printf '%s\n1,S,+,K,+,4,6,4,4,1,1,1,0,160,"A;B;b;N"\n' "$MOD_HDR" > "$OUT/gap.bubbles.csv"
gap_log=$("$BIN" call -i "$OUT/gap.gfa" -c "$OUT/gap.bubbles.csv" -r ref -o "$OUT/gap" 2>&1 >/dev/null)
gt_of() { awk -F'\t' -v s="$2" -v t="$3" '/^#CHROM/{for(i=10;i<=NF;i++) if($i ~ s) c=i}
            c && !/^#/ && $8 ~ ("SVTYPE="t";") {split($c,g,":"); print g[1]; exit}' "$1"; }   # no column: nothing
VCF="$OUT/gap.region.vcf"
check "a haplotype crossing an assembly gap: GT on the substitution, typed INS by its net change" \
  "$(gt_of "$VCF" hgap INS):$(info INS SVLEN):$(grep -v '^#' "$VCF" | grep -c SUBSTITUTION)" ".:60:1"
check "the gap haplotype is out of AN; the case variant joins its allele (INS AN:AC)" \
  "$(info INS AN):$(info INS AC)" "2:2"
# Aligned, the replacement is trimmed of the bases A and B share at either end, so ALT, after its base at
# POS, is part of B.
ins=$(field INS 5 | cut -c2-)
ins_ok=no
[ -n "$ins" ] && [ "$ins" = "$(printf '%s' "$ins" | tr a-z A-Z)" ] && [[ "$GB" == *"$ins"* ]] && ins_ok=yes
check "the inserted bases (ALT) are upper case, and B's" "$ins_ok" "yes"
grep -v '^#' "$VCF" | grep -q 'NNNNNNNNNN' && bad "a record spells the assembly gap's N bases" \
  || ok "no record spells the gap's N bases"
printf '%s\n' "$gap_log" | grep -q 'GT . across an assembly gap: 1' \
  && ok "the run summary counts the haplotype-site set missing" \
  || bad "the run summary does not report the assembly gap"
VCF="$OUT/cc.region.vcf"

# ---------------------------------------------------------------- coalescing in a reverse bubble
# Reference F S A M C K G; hcut deletes A and C, two 40 bp pieces 30 bp apart. Each is under
# --min-sv-bp 50; coalesced they are one 80 bp IMPRECISE deletion. The bubble given as (K,S), against
# the reference, walks the pieces in decreasing coordinate order, which the gap test read as
# -(40+30+40) bp apart and never joined. Coalescing belongs to the node comparison (--no-align-alleles):
# aligned at base level they are two 40 bp deletions, under --min-sv-bp.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%s\nS\tA\t%s\nS\tM\t%s\nS\tC\t%s\nS\tK\t%s\nS\tG\t%s\n' "$(rseq 300 81)" \
    "$(rseq 40 82)" "$(rseq 40 83)" "$(rseq 30 84)" "$(rseq 40 85)" "$(rseq 40 86)" "$(rseq 300 87)"
  for e in "F S" "S A" "A M" "M C" "C K" "S M" "M K" "K G"; do
    set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref#0#chr1:1-790\tF+,S+,A+,M+,C+,K+,G+\t*\n'
  printf 'P\thcut#1#chr1\tF+,S+,M+,K+,G+\t*\n'; } > "$OUT/coal.gfa"
printf '%s\n1,S,+,K,+,3,5,2,2,1,1,1,0,110,"A;M;C"\n' "$MOD_HDR" > "$OUT/coal_fwd.bubbles.csv"
printf '%s\n1,K,+,S,+,3,5,2,2,1,1,1,0,110,"A;M;C"\n' "$MOD_HDR" > "$OUT/coal_rev.bubbles.csv"
for t in fwd rev; do
  "$BIN" call -i "$OUT/coal.gfa" -c "$OUT/coal_$t.bubbles.csv" -r ref -o "$OUT/coal_$t" --no-align-alleles -q >/dev/null 2>&1
  VCF="$OUT/coal_$t.region.vcf"
  check "two deletion pieces 30 bp apart coalesce ($t bubble): POS:END:SVLEN" \
    "$(field DEL 2):$(info DEL END):$(info DEL SVLEN)" "340:450:-80"
done
VCF="$OUT/cc.region.vcf"

# ---------------------------------------------------------------- an insertion beside a counted array
# Reference F S U U K G with U a 60 bp loop unit; hrep carries U x4 (a UNIT gain) and an unrelated
# 80 bp insertion I between the array and K. The array's copies leave the walks, so I is anchored
# next to the array: at its start (POS 340, the UNIT record's POS) in either bubble orientation. A
# reverse bubble used to put it at the array's end, 120 bp away.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%s\nS\tU\t%s\nS\tI\t%s\nS\tK\t%s\nS\tG\t%s\n' "$(rseq 300 91)" \
    "$(rseq 40 92)" "$(rseq 60 93)" "$(rseq 80 94)" "$(rseq 40 95)" "$(rseq 300 96)"
  for e in "F S" "S U" "U U" "U K" "U I" "I K" "K G"; do
    set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref#0#chr1:1-800\tF+,S+,U+,U+,K+,G+\t*\n'
  printf 'P\threp#1#chr1\tF+,S+,U+,U+,U+,U+,I+,K+,G+\t*\n'; } > "$OUT/arr.gfa"
printf '%s\n1,S,+,K,+,2,4,2,2,1,1,1,0,320,"U;I"\n' "$MOD_HDR" > "$OUT/arr_fwd.bubbles.csv"
printf '%s\n1,K,+,S,+,2,4,2,2,1,1,1,0,320,"U;I"\n' "$MOD_HDR" > "$OUT/arr_rev.bubbles.csv"
for t in fwd rev; do
  "$BIN" call -i "$OUT/arr.gfa" -c "$OUT/arr_$t.bubbles.csv" -r ref -o "$OUT/arr_$t" -q >/dev/null 2>&1
  VCF="$OUT/arr_$t.region.vcf"
  check "an insertion beside a UNIT array sits at the array's start ($t bubble): INS POS:UNIT POS" \
    "$(field INS 2):$(field CNV 2)" "340:340"
done
VCF="$OUT/cc.region.vcf"

# ---------------------------------------------------------------- module records only claim what they keep
# The module graph above with a 60 bp node Y after S carrying a self-link that every path walks once.
# The loop never varies, so it has no record and must not switch the module route off: htri's three
# copies are a MODULE DUP. It used to block the route on length alone.
modY_gfa() {   # modY_gfa <out.gfa> <Y-loop: 1|0> <path line>...
  local out=$1 yloop=$2; shift 2
  local -a haps=("$@")
  { printf 'H\tVN:Z:1.0\n'
    printf 'S\tF\t%s\nS\tS\t%s\nS\tY\t%s\nS\tP\t%s\nS\tQ\t%s\nS\tZ\t%s\nS\tK\t%s\nS\tG\t%s\n' "$(rseq 500 51)" \
      "$(rseq 40 52)" "$(rseq 60 58)" "$(rseq 150 53)" "$(rseq 150 54)" "$(rseq 30 59)" "$(rseq 40 56)" "$(rseq 500 57)"
    for e in "F S" "S Y" "Y P" "P Q" "Q P" "Q Z" "Z Z" "Z K" "K G"; do
      set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
    [ "$yloop" = 1 ] && printf 'L\tY\t+\tY\t+\t0M\n'
    printf 'P\tref#0#chr1:1-1770\tF+,S+,Y+,P+,Q+,P+,Q+,Z+,K+,G+\t*\n'
    for h in "${haps[@]}"; do printf 'P\t%s\t*\n' "$h"; done; } > "$out"
}
modY_gfa "$OUT/my.gfa" 1 $'htri#1#chr1\tF+,S+,Y+,P+,Q+,P+,Q+,P+,Q+,Z+,K+,G+'
printf '%s\n1,S,+,K,+,4,6,2,2,1,1,1,0,1000,"Y;P;Q;Z"\n' "$MOD_HDR" > "$OUT/my.bubbles.csv"
"$BIN" call -i "$OUT/my.gfa" -c "$OUT/my.bubbles.csv" -r ref -o "$OUT/my" --cn -q >/dev/null 2>&1
check "an SV-sized loop every path walks once does not block the module (htri CN_METHOD)" \
  "$(grep -v '^#' "$OUT/my.region.vcf" | grep -o 'CN_METHOD=[A-Z_]*' | head -1)" "CN_METHOD=MODULE"
# The 30 bp loop Z: hz1 and hz2 walk it 10 times (a 270 bp UNIT gain) at the reference's module copy
# number; htri gains a module copy alone. With --min-haplotypes 2 the module record (1 carrier)
# is not kept, so it must not claim the bubble: the UNIT record (2 carriers) is written. Claiming
# first and dropping later left neither.
zz='Z+,Z+,Z+,Z+,Z+,Z+,Z+,Z+,Z+,Z+'
modY_gfa "$OUT/mz.gfa" 0 $'htri#1#chr1\tF+,S+,Y+,P+,Q+,P+,Q+,P+,Q+,Z+,K+,G+' \
  "hz1#1#chr1	F+,S+,Y+,P+,Q+,P+,Q+,$zz,K+,G+" "hz2#1#chr1	F+,S+,Y+,P+,Q+,P+,Q+,$zz,K+,G+"
printf '%s\n1,S,+,K,+,4,6,4,3,1,2,1,0,1300,"Y;P;Q;Z"\n' "$MOD_HDR" > "$OUT/mz.bubbles.csv"
"$BIN" call -i "$OUT/mz.gfa" -c "$OUT/mz.bubbles.csv" -r ref -o "$OUT/mz" --cn --min-haplotypes 2 -q >/dev/null 2>&1
check "a module record dropped by --min-haplotypes leaves the bubble to UNIT (CN_METHOD:AC)" \
  "$(grep -v '^#' "$OUT/mz.region.vcf" | awk -F'\t' '$8~/SVTYPE=CNV/{match($8,/CN_METHOD=[A-Z_]+/); m=substr($8,RSTART+10,RLENGTH-10); match($8,/;AC=[0-9]+/); print m":"substr($8,RSTART+4,RLENGTH-4); exit}')" "UNIT:2"

# ---------------------------------------------------------------- a stretch too large to align
# Between the shared nodes S and K the reference walks N one-base nodes and hbig N others; hdel1 and
# hdel2 delete the whole stretch. At N = 2001 (4M cells, the old cap) the stretch is compared, and
# hbig's replacement of it is a substitution, not the deletion. At N = 8001 it is over the 64M cap and not
# compared: hbig found no deletion only because it was never compared, so its GT is unknown, not 0,
# it leaves AN, and the run summary names the bubble.
big_gfa() {   # big_gfa <out.gfa> <N>
  { printf 'H\tVN:Z:1.0\n'
    printf 'S\tF\t%s\nS\tS\t%s\nS\tK\t%s\nS\tG\t%s\n' "$(rseq 300 91)" "$(rseq 40 92)" "$(rseq 40 93)" "$(rseq 300 94)"
    # one base per node, drawn at random: a poly-A reference would be a tandem repeat, compared as
    # whole alleles rather than aligned
    awk -v n="$2" 'BEGIN{b="ACGT"; x=12345; for(i=1;i<=n;i++){ x=(x*16807)%2147483647; r=substr(b, int(x*4/2147483647)+1, 1)
         x=(x*16807)%2147483647; printf "S\tr%d\t%s\nS\ta%d\t%s\n", i, r, i, substr(b, int(x*4/2147483647)+1, 1)}}'
    printf 'L\tF\t+\tS\t+\t0M\nL\tS\t+\tK\t+\t0M\nL\tK\t+\tG\t+\t0M\n'
    awk -v n="$2" 'BEGIN{for(p=1;p<=2;p++){c=p==1?"r":"a"; printf "L\tS\t+\t%s1\t+\t0M\n", c
         for(i=1;i<n;i++) printf "L\t%s%d\t+\t%s%d\t+\t0M\n", c, i, c, i+1
         printf "L\t%s%d\t+\tK\t+\t0M\n", c, n}}'
    walk() { awk -v c="$1" -v n="$2" 'BEGIN{printf "F+,S+"; for(i=1;i<=n;i++) printf ",%s%d+", c, i; printf ",K+,G+"}'; }
    printf 'P\tref#0#chr1:1-%d\t%s\t*\n' $(($2 + 680)) "$(walk r "$2")"
    printf 'P\thbig#1#chr1\t%s\t*\n' "$(walk a "$2")"
    printf 'P\thdel1#1#chr1\tF+,S+,K+,G+\t*\nP\thdel2#1#chr1\tF+,S+,K+,G+\t*\n'; } > "$1"
  printf '%s\n1,S,+,K,+,%d,%d,4,3,1,2,1,0,%d,"r1;a1"\n' "$MOD_HDR" $((2 * $2)) $((2 * $2 + 2)) "$2" \
    > "${1%.gfa}.bubbles.csv"
}
big_gfa "$OUT/big4.gfa" 2001
"$BIN" call -i "$OUT/big4.gfa" -c "$OUT/big4.bubbles.csv" -r ref -o "$OUT/big4" -q >/dev/null 2>&1
VCF="$OUT/big4.region.vcf"
check "a 2001 x 2001-node stretch is compared: hbig carries its SUB, not the DEL (GT SUB:GT DEL:AN:AC)" \
  "$(gt_of "$VCF" hbig SUB):$(gt_of "$VCF" hbig DEL):$(info DEL AN):$(info DEL AC)" "1:0:3:2"
big_gfa "$OUT/big.gfa" 8001
big_log=$("$BIN" call -i "$OUT/big.gfa" -c "$OUT/big.bubbles.csv" -r ref -o "$OUT/big" 2>&1 >/dev/null)
VCF="$OUT/big.region.vcf"
check "a haplotype with an unaligned stretch is unknown on the bubble's DEL (GT:AN:AC)" \
  "$(gt_of "$VCF" hbig DEL):$(info DEL AN):$(info DEL AC)" ".:2:2"
printf '%s\n' "$big_log" | grep -q 'blocks too large to align: [0-9]* (sites 1); GT . there: 1' \
  && ok "the run summary names the bubble with the unaligned stretch" \
  || bad "the run summary does not name the unaligned bubble: $(printf '%s' "$big_log" | grep -i align)"
VCF="$OUT/cc.region.vcf"

# ---------------------------------------------------------------- the representative of a merged record
# Reference F S A B K G with A 80 bp and B 20 bp. hd1-hd3 delete A (80 bp), hd4 deletes A and B (100 bp);
# node Jaccard 0.8 merges them into one record. It reports the event most carriers hold, -80, not the
# largest member's -100, and SVLEN_RANGE still shows both.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%s\nS\tA\t%s\nS\tB\t%s\nS\tK\t%s\nS\tG\t%s\n' "$(rseq 300 101)" "$(rseq 40 102)" \
    "$(rseq 80 103)" "$(rseq 20 104)" "$(rseq 40 105)" "$(rseq 300 106)"
  for e in "F S" "S A" "A B" "B K" "S B" "S K" "K G"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref#0#chr1:1-780\tF+,S+,A+,B+,K+,G+\t*\n'
  for h in hd1 hd2 hd3; do printf 'P\t%s#1#chr1\tF+,S+,B+,K+,G+\t*\n' $h; done
  printf 'P\thd4#1#chr1\tF+,S+,K+,G+\t*\n'; } > "$OUT/rep.gfa"
printf '%s\n1,S,+,K,+,2,4,5,3,1,3,1,0,100,"A;B"\n' "$MOD_HDR" > "$OUT/rep.bubbles.csv"
"$BIN" call -i "$OUT/rep.gfa" -c "$OUT/rep.bubbles.csv" -r ref -o "$OUT/rep" -q >/dev/null 2>&1
VCF="$OUT/rep.region.vcf"
check "a merged DEL reports its most-carried member (SVLEN:SVLEN_RANGE:AC)" \
  "$(info DEL SVLEN):$(info DEL SVLEN_RANGE):$(info DEL AC)" "-80:-100,-80:4"
grep -q 'the exact event carried by the most haplotypes' "$VCF" \
  && ok "the header states the representative rule" || bad "the header does not state the representative rule"
VCF="$OUT/cc.region.vcf"

# ---------------------------------------------------------------- left-alignment
# S ends in G then ACACAC; hrep inserts (AC)x30 after S, where the graph puts it (POS 340). The same
# haplotype results with the insertion anywhere along the ACACAC run, so it is left-aligned to POS 334
# (before the run) with HOMLEN 6; K starts with G, so it cannot slide past S.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%sGACACAC\nS\tI\t%s\nS\tK\tG%s\nS\tG\t%s\n' "$(rseq 300 111)" "$(rseq 33 112)" \
    "$(printf 'AC%.0s' $(seq 30))" "$(rseq 39 113)" "$(rseq 300 114)"
  for e in "F S" "S I" "I K" "S K" "K G"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref#0#chr1:1-680\tF+,S+,K+,G+\t*\n'
  printf 'P\threp#1#chr1\tF+,S+,I+,K+,G+\t*\nP\thnon#1#chr1\tF+,S+,K+,G+\t*\n'; } > "$OUT/la.gfa"
printf '%s\n1,S,+,K,+,1,3,3,2,2,1,1,0,60,"I"\n' "$MOD_HDR" > "$OUT/la.bubbles.csv"
"$BIN" call -i "$OUT/la.gfa" -c "$OUT/la.bubbles.csv" -r ref -o "$OUT/la" -q >/dev/null 2>&1
VCF="$OUT/la.region.vcf"
check "an INS inside a repeat is left-aligned (POS:HOMLEN:CIPOS)" \
  "$(field INS 2):$(info INS HOMLEN):$(info INS CIPOS)" "334:6:0,6"
VCF="$OUT/cc.region.vcf"

# ---------------------------------------------------------------- --diploid columns
# The representative fixture's walks as two samples: s1#1 and s1#2 delete A (and B), s2#1 is
# reference-like, s2#2 deletes A. --diploid writes one column per sample with a phased GT and leaves
# the reference path out of the columns and of AN.
sed -e 's/^P\thd1#1#chr1/P\ts1#1#chr1/' -e 's/^P\thd4#1#chr1/P\ts1#2#chr1/' -e 's/^P\thd2#1#chr1/P\ts2#2#chr1/' \
    -e 's/^P\thd3#1#chr1\t.*/P\ts2#1#chr1\tF+,S+,A+,B+,K+,G+\t*/' "$OUT/rep.gfa" > "$OUT/dip.gfa"
"$BIN" call -i "$OUT/dip.gfa" -c "$OUT/rep.bubbles.csv" -r ref -o "$OUT/dip" --diploid -q >/dev/null 2>&1
VCF="$OUT/dip.region.vcf"
check "--diploid: columns, phased GTs and AN" \
  "$(awk -F'\t' '/^#CHROM/{print $10","$11} !/^#/{split($10,a,":"); split($11,b,":"); print a[1]"/"b[1]; exit}' "$VCF" | paste -sd: -):$(info DEL AN)" "s1,s2:1|1/0|1:4"
check "--diploid: HCN only on a copy-number record (DEL FORMAT)" "$(field DEL 9)" "GT:CN"
VCF="$OUT/cc.region.vcf"

# ---------------------------------------------------------------- substitutions
# Reference F S A K G with A 200 bp. hb replaces A by B (40 bp): a net loss of 160 bp, written as a DEL
# of -160 that replaces POS+1..END (200 bp) by B and is flagged SUBSTITUTION, its two sides spelled in REF
# and ALT. hc replaces
# A by C (210 bp): a net +10, under --min-sv-bp, so a SUB spelling both sequences. ht walks A2, a
# second node with A's sequence: no change, so no record.
SA=$(rseq 200 121); SB=$(rseq 40 122); SC=$(rseq 210 123)
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%s\nS\tA\t%s\nS\tA2\t%s\nS\tB\t%s\nS\tC\t%s\nS\tK\t%s\nS\tG\t%s\n' "$(rseq 300 124)" \
    "$(rseq 40 125)" "$SA" "$SA" "$SB" "$SC" "$(rseq 40 126)" "$(rseq 300 127)"
  for e in "F S" "S A" "A K" "S A2" "A2 K" "S B" "B K" "S C" "C K" "K G"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref#0#chr1:1-880\tF+,S+,A+,K+,G+\t*\n'
  printf 'P\thb#1#chr1\tF+,S+,B+,K+,G+\t*\nP\thc#1#chr1\tF+,S+,C+,K+,G+\t*\nP\tht#1#chr1\tF+,S+,A2+,K+,G+\t*\n'; } > "$OUT/sub.gfa"
printf '%s\n1,S,+,K,+,4,6,4,4,1,1,1,0,210,"A;A2;B;C"\n' "$MOD_HDR" > "$OUT/sub.bubbles.csv"
"$BIN" call -i "$OUT/sub.gfa" -c "$OUT/sub.bubbles.csv" -r ref -o "$OUT/sub" -q >/dev/null 2>&1
VCF="$OUT/sub.region.vcf"
check "a replacement losing 160 bp: one DEL spelled in REF/ALT (SVLEN:END-POS:REF is A:ALT is B:flag)" \
  "$(info DEL SVLEN):$(( $(info DEL END) - $(field DEL 2) )):$([ "$(field DEL 4)" = "$(field DEL 4 | cut -c1)$SA" ] && echo A):$([ "$(field DEL 5)" = "$(field DEL 4 | cut -c1)$SB" ] && echo B):$(grep -v '^#' "$VCF" | awk -F'\t' '$8~/SVTYPE=DEL;/' | grep -c ';SUBSTITUTION')" \
  "-160:200:A:B:1"
# REF and ALT are A and C less the bases they share at either end, after the base before.
trim_pair() { "${PYTHON:-python3}" -c 'import sys
a, b = sys.argv[1:3]; p = 0
while p < min(len(a), len(b)) and a[p] == b[p]: p += 1
s = 0
while s < min(len(a), len(b)) - p and a[len(a) - 1 - s] == b[len(b) - 1 - s]: s += 1
print(a[p:len(a) - s], b[p:len(b) - s])' "$1" "$2"; }
read -r TA TC <<< "$(trim_pair "$SA" "$SC")"
check "a replacement changing 10 bp: one SUB spelling both sequences (SVLEN:REF:ALT)" \
  "$(info SUB SVLEN):$([ "$(field SUB 4)" = "$(field SUB 4 | cut -c1)$TA" ] && echo A):$([ "$(field SUB 5)" = "$(field SUB 4 | cut -c1)$TC" ] && echo C)" "10:A:C"
check "a walk through another node with the reference's sequence gives no record (records:ht carries)" \
  "$(grep -vc '^#' "$VCF"):$(awk -F'\t' '/^#CHROM/{for(i=10;i<=NF;i++) if($i ~ "^ht#") c=i} !/^#/{split($c,g,":"); if(g[1]=="1") n++} END{print n+0}' "$VCF")" "2:0"
VCF="$OUT/cc.region.vcf"

# ---------------------------------------------------------------- copies of a repeat are counted when merging
# Reference F S K G. h2a/h2b insert the 60 bp cycle R1 R2 twice (120 bp), h5a/h5b five times (300 bp).
# Both walk the same nodes, so over node SETS they were one record (Jaccard 1) reporting one size for
# all four; counting each node as often as it is walked keeps 2 and 5 copies apart (Jaccard 0.4).
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%s\nS\tR1\t%s\nS\tR2\t%s\nS\tK\t%s\nS\tG\t%s\n' "$(rseq 300 131)" "$(rseq 40 132)" \
    "$(rseq 30 133)" "$(rseq 30 134)" "$(rseq 40 135)" "$(rseq 300 136)"
  for e in "F S" "S K" "S R1" "R1 R2" "R2 R1" "R2 K" "K G"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  r2='R1+,R2+,R1+,R2+'; r5='R1+,R2+,R1+,R2+,R1+,R2+,R1+,R2+,R1+,R2+'
  printf 'P\tref#0#chr1:1-680\tF+,S+,K+,G+\t*\n'
  printf 'P\th2a#1#chr1\tF+,S+,%s,K+,G+\t*\nP\th2b#1#chr1\tF+,S+,%s,K+,G+\t*\n' "$r2" "$r2"
  printf 'P\th5a#1#chr1\tF+,S+,%s,K+,G+\t*\nP\th5b#1#chr1\tF+,S+,%s,K+,G+\t*\n' "$r5" "$r5"; } > "$OUT/mset.gfa"
printf '%s\n1,S,+,K,+,2,4,5,3,1,2,2,0,300,"R1;R2"\n' "$MOD_HDR" > "$OUT/mset.bubbles.csv"
"$BIN" call -i "$OUT/mset.gfa" -c "$OUT/mset.bubbles.csv" -r ref -o "$OUT/mset" -q >/dev/null 2>&1
check "2 and 5 copies of one repeat are two INS records (SVLEN:AC per record)" \
  "$(grep -v '^#' "$OUT/mset.region.vcf" | awk -F'\t' '$8~/SVTYPE=INS;/{match($8,/SVLEN=[0-9]+/); s=substr($8,RSTART+6,RLENGTH-6); match($8,/;AC=[0-9]+/); print s"/"substr($8,RSTART+4,RLENGTH-4)}' | sort -n | paste -sd, -)" "120/2,300/2"

# ---------------------------------------------------------------- repeat expansions compared at base level
# The reference crosses the site F..G through an 80 bp (AC)n microsatellite. h60a/h60b carry 70 copies
# (+60 bp), h64 72 (+64), hsw the same length of (TG)n (every base differs), hm4 38 copies (-4 bp). Each
# allele is one node, so node by node they share nothing. Aligned at base level: +60 and +64 are 60 and
# 64 bp of AC inserted, which give the same haplotype to within 4 bp and so are one INS, left-aligned to
# the repeat's start (after F, POS 300); the (TG)n swap changes no length but 80 bp of sequence, a SUB
# spelling both; -4 bp is under --min-sv-bp and gives nothing.
ac() { local o="" i; for ((i = 0; i < $2; i++)); do o="$o$1"; done; printf '%s' "$o"; }
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tG\t%s\n' "$(rseq 300 141)" "$(rseq 300 142)"
  printf 'S\tR\t%s\nS\tX60\t%s\nS\tX64\t%s\nS\tXsw\t%s\nS\tXm4\t%s\n' "$(ac AC 40)" "$(ac AC 70)" "$(ac AC 72)" "$(ac TG 40)" "$(ac AC 38)"
  for n in R X60 X64 Xsw Xm4; do printf 'L\tF\t+\t%s\t+\t0M\nL\t%s\t+\tG\t+\t0M\n' "$n" "$n"; done
  printf 'P\tref#0#chr1:1-680\tF+,R+,G+\t*\n'
  printf 'P\th60a#1#chr1\tF+,X60+,G+\t*\nP\th60b#1#chr1\tF+,X60+,G+\t*\nP\th64#1#chr1\tF+,X64+,G+\t*\n'
  printf 'P\thsw#1#chr1\tF+,Xsw+,G+\t*\nP\thm4#1#chr1\tF+,Xm4+,G+\t*\n'; } > "$OUT/ra.gfa"
printf '%s\n1,F,+,G,+,5,7,6,5,1,2,1,76,144,"R;X60;X64;Xsw;Xm4"\n' "$MOD_HDR" > "$OUT/ra.bubbles.csv"
"$BIN" call -i "$OUT/ra.gfa" -c "$OUT/ra.bubbles.csv" -r ref -o "$OUT/ra" -q >/dev/null 2>&1
"$BIN" call -i "$OUT/ra.gfa" -c "$OUT/ra.bubbles.csv" -r ref -o "$OUT/ran" --no-align-alleles -q >/dev/null 2>&1
VCF="$OUT/ra.region.vcf"
check "a +60 and a +64 bp expansion are one INS, left-aligned (records:POS:SVLEN:AC:INSSEQ is AC)" \
  "$(grep -vc '^#' "$VCF"):$(field INS 2):$(info INS SVLEN):$(info INS AC):$([ "$(info INS INSSEQ)" = "$(ac AC 30)" ] && echo AC)" \
  "2:300:60:3:AC"
check "a same-length motif swap is a SUB spelling both alleles (SVLEN:REF:ALT)" \
  "$(info SUB SVLEN):$([ "$(field SUB 4)" = "$(field SUB 4 | cut -c1)$(ac AC 40)" ] && echo AC):$([ "$(field SUB 5)" = "$(field SUB 4 | cut -c1)$(ac TG 40)" ] && echo TG)" \
  "0:AC:TG"
# Node by node, each allele replaces the whole 80 bp node, but a substitution is sized by what changed, its
# edit distance: the 4 bp contraction is no 80 bp event, and node by node gives the same two records.
check "--no-align-alleles: the 4 bp contraction is no 80 bp SUB, the same two records (records)" \
  "$(grep -vc '^#' "$OUT/ran.region.vcf")" "2"

# ---------------------------------------------------------------- alleles aligned at base level
# hds's allele is the reference's with a 100 bp deletion and, 12 bp before it, one base changed, all on
# one node of its own. Node by node that is one replacement, written as a DEL flagged SUBSTITUTION that
# carries the base change; aligned, the deletion is the change and the base change a small variant.
SX=$(rseq 30 151); SY=$(rseq 12 152); SD=$(rseq 100 153); SZ=$(rseq 30 154)
snv=$([ "${SX: -1}" = "A" ] && echo C || echo A)
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tR\t%s\nS\tH\t%s\nS\tG\t%s\n' "$(rseq 300 155)" "$SX$SY$SD$SZ" "${SX%?}$snv$SY$SZ" "$(rseq 300 156)"
  printf 'L\tF\t+\tR\t+\t0M\nL\tR\t+\tG\t+\t0M\nL\tF\t+\tH\t+\t0M\nL\tH\t+\tG\t+\t0M\n'
  printf 'P\tref#0#chr1:1-772\tF+,R+,G+\t*\nP\thds#1#chr1\tF+,H+,G+\t*\n'; } > "$OUT/ds.gfa"
printf '%s\n1,F,+,G,+,2,4,2,2,1,1,1,72,172,"R;H"\n' "$MOD_HDR" > "$OUT/ds.bubbles.csv"
"$BIN" call -i "$OUT/ds.gfa" -c "$OUT/ds.bubbles.csv" -r ref -o "$OUT/ds" -q >/dev/null 2>&1
"$BIN" call -i "$OUT/ds.gfa" -c "$OUT/ds.bubbles.csv" -r ref -o "$OUT/dsn" --no-align-alleles -q >/dev/null 2>&1
check "a 100 bp deletion beside a base change: aligned, a plain DEL (records:SVLEN:SUBSTITUTION)" \
  "$(grep -vc '^#' "$OUT/ds.region.vcf"):$(grep -v '^#' "$OUT/ds.region.vcf" | grep -o 'SVLEN=-*[0-9]*' | head -1):$(grep -v '^#' "$OUT/ds.region.vcf" | grep -c SUBSTITUTION)" \
  "1:SVLEN=-100:0"
check "  ... node by node, one replacement flagged SUBSTITUTION" "$(grep -v '^#' "$OUT/dsn.region.vcf" | grep -c SUBSTITUTION)" "1"

# Four expansions of a 4 bp unit, +60, +68, +76 and +84 bp, one haplotype each. Single linkage chains
# them into one record (each within a fifth of the next); merged around a representative, every member
# must match it: the +84 representative takes +76 and +68, and +60 (under 0.8 of 84) is a record of its own.
au() { local o="" i; for ((i = 0; i < $1; i++)); do o="${o}AAAG"; done; printf '%s' "$o"; }
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tG\t%s\nS\tR\t%s\n' "$(rseq 300 161)" "$(rseq 300 162)" "$(au 10)"
  for k in 25 27 29 31; do printf 'S\tE%s\t%s\nL\tF\t+\tE%s\t+\t0M\nL\tE%s\t+\tG\t+\t0M\n' "$k" "$(au $k)" "$k" "$k"; done
  printf 'L\tF\t+\tR\t+\t0M\nL\tR\t+\tG\t+\t0M\nP\tref#0#chr1:1-640\tF+,R+,G+\t*\n'
  for k in 25 27 29 31; do printf 'P\th%s#1#chr1\tF+,E%s+,G+\t*\n' "$k" "$k"; done; } > "$OUT/lad.gfa"
printf '%s\n1,F,+,G,+,5,7,5,5,1,1,1,40,124,"R;E25;E27;E29;E31"\n' "$MOD_HDR" > "$OUT/lad.bubbles.csv"
"$BIN" call -i "$OUT/lad.gfa" -c "$OUT/lad.bubbles.csv" -r ref -o "$OUT/lad" -q >/dev/null 2>&1
check "repeat-length alleles merge around a representative, not by chaining (SVLEN/AC per INS record)" \
  "$(grep -v '^#' "$OUT/lad.region.vcf" | awk -F'\t' '$8~/SVTYPE=INS;/{match($8,/SVLEN=[0-9]+/); s=substr($8,RSTART+6,RLENGTH-6); match($8,/;AC=[0-9]+/); print s"/"substr($8,RSTART+4,RLENGTH-4)}' | sort -n | paste -sd, -)" \
  "60/1,84/3"

# ---------------------------------------------------------------- aligned alleles: fixes from the module check
# One generator writes several single-bubble cases. Each case's interior nodes, and the haplotypes as
# lists of interior sequences; a sequence used twice is one node.
"${PYTHON:-python3}" - "$OUT" <<'PYEOF2'
import random, sys
out = sys.argv[1]
def rnd(n, seed):
    r = random.Random(seed); return ''.join(r.choice('ACGT') for _ in range(n))
HDR = ('bubble_id,source,source_orient,sink,sink_orient,inside_node_count,total_node_count,path_support,'
       'distinct_alleles,ref_allele_support,alt_allele_support_max,alt_allele_support_min,min_inside_bp,'
       'max_inside_bp,inside_nodes')
def write(name, ref, haps, seed):
    F, G = rnd(300, seed), rnd(300, seed + 1)
    nodes, links, paths = {'F': F, 'G': G}, set(), {}
    ids = {}
    def nid(s):
        if s not in ids:
            ids[s] = 'n%d' % (len(ids) + 1); nodes[ids[s]] = s
        return ids[s]
    total = 300 + sum(map(len, ref)) + 300
    for name_p, segs in [('ref#0#chr1:1-%d' % total, ref)] + list(haps.items()):
        walk = ['F'] + [nid(x) for x in segs] + ['G']
        for a, b in zip(walk, walk[1:]): links.add((a, b))
        paths[name_p] = walk
    with open(f'{out}/{name}.gfa', 'w') as o:
        o.write('H\tVN:Z:1.0\n')
        for k, v in nodes.items(): o.write(f'S\t{k}\t{v}\n')
        for a, b in sorted(links): o.write(f'L\t{a}\t+\t{b}\t+\t0M\n')
        for k, v in paths.items(): o.write(f'P\t{k}\t' + ','.join(x + '+' for x in v) + '\t*\n')
    inside = sorted(set(n for w in paths.values() for n in w) - {'F', 'G'})
    lens = [sum(len(nodes[n]) for n in w[1:-1]) for w in paths.values()]
    with open(f'{out}/{name}.bubbles.csv', 'w') as o:
        o.write(HDR + '\n')
        o.write('1,F,+,G,+,%d,%d,%d,%d,1,1,1,%d,%d,"%s"\n' % (len(inside), len(inside) + 2, len(paths),
                len(paths), min(lens), max(lens), ';'.join(inside)))
A, B, C, D = rnd(400, 3), rnd(300, 4), rnd(400, 5), rnd(400, 6)
L, I = rnd(3950, 7), rnd(60, 8)
# a 3950 bp deletion carried by three haplotypes; h3 also inserts 60 bp, on a node of its own, so its
# core is too large for the banded alignment and goes to minimap2
write('mm', [A, L, B, C], {'h1#1#c': [A, B, C], 'h2#1#c': [A, B, C], 'h3#1#c': [A, B[:100] + I + B[100:], C]}, 11)
# one node holding a 500 bp insertion and, 1.5 kb further, a 500 bp deletion: no length change overall
B2, C2, X = rnd(1500, 12), rnd(1000, 13), rnd(500, 14)
write('band', [A, B2, C2, D], {'h1#1#c': [A + X + B2 + C2[500:] + D]}, 15)
# a 20 kb insertion with 100 bp shared either side, split over nodes in h1 and on one node in h2
L2 = rnd(20000, 16)
write('onesided', [A[:100], C[:100]], {'h1#1#c': [A[:100], L2, C[:100]], 'h2#1#c': [A[:100] + L2 + C[:100]]}, 17)
# deletions of 27.7 and 32.2 kb: their haplotypes differ by a 4.5 kb segment, so they are two alleles
D1, D2 = rnd(27700, 18), rnd(4500, 19)
write('bigdel', [A, D1, D2, B], {'h1#1#c': [A, D2, B], 'h2#1#c': [A, B], 'h3#1#c': [A, B]}, 20)
# --diploid: S1 has numbered haplotypes 1 and 2; S2 two unassigned contigs (#U#), one deleting B;
# S3 has haplotype 1 in two contigs that both cross the site, one of them deleting B
write('dip', [A, B, C], {'S1#1#a': [A, B, C], 'S1#2#b': [A, B, C], 'S2#U#c': [A, C], 'S2#U#d': [A, B, C],
                        'S3#1#e': [A, C], 'S3#1#f': [A, B, C], 'S3#2#g': [A, C]}, 21)
# the same 20 kb insertion in four haplotypes, each with 30 base changes of its own: two of them differ
# in 60 bases but not in length, so they are one allele
def snv(x, k, seed):
    r = random.Random(seed); x = list(x)
    for p in r.sample(range(len(x)), k): x[p] = 'ACGT'[('ACGT'.index(x[p]) + r.randint(1, 3)) % 4]
    return ''.join(x)
write('snvins', [A, C], {'h%d#1#c' % i: [A, snv(L2, 30, 30 + i), C] for i in range(1, 5)}, 22)
# --diploid with an allele-VCF record: S1 deletes B on haplotype 1 and inserts X on haplotype 2,
# S2 is reference-like on 1 and deletes B on 2, S3 has two contigs on haplotype 1 (unknown) and X on 2
X2 = rnd(200, 23)
# a 250 kb deletion, its end 10 bp apart in h1 and h2: one allele, past any cap on the window compared.
# h3 deletes what h1 does, on a flank split into two nodes of its own: the same event by other nodes
D3 = rnd(250000, 25)
write('hugedel', [A, D3[:10], D3[10:], B], {'h1#1#c': [A, B], 'h2#1#c': [A, D3[:10], B],
                                           'h3#1#c': [A[:200], A[200:], B]}, 26)
# two different 1 kb deletions 300 bp apart in unique sequence (h1 and h2 delete P1 P2, h3 P2 P3): the
# same size, and the haplotypes alike over most of the event, but where they differ, they share nothing
P1, P2, P3, P4 = rnd(300, 27), rnd(700, 28), rnd(300, 29), rnd(700, 30)
write('offdel', [A, P1, P2, P3, P4, B], {'h1#1#c': [A, P3, P4, B], 'h2#1#c': [A, P3, P4, B], 'h3#1#c': [A, P1, P4, B]}, 31)
# the 20 kb insertion again, each carrier's 30 base changes now within 200 bp of each other: the same
# event, its carriers differing in one patch (one record)
def snv_in(x, k, lo, hi, seed):
    r = random.Random(seed); x = list(x)
    for p in r.sample(range(lo, hi), k): x[p] = 'ACGT'[('ACGT'.index(x[p]) + r.randint(1, 3)) % 4]
    return ''.join(x)
write('snvpatch', [A, C], {'h%d#1#c' % i: [A, snv_in(L2, 30, 9000, 9200, 40 + i), C] for i in range(1, 5)}, 32)
# one 2 kb sequence inserted before a 200 bp stretch of unique sequence (h1, h2) and after it (h3): two sites
I2, Q = rnd(2000, 33), rnd(200, 34)
write('twoins', [A, Q, B], {'h1#1#c': [A, I2, Q, B], 'h2#1#c': [A, I2, Q, B], 'h3#1#c': [A, Q, I2, B]}, 35)
PYEOF2
cnt() { grep -v '^#' "$1" | grep -c "$2"; }
"$BIN" call -i "$OUT/mm.gfa" -c "$OUT/mm.bubbles.csv" -r ref -o "$OUT/mm" -q >/dev/null 2>&1
check "a deletion found by minimap2 and by the banded alignment is one record (DEL records:AC)" \
  "$(cnt "$OUT/mm.region.vcf" 'SVTYPE=DEL;'):$(grep -v '^#' "$OUT/mm.region.vcf" | grep 'SVTYPE=DEL;' | grep -o ';AC=[0-9]*' | cut -c5-)" "1:3"
"$BIN" call -i "$OUT/band.gfa" -c "$OUT/band.bubbles.csv" -r ref -o "$OUT/band" -q >/dev/null 2>&1
check "an insertion and a deletion 1.5 kb apart on one node are found as such (SVLEN of INS,DEL)" \
  "$(grep -v '^#' "$OUT/band.region.vcf" | grep -o 'SVTYPE=[A-Z]*;SVLEN=-*[0-9]*' | sed 's/SVTYPE=//;s/;SVLEN=/:/' | sort | paste -sd, -)" "DEL:-500,INS:500"
"$BIN" call -i "$OUT/onesided.gfa" -c "$OUT/onesided.bubbles.csv" -r ref -o "$OUT/onesided" -q >/dev/null 2>&1
VCF="$OUT/onesided.region.vcf"
check "a 20 kb insertion with short shared flanks is one INS record (records:SVLEN:AC)" \
  "$(grep -vc '^#' "$OUT/onesided.region.vcf"):$(info INS SVLEN):$(grep -v '^#' "$OUT/onesided.region.vcf" | grep -o ';AC=[0-9]*' | cut -c5-)" "1:20000:2"
"$BIN" call -i "$OUT/bigdel.gfa" -c "$OUT/bigdel.bubbles.csv" -r ref -o "$OUT/bigdel" -q >/dev/null 2>&1
check "deletions of 27.7 and 32.2 kb stay two records (SVLEN/AC)" \
  "$(grep -v '^#' "$OUT/bigdel.region.vcf" | awk -F'\t' '{match($8,/SVLEN=-?[0-9]+/); s=substr($8,RSTART+6,RLENGTH-6); match($8,/;AC=[0-9]+/); print s"/"substr($8,RSTART+4,RLENGTH-4)}' | sort -n | paste -sd, -)" "-32200/2,-27700/1"
"$BIN" call -i "$OUT/dip.gfa" -c "$OUT/dip.bubbles.csv" -r ref -o "$OUT/dip" --diploid -q >/dev/null 2>&1
VCF="$OUT/dip.region.vcf"
dipgt() { awk -F'\t' -v s="$1" '/^#CHROM/{for(i=10;i<=NF;i++) if($i==s) c=i} !/^#/{split($c,g,":"); print g[1]; exit}' "$OUT/dip.region.vcf"; }
check "--diploid: unassigned contigs are two unphased haplotypes; of two contigs crossing, the carrying one makes the slot a carrier (S1:S2:S3)" \
  "$(dipgt S1):$(dipgt S2):$(dipgt S3)" "0|0:1/0:1|1"
check "--diploid: AN and AC count the alleles written" "$(info DEL AN):$(info DEL AC)" "6:3"
"$BIN" call -i "$OUT/hugedel.gfa" -c "$OUT/hugedel.bubbles.csv" -r ref -o "$OUT/hugedel" -q >/dev/null 2>&1
check "a 250 kb deletion whose end moves by 10 bp is one record, with an identical one on other nodes (records:AC)" \
  "$(grep -vc '^#' "$OUT/hugedel.region.vcf"):$(grep -v '^#' "$OUT/hugedel.region.vcf" | grep -o ';AC=[0-9]*' | cut -c5-)" "1:3"
"$BIN" call -i "$OUT/offdel.gfa" -c "$OUT/offdel.bubbles.csv" -r ref -o "$OUT/offdel" -q >/dev/null 2>&1
check "two same-size deletions 300 bp apart in unique sequence are two records (POS/SVLEN/AC)" \
  "$(grep -v '^#' "$OUT/offdel.region.vcf" | awk -F'\t' '{match($8,/SVLEN=-?[0-9]+/); s=substr($8,RSTART+6,RLENGTH-6); match($8,/;AC=[0-9]+/); print $2"/"s"/"substr($8,RSTART+4,RLENGTH-4)}' | paste -sd, -)" \
  "700/-1000/2,1000/-1000/1"
"$BIN" call -i "$OUT/snvpatch.gfa" -c "$OUT/snvpatch.bubbles.csv" -r ref -o "$OUT/snvpatch" -q >/dev/null 2>&1
check "one 20 kb insertion whose carriers differ in a 200 bp patch of base changes is one record (records:AC)" \
  "$(grep -vc '^#' "$OUT/snvpatch.region.vcf"):$(grep -v '^#' "$OUT/snvpatch.region.vcf" | grep -o ';AC=[0-9]*' | cut -c5-)" "1:4"
"$BIN" call -i "$OUT/twoins.gfa" -c "$OUT/twoins.bubbles.csv" -r ref -o "$OUT/twoins" -q >/dev/null 2>&1
check "one sequence inserted at two sites 200 bp apart in unique sequence is two records (POS/SVLEN/AC)" \
  "$(grep -v '^#' "$OUT/twoins.region.vcf" | awk -F'\t' '{match($8,/SVLEN=-?[0-9]+/); s=substr($8,RSTART+6,RLENGTH-6); match($8,/;AC=[0-9]+/); print $2"/"s"/"substr($8,RSTART+4,RLENGTH-4)}' | paste -sd, -)" \
  "700/2000/2,900/2000/1"
"$BIN" call -i "$OUT/snvins.gfa" -c "$OUT/snvins.bubbles.csv" -r ref -o "$OUT/snvins" -q >/dev/null 2>&1
check "one 20 kb insertion whose carriers differ in 30 base changes each is one record (records:AC)" \
  "$(grep -vc '^#' "$OUT/snvins.region.vcf"):$(grep -v '^#' "$OUT/snvins.region.vcf" | grep -o ';AC=[0-9]*' | cut -c5-)" "1:4"

# The UNIT loop counts copies only when its record survives: with --min-haplotypes 2 the UNIT record
# (one carrier) is dropped, the loop stays in the walks, and the 70 bp deletion keeps both carriers.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%s\nS\tA\t%s\nS\tD\t%s\nS\tL\t%s\nS\tC\t%s\nS\tK\t%s\nS\tG\t%s\n' "$(rseq 300 171)" "$(rseq 40 172)" \
    "$(rseq 40 173)" "$(rseq 30 174)" "$(rseq 20 175)" "$(rseq 40 176)" "$(rseq 40 177)" "$(rseq 300 178)"
  for e in "F S" "S A" "A D" "D L" "L L" "L C" "A C" "C K" "K G"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref#0#chr1:1-830\tF+,S+,A+,D+,L+,L+,C+,K+,G+\t*\n'
  printf 'P\th1#1#chr1\tF+,S+,A+,D+,L+,L+,L+,L+,L+,C+,K+,G+\t*\nP\th2#1#chr1\tF+,S+,A+,D+,L+,L+,L+,C+,K+,G+\t*\n'
  printf 'P\th3#1#chr1\tF+,S+,A+,D+,L+,L+,L+,C+,K+,G+\t*\nP\th4#1#chr1\tF+,S+,A+,C+,K+,G+\t*\n'
  printf 'P\th5#1#chr1\tF+,S+,A+,C+,K+,G+\t*\nP\th6#1#chr1\tF+,S+,A+,D+,L+,L+,C+,K+,G+\t*\n'; } > "$OUT/rps.gfa"
printf '%s\n1,S,+,K,+,4,6,7,4,2,2,1,40,170,"A;D;L;C"\n' "$MOD_HDR" > "$OUT/rps.bubbles.csv"
"$BIN" call -i "$OUT/rps.gfa" -c "$OUT/rps.bubbles.csv" -r ref -o "$OUT/rps" --min-haplotypes 2 -q >/dev/null 2>&1
VCF="$OUT/rps.region.vcf"
check "with --min-haplotypes 2 a dropped UNIT record leaves its loop in the walks (DEL SVLEN:AC)" \
  "$(info DEL SVLEN):$(info DEL AC)" "-70:2"

# A locus with no sites: bubble and panphorte write a header-only CSV, and call a header-only VCF.
printf 'H\tVN:Z:1.0\nS\t1\tACGTACGTACGT\nP\tref\t1+\t*\nP\th1\t1+\t*\n' > "$OUT/nosite.gfa"
"$BIN" bubble -i "$OUT/nosite.gfa" -r ref -o "$OUT/nosite_b" -q >/dev/null 2>&1
"$BIN" call -i "$OUT/nosite_b.sorted.gfa" -c "$OUT/nosite_b.bubbles.csv" -r ref -o "$OUT/nosite" -q >/dev/null 2>&1
check "a locus with no sites gives a header-only VCF and exit 0 (exit:records:header)" \
  "$?:$(grep -vc '^#' "$OUT/nosite.region.vcf" 2>/dev/null):$(grep -c '^#CHROM' "$OUT/nosite.region.vcf" 2>/dev/null)" "0:0:1"

# An assembly gap on the reference's own allele: its sequence there is unknown, so the site is skipped
# and counted, rather than the haplotypes skipping the gap being called a deletion of N bases.
NS=$(printf 'N%.0s' $(seq 200))
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%s\nS\tN\t%s\nS\tK\t%s\nS\tG\t%s\n' "$(rseq 300 181)" "$(rseq 40 182)" "$NS" "$(rseq 40 183)" "$(rseq 300 184)"
  for e in "F S" "S N" "N K" "S K" "K G"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref#0#chr1:1-880\tF+,S+,N+,K+,G+\t*\nP\th1#1#chr1\tF+,S+,N+,K+,G+\t*\n'
  printf 'P\th2#1#chr1\tF+,S+,K+,G+\t*\nP\th3#1#chr1\tF+,S+,K+,G+\t*\n'; } > "$OUT/refgap.gfa"
printf '%s\n1,S,+,K,+,1,3,4,2,2,2,2,0,200,"N"\n' "$MOD_HDR" > "$OUT/refgap.bubbles.csv"
refgap_log=$("$BIN" call -i "$OUT/refgap.gfa" -c "$OUT/refgap.bubbles.csv" -r ref -o "$OUT/refgap" 2>&1 >/dev/null)
check "a site whose reference allele is an assembly gap is skipped and counted (records:logged)" \
  "$(grep -vc '^#' "$OUT/refgap.region.vcf"):$(printf '%s\n' "$refgap_log" | grep -c 'sites skipped, assembly gap in the reference: 1 (1)')" "0:1"

# An assembly gap in the site's source and sink: every path crossing the site shares them, so they say
# nothing about the allele, and the 120 bp deletion inside is called with its carriers.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%s%s%s\nS\tA\t%s\nS\tK\t%s%s%s\nS\tG\t%s\n' "$(rseq 300 191)" "$(rseq 40 192)" "$NS" \
    "$(rseq 20 193)" "$(rseq 120 194)" "$(rseq 20 195)" "$NS" "$(rseq 40 196)" "$(rseq 300 197)"
  for e in "F S" "S A" "A K" "S K" "K G"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref#0#chr1:1-1240\tF+,S+,A+,K+,G+\t*\nP\th1#1#chr1\tF+,S+,A+,K+,G+\t*\n'
  printf 'P\th2#1#chr1\tF+,S+,K+,G+\t*\nP\th3#1#chr1\tF+,S+,K+,G+\t*\n'; } > "$OUT/bgap.gfa"
printf '%s\n1,S,+,K,+,1,3,4,2,2,2,2,0,120,"A"\n' "$MOD_HDR" > "$OUT/bgap.bubbles.csv"
bgap_log=$("$BIN" call -i "$OUT/bgap.gfa" -c "$OUT/bgap.bubbles.csv" -r ref -o "$OUT/bgap" 2>&1 >/dev/null)
VCF="$OUT/bgap.region.vcf"
check "an assembly gap in a site's source and sink leaves its deletion called (SVLEN:AN:AC:set missing)" \
  "$(info DEL SVLEN):$(info DEL AN):$(info DEL AC):$(printf '%s\n' "$bgap_log" | grep -c 'GT . across an assembly gap')" "-120:3:2:0"

# ---- a folded repeat as a unit graph (panphorte's default approximate fold) ----------------------
# Copies of a 64 bp unit, two with a SNV of their own: the unit graph keeps every copy's sequence, so no node
# loops and call used to see the copies as INS/DEL. Read from panphorte's unit_graph.tsv, the unit graph is
# one UNIT unit counted at its anchor: ref 2, h1 2, h2 3, h0 0 copies.
UG_UNIT="TATACGCCATTCAATAACAAGCCAATTCGCTGGGTCAACTCCCAGCCAGACGCCCCGTTGGCGC"
UG_A="$(blk 2000 101)"; UG_B="$(blk 2000 202)"
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t5\t%s\nS\t4\t%s\n' "$UG_A" "$UG_UNIT" "${UG_UNIT:0:30}A${UG_UNIT:31}" \
    "${UG_UNIT:0:40}$([ "${UG_UNIT:40:1}" = C ] && echo G || echo C)${UG_UNIT:41}" "$UG_B"
  for e in "1 2" "2 2" "2 3" "3 2" "3 4" "2 4" "1 4" "2 5" "5 2"; do set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tref\t1+,2+,3+,4+\t*\nP\th1\t1+,2+,2+,4+\t*\nP\th2\t1+,2+,5+,2+,4+\t*\nP\th0\t1+,4+\t*\n'; } > "$OUT/ug.gfa"
"$BIN" bubble -i "$OUT/ug.gfa" -r ref -o "$OUT/ug_b" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/ug_b.sorted.gfa" -b "$OUT/ug_b" -o "$OUT/ug_p" -r ref --min-similarity 0.95 -q >/dev/null 2>&1
"$BIN" call -i "$OUT/ug_p.normalized.sorted.gfa" -b "$OUT/ug_p" -r ref -o "$OUT/ug_c" --no-per-bubble-vcf -q >/dev/null 2>&1
ug_cn() { awk -F'\t' '/^#CHROM/{for(i=10;i<=NF;i++) s[i]=$i; next} !/^#/ && $8 ~ /CN_METHOD=UNIT/ {
            n=split($9,f,":"); for(k=1;k<=n;k++) if(f[k]=="CN") c=k
            for(i=10;i<=NF;i++){ split($i,v,":"); printf "%s %s ", s[i], v[c] } }' "$1"; }
check "unit graph: copies counted at the anchor as one UNIT unit (ref 2, h1 2, h2 3, h0 0), no copy INS/DEL" \
  "$(ug_cn "$OUT/ug_c.region.vcf")|$(grep -v '^#' "$OUT/ug_c.region.vcf" | grep -cE 'SVTYPE=(INS|DEL)')" "ref 2 h1 2 h2 3 h0 0 |0"
# The record spans the reference's two copies from the first unit graph node, not the segment the unit is keyed
# by (on LPA, 29 bp of a 5547 bp copy).
check "unit graph: the record spans the reference's copies (END-POS = REF_CN x RU_LEN = 128)" \
  "$(grep -v '^#' "$OUT/ug_c.region.vcf" | awk -F'\t' '/CN_METHOD=UNIT/{match($8,/(^|;)END=[0-9]+/); e=$8; e=substr(e,RSTART,RLENGTH); sub(/.*END=/,"",e); match($8,/RU_LEN=[0-9]+/); print e-$2":"substr($8,RSTART+7,RLENGTH-7)}')" "128:64"

# Counted copies leave the walks before they are compared, so what lay on either side of them comes together.
# Reference F S U J U J U T K G: three copies of the 300 bp unit U with a 1 bp junction J between them, then
# T (60 bp). hX has one copy and hY three, both with V (T less 4 bp) for T. hX's walk F S U V K G met the
# reference's F S J J T K G, and J J T against V was one substitution at the first junction spelling all three
# places (ANKRD36C altref, 11 kb). The block is split where copies left: T against V after the last copy, one
# record carried by hX and hY, its REF the reference's bases at POS; each junction alone, too small to report.
# V replaces 60 bp in the middle of T (120 bp) with 56 others, a change the alignment sees too; the copies
# that left the walks lie before it, so it is aligned against the piece of reference after them.
JU=$(blk 300 211); JT=$(blk 120 212); JT2="${JT:0:30}$(blk 56 217)${JT:90}"
J_R='F+,S+,U+,J+,U+,J+,U+,T+,K+,G+'; J_X='F+,S+,U+,V+,K+,G+'; J_Y='F+,S+,U+,J+,U+,J+,U+,V+,K+,G+'; J_1='F+,S+,U+,J+,U+,T+,K+,G+'
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%s\nS\tU\t%s\nS\tJ\tA\nS\tT\t%s\nS\tV\t%s\nS\tK\t%s\nS\tG\t%s\n' "$(blk 500 213)" "$(blk 40 214)" \
    "$JU" "$JT" "$JT2" "$(blk 40 215)" "$(blk 500 216)"
  for p in "$J_R" "$J_X" "$J_Y" "$J_1"; do printf '%s\n' "$p" | tr ',' '\n' | awk 'NR > 1 { print prev "\t" $0 } { prev = $0 }'; done \
    | sort -u | awk -F'\t' '{ printf "L\t%s\t%s\t%s\t%s\t0M\n", substr($1, 1, length($1) - 1), substr($1, length($1)),
                               substr($2, 1, length($2) - 1), substr($2, length($2)) }'
  printf 'P\tref#0#chr1:1-2102\t%s\t*\nP\thX#1#chr1\t%s\t*\nP\thY#1#chr1\t%s\t*\nP\th1#1#chr1\t%s\t*\n' "$J_R" "$J_X" "$J_Y" "$J_1"; } > "$OUT/junc.gfa"
printf '%s\n1,S,+,K,+,4,6,4,4,1,1,1,356,962,"U;J;T;V"\n' "$MOD_HDR" > "$OUT/junc.bubbles.csv"
printf 'input_bubble_id\tunit_bp\tcanonical_motif\tcopy_quantum\tunit\tanchor_pos\tanchor_nodes\tn_nodes\tnodes\n1\t300\t%s\t1\t%s\t1\tU\t1\tU\n' \
  "$JU" "$JU" > "$OUT/junc.panphorte.unit_graph.tsv"
"$BIN" call -i "$OUT/junc.gfa" -c "$OUT/junc.bubbles.csv" --panphorte-prefix "$OUT/junc" -r ref -o "$OUT/junc" \
  --no-per-bubble-vcf -q >/dev/null 2>&1
"$BIN" call -i "$OUT/junc.gfa" -c "$OUT/junc.bubbles.csv" --panphorte-prefix "$OUT/junc" -r ref -o "$OUT/juncn" \
  --no-per-bubble-vcf --no-align-alleles -q >/dev/null 2>&1
J_REF=$(printf '%s' "$(blk 500 213)$(blk 40 214)${JU}A${JU}A${JU}")
j_sub() { grep -v '^#' "$1" | awk -F'\t' -v r="$J_REF$JT" '$8 ~ /SVTYPE=SUB/ { match($8, /;AC=[0-9]+/)
      printf "%s %s %s ", ($2 >= 1442 && $2 < 1442 + 120 ? "after" : $2), substr($8, RSTART + 4, RLENGTH - 4),
             (substr(r, $2, length($4)) == $4 ? "yes" : "no") }'; }
check "copies removed between pieces: one SUB in T after the last copy, REF at POS (records POS AC REF-matches)" \
  "$(j_sub "$OUT/junc.region.vcf")" "after 2 yes "
check "copies removed between pieces, node by node: the same (records POS AC REF-matches)" \
  "$(j_sub "$OUT/juncn.region.vcf")" "after 2 yes "

# Two phase-rotated units of one motif at one site are two one-node unit graphs (an unsplit node cannot hold
# both rotations); hX walks one copy of each. Counted apart, hX read as a loss on the first and a gain on the
# second; the table's site and motif make them one repeat: ref 2, h1 2, h2 3, hX 2.
PCC="$(cd "$(dirname "$0")/.." && pwd)/synthetic_data/panphorte_cases"
"$BIN" bubble -i "$PCC/refine_rep_no_loop.gfa" -r ref -o "$OUT/rot_b" --min-variant-bp 0 --merge-nearby-bp 0 -q >/dev/null 2>&1
"$BIN" panphorte -i "$OUT/rot_b.sorted.gfa" -b "$OUT/rot_b" -o "$OUT/rot_p" -r ref --threads 1 -q >/dev/null 2>&1
"$BIN" call -i "$OUT/rot_p.normalized.sorted.gfa" -b "$OUT/rot_p" -r ref -o "$OUT/rot_c" --no-per-bubble-vcf -q >/dev/null 2>&1
check "two rotations of one motif at one site are one unit (CN per haplotype | CNV records)" \
  "$(ug_cn "$OUT/rot_c.region.vcf")|$(grep -v '^#' "$OUT/rot_c.region.vcf" | grep -c 'SVTYPE=CNV')" "ref 2 h1 2 h2 3 hX 2 |1"

# A site every path crosses twice, the second time inverted: the scored walk holds the first crossing only.
# h2's 100 bp deletion is in the second, so no record can describe it; it is counted and flagged
# (MULTI_CROSSING_UNSCORED) rather than lost silently. h1's, in the first crossing, is called.
BC="$(cd "$(dirname "$0")/.." && pwd)/synthetic_data/bubble_cases"
"$BIN" call -i "$BC/inverted_second_crossing.gfa" -c "$BC/inverted_second_crossing.csv" -r ref -o "$OUT/isc" \
  --no-per-bubble-vcf > /dev/null 2> "$OUT/isc.err"
check "twice-crossed site: an SV-sized change outside the scored walk is flagged (UNSCORED | warning | DEL AC)" \
  "$(grep -o 'MULTI_CROSSING_UNSCORED=[0-9]*' "$OUT/isc.region.vcf" | sort -u | tr '\n' ' ')|$(grep -c 'MULTI_CROSSING_UNSCORED' "$OUT/isc.err")|$(grep -v '^#' "$OUT/isc.region.vcf" | grep -o 'SVTYPE=DEL;SVLEN=-100' | wc -l | tr -d ' ')" \
  "MULTI_CROSSING_UNSCORED=1 |1|1"

# hapA duplicates X Y, hapB Y Z, one carrier each: two different tandem duplications, two INS records,
# the same whichever path comes first. (The PEAK route labelled one record by a gaining allele's busiest
# node, and took that allele from the P-line order.)
PS=$(blk 100 41); PX=$(blk 60 43); PY=$(blk 60 47); PZ=$(blk 60 53); PT=$(blk 100 59)
pk_gfa() {  # pk_gfa <out> <order: ab|ba>
  { printf 'H\tVN:Z:1.0\nS\ts\t%s\nS\tX\t%s\nS\tY\t%s\nS\tZ\t%s\nS\tt\t%s\n' "$PS" "$PX" "$PY" "$PZ" "$PT"
    printf 'L\ts\t+\tX\t+\t0M\nL\tX\t+\tY\t+\t0M\nL\tY\t+\tZ\t+\t0M\nL\tZ\t+\tt\t+\t0M\nL\tY\t+\tX\t+\t0M\nL\tZ\t+\tY\t+\t0M\n'
    printf 'P\tref#0#chr1:1-380\ts+,X+,Y+,Z+,t+\t*\n'
    a='P\thapA#1#c1\ts+,X+,Y+,X+,Y+,Z+,t+\t*\n'; b='P\thapB#1#c1\ts+,X+,Y+,Z+,Y+,Z+,t+\t*\n'
    if [ "$2" = ab ]; then printf "$a$b"; else printf "$b$a"; fi; } > "$1"; }
pk_gfa "$OUT/pk_ab.gfa" ab; pk_gfa "$OUT/pk_ba.gfa" ba
printf 'bubble_id,source,source_orient,sink,sink_orient,inside_node_count,total_node_count,path_support,distinct_alleles,ref_allele_support,alt_allele_support_max,alt_allele_support_min,min_inside_bp,max_inside_bp,inside_nodes\n1,s,+,t,+,3,5,3,3,1,1,1,180,300,"X;Y;Z"\n' > "$OUT/pk.csv"
for o in ab ba; do "$BIN" call -i "$OUT/pk_$o.gfa" -c "$OUT/pk.csv" -r 'ref#0#chr1:1-380' -o "$OUT/pk_$o" --cn --no-per-bubble-vcf -q >/dev/null 2>&1; done
pk_lab() { grep -v '^#' "$1" | awk -F'\t' '/^#CHROM/{next} {print $3}' | paste -sd, -; }
check "two tandem duplications are two INS records whatever the P-line order (IDs, paths A,B | B,A)" \
  "$(pk_lab "$OUT/pk_ab.region.vcf")|$(pk_lab "$OUT/pk_ba.region.vcf")" "bubble1_INS_X,bubble1_INS_Y|bubble1_INS_X,bubble1_INS_Y"

# --diploid with a haplotype in two contigs that both cross the site. Agreeing (both delete d), they are the
# haplotype there: S1 carries the DEL (AC 1); read as unknown, the record was written with AC=0. Disagreeing
# (one deletes d), the contig carrying it makes S1 a carrier, as AC counts it per path.
dip_gfa() {  # dip_gfa <out> <S1 contig 2's walk>
  { printf 'H\tVN:Z:1.0\nS\ta\t%s\nS\td\t%s\nS\tb\t%s\n' "$(blk 200 61)" "$(blk 120 62)" "$(blk 200 63)"
    printf 'L\ta\t+\td\t+\t0M\nL\td\t+\tb\t+\t0M\nL\ta\t+\tb\t+\t0M\n'
    printf 'P\tref#0#chr1:1-520\ta+,d+,b+\t*\nP\tS1#1#c1\ta+,b+\t*\nP\tS1#1#c2\t%s\t*\n' "$2"
    printf 'P\tS2#1#c1\ta+,d+,b+\t*\nP\tS2#2#c1\ta+,d+,b+\t*\n'; } > "$1"; }
dip_gfa "$OUT/dip_agree.gfa" "a+,b+"; dip_gfa "$OUT/dip_differ.gfa" "a+,d+,b+"
printf 'bubble_id,source,source_orient,sink,sink_orient,inside_node_count,total_node_count,path_support,distinct_alleles,ref_allele_support,alt_allele_support_max,alt_allele_support_min,min_inside_bp,max_inside_bp,inside_nodes\n1,a,+,b,+,1,3,5,2,3,2,2,0,120,"d"\n' > "$OUT/dip.csv"
for o in agree differ; do "$BIN" call -i "$OUT/dip_$o.gfa" -c "$OUT/dip.csv" -r 'ref#0#chr1:1-520' -o "$OUT/dip_$o" --diploid \
  --no-per-bubble-vcf -q >/dev/null 2>&1; done
check "--diploid: agreeing contigs carry the record, and so does a carrying contig beside one that is not (AC agree | AC, S1 GT differ)" \
  "$(grep -v '^#' "$OUT/dip_agree.region.vcf" | grep -o ';AC=[0-9]*' | head -1)|$(grep -v '^#' "$OUT/dip_differ.region.vcf" | grep -o ';AC=[0-9]*' | head -1) $(awk -F'\t' '/^#CHROM/{for(i=10;i<=NF;i++) if($i=="S1") c=i} !/^#/{split($c,g,":"); print g[1]; exit}' "$OUT/dip_differ.region.vcf")" ";AC=1|;AC=1 1"
# Per path the same counting: AN and AC over haplotypes, the reference path left out. S1#1's two contigs are
# one haplotype (AN 3: S1#1, S2#1, S2#2), a carrier whether both or one of them carry the record.
for o in agree differ; do "$BIN" call -i "$OUT/dip_$o.gfa" -c "$OUT/dip.csv" -r 'ref#0#chr1:1-520' -o "$OUT/pp_$o" \
  --no-per-bubble-vcf -q >/dev/null 2>&1; done
check "per path: a haplotype in two contigs counts once, the reference not at all (AN;AC agree | differ)" \
  "$(grep -v '^#' "$OUT/pp_agree.region.vcf" | grep -o 'AN=[0-9]*;AC=[0-9]*' | head -1)|$(grep -v '^#' "$OUT/pp_differ.region.vcf" | grep -o 'AN=[0-9]*;AC=[0-9]*' | head -1)" \
  "AN=3;AC=1|AN=3;AC=1"

# Events right after an N run in the reference: the source node ends in 200 N, A (120 bp) is deleted by h2/h3 and
# I (80 bp) inserted before A by h4/h5. Each breakpoint is known only to within the run, and both records were
# written as precise with N as their anchor base.
NN=$(printf 'N%.0s' $(seq 1 200))
{ printf 'H\tVN:Z:1.0\nS\ts\t%s%s\nS\tA\t%s\nS\tI\t%s\nS\tt\t%s\n' "$(blk 40 71)" "$NN" "$(blk 120 72)" "$(blk 80 73)" "$(blk 200 74)"
  printf 'L\ts\t+\tA\t+\t0M\nL\tA\t+\tt\t+\t0M\nL\ts\t+\tt\t+\t0M\nL\ts\t+\tI\t+\t0M\nL\tI\t+\tA\t+\t0M\n'
  printf 'P\tref#0#chr1:1-560\ts+,A+,t+\t*\nP\th1#1#c\ts+,A+,t+\t*\nP\th2#1#c\ts+,t+\t*\nP\th3#1#c\ts+,t+\t*\n'
  printf 'P\th4#1#c\ts+,I+,A+,t+\t*\nP\th5#1#c\ts+,I+,A+,t+\t*\n'; } > "$OUT/nadj.gfa"
printf 'bubble_id,source,source_orient,sink,sink_orient,inside_node_count,total_node_count,path_support,distinct_alleles,ref_allele_support,alt_allele_support_max,alt_allele_support_min,min_inside_bp,max_inside_bp,inside_nodes\n1,s,+,t,+,2,4,5,3,2,2,2,0,200,"A;I"\n' > "$OUT/nadj.csv"
"$BIN" call -i "$OUT/nadj.gfa" -c "$OUT/nadj.csv" -r 'ref#0#chr1:1-560' -o "$OUT/nadj" --no-per-bubble-vcf -q >/dev/null 2>&1
check "events against an N run in the reference are IMPRECISE (DEL | INS)" \
  "$(grep -v '^#' "$OUT/nadj.region.vcf" | awk -F'\t' '/SVTYPE=DEL/{print ($8 ~ /IMPRECISE/) ? "imprecise" : "precise"}')|$(grep -v '^#' "$OUT/nadj.region.vcf" | awk -F'\t' '/SVTYPE=INS/{print ($8 ~ /IMPRECISE/) ? "imprecise" : "precise"}')" "imprecise|imprecise"
check "and their start can lie anywhere in the 200 N run (DEL CIPOS CIEND | INS CIPOS)" \
  "$(grep -v '^#' "$OUT/nadj.region.vcf" | awk -F'\t' '/SVTYPE=DEL/{n=split($8,a,";"); for(i=1;i<=n;i++){split(a[i],kv,"="); if(kv[1]=="CIPOS") p=kv[2]; if(kv[1]=="CIEND") e=kv[2]} print p" "e}')|$(grep -v '^#' "$OUT/nadj.region.vcf" | awk -F'\t' '/SVTYPE=INS/{n=split($8,a,";"); for(i=1;i<=n;i++){split(a[i],kv,"="); if(kv[1]=="CIPOS") p=kv[2]} print p}')" "-200,0 0,0|-200,0"

# ---------------------------------------------------------------- merging next to copy-number records
# Fixtures from a line spec (tests/tools/spec_gfa.py). th prints one item per record, in POS order:
# SVTYPE, SVLEN, /CN_METHOD on a CNV, then each haplotype's GT (first character), in path order.
TOOLS="$(cd "$(dirname "$0")/.." && pwd)/tools"
th() {   # th <name> <spec-file> [call options...]
  local d="$OUT/th_$1" spec="$2"; shift 2
  python3 -I "$TOOLS/spec_gfa.py" "$spec" "$d"
  "$BIN" call -i "$d/g.gfa" -c "$d/g.bubbles.csv" -r "$(awk '$1 == "P" { print $2; exit }' "$spec")" -o "$d/v" \
    --no-per-bubble-vcf -q "$@" >/dev/null 2>&1
  grep -v '^#' "$d/v.region.vcf" | awk -F'\t' '{ t = l = m = ""; n = split($8, a, ";")
      for (i = 1; i <= n; i++) { split(a[i], kv, "="); if (kv[1] == "SVTYPE") t = kv[2]; if (kv[1] == "SVLEN") l = kv[2]
                                 if (kv[1] == "CN_METHOD") m = "/" kv[2] }
      g = ""; for (i = 10; i <= NF; i++) g = g substr($i, 1, 1); printf "%s%s%s:%s ", t, l, m, g }'; }
# Reference a(400) g(60) b(60) C L; hA loses a and b, hB and hC a, hD a and g, hE all three; h5 carries a
# 100 bp expansion of the 20 bp repeat L, a UNIT record. That record used to put the whole site on node
# comparison, whose rules have no length gate: hA's two deletions coalesced across g into one of 460 bp,
# hD's 460 bp deletion merged into the 400 bp one, and hA's 60 bp deletion was lost. Aligned now, as at a
# site without a copy-number record. And a haplotype whose own larger deletion takes in a record's bases
# has no allele there: hD at the 400 bp deletion, hE at all three, GT '.', not 0.
cat > "$OUT/th_nest.spec" <<'SPEC'
N S 2000 1
N a 400 2
N g 60 7
N b 60 3
N C 200 5
N L 20 4
N K 2000 6
P ref#0#chr1:0-4740 S+,a+,g+,b+,C+,L+,K+
P h1#1#chr1 S+,a+,g+,b+,C+,L+,K+
P hA#1#chr1 S+,g+,C+,L+,K+
P hB#1#chr1 S+,g+,b+,C+,L+,K+
P hC#1#chr1 S+,g+,b+,C+,L+,K+
P hD#1#chr1 S+,b+,C+,L+,K+
P hE#1#chr1 S+,C+,L+,K+
P h5#1#chr1 S+,a+,g+,b+,C+,L+,L+,L+,L+,L+,L+,K+
B S K a;g;b;C;L
SPEC
check "beside a UNIT record, deletions are aligned; a haplotype inside a larger deletion is '.' (records)" \
  "$(th nest "$OUT/th_nest.spec")" \
  "DEL-400:00111..0 DEL-460:000001.0 DEL-520:00000010 DEL-60:001000.0 CNV100/UNIT:00000001 "
# The same insertion after a and after b, 300 bp apart: two records, node by node too, where sharing their
# one node made them one.
cat > "$OUT/th_two.spec" <<'SPEC'
N S 2000 1
N a 200 2
N b 300 3
N c 200 4
N I 300 9
N K 2000 6
P ref#0#chr1:0-4700 S+,a+,b+,c+,K+
P hX#1#chr1 S+,a+,I+,b+,c+,K+
P hX2#1#chr1 S+,a+,I+,b+,c+,K+
P hY#1#chr1 S+,a+,b+,I+,c+,K+
B S K a;b;c;I
SPEC
check "one inserted node at two places is two INS records, aligned and node by node (records | records)" \
  "$(th two "$OUT/th_two.spec")| $(th twon "$OUT/th_two.spec" --no-align-alleles)" \
  "INS300:0110 INS300:0001 | INS300:0110 INS300:0001 "
# 300 bp deleted, and the same 300 bp replaced by 4: the haplotypes they make are 99% the same and their net
# changes differ by 4 bp, so they are one record, as a deletion and a substitution never were.
cat > "$OUT/th_subdel.spec" <<'SPEC'
N S 2000 1
N a 200 2
N D 300 3
N X 4 8
N b 200 4
N K 2000 6
P ref#0#chr1:0-4700 S+,a+,D+,b+,K+
P h1#1#chr1 S+,a+,b+,K+
P h2#1#chr1 S+,a+,b+,K+
P h3#1#chr1 S+,a+,X+,b+,K+
P h4#1#chr1 S+,a+,X+,b+,K+
B S K a;D;X;b
SPEC
check "a deletion and a substitution of the same bases are one record (records)" "$(th subdel "$OUT/th_subdel.spec")" "DEL-300:01111 "
# One base changed on a 100 bp node (Mp for M), beside a UNIT record (L): a 1 bp change, not a 100 bp SUB.
cat > "$OUT/th_sub1.spec" <<'SPEC'
N S 2000 1
N A 200 2
N M =TGGCCAGTAGATCTTCCCAACATAGCCTAGCTGGACATATTCACTAAACCGAACAATCTATCACCAAGCGAATCCAGAGAGTCTCATGATACCTGGAGGA
N Mp =TGGCCAGTAGATCTTCCCAACATAGCCTAGCTGGACATATCACTAAACCGAACAATCTATCACCAAGCGAATCCAGAGAGTCTCATGATACCTGGAGGA
N C 200 5
N L 20 4
N K 2000 6
P ref#0#chr1:0-4520 S+,A+,M+,C+,L+,K+
P h1#1#chr1 S+,A+,M+,C+,L+,K+
P h2#1#chr1 S+,A+,Mp+,C+,L+,K+
P h5#1#chr1 S+,A+,M+,C+,L+,L+,L+,L+,L+,L+,K+
B S K A;M;Mp;C;L
SPEC
check "a 1 bp change on a 100 bp node is no SV record, aligned or node by node (records | records)" \
  "$(th sub1 "$OUT/th_sub1.spec")| $(th sub1n "$OUT/th_sub1.spec" --no-align-alleles)" "CNV100/UNIT:0001 | CNV100/UNIT:0001 "

# ---------------------------------------------------------------- the module route (--cn)
# f1: module A (the reference visits it twice), h2 one more copy (U1 A); a 20 bp repeat L outside the module,
# expanded by 100 bp in h3 and h4. The module record no longer switches L's UNIT record off.
cat > "$OUT/th_u1.spec" <<'SPEC'
N S 100 1
N A 200 2
N U1 100 3
N L 20 4
N U2 100 5
N K 100 6
P ref#0#chr1:0-820 S+,A+,U1+,A+,L+,U2+,K+
P h1#1#chr1 S+,A+,U1+,A+,L+,U2+,K+
P h2#1#chr1 S+,A+,U1+,A+,U1+,A+,L+,U2+,K+
P h3#1#chr1 S+,A+,U1+,A+,L+,L+,L+,L+,L+,L+,U2+,K+
P h4#1#chr1 S+,A+,U1+,A+,L+,L+,L+,L+,L+,L+,U2+,K+
B S K A;U1;L;U2
SPEC
check "--cn module: a module record beside a repeat outside it, which keeps its UNIT record (records)" \
  "$(th u1 "$OUT/th_u1.spec" --cn)" "CNV/MODULE:00100 CNV100/UNIT:00011 "
# f2: copies A1 A2 or A1 B2 (a variant copy), U1 between copies. h2 holds three copies on B2, h3 four, h4 three
# on A2. MODULE_BP counted only the reference's revisited nodes and read h2 as 2 copies, its gain then a 300 bp
# INS; counted at the anchor (A1, in every copy) h2 has 3.
cat > "$OUT/th_u2.spec" <<'SPEC'
N S 2000 1
N A1 150 2
N A2 50 3
N B2 50 7
N U1 100 4
N K 2000 6
P ref#0#chr1:0-4500 S+,A1+,A2+,U1+,A1+,A2+,K+
P h1#1#chr1 S+,A1+,A2+,U1+,A1+,A2+,K+
P h2#1#chr1 S+,A1+,B2+,U1+,A1+,B2+,U1+,A1+,B2+,K+
P h3#1#chr1 S+,A1+,B2+,U1+,A1+,B2+,U1+,A1+,B2+,U1+,A1+,B2+,K+
P h4#1#chr1 S+,A1+,A2+,U1+,A1+,A2+,U1+,A1+,A2+,K+
B S K A1;A2;B2;U1
SPEC
u2_cn() { grep -v '^#' "$OUT/th_u2/v.region.vcf" | awk -F'\t' '/MODULE/ { for (i = 10; i <= NF; i++) { split($i, f, ":"); printf "%s", f[2] } }'; }
check "--cn module: copies of either type are counted at the anchor (records | CN per path)" \
  "$(th u2 "$OUT/th_u2.spec" --cn)| $(u2_cn)" "CNV/MODULE:00111 | 22343"
# f3: X Y duplicated in h2 (40 bp, under --min-sv-bp) and a novel 500 bp NN in h2 and h3. PEAK counted the novel
# node as duplicated content and took the insertion from h2; counted, no module record, and both keep it.
cat > "$OUT/th_u3.spec" <<'SPEC'
N S 2000 1
N X 20 2
N Y 20 3
N Z 100 4
N NN 500 5
N K 2000 6
P ref#0#chr1:0-4140 S+,X+,Y+,Z+,K+
P h1#1#chr1 S+,X+,Y+,Z+,K+
P h2#1#chr1 S+,X+,Y+,X+,Y+,Z+,NN+,K+
P h3#1#chr1 S+,X+,Y+,Z+,NN+,K+
B S K X;Y;Z;NN
SPEC
check "--cn module: a novel insertion is not duplicated content (records)" "$(th u3 "$OUT/th_u3.spec" --cn)" "INS500:0011 "
# f7: copies A L B, U between them; L is a 60 bp repeat inside each copy. hv expands L in one copy, hd and hd2
# carry a third copy. L's UNIT record used to take the bubble: a copy gain read as +60 bp of L and two INS.
cat > "$OUT/th_u7.spec" <<'SPEC'
N S 2000 1
N A 300 2
N L 60 3
N B 300 4
N U 100 5
N K 2000 6
P ref#0#chr1:0-5420 S+,A+,L+,B+,U+,A+,L+,B+,K+
P h1#1#chr1 S+,A+,L+,B+,U+,A+,L+,B+,K+
P hv#1#chr1 S+,A+,L+,L+,B+,U+,A+,L+,B+,K+
P hd#1#chr1 S+,A+,L+,B+,U+,A+,L+,B+,U+,A+,L+,B+,K+
P hd2#1#chr1 S+,A+,L+,B+,U+,A+,L+,B+,U+,A+,L+,B+,K+
B S K A;L;B;U
SPEC
check "--cn module: a copy gain is the module's, a repeat change inside one copy an INS (records)" \
  "$(th u7 "$OUT/th_u7.spec" --cn)" "CNV/MODULE:00011 INS60:00100 "

# ---------------------------------------------------------------- the independent re-check of call (findings fixed)
# A copy gain whose copy is of another type, beside a change inside a copy: hX holds long + short + short copies
# (A1 H A2 / A1 A2) against two long ones. Its copy INS (2,020 bp) restates its CN gain and goes; its loss of H
# (600 bp, on the module's nodes) is part of its module change and stays as the DEL it is. Measuring the
# module's change without taking H off it kept the copy INS beside the CNV (C4: 26.4 kb on 7 haplotypes).
cat > "$OUT/th_c4like.spec" <<'SPEC'
N S 2000 1
N A1 1000 2
N H 600 3
N A2 1000 4
N m 20 5
N K 2000 7
P ref#0#chr1:0-9220 S+,A1+,H+,A2+,m+,A1+,H+,A2+,K+
P h1#1#chr1 S+,A1+,H+,A2+,m+,A1+,H+,A2+,K+
P hS#1#chr1 S+,A1+,H+,A2+,m+,A1+,A2+,K+
P hX#1#chr1 S+,A1+,H+,A2+,m+,A1+,A2+,m+,A1+,A2+,K+
P hY#1#chr1 S+,A1+,H+,A2+,m+,A1+,H+,A2+,m+,A1+,H+,A2+,K+
B S K A1;H;A2;m
SPEC
R6_C4LIKE="$(th c4like "$OUT/th_c4like.spec" --cn)"
check "a copy of another type beside a change inside a copy: the copy INS goes, the DEL stays (records)" "$R6_C4LIKE" "CNV/MODULE:00011 DEL-600:00110 "
# Two duplications in one site: hX gains a copy of A (counted, the anchor) and, apart from it, one of B. Only
# an event holding a copy of the anchor restates the copy change: A's INS goes, B's 1,520 bp INS stays.
cat > "$OUT/th_twomod.spec" <<'SPEC'
N S 2000 1
N A 1000 2
N m 20 3
N n 30 4
N B 1500 5
N p 20 6
N K 2000 7
P ref#0#chr1:0-7550 S+,A+,m+,A+,n+,B+,K+
P h1#1#chr1 S+,A+,m+,A+,n+,B+,K+
P h2#1#chr1 S+,A+,m+,A+,n+,B+,K+
P hX#1#chr1 S+,A+,m+,A+,m+,A+,n+,B+,p+,B+,K+
B S K A;m;n;B;p
SPEC
R6_TWOMOD="$(th twomod "$OUT/th_twomod.spec" --cn)"
check "two duplications in one site: only the anchor's copy is restated; the other INS stays (records)" "$R6_TWOMOD" "CNV/MODULE:0001 INS1520:0001 "
# B is duplicated in every path (two copies each), A varies (1, 2 or 3 copies). A duplication the whole panel
# carries says nothing about copy number and is not the anchor; counted at it, A's changes had no record.
cat > "$OUT/th_fixdup.spec" <<'SPEC'
N S 2000 1
N A 1000 2
N m 20 3
N n 30 4
N B 1500 5
N p 20 6
N K 2000 7
P ref#0#chr1:0-9070 S+,A+,m+,A+,n+,B+,p+,B+,K+
P h1#1#chr1 S+,A+,m+,A+,n+,B+,p+,B+,K+
P h2#1#chr1 S+,A+,m+,A+,n+,B+,p+,B+,K+
P h3#1#chr1 S+,A+,m+,A+,n+,B+,p+,B+,K+
P h4#1#chr1 S+,A+,m+,A+,m+,A+,n+,B+,p+,B+,K+
P h5#1#chr1 S+,A+,m+,A+,m+,A+,n+,B+,p+,B+,K+
P h6#1#chr1 S+,A+,n+,B+,p+,B+,K+
B S K A;m;n;B;p
SPEC
R6_FIXDUP="$(th fixdup "$OUT/th_fixdup.spec" --cn)"
check "a duplication fixed in the panel is not the anchor; the varying module is counted (records)" "$R6_FIXDUP" "CNV/MODULE:0000111 "
# The reference holds no copy of the module A1 L L A2 (L a 30 bp loop): L is inside the copies, read from the
# haplotypes' copies, so it has no UNIT record of its own beside the module's (it counted the gain twice).
cat > "$OUT/th_ref0.spec" <<'SPEC'
N S 2000 1
N A1 600 2
N L 30 3
N A2 600 4
N m 20 5
N K 2000 7
P ref#0#chr1:0-4000 S+,K+
P h1#1#chr1 S+,A1+,L+,L+,A2+,K+
P h2#1#chr1 S+,A1+,L+,L+,A2+,m+,A1+,L+,L+,A2+,K+
P h3#1#chr1 S+,A1+,L+,L+,A2+,m+,A1+,L+,L+,A2+,K+
P h4#1#chr1 S+,K+
B S K A1;L;A2;m
SPEC
R6_REF0="$(th ref0 "$OUT/th_ref0.spec" --cn)"
check "no reference copy: a repeat inside the copies has no UNIT record of its own (records)" "$R6_REF0" "CNV/MODULE:01110 "
# A module (M) and, apart from it, a UNIT array of three 300 bp copies with 40 bp junctions; hX keeps one copy.
# The gap left where its copies leave the walks is split there as without a module record: no 80 bp DEL of the
# two junctions (IMPRECISE) beside the UNIT record.
cat > "$OUT/th_split.spec" <<'SPEC'
N S 2000 1
N M 1000 2
N n 20 3
N x 100 4
N U 300 5
N j1 40 6
N j2 40 8
N y 100 9
N K 2000 7
P ref#0#chr1:0-7480 S+,M+,n+,M+,x+,U+,j1+,U+,j2+,U+,y+,K+
P h1#1#chr1 S+,M+,n+,M+,x+,U+,j1+,U+,j2+,U+,y+,K+
P h2#1#chr1 S+,M+,n+,M+,n+,M+,x+,U+,j1+,U+,j2+,U+,y+,K+
P h3#1#chr1 S+,M+,n+,M+,n+,M+,x+,U+,j1+,U+,j2+,U+,y+,K+
P hX#1#chr1 S+,M+,n+,M+,x+,U+,y+,K+
B S K M;n;x;U;j1;j2;y
SPEC
cat > "$OUT/th_split.panphorte.unit_graph.tsv" <<'TSV'
input_bubble_id	unit_bp	anchor_pos	anchor_nodes	n_nodes	nodes
1	340	0	U	1	U
TSV
R6_SPLIT="$(th split "$OUT/th_split.spec" --cn --panphorte-prefix "$OUT/th_split")"
check "a module record beside a UNIT array: no DEL of the junctions between removed copies (records)" "$R6_SPLIT" "CNV/MODULE:00110 CNV-680/UNIT:00001 "
# Two arrays of one motif with 500 bp between them are two repeats: hX gains a copy on the first and loses one
# on the second. Merged as one unit (one site, one motif) the two cancelled and hX read as unchanged.
cat > "$OUT/th_motif.spec" <<'SPEC'
N S 2000 1
N a1 40 2
N b1 20 3
N M 500 4
N a2 40 5
N b2 20 6
N K 2000 7
P ref#0#chr1:0-4740 S+,a1+,b1+,a1+,b1+,M+,a2+,b2+,a2+,b2+,K+
P h1#1#chr1 S+,a1+,b1+,a1+,b1+,M+,a2+,b2+,a2+,b2+,K+
P hX#1#chr1 S+,a1+,b1+,a1+,b1+,a1+,b1+,M+,a2+,b2+,K+
P hY#1#chr1 S+,a1+,b1+,a1+,b1+,a1+,b1+,M+,a2+,b2+,a2+,b2+,K+
B S K a1;b1;M;a2;b2
SPEC
cat > "$OUT/th_motif.panphorte.unit_graph.tsv" <<'TSV'
input_bubble_id	unit_bp	anchor_pos	anchor_nodes	n_nodes	nodes	canonical_motif
1	60	0	a1	2	a1;b1	ACGT
1	60	0	a2	2	a2;b2	ACGT
TSV
R6_MOTIF="$(th motif "$OUT/th_motif.spec" --panphorte-prefix "$OUT/th_motif")"
check "two arrays of one motif apart are two repeats (records)" "$R6_MOTIF" "CNV60/UNIT:0011 CNV-60/UNIT:0010 "
# A reverse-oriented site compared node by node: h2 deletes a and b, taking in the point where h1 inserts X.
# h2 has no allele at the insertion ('.'), as in a forward site.
cat > "$OUT/th_insr.spec" <<'SPEC'
N S 2000 1
N a 300 2
N X 120 3
N b 300 4
N c 300 5
N K 2000 7
P ref#0#chr1:0-4900 K-,c-,b-,a-,S-
P h1#1#chr1 K-,c-,b-,X-,a-,S-
P h2#1#chr1 K-,c-,S-
P h3#1#chr1 K-,c-,b-,a-,S-
B S K a;X;b;c
SPEC
R6_INSR="$(th insr "$OUT/th_insr.spec" --no-align-alleles)"
check "a reverse site node by node: a haplotype deleting across an insertion point has no allele there (records)" "$R6_INSR" "DEL-600:0010 INS120:01.0 "
# A module spanning more than 80% of the reference is no copy number: this one spans 4,040 of 4,240 reference
# bp (95%), so it gets no CNV, its copy changes are the INS and DEL, and a warning names the site even with -q.
cat > "$OUT/th_oversz.spec" <<'SPEC'
N S 100 1
N A 2000 2
N m 20 3
N K 100 7
P ref#0#chr1:0-4240 S+,A+,m+,A+,m+,K+
P h1#1#chr1 S+,A+,m+,A+,m+,K+
P h2#1#chr1 S+,A+,m+,K+
P h3#1#chr1 S+,A+,m+,A+,m+,A+,m+,K+
B S K A;m
SPEC
R6_OVERSZ="$(th oversz "$OUT/th_oversz.spec" --cn)"
check "a module over 80% of the reference: no CNV over 95% of it (records)" "$R6_OVERSZ" "INS2020:0001 DEL-2020:0010 "
oversz_log=$("$BIN" call -i "$OUT/th_oversz/g.gfa" -c "$OUT/th_oversz/g.bubbles.csv" -r ref -o "$OUT/th_oversz/w" --cn -q 2>&1 >/dev/null)
printf '%s' "$oversz_log" | grep -q "no module copy number at site(s) 1: .*panvar rebuild" \
  && ok "a module over 80% of the reference is named in a warning pointing at rebuild, even with -q" \
  || bad "oversized module: expected a warning naming site 1, got: $oversz_log"
# Node by node, a deletion and an insertion 3 bp apart whose length changes cancel are one substitution. h1
# and h2 carry X again 3 bp later with 4 bases changed (Y): the same haplotype bar a few bases, which was a
# 475 bp DEL plus a 475 bp INS. h3 carries unrelated sequence there (Z): one SUB of the two, no length change.
cat > "$OUT/th_join.spec" <<'SPEC'
N S 100 1
N X 475 2
N t 3 3
N Y @X:4:5
N Z 475 6
N K 100 4
P ref#0#chr1:0-678 S+,X+,t+,K+
P h1#1#chr1 S+,t+,Y+,K+
P h2#1#chr1 S+,t+,Y+,K+
P h3#1#chr1 S+,t+,Z+,K+
B S K X;t;Y;Z
SPEC
check "a deletion and an insertion 3 bp apart that cancel are one substitution, node by node (records)" \
  "$(th join "$OUT/th_join.spec" --no-align-alleles)" "SUB0:0001 "
# Pieces are joined only where they lie in one place on both walks. Here the UNIT copies of U left the walks
# between D's deletion and I's insertion, so the tokens t U u look adjacent although a whole array lies
# between them: joined, they made a substitution whose REF spelled bases the reference does not have there.
cat > "$OUT/th_joinu.spec" <<'SPEC'
N a 200 1
N D 500 2
N t 10 3
N U 100 4
N u 10 5
N I 500 6
N b 200 7
P ref a+,D+,t+,U+,U+,U+,u+,b+
P h1 a+,t+,U+,u+,I+,b+
P h2 a+,t+,U+,u+,I+,b+
P h3 a+,D+,t+,U+,U+,U+,u+,b+
B a b D;t;U;u;I
SPEC
check "no join across copies a UNIT record took out of the walks (records)" \
  "$(th joinu "$OUT/th_joinu.spec")" "DEL-500:0110 CNV-200/UNIT:0110 INS500:0110 "
# A small piece does not join a large one: h1 and h2 delete D and insert 5 bp 10 bp further on, h3 and h4
# only delete D. One DEL for all four; joined, h1 and h2 had a substitution of their own and GT . on the DEL.
cat > "$OUT/th_joins.spec" <<'SPEC'
N a 200 1
N D 500 2
N t 10 3
N s 5 8
N u 30 5
N b 200 7
P ref a+,D+,t+,u+,b+
P h1 a+,t+,s+,u+,b+
P h2 a+,t+,s+,u+,b+
P h3 a+,t+,u+,b+
P h4 a+,t+,u+,b+
B a b D;t;s;u
SPEC
check "a 5 bp insertion beside a 500 bp deletion leaves one DEL for all carriers, node by node (records)" \
  "$(th joins "$OUT/th_joins.spec" --no-align-alleles)" "DEL-500:01111 "
# Pairs, not chains: D replaced by I is one substitution for all four, whether or not E is also lost.
cat > "$OUT/th_joinp.spec" <<'SPEC'
N a 200 1
N D 300 2
N t 10 3
N I 300 6
N u 10 5
N E 600 8
N b 200 7
P ref a+,D+,t+,u+,E+,b+
P h1 a+,t+,I+,u+,b+
P h2 a+,t+,I+,u+,b+
P h3 a+,t+,I+,u+,E+,b+
P h4 a+,t+,I+,u+,E+,b+
B a b D;t;I;u;E
SPEC
check "a replacement next to a deletion is one substitution for every carrier, node by node (records)" \
  "$(th joinp "$OUT/th_joinp.spec" --no-align-alleles)" "SUB0:01111 DEL-600:01100 "


# ---------------------------------------------------------------- per-gene copy number (--gtf), 2026-10-31
# gs <name> <spec> <GENE:node>...: the spec's reference renamed ref#0#chr9:1-LEN, a GTF with one gene (gene and
# CDS) over each named node's first reference occurrence, then call --cn --gtf. Prints dup_gene_cn rows as
# variant:sample:gene:cn:reliable, one per row.
gs() {
  local d="$OUT/gs_$1" spec="$2"; shift 2
  python3 -I "$TOOLS/spec_gfa.py" "$spec" "$d" >/dev/null
  python3 -I - "$d" "$@" <<'GPY'
import os, sys
d, genes = sys.argv[1], sys.argv[2:]
lines = open(os.path.join(d, 'g.gfa')).read().split('\n')
seg, ref = {}, None
for l in lines:
    f = l.split('\t')
    if f[0] == 'S': seg[f[1]] = f[2]
    if f[0] == 'P' and f[1] == 'ref': ref = f[2].split(',')
name = 'ref#0#chr9:1-%d' % sum(len(seg[s[:-1]]) for s in ref)
open(os.path.join(d, 'g.gfa'), 'w').write('\n'.join(l.replace('P\tref\t', 'P\t%s\t' % name, 1) for l in lines))
pos, at = {}, 0
for s in ref:
    pos.setdefault(s[:-1], at + 1); at += len(seg[s[:-1]])
with open(os.path.join(d, 'g.gtf'), 'w') as f:
    for g in genes:
        nm, node = g.split(':'); a = pos[node]; b = a + len(seg[node]) - 1
        for kind in ('gene', 'CDS'):
            f.write('chr9\tt\t%s\t%d\t%d\t.\t+\t.\tgene_id "%s"; gene_name "%s";\n' % (kind, a, b, nm, nm))
GPY
  "$BIN" call -i "$d/g.gfa" -c "$d/g.bubbles.csv" -r ref -o "$d/v" --cn --gtf "$d/g.gtf" --no-per-bubble-vcf -q >/dev/null 2>&1
  tail -n +2 "$d/v.dup_gene_cn.tsv" 2>/dev/null | awk -F'\t' '{split($3, s, "#"); printf "%s:%s:%s:%s:%s ", $2, s[1], $4, $5, $6}'
}
# A gene beside a module (GENEZ, single-copy, before two-copy A-l1-C-lk) is not the module's copies: it
# was given the module's 2/1/3/2 copies, marked reliable. No row now.
cat > "$OUT/gs_single.spec" <<'SPEC'
N f1 4000 1
N gz 1500 21
N sp 200 22
N gA 1200 2
N l1 300 4
N gC 1200 8
N lk 300 6
N f2 4000 7
P ref f1+,gz+,sp+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,f2+
P h1 f1+,gz+,sp+,gA+,l1+,gC+,lk+,f2+
P h2 f1+,gz+,sp+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,f2+
P h3 f1+,gz+,sp+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,f2+
B f1 f2 gz;sp;gA;l1;gC;lk
SPEC
check "a gene beside a module's copies gets no copy number from it" "$(gs single "$OUT/gs_single.spec" GENEZ:gz)" ""
# A near-identical pair, one copy each, with a 60 bp repeat between them (a UNIT record, 6 units on the
# reference). The pair was split by the repeat's unit count (h3 5/4). Each gene reads its own k-mers: 1.
cat > "$OUT/gs_unitpair.spec" <<'SPEC'
N f1 500 1
N gA 1500 2
N gB @gA:8:3
N l1 300 4
N u 60 5
N l2 300 6
N f2 500 7
P ref f1+,gA+,l1+,u+,u+,u+,u+,u+,u+,l2+,gB+,f2+
P h1 f1+,gA+,l1+,u+,u+,u+,l2+,gB+,f2+
P h2 f1+,gA+,l1+,u+,u+,u+,u+,l2+,gB+,f2+
P h3 f1+,gA+,l1+,u+,u+,u+,u+,u+,u+,u+,u+,u+,l2+,gB+,f2+
P h4 f1+,gA+,l1+,u+,u+,u+,u+,u+,u+,l2+,gB+,f2+
B f1 f2 gA;gB;l1;u;l2
SPEC
GSU="$(gs unitpair "$OUT/gs_unitpair.spec" GENEA:gA GENEB:gB)"
check "a pair beside a repeat is not split by the repeat's unit count (h3 GENEA, GENEB)" \
      "$(printf '%s' "$GSU" | tr ' ' '\n' | awk -F: '$2=="h3"{print $3"="$4}' | paste -sd' ' -)" "GENEA=1 GENEB=1"
# A tandem copy holding a whole gene, folded into a repeat unit, does count the gene: the unit is its copy.
cat > "$OUT/gs_unitgene.spec" <<'SPEC'
N f1 500 1
N g 900 2
N f2 500 7
P ref f1+,g+,g+,f2+
P h1 f1+,g+,f2+
P h3 f1+,g+,g+,g+,f2+
B f1 f2 g
SPEC
check "a gene inside a repeat unit is counted by it (ref, h1, h3)" \
      "$(gs unitgene "$OUT/gs_unitgene.spec" GENEG:g | tr ' ' '\n' | awk -F: 'NF{print $2"="$4}' | paste -sd' ' -)" "ref=2 h1=1 h3=3"

# call removes the optional outputs an earlier run at the same prefix wrote and this one does not.
"$BIN" call -i "$OUT/th_nest/g.gfa" -c "$OUT/th_nest/g.bubbles.csv" -r "$(awk '$1 == "P" { print $2; exit }' "$OUT/th_nest.spec")" \
       -o "$OUT/stale" --no-per-bubble-vcf -q >/dev/null 2>&1
"$BIN" call -i "$OUT/th_nest/g.gfa" -c "$OUT/th_nest/g.bubbles.csv" -r "$(awk '$1 == "P" { print $2; exit }' "$OUT/th_nest.spec")" \
       -o "$OUT/stale" --no-variant-nodes --no-per-bubble-vcf -q >/dev/null 2>&1
check "outputs of an earlier run with other options are removed (variant nodes, carrier events)" \
      "$(ls "$OUT"/stale.* | sed 's|.*/||' | paste -sd' ' -)" "stale.region.vcf"

# ---------------------------------------------------------------- inversion sites, 2026-11-01
# A large inversion between inverted repeats: every path walks the repeat nodes in both orientations, so no
# snarl holds it (the snarls at the breakpoints have no interior a path crosses). bubble makes it a site
# from the walks -- kind 'inversion', between flanks every path crosses once forward -- and call types it
# there: a run of reference-unique nodes walked reversed, between two runs walked forward. Two graph forms:
# the repeat collapsed into one node the reference walks twice (pggb), and its two copies kept apart
# (rebuild). hR is stored on the other strand: reversed everywhere, no inversion.
thb() {   # thb <name> <spec-file> [call options...]: as th, with the sites from panvar bubble
  local d="$OUT/th_$1" spec="$2" ref; shift 2
  ref="$(awk '$1 == "P" { print $2; exit }' "$spec")"
  python3 -I "$TOOLS/spec_gfa.py" "$spec" "$d"
  "$BIN" bubble -i "$d/g.gfa" -r "$ref" -o "$d/b" -q >/dev/null 2>&1
  "$BIN" call -i "$d/b.sorted.gfa" -b "$d/b" -r "$ref" -o "$d/v" --no-per-bubble-vcf -q "$@" >/dev/null 2>&1
  grep -v '^#' "$d/v.region.vcf" | awk -F'\t' '{ t = l = ""; n = split($8, a, ";")
      for (i = 1; i <= n; i++) { split(a[i], kv, "="); if (kv[1] == "SVTYPE") t = kv[2]; if (kv[1] == "SVLEN") l = kv[2] }
      g = ""; for (i = 10; i <= NF; i++) g = g substr($i, 1, 1); printf "%s%s:%s ", t, l, g }'; }
cat > "$OUT/th_invcol.spec" <<'SPEC'
N S 2000 1
N I 300 2
N M1 3000 3
N M2 3000 4
N M3 3000 5
N T 2000 6
P ref#0#chr1:0-13600 S+,I+,M1+,M2+,M3+,I-,T+
P h1#1#chr1 S+,I+,M1+,M2+,M3+,I-,T+
P hI#1#chr1 S+,I+,M3-,M2-,M1-,I-,T+
P hJ#1#chr1 S+,I+,M3-,M2-,M1-,I-,T+
P hR#1#chr1 T-,I+,M3-,M2-,M1-,I-,S-
SPEC
cat > "$OUT/th_invsep.spec" <<'SPEC'
N S 2000 1
N A 300 2
N M1 3000 3
N M2 3000 4
N M3 3000 5
N B 300 7
N T 2000 6
P ref#0#chr1:0-13600 S+,A+,M1+,M2+,M3+,B+,T+
P h1#1#chr1 S+,A+,M1+,M2+,M3+,B+,T+
P hI#1#chr1 S+,B-,M3-,M2-,M1-,B+,T+
P hR#1#chr1 T-,B-,M3-,M2-,M1-,A-,S-
SPEC
check "an inversion between collapsed inverted repeats is an inversion site, typed from the walks (records)" \
  "$(thb invcol "$OUT/th_invcol.spec")" "INV9600:00110 "
check "its site runs between the flanks, kind 'inversion'" \
  "$(awk -F, 'NR > 1 { print $2 ".." $4 ":" $16 }' "$OUT/th_invcol/b.bubbles.csv" | paste -sd' ' -)" "1..6:inversion"
check "and between inverted repeats kept as two copies (records)" "$(thb invsep "$OUT/th_invsep.spec")" "INV9600:0010 "
check "its breakpoints are the repeats: IMPRECISE, CIPOS 0,300, CIEND -300,0, at its site" \
  "$(grep -v '^#' "$OUT/th_invcol/v.region.vcf" | grep -o 'PATH_INVERSION;BUBBLE_ID=1;IMPRECISE;CIPOS=0,300;CIEND=-300,0' | head -1)" \
  "PATH_INVERSION;BUBBLE_ID=1;IMPRECISE;CIPOS=0,300;CIEND=-300,0"
check "its variant nodes are the site's own" "$(tail -n +2 "$OUT/th_invcol/v.variant_nodes.tsv" | cut -f1-4)" \
  "$(printf 'bubble1_INV_1\t1\tINV\t2,3,4,5')"
# A site inside the inversion keeps its own variant and names the inversion site as its parent: h1 and
# hI lose X (60 bp) inside M2; hI also carries the inversion. The ordinary site keeps id 1 (node ids:
# bubble's sorted graph).
cat > "$OUT/th_invnest.spec" <<'SPEC'
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
check "a site inside an inversion site keeps its variant; both are typed (records)" \
  "$(thb invnest "$OUT/th_invnest.spec")" "INV9660:00110 DEL-60:01100 "
check "the inner site names the inversion site as its parent and keeps its id" \
  "$(awk -F, 'NR > 1 { print $1 ":" $2 ".." $4 ":" $16 ":" $17 }' "$OUT/th_invnest/b.bubbles.csv" | paste -sd' ' -)" \
  "1:4..6:site:2 2:1..8:inversion:"
# An inversion a site holds is typed there, and not written again from the walks.
cat > "$OUT/th_invsite.spec" <<'SPEC'
N S 2000 1
N M 3000 3
N T 2000 6
P ref#0#chr1:0-7000 S+,M+,T+
P h1#1#chr1 S+,M+,T+
P hI#1#chr1 S+,M-,T+
B S T M
SPEC
check "an inversion a site holds is one record, not two (records)" "$(thb invsite "$OUT/th_invsite.spec")" "INV3000:001 "
check "and no inversion site is made for it" "$(awk -F, 'NR > 1 { print $16 }' "$OUT/th_invsite/b.bubbles.csv" | paste -sd' ' -)" "site"

# ---------------------------------------------------------------------------------------------
# A copy loss beside a novel insertion at a collapsed module. h1 and h2 lose one gA..lk copy and carry X,
# h5 carries X only. The node comparison sees the lost copy and X in one gap; written as one replacement,
# the copy change took it whole and X was lost for h1 and h2. Split, the copy goes to the CNV and X stays.
cat > "$OUT/th_subcn.spec" <<'SPEC'
N f1 2000 1
N t 10 9
N X 1000 10
N gA 1200 2
N l1 300 4
N gC 1200 8
N lk 300 6
N sp 200 11
N f2 2000 7
P ref#0#chr1:0-9010 f1+,t+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,sp+,f2+
P h1#1#chr1 f1+,t+,X+,gA+,l1+,gC+,lk+,sp+,f2+
P h2#1#chr1 f1+,t+,X+,gA+,l1+,gC+,lk+,sp+,f2+
P h3#1#chr1 f1+,t+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,sp+,f2+
P h4#1#chr1 f1+,t+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,sp+,f2+
P h5#1#chr1 f1+,t+,X+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,sp+,f2+
B f1 f2 t;X;gA;l1;gC;lk;sp
SPEC
check "a copy loss beside a novel insertion: the insertion keeps every carrier (records)" \
  "$(th subcn "$OUT/th_subcn.spec" --cn)" "CNV/MODULE:011010 INS1000:011001 "

# Two deletions 500 bp apart on the reference, whose nodes the haplotype walks side by side elsewhere (the
# reference holds A and B twice): two
# records, also node by node. Joined through the haplotype's own sequence, they were one 120 bp deletion.
cat > "$OUT/th_farcoal.spec" <<'SPEC'
N S 100 1
N A 60 2
N M 500 3
N B 60 4
N N 500 5
N K 100 6
P ref#0#chr1:1-1440 S+,A+,M+,B+,N+,A+,B+,K+
P h1#1#chr1 S+,M+,N+,A+,B+,K+
B S K A;M;B;N
SPEC
check "two deletions far apart on the reference stay two, node by node (records)" \
  "$(th farcoal "$OUT/th_farcoal.spec" --no-align-alleles)" "DEL-60:01 DEL-60:01 "

# Equal events on different nodes: the representative is the same whatever the nodes are numbered.
for v in a b; do
  if [ "$v" = a ]; then n2=2; n3=3; else n2=30; n3=20; fi
  cat > "$OUT/th_renum_$v.spec" <<SPEC
N 1 100 1
N $n2 60 2
N $n3 @$n2:0:3
N 5 100 5
P ref#0#chr1:1-200 1+,5+
P hapX#1#chr1 1+,$n2+,5+
P hapY#1#chr1 1+,$n3+,5+
B 1 5 $n2;$n3
SPEC
  th "renum_$v" "$OUT/th_renum_$v.spec" >/dev/null
done
check "the representative of equal events does not depend on node numbering (EVENT_NODES)" \
  "$(for v in a b; do grep -v '^#' "$OUT/th_renum_$v/v.region.vcf" | grep -o 'EVENT_NODES=[^;]*'; done | paste -sd' ' -)" \
  "EVENT_NODES=2 EVENT_NODES=30"

# GENES from the record's reference interval: an insertion inside a gene names it (it has no reference
# node of its own), one between two genes names neither, and a name with ';' is percent-encoded.
cat > "$OUT/tg.spec" <<'SPEC'
N f1 500 1
N X 200 6
N m 500 2
N Y 200 5
N f2 500 7
P ref#0#chr9:1-1500 f1+,m+,f2+
P h1#1#chr9 f1+,X+,m+,f2+
P h2#1#chr9 f1+,m+,Y+,f2+
P h3#1#chr9 f1+,m+,f2+
B f1 f2 X;m;Y
SPEC
python3 -I "$TOOLS/spec_gfa.py" "$OUT/tg.spec" "$OUT/tg"
printf 'chr9\tt\tgene\t300\t700\t.\t+\t.\tgene_id "G1"; gene_name "IN;SIDE";\nchr9\tt\tgene\t1300\t1400\t.\t+\t.\tgene_id "G2"; gene_name "LATER";\n' > "$OUT/tg/g.gtf"
"$BIN" call -i "$OUT/tg/g.gfa" -c "$OUT/tg/g.bubbles.csv" -r 'ref#0#chr9:1-1500' -o "$OUT/tg/v" --gtf "$OUT/tg/g.gtf" \
  --no-per-bubble-vcf -q >/dev/null 2>&1
check "GENES of an insertion inside a gene, and of one between genes (POS:GENES)" \
  "$(grep -v '^#' "$OUT/tg/v.region.vcf" | awk -F'\t' '{ g = "-"; if (match($8, /GENES=[^;]*/)) g = substr($8, RSTART + 6, RLENGTH - 6); printf "%s:%s ", $2, g }')" \
  "500:IN%3BSIDE 1000:- "

# An input is never removed as a stale output: a GTF named <prefix>.node_genes.tsv with no gene on the
# reference's chromosome stays.
cp "$OUT/tg/g.gtf" "$OUT/tg/w.node_genes.tsv"
sed -i.bak 's/^chr9/chr2/' "$OUT/tg/w.node_genes.tsv"
"$BIN" call -i "$OUT/tg/g.gfa" -c "$OUT/tg/g.bubbles.csv" -r 'ref#0#chr9:1-1500' -o "$OUT/tg/w" --gtf "$OUT/tg/w.node_genes.tsv" \
  --no-per-bubble-vcf -q >/dev/null 2>&1
check "a GTF named like a stale gene table is not removed" "$([ -s "$OUT/tg/w.node_genes.tsv" ] && echo kept || echo removed)" "kept"

# An inversion-site record passes the same carrier filters as every other record.
check "an inversion-site record carried by one haplotype, with --min-haplotypes 2 (records)" \
  "$(thb invsep2 "$OUT/th_invsep.spec" --min-haplotypes 2)" ""

# Hybrid copies are read on each copy's own strand: with GENEB's copies walked reverse-complemented, the
# hybrid fixture reads the same.
sed -e "s/    return st + \['lk+'\]/    if word == 'bbbbbb': st = [x[:-1] + '-' for x in reversed(st)]\n    return st + ['lk+']/" \
    -e "s/^seq = lambda st: .*/rc = lambda q: q[::-1].translate(str.maketrans('ACGT', 'TGCA'))\nseq = lambda st: ''.join(segs[s[:-1]] if s[-1] == '+' else rc(segs[s[:-1]]) for s in st)/" \
    -e "s/f.write('L\\\\t%s\\\\t+\\\\t%s\\\\t+\\\\t0M\\\\n' % (x\[:-1\], y\[:-1\]))/f.write('L\\\\t%s\\\\t%s\\\\t%s\\\\t%s\\\\t0M\\\\n' % (x[:-1], x[-1], y[:-1], y[-1]))/" \
    "$OUT/hyb_gen.py" > "$OUT/hyb_inv_gen.py"
python3 "$OUT/hyb_inv_gen.py" "$OUT/hybi"
"$BIN" call -i "$OUT/hybi.gfa" -c "$OUT/hybi.csv" -r 'ref#0#chr9' -o "$OUT/hybi" --gtf "$OUT/hybi.gtf" --cn -q >/dev/null 2>&1
hybi_of() { awk -F'\t' -v h="$1" -v g="$2" 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "hybrid_copies") c = i }
  $3 ~ ("^" h "#") && $4 == g { print $5 "/" $c; exit }' "$OUT/hybi.dup_gene_cn.tsv"; }
check "hybrid copies with the GENEB copies reversed: cn/hybrids for h2, hhyb, hconv" \
  "$(hybi_of h2 GENEA) $(hybi_of h2 GENEB) $(hybi_of hhyb GENEA) $(hybi_of hhyb GENEB) $(hybi_of hconv GENEA) $(hybi_of hconv GENEB)" \
  "2/0 1/0 2/1 1/1 2/0 1/0"

# The tangle gate counts hubs as rebuild's gate does: fifty neighbours or more on one end. Ten interior
# nodes each linked to 50 leaves on one end make a tangle: no module record, and a warning even with -q.
# With 25 leaves on each end they are clean branches, and the module is counted.
python3 -I - "$OUT/tangle" <<'PY'
import os, random, sys
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
rs = lambda n, s: ''.join(random.Random(s).choice('ACGT') for _ in range(n))
for name, per_end in (('one', (50, 0)), ('two', (25, 25))):
    seg = {'S': rs(500, 1), 'M': rs(300, 2), 'x': rs(100, 3), 'T': rs(500, 4)}
    links = set()
    hubs = ['H%d' % i for i in range(10)]
    for i, h in enumerate(hubs):
        seg[h] = rs(40, 10 + i)
        links.add((h, '+', 'x', '+'))
        for side, n in zip('+-', per_end):
            for j in range(n):
                leaf = 'L%d_%s%d' % (i, side, j); seg[leaf] = rs(5, 1000 + 100 * i + j + (50 if side == '-' else 0))
                links.add((h, side, leaf, '+'))
    paths = {'ref#0#chr1:1-1700': 'S M x M T', 'h1#1#chr1': 'S M x M x M T', 'h2#1#chr1': 'S M x M x M T',
             'h3#1#chr1': 'S M T', 'h4#1#chr1': 'S M x M T'}
    for st in paths.values():
        w = st.split()
        links.update((a, '+', b, '+') for a, b in zip(w, w[1:]))
    with open(os.path.join(out, name + '.gfa'), 'w') as f:
        for k, v in seg.items(): f.write('S\t%s\t%s\n' % (k, v))
        for a, ao, b, bo in sorted(links): f.write('L\t%s\t%s\t%s\t%s\t0M\n' % (a, ao, b, bo))
        for n, st in paths.items(): f.write('P\t%s\t%s\t*\n' % (n, ','.join(s + '+' for s in st.split())))
    with open(os.path.join(out, name + '.csv'), 'w') as f:
        f.write('bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n')
        f.write('1,S,+,T,+,5,300,1100,"%s"\n' % ';'.join(['M', 'x'] + hubs))
PY
for v in one two; do
  "$BIN" call -i "$OUT/tangle/$v.gfa" -c "$OUT/tangle/$v.csv" -r 'ref#0#chr1:1-1700' -o "$OUT/tangle/$v" --cn \
    --no-per-bubble-vcf -q > /dev/null 2> "$OUT/tangle/$v.err"
done
check "hubs on one end make a tangle, hubs spread over both ends do not (module records, warnings)" \
  "$(grep -v '^#' "$OUT/tangle/one.region.vcf" | grep -c 'CN_METHOD=MODULE') $(grep -c 'tangle' "$OUT/tangle/one.err") $(grep -v '^#' "$OUT/tangle/two.region.vcf" | grep -c 'CN_METHOD=MODULE') $(grep -c 'tangle' "$OUT/tangle/two.err")" \
  "0 1 1 0"

# ---------------------------------------------------------------------------------------------
# The split of a copy change from a new insertion beside it, cut where the two walks share a node: h1 and h2
# lose a copy and carry X after the 20 bp node s, which a join had taken into one replacement with the lost
# copy. The insertion is h5's X, at the same place, for all three.
cat > "$OUT/th_splitj.spec" <<'SPEC'
N f1 2000 1
N t 10 9
N X 1000 10
N s 20 14
N gA 1200 2
N l1 300 4
N gC 1200 8
N lk 300 6
N sp 200 11
N f2 2000 7
P ref#0#chr1:0-9030 f1+,t+,gA+,l1+,gC+,lk+,s+,gA+,l1+,gC+,lk+,sp+,f2+
P h1#1#chr1 f1+,t+,s+,X+,gA+,l1+,gC+,lk+,sp+,f2+
P h2#1#chr1 f1+,t+,s+,X+,gA+,l1+,gC+,lk+,sp+,f2+
P h3#1#chr1 f1+,t+,gA+,l1+,gC+,lk+,s+,gA+,l1+,gC+,lk+,sp+,f2+
P h4#1#chr1 f1+,t+,gA+,l1+,gC+,lk+,s+,gA+,l1+,gC+,lk+,s+,gA+,l1+,gC+,lk+,sp+,f2+
P h5#1#chr1 f1+,t+,gA+,l1+,gC+,lk+,s+,X+,gA+,l1+,gC+,lk+,sp+,f2+
B f1 f2 t;X;s;gA;l1;gC;lk;sp
SPEC
check "a copy loss joined to an insertion across a shared node: one insertion, where h5's is (records)" \
  "$(th splitj "$OUT/th_splitj.spec" --cn)" "CNV/MODULE:011010 INS1000:011001 "
# A copy replaced by a diverged version of itself, 600 bp shorter, is no copy change beside an insertion: it is
# not split into a deletion of the copy and an insertion of the diverged one.
cat > "$OUT/th_splitd.spec" <<'SPEC'
N f1 2000 1
N t 10 9
N gA 1200 2
N l1 300 4
N gC 1200 8
N lk 300 6
N dA @gA:12:3
N dC @gC:12:5
N sp 200 11
N f2 2000 7
P ref#0#chr1:0-9010 f1+,t+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,sp+,f2+
P h1#1#chr1 f1+,t+,gA+,l1+,gC+,lk+,dA+,dC+,sp+,f2+
P h2#1#chr1 f1+,t+,gA+,l1+,gC+,lk+,dA+,dC+,sp+,f2+
P h3#1#chr1 f1+,t+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,sp+,f2+
P h4#1#chr1 f1+,t+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,sp+,f2+
B f1 f2 t;gA;l1;gC;lk;dA;dC;sp
SPEC
SPLITD="$(th splitd "$OUT/th_splitd.spec" --cn)"
check "a diverged copy is not split into a copy deletion and an insertion (records)" \
  "$(printf '%s' "$SPLITD" | grep -q 'INS2400' && echo split || echo kept)" "kept"
# A copy loss beside a larger novel insertion: the copy goes to the CNV, the insertion keeps its full size.
cat > "$OUT/th_splitb.spec" <<'SPEC'
N f1 2000 1
N t 10 9
N X 4000 10
N gA 1200 2
N l1 300 4
N gC 1200 8
N lk 300 6
N sp 200 11
N f2 2000 7
P ref#0#chr1:0-9010 f1+,t+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,sp+,f2+
P h1#1#chr1 f1+,t+,X+,gA+,l1+,gC+,lk+,sp+,f2+
P h2#1#chr1 f1+,t+,X+,gA+,l1+,gC+,lk+,sp+,f2+
P h3#1#chr1 f1+,t+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,sp+,f2+
P h4#1#chr1 f1+,t+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,sp+,f2+
P h5#1#chr1 f1+,t+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,sp+,f2+
B f1 f2 t;X;gA;l1;gC;lk;sp
SPEC
check "a copy loss beside a larger novel insertion: the whole insertion stays (records)" \
  "$(th splitb "$OUT/th_splitb.spec" --cn)" "CNV/MODULE:011010 INS4000:011000 "
# Hybrid copies: a partial copy after a hybrid one is a copy of its own. hpa carries GENEA's last two exons
# after its hybrid copy; read into the hybrid, they made it a GENEA copy again.
cat > "$OUT/hybp_gen.py" <<'PYEOF'
import random, sys
out = sys.argv[1]; rnd = random.Random(11)
rs = lambda n: ''.join(rnd.choice('ACGT') for _ in range(n))
mut = lambda s: ''.join({'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[c] if i % 24 == 5 else c for i, c in enumerate(s))
segs = {'f1': rs(800), 'f2': rs(800), 'lk': rs(400)}
for i in range(6):
    segs['ea%d' % i] = rs(150); segs['eb%d' % i] = mut(segs['ea%d' % i])
    if i < 5: segs['in%d' % i] = rs(250)
def copy(word):   # word: per exon 'a' or 'b' or '-' (absent)
    st = []
    for i, c in enumerate(word):
        if c != '-': st.append('e%s%d+' % (c, i))
        if i < 5 and (c != '-' or (i + 1 < len(word) and word[i+1] != '-')) and any(x != '-' for x in word[:i+1]): st.append('in%d+' % i)
    return st + ['lk+']
haps = {'ref': ['aaaaaa', 'bbbbbb'], 'h1': ['aaaaaa', 'bbbbbb'], 'h2': ['aaaaaa', 'aaaaaa', 'bbbbbb'],
        'h3': ['aaaaaa', 'bbbbbb', 'bbbbbb'], 'hhyb': ['aaaaaa', 'aaaaab', 'bbbbbb'],
        'hconv': ['aaaaaa', 'aabaaa', 'bbbbbb'], 'hpa': ['aaaaaa', 'aaaabb', '----aa', 'bbbbbb'], 'hpb': ['aaaaaa', 'aaaabb', 'bbbbbb']}
paths = {n: ['f1+'] + [s for w in ws for s in copy(w)] + ['f2+'] for n, ws in haps.items()}
edges = set((x, y) for st in paths.values() for x, y in zip(st, st[1:]))
seq = lambda st: ''.join(segs[s[:-1]] for s in st)
reflen = len(seq(paths['ref']))
with open(out + '.gfa', 'w') as f:
    f.write('H\tVN:Z:1.0\n')
    for k, v in segs.items(): f.write('S\t%s\t%s\n' % (k, v))
    for x, y in sorted(edges): f.write('L\t%s\t+\t%s\t+\t0M\n' % (x[:-1], y[:-1]))
    for n, st in paths.items():
        f.write('P\t%s\t%s\t*\n' % ('ref#0#chr9:1-%d' % reflen if n == 'ref' else n + '#1#chr9', ','.join(st)))
with open(out + '.csv', 'w') as f:
    f.write('bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n')
    f.write('1,f1,+,f2,+,%d,0,%d,"%s"\n' % (len(paths), reflen, ';'.join(sorted(k for k in segs if k not in ('f1', 'f2')))))
with open(out + '.gtf', 'w') as f:
    at = len(segs['f1']) + 1
    for name, c in (('GENEA', 'a'), ('GENEB', 'b')):
        lo = at
        for i in range(6):
            hi = at + 149
            f.write('chr9\tt\tCDS\t%d\t%d\t.\t+\t.\tgene_id "%s"; gene_name "%s";\n' % (at, hi, name, name))
            at = hi + 1 + (250 if i < 5 else 0)
        f.write('chr9\tt\tgene\t%d\t%d\t.\t+\t.\tgene_id "%s"; gene_name "%s";\n' % (lo, at - 1, name, name))
        at += len(segs['lk'])
PYEOF
python3 "$OUT/hybp_gen.py" "$OUT/hybp"
"$BIN" call -i "$OUT/hybp.gfa" -c "$OUT/hybp.csv" -r 'ref#0#chr9' -o "$OUT/hybp" --gtf "$OUT/hybp.gtf" --cn -q >/dev/null 2>&1
hybp_of() { awk -F'\t' -v h="$1" -v g="$2" 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "hybrid_copies") c = i }
  $3 ~ ("^" h "#") && $4 == g { print $c; exit }' "$OUT/hybp.dup_gene_cn.tsv"; }
check "a partial copy after a hybrid copy does not hide it (hybrids hpa GENEA, GENEB; hpb GENEA, GENEB)" \
  "$(hybp_of hpa GENEA) $(hybp_of hpa GENEB) $(hybp_of hpb GENEA) $(hybp_of hpb GENEB)" "1 1 1 1"
# Private k-mers are absent from every gene of the site, also one beside the copies: GENEC, beside them, shares
# GENEA's last three exons, so hfalse's copy (GENEB's first three exons, then GENEC's) is no hybrid of the pair,
# and with those exons out the pair has too few left to be assessed.
cat > "$OUT/hybc_gen.py" <<'PYEOF'
import random, sys

out = sys.argv[1]
rnd = random.Random(17)
rs = lambda n: ''.join(rnd.choice('ACGT') for _ in range(n))
mut = lambda s: ''.join({'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[c] if i % 24 == 5 else c
                        for i, c in enumerate(s))
segs = {'f1': rs(800), 'f2': rs(800), 'lk': rs(400)}
for i in range(6):
    segs['ea%d' % i] = rs(150)
    segs['eb%d' % i] = mut(segs['ea%d' % i])
    if i < 5:
        segs['in%d' % i] = rs(250)
for i in range(3, 6):
    segs['ec%d' % i] = segs['ea%d' % i]
    if i < 5:
        segs['ci%d' % i] = rs(250)

def full(letter):
    st = []
    for i in range(6):
        st.append('e%s%d+' % (letter, i))
        if i < 5:
            st.append('in%d+' % i)
    return st + ['lk+']

def ctail():
    st = []
    for i in range(3, 6):
        st.append('ec%d+' % i)
        if i < 5:
            st.append('ci%d+' % i)
    return st

def bhead():
    st = []
    for i in range(3):
        st.append('eb%d+' % i)
        st.append('in%d+' % i)
    return st

paths = {'ref': ['f1+'] + full('a') + full('b') + ctail() + ['f2+']}
for i in range(10):
    paths['h%d' % i] = ['f1+'] + full('a') + full('b') + ['f2+']
paths['hfalse'] = ['f1+'] + full('a') + bhead() + ctail() + ['lk+', 'f2+']
edges = set((x, y) for st in paths.values() for x, y in zip(st, st[1:]))

with open(out + '.gfa', 'w') as f:
    f.write('H\tVN:Z:1.0\n')
    for name, seq in segs.items():
        f.write('S\t%s\t%s\n' % (name, seq))
    for x, y in sorted(edges):
        f.write('L\t%s\t+\t%s\t+\t0M\n' % (x[:-1], y[:-1]))
    reflen = sum(len(segs[x[:-1]]) for x in paths['ref'])
    for name, st in paths.items():
        pname = 'ref#0#chr9:1-%d' % reflen if name == 'ref' else name + '#1#chr9'
        f.write('P\t%s\t%s\t*\n' % (pname, ','.join(st)))

inside = sorted(k for k in segs if k not in ('f1', 'f2'))
with open(out + '.csv', 'w') as f:
    f.write('bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n')
    f.write('1,f1,+,f2,+,%d,0,10000,"%s"\n' % (len(paths), ';'.join(inside)))

# Coordinates on the reference path.
at = 1
coords = {}
for token in paths['ref']:
    name = token[:-1]
    coords.setdefault(name, []).append((at, at + len(segs[name]) - 1))
    at += len(segs[name])
with open(out + '.gtf', 'w') as f:
    for gene, letter, occurrence in (('GENEA', 'a', 0), ('GENEB', 'b', 0)):
        ivs = []
        for i in range(6):
            lo, hi = coords['e%s%d' % (letter, i)][occurrence]
            ivs.append((lo, hi))
            f.write('chr9\tt\tCDS\t%d\t%d\t.\t+\t.\tgene_id "%s"; gene_name "%s";\n' %
                    (lo, hi, gene, gene))
        f.write('chr9\tt\tgene\t%d\t%d\t.\t+\t.\tgene_id "%s"; gene_name "%s";\n' %
                (ivs[0][0], ivs[-1][1], gene, gene))
    ivs = []
    for i in range(3, 6):
        lo, hi = coords['ec%d' % i][0]
        ivs.append((lo, hi))
        f.write('chr9\tt\tCDS\t%d\t%d\t.\t+\t.\tgene_id "GENEC"; gene_name "GENEC";\n' % (lo, hi))
    f.write('chr9\tt\tgene\t%d\t%d\t.\t+\t.\tgene_id "GENEC"; gene_name "GENEC";\n' %
            (ivs[0][0], ivs[-1][1]))
PYEOF
python3 "$OUT/hybc_gen.py" "$OUT/hybc"
"$BIN" call -i "$OUT/hybc.gfa" -c "$OUT/hybc.csv" -r "$(awk -F'\t' '$1 == "P" && $2 ~ /^ref#/ { print $2 }' "$OUT/hybc.gfa")" \
  -o "$OUT/hybc" --gtf "$OUT/hybc.gtf" --cn -q >/dev/null 2>&1
check "a gene beside the copies keeps its k-mers out of the pair's (hfalse GENEA, GENEB hybrid_copies)" \
  "$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "hybrid_copies") c = i }
      $3 ~ /^hfalse#/ && ($4 == "GENEA" || $4 == "GENEB") { printf "%s ", $c }' "$OUT/hybc.dup_gene_cn.tsv")" ". . "
# A module record under the carrier filters still gives its genes' copy numbers, whether it fell under them when
# made or after giving up a carrier: rows with variant_id '.'.
cat > "$OUT/th_ncf.spec" <<'SPEC'
N f1 2000 1
N t 10 9
N U 5000 13
N gA 1200 2
N l1 300 4
N gC 1200 8
N lk 300 6
N sp 200 11
N f2 2000 7
P ref#0#chr1:0-14210 f1+,t+,U+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,sp+,f2+
P h1#1#chr1 f1+,t+,U+,gA+,l1+,gC+,lk+,sp+,f2+
P h2#1#chr1 f1+,t+,U+,gA+,l1+,gC+,lk+,sp+,f2+
P h4#1#chr1 f1+,t+,U+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,sp+,f2+
B f1 f2 t;gA;l1;gC;lk;U;sp
SPEC
python3 -I "$TOOLS/spec_gfa.py" "$OUT/th_ncf.spec" "$OUT/th_ncf"
printf 'chr1\tt\tgene\t7011\t8210\t.\t+\t.\tgene_id "GA"; gene_name "GA";\nchr1\tt\tCDS\t7011\t8210\t.\t+\t.\tgene_id "GA"; gene_name "GA";\n' > "$OUT/th_ncf/g.gtf"
"$BIN" call -i "$OUT/th_ncf/g.gfa" -c "$OUT/th_ncf/g.bubbles.csv" -r 'ref#0#chr1:0-14210' -o "$OUT/th_ncf/v" --cn \
  --gtf "$OUT/th_ncf/g.gtf" --min-haplotypes 3 --no-per-bubble-vcf -q >/dev/null 2>&1
check "a module record under --min-haplotypes still writes its gene rows (variant_id sample:cn)" \
  "$(awk -F'\t' 'NR > 1 { split($3, s, "#"); printf "%s %s:%s ", $2, s[1], $5 }' "$OUT/th_ncf/v.dup_gene_cn.tsv")" \
  ". ref:2 . h1:1 . h2:1 . h4:2 "

# A gene in the copies has no more copies than the module. h2 carries a stray GA copy (gA2) outside the
# copies: its GA dosage nears 2 while the module total is 1, so GA reads 1, unreliable.
cat > "$OUT/gs_cap.spec" <<'SPEC'
N f1 4000 1
N sp 200 22
N gA 1200 2
N l1 300 4
N gC 1200 8
N lk 300 6
N f2 4000 7
N f3 3000 9
N gA2 @gA:2:5
P ref f1+,sp+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,f2+,f3+
P h1 f1+,sp+,gA+,l1+,gC+,lk+,f2+,f3+
P h2 f1+,gA2+,sp+,gA+,l1+,gC+,lk+,f2+,f3+
P h3 f1+,sp+,gA+,l1+,gC+,lk+,gA+,l1+,gC+,lk+,f2+,f3+
B f1 f2 sp;gA;l1;gC;lk
SPEC
check "a gene's copies are capped at the module's: h2 GA cn:reliable" \
  "$(gs cap "$OUT/gs_cap.spec" GA:gA GC:gC | tr ' ' '\n' | grep ':h2:GA:' | cut -d: -f4,5)" "1:0"
# Hybrids common in the panel, or beside a partial copy, are still read: their crossover exons stay in.
sed "s/'hconv': \['aaaaaa', 'aabaaa', 'bbbbbb'\]}/'hconv': ['aaaaaa', 'aabaaa', 'bbbbbb'], 'hyb2': ['aaaaaa', 'aaaaab', 'bbbbbb'], 'hyb3': ['aaaaaa', 'aaaaab', 'bbbbbb'], 'hpa': ['aaaaaa', 'aaa', 'bbbbbb']}/" \
    "$OUT/hyb_gen.py" > "$OUT/hybm_gen.py"
python3 "$OUT/hybm_gen.py" "$OUT/hybm"
"$BIN" call -i "$OUT/hybm.gfa" -c "$OUT/hybm.csv" -r 'ref#0#chr9' -o "$OUT/hybm" --gtf "$OUT/hybm.gtf" --cn -q >/dev/null 2>&1
hybm_of() { awk -F'\t' -v h="$1" -v g="$2" 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "hybrid_copies") c = i }
  $3 ~ ("^" h "#") && $4 == g { print $c; exit }' "$OUT/hybm.dup_gene_cn.tsv"; }
check "three hybrid haplotypes and a partial copy: hybrids of hhyb, hyb2, hyb3 (GENEA), hconv, hpa" \
  "$(hybm_of hhyb GENEA) $(hybm_of hyb2 GENEA) $(hybm_of hyb3 GENEA) $(hybm_of hconv GENEA) $(hybm_of hpa GENEA)" "1 1 1 0 0"
# A gene in two separate hybrid copies, one with each of two paralogs, has two; one copy seen by two pairs
# stays one.
cat > "$OUT/hyb3_gen.py" <<'PYEOF'
import random, sys
out = sys.argv[1]; rnd = random.Random(11)
rs = lambda n: ''.join(rnd.choice('ACGT') for _ in range(n))
sub = lambda s, m: ''.join(m[c] if i % 24 == 5 else c for i, c in enumerate(s))
segs = {'f1': rs(2000), 'f2': rs(2000), 'lk': rs(400)}
for i in range(6):
    segs['ea%d' % i] = rs(150)
    segs['eb%d' % i] = sub(segs['ea%d' % i], {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'})
    segs['ec%d' % i] = sub(segs['ea%d' % i], {'A': 'G', 'C': 'T', 'G': 'A', 'T': 'C'})
    if i < 5: segs['in%d' % i] = rs(250)
def copy(word):
    st = []
    for i, c in enumerate(word):
        st.append('e%s%d+' % (c, i))
        if i < 5: st.append('in%d+' % i)
    return st + ['lk+']
haps = {'ref': ['aaaaaa', 'bbbbbb', 'cccccc'], **{'c%d' % i: ['aaaaaa', 'bbbbbb', 'cccccc'] for i in range(12)},
        'target': ['aaaaaa', 'aaaaab', 'aaaaac', 'bbbbbb', 'cccccc']}
paths = {n: ['f1+'] + [s for w in ws for s in copy(w)] + ['f2+'] for n, ws in haps.items()}
edges = set((x, y) for st in paths.values() for x, y in zip(st, st[1:]))
reflen = len(''.join(segs[s[:-1]] for s in paths['ref']))
with open(out + '.gfa', 'w') as f:
    f.write('H\tVN:Z:1.0\n')
    for k, v in segs.items(): f.write('S\t%s\t%s\n' % (k, v))
    for x, y in sorted(edges): f.write('L\t%s\t+\t%s\t+\t0M\n' % (x[:-1], y[:-1]))
    for n, st in paths.items():
        f.write('P\t%s\t%s\t*\n' % ('ref#0#chr9:1-%d' % reflen if n == 'ref' else n + '#1#chr9', ','.join(st)))
with open(out + '.csv', 'w') as f:
    f.write('bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n')
    f.write('1,f1,+,f2,+,%d,0,%d,"%s"\n' % (len(paths), reflen, ';'.join(sorted(k for k in segs if k not in ('f1', 'f2')))))
with open(out + '.gtf', 'w') as f:
    at = len(segs['f1']) + 1
    for name in ('GENEA', 'GENEB', 'GENEC'):
        lo = at
        for i in range(6):
            hi = at + 149
            f.write('chr9\tt\tCDS\t%d\t%d\t.\t+\t.\tgene_id "%s"; gene_name "%s";\n' % (at, hi, name, name))
            at = hi + 1 + (250 if i < 5 else 0)
        f.write('chr9\tt\tgene\t%d\t%d\t.\t+\t.\tgene_id "%s"; gene_name "%s";\n' % (lo, at - 1, name, name))
        at += len(segs['lk'])
PYEOF
python3 "$OUT/hyb3_gen.py" "$OUT/hyb3"
"$BIN" call -i "$OUT/hyb3.gfa" -c "$OUT/hyb3.csv" -r 'ref#0#chr9' -o "$OUT/hyb3" --gtf "$OUT/hyb3.gtf" --cn -q >/dev/null 2>&1
check "two hybrid copies of GENEA (with GENEB, with GENEC): hybrids of GENEA GENEB GENEC" \
  "$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "hybrid_copies") c = i }
     $3 ~ /^target#/ { printf "%s ", $c }' "$OUT/hyb3.dup_gene_cn.tsv")" "2 1 1 "

# A W line may give a start and no end ('*'): the start still places the records.
for wend in 1260 '*'; do
  d="$OUT/wend_${wend/\*/star}"; mkdir -p "$d"
  { printf 'H\tVN:Z:1.0\nS\tf1\t%s\nS\tx\t%s\nS\tf2\t%s\n' "$(printf 'ACGTTGCA%.0s' $(seq 1 13))" \
      "$(printf 'GGATCCAT%.0s' $(seq 1 8))" "$(printf 'TTGACCAG%.0s' $(seq 1 13))"
    printf 'L\tf1\t+\tx\t+\t0M\nL\tx\t+\tf2\t+\t0M\nL\tf1\t+\tf2\t+\t0M\n'
    printf 'W\tref\t0\tchr9\t1000\t%s\t>f1>x>f2\nW\th1\t1\tchr9\t0\t*\t>f1>f2\n' "$wend"; } > "$d/g.gfa"
  printf 'bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,f1,+,f2,+,2,0,64,x\n' > "$d/g.csv"
  "$BIN" call -i "$d/g.gfa" -c "$d/g.csv" -r ref -o "$d/v" -q >/dev/null 2>&1
done
check "a W start with an unknown end places the record as with a known end (POS)" \
  "$(grep -v '^#' "$OUT/wend_star/v.region.vcf" | cut -f2 | head -1)" "$(grep -v '^#' "$OUT/wend_1260/v.region.vcf" | cut -f2 | head -1)"

# ---------------------------------------------------------------------------------------------
# A folded unit the reference does not carry (REF_CN=0) has no reference orientation. Its direct
# orientation is the one most haplotypes' copies take, whichever way the node is stored: u- in one
# graph is u+ in the same graph with u stored reverse-complemented. h1 and h2 carry direct copies,
# h3 one direct and one inverted copy, so CNINV is 0 0 1 and INVDUP is set either way.
# ---------------------------------------------------------------------------------------------
U=$(blk 60 41); UR=$(printf '%s' "$U" | rev | tr ACGT TGCA)
for st in minus plus; do
  d="$OUT/orient_$st"; mkdir -p "$d"
  if [ "$st" = minus ]; then us=$U; a=-; b=+; else us=$UR; a=+; b=-; fi
  { printf 'H\tVN:Z:1.0\nS\tf1\t%s\nS\tu\t%s\nS\tf2\t%s\n' "$(blk 100 43)" "$us" "$(blk 100 47)"
    printf 'L\tf1\t+\tf2\t+\t0M\nL\tf1\t+\tu\t%s\t0M\nL\tu\t%s\tu\t%s\t0M\nL\tu\t%s\tu\t%s\t0M\n' "$a" "$a" "$a" "$a" "$b"
    printf 'L\tu\t%s\tf2\t+\t0M\nL\tu\t%s\tf2\t+\t0M\n' "$a" "$b"
    printf 'P\tref#0#chr1:1-200\tf1+,f2+\t*\n'
    printf 'P\th1#1#c\tf1+,u%s,u%s,f2+\t*\nP\th2#1#c\tf1+,u%s,f2+\t*\nP\th3#1#c\tf1+,u%s,u%s,f2+\t*\n' "$a" "$a" "$a" "$a" "$b"
  } > "$d/g.gfa"
  printf 'bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,f1,+,f2,+,4,0,120,u\n' > "$d/g.csv"
  "$BIN" call -i "$d/g.gfa" -c "$d/g.csv" -r 'ref#0#chr1:1-200' -o "$d/v" -q >/dev/null 2>&1
  check "REF_CN=0 unit walked $st: CNINV of h1 h2 h3, INVDUP" \
    "$(awk -F'\t' '/^#CHROM/ { for (i = 10; i <= NF; i++) col[$i] = i }
         !/^#/ && $8 ~ /SVTYPE=CNV/ { n = split($9, k, ":"); for (j = 1; j <= n; j++) if (k[j] == "CNINV") c = j
           for (h = 1; h <= 3; h++) { split($(col["h" h "#1#c"]), v, ":"); printf "%s ", v[c] }
           printf "%s", ($8 ~ /(^|;)INVDUP(;|$)/ ? "INVDUP" : "-") }' "$d/v.region.vcf")" "0 0 1 INVDUP"
done

# ---------------------------------------------------------------------------------------------
# A substitution merged with deletions: every member counts by its net change, as SVLEN does. Node b
# (100 bp) is deleted by hD haplotypes and replaced by 15 unrelated bases (net -85) by hS ones; the
# record's range is -100..-85 whichever member represents it.
# ---------------------------------------------------------------------------------------------
B100=$(blk 100 53 | tr GT AC); S15=$(blk 15 89 | tr AC GT)   # no base in common: one change, no chance matches
for mix in "hD1 hD2 hS1" "hS1 hS2 hD1"; do
  d="$OUT/submix_${mix// /_}"; mkdir -p "$d"
  { printf 'H\tVN:Z:1.0\nS\tf1\t%s\nS\tb\t%s\nS\ts\t%s\nS\tf2\t%s\n' "$(blk 100 53)" "$B100" "$S15" "$(blk 100 59)"
    printf 'L\tf1\t+\tb\t+\t0M\nL\tb\t+\tf2\t+\t0M\nL\tf1\t+\tf2\t+\t0M\nL\tf1\t+\ts\t+\t0M\nL\ts\t+\tf2\t+\t0M\n'
    printf 'P\tref#0#chr1:1-300\tf1+,b+,f2+\t*\nP\thR#1#c\tf1+,b+,f2+\t*\n'
    for h in $mix; do
      case $h in hD*) printf 'P\t%s#1#c\tf1+,f2+\t*\n' "$h";; hS*) printf 'P\t%s#1#c\tf1+,s+,f2+\t*\n' "$h";; esac
    done; } > "$d/g.gfa"
  printf 'bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,f1,+,f2,+,5,0,100,b;s\n' > "$d/g.csv"
  "$BIN" call -i "$d/g.gfa" -c "$d/g.csv" -r 'ref#0#chr1:1-300' -o "$d/v" -q >/dev/null 2>&1
  check "deletion and substitution merged ($mix): records, SVLEN_RANGE" \
    "$(awk -F'\t' '!/^#/ { n++; r = "-"; m = "-"; k = split($8, a, ";")
         for (i = 1; i <= k; i++) if (a[i] ~ /^SVLEN_RANGE=/) r = substr(a[i], 13) }
         END { printf "%d %s", n, r }' "$d/v.region.vcf")" "1 -100,-85"
done

# ---------------------------------------------------------------------------------------------
# Two replacements of node b (100 bp) merged into one record: hS1 puts 170 bases there, hS2 the same and
# 10 more (net +70 and +80). Members' sizes in SVLEN_RANGE are their net changes; a replacement's larger
# side (170, 180) is no size SVLEN gives.
# ---------------------------------------------------------------------------------------------
d="$OUT/subsub"; mkdir -p "$d"
S170=$(blk 170 97 | tr AC GT); S180="$S170$(blk 10 101 | tr AC GT)"
{ printf 'H\tVN:Z:1.0\nS\tf1\t%s\nS\tb\t%s\nS\ts1\t%s\nS\ts2\t%s\nS\tf2\t%s\n' "$(blk 100 53)" "$B100" "$S170" "$S180" "$(blk 100 59)"
  printf 'L\tf1\t+\tb\t+\t0M\nL\tb\t+\tf2\t+\t0M\nL\tf1\t+\ts1\t+\t0M\nL\ts1\t+\tf2\t+\t0M\nL\tf1\t+\ts2\t+\t0M\nL\ts2\t+\tf2\t+\t0M\n'
  printf 'P\tref#0#chr1:1-300\tf1+,b+,f2+\t*\nP\thR#1#c\tf1+,b+,f2+\t*\n'
  printf 'P\thS1#1#c\tf1+,s1+,f2+\t*\nP\thS2#1#c\tf1+,s2+,f2+\t*\n'; } > "$d/g.gfa"
printf 'bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,f1,+,f2,+,4,100,180,b;s1;s2\n' > "$d/g.csv"
"$BIN" call -i "$d/g.gfa" -c "$d/g.csv" -r 'ref#0#chr1:1-300' -o "$d/v" -q >/dev/null 2>&1
check "two replacements merged: records, SVLEN_RANGE by net change" \
  "$(awk -F'\t' '!/^#/ { n++; r = "-"; k = split($8, a, ";")
       for (i = 1; i <= k; i++) if (a[i] ~ /^SVLEN_RANGE=/) r = substr(a[i], 13) }
       END { printf "%d %s", n, r }' "$d/v.region.vcf")" "1 70,80"

# ---------------------------------------------------------------------------------------------
# --min-alt-af is tested on the AF the record is written with. Per path: hBig deletes a, b and c, so it
# has no allele at hSmall's deletion of b (GT .), whose AF is 1/2; hBig's own is 1/3 (the reference path is
# out of AN). Per sample (--diploid): S1's deletion is 1 of 4 haplotypes.
# ---------------------------------------------------------------------------------------------
d="$OUT/maf"; mkdir -p "$d"
{ printf 'H\tVN:Z:1.0\nS\tf1\t%s\nS\ta\t%s\nS\tb\t%s\nS\tc\t%s\nS\tf2\t%s\n' \
    "$(blk 100 61)" "$(blk 100 67)" "$(blk 100 71)" "$(blk 100 73)" "$(blk 100 79)"
  printf 'L\tf1\t+\ta\t+\t0M\nL\ta\t+\tb\t+\t0M\nL\tb\t+\tc\t+\t0M\nL\tc\t+\tf2\t+\t0M\nL\tf1\t+\tf2\t+\t0M\nL\ta\t+\tc\t+\t0M\n'
  printf 'P\tref#0#chr1:1-500\tf1+,a+,b+,c+,f2+\t*\nP\thR#1#c\tf1+,a+,b+,c+,f2+\t*\n'
  printf 'P\thBig#1#c\tf1+,f2+\t*\nP\thSmall#1#c\tf1+,a+,c+,f2+\t*\n'; } > "$d/g.gfa"
printf 'bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,f1,+,f2,+,4,0,300,a;b;c\n' > "$d/g.csv"
svlens() { awk -F'\t' '!/^#/ { k = split($8, a, ";"); for (i = 1; i <= k; i++) if (a[i] ~ /^(SVLEN|AF)=/) printf "%s ", a[i] }' "$1"; }
# A run that drops every record still succeeds and writes a VCF of its header, so each run's exit status
# and file are checked: an empty result alone would also pass for a run that failed. The log counts only
# the '.' genotypes of records written.
maf_runs=""
for af in 0 0.4 0.51; do
  "$BIN" call -i "$d/g.gfa" -c "$d/g.csv" -r 'ref#0#chr1:1-500' -o "$d/af$af" --min-alt-af "$af" > "$d/af$af.log" 2>&1
  maf_runs="$maf_runs$? $([ -f "$d/af$af.region.vcf" ] && echo vcf || echo none) "
done
check "--min-alt-af 0 / 0.4 / 0.51: exit status and VCF" "$maf_runs" "0 vcf 0 vcf 0 vcf "
check "--min-alt-af 0 / 0.4 / 0.51 keep the records written with AF at or above it" \
  "$(svlens "$d/af0.region.vcf")| $(svlens "$d/af0.4.region.vcf")| $(svlens "$d/af0.51.region.vcf")" \
  "SVLEN=-300 AF=0.3333 SVLEN=-100 AF=0.5 | SVLEN=-100 AF=0.5 | "
check "--min-alt-af: the log's 'GT . within a larger deletion' counts written records only (0 / 0.4 / 0.51)" \
  "$(for af in 0 0.4 0.51; do printf '%s ' "$(grep -c 'GT . within a larger deletion' "$d/af$af.log")"; done)" "1 1 0 "
dd="$OUT/maf_diploid"; mkdir -p "$dd"
{ grep -v '^P' "$d/g.gfa"
  printf 'P\tref#0#chr1:1-500\tf1+,a+,b+,c+,f2+\t*\nP\tS1#1#c\tf1+,b+,c+,f2+\t*\nP\tS1#2#c\tf1+,a+,b+,c+,f2+\t*\n'
  printf 'P\tS2#1#c\tf1+,a+,b+,c+,f2+\t*\nP\tS2#2#c\tf1+,a+,b+,c+,f2+\t*\n'; } > "$dd/g.gfa"
printf 'L\tf1\t+\tb\t+\t0M\n' >> "$dd/g.gfa"
dip_runs=""
for af in 0.25 0.26; do
  "$BIN" call -i "$dd/g.gfa" -c "$d/g.csv" -r 'ref#0#chr1:1-500' -o "$dd/af$af" --diploid --min-alt-af "$af" -q >/dev/null 2>&1
  dip_runs="$dip_runs$? $([ -f "$dd/af$af.region.vcf" ] && echo vcf || echo none) "
done
check "--diploid --min-alt-af 0.25 / 0.26: exit status and VCF" "$dip_runs" "0 vcf 0 vcf "
check "--diploid --min-alt-af 0.25 / 0.26 against the written AF of 1 in 4" \
  "$(svlens "$dd/af0.25.region.vcf")| $(svlens "$dd/af0.26.region.vcf")" "SVLEN=-100 AF=0.25 | "
check "a record --min-alt-af drops is in neither the per-site VCF nor the node table (records in each)" \
  "$(grep -vc '^#' "$d"/af0.4.bubble_1.vcf) $(tail -n +2 "$d/af0.4.variant_nodes.tsv" | wc -l | tr -d ' ')" "1 1"

# ---------------------------------------------------------------------------------------------
# Site ids start at 1: 0 is what "no site" is written as, so a table holding one is refused.
# ---------------------------------------------------------------------------------------------
sed '2s/^1,/0,/' "$OUT/cb.bubbles.csv" > "$OUT/id0.bubbles.csv"
refuses "a bubbles CSV with bubble_id 0 is refused" "bubble_id 0" \
  "$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/id0.bubbles.csv" -r "$REF_NAME" -q

# A file named like a per-site VCF with an id too large for one is not call's: the run neither fails
# on it nor removes it.
d="$OUT/bigid"; mkdir -p "$d"; : > "$d/v.bubble_123456789012345678901234.vcf"
"$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" -o "$d/v" -q >/dev/null 2>&1; rc=$?
check "an over-long per-site VCF name: exit status, file kept" "$rc $(ls "$d"/v.bubble_123456789012345678901234.vcf 2>/dev/null | wc -l | tr -d ' ')" "0 1"

# ---------------------------------------------------------------------------------------------
# A block too large for the node comparison that base alignment then compares in full is not
# reported as skipped: here 8100 one-base nodes, the haplotype on nodes of its own with 100 bases
# deleted, which the aligned allele calls.
# ---------------------------------------------------------------------------------------------
d="$OUT/bigblock"; mkdir -p "$d"
awk -v L=8100 'BEGIN { b = "ACGT"; x = 9 * 7919 + 1
    printf "H\tVN:Z:1.0\nS\tS\t%s\nS\tK\t%s\n", "ACGTTGCAAG", "TTGACCAGTC"
    for (i = 0; i < L; i++) { x = (x * 16807) % 2147483647; s[i] = substr(b, int(x * 4 / 2147483647) + 1, 1)
      printf "S\tr%d\t%s\n", i, s[i]; if (i < 4000 || i >= 4100) printf "S\th%d\t%s\n", i, s[i] }
    r = "S+"; h = "S+"; pr = "S"; ph = "S"
    for (i = 0; i < L; i++) { r = r ",r" i "+"; printf "L\t%s\t+\tr%d\t+\t0M\n", pr, i; pr = "r" i
      if (i < 4000 || i >= 4100) { h = h ",h" i "+"; printf "L\t%s\t+\th%d\t+\t0M\n", ph, i; ph = "h" i } }
    printf "L\t%s\t+\tK\t+\t0M\nL\t%s\t+\tK\t+\t0M\n", pr, ph
    printf "P\tref#0#chr1:1-8120\t%s,K+\t*\nP\thD#1#c\t%s,K+\t*\nP\thR#1#c\t%s,K+\t*\n", r, h, r }' > "$d/g.gfa"
{ printf 'bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,S,+,K,+,3,8000,8100,"'
  awk -F'\t' '$1 == "S" && $2 != "S" && $2 != "K" { printf "%s%s", (n++ ? ";" : ""), $2 }' "$d/g.gfa"; printf '"\n'; } > "$d/g.csv"
"$BIN" call -i "$d/g.gfa" -c "$d/g.csv" -r 'ref#0#chr1:1-8120' -o "$d/v" > "$d/log" 2>&1
check "a block the aligned allele covers: DEL SVLEN and GT, 'too large to align' log lines" \
  "$(awk -F'\t' '/^#CHROM/ { for (i = 10; i <= NF; i++) if ($i == "hD#1#c") c = i }
       !/^#/ && $8 ~ /SVTYPE=DEL/ { k = split($8, a, ";"); for (i = 1; i <= k; i++) if (a[i] ~ /^SVLEN=/) printf "%s ", a[i]; split($c, g, ":"); printf "%s ", g[1] }' "$d/v.region.vcf")$(grep -c 'too large to align' "$d/log")" \
  "SVLEN=-100 1 0"

# ---------------------------------------------------------------------------------------------
# One row per gene and haplotype in dup_gene_cn.tsv when a site holds two folded repeats and no module:
# repeat ra carries gene GA, repeat rb gene GB, each written once, under the record whose copies hold it.
# ---------------------------------------------------------------------------------------------
d="$OUT/tworep"; mkdir -p "$d"
{ printf 'H\tVN:Z:1.0\nS\t1\t%s\nS\t2\t%s\nS\t3\t%s\nS\t4\t%s\nS\t5\t%s\n' \
    "$(blk 400 101)" "$(blk 1200 103)" "$(blk 300 107)" "$(blk 1100 109)" "$(blk 400 113)"
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\nL\t3\t+\t4\t+\t0M\nL\t4\t+\t4\t+\t0M\nL\t4\t+\t5\t+\t0M\nL\t1\t+\t5\t+\t0M\n'
  printf 'P\tgref#0#chr9:1-5700\t1+,2+,2+,3+,4+,4+,5+\t*\nP\th1#1#chr9\t1+,2+,2+,2+,3+,4+,4+,5+\t*\n'
  printf 'P\th2#1#chr9\t1+,2+,2+,3+,4+,4+,4+,5+\t*\nP\thD#1#chr9\t1+,5+\t*\n'; } > "$d/g.gfa"
printf 'bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n1,1,+,5,+,4,0,6100,2;3;4\n' > "$d/g.csv"
for g in "GA 401 1600" "GB 3101 4200"; do
  set -- $g
  printf 'chr9\tt\tgene\t%s\t%s\t.\t+\t.\tgene_id "%s"; gene_name "%s";\n' "$2" "$3" "$1" "$1"
  printf 'chr9\tt\tCDS\t%s\t%s\t.\t+\t0\tgene_id "%s"; gene_name "%s";\n' "$2" "$3" "$1" "$1"
done > "$d/g.gtf"
"$BIN" call -i "$d/g.gfa" -c "$d/g.csv" -r 'gref#0#chr9:1-5700' --gtf "$d/g.gtf" -o "$d/v" -q >/dev/null 2>&1
check "two repeats, one gene each: rows, and the records GA and GB are written under" \
  "$(tail -n +2 "$d/v.dup_gene_cn.tsv" | wc -l | tr -d ' ') $(tail -n +2 "$d/v.dup_gene_cn.tsv" | cut -f2,4 | sort -u | tr '\t\n' ': ')" \
  "8 bubble1_CNV_2:GA bubble1_CNV_4:GB "

# ---------------------------------------------------------------------------------------------
# Copy orientation where the reference has no copy of the unit u. cnv_inv <name> <spec> [call options]
# prints INVDUP (or -) and each column's CN:CNINV.
# ---------------------------------------------------------------------------------------------
cnv_inv() {
  local d="$OUT/oi_$1" spec="$2"; shift 2
  python3 -I "$TOOLS/spec_gfa.py" "$spec" "$d"
  "$BIN" call -i "$d/g.gfa" -c "$d/g.bubbles.csv" -r "$(awk '$1 == "P" { print $2; exit }' "$spec")" -o "$d/v" \
    --no-per-bubble-vcf -q "$@" >/dev/null 2>&1
  awk -F'\t' '!/^#/ && $8 ~ /SVTYPE=CNV/ { printf "%s", ($8 ~ /(^|;)INVDUP(;|$)/ ? "INVDUP" : "-")
      n = split($9, f, ":"); for (k = 1; k <= n; k++) { if (f[k] == "CN") c = k; if (f[k] == "CNINV") v = k }
      for (i = 10; i <= NF; i++) { split($i, x, ":"); printf " %s:%s", x[c], x[v] } }' "$d/v.region.vcf"; }
# h1 holds a copy and its reverse complement, a fold-back: an inverted duplication whichever orientation is
# called direct, even where the panel's copies are evenly split between the two orientations.
cat > "$OUT/oi_fold.spec" <<'SPEC'
N f1 100 301
N u 60 302
N f2 100 303
P ref#0#chr1:1-200 f1+,f2+
P h1#1#c f1+,u+,u-,f2+
P h2#1#c f1+,f2+
P h3#1#c f1+,f2+
B f1 f2 u
SPEC
check "no reference copy, a fold-back carrier: INVDUP, CN:CNINV" "$(cnv_inv fold "$OUT/oi_fold.spec")" "INVDUP 0:0 2:1 0:0 0:0"
# h4's one copy runs against the other haplotypes' and duplicates nothing: CNINV 1, no INVDUP.
cat > "$OUT/oi_minor.spec" <<'SPEC'
N f1 100 301
N u 60 302
N f2 100 303
P ref#0#chr1:1-200 f1+,f2+
P h1#1#c f1+,u+,f2+
P h2#1#c f1+,u+,f2+
P h3#1#c f1+,u+,u+,f2+
P h4#1#c f1+,u-,f2+
B f1 f2 u
SPEC
check "no reference copy, one copy against the panel's orientation: no INVDUP, CN:CNINV" "$(cnv_inv minor "$OUT/oi_minor.spec")" "- 0:0 1:0 1:0 2:0 1:1"
# Evenly split, each haplotype in one orientation: no copy reverse-complements another of its haplotype.
cat > "$OUT/oi_tie.spec" <<'SPEC'
N f1 100 301
N u 60 302
N f2 100 303
P ref#0#chr1:1-200 f1+,f2+
P h1#1#c f1+,u+,u+,f2+
P h2#1#c f1+,u-,u-,f2+
B f1 f2 u
SPEC
check "no reference copy, an even split across haplotypes: no INVDUP, CN:CNINV" "$(cnv_inv tie "$OUT/oi_tie.spec")" "- 0:0 2:0 2:0"
# A unit graph whose two anchors are alternative forward copies (V is U with one change): h1 holds two V
# copies. Stored on the other strand it is the same haplotype, and its copies are not inverted.
for o in fwd rev; do
  { printf 'N F 500 311\nN S 100 312\nN U 100 313\nN V @U:1:9\nN K 100 315\nN G 500 316\n'
    printf 'P ref#0#chr1:1-1300 F+,S+,U+,K+,G+\n'
    [ $o = fwd ] && printf 'P h1#1#chr1 F+,S+,V+,V+,K+,G+\n' || printf 'P h1#1#chr1 G-,K-,V-,V-,S-,F-\n'
    printf 'B S K U;V\n'; } > "$OUT/oi_ug$o.spec"
  printf 'input_bubble_id\tunit_bp\tanchor_pos\tanchor_nodes\tn_nodes\tnodes\n1\t100\t0\tU;V\t2\tU;V\n' > "$OUT/oi_ug$o.panphorte.unit_graph.tsv"
done
check "a unit graph's copies, haplotype stored on either strand: INVDUP, CN:CNINV (forward | reverse)" \
  "$(cnv_inv ugfwd "$OUT/oi_ugfwd.spec" --panphorte-prefix "$OUT/oi_ugfwd")|$(cnv_inv ugrev "$OUT/oi_ugrev.spec" --panphorte-prefix "$OUT/oi_ugrev")" \
  "- 1:0 2:0|- 1:0 2:0"
# h1 keeps U and gains V inverted: an inverted duplication on either strand. Vr is V stored on the other
# strand; walked Vr- its copy is direct, walked Vr+ inverted. The sequence, not the walks, says how each
# alternative anchor is stored, so a haplotype's only copy on it can be inverted.
ugv=$(awk -F'\t' '$1 == "S" && $2 == "V" { print $3 }' "$OUT/oi_ugfwd/g.gfa")
ugvr=$(printf '%s' "$ugv" | rev | tr ACGT TGCA)
ug_inv=""
for c in "fwdinv V F+,S+,U+,V-,K+,G+" "revinv V G-,K-,V+,U-,S-,F-" "flipdir Vr F+,S+,U+,Vr-,K+,G+" "flipdirrev Vr G-,K-,Vr+,U-,S-,F-" \
         "flipinv Vr F+,S+,U+,Vr+,K+,G+" "flipinvrev Vr G-,K-,Vr-,U-,S-,F-"; do
  read -r name alt walk <<< "$c"
  { printf 'N F 500 311\nN S 100 312\nN U 100 313\n'
    [ "$alt" = V ] && printf 'N V =%s\n' "$ugv" || printf 'N Vr =%s\n' "$ugvr"
    printf 'N K 100 315\nN G 500 316\nP ref#0#chr1:1-1300 F+,S+,U+,K+,G+\nP h1#1#chr1 %s\nB S K U;%s\n' "$walk" "$alt"; } \
    > "$OUT/oi_ug$name.spec"
  printf 'input_bubble_id\tunit_bp\tanchor_pos\tanchor_nodes\tn_nodes\tnodes\n1\t100\t0\tU;%s\t2\tU;%s\n' "$alt" "$alt" \
    > "$OUT/oi_ug$name.panphorte.unit_graph.tsv"
  ug_inv="$ug_inv$(cnv_inv "ug$name" "$OUT/oi_ug$name.spec" --panphorte-prefix "$OUT/oi_ug$name")|"
done
check "U and an inverted V, either strand; V stored flipped, walked direct then inverted (CN:CNINV)" "$ug_inv" \
  "INVDUP 1:0 2:1|INVDUP 1:0 2:1|- 1:0 2:0|- 1:0 2:0|INVDUP 1:0 2:1|INVDUP 1:0 2:1|"
# Short alternative anchors (a 1 bp C against G, A against T, a 3 bp GAT against C) match the other strand by
# chance; every copy here is direct, so no haplotype has an inverted copy.
ug_short=""
for c in "cg C G" "at A T" "gat C GAT"; do
  read -r name u v <<< "$c"
  { printf 'N F 500 311\nN S 100 312\nN A 49 320\nN U =%s\nN V =%s\nN B 50 321\nN K 100 315\nN G 500 316\n' "$u" "$v"
    printf 'P ref#0#chr1:1-1300 F+,S+,A+,U+,B+,K+,G+\nP h1#1#chr1 F+,S+,A+,U+,B+,A+,V+,B+,K+,G+\n'
    printf 'P h2#1#chr1 F+,S+,A+,V+,B+,K+,G+\nP h3#1#chr1 F+,S+,A+,U+,B+,K+,G+\nB S K A;U;V;B\n'; } > "$OUT/oi_ugs$name.spec"
  printf 'input_bubble_id\tunit_bp\tanchor_pos\tanchor_nodes\tn_nodes\tnodes\n1\t100\t50\tU;V\t4\tA;U;V;B\n' \
    > "$OUT/oi_ugs$name.panphorte.unit_graph.tsv"
  ug_short="$ug_short$(cnv_inv "ugs$name" "$OUT/oi_ugs$name.spec" --panphorte-prefix "$OUT/oi_ugs$name")|"
done
check "1-3 bp alternative anchors, every copy direct: no INVDUP (CN:CNINV; C/G, A/T, C/GAT)" "$ug_short" \
  "- 1:0 2:0 1:0 1:0|- 1:0 2:0 1:0 1:0|- 1:0 2:0 1:0 1:0|"
# Two one-node rows of one motif joined into one repeat: V is U rotated, and Vr that stored on the other
# strand. h1 keeps U and gains V inverted (either strand of h1), V direct, Vr walked direct, Vr inverted.
ju=$(blk 100 361); jv="${ju:17}${ju:0:17}"; jvr=$(printf '%s' "$jv" | rev | tr ACGT TGCA)
ug_join=""
for c in "inv $jv S+,U+,V-,K+" "invrev $jv K-,V+,U-,S-" "dir $jv S+,U+,V+,K+" "rdir $jvr S+,U+,V-,K+" "rinv $jvr S+,U+,V+,K+"; do
  read -r name vs walk <<< "$c"
  printf 'N S 100 312\nN U =%s\nN V =%s\nN K 100 315\nP ref#0#chr1:1-300 S+,U+,K+\nP h1#1#chr1 %s\nB S K U;V\n' \
    "$ju" "$vs" "$walk" > "$OUT/oi_ugj$name.spec"
  { printf 'input_bubble_id\tunit_bp\tcanonical_motif\tunit\tanchor_pos\tanchor_nodes\tn_nodes\tnodes\n'
    printf '1\t100\t%s\t%s\t1\tU\t1\tU\n1\t100\t%s\t%s\t1\tV\t1\tV\n' "$ju" "$ju" "$ju" "$vs"; } \
    > "$OUT/oi_ugj$name.panphorte.unit_graph.tsv"
  ug_join="$ug_join$(cnv_inv "ugj$name" "$OUT/oi_ugj$name.spec" --panphorte-prefix "$OUT/oi_ugj$name")|"
done
check "joined rows, a rotated copy: V inverted (two strands), V direct, Vr direct, Vr inverted (CN:CNINV)" "$ug_join" \
  "INVDUP 1:0 2:1|INVDUP 1:0 2:1|- 1:0 2:0|- 1:0 2:0|INVDUP 1:0 2:1|"

# ---------------------------------------------------------------------------------------------
# --diploid --min-alt-af on a copy-number record: S1's two contigs cross the site with different copy
# numbers, so its slot is '.' and leaves AN. The record is written with AF 1/2 and kept at 0.5.
# ---------------------------------------------------------------------------------------------
{ printf 'N a 60 321\nN L 30 322\nN b 60 323\n'
  printf 'P REF#0#chr1:1-180 a+,L+,L+,b+\nP S1#1#a a+,L+,L+,b+\nP S1#1#b a+,L+,L+,L+,b+\n'
  printf 'P S2#1#c a+,L+,L+,L+,L+,b+\nP S3#1#c a+,L+,L+,b+\nB a b L\n'; } > "$OUT/dcn.spec"
python3 -I "$TOOLS/spec_gfa.py" "$OUT/dcn.spec" "$OUT/dcn"
dcn=""
for af in 0.5 0.51; do
  "$BIN" call -i "$OUT/dcn/g.gfa" -c "$OUT/dcn/g.bubbles.csv" -r 'REF#0#chr1:1-180' -o "$OUT/dcn/af$af" --diploid \
    --min-alt-af "$af" --no-per-bubble-vcf -q >/dev/null 2>&1
  dcn="$dcn$? $(grep -vc '^#' "$OUT/dcn/af$af.region.vcf") $(grep -v '^#' "$OUT/dcn/af$af.region.vcf" | grep -o 'AN=[0-9]*;AC=[0-9]*')|"
done
check "--diploid --min-alt-af 0.5 / 0.51 on a CNV written with AN=2, AC=1 (exit, records, AN;AC)" "$dcn" "0 1 AN=2;AC=1|0 0 |"

# ---------------------------------------------------------------------------------------------
# 'with records' counts ordinary sites, as 'sites' does, also after --min-alt-af drops a record: here an
# inversion site holds site 1, whose deletion (AF 1/8) goes at 0.2, and the inversion record carries the
# inversion site's id.
# ---------------------------------------------------------------------------------------------
cat > "$OUT/wr.spec" <<'SPEC'
N L 1500 11
N R 400 12
N A 2500 13
N D 80 14
N E 90 15
N B 2500 16
N Rt 1500 17
P ref#0#chr1 L+,R+,A+,D+,B+,R-,Rt+
P s1#1#chr1 L+,R+,A+,B+,R-,Rt+
P s1#2#chr1 L+,R+,A+,D+,E+,B+,R-,Rt+
P s2#1#chr1 L+,R+,A+,D+,E+,B+,R-,Rt+
P s2#2#chr1 L+,R+,B-,E-,D-,A-,R-,Rt+
P s3#1#chr1 L+,R+,B-,E-,D-,A-,R-,Rt+
P s3#2#chr1 L+,R+,A+,D+,B+,R-,Rt+
P s4#1#chr1 L+,R+,B-,D-,A-,R-,Rt+
SPEC
python3 -I "$TOOLS/spec_gfa.py" "$OUT/wr.spec" "$OUT/wr"
"$BIN" bubble -i "$OUT/wr/g.gfa" -r ref -o "$OUT/wr/b" -q >/dev/null 2>&1
wr=""
for af in 0 0.2; do
  wr="$wr$("$BIN" call -i "$OUT/wr/b.sorted.gfa" -c "$OUT/wr/b.bubbles.csv" -r ref -o "$OUT/wr/af$af" --min-alt-af "$af" \
             --no-per-bubble-vcf 2>&1 | grep -E '\] (sites|records): ' | sed 's/^\[call [0-9:]*\] //' | cut -d';' -f1,2 | tr '\n' ' ')| "
done
check "an inversion site holding a site: the log's sites and records at --min-alt-af 0 | 0.2" "$wr" \
  "sites: 1; with records: 1 records: 3; DEL: 1 | sites: 1; with records: 1 records: 2; DEL: 0 | "

# ---------------------------------------------------------------------------------------------
# Many sites under a low open-file limit: each per-site VCF is written and closed before the next is
# opened, so the run does not need a descriptor per site. The log counts the files it wrote.
# ---------------------------------------------------------------------------------------------
d="$OUT/many"; mkdir -p "$d"
awk -v n=80 'BEGIN { b = "ACGT"; x = 4242
    printf "H\tVN:Z:1.0\n"
    for (i = 0; i <= n; i++) { s = ""; for (k = 0; k < 100; k++) { x = (x * 16807) % 2147483647; s = s substr(b, int(x * 4 / 2147483647) + 1, 1) }
      printf "S\tf%d\t%s\n", i, s }
    for (i = 0; i < n; i++) { s = ""; for (k = 0; k < 60; k++) { x = (x * 16807) % 2147483647; s = s substr(b, int(x * 4 / 2147483647) + 1, 1) }
      printf "S\ti%d\t%s\nL\tf%d\t+\tf%d\t+\t0M\nL\tf%d\t+\ti%d\t+\t0M\nL\ti%d\t+\tf%d\t+\t0M\n", i, s, i, i + 1, i, i, i, i + 1
      r = r "f" i "+,"; a = a "f" i "+,i" i "+," }
    printf "P\tref#0#chr1:1-%d\t%sf%d+\t*\nP\th1#1#c\t%sf%d+\t*\nP\th2#1#c\t%sf%d+\t*\n", (n + 1) * 100, r, n, a, n, r, n }' > "$d/g.gfa"
{ printf 'bubble_id,source,source_orient,sink,sink_orient,path_support,min_inside_bp,max_inside_bp,inside_nodes\n'
  for i in $(seq 0 79); do printf '%d,f%d,+,f%d,+,3,0,60,i%d\n' $((i + 1)) "$i" $((i + 1)) "$i"; done; } > "$d/g.csv"
( ulimit -n 40; "$BIN" call -i "$d/g.gfa" -c "$d/g.csv" -r 'ref#0#chr1:1-8100' -o "$d/v" --threads 2 > "$d/log" 2>&1 ); rc=$?
check "80 sites with records under an open-file limit of 40: exit status, per-site VCFs, files the log counts, files written" \
  "$rc $(ls "$d" | grep -c '^v\.bubble_.*\.vcf$') $(grep -o 'wrote [0-9]* files' "$d/log") $(ls "$d" | grep -c '^v\.')" \
  "0 80 wrote 83 files 83"
"$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" -o "$OUT/short_t" -t 1 -q >/dev/null 2>&1
check "-t is --threads (exit status, records)" "$? $(grep -vc '^#' "$OUT/short_t.region.vcf" 2>/dev/null)" "0 3"
check "the log names no per-site VCF where none was written (records dropped: --min-alt-af 0.51)" \
  "$(grep -c 'bubble_<id>' "$OUT/maf/af0.51.log") $(grep -o 'wrote [0-9]* files' "$OUT/maf/af0.51.log")" "0 wrote 3 files"

# ---------------------------------------------------------------------------------------------
# A unit-graph table from another run, whose node ids all exist here, is refused: its one-node unit does
# not spell the node's sequence, so it would turn a deletion into copies of a unit of another length.
# ---------------------------------------------------------------------------------------------
del_node=$(awk -F'\t' '$1 == "S" && $3 == "'"$N2"'" { print $2 }' "$OUT/cb.sorted.gfa")
flank_node=$(awk -F'\t' '$1 == "S" && $3 == "'"$N3"'" { print $2 }' "$OUT/cb.sorted.gfa")
printf 'input_bubble_id\tunit_bp\tcanonical_motif\tcopy_quantum\tunit\tanchor_pos\tanchor_nodes\tn_nodes\tnodes\n1\t60\t%s\t1\t%s\t0\t%s\t1\t%s\n' \
  "$(blk 60 41)" "$(blk 60 41)" "$del_node" "$del_node" > "$OUT/stale.panphorte.unit_graph.tsv"
refuses "a unit-graph table that does not describe the graph is refused" "does not spell its unit" \
  "$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" --panphorte-prefix "$OUT/stale" -r "$REF_NAME" -q
printf 'input_bubble_id\tunit_bp\tanchor_pos\tanchor_nodes\tn_nodes\tnodes\n1\t70\t0\t%s\t1\t%s\n' "$flank_node" "$del_node" \
  > "$OUT/stale2.panphorte.unit_graph.tsv"
refuses "a unit-graph table whose anchor is not among its nodes is refused" "is not among its nodes" \
  "$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" --panphorte-prefix "$OUT/stale2" -r "$REF_NAME" -q

# ---------------------------------------------------------------------------------------------
# A gene between two copies of a repeat that are not adjacent is not in the copies: U is walked twice by
# the reference, around gene G on A, and the haplotype gains a U. G is the site's only gene and lies beside
# the copies, so it has no row, not the repeat's count.
# ---------------------------------------------------------------------------------------------
cat > "$OUT/gbc.spec" <<'SPEC'
N S 100 331
N U 100 332
N A 1000 333
N K 100 334
P ref#1#chr1:1-1400 S+,U+,A+,U+,K+
P hap#1#chr1 S+,U+,U+,A+,U+,K+
B S K U;A
SPEC
python3 -I "$TOOLS/spec_gfa.py" "$OUT/gbc.spec" "$OUT/gbc"
printf 'chr1\tt\tgene\t201\t1200\t.\t+\t.\tgene_id "G"; gene_name "G";\n' > "$OUT/gbc/g.gtf"
"$BIN" call -i "$OUT/gbc/g.gfa" -c "$OUT/gbc/g.bubbles.csv" -r 'ref#1#chr1:1-1400' --gtf "$OUT/gbc/g.gtf" -o "$OUT/gbc/v" \
  --no-per-bubble-vcf -q >/dev/null 2>&1
check "a gene between non-adjacent repeat copies: the CNV record, and dup_gene_cn rows" \
  "$(grep -v '^#' "$OUT/gbc/v.region.vcf" | grep -c 'CN_METHOD=UNIT') $(tail -n +2 "$OUT/gbc/v.dup_gene_cn.tsv" | wc -l | tr -d ' ')" "1 0"

# ---------------------------------------------------------------------------------------------
# An assembly gap the graph cut into nodes of 9 N each is a gap: hA crosses it forward, hR stored on the
# other strand, both '.' and not a 90 bp insertion of N. hC's 9 N on a node of their own are no gap.
# ---------------------------------------------------------------------------------------------
{ printf 'N S 100 341\nN K 100 342\nN N0 =NNNNNNNNN\n'
  for i in 1 2 3 4 5 6 7 8 9 10; do printf 'N N%d =NNNNNNNNN\n' "$i"; done
  printf 'P ref#0#chr1:1-200 S+,K+\nP hA#1#chr1 S+,N1+,N2+,N3+,N4+,N5+,N6+,N7+,N8+,N9+,N10+,K+\n'
  printf 'P hR#1#chr1 K-,N10-,N9-,N8-,N7-,N6-,N5-,N4-,N3-,N2-,N1-,S-\nP hB#1#chr1 S+,K+\nP hC#1#chr1 S+,N0+,K+\n'
  printf 'B S K N0;N1;N2;N3;N4;N5;N6;N7;N8;N9;N10\n'; } > "$OUT/sgap.spec"
python3 -I "$TOOLS/spec_gfa.py" "$OUT/sgap.spec" "$OUT/sgap"
"$BIN" call -i "$OUT/sgap/g.gfa" -c "$OUT/sgap/g.bubbles.csv" -r 'ref#0#chr1:1-200' -o "$OUT/sgap/v" \
  --no-per-bubble-vcf > "$OUT/sgap/log" 2>&1
check "an assembly gap split into nodes under 10 N: records, haplotypes set '.' across a gap" \
  "$(grep -vc '^#' "$OUT/sgap/v.region.vcf") $(grep -o 'GT . across an assembly gap: [0-9]*' "$OUT/sgap/log")" \
  "0 GT . across an assembly gap: 2"
# Gaps are read per walk. gs3: hA's N1+N2 (18 N) is a gap, but the reference's and hC's N1 alone (9 N) is
# not, so the site is typed and hB's deletion of D called. gs2: hC's 9 N on N1 keep it callable although hA
# walks N1 inside its gap.
gsplit() { # gsplit <name> <ref steps> <hA> <hB> <hC> [hD]: hA..hD's GT on the DEL, AN, gap count
  { printf 'N S 100 341\nN D 300 345\nN K 100 342\nN N1 =NNNNNNNNN\nN N2 =NNNNNNNNN\n'
    printf 'P ref#0#chr1:1-600 %s\nP hA#1#chr1 %s\nP hB#1#chr1 %s\nP hC#1#chr1 %s\n' "$2" "$3" "$4" "$5"
    [ -n "${6:-}" ] && printf 'P hD#1#chr1 %s\n' "$6"
    printf 'B S K D;N1;N2\n'; } > "$OUT/$1.spec"
  python3 -I "$TOOLS/spec_gfa.py" "$OUT/$1.spec" "$OUT/$1"
  "$BIN" call -i "$OUT/$1/g.gfa" -c "$OUT/$1/g.bubbles.csv" -r 'ref#0#chr1:1-600' -o "$OUT/$1/v" --no-per-bubble-vcf \
    > "$OUT/$1/log" 2>&1
  local h out=""
  for h in hA hB hC hD; do out="$out$(gt_of "$OUT/$1/v.region.vcf" "$h#" DEL)"; done
  VCF="$OUT/$1/v.region.vcf"
  printf '%s %s %s' "${out:--}" "$(info DEL AN)" "$(grep -c 'assembly gap' "$OUT/$1/log")"; }
check "a 9 N reference beside an 18 N haplotype on shared nodes: hA..hC GT, AN, gap lines" \
  "$(gsplit gs3 S+,N1+,D+,K+ S+,N1+,N2+,D+,K+ S+,N1+,K+ S+,N1+,D+,K+)" ".10 2 1"
check "a 9 N haplotype on a node of an 18 N one: hA..hD GT, AN, gap lines" \
  "$(gsplit gs2 S+,D+,K+ S+,N1+,N2+,K+ S+,K+ S+,N1+,D+,K+ S+,D+,K+)" ".100 3 1"
# A gap split across two partly-N nodes (X ends in 5 N, Y starts with 5): h1, and h2 stored on the other
# strand, cross it and are '.'; h3's 5 N on X alone are no gap, and its 65 bp insertion is called.
{ printf 'N S 100 341\nN K 100 342\nN X =%sNNNNN\nN Y =NNNNN%s\n' "$(blk 60 351)" "$(blk 60 352)"
  printf 'P ref#0#chr1:1-200 S+,K+\nP h1#1#chr1 S+,X+,Y+,K+\nP h2#1#chr1 K-,Y-,X-,S-\nP h3#1#chr1 S+,X+,K+\nB S K X;Y\n'; } \
  > "$OUT/pgap.spec"
python3 -I "$TOOLS/spec_gfa.py" "$OUT/pgap.spec" "$OUT/pgap"
"$BIN" call -i "$OUT/pgap/g.gfa" -c "$OUT/pgap/g.bubbles.csv" -r 'ref#0#chr1:1-200' -o "$OUT/pgap/v" \
  --no-per-bubble-vcf > "$OUT/pgap/log" 2>&1
VCF="$OUT/pgap/v.region.vcf"
check "a gap split across two partly-N nodes: records, h1 h2 h3 GT, SVLEN, AN, gap line" \
  "$(grep -vc '^#' "$VCF") $(gt_of "$VCF" h1 INS)$(gt_of "$VCF" h2 INS)$(gt_of "$VCF" h3 INS) $(info INS SVLEN) $(info INS AN) $(grep -o 'GT . across an assembly gap: [0-9]*' "$OUT/pgap/log")" \
  "1 ..1 65 1 GT . across an assembly gap: 2"
"$BIN" call -i "$OUT/pgap/g.gfa" -c "$OUT/pgap/g.bubbles.csv" -r 'ref#0#chr1:1-200' -o "$OUT/pgap/one" \
  --no-per-bubble-vcf --no-variant-nodes > "$OUT/pgap/one.log" 2>&1
check "one file written: the log says 'file'" "$(grep -o 'wrote 1 file[s]* in' "$OUT/pgap/one.log")" "wrote 1 file in"

# Node comparison: h1 walks node b in place of the reference's a, both the single base G, and skips the 100 bp
# node D after a. The two sides of the replacement share that G, so it is written as its minimal change: a plain
# 100 bp deletion, left-aligned (D ends in G, so POS 300 and END 400), with no SUBSTITUTION flag.
printf 'N L 300 91\nN a =G\nN b =G\nN D 100 92\nN R 300 93\nP ref#0#chr1:0-701 L+,a+,D+,R+\nP h1#1#chr1 L+,b+,R+\nP h2#1#chr1 L+,a+,D+,R+\nB L R a;b;D\n' \
  > "$OUT/ntrim.spec"
python3 -I "$TOOLS/spec_gfa.py" "$OUT/ntrim.spec" "$OUT/ntrim"
"$BIN" call -i "$OUT/ntrim/g.gfa" -c "$OUT/ntrim/g.bubbles.csv" -r 'ref#0#chr1:0-701' -o "$OUT/ntrim/v" \
  --no-align-alleles --no-per-bubble-vcf -q >/dev/null 2>&1
check "node replacement sharing a base with the reference: a plain deletion (POS END SVLEN flag)" \
  "$(awk -F'\t' '!/^#/{n=split($8,f,";"); e=""; l=""; for(i=1;i<=n;i++){split(f[i],kv,"=");
       if(kv[1]=="END") e=kv[2]; if(kv[1]=="SVLEN") l=kv[2]}
       print $2, e, l, ($8 ~ /SUBSTITUTION/ ? "SUBSTITUTION" : "plain")}' "$OUT/ntrim/v.region.vcf")" "300 400 -100 plain"
# rinfo_of <vcf> <key>: the first record's INFO value for key
rinfo_of() { awk -F'\t' -v k="$2" '!/^#/{n=split($8,f,";"); for(i=1;i<=n;i++){split(f[i],kv,"="); if(kv[1]==k){print kv[2]; exit}}}' "$1"; }
D=$(awk -F'\t' '$1=="S" && $2=="D"{print $3}' "$OUT/ntrim/g.gfa")
check "the deletion slides one base left (a's G, then D less its last base) and 4 bp in all (HOMLEN CIPOS CIEND DELSEQ)" \
  "$(rinfo_of "$OUT/ntrim/v.region.vcf" HOMLEN) $(rinfo_of "$OUT/ntrim/v.region.vcf" CIPOS) $(rinfo_of "$OUT/ntrim/v.region.vcf" CIEND) $(rinfo_of "$OUT/ntrim/v.region.vcf" DELSEQ)" \
  "4 0,4 0,4 G${D%?}"
# h2 deletes D through the reference's own a; h1's replacement is the same deletion, so one record holds both.
sed 's/^P h2#1#chr1 .*/P h2#1#chr1 L+,a+,R+/' "$OUT/ntrim.spec" > "$OUT/ndup.spec"
python3 -I "$TOOLS/spec_gfa.py" "$OUT/ndup.spec" "$OUT/ndup"
"$BIN" call -i "$OUT/ndup/g.gfa" -c "$OUT/ndup/g.bubbles.csv" -r 'ref#0#chr1:0-701' -o "$OUT/ndup/v" \
  --no-align-alleles --no-per-bubble-vcf -q >/dev/null 2>&1
check "a node replacement that is a plain deletion merges with the same deletion (records, h1 h2 GT)" \
  "$(grep -vc '^#' "$OUT/ndup/v.region.vcf") $(awk -F'\t' '/^#CHROM/{for(i=10;i<=NF;i++) c[$i]=i} !/^#/{split($c["h1#1#chr1"],a,":"); split($c["h2#1#chr1"],b,":"); print a[1] b[1]; exit}' "$OUT/ndup/v.region.vcf")" "1 11"
# h1 walks b, a's G and 97 bp more, for a: a plain 97 bp insertion, classified as one (INS_SUBTYPE).
printf 'N L 300 91\nN a =G\nN I 90 95\nN R 300 93\nP ref#0#chr1:0-601 L+,a+,R+\nP h2#1#chr1 L+,a+,I+,R+\n' > "$OUT/ncls.spec"
printf 'N b =G%s\nP h1#1#chr1 L+,b+,R+\nB L R a;b;I\n' "$(blk 97 96)" >> "$OUT/ncls.spec"
python3 -I "$TOOLS/spec_gfa.py" "$OUT/ncls.spec" "$OUT/ncls"
"$BIN" call -i "$OUT/ncls/g.gfa" -c "$OUT/ncls/g.bubbles.csv" -r 'ref#0#chr1:0-601' -o "$OUT/ncls/v" \
  --no-align-alleles --classify-ins --no-per-bubble-vcf -q >/dev/null 2>&1
check "--classify-ins: both insertions, the node replacement's too, are classified (SVLEN:INS_SUBTYPE)" \
  "$(awk -F'\t' '!/^#/{n=split($8,f,";"); l=""; t="-"; for(i=1;i<=n;i++){split(f[i],kv,"="); if(kv[1]=="SVLEN") l=kv[2]; if(kv[1]=="INS_SUBTYPE") t=kv[2]} printf "%s:%s ", l, t}' "$OUT/ncls/v.region.vcf")" \
  "90:NOVEL 97:NOVEL "
# h1 walks x, a copy of R, in place of the reference's 20 bp node a: a replacement of net +280, written as an INS
# flagged SUBSTITUTION. Its inserted bases are classified too, and are a duplication of R.
printf 'N L 300 91\nN a 20 97\nN R 300 93\nN x @R:0:1\nP ref#0#chr1:0-621 L+,a+,R+\nP h1#1#chr1 L+,x+,R+\nP h2#1#chr1 L+,a+,R+\nB L R a;x\n' \
  > "$OUT/nsub.spec"
python3 -I "$TOOLS/spec_gfa.py" "$OUT/nsub.spec" "$OUT/nsub"
"$BIN" call -i "$OUT/nsub/g.gfa" -c "$OUT/nsub/g.bubbles.csv" -r 'ref#0#chr1:0-621' -o "$OUT/nsub/v" \
  --classify-ins --no-per-bubble-vcf -q >/dev/null 2>&1
check "--classify-ins: an INS flagged SUBSTITUTION is classified from its inserted bases (SVTYPE:flag:INS_SUBTYPE)" \
  "$(awk -F'\t' '!/^#/{n=split($8,f,";"); t="-"; for(i=1;i<=n;i++){split(f[i],kv,"="); if(kv[1]=="INS_SUBTYPE") t=kv[2]}
       printf "%s:%s:%s ", ($8 ~ /SVTYPE=INS;/ ? "INS" : "other"), ($8 ~ /SUBSTITUTION/ ? "SUBSTITUTION" : "plain"), t}' "$OUT/nsub/v.region.vcf")" \
  "INS:SUBSTITUTION:DUP "
# The reference runs X1 D X2 with X1 and X2 the same 30 bases. Node by node, h1a/h1b delete D X2 and h2a/h2b
# delete X1 D: two placements of one deletion, written the same once left-aligned (130 bp after F). h3a..h3c
# delete D less its first 10 bases and X2 (120 bp). The representative is the event most carriers hold as
# written, four, not the three that are one placement alone.
{ printf 'N F 300 111\nN X1 =%s\nN D 100 113\nN X2 =%s\nN D3 @D:0:1\nN R 300 115\n' "$(blk 30 112)" "$(blk 30 112)"
  printf 'P ref#0#chr1:0-760 F+,X1+,D+,X2+,R+\nP hr#1#c F+,X1+,D+,X2+,R+\n'
  printf 'P h1a#1#c F+,X1+,R+\nP h1b#1#c F+,X1+,R+\nP h2a#1#c F+,X2+,R+\nP h2b#1#c F+,X2+,R+\n'
  printf 'P h3a#1#c F+,X1+,D3+,R+\nP h3b#1#c F+,X1+,D3+,R+\nP h3c#1#c F+,X1+,D3+,R+\nB F R X1;D;X2;D3\n'; } > "$OUT/wrep.spec"
python3 -I "$TOOLS/spec_gfa.py" "$OUT/wrep.spec" "$OUT/wrep"
D10=$(awk -F'\t' '$1=="S" && $2=="D"{print substr($3, 1, 10)}' "$OUT/wrep/g.gfa")
awk -F'\t' -v OFS='\t' -v d="$D10" '$1=="S" && $2=="D3"{$3=d} {print}' "$OUT/wrep/g.gfa" > "$OUT/wrep/g10.gfa"
"$BIN" call -i "$OUT/wrep/g10.gfa" -c "$OUT/wrep/g.bubbles.csv" -r 'ref#0#chr1:0-760' -o "$OUT/wrep/v" \
  --no-align-alleles --no-per-bubble-vcf -q >/dev/null 2>&1
check "the representative is the event most carriers hold as written (SVLEN:SVLEN_RANGE:AC)" \
  "$(awk -F'\t' '!/^#/{n=split($8,f,";"); for(i=1;i<=n;i++){split(f[i],kv,"="); v[kv[1]]=kv[2]} print v["SVLEN"]":"v["SVLEN_RANGE"]":"v["AC"]}' "$OUT/wrep/v.region.vcf")" \
  "-130:-130,-120:7"
# A haplotype in two contigs at a module site: s1#1#c1 deletes the second A2 and X, s1#1#c2 replaces that A2 by
# Q, a 2-base variant. Only c1 carries the deletion; per path the haplotype counts as a carrier (no AC=0 record),
# and with --diploid sample s1 carries it with c1's values.
{ printf 'N L 2000 81\nN A3 20 82\nN A2 40 83\nN Q @A2:2:7\nN S 40 84\nN X 60 85\nN R 2000 86\n'
  printf 'P ref#0#chr1:0-4220 L+,A3+,A2+,S+,A3+,A2+,X+,R+\nP hfull#1#chr1 L+,A3+,A2+,X+,R+\n'
  printf 'P hgain#1#chr1 L+,A3+,A2+,S+,A3+,A2+,S+,A3+,A2+,X+,R+\nP s1#1#c1 L+,A3+,A2+,S+,A3+,R+\nP s1#1#c2 L+,A3+,A2+,S+,A3+,Q+,X+,R+\n'
  for i in 1 2 3 4 5 6 7 8 9 10; do printf 'P hr%s#1#chr1 L+,A3+,A2+,S+,A3+,A2+,X+,R+\n' "$i"; done
  printf 'B L R A3;A2;S;X;Q\n'; } > "$OUT/acsplit.spec"
python3 -I "$TOOLS/spec_gfa.py" "$OUT/acsplit.spec" "$OUT/acsplit"
for m in path diploid; do x=""; [ "$m" = diploid ] && x="--diploid"
  # shellcheck disable=SC2086
  "$BIN" call -i "$OUT/acsplit/g.gfa" -c "$OUT/acsplit/g.bubbles.csv" -r 'ref#0#chr1:0-4220' -o "$OUT/acsplit_$m" --cn $x \
    --no-per-bubble-vcf -q >/dev/null 2>&1; done
check "a haplotype one of whose contigs carries a record is its carrier (DEL_A2 AN;AC per path | per sample, s1 GT)" \
  "$(grep -v '^#' "$OUT/acsplit_path.region.vcf" | awk -F'\t' '$3=="bubble1_DEL_A2"' | grep -o 'AN=[0-9]*;AC=[0-9]*')|$(awk -F'\t' '/^#CHROM/{for(i=10;i<=NF;i++) if($i=="s1") c=i} !/^#/ && $3=="bubble1_DEL_A2"{match($8,/AN=[0-9]*;AC=[0-9]*/); split($c,g,":"); print substr($8,RSTART,RLENGTH), g[1]}' "$OUT/acsplit_diploid.region.vcf")" \
  "AN=13;AC=1|AN=13;AC=1 1"
check "no record is written with AC=0 (per path, per sample)" \
  "$(cat "$OUT/acsplit_path.region.vcf" "$OUT/acsplit_diploid.region.vcf" | grep -v '^#' | grep -c ';AC=0;')" "0"

# --min-alt-af on a copy-number record, tested before the record is written, counts as the writer does: the
# reference path is out of AN. hgain gains a copy and hr does not: AF 1/2, so 0.5 keeps the CNV.
printf 'N L 2000 81\nN A3 20 82\nN A2 40 83\nN S 40 84\nN X 60 85\nN R 2000 86\nP ref#0#chr1:0-4220 L+,A3+,A2+,S+,A3+,A2+,X+,R+\nP hgain#1#chr1 L+,A3+,A2+,S+,A3+,A2+,S+,A3+,A2+,X+,R+\nP hr#1#chr1 L+,A3+,A2+,S+,A3+,A2+,X+,R+\n' > "$OUT/cnaf.spec"
python3 -I "$TOOLS/spec_gfa.py" "$OUT/cnaf.spec" "$OUT/cnaf"
"$BIN" bubble -i "$OUT/cnaf/g.gfa" -r ref --min-variant-bp 0 -o "$OUT/cnaf.b" --merge-nearby-bp 0 -q >/dev/null 2>&1
"$BIN" call -i "$OUT/cnaf.b.sorted.gfa" --bubble-prefix "$OUT/cnaf.b" -r ref -o "$OUT/cnaf" --cn --min-alt-af 0.5 -q >/dev/null 2>&1
check "--min-alt-af 0.5 keeps a CNV written with AF 1/2 (records: SVTYPE AF)" \
  "$(grep -v '^#' "$OUT/cnaf.region.vcf" | awk -F'\t' '{n=split($8,f,";"); for(i=1;i<=n;i++){split(f[i],kv,"="); v[kv[1]]=kv[2]} printf "%s %s ", v["SVTYPE"], v["AF"]}')" "CNV 0.5 "

# A 52 bp unit U, six copies on the reference, then D, a node spelling U, and E. h2 skips D; h1 skips it too but
# walks F, a node spelling E, in place of E. h1's replacement is the deletion of D, which slides along the copies,
# but no further than where the copies left the walks: it is h2's deletion, one record.
U52=$(blk 52 103)
{ printf 'N L 300 101\nN A 100 102\nN U =%s\nN D =%s\nN E =%s\nN F =%s\nN R 300 105\n' "$U52" "$U52" "$(blk 30 108)" "$(blk 30 108)"
  M="U+,U+,U+,U+,U+,U+"
  printf 'P ref#0#chr1:0-1094 L+,A+,%s,D+,E+,R+\nP h1#1#chr1 L+,A+,%s,F+,R+\nP h2#1#chr1 L+,A+,%s,E+,R+\n' "$M" "$M" "$M"
  printf 'P h3#1#chr1 L+,A+,%s,U+,U+,D+,E+,R+\nB L R A;U;D;E;F\n' "$M"; } > "$OUT/nfloor.spec"
python3 -I "$TOOLS/spec_gfa.py" "$OUT/nfloor.spec" "$OUT/nfloor"
"$BIN" call -i "$OUT/nfloor/g.gfa" -c "$OUT/nfloor/g.bubbles.csv" -r 'ref#0#chr1:0-1094' -o "$OUT/nfloor/v" \
  --no-align-alleles --no-per-bubble-vcf -q >/dev/null 2>&1
check "a node replacement's deletion beside counted copies merges with the same deletion (DEL records, AC)" \
  "$(grep -v '^#' "$OUT/nfloor/v.region.vcf" | grep -c 'SVTYPE=DEL;') $(awk -F'\t' '!/^#/ && /SVTYPE=DEL;/{n=split($8,f,";"); for(i=1;i<=n;i++) if(f[i] ~ /^AC=/) print f[i]}' "$OUT/nfloor/v.region.vcf")" "1 AC=2"
# The same with the walks stored on the other strand (the reference walked backwards): the copies' side is then
# the walks' later side, and the slide stops there all the same.
{ printf 'N L 300 101\nN A 100 102\nN U =%s\nN D =%s\nN E =%s\nN F =%s\nN R 300 105\n' "$U52" "$U52" "$(blk 30 108)" "$(blk 30 108)"
  M="U-,U-,U-,U-,U-,U-"
  printf 'P ref#0#chr1:0-1094 R-,%s,D-,E-,A-,L-\nP h1#1#chr1 R-,%s,F-,A-,L-\nP h2#1#chr1 R-,%s,E-,A-,L-\n' "$M" "$M" "$M"
  printf 'P h3#1#chr1 R-,U-,U-,%s,D-,E-,A-,L-\nB L R A;U;D;E;F\n' "$M"; } > "$OUT/nfloorr.spec"
python3 -I "$TOOLS/spec_gfa.py" "$OUT/nfloorr.spec" "$OUT/nfloorr"
"$BIN" call -i "$OUT/nfloorr/g.gfa" -c "$OUT/nfloorr/g.bubbles.csv" -r 'ref#0#chr1:0-1094' -o "$OUT/nfloorr/v" \
  --no-align-alleles --no-per-bubble-vcf -q >/dev/null 2>&1
check "the same, the reference walked backwards (DEL records, AC)" \
  "$(grep -v '^#' "$OUT/nfloorr/v.region.vcf" | grep -c 'SVTYPE=DEL;') $(awk -F'\t' '!/^#/ && /SVTYPE=DEL;/{n=split($8,f,";"); for(i=1;i<=n;i++) if(f[i] ~ /^AC=/) print f[i]}' "$OUT/nfloorr/v.region.vcf")" "1 AC=2"

# Aligned alleles: the band is set from the node comparison's drift, which can understate the change where the
# walks are matched at short nodes out of place. The best placement is a 123 bp deletion and an 86 bp insertion
# (net -37); a band of the node drift alone holds only a worse one, a 71 -> 78 bp replacement.
python3 -I "$TOOLS/spec_gfa.py" "$(cd "$(dirname "$0")/.." && pwd)/synthetic_data/call_cases/narrow_band.spec" "$OUT/nband"
"$BIN" call -i "$OUT/nband/g.gfa" -c "$OUT/nband/g.bubbles.csv" -r 'ref#0#chr1:0-1186' -o "$OUT/nband/v" \
  --no-per-bubble-vcf -q >/dev/null 2>&1
check "aligned alleles: the placement a wider band finds (records SVTYPE:SVLEN)" \
  "$(awk -F'\t' '!/^#/{n=split($8,f,";"); for(i=1;i<=n;i++) if(f[i] ~ /^SVTYPE=|^SVLEN=/) s=s f[i] " "; } END{print s}' "$OUT/nband/v.region.vcf")" \
  "SVTYPE=DEL SVLEN=-123 SVTYPE=INS SVLEN=86 "

printf "\n"
if [ "$fails" -eq 0 ]; then printf "call_stats: all assertions passed\n"; exit 0; fi
# repeated last, shortened: a long output can be cut in the middle by ctest
printf "call_stats: %d assertion(s) failed:\n" "$fails"; printf "  FAIL %s\n" "${failed[@]}"; exit 1
