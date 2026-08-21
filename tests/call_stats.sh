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
bad() { printf "  FAIL %s\n" "$1"; fails=$((fails + 1)); }

# Deterministic sequence, distinct per node: repeat a fixed 40-mer rotated by $2 so no two nodes
# share content (identical flanks would give the snarl finder a different graph than intended).
blk() { local n=$1 r=${2:-0} u="ACGTTGCAATTCCGGATCAGGTCAAGCTTGACCTAGGACT" out=""
        u="${u:$r}${u:0:$r}"
        while [ ${#out} -lt "$n" ]; do out="$out$u"; done; printf '%s' "${out:0:$n}"; }

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

"$BIN" bubble -i "$OUT/coord.gfa" -r "$REF_NAME" -o "$OUT/cb" -q >/dev/null 2>&1
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
# against everything and silently disabled the guard it was meant to set.
for bad in nan 2 -0.5; do
  refuses "--max-dup-region-frac $bad is refused" "must be in [0,1]" \
    "$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" --max-dup-region-frac "$bad" -q
done
refuses "--minimap-preset must be one of the advertised three" "must be asm5" \
  "$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" --minimap-preset bogus -q

# --min-alt-af is the honest name; --min-maf stays accepted so existing pipelines do not break.
"$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" -o "$OUT/af1" --min-alt-af 0.1 -q >/dev/null 2>&1
"$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" -o "$OUT/af2" --min-maf 0.1 -q >/dev/null 2>&1
if [ -s "$OUT/af1.region.vcf" ] && cmp -s "$OUT/af1.region.vcf" "$OUT/af2.region.vcf"; then
  ok "--min-alt-af and its --min-maf alias agree"
else bad "--min-alt-af and --min-maf disagree (or one failed)"; fi

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
"$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" -o "$OUT/t1" --threads 1 --allele-vcf -q >/dev/null 2>&1
"$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" -o "$OUT/t8" --threads 8 --allele-vcf -q >/dev/null 2>&1
if cmp -s "$OUT/t1.region.vcf" "$OUT/t8.region.vcf" && cmp -s "$OUT/t1.alleles.vcf" "$OUT/t8.alleles.vcf"; then
  ok "1 thread and 8 threads give byte-identical region and allele VCFs"
else bad "output differs between 1 and 8 threads"; fi

# The allele VCF is the lossless record and must not depend on anything being interpretable. Raising
# --min-sv-bp above every event empties the region VCF; the allele VCF must still describe the site.
"$BIN" call -i "$OUT/cb.sorted.gfa" -c "$OUT/cb.bubbles.csv" -r "$REF_NAME" -o "$OUT/nolo" \
      --min-sv-bp 100000 --allele-vcf -q >/dev/null 2>&1
nolo_region=$(grep -v "^#" "$OUT/nolo.region.vcf" 2>/dev/null | wc -l | tr -d " ")
nolo_alleles=$(grep -v "^#" "$OUT/nolo.alleles.vcf" 2>/dev/null | wc -l | tr -d " ")
[ "$nolo_region" = "0" ] && [ "$nolo_alleles" -gt 0 ] \
  && ok "the allele VCF survives a region VCF with no interpreted calls ($nolo_alleles records)" \
  || bad "allele-VCF independence: region=$nolo_region alleles=$nolo_alleles (want region 0, alleles > 0)"

# ---------------------------------------------------------------------------------------------
# Per-gene copy number: a haplotype with no evidence must not be published as a confident zero.
#
# A self-loop REP node carrying two near-identical genes fires the REP DUP route, and --gtf then
# splits it per gene. hapD deletes the module entirely, so it has no hits at any divergent site and
# no basis for a split. It used to be reported CN=0 reliable=1 -- indistinguishable from a measured
# absence -- because separability was decided from the pair's marker sets, which exist for every
# haplotype, rather than from whether THIS haplotype resolved.
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

"$BIN" bubble -i "$OUT/gene.gfa" -r 'gref#0#chr9' -o "$OUT/gb" -q >/dev/null 2>&1
"$BIN" call -i "$OUT/gb.sorted.gfa" -c "$OUT/gb.bubbles.csv" -r 'gref#0#chr9' -o "$OUT/gc" \
      --cn --gtf "$OUT/gene.gtf" -q >/dev/null 2>&1
GTSV="$OUT/gc.dup_gene_cn.tsv"
if [ ! -s "$GTSV" ]; then
  bad "the two-gene REP fixture produced no dup_gene_cn table"
else
  # a haplotype that DOES carry the module resolves and is reliable
  carrier_rel=$(awk -F'\t' '$3 ~ /ghapB/ && $4=="GENEA" {print $6; exit}' "$GTSV")
  [ "$carrier_rel" = "1" ] && ok "a haplotype carrying the module reports a reliable per-gene split" \
    || bad "carrier per-gene split: expected reliable=1, got ${carrier_rel:-<missing>}"
  # the haplotype that lost the module must NOT be reported as a confident zero
  lost_rel=$(awk -F'\t' '$3 ~ /ghapD/ && $4=="GENEA" {print $6; exit}' "$GTSV")
  lost_dos=$(awk -F'\t' '$3 ~ /ghapD/ && $4=="GENEA" {print $7; exit}' "$GTSV")
  [ "$lost_rel" = "0" ] \
    && ok "a haplotype with no per-site evidence falls back (reliable=0), not a confident zero" \
    || bad "module-loss haplotype: expected reliable=0, got ${lost_rel:-<missing>} (dosage ${lost_dos:-?})"
fi

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
RF1=$(blk 100 0); RR=$(blk 50 7); RF2=$(blk 100 13); RA=$(blk 70 21); RF3=$(blk 100 29)
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

"$BIN" bubble -i "$OUT/rep.gfa" -r "$RREF_NAME" -o "$OUT/rb" -q >/dev/null 2>&1
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
"$BIN" bubble -i "$OUT/rev.gfa" -r 'rvref#0#chr1:1-470' -o "$OUT/rvb" --no-flip -q >/dev/null 2>&1
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
# The PEAK copy-number route, which had no test of any kind. It cannot be checked against the
# pangene truth the other routes use: every PEAK record on the reference loci describes an 11-92 bp
# micro-module, and gene-level truth counts GENE copies, so the two measure different things. A
# fixture with a known micro-module copy number is the only way to test it.
#
# The reference visits u once; hapdup visits it three times, NON-consecutively so the REP self-loop
# route does not claim it first. u is 60 bp, so a carrier gains 120 bp and clears --min-sv-bp.
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
"$BIN" bubble -i "$OUT/peak.gfa" -r 'pkref#0#chr8:1-860' -o "$OUT/pkb" -q >/dev/null 2>&1
"$BIN" call -i "$OUT/pkb.sorted.gfa" -c "$OUT/pkb.bubbles.csv" -r 'pkref#0#chr8:1-860' \
      -o "$OUT/pkc" --cn -q >/dev/null 2>&1
PKV="$OUT/pkc.region.vcf"
pk_method=$(awk -F'\t' '$0!~/^#/ && $8~/SVTYPE=DUP/{n=split($8,a,";");
              for(i=1;i<=n;i++){split(a[i],kv,"="); if(kv[1]=="CN_METHOD"){print kv[2]; exit}}}' "$PKV")
[ "$pk_method" = "PEAK" ] \
  && ok "a non-consecutive micro-module duplication is called by the PEAK route" \
  || bad "expected CN_METHOD=PEAK, got '${pk_method:-<no DUP record>}'"
# the duplicating haplotype must read 3 copies, the flat one 1 -- a positive AND a negative
pk_cn() { awk -F'\t' -v s="$1" '/^#CHROM/{for(i=10;i<=NF;i++) if($i ~ s) c=i}
          $0!~/^#/ && $8~/SVTYPE=DUP/{split($9,k,":"); for(j in k) if(k[j]=="CN") ci=j
            split($c,v,":"); print v[ci]; exit}' "$PKV"; }
[ "$(pk_cn pkdup)" = "3" ] \
  && ok "PEAK reports 3 copies for the haplotype that visits the unit three times" \
  || bad "PEAK CN for the duplicating haplotype: expected 3, got '$(pk_cn pkdup)'"
[ "$(pk_cn pkflat)" = "1" ] \
  && ok "PEAK reports 1 copy for the haplotype that does not duplicate" \
  || bad "PEAK CN for the non-duplicating haplotype: expected 1, got '$(pk_cn pkflat)'"
# it must still be labelled heuristic: nothing here validates it against external truth
grep -q "CN_CONFIDENCE=HEURISTIC" "$PKV" \
  && ok "the PEAK record still declares itself heuristic" \
  || bad "a PEAK record lost its CN_CONFIDENCE=HEURISTIC label"

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
  printf 'P\tcref#0#chr9:1-570\tc1+,cR+,c2+,cR+,c3+,d1+,ck+,d2+,c1x+\t*\n'
  printf 'P\tchapref#1#chr9\tc1+,cR+,c2+,cR+,c3+,d1+,ck+,d2+,c1x+\t*\n'
  printf 'P\tchapins#1#chr9\tc1+,cR+,i1+,c2+,cR+,i2+,c3+,d1+,ck+,d2+,c1x+\t*\n'
  printf 'P\tchapdel#1#chr9\tc1+,cR+,c2+,cR+,c3+,ck+,c1x+\t*\n'
  # a haplotype that removes the whole span, so c3..c1x is ONE snarl containing d1, ck and d2 --
  # otherwise ck is a cut vertex, the two deletions land in different bubbles, and within-haplotype
  # coalescing (which is per bubble) never gets to see them
  printf 'P\tchapwipe#1#chr9\tc1+,cR+,c2+,cR+,c3+,c1x+\t*\n'; } > "$OUT/coal.gfa"
"$BIN" bubble -i "$OUT/coal.gfa" -r 'cref#0#chr9:1-570' -o "$OUT/cob" -q >/dev/null 2>&1
"$BIN" call -i "$OUT/cob.sorted.gfa" -c "$OUT/cob.bubbles.csv" -r 'cref#0#chr9:1-570' \
      -o "$OUT/coc" -q >/dev/null 2>&1
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
  fi
fi


# ---------------------------------------------------------------- repeated module boundary
# CN, CNBP, CN_MODULE_REF_BP, POS and END share one occurrence-aware resolver, which deliberately
# selects the WIDEST span (first source to last sink) when a boundary recurs. No reference locus
# produces a CN record with a repeated boundary, so nothing exercised that choice.
#
# The reference here visits S and K TWICE, with a full module between each pair. Nothing is
# self-adjacent, so the REP route cannot fire and the collapsed-module route is used.
{ printf 'H\tVN:Z:1.0\n'
  printf 'S\tF\t%s\nS\tS\t%s\nS\tP\t%s\nS\tQ\t%s\nS\tK\t%s\nS\tG\t%s\n' \
    "$(gblk 150 21)" "$(gblk 50 22)" "$(gblk 150 23)" "$(gblk 150 24)" "$(gblk 50 25)" "$(gblk 150 26)"
  for e in "F S" "S P" "P Q" "Q K" "K S" "K G" "Q P" "P K"; do
    set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"; done
  printf 'P\tgrch38#0#chr1:1-3000\tF+,S+,P+,Q+,K+,S+,P+,Q+,K+,G+\t*\n'
  printf 'P\th1#1#chr1\tF+,S+,P+,Q+,P+,Q+,K+,S+,P+,Q+,K+,G+\t*\n'
  printf 'P\th2#1#chr1\tF+,S+,P+,K+,S+,P+,Q+,K+,G+\t*\n'; } > "$OUT/modrep.gfa"
printf 'bubble_id,source,source_orient,sink,sink_orient,inside_node_count,total_node_count,path_support,distinct_alleles,ref_allele_support,alt_allele_support_max,alt_allele_support_min,min_inside_bp,max_inside_bp,inside_nodes\n1,S,+,K,+,2,4,3,3,1,1,1,150,450,"P;Q"\n' > "$OUT/modrep.bubbles.csv"
"$BIN" call -i "$OUT/modrep.gfa" -c "$OUT/modrep.bubbles.csv" \
   --reference-path "grch38#0#chr1:1-3000" -o "$OUT/modrep" --cn --min-sv-bp 20 -q >/dev/null 2>&1
MRV="$OUT/modrep.region.vcf"
inf() { awk -F'\t' -v k="$1" '$0!~/^#/ && $8~/SVTYPE=DUP/{n=split($8,a,";");
          for(i=1;i<=n;i++){split(a[i],kv,"="); if(kv[1]==k){print kv[2]; exit}}}' "$MRV"; }
fmtv() { awk -F'\t' -v s="$1" -v k="$2" '$0~/^#CHROM/{for(i=10;i<=NF;i++) if($i==s) col=i}
          $0!~/^#/ && $8~/SVTYPE=DUP/{split($9,f,":"); for(j in f) if(f[j]==k) ki=j
            split($col,v,":"); print v[ki]; exit}' "$MRV"; }
[ "$(inf CN_METHOD)" = "MODULE_BP" ] \
  && ok "a repeated module boundary still routes to MODULE_BP" \
  || bad "CN_METHOD=$(inf CN_METHOD), expected MODULE_BP"
[ "$(inf CN_SPAN_AMBIGUOUS)" = "3" ] \
  && ok "CN_SPAN_AMBIGUOUS counts every traverser with a repeated boundary (3)" \
  || bad "CN_SPAN_AMBIGUOUS=$(inf CN_SPAN_AMBIGUOUS), expected 3"
[ "$(inf CN_MODULE_REF_BP)" = "600" ] \
  && ok "CN_MODULE_REF_BP sums the module nodes over both visits (600)" \
  || bad "CN_MODULE_REF_BP=$(inf CN_MODULE_REF_BP), expected 600"
[ "$(fmtv h1#1#chr1 CNBP)" = "300" ] && [ "$(fmtv h2#1#chr1 CNBP)" = "-150" ] \
  && ok "CNBP is measured over the same node set (+300 / -150)" \
  || bad "CNBP h1=$(fmtv h1#1#chr1 CNBP) h2=$(fmtv h2#1#chr1 CNBP), expected 300 / -150"
# POS/END are the widest SPAN, which with a repeated boundary encloses the intervening K and S. So
# END-POS EXCEEDS CN_MODULE_REF_BP by exactly that intervening sequence. The identity
# END-POS == CN_MODULE_REF_BP holds only when every boundary is visited once, and CN_SPAN_AMBIGUOUS
# is precisely the flag that says it does not.
mr_pos=$(awk -F'\t' '$0!~/^#/ && $8~/SVTYPE=DUP/{print $2; exit}' "$MRV")
mr_span=$(( $(inf END) - mr_pos ))
[ "$mr_span" = "700" ] \
  && ok "END-POS is the widest span (700), not the node sum" \
  || bad "END-POS=$mr_span, expected 700"
[ "$(( mr_span - $(inf CN_MODULE_REF_BP) ))" = "100" ] \
  && ok "and the excess is exactly the boundary visits between the outermost ones (K+S=100)" \
  || bad "span excess $(( mr_span - $(inf CN_MODULE_REF_BP) )), expected 100"

printf "\n"
if [ "$fails" -eq 0 ]; then printf "call_stats: all assertions passed\n"; exit 0; fi
printf "call_stats: %d assertion(s) failed\n" "$fails"; exit 1
