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

printf "\n"
if [ "$fails" -eq 0 ]; then printf "call_stats: all assertions passed\n"; exit 0; fi
printf "call_stats: %d assertion(s) failed\n" "$fails"; exit 1
