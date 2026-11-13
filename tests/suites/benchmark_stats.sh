#!/usr/bin/env bash
# benchmark_stats.sh - contract assertions for `panvar benchmark`.
#
# benchmark applies to the reference the events call assigned to each haplotype (each carrier's own,
# from <call>.carrier_events.tsv), aligns the result to the haplotype, and lists the large differences
# left with a cause. Each fixture below has a hand-countable answer: what the haplotypes differ by, and
# which of those call accounts for. The error cases edit call's output consistently -- the VCF and the
# carrier events together -- to stand for a call that went wrong.
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
# bm <benchmark args>: benchmark with its log kept at <-o prefix>.log (-q dropped), for the counts it logs.
bm() {
  local a=() o="" prev=""
  for x in "$@"; do [ "$prev" = "-o" ] && o=$x; [ "$x" = "-q" ] || a+=("$x"); prev=$x; done
  "$BIN" benchmark "${a[@]}" 2> "$o.log"
}
# lg <-o name> <sed pattern capturing one number>: that number from the run's log, 0 when the line is absent.
lg() { local v; v=$(sed -n "s/$2/\1/p" "$OUT/$1.log" | head -1); echo "${v:-0}"; }
APPLIED_CALL='.*events applied: call: \([0-9]*\);.*'
CONFLICTS_CALL='.*overlapping another: call: \([0-9]*\);.*'
INSIDE_CN='.*own copy stretch: \([0-9]*\);.*'
OUTSIDE='.*past the scored stretch: \([0-9]*\)$'
UNRESOLVED='.*warning: \([0-9]*\) carrier event(s) call could not place.*'
UNSPELLED='.*warning: \([0-9]*\) VCF record(s) that could not be spelled.*'
UNCHECKED='.*warning: \([0-9]*\) carrier check(s) skipped.*'
R='ref#0#chr1'

# Fixture graphs. Every node sequence is random (seeded by the fixture name), so a difference is only
# where a walk puts one.
cat > "$OUT/gen.py" <<'PYEOF'
import sys, random
name, out = sys.argv[1], sys.argv[2]
rng = random.Random(sum(map(ord, name)))
def seq(n): return ''.join(rng.choice('ACGT') for _ in range(n))
def snp(s, k=1):
    s = list(s)
    for i in rng.sample(range(5, len(s) - 5), k):
        s[i] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[s[i]]
    return ''.join(s)
def rc(s): return s[::-1].translate(str.maketrans('ACGT', 'TGCA'))
def write(nodes, paths, fn):
    links = set()
    for _, steps in paths:
        for a, b in zip(steps, steps[1:]):
            links.add((a[:-1], a[-1], b[:-1], b[-1]))
    with open(fn, 'w') as f:
        f.write('H\tVN:Z:1.0\n')
        for k, v in nodes.items(): f.write(f'S\t{k}\t{v}\n')
        for a, ao, b, bo in sorted(links): f.write(f'L\t{a}\t{ao}\t{b}\t{bo}\t0M\n')
        for pn, steps in paths:
            if pn == REF:   # the reference's interval, 1-based inclusive, ends where its walk does
                pn = 'ref#0#chr1:1001-%d' % (1000 + sum(len(nodes[st[:-1]]) for st in steps))
            f.write(f'P\t{pn}\t{",".join(steps)}\t*\n')
REF = 'ref#0#chr1:1001-2000'
if name in ('snpswap', 'snponly'):
    # one bubble: ref f1 A f2; hapsnp takes A' (A with 1 SNP, 60 bp); hapdup takes A2 (same sequence as A,
    # another node); hapdel deletes A (a real 60 bp deletion, so the bubble has a call)
    A = seq(60)
    n = {'f1': seq(100), 'A': A, 'Asnp': snp(A), 'Adup': A, 'f2': seq(100)}
    p = [(REF, ['f1+', 'A+', 'f2+']), ('hapref#1#c', ['f1+', 'A+', 'f2+']),
         ('hapsnp#1#c', ['f1+', 'Asnp+', 'f2+']), ('hapdup#1#c', ['f1+', 'Adup+', 'f2+']),
         ('hapdel#1#c', ['f1+', 'f2+'])]
    if name == 'snponly': p = p[:4]
elif name == 'frag':
    # a 75 bp deletion the graph holds as three 25 bp nodes separated by 5 bp shared nodes; hapx replaces
    # the whole interior so all of it is one bubble.
    n = {'s': seq(100), 'd1': seq(25), 'm1': seq(5), 'd2': seq(25), 'm2': seq(5), 'd3': seq(25), 'x': seq(30), 't': seq(100)}
    full = ['s+', 'd1+', 'm1+', 'd2+', 'm2+', 'd3+', 't+']
    p = [(REF, full), ('hapref#1#c', full), ('hapfrag#1#c', ['s+', 'm1+', 'm2+', 't+']),
         ('hapx#1#c', ['s+', 'x+', 't+'])]
elif name in ('trace', 'trace1'):
    # worked-trace bubble 1: ref s,n2(60),n3,n4(20),t ; H: s,n3,n4b(20),t
    n = {'s': seq(100), 'n2': seq(60), 'n3': seq(40), 'n4': seq(20), 'n4b': seq(20), 't': seq(100)}
    p = [(REF, ['s+', 'n2+', 'n3+', 'n4+', 't+']), ('H#1#c', ['s+', 'n3+', 'n4b+', 't+']),
         ('hapref#1#c', ['s+', 'n2+', 'n3+', 'n4+', 't+'])]
    if name == 'trace1':
        n['x'] = seq(30); p.append(('hapx#1#c', ['s+', 'x+', 't+']))
elif name in ('strand', 'strandrev', 'refrev'):
    # forward reference, a 60 bp deletion and a 70 bp insertion; 'strandrev' writes the haplotypes
    # (and optionally) as reverse-complement walks
    n = {'f1': seq(100), 'D': seq(60), 'f2': seq(100), 'I': seq(70), 'f3': seq(100)}
    ref = ['f1+', 'D+', 'f2+', 'f3+']
    hd = ['f1+', 'f2+', 'f3+']; hi = ['f1+', 'D+', 'f2+', 'I+', 'f3+']
    def rev(st): return [s[:-1] + ('-' if s[-1] == '+' else '+') for s in reversed(st)]
    if name == 'strandrev': hd, hi = rev(hd), rev(hi)
    p = [(REF, ref), ('hapref#1#c', ref), ('hapdel#1#c', hd), ('hapins#1#c', hi)]
    if name == 'refrev': p[0] = (REF, rev(ref))
elif name == 'shift':
    # hap H spells exactly the reference sequence with a different node partition around a shared node m
    m = seq(10); A0 = seq(50); B = seq(20)
    n = {'s': seq(100), 'A': A0 + m, 'm': m, 'B': B, 'A1': A0, 'B1': m + B, 'x': seq(30), 't': seq(100)}
    p = [(REF, ['s+', 'A+', 'm+', 'B+', 't+']), ('hapref#1#c', ['s+', 'A+', 'm+', 'B+', 't+']),
         ('H#1#c', ['s+', 'A1+', 'm+', 'B1+', 't+']), ('hapdel#1#c', ['s+', 'm+', 'B+', 't+']),
         ('hapx#1#c', ['s+', 'x+', 't+'])]
elif name == 'twice':
    # hapdup crosses bubble (s,t) twice: s a t u s b t u2 ; ref s a t u u2... second visit takes b
    n = {'s': seq(100), 'a': seq(30), 'b': seq(30), 't': seq(100), 'u': seq(200), 'u2': seq(100)}
    p = [(REF, ['s+', 'a+', 't+', 'u+', 'u2+']), ('hapref#1#c', ['s+', 'a+', 't+', 'u+', 'u2+']),
         ('hapb#1#c', ['s+', 'b+', 't+', 'u+', 'u2+']),
         ('hapdup#1#c', ['s+', 'a+', 't+', 'u+', 's+', 'b+', 't+', 'u+', 'u2+'])]
elif name == 'rep':
    # a 60 bp tandem unit U (self loop): the reference holds 2 copies, hap4 4, hap1 1, hap2 2
    U = seq(60)
    n = {'s': seq(100), 'U': U, 't': seq(100)}
    p = [(REF, ['s+', 'U+', 'U+', 't+']), ('hap2#1#c', ['s+', 'U+', 'U+', 't+']),
         ('hap4#1#c', ['s+', 'U+', 'U+', 'U+', 'U+', 't+']), ('hap4b#1#c', ['s+', 'U+', 'U+', 'U+', 'U+', 't+']),
         ('hap1#1#c', ['s+', 'U+', 't+'])]
elif name == 'motif':
    # 'rep' with a 600 bp unit holding one 40 bp motif twice, 480 bp apart: k-mers recur near 480 as well
    M = seq(40)
    U = seq(20) + M + seq(440) + M + seq(60)
    n = {'s': seq(100), 'U': U, 't': seq(100)}
    p = [(REF, ['s+', 'U+', 'U+', 't+']), ('hap2#1#c', ['s+', 'U+', 'U+', 't+']),
         ('hap4#1#c', ['s+', 'U+', 'U+', 'U+', 'U+', 't+']), ('hap1#1#c', ['s+', 'U+', 't+'])]
elif name == 'partial':
    # the strand fixture's deletion, with one haplotype starting after the first flank
    n = {'f0': seq(100), 'f1': seq(100), 'D': seq(60), 'f2': seq(100), 'f3': seq(100)}
    ref = ['f0+', 'f1+', 'D+', 'f2+', 'f3+']
    p = [(REF, ref), ('hapref#1#c', ref), ('hapdel#1#c', ['f0+', 'f1+', 'f2+', 'f3+']),
         ('hapcut#1#c', ['f1+', 'f2+', 'f3+'])]
write(n, p, f'{out}/{name}.gfa')
open(f'{out}/{name}.seqs', 'w').write('\n'.join(f'{k}\t{v}' for k, v in n.items()) + '\n')
PYEOF

# pipeline <fixture> <label> [bubble options] -- [call options]: bubble, call, benchmark
pipeline() {
  local fx=$1 lb=$2; shift 2
  local bopt=() copt=() into=b
  for a in "$@"; do if [ "$a" = "--" ]; then into=c; elif [ $into = b ]; then bopt+=("$a"); else copt+=("$a"); fi; done
  "${PYTHON:-python3}" -I "$OUT/gen.py" "$fx" "$OUT" >/dev/null
  "$BIN" bubble -i "$OUT/$fx.gfa" -r "$R" -o "$OUT/${lb}_b" --merge-nearby-bp 0 --min-variant-bp 1 -q ${bopt[@]+"${bopt[@]}"} >/dev/null 2>&1
  "$BIN" call -i "$OUT/${lb}_b.sorted.gfa" -c "$OUT/${lb}_b.bubbles.csv" -r "$R" -o "$OUT/${lb}_c" \
         --no-per-bubble-vcf -q ${copt[@]+"${copt[@]}"} >/dev/null 2>&1
  bench "$lb"
}
bench() { bm -i "$OUT/$1_b.sorted.gfa" -c "$OUT/$1_b.bubbles.csv" -r "$R" --calls "$OUT/$1_c" \
                 -o "$OUT/$1_B" -q >/dev/null 2>&1; }
val() { awk -F'\t' -v k="$2" '$1==k{print $2}' "$OUT/$1_B.summary.tsv"; }
# present:left:left_bp for a label
pl() { echo "$(val "$1" large_present):$(val "$1" large_left_call):$(val "$1" large_left_call_bp)"; }
# call's leftovers: type size reason bubbles records
left_rows() { tail -n +2 "$OUT/$1_B.events.tsv" | awk -F'\t' '$2=="call"{print $5" "$6" "$9" "$10" "$11}' | paste -sd';' -; }
# the VCF's leftovers with a VCF reason (lost in the VCF): haplotype type size reason
lost_rows() { tail -n +2 "$OUT/$1_B.events.tsv" | awk -F'\t' '$2=="vcf" && $9 ~ /^(no_sequence|copy_count|representative|other)$/{print substr($1,1,2)" "$5" "$6" "$9}' | paste -sd';' -; }

# ---- what call accounts for leaves nothing --------------------------------------------------------
# snpswap: a 60 bp deletion (hapdel), a 1 bp SNP in a 60 bp node (hapsnp), and another node spelling the
# reference's own 60 bp (hapdup). One large event; the SNP is small and the identical node is nothing.
pipeline snpswap snpswap
check "snpswap: large events present:left:left_bp" "$(pl snpswap)" "1:0:0"
check "snpswap: the SNP is the only small difference, before and after" "$(val snpswap snv_present):$(val snpswap snv_left_call)" "1:1"
# frag: a 75 bp deletion the graph holds as three 25 bp nodes, and a 30-for-80 bp replacement
pipeline frag frag
check "frag: a deletion split over three nodes is one event, and is recovered" "$(pl frag)" "2:0:0"
# strand / strandrev: a 60 bp deletion and a 70 bp insertion, the haplotypes walked forward or reversed
pipeline strand strand
pipeline strandrev strandrev
check "strand: deletion and insertion recovered" "$(pl strand)" "2:0:0"
check "strandrev: the same with the haplotypes walked on the other strand" "$(pl strandrev)" "2:0:0"
# twice: a haplotype crossing the bubble twice, with a node between the crossings
pipeline twice twice
check "twice: a haplotype crossing a bubble twice is scored and recovered" "$(val twice haplotypes_scored):$(pl twice)" "3:1:0:0"
# shift: one haplotype spells the reference with its nodes cut differently
pipeline shift shift
check "shift: a different cut of the same sequence is no difference" "$(pl shift)" "2:0:0"
# rep: a 60 bp tandem unit at 2 copies on the reference, 4 on two haplotypes and 1 on another
pipeline rep rep
check "rep: copy-number changes recovered (+120, +120, -60)" "$(pl rep):$(val rep large_present_bp)" "3:0:0:300"
check "rep: each carrier's CN matches the copies its walk holds" "$(val rep carriers_checked):$(val rep carrier_mismatches)" "3:0"
# partial: a haplotype starting one flank late is cut to what it covers, not charged for the rest
pipeline partial partial
check "partial: a haplotype covering part of the region is scored over that part" \
      "$(awk -F'\t' '$1=="hapcut#1#c"{print $2":"$3":"$4}' "$OUT/partial_B.haplotypes.tsv")" "scored:1101:1460"
check "partial: and nothing is left" "$(pl partial)" "2:0:0"

# ---- what call does not account for, by cause --------------------------------------------------------
pipeline strand norec -- --min-sv-bp 65
check "a 60 bp deletion call did not report (call at 65 bp) is not_reported, in its site" "$(left_rows norec)" "DEL 60 not_reported 1 ."
pipeline strand outside --min-variant-bp 65
check "one in no bubble at all is not_reported, in no site" "$(left_rows outside)" "DEL 60 not_reported . ."

# A missed carrier: hapdel's genotype for the deletion is 0, in the VCF and the carrier events alike.
pipeline strand missed
del_id=$(awk -F'\t' '!/^#/ && $8 ~ /SVTYPE=DEL/{print $3}' "$OUT/missed_c.region.vcf")
awk -F'\t' -v OFS='\t' -v id="$del_id" '/^#CHROM/{for(i=10;i<=NF;i++) if($i=="hapdel#1#c") c=i} !/^#/ && $3==id{split($c,f,":"); f[1]="0"; $c=f[1]; for(j=2;j in f;j++) $c=$c":"f[j]} {print}' \
  "$OUT/missed_c.region.vcf" > "$OUT/x.vcf" && mv "$OUT/x.vcf" "$OUT/missed_c.region.vcf"
awk -F'\t' -v id="$del_id" '$1!=id' "$OUT/missed_c.carrier_events.tsv" > "$OUT/x.tsv" && mv "$OUT/x.tsv" "$OUT/missed_c.carrier_events.tsv"
bench missed
check "a carrier call missed is not_carrier, naming the record" "$(left_rows missed | cut -d' ' -f1-3,5)" "DEL 60 not_carrier $del_id"
vcf_rows() { tail -n +2 "$OUT/$1_B.events.tsv" | awk -F'\t' '$2=="vcf"{print $5" "$6" "$9}' | paste -sd';' -; }
check "the VCF's rebuild leaves that miss again, whole, with its cause, not as its own (vcf:vcf_only)" \
      "$(vcf_rows missed):$(val missed large_left_vcf):$(val missed large_left_vcf_only)" "DEL 60 not_carrier:1:0"

# An event call got wrong: hapins's own insertion cut from 70 bp to 10. What is left is inside the event
# applied to it (incomplete), and its own event no longer matches the record's SVLEN.
pipeline strand cut
awk -F'\t' -v OFS='\t' '$3=="seq" || $3=="ins" { if ($9 == 70) { $8 = 10; $9 = 10; $10 = substr($10, 1, 10) } } {print}' \
  "$OUT/cut_c.carrier_events.tsv" > "$OUT/x.tsv" && mv "$OUT/x.tsv" "$OUT/cut_c.carrier_events.tsv"
bench cut
check "an event shorter than the haplotype's leaves the rest as incomplete" "$(left_rows cut | cut -d' ' -f1-3)" "INS 60 incomplete"
check "and its carrier is flagged against the record's SVLEN (own 10 bp, SVLEN 70)" \
      "$(tail -n +2 "$OUT/cut_B.carrier_checks.tsv" | cut -f2,6)" "$(printf 'hapins#1#c\t-60')"
check "the VCF's record spells the whole insertion: its rebuild leaves nothing (call:vcf)" \
      "$(val cut large_left_call):$(val cut large_left_vcf):$(vcf_rows cut)" "1:0:"

# A wrong copy number: hap4 reported at CN 3 while its walk holds 4 copies.
pipeline rep wrongcn
awk -F'\t' -v OFS='\t' '/^#CHROM/{for(i=10;i<=NF;i++) if($i=="hap4#1#c") c=i}
  !/^#/ && $8 ~ /SVTYPE=CNV/ { n=split($9,fm,":"); for(j=1;j<=n;j++) if(fm[j]=="CN") k=j; split($c,f,":"); f[k]="3"; $c=f[1]; for(j=2;j<=n;j++) $c=$c":"f[j] } {print}' \
  "$OUT/wrongcn_c.region.vcf" > "$OUT/x.vcf" && mv "$OUT/x.vcf" "$OUT/wrongcn_c.region.vcf"
bench wrongcn
check "a CN that disagrees with the copies in the walk is flagged (CN 3, walk 4)" \
      "$(tail -n +2 "$OUT/wrongcn_B.carrier_checks.tsv" | awk -F'\t' '$3=="cn"' | cut -f2,4,5)" "$(printf 'hap4#1#c\tCN=3\t4 copies in its walk')"
check "and its copy stretch is 60 bp longer than CN 3 allows" \
      "$(tail -n +2 "$OUT/wrongcn_B.carrier_checks.tsv" | awk -F'\t' '$3=="cn_stretch"' | cut -f2,6)" "$(printf 'hap4#1#c\t60')"

# A copy is the reference's period, not the record's unit: here the unit is given as 480 bp (a share of the
# 600 bp copy), where a motif in each copy also recurs. Copies of the 600 bp period rebuild every haplotype.
pipeline motif motif
for f in "$OUT"/motif_b.* "$OUT"/motif_c.*; do cp "$f" "${f/motif_/motifs_}"; done
sed -i.bak 's/RU_LEN=600/RU_LEN=480/' "$OUT/motifs_c.region.vcf" && rm -f "$OUT/motifs_c.region.vcf.bak"
bench motifs
check "a unit shorter than the copies still rebuilds them from the period (left by the VCF only)" \
      "$(val motif large_left_vcf_only):$(val motifs large_left_vcf_only)" "0:0"

# ---- copy-number stretches, and the VCF mode (spec fixtures: tests/tools/spec_gfa.py) --------------------
SPEC_GFA="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/tools/spec_gfa.py"
# spec <name> [call options]: build the fixture from $OUT/<name>.spec (plus <name>.ug.tsv), call, benchmark
spec() {
  local n=$1; shift
  "${PYTHON:-python3}" -I "$SPEC_GFA" "$OUT/$n.spec" "$OUT/$n" >/dev/null
  [ -f "$OUT/$n.ug.tsv" ] && cp "$OUT/$n.ug.tsv" "$OUT/$n/g.panphorte.unit_graph.tsv"
  local r; r=$(awk '$1=="P"{print $2; exit}' "$OUT/$n.spec")
  if [ -n "${BUBBLE:-}" ]; then   # the sites from panvar bubble (and its sorted graph), not the spec's
    "$BIN" bubble -i "$OUT/$n/g.gfa" -r "$r" -o "$OUT/$n/b" -q >/dev/null 2>&1
    mv "$OUT/$n/b.sorted.gfa" "$OUT/$n/g.gfa"; mv "$OUT/$n/b.bubbles.csv" "$OUT/$n/g.bubbles.csv"
  fi
  "$BIN" call -i "$OUT/$n/g.gfa" -c "$OUT/$n/g.bubbles.csv" -r "$r" -o "$OUT/${n}_c" --no-per-bubble-vcf -q "$@" >/dev/null 2>&1
  bm -i "$OUT/$n/g.gfa" -c "$OUT/$n/g.bubbles.csv" -r "$r" --calls "$OUT/${n}_c" -o "$OUT/${n}_B" -q >/dev/null 2>&1
}
sv() { awk -F'\t' -v k="$2" '$1==k{print $2}' "$OUT/$1.summary.tsv"; }

# noforms: a VCF without the copies its module records name (CN_FORMS, CNFORM), as an older call wrote it, so
# benchmark finds each carrier's copy among the reference's
noforms() { awk -F'\t' -v OFS='\t' '/^#/ || $9 !~ /:CNFORM$/ {print; next}
  { sub(/;CN_FORMS=[^;]*/, "", $8); n = split($9, k, ":"); sub(/:CNFORM$/, "", $9)
    for (i = 10; i <= NF; i++) if (split($i, v, ":") == n) sub(/:[^:]*$/, "", $i); print }'; }
# vleft <name>: large events left by the VCF's records
vleft() { sv "$1_B" large_left_vcf; }
# A 60 bp repeat unit U, two copies on the reference; hX holds an 80 bp X between its first and second
# copy. Its copy stretch holds X, and call's insertion of X is inside it: applied twice, it left 80 bp.
cat > "$OUT/ubetween.spec" <<'EOF'
N S 300 1
N U 60 2
N X 80 3
N T 300 4
P ref#0#chr1:0-720 S+,U+,U+,T+
P h1#1#chr1 S+,U+,U+,T+
P hX#1#chr1 S+,U+,X+,U+,U+,T+
P h4#1#chr1 S+,U+,U+,U+,U+,T+
B S T U;X
EOF
spec ubetween
check "an insertion between a carrier's copies is part of its copy stretch, not applied twice (left:inside:mismatches)" \
      "$(sv ubetween_B large_left_call):$(lg ubetween_B "$INSIDE_CN"):$(sv ubetween_B carrier_mismatches)" "0:1:0"
# An insertion that repeats its neighbours can go in anywhere over its CIPOS, and its listed place can fall
# on the wrong side of a stretch's edge; the stretch's length says which side it is on. ubedge: hX's
# stretch holds X, listed past the stretch. uedge: X goes in before the copies (an aligned allele's stretch
# is its copies alone), listed at the stretch's first base.
sed 's/^P hX.*/P hX#1#chr1 S+,X+,U+,U+,U+,T+/' "$OUT/ubetween.spec" > "$OUT/uedge.spec"
spec uedge
awk -F'\t' 'BEGIN{OFS="\t"} $1=="bubble1_INS_X"{$4=422; $5=421} 1' "$OUT/ubetween_c.carrier_events.tsv" > "$OUT/ubedge_c.carrier_events.tsv"
sed 's/SVTYPE=INS;SVLEN=80;/SVTYPE=INS;SVLEN=80;CIPOS=0,121;/' "$OUT/ubetween_c.region.vcf" > "$OUT/ubedge_c.region.vcf"
sed 's/SVTYPE=INS;SVLEN=80;/SVTYPE=INS;SVLEN=80;CIPOS=-5,0;/' "$OUT/uedge_c.region.vcf" > "$OUT/x.vcf" && mv "$OUT/x.vcf" "$OUT/uedge_c.region.vcf"
for n in ubetween:ubedge uedge:uedge; do
  bm -i "$OUT/${n%%:*}/g.gfa" -c "$OUT/${n%%:*}/g.bubbles.csv" -r 'ref#0#chr1:0-720' --calls "$OUT/${n#*:}_c" \
         -o "$OUT/${n#*:}_B" -q >/dev/null 2>&1
done
check "an ambiguous insertion listed past a stretch that holds it is inside (left:inside:mismatches)" \
      "$(sv ubedge_B large_left_call):$(lg ubedge_B "$INSIDE_CN"):$(sv ubedge_B carrier_mismatches)" "0:1:0"
check "an ambiguous insertion listed at a stretch that lacks it is applied (left:inside:mismatches)" \
      "$(sv uedge_B large_left_call):$(lg uedge_B "$INSIDE_CN"):$(sv uedge_B carrier_mismatches)" "0:0:0"
# A unit U, three copies on the reference. h1 deletes the last 3 bp of A (node P), keeps the 8 bp node B and
# inserts a 180 bp X before its one copy. Aligned, the deletion and X are one replacement starting in A,
# before the copies; its copy stretch is its copies alone, so X is applied once (uflank). In uedge2, h1 keeps
# A whole: X goes in at the stretch's first base, and the stretch's length says it is beside the copies.
cat > "$OUT/uflank.spec" <<'EOF'
N L 300 101
N A =CGCCTGTACTTCCAACTCCTTTACAGTCTTTCCGCACGTCTGGTTCTCCCGATATCGATGTCCGTCCCTCCTCAACTGCAACATGGGTTCTCGCATCTCC
N P =CGCCTGTACTTCCAACTCCTTTACAGTCTTTCCGCACGTCTGGTTCTCCCGATATCGATGTCCGTCCCTCCTCAACTGCAACATGGGTTCTCGCATC
N B 8 106
N U 52 103
N X 180 104
N R 300 105
P ref#0#chr1:0-864 L+,A+,B+,U+,U+,U+,R+
P h1#1#chr1 L+,P+,B+,X+,U+,R+
P h2#1#chr1 L+,A+,B+,U+,U+,U+,R+
P h3#1#chr1 L+,A+,B+,U+,U+,U+,U+,R+
B L R A;P;B;U;X
EOF
spec uflank
sed 's/^P h1#1#chr1 L+,P+,/P h1#1#chr1 L+,A+,/' "$OUT/uflank.spec" > "$OUT/uedge2.spec"
spec uedge2
for n in uflank uedge2; do
  check "$n: an aligned insertion beside the copies is applied once (left:stretch bp of h1)" \
        "$(sv ${n}_B large_left_call):$(awk -F'\t' '$3=="cn" && $11 ~ /h1#/{print $9}' "$OUT/${n}_c.carrier_events.tsv")" "0:52"
done
# A 120 bp spacer Z between the copies on one walk only: h1 holds its one copy before Z (sp2) or Z between its
# copies (sp3), or no copy (sp4). The aligned walks, copies left out, hold Z alike, so no event of h1's has it:
# its copy stretch must, running between the steps the two walks share, or Z is lost or applied twice.
cat > "$OUT/sp2.spec" <<'EOF'
N L 300 101
N A 100 102
N U 52 103
N Z 120 107
N R 300 105
P ref#0#chr1:0-976 L+,A+,U+,U+,Z+,U+,R+
P h1#1#chr1 L+,A+,U+,Z+,R+
P h2#1#chr1 L+,A+,U+,U+,Z+,U+,R+
P h3#1#chr1 L+,A+,U+,U+,U+,Z+,U+,R+
B L R A;U;Z
EOF
spec sp2
sed -e 's/^P ref#0#chr1:0-976 .*/P ref#0#chr1:0-976 L+,A+,Z+,U+,U+,U+,R+/' -e 's/^P h1#1#chr1 .*/P h1#1#chr1 L+,A+,U+,Z+,U+,R+/' \
    -e 's/^P h2#1#chr1 .*/P h2#1#chr1 L+,A+,Z+,U+,U+,U+,R+/' -e 's/^P h3#1#chr1 .*/P h3#1#chr1 L+,A+,Z+,U+,U+,U+,U+,R+/' \
    "$OUT/sp2.spec" > "$OUT/sp3.spec"
spec sp3
sed 's/^P h1#1#chr1 .*/P h1#1#chr1 L+,A+,Z+,R+/' "$OUT/sp2.spec" > "$OUT/sp4.spec"
spec sp4
# dd: h1 deletes A's last 3 bp, keeps B, inserts X before its first copy (an aligned event, reported), and holds
# Y, Z with one base changed, between its two copies: Y pairs with the reference's Z by place, so the stretch
# stays its copies and Y, and X is applied once. ddsame: the same with Z itself.
cat > "$OUT/dd.spec" <<'EOF'
N L 300 101
N A =CGCCTGTACTTCCAACTCCTTTACAGTCTTTCCGCACGTCTGGTTCTCCCGATATCGATGTCCGTCCCTCCTCAACTGCAACATGGGTTCTCGCATCTCC
N P =CGCCTGTACTTCCAACTCCTTTACAGTCTTTCCGCACGTCTGGTTCTCCCGATATCGATGTCCGTCCCTCCTCAACTGCAACATGGGTTCTCGCATC
N B 8 106
N U 52 103
N Z 120 107
N Y @Z:1:5
N X 180 104
N R 300 105
P ref#0#chr1:0-1036 L+,A+,B+,U+,Z+,U+,U+,R+
P h1#1#chr1 L+,P+,B+,X+,U+,Y+,U+,R+
P h2#1#chr1 L+,A+,B+,U+,Z+,U+,U+,R+
P h3#1#chr1 L+,A+,B+,U+,Z+,U+,U+,U+,R+
B L R A;P;B;U;X;Z;Y
EOF
spec dd
sed 's/^P h1#1#chr1 .*/P h1#1#chr1 L+,P+,B+,X+,U+,Z+,U+,R+/' "$OUT/dd.spec" > "$OUT/ddsame.spec"
spec ddsame
for n in sp2:172 sp3:224 sp4:120 dd:224 ddsame:224; do
  check "${n%%:*}: a spacer between copies on one walk is kept once (left:stretch bp of h1)" \
        "$(sv ${n%%:*}_B large_left_call):$(awk -F'\t' '$3=="cn" && $11 ~ /h1#/{print $9}' "$OUT/${n%%:*}_c.carrier_events.tsv")" "0:${n#*:}"
done
# An 80 bp X beside hX's copies, after them (eafter, three copies; eone, one) or before them (ebefore), and the
# same stored on the other strand (rebefore, reafter). Aligned, X goes in at the copies' start on the reference;
# its bases must stay on their side of hX's copies, so the stretch holds them where they follow the copies.
for n in eafter:S+,U+,U+,U+,X+,T+ eone:S+,U+,X+,T+ ebefore:S+,X+,U+,U+,U+,T+; do
  sed "s/^P hX#1#chr1 .*/P hX#1#chr1 ${n#*:}/" "$OUT/ubetween.spec" > "$OUT/${n%%:*}.spec"
  spec "${n%%:*}"
done
for n in rebefore:T-,U-,U-,U-,X-,S- reafter:T-,X-,U-,U-,U-,S-; do
  sed -e "s/^P hX#1#chr1 .*/P hX#1#chr1 ${n#*:}/" -e 's/^P ref#0#chr1:0-720 .*/P ref#0#chr1:0-720 T-,U-,U-,S-/' \
      -e 's/^P h1#1#chr1 .*/P h1#1#chr1 T-,U-,U-,S-/' -e 's/^P h4#1#chr1 .*/P h4#1#chr1 T-,U-,U-,U-,U-,S-/' \
      "$OUT/ubetween.spec" > "$OUT/${n%%:*}.spec"
  spec "${n%%:*}"
done
check "an insertion beside the copies keeps its side of them (left: eafter eone ebefore rebefore reafter)" \
      "$(sv eafter_B large_left_call) $(sv eone_B large_left_call) $(sv ebefore_B large_left_call) $(sv rebefore_B large_left_call) $(sv reafter_B large_left_call)" "0 0 0 0 0"
# No reference copy: X after the copies (hX), or before them (hY, whose copies call could not place).
cat > "$OUT/uzero.spec" <<'EOF'
N S 300 1
N U 60 2
N X 80 3
N T 300 4
P ref#0#chr1:0-600 S+,T+
P h1#1#chr1 S+,T+
P hX#1#chr1 S+,U+,U+,X+,T+
P hY#1#chr1 S+,X+,U+,U+,T+
P h2#1#chr1 S+,U+,U+,T+
B S T U;X
EOF
spec uzero
check "with no reference copy, the copies go where the walk has them, either side of an insertion (left:unresolved)" \
      "$(sv uzero_B large_left_call):$(lg uzero_B "$UNRESOLVED")" "0:0"
check "the VCF cannot spell copies with REF_CN 0: counted, and what they leave is lost (no_sequence)" \
      "$(lg uzero_B "$UNSPELLED"):$(lost_rows uzero | tr ';' '\n' | awk '$1=="h2"{print $4}')" "1:no_sequence"
check "the warning names the cause: no reference copy to paste, not a sequence call left out (lines)" \
      "$(grep -c 'no reference copy to paste' "$OUT/uzero_B.log"):$(grep -c 'max-seq-bp' "$OUT/uzero_B.log")" "1:0"
# A twin of the loop (V: U's sequence, same neighbours, no loop) and a unit graph whose copy is
# (a1|a2) b: a carrier's copies are counted as the record counts them, not on the key node alone.
cat > "$OUT/twin.spec" <<'EOF'
N S 300 1
N U 60 2
N V @U:0:1
N T 300 4
P ref#0#chr1:0-720 S+,U+,U+,T+
P h1#1#chr1 S+,U+,U+,T+
P hV#1#chr1 S+,V+,T+
P h4#1#chr1 S+,U+,U+,U+,U+,T+
P h0#1#chr1 S+,T+
B S T U;V
EOF
cat > "$OUT/ug.spec" <<'EOF'
N S 300 1
N a1 40 2
N a2 @a1:2:5
N b 60 3
N T 300 4
P ref#0#chr1:0-800 S+,a1+,b+,a1+,b+,T+
P h2#1#chr1 S+,a1+,b+,a1+,b+,T+
P h3#1#chr1 S+,a1+,b+,a2+,b+,a1+,b+,T+
P h1#1#chr1 S+,a2+,b+,T+
P h4#1#chr1 S+,a2+,b+,a2+,b+,a1+,b+,a1+,b+,T+
B S T a1;a2;b
EOF
printf 'unit_bp\tanchor_nodes\tnodes\n100\ta1;a2\ta1;a2;b\n' > "$OUT/ug.ug.tsv"
spec twin; spec ug
check "a twin's copy and a unit graph's other anchor count as copies (checked:mismatches)" \
      "$(sv twin_B carriers_checked):$(sv twin_B carrier_mismatches) $(sv ug_B carriers_checked):$(sv ug_B carrier_mismatches)" "3:0 3:0"
check "the VCF: tandem copies of the reference's copy rebuild a repeat array exactly" "$(vleft twin) $(vleft ug)" "0 0"
# --cn MODULE: A B (600 bp), one reference copy; h2, hX and hZ hold two. hX also holds an 80 bp X before
# its copies, hZ after them: each X is a record of its own, apart from the copies the module record names.
cat > "$OUT/modins.spec" <<'EOF'
N S 2000 1
N A 300 2
N B 300 3
N X 80 4
N T 2000 5
P ref#0#chr1:0-4600 S+,A+,B+,T+
P h1#1#chr1 S+,A+,B+,T+
P h2#1#chr1 S+,A+,B+,A+,B+,T+
P hX#1#chr1 S+,X+,A+,B+,A+,B+,T+
P hZ#1#chr1 S+,A+,B+,A+,B+,X+,T+
B S T A;B;X
EOF
spec modins --cn
check "MODULE carriers are checked: their copies by the anchor rule, and the two X records (checked:unchecked)" \
      "$(sv modins_B carriers_checked):$(lg modins_B "$UNCHECKED")" "5:0"
check "no check flags a carrier, its copies and its X making its stretch (rows; cn_exact checked; large_left_call)" \
      "$(tail -n +2 "$OUT/modins_B.carrier_checks.tsv" | wc -l | tr -d ' ') $(sv modins_B cn_exact_checked) $(sv modins_B large_left_call)" "0 3 0"
check "the VCF: hX's copy and its insertion at one point are both applied, and nothing is lost" \
      "$(lost_rows modins)" ""
# Two variants of call's output, kept in tests/synthetic_data/benchmark_cases:
# - modins_short_unit: CN_UNIT_BP is 300, half a copy; with one reference copy, the copy is the record's
#   span, so the copies still rebuild and only hX's 80 bp is lost.
# - modmid_missed_in_stretch: on a graph where hX's X lies between its copies, X is cut from hX's copy
#   stretch and its own record dropped, so call misses it inside the stretch. The miss is counted once: the
#   VCF stage's leftover lies on it, wherever in the stretch.
CASES="$(cd "$(dirname "$0")/.." && pwd)/synthetic_data/benchmark_cases"
bm -i "$OUT/modins/g.gfa" -c "$OUT/modins/g.bubbles.csv" -r ref --calls "$CASES/modins_short_unit" \
  -o "$OUT/modins_short_unit_B" -q >/dev/null 2>&1
check "one reference copy: the copy is the record's span, not CN_UNIT_BP (large_left_vcf_only_bp, carrier_mismatches)" \
      "$(sv modins_short_unit_B large_left_vcf_only_bp) $(sv modins_short_unit_B carrier_mismatches)" "80 1"
# With a unit only 50 bp short of the span, both are within the forms' tolerance of a carrier whose CNBP says
# 550 bp: the nearer, the unit, is pasted, and the 50 bp the haplotype really holds beyond it are lost.
sed -e 's/CN_UNIT_BP=300/CN_UNIT_BP=550/' -e 's/1:2:600:0\t1:2:600:0\t1:2:600:0$/1:2:550:0\t1:2:600:0\t1:2:600:0/' \
  "$CASES/modins_short_unit.region.vcf" > "$OUT/modins_near_unit.region.vcf"
cp "$CASES/modins_short_unit.carrier_events.tsv" "$OUT/modins_near_unit.carrier_events.tsv"
bm -i "$OUT/modins/g.gfa" -c "$OUT/modins/g.bubbles.csv" -r ref --calls "$OUT/modins_near_unit" \
  -o "$OUT/modins_near_unit_B" -q >/dev/null 2>&1
check "one reference copy, unit near the span: the form nearer the carrier's CNBP (large_left_vcf_only_bp)" \
      "$(sv modins_near_unit_B large_left_vcf_only_bp)" "130"
# A VCF with Windows line endings scores the same; an insertion runs from the base before it to the base after.
cp "$OUT/modins_c.carrier_events.tsv" "$OUT/modins_crlf.carrier_events.tsv"
sed 's/$/\r/' "$OUT/modins_c.region.vcf" > "$OUT/modins_crlf.region.vcf"
bm -i "$OUT/modins/g.gfa" -c "$OUT/modins/g.bubbles.csv" -r ref --calls "$OUT/modins_crlf" \
  -o "$OUT/modins_crlf_B" -q >/dev/null 2>&1
check "CRLF line endings: the same summary" "$(cmp -s "$OUT/modins_B.summary.tsv" "$OUT/modins_crlf_B.summary.tsv" && echo same)" "same"
check "present.tsv: a pure insertion runs from the base before it to the base after (span 1)" \
      "$(awk -F'\t' '$4=="INS"{print $3-$2}' "$OUT/modins_B.present.tsv" | sort -n | head -n 1)" "1"
# A gained copy of a shorter form: two reference copies of A H B; hS holds a third copy, before them, without
# H, and hD deletes H from its first copy, so the VCF has a record of that deletion. The copy-number record says
# hS's change is 200 bp short of a whole copy (FORMAT:CNBP), which that record accounts for: the VCF rebuild
# pastes the copy without H where it inserts copies, and nothing is lost.
cat > "$OUT/shortcopy.spec" <<'EOF'
N S 2000 1
N A 300 2
N H 200 3
N B 300 4
N T 2000 5
P ref#0#chr1:0-5600 S+,A+,H+,B+,A+,H+,B+,T+
P h1#1#chr1 S+,A+,H+,B+,A+,H+,B+,T+
P hS#1#chr1 S+,A+,B+,A+,H+,B+,A+,H+,B+,T+
P hD#1#chr1 S+,A+,B+,A+,H+,B+,T+
P hG#1#chr1 S+,A+,H+,B+,A+,H+,B+,A+,H+,B+,T+
B S T A;H;B
EOF
spec shortcopy --cn
check "a gained copy of a shorter form, which a deletion record accounts for, is rebuilt (large_left_vcf_only_bp)" \
      "$(sv shortcopy_B large_left_vcf_only_bp)" "0"
# Copies and the other events of each haplotype: the module record names the copies each carrier gained or lost
# (CN_FORMS, CNFORM), and its other records the rest, so the VCF rebuilds every haplotype. The reference runs
# A3 A2 twice, S between the copies and X1 X2 X3 after them (see call_stats.sh, the same haplotypes).
{ printf 'N L 2000 81\nN A3 60 82\nN A2 100 83\nN S 40 84\nN X1 80 85\nN X2 80 86\nN X3 240 87\nN Y 70 89\nN R 2000 88\n'
  printf 'P ref#0#chr1:0-4680 L+,A3+,A2+,S+,A3+,A2+,X1+,X2+,X3+,R+\n'
  printf 'P hgain#1#chr1 L+,A3+,A2+,S+,A3+,A2+,S+,A3+,A2+,X1+,X2+,X3+,R+\nP hloss#1#chr1 L+,A3+,A2+,X1+,X2+,X3+,R+\n'
  printf 'P hhalf#1#chr1 L+,A3+,A2+,S+,X2+,X3+,R+\nP heq#1#chr1 L+,A3+,A2+,S+,X3+,R+\nP hmore#1#chr1 L+,A3+,A2+,S+,R+\n'
  printf 'P hdir#1#chr1 L+,A3+,A2+,S+,A3+,A2+,S+,A3+,A2+,R+\nP hinv#1#chr1 L+,A3+,A2+,S+,A3+,A2+,S-,A2-,A3-,X1+,X2+,X3+,R+\n'
  printf 'P hyy#1#chr1 L+,A3+,Y+,A2+,S+,A3+,Y+,A2+,X1+,X2+,X3+,R+\n'
  for i in 1 2 3 4 5 6 7 8; do printf 'P hr%s#1#chr1 L+,A3+,A2+,S+,A3+,A2+,X1+,X2+,X3+,R+\n' "$i"; done; } > "$OUT/sepb.spec"
BUBBLE=1 spec sepb --cn
# cn_exact <name>: carriers checked by cn_exact, and those that differ
cn_exact() { echo "$(awk -F'\t' '$1=="cn_exact_checked"{print $2}' "$OUT/$1_B.summary.tsv"):$(tail -n +2 "$OUT/$1_B.carrier_checks.tsv" | awk -F'\t' '$3=="cn_exact"' | wc -l | tr -d ' ')"; }
check "copies: the VCF rebuilds every haplotype (large_left_vcf_bp; cn_exact checked:differing)" \
  "$(sv sepb_B large_left_vcf_bp) $(cn_exact sepb)" "0 6:0"
cp "$OUT/sepb.spec" "$OUT/sepbs.spec"
BUBBLE=1 spec sepbs --cn --max-seq-bp 100
check "copies, --max-seq-bp under a copy: gained copies named by the reference stretches they repeat still rebuild (large_left_vcf_bp; cn_exact)" \
  "$(sv sepbs_B large_left_vcf_bp) $(cn_exact sepbs)" "0 6:0"
awk '$1=="P"{n=split($3,s,","); o=""; for(i=n;i>=1;i--){x=s[i]; c=substr(x,length(x)); x=substr(x,1,length(x)-1) (c=="+"?"-":"+")
            o=o (o==""?"":",") x}; print $1, $2, o; next} {print}' "$OUT/sepb.spec" > "$OUT/sepbr.spec"
BUBBLE=1 spec sepbr --cn
check "copies, the reference walked backwards: the VCF rebuilds every haplotype (large_left_vcf_bp; cn_exact)" \
  "$(sv sepbr_B large_left_vcf_bp) $(cn_exact sepbr)" "0 6:0"
# No reference copy: the gained copies are the record's own sequence.
{ printf 'N L 2000 91\nN A 150 92\nN B 150 93\nN C 30 94\nN R 2000 95\nP ref#0#chr1:0-4030 L+,C+,R+\n'
  for i in 1 2 3 4; do printf 'P hone%s#1#chr1 L+,A+,B+,C+,R+\n' "$i"; done
  for i in 1 2; do printf 'P htwo%s#1#chr1 L+,A+,B+,A+,B+,C+,R+\n' "$i"; done
  for i in 1 2 3 4 5 6; do printf 'P hr%s#1#chr1 L+,C+,R+\n' "$i"; done; } > "$OUT/sepb0.spec"
BUBBLE=1 spec sepb0 --cn
check "copies, no reference copy: the VCF rebuilds the gains (large_left_vcf_bp, no_sequence; cn_exact)" \
  "$(sv sepb0_B large_left_vcf_bp) $(sv sepb0_B large_left_vcf_no_sequence) $(cn_exact sepb0)" "0 0 6:0"

sed 's/^P hX#1#chr1 .*/P hX#1#chr1 S+,A+,B+,A+,X+,B+,T+/' "$OUT/modins.spec" > "$OUT/modmid.spec"
spec modmid --cn
bm -i "$OUT/modmid/g.gfa" -c "$OUT/modmid/g.bubbles.csv" -r ref --calls "$CASES/modmid_missed_in_stretch" \
  -o "$OUT/modmid_missed_B" -q >/dev/null 2>&1
check "a miss inside a copy stretch is counted once (large_left_call_bp:large_left_vcf_only_bp)" \
      "$(sv modmid_missed_B large_left_call_bp):$(sv modmid_missed_B large_left_vcf_only_bp)" "80:0"
# Two decoy deletion records nobody carries, 40 and 80 bp before H and 20 and 40 bp shorter: each explains
# hS's short copy almost as well, and the first and the last lie 80 bp apart, so the VCF does not say which
# stretch is missing and none is applied, whatever the records' order (large_left_vcf_only_bp: H's 200 bp).
awk -F'\t' -v OFS='\t' '!/^#/ && $3 ~ /_DEL_H$/ {
    n=split($0, f, "\t")
    for (k = 10; k <= n; k++) g = g "\t0"
    print $1, $2 - 80, "decoy1", $4, "<DEL>", ".", ".", "END=" $2 + 80 ";SVTYPE=DEL;SVLEN=-160", "GT" g
    print $1, $2 - 40, "decoy2", $4, "<DEL>", ".", ".", "END=" $2 + 140 ";SVTYPE=DEL;SVLEN=-180", "GT" g
  } {print}' "$OUT/shortcopy_c.region.vcf" | noforms > "$OUT/decoys_c.region.vcf"
cp "$OUT/shortcopy_c.carrier_events.tsv" "$OUT/decoys_c.carrier_events.tsv"
bm -i "$OUT/shortcopy/g.gfa" -c "$OUT/shortcopy/g.bubbles.csv" -r "ref#0#chr1:0-5600" --calls "$OUT/decoys_c" \
  -o "$OUT/decoys_B" -q >/dev/null 2>&1
check "deletions that would explain a short copy, a chain of small shifts reaching min_sv apart: none applied (large_left_vcf_only_bp)" \
      "$(sv decoys_B large_left_vcf_only_bp)" "200"
# An inversion is mismatches only: one call misses is a large event left, not 200 bp of SNVs.
cat > "$OUT/inv.spec" <<'EOF'
N S 300 1
N A 200 2
N C 100 3
N T 300 4
P ref#0#chr1:0-900 S+,A+,C+,T+
P h1#1#chr1 S+,A+,C+,T+
P hI#1#chr1 S+,A-,C+,T+
P hK#1#chr1 S+,C+,T+
B S T A;C
EOF
spec inv
check "an inversion counts as a large event present, and call's INV leaves nothing" \
      "$(sv inv_B large_present):$(sv inv_B large_left_call)" "2:0"
grep -v 'INV_' "$OUT/inv_c.region.vcf" > "$OUT/invm_c.region.vcf"
awk -F'\t' '$1 !~ /INV_/' "$OUT/inv_c.carrier_events.tsv" > "$OUT/invm_c.carrier_events.tsv"
bm -i "$OUT/inv/g.gfa" -c "$OUT/inv/g.bubbles.csv" -r 'ref#0#chr1:0-900' --calls "$OUT/invm_c" -o "$OUT/invm_B" -q >/dev/null 2>&1
check "an inversion call missed is a large event left" "$(sv invm_B large_left_call)" "1"
bm -i "$OUT/inv/g.gfa" -c "$OUT/inv/g.bubbles.csv" -r 'ref#0#chr1:0-900' --calls "$OUT/inv_c" \
       --vcf "$OUT/inv_c.region.vcf" -o "$OUT/x" -q >/dev/null 2>&1
check "--vcf is not an option: both stages come from --calls" "$?" "2"
# A 9 kb inversion between collapsed inverted repeats, longer than its two flanks together. A haplotype's
# direction comes from its ends: read the way most of its anchors run, the inverted ones were turned
# around and scored over the inverted middle alone, and the inversion was neither present nor left.
# No snarl holds it; bubble makes it an inversion site, which call types.
cat > "$OUT/binv.spec" <<'EOF'
N S 2000 1
N I 300 2
N M1 3000 3
N M2 3000 4
N M3 3000 5
N T 2000 6
P ref#0#chr1:0-13600 S+,I+,M1+,M2+,M3+,I-,T+
P h1#1#chr1 S+,I+,M1+,M2+,M3+,I-,T+
P hI#1#chr1 S+,I+,M3-,M2-,M1-,I-,T+
EOF
BUBBLE=1 spec binv
# (an inversion aligns as dense mismatches and comes out as several large events, so counts are > 0, not 1)
check "an inversion longer than its flanks is scored: present, and left by nothing once call types its site" \
      "$([ "$(sv binv_B large_present)" -gt 0 ] && echo present):$(sv binv_B large_left_call) $(awk -F'\t' '$1 ~ /^hI/{print $3"-"$4}' "$OUT/binv_B.haplotypes.tsv")" "present:0 1-13600"
grep -v 'INV' "$OUT/binv_c.region.vcf" > "$OUT/binvm_c.region.vcf"
awk -F'\t' 'NR == 1 || $3 != "inv"' "$OUT/binv_c.carrier_events.tsv" > "$OUT/binvm_c.carrier_events.tsv"
bm -i "$OUT/binv/g.gfa" -c "$OUT/binv/g.bubbles.csv" -r 'ref#0#chr1:0-13600' --calls "$OUT/binvm_c" -o "$OUT/binvm_B" -q >/dev/null 2>&1
check "and left over when no record has it" "$([ "$(sv binvm_B large_left_call)" -gt 0 ] && echo left)" "left"
# The same inversion with a 200 bp D in its middle, which hID deletes and inverts (and hD only deletes):
# the deletion is applied inside the inversion, not left out as a conflict with it.
cat > "$OUT/binvd.spec" <<'EOF'
N S 2000 1
N I 300 2
N M1 3000 3
N M2 1400 4
N D 200 7
N M4 1400 8
N M3 3000 5
N T 2000 6
P ref#0#chr1:0-13600 S+,I+,M1+,M2+,D+,M4+,M3+,I-,T+
P h1#1#chr1 S+,I+,M1+,M2+,D+,M4+,M3+,I-,T+
P hI#1#chr1 S+,I+,M3-,M4-,D-,M2-,M1-,I-,T+
P hID#1#chr1 S+,I+,M3-,M4-,M2-,M1-,I-,T+
P hD#1#chr1 S+,I+,M1+,M2+,M4+,M3+,I-,T+
EOF
BUBBLE=1 spec binvd
check "an event inside an inversion is applied inside it, in both stages (left:lost:conflicts:clean)" \
      "$(sv binvd_B large_left_call):$(sv binvd_B large_left_vcf_only):$(lg binvd_B "$CONFLICTS_CALL"):$(sv binvd_B haplotypes_clean_vcf)" "0:0:0:4"

# ---- sizes, alignment and attribution -------------------------------------------------------------------
# An inversion is sized by its length: seven of 60 to 300 bp, one per haplotype. The bases at its ends that
# happen to match the reference are not part of the event, so the total is a little under 1,010.
cat > "$OUT/invs.spec" <<'EOF'
N S 400 1
N A60 60 60
N F60 400 67
N A80 80 80
N F80 400 87
N A100 100 100
N F100 400 107
N A120 120 120
N F120 400 127
N A150 150 150
N F150 400 157
N A200 200 200
N F200 400 207
N A300 300 300
N F300 400 307
N T 400 2
P ref#0#chr1 S+,A60+,F60+,A80+,F80+,A100+,F100+,A120+,F120+,A150+,F150+,A200+,F200+,A300+,F300+,T+
P h60#1#chr1 S+,A60-,F60+,A80+,F80+,A100+,F100+,A120+,F120+,A150+,F150+,A200+,F200+,A300+,F300+,T+
P h80#1#chr1 S+,A60+,F60+,A80-,F80+,A100+,F100+,A120+,F120+,A150+,F150+,A200+,F200+,A300+,F300+,T+
P h100#1#chr1 S+,A60+,F60+,A80+,F80+,A100-,F100+,A120+,F120+,A150+,F150+,A200+,F200+,A300+,F300+,T+
P h120#1#chr1 S+,A60+,F60+,A80+,F80+,A100+,F100+,A120-,F120+,A150+,F150+,A200+,F200+,A300+,F300+,T+
P h150#1#chr1 S+,A60+,F60+,A80+,F80+,A100+,F100+,A120+,F120+,A150-,F150+,A200+,F200+,A300+,F300+,T+
P h200#1#chr1 S+,A60+,F60+,A80+,F80+,A100+,F100+,A120+,F120+,A150+,F150+,A200-,F200+,A300+,F300+,T+
P h300#1#chr1 S+,A60+,F60+,A80+,F80+,A100+,F100+,A120+,F120+,A150+,F150+,A200+,F200+,A300-,F300+,T+
B S F60 A60
B F60 F80 A80
B F80 F100 A100
B F100 F120 A120
B F120 F150 A150
B F150 F200 A200
B F200 F300 A300
EOF
spec invs
check "seven inversions are seven large events, present, and none left once called" \
      "$(sv invs_B large_present):$([ "$(sv invs_B large_present_bp)" -ge 900 ] && echo sized):$(sv invs_B large_left_call)" "7:sized:0"
# A 60 bp replacement by 45 other bases is a large event, sized 60, in each of ten haplotypes.
cat > "$OUT/sub.spec" <<'EOF'
N S 300 1
N A 60 2
N B1 45 101
N B2 45 102
N B3 45 103
N B4 45 104
N B5 45 105
N B6 45 106
N B7 45 107
N B8 45 108
N B9 45 109
N B10 45 110
N T 300 3
P ref#0#chr1:0-660 S+,A+,T+
P h1#1#chr1 S+,B1+,T+
P h2#1#chr1 S+,B2+,T+
P h3#1#chr1 S+,B3+,T+
P h4#1#chr1 S+,B4+,T+
P h5#1#chr1 S+,B5+,T+
P h6#1#chr1 S+,B6+,T+
P h7#1#chr1 S+,B7+,T+
P h8#1#chr1 S+,B8+,T+
P h9#1#chr1 S+,B9+,T+
P h10#1#chr1 S+,B10+,T+
B S T A;B1;B2;B3;B4;B5;B6;B7;B8;B9;B10;
EOF
spec sub
check "ten 60-for-45 bp replacements are large" "$(sv sub_B large_present):$([ "$(sv sub_B large_present_bp)" -ge 580 ] && echo sized)" "10:sized"
# A replacement whose new sequence the VCF does not hold cannot be rebuilt from it.
awk -F'\t' -v OFS='\t' '$3 == "bubble1_SUB_A" { $4 = substr($4, 1, 1); $5 = "<DEL>"; $8 = $8 ";SUBSTITUTION" } 1' \
    "$OUT/sub_c.region.vcf" > "$OUT/subn_c.region.vcf"
cp "$OUT/sub_c.carrier_events.tsv" "$OUT/subn_c.carrier_events.tsv"
bm -i "$OUT/sub/g.gfa" -c "$OUT/sub/g.bubbles.csv" -r 'ref#0#chr1:0-660' --calls "$OUT/subn_c" -o "$OUT/subn_B" -q >/dev/null 2>&1
check "a replacement with no INSSEQ is not rebuilt from the VCF, and what it leaves is lost (no_sequence)" \
      "$([ "$(lg subn_B "$UNSPELLED")" -ge 1 ] && echo counted):$(lost_rows subn | cut -d' ' -f4 | sort -u | paste -sd, -)" "counted:no_sequence"
check "the warning names the cause: a sequence call left out above --max-seq-bp (lines)" \
      "$(grep -c 'call left out above its --max-seq-bp' "$OUT/subn_B.log"):$(grep -c 'no reference copy' "$OUT/subn_B.log")" "1:0"
check "every haplotype is clean once its replacement is applied" "$(sv sub_B haplotypes_clean_call):$(sv sub_B haplotypes_clean_vcf):$(sv sub_B haplotypes_scored)" "10:10:10"
# Two similar segments in the other order: 31 differences each, not a 600 bp deletion and insertion.
cat > "$OUT/swap.spec" <<'EOF'
N S 500 1
N A 600 2
N B @A:31:7
N X @B:0:1
N Y @A:0:1
N T 500 3
P ref#0#chr1:0-2200 S+,A+,B+,T+
P href#1#chr1 S+,A+,B+,T+
P hswap#1#chr1 S+,X+,Y+,T+
B S T A;B;X;Y
EOF
spec swap
check "swapped similar segments are small differences, not a large deletion and insertion" \
      "$(sv swap_B large_present):$([ "$(sv swap_B snv_present)" -ge 50 ] && echo snvs)" "0:snvs"
# A 1 kb segment moved past a 5 kb one: the exact blocks are out of order, the cut alignment would
# make it a 5 kb deletion and insertion, and the exact one is 1 kb deleted and 1 kb inserted.
cat > "$OUT/tp2.spec" <<'EOF'
N S 2000 1
N A 2000 2
N X 1000 32
N B 5000 12
N C 2000 5
N T 2000 6
P ref#0#chr1 S+,A+,X+,B+,C+,T+
P h1#1#chr1 S+,A+,X+,B+,C+,T+
P hm#1#chr1 S+,A+,B+,X+,C+,T+
EOF
spec tp2
check "a moved segment is its own size deleted and inserted, not its neighbour's" "$(sv tp2_B large_present_bp)" "2000"
# Lower-case sequence: call writes REF upper-case; --vcf compares it ignoring case.
cat > "$OUT/lc.spec" <<'EOF'
N S =GGATCACAGTCTACACTGCTCACTCCAACCCCGGCCCCTGAGTCCGAGGAGAGGGTGCTTCAGAGTATGTATACCACTGGGTAGGATACGGCGGAGGGCACGTCAATACGGTTCAATGCCCTACTGCATGCTCTTGTGGTTCATCTGCATGGAGAGGGTGGGCATGGGTGGGGGTGCTGGCCCGTGATCTGGACCTCCCATCCACAGCTCATTGTACCGAGTGTAGAGAGGGGCTTGTCCTTCCAGATAGCGTTTCTGTTTCGGTGTAGGTGCTAATCGACTATGCTACTGCGGTTAACG
N A =gggatggcaagtacattttttcgtagatgtgccttgctaacgaaagtattaaacacgtcc
N B =ctcacaatagaatcatagttggacgcgcgacggccgttccagaaaatctttgaat
N T =ACTCAATCCTGCGGGTTCGGTGACCTAAAACCCATTGATTGTGTTACCCAGTTCGAGCGCATAGGGAATTCAGGTCCACACATGGCTGGATCCCCATGATATTCAAGAACTATACATTAAGTTGAACCTCCAGAACACATGTTTCAGTCACGTAGTGCCATCATCGATCACGGAATGTAGCATCAATGATCGAGCCGTGGAAAAAACGTGACTCGCGGACCAGCCTTTAGGTCTTCTACTTAACTACAACTGTTCCGCGGCGGCATTGCCCTTAACTAGCGTTACTAACTAGAGTTTTAC
P ref#0#chr1 S+,A+,T+
P h1#1#chr1 S+,A+,T+
P h2#1#chr1 S+,B+,T+
B S T A;B
EOF
spec lc
check "the VCF stage runs on a soft-masked graph" "$(sv lc_B haplotypes_scored):$(sv lc_B large_left_vcf_only)" "2:0"
# An insertion just past a haplotype's shared stretch is outside it, not applied.
cat > "$OUT/edge.spec" <<'EOF'
N S 300 1
N A 300 2
N X 100 3
N T 300 4
N U 300 5
P ref#0#chr1:0-1200 S+,A+,T+,U+
P h1#1#chr1 S+,A+,T+,U+
P hX#1#chr1 S+,A+,X+,T+,U+
P hT#1#chr1 S+,A+,X+,T+,T+
P hP#1#chr1 S+,A+,X+
P hQ#1#chr1 X+,T+,U+
B A T X
EOF
spec edge
check "and the VCF stage leaves it out too: nothing lost" "$(sv edge_B large_left_vcf_only)" "0"
check "an insertion past the stretch a haplotype shares is outside it, and leaves nothing" \
      "$([ "$(lg edge_B "$OUTSIDE")" -ge 1 ] && echo outside):$(awk -F'\t' '$1 ~ /^hT/{print $7}' "$OUT/edge_B.haplotypes.tsv")" "outside:0"
# A haplotype that shares no node with the reference is not scored, and the log's count says so.
cat > "$OUT/noanc.spec" <<'EOF'
N S 300 1
N D 60 2
N F 300 3
N X 400 4
P ref#0#chr1:0-660 S+,D+,F+
P h1#1#chr1 S+,D+,F+
P hD#1#chr1 S+,F+
P hX#1#chr1 X+
B S F D
EOF
spec noanc
"$BIN" benchmark -i "$OUT/noanc/g.gfa" -c "$OUT/noanc/g.bubbles.csv" -r 'ref#0#chr1:0-660' --calls "$OUT/noanc_c" \
       -o "$OUT/noanc_L" 2> "$OUT/noanc_L.log" >/dev/null
check "the log counts only the haplotypes scored, and warns of the one skipped (log:summary:skipped)" \
      "$(sed -n 's/.*haplotypes scored: \([0-9]*\);.*/\1/p' "$OUT/noanc_L.log"):$(sv noanc_L haplotypes_scored):$(grep -c 'warning: 1 haplotype(s) share no anchor' "$OUT/noanc_L.log")" \
      "2:2:1"
# A miss inside an inversion site is in that site, not outside every site.
check "a miss inside an inversion site is not_reported there" \
      "$(tail -n +2 "$OUT/binvm_B.events.tsv" | awk -F'\t' '$2=="call"{print $9, ($10 == "." ? "nowhere" : "in a site")}' | sort -u | paste -sd, -)" "not_reported in a site"
# Three overlapping events of one haplotype, under its record: A over P..R (its own sequence there), C the
# deletion of Q inside A, B over R..W (no change) and longer than A. B replaces A; C does not overlap B, so
# it is applied, and nothing is left.
cat > "$OUT/three.spec" <<'EOF'
N S 300 1
N P 100 2
N Q 60 3
N R 100 4
N W 200 5
N T 300 6
P ref#0#chr1:0-1060 S+,P+,Q+,R+,W+,T+
P h1#1#chr1 S+,P+,Q+,R+,W+,T+
P hD#1#chr1 S+,P+,R+,W+,T+
B P R Q
EOF
spec three
nseq() { awk -v n="$2" '$1=="S" && $2==n{print $3}' "$OUT/$1/g.gfa"; }
{ printf 'bubble1_DEL_Q\t1\tseq\t301\t560\t.\t.\t-60\t200\t%s\thD#1#chr1\n' "$(nseq three P)$(nseq three R)"
  printf 'bubble1_DEL_Q\t1\tseq\t461\t760\t.\t.\t0\t300\t%s\thD#1#chr1\n' "$(nseq three R)$(nseq three W)"
} >> "$OUT/three_c.carrier_events.tsv"
bm -i "$OUT/three/g.gfa" -c "$OUT/three/g.bubbles.csv" -r 'ref#0#chr1:0-1060' --calls "$OUT/three_c" -o "$OUT/three_B" -q >/dev/null 2>&1
check "an event left out for one later replaced is applied when it overlaps no kept one (left:applied:conflicts)" \
      "$(sv three_B large_left_call):$(lg three_B "$APPLIED_CALL"):$(lg three_B "$CONFLICTS_CALL")" "0:2:1"
# Longest first, whatever lies around: with D over W..T (no change, longer than B) and an empty insertion
# inside A too, D and A are kept and C, the insertion and B are left out, in either order of the rows.
rseq() { awk '$1=="S"{printf "%s", $3}' "$OUT/three/g.gfa" | cut -c"$1"-"$2"; }   # S P Q R W T in order
qrow() { awk -F'\t' 'NR==1 || ($4==401 && $5==460)' "$OUT/three_c.carrier_events.tsv"; }   # call's own: Q deleted
row() { printf 'bubble1_DEL_Q\t1\t%s\t%s\t%s\t.\t.\t%s\t%s\t%s\thD#1#chr1\n' "$@"; }
bench3() { cp "$OUT/three_c.region.vcf" "$OUT/$1_c.region.vcf"
           bm -i "$OUT/three/g.gfa" -c "$OUT/three/g.bubbles.csv" -r 'ref#0#chr1:0-1060' --calls "$OUT/$1_c" \
                  -o "$OUT/$1_B" -q >/dev/null 2>&1; }
lac() { echo "$(sv "$1_B" large_left_call):$(lg "$1_B" "$APPLIED_CALL"):$(lg "$1_B" "$CONFLICTS_CALL")"; }
{ cat "$OUT/three_c.carrier_events.tsv"; row seq 450 449 0 0 .; row seq 561 1000 0 440 "$(rseq 561 1000)"; } \
  > "$OUT/four_c.carrier_events.tsv"
awk 'NR==1{print; next} {a[NR]=$0} END{for (i=NR; i>1; i--) print a[i]}' "$OUT/four_c.carrier_events.tsv" \
  > "$OUT/fourrev_c.carrier_events.tsv"
# The same inside W with a 2 bp inversion that turns the reference: the 15 bp event over it and the 21 bp
# one past it are kept; the inversion and the events inside the 15 bp one are left out. All but the
# inversion spell the reference, so applying the inversion would leave its 2 bp as mismatches.
{ qrow
  for iv in 591:605 598:605 603:604 605:622 606:626; do
    row seq "${iv%:*}" "${iv#*:}" 0 "$(( ${iv#*:} - ${iv%:*} + 1 ))" "$(rseq "${iv%:*}" "${iv#*:}")"
  done
  row inv 604 605 2 2 "$(rseq 604 605 | rev | tr ACGT TGCA)"
} > "$OUT/inv6_c.carrier_events.tsv"
# Ties: of two crossing intervals of one length the first is kept (E deletes Q, F after it does not; F's
# row comes first), and of two over one interval the longer replacement (G spells Q, so Q's deletion is left).
{ qrow | head -1; row seq 401 560 0 160 "$(rseq 401 560)"; row seq 301 460 -60 100 "$(rseq 301 400)"; qrow | tail -n +2; } \
  > "$OUT/tie_c.carrier_events.tsv"
{ qrow; row seq 401 460 0 60 "$(rseq 401 460)"; } > "$OUT/tieg_c.carrier_events.tsv"
for p in four fourrev inv6 tie tieg; do bench3 "$p"; done
check "of overlapping events the longer is kept, in either row order (left:applied:conflicts)" \
      "$(lac four)/$(lac fourrev)" "0:2:3/0:2:3"
check "an event longer than an inversion it holds is kept over it (left:applied:conflicts:snvs left)" \
      "$(lac inv6):$(sv inv6_B snv_left_call)" "0:3:4:0"
check "ties: the first of two crossing intervals, the longer of two replacements (left:applied:conflicts)" \
      "$(lac tie)/$(lac tieg)" "0:1:2/1:1:1"
# Why the VCF loses what call found, from one 200 bp deletion of Q (before a 300 bp tandem repeat R) with
# its record edited: deleting only Q's last 20 bp leaves the record's own allele short (representative),
# the same as IMPRECISE (other), and deleting 60 bp of R too leaves a gap the alignment places at R's
# far end; it can slide back along R to the record, so it is that record's allele (representative).
cat > "$OUT/lostr.spec" <<EOF
N S 300 1
N P 100 2
N Q 200 3
N R =$(printf 'ACGTTGCAAG%.0s' {1..30})
N W 200 5
N T 300 6
P ref#0#chr1:0-1400 S+,P+,Q+,R+,W+,T+
P h1#1#chr1 S+,P+,Q+,R+,W+,T+
P hD#1#chr1 S+,P+,R+,W+,T+
B P R Q
EOF
spec lostr
for e in rep:'$2=580' imp:'$2=580; $8=$8";IMPRECISE"' oth:'sub(/END=[0-9]+/, "END=660", $8)'; do
  awk -F'\t' -v OFS='\t' "!/^#/{${e#*:}} 1" "$OUT/lostr_c.region.vcf" > "$OUT/${e%%:*}_c.region.vcf"
  cp "$OUT/lostr_c.carrier_events.tsv" "$OUT/${e%%:*}_c.carrier_events.tsv"
  bm -i "$OUT/lostr/g.gfa" -c "$OUT/lostr/g.bubbles.csv" -r 'ref#0#chr1:0-1400' --calls "$OUT/${e%%:*}_c" \
         -o "$OUT/${e%%:*}_B" -q >/dev/null 2>&1
done
check "lost in the VCF: a carried record's other allele, an IMPRECISE one, a gap sliding back to the record (reasons)" \
      "$(lost_rows rep | cut -d' ' -f4):$(lost_rows imp | cut -d' ' -f4):$(lost_rows oth | cut -d' ' -f4)" \
      "representative:other:representative"
# Records are named by unique IDs.
awk -F'\t' -v OFS='\t' '!/^#/ { $3 = "." } 1' "$OUT/strand_c.region.vcf" > "$OUT/dot_c.region.vcf"
cp "$OUT/strand_c.carrier_events.tsv" "$OUT/dot_c.carrier_events.tsv"
bm -i "$OUT/strand_b.sorted.gfa" -c "$OUT/strand_b.bubbles.csv" -r "$R" --calls "$OUT/dot_c" -o "$OUT/dot" -q >/dev/null 2>&1
check "a VCF record with no ID is refused" "$?" "2"

# ---- refusals --------------------------------------------------------------------------------------------
# The VCF and the carrier events must describe one run.
pipeline strand stale
awk -F'\t' 'NR==1 || $3!="seq" || $9!=70' "$OUT/stale_c.carrier_events.tsv" > "$OUT/x.tsv" && mv "$OUT/x.tsv" "$OUT/stale_c.carrier_events.tsv"
bench stale
check "carrier events and a VCF from different runs are refused" "$?" "2"
bm -i "$OUT/strand_b.sorted.gfa" -c "$OUT/strand_b.bubbles.csv" -r "$R" --calls "$OUT/strand_c" \
       --variant-nodes "$OUT/strand_c.variant_nodes.tsv" -o "$OUT/old" -q >/dev/null 2>&1
check "the previous benchmark's flags are refused, pointing to --calls" "$?" "2"
# The site table must be the graph's: one naming nodes the graph lacks, or one id twice, is refused.
awk -F, -v OFS=, 'NR > 1 { $2 = "X" $2; $4 = "X" $4 } 1' "$OUT/strand_b.bubbles.csv" > "$OUT/foreign.csv"
{ cat "$OUT/strand_b.bubbles.csv"; sed -n 2p "$OUT/strand_b.bubbles.csv"; } > "$OUT/twice.csv"
for t in foreign twice; do
  bm -i "$OUT/strand_b.sorted.gfa" -c "$OUT/$t.csv" -r "$R" --calls "$OUT/strand_c" -o "$OUT/$t" -q >/dev/null 2>&1
  eval "rc_$t=\$?"
done
check "a site table naming nodes the graph lacks, or an id twice, is refused" "$rc_foreign:$rc_twice" "2:2"
# A --diploid VCF has a column per sample: refused as such, even where each sample has one haplotype.
"$BIN" call -i "$OUT/strand_b.sorted.gfa" -c "$OUT/strand_b.bubbles.csv" -r "$R" -o "$OUT/dip_c" --no-per-bubble-vcf \
       --diploid -q >/dev/null 2>&1
"$BIN" benchmark -i "$OUT/strand_b.sorted.gfa" -c "$OUT/strand_b.bubbles.csv" -r "$R" --calls "$OUT/dip_c" -o "$OUT/dip" -q \
       >/dev/null 2> "$OUT/dip.err"
check "a --diploid VCF is refused as one (rc:says so)" "$?:$(grep -c 'without --diploid' "$OUT/dip.err")" "2:1"

printf "%d assertion(s) failed\n" "$fails"
[ "$fails" -eq 0 ]
