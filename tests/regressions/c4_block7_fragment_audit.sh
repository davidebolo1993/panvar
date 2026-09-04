#!/usr/bin/env bash
# c4_block7_fragment_audit.sh - which fragments actually decide c4/NA19240, under COMPLETE in-band
# search rather than the accelerated placer.
#
#   c4_block7_fragment_audit.sh <out-dir> <repo> <panvar-binary>
#
# THE QUESTION. The accelerated production scorer puts the truth-equivalent pair 499.17 nats behind
# a competitor, and that total is carried by 76 of 23953 fragments: 23 the competitor alone explains
# (-3618.56) against 53 only truth explains (+2384.39). Those labels come from the accelerated
# placer. This asks whether they survive a search that is complete within the divergence band.
#
# WHAT IT CAN AND CANNOT RECONCILE. Searching only the 76 selected fragments reconciles to the
# AUDIT-SUBSET delta (-1234.17), NOT to the whole-dataset -499.17. The other 23877 fragments are not
# assumed unchanged: complete search can find additional or better placements among fragments
# previously called `both` or `neither` too. Reconciling the full -499.17 needs a bounded-only run
# over all 23953 fragments, which is a separate exercise. The column is named audit_subset_delta for
# that reason.
#
# PROVENANCE IS PINNED BEFORE ANY MEASUREMENT, because the inputs are regenerated rather than stored:
# read files, the production dump, candidate names with their authoritative walk md5s, the fragment
# name list, and the parameters. The scratchpad is ephemeral -- an earlier run of this audit lost
# every input to a cleanup, which is why the script is the artifact and the data is not.
set -uo pipefail
OUT="${1:?out}"; REPO="${2:?repo}"; BIN="${3:?bin}"
mkdir -p "$OUT"
PY="${PYTHON:-python3}"
md5of() { if command -v md5 >/dev/null 2>&1; then md5 -q "$1"; else md5sum "$1" | cut -d' ' -f1; fi; }
fails=0
ok()  { printf '  ok   %s\n' "$*"; }
bad() { printf '  FAIL %s\n' "$*"; fails=$((fails+1)); }

G="$REPO/results/real_data/c4/bubble/bubble.sorted.gfa"
P="$REPO/results/real_data/c4/bubble/bubble"
T1="NA19240#1#haplotype1-0000009:31981961-32208331"
TE="NA19238#1#haplotype1-0000026:31825469-32078164"
WC="HG00171#1#haplotype1-0000033:140888062-141114404"
# C4 PARAMETERS, frozen. 150 bp mates at 5% give mate_band_edits = floor(0.05*150)+1 = 8, and the
# shared insert support is [max(r1+r2), mean+4sd] = [300,550].
SEED=7; COV=30; DIV=0.05; FLEN=350; FSD=50; ERR=0.001
EXP_D=8; EXP_ILO=300; EXP_IHI=550
EXP_COMP=23; EXP_TRUTH=53; EXP_SUBSET=-1234.17

# ---- truth haplotypes, spelled from the graph walk ---------------------------------------------
"$BIN" genotype-frag -i "$G" -b "$P" -o "$OUT/dump" --dump-scored-sequences "$OUT/dump" -q \
  >/dev/null 2>&1
if [ ! -s "$OUT/dump.scored_sequences.fa" ]; then
  echo "FATAL: could not spell the c4 panel"; exit 1
fi
"$PY" - "$OUT/dump.scored_sequences.fa" "$OUT/dump.scored_sequences.tsv" "$OUT" \
       "$T1" "$TE" "$WC" > "$OUT/candidates.tsv" <<'PYEOF'
import sys, hashlib
fa,tsv,out,t1,te,wc=sys.argv[1:7]
frame={};h=None
for l in open(tsv):
    f=l.rstrip('\n').split('\t')
    if h is None: h=f; continue
    d=dict(zip(h,f)); frame[d['name']]=d.get('frame','fwd')
seqs={};n=None;b=[]
def flush():
    if n is not None: seqs[n]=''.join(b)
for l in open(fa):
    if l[0]=='>': flush(); n=l[1:].split()[0]; b=[]
    else: b.append(l.strip())
flush()
def rc(s): return s[::-1].translate(str.maketrans('ACGTN','TGCAN'))
print("role\tname\twalk_md5\tbp\tframe")
for role,nm in (("truth_h1",t1),("truth_equiv_h2",te),("competitor",wc)):
    if nm not in seqs:
        sys.stderr.write("missing candidate: %s\n"%nm); sys.exit(3)
    s=seqs[nm].upper()
    o=rc(s) if frame.get(nm)=='rc' else s
    open(f"{out}/{role}.fa","w").write(">%s\n%s\n"%(role,o))
    print("%s\t%s\t%s\t%d\t%s"%(role,nm,hashlib.md5(o.encode()).hexdigest(),len(o),frame.get(nm,'fwd')))
PYEOF
[ -s "$OUT/candidates.tsv" ] || { echo "FATAL: candidate spelling failed"; exit 1; }
cat "$OUT/truth_h1.fa" "$OUT/truth_equiv_h2.fa" > "$OUT/truth.fa"
L=$(awk '/^>/{next}{n+=length($0)}END{print n}' "$OUT/truth.fa")
N=$(( COV * L / 600 ))
wgsim -N "$N" -1 150 -2 150 -d "$FLEN" -s "$FSD" -e "$ERR" -r 0 -R 0 -X 0 -S "$SEED" \
  "$OUT/truth.fa" "$OUT/r1.fq" "$OUT/r2.fq" >/dev/null 2>&1

{ echo "commit            $(cd "$REPO" && git rev-parse --short HEAD 2>/dev/null || echo unknown)"
  echo "binary_md5        $(md5of "$BIN")"
  echo "graph_md5         $(md5of "$G")"
  echo "bubbles_md5       $(md5of "$P.bubbles.csv")"
  echo "truth_fa_md5      $(md5of "$OUT/truth.fa")"
  echo "r1_md5            $(md5of "$OUT/r1.fq")"
  echo "r2_md5            $(md5of "$OUT/r2.fq")"
  echo "seed              $SEED"
  echo "coverage          ${COV}x"
  echo "read_pairs        $N"
  echo "params            div=$DIV flen=$FLEN fsd=$FSD err=$ERR"
  echo "expected_band     d=$EXP_D per 150 bp mate"
  echo "expected_insert   [$EXP_ILO,$EXP_IHI]"
  echo "date              $(date -u +%Y-%m-%dT%H:%M:%SZ)"; } > "$OUT/PROVENANCE.txt"
cat "$OUT/PROVENANCE.txt"; echo; cat "$OUT/candidates.tsv"; echo

# ---- REPRODUCE THE ORIGINAL PARTITION, on the current default-off binary ------------------------
# The audit is only meaningful if the thing being explained still happens. Production scoring has
# been touched since the original measurement (the FR downstream check was factored into a shared
# helper), so this is a real check, not a formality.
for TAG in TE WC; do
  case $TAG in TE) P2="$TE";; WC) P2="$WC";; esac
  "$BIN" genotype-frag -i "$G" -b "$P" -o "$OUT/$TAG" -R "$OUT/r1.fq" -R "$OUT/r2.fq" \
    --haplotype-mode --hamming-emission --max-divergence "$DIV" --fragment-len "$FLEN" \
    --fragment-sd "$FSD" --error-rate "$ERR" --max-anchor-occ 32 --placement-topk 8 \
    --top-pairs 4000 --force-haplotypes "$T1,$TE,$WC" \
    --dump-fragment-mass "$OUT/$TAG.mass.tsv" --dump-mass-pair "$T1,$P2" -q >/dev/null 2>&1
done
if [ ! -s "$OUT/TE.mass.tsv" ] || [ ! -s "$OUT/WC.mass.tsv" ]; then
  bad "the production mass dumps were not written; the partition cannot be reproduced"
else
  "$PY" - "$OUT/TE.mass.tsv" "$OUT/WC.mass.tsv" "$OUT/partition.tsv" "$OUT/deciders.txt" \
    > "$OUT/partition_summary.txt" <<'PYEOF'
import sys
te,wc,outp,outd=sys.argv[1:5]
def load(p):
    rows={};hdr=None
    for l in open(p):
        if l.startswith('#'): continue
        f=l.rstrip('\n').split('\t')
        if hdr is None: hdr=f; continue
        d=dict(zip(hdr,f)); rows[d['fragment']]=d
    return rows
A=load(te); B=load(wc)
common=sorted(set(A)&set(B))
cls={}; dec=[]
for k in common:
    a,b=A[k],B[k]
    ea=float(a['ll_b'])>float(a['floor'])+1e-9
    eb=float(b['ll_b'])>float(b['floor'])+1e-9
    c=('both' if (ea and eb) else 'truth_only' if ea else 'comp_only' if eb else 'neither')
    d=float(a['contrib'])-float(b['contrib'])
    e=cls.setdefault(c,[0,0.0]); e[0]+=1; e[1]+=d
    if c in ('truth_only','comp_only'): dec.append((k,c,d))
with open(outp,'w') as f:
    f.write("class\tn\tsum_delta\n")
    for c in ('comp_only','truth_only','both','neither'):
        n,s=cls.get(c,[0,0.0]); f.write("%s\t%d\t%.4f\n"%(c,n,s))
with open(outd,'w') as f:
    for k,c,d in dec: f.write("%s\t%s\t%.4f\n"%(k,c,d))
sub=cls.get('comp_only',[0,0])[1]+cls.get('truth_only',[0,0])[1]
tot=sum(v[1] for v in cls.values())
print("fragments %d" % len(common))
print("comp_only %d %.2f" % tuple(cls.get('comp_only',[0,0.0])))
print("truth_only %d %.2f" % tuple(cls.get('truth_only',[0,0.0])))
print("audit_subset_delta %.2f" % sub)
print("whole_dataset_delta %.2f" % tot)
PYEOF
  NC=$(awk '/^comp_only/{print $2}' "$OUT/partition_summary.txt")
  NT=$(awk '/^truth_only/{print $2}' "$OUT/partition_summary.txt")
  SUB=$(awk '/^audit_subset_delta/{print $2}' "$OUT/partition_summary.txt")
  TOT=$(awk '/^whole_dataset_delta/{print $2}' "$OUT/partition_summary.txt")
  NF2=$(awk '/^fragments/{print $2}' "$OUT/partition_summary.txt")
  echo "  reproduced: fragments=$NF2 comp_only=$NC truth_only=$NT subset=$SUB whole=$TOT"
  [ "${NC:-0}" = "$EXP_COMP" ] && ok "competitor-only count reproduces ($NC)" \
                               || bad "competitor-only is $NC, original was $EXP_COMP"
  [ "${NT:-0}" = "$EXP_TRUTH" ] && ok "truth-only count reproduces ($NT)" \
                                || bad "truth-only is $NT, original was $EXP_TRUTH"
  "$PY" -c "import sys; sys.exit(0 if abs(float('$SUB')-($EXP_SUBSET))<0.5 else 1)" 2>/dev/null \
    && ok "audit_subset_delta reproduces ($SUB vs $EXP_SUBSET)" \
    || bad "audit_subset_delta is $SUB, original was $EXP_SUBSET"
  NDEC=$(wc -l < "$OUT/deciders.txt" | tr -d ' ')
  [ "${NDEC:-0}" = 76 ] && ok "exactly 76 deciding fragments selected by name" \
                        || bad "$NDEC deciding fragments, expected 76"
fi

if [ "$fails" -ne 0 ]; then
  echo; echo "c4 audit preconditions: $fails FAILED -- the audit must not run"; exit "$fails"
fi
echo "  preconditions pinned and reproduced"; echo

# ---- THE AUDIT: reclassify the 76 under complete in-band search ---------------------------------
# Restricted to the three named haplotypes. Searching all 131 panel paths is unnecessary to classify
# two named pairs, and the exhaustive half of the A/B scales with the panel.
{ printf 'H\tVN:Z:1.0\n'
  i=0
  for R in truth_h1 truth_equiv_h2 competitor; do
    i=$((i+1))
    printf 'S\t%d\t%s\n' "$i" "$(awk 'NR>1' "$OUT/$R.fa" | tr -d '\n')"
  done
  printf 'L\t1\t+\t2\t+\t0M\nL\t2\t+\t3\t+\t0M\n'
  printf 'P\ttruth_h1\t1+\t*\n'; printf 'P\ttruth_equiv_h2\t2+\t*\n'
  printf 'P\tcompetitor\t3+\t*\n'; } > "$OUT/three.gfa"
"$BIN" bubble -i "$OUT/three.gfa" -r truth_h1 -o "$OUT/three" --min-variant-bp 0 -q >/dev/null 2>&1

# The 76 fragments BY NAME from the production dump, never rediscovered under the new search.
cut -f1 "$OUT/deciders.txt" | sort > "$OUT/want76.txt"
"$PY" - "$OUT/r1.fq" "$OUT/r2.fq" "$OUT/want76.txt" "$OUT/d1.fq" "$OUT/d2.fq" > "$OUT/recover.txt" <<'PYEOF'
import sys
r1,r2,want,o1,o2=sys.argv[1:6]
w=set(l.strip() for l in open(want))
def filt(src,dst):
    got=set()
    with open(src) as f, open(dst,'w') as o:
        while True:
            h=f.readline()
            if not h: break
            s=f.readline(); p=f.readline(); q=f.readline()
            nm=h[1:].strip().split('/')[0]
            if nm in w: o.write(h+s+p+q); got.add(nm)
    return got
g1=filt(r1,o1); g2=filt(r2,o2)
print("wanted %d mate1 %d mate2 %d both %d" % (len(w),len(g1),len(g2),len(g1&g2)))
PYEOF
read -r _ NW _ NM1 _ NM2 _ NB <<<"$(cat "$OUT/recover.txt")"
{ [ "$NW" = 76 ] && [ "$NB" = 76 ]; }   && ok "all 76 selected fragments recovered with BOTH mates from the pinned reads"   || bad "recovered $NB of $NW fragments with both mates (m1=$NM1 m2=$NM2)"

"$BIN" genotype-frag -i "$OUT/three.gfa" -b "$OUT/three" -o "$OUT/bs"   -R "$OUT/d1.fq" -R "$OUT/d2.fq" --max-divergence "$DIV" --fragment-len "$FLEN"   --fragment-sd "$FSD" --bounded-search "$OUT/bs.tsv" -q >/dev/null 2>&1
if [ ! -s "$OUT/bs.tsv.states.tsv" ]; then
  bad "the bounded search produced no fragment states for the 76"
else
  MD=$(awk -F'\t' 'NR>1{print $9; exit}' "$OUT/bs.tsv" 2>/dev/null)
  ILO=$(awk -F'\t' 'NR==2{print $8}' "$OUT/bs.tsv.states.tsv")
  IHI=$(awk -F'\t' 'NR==2{print $9}' "$OUT/bs.tsv.states.tsv")
  [ "${MD:-0}" = "$EXP_D" ] && ok "the audit ran at the production band (d=$MD per 150 bp mate)"                             || bad "band is ${MD:-?}, expected $EXP_D"
  { [ "${ILO:-0}" = "$EXP_ILO" ] && [ "${IHI:-0}" = "$EXP_IHI" ]; }     && ok "and the production insert support ([$ILO,$IHI])"     || bad "insert support [${ILO:-?},${IHI:-?}], expected [$EXP_ILO,$EXP_IHI]"
  "$PY" - "$OUT/bs.tsv.states.tsv" "$OUT/deciders.txt" > "$OUT/reclass.tsv" <<'PYEOF'
import sys, collections
st,dec=sys.argv[1:3]
have=collections.defaultdict(set)
h=None
for l in open(st):
    f=l.rstrip(chr(10)).split(chr(9))
    if h is None: h=f; continue
    d=dict(zip(h,f))
    if int(d['bounded_states'])>0: have[d['fragment']].add(d['haplotype'])
old={}
for l in open(dec):
    k,c,v=l.rstrip(chr(10)).split(chr(9)); old[k]=(c,float(v))
tab=collections.Counter()
print(chr(9).join(["fragment","old_class","new_class","old_delta"]))
for k,(c,v) in sorted(old.items()):
    hs=have.get(k,set())
    te='truth_equiv_h2' in hs; wc='competitor' in hs
    new=('both' if (te and wc) else 'truth_only' if te else 'comp_only' if wc else 'neither')
    tab[(c,new)]+=1
    print(chr(9).join([k,c,new,"%.4f"%v]))
PYEOF
  echo; echo "  old class -> bounded-complete class:"
  "$PY" - "$OUT/reclass.tsv" <<'PYEOF'
import sys, collections
t=collections.Counter(); h=None
for l in open(sys.argv[1]):
    f=l.rstrip(chr(10)).split(chr(9))
    if h is None: h=f; continue
    t[(f[1],f[2])]+=1
for (a,b),n in sorted(t.items()):
    print("    %-11s -> %-11s %3d%s" % (a,b,n," (unchanged)" if a==b else ""))
PYEOF
fi

echo
if [ "$fails" -eq 0 ]; then echo "c4 block-7 audit: complete"; else
  echo "c4 block-7 audit: $fails assertion(s) failed"; fi
exit "$fails"
