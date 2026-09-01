#!/usr/bin/env bash
# genotype_path_blocks.sh - the authoritative path->block projection, on REAL antiparallel paths.
#
#   genotype_path_blocks.sh <panvar-binary> <out-dir>
#
# The projection is the coordinate foundation every per-block measurement stands on, so its gates
# run against real loci rather than synthetic fixtures. c4 and cyp2d6 are used because 60 of 131 and
# 59 of 127 of their paths are ANTIPARALLEL to the block frame -- the case a previous Python
# projector got wrong by assuming source->sink order in walk coordinates, returning ABSENT for every
# block of every reverse path. A fixture with one hand-built reverse path would not have caught the
# scale of that, and the donor validated first (HG00096) is forward at both loci, so a donor-level
# check would not have caught it at all.
set -uo pipefail
BIN="${1:?usage: genotype_path_blocks.sh <panvar> <outdir>}"
OUT="${2:?}"; mkdir -p "$OUT"; OUT="$OUT/run.$$"; rm -rf "$OUT"; mkdir -p "$OUT"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/.." && pwd)"
PY="${PYTHON:-python3}"
fails=0
ok()  { printf '  ok   %s\n' "$*"; }
bad() { printf '  FAIL %s\n' "$*"; fails=$((fails+1)); }

for LOCUS in c4 cyp2d6; do
  G="$REPO/results/real_data/$LOCUS/bubble/bubble.sorted.gfa"
  P="$REPO/results/real_data/$LOCUS/bubble/bubble"
  if [ ! -f "$G" ]; then printf '  skip %s (no substrate)\n' "$LOCUS"; continue; fi
  "$BIN" genotype-frag -i "$G" -b "$P" -o "$OUT/$LOCUS" --dump-scored-sequences "$OUT/$LOCUS" -q \
    >/dev/null 2>&1
  if [ ! -s "$OUT/$LOCUS.path_blocks.tsv" ]; then bad "$LOCUS: no path_blocks.tsv"; continue; fi

  # NON-VACUITY: the antiparallel case must actually be present, or every gate below is about
  # forward paths only.
  NRC=$(awk -F'\t' 'NR>1 && $3=="rc"{n++} END{print n+0}' "$OUT/$LOCUS.path_blocks.tsv")
  [ "$NRC" -gt 0 ] && ok "$LOCUS: $NRC antiparallel slices present, so the reverse case is exercised" \
                   || bad "$LOCUS: no antiparallel slices; these gates would only cover forward paths"

  # GATE: concatenating mapped slices -- reverse-complemented back to walk orientation where the
  # strand is '-' -- plus the explicit unmapped intervals reproduces the authoritative walk BYTE FOR
  # BYTE. This is the contract the whole surface rests on.
  R=$("$PY" - "$OUT/$LOCUS" <<'PY'
import sys
S=sys.argv[1]
def rc(s): return s[::-1].translate(str.maketrans('ACGTN','TGCAN'))
def fa(p):
    d={};n=None;b=[]
    for l in open(p):
        if l[0]=='>':
            if n: d[n]=''.join(b)
            n=l[1:].split()[0]; b=[]
        else: b.append(l.strip())
    if n: d[n]=''.join(b)
    return d
walks=fa(f"{S}.scored_sequences.fa")
sl={};hdr=None;b=[]
for l in open(f"{S}.path_blocks.fa"):
    if l[0]=='>':
        if hdr: sl[hdr]=''.join(b)
        p=l[1:].rstrip('\n').split()
        wb=int([x for x in p if x.startswith('walk=')][0][5:].split('-')[0])
        hdr=(p[0],wb); b=[]
    else: b.append(l.strip())
if hdr: sl[hdr]=''.join(b)
rows={};h=None
for l in open(f"{S}.path_blocks.tsv"):
    f=l.rstrip('\n').split('\t')
    if h is None: h=f; continue
    d=dict(zip(h,f)); rows.setdefault(d['path'],[]).append(d)
ok=bad=0
for path,rs in rows.items():
    if rs[0]['projection_status']=='unprojectable': continue
    w=walks.get(path)
    if w is None: continue
    parts=[]
    for d in rs:
        lo=int(d['walk_begin'])
        if d['projection_status']=='unmapped':
            parts.append((lo,w[lo:int(d['walk_end'])])); continue
        s=sl.get((path,lo),'')
        parts.append((lo, rc(s) if d['walk_strand']=='-' else s))
    parts.sort()
    if ''.join(p[1] for p in parts).upper()==w.upper(): ok+=1
    else: bad+=1
print("%d %d" % (ok,bad))
PY
)
  read -r NOK NBAD <<<"$R"
  { [ "${NBAD:-1}" = 0 ] && [ "${NOK:-0}" -gt 0 ]; } \
    && ok "$LOCUS: $NOK/$NOK paths reconstruct the walk byte-for-byte from their block slices" \
    || bad "$LOCUS: $NBAD path(s) do not reconstruct; the coordinate foundation is wrong"

  # GATE: a retained path is representable in its own catalogue, by construction. This is the
  # control -- if it ever fails, the sequence-equality match is broken and every held-out number
  # below is meaningless.
  PN=$(awk -F'\t' 'NR>1 && $1=="panel" && $4=="complete"{n++} END{print n+0}' "$OUT/$LOCUS.path_blocks.tsv")
  PR=$(awk -F'\t' 'NR>1 && $1=="panel" && $4=="complete" && $15==1{n++} END{print n+0}' "$OUT/$LOCUS.path_blocks.tsv")
  { [ "$PN" -gt 0 ] && [ "$PN" = "$PR" ]; } \
    && ok "$LOCUS: every retained slice ($PR/$PN) is representable in the calling catalogue" \
    || bad "$LOCUS: only $PR of $PN retained slices are representable"

  # GATE: a slice carrying a catalogue index must be representable, and one without must not be.
  INC=$(awk -F'\t' 'NR>1 && $4=="complete" && (($14=="NA" && $15==1) || ($14!="NA" && $15==0)){n++} END{print n+0}' "$OUT/$LOCUS.path_blocks.tsv")
  [ "$INC" = 0 ] \
    && ok "$LOCUS: catalogue_allele and catalogue_representable never disagree" \
    || bad "$LOCUS: $INC slice(s) have an index without representability, or the reverse"
done

# GATE: HELD-OUT REPRESENTABILITY IS DECIDED BY THE REDUCED PANEL. Projecting a held-out path with
# its own block set and then looking its allele up there answers a question nobody asked: it is
# representable by construction. Holding out a real donor must produce BOTH outcomes -- alleles a
# retained path still carries, and alleles unique to the held-out pair.
G="$REPO/results/real_data/c4/bubble/bubble.sorted.gfa"
P="$REPO/results/real_data/c4/bubble/bubble"
if [ -f "$G" ]; then
  T1="HG00096#1#haplotype1-0000024:31848049-32074186"
  T2="HG00096#2#haplotype2-0000117:31856177-32049825"
  "$BIN" genotype-frag -i "$G" -b "$P" -o "$OUT/held" --dump-scored-sequences "$OUT/held" \
    --exclude-haplotypes "$T1,$T2" -q >/dev/null 2>&1
  HN=$(awk -F'\t' 'NR>1 && $1=="held_out" && $4=="complete"{n++} END{print n+0}' "$OUT/held.path_blocks.tsv")
  H1=$(awk -F'\t' 'NR>1 && $1=="held_out" && $4=="complete" && $15==1{n++} END{print n+0}' "$OUT/held.path_blocks.tsv")
  H0=$(awk -F'\t' 'NR>1 && $1=="held_out" && $4=="complete" && $15==0{n++} END{print n+0}' "$OUT/held.path_blocks.tsv")
  [ "${HN:-0}" -gt 0 ] \
    && ok "held-out projection produced $HN slices against the reduced panel" \
    || bad "held-out projection produced no slices"
  [ "${H1:-0}" -gt 0 ] \
    && ok "a held-out allele still carried by a retained path is representable ($H1)" \
    || bad "no held-out allele was representable; the reduced catalogue is not being consulted"
  [ "${H0:-0}" -gt 0 ] \
    && ok "a held-out allele carried by nobody else is NOT representable ($H0)" \
    || bad "every held-out allele was representable -- the tautology this gate exists to catch"
  UNIQ=$(awk -F'\t' 'NR>1 && $1=="held_out" && $4=="complete" && $15==0 && $14!="NA"{n++} END{print n+0}' "$OUT/held.path_blocks.tsv")
  [ "$UNIQ" = 0 ] \
    && ok "an unrepresentable held-out allele carries NA, not a borrowed index" \
    || bad "$UNIQ unrepresentable held-out slice(s) carry a catalogue index"
fi

echo
if [ "$fails" -eq 0 ]; then echo "path blocks: all assertions passed"; else
  echo "path blocks: $fails assertion(s) failed"; fi
exit "$fails"
