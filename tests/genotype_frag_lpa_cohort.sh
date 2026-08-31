#!/usr/bin/env bash
# genotype_frag_lpa_cohort.sh - the LPA development cohort, reproducibly.
#
#   genotype_frag_lpa_cohort.sh <out-dir> [arms] [donors]
#
# WHAT THIS IS FOR. Every LPA accuracy number on this branch was at one time wrong for a reason that
# had nothing to do with the genotyper: the panel was spelled from a different pipeline stage than
# the reads, the call table was spelled against a different allele catalogue, the binary changed
# mid-run, the floor was computed against a candidate list that included the held-out truth itself.
# Each produced numbers that were plausible and wrong. So this harness freezes the things that were
# silently varying, and REFUSES rather than reporting when one of them does not hold.
#
# WHAT IT MEASURES, per donor, leave-one-out:
#   * the drift floor       -- truth spelled from THIS graph against the reads' own FASTA; must be 0
#   * the certified floor   -- min edit distance from any remaining panel haplotype, exact, banded
#   * the called distance   -- the caller's pair, spelled BY NAME, against truth
#   * excess                -- called minus floor, the reconstruction headroom
#
# SUBSTRATE. The bubble stage, always. panphorte folds KIV-2 approximately (LPA runs
# --min-similarity 0.95) and that rewrites ~2100 edits per haplotype; the reads and truth come from
# the bubble stage, so scoring against a later stage measures the folding. See
# scripts/spell_paths.py, which warns about exactly this.
#
# DONORS. Two kinds, and they behave the same:
#   * experiments/cosigt_mapping_work/pair8..pair23 -- SYNTHETIC pairs joining haplotypes of
#     unrelated individuals, except pair8 which is HG00232's own diploid;
#   * assembly-matched diploids built here, both homologues from one sample.
# These are ASSEMBLY-MATCHED SIMULATED reads (wgsim from the assemblies' spelled sequence), not real
# reads. No real reads for this locus exist in the repository.
#
# SCOPE. One locus. Nothing here licenses a claim about the caller in general.
set -uo pipefail

OUT="${1:?usage: genotype_frag_lpa_cohort.sh <out-dir> [arms] [donors]}"
ARMS="${2:-base bandfloor}"
DONORS="${3:-8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/.." && pwd)"
BIN="${PANVAR_BIN:-$REPO/build_genotype_review/panvar}"
W="$REPO/experiments/cosigt_mapping_work"
G="$REPO/results/real_data/lpa/bubble/bubble.sorted.gfa"
PFX="$REPO/results/real_data/lpa/bubble/bubble"
BAND="${BAND:-2048}"
mkdir -p "$OUT"

md5of() { if command -v md5 >/dev/null 2>&1; then md5 -q "$1"; else md5sum "$1" | cut -d' ' -f1; fi; }

# ---- provenance, recorded before anything is measured -------------------------------------------
# A run whose binary or substrate is not the one reported is worse than no run: it looks like
# evidence. Measured the hard way -- a cohort and a ctest run were both invalidated by rebuilding
# underneath them.
for f in "$BIN" "$G" "$PFX.bubbles.csv"; do
  [ -f "$f" ] || { echo "FATAL: missing $f" >&2; exit 1; }
done
BIN_MD5="$(md5of "$BIN")"
{ echo "binary          $BIN_MD5"
  echo "graph           $(md5of "$G")"
  echo "bubbles         $(md5of "$PFX.bubbles.csv")"
  echo "band            $BAND"
  echo "arms            $ARMS"
  echo "donors          $DONORS"
  echo "date            $(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$OUT/PROVENANCE.txt"
cat "$OUT/PROVENANCE.txt"
# Copy the binary: a rebuild during the run must not change what is being measured.
cp "$BIN" "$OUT/panvar.frozen"; BIN="$OUT/panvar.frozen"

"$BIN" genotype-frag --help 2>/dev/null | grep -q -- --dump-scored-sequences || {
  echo "FATAL: this binary has no --dump-scored-sequences; it predates the cohort" >&2; exit 1; }

printf 'donor\tarm\tdrift\tfloor\tcalled\texcess\thap1\thap2\tmargin\n' > "$OUT/cohort.tsv"
fails=0

for P in $DONORS; do
  D="$W/pair$P"
  [ -f "$D/truth_paths.tsv" ] || { echo "  pair$P: no truth_paths.tsv, skipped"; continue; }
  T1=$(awk -F'\t' 'NR==2{print $3}' "$D/truth_paths.tsv")
  T2=$(awk -F'\t' 'NR==3{print $3}' "$D/truth_paths.tsv")
  for i in 1 2; do
    python3 -c "
b=[l.strip() for l in open('$D/truth$i.fa') if l[0]!='>']
open('$OUT/p$P.t$i.fa','w').write('>t\n'+''.join(b).upper()+'\n')"
  done

  # ---- STALE-OUTPUT GUARD ----------------------------------------------------------------------
  # Any previous run's files for this donor are removed. A partially written donor directory once
  # produced a NOFLOOR that looked like a real measurement.
  rm -f "$OUT/p$P".*.hap_pairs.tsv "$OUT/p$P".*.called.fa "$OUT/p$P.dmp.scored_sequences".*

  # ---- DRIFT CHECK, before any accuracy number ------------------------------------------------
  # Spell the truth haplotypes from THIS graph and compare to the FASTA the reads were drawn from.
  # Anything but 0 means the substrate is wrong and every number below it would be charging a
  # representation difference to the genotyper.
  "$BIN" genotype-frag -i "$G" -b "$PFX" -o "$OUT/p$P.dmp" \
    --dump-scored-sequences "$OUT/p$P.dmp" --exclude-haplotypes "$T1,$T2" -q >/dev/null 2>&1
  DRIFT=$(python3 - "$OUT/p$P.dmp.scored_sequences.fa" "$OUT/p$P" "$T1" "$T2" "$BIN" <<'PY'
import sys, subprocess
fa,pfx,n1,n2,BIN=sys.argv[1:6]
want={n1:1,n2:2}; n=None;buf=[];out={}
def flush():
    if n in want: out[want[n]]="".join(buf)
for l in open(fa):
    if l[0]=='>': flush(); n=l[1:].split()[0]; buf=[]
    else: buf.append(l.strip())
flush()
tot=0
for i in (1,2):
    if i not in out: print("NOSPELL"); sys.exit()
    open(f"{pfx}.g{i}.fa","w").write(">g\n"+out[i].upper()+"\n")
    r=subprocess.run([BIN,"genotype-frag","--exact-distance",f"{pfx}.g{i}.fa",f"{pfx}.t{i}.fa",
                      "--distance-band","4096"],capture_output=True,text=True)
    v=r.stdout.strip()
    if not v or v.startswith('>'): print("DRIFTED"); sys.exit()
    tot+=int(v)
print(tot)
PY
)
  if [ "$DRIFT" != "0" ]; then
    echo "  pair$P: DRIFT CHECK FAILED ($DRIFT) -- wrong substrate, donor skipped"
    fails=$((fails+1)); continue
  fi

  # ---- certified floor: every remaining panel haplotype, exact, banded -------------------------
  # The dump carries BOTH the panel and held_out groups, and held_out IS the truth. Scoring it makes
  # every floor 0 -- measured, and uniform enough to look right. Only the panel group is a candidate.
  FLOOR=$(python3 - "$OUT/p$P.dmp.scored_sequences.fa" "$OUT/p$P" "$BIN" "$BAND" <<'PY'
import sys, subprocess
fa,pfx,BIN,band=sys.argv[1:5]
n=None;grp=None;buf=[];names=[];seqs=[]
def flush():
    if n and grp=="panel": names.append(n); seqs.append("".join(buf))
for l in open(fa):
    if l[0]=='>':
        flush(); p=l[1:].split(); n=p[0]; grp=p[1] if len(p)>1 else "panel"; buf=[]
    else: buf.append(l.strip())
flush()
best=[10**9,10**9]
for s in seqs:
    open(f"{pfx}.c.fa","w").write(">c\n"+s.upper()+"\n")
    for h in (0,1):
        r=subprocess.run([BIN,"genotype-frag","--exact-distance",f"{pfx}.c.fa",f"{pfx}.t{h+1}.fa",
                          "--distance-band",band],capture_output=True,text=True)
        v=r.stdout.strip()
        if v and not v.startswith('>'): best[h]=min(best[h],int(v))
print("NOFLOOR" if 10**9 in best else best[0]+best[1])
PY
)
  rm -f "$OUT/p$P.dmp.scored_sequences.fa" "$OUT/p$P.dmp.scored_sequences.tsv"

  for ARM in $ARMS; do
    EXW=()
    case "$ARM" in
      base)      EXW=() ;;
      bandfloor) EXW=(--band-floor) ;;
      *) echo "  unknown arm $ARM" >&2; continue ;;
    esac
    "$BIN" genotype-frag -i "$G" -b "$PFX" -o "$OUT/p$P.$ARM" \
      -R "$D/reads_1.fq.gz" -R "$D/reads_2.fq.gz" \
      --exclude-haplotypes "$T1,$T2" \
      --haplotype-mode --hamming-emission --max-divergence 0.05 \
      --fragment-len 350 --fragment-sd 50 --error-rate 0.001 --haploid-depth 0.05 \
      --max-anchor-occ 32 --placement-topk 8 \
      --top-pairs 20 -t "${THREADS:-4}" ${EXW[@]+"${EXW[@]}"} -q >/dev/null 2>&1
    # Spell BY NAME: allele indices are only valid for one catalogue, and spelling across a
    # mismatch is silent. Path names are stable.
    "$BIN" genotype-frag -i "$G" -b "$PFX" -o "$OUT/p$P.$ARM.sp" \
      --spell-pair "$OUT/p$P.$ARM.hap_pairs.tsv" -q >/dev/null 2>&1
    CALLED=$(python3 - "$OUT/p$P.$ARM.sp.called.fa" "$OUT/p$P.t1.fa" "$OUT/p$P.t2.fa" \
                       "$OUT/p$P.$ARM" "$BIN" <<'PY'
import sys, subprocess
def recs(p):
    out=[];n=None;b=[]
    for l in open(p):
        if l[0]=='>':
            if n is not None: out.append("".join(b))
            n=1;b=[]
        else: b.append(l.strip())
    if n is not None: out.append("".join(b))
    return out
try: c=recs(sys.argv[1])
except Exception: print("NA"); sys.exit()
if len(c)<2: print("NA"); sys.exit()
pfx,BIN=sys.argv[4],sys.argv[5]
def d(a,b):
    open(f"{pfx}.q.fa","w").write(">q\n"+a.upper()+"\n")
    r=subprocess.run([BIN,"genotype-frag","--exact-distance",f"{pfx}.q.fa",b,
                      "--distance-band","262144"],capture_output=True,text=True)
    v=r.stdout.strip()
    return int(v) if v and not v.startswith('>') else 10**9
print(min(d(c[0],sys.argv[2])+d(c[1],sys.argv[3]), d(c[0],sys.argv[3])+d(c[1],sys.argv[2])))
PY
)
    H1=$(sed -n 2p "$OUT/p$P.$ARM.hap_pairs.tsv" 2>/dev/null | cut -f2)
    H2=$(sed -n 2p "$OUT/p$P.$ARM.hap_pairs.tsv" 2>/dev/null | cut -f3)
    MG=$(awk -F'\t' 'NR==2{print $1}' "$OUT/p$P.$ARM.equivalence.tsv" 2>/dev/null)
    EX="NA"
    case "$FLOOR$CALLED" in *NA*|*NOFLOOR*) ;; *) EX=$((CALLED-FLOOR)) ;; esac
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$P" "$ARM" "$DRIFT" "$FLOOR" "$CALLED" "$EX" "${H1:-NA}" "${H2:-NA}" "${MG:-NA}" >> "$OUT/cohort.tsv"
    echo "  pair$P $ARM: drift $DRIFT  floor $FLOOR  called $CALLED  excess $EX"
  done
done

echo
echo "binary $BIN_MD5"
awk -F'\t' 'NR>1 && $6!="NA"{n++; s+=$6; if($5<1000) g++} END{
  if(n) printf "%d donors: median-free mean excess %.0f; within 1000 edits %d/%d\n", n, s/n, g+0, n}' "$OUT/cohort.tsv"
[ "$fails" -eq 0 ] || echo "$fails donor(s) failed the drift check and were skipped"
exit 0
