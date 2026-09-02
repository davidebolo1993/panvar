#!/usr/bin/env bash
# C4/NA19240 block-7 contrast: signed delta between two NAMED pairs, over seeds and depths.
#
# Delta = score(truth-equivalent) - score(wrong competitor). Positive means the model prefers truth.
# Reporting only "truth won" censors the contrast -- three seeds showed 461.55, 88.32 and 0.00, where
# the 0.00 was simply "it won" and carried no magnitude. Both pairs are scored at every cell.
#
# --force-haplotypes ADDS to the shortlist rather than restricting it (measured), which is exactly
# what is needed here: it guarantees all three haplotypes are scored so both pairs appear in
# hap_pairs.tsv regardless of rank.
set -uo pipefail
D="${1:?out}"; REPO="${2:?repo}"; BIN="${3:?bin}"
DEPTHS="${DEPTHS:-10 20 30 60}"; SEEDS="${SEEDS:-1 2 3 4 5 6 7 8 9 10 11 12}"
mkdir -p "$D"
md5of() { if command -v md5 >/dev/null 2>&1; then md5 -q "$1"; else md5sum "$1" | cut -d' ' -f1; fi; }
{ echo "commit   $(cd "$REPO" && git rev-parse --short HEAD 2>/dev/null || echo unknown)"
  echo "binary   $(md5of "$BIN")"
  echo "depths   $DEPTHS"
  echo "seeds    $SEEDS"
  echo "date     $(date -u +%Y-%m-%dT%H:%M:%SZ)"; } > "$D/PROVENANCE.txt"
cat "$D/PROVENANCE.txt"
T1="NA19240#1#haplotype1-0000009:31981961-32208331"   # truth homologue 1
TE="NA19238#1#haplotype1-0000026:31825469-32078164"   # sequence-identical twin of truth homologue 2
WC="HG00171#1#haplotype1-0000033:140888062-141114404" # the competitor that keeps winning
G="$REPO/results/real_data/c4/bubble/bubble.sorted.gfa"
P="$REPO/results/real_data/c4/bubble/bubble"
TRUTH="$D/truth.fa"
# Distinct record names. Both were ">s", which wgsim tolerates but makes any per-record accounting
# -- including "which homologue did this fragment come from" -- impossible to reconstruct later.
awk '/^>/{ printf ">hap%d\n", ++i; next } {print}' "$TRUTH" > "$TRUTH.named" && mv "$TRUTH.named" "$TRUTH"
L=$(awk '/^>/{next}{n+=length($0)}END{print n}' "$TRUTH")
printf 'depth\tseed\tn_pairs_requested\tn_fragments_loaded\tscore_truth_equiv\tscore_competitor\tdelta\ttruth_wins\twinner\n' > "$D/contrast.tsv"
cells=0; expected=0
for COV in $DEPTHS; do
  N=$(( COV * L / 600 ))
  for S in $SEEDS; do
    expected=$((expected+1))
    wgsim -N "$N" -1 150 -2 150 -d 350 -s 50 -e 0.001 -r 0 -R 0 -X 0 -S "$S" \
      "$TRUTH" "$D/r1.fq" "$D/r2.fq" >/dev/null 2>&1 || continue
    "$BIN" genotype-frag -i "$G" -b "$P" -o "$D/run" -R "$D/r1.fq" -R "$D/r2.fq" \
      --haplotype-mode --hamming-emission --max-divergence 0.05 --fragment-len 350 \
      --fragment-sd 50 --error-rate 0.001 --max-anchor-occ 32 \
      --placement-topk 8 --top-pairs 4000 --force-haplotypes "$T1,$TE,$WC" 2>"$D/run.err" >/dev/null
    # LOADED, not requested. wgsim writes what it writes and the caller drops what it drops; a cell
    # that silently loaded half its fragments would look like a weak signal rather than a bad cell.
    NLOAD=$(sed -nE 's/.*-> ([0-9]+) fragments.*/\1/p' "$D/run.err" | head -1)
    python3 - "$D/run.hap_pairs.tsv" "$T1" "$TE" "$WC" "$COV" "$S" "$N" "${NLOAD:-NA}" \
      >> "$D/contrast.tsv" <<'PY'
import sys
f,t1,te,wc,cov,seed,n,nload=sys.argv[1:9]
rows=[l.rstrip('\n').split('\t') for l in open(f)][1:]
def find(a,b):
    for r in rows:
        if {r[1],r[2]}=={a,b}: return float(r[3])
    return None
st=find(t1,te); sc=find(t1,wc)
win=rows[0][1].split('#')[0]+"#"+rows[0][1].split('#')[1]+"/"+rows[0][2].split('#')[0]+"#"+rows[0][2].split('#')[1] if rows else "NA"
if st is None or sc is None:
    print("%s\t%s\t%s\t%s\tNA\tNA\tNA\tNA\t%s"%(cov,seed,n,nload,win))
else:
    d=st-sc
    print("%s\t%s\t%s\t%s\t%.4f\t%.4f\t%.4f\t%s\t%s"%(cov,seed,n,nload,st,sc,d,"yes" if d>0 else "no",win))
PY
  done
done

# EVERY CELL MUST SUCCEED. The first version continued past a failed arm, so a partial ladder would
# have been averaged as if complete -- and a missing 60x cell is exactly the one that changes the
# conclusion from "insensitive" to "insensitive at 30x".
cells=$(( $(wc -l < "$D/contrast.tsv") - 1 ))
bad=$(awk -F'\t' 'NR>1 && ($7=="NA"||$4=="NA"){n++} END{print n+0}' "$D/contrast.tsv")
echo
awk -F'\t' 'NR>1 && $7!="NA"{n[$1]++; s[$1]+=$7; ss[$1]+=$7*$7; if($7>0) w[$1]++}
  END{printf "  %-6s %10s %8s %10s\n","depth","meanD","wins","sd";
      for(d in n){m=s[d]/n[d]; v=ss[d]/n[d]-m*m
        printf "  %-6s %10.0f %5d/%-2d %10.0f\n", d, m, w[d]+0, n[d], (v>0?sqrt(v):0)}}' \
  "$D/contrast.tsv" | sort -n
echo
if [ "$cells" != "$expected" ] || [ "$bad" != 0 ]; then
  echo "FAIL: $cells of $expected cells, $bad unusable -- an incomplete ladder must not be averaged"
  exit 1
fi
echo "c4 block-7 contrast: $cells/$expected cells, all usable"
