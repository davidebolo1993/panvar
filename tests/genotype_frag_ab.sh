#!/usr/bin/env bash
# genotype_frag_ab.sh - run the production marker genotyper and the fragment prototype over the SAME
# simulated reads, and report the two side by side, per stratum.
#
#   genotype_frag_ab.sh <locus> [n_donors] [depth] [error]
#
# Both arms see identical reads, an identical panel and an identical leave-one-out exclusion, so a
# difference between them is the evidence model and nothing else. Every figure is PAIRED per donor
# per block: an aggregate over donors hides which donors moved, and this project has twice drawn a
# conclusion from a comparison that was selective in a way nobody noticed at the time.
#
# Env: LOO=1 (default) leave-one-out; LOO=0 leave-zero-out, where an exact call is reachable and a
#      failure is an implementation defect rather than a panel limit -- run it first, always.
#      ARMS="prod hap block" which arms to run (default "prod hap")
#      PANVAR_BIN, SEED, OUT, FRAG_EXTRA (extra flags for genotype-frag)
set -uo pipefail

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  sed -n '2,16p' "$0" | sed 's/^# \{0,1\}//'
  exit 0
fi

LOCUS="${1:?usage: genotype_frag_ab.sh <locus> [n_donors] [depth] [error]}"
DONORS_N="${2:-5}"
DEPTH="${3:-30}"
ERR="${4:-0.001}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
BIN="${PANVAR_BIN:-$REPO/build_genotype_review/panvar}"
PY="${PYTHON:-python3}"
OUT="${OUT:-${TMPDIR:-/tmp}/panvar_frag_ab}/$LOCUS"
SEED="${SEED:-42}"
LOO="${LOO:-1}"
ARMS="${ARMS:-prod hap}"

G="$REPO/results/real_data/$LOCUS/bubble/bubble.sorted.gfa"
PFX="$REPO/results/real_data/$LOCUS/bubble/bubble"
[[ -s "$G" ]] || { echo "no bubble-stage graph for $LOCUS ($G)"; exit 1; }
[[ -x "$BIN" ]] || { echo "no binary at $BIN (build with -DPANVAR_ENABLE_EXPERIMENTAL_GENOTYPE=ON)"; exit 1; }
mkdir -p "$OUT"

REF="$(gunzip -c "$REPO/tests/real_data/$LOCUS.gfa.gz" 2>/dev/null | awk -F'\t' \
  '($1=="P"||$1=="W"){n=$2; if(n~/[Gg][Rr][Cc]h38/){print n;exit} if(!f)f=n} END{if(f&&!d)print f}' | head -1)"

NAMES=(); while IFS= read -r l; do NAMES+=("$l"); done < <("$PY" "$REPO/scripts/spell_paths.py" -i "$G" --list | cut -f1)
DONORS=(); while IFS= read -r l; do DONORS+=("$l"); done < <(
  printf '%s\n' "${NAMES[@]}" | awk -F'#' '{c[$1]++; if(c[$1]==1) first[$1]=$0; else if(c[$1]==2) second[$1]=$0}
    END{for(s in c) if(c[s]>=2) print first[s]"\t"second[s]}' | sort)
ND=${#DONORS[@]}
(( ND > 0 )) || { echo "no donor has two haplotypes in $G"; exit 1; }

# One row per (donor, block, arm). Everything below is computed from this file, so any claim in the
# summary can be traced to the blocks it came from.
ROWS="$OUT/blocks.tsv"
printf 'locus\tloo\tdonor\tarm\tblock\tkind\tn_alleles\trepresentable\texact\n' > "$ROWS"

echo "locus $LOCUS: ${#NAMES[@]} haplotypes, $ND donors; $DONORS_N donors at ${DEPTH}x, error $ERR, LOO=$LOO"
echo "arms: $ARMS"

for ((p=0; p<DONORS_N && p<ND; p++)); do
  d=$(( (SEED + p * 7919) % ND ))
  H1="${DONORS[$d]%%$'\t'*}"; H2="${DONORS[$d]##*$'\t'}"
  DONOR="${H1%%#*}"
  rm -rf "$OUT/fa"; "$PY" "$REPO/scripts/spell_paths.py" -i "$G" -o "$OUT/fa" --paths "$H1,$H2" >/dev/null
  rm -f "$OUT/r_1.fq" "$OUT/r_2.fq"; hidx=0
  for f in "$OUT"/fa/*.fa; do
    L=$(awk 'NR>1{n+=length($0)} END{print n}' "$f")
    NRD=$(( DEPTH * L / 2 / 300 ))
    wgsim -N "$NRD" -1 150 -2 150 -d 350 -s 50 -e "$ERR" -r 0 -R 0 -X 0 -S $((SEED + p*97 + hidx)) \
      "$f" "$OUT/a_1.fq" "$OUT/a_2.fq" >/dev/null 2>&1
    cat "$OUT/a_1.fq" >> "$OUT/r_1.fq"; cat "$OUT/a_2.fq" >> "$OUT/r_2.fq"; hidx=$((hidx+1))
  done
  gzip -f "$OUT/r_1.fq" "$OUT/r_2.fq"
  EXCL=""; [[ "$LOO" != "0" ]] && EXCL="--exclude-haplotypes $H1,$H2"

  for arm in $ARMS; do
    case "$arm" in
      prod)
        "$BIN" genotype -i "$G" --bubble-prefix-in "$PFX" -r "$REF" -o "$OUT/prod" \
          -R "$OUT/r_1.fq.gz" -R "$OUT/r_2.fq.gz" --truth-haplotypes "$H1,$H2" $EXCL \
          --dump-haplotype-alleles "$OUT/alleles.tsv" -q >/dev/null 2>&1
        awk -F'\t' -v L="$LOCUS" -v LO="$LOO" -v D="$DONOR" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
          { rep = ($(c["truth1"])>=0 && $(c["truth2"])>=0) ? 1 : 0
            ex="NA"
            if (rep) {
              lo=($(c["allele1"])<$(c["allele2"])?$(c["allele1"]):$(c["allele2"]))
              hi=($(c["allele1"])<$(c["allele2"])?$(c["allele2"]):$(c["allele1"]))
              tl=($(c["truth1"])<$(c["truth2"])?$(c["truth1"]):$(c["truth2"]))
              th=($(c["truth1"])<$(c["truth2"])?$(c["truth2"]):$(c["truth1"]))
              ex=(lo==tl&&hi==th)?1:0 }
            printf "%s\t%s\t%s\tprod\t%s\t%s\t%s\t%s\t%s\n", L,LO,D,$(c["block_index"]),$(c["block_kind"]),$(c["n_alleles"]),rep,ex }' \
          "$OUT/prod.genotypes.tsv" >> "$ROWS" 2>/dev/null
        ;;
      hap)
        "$BIN" genotype-frag -i "$G" -b "$PFX" -o "$OUT/hap" -R "$OUT/r_1.fq.gz" -R "$OUT/r_2.fq.gz" \
          --haplotype-mode --truth-haplotypes "$H1,$H2" $EXCL ${FRAG_EXTRA:-} -q >/dev/null 2>&1
        awk -F'\t' -v L="$LOCUS" -v LO="$LOO" -v D="$DONOR" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
          { printf "%s\t%s\t%s\thap\t%s\t%s\t%s\t%s\t%s\n", L,LO,D,$(c["block"]),$(c["kind"]),$(c["n_alleles"]),$(c["representable"]),$(c["exact"]) }' \
          "$OUT/hap.hap_blocks.tsv" >> "$ROWS" 2>/dev/null
        ;;
      ceiling)
        # The best ANY pair of complete panel haplotypes could do. Not the panel's ceiling: a mosaic
        # model may take a different haplotype at every block. The gap between this and production is
        # what a mosaic layer is worth, so it answers "is the evidence model losing, or is the
        # missing mosaic layer losing?" with a number instead of an argument.
        [[ -s "$OUT/alleles.tsv" && -s "$OUT/prod.genotypes.tsv" ]] || continue
        "$PY" "$REPO/scripts/genotype_pair_ceiling.py" "$OUT/alleles.tsv" "$OUT/prod.genotypes.tsv" \
          "$LOCUS" "$LOO" "$DONOR" >> "$ROWS"
        ;;
      block)
        "$BIN" genotype-frag -i "$G" -b "$PFX" -o "$OUT/blk" -R "$OUT/r_1.fq.gz" -R "$OUT/r_2.fq.gz" \
          --all-blocks --truth-haplotypes "$H1,$H2" $EXCL ${FRAG_EXTRA:-} -q >/dev/null 2>&1
        awk -F'\t' -v L="$LOCUS" -v LO="$LOO" -v D="$DONOR" 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
          { rep = ($(c["truth_rank"])!=-1) ? 1 : 0
            printf "%s\t%s\t%s\tblock\t%s\t%s\t%s\t%s\t%s\n", L,LO,D,$(c["block"]),$(c["kind"]),$(c["n_alleles"]),rep,$(c["exact"]) }' \
          "$OUT/blk.frag_blocks.tsv" >> "$ROWS" 2>/dev/null
        ;;
    esac
  done
  echo "  donor $DONOR done"
done

echo
"$PY" - "$ROWS" <<'PYEOF'
import sys, collections
rows=[l.rstrip('\n').split('\t') for l in open(sys.argv[1])][1:]
if not rows: sys.exit("no rows")
arms=sorted({r[3] for r in rows}, key=lambda a:{'prod':0,'hap':1,'block':2}.get(a,9))
# Paired on (donor, block): only blocks every arm scored and where the truth is representable in
# BOTH. An arm that silently skips a block would otherwise get credit for an easier denominator.
by=collections.defaultdict(dict)
for r in rows:
    by[(r[2],r[4])][r[3]]=(r[5], r[7], r[8])
common=[k for k,v in by.items() if all(a in v for a in arms) and all(v[a][1]=='1' for a in arms)]
print(f"paired blocks with a representable truth in every arm: {len(common)}")
def rate(arm, pred=lambda kind: True):
    sel=[k for k in common if pred(by[k][arm][0])]
    ok=sum(1 for k in sel if by[k][arm][2]=='1')
    return ok, len(sel)
print(f"\n{'stratum':<22}" + "".join(f"{a:>16}" for a in arms))
for label, pred in [("all blocks", lambda k: True),
                    ("bubble", lambda k: k=="bubble"),
                    ("backbone", lambda k: k=="backbone")]:
    line=f"{label:<22}"
    for a in arms:
        ok,n=rate(a,pred)
        line += f"{(str(ok)+'/'+str(n)) if n else '-':>10}" + (f"{100*ok/n:>5.0f}%" if n else "     ")
    print(line)
if len(arms)>1:
    base=arms[0]
    for a in arms[1:]:
        b=sum(1 for k in common if by[k][a][2]=='1' and by[k][base][2]!='1')
        w=sum(1 for k in common if by[k][a][2]!='1' and by[k][base][2]=='1')
        print(f"\n{a} vs {base}: {b} blocks fixed, {w} blocks broken")
        # Per donor, so a net gain carried by one donor cannot pass as a general one.
        pd=collections.Counter()
        for k in common:
            if by[k][a][2]=='1' and by[k][base][2]!='1': pd[k[0]]+=1
            elif by[k][a][2]!='1' and by[k][base][2]=='1': pd[k[0]]-=1
        print("  net per donor: " + ", ".join(f"{d}{'+' if v>=0 else ''}{v}" for d,v in sorted(pd.items())))
PYEOF
echo
echo "per-block rows: $ROWS"
