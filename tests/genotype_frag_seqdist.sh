#!/usr/bin/env bash
# genotype_frag_seqdist.sh - does the prototype's pair beat the label-ceiling pair IN SEQUENCE?
#
#   genotype_frag_seqdist.sh <locus> [n_donors] [depth] [error]
#
# The non-mosaic ceiling maximises exact-label block agreement. A fragment likelihood optimises how
# well a sequence explains the reads. Under leave-one-out the sample is off-panel, so those are
# different objectives and the label ceiling is not automatically a target the likelihood should
# reach. This measures the other objective for the same donors: divergence from the sample's own two
# haplotypes, for the ceiling pair, the prototype's pair and production's own best panel pair.
set -uo pipefail

LOCUS="${1:?usage: genotype_frag_seqdist.sh <locus> [n_donors] [depth] [error]}"
DONORS_N="${2:-10}"; DEPTH="${3:-30}"; ERR="${4:-0.001}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/.." && pwd)"
BIN="${PANVAR_BIN:-$REPO/build_genotype_review/panvar}"; PY="${PYTHON:-python3}"
OUT="${OUT:-${TMPDIR:-/tmp}/panvar_frag_seq}/$LOCUS"; SEED="${SEED:-42}"; MAXHAP="${MAXHAP:-48}"
G="$REPO/results/real_data/$LOCUS/bubble/bubble.sorted.gfa"
PFX="$REPO/results/real_data/$LOCUS/bubble/bubble"
mkdir -p "$OUT"
REF="$("$PY" "$REPO/scripts/spell_paths.py" -i "$G" --list | cut -f1 | grep -i grch38 | head -1)"

NAMES=(); while IFS= read -r l; do NAMES+=("$l"); done < <("$PY" "$REPO/scripts/spell_paths.py" -i "$G" --list | cut -f1)
DONORS=(); while IFS= read -r l; do DONORS+=("$l"); done < <(
  printf '%s\n' "${NAMES[@]}" | awk -F'#' '{c[$1]++; if(c[$1]==1) first[$1]=$0; else if(c[$1]==2) second[$1]=$0}
    END{for(s in c) if(c[s]>=2) print first[s]"\t"second[s]}' | sort)
ND=${#DONORS[@]}

# The panel is the same for every donor except which two are held out, so it is spelled ONCE and
# filtered per donor. Spelling it per donor would dominate the runtime at the larger loci.
PANELDIR="$OUT/panel_all"
if [[ ! -d "$PANELDIR" ]]; then
  "$PY" "$REPO/scripts/spell_paths.py" -i "$G" -o "$PANELDIR" >/dev/null
fi

RES="$OUT/seqdist.tsv"
printf 'locus\tdonor\tarm\ttotal\tnm\tunaligned\ttotal1\ttotal2\taligned1\taligned2\n' > "$RES"
echo "locus $LOCUS: $DONORS_N donors, leave-one-out; divergence from the donor's own haplotypes"

for ((p=0; p<DONORS_N && p<ND; p++)); do
  d=$(( (SEED + p * 7919) % ND ))
  H1="${DONORS[$d]%%$'\t'*}"; H2="${DONORS[$d]##*$'\t'}"; DONOR="${H1%%#*}"
  rm -rf "$OUT/fa"; "$PY" "$REPO/scripts/spell_paths.py" -i "$G" -o "$OUT/fa" --paths "$H1,$H2" >/dev/null
  T1=$(ls "$OUT"/fa/*.fa | head -1); T2=$(ls "$OUT"/fa/*.fa | tail -1)
  rm -f "$OUT/r_1.fq" "$OUT/r_2.fq"; hidx=0
  for f in "$OUT"/fa/*.fa; do
    L=$(awk 'NR>1{n+=length($0)} END{print n}' "$f")
    wgsim -N $(( DEPTH * L / 2 / 300 )) -1 150 -2 150 -d 350 -s 50 -e "$ERR" -r 0 -R 0 -X 0 \
      -S $((SEED + p*97 + hidx)) "$f" "$OUT/a_1.fq" "$OUT/a_2.fq" >/dev/null 2>&1
    cat "$OUT/a_1.fq" >> "$OUT/r_1.fq"; cat "$OUT/a_2.fq" >> "$OUT/r_2.fq"; hidx=$((hidx+1))
  done
  gzip -f "$OUT/r_1.fq" "$OUT/r_2.fq"

  "$BIN" genotype -i "$G" --bubble-prefix-in "$PFX" -r "$REF" -o "$OUT/prod" \
    -R "$OUT/r_1.fq.gz" -R "$OUT/r_2.fq.gz" --truth-haplotypes "$H1,$H2" \
    --exclude-haplotypes "$H1,$H2" --dump-haplotype-alleles "$OUT/alleles.tsv" -q >/dev/null 2>&1
  [[ -s "$OUT/alleles.tsv" ]] || { echo "  $DONOR skipped"; continue; }
  "$PY" "$REPO/scripts/genotype_pair_ceiling.py" "$OUT/alleles.tsv" "$OUT/prod.genotypes.tsv" \
    "$LOCUS" 1 "$DONOR" "$OUT/optimum.tsv" > /dev/null
  read -r _ C1 C2 _ < <(sed -n 2p "$OUT/optimum.tsv")

  "$BIN" genotype-frag -i "$G" -b "$PFX" -o "$OUT/frag" -R "$OUT/r_1.fq.gz" -R "$OUT/r_2.fq.gz" \
    --haplotype-mode --exclude-haplotypes "$H1,$H2" --max-haplotypes "$MAXHAP" -q >/dev/null 2>&1
  read -r _ P1 P2 _ < <(sed -n 2p "$OUT/frag.hap_pairs.tsv")

  emit() {  # arm name, two haplotype names
    rm -rf "$OUT/cand"; "$PY" "$REPO/scripts/spell_paths.py" -i "$G" -o "$OUT/cand" --paths "$2,$3" >/dev/null
    local A B; A=$(ls "$OUT"/cand/*.fa | head -1); B=$(ls "$OUT"/cand/*.fa | tail -1)
    printf '%s\t%s\t' "$LOCUS" "$DONOR" >> "$RES"
    "$PY" "$REPO/scripts/genotype_pair_sequence_distance.py" "$T1" "$T2" "$A" "$B" "$1" >> "$RES"
  }
  # E*_panel: the best ANY panel pair could reconstruct. Separates panel limitation from inference
  # failure, which the exact-label ceiling could not.
  cat $(ls "$PANELDIR"/*.fa | grep -v "$DONOR") > "$OUT/panel_loo.fa"
  printf '%s\t%s\t' "$LOCUS" "$DONOR" >> "$RES"
  "$PY" "$REPO/scripts/genotype_panel_floor.py" "$T1" "$T2" "$OUT/panel_loo.fa" panel_floor \
    | awk -F'\t' '{printf "%s\t%s\t0\t0\t%s\t%s\t0\t0\n",$1,$2,$3,$4}' >> "$RES"

  emit ceiling "$C1" "$C2"
  emit prototype "$P1" "$P2"

  # Both callers, spelled from their own per-block call tables, so the comparison is on the objective
  # the release goal is stated in rather than on allele-index agreement.
  spell() {  # arm name, call table
    "$BIN" genotype-frag -i "$G" -b "$PFX" -o "$OUT/$1" --exclude-haplotypes "$H1,$H2" \
      --spell-calls "$2" -q >/dev/null 2>&1 || return
    [[ -s "$OUT/$1.called.fa" ]] || return
    "$PY" - "$OUT/$1.called.fa" "$OUT/$1" <<'PYS'
import sys
seqs={}; name=None
for line in open(sys.argv[1]):
    if line[0]=='>': name=line[1:].strip(); seqs[name]=[]
    else: seqs[name].append(line.strip())
for i,(n,v) in enumerate(seqs.items(), 1):
    with open(f"{sys.argv[2]}_{i}.fa","w") as fh:
        fh.write(f">{n}\n{''.join(v)}\n")
PYS
    printf '%s\t%s\t' "$LOCUS" "$DONOR" >> "$RES"
    "$PY" "$REPO/scripts/genotype_pair_sequence_distance.py" "$T1" "$T2" \
      "$OUT/${1}_1.fa" "$OUT/${1}_2.fa" "$1" >> "$RES"
  }
  spell production "$OUT/prod.genotypes.tsv"
  spell prototype_blocks "$OUT/frag.hap_blocks.tsv"
  echo "  $DONOR done"
done

echo
"$PY" - "$RES" <<'PYEOF'
import sys, collections
rows=[l.rstrip('\n').split('\t') for l in open(sys.argv[1])][1:]
by=collections.defaultdict(dict)
for r in rows: by[r[1]][r[2]]=int(r[3])
arms=[a for a in ('panel_floor','ceiling','prototype','production','prototype_blocks')
      if any(a in v for v in by.values())]
both=[d for d,v in by.items() if all(a in v for a in arms)]
if not both: sys.exit("no paired donors")
print("total distance from the donor's own two haplotypes (aligned NM + unaligned bp), lower is better")
unal={}
for r in rows: unal.setdefault(r[1],{})[r[2]]=int(r[5])
print(f"\n  {'donor':<10}" + "".join(f"{a:>19}" for a in arms))
for d in sorted(both):
    print(f"  {d:<10}" + "".join(f"{by[d][a]:>19d}" for a in arms))
med={a: sorted(by[d][a] for d in both)[len(both)//2] for a in arms}
print(f"  {'MEDIAN':<10}" + "".join(f"{med[a]:>19d}" for a in arms))
if 'panel_floor' in arms:
    print("\nEXCESS over the panel floor -- the part a better CALLER could still recover.")
    print("(the floor itself is the part only a better PANEL could)")
    print(f"\n  {'donor':<10}" + "".join(f"{a:>19}" for a in arms if a!='panel_floor'))
    for d in sorted(both):
        print(f"  {d:<10}" + "".join(f"{by[d][a]-by[d]['panel_floor']:>19d}"
                                     for a in arms if a!='panel_floor'))
    mex={a: sorted(by[d][a]-by[d]['panel_floor'] for d in both)[len(both)//2]
         for a in arms if a!='panel_floor'}
    print(f"  {'MEDIAN':<10}" + "".join(f"{mex[a]:>19d}" for a in arms if a!='panel_floor'))
if 'production' in arms:
    w=sum(1 for d in both if by[d]['prototype_blocks'] < by[d]['production'])
    print(f"\n  prototype closer than PRODUCTION in sequence: {w}/{len(both)} donors")
print()
win=sum(1 for d in both if by[d]['prototype'] < by[d]['ceiling'])
print(f"donors: {len(both)}")
ties=sum(1 for d in both if by[d]['prototype']==by[d]['ceiling'])
print(f"  prototype's pair is CLOSER IN SEQUENCE than the label-ceiling pair: {win}/{len(both)}"
      f"  ({ties} exact ties)")
mc=sorted(by[d]['ceiling'] for d in both); mp=sorted(by[d]['prototype'] for d in both)
print(f"  median total edits -- ceiling {mc[len(mc)//2]}, prototype {mp[len(mp)//2]}")
print()
print(f"  {'donor':<10} {'ceiling':>9} {'prototype':>10} {'delta':>8}   closer")
for d in sorted(both):
    c,p=by[d]['ceiling'],by[d]['prototype']
    print(f"  {d:<10} {c:>9d} {p:>10d} {p-c:>+8d}   {'prototype' if p<c else ('tie' if p==c else 'ceiling')}")
PYEOF
echo
echo "rows: $RES"
