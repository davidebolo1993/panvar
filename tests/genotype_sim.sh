#!/usr/bin/env bash
# genotype_sim.sh - simulate a diploid from two panel haplotypes and score the genotyper against it.
#
# Reads are simulated from the pre-panphorte `bubble` stage, where paths still spell their haplotypes
# exactly; simulating from a folded graph would test the genotyper against the collapsed sequence
# rather than the real one.
#
#   genotype_sim.sh <locus> [n_pairs] [depth] [error]
#
# Env: LOO=1 for leave-one-out (the sample's own haplotypes dropped from the panel)
# Env: PANVAR_BIN (default build/panvar), PYTHON, OUT (default a scratch dir), SEED,
#      GENOTYPE_EXTRA (extra flags passed through to `panvar genotype`)
set -uo pipefail

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'
  exit 0
fi

LOCUS="${1:?usage: genotype_sim.sh <locus> [n_pairs] [depth] [error]}"
PAIRS="${2:-5}"
DEPTH="${3:-30}"
ERR="${4:-0.001}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"
BIN="${PANVAR_BIN:-$REPO/build/panvar}"
PY="${PYTHON:-python3}"
OUT="${OUT:-${TMPDIR:-/tmp}/panvar_genotype_sim}/$LOCUS"
SEED="${SEED:-42}"

G="$REPO/results/real_data/$LOCUS/bubble/bubble.sorted.gfa"
PFX="$REPO/results/real_data/$LOCUS/bubble/bubble"
[[ -s "$G" ]] || { echo "no bubble-stage graph for $LOCUS ($G); run scripts/regen_results.sh first"; exit 1; }
mkdir -p "$OUT"

REF="$(gzcat "$REPO/tests/real_data/$LOCUS.gfa.gz" 2>/dev/null | awk -F'\t' \
  '($1=="P"||$1=="W"){n=$2; if(n~/[Gg][Rr][Cc]h38/){print n;exit} if(!f)f=n} END{if(f&&!d)print f}' | head -1)"

# bash 3.2 (macOS default) has no mapfile
NAMES=()
while IFS= read -r line; do NAMES+=("$line"); done < <("$PY" "$REPO/scripts/spell_paths.py" -i "$G" --list | cut -f1)
N=${#NAMES[@]}
echo "locus $LOCUS: $N panel haplotypes; $PAIRS pairs at ${DEPTH}x, error $ERR"

# Per-call records for the calibration plot: one row per scored block across every pair.
CALLS="${CALLS:-$OUT/calls.tsv}"
# Truncate unless the caller is deliberately accumulating across loci in one sweep (CALLS_APPEND=1).
# Appending by default silently mixed rows from earlier runs into the calibration plot, so a figure
# could be built from two different builds of the binary.
if [[ "${CALLS_APPEND:-0}" == "0" || ! -s "$CALLS" ]]; then
  printf 'locus\tdepth\terror\tloo\tpair\tblock_kind\tn_alleles\tn_markers\tgq\texplained\tcorrect\tfilter\n' > "$CALLS"
fi

exact=0; partial=0; wrong=0; total=0; bex=0; btot=0
for ((p=0; p<PAIRS; p++)); do
  i=$(( (SEED + p * 7919) % N ))
  j=$(( (SEED + p * 104729 + 13) % N ))
  [[ "$i" == "$j" ]] && j=$(( (j + 1) % N ))
  H1="${NAMES[$i]}"; H2="${NAMES[$j]}"
  rm -rf "$OUT/fa"; "$PY" "$REPO/scripts/spell_paths.py" -i "$G" -o "$OUT/fa" --paths "$H1,$H2" >/dev/null
  rm -f "$OUT/r_1.fq" "$OUT/r_2.fq"; hidx=0
  for f in "$OUT"/fa/*.fa; do
    L=$(awk 'NR>1{n+=length($0)} END{print n}' "$f")
    NR_=$(( DEPTH * L / 2 / 300 ))          # per haplotype: half the diploid depth
    # A distinct seed per homolog. Sharing one across both makes their coverage fluctuations identical,
    # which is the one thing a diploid simulation must not do -- it hides exactly the per-haplotype
    # depth noise that decides het against hom.
    hseed=$((SEED + p*97 + hidx)); hidx=$((hidx+1))
    wgsim -N "$NR_" -1 150 -2 150 -d 350 -s 50 -e "$ERR" -r 0 -R 0 -X 0 -S "$hseed" \
      "$f" "$OUT/a_1.fq" "$OUT/a_2.fq" >/dev/null 2>&1
    cat "$OUT/a_1.fq" >> "$OUT/r_1.fq"; cat "$OUT/a_2.fq" >> "$OUT/r_2.fq"
  done
  gzip -f "$OUT/r_1.fq" "$OUT/r_2.fq"
  line=$("$BIN" genotype -i "$G" --bubble-prefix-in "$PFX" -r "$REF" -o "$OUT/gt" \
    -R "$OUT/r_1.fq.gz" -R "$OUT/r_2.fq.gz" --truth-haplotypes "$H1,$H2" \
    $( [ "${LOO:-0}" != "0" ] && printf -- "--exclude-haplotypes %s,%s" "$H1" "$H2" ) ${GENOTYPE_EXTRA:-} 2>&1 | grep -E "truth check|unrepresentable")
  unrep=$(sed -nE 's/.*unrepresentable blocks[^0-9]*([0-9]+)\/.*/\1/p' <<<"$line" | head -1)
  line=$(grep "truth check" <<<"$line")
  e=$(sed -E 's/.*check: ([0-9]+)\/([0-9]+).*/\1/' <<<"$line")
  t=$(sed -E 's/.*check: ([0-9]+)\/([0-9]+).*/\2/' <<<"$line")
  pa=$(sed -E 's/.*\(([0-9]+) one allele right, ([0-9]+) both wrong\).*/\1/' <<<"$line")
  wr=$(sed -E 's/.*\(([0-9]+) one allele right, ([0-9]+) both wrong\).*/\2/' <<<"$line")
  be=$(sed -E 's/.*bubble blocks ([0-9]+)\/([0-9]+).*/\1/' <<<"$line")
  bt=$(sed -E 's/.*bubble blocks ([0-9]+)\/([0-9]+).*/\2/' <<<"$line")
  # Columns are looked up BY NAME. Fixed indices silently drift the moment a column is added, and did:
  # truth moved from 17-18 to 25-26 and filter from 16 to 24 as length estimates were inserted, so every
  # calibration figure this script produced after that was scored against called_bp and mass_bp.
  awk -F'\t' -v L="$LOCUS" -v D="$DEPTH" -v E="$ERR" -v LO="${LOO:-0}" -v P="$p" '
    NR==1 { for (i=1;i<=NF;i++) c[$i]=i; next }
    $(c["truth1"])>=0 {
      lo=($(c["allele1"])<$(c["allele2"])?$(c["allele1"]):$(c["allele2"]));
      hi=($(c["allele1"])<$(c["allele2"])?$(c["allele2"]):$(c["allele1"]));
      tl=($(c["truth1"])<$(c["truth2"])?$(c["truth1"]):$(c["truth2"]));
      th=($(c["truth1"])<$(c["truth2"])?$(c["truth2"]):$(c["truth1"]));
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\n", L,D,E,LO,P,
             $(c["block_kind"]),$(c["n_alleles"]),$(c["n_markers"]),$(c["gq"]),$(c["explained"]),
             (lo==tl&&hi==th)?1:0,$(c["filter"])
    }' "$OUT/gt.genotypes.tsv" >> "$CALLS" 2>/dev/null
  printf "  pair %d: %s/%s blocks exact, bubbles %s/%s%s\n" "$p" "$e" "$t" "$be" "$bt" "${unrep:+, $unrep unrepresentable}"
  exact=$((exact+e)); total=$((total+t)); partial=$((partial+pa)); wrong=$((wrong+wr))
  bex=$((bex+be)); btot=$((btot+bt))
done
echo "TOTAL $LOCUS: $exact/$total blocks exact ($partial one-allele, $wrong both-wrong); bubbles $bex/$btot"
# Accuracy over ALL scored blocks hides the filter doing its job: a call the model declined counts the
# same as one it got wrong. Report the PASS subset and its call rate alongside.
awk -F'\t' -v L="$LOCUS" 'NR>1 && $1==L {
    n++; c+=$11;
    if ($12=="PASS") { pn++; pc+=$11 }
  } END {
    if (n) printf "  all calls: %d/%d exact (%.1f%%); PASS only: %d/%d exact (%.1f%%), call rate %.1f%%\n",
      c, n, 100*c/n, pc, pn, (pn?100*pc/pn:0), 100*pn/n
  }' "$CALLS"
