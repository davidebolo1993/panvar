#!/usr/bin/env bash
# Shared helpers for the per-gene drivers (scripts/genes/<gene>.sh); sourced, not run directly.
# Each gene runs the data pipeline (bubble -> inspect -> panphorte -> inspect -> call -> describe);
# plot commands sit commented at the bottom of each. Env: PANVAR_BIN, PYTHON, RSCRIPT, THREADS, GTF.

set -uo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
BIN="${PANVAR_BIN:-$REPO/build/panvar}"
PY="${PYTHON:-python3}"
RS="${RSCRIPT:-Rscript}"
DATA="$REPO/tests/real_data"
OUT="$REPO/results/real_data"
THREADS="${THREADS:-0}"
GTF="${GTF:-$DATA/Homo_sapiens.GRCh38.116.gtf.gz}"
GTFOPT=(); [[ -f "$GTF" ]] && GTFOPT=(--gtf "$GTF")

# pick a reference path: prefer GRCh38, then the first path in the GFA
ref_of() {
  gzcat "$1" 2>/dev/null | awk -F'\t' '($1=="P"||$1=="W"){n=$2;
    if(n~/[Gg][Rr][Cc]h38/){print n;exit} if(!first)first=n} END{if(first)print first}' | head -1
}

# largest-DUP bubble id from a region VCF (the copy-number module to inspect/plot)
main_dup_bubble() {
  awk -F'\t' '/SVTYPE=DUP/{ sv=$8; sub(/.*SVLEN=/,"",sv); sub(/;.*/,"",sv); if(sv<0)sv=-sv;
      b=$8; sub(/.*BUBBLE_ID=/,"",b); sub(/;.*/,"",b);
      if(sv+0>m){m=sv+0;bb=b}} END{print bb}' "$1"
}

# largest repeat bubble id from a bubble.bubbles.csv (column 8 = size-like metric)
main_bubble_csv() { awk -F',' 'NR>1 && $8+0>m{m=$8+0; b=$1} END{print b}' "$1"; }

# run_gene_data <region> <gfa.gz> <panphorte_min_sim> <call_graph: bubble|panphorte> <call_extra> [cluster_sim]
run_gene_data() {
  local region="$1" gfa="$2" pan_sim="$3" call_graph="$4" call_extra="$5" clu="${6:-0.97}"
  local d="$OUT/$region"
  echo "############ $region (data only) ############"
  mkdir -p "$d"/{bubble,panphorte,call,describe,inspect}
  local ref; ref="$(ref_of "$gfa")"
  echo "[$region] reference: $ref ; call graph: $call_graph ; panphorte sim: $pan_sim ; inspect cluster: $clu"

  # 1) bubble (reference-oriented sort + cactus snarls)
  "$BIN" bubble -i "$gfa" -o "$d/bubble/bubble" -r "$ref" "${GTFOPT[@]}" --quiet || return 1
  local sgfa="$d/bubble/bubble.sorted.gfa"

  # 2) inspect the largest repeat bubble on the BUBBLE graph (pre-panphorte), with walk clustering
  local bbub; bbub="$(main_bubble_csv "$d/bubble/bubble.bubbles.csv")"
  [[ -n "$bbub" ]] && "$BIN" inspect -i "$sgfa" --bubble-prefix-in "$d/bubble/bubble" \
      --bubble-id "$bbub" --cluster --cluster-similarity "$clu" -o "$d/inspect/inspect_bubblestage" --quiet || true

  # 3) panphorte (normalized + re-sorted)
  "$BIN" panphorte -i "$sgfa" --bubble-prefix-in "$d/bubble/bubble" -o "$d/panphorte/panphorte" \
    --reference-path "$ref" --min-similarity "$pan_sim" --threads "$THREADS" "${GTFOPT[@]}" --quiet || return 1
  local pgfa="$d/panphorte/panphorte.normalized.sorted.gfa"

  # graph + bubble prefix that call/describe read, per topology
  local cgfa cpfx
  if [[ "$call_graph" == "bubble" ]]; then cgfa="$sgfa"; cpfx="$d/bubble/bubble";
  else cgfa="$pgfa"; cpfx="$d/panphorte/panphorte"; fi

  # 4) call
  "$BIN" call -i "$cgfa" --bubble-prefix-in "$cpfx" --reference-path "$ref" \
    -o "$d/call/call" --threads "$THREADS" "${GTFOPT[@]}" --quiet $call_extra || return 1

  # 5) inspect the main copy-number bubble on the CALL graph (node ids line up 1:1 with the VCF)
  local cbub; cbub="$(main_dup_bubble "$d/call/call.region.vcf")"
  [[ -n "$cbub" ]] && "$BIN" inspect -i "$cgfa" --bubble-prefix-in "$cpfx" --bubble-id "$cbub" \
      --cluster --cluster-similarity "$clu" -o "$d/inspect/inspect" --quiet || true

  # 6) describe (variant-restricted markers; BIMBAM dosage + per-bubble feature tables)
  "$BIN" describe -i "$cgfa" --bubble-prefix-in "$cpfx" --out-dir "$d/describe" \
    --variant-nodes "$d/call/call.variant_nodes.tsv" --no-wide-matrix --threads "$THREADS" --quiet || return 1

  echo "[$region] data done -> $d  (bubble bubble=$bbub ; call DUP bubble=$cbub)"
}
