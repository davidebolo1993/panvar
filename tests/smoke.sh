#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <panvar_bin> <graph.gfa> <out_dir>" >&2
  exit 2
fi

PANVAR_BIN="$1"
GFA="$2"
OUT_DIR="$3"

if [[ ! -x "$PANVAR_BIN" ]]; then
  echo "error: panvar binary not executable: $PANVAR_BIN" >&2
  exit 1
fi
if [[ ! -f "$GFA" ]]; then
  echo "error: missing GFA: $GFA" >&2
  exit 1
fi

SNARLS="${GFA%.gfa}.snarls.jsonl"
if [[ ! -f "$SNARLS" ]]; then
  echo "error: missing snarls JSONL next to GFA: $SNARLS" >&2
  exit 1
fi

REF_PATH="$(awk -F '\t' '$1=="P" || $1=="W" { print $2; exit }' "$GFA")"
if [[ -z "$REF_PATH" ]]; then
  echo "error: could not infer a reference path from $GFA" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

BUBBLE_PREFIX="$OUT_DIR/bubble"
ALLELE_PREFIX="$OUT_DIR/allele"
CALL_PREFIX="$OUT_DIR/call"

"$PANVAR_BIN" bubble \
  -i "$GFA" \
  -o "$BUBBLE_PREFIX" \
  --snarls-in "$SNARLS"

"$PANVAR_BIN" allele \
  -i "$GFA" \
  -o "$ALLELE_PREFIX" \
  --bubble-prefix-in "$BUBBLE_PREFIX" \
  --quiet

"$PANVAR_BIN" call \
  -i "$GFA" \
  -o "$CALL_PREFIX" \
  --bubble-prefix-in "$BUBBLE_PREFIX" \
  --allele-prefix-in "$ALLELE_PREFIX" \
  --reference-path "$REF_PATH" \
  --quiet

test -s "$BUBBLE_PREFIX.bubbles.csv"
test -s "$ALLELE_PREFIX.allele_clusters.csv"
test -s "$ALLELE_PREFIX.allele_assignments.csv"
test -s "$CALL_PREFIX.region.vcf"
