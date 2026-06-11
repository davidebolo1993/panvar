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
CALL_PREFIX="$OUT_DIR/call"
DESCRIBE_DIR="$OUT_DIR/describe"

"$PANVAR_BIN" bubble \
  -i "$GFA" \
  -o "$BUBBLE_PREFIX" \
  --snarls-in "$SNARLS"

FIRST_BUBBLE_ID="$(awk -F ',' 'NR==2 { print $1; exit }' "$BUBBLE_PREFIX.bubbles.csv")"
if [[ -n "$FIRST_BUBBLE_ID" ]]; then
  INSPECT_PREFIX="$OUT_DIR/inspect/bubble_${FIRST_BUBBLE_ID}"
  INSPECT_ALL_PREFIX="$OUT_DIR/inspect/all"
  "$PANVAR_BIN" inspect \
    -i "$GFA" \
    --bubble-prefix-in "$BUBBLE_PREFIX" \
    --bubble-id "$FIRST_BUBBLE_ID" \
    -o "$INSPECT_PREFIX" >/dev/null
  "$PANVAR_BIN" inspect \
    -i "$GFA" \
    --bubble-prefix-in "$BUBBLE_PREFIX" \
    -o "$INSPECT_ALL_PREFIX" >/dev/null

  "$PANVAR_BIN" describe \
    -i "$GFA" \
    --bubble-prefix-in "$BUBBLE_PREFIX" \
    --bubble-id "$FIRST_BUBBLE_ID" \
    --out-dir "$DESCRIBE_DIR" \
    --kmer-size 21 \
    --max-wide-features 0 \
    --quiet >/dev/null
fi

"$PANVAR_BIN" call \
  -i "$GFA" \
  -o "$CALL_PREFIX" \
  --bubble-prefix-in "$BUBBLE_PREFIX" \
  --reference-path "$REF_PATH" \
  --quiet

test -s "$BUBBLE_PREFIX.bubbles.csv"
if [[ -n "${FIRST_BUBBLE_ID:-}" ]]; then
  test -s "$INSPECT_PREFIX.bubble_${FIRST_BUBBLE_ID}.paths.fa.gz"
  test -s "$INSPECT_PREFIX.bubble_${FIRST_BUBBLE_ID}.node_counts.tsv"
  gzip -t "$INSPECT_PREFIX.bubble_${FIRST_BUBBLE_ID}.paths.fa.gz"
  test -s "$INSPECT_ALL_PREFIX.bubble_${FIRST_BUBBLE_ID}.paths.fa.gz"
  test -s "$INSPECT_ALL_PREFIX.bubble_${FIRST_BUBBLE_ID}.node_counts.tsv"
  gzip -t "$INSPECT_ALL_PREFIX.bubble_${FIRST_BUBBLE_ID}.paths.fa.gz"
  test -s "$DESCRIBE_DIR/describe.index.tsv"
  DESCRIBE_BUBBLE_DIR="$DESCRIBE_DIR/bubble_${FIRST_BUBBLE_ID}"
  test -s "$DESCRIBE_BUBBLE_DIR/kmer_features.tsv.gz"
  test -s "$DESCRIBE_BUBBLE_DIR/kmer_counts.jsonl.gz"
  test -s "$DESCRIBE_BUBBLE_DIR/kmer_matrix.tsv.gz"
  gzip -t "$DESCRIBE_BUBBLE_DIR/kmer_features.tsv.gz"
  gzip -t "$DESCRIBE_BUBBLE_DIR/kmer_counts.jsonl.gz"
  gzip -t "$DESCRIBE_BUBBLE_DIR/kmer_matrix.tsv.gz"
fi
test -s "$CALL_PREFIX.region.vcf"
