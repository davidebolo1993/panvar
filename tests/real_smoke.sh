#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <panvar_bin> <graph.gfa[.gz]> <out_dir>" >&2
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

# Read the GFA whether plain or gzipped. Print the first path name without an early
# awk `exit` (which would SIGPIPE the upstream reader under `set -o pipefail`).
gfa_cat() { case "$GFA" in *.gz) gzip -dc "$GFA";; *) cat "$GFA";; esac; }
REF_PATH="$(gfa_cat | awk -F '\t' '($1=="P" || $1=="W") && !f { print $2; f=1 }')"
if [[ -z "$REF_PATH" ]]; then
  echo "error: could not infer a reference path from $GFA" >&2
  exit 1
fi

# Optional: only used to validate the final VCF is bgzip/tabix-indexable.
BCFTOOLS_BIN="${BCFTOOLS_BIN:-$(command -v bcftools || true)}"

mkdir -p "$OUT_DIR"
BUBBLE_PREFIX="$OUT_DIR/bubble"
DESCRIBE_DIR="$OUT_DIR/describe"
PANPHORTE_PREFIX="$OUT_DIR/panphorte"
CALL_PREFIX="$OUT_DIR/call"

# 1) bubble: internal sort+flip + cactus snarls (no vg/odgi). Writes <prefix>.sorted.gfa.
"$PANVAR_BIN" bubble -i "$GFA" --reference-path "$REF_PATH" -o "$BUBBLE_PREFIX" --quiet
SORTED_GFA="$BUBBLE_PREFIX.sorted.gfa"
test -s "$SORTED_GFA"

FIRST_BUBBLE_ID="$(awk -F ',' 'NR==2 { print $1; exit }' "$BUBBLE_PREFIX.bubbles.csv")"

# 2) inspect (single bubble + all) on the sorted graph
if [[ -n "$FIRST_BUBBLE_ID" ]]; then
  INSPECT_PREFIX="$OUT_DIR/inspect/bubble_${FIRST_BUBBLE_ID}"
  INSPECT_ALL_PREFIX="$OUT_DIR/inspect/all"
  "$PANVAR_BIN" inspect -i "$SORTED_GFA" --bubble-prefix-in "$BUBBLE_PREFIX" \
    --bubble-id "$FIRST_BUBBLE_ID" -o "$INSPECT_PREFIX" --cluster --quiet >/dev/null
  "$PANVAR_BIN" inspect -i "$SORTED_GFA" --bubble-prefix-in "$BUBBLE_PREFIX" \
    -o "$INSPECT_ALL_PREFIX" --quiet >/dev/null

  # 3) describe
  "$PANVAR_BIN" describe -i "$SORTED_GFA" --bubble-prefix-in "$BUBBLE_PREFIX" \
    --bubble-id "$FIRST_BUBBLE_ID" --out-dir "$DESCRIBE_DIR" \
    --kmer-size 21 --max-wide-features 0 --quiet >/dev/null
fi

# 4) panphorte: normalize tandems AND internally re-sort + re-snarl (no vg/odgi).
#    Writes <prefix>.normalized.sorted.gfa + <prefix>.bubbles.csv, ready for call.
"$PANVAR_BIN" panphorte -i "$SORTED_GFA" --bubble-prefix-in "$BUBBLE_PREFIX" \
  --reference-path "$REF_PATH" -o "$PANPHORTE_PREFIX" --min-similarity 0.90 --quiet
NORM_GFA="$PANPHORTE_PREFIX.normalized.sorted.gfa"
test -s "$NORM_GFA"

# 5) call on the panphorte sorted output
"$PANVAR_BIN" call -i "$NORM_GFA" --bubble-prefix-in "$PANPHORTE_PREFIX" \
  --reference-path "$REF_PATH" -o "$CALL_PREFIX" --cn --quiet

# ---- assertions ----
test -s "$BUBBLE_PREFIX.bubbles.csv"
if [[ -n "${FIRST_BUBBLE_ID:-}" ]]; then
  test -s "$INSPECT_PREFIX.bubble_${FIRST_BUBBLE_ID}.paths.fa.gz"
  test -s "$INSPECT_PREFIX.bubble_${FIRST_BUBBLE_ID}.node_counts.tsv"
  test -s "$INSPECT_PREFIX.bubble_${FIRST_BUBBLE_ID}.edge_counts.tsv"
  test -s "$INSPECT_PREFIX.bubble_${FIRST_BUBBLE_ID}.node_lengths.tsv"
  test -s "$INSPECT_PREFIX.bubble_${FIRST_BUBBLE_ID}.clusters.tsv"
  gzip -t "$INSPECT_PREFIX.bubble_${FIRST_BUBBLE_ID}.paths.fa.gz"
  test -s "$INSPECT_ALL_PREFIX.bubble_${FIRST_BUBBLE_ID}.node_counts.tsv"
  test -s "$INSPECT_ALL_PREFIX.bubble_${FIRST_BUBBLE_ID}.node_lengths.tsv"
  test -s "$DESCRIBE_DIR/describe.index.tsv"
  DESCRIBE_BUBBLE_DIR="$DESCRIBE_DIR/bubble_${FIRST_BUBBLE_ID}"
  test -s "$DESCRIBE_BUBBLE_DIR/kmer_features.tsv.gz"
  gzip -t "$DESCRIBE_BUBBLE_DIR/kmer_features.tsv.gz"
fi
test -s "$PANPHORTE_PREFIX.normalized.sorted.gfa"   # reference given -> sorted graph (no plain .normalized.gfa)
test -s "$PANPHORTE_PREFIX.panphorte.report.tsv"
test -s "$CALL_PREFIX.region.vcf"
test -s "$CALL_PREFIX.variant_nodes.tsv"   # the describe --variant-nodes and benchmark handoff

# The region VCF must be coordinate-sorted (POS non-decreasing) and have unique IDs.
awk -F '\t' '!/^#/ { if ($2 < prev) { print "VCF not sorted at " $2 " < " prev > "/dev/stderr"; exit 1 } prev=$2 }' \
  "$CALL_PREFIX.region.vcf"
if awk -F '\t' '!/^#/ { print $3 }' "$CALL_PREFIX.region.vcf" | sort | uniq -d | grep -q .; then
  echo "error: duplicate VCF IDs in $CALL_PREFIX.region.vcf" >&2
  exit 1
fi

# Optional: validate the VCF parses and is bgzip/tabix-indexable.
if [[ -n "$BCFTOOLS_BIN" ]]; then
  "$BCFTOOLS_BIN" view "$CALL_PREFIX.region.vcf" >/dev/null
fi
if command -v bgzip >/dev/null && command -v tabix >/dev/null; then
  bgzip -f -c "$CALL_PREFIX.region.vcf" > "$CALL_PREFIX.region.vcf.gz"
  tabix -f -p vcf "$CALL_PREFIX.region.vcf.gz"
fi

echo "smoke: OK"
