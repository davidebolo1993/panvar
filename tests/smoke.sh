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

# Snarls live next to the GFA: strip .gz then .gfa, append .snarls.jsonl.
BASE="${GFA%.gz}"; BASE="${BASE%.gfa}"
SNARLS="${BASE}.snarls.jsonl"
if [[ ! -f "$SNARLS" ]]; then
  echo "error: missing snarls JSONL next to GFA: $SNARLS" >&2
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

# Optional external tools (full pipeline incl. re-snarl). CI without vg still runs
# the core pipeline (call on the original graph).
VG_BIN="${VG_BIN:-$(command -v vg || true)}"
BCFTOOLS_BIN="${BCFTOOLS_BIN:-$(command -v bcftools || true)}"
ODGI_BIN="${ODGI_BIN:-$(command -v odgi || true)}"

mkdir -p "$OUT_DIR"
BUBBLE_PREFIX="$OUT_DIR/bubble"
DESCRIBE_DIR="$OUT_DIR/describe"
PANPHORTE_PREFIX="$OUT_DIR/panphorte"
CALL_PREFIX="$OUT_DIR/call"

# 1) bubble
"$PANVAR_BIN" bubble -i "$GFA" -o "$BUBBLE_PREFIX" --snarls-in "$SNARLS" --quiet

FIRST_BUBBLE_ID="$(awk -F ',' 'NR==2 { print $1; exit }' "$BUBBLE_PREFIX.bubbles.csv")"

# 2) inspect (single bubble + all)
if [[ -n "$FIRST_BUBBLE_ID" ]]; then
  INSPECT_PREFIX="$OUT_DIR/inspect/bubble_${FIRST_BUBBLE_ID}"
  INSPECT_ALL_PREFIX="$OUT_DIR/inspect/all"
  "$PANVAR_BIN" inspect -i "$GFA" --bubble-prefix-in "$BUBBLE_PREFIX" \
    --bubble-id "$FIRST_BUBBLE_ID" -o "$INSPECT_PREFIX" --cluster --quiet >/dev/null
  "$PANVAR_BIN" inspect -i "$GFA" --bubble-prefix-in "$BUBBLE_PREFIX" \
    -o "$INSPECT_ALL_PREFIX" --quiet >/dev/null

  # 3) describe
  "$PANVAR_BIN" describe -i "$GFA" --bubble-prefix-in "$BUBBLE_PREFIX" \
    --bubble-id "$FIRST_BUBBLE_ID" --out-dir "$DESCRIBE_DIR" \
    --kmer-size 21 --max-wide-features 0 --quiet >/dev/null
fi

# 4) panphorte (normalize tandem-repeat bubbles)
"$PANVAR_BIN" panphorte -i "$GFA" --bubble-prefix-in "$BUBBLE_PREFIX" \
  -o "$PANPHORTE_PREFIX" --min-similarity 0.90 --quiet

# 5) call. With vg, run the full intended pipeline (re-snarl the normalized graph,
#    re-bubble, then call on the normalized graph). Without vg, call on the original.
if [[ -n "$VG_BIN" ]]; then
  NORM_GFA="$PANPHORTE_PREFIX.normalized.gfa"
  NORM_SNARLS="$PANPHORTE_PREFIX.normalized.snarls.jsonl"
  NORM_BUBBLE="$OUT_DIR/bubble_norm"

  # panphorte appends new REP nodes at the END of the node list, so numeric node id
  # no longer follows the reference. odgi sort repositions them topologically along
  # the reference, restoring the "numeric id == reference order" invariant that
  # call's node_track.tsv / plotting assume. Optional: skipped when odgi is absent.
  if [[ -n "$ODGI_BIN" ]]; then
    SORTED_GFA="$PANPHORTE_PREFIX.normalized.sorted.gfa"
    "$ODGI_BIN" build -g "$NORM_GFA" -o - | "$ODGI_BIN" sort -i - -o - -p Ygs -P \
      | "$ODGI_BIN" view -i - -g > "$SORTED_GFA"
    NORM_GFA="$SORTED_GFA"
  else
    echo "note: odgi not found; calling on the unsorted normalized graph (REP nodes appended at end)"
  fi

  "$VG_BIN" snarls -A integrated "$NORM_GFA" | "$VG_BIN" view -R -j - > "$NORM_SNARLS"
  "$PANVAR_BIN" bubble -i "$NORM_GFA" -o "$NORM_BUBBLE" --snarls-in "$NORM_SNARLS" --quiet
  "$PANVAR_BIN" call -i "$NORM_GFA" --bubble-prefix-in "$NORM_BUBBLE" \
    --reference-path "$REF_PATH" -o "$CALL_PREFIX" --classify-ins --cn-from-multiplicity --quiet
else
  echo "note: vg not found; calling on the original graph (skipping re-snarl stage)"
  "$PANVAR_BIN" call -i "$GFA" --bubble-prefix-in "$BUBBLE_PREFIX" \
    --reference-path "$REF_PATH" -o "$CALL_PREFIX" --cn-from-multiplicity --quiet
fi

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
test -s "$PANPHORTE_PREFIX.normalized.gfa"
test -s "$PANPHORTE_PREFIX.panphorte.report.tsv"
test -s "$CALL_PREFIX.region.vcf"
test -s "$CALL_PREFIX.variant_paths.tsv"
test -s "$CALL_PREFIX.node_track.tsv"

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
