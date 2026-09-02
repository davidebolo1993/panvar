#!/usr/bin/env bash
# ab_manifest.sh - a normalised, byte-identity manifest of the whole pipeline's output.
#
# Purpose: prove that a refactor changed NOTHING. Capture a manifest before the change, capture one
# after, diff them. An empty diff is the proof; any line in the diff is either a bug the refactor
# introduced or a difference we argued for in advance.
#
# It is deliberately separate from real_smoke.sh. That script asserts that outputs EXIST and satisfy
# a handful of contracts; this one asserts that their CONTENT has not moved. The two answer different
# questions and a refactor needs the second one.
#
#   ab_manifest.sh <panvar-binary> <out-dir> [locus ...]
#
# Loci default to the fast set (c4 f7 cyp2d6); pass names, or ALL, to widen. A locus name is the
# basename of tests/real_data/<name>.gfa.gz.
#
# Env:
#   PANVAR_AB_THREADS  passed to every module that takes --threads (unset = module default).
#                      Set to 1 and diff against the default run to test thread invariance.
#   PANVAR_AB_GTF      1 (default) to include --gtf, exercising gtf.cpp + gene_cn_kmer.cpp.
#   PANVAR_AB_ASSOC    1 (default) to run associate on the committed LPA GWAS fixtures when the
#                      `lpa` locus is in the set.
#
# WHY THE NORMALISATION BELOW IS NEEDED (each rule is a real source of run-to-run difference, not a
# precaution): gzip members carry an mtime in their header, so two identical payloads hash
# differently; describe.index.tsv and describe.params.json embed input/output PATHS, which vary with
# where the run was rooted; describe.index.tsv currently embeds the staging directory name, which
# carries the PID; and describe.params.json records the thread count it was INVOKED with, which is a
# record of the call rather than of the result -- normalising it is what lets a 1-thread run be diffed
# against a default run to test thread invariance directly.
# Normalising these is what makes "identical output" a testable claim at all.
set -uo pipefail

BIN="${1:?usage: ab_manifest.sh <panvar-binary> <out-dir> [locus ...]}"
OUT="${2:?usage: ab_manifest.sh <panvar-binary> <out-dir> [locus ...]}"
shift 2

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="$HERE/real_data"
GTF_FILE="$DATA/Homo_sapiens.GRCh38.116.gtf.gz"
GWAS="$HERE/gwas/lpa"

FAST_LOCI=(c4 f7 cyp2d6)
# myom2 is deliberately absent: panphorte refuses it at HEAD ("2 rewritten path(s) can still be
# walked along their ORIGINAL route", bubble 7), which is a known unimplemented case -- the error
# says so -- and `rebuild` does not clear it, because rebuild REJECTS its own output there (313
# haplotypes below the identity bound) and passes the original through. A harness that always
# reports one expected failure teaches its readers to ignore failures, so the locus is named here
# rather than left to fail on every run. Pass it explicitly to reproduce.
ALL_LOCI=(c4 f7 cyp2d6 gstm1 acot lpa ankrd36c)

if [ "$#" -eq 0 ]; then
  LOCI=("${FAST_LOCI[@]}")
elif [ "$1" = "ALL" ]; then
  LOCI=("${ALL_LOCI[@]}")
else
  LOCI=("$@")
fi

[ -x "$BIN" ] || { echo "error: not executable: $BIN" >&2; exit 2; }
BIN="$(cd "$(dirname "$BIN")" && pwd)/$(basename "$BIN")"

rm -rf "$OUT"; mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"
RUN="$OUT/run"; mkdir -p "$RUN"
MANIFEST="$OUT/manifest.tsv"
LOG="$OUT/commands.log"
: > "$LOG"

THREAD_OPT=()
[ -n "${PANVAR_AB_THREADS:-}" ] && THREAD_OPT=(--threads "$PANVAR_AB_THREADS")

WANT_GTF="${PANVAR_AB_GTF:-1}"
WANT_ASSOC="${PANVAR_AB_ASSOC:-1}"

SHA=""
for c in sha256sum "shasum -a 256"; do
  if command -v ${c%% *} >/dev/null 2>&1; then SHA="$c"; break; fi
done
[ -n "$SHA" ] || { echo "error: no sha256sum or shasum on PATH" >&2; exit 2; }

fails=0
note() { printf '  %s\n' "$1"; }
run() {  # run <label> <cmd...>  -- record, execute, keep going on failure so the manifest is complete
  local label="$1"; shift
  printf '### %s\n%s\n' "$label" "$*" >> "$LOG"
  if ! "$@" >> "$LOG" 2>&1; then
    printf '  FAIL %s (see %s)\n' "$label" "$LOG"
    fails=$((fails + 1))
    return 1
  fi
  return 0
}

# The reference path of a locus: prefer a GRCh38-looking name so the choice does not depend on file
# order, else the first P/W line. Same rule real_regressions.sh uses.
ref_of() {
  gunzip -c "$1" 2>/dev/null | awk -F'\t' '($1=="P"||$1=="W"){n=$2;
    if(n~/[Gg][Rr][Cc]h38/){print n; exit} if(!f)f=n} END{if(f)print f}' | head -1
}

# ---------------------------------------------------------------- the pipeline, per locus
for L in "${LOCI[@]}"; do
  G="$DATA/$L.gfa.gz"
  if [ ! -f "$G" ]; then note "skip $L (no $G)"; continue; fi
  R="$(ref_of "$G")"
  if [ -z "$R" ]; then note "skip $L (no P/W path)"; continue; fi
  D="$RUN/$L"; mkdir -p "$D"
  printf '== %s (ref %s)\n' "$L" "$R"

  # rebuild first: it is the documented pre-step to bubble, and it is a pass-through on a healthy
  # graph, so including it costs ~0.5s per locus and is the only coverage this module gets here.
  # Its output feeds bubble, which is how the pipeline is meant to run.
  run "$L rebuild" "$BIN" rebuild -i "$G" -o "$D/rebuilt.gfa" --quiet
  IN="$G"
  [ -s "$D/rebuilt.gfa" ] && IN="$D/rebuilt.gfa"

  run "$L bubble" "$BIN" bubble -i "$IN" -r "$R" -o "$D/bub" ${THREAD_OPT[@]+"${THREAD_OPT[@]}"} --quiet || continue
  SORTED="$D/bub.sorted.gfa"

  # inspect: one bubble (with clustering) and the all-bubbles sweep.
  FIRST="$(awk -F',' 'NR==2{print $1; exit}' "$D/bub.bubbles.csv" 2>/dev/null)"
  if [ -n "$FIRST" ]; then
    mkdir -p "$D/inspect"
    run "$L inspect one" "$BIN" inspect -i "$SORTED" --bubble-prefix-in "$D/bub" \
      --bubble-id "$FIRST" -o "$D/inspect/one" --cluster --quiet
    run "$L inspect all" "$BIN" inspect -i "$SORTED" --bubble-prefix-in "$D/bub" \
      -o "$D/inspect/all" --quiet
  fi

  run "$L panphorte" "$BIN" panphorte -i "$SORTED" --bubble-prefix-in "$D/bub" \
    -r "$R" -o "$D/pan" --min-similarity 0.90 ${THREAD_OPT[@]+"${THREAD_OPT[@]}"} --quiet || continue
  NORM="$D/pan.normalized.sorted.gfa"

  run "$L refine" "$BIN" refine -i "$NORM" --bubble-prefix-in "$D/pan" \
    -r "$R" -o "$D/ref" --quiet
  if [ -s "$D/ref.normalized.sorted.gfa" ]; then
    run "$L call(refined)" "$BIN" call -i "$D/ref.normalized.sorted.gfa" \
      --bubble-prefix-in "$D/ref" -r "$R" -o "$D/refcall" --cn ${THREAD_OPT[@]+"${THREAD_OPT[@]}"} --quiet
  fi

  # call on the panphorte graph, with the GTF when asked (this is the only path that reaches
  # gtf.cpp and gene_cn_kmer.cpp).
  GTF_OPT=()
  [ "$WANT_GTF" = "1" ] && [ -f "$GTF_FILE" ] && GTF_OPT=(--gtf "$GTF_FILE")
  run "$L call" "$BIN" call -i "$NORM" --bubble-prefix-in "$D/pan" -r "$R" -o "$D/call" \
    --cn ${GTF_OPT[@]+"${GTF_OPT[@]}"} ${THREAD_OPT[@]+"${THREAD_OPT[@]}"} --quiet || continue

  run "$L benchmark" "$BIN" benchmark -i "$NORM" --bubble-prefix-in "$D/pan" -r "$R" \
    --variant-nodes "$D/call.variant_nodes.tsv" -o "$D/bench" ${THREAD_OPT[@]+"${THREAD_OPT[@]}"} --quiet

  # describe: all bubbles, all three substrates, BIMBAM on. The wide matrix is off because its size
  # is quadratic in features and it adds no coverage the sparse tables do not already give.
  run "$L describe" "$BIN" describe -i "$NORM" --bubble-prefix-in "$D/pan" --out-dir "$D/desc" \
    --kmer-size 31 --no-wide-matrix --variant-nodes "$D/call.variant_nodes.tsv" \
    --variant-vcf "$D/call.region.vcf" ${THREAD_OPT[@]+"${THREAD_OPT[@]}"} --quiet

  # associate, on this locus's own haplotype-level BIMBAM.
  #
  # The phenotype is synthesised from the haplotype NAMES with a fixed FNV-1a hash, not drawn at
  # random: associate has to be in the FAST set or a change to it is only ever checked by the slow
  # run, and a random cohort would make the manifest unstable and useless as a baseline. Both a
  # quantitative and a binary phenotype are emitted so the linear path and the score/SPA path are
  # both exercised, each with two covariates so the covariate-adjusted code is reached too.
  if [ "$WANT_ASSOC" = "1" ]; then
    for SUB in kmers graph variant; do
      SD="$D/desc/haplotype/$SUB"
      GENO="$SD/bimbam_$SUB.bimbam.gz"
      [ -f "$GENO" ] && [ -f "$SD/samples.txt.gz" ] || continue
      PHENO="$D/pheno.$SUB.tsv"
      gunzip -c "$SD/samples.txt.gz" | awk '
        BEGIN { print "sample\tphenotype\tbin\tcov1\tcov2" }
        {
          # A fixed positional hash over the name. Plain arithmetic rather than xor(), which is a
          # gawk extension the BSD awk on macOS does not have.
          v = 2166136261
          n = length($0)
          for (i = 1; i <= n; i++) { c = index(ORDER, substr($0, i, 1)); v = (v + c * i * 131) % 1000003 }
          q = (v % 10000) / 1000.0
          b = (v % 2)
          c1 = (v % 7) - 3
          c2 = ((int(v / 7)) % 5) - 2
          printf "%s\t%.4f\t%d\t%d\t%d\n", $0, q, b, c1, c2
        }' ORDER=" ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789#-_.:" > "$PHENO"
      # Quantitative (linear + Student-t tail) and binary (Rao score + saddlepoint) in turn.
      run "$L associate $SUB quant" "$BIN" associate --genotypes "$GENO" \
        --samples "$SD/samples.txt.gz" --feature-annot "$SD/feature_annot.$SUB.tsv.gz" \
        --phenotype "$PHENO" -o "$D/assoc.$SUB.quant" --quiet
      awk -F'\t' 'NR==1{print "sample\tphenotype\tcov1\tcov2"; next}
                  {print $1"\t"$3"\t"$4"\t"$5}' "$PHENO" > "$PHENO.bin"
      run "$L associate $SUB binary" "$BIN" associate --genotypes "$GENO" \
        --samples "$SD/samples.txt.gz" --feature-annot "$SD/feature_annot.$SUB.tsv.gz" \
        --phenotype "$PHENO.bin" -o "$D/assoc.$SUB.binary" --quiet
    done
  fi

  # associate, on real describe output with the committed (deterministic) LPA cohort. Adds the
  # SAMPLE-level (diploid) substrate and a 6,000-individual cohort, which the per-locus block above
  # cannot reach; only runs when the lpa locus is in the set.
  if [ "$L" = "lpa" ] && [ "$WANT_ASSOC" = "1" ] && [ -f "$GWAS/samples.tsv" ]; then
    run "$L describe(samples)" "$BIN" describe -i "$NORM" --bubble-prefix-in "$D/pan" \
      --out-dir "$D/desc_s" --kmer-size 31 --no-wide-matrix \
      --variant-nodes "$D/call.variant_nodes.tsv" --variant-vcf "$D/call.region.vcf" \
      --samples "$GWAS/samples.tsv" ${THREAD_OPT[@]+"${THREAD_OPT[@]}"} --quiet
    for SUB in variant graph; do
      SD="$D/desc_s/sample/$SUB"
      [ -f "$SD/bimbam.gz" ] || continue
      for PH in pheno.quant pheno.binary; do
        [ -f "$GWAS/$PH.tsv" ] || continue
        run "$L associate $SUB $PH" "$BIN" associate \
          --genotypes "$SD/bimbam.gz" --samples "$SD/samples.txt.gz" \
          --feature-annot "$SD/feature_annot.tsv.gz" \
          --phenotype "$GWAS/$PH.tsv" -o "$D/assoc.$SUB.$PH" --quiet
      done
    done
  fi
done

# ---------------------------------------------------------------- normalise + hash
#
# Every file produced, hashed after normalisation. Directories are walked in sorted order so the
# manifest is stable regardless of filesystem enumeration order.
normalise_and_hash() {  # normalise_and_hash <abs-path> -> sha256 hex on stdout
  local f="$1"
  local decoded
  case "$f" in
    *.gz)
      # gzip headers carry an mtime, so identical payloads hash differently. Decompress first.
      # A non-gzip file that merely ends in .gz falls back to raw bytes.
      decoded="$(gunzip -c "$f" 2>/dev/null)" || decoded="$(cat "$f")"
      ;;
    *) decoded="$(cat "$f" 2>/dev/null)" ;;
  esac
  printf '%s' "$decoded" | sed \
    -e "s#$RUN#@RUN@#g" \
    -e "s#$OUT#@OUT@#g" \
    -e "s#$HERE#@TESTS@#g" \
    -e 's#\.describe-staging\.[0-9][0-9]*#.describe-staging.@PID@#g' \
    -e 's#-tmp\.[0-9][0-9]*\.[0-9][0-9]*#-tmp.@STAGE@#g' \
    -e 's#"threads": [0-9][0-9]*#"threads": @THREADS@#g' \
    | $SHA | awk '{print $1}'
}

: > "$MANIFEST"
find "$RUN" -type f -print0 \
  | LC_ALL=C sort -z \
  | while IFS= read -r -d '' f; do
      printf '%s\t%s\n' "${f#$RUN/}" "$(normalise_and_hash "$f")" >> "$MANIFEST"
    done

N="$(wc -l < "$MANIFEST" | tr -d ' ')"
printf '\nmanifest: %s (%s files)\n' "$MANIFEST" "$N"
if [ "$fails" -gt 0 ]; then
  printf 'ab_manifest: %d command(s) FAILED - see %s\n' "$fails" "$LOG"
  exit 1
fi
[ "$N" -gt 0 ] || { echo "ab_manifest: manifest is empty" >&2; exit 1; }
echo "ab_manifest: OK"
