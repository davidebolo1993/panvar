#!/usr/bin/env bash
# legacy_ab.sh - prove that a change leaves every LEGACY output byte-identical.
#
#   tests/legacy_ab.sh <binary-A> <binary-B> <workdir>
#
# WHY THIS EXISTS IN THE REPO. The previous harness lived in a scratch directory and compared a
# freshly generated manifest against a COPIED one whose digests could not be reproduced from the
# files they named -- by sha256, sha3-256, blake2s, sorted content, CR-stripped content or
# path-normalised content. A manifest nobody can recompute is not a baseline; it is a number that
# agrees with itself. So this harness:
#
#   1. builds BOTH sides from binaries the caller pins, into fresh directories;
#   2. generates its own phenotype inputs DETERMINISTICALLY (they are inputs, not outputs -- the
#      old harness produced them out of band, which is why a replay of its command log could not
#      run the association steps at all);
#   3. hashes with one documented command: `shasum -a 256 <file>`;
#   4. RECOMPUTES AND VERIFIES its own manifest against the tree it just wrote, so a manifest that
#      does not describe its own files is caught here rather than by whoever trusts it next;
#   5. normalises ONLY the fields that are legitimately volatile -- the absolute output root, which
#      `describe` embeds in describe.index.tsv and describe.params.json -- and nothing else.
set -u
A=${1:?binary A (the pinned baseline)}
B=${2:?binary B (the candidate)}
W=${3:?workdir}
A=$(cd "$(dirname "$A")" && pwd)/$(basename "$A")
B=$(cd "$(dirname "$B")" && pwd)/$(basename "$B")
REPO=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$W"; W=$(cd "$W" && pwd)

REFS_c4='grch38#1#chr6:31891045-32123783'
REFS_f7='GRCh38#0#chr13:113105788-113121779'
REFS_cyp2d6='grch38#1#chr22:42031864-42245566'

# Deterministic synthetic phenotypes, from the sample name alone. Values are arbitrary; what
# matters is that both arms get the SAME ones, which a name-derived hash guarantees without any
# stored file to drift.
gen_pheno () {  # <gfa.gz> <outdir>
  local gfa=$1 out=$2
  gzcat "$gfa" 2>/dev/null | awk -F'\t' '$1=="P"{print $2}' | sort -u > "$out/.samples"
  for kind in kmers graph variant; do
    awk -v k="$kind" 'BEGIN{OFS="\t"; print "sample","phenotype","bin","cov1","cov2"}
      { h=0; s=$0 k; for(i=1;i<=length(s);i++) h=(h*31+index("ACGTacgt0123456789#-:.",substr(s,i,1))+i)%100003;
        printf "%s\t%.4f\t%d\t%d\t%d\n", $0, 1.0+(h%1000)/1000.0, h%2, h%5, (h%7)-3 }' \
      "$out/.samples" > "$out/pheno.$kind.tsv"
    awk 'BEGIN{OFS="\t"} NR==1{print "sample","phenotype","cov1","cov2"; next}
         {print $1,$3,$4,$5}' "$out/pheno.$kind.tsv" > "$out/pheno.$kind.tsv.bin"
  done
  rm -f "$out/.samples"
}

run_arm () {  # <binary> <root>
  local P=$1 R=$2
  rm -rf "$R"; mkdir -p "$R"
  for L in c4 f7 cyp2d6; do
    local ref; eval "ref=\$REFS_$L"
    local G=$REPO/tests/real_data/$L.gfa.gz
    local D=$R/$L; mkdir -p "$D/inspect"
    "$P" rebuild   -i "$G" -o "$D/rebuilt.gfa" --quiet || return 1
    "$P" bubble    -i "$D/rebuilt.gfa" -r "$ref" -o "$D/bub" --quiet || return 1
    "$P" inspect   -i "$D/bub.sorted.gfa" --bubble-prefix-in "$D/bub" --bubble-id 1 \
                   -o "$D/inspect/one" --cluster --quiet || return 1
    "$P" inspect   -i "$D/bub.sorted.gfa" --bubble-prefix-in "$D/bub" \
                   -o "$D/inspect/all" --quiet || return 1
    "$P" panphorte -i "$D/bub.sorted.gfa" --bubble-prefix-in "$D/bub" -r "$ref" \
                   -o "$D/pan" --min-similarity 0.90 --quiet || return 1
    "$P" refine    -i "$D/pan.normalized.sorted.gfa" --bubble-prefix-in "$D/pan" -r "$ref" \
                   -o "$D/ref" --quiet || return 1
    "$P" call      -i "$D/ref.normalized.sorted.gfa" --bubble-prefix-in "$D/ref" -r "$ref" \
                   -o "$D/refcall" --cn --quiet || return 1
    "$P" call      -i "$D/pan.normalized.sorted.gfa" --bubble-prefix-in "$D/pan" -r "$ref" \
                   -o "$D/call" --cn --gtf "$REPO/tests/real_data/Homo_sapiens.GRCh38.116.gtf.gz" \
                   --quiet || return 1
    "$P" benchmark -i "$D/pan.normalized.sorted.gfa" --bubble-prefix-in "$D/pan" -r "$ref" \
                   --variant-nodes "$D/call.variant_nodes.tsv" -o "$D/bench" --quiet || return 1
    "$P" describe  -i "$D/pan.normalized.sorted.gfa" --bubble-prefix-in "$D/pan" \
                   --out-dir "$D/desc" --kmer-size 31 --no-wide-matrix \
                   --variant-nodes "$D/call.variant_nodes.tsv" \
                   --variant-vcf "$D/call.region.vcf" --quiet || return 1
    gen_pheno "$G" "$D"
    for kind in kmers graph variant; do
      for suf in "" ".bin"; do
        local tag=quant; [ -n "$suf" ] && tag=binary
        "$P" associate \
          --genotypes    "$D/desc/haplotype/$kind/bimbam_$kind.bimbam.gz" \
          --samples      "$D/desc/haplotype/$kind/samples.txt.gz" \
          --feature-annot "$D/desc/haplotype/$kind/feature_annot.$kind.tsv.gz" \
          --phenotype    "$D/pheno.$kind.tsv$suf" \
          -o "$D/assoc.$kind.$tag" --quiet || return 1
      done
    done
  done
}

# ONE documented digest command, over content with the output root normalised away.
manifest () {  # <root> <outfile>
  local R=$1 M=$2
  ( cd "$R" && find . -type f | LC_ALL=C sort | while read -r f; do
      printf '%s\t%s\n' "${f#./}" \
        "$(sed "s#$R#<ROOT>#g" "$f" | shasum -a 256 | cut -d' ' -f1)"
    done ) > "$M"
}

verify_manifest () {  # <root> <manifest>  -- the harness checking its own output
  local R=$1 M=$2 tmp; tmp=$(mktemp)
  manifest "$R" "$tmp"
  if cmp -s "$M" "$tmp"; then rm -f "$tmp"; return 0; fi
  echo "MANIFEST DOES NOT DESCRIBE ITS OWN TREE ($R)"; diff "$M" "$tmp" | head -5; rm -f "$tmp"; return 1
}

echo "A: $A"; echo "B: $B"
run_arm "$A" "$W/A" || { echo "arm A FAILED"; exit 1; }
run_arm "$B" "$W/B" || { echo "arm B FAILED"; exit 1; }
manifest "$W/A" "$W/A.manifest.tsv"
manifest "$W/B" "$W/B.manifest.tsv"
verify_manifest "$W/A" "$W/A.manifest.tsv" || exit 1
verify_manifest "$W/B" "$W/B.manifest.tsv" || exit 1
echo "manifests verify against their own trees ($(wc -l < "$W/A.manifest.tsv" | tr -d ' ') files each)"
if diff "$W/A.manifest.tsv" "$W/B.manifest.tsv" > "$W/diff.txt"; then
  echo "LEGACY A/B IDENTICAL ($(wc -l < "$W/A.manifest.tsv" | tr -d ' ') files)"
else
  echo "LEGACY A/B DIFFERS ($(grep -c '^[<>]' "$W/diff.txt") lines):"; head -20 "$W/diff.txt"; exit 1
fi
