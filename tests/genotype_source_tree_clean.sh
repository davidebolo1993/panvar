#!/usr/bin/env bash
# Gate: the genotype tests and developer probes leave no generated files in the source tree. Generated graphs,
# reads, logs and run tables belong under the build directory or the gitignored tests/results. An untracked
# genotype-named graph, FASTA/FASTQ, log or timing file under the tracked source directories fails this test, as
# does a run-output table (mosaic_*, genotypes, block_genotypes, resources, eval, trace). Hand-written source such as
# a settings or manifest .tsv that is simply not committed yet does not.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
stray="$(git -C "$REPO" status --porcelain --untracked-files=all -- tests src include docs scripts \
  | awk '$1=="??"{print $2}' | grep -i 'genotype' \
  | grep -E '\.(gfa|fa|fasta|fq|fastq|log|time)$|(\.mosaic_[a-z_]+|\.block_genotypes|\.genotypes|\.resources|\.eval\.[a-z_.]+|\.trace(\.summary)?)\.tsv$' || true)"
if [ -n "$stray" ]; then
  echo "generated genotype files in the source tree:" >&2
  echo "$stray" >&2
  exit 1
fi
echo "genotype source tree: no generated files"
