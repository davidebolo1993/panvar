#!/usr/bin/env bash
# Gate: the genotype tests and developer probes leave no generated files in the source tree. Generated graphs,
# tables, reads and logs belong under the build directory or the gitignored tests/results. An untracked
# genotype-named file with a generated-output extension under the tracked source directories fails this test;
# new source files that are simply not committed yet do not.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
stray="$(git -C "$REPO" status --porcelain --untracked-files=all -- tests src include docs scripts \
  | awk '$1=="??"{print $2}' | grep -i 'genotype' | grep -E '\.(gfa|tsv|fa|fasta|fq|fastq|log|time)$' || true)"
if [ -n "$stray" ]; then
  echo "generated genotype files in the source tree:" >&2
  echo "$stray" >&2
  exit 1
fi
echo "genotype source tree: no generated files"
