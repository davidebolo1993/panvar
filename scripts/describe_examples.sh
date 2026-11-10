#!/usr/bin/env bash
# Describe the graphs call_examples.sh called: k-mer, node/edge and variant features of every site with a record,
# per haplotype and per sample (samples from the PanSN path names).
#   scripts/describe_examples.sh [panvar binary] [output dir] [refine_examples.sh output dir] [call_examples.sh output dir]
set -euo pipefail
here=$PWD
cd "$(dirname "$0")/.."
# Relative paths given as arguments are read from the calling directory.
arg() { if [ -z "$1" ]; then echo "$2"; elif [ "${1#/}" != "$1" ]; then echo "$1"; else echo "$here/$1"; fi; }
# Run a step with stderr to a log; on failure show the log's end and stop.
run() { local log=$1; shift; "$@" 2> "$log" || { echo "failed: $* (log: $log)"; tail -n 3 "$log"; exit 1; }; }
command -v python3 > /dev/null || { echo "python3 not found on PATH"; exit 1; }

PANVAR=$(arg "${1:-}" build/panvar)
OUT=$(arg "${2:-}" tests/examples/describe)
REFINED=$(arg "${3:-}" tests/examples/refine)
CALLED=$(arg "${4:-}" tests/examples/call)

# <graph.gfa> <table>: samples with both haplotypes, from PanSN or sample_hap_contig P lines and from W lines;
# a haplotype in several contigs is one comma-separated cell
samples() {
  python3 - "$1" "$2" <<'EOF'
import collections, sys
haps = collections.defaultdict(lambda: collections.defaultdict(list))
numeric = lambda t: 0 < len(t) < 19 and all(c in '0123456789' for c in t)
with open(sys.argv[1]) as fh:
    for line in fh:
        if line.startswith('P\t'):
            name = line.split('\t')[1]
            parts = name.split('#') if '#' in name else name.split('_', 2)
            if len(parts) >= 3:
                haps[parts[0]][parts[1]].append(name)
        elif line.startswith('W\t'):
            sample, hap, seqid, start, end = line.rstrip('\n').split('\t')[1:6]
            name = '%s#%s#%s' % (sample, hap, seqid)
            if numeric(start):
                name += ':%d-%s' % (int(start) + 1, end if numeric(end) else '*')
            haps[sample][hap].append(name)
with open(sys.argv[2], 'w') as out:
    out.write('sample\thap1\thap2\n')
    for s in sorted(haps):
        if '1' in haps[s] and '2' in haps[s]:
            out.write('%s\t%s\t%s\n' % (s, ','.join(haps[s]['1']), ','.join(haps[s]['2'])))
EOF
}

describe() {  # <name>
  local name=$1 dir=$OUT/$1
  local gfa=$REFINED/$name/$name.normalized.sorted.gfa
  [ -f "$CALLED/$name/$name.region.vcf" ] ||
    { echo "no $CALLED/$name/$name.region.vcf: run scripts/call_examples.sh first"; exit 1; }
  mkdir -p "$dir"
  samples "$gfa" "$dir/$name.samples.tsv"
  run "$dir/$name.describe.log" "$PANVAR" describe -i "$gfa" -b "$REFINED/$name/$name" -r GRCh38 -o "$dir" \
    --variant-vcf "$CALLED/$name/$name.region.vcf" --variant-nodes "$CALLED/$name/$name.variant_nodes.tsv" \
    --samples "$dir/$name.samples.tsv" --threads 2
  echo "$name: $(grep -E '\] (sites|k-mer features kept|node and edge features kept|variant records|samples): ' \
    "$dir/$name.describe.log" | sed 's/^\[describe [0-9:]*\] //' | paste -sd ';' - | sed 's/;/; /g')"
}

for name in lpa lpa_lossy c4 gstm cyp2d acot myom2_rebuilt; do
  describe "$name"
done
