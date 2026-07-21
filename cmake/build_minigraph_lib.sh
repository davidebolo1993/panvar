#!/usr/bin/env bash
# Build libminigraph.a for panvar.
#
# minimap2 vendored the same lh3 chaining/util code as minigraph, so the two export ~24 identical symbols
# (mg_lchain_dp, gfa_aux_*, seq_nt4_table, ksort/krmq instantiations) and cannot be linked side by side.
# We partial-link (ld -r) every minigraph object except main.o into one relocatable object, keeping only
# the entry points `rebuild` calls global and localizing the rest: minigraph's internal calls then bind to
# its own copies and nothing collides with minimap2.
#
# Usage: build_minigraph_lib.sh <minigraph_dir> <exports_file>
set -euo pipefail

dir="$1"
exports="$2"
cd "$dir"

objs=$(ls ./*.o | grep -v '^\./main\.o$')
if [ -z "$objs" ]; then
    echo "build_minigraph_lib.sh: no objects in $dir (did 'make minigraph' run?)" >&2
    exit 1
fi

rm -f libminigraph.a minigraph_merged.o

# shellcheck disable=SC2086  # word splitting of the object list is intended
if [ "$(uname -s)" = "Darwin" ]; then
    ld -r -o minigraph_merged.o $objs -exported_symbols_list "$exports"
else
    ld -r -o minigraph_merged.o $objs
    objcopy --keep-global-symbols="$exports" minigraph_merged.o
fi

ar rcs libminigraph.a minigraph_merged.o
