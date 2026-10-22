#!/usr/bin/env bash
# M15-K: the build-time source provenance (cmake/source_provenance.cmake), in a temporary git repository: a clean commit, a
# tracked edit (dirty), a new commit (a new revision with no reconfigure), an unchanged revision (the header is not
# rewritten, so nothing recompiles), and a directory outside git (unavailable / unknown).
#   source_provenance_test.sh <cmake> <source_provenance.cmake> <scratch dir>
set -uo pipefail
CMAKE="${1:?}"; SCRIPT="${2:?}"; OUT="${3:?}/run.$$"
rm -rf "$OUT"; mkdir -p "$OUT/repo" "$OUT/nogit"
fails=0
ok()  { printf "  ok   %s\n" "$1"; }
bad() { printf "  FAIL %s\n" "$1"; fails=$((fails + 1)); }
gen() { "$CMAKE" -DSOURCE_DIR="$1" -DOUTPUT="$2" -P "$SCRIPT" >/dev/null 2>&1; }
field() { sed -n "s/^#define $2 \"\\(.*\\)\"$/\\1/p" "$1"; }
cd "$OUT/repo" || exit 1
git init -q . && git config user.email t@example.org && git config user.name t
echo one > tracked.txt && git add tracked.txt && git commit -q -m one
H="$OUT/h1.h"
gen "$OUT/repo" "$H"
[ "$(field "$H" PANVAR_SOURCE_COMMIT)" = "$(git rev-parse HEAD)" ] && [ "$(field "$H" PANVAR_SOURCE_TREE)" = "clean" ] \
  && ok "a clean checkout records its HEAD and 'clean'" || bad "clean checkout: $(cat "$H")"
first=$(git rev-parse HEAD)
echo two > tracked.txt
gen "$OUT/repo" "$H"
[ "$(field "$H" PANVAR_SOURCE_COMMIT)" = "$first" ] && [ "$(field "$H" PANVAR_SOURCE_TREE)" = "dirty" ] \
  && ok "a tracked edit records 'dirty' beside the same commit" || bad "dirty checkout: $(cat "$H")"
echo untracked > new.txt
git commit -q -am two
gen "$OUT/repo" "$H"
[ "$(field "$H" PANVAR_SOURCE_COMMIT)" = "$(git rev-parse HEAD)" ] && [ "$(field "$H" PANVAR_SOURCE_COMMIT)" != "$first" ] &&
  [ "$(field "$H" PANVAR_SOURCE_TREE)" = "clean" ] \
  && ok "a new commit is recorded on the next run, with no reconfigure (untracked files ignored)" || bad "new commit: $(cat "$H")"
before=$(ls -li "$H" | awk '{print $1, $7, $8, $9}'); sum=$(cksum < "$H")
sleep 1
gen "$OUT/repo" "$H"
[ "$(ls -li "$H" | awk '{print $1, $7, $8, $9}')" = "$before" ] && [ "$(cksum < "$H")" = "$sum" ] \
  && ok "an unchanged revision leaves the header untouched (nothing recompiles)" || bad "the header was rewritten"
H2="$OUT/h2.h"
# The scratch directory may itself lie inside a git checkout (the build tree); stop git's upward search at it.
GIT_CEILING_DIRECTORIES="$OUT" gen "$OUT/nogit" "$H2"
[ "$(field "$H2" PANVAR_SOURCE_COMMIT)" = "unavailable" ] && [ "$(field "$H2" PANVAR_SOURCE_TREE)" = "unknown" ] \
  && ok "outside a git checkout: unavailable / unknown, never a guess" || bad "no git: $(cat "$H2")"
echo
if [ "$fails" -eq 0 ]; then echo "source_provenance: all assertions passed"; exit 0; fi
echo "source_provenance: $fails assertion(s) FAILED"; exit 1
