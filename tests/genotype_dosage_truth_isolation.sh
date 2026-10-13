#!/usr/bin/env bash
# Truth isolation gate (docs/algorithms/genotype.md, stage 11): the dosage audit's sources and headers read no
# truth or evaluation input. They take no path argument other than the run prefix they audit and the prefix they write;
# no truth GFA, evaluation table or truth-named field appears in them. Usage: <source root>
# The scan covers executable source and string literals; C/C++ comments are removed first (line numbers kept), so a
# comment that merely mentions an evaluator or a known-truth gate is not a finding, while any code or literal is.
set -uo pipefail
ROOT=${1:-.}
files=$(ls "$ROOT"/src/genotype_dosage*.cpp "$ROOT"/include/panvar/genotype_dosage*.hpp)
fail=0
# Prints the file with // and /* */ comments blanked out, string and character literals kept, and every newline kept.
strip_comments() {
  python3 - "$1" <<'PY'
import sys
s = open(sys.argv[1], encoding="utf-8", errors="replace").read()
out, i, n = [], 0, len(s)
while i < n:
    c = s[i]
    if c in "\"'":
        j = i + 1
        while j < n and s[j] != c:
            j += 2 if s[j] == "\\" else 1
        out.append(s[i:j + 1]); i = j + 1
    elif s.startswith("//", i):
        while i < n and s[i] != "\n": i += 1
    elif s.startswith("/*", i):
        j = s.find("*/", i + 2); j = n if j < 0 else j + 2
        out.append("".join(ch if ch == "\n" else " " for ch in s[i:j])); i = j
    else:
        out.append(c); i += 1
sys.stdout.write("".join(out))
PY
}
for f in $files; do
  hits=$(strip_comments "$f" | grep -niE 'truth|\.eval\.|eval\.blocks|evaluator|eval_summary|truth_gfa|posthoc|phased' | grep -viE 'truth-free|never.*truth|no truth|truth is|evaluated in this order|evaluated in parallel|is evaluated' || true)
  if [ -n "$hits" ]; then echo "FAIL $f:"; echo "$hits"; fail=1; fi
done
# the only files the audit opens are the run's own outputs
opens=$(for f in $files; do strip_comments "$f"; done | grep -hoE '(run_prefix|out_prefix) \+ "\.[a-z_.]+' | sed 's/.*"//' | sort -u | tr '\n' ' ')
echo "files the audit reads or writes (by suffix): $opens"
echo "$opens" | grep -qE 'eval|truth' && { echo "FAIL a truth or evaluation file is referenced"; fail=1; }
[ $fail -eq 0 ] && echo "PASS truth isolation ($(echo $files | wc -w | tr -d ' ') files)" || echo FAILED
exit $fail
