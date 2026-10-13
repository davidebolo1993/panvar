#!/usr/bin/env bash
# Truth isolation gate (docs/algorithms/genotype.md, stage 11): the dosage audit's sources and headers read no
# truth or evaluation input. They take no path argument other than the run prefix they audit and the prefix they write;
# no truth GFA, evaluation table or truth-named field appears in them. Usage: <source root>
set -uo pipefail
ROOT=${1:-.}
files=$(ls "$ROOT"/src/genotype_dosage*.cpp "$ROOT"/include/panvar/genotype_dosage*.hpp)
fail=0
for f in $files; do
  hits=$(grep -niE 'truth|\.eval\.|eval\.blocks|evaluator|eval_summary|truth_gfa|posthoc|phased' "$f" | grep -viE 'truth-free|never.*truth|no truth|truth is|evaluated in this order|evaluated in parallel|is evaluated' || true)
  if [ -n "$hits" ]; then echo "FAIL $f:"; echo "$hits"; fail=1; fi
done
# the only files the audit opens are the run's own outputs
opens=$(grep -hoE '(run_prefix|out_prefix) \+ "\.[a-z_.]+' $files | sed 's/.*"//' | sort -u | tr '\n' ' ')
echo "files the audit reads or writes (by suffix): $opens"
echo "$opens" | grep -qE 'eval|truth' && { echo "FAIL a truth or evaluation file is referenced"; fail=1; }
[ $fail -eq 0 ] && echo "PASS truth isolation ($(echo $files | wc -w | tr -d ' ') files)" || echo FAILED
exit $fail
