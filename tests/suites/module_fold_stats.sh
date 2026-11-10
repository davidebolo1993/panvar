#!/usr/bin/env bash
# module_fold_stats.sh - rebuild's module fold (--mode folded, the default).
#
# A synthetic module first: an 8 kb unit U and a near-identical copy U2 between unique flanks; the reference
# holds two copies, h1 one, h2 three, h3 two of U2. Then fixtures in tests/synthetic_data/module_fold_cases/: both module
# boundaries in one node (f1), a haplotype ending inside the module (f2), reverse-stored paths and an inverted
# copy (f3), a truncated copy (f4b), long differences inside copies (f5), a one-copy reference (f6), and an
# indel just before the module's start (f7).
#
#   module_fold_stats.sh <panvar-binary> <out-dir>
set -uo pipefail
BIN="${1:?usage: module_fold_stats.sh <panvar> <outdir>}"
OUT="${2:?}"
TOOLS="$(cd "$(dirname "$0")/.." && pwd)/tools"
mkdir -p "$OUT"; OUT="$OUT/run.$$.$(date +%s)"; rm -rf "$OUT"; mkdir -p "$OUT"
fails=0
ok()  { printf "  ok   %s\n" "$1"; }
bad() { printf "  FAIL %s\n" "$1"; fails=$((fails + 1)); }
check() { [ "$2" = "$3" ] && ok "$1 = $3" || bad "$1: expected $3, got ${2:-<missing>}"; }

cat > "$OUT/fold.spec" <<'SPEC'
N F1 10000 1
N U 8000 2
N U2 @U:100:3
N F2 10000 5
P ref#0#chr1:0-36000 F1+,U+,U2+,F2+
P h1#1#chr1 F1+,U+,F2+
P h2#1#chr1 F1+,U+,U2+,U+,F2+
P h3#1#chr1 F1+,U2+,U2+,F2+
P h4#1#chr1 F1+,U+,U2+,F2+
SPEC
"${PYTHON:-python3}" -I "$TOOLS/spec_gfa.py" "$OUT/fold.spec" "$OUT/in" >/dev/null
REF='ref#0#chr1:0-36000'
"$BIN" rebuild -i "$OUT/in/g.gfa" -o "$OUT/folded.gfa" -r "$REF" --force --allow-loss 2>"$OUT/folded.err"
"$BIN" rebuild -i "$OUT/in/g.gfa" -o "$OUT/plain.gfa" -r "$REF" --force --allow-loss --mode minigraph -q 2>/dev/null

# Per path: does it spell what it did in the input, and how often does it pass the unit's first node.
report() {   # report <input.gfa> <output.gfa>
  "${PYTHON:-python3}" -I - "$1" "$2" <<'PY'
import sys
def rd(p):
    seq, paths, tag = {}, {}, {}
    for l in open(p):
        x = l.rstrip('\n').split('\t')
        if x[0] == 'S':
            seq[x[1]] = x[2]; tag[x[1]] = ' '.join(x[3:])
        elif x[0] == 'P':
            paths[x[1]] = [(s[:-1], s[-1] == '-') for s in x[2].split(',')]
    return seq, paths, tag
C = str.maketrans('ACGT', 'TGCA')
sp = lambda seq, st: ''.join(seq[n].translate(C)[::-1] if r else seq[n] for n, r in st)
a, b = rd(sys.argv[1]), rd(sys.argv[2])
u0 = [n for n, t in b[2].items() if 'mn:Z:fold_unit:0-' in t]
for p in a[1]:
    exact = sp(a[0], a[1][p]) == sp(b[0], b[1][p])
    passes = sum(1 for n, _ in b[1][p] if n in u0) if u0 else 0
    print(f"{p.split('#')[0]}:{'exact' if exact else 'differs'}:{passes}", end=' ')
PY
}
check "folded: every path spells its own sequence, passing the unit once per copy" \
      "$(report "$OUT/in/g.gfa" "$OUT/folded.gfa")" "ref:exact:2 h1:exact:1 h2:exact:3 h3:exact:2 h4:exact:2 "
check "the run reports the module folded, with its copies" \
      "$(grep -o 'folded copy-number modules: 1 (copies: 1-3, reference: 2)' "$OUT/folded.err" | head -1)" "folded copy-number modules: 1 (copies: 1-3, reference: 2)"
check "--mode minigraph leaves the copies apart: no unit nodes, h1's lost copy is written in as a deletion, h3 is spelled A-like" \
      "$(report "$OUT/in/g.gfa" "$OUT/plain.gfa")" "ref:exact:0 h1:exact:0 h2:exact:0 h3:differs:0 h4:exact:0 "
check "the report records the mode and the module, folded with its copies" \
      "$(grep '^#mode\|^#module' "$OUT/folded.rebuild.report.tsv" | cut -f1,2,6,7,8 | tr '\t\n' '  ')" \
      "#mode folded #module 9920-26080 folded copies=1-3 reference_copies=2 "
# Acceptance is judged on what is written: h1's rebuilt walk keeps a lost copy, its folded path does not.
"$BIN" rebuild -i "$OUT/in/g.gfa" -o "$OUT/strict.gfa" -r "$REF" --force -q 2>/dev/null
check "a rebuild whose only failing walk the fold makes exact is accepted without --allow-loss" \
      "$(grep '^#verdict' "$OUT/strict.rebuild.report.tsv" | cut -f2) $(awk -F'\t' '$1 ~ /^h1/ { print $(NF-1) }' "$OUT/strict.rebuild.report.tsv")" \
      "accepted exact"

# The fixtures. Each is rebuilt folded; the module line of its table and every path's spelling are checked.
FIX="$(cd "$(dirname "$0")/.." && pwd)/synthetic_data/module_fold_cases"
fixture() {   # fixture <name>: rebuild it, print its module line's counts and the per-path report
  local d="$OUT/$1"
  "${PYTHON:-python3}" -I "$TOOLS/spec_gfa.py" "$FIX/$1.spec" "$d" >/dev/null
  "$BIN" rebuild -i "$d/g.gfa" -o "$d/out.gfa" -r ref --force --allow-loss -t 1 2>"$d/err"
  report "$d/g.gfa" "$d/out.gfa" >"$d/report"
}
for f in f1 f2 f3 f4b f5 f6 f7; do fixture "$f"; done
# the module's span, copies, reference copies, copy identity, paths folded and those from their own sequence
modline() { awk -F'\t' '$1=="#module"{print $2, $7, $8, $9, $10, $11}' "$OUT/$1/out.rebuild.report.tsv" | paste -sd' ' -; }
check "f1: both boundaries in one node -- the module, not the path's end, is folded; all from own sequence" \
      "$(modline f1)" "14940-33060 copies=3-3 reference_copies=3 copy_identity=0.95 paths_folded=3/3 exact_stretches=3/3"
check "f2: a haplotype ending inside the module keeps its walk; the others fold" \
      "$(modline f2)" "14920-31080 copies=1-3 reference_copies=2 copy_identity=0.95 paths_folded=3/4 exact_stretches=3/3"
check "f2: every path spells its own sequence" "$(cat "$OUT/f2/report")" "ref:exact:2 h1:exact:1 h2:exact:3 h3:exact:0 "
check "f4b: a truncated copy covering most of the unit is a copy of its own, not part of the next" \
      "$(modline f4b)" "14950-22050 copies=1-3 reference_copies=2 copy_identity=0.95 paths_folded=4/4 exact_stretches=4/4"
# Folded, with every spanning path folded and taken from its own sequence; a not_folded line carries a reason in
# those fields instead, and fails.
folded_all() { awk '{split($5,a,"[=/]"); split($6,b,"[=/]")
  print (a[1]=="paths_folded" && b[1]=="exact_stretches" && a[2]>0 && a[2]==a[3] && b[2]==b[3]) ? "all" : $0}'; }
check "folded_all fails a module that was not folded" \
      "$(printf '1-2 no haplotype holds two copies\n' | folded_all)" "1-2 no haplotype holds two copies"
for f in f3 f5 f6; do
  check "$f: every path folded from its own sequence" "$(modline "$f" | folded_all)" "all"
done
check "f7: an indel just before the module's start -- every junction still placed in the haplotype's own sequence" \
      "$(modline f7)" "14920-31080 copies=1-3 reference_copies=2 copy_identity=0.95 paths_folded=4/4 exact_stretches=4/4"
# f7's indel lies just before the module's copies: the span, widened by 1% of the unit, takes it in, so the fold
# handles it in one pass.
check "f7: accepted in one pass, the indel inside the widened module" \
      "$(grep -c 'written in on a second pass' "$OUT/f7/err"):$(awk -F'\t' '$1=="#verdict"{print $2}' "$OUT/f7/out.rebuild.report.tsv")" "0:accepted"
"$BIN" rebuild -i "$OUT/in/g.gfa" -o "$OUT/x.gfa" --mode fold >/dev/null 2>&1
check "an unknown --mode is refused" "$?" "2"
"$BIN" rebuild -i "$OUT/in/g.gfa" -o "$OUT/x.txt" -q >/dev/null 2>&1
check "an output not named .gfa or .gfa.gz is refused" "$?" "2"

# Cases for the repair and the fold's discovery, generated by tests/tools/rebuild_cases.py (see its header).
"${PYTHON:-python3}" -I "$TOOLS/rebuild_cases.py" "$OUT/cases" >/dev/null
rb() {   # rb <case> <label> [options]: rebuild a case into $OUT/cases/<label>
  local c=$1 l=$2; shift 2
  "$BIN" rebuild -i "$OUT/cases/$c.gfa" -o "$OUT/cases/$l.gfa" -r ref --force -t 1 "$@" 2>"$OUT/cases/$l.err"
}
rep() { awk -F'\t' -v k="$2" '$1==k' "$OUT/cases/$1.rebuild.report.tsv" 2>/dev/null; }
modules() { rep "$1" '#module' | awk -F'\t' '$6=="folded"{print $7, $8}' | paste -sd';' -; }
rb unit2999 u2999; rb unit3000 u3000
check "a tandem unit under 3 kb is left to panphorte, not folded at twice its period" "$(modules u2999)" ""
check "a 3 kb unit folds with its own copy counts" "$(modules u3000)" "copies=1-3 reference_copies=2"
rb equal0 e0; rb equal1000 e1000
check "two arrays changing in opposite directions are both found" \
      "$(modules e0)" "copies=1-3 reference_copies=3;copies=1-3 reference_copies=1"
check "an array whose extra copies sit right at its boundary is found and folded" \
      "$(modules e1000)" "copies=1-3 reference_copies=3;copies=1-3 reference_copies=1"
rb adjacent adj
check "an array next to one whose cut runs into its first copy still folds" \
      "$(modules adj)" "copies=1-3 reference_copies=2;copies=1-3 reference_copies=2"
check "  ... and the first module's span ends where the second's begins" \
      "$(rep adj '#module' | awk -F'\t' '$6=="folded"{split($2,s,"-"); if (n++ && s[1] < e) bad=1; e=s[2]} END{print (n==2 && !bad) ? "apart" : "overlap"}')" "apart"
check "  ... and a node the second cut split carries its own range in its mn tag" \
      "$(awk -F'\t' '$1=="S" && match($4, /:[0-9]+-[0-9]+$/) {split(substr($4, RSTART+1), r, "-"); if (r[2]-r[1] != length($3)) print $4}' "$OUT/cases/adj.gfa" | head -1)" ""
rb adjacent_snp adjs
check "  ... also with differences in the flank stretch the second cut takes in" \
      "$(modules adjs)" "copies=1-3 reference_copies=2;copies=1-3 reference_copies=2"
check "  ... every path spelling its own sequence" \
      "$(report "$OUT/cases/adjacent_snp.gfa" "$OUT/cases/adjs.gfa" | sed 's/:[0-9]* /:/g')" "ref:exact:h1:exact:h2:exact:h3:exact:h4:exact:"
rb adjacent600 a600
check "an array whose cut would take in a folded array's copies is refused" \
      "$(rep a600 '#module' | cut -f6,7 | paste -sd';' -)" "$(printf 'folded\tcopies=1-3;not_folded\toverlaps a module already folded')"
check "  ... and the folded array keeps its unit, each path passing it once per copy" \
      "$(report "$OUT/cases/adjacent600.gfa" "$OUT/cases/a600.gfa")" "ref:exact:2 h1:exact:1 h2:exact:3 h3:exact:2 h4:exact:3 "
for c in fragment:"37% of a unit" fragment_half:"52% of a unit" fragment_full:"a full unit"; do
  rb "${c%%:*}" "${c%%:*}.out"
  check "a copy of ${c#*:} far past its array does not stretch the module over the array between" \
        "$(modules "${c%%:*}.out")" "copies=1-3 reference_copies=1;copies=1-3 reference_copies=2"
done
rb fragment_up frag_up
check "  ... nor does one before it" "$(modules frag_up)" "copies=1-3 reference_copies=1;copies=1-3 reference_copies=2"
rb fragment_long frag_long
check "  ... nor one closer than a haplotype's long array is long" \
      "$(modules frag_long)" "copies=1-10 reference_copies=1;copies=1-3 reference_copies=2"
rb fragment_pair frag_pair
check "  ... and two copies far past it do not take the module from the array" \
      "$(modules frag_pair)" "copies=1-3 reference_copies=1;copies=1-3 reference_copies=2"
rb diverged div
check "reference copies at 92.5% identity are one module" "$(modules div)" "copies=1-3 reference_copies=2"
rb second_pass sp
check "a copy loss the fold does not take is written in on a second pass, and the rebuild is accepted" \
      "$(grep -c 'written in on a second pass' "$OUT/cases/sp.err"):$(rep sp '#verdict' | cut -f2)" "1:accepted"
rb second_pass cl --mode minigraph
check "--mode minigraph: a lost copy is one edge past it, not a node copying graph sequence" \
      "$(grep -c 'mn:Z:insertion' "$OUT/cases/cl.gfa"):$(rep cl '#verdict' | cut -f2)" "0:accepted"
rb inflate inf --mode minigraph --min-matched-cover 0.999
check "a repair counts only the bases it makes match: 100 collapsed differences keep matched cover under 0.999" \
      "$(rep inf alt | cut -f5,17):$(rep inf '#verdict' | cut -f2)" "$(printf '0.9979381443\tlow_cover'):rejected"
rb endins ei
check "insertions the generator left out near the end are written in, though the generated graph has no variation" \
      "$(rep ei '#indels_added' | cut -f2):$(rep ei '#verdict' | cut -f2)" "2:accepted"
rb overrun_nospan on
check "a fold elsewhere does not pass a walk that runs past its haplotype's end" \
      "$(rep on h1 | cut -f16,17):$(rep on '#verdict' | cut -f2)" "$(printf '.\tend_overrun'):rejected"
rb del_flat df
check "deletions written into a graph with no other variation are kept, not passed through" \
      "$(rep df '#indels_added' | cut -f2):$(rep df '#verdict' | cut -f2)" "2:accepted"
"$BIN" rebuild -i "$OUT/in/g.gfa" -o "$OUT/gz.gfa.gz" -r "$REF" --force -q 2>/dev/null
check "the report is named from the output's prefix" \
      "$([ -s "$OUT/gz.rebuild.report.tsv" ] && [ ! -e "$OUT/gz.gfa.gz.rebuild.report.tsv" ] && echo yes || echo no)" "yes"

echo
if [ "$fails" -eq 0 ]; then echo "module_fold_stats: all assertions passed"; exit 0; fi
echo "module_fold_stats: $fails assertion(s) FAILED"; exit 1
