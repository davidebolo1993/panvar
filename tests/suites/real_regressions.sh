#!/usr/bin/env bash
# real_regressions.sh - opt-in regressions against the committed real loci.
#
# The per-module *_stats.sh suites use hand-sized fixtures, which is right: every expected number is
# derivable rather than recorded. But several contracts are only reachable on real data, and the
# review ledger names them one by one:
#
#   * ANKRD36C is the locus where Bubble used to emit overlapping sites and Panphorte then refused the
#     whole run. Nothing synthetic reproduces a ten-site nested tangle.
#   * LPA KIV-2 is the folding case the whole project is built around, and no synthetic fixture has its
#     466 haplotypes or its 5,547 bp unit.
#   * The POA guards are resource bounds. Whether they SKIP anything is a property of real sequence,
#     and "no region was skipped" is a different claim from "the guards are implemented".
#   * An accepted rebuilt GFA has to be re-readable by the next module. In-memory acceptance and
#     emitted-file semantics are two different things.
#
# Registered under the `real` label and PANVAR_SLOW_TESTS, because it runs the real pipeline.
#
#   real_regressions.sh <panvar-binary> <out-dir>
set -uo pipefail

BIN="${1:?usage: real_regressions.sh <panvar> <outdir>}"
OUT="${2:?}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA="$HERE/real_data"
mkdir -p "$OUT"; OUT="$OUT/run.$$.$(date +%s)"; rm -rf "$OUT"; mkdir -p "$OUT"
fails=0
ok()  { printf "  ok   %s\n" "$1"; }
bad() { printf "  FAIL %s\n" "$1"; fails=$((fails + 1)); }
skip(){ printf "  --   %s (skipped)\n" "$1"; }

ref_of() { gunzip -c "$1" 2>/dev/null | awk -F'\t' '($1=="P"||$1=="W"){n=$2;
             if(n~/[Gg][Rr][Cc]h38/){print n;exit} if(!f)f=n} END{if(f&&!d)print f}' | head -1; }
# every interior node id claimed by an emitted site, one per line. An inversion site holds the sites
# nested in it (their parent_id), so it is left out here and checked by outside_parent instead.
interiors() { python3 -c '
import csv, sys
for r in csv.DictReader(open(sys.argv[1])):
    if r["kind"] != "inversion": print("\n".join(n for n in r["inside_nodes"].split(";") if n))' "$1"; }
# number of sites that overlap an inversion site without one being nested in the other: a node of
# either, boundary or interior, in the other's interior. Only a shared boundary is allowed.
# interiors() leaves inversion sites out, so an overlap with one is only caught here.
inversion_overlaps() { python3 -c '
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))
inside = {r["bubble_id"]: set(filter(None, r["inside_nodes"].split(";"))) for r in rows}
ends = {r["bubble_id"]: {r["source"], r["sink"]} for r in rows}
parent = {r["bubble_id"]: r["parent_id"] for r in rows}
def ancestors(b):
    seen = set()
    while parent.get(b, "") not in ("", "0") and b not in seen:
        seen.add(b); b = parent[b]
        yield b
n = 0
for inv in (r["bubble_id"] for r in rows if r["kind"] == "inversion"):
    for b in inside:
        if b == inv or inv in ancestors(b) or b in ancestors(inv): continue
        if (inside[b] | ends[b]) & inside[inv] or (inside[inv] | ends[inv]) & inside[b]:
            n += 1
print(n)' "$1"; }
# number of nested sites reaching outside their parent: an interior node outside the parent's
# interior, or a boundary outside the parent's interior and boundaries
outside_parent() { python3 -c '
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))
inside = {r["bubble_id"]: set(filter(None, r["inside_nodes"].split(";"))) for r in rows}
ends = {r["bubble_id"]: {r["source"], r["sink"]} for r in rows}
def outside(b, p):
    return not inside[b] <= inside.get(p, set()) or not ends[b] <= inside.get(p, set()) | ends.get(p, set())
print(sum(1 for r in rows if r["parent_id"] not in ("", "0") and outside(r["bubble_id"], r["parent_id"])))' "$1"; }

# The checks above on hand-made tables, so a check that cannot fail is caught: an inversion 1..5, and a
# site sharing its sink (allowed), starting inside it (an overlap), or nested in it but ending outside.
H=bubble_id,source,sink,source_orient,sink_orient,inside_nodes,kind,parent_id
printf '%s\n' "$H" '1,1,5,+,+,2;3;4,inversion,0' '2,5,7,+,+,6,site,0' > "$OUT/check_shared.csv"
printf '%s\n' "$H" '1,1,5,+,+,2;3;4,inversion,0' '2,4,6,+,+,5,site,0' > "$OUT/check_overlap.csv"
printf '%s\n' "$H" '1,1,5,+,+,2;3;4,inversion,0' '2,4,6,+,+,,site,1' > "$OUT/check_nested.csv"
got="$(inversion_overlaps "$OUT/check_shared.csv")$(outside_parent "$OUT/check_shared.csv")"
got="$got $(inversion_overlaps "$OUT/check_overlap.csv")$(outside_parent "$OUT/check_overlap.csv")"
got="$got $(inversion_overlaps "$OUT/check_nested.csv")$(outside_parent "$OUT/check_nested.csv")"
[ "$got" = "00 10 01" ] && ok "the overlap checks pass a shared boundary and catch a boundary inside an inversion" \
                       || bad "the overlap checks give '$got' on hand-made tables (expected '00 10 01')"

# ---------------------------------------------------------------- ANKRD36C: disjoint sites
# Bubble emitted a site enclosing all ten others here, folding the same 5616 bp unit as the site
# inside it, and Panphorte failed the pipeline outright at "bubbles 1 and 2 both claim interior node
# 7467". The conflict pass fixed it; this is the only locus that exercises it.
G="$DATA/inversion/ankrd36c.gfa.gz"
if [ ! -f "$G" ]; then
  skip "ANKRD36C disjointness (no test graph)"
else
  R=$(ref_of "$G")
  "$BIN" bubble -i "$G" -r "$R" -o "$OUT/ank" --quiet >/dev/null 2>&1
  if [ ! -s "$OUT/ank.bubbles.csv" ]; then
    bad "ANKRD36C produced no bubbles CSV"
  else
    n=$(($(wc -l < "$OUT/ank.bubbles.csv") - 1))
    dup=$(interiors "$OUT/ank.bubbles.csv" | sort | uniq -d | wc -l | tr -d ' ')
    [ "$dup" = "0" ] && ok "ANKRD36C emits pairwise-disjoint interiors across $n sites" \
                     || bad "$dup interior node(s) claimed by more than one ANKRD36C site"
    out=$(outside_parent "$OUT/ank.bubbles.csv")
    [ "$out" = "0" ] && ok "and every nested site lies inside its parent" \
                     || bad "$out nested ANKRD36C site(s) reach outside their parent"
    out=$(inversion_overlaps "$OUT/ank.bubbles.csv")
    [ "$out" = "0" ] && ok "and every site touching an inversion site's interior is nested in it" \
                     || bad "$out ANKRD36C site(s) partly overlap an inversion site"
    # The measured policy: ten well-localised sites reconstruct far better than one enclosing snarl
    # (89.54% against 3.98% of the gap closed), so the smaller set is the one to keep.
    [ "$n" -ge 2 ] && ok "and keeps the smaller sites rather than one enclosing snarl ($n)" \
                   || bad "ANKRD36C collapsed to $n site(s); the enclosing snarl was retained"
    # GRCh38 walks an inverted duplication here, so the snarl decomposition roots in a cycle and one
    # top-level snarl closes it, wrapping the locus ends: every site upstream was nested in it and none
    # was reported (the first site started 26% into the reference). Its contents are now emitted in
    # its place (first site at 11%).
    first=$(awk -F'\t' -v r="$R" 'NR==FNR{if($1=="P" && $2==r){n=split($3,a,","); for(i=1;i<=n;i++){id=substr(a[i],1,length(a[i])-1); if(!(id in p)) p[id]=i/n}} next}
                 FNR>1{split($0,f,","); if((f[2] in p) && (m=="" || p[f[2]]<m)) m=p[f[2]]} END{printf "%.2f", m}' \
                 "$OUT/ank.sorted.gfa" "$OUT/ank.bubbles.csv")
    awk -v f="$first" 'BEGIN{exit !(f < 0.20)}' \
      && ok "and reports sites upstream too (first site at $first of the reference)" \
      || bad "ANKRD36C's first site lies at $first of the reference; the region upstream has none"
    # The same inverted duplication makes the locus's last 81 kb a snarl every path enters through
    # the reverse visit of a shared node and ends in. None crosses it, so it is no site, and the sites
    # nested in it were lost with it (the last site sat at 0.88 of the walk, before a 131 bp deletion that 411 of
    # 464 haplotypes carry at 0.86). A snarl holding an end of most paths is now descended into too.
    last=$(awk -F'\t' -v r="$R" 'NR==FNR{if($1=="P" && $2==r){n=split($3,a,","); for(i=1;i<=n;i++){id=substr(a[i],1,length(a[i])-1); p[id]=i/n}} next}
                 FNR>1{split($0,f,","); if((f[2] in p) && p[f[2]]>m) m=p[f[2]]} END{printf "%.2f", m}' \
                 "$OUT/ank.sorted.gfa" "$OUT/ank.bubbles.csv")
    awk -v f="$last" 'BEGIN{exit !(f > 0.90)}' \
      && ok "and reports sites at the locus end too (last site at $last of the reference)" \
      || bad "ANKRD36C's last site lies at $last of the reference; the locus end has none"
    # Panphorte's preflight is what used to fail. It must now accept the same set.
    "$BIN" panphorte -i "$OUT/ank.sorted.gfa" --bubble-prefix "$OUT/ank" \
       -o "$OUT/ank_pan" --reference-path "$R" --min-similarity 0.95 --quiet >/dev/null 2>&1
    [ "$?" -eq 0 ] && [ -s "$OUT/ank_pan.normalized.sorted.gfa" ] \
      && ok "Panphorte's preflight accepts the ANKRD36C site set" \
      || bad "Panphorte refused the ANKRD36C site set; the disjointness fix has regressed"
  fi
fi

# ---------------------------------------------------------------- LPA: the KIV-2 fold
# The defining folding case. Pinned: the unit length, that every haplotype carries it, the copy-number
# range, and that the fold is recorded in the provenance table rather than only in a log line.
G="$DATA/tandem_repeat/lpa.gfa.gz"
if [ ! -f "$G" ]; then
  skip "LPA KIV-2 folding (no test graph)"
else
  R=$(ref_of "$G")
  "$BIN" bubble -i "$G" -r "$R" -o "$OUT/lpa" --quiet >/dev/null 2>&1
  "$BIN" panphorte -i "$OUT/lpa.sorted.gfa" --bubble-prefix "$OUT/lpa" -o "$OUT/lpa_pan" \
     --reference-path "$R" --min-similarity 0.95 --quiet >/dev/null 2>&1
  REP="$OUT/lpa_pan.panphorte.report.tsv"
  if [ ! -s "$REP" ]; then
    bad "LPA panphorte wrote no report"
  else
    read -r unit paths mn mx coll <<EOF
$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next} $(c["normalized"])=="yes"{
    print $(c["unit_bp"]), $(c["paths_normalized"]), $(c["min_copies"]), $(c["max_copies"]), $(c["nodes_collapsed"]); exit}' "$REP")
EOF
    [ "${unit:-0}" = "5547" ] && ok "LPA folds the 5,547 bp KIV-2 unit" \
                             || bad "LPA unit_bp=${unit:-none}, expected 5547"
    [ "${paths:-0}" = "466" ] && ok "and every one of the 466 haplotypes is normalized" \
                             || bad "paths_normalized=${paths:-none}, expected 466"
    [ "${mn:-0}" = "1" ] && [ "${mx:-0}" = "32" ] \
      && ok "with copy numbers spanning 1-32" \
      || bad "copies ${mn:-?}-${mx:-?}, expected 1-32"
    # nodes_collapsed counts nodes the fold removed from the graph. It once counted every node a folded
    # copy range touched, with repeats, and read 3485 here.
    [ "${coll:-0}" = "1508" ] && ok "collapsing ${coll} nodes onto one unit graph" \
                             || bad "nodes_collapsed=${coll:-none}, expected 1508"
    PROV="$OUT/lpa_pan.panphorte.unit_graph.tsv"
    [ -s "$PROV" ] && [ "$(($(wc -l < "$PROV") - 1))" -ge 1 ] \
      && ok "the fold is recorded in the unit graph table" \
      || bad "no unit graph row for the LPA fold"
    # Folding must not change what the haplotypes spell beyond the declared similarity.
    [ -s "$OUT/lpa_pan.panphorte.copies.tsv" ] \
      && [ "$(($(wc -l < "$OUT/lpa_pan.panphorte.copies.tsv") - 1))" -ge 466 ] \
      && ok "and a per-haplotype copy count is emitted for every path" \
      || bad "the copies table does not cover all 466 haplotypes"
  fi
fi

# ---------------------------------------------------------------- Refine: POA guard coverage
# The guards are resource bounds, and "safer accounting" must not be mistaken for unchanged coverage.
# This REPORTS what each guard skipped per locus and fails only if refinement did nothing at all
# anywhere -- the number itself is the deliverable, not a threshold.
echo "  ---- refine POA guard coverage (reported, not thresholded) ----"
refined_any=0
for loc in lpa c4 gstm cyp2d acot ankrd36c; do
  G=$(echo "$DATA"/*/"$loc".gfa.gz); [ -f "$G" ] || continue
  R=$(ref_of "$G")
  "$BIN" bubble -i "$G" -r "$R" -o "$OUT/${loc}_b" --quiet >/dev/null 2>&1 || continue
  "$BIN" panphorte -i "$OUT/${loc}_b.sorted.gfa" --bubble-prefix "$OUT/${loc}_b" \
     -o "$OUT/${loc}_p" --reference-path "$R" --min-similarity 0.95 --quiet >/dev/null 2>&1 || continue
  "$BIN" refine -i "$OUT/${loc}_p.normalized.sorted.gfa" --bubble-prefix "$OUT/${loc}_p" \
     --reference-path "$R" -o "$OUT/${loc}_r" --quiet >/dev/null 2>&1 || continue
  RR="$OUT/${loc}_r.refine.report.tsv"; [ -s "$RR" ] || continue
  printf '       %-9s %s\n' "$loc" "$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
      { d=$(c["decision"]); n[d]++ } END { for (k in n) printf "%s=%d ", k, n[k] }' "$RR")"
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next} $(c["decision"]) ~ /REBUILT|rebuilt/ {found=1}
      END{exit !found}' "$RR" && refined_any=1
done
[ "$refined_any" = "1" ] && ok "refinement rebuilds at least one region on the real loci" \
                         || bad "no locus had a single region rebuilt; the guards may be skipping everything"

# ---------------------------------------------------------------- Rebuild: downstream round-trip
# In-memory acceptance and emitted-file semantics are different things. An accepted rebuilt GFA must
# parse again and decompose, or the contract only holds inside the process that wrote it.
G="$DATA/paralogs/c4.gfa.gz"
if [ ! -f "$G" ]; then
  skip "rebuild round-trip (no C4 graph)"
else
  R=$(ref_of "$G")
  gunzip -c "$G" > "$OUT/c4.gfa"
  # --force, because this graph passes the pathology gate. Without it rebuild copies the input
  # through unchanged and the round trip below would prove only that bubble can read its own input.
  # --base-align and -r: without base alignment one haplotype misses a 5.6 kb indel and the length
  # bound refuses the rebuild; with -r the reference seeds it and comes back exactly.
  "$BIN" rebuild -i "$OUT/c4.gfa" -o "$OUT/c4_rebuilt.gfa" --force --base-align -r "$R" \
     --quiet >/dev/null 2>&1
  AUD="$OUT/c4_rebuilt.rebuild.report.tsv"
  if [ ! -s "$OUT/c4_rebuilt.gfa" ]; then
    bad "rebuild produced no output"
  else
    # The report is what proves a rebuild actually happened and was accepted, rather than the input
    # being passed through: a pass-through writes no report at all.
    if [ ! -s "$AUD" ]; then
      bad "no rebuild report was written, so the graph was passed through rather than rebuilt"
    else
      verdict=$(awk -F'\t' '$1=="#verdict"{print $2}' "$AUD")
      nrows=$(awk -F'\t' '$0!~/^#/ && NR>1' "$AUD" | wc -l | tr -d ' ')
      nok=$(awk -F'\t' '$0!~/^#/ && NR>1 && $NF=="ok"' "$AUD" | wc -l | tr -d ' ')
      [ "$verdict" = "accepted" ] && ok "the rebuild is accepted by its own contract" \
                                  || bad "rebuild verdict is '${verdict:-none}', expected accepted"
      [ "$nrows" -gt 0 ] && [ "$nok" = "$nrows" ] \
        && ok "and every one of the $nrows haplotypes is recovered within bounds" \
        || bad "$nok of $nrows haplotypes recovered ok"
    fi
    # ...and the emitted graph must genuinely differ from the input, or nothing was rebuilt.
    in_nodes=$(grep -c "^S" "$OUT/c4.gfa"); out_nodes=$(grep -c "^S" "$OUT/c4_rebuilt.gfa")
    [ "$in_nodes" != "$out_nodes" ] \
      && ok "the emitted graph differs from the input ($in_nodes -> $out_nodes nodes)" \
      || bad "the emitted graph has the same node count as the input; it was passed through"
    "$BIN" bubble -i "$OUT/c4_rebuilt.gfa" -r "$R" -o "$OUT/c4_rb_bub" --quiet >/dev/null 2>&1
    rc=$?
    [ "$rc" -eq 0 ] && [ -s "$OUT/c4_rb_bub.bubbles.csv" ] \
      && ok "and the next module parses and decomposes it (exit 0)" \
      || bad "bubble could not consume the rebuilt graph (exit $rc)"
    dup=$(interiors "$OUT/c4_rb_bub.bubbles.csv" 2>/dev/null | sort | uniq -d | wc -l | tr -d ' ')
    [ "${dup:-0}" = "0" ] && ok "with pairwise-disjoint interiors on the round-tripped graph" \
                          || bad "$dup overlapping interior(s) after the round trip"
    out=$(inversion_overlaps "$OUT/c4_rb_bub.bubbles.csv" 2>/dev/null)
    [ "${out:-0}" = "0" ] && ok "and no site partly overlaps an inversion site there" \
                          || bad "$out site(s) partly overlap an inversion site after the round trip"
  fi
fi

echo
if [ "$fails" -eq 0 ]; then echo "real_regressions: all assertions passed"; exit 0; fi
echo "real_regressions: $fails assertion(s) FAILED"; exit 1
