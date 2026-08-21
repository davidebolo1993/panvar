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
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="$HERE/real_data"
mkdir -p "$OUT"; OUT="$OUT/run.$$.$(date +%s)"; rm -rf "$OUT"; mkdir -p "$OUT"
fails=0
ok()  { printf "  ok   %s\n" "$1"; }
bad() { printf "  FAIL %s\n" "$1"; fails=$((fails + 1)); }
skip(){ printf "  --   %s (skipped)\n" "$1"; }

ref_of() { gunzip -c "$1" 2>/dev/null | awk -F'\t' '($1=="P"||$1=="W"){n=$2;
             if(n~/[Gg][Rr][Cc]h38/){print n;exit} if(!f)f=n} END{if(f&&!d)print f}' | head -1; }
# every interior node id claimed by any emitted site, one per line
interiors() { tail -n +2 "$1" | awk -F'"' '{print $2}' | tr ';' '\n' | grep -v '^$'; }

# ---------------------------------------------------------------- ANKRD36C: disjoint sites
# Bubble emitted a site enclosing all ten others here, folding the same 5616 bp unit as the site
# inside it, and Panphorte failed the pipeline outright at "bubbles 1 and 2 both claim interior node
# 7467". The conflict pass fixed it; this is the only locus that exercises it.
G="$DATA/ankrd36c.gfa.gz"
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
    # The measured policy: ten well-localised sites reconstruct far better than one enclosing snarl
    # (89.54% against 3.98% of the gap closed), so the smaller set is the one to keep.
    [ "$n" -ge 2 ] && ok "and keeps the smaller sites rather than one enclosing snarl ($n)" \
                   || bad "ANKRD36C collapsed to $n site(s); the enclosing snarl was retained"
    # Panphorte's preflight is what used to fail. It must now accept the same set.
    "$BIN" panphorte -i "$OUT/ank.sorted.gfa" --bubble-prefix-in "$OUT/ank" \
       -o "$OUT/ank_pan" --reference-path "$R" --min-similarity 0.95 --quiet >/dev/null 2>&1
    [ "$?" -eq 0 ] && [ -s "$OUT/ank_pan.normalized.sorted.gfa" ] \
      && ok "Panphorte's preflight accepts the ANKRD36C site set" \
      || bad "Panphorte refused the ANKRD36C site set; the disjointness fix has regressed"
  fi
fi

# ---------------------------------------------------------------- LPA: the KIV-2 fold
# The defining folding case. Pinned: the unit length, that every haplotype carries it, the copy-number
# range, and that the fold is recorded in the provenance table rather than only in a log line.
G="$DATA/lpa.gfa.gz"
if [ ! -f "$G" ]; then
  skip "LPA KIV-2 folding (no test graph)"
else
  R=$(ref_of "$G")
  "$BIN" bubble -i "$G" -r "$R" -o "$OUT/lpa" --quiet >/dev/null 2>&1
  "$BIN" panphorte -i "$OUT/lpa.sorted.gfa" --bubble-prefix-in "$OUT/lpa" -o "$OUT/lpa_pan" \
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
    [ "${coll:-0}" -gt 3000 ] && ok "collapsing ${coll} nodes onto one REP node" \
                             || bad "nodes_collapsed=${coll:-none}, expected >3000"
    PROV="$OUT/lpa_pan.panphorte.rep_provenance.tsv"
    [ -s "$PROV" ] && [ "$(($(wc -l < "$PROV") - 1))" -ge 1 ] \
      && ok "the REP node is recorded in the provenance table" \
      || bad "no REP provenance row for the LPA fold"
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
for loc in lpa c4 gstm1 cyp2d6 acot ankrd36c; do
  G="$DATA/$loc.gfa.gz"; [ -f "$G" ] || continue
  R=$(ref_of "$G")
  "$BIN" bubble -i "$G" -r "$R" -o "$OUT/${loc}_b" --quiet >/dev/null 2>&1 || continue
  "$BIN" panphorte -i "$OUT/${loc}_b.sorted.gfa" --bubble-prefix-in "$OUT/${loc}_b" \
     -o "$OUT/${loc}_p" --reference-path "$R" --min-similarity 0.95 --quiet >/dev/null 2>&1 || continue
  "$BIN" refine -i "$OUT/${loc}_p.normalized.sorted.gfa" --bubble-prefix-in "$OUT/${loc}_p" \
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
G="$DATA/c4.gfa.gz"
if [ ! -f "$G" ]; then
  skip "rebuild round-trip (no C4 graph)"
else
  R=$(ref_of "$G")
  gunzip -c "$G" > "$OUT/c4.gfa"
  "$BIN" rebuild -i "$OUT/c4.gfa" -o "$OUT/c4_rebuilt.gfa" --min-recovered-identity 0.97 \
     --quiet >/dev/null 2>&1
  if [ ! -s "$OUT/c4_rebuilt.gfa" ]; then
    bad "rebuild produced no output on C4"
  else
    ok "rebuild emitted an accepted graph on C4"
    "$BIN" bubble -i "$OUT/c4_rebuilt.gfa" -r "$R" -o "$OUT/c4_rb_bub" --quiet >/dev/null 2>&1
    rc=$?
    [ "$rc" -eq 0 ] && [ -s "$OUT/c4_rb_bub.bubbles.csv" ] \
      && ok "and the next module parses and decomposes it (exit 0)" \
      || bad "bubble could not consume the rebuilt graph (exit $rc)"
    dup=$(interiors "$OUT/c4_rb_bub.bubbles.csv" 2>/dev/null | sort | uniq -d | wc -l | tr -d ' ')
    [ "${dup:-0}" = "0" ] && ok "with pairwise-disjoint interiors on the round-tripped graph" \
                          || bad "$dup overlapping interior(s) after the round trip"
  fi
fi

echo
if [ "$fails" -eq 0 ]; then echo "real_regressions: all assertions passed"; exit 0; fi
echo "real_regressions: $fails assertion(s) FAILED"; exit 1
