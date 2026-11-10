#!/usr/bin/env bash
# call_regress.sh - ORIGINAL vs CURRENT `panvar call` on the committed real loci, field by field.
#
# A byte comparison against the previous snapshot is what hid the C4 coverage-route regression: the
# snapshot already carried the change. This compares against a pinned ORIGINAL binary instead, and
# only `call` differs between the two sides: the inputs are built ONCE, with the current binary
# (bubble -> panphorte --min-similarity 0.95 on the GRCh38 path, as real_regressions.sh does, then
# refine for the second input set), and both binaries call the same files.
#
#   inputs   base     <loc>_p: bubble + panphorte
#            refined  <loc>_r: refine on base (refine is part of the pipeline, and it renumbers)
#            altref   <loc>_a: bubble + panphorte along the FIRST P path of the .gfa.gz (W if no P),
#                     not GRCh38, still called with -r GRCh38: the reference runs through bubbles
#                     that were opened on another path, often on the minus strand, which the
#                     GRCh38-built sets never exercise
#   configs  default | cn (--cn)
#            classify (--cn --classify-ins); a run that aborts has its exit status and stated reason
#            compared, its records when both run
#
# Per run: tools/vcfdiff.py on the region VCF, tools/nodesdiff.py on
# variant_nodes.tsv, tools/cnrecords.py on the CURRENT region VCF (no bubble may carry more than one
# DUP record), then the CURRENT binary's `describe --only-variant` on both VCFs and
# tools/bimbamdiff.py on the dosages. Per input graph: tools/crossings.py, which must find no path
# with two traversals of one bubble (source..sink or sink..source with only interior nodes between;
# CNACCT-8).
#
# Every difference must be listed in call_regress.expected.tsv, with the pin of its exact value:
# 'approved', or 'pending' (printed under PENDING USER APPROVAL, does not fail). Anything else fails,
# and so does a listed difference whose value moved on. Header-line changes and ID renames are
# listed and counted and never fail by themselves. tools/regress_report.py does the matching and
# prints the report; `regress_report.py <out>/run.*/items call_regress.expected.tsv --emit-expected`
# prints the file back with this run's pins and a pending row per new difference.
#
# The reference path is pinned: the P-line named GRCh38#0..., else the first P/W name matching grch38
# case-insensitively, read from each set's own graph. Agents that picked different paths reported
# different aborts.
#
#   call_regress.sh <orig-panvar> <cur-panvar> <out-dir>
#
# Environment:
#   LOCI       default "c4 cyp2d ankrd36c gstm lpa acot" (myom2 is excluded on purpose)
#   CONFIGS    default "default cn classify"
#   SETS       default "base refined altref"
#   INPUT_DIR  reuse prepared base inputs <INPUT_DIR>/<loc><INPUT_SUFFIX>{.normalized.sorted.gfa,
#              .bubbles.csv} instead of building them (INPUT_SUFFIX default _grch38_p0.95)
#   REFINED_DIR / REFINED_SUFFIX  likewise for the refined set (default: refine the base inputs)
#   ALTREF_DIR / ALTREF_SUFFIX  likewise for the altref set (default suffix _first_p0.95)
#   CHECK_INPUTS_DIR / CHECK_INPUTS_SUFFIX / CHECK_ALTREF_SUFFIX  when base or altref inputs are built
#              here, cmp them against <dir>/<loc><suffix>.* (suffixes default _grch38_p0.95 and
#              _first_p0.95)
#   EXPECTED   expected-diff file (default: call_regress.expected.tsv next to this script)
#   DESCRIBE   1 (default) runs the BIMBAM comparison; 0 skips it
#   PARALLEL   1 (default) runs the loci concurrently; 0 runs them one after another
set -uo pipefail

ORIG="${1:?usage: call_regress.sh <orig-panvar> <cur-panvar> <outdir>}"
CUR="${2:?usage: call_regress.sh <orig-panvar> <cur-panvar> <outdir>}"
OUT="${3:?usage: call_regress.sh <orig-panvar> <cur-panvar> <outdir>}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TOOLS="$HERE/tools"
DATA="$HERE/real_data"
# locus_gfa <loc>: the locus graph, from whichever sub-folder of real_data holds it.
locus_gfa() { set -- "$DATA"/*/"$1".gfa.gz; echo "$1"; }
LOCI="${LOCI:-c4 cyp2d ankrd36c gstm lpa acot}"
CONFIGS="${CONFIGS:-default cn classify}"
SETS="${SETS:-base refined altref}"
INPUT_SUFFIX="${INPUT_SUFFIX:-_grch38_p0.95}"
REFINED_SUFFIX="${REFINED_SUFFIX:-_new3}"
ALTREF_SUFFIX="${ALTREF_SUFFIX:-_first_p0.95}"
CHECK_INPUTS_SUFFIX="${CHECK_INPUTS_SUFFIX:-_grch38_p0.95}"
CHECK_ALTREF_SUFFIX="${CHECK_ALTREF_SUFFIX:-_first_p0.95}"
EXPECTED="${EXPECTED:-$HERE/suites/call_regress.expected.tsv}"
DESCRIBE="${DESCRIBE:-1}"
PARALLEL="${PARALLEL:-1}"
PY="python3 -I"

for b in "$ORIG" "$CUR"; do
  [ -x "$b" ] || { echo "call_regress: not an executable: $b" >&2; exit 2; }
done
[ -f "$EXPECTED" ] || { echo "call_regress: no expected-diff file at $EXPECTED" >&2; exit 2; }
case "$LOCI" in *myom2*) echo "call_regress: myom2 is excluded until it is rebuilt" >&2; exit 2;; esac
for set in $SETS; do
  case "$set" in base|refined|altref) ;; *) echo "call_regress: unknown input set $set" >&2; exit 2;; esac
done
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)/run.$$.$(date +%s)"; rm -rf "$OUT"
IN="$OUT/inputs"; RUNS="$OUT/runs"; ITEMS="$OUT/items"
mkdir -p "$IN" "$RUNS" "$ITEMS"
t0=$(date +%s)

# ref_of <gfa or gfa.gz>: GRCh38#0... if present, else the first P/W name matching grch38
ref_of() {
  local cat=cat; case "$1" in *.gz) cat="gunzip -c";; esac
  $cat "$1" | awk -F'\t' '$1=="P"{print $2} $1=="W"{print $2"#"$3"#"$4}' \
    | awk '/^GRCh38#0/{print; f=1; exit} tolower($0)~/grch38/ && g==""{g=$0}
           END{if(!f && g!="") print g}' | head -1
}

# first_of <gfa.gz>: the first P name, else the first W name (the altref set is built along it)
first_of() {
  gunzip -c "$1" | awk -F'\t' '$1=="P"{print $2; exit} $1=="W" && w==""{w=$2"#"$3"#"$4}
                              END{if(w!="") print w}' | head -1
}

case_args() {
  case "$1" in
    default)   echo "";;
    cn)        echo "--cn";;
    classify)  echo "--cn --classify-ins";;
    *) return 1;;
  esac
}
for cfg in $CONFIGS; do case_args "$cfg" >/dev/null || { echo "call_regress: unknown config $cfg" >&2; exit 2; }; done

# the stated reason of a failed run: its Error lines, with any [module hh:mm:ss] stamp removed
reason_of() { grep -E '^(Error|error|terminate|libc\+\+abi)' "$1" 2>/dev/null | sed -E 's/\[[a-z]+ [0-9:]+\] //' | head -3; }

# check_built <loc> <set> <built prefix> <reference prefix>: cmp freshly built inputs against a
# prepared copy, so a reused input directory is known to be what this script would build
check_built() {
  local x c
  for x in normalized.sorted.gfa bubbles.csv; do
    cmp -s "$3.$x" "$4.$x" && c=same || c=DIFFERENT
    echo "input check: $1 $2 .$x vs $4: $c" >> "$IN/inputs_check.txt.$1"
  done
}

# build_bp <loc> <path> <prefix>: bubble along <path>, then panphorte 0.95 on the same path
build_bp() {
  "$CUR" bubble -i "$(locus_gfa "$1")" -r "$2" -o "$3_b" -q >> "$IN/$1.prep.log" 2>&1 \
    && "$CUR" panphorte -i "$3_b.sorted.gfa" --bubble-prefix "$3_b" -o "$3_p" \
         --reference-path "$2" --min-similarity 0.95 -q >> "$IN/$1.prep.log" 2>&1
}

# prepare <loc>: writes $IN/<loc>.env with G_<set>/P_<set>/REF_<set> for every set in SETS. REF is
# GRCh38 for every set, read from that set's graph.
prepare() {
  local loc=$1 gz="$(locus_gfa "$1")" env="$IN/$1.env" p first
  : > "$env"
  if [[ " $SETS " = *" base "* || ( " $SETS " = *" refined "* && -z "${REFINED_DIR:-}" ) ]]; then
    if [ -n "${INPUT_DIR:-}" ]; then
      p="$INPUT_DIR/$loc$INPUT_SUFFIX"
      [ -s "$p.normalized.sorted.gfa" ] || { echo "# PREP_FAIL base: missing $p.normalized.sorted.gfa" >> "$env"; return; }
    else
      [ -s "$gz" ] || { echo "# PREP_FAIL base: missing $gz" >> "$env"; return; }
      p="$IN/${loc}"
      build_bp "$loc" "$(ref_of "$gz")" "$p" \
        || { echo "# PREP_FAIL base: bubble/panphorte, see $IN/$loc.prep.log" >> "$env"; return; }
      p="$IN/${loc}_p"
      [ -n "${CHECK_INPUTS_DIR:-}" ] && check_built "$loc" base "$p" "$CHECK_INPUTS_DIR/$loc$CHECK_INPUTS_SUFFIX"
    fi
    printf 'G_base=%q\nP_base=%q\nREF_base=%q\n' "$p.normalized.sorted.gfa" "$p" \
      "$(ref_of "$p.normalized.sorted.gfa")" >> "$env"
  fi
  if [[ " $SETS " = *" refined "* ]]; then
    if [ -n "${REFINED_DIR:-}" ]; then
      p="$REFINED_DIR/$loc$REFINED_SUFFIX"
    else
      # refine keeps the base set's reference
      ( . "$env"; "$CUR" refine -i "$G_base" --bubble-prefix "$P_base" -r "$REF_base" -o "$IN/${loc}_r" -q ) \
        >> "$IN/$loc.prep.log" 2>&1 || { echo "# PREP_FAIL refined: refine, see $IN/$loc.prep.log" >> "$env"; return; }
      p="$IN/${loc}_r"
    fi
    printf 'G_refined=%q\nP_refined=%q\nREF_refined=%q\n' "$p.normalized.sorted.gfa" "$p" \
      "$(ref_of "$p.normalized.sorted.gfa")" >> "$env"
  fi
  if [[ " $SETS " = *" altref "* ]]; then
    if [ -n "${ALTREF_DIR:-}" ]; then
      p="$ALTREF_DIR/$loc$ALTREF_SUFFIX"
      [ -s "$p.normalized.sorted.gfa" ] || { echo "# PREP_FAIL altref: missing $p.normalized.sorted.gfa" >> "$env"; return; }
    else
      [ -s "$gz" ] || { echo "# PREP_FAIL altref: missing $gz" >> "$env"; return; }
      first=$(first_of "$gz")
      build_bp "$loc" "$first" "$IN/${loc}_alt" \
        || { echo "# PREP_FAIL altref: bubble/panphorte on $first, see $IN/$loc.prep.log" >> "$env"; return; }
      p="$IN/${loc}_alt_p"
      echo "ALTREF_PATH=$(printf %q "$first")" >> "$env"
      [ -n "${CHECK_INPUTS_DIR:-}" ] && check_built "$loc" altref "$p" "$CHECK_INPUTS_DIR/$loc$CHECK_ALTREF_SUFFIX"
    fi
    printf 'G_altref=%q\nP_altref=%q\nREF_altref=%q\n' "$p.normalized.sorted.gfa" "$p" \
      "$(ref_of "$p.normalized.sorted.gfa")" >> "$env"
  fi
}

# item <set> <loc> <cfg> <file> <kind> <fields> <detail>: a difference call_regress finds itself.
# Its pin is the hash of the detail, so the run directory is taken out of it.
item() { printf '%s\t-\t-\t-\t%s\t1\t%s\n' "$5" "$6" "${7//$OUT/<out>}" >> "$ITEMS/$1__$2__$3__$4.items.tsv"; }

run_locus() {
  local loc=$1 env="$IN/$1.env" set cfg tag bin o rc_o rc_c why_o why_c d x key G P REF
  prepare "$loc"
  if grep -q PREP_FAIL "$env"; then
    set=$(grep PREP_FAIL "$env" | head -1 | sed -E 's/^# PREP_FAIL ([a-z]+):.*/\1/')
    item "$set" "$loc" - inputs PREP - "$(grep PREP_FAIL "$env" | head -1)"
    echo "inputs failed" > "$ITEMS/${set}__${loc}__-__inputs.summary"
  fi
  . "$env"
  for set in $SETS; do
    eval "G=\${G_$set:-}; P=\${P_$set:-}; REF=\${REF_$set:-}"
    if [ -z "$G" ]; then   # its own failure is reported above; a set after it was never prepared
      grep -q "PREP_FAIL $set:" "$env" || item "$set" "$loc" - inputs PREP - "not prepared: an earlier set failed"
      continue
    fi
    [ -s "$G" ] && [ -s "$P.bubbles.csv" ] || { item "$set" "$loc" - inputs PREP - "missing $G"; continue; }
    [ -n "$REF" ] || { item "$set" "$loc" - inputs PREP - "no GRCh38 path in $G"; continue; }
    key="${set}__${loc}__-__crossing"
    $PY "$TOOLS/crossings.py" "$G" "$P.bubbles.csv" --label "$set/$loc" --items "$ITEMS/$key.items.tsv" \
      > "$RUNS/$key.txt" 2>&1
    head -1 "$RUNS/$key.txt" | sed 's/^[^:]*: //' > "$ITEMS/$key.summary"
    for cfg in $CONFIGS; do
      d="$RUNS/$set/$loc/$cfg"; mkdir -p "$d"; key="${set}__${loc}__${cfg}"
      for tag in orig cur; do
        bin=$ORIG; [ "$tag" = cur ] && bin=$CUR
        # -b: the original binary names the long option --bubble-prefix-in
        # shellcheck disable=SC2046
        "$bin" call -i "$G" -b "$P" -r "$REF" -o "$d/$tag" --no-per-bubble-vcf -q \
          $(case_args "$cfg") > "$d/$tag.log" 2>&1
        echo $? > "$d/$tag.rc"
      done
      rc_o=$(cat "$d/orig.rc"); rc_c=$(cat "$d/cur.rc")
      why_o=$(reason_of "$d/orig.log"); why_c=$(reason_of "$d/cur.log")
      if [ "$rc_o" != "$rc_c" ]; then
        item "$set" "$loc" "$cfg" exit EXIT rc "orig exit $rc_o, cur exit $rc_c: ${why_o:-$why_c}"
      elif [ "$rc_o" != 0 ] && [ "$why_o" != "$why_c" ]; then
        item "$set" "$loc" "$cfg" exit EXIT reason "both exit $rc_o; orig: $why_o | cur: $why_c"
      fi
      if [ "$rc_o" != 0 ] || [ "$rc_c" != 0 ]; then
        echo "orig exit $rc_o, cur exit $rc_c: $(printf '%s' "${why_o:-$why_c}" | head -1 | cut -c1-140)" \
          > "$ITEMS/${key}__exit.summary"
        continue
      fi
      for x in region.vcf variant_nodes.tsv; do
        f=${x%%.*}; [ "$f" = variant_nodes ] && f=nodes
        [ -f "$d/orig.$x" ] || [ -f "$d/cur.$x" ] || continue
        if [ ! -f "$d/orig.$x" ] || [ ! -f "$d/cur.$x" ]; then
          item "$set" "$loc" "$cfg" "$f" FILE presence "$x written by only one binary"; continue
        fi
        tool=vcfdiff; [ "$f" = nodes ] && tool=nodesdiff
        $PY "$TOOLS/$tool.py" "$d/orig.$x" "$d/cur.$x" --label "$set/$loc/$cfg/$f" \
          --items "$ITEMS/${key}__$f.items.tsv" > "$d/$f.diff.txt" 2>&1
        [ $? -le 1 ] || item "$set" "$loc" "$cfg" "$f" TOOL "$tool" "$tool.py crashed, see $d/$f.diff.txt"
        head -1 "$d/$f.diff.txt" | sed 's/^[^:]*: //' > "$ITEMS/${key}__$f.summary"
      done
      # one copy-number record per bubble: a second DUP on a bubble counts its copies twice
      if [ -f "$d/cur.region.vcf" ]; then
        $PY "$TOOLS/cnrecords.py" "$d/cur.region.vcf" --orig "$d/orig.region.vcf" --label "$set/$loc/$cfg" \
          --items "$ITEMS/${key}__cnrec.items.tsv" > "$d/cnrec.txt" 2>&1
        [ $? -le 1 ] || item "$set" "$loc" "$cfg" cnrec TOOL cnrecords "cnrecords.py crashed, see $d/cnrec.txt"
        head -1 "$d/cnrec.txt" | sed 's/^[^:]*: //' > "$ITEMS/${key}__cnrec.summary"
      fi
      [ "$DESCRIBE" = 1 ] || continue
      for tag in orig cur; do
        "$CUR" describe --only-variant --variant-vcf "$d/$tag.region.vcf" -o "$d/$tag.dsc" -q \
          > "$d/$tag.dsc.log" 2>&1
        echo $? > "$d/$tag.dsc.rc"
      done
      rc_o=$(cat "$d/orig.dsc.rc"); rc_c=$(cat "$d/cur.dsc.rc")
      if [ "$rc_o" != 0 ] || [ "$rc_c" != 0 ]; then
        # both failing is still an item: the dosages of this run were not compared at all
        if [ "$rc_o" = "$rc_c" ]; then
          item "$set" "$loc" "$cfg" bimbam EXIT describe \
            "describe failed on both sides (exit $rc_o), dosages not compared: $(reason_of "$d/cur.dsc.log")"
        else
          item "$set" "$loc" "$cfg" bimbam EXIT rc \
            "describe exit orig $rc_o, cur $rc_c: $(reason_of "$d/orig.dsc.log")$(reason_of "$d/cur.dsc.log")"
        fi
        echo "describe exit orig $rc_o, cur $rc_c" > "$ITEMS/${key}__bimbam.summary"
        continue
      fi
      $PY "$TOOLS/bimbamdiff.py" "$d/orig.dsc" "$d/cur.dsc" --label "$set/$loc/$cfg" \
        --items "$ITEMS/${key}__bimbam.items.tsv" > "$d/bimbam.diff.txt" 2>&1
      [ $? -le 1 ] || item "$set" "$loc" "$cfg" bimbam TOOL bimbamdiff "bimbamdiff.py crashed, see $d/bimbam.diff.txt"
      head -1 "$d/bimbam.diff.txt" | sed 's/^[^:]*: //' > "$ITEMS/${key}__bimbam.summary"
    done
  done
}

echo "call_regress: orig=$ORIG"
echo "              cur=$CUR"
echo "              loci=[$LOCI] configs=[$CONFIGS] sets=[$SETS] out=$OUT"
[ -n "${INPUT_DIR:-}" ] && echo "              base inputs reused from $INPUT_DIR (*$INPUT_SUFFIX)"
[ -n "${REFINED_DIR:-}" ] && echo "              refined inputs reused from $REFINED_DIR (*$REFINED_SUFFIX)"
[ -n "${ALTREF_DIR:-}" ] && echo "              altref inputs reused from $ALTREF_DIR (*$ALTREF_SUFFIX)"
if [ "$PARALLEL" = 1 ]; then
  for loc in $LOCI; do run_locus "$loc" & done
  wait
else
  for loc in $LOCI; do run_locus "$loc"; done
fi
for loc in $LOCI; do
  [ -f "$IN/$loc.env" ] || continue
  for set in $SETS; do
    printf '  reference %-8s %-9s %s\n' "$set" "$loc" "$(. "$IN/$loc.env" 2>/dev/null; eval "echo \"\${REF_$set:-?}\"")"
  done
  x=$(. "$IN/$loc.env" 2>/dev/null; echo "${ALTREF_PATH:-}")
  [ -n "$x" ] && printf '  %-19s%-9s %s\n' 'altref built along' "$loc" "$x"
done
for loc in $LOCI; do [ -f "$IN/inputs_check.txt.$loc" ] && sed 's/^/  /' "$IN/inputs_check.txt.$loc"; done
echo
$PY "$TOOLS/regress_report.py" "$ITEMS" "$EXPECTED"
rc=$?
echo "call_regress: $(( $(date +%s) - t0 )) s; per-run diffs under $RUNS"
exit $rc
