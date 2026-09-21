#!/usr/bin/env bash
# The C4 synthetic-mosaic experiment: can proposal and scoring recover a mosaic whose two haplotypes are
# in no panel row?  Round 14 ran this ad hoc and drew a conclusion the measurement could not support.
# This harness exists so the experiment is reproducible, and so the claim it licenses is bounded by what
# it actually measures.
#
#   BUILD=<build dir> tests/genotype_c4_mosaic_experiment.sh <out-dir>
#
# Stages, each skipped when its outputs already exist:
#   1. stage   -- build the five mosaic truths and their staged graphs
#   2. gate    -- the curated staged model must reproduce ORDINARY C4 exactly
#   3. reads   -- simulate 30x from each staged mosaic
#   4. call    -- genotype against the staged graph with the mosaics curated out, K = 0 and 8
#   5. trace   -- proposal-only diagnostic at state limits 1024 and 4096
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${BUILD:-$REPO/build}"
OUT="${1:?usage: BUILD=<build> tests/genotype_c4_mosaic_experiment.sh <out-dir>}"
PROBE="$BUILD/genotype_probe"
STAGES="${STAGES:-stage gate reads call eval concord trace}"

# ---- frozen settings -------------------------------------------------------------------------------
REFERENCE='HG00096#1#haplotype1-0000024:31848049-32074186'
DONOR_A=HG02666           # both haplotypes forward-oriented; see the orientation note below
DONOR_B=NA19240
SEED=20260916             # salted per plan, so no two plans draw the same reads
DEPTH=30
READ_LENGTH=100
INSERT_MEAN=350
INSERT_SD=40
SIGMAS=4
SUBSTITUTION_ERROR=0.002
OUTLIER_WEIGHT=0.001
# Overridable for a targeted diagnostic run; the committed defaults are what the frozen experiment uses.
# K=128 in particular is a development diagnostic for every_block, NOT a new default.
STATE_LIMITS="${STATE_LIMITS:-1024 4096}"
CALL_PREFIXES="${CALL_PREFIXES:-0,8}"
TRACE_CANDIDATES=256

# Ordinary C4, for the model-equivalence gate: what the model must still be after two synthetic mosaic
# paths are added to the graph and then curated out of the panel.  Obtained from a plain-C4 run at
# commit 9f80baa (tests/results-style bench, reference below, min_variant_bp 50).  If the graph, the
# bubble caller or the block decomposition changes, this constant changes with it -- and the gate
# failing is then the correct outcome, not a nuisance to be updated away.
PLAIN_C4_FINGERPRINT=d2b8c6c2e0b4079d
PLAIN_C4_DIPLOTYPES=8515
PLAIN_C4_VARIANT_BLOCKS=5
PLAIN_C4_CONTEXT_BLOCKS=6

# C4's chain is 11 blocks: ctx:0 sv:1 ctx:2 sv:3 ctx:4 sv:5 ctx:6 sv:7 ctx:8 sv:9 ctx:10.
# A plan names the block indices at which the source donor flips, or "every".
#
# bridgeable/unbridgeable refer to whether a read pair can span the switch: sv:3 (166-293 bp) and
# sv:5 (0-117 bp) sit within a 350 bp insert of their neighbours, while sv:7 is 26-85 kb and cannot be
# crossed by any pair at any depth.
PLANS=(
  "one_switch:5"
  "bridgeable:3,5"
  "unbridgeable:8"
  "three_switches:3,6,9"
  "every_block:every"
)

# Orientation: C4 is 71 forward to 60 reverse templates, and mosaic-truth refuses to splice across
# orientations rather than silently joining walks that traverse the locus in opposite directions.  Of the
# eight frozen development donors, HG00733, HG03371 and NA20509 are reverse on both haplotypes and cannot
# be used until canonical-orientation splicing exists.  HG02666, HG02818, HG04036, NA19036 and NA19240 can.

mkdir -p "$OUT"
C4="$OUT/c4.gfa"
[ -f "$C4" ] || gzip -dc "$REPO/tests/real_data/c4.gfa.gz" > "$C4"
CURATION="$OUT/mosaic_exclude.tsv"
cat > "$CURATION" <<'EOF'
# The two synthetic mosaics are removed from the panel: the experiment requires that neither truth
# haplotype is reachable without recombining real templates.  They remain in the graph so reads can be
# simulated from them and so truth can be projected -- curation removes them from the PANEL, not the file.
MOSAIC#1#synthetic:0-0
MOSAIC#2#synthetic:0-0
EOF

has_stage() { case " $STAGES " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }

# One canonical signature over everything a call consumes: both read files, the library, the staged
# graph, the declared simulation parameters with the plan's salt, and the simulator binary itself.
# ALWAYS recomputed from the files on disk -- a stored signature is only ever compared against, never
# trusted as a description of what is there now.  That is what makes a changed reads_2.fq, a changed
# library.txt or a changed simulation parameter detectable even when the staged graph is untouched.
#
# Hashed by CONTENT (`shasum < file`), never `shasum file`, whose output embeds the path and would make
# a copied tree look tampered with.  A missing input yields MISSING:<what>, which cannot equal any
# stored signature, so absence refuses instead of silently comparing equal.
input_bundle_signature() {
  local d="$1" name="$2" f
  for f in "$d/staged.gfa" "$d/source/reads_1.fq" "$d/source/reads_2.fq" "$d/source/library.txt"; do
    [ -f "$f" ] || { printf 'MISSING:%s\n' "$(basename "$f")"; return 0; }
  done
  [ -f "$PROBE" ] || { printf 'MISSING:probe\n'; return 0; }
  {
    shasum -a 256 < "$d/staged.gfa"
    shasum -a 256 < "$d/source/reads_1.fq"
    shasum -a 256 < "$d/source/reads_2.fq"
    shasum -a 256 < "$d/source/library.txt"
    shasum -a 256 < "$PROBE"
    printf 'salt=%s seed=%s depth=%s len=%s ins=%s sd=%s sig=%s err=%s out=%s\n' \
      "$name" "$SEED" "$DEPTH" "$READ_LENGTH" "$INSERT_MEAN" "$INSERT_SD" "$SIGMAS" \
      "$SUBSTITUTION_ERROR" "$OUTLIER_WEIGHT"
  } | shasum -a 256 | cut -d' ' -f1
}

# ---- provenance ------------------------------------------------------------------------------------
{
  printf 'key\tvalue\n'
  printf 'git_commit\t%s\n' "$(git -C "$REPO" rev-parse HEAD)"
  printf 'git_tree_dirty\t%s\n' "$(git -C "$REPO" status --porcelain | grep -q . && echo true || echo false)"
  printf 'probe_sha256\t%s\n' "$(shasum -a 256 "$PROBE" | cut -d' ' -f1)"
  printf 'script_sha256\t%s\n' "$(shasum -a 256 "${BASH_SOURCE[0]}" | cut -d' ' -f1)"
  printf 'c4_gfa_sha256\t%s\n' "$(shasum -a 256 "$REPO/tests/real_data/c4.gfa.gz" | cut -d' ' -f1)"
  printf 'curation_sha256\t%s\n' "$(shasum -a 256 "$CURATION" | cut -d' ' -f1)"
  printf 'reference\t%s\n' "$REFERENCE"
  printf 'donor_a\t%s\ndonor_b\t%s\n' "$DONOR_A" "$DONOR_B"
  printf 'seed\t%s\ndepth\t%s\nread_length\t%s\n' "$SEED" "$DEPTH" "$READ_LENGTH"
  printf 'insert_mean_bp\t%s\ninsert_sd_bp\t%s\nsupport_sigmas\t%s\n' "$INSERT_MEAN" "$INSERT_SD" "$SIGMAS"
  printf 'substitution_error\t%s\noutlier_weight\t%s\n' "$SUBSTITUTION_ERROR" "$OUTLIER_WEIGHT"
  printf 'state_limits\t%s\ncall_prefixes\t%s\ntrace_candidates\t%s\n' "$STATE_LIMITS" "$CALL_PREFIXES" "$TRACE_CANDIDATES"
  for spec in "${PLANS[@]}"; do printf 'plan\t%s\n' "$spec"; done
} > "$OUT/manifest.tsv"
echo "manifest: $OUT/manifest.tsv"

# ---- 1. stage --------------------------------------------------------------------------------------
if has_stage stage; then
  echo "=== staging mosaic truths ==="
  for spec in "${PLANS[@]}"; do
    name="${spec%%:*}"; switches="${spec##*:}"
    d="$OUT/$name"; mkdir -p "$d"
    "$PROBE" mosaic-truth --gfa "$C4" --reference "$REFERENCE" \
      --donor-a "$DONOR_A" --donor-b "$DONOR_B" --switches "$switches" \
      --out-gfa "$d/staged.gfa" --out "$d/truth.tsv" > "$d/stage.log" 2>&1
    # The signature every later stage binds to.  Reads, calls and evaluations are only meaningful
    # against the graph that produced the truth they are scored on.
    shasum -a 256 "$d/staged.gfa" | cut -d' ' -f1 > "$d/staged.gfa.sha"
    printf '  %-16s %s\n' "$name" "$(tail -1 "$d/stage.log")"
  done
fi

# ---- 2. gate: the curated staged model must BE ordinary C4 -----------------------------------------
# Adding two P lines renumbers nodes wholesale (2,694 of 4,740 on C4), so the staged graph must be shown
# to leave the MODEL untouched once the mosaics are curated out -- not assumed to.
gate_models() {
  echo "=== model-equivalence gate ==="
  local bad=0 checked=0
  for spec in "${PLANS[@]}"; do
    local name="${spec%%:*}"
    local s="$OUT/$name/call.K0.mosaic_status.tsv"
    [ -f "$s" ] || { printf '  %-16s NO CALL: cannot gate a model that was never built\n' "$name"; bad=$((bad+1)); continue; }
    local fp cur dip vb cb
    fp="$(awk -F'\t' '$1=="model_fingerprint"{print $2}' "$s")"
    cur="$(awk -F'\t' '$1=="templates_curated_out"{print $2}' "$s")"
    dip="$(awk -F'\t' '$1=="panel_diplotypes"{print $2}' "$s")"
    vb="$(awk -F'\t' '$1=="variant_blocks"{print $2}' "$s")"
    cb="$(awk -F'\t' '$1=="context_blocks"{print $2}' "$s")"
    if [ "$fp" != "$PLAIN_C4_FINGERPRINT" ] || [ "$dip" != "$PLAIN_C4_DIPLOTYPES" ] \
       || [ "$vb" != "$PLAIN_C4_VARIANT_BLOCKS" ] || [ "$cb" != "$PLAIN_C4_CONTEXT_BLOCKS" ] \
       || [ "$cur" != "2" ]; then
      printf '  %-16s DIFFERS from ordinary C4: fp=%s (want %s) diplotypes=%s (want %s) blocks=%s/%s curated=%s\n' \
        "$name" "$fp" "$PLAIN_C4_FINGERPRINT" "$dip" "$PLAIN_C4_DIPLOTYPES" "$vb" "$cb" "$cur"
      bad=$((bad+1))
    else
      printf '  %-16s identical to ordinary C4: fp=%s diplotypes=%s blocks=%s/%s curated=%s\n' \
        "$name" "$fp" "$dip" "$vb" "$cb" "$cur"
      checked=$((checked+1))
    fi
  done
  if [ "$bad" -ne 0 ]; then
    echo "model-equivalence gate FAILED: staging the mosaics changed the model, so any result is void"
    exit 1
  fi
  echo "  $checked plans: adding two mosaic P lines and curating them out leaves the model untouched"
}

# ---- 3. reads --------------------------------------------------------------------------------------
if has_stage reads; then
  echo "=== simulating reads ==="
  for spec in "${PLANS[@]}"; do
    name="${spec%%:*}"; d="$OUT/$name"
    # Reads may be reused ONLY when they were drawn from the staged graph that is there now.  Reusing on
    # mere existence is what produced round 14's false conclusion: the graph was regenerated after a
    # fix, the FASTQs were not, and every later call scored reads from one pair of haplotypes against a
    # truth projected from another.  Silence was the whole problem, so a mismatch is loud and re-draws.
    want="$(input_bundle_signature "$d" "$name")"
    have="$(cat "$d/source/bundle.sha" 2>/dev/null || true)"
    if [ -f "$d/source/reads_1.fq" ] && [ -n "$have" ] && [ "$have" = "$want" ]; then
      printf '  %-16s reused (input bundle %s)\n' "$name" "$(echo "$want" | cut -c1-12)"
      continue
    fi
    if [ -f "$d/source/reads_1.fq" ]; then
      printf '  %-16s STALE INPUTS: bundle was %s, is now %s -- re-simulating\n' \
        "$name" "$(echo "${have:-MISSING}" | cut -c1-12)" "$(echo "$want" | cut -c1-12)"
    fi
    mkdir -p "$d/source"
    "$PROBE" simulate --gfa "$d/staged.gfa" --reference "$REFERENCE" --donor MOSAIC --depth "$DEPTH" \
      --read-length "$READ_LENGTH" --insert-mean "$INSERT_MEAN" --insert-sd "$INSERT_SD" \
      --sigmas "$SIGMAS" --error "$SUBSTITUTION_ERROR" --outlier "$OUTLIER_WEIGHT" \
      --seed "$SEED" --seed-salt "$name" --out "$d/source" > "$d/sim.log" 2>&1
    # graph.sha is informational provenance only; bundle.sha is what every later stage checks.
    cp "$d/staged.gfa.sha" "$d/source/graph.sha"
    input_bundle_signature "$d" "$name" > "$d/source/bundle.sha"
    printf '  %-16s %s fragments\n' "$name" "$(awk -F'\t' '$1=="fragments_written"{print $2}' "$d/source/simulation.tsv")"
  done
fi

# ---- 4. call ---------------------------------------------------------------------------------------
# The probe exits 0 even when a run refuses, so the STATUS FILE is the signal, never the exit code.
if has_stage call; then
  echo "=== calling (staged graph, mosaics curated out) ==="
  for spec in "${PLANS[@]}"; do
    name="${spec%%:*}"; d="$OUT/$name"
    # Written under a temporary prefix and published only once a new output set exists.  A failed rerun
    # must never leave last run's outputs beside a signature claiming they match the current inputs:
    # that combination is indistinguishable from success and is precisely what the signatures exist to
    # prevent.  The signature is written after publication, never before, and never on the failure path.
    rm -f "$d"/pending.*
    "$PROBE" genotype --gfa "$d/staged.gfa" --reference "$REFERENCE" \
      --reads1 "$d/source/reads_1.fq" --reads2 "$d/source/reads_2.fq" --library "$d/source/library.txt" \
      --out-prefix "$d/pending" --state-limit 1024 --prefixes "$CALL_PREFIXES" --phase-audit 0 --threads 2 \
      --panel-curation "$CURATION" > "$d/pending.log" 2>&1 || true
    if [ ! -f "$d/pending.K0.mosaic_status.tsv" ]; then
      printf '  %-16s CALL PRODUCED NO PER-K STATUS: previous outputs and signature left untouched: %s\n' \
        "$name" "$(awk -F'\t' '$1=="note"{print substr($2,1,90)}' "$d/pending.mosaic_status.tsv" 2>/dev/null)"
      rm -f "$d"/pending.*
      continue
    fi
    rm -f "$d"/call.*
    for f in "$d"/pending.*; do mv "$f" "$d/call.${f##*/pending.}"; done
    input_bundle_signature "$d" "$name" > "$d/call.bundle.sha"
    printf '  %-16s K0=%s K8=%s\n' "$name" \
      "$(awk -F'\t' '$1=="status"{print $2}' "$d/call.K0.mosaic_status.tsv")" \
      "$(awk -F'\t' '$1=="status"{print $2}' "$d/call.K8.mosaic_status.tsv" 2>/dev/null)"
  done
fi

# The gate runs AFTER the call stage: it compares the model each call actually built, so it has nothing
# to check until those statuses exist.
if has_stage gate; then gate_models; fi

# ---- 4b. eval: accuracy against the truth the reads were actually drawn from -----------------------
# Every stage below re-checks the read/graph binding.  A call is only interpretable against the graph its
# reads came from, and that is the invariant round 14 broke.
check_binding() {
  local d="$1" name="$2" now stored called
  now="$(input_bundle_signature "$d" "$name")"
  case "$now" in
    MISSING:*) printf '  %-16s INPUT MISSING (%s) -- rerun the reads stage\n' "$name" "$now"; return 1 ;;
  esac
  stored="$(cat "$d/source/bundle.sha" 2>/dev/null || true)"
  # A missing signature is a refusal, not compatibility: unsigned inputs are exactly the state whose
  # interpretation round 14 got wrong, and defaulting to "probably fine" is how that happened.
  if [ -z "$stored" ]; then
    printf '  %-16s NO INPUT SIGNATURE -- rerun the reads stage\n' "$name"; return 1
  fi
  if [ "$now" != "$stored" ]; then
    printf '  %-16s INPUTS CHANGED since simulation (%s -> %s) -- rerun the reads stage\n' \
      "$name" "$(echo "$stored" | cut -c1-12)" "$(echo "$now" | cut -c1-12)"; return 1
  fi
  called="$(cat "$d/call.bundle.sha" 2>/dev/null || true)"
  if [ -z "$called" ]; then
    printf '  %-16s NO CALL SIGNATURE -- rerun the call stage\n' "$name"; return 1
  fi
  if [ "$called" != "$now" ]; then
    printf '  %-16s CALLS MADE ON OTHER INPUTS (%s -> %s) -- rerun the call stage\n' \
      "$name" "$(echo "$called" | cut -c1-12)" "$(echo "$now" | cut -c1-12)"; return 1
  fi
  return 0
}

if has_stage eval; then
  echo "=== accuracy against truth ==="
  # Truth is evaluation only and never enters consensus: blocks_exact below is scored AFTER the call,
  # against a truth the caller never saw.
  printf '%-15s %-2s %-18s %-4s %-7s %-8s %-7s %-12s %s\n' \
    plan K status tied called nocalled exact phase uniq_hap
  for spec in "${PLANS[@]}"; do
    name="${spec%%:*}"; d="$OUT/$name"
    check_binding "$d" "$name" || continue
    for k in ${CALL_PREFIXES//,/ }; do
      s="$d/call.K$k.mosaic_status.tsv"
      [ -f "$s" ] || continue
      st="$(awk -F'\t' '$1=="status"{print $2}' "$s")"
      ws="$(awk -F'\t' '$1=="winner_set_size"{print $2}' "$s")"
      cc="$(awk -F'\t' '$1=="consensus_blocks_called"{print $2}' "$s")"
      cn="$(awk -F'\t' '$1=="consensus_blocks_no_called"{print $2}' "$s")"
      ph="$(awk -F'\t' '$1=="phase_support"{print $2}' "$s")"
      fa=no; [ -f "$d/call.K$k.mosaic_haplotypes.fa" ] && fa=yes
      bx=NA
      if "$PROBE" evaluate --run-prefix "$d/call.K$k" --out "$d/eval.K$k" \
           --scheme "c4-mosaic-$name-K$k" --donor MOSAIC --reference "$REFERENCE" \
           --model-gfa "$d/staged.gfa" --truth-gfa "$d/staged.gfa" \
           --panel-curation "$CURATION" > "$d/eval.K$k.log" 2>&1; then
        bx="$(awk -F'\t' '$1=="blocks_exact"{print $2}' "$d/eval.K$k.summary.tsv")"
      fi
      printf '%-15s %-2s %-18s %-4s %-7s %-8s %-7s %-12s %s\n' \
        "$name" "$k" "$st" "${ws:-NA}" "${cc:-NA}" "${cn:-NA}" "${bx:-NA}" "${ph:-NA}" "$fa"
    done
  done
fi

# ---- 4c. concord: do the tied winners agree block by block? ----------------------------------------
# A refusal declines a whole-locus winner.  Whether the tied candidates nonetheless agree on each block
# genotype is the measurement that decides whether a block-consensus call would be sound -- so the
# harness produces it rather than leaving it to a one-off script.
if has_stage concord; then
  echo "=== winner-set block concordance ==="
  for spec in "${PLANS[@]}"; do
    name="${spec%%:*}"; d="$OUT/$name"
    check_binding "$d" "$name" || continue
    for k in ${CALL_PREFIXES//,/ }; do
      python3 - "$d" "$name" "$k" <<'PY'
import csv, os, sys
d, name, k = sys.argv[1], sys.argv[2], sys.argv[3]
sc, cd = f"{d}/call.K{k}.mosaic_scores.tsv", f"{d}/call.K{k}.mosaic_candidates.tsv"
if not (os.path.exists(sc) and os.path.exists(cd)):
    print(f"  {name:<16} K{k:<3} no scores/candidates written"); raise SystemExit
rows = lambda p: list(csv.DictReader(open(p), delimiter='\t'))
win = {r["representative_beam_candidate"] for r in rows(sc) if r.get("in_winner_set") == "true"}
by_id = {c["candidate_id"]: c for c in rows(cd)}
sel = [by_id[w] for w in win if w in by_id]
if not sel:
    print(f"  {name:<16} K{k:<3} winner set {len(win)}, none resolvable"); raise SystemExit
al = lambda s: [int(x) for x in s.split(",") if x != ""]
nb = len(al(sel[0]["hap1_alleles"]))
agree = [b for b in range(nb)
         if len({tuple(sorted((al(c["hap1_alleles"])[b], al(c["hap2_alleles"])[b]))) for c in sel}) == 1]
note = "trivial: one winner" if len(sel) == 1 else ("unanimous" if len(agree) == nb else "disagree somewhere")
print(f"  {name:<16} K{k:<3} winners={len(sel):<3} blocks_agreeing={len(agree)}/{nb}  {note}")
PY
    done
  done
fi

# ---- 5. trace: the proposal-only diagnostic --------------------------------------------------------
# Searches the COMPLETE final beam, not the ranked prefix written to mosaic_candidates.tsv.  An absent
# verdict from the prefix cannot distinguish ranked-beyond-K, pruned-at-a-layer and never-generated;
# truth_in_full_beam and first_block_truth_not_retained can.
if has_stage trace; then
  echo "=== proposal-only trace ==="
  printf '%-16s %-7s %-14s %-12s %-8s %-10s %-10s %s\n' \
    plan limit truth_in_beam full_beam_rank found agrees_id agrees_score first_loss
  for spec in "${PLANS[@]}"; do
    name="${spec%%:*}"; d="$OUT/$name"
    [ -f "$d/call.K8.mosaic_status.tsv" ] || continue
    for limit in $STATE_LIMITS; do
      "$PROBE" trace --gfa "$d/staged.gfa" --reference "$REFERENCE" --donor MOSAIC \
        --truth-gfa "$d/staged.gfa" --panel-curation "$CURATION" --grouped 1 \
        --run "$d/call.K8" --reads1 "$d/source/reads_1.fq" --reads2 "$d/source/reads_2.fq" \
        --library "$d/source/library.txt" --state-limit "$limit" --candidates "$TRACE_CANDIDATES" \
        --out "$d/trace.$limit" > "$d/trace.$limit.log" 2>&1
      f="$d/trace.$limit.trace.summary.tsv"
      # Crowding, over the layers: the largest share of the width taken by one diplotype prefix, and
      # the fewest distinct prefixes any layer held.  Columns are looked up by NAME so the report does
      # not silently follow a column that moved.
      t="$d/trace.$limit.trace.tsv"
      printf '%-16s %-7s %-14s %-12s %-8s %-10s %-10s %s\n' "$name" "$limit" \
        "$(awk -F'\t' '$1=="truth_in_full_beam"{print $2}' "$f")" \
        "$(awk -F'\t' '$1=="truth_full_beam_rank"{print $2}' "$f")" \
        "$(awk -F'\t' '$1=="unique_candidates_found"{print $2}' "$f")" \
        "$(awk -F'\t' '$1=="proposal_matches_run"{print $2}' "$f")" \
        "$(awk -F'\t' '$1=="proposal_matches_run_score"{print $2}' "$f")" \
        "$(awk -F'\t' '$1=="first_block_truth_not_retained"{print $2}' "$f")"
    done
  done
fi
echo "done: $OUT"
