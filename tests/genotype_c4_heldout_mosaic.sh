#!/usr/bin/env bash
# The HELD-OUT C4 synthetic-mosaic validation.
#
#   tests/genotype_c4_heldout_mosaic.sh --check-manifest          # the manifest gate alone (ctest)
#   BUILD=<build dir> tests/genotype_c4_heldout_mosaic.sh <out-dir>
#
# THE DEFAULT STAGES DO NOT SIMULATE, CALL OR EVALUATE ANYTHING. They are "manifest stage verify":
# check that the committed plan table is what the declared draw rule produces, build the ten mosaic
# truths, and prove each one is admissible against the frozen ordinary-C4 model. Simulation waits for
# review, and is reached only by asking for it explicitly:
#
#   STAGES="reads nested call complete eval concord trace" BUILD=<build> \
#     tests/genotype_c4_heldout_mosaic.sh <out-dir>
#
# That split is the point of this harness. A protocol that can be half-run by accident is not frozen,
# and the whole value of held-out data is that nothing touched it before the protocol was fixed.
#
# Everything the run obeys comes from tests/genotype_c4_heldout_mosaic_manifest.tsv, which was committed
# before any of this ran. Nothing here re-selects a caller setting: width, K, the diagnostic comparison
# width and the phase audit are read from the manifest and asserted to be the frozen production values.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
HARNESS="${BASH_SOURCE[0]}"
MANIFEST="$REPO/tests/genotype_c4_heldout_mosaic_manifest.tsv"
ORIENTATION="$REPO/tests/genotype_c4_template_orientation.tsv"
PANEL_ALLELES="$REPO/tests/genotype_c4_panel_alleles.tsv"
DONORS="$REPO/tests/genotype_development_donors.tsv"
DRAW="$REPO/scripts/draw_heldout_mosaic_plans.py"
GRAPH_GZ="$REPO/tests/real_data/c4.gfa.gz"

setting() { awk -F'\t' -v key="$1" '$1==key{print $2; exit}' "$MANIFEST"; }
plan_rows() { awk -F'\t' '$1=="plan" && $2!="name"' "$MANIFEST"; }

# Every loop below reads plan rows as pipe-separated, NOT tab-separated, and it matters. Tab is IFS
# WHITESPACE, so `IFS=$'\t' read` collapses consecutive tabs into one delimiter and drops empty fields.
# A one_homologue_once row has an empty switches_hap2, so reading it that way silently shifted
# require_off_panel into the switch spec and handed the probe "one" as a block index, which surfaced as
# "stoul: no conversion" from deep inside mosaic-truth. A non-whitespace delimiter preserves empties.
# Donor names, switch specs and class names never contain '|'.

# ---- the manifest gate ------------------------------------------------------------------------------
# Re-derives the plan table from the rule and refuses on any difference. A manifest that has drifted from
# its own stated rule is no longer a rule: it is a set of choices wearing one, and it would be impossible
# to tell afterwards which rows were drawn and which were edited.
check_manifest() {
  echo "=== manifest gate ==="
  local bad=0

  local derived committed
  derived="$(python3 "$DRAW" "$ORIENTATION" "$DONORS" "$PANEL_ALLELES" | awk -F'\t' '$1=="plan" && $2!="name"')"
  committed="$(plan_rows)"
  if [ "$derived" != "$committed" ]; then
    echo "  DRAW MISMATCH: the committed plan table is not what the rule produces"
    diff <(printf '%s\n' "$committed") <(printf '%s\n' "$derived") | sed 's/^/    /' || true
    bad=$((bad + 1))
  else
    echo "  plan table re-derives from the rule: $(printf '%s\n' "$committed" | grep -c .) plans"
  fi

  # No development donor may appear. The draw excludes them by role, so a hit here means either the
  # role table or the draw changed, and both invalidate the held-out claim.
  local dev_hits=0 donor
  for donor in $(awk -F'\t' '$2=="development"{print $1}' "$DONORS"); do
    if plan_rows | cut -f4,5 | grep -qw "$donor"; then
      echo "  DEVELOPMENT DONOR $donor appears in the held-out draw"
      dev_hits=$((dev_hits + 1))
    fi
  done
  [ "$dev_hits" -eq 0 ] && echo "  no development donor appears in any plan" || bad=$((bad + 1))

  # No plan may repeat a development switch plan, in either homologue's spec.
  local repeat=0 spec
  while IFS= read -r spec; do
    [ -n "$spec" ] || continue
    case "$(printf '%s' "$spec" | tr ',' ' ' | tr ' ' '\n' | sort -n | tr '\n' ',' )" in
      "5,"|"3,5,"|"8,"|"3,6,9,") echo "  DEVELOPMENT SWITCH PLAN reused: $spec"; repeat=$((repeat + 1)) ;;
      "every,") echo "  DEVELOPMENT 'every' plan reused"; repeat=$((repeat + 1)) ;;
    esac
  done < <(plan_rows | cut -f7,8 | tr '\t' '\n')
  [ "$repeat" -eq 0 ] && echo "  no plan repeats a development switch plan" || bad=$((bad + 1))

  # Donors are used at most once, so no two plans share a haplotype source.
  local dup
  dup="$(plan_rows | cut -f4,5 | tr '\t' '\n' | sort | uniq -d | tr '\n' ' ')"
  if [ -n "$dup" ]; then echo "  DONOR REUSED across plans: $dup"; bad=$((bad + 1));
  else echo "  every plan draws two donors used by no other plan"; fi

  # Classes: five, two plans each.
  local classes
  classes="$(plan_rows | cut -f3 | sort | uniq -c | awk '{print $2"="$1}' | tr '\n' ' ')"
  if [ "$(plan_rows | cut -f3 | sort -u | wc -l | tr -d ' ')" != "5" ] ||
     [ -n "$(plan_rows | cut -f3 | sort | uniq -c | awk '$1!=2')" ]; then
    echo "  CLASS SHAPE WRONG: $classes"; bad=$((bad + 1))
  else
    echo "  five classes, two plans each: $classes"
  fi

  # The frozen caller settings, asserted rather than trusted: this protocol consumes the production
  # settings and may not quietly run at some other width or K.
  local width k diag audit
  width="$(setting proposal_state_limit)"; k="$(setting max_mosaic_candidates)"
  diag="$(setting diagnostic_comparison_width)"; audit="$(setting phase_audit_candidates)"
  if [ "$width" = "256" ] && [ "$k" = "128" ] && [ "$diag" = "1024" ] && [ "$audit" = "0" ]; then
    echo "  frozen settings: width $width, K $k, diagnostic comparison width $diag, phase audit $audit"
  else
    echo "  FROZEN SETTINGS CHANGED: width=$width K=$k diagnostic=$diag audit=$audit"; bad=$((bad + 1))
  fi

  # Four coverages, nested from one simulation at the source coverage, and the source must be among them.
  local covs ncov
  covs="$(setting coverages)"; ncov="$(printf '%s' "$covs" | tr ',' '\n' | grep -c .)"
  if [ "$ncov" = "4" ] && printf '%s' "$covs" | tr ',' '\n' | grep -qx "$(setting source_coverage)"; then
    echo "  four coverages, source included: $covs"
  else
    echo "  COVERAGE SHAPE WRONG: $covs (source $(setting source_coverage))"; bad=$((bad + 1))
  fi

  # Completeness is declared here, not inferred later from whatever succeeded.
  local plans_n want_runs want_status
  plans_n="$(plan_rows | grep -c .)"
  want_runs="$(setting expected_shared_runs)"; want_status="$(setting expected_k_arm_status_records)"
  if [ "$want_runs" = "$((plans_n * ncov))" ] &&
     [ "$want_status" = "$((plans_n * ncov * $(printf '%s' "$(setting call_prefixes)" | tr ',' '\n' | grep -c .)))" ]; then
    echo "  completeness declared: $want_runs shared runs, $want_status K-arm status records"
  else
    echo "  COMPLETENESS DECLARATION INCONSISTENT: $want_runs runs / $want_status records for"
    echo "    $plans_n plans x $ncov coverages x $(printf '%s' "$(setting call_prefixes)" | tr ',' '\n' | grep -c .) arms"
    bad=$((bad + 1))
  fi

  # Every metric and diagnostic is declared before any result exists, each naming its source key.
  local metrics
  metrics="$(awk -F'\t' '$1=="metric"' "$MANIFEST" | wc -l | tr -d ' ')"
  if [ "$metrics" -ge 30 ]; then echo "  $metrics metrics and diagnostics frozen, each naming its source key"
  else echo "  METRICS MISSING: only $metrics declared"; bad=$((bad + 1)); fi

  if [ "$bad" -ne 0 ]; then echo "manifest gate FAILED"; exit 1; fi
  echo "manifest gate passed"
}

if [ "${1:-}" = "--check-manifest" ]; then check_manifest; exit 0; fi

OUT="${1:?usage: BUILD=<build> tests/genotype_c4_heldout_mosaic.sh <out-dir> | --check-manifest}"
BUILD="${BUILD:-$REPO/build}"
PROBE="$BUILD/genotype_probe"
STAGES="${STAGES:-manifest stage verify}"
has_stage() { case " $STAGES " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }

REFERENCE="$(setting reference)"
SEED="$(setting seed)"
SOURCE_COVERAGE="$(setting source_coverage)"
COVERAGES="$(setting coverages | tr ',' ' ')"
READ_LENGTH="$(setting read_length)"
INSERT_MEAN="$(setting insert_mean_bp)"
INSERT_SD="$(setting insert_sd_bp)"
SIGMAS="$(setting support_sigmas)"
SUBSTITUTION_ERROR="$(setting substitution_error)"
OUTLIER_WEIGHT="$(setting outlier_weight)"
STATE_LIMIT="$(setting proposal_state_limit)"
CALL_PREFIXES="$(setting call_prefixes)"
DIAGNOSTIC_WIDTH="$(setting diagnostic_comparison_width)"
PHASE_AUDIT="$(setting phase_audit_candidates)"
TRACE_CANDIDATES="$(setting trace_candidates)"
MIN_VARIANT_BP="$(setting min_variant_bp)"
EXPECTED_RUNS="$(setting expected_shared_runs)"
EXPECTED_STATUS="$(setting expected_k_arm_status_records)"

mkdir -p "$OUT"
C4="$OUT/c4.gfa"
[ -f "$C4" ] || gzip -dc "$GRAPH_GZ" > "$C4"
CURATION="$OUT/mosaic_exclude.tsv"
cat > "$CURATION" <<'EOF'
# The two synthetic mosaics are removed from the PANEL, not from the file: reads are simulated from them
# and truth is projected through them, but no candidate may reach them without recombining real panel
# templates. For the panel_control plans the two paths are copies of real panel rows, and curating them
# out leaves those rows in the panel -- which is exactly what makes a control a control.
MOSAIC#1#synthetic:0-0
MOSAIC#2#synthetic:0-0
EOF

# =====================================================================================================
# Signatures.
#
# The old provenance file was written fresh on every invocation, so the record of what a result came
# from was destroyed by the next run that touched the directory -- including a run with different
# inputs. These are content-addressed instead: the signature IS the hash of the bound content, it is
# written to a file named by that hash, and an existing file is never rewritten. A changed input
# produces a different id and a new file beside the old one, so the history of what was run cannot be
# overwritten, only added to.
# =====================================================================================================

content_sha() {
  local file="$1"
  [ -f "$file" ] || { printf 'MISSING:%s\n' "$(basename "$file")"; return 0; }
  shasum -a 256 < "$file" | cut -d' ' -f1
}

# The caller settings a result depends on. A run at another width or K is a different experiment, and
# must not be able to masquerade as this one.
caller_settings_blob() {
  printf 'cov=%s src=%s seed=%s len=%s ins=%s sd=%s sig=%s err=%s out=%s width=%s prefixes=%s audit=%s minvar=%s\n' \
    "$1" "$SOURCE_COVERAGE" "$SEED" "$READ_LENGTH" "$INSERT_MEAN" "$INSERT_SD" "$SIGMAS" \
    "$SUBSTITUTION_ERROR" "$OUTLIER_WEIGHT" "$STATE_LIMIT" "$CALL_PREFIXES" "$PHASE_AUDIT" "$MIN_VARIANT_BP"
}

# Everything a result at (plan, coverage) is conditional on, component by component. ALWAYS recomputed
# from the files on disk: a stored signature is only ever compared against, never trusted as a
# description of what is there now.
experiment_signature() {
  # Split deliberately. `local a="$1" d="$OUT/$a"` expands every assignment word BEFORE performing any
  # of the assignments, so $a is still unset there and `set -u` aborts the whole run. None of these
  # functions runs in the default stages, so the fault would have surfaced at the first simulated run.
  local name="$1" cov="$2"
  local d="$OUT/$name"
  printf 'manifest\t%s\n' "$(content_sha "$MANIFEST")"
  printf 'harness\t%s\n' "$(content_sha "$HARNESS")"
  printf 'draw_script\t%s\n' "$(content_sha "$DRAW")"
  printf 'orientation_table\t%s\n' "$(content_sha "$ORIENTATION")"
  printf 'panel_allele_table\t%s\n' "$(content_sha "$PANEL_ALLELES")"
  printf 'graph\t%s\n' "$(content_sha "$GRAPH_GZ")"
  printf 'curation\t%s\n' "$(content_sha "$CURATION")"
  printf 'plan_row\t%s\n' "$(plan_rows | awk -F'\t' -v n="$name" '$2==n' | shasum -a 256 | cut -d' ' -f1)"
  printf 'staged_graph\t%s\n' "$(content_sha "$d/staged.gfa")"
  printf 'staged_truth\t%s\n' "$(content_sha "$d/truth.tsv")"
  printf 'reads_1\t%s\n' "$(content_sha "$d/cov$cov/reads_1.fq")"
  printf 'reads_2\t%s\n' "$(content_sha "$d/cov$cov/reads_2.fq")"
  printf 'library\t%s\n' "$(content_sha "$d/cov$cov/library.txt")"
  printf 'binary\t%s\n' "$(content_sha "$PROBE")"
  printf 'caller_settings\t%s\n' "$(caller_settings_blob "$cov" | shasum -a 256 | cut -d' ' -f1)"
}

experiment_signature_id() { experiment_signature "$1" "$2" | shasum -a 256 | cut -d' ' -f1; }

# A signature naming a missing component can never equal a stored one, but saying so plainly beats
# leaving the caller to notice that two opaque hashes differ.
signature_incomplete() {
  # Deliberately NOT `experiment_signature ... | grep -q`. Under `set -o pipefail` grep -q exits at the
  # first match, the producer dies of SIGPIPE, and the pipeline reports 141 -- so a signature naming a
  # MISSING component came back "complete" and the call stage's input-missing refusal never fired.
  case "$(experiment_signature "$1" "$2")" in
    *MISSING:*) return 0 ;;
    *) return 1 ;;
  esac
}

record_experiment_signature() {
  local name="$1" cov="$2" components id file
  components="$(experiment_signature "$name" "$cov")"
  id="$(printf '%s\n' "$components" | shasum -a 256 | cut -d' ' -f1)"
  mkdir -p "$OUT/signatures"
  file="$OUT/signatures/$id.tsv"
  if [ ! -f "$file" ]; then
    {
      printf 'key\tvalue\n'
      printf 'signature_id\t%s\n' "$id"
      printf 'plan\t%s\ncoverage\t%s\n' "$name" "$cov"
      printf 'recorded_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
      printf 'git_commit\t%s\n' "$(git -C "$REPO" rev-parse HEAD 2>/dev/null || echo unknown)"
      printf 'git_tree_dirty\t%s\n' "$(git -C "$REPO" status --porcelain 2>/dev/null | grep -q . && echo true || echo false)"
      printf '%s\n' "$components"
    } > "$file.tmp"
    mv "$file.tmp" "$file"
    # Advisory only -- it does not stop root or a -f removal, but it does stop the casual overwrite
    # that destroyed the old provenance file.
    chmod a-w "$file" 2>/dev/null || true
    printf '%s\t%s\t%s\t%s\n' "$name" "$cov" "$id" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$OUT/signatures/index.tsv"
  fi
  printf '%s\n' "$id"
}

# The call signature binds the experiment signature to the outputs that call produced. eval, concord and
# trace all verify it before reading anything: a changed read, curation file, caller setting or binary
# moves the experiment signature, and a changed call output moves this one, and either refuses
# interpretation rather than quietly reinterpreting one run's outputs against another run's inputs.
call_signature() {
  local name="$1" cov="$2" k
  local d="$OUT/$name"
  {
    printf 'experiment\t%s\n' "$(experiment_signature_id "$name" "$cov")"
    for k in ${CALL_PREFIXES//,/ }; do
      printf 'status_K%s\t%s\n' "$k" "$(content_sha "$d/call.$cov.K$k.mosaic_status.tsv")"
      printf 'scores_K%s\t%s\n' "$k" "$(content_sha "$d/call.$cov.K$k.mosaic_scores.tsv")"
      printf 'candidates_K%s\t%s\n' "$k" "$(content_sha "$d/call.$cov.K$k.mosaic_candidates.tsv")"
      printf 'blocks_K%s\t%s\n' "$k" "$(content_sha "$d/call.$cov.K$k.block_genotypes.tsv")"
    done
  } | shasum -a 256 | cut -d' ' -f1
}

require_call_signature() {
  local name="$1" cov="$2" stored now
  local d="$OUT/$name"
  stored="$(cat "$d/call.$cov.signature" 2>/dev/null || true)"
  if [ -z "$stored" ]; then
    printf '  %-22s cov%-3s NO CALL SIGNATURE -- rerun the call stage\n' "$name" "$cov"; return 1
  fi
  now="$(call_signature "$name" "$cov")"
  if [ "$now" != "$stored" ]; then
    printf '  %-22s cov%-3s REFUSED: inputs or outputs changed since calling (%s -> %s)\n' \
      "$name" "$cov" "$(printf '%s' "$stored" | cut -c1-12)" "$(printf '%s' "$now" | cut -c1-12)"
    return 1
  fi
  return 0
}

if has_stage manifest; then check_manifest; fi

# ---- stage: build the ten mosaic truths -------------------------------------------------------------
if has_stage stage; then
  echo "=== staging held-out mosaic truths ==="
  while IFS='|' read -r _ name class donor_a donor_b swap hap1 hap2 off_panel; do
    d="$OUT/$name"; mkdir -p "$d"
    # A column shift should say so here, not arrive as a parse error from inside the probe.
    case "$off_panel" in
      both|one|none) ;;
      *) echo "  $name: MALFORMED PLAN ROW (require_off_panel='$off_panel'): the manifest columns did not parse"
         exit 1 ;;
    esac
    "$PROBE" mosaic-truth --gfa "$C4" --reference "$REFERENCE" \
      --donor-a "$donor_a" --donor-b "$donor_b" \
      --switches-hap1 "$hap1" --switches-hap2 "$hap2" \
      --swap-hap-b "$swap" --require-off-panel "$off_panel" \
      --out-gfa "$d/staged.gfa" --out "$d/truth.tsv" > "$d/stage.log" 2>&1
    printf '  %-22s %-20s %s\n' "$name" "$class" "$(tail -1 "$d/stage.log")"
  done < <(plan_rows | tr '\t' '|')
fi

# ---- verify: is every staged model admissible? ------------------------------------------------------
# Everything here runs BEFORE any read exists, and refuses rather than reports.
if has_stage verify; then
  echo "=== verify (no reads, no calls) ==="
  bad=0

  # 1. Staging added the two mosaic paths and NOTHING else. Each staged graph is the reference-sorted
  # C4 plus two P lines, so stripping those lines must leave every plan holding byte-identical graphs.
  base_sha=""; base_plan=""
  while IFS='|' read -r _ name _ _ _ _ _ _ _; do
    sha="$(grep -v $'^P\tMOSAIC#' "$OUT/$name/staged.gfa" | shasum -a 256 | cut -d' ' -f1)"
    if [ -z "$base_sha" ]; then base_sha="$sha"; base_plan="$name"; fi
    if [ "$sha" != "$base_sha" ]; then
      printf '  %-22s STAGED GRAPH DIFFERS from %s outside the two mosaic paths\n' "$name" "$base_plan"
      bad=$((bad + 1))
    fi
    added="$(grep -c $'^P\tMOSAIC#' "$OUT/$name/staged.gfa" || true)"
    [ "$added" = "2" ] || { printf '  %-22s ADDED %s mosaic paths, expected 2\n' "$name" "$added"; bad=$((bad + 1)); }
  done < <(plan_rows | tr '\t' '|')
  [ "$bad" -eq 0 ] && echo "  all plans: staged graph is one shared sorted C4 plus exactly two mosaic paths"

  # 2. THE MODEL-EQUIVALENCE GATE, before simulation. Adding two P lines renumbers nodes wholesale, so
  # the curated staged model has to be SHOWN to still be ordinary C4. preflight now takes the curation
  # manifest, so this runs with no reads and no call: it builds each staged model exactly as a run
  # would, and compares fingerprint, panel diplotypes, block counts and the number of templates actually
  # curated out against the frozen constants.
  want_fp="$(setting plain_c4_model_fingerprint)"; want_dip="$(setting plain_c4_panel_diplotypes)"
  want_blocks="$(setting plain_c4_blocks)"; want_vb="$(setting plain_c4_variant_blocks)"
  want_cb="$(setting plain_c4_context_blocks)"; want_cur="$(setting plain_c4_templates_curated_out)"
  while IFS='|' read -r _ name _ _ _ _ _ _ _; do
    d="$OUT/$name"
    "$PROBE" preflight --gfa "$d/staged.gfa" --reference "$REFERENCE" \
      --panel-curation "$CURATION" --min-variant-bp "$MIN_VARIANT_BP" \
      --state-limit "$STATE_LIMIT" --max-mosaic-candidates "${CALL_PREFIXES##*,}" \
      --phase-audit-candidates "$PHASE_AUDIT" --out "$d/preflight.tsv" > "$d/preflight.log" 2>&1 || true
    if [ ! -f "$d/preflight.tsv" ]; then
      printf '  %-22s PREFLIGHT PRODUCED NOTHING: %s\n' "$name" "$(tail -1 "$d/preflight.log")"
      bad=$((bad + 1)); continue
    fi
    p() { awk -F'\t' -v key="$1" '$1==key{print $2; exit}' "$d/preflight.tsv"; }
    fp="$(p model_fingerprint)"; dip="$(p panel_diplotypes)"; blocks="$(p blocks)"
    vb="$(p variant_blocks)"; cb="$(p context_blocks)"; cur="$(p templates_curated_out)"
    if [ "$fp" = "$want_fp" ] && [ "$dip" = "$want_dip" ] && [ "$blocks" = "$want_blocks" ] &&
       [ "$vb" = "$want_vb" ] && [ "$cb" = "$want_cb" ] && [ "$cur" = "$want_cur" ]; then
      printf '  %-22s curated staged model IS ordinary C4: fp=%s diplotypes=%s blocks=%s/%s curated=%s\n' \
        "$name" "$fp" "$dip" "$vb" "$cb" "$cur"
    else
      printf '  %-22s MODEL DIFFERS: fp=%s(want %s) dip=%s(want %s) blocks=%s(want %s) %s/%s(want %s/%s) curated=%s(want %s)\n' \
        "$name" "$fp" "$want_fp" "$dip" "$want_dip" "$blocks" "$want_blocks" "$vb" "$cb" "$want_vb" "$want_cb" \
        "$cur" "$want_cur"
      bad=$((bad + 1))
    fi
  done < <(plan_rows | tr '\t' '|')

  # 3. The declared off-panel requirement, and how many homologues actually are panel rows. mosaic-truth
  # refuses the violating cases itself; this asserts the EXACT count each class implies, which is
  # stricter than what the probe enforces and catches a class whose draw degenerated.
  while IFS='|' read -r _ name class _ _ _ _ _ off_panel; do
    on="$(sed -n 's/.*, \([0-9]*\) homologue(s) on panel.*/\1/p' "$OUT/$name/stage.log")"
    case "$off_panel" in
      both) want=0 ;;
      one)  want=1 ;;
      none) want=2 ;;
    esac
    if [ "$on" != "$want" ]; then
      printf '  %-22s OFF-PANEL WRONG: %s homologue(s) on panel, class %s requires %s\n' \
        "$name" "${on:-NA}" "$class" "$want"
      bad=$((bad + 1))
    else
      printf '  %-22s %-20s %s of 2 homologues on panel, as %s requires\n' "$name" "$class" "$on" "$off_panel"
    fi
  done < <(plan_rows | tr '\t' '|')

  # 4. Representability: every block allele of both truth haplotypes must be carried by at least one
  # template of the curated panel. Carrier counts come from mosaic-truth, which computed them over the
  # plain C4 panel -- the two MOSAIC paths are not in that model, so these counts are already the
  # curated-panel counts. A truth the panel cannot spell would measure representation, not genotyping.
  while IFS='|' read -r _ name _ _ _ _ _ _ _; do
    # +0 forces numeric comparison: split() results compare as strings otherwise, and "9" < "10" is
    # false as a string, which would report a wrong minimum carrier count.
    min="$(awk -F'\t' 'NR>1{split($7,a,":"); split($8,b,":"); x=a[1]+0; y=b[1]+0;
                        m=(x<y?x:y); if (n=="" || m<n+0) n=m} END{print n+0}' "$OUT/$name/truth.tsv")"
    blocks="$(awk -F'\t' 'NR>1' "$OUT/$name/truth.tsv" | wc -l | tr -d ' ')"
    if [ -z "$min" ] || [ "$min" -lt 1 ]; then
      printf '  %-22s UNREPRESENTABLE: a block allele has %s panel carriers\n' "$name" "${min:-NA}"
      bad=$((bad + 1))
    else
      printf '  %-22s %s blocks, every truth allele carried by >= %s panel template(s)\n' "$name" "$blocks" "$min"
    fi
  done < <(plan_rows | tr '\t' '|')

  if [ "$bad" -ne 0 ]; then echo "verify FAILED: $bad problem(s); nothing may be simulated"; exit 1; fi
  echo "verify passed: every staged model is ordinary C4 and every truth is admissible"
fi

# =====================================================================================================
# Below this line nothing runs by default. These stages consume held-out data and are reached only by
# naming them in STAGES, after the protocol above has been reviewed.
# =====================================================================================================

# ---- reads: one simulation per plan, then NESTED subsets --------------------------------------------
if has_stage reads; then
  echo "=== simulating reads (held-out data) ==="
  while IFS='|' read -r _ name _ _ _ _ _ _ _; do
    d="$OUT/$name"
    mkdir -p "$d/cov$SOURCE_COVERAGE"
    "$PROBE" simulate --gfa "$d/staged.gfa" --reference "$REFERENCE" --donor MOSAIC \
      --depth "$SOURCE_COVERAGE" --read-length "$READ_LENGTH" --insert-mean "$INSERT_MEAN" \
      --insert-sd "$INSERT_SD" --sigmas "$SIGMAS" --error "$SUBSTITUTION_ERROR" \
      --outlier "$OUTLIER_WEIGHT" --seed "$SEED" --seed-salt "$name" \
      --out "$d/cov$SOURCE_COVERAGE" > "$d/sim.log" 2>&1
    for cov in $COVERAGES; do
      [ "$cov" = "$SOURCE_COVERAGE" ] && continue
      mkdir -p "$d/cov$cov"
      "$PROBE" subsample --reads1 "$d/cov$SOURCE_COVERAGE/reads_1.fq" \
        --reads2 "$d/cov$SOURCE_COVERAGE/reads_2.fq" --library "$d/cov$SOURCE_COVERAGE/library.txt" \
        --source-coverage "$SOURCE_COVERAGE" --coverage "$cov" --name-salt "$name" \
        --out "$d/cov$cov" > "$d/subsample.$cov.log" 2>&1
    done
    for cov in $COVERAGES; do record_experiment_signature "$name" "$cov" > /dev/null; done
    printf '  %-22s %s fragments at %sx, nested to %s\n' "$name" \
      "$(awk -F'\t' '$1=="fragments_written"{print $2}' "$d/cov$SOURCE_COVERAGE/simulation.tsv")" \
      "$SOURCE_COVERAGE" "$(echo $COVERAGES | tr ' ' ',')"
  done < <(plan_rows | tr '\t' '|')
fi

# ---- nested: prove the subsets really are nested, by read identity ----------------------------------
# The subsampler's nesting is a property of its selector, and the selector is tested on stubs elsewhere.
# This checks the ACTUAL read files this experiment will score: mates in step, and every fragment kept
# at a coverage kept at every higher one. A broken nested subset would make a coverage comparison a
# comparison of two different draws, which is exactly the confound the nesting exists to remove.
if has_stage nested; then
  echo "=== post-subsampling gate: mate synchronisation and nesting ==="
  nested_bad=0
  while IFS='|' read -r _ name _ _ _ _ _ _ _; do
    if ! python3 - "$OUT/$name" "$name" "$COVERAGES" <<'PY'
import sys
d, name, covs = sys.argv[1], sys.argv[2], sys.argv[3].split()
covs = sorted((int(c) for c in covs))

def names(path):
    out = []
    with open(path) as handle:
        for index, line in enumerate(handle):
            if index % 4 == 0:
                label = line[1:].strip()
                out.append(label.rsplit("/", 1)[0] if "/" in label else label)
    return out

problems = []
kept = {}
for cov in covs:
    first, second = names(f"{d}/cov{cov}/reads_1.fq"), names(f"{d}/cov{cov}/reads_2.fq")
    if len(first) != len(second):
        problems.append(f"cov{cov}: {len(first)} vs {len(second)} records -- mates out of step")
        continue
    mismatched = [i for i, (a, b) in enumerate(zip(first, second)) if a != b]
    if mismatched:
        problems.append(f"cov{cov}: {len(mismatched)} mate pairs disagree, first at record {mismatched[0]}")
    if len(set(first)) != len(first):
        problems.append(f"cov{cov}: duplicate fragment names")
    kept[cov] = set(first)

for lower, higher in zip(covs, covs[1:]):
    if lower in kept and higher in kept:
        escaped = kept[lower] - kept[higher]
        if escaped:
            problems.append(f"{lower}x is NOT a subset of {higher}x: {len(escaped)} fragments escaped")

if problems:
    print(f"  {name:<22} NESTING BROKEN")
    for problem in problems:
        print(f"    {problem}")
    raise SystemExit(1)
chain = " subset of ".join(f"{c}x({len(kept[c])})" for c in covs)
print(f"  {name:<22} mates in step; {chain}")
PY
    then nested_bad=$((nested_bad + 1)); fi
  done < <(plan_rows | tr '\t' '|')
  if [ "$nested_bad" -ne 0 ]; then
    echo "post-subsampling gate FAILED on $nested_bad plan(s); no coverage comparison is interpretable"
    exit 1
  fi
  echo "post-subsampling gate passed"
fi

# ---- call: one shared exact-scoring pass per plan and coverage --------------------------------------
if has_stage call; then
  echo "=== calling (staged graph, mosaics curated out) ==="
  printf 'plan\tcoverage\twall_s\tmax_rss_bytes\n' > "$OUT/runtime.tsv"
  if /usr/bin/time -l true >/dev/null 2>&1; then TIME_FLAG="-l"; else TIME_FLAG="-v"; fi
  call_bad=0
  while IFS='|' read -r _ name _ _ _ _ _ _ _; do
    d="$OUT/$name"
    for cov in $COVERAGES; do
      if signature_incomplete "$name" "$cov"; then
        printf '  %-22s cov%-3s INPUT MISSING -- rerun the reads stage\n' "$name" "$cov"
        call_bad=$((call_bad + 1)); continue
      fi
      # The inputs must still be what they were when the signature was recorded.
      if [ ! -f "$OUT/signatures/$(experiment_signature_id "$name" "$cov").tsv" ]; then
        printf '  %-22s cov%-3s INPUTS CHANGED since simulation -- no recorded signature matches\n' "$name" "$cov"
        call_bad=$((call_bad + 1)); continue
      fi
      rm -f "$d"/pending.*
      /usr/bin/time $TIME_FLAG -o "$d/time.$cov" \
        "$PROBE" genotype --gfa "$d/staged.gfa" --reference "$REFERENCE" \
          --reads1 "$d/cov$cov/reads_1.fq" --reads2 "$d/cov$cov/reads_2.fq" \
          --library "$d/cov$cov/library.txt" --out-prefix "$d/pending" \
          --state-limit "$STATE_LIMIT" --prefixes "$CALL_PREFIXES" --phase-audit "$PHASE_AUDIT" \
          --threads 2 --panel-curation "$CURATION" > "$d/pending.log" 2>&1 || true
      if [ ! -f "$d/pending.K0.mosaic_status.tsv" ]; then
        printf '  %-22s cov%-3s NO PER-K STATUS; previous outputs untouched\n' "$name" "$cov"
        rm -f "$d"/pending.*; call_bad=$((call_bad + 1)); continue
      fi
      rm -f "$d"/call.$cov.*
      for f in "$d"/pending.*; do mv "$f" "$d/call.$cov.${f##*/pending.}"; done
      # Written only after publication, and never on the failure path: last run's outputs must never sit
      # beside a signature claiming they match the current inputs.
      call_signature "$name" "$cov" > "$d/call.$cov.signature"
      if [ "$TIME_FLAG" = "-l" ]; then
        w="$(awk '/ real /{print $1; exit}' "$d/time.$cov")"
        r="$(awk '/maximum resident set size/{print $1}' "$d/time.$cov")"
      else
        w="$(awk -F': ' '/Elapsed \(wall clock\)/{n=split($2,p,":"); s=0; for(i=1;i<=n;i++) s=s*60+p[i]; print s}' "$d/time.$cov")"
        r="$(awk -F: '/Maximum resident set size/{gsub(/ /,"",$2); print $2*1024}' "$d/time.$cov")"
      fi
      printf '%s\t%s\t%s\t%s\n' "$name" "$cov" "${w:-NA}" "${r:-NA}" >> "$OUT/runtime.tsv"
      printf '  %-22s cov%-3s %s\n' "$name" "$cov" \
        "$(awk -F'\t' '$1=="status"{print $2}' "$d/call.$cov.K0.mosaic_status.tsv")"
    done
  done < <(plan_rows | tr '\t' '|')
  [ "$call_bad" -eq 0 ] || { echo "call stage FAILED on $call_bad run(s)"; exit 1; }
fi

# ---- complete: every declared run and output exists -------------------------------------------------
# A rate computed over whatever happened to succeed has a denominator nobody declared. The experiment
# fails on a missing run or output instead.
if has_stage complete; then
  echo "=== completeness gate ==="
  runs=0; records=0; missing=0
  while IFS='|' read -r _ name _ _ _ _ _ _ _; do
    d="$OUT/$name"
    for cov in $COVERAGES; do
      if [ -f "$d/call.$cov.signature" ]; then runs=$((runs + 1)); else
        printf '  %-22s cov%-3s MISSING RUN (no call signature)\n' "$name" "$cov"; missing=$((missing + 1))
      fi
      for k in ${CALL_PREFIXES//,/ }; do
        for suffix in mosaic_status.tsv mosaic_scores.tsv mosaic_candidates.tsv block_genotypes.tsv; do
          f="$d/call.$cov.K$k.$suffix"
          if [ -f "$f" ]; then
            [ "$suffix" = "mosaic_status.tsv" ] && records=$((records + 1))
          else
            printf '  %-22s cov%-3s K%-3s MISSING OUTPUT %s\n' "$name" "$cov" "$k" "$suffix"
            missing=$((missing + 1))
          fi
        done
      done
    done
  done < <(plan_rows | tr '\t' '|')
  printf '  %s shared runs (want %s), %s K-arm status records (want %s), %s missing output(s)\n' \
    "$runs" "$EXPECTED_RUNS" "$records" "$EXPECTED_STATUS" "$missing"
  if [ "$runs" != "$EXPECTED_RUNS" ] || [ "$records" != "$EXPECTED_STATUS" ] || [ "$missing" -ne 0 ]; then
    echo "completeness gate FAILED: the experiment is incomplete and no rate over it is interpretable"
    exit 1
  fi
  echo "completeness gate passed"
fi

# ---- eval: accuracy against a truth the caller never saw --------------------------------------------
if has_stage eval; then
  echo "=== CERTIFIED CALLS: accuracy against truth ==="
  printf 'plan\tcoverage\tK\tmin_gq\tblocks_reportable\tgenotype_support\tcandidate_scope\tfinal_candidate_scope\tsearch_exhaustive\tseparation_gap\tfragments_total\n' \
    > "$OUT/diagnostics.tsv"
  printf '%-22s %-4s %-4s %-18s %-7s %-7s %-9s %-12s %s\n' \
    plan cov K status exact dist truth_rank phase uniq_hap
  eval_bad=0
  while IFS='|' read -r _ name _ _ _ _ _ _ _; do
    d="$OUT/$name"
    for cov in $COVERAGES; do
      require_call_signature "$name" "$cov" || { eval_bad=$((eval_bad + 1)); continue; }
      for k in ${CALL_PREFIXES//,/ }; do
        s="$d/call.$cov.K$k.mosaic_status.tsv"
        b="$d/call.$cov.K$k.block_genotypes.tsv"
        sc="$d/call.$cov.K$k.mosaic_scores.tsv"
        "$PROBE" evaluate --run-prefix "$d/call.$cov.K$k" --out "$d/eval.$cov.K$k" \
          --scheme "c4-heldout-$name-cov$cov-K$k" --donor MOSAIC --reference "$REFERENCE" \
          --model-gfa "$d/staged.gfa" --truth-gfa "$d/staged.gfa" \
          --panel-curation "$CURATION" > "$d/eval.$cov.K$k.log" 2>&1 || true
        e="$d/eval.$cov.K$k.summary.tsv"
        get() { awk -F'\t' -v key="$1" '$1==key{print $2; exit}' "$2" 2>/dev/null; }
        col() { awk -F'\t' -v want="$1" 'NR==1{for(i=1;i<=NF;i++) if($i==want) c=i; next} c{print $c}' "$2" 2>/dev/null; }
        fa=no; [ -f "$d/call.$cov.K$k.mosaic_haplotypes.fa" ] && fa=yes
        printf '%-22s %-4s %-4s %-18s %-7s %-7s %-9s %-12s %s\n' \
          "$name" "$cov" "$k" "$(get status "$s")" "$(get blocks_exact "$e")" \
          "$(get unordered_distance_total "$e")" "$(get truth_fragment_score_rank "$e")" \
          "$(get phase_support "$s")" "$fa"
        # Diagnostics, exactly as frozen in the manifest: reported beside the accuracy above, never
        # merged into it, and with no threshold applied to any of them.
        min_gq="$(paste <(col reportable "$b") <(col conditional_GQ "$b") |
                  awk -F'\t' '$1=="true"{ if (n=="" || $2+0 < n+0) n=$2 } END{print (n==""?"NA":n)}')"
        reportable="$(col reportable "$b" | grep -c '^true$' || true)"
        support="$(col genotype_support "$b" | sort | uniq -c | awk '{printf "%s%s:%s", (NR>1?",":""), $2, $1}')"
        frags="$(col fragments_total "$sc" | head -1)"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
          "$name" "$cov" "$k" "${min_gq:-NA}" "${reportable:-0}" "${support:-NA}" \
          "$(get candidate_scope "$s")" "$(get final_candidate_scope "$s")" \
          "$(get proposal_search_exhaustive "$s")" "$(get numeric_separation_gap "$s")" \
          "${frags:-NA}" >> "$OUT/diagnostics.tsv"
      done
    done
  done < <(plan_rows | tr '\t' '|')
  [ "$eval_bad" -eq 0 ] || { echo "eval FAILED: $eval_bad run(s) refused interpretation"; exit 1; }
  echo "diagnostics: $OUT/diagnostics.tsv"
fi

# ---- concord: CONDITIONAL tied-winner block concordance ---------------------------------------------
# Reported separately from the certified calls above and never merged with them. This is a measurement
# of what tied candidates agree on; it is not a call, is not counted as one, and does not enter any
# accuracy rate. The strict policy is unchanged: a pruned search must not produce certified block
# consensus, so on C4 this stays a diagnostic.
if has_stage concord; then
  echo "=== CONDITIONAL (NOT CALLS): tied-winner block concordance ==="
  concord_bad=0
  while IFS='|' read -r _ name _ _ _ _ _ _ _; do
    d="$OUT/$name"
    for cov in $COVERAGES; do
      require_call_signature "$name" "$cov" || { concord_bad=$((concord_bad + 1)); continue; }
      for k in ${CALL_PREFIXES//,/ }; do
        python3 - "$d" "$name" "$cov" "$k" <<'PY'
import csv, os, sys
d, name, cov, k = sys.argv[1:5]
sc, cd = f"{d}/call.{cov}.K{k}.mosaic_scores.tsv", f"{d}/call.{cov}.K{k}.mosaic_candidates.tsv"
if not (os.path.exists(sc) and os.path.exists(cd)):
    print(f"  {name:<22} cov{cov:<3} K{k:<3} no scores/candidates"); raise SystemExit
rows = lambda p: list(csv.DictReader(open(p), delimiter='\t'))
win = {r["representative_beam_candidate"] for r in rows(sc) if r.get("in_winner_set") == "true"}
by_id = {c["candidate_id"]: c for c in rows(cd)}
sel = [by_id[w] for w in win if w in by_id]
if not sel:
    print(f"  {name:<22} cov{cov:<3} K{k:<3} winner set {len(win)}, none resolvable"); raise SystemExit
al = lambda s: [int(x) for x in s.split(",") if x != ""]
nb = len(al(sel[0]["hap1_alleles"]))
agree = [b for b in range(nb)
         if len({tuple(sorted((al(c["hap1_alleles"])[b], al(c["hap2_alleles"])[b]))) for c in sel}) == 1]
note = "trivial: one winner" if len(sel) == 1 else ("unanimous" if len(agree) == nb else "disagree somewhere")
print(f"  {name:<22} cov{cov:<3} K{k:<3} winners={len(sel):<3} blocks_agreeing={len(agree)}/{nb}  {note}")
PY
      done
    done
  done < <(plan_rows | tr '\t' '|')
  [ "$concord_bad" -eq 0 ] || { echo "concord FAILED: $concord_bad run(s) refused interpretation"; exit 1; }
fi

# ---- trace: proposal-only diagnostic, at the frozen width and the comparison width ------------------
if has_stage trace; then
  echo "=== proposal-only trace ==="
  printf '%-22s %-6s %-14s %-12s %s\n' plan width truth_in_beam full_beam_rank first_loss
  trace_bad=0
  while IFS='|' read -r _ name _ _ _ _ _ _ _; do
    d="$OUT/$name"
    require_call_signature "$name" "$SOURCE_COVERAGE" || { trace_bad=$((trace_bad + 1)); continue; }
    for width in "$STATE_LIMIT" "$DIAGNOSTIC_WIDTH"; do
      "$PROBE" trace --gfa "$d/staged.gfa" --reference "$REFERENCE" --donor MOSAIC \
        --truth-gfa "$d/staged.gfa" --panel-curation "$CURATION" --grouped 1 \
        --run "$d/call.$SOURCE_COVERAGE.K${CALL_PREFIXES##*,}" \
        --reads1 "$d/cov$SOURCE_COVERAGE/reads_1.fq" --reads2 "$d/cov$SOURCE_COVERAGE/reads_2.fq" \
        --library "$d/cov$SOURCE_COVERAGE/library.txt" --state-limit "$width" \
        --candidates "$TRACE_CANDIDATES" --out "$d/trace.$width" > "$d/trace.$width.log" 2>&1 || true
      f="$d/trace.$width.trace.summary.tsv"
      printf '%-22s %-6s %-14s %-12s %s\n' "$name" "$width" \
        "$(awk -F'\t' '$1=="truth_in_full_beam"{print $2}' "$f" 2>/dev/null)" \
        "$(awk -F'\t' '$1=="truth_full_beam_rank"{print $2}' "$f" 2>/dev/null)" \
        "$(awk -F'\t' '$1=="first_block_truth_not_retained"{print $2}' "$f" 2>/dev/null)"
    done
  done < <(plan_rows | tr '\t' '|')
  [ "$trace_bad" -eq 0 ] || { echo "trace FAILED: $trace_bad plan(s) refused interpretation"; exit 1; }
fi

echo "done: $OUT"
