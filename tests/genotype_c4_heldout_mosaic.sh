#!/usr/bin/env bash
# The HELD-OUT C4 synthetic-mosaic validation.
#
#   tests/genotype_c4_heldout_mosaic.sh --check-manifest          # the manifest gate alone (ctest)
#   BUILD=<build dir> tests/genotype_c4_heldout_mosaic.sh <out-dir>
#
# THE DEFAULT STAGES DO NOT SIMULATE, CALL OR EVALUATE ANYTHING. They are "manifest stage verify":
# check that the committed plan table is what the declared draw rule produces, build the ten mosaic
# truths, and prove each one is admissible. Simulation waits for review, and is reached only by asking
# for it explicitly:
#
#   STAGES="reads call eval concord trace" BUILD=<build> tests/genotype_c4_heldout_mosaic.sh <out-dir>
#
# That split is the point of this harness. A protocol that can be half-run by accident is not frozen,
# and the whole value of held-out data is that nothing touched it before the protocol was fixed.
#
# Everything the run obeys comes from tests/genotype_c4_heldout_mosaic_manifest.tsv, which was committed
# before any of this ran. Nothing here re-selects a caller setting: width, K, the diagnostic comparison
# width and the phase audit are read from the manifest and asserted to be the frozen production values.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MANIFEST="$REPO/tests/genotype_c4_heldout_mosaic_manifest.tsv"
ORIENTATION="$REPO/tests/genotype_c4_template_orientation.tsv"
PANEL_ALLELES="$REPO/tests/genotype_c4_panel_alleles.tsv"
DONORS="$REPO/tests/genotype_development_donors.tsv"
DRAW="$REPO/scripts/draw_heldout_mosaic_plans.py"

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

  # Every metric is declared before any result exists, and each names the file and key it comes from.
  local metrics
  metrics="$(awk -F'\t' '$1=="metric"' "$MANIFEST" | wc -l | tr -d ' ')"
  if [ "$metrics" -ge 20 ]; then echo "  $metrics metrics frozen, each naming its source key"
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

mkdir -p "$OUT"
C4="$OUT/c4.gfa"
[ -f "$C4" ] || gzip -dc "$REPO/tests/real_data/c4.gfa.gz" > "$C4"
CURATION="$OUT/mosaic_exclude.tsv"
cat > "$CURATION" <<'EOF'
# The two synthetic mosaics are removed from the PANEL, not from the file: reads are simulated from them
# and truth is projected through them, but no candidate may reach them without recombining real panel
# templates. For the panel_control plans the two paths are copies of real panel rows, and curating them
# out leaves those rows in the panel -- which is exactly what makes a control a control.
MOSAIC#1#synthetic:0-0
MOSAIC#2#synthetic:0-0
EOF

# One signature over everything a call consumes, always recomputed from the files on disk. A stored
# signature is only ever compared against, never trusted as a description of what is there now.
input_bundle_signature() {
  local d="$1" name="$2" cov="$3" f
  for f in "$d/staged.gfa" "$d/cov$cov/reads_1.fq" "$d/cov$cov/reads_2.fq" "$d/cov$cov/library.txt"; do
    [ -f "$f" ] || { printf 'MISSING:%s\n' "$(basename "$f")"; return 0; }
  done
  [ -f "$PROBE" ] || { printf 'MISSING:probe\n'; return 0; }
  {
    shasum -a 256 < "$d/staged.gfa"
    shasum -a 256 < "$d/cov$cov/reads_1.fq"
    shasum -a 256 < "$d/cov$cov/reads_2.fq"
    shasum -a 256 < "$d/cov$cov/library.txt"
    shasum -a 256 < "$PROBE"
    printf 'salt=%s seed=%s cov=%s src=%s len=%s ins=%s sd=%s sig=%s err=%s out=%s\n' \
      "$name" "$SEED" "$cov" "$SOURCE_COVERAGE" "$READ_LENGTH" "$INSERT_MEAN" "$INSERT_SD" \
      "$SIGMAS" "$SUBSTITUTION_ERROR" "$OUTLIER_WEIGHT"
  } | shasum -a 256 | cut -d' ' -f1
}

{
  printf 'key\tvalue\n'
  printf 'git_commit\t%s\n' "$(git -C "$REPO" rev-parse HEAD)"
  printf 'git_tree_dirty\t%s\n' "$(git -C "$REPO" status --porcelain | grep -q . && echo true || echo false)"
  printf 'probe_sha256\t%s\n' "$([ -f "$PROBE" ] && shasum -a 256 "$PROBE" | cut -d' ' -f1 || echo MISSING)"
  printf 'script_sha256\t%s\n' "$(shasum -a 256 "${BASH_SOURCE[0]}" | cut -d' ' -f1)"
  printf 'manifest_sha256\t%s\n' "$(shasum -a 256 "$MANIFEST" | cut -d' ' -f1)"
  printf 'draw_script_sha256\t%s\n' "$(shasum -a 256 "$DRAW" | cut -d' ' -f1)"
  printf 'orientation_sha256\t%s\n' "$(shasum -a 256 "$ORIENTATION" | cut -d' ' -f1)"
  printf 'panel_alleles_sha256\t%s\n' "$(shasum -a 256 "$PANEL_ALLELES" | cut -d' ' -f1)"
  printf 'c4_gfa_sha256\t%s\n' "$(shasum -a 256 "$REPO/tests/real_data/c4.gfa.gz" | cut -d' ' -f1)"
  printf 'curation_sha256\t%s\n' "$(shasum -a 256 "$CURATION" | cut -d' ' -f1)"
  printf 'stages_run\t%s\n' "$STAGES"
} > "$OUT/signature.tsv"
echo "signature: $OUT/signature.tsv"

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
    shasum -a 256 "$d/staged.gfa" | cut -d' ' -f1 > "$d/staged.gfa.sha"
    printf '  %-22s %-20s %s\n' "$name" "$class" "$(tail -1 "$d/stage.log")"
  done < <(plan_rows | tr '\t' '|')
fi

# ---- verify: is every staged truth admissible? ------------------------------------------------------
# Everything here runs BEFORE any read exists, and refuses rather than reports. What it cannot check
# yet is the full model-equivalence fingerprint: genotype_probe preflight takes no --panel-curation, so
# the curated staged model cannot be built without a call. That gate therefore stays where the
# development harness has it -- at the call stage, before any result is read -- and is NOT claimed here.
if has_stage verify; then
  echo "=== verify (no reads, no calls) ==="
  bad=0

  # 1. Staging added the two mosaic paths and NOTHING else. Each staged graph is the reference-sorted
  # C4 plus two P lines, so stripping those lines must leave every plan holding byte-identical graphs.
  # A donor-dependent difference here would mean the model each plan is scored against is not the same
  # model, which would make the ten results incomparable.
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

  # 2. The declared off-panel requirement, and how many homologues actually are panel rows. mosaic-truth
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

  # 3. Representability: every block allele of both truth haplotypes must be carried by at least one
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
  echo "verify passed: every staged truth is admissible"
fi

# =====================================================================================================
# Below this line nothing runs by default. These stages consume held-out data and are reached only by
# naming them in STAGES, after the protocol above has been reviewed.
# =====================================================================================================

# ---- reads: one simulation per plan, then NESTED subsets ---------------------------------------------
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
    for cov in $COVERAGES; do input_bundle_signature "$d" "$name" "$cov" > "$d/cov$cov/bundle.sha"; done
    printf '  %-22s %s fragments at %sx, nested to %s\n' "$name" \
      "$(awk -F'\t' '$1=="fragments_written"{print $2}' "$d/cov$SOURCE_COVERAGE/simulation.tsv")" \
      "$SOURCE_COVERAGE" "$(echo $COVERAGES | tr ' ' ',')"
  done < <(plan_rows | tr '\t' '|')
fi

check_binding() {
  local d="$1" name="$2" cov="$3" now stored
  now="$(input_bundle_signature "$d" "$name" "$cov")"
  case "$now" in MISSING:*) printf '  %-22s INPUT MISSING (%s)\n' "$name" "$now"; return 1 ;; esac
  stored="$(cat "$d/cov$cov/bundle.sha" 2>/dev/null || true)"
  [ -n "$stored" ] || { printf '  %-22s NO INPUT SIGNATURE\n' "$name"; return 1; }
  [ "$now" = "$stored" ] || { printf '  %-22s INPUTS CHANGED since simulation\n' "$name"; return 1; }
  return 0
}

# ---- call: one shared exact-scoring pass per plan and coverage ---------------------------------------
if has_stage call; then
  echo "=== calling (staged graph, mosaics curated out) ==="
  : > "$OUT/runtime.tsv"; printf 'plan\tcoverage\twall_s\tmax_rss_bytes\n' > "$OUT/runtime.tsv"
  if /usr/bin/time -l true >/dev/null 2>&1; then TIME_FLAG="-l"; else TIME_FLAG="-v"; fi
  while IFS='|' read -r _ name _ _ _ _ _ _ _; do
    d="$OUT/$name"
    for cov in $COVERAGES; do
      check_binding "$d" "$name" "$cov" || continue
      rm -f "$d"/pending.*
      /usr/bin/time $TIME_FLAG -o "$d/time.$cov" \
        "$PROBE" genotype --gfa "$d/staged.gfa" --reference "$REFERENCE" \
          --reads1 "$d/cov$cov/reads_1.fq" --reads2 "$d/cov$cov/reads_2.fq" \
          --library "$d/cov$cov/library.txt" --out-prefix "$d/pending" \
          --state-limit "$STATE_LIMIT" --prefixes "$CALL_PREFIXES" --phase-audit "$PHASE_AUDIT" \
          --threads 2 --panel-curation "$CURATION" > "$d/pending.log" 2>&1 || true
      if [ ! -f "$d/pending.K0.mosaic_status.tsv" ]; then
        printf '  %-22s cov%-3s NO PER-K STATUS; previous outputs untouched\n' "$name" "$cov"
        rm -f "$d"/pending.*; continue
      fi
      rm -f "$d"/call.$cov.*
      for f in "$d"/pending.*; do mv "$f" "$d/call.$cov.${f##*/pending.}"; done
      input_bundle_signature "$d" "$name" "$cov" > "$d/call.$cov.bundle.sha"
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
fi

# ---- eval: accuracy against a truth the caller never saw --------------------------------------------
if has_stage eval; then
  echo "=== CERTIFIED CALLS: accuracy against truth ==="
  printf '%-22s %-4s %-2s %-18s %-7s %-7s %-9s %-12s %s\n' \
    plan cov K status exact dist truth_rank phase uniq_hap
  while IFS='|' read -r _ name _ _ _ _ _ _ _; do
    d="$OUT/$name"
    for cov in $COVERAGES; do
      for k in ${CALL_PREFIXES//,/ }; do
        s="$d/call.$cov.K$k.mosaic_status.tsv"; [ -f "$s" ] || continue
        "$PROBE" evaluate --run-prefix "$d/call.$cov.K$k" --out "$d/eval.$cov.K$k" \
          --scheme "c4-heldout-$name-cov$cov-K$k" --donor MOSAIC --reference "$REFERENCE" \
          --model-gfa "$d/staged.gfa" --truth-gfa "$d/staged.gfa" \
          --panel-curation "$CURATION" > "$d/eval.$cov.K$k.log" 2>&1 || true
        e="$d/eval.$cov.K$k.summary.tsv"
        get() { awk -F'\t' -v key="$1" '$1==key{print $2; exit}' "$2" 2>/dev/null; }
        fa=no; [ -f "$d/call.$cov.K$k.mosaic_haplotypes.fa" ] && fa=yes
        printf '%-22s %-4s %-2s %-18s %-7s %-7s %-9s %-12s %s\n' \
          "$name" "$cov" "$k" "$(get status "$s")" "$(get blocks_exact "$e")" \
          "$(get unordered_distance_total "$e")" "$(get truth_fragment_score_rank "$e")" \
          "$(get phase_support "$s")" "$fa"
      done
    done
  done < <(plan_rows | tr '\t' '|')
fi

# ---- concord: CONDITIONAL tied-winner block concordance ---------------------------------------------
# Reported separately from the certified calls above and never merged with them. This is a measurement
# of what tied candidates agree on; it is not a call, is not counted as one, and does not enter any
# accuracy rate. The strict policy is unchanged: a pruned search must not produce certified block
# consensus, so on C4 this stays a diagnostic.
if has_stage concord; then
  echo "=== CONDITIONAL (NOT CALLS): tied-winner block concordance ==="
  while IFS='|' read -r _ name _ _ _ _ _ _ _; do
    d="$OUT/$name"
    for cov in $COVERAGES; do
      for k in ${CALL_PREFIXES//,/ }; do
        [ -f "$d/call.$cov.K$k.mosaic_scores.tsv" ] || continue
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
fi

# ---- trace: proposal-only diagnostic, at the frozen width and the comparison width -------------------
if has_stage trace; then
  echo "=== proposal-only trace ==="
  printf '%-22s %-6s %-14s %-12s %s\n' plan width truth_in_beam full_beam_rank first_loss
  while IFS='|' read -r _ name _ _ _ _ _ _ _; do
    d="$OUT/$name"
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
fi

echo "done: $OUT"
