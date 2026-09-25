#!/usr/bin/env bash
# Gate for the prepared-not-launched final C4 campaign: the plan is derived from the committed manifests alone
# (484 jobs: 224 LZO, 220 LOO, 40 synthetic mosaic, each cohort at 2x, 5x, 10x and 30x), the reference donor
# HG00096 never enters LOO, every call names the v6 campaign-local profile and no LOO call names its truth, the
# scientific provenance is the v4 profile through the v5/v6 profile gates, the plan is bound to a clean build of
# the production commit (or of a descendant with identical production sources) and to a fresh directory named
# from it, legitimate no-calls are completed outcomes, resuming refuses a changed binary, profile or scheduler,
# and launching without the reviewed plan's hash is refused. Nothing is simulated, called or evaluated.
#
# The build under test is a throwaway clone of this repository at the production commit with stand-in binaries,
# so the provenance checks are exercised both ways whatever state the tree running the suite is in.
#
#   tests/genotype_production_campaign_gate.sh <scratch dir>
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${1:?usage: genotype_production_campaign_gate.sh <scratch dir>}"
if [ ! -d "$REPO/tests/results/genotype_coverage_c4_held_out_abb35a0" ]; then
  echo "skipped: the frozen held-out read package under tests/results is absent (gitignored)."; exit 0
fi
HARNESS="$REPO/tests/genotype_production_c4_campaign.sh"
PRODUCTION_COMMIT="$(awk -F= '/^PRODUCTION_COMMIT=/{print $2}' "$HARNESS")"
if ! git -C "$REPO" cat-file -e "$PRODUCTION_COMMIT^{commit}" 2> /dev/null; then
  echo "skipped: the production commit $PRODUCTION_COMMIT is not in this repository's history."; exit 0
fi
rm -rf "$WORK"; mkdir -p "$WORK"
fail=0
same() { if [ "$2" = "$3" ]; then echo "  ok: $1 = $2"; else echo "  GATE FAILED: $1 is '$2', expected '$3'"; fail=1; fi; }
key() { awk -F'\t' -v k="$2" '$1==k{print $2; exit}' "$1"; }
sha() { shasum -a 256 "$1" | cut -d' ' -f1; }
G=(git -c user.name=gate -c user.email=gate@invalid -c advice.detachedHead=false)

# ---- stand-in production builds: clean clones with a build/ holding executable placeholders
fake_build() {   # fake_build <dir> <commit> <tag>
  "${G[@]}" clone -q --shared "$REPO" "$1" && "${G[@]}" -C "$1" checkout -q --detach "$2" || return 1
  mkdir -p "$1/build"
  printf '#!/bin/sh\n# stand-in %s\nexit 99\n' "$3" > "$1/build/panvar"
  printf '#!/bin/sh\n# stand-in evaluator %s\nexit 99\n' "$3" > "$1/build/genotype_probe"
  chmod +x "$1/build/panvar" "$1/build/genotype_probe"
}
P7="${PRODUCTION_COMMIT:0:7}"
fake_build "$WORK/prod" "$PRODUCTION_COMMIT" production
PLAN="$WORK/genotype_c4_production_$P7"

echo "== plan"
bash "$HARNESS" plan "$PLAN" "$WORK/prod/build" > "$WORK/plan.log" 2>&1
same "plan exit" "$?" "0"
S="$PLAN/settings.tsv"
same "jobs" "$(key "$S" jobs)" "484"
same "LZO jobs" "$(key "$S" jobs_lzo)" "224"
same "LOO jobs" "$(key "$S" jobs_loo)" "220"
same "mosaic jobs" "$(key "$S" jobs_mosaic)" "40"
same "jobs.tsv rows" "$(($(wc -l < "$PLAN/jobs.tsv") - 1))" "484"
same "jobs_sha256 is jobs.tsv's" "$(key "$S" jobs_sha256)" "$(sha "$PLAN/jobs.tsv")"
same "launched" "$(key "$S" launched)" "false"
same "HG00096 LOO jobs" "$(awk -F'\t' '$2=="loo" && $3=="HG00096"' "$PLAN/jobs.tsv" | wc -l | tr -d ' ')" "0"
same "coverages per subject" "$(awk -F'\t' 'NR>1{print $4}' "$PLAN/jobs.tsv" | sort -un | tr '\n' ' ')" "2 5 10 30 "
# Path columns: model, reads, library, curation, truth, the three sources, and the output prefix.
same "every path in jobs.tsv is absolute" "$(awk -F'\t' 'NR>1{for(i=7;i<=18;i++) if(i!=16 && i!=17 && $i!="none" && $i!="harness" && substr($i,1,1)!="/") n++} END{print n+0}' "$PLAN/jobs.tsv")" "0"
(cd "$WORK" && bash "$HARNESS" plan "rel/genotype_c4_production_$P7" "$WORK/prod/build" > "$WORK/plan_rel.log" 2>&1)
same "a relative plan directory yields absolute job paths" "$(awk -F'\t' 'NR==2{print substr($18,1,1)}' "$WORK/rel/genotype_c4_production_$P7/jobs.tsv")" "/"
same "every job writes inside the plan directory" "$(awk -F'\t' -v p="$PLAN/" 'NR>1 && index($18, p)!=1' "$PLAN/jobs.tsv" | wc -l | tr -d ' ')" "0"

echo "== profiles and provenance"
same "production profile" "$(key "$S" production_profile)" "profiles/genotype_c4_v6_campaign_local.tsv"
same "production profile sha256" "$(key "$S" production_profile_sha256)" "$(sha "$REPO/profiles/genotype_c4_v6_campaign_local.tsv")"
same "production profile name" "$(key "$S" production_profile_name)" "c4-panel-anchored-v6-campaign-local"
same "scientific profile" "$(key "$S" science_profile)" "profiles/genotype_c4_v4.tsv"
same "scientific profile sha256" "$(key "$S" science_profile_sha256)" "$(sha "$REPO/profiles/genotype_c4_v4.tsv")"
same "no v1 profile bound in settings" "$(grep -v '^estimate_provenance' "$S" | grep -c -E 'genotype_c4(_campaign_local)?\.tsv')" "0"
same "production commit" "$(key "$S" production_commit)" "$PRODUCTION_COMMIT"
same "build commit" "$(key "$S" build_commit)" "$PRODUCTION_COMMIT"
same "build tree clean" "$(key "$S" build_tree_clean)" "true"
same "production binary sha256" "$(key "$S" production_binary_sha256)" "$(sha "$WORK/prod/build/panvar")"
same "evaluator binary sha256" "$(key "$S" evaluator_binary_sha256)" "$(sha "$WORK/prod/build/genotype_probe")"
for f in tests/genotype_production_c4_campaign.sh scripts/c4_production_campaign.py scripts/consolidate_c4_production.py \
         tests/genotype_c4_heldout_mosaic.sh tests/genotype_v5_profile_gate.sh tests/genotype_v6_profile_gate.sh; do
  same "$f hash recorded" "$(key "$S" "${f}_sha256")" "$(sha "$REPO/$f")"
done
same "concurrent jobs" "$(key "$S" concurrent_jobs)" "1"
same "threads per call" "$(key "$S" threads_per_call)" "8"
same "estimate labelled historical" "$(head -1 "$PLAN/estimate.tsv" | grep -c 'HISTORICAL, NOT A v6 PREDICTION')" "1"
same "estimate provenance recorded" "$(key "$S" estimate_provenance | cut -c1-10)" "HISTORICAL"
# The scientific chain, each step exactly its reviewed differences (the v5/v6 gates define them).
science() { grep -v '^#' "$1" | grep -v '^$' | grep -v -E '^(profile|threads|footprint_limit_bytes|time_limit_seconds)'$'\t' | sort; }
same "v4 campaign-local carries the science of v4" "$(diff <(science "$REPO/profiles/genotype_c4_v4.tsv") <(science "$REPO/profiles/genotype_c4_v4_campaign_local.tsv") > /dev/null && echo identical || echo DIFFERS)" "identical"
bash "$REPO/tests/genotype_v5_profile_gate.sh" > "$WORK/v5.log" 2>&1; same "v5 profile gate" "$?" "0"
bash "$REPO/tests/genotype_v6_profile_gate.sh" > "$WORK/v6.log" 2>&1; same "v6 profile gate" "$?" "0"
for kv in "repair_tie_branching 1" "repair_tie_max_branches 16" "representation_guard 1" "repair_rounds 32"; do
  for p in genotype_c4_v4 genotype_c4_v6_campaign_local; do
    same "$p ${kv% *}" "$(key "$REPO/profiles/$p.tsv" "${kv% *}")" "${kv#* }"
  done
done
# The historical chain v1 -> v2 -> v3 -> v4 is unchanged (the frozen v1 profiles are never reinterpreted).
v2() { grep -v '^#' "$1" | grep -v -E '^(profile|repair_tie_branching|repair_tie_max_branches)'$'\t' | sort; }
v3() { grep -v '^#' "$1" | grep -v -E '^(profile|representation_guard)'$'\t' | sort; }
v4() { grep -v '^#' "$1" | grep -v -E '^(profile|repair_rounds)'$'\t' | sort; }
for v in genotype_c4 genotype_c4_campaign_local; do
  p1="$REPO/profiles/$v.tsv"; p2="$REPO/profiles/${v/genotype_c4/genotype_c4_v2}.tsv"
  p3="$REPO/profiles/${v/genotype_c4/genotype_c4_v3}.tsv"; p4="$REPO/profiles/${v/genotype_c4/genotype_c4_v4}.tsv"
  same "$v v2 = v1 apart from tie branching" "$(diff <(v2 "$p1") <(v2 "$p2") > /dev/null && echo identical || echo DIFFERS)" "identical"
  same "$v v3 = v2 apart from the guard" "$(diff <(v3 "$p2") <(v3 "$p3") > /dev/null && echo identical || echo DIFFERS)" "identical"
  same "$v v4 = v3 apart from the round limit" "$(diff <(v4 "$p3") <(v4 "$p4") > /dev/null && echo identical || echo DIFFERS)" "identical"
  same "$v v1 sets no tie branching or guard" "$(grep -c -E '^(repair_tie_branching|representation_guard)'$'\t' "$p1")" "0"
done
same "v1 campaign-local footprint limit (frozen)" "$(key "$REPO/profiles/genotype_c4_campaign_local.tsv" footprint_limit_bytes)" "3758096384"

echo "== call commands and outcomes (the runner itself)"
python3 - "$REPO" "$PLAN" "$WORK" <<'EOF' || fail=1
import csv, importlib.util, pathlib, sys
repo, plan, work = map(pathlib.Path, sys.argv[1:4])
spec = importlib.util.spec_from_file_location("runner", repo / "scripts" / "c4_production_campaign.py")
runner = importlib.util.module_from_spec(spec); spec.loader.exec_module(runner)
bad = []
def same(what, got, want):
    print(("  ok: " if got == want else "  GATE FAILED: ") + f"{what} = {got!r}" + ("" if got == want else f", expected {want!r}"))
    if got != want: bad.append(what)
jobs = runner.read_jobs(plan)
build = work / "prod" / "build"
commands = {j["job"]: list(map(str, runner.call_command(j, build, 8))) for j in jobs}
profile = str(repo / "profiles" / "genotype_c4_v6_campaign_local.tsv")
same("every call names the v6 campaign-local profile", sum(c[c.index("--profile") + 1] == profile for c in commands.values()), len(jobs))
loo = [j for j in jobs if j["campaign"] == "loo"]
same("LOO calls naming their truth graph", sum(j["truth_gfa"] in commands[j["job"]] for j in loo), 0)
same("LOO model is the donor's leave-one-out graph", sum(j["model_gfa"].endswith(f"/{j['subject']}.loo.gfa") for j in loo), len(loo))
leaky = dict(loo[0]); leaky["model_gfa"] = leaky["truth_gfa"]
try:
    runner.call_command(leaky, build, 8); same("a LOO call naming its truth is refused", "accepted", "refused")
except runner.Defect:
    same("a LOO call naming its truth is refused", "refused", "refused")
# Legitimate no-calls are completed outcomes with their evidence, and never carry genotype output.
def job_with(name, status, files):
    prefix = work / "outcomes" / name
    prefix.parent.mkdir(parents=True, exist_ok=True)
    pathlib.Path(str(prefix) + ".mosaic_status.tsv").write_text("key\tvalue\n" + "".join(f"{k}\t{v}\n" for k, v in status.items()))
    for suffix in files: pathlib.Path(str(prefix) + suffix).write_text("x\n")
    return {"job": name, "campaign": "loo", "subject": "X", "coverage": "30", "out_prefix": str(prefix)}
tie_files = [".repair_rounds.tsv", ".tie_branches.tsv", ".tie_terminal_ties.tsv"]
cases = {
    "strict": ({"status": "COMPLETE", "call_class": "strict_call", "no_call_reason": "none", "locus_filter": "PASS"},
               [".genotypes.tsv", ".block_genotypes.tsv", ".mosaic_haplotypes.fa", ".repair_rounds.tsv"], "strict_call"),
    "equivalent": ({"status": "CONDITIONAL_GENOTYPE_EQUIVALENT", "call_class": "genotype_equivalent_call", "no_call_reason": "none",
                    "search_conditional": "tie_branched", "locus_filter": "PASS"},
                   [".genotypes.tsv", ".block_genotypes.tsv", ".equivalence_set.tsv", ".repair_rounds.tsv"] + tie_files[1:], "genotype_equivalent_call"),
    "branch_cap": ({"status": "REPAIR_INCOMPLETE", "call_class": "no_call", "no_call_reason": "tie_branching_incomplete",
                    "search_conditional": "tie_branched", "locus_filter": "NA"}, tie_files, "no_call"),
    "divergent": ({"status": "RESCORE_INCOMPLETE", "call_class": "no_call", "no_call_reason": "divergent_tie",
                   "search_conditional": "tie_branched", "locus_filter": "NA"}, tie_files + [".equivalence_set.tsv"], "no_call"),
    "branch_cap_with_genotypes": ({"status": "REPAIR_INCOMPLETE", "call_class": "no_call", "no_call_reason": "tie_branching_incomplete",
                                   "search_conditional": "tie_branched", "locus_filter": "NA"}, tie_files + [".genotypes.tsv"], "defect"),
    "branch_cap_with_fasta": ({"status": "REPAIR_INCOMPLETE", "call_class": "no_call", "no_call_reason": "tie_branching_incomplete",
                               "search_conditional": "tie_branched", "locus_filter": "NA"}, tie_files + [".mosaic_haplotypes.fa"], "defect"),
    "divergent_with_block_genotypes": ({"status": "RESCORE_INCOMPLETE", "call_class": "no_call", "no_call_reason": "divergent_tie",
                                        "search_conditional": "tie_branched", "locus_filter": "NA"},
                                       tie_files + [".equivalence_set.tsv", ".block_genotypes.tsv"], "defect"),
    "branch_cap_without_tie_tables": ({"status": "REPAIR_INCOMPLETE", "call_class": "no_call", "no_call_reason": "tie_branching_incomplete",
                                       "search_conditional": "tie_branched", "locus_filter": "NA"}, [".repair_rounds.tsv"], "defect"),
    "divergent_not_tie_branched": ({"status": "RESCORE_INCOMPLETE", "call_class": "no_call", "no_call_reason": "divergent_tie",
                                    "search_conditional": "none", "locus_filter": "NA"}, tie_files + [".equivalence_set.tsv"], "defect"),
}
made = {}
for name, (status, files, want) in cases.items():
    made[name] = job_with(name, status, files)
    same(f"classify {name}", runner.classify(made[name])[0], want)
same("tie_branching_incomplete detail", runner.classify(made["branch_cap"])[1], "tie_branching_incomplete")
same("divergent_tie detail", runner.classify(made["divergent"])[1], "divergent_tie")
# Resumability: a recorded no-call is finished (never rerun or retried); a pending or failed job is not; the
# ledger reloads what a previous launch recorded.
ledger = runner.Ledger(work / "outcomes", list(made.values()))
for name in ("strict", "equivalent", "branch_cap", "divergent"):
    ledger.record(made[name], *runner.classify(made[name]), 1)
ledger.record(made["branch_cap_with_genotypes"], "failed", "exit 139", 1)
again = runner.Ledger(work / "outcomes", list(made.values()))
same("finished after reload", sorted(n for n in made if again.finished(n)), ["branch_cap", "divergent", "equivalent", "strict"])
same("an external failure is resumed, not finished", again.finished("branch_cap_with_genotypes"), False)
same("a never-run job is pending", again.rows["branch_cap_without_tie_tables"]["outcome"], "pending")
sys.exit(1 if bad else 0)
EOF

echo "== launch refusals"
bash "$HARNESS" launch "$WORK/prod/build" "$PLAN" > "$WORK/launch.log" 2>&1
same "launch without approval" "$?" "2"
APPROVED="$(sha "$PLAN/jobs.tsv")"
GENOTYPE_C4_CAMPAIGN_APPROVED=0000 bash "$HARNESS" launch "$WORK/prod/build" "$PLAN" > "$WORK/launch_wrong.log" 2>&1
same "launch with another plan's hash" "$?" "2"
fake_build "$WORK/other" "$PRODUCTION_COMMIT" other
GENOTYPE_C4_CAMPAIGN_APPROVED="$APPROVED" bash "$HARNESS" launch "$WORK/other/build" "$PLAN" > "$WORK/launch_binary.log" 2>&1
same "launch with a binary that is not the plan's" "$?" "2"
same "  refused for the binary" "$(grep -c "production_binary_sha256" "$WORK/launch_binary.log")" "1"
# A resume whose previous launch ran another profile is refused before anything is staged.
printf "key\tvalue\npanvar_sha256\t%s\nprofile\tprofiles/genotype_c4_v6_campaign_local.tsv\nprofile_sha256\tstale\njobs\t1\nthreads\t8\njobs_sha256\t%s\n" \
  "$(sha "$WORK/prod/build/panvar")" "$APPROVED" > "$PLAN/launch.tsv"
GENOTYPE_C4_CAMPAIGN_APPROVED="$APPROVED" bash "$HARNESS" launch "$WORK/prod/build" "$PLAN" > "$WORK/launch_resume.log" 2>&1
same "resume after a profile change" "$?" "2"
same "  refused as a resume" "$(grep -c "resume refused: profile_sha256" "$WORK/launch_resume.log")" "1"
rm -f "$PLAN/launch.tsv"
same "nothing launched or staged" "$(find "$PLAN" \( -name '*.mosaic_status.tsv' -o -name 'c4.gfa' -o -name 'inputs.tsv' \) | wc -l | tr -d ' ')" "0"

echo "== plan refusals"
refused() {   # refused <what> <log> <command...>
  local what="$1" log="$2"; shift 2
  "$@" > "$log" 2>&1; local code=$?
  same "$what refused" "$([ "$code" -ne 0 ] && grep -q REFUSED "$log" && echo refused || echo "accepted (exit $code)")" "refused"
}
bash "$HARNESS" plan "$WORK/x/genotype_c4_production_$P7" > "$WORK/r_nobuild.log" 2>&1
same "a plan without a build refused" "$([ $? -ne 0 ] && [ ! -e "$WORK/x" ] && echo refused || echo accepted)" "refused"
refused "a plan named for another build" "$WORK/r_name.log" bash "$HARNESS" plan "$WORK/genotype_c4_production_bc5dd97" "$WORK/prod/build"
refused "a plan into the frozen v1 campaign" "$WORK/r_old.log" bash "$HARNESS" plan "$REPO/tests/results/genotype_c4_production_bc5dd97/campaign" "$WORK/prod/build"
mkdir -p "$WORK/used/genotype_c4_production_$P7"; : > "$WORK/used/genotype_c4_production_$P7/completeness.tsv"
refused "a plan into a directory a campaign ran in" "$WORK/r_used.log" bash "$HARNESS" plan "$WORK/used/genotype_c4_production_$P7" "$WORK/prod/build"
refused "two concurrent jobs without a calibration" "$WORK/r_jobs.log" env CONCURRENT_JOBS=2 bash "$HARNESS" plan "$WORK/j/genotype_c4_production_$P7" "$WORK/prod/build"
echo "stray" > "$WORK/other/stray.txt"
refused "a dirty build tree" "$WORK/r_dirty.log" bash "$HARNESS" plan "$WORK/d/genotype_c4_production_$P7" "$WORK/other/build"
rm -f "$WORK/other/stray.txt"
# A descendant: accepted when only non-production files changed, refused once a production source changes.
echo "# gate" >> "$WORK/other/tests/genotype_production_c4_timing.tsv"
"${G[@]}" -C "$WORK/other" commit -q -am "gate: a non-production change"
D="$("${G[@]}" -C "$WORK/other" rev-parse HEAD)"
bash "$HARNESS" plan "$WORK/desc/genotype_c4_production_${D:0:7}" "$WORK/other/build" > "$WORK/r_desc.log" 2>&1
same "a descendant with identical production sources is accepted" "$?" "0"
same "  its build commit is recorded" "$(key "$WORK/desc/genotype_c4_production_${D:0:7}/settings.tsv" build_commit)" "$D"
echo "// gate" >> "$WORK/other/src/genotype_command.cpp"
"${G[@]}" -C "$WORK/other" commit -q -am "gate: a production change"
D2="$("${G[@]}" -C "$WORK/other" rev-parse HEAD)"
refused "a descendant whose production sources differ" "$WORK/r_src.log" bash "$HARNESS" plan "$WORK/src/genotype_c4_production_${D2:0:7}" "$WORK/other/build"
PARENT="$(git -C "$REPO" rev-parse "$PRODUCTION_COMMIT^")"
"${G[@]}" -C "$WORK/other" checkout -q --detach "$PARENT"
refused "a build before the production commit" "$WORK/r_old_commit.log" bash "$HARNESS" plan "$WORK/p/genotype_c4_production_${PARENT:0:7}" "$WORK/other/build"

rm -rf "$WORK/prod" "$WORK/other"
[ "$fail" -eq 0 ] && echo "production campaign gate: PASS" || { echo "production campaign gate: FAILED"; exit 1; }
