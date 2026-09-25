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

# A mirror of the frozen held-out package, file for file by symbolic link, so a source can be mutated in the mirror
# without ever touching the package.
MIRROR="$WORK/heldout_mirror"
python3 - "$REPO/tests/results/genotype_coverage_c4_held_out_abb35a0" "$MIRROR" <<'PY'
import os, sys
src, dst = sys.argv[1], sys.argv[2]
for root, _, files in os.walk(src):
    target = os.path.join(dst, os.path.relpath(root, src))
    os.makedirs(target, exist_ok=True)
    for name in files:
        os.symlink(os.path.join(root, name), os.path.join(target, name))
PY

echo "== plan"
# A copy of the frozen mosaic signatures (41 small files), so a record can be mutated without touching the frozen one.
SIGS="$WORK/signatures_mirror"
cp -R "$REPO/tests/results/genotype_c4_mosaic_c8fe2a4_raw/signatures" "$SIGS"
chmod -R u+w "$SIGS"   # the frozen records are read-only; only the copy is mutated
GENOTYPE_C4_HELDOUT_READS="$MIRROR" GENOTYPE_C4_MOSAIC_SIGNATURES="$SIGS" \
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
same "estimate labelled historical" "$(head -1 "$PLAN/estimate.tsv" | grep -c 'HISTORICAL, NOT A PREDICTION FOR THE PRODUCTION PROFILE')" "1"
same "estimate provenance recorded" "$(key "$S" estimate_provenance | cut -c1-10)" "HISTORICAL"
# The scientific chain, each step exactly its reviewed differences (the v5/v6 gates define them).
science() { grep -v '^#' "$1" | grep -v '^$' | grep -v -E '^(profile|threads|footprint_limit_bytes|time_limit_seconds)'$'\t' | sort; }
same "v4 campaign-local carries the science of v4" "$(diff <(science "$REPO/profiles/genotype_c4_v4.tsv") <(science "$REPO/profiles/genotype_c4_v4_campaign_local.tsv") > /dev/null && echo identical || echo DIFFERS)" "identical"
bash "$REPO/tests/genotype_v5_profile_gate.sh" > "$WORK/v5.log" 2>&1; same "v5 profile gate" "$?" "0"
bash "$REPO/tests/genotype_v6_profile_gate.sh" > "$WORK/v6.log" 2>&1; same "v6 profile gate" "$?" "0"
same "the production profile sets no seed cache (measured, not adopted)" "$(grep -c '^repair_tie_seed_cache_bytes' "$REPO/profiles/genotype_c4_v6_campaign_local.tsv")" "0"
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

echo "== plan manifest and frozen sources"
same "the plan records its held-out root" "$(key "$S" heldout_reads_root)" "$MIRROR"
same "  as a non-canonical (mirror) root" "$(key "$S" heldout_reads_root_canonical)" "false"
role_count() { awk -F'\t' -v r="$1" 'NR>1 && $4==r' "$PLAN/source_inputs.tsv" | wc -l | tr -d ' '; }
same "source inputs: every unique held-out reads_1, reads_2 and library" \
  "$(awk -F'\t' 'NR>1 && ($4=="reads_1" || $4=="reads_2" || $4=="library")' "$PLAN/source_inputs.tsv" | wc -l | tr -d ' ')" \
  "$(awk -F'\t' 'NR>1 && $13!="harness"{print $13; print $14; print $15}' "$PLAN/jobs.tsv" | sort -u | wc -l | tr -d ' ')"
same "  56 donors x 3 files" "$(($(role_count reads_1) + $(role_count reads_2) + $(role_count library)))" "168"
for r in graph_gz donor_manifest eligibility_manifest mosaic_manifest mosaic_signature_index; do
  same "  one $r" "$(role_count $r)" "1"
done
same "  graph gzip at its repository path" "$(awk -F'\t' '$4=="graph_gz"{print $1}' "$PLAN/source_inputs.tsv")" "$REPO/tests/real_data/c4.gfa.gz"
same "  40 mosaic signatures" "$(role_count mosaic_signature)" "40"
same "  keyed by exactly the planned mosaic (plan, coverage)" \
  "$(awk -F'\t' '$4=="mosaic_signature"{print $5}' "$PLAN/source_inputs.tsv" | sort | tr '\n' ' ')" \
  "$(awk -F'\t' 'NR>1 && $2=="mosaic"{print $3":"$4}' "$PLAN/jobs.tsv" | sort | tr '\n' ' ')"
same "  213 frozen inputs in all" "$(($(wc -l < "$PLAN/source_inputs.tsv") - 1))" "213"
same "  the mirror roots are recorded" "$(key "$S" mosaic_signatures_root)|$(key "$S" mosaic_signatures_root_canonical)" "$SIGS|false"
same "  every source path absolute" "$(awk -F'\t' 'NR>1 && substr($1,1,1)!="/"' "$PLAN/source_inputs.tsv" | wc -l | tr -d ' ')" "0"
same "  hashes and sizes as the files are now" "$(awk -F'\t' 'NR>1{print $1"\t"$2"\t"$3}' "$PLAN/source_inputs.tsv" | while IFS=$'\t' read -r f h b; do [ "$(sha "$f")" = "$h" ] && [ "$(wc -c < "$f" | tr -d ' ')" = "$b" ] || echo bad; done | wc -l | tr -d ' ')" "0"
same "manifest components" "$(awk -F'\t' 'NR>1{print $1}' "$PLAN/plan_manifest.tsv" | tr '\n' ' ')" "jobs.tsv settings.tsv source_inputs.tsv "
for c in jobs.tsv settings.tsv source_inputs.tsv; do
  same "  $c hash in the manifest" "$(key "$PLAN/plan_manifest.tsv" "$c")" "$(sha "$PLAN/$c")"
done
same "the manifest is not named by the settings (non-circular)" "$(grep -c plan_manifest "$S")" "0"
APPROVED="$(sha "$PLAN/plan_manifest.tsv")"
JOBS_SHA="$(sha "$PLAN/jobs.tsv")"
staged() { find "$PLAN" \( -name '*.mosaic_status.tsv' -o -name 'c4.gfa' -o -name 'inputs.tsv' -o -name 'launch.tsv' -o -name '*.command.txt' \) | wc -l | tr -d ' '; }
launch_as() {   # launch_as <log> <approval> [build]
  GENOTYPE_C4_CAMPAIGN_APPROVED="$2" bash "$HARNESS" launch "${3:-$WORK/prod/build}" "$PLAN" > "$WORK/$1" 2>&1
}
refuses() {   # refuses <what> <log> <expected message>
  local code=$?
  same "$1 refused" "$code" "2"
  same "  for the intended reason" "$(grep -c -F "$3" "$WORK/$2")" "1"
  same "  before anything is staged or called" "$(staged)" "0"
}
reseal() {   # the plan manifest recomputed over the current components: a self-consistent tampered plan
  { printf "component\tsha256\tbytes\n"; for c in jobs.tsv settings.tsv source_inputs.tsv; do
      printf "%s\t%s\t%s\n" "$c" "$(sha "$PLAN/$c")" "$(wc -c < "$PLAN/$c" | tr -d ' ')"; done; } > "$PLAN/plan_manifest.tsv"
  sha "$PLAN/plan_manifest.tsv"
}
keep() { cp "$PLAN/$1" "$WORK/reviewed.$1"; }
restore() { cp "$WORK/reviewed.$1" "$PLAN/$1"; }
for c in jobs.tsv settings.tsv source_inputs.tsv plan_manifest.tsv; do keep "$c"; done

echo "== launch refusals"
bash "$HARNESS" launch "$WORK/prod/build" "$PLAN" > "$WORK/launch.log" 2>&1
refuses "launch without approval" launch.log "requires review of"
launch_as launch_wrong.log 0000
refuses "launch with another plan's approval" launch_wrong.log "requires review of"
launch_as launch_jobs_token.log "$JOBS_SHA"
refuses "launch approved by the jobs.tsv hash alone" launch_jobs_token.log "requires review of"
# settings.tsv changed, jobs.tsv untouched: refused by the manifest.
printf "tampered\ttrue\n" >> "$PLAN/settings.tsv"
launch_as launch_settings.log "$APPROVED"
refuses "launch after settings.tsv changed (jobs.tsv unchanged)" launch_settings.log "settings.tsv is not the component the plan manifest reviewed"
restore settings.tsv
# source_inputs.tsv changed: refused by the manifest.
sed -i.bak -e '2s/\t[0-9a-f]\{64\}\t/\t0000000000000000000000000000000000000000000000000000000000000000\t/' "$PLAN/source_inputs.tsv" && rm -f "$PLAN/source_inputs.tsv.bak"
launch_as launch_sources_file.log "$APPROVED"
refuses "launch after source_inputs.tsv changed" launch_sources_file.log "source_inputs.tsv is not the component the plan manifest reviewed"
restore source_inputs.tsv
# One held-out source FASTQ changed (in the mirror only).
FQ="$(awk -F'\t' 'NR>1 && $4=="reads_1"{print $1; exit}' "$PLAN/source_inputs.tsv")"
rm "$FQ"; { gzip -dc "$(readlink -f "$REPO/tests/results/genotype_coverage_c4_held_out_abb35a0/${FQ#$MIRROR/}")"; printf '@x\nA\n+\nI\n'; } | gzip -c > "$FQ"
launch_as launch_fastq.log "$APPROVED"
refuses "launch after a source FASTQ changed" launch_fastq.log "the frozen source $FQ is not the one the plan bound"
rm "$FQ"; ln -s "$REPO/tests/results/genotype_coverage_c4_held_out_abb35a0/${FQ#$MIRROR/}" "$FQ"
# One held-out source library changed (in the mirror only).
LIB="$(awk -F'\t' 'NR>1 && $4=="library"{print $1; exit}' "$PLAN/source_inputs.tsv")"
rm "$LIB"; { cat "$REPO/tests/results/genotype_coverage_c4_held_out_abb35a0/${LIB#$MIRROR/}"; echo "# changed"; } > "$LIB"
launch_as launch_library.log "$APPROVED"
refuses "launch after a source library changed" launch_library.log "the frozen source $LIB is not the one the plan bound"
rm "$LIB"; ln -s "$REPO/tests/results/genotype_coverage_c4_held_out_abb35a0/${LIB#$MIRROR/}" "$LIB"
# A binary that is not the plan's.
fake_build "$WORK/other" "$PRODUCTION_COMMIT" other
launch_as launch_binary.log "$APPROVED" "$WORK/other/build"
refuses "launch with a binary that is not the plan's" launch_binary.log "production_binary_sha256"
# A campaign file changed after review: a self-consistent plan whose settings name another runner hash.
sed -i.bak -e 's/^\(scripts\/c4_production_campaign.py_sha256\t\).*/\1changed/' "$PLAN/settings.tsv" && rm -f "$PLAN/settings.tsv.bak"
RESEALED="$(reseal)"
launch_as launch_file.log "$RESEALED"
refuses "launch after a campaign file changed" launch_file.log "scripts/c4_production_campaign.py is not the file the plan reviewed"
for c in settings.tsv plan_manifest.tsv; do restore "$c"; done
same "the reviewed manifest is restored" "$(sha "$PLAN/plan_manifest.tsv")" "$APPROVED"
echo "== frozen repository inputs and mosaic records"
# A self-consistent tampered plan: one source_inputs.tsv row edited by awk, then the manifest resealed.
tamper() {   # tamper <awk program over source_inputs.tsv>
  awk -F'\t' -v OFS='\t' "$1" "$WORK/reviewed.source_inputs.tsv" > "$PLAN/source_inputs.tsv"; reseal
}
untamper() { restore source_inputs.tsv; restore plan_manifest.tsv; }
# The tracked repository inputs are never modified: their bound record is changed instead, which a launch must
# refuse exactly as it would refuse the file changing.
for r in graph_gz donor_manifest eligibility_manifest mosaic_manifest; do
  F="$(awk -F'\t' -v r=$r '$4==r{print $1}' "$WORK/reviewed.source_inputs.tsv")"
  T="$(tamper "\$4==\"$r\"{\$2=\"0000000000000000000000000000000000000000000000000000000000000000\"} {print}")"
  launch_as "launch_$r.log" "$T"
  refuses "launch with the $r not as bound" "launch_$r.log" "the frozen source $F is not the one the plan bound"
  untamper
done
# A repository input bound at another path (a copy, correct hash) is not the input.
cp "$REPO/tests/genotype_development_donors.tsv" "$WORK/donors_copy.tsv"
T="$(tamper "\$4==\"donor_manifest\"{\$1=\"$WORK/donors_copy.tsv\"} {print}")"
launch_as launch_donor_path.log "$T"
refuses "launch with the donor manifest bound at another path" launch_donor_path.log "must bind $REPO/tests/genotype_development_donors.tsv once as donor_manifest"
untamper
# The signature index and one signature changed, and one signature missing (in the copy only).
cp "$SIGS/index.tsv" "$WORK/index.reviewed"; printf '# changed\n' >> "$SIGS/index.tsv"
same "  (the index mutation took effect)" "$(cmp -s "$SIGS/index.tsv" "$WORK/index.reviewed" && echo no || echo yes)" "yes"
launch_as launch_index.log "$APPROVED"
refuses "launch after the signature index changed" launch_index.log "the frozen source $SIGS/index.tsv is not the one the plan bound"
cp "$WORK/index.reviewed" "$SIGS/index.tsv"
SIG="$(awk -F'\t' '$4=="mosaic_signature"{print $1; exit}' "$PLAN/source_inputs.tsv")"
cp "$SIG" "$WORK/sig.reviewed"; printf 'extra\tline\n' >> "$SIG"
same "  (the signature mutation took effect)" "$(cmp -s "$SIG" "$WORK/sig.reviewed" && echo no || echo yes)" "yes"
launch_as launch_signature.log "$APPROVED"
refuses "launch after one signature changed" launch_signature.log "the frozen source $SIG is not the one the plan bound"
rm "$SIG"
launch_as launch_signature_missing.log "$APPROVED"
refuses "launch with one signature missing" launch_signature_missing.log "the frozen source $SIG is not the one the plan bound"
cp "$WORK/sig.reviewed" "$SIG"
# A duplicate and a missing (plan, coverage) mapping.
K1="$(awk -F'\t' '$4=="mosaic_signature"{print $5; exit}' "$WORK/reviewed.source_inputs.tsv")"
K2="$(awk -F'\t' '$4=="mosaic_signature"{n++; if (n==2) {print $5; exit}}' "$WORK/reviewed.source_inputs.tsv")"
T="$(tamper "\$4==\"mosaic_signature\" && \$5==\"$K2\"{\$5=\"$K1\"} {print}")"
launch_as launch_duplicate.log "$T"
refuses "launch with a duplicate (plan, coverage)" launch_duplicate.log "binds (plan, coverage) ['$K1'] more than once"
untamper
T="$(tamper "!(\$4==\"mosaic_signature\" && \$5==\"$K2\")")"
launch_as launch_missing_key.log "$T"
refuses "launch with a (plan, coverage) unmapped" launch_missing_key.log "do not match the planned mosaic jobs exactly"
untamper
same "the reviewed plan is intact after the mutations" "$(sha "$PLAN/plan_manifest.tsv")|$(sha "$PLAN/source_inputs.tsv")" "$APPROVED|$(sha "$WORK/reviewed.source_inputs.tsv")"
# The comparison after regeneration, and an unrelated file beside the bound records.
python3 - "$REPO" "$PLAN" "$SIGS" <<'PY' || fail=1
import importlib.util, pathlib, sys, shutil
repo, plan, sigs = map(pathlib.Path, sys.argv[1:4])
spec = importlib.util.spec_from_file_location("runner", repo / "scripts" / "c4_production_campaign.py")
runner = importlib.util.module_from_spec(spec); spec.loader.exec_module(runner)
bad = []
def same(what, got, want):
    print(("  ok: " if got == want else "  GATE FAILED: ") + f"{what} = {got!r}" + ("" if got == want else f", expected {want!r}"))
    if got != want: bad.append(what)
def refused(what, call, text):
    try:
        call(); same(what, "accepted", "refused")
    except runner.Defect as defect:
        same(what, "refused" if text in str(defect) else f"refused for another reason: {defect}", "refused")
jobs = runner.read_jobs(plan)
counts, records = runner.verify_sources(plan, jobs)
same("40 frozen mosaic records", len(records), 40)
same("  exactly the planned mosaic jobs", sorted(records), sorted((j["subject"], j["coverage"]) for j in jobs if j["campaign"] == "mosaic"))
# An unrelated record in the signatures directory, even one claiming a planned key, never overrides a bound one.
k = sorted(records)[0]
fake = dict(records[k]); fake["reads_1"] = "0" * 64
(sigs / ("f" * 64 + ".tsv")).write_text("key\tvalue\n" + "".join(f"{a}\t{b}\n" for a, b in fake.items()))
_, again = runner.verify_sources(plan, jobs)
same("an unrelated signature file beside the bound set is ignored", again[k]["reads_1"], records[k]["reads_1"])
(sigs / ("f" * 64 + ".tsv")).unlink()
regen = {key: {n: r[n] for n in runner.MOSAIC_COMPARED} for key, r in records.items()}
same("every regenerated cell identical to its record", set(runner.compare_mosaic(regen, records).values()), {"identical"})
extra = dict(regen); extra[("dense_switching_1", "15")] = regen[k]
refused("a regenerated cell with no frozen record", lambda: runner.compare_mosaic(extra, records), "no frozen record")
short = dict(regen); del short[k]
refused("a frozen record with no regenerated cell", lambda: runner.compare_mosaic(short, records), "no regenerated cell")
changed = {key: dict(v) for key, v in regen.items()}; changed[k]["library"] = "0" * 64
same("a regenerated cell that differs is DIFFERS", runner.compare_mosaic(changed, records)[k], "DIFFERS")
# A resume over a prepared.tsv that is not the bound set, all identical, is refused before anything runs.
out = plan.parent / "resume_probe"; (out / "mosaic_inputs").mkdir(parents=True)
rows = ["plan\tcoverage\treads_1_sha256\treads_2_sha256\tlibrary_sha256\tstaged_graph_sha256\tvs_c8fe2a4_record"]
rows += [f"{p}\t{c}\t-\t-\t-\t-\t{'no_record' if (p, c) == k else 'identical'}" for p, c in sorted(records)]
(out / "mosaic_inputs" / "prepared.tsv").write_text("\n".join(rows) + "\n")
refused("a resume whose prepared mosaic set holds a no_record cell", lambda: runner.prepare_mosaic(out, out, None, out / "log", records), "every cell identical")
shutil.rmtree(out)
sys.exit(1 if bad else 0)
PY

# A resume whose previous launch ran another profile: every plan and source check passes, then the resume refuses.
printf "key\tvalue\npanvar_sha256\t%s\nprofile\tprofiles/genotype_c4_v6_campaign_local.tsv\nprofile_sha256\tstale\njobs\t1\nthreads\t8\njobs_sha256\t%s\n" \
  "$(sha "$WORK/prod/build/panvar")" "$JOBS_SHA" > "$PLAN/launch.tsv"
launch_as launch_resume.log "$APPROVED"
same "resume after a profile change" "$?" "2"
same "  every source verified first" "$(grep -c "213 frozen source inputs verified" "$WORK/launch_resume.log")" "1"
same "  refused as a resume" "$(grep -c "resume refused: profile_sha256" "$WORK/launch_resume.log")" "1"
rm -f "$PLAN/launch.tsv"
same "nothing launched or staged" "$(staged)" "0"

echo "== worker failures (the runner's dispatcher)"
python3 - "$REPO" "$WORK" <<'PY' || fail=1
import importlib.util, pathlib, sys, os
repo, work = map(pathlib.Path, sys.argv[1:3])
spec = importlib.util.spec_from_file_location("runner", repo / "scripts" / "c4_production_campaign.py")
runner = importlib.util.module_from_spec(spec); spec.loader.exec_module(runner)
bad = []
def same(what, got, want):
    print(("  ok: " if got == want else "  GATE FAILED: ") + f"{what} = {got!r}" + ("" if got == want else f", expected {want!r}"))
    if got != want: bad.append(what)
root = work / "workers"; root.mkdir()
# A real process that exits 0 without writing a status: an external failure, never a finished job.
build = root / "build"; build.mkdir()
(build / "panvar").write_text("#!/bin/sh\nexit 0\n"); (build / "panvar").chmod(0o755)
job = {"job": "loo.X.10", "campaign": "loo", "subject": "X", "coverage": "10", "model_gfa": str(root / "m.gfa"),
       "reads1": str(root / "r1"), "reads2": str(root / "r2"), "library": str(root / "lib"), "curation": "none",
       "truth_gfa": str(root / "t.gfa"), "out_prefix": str(root / "loo" / "X" / "call.10")}
outcome, detail = runner.run_job(job, build, build / "genotype_probe", 8)
same("a process exiting without a status is an external failure", outcome, "failed")
same("  and says so", "without a status" in detail, True)
def jobs(n): return [dict(job, job=f"loo.J{i}.10", out_prefix=str(root / f"J{i}" / "call.10")) for i in range(n)]
def ledger_for(js, name):
    d = root / name; d.mkdir(); return runner.Ledger(d, js)
quiet = lambda message: None
# Retry after a missing status succeeds.
js = jobs(1); led = ledger_for(js, "retry"); calls = []
def flaky(j):
    calls.append(j["job"])
    return ("failed", "the process exited (0) without a status") if len(calls) == 1 else ("strict_call", "PRIOR_ONLY")
reasons = runner.dispatch(js, led, 1, flaky, log=quiet)
same("retry after a missing status: no stop", reasons, [])
same("  finished as a call on attempt 2", (led.rows[js[0]["job"]]["outcome"], led.rows[js[0]["job"]]["attempts"]), ("strict_call", "2"))
# A repeated missing status stops the campaign, and nothing after it runs.
js = jobs(3); led = ledger_for(js, "repeat"); ran = []
def missing(j):
    ran.append(j["job"]); return ("failed", "the process exited (0) without a status")
reasons = runner.dispatch(js, led, 1, missing, log=quiet)
same("repeated missing status stops the campaign", reasons, [f"{js[0]['job']}: repeated external failure"])
same("  after exactly MAX_ATTEMPTS attempts and no other job", ran, [js[0]["job"]] * runner.MAX_ATTEMPTS)
same("  the job is recorded failed, not finished", (led.rows[js[0]["job"]]["outcome"], led.finished(js[0]["job"])), ("failed", False))
# A worker that raises: recorded, stop, nonzero -- and dispatch does not carry on.
for workers in (1, 2):
    js = jobs(4); led = ledger_for(js, f"raise{workers}"); ran = []
    def explode(j):
        ran.append(j["job"])
        if j["job"] == js[0]["job"]: raise RuntimeError("scripted worker bug")
        return ("strict_call", "PRIOR_ONLY")
    reasons = runner.dispatch(js, led, workers, explode, log=quiet)
    same(f"a raising worker stops the campaign ({workers} worker(s))", any("worker raised RuntimeError: scripted worker bug" in r for r in reasons), True)
    same("  and is recorded as a defect", led.rows[js[0]["job"]]["outcome"], "defect")
    same("  and dispatch does not continue past it", len(ran) <= workers, True)
# An unknown outcome label fails closed.
js = jobs(2); led = ledger_for(js, "unknown")
reasons = runner.dispatch(js, led, 1, lambda j: ("half_call", "?"), log=quiet)
same("an unknown outcome stops the campaign", len(reasons) == 1 and "unknown outcome 'half_call'" in reasons[0], True)
same("  never a finished job", led.finished(js[0]["job"]), False)
# `pending` returned by an attempt is not a finished job either.
js = jobs(1); led = ledger_for(js, "pending")
reasons = runner.dispatch(js, led, 1, lambda j: ("pending", "no status"), log=quiet)
same("a pending attempt fails closed", len(reasons) == 1 and led.finished(js[0]["job"]) is False, True)
# Input preparation that raises stops before the job is submitted.
js = jobs(2); led = ledger_for(js, "prepare"); ran = []
def prepare(j): raise runner.Defect("input moved")
reasons = runner.dispatch(js, led, 1, lambda j: (ran.append(j) or ("strict_call", "x")), prepare=prepare, log=quiet)
same("a failed preparation stops before any call", (len(reasons), ran), (1, []))
# Every finished job is finished; a clean run returns no reasons.
js = jobs(3); led = ledger_for(js, "clean")
reasons = runner.dispatch(js, led, 2, lambda j: ("no_call", "tie_branching_incomplete"), log=quiet)
same("official no-calls finish their jobs and never stop", (reasons, all(led.finished(j["job"]) for j in js)), ([], True))
sys.exit(1 if bad else 0)
PY

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
# A signature index that names a (plan, coverage) twice or not at all, or points at a missing record: the plan refuses.
bad_sigs() { rm -rf "$WORK/bad_sigs"; cp -R "$REPO/tests/results/genotype_c4_mosaic_c8fe2a4_raw/signatures" "$WORK/bad_sigs"; chmod -R u+w "$WORK/bad_sigs"; }
bad_sigs; head -1 "$WORK/bad_sigs/index.tsv" >> "$WORK/bad_sigs/index.tsv"
refused "a plan whose signature index names a key twice" "$WORK/r_sig_dup.log" env GENOTYPE_C4_MOSAIC_SIGNATURES="$WORK/bad_sigs" bash "$HARNESS" plan "$WORK/sd/genotype_c4_production_$P7" "$WORK/prod/build"
bad_sigs; sed -i.bak -e '1d' "$WORK/bad_sigs/index.tsv"
refused "a plan whose signature index misses a planned key" "$WORK/r_sig_missing.log" env GENOTYPE_C4_MOSAIC_SIGNATURES="$WORK/bad_sigs" bash "$HARNESS" plan "$WORK/sm/genotype_c4_production_$P7" "$WORK/prod/build"
bad_sigs; rm "$WORK/bad_sigs/$(awk -F'\t' 'NR==1{print $3}' "$WORK/bad_sigs/index.tsv").tsv"
refused "a plan whose indexed signature is missing" "$WORK/r_sig_file.log" env GENOTYPE_C4_MOSAIC_SIGNATURES="$WORK/bad_sigs" bash "$HARNESS" plan "$WORK/sf/genotype_c4_production_$P7" "$WORK/prod/build"
PARENT="$(git -C "$REPO" rev-parse "$PRODUCTION_COMMIT^")"
"${G[@]}" -C "$WORK/other" checkout -q --detach "$PARENT"
refused "a build before the production commit" "$WORK/r_old_commit.log" bash "$HARNESS" plan "$WORK/p/genotype_c4_production_${PARENT:0:7}" "$WORK/other/build"

rm -rf "$WORK/prod" "$WORK/other"
[ "$fail" -eq 0 ] && echo "production campaign gate: PASS" || { echo "production campaign gate: FAILED"; exit 1; }
