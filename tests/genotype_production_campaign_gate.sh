#!/usr/bin/env bash
# Gate for the prepared-not-launched final C4 campaigns: the plan is derived from the committed
# manifests alone (484 jobs: 224 LZO, 220 LOO, 40 synthetic mosaic, each cohort at 2x, 5x, 10x and 30x),
# the reference donor HG00096 never enters LOO, every call names the committed profile, and launching
# without the reviewed plan's hash is refused. Nothing is simulated, called or evaluated.
#
#   tests/genotype_production_campaign_gate.sh <scratch dir>
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="${1:?usage: genotype_production_campaign_gate.sh <scratch dir>}"
if [ ! -d "$REPO/tests/results/genotype_coverage_c4_held_out_abb35a0" ]; then
  echo "skipped: the frozen held-out read package under tests/results is absent (gitignored)."; exit 0
fi
rm -rf "$WORK"; mkdir -p "$WORK"
fail=0
same() { if [ "$2" = "$3" ]; then echo "  ok: $1 = $2"; else echo "  GATE FAILED: $1 is '$2', expected '$3'"; fail=1; fi; }
key() { awk -F'\t' -v k="$2" '$1==k{print $2; exit}' "$1"; }
bash "$REPO/tests/genotype_production_c4_campaign.sh" plan "$WORK/plan" > "$WORK/plan.log" 2>&1
same "plan exit" "$?" "0"
S="$WORK/plan/settings.tsv"
same "jobs" "$(key "$S" jobs)" "484"
same "LZO jobs" "$(key "$S" jobs_lzo)" "224"
same "LOO jobs" "$(key "$S" jobs_loo)" "220"
same "mosaic jobs" "$(key "$S" jobs_mosaic)" "40"
same "launched" "$(key "$S" launched)" "false"
same "HG00096 LOO jobs" "$(awk -F'\t' '$2=="loo" && $3=="HG00096"' "$WORK/plan/jobs.tsv" | wc -l | tr -d ' ')" "0"
same "coverages per subject" "$(awk -F'\t' 'NR>1{print $4}' "$WORK/plan/jobs.tsv" | sort -un | tr '\n' ' ')" "2 5 10 30 "
# The local campaign profile carries the production profile's science, key for key; only the process
# settings (its name, the footprint and time limits, the thread count) may differ.
science() { grep -v '^#' "$1" | grep -v -E '^(profile|threads|footprint_limit_bytes|time_limit_seconds)'$'\t' | sort; }
same "campaign profile science" "$(diff <(science "$REPO/profiles/genotype_c4.tsv") <(science "$REPO/profiles/genotype_c4_campaign_local.tsv") > /dev/null && echo identical || echo DIFFERS)" "identical"
same "campaign profile footprint limit" "$(awk -F'\t' '$1=="footprint_limit_bytes"{print $2}' "$REPO/profiles/genotype_c4_campaign_local.tsv")" "3758096384"
same "campaign profile leaves threads to the scheduler" "$(grep -c '^threads'$'\t' "$REPO/profiles/genotype_c4_campaign_local.tsv")" "0"
# The v2 profiles are their v1 profiles, key for key, plus exactly the two tie-branching keys.
v2() { grep -v '^#' "$1" | grep -v -E '^(profile|repair_tie_branching|repair_tie_max_branches)'$'\t' | sort; }
for v in genotype_c4 genotype_c4_campaign_local; do
  same "$v v2 = v1 apart from tie branching" "$(diff <(v2 "$REPO/profiles/$v.tsv") <(v2 "$REPO/profiles/${v/genotype_c4/genotype_c4_v2}.tsv") > /dev/null && echo identical || echo DIFFERS)" "identical"
  same "$v v2 enables tie branching" "$(grep -c -E '^(repair_tie_branching'$'\t''1|repair_tie_max_branches'$'\t''16)$' "$REPO/profiles/${v/genotype_c4/genotype_c4_v2}.tsv")" "2"
done
# The v3 profiles are their v2 profiles, key for key, plus exactly representation_guard 1; v1 and v2 never set it.
v3() { grep -v '^#' "$1" | grep -v -E '^(profile|representation_guard)'$'\t' | sort; }
for v in genotype_c4 genotype_c4_campaign_local; do
  v2p="$REPO/profiles/${v/genotype_c4/genotype_c4_v2}.tsv"; v3p="$REPO/profiles/${v/genotype_c4/genotype_c4_v3}.tsv"
  same "$v v3 = v2 apart from the guard" "$(diff <(v3 "$v2p") <(v3 "$v3p") > /dev/null && echo identical || echo DIFFERS)" "identical"
  same "$v v3 binds the guard" "$(grep -c -E '^representation_guard'$'\t''1$' "$v3p")" "1"
  same "$v v3 names itself v3" "$(awk -F'\t' '$1=="profile"{print $2}' "$v3p")" "$(awk -F'\t' '$1=="profile"{print $2}' "$v2p" | sed 's/-v2/-v3/')"
  same "$v v1 and v2 never set the guard" "$(cat "$REPO/profiles/$v.tsv" "$v2p" | grep -c '^representation_guard')" "0"
done
# The v4 profiles are their v3 profiles, key for key, except the name and repair_rounds 11 -> 32.
v4() { grep -v '^#' "$1" | grep -v -E '^(profile|repair_rounds)'$'\t' | sort; }
for v in genotype_c4 genotype_c4_campaign_local; do
  v3p="$REPO/profiles/${v/genotype_c4/genotype_c4_v3}.tsv"; v4p="$REPO/profiles/${v/genotype_c4/genotype_c4_v4}.tsv"
  same "$v v4 = v3 apart from the round limit" "$(diff <(v4 "$v3p") <(v4 "$v4p") > /dev/null && echo identical || echo DIFFERS)" "identical"
  same "$v v3 repair_rounds" "$(awk -F'\t' '$1=="repair_rounds"{print $2}' "$v3p")" "11"
  same "$v v4 repair_rounds" "$(awk -F'\t' '$1=="repair_rounds"{print $2}' "$v4p")" "32"
  same "$v v4 names itself v4" "$(awk -F'\t' '$1=="profile"{print $2}' "$v4p")" "$(awk -F'\t' '$1=="profile"{print $2}' "$v3p" | sed 's/-v3/-v4/')"
done
bash "$REPO/tests/genotype_production_c4_campaign.sh" launch /nonexistent "$WORK/plan" > "$WORK/launch.log" 2>&1
same "launch without approval" "$?" "2"
same "nothing launched" "$(find "$WORK/plan" -name '*.mosaic_status.tsv' | wc -l | tr -d ' ')" "0"
[ "$fail" -eq 0 ] && echo "production campaign gate: PASS" || { echo "production campaign gate: FAILED"; exit 1; }
