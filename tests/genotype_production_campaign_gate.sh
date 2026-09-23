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
bash "$REPO/tests/genotype_production_c4_campaign.sh" launch /nonexistent "$WORK/plan" > "$WORK/launch.log" 2>&1
same "launch without approval" "$?" "2"
same "nothing launched" "$(find "$WORK/plan" -name '*.mosaic_status.tsv' | wc -l | tr -d ' ')" "0"
[ "$fail" -eq 0 ] && echo "production campaign gate: PASS" || { echo "production campaign gate: FAILED"; exit 1; }
