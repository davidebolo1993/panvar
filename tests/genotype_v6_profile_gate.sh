#!/usr/bin/env bash
# The v6 campaign-local profile is its v5 profile, key for key and value for value, except its name and
# repair_origin_cache_bytes (384 MiB -> 768 MiB). Nothing is run.
#
#   tests/genotype_v6_profile_gate.sh
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
V5="$REPO/profiles/genotype_c4_v5_campaign_local.tsv"
V6="$REPO/profiles/genotype_c4_v6_campaign_local.tsv"
fail=0
same() { if [ "$2" = "$3" ]; then echo "  ok: $1"; else echo "  GATE FAILED: $1 ('$2', expected '$3')"; fail=1; fi; }
settings() { grep -v '^#' "$1" | grep -v '^$' | sort; }
INTENDED='^(profile|repair_origin_cache_bytes)'$'\t'
same "v6 = v5 apart from its name and the repair cache budget" \
  "$(diff <(settings "$V5" | grep -v -E "$INTENDED") <(settings "$V6" | grep -v -E "$INTENDED") > /dev/null && echo identical || echo DIFFERS)" "identical"
same "every v6 key appears once" "$(settings "$V6" | cut -f1 | sort | uniq -d | wc -l | tr -d ' ')" "0"
same "v6 profile" "$(awk -F'\t' '$1=="profile"{print $2}' "$V6")" "c4-panel-anchored-v6-campaign-local"
same "v6 repair_origin_cache_bytes" "$(awk -F'\t' '$1=="repair_origin_cache_bytes"{print $2}' "$V6")" "805306368"
same "v5 repair_origin_cache_bytes (frozen)" "$(awk -F'\t' '$1=="repair_origin_cache_bytes"{print $2}' "$V5")" "402653184"
same "v6 footprint limit = v5's" "$(awk -F'\t' '$1=="footprint_limit_bytes"{print $2}' "$V6")" "4294967296"
MUTANT="$(mktemp)"
sed -e 's/^repair_rounds\t32$/repair_rounds\t33/' "$V6" > "$MUTANT"
same "a v6 copy with another changed setting is caught" \
  "$(diff <(settings "$V5" | grep -v -E "$INTENDED") <(settings "$MUTANT" | grep -v -E "$INTENDED") > /dev/null && echo identical || echo DIFFERS)" "DIFFERS"
rm -f "$MUTANT"
[ "$fail" -eq 0 ] && echo "v6 profile gate: PASS" || { echo "v6 profile gate: FAILED"; exit 1; }
