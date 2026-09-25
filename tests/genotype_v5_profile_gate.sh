#!/usr/bin/env bash
# The v5 campaign-local profile is its v4 profile, key for key and value for value, except exactly the reviewed
# speed and memory settings (and its name); v4 sets none of the new keys. Nothing is run.
#
#   tests/genotype_v5_profile_gate.sh
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
V4="$REPO/profiles/genotype_c4_v4_campaign_local.tsv"
V5="$REPO/profiles/genotype_c4_v5_campaign_local.tsv"
fail=0
same() { if [ "$2" = "$3" ]; then echo "  ok: $1"; else echo "  GATE FAILED: $1 ('$2', expected '$3')"; fail=1; fi; }
settings() { grep -v '^#' "$1" | grep -v '^$' | sort; }
INTENDED='^(profile|footprint_limit_bytes|union_placement_cache_bytes|max_union_cached_placements_consumed|union_origin_cache_bytes|max_union_cached_origins_consumed|max_union_cached_states_consumed|union_prefilter_cache_bytes|max_union_cached_prefilters_consumed|start_sort_scratch_bytes|emission_lookup_timing)'$'\t'
same "v5 = v4 apart from the intended keys" \
  "$(diff <(settings "$V4" | grep -v -E "$INTENDED") <(settings "$V5" | grep -v -E "$INTENDED") > /dev/null && echo identical || echo DIFFERS)" "identical"
same "every v5 key appears once" "$(settings "$V5" | cut -f1 | sort | uniq -d | wc -l | tr -d ' ')" "0"
expect() { same "v5 $1 = $2" "$(awk -F'\t' -v k="$1" '$1==k{print $2}' "$V5")" "$2"; }
expect profile c4-panel-anchored-v5-campaign-local
expect footprint_limit_bytes 4294967296
expect union_placement_cache_bytes 536870912
expect max_union_cached_placements_consumed 2000000000
expect union_origin_cache_bytes 1073741824
expect max_union_cached_origins_consumed 2000000000
expect max_union_cached_states_consumed 2000000000
expect union_prefilter_cache_bytes 134217728
expect max_union_cached_prefilters_consumed 2000000000
expect start_sort_scratch_bytes 335544320
expect emission_lookup_timing 0
# Unchanged from v4, named because the review named them.
for k in repair_rounds representation_guard repair_tie_branching repair_tie_max_branches max_physical_origins \
         max_state_classes emission_cache_bytes time_limit_seconds; do
  same "v5 $k = v4 $k" "$(awk -F'\t' -v k="$k" '$1==k{print $2}' "$V5")" "$(awk -F'\t' -v k="$k" '$1==k{print $2}' "$V4")"
done
expect repair_rounds 32
expect representation_guard 1
same "v4 sets none of the new keys" \
  "$(settings "$V4" | grep -c -E '^(union_placement_cache_bytes|max_union_cached_placements_consumed|union_origin_cache_bytes|max_union_cached_origins_consumed|max_union_cached_states_consumed|union_prefilter_cache_bytes|max_union_cached_prefilters_consumed|start_sort_scratch_bytes|emission_lookup_timing)'$'\t')" "0"
same "v4 footprint limit" "$(awk -F'\t' '$1=="footprint_limit_bytes"{print $2}' "$V4")" "3758096384"
same "the v5 profile leaves threads to the scheduler" "$(grep -c '^threads'$'\t' "$V5")" "0"
# The comparison can fail: a copy of v5 with one more changed setting must be caught.
MUTANT="$(mktemp)"
sed -e 's/^max_physical_origins\t700000000$/max_physical_origins\t700000001/' "$V5" > "$MUTANT"
same "a v5 copy with another changed setting is caught" \
  "$(diff <(settings "$V4" | grep -v -E "$INTENDED") <(settings "$MUTANT" | grep -v -E "$INTENDED") > /dev/null && echo identical || echo DIFFERS)" "DIFFERS"
rm -f "$MUTANT"
[ "$fail" -eq 0 ] && echo "v5 profile gate: PASS" || { echo "v5 profile gate: FAILED"; exit 1; }
