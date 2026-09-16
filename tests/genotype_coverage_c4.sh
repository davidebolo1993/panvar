#!/usr/bin/env bash
# Developer harness, not a regression test: the C4 nested-coverage leave-zero-out experiment (review round 8).
#
#   BUILD=<clean build dir> ROLE=development|held_out tests/genotype_coverage_c4.sh <out-dir>
#
# Settings come only from tests/genotype_coverage_settings.tsv, frozen before any read of this experiment was
# simulated. Each donor is simulated ONCE at the source coverage with the frozen seed; every lower coverage is a
# nested subset of that simulation (genotype_probe subsample). At each coverage both arms -- panel (K = 0) and
# mosaic (K = 8) -- are ranked from one proposal and one scoring pass, with the frozen phase audit, and each is
# evaluated against the donor's truth. Writes one row per donor, arm, coverage and block, plus a run-level table.
# Development donors and held-out donors are run separately through ROLE; a role never writes another's directory.
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${BUILD:-$REPO/build}"
OUT="${1:?usage: BUILD=<build> ROLE=<development|held_out> tests/genotype_coverage_c4.sh <out-dir>}"
ROLE="${ROLE:?ROLE must be development or held_out}"
[ "$ROLE" = development ] || [ "$ROLE" = held_out ] || { echo "ROLE must be development or held_out" >&2; exit 1; }
MANIFEST="$REPO/tests/genotype_development_donors.tsv"
FROZEN="$REPO/tests/genotype_coverage_settings.tsv"
setting() { awk -F'\t' -v key="$1" '$1==key{print $2}' "$FROZEN"; }
SEED="$(setting seed)"
SOURCE_COVERAGE="$(setting source_coverage)"
COVERAGES="$(setting coverages | tr ',' ' ')"
STATE_LIMIT="$(setting state_limit)"
PANEL_K="$(setting panel_arm_mosaic_candidates)"
MOSAIC_K="$(setting mosaic_arm_mosaic_candidates)"
AUDIT="$(setting phase_audit_candidates)"
READ_LENGTH="$(setting read_length)"
[ -n "$SEED" ] && [ -n "$SOURCE_COVERAGE" ] && [ -n "$COVERAGES" ] && [ -n "$STATE_LIMIT" ] && [ "$PANEL_K" = 0 ] &&
  [ -n "$MOSAIC_K" ] && [ -n "$AUDIT" ] || { echo "incomplete $FROZEN" >&2; exit 1; }
DONORS="${DONORS:-$(awk -F'\t' -v role="$ROLE" '$2==role{print $1}' "$MANIFEST")}"
JOBS="${JOBS:-4}"
REF="HG00096#1#haplotype1-0000024:31848049-32074186"
PROBE="$BUILD/genotype_probe"

mkdir -p "$OUT"
[ -z "$(ls -A "$OUT" 2>/dev/null)" ] || { echo "$OUT is not empty; this experiment writes a new directory" >&2; exit 1; }
gzip -dc "$REPO/tests/real_data/c4.gfa.gz" > "$OUT/c4.gfa"

if /usr/bin/time -l true >/dev/null 2>&1; then TIME_FLAG="-l"; else TIME_FLAG="-v"; fi
rss_of() {
  if [ "$TIME_FLAG" = "-l" ]; then awk '/maximum resident set size/{print $1}' "$1"
  else awk -F: '/Maximum resident set size/{gsub(/ /,"",$2); print $2*1024}' "$1"; fi
}
wall_of() {
  if [ "$TIME_FLAG" = "-l" ]; then awk '/ real /{print $1; exit}' "$1"
  else awk -F': ' '/Elapsed \(wall clock\)/{n=split($2,p,":"); s=0; for(i=1;i<=n;i++) s=s*60+p[i]; print s}' "$1"; fi
}

{
  printf 'key\tvalue\n'
  printf 'git_commit\t%s\n' "$(git -C "$REPO" rev-parse HEAD)"
  printf 'git_tree_dirty\t%s\n' "$(git -C "$REPO" status --porcelain | grep -q . && echo true || echo false)"
  printf 'probe_sha256\t%s\n' "$(shasum -a 256 "$PROBE" | cut -d' ' -f1)"
  printf 'script_sha256\t%s\n' "$(shasum -a 256 "${BASH_SOURCE[0]}" | cut -d' ' -f1)"
  printf 'settings_sha256\t%s\n' "$(shasum -a 256 "$FROZEN" | cut -d' ' -f1)"
  printf 'plot_script_sha256\t%s\n' "$(shasum -a 256 "$REPO/scripts/plot_genotype_coverage.py" | cut -d' ' -f1)"
  printf 'manifest_sha256\t%s\n' "$(shasum -a 256 "$MANIFEST" | cut -d' ' -f1)"
  printf 'graph_sha256\t%s\n' "$(shasum -a 256 "$OUT/c4.gfa" | cut -d' ' -f1)"
  printf 'role\t%s\nseed\t%s\nsource_coverage\t%s\ncoverages\t%s\nstate_limit\t%s\npanel_arm_K\t%s\nmosaic_arm_K\t%s\nphase_audit\t%s\nreference\t%s\ndonors\t%s\n' \
    "$ROLE" "$SEED" "$SOURCE_COVERAGE" "$COVERAGES" "$STATE_LIMIT" "$PANEL_K" "$MOSAIC_K" "$AUDIT" "$REF" "$(echo $DONORS)"
} > "$OUT/signature.tsv"

run_donor() {
  local donor="$1" dir="$OUT/$1"
  mkdir -p "$dir/source"
  "$PROBE" simulate --gfa "$OUT/c4.gfa" --reference "$REF" --donor "$donor" --depth "$SOURCE_COVERAGE" \
    --read-length "$READ_LENGTH" --seed "$SEED" --out "$dir/source" > "$dir/source.log" 2>&1
  for coverage in $COVERAGES; do
    local reads="$dir/cov$coverage"
    mkdir -p "$reads"
    "$PROBE" subsample --reads1 "$dir/source/reads_1.fq" --reads2 "$dir/source/reads_2.fq" \
      --library "$dir/source/library.txt" --source-coverage "$SOURCE_COVERAGE" --coverage "$coverage" \
      --out "$reads" >> "$dir/source.log" 2>&1
    local prefix="$reads/lzo"
    /usr/bin/time "$TIME_FLAG" "$PROBE" genotype --gfa "$OUT/c4.gfa" --reference "$REF" \
      --reads1 "$reads/reads_1.fq" --reads2 "$reads/reads_2.fq" --library "$reads/library.txt" \
      --out-prefix "$prefix" --state-limit "$STATE_LIMIT" --prefixes "$PANEL_K,$MOSAIC_K" \
      --phase-audit "$AUDIT" > "$prefix.log" 2> "$prefix.time"
    printf 'peak_rss_bytes\t%s\nwall_seconds\t%s\n' "$(rss_of "$prefix.time")" "$(wall_of "$prefix.time")" \
      > "$prefix.resources.tsv"
    for k in "$PANEL_K" "$MOSAIC_K"; do
      if [ -f "$prefix.K$k.mosaic_status.tsv" ]; then
        "$PROBE" evaluate --scheme lzo --donor "$donor" --reference "$REF" --run-prefix "$prefix.K$k" \
          --model-gfa "$OUT/c4.gfa" --truth-gfa "$OUT/c4.gfa" --out "$prefix.K$k.eval" >> "$prefix.log" 2>&1
      fi
    done
  done
  echo "done $donor"
}
export -f run_donor rss_of wall_of
export OUT PROBE REF SEED SOURCE_COVERAGE COVERAGES STATE_LIMIT PANEL_K MOSAIC_K AUDIT READ_LENGTH TIME_FLAG
printf '%s\n' $DONORS | xargs -P "$JOBS" -I{} bash -c 'run_donor "$@"' _ {}

# ---- collect: one row per donor, arm, coverage and block; and one per run --------------------------------
OUT="$OUT" DONORS="$DONORS" COVERAGES="$COVERAGES" PANEL_K="$PANEL_K" MOSAIC_K="$MOSAIC_K" ROLE="$ROLE" python3 - <<'COLLECT'
import csv, os
out = os.environ["OUT"]; donors = os.environ["DONORS"].split(); coverages = os.environ["COVERAGES"].split()
arms = {"panel": os.environ["PANEL_K"], "mosaic": os.environ["MOSAIC_K"]}; role = os.environ["ROLE"]

def read_keys(path):
    if not os.path.exists(path):
        return {}
    values = {}
    for line in open(path):
        if "\t" in line:
            key, value = line.rstrip("\n").split("\t", 1)
            values[key] = value
    return values

blocks = open(os.path.join(out, "blocks.tsv"), "w")
blocks.write("role\tdonor\tarm\tcoverage\tblock\tblock_type\tunordered_diploid_edlib\texact\tconditional_GQ\t"
             "genotype_support\tphase_support\n")
runs = open(os.path.join(out, "runs.tsv"), "w")
runs.write("role\tdonor\tarm\tcoverage\tstatus\trun_no_call\tnote\tmarker_termination\tpanel_baseline_complete\t"
           "mosaic_candidates\tmosaic_target_reached\tcandidate_scope\ttruth_in_candidates\ttruth_source\t"
           "blocks_exact\tsv_blocks_exact\tunordered_distance\texcess\tphase_within_candidates\tphase_support\t"
           "phase_audit_status\tmin_block_GQ\tkept_fragments\twall_seconds\tpeak_rss_bytes\n")
for donor in donors:
    for coverage in coverages:
        reads = os.path.join(out, donor, "cov" + coverage)
        subsample = read_keys(os.path.join(reads, "subsample.tsv"))
        resources = read_keys(os.path.join(reads, "lzo.resources.tsv"))
        for arm, k in arms.items():
            prefix = os.path.join(reads, "lzo.K%s" % k)
            status = read_keys(prefix + ".mosaic_status.tsv") or read_keys(os.path.join(reads, "lzo.mosaic_status.tsv"))
            summary = read_keys(prefix + ".eval.summary.tsv")
            complete = status.get("status") == "COMPLETE"
            distances = {}
            if os.path.exists(prefix + ".eval.blocks.tsv"):
                for row in csv.DictReader(open(prefix + ".eval.blocks.tsv"), delimiter="\t"):
                    distances[row["block_id"]] = row["unordered_distance"]
            lowest = ""
            if os.path.exists(prefix + ".block_genotypes.tsv"):
                for row in csv.DictReader(open(prefix + ".block_genotypes.tsv"), delimiter="\t"):
                    distance = distances.get(row["block_id"], "NA")
                    exact = "NA" if distance == "NA" else ("true" if int(distance) == 0 else "false")
                    quality = row["conditional_GQ"]
                    if quality != "inf" and (lowest == "" or float(quality) < float(lowest)):
                        lowest = quality
                    blocks.write("\t".join([role, donor, arm, coverage, row["block_id"], row["kind"], distance,
                                            exact, quality, row["genotype_support"], row["phase_support"]]) + "\n")
            runs.write("\t".join([role, donor, arm, coverage, status.get("status", "NA"),
                                  "false" if complete else "true", status.get("note", "").replace("\t", " ")[:90] or "NA",
                                  status.get("marker_depth_termination", "NA"),
                                  status.get("panel_baseline_complete", "NA"), status.get("mosaic_candidates", "NA"),
                                  status.get("mosaic_target_reached", "NA"), status.get("candidate_scope", "NA"),
                                  summary.get("truth_in_candidates", "NA"), summary.get("truth_candidate_source", "NA"),
                                  summary.get("blocks_exact", "NA"), summary.get("sv_blocks_exact", "NA"),
                                  summary.get("unordered_distance_total", "NA"),
                                  summary.get("excess_over_floor_total", "NA"),
                                  status.get("phase_support_within_candidates", "NA"),
                                  status.get("phase_support", "NA"), status.get("phase_audit_status", "NA"),
                                  lowest or "NA", subsample.get("kept_fragments", "NA"),
                                  resources.get("wall_seconds", "NA"), resources.get("peak_rss_bytes", "NA")]) + "\n")
blocks.close(); runs.close()
COLLECT

# ---- per arm and coverage: no-calls, exact rates, distance, phase and GQ ---------------------------------
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
  { key=$c["arm"]"\t"$c["coverage"]; runs[key]++
    if ($c["run_no_call"]=="true") { nocall[key]++; next }
    called[key]++; exact[key]+=$c["blocks_exact"]; sv[key]+=$c["sv_blocks_exact"]; dist[key]+=$c["unordered_distance"]
    if ($c["unordered_distance"]+0==0) donorexact[key]++
    phase[key"\t"$c["phase_support"]]++ }
  END { printf "arm\tcoverage\truns\trun_no_calls\tdonors_fully_exact\tblocks_exact\tsv_blocks_exact\tdistance_total\n"
        for (k in runs) printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", k, runs[k], nocall[k], donorexact[k], exact[k], sv[k], dist[k] }' \
  "$OUT/runs.tsv" | (read -r header; echo "$header"; sort -t$'\t' -k1,1 -k2,2n) > "$OUT/summary.tsv"
column -t -s $'\t' "$OUT/summary.tsv"
