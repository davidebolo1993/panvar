#!/usr/bin/env bash
# genotype_baseline_cohort.sh - the frozen pre-estimator baseline, six loci, LZO and LOO.
#
#   genotype_baseline_cohort.sh <out-dir> [config.tsv] [loci] [donors]
#
# WHAT THIS IS. The regression baseline the bounded pigeonhole estimator and any later block-level
# caller are measured against. It is CONFIGURATION-DRIVEN: one table of loci, one code path. Six
# copied locus branches would drift, and a drifted branch reports a number for the wrong substrate.
#
# NO TUNING. Every caller parameter is fixed in CALLER_ARGS below. A baseline that was tuned is not
# a baseline. If a parameter must change, that is a NEW baseline with a new commit, not an edit.
#
# TWO TABLES, because this is meant to become a BLOCK-level genotyper:
#   runs.tsv    one row per (locus, donor, arm) -- the whole-locus summary
#   blocks.tsv  one row per block
# A whole-locus median hides a catastrophic KIV-2 block, and one large array hides otherwise exact
# ordinary blocks. Neither table is interpretable without the other.
#
# CONTRACTS, each of which has already been violated at least once on this branch:
#
#  1. THE WALK IS THE SEQUENCE. Truth, candidates and calls are all spelled from the graph walk
#     (--dump-scored-sequences, --spell-pair), never rebuilt from block alleles and never taken from
#     a different pipeline stage. `metric` records exact vs bounded, and the band.
#  2. called_distance >= panel_floor, ALWAYS. A negative excess is not a good result, it is a metric
#     or substrate mismatch, and the run FAILS rather than reporting it.
#  3. THE DECOMPOSITION IS FIXED. LOO excludes candidates; it never recomputes bubbles. Recomputing
#     would change the block indices that blocks.tsv is keyed on, between arms of the same donor.
#  4. PROVENANCE IS RECORDED AND CHECKED: graph, bubble catalogue, truth FASTA, reads, binary. The
#     LPA numbers were once measuring panphorte's KIV-2 folding rather than the genotyper.
#  5. SEED AND READ COUNT ARE RECORDED. This is a DETERMINISTIC regression baseline, one seed. It is
#     not a scientific estimate; a multi-seed study is a separate thing and must be labelled so.
#  6. TRUTH IS AN EQUIVALENCE CLASS. Sequence-identical and reverse-complement-equivalent paths all
#     count as correct; requiring one arbitrary label to win measures path naming, not genotyping.
#  7. FIT IS NORMALISED PER FRAGMENT and split into fragment and exposure/dosage components, and is
#     comparable only within a locus and read regime. Raw log-likelihoods across loci rank by depth.
#  8. REFUSE, DO NOT DEGRADE. Any failed step means NO ROW. A plausible partial row is worse than a
#     missing one because it will be averaged.
#
# NOT COMPUTED, and left blank rather than zero-filled:
#   * per-block edit distance and per-block floor need block-boundary offsets in walk coordinates,
#     which no command exposes today. blocks.tsv reports allele identity and availability instead,
#     which is what determines PASS/OFF_PANEL. A zero here would be a number from an instrument that
#     never ran -- this cohort has produced four of those already.
set -uo pipefail

OUT="${1:?usage: genotype_baseline_cohort.sh <out-dir> [config.tsv] [loci] [donors]}"
CFG="${2:-}"
ONLY_LOCI="${3:-}"
ONLY_DONORS="${4:-}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; REPO="$(cd "$HERE/.." && pwd)"
BIN="${PANVAR_BIN:-$REPO/build_genotype_review/panvar}"
PY="${PYTHON:-python3}"
BAND="${BAND:-4096}"
SEED="${SEED:-7}"
COVERAGE="${COVERAGE:-30}"
mkdir -p "$OUT"

md5of() { if command -v md5 >/dev/null 2>&1; then md5 -q "$1"; else md5sum "$1" | cut -d' ' -f1; fi; }
say() { printf '%s\n' "$*"; }
fails=0; rows=0

# ---- the caller, FROZEN --------------------------------------------------------------------------
CALLER_ARGS=(--haplotype-mode --hamming-emission --max-divergence 0.05
             --fragment-len 350 --fragment-sd 50 --error-rate 0.001 --haploid-depth 0.05
             --max-anchor-occ 32 --placement-topk 8 --top-pairs 20)

# ---- configuration -------------------------------------------------------------------------------
# locus <TAB> graph <TAB> bubble-prefix <TAB> reference-path <TAB> donors (comma-separated)
# The frozen cohort. Written here, in the file, so a baseline is reproducible from the commit alone.
# These six are present with both homologues at every locus in results/real_data.
# THE FOUR DONORS PRESENT WITH BOTH HOMOLOGUES IN ALL SIX PANELS. The six locus graphs were built
# from DIFFERENT sample cohorts -- 232/231/231/231 complete diploid donors at acot/gstm1/lpa/
# ankrd36c against 64 at c4 and 59 at cyp2d6 -- so a donor chosen from the pharmacogene panels is
# usually absent entirely from the others. HG00096/HG00171/HG00268 are in c4 and cyp2d6 only.
# Frozen as a SET, so the cohort is not selected on results. Development/regression donors, not a
# holdout.
FROZEN_DONORS="${FROZEN_DONORS:-HG00733,HG02818,NA19036,NA19240}"
if [ -z "$CFG" ]; then
  CFG="$OUT/config.tsv"
  { for L in acot gstm1 lpa ankrd36c c4 cyp2d6; do
      printf '%s\t%s\t%s\t%s\t%s\n' "$L" \
        "$REPO/results/real_data/$L/bubble/bubble.sorted.gfa" \
        "$REPO/results/real_data/$L/bubble/bubble" "ref" "$FROZEN_DONORS"
    done; } > "$CFG"
fi

[ -x "$BIN" ] || { say "FATAL: binary not executable: $BIN"; exit 1; }
BIN_MD5="$(md5of "$BIN")"
COMMIT="$(cd "$REPO" && git rev-parse --short HEAD 2>/dev/null || echo unknown)"
# Tracked modifications and untracked files are DIFFERENT facts. Counting them together made an
# untracked test fixture look like a source change and reported tracked_dirty 1 on a clean tree.
DIRTY="$(cd "$REPO" && git status --porcelain --untracked-files=no -- src include tests CMakeLists.txt 2>/dev/null | wc -l | tr -d ' ')"
UNTRACKED="$(cd "$REPO" && git ls-files --others --exclude-standard -- src include tests CMakeLists.txt 2>/dev/null | wc -l | tr -d ' ')"
# The binary is COPIED. A rebuild during the run must not change what is being measured -- that has
# invalidated a ctest run and a 16-donor cohort on this branch.
cp "$BIN" "$OUT/panvar.frozen"; BIN="$OUT/panvar.frozen"

{ echo "commit          $COMMIT"
  echo "tracked_dirty   $DIRTY files"
  echo "untracked       $UNTRACKED files (not part of the build)"
  echo "binary_md5      $BIN_MD5"
  echo "band            $BAND"
  echo "seed            $SEED"
  echo "coverage        ${COVERAGE}x"
  echo "caller_args     ${CALLER_ARGS[*]}"
  echo "date            $(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$OUT/PROVENANCE.txt"
cat "$OUT/PROVENANCE.txt"
# A FROZEN run ABORTS on tracked modifications. A baseline that cannot be rebuilt from its recorded
# commit is not a baseline, and a warning is too easy to scroll past. FROZEN=0 for diagnostic runs.
if [ "$DIRTY" != "0" ]; then
  if [ "${FROZEN:-1}" = "1" ]; then
    say "FATAL: $DIRTY tracked source file(s) modified; this baseline would not be reproducible from $COMMIT."
    say "       Commit them, or re-run with FROZEN=0 to take a diagnostic (non-baseline) measurement."
    exit 1
  fi
  say "WARNING: $DIRTY tracked file(s) modified; FROZEN=0, so this is a diagnostic run, NOT a baseline"
fi

# ---- PREFLIGHT ------------------------------------------------------------------------------
# The WHOLE configuration is validated before a single read is simulated or a single row written.
# Discovering at locus five that a donor is absent leaves four loci of results already on disk that
# nobody will delete, and they will be read as "the baseline" with a footnote nobody applies. The
# preflight reports h1_found/h2_found per requested cell and refuses the entire run if any is
# missing, before any work.
PRE="$OUT/preflight.tsv"
printf 'locus	donor	h1_found	h2_found	paths_for_donor	status
' > "$PRE"
pre_fail=0
while IFS=$'	' read -r LOCUS GFA PFX REFNAME DONOR_LIST; do
  [ -z "${LOCUS:-}" ] && continue
  case "$LOCUS" in \#*) continue ;; esac
  if [ -n "$ONLY_LOCI" ] && ! printf '%s' " $ONLY_LOCI " | grep -q " $LOCUS "; then continue; fi
  if [ ! -f "$GFA" ] || [ ! -f "$PFX.bubbles.csv" ]; then
    printf '%s	-	-	-	-	NO_SUBSTRATE
' "$LOCUS" >> "$PRE"; pre_fail=$((pre_fail+1)); continue
  fi
  PD="$OUT/.pre_$LOCUS"; mkdir -p "$PD"
  "$BIN" genotype-frag -i "$GFA" -b "$PFX" -o "$PD/p" --dump-scored-sequences "$PD/p" -q     >/dev/null 2>&1
  if [ ! -s "$PD/p.scored_sequences.tsv" ]; then
    printf '%s	-	-	-	-	NO_PANEL
' "$LOCUS" >> "$PRE"; pre_fail=$((pre_fail+1)); continue
  fi
  DL="${ONLY_DONORS:-$DONOR_LIST}"
  DL="$(printf '%s' "$DL" | tr ',' ' ')"
  if [ -z "$DL" ]; then
    printf '%s	-	-	-	-	NO_MANIFEST
' "$LOCUS" >> "$PRE"; pre_fail=$((pre_fail+1)); continue
  fi
  for D in $DL; do
    H1=$(awk -F'	' -v d="$D" 'NR>1 && index($2,d"#1#")==1{print "yes"; exit}' "$PD/p.scored_sequences.tsv")
    H2=$(awk -F'	' -v d="$D" 'NR>1 && index($2,d"#2#")==1{print "yes"; exit}' "$PD/p.scored_sequences.tsv")
    NP=$(awk -F'	' -v d="$D" 'NR>1 && index($2,d"#")==1{n++} END{print n+0}' "$PD/p.scored_sequences.tsv")
    ST=OK; [ -z "$H1" ] || [ -z "$H2" ] && { ST=MISSING; pre_fail=$((pre_fail+1)); }
    printf '%s	%s	%s	%s	%s	%s
' "$LOCUS" "$D" "${H1:-no}" "${H2:-no}" "$NP" "$ST" >> "$PRE"
  done
  rm -rf "$PD"
done < "$CFG"
column -t "$PRE" 2>/dev/null || cat "$PRE"
if [ "$pre_fail" != 0 ]; then
  say ""
  say "FATAL: preflight failed for $pre_fail cell(s). Nothing was simulated and no rows were written."
  say "       The six locus graphs come from DIFFERENT cohorts -- 232/231/231/231 complete diploid"
  say "       donors at acot/gstm1/lpa/ankrd36c against 64 at c4 and 59 at cyp2d6 -- so a donor"
  say "       present in one panel is often absent from another entirely."
  exit 1
fi
say "preflight: every requested locus/donor cell has a complete diploid"
say ""

RUNS="$OUT/runs.tsv"; BLOCKS="$OUT/blocks.tsv"
printf 'commit\tbinary_md5\tlocus\tdonor\tarm\treads\tseed\tn_read_pairs\tgraph_md5\tbubbles_md5\ttruth_md5\treads_md5\tcatalogue_md5\tshortlist_md5\ttruth1\ttruth2\tpanel_floor\tfloor_pair_survived\tcalled1\tcalled2\tcalled_distance\texcess\tmetric\tequiv_set_size\tfit_per_frag\tfit_frag\tfit_exposure\tn_blocks\tn_projectable\tn_determined\tn_exact\tn_wrong\truntime_s\tpeak_rss_bytes\tblock_sum_distance\tblock_sum_residual\tirreducible_bp\trecoverable_bp\tstatus\n' > "$RUNS"
printf 'commit\tlocus\tdonor\tarm\tblock\tkind\tbubble_id\tn_alleles\ttruth_a1\ttruth_a2\tcalled_a1\tcalled_a2\tprojectable\tdetermined\tblock_equivalence_size\tstatus\tblock_called_distance\tblock_floor\tblock_excess\ttruth_bp\tcalled_bp\tmetric\n' > "$BLOCKS"

# ---- peak RSS, units made explicit ---------------------------------------------------------------
# Darwin /usr/bin/time -l reports BYTES; GNU time -v reports KILOBYTES. Recording the number without
# the unit makes a 1000x difference invisible in a table that will be compared across machines.
peak_rss_bytes() {   # peak_rss_bytes <timefile>
  local f="$1" v
  v=$(awk '/maximum resident set size/ {print $1; exit}' "$f" 2>/dev/null)
  if [ -n "$v" ]; then echo "$v"; return; fi
  v=$(awk -F': *' '/Maximum resident set size/ {print $2; exit}' "$f" 2>/dev/null)
  if [ -n "$v" ]; then echo $(( v * 1024 )); return; fi
  echo ""
}
# run_timed <timefile> <cmd...> -> propagates the command's exit code
run_timed() {
  local tf="$1"; shift
  /usr/bin/time -l "$@" 2>"$tf" >/dev/null
  return $?
}

# ---- helpers -------------------------------------------------------------------------------------
# Write one named record out of a scored_sequences.fa, uppercased and in REFERENCE ORIENTATION.
#
# scored_sequences.fa carries WALK bytes, and a path whose frame is `rc` runs antiparallel to the
# reference/block frame -- its walk is the reverse complement of every forward path's. Comparing
# walk bytes across frames therefore measures STRAND, not genotype: c4/HG00171 leave-ZERO-out, with
# its own truth in the panel and all 11 blocks matching, reported a distance of 240094. Every
# distance here (truth, candidate, called) is oriented by the frame column the binary itself
# reports, which is the same convention the block projection emits.
extract_fa() {   # extract_fa <fa> <tsv> <name> <out>
  "$PY" - "$1" "$2" "$3" "$4" <<'PY'
import sys
fa,tsv,want,out=sys.argv[1:5]
frame={}
h=None
for l in open(tsv):
    f=l.rstrip('\n').split('\t')
    if h is None: h=f; continue
    d=dict(zip(h,f)); frame[d['name']]=d.get('frame','fwd')
n=None;buf=[];hit=None
def flush():
    global hit
    if n is not None and n==want: hit="".join(buf)
for l in open(fa):
    if l[0]=='>':
        flush()
        if hit is not None: break
        n=l[1:].split()[0]; buf=[]
    else: buf.append(l.strip())
flush()
if hit is None: sys.exit(3)
seq=hit.upper()
if frame.get(want,'fwd')=='rc':
    seq=seq[::-1].translate(str.maketrans('ACGTN','TGCAN'))
open(out,'w').write(">s\n"+seq+"\n")
PY
}
# Exact edit distance between two one-record FASTAs, via the SAME command every other distance uses.
# Escalates the band until the distance is EXACT, and reports the band it needed. A distance that
# fell outside the band is not a large distance, it is no distance -- and silently treating the cap
# as the answer would make a badly wrong call look like a near miss.
#
# The band actually used is written to a FILE, not a variable. dist() runs inside a command
# substitution -- a subshell -- so an assignment there never reaches the caller: the metric column
# reported band=4096 for a distance of 93723, which is arithmetically impossible and went unnoticed
# because the DISTANCES were right. Only the certification was wrong.
BANDFILE=""
dist() {   # dist <a.fa> <b.fa> -> integer, or empty
  local v b
  for b in $BAND 16384 65536 262144; do
    v=$("$BIN" genotype-frag --exact-distance "$1" "$2" --distance-band "$b" 2>/dev/null | tr -d ' ')
    case "$v" in ''|*[!0-9]*) continue ;; esac
    if [ -n "$BANDFILE" ]; then
      prev=$(cat "$BANDFILE" 2>/dev/null || echo 0)
      [ "$b" -gt "${prev:-0}" ] && printf '%s' "$b" > "$BANDFILE"
    fi
    echo "$v"; return
  done
  echo ""
}

while IFS=$'\t' read -r LOCUS GFA PFX REFNAME DONOR_LIST; do
  [ -z "${LOCUS:-}" ] && continue
  case "$LOCUS" in \#*) continue ;; esac
  if [ -n "$ONLY_LOCI" ] && ! printf '%s' " $ONLY_LOCI " | grep -q " $LOCUS "; then continue; fi
  if [ ! -f "$GFA" ] || [ ! -f "$PFX.bubbles.csv" ]; then
    say "REFUSE $LOCUS: missing substrate ($GFA / $PFX.bubbles.csv)"; fails=$((fails+1)); continue
  fi
  GFA_MD5="$(md5of "$GFA")"; BUB_MD5="$(md5of "$PFX.bubbles.csv")"
  LD="$OUT/$LOCUS"; rm -rf "$LD"; mkdir -p "$LD"

  # Panel inventory, from the graph walks. This is also the candidate catalogue for LZO.
  "$BIN" genotype-frag -i "$GFA" -b "$PFX" -o "$LD/all" --dump-scored-sequences "$LD/all" -q \
    >/dev/null 2>&1
  if [ ! -s "$LD/all.scored_sequences.fa" ]; then
    say "REFUSE $LOCUS: could not spell the panel from the graph"; fails=$((fails+1)); continue
  fi
  # DONORS ARE FROZEN IN THE CONFIGURATION, field 5. `sort | head -N` is not a cohort: it silently
  # changes when a path is added to or removed from the panel, and a baseline whose membership moves
  # is not a baseline. A locus with no manifest is REFUSED rather than filled in by discovery.
  DONORS="${DONOR_LIST:-}"
  [ -n "$ONLY_DONORS" ] && DONORS="$ONLY_DONORS"
  if [ -z "$DONORS" ]; then
    say "REFUSE $LOCUS: no frozen donor manifest in the configuration (field 5)"
    say "       candidates with two homologues: $(awk -F'\t' 'NR>1{split($2,a,"#"); if(a[2]=="1"||a[2]=="2") c[a[1]]=c[a[1]] a[2]} END{for(s in c) if(length(c[s])==2) printf "%s ", s}' "$LD/all.scored_sequences.tsv" | tr ' ' '\n' | sort | tr '\n' ' ')"
    fails=$((fails+1)); continue
  fi
  DONORS="$(printf '%s' "$DONORS" | tr ',' ' ')"

  for DONOR in $DONORS; do
    T1=$(awk -F'\t' -v d="$DONOR" 'NR>1 && index($2,d"#1#")==1 {print $2; exit}' "$LD/all.scored_sequences.tsv")
    T2=$(awk -F'\t' -v d="$DONOR" 'NR>1 && index($2,d"#2#")==1 {print $2; exit}' "$LD/all.scored_sequences.tsv")
    if [ -z "$T1" ] || [ -z "$T2" ]; then
      # "does not have two homologues" was misleading: most often the donor is absent from this
      # panel ENTIRELY -- zero paths, not one -- because the six graphs come from different cohorts.
      NFOUND=$(awk -F'\t' -v d="$DONOR" 'NR>1 && index($2,d"#")==1{n++} END{print n+0}' "$LD/all.scored_sequences.tsv")
      say "REFUSE $LOCUS/$DONOR: $NFOUND path(s) for this donor (h1=$([ -n "$T1" ] && echo yes || echo no) h2=$([ -n "$T2" ] && echo yes || echo no)); the panel does not carry a complete diploid"
      fails=$((fails+1)); continue
    fi
    DD="$LD/$DONOR"; mkdir -p "$DD"
    extract_fa "$LD/all.scored_sequences.fa" "$LD/all.scored_sequences.tsv" "$T1" "$DD/t1.fa" || { say "REFUSE $LOCUS/$DONOR: cannot spell truth 1"; fails=$((fails+1)); continue; }
    extract_fa "$LD/all.scored_sequences.fa" "$LD/all.scored_sequences.tsv" "$T2" "$DD/t2.fa" || { say "REFUSE $LOCUS/$DONOR: cannot spell truth 2"; fails=$((fails+1)); continue; }
    cat "$DD/t1.fa" "$DD/t2.fa" > "$DD/truth.fa"; TRUTH_MD5="$(md5of "$DD/truth.fa")"

    # TRUTH EQUIVALENCE CLASS: every panel path whose walk is sequence-identical to a truth
    # haplotype, or its reverse complement. Requiring one arbitrary label to win measures naming.
    "$PY" - "$LD/all.scored_sequences.fa" "$DD/t1.fa" "$DD/t2.fa" > "$DD/truth_class.txt" <<'PY'
import sys
def one(p):
    return "".join(l.strip() for l in open(p) if l[0]!='>').upper()
def rc(s): return s[::-1].translate(str.maketrans('ACGT','TGCA'))
fa,p1,p2=sys.argv[1:4]
t=[one(p1),one(p2)]; t+= [rc(x) for x in t]
n=None;buf=[]
def flush():
    if n is not None and "".join(buf).upper() in t: print(n)
for l in open(fa):
    if l[0]=='>': flush(); n=l[1:].split()[0]; buf=[]
    else: buf.append(l.strip())
flush()
PY

    # READS: simulated, one fixed seed, count recorded. A deterministic regression baseline.
    # 600, not 300. truth.fa is DIPLOID, so lambda = N/L_diploid. With N = C*L/300 that is
    # lambda = 0.10 while the caller is told --haploid-depth 0.05: the reads carried twice the depth
    # the model assumed, and every dosage/exposure term was fitted against the wrong rate.
    # C*L/600 gives lambda = 0.05. Verified arithmetically, not assumed.
    NPAIRS=$(( COVERAGE * $(awk '/^>/{next}{n+=length($0)}END{print n}' "$DD/truth.fa") / 600 ))
    if ! wgsim -N "$NPAIRS" -1 150 -2 150 -d 350 -s 50 -e 0.001 -r 0 -R 0 -X 0 -S "$SEED" \
         "$DD/truth.fa" "$DD/r1.fq" "$DD/r2.fq" >/dev/null 2>&1; then
      say "REFUSE $LOCUS/$DONOR: wgsim failed"; fails=$((fails+1)); continue
    fi
    # BOTH mates. Hashing R1 alone leaves R2 unprovenanced, and R2 is half the evidence.
    READS_MD5="$(cat "$DD/r1.fq" "$DD/r2.fq" | { if command -v md5 >/dev/null 2>&1; then md5 -q; else md5sum | cut -d' ' -f1; fi; })"

    for ARM in LZO LOO; do
      AD="$DD/$ARM"; rm -rf "$AD"; mkdir -p "$AD"
      EXCL=(); [ "$ARM" = LOO ] && EXCL=(--exclude-haplotypes "$T1,$T2")

      # The candidate catalogue AFTER exclusion. The floor must be computed against exactly this
      # set: including the held-out truth makes every floor 0, which is uniform enough to look right.
      "$BIN" genotype-frag -i "$GFA" -b "$PFX" -o "$AD/cat" --dump-scored-sequences "$AD/cat" \
        ${EXCL[@]+"${EXCL[@]}"} -q >/dev/null 2>&1
      if [ ! -s "$AD/cat.scored_sequences.fa" ]; then
        say "REFUSE $LOCUS/$DONOR/$ARM: no candidate catalogue"; fails=$((fails+1)); continue
      fi
      CAT_MD5="$(md5of "$AD/cat.scored_sequences.fa")"

      # PANEL FLOOR: the best any remaining candidate PAIR can do, exact, same metric as everything
      # else. Filtered on the dump's group column -- the held_out group IS the truth.
      FLOORLINE=$("$PY" - "$AD/cat.scored_sequences.fa" "$AD" "$BIN" "$BAND" "$DD/t1.fa" "$DD/t2.fa" \
                       "$AD/cat.scored_sequences.tsv" <<'PY'
import sys, subprocess
fa,ad,BIN,band,t1,t2,tsv=sys.argv[1:8]
frame={};h=None
for l in open(tsv):
    f=l.rstrip('\n').split('\t')
    if h is None: h=f; continue
    d=dict(zip(h,f)); frame[d['name']]=d.get('frame','fwd')
def orient(nm,s):
    s=s.upper()
    return s[::-1].translate(str.maketrans('ACGTN','TGCAN')) if frame.get(nm,'fwd')=='rc' else s
names=[];seqs=[];n=None;grp=None;buf=[]
def flush():
    if n is not None and grp!="held_out": names.append(n); seqs.append("".join(buf))
for l in open(fa):
    if l[0]=='>':
        flush(); p=l[1:].split(); n=p[0]; grp=p[1] if len(p)>1 else "panel"; buf=[]
    else: buf.append(l.strip())
flush()
if not seqs: print("NOFLOOR"); sys.exit()
best=[10**9,10**9]; who=[None,None]
for nm,s in zip(names,seqs):
    open(f"{ad}/c.fa","w").write(">c\n"+orient(nm,s)+"\n")
    for h,t in ((0,t1),(1,t2)):
        # ESCALATE, like the called distance does. A floor that refuses at the default band aborts
        # the whole donor -- ankrd36c/NA19240/LOO was lost that way -- and a candidate outside the
        # band is not a distant candidate, it is an unmeasured one.
        v=None
        for bd in (band,'16384','65536','262144'):
            r=subprocess.run([BIN,"genotype-frag","--exact-distance",f"{ad}/c.fa",t,
                              "--distance-band",bd],capture_output=True,text=True)
            o=r.stdout.strip()
            if o.isdigit(): v=int(o); break
        if v is not None and v<best[h]: best[h]=v; who[h]=nm
if 10**9 in best: print("NOFLOOR")
else: print("%d\t%s\t%s" % (best[0]+best[1], who[0], who[1]))
PY
)
      case "$FLOORLINE" in
        NOFLOOR|"") say "REFUSE $LOCUS/$DONOR/$ARM: floor not computable within band $BAND"
                    fails=$((fails+1)); continue ;;
      esac
      read -r PANEL_FLOOR FLOOR_H1 FLOOR_H2 <<<"$FLOORLINE"

      # THE CALLER, timed, into a directory guaranteed empty. A stale output from a previous run
      # read back as this run's result is a silent way to report the wrong number.
      rm -f "$AD/call".* 
      if ! run_timed "$AD/time.txt" \
           "$BIN" genotype-frag -i "$GFA" -b "$PFX" -o "$AD/call" \
             -R "$DD/r1.fq" -R "$DD/r2.fq" ${EXCL[@]+"${EXCL[@]}"} \
             "${CALLER_ARGS[@]}" -t "${THREADS:-4}" -q; then
        say "REFUSE $LOCUS/$DONOR/$ARM: caller exited nonzero"; fails=$((fails+1)); continue
      fi
      if [ ! -s "$AD/call.hap_pairs.tsv" ]; then
        say "REFUSE $LOCUS/$DONOR/$ARM: caller produced no hap_pairs.tsv"; fails=$((fails+1)); continue
      fi
      # Darwin `time -l` prints "        0.34 real         0.30 user" -- `real` is field 2 and the
      # number is field 1. Matching /^ *real/ silently yields an empty column.
      RUNTIME=$(awk '$2=="real" {print $1; exit}' "$AD/time.txt" 2>/dev/null)
      [ -z "$RUNTIME" ] && RUNTIME=$(awk -F': *' '/Elapsed \(wall clock\)/ {print $2; exit}' "$AD/time.txt")
      PEAK=$(peak_rss_bytes "$AD/time.txt")

      # SPELL BY NAME, from the walks. Allele indices are valid for one catalogue only and spelling
      # across a mismatch is silent.
      "$BIN" genotype-frag -i "$GFA" -b "$PFX" -o "$AD/sp" --spell-pair "$AD/call.hap_pairs.tsv" \
        -q >/dev/null 2>&1
      if [ ! -s "$AD/sp.called.fa" ]; then
        say "REFUSE $LOCUS/$DONOR/$ARM: could not spell the called pair by name"; fails=$((fails+1)); continue
      fi
      C1=$(awk -F'\t' 'NR==2{print $2}' "$AD/call.hap_pairs.tsv")
      C2=$(awk -F'\t' 'NR==2{print $3}' "$AD/call.hap_pairs.tsv")
      # --spell-pair ALREADY emits reference orientation -- verified: for the antiparallel
      # c4/HG00171 pair, its record equals the frame-corrected truth byte for byte. Only
      # scored_sequences.fa carries raw walk bytes and needs orienting. Applying the correction here
      # too double-flips the called pair, which leaves the distance wrong by exactly as much as
      # doing nothing did, and looks like the fix having no effect.
      "$PY" - "$AD/sp.called.fa" "$AD/c1.fa" "$AD/c2.fa" <<'PY'
import sys
fa,o1,o2=sys.argv[1:4]
recs=[];n=None;buf=[]
for l in open(fa):
    if l[0]=='>':
        if n is not None: recs.append("".join(buf))
        n=1;buf=[]
    else: buf.append(l.strip())
if n is not None: recs.append("".join(buf))
if len(recs)<2: sys.exit(3)
open(o1,'w').write(">c1\n"+recs[0].upper()+"\n")
open(o2,'w').write(">c2\n"+recs[1].upper()+"\n")
PY
      [ -s "$AD/c1.fa" ] && [ -s "$AD/c2.fa" ] || { say "REFUSE $LOCUS/$DONOR/$ARM: called pair not spelled as two records"; fails=$((fails+1)); continue; }
      # best of the two orientation assignments -- the pair is unordered
      BANDFILE="$AD/used_band"; printf '0' > "$BANDFILE"
      D11=$(dist "$AD/c1.fa" "$DD/t1.fa"); D22=$(dist "$AD/c2.fa" "$DD/t2.fa")
      D12=$(dist "$AD/c1.fa" "$DD/t2.fa"); D21=$(dist "$AD/c2.fa" "$DD/t1.fa")
      CALLED_D=$("$PY" -c "
import sys
v=[]
for a,b in (('$D11','$D22'),('$D12','$D21')):
    if a.isdigit() and b.isdigit(): v.append(int(a)+int(b))
print(min(v) if v else '')" 2>/dev/null)
      if [ -z "$CALLED_D" ]; then
        say "REFUSE $LOCUS/$DONOR/$ARM: called distance not computable within band $BAND"
        fails=$((fails+1)); continue
      fi
      # CONTRACT: NAME IDENTITY IMPLIES SEQUENCE IDENTITY. If the caller returned the truth pair by
      # name, the distance to truth must be exactly 0, whatever the frames involved. This is a check
      # on the METRIC, not on the caller, and it is what finally caught the orientation defect:
      # c4/HG00171 leave-zero-out returned its own truth pair and scored it 240094, because
      # scored_sequences.fa carries walk bytes and both truth paths are antiparallel while
      # --spell-pair emits reference orientation. The floor was 0 either way, so no floor-based
      # invariant could see it, and only a donor with an rc-frame haplotype exposes it at all.
      if "$PY" -c "
import sys
t={'$T1','$T2'}; c={'$C1','$C2'}
sys.exit(0 if t==c else 1)" 2>/dev/null && [ "$CALLED_D" != 0 ]; then
        say "REFUSE $LOCUS/$DONOR/$ARM: the called pair IS the truth pair by name, but distance is $CALLED_D -- the metric is comparing different orientations, not different genotypes"
        fails=$((fails+1)); continue
      fi

      # CONTRACT 2. A negative excess is a metric or substrate mismatch, never a good result.
      if [ "$CALLED_D" -lt "$PANEL_FLOOR" ]; then
        say "REFUSE $LOCUS/$DONOR/$ARM: called_distance $CALLED_D < panel_floor $PANEL_FLOOR -- metric or substrate mismatch"
        fails=$((fails+1)); continue
      fi
      EXCESS=$(( CALLED_D - PANEL_FLOOR ))
      SHORTLIST_MD5=$(awk -F'\t' 'NR>1{print $1}' "$AD/call.hap_scores.tsv" 2>/dev/null | sort | \
                      { if command -v md5 >/dev/null 2>&1; then md5 -q; else md5sum | cut -d' ' -f1; fi; })
      # floor pair survives into the SHORTLIST when both floor haplotypes are scored there
      # EXACT field match. grep -F on a whole path name still matches a longer name that contains
      # it, and these names share long prefixes, so a substring hit would report a floor pair as
      # surviving when a different haplotype was in the shortlist.
      FSURV=yes
      for FH in "$FLOOR_H1" "$FLOOR_H2"; do
        awk -F'\t' -v n="$FH" 'NR>1 && $1==n {found=1} END{exit !found}' \
          "$AD/call.hap_scores.tsv" 2>/dev/null || FSURV=no
      done
      EQS=$(awk -F'\t' 'NR==2{print $2}' "$AD/call.equivalence.tsv" 2>/dev/null)

      # NORMALISED FIT. Definition, since the denominator is a real choice and not a detail:
      #
      #     fit = mean_f [ log P(f | called pair) - log P_bg(f) ]
      #
      # PURE background as the denominator, not the mixture-weighted eta*P_bg. The two differ by a
      # constant -log(eta) per fragment, which is arbitrary and would ride on every locus's number
      # while looking like signal. A mixture-weighted variant may be kept, but must be NAMED as such
      # rather than reported under the same column.
      #
      # NOT COMPUTED YET, but no longer blocked: --dump-fragment-mass works on the production path
      # (its own header says "production path (no --joint-depth)"), so the plan doc's claim that it
      # requires --joint-depth is stale. What remains is to wire the dump into this harness and
      # calibrate the distribution separately per locus and depth against the LZO controls -- an
      # absolute fit compared across loci or depths without that calibration is not comparable.
      # OLD NOTE, retained because the reasoning still applies to the previous attempt:
      # The previous version summed the two haplotypes' solo_ll and called the remainder
      # "exposure". That is invalid: the diploid per-fragment term is
      # log[(1-eta)*lambda*(M_a+M_b) + eta*P_bg], which is not the sum of two haploid terms -- the
      # background is counted twice and the fragment mass is combined before the log, not after.
      # A valid split needs the scorer's own terms (fragment_sum, coverage_a, coverage_b, dosage,
      # total_score) normalised by the LOADED fragment count, and --dump-fragment-mass reaches them
      # only under --joint-depth, not on the production path. That is the known instrumentation gap.
      # Emitting NA rather than an invalid number: this cohort has produced five figures from
      # instruments that never ran, and a plausible-looking fit would be the sixth.
      FIT=NA; FITF=NA; FITE=NA

      # ---- BLOCKS, from the binary's authoritative projection ------------------------------------
      # $AD/cat.path_blocks.tsv comes from the SAME --dump-scored-sequences run that produced the
      # candidate catalogue, with this arm's exclusion already applied. So truth arrives as group
      # held_out, the called pair as group panel, and both are projected against the reduced calling
      # catalogue. No second invocation, and no second projector.
      #
      # The previous implementation used a Python projector that assumed source->sink order in walk
      # coordinates and returned ABSENT for every block of every ANTIPARALLEL path -- 60 of 131 at
      # c4, 59 of 127 at cyp2d6. HG00096 is forward at both loci, which is the only reason the first
      # validated donor looked right. Comparison is on block_md5, which is emitted in
      # reference/block orientation and so is directly comparable across a forward and a reverse
      # path; canonical_md5 exists for identity across whole-path orientation and is not used here.
      NB=""; NPROJ=""; NDET=""; NEX=""; NWR=""
      "$PY" - "$AD/cat.path_blocks.tsv" "$AD/cat.path_blocks.fa" "$AD/call.hap_blocks.tsv" \
             "$COMMIT" "$LOCUS" "$DONOR" "$ARM" "${EQS:-1}" "$T1" "$T2" "$C1" "$C2" \
             "$BIN" "$AD" "$BAND" \
             >> "$BLOCKS" 2>"$AD/bstat.txt" <<'PY'
import sys, hashlib, subprocess, os
tsv,fasta,armcall,commit,locus,donor,arm,eqs,t1,t2,c1,c2,BIN,AD,band=sys.argv[1:16]
want=[t1,t2,c1,c2]
def die(msg):
    sys.stderr.write("0 0 0 0 0\n"); raise SystemExit("block projection: "+msg)
rows={};hdr=None
for l in open(tsv):
    f=l.rstrip('\n').split('\t')
    if hdr is None: hdr=f; continue
    d=dict(zip(hdr,f))
    # EXACT path name. A substring or prefix match would pick up a different haplotype of the same
    # sample, and these names share long prefixes.
    if d['path'] in want: rows.setdefault(d['path'],[]).append(d)
per={}
for nm in want:
    rs=rows.get(nm)
    if not rs: die("no rows for "+nm)
    if any(r['projection_status']=='unprojectable' for r in rs): die("unprojectable: "+nm)
    m={}
    for r in rs:
        if r['projection_status']=='unmapped': continue
        b=int(r['block'])
        if b in m: die("duplicate row for %s block %d" % (nm,b))
        m[b]=r
    if not m: die("no mapped blocks for "+nm)
    per[nm]=m
sets=[set(per[nm]) for nm in want]
if any(s!=sets[0] for s in sets):
    die("block sets differ between requested paths")
# FASTA must agree with the TSV: length and md5, per slice. A dump whose two halves disagree is
# not a dump, and nothing downstream would notice.
seqs={};h=None;buf=[]
for l in open(fasta):
    if l[0]=='>':
        if h: seqs[h]=''.join(buf)
        p=l[1:].rstrip('\n').split()
        h=(p[0], int([x for x in p if x.startswith('block=')][0][6:])); buf=[]
    else: buf.append(l.strip())
if h: seqs[h]=''.join(buf)
for nm in want:
    for b,r in per[nm].items():
        s=seqs.get((nm,b))
        if s is None: die("no FASTA record for %s block %d" % (nm,b))
        if len(s)!=int(r['block_bp']): die("FASTA length disagrees with TSV at %s block %d" % (nm,b))
        if hashlib.md5(s.encode()).hexdigest()!=r['block_md5']:
            die("FASTA md5 disagrees with TSV at %s block %d" % (nm,b))
A={};ah=None
for l in open(armcall):
    if l.startswith('#'): continue
    f=l.rstrip('\n').split('\t')
    if ah is None: ah=f; continue
    d=dict(zip(ah,f)); A[int(d['block'])]=d
# PER-BLOCK EXACT DISTANCE AND FLOOR, under the SAME metric and the same unordered diploid rule as
# the whole-locus number: min[d(t1,c1)+d(t2,c2), d(t1,c2)+d(t2,c1)].
#
# The floor is over PANEL ALLELES at the block, chosen independently per homologue. That is exact
# rather than a lower bound, because at a single block any two panel alleles are jointly realisable:
# they are carried by some pair of candidates, and nothing at this level couples them. Linkage
# between blocks is precisely what the whole-locus floor measures and this one does not -- which is
# why the two are reported separately and their difference is not an error.
alleles={}   # block -> {md5: seq} over every PANEL path, from the same dump
for l in open(tsv):
    pass
seqs_all={}
h=None;buf=[];key=None
for l in open(fasta):
    if l[0]=='>':
        if key: seqs_all[key]=''.join(buf)
        p=l[1:].rstrip('\n').split()
        key=(p[0], int([x for x in p if x.startswith('block=')][0][6:])); buf=[]
    else: buf.append(l.strip())
if key: seqs_all[key]=''.join(buf)
grp={}
hh=None
for l in open(tsv):
    f=l.rstrip('\n').split('\t')
    if hh is None: hh=f; continue
    d=dict(zip(hh,f))
    if d['projection_status']!='complete' and d['projection_status']!='partial': continue
    if d['group']!='panel': continue
    try: b=int(d['block'])
    except ValueError: continue
    sq=seqs_all.get((d['path'],b))
    if sq is not None: alleles.setdefault(b,{})[d['block_md5']]=sq
_dcache={}
def dist(a,b):
    if a==b: return 0
    k=(hashlib.md5(a.encode()).hexdigest(),hashlib.md5(b.encode()).hexdigest())
    if k in _dcache: return _dcache[k]
    fa=os.path.join(AD,'_da.fa'); fb=os.path.join(AD,'_db.fa')
    open(fa,'w').write(">a\n"+a+"\n"); open(fb,'w').write(">b\n"+b+"\n")
    v=None
    for bd in (band,'16384','65536','262144'):
        r=subprocess.run([BIN,"genotype-frag","--exact-distance",fa,fb,"--distance-band",bd],
                         capture_output=True,text=True)
        t=r.stdout.strip()
        if t.isdigit(): v=int(t); break
    _dcache[k]=v
    return v
def pair_dist(x1,x2,y1,y2):
    a=[dist(x1,y1),dist(x2,y2)]; b=[dist(x1,y2),dist(x2,y1)]
    v=[]
    if None not in a: v.append(a[0]+a[1])
    if None not in b: v.append(b[0]+b[1])
    return min(v) if v else None
def na(x): return "NA" if x is None else str(x)
nb=nproj=ndet=nex=nwr=0
sum_bd=0; sum_ok=True; sum_bf=0; sum_bx=0
for b in sorted(sets[0]):
    a=A.get(b,{})
    r1,r2,q1,q2=(per[t1][b],per[t2][b],per[c1][b],per[c2][b])
    # A RETAINED path must have index and representability agreeing. A HELD-OUT path may legitimately
    # be unrepresentable -- that is the interesting state, not an error.
    for r in (q1,q2):
        if (r['catalogue_allele']=='NA') != (r['catalogue_representable']=='0'):
            die("called path block %d: index and representability disagree" % b)
    nb+=1
    pj=a.get('projectable',''); dt=a.get('determined','')
    if pj=='yes': nproj+=1
    if dt=='1': ndet+=1
    tset={r1['block_md5'],r2['block_md5']}; cset={q1['block_md5'],q2['block_md5']}
    st='PASS' if tset==cset else 'WRONG'
    if st=='PASS': nex+=1
    else: nwr+=1
    ts1=seqs_all.get((t1,b),''); ts2=seqs_all.get((t2,b),'')
    cs1=seqs_all.get((c1,b),''); cs2=seqs_all.get((c2,b),'')
    bd=pair_dist(ts1,ts2,cs1,cs2)
    # Floor: best panel allele for each truth homologue, independently.
    cand=list(alleles.get(b,{}).values())
    if cand:
        d1=[dist(ts1,x) for x in cand]; d2=[dist(ts2,x) for x in cand]
        d1=[x for x in d1 if x is not None]; d2=[x for x in d2 if x is not None]
        bf=(min(d1)+min(d2)) if (d1 and d2) else None
    else:
        bf=None
    bx=None if (bd is None or bf is None) else bd-bf
    if bd is None: sum_ok=False
    else: sum_bd+=bd
    if bf is not None: sum_bf+=bf
    if bx is not None: sum_bx+=bx
    if bx is not None and bx<0:
        sys.stderr.write("0 0 0 0 0\n")
        raise SystemExit("block %d: called distance %d below floor %d" % (b,bd,bf))
    eqb=a.get('block_equivalence_size','NA')
    sys.stdout.write('\t'.join([commit,locus,donor,arm,str(b),
        a.get('kind',r1.get('kind','')),r1.get('bubble_id',''),a.get('n_alleles',''),
        r1['block_md5'][:8],r2['block_md5'][:8],q1['block_md5'][:8],q2['block_md5'][:8],
        pj,dt,eqb,st,na(bd),na(bf),na(bx),str(len(ts1))+"/"+str(len(ts2)),
        str(len(cs1))+"/"+str(len(cs2)),"exact"])+'\n')
# THE DECOMPOSITION, PER EDIT AND NOT PER BLOCK. Classifying a whole block by whether its floor is
# zero assigns the ENTIRE called distance to "unrepresentable" whenever the floor is even 1 edit, so
# a block whose nearest panel allele is 10 edits away and whose call is 10000 away contributes 10000
# irreducible instead of 10. That inverted the headline of the first sweep: lpa looked like a 94835
# edit coverage dead end when it is ~1075 irreducible and ~94k of recoverable headroom.
#   irreducible = sum(block_floor)      recoverable = sum(block_excess)
sys.stderr.write("%d %d %d %d %d %s %s %s\n" % (nb,nproj,ndet,nex,nwr,
                 str(sum_bd) if sum_ok else "NA", str(sum_bf), str(sum_bx)))
PY
      if [ ! -s "$AD/bstat.txt" ]; then
        say "REFUSE $LOCUS/$DONOR/$ARM: block projection produced no status"; fails=$((fails+1)); continue
      fi
      read -r NB NPROJ NDET NEX NWR BSUM BIRR BREC <<<"$(head -1 "$AD/bstat.txt")"
      # An aborted projection writes 0 0 0 0 0 and a reason. It must REFUSE the donor, not leave the
      # block columns blank and carry on -- a run row with no blocks would still be averaged.
      if [ "${NB:-0}" = 0 ]; then
        say "REFUSE $LOCUS/$DONOR/$ARM: $(tail -1 "$AD/bstat.txt")"; fails=$((fails+1)); continue
      fi

      # BLOCK-SUM RESIDUAL, reported and NOT enforced. The sum of per-block global alignments need
      # not equal the whole-locus global alignment: a locus alignment can place gaps across a block
      # boundary that no per-block alignment is allowed to. A nonzero residual is information about
      # how much of the distance is boundary-crossing, not a discrepancy to fix. Requiring equality
      # would be wrong and would fail on every indel spanning a bubble edge.
      if [ "${BSUM:-NA}" != "NA" ]; then BRESID=$(( CALLED_D - BSUM )); else BRESID=NA; fi
      printf '%s\t%s\t%s\t%s\t%s\tsim\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\texact(band=%s)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tOK\n' \
        "$COMMIT" "$BIN_MD5" "$LOCUS" "$DONOR" "$ARM" "$SEED" "$NPAIRS" \
        "$GFA_MD5" "$BUB_MD5" "$TRUTH_MD5" "$READS_MD5" "$CAT_MD5" "${SHORTLIST_MD5:-}" \
        "$T1" "$T2" "$PANEL_FLOOR" "$FSURV" "$C1" "$C2" "$CALLED_D" "$EXCESS" "$(cat "$AD/used_band" 2>/dev/null || echo "$BAND")" \
        "${EQS:-}" "$FIT" "$FITF" "$FITE" \
        "${NB:-}" "${NPROJ:-}" "${NDET:-}" "${NEX:-}" "${NWR:-}" \
        "${RUNTIME:-}" "${PEAK:-}" "${BSUM:-NA}" "${BRESID:-NA}" "${BIRR:-NA}" "${BREC:-NA}" >> "$RUNS"
      rows=$((rows+1))
      say "  $LOCUS/$DONOR/$ARM floor=$PANEL_FLOOR called=$CALLED_D excess=$EXCESS blocks=${NEX:-?}/${NB:-?} exact"
    done
  done
done < "$CFG"

say ""
say "rows written: $rows    refusals: $fails"
say "  $RUNS"
say "  $BLOCKS"
[ "$rows" -gt 0 ] || { say "NO ROWS -- the baseline is empty, which is a failure, not an empty result"; exit 1; }
exit $(( fails > 0 ? 1 : 0 ))
