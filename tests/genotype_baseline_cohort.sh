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
FROZEN_DONORS="${FROZEN_DONORS:-HG00096,HG00171,HG00268}"
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
DIRTY="$(cd "$REPO" && git status --porcelain -- src include tests CMakeLists.txt 2>/dev/null | wc -l | tr -d ' ')"
# The binary is COPIED. A rebuild during the run must not change what is being measured -- that has
# invalidated a ctest run and a 16-donor cohort on this branch.
cp "$BIN" "$OUT/panvar.frozen"; BIN="$OUT/panvar.frozen"

{ echo "commit          $COMMIT"
  echo "tracked_dirty   $DIRTY files"
  echo "binary_md5      $BIN_MD5"
  echo "band            $BAND"
  echo "seed            $SEED"
  echo "coverage        ${COVERAGE}x"
  echo "caller_args     ${CALLER_ARGS[*]}"
  echo "date            $(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$OUT/PROVENANCE.txt"
cat "$OUT/PROVENANCE.txt"
[ "$DIRTY" != "0" ] && say "WARNING: $DIRTY tracked source file(s) modified; this baseline is not reproducible from $COMMIT alone"

RUNS="$OUT/runs.tsv"; BLOCKS="$OUT/blocks.tsv"
printf 'commit\tbinary_md5\tlocus\tdonor\tarm\treads\tseed\tn_read_pairs\tgraph_md5\tbubbles_md5\ttruth_md5\treads_md5\tcatalogue_md5\tshortlist_md5\ttruth1\ttruth2\tpanel_floor\tfloor_pair_survived\tcalled1\tcalled2\tcalled_distance\texcess\tmetric\tequiv_set_size\tfit_per_frag\tfit_frag\tfit_exposure\tn_blocks\tn_projectable\tn_determined\tn_exact\tn_wrong\truntime_s\tpeak_rss_bytes\tstatus\n' > "$RUNS"
printf 'commit\tlocus\tdonor\tarm\tblock\tkind\tbubble_id\tn_alleles\ttruth_a1\ttruth_a2\tcalled_a1\tcalled_a2\tprojectable\tdetermined\tequiv_set_size\tstatus\n' > "$BLOCKS"

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
# Write one named record out of a scored_sequences.fa, uppercased.
extract_fa() {   # extract_fa <fa> <name> <out>
  "$PY" - "$1" "$2" "$3" <<'PY'
import sys
fa,want,out=sys.argv[1:4]
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
open(out,'w').write(">s\n"+hit.upper()+"\n")
PY
}
# Exact edit distance between two one-record FASTAs, via the SAME command every other distance uses.
# Escalates the band until the distance is EXACT, and reports the band it needed. A distance that
# fell outside the band is not a large distance, it is no distance -- and silently treating the cap
# as the answer would make a badly wrong call look like a near miss. USED_BAND is set as a side
# effect so the row can record what the number actually cost.
USED_BAND=""
dist() {   # dist <a.fa> <b.fa> -> integer, or empty
  local v b
  for b in $BAND 16384 65536 262144; do
    v=$("$BIN" genotype-frag --exact-distance "$1" "$2" --distance-band "$b" 2>/dev/null | tr -d ' ')
    case "$v" in ''|*[!0-9]*) continue ;; esac
    if [ -z "$USED_BAND" ] || [ "$b" -gt "$USED_BAND" ]; then USED_BAND="$b"; fi
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
      say "REFUSE $LOCUS/$DONOR: donor does not have two homologues in the panel"; fails=$((fails+1)); continue
    fi
    DD="$LD/$DONOR"; mkdir -p "$DD"
    extract_fa "$LD/all.scored_sequences.fa" "$T1" "$DD/t1.fa" || { say "REFUSE $LOCUS/$DONOR: cannot spell truth 1"; fails=$((fails+1)); continue; }
    extract_fa "$LD/all.scored_sequences.fa" "$T2" "$DD/t2.fa" || { say "REFUSE $LOCUS/$DONOR: cannot spell truth 2"; fails=$((fails+1)); continue; }
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
      FLOORLINE=$("$PY" - "$AD/cat.scored_sequences.fa" "$AD" "$BIN" "$BAND" "$DD/t1.fa" "$DD/t2.fa" <<'PY'
import sys, subprocess
fa,ad,BIN,band,t1,t2=sys.argv[1:7]
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
    open(f"{ad}/c.fa","w").write(">c\n"+s.upper()+"\n")
    for h,t in ((0,t1),(1,t2)):
        r=subprocess.run([BIN,"genotype-frag","--exact-distance",f"{ad}/c.fa",t,
                          "--distance-band",band],capture_output=True,text=True)
        v=r.stdout.strip()
        if v.isdigit() and int(v)<best[h]: best[h]=int(v); who[h]=nm
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
open(o1,'w').write(">c1\n"+recs[0].upper()+"\n"); open(o2,'w').write(">c2\n"+recs[1].upper()+"\n")
PY
      [ -s "$AD/c1.fa" ] && [ -s "$AD/c2.fa" ] || { say "REFUSE $LOCUS/$DONOR/$ARM: called pair not spelled as two records"; fails=$((fails+1)); continue; }
      # best of the two orientation assignments -- the pair is unordered
      USED_BAND=""
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

      # NORMALISED FIT: NOT COMPUTED, deliberately.
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
             >> "$BLOCKS" 2>"$AD/bstat.txt" <<'PY'
import sys, hashlib
tsv,fasta,armcall,commit,locus,donor,arm,eqs,t1,t2,c1,c2=sys.argv[1:13]
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
nb=nproj=ndet=nex=nwr=0
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
    sys.stdout.write('\t'.join([commit,locus,donor,arm,str(b),
        a.get('kind',r1.get('kind','')),r1.get('bubble_id',''),a.get('n_alleles',''),
        r1['block_md5'][:8],r2['block_md5'][:8],q1['block_md5'][:8],q2['block_md5'][:8],
        pj,dt,eqs,st])+'\n')
sys.stderr.write("%d %d %d %d %d\n" % (nb,nproj,ndet,nex,nwr))
PY
      if [ ! -s "$AD/bstat.txt" ]; then
        say "REFUSE $LOCUS/$DONOR/$ARM: block projection produced no status"; fails=$((fails+1)); continue
      fi
      read -r NB NPROJ NDET NEX NWR <<<"$(head -1 "$AD/bstat.txt")"
      # An aborted projection writes 0 0 0 0 0 and a reason. It must REFUSE the donor, not leave the
      # block columns blank and carry on -- a run row with no blocks would still be averaged.
      if [ "${NB:-0}" = 0 ]; then
        say "REFUSE $LOCUS/$DONOR/$ARM: $(tail -1 "$AD/bstat.txt")"; fails=$((fails+1)); continue
      fi

      printf '%s\t%s\t%s\t%s\t%s\tsim\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\texact(band=%s)\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tOK\n' \
        "$COMMIT" "$BIN_MD5" "$LOCUS" "$DONOR" "$ARM" "$SEED" "$NPAIRS" \
        "$GFA_MD5" "$BUB_MD5" "$TRUTH_MD5" "$READS_MD5" "$CAT_MD5" "${SHORTLIST_MD5:-}" \
        "$T1" "$T2" "$PANEL_FLOOR" "$FSURV" "$C1" "$C2" "$CALLED_D" "$EXCESS" "${USED_BAND:-$BAND}" \
        "${EQS:-}" "$FIT" "$FITF" "$FITE" \
        "${NB:-}" "${NPROJ:-}" "${NDET:-}" "${NEX:-}" "${NWR:-}" \
        "${RUNTIME:-}" "${PEAK:-}" >> "$RUNS"
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
