#!/usr/bin/env bash
# genotype_bounded_search.sh - acceptance gates for the BOUNDED-COMPLETE placement search.
#
#   genotype_bounded_search.sh <panvar-binary> <out-dir>
#
# WRITTEN BEFORE THE IMPLEMENTATION, deliberately. Every acceleration on this branch that was
# measured only after it existed turned out to be measuring something else: recruitment labels were
# not valid factor scopes, a projector assumed source->sink order, a "projection" of a named pair
# returned the called pair. The conditions below are frozen first so the search is built against
# them rather than described by them afterwards.
#
# WHY IT IS NEEDED, concretely: --reference-score enumerates every fragment start on both
# haplotypes and does not finish in ten minutes on c4's 226kb/252kb pair even for 76 fragments. The
# exact audit of the 76 fragments that decide c4/NA19240 is blocked on this, as is any trustworthy
# fragment linkage factor.
#
# THE CONTRACT:
#   1. COMPLETENESS. Every mate placement within the declared edit band is found. No occurrence cap,
#      no top-k. A cap makes completeness a function of a tuning parameter.
#   2. STATE FIDELITY. Strand, start, end, insert length and haplotype identity are preserved. A
#      placement is a state, not a position: the fragment model's state is (start, insert,
#      orientation) and collapsing any component changes the likelihood.
#   3. VALID FR ONLY, which is three conditions and not one: forward/reverse orientation, the
#      reverse mate DOWNSTREAM of the forward one, and an insert length inside the prior's support.
#      Orientation alone admits pairs the model assigns zero probability.
#   4. MULTIPLICITY PRESERVED MATHEMATICALLY, not as a record count.
#        * the same placement discovered through several anchors is ONE state;
#        * different repeat-copy coordinates are DISTINCT origins;
#        * identical-likelihood origins MAY be compressed to one representative plus
#          log(multiplicity), provided placement mass and exposure both stay exact.
#      Requiring a particular number of stored records would forbid precisely the compression LPA
#      needs. The gate compares LIKELIHOOD and EXPOSURE against exhaustive enumeration, never
#      storage. Collapsing distinct origins destroys copy-number information in a marginal placement
#      model and must change neither placement mass nor exposure.
#   5. EXHAUSTIVE AGREEMENT on small fixtures: substitution, indel, reverse strand, repeat, junction,
#      and unseeded -- meaning no usable CURRENT syncmer seed. Not "no exact seed anywhere": under
#      the pigeonhole construction a read within d edits must contain an exact piece, ambiguous
#      bases aside, so a fixture claiming otherwise would be testing an impossibility.
#   6. REPORTED WORK: candidate starts, verified starts, placements kept, runtime, and the
#      omitted-mass bound. A search that cannot say what it skipped cannot be trusted not to skip.
#
# THE ACCEPTANCE EQUALITY, which subsumes the rest:
#
#     bounded likelihood mass == exhaustive likelihood mass
#     bounded exposure        == exhaustive exposure
#
# for EVERY candidate and every diploid pair on the small fixtures, with exact multiplicity-aware
# compression permitted internally. Anything that changes either quantity is a different model, not
# a faster one.
#
# STATUS: the search is NOT IMPLEMENTED. This file exits 77 (skip) until the entry point exists, so
# it can be registered now and start gating the moment there is something to gate. It must never be
# made to pass by weakening a condition.
set -uo pipefail
BIN="${1:?usage: genotype_bounded_search.sh <panvar> <outdir>}"
OUT="${2:?}"; mkdir -p "$OUT"

if ! "$BIN" genotype-frag --help 2>&1 | grep -q -- "--bounded-search"; then
  echo "SKIP: --bounded-search is not implemented yet; gates are frozen and waiting"
  echo "      acceptance: bounded mass == exhaustive mass AND bounded exposure == exhaustive"
  echo "                  exposure, per candidate and per diploid pair, on 6 fixture classes"
  echo "                  (substitution, indel, reverse strand, repeat, junction, unseeded);"
  echo "                  completeness within band, no cap/top-k, valid-FR join, exact"
  echo "                  multiplicity-aware compression allowed, work and omitted-mass reported"
  exit 77
fi

echo "FAIL: --bounded-search exists but its gates have not been written."
echo "      Implementing the search without implementing these checks is the failure mode this"
echo "      file exists to prevent. Write the assertions, do not delete the file."
exit 1
