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
#   5. EXHAUSTIVE AGREEMENT on small fixtures. "Unseeded" means no usable CURRENT syncmer seed --
#      NOT "no exact seed anywhere": under the pigeonhole construction a read within d edits must
#      contain an exact piece, ambiguous bases aside, so the latter would test an impossibility.
#
#      TWO STAGES, because the emissions differ and must not be conflated:
#
#      STAGE 1 -- HAMMING-COMPLETE. Pigeonhole recruitment plus exhaustive verification of candidate
#      starts, certified against the current reference. reference_emission is FIXED-POSITION
#      HAMMING: it counts mismatches at one offset and has no gap model (src/genotype_fragments.cpp,
#      reference_emission). It can therefore certify substitution, reverse strand, repeat, junction
#      and unseeded fixtures, and nothing else. This stage alone unblocks the c4 76-fragment audit,
#      which was run with --hamming-emission.
#
#      STAGE 2 -- INDEL-AWARE, and NOT gated by this file until it has its own oracle. With indels,
#      start, reference consumption, end coordinate and insert length all need a defined alignment
#      contract, and an exhaustive edit-distance oracle must be written against that same
#      definition. Activating an indel fixture against the Hamming reference would compare two
#      DIFFERENT EMISSIONS and call the difference an acceleration defect -- the exact
#      reference-versus-accelerated confusion this branch has already paid for twice.
#   6. REPORTED WORK: candidate starts, verified starts, placements kept, runtime, and the
#      omitted-mass bound. A search that cannot say what it skipped cannot be trusted not to skip.
#
# THE ACCEPTANCE CONDITIONS. An earlier draft demanded
#
#     bounded likelihood mass == exhaustive likelihood mass
#
# which is IMPOSSIBLE and would have made this gate unsatisfiable. The exhaustive reference
# integrates every start and every insert length and does NOT truncate at max_divergence --
# reference_emission is finite at every position, which is exactly why a union-of-spans scope came
# out locus-wide in the oracle work and had to be replaced by a counterfactual test. So out-of-band
# states carry small but NONZERO probability, and a search complete only within the band cannot
# equal the total. Requiring it would force the bounded search to evaluate the outside-band tail
# exactly, defeating its purpose, or would need a separate hard-band reference -- a DIFFERENT MODEL
# from the current exhaustive one, not a faster implementation of it.
#
# What must hold instead:
#
#   A. bounded IN-BAND mass == exhaustive IN-BAND mass        (exact, every in-band placement found)
#   B. bounded exposure     == exhaustive exposure            (exact, analytic)
#   C. the omitted mass is BOUNDED, and the bound is valid:
#
#          M_found <= M_exhaustive <= M_found + M_omitted_bound
#
#      in log space, log M_exhaustive lies in
#          [ log M_found , log_add(log M_found, log M_omitted_bound) ]
#
#   D. per fragment, per candidate AND per diploid pair -- a per-fragment bound does not imply the
#      pair bound, the same reason the scope oracle's guarantee had to be made diploid;
#   E. the exhaustive reference PAIR SCORE lies inside the reported interval;
#   F. TWO TOLERANCES, recorded separately, because they are different kinds of quantity:
#        * NUMERICAL, for the exact in-band equalities A and B -- floating-point agreement;
#        * APPROXIMATION, for the omitted mass in C -- expressed as NATS PER FRAGMENT or as an
#          omitted-mass fraction, never as one raw-likelihood threshold. A shared absolute threshold
#          is meaningless across loci whose totals differ by orders of magnitude: c4 pair scores run
#          near -1.1e6 while the deciding margin is ~500 nats.
#
#   G. CERTIFICATION USES THE INTERVALS, not the width alone:
#
#          certified(g)  iff  L(g) > max over other genotypes of U(other)
#
#      A narrow interval does not certify a winner whose interval still overlaps a competitor's.
#      Genotypes whose intervals overlap are reported as an equivalence set. At c4/NA19240 the
#      winner already changes between the pessimistic and optimistic tail treatments, so it fails
#      even the weaker same-winner test and "unresolved" is the correct output there.
#
# This is not a weaker gate. It stops the test both from demanding an impossible equality and from
# passing merely because omitted likelihood rounds to zero on fixtures small enough to hide it.
# Exact multiplicity-aware compression stays permitted internally: it must change neither the
# in-band mass nor the exposure.
#
# BUILD ORDER for stage 1, frozen with the contract so the sequence is not re-litigated mid-build:
#   1. complete SINGLE-MATE Hamming search -- split into d+1 pieces, index EVERY occurrence with no
#      cap, infer the fixed start from each occurrence, deduplicate starts, verify Hamming <= d;
#   2. both strands, and reads with no usable current syncmer seed;
#   3. combine mate placements into valid FR fragment states (orientation, reverse mate downstream,
#      insert inside the prior's support);
#   4. preserve repeat-origin multiplicity while deduplicating rediscoveries of the SAME state --
#      the same placement reached through several pieces is one state, two repeat copies are two;
#   5. exact in-band mass, analytic exposure, and the outside-band bound;
#   6. propagate those bounds through haplotypes, diploid pairs, and the background mixture -- the
#      bound must survive the mixture, since log[(1-eta)*lambda*M + eta*P_bg] is where it is used;
#   7. apply certification (G) and emit the equivalence set.
# Mutation-test every fixture BEFORE running c4. A fixture that cannot fail is not evidence, and
# five vacuous checks have already been caught on this branch by asking that question late.
#
# TWO DIAGNOSTIC QUANTITIES, kept in separate columns and never merged:
#   * BEST in-band placement -- what the current max-placement caller actually uses, and the
#     quantity that reclassifies c4's 23 competitor-only and 53 truth-only fragments;
#   * SUMMED in-band placement mass -- what the reference model requires and what any future
#     linkage factor consumes.
# They answer different questions. Sharing one column would make a max-placement result look like
# evidence about the marginal model, which is precisely the confusion that produced the retracted
# "collapsed multiplicity caused the c4 failures" claim.
#
# FIRST SCIENTIFIC RESULT, in this order: the c4 table of which of the 23 competitor-only and 53
# truth-only fragments remain EXCLUSIVE under complete in-band search, then the resulting lower and
# upper genotype intervals. Not "the search runs".
#
# STATUS: the search is NOT IMPLEMENTED. This file exits 77 (skip) until the entry point exists, so
# it can be registered now and start gating the moment there is something to gate. It must never be
# made to pass by weakening a condition.
set -uo pipefail
BIN="${1:?usage: genotype_bounded_search.sh <panvar> <outdir>}"
OUT="${2:?}"; mkdir -p "$OUT"

if ! "$BIN" genotype-frag --help 2>&1 | grep -q -- "--bounded-search"; then
  echo "SKIP: --bounded-search is not implemented yet; gates are frozen and waiting"
  echo "      acceptance: in-band mass and exposure EXACT; omitted mass validly BOUNDED, with"
  echo "                  log M_exhaustive in [log M_found, logadd(log M_found, log M_bound)],"
  echo "                  per fragment/candidate/pair; reference score inside the interval;"
  echo "                  interval width within tolerance before certifying. 6 fixture classes"
  echo "                  (substitution, indel, reverse strand, repeat, junction, unseeded);"
  echo "                  completeness within band, no cap/top-k, valid-FR join, exact"
  echo "                  multiplicity-aware compression allowed, work and omitted-mass reported"
  exit 77
fi

echo "FAIL: --bounded-search exists but its gates have not been written."
echo "      Implementing the search without implementing these checks is the failure mode this"
echo "      file exists to prevent. Write the assertions, do not delete the file."
exit 1
