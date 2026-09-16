#!/usr/bin/env bash
# Gate: genotype_probe subsample produces NESTED lower-coverage read sets from one simulation.
#
#   tests/genotype_subsample_nested.sh <build dir> <scratch dir>
#
# Every fragment kept at a coverage must be kept at every higher coverage, the choice must depend only on the salt
# and the fragment name (so a rerun is byte-identical and a different source coverage with the same ratio keeps the
# same fragments), two different salts must give different masks while each stays nested, mates must stay paired
# and in order, the kept fraction must be close to the requested one, and the declared fragment-start intensity
# must follow the coverage.
set -euo pipefail
BUILD="${1:?usage: tests/genotype_subsample_nested.sh <build dir> <scratch dir>}"
SCRATCH="${2:?usage: tests/genotype_subsample_nested.sh <build dir> <scratch dir>}"
PROBE="$BUILD/genotype_probe"
rm -rf "$SCRATCH"; mkdir -p "$SCRATCH/source"

python3 - "$SCRATCH/source" <<'FIXTURE'
import random, sys
out = sys.argv[1]
rng = random.Random(20260916)
bases = "ACGT"
with open(out + "/reads_1.fq", "w") as r1, open(out + "/reads_2.fq", "w") as r2:
    for i in range(1, 4001):
        left = "".join(rng.choice(bases) for _ in range(100))
        right = "".join(rng.choice(bases) for _ in range(100))
        r1.write("@probe_%d/1\n%s\n+\n%s\n" % (i, left, "I" * 100))
        r2.write("@probe_%d/2\n%s\n+\n%s\n" % (i, right, "I" * 100))
with open(out + "/library.txt", "w") as lib:
    lib.write("# fixture\ninsert_mean_bp 350\ninsert_sd_bp 40\nsupport_sigmas 4\n"
              "substitution_error 0.002\noutlier_weight 0.001\nfragment_start_intensity 0.08\n")
FIXTURE

for coverage in 30 20 15 10 5 2; do
  mkdir -p "$SCRATCH/cov$coverage"
  "$PROBE" subsample --reads1 "$SCRATCH/source/reads_1.fq" --reads2 "$SCRATCH/source/reads_2.fq" \
    --library "$SCRATCH/source/library.txt" --source-coverage 30 --coverage "$coverage" \
    --out "$SCRATCH/cov$coverage" > /dev/null
done
# Determinism, and independence from how the ratio is expressed: 10/30 must keep exactly what 1/3 keeps.
mkdir -p "$SCRATCH/repeat" "$SCRATCH/ratio"
"$PROBE" subsample --reads1 "$SCRATCH/source/reads_1.fq" --reads2 "$SCRATCH/source/reads_2.fq" \
  --library "$SCRATCH/source/library.txt" --source-coverage 30 --coverage 10 --out "$SCRATCH/repeat" > /dev/null
"$PROBE" subsample --reads1 "$SCRATCH/source/reads_1.fq" --reads2 "$SCRATCH/source/reads_2.fq" \
  --library "$SCRATCH/source/library.txt" --source-coverage 3 --coverage 1 --out "$SCRATCH/ratio" > /dev/null
# Two donors, same source reads: different salts must give different masks, each nested within itself.
for salt in donorA donorB; do
  for coverage in 10 5; do
    mkdir -p "$SCRATCH/$salt.cov$coverage"
    "$PROBE" subsample --reads1 "$SCRATCH/source/reads_1.fq" --reads2 "$SCRATCH/source/reads_2.fq" \
      --library "$SCRATCH/source/library.txt" --source-coverage 30 --coverage "$coverage" --name-salt "$salt" \
      --out "$SCRATCH/$salt.cov$coverage" > /dev/null
  done
done

python3 - "$SCRATCH" <<'CHECK'
import os, sys
scratch = sys.argv[1]
coverages = [30, 20, 15, 10, 5, 2]

def names(path):
    out = []
    with open(path) as handle:
        for index, line in enumerate(handle):
            if index % 4 == 0:
                out.append(line.strip()[1:])
    return out

def fail(message):
    sys.exit("genotype subsample nesting: FAIL: " + message)

kept = {}
for coverage in coverages:
    first = names("%s/cov%d/reads_1.fq" % (scratch, coverage))
    second = names("%s/cov%d/reads_2.fq" % (scratch, coverage))
    if [n[:-2] for n in first] != [n[:-2] for n in second]:
        fail("mates are not paired and in order at %dx" % coverage)
    if any(not n.endswith("/1") for n in first) or any(not n.endswith("/2") for n in second):
        fail("mate suffixes are not preserved at %dx" % coverage)
    kept[coverage] = [n[:-2] for n in first]

source = names(scratch + "/source/reads_1.fq")
if kept[30] != [n[:-2] for n in source]:
    fail("subsampling to the source coverage must keep every fragment in order")
for low, high in zip(coverages[1:], coverages[:-1]):
    if not set(kept[low]).issubset(set(kept[high])):
        fail("the %dx set is not a subset of the %dx set" % (low, high))
    if kept[low] != [n for n in kept[high] if n in set(kept[low])]:
        fail("the %dx set does not preserve the order of the %dx set" % (low, high))
    if len(kept[low]) >= len(kept[high]):
        fail("the %dx set is not smaller than the %dx set" % (low, high))

total = len(source)
for coverage in coverages:
    expected = total * coverage / 30.0
    if abs(len(kept[coverage]) - expected) > 4.0 * (expected ** 0.5) + 1.0:
        fail("the %dx set keeps %d fragments, far from the requested %.0f" % (coverage, len(kept[coverage]), expected))

if kept[10] != [n[:-2] for n in names(scratch + "/repeat/reads_1.fq")]:
    fail("a repeated subsample is not identical")
if kept[10] != [n[:-2] for n in names(scratch + "/ratio/reads_1.fq")]:
    fail("the same ratio expressed differently keeps different fragments")

salted = {}
for salt in ("donorA", "donorB"):
    for coverage in (10, 5):
        salted[(salt, coverage)] = [n[:-2] for n in names("%s/%s.cov%d/reads_1.fq" % (scratch, salt, coverage))]
    if not set(salted[(salt, 5)]).issubset(set(salted[(salt, 10)])):
        fail("the %s 5x set is not nested inside its own 10x set" % salt)
for coverage in (10, 5):
    a, b = salted[("donorA", coverage)], salted[("donorB", coverage)]
    if a == b:
        fail("two salts produced the same %dx mask" % coverage)
    shared = len(set(a) & set(b)) / max(1, len(set(a) | set(b)))
    if shared > 0.5:
        fail("two salts share %.0f%% of the %dx mask; they should be independent" % (100 * shared, coverage))
    if a == kept[coverage] or b == kept[coverage]:
        fail("a salted %dx mask equals the unsalted one" % coverage)

def lambda_of(path):
    for line in open(path):
        if line.startswith("fragment_start_intensity"):
            return float(line.split()[1])
    fail("no fragment_start_intensity in " + path)

for coverage in coverages:
    got = lambda_of("%s/cov%d/library.txt" % (scratch, coverage))
    want = 0.08 * coverage / 30.0
    if abs(got - want) > 1e-12:
        fail("the declared intensity at %dx is %g, not %g" % (coverage, got, want))

print("genotype subsample nesting: PASS (%s kept of %d, nested and deterministic; two salts differ, each nested)"
      % ("/".join(str(len(kept[c])) for c in coverages), total))
CHECK
