#!/usr/bin/env python3
"""Draw the held-out C4 synthetic-mosaic plans by a deterministic rule.

The manifest this prints is committed as tests/genotype_c4_heldout_mosaic_manifest.tsv and is frozen
BEFORE anything is simulated. This script exists so the draw is a rule rather than a choice: anyone can
re-run it and get the committed manifest back, byte for byte, and the harness gate does exactly that.

The rule, in full:

  pool        every donor whose role in tests/genotype_development_donors.tsv is held_out, that has
              exactly two C4 panel rows, and whose BOTH rows are forward-oriented. Orientation comes
              from `genotype_probe orientations`, never from a re-derivation here: mosaic-truth refuses
              to splice a reverse walk against a forward one, so a reverse donor cannot source a mosaic
              at all and must be excluded before the draw rather than discovered when staging throws.
              The C4 REFERENCE sample is excluded on top of that -- using it as a mosaic source would
              make one homologue the chain frame itself, which is a different experiment.
  order       ascending sha256(salt + ":" + donor).
  pairs       consecutive disjoint pairs of that order: (0,1), (2,3), ... Ten pairs are used; any
              remainder is named in the manifest and left unused, so the pool size does not silently
              decide how many plans there are.
  classes     five frozen classes, two pairs each, in the fixed order below.
  haplotypes  swap_hap_b is bit 0 of sha256(salt + ":pairing:" + a + "x" + b). It decides which of
              donor B's rows pairs with donor A's row h, so a plan varies the donor/HAPLOTYPE
              combination and not only the switch points.
  switches    drawn from sha256(salt + ":" + class + ":" + a + "x" + b + ":" + counter), taking decimal
              digits as block indices in 1..10 until enough distinct ones are collected. Block 0 is
              never a switch point: flipping at the first block yields no switch, only a different
              source. The counter advances until the drawn plan is ADMISSIBLE, and the first admissible
              counter is the plan.

  admissible  two conditions, both checked here against the frozen panel allele matrix rather than
              discovered when staging throws:

              1. no emitted homologue spec repeats a DEVELOPMENT switch plan. The check is per
                 homologue, on the spec that is actually emitted -- checking the pooled set of drawn
                 points instead once let both_homologues emit hap1="5", which is exactly the
                 development one_switch plan.

              2. the homologues are on-panel exactly as often as the class requires: 0 of 2 for the
                 classes that recombine both, 1 of 2 for a mosaic beside an intact homologue, 2 of 2
                 for the controls. This is the same test mosaic-truth enforces, applied to the same
                 array (panel().template_alleles), so a plan that would be refused at staging is never
                 committed. It cannot be skipped: a switch point next to a chain end takes only a block
                 or two from donor A, and where donor B's row carries the same alleles there, the
                 "mosaic" IS donor B's row. That is a property of the data, invisible to a draw over
                 donor names, and it is what the first version of this rule kept proposing.

Development exclusions, which this rule enforces rather than assumes: the development donor pair
(HG02666 x NA19240) and every other development donor are outside the pool because their role is not
held_out, and the development switch plans {5}, {3,5}, {8}, {3,6,9} are rejected by condition 1.

Usage:
  scripts/draw_heldout_mosaic_plans.py <orientation.tsv> <donor manifest.tsv> <panel alleles.tsv>
"""

import collections
import hashlib
import sys

SALT = "panvar-genotype-mosaic-heldout-c4-2026-09-21"
REFERENCE_SAMPLE = "HG00096"
PAIRS_USED = 10
CLASSES = [
    "one_homologue_once",
    "both_homologues",
    "three_switches",
    "dense_switching",
    "panel_control",
]
# The five development plans, as switch-point SETS. "every" is not listed: the dense class below draws
# eight points rather than alternating at every block, so it cannot collide with it.
DEVELOPMENT_SWITCH_SETS = {
    frozenset([5]),
    frozenset([3, 5]),
    frozenset([8]),
    frozenset([3, 6, 9]),
}
BLOCKS = 11            # C4: ctx:0 sv:1 ctx:2 sv:3 ctx:4 sv:5 ctx:6 sv:7 ctx:8 sv:9 ctx:10
DENSE_POINTS = 8       # dense, but deliberately not the development "every" plan
# How many of the two homologues a class requires to be panel rows.
ON_PANEL_REQUIRED = {
    "one_homologue_once": 1,
    "both_homologues": 0,
    "three_switches": 0,
    "dense_switching": 0,
    "panel_control": 2,
}
OFF_PANEL_DECLARED = {
    "one_homologue_once": "one",
    "both_homologues": "both",
    "three_switches": "both",
    "dense_switching": "both",
    "panel_control": "none",
}


def digest(text):
    return hashlib.sha256((SALT + ":" + text).encode()).hexdigest()


def load_panel(path):
    """The panel allele matrix, in panel order -- the order mosaic-truth's rows_of() sees."""
    names, alleles = [], []
    for line in open(path):
        if line.startswith("template\t") or line.startswith("#"):
            continue
        name, spec = line.rstrip("\n").split("\t")
        names.append(name)
        alleles.append(tuple(int(value) for value in spec.split(",")))
    return names, alleles


def eligible_pool(orientation_path, manifest_path):
    orientations = collections.defaultdict(list)
    for line in open(orientation_path):
        if line.startswith("template\t") or line.startswith("#"):
            continue
        name, orientation = line.rstrip("\n").split("\t")
        orientations[name.split("#")[0]].append(orientation)
    role = {}
    for line in open(manifest_path):
        if line.startswith("#") or line.startswith("donor\t"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 2:
            role[fields[0]] = fields[1]
    pool = [
        donor
        for donor, rows in orientations.items()
        if role.get(donor) == "held_out"
        and len(rows) == 2
        and set(rows) == {"forward"}
        and donor != REFERENCE_SAMPLE
    ]
    return sorted(pool, key=digest)


def rows_of(names, donor):
    rows = [index for index, name in enumerate(names) if name.startswith(donor + "#")]
    if len(rows) != 2:
        raise SystemExit("donor %s has %d panel rows, expected 2" % (donor, len(rows)))
    return rows


def donor_of_block(spec, blocks):
    """Which donor (0 or 1) supplies each block, from one homologue's switch spec."""
    if spec == "every":
        return [block % 2 for block in range(blocks)]
    plan, current = [0] * blocks, 0
    if spec:
        flips = {int(point) for point in spec.split(",") if point != ""}
        for block in range(blocks):
            if block in flips:
                current = 1 - current
            plan[block] = current
    return plan


def haplotype_paths(names, alleles, donor_a, donor_b, swap, specs):
    """Both homologues' allele paths, built exactly as mosaic-truth builds them."""
    blocks = len(alleles[0])
    rows_a, rows_b = rows_of(names, donor_a), rows_of(names, donor_b)
    paths = []
    for hap in (0, 1):
        plan = donor_of_block(specs[hap], blocks)
        row_of = lambda donor, hap=hap: (
            rows_a[hap] if donor == 0 else rows_b[1 - hap if swap else hap]
        )
        paths.append(tuple(alleles[row_of(plan[block])][block] for block in range(blocks)))
    return paths


def draw_points(tag, count, counter):
    """Distinct block indices in 1..10 from one counter. Admissibility is decided by the caller."""
    value = int(digest("%s:%d" % (tag, counter)), 16)
    points = []
    while len(points) < count and value > 0:
        candidate = 1 + (value % (BLOCKS - 1))
        value //= (BLOCKS - 1)
        if candidate not in points:
            points.append(candidate)
    return sorted(points) if len(points) == count else None


def spec_set(spec):
    return frozenset(int(point) for point in spec.split(",") if point != "")


def candidate_specs(klass, tag, counter):
    if klass == "one_homologue_once":
        points = draw_points(tag, 1, counter)
        return (",".join(map(str, points)), "") if points else None
    if klass == "both_homologues":
        points = draw_points(tag, 2, counter)
        return (str(points[0]), str(points[1])) if points else None
    if klass == "three_switches":
        points = draw_points(tag, 3, counter)
        return ((",".join(map(str, points)),) * 2) if points else None
    if klass == "dense_switching":
        points = draw_points(tag, DENSE_POINTS, counter)
        return ((",".join(map(str, points)),) * 2) if points else None
    raise SystemExit("unknown class " + klass)


def class_specs(klass, tag, names, alleles, donor_a, donor_b, swap):
    """The first admissible draw for this class and donor pair."""
    panel = set(alleles)
    required = ON_PANEL_REQUIRED[klass]
    if klass == "panel_control":
        specs = ("", "")
    else:
        specs = None
        for counter in range(1000):
            candidate = candidate_specs(klass, tag, counter)
            if not candidate:
                continue
            if any(spec_set(spec) in DEVELOPMENT_SWITCH_SETS for spec in candidate if spec != ""):
                continue
            paths = haplotype_paths(names, alleles, donor_a, donor_b, swap, candidate)
            if sum(1 for path in paths if path in panel) != required:
                continue
            specs = candidate
            break
        if specs is None:
            raise SystemExit("no admissible draw for " + tag)
    # The controls are drawn from nothing, but still have to BE controls.
    paths = haplotype_paths(names, alleles, donor_a, donor_b, swap, specs)
    if sum(1 for path in paths if path in panel) != required:
        raise SystemExit("%s: %d of 2 homologues on panel, class requires %d" %
                         (tag, sum(1 for path in paths if path in panel), required))
    return specs


def plans(pool, names, alleles):
    rows = []
    for index in range(PAIRS_USED):
        donor_a, donor_b = pool[2 * index], pool[2 * index + 1]
        klass = CLASSES[index // 2]
        swap = 0 if klass == "panel_control" else \
            int(digest("pairing:%sx%s" % (donor_a, donor_b)), 16) & 1
        hap1, hap2 = class_specs(
            klass, "%s:%sx%s" % (klass, donor_a, donor_b), names, alleles, donor_a, donor_b, swap
        )
        rows.append(
            ["plan", "%s_%d" % (klass, index % 2 + 1), klass, donor_a, donor_b,
             str(swap), hap1, hap2, OFF_PANEL_DECLARED[klass]]
        )
    return rows


def main():
    if len(sys.argv) != 4:
        raise SystemExit(__doc__.strip().splitlines()[-1].strip())
    pool = eligible_pool(sys.argv[1], sys.argv[2])
    names, alleles = load_panel(sys.argv[3])
    if len(pool) < 2 * PAIRS_USED:
        raise SystemExit("pool of %d cannot supply %d disjoint pairs" % (len(pool), PAIRS_USED))
    print("# salt\t%s" % SALT)
    print("# pool_size\t%d" % len(pool))
    print("# pool\t%s" % ",".join(pool))
    print("# unused_remainder\t%s" % (",".join(pool[2 * PAIRS_USED:]) or "none"))
    print("plan\tname\tclass\tdonor_a\tdonor_b\tswap_hap_b\tswitches_hap1\tswitches_hap2\trequire_off_panel")
    for row in plans(pool, names, alleles):
        print("\t".join(row))


if __name__ == "__main__":
    main()
