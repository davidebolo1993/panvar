#!/usr/bin/env python3
"""M1-D: representation reachability of truth's sv:7 units from the FROZEN panel-only inventory (diagnostic evaluator).

  reach.py pairs  --inventory PREFIX --frozen FROZEN.tsv --truth TRUTH.fa --alleles ALLELES.fa --out OUT_PREFIX
      Verifies the inventory against its committed freeze record, partitions each truth block with the inventory's own
      phase-A anchor (the same rule as the inventory; the anchor was chosen from panel bytes only), and writes the pairs
      for an exhaustive edit-distance pass: every truth component against every unique inventory component of the same
      type, and every truth block against every carried whole allele.
  reach.py select --inventory PREFIX --truth TRUTH.fa --alleles ALLELES.fa --distances D.tsv --components TRUTH_COMPONENTS.fa
                  --out OUT_PREFIX
      Reports exact availability and the nearest component per truth component (ties: smallest distance, then record
      name), the nearest whole allele per truth block, and composes, from existing components only (no new bases, no
      consensus, no within-unit patching):
        compose.truthH       nearest prefix + nearest period(s) + nearest tail, joined at the anchor
        periodonly.truthH    the nearest whole allele W(H) with only its period(s) replaced by the nearest period(s)
                             (only when W(H) has truth H's unit count)
        whole.truthH         W(H) itself
      A component-wise distance sum bounds a composition's distance; its actual distance is measured later. Nearest
      components establish nearest-unit distances, not the best whole-diplotype likelihood.
"""
import argparse
import csv
import hashlib
import sys

K = 31


def fasta(path):
    out, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name:
            out[name].append(line)
    return {k: "".join(v) for k, v in out.items()}


def sha_file(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def positions(s, k):
    out, i = [], s.find(k)
    while i >= 0:
        out.append(i)
        i = s.find(k, i + 1)
    return out


def partition(seq, anchor):
    p = positions(seq, anchor)
    assert p, "the inventory anchor does not occur in a truth block"
    parts = [("prefix", seq[:p[0]])] + [(f"period{i + 1}", seq[p[i]:p[i + 1]]) for i in range(len(p) - 1)] + [("tail", seq[p[-1]:])]
    assert "".join(x for _, x in parts) == seq
    return parts, p


def kind(part):
    return "period" if part.startswith("period") else part


def pairs(a):
    frozen = dict(line.rstrip("\n").split("\t", 1) for line in open(a.frozen) if "\t" in line)
    for suffix in (".tsv", ".fa", ".anchors.tsv"):
        got = sha_file(a.inventory + suffix)
        assert frozen.get("sha256" + suffix) == got, f"inventory{suffix} differs from its freeze record"
    anchor = next(r["anchor"] for r in csv.DictReader(open(a.inventory + ".anchors.tsv"), delimiter="\t") if r["phase"] == "A")
    inv = fasta(a.inventory + ".fa")
    truth = fasta(a.truth)
    carried = sorted({r["source_allele"] for r in csv.DictReader(open(a.inventory + ".tsv"), delimiter="\t")},
                     key=lambda n: int(n.split("_")[1]))
    out, rows = [], []
    with open(a.out + ".truth_components.fa", "w") as fa:
        for t in sorted(truth):
            parts, p = partition(truth[t], anchor)
            rows.append(f"{t}\t{len(p)}\t{','.join(map(str, p))}\t" + ",".join(f"{n}:{len(x)}" for n, x in parts))
            for part, piece in parts:
                fa.write(f">{t}.{part}\n{piece}\n")
                for rec in inv:
                    if rec.startswith(kind(part) + "_"):
                        out.append(f"comp:{t}.{part}:{rec}\t{a.out}.truth_components.fa\t{t}.{part}\t{a.inventory}.fa\t{rec}")
            for allele in carried:
                out.append(f"whole:{t}:{allele}\t{a.truth}\t{t}\t{a.alleles}\t{allele}")
    open(a.out + ".pairs.tsv", "w").write("\n".join(out) + "\n")
    open(a.out + ".truth_units.tsv", "w").write("truth\tunits\tanchor_positions\tcomponents\n" + "\n".join(rows) + "\n")
    print(open(a.out + ".truth_units.tsv").read(), end="")
    print(f"{len(out)} pairs")


def select(a):
    inv = fasta(a.inventory + ".fa")
    alleles = fasta(a.alleles)
    truth = fasta(a.truth)
    comps = fasta(a.components)
    meta = list(csv.DictReader(open(a.inventory + ".tsv"), delimiter="\t"))
    units_of = {r["source_allele"]: int(r["source_units"]) for r in meta}
    parts_of = {}
    for r in meta:
        parts_of.setdefault(r["source_allele"], []).append(r)
    sources = {}
    for r in meta:
        sources.setdefault(r["unique_record"], []).append(r["component"])
    dist = list(csv.DictReader(open(a.distances), delimiter="\t"))
    best, whole = {}, {}
    for r in dist:
        label, d = r["label"], int(r["distance"])
        if label.startswith("comp:"):
            tc, rec = label.split(":")[1], label.split(":")[2]
            best.setdefault(tc, []).append((d, rec))
        else:
            t, allele = label.split(":")[1], label.split(":")[2]
            whole.setdefault(t, []).append((d, int(allele.split("_")[1]), allele))
    report = ["truth_component\tbp\texact_available\tnearest_record\tnearest_distance\ttied_at_nearest\tsecond_distance\t"
              "nearest_sources"]
    nearest = {}
    for tc in sorted(best):
        ranked = sorted(best[tc])
        d0, rec = ranked[0]
        ties = sum(1 for d, _ in ranked if d == d0)
        second = next((d for d, _ in ranked if d > d0), "NA")
        nearest[tc] = (rec, d0)
        report.append(f"{tc}\t{len(comps[tc])}\t{d0 == 0}\t{rec}\t{d0}\t{ties}\t{second}\t"
                      f"{';'.join(sources[rec][:6])}{' ...' if len(sources[rec]) > 6 else ''}")
    composed, structure = {}, []
    lines = ["truth\tunits\tnearest_whole\twhole_distance\twhole_units\tcomponent_distance_sum\tcomposition"]
    for t in sorted(truth):
        names = sorted([tc for tc in nearest if tc.startswith(t + ".")],
                       key=lambda n: (0 if n.endswith("prefix") else 2 if n.endswith("tail") else 1, n))
        n_units = sum(1 for n in names if ".period" in n) + 1
        pieces = [nearest[n][0] for n in names]
        composed[f"compose.{t}"] = "".join(inv[p] for p in pieces)
        structure.append(f"compose.{t}\t" + "+".join(pieces))
        wd, _, w = sorted(whole[t])[0]
        composed[f"whole.{t}"] = alleles[w]
        structure.append(f"whole.{t}\t{w}")
        if units_of[w] == n_units:
            wp = parts_of[w]
            pre = next(inv[r["unique_record"]] for r in wp if r["type"] == "prefix")
            tail = next(inv[r["unique_record"]] for r in wp if r["type"] == "tail")
            periods = [inv[nearest[n][0]] for n in names if ".period" in n]
            composed[f"periodonly.{t}"] = pre + "".join(periods) + tail
            structure.append(f"periodonly.{t}\t{w}.prefix+" + "+".join(nearest[n][0] for n in names if ".period" in n) + f"+{w}.tail")
        lines.append(f"{t}\t{n_units}\t{w}\t{wd}\t{units_of[w]}\t{sum(nearest[n][1] for n in names)}\t{'+'.join(pieces)}")
    with open(a.out + ".composed.fa", "w") as f:
        for name, seq in composed.items():
            f.write(f">{name}\n{seq}\n")
    open(a.out + ".composed_structure.tsv", "w").write("record\tcomponents\n" + "\n".join(structure) + "\n")
    open(a.out + ".reachability.tsv", "w").write("\n".join(report) + "\n")
    open(a.out + ".blocks.tsv", "w").write("\n".join(lines) + "\n")
    print("\n".join(report))
    print("\n".join(lines))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=("pairs", "select"))
    for x in ("inventory", "truth", "alleles", "out"):
        ap.add_argument("--" + x, required=True)
    ap.add_argument("--frozen")
    ap.add_argument("--distances")
    ap.add_argument("--components")
    args = ap.parse_args()
    sys.exit(pairs(args) if args.mode == "pairs" else select(args))
