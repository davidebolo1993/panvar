#!/usr/bin/env python3
"""M1-D: the panel-only repeat-unit inventory of C4 sv:7 (no truth, no NA19650 path, no other graph sequence).

Source: `genotype_probe block-sequences` run WITHOUT --truth-gfa on the LOO model (NA19650 paths are absent from it), i.e.
the sv:7 allele bytes of the active panel in the model frame, with the templates carrying each allele. Only alleles
carried by at least one panel template are used.

Unit boundaries (declared before execution; panel bytes only):
  * S = the shortest carried allele of length >= 10,000 bp (ties: lowest allele index).
  * An anchor is a 31-mer starting in S that occurs exactly once in S, occurs at least once in EVERY carried allele, has
    no N, and whose consecutive occurrences in every carried allele are >= 20,000 bp apart.
  * Phase A = the qualifying anchor with the smallest start in S; phase B = the first qualifying anchor starting at least
    len(S)/2 after phase A's. Boundaries are RELIABLE only if both phases give every carried allele the same unit count
    and every phase-A spacing lies in [20,000, 40,000] bp; otherwise the limitation is reported and no inventory is used.
  * With phase-A anchor starts p1 < ... < pn in allele X (n = X's unit count): prefix = X[0:p1], period_i = X[p_i:p_(i+1)]
    for i < n, tail = X[p_n:]. prefix + periods + tail must reconstruct X byte for byte.
Every component records its source allele, coordinates, orientation (the model frame, as the block is spelled), length,
SHA-256 and carrying templates. Identical bytes of one component type are listed once in the FASTA with every source.

  inventory.py --sequences PREFIX --block ID --out OUT_PREFIX
"""
import argparse
import hashlib

K, MIN_SOURCE, MIN_GAP, MAX_GAP = 31, 10000, 20000, 40000


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


def sha(s):
    return hashlib.sha256(s.encode()).hexdigest()


def positions(s, k):
    out, i = [], s.find(k)
    while i >= 0:
        out.append(i)
        i = s.find(k, i + 1)
    return out


def qualifies(k, source, alleles):
    if "N" in k or len(positions(source, k)) != 1:
        return None
    found = {}
    for name, seq in alleles.items():
        p = positions(seq, k)
        if not p or any(b - a < MIN_GAP for a, b in zip(p, p[1:])):
            return None
        found[name] = p
    return found


def main():
    ap = argparse.ArgumentParser()
    for x in ("sequences", "block", "out"):
        ap.add_argument("--" + x, required=True)
    a = ap.parse_args()
    raw = fasta(f"{a.sequences}.{a.block}.alleles.fa")
    carriers = {}
    for line in open(f"{a.sequences}.block_sequences.tsv"):
        f = line.rstrip("\n").split("\t")
        if f[0] == a.block and f[1] == "panel" and f[5] != "-":
            carriers[f[3]] = f[5].split(",")
        assert not (f[0] == a.block and f[1] == "truth"), "the inventory source must not contain truth"
    assert not any("NA19650" in t for ts in carriers.values() for t in ts), "an excluded NA19650 path carries an allele"
    alleles = {n: raw[n] for n in raw if n in carriers}
    order = sorted(alleles, key=lambda n: int(n.split("_")[1]))
    source = min((n for n in order if len(alleles[n]) >= MIN_SOURCE), key=lambda n: (len(alleles[n]), int(n.split("_")[1])))
    s = alleles[source]
    phase_a = phase_b = None
    for i in range(len(s) - K + 1):
        if phase_a is not None and i < phase_a[0] + len(s) // 2:
            continue
        hit = qualifies(s[i:i + K], s, alleles)
        if hit is None:
            continue
        if phase_a is None:
            phase_a = (i, s[i:i + K], hit)
        else:
            phase_b = (i, s[i:i + K], hit)
            break
    log = [f"block\t{a.block}", f"carried_alleles\t{len(alleles)}", f"uncarried_alleles_excluded\t{len(raw) - len(alleles)}",
           f"unit_source\t{source} ({len(s)} bp; {','.join(carriers[source])})"]
    reliable = phase_a is not None and phase_b is not None
    if reliable:
        counts_a = {n: len(p) for n, p in phase_a[2].items()}
        counts_b = {n: len(p) for n, p in phase_b[2].items()}
        disagree = [n for n in order if counts_a[n] != counts_b[n]]
        bad_gap = [n for n in order for x, y in zip(phase_a[2][n], phase_a[2][n][1:]) if not MIN_GAP <= y - x <= MAX_GAP]
        reliable = not disagree and not bad_gap
        log += [f"phase_A_anchor\t{phase_a[1]} (start {phase_a[0]} in {source})",
                f"phase_B_anchor\t{phase_b[1]} (start {phase_b[0]} in {source})",
                f"unit_count_disagreements_between_phases\t{len(disagree)} {','.join(disagree)}",
                f"spacings_outside_[{MIN_GAP},{MAX_GAP}]\t{len(bad_gap)}"]
        hist = {}
        for n in order:
            hist[counts_a[n]] = hist.get(counts_a[n], 0) + 1
        log.append("unit_count_histogram\t" + ", ".join(f"{u} units: {c}" for u, c in sorted(hist.items())))
    else:
        log.append("anchors\tNOT FOUND for both phases")
    log.append(f"boundaries_reliable\t{reliable}")
    rows, unique = [], {}
    if reliable:
        for n in order:
            x, p = alleles[n], phase_a[2][n]
            parts = [("prefix", 0, p[0])] + [(f"period{i + 1}", p[i], p[i + 1]) for i in range(len(p) - 1)] + [("tail", p[-1], len(x))]
            assert "".join(x[b:e] for _, b, e in parts) == x, f"{n}: components do not reconstruct the allele"
            for part, b, e in parts:
                kind = "period" if part.startswith("period") else part
                piece = x[b:e]
                key = (kind, piece)
                if key not in unique:
                    unique[key] = f"{kind}_u{sum(1 for k in unique if k[0] == kind) + 1}"
                rows.append((f"{n}.{part}", kind, n, len(p), b, e, e - b, "model_frame_forward", sha(piece), unique[key],
                             len(carriers[n]), ",".join(carriers[n])))
        with open(a.out + ".fa", "w") as fa:
            for (kind, piece), name in unique.items():
                fa.write(f">{name}\n{piece}\n")
        log.append("unique_components\t" + ", ".join(f"{k}: {sum(1 for u in unique if u[0] == k)}" for k in ("prefix", "period", "tail")))
        log.append(f"components\t{len(rows)} (every carried allele reconstructed byte for byte)")
    with open(a.out + ".tsv", "w") as f:
        f.write("component\ttype\tsource_allele\tsource_units\tbegin\tend\tbp\torientation\tsha256\tunique_record\t"
                "carrying_templates\ttemplates\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")
    with open(a.out + ".anchors.tsv", "w") as f:
        f.write("phase\tanchor\n")
        if reliable:
            f.write(f"A\t{phase_a[1]}\nB\t{phase_b[1]}\n")
    open(a.out + ".log.tsv", "w").write("\n".join(log) + "\n")
    print("\n".join(log))


if __name__ == "__main__":
    main()
