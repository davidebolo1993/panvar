#!/usr/bin/env python3
"""M1-B diagnostic evaluator: repeat-unit copy number per haplotype at a structural block, from bytes (truth allowed here,
never in a caller).

Inputs: the block-sequences output of one case (<prefix>.<block>.alleles.fa, <prefix>.<block>.truth.fa,
<prefix>.block_sequences.tsv) and the historical and fresh evaluator block tables (for the called alleles).

Unit copies (declared before execution; no length-based assumption):
  * the unit source S0 is the shortest panel allele at the block of length >= 10,000 bp that some panel template carries;
  * anchors are the canonical 31-mers starting every 25 bp of S0 that occur exactly once in S0 and at least once in
    >= 90 % of the template-carried alleles (conserved);
  * an anchor is informative when its occurrence count differs between carried alleles (max count >= 2);
  * a sequence's unit copies = the median count over informative anchors; support = the share of informative anchors with
    exactly that count. Fewer than 50 informative anchors is reported LOW_INFORMATION, not a count.
  * the informative span = (last - first informative anchor position in S0) + 31: the length of S0 that varies in copy.

Total dosage vs allocation: total = copies(hap1) + copies(hap2); allocation = the unordered pair. A call with the true
total and a different pair is an allocation error; a different total is a dosage error.

  structure.py --case NAME --block ID --sequences PREFIX --hist-eval PATH --fresh-eval PATH --out OUT_PREFIX
"""
import argparse
import csv
import statistics

K, STRIDE, MIN_SOURCE_BP, CONSERVED, MIN_INFORMATIVE = 31, 25, 10000, 0.90, 50
COMP = str.maketrans("ACGTN", "TGCAN")


def canon(s):
    r = s.translate(COMP)[::-1]
    return min(s, r)


def fasta(path):
    out, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].split()[0]
            out[name] = []
        elif name:
            out[name].append(line)
    return {k: "".join(v).upper() for k, v in out.items()}


def counts(seq, wanted):
    c = dict.fromkeys(wanted, 0)
    for i in range(len(seq) - K + 1):
        k = canon(seq[i:i + K])
        if k in c:
            c[k] += 1
    return c


def eval_called(path, block):
    for r in csv.DictReader(open(path), delimiter="\t"):
        if r["block_id"] == block:
            return r
    return None


def main():
    ap = argparse.ArgumentParser()
    for a in ("case", "block", "sequences", "hist_eval", "fresh_eval", "out"):
        ap.add_argument("--" + a.replace("_", "-"), required=True)
    a = ap.parse_args()
    alleles = fasta(f"{a.sequences}.{a.block}.alleles.fa")
    truth = fasta(f"{a.sequences}.{a.block}.truth.fa")
    carried = {}
    for r in csv.DictReader(open(f"{a.sequences}.block_sequences.tsv"), delimiter="\t"):
        if r["block_id"] == a.block and r["source"] == "panel":
            carried[r["name"]] = r["templates"]
    used = [n for n in alleles if carried.get(n, "-") != "-"]
    source = min((n for n in used if len(alleles[n]) >= MIN_SOURCE_BP), key=lambda n: (len(alleles[n]), n))
    s0 = alleles[source]
    pos_of, seen = {}, {}
    for i in range(len(s0) - K + 1):
        k = canon(s0[i:i + K])
        seen[k] = seen.get(k, 0) + 1
        if i % STRIDE == 0:
            pos_of.setdefault(k, i)
    anchors = {k: p for k, p in pos_of.items() if seen[k] == 1 and "N" not in k}
    table = {n: counts(alleles[n], anchors) for n in alleles}
    table.update({n: counts(s, anchors) for n, s in truth.items()})
    conserved = [k for k in anchors if sum(table[n][k] >= 1 for n in used) >= CONSERVED * len(used)]
    informative = sorted((k for k in conserved if len({table[n][k] for n in used}) > 1 and max(table[n][k] for n in used) >= 2),
                         key=lambda k: anchors[k])
    span = (anchors[informative[-1]] - anchors[informative[0]] + K) if informative else 0

    def copies(name):
        if len(informative) < MIN_INFORMATIVE:
            return "LOW_INFORMATION", 0.0
        vals = [table[name][k] for k in informative]
        m = statistics.median_low(vals)
        return m, sum(v == m for v in vals) / len(vals)

    rows = []
    for n in list(alleles) + list(truth):
        c, s = copies(n)
        rows.append((n, "truth" if n in truth else "panel", len(alleles.get(n, truth.get(n, ""))), c, f"{s:.3f}",
                     carried.get(n, "-").count(",") + 1 if carried.get(n, "-") != "-" else 0))
    with open(f"{a.out}.units.tsv", "w") as f:
        f.write("name\tsource\tbp\tunit_copies\tsupport\tcarrying_templates\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    unit = {r[0]: r[3] for r in rows}
    t1, t2 = (unit[n] for n in sorted(truth))

    def call(path):
        r = eval_called(path, a.block)
        if r is None:
            return None
        x, y = f"allele_{r['called_allele1']}", f"allele_{r['called_allele2']}"
        return r["called_allele1"], r["called_allele2"], unit.get(x, "NA"), unit.get(y, "NA"), r["excess"]

    def verdict(c):
        if c is None or "LOW_INFORMATION" in (t1, t2) or "NA" in (c[2], c[3]):
            return "NOT_ASSESSED"
        if sorted([c[2], c[3]]) == sorted([t1, t2]):
            return "UNIT_PAIR_CORRECT"
        return "ALLOCATION_ERROR" if c[2] + c[3] == t1 + t2 else "DOSAGE_ERROR"

    with open(f"{a.out}.structure.tsv", "w") as f:
        f.write("key\tvalue\n")
        f.write(f"case\t{a.case}\nblock\t{a.block}\nunit_source\t{source} ({len(s0)} bp; carried by {carried[source]})\n"
                f"anchors_unique_in_source\t{len(anchors)}\nanchors_conserved\t{len(conserved)}\n"
                f"anchors_informative\t{len(informative)}\ninformative_span_bp\t{span}\n")
        for n in sorted(truth):
            f.write(f"truth_{n}\t{len(truth[n])} bp, {unit[n]} unit copies\n")
        f.write(f"truth_total_units\t{t1 + t2 if 'LOW_INFORMATION' not in (t1, t2) else 'NA'}\n")
        for label, path in (("historical", a.hist_eval), ("fresh", a.fresh_eval)):
            c = call(path)
            if c is None:
                f.write(f"{label}_call\tNA (no evaluator row)\n")
                continue
            f.write(f"{label}_called_alleles\t{c[0]},{c[1]}\n{label}_called_units\t{c[2]}+{c[3]}\n"
                    f"{label}_target_excess\t{c[4]}\n{label}_structure_verdict\t{verdict(c)}\n")
    print(open(f"{a.out}.structure.tsv").read(), end="")


if __name__ == "__main__":
    main()
