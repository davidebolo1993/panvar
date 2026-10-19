#!/usr/bin/env python3
"""M1-B Amendment 1: declare <= 12 alternatives for one persistent case by a frozen rule, and report whether the fresh caller
proposed and evaluated a structurally appropriate alternative (saved outputs only; truth enters through the diagnostic
structure table and the evaluator's closest-allele columns, never a caller).

Rule (in order, capped at 12):
  1. `reported`: the rank-1 row of the fresh run's equivalence_set.tsv (full-precision saved total).
  2. `saved_<a>+<b>`: for each unit-copy class at the target block (unordered pair of per-haplotype unit copies from
     structure.py), the best-scoring saved scored sequence candidate (mosaic_scores.tsv, its representative beam
     candidate's allele paths), classes other than the reported candidate's own, best score first, at most 7 classes;
     expected total = the saved score (printed to 1e-6).
  3. `sub_truth_a`, `sub_truth_b`: the reported diplotype with ONLY the target block replaced by the evaluator's closest
     panel alleles to truth hap1 and hap2, in both homologue assignments (no expectation; never a whole-locus oracle).
  4. `cross_<case>`: the reported diplotype of a declared twin run (CYP2D6 5x only: the fresh 30x report), no expectation.
Decomposition pair: reported vs the best saved true-class candidate if declared, else sub_truth_a.
Membership: number of saved scored sequence candidates per class, the best rank and score gap of the true class.

  alternatives.py --run PREFIX --units UNITS.tsv --eval-blocks EVAL.tsv --block ID --out OUT_PREFIX [--twin PREFIX --twin-label L]
"""
import argparse
import csv


def read(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    for a in ("run", "units", "eval_blocks", "block", "out"):
        ap.add_argument("--" + a.replace("_", "-"), required=True)
    ap.add_argument("--twin")
    ap.add_argument("--twin-label")
    a = ap.parse_args()
    ev = read(a.eval_blocks)
    t = [r["block_id"] for r in ev].index(a.block)
    target = ev[t]
    units = {r["name"]: r["unit_copies"] for r in read(a.units)}
    true_class = tuple(sorted(int(units[n]) for n in ("truth1", "truth2")))

    def cls(h1, h2):
        return tuple(sorted(int(units[f"allele_{p.split(',')[t]}"]) for p in (h1, h2)))

    def first_alt(prefix):
        r = next(r for r in read(prefix + ".equivalence_set.tsv") if r["rank"] == "1")
        return r["hap1_alleles"], r["hap2_alleles"], r["score"], r["sequence_candidate_id"]

    cand = {r["candidate_id"]: r for r in read(a.run + ".mosaic_candidates.tsv")}
    scores = sorted(read(a.run + ".mosaic_scores.tsv"), key=lambda r: -float(r["score"]))
    rh1, rh2, rscore, rid = first_alt(a.run)
    rcls = cls(rh1, rh2)
    best, count, ranks = {}, {}, {}
    for rank, s in enumerate(scores, 1):
        c = cand[s["representative_beam_candidate"]]
        k = cls(c["hap1_alleles"], c["hap2_alleles"])
        count[k] = count.get(k, 0) + 1
        if k not in best:
            best[k] = (s, c)
            ranks[k] = rank
    out = [("reported", rh1, rh2, rscore)]
    for k in sorted((k for k in best if k != rcls), key=lambda k: -float(best[k][0]["score"]))[:7]:
        s, c = best[k]
        out.append((f"saved_{k[0]}+{k[1]}", c["hap1_alleles"], c["hap2_alleles"], s["score"]))

    def sub(path, allele):
        p = path.split(",")
        p[t] = allele
        return ",".join(p)

    x, y = target["truth1_closest_allele"], target["truth2_closest_allele"]
    out.append(("sub_truth_a", sub(rh1, x), sub(rh2, y), "NA"))
    out.append(("sub_truth_b", sub(rh1, y), sub(rh2, x), "NA"))
    if a.twin:
        th1, th2, _, _ = first_alt(a.twin)
        out.append((f"cross_{a.twin_label}", th1, th2, "NA"))
    out = out[:12]
    with open(a.out + ".candidates.tsv", "w") as f:
        f.write(f"# M1-B Amendment 1 rule; target {a.block} (block index {t}); true unit class {true_class[0]}+{true_class[1]}\n")
        for label, h1, h2, e in out:
            f.write(f"{label}\t{h1}\t{h2}\t{e}\n")
    with open(a.out + ".membership.tsv", "w") as f:
        f.write("unit_class\tsaved_scored_sequence_candidates\tbest_rank\tbest_score\tgap_to_reported\tis_true_class\tis_reported_class\n")
        for k in sorted(count, key=lambda k: ranks[k]):
            f.write(f"{k[0]}+{k[1]}\t{count[k]}\t{ranks[k]}\t{best[k][0]['score']}\t"
                    f"{float(best[k][0]['score']) - float(rscore):.6f}\t{k == true_class}\t{k == rcls}\n")
    labels = [o[0] for o in out]
    b = f"saved_{true_class[0]}+{true_class[1]}"
    open(a.out + ".decompose.txt", "w").write(f"reported,{b if b in labels else 'sub_truth_a'}\n")
    print(open(a.out + ".membership.tsv").read(), end="")
    print(f"declared {len(out)} alternatives; truth closest alleles {x} ({units.get('allele_' + x)} units), "
          f"{y} ({units.get('allele_' + y)} units)")


if __name__ == "__main__":
    main()
