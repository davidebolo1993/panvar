#!/usr/bin/env python3
"""M3-F POST-HOC diagnostic (declared after the frozen attribution; truth-guided; changes no class, prediction or setting).

The frozen attribution compares PHASED paths (exact spelled pair). For an unphased (genotype-equivalent) call this mixes
phase with genotype. Here, per mosaic error cell with a diagnostic catalogue: the incumbent's unordered genotype errors
against the truth, whether any catalogue member carries the true unordered genotype at every block (and its rank and
F/E/P deltas), and whether each wrong block needs one or both homologues changed.

  posthoc_genotype_level.py <runs dir> cell...
"""
import csv
import sys

runs, cells = sys.argv[1], sys.argv[2:]
for cell in cells:
    tb = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(f"{runs}/tbeval.{cell}/eval.eval_summary.tsv")}
    t = [list(map(int, tb["truth_hap1_path"].split(","))), list(map(int, tb["truth_hap2_path"].split(",")))]
    tg = [tuple(sorted(p)) for p in zip(*t)]
    rows = sorted(csv.DictReader(open(f"{runs}/lik.{cell}/lik.scores.tsv"), delimiter="\t"), key=lambda r: -float(r["total"]))
    g = lambda r: [tuple(sorted(p)) for p in zip(map(int, r["hap1_alleles"].split(",")), map(int, r["hap2_alleles"].split(",")))]
    inc = next(r for r in rows if r["label"] == "incumbent")
    ig = g(inc)
    wrong = [b for b in range(len(tg)) if ig[b] != tg[b]]
    print(f"{cell}: catalogue {len(rows)} (incl. incumbent); incumbent rank {rows.index(inc) + 1}")
    for b in wrong:
        shared = len(set(ig[b]) & set(tg[b])) if ig[b][0] != ig[b][1] or tg[b][0] != tg[b][1] else (2 if ig[b] == tg[b] else 0)
        print(f"  block index {b}: called {ig[b]} truth {tg[b]}: {'one homologue' if shared >= 1 else 'BOTH homologues'} must change")
    ok = [r for r in rows if g(r) == tg]
    print(f"  catalogue members with the true unordered genotype at every block: {len(ok)}")
    for r in ok[:3]:
        print(f"    {r['label']} rank {rows.index(r) + 1}: total-inc {float(r['total']) - float(inc['total']):.4f} "
              f"(F {float(r['fragment_term']) - float(inc['fragment_term']):.4f}, E {float(r['exposure_term']) - float(inc['exposure_term']):.4f}, "
              f"P {float(r['prior']) - float(inc['prior']):.4f})")
