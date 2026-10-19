#!/usr/bin/env python3
"""M5-A summary (after the freeze): joins the frozen combination with the truth evaluation (score-set distances/excess per
candidate and block), the score components (F, E, prior) of every original and terminal, and the ledger costs.
Recommendation rule (frozen):
  BENEFIT  a cell whose combined decision changes the unordered genotype and lowers excess, no control changed or harmed
  HARM     a combined decision that changes the genotype and raises excess, or any control whose genotype changes
  -> INTEGRATE if >= 1 BENEFIT and 0 HARM; INVESTIGATE if any HARM or a score-improving genotype change with equal excess;
     PARK if no cell's combined genotype changes (no scientific benefit to weigh against the cost)

  summary.py <runs dir> <ledger> <out prefix> cell...
"""
import csv
import sys

runs, ledger, prefix, cells = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
T = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
cellrows = {r["cell"]: r for r in T(f"{runs}/combine/m5a.cells.tsv")}
led = {r["unit"]: r for r in T(ledger)}
lines = ["cell\tlabel\ttotal\tfragment\texposure\tprior\texpected_reproduced\tdistance\tfloor\texcess\tin_combined_top_set"]
verdicts = []
for c in cells:
    ev = {r["label"]: r for r in T(f"{runs}/eval.{c}/ev.scores.tsv")}
    top = set(cellrows[c]["top_set"].split(","))
    for k, r in ev.items():
        rep = "NA" if r["expected_total"] == "NA" else str(abs(float(r["total"]) - float(r["expected_total"])) <= 1e-6)
        lines.append("\t".join([c, k, r["total"], r["fragment_term"], r["exposure_term"], r["prior"], rep, r["distance"], r["floor"], r["excess"],
                                str(k in top)]))
    orig = min(int(ev[k]["excess"]) for k in ev if k.startswith("orig_"))
    comb = sorted({int(ev[k]["excess"]) for k in top if k in ev})
    changed = cellrows[c]["combined_genotype_equals_original"] != "True"
    control = c != "c4_three_switches_1_2"
    if not changed:
        v = "NO_GENOTYPE_CHANGE"
    elif control:
        v = "HARM(control genotype changed)"
    elif max(comb) < orig:
        v = "BENEFIT"
    elif min(comb) > orig:
        v = "HARM"
    else:
        v = "EQUAL_OR_MIXED"
    verdicts.append((c, v, orig, comb, cellrows[c]["combined_state"], cellrows[c]["best_minus_original"]))
lines += ["", "cell\tverdict\toriginal_excess\tcombined_top_excess\tcombined_state\tbest_minus_original_total"]
lines += ["\t".join(map(str, v)) for v in verdicts]
vs = [v[1] for v in verdicts]
rec = ("INTEGRATE" if vs.count("BENEFIT") >= 1 and not any(x.startswith("HARM") for x in vs) else
       "INVESTIGATE" if any(x.startswith("HARM") or x == "EQUAL_OR_MIXED" for x in vs) else "PARK")
lines += ["", "unit\twall_seconds\tpeak_rss_gib\texit"] + [f"{u}\t{r['wall_seconds']}\t{r['peak_rss_gib']}\t{r['exit']}" for u, r in led.items() if u.startswith("leg.")]
lines += ["", f"RECOMMENDATION\t{rec}"]
open(prefix + ".tsv", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
