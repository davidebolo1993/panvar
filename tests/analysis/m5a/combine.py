#!/usr/bin/env python3
"""M5-A combination of the pooled-winner continuation legs (TRUTH-FREE; frozen before execution).

Per root: the leg's run state, rounds, moves, stop reason, whether round 1 moved (the original call had an improving
neighbour), ties (stop reason rescore_incomplete or tie_at_best: the engine never breaks a tie), and the terminal incumbent
(the last round's incumbent; a leg that ends on a moved round hit its round bound and is not converged).
Per cell, the diagnostic candidates are every ORIGINAL root (the reported prediction is retained) plus every completed
leg's terminal. The combined top set holds every candidate within 1e-6 of the best total; if it spans more than one
unordered block genotype the combined decision is AMBIGUOUS (no representative is chosen); several sequence pairs with
one unordered genotype are GENOTYPE_EQUIVALENT. An incomplete leg (refusal, time limit, failure) makes the cell INCOMPLETE:
its ambiguity is never discarded. Writes per-cell candidate files for the truth evaluation and a selection digest.

  combine.py <runs dir> <ROOTS.tsv> <out prefix>
"""
import csv
import hashlib
import os
import sys

runs, roots_path, prefix = sys.argv[1:4]
roots = list(csv.DictReader(open(roots_path), delimiter="\t"))
KV = lambda p: {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(p) if "\t" in l} if os.path.exists(p) else {}
geno = lambda h1, h2: tuple(tuple(sorted(p)) for p in zip(map(int, h1.split(",")), map(int, h2.split(","))))
pair = lambda h1, h2: tuple(sorted((h1, h2)))
root_rows = ["cell\troot\tleg_state\trounds\tmoves\tstop_reason\toriginal_had_improving_neighbour\ttie\tconverged\tstart_total\tterminal_total\t"
             "delta_total\tdelta_fragment\tdelta_exposure\tdelta_prior\tterminal_hap1\tterminal_hap2\tblocks_changed_unordered\tphase_only_change\t"
             "footprint_peak_bytes"]
cell_rows = ["cell\tcandidates\tcombined_state\ttop_set\ttop_genotypes\ttop_sequence_pairs\tcombined_best_total\toriginal_total\tbest_minus_original\t"
             "combined_genotype_equals_original"]
by_cell = {}
for r in roots:
    by_cell.setdefault(r["cell"], []).append(r)
for cell, rs in by_cell.items():
    cands, incomplete = [], False
    for r in rs:
        cands.append((f"orig_{r['root']}", r["hap1"], r["hap2"], float(r["score"])))
        d = f"{runs}/leg.{cell}.{r['root']}"
        st, sm = KV(f"{d}/leg.run_status.tsv"), KV(f"{d}/leg.summary.tsv")
        state = st.get("state", st.get("status", "NA"))
        rounds = list(csv.DictReader(open(f"{d}/leg.rounds.tsv"), delimiter="\t")) if os.path.exists(f"{d}/leg.rounds.tsv") else []
        if not rounds or not sm:
            incomplete = True
            root_rows.append("\t".join([cell, r["root"], state] + ["NA"] * 17))
            continue
        last = rounds[-1]
        converged = last["moved"] == "false"
        th1, th2 = last["hap1_alleles"], last["hap2_alleles"]
        tt = float(last["start_total"]) if converged else float(last["total"])
        moves = sum(1 for x in rounds if x["moved"] == "true")
        first = rounds[0]
        # start components from round 1 (the incumbent's row is the start), terminal components from the last row
        s_frag = s_exp = s_prior = None
        g0, g1 = geno(r["hap1"], r["hap2"]), geno(th1, th2)
        changed = [b for b in range(len(g0)) if g0[b] != g1[b]]
        phase_only = (not changed) and pair(th1, th2) != pair(r["hap1"], r["hap2"])
        tie = sm.get("stop_reason", "") in ("rescore_incomplete", "tie_at_best")
        root_rows.append("\t".join(map(str, [cell, r["root"], state, len(rounds), moves, sm.get("stop_reason", "NA"), first["moved"] == "true", tie,
                                             converged, r["score"], f"{tt:.10f}", f"{tt - float(r['score']):.6f}", "see eval", "see eval", "see eval",
                                             th1, th2, ",".join(map(str, changed)) or "-", phase_only, sm.get("phys_footprint_peak_bytes", "NA")])))
        if not converged:
            incomplete = True   # a leg at its round bound is reported, never treated as a finished optimum
        cands.append((f"term_{r['root']}", th1, th2, tt))
    best = max(t for *_, t in cands)
    top = [x for x in cands if x[3] >= best - 1e-6]
    genos = {geno(x[1], x[2]) for x in top}
    seqs = {pair(x[1], x[2]) for x in top}
    orig = max(float(r["score"]) for r in rs)
    orig_geno = {geno(r["hap1"], r["hap2"]) for r in rs}
    state = "INCOMPLETE" if incomplete else ("AMBIGUOUS" if len(genos) > 1 else ("GENOTYPE_EQUIVALENT" if len(seqs) > 1 else "UNIQUE"))
    cell_rows.append("\t".join(map(str, [cell, len(cands), state, ",".join(x[0] for x in top), len(genos), len(seqs), f"{best:.10f}", f"{orig:.10f}",
                                         f"{best - orig:.6f}", genos == orig_geno])))
    with open(f"{prefix}.candidates.{cell}.tsv", "w") as f:
        seen = set()
        for label, h1, h2, t in cands:
            key = pair(h1, h2)
            if key in seen:
                continue
            seen.add(key)
            f.write(f"{label}\t{h1}\t{h2}\t{t!r}\n")
open(prefix + ".roots.tsv", "w").write("\n".join(root_rows) + "\n")
open(prefix + ".cells.tsv", "w").write("\n".join(cell_rows) + "\n")
h = hashlib.sha256(open(prefix + ".roots.tsv", "rb").read() + open(prefix + ".cells.tsv", "rb").read()).hexdigest()
open(prefix + ".selection.sha256", "w").write(h + "\n")
print("\n".join(root_rows))
print("\n".join(cell_rows))
print("selection digest", h)
