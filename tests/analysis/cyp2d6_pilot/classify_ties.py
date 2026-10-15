#!/usr/bin/env python3
"""CYP2D6/D7 pilot (AMENDMENT_5): classify each cell's outcome from call outputs and the model alone (no truth).

The caller's status is kept unchanged. A no-call is additionally classified
SEQUENCE_IDENTICAL_BLOCK_LABEL_AMBIGUITY only when every member of its tied set (every row of the caller's
equivalence_set.tsv for a tie no-call; a divergent_tie ranks its tied members 1..n although their totals are equal):
  1. shares one sequence_candidate_id;
  2. has identical score (total), fragment_term, exposure_term and sequence_pair_prior (exact text);
  3. spells the same unordered diploid sequence (each haplotype spelled from the call model's block alleles);
  4. differs from the first member only at blocks lying inside a model-derived degeneracy interval
     (degeneracy/degeneracy_summary.tsv), i.e. only in the allocation of allele labels among blocks;
  5. and no distinct sequence is present in the top set (implied by 1 and 3, checked separately).
Every tie no-call (caller reasons distinct_tie and divergent_tie) is examined; a tie failing any criterion is
DISTINCT_SEQUENCE_TIE. Other no-calls keep their reason; calls are CALL_<call_class>.
The classification picks no canonical allele combination and invents no block genotype.

Usage: classify_ties.py <pilot package> <runs dir name> <out tsv> <job>...
"""
import csv
import os
import sys


def kv(path):
    out = {}
    for line in open(path).read().splitlines():
        p = line.split("\t")
        if len(p) >= 2:
            out[p[0]] = p[1]
    return out


def model_alleles(pkg, model):
    blocks, ids = [], []
    with open(os.path.join(pkg, "degeneracy", "models", f"{model}.block_alleles.tsv")) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            b = int(r["block_index"])
            while len(blocks) <= b:
                blocks.append([])
                ids.append(r["block_id"])
            blocks[b].append(r["sequence"])
    return blocks, ids


def intervals(pkg, model, ids):
    with open(os.path.join(pkg, "degeneracy", "degeneracy_summary.tsv")) as f:
        row = next(r for r in csv.DictReader(f, delimiter="\t") if r["model"] == model)
    out = []
    if row["coalescence_intervals"] != "none":
        for iv in row["coalescence_intervals"].split(";"):
            lo, hi = iv.split("..")
            out.append((ids.index(lo), ids.index(hi)))
    return out


def classify(pkg, runs, job, jobrow):
    prefix = os.path.join(pkg, "data", runs, job)
    rec = dict(zip(*[l.split("\t") for l in open(prefix + ".run.tsv").read().splitlines()[:2]]))
    st = kv(prefix + ".mosaic_status.tsv")
    model = "original" if jobrow["arm"] in ("lzo", "mosaic") else f"loo.{jobrow['sample']}"
    exit_code, exit_source = rec["exit"], "runner"
    if exit_code not in ("0", "1") and st.get("status") not in (None, "RUNNING") and "seconds_total" in st:
        # the runner lost the process (STOP_4: the backstop killed only the time wrapper) but the caller finished on
        # its own and wrote a final status; its outcome is read from that status
        exit_code = "0" if st.get("call_class") in ("strict_call", "genotype_equivalent_call") else "1"
        exit_source = f"caller_status (runner record: {rec['exit']})"
    out = {"job": job, "arm": jobrow["arm"], "sample": jobrow["sample"], "coverage": jobrow["coverage"],
           "exit": exit_code, "exit_source": exit_source, "caller_status": st.get("status", "NA"), "call_class": st.get("call_class", "NA"),
           "no_call_reason": st.get("no_call_reason", "NA"), "model": model, "top_set_members": "NA",
           "sequence_candidate_ids": "NA", "scores_identical": "NA", "unordered_sequence_identical": "NA",
           "differing_blocks": "NA", "differences_within_degeneracy_interval": "NA", "classification": "NA"}
    if exit_code == "0":
        out["classification"] = f"CALL_{st.get('call_class', 'NA')}"
        return out
    if exit_code != "1":
        out["classification"] = f"NOT_COMPLETED_exit_{exit_code}"
        return out
    if st.get("no_call_reason") not in ("distinct_tie", "divergent_tie"):
        reason = st.get("no_call_reason", "NA")
        out["classification"] = f"NO_CALL_{reason if reason != 'NA' else st.get('status', 'NA')}"
        return out
    with open(prefix + ".equivalence_set.tsv") as f:
        top = list(csv.DictReader(f, delimiter="\t"))
    blocks, ids = model_alleles(pkg, model)
    ivs = intervals(pkg, model, ids)
    spell = lambda alleles: "".join(blocks[b][int(a)] for b, a in enumerate(alleles.split(",")))  # noqa: E731
    seq_ids = sorted({r["sequence_candidate_id"] for r in top})
    scores_same = all(len({r[k] for r in top}) == 1
                      for k in ("score", "fragment_term", "exposure_term", "sequence_pair_prior"))
    dip = [tuple(sorted((spell(r["hap1_alleles"]), spell(r["hap2_alleles"])))) for r in top]
    seq_same = len(set(dip)) == 1
    first = [top[0]["hap1_alleles"].split(","), top[0]["hap2_alleles"].split(",")]
    differing, inside = set(), True
    for r in top[1:]:
        other = [r["hap1_alleles"].split(","), r["hap2_alleles"].split(",")]
        # the haplotype matching under which the sequences agree (either assignment)
        best = None
        for perm in ((0, 1), (1, 0)):
            if all(spell(",".join(first[h])) == spell(",".join(other[perm[h]])) for h in (0, 1)):
                best = perm
                break
        if best is None:
            inside = False
            continue
        for h in (0, 1):
            for b, (x, y) in enumerate(zip(first[h], other[best[h]])):
                if x != y:
                    differing.add(b)
                    if not any(lo <= b <= hi for lo, hi in ivs):
                        inside = False
    out.update({"top_set_members": len(top), "sequence_candidate_ids": ";".join(seq_ids),
                "scores_identical": str(scores_same).lower(), "unordered_sequence_identical": str(seq_same).lower(),
                "differing_blocks": ",".join(ids[b] for b in sorted(differing)) or "none",
                "differences_within_degeneracy_interval": str(inside and bool(differing)).lower()})
    ok = len(top) >= 2 and len(seq_ids) == 1 and scores_same and seq_same and inside and bool(differing)
    out["classification"] = "SEQUENCE_IDENTICAL_BLOCK_LABEL_AMBIGUITY" if ok else "DISTINCT_SEQUENCE_TIE"
    return out


def main(pkg, runs, out_path, jobs):
    with open(os.path.join(pkg, "jobs.tsv")) as f:
        table = {r["job"]: r for r in csv.DictReader(f, delimiter="\t")}
    rows = [classify(pkg, runs, j, table[j]) for j in jobs]
    with open(out_path, "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(f"{r['job']}\t{r['classification']}\t{r['differing_blocks']}")


if __name__ == "__main__":
    main(os.path.abspath(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4:])
