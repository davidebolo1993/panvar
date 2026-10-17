#!/usr/bin/env python3
"""S1-A: reporting-only projection of one CURRENT caller run onto the proposed sequence-first outputs
(tests/genotype_sequence_label_reporting_design.md section 7). Reads the run's own tables and the model dump; writes
nothing into the run. No score, prior, posterior or calibration is computed here: every number shown comes from the
caller's tables, and every proposed field is a deterministic function of them.

Projection coverage: block projections are COMPLETE over the model's allele tuples for every explored sequence group
(an exhaustive spelling search, fixture-scale); sequence groups are the EXPLORED pool only (mosaic_scores.tsv).
"""
import csv
import hashlib
import os

HARD_FAILURES = {"CANDIDATE_INCOMPLETE", "INPUT_INVALID", "MODEL_REFUSED", "FAILED", "RESOURCE_REFUSAL", "TIME_LIMIT"}


def rows(path):
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"))


def kv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path).read().splitlines():
            p = line.split("\t")
            if len(p) >= 2:
                out.setdefault(p[0], p[1])
    return out


def h16(s):
    return hashlib.sha256(s.encode()).hexdigest()[:16]


class Model:
    def __init__(self, dump_prefix):
        self.blocks, self.ids = [], []
        for r in rows(dump_prefix + ".block_alleles.tsv"):
            b = int(r["block_index"])
            while len(self.blocks) <= b:
                self.blocks.append([])
                self.ids.append(None)
            self.blocks[b].append(r["sequence"])
            self.ids[b] = r["block_id"]

    def spell(self, labels):
        return "".join(self.blocks[b][int(a)] for b, a in enumerate(labels.split(",")))

    def projections(self, sequence):
        """Every allele tuple (one per block, in order) whose concatenation is `sequence` -- complete, by DFS."""
        out, n = [], len(self.blocks)

        def walk(b, at, path):
            if b == n:
                if at == len(sequence):
                    out.append(tuple(path))
                return
            for i, a in enumerate(self.blocks[b]):
                if sequence.startswith(a, at):
                    walk(b + 1, at + len(a), path + [i])
        walk(0, 0, [])
        return out


def project(prefix, dump_prefix):
    model = Model(dump_prefix)
    st = kv(prefix + ".mosaic_status.tsv")
    status = st.get("status", "NA")
    current = {"status": status, "call_class": st.get("call_class", "NA"), "no_call_reason": st.get("no_call_reason", "NA"),
               "fasta_records": 0, "block_rows": 0, "identified_block_rows": 0}
    if os.path.exists(prefix + ".mosaic_haplotypes.fa"):
        current["fasta_records"] = sum(1 for l in open(prefix + ".mosaic_haplotypes.fa") if l.startswith(">"))
    current_gq = {}
    for r in rows(prefix + ".block_genotypes.tsv"):
        current["block_rows"] += 1
        current["identified_block_rows"] += 1 if r.get("allele1", "NA") != "NA" else 0
        current_gq[r["block_id"]] = r.get("conditional_GQ", "NA")
    # explored labellings per sequence group, and the groups' scores, probabilities and decision-set membership
    labellings = {}
    for r in rows(prefix + ".mosaic_candidates.tsv"):
        labellings.setdefault(r["sequence_candidate_id"], []).append((r["hap1_alleles"], r["hap2_alleles"]))
    for r in rows(prefix + ".equivalence_set.tsv"):
        labellings.setdefault(r["sequence_candidate_id"], []).append((r["hap1_alleles"], r["hap2_alleles"]))
    groups = {}
    for r in rows(prefix + ".mosaic_scores.tsv"):
        g = r["sequence_candidate_id"]
        if g not in labellings:
            continue
        a, b = labellings[g][0]
        s1, s2 = model.spell(a), model.spell(b)
        pair = tuple(sorted((s1, s2)))
        geno = [set() for _ in model.blocks]   # every block genotype over every projection of the pair (complete)
        for p1 in model.projections(pair[0]):
            for p2 in model.projections(pair[1]):
                for blk in range(len(model.blocks)):
                    geno[blk].add(tuple(sorted((p1[blk], p2[blk]))))
        explored = [set() for _ in model.blocks]   # what the caller explored (labellings it holds)
        for x, y in labellings[g]:
            for blk, (i, j) in enumerate(zip(x.split(","), y.split(","))):
                explored[blk].add(tuple(sorted((int(i), int(j)))))
        groups[g] = {"score": float(r["score"]), "prob": float(r["candidate_set_probability"]), "top": r["in_winner_set"] == "true",
                     "pair": pair, "hashes": tuple(sorted((h16(pair[0]), h16(pair[1])))), "geno": geno, "explored": explored,
                     "ambiguous_flag": r["projection_ambiguous"] == "true"}
    proposed = {"prediction": "NO_PREDICTION", "sequence_status": "NA", "provisional": "no", "representative": "none",
                "fasta": "none", "blocks": [], "intervals": [], "alternatives": [], "pool_conditional_posterior": "NA",
                "structure_confidence": "NA (UNCALIBRATED)", "coverage": "projections complete over model allele tuples; groups = explored pool"}
    provisional = status == "REPAIR_INCOMPLETE" and st.get("repair_final_incumbent_hap1", "NA") not in ("NA", "")
    if status in HARD_FAILURES or (not groups and not provisional):
        proposed["reason"] = "hard failure or no admissible evaluated candidate"
        return current, proposed
    if not groups and provisional:   # the search stopped: its last committed incumbent is a valid, exactly scored candidate
        s1, s2 = model.spell(st["repair_final_incumbent_hap1"]), model.spell(st["repair_final_incumbent_hap2"])
        proposed.update(prediction="PREDICTION", sequence_status="NOT_DETERMINED", provisional="yes",
                        representative="/".join(sorted((h16(s1), h16(s2)))), fasta="/".join(sorted((h16(s1), h16(s2)))),
                        reason="search incomplete (" + st.get("repair_stop_reason", "NA") + "); the last committed incumbent, "
                               "exactly scored; decision set not determined; confidence not defined")
        return current, proposed
    top = {g: v for g, v in groups.items() if v["top"]}
    # canonical representative: best score, then the smallest canonical sequence-hash pair (never candidate order)
    rep_id = sorted(top, key=lambda g: (-top[g]["score"], top[g]["hashes"]))[0]
    rep = top[rep_id]
    n_blocks = len(model.blocks)
    union = [set().union(*(top[g]["geno"][b] for g in top)) for b in range(n_blocks)]
    if len(top) == 1:
        sstatus = "RESOLVED"
    elif all(len(u) == 1 for u in union):
        sstatus = "GENOTYPE_RESOLVED_SEQUENCE_SET"
    else:
        sstatus = "SEQUENCE_TIED"
    blocks, run = [], []
    for b in range(n_blocks):
        if len(union[b]) == 1:
            g = next(iter(union[b]))
            # projection-aware error-mass bounds from the explored pool: lower = groups whose every projection differs,
            # upper = groups with any projection that differs; one weight per sequence group; no point estimate invented
            lo = sum(v["prob"] for v in groups.values() if all(x != g for x in v["geno"][b]))
            hi = sum(v["prob"] for v in groups.values() if any(x != g for x in v["geno"][b]))
            blocks.append((model.ids[b], "IDENTIFIED", f"{g[0]}/{g[1]}", f"{lo:.3g}", f"{hi:.3g}", current_gq.get(model.ids[b], "NA")))
        elif len(top) == 1:
            blocks.append((model.ids[b], "REPRESENTATION_AMBIGUOUS", ";".join(f"{a}/{c}" for a, c in sorted(union[b])), "NA", "NA",
                           current_gq.get(model.ids[b], "NA")))
        else:
            blocks.append((model.ids[b], "GENOTYPE_UNRESOLVED", ";".join(f"{a}/{c}" for a, c in sorted(union[b])), "NA", "NA",
                           current_gq.get(model.ids[b], "NA")))
    for b, row in enumerate(blocks):   # representation intervals: maximal runs of representation-ambiguous blocks
        if row[1] == "REPRESENTATION_AMBIGUOUS":
            run.append(b)
        if run and (row[1] != "REPRESENTATION_AMBIGUOUS" or b == n_blocks - 1):
            proposed["intervals"].append(f"{model.ids[run[0]]}..{model.ids[run[-1]]}")
            run = []
    alternatives = sorted(groups, key=lambda g: (-groups[g]["score"], groups[g]["hashes"]))[:5]
    for g in alternatives:
        v = groups[g]
        differing = sum(1 for b in range(n_blocks) if v["geno"][b] != rep["geno"][b])
        proposed["alternatives"].append(f"{'/'.join(v['hashes'])}:p={v['prob']:.3g}:{'top' if v['top'] else 'below'}:"
                                        f"{'same_sequence' if v['pair'] == rep['pair'] else 'distinct'}:blocks_differing={differing}")
    proposed.update(prediction="PREDICTION", sequence_status=sstatus, representative="/".join(rep["hashes"]),
                    fasta="/".join(rep["hashes"]), blocks=blocks, pool_conditional_posterior=f"{rep['prob']:.6g}",
                    provisional="yes" if provisional else "no",
                    reason="completed ranking" if not provisional else "search incomplete; provisional")
    lower_ambiguous = [g for g, v in groups.items() if not v["top"] and v["prob"] > 0 and any(len(x) > 1 for x in v["geno"])]
    proposed["lower_projection_ambiguous_groups"] = len(lower_ambiguous)
    unexplored = sum(1 for v in groups.values() for b in range(n_blocks) if v["geno"][b] - v["explored"][b])
    proposed["projections_unexplored_by_caller"] = unexplored
    return current, proposed
