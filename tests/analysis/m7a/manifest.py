#!/usr/bin/env python3
"""M7-A frozen candidate manifest (read-only over saved records; no scoring, no candidate generation from truth).

Tiers:
  D  deployable: (1) the run's reported call (equivalence-set rank 1) and other equivalence-set members;
     (2) the caller's own best saved scored candidate per target-block unit class (M1-B Amendment 1 rule, steps 1-2,
     which are truth-free; recomputed here for the LZO controls, read from M1-B for the LOO cells);
     (3) M6-A's frozen count-model evaluation selection (truth-free: caller + best count group + first members of
     groups 2-5).
  X  cross-run: a prediction from another run of the same donor (CYP2D6 LOO: M1-B's `cross_30x`); never pooled with D.
  O  oracle diagnostics (truth-selected or truth-patched; never deployable): M1-B `sub_truth_a/b`, M1-D constructions,
     M1-E constructions (incl. the complete true diplotype as a reference).
A candidate reached by both a truth-free and a truth-selected rule stays in D, with every provenance recorded.

  manifest.py <out dir>
"""
import csv
import hashlib
import os
import sys
from collections import OrderedDict

P = "/Users/davide.bolognini/panvar/tests/results"
out = sys.argv[1]
T = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
sha = lambda p: hashlib.sha256(open(p, "rb").read()).hexdigest()
CELLS = OrderedDict([
    ("c4_loo_na19650_10", dict(locus="C4", block="sv:7", t=7,
        model=f"{P}/genotype_c4_production_11fe681/graphs/NA19650.loo.gfa", reads=f"{P}/genotype_c4_production_11fe681/reads/NA19650/cov10",
        run=f"{P}/m1b_structural_20261019/runs/call.c4_10/call", units=f"{P}/m1b_structural_20261019/runs/structure.c4_10/c4_10.units.tsv",
        alts=f"{P}/m1b_structural_20261019/runs/alts.c4_10/c4_10.candidates.tsv", m6a="c4_loo_na19650_10",
        m1d=f"{P}/m1d_unit_reachability_20261019/runs/construct.c4_10.a2", m1e=f"{P}/m1e_residual_sensitivity_20261019/runs/construct.c4_10")),
    ("cyp_loo_na19650_5", dict(locus="CYP2D6", block="sv:9", t=9,
        model=f"{P}/cyp2d6_shortread_pilot_v1/data/graphs/loo.NA19650.gfa", reads=f"{P}/cyp2d6_shortread_pilot_v1/data/reads/NA19650/cov5",
        run=f"{P}/m1b_structural_20261019/runs/call.cyp_5/call", units=f"{P}/m1b_structural_20261019/runs/structure.cyp_5/cyp_5.units.tsv",
        alts=f"{P}/m1b_structural_20261019/runs/alts.cyp_5/cyp_5.candidates.tsv", m6a="cyp_loo_na19650_5", m1d=None, m1e=None)),
    ("c4_lzo_hg00171_5", dict(locus="C4", block="sv:7", t=7,
        model=f"{P}/genotype_c4_production_11fe681/graphs/c4.gfa", reads=f"{P}/genotype_c4_production_11fe681/reads/HG00171/cov5",
        run=f"{P}/m4b_distance_prior_20261019/runs/run.c4_lzo_hg00171_5.const/call", units=f"{P}/m6a_count_model_20261019/runs/units.c4_lzo_hg00171_5/units.tsv",
        alts=None, m6a="c4_lzo_hg00171_5", m1d=None, m1e=None)),
    ("cyp_lzo_hg02769_5", dict(locus="CYP2D6", block="sv:9", t=9,
        model=f"{P}/cyp2d6_shortread_pilot_v1/data/graphs/cyp2d6.gfa", reads=f"{P}/cyp2d6_shortread_pilot_v1/data/reads/HG02769/cov5",
        run=f"{P}/m4b_distance_prior_20261019/runs/run.cyp_lzo_hg02769_5.const/call", units=f"{P}/m6a_count_model_20261019/runs/units.cyp_lzo_hg02769_5/units.tsv",
        alts=None, m6a="cyp_lzo_hg02769_5", m1d=None, m1e=None)),
])
CURATION = "/Users/davide.bolognini/panvar/tests/genotype_cyp2d6_curation.tsv"
bindings, rows = [], []


def bind(cell, role, path):
    bindings.append((cell, role, path, sha(path)))


for cell, c in CELLS.items():
    for role, p in [("model", c["model"]), ("reads_1", c["reads"] + "/reads_1.fq"), ("reads_2", c["reads"] + "/reads_2.fq"), ("library", c["reads"] + "/library.txt"),
                    ("equivalence_set", c["run"] + ".equivalence_set.tsv"), ("mosaic_scores", c["run"] + ".mosaic_scores.tsv"),
                    ("mosaic_candidates", c["run"] + ".mosaic_candidates.tsv"), ("mosaic_status", c["run"] + ".mosaic_status.tsv"), ("units", c["units"]),
                    ("m6a_selection", f"{P}/m6a_count_model_20261019/runs/freeze/eval_candidates.{c['m6a']}.tsv"),
                    ("m6a_scores", f"{P}/m6a_count_model_20261019/runs/fragscore.{c['m6a']}/fs.scores_truth_free.tsv")]:
        bind(cell, role, p)
    if c["locus"] == "CYP2D6":
        bind(cell, "curation", CURATION)
    if c["alts"]:
        bind(cell, "m1b_alternatives", c["alts"])
    for k in ("m1d", "m1e"):
        if c[k]:
            bind(cell, f"{k}_construction", c[k] + "/candidates.construction.tsv")
            bind(cell, f"{k}_explicit", c[k] + "/candidates.explicit.tsv")
            bind(cell, f"{k}_fasta", c[k] + "/candidates.fa")

    units = {r["name"]: int(r["unit_copies"]) for r in T(c["units"])}
    t = c["t"]
    ucopy = lambda path: units.get(f"allele_{path.split(',')[t]}")
    cand = {r["candidate_id"]: r for r in T(c["run"] + ".mosaic_candidates.tsv")}
    scored = T(c["run"] + ".mosaic_scores.tsv")
    scored_pairs = set()
    for s in scored:
        r = cand[s["representative_beam_candidate"]]
        scored_pairs.add(frozenset([r["hap1_alleles"], r["hap2_alleles"]]) if r["hap1_alleles"] != r["hap2_alleles"] else frozenset([r["hap1_alleles"]]))
    prior_of = {}
    for r in cand.values():
        prior_of[frozenset([r["hap1_alleles"], r["hap2_alleles"]])] = r["exact_path_prior"]
    entries = OrderedDict()

    def add(tier, label, h1, h2, source, expected="NA", F="NA", E="NA", prior="NA"):
        key = frozenset([h1, h2]) if h1 != h2 else frozenset([h1])
        if key in entries:
            e = entries[key]
            e["labels"].append(label)
            e["sources"].append(source)
            if tier == "D":
                e["tier"] = "D"
            for k, v in (("expected_total", expected), ("fragment_term", F), ("exposure_term", E), ("prior", prior)):
                if e[k] == "NA" and v != "NA":
                    e[k] = v
            return
        entries[key] = dict(tier=tier, labels=[label], sources=[source], hap1_alleles=h1, hap2_alleles=h2, expected_total=expected,
                            fragment_term=F, exposure_term=E, prior=prior)

    eq = T(c["run"] + ".equivalence_set.tsv")
    for r in eq:
        add("D", "reported" if r["rank"] == "1" else f"equivalence_rank{r['rank']}", r["hap1_alleles"], r["hap2_alleles"], "caller equivalence set",
            r["score"], r["fragment_term"], r["exposure_term"], r["sequence_pair_prior"])
    if c["alts"]:   # M1-B's saved file: saved_* are truth-free (rule steps 1-2); sub_truth_* are truth-selected; cross_* is another run
        for line in open(c["alts"]):
            if line.startswith("#"):
                continue
            label, h1, h2, e = line.rstrip("\n").split("\t")
            if label == "reported":
                continue
            tier = "D" if label.startswith("saved_") else ("O" if label.startswith("sub_truth") else "X")
            add(tier, label, h1, h2, "M1-B alternatives (" + ("caller best saved per unit class" if tier == "D" else "truth-closest substitution" if tier == "O" else "30x report of the same donor") + ")", e)
    else:   # the same truth-free rule (M1-B Amendment 1 steps 1-2) on the control's saved scored candidates
        rh1, rh2 = eq[0]["hap1_alleles"], eq[0]["hap2_alleles"]
        cls = lambda h1, h2: tuple(sorted((ucopy(h1), ucopy(h2))))
        best = {}
        for s in sorted(scored, key=lambda s: -float(s["score"])):
            r = cand[s["representative_beam_candidate"]]
            k = cls(r["hap1_alleles"], r["hap2_alleles"])
            if k not in best:
                best[k] = (s, r)
        for k in sorted((k for k in best if k != cls(rh1, rh2)), key=lambda k: -float(best[k][0]["score"]))[:7]:
            s, r = best[k]
            add("D", f"saved_{k[0]}+{k[1]}", r["hap1_alleles"], r["hap2_alleles"], "caller best saved per unit class (M1-B rule steps 1-2, truth-free)",
                s["score"], s["fragment_term"], s["exposure_term"], s["sequence_pair_prior"])
    fs = {r["label"]: r for r in T(f"{P}/m6a_count_model_20261019/runs/fragscore.{c['m6a']}/fs.scores_truth_free.tsv")}
    for r in csv.reader(open(f"{P}/m6a_count_model_20261019/runs/freeze/eval_candidates.{c['m6a']}.tsv"), delimiter="\t"):
        label, h1, h2 = r[0], r[1], r[2]
        s = fs.get(label, {})
        add("D", f"m6a_{label}", h1, h2, "M6-A count-model frozen selection (truth-free)", s.get("total", "NA"), s.get("fragment_term", "NA"),
            s.get("exposure_term", "NA"), s.get("prior", "NA"))
    for e in entries.values():
        if True:
            u1, u2 = ucopy(e["hap1_alleles"]), ucopy(e["hap2_alleles"])
            key = frozenset([e["hap1_alleles"], e["hap2_alleles"]]) if e["hap1_alleles"] != e["hap2_alleles"] else frozenset([e["hap1_alleles"]])
            rows.append(dict(cell=cell, locus=c["locus"], target_block=c["block"], tier=e["tier"], labels="|".join(e["labels"]), sources="|".join(sorted(set(e["sources"]))),
                             representation="panel_allele_paths", hap1_alleles=e["hap1_alleles"], hap2_alleles=e["hap2_alleles"], fasta="-", records="-",
                             target_alleles=f"{e['hap1_alleles'].split(',')[t]}|{e['hap2_alleles'].split(',')[t]}",
                             unit_copies=f"{u1}+{u2}", unit_total=(u1 + u2) if None not in (u1, u2) else "NA", unit_class="+".join(str(x) for x in sorted((u1, u2))) if None not in (u1, u2) else "NA",
                             reachability="caller_scored" if key in scored_pairs else "catalogue_member_not_scored_by_caller",
                             copying_prior_defined="yes", saved_expected_total=e["expected_total"], saved_fragment_term=e["fragment_term"],
                             saved_exposure_term=e["exposure_term"], saved_prior=e["prior"]))
    # oracle constructions (explicit sequences; off-panel; no copying prior)
    for k in ("m1d", "m1e"):
        if not c[k]:
            continue
        seqs = {}
        name = None
        for line in open(c[k] + "/candidates.fa"):
            if line.startswith(">"):
                name = line[1:].strip()
                seqs[name] = []
            else:
                seqs[name].append(line.strip())
        seqs = {n: "".join(v) for n, v in seqs.items()}
        explicit = [dict(zip(["label", "fasta", "record1", "record2"], l.rstrip("\n").split("\t")[:4]))
                    for l in open(c[k] + "/candidates.explicit.tsv") if l.strip() and not l.startswith("#")]
        for r in explicit:
            if r["label"] in ("reported", "saved_2+2"):
                continue   # already in D through the caller's own records
            h = [hashlib.sha256(seqs[r[x]].encode()).hexdigest()[:16] for x in ("record1", "record2")]
            rows.append(dict(cell=cell, locus=c["locus"], target_block=c["block"], tier="O", labels=f"{k}_{r['label']}",
                             sources=f"{k.upper()} truth-selected/patched construction" + (" (complete true diplotype; reference only)" if r["label"] == "truth" else ""),
                             representation="explicit_sequences", hap1_alleles="-", hap2_alleles="-", fasta=r["fasta"], records=f"{r['record1']}|{r['record2']}",
                             target_alleles="-", unit_copies="see construction record", unit_total="NA", unit_class="NA", reachability="off_panel",
                             copying_prior_defined="no (prior_support none)", saved_expected_total="NA", saved_fragment_term="NA", saved_exposure_term="NA",
                             saved_prior="NA", sequence_sha256_prefix=f"{h[0]}|{h[1]}"))

os.makedirs(out, exist_ok=True)
keys = list(dict.fromkeys(k for r in rows for k in r))
with open(f"{out}/CANDIDATES.tsv", "w") as f:
    w = csv.DictWriter(f, fieldnames=keys, delimiter="\t", lineterminator="\n", restval="-")
    w.writeheader()
    w.writerows(rows)
with open(f"{out}/BINDINGS.tsv", "w") as f:
    f.write("cell\trole\tfile\tsha256\n")
    for b in bindings:
        f.write("\t".join(b) + "\n")
from collections import Counter
print(Counter((r["cell"], r["tier"]) for r in rows))
for r in rows:
    if r["tier"] != "O":
        print(r["cell"], r["tier"], r["labels"][:60], r["unit_copies"], r["reachability"], r["saved_expected_total"][:14])
