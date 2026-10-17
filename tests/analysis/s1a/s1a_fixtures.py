#!/usr/bin/env python3
"""S1-A executable fixture evidence: CURRENT caller outputs against the PROPOSED sequence-first outputs.

Reuses the representation fixtures F1-F4 (tests/analysis/representation_fixtures/fixtures.py; graphs rebuilt
bit-identically from the same seeds, reads reused read-only from tests/results/representation_fixtures_20261015) and
adds only the missing cases, each on DETERMINISTIC error-free tiled reads (fixture construction; no simulator, no
randomness in the reads):
  F5 partial resolution: bubble C with three equal-length alleles that no fragment reaches (reads omitted around C),
     symmetric panel, plus an observed heterozygous bubble E -> a genuine three-way tie of distinct sequences that differ
     at C, with E identifiable.
  F7 provisional: three linked bubbles, a two-haplotype panel and an off-panel donor (phase switched at the middle
     bubble), union = panel only, repair_rounds 1 -> the search stops at its round cap with a valid, exactly scored
     incumbent. F7 'full' is the same input with the default round cap (a completed control).
F2 already holds a lower-ranked projection-ambiguous sequence group with substantial probability (0.336 at the 2026-10-15
baseline), so no separate fixture is added for that case. Variants rev / ren / sam (P-line order, node ids, sample names)
are run for every fixture; the proposed outputs must be identical across them.
  s1a_fixtures.py <out dir> <panvar> <cyp2d6_block_alleles>
"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "representation_fixtures"))
sys.path.insert(0, os.path.join(HERE, "..", "cyp2d6_pilot"))
sys.path.insert(0, HERE)
import fixtures  # noqa: E402  (module-level sequences come from its fixed seed: identical graphs)
import pilot  # noqa: E402
import random  # noqa: E402
import s1a_project  # noqa: E402

REF = fixtures.REF
OLD = "/Users/davide.bolognini/panvar/tests/results/representation_fixtures_20261015"
rng = random.Random(20261017)


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


L5, C1, C2, C3, SP5, E1, E2, RF5 = seq(1500), seq(80), seq(80), seq(80), seq(3000), seq(80), seq(80), seq(1500)
L7, X1, X2, S71, Y1, Y2, S72, Z1, Z2, RF7 = seq(1500), seq(80), seq(80), seq(150), seq(80), seq(80), seq(150), seq(80), seq(80), seq(1500)


def build(kind):
    nodes = {}

    def node(s):
        nid = str(len(nodes) + 1)
        nodes[nid] = s
        return nid
    if kind == "F5":
        nl, c = node(L5), [node(C1), node(C2), node(C3)]
        nsp, e, nrf = node(SP5), [node(E1), node(E2)], node(RF5)
        hap = lambda i, j: [nl, c[i - 1], nsp, e[j - 1], nrf]  # noqa: E731
        rows = [(REF, hap(1, 1)), ("S1#1#c", hap(2, 2)), ("S1#2#c", hap(3, 1)), ("S2#1#c", hap(1, 2)),
                ("S2#2#c", hap(2, 1)), ("S3#1#c", hap(3, 2))]
        donor = [hap(1, 1), hap(2, 2)]
        gap = (len(L5), len(L5) + 80)   # C on every haplotype
    else:   # F7
        nl, x, s1, y, s2, z, nrf = node(L7), [node(X1), node(X2)], node(S71), [node(Y1), node(Y2)], node(S72), [node(Z1), node(Z2)], node(RF7)
        hap = lambda a, b, cc: [nl, x[a - 1], s1, y[b - 1], s2, z[cc - 1], nrf]  # noqa: E731
        rows = [(REF, hap(1, 1, 1)), ("S1#1#c", hap(2, 2, 2)), ("S1#2#c", hap(1, 1, 1)), ("S2#1#c", hap(2, 2, 2)),
                ("S2#2#c", hap(1, 1, 1)), ("S3#1#c", hap(2, 2, 2))]
        donor = [hap(1, 2, 1), hap(2, 1, 2)]
        gap = None
    return nodes, rows, donor, gap


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def tiled_reads(nodes, donor, gap, out, step=13, insert=350, length=100):
    """Deterministic, error-free FR pairs every `step` bp on each donor haplotype; fragments touching `gap` are omitted."""
    os.makedirs(out, exist_ok=True)
    n = 0
    with open(os.path.join(out, "reads_1.fq"), "w") as r1, open(os.path.join(out, "reads_2.fq"), "w") as r2:
        for h, steps in enumerate(donor):
            s = "".join(nodes[x] for x in steps)
            for start in range(0, len(s) - insert + 1, step):
                if gap and start < gap[1] + 20 and start + insert > gap[0] - 20:
                    continue
                a, b = s[start:start + length], revcomp(s[start + insert - length:start + insert])
                r1.write(f"@t{n}/1\n{a}\n+\n{'I' * length}\n")
                r2.write(f"@t{n}/2\n{b}\n+\n{'I' * length}\n")
                n += 1
    open(os.path.join(out, "library.txt"), "w").write(
        "# S1-A deterministic tiled reads (error-free, every 13 bp per haplotype, insert 350)\ninsert_mean_bp 350\n"
        "insert_sd_bp 40\nsupport_sigmas 4\nsubstitution_error 0.002\noutlier_weight 0.001\n"
        f"fragment_start_intensity {1.0 / step:.7g}\nmarker_background 1e-06\n")
    return n


def main(out, panvar, dumper):
    os.makedirs(out, exist_ok=True)
    base_profile = "profile\tfixture-bounds\ntime_limit_seconds\t300\nfootprint_limit_bytes\t4294967296\n"
    open(os.path.join(out, "profile.tsv"), "w").write(base_profile)
    open(os.path.join(out, "profile_round1.tsv"), "w").write(base_profile + "repair_rounds\t1\nmax_mosaic_candidates\t0\n")
    open(os.path.join(out, "profile_k0.tsv"), "w").write(base_profile + "max_mosaic_candidates\t0\n")
    summary, block_rows = [], []
    plan = [("F1", None), ("F2", None), ("F3", None), ("F4", None), ("F5", None), ("F7", "profile_round1.tsv"), ("F7full", "profile_k0.tsv")]
    for kind, profile in plan:
        d = os.path.join(out, kind)
        os.makedirs(d, exist_ok=True)
        if kind in ("F1", "F2", "F3", "F4"):
            nodes, rows = fixtures.build(kind)
            reads = os.path.join(OLD, kind, "reads")   # reused read-only
        else:
            nodes, rows, donor, gap = build("F7" if kind.startswith("F7") else kind)
            reads = os.path.join(d, "reads")
            tiled_reads(nodes, donor, gap, reads)
        perm = list(nodes)
        random.Random(7).shuffle(perm)
        variants = {"base": {}, "rev": {"reverse": True}, "ren": {"rename": {n: str(100 + i) for i, n in enumerate(perm)}},
                    "sam": {"samples": {"S1": "Z1", "S2": "A9", "S3": "M5", "S4": "B2"}}}
        canon = {}
        for name, opts in variants.items():
            gfa = os.path.join(d, f"{name}.gfa")
            fixtures.write_gfa(gfa, nodes, rows, **opts)
            prefix = os.path.join(d, f"{name}.call")
            cmd = [panvar, "genotype", "--preset", "short-read", "--profile", os.path.join(out, profile or "profile.tsv"),
                   "--threads", "2", "--graph", gfa, "--reference-path", REF, "--reads1", os.path.join(reads, "reads_1.fq"),
                   "--reads2", os.path.join(reads, "reads_2.fq"), "--library", os.path.join(reads, "library.txt"),
                   "--dosage-mode", "repair", "--out", prefix, "--quiet"]
            with open(prefix + ".log", "w") as log:
                code, _ = pilot.run_timed(cmd, log, 300)
            fixtures.run([dumper, gfa, REF, "none", os.path.join(d, f"{name}.model")], os.path.join(d, f"{name}.model.log"))
            cur, prop = s1a_project.project(prefix, os.path.join(d, f"{name}.model"))
            key = (cur["status"], prop["prediction"], prop["sequence_status"], prop["provisional"], prop["representative"],
                   tuple((b[0], b[1], b[2]) for b in prop["blocks"]), tuple(prop["intervals"]))
            canon[name] = key
            summary.append({"fixture": kind, "variant": name, "exit": code, "current_status": cur["status"],
                            "current_call_class": cur["call_class"], "current_no_call_reason": cur["no_call_reason"],
                            "current_fasta_records": cur["fasta_records"], "current_identified_block_rows": cur["identified_block_rows"],
                            "proposed_prediction": prop["prediction"], "proposed_sequence_status": prop["sequence_status"],
                            "proposed_provisional": prop["provisional"], "proposed_representative": prop["representative"],
                            "proposed_fasta": prop["fasta"], "proposed_identified_blocks": sum(1 for b in prop["blocks"] if b[1] == "IDENTIFIED"),
                            "proposed_ambiguous_blocks": sum(1 for b in prop["blocks"] if b[1] == "REPRESENTATION_AMBIGUOUS"),
                            "proposed_unresolved_blocks": sum(1 for b in prop["blocks"] if b[1] == "GENOTYPE_UNRESOLVED"),
                            "proposed_intervals": ";".join(prop["intervals"]) or "none",
                            "pool_conditional_posterior": prop["pool_conditional_posterior"],
                            "structure_confidence": prop["structure_confidence"],
                            "lower_projection_ambiguous_groups": prop.get("lower_projection_ambiguous_groups", "NA"),
                            "projections_unexplored_by_caller": prop.get("projections_unexplored_by_caller", "NA"),
                            "alternatives": " | ".join(prop["alternatives"]) or "none",
                            "reason": prop.get("reason", ""),
                            "order_independent": "reference" if name == "base" else ("yes" if canon[name] == canon["base"] else "NO")})
            if name == "base":
                for b in prop["blocks"]:
                    block_rows.append({"fixture": kind, "block": b[0], "proposed_status": b[1], "genotype_or_alternatives": b[2],
                                       "error_mass_lower": b[3], "error_mass_upper": b[4], "current_conditional_GQ": b[5]})
            print(kind, name, cur["status"], "->", prop["prediction"], prop["sequence_status"], prop["provisional"],
                  summary[-1]["order_independent"], flush=True)
    for name, data in (("current_vs_proposed.tsv", summary), ("proposed_block_rows.tsv", block_rows)):
        with open(os.path.join(out, name), "w") as f:
            w = csv.DictWriter(f, fieldnames=list(data[0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(data)


if __name__ == "__main__":
    main(os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), os.path.abspath(sys.argv[3]))
