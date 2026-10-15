#!/usr/bin/env python3
"""Adversarial sequence-versus-label fixtures (tests/genotype_sequence_label_reporting_design.md section 4).

Builds four small synthetic loci and permuted variants of each, simulates the donor's reads once per locus, runs the
production caller and records its CURRENT behaviour (baseline_current.tsv) plus an order-independence verdict for each
permuted variant (identical canonicalised outputs required). Nothing here changes the caller.

  F1 pure label ambiguity: unit R (60 bp); bubble A (- / R), context C = R, bubble B (- / R). (R,-) and (-,R) spell the
     same haplotype; both labellings are panel rows. Donor (R,-)/(-,-).
  F2 F1 plus a distant heterozygous bubble E with distinct 80 bp alleles X / Y.
  F3 phase tie: bubbles P (X1 / X2) and Q (Y1 / Y2) 3 kb apart (beyond any fragment); the four phased haplotypes
     appear equally often in the panel. Donor (X1,Y1)/(X2,Y2).
  F4 F1 combined with F3.
Variants: rev (P lines reversed), ren (node ids renumbered), sam (sample names permuted, so their sort order changes).
Canonical comparison: exit, caller status, call class, no-call reason, winner total and fragment term, FASTA sequences
(as a set), block genotypes as (block id, sorted allele lengths, genotype filter), equivalence set as the set of spelled
unordered sequence pairs (spelled from each variant's own model dump).

Every fixture's library file declares marker_background 1e-06 (the value the C4 runs estimate): with random flanks
every invariant anchor has one copy, so background and per-copy rate are otherwise not separately identifiable and
the caller refuses (CANDIDATE_INCOMPLETE). The declared value is identical for every variant. The degeneracy
inventory of every fixture model is recorded (fixture_degeneracies.tsv), so the intended label equivalences are shown
to exist in the caller's own block decomposition rather than assumed.

Usage: fixtures.py <out dir> <panvar> <genotype_probe> <cyp2d6_block_alleles>
"""
import csv
import hashlib
import os
import random
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cyp2d6_pilot"))
import pilot  # noqa: E402

REF = "ref#0#chain"
rng = random.Random(20261015)


def seq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


R, L, RF, SP1, SP2, RF2 = seq(60), seq(1500), seq(1500), seq(3000), seq(3000), seq(1500)
X, Y, X1, X2, Y1, Y2 = seq(80), seq(80), seq(80), seq(80), seq(80), seq(80)


def build(kind):
    """Returns (nodes {id: seq}, paths [(name, [node ids])]); every step '+'."""
    nodes, order = {}, []

    def node(s):
        nid = str(len(nodes) + 1)
        nodes[nid] = s
        return nid

    if kind in ("F1", "F2", "F4"):
        n_l, n_a, n_c, n_b = node(L), node(R), node(R), node(R)
        amb = {"--": [n_l, n_c], "R-": [n_l, n_a, n_c], "-R": [n_l, n_c, n_b], "RR": [n_l, n_a, n_c, n_b]}
    else:
        n_l = node(L)
        amb = {"--": [n_l]}
    if kind in ("F3", "F4"):
        n_sp, n_x1, n_x2, n_sp2, n_y1, n_y2 = node(SP1), node(X1), node(X2), node(SP2), node(Y1), node(Y2)
        ph = {"11": [n_sp, n_x1, n_sp2, n_y1], "22": [n_sp, n_x2, n_sp2, n_y2], "12": [n_sp, n_x1, n_sp2, n_y2],
              "21": [n_sp, n_x2, n_sp2, n_y1]}
    n_rf = node(RF)
    if kind == "F2":
        n_x, n_y, n_rf2 = node(X), node(Y), node(RF2)
        e = {"X": [n_rf, n_x, n_rf2], "Y": [n_rf, n_y, n_rf2]}

    def hap(a="--", p=None, ee=None):
        steps = list(amb[a])
        if p is not None:
            steps += ph[p]
        steps += e[ee] if ee else [n_rf]
        return steps

    if kind == "F1":
        rows = [(REF, hap("--")), ("S1#1#c", hap("R-")), ("S1#2#c", hap("--")), ("S2#1#c", hap("-R")),
                ("S2#2#c", hap("RR")), ("S3#1#c", hap("-R")), ("S3#2#c", hap("--")), ("S4#1#c", hap("R-")),
                ("S4#2#c", hap("RR"))]
    elif kind == "F2":
        rows = [(REF, hap("--", ee="X")), ("S1#1#c", hap("R-", ee="X")), ("S1#2#c", hap("--", ee="Y")),
                ("S2#1#c", hap("-R", ee="X")), ("S2#2#c", hap("RR", ee="Y")), ("S3#1#c", hap("-R", ee="Y")),
                ("S3#2#c", hap("--", ee="X")), ("S4#1#c", hap("R-", ee="Y")), ("S4#2#c", hap("RR", ee="X"))]
    elif kind == "F3":
        rows = [(REF, hap(p="11")), ("S1#1#c", hap(p="11")), ("S1#2#c", hap(p="22")), ("S2#1#c", hap(p="12")),
                ("S2#2#c", hap(p="21")), ("S3#1#c", hap(p="22")), ("S3#2#c", hap(p="12")), ("S4#1#c", hap(p="21"))]
    else:
        rows = [(REF, hap("--", p="11")), ("S1#1#c", hap("R-", p="11")), ("S1#2#c", hap("--", p="22")),
                ("S2#1#c", hap("-R", p="12")), ("S2#2#c", hap("RR", p="21")), ("S3#1#c", hap("-R", p="22")),
                ("S3#2#c", hap("--", p="12")), ("S4#1#c", hap("R-", p="21")), ("S4#2#c", hap("RR", p="11"))]
    return nodes, rows


def write_gfa(path, nodes, rows, rename=None, reverse=False, samples=None):
    rename = rename or {n: n for n in nodes}
    samples = samples or {}
    edges = set()
    for _, steps in rows:
        edges.update(zip(steps, steps[1:]))
    out = ["H\tVN:Z:1.0"] + [f"S\t{rename[n]}\t{s}" for n, s in sorted(nodes.items(), key=lambda t: int(rename[t[0]]))]
    out += [f"L\t{rename[a]}\t+\t{rename[b]}\t+\t0M" for a, b in sorted(edges, key=lambda e: (int(rename[e[0]]), int(rename[e[1]])))]
    plines = []
    for name, steps in rows:
        s = name.split("#")[0]
        new = name.replace(s + "#", samples.get(s, s) + "#", 1) if name != REF else name
        plines.append(f"P\t{new}\t{','.join(rename[n] + '+' for n in steps)}\t*")
    out += list(reversed(plines)) if reverse else plines
    open(path, "w").write("\n".join(out) + "\n")


def kv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path).read().splitlines():
            p = line.split("\t")
            if len(p) >= 2:
                out.setdefault(p[0], p[1])
    return out


def canonical(prefix, dump_prefix, code):
    st = kv(prefix + ".mosaic_status.tsv")
    blocks = {}
    if os.path.exists(dump_prefix + ".block_alleles.tsv"):
        with open(dump_prefix + ".block_alleles.tsv") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                blocks.setdefault(int(r["block_index"]), []).append(r["sequence"])
    spell = lambda a: "".join(blocks[b][int(x)] for b, x in enumerate(a.split(",")))  # noqa: E731
    eq = set()
    if os.path.exists(prefix + ".equivalence_set.tsv"):
        with open(prefix + ".equivalence_set.tsv") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                pair = sorted(hashlib.sha256(spell(r[h]).encode()).hexdigest()[:16] for h in ("hap1_alleles", "hap2_alleles"))
                eq.add(tuple(pair))
    bl = []
    if os.path.exists(prefix + ".block_genotypes.tsv"):
        with open(prefix + ".block_genotypes.tsv") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                bl.append((r["block_id"], tuple(sorted((r["allele1_bp"], r["allele2_bp"]))), r.get("genotype_filter", "")))
    fa = []
    if os.path.exists(prefix + ".mosaic_haplotypes.fa"):
        fa = sorted(hashlib.sha256(l.strip().encode()).hexdigest()[:16]
                    for l in open(prefix + ".mosaic_haplotypes.fa") if l.strip() and not l.startswith(">"))
    return {"exit": str(code), "status": st.get("status", "NA"), "call_class": st.get("call_class", "NA"),
            "no_call_reason": st.get("no_call_reason", "NA"), "note": st.get("note", "")[:140],
            "blocks": str(len(blocks)), "equivalence_set": ";".join("/".join(p) for p in sorted(eq)) or "none",
            "block_genotypes": ";".join(f"{b}:{'|'.join(a)}:{g}" for b, a, g in bl) or "none",
            "fasta": ";".join(fa) or "none"}


def run(cmd, log):
    with open(log, "w") as f:
        return subprocess.run([str(c) for c in cmd], stdout=f, stderr=subprocess.STDOUT).returncode


def main(out, panvar, probe, dumper):
    os.makedirs(out, exist_ok=True)
    profile = os.path.join(out, "profile.tsv")
    open(profile, "w").write("profile\tfixture-bounds\ntime_limit_seconds\t300\nfootprint_limit_bytes\t4294967296\n")
    rows_out = []
    for kind in ("F1", "F2", "F3", "F4"):
        nodes, rows = build(kind)
        d = os.path.join(out, kind)
        os.makedirs(d, exist_ok=True)
        base = os.path.join(d, "base.gfa")
        write_gfa(base, nodes, rows)
        reads = os.path.join(d, "reads")
        os.makedirs(reads, exist_ok=True)
        run([probe, "simulate", "--gfa", base, "--reference", REF, "--donor", "S1", "--depth", "30", "--read-length", "100",
             "--insert-mean", "350", "--insert-sd", "40", "--sigmas", "4", "--error", "0.002", "--outlier", "0.001",
             "--seed", "20261015", "--seed-salt", kind, "--out", reads], os.path.join(d, "simulate.log"))
        with open(os.path.join(reads, "library.txt"), "a") as lib:
            lib.write("marker_background 1e-06\n")
        perm = list(nodes)
        random.Random(7).shuffle(perm)
        variants = {"base": {}, "rev": {"reverse": True},
                    "ren": {"rename": {n: str(100 + i) for i, n in enumerate(perm)}},
                    "sam": {"samples": {"S1": "Z1", "S2": "A9", "S3": "M5", "S4": "B2"}}}
        results = {}
        for name, opts in variants.items():
            gfa = os.path.join(d, f"{name}.gfa")
            if name != "base":
                write_gfa(gfa, nodes, rows, **opts)
            prefix = os.path.join(d, f"{name}.call")
            cmd = [panvar, "genotype", "--preset", "short-read", "--profile", profile, "--threads", "2", "--graph", gfa,
                   "--reference-path", REF, "--reads1", os.path.join(reads, "reads_1.fq"),
                   "--reads2", os.path.join(reads, "reads_2.fq"), "--library", os.path.join(reads, "library.txt"),
                   "--dosage-mode", "repair", "--out", prefix, "--quiet"]
            with open(prefix + ".log", "w") as log:
                code, _ = pilot.run_timed(cmd, log, 360)
            run([dumper, gfa, REF, "none", os.path.join(d, f"{name}.model")], os.path.join(d, f"{name}.model.log"))
            results[name] = canonical(prefix, os.path.join(d, f"{name}.model"), code)
            fp = kv(os.path.join(d, f"{name}.model.summary.tsv")).get("model_fingerprint", "NA")
            row = {"fixture": kind, "variant": name, "model_fingerprint": fp, **results[name]}
            if name != "base":
                same = {k: results[name][k] == results["base"][k] for k in results["base"] if k != "note"}
                row["order_independent"] = "yes" if all(same.values()) else "NO: " + ",".join(k for k, v in same.items() if not v)
            else:
                row["order_independent"] = "reference"
            rows_out.append(row)
            print(kind, name, row["exit"], row["status"], row["call_class"], row["no_call_reason"], row["order_independent"],
                  flush=True)
    inventory = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "cyp2d6_pilot", "degeneracy_inventory.py")
    run([sys.executable, "-I", inventory, os.path.join(out, "fixture_degeneracies.tsv"),
         os.path.join(out, "fixture_degeneracy_summary.tsv")]
        + [f"{k}={os.path.join(out, k, 'base.model')}" for k in ("F1", "F2", "F3", "F4")], os.path.join(out, "inventory.log"))
    with open(os.path.join(out, "baseline_current.tsv"), "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows_out[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows_out)


if __name__ == "__main__":
    main(os.path.abspath(sys.argv[1]), *[os.path.abspath(a) for a in sys.argv[2:5]])
