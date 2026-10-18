#!/usr/bin/env python3
"""S1-D fixture gates for output_schema sequence_first_v2 (caller level, deterministic tiled reads).

Fixtures (graphs built here or by the S1-A/B builders; reads error-free and tiled, as in S1-A):
  F3     two panel haplotypes (S1-A F3; reads reused read-only)
  F9     homozygous panel pair (S1-B F9)
  F1     sequence-equivalent aliases and alternative labellings of one sequence (S1-A F1; reads reused)
  F7full two off-panel haplotypes (S1-A F7 graph, off-panel donor, k0)
  F7     provisional (round cap 1), off-panel
  F11    three bubbles; donor (1,1,1)/(2,1,2); panel rows (1,1,1), (2,2,2), (1,2,2), (2,1,2), (1,1,1).
         variant 'uncurated': both haplotypes panel-represented (aliases REF and S2#2 for (1,1,1)).
         variant 'curated' (S2#1#c, the only (2,1,2) template, curated out): one panel, one off-panel -> FASTA.
  F12    F11 'uncurated' with S2#1#c written in REVERSE orientation (steps reversed, '-'): the panel match is reported
         reverse_complement and the extraction must reverse-complement the path.
Per case: membership is checked independently (a haplotype is PANEL iff an active path's sequence, or its reverse
complement, equals it -- extracted from the GFA, curated paths excluded); every PANEL row is extracted with
extract_panel_haplotype.py --check; the FASTA decision is checked; legacy identity (unset schema, --dosage-mode audit)
is checked on F11; F11 'uncurated' is also run through a persisted prepared-locus index (write-panel, loaded) and must
equal its cold outputs; one prefix runs F11 curated (off-panel) -> F11 uncurated (panel only) -> a watchdog stop.
  s1d_fixtures.py <out dir> <panvar> <genotype_prepared_locus_driver> [F3,F9,...]
"""
import csv
import hashlib
import os
import random
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
for sub in ("../representation_fixtures", "../cyp2d6_pilot", "../s1a", "../s1b", "."):
    sys.path.insert(0, os.path.join(HERE, sub))
import fixtures  # noqa: E402
import pilot  # noqa: E402
import s1a_fixtures  # noqa: E402
import extract_panel_haplotype as ex  # noqa: E402

REF = fixtures.REF
OLD = "/Users/davide.bolognini/panvar/tests/results/representation_fixtures_20261015"
rng = random.Random(20261018)
SEQS = ["".join(rng.choice("ACGT") for _ in range(n)) for n in (1500, 80, 80, 150, 80, 80, 150, 80, 80, 1500)]
L, X1, X2, S1, Y1, Y2, S2, Z1, Z2, RF = SEQS


def f11(reverse_s21=False):
    nodes = {str(i + 1): s for i, s in enumerate(SEQS)}
    nl, x, s1, y, s2, z, nrf = "1", ["2", "3"], "4", ["5", "6"], "7", ["8", "9"], "10"
    hap = lambda a, b, c: [nl, x[a - 1], s1, y[b - 1], s2, z[c - 1], nrf]  # noqa: E731
    rows = [(REF, hap(1, 1, 1)), ("S1#1#c", hap(2, 2, 2)), ("S1#2#c", hap(1, 2, 2)), ("S2#1#c", hap(2, 1, 2)),
            ("S2#2#c", hap(1, 1, 1))]
    donor = [hap(1, 1, 1), hap(2, 1, 2)]
    return nodes, rows, donor, {"S2#1#c"} if reverse_s21 else set()


def write_gfa(path, nodes, rows, reversed_paths=()):
    edges = set()
    plines = []
    for name, steps in rows:
        if name in reversed_paths:
            st = [(n, "-") for n in reversed(steps)]
        else:
            st = [(n, "+") for n in steps]
        for (a, oa), (b, ob) in zip(st, st[1:]):
            edges.add((a, oa, b, ob))
        plines.append(f"P\t{name}\t{','.join(n + o for n, o in st)}\t*")
    out = ["H\tVN:Z:1.0"] + [f"S\t{n}\t{s}" for n, s in sorted(nodes.items(), key=lambda t: int(t[0]))]
    out += [f"L\t{a}\t{oa}\t{b}\t{ob}\t0M" for a, oa, b, ob in sorted(edges)]
    open(path, "w").write("\n".join(out + plines) + "\n")


def run(cmd, log, timeout=300):
    with open(log, "w") as f:
        f.write("$ " + " ".join(cmd) + "\n")
        f.flush()
        code, _ = pilot.run_timed(cmd, f, timeout)
    return code


def kv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path).read().splitlines():
            p = line.split("\t", 1)
            if len(p) == 2:
                out[p[0]] = p[1]
    return out


def membership_check(gfa, hap_table, curated):
    """Independent: PANEL iff an active path (curated excluded) spells the haplotype in either orientation."""
    segments, paths = ex.load(gfa)
    active = {}
    for name in paths:
        if name in curated:
            continue
        seq = ex.path_sequence(segments, paths, name)
        for s in (seq, ex.revcomp(seq)):
            active.setdefault(hashlib.sha256(s.encode()).hexdigest(), set()).add(name)
    rows = list(csv.DictReader(open(hap_table), delimiter="\t"))
    problems = []
    for r in rows:
        want = sorted(active.get(r["sha256"], set()))
        got = [] if r["panel_paths"] == "NA" else r["panel_paths"].split(",")
        if (r["origin"] == "PANEL") != bool(want) or sorted(got) != want:
            problems.append(f"{r['haplotype']}: {r['origin']} {got} vs independent {want}")
    return rows, problems


def main(out, panvar, driver, only=None):
    os.makedirs(out, exist_ok=True)
    base = "profile\tfixture-bounds\ntime_limit_seconds\t300\nfootprint_limit_bytes\t4294967296\n"
    profiles = {"default": base, "round1": base + "repair_rounds\t1\nmax_mosaic_candidates\t0\n",
                "k0": base + "max_mosaic_candidates\t0\n",
                "watchdog": "profile\ts1d-watchdog\ntime_limit_seconds\t0.000001\nfootprint_limit_bytes\t4294967296\n"}
    for name, text in profiles.items():
        open(os.path.join(out, f"profile_{name}.tsv"), "w").write(text)
    open(os.path.join(out, "curation_F11.txt"), "w").write("# S1-D: the only (2,1,2) template\nS2#1#c\n")
    plan = [("F3", "default", None, "PANEL_REPRESENTED_PAIR"), ("F9", "k0", None, "PANEL_REPRESENTED_PAIR"),
            ("F1", "default", None, "PANEL_REPRESENTED_PAIR"), ("F7full", "k0", None, "OFF_PANEL_PAIR"),
            ("F7", "round1", None, "OFF_PANEL_PAIR"), ("F11u", "k0", None, "PANEL_REPRESENTED_PAIR"),
            ("F11c", "k0", "curation_F11.txt", "OFF_PANEL_PAIR"), ("F12", "k0", None, "PANEL_REPRESENTED_PAIR")]
    if only:
        plan = [p for p in plan if p[0] in only]
    summary = []
    for kind, profile, curation, expect in plan:
        d = os.path.join(out, kind)
        os.makedirs(d, exist_ok=True)
        gfa = os.path.join(d, "graph.gfa")
        reversed_paths = set()
        if kind in ("F1", "F3"):
            nodes, rows = fixtures.build(kind)
            fixtures.write_gfa(gfa, nodes, rows)
            reads = os.path.join(OLD, kind, "reads")
        elif kind == "F9":
            import s1b_fixtures
            nodes, rows, donor, gap, _ = s1b_fixtures.build("F9")
            fixtures.write_gfa(gfa, nodes, rows)
            reads = os.path.join(d, "reads")
            s1a_fixtures.tiled_reads(nodes, donor, gap, reads)
        elif kind.startswith("F7"):
            nodes, rows, donor, gap = s1a_fixtures.build("F7")
            fixtures.write_gfa(gfa, nodes, rows)
            reads = os.path.join(d, "reads")
            s1a_fixtures.tiled_reads(nodes, donor, gap, reads)
        else:
            nodes, rows, donor, reversed_paths = f11(reverse_s21=(kind == "F12"))
            write_gfa(gfa, nodes, rows, reversed_paths)
            reads = os.path.join(d, "reads")
            s1a_fixtures.tiled_reads(nodes, donor, None, reads)
        args = ["--preset", "short-read", "--profile", os.path.join(out, f"profile_{profile}.tsv"), "--threads", "2",
                "--graph", gfa, "--reference-path", REF, "--reads1", os.path.join(reads, "reads_1.fq"), "--reads2",
                os.path.join(reads, "reads_2.fq"), "--library", os.path.join(reads, "library.txt"), "--quiet"]
        if curation:
            args[args.index("--reads1"):args.index("--reads1")] = ["--panel-curation", os.path.join(out, curation)]
        v2 = os.path.join(d, "v2")
        os.makedirs(v2, exist_ok=True)
        code = run([panvar, "genotype", *args, "--output-schema", "sequence_first_v2", "--out", f"{v2}/call"],
                   os.path.join(d, "v2.log"))
        loc = kv(f"{v2}/call.sequence_locus.tsv")
        hap_rows, problems = membership_check(gfa, f"{v2}/call.sequence_haplotypes.tsv",
                                              {"S2#1#c"} if curation else set())
        check = subprocess.run([sys.executable, os.path.join(HERE, "extract_panel_haplotype.py"), "--check", gfa,
                                f"{v2}/call.sequence_haplotypes.tsv"], capture_output=True, text=True)
        extraction = check.stdout.strip().splitlines()[-1] if check.stdout.strip() else "no output"
        fasta_ok = (loc.get("fasta_reason") == expect and
                    os.path.exists(f"{v2}/call.sequence_haplotypes.fa") == (expect == "OFF_PANEL_PAIR"))
        orientations = ";".join(sorted({o for r in hap_rows for o in r["path_orientations"].split(",") if o != "NA"}))
        row = {"fixture": kind, "exit": code, "prediction": loc.get("prediction"), "sequence_status": loc.get("sequence_status"),
               "hap_origins": f"{loc.get('hap1_origin')}/{loc.get('hap2_origin')}", "fasta_reason": loc.get("fasta_reason"),
               "expected_fasta_reason": expect, "fasta_decision_ok": "yes" if fasta_ok else "NO",
               "membership_independent": "agree" if not problems else " | ".join(problems),
               "extraction": extraction if check.returncode == 0 else "FAILED: " + extraction,
               "orientations": orientations or "NA", "dosage_mode": loc.get("dosage_mode"),
               "dosage_mode_source": loc.get("dosage_mode_source"), "extra": "-"}
        if kind == "F11u":
            # legacy identity: unset schema with the same (audit) dosage mode
            leg = os.path.join(d, "legacy")
            os.makedirs(leg, exist_ok=True)
            run([panvar, "genotype", *args, "--dosage-mode", "audit", "--out", f"{leg}/call"], os.path.join(d, "legacy.log"))
            diffs = []
            v2_files = {f for f in os.listdir(v2) if f.startswith("call.")}
            leg_files = {f for f in os.listdir(leg) if f.startswith("call.")}
            seq_files = {f for f in v2_files if ".sequence_" in f or ".representation_intervals" in f}
            if v2_files - seq_files != leg_files:
                diffs.append(f"file sets {sorted((v2_files - seq_files) ^ leg_files)}")
            declared = {"output_schema", "sequence_prediction", "sequence_status", "sequence_completion", "sequence_exit_code",
                        "sequence_publication", "sequence_fasta", "sequence_fasta_reason"}
            for f in sorted(leg_files & v2_files):
                a, b = open(f"{leg}/{f}", errors="replace").read(), open(f"{v2}/{f}", errors="replace").read()
                a, b = a.replace(leg, "<D>"), b.replace(v2, "<D>")
                if f.endswith("_status.tsv"):
                    ka, kb = kv(f"{leg}/{f}"), kv(f"{v2}/{f}")
                    for k in declared:
                        kb.pop(k, None)
                    vol = ("seconds", "footprint", "rss", "peak", "_bytes", "worker", "threads_", "profile_sha256")
                    ka = {k: v.replace(leg, "<D>") for k, v in ka.items() if not any(x in k for x in vol)}
                    kb = {k: v.replace(v2, "<D>") for k, v in kb.items() if not any(x in k for x in vol)}
                    if ka != kb:
                        diffs.append(f"{f}: {sorted(set(ka.items()) ^ set(kb.items()))[:3]}")
                elif f == "call.settings.tsv":
                    if b.replace("output_schema\tsequence_first_v2\tflag --output-schema\n", "", 1) != a:
                        diffs.append("settings")
                elif f in ("call.repair_memory.tsv",):
                    continue
                elif f == "call.repair_rounds.tsv":
                    cut = lambda t: [l.split("\t")[:28] for l in t.splitlines()]  # noqa: E731
                    if cut(a) != cut(b):
                        diffs.append("trajectory")
                elif a != b:
                    diffs.append(f"{f} differs")
            # loaded prepared-locus route (persisted index) against the cold v2 outputs
            idx = os.path.join(d, "F11u.panel.pvpl")
            if os.path.exists(idx):
                os.remove(idx)
            ld = os.path.join(d, "v2.loaded")
            os.makedirs(ld, exist_ok=True)
            argfile = os.path.join(d, "v2.loaded.args.txt")
            open(argfile, "w").write("".join(a + "\n" for a in args + ["--output-schema", "sequence_first_v2",
                                                                         "--out", f"{ld}/call"]))
            w = run([driver, "write-panel", str(4 << 30), "0", os.path.join(d, "drv_write"), idx, argfile],
                    os.path.join(d, "index_write.log"))
            lc = run([driver, "loaded", str(4 << 30), "0", os.path.join(d, "drv_load"), idx, argfile],
                     os.path.join(d, "index_load.log"))
            same = all(open(f"{v2}/{f}").read().replace(v2, "<D>") == open(f"{ld}/{f}").read().replace(ld, "<D>")
                       for f in ("call.sequence_locus.tsv", "call.sequence_haplotypes.tsv", "call.sequence_blocks.tsv",
                                 "call.representation_intervals.tsv", "call.sequence_alternatives.tsv"))
            row["extra"] = (f"legacy identity: {'identical' if not diffs else ' | '.join(diffs)}; "
                            f"loaded index: write {w} load {lc} {'identical' if same else 'DIFFERENT'}")
        summary.append(row)
        print("\t".join(str(v) for v in row.values()), flush=True)
    if not only or "F11c" in only:
        # one prefix: F11 curated (off-panel) -> F11 uncurated (panel only) -> watchdog stop
        nodes, rows, donor, _ = f11()
        d = os.path.join(out, "lifecycle")
        os.makedirs(d, exist_ok=True)
        gfa = os.path.join(d, "graph.gfa")
        write_gfa(gfa, nodes, rows)
        reads = os.path.join(d, "reads")
        s1a_fixtures.tiled_reads(nodes, donor, None, reads)
        prefix = os.path.join(d, "call")
        steps = []
        for label, extra, profile in (("off_panel", ["--panel-curation", os.path.join(out, "curation_F11.txt")], "k0"),
                                      ("panel_only", [], "k0"), ("watchdog", [], "watchdog")):
            args = ["--preset", "short-read", "--profile", os.path.join(out, f"profile_{profile}.tsv"), "--threads", "2",
                    "--graph", gfa, "--reference-path", REF, *extra, "--reads1", os.path.join(reads, "reads_1.fq"),
                    "--reads2", os.path.join(reads, "reads_2.fq"), "--library", os.path.join(reads, "library.txt"),
                    "--output-schema", "sequence_first_v2", "--out", prefix, "--quiet"]
            code = run([panvar, "genotype", *args], os.path.join(d, f"{label}.log"))
            present = sorted(f.replace("call", "") for f in os.listdir(d) if f.startswith("call.sequence") or
                             f.startswith("call.representation"))
            steps.append(f"{label}: exit {code}, fasta {'present' if os.path.exists(prefix + '.sequence_haplotypes.fa') else 'absent'}, "
                         f"status {kv(prefix + '.mosaic_status.tsv').get('status')}, sequence files {len(present)}")
        ok = ("fasta present" in steps[0] and "fasta absent" in steps[1] and steps[2].endswith("sequence files 0") and
              "TIME_LIMIT" in steps[2])
        print("lifecycle\t" + " | ".join(steps) + f"\t{'OK' if ok else 'FAILED'}", flush=True)
        summary.append({"fixture": "lifecycle", "exit": "-", "prediction": "-", "sequence_status": "-", "hap_origins": "-",
                        "fasta_reason": "-", "expected_fasta_reason": "-", "fasta_decision_ok": "yes" if ok else "NO",
                        "membership_independent": "-", "extraction": "-", "orientations": "-", "dosage_mode": "-",
                        "dosage_mode_source": "-", "extra": " | ".join(steps)})
    path = os.path.join(out, "fixtures.tsv")
    kept = [r for r in csv.DictReader(open(path), delimiter="\t")] if os.path.exists(path) else []
    kept = [r for r in kept if r["fixture"] not in {s["fixture"] for s in summary}]
    with open(path, "w") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(kept + summary)
    bad = [r for r in summary if r["fasta_decision_ok"] != "yes" or (r["membership_independent"] not in ("agree", "-")) or
           r["extraction"].startswith("FAILED") or ("DIFFERENT" in r["extra"]) or
           ("legacy identity" in r["extra"] and "legacy identity: identical" not in r["extra"])]
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]), os.path.abspath(sys.argv[3]),
                  sys.argv[4].split(",") if len(sys.argv) > 4 else None))
