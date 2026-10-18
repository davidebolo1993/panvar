#!/usr/bin/env python3
"""S1-C: Sol's two lifecycle reproductions through the CLI, with the S1-B binary (before) and the S1-C binary (after),
on fixture F1 (graph and reads from S1-B's fixture directory; read only). Per build, at ONE prefix:
  A  sequence_first_v1 run (PREDICTION) -> legacy run          : which sequence reports survive beside the legacy call
  B  sequence_first_v1 run (PREDICTION) -> watchdog stop          : exit, status, and which sequence reports survive
     (profile time_limit_seconds 0.000001, as in Sol's watchdog.tsv)
  repro_lifecycle.py <out> <before panvar> <after panvar> <F1 fixture dir> <reads dir>
"""
import os
import subprocess
import sys

OWNED = [".sequence_locus.tsv", ".sequence_blocks.tsv", ".sequence_haplotypes.fa", ".representation_intervals.tsv",
         ".sequence_alternatives.tsv"]


def run(panvar, gfa, reads, profile, prefix, v1):
    cmd = [panvar, "genotype", "--preset", "short-read", "--profile", profile, "--threads", "2", "--graph", gfa,
           "--reference-path", "ref#0#chain", "--reads1", f"{reads}/reads_1.fq", "--reads2", f"{reads}/reads_2.fq",
           "--library", f"{reads}/library.txt", "--dosage-mode", "repair", "--out", prefix, "--quiet"]
    if v1:
        cmd += ["--output-schema", "sequence_first_v1"]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=120).returncode


def status(prefix):
    for line in open(prefix + ".mosaic_status.tsv").read().splitlines():
        if line.startswith("status\t"):
            return line.split("\t")[1]
    return "NA"


def main(out, before, after, fixture, reads):
    os.makedirs(out, exist_ok=True)
    normal = os.path.join(out, "profile.tsv")
    watchdog = os.path.join(out, "watchdog.tsv")
    open(normal, "w").write("profile\tfixture-bounds\ntime_limit_seconds\t300\nfootprint_limit_bytes\t4294967296\n")
    open(watchdog, "w").write("profile\ts1c-watchdog\ntime_limit_seconds\t0.000001\nfootprint_limit_bytes\t4294967296\n")
    gfa = os.path.join(fixture, "base.gfa")
    rows = ["build\tsequence\tfirst_exit\tsecond_exit\tsecond_status\tsurviving_sequence_reports\tlocus_prediction_left"]
    bad = 0
    for label, panvar in (("before_S1B", before), ("after_S1C", after)):
        for seq, second_v1, second_profile in (("A_prediction_then_legacy", False, normal),
                                               ("B_prediction_then_watchdog", True, watchdog)):
            d = os.path.join(out, label, seq)
            os.makedirs(d, exist_ok=True)
            prefix = os.path.join(d, "call")
            first = run(panvar, gfa, reads, normal, prefix, True)
            second = run(panvar, gfa, reads, second_profile, prefix, second_v1)
            survivors = [s for s in OWNED if os.path.exists(prefix + s)]
            left = "none"
            if os.path.exists(prefix + ".sequence_locus.tsv"):
                for line in open(prefix + ".sequence_locus.tsv").read().splitlines():
                    if line.startswith("prediction\t"):
                        left = line.split("\t")[1]
            rows.append(f"{label}\t{seq}\t{first}\t{second}\t{status(prefix)}\t{','.join(survivors) or 'none'}\t{left}")
            if label == "after_S1C" and survivors:
                bad += 1
    open(os.path.join(out, "repro_lifecycle.tsv"), "w").write("\n".join(rows) + "\n")
    print("\n".join(rows))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:6]))
