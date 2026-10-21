#!/usr/bin/env python3
"""M15-E: the public panel-index flags on the development cells (C4 NA19705, CYP2D6 HG02769, ACOT HG03209; a second C4
sample NA19650). No reserved family. Each example asserts the frozen expectations (DESIGN.md section 2) and records
times and sizes, never checking them against a speedup. Usage: examples.py <example> <unit dir>"""
import csv
import hashlib
import os
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "m15d"))
import examples as E  # noqa: E402

R = E.R
OUT = f"{R}/m15e_panel_index_20261021/runs"
C4B = f"{R}/genotype_c4_production_11fe681/reads/NA19650/cov5"
CHECKS = []
# Science outputs compared byte for byte (prefix-derived names normalised); the status and settings carry declared
# provenance and timing differences only (compared below with named allowances).
SCIENCE = [".sequence_haplotypes.tsv", ".sequence_blocks.tsv", ".equivalence_set.tsv", ".block_genotypes.tsv", ".genotypes.tsv",
           ".mosaic_scores.tsv", ".mosaic_candidates.tsv", ".representation_intervals.tsv", ".sequence_alternatives.tsv",
           ".sequence_compound.tsv", ".library_sensitivity_blocks.tsv"]


def check(ok, what):
    CHECKS.append((ok, what))
    print(("PASS " if ok else "FAIL ") + what)


def run(a):
    print("$", " ".join(a).replace(R, "<results>").replace(E.WT, "<worktree>"))
    t = time.time()
    p = subprocess.run(a, capture_output=True, text=True)
    p.wall = time.time() - t
    print(f"exit {p.returncode}  wall {p.wall:.2f} s")
    if p.stderr.strip():
        print("stderr:", p.stderr.strip()[-900:])
    return p


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def named(prefix, suffix):
    """A file's text with the run's own prefix name replaced (FASTA record names, the locus table's sibling links)."""
    f = prefix + suffix
    if not os.path.exists(f):
        return None
    return open(f).read().replace(prefix, "<RUN>").replace(os.path.basename(prefix), "<RUN>")


def index_path(cell):
    return f"{OUT}/ex.build_{cell}/index/{cell}.pvpl"


def locus_cmd(cell, extra=()):
    c = E.CELLS[cell]
    a = [E.PANVAR, "genotype", "--graph", c["graph"], "--reference-path", c["ref"]]
    if c["cur"]:
        a += ["--panel-curation", c["cur"]]
    return a + list(extra)


def build(cell, d):
    d = f"{d}/index"   # the unit directory also holds the runner's log
    os.makedirs(d, exist_ok=True)
    idx = f"{d}/{cell}.pvpl"
    p = run(locus_cmd(cell, ["--build-panel-index", idx]))
    open(f"{d}/{cell}.build.tsv", "w").write(p.stdout)
    rec = dict(l.split("\t", 1) for l in p.stdout.strip().split("\n")[1:] if "\t" in l)
    rec = {k: v.strip() for k, v in rec.items()}
    check(p.returncode == 0 and os.path.exists(idx), f"{cell}: built without reads or a library")
    check(sorted(os.listdir(d)) == sorted([f"{cell}.pvpl", f"{cell}.build.tsv"]), f"{cell}: the build wrote the index only ({sorted(os.listdir(d))})")
    check(all(rec.get(k) == "1" for k in ("schema_version", "preparation_algorithm_version", "panel_extension_version",
                                          "marker_panel_algorithm_version")) and rec.get("identity_gfa_sha256") == sha(E.CELLS[cell]["graph"]),
          f"{cell}: the record states the versions and the graph identity")
    check(rec.get("file_sha256") == sha(idx) and int(rec.get("file_bytes", "0")) == os.path.getsize(idx), f"{cell}: the record states the file's size and SHA-256")
    print(f"index {os.path.getsize(idx)} bytes; prepare {rec.get('prepare_seconds')} s, write {rec.get('write_seconds')} s; "
          f"curation {rec.get('identity_curation')}; fingerprint {rec.get('model_fingerprint')}")
    before = sha(idx)
    p2 = run(locus_cmd(cell, ["--build-panel-index", idx]))
    check(p2.returncode == 2 and "(exists)" in p2.stderr and sha(idx) == before, f"{cell}: an existing destination is refused, byte-identical")
    return 0


def ex_build_c4(d):
    return build("c4", d)


def ex_build_cyp2d6(d):
    return build("cyp2d6", d)


def ex_build_acot(d):
    return build("acot", d)


def compare(cell, cold, idx, what, dosage=False):
    for suffix in SCIENCE + ([".dosage_blocks.tsv", ".dosage_census.tsv"] if dosage else []):
        a, b = named(cold, suffix), named(idx, suffix)
        check(a == b, f"{what}: {suffix} identical ({'present' if a is not None else 'absent in both' if b is None else 'absent in one'})")
    check(E.fasta_key(cold + ".sequence_haplotypes.fa") == E.fasta_key(idx + ".sequence_haplotypes.fa"), f"{what}: FASTA identical (record names aside)")
    check(E.untimed(cold + ".library_estimate.tsv") == E.untimed(idx + ".library_estimate.tsv"), f"{what}: library estimate identical (timing aside)")
    check(E.prediction_key(cold) == E.prediction_key(idx), f"{what}: prediction identical (score, paths, haplotype hashes, status, fingerprint)")
    check(named(cold, ".locus_genotype.tsv") == named(idx, ".locus_genotype.tsv"), f"{what}: locus table identical (prefix names aside)")
    sc, si = E.settings(cold), E.settings(idx)
    check({k: v[0] for k, v in sc.items()} == {k: v[0] for k, v in si.items()}, f"{what}: every effective setting identical")
    st = E.KV(idx + ".mosaic_status.tsv")
    check(st.get("prepared_locus_origin") == "loaded" and st.get("prepared_locus_panel_source") == "loaded" and
          st.get("panel_index_panel_rebuild_reason") == "none", f"{what}: the indexed run loaded the model AND the panel")
    if dosage:
        x, y = E.untimed(cold + ".dosage_status.tsv"), E.untimed(idx + ".dosage_status.tsv")
        check(x is not None and x == y, f"{what}: dosage status identical (timing aside)")
    tc, ti = E.KV(cold + ".mosaic_status.tsv"), st
    print(f"{what}: seconds_total cold {tc.get('seconds_total')} vs indexed {ti.get('seconds_total')}; index load "
          f"{ti.get('panel_index_load_seconds')} s ({ti.get('panel_index_file_bytes')} bytes); peak footprint cold "
          f"{tc.get('resources_phys_footprint_peak_bytes')} vs indexed {ti.get('resources_phys_footprint_peak_bytes')}")


def cold_vs_index(cell, d, library=None, extra=(), reads=None, tag=""):
    cold, idx = f"{d}/{cell}{tag}.cold", f"{d}/{cell}{tag}.index"
    p1 = run(E.cmd(cell, cold, library, list(extra), reads=reads))
    p2 = run(E.cmd(cell, idx, library, list(extra) + ["--panel-index", index_path(cell)], reads=reads))
    check(p1.returncode == p2.returncode and p1.returncode in (0, 1, 3), f"{cell}{tag}: cold and indexed exits agree ({p1.returncode}, {p2.returncode})")
    print(f"{cell}{tag}: wall cold {p1.wall:.1f} s, indexed {p2.wall:.1f} s")
    return cold, idx


def ex_index_c4(d):
    cold, idx = cold_vs_index("c4", d)
    compare("c4", cold, idx, "c4 (auto library)")
    return 0


def ex_index_cyp2d6(d):
    cold, idx = cold_vs_index("cyp2d6", d)
    compare("cyp2d6", cold, idx, "cyp2d6 (auto library)")
    return 0


def ex_index_acot(d):
    cold, idx = cold_vs_index("acot", d)
    compare("acot", cold, idx, "acot (auto library)")
    return 0


def ex_dosage_check_c4(d):
    cold, idx = cold_vs_index("c4", d, "declared", ["--dosage-check"], tag=".check")
    compare("c4", cold, idx, "c4 (declared library, --dosage-check)", dosage=True)
    return 0


def ex_two_samples_c4(d):
    """A second C4 sample from the same index; the index never changes."""
    before = sha(index_path("c4"))
    cold, idx = cold_vs_index("c4", d, reads=C4B, tag=".NA19650")
    compare("c4", cold, idx, "c4 second sample NA19650 (auto library)")
    a, b = E.prediction_key(idx), E.prediction_key(f"{OUT}/ex.index_c4/c4.index")
    print(f"recorded: the two samples' indexed haplotype pairs {'differ' if a[3] != b[3] else 'are equal'}")
    rec = dict(l.rstrip("\n").split("\t", 1) for l in open(f"{OUT}/ex.build_c4/index/c4.build.tsv") if "\t" in l)
    check(sha(index_path("c4")) == before == rec.get("file_sha256"), "c4: the index is byte-identical to its build after every use")
    return 0


def ex_refusals_c4(d):
    """Incompatible, missing and damaged indexes are refused with a named reason before the prefix is touched."""
    idx = index_path("c4")
    graph_paths = [l.split("\t")[1] for l in open(E.CELLS["c4"]["graph"]) if l.startswith("P\t")]
    other_ref = next(p for p in graph_paths if p != E.CELLS["c4"]["ref"])
    cur = f"{d}/manifest_curation.txt"   # A1: a name outside every refusal prefix (curation.txt was counted at prefix "curation")
    open(cur, "w").write(other_ref + "\n")
    data = open(idx, "rb").read()
    open(f"{d}/half_index.pvpl", "wb").write(data[: len(data) // 2])
    flipped = bytearray(data)
    flipped[len(flipped) - 17] ^= 0x5A
    open(f"{d}/flipped_index.pvpl", "wb").write(bytes(flipped))

    def refused(tag, a, needle):
        prefix = f"{d}/{tag}"
        sentinels = [".mosaic_status.tsv", ".sequence_locus.tsv", ".settings.tsv", ".sequence_haplotypes.fa"]
        for s in sentinels:
            open(prefix + s, "w").write("sentinel\n")
        p = run(a)
        kept = all(open(prefix + s).read() == "sentinel\n" for s in sentinels)
        owned = sorted(f for f in os.listdir(d) if f.startswith(tag + "."))
        check(p.returncode == 2 and needle in p.stderr and kept and len(owned) == len(sentinels),
              f"refused ({tag}): exit 2, '{needle}' named, the prefix untouched ({owned})")

    base = lambda tag, idx_file=idx, cell="c4": E.cmd(cell, f"{d}/{tag}", "declared", ["--panel-index", idx_file])
    refused("curation", base("curation") + ["--panel-curation", cur], "panel curation")
    a = base("kmer") + ["--kmer-size", "27"]
    refused("kmer", a, "kmer_size")
    a = base("reference")
    a[a.index("--reference-path") + 1] = other_ref
    refused("reference", a, "reference path")
    a = E.cmd("cyp2d6", f"{d}/graph", "declared", ["--panel-index", idx])
    refused("graph", a, "graph content")
    refused("missing", base("missing", f"{d}/no_such.pvpl"), "(io)")
    refused("truncated", base("truncated", f"{d}/half_index.pvpl"), "(truncated)")
    refused("corrupted", base("corrupted", f"{d}/flipped_index.pvpl"), "(checksum_mismatch)")
    p = run(locus_cmd("c4", ["--build-panel-index", f"{d}/x.pvpl", "--panel-index", idx]))
    check(p.returncode == 2 and "cannot be combined" in p.stderr and not os.path.exists(f"{d}/x.pvpl"), "build + load refused, nothing built")
    p = run(locus_cmd("c4", ["--build-panel-index", f"{d}/y.pvpl", "--reads1", "r1.fq"]))
    check(p.returncode == 2 and "is a sample argument" in p.stderr and not os.path.exists(f"{d}/y.pvpl"), "a build with reads refused, nothing built")
    return 0


if __name__ == "__main__":
    name, d = sys.argv[1], sys.argv[2]
    os.makedirs(d, exist_ok=True)
    globals()["ex_" + name](d)
    bad = [w for ok, w in CHECKS if not ok]
    print(f"{len(CHECKS) - len(bad)} of {len(CHECKS)} checks passed" + (f"; FAILED: {bad}" if bad else ""))
    sys.exit(1 if bad else 0)
