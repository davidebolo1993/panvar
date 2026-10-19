#!/usr/bin/env python3
"""M1-A: the frozen computational development baseline (accepted in S2-L), as a compact manifest built ONLY from the
accepted records, the code and the binaries -- no genotyping.
  build_baseline.py <out dir>
Writes: baseline_manifest.tsv (code, binaries, envelope, records), baseline_cells.tsv (inputs and identities per cell),
baseline_settings.tsv (every effective setting value and its source, per cell, from the S2-L combined runs; the two
experimental keys' explicit defaults for the control arm documented separately).
"""
import hashlib
import os
import subprocess
import sys

WT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", ".."))
S2L = "/Users/davide.bolognini/panvar/tests/results/s2l_cross_locus_20261019"
CELLS = ["c4_5", "cyp_5", "acot_5", "c4_30", "cyp_30", "acot_30"]
RECORDS = [
    ("s2e_fragment_blocked_20261018", "blocked union bit-identical to bulk at 5x"),
    ("s2f_blocked_fine_aggregates_20261018", "blocked+fine bit-identical; ACOT 30x union completes"),
    ("s2g_blocked_integration_20261018", "ACOT 30x union reference (1 B caps) bit-identical, 108,939 candidates; caller 5x identical"),
    ("s2j_origin_admission_20261018", "round-1 reuse-only admission; ACOT 30x full caller completes"),
    ("s2k_post_refinement_20261019", "guard lifecycle; 6 GiB-envelope ACOT 30x validation"),
    ("s2l_cross_locus_20261019", "cross-locus control vs combined: science identical on every completed pair; ACOT 30x repeat identical"),
]


def sha(path, n=None):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()[:n] if n else h.hexdigest()


def kv(path):
    out = {}
    for line in open(path, errors="replace"):
        k, _, v = line.rstrip("\n").partition("\t")
        out[k] = v
    return out


def git(*args):
    return subprocess.run(["git", "-C", WT, *args], capture_output=True, text=True).stdout.strip()


def main(out):
    os.makedirs(out, exist_ok=True)
    rows = [("baseline", "S2-L accepted computational development baseline (not an accuracy, calibration or production-default claim)"),
            ("code_revision", git("rev-parse", "HEAD")),
            ("cpp_last_changed", git("log", "-1", "--format=%H", "--", "src", "include", "CMakeLists.txt")),
            ("source_tree", git("rev-parse", "HEAD^{tree}")),
            ("binary_panvar_sha256", sha(f"{WT}/build/panvar")),
            ("binary_libpanvarlib_sha256", sha(f"{WT}/build/libpanvarlib.a")),
            ("binary_libedlib_sha256", sha(f"{WT}/build/libedlib.a")),
            ("envelope_footprint_limit_bytes", "6442450944 (declared, every run; a sampled guard, not an allocation ceiling)"),
            ("threads", "4 (one process at a time)"),
            ("combined_settings", "P1 'opt' profile of the cell + union_fragment_blocks fine_core_4096_experimental + "
                                  "repair_round1_origin_admission reuse_only_experimental + early_incumbent experimental_v1 + "
                                  "output_schema sequence_first_v2 + footprint_limit_bytes 6442450944; time limit = unit timeout - 20 s"),
            ("control_settings_note", "S2-L control settings files OMIT union_fragment_blocks and repair_round1_origin_admission "
                                      "(54 vs 56 keys): the control ran at their code defaults, verified in "
                                      "include/panvar/genotype_command.hpp -- union_fragment_blocks \"off\", repair_round1_origin_admission "
                                      "\"existing\". The key sets are therefore NOT identical; the effective values differ only in those two."),
            ("default_union_fragment_blocks", "off"),
            ("default_repair_round1_origin_admission", "existing"),
            ("default_early_incumbent", "off (the baseline sets experimental_v1)"),
            ("default_footprint_limit_bytes", "0 (none; the baseline declares 6442450944)")]
    for name, what in RECORDS:
        rows.append((f"record:{name}", what))
    with open(f"{out}/baseline_manifest.tsv", "w") as o:
        o.write("key\tvalue\n" + "".join(f"{k}\t{v}\n" for k, v in rows))
    sys.path.insert(0, os.path.join(WT, "tests", "analysis", "p1"))
    import p1_run  # noqa: E402
    cells = p1_run.cells()
    with open(f"{out}/baseline_cells.tsv", "w") as o:
        o.write("cell\tgraph\treference\treads1\treads2\tlibrary\tcuration\tbase_profile\tbase_profile_sha256\tmodel_fingerprint\t"
                "sha256_gfa\tsha256_reads1\tsha256_reads2\tsha256_library\tprepared_locus\tsettings_source_run\n")
        for c in CELLS:
            cell = cells[c]
            st = kv(f"{S2L}/runs/combined.{c}/call.mosaic_status.tsv")
            o.write("\t".join([c, cell["graph"], cell["reference"], cell["reads1"], cell["reads2"], cell["library"], cell["curation"],
                               cell["base_profile"], sha(cell["base_profile"]), st.get("model_fingerprint", "NA"),
                               st.get("sha256_gfa", "NA"), st.get("sha256_reads1", "NA"), st.get("sha256_reads2", "NA"),
                               st.get("sha256_library", "NA"), "none", f"runs/combined.{c}"]) + "\n")
    with open(f"{out}/baseline_settings.tsv", "w") as o:
        o.write("cell\tsetting\tvalue\tsource_kind\n")
        for c in CELLS:
            for k, v in kv(f"{S2L}/runs/combined.{c}/call.settings.tsv").items():
                if k == "key":
                    continue
                value, _, source = v.partition("\t")
                kind = source.split(" ")[0] if source else "NA"
                o.write(f"{c}\t{k}\t{value}\t{kind}\n")
    print(open(f"{out}/baseline_manifest.tsv").read())
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
