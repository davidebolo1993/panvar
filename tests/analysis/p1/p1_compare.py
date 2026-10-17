#!/usr/bin/env python3
"""P1: compare two completed runs of one cell and report EVERY difference, strict first.

Two verdicts per pair, both printed:
  * strict: tests/analysis/perf/perf_compare.py exactly as committed (its own DIAGNOSTIC class only: time, memory,
    cache and executed-work keys), with the COMPLETE problem list (perf_compare's CLI prints only twelve);
  * qualified: tests/analysis/i4c/i4c_compare.py --profile i4b (the reviewed I-series chain, which additionally allows the
    backend/provenance keys each experiment declared, and the tie tables' timing and memory columns).
A pair PASSES only if the qualified verdict passes AND no strict problem touches a scientific output (SCIENCE below) or
the repair trajectory. Every strict problem that remains (settings rows, status keys) is listed with its class:
key present on one side only ("added"/"removed") or a value change ("changed"). Nothing is hidden.

  p1_compare.py KIND DIR_A DIR_B [PREFIX]      -> one TSV row on stdout; exit 0 iff PASS
"""
import importlib.util
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.join(HERE, "..")
# outputs that carry the science: genotypes, phase/support and filters (block and locus tables), scores and components,
# candidates, ranking and ties (equivalence sets), FASTA, dosage, library sensitivity, and the trajectory
SCIENCE = (".block_genotypes.tsv", ".genotypes.tsv", ".locus_genotype.tsv", ".mosaic_scores.tsv", ".mosaic_candidates.tsv",
           ".equivalence_set.tsv", ".mosaic_haplotypes.fa", ".dosage", ".library_sensitivity.tsv", ".repair_rounds.tsv",
           ".original")
# Amendment 1 (declared after the first comparison, before any further execution): (a) a status file -- the
# .mosaic_status.tsv / .dosage_status.tsv of a run, including their .original copies -- is NOT a scientific output by its
# file name; the first version's suffix rule matched ".original" and ".dosage" inside those names and so classified
# executed-work counters of the status copy as science. Status differences are judged key by key instead: (b) the keys
# below must be strictly identical (this is stricter than the qualified chain alone), and every other status key
# difference is listed and must pass the qualified chain.
SCIENCE_STATUS = re.compile(r"^(status|call_class|no_call_reason|phase_support|score_gap_to_runner_up|"
                            r"best_other_genotype_score|alleles_scored|every_allele_scored|blocks_genotype\w*|"
                            r"quality_\w+|dosage_audit\w*|dosage_status_\w+|dosage_repair\w*|dosage_allocation\w*|"
                            r"dosage_policy|dosage_census_sensitive_blocks|dosage_molecule_span_bp|\w*genotype\w*)$")


def strict_compare(a, b, name):
    spec = importlib.util.spec_from_file_location("perf_compare_strict", os.path.join(ANALYSIS, "perf", "perf_compare.py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)   # a fresh module: the I-series chain extends the class of a shared one
    return mod.compare(a, b, name)


def qualified_compare(a, b, name):
    r = subprocess.run([sys.executable, os.path.join(ANALYSIS, "i4c", "i4c_compare.py"), a, b, name, "--profile", "i4b"],
                       capture_output=True, text=True)
    return r.returncode == 0, r.stdout.strip().replace("\t", " ")


def classify(problem):
    if ".mosaic_status.tsv: " in problem or ".dosage_status.tsv: " in problem:
        _, key, rest = problem.split(": ", 2)
        left, right = rest.split(" vs ", 1)
        kind = "added" if left == "None" else "removed" if right == "None" else "changed"
        return f"status:{key}:{kind}"
    return problem


def main(argv):
    kind, a, b = argv[0], argv[1], argv[2]
    name = argv[3] if len(argv) > 3 else "call"
    s_ok, problems, diagnostics, compared = strict_compare(a, b, name)
    q_ok, q_line = qualified_compare(a, b, name)
    def is_status(p):
        return p.split(":")[0].endswith("_status.tsv")
    science = [p for p in problems if (not is_status(p) and any(x in p.split(":")[0] for x in SCIENCE))
               or "file sets differ" in p]
    status_science = [p for p in problems if is_status(p) and len(p.split(": ")) > 2 and
                      SCIENCE_STATUS.match(p.split(": ")[1])]
    ok = q_ok and not science and not status_science
    row = [kind, os.path.basename(a.rstrip("/")) + "/" + name, os.path.basename(b.rstrip("/")) + "/" + name,
           "PASS" if ok else "FAIL", "IDENTICAL" if s_ok else "DIFFERENT", str(compared), str(len(diagnostics)),
           str(len(problems)), " | ".join(classify(p) for p in problems) or "none",
           " | ".join(science + status_science) or "none", q_line]
    print("\t".join(row))
    return 0 if ok else 1


HEADER = ["kind", "a", "b", "verdict", "strict", "files_compared", "diagnostic_keys_changed", "strict_problems",
          "strict_problem_list", "science_problems", "qualified_i4c_i4b"]

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
