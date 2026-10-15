#!/usr/bin/env python3
"""Is a performance-only setting scientifically invisible? Compares two runs of one cell (no truth is read).

Required identical (after replacing each run's directory with one placeholder):
  * byte-identical: block/locus genotype tables, genotypes, equivalence sets, candidate scores and candidates,
    FASTA, dosage tables, library sensitivity (each also as its .original copy where written), settings apart from the
    declared budget keys and the profile path/identity;
  * status (key/value): every key, except time, memory, cache, executed-work and profile/budget keys (DIAGNOSTIC);
  * repair trajectory (.repair_rounds.tsv): the scientific columns (SCIENTIFIC_ROUND_COLUMNS);
Not compared, listed: the repair memory ledger (.repair_memory.tsv) and the runner's own log/record.
Both runs must have completed. Differences in DIAGNOSTIC keys are reported (they are the measured effect).

Usage: perf_compare.py <run dir A> <run dir B> [<prefix name, default call>]  -> prints a verdict; exit 0 if identical
"""
import os
import re
import sys

DIAGNOSTIC = re.compile(
    r"(^seconds|_seconds$|seconds_|wall|rss|footprint|phys_|_gib|timing|cache|cached|scratch|_hits$|_misses$|"
    r"^profile$|^profile_name$|^profile_sha256$|^physical_origins$|^statistical_state_classes$|^merged_state_classes$|"
    r"^piece_hits$|^reused_|^scorer_|^threads_peak|^resources_|^repair_cumulative_|^repair_margin_|_peak_bytes$|"
    r"^max_union_|^union_|^emission_|^start_sort|^repair_origin_cache|^tie_seed_|^tie_branch_physical|"
    r"^representation_guard_physical_origins$|^marker_lookup|^seeding_|"
    # executed-work counters: what was computed rather than replayed from a cache (they are the measured effect)
    r"^candidate_bases_scanned$|^dp_evaluations$|^full_read_verifications$|^prefilter_rejections$|^proposed_starts$|"
    r"^unique_starts$|^valid_fr_states$|^threads_start_sort_)")
SCIENTIFIC_ROUND_COLUMNS = ["round", "start_total", "moved", "repaired_block", "local_pair_rank", "assignment", "allele_x",
                            "allele_y", "hap1_alleles", "hap2_alleles", "fragment_term", "exposure_term", "prior", "total",
                            "delta_total", "status", "unique_best", "separation", "neighbours_generated",
                            "neighbours_scored", "duplicate", "suppressed_visited", "self_skipped",
                            "rejected_infinite_prior", "unique_haplotypes", "scored_pairs", "order_invariant",
                            "swap_invariant"]
NOT_COMPARED = (".repair_memory.tsv", ".call.log", ".run.tsv", "profile.tsv")


def recorded_dirs(d, name):
    """The run's own directory as given and as the caller recorded it (the status `profile` path, whose file lives in
    the run directory, and the locus table's haplotype ids), longest first."""
    dirs = {os.path.abspath(d)}
    st = os.path.join(d, name + ".mosaic_status.tsv")
    if os.path.exists(st):
        for line in open(st).read().split("\n"):
            p = line.split("\t")
            if p[0] == "profile" and len(p) > 1:
                dirs.add(os.path.dirname(p[1]))
    lg = os.path.join(d, name + ".locus_genotype.tsv")
    if os.path.exists(lg):
        lines = open(lg).read().split("\n")
        head = lines[0].split("\t")
        if "hap1_id" in head and len(lines) > 1 and "/" in lines[1].split("\t")[head.index("hap1_id")]:
            dirs.add(os.path.dirname(lines[1].split("\t")[head.index("hap1_id")]))
    return sorted((x for x in dirs if len(x) > 1), key=len, reverse=True)


def read(path, mine, placeholder="<RUN>"):
    text = open(path, "rb").read().decode("utf-8", "surrogateescape")
    for m in mine:
        text = text.replace(m, placeholder)
    return text


def kv(text):
    out = {}
    for line in text.split("\n"):
        p = line.split("\t")
        if len(p) >= 2:
            out.setdefault(p[0], "\t".join(p[1:]))
    return out


def compare(dir_a, dir_b, name="call"):
    fa = sorted(f for f in os.listdir(dir_a) if f.startswith(name + ".") or f == "profile.tsv")
    fb = sorted(f for f in os.listdir(dir_b) if f.startswith(name + ".") or f == "profile.tsv")
    problems, diagnostics, compared = [], set(), 0
    mine_a, mine_b = recorded_dirs(dir_a, name), recorded_dirs(dir_b, name)
    if fa != fb:
        problems.append(f"file sets differ: {sorted(set(fa) ^ set(fb))}")
    for f in sorted(set(fa) & set(fb)):
        if f.endswith(NOT_COMPARED):
            continue
        a, b = read(os.path.join(dir_a, f), mine_a), read(os.path.join(dir_b, f), mine_b)
        compared += 1
        if f.endswith(".mosaic_status.tsv") or f.endswith(".dosage_status.tsv"):
            ka, kb = kv(a), kv(b)
            for k in sorted(set(ka) | set(kb)):
                if ka.get(k) != kb.get(k):
                    (diagnostics.add(k) if DIAGNOSTIC.search(k) else problems.append(f"{f}: {k}: {ka.get(k)} vs {kb.get(k)}"))
        elif f.endswith(".settings.tsv"):
            la = [l for l in a.split("\n") if not DIAGNOSTIC.search(l.split("\t")[0])]
            lb = [l for l in b.split("\n") if not DIAGNOSTIC.search(l.split("\t")[0])]
            strip = lambda ls: [re.sub(r"profile \S+", "profile <P>", l) for l in ls]  # noqa: E731
            if strip(la) != strip(lb):
                problems.append(f"{f}: non-budget settings differ")
        elif f.endswith(".repair_rounds.tsv"):
            ra, rb = [l.split("\t") for l in a.split("\n") if l], [l.split("\t") for l in b.split("\n") if l]
            ha, hb = ra[0], rb[0]
            cols = [c for c in SCIENTIFIC_ROUND_COLUMNS if c in ha]
            pa = [[r[ha.index(c)] for c in cols] for r in ra[1:]]
            pb = [[r[hb.index(c)] for c in cols] for r in rb[1:]]
            if pa != pb:
                problems.append(f"{f}: repair trajectory differs")
        elif a != b:
            problems.append(f"{f}: differs")
    return not problems, problems, sorted(diagnostics), compared


if __name__ == "__main__":
    ok, problems, diagnostics, compared = compare(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "call")
    print(("IDENTICAL" if ok else "DIFFERENT") + f": {compared} files compared; diagnostic keys that changed: "
          f"{len(diagnostics)}" + ("" if ok else "; " + " | ".join(problems[:12])))
    sys.exit(0 if ok else 1)
