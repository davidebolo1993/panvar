#!/usr/bin/env python3
"""CYP2D6/D7 pilot: are two runs of one cell scientifically identical? (call outputs only; no truth)

Every output file of the cell (all <job>.* files except the runner's own .call.log and .run.tsv, whose exit codes are
compared separately) must exist in both runs and be byte-identical after removing ONLY the declared differences:
  * the output path: each run's directory name is replaced by one placeholder;
  * the profile path, and the profile's own identity (profile, profile_name, profile_sha256 rows of status files and
    the profile path in settings.tsv's source column);
  * the declared settings (--declared, e.g. max_state_classes): their rows in status and settings files are removed,
    after checking that each run's value equals what its own profile declares (or the built-in default it reports);
  * time and memory measurements: keys or columns named by TIMING_MEMORY below, and nothing else.
Everything else (genotype and equivalence set, scores and ranking, block and locus tables, filters, repair trajectory
and its scientific counters, FASTA) must be byte-identical. A file whose normalised content differs fails the gate.

Usage: identity_gate.py <out tsv> <declared keys, comma list> <run dir A> <profile A> <run dir B> <profile B> <job>...
"""
import os
import re
import sys

TIMING_MEMORY = re.compile(
    r"^(seconds(_.+)?|.+_seconds|.+_seconds_all|round_wall_seconds|peak_rss_gib|rss_current|rss_peak|rss_after_relief|"
    r"rss_minus_capacity|phys_footprint|phys_footprint_peak|resources_rss_peak_bytes|resources_phys_footprint_peak_bytes|"
    r"resources_footprint_at_repair_start_bytes)$")
PROFILE_IDENTITY = {"profile", "profile_name", "profile_sha256"}


def declared_value(profile, key):
    for line in open(profile):
        p = line.rstrip("\n").split("\t")
        if p[0] == key:
            return p[1]
    return None


def normalise(path, run_dir, profile, declared):
    text = open(path, "rb").read().decode("utf-8", errors="surrogateescape")
    text = text.replace(os.path.abspath(profile), "<PROFILE>").replace(run_dir, "<RUN>")
    lines = text.split("\n")
    name = os.path.basename(path)
    head = lines[0].split("\t") if lines else []
    problems = []
    if name.endswith(".settings.tsv") or head[:2] == ["key", "value"]:
        kept = []
        for line in lines:
            key = line.split("\t")[0]
            if key in declared:
                want = declared_value(profile, key)
                got = line.split("\t")[1] if "\t" in line else None
                if want is not None and got != want:
                    problems.append(f"{name}: {key} = {got}, profile declares {want}")
                continue
            if TIMING_MEMORY.match(key) or key in PROFILE_IDENTITY:
                continue
            kept.append(line)
        return "\n".join(kept), problems
    if len(head) > 1 and not name.endswith((".fa", ".fasta")):
        keep = [i for i, h in enumerate(head) if not TIMING_MEMORY.match(h)]
        return "\n".join("\t".join(p[i] for i in keep if i < len(p)) for p in (l.split("\t") for l in lines)), problems
    return text, problems


def files_of(run_dir, job):
    skip = {job + ".call.log", job + ".run.tsv"}
    return sorted(f for f in os.listdir(run_dir) if f.startswith(job + ".") and f not in skip)


def exit_of(run_dir, job):
    lines = open(os.path.join(run_dir, job + ".run.tsv")).read().splitlines()
    return dict(zip(lines[0].split("\t"), lines[1].split("\t")))["exit"]


def compare(job, dir_a, profile_a, dir_b, profile_b, declared):
    """Returns (ok, detail)."""
    fa, fb = files_of(dir_a, job), files_of(dir_b, job)
    if fa != fb:
        return False, f"file sets differ: {sorted(set(fa) ^ set(fb))}"
    ea, eb = exit_of(dir_a, job), exit_of(dir_b, job)
    if ea != eb:
        return False, f"exit codes differ: {ea} vs {eb}"
    identical, normalised_equal, differing, problems = 0, 0, [], []
    for f in fa:
        pa, pb = os.path.join(dir_a, f), os.path.join(dir_b, f)
        if open(pa, "rb").read() == open(pb, "rb").read():
            identical += 1
            continue
        na, qa = normalise(pa, os.path.basename(os.path.normpath(dir_a)), profile_a, declared)
        nb, qb = normalise(pb, os.path.basename(os.path.normpath(dir_b)), profile_b, declared)
        problems += qa + qb
        if na == nb:
            normalised_equal += 1
        else:
            differing.append(f)
    ok = not differing and not problems
    return ok, (f"exit {ea}; {len(fa)} files: {identical} byte-identical, {normalised_equal} identical after removing "
                f"declared differences" + (f"; DIFFER: {differing}" if differing else "")
                + (f"; declared-value problems: {problems}" if problems else ""))


def main():
    out, declared, dir_a, profile_a, dir_b, profile_b = sys.argv[1:7]
    declared = {k for k in declared.split(",") if k}
    rows = []
    for job in sys.argv[7:]:
        ok, detail = compare(job, dir_a, profile_a, dir_b, profile_b, declared)
        rows.append((job, "PASS" if ok else "FAIL", detail))
        print("\t".join(rows[-1]))
    with open(out, "w") as f:
        f.write(f"# identity of {dir_a} ({profile_a}) and {dir_b} ({profile_b}); declared: {sorted(declared)}\n"
                "job\tstatus\tdetail\n" + "".join("\t".join(r) + "\n" for r in rows))
    failed = sum(r[1] != "PASS" for r in rows)
    print(f"{len(rows)} cells, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
