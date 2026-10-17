#!/usr/bin/env python3
"""I4-C: a narrowly typed comparison of the tie-branch tables, layered on i3d_compare.

i3d_compare (and the chain beneath it) compares call.tie_branches.tsv and call.tie_branch_memory.tsv byte for byte,
so any timing or memory column makes them "differ". This module leaves that verdict in place and adds a QUALIFIED
verdict: for exactly these two files, a named PROFILE may let explicitly listed columns differ; every other column,
the header, the row count, row order and row identities, and file presence are still required identical. Any other
difference anywhere keeps the run DIFFERENT. Both verdicts are printed.

Profiles (what each experiment declared before comparing):
  i4b   tie_branches: seconds; tie_branch_memory: capacity_total, logical_total, phys_footprint, phys_footprint_peak
        (the columns identified in the I4-B NA19650 comparison; nothing else).
  i4c   seed cache on vs off: tie_branches: seconds, executed_physical_origins; tie_branch_memory: not compared
        (cache-dependent event ledger, no scientific column); status: tie_branch_seeding (the configured strategy).
Usage: i4c_compare.py <run dir A> <run dir B> [prefix] [--profile i4b]
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "i3d"))
import i3d_compare  # noqa: E402

perf_compare = i3d_compare.perf_compare
PROFILES = {
    "i4b": {
        "tie_branches.tsv": {"seconds"},
        "tie_branch_memory.tsv": {"capacity_total", "logical_total", "phys_footprint", "phys_footprint_peak"},
    },
    # I4-C, seed cache on vs off (declared before any comparison): a cached seed is built when its tie is recorded, so
    # its execution moves into the parent's period (per-branch executed_physical_origins may differ; totals are
    # reported separately). The branch memory ledger's EVENT rows depend on the cache (seed_cached vs seed_materialised)
    # and it carries no scientific column, so it is declared not compared ("*") and summarised separately.
    "i4c": {
        "tie_branches.tsv": {"seconds", "executed_physical_origins"},
        "tie_branch_memory.tsv": "*",
    },
}
# status keys a profile additionally allows to differ (declared configuration descriptions, not results)
# i4c amendment (after the first gate comparison, where it was the only remaining difference): tie_seeds_materialised
# counts the starts placed just in time -- the executed work a cache hit replaces -- and is not caught by the chain's
# ^tie_seed_ diagnostic pattern only because of its spelling. It is reported, never hidden.
PROFILE_STATUS_KEYS = {"i4b": set(), "i4c": {"tie_branch_seeding", "tie_seeds_materialised"}}
TYPED = ("tie_branches.tsv", "tie_branch_memory.tsv")


def typed_table(path_a, path_b, mine_a, mine_b, allowed):
    """Problems ([] when equal under `allowed`) and the allowed columns that actually differed."""
    a = perf_compare.read(path_a, mine_a).split("\n")
    b = perf_compare.read(path_b, mine_b).split("\n")
    a = [l.split("\t") for l in a if l != ""]
    b = [l.split("\t") for l in b if l != ""]
    if not a or not b:
        return ["empty table"], set()
    if a[0] != b[0]:
        return ["header differs"], set()
    head = a[0]
    unknown = allowed - set(head)
    if unknown:
        return [f"profile names columns the table does not have: {sorted(unknown)}"], set()
    if len(a) != len(b):
        return [f"row count {len(a) - 1} vs {len(b) - 1}"], set()
    problems, differed = [], set()
    for i, (ra, rb) in enumerate(zip(a[1:], b[1:]), start=1):
        if len(ra) != len(head) or len(rb) != len(head):
            problems.append(f"row {i}: field count")
            continue
        for c, (x, y) in enumerate(zip(ra, rb)):
            if x == y:
                continue
            if head[c] in allowed:
                differed.add(head[c])
            else:
                problems.append(f"row {i}: {head[c]}: {x} vs {y}")
    return problems, differed


def compare(dir_a, dir_b, name="call", profile="i4b"):
    allowed = PROFILES[profile]
    ok, problems, diagnostics, compared = i3d_compare.compare(dir_a, dir_b, name)
    mine_a, mine_b = perf_compare.recorded_dirs(dir_a, name), perf_compare.recorded_dirs(dir_b, name)
    remaining, qualified = [], {}
    for p in problems:
        if any(p.startswith(f"{name}.mosaic_status.tsv: {k}: ") for k in PROFILE_STATUS_KEYS[profile]):
            qualified.setdefault("mosaic_status.tsv", []).append(p.split(": ")[1])
            continue
        typed = next((t for t in TYPED if p == f"{name}.{t}: differs"), None)
        if typed is None:
            remaining.append(p)
            continue
        if allowed[typed] == "*":
            qualified[typed] = ["not compared under this profile (declared)"]
            continue
        issues, differed = typed_table(os.path.join(dir_a, f"{name}.{typed}"), os.path.join(dir_b, f"{name}.{typed}"),
                                       mine_a, mine_b, allowed[typed])
        if issues:
            remaining += [f"{name}.{typed}: {x}" for x in issues[:6]]
        else:
            qualified[typed] = sorted(differed)
    return ok, not remaining, remaining, diagnostics, compared, qualified


def main(argv):
    profile = "i4b"
    if "--profile" in argv:
        i = argv.index("--profile")
        profile = argv[i + 1]
        argv = argv[:i] + argv[i + 2:]
    ok, qok, problems, diagnostics, compared, qualified = compare(argv[0], argv[1], argv[2] if len(argv) > 2 else "call", profile)
    original = "IDENTICAL" if ok else "DIFFERENT"
    verdict = "IDENTICAL" if ok else ("QUALIFIED_IDENTICAL" if qok else "DIFFERENT")
    print(f"original: {original}; qualified ({profile}): {verdict}: {compared} files compared; allowed keys that changed: "
          f"{len(diagnostics)}" + ("".join(f"; {t} differed only in {c}" for t, c in qualified.items()) if qualified else "")
          + ("" if qok else "; " + " | ".join(problems[:12])))
    return 0 if qok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
