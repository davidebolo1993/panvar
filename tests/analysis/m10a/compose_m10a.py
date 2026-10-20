#!/usr/bin/env python3
"""M10-A: compose ONE profile per run from the three layers of genotype-dev-m10a-v1, so a run records a single profile
path and SHA-256 that contains every setting and lists every input with its SHA-256.

  layer 3   the cell's base profile: RESOURCE requirements only (RESOURCE_KEYS; anything else is refused)
  layer 1a  profiles/genotype_dev_performance_v1.tsv through tests/analysis/p1/compose_profile.py (full | reference)
  layer 1b  profiles/genotype_dev_m10a_v1.tsv: the sections the arm selects (ARMS); keys outside SECTIONS are refused
  layer 2   the model's compound_interval_declarations_v1 artifact (only with the annotation section)
Operational overrides (time / footprint guards only) follow compose_profile's rules. No key may come from two layers.

Arms (the validation plan's comparisons, DESIGN.md section 3):
  REF  reference backends (compose_profile 'reference'), no route, no JP1.1        -- backend reference
  O    full backends, no route, no JP1.1                                           -- P1 'opt' route
  D0   full backends + retained route, no JP1.1                                     -- the configuration without JP1.1
  D    full backends + retained route + JP1.1  = genotype-dev-m10a-v1               -- the named configuration
Every arm carries the reporting contract and annotation sections, so the annotation itself is compared between arms.

  compose_m10a.py OUT NAME ARM BASE_PROFILE DECLARATIONS|none [KEY=VALUE:REASON ...]
"""
import hashlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "p1"))
import compose_profile  # noqa: E402

M10A = os.path.abspath(os.path.join(HERE, "..", "..", "..", "profiles", "genotype_dev_m10a_v1.tsv"))
RESOURCE_KEYS = {"time_limit_seconds", "footprint_limit_bytes", "max_state_classes", "repair_max_neighbours"}
SECTIONS = {
    "retained route": ["early_incumbent", "union_fragment_blocks", "repair_round1_origin_admission"],
    "experimental scientific switch": ["repair_junction_proposals"],
    "reporting contract and additive annotation": ["output_schema", "compound_annotation"],
    "stated explicitly, equal to the default": ["repair_joint_search"],
}
COMMON = ("reporting contract and additive annotation", "stated explicitly, equal to the default")
ARMS = {"REF": ("reference", COMMON), "O": ("full", COMMON), "D0": ("full", COMMON + ("retained route",)),
        "D": ("full", COMMON + ("retained route", "experimental scientific switch"))}


def sha256(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def m10a_rows():
    rows = dict(compose_profile.read_profile(M10A))
    if rows.pop("profile", None) != "genotype-dev-m10a-v1":
        compose_profile.refuse("layer 1b must name itself genotype-dev-m10a-v1")
    allowed = {k for keys in SECTIONS.values() for k in keys}
    if set(rows) != allowed:
        compose_profile.refuse(f"layer 1b keys differ from its declared sections: {sorted(set(rows) ^ allowed)}")
    return rows


def compose(out, name, arm, base, declarations, overrides):
    if arm not in ARMS:
        compose_profile.refuse(f"arm must be one of {sorted(ARMS)}: {arm!r}")
    extra = {k for k, _ in compose_profile.read_profile(base)} - RESOURCE_KEYS - {"profile"}
    if extra:
        compose_profile.refuse(f"{base}: layer 3 holds resource requirements only; refused keys {sorted(extra)}")
    perf_mode, sections = ARMS[arm]
    tmp = out + ".layers13"
    compose_profile.compose(tmp, name, base, perf_mode, overrides)
    text = open(tmp).read()
    os.remove(tmp)
    present = {line.split("\t", 1)[0] for line in text.split("\n") if "\t" in line and not line.startswith("#")}
    rows = m10a_rows()
    lines = [f"# M10-A layer 1b: {M10A} sha256 {sha256(M10A)}; arm {arm}: {', '.join(sections)}"]
    for section, keys in SECTIONS.items():
        if section not in sections:
            continue
        lines.append(f"## m10a: {section}")
        for k in keys:
            if k in present:
                compose_profile.refuse(f"{k} is set by two layers")
            lines.append(f"{k}\t{rows[k]}")
    if declarations != "none":
        lines += ["## layer 2: model-bound declarations", f"# {declarations} sha256 {sha256(declarations)}",
                  f"compound_interval_declarations\t{declarations}"]
    open(out, "w").write(text + "\n".join(lines) + "\n")
    compose_profile.read_profile(out)   # no key twice in the result
    compose_profile.check_composed(out)
    return out


if __name__ == "__main__":
    compose(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6:])
