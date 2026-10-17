#!/usr/bin/env python3
"""P1.1: positive and negative controls for compose_profile.py (no panvar run; temporary files only).

Required refusals: scientific/model overrides (kmer_size, syncmer_s, min_variant_bp, the copying prior
switch_probability, the emission model mismatch_band), work-cap overrides (max_state_classes, repair_max_neighbours),
an override of a key the base does not set, a missing/empty reason or value, a malformed override, a duplicated
override, an invalid mode, a duplicate key in the base, a key set by both base and overlay, and an audited profile whose
override section holds a scientific key or lacks a reason. Required acceptances: the three modes with the two allowed
guards, each output passing the audit and carrying exactly the declared sections.
  test_compose_profile.py        -> prints one line per case; exit 0 iff every case behaves as required
"""
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compose_profile as cp  # noqa: E402

BASE = """profile\tfixture-base
time_limit_seconds\t900
footprint_limit_bytes\t6442450944
max_state_classes\t600000000
repair_max_neighbours\t958
kmer_size\t31
syncmer_s\t11
min_variant_bp\t50
switch_probability\t0.01
mismatch_band\t6
"""
OK_OVERRIDES = ["time_limit_seconds=280:declared unit timeout minus 20 s",
                "footprint_limit_bytes=4294967296:4 GiB local test envelope"]
failures = []


def case(name, fn, refused):
    try:
        fn()
        got = False
    except ValueError as e:
        got, msg = True, str(e)
    ok = got == refused
    print(f"{'ok' if ok else 'FAILED'}: {name}: {'refused' if got else 'accepted'}" + (f" ({msg[:90]})" if got else ""))
    if not ok:
        failures.append(name)


def main():
    d = tempfile.mkdtemp(prefix="p11_compose_")
    base = os.path.join(d, "base.tsv")
    open(base, "w").write(BASE)
    out = os.path.join(d, "out.tsv")

    def compose(mode="full", overrides=OK_OVERRIDES, base_path=base):
        return lambda: cp.compose(out, "fixture", base_path, mode, overrides)

    for key, value in [("kmer_size", "23"), ("syncmer_s", "9"), ("min_variant_bp", "40"),
                       ("switch_probability", "0.05"), ("mismatch_band", "8"), ("max_state_classes", "900000000"),
                       ("repair_max_neighbours", "2000")]:
        case(f"scientific/cap override {key} refused", compose(overrides=[f"{key}={value}:a reason"]), True)
    nobase = os.path.join(d, "nobase.tsv")
    open(nobase, "w").write("profile\tx\nmax_state_classes\t1\n")
    case("override of a key the base does not set refused", compose(base_path=nobase), True)
    case("missing reason refused", compose(overrides=["time_limit_seconds=280"]), True)
    case("empty reason refused", compose(overrides=["time_limit_seconds=280:  "]), True)
    case("empty value refused", compose(overrides=["time_limit_seconds=:reason"]), True)
    case("malformed override refused", compose(overrides=["time_limit_seconds:280"]), True)
    case("duplicated override refused", compose(overrides=OK_OVERRIDES + [OK_OVERRIDES[0]]), True)
    for mode in ["", "Full", "fulll", "ref", "reference "]:
        case(f"invalid mode {mode!r} refused", compose(mode=mode), True)
    dup = os.path.join(d, "dup.tsv")
    open(dup, "w").write(BASE + "kmer_size\t31\n")
    case("duplicate key in the base refused", compose(base_path=dup), True)
    clash = os.path.join(d, "clash.tsv")
    open(clash, "w").write(BASE + "emission_cache_bytes\t1\n")
    case("key in both base and overlay refused (full)", compose(base_path=clash), True)
    case("key in both base and overlay refused (reference)", compose(mode="reference", base_path=clash), True)
    case("the same base composes with no overlay (none)", compose(mode="none", base_path=clash), False)

    for mode in cp.MODES:
        case(f"mode {mode} with the two allowed guards accepted", compose(mode=mode), False)
        rows = dict(cp.read_profile(out))
        text = open(out).read()
        if rows["time_limit_seconds"] != "280" or rows["kmer_size"] != "31" or rows["switch_probability"] != "0.01":
            failures.append(f"{mode}: values")
        overlay_keys = set(rows) & set(k for ks in cp.OVERLAY_SECTIONS.values() for k in ks)
        expect = (set(k for ks in cp.OVERLAY_SECTIONS.values() for k in ks) if mode == "full"
                  else set(cp.REFERENCE_ROUTE) if mode == "reference" else set())
        if overlay_keys != expect or cp.check_composed(out) != ["time_limit_seconds", "footprint_limit_bytes"]:
            failures.append(f"{mode}: sections")
        print(f"{'ok' if f'{mode}: sections' not in failures else 'FAILED'}: mode {mode} overlay keys = expected "
              f"({len(expect)}), audit passes ({'## operational overrides' in text})")
    case("no override at all accepted", compose(overrides=[]), False)

    bad = os.path.join(d, "bad.tsv")
    open(bad, "w").write("profile\tx\n## operational overrides (declared; resource guards only)\n# kmer_size: smuggled\n"
                         "kmer_size\t23\n")
    case("audit refuses a scientific key in an override section", lambda: cp.check_composed(bad), True)
    noreason = os.path.join(d, "noreason.tsv")
    open(noreason, "w").write("profile\tx\n## operational overrides (declared; resource guards only)\n"
                              "time_limit_seconds\t10\n")
    case("audit refuses an override without its reason", lambda: cp.check_composed(noreason), True)
    print(f"{'ALL PASSED' if not failures else 'FAILURES: ' + ', '.join(failures)}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
