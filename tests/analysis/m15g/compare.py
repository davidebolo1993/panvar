#!/usr/bin/env python3
"""M15-G: the C++ assessment agrees with the executable specification, byte for byte, on every fixture (and on any
extra prefixes given). Usage: compare.py <genotype_assessment_test> <scratch dir> [prefix ...]"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import assessment_spec as S  # noqa: E402

binary, scratch = sys.argv[1], sys.argv[2]
fixtures = S.write_fixtures(os.path.join(scratch, "fixtures"))
prefixes = [(n, p) for n, p in fixtures] + [(os.path.basename(p), p) for p in sys.argv[3:]]
bad = 0
for name, prefix in prefixes:
    want = S.render(prefix)
    got = subprocess.run([binary, "show", prefix], capture_output=True, text=True).stdout
    if got != want:
        bad += 1
        print(f"DIFF {name}")
        for a, b in zip(want.splitlines(), got.splitlines()):
            if a != b:
                print(f"  spec: {a}\n  c++ : {b}")
                break
        else:
            print(f"  line counts: spec {len(want.splitlines())}, c++ {len(got.splitlines())}")
print(f"assessment C++ vs specification: {len(prefixes) - bad} of {len(prefixes)} identical")
fx_bad = 0
for name, prefix in fixtures:
    if name in S.EXPECT:
        L, Bk = S.assess(S.read(prefix))
        fx_bad += not S.EXPECT[name](L, Bk)
print(f"counterexample expectations: {len(S.EXPECT) - fx_bad} of {len(S.EXPECT)} hold")
sys.exit(1 if bad or fx_bad else 0)
