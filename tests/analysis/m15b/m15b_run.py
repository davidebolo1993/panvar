#!/usr/bin/env python3
"""M15-B runner: budget 7,200 s including preparation, tests, failures and evaluation; 4 threads; one process at a time; 6 GiB.
Order: build (+ frozen-binary assertion) -> full enabled ctest -> prepare (manifest, lineage, reads) -> 18 calls -> freeze
(prediction hashes) -> 18 evaluations + truth projections -> analyse. A failed regression stops reserved execution (prepare
needs ctest.full). Calls are admitted only if spent + 300 + the evaluation reserve fits. Usage: m15b_run.py manifest | next | <unit>."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "m11e"))
sys.path.insert(0, os.path.join(HERE, "..", "m10a"))
import m15b  # noqa: E402
import runner  # noqa: E402
import compose_m10a  # noqa: E402

WT, OUT = m15b.WT, m15b.OUT
BUDGET, PRE = 7200, 90.0   # disclosed: inventory/curation checks, archived-hash lookups, script writing (no caller run)
RC_PANVAR = "47042a5cbae5b6418d4c02f147c27772ec1a85791685dcf7c8f11bd8fa2af568"
PERF = "/Users/davide.bolognini/panvar/tests/results/genotype_perf_bounded_20261015/inputs/profiles"
WATCHDOG, GUARD = 300, 280
DONORS = [d for _, d in m15b.families()]
UNITS = [("build", 300), ("ctest.full", 2400), ("prepare", 600)] + [(f"call.{d}", WATCHDOG) for d in DONORS] + [("freeze", 60)]
for d in DONORS:
    UNITS += [(f"eval.{d}", 120), (f"truth.{d}", 60)]
UNITS += [("analyse", 300)]
RESERVE = 60 + 18 * (120 + 60) + 300


def needs(unit):
    if unit == "ctest.full":
        return ["build"]
    if unit == "prepare":
        return ["build", "ctest.full"]
    if unit.startswith("call."):
        return ["prepare"]
    if unit in ("freeze", "analyse") or unit.startswith(("eval.", "truth.")):
        return ["freeze"] if unit != "freeze" else ["prepare"]
    return []


def command(unit, d):
    py = [sys.executable, os.path.join(HERE, "m15b.py")]
    if unit == "build":
        check = (f"import hashlib,sys; h=hashlib.sha256(open('{WT}/build/panvar','rb').read()).hexdigest(); "
                 f"print('panvar', h); sys.exit(0 if h=='{RC_PANVAR}' else 1)")
        return ["sh", "-c", f"cmake --build {WT}/build -j 4 && {sys.executable} -c \"{check}\""]
    if unit == "ctest.full":
        return ["ctest", "--test-dir", f"{WT}/build", "-j", "1", "--output-on-failure", "--timeout", "900"]
    if unit in ("prepare", "freeze", "analyse"):
        return py + [unit, d]
    kind, donor = unit.split(".", 1)
    if kind == "truth":
        return py + ["truth", donor, d]
    if kind == "eval":
        return [f"{WT}/build/genotype_probe", "evaluate", "--scheme", "loo", "--donor", donor, "--reference", m15b.ACOT_REF,
                "--run-prefix", f"{OUT}/runs/call.{donor}/call", "--model-gfa", m15b.ACOT_GFA, "--truth-gfa", m15b.ACOT_GFA,
                "--out", f"{d}/eval", "--panel-curation", m15b.CURATION]
    # call: the named experimental profile (M10-A arm D, acot_5 base), the declared time guard, block diagnostics (diagnostic)
    prof = compose_m10a.compose(f"{d}/profile.tsv", f"m15b-rc-{donor}", "D", f"{PERF}/acot_5.tsv", "none",
                                [f"time_limit_seconds={GUARD}:declared watchdog {WATCHDOG} s minus 20 s (the caller records TIME_LIMIT first)"])
    diag = f"{d}/profile.block_diagnostics.tsv"
    open(diag, "w").write("profile\tm15b-block-diagnostics\nblock_diagnostics\tdiagnostic_v1\n")
    r = m15b.reads_dir(donor)
    return ["nice", f"{WT}/build/panvar", "genotype", "--preset", "short-read", "--profile", prof, "--threads", "4", "--graph", m15b.ACOT_GFA,
            "--reference-path", m15b.ACOT_REF, "--panel-curation", m15b.CURATION, "--reads1", f"{r}/reads_1.fq", "--reads2", f"{r}/reads_2.fq",
            "--library", f"{r}/library.txt", "--dosage-mode", "audit", "--out", f"{d}/call", "--quiet", "--profile", diag]


class ReservingRunner(runner.Runner):
    def precondition(self, unit):
        why = super().precondition(unit)
        if why or not unit.startswith("call."):
            return why
        spent, timeout = self.spent(), dict(self.units)[unit]
        if spent + timeout + RESERVE > self.budget:
            return f"budget: {spent:.1f} + {timeout} + evaluation reserve {RESERVE} > {self.budget} s"
        return None


R = ReservingRunner(OUT, BUDGET, PRE, UNITS, command, needs, revision=lambda: runner.git_revision(WT), require_clean=runner.git_clean(WT))

if __name__ == "__main__":
    if sys.argv[1] == "manifest":
        os.makedirs(f"{OUT}/runs", exist_ok=True)
        lines = [f"# executable\tpanvar\t{m15b.sha(m15b.PANVAR)}\t(release candidate {RC_PANVAR})", f"# executable\tgenotype_probe\t{m15b.sha(m15b.PROBE)}",
                 f"# graph\t{m15b.ACOT_GFA}\t{m15b.sha(m15b.ACOT_GFA)}", f"# curation\t{m15b.CURATION}\t{m15b.sha(m15b.CURATION)}",
                 f"# profiles\tbase {PERF}/acot_5.tsv {m15b.sha(PERF + '/acot_5.tsv')}\tperformance {m15b.sha(WT + '/profiles/genotype_dev_performance_v1.tsv')}"
                 f"\tm10a {m15b.sha(WT + '/profiles/genotype_dev_m10a_v1.tsv')}"]
        lines += [f"# family\t{f}\t{d}" for f, d in m15b.families()]
        lines += ["unit\tdeclared_timeout_seconds\tneeds"] + [f"{u}\t{t}\t{','.join(needs(u)) or '-'}" for u, t in UNITS]
        text = "\n".join(lines) + "\n"
        p = f"{OUT}/MANIFEST.tsv"
        assert not os.path.exists(p) or open(p).read() == text
        open(p, "w").write(text)
        print(text[:3000])
    elif sys.argv[1] == "next":
        assert os.path.exists(f"{OUT}/MANIFEST.tsv"), "manifest first"
        sys.exit(R.next())
    else:
        sys.exit(R.run_unit(sys.argv[1]))
