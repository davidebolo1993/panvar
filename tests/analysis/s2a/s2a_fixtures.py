#!/usr/bin/env python3
"""S2-A fixture gates: early_incumbent experimental_v1 on the S1-A fixtures (inputs reused read-only from the S1-C
results tree, commands replayed from its logs with this library, sequence_first_v2, dosage by the v2 default).

Cases: F1 (sequence-equivalent labellings), F5 (a genuine three-way tie of distinct sequences), F7full (completed
control) and F7 (round cap 1: a provisional run). Per case:
  off  -- sequence_first_v2 without the early incumbent;
  on   -- the same plus `early_incumbent experimental_v1`: a completed run must equal off in every output except the
          declared status/settings rows (completed-run mathematics unchanged); an incomplete run is compared and its
          differences recorded (the provisional incumbent is the best fully evaluated candidate, not the search's);
  stop -- on plus `test_interrupt_point adaptive_refinement`: the kept checkpoint's stage, ties (union tied groups),
          projection ambiguity and explored labellings are recorded; F5 must report a genuine tie (>= 2 tied groups)
          and F1 several explored labellings when its checkpoint is the union's best.
  s2a_fixtures.py <out dir> <panvar>
"""
import csv
import os
import shlex
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "cyp2d6_pilot"))
import pilot  # noqa: E402

S1C = "/Users/davide.bolognini/panvar/tests/results/s1c_reporting_integrity_20261018/fixtures"
VOLATILE = ("seconds", "footprint", "rss", "peak", "_bytes", "worker", "threads_", "profile_sha256", "timing", "elapsed")
# the on/stop arms pass one more --profile (the early_incumbent key): the status names the last profile file
VOLATILE_EXACT = ("profile", "profile_name")
DECLARED = ("early_incumbent", "checkpoint_")


def kv(path):
    out = {}
    if os.path.exists(path):
        for line in open(path, errors="replace").read().splitlines():
            p = line.split("\t", 1)
            if len(p) == 2:
                out[p[0]] = p[1]
    return out


def run(cmd, log, timeout=240):
    with open(log, "w") as f:
        f.write("$ " + " ".join(cmd) + "\n")
        f.flush()
        code, _ = pilot.run_timed(cmd, f, timeout)
    return code


def compare(a, b):
    diffs = []
    fa = {f for f in os.listdir(a) if f.startswith("call.")}
    fb = {f for f in os.listdir(b) if f.startswith("call.")}
    if fa != fb:
        diffs.append(f"file sets {sorted(fa ^ fb)}")
    for f in sorted(fa & fb):
        if f in ("call.repair_memory.tsv",):
            continue
        ta = open(os.path.join(a, f), errors="replace").read().replace(a, "<D>")
        tb = open(os.path.join(b, f), errors="replace").read().replace(b, "<D>")
        if f.endswith("_status.tsv"):
            ka = {k: v.replace(a, "<D>") for k, v in kv(os.path.join(a, f)).items()
                  if not any(x in k for x in VOLATILE) and k not in VOLATILE_EXACT}
            kb = {k: v.replace(b, "<D>") for k, v in kv(os.path.join(b, f)).items()
                  if not any(x in k for x in VOLATILE) and k not in VOLATILE_EXACT and not k.startswith(DECLARED)}
            if ka != kb:
                diffs.append(f"{f}: {sorted(set(ka.items()) ^ set(kb.items()))[:2]}")
        elif f == "call.settings.tsv":
            if tb.replace("early_incumbent\texperimental_v1\t", "\x00").split("\x00")[0] != ta.split("early_incumbent")[0] or \
                    len(tb.splitlines()) != len(ta.splitlines()) + 1:
                diffs.append("settings")
        elif f == "call.repair_rounds.tsv":
            cut = lambda t: [line.split("\t")[:28] for line in t.splitlines()]  # noqa: E731
            if cut(ta) != cut(tb):
                diffs.append("trajectory")
        elif ta != tb:
            if f == "call.sequence_locus.tsv":
                ka, kb = kv(os.path.join(a, f)), kv(os.path.join(b, f))
                keys = sorted(k for k in set(ka) | set(kb) if ka.get(k, "").replace(a, "<D>") != kb.get(k, "").replace(b, "<D>"))
                diffs.append(f"{f} keys {','.join(keys[:6])}")
            else:
                diffs.append(f"{f} differs")
    return diffs


def hap_sha(d):
    p = os.path.join(d, "call.sequence_haplotypes.tsv")
    rows = list(csv.DictReader(open(p), delimiter="\t")) if os.path.exists(p) else []
    return "/".join(r["sha256"][:12] for r in rows) or "none"


def main(out, panvar):
    os.makedirs(out, exist_ok=True)
    on_profile = os.path.join(out, "profile_early.tsv")
    stop_profile = os.path.join(out, "profile_early_stop.tsv")
    open(on_profile, "w").write("profile\ts2a-early\nearly_incumbent\texperimental_v1\n")
    open(stop_profile, "w").write("profile\ts2a-early-stop\nearly_incumbent\texperimental_v1\n"
                                  "test_interrupt_point\tadaptive_refinement\n")
    rows = []
    for case in ("F1", "F5", "F7full", "F7"):
        base = shlex.split(open(os.path.join(S1C, case, "base.v1.log")).readline()[2:].strip())
        base[0] = panvar
        i = base.index("--dosage-mode")
        del base[i:i + 2]   # sequence_first_v2 refuses repair; its default is audit
        base[base.index("--output-schema") + 1] = "sequence_first_v2"
        o = base.index("--out")
        results = {}
        for arm, extra in (("off", []), ("on", ["--profile", on_profile]), ("stop", ["--profile", stop_profile])):
            d = os.path.join(out, case, arm)
            os.makedirs(d, exist_ok=True)
            for f in os.listdir(d):
                os.remove(os.path.join(d, f))
            cmd = list(base)
            cmd[o + 1] = os.path.join(d, "call")
            cmd[o:o] = extra
            results[arm] = run(cmd, os.path.join(out, case, f"{arm}.log"))
        off, on, stop = (os.path.join(out, case, a) for a in ("off", "on", "stop"))
        diffs = compare(off, on)
        loff, lon, lstop = (kv(os.path.join(d, "call.sequence_locus.tsv")) for d in (off, on, stop))
        son = kv(os.path.join(on, "call.mosaic_status.tsv"))
        completed = results["off"] == 0
        problems = []
        if completed and diffs:
            problems.append("completed run changed: " + " | ".join(diffs[:3]))
        if son.get("early_incumbent_union_check") not in ("reproduced",):
            problems.append(f"union check {son.get('early_incumbent_union_check')}")
        stage = lstop.get("checkpoint_stage")
        if results["stop"] != 3 or lstop.get("prediction") != "PROVISIONAL":
            problems.append(f"stop exit {results['stop']} {lstop.get('prediction')}")
        if case == "F5" and stage == "union_best" and int(lstop.get("checkpoint_tied_groups", "0")) < 2:
            problems.append("F5: the genuine tie is not reported")
        if case == "F1" and stage == "union_best" and int(lstop.get("checkpoint_explored_labellings", "0")) < 2:
            problems.append("F1: equivalent labellings not reported")
        rows.append({"case": case, "exit_off": results["off"], "exit_on": results["on"], "exit_stop": results["stop"],
                     "completed": "yes" if completed else "no (provisional control)",
                     "off_vs_on": "identical" if not diffs else " | ".join(diffs[:4]),
                     "prediction_off": loff.get("prediction"), "prediction_on": lon.get("prediction"),
                     "on_hap_sha": hap_sha(on), "off_hap_sha": hap_sha(off), "stop_hap_sha": hap_sha(stop),
                     "early_status": son.get("early_incumbent_status"), "union_check": son.get("early_incumbent_union_check"),
                     "checkpoint_commits": son.get("checkpoint_commits"), "stop_stage": stage,
                     "stop_tied_groups": lstop.get("checkpoint_tied_groups"),
                     "stop_projection_ambiguous": lstop.get("checkpoint_projection_ambiguous"),
                     "stop_explored_labellings": lstop.get("checkpoint_explored_labellings"),
                     "stop_evaluated": lstop.get("checkpoint_evaluated_candidates"),
                     "result": "PASS" if not problems else "FAIL: " + " | ".join(problems)})
        print("\t".join(str(v) for v in rows[-1].values()), flush=True)
    with open(os.path.join(out, "fixtures.tsv"), "w") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return 0 if all(r["result"] == "PASS" for r in rows) else 1


if __name__ == "__main__":
    sys.exit(main(os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2])))
