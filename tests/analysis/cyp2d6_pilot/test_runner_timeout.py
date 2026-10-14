#!/usr/bin/env python3
"""Synthetic checks of pilot.run_timed (the STOP_4 runner defect): a timed-out cell leaves no process behind.

A wrapper (standing in for /usr/bin/time) spawns a long-sleeping child (standing in for panvar) and writes the child's
pid. Checks:
  1. run_timed times out, and afterwards neither the wrapper nor the child exists;
  2. a wrapper that ignores SIGTERM is still removed (SIGKILL after the grace period), with its child;
  3. a cell that finishes normally returns its exit code and leaves nothing behind;
  4. negative control: the old wrapper-only behaviour (subprocess.run(..., timeout=...)) leaves the child alive.
Usage: test_runner_timeout.py
"""
import os
import signal
import subprocess
import sys
import tempfile
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pilot  # noqa: E402

WRAPPER = r'''
import signal, subprocess, sys, time
if sys.argv[2] == "ignore_term":
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(600)"])
open(sys.argv[1], "w").write(str(child.pid))
if sys.argv[2] == "quick":
    child.kill(); child.wait(); sys.exit(3)
time.sleep(600)
'''


def alive(pid):
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    # a zombie still answers kill(0); count it as gone once its state is Z
    state = subprocess.run(["ps", "-o", "stat=", "-p", str(pid)], capture_output=True, text=True).stdout.strip()
    return bool(state) and not state.startswith("Z")


def gone_within(pid, seconds):
    deadline = time.monotonic() + seconds
    while alive(pid):
        if time.monotonic() > deadline:
            return False
        time.sleep(0.1)
    return True


def wait_for(path):
    deadline = time.monotonic() + 10
    while not (os.path.exists(path) and open(path).read().strip()):
        if time.monotonic() > deadline:
            raise RuntimeError("wrapper did not report its child")
        time.sleep(0.05)
    return int(open(path).read())


failures = []


def check(ok, what):
    print(("PASS " if ok else "FAIL ") + what)
    if not ok:
        failures.append(what)


def main():
    tmp = tempfile.mkdtemp(prefix="runner_timeout_")
    wrapper = os.path.join(tmp, "wrapper.py")
    open(wrapper, "w").write(WRAPPER)
    with open(os.path.join(tmp, "log"), "w") as log:
        for mode in ("plain", "ignore_term"):
            pidfile = os.path.join(tmp, f"child_{mode}.pid")
            t0 = time.monotonic()
            code, timed_out = pilot.run_timed([sys.executable, wrapper, pidfile, mode], log, timeout=2, term_grace=1)
            child = wait_for(pidfile)
            check(timed_out and code == "timeout", f"[{mode}] the cell times out ({time.monotonic() - t0:.1f} s)")
            check(gone_within(child, 5), f"[{mode}] the child (pid {child}) is gone after run_timed returns")
        pidfile = os.path.join(tmp, "child_quick.pid")
        code, timed_out = pilot.run_timed([sys.executable, wrapper, pidfile, "quick"], log, timeout=30)
        child = wait_for(pidfile)
        check(code == 3 and not timed_out and not alive(child), "a normally finishing cell returns its exit code, nothing left")
        # negative control: the old runner killed only the direct child (the wrapper)
        pidfile = os.path.join(tmp, "child_old.pid")
        try:
            subprocess.run([sys.executable, wrapper, pidfile, "plain"], stdout=log, stderr=subprocess.STDOUT, timeout=2)
        except subprocess.TimeoutExpired:
            pass
        child = wait_for(pidfile)
        time.sleep(1)
        orphaned = alive(child)
        check(orphaned, f"negative control: the old wrapper-only kill leaves the child (pid {child}) alive")
        if orphaned:
            os.kill(child, signal.SIGKILL)
    print(f"{len(failures)} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
