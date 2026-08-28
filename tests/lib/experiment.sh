# experiment.sh - contract for experimental harnesses. Source it, do not execute it.
#
#   source "$(dirname "${BASH_SOURCE[0]}")/lib/experiment.sh"
#   exp_init genotype-frag-multiseed "$BIN"
#   exp_run  "$OUT/arm1" "$OUT/arm1.hap_pairs.tsv" -- "$BIN" genotype-frag ...
#
# WHY THIS EXISTS. The same failure has now occurred twice in this project and both times it produced
# a confident null result that was wrong:
#
#   * a run failed, its output file from a PREVIOUS arm was still on disk, and the stale file was read
#     back as that arm's result -- so a live flag looked inert and "the channel does nothing" was
#     nearly recorded as a finding;
#   * an audit column was read from an instrument that had never run, and reported 0 rather than NA.
#
# A null result is only evidence if the arm demonstrably ran. So: every arm gets a unique prefix, that
# prefix is cleared before the run, the exit status is checked, and every expected output is asserted
# to exist and be non-empty before anything reads it. The binary hash and full argument list are
# recorded so a number can be traced to what produced it.
set -euo pipefail

EXP_NAME=""
EXP_LOG=""
EXP_BIN=""

exp_init() {   # <experiment-name> <binary> [logfile]
    EXP_NAME="${1:?exp_init needs a name}"
    EXP_BIN="${2:?exp_init needs the binary under test}"
    EXP_LOG="${3:-${OUT:-.}/provenance.txt}"
    mkdir -p "$(dirname "$EXP_LOG")"
    {
        printf 'experiment\t%s\n' "$EXP_NAME"
        printf 'date\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'binary\t%s\n' "$EXP_BIN"
        printf 'binary_md5\t%s\n' "$(md5 -q "$EXP_BIN" 2>/dev/null || md5sum "$EXP_BIN" | cut -d' ' -f1)"
        printf 'commit\t%s\n' "$(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
        printf 'dirty\t%s\n' "$(git status --porcelain 2>/dev/null | wc -l | tr -d ' ')"
    } > "$EXP_LOG"
    printf 'binary %s  commit %s\n' \
        "$(md5 -q "$EXP_BIN" 2>/dev/null | cut -c1-8 || echo ?)" \
        "$(git rev-parse --short HEAD 2>/dev/null || echo ?)" >&2
}

# exp_run <prefix> <expected-output>... -- <command>...
# Clears the prefix, runs the command, and refuses to continue unless it exited 0 and every named
# output exists and is non-empty. Nothing downstream may read a file this did not certify.
exp_run() {
    local prefix="${1:?exp_run needs an output prefix}"; shift
    local -a expected=()
    while [[ $# -gt 0 && "$1" != "--" ]]; do expected+=("$1"); shift; done
    [[ "${1:-}" == "--" ]] || { echo "exp_run: missing -- before the command" >&2; return 2; }
    shift
    rm -f "$prefix".* 2>/dev/null || true
    for f in "${expected[@]}"; do rm -f "$f" 2>/dev/null || true; done
    printf 'cmd\t%s\n' "$*" >> "$EXP_LOG"
    local rc=0
    "$@" >/dev/null 2>>"$EXP_LOG" || rc=$?
    if [[ $rc -ne 0 ]]; then
        printf '  FAIL rc=%d: %s\n' "$rc" "$*" >&2
        return "$rc"
    fi
    for f in "${expected[@]}"; do
        if [[ ! -s "$f" ]]; then
            printf '  FAIL missing or empty output: %s\n' "$f" >&2
            return 3
        fi
    done
    return 0
}

# A value that must exist. Empty or literal "NA" is a hard stop, not a 0 to be averaged in later.
exp_require_value() {   # <label> <value>
    if [[ -z "${2:-}" || "${2}" == "NA" ]]; then
        printf '  FAIL %s is empty/NA; refusing to record it as a number\n' "$1" >&2
        return 4
    fi
    printf '%s' "$2"
}
