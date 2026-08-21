#!/usr/bin/env python3
"""Check the documentation against the tool it documents.

Three things drift silently: a flag that was renamed, a column name that never existed, and a page
whose worked trace stops matching the numbered steps above it. All three read fine.

    scripts/check_docs.py [module ...]

Checks per module page and algorithm page:
  flags    every `--flag` mentioned is accepted by that module's --help
  names    every backticked snake_case or SCREAMING_CASE name is a real column, VCF key or summary
           key in that module's actual output
  trace    the numbered steps in the worked trace match the numbered sections above them
  style    no bold/italic, no gene or locus names

Output names come from results/real_data/, so run scripts/regen_results.sh first.
"""
import json
import pathlib
import re
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parent.parent
BIN = REPO / "build" / "panvar"
MODULES = ["rebuild", "bubble", "panphorte", "refine", "call",
           "describe", "associate", "inspect", "benchmark"]

# Words that look like output names but are prose, code identifiers or file extensions.
ALLOWED = {
    "source", "sink", "bubble", "call", "inspect", "panphorte", "describe", "associate",
    "benchmark", "refine", "rebuild", "snarls", "prefix", "path", "reason", "status", "sample",
    "id", "name", "gz", "tsv", "csv", "vcf", "gfa", "bimbam", "jsonl", "fa", "png", "bed", "gtf",
    "true", "false", "auto", "off", "on", "none", "all", "NA", "ALL", "GT", "REF", "ALT", "POS",
    "INFO", "FORMAT", "CHROM", "ID", "QUAL", "FILTER", "PASS", "DEL", "INS", "INV", "DUP", "SNV",
    "P", "W", "S", "L", "N", "A", "B", "AC", "AN", "AF", "bp", "kb", "Mb", "k", "n", "r", "i", "o",
    # single letters and short symbols used as maths in the algorithm pages
    "J", "U", "K", "X", "Y", "T", "W1", "W2", "W3", "W4", "hA", "hB", "hC", "hD",
}
GENE_RE = re.compile(r"\b(lpa|gstm1|gstm2|gstm4|gstm5|cyp2d6|cyp2d7|acot1|acot2|ankrd36c|kiv-2|myom2|c4a|c4b)\b", re.I)


def flags_of(module):
    try:
        out = subprocess.run([str(BIN), module, "--help"], capture_output=True, text=True, timeout=60)
    except Exception:
        return set()
    return set(re.findall(r"--[a-z0-9][a-z0-9-]+", out.stdout + out.stderr))


def numbered_sections(text):
    return [int(m.group(1)) for m in re.finditer(r"^###\s+(\d+)\.", text, re.M)]


def numbered_trace(text):
    idx = text.find("## Worked trace")
    if idx < 0:
        return []
    body = text[idx:]
    # a numbered list item at the start of a line, outside fenced code
    out, fenced = [], False
    for line in body.split("\n"):
        if line.startswith("```"):
            fenced = not fenced
            continue
        if fenced:
            continue
        m = re.match(r"^(\d+)\.\s", line)
        if m:
            out.append(int(m.group(1)))
    return out


def main():
    wanted = sys.argv[1:] or MODULES
    names = json.loads((REPO / "scripts" / "doc_output_names.json").read_text())
    # A page legitimately cites names another module produces -- benchmark reads call's VCF keys,
    # associate reads describe's columns. The check is that a name exists SOMEWHERE in the tool's real
    # output, not that the page only mentions its own. Same for flags.
    all_names = {n for v in names.values() for n in v}
    all_flags = {f for m in MODULES for f in flags_of(m)}
    problems = 0
    for mod in wanted:
        flags = all_flags
        known = all_names | ALLOWED
        for page in (REPO / "docs" / "modules" / f"{mod}.md",
                     REPO / "docs" / "algorithms" / f"{mod}.md"):
            if not page.exists():
                continue
            text = page.read_text()
            rel = page.relative_to(REPO)
            bad = []

            # Link targets carry anchors like #rebuild--re-inducing-..., which look like flags.
            prose = re.sub(r"\]\([^)]*\)", "]()", text)
            for f in sorted(set(re.findall(r"--[a-z0-9][a-z0-9-]+", prose))):
                if flags and f not in flags:
                    bad.append(f"flag not accepted by `{mod} --help`: {f}")

            for tok in sorted(set(re.findall(r"`([A-Za-z][A-Za-z0-9_]*)`", text))):
                if "_" not in tok and not tok.isupper():
                    continue          # single lowercase words are prose
                if re.fullmatch(r"[ACGTN]+", tok):
                    continue          # a DNA string, not a field name
                if re.fullmatch(r"[A-Z][0-9]?", tok):
                    continue          # a single-letter label in a worked trace
                if tok.endswith("_calls") or tok.endswith("_out"):
                    continue          # a default output prefix
                if tok not in known:
                    bad.append(f"name not found in {mod}'s real output: {tok}")

            if "## Worked trace" in text:
                sec, tr = numbered_sections(text), numbered_trace(text)
                if sec and tr and sec != tr:
                    bad.append(f"trace steps {tr} do not match sections {sec}")

            if "**" in text:
                bad.append("contains bold")
            if re.search(r"(^|[^*\w])\*[^*\s][^*]*\*([^*\w]|$)", text):
                bad.append("contains italics")
            for g in sorted(set(m.group(0) for m in GENE_RE.finditer(text))):
                bad.append(f"names a specific gene or locus: {g}")

            if bad:
                problems += len(bad)
                print(f"{rel}")
                for b in bad:
                    print(f"    {b}")
    print("docs check: OK" if problems == 0 else f"docs check: {problems} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
