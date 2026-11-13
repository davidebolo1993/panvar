#!/usr/bin/env python3
"""Precision, recall and F1 of call and svim-asm against benchmark's <prefix>.present.tsv: every difference
of --min-sv-bp or more between the reference and each haplotype, aligned directly.

Each caller is scored at two levels, as benchmark scores call:
  own      what it found in each haplotype: call's carrier_events.tsv; svim-asm's per-haplotype calls
           (assembly_variants.py --sv-tables: <prefix>.haplotype_svs.tsv)
  merged   what its multi-sample output says each haplotype carries: call's region.vcf (GT, and CN against
           REF_CN on copy-number records); svim-asm after truvari collapse (<prefix>.merged_svs.tsv)

Matching is implemented here (truvari bench cannot read present.tsv, per-haplotype carriers or copy
numbers), following truvari bench's rules and defaults but without its sequence comparison, since in a
repeat a copy can sit at several offsets that spell nearly the same sequence. A call matches a truth event
of the same haplotype when:
  - their intervals overlap, or lie within --refdist;
  - their sizes (the larger of bases removed and inserted) agree to --pctsize; a replacement may also match
    by its net change, and a deletion under 1 kb by the reference it spans;
  - their directions agree: a loss with a loss, a gain with a gain; a replacement or inversion with either.
Pairs are matched one to one first, closest first. Then, like truvari's multi-match, an event claims a
consecutive run of unmatched events of the other side (same haplotype and direction, near it) once together
they reach its size: in a tandem array one copy-number change can be several truth events. Only the
haplotypes and stretches benchmark scored count. In call's own events, a carrier's events inside its copy
stretch belong to the stretch; in the VCF a copy-number record is a count, so records inside it stay calls.

  compare_sv_callers.py --label C4 --benchmark bench/c4 --calls call/c4/c4 --svim svim/c4 -o c4.callers.tsv
                        [--details c4.callers.events.tsv] [--refdist 500] [--pctsize 0.7] [--min-sv-bp 50]

Output: one row per caller, level and stratum (all, and by size: of the truth event for recall, of the call
for precision: 50-999 bp, 1-10 kb, over 10 kb): truth events and those matched, calls and those matched,
precision, recall, F1.
"""
import argparse
import collections
import csv
import sys

SIZE_BINS = [("50-999 bp", 0, 1000), ("1-10 kb", 1000, 10000), (">10 kb", 10000, float("inf"))]


class Event:
    """One SV of one haplotype. `size` is the larger of the bases removed and inserted; `sizes` also holds
    the net change where that can be the event's real size (a replacement, or an event written with the
    bases around it as context), and any size matches."""
    __slots__ = ("hap", "lo", "hi", "size", "sizes", "net", "kind", "matched", "label", "index", "point")

    def __init__(self, hap, lo, hi, size, net, kind, label="", point=False):
        self.hap, self.lo, self.hi = hap, min(lo, hi), max(lo, hi)
        self.size, self.net, self.kind, self.label = size, net, kind, label
        self.point = point   # an insertion: goes in before `lo`, replacing no reference base
        self.sizes = {size} | ({abs(net)} if kind == "any" and net else set())
        self.matched = False
        self.index = 0


def kind_of(removed, inserted):
    if inserted == 0:
        return "loss"
    if removed == 0:
        return "gain"
    return "any"


def read_tsv(path):
    csv.field_size_limit(sys.maxsize)   # call's events hold whole copy stretches
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def truth_events(prefix, min_sv):
    stretch = {}
    for r in read_tsv(prefix + ".haplotypes.tsv"):
        if r["status"] == "scored":
            stretch[r["haplotype"]] = (int(r["ref_start"]), int(r["ref_end"]))
    events = []
    for r in read_tsv(prefix + ".present.tsv"):
        if r["haplotype"] not in stretch:
            continue
        d, i = int(r["del_bp"]), int(r["ins_bp"])
        if max(d, i) < min_sv:
            continue
        kind = {"DEL": "loss", "INS": "gain"}.get(r["type"], "any")
        # present.tsv spans the bases around the change; here an event is its first and last changed base,
        # an insertion the base it goes before
        g0, g1 = int(r["ref_start"]), int(r["ref_end"])
        lo, hi = (g0 + 1, g1 - 1) if g1 - 1 >= g0 + 1 else (g0 + 1, g0 + 1)
        e = Event(r["haplotype"], lo, hi, max(d, i), i - d, kind, r["type"])
        # A deletion in a short repeat can align with chance matches inside its gap, so its reference span is
        # also its size (as truvari sizes deletions); only under 1 kb and up to 2.5 times its gap bases.
        span = e.hi - e.lo + 1
        if r["type"] == "DEL" and e.size < 1000 and e.size < span <= 2.5 * e.size:
            e.sizes.add(span)
        events.append(e)
    return stretch, events


def inside_at(c, e):
    """benchmark's rule for an event inside a copy stretch: an interval within it, an insertion going in
    anywhere from before its first base to before the base after its last. By the listed place only:
    benchmark also counts an insertion whose CIPOS crosses the stretch's edge as inside when that brings the
    stretch's length nearer its own, which needs the copy forms it reads from the reference."""
    if c.point:
        return e.point and e.lo == c.lo
    if e.point:
        return c.lo <= e.lo <= c.hi + 1
    return c.lo <= e.lo and e.hi <= c.hi


def without_inner(events, min_sv):
    """A carrier's events inside its own copy-number stretch are part of it: the stretch's net change already
    holds them (benchmark treats them so too), so they are not separate calls. A stretch whose net change is
    under min_sv still holds them, but is no call itself."""
    stretches = collections.defaultdict(list)
    for e in events:
        if e.label.startswith("cn:"):
            stretches[e.hap].append(e)
    out = []
    for e in events:
        if e.label.startswith("cn:"):
            if abs(e.net) >= min_sv:
                out.append(e)
        elif not any(inside_at(c, e) for c in stretches[e.hap]):
            out.append(e)
    for e in out:
        e.label = e.label[3:] if e.label.startswith("cn:") else e.label
    return out


def call_own(prefix, min_sv):
    events = []
    for r in read_tsv(prefix + ".carrier_events.tsv"):
        if r["kind"] == "unresolved" or r["event_bp"] in (".", ""):
            continue   # an event call could not place: no position to match
        lo, hi, net = int(r["ref_start"]), int(r["ref_end"]), int(r["event_bp"])
        point = hi < lo
        tag = ""
        if r["kind"] == "cn":   # a copy stretch: the copies gained or lost, anywhere over it
            size, kind, tag = abs(net), ("gain" if net > 0 else "loss"), "cn:"
        elif r["kind"] == "inv":
            size, kind = hi - lo + 1, "any"
        else:
            # The replaced interval can carry unchanged bases as context, so the net change is also a size.
            removed = max(0, hi - lo + 1)
            inserted = int(r["alt_bp"]) if r["alt_bp"] not in (".", "") else 0
            size, kind = max(removed, inserted), kind_of(removed, inserted)
            if kind != "any" and abs(net) != size:
                kind = "any"
        if max(size, abs(net)) < min_sv and not tag:
            continue   # a copy stretch stays until the events inside it are taken (without_inner)
        if hi < lo:
            hi = lo
        for hap in r["haplotypes"].split(","):
            if hap:
                events.append(Event(hap, lo, hi, size, net, kind, tag + r["variant_id"], point))
    return events


def info_dict(field):
    out = {}
    for kv in field.split(";"):
        k, _, v = kv.partition("=")
        out[k] = v
    return out


def call_merged(prefix, min_sv):
    events = []
    samples = []
    with open(prefix + ".region.vcf") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line or line.startswith("##"):
                continue
            f = line.split("\t")
            if line.startswith("#"):
                samples = f[9:]
                continue
            inf = info_dict(f[7])
            t = inf.get("SVTYPE", "")
            pos = int(f[1])
            end = int(inf["END"]) if "END" in inf else pos + len(f[3]) - 1
            alts = f[4].split(",")
            svlens = [int(x) for x in inf["SVLEN"].split(",")] if inf.get("SVLEN") else []
            fmt = f[8].split(":")
            gi = fmt.index("GT")
            ci = fmt.index("CN") if "CN" in fmt else -1
            bi = fmt.index("CNBP") if "CNBP" in fmt else -1
            ref_cn = int(inf["REF_CN"]) if inf.get("REF_CN", "").lstrip("-").isdigit() else None
            unit = next((int(inf[k]) for k in ("CN_UNIT_BP", "RU_LEN") if inf.get(k, "").isdigit()), 0)
            span = end - pos   # a CNV covers the reference's copies, POS+1..END
            for hap, cell in zip(samples, f[9:]):
                v = cell.split(":")
                if t == "CNV":   # the copies changed, anywhere over the reference's copies
                    if ci < 0 or ref_cn is None or len(v) <= ci or v[ci] in (".", "") or int(v[ci]) == ref_cn:
                        continue
                    cn = int(v[ci])
                    if 0 <= bi < len(v) and v[bi] not in (".", ""):
                        net = int(v[bi])
                    else:   # no CNBP: the copies counted, each a unit long
                        net = (cn - ref_cn) * unit
                    events.append(Event(hap, pos + 1, pos + span, abs(net), net, "gain" if cn > ref_cn else "loss", "cn:" + f[2]))
                    continue
                a = v[gi]
                if not a.isdigit() or int(a) == 0 or int(a) > len(alts):
                    continue
                k = int(a) - 1
                svlen = svlens[k] if k < len(svlens) else None
                alt = alts[k]
                if t == "DEL":
                    size = abs(svlen) if svlen is not None else end - pos
                    events.append(Event(hap, pos + 1, max(end, pos + 1), size, -size, "loss", f[2]))
                elif t == "INS":
                    size = abs(svlen) if svlen is not None else len(alt) - len(f[3])
                    events.append(Event(hap, pos + 1, pos + 1, size, size, "gain", f[2]))
                elif t == "INV":
                    events.append(Event(hap, pos + 1, max(end, pos + 1), end - pos, 0, "any", f[2]))
                else:   # SUB, written as REF and ALT; the leading bases they share are padding
                    ref = f[3]
                    pad = 0
                    while pad < min(len(ref), len(alt)) and ref[pad].upper() == alt[pad].upper():
                        pad += 1
                    removed, inserted = len(ref) - pad, len(alt) - pad
                    lo = pos + pad
                    events.append(Event(hap, lo, max(lo, pos + len(ref) - 1), max(removed, inserted),
                                        inserted - removed, "any", f[2], removed == 0))
    return [e for e in events if max(e.sizes) >= min_sv]


def svim_events(rows, min_sv, haps_of):
    events = []
    for r in rows:
        t, size = r["type"], int(r["length"])
        if size < min_sv:
            continue
        lo, hi = int(r["pos"]), int(r["end"])
        kind = {"DEL": "loss", "INS": "gain", "DUP": "gain"}.get(t, "any")
        net = -size if kind == "loss" else (size if kind == "gain" else 0)
        for hap in haps_of(r):
            events.append(Event(hap, lo + 1, max(hi, lo + 1) if t != "INS" else lo + 1, size, net, kind, r.get("id", t)))
    return events


def compatible(t, c, refdist, pctsize):
    if t.hap != c.hap:
        return None
    gap = max(0, max(t.lo, c.lo) - min(t.hi, c.hi))
    if gap > refdist:
        return None
    if t.kind != "any" and c.kind != "any" and t.kind != c.kind:
        return None
    sim = max((min(a, b) / max(a, b) for a in t.sizes for b in c.sizes if a > 0 and b > 0), default=0.0)
    if sim < pctsize:
        return None
    return (gap, 1 - sim)


def match(truth, calls, refdist, pctsize):
    """The most matches one to one (augmenting paths), each truth event trying its closest calls first; ties
    broken by file order, so a run is reproducible."""
    for i, e in enumerate(truth):
        e.matched, e.index = False, i
    for i, e in enumerate(calls):
        e.matched, e.index = False, i
    by_hap = collections.defaultdict(list)
    for c in calls:
        by_hap[c.hap].append(c)
    edges = {}
    for t in truth:
        opts = []
        for c in by_hap.get(t.hap, ()):
            key = compatible(t, c, refdist, pctsize)
            if key is not None:
                opts.append((key, c.index, c))
        opts.sort(key=lambda o: (o[0], o[1]))
        edges[t.index] = [c for _, _, c in opts]
    owner = {}   # call index -> truth event
    order = sorted((t for t in truth if edges[t.index]),
                   key=lambda t: (compatible(t, edges[t.index][0], refdist, pctsize), t.index))

    def augment(t, seen):
        for c in edges[t.index]:
            if c.index in seen:
                continue
            seen.add(c.index)
            if c.index not in owner or augment(owner[c.index], seen):
                owner[c.index] = t
                return True
        return False

    for t in order:
        augment(t, set())
    for ci, t in owner.items():
        t.matched = True
        calls[ci].matched = True
    group_match(truth, calls, refdist, pctsize)


def group_match(truth, calls, refdist, pctsize):
    """Inside a tandem array the direct alignment places each gained or lost copy on its own, often between
    reference copies, so one copy-number change is several truth events (and a caller can split one event
    too). What one-to-one matching left is matched one to several: an event claims unmatched events of the
    other side, same haplotype and direction, within --refdist plus its own size of it (a point insertion
    says nothing of where in the array its copies sit), lying together, the run nearest it first, once
    together they reach its size within --pctsize, and no further."""
    def claim(one, others):
        size = one.size   # its own size: the alternatives are for pairing one to one
        reach = refdist + size
        gap_to = lambda lo, hi: max(0, max(one.lo, lo) - min(one.hi, hi))
        for direction in ("loss", "gain"):
            if one.kind not in ("any", direction):
                continue
            # The pieces: of this direction, each no larger than the event allows, near it.
            near = [o for o in others if not o.matched and o.hap == one.hap and o.kind in (direction, "any")
                    and o.size <= size / pctsize and gap_to(o.lo, o.hi) <= reach]
            if len(near) < 2:
                continue
            # The pieces of one change lie together: a chain along the reference, each within --refdist or its
            # own size of what the chain covers so far.
            near.sort(key=lambda o: (o.lo, o.hi, -o.size))
            chains, cur, cur_hi = [], [near[0]], near[0].hi
            for o in near[1:]:
                if o.lo - cur_hi <= max(refdist, cur[-1].size, o.size):
                    cur.append(o)
                    cur_hi = max(cur_hi, o.hi)
                else:
                    chains.append((cur, cur_hi))
                    cur, cur_hi = [o], o.hi
            chains.append((cur, cur_hi))
            # Consecutive pieces of a chain, the run nearest the event first, until they reach its size.
            runs = []
            for chain, _ in chains:
                for a in range(len(chain)):
                    total, hi = 0, chain[a].hi
                    for b in range(a, len(chain)):
                        total += chain[b].size
                        hi = max(hi, chain[b].hi)
                        if total > size / pctsize:
                            break
                        if b > a and total >= size * pctsize:
                            runs.append((gap_to(chain[a].lo, hi), chain[a].lo, chain[a:b + 1]))
                            break
            if runs:
                picked = min(runs, key=lambda r: (r[0], r[1]))[2]
                one.matched = True
                for o in picked:
                    o.matched = True
                return

    by_hap_t, by_hap_c = collections.defaultdict(list), collections.defaultdict(list)
    for t in truth:
        by_hap_t[t.hap].append(t)
    for c in calls:
        by_hap_c[c.hap].append(c)
    for c in sorted(calls, key=lambda e: (-e.size, e.lo, e.hi)):   # the largest calls first
        if not c.matched:
            claim(c, by_hap_t[c.hap])
    for t in sorted(truth, key=lambda e: (-e.size, e.lo, e.hi)):
        if not t.matched:
            claim(t, by_hap_c[t.hap])


def in_stretch(events, stretch):
    """Inside the stretch benchmark scored, as it decides: an interval within it, a point past its first base."""
    out = []
    for e in events:
        s = stretch.get(e.hap)
        if s is None:
            continue
        # an insertion (a point: the base it goes before) past the first base; any other event within it
        if (e.lo > s[0] and e.lo <= s[1]) if e.lo == e.hi else (e.lo >= s[0] and e.hi <= s[1]):
            out.append(e)
    return out


def metrics(truth, calls):
    nt, mt = len(truth), sum(t.matched for t in truth)
    nc, mc = len(calls), sum(c.matched for c in calls)
    p = mc / nc if nc else float("nan")
    r = mt / nt if nt else float("nan")
    f1 = (2 * p * r / (p + r) if (p + r) > 0 else 0.0) if nc and nt else float("nan")
    return nt, mt, nc, mc, p, r, f1


def main():
    ap = argparse.ArgumentParser(description="Precision, recall and F1 of call and svim-asm against benchmark's present.tsv.",
                                 formatter_class=lambda prog: argparse.HelpFormatter(prog, max_help_position=32, width=100))
    ap.add_argument("--label", required=True, help="locus name in the table")
    ap.add_argument("--benchmark", required=True, metavar="PREFIX", help="benchmark's -o prefix")
    ap.add_argument("--calls", required=True, metavar="PREFIX", help="call's -o prefix")
    ap.add_argument("--svim", required=True, metavar="PREFIX", help="assembly_variants.py --sv-tables prefix")
    ap.add_argument("-o", "--out", required=True, metavar="TSV", help="output table")
    ap.add_argument("--details", metavar="TSV", help="also every event and whether it matched")
    ap.add_argument("--refdist", type=int, default=500, metavar="BP", help="largest distance to match (default 500)")
    ap.add_argument("--pctsize", type=float, default=0.7, metavar="F", help="smallest size ratio to match (default 0.7)")
    ap.add_argument("--min-sv-bp", type=int, default=50, metavar="N", help="smallest event (default 50)")
    a = ap.parse_args()

    stretch, truth = truth_events(a.benchmark, a.min_sv_bp)
    if not stretch:
        sys.exit("compare_sv_callers: no scored haplotype in " + a.benchmark + ".haplotypes.tsv")
    merged_rows = read_tsv(a.svim + ".merged_svs.tsv")
    callers = {
        ("panvar call", "own"): call_own(a.calls, a.min_sv_bp),
        ("panvar call", "merged"): call_merged(a.calls, a.min_sv_bp),
        ("svim-asm", "own"): svim_events(read_tsv(a.svim + ".haplotype_svs.tsv"), a.min_sv_bp, lambda r: [r["haplotype"]]),
        ("svim-asm", "merged"): svim_events(merged_rows, a.min_sv_bp, lambda r: [h for h in r["carriers"].split(",") if h]),
    }
    rows, details = [], []
    for (caller, level), calls in callers.items():
        calls = in_stretch(calls, stretch)
        if (caller, level) == ("panvar call", "own"):
            # a copy stretch's own sequence holds the carrier's events inside it; a merged record holds a count,
            # so there the records inside it stay calls of their own
            calls = without_inner(calls, a.min_sv_bp)
        for c in calls:
            c.label = c.label[3:] if c.label.startswith("cn:") else c.label
        match(truth, calls, a.refdist, a.pctsize)
        strata = [("all", truth, calls)]
        for name, lo, hi in SIZE_BINS:
            strata.append((name, [t for t in truth if lo <= t.size < hi], [c for c in calls if lo <= c.size < hi]))
        for name, ts, cs in strata:
            nt, mt, nc, mc, p, r, f1 = metrics(ts, cs)
            rows.append([a.label, caller, level, name, nt, mt, nc, mc] + ["%.4f" % x if x == x else "NA" for x in (p, r, f1)])
        if a.details:
            for t in truth:
                details.append([a.label, caller, level, "truth", t.hap, t.lo, t.hi, t.label, t.kind, t.size, int(t.matched)])
            for c in calls:
                details.append([a.label, caller, level, "call", c.hap, c.lo, c.hi, c.label, c.kind, c.size, int(c.matched)])
    with open(a.out, "w") as fh:
        fh.write("locus\tcaller\tlevel\tstratum\ttruth_events\ttruth_matched\tcalls\tcalls_matched\tprecision\trecall\tf1\n")
        for r in rows:
            fh.write("\t".join(map(str, r)) + "\n")
    if a.details:
        with open(a.details, "w") as fh:
            fh.write("locus\tcaller\tlevel\tside\thaplotype\tref_start\tref_end\tevent\tdirection\tsize_bp\tmatched\n")
            for r in details:
                fh.write("\t".join(map(str, r)) + "\n")
    for r in rows:
        if r[3] == "all":
            print("%s %-12s %-6s truth %4d (%4d found)  calls %4d (%4d right)  precision %s  recall %s  F1 %s" % tuple(r[0:3] + r[4:]))


if __name__ == "__main__":
    main()
