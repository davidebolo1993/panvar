#!/usr/bin/env python3
"""M10-B.1 structural-projection prototype, version 2 (diagnostic only; project.py, version 1, is kept unchanged).

Changes from version 1:
  1. Missing metadata is never structure. An allele whose recorded walks are fewer than its walk source records, or that
     has no walk record and no valid empty record ('empty' / 'bypass'), is MISSING; any projection using it is
     INCOMPLETE. A valid empty record contributes an empty walk. (Version 1 turned an absent walk into a RESOLVED skeleton
     holding the token NO_WALK, so unrelated sequences could compare EQUAL.)
  2. Class identity is the canonical node SEQUENCE BYTES. The short display id (a hash prefix) is assigned by a registry
     that disambiguates any prefix collision; the hash function is injectable for the forced-collision gate.
  3. Coordinates are reconciled with the implementation's boundary transformation (genotype_locus.cpp):
       - a variant block's recorded walk includes its bubble boundary nodes; the divergence trim then removes the bytes
         common to all its source alleles: normalized = spell(walk)[P : len - S], with P / S the report's
         trimmed_prefix_bp / trimmed_suffix_bp (no trim in a block carrying a bypass);
       - the removed bytes move to the neighbouring context blocks: context normalized = S(prev) + spell(walk) + P(next);
       - an abutting bubble drops the shared boundary node from one side, so consecutive walks share no node.
     The relation is VERIFIED per block (audit()) for every allele and walk before any projection uses that block.
     A scope f..l is projected in normalized coordinates as the concatenated walk W clipped by P(f) (if f is a variant
     block) and S(l) (if l is a variant block). The scope is UNSUPPORTED when
       - bytes enter it from outside (f is a context block after a trimmed variant block, or l a context block before
         one): their nodes belong to a neighbour's walk, which the scope does not determine;
       - a skeleton node (long, or declared) straddles a clip boundary;
       - a block in the scope failed the audit.
  4. Comparison categories: UNIQUE_AGREEMENT, SET_IDENTICAL (equal multisets, some projection set-valued), COMPATIBLE,
     DISAGREEMENT (with the DOSAGE_EQUAL_ALLOCATION_DIFFERS flag), UNSUPPORTED (any side UNKNOWN, INCOMPLETE,
     UNSUPPORTED or OVERFLOW). Only UNIQUE_AGREEMENT establishes one structural genotype.
Unchanged: pi is defined on sequences (the set of skeletons over every exact parse and every walk); ANCHOR_BP = 50 is the
configured min_variant_bp used as an experimental convention; class-defining sites are an explicit, separate declaration.
"""
import collections
import hashlib
import itertools
import os

ANCHOR_BP = 50
MAX_PARSES = 100000
COMP = str.maketrans("ACGTN", "TGCAN")
FAILING = ("UNKNOWN", "INCOMPLETE", "UNSUPPORTED", "OVERFLOW")


def rc(s):
    return s.translate(COMP)[::-1]


def T(p):
    import csv
    return list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t")) if os.path.exists(p) else []


def default_hash(b):
    return hashlib.sha256(b.encode()).hexdigest()[:10]


class Model:
    def __init__(self, prefix, class_defining_sequences=(), hasher=default_hash):
        self.prefix = prefix
        self.nodes = {r["node"]: r["sequence"] for r in T(prefix + ".nodes.tsv")}
        self.seq, self.block_ids, self.kind = {}, [], []
        for r in T(prefix + ".block_alleles.tsv"):
            b = int(r["block_index"])
            self.seq[(b, int(r["allele_index"]))] = r["sequence"]
            if len(self.block_ids) <= b:
                self.block_ids.append(r["block_id"])
                self.kind.append(r["kind"])
        self.walks = collections.defaultdict(list)
        for r in T(prefix + ".walks.tsv"):
            self.walks[(int(r["block_index"]), int(r["allele_index"]))].append(r["walk"].split(",") if r["walk"] else [])
        self.sources = {(int(r["block_index"]), int(r["allele_index"])): {k: int(r[k]) for k in ("walk_sources", "empty_sources",
                        "bypass_sources", "other_sources")} for r in T(prefix + ".sources.tsv")}
        self.trim = {int(r["block_index"]): (int(r["trimmed_prefix_bp"]), int(r["trimmed_suffix_bp"])) for r in T(prefix + ".blocks.tsv")
                     if r["trimmed_prefix_bp"] not in ("NA", "")}
        self.alleles = collections.defaultdict(list)
        for (b, a), s in sorted(self.seq.items()):
            self.alleles[b].append((a, s))
        self.panel = {r["template"]: [int(x) for x in r["alleles"].split(",")] for r in T(prefix + ".panel.tsv")}
        self.fingerprint = {r["key"]: r["value"] for r in T(prefix + ".summary.tsv")}.get("model_fingerprint")
        self.hasher = hasher
        self.registry = {}                     # canonical bytes -> display id (identity is the bytes)
        self.used_ids = {}
        self.class_defining = {min(s, rc(s)) for s in class_defining_sequences}
        self._audit = None

    # ---- class identity (byte-verified)
    def class_of(self, s):
        canon = min(s, rc(s))
        if canon not in self.registry:
            base = "U" + self.hasher(canon)
            cid, n = base, 1
            while cid in self.used_ids and self.used_ids[cid] != canon:   # a prefix collision: disambiguate, never merge
                n += 1
                cid = f"{base}~{n}"
            self.registry[canon] = cid
            self.used_ids[cid] = canon
        return canon, self.registry[canon], ("0" if s == rc(s) else "+" if s == canon else "-")

    def step(self, st):
        s = self.nodes[st[:-1]]
        return s if st[-1] == "+" else rc(s)

    def keep(self, s):
        return len(s) >= ANCHOR_BP or min(s, rc(s)) in self.class_defining

    # ---- metadata status
    def allele_walks(self, b, a):
        """(status, walks): OK (recorded walks), VALID_EMPTY ([[]]), or MISSING."""
        src = self.sources.get((b, a))
        rec = self.walks.get((b, a), [])
        if src is None:
            return ("MISSING", []) if not rec else ("MISSING", rec)   # no source record: metadata absent
        if src["walk_sources"] > 0:
            return ("OK", rec) if len(rec) == src["walk_sources"] else ("MISSING", rec)
        if src["empty_sources"] + src["bypass_sources"] > 0 and not rec:
            return ("VALID_EMPTY", [[]])
        return ("MISSING", rec)

    # ---- boundary audit
    def flank(self, b, side):
        """The bytes the trim moved out of variant block b (identical for all its alleles), or None if inconsistent."""
        P, S = self.trim.get(b, (0, 0))
        got = set()
        for (bb, a), _ in self.seq.items():
            if bb != b:
                continue
            for w in self.allele_walks(b, a)[1]:
                sp = "".join(self.step(s) for s in w)
                got.add(sp[:P] if side == "prefix" else (sp[len(sp) - S:] if S else ""))
        return got.pop() if len(got) == 1 else None

    def audit(self):
        if self._audit is not None:
            return self._audit
        res = {}
        nb = len(self.block_ids)
        flanks = {b: (self.flank(b, "prefix"), self.flank(b, "suffix")) for b in range(nb) if self.kind[b] == "SV"}
        for b in range(nb):
            ok = bad = missing = 0
            for a, s in self.alleles[b]:
                status, ws = self.allele_walks(b, a)
                if status == "MISSING":
                    missing += 1
                    continue
                for w in ws:
                    sp = "".join(self.step(x) for x in w)
                    if self.kind[b] == "SV":
                        P, S = self.trim.get(b, (0, 0))
                        expect = sp[P:len(sp) - S]
                    else:
                        left = flanks.get(b - 1, (None, ""))[1] if b > 0 and self.kind[b - 1] == "SV" else ""
                        right = flanks.get(b + 1, ("", None))[0] if b + 1 < nb and self.kind[b + 1] == "SV" else ""
                        expect = None if left is None or right is None else left + sp + right
                    if expect == s:
                        ok += 1
                    else:
                        bad += 1
            res[b] = dict(block=self.block_ids[b], kind=self.kind[b], trimmed_prefix=self.trim.get(b, (0, 0))[0],
                          trimmed_suffix=self.trim.get(b, (0, 0))[1], walks_verified=ok, walks_failed=bad, alleles_missing=missing,
                          verified=bad == 0)   # missing metadata makes projections using that allele INCOMPLETE, not the block
        self._audit = res
        return res

    # ---- projection
    def block(self, block_id):
        return self.block_ids.index(block_id)

    def spell(self, f, alleles):
        return "".join(self.seq[(f + i, a)] for i, a in enumerate(alleles))

    def parses(self, H, f, l):
        memo = {}

        def go(b, o):
            if b > l:
                return [()] if o == len(H) else []
            if (b, o) not in memo:
                res = []
                for a, s in self.alleles[b]:
                    if H.startswith(s, o):
                        for rest in go(b + 1, o + len(s)):
                            res.append((a,) + rest)
                            if len(res) > MAX_PARSES:
                                raise OverflowError
                memo[(b, o)] = res
            return memo[(b, o)]
        try:
            return go(f, 0)
        except OverflowError:
            return "OVERFLOW"

    def scope(self, f, l):
        """(supported, clip_left, clip_right, reason) from the audit and the trims."""
        aud = self.audit()
        bad = [self.block_ids[b] for b in range(f, l + 1) if not aud[b]["verified"]]
        if bad:
            return False, 0, 0, "audit failed in " + ",".join(bad)
        nb = len(self.block_ids)
        if self.kind[f] != "SV" and f > 0 and self.kind[f - 1] == "SV" and self.trim.get(f - 1, (0, 0))[1] > 0:
            return False, 0, 0, f"{self.trim[f - 1][1]} bytes enter from {self.block_ids[f - 1]} (neighbour walk not determined)"
        if self.kind[l] != "SV" and l + 1 < nb and self.kind[l + 1] == "SV" and self.trim.get(l + 1, (0, 0))[0] > 0:
            return False, 0, 0, f"{self.trim[l + 1][0]} bytes enter from {self.block_ids[l + 1]} (neighbour walk not determined)"
        cl = self.trim.get(f, (0, 0))[0] if self.kind[f] == "SV" else 0
        cr = self.trim.get(l, (0, 0))[1] if self.kind[l] == "SV" else 0
        return True, cl, cr, ""

    def project(self, H, f, l):
        """(status, value, parses, note). status: RESOLVED | AMBIGUOUS | UNKNOWN | INCOMPLETE | UNSUPPORTED | OVERFLOW."""
        ok, cl, cr, why = self.scope(f, l)
        if not ok:
            return ("UNSUPPORTED", why, 0, "")
        ps = self.parses(H, f, l)
        if ps == "OVERFLOW":
            return ("OVERFLOW", f"more than {MAX_PARSES} parses", 0, "")
        if not ps:
            return ("UNKNOWN", "no exact parse into model alleles (off-model sequence)", 0, "")
        sks, empty_used = set(), False
        for p in ps:
            per_block = []
            for i, a in enumerate(p):
                status, ws = self.allele_walks(f + i, a)
                if status == "MISSING":
                    return ("INCOMPLETE", f"missing walk metadata for {self.block_ids[f + i]} allele {a}", len(ps), "")
                empty_used |= status == "VALID_EMPTY"
                per_block.append(ws)
            for combo in itertools.product(*per_block):
                W = [st for w in combo for st in w]
                total = sum(len(self.step(st)) for st in W)
                if total != len(H) + cl + cr:
                    return ("UNSUPPORTED", f"coordinate mismatch: walks spell {total}, expected {len(H) + cl + cr}", len(ps), "")
                pos, out = 0, []
                for st in W:
                    s = self.step(st)
                    start, end = pos, pos + len(s)
                    pos = end
                    inside = start >= cl and end <= total - cr
                    outside = end <= cl or start >= total - cr
                    if not self.keep(s):
                        continue
                    if not inside and not outside:
                        return ("UNSUPPORTED", f"a skeleton node ({len(s)} bp) straddles the scope clip", len(ps), "")
                    if inside:
                        canon, cid, orient = self.class_of(s)
                        out.append((canon, orient))
                sks.add(tuple(out))
                if len(sks) > MAX_PARSES:
                    return ("OVERFLOW", f"more than {MAX_PARSES} skeletons", len(ps), "")
        sks = frozenset(sks)
        return ("RESOLVED" if len(sks) == 1 else "AMBIGUOUS", sks, len(ps), "valid empty allele used" if empty_used else "")

    def display(self, skeleton):
        return ",".join(self.registry[c] + o for c, o in skeleton)


def profile(skeleton):
    c = collections.Counter(cls for cls, _ in skeleton)
    rep = sorted((n for n in c.values() if n > 1), reverse=True)
    return {"steps": len(skeleton), "max_multiplicity": max(c.values()) if c else 0,
            "multiplicity_signature": ",".join(map(str, collections.Counter(rep).most_common())) if rep else "none"}


def compare(pa, pb):
    if any(p[0] in FAILING for p in pa + pb):
        return "UNSUPPORTED", ";".join(f"{p[0]}:{p[1]}" for p in pa + pb if p[0] in FAILING)
    A, B = [p[1] for p in pa], [p[1] for p in pb]
    setv = any(len(x) > 1 for x in A + B)
    if collections.Counter(A) == collections.Counter(B):
        return ("SET_IDENTICAL" if setv else "UNIQUE_AGREEMENT"), ""
    for x, y in ((B[0], B[1]), (B[1], B[0])):
        if A[0] & x and A[1] & y:
            return "COMPATIBLE", ""

    def union(S):
        out = collections.Counter()
        for proj in S:
            out += collections.Counter(cls for cls, _ in next(iter(proj)))
        return out
    return "DISAGREEMENT", "DOSAGE_EQUAL_ALLOCATION_DIFFERS" if not setv and union(A) == union(B) else ""
