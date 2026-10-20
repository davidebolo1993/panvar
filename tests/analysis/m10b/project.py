#!/usr/bin/env python3
"""M10-B offline structural projection prototype (diagnostic only; no caller integration).

Input: a model dump from m10b_walk_dump (<prefix>.block_alleles.tsv, .walks.tsv, .nodes.tsv, .panel.tsv, .summary.tsv).

pi(H) for a haplotype SEQUENCE H over consecutive blocks f..l (never for a labelling):
  1. PARSES: every tuple of model alleles (one per block f..l) whose concatenated sequences equal H exactly. Alleles are
     sequence-unique per block (the model canonicalises them), so the only freedom is where block boundaries fall.
  2. SKELETON of a parse: the concatenation of each chosen allele's oriented node walk (EVERY recorded walk of the allele:
     the graph can spell one sequence through different nodes), keeping only
       - steps whose node is at least ANCHOR_BP long (EXPERIMENTAL DIAGNOSTIC CONVENTION: the configured min_variant_bp,
         50, reused as is -- not an approved production semantic), and
       - steps whose node sequence class is declared class-defining (an explicit, separate rule; empty unless declared).
     A step is written as its SEQUENCE CLASS (sha256 of the canonical strand of the node sequence, first 10 hex) and its
     orientation relative to that canonical strand (+, -, or 0 for a palindrome). Node ids never appear, so two nodes with
     identical sequence are the same class.
  3. pi(H) = the SET of skeletons over all parses (and all walk choices). Identical sequences therefore have identical
     projections by construction. A set of size > 1 is AMBIGUOUS (the model places the same sequence in more than one
     structural way). No parse: UNKNOWN (no exact parse into model alleles; never equated with anything).
Comparison of two diploid genotypes {pi(a1), pi(a2)} vs {pi(b1), pi(b2)} (unordered; per homologue, so allocation counts):
  UNKNOWN     any side UNKNOWN
  EQUAL       equal as multisets of projections (EQUAL_SET_VALUED when a projection is a set of size > 1)
  COMPATIBLE  not equal, but some pairing has intersecting projections on both homologues
  DISJOINT    otherwise; flagged DOSAGE_EQUAL_ALLOCATION_DIFFERS when the multiset union of skeleton classes over both
              homologues is equal (dosage alone never makes EQUAL)
No pairwise tolerance and no clustering anywhere.
"""
import collections
import csv
import hashlib
import itertools

ANCHOR_BP = 50          # = configured min_variant_bp; experimental diagnostic convention only
MAX_PARSES = 100000     # computational guard: exceeding it is reported as OVERFLOW, never truncated silently
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def T(p):
    return list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))


class Model:
    def __init__(self, prefix, class_defining=()):
        self.prefix = prefix
        self.nodes = {r["node"]: r["sequence"] for r in T(prefix + ".nodes.tsv")}
        self.seq = {}
        self.block_ids = []
        for r in T(prefix + ".block_alleles.tsv"):
            b = int(r["block_index"])
            self.seq[(b, int(r["allele_index"]))] = r["sequence"]
            if len(self.block_ids) <= b:
                self.block_ids.append(r["block_id"])
        self.walks = collections.defaultdict(list)
        for r in T(prefix + ".walks.tsv"):
            self.walks[(int(r["block_index"]), int(r["allele_index"]))].append(r["walk"].split(",") if r["walk"] else [])
        self.alleles = collections.defaultdict(list)
        for (b, a), s in sorted(self.seq.items()):
            self.alleles[b].append((a, s))
        self.class_defining = set(class_defining)   # sequence classes declared class-defining (explicit rule)
        self.panel = {r["template"]: [int(x) for x in r["alleles"].split(",")] for r in T(prefix + ".panel.tsv")}
        self.fingerprint = {r["key"]: r["value"] for r in T(prefix + ".summary.tsv")}.get("model_fingerprint")
        self._cls = {}

    def block(self, block_id):
        return self.block_ids.index(block_id)

    def step_class(self, step):
        """(class, orientation) of an oriented step, or None if it is not a skeleton step."""
        if step not in self._cls:
            s = self.nodes[step[:-1]]
            s = s if step[-1] == "+" else rc(s)
            canon = min(s, rc(s))
            cid = "U" + hashlib.sha256(canon.encode()).hexdigest()[:10]
            orient = "0" if s == rc(s) else ("+" if s == canon else "-")
            keep = len(s) >= ANCHOR_BP or cid in self.class_defining
            self._cls[step] = (cid, orient, len(s)) if keep else None
        return self._cls[step]

    def block_skeletons(self, b, a):
        out = set()
        for w in self.walks[(b, a)]:
            out.add(tuple((c[0], c[1]) for c in (self.step_class(s) for s in w) if c))
        return out or {("NO_WALK",)}

    def spell(self, f, alleles):
        return "".join(self.seq[(f + i, a)] for i, a in enumerate(alleles))

    def parses(self, H, f, l):
        """Every allele tuple over blocks f..l spelling H exactly (list of tuples), or 'OVERFLOW'."""
        memo = {}

        def go(b, o):
            if b > l:
                return [()] if o == len(H) else []
            key = (b, o)
            if key not in memo:
                res = []
                for a, s in self.alleles[b]:
                    if H.startswith(s, o):
                        for rest in go(b + 1, o + len(s)):
                            res.append((a,) + rest)
                            if len(res) > MAX_PARSES:
                                raise OverflowError
                memo[key] = res
            return memo[key]
        try:
            return go(f, 0)
        except OverflowError:
            return "OVERFLOW"

    def project(self, H, f, l):
        """('RESOLVED'|'AMBIGUOUS', frozenset of skeletons, parses) or ('UNKNOWN'|'OVERFLOW', reason, 0)."""
        ps = self.parses(H, f, l)
        if ps == "OVERFLOW":
            return ("OVERFLOW", f"more than {MAX_PARSES} parses", 0)
        if not ps:
            return ("UNKNOWN", "no exact parse into model alleles (off-model sequence)", 0)
        sks = set()
        for p in ps:
            per_block = [self.block_skeletons(f + i, a) for i, a in enumerate(p)]
            for combo in itertools.product(*per_block):
                sks.add(tuple(step for part in combo for step in part))
                if len(sks) > MAX_PARSES:
                    return ("OVERFLOW", f"more than {MAX_PARSES} skeletons", len(ps))
        sks = frozenset(sks)
        return ("RESOLVED" if len(sks) == 1 else "AMBIGUOUS", sks, len(ps))


def profile(skeleton):
    """Descriptive copy profile of one skeleton: steps, repeated classes and their multiplicities."""
    c = collections.Counter(cls for cls, _ in skeleton)
    rep = sorted((n for n in c.values() if n > 1), reverse=True)
    return {"steps": len(skeleton), "classes": len(c), "repeated_classes": len(rep), "max_multiplicity": max(c.values()) if c else 0,
            "multiplicity_signature": ",".join(map(str, collections.Counter(rep).most_common())) if rep else "none"}


def compare(pa, pb):
    """pa, pb: pairs of projections (status, value, parses) for the two homologues. Returns (verdict, detail)."""
    if any(p[0] in ("UNKNOWN", "OVERFLOW") for p in pa + pb):
        return "UNKNOWN", ";".join(f"{p[0]}:{p[1]}" for p in pa + pb if p[0] in ("UNKNOWN", "OVERFLOW"))
    A, B = [p[1] for p in pa], [p[1] for p in pb]
    setv = any(len(x) > 1 for x in A + B)
    if collections.Counter(A) == collections.Counter(B):
        return ("EQUAL_SET_VALUED" if setv else "EQUAL"), ""
    for x, y in ((B[0], B[1]), (B[1], B[0])):
        if A[0] & x and A[1] & y:
            return "COMPATIBLE", "some pairing has intersecting projections on both homologues"

    def union(S):
        out = collections.Counter()
        for proj in S:
            out += collections.Counter(cls for cls, _ in min(proj))   # canonical member for the dosage diagnostic
        return out
    flag = "DOSAGE_EQUAL_ALLOCATION_DIFFERS" if not setv and union(A) == union(B) else ""
    return "DISJOINT", flag
