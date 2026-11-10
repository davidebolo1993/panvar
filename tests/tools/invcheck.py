#!/usr/bin/env python3
"""Record-level assertions for tests/suites/call_invariants.sh, on top of vcfdiff.py's parser and pairing.

  invcheck.py same A B --fields F[,F...] [--only COND] [--exclude COND]
      Pair A's records with B's (vcfdiff.py's pairing passes) and compare the listed fields. A field is
      POS, REF, ALT, ID, an INFO key (END, INSSEQ, ...), GT (per-sample, by sample name) or ALL (every
      column, INFO as a map with EVENT_NODES as a set). Records of A selected by --only/--exclude
      must all pair and agree; selected records of B left unpaired also count. Exit 1 on any mismatch,
      and when no record of either is selected (nothing was compared), unless --allow-empty.
  invcheck.py carriers A B [--allow-empty]
      The multiset of (SVTYPE, SVLEN, carrier set) must be equal, whatever the record ids and POS. Exit 1
      when both hold no record, unless --allow-empty.
  invcheck.py count VCF [COND]
      Print the number of records matching COND.
  invcheck.py get VCF COND FIELD
      Print FIELD of the single record matching COND; FIELD is as for 'same', or FMT:<key>:<sample>.

COND is a comma-separated conjunction of KEY=GLOB (KEY is ID, SVTYPE or an INFO key), KEY (INFO key
present) and !KEY (absent).
"""
import argparse
import collections
import fnmatch
import os
import sys

sys.dont_write_bytecode = True   # no __pycache__ in the source tree
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import vcfdiff  # noqa: E402


def load(path):
    _, samples, recs = vcfdiff.parse(path)
    return samples, recs


def cond_ok(r, cond):
    if not cond:
        return True
    for t in cond.split(','):
        t = t.strip()
        if t.startswith('!'):
            if t[1:] in r['info']:
                return False
            continue
        k, eq, v = t.partition('=')
        if not eq:
            if k not in r['info']:
                return False
            continue
        val = r['id'] if k == 'ID' else (vcfdiff.svtype(r) if k == 'SVTYPE' else r['info'].get(k))
        if val is None or not fnmatch.fnmatchcase(str(val), v):
            return False
    return True


def field(r, f):
    if f == 'POS':
        return r['pos']
    if f in ('REF', 'ALT', 'ID'):
        return r[f.lower()]
    if f == 'GT':
        return tuple(sorted((s, d.get('GT')) for s, d in r['per'].items()))
    if f == 'ALL':
        info = dict(r['info'])
        if 'EVENT_NODES' in info:
            info['EVENT_NODES'] = tuple(sorted(str(info['EVENT_NODES']).split(',')))
        per = tuple(sorted((s, tuple(sorted(d.items()))) for s, d in r['per'].items()))
        return (r['chrom'], r['pos'], r['id'], r['ref'], r['alt'], r['qual'], r['filter'],
                tuple(sorted(info.items())), tuple(r['fmt']), per)
    if f.startswith('FMT:'):
        _, k, s = f.split(':', 2)
        return r['per'].get(s, {}).get(k)
    return r['info'].get(f, '<absent>')


def short(v):
    return vcfdiff.short(v, 70)


def cmd_same(a):
    _, A = load(a.a)
    _, B = load(a.b)
    sel = lambda r: cond_ok(r, a.only) and not (a.exclude and cond_ok(r, a.exclude))
    if not a.allow_empty and not any(sel(r) for r in A) and not any(sel(r) for r in B):
        print('    no record selected in either VCF: nothing compared (--allow-empty if that is expected)')
        return 1
    pairs, ua, ub = vcfdiff.pair(A, B)
    bad = 0
    fields = a.fields.split(',')
    for i, j, how in pairs:
        if not sel(A[i]):
            continue
        for f in fields:
            x, y = field(A[i], f), field(B[j], f)
            if x != y:
                bad += 1
                if f == 'ALL':
                    diff = [k for k in ('POS', 'REF', 'ALT', 'ID') if field(A[i], k) != field(B[j], k)]
                    diff += [f'INFO.{k}' for k in set(A[i]['info']) | set(B[j]['info'])
                             if A[i]['info'].get(k) != B[j]['info'].get(k)]
                    print(f'    {A[i]["id"]} vs {B[j]["id"]}: differs in {",".join(sorted(diff)) or "FORMAT"}')
                else:
                    print(f'    {A[i]["id"]} vs {B[j]["id"]}: {f} {short(x)} != {short(y)}')
    for i in ua:
        if sel(A[i]):
            bad += 1
            print(f'    {A[i]["id"]}: no counterpart in {os.path.basename(a.b)}')
    for j in ub:
        if sel(B[j]):
            bad += 1
            print(f'    {B[j]["id"]}: no counterpart in {os.path.basename(a.a)}')
    return 1 if bad else 0


def cmd_carriers(a):
    def ms(path):
        _, R = load(path)
        return collections.Counter((vcfdiff.svtype(r), r['info'].get('SVLEN'), frozenset(vcfdiff.carriers(r)))
                                   for r in R)
    x, y = ms(a.a), ms(a.b)
    if not a.allow_empty and not x and not y:
        print('    no record in either VCF: nothing compared (--allow-empty if that is expected)')
        return 1
    for k in (x - y):
        print(f'    only in A: {k[0]} SVLEN={k[1]} carriers={sorted(k[2])}')
    for k in (y - x):
        print(f'    only in B: {k[0]} SVLEN={k[1]} carriers={sorted(k[2])}')
    return 1 if x != y else 0


def cmd_count(a):
    _, R = load(a.vcf)
    print(sum(1 for r in R if cond_ok(r, a.cond)))
    return 0


def cmd_get(a):
    _, R = load(a.vcf)
    m = [r for r in R if cond_ok(r, a.cond)]
    if len(m) != 1:
        print(f'<{len(m)} records match>')
        return 1
    print(field(m[0], a.field))
    return 0


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    p = sp.add_parser('same'); p.add_argument('a'); p.add_argument('b')
    p.add_argument('--fields', required=True); p.add_argument('--only'); p.add_argument('--exclude')
    p.add_argument('--allow-empty', action='store_true')
    p = sp.add_parser('carriers'); p.add_argument('a'); p.add_argument('b')
    p.add_argument('--allow-empty', action='store_true')
    p = sp.add_parser('count'); p.add_argument('vcf'); p.add_argument('cond', nargs='?', default='')
    p = sp.add_parser('get'); p.add_argument('vcf'); p.add_argument('cond'); p.add_argument('field')
    a = ap.parse_args()
    return {'same': cmd_same, 'carriers': cmd_carriers, 'count': cmd_count,
            'get': cmd_get}[a.cmd](a)


if __name__ == '__main__':
    sys.exit(main())
