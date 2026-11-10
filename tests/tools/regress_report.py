#!/usr/bin/env python3
"""Check the differences found by call_regress.sh against the expected-diff file and print the report.

  regress_report.py <items-dir> <expected.tsv>
  regress_report.py <items-dir> <expected.tsv> --emit-expected > new.tsv

<items-dir> holds one <set>__<locus>__<config>__<file>.items.tsv per comparison (written by
vcfdiff.py, bimbamdiff.py, nodesdiff.py, cnrecords.py, crossings.py and call_regress.sh itself) plus
a matching .summary line. Each item row is: kind, orig_id, cur_id, svtype, fields, critical, detail
and an optional pin, a hash of the complete difference (the detail text is shortened, so it is not
enough to pin a value); without one, the pin is the hash of the detail.

Expected-diff rows (tab-separated, '#' starts a comment; note is optional):

  status  locus  inputs  config  file  kind  id  fields  pin  note

status is 'approved' or 'pending'. locus/inputs/config/file/kind are globs ('*' = any; each may be a
comma list of globs). id is a glob matched against the original OR the current id. fields applies
to CHANGED items: every changed field must match one of its comma-separated globs, so a row that
approves a POS move does not also approve a GT change on the same record. pin is a comma list of
the pins this row accepts: a matching item whose pin is not listed changed VALUE since the row was
written (one more haplotype's GT, a further POS shift), and fails, pending or approved.

Every difference must be covered by a row. An uncovered one is a FAIL (exit 1). A 'pending' one is
listed under PENDING USER APPROVAL and does not fail. Header-line (##) changes and ID renames are
listed and counted but never fail on their own: a renamed record whose content also changed has a
CHANGED item that does. Rows that matched nothing are listed as a warning (they are stale, or the
run was a subset).

--emit-expected prints the expected-diff file back with every exercised row's pin list set to what
this run observed, and appends a 'pending' row (note TODO) for each uncovered difference. A row
whose pins changed is printed as 'pending' whatever its status was, and is named on stderr. Review
the diff against the committed file before replacing it; only the user moves a row to approved.
"""
import argparse
import collections
import fnmatch
import glob
import hashlib
import os
import re
import sys

NO_FAIL_KINDS = ('HEADER', 'RENAME')
COLUMNS = ('status', 'locus', 'inputs', 'config', 'file', 'kind', 'id', 'fields', 'pin', 'note')
PIN_RE = re.compile(r'^[0-9a-f]{12}$')
CONFIG_ORDER = ('default', 'cn', 'allele', 'classify', 'cnspacing', '-')


def load_expected(path):
    rows = []
    with open(path) as fh:
        for n, line in enumerate(fh, 1):
            line = line.rstrip('\n')
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            f = line.split('\t')
            if f[0] == 'status':
                continue
            if len(f) < 9 or f[0] not in ('approved', 'pending'):
                sys.exit(f'{path}:{n}: malformed expected-diff row (need 9 tab-separated columns, status '
                         f'approved|pending, then locus inputs config file kind id fields pin; a 10th, '
                         f'note, is optional): {line}')
            pins = [p.strip() for p in f[8].split(',')]
            bad = [p for p in pins if not PIN_RE.match(p)]
            if bad:
                sys.exit(f'{path}:{n}: malformed pin {bad[0]!r} (a comma list of 12-hex-digit pins, as '
                         f'--emit-expected prints them): {line}')
            rows.append(dict(n=n, status=f[0], locus=f[1], inputs=f[2], config=f[3], file=f[4], kind=f[5],
                             id=f[6], fields=f[7], pins=pins, note=f[9] if len(f) > 9 else ''))
    return rows


def gm(value, pattern):
    return any(fnmatch.fnmatchcase(value, p.strip()) for p in pattern.split(','))


def matches(row, it):
    """The row covers this kind of difference on this record (its value is checked by the pin)."""
    if not (gm(it['locus'], row['locus']) and gm(it['inputs'], row['inputs'])
            and gm(it['config'], row['config']) and gm(it['file'], row['file'])
            and gm(it['kind'], row['kind'])):
        return False
    if not (fnmatch.fnmatchcase(it['orig_id'], row['id']) or fnmatch.fnmatchcase(it['cur_id'], row['id'])):
        return False
    if it['kind'] == 'CHANGED':
        pats = [p.strip() for p in row['fields'].split(',')]
        return all(any(fnmatch.fnmatchcase(f, p) for p in pats) for f in it['fields'].split(','))
    return True


def brief(fil, line):
    """One token per comparison for the per-run line; the full summary stays in the run's files."""
    m = re.search(r'orig=(\d+) cur=(\d+)', line)
    if fil in ('region', 'alleles', 'bimbam', 'nodes') and m:
        return f'{fil} {m.group(1)}/{m.group(2)}'
    m = re.search(r'traversals=(\d+) .*paths_multi_crossing=(\d+)', line)
    if fil == 'crossing' and m:
        return f'traversals {m.group(1)}, multi-crossing paths={m.group(2)}'
    m = re.search(r'bubbles_multi_cn cur=(\d+)', line)
    if fil == 'cnrec' and m:
        return f'multi-CN bubbles={m.group(1)}'
    return f'{fil}: ' + line.replace('Error: call: ', '')[:110]


def load_items(d):
    items, summaries = [], collections.OrderedDict()
    for p in sorted(glob.glob(os.path.join(d, '*.items.tsv'))):
        base = os.path.basename(p)[:-len('.items.tsv')]
        st, loc, cfg, fil = base.split('__')
        with open(p) as fh:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) < 7:
                    continue
                pin = f[7] if len(f) > 7 and f[7] else hashlib.sha1(f[6].encode()).hexdigest()[:12]
                items.append(dict(inputs=st, locus=loc, config=cfg, file=fil, kind=f[0], orig_id=f[1],
                                  cur_id=f[2], svtype=f[3], fields=f[4], crit=f[5] == '1', detail=f[6], pin=pin))
    for p in sorted(glob.glob(os.path.join(d, '*.summary'))):
        base = os.path.basename(p)[:-len('.summary')]
        with open(p) as fh:
            summaries[base] = fh.read().strip()
    return items, summaries


def cfg_key(c):
    return (CONFIG_ORDER.index(c) if c in CONFIG_ORDER else len(CONFIG_ORDER), c)


def emit_expected(exp, rows, items):
    """Print the expected-diff file with observed pins; see the module docstring."""
    by_row, new = collections.defaultdict(list), collections.OrderedDict()
    for it in items:
        if it['kind'] in NO_FAIL_KINDS:
            continue
        row = next((r for r in rows if matches(r, it)), None)
        if row:
            by_row[row['n']].append(it)
        else:
            rid = it['cur_id'] if it['cur_id'] != '-' else it['orig_id']
            new.setdefault((it['locus'], it['inputs'], it['file'], it['kind'], rid), []).append(it)
    rown = {r['n']: r for r in rows}
    with open(exp) as fh:
        for n, line in enumerate(fh, 1):
            line = line.rstrip('\n')
            if line.startswith('status\t'):
                print('\t'.join(COLUMNS))
            elif n in rown and n in by_row:
                r = rown[n]
                pins = sorted({it['pin'] for it in by_row[n]})
                status = r['status']
                if not set(pins) <= set(r['pins']):
                    status = 'pending'
                    print(f'regress_report: row {n} ({r["note"][:60]}): pins {",".join(r["pins"])} -> '
                          f'{",".join(pins)}, printed as pending', file=sys.stderr)
                f = line.split('\t')
                f[0], f[8] = status, ','.join(pins)
                print('\t'.join(f))
            else:
                print(line)
    if new:
        print('# New differences (emitted by regress_report.py --emit-expected): replace each TODO with what')
        print('# the change is and why, and leave the row pending until the user approves it.')
        for (loc, st, fil, kind, rid), its in new.items():
            cfgs = sorted({it['config'] for it in its}, key=cfg_key)
            fields = list(dict.fromkeys(f for it in its for f in it['fields'].split(','))) \
                if kind == 'CHANGED' else ['-']
            pins = sorted({it['pin'] for it in its})
            print('\t'.join(('pending', loc, st, ','.join(cfgs), fil, kind, rid, ','.join(fields), ','.join(pins),
                             f'TODO {its[0]["detail"][:160]}')))
        print(f'regress_report: {len(new)} new row(s) appended as pending', file=sys.stderr)
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('items_dir')
    ap.add_argument('expected')
    ap.add_argument('--emit-expected', action='store_true',
                    help='print the expected-diff file with the pins and rows this run observed')
    a = ap.parse_args()
    rows = load_expected(a.expected)
    items, summaries = load_items(a.items_dir)
    if a.emit_expected:
        return emit_expected(a.expected, rows, items)

    used = collections.Counter()
    unexpected, pending, approved = [], collections.defaultdict(list), collections.defaultdict(list)
    headers, renames = collections.defaultdict(list), collections.defaultdict(list)
    per_run = collections.defaultdict(collections.Counter)
    for it in items:
        run = f'{it["inputs"]}/{it["locus"]}/{it["config"]}'
        if it['kind'] == 'HEADER':
            headers[it['detail']].append(f'{run}:{it["file"]}')
            per_run[run]['header'] += 1
            continue
        if it['kind'] == 'RENAME':
            renames[(it['file'], it['svtype'], it['orig_id'], it['cur_id'], it['detail'])].append(run)
            per_run[run]['rename'] += 1
            continue
        cover = [r for r in rows if matches(r, it)]
        row = next((r for r in cover if it['pin'] in r['pins']), None)
        if row is None:
            # a row that names this difference but pins another value is not a cover
            it['value_of'] = cover[0] if cover else None
            unexpected.append(it)
            per_run[run]['UNEXPECTED'] += 1
        else:
            used[row['n']] += 1
            (pending if row['status'] == 'pending' else approved)[row['n']].append(it)
            per_run[run][row['status']] += 1
    for it in unexpected:
        if it['value_of']:
            used[it['value_of']['n']] += 1

    def fmt(it):
        ids = it['orig_id'] if it['cur_id'] in ('-', it['orig_id']) else (
            it['cur_id'] if it['orig_id'] == '-' else f'{it["orig_id"]}->{it["cur_id"]}')
        fl = f' [{it["fields"]}]' if it['kind'] == 'CHANGED' else ''
        det = it['detail'] if len(it['detail']) <= 240 else it['detail'][:240] + ' ...'
        return f'{it["inputs"]}/{it["locus"]}/{it["config"]} {it["file"]} {it["kind"]} {it["svtype"]} {ids}{fl}: {det}'

    print('== per-run summary (inputs/locus/config: verdict, record counts orig/cur) ==')
    runs = collections.OrderedDict()
    for k, v in summaries.items():
        st, loc, cfg, fil = k.split('__')
        runs.setdefault(f'{st}/{loc}/{cfg}', []).append((fil, v))
    for run, parts in runs.items():
        c = per_run.get(run, collections.Counter())
        verdict = 'identical' if not (c['UNEXPECTED'] or c['pending'] or c['approved']) else \
            ', '.join(f'{n} {k}' for k, n in (('UNEXPECTED', c['UNEXPECTED']), ('pending', c['pending']),
                                                ('approved', c['approved'])) if n)
        extra = [f'+{c[k]} {w}' for k, w in (('rename', 'rename(s)'), ('header', 'header line change(s)')) if c[k]]
        if extra:
            verdict += f' ({", ".join(extra)})'
        print(f'  {run:<26} {verdict:<22} {"  ".join(brief(f, v) for f, v in sorted(parts))}')

    if headers:
        print('\n== HEADER LINE CHANGES (reported, not failures) ==')
        for l, where in headers.items():
            print(f'  {l}   [{len(where)} file(s), e.g. {where[0]}]')

    if renames:
        print(f'\n== ID RENAMES (reported, not failures; {sum(map(len, renames.values()))} item(s), '
              f'{len(renames)} distinct) ==')
        for (fil, svt, a_id, b_id, det), where in renames.items():
            print(f'  {fil} {svt} {a_id}->{b_id}: {det}   [{len(where)} run(s), e.g. {where[0]}]')

    if approved:
        print('\n== APPROVED DIFFERENCES ==')
        for n, its in sorted(approved.items()):
            r = next(r for r in rows if r['n'] == n)
            print(f'  row {n}: {r["note"]} -- {len(its)} item(s), pin(s) {",".join(sorted({i["pin"] for i in its}))}')

    if pending:
        print('\n== PENDING USER APPROVAL (not failures until approved or rejected) ==')
        for n, its in sorted(pending.items()):
            r = next(r for r in rows if r['n'] == n)
            print(f'  row {n}: {r["note"]}')
            # one line per distinct difference, with the runs it was seen in
            groups = collections.OrderedDict()
            for it in its:
                k = (it['file'], it['kind'], it['svtype'], it['orig_id'], it['cur_id'], it['fields'], it['detail'],
                     it['pin'])
                groups.setdefault(k, []).append(f'{it["inputs"]}/{it["locus"]}/{it["config"]}')
            for k, where in groups.items():
                it = dict(zip(('file', 'kind', 'svtype', 'orig_id', 'cur_id', 'fields', 'detail', 'pin'), k))
                it.update(inputs='', locus='', config='')
                print(f'      {fmt(it).lstrip("/ ")}')
                print(f'        pin {it["pin"]}, seen in {len(where)} run(s): {", ".join(sorted(where))}')

    stale = [r for r in rows if not used[r['n']]]
    if stale:
        print('\n== expected-diff rows that matched nothing (stale, or not exercised by this run) ==')
        for r in stale:
            print(f'  row {r["n"]} [{r["status"]}] {r["locus"]}/{r["inputs"]}/{r["config"]} {r["file"]} '
                  f'{r["kind"]} {r["id"]}: {r["note"]}')

    if unexpected:
        unexpected.sort(key=lambda it: (not it['crit'], it['value_of'] is None))
        nval = sum(1 for it in unexpected if it['value_of'])
        print(f'\n== UNEXPECTED DIFFERENCES: FAIL ({len(unexpected)}, CN-critical first; {nval} of them a '
              f'changed value on a listed difference) ==')
        for it in unexpected:
            print(f'  {"CRIT " if it["crit"] else ""}{fmt(it)}')
            r = it['value_of']
            if r:
                print(f'      VALUE DIFFERS from row {r["n"]} [{r["status"]}] ({r["note"][:80]}): '
                      f'pin {it["pin"]}, row accepts {",".join(r["pins"])}')
        print(f'\ncall_regress: FAIL, {len(unexpected)} unexpected difference(s); '
              f'{sum(map(len, pending.values()))} pending, {sum(map(len, approved.values()))} approved, '
              f'{sum(map(len, renames.values()))} rename(s)')
        return 1
    print(f'\ncall_regress: PASS, no unexpected difference; {sum(map(len, pending.values()))} item(s) '
          f'pending user approval, {sum(map(len, approved.values()))} approved, '
          f'{sum(map(len, renames.values()))} rename(s)')
    return 0


if __name__ == '__main__':
    sys.exit(main())
