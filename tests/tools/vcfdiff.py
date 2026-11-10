#!/usr/bin/env python3
"""Field-by-field diff of two panvar call VCFs (region or allele VCF), ORIGINAL vs CURRENT.

Record IDs are not a safe join key: the _2/_3 suffixes are assigned in emission order and the
START_NODE part renumbers whenever the graph is rebuilt, so a plain ID diff reports a changed record
set as carrier flips. Records are therefore paired in five passes, each over what is still unpaired:

  1. content    CHROM, POS, BUBBLE_ID, SVTYPE, END, SVLEN, REF, ALT and the event sequence
                (INSSEQ/DELSEQ/INVSEQ), hashed; same-ID candidates are taken first
  2. carriers   (BUBBLE_ID, SVTYPE, SVLEN, carrier set): the same event with the same carriers is
                one record whatever its ID and position -- a representative that moved (a P-line
                order tie), or two records that swapped their _2/_3 suffixes while their position
                also changed. Before the ID pass, which would pair each with the other's old record
                and report the swap as carrier flips
  3. id         the record ID
  4. site       (BUBBLE_ID, SVTYPE, POS, SVLEN)
  5. nodes      (BUBBLE_ID, SVTYPE, sorted EVENT_NODES)

A pair whose IDs differ is reported as an ID rename, separately from any content change: a rename
alone is not a difference, a rename with a content change still has its CHANGED item. Every
INFO key, every FORMAT key and every per-haplotype FORMAT value is compared. CN-critical changes
come first: records added or removed, per-haplotype GT/CN/CNBP on a DUP, and the carrier set
(GT) of any other record.

stdout: a one-line summary, then the details. --items writes one machine-readable row per
difference (kind, orig_id, cur_id, svtype, fields, critical, detail, pin) for the expected-diff
check. pin is a hash of the complete difference (every changed value, every haplotype; the detail
text is shortened), so an expected-diff row can pin the exact change and not only its field names.
Exit 0 when there is no record-level difference, 1 otherwise (header lines and renames alone do not
count).
"""
import argparse
import collections
import hashlib
import sys

SEQ_KEYS = ('INSSEQ', 'DELSEQ', 'INVSEQ')
CN_KEYS = ('GT', 'CN', 'CNBP')


def h(s):
    return hashlib.sha1(str(s).encode()).hexdigest()[:12]


def parse(path):
    headers, samples, recs = [], [], []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('##'):
                headers.append(line)
                continue
            if line.startswith('#CHROM'):
                samples = line.split('\t')[9:]
                continue
            if not line:
                continue
            f = line.split('\t')
            info = collections.OrderedDict()
            if f[7] != '.':
                for kv in f[7].split(';'):
                    k, eq, v = kv.partition('=')
                    info[k] = v if eq else True
            fmt = f[8].split(':') if len(f) > 8 else []
            per = {}
            for s, col in zip(samples, f[9:]):
                vals = col.split(':')
                per[s] = {k: (vals[i] if i < len(vals) else '.') for i, k in enumerate(fmt)}
            recs.append(dict(chrom=f[0], pos=int(f[1]), id=f[2], ref=f[3], alt=f[4], qual=f[5],
                             filter=f[6], info=info, fmt=fmt, per=per))
    return headers, samples, recs


def svtype(r):
    if 'SVTYPE' in r['info']:
        return r['info']['SVTYPE']
    # the --multiallelic-loci record and the allele VCF's per-bubble record carry no SVTYPE
    return 'MULTI' if '_MULTI_' in r['id'] else 'ALLELES' if r['id'].endswith('_ALLELES') else '?'


def carriers(r):
    # any non-reference, non-missing allele index in GT
    out = set()
    for s, d in r['per'].items():
        gt = d.get('GT', '.')
        if any(a not in ('0', '.', '') for a in gt.replace('|', '/').split('/')):
            out.add(s)
    return out


def k_content(r):
    i = r['info']
    seq = next((i[k] for k in SEQ_KEYS if k in i), '')
    return (r['chrom'], r['pos'], i.get('BUBBLE_ID'), svtype(r), i.get('END'), i.get('SVLEN'),
            h(r['ref']), h(r['alt']), h(seq))


def k_id(r):
    return r['id']


def k_site(r):
    i = r['info']
    return (i.get('BUBBLE_ID'), svtype(r), r['pos'], i.get('SVLEN'))


def k_nodes(r):
    i = r['info']
    return (i.get('BUBBLE_ID'), svtype(r), tuple(sorted(str(i.get('EVENT_NODES', '')).split(','))))


def k_carriers(r):
    i = r['info']
    return (i.get('BUBBLE_ID'), svtype(r), i.get('SVLEN'), frozenset(carriers(r)))


def whole(r):
    """The record minus its ID, samples in name order: what an ADDED/REMOVED pin covers."""
    per = ';'.join(f'{s}={":".join(r["per"][s].get(k, ".") for k in r["fmt"])}' for s in sorted(r['per']))
    return '\t'.join(map(str, (r['chrom'], r['pos'], r['ref'], r['alt'], r['qual'], r['filter'],
                                ';'.join(f'{k}={v}' for k, v in r['info'].items()), ':'.join(r['fmt']), per)))


def pair(A, B):
    """Return [(ia, ib, pass_name)], unpaired A indices, unpaired B indices."""
    ua, ub = list(range(len(A))), set(range(len(B)))
    pairs = []
    for name, keyf in (('content', k_content), ('carriers', k_carriers), ('id', k_id), ('site', k_site),
                       ('nodes', k_nodes)):
        idx = collections.defaultdict(list)
        for j in sorted(ub):
            idx[keyf(B[j])].append(j)
        rest = []
        for i in ua:
            cand = idx.get(keyf(A[i]))
            if not cand:
                rest.append(i)
                continue
            # among equal keys, prefer the candidate that keeps the ID, then emission order
            same = [j for j in cand if B[j]['id'] == A[i]['id']]
            j = same[0] if same else cand[0]
            cand.remove(j)
            ub.discard(j)
            pairs.append((i, j, name))
        ua = rest
    return pairs, ua, sorted(ub)


def short(v, n=60):
    s = str(v)
    return s if len(s) <= n else f'{s[:24]}...({len(s)}ch,{h(s)[:6]})'


def hdr(line):
    # a structured header line is named by its ID; the rest is a description, hashed
    if '<ID=' in line:
        return line.split(',', 1)[0] + (f',...> ({h(line)[:6]})' if ',' in line else '')
    return short(line, 100)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('orig')
    ap.add_argument('cur')
    ap.add_argument('--label', default='')
    ap.add_argument('--items', help='write one TSV row per difference here')
    ap.add_argument('--max-lines', type=int, default=400)
    a = ap.parse_args()

    ha, sa, A = parse(a.orig)
    hb, sb, B = parse(a.cur)
    items = []   # (kind, orig_id, cur_id, svtype, fields, critical, detail, pin)

    if sa != sb:
        items.append(('SAMPLES', '-', '-', '-', 'samples', int(set(sa) != set(sb)),
                      f'orig {len(sa)} cur {len(sb)}; only orig {len(set(sa) - set(sb))}, '
                      f'only cur {len(set(sb) - set(sa))}'
                      + ('' if set(sa) != set(sb) else ' (same set, different order)'),
                      h(f'{sa}|{sb}')))
    # header lines: reported, never a failure (the expected-diff check ignores kind HEADER)
    ca, cb = collections.Counter(ha), collections.Counter(hb)
    for l in sorted((ca - cb).elements()):
        items.append(('HEADER', '-', '-', '-', 'header', 0, '- ' + hdr(l), h('- ' + l)))
    for l in sorted((cb - ca).elements()):
        items.append(('HEADER', '-', '-', '-', 'header', 0, '+ ' + hdr(l), h('+ ' + l)))

    pairs, ua, ub = pair(A, B)
    for i in ua:
        r = A[i]
        items.append(('REMOVED', r['id'], '-', svtype(r), '-', 1,
                      f'pos={r["pos"]} SVLEN={r["info"].get("SVLEN")} AC={r["info"].get("AC")}', h(whole(r))))
    for j in ub:
        r = B[j]
        items.append(('ADDED', '-', r['id'], svtype(r), '-', 1,
                      f'pos={r["pos"]} SVLEN={r["info"].get("SVLEN")} AC={r["info"].get("AC")}', h(whole(r))))

    changed, renamed, passes = 0, 0, collections.Counter()
    common = [s for s in sa if s in set(sb)]
    for i, j, how in pairs:
        ra, rb = A[i], B[j]
        passes[how] += 1
        st = svtype(ra)
        if ra['id'] != rb['id']:
            renamed += 1
            items.append(('RENAME', ra['id'], rb['id'], st, 'ID', 0, f'paired by {how}',
                          h(f'{ra["id"]}->{rb["id"]}')))
        # full: every changed value unshortened, for the pin; det is the readable, shortened form
        fields, det, full, crit = [], [], [], 0
        for col in ('chrom', 'pos', 'ref', 'alt', 'qual', 'filter'):
            if ra[col] != rb[col]:
                fields.append(col.upper())
                det.append(f'{col.upper()} {short(ra[col])}->{short(rb[col])}')
                full.append(f'{col.upper()} {ra[col]}->{rb[col]}')
        for k in dict.fromkeys(list(ra['info']) + list(rb['info'])):
            va, vb = ra['info'].get(k, '<absent>'), rb['info'].get(k, '<absent>')
            if va != vb:
                fields.append(f'INFO.{k}')
                det.append(f'INFO.{k} {short(va)}->{short(vb)}')
                full.append(f'INFO.{k} {va}->{vb}')
        if ra['fmt'] != rb['fmt']:
            fields.append('FORMAT')
            det.append(f'FORMAT {":".join(ra["fmt"])}->{":".join(rb["fmt"])}')
            full.append(det[-1])
        fdiff, cnd = collections.Counter(), collections.defaultdict(list)
        for s in common:
            pa, pb = ra['per'].get(s, {}), rb['per'].get(s, {})
            for k in dict.fromkeys(list(pa) + list(pb)):
                x, y = pa.get(k, '<absent>'), pb.get(k, '<absent>')
                if x != y:
                    fdiff[k] += 1
                    full.append(f'{s}:{k}:{x}->{y}')
                    if st == 'DUP' and k in CN_KEYS:
                        cnd[k].append(f'{s}:{x}->{y}')
        for k, n in fdiff.items():
            fields.append(f'FORMAT.{k}')
            det.append(f'FORMAT.{k} differs in {n} haplotypes')
        if st == 'DUP' and cnd:
            crit = 1
            for k, l in cnd.items():
                det.append(f'{k}: ' + ', '.join(l[:12]) + (f' ... +{len(l) - 12}' if len(l) > 12 else ''))
        elif st != 'DUP':
            c1, c2 = carriers(ra), carriers(rb)
            if c1 != c2:
                crit = 1
                det.append(f'carriers {len(c1)}->{len(c2)} (-{len(c1 - c2)} +{len(c2 - c1)})')
            elif fdiff.get('GT'):
                crit = 1   # same carriers, different allele index or phasing
        if fields:
            changed += 1
            items.append(('CHANGED', ra['id'], rb['id'], st, ','.join(fields), crit, '; '.join(det),
                          h('\n'.join(sorted(full)))))

    rec_items = [x for x in items if x[0] not in ('HEADER', 'RENAME')]
    ncrit = sum(1 for x in rec_items if x[5])
    nhdr = sum(1 for x in items if x[0] == 'HEADER')
    print(f'{a.label}: orig={len(A)} cur={len(B)} paired={len(pairs)} '
          f'(content={passes["content"]} id={passes["id"]} site={passes["site"]} nodes={passes["nodes"]} '
          f'carriers={passes["carriers"]}) '
          f'changed={changed} renamed={renamed} removed={len(ua)} added={len(ub)} '
          f'cn_critical={ncrit} header_lines={nhdr}')
    order = sorted(items, key=lambda x: (x[0] == 'HEADER', -x[5], x[0] == 'RENAME'))
    lines = [f'{"CRIT " if x[5] else ""}{x[0]} {x[3]} {x[1]}'
             + (f' -> {x[2]}' if x[2] not in ('-', x[1]) else '')
             + (f' [{x[4]}]' if x[0] == 'CHANGED' else '') + f': {x[6]}' for x in order]
    for l in lines[:a.max_lines]:
        print('  ' + l)
    if len(lines) > a.max_lines:
        print(f'  ... {len(lines) - a.max_lines} more')
    if a.items:
        with open(a.items, 'w') as fo:
            for x in order:
                fo.write('\t'.join(str(v).replace('\t', ' ') for v in x) + '\n')
    sys.exit(1 if rec_items else 0)


if __name__ == '__main__':
    main()
