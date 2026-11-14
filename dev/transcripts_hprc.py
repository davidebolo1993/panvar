#!/usr/bin/env python3
"""`panvar transcripts` against HPRC's own annotations of the same assemblies (CAT and Liftoff, release 2):
copies per gene and haplotype, and each annotated transcript against our HST of it.

  transcripts_hprc.py select <dir> <locus> <prefix>.hst.tsv <cat.csv> <liftoff.csv> <gene> <copies=haplotypes,...> [extra gene ...]
      haplotypes with a CAT annotation, drawn by their copies of <gene> (e.g. GSTM1 0=8,1=9,2=2), and two
      with other than one copy of each extra gene; appended to <dir>/selected.tsv
  transcripts_hprc.py fetch <dir>
      each selected haplotype's CAT and Liftoff records over its path, into <dir>/cat and <dir>/lo
  transcripts_hprc.py compare <dir> <locus> <prefix>.hst.tsv <graph.gfa.gz> <genes.gtf.gz>
      <dir>/<locus>.pairs.tsv and .copies.tsv, then the summary on stdout

cat.csv and liftoff.csv are HPRC's release-2 annotation index files (sample_id, haplotype, assembly_name,
location). An annotated transcript is matched to our HST by its Ensembl id (CAT), MANE tag (Liftoff), same
exons, or same introns; the summary keeps each source transcript's best alignment, and classes the
disagreements (transcript version, novel isoforms, copies on another paralog, our statuses).
"""
import collections
import csv
import gzip
import os
import random
import re
import subprocess
import sys
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from transcripts_common import edit_distances, rc, walks  # noqa: E402,F401

HTTPS = 'https://s3-us-west-2.amazonaws.com/human-pangenomics'


def select(D, locus, hst_tsv, cat_csv, lo_csv, gene, plan, extra):
    cat = {(r['sample_id'], r['haplotype']): r['location'] for r in csv.DictReader(open(cat_csv))}
    lo = {(r['sample_id'], r['haplotype']): r['location'] for r in csv.DictReader(open(lo_csv))}
    d = collections.defaultdict(dict)
    with open(hst_tsv) as f:
        h = f.readline().rstrip('\n').split('\t')
        for l in f:
            t = dict(zip(h, l.rstrip('\n').split('\t')))
            if t['gene_name'] in [gene] + extra:
                d[t['haplotype']][t['gene_name']] = int(t['copies'])
    ok = [hp for hp in d if tuple(hp.split('#')[:2]) in cat]
    random.seed(7)
    picked = []
    for k, n in plan.items():
        pool = sorted(x for x in ok if d[x].get(gene, 0) == k)
        random.shuffle(pool)
        picked += pool[:n]
    for e in extra:
        picked += sorted(x for x in ok if d[x].get(e, 1) != 1 and x not in picked)[:2]
    with open(D + '/selected.tsv', 'a') as out:
        for hp in picked:
            s, h = hp.split('#')[:2]
            out.write('\t'.join([locus, hp, s, h, cat[(s, h)], lo.get((s, h), '.'),
                                 ','.join('%s=%d' % kv for kv in sorted(d[hp].items()))]) + '\n')


def fetch_cat(p, u, out):
    """A CAT GFF3's lines over path p's interval (the whole file is streamed)."""
    c = p.rsplit('#', 1)[1]
    c, iv = c.rsplit(':', 1)
    s, e = iv.split('-')
    cmd = ("curl -sf '%s' | gzip -dc | LC_ALL=C grep -F '%s\t' | awk -F'\\t' -v c='%s' -v s=%s -v e=%s "
           "'$1==c && $5>s && $4<=e'" % (u.replace('s3://human-pangenomics', HTTPS), c, c, s, e))
    with open(out + '.tmp', 'w') as o:
        subprocess.run(cmd, shell=True, stdout=o, check=True)
    os.rename(out + '.tmp', out)


def compare(D, locus, hst_tsv, gfa, gtf):
    """<locus>.pairs.tsv (one row per annotated transcript, or our HST without a partner) and
    <locus>.copies.tsv (copies per gene and haplotype on each side) in D; returns the pairs to edit-compare."""
    sel = [l.rstrip('\n').split('\t') for l in open(D + '/selected.tsv') if l.startswith(locus + '\t')]
    haps = [r[1] for r in sel]
    want = set(haps)
    walk = walks(gfa, want)

    # reference annotation: gene id <-> name, transcript -> gene, ref exon lengths
    tv = {}; mane = set(); tx_gene = {}; gene_name = {}; ref_ex = collections.defaultdict(list); tx_strand = {}
    with gzip.open(gtf, 'rt') as f:
        for l in f:
            if l[0] == '#': continue
            t = l.split('\t')
            if t[2] != 'exon': continue
            a = dict(re.findall(r'(\S+) "([^"]*)"', t[8]))
            tx_gene[a['transcript_id']] = a['gene_id']
            gene_name[a['gene_id']] = a.get('gene_name', a['gene_id'])
            ref_ex[a['transcript_id']].append((int(t[3]), int(t[4])))
            tx_strand[a['transcript_id']] = t[6]
            tv[a['transcript_id']] = a.get('transcript_version', '?')
            if 'tag "MANE_Select"' in t[8]: mane.add(a['transcript_id'])
    name_gene = {v: k for k, v in gene_name.items()}

    # our HSTs
    ours = collections.defaultdict(list)   # hap -> rows
    ours_copies = collections.defaultdict(dict)  # hap -> gene_name -> copies
    ours_genes = set()
    with open(hst_tsv) as f:
        h = f.readline().rstrip('\n').split('\t')
        for l in f:
            t = l.rstrip('\n').split('\t', 7)
            if t[6] not in want: continue
            t = dict(zip(h, l.rstrip('\n').split('\t')))
            ours_genes.add(t['gene_name'])
            ours_copies[t['haplotype']][t['gene_name']] = int(t['copies']) if t['copy'] != '0' else 0
            if t['copy'] == '0' or t['exon_starts'] == '.':
                ours[t['haplotype']].append(t); continue
            st = list(map(int, t['exon_starts'].split(','))); ln = list(map(int, t['exon_lengths'].split(',')))
            t['ex'] = sorted((s, s + n) for s, n in zip(st, ln))
            ours[t['haplotype']].append(t)

    def spliced(w, ex, strand):
        s = ''.join(w[a:b] for a, b in sorted(ex))
        return rc(s) if strand == '-' else s

    def parse_gff(fn, s0, e0, kind):
        tx = {}; genes = {}
        for l in open(fn):
            t = l.rstrip('\n').split('\t')
            if len(t) < 9: continue
            a = dict(x.split('=', 1) for x in t[8].split(';') if '=' in x)
            if t[2] == 'gene':
                genes[a['ID']] = dict(a, start=int(t[3]), end=int(t[4]), strand=t[6])
            elif t[2] in ('transcript', 'mRNA', 'ncRNA', 'lnc_RNA', 'snRNA', 'miRNA', 'primary_transcript', 'snoRNA', 'rRNA', 'tRNA', 'V_gene_segment', 'C_gene_segment'):
                tx[a['ID']] = dict(a, start=int(t[3]), end=int(t[4]), strand=t[6], ex=[], ftype=t[2])
            elif t[2] == 'exon':
                for p in a['Parent'].split(','):
                    tx.setdefault(p, {'ex': [], 'orphan': True})['ex'].append((int(t[3]) - 1 - s0, int(t[4]) - s0))
        out = []
        for k, v in tx.items():
            if 'start' not in v or not v['ex']: continue
            v['ex'] = sorted(set(v['ex']))
            v['inside'] = v['ex'][0][0] >= 0 and v['ex'][-1][1] <= e0 - s0
            v['tid'] = k
            out.append(v)
        return out, genes

    def ivl_ov(a, b): return min(a[1], b[1]) - max(a[0], b[0])

    pairs = open('%s/%s.pairs.tsv' % (D, locus), 'w')
    pairs.write('\t'.join(['annot', 'haplotype', 'side', 'gene', 'annot_gene', 'transcript', 'annot_tx', 'annot_class', 'our_copy', 'our_status',
                           'match', 'same_exons', 'same_introns', 'annot_len', 'our_len', 'ed_key', 'note']) + '\n')
    cop = open('%s/%s.copies.tsv' % (D, locus), 'w')
    cop.write('annot\thaplotype\tgene\tours\tannot_copies\tannot_on_other_gene\tnote\n')
    edin = []
    key = 0

    for r in sel:
        hp = r[1]; s0, e0 = map(int, hp.rsplit(':', 1)[1].split('-'))
        w = walk[hp]
        assert len(w) == e0 - s0, (hp, len(w), e0 - s0)
        rows = [t for t in ours[hp] if 'ex' in t]
        # our copy intervals per gene
        ourcopy = collections.defaultdict(list)
        for t in rows:
            ourcopy[(t['gene_name'], t['copy'])].append(t)
        cint = {k: (min(t['ex'][0][0] for t in v), max(t['ex'][-1][1] for t in v)) for k, v in ourcopy.items()}
        for annot in ('cat', 'lo'):
            fn = '%s/%s/%s.%s.gff3' % (D, annot, hp.replace('#', '_').replace(':', '_'), annot)
            if not os.path.exists(fn): continue
            txs, genes = parse_gff(fn, s0, e0, annot)
            # gene label of each annotated transcript
            for v in txs:
                if annot == 'cat':
                    st = v.get('source_transcript', 'N/A').split('.')[0]
                    v['enst'] = st if st.startswith('ENST') else None
                    gid = v.get('source_gene', '').split('.')[0]
                    v['gname'] = gene_name.get(gid, v.get('source_gene_common_name', v.get('gene_name')))
                    v['cls'] = v.get('transcript_class', '.') + '/' + v.get('alignment_mode', '.')
                    v['gkey'] = v['Parent']
                else:
                    v['enst'] = None
                    v['gname'] = v.get('gene', v.get('gene_name', '.'))

                    v['cls'] = ('MANE/' if 'MANE Select' in v.get('tag', '') else '') + v.get('ID', '.').split('_')[0][:2] + ('/copy' if 'extra_copy_number' in genes.get(v.get('Parent', ''), {}) and genes[v['Parent']]['extra_copy_number'] != '0' else '')
                    v['gkey'] = v.get('Parent', '.')
            loc_tx = [v for v in txs if v['gname'] in ours_genes and v['inside']]
            # copies per gene: annotated gene records holding at least one inside transcript
            agc = collections.defaultdict(set)
            for v in loc_tx: agc[v['gname']].add(v['gkey'])
            for g in sorted(ours_genes):
                n_our = ours_copies[hp].get(g, 0)
                n_an = len(agc.get(g, ()))
                other = []
                for gk in sorted(agc.get(g, ())):
                    vs = [v for v in loc_tx if v['gkey'] == gk]
                    iv = (min(v['ex'][0][0] for v in vs), max(v['ex'][-1][1] for v in vs))
                    hit = [k for k, c in cint.items() if k[0] == g and ivl_ov(c, iv) > 0]
                    if not hit:
                        oth = sorted(set(k[0] for k, c in cint.items() if ivl_ov(c, iv) > 0.3 * (iv[1] - iv[0])))
                        other.append('%d-%d:%s' % (iv[0], iv[1], '/'.join(oth) or 'none'))
                cop.write('\t'.join([annot, hp, g, str(n_our), str(n_an), str(len(other)), ';'.join(other) or '.']) + '\n')
            # match annotated transcripts to our HSTs
            used = set()
            for v in loc_tx:
                ex = v['ex']; iv = (ex[0][0], ex[-1][1])
                strand = v['strand']
                aseq = spliced(w, ex, strand)
                cands = [t for t in rows if t['gene_name'] == v['gname'] and ivl_ov((t['ex'][0][0], t['ex'][-1][1]), iv) > 0]
                m = None; how = 'none'
                if v['enst']:
                    c = [t for t in cands if t['transcript_id'] == v['enst']]
                    if c: m = c[0]; how = 'id'
                if m is None and annot == 'lo' and 'MANE Select' in v.get('tag', ''):
                    c = [t for t in cands if t['transcript_id'] in mane]
                    if c: m = c[0]; how = 'mane'
                if m is None:
                    c = [t for t in cands if t['ex'] == ex]
                    if c: m = c[0]; how = 'exons'
                if m is None:
                    introns = [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]
                    c = [] if v['enst'] else [t for t in cands if len(introns) > 0 and [(t['ex'][i][1], t['ex'][i + 1][0]) for i in range(len(t['ex']) - 1)] == introns]
                    if c: m = c[0]; how = 'introns'
                note = []
                if m is None:
                    # any HST of any gene on this stretch with the same intron chain or exons?
                    oth = [t for t in rows if t['gene_name'] != v['gname'] and t['ex'] == ex]
                    if oth: note.append('same_exons_as:%s:%s' % (oth[0]['gene_name'], oth[0]['transcript_id']))
                    ovg = sorted(set(t['gene_name'] for t in rows if ivl_ov((t['ex'][0][0], t['ex'][-1][1]), iv) > 0.5 * (iv[1] - iv[0])))
                    note.append('overlaps:' + ('/'.join(ovg) or 'none'))
                    if v['enst']:
                        # do we have this transcript at all on this haplotype (other copy) or as a non-written row?
                        allr = [t for t in ours[hp] if t['transcript_id'] == v['enst']]
                        note.append('ours_rows:' + ('/'.join('%s:%s' % (t['copy'], t['status']) for t in allr) or 'none'))
                        nc = [t['detail'] for t in allr if t['copy'] == '0']
                        if nc: note.append('ours_detail:' + nc[0].replace(';', ','))
                        if v['enst'] not in tx_gene: note.append('not_in_gtf')
                        elif gene_name[tx_gene[v['enst']]] != v['gname']: note.append('gtf_gene:' + gene_name[tx_gene[v['enst']]])
                    pairs.write('\t'.join([annot, hp, 'annot_only', v['gname'], v['gkey'], v.get('enst') or '.', v['tid'], v['cls'], '.', '.', how,
                                           '.', '.', str(len(aseq)), '.', '.', ';'.join(note)]) + '\n')
                    continue
                used.add(id(m))
                oseq = spliced(w, m['ex'], strand)
                same_ex = m['ex'] == ex
                mi = [(m['ex'][i][1], m['ex'][i + 1][0]) for i in range(len(m['ex']) - 1)]
                ai = [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]
                k = '.'
                if aseq != oseq:
                    key += 1; k = 'k%d' % key
                    edin.append((k, aseq, oseq))
                if not same_ex:
                    # describe exon differences
                    d = []
                    if len(m['ex']) != len(ex): d.append('nexons:%d/%d' % (len(ex), len(m['ex'])))
                    else:
                        for i, (a, b) in enumerate(zip(ex, m['ex'])):
                            if a != b: d.append('x%d:%+d,%+d' % (i, b[0] - a[0], b[1] - a[1]))
                    note.append('diff(ours-annot):' + ','.join(d[:6]))
                if v['enst']:
                    sv = v.get('source_transcript', '').split('.')
                    note.append('ver:%s/%s' % (sv[1] if len(sv) > 1 else '?', tv.get(v['enst'], '?')))
                if m['transcript_id'] != (v['enst'] or m['transcript_id']): note.append('ours_tx:' + m['transcript_id'])
                if annot == 'lo' or how != 'id': note.append('ours_tx:' + m['transcript_id'])
                pairs.write('\t'.join([annot, hp, 'both', v['gname'], v['gkey'], v.get('enst') or '.', v['tid'], v['cls'], m['copy'], m['status'], how,
                                       'yes' if same_ex else 'no', 'yes' if mi == ai else 'no', str(len(aseq)), str(len(oseq)), k, ';'.join(note) or '.']) + '\n')
            if annot == 'cat':
                # our HSTs with no CAT partner
                for t in ours[hp]:
                    if 'ex' not in t or id(t) in used: continue
                    cat_has = any(v['enst'] == t['transcript_id'] for v in txs)
                    cat_in = any(v['enst'] == t['transcript_id'] for v in loc_tx)
                    st = t['status'] if t['copy'] != '0' else 'no_copy'
                    pairs.write('\t'.join(['cat', hp, 'ours_only', t['gene_name'], '.', t['transcript_id'], '.', '.', t['copy'], st, 'none',
                                           '.', '.', '.', t.get('length', '.'), '.', 'cat_has_tx_elsewhere' if cat_has else ('cat_lacks_tx' if not cat_in else 'cat_tx_other_copy')]) + '\n')
    pairs.close(); cop.close()
    return edin


def summary(D, locus, ed):
    """Copies and transcripts compared, per annotation, from the tables compare() wrote."""
    rows = list(csv.DictReader(open('%s/%s.pairs.tsv' % (D, locus)), delimiter='\t'))
    cops = list(csv.DictReader(open('%s/%s.copies.tsv' % (D, locus)), delimiter='\t'))
    short = lambda h: h.split(':')[0].rsplit('#', 1)[0]

    for annot in ('cat', 'lo'):
        R = [r for r in rows if r['annot'] == annot]
        haps = sorted(set(r['haplotype'] for r in R if r['side'] != 'ours_only') | set(c['haplotype'] for c in cops if c['annot'] == annot))
        if not haps: continue
        print('==== %s %s: %d haplotypes' % (locus, annot, len(haps)))
        # copies
        C = [c for c in cops if c['annot'] == annot and c['haplotype'] in haps]
        per = collections.defaultdict(lambda: [0, 0, 0, []])
        for c in C:
            o, a, oth = int(c['ours']), int(c['annot_copies']), int(c['annot_on_other_gene'])
            if o == 0 and a == 0: continue
            p = per[c['gene']]
            p[0] += 1
            if o == a: p[1] += 1
            if o == a - oth: p[2] += 1
            if o != a: p[3].append('%s:%d/%d%s' % (short(c['haplotype']), o, a, ('[' + c['note'] + ']') if c['note'] != '.' else ''))
        print('gene\thaps\tagree_raw\tagree_excl_overlaid\tdisagreements(ours/annot)')
        for g, p in per.items():
            print('%s\t%d\t%d\t%d\t%s' % (g, p[0], p[1], p[2], ' '.join(p[3][:4]) + (' ...+%d' % (len(p[3]) - 4) if len(p[3]) > 4 else '')))
        # transcripts: group annotated records by (hap, gene record, enst or tid) and keep the best
        groups = collections.defaultdict(list)
        for r in R:
            if r['side'] == 'ours_only': continue
            k = (r['haplotype'], r['annot_gene'], r['transcript'] if r['transcript'] != '.' else r['annot_tx'])
            groups[k].append(r)
        def score(r):
            if r['side'] != 'both': return 10 ** 9
            return 0 if r['ed_key'] == '.' else ed.get(r['ed_key'], 10 ** 6)
        cls = collections.Counter(); ex = collections.defaultdict(list)
        n = ident = 0; alt_dis = 0
        for k, g in groups.items():
            g.sort(key=score); b = g[0]; n += 1
            if b['side'] == 'both' and b['ed_key'] == '.':
                ident += 1
                if any(score(x) > 0 for x in g[1:]):
                    alt_dis += 1; cls['cat_alt_alignment_differs (best identical)'] += 1
                    if len(ex['cat_alt_alignment_differs (best identical)']) < 3:
                        x = [x for x in g if score(x) > 0][0]
                        ex['cat_alt_alignment_differs (best identical)'].append('%s %s %s %s ed=%s' % (short(x['haplotype']), x['gene'], x['transcript'], x['annot_class'], ed.get(x['ed_key'])))
                continue
            # classify the disagreement
            c = None
            note = b['note']
            if b['side'] == 'annot_only':
                if 'novel' in b['annot_class'] or (annot == 'cat' and b['transcript'] == '.'):
                    c = 'annot novel isoform (no ref transcript)'
                elif annot == 'lo':
                    c = 'annot transcript absent from our GTF (RefSeq isoform)'
                elif 'gtf_gene:' in note:
                    c = 'annotation naming (ENST under another gene)'
                elif 'not_in_gtf' in note:
                    c = 'annot transcript absent from our GTF (release)'
                elif 'ours_rows:' in note and ':ok' in note.split('ours_rows:')[1] or ':boundary_moved' in note or ':splice_site' in note:
                    ov = note.split('overlaps:')[1].split(';')[0]
                    c = 'annot extra copy, on our %s' % ov
                elif 'ours_rows:none' in note and annot == 'cat':
                    c = 'not projected by us (transcript reaches past the GRCh38 path)'
                elif 'ours_detail:' in note:
                    d = note.split('ours_detail:')[1]
                    if 'taken_by:' in d: c = 'copy assigned to a different paralog (ours: taken_by %s)' % d.split('taken_by:')[1].split(',')[0]
                    elif 'inside:' in d: c = 'copy nested in another gene (ours: inside %s)' % d.split('inside:')[1].split(',')[0]
                    else: c = 'partial copy below --min-exon-share (ours: no_copy, %s)' % d.split(',')[0]
                elif 'ours_rows:none' in note or 'no_copy' in note:
                    ov = note.split('overlaps:')[1].split(';')[0]
                    c = 'annot copy on our other paralog (%s)' % ov if ov != 'none' and b['gene'] not in ov.split('/') else 'missing on our side (no copy)'
                elif 'missing_exon' in note:
                    c = 'missing on our side (missing_exon)'
                else:
                    c = 'missing on our side (other: %s)' % note[:60]
            else:
                enst = b['transcript']
                cver = None
                if 'novel' in b['annot_class']:
                    c = 'annot novel isoform (matched by exons/introns)'
                elif annot == 'lo' and b['match'] == 'mane':
                    c = 'MANE: RefSeq vs Ensembl MANE differ (same introns)' if b['same_introns'] == 'yes' else 'MANE: different intron chain (inspect)'
                elif annot == 'lo':
                    c = 'RefSeq isoform differs at ends (same introns)' if b['same_introns'] == 'yes' else 'RefSeq isoform, different intron chain'
                elif 'ver:' in note and note.split('ver:')[1].split(';')[0].split('/')[0] != note.split('ver:')[1].split(';')[0].split('/')[1]:
                    c = 'annotation release (transcript version differs)'
                elif b['our_status'] == 'boundary_moved':
                    c = 'boundary we moved'
                else:
                    c = 'exons differ, same ref version (inspect)'
                b['_ed'] = ed.get(b['ed_key'])
            cls[c] += 1
            if len(ex[c]) < 3:
                ex[c].append('%s %s %s %s %s ed=%s %s' % (short(b['haplotype']), b['gene'], b['transcript'], b['annot_class'], b['our_status'], b.get('_ed', '.'), note[:140]))
        print('transcripts compared (annot records, best alignment per source transcript): %d; identical: %d (%.1f%%)' % (n, ident, 100.0 * ident / max(n, 1)))
        if annot == 'lo':
            mg = [sorted(g, key=score)[0] for g in groups.values() if g[0]['annot_class'].startswith('MANE')]
            mi = sum(1 for b in mg if b['side'] == 'both' and b['ed_key'] == '.')
            print('  MANE Select records: %d; matched to our Ensembl MANE: %d; identical: %d' % (len(mg), sum(1 for b in mg if b['side'] == 'both'), mi))
            eds = sorted(ed.get(b['ed_key'], -1) for b in mg if b['side'] == 'both' and b['ed_key'] != '.')
            print('  MANE differing edit distances:', eds)
        for c, k in cls.most_common():
            print('  %5d  %s' % (k, c))
            for e in ex[c]: print('           e.g. ' + e)
        if annot == 'cat':
            oo = [r for r in R if r['side'] == 'ours_only']
            print('our HSTs without a CAT partner: %d' % len(oo), dict(collections.Counter(r['note'] for r in oo)))


def fetch_liftoff(p, u, out):
    """A Liftoff GFF3's lines over path p's interval, by binary search over the sorted remote file."""
    u = u.replace('s3://human-pangenomics', 'https://s3-us-west-2.amazonaws.com/human-pangenomics')
    c = p.rsplit(':', 1)[0]  # sample#hap#contig
    s, e = map(int, p.rsplit(':', 1)[1].split('-'))
    def get(a, b):
        r = urllib.request.Request(u, headers={'Range': 'bytes=%d-%d' % (a, b)})
        return urllib.request.urlopen(r, timeout=60).read()
    size = int(urllib.request.urlopen(urllib.request.Request(u, method='HEAD'), timeout=60).headers['Content-Length'])
    def key_at(off):
        b = get(off, off + 20000).split(b'\n')
        for l in b[1:-1]:
            if l and not l.startswith(b'#'):
                t = l.split(b'\t')
                return (t[0].decode(), int(t[3]))
        return None
    target = (c, s - 3000000)
    lo, hi = 0, size
    while hi - lo > 200000:
        mid = (lo + hi) // 2
        k = key_at(mid)
        if k is None or k < target: lo = mid
        else: hi = mid
    # stream from lo
    r = urllib.request.Request(u, headers={'Range': 'bytes=%d-' % lo})
    f = urllib.request.urlopen(r, timeout=120)
    buf = b''; n = 0; seen = False
    with open(out + '.tmp', 'w') as o:
        first = True
        while True:
            ch = f.read(1 << 20)
            if not ch: break
            buf += ch
            lines = buf.split(b'\n'); buf = lines.pop()
            if first: lines = lines[1:]; first = False
            stop = False
            for l in lines:
                if not l or l.startswith(b'#'): continue
                t = l.decode().split('\t')
                if t[0] == c:
                    seen = True
                    a, b = int(t[3]), int(t[4])
                    if b > s and a <= e: o.write(l.decode() + '\n'); n += 1
                    if a > e + 3000000: stop = True; break
                elif seen or t[0] > c: stop = True; break
            if stop: break
    os.rename(out + '.tmp', out)
    print(p, n, 'seen' if seen else 'NOT_SEEN')


def main():
    cmd, D = sys.argv[1], sys.argv[2]
    if cmd == 'select':
        locus, hst_tsv, cat_csv, lo_csv, gene, plan = sys.argv[3:9]
        select(D, locus, hst_tsv, cat_csv, lo_csv, gene,
               {int(k): int(v) for k, v in (x.split('=') for x in plan.split(','))}, sys.argv[9:])
    elif cmd == 'fetch':
        for annot in ('cat', 'lo'):
            os.makedirs(os.path.join(D, annot), exist_ok=True)
        for r in (l.rstrip('\n').split('\t') for l in open(D + '/selected.tsv')):
            hp, name = r[1], r[1].replace('#', '_').replace(':', '_')
            for annot, url in (('cat', r[4]), ('lo', r[5])):
                out = '%s/%s/%s.%s.gff3' % (D, annot, name, annot)
                if url != '.' and not os.path.exists(out):
                    (fetch_cat if annot == 'cat' else fetch_liftoff)(hp, url, out)
    elif cmd == 'compare':
        locus, hst_tsv, gfa, gtf = sys.argv[3:7]
        ed = edit_distances(compare(D, locus, hst_tsv, gfa, gtf))
        summary(D, locus, ed)


if __name__ == '__main__':
    main()
