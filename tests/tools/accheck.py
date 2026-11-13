#!/usr/bin/env python3
"""accheck.py [--min-af X] <vcf>...: every record's AC is the number of haplotypes with a carrying column, and none is 0.

Per path, a haplotype is a PanSN sample#haplotype (numbered), or a path of its own; the reference (##reference)
is none. Per sample (--diploid), each haplotype of a column's GT that names an allele. With --min-af, every record's
AF must also reach X. Prints each violation; exit 1 if any."""
import sys

def check(path, min_af):
    bad, ref, diploid, samples = [], None, False, []
    for line in open(path):
        if line.startswith('##reference='):
            ref = line.rstrip('\n').split('=', 1)[1]
        elif line.startswith('##panvar_columns=diploid'):
            diploid = True
        elif line.startswith('#CHROM'):
            samples = line.rstrip('\n').split('\t')[9:]
        elif not line.startswith('#'):
            f = line.rstrip('\n').split('\t')
            info = dict(kv.split('=', 1) for kv in f[7].split(';') if '=' in kv)
            gi = f[8].split(':').index('GT')
            carriers = set()
            for name, cell in zip(samples, f[9:]):
                if name == ref:
                    continue
                gt = cell.split(':')[gi].replace('|', '/').split('/')
                if not any(a not in ('0', '.', '') for a in gt):
                    continue
                if diploid:
                    carriers.update(f'{name}\x00{k}' for k, a in enumerate(gt) if a not in ('0', '.', ''))
                else:
                    p = name.split('#')
                    carriers.add('#'.join(p[:2]) if len(p) > 1 and p[1].isdigit() else '\x01' + name)
            ac = int(info.get('AC', '-1'))
            if ac != len(carriers):
                bad.append(f'{path}: {f[2]} AC={ac}, carrying haplotypes {len(carriers)}')
            if ac == 0:
                bad.append(f'{path}: {f[2]} written with AC=0')
            if min_af is not None and float(info.get('AF', '0')) < min_af:
                bad.append(f'{path}: {f[2]} AF={info.get("AF")} under --min-alt-af {min_af}')
    return bad

args = sys.argv[1:]
min_af = None
if args[:1] == ['--min-af']:
    min_af = float(args[1]); args = args[2:]
problems = [b for p in args for b in check(p, min_af)]
for b in problems[:20]:
    print(b)
sys.exit(1 if problems else 0)
