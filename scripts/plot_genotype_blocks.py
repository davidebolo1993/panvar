#!/usr/bin/env python3
"""Per-block genotyping error, one figure per locus.

Two rows: leave-zero-out above (the control -- the donor's own haplotypes are IN the panel, so a
perfect model reads 0), leave-one-out below (the donor is held out, so part of the error is a
haplotype the panel no longer contains). One violin pair per block, real reads against simulated.

Metric is normalized edit distance to truth, 1 - identity, so blocks of very different length share
an axis. The same metric is used in both rows on purpose: mixing "distance to truth" above with
"distance to the best available pair" below would make the two rows incomparable.

The panel floor (1 - best_identity) is drawn as a black diamond. Under leave-one-out the violin
ABOVE that diamond is the model's own error; the diamond is what removing the donor costs. Plotting
the violin alone would charge the model for a haplotype nobody could have picked.

NOTE `best_identity` is computed over a top-16 syncmer-Jaccard shortlist, so it is a LOWER bound on
the true panel ceiling: the real floor may be lower and the model's own share correspondingly larger.

  plot_genotype_blocks.py --table all_blocks.tsv --out DIR [--loci c4,lpa] [--dpi 150]
"""
import argparse, os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

PAL = {'real': '#B2182B', 'simulated': '#2166AC'}
ROWS = [('LZO', 'leave-zero-out  (donor kept in panel)'),
        ('LOO', 'leave-one-out  (donor removed)')]

ap = argparse.ArgumentParser(description=__doc__,
                             formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--table', required=True, help='per-block cohort table')
ap.add_argument('--out', required=True, help='output directory')
ap.add_argument('--loci', default=None, help='comma-separated subset')
ap.add_argument('--dpi', type=int, default=150)
a = ap.parse_args()

d = pd.read_csv(a.table, sep='\t')
need = {'locus', 'read_source', 'regime', 'block_index', 'block_kind',
        'sample', 'identity', 'best_identity'}
miss = need - set(d.columns)
if miss:
    sys.exit(f'table lacks columns: {", ".join(sorted(miss))}')

d['identity'] = pd.to_numeric(d['identity'], errors='coerce')
d['best_identity'] = pd.to_numeric(d['best_identity'], errors='coerce')
d = d[np.isfinite(d['identity'])].copy()
d['err'] = 1.0 - d['identity']
d['floor'] = 1.0 - d['best_identity']
if 'exact' in d.columns:
    d['exact'] = pd.to_numeric(d['exact'], errors='coerce')
else:
    d['exact'] = np.nan
if a.loci:
    d = d[d['locus'].isin(a.loci.split(','))]
if d.empty:
    sys.exit('no rows to plot')

os.makedirs(a.out, exist_ok=True)

# Per-block summary written alongside, since "what percent did we genotype properly" is a number
# people want to quote and reading it off a violin is guesswork.
rows = []
for (loc, reg, blk, kind, src), g in d.groupby(
        ['locus', 'regime', 'block_index', 'block_kind', 'read_source'], sort=True):
    rows.append(dict(locus=loc, regime=reg, block_index=blk, block_kind=kind, read_source=src,
                     n=len(g), exact_pct=100.0 * g['exact'].mean(skipna=True),
                     mean_err=g['err'].mean(), median_err=g['err'].median(),
                     max_err=g['err'].max(), panel_floor=g['floor'].mean(skipna=True)))
summary = pd.DataFrame(rows).sort_values(['locus', 'regime', 'block_index', 'read_source'])
summary.to_csv(os.path.join(a.out, 'summary.tsv'), sep='\t', index=False,
               float_format='%.6f')

for loc in sorted(d['locus'].unique()):
    sub = d[d['locus'] == loc]
    blocks = sorted(sub['block_index'].unique())
    kind = {b: str(sub.loc[sub['block_index'] == b, 'block_kind'].iloc[0]) for b in blocks}
    fig, axes = plt.subplots(2, 1, figsize=(max(7.0, 0.46 * len(blocks) + 2.6), 7.2), sharex=True)

    for ax, (reg, title) in zip(axes, ROWS):
        r = sub[sub['regime'] == reg]
        for off, src in ((-0.19, 'real'), (0.19, 'simulated')):
            xs, data = [], []
            for i, b in enumerate(blocks):
                v = r.loc[(r['block_index'] == b) & (r['read_source'] == src), 'err'].dropna()
                if len(v) == 0:
                    continue
                # A block where every sample is exact has zero variance; violin() fails on that, so
                # those are drawn as a flat marker rather than silently dropped.
                if len(v) > 1 and v.nunique() > 1:
                    xs.append(i + off); data.append(v.values)
                else:
                    ax.plot([i + off - 0.13, i + off + 0.13], [v.iloc[0]] * 2,
                            color=PAL[src], lw=1.6, solid_capstyle='butt', zorder=3)
                jit = (np.random.default_rng(42).random(len(v)) - 0.5) * 0.16
                ax.scatter(i + off + jit, v.values, s=3, color=PAL[src], alpha=0.45,
                           linewidths=0, zorder=4)
            if data:
                parts = ax.violinplot(data, positions=xs, widths=0.34, showextrema=False)
                for pc in parts['bodies']:
                    pc.set_facecolor(PAL[src]); pc.set_alpha(0.45); pc.set_edgecolor('none')
        fl = [r.loc[r['block_index'] == b, 'floor'].mean() for b in blocks]
        ax.scatter(range(len(blocks)), fl, marker='D', s=14, color='black', zorder=5,
                   label='panel floor (1 - best_identity)')
        ax.set_title(title, fontsize=10, loc='left', fontweight='bold')
        ax.set_ylabel('1 - identity')
        ax.grid(axis='y', lw=0.4, alpha=0.35)
        ax.set_axisbelow(True)
        ax.margins(x=0.01)

    axes[1].set_xticks(range(len(blocks)))
    axes[1].set_xticklabels([f'{b}\n{kind[b][:4]}' for b in blocks], fontsize=7)
    axes[1].set_xlabel('block (chain order; kind abbreviated)')
    handles = [Patch(facecolor=PAL['real'], alpha=0.6, label='real reads'),
               Patch(facecolor=PAL['simulated'], alpha=0.6, label='simulated reads'),
               plt.Line2D([], [], marker='D', color='black', ls='none', ms=4,
                          label='panel floor (unreachable under LOO)')]
    axes[0].legend(handles=handles, loc='upper right', fontsize=7, framealpha=0.9)
    ns = sub.groupby('read_source')['sample'].nunique().to_dict()
    fig.suptitle(f'{loc} - per-block genotyping error   '
                 f'({ns.get("real", 0)} real / {ns.get("simulated", 0)} simulated donors)',
                 fontsize=12, fontweight='bold')
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    p = os.path.join(a.out, f'{loc}.png')
    fig.savefig(p, dpi=a.dpi); plt.close(fig)
    print(f'wrote {p}  ({len(blocks)} blocks, {len(sub)} observations)')

print(f'wrote {os.path.join(a.out, "summary.tsv")}')
