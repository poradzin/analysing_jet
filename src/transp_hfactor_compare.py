#!/usr/bin/env python
"""
Compare confinement times and H-factors across TRANSP runs of different
pulses, e.g. D, T and DT discharges.

Each entry is PULSE/RUNID, optionally followed by its own time window in JET
time, PULSE/RUNID:T1-T2. Without a window the whole run is averaged. Flat-top
times differ between pulses, so give a window for each entry.

For every entry the script prints mean +/- std of TAUEA, TAUA1, TAUEE and of
the chosen H-factors over its window. --plot draws a summary (mean +/- std per
entry); --traces adds the time traces against JET time, with each averaging
window drawn bold (--relative-time: against time from window start instead).
Entries sharing a --labels value (e.g. several D pulses) share a colour;
--colors sets colours, one per entry or one per distinct label.

Definitions and signal reading come from transp_hfactor.py, which must be in
the same directory:
  H98Y*  = TAUEA / tau_IPB98(y | y,1 | y,2)
  TAUEA  = W_th / (P_loss - dW_th/dt)       thermal, dW/dt-corrected
  TAUA1  = (W_th + W_fast) / P_loss         total, NO dW/dt correction
Scalings IPB98(y), IPB98(y,1), IPB98(y,2): ITER Physics Basis, Chapter 2,
Nucl. Fusion 39, 2175 (1999). Note that the IPB98 scalings contain the isotope
mass (M^0.19 in IPB98(y,2)), so H-factors already include the expected
isotope effect while the raw confinement times do not.

Examples:
  python transp_hfactor_compare.py 99633/M10:48.0-49.5 99950/T05:47.5-49.0 \\
      104614/M30:52.5-54.0 --labels D T DT --plot
  python transp_hfactor_compare.py 99633/M10 99634/M02 99950/T05 \\
      --labels D D T --hfactors H98Y2 H98Y1 --plot --traces --csv out.csv
  python transp_hfactor_compare.py 99633/M10 99634/M02 99950/T05 104614/M30 \\
      --labels D D T DT --colors tab:blue tab:red '#2ca02c' --plot --traces
"""
import argparse
import csv
import os
import re
import sys

import numpy as np
from netCDF4 import Dataset

from transp_hfactor import (TRANSP_PATH, HFACTORS, SCALING_NAMES, TAUS,
                            cdf_path, read_signal, stats)

ENTRY_RE = re.compile(r'^(\d+)/(\w+)(?::([\d.]+)-([\d.]+))?$')
TAU_NAMES = {'TAUEA': r'$\tau_{E,th}$ (TAUEA)',
             'TAUA1': r'$\tau_{E,tot}$ (TAUA1)',
             'TAUEE': r'$\tau_{E,e}$ (TAUEE)'}


def parse_entry(s):
    m = ENTRY_RE.match(s)
    if not m:
        raise argparse.ArgumentTypeError(
            "bad entry '%s'; use PULSE/RUNID or PULSE/RUNID:T1-T2" % s)
    pulse, runid, t1, t2 = m.groups()
    win = tuple(sorted((float(t1), float(t2)))) if t1 else None
    return int(pulse), runid, win


def main():
    p = argparse.ArgumentParser(
        description='Compare TRANSP confinement times and H-factors across pulses.')
    p.add_argument('entries', nargs='+', type=parse_entry, metavar='PULSE/RUNID[:T1-T2]',
                   help='e.g. 104614/M30:52.5-54.0 (window in JET time, optional)')
    p.add_argument('--labels', nargs='+',
                   help='one label per entry, e.g. D T DT; equal labels share a colour')
    p.add_argument('--toffset', type=float, default=40.0,
                   help='offset added to TRANSP time to get JET time [s] (default 40)')
    p.add_argument('--hfactors', nargs='+', default=['H98Y2'], choices=HFACTORS,
                   help='H-factors to report/plot (default: H98Y2)')
    p.add_argument('--plot', action='store_true', help='summary plot of means +/- std')
    p.add_argument('--traces', action='store_true',
                   help='also plot time traces against time from window start')
    p.add_argument('--relative-time', action='store_true',
                   help='plot traces against time from window start instead of JET time')
    p.add_argument('--colors', nargs='+', metavar='COLOR',
                   help='matplotlib colours (names or hex, e.g. tab:blue red #2ca02c): '
                        'one per entry, or one per distinct label')
    p.add_argument('--save', metavar='FILE',
                   help='save the summary plot to FILE (traces go to FILE_traces)')
    p.add_argument('--csv', metavar='FILE', help='write the table to a CSV file')
    p.add_argument('--root', default=TRANSP_PATH, help='TRANSP results root directory')
    args = p.parse_args()

    if args.labels and len(args.labels) != len(args.entries):
        sys.exit('--labels needs one label per entry (%d given, %d entries)'
                 % (len(args.labels), len(args.entries)))
    labels = args.labels or [''] * len(args.entries)
    if args.colors:
        from matplotlib.colors import is_color_like
        bad = [c for c in args.colors if not is_color_like(c)]
        if bad:
            sys.exit('Not a valid colour: %s' % ', '.join(bad))
        n_groups = len(set(labels)) if args.labels else None
        if len(args.colors) != len(args.entries) and len(args.colors) != n_groups:
            sys.exit('--colors needs one colour per entry (%d) or per distinct label (%s); '
                     'got %d' % (len(args.entries), n_groups or 'no --labels', len(args.colors)))
    keys = list(TAUS) + list(args.hfactors)

    # ---- read and average -------------------------------------------------
    rows = []
    for (pulse, runid, win), lab in zip(args.entries, labels):
        path = cdf_path(pulse, runid, args.root)
        if not os.path.isfile(path):
            print('WARNING: CDF not found, skipping %d/%s: %s' % (pulse, runid, path))
            continue
        with Dataset(path) as f:
            data = {}
            for k in keys:
                t, y = read_signal(f, k)
                data[k] = (t + args.toffset, y)
        ta = data['TAUEA'][0]
        a, b = np.nanmin(ta), np.nanmax(ta)
        t1, t2 = win if win else (a, b)
        if b < t1 or a > t2:
            sys.exit('%d/%s covers %.3f-%.3f s, no samples in window %.3f-%.3f s. '
                     'Check the window and --toffset.' % (pulse, runid, a, b, t1, t2))
        if win and (a > t1 or b < t2):
            print('WARNING: %d/%s covers only part of its window' % (pulse, runid))
        res = {}
        for k in keys:
            t, y = data[k]
            res[k] = stats(y[(t >= t1) & (t <= t2)])
        name = '%d/%s' % (pulse, runid)
        rows.append(dict(name=name, label=lab, t1=t1, t2=t2, data=data, res=res,
                         tag=(lab + ' ' + name).strip()))
    if not rows:
        sys.exit('No CDFs found.')

    # ---- print table -------------------------------------------------------
    print('Mean +/- std over each window (JET time, toffset = %g s)' % args.toffset)
    hdr = '%-6s %-14s %-17s' % ('label', 'pulse/run', 'window [s]') + \
          ''.join(' %17s' % (k + (' [s]' if k.startswith('TAU') else '')) for k in keys)
    print('-' * len(hdr))
    print(hdr)
    print('-' * len(hdr))
    for r in rows:
        line = '%-6s %-14s %7.3f - %7.3f' % (r['label'], r['name'], r['t1'], r['t2'])
        for k in keys:
            line += ' %17s' % ('%.4f +/- %.4f' % r['res'][k][:2])
        print(line)

    # group means if labels repeat
    groups = list(dict.fromkeys(r['label'] for r in rows if r['label']))
    if any(sum(r['label'] == g for r in rows) > 1 for g in groups):
        print('-' * len(hdr))
        print('Group mean +/- std of the entry means')
        for g in groups:
            sel = [r for r in rows if r['label'] == g]
            line = '%-6s %-14s %17s' % (g, 'n = %d' % len(sel), '')
            for k in keys:
                v = np.array([r['res'][k][0] for r in sel])
                cell = '%.4f +/- %.4f' % (np.nanmean(v), np.nanstd(v)) if len(sel) > 1 \
                    else '%.4f         ' % v[0]
                line += ' %17s' % cell
            print(line)

    if args.csv:
        with open(args.csv, 'w', newline='') as fh:
            w = csv.writer(fh)
            w.writerow(['label', 'pulse', 'runid', 't1', 't2'] +
                       ['%s_%s' % (k, s) for k in keys for s in ('mean', 'std', 'median', 'N')])
            for r in rows:
                pulse, runid = r['name'].split('/')
                w.writerow([r['label'], pulse, runid, '%.4f' % r['t1'], '%.4f' % r['t2']] +
                           ['%.6g' % x for k in keys for x in r['res'][k]])
        print('Table written to %s' % args.csv)

    if not (args.plot or args.traces):
        return

    # ---- plots -------------------------------------------------------------
    import matplotlib
    headless = bool(args.save) or not os.environ.get('DISPLAY')
    if headless:
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # colours: --colors per entry, or per label group, else default cycle by group
    keyed = [r['label'] or r['name'] for r in rows]
    groups_k = list(dict.fromkeys(keyed))
    if args.colors and len(args.colors) == len(args.entries):
        given = dict(zip(['%d/%s' % e[:2] for e in args.entries], args.colors))
        colors = [given[r['name']] for r in rows]
    elif args.colors:
        given = dict(zip(dict.fromkeys(labels), args.colors))
        colors = [given[g] for g in keyed]
    else:
        palette = plt.rcParams['axes.prop_cycle'].by_key()['color']
        cmap = {g: palette[i % len(palette)] for i, g in enumerate(groups_k)}
        colors = [cmap[g] for g in keyed]
    markers = {'H98Y2': 'o', 'H98Y1': 's', 'H98Y': '^',
               'TAUEA': 'o', 'TAUA1': 's', 'TAUEE': '^'}
    x = np.arange(len(rows))
    figs = []

    if args.plot:
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True,
                                       figsize=(max(6, 1.1 * len(rows) + 3), 7))
        for ax, ks, names in ((ax1, args.hfactors, SCALING_NAMES), (ax2, TAUS, TAU_NAMES)):
            n = len(ks)
            for j, k in enumerate(ks):
                dx = (j - (n - 1) / 2) * 0.18
                for i, r in enumerate(rows):
                    m, s = r['res'][k][:2]
                    ax.errorbar(x[i] + dx, m, yerr=s, fmt=markers[k], color=colors[i],
                                ms=7, capsize=3)
                ax.plot([], [], markers[k], color='0.4', label=names[k])  # legend key
            ax.legend(fontsize='small')
            ax.grid(alpha=0.3)
        ax1.set_ylabel(r'$H = \tau_{E,th}/\tau_E^{\mathrm{scaling}}$')
        ax2.set_ylabel(r'$\tau_E$ [s]')
        ax2.set_xticks(x)
        ax2.set_xticklabels(['%s\n%s' % (r['label'], r['name']) if r['label'] else r['name']
                             for r in rows], fontsize='small')
        for c, tick in zip(colors, ax2.get_xticklabels()):
            tick.set_color(c)
        fig.suptitle('Confinement comparison (mean $\\pm$ std over each window)')
        fig.tight_layout()
        figs.append((fig, args.save or 'hfactor_compare.png'))

    if args.traces:
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(9, 7))
        styles = {'H98Y2': '-', 'H98Y1': '--', 'H98Y': ':'}

        def trace(ax, r, k, c, ls, lab):
            t, y = r['data'][k]
            m = (t >= r['t1']) & (t <= r['t2'])
            if args.relative_time:          # window only, time from window start
                ax.plot(t[m] - r['t1'], y[m], color=c, ls=ls, label=lab)
            else:                           # JET time: full run faint, window bold
                ax.plot(t, y, color=c, ls=ls, lw=0.8, alpha=0.3)
                ax.plot(np.where(m, t, np.nan), np.where(m, y, np.nan),
                        color=c, ls=ls, lw=1.5, label=lab)

        for i, r in enumerate(rows):
            for k in args.hfactors:
                lab = r['tag'] if len(args.hfactors) == 1 else '%s  %s' % (r['tag'], SCALING_NAMES[k])
                trace(ax1, r, k, colors[i], styles[k], lab)
            trace(ax2, r, 'TAUEA', colors[i], '-', r['tag'])
        for ax in (ax1, ax2):
            ax.grid(alpha=0.3)
            ax.legend(fontsize='x-small')
        ax1.set_ylabel(r'$H = \tau_{E,th}/\tau_E^{\mathrm{scaling}}$')
        ax2.set_ylabel(r'$\tau_{E,th}$ (TAUEA) [s]')
        ax2.set_xlabel('time from window start [s]  (each trace shown over its own window)'
                       if args.relative_time else
                       'JET time [s]  (bold: averaging window, faint: rest of run)')
        fig.tight_layout()
        base, ext = os.path.splitext(args.save or 'hfactor_compare.png')
        figs.append((fig, base + '_traces' + (ext or '.png')))

    if headless:
        for fig, out in figs:
            fig.savefig(out, dpi=150)
            print('Plot saved to %s' % out)
    else:
        plt.show()


if __name__ == '__main__':
    main()
