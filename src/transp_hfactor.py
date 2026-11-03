#!/usr/bin/env python
"""
Confinement time and H-factor summary from a TRANSP run.

Prints the mean, standard deviation and median of the thermal (TAUEA), total
(TAUA1) and electron (TAUEE) energy confinement times and of the chosen
H-factors, over a time window or over the whole run.

Definitions (verified on 103073/M05 and 104614/M30):
  H98Y*  = TAUEA / tau_IPB98(y | y,1 | y,2)
  TAUEA  = W_th / (P_loss - dW_th/dt)       thermal, dW/dt-corrected
  TAUA1  = (W_th + W_fast) / P_loss         total, NO dW/dt correction

Scaling laws IPB98(y), IPB98(y,1) and IPB98(y,2) (TRANSP: H98Y, H98Y1, H98Y2):
  ITER Physics Expert Groups on Confinement and Transport and Confinement
  Modelling and Database, ITER Physics Basis Editors,
  "Chapter 2: Plasma confinement and transport",
  Nucl. Fusion 39, 2175 (1999). IPB98(y,2) is the standard H98 reference.

Several run IDs compare runs of the same pulse: a side-by-side table of
mean +/- std, the change of each mean relative to the first run, and (with
--plot) one colour per run. Without --window the statistics are taken over the
time range common to all runs, so the averages are directly comparable.

Examples:
  python transp_hfactor.py 104614 M30
  python transp_hfactor.py 104614 M30 M31 M32 --window 52.5 54.0 --plot
  python transp_hfactor.py 104614 M30 --window 52.5 54.0
  python transp_hfactor.py 104614 M30 --window 52.5 54.0 --plot
  python transp_hfactor.py 104614 M30 --hfactors H98Y H98Y1 H98Y2 --plot --save h.png
"""
import argparse
import os
import sys

import numpy as np
from netCDF4 import Dataset

TRANSP_PATH = '/common/transp_shared/Data/result/JET'
HFACTORS = ('H98Y', 'H98Y1', 'H98Y2')
# Official scaling-law names (ITER Physics Basis, Nucl. Fusion 39, 2175, 1999)
SCALING_NAMES = {'H98Y': 'IPB98(y)', 'H98Y1': 'IPB98(y,1)', 'H98Y2': 'IPB98(y,2)'}
TAUS = ('TAUEA', 'TAUA1', 'TAUEE')


def cdf_path(pulse, runid, root=TRANSP_PATH):
    pulse = str(pulse)
    return os.path.join(root, pulse, runid, pulse + runid + '.CDF')


def read_signal(f, key):
    """Return (time, values) for a 1-D time trace, NaN where masked."""
    if key not in f.variables:
        raise KeyError('Signal %s not found in CDF' % key)
    var = f.variables[key]
    tname = var.dimensions[0]  # TIME or TIME3 depending on the signal
    t = np.ma.filled(np.ma.asarray(f.variables[tname][:], dtype=float), np.nan)
    y = np.ma.filled(np.ma.asarray(var[:], dtype=float), np.nan)
    return t, y


def stats(y):
    y = y[np.isfinite(y)]
    if y.size == 0:
        return np.nan, np.nan, np.nan, 0
    return np.mean(y), np.std(y), np.median(y), y.size


def load_run(pulse, runid, keys, root, toffset):
    """Read the requested signals of one run; None if the CDF is missing."""
    path = cdf_path(pulse, runid, root)
    if not os.path.isfile(path):
        print('WARNING: CDF not found, skipping %s: %s' % (runid, path))
        return None
    data = {}
    with Dataset(path) as f:
        for k in keys:
            t, y = read_signal(f, k)
            data[k] = (t + toffset, y)
    return {'path': path, 'data': data}


def main():
    p = argparse.ArgumentParser(
        description='TRANSP confinement time and H-factor summary; '
                    'give several run IDs to compare runs of the same pulse.')
    p.add_argument('pulse', type=int, help='JET pulse number, e.g. 104614')
    p.add_argument('runids', nargs='+', metavar='runid',
                   help='one or more TRANSP run IDs, e.g. M30 M31 M32')
    p.add_argument('--window', nargs=2, type=float, metavar=('T1', 'T2'),
                   help='time window in JET time [s] (TRANSP time + toffset); '
                        'default: whole run (common time range if several runs)')
    p.add_argument('--toffset', type=float, default=40.0,
                   help='offset added to TRANSP time to get JET time [s] (default 40)')
    p.add_argument('--hfactors', nargs='+', default=['H98Y2'], choices=HFACTORS,
                   help='H-factors to report/plot (default: H98Y2)')
    p.add_argument('--plot', action='store_true', help='plot H-factors and TAUEA')
    p.add_argument('--save', metavar='FILE',
                   help='save the plot to FILE instead of showing it')
    p.add_argument('--root', default=TRANSP_PATH, help='TRANSP results root directory')
    args = p.parse_args()

    keys = list(TAUS) + list(args.hfactors)
    runs = {}
    for r in dict.fromkeys(args.runids):          # drop duplicates, keep order
        run = load_run(args.pulse, r, keys, args.root, args.toffset)
        if run is not None:
            runs[r] = run
    if not runs:
        sys.exit('No CDFs found for pulse %d.' % args.pulse)

    # time range of each run and the window to average over
    ranges = {r: (np.nanmin(run['data']['TAUEA'][0]), np.nanmax(run['data']['TAUEA'][0]))
              for r, run in runs.items()}
    if args.window:
        t1, t2 = sorted(args.window)
        win_txt = '%.3f - %.3f s' % (t1, t2)
    else:
        t1 = max(a for a, b in ranges.values())
        t2 = min(b for a, b in ranges.values())
        if t1 >= t2:
            sys.exit('Runs do not overlap in time; give --window explicitly.')
        win_txt = ('whole run' if len(runs) == 1 else 'common time range') + \
                  ' (%.3f - %.3f s)' % (t1, t2)

    print('Pulse %d' % args.pulse)
    for r, run in runs.items():
        print('  %-6s %s  (%.3f - %.3f s)' % ((r, run['path']) + ranges[r]))
    print('Window %s  (JET time, toffset = %g s)' % (win_txt, args.toffset))

    for r, (a, b) in ranges.items():
        if b < t1 or a > t2:
            sys.exit('Run %s (%.3f-%.3f s) has no samples in window %.3f-%.3f s. '
                     'Check --window and --toffset.' % (r, a, b, t1, t2))
        if args.window and (a > t1 or b < t2):
            print('WARNING: run %s covers only part of the window' % r)

    # statistics: res[run][signal] = (mean, std, median, N)
    res = {}
    for r, run in runs.items():
        res[r] = {}
        for k in keys:
            t, y = run['data'][k]
            res[r][k] = stats(y[(t >= t1) & (t <= t2)])

    if len(runs) == 1:
        r = next(iter(runs))
        print('-' * 64)
        print('%-8s %12s %12s %12s %6s' % ('signal', 'mean', 'std', 'median', 'N'))
        for k in keys:
            mean, std, med, n = res[r][k]
            unit = ' s' if k.startswith('TAU') else ''
            print('%-8s %12.4f %12.4f %12.4f %6d%s' % (k, mean, std, med, n, unit))
    else:
        ids = list(runs)
        w = 18
        line = '-' * (10 + (w + 1) * len(ids))
        print(line)
        print('Mean +/- std over the window')
        print('%-10s' % 'signal' + ''.join(' %*s' % (w, r) for r in ids))
        for k in keys:
            cells = ['%.4f +/- %.4f' % res[r][k][:2] for r in ids]
            print('%-10s' % (k + (' [s]' if k.startswith('TAU') else ''))
                  + ''.join(' %*s' % (w, c) for c in cells))
        print(line)
        print('Mean relative to %s [%%]' % ids[0])
        print('%-10s' % 'signal' + ''.join(' %*s' % (w, r) for r in ids[1:]))
        for k in keys:
            ref = res[ids[0]][k][0]
            cells = ['%+.2f' % (100.0 * (res[r][k][0] / ref - 1.0)) for r in ids[1:]]
            print('%-10s' % k + ''.join(' %*s' % (w, c) for c in cells))

    if args.plot:
        import matplotlib
        if args.save or not os.environ.get('DISPLAY'):
            matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        # colour = run, line style = scaling law
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        styles = {'H98Y2': '-', 'H98Y1': '--', 'H98Y': ':'}
        multi = len(runs) > 1

        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(9, 7))
        for i, (r, run) in enumerate(runs.items()):
            c = colors[i % len(colors)] if multi else None
            for k in args.hfactors:
                t, y = run['data'][k]
                if multi:
                    label = r if len(args.hfactors) == 1 else '%s  %s' % (r, SCALING_NAMES[k])
                    ax1.plot(t, y, color=c, ls=styles[k], label=label)
                else:
                    ax1.plot(t, y, label=SCALING_NAMES[k])
            t, y = run['data']['TAUEA']
            ax2.plot(t, y, color=c if multi else 'k',
                     label=r if multi else 'TAUEA (thermal)')
        ylab = r'$H = \tau_{E,th}/\tau_E^{\mathrm{scaling}}$'
        if multi and len(args.hfactors) == 1:
            ylab = SCALING_NAMES[args.hfactors[0]] + '  ' + ylab
        ax1.set_ylabel(ylab)
        ax1.legend(fontsize='small')
        ax2.set_ylabel(r'$\tau_{E,th}$ (TAUEA) [s]')
        ax2.set_xlabel('JET time [s]')
        ax2.legend(fontsize='small')
        for ax in (ax1, ax2):
            ax.grid(alpha=0.3)
            if args.window or multi:
                ax.axvspan(t1, t2, color='0.9', zorder=0)
        fig.suptitle('%d / %s' % (args.pulse, ', '.join(runs)))
        fig.tight_layout()

        if args.save or not os.environ.get('DISPLAY'):
            out = args.save or '%d%s_hfactor.png' % (args.pulse, '_'.join(runs))
            fig.savefig(out, dpi=150)
            print('Plot saved to %s' % out)
        else:
            plt.show()


if __name__ == '__main__':
    main()
