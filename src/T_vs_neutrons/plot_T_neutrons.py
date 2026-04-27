#!/usr/bin/env python
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import plotWindow as pw
import matplotlib.pyplot as plt
import numpy as np
import profiles as ps
import argparse

parser = argparse.ArgumentParser(
    description='Check how the tritium profile affects neutron rate and TH/BT ratio.')
parser.add_argument('pulse', type=int)
parser.add_argument('runid', type=str)
parser.add_argument('-t', '--time', type=float, default=None,
                    help='Time slice for profile plot [s]. Default: last time slice.')
parser.add_argument('-nf', '--no_fast', action='store_true',
                    help='Also plot NT/(NT+ND) profile without fast ion density.')
parser.add_argument('-c', '--compare', nargs='+', default=None, metavar='RUNID',
                    help='Run IDs to compare, e.g. -c M27 or -c M27 M28.')
parser.add_argument('-hs', '--hstrap', nargs=2, type=float, metavar=('YMIN', 'YMAX'),
                    help='Shade a horizontal band on the TH/BT ratio plot (diagnostic measurement).')
parser.add_argument('--strap_color', default='lightgrey',
                    help='Colour of the strap (default: lightgrey).')
parser.add_argument('--strap_alpha', type=float, default=0.5,
                    help='Transparency of the strap (default: 0.5).')
args = parser.parse_args()

pulse = args.pulse
runid = args.runid

COLORS = ['k', 'royalblue', 'tomato', 'forestgreen', 'darkorange', 'purple']
BND    = 0.2   # inner-region rho boundary for TH/BT(rho < BND)

def safe_divide(x, y):
    result = np.zeros_like(x, dtype=float)
    mask = y > 0
    result[mask] = x[mask] / y[mask]
    return result

def thbt_inner(run, bnd):
    """TH/BT ratio cumulative-integrated up to rho = bnd, at each time step."""
    th_c = np.cumsum(run.get_key('THNTX') * run.dvol, axis=1)
    bt_c = np.cumsum(run.get_key('BTNTX') * run.dvol, axis=1)
    ratio = safe_divide(th_c, bt_c)
    result = np.zeros(run.t.shape)
    for i in range(len(run.t)):
        i_bnd = np.abs(run.x[i, :] - bnd).argmin()
        result[i] = ratio[i, i_bnd]
    return result

def get_nt_profile(run, tidx):
    nt = run._transp['NT'][tidx]
    nd = run._transp['ND'][tidx] if 'ND' in run._transp else np.zeros_like(nt)
    bd = run._transp['BDENS'][tidx] if 'BDENS' in run._transp else np.zeros_like(nt)
    return nt, nd, bd

# --- load runs ---
main_run = ps.Neutrons(pulse, runid)
main_run.get_transp_neutrons()
main_run.get_exp()
main_run.add_data('BDENS')

compare_runs = []
if args.compare:
    for cid in args.compare:
        r = ps.Neutrons(pulse, cid)
        r.get_transp_neutrons()
        r.add_data('BDENS')
        compare_runs.append(r)
all_runs = [main_run] + compare_runs

rnt_time, rnt = main_run.rnt
rnt_unc = 0.1

win = pw.plotWindow()

# =========================================================
# Tab 1: th vs. beam
# =========================================================
fig = plt.figure()
fig.suptitle('Thermal vs. beam', fontsize=13)
ax = fig.add_subplot(111)

ax.plot(rnt_time, rnt, color='r', linewidth=2, label='Measured')
ax.fill_between(rnt_time,
                (1 - rnt_unc) * rnt,
                (1 + rnt_unc) * rnt,
                alpha=0.4, facecolor='r', edgecolor='none')

for run, color in zip(all_runs, COLORS):
    rid = run.runid
    ax.plot(run.transp_time + 40, run.transp('NEUTT'),
            color=color, linewidth=2, label=f'{rid} total')
    if run.transp('NEUTX') is not None:
        ax.plot(run.transp_time + 40, run.transp('NEUTX'),
                color=color, linewidth=2, linestyle='dashed', label=f'{rid} thermal')
    if run.transp('BTNTS') is not None:
        ax.plot(run.transp_time + 40, run.transp('BTNTS'),
                color=color, linewidth=2, linestyle='dotted', label=f'{rid} beam-target')

ax.set_xlabel('Time [s]')
ax.set_ylabel('Neutron rate [#/s]')
ax.set_xlim(main_run.transp_time[0] + 40, main_run.transp_time[-1] + 40)
leg = ax.legend(fontsize='x-small', title='Neutron rate')
leg.set(draggable=True)
win.addPlot('th vs. beam', fig)

# =========================================================
# Tab 2: TH/BT ratio
# =========================================================
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_title(f'{main_run.transpcid} TH/BT ratio')

if args.hstrap:
    ymin, ymax = args.hstrap
    ax.axhspan(ymin, ymax, color=args.strap_color, alpha=args.strap_alpha)

for run, color in zip(all_runs, COLORS):
    rid = run.runid
    if run.transp('BTNTS') is not None:
        ax.plot(run.transp_time + 40,
                safe_divide(run.transp('NEUTX'), run.transp('BTNTS')),
                color=color, linewidth=2, label=f'{rid} average')
    if run.check_signal('BTNTX'):
        ax.plot(run.transp_time + 40, thbt_inner(run, BND),
                color=color, linewidth=2, linestyle='dashed',
                label=rf'{rid} ($\rho<${BND})')

ax.set_xlabel('Time [s]')
ax.set_xlim(main_run.transp_time[0] + 40, main_run.transp_time[-1] + 40)
ax.set_ylim(0, 1.05)
leg = ax.legend()
leg.set_draggable(True)
win.addPlot('TH/BT ratio', fig)

# =========================================================
# Tab 3: NT/(NT+ND+BDENS) profile
# =========================================================
if main_run.signal('NT'):
    t_wall  = main_run.t + 40.
    t_idx   = int(np.argmin(np.abs(t_wall - args.time))) if args.time is not None else -1
    t_slice = t_wall[t_idx]

    show_legend = args.no_fast or len(all_runs) > 1
    if len(all_runs) > 1:
        solid_label = main_run.runid
        dash_label  = main_run.runid + ' (no beam)'
    else:
        solid_label = r'$n_T/(n_T+n_D+n_{beam})$' if args.no_fast else None
        dash_label  = r'$n_T/(n_T+n_D)$'

    nt_p, nd_p, bd_p = get_nt_profile(main_run, t_idx)
    frac_prof = np.where(nt_p + nd_p + bd_p > 0, nt_p / (nt_p + nd_p + bd_p), 0.0)

    fig = plt.figure()
    fig.suptitle('NT/(NT+ND+BDENS) profile', fontsize=13)
    ax = fig.add_subplot(111)
    ax.set_title(f'{main_run.transpcid}  t = {t_slice:.3f} s')
    ax.plot(main_run.x[t_idx], frac_prof * 100, color='k', linewidth=2, label=solid_label)
    if args.no_fast:
        frac_th = np.where(nt_p + nd_p > 0, nt_p / (nt_p + nd_p), 0.0)
        ax.plot(main_run.x[t_idx], frac_th * 100, color='k', linewidth=2,
                linestyle='dashed', label=dash_label)

    for run, color in zip(compare_runs, COLORS[1:]):
        if 'NT' not in run._transp:
            continue
        cp_tidx = int(np.argmin(np.abs(run.t + 40. - t_slice)))
        nt_c, nd_c, bd_c = get_nt_profile(run, cp_tidx)
        frac_c = np.where(nt_c + nd_c + bd_c > 0, nt_c / (nt_c + nd_c + bd_c), 0.0)
        ax.plot(run.x[cp_tidx], frac_c * 100, color=color, linewidth=2, label=run.runid)
        if args.no_fast:
            frac_th_c = np.where(nt_c + nd_c > 0, nt_c / (nt_c + nd_c), 0.0)
            ax.plot(run.x[cp_tidx], frac_th_c * 100, color=color, linewidth=2,
                    linestyle='dashed', label=run.runid + ' (no beam)')

    if show_legend:
        ax.legend()
    ax.set_xlabel(r'$\rho$')
    ax.set_ylabel(r'$n_T\,/\,(n_T+n_D+n_{beam})$ [%]')
    win.addPlot('NT profile', fig)

win.show()
