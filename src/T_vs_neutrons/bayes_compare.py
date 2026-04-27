#!/usr/bin/env python
"""Bayesian comparison of TRANSP runs against RNT and TH/BT measurements.

For each run, time-averages NEUTT and the TH/BT ratio (NEUTX/BTNTS) over a
window, compares against the measured values with Gaussian likelihoods, and
returns posterior probabilities under a uniform prior. Assumes the two
measurements are independent.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import argparse
import profiles as ps

parser = argparse.ArgumentParser(
    description='Bayesian comparison of TRANSP runs based on total neutron rate '
                'and TH/BT ratio measurements.')
parser.add_argument('pulse', type=int)
parser.add_argument('runs', nargs='+',
                    help='TRANSP run IDs to compare (e.g. M26 M27 M28)')
parser.add_argument('-tw', '--time_window', nargs=2, type=float, required=True,
                    metavar=('T1', 'T2'),
                    help='Wall-time window [s] over which to average (e.g. -tw 51.5 52.5)')
parser.add_argument('--thbt', type=float, required=True,
                    help='Measured TH/BT ratio (central value)')
parser.add_argument('--thbt_sigma', type=float, default=0.06,
                    help='Absolute 1-sigma uncertainty on TH/BT (default: 0.06)')
parser.add_argument('--rnt_unc', type=float, default=0.10,
                    help='Relative 1-sigma uncertainty on RNT (default: 0.10)')
args = parser.parse_args()

t1, t2 = args.time_window

def window_mean(t, values, dt=None):
    """Mean of values over [t1, t2]; weighted by dt if provided. NaN if window empty."""
    mask = (t >= t1) & (t <= t2)
    if not np.any(mask):
        return np.nan
    if dt is not None:
        return np.average(values[mask], weights=dt[mask])
    return np.mean(values[mask])

print(f'\nPulse {args.pulse}, time window [{t1:.3f}, {t2:.3f}] s\n')

# --- load runs and collect time-averaged predictions ---
runs_data = []
rnt_meas = None
for rid in args.runs:
    n = ps.Neutrons(args.pulse, rid)
    n.get_transp_neutrons()
    if rnt_meas is None:
        n.get_exp()
        rnt_t, rnt = n.rnt
        rnt_meas = window_mean(rnt_t, rnt)

    twall = n.transp_time + 40.
    neutt = n.transp('NEUTT')
    if neutt is None:
        print(f'WARNING: {rid} has no NEUTT - skipping')
        continue
    neutt_pred = window_mean(twall, neutt, dt=n.dt)

    neutx = n.transp('NEUTX')
    btnts = n.transp('BTNTS')
    if neutx is not None and btnts is not None:
        neutx_avg = window_mean(twall, neutx, dt=n.dt)
        btnts_avg = window_mean(twall, btnts, dt=n.dt)
        thbt_pred = neutx_avg / btnts_avg if btnts_avg > 0 else np.nan
    else:
        thbt_pred = np.nan
    runs_data.append((rid, neutt_pred, thbt_pred))

if not runs_data or rnt_meas is None or np.isnan(rnt_meas):
    print('\nERROR: no usable runs or no RNT data in the window. Aborting.')
    sys.exit(1)

# --- likelihoods ---
sigma_rnt  = args.rnt_unc * rnt_meas
sigma_thbt = args.thbt_sigma

print('Measurement:')
print(f'  <RNT>  = {rnt_meas:.3e} n/s   (sigma = {sigma_rnt:.3e}, {args.rnt_unc*100:.1f}%)')
print(f'  TH/BT  = {args.thbt:.3f}        (sigma = {sigma_thbt:.3f})\n')

results = []
for rid, neutt_pred, thbt_pred in runs_data:
    c2_rnt  = ((neutt_pred - rnt_meas) / sigma_rnt)**2 if not np.isnan(neutt_pred) else 0.0
    c2_thbt = ((thbt_pred - args.thbt) / sigma_thbt)**2 if not np.isnan(thbt_pred) else 0.0
    log_L   = -0.5 * (c2_rnt + c2_thbt)
    results.append((rid, neutt_pred, thbt_pred, c2_rnt, c2_thbt, log_L))

# --- posterior under uniform prior ---
log_Ls = np.array([r[5] for r in results])
posteriors = np.exp(log_Ls - log_Ls.max())
posteriors /= posteriors.sum()

# --- print table ---
header = (f'{"Run":<8}{"NEUTT [n/s]":<16}{"TH/BT":<10}'
          f'{"chi2(RNT)":<12}{"chi2(TH/BT)":<14}{"log L":<12}{"P(M|d)":>10}')
print(header)
print('-' * len(header))
for (rid, neutt, thbt, c2r, c2t, lnL), post in zip(results, posteriors):
    thbt_str = f'{thbt:.3f}' if not np.isnan(thbt) else 'n/a'
    print(f'{rid:<8}{neutt:<16.3e}{thbt_str:<10}'
          f'{c2r:<12.3f}{c2t:<14.3f}{lnL:<12.3f}{post*100:>9.2f}%')

best = int(np.argmax(posteriors))
print(f'\nMost likely: {results[best][0]}  (P = {posteriors[best]*100:.1f}%)')
print('\nAssumes RNT and TH/BT measurements are independent.')
print('Uniform prior across the supplied runs.')
