#!/usr/bin/env python
"""make_los_weights.py - time-resolved LoS weight functions for KM9 / KM14.

For each TRANSP TIME3 slice in a run, compute the per-flux-shell LoS detector
coupling weight ``f_det(rhot)`` for one of the JET neutron diagnostics
(KM9 / MPRu or KM14 diamond), and dump the time-resolved arrays to one .npz
file. Geometry-only / channel-independent: the same weight applies to TH and
to (~isotropic) BT emissivities.

The intended consumer is the ``~/jet/jet_tritium`` repository. With this file
it can replace its current crude proxies (KM14: volume integral from 0 to 0.2;
KM9: full-volume integral) with a proper LoS-weighted integral

    Rate_det(t) = sum_i  eps_i(rhot_i, t) * f_det_i(t) * DVOL_i(t)              [n/s]

where ``eps`` is any volumetric emissivity (TH, BT, total, ...) sampled on the
TRANSP rhot grid ``X(t)`` and ``i`` runs over the TRANSP flux shells.
Closure of the construction (verified per slice):
``sum(THNTX * DVOL * f_det) == Rate_det = sum_cells eps * C_cell``.

Output schema (``.npz`` written via ``np.savez_compressed``)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
2-D arrays indexed ``[it, irho]``, ``it`` over selected TIME3 slices, ``irho``
over the TRANSP rhot grid::

    rhot           TRANSP rhot grid X (sorted ascending)
    dvol_m3        Flux-shell volume DVOL [m^3]
    f_det          Detector coupling weight C_bin / DVOL  [dimensionless]
    thntx_si       Reference THNTX [n/m^3/s] for the saved channel
    thlos_det_si   THNTX * f_det  [n/m^3/s]   (closure: sum(.*DVOL) = Rate_det)

1-D arrays indexed ``[it]``::

    times_jet      JET time [s]                (= TRANSP TIME3 + 40 s)
    trinds         original TIME3 index in the main CDF
    rate_det       LoS-weighted thermal rate [n/s]
    rate_chord     bare real-chord emission sum(eps*V) [n/s] (no solid angle)
    rate_closure   sum(THNTX * f_det * DVOL) [n/s] (== rate_det to ~1e-12)
    rho_50         signal-median radius (50 % of rate_det enclosed)
    rhot_min       chord closest-approach rhot (KM9; NaN for vertical KM14)
    rhot_crit      enclosed-shell rhot (KM14; NaN for horizontal KM9)

0-D string / scalar metadata::

    pulse, runid, diagnostic, channel, th_var, eq_source, chord_kind, los_file

Usage
-----
    # KM9, every TIME3 slice, self-contained CDF equilibrium (default)
    python make_los_weights.py 104614 M29 --diag km9

    # KM14, time window 53.0-54.0 s, every 5th slice
    python make_los_weights.py 104614 M29 --diag km14 -t 53.0 54.0 --stride 5

    # Custom output path
    python make_los_weights.py 104614 M29 --diag km9 --output /tmp/foo.npz

Reader sketch (drop into jet_tritium)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    import numpy as np
    from scipy.interpolate import PchipInterpolator

    z = np.load('104614M29_KM9_los_weights.npz')
    t  = z['times_jet']      # (n_t,)
    rhot = z['rhot']         # (n_t, n_rho)
    fdet = z['f_det']        # (n_t, n_rho)
    dvol = z['dvol_m3']      # (n_t, n_rho)

    # If your emissivity eps(rhot, t) lives on a different rhot grid, PCHIP
    # it onto rhot[it] per slice; if you produce it directly on the TRANSP X
    # grid (e.g. straight from the same CDF) you can skip the interp.
    rate_t = np.array([
        np.sum(eps_on_rhot_i * fdet[i] * dvol[i]) for i, eps_on_rhot_i in ...
    ])

Notes
-----
* The weight is *channel-independent*; ``--channel`` only controls which
  reference THNTX is co-stored for self-consistency checking.
* KM9's BT was treated as isotropic in v1, so the same f_det is appropriate
  for BT in jet_tritium (see CONTEXT.md). For KM14 the BT anisotropy factor
  is g ~ 1.00 so f_det is equally appropriate.
* The equilibrium ``cubic griddata`` is the per-slice bottleneck. For
  ``--diag km9`` (real-cell path only) expect ~3-6 s per slice; --stride
  thins the time grid if you do not need every slice.
"""

from __future__ import annotations

import argparse
import os
import sys
import time as _time

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_COMMON_DIR = os.path.join(_HERE, "common")
if _COMMON_DIR not in sys.path:
    sys.path.insert(0, _COMMON_DIR)

from los_common import (                                                # noqa: E402
    find_run_dir,
    read_time_grid,
    read_thermal_slice,
    read_los_file,
    EqCDF,
    EqPPF,
    los_file_detector_rate,
)


DIAGS = {
    "km9":  {"chord_kind": "horizontal",
             "los_file":   os.path.join(_HERE, "km9",  "KM9.los")},
    "km14": {"chord_kind": "vertical",
             "los_file":   os.path.join(_HERE, "km14", "KM3.los")},
}

CHANNELS = {
    "total": ("THNTX",    "total (DD+DT+...)"),
    "dd":    ("THNTX_DD", "DD (2.45 MeV)"),
    "dt":    ("THNTX_DT", "DT (14 MeV)"),
}


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Build a time-resolved LoS weight .npz for KM9 / KM14.")
    p.add_argument('pulse', type=int, help='JET pulse number')
    p.add_argument('runid', type=str, help='TRANSP run suffix (e.g. M29)')
    p.add_argument('--diag', choices=list(DIAGS), required=True,
                   help='Diagnostic to compute the LoS weight for.')
    p.add_argument('--channel', choices=list(CHANNELS), default='total',
                   help='Reference thermal channel to co-store (default: '
                        'total). The weight itself is channel-independent.')
    p.add_argument('--eq-source', choices=['cdf', 'ppf'], default='cdf',
                   help='Equilibrium source: TRANSP CDF (default) or JET PPF.')
    p.add_argument('--dda', default='eftp', help='[ppf] PPF DDA (default eftp)')
    p.add_argument('--uid', default='jetppf', help='[ppf] PPF UID')
    p.add_argument('--seq', type=int, default=0, help='[ppf] PPF sequence')
    p.add_argument('--data-dir', default=None,
                   help='Override the TRANSP data base dir.')
    p.add_argument('--los-file', default=None,
                   help='Override the diagnostic-default LOS cell file.')
    p.add_argument('-t', '--time', type=float, nargs='+', default=None,
                   metavar='T',
                   help='Restrict TIME3 slices. One value: nearest slice. '
                        'Two values t1 t2 (JET s): inclusive window.')
    p.add_argument('--stride', type=int, default=1,
                   help='Take every Nth selected TIME3 slice (default 1).')
    p.add_argument('--chord-kind', choices=['auto', 'vertical', 'horizontal'],
                   default=None,
                   help='Override the diagnostic-default chord kind.')
    p.add_argument('--no-subgrid', action='store_true',
                   help='Point binning instead of anti-aliased subgrid.')
    p.add_argument('--no-median', action='store_true',
                   help='Disable cosmetic running-median de-speckle of f_det.')
    p.add_argument('--output', default=None,
                   help='Output .npz path (default: '
                        '<repo>/src/tmp/<run_id>_<DIAG>_los_weights.npz).')
    return p.parse_args(argv)


def _select_trinds(t_jet_grid, times, stride):
    if times is None:
        sel = np.ones_like(t_jet_grid, dtype=bool)
    elif len(times) == 1:
        sel = np.zeros_like(t_jet_grid, dtype=bool)
        sel[int(np.argmin(np.abs(t_jet_grid - float(times[0]))))] = True
    elif len(times) == 2:
        t1, t2 = sorted(float(t) for t in times)
        sel = (t_jet_grid >= t1) & (t_jet_grid <= t2)
    else:
        raise SystemExit('-t/--time takes 1 or 2 values')
    trinds = np.where(sel)[0]
    if stride > 1:
        trinds = trinds[::stride]
    return trinds


def main(argv=None):
    args = parse_args(argv)

    cfg = DIAGS[args.diag]
    los_file   = args.los_file   or cfg['los_file']
    chord_kind = args.chord_kind or cfg['chord_kind']

    run_id = f'{args.pulse}{args.runid}'
    run_dir = find_run_dir(args.pulse, args.runid, args.data_dir)
    cdf_path = run_dir / f'{run_id}.CDF'
    th_var, chan_label = CHANNELS[args.channel]

    print(f'Run dir   : {run_dir}')
    print(f'Diagnostic: {args.diag.upper()}  (chord_kind = {chord_kind})')
    print(f'LOS file  : {los_file}')
    print(f'Channel   : {chan_label}  ({th_var})')

    t_jet_grid = read_time_grid(cdf_path) + 40.0
    trinds = _select_trinds(t_jet_grid, args.time, args.stride)
    if trinds.size == 0:
        raise SystemExit('No TIME3 slices selected.')
    print(f'Slices    : {trinds.size}  '
          f'(t_JET = {t_jet_grid[trinds[0]]:.3f} .. '
          f'{t_jet_grid[trinds[-1]]:.3f} s, stride={args.stride})')

    cells = read_los_file(los_file)
    print(f'LOS cells : {cells["R"].size} from {los_file}')

    # Probe first slice to size the rhot grid.
    (_, x0, _, thntx_to_si, dvol_to_m3, unit_th) = read_thermal_slice(
        cdf_path, th_var, int(trinds[0]))
    n_rho = int(x0.size)
    n_t   = int(trinds.size)
    print(f'rhot grid : {n_rho} bins   THNTX units = [{unit_th}]')
    print()

    rhot_arr   = np.zeros((n_t, n_rho))
    dvol_arr   = np.zeros((n_t, n_rho))
    fdet_arr   = np.zeros((n_t, n_rho))
    thntx_arr  = np.zeros((n_t, n_rho))
    thlos_arr  = np.zeros((n_t, n_rho))
    times_arr  = np.zeros(n_t)
    rate_det   = np.zeros(n_t)
    rate_chord = np.zeros(n_t)
    rate_clos  = np.zeros(n_t)
    rho_50     = np.full(n_t, np.nan)
    rhot_min   = np.full(n_t, np.nan)
    rhot_crit  = np.full(n_t, np.nan)
    ok         = np.zeros(n_t, dtype=bool)

    t0 = _time.time()
    for i, ti in enumerate(trinds):
        ti = int(ti)
        t_jet = float(t_jet_grid[ti])
        try:
            (thntx_x, x_rhot, dvol_x, _, _, _) = read_thermal_slice(
                cdf_path, th_var, ti)
            if args.eq_source == 'ppf':
                eqs = EqPPF(args.pulse, args.dda, args.uid, args.seq, t_jet)
            else:
                eqs = EqCDF(cdf_path, t_jet)

            order = np.argsort(x_rhot)
            xs           = np.asarray(x_rhot)[order]
            thntx_sorted = np.asarray(thntx_x)[order]
            dvol_m3      = np.asarray(dvol_x)[order] * dvol_to_m3

            losf = los_file_detector_rate(
                cells, eqs, xs, thntx_sorted, thntx_to_si, dvol_m3,
                subgrid=not args.no_subgrid, chord_kind=chord_kind,
                median=not args.no_median,
            )

            rhot_arr[i]  = xs
            dvol_arr[i]  = dvol_m3
            fdet_arr[i]  = losf['f_det']
            thntx_arr[i] = thntx_sorted * thntx_to_si
            thlos_arr[i] = losf['thkm14_det_si']   # = THNTX * f_det in SI
            times_arr[i] = t_jet
            rate_det[i]   = losf['rate_det']
            rate_chord[i] = losf['rate_chord']
            rate_clos[i]  = losf['rate_from_profile']
            rho_50[i]     = losf['rho_med']
            if losf.get('rhot_min')  is not None: rhot_min[i]  = losf['rhot_min']
            if losf.get('rhot_crit') is not None: rhot_crit[i] = losf['rhot_crit']
            ok[i] = True

            elapsed = _time.time() - t0
            eta = elapsed / (i + 1) * (n_t - i - 1)
            print(f'  [{i + 1:3d}/{n_t:3d}] t_JET={t_jet:7.4f}s  '
                  f'Rate_det={rate_det[i]:.3e} n/s  '
                  f'rho_50={rho_50[i]:.4f}  '
                  f'(ETA {eta:5.0f}s)')

        except Exception as exc:
            times_arr[i] = t_jet
            print(f'  [{i + 1:3d}/{n_t:3d}] t_JET={t_jet:7.4f}s  '
                  f'SKIPPED ({type(exc).__name__}: {exc})')

    if not ok.any():
        raise SystemExit('All slices failed; no output written.')

    if args.output:
        out_path = args.output
    else:
        out_dir = os.path.abspath(os.path.join(_HERE, '..', 'tmp'))
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(
            out_dir, f'{run_id}_{args.diag.upper()}_los_weights.npz')

    np.savez_compressed(
        out_path,
        rhot=rhot_arr,
        dvol_m3=dvol_arr,
        f_det=fdet_arr,
        thntx_si=thntx_arr,
        thlos_det_si=thlos_arr,
        times_jet=times_arr,
        trinds=np.asarray(trinds, dtype=np.int32),
        rate_det=rate_det,
        rate_chord=rate_chord,
        rate_closure=rate_clos,
        rho_50=rho_50,
        rhot_min=rhot_min,
        rhot_crit=rhot_crit,
        ok=ok,
        pulse=np.int32(args.pulse),
        runid=np.array(args.runid),
        diagnostic=np.array(args.diag),
        channel=np.array(args.channel),
        th_var=np.array(th_var),
        eq_source=np.array(args.eq_source),
        chord_kind=np.array(chord_kind),
        los_file=np.array(los_file),
    )

    print()
    print(f'Wrote LoS weight file -> {out_path}')
    print(f'  {ok.sum()}/{n_t} slices succeeded.')
    print('  Apply in jet_tritium:')
    print('    z = np.load(path)')
    print('    rate_t[i] = np.sum(eps_on_rhot_i * z["f_det"][i] * z["dvol_m3"][i])')
    return 0


if __name__ == '__main__':
    sys.exit(main())
