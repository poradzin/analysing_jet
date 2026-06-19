#!/usr/bin/env python
"""
los_thermal_rate.py  (KM9 / MPRu)
---------------------------------

Compute the thermal-neutron rate reaching the KM9 (MPRu) detector from the
volume swept by its real line of sight (LOS).

Geometry
~~~~~~~~
KM9 is a **horizontal** chord that crosses the JET vessel twice (the LOS
file ``KM9.los`` has |x| up to ~3.9 m in the tangential/toroidal direction,
versor ``w ~ 0.05``). This is fundamentally different from KM14's near-vertical
pencil chord, so the idealised "vertical box" framework used by
``km14/los_thermal_rate.py`` (fixed ``R in [Rmin, Rmax]``, full Z extent,
toroidal width ``w_tor``, equivalent ``rho_bnd``, A/B/C decomposition) has no
physical analogue here.

This script therefore uses **only the real-cell detector-rate path** shared in
``los_common.los_file_detector_rate``, which is geometry-agnostic. Each LOS
cell carries its own volume ``V`` and etendue weight ``C`` (with
``C = V * Omega/4pi``, ``Omega`` the detector solid angle seen from the cell),
so for the isotropic thermal source

    Rate_chord = sum_cells eps * V      (bare real-chord emission)
    Rate_det   = sum_cells eps * C      <== thermal rate reaching the detector

The detector response is binned onto the TRANSP flux shells to give the
LOS detector-coupling weight ``f_det(rhot) = C_bin / DVOL`` and the
LOS-weighted emissivity ``THLOS_det = THNTX * f_det`` with closure
``sum(THNTX * DVOL * f_det) == Rate_det``. The signal-median radius
``rho_50`` (50 % of ``Rate_det`` enclosed) is the detector-weighted analogue
of KM14's ``rho_bnd``.

The chord is classified as **horizontal** (``chord_kind="auto"`` infers it from
the cell versor ``median(|w|) <~ 0.05``), so the enclosed-shell flattening and
cosmetic median used for KM14's vertical chord are skipped -- a horizontal
chord that wraps the torus tangentially never fully sweeps a flux surface, so
there is no "enclosed shell" plateau.

Equilibrium source
~~~~~~~~~~~~~~~~~~
Same dual interface as the KM14 script: ``--eq-source cdf`` (default,
self-contained on the TRANSP main CDF, no ``ppf``) or ``--eq-source ppf``
(JET PPF EFTP, imports ``ppf``/``profiles``/``change_rho`` lazily).

Usage
-----
    # default: self-contained CDF equilibrium, real KM9 LOS cells
    python los_thermal_rate.py 104614 M29 -t 53.5268 --plot

    # PPF equilibrium (Heimdall / freia)
    python los_thermal_rate.py 104614 M29 --eq-source ppf \
            --dda eftp --uid gszepesi --seq 405 -t 53.5268

The convenience ``key=value`` form (e.g. ``dda='eftp'``) is also accepted.
"""

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# Shared LOS / equilibrium / binning library lives in src/neutron/common/, one
# level up from km9/; inject it onto sys.path so this script can be run
# directly (no package install).
_COMMON_DIR = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "../common"))
if _COMMON_DIR not in sys.path:
    sys.path.insert(0, _COMMON_DIR)

from los_common import (
    find_run_dir,
    read_thermal_slice,
    read_los_file,
    resolve_time_selection,
    EqCDF,
    EqPPF,
    los_file_detector_rate,
)

# Default real KM9 LOS cell file sits next to this script.
DEFAULT_LOS_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "KM9.los")

# Thermal-emissivity CDF variable per channel.
CHANNELS = {
    "total": ("THNTX",     "total (DD+DT+...)"),
    "dd":    ("THNTX_DD",  "DD (2.45 MeV)"),
    "dt":    ("THNTX_DT",  "DT (14 MeV)"),
}


# -----------------------------------------------------------------------------
# CLI helpers
# -----------------------------------------------------------------------------
def _preprocess_argv(argv):
    """Rewrite ``key=value`` tokens as ``--key value`` for argparse."""
    out = []
    for tok in argv:
        if tok.startswith('-') or '=' not in tok:
            out.append(tok)
            continue
        k, v = tok.split('=', 1)
        v = v.strip().strip("'").strip('"')
        out.extend([f'--{k}', v])
    return out


def parse_args(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    argv = _preprocess_argv(argv)

    parser = argparse.ArgumentParser(
        description='Thermal-neutron rate reaching the KM9 (MPRu) detector '
                    'from its real line-of-sight cell cloud.',
    )
    parser.add_argument('pulse', type=int, help='JET pulse number')
    parser.add_argument('runid', type=str, help='TRANSP run suffix (e.g. M29)')
    parser.add_argument('--eq-source', choices=['cdf', 'ppf'], default='cdf',
                        help='Equilibrium source: self-contained TRANSP CDF '
                             '(default) or JET PPF (lazily imports ppf).')
    parser.add_argument('--channel', choices=list(CHANNELS), default='total',
                        help='Thermal channel: total=THNTX (default), '
                             'dd=THNTX_DD, dt=THNTX_DT.')
    parser.add_argument('--data-dir', default=None,
                        help='Base data dir for the run CDFs '
                             '(<base>/<pulse>/<runid>). Falls back to '
                             '~/jet/data then the heimdall results tree.')
    parser.add_argument('--dda', type=str, default='eftp',
                        help='[ppf mode] Equilibrium PPF DDA (default: eftp)')
    parser.add_argument('--uid', type=str, default='jetppf',
                        help='[ppf mode] Equilibrium PPF UID (default: jetppf)')
    parser.add_argument('--seq', type=int, default=0,
                        help='[ppf mode] Equilibrium PPF sequence (default: 0)')
    parser.add_argument('-t', '--time', type=float, nargs='+', default=None,
                        metavar='T',
                        help='Time in JET convention (s, t>40). One value: '
                             'snap to the nearest TRANSP output slice. Two '
                             'values t1 t2: average the TRANSP signals over '
                             '[t1, t2] (eq mapped to the window midpoint). '
                             'If omitted, the run midpoint is used. Mutually '
                             'exclusive with --idx.')
    parser.add_argument('--idx', type=int, default=None,
                        help='Fast-ion output index (1-based, as in '
                             '_fi_<idx>.cdf). Averages the TRANSP signals over '
                             'the namelist window [OUTTIM(idx)-AVGTIM, '
                             'OUTTIM(idx)] and maps the equilibrium to '
                             'OUTTIM(idx). Mutually exclusive with -t.')
    parser.add_argument('--los-file', default=DEFAULT_LOS_FILE,
                        help='Path to the real KM9 LOS cell file '
                             f'(default: {DEFAULT_LOS_FILE}).')
    parser.add_argument('--chord-kind', choices=['auto', 'horizontal', 'vertical'],
                        default='auto',
                        help='Chord classification controlling the f_det '
                             'enclosed-shell flattening (default: auto -> '
                             'horizontal for KM9 cells).')
    parser.add_argument('--no-subgrid', action='store_true',
                        help='Disable subgrid (anti-aliased) cell->shell '
                             'binning of f_det; use point binning instead.')
    parser.add_argument('--no-median', action='store_true',
                        help='Disable the cosmetic 5-bin running-median '
                             'de-speckle of f_det (show the raw binned f_det). '
                             'Rate_det/rho_50 are unchanged either way.')
    parser.add_argument('--plot', action='store_true',
                        help='Show diagnostic plots.')
    parser.add_argument('--plot-los', action='store_true',
                        help='Show the LOS geometry figure (plot_LoS.py): '
                             'poloidal, top, and side-elevation views of the '
                             'real LOS cells.')
    parser.add_argument('--save', action='store_true',
                        help='Save (rhot, THNTX, f_det, THLOS_det, DVOL) '
                             'profile to src/tmp/.')
    args = parser.parse_args(argv)
    if args.time is not None and len(args.time) not in (1, 2):
        parser.error('-t/--time takes 1 value (nearest slice) or 2 values '
                     '(t1 t2 averaging window)')
    if args.time is not None and args.idx is not None:
        parser.error('-t/--time and --idx are mutually exclusive')
    return args


# -----------------------------------------------------------------------------
# Diagnostics  (KM9-tailored 1x3; LOS geometry lives in plot_LoS.py / --plot-los)
# -----------------------------------------------------------------------------
def diagnostic_plots(losf, xs, thntx_si_prof, dvol_m3,
                     pulse, runid, time_jet, t_eq_jet, eq_label, chan_label):
    """KM9 1x3 analysis diagnostic (no geometry panels -- those are in
    ``plot_LoS.py``, shown with ``--plot-los``):
        (0) f_det(rhot) detector-coupling weight (single curve)
        (1) per-shell emission rate: THNTX*DVOL (whole plasma) vs
            THNTX*f_det*DVOL (detector-weighted; integral = Rate_det)
        (2) cumulative detector-signal fraction with rho_50 marked
    """
    fig, axes = plt.subplots(1, 3, figsize=(18.0, 5.6))
    fig.suptitle(
        f'KM9 LOS thermal-neutron rate  -  pulse {pulse}  TRANSP {runid}  '
        f'{chan_label}\n'
        f't_JET={time_jet:.3f}s  t_EQ={t_eq_jet:.3f}s  ({eq_label})'
    )

    # ---- (0) f_det(rhot) -------------------------------------------------
    ax1 = axes[0]
    fd = np.asarray(losf['f_det'])
    ax1.plot(xs, fd, 'm-', lw=1.4, label='f_det (real LOS)')
    ax1.fill_between(xs, 0.0, fd, color='m', alpha=0.15)
    rmin = losf.get('rhot_min')
    if rmin is not None and np.isfinite(rmin):
        ax1.axvline(rmin, color='c', ls=':', lw=1.0,
                    label=f'min sampled rhot (resolution floor) = {rmin:.4f}')
    if np.isfinite(losf['rho_med']):
        ax1.axvline(losf['rho_med'], color='k', ls='-.', lw=1.0,
                    label=f'rho_50 = {losf["rho_med"]:.4f}')
    ax1.set_xlabel('rhot')
    ax1.set_ylabel(r'$f_{det}(\rho_t) = C_{bin}/DVOL$')
    ax1.set_title('Detector-coupling weight')
    ax1.set_xlim(0.0, 1.0)
    ax1.grid(True, ls=':', lw=0.5)
    ax1.legend(loc='best', fontsize=8)

    # ---- (1) per-shell emission rate: THNTX*DVOL vs THNTX*f_det*DVOL ------
    # THNTX*DVOL is the whole-plasma thermal emission per flux shell [n/s];
    # THNTX*f_det*DVOL is the detector-weighted per-shell contribution whose
    # integral over rhot equals Rate_det. The detector curve is ~f_det (~1e-8)
    # times smaller, so it gets its own right-hand axis (unlike the KM14 box
    # where f<=1 keeps the two comparable on one axis); the shapes show how the
    # LOS+etendue reweight the plasma emission radially.
    ax3 = axes[1]
    dvol = np.asarray(dvol_m3)
    fd_p = np.asarray(losf['f_det'])
    emis_plasma = thntx_si_prof * dvol             # n/s per shell (whole plasma)
    emis_det = thntx_si_prof * fd_p * dvol         # n/s per shell (detector)
    l_a = ax3.plot(xs, emis_plasma, 'b-', lw=1.4,
                   label=r'THNTX$\cdot$DVOL (whole plasma)')
    ax3.fill_between(xs, 0.0, emis_plasma, color='b', alpha=0.10)
    ax3.set_xlabel('rhot')
    ax3.set_ylabel(r'THNTX$\cdot$DVOL  [n/s per shell]', color='b')
    ax3.tick_params(axis='y', labelcolor='b')
    ax3.set_xlim(0.0, 1.0)
    ax3.grid(True, ls=':', lw=0.5)
    ax3b = ax3.twinx()
    l_b = ax3b.plot(xs, emis_det, 'm-', lw=1.4,
                    label=r'THNTX$\cdot f_{det}\cdot$DVOL (detector)')
    ax3b.fill_between(xs, 0.0, emis_det, color='m', alpha=0.18)
    ax3b.set_ylabel(r'THNTX$\cdot f_{det}\cdot$DVOL  [n/s per shell]', color='m')
    ax3b.tick_params(axis='y', labelcolor='m')
    ax3b.set_ylim(bottom=0.0)
    ax3.set_ylim(bottom=0.0)
    lns = list(l_a) + list(l_b)
    if np.isfinite(losf['rho_med']):
        vl = ax3.axvline(losf['rho_med'], color='k', ls='-.', lw=1.0,
                         label=f'rho_50 = {losf["rho_med"]:.4f}')
        lns += [vl]
    ax3.set_title(r'Per-shell emission   ($\int$detector = Rate$_{det}$ = '
                  f'{losf["rate_det"]:.2e} n/s)')
    ax3.legend(lns, [l.get_label() for l in lns], loc='upper right', fontsize=7)

    # ---- (2) cumulative detector-signal fraction -------------------------
    ax2 = axes[2]
    ax2.plot(losf['cum_rhot'], losf['cum_frac'], 'b-', lw=1.4,
             label='cum. detector signal')
    ax2.axhline(0.5, color='0.6', ls=':', lw=0.8)
    if np.isfinite(losf['rho_med']):
        ax2.axvline(losf['rho_med'], color='k', ls='-.', lw=1.0,
                    label=f'rho_50 = {losf["rho_med"]:.4f}')
    ax2.set_xlabel('rhot')
    ax2.set_ylabel(r'cumulative fraction of Rate$_{det}$')
    ax2.set_title('Where KM9 picks up its thermal signal')
    ax2.set_xlim(0.0, 1.0)
    ax2.set_ylim(0.0, 1.02)
    ax2.grid(True, ls=':', lw=0.5)
    ax2.legend(loc='best', fontsize=8)

    plt.tight_layout()
    plt.show()


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main(argv=None):
    args = parse_args(argv)

    run_id = f"{args.pulse}{args.runid}"
    run_dir = find_run_dir(args.pulse, args.runid, args.data_dir)
    cdf_path = run_dir / f"{run_id}.CDF"
    th_var, chan_label = CHANNELS[args.channel]
    print(f'Run dir : {run_dir}')
    print(f'Channel : {chan_label}  (thermal var = {th_var})')

    # ------- Time selection (-t / --idx) + thermal slice ----------------
    trinds, eq_ref_jet, time_jet, time_desc = resolve_time_selection(
        cdf_path, run_dir, run_id, times=args.time, idx=args.idx)
    print(f'Time selection: {time_desc}')

    (thntx_x, x_rhot, dvol_x, thntx_to_si,
     dvol_to_m3, unit_th) = read_thermal_slice(cdf_path, th_var, trinds)

    print(f'Averaged TRANSP slice(s): {trinds.tolist()}')
    print(f'  {th_var} units = [{unit_th}],  conversion to n/m3/s = '
          f'{thntx_to_si:.1e}')
    print(f'  {th_var}(rhot=0) = {thntx_x[0]:.3e} [{unit_th}]  '
          f'= {thntx_x[0] * thntx_to_si:.3e} n/m3/s')

    # ------- Equilibrium (single slice at eq_ref_jet) -------------------
    if args.eq_source == 'ppf':
        print(f'Equilibrium: PPF {args.dda}/{args.uid}/{args.seq} '
              f'@ {eq_ref_jet:.3f}s (importing ppf) ...')
        eqs = EqPPF(args.pulse, args.dda, args.uid, args.seq, eq_ref_jet)
    else:
        print(f'Equilibrium: TRANSP CDF @ {eq_ref_jet:.3f}s '
              f'(self-contained, no ppf) ...')
        eqs = EqCDF(cdf_path, eq_ref_jet)
    Rmag, Zmag = eqs.Rmag, eqs.Zmag
    print(f'Equilibrium slice: t = {eqs.t_eq_jet:.3f} s   [{eqs.label}]')
    print(f'  Rmag = {Rmag:.3f} m, Zmag = {Zmag:.3f} m')

    # ------- Sorted TRANSP profile arrays -------------------------------
    order = np.argsort(x_rhot)
    xs_sorted = np.asarray(x_rhot)[order]
    thntx_sorted = np.asarray(thntx_x)[order]
    dvol_m3 = np.asarray(dvol_x)[order] * dvol_to_m3
    thntx_si_profile = thntx_sorted * thntx_to_si      # n/m3/s

    # ------- Real KM9 LOS cell file (detector-weighted) -----------------
    print()
    print('=============  Real KM9 LOS cell file  =============')
    cells = read_los_file(args.los_file)
    n_cells = cells["R"].size
    print(f'  Loaded {n_cells} cells from {args.los_file}')
    print(f'  Cell R in [{cells["R"].min():.3f}, {cells["R"].max():.3f}] m, '
          f'Z in [{cells["Z"].min():.3f}, {cells["Z"].max():.3f}] m')
    print(f'  |x| <= {np.abs(cells["x"]).max():.3f} m (tangential), '
          f'median |w| = {np.median(np.abs(cells["w"])):.3f} '
          f'(versor vertical component)')

    losf = los_file_detector_rate(
        cells, eqs, xs_sorted, thntx_sorted, thntx_to_si, dvol_m3,
        subgrid=not args.no_subgrid, chord_kind=args.chord_kind,
        median=not args.no_median,
    )
    rmin = losf["rhot_min"]
    rmin_str = f'{rmin:.4f}' if rmin is not None else 'n/a'
    print(f'  chord_kind resolved -> {losf["chord_kind"]}  '
          f'(rhot_crit = {losf["rhot_crit"]}, '
          f'closest-approach rhot_min = {rmin_str})')

    # Open-question check (CONTEXT Stage B): the KM9 LOS reaches well outside
    # the LCFS and even outside the equilibrium computational box, so the
    # scattered griddata rhot may return NaN for far-out cells before the
    # inside-LCFS mask zeros them. Report the fraction.
    rhot_cells = np.asarray(losf['rhot_cells'])
    n_nan = int(np.count_nonzero(~np.isfinite(rhot_cells)))
    nan_inside = int(np.count_nonzero(
        ~np.isfinite(rhot_cells) & np.asarray(losf['inside_cells'])))
    print(f'  rhot NaN cells: {n_nan} / {n_cells} '
          f'({100.0 * n_nan / n_cells:.2f} %)  '
          f'-- of which inside LCFS: {nan_inside} (should be 0)')

    print()
    print(f'  Cells inside LCFS    : {losf["n_inside"]} / {n_cells} '
          f'({100.0 * losf["n_inside"] / n_cells:.1f} %)')
    print(f'  LOS volume (sum V)   = {losf["vol_tot"]:.4e} m^3  '
          f'(in-LCFS cells)')
    print(f'  LOS etendue (sum C)  = {losf["c_tot"]:.4e} m^3  '
          f'(= sum V*Omega/4pi)')
    print()
    print(f'  Rate_chord = sum eps*V = {losf["rate_chord"]:.4e} n/s  '
          f'(real chord, no solid angle)')
    print(f'  Rate_det   = sum eps*C = {losf["rate_det"]:.4e} n/s  '
          f'<== thermal rate reaching the KM9 detector')
    rel = (abs(losf["rate_from_profile"] - losf["rate_det"]) / losf["rate_det"]
           if losf["rate_det"] > 0 else float('nan'))
    print(f'  closure sum(THNTX*DVOL*f_det) = {losf["rate_from_profile"]:.4e} '
          f'n/s  (rel. residual {rel:.2e})   <== acceptance test')
    print()
    print('  --- where does KM9 see its thermal signal? ---')
    print(f'  signal-median radius rho_50 (50% of Rate_det inside) = '
          f'{losf["rho_med"]:.4f}')
    print('====================================================')

    # ------- Save -------------------------------------------------------
    if args.save:
        # Repo-local tmp/ (src/tmp), independent of where the CDFs live.
        src_dir = os.path.abspath(
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
        outdir = os.path.join(src_dir, 'tmp')
        os.makedirs(outdir, exist_ok=True)
        out_path = os.path.join(
            outdir,
            f'{run_id}_KM9_LOS_profile_{args.channel}_{args.eq_source}'
            f'_t{time_jet:.3f}s.txt',
        )
        cols = [xs_sorted, thntx_si_profile, np.asarray(losf['f_det']),
                np.asarray(losf['thkm14_det_si']), dvol_m3]
        col_hdr = ('    rhot           THNTX [n/m3/s]   f_det           '
                   'THLOS_det [n/m3/s]  DVOL [m3]')
        data = np.column_stack(cols)
        np.savetxt(
            out_path, data,
            header=(f'pulse {args.pulse}  TRANSP {args.runid}  '
                    f'channel {args.channel} ({th_var})  '
                    f't_JET = {time_jet:.3f} s\n'
                    f'equilibrium: {eqs.label}  t_EQ = {eqs.t_eq_jet:.3f} s\n'
                    f'real LOS (KM9.los): chord_kind = {losf["chord_kind"]}   '
                    f'Rate_det = {losf["rate_det"]:.4e} n/s   '
                    f'Rate_chord = {losf["rate_chord"]:.4e} n/s   '
                    f'rho_50 = {losf["rho_med"]:.4f}\n'
                    f'{col_hdr}'),
            fmt='%.8e',
        )
        print(f'Saved LOS profile to {out_path}')

    # ------- Plot -------------------------------------------------------
    if args.plot_los:
        import plot_LoS
        Rb, Zb = eqs.lcfs()
        plot_LoS.plot_los_geometry(
            cells, losf['inside_cells'], Rb, Zb, Rmag, Zmag,
            title=f'KM9 LOS geometry  -  pulse {args.pulse}  TRANSP {args.runid}'
                  f'   t_EQ={eqs.t_eq_jet:.3f}s',
            rhot_min=losf['rhot_min'])
    if args.plot:
        diagnostic_plots(losf, xs_sorted, thntx_si_profile, dvol_m3,
                         args.pulse, args.runid, time_jet,
                         eqs.t_eq_jet, eqs.label, chan_label)

    return 0


if __name__ == '__main__':
    sys.exit(main())
