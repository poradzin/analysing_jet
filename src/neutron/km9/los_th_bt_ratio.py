#!/usr/bin/env python
"""
los_th_bt_ratio.py  (KM9 / MPRu)
--------------------------------

Thermal-to-beam-target (TH/BT) neutron ratio as seen by the KM9 (MPRu)
line of sight, using the **real LOS cell cloud** (``KM9.los``).

This is the KM9 sibling of ``km14/los_th_bt_ratio.py``. KM9 is a *horizontal*
chord that crosses the vessel twice, so the idealised vertical "box chord"
machinery of the KM14 script (``Rmin/Rmax/wtor`` rates, ``chord_integral``, and
the **point-detector** BT-anisotropy approximation) has no analogue here and is
dropped -- exactly as ``km9/los_thermal_rate.py`` drops the box from its KM14
counterpart. What remains is the geometry-agnostic **real-cell detector path**:

  * f_det(rhot) = C_bin / DVOL -- the per-shell detector coupling from the real
    KM9 cells (purely geometric, emissivity-independent), via
    ``los_common.los_file_detector_rate``.
  * Apply the **same** f_det to both the thermal (THNTX) and the beam-target
    (flux-averaged BTN) emissivity (f_det is the etendue/geometry coupling; in
    v1 the BT emission is treated as **isotropic**, so TH and BT share it), and
    form the detector-weighted cumulative

        (TH/BT)_LOS = sum(TH * f_det * DVOL) / sum(BT * f_det * DVOL)
                    = Rate_det(TH) / Rate_det(BT)

    alongside the local ratio TH(rhot)/BT(rhot) and the whole-plasma 0D TH/BT.

BT emission anisotropy (the KM14 ``--bt-aniso`` directional g-factor) is **not**
included in v1 -- for a horizontal chord it is a separate, heavier calculation
(the per-cell ``(u,v,w)`` versors are near-horizontal). Deferred.

Time selection follows the fast-ion output: ``--idx X`` picks the
``_fi_X``/``_neut_X`` file written at ``OUTTIM(X)``; the thermal profile and
equilibrium are taken at the matching TIME3 slice.

Usage
-----
    python los_th_bt_ratio.py 104614 M29 --idx 2 --plot
    python los_th_bt_ratio.py 104614 M29 --idx 2 --channel dd --save out.png

The ``key=value`` form (e.g. ``channel=dd``) is also accepted.
"""

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# Shared LOS / equilibrium / binning library (src/neutron/common/).
_COMMON_DIR = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "../common"))
if _COMMON_DIR not in sys.path:
    sys.path.insert(0, _COMMON_DIR)
# bt_zone_integrator (the NUBEAM _fi/_neut reader) lives in km14/; it is shared
# NUBEAM infrastructure, so add that directory too. (A future refactor could
# move it into common/ -- see CONTEXT Stage B follow-ups.)
_KM14_DIR = os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "../km14"))
if _KM14_DIR not in sys.path:
    sys.path.insert(0, _KM14_DIR)

from los_common import (
    find_run_dir,
    read_thermal_slice,
    read_scalar_totals,
    read_los_file,
    flux_avg_profile,
    interp_flux_to,
    EqCDF,
    los_file_detector_rate,
)
import bt_zone_integrator as bzi

# Default real KM9 LOS cell file sits next to this script.
DEFAULT_LOS_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "KM9.los")

# (thermal profile var, list of beam-target _neut keys to sum, label).
# Empty bt-key list => sum every BT component present in _neut.
CHANNELS = {
    "total": ("THNTX",    [],      "total (DD+DT+...)"),
    "dd":    ("THNTX_DD", ["DD"],  "DD (2.45 MeV)"),
    "dt":    ("THNTX_DT", ["DT"],  "DT (14 MeV)"),
}


# -----------------------------------------------------------------------------
# CLI helpers
# -----------------------------------------------------------------------------
def _preprocess_argv(argv):
    out = []
    for tok in argv:
        if tok.startswith('-') or '=' not in tok:
            out.append(tok)
            continue
        k, v = tok.split('=', 1)
        out.extend([f'--{k}', v.strip().strip("'").strip('"')])
    return out


def parse_args(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    argv = _preprocess_argv(argv)
    p = argparse.ArgumentParser(
        description='KM9 (MPRu) real-LOS thermal/beam-target neutron ratio.')
    p.add_argument('pulse', type=int, help='JET pulse number')
    p.add_argument('runid', type=str, help='TRANSP run suffix (e.g. M29)')
    p.add_argument('--idx', type=int, default=None,
                   help='Fast-ion output index (1-based; selects '
                        '_fi_<idx>/_neut_<idx>). Default: first available.')
    p.add_argument('--data-dir', default=None,
                   help='Base data dir (<base>/<pulse>/<runid>); falls back to '
                        '~/jet/data then the heimdall results tree.')
    p.add_argument('--channel', choices=list(CHANNELS), default='total',
                   help='Channel: total (default), dd, dt.')
    p.add_argument('--los-file', default=DEFAULT_LOS_FILE,
                   help=f'Real KM9 LOS cell file (default: {DEFAULT_LOS_FILE}).')
    p.add_argument('--chord-kind', choices=['auto', 'horizontal', 'vertical'],
                   default='auto',
                   help='Chord classification for f_det (default auto -> '
                        'horizontal for KM9).')
    p.add_argument('--no-subgrid', action='store_true',
                   help='Disable subgrid (anti-aliased) cell->shell binning.')
    p.add_argument('--no-median', action='store_true',
                   help='Disable the cosmetic running-median de-speckle of f_det.')
    p.add_argument('--plot', action='store_true', help='Show diagnostic plots.')
    p.add_argument('--save', default=None, help='Save the figure to this path.')
    p.add_argument('--no-plot', action='store_true',
                   help='Suppress the figure even with --save.')
    return p.parse_args(argv)


# -----------------------------------------------------------------------------
# Diagnostics (2x3)
# -----------------------------------------------------------------------------
def diagnostic_plots(res, cells, save=None):
    """KM9 TH/BT 2x3 diagnostic:
        (0,0) poloidal (R, Z) LOS cells coloured by log10 C over the LCFS
        (0,1) top-of-machine (x-y) view vs the R=Rmag axis circle
        (0,2) f_det(rhot) detector coupling
        (1,0) local TH/BT(rhot) and cumulative TH/BT (C-weighted vs full plasma)
        (1,1) per-shell detector-weighted emission TH*f_det*DVOL & BT*f_det*DVOL
        (1,2) unused (removed)
    """
    xs = res['xs']
    fig, axes = plt.subplots(2, 3, figsize=(18.5, 10.5))
    fig.suptitle(
        f'KM9 LOS TH/BT  -  pulse {res["pulse"]}  TRANSP {res["runid"]}  '
        f'{res["chan_label"]}\n'
        f'idx {res["idx"]}  t_TRANSP={res["t_transp"]:.3f}s  '
        f't_EQ={res["t_eq_jet"]:.3f}s  ({res["eq_label"]})'
    )

    Rc = np.asarray(res['R_cells']); Zc = np.asarray(res['Z_cells'])
    Cc = np.asarray(res['C_cells']); ins = np.asarray(res['inside_cells'])
    xcell = np.asarray(cells['x']); ycell = np.asarray(cells['y'])
    Rb, Zb = res['Rb'], res['Zb']
    Rmag, Zmag = res['Rmag'], res['Zmag']
    idxc = np.where(Cc > 0)[0]
    if idxc.size > 40000:
        idxc = idxc[np.linspace(0, idxc.size - 1, 40000).astype(int)]
    out = idxc[~ins[idxc]]
    inn = idxc[ins[idxc]]

    # ---- (0,0) poloidal ---------------------------------------------------
    ax = axes[0, 0]
    ax.plot(Rb, Zb, 'b-', lw=1.4, label='LCFS')
    ax.scatter(Rc[out], Zc[out], c='0.8', s=2, alpha=0.4, linewidths=0)
    sc = ax.scatter(Rc[inn], Zc[inn], c=np.log10(Cc[inn]), s=4, cmap='plasma',
                    alpha=0.7, linewidths=0, label='cells inside LCFS')
    fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04).set_label(
        r'$\log_{10}$ C  [m$^3$]')
    ax.plot([Rmag], [Zmag], 'r+', ms=10, label=f'axis ({Rmag:.2f},{Zmag:.2f})')
    ax.set_xlabel('R [m]'); ax.set_ylabel('Z [m]')
    ax.set_title('Poloidal (R, Z): KM9 LOS cells + LCFS')
    ax.legend(loc='best', fontsize=8); ax.grid(True, ls=':', lw=0.5)
    ax.set_aspect('equal', adjustable='box')

    # ---- (0,1) top view ---------------------------------------------------
    ax = axes[0, 1]
    th = np.linspace(0.0, 2.0 * np.pi, 361)
    ax.scatter(xcell[out], ycell[out], c='0.8', s=2, alpha=0.4, linewidths=0)
    sc = ax.scatter(xcell[inn], ycell[inn], c=np.log10(Cc[inn]), s=4,
                    cmap='plasma', alpha=0.7, linewidths=0)
    fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04).set_label(
        r'$\log_{10}$ C  [m$^3$]')
    ax.plot(Rmag * np.cos(th), Rmag * np.sin(th), 'r-', lw=1.4,
            label=f'magnetic axis R={Rmag:.2f} m')
    ax.plot(Rb.min() * np.cos(th), Rb.min() * np.sin(th), 'b--', lw=0.9,
            label=f'LCFS R_in={Rb.min():.2f}')
    ax.plot(Rb.max() * np.cos(th), Rb.max() * np.sin(th), 'b--', lw=0.9,
            label=f'LCFS R_out={Rb.max():.2f}')
    ax.plot([0.0], [0.0], 'k+', ms=9, label='machine axis')
    ax.set_xlabel('x [m] (toroidal)'); ax.set_ylabel('y [m] (radial)')
    ax.set_title('Top view (x-y): LOS vs magnetic-axis circle')
    ax.legend(loc='upper right', fontsize=7); ax.grid(True, ls=':', lw=0.5)
    ax.set_aspect('equal', adjustable='box')

    # ---- (0,2) f_det ------------------------------------------------------
    ax = axes[0, 2]
    ax.plot(xs, res['f_det'], 'm-', lw=1.4, label='f_det (real LOS)')
    ax.fill_between(xs, 0.0, res['f_det'], color='m', alpha=0.15)
    if res['rhot_min'] is not None and np.isfinite(res['rhot_min']):
        ax.axvline(res['rhot_min'], color='c', ls=':', lw=1.0,
                   label=f'min sampled rhot (resolution floor) = {res["rhot_min"]:.4f}')
    ax.set_xlabel('rhot'); ax.set_ylabel(r'$f_{det}=C_{bin}/DVOL$')
    ax.set_title('Detector-coupling weight'); ax.set_xlim(0, 1)
    ax.grid(True, ls=':', lw=0.5); ax.legend(loc='best', fontsize=8)

    # ---- (1,0) TH/BT vs rhot ---------------------------------------------
    ax = axes[1, 0]
    ax.plot(xs, res['ratio_local'], 'k-', lw=1.4, label='local TH/BT(rhot)')
    ax.plot(xs, res['ratio_cum_det'], 'm-', lw=1.6,
            label=f'cum. C-weighted -> {res["thbt_los"]:.3f}')
    ax.plot(xs, res['ratio_cum_plasma'], color='0.5', ls='--', lw=1.3,
            label=f'cum. full plasma -> {res["thbt_plasma"]:.3f}')
    if np.isfinite(res['rho_med_th']):
        ax.axvline(res['rho_med_th'], color='c', ls='-.', lw=0.9,
                   label=f'TH rho_50={res["rho_med_th"]:.3f}')
    ax.set_xlabel('rhot'); ax.set_ylabel('TH / BT')
    ax.set_title('TH/BT vs rhot'); ax.set_xlim(0, 1)
    ax.set_ylim(bottom=0.0)
    ax.grid(True, ls=':', lw=0.5); ax.legend(loc='best', fontsize=8)

    # ---- (1,1) per-shell detector-weighted emission -----------------------
    ax = axes[1, 1]
    ax.plot(xs, res['shell_th_det'], 'b-', lw=1.4,
            label=r'TH$\cdot f_{det}\cdot$DVOL')
    ax.fill_between(xs, 0.0, res['shell_th_det'], color='b', alpha=0.12)
    ax.plot(xs, res['shell_bt_det'], 'r-', lw=1.4,
            label=r'BT$\cdot f_{det}\cdot$DVOL')
    ax.fill_between(xs, 0.0, res['shell_bt_det'], color='r', alpha=0.12)
    ax.set_xlabel('rhot'); ax.set_ylabel('detector-weighted rate [n/s per shell]')
    ax.set_title(r'Per-shell detector signal  ($\int$: TH=%.2e, BT=%.2e n/s)'
                 % (res['rate_det_th'], res['rate_det_bt']))
    ax.set_xlim(0, 1); ax.set_ylim(bottom=0.0)
    ax.grid(True, ls=':', lw=0.5); ax.legend(loc='best', fontsize=8)

    fig.delaxes(axes[1, 2])
    plt.tight_layout()
    if save:
        fig.savefig(save, dpi=130)
        print(f'Saved figure to {save}')
    plt.show()


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main(argv=None):
    args = parse_args(argv)
    run_id = f"{args.pulse}{args.runid}"
    run_dir = find_run_dir(args.pulse, args.runid, args.data_dir)
    cdf_path = run_dir / f"{run_id}.CDF"
    th_var, bt_keys_req, chan_label = CHANNELS[args.channel]

    idxs = bzi.list_fbm_indices(run_dir, run_id)
    if not idxs:
        print(f"No {run_id}_fi_*.cdf in {run_dir}", file=sys.stderr)
        return 1
    idx = args.idx if args.idx is not None else idxs[0]
    fi_path = run_dir / f"{run_id}_fi_{idx}.cdf"
    neut_path = run_dir / f"{run_id}_neut_{idx}.cdf"

    print(f'Run dir : {run_dir}')
    print(f'Channel : {chan_label}   thermal={th_var}')

    # ------- BT zone data + flux-averaged BT profile --------------------
    fi = bzi.read_fi_distribution(fi_path)
    neut = bzi.read_neut_rates(neut_path)
    available = [k for k in ("DD", "DT", "TT", "TD") if k in neut]
    bt_keys = bt_keys_req if bt_keys_req else available
    missing = [k for k in bt_keys if k not in neut]
    if missing:
        print(f"Beam-target key(s) {missing} not in {neut_path}", file=sys.stderr)
        return 1
    bt_zone_sum = sum(neut[k] for k in bt_keys)
    x_bt, bt_prof = flux_avg_profile(bt_zone_sum, fi["x2d"], fi["bmvol"])

    # ------- Equilibrium (CDF) at the fast-ion output time --------------
    eqs = EqCDF(cdf_path, fi["time"] + 40.0)   # EqCDF wants JET time
    Rb, Zb = eqs.lcfs()
    print(f'FBM idx : {idx}   t_TRANSP = {fi["time"]:.4f} s   '
          f'(TIME3[{eqs.tind}] = {eqs.t_eq_jet - 40.0:.4f} s)')
    print(f'BT keys : BTN[{"+".join(bt_keys)}]   '
          f'Rmag={eqs.Rmag:.3f} m, Zmag={eqs.Zmag:.3f} m')

    # ------- Thermal profile + DVOL on the TRANSP X grid ----------------
    (thntx_x, x_rhot, dvol_x, th_to_si,
     dvol_to_m3, unit_th) = read_thermal_slice(cdf_path, th_var, eqs.tind)
    order = np.argsort(x_rhot)
    xs = np.asarray(x_rhot)[order]
    th_native = np.asarray(thntx_x)[order]                 # native CGS
    dvol_m3 = np.asarray(dvol_x)[order] * dvol_to_m3
    th_si = th_native * th_to_si                           # n/m3/s
    bt_native = interp_flux_to(x_bt, bt_prof, xs)          # 1/cm3/s on xs
    bt_si = bt_native * 1.0e6                              # n/m3/s

    # ------- Real KM9 LOS cells: f_det (geometric) + TH/BT detector rates --
    cells = read_los_file(args.los_file)
    print(f'Loaded {cells["R"].size} cells from {args.los_file}')
    common = dict(subgrid=not args.no_subgrid, chord_kind=args.chord_kind,
                  median=not args.no_median)
    fTH = los_file_detector_rate(cells, eqs, xs, th_native, th_to_si,
                                 dvol_m3, **common)
    fBT = los_file_detector_rate(cells, eqs, xs, bt_native, 1.0e6,
                                 dvol_m3, **common)
    f_det = fTH['f_det']
    # f_det is purely geometric, so the TH- and BT-driven calls must agree.
    assert np.allclose(f_det, fBT['f_det']), 'f_det mismatch (should be geometric)'

    # ------- Ratios -----------------------------------------------------
    shell_th_det = th_si * f_det * dvol_m3                 # n/s per shell
    shell_bt_det = bt_si * f_det * dvol_m3
    cum_th_det = np.cumsum(shell_th_det)
    cum_bt_det = np.cumsum(shell_bt_det)
    cum_th_plasma = np.cumsum(th_si * dvol_m3)
    cum_bt_plasma = np.cumsum(bt_si * dvol_m3)
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio_local = np.where(bt_si > 0, th_si / bt_si, np.nan)
        ratio_cum_det = np.where(cum_bt_det > 0, cum_th_det / cum_bt_det, np.nan)
        ratio_cum_plasma = np.where(cum_bt_plasma > 0,
                                    cum_th_plasma / cum_bt_plasma, np.nan)
    thbt_los = fTH['rate_det'] / fBT['rate_det']           # = ratio_cum_det[-1]
    thbt_plasma = float(ratio_cum_plasma[-1])

    # 0D whole-plasma totals for cross-check.
    th_plasma = float(np.sum(th_si * dvol_m3))
    bt_plasma = float(np.sum(bt_zone_sum * fi["bmvol"]))
    totals = read_scalar_totals(cdf_path, eqs.tind)
    btnts = sum(totals.get('BTNTS_' + k, 0.0) for k in bt_keys)

    # ------- Report -----------------------------------------------------
    print()
    print('================  KM9 real-LOS TH / BT  ==================')
    print(f'  cells inside LCFS : {fTH["n_inside"]} / {cells["R"].size} '
          f'({100.0 * fTH["n_inside"] / cells["R"].size:.1f} %)')
    print(f'  chord_kind        : {fTH["chord_kind"]}  '
          f'(min sampled rhot = '
          f'{fTH["rhot_min"]:.4f})' if fTH["rhot_min"] is not None
          else f'  chord_kind        : {fTH["chord_kind"]}')
    print(f'  Rate_det(TH)      = {fTH["rate_det"]:.4e} n/s   '
          f'(closure {abs(fTH["rate_from_profile"]-fTH["rate_det"])/fTH["rate_det"]:.1e})')
    print(f'  Rate_det(BT)      = {fBT["rate_det"]:.4e} n/s   '
          f'(closure {abs(fBT["rate_from_profile"]-fBT["rate_det"])/fBT["rate_det"]:.1e})')
    print(f'  >> (TH/BT)_LOS    = {thbt_los:.4f}    (BT/TH = {1.0/thbt_los:.4f})')
    print(f'     cumulative profile endpoint check = {ratio_cum_det[-1]:.4f}')
    print('  ---- cross-checks (chord geometry cancels) ----')
    print(f'  (TH/BT) whole plasma (0D) = {th_plasma / bt_plasma:.4f}')
    print(f'     core enhancement LOS vs plasma = '
          f'{thbt_los / (th_plasma / bt_plasma):.3f}x')
    print(f'  TH signal-median rho_50 = {fTH["rho_med"]:.4f}   '
          f'BT signal-median rho_50 = {fBT["rho_med"]:.4f}')
    print('==========================================================')
    print(f'  TH whole-plasma = {th_plasma:.4e} n/s')
    print(f'  BT whole-plasma = {bt_plasma:.4e} n/s   '
          f'(0D BTNTS[{"+".join(bt_keys)}] = {btnts:.4e})')

    res = dict(
        pulse=args.pulse, runid=args.runid, idx=idx, chan_label=chan_label,
        t_transp=fi["time"], t_eq_jet=eqs.t_eq_jet, eq_label=eqs.label,
        xs=xs, f_det=f_det, rhot_min=fTH['rhot_min'],
        ratio_local=ratio_local, ratio_cum_det=ratio_cum_det,
        ratio_cum_plasma=ratio_cum_plasma,
        shell_th_det=shell_th_det, shell_bt_det=shell_bt_det,
        rate_det_th=fTH['rate_det'], rate_det_bt=fBT['rate_det'],
        thbt_los=thbt_los, thbt_plasma=thbt_plasma,
        rho_med_th=fTH['rho_med'], rho_med_bt=fBT['rho_med'],
        R_cells=fTH['R_cells'], Z_cells=fTH['Z_cells'],
        C_cells=fTH['C_cells'], inside_cells=fTH['inside_cells'],
        Rb=Rb, Zb=Zb, Rmag=eqs.Rmag, Zmag=eqs.Zmag,
    )

    if (args.plot or args.save) and not args.no_plot:
        diagnostic_plots(res, cells, save=args.save)
    return 0


if __name__ == '__main__':
    sys.exit(main())
