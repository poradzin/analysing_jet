#!/usr/bin/env python
"""
los_thermal_rate.py
-------------------

Compute the total thermal-neutron rate emitted from the volume swept by the
KM14 vertical line of sight (LOS).

Geometry
~~~~~~~~
The KM14 LOS is a vertical chord (constant toroidal angle) located above
the JET vessel. In the (R, Z) poloidal plane it spans

    R in [R_MIN, R_MAX]   (default: 2.7 m .. 3.1 m)
    Z covering the plasma vertical extent (default: equilibrium Z box)

A fixed effective toroidal width ``W_TOR`` (default 0.4 m) closes the
chord cross section so the volume element is

    dV = W_TOR * dR * dZ

(this is the same convention used by ``estimate_outer_th_fraction.py``).

Equilibrium source
~~~~~~~~~~~~~~~~~~
By default the equilibrium (rhot(R, Z) map + LCFS) is taken **self-contained
from the TRANSP main CDF** (``--eq-source cdf``), exactly like
``los_th_bt_ratio.py``: psi(R, Z) from ``PSIRZ``, normalized by
``PSI0_TR``/``PLFLXA``, with rhot(psi_n) inverted from ``PLFLX`` vs ``XB``;
the LCFS is reconstructed from the time-resolved asymmetric boundary Fourier
moments (``RMCB*``/``RMSB*``/``YMCB*``/``YMSB*``). This needs no ``ppf`` and
therefore runs in the WSL dev env, on freia and on heimdall alike.

Passing ``--eq-source ppf`` instead uses a JET PPF equilibrium (default DDA
``eftp``) via ``profiles.Eq`` and ``change_rho`` -- those modules are imported
**lazily**, only when this option is selected, so the default CDF path does not
require ``ppf`` to be importable.

The thermal emissivity profile and DVOL are read directly from the TRANSP main
CDF in both modes.

Algorithm
~~~~~~~~~
1.  Map a dense (R, Z) chord grid to ``rhot`` using the selected equilibrium.
2.  Build the inside-LCFS mask as a point-in-polygon test against the closed
    LCFS boundary. This excludes the private flux region below the X-point,
    where ``psin < 1`` but TRANSP does not model emission.
3.  PCHIP-interpolate the TRANSP ``THNTX[_<chan>](X)`` profile (X = rhot) onto
    the dense rhot map; set ``THNTX = 0`` outside the LCFS.
4.  Convert ``THNTX`` from TRANSP CGS units (N/CM3/SEC) to SI (n/m3/s) and
    integrate using the trapezoidal rule:

        Rate_LOS [n/s] = W_TOR * integral over R, Z of THNTX(R, Z) dR dZ

Usage
-----
    # default: self-contained CDF equilibrium
    python los_thermal_rate.py 104614 M30 -t 53.5268 --plot

    # PPF equilibrium (Heimdall / freia, imports ppf lazily)
    python los_thermal_rate.py 104614 M29 --eq-source ppf \
            --dda eftp --uid gszepesi --seq 405 -t 53.5268 --plot

The convenience ``key=value`` form (e.g. ``dda='eftp'``) is also accepted.
"""

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt

# Shared LOS / equilibrium / binning library lives in src/neutron/common/, one
# level up from km14/; inject it onto sys.path so this script can be run
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
    thntx_on_grid,
    find_rho_bnd,
    los_shell_fraction,
    los_file_detector_rate,
)


# -----------------------------------------------------------------------------
# Geometric defaults for the KM14 LOS
# -----------------------------------------------------------------------------
R_MIN_DEFAULT = 2.60   # inner R limit of the chord [m] (real LOS footprint)
R_MAX_DEFAULT = 3.16   # outer R limit of the chord [m] (real LOS footprint)
W_TOR_DEFAULT = 0.40   # effective toroidal width of the chord [m]
NR_DEFAULT = 100
NZ_DEFAULT = 100
Z_MARGIN_DEFAULT = 0.0  # extra Z padding beyond LCFS [m]

# Thermal-emissivity CDF variable per channel.
CHANNELS = {
    "total": ("THNTX",    "total (DD+DT+...)"),
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
        description='Total thermal-neutron rate emitted from the KM14 '
                    'vertical LOS volume.',
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
    parser.add_argument('--Rmin', type=float, default=R_MIN_DEFAULT,
                        help=f'LOS inner R [m] (default {R_MIN_DEFAULT})')
    parser.add_argument('--Rmax', type=float, default=R_MAX_DEFAULT,
                        help=f'LOS outer R [m] (default {R_MAX_DEFAULT})')
    parser.add_argument('--wtor', type=float, default=W_TOR_DEFAULT,
                        help=f'Effective toroidal width of the LOS [m] '
                             f'(default {W_TOR_DEFAULT})')
    parser.add_argument('--nR', type=int, default=NR_DEFAULT,
                        help=f'Number of R samples (default {NR_DEFAULT})')
    parser.add_argument('--nZ', type=int, default=NZ_DEFAULT,
                        help=f'Number of Z samples (default {NZ_DEFAULT})')
    parser.add_argument('--zmargin', type=float, default=Z_MARGIN_DEFAULT,
                        help='Extra Z padding above/below the LCFS extent [m] '
                             f'(default {Z_MARGIN_DEFAULT})')
    parser.add_argument('--los-file', default=None,
                        help='Path to the real KM14 LOS cell file (_KM3.los). '
                             'When given, the detector-weighted thermal rate '
                             '(Rate_det = sum eps*C) and the C-weighted LOS '
                             'response profile are computed from the actual '
                             'cells, in addition to the idealised-chord box.')
    parser.add_argument('--plot', action='store_true',
                        help='Show diagnostic plots.')
    parser.add_argument('--save', action='store_true',
                        help='Save (rhot, THNTX, THKM14, f) profile to tmp/.')
    parser.add_argument('--no-subgrid', action='store_true',
                        help='Disable subgrid (R,Z)-cell rhot distribution '
                             'for f(rhot); use point binning instead. '
                             'Subgrid is on by default to suppress aliasing.')
    parser.add_argument('--no-median', action='store_true',
                        help='Disable the cosmetic running-median de-speckle of '
                             'the real-LOS f_det (--los-file); show raw f_det. '
                             'Rate_det/rho_50 are unchanged either way.')
    args = parser.parse_args(argv)
    if args.time is not None and len(args.time) not in (1, 2):
        parser.error('-t/--time takes 1 value (nearest slice) or 2 values '
                     '(t1 t2 averaging window)')
    if args.time is not None and args.idx is not None:
        parser.error('-t/--time and --idx are mutually exclusive')
    return args


# (read_time_grid / read_thermal_slice / read_los_file / equilibrium classes
#  EqCDF / EqPPF and their helpers _boundary_from_moments, _rhot_scatter,
#  _lcfs_contour all live in src/neutron/common/los_common.py; see the import
#  block at the top of this file.)


# -----------------------------------------------------------------------------
# LOS grid + emissivity (source-agnostic)
# -----------------------------------------------------------------------------
def build_rz_grid(z_bot, z_top, Rmin, Rmax, nR, nZ, z_margin):
    """Build a (nR x nZ) Cartesian grid covering the LOS region.

    The Z range is taken from the equilibrium computational box (the whole
    vessel), not from the LCFS extent: vacuum points are zeroed out by the
    inside-LCFS polygon mask.
    """
    R = np.linspace(Rmin, Rmax, nR)
    Z = np.linspace(z_bot - z_margin, z_top + z_margin, nZ)
    return R, Z


def integrate_rate(thntx_grid_si, R, Z, w_tor):
    """Trapezoidal integral of THNTX (in n/m3/s) * w_tor over (R, Z).

    Returns the rate in n/s and the dR-integrated linear emissivity
    (n/s per meter of Z) for diagnostics.
    """
    inner_R = np.trapezoid(thntx_grid_si, R, axis=0)        # over R -> shape (nZ,)
    rate = w_tor * np.trapezoid(inner_R, Z)                 # over Z -> scalar
    return float(rate), inner_R


def toroidal_rate(thntx_grid_si, R, Z, mask=None):
    """Toroidally-integrated emission rate over the LOS (R, Z) area.

    Replaces the linear toroidal width with the proper toroidal volume
    element 2*pi*R, so the result is the full-torus emission from the
    (R, Z) area covered by the LOS:

        Rate_tor [n/s] = integral over R, Z of  2*pi*R * THNTX(R, Z) dR dZ

    If ``mask`` is given (boolean, same shape as ``thntx_grid_si``), the
    integrand is zeroed outside the mask -- used for sub-region integrals
    like "inside LOS AND inside rho_bnd".
    """
    Rg, _ = np.meshgrid(R, Z, indexing='ij')
    integrand = thntx_grid_si * (2.0 * np.pi) * Rg
    if mask is not None:
        integrand = integrand * np.asarray(mask, dtype=float)
    inner = np.trapezoid(integrand, R, axis=0)
    return float(np.trapezoid(inner, Z))


# -----------------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------------
def diagnostic_plots(R, Z, rhot, inside, psin_dense, thntx_grid_si, inner_R,
                     Rb, Zb, native, pulse, runid, time_jet, t_eq_jet,
                     eq_label, chan_label, Rmag, Zmag, rho_bnd=None,
                     xs=None, thntx_si_prof=None,
                     thkm14_si_prof=None, f_prof=None, losf=None):
    """Six-panel (2x3) diagnostic:
        (0,0) rhot(R,Z) in LOS region + LCFS + rho_bnd
        (0,1) THNTX(R,Z) in LOS region + LCFS + rho_bnd
        (0,2) THNTX(rhot) and THKM14(rhot) profiles
        (1,0) Full poloidal LCFS with LOS box and rho_bnd surface
        (1,1) R-integrated emissivity vs Z
        (1,2) LOS weight function f(rhot)
    """
    Rg, Zg = np.meshgrid(R, Z, indexing='ij')
    rhot_plot = np.where(inside, rhot, np.nan)

    Rg_n, Zg_n, rhot_native = native

    rho_bnd_lbl = (f'rho_bnd = {rho_bnd:.4f}'
                   if rho_bnd is not None and np.isfinite(rho_bnd)
                   else 'rho_bnd: N/A')

    fig, axes = plt.subplots(2, 3, figsize=(17.0, 10.0))
    fig.suptitle(
        f'KM14 LOS thermal-neutron rate  -  pulse {pulse}  TRANSP {runid}  '
        f'{chan_label}\n'
        f't_JET={time_jet:.3f}s  t_EQ={t_eq_jet:.3f}s  ({eq_label})'
    )

    ax1 = axes[0, 0]
    pc1 = ax1.pcolormesh(Rg, Zg, rhot_plot, shading='auto', cmap='viridis')
    fig.colorbar(pc1, ax=ax1, label='rhot')
    ax1.plot(Rb, Zb, 'w-', lw=1.2, label='LCFS')
    if psin_dense is not None:
        ax1.contour(Rg, Zg, psin_dense, levels=[1.0],
                    colors='white', linewidths=0.8, linestyles='--')
    if rho_bnd is not None and np.isfinite(rho_bnd):
        ax1.contour(Rg, Zg, rhot, levels=[rho_bnd],
                    colors='magenta', linewidths=1.4)
    ax1.plot([Rmag], [Zmag], 'r+', ms=10, label='axis')
    ax1.set_xlim(np.nanmin(Rb)-0.1, np.nanmax(Rb)+0.1)
    ax1.set_ylim(np.nanmin(Zb)-0.1, np.nanmax(Zb)+0.1)
    ax1.set_xlabel('R [m]')
    ax1.set_ylabel('Z [m]')
    ax1.set_title('rhot on LOS grid')
    ax1.legend(loc='upper right', fontsize=8)
    ax1.set_aspect('equal', adjustable='box')

    ax2 = axes[0, 1]
    pc2 = ax2.pcolormesh(Rg, Zg, thntx_grid_si, shading='auto', cmap='inferno')
    fig.colorbar(pc2, ax=ax2, label='THNTX [n/m3/s]')
    ax2.plot(Rb, Zb, 'c-', lw=1.2, label='LCFS')
    if rho_bnd is not None and np.isfinite(rho_bnd):
        ax2.contour(Rg, Zg, rhot, levels=[rho_bnd],
                    colors='magenta', linewidths=1.4)
    ax2.plot([Rmag], [Zmag], 'w+', ms=10, label='axis')
    ax2.set_xlim(np.nanmin(Rb)-0.1, np.nanmax(Rb)+0.1)
    ax2.set_ylim(np.nanmin(Zb)-0.1, np.nanmax(Zb)+0.1)
    ax2.set_xlabel('R [m]')
    ax2.set_ylabel('Z [m]')
    ax2.set_title(f'THNTX(R, Z) inside LOS   ({rho_bnd_lbl})')
    ax2.legend(loc='upper right', fontsize=8)
    ax2.set_aspect('equal', adjustable='box')

    ax3 = axes[1, 0]
    ax3.plot(Rb, Zb, 'b-', lw=1.4, label='LCFS')
    # Real KM14 LOS cells (when --los-file given): scatter the actual footprint
    # coloured by etendue C (detector sensitivity), to compare against the
    # simplified box rectangle drawn below. Subsample for responsiveness.
    if losf is not None and 'R_cells' in losf:
        Rc = np.asarray(losf['R_cells'])
        Zc = np.asarray(losf['Z_cells'])
        Cc = np.asarray(losf['C_cells'])
        ins = np.asarray(losf['inside_cells'])
        idx = np.where(ins & (Cc > 0))[0]
        if idx.size > 25000:
            idx = idx[np.linspace(0, idx.size - 1, 25000).astype(int)]
        sc = ax3.scatter(Rc[idx], Zc[idx], c=np.log10(Cc[idx]), s=3,
                         cmap='plasma', alpha=0.55, linewidths=0,
                         label='real LOS cells')
        cb = fig.colorbar(sc, ax=ax3, fraction=0.046, pad=0.04)
        cb.set_label(r'$\log_{10}$ C  [m$^3$]  (etendue)')
    if rho_bnd is not None and np.isfinite(rho_bnd):
        ax3.contour(Rg_n, Zg_n, rhot_native, levels=[rho_bnd],
                    colors='magenta', linewidths=1.4)
        # Proxy handle for the legend: QuadContourSet.collections was removed
        # in matplotlib 3.8, so attach the label via an empty line instead.
        ax3.plot([], [], color='magenta', lw=1.4, label=rho_bnd_lbl)
    ax3.plot([Rmag], [Zmag], 'r+', ms=10, label=f'axis ({Rmag:.3f}, {Zmag:.3f})')
    los_box_R = [R[0], R[-1], R[-1], R[0], R[0]]
    los_box_Z = [Z[0], Z[0], Z[-1], Z[-1], Z[0]]
    ax3.plot(los_box_R, los_box_Z, 'r-', lw=1.2,
             label=f'simplified box R=[{R[0]:.2f},{R[-1]:.2f}]')
    ax3.axvline(0.5 * (R[0] + R[-1]), color='r', ls=':', lw=0.8,
                label=f'box centre R={0.5*(R[0]+R[-1]):.2f} m')
    ax3.set_xlabel('R [m]')
    ax3.set_ylabel('Z [m]')
    ax3.set_title('Poloidal cross-section: real LOS cells vs simplified box'
                  if losf is not None else
                  'Poloidal cross-section: LCFS + KM14 LOS')
    ax3.legend(loc='best', fontsize=8)
    ax3.grid(True, ls=':', lw=0.5)
    ax3.set_aspect('equal', adjustable='box')

    ax4 = axes[1, 1]
    ax4.plot(inner_R, Z, lw=1.4)
    ax4.set_xlabel(r'$\int$ THNTX dR  [n/m2/s]')
    ax4.set_ylabel('Z [m]')
    ax4.set_title('R-integrated emissivity vs Z')
    ax4.axhline(Zmag, color='r', ls='--', lw=0.8, label=f'Zmag={Zmag:.3f} m')
    ax4.grid(True, ls=':', lw=0.5)
    ax4.legend(loc='best', fontsize=8)

    # ---- (0, 2): emissivity profiles --------------------------------
    ax5 = axes[0, 2]
    if xs is not None and thntx_si_prof is not None:
        ax5.plot(xs, thntx_si_prof, 'b-', lw=1.4,
                 label='THNTX (TRANSP)')
        ax5.plot(xs, thkm14_si_prof, 'r-', lw=1.4,
                 label=r'THKM14 = THNTX$\cdot f$ (box)')
        if losf is not None:
            # Detector-weighted profile rescaled to the geometric THKM14 peak so
            # its radial *shape* is comparable on the same axis (f_det carries
            # the tiny Omega/4pi factor, so the raw curve would be invisible).
            tk = np.asarray(losf['thkm14_det_si'])
            scale = (np.nanmax(thkm14_si_prof) / np.nanmax(tk)
                     if np.nanmax(tk) > 0 else 1.0)
            ax5.plot(xs, tk * scale, 'm--', lw=1.4,
                     label=f'THKM14_det (x{scale:.1e}, real LOS)')
        ax5.set_xlabel('rhot')
        ax5.set_ylabel('Emissivity [n/m3/s]')
        ax5.set_title('Emissivity profiles vs rhot')
        if rho_bnd is not None and np.isfinite(rho_bnd):
            ax5.axvline(rho_bnd, color='m', ls=':', lw=1.0,
                        label=f'rho_bnd = {rho_bnd:.4f}')
        ax5.set_xlim(0.0, 1.0)
        ax5.grid(True, ls=':', lw=0.5)
        ax5.legend(loc='upper right', fontsize=8)

    # ---- (1, 2): LOS weight function f(rhot) ------------------------
    ax6 = axes[1, 2]
    if xs is not None and f_prof is not None:
        ax6.plot(xs, f_prof, 'g-', lw=1.4, label='f (box, geometric)')
        ax6.fill_between(xs, 0.0, f_prof, color='g', alpha=0.15)
        ax6.set_xlabel('rhot')
        ax6.set_ylabel(r'$f(\rho_t)$ = LOS shell fraction')
        ax6.set_title('LOS weight function')
        if rho_bnd is not None and np.isfinite(rho_bnd):
            ax6.axvline(rho_bnd, color='m', ls=':', lw=1.0,
                        label=f'rho_bnd = {rho_bnd:.4f}')
        if losf is not None:
            # Detector-coupling weight, normalised to its own peak so its shape
            # overlays the [0,1] geometric f on the same axis.
            fd = np.asarray(losf['f_det'])
            fdn = fd / np.nanmax(fd) if np.nanmax(fd) > 0 else fd
            ax6.plot(xs, fdn, 'm--', lw=1.4, label='f_det (real LOS, norm.)')
            rc = losf.get('rhot_crit')
            if rc is not None and np.isfinite(rc):
                ax6.axvline(rc, color='c', ls=':', lw=1.0,
                            label=f'rhot_crit = {rc:.3f} (enclosed)')
            if np.isfinite(losf['rho_med']):
                ax6.axvline(losf['rho_med'], color='k', ls='-.', lw=1.0,
                            label=f'rho_50 = {losf["rho_med"]:.4f}')
        ax6.legend(loc='best', fontsize=8)
        ax6.set_xlim(0.0, 1.0)
        ax6.set_ylim(0.0, 1.05)
        ax6.grid(True, ls=':', lw=0.5)

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

    # ------- LOS grid ---------------------------------------------------
    z_bot, z_top = eqs.z_extent()
    R, Z = build_rz_grid(z_bot, z_top, args.Rmin, args.Rmax,
                         args.nR, args.nZ, args.zmargin)
    # Insert the magnetic axis into the sampling arrays so a chord sample sits
    # exactly at (Rmag, Zmag). Without it the centremost f(rhot) bin is
    # under-sampled -- the innermost flux shell is smaller than one (R, Z) cell
    # -- so THKM14 stays ~0 on axis even once psi_n is pinned. Source-agnostic.
    if R[0] <= Rmag <= R[-1]:
        R = np.unique(np.concatenate([R, [Rmag]]))
    if Z[0] <= Zmag <= Z[-1]:
        Z = np.unique(np.concatenate([Z, [Zmag]]))
    print(f'  LOS R in [{R[0]:.3f}, {R[-1]:.3f}] m  ({len(R)} samples)')
    print(f'  LOS Z in [{Z[0]:.3f}, {Z[-1]:.3f}] m  ({len(Z)} samples) '
          f'[box Z extent: {z_bot:.3f} .. {z_top:.3f}]')
    print(f'  LOS centre R_c = {0.5 * (R[0] + R[-1]):.3f} m, '
          f'Rmag - R_c = {Rmag - 0.5 * (R[0] + R[-1]):+.3f} m')

    # ------- (R, Z) -> rhot --------------------------------------------
    rhot, inside, psin_dense = eqs.rhot_on_grid(R, Z)
    Rb, Zb = eqs.lcfs()
    n_inside = int(np.count_nonzero(inside))
    print(f'  Grid points inside LCFS: {n_inside} / {inside.size} '
          f'({100.0 * n_inside / inside.size:.1f} %)')

    # ------- THNTX(R, Z) and integral ----------------------------------
    thntx_grid_native = thntx_on_grid(x_rhot, thntx_x, rhot, inside)
    thntx_grid_si = thntx_grid_native * thntx_to_si  # n/m3/s

    rate, inner_R = integrate_rate(thntx_grid_si, R, Z, args.wtor)
    rate_tor = toroidal_rate(thntx_grid_si, R, Z)

    # Peak emissivity for a sanity check.
    peak = float(np.nanmax(thntx_grid_si))
    peak_idx = np.unravel_index(int(np.nanargmax(thntx_grid_si)),
                                thntx_grid_si.shape)
    peak_R = R[peak_idx[0]]
    peak_Z = Z[peak_idx[1]]

    # Reference: volumetric average emissivity over the LOS box.
    box_volume = args.wtor * (R[-1] - R[0]) * (Z[-1] - Z[0])
    mean_emis = rate / box_volume if box_volume > 0 else float('nan')

    # ------- Equivalent rho_bnd ----------------------------------------
    # Match Rate_tor (the toroidally-integrated LOS-region rate) against the
    # TRANSP cumulative flux-surface integral cumsum(THNTX*DVOL).
    order = np.argsort(x_rhot)
    xs_sorted = np.asarray(x_rhot)[order]
    thntx_sorted = np.asarray(thntx_x)[order]
    dvol_sorted = np.asarray(dvol_x)[order]
    cum_emis = np.cumsum(thntx_sorted * dvol_sorted)
    rho_bnd, total_plasma = find_rho_bnd(rate_tor, xs_sorted, cum_emis)

    # ------- LOS-weighted profile THKM14(rhot) -------------------------
    dvol_m3 = dvol_sorted * dvol_to_m3
    f_profile, los_vol_bin = los_shell_fraction(
        rhot, inside, R, Z, xs_sorted, dvol_m3,
        subgrid=not args.no_subgrid, rmag=Rmag,
    )
    thntx_si_profile = thntx_sorted * thntx_to_si      # n/m3/s
    thkm14_si_profile = thntx_si_profile * f_profile   # n/m3/s
    # Sanity: bin-summed rate from LOS-weighted profile vs Rate_tor.
    rate_from_profile = float(
        np.sum(thntx_sorted * dvol_sorted * f_profile)
    )

    # ------- Report -----------------------------------------------------
    print()
    print('-------------------- LOS volume integral --------------------')
    print(f'  Toroidal width w_tor   = {args.wtor:.3f} m')
    print(f'  LOS box volume         = {box_volume:.3e} m^3   '
          f'(= w_tor * dR * dZ)')
    print(f'  Peak THNTX in box      = {peak:.3e} n/m3/s '
          f'at (R={peak_R:.3f} m, Z={peak_Z:.3f} m)')
    print(f'  Mean THNTX in box      = {mean_emis:.3e} n/m3/s')
    print()
    print('====================  Rate_LOS  ===========================')
    print(f'  Rate_LOS (chord, w_tor={args.wtor:.2f} m) = {rate:.4e} n/s')
    print(f'  Rate_tor (full torus, same R,Z)   = {rate_tor:.4e} n/s')
    print(f'  Rate_tor / Rate_LOS               = {rate_tor / rate:.2f}  '
          f'(expect ~2*pi*R_c/w_tor = {2 * np.pi * 0.5 * (R[0] + R[-1]) / args.wtor:.2f})')
    print('============================================================')
    print()
    print('-------------- Equivalent rho_bnd (TRANSP cumvol) --------------')
    print(f'  Total plasma cumsum(THNTX*DVOL) = {total_plasma:.4e} n/s')
    print(f'  Rate_tor / total_plasma         = {rate_tor / total_plasma:.4f}'
          if total_plasma > 0 else
          '  total_plasma is zero, cannot compute ratio')
    if np.isnan(rho_bnd):
        over = rate_tor / total_plasma if total_plasma > 0 else float('nan')
        print(f'  Rate_tor exceeds the full-plasma integral '
              f'(ratio = {over:.4f});')
        print(f'  no rho_bnd inside the LCFS reproduces it.')
    else:
        print(f'  rho_bnd such that cumsum(THNTX*DVOL)|_rho_bnd = Rate_tor: '
              f'{rho_bnd:.4f}')

        # ------- A / B / C decomposition --------------------------------
        # A = inside LOS  AND  inside rho_bnd       (LOS sees the core)
        # B = inside LOS  AND  outside rho_bnd      (LOS sees the wings)
        # C = inside rho_bnd  AND  outside LOS R-band (core LOS misses)
        # By construction rate(B) = rate(C) since Rate_tor = cumsum|_rho_bnd.
        mask_inner = inside & (rhot <= rho_bnd)
        mask_outer = inside & (rhot > rho_bnd)
        rate_A = toroidal_rate(thntx_grid_si, R, Z, mask=mask_inner)
        rate_B = toroidal_rate(thntx_grid_si, R, Z, mask=mask_outer)
        rate_C = rate_tor - rate_A  # = rate_B analytically; numerical check
        print()
        print('---- A / B / C decomposition (toroidally integrated) ----')
        print(f'  A: inside LOS  AND  inside rho_bnd  = {rate_A:.4e} n/s'
              f'  ({100.0 * rate_A / rate_tor:.2f} %)')
        print(f'  B: inside LOS  AND  outside rho_bnd = {rate_B:.4e} n/s'
              f'  ({100.0 * rate_B / rate_tor:.2f} %)')
        print(f'  C: outside LOS AND  inside rho_bnd  = {rate_C:.4e} n/s'
              f'  (= rate_tor - A)')
        print(f'  Sanity: A + B           = {rate_A + rate_B:.4e} n/s  '
              f'(should equal Rate_tor = {rate_tor:.4e})')
        print(f'  Sanity: |B - C| / B     = '
              f'{abs(rate_B - rate_C) / rate_B if rate_B > 0 else float("nan"):.2e}'
              f'  (analytic zero; trapz/PCHIP residual)')

    print()
    print('--------------- THKM14(rhot) (LOS-weighted) ----------------')
    print(f'  binning mode = '
          f'{"subgrid (anti-aliased)" if not args.no_subgrid else "point (raw histogram)"}')
    print(f'  f(rhot) = LOS shell volume / TRANSP DVOL  '
          f'(max f = {f_profile.max():.4f})')
    print(f'  sum(THNTX * DVOL * f) = {rate_from_profile:.4e} n/s  '
          f'(Rate_tor = {rate_tor:.4e})')
    print(f'  relative residual     = '
          f'{abs(rate_from_profile - rate_tor) / rate_tor:.2e}')
    # Volume-conservation diagnostic: total binned LOS volume vs cumulative
    # LOS box volume should agree (the subgrid splat preserves total mass).
    Rg_chk, _ = np.meshgrid(R, Z, indexing='ij')
    dR_chk = np.gradient(R)
    dZ_chk = np.gradient(Z)
    total_cell_vol = float(np.where(
        inside,
        2.0 * np.pi * Rg_chk * dR_chk[:, None] * dZ_chk[None, :], 0.0,
    ).sum())
    print(f'  total LOS volume (cell sum) = {total_cell_vol:.4e} m^3')
    print(f'  total LOS volume (binned)   = {los_vol_bin.sum():.4e} m^3   '
          f'(consistency check)')

    # ------- Real KM14 LOS cell file (detector-weighted) ----------------
    losf = None
    if args.los_file:
        print()
        print('============  Real KM14 LOS cell file (_KM3.los)  ============')
        cells = read_los_file(args.los_file)
        print(f'  Loaded {cells["R"].size} cells from {args.los_file}')
        print(f'  Cell R in [{cells["R"].min():.3f}, {cells["R"].max():.3f}] m, '
              f'Z in [{cells["Z"].min():.3f}, {cells["Z"].max():.3f}] m, '
              f'|x| <= {np.abs(cells["x"]).max():.3f} m (tangential)')
        losf = los_file_detector_rate(
            cells, eqs, xs_sorted, thntx_sorted, thntx_to_si, dvol_m3,
            subgrid=not args.no_subgrid, median=not args.no_median,
        )
        print(f'  Cells inside LCFS    : {losf["n_inside"]} / '
              f'{cells["R"].size} '
              f'({100.0 * losf["n_inside"] / cells["R"].size:.1f} %)')
        print(f'  LOS volume (sum V)   = {losf["vol_tot"]:.4e} m^3  '
              f'(real narrow chord, cf. box {box_volume:.4e})')
        print(f'  LOS etendue (sum C)  = {losf["c_tot"]:.4e} m^3  '
              f'(= sum V*Omega/4pi)')
        print()
        print(f'  Rate_chord = sum eps*V = {losf["rate_chord"]:.4e} n/s  '
              f'(real chord, no solid angle; cf. box Rate_LOS = {rate:.4e})')
        print(f'  Rate_det   = sum eps*C = {losf["rate_det"]:.4e} n/s  '
              f'<== thermal rate reaching the KM14 detector')
        print(f'  closure sum(THNTX*DVOL*f_det) = '
              f'{losf["rate_from_profile"]:.4e} n/s  '
              f'(rel. residual {abs(losf["rate_from_profile"] - losf["rate_det"]) / losf["rate_det"]:.2e})')
        print()
        print('  --- where does KM14 see its thermal signal? ---')
        print(f'  signal-median radius rho_50 (50% of Rate_det inside) = '
              f'{losf["rho_med"]:.4f}')
        if np.isfinite(rho_bnd):
            print(f'  (cf. idealised-box equivalent rho_bnd = {rho_bnd:.4f})')
        print('==============================================================')

    if args.save:
        import os
        # Repo-local tmp/ (src/tmp), independent of where the CDFs live, so the
        # output is easy to find regardless of --data-dir / the data tree.
        src_dir = os.path.abspath(
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
        outdir = os.path.join(src_dir, 'tmp')
        os.makedirs(outdir, exist_ok=True)
        # Tag the eq source so cdf and ppf runs don't overwrite each other.
        out_path = os.path.join(
            outdir,
            f'{run_id}_KM14_LOS_profile_{args.channel}_{args.eq_source}'
            f'_t{time_jet:.3f}s.txt',
        )
        cols = [xs_sorted, thntx_si_profile, thkm14_si_profile, f_profile,
                dvol_m3]
        col_hdr = ('    rhot           THNTX [n/m3/s]   '
                   'THKM14 [n/m3/s]  f               DVOL [m3]')
        losf_hdr = ''
        if losf is not None:
            cols += [losf['f_det'], losf['thkm14_det_si']]
            col_hdr += '       f_det           THKM14_det [n/m3/s]'
            losf_hdr = (f'real LOS (_KM3.los): Rate_det = {losf["rate_det"]:.4e} '
                        f'n/s   Rate_chord = {losf["rate_chord"]:.4e} n/s   '
                        f'rho_50 = {losf["rho_med"]:.4f}\n')
        data = np.column_stack(cols)
        np.savetxt(
            out_path, data,
            header=(f'pulse {args.pulse}  TRANSP {args.runid}  '
                    f'channel {args.channel} ({th_var})  '
                    f't_JET = {time_jet:.3f} s\n'
                    f'equilibrium: {eqs.label}  t_EQ = {eqs.t_eq_jet:.3f} s\n'
                    f'LOS R in [{R[0]:.3f}, {R[-1]:.3f}] m  '
                    f'(w_tor = {args.wtor:.3f} m)\n'
                    f'Rate_LOS = {rate:.4e} n/s   '
                    f'Rate_tor = {rate_tor:.4e} n/s   '
                    f'rho_bnd = {rho_bnd}\n'
                    f'{losf_hdr}'
                    f'{col_hdr}'),
            fmt='%.8e',
        )
        print(f'Saved LOS profile to {out_path}')

    if args.plot:
        diagnostic_plots(R, Z, rhot, inside, psin_dense, thntx_grid_si,
                         inner_R, Rb, Zb, eqs.native_rhot(),
                         args.pulse, args.runid, time_jet, eqs.t_eq_jet,
                         eqs.label, chan_label, Rmag, Zmag,
                         rho_bnd=rho_bnd,
                         xs=xs_sorted,
                         thntx_si_prof=thntx_si_profile,
                         thkm14_si_prof=thkm14_si_profile,
                         f_prof=f_profile, losf=losf)

    return 0


if __name__ == '__main__':
    sys.exit(main())
