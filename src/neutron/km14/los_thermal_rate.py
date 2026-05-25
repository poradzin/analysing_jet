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
    Z covering the plasma vertical extent (default: LCFS Z range +/- margin)

A fixed effective toroidal width ``W_TOR`` (default 0.4 m) closes the
chord cross section so the volume element is

    dV = W_TOR * dR * dZ

(this is the same convention used by ``estimate_outer_th_fraction.py``).

Algorithm
~~~~~~~~~
1.  Take the Eq class's pre-computed ``eq._psirzmg`` and
    ``eq._psi_norm[tind]`` (already normalized PSIN on the native EFIT
    grid). They share the same C-order flat layout, so ravel'd
    coordinates and values pair correctly with no reshape ambiguity.
2.  Append the magnetic-axis point ``(RMAG, ZMAG, PSIN=0)`` to anchor
    the centre.
3.  ``griddata`` cubic from this scattered set onto a dense (R, Z)
    meshgrid (default 100 x 100) covering the LOS extent; linear
    fallback for edge NaNs.
4.  Mask ``PSIN > 1`` as outside the LCFS; clip the rest to [0, 1] and
    map to ``rhot = sqrt(normalized toroidal flux)`` via
    ``psin_to_sqrt_ftor_norm``.
5.  PCHIP-interpolate the TRANSP ``THNTX(X)`` profile (X = rhot) onto the
    dense rhot map; set ``THNTX = 0`` outside the LCFS.
6.  Convert ``THNTX`` from TRANSP CGS units (N/CM3/SEC) to SI (n/m3/s)
    and integrate using the trapezoidal rule:

        Rate_LOS [n/s] = W_TOR * integral over R, Z of THNTX(R, Z) dR dZ

Usage
-----
    python los_thermal_rate.py 104614 M29 \
            --dda eftp --uid gszepesi --seq 405 -t 53.5268 --plot

The convenience ``key=value`` form (e.g. ``dda='eftp'``) is also accepted.
"""

import argparse
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator, griddata

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, "../.."))
sys.path.append(SRC_DIR)

import profiles as ps
from change_rho import psin_to_sqrt_ftor_norm


# -----------------------------------------------------------------------------
# Geometric defaults for the KM14 LOS
# -----------------------------------------------------------------------------
R_MIN_DEFAULT = 2.70   # inner R limit of the chord [m]
R_MAX_DEFAULT = 3.10   # outer R limit of the chord [m]
W_TOR_DEFAULT = 0.40   # effective toroidal width of the chord [m]
NR_DEFAULT = 100
NZ_DEFAULT = 100
Z_MARGIN_DEFAULT = 0.0  # extra Z padding beyond LCFS [m]


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
    parser.add_argument('runid', type=str, help='TRANSP runid (e.g. M29)')
    parser.add_argument('--dda', type=str, default='eftp',
                        help='Equilibrium PPF DDA (default: eftp)')
    parser.add_argument('--uid', type=str, default='jetppf',
                        help='Equilibrium PPF UID (default: jetppf)')
    parser.add_argument('--seq', type=int, default=0,
                        help='Equilibrium PPF sequence (default: 0)')
    parser.add_argument('-t', '--time', type=float, default=None,
                        help='Time in JET convention (s, t>40). '
                             'If omitted, the midpoint of the TRANSP run is used.')
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
    parser.add_argument('--plot', action='store_true',
                        help='Show diagnostic plots.')
    return parser.parse_args(argv)


# -----------------------------------------------------------------------------
# TRANSP / equilibrium loading
# -----------------------------------------------------------------------------
def get_transp_slice(tr, time_jet):
    """Return (THNTX, rhot, trind) at the JET time *time_jet*."""
    t_jet_grid = tr.t + 40.0
    trind = int(np.abs(t_jet_grid - time_jet).argmin())

    thntx_all = tr.profile('THNTX')
    x_all = tr.x

    thntx = thntx_all[trind] if thntx_all.ndim == 2 else thntx_all
    x_rhot = x_all[trind] if x_all.ndim == 2 else x_all
    return thntx, x_rhot, trind


def build_rz_grid(eq, tind_eq, Rmin, Rmax, nR, nZ, z_margin):
    """Build a (nR x nZ) Cartesian grid covering the LOS region.

    The Z range is taken from the EFIT computational box ``eq._psiz`` (the
    whole vessel), not from ZBND. We saw on pulse 104614 that ``_Zbnd``
    appears to be transposed relative to its reshape parameter and slicing
    with ``[:, tind_eq]`` returns a single boundary point sampled over the
    full discharge -- giving a spuriously narrow Z range. Using the PSI
    grid is robust: vacuum points get zeroed out by the rhot >= 1 mask.
    """
    PSIZ = np.asarray(eq._psiz)
    z_bot = float(np.nanmin(PSIZ))
    z_top = float(np.nanmax(PSIZ))

    R = np.linspace(Rmin, Rmax, nR)
    Z = np.linspace(z_bot - z_margin, z_top + z_margin, nZ)
    return R, Z, z_bot, z_top


def map_grid_to_rhot(R, Z, eq, tind_eq, time_jet):
    """Return rhot(R, Z), inside-LCFS mask, and the raw psin map.

    Uses the Eq class's pre-computed structures directly so there is no
    chance of mismatching coordinate ordering with values:
        * ``eq._psirzmg`` -- meshgrid of (psir, psiz) with default
          ``indexing='xy'`` so axis 0 = Z, axis 1 = R.
        * ``eq._psi_norm[tind]`` -- flat 1D array stored in the *same*
          C-order, so ravel'd coordinates and values pair correctly.

    Algorithm:
        1. Build scattered (R, Z, psin) from psirzmg + psi_norm; append
           the magnetic axis (Rmag, Zmag, psin = 0) to anchor the centre.
        2. griddata cubic onto the dense (R, Z) meshgrid; linear fallback
           for edge NaNs.
        3. psin > 1 -> outside LCFS; clip and map to rhot via
           psin_to_sqrt_ftor_norm.
    """
    Rg_n = np.asarray(eq._psirzmg[0])
    Zg_n = np.asarray(eq._psirzmg[1])
    psin_flat = np.asarray(eq._psi_norm[tind_eq]).ravel()

    pts_R = np.concatenate([Rg_n.ravel(), [float(eq._Rmag[tind_eq])]])
    pts_Z = np.concatenate([Zg_n.ravel(), [float(eq._Zmag[tind_eq])]])
    pts_psin = np.concatenate([psin_flat, [0.0]])

    Rg, Zg = np.meshgrid(R, Z, indexing='ij')  # (nR, nZ)
    query = np.column_stack([Rg.ravel(), Zg.ravel()])

    psin_dense = griddata((pts_R, pts_Z), pts_psin, query, method='cubic')
    if np.any(np.isnan(psin_dense)):
        psin_lin = griddata((pts_R, pts_Z), pts_psin, query, method='linear')
        psin_dense = np.where(np.isnan(psin_dense), psin_lin, psin_dense)
    psin_dense = psin_dense.reshape(Rg.shape)

    inside = (psin_dense >= 0.0) & (psin_dense <= 1.0) & ~np.isnan(psin_dense)
    psin_clipped = np.clip(np.where(np.isnan(psin_dense), 1.0, psin_dense),
                           0.0, 1.0)

    rhot = psin_to_sqrt_ftor_norm(psin_clipped.ravel(), eq, time_jet)
    rhot = np.asarray(rhot).reshape(Rg.shape)
    return rhot, inside, psin_dense, Rg, Zg


def thntx_on_grid(x_rhot, thntx_x, rhot_grid, inside_mask):
    """PCHIP-interpolate THNTX(X) onto rhot(R, Z); zero outside the LCFS."""
    order = np.argsort(x_rhot)
    xs = np.asarray(x_rhot)[order]
    ys = np.asarray(thntx_x)[order]
    uniq = np.concatenate(([True], np.diff(xs) > 0))
    xs = xs[uniq]
    ys = ys[uniq]

    # Pad so the interpolator covers [0, 1]:
    if xs[0] > 0.0:
        xs = np.concatenate(([0.0], xs))
        ys = np.concatenate(([ys[0]], ys))
    if xs[-1] < 1.0:
        xs = np.concatenate((xs, [1.0]))
        ys = np.concatenate((ys, [0.0]))

    f = PchipInterpolator(xs, ys, extrapolate=False)
    out = f(rhot_grid)
    out = np.where(np.isnan(out), 0.0, out)
    out = np.where(inside_mask, out, 0.0)
    return out


def integrate_rate(thntx_grid_si, R, Z, w_tor):
    """Trapezoidal integral of THNTX (in n/m3/s) * w_tor over (R, Z).

    Returns the rate in n/s and the dR-integrated linear emissivity
    (n/s per meter of Z) for diagnostics.
    """
    inner_R = np.trapz(thntx_grid_si, R, axis=0)        # over R -> shape (nZ,)
    rate = w_tor * np.trapz(inner_R, Z)                 # over Z -> scalar
    return float(rate), inner_R


# -----------------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------------
def _lcfs_contour(eq, tind_eq):
    """Return (Rbnd, Zbnd) at tind_eq using the correct time-first indexing.

    profiles.py reshapes RBND/ZBND with ``(-1, n_t)`` but the array actually
    ends up shaped ``(n_t, n_bnd)`` (confirmed empirically: on pulse 104614
    only ``_Zbnd[tind, :]`` gives the expected ~[-1.37, 1.68] range, while
    ``_Zbnd[:, tind]`` returns a single bnd point sampled across the
    discharge). The LCFS may be padded with zeros at unused slots, so we
    drop the trailing (0, 0) entries.
    """
    Rb = np.asarray(eq._Rbnd[tind_eq, :], dtype=float)
    Zb = np.asarray(eq._Zbnd[tind_eq, :], dtype=float)
    valid = ~((Rb == 0.0) & (Zb == 0.0)) & ~np.isnan(Rb) & ~np.isnan(Zb)
    return Rb[valid], Zb[valid]


def diagnostic_plots(R, Z, rhot, inside, psin_dense, thntx_grid_si, inner_R,
                     eq, tind_eq, pulse, runid, time_jet,
                     Rmag, Zmag):
    """Four-panel diagnostic:
        (0,0) rhot(R,Z) in LOS region + LCFS overlay
        (0,1) THNTX(R,Z) in LOS region + LCFS overlay
        (1,0) Full poloidal LCFS shape with LOS box marked
        (1,1) R-integrated emissivity vs Z
    """
    Rg, Zg = np.meshgrid(R, Z, indexing='ij')
    rhot_plot = np.where(inside, rhot, np.nan)

    Rb, Zb = _lcfs_contour(eq, tind_eq)

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 11.0))
    fig.suptitle(
        f'KM14 LOS thermal-neutron rate  -  pulse {pulse}  TRANSP {runid}  '
        f't_JET={time_jet:.3f}s  t_EQ={eq.t[tind_eq]:.3f}s'
    )

    ax1 = axes[0, 0]
    pc1 = ax1.pcolormesh(Rg, Zg, rhot_plot, shading='auto', cmap='viridis')
    fig.colorbar(pc1, ax=ax1, label='rhot')
    ax1.plot(Rb, Zb, 'w-', lw=1.2, label='LCFS (RBND,ZBND)')
    ax1.contour(Rg, Zg, psin_dense, levels=[1.0],
                colors='white', linewidths=0.8, linestyles='--')
    ax1.plot([Rmag], [Zmag], 'r+', ms=10, label='axis')
    ax1.set_xlim(R[0], R[-1])
    ax1.set_ylim(Z[0], Z[-1])
    ax1.set_xlabel('R [m]')
    ax1.set_ylabel('Z [m]')
    ax1.set_title('rhot on LOS grid')
    ax1.legend(loc='upper right', fontsize=8)
    ax1.set_aspect('equal', adjustable='box')

    ax2 = axes[0, 1]
    pc2 = ax2.pcolormesh(Rg, Zg, thntx_grid_si, shading='auto', cmap='inferno')
    fig.colorbar(pc2, ax=ax2, label='THNTX [n/m3/s]')
    ax2.plot(Rb, Zb, 'c-', lw=1.2)
    ax2.plot([Rmag], [Zmag], 'w+', ms=10)
    ax2.set_xlim(R[0], R[-1])
    ax2.set_ylim(Z[0], Z[-1])
    ax2.set_xlabel('R [m]')
    ax2.set_ylabel('Z [m]')
    ax2.set_title('THNTX(R, Z) inside LOS')
    ax2.set_aspect('equal', adjustable='box')

    ax3 = axes[1, 0]
    ax3.plot(Rb, Zb, 'b-', lw=1.4, label='LCFS')
    ax3.plot([Rmag], [Zmag], 'r+', ms=10, label=f'axis ({Rmag:.3f}, {Zmag:.3f})')
    los_box_R = [R[0], R[-1], R[-1], R[0], R[0]]
    los_box_Z = [Z[0], Z[0], Z[-1], Z[-1], Z[0]]
    ax3.plot(los_box_R, los_box_Z, 'r-', lw=1.0,
             label=f'LOS box R=[{R[0]:.2f},{R[-1]:.2f}]')
    ax3.axvline(0.5 * (R[0] + R[-1]), color='r', ls=':', lw=0.8,
                label=f'LOS centre R={0.5*(R[0]+R[-1]):.2f} m')
    ax3.set_xlabel('R [m]')
    ax3.set_ylabel('Z [m]')
    ax3.set_title('Poloidal cross-section: LCFS + KM14 LOS')
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

    plt.tight_layout()
    plt.show()


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------
def main(argv=None):
    args = parse_args(argv)

    # ------- TRANSP -----------------------------------------------------
    print(f'Loading TRANSP pulse {args.pulse} run {args.runid} ...')
    tr = ps.Transp(args.pulse, args.runid)

    t_jet_grid = tr.t + 40.0
    if args.time is None:
        time_jet = float(t_jet_grid[len(t_jet_grid) // 2])
        print(f'No time supplied; using midpoint t_JET = {time_jet:.3f} s.')
    else:
        time_jet = float(args.time)

    thntx_x, x_rhot, trind = get_transp_slice(tr, time_jet)

    tr.add_data('THNTX')
    unit_th = tr.units('THNTX') or ''
    # TRANSP THNTX is in N/CM3/SEC; convert to N/M3/SEC for the integral.
    thntx_to_si = 1.0e6 if 'CM' in unit_th.upper() else 1.0

    print(f'TRANSP slice [{trind}]: t_JET = {t_jet_grid[trind]:.3f} s '
          f'(t_TRANSP = {tr.t[trind]:.3f} s)')
    print(f'  THNTX units = [{unit_th}],  conversion to n/m3/s = {thntx_to_si:.1e}')
    print(f'  THNTX(rhot=0) = {thntx_x[0]:.3e} [{unit_th}]  '
          f'= {thntx_x[0] * thntx_to_si:.3e} n/m3/s')

    # ------- Equilibrium ------------------------------------------------
    print(f'Loading equilibrium PPF {args.dda}/{args.uid}/{args.seq} ...')
    eq = ps.Eq(args.pulse, dda=args.dda, uid=args.uid, seq=args.seq)
    tind_eq = int(np.abs(eq.t - time_jet).argmin())
    Rmag = float(eq._Rmag[tind_eq])
    Zmag = float(eq._Zmag[tind_eq])
    print(f'Equilibrium slice [{tind_eq}]: t = {eq.t[tind_eq]:.3f} s')
    print(f'  Rmag = {Rmag:.3f} m, Zmag = {Zmag:.3f} m')

    # ------- LOS grid ---------------------------------------------------
    R, Z, z_bot, z_top = build_rz_grid(
        eq, tind_eq, args.Rmin, args.Rmax, args.nR, args.nZ, args.zmargin
    )
    print(f'  LOS R in [{R[0]:.3f}, {R[-1]:.3f}] m  ({args.nR} samples)')
    print(f'  LOS Z in [{Z[0]:.3f}, {Z[-1]:.3f}] m  ({args.nZ} samples) '
          f'[EFIT PSIZ extent: {z_bot:.3f} .. {z_top:.3f}]')
    print(f'  LOS centre R_c = {0.5 * (R[0] + R[-1]):.3f} m, '
          f'Rmag - R_c = {Rmag - 0.5 * (R[0] + R[-1]):+.3f} m')

    # ------- (R, Z) -> rhot --------------------------------------------
    rhot, inside, psin_dense, Rg, Zg = map_grid_to_rhot(
        R, Z, eq, tind_eq, time_jet
    )
    n_inside = int(np.count_nonzero(inside))
    print(f'  Grid points inside LCFS: {n_inside} / {inside.size} '
          f'({100.0 * n_inside / inside.size:.1f} %)')

    # ------- THNTX(R, Z) and integral ----------------------------------
    thntx_grid_native = thntx_on_grid(x_rhot, thntx_x, rhot, inside)
    thntx_grid_si = thntx_grid_native * thntx_to_si  # n/m3/s

    rate, inner_R = integrate_rate(thntx_grid_si, R, Z, args.wtor)

    # Peak emissivity for a sanity check.
    peak = float(np.nanmax(thntx_grid_si))
    peak_idx = np.unravel_index(int(np.nanargmax(thntx_grid_si)),
                                thntx_grid_si.shape)
    peak_R = R[peak_idx[0]]
    peak_Z = Z[peak_idx[1]]

    # Reference: volumetric average emissivity over the LOS box.
    box_volume = args.wtor * (R[-1] - R[0]) * (Z[-1] - Z[0])
    mean_emis = rate / box_volume if box_volume > 0 else float('nan')

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
    print(f'  Rate_LOS = {rate:.4e} n/s')
    print('============================================================')

    if args.plot:
        diagnostic_plots(R, Z, rhot, inside, psin_dense, thntx_grid_si,
                         inner_R, eq, tind_eq, args.pulse, args.runid,
                         time_jet, Rmag, Zmag)

    return 0


if __name__ == '__main__':
    sys.exit(main())
