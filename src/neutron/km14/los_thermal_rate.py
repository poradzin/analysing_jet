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
4.  Build the inside-LCFS mask as a point-in-polygon test against
    ``(RBND, ZBND)`` (closed boundary). This excludes the private flux
    region below the X-point, which is where ``psin < 1`` but TRANSP
    does not model emission. Then clip ``psin`` to [0, 1] and map to
    ``rhot = sqrt(normalized toroidal flux)`` via
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
from matplotlib.path import Path as MplPath
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
    parser.add_argument('--save', action='store_true',
                        help='Save (rhot, THNTX, THKM14, f) profile to tmp/.')
    parser.add_argument('--no-subgrid', action='store_true',
                        help='Disable subgrid (R,Z)-cell rhot distribution '
                             'for f(rhot); use point binning instead. '
                             'Subgrid is on by default to suppress aliasing.')
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

    The ``inside`` mask is the closed LCFS polygon test on (RBND, ZBND),
    NOT ``psin <= 1``. The latter would include the private flux region
    below the X-point -- ``psin`` is < 1 there but TRANSP does not model
    emission outside the confined plasma, so including the PFR would
    bias the integral upward by a small but spurious amount.
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

    # Inside-LCFS mask from the closed boundary polygon -- excludes the PFR.
    Rb, Zb = _lcfs_contour(eq, tind_eq)
    lcfs_poly = MplPath(np.column_stack([Rb, Zb]))
    inside = lcfs_poly.contains_points(query).reshape(Rg.shape)

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


def cumulative_emission(x_rhot, thntx_x, dvol_x):
    """Return (xs_sorted, cumsum(THNTX*DVOL)) on the TRANSP rhot grid [n/s].

    Note: THNTX * DVOL has units of n/s regardless of CGS vs SI
    (N/CM3/SEC * CM**3 = N/SEC), so no conversion is needed.
    """
    order = np.argsort(x_rhot)
    xs = np.asarray(x_rhot)[order]
    th = np.asarray(thntx_x)[order]
    dv = np.asarray(dvol_x)[order]
    return xs, np.cumsum(th * dv)


def find_rho_bnd(target, xs, cum_emis):
    """Find rho_bnd such that cumsum(THNTX*DVOL)|_{rho_bnd} = target.

    Returns (rho_bnd, total_plasma). rho_bnd is NaN if target exceeds the
    full-plasma cumulative emission cum_emis[-1].
    """
    total = float(cum_emis[-1])
    if target <= cum_emis[0]:
        return float(xs[0]), total
    if target >= total:
        return float('nan'), total
    keep = np.concatenate(([True], np.diff(cum_emis) > 0))
    f_inv = PchipInterpolator(cum_emis[keep], xs[keep], extrapolate=False)
    return float(f_inv(target)), total


def los_shell_fraction(rhot_los, inside, R, Z, x_rhot, dvol_x_m3,
                       subgrid=True):
    """LOS-weighted fraction f(rhot) of each TRANSP flux shell.

    For every (R, Z) cell of the LOS grid we accumulate the toroidal
    volume element ``dV = 2*pi*R*dR*dZ`` (zero outside the LCFS) into
    bins centred on the TRANSP rhot grid, then

        f(rhot_i) = LOS_vol_bin(i) / DVOL_TRANSP_bin(i)

    so ``f`` lies in [0, 1] and ``sum(THNTX * DVOL * f) == Rate_tor`` by
    construction (modulo discretisation).

    Binning modes
    -------------
    * ``subgrid=False`` -- point binning: each cell's volume goes to a
      single bin chosen by its centre rhot. Suffers from aliasing when
      the rhot change across one (R, Z) cell exceeds one TRANSP bin
      width (typical for the Z direction at modest nZ).
    * ``subgrid=True`` (default) -- anti-aliased: each cell spans a rhot
      range estimated from ``|grad rhot| * cell_size`` (L1 corner half-
      range). Its volume is distributed uniformly across that range and
      added to every TRANSP bin proportional to overlap, vectorised via
      ``np.add.at``. Eliminates the period-N oscillation in f(rhot) and
      correctly accounts for the mismatch between a rectangular (R, Z)
      grid and curved flux surfaces.

    Parameters
    ----------
    rhot_los  : 2D array, rhot on the LOS (R, Z) grid
    inside    : 2D bool array, inside-LCFS mask (same shape)
    R, Z      : 1D LOS grid axes [m]
    x_rhot    : 1D TRANSP rhot grid (sorted ascending)
    dvol_x_m3 : 1D TRANSP DVOL converted to m^3 (same shape as x_rhot)
    subgrid   : bool, see above

    Returns
    -------
    f            : 1D weight, shape (len(x_rhot),)
    los_vol_bin  : 1D LOS toroidal volume per TRANSP shell [m^3]
    """
    R_arr = np.asarray(R, dtype=float)
    Z_arr = np.asarray(Z, dtype=float)
    dR = np.gradient(R_arr)
    dZ = np.gradient(Z_arr)
    Rg, _ = np.meshgrid(R_arr, Z_arr, indexing='ij')
    cell_vol = np.where(
        inside,
        2.0 * np.pi * Rg * dR[:, None] * dZ[None, :],
        0.0,
    )

    xs = np.asarray(x_rhot, dtype=float)
    edges = np.concatenate(([0.0], 0.5 * (xs[:-1] + xs[1:]), [1.0]))
    n_bins = len(edges) - 1

    if not subgrid:
        los_vol_bin, _ = np.histogram(
            rhot_los.ravel(), bins=edges, weights=cell_vol.ravel()
        )
    else:
        # rhot half-extent across the cell (L1 corner half-range)
        grad_R, grad_Z = np.gradient(rhot_los, R_arr, Z_arr)
        h_R = 0.5 * np.abs(grad_R) * dR[:, None]
        h_Z = 0.5 * np.abs(grad_Z) * dZ[None, :]
        half = h_R + h_Z

        rhot_c = rhot_los.ravel()
        vol_flat = cell_vol.ravel()
        half_flat = half.ravel()

        keep = vol_flat > 0.0
        rhot_c = rhot_c[keep]
        vol = vol_flat[keep]
        half_k = half_flat[keep]

        # Clip to [0, 1] so volume isn't spilled outside the physical range.
        rhot_lo = np.clip(rhot_c - half_k, 0.0, 1.0)
        rhot_hi = np.clip(rhot_c + half_k, 0.0, 1.0)
        span = rhot_hi - rhot_lo

        # Which bins each cell touches.
        i_lo = np.clip(
            np.searchsorted(edges, rhot_lo, side='right') - 1,
            0, n_bins - 1,
        )
        i_hi = np.clip(
            np.searchsorted(edges, rhot_hi, side='right') - 1,
            0, n_bins - 1,
        )
        n_per_cell = i_hi - i_lo + 1

        # Build flat (cell, bin) assignment vectors.
        total = int(n_per_cell.sum())
        flat_cell = np.repeat(np.arange(len(vol)), n_per_cell)
        starts = np.concatenate(([0], np.cumsum(n_per_cell[:-1])))
        within_offset = np.arange(total) - np.repeat(starts, n_per_cell)
        flat_bin = np.repeat(i_lo, n_per_cell) + within_offset

        cell_lo_b = rhot_lo[flat_cell]
        cell_hi_b = rhot_hi[flat_cell]
        cell_span_b = span[flat_cell]
        cell_vol_b = vol[flat_cell]
        bin_lo = edges[flat_bin]
        bin_hi = edges[flat_bin + 1]
        ovl = np.maximum(
            0.0,
            np.minimum(cell_hi_b, bin_hi) - np.maximum(cell_lo_b, bin_lo),
        )
        # Cells with zero rhot span -> all volume to their single bin.
        zero_span = cell_span_b == 0
        weight = np.where(
            zero_span,
            cell_vol_b,
            cell_vol_b * ovl
            / np.where(cell_span_b > 0, cell_span_b, 1.0),
        )

        los_vol_bin = np.zeros(n_bins)
        np.add.at(los_vol_bin, flat_bin, weight)

    dvol = np.asarray(dvol_x_m3, dtype=float)
    with np.errstate(invalid='ignore', divide='ignore'):
        f = np.where(dvol > 0.0, los_vol_bin / dvol, 0.0)
    f = np.clip(f, 0.0, 1.0)
    return f, los_vol_bin


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
                     Rmag, Zmag, rho_bnd=None,
                     xs=None, thntx_si_prof=None,
                     thkm14_si_prof=None, f_prof=None):
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

    Rb, Zb = _lcfs_contour(eq, tind_eq)

    # Vessel-wide rhot map for the poloidal overview panel (used to draw
    # the rho_bnd surface across the whole plasma). Masked outside the LCFS
    # polygon so we don't trace spurious contours in the PFR or vacuum.
    Rg_n = np.asarray(eq._psirzmg[0])
    Zg_n = np.asarray(eq._psirzmg[1])
    psin_native = np.asarray(eq._psi_norm[tind_eq]).reshape(Rg_n.shape)
    psin_native_clip = np.clip(
        np.where(np.isnan(psin_native), 1.0, psin_native), 0.0, 1.0
    )
    rhot_native = psin_to_sqrt_ftor_norm(
        psin_native_clip.ravel(), eq, time_jet
    ).reshape(Rg_n.shape)
    lcfs_poly = MplPath(np.column_stack([Rb, Zb]))
    native_mask = lcfs_poly.contains_points(
        np.column_stack([Rg_n.ravel(), Zg_n.ravel()])
    ).reshape(Rg_n.shape)
    rhot_native = np.where(native_mask, rhot_native, np.nan)

    rho_bnd_lbl = (f'rho_bnd = {rho_bnd:.4f}'
                   if rho_bnd is not None and np.isfinite(rho_bnd)
                   else 'rho_bnd: N/A')

    fig, axes = plt.subplots(2, 3, figsize=(17.0, 10.0))
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
    if rho_bnd is not None and np.isfinite(rho_bnd):
        cs = ax3.contour(Rg_n, Zg_n, rhot_native, levels=[rho_bnd],
                         colors='magenta', linewidths=1.4)
        if cs.collections:
            cs.collections[0].set_label(rho_bnd_lbl)
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

    # ---- (0, 2): emissivity profiles --------------------------------
    ax5 = axes[0, 2]
    if xs is not None and thntx_si_prof is not None:
        ax5.plot(xs, thntx_si_prof, 'b-', lw=1.4,
                 label='THNTX (TRANSP)')
        ax5.plot(xs, thkm14_si_prof, 'r-', lw=1.4,
                 label=r'THKM14 = THNTX$\cdot f$')
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
        ax6.plot(xs, f_prof, 'g-', lw=1.4)
        ax6.fill_between(xs, 0.0, f_prof, color='g', alpha=0.15)
        ax6.set_xlabel('rhot')
        ax6.set_ylabel(r'$f(\rho_t)$ = LOS shell fraction')
        ax6.set_title('LOS weight function')
        if rho_bnd is not None and np.isfinite(rho_bnd):
            ax6.axvline(rho_bnd, color='m', ls=':', lw=1.0,
                        label=f'rho_bnd = {rho_bnd:.4f}')
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
    tr.add_data('DVOL')
    unit_dv = tr.units('DVOL') or ''
    dvol_to_m3 = 1.0e-6 if 'CM' in unit_dv.upper() else 1.0

    dvol_all = tr.profile('DVOL')
    dvol_x = dvol_all[trind] if dvol_all.ndim == 2 else dvol_all

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
        subgrid=not args.no_subgrid,
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

    if args.save:
        outdir = os.path.join(SRC_DIR, 'tmp')
        os.makedirs(outdir, exist_ok=True)
        out_path = os.path.join(
            outdir,
            f'{args.pulse}_{args.runid}_KM14_LOS_profile_t{time_jet:.3f}s.txt',
        )
        data = np.column_stack([
            xs_sorted,
            thntx_si_profile,
            thkm14_si_profile,
            f_profile,
            dvol_m3,
        ])
        np.savetxt(
            out_path, data,
            header=(f'pulse {args.pulse}  TRANSP {args.runid}  '
                    f't_JET = {time_jet:.3f} s\n'
                    f'eq = {args.dda}/{args.uid}/{args.seq}  '
                    f't_EQ = {eq.t[tind_eq]:.3f} s\n'
                    f'LOS R in [{R[0]:.3f}, {R[-1]:.3f}] m  '
                    f'(w_tor = {args.wtor:.3f} m)\n'
                    f'Rate_LOS = {rate:.4e} n/s   '
                    f'Rate_tor = {rate_tor:.4e} n/s   '
                    f'rho_bnd = {rho_bnd}\n'
                    f'    rhot           THNTX [n/m3/s]   '
                    f'THKM14 [n/m3/s]  f               DVOL [m3]'),
            fmt='%.8e',
        )
        print(f'Saved LOS profile to {out_path}')

    if args.plot:
        diagnostic_plots(R, Z, rhot, inside, psin_dense, thntx_grid_si,
                         inner_R, eq, tind_eq, args.pulse, args.runid,
                         time_jet, Rmag, Zmag,
                         rho_bnd=rho_bnd,
                         xs=xs_sorted,
                         thntx_si_prof=thntx_si_profile,
                         thkm14_si_prof=thkm14_si_profile,
                         f_prof=f_profile)

    return 0


if __name__ == '__main__':
    sys.exit(main())
