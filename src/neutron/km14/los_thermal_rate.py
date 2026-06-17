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
import sys

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib.path import Path as MplPath
from scipy.interpolate import PchipInterpolator

import bt_zone_integrator as bzi
from los_th_bt_ratio import CdfEquilibrium


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
    return parser.parse_args(argv)


# -----------------------------------------------------------------------------
# Thermal profiles from the TRANSP main CDF
# -----------------------------------------------------------------------------
def read_time_grid(cdf_path):
    """Return TIME3 (TRANSP seconds) from the main CDF."""
    with Dataset(cdf_path, "r") as d:
        return np.array(d["TIME3"][:])


def read_thermal_slice(cdf_path, th_var, trind):
    """Return (thntx, x_rhot, dvol, thntx_to_si, dvol_to_m3, units) at *trind*."""
    with Dataset(cdf_path, "r") as d:
        if th_var not in d.variables:
            raise KeyError(f"{th_var} not present in {cdf_path}")
        thntx = np.array(d[th_var][trind])
        x_rhot = np.array(d["X"][trind])
        dvol = np.array(d["DVOL"][trind])
        unit_th = getattr(d[th_var], "units", "") or ""
        unit_dv = getattr(d["DVOL"], "units", "") or ""
    thntx_to_si = 1.0e6 if 'CM' in unit_th.upper() else 1.0   # N/CM3/SEC -> n/m3/s
    dvol_to_m3 = 1.0e-6 if 'CM' in unit_dv.upper() else 1.0   # CM**3 -> m^3
    return thntx, x_rhot, dvol, thntx_to_si, dvol_to_m3, unit_th


# -----------------------------------------------------------------------------
# Real KM14 line-of-sight cell file (_KM3.los, supplied by M. Nocente)
# -----------------------------------------------------------------------------
def read_los_file(path):
    """Read the KM14 LOS cell file (``_KM3.los``).

    Each row is one LOS cell with 8 whitespace-separated columns

        x, y, z, C, V, u, v, w

    (see ``KM3_LoS_readme.txt``). All in a fixed, discharge-/time-independent
    right-handed Cartesian frame: the poloidal plane is x-z, the toroidal plane
    is x-y, so

      * x       -- tangential (toroidal) offset from the detector poloidal
                   plane [m] (small, |x| < 0.1 m)
      * y       -- horizontal in-plane coordinate [m] -- essentially the major
                   radius; the cell major radius is R = sqrt(x^2 + y^2)
      * z       -- vertical coordinate Z [m]
      * C       -- etendue weight [m^3]. The detector solid angle seen from the
                   cell is Omega = 4*pi*C/V, so for an *isotropic* emitter of
                   volumetric emissivity eps the rate of neutrons that actually
                   reach the detector is eps * V * Omega/(4*pi) = eps * C.
      * V       -- cell volume [m^3]
      * u, v, w -- unit vector of the emission direction a neutron born in the
                   cell must have to reach the detector. Needed for the
                   anisotropic beam-target channel; unused for the isotropic
                   thermal channel handled here, but read and returned so the
                   sister script can reuse this loader.

    Returns
    -------
    dict with 1-D arrays ``x, y, z, C, V, u, v, w`` plus derived
    ``R = sqrt(x^2 + y^2)`` and ``Z = z``.
    """
    raw = np.loadtxt(path)
    if raw.ndim != 2 or raw.shape[1] != 8:
        raise ValueError(
            f"{path}: expected an (N, 8) table (x,y,z,C,V,u,v,w), "
            f"got shape {raw.shape}"
        )
    x, y, z, C, V, u, v, w = raw.T
    return {
        "x": x, "y": y, "z": z, "C": C, "V": V, "u": u, "v": v, "w": w,
        "R": np.hypot(x, y), "Z": z.copy(),
    }


# -----------------------------------------------------------------------------
# Equilibrium sources -- common interface used by main()
#
#   .tind        thermal/eq time index on the relevant time grid (int)
#   .t_eq_jet    equilibrium time in JET convention [s]
#   .label       human-readable source description
#   .Rmag, .Zmag magnetic axis [m]
#   .z_extent()  -> (z_bot, z_top) computational-box Z range [m]
#   .lcfs()      -> (Rb, Zb) closed LCFS polygon [m]
#   .rhot_on_grid(R, Z) -> (rhot, inside, psin_dense) on the (nR, nZ) grid
#   .native_rhot()      -> (Rg_n, Zg_n, rhot_native) masked to the LCFS, or None
# -----------------------------------------------------------------------------
def _boundary_from_moments(d, tind, ntheta=257):
    """Reconstruct the LCFS (R, Z) [m] from TRANSP asymmetric boundary moments.

    TRANSP stores the plasma boundary as a Fourier series in a poloidal angle:

        R(theta) = RMCB0 + sum_{n>=1} [RMCBn cos(n th) + RMSBn sin(n th)]
        Z(theta) = YMCB0 + sum_{n>=1} [YMCBn cos(n th) + YMSBn sin(n th)]

    (units cm in the CDF). Verified to reproduce the ``_fi`` RSURF/ZSURF
    outermost surface to ~1e-4 cm on 104614 M30.
    """
    th = np.linspace(0.0, 2.0 * np.pi, ntheta)
    R = np.zeros_like(th)
    Z = np.zeros_like(th)
    n = 0
    while f"RMCB{n}" in d.variables:
        R += float(d[f"RMCB{n}"][tind]) * np.cos(n * th)
        Z += float(d[f"YMCB{n}"][tind]) * np.cos(n * th)
        if n >= 1 and f"RMSB{n}" in d.variables:
            R += float(d[f"RMSB{n}"][tind]) * np.sin(n * th)
            Z += float(d[f"YMSB{n}"][tind]) * np.sin(n * th)
        n += 1
    return R / 100.0, Z / 100.0


def _rhot_scatter(Rq, Zq, pts_R, pts_Z, pts_psin, lcfs_RZ, psin_to_rhot):
    """Evaluate (rhot, inside, psin) at scattered query points (Rq, Zq).

    Shared core of ``EqCDF``/``EqPPF`` -- griddata-cubic of the psi_n node set
    (PSIRZ grid + pinned axis at psi_n=0 + pinned LCFS at psi_n=1), linear
    fallback for NaNs, point-in-LCFS-polygon mask, then psi_n -> rhot via the
    source-specific ``psin_to_rhot`` callable. Works for any-shape flat inputs,
    so it serves both the dense (R, Z) box grid (via ``rhot_on_grid``) and the
    irregular cloud of real LOS cells (via ``rhot_at_points``).
    """
    from scipy.interpolate import griddata
    Rq = np.asarray(Rq, dtype=float).ravel()
    Zq = np.asarray(Zq, dtype=float).ravel()
    query = np.column_stack([Rq, Zq])
    psin = griddata((pts_R, pts_Z), pts_psin, query, method='cubic')
    if np.any(np.isnan(psin)):
        psin_lin = griddata((pts_R, pts_Z), pts_psin, query, method='linear')
        psin = np.where(np.isnan(psin), psin_lin, psin)
    inside = MplPath(lcfs_RZ).contains_points(query)
    psin_clip = np.clip(np.where(np.isnan(psin), 1.0, psin), 0.0, 1.0)
    rhot = np.clip(np.asarray(psin_to_rhot(psin_clip)), 0.0, 1.0)
    return rhot, inside, psin


class EqCDF:
    """Self-contained equilibrium from the TRANSP main CDF (no ppf)."""

    label = "TRANSP CDF (PSIRZ + boundary moments)"

    def __init__(self, cdf_path, time_jet):
        self.ceq = CdfEquilibrium(cdf_path, time_jet - 40.0)  # TRANSP time
        self.tind = self.ceq.tind
        self.t_eq_jet = self.ceq.t_used + 40.0
        self.RG = self.ceq.RG
        self.ZG = self.ceq.ZG
        with Dataset(cdf_path, "r") as d:
            self.Rmag = float(d["RAXIS"][self.tind]) / 100.0
            self.Zmag = float(d["YAXIS"][self.tind]) / 100.0
            self._Rb, self._Zb = _boundary_from_moments(d, self.tind)

    def z_extent(self):
        return float(self.ZG.min()), float(self.ZG.max())

    def lcfs(self):
        return self._Rb, self._Zb

    def _scatter_nodes(self):
        """psi_n node set for griddata: PSIRZ grid + pinned axis + pinned LCFS.

        The magnetic axis is pinned at psi_n = 0 because plain interpolation of
        the coarse PSIRZ grid (dR ~ 2 cm) cannot recover the true paraboloidal
        axis minimum (it floors psi_n ~ 3e-4, i.e. rhot ~ 0.02, zeroing the
        innermost f(rhot) bins). The LCFS polygon is pinned at psi_n = 1 because
        the cubic interpolation otherwise floors short of 1 inside the mask
        (max psi_n ~ 0.999), starving the outermost flux shell of LOS cells and
        cratering f(rhot) in the last 1-2 bins.
        """
        Rg_n, Zg_n = np.meshgrid(self.RG, self.ZG)  # (nZ, nR), matches psin[iz,ir]
        Rb, Zb = self.lcfs()
        pts_R = np.concatenate([Rg_n.ravel(), [self.Rmag], Rb])
        pts_Z = np.concatenate([Zg_n.ravel(), [self.Zmag], Zb])
        pts_psin = np.concatenate([self.ceq.psin.ravel(), [0.0],
                                   np.ones(Rb.size)])
        return pts_R, pts_Z, pts_psin

    def rhot_at_points(self, Rq, Zq):
        """(rhot, inside, psin) at scattered query points (flat arrays)."""
        pts_R, pts_Z, pts_psin = self._scatter_nodes()
        Rb, Zb = self.lcfs()
        return _rhot_scatter(
            Rq, Zq, pts_R, pts_Z, pts_psin, np.column_stack([Rb, Zb]),
            lambda p: np.interp(p, self.ceq._psin_b, self.ceq._rhot_b),
        )

    def rhot_on_grid(self, R, Z):
        Rg, Zg = np.meshgrid(R, Z, indexing='ij')
        rhot, inside, psin = self.rhot_at_points(Rg.ravel(), Zg.ravel())
        return (rhot.reshape(Rg.shape), inside.reshape(Rg.shape),
                psin.reshape(Rg.shape))

    def native_rhot(self):
        Rg_n, Zg_n = np.meshgrid(self.RG, self.ZG)  # (nZ, nR), matches psin[iz,ir]
        psin_clip = np.clip(self.ceq.psin, 0.0, 1.0)
        rhot = np.clip(np.interp(psin_clip.ravel(), self.ceq._psin_b,
                                 self.ceq._rhot_b), 0.0, 1.0).reshape(Rg_n.shape)
        mask = MplPath(np.column_stack(self.lcfs())).contains_points(
            np.column_stack([Rg_n.ravel(), Zg_n.ravel()])).reshape(Rg_n.shape)
        return Rg_n, Zg_n, np.where(mask, rhot, np.nan)


class EqPPF:
    """JET PPF equilibrium via profiles.Eq (imports ppf lazily)."""

    label = "PPF equilibrium"

    def __init__(self, pulse, dda, uid, seq, time_jet):
        # profiles.py / change_rho.py live in src/ (two levels up from this
        # km14/ dir), which is not on sys.path when the script is run directly
        # -- only the script's own dir is. Add src/ before the lazy import.
        import os
        src_dir = os.path.abspath(
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
        if src_dir not in sys.path:
            sys.path.append(src_dir)
        import profiles as ps
        from change_rho import psin_to_sqrt_ftor_norm
        self._psin_to_rhot = psin_to_sqrt_ftor_norm
        self.eq = ps.Eq(pulse, dda=dda, uid=uid, seq=seq)
        self.time_jet = time_jet
        self.tind = int(np.abs(self.eq.t - time_jet).argmin())
        self.t_eq_jet = float(self.eq.t[self.tind])
        self.Rmag = float(self.eq._Rmag[self.tind])
        self.Zmag = float(self.eq._Zmag[self.tind])
        self.label = f"PPF equilibrium {dda}/{uid}/{seq}"

    def z_extent(self):
        PSIZ = np.asarray(self.eq._psiz)
        return float(np.nanmin(PSIZ)), float(np.nanmax(PSIZ))

    def lcfs(self):
        return _lcfs_contour(self.eq, self.tind)

    def _scatter_nodes(self):
        """psi_n node set: PSIRZ grid + pinned axis (psi_n=0) + pinned LCFS (1).

        Same pinning rationale as ``EqCDF._scatter_nodes`` (axis floor and edge
        crater), so both equilibrium sources share one scattered evaluator.
        """
        eq, tind_eq = self.eq, self.tind
        Rg_n = np.asarray(eq._psirzmg[0])
        Zg_n = np.asarray(eq._psirzmg[1])
        psin_flat = np.asarray(eq._psi_norm[tind_eq]).ravel()
        Rb, Zb = self.lcfs()
        pts_R = np.concatenate([Rg_n.ravel(), [self.Rmag], Rb])
        pts_Z = np.concatenate([Zg_n.ravel(), [self.Zmag], Zb])
        pts_psin = np.concatenate([psin_flat, [0.0], np.ones(Rb.size)])
        return pts_R, pts_Z, pts_psin

    def rhot_at_points(self, Rq, Zq):
        """(rhot, inside, psin) at scattered query points (flat arrays)."""
        pts_R, pts_Z, pts_psin = self._scatter_nodes()
        Rb, Zb = self.lcfs()
        return _rhot_scatter(
            Rq, Zq, pts_R, pts_Z, pts_psin, np.column_stack([Rb, Zb]),
            lambda p: self._psin_to_rhot(p, self.eq, self.time_jet),
        )

    def rhot_on_grid(self, R, Z):
        Rg, Zg = np.meshgrid(R, Z, indexing='ij')
        rhot, inside, psin = self.rhot_at_points(Rg.ravel(), Zg.ravel())
        return (rhot.reshape(Rg.shape), inside.reshape(Rg.shape),
                psin.reshape(Rg.shape))

    def native_rhot(self):
        eq, tind_eq = self.eq, self.tind
        Rg_n = np.asarray(eq._psirzmg[0])
        Zg_n = np.asarray(eq._psirzmg[1])
        psin_native = np.asarray(eq._psi_norm[tind_eq]).reshape(Rg_n.shape)
        psin_clip = np.clip(np.where(np.isnan(psin_native), 1.0, psin_native),
                            0.0, 1.0)
        rhot = self._psin_to_rhot(
            psin_clip.ravel(), eq, self.time_jet).reshape(Rg_n.shape)
        Rb, Zb = self.lcfs()
        mask = MplPath(np.column_stack([Rb, Zb])).contains_points(
            np.column_stack([Rg_n.ravel(), Zg_n.ravel()])).reshape(Rg_n.shape)
        return Rg_n, Zg_n, np.where(mask, rhot, np.nan)


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


def _running_median3(a):
    """3-point running median; endpoints unchanged. Cosmetic de-speckle."""
    a = np.asarray(a, dtype=float)
    if a.size < 3:
        return a.copy()
    out = a.copy()
    out[1:-1] = np.median(np.stack([a[:-2], a[1:-1], a[2:]]), axis=0)
    return out


def _subgrid_bin(rhot_c, half, weight, edges, density=None):
    """Anti-aliased binning: spread each sample over [rhot-half, rhot+half].

    Shared by the box LOS path (:func:`los_shell_fraction`) and the real-cell
    detector path (:func:`los_file_detector_rate`). A sample at ``rhot_c`` with
    half-range ``half`` distributes its ``weight`` across the rhot interval and
    adds to each bin in proportion to its overlap with that interval. Total
    weight is conserved (zero-half samples drop entirely into their containing
    bin). Vectorised with ``np.add.at``.

    Parameters
    ----------
    density : 1-D array (len = n_bins) or None
        Optional per-bin density for a **Jacobian-weighted** distribution.
        With ``density=None`` (default) the sample is spread *uniformly in
        rhot* (overlap-weighted) -- the original behaviour, used by the box
        path. When ``density`` is given the within-sample weight goes as
        ``overlap * density_bin`` (re-normalised per sample so the total is
        unchanged). Passing ``density = DVOL`` makes the split follow the
        flux-shell volume, which is the physically correct distribution of a
        Cartesian cell's volume across shells near the axis (where
        ``dV/drhot -> 0`` so a cell holds more volume in its outer-rhot part).
        This removes the near-axis over-fill/starve artifact in ``f_det``
        without the geometric path's "pin enclosed shells to 1" trick (which
        is inapplicable because ``f_det`` is not unity there).

    Returns
    -------
    binned weight array of shape ``(len(edges) - 1,)``.
    """
    rhot_c = np.asarray(rhot_c, dtype=float)
    half = np.asarray(half, dtype=float)
    weight = np.asarray(weight, dtype=float)
    n_bins = len(edges) - 1
    out = np.zeros(n_bins)
    if rhot_c.size == 0:
        return out

    # Clip to [0, 1] so weight isn't spilled outside the physical range.
    rhot_lo = np.clip(rhot_c - half, 0.0, 1.0)
    rhot_hi = np.clip(rhot_c + half, 0.0, 1.0)
    span = rhot_hi - rhot_lo

    i_lo = np.clip(np.searchsorted(edges, rhot_lo, side='right') - 1,
                   0, n_bins - 1)
    i_hi = np.clip(np.searchsorted(edges, rhot_hi, side='right') - 1,
                   0, n_bins - 1)
    n_per_cell = i_hi - i_lo + 1

    total = int(n_per_cell.sum())
    flat_cell = np.repeat(np.arange(rhot_c.size), n_per_cell)
    starts = np.concatenate(([0], np.cumsum(n_per_cell[:-1])))
    within_offset = np.arange(total) - np.repeat(starts, n_per_cell)
    flat_bin = np.repeat(i_lo, n_per_cell) + within_offset

    cell_lo_b = rhot_lo[flat_cell]
    cell_hi_b = rhot_hi[flat_cell]
    cell_span_b = span[flat_cell]
    cell_w_b = weight[flat_cell]
    bin_lo = edges[flat_bin]
    bin_hi = edges[flat_bin + 1]
    ovl = np.maximum(
        0.0, np.minimum(cell_hi_b, bin_hi) - np.maximum(cell_lo_b, bin_lo))

    # Per-(cell, bin) contribution before normalisation.
    if density is None:
        contrib = ovl
    else:
        contrib = ovl * np.asarray(density, dtype=float)[flat_bin]

    # Per-cell normalisation so each sample's total weight is preserved. For
    # density=None this norm equals the (clipped) span, recovering the original
    # ovl/span formula exactly.
    norm = np.zeros(rhot_c.size)
    np.add.at(norm, flat_cell, contrib)
    norm_b = norm[flat_cell]

    zero_span = cell_span_b == 0
    good = (~zero_span) & (norm_b > 0.0)
    # Default: zero-span samples deposit their full weight in their one bin.
    w = np.where(zero_span, cell_w_b, 0.0)
    w = np.where(good, cell_w_b * contrib / np.where(norm_b > 0.0, norm_b, 1.0), w)
    # Degenerate span>0 but zero density over the touched bins: fall back to
    # the uniform (overlap/span) split so weight is never lost.
    bad = (~zero_span) & (norm_b <= 0.0)
    if np.any(bad):
        w = np.where(
            bad,
            cell_w_b * ovl / np.where(cell_span_b > 0.0, cell_span_b, 1.0),
            w,
        )
    np.add.at(out, flat_bin, w)
    return out


def los_shell_fraction(rhot_los, inside, R, Z, x_rhot, dvol_x_m3,
                       subgrid=True, rmag=None):
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
    rmag      : float or None. Magnetic-axis R [m]. When given and inside the
                LOS R-band, every flux shell whose full poloidal extent stays
                within [R[0], R[-1]] is swept in its entirety by the chord, so
                ``f == 1`` exactly there; those inner bins are set to 1
                analytically (see "Enclosed-shell correction" below).

    Enclosed-shell correction
    -------------------------
    Near the axis the flux-shell volume Jacobian dV/drhot -> 0, while a single
    Cartesian (R, Z) cell spans many rhot bins. Both binning modes then
    mis-share that volume among the innermost bins (over-filling the centre,
    starving the next ring), producing spurious f < 1 dips even though those
    shells are entirely inside the R-band and therefore fully seen. We compute
    ``rhot_crit`` = the rhot of the innermost flux surface that reaches either
    R boundary (min over Z, inside the LCFS, of rhot at the Rmin/Rmax columns)
    and pin ``f = 1`` for every bin lying fully inside it.

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
        half = (0.5 * np.abs(grad_R) * dR[:, None]
                + 0.5 * np.abs(grad_Z) * dZ[None, :])
        keep = cell_vol.ravel() > 0.0
        los_vol_bin = _subgrid_bin(
            rhot_los.ravel()[keep], half.ravel()[keep],
            cell_vol.ravel()[keep], edges,
        )

    dvol = np.asarray(dvol_x_m3, dtype=float)
    with np.errstate(invalid='ignore', divide='ignore'):
        f = np.where(dvol > 0.0, los_vol_bin / dvol, 0.0)
    f = np.clip(f, 0.0, 1.0)

    # Enclosed-shell correction: pin f = 1 for flux shells fully inside the
    # R-band (see docstring). Only meaningful when the axis sits in the band;
    # otherwise the innermost shells lie outside the chord (f -> 0, not 1).
    if rmag is not None and R_arr[0] <= rmag <= R_arr[-1]:
        col_lo = rhot_los[0, :][inside[0, :]]    # Rmin column, inside LCFS
        col_hi = rhot_los[-1, :][inside[-1, :]]  # Rmax column, inside LCFS
        cols = np.concatenate([col_lo, col_hi])
        if cols.size:
            rhot_crit = float(cols.min())
            enclosed = edges[1:] <= rhot_crit    # bin upper edge inside band
            f = np.where(enclosed, 1.0, f)
    return f, los_vol_bin


# -----------------------------------------------------------------------------
# Detector-weighted thermal rate from the real KM14 LOS cell file
# -----------------------------------------------------------------------------
def los_file_detector_rate(cells, eqs, xs_sorted, thntx_sorted,
                           thntx_to_si, dvol_m3, subgrid=True):
    """Thermal-neutron quantities from the real KM14 LOS cell cloud.

    Each cell carries its own volume ``V`` and etendue weight ``C`` (with
    ``C = V * Omega/4pi``, ``Omega`` = detector solid angle from the cell). For
    the isotropic thermal source the *detector-reaching* rate is therefore

        Rate_det   = sum_cells  eps(rhot_cell) * C_cell      [n/s]

    and the bare real-chord emission (no solid-angle weighting, for comparison
    with the idealised ``Rate_LOS``) is

        Rate_chord = sum_cells  eps(rhot_cell) * V_cell      [n/s]

    The detector response is also binned onto the TRANSP flux shells to give a
    LOS detector-coupling weight ``f_det(rhot) = C_bin / DVOL`` (dimensionless,
    ~Omega/4pi times the swept volume fraction) and the LOS-weighted emissivity
    ``THKM14_det = THNTX * f_det`` with the same closure as the geometric path:
    ``sum(THNTX * DVOL * f_det) == Rate_det``.

    Parameters
    ----------
    cells       : dict from :func:`read_los_file`
    eqs         : equilibrium object exposing ``rhot_at_points``
    xs_sorted   : TRANSP rhot grid (ascending)
    thntx_sorted: THNTX on ``xs_sorted`` (native CGS units)
    thntx_to_si : CGS->SI factor for THNTX
    dvol_m3     : TRANSP DVOL on ``xs_sorted`` in m^3

    Returns
    -------
    dict with keys
        rate_det, rate_chord, rhot_cells, inside_cells, eps_si,
        c_bin, v_bin, f_det, thkm14_det_si, rate_from_profile,
        cum_rhot, cum_frac, rho_med, vol_tot, c_tot, n_inside
    """
    R = cells["R"]
    Z = cells["Z"]
    C = np.asarray(cells["C"], dtype=float)
    V = np.asarray(cells["V"], dtype=float)

    rhot_c, inside_c, _ = eqs.rhot_at_points(R, Z)

    # Emissivity at each cell (n/m3/s), zeroed outside the LCFS. Reuse the same
    # padded-PCHIP interpolation used for the box grid (works on any shape).
    eps_si = thntx_on_grid(xs_sorted, thntx_sorted, rhot_c, inside_c) * thntx_to_si

    sig = eps_si * C                       # per-cell detector-reaching rate
    rate_det = float(np.sum(sig))
    rate_chord = float(np.sum(eps_si * V))

    # Bin the etendue / volume onto the TRANSP flux shells (same edges as
    # los_shell_fraction: midpoints of xs, plus 0 and 1).
    edges = np.concatenate(([0.0], 0.5 * (xs_sorted[:-1] + xs_sorted[1:]), [1.0]))
    w_in = inside_c.astype(float)
    rhot_crit = None
    if not subgrid:
        c_bin, _ = np.histogram(rhot_c, bins=edges, weights=C * w_in)
        v_bin, _ = np.histogram(rhot_c, bins=edges, weights=V * w_in)
    else:
        # Anti-aliased binning. The LOS cells are Cartesian while the flux
        # surfaces are ~cylindrical, so one cell spans several TRANSP rhot bins
        # (mostly in Z) -> point binning aliases f_det (CONTEXT.md). Give each
        # cell a rhot half-span 0.5*(|d rhot/dR|*dR + |d rhot/dZ|*dZ) and spread
        # its weight over that range (same scheme as the box subgrid path).
        # Local |grad rhot| comes from a dense structured rhot field; per-cell
        # sizes are dZ = z-slice spacing and dR = sqrt(V/dZ) (a cell is ~one
        # slice thick, so V ~= dR*dx*dZ with dR ~= dx in the poloidal plane).
        from scipy.interpolate import RegularGridInterpolator
        nrg, nzg = 200, 600
        Rg1 = np.linspace(float(R.min()), float(R.max()), nrg)
        Zg1 = np.linspace(float(Z.min()), float(Z.max()), nzg)
        rhot_g, inside_g, _ = eqs.rhot_on_grid(Rg1, Zg1)   # (nrg, nzg), 'ij'
        gR, gZ = np.gradient(rhot_g, Rg1, Zg1)
        igR = RegularGridInterpolator((Rg1, Zg1), np.abs(gR),
                                      bounds_error=False, fill_value=0.0)
        igZ = RegularGridInterpolator((Rg1, Zg1), np.abs(gZ),
                                      bounds_error=False, fill_value=0.0)
        keep = inside_c & (C > 0.0)
        pts = np.column_stack([R[keep], Z[keep]])
        dZc = float(np.median(np.diff(np.unique(Z))))
        dRc = np.sqrt(np.maximum(V[keep], 0.0) / dZc)
        half = 0.5 * (igR(pts) * dRc + igZ(pts) * dZc)
        # DVOL-weighted (Jacobian) split: a Cartesian cell's volume is spread
        # across the shells it spans in proportion to each shell's volume, not
        # uniformly in rhot. This is the physically correct distribution and
        # removes the near-axis over-fill/starve dip in f_det (the geometric
        # path's "pin enclosed shells to 1" fix does not apply -- f_det is not
        # unity on the plateau). Total C / V are conserved.
        c_bin = _subgrid_bin(rhot_c[keep], half, C[keep], edges, density=dvol_m3)
        v_bin = _subgrid_bin(rhot_c[keep], half, V[keep], edges, density=dvol_m3)

        # rhot_crit = innermost flux surface that reaches an R-boundary of the
        # cell footprint (same construction as the geometric enclosed-shell
        # correction). Shells inside it are fully swept by the chord.
        if Rg1[0] <= eqs.Rmag <= Rg1[-1]:
            col_lo = rhot_g[0, :][inside_g[0, :]]
            col_hi = rhot_g[-1, :][inside_g[-1, :]]
            cols = np.concatenate([col_lo, col_hi])
            if cols.size:
                rhot_crit = float(cols.min())

    with np.errstate(invalid="ignore", divide="ignore"):
        f_det = np.where(dvol_m3 > 0.0, c_bin / dvol_m3, 0.0)

    # Enclosed-shell flattening. For shells fully inside the chord footprint
    # (rhot <= rhot_crit) the wide R-band captures the whole poloidal cross
    # section, so the chord samples a constant toroidal+solid-angle fraction and
    # f_det is physically *flat* (verified: 0.94-1.04 x plateau over
    # [0.018, rhot_crit]). The few innermost bins are nonetheless under-sampled
    # (the ~13 mm cells barely resolve the tiny axis tube, where DVOL -> 0), so
    # f_det dips/rises there. Pin the enclosed region to its DVOL-weighted mean
    # -- the detector analogue of the geometric path pinning f = 1 there, but to
    # the real plateau value since f_det != 1. The volume weighting naturally
    # down-weights the corrupt tiny-DVOL inner bins, and Sum(DVOL*f_det) over
    # the region is unchanged so Rate_det closure is preserved.
    if rhot_crit is not None:
        enclosed = edges[1:] <= rhot_crit
        dsum = float(dvol_m3[enclosed].sum())
        if enclosed.any() and dsum > 0.0:
            f_det = np.where(enclosed, float(c_bin[enclosed].sum()) / dsum, f_det)
        # Light 3-bin running median *outside* rhot_crit only, to suppress the
        # ~1-2% finite-cell sampling wiggle (the physical bump/roll-off survives;
        # the flat enclosed plateau is left untouched). Cosmetic -- Rate_det
        # (= sum eps*C) is computed directly from the cells and is unchanged.
        if subgrid:
            outside = edges[1:] > rhot_crit
            f_det = np.where(outside, _running_median3(f_det), f_det)

    thntx_si = thntx_sorted * thntx_to_si
    thkm14_det_si = thntx_si * f_det
    # Closure check: profile-reconstructed rate should match Rate_det.
    rate_from_profile = float(np.sum(thntx_si * dvol_m3 * f_det))

    # Cumulative detector-signal fraction vs rhot -> signal-median radius.
    order = np.argsort(rhot_c)
    cum_rhot = rhot_c[order]
    cum_sig = np.cumsum(sig[order])
    tot = cum_sig[-1]
    cum_frac = cum_sig / tot if tot > 0 else np.zeros_like(cum_sig)
    rho_med = float(np.interp(0.5, cum_frac, cum_rhot)) if tot > 0 else float("nan")

    return {
        "rate_det": rate_det,
        "rate_chord": rate_chord,
        "R_cells": R,
        "Z_cells": Z,
        "C_cells": C,
        "rhot_cells": rhot_c,
        "inside_cells": inside_c,
        "eps_si": eps_si,
        "c_bin": c_bin,
        "v_bin": v_bin,
        "f_det": f_det,
        "thkm14_det_si": thkm14_det_si,
        "rate_from_profile": rate_from_profile,
        "cum_rhot": cum_rhot,
        "cum_frac": cum_frac,
        "rho_med": rho_med,
        "rhot_crit": rhot_crit,
        "vol_tot": float(np.sum(V * w_in)),
        "c_tot": float(np.sum(C * w_in)),
        "n_inside": int(np.count_nonzero(inside_c)),
    }


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
    run_dir = bzi.find_run_dir(args.pulse, args.runid, args.data_dir)
    cdf_path = run_dir / f"{run_id}.CDF"
    th_var, chan_label = CHANNELS[args.channel]
    print(f'Run dir : {run_dir}')
    print(f'Channel : {chan_label}  (thermal var = {th_var})')

    # ------- Time / thermal slice from the main CDF ---------------------
    t_jet_grid = read_time_grid(cdf_path) + 40.0
    if args.time is None:
        time_jet = float(t_jet_grid[len(t_jet_grid) // 2])
        print(f'No time supplied; using midpoint t_JET = {time_jet:.3f} s.')
    else:
        time_jet = float(args.time)
    trind = int(np.abs(t_jet_grid - time_jet).argmin())

    (thntx_x, x_rhot, dvol_x, thntx_to_si,
     dvol_to_m3, unit_th) = read_thermal_slice(cdf_path, th_var, trind)

    print(f'TRANSP slice [{trind}]: t_JET = {t_jet_grid[trind]:.3f} s '
          f'(t_TRANSP = {t_jet_grid[trind] - 40.0:.3f} s)')
    print(f'  {th_var} units = [{unit_th}],  conversion to n/m3/s = '
          f'{thntx_to_si:.1e}')
    print(f'  {th_var}(rhot=0) = {thntx_x[0]:.3e} [{unit_th}]  '
          f'= {thntx_x[0] * thntx_to_si:.3e} n/m3/s')

    # ------- Equilibrium ------------------------------------------------
    if args.eq_source == 'ppf':
        print(f'Equilibrium: PPF {args.dda}/{args.uid}/{args.seq} '
              f'(importing ppf) ...')
        eqs = EqPPF(args.pulse, args.dda, args.uid, args.seq, time_jet)
    else:
        print('Equilibrium: TRANSP CDF (self-contained, no ppf) ...')
        eqs = EqCDF(cdf_path, time_jet)
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
            subgrid=not args.no_subgrid,
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
