#!/usr/bin/env python3
"""los_th_bt_ratio.py -- KM14 line-of-sight thermal/beam-target neutron ratio.

Computes the TH/BT neutron ratio that the KM14 vertical camera sees, for a
given TRANSP run, so it can be compared against the spectroscopically-measured
KM14 TH/BT. This is the beam-target analogue / extension of
``los_thermal_rate.py``, but **self-contained on the TRANSP CDFs** (no ppf /
profiles.Eq dependency): the equilibrium, thermal profile and beam-target
emissivity all come from the run's own CDF + _fi/_neut files, so it runs in the
local dev environment as well as on freia.

What it does
------------
1. Builds the KM14 chord (R, Z) grid (vertical chord, R in [Rmin, Rmax]).
2. Maps (R, Z) -> rhot from the TRANSP equilibrium: psi(R, Z) from ``PSIRZ``
   (C-order (Z, R)), normalized by the axis/boundary poloidal flux, then
   rhot(psi_n) inverted from ``PLFLX`` vs ``XB``.
3. Thermal emissivity eps_TH(rhot): the flux-function ``THNTX_<chan>`` profile.
4. Beam-target emissivity eps_BT(rhot): flux-surface average (BMVOL-weighted
   per x-row) of the per-zone ``BTN4``/``BTN1`` from ``_neut``. (Commit 3
   showed the *emission direction* is isotropic to <0.5% for KM14, so the 4*pi
   per-zone rate is the right thing to integrate; the remaining spatial
   poloidal asymmetry is optionally checked with ``--bt-mode zone``.)
5. Integrates both along the chord. For the **ratio** the chord geometry
   (toroidal width, solid angle, 2*pi*R) cancels, so

       (TH/BT)_LOS = integral_chord eps_TH dR dZ / integral_chord eps_BT dR dZ

   is robust. Reported alongside the full-torus (volume) ratio, whose value
   matches the 0D ``THNTX*DVOL`` / ``BTN*BMVOL`` totals as a cross-check.

Channels
--------
``--channel total`` (default): the unseparated thermal ``THNTX`` (DD+DT+...)
versus the sum of every BT component present in ``_neut`` (DD+DT+TT+TD,
whichever exist) -- the unfolded total neutron rate KM14 sees if no
spectroscopic separation is applied.

``--channel dd`` and ``--channel dt`` keep the spectroscopically-separated
DD-only (THNTX_DD vs BTN4, 2.45 MeV) and DT-only (THNTX_DT vs BTN1, 14 MeV)
options for comparison against KM14's separated yields.

Beam-target emission anisotropy (``--bt-aniso``)
------------------------------------------------
By default BT emission is treated as isotropic (the 4*pi NUBEAM ``BTN*`` rate),
which commit 3 (``bt_los_emissivity.py``) showed is good to <0.5% for a
perfectly vertical KM14 view. ``--bt-aniso`` instead accounts for the real
detector geometry: KM14 is a point ~10 m above ``Zmag`` on the chord axis, so
neutrons reach it within a narrow upward cone whose axis tilts a couple of
degrees off vertical across the chord. Each zone's BTN rate is then weighted by
the directional factor ``g(zone) = (dEps/dOmega toward the detector)/(Eps/4pi)``,
computed by MC-sampling the fast-ion distribution + CM->lab boost (reusing the
``bt_los_emissivity`` machinery) with a per-zone line-of-sight pointing at the
detector. The 4*pi *magnitude* still comes from NUBEAM's validated ``BTN*``;
only the dimensionless ``g`` is taken from the MC. In the ratio the common
detector solid angle cancels, so the net effect is exactly
``(TH/BT)_LOS -> (TH/BT)_LOS / <g>_BT``. Detector at ``--detector-R`` (default
chord centre 2.90 m), ``--detector-height`` above Zmag (default 10 m).

``--los-file _KM3.los`` replaces that point-detector direction model with the
**exact** per-cell emission versors ``(u, v, w)`` from M. Nocente's KM14 LOS
file: the versor (detector Cartesian x=toroidal, y=radial, z=vertical) maps to
the ``(R, phi, Z)`` MC basis as ``(v, u, w)`` and is interpolated to each zone
in the LOS footprint (others keep the point-detector estimate). On 104614 M29
the real directions differ from the point-detector model by only ~0.6 deg
(median), so ``<g>_BT`` is essentially unchanged (1.0044 -> 1.0042) -- the input
is now exact rather than approximated.

``--los-file`` also (independently of ``--bt-aniso``) adds the real-LOS detector
weighting to the plot: the per-shell etendue coupling ``f_det = C_bin/DVOL``
(reused from ``los_thermal_rate``, emissivity-independent) is applied to both TH
and BT for a C-weighted cumulative ``Sum(TH f_det DVOL)/Sum(BT f_det DVOL)``,
overlaid on the geometric-LOS and full-plasma cumulatives, plus the real LOS
cell footprint over the TH/BT(local) map and ``f_det`` on the weight panel.

CLI
---
    python los_th_bt_ratio.py 104614 M30 --idx 1 --plot          # total (default)
    python los_th_bt_ratio.py 104614 M30 --channel dd --bt-mode zone
    python los_th_bt_ratio.py 104614 M29 --channel dd --bt-aniso # finite-height view

Flags: --idx --data-dir --channel {total,dd,dt} --bt-mode {flux,zone}
       --Rmin --Rmax --wtor --nR --nZ --plot --save --no-plot
       --bt-aniso --detector-R --detector-height --los-file --cone-deg --nsamp
       --fast-norm {bdens,ntot} --no-rotation --seed
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
from netCDF4 import Dataset
from matplotlib.path import Path as MplPath
from scipy.interpolate import PchipInterpolator, griddata

import bt_zone_integrator as bzi

R_MIN_DEFAULT = 2.70
R_MAX_DEFAULT = 3.10
W_TOR_DEFAULT = 0.40
NR_DEFAULT = 120
NZ_DEFAULT = 200

# KM14 collimator geometry for the optional BT-anisotropy treatment: a point
# detector above the plasma, on the vertical chord axis. Neutrons reach it
# travelling up a narrow cone (a slightly different axis from each emission
# point, since the detector is at a finite height).
DET_R_DEFAULT = 2.90          # m, above the chord centre (R in [2.70, 3.10])
DET_HEIGHT_DEFAULT = 10.0     # m above Zmag

# (thermal profile var, list of beam-target _neut keys to sum, label).
# An empty bt-key list means "sum every BT component present in _neut".
CHANNELS = {
    "total": ("THNTX",    [],            "total (DD+DT+TT+TD)"),
    "dd":    ("THNTX_DD", ["DD"],        "DD (2.45 MeV)"),
    "dt":    ("THNTX_DT", ["DT"],        "DT (14 MeV)"),
}


# ----------------------------------------------------------------------
# Equilibrium from the TRANSP CDF (psi(R,Z) -> rhot), no ppf
# ----------------------------------------------------------------------

class CdfEquilibrium:
    """rhot(R, Z) from the TRANSP main CDF at a chosen time."""

    def __init__(self, cdf_path, time):
        with Dataset(cdf_path, "r") as d:
            t3 = np.array(d["TIME3"][:])
            ti = int(np.argmin(np.abs(t3 - time)))
            self.tind = ti
            self.t_used = float(t3[ti])
            self.RG = np.array(d["RGRID"][ti]) / 100.0      # m
            self.ZG = np.array(d["ZGRID"][ti]) / 100.0
            nR, nZ = len(self.RG), len(self.ZG)
            psi = np.array(d["PSIRZ"][ti]).reshape(nZ, nR)  # [iz, ir] Wb/rad
            psi0 = float(d["PSI0_TR"][ti])
            psibnd = float(d["PLFLXA"][ti])
            XB = np.array(d["XB"][ti])                       # rhot boundaries
            PLFLX = np.array(d["PLFLX"][ti])                 # poloidal flux on XB
            self.Rmag = float(d["RAXIS"][ti]) / 100.0        # magnetic axis [m]
            self.Zmag = float(d["YAXIS"][ti]) / 100.0
        self.psin = (psi - psi0) / (psibnd - psi0)           # normalized psi(R,Z)
        # rhot(psi_n): psi_n on the XB grid is monotonic 0->1; XB is rhot
        psin_b = (PLFLX - PLFLX[0]) / (PLFLX[-1] - PLFLX[0])
        self._psin_b = psin_b
        self._rhot_b = XB

    def rhot(self, R, Z):
        """Bilinear psi_n at (R, Z) [m] -> rhot (sqrt tor flux). Clipped [0,1].

        NB: bilinear floors psi_n at the coarse PSIRZ axis node (~3e-4 -> rhot
        ~0.02), which zeroes the innermost LOS flux shells. Use ``rhot_pinned``
        for the chord grid; this plain version is kept for ad-hoc point lookups.
        """
        RG, ZG = self.RG, self.ZG
        ir = np.clip(np.searchsorted(RG, R) - 1, 0, len(RG) - 2)
        iz = np.clip(np.searchsorted(ZG, Z) - 1, 0, len(ZG) - 2)
        tr = (R - RG[ir]) / (RG[ir + 1] - RG[ir])
        tz = (Z - ZG[iz]) / (ZG[iz + 1] - ZG[iz])
        P = self.psin
        pn = (P[iz, ir] * (1 - tr) * (1 - tz) + P[iz, ir + 1] * tr * (1 - tz)
              + P[iz + 1, ir] * (1 - tr) * tz + P[iz + 1, ir + 1] * tr * tz)
        pn = np.clip(pn, 0.0, 1.0)
        return np.clip(np.interp(pn, self._psin_b, self._rhot_b), 0.0, 1.0)

    def rhot_pinned(self, R, Z, Rb=None, Zb=None):
        """rhot at (R, Z) [m] via cubic griddata of psi_n with the magnetic axis
        pinned at psi_n = 0 and (when ``Rb``/``Zb`` given) the LCFS at psi_n = 1.

        Bilinear interpolation of the coarse PSIRZ grid (dR ~ 2 cm) cannot reach
        the true axis minimum (psi is paraboloidal with its vertex between
        nodes), flooring psi_n at ~3e-4, i.e. rhot ~ 0.02 -- so the innermost
        f(rhot) bins get no LOS volume. Injecting the axis as a scattered node
        removes the floor. Symmetrically, the cubic fit also floors *short* of
        psi_n = 1 inside the LCFS (max ~0.999), starving the outermost f(rhot)
        bins so the weight craters near rhot = 1; passing the LCFS polygon
        (``Rb``, ``Zb`` [m], the same one used for the inside mask) pins it at
        psi_n = 1 and removes that crater. Mirrors
        ``los_thermal_rate.EqCDF.rhot_on_grid``.
        """
        from scipy.interpolate import griddata
        Rg_n, Zg_n = np.meshgrid(self.RG, self.ZG)  # (nZ, nR), matches psin[iz,ir]
        pin_R = [self.Rmag]
        pin_Z = [self.Zmag]
        pin_psin = [0.0]
        if Rb is not None and Zb is not None:
            Rb = np.asarray(Rb).ravel()
            Zb = np.asarray(Zb).ravel()
            pin_R = np.concatenate([pin_R, Rb])
            pin_Z = np.concatenate([pin_Z, Zb])
            pin_psin = np.concatenate([pin_psin, np.ones(Rb.size)])
        pts_R = np.concatenate([Rg_n.ravel(), pin_R])
        pts_Z = np.concatenate([Zg_n.ravel(), pin_Z])
        pts_psin = np.concatenate([self.psin.ravel(), pin_psin])
        query = np.column_stack([np.asarray(R).ravel(), np.asarray(Z).ravel()])
        pd = griddata((pts_R, pts_Z), pts_psin, query, method="cubic")
        if np.any(np.isnan(pd)):
            pl = griddata((pts_R, pts_Z), pts_psin, query, method="linear")
            pd = np.where(np.isnan(pd), pl, pd)
        pc = np.clip(np.where(np.isnan(pd), 1.0, pd), 0.0, 1.0)
        return np.clip(np.interp(pc, self._psin_b, self._rhot_b), 0.0, 1.0)


def read_lcfs(fi_path):
    """Outermost flux surface (LCFS) polygon (R, Z) [m] from _fi RSURF/ZSURF."""
    with Dataset(fi_path, "r") as d:
        R = np.array(d["RSURF"][:])[-1, :] / 100.0
        Z = np.array(d["ZSURF"][:])[-1, :] / 100.0
    return R, Z


def read_thermal_profile(cdf_path, var, tind):
    """Return (rhot grid X, profile) for a THNTX_* thermal emissivity [N/CM3/SEC]."""
    with Dataset(cdf_path, "r") as d:
        X = np.array(d["X"][tind])
        prof = np.array(d[var][tind])
    return X, prof


def read_scalar_totals(cdf_path, tind):
    with Dataset(cdf_path, "r") as d:
        out = {}
        for k in ("BTNTS_DD", "BTNTS_DT", "BTNTS_TT", "BTNTS_TD"):
            if k in d.variables:
                out[k] = float(d[k][tind])
    return out


# ----------------------------------------------------------------------
# Optional beam-target emission anisotropy (KM14 finite-height detector)
# ----------------------------------------------------------------------

# _neut beam-target key -> (bt_kinematics reaction, thermal-target mass amu,
# thermal-density key in bzi.read_thermal_profiles). Keys absent here (e.g.
# "TT", which has no Bosch-Hale sigma) are treated as isotropic (g = 1).
def _bt_reaction_map():
    import bt_kinematics as btk
    return {
        "DD": ("DD_n", btk.M_D_AMU, "nd"),   # beam-D + thermal-D  (2.45 MeV)
        "DT": ("DT",   btk.M_T_AMU, "nt"),   # beam-D + thermal-T  (14 MeV)
    }


def los_directions_for_zones(los_path, Rz, Zz, max_cells=30000, seed=0):
    """Per-zone emission direction n_los (R, phi, Z) from the real KM14 LOS file.

    Reads ``_KM3.los`` and interpolates its per-cell emission versor onto each
    NUBEAM zone position ``(Rz, Zz)`` [m]. The file's versor ``(u, v, w)`` is in
    the detector Cartesian frame (x = toroidal, y = radial since y ~ R,
    z = vertical), so in the anisotropy MC's ``(R, phi, Z)`` basis it maps to
    ``(v, u, w)`` -- this also restores the small *toroidal* tilt (u) that the
    point-detector model zeroes.

    The cell cloud is huge and many cells share nearly the same (R, Z) at
    different toroidal x; ``griddata`` linear over a random subsample averages
    those naturally. Zones outside the LOS footprint interpolate to NaN.

    Returns ``(n_los_zone, valid)``: ``n_los_zone`` is (n_zone, 3) unit vectors
    (NaN rows where invalid) and ``valid`` is the boolean in-footprint mask; the
    caller falls back to the point-detector direction for invalid zones.
    """
    from los_thermal_rate import read_los_file
    cells = read_los_file(los_path)
    Rc, Zc = cells["R"], cells["Z"]
    u, v, w = cells["u"], cells["v"], cells["w"]
    if Rc.size > max_cells:
        sel = np.random.default_rng(seed).choice(Rc.size, max_cells, replace=False)
        Rc, Zc, u, v, w = Rc[sel], Zc[sel], u[sel], v[sel], w[sel]
    pts = np.column_stack([Rc, Zc])
    q = np.column_stack([np.asarray(Rz, float).ravel(),
                         np.asarray(Zz, float).ravel()])
    nr = griddata(pts, v, q, method="linear")    # R-component  <- y-versor v
    nph = griddata(pts, u, q, method="linear")   # phi-component <- x-versor u
    nz = griddata(pts, w, q, method="linear")    # Z-component  <- z-versor w
    valid = ~(np.isnan(nr) | np.isnan(nph) | np.isnan(nz))
    n_los = np.column_stack([nr, nph, nz])
    nrm = np.sqrt(np.nansum(n_los ** 2, axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        n_los = n_los / np.where(nrm > 0, nrm, np.nan)[:, None]
    return n_los, valid


def bt_anisotropy_factors(fi, cdf_path, bt_keys, det_R, det_Z, args):
    """Per-zone anisotropy factor g(zone) = (dEps/dOmega)/(Eps_4pi/4pi) along the
    line from each zone up to the KM14 detector, for each handled BT component.

    The detector is a point at (``det_R``, ``det_Z``) [m] on the chord toroidal
    plane, so the LOS from a zone at (R, Z) is ``(det_R - R, 0, det_Z - Z)``
    (zero toroidal component) in the (R, phi, Z) basis -- a narrow upward cone
    whose axis tilts a couple of degrees off vertical across the chord. The 4*pi
    magnitude is *not* taken from this MC (NUBEAM's BTN* is the validated
    reactivity); only the dimensionless directional factor g is used, so the
    caller multiplies the per-zone BTN rate by g.

    Returns ``(g_map, info)`` where ``g_map[key]`` is a per-zone array (1.0 for
    keys without a kinematics mapping) and ``info`` holds reporting diagnostics.
    """
    import bt_kinematics as btk
    import bt_los_emissivity as ble

    therm = bzi.read_thermal_profiles(cdf_path, fi["time"], fi["x2d"])
    bfield = ble.BField(cdf_path, fi["time"])
    n_fast = bzi.fast_ion_density(fi)
    if args.fast_norm == "bdens":
        n_fast, _ = bzi.renormalize_to_bdens(n_fast, fi["x2d"], fi["bmvol"],
                                             therm["bdens_d"])

    Rz = fi["r2d"] / 100.0
    Zz = fi["z2d"] / 100.0
    bhat_zone = bfield.bhat(Rz, Zz)                       # (n_zone, 3)
    if args.no_rotation:
        vrot_zone = np.zeros_like(Rz)
    else:
        with Dataset(str(cdf_path), "r") as d:
            t3 = np.array(d["TIME3"][:])
            ti = int(np.argmin(np.abs(t3 - fi["time"])))
            X = np.array(d["X"][ti]); OM = np.array(d["OMEGA"][ti])
        vrot_zone = np.interp(fi["x2d"], X, OM) * Rz       # m/s, toroidal (+phi)

    # per-zone LOS direction (R, phi, Z). Default: point-detector model -- from
    # the zone straight up to the detector point (zero toroidal component).
    dR = det_R - Rz
    dZ = det_Z - Zz
    nrm = np.sqrt(dR * dR + dZ * dZ)
    n_los_zone = np.column_stack([dR / nrm, np.zeros_like(dR), dZ / nrm])

    # If the real KM14 LOS cell file is supplied, replace the point-detector
    # direction by the file's exact (u,v,w) versor wherever the zone falls in
    # the LOS footprint (others keep the point-detector estimate).
    los_info = None
    los_file = getattr(args, "los_file", None)
    if los_file:
        n_real, valid = los_directions_for_zones(los_file, Rz, Zz)
        cosd = np.clip(np.abs(np.einsum("zi,zi->z", n_los_zone,
                                        np.nan_to_num(n_real))), 0.0, 1.0)
        dang = np.degrees(np.arccos(cosd))
        n_los_zone = np.where(valid[:, None], n_real, n_los_zone)
        los_info = dict(
            n_used=int(valid.sum()), n_zone=int(valid.size),
            dang_med=float(np.median(dang[valid])) if valid.any() else float("nan"),
            dang_max=float(np.max(dang[valid])) if valid.any() else float("nan"),
        )

    tilt = np.degrees(np.arctan2(
        np.abs(n_los_zone[:, 0]), np.abs(n_los_zone[:, 2])))   # off-vertical tilt
    ang_bL = np.degrees(np.arccos(np.clip(
        np.abs(np.einsum("zi,zi->z", bhat_zone, n_los_zone)), 0, 1)))

    cos_alpha = np.cos(np.radians(args.cone_deg))
    rng = np.random.default_rng(args.seed)
    rxn = _bt_reaction_map()

    g_map = {}
    gbar = {}
    consistency = {}
    n_zone = fi["F"].shape[0]
    for k in bt_keys:
        if k not in rxn:
            g_map[k] = np.ones(n_zone)            # no sigma (e.g. TT) -> isotropic
            continue
        reaction, m_th, dkey = rxn[k]
        n_th = therm[dkey]
        eps4pi, dedom = ble.zone_emissivity(
            fi, n_fast, n_th, therm["ti"], reaction, btk.M_D_AMU, m_th,
            bhat_zone, vrot_zone, n_los_zone, args.nsamp, cos_alpha, rng)
        iso = eps4pi / (4.0 * np.pi)
        g = np.divide(dedom, iso, out=np.ones_like(dedom), where=iso > 0)
        g_map[k] = g
        w = eps4pi * fi["bmvol"]
        gbar[k] = float(np.sum(g * w) / np.sum(w)) if np.sum(w) > 0 else 1.0
        ref = fi.get("_neut_ref", {}).get(k)  # optional NUBEAM cross-check
        if ref is not None and np.sum(ref * fi["bmvol"]) > 0:
            consistency[k] = float(np.sum(eps4pi * fi["bmvol"])
                                   / np.sum(ref * fi["bmvol"]))
    info = dict(gbar=gbar, consistency=consistency,
                tilt_max=float(tilt.max()),
                ang_bL=(float(ang_bL.min()), float(np.median(ang_bL)),
                        float(ang_bL.max())),
                det_R=det_R, det_Z=det_Z, los=los_info)
    return g_map, info


# ----------------------------------------------------------------------
# Emissivity on the chord grid
# ----------------------------------------------------------------------

def flux_avg_profile(values_zone, x2d, bmvol):
    """BMVOL-weighted flux-surface average per x-row -> (x_rows, profile)."""
    xu = np.unique(x2d)
    prof = np.array([np.sum(values_zone[x2d == xv] * bmvol[x2d == xv])
                     / np.sum(bmvol[x2d == xv]) for xv in xu])
    return xu, prof


def profile_on_grid(x_prof, val_prof, rhot_grid, inside):
    """PCHIP-interpolate a flux profile onto rhot(R, Z); zero outside the LCFS."""
    order = np.argsort(x_prof)
    xs = np.asarray(x_prof)[order]
    ys = np.asarray(val_prof)[order]
    uniq = np.concatenate(([True], np.diff(xs) > 0))
    xs, ys = xs[uniq], ys[uniq]
    if xs[0] > 0.0:
        xs = np.concatenate(([0.0], xs)); ys = np.concatenate(([ys[0]], ys))
    if xs[-1] < 1.0:
        xs = np.concatenate((xs, [1.0])); ys = np.concatenate((ys, [0.0]))
    f = PchipInterpolator(xs, ys, extrapolate=False)
    out = f(rhot_grid)
    out = np.where(np.isnan(out), 0.0, out)
    return np.where(inside, out, 0.0)


def zone_emis_on_grid(values_zone, r2d_m, z2d_m, Rg, Zg, inside):
    """Per-zone scattered emissivity -> chord grid via cubic griddata (asymmetry)."""
    query = np.column_stack([Rg.ravel(), Zg.ravel()])
    g = griddata((r2d_m, z2d_m), values_zone, query, method="cubic")
    g_lin = griddata((r2d_m, z2d_m), values_zone, query, method="linear")
    g = np.where(np.isnan(g), g_lin, g)
    g = np.where(np.isnan(g), 0.0, g).reshape(Rg.shape)
    return np.where(inside & (g > 0), g, 0.0)


def chord_integral(grid_si, R, Z, w_tor):
    """w_tor * integral eps dR dZ  [n/s] (chord), and the R-integral vs Z."""
    inner = np.trapezoid(grid_si, R, axis=0)
    return float(w_tor * np.trapezoid(inner, Z)), inner


def toroidal_integral(grid_si, R, Z):
    """integral 2*pi*R eps dR dZ  [n/s] (full torus over the chord R,Z area)."""
    Rg, _ = np.meshgrid(R, Z, indexing="ij")
    inner = np.trapezoid(grid_si * 2.0 * np.pi * Rg, R, axis=0)
    return float(np.trapezoid(inner, Z))


# ----------------------------------------------------------------------
# LOS weight function f(rhot) and LOS-weighted TH/BT profiles vs rhot
# ----------------------------------------------------------------------

def _interp_flux_to(x_src, y_src, x_dst):
    """PCHIP a flux profile y_src(x_src) onto x_dst, padded to cover [0, 1].

    Same padding convention as ``profile_on_grid``: rhot=0 takes the on-axis
    value, rhot=1 is forced to zero, and points outside the source support are
    set to 0 (never negative).
    """
    order = np.argsort(x_src)
    xs = np.asarray(x_src, float)[order]
    ys = np.asarray(y_src, float)[order]
    uniq = np.concatenate(([True], np.diff(xs) > 0))
    xs, ys = xs[uniq], ys[uniq]
    if xs[0] > 0.0:
        xs = np.concatenate(([0.0], xs)); ys = np.concatenate(([ys[0]], ys))
    if xs[-1] < 1.0:
        xs = np.concatenate((xs, [1.0])); ys = np.concatenate((ys, [0.0]))
    out = PchipInterpolator(xs, ys, extrapolate=False)(np.asarray(x_dst, float))
    return np.clip(np.where(np.isnan(out), 0.0, out), 0.0, None)


def rhot_weight_profiles(eq, rhot, inside, R, Z, x_th, th_prof, x_bt, bt_prof,
                         cdf_path):
    """LOS shell-fraction weight f(rhot) and LOS-weighted TH/BT vs rhot.

    Mirrors ``los_thermal_rate.py``: f(rhot) is the *geometric* fraction of
    each TRANSP flux shell's toroidal volume that falls inside the chord R-band
    (computed once, shared by TH and BT). Everything is put on the TRANSP ``X``
    grid (where THNTX and DVOL live); the BT flux profile is PCHIP-interpolated
    onto it.

    Returns a dict with, on the sorted ``X`` grid ``xs``:
      f            geometric LOS shell fraction
      th_si, bt_si flux-profile emissivities                 [n/m3/s]
      thkm14,btkm14 LOS-weighted emissivities = eps * f       [n/m3/s]
      dvol_m3      TRANSP shell volume                        [m3]
      shell_th, shell_bt  per-shell LOS rate = eps*f*DVOL     [n/s]
      cum_th, cum_bt      cumulative per-shell LOS rate       [n/s]
      ratio_local  th_si/bt_si  (= thkm14/btkm14; f, DVOL cancel)
      ratio_cum    cum_th/cum_bt  (-> (TH/BT)_LOS at rhot=1)
      ratio_cum_plasma  cumSum(TH*DVOL)/cumSum(BT*DVOL), no LOS weight
                        (-> whole-plasma 0D TH/BT at rhot=1)
      sum_th, sum_bt      total per-shell LOS rate (= cum[-1])
    """
    from los_thermal_rate import los_shell_fraction  # lazy: avoids import cycle
    with Dataset(str(cdf_path), "r") as d:
        dvol = np.array(d["DVOL"][eq.tind]) * 1.0e-6   # CM**3 -> m^3

    order = np.argsort(x_th)
    xs = np.asarray(x_th, float)[order]
    dvol_m3 = np.asarray(dvol, float)[order]

    f, _ = los_shell_fraction(rhot, inside, R, Z, xs, dvol_m3, subgrid=True,
                              rmag=eq.Rmag)

    th_si = np.asarray(th_prof, float)[order] * 1.0e6        # n/m3/s
    bt_si = _interp_flux_to(x_bt, bt_prof, xs) * 1.0e6       # n/m3/s

    thkm14 = th_si * f
    btkm14 = bt_si * f
    shell_th = thkm14 * dvol_m3                              # n/s per shell
    shell_bt = btkm14 * dvol_m3
    cum_th = np.cumsum(shell_th)
    cum_bt = np.cumsum(shell_bt)

    # full-plasma (no-LOS) cumulative: same emissivities and DVOL but without
    # the geometric weight f, i.e. cumSum(TH*DVOL)/cumSum(BT*DVOL). Its rhot=1
    # endpoint is the whole-plasma 0D TH/BT; comparing it to ratio_cum shows how
    # narrowing the emission to the KM14 chord reweights the ratio.
    cum_th_plasma = np.cumsum(th_si * dvol_m3)
    cum_bt_plasma = np.cumsum(bt_si * dvol_m3)

    with np.errstate(invalid="ignore", divide="ignore"):
        ratio_local = np.where(bt_si > 0, th_si / bt_si, np.nan)
        ratio_cum = np.where(cum_bt > 0, cum_th / cum_bt, np.nan)
        ratio_cum_plasma = np.where(cum_bt_plasma > 0,
                                    cum_th_plasma / cum_bt_plasma, np.nan)

    return dict(xs=xs, f=f, th_si=th_si, bt_si=bt_si,
                thkm14=thkm14, btkm14=btkm14, dvol_m3=dvol_m3,
                shell_th=shell_th, shell_bt=shell_bt,
                cum_th=cum_th, cum_bt=cum_bt,
                ratio_local=ratio_local, ratio_cum=ratio_cum,
                ratio_cum_plasma=ratio_cum_plasma,
                sum_th=float(cum_th[-1]), sum_bt=float(cum_bt[-1]))


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("pulse", type=int)
    p.add_argument("run_suffix")
    p.add_argument("--idx", type=int, default=None)
    p.add_argument("--data-dir", default=None)
    p.add_argument("--channel", choices=list(CHANNELS), default="total")
    p.add_argument("--bt-mode", choices=["flux", "zone"], default="flux",
                   help="beam-target emissivity: flux-surface average (default) "
                        "or per-zone griddata (keeps poloidal asymmetry)")
    # --- optional BT emission-anisotropy treatment (KM14 finite-height view) ---
    p.add_argument("--bt-aniso", action="store_true",
                   help="weight the beam-target emissivity by the directional "
                        "factor g(R,Z) for neutrons emitted up toward the KM14 "
                        "detector (default: isotropic 4*pi BTN rates)")
    p.add_argument("--detector-R", type=float, default=DET_R_DEFAULT,
                   help="R of the KM14 detector point [m] (default chord centre "
                        f"{DET_R_DEFAULT})")
    p.add_argument("--detector-height", type=float, default=DET_HEIGHT_DEFAULT,
                   help=f"detector height above Zmag [m] (default {DET_HEIGHT_DEFAULT})")
    p.add_argument("--los-file", default=None,
                   help="real KM14 LOS cell file (_KM3.los). With --bt-aniso, "
                        "the per-zone emission direction toward the detector is "
                        "taken from the file's exact (u,v,w) versors instead of "
                        "the point-detector approximation (--detector-R/-height).")
    p.add_argument("--cone-deg", type=float, default=20.0,
                   help="half-angle of the LOS acceptance cone for the g MC estimator")
    p.add_argument("--nsamp", type=int, default=60000,
                   help="MC samples per zone for the anisotropy factor")
    p.add_argument("--fast-norm", choices=["bdens", "ntot"], default="bdens",
                   help="fast-ion density normalization for the anisotropy MC")
    p.add_argument("--no-rotation", action="store_true",
                   help="disable thermal toroidal rotation in the anisotropy MC")
    p.add_argument("--seed", type=int, default=12345)
    p.add_argument("--Rmin", type=float, default=R_MIN_DEFAULT)
    p.add_argument("--Rmax", type=float, default=R_MAX_DEFAULT)
    p.add_argument("--wtor", type=float, default=W_TOR_DEFAULT)
    p.add_argument("--nR", type=int, default=NR_DEFAULT)
    p.add_argument("--nZ", type=int, default=NZ_DEFAULT)
    p.add_argument("--plot", action="store_true")
    p.add_argument("--save", default=None)
    p.add_argument("--no-plot", action="store_true")
    args = p.parse_args(argv)

    run_id = f"{args.pulse}{args.run_suffix}"
    run_dir = bzi.find_run_dir(args.pulse, args.run_suffix, args.data_dir)
    idxs = bzi.list_fbm_indices(run_dir, run_id)
    if not idxs:
        print(f"No {run_id}_fi_*.cdf in {run_dir}", file=sys.stderr)
        return 1
    idx = args.idx if args.idx is not None else idxs[0]
    fi_path = run_dir / f"{run_id}_fi_{idx}.cdf"
    neut_path = run_dir / f"{run_id}_neut_{idx}.cdf"
    cdf_path = run_dir / f"{run_id}.CDF"

    th_var, bt_keys_req, chan_label = CHANNELS[args.channel]

    fi = bzi.read_fi_distribution(fi_path)
    neut = bzi.read_neut_rates(neut_path)
    # Resolve the BT components to sum: explicit list, or (for "total")
    # every BT key the _neut file actually provides.
    available = [k for k in ("DD", "DT", "TT", "TD") if k in neut]
    bt_keys = bt_keys_req if bt_keys_req else available
    missing = [k for k in bt_keys if k not in neut]
    if missing:
        print(f"Beam-target key(s) {missing} not in {neut_path}", file=sys.stderr)
        return 1
    bt_zone_sum = sum(neut[k] for k in bt_keys)
    eq = CdfEquilibrium(cdf_path, fi["time"])
    Rb, Zb = read_lcfs(fi_path)
    x_th, th_prof = read_thermal_profile(cdf_path, th_var, eq.tind)
    totals = read_scalar_totals(cdf_path, eq.tind)

    print(f"Run dir : {run_dir}")
    print(f"FBM idx : {idx}   t = {fi['time']:.4f} s   "
          f"(equilibrium/thermal slice TIME3[{eq.tind}] = {eq.t_used:.4f} s)")
    print(f"Channel : {chan_label}   thermal={th_var}  "
          f"beam-target=BTN[{'+'.join(bt_keys)}]")
    print(f"BT mode : {args.bt_mode}")

    # chord grid. Insert the magnetic axis into the R, Z sampling so a chord
    # sample sits exactly at (Rmag, Zmag) -- needed for the centremost f(rhot)
    # bin (the innermost flux shell is smaller than one cell). Source-agnostic.
    R = np.linspace(args.Rmin, args.Rmax, args.nR)
    Z = np.linspace(float(eq.ZG.min()), float(eq.ZG.max()), args.nZ)
    if R[0] <= eq.Rmag <= R[-1]:
        R = np.unique(np.concatenate([R, [eq.Rmag]]))
    if Z[0] <= eq.Zmag <= Z[-1]:
        Z = np.unique(np.concatenate([Z, [eq.Zmag]]))
    Rg, Zg = np.meshgrid(R, Z, indexing="ij")
    # rhot_pinned (griddata + pinned axis + pinned LCFS) instead of plain
    # bilinear eq.rhot, so f(rhot) reaches ~1 on axis and is not starved near
    # rhot = 1 (see rhot_pinned docstring). Pin the LCFS with the same (Rb, Zb)
    # used for the inside mask so field and mask stay consistent at the edge.
    rhot = eq.rhot_pinned(Rg.ravel(), Zg.ravel(), Rb=Rb, Zb=Zb).reshape(Rg.shape)
    inside = MplPath(np.column_stack([Rb, Zb])).contains_points(
        np.column_stack([Rg.ravel(), Zg.ravel()])).reshape(Rg.shape)
    print(f"chord   : R[{R[0]:.2f},{R[-1]:.2f}] x Z[{Z[0]:.2f},{Z[-1]:.2f}] m, "
          f"{int(inside.sum())}/{inside.size} pts inside LCFS")

    # thermal emissivity on grid (N/CM3/SEC -> n/m3/s) + chord/torus integrals
    th_grid = profile_on_grid(x_th, th_prof, rhot, inside) * 1.0e6
    th_los, th_innerZ = chord_integral(th_grid, R, Z, args.wtor)
    th_tor = toroidal_integral(th_grid, R, Z)

    # beam-target zone reactivity. With --bt-aniso each zone's NUBEAM BTN rate is
    # weighted by the directional factor g(zone) for neutrons emitted up toward
    # the KM14 detector (a point detector_height above Zmag on the chord axis);
    # otherwise the isotropic 4*pi BTN rate is used. The 0D cross-checks further
    # down keep the isotropic bt_zone_sum unchanged.
    bt_zone_eff = bt_zone_sum
    aniso = None
    if args.bt_aniso:
        det_Z = eq.Zmag + args.detector_height
        fi["_neut_ref"] = neut
        g_map, aniso = bt_anisotropy_factors(
            fi, cdf_path, bt_keys, args.detector_R, det_Z, args)
        bt_zone_eff = sum(neut[k] * g_map[k] for k in bt_keys)
        amn, amed, amx = aniso["ang_bL"]
        if aniso.get("los"):
            li = aniso["los"]
            print(f"BT aniso: LOS direction from {args.los_file} "
                  f"({li['n_used']}/{li['n_zone']} zones in footprint; "
                  f"vs point-detector model med/max "
                  f"{li['dang_med']:.1f}/{li['dang_max']:.1f} deg)")
        else:
            print(f"BT aniso: detector (R,Z)=({aniso['det_R']:.2f},"
                  f"{aniso['det_Z']:.2f}) m (point detector), ", end="")
        print(f"max cone tilt {aniso['tilt_max']:.1f} deg; "
              f"angle(B,LOS) min/med/max {amn:.0f}/{amed:.0f}/{amx:.0f} deg "
              f"(nsamp={args.nsamp}, cone={args.cone_deg:g} deg)")

    # beam-target emissivity on grid. The flux-surface average is computed
    # unconditionally -- it feeds the rhot weight function below (a flux-surface
    # quantity) even when --bt-mode zone is used for the 2-D maps.
    def _bt_chord(zone_arr):
        xb, pb = flux_avg_profile(zone_arr, fi["x2d"], fi["bmvol"])
        if args.bt_mode == "flux":
            grid = profile_on_grid(xb, pb, rhot, inside) * 1.0e6
        else:
            grid = zone_emis_on_grid(zone_arr, fi["r2d"] / 100.0,
                                     fi["z2d"] / 100.0, Rg, Zg, inside) * 1.0e6
        los, innerZ = chord_integral(grid, R, Z, args.wtor)
        return grid, los, innerZ, xb, pb

    bt_grid, bt_los, bt_innerZ, x_bt, bt_prof = _bt_chord(bt_zone_eff)
    bt_tor = toroidal_integral(bt_grid, R, Z)
    # isotropic chord rate for the anisotropy cross-check (same geometry)
    bt_los_iso = _bt_chord(bt_zone_sum)[1] if args.bt_aniso else bt_los

    # 0D plasma totals for cross-check
    with Dataset(cdf_path, "r") as d:
        dvol = np.array(d["DVOL"][eq.tind]) * 1.0e-6   # m^3
        th_all = np.array(d[th_var][eq.tind]) * 1.0e6  # n/m3/s
    th_plasma = float(np.sum(th_all * dvol))
    bt_plasma = float(np.sum(bt_zone_sum * fi["bmvol"]))

    print("\n==================  KM14 LOS TH / BT  ====================")
    print(f"  TH chord rate   = {th_los:.4e} n/s")
    print(f"  BT chord rate   = {bt_los:.4e} n/s")
    print(f"  >> (TH/BT)_LOS  = {th_los / bt_los:.4f}    "
          f"(BT/TH = {bt_los / th_los:.4f})")
    if args.bt_aniso:
        gfac = bt_los / bt_los_iso                  # chord-integrated <g>_BT
        print(f"     [BT anisotropy ON]  isotropic (TH/BT)_LOS = "
              f"{th_los / bt_los_iso:.4f}")
        print(f"     chord <g>_BT = {gfac:.4f}  => (TH/BT)_LOS shifts "
              f"{(1.0 / gfac - 1.0) * 100:+.2f}% vs isotropic")
        for k in bt_keys:
            if k in aniso["gbar"]:
                cons = aniso["consistency"].get(k)
                cons_s = f"  [vector Eps_4pi/NUBEAM = {cons:.3f}]" if cons else ""
                print(f"       reactivity-wtd g_{k} = {aniso['gbar'][k]:.4f}{cons_s}")
    print("  ---- cross-checks (ratios; chord geometry cancels) ----")
    print(f"  (TH/BT) full-torus over chord area = {th_tor / bt_tor:.4f}")
    print(f"  (TH/BT) whole plasma (0D)          = {th_plasma / bt_plasma:.4f}")
    print(f"     core enhancement LOS vs plasma  = "
          f"{(th_los / bt_los) / (th_plasma / bt_plasma):.3f}x")
    print("==========================================================")
    print(f"  TH whole-plasma   = {th_plasma:.4e} n/s")
    btnts_parts = [(k, totals['BTNTS_' + k]) for k in bt_keys
                   if ('BTNTS_' + k) in totals]
    btnts_sum = sum(v for _, v in btnts_parts) if btnts_parts else float('nan')
    btnts_str = '+'.join(f"BTNTS_{k}" for k, _ in btnts_parts) or 'n/a'
    print(f"  BT whole-plasma   = {bt_plasma:.4e} n/s   "
          f"(0D {btnts_str} = {btnts_sum:.4e})")
    if len(btnts_parts) > 1:
        for k, v in btnts_parts:
            print(f"     BTNTS_{k:<2} = {v:.4e} n/s")

    # ----- LOS weight function f(rhot) and TH/BT vs rhot --------------------
    wp = rhot_weight_profiles(eq, rhot, inside, R, Z, x_th, th_prof,
                              x_bt, bt_prof, cdf_path)
    print("\n----------- LOS weight function & TH/BT vs rhot -----------")
    print(f"  f(rhot) = LOS shell volume / TRANSP DVOL   (max f = {wp['f'].max():.4f})")
    print(f"  Sum(TH*f*DVOL) = {wp['sum_th']:.4e} n/s   "
          f"(th_tor = {th_tor:.4e}, resid {abs(wp['sum_th']-th_tor)/th_tor:.2e})")
    print(f"  Sum(BT*f*DVOL) = {wp['sum_bt']:.4e} n/s   "
          f"(bt_tor = {bt_tor:.4e}, resid {abs(wp['sum_bt']-bt_tor)/bt_tor:.2e})")
    if args.bt_mode == "flux":
        print(f"  cumulative (TH/BT)(rhot=1) = {wp['ratio_cum'][-1]:.4f}   "
              f"(should match (TH/BT)_LOS = {th_los / bt_los:.4f})")
    else:
        # The rhot weight function always uses the flux-averaged BT, so its
        # endpoint matches the flux-mode LOS ratio, not the asymmetric zone one.
        print(f"  cumulative (TH/BT)(rhot=1) = {wp['ratio_cum'][-1]:.4f}   "
              f"(flux-averaged BT; differs from zone (TH/BT)_LOS "
              f"= {th_los / bt_los:.4f} by the poloidal asymmetry)")

    # ----- optional real-LOS (etendue C) detector weighting ---------------
    # f_det(rhot) = C_bin/DVOL is the per-shell detector coupling from the real
    # KM14 LOS cell file (purely geometric, emissivity-independent), reused from
    # los_thermal_rate. Apply it to BOTH TH and BT (BT ~ isotropic for KM14, so
    # the etendue coupling is the same) to get the real-detector cumulative
    # TH/BT, alongside the geometric-LOS and full-plasma cumulatives.
    los_det = None
    if args.los_file:
        from los_thermal_rate import (EqCDF, los_file_detector_rate,
                                       read_los_file)
        eqcdf = EqCDF(cdf_path, fi["time"] + 40.0)   # EqCDF wants JET time
        cells_los = read_los_file(args.los_file)
        fres = los_file_detector_rate(
            cells_los, eqcdf, wp["xs"], wp["th_si"] / 1.0e6, 1.0e6,
            wp["dvol_m3"])
        f_det = fres["f_det"]
        cum_th_d = np.cumsum(wp["th_si"] * f_det * wp["dvol_m3"])
        cum_bt_d = np.cumsum(wp["bt_si"] * f_det * wp["dvol_m3"])
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio_cum_det = np.where(cum_bt_d > 0, cum_th_d / cum_bt_d, np.nan)
        los_det = dict(
            f_det=f_det, ratio_cum_det=ratio_cum_det,
            rhot_crit=fres["rhot_crit"], thbt_det=float(ratio_cum_det[-1]),
            R_cells=fres["R_cells"], Z_cells=fres["Z_cells"],
            C_cells=fres["C_cells"], inside_cells=fres["inside_cells"])
        print(f"  real-LOS (C-weighted) cumulative (TH/BT)(rhot=1) = "
              f"{los_det['thbt_det']:.4f}   "
              f"(geometric-LOS {wp['ratio_cum'][-1]:.4f}, "
              f"full plasma {wp['ratio_cum_plasma'][-1]:.4f})")

    plot_label = chan_label + (" [BT aniso]" if args.bt_aniso else "")
    if args.save:
        _save_weight_profile(wp, run_id, idx, chan_label, th_var, bt_keys,
                             fi["time"], th_los / bt_los, aniso=args.bt_aniso)

    if (args.plot or args.save) and not args.no_plot:
        make_plot(R, Z, rhot, inside, th_grid, bt_grid, th_innerZ, bt_innerZ,
                  Rb, Zb, run_id, idx, plot_label, th_los, bt_los, wp,
                  save=args.save, los_det=los_det)
    return 0


def _save_weight_profile(wp, run_id, idx, chan_label, th_var, bt_keys,
                         time_s, th_bt_los, aniso=False):
    """Write the rhot weight-function / TH-BT profile to src/tmp/ as text."""
    import os
    src_dir = os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
    outdir = os.path.join(src_dir, "tmp")
    os.makedirs(outdir, exist_ok=True)
    tag = "_aniso" if aniso else ""
    out_path = os.path.join(
        outdir, f"{run_id}_KM14_LOS_THBT_weight_idx{idx}_t{time_s:.3f}s{tag}.txt")
    data = np.column_stack([
        wp["xs"], wp["f"], wp["th_si"], wp["bt_si"],
        wp["thkm14"], wp["btkm14"], wp["dvol_m3"],
        wp["shell_th"], wp["shell_bt"], wp["ratio_local"], wp["ratio_cum"],
    ])
    np.savetxt(
        out_path, data,
        header=(f"{run_id}  idx {idx}  channel {chan_label}"
                f"{'  [BT anisotropy: dEps/dOmega toward KM14]' if aniso else ''}  "
                f"TH={th_var}  BT=BTN[{'+'.join(bt_keys)}]  t = {time_s:.3f} s\n"
                f"(TH/BT)_LOS = {th_bt_los:.6f}  "
                f"[cum (TH/BT)(rhot=1) must match this]\n"
                f"  rhot          f             TH[n/m3/s]    BT[n/m3/s]    "
                f"THKM14        BTKM14        DVOL[m3]      "
                f"TH*f*DVOL[n/s] BT*f*DVOL[n/s] TH/BT_local   TH/BT_cum"),
        fmt="%.8e",
    )
    print(f"  saved weight profile -> {out_path}")
    return out_path


def make_plot(R, Z, rhot, inside, th_grid, bt_grid, th_innerZ, bt_innerZ,
              Rb, Zb, run_id, idx, chan_label, th_los, bt_los, wp, save=None,
              los_det=None):
    import matplotlib
    if save is not None:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    Rg, Zg = np.meshgrid(R, Z, indexing="ij")
    # Local TH/BT map: only where BT is a meaningful fraction of its peak. Near
    # the LCFS the BT flux profile -> 0 faster than TH, so the bare ratio spikes
    # in a few edge cells and autoscaling washes the whole map to one flat colour
    # -- mask those out and clip the colour range to robust percentiles so the
    # physical core/edge gradient is visible.
    bt_peak = np.nanmax(np.where(inside, bt_grid, np.nan))
    valid = inside & (bt_grid > 1e-3 * bt_peak) & (th_grid > 0)
    ratio_grid = np.divide(th_grid, bt_grid, out=np.full_like(th_grid, np.nan),
                           where=valid)
    rin = ratio_grid[np.isfinite(ratio_grid)]
    rlo, rhi = (float(np.percentile(rin, 2)), float(np.percentile(rin, 98))) \
        if rin.size else (0.0, 1.0)
    fig, ax = plt.subplots(2, 4, figsize=(20, 11))
    # ---- top row: (R, Z) emissivity maps + R-integrated vs Z ----
    for a, grid, ttl, cmap, lims in [
            (ax[0, 0], np.where(inside, th_grid, np.nan), "eps_TH [n/m3/s]", "inferno", (None, None)),
            (ax[0, 1], np.where(inside, bt_grid, np.nan), "eps_BT [n/m3/s]", "inferno", (None, None)),
            (ax[0, 2], ratio_grid, "TH/BT (local)", "coolwarm", (rlo, rhi))]:
        pc = a.pcolormesh(Rg, Zg, grid, shading="auto", cmap=cmap,
                          vmin=lims[0], vmax=lims[1])
        a.plot(Rb, Zb, "c-", lw=1)
        fig.colorbar(pc, ax=a)
        a.set_aspect("equal"); a.set_xlabel("R [m]"); a.set_ylabel("Z [m]")
        a.set_title(ttl)
    # Real KM14 LOS footprint overlaid on the TH/BT(local) map: shows which part
    # of the ratio field the detector actually samples (slanted, narrow band).
    if los_det is not None:
        ins = np.asarray(los_det["inside_cells"])
        Rc = np.asarray(los_det["R_cells"]); Zc = np.asarray(los_det["Z_cells"])
        sel = np.where(ins)[0]
        if sel.size > 15000:
            sel = sel[np.linspace(0, sel.size - 1, 15000).astype(int)]
        ax[0, 2].scatter(Rc[sel], Zc[sel], s=1.0, c="lime", alpha=0.12,
                         linewidths=0)
        ax[0, 2].set_title("TH/BT (local) + real LOS cells")
    ax[0, 3].plot(th_innerZ, Z, "r-", label="TH")
    ax[0, 3].plot(bt_innerZ, Z, "b-", label="BT")
    ax[0, 3].set_xlabel(r"$\int \epsilon\, dR$ [n/m2/s]"); ax[0, 3].set_ylabel("Z [m]")
    ax[0, 3].set_title(f"R-integrated vs Z\n(TH/BT)_LOS={th_los/bt_los:.3f}")
    ax[0, 3].legend(); ax[0, 3].grid(True, ls=":")

    # ---- bottom row: LOS weight function & TH/BT vs rhot ----
    xs = wp["xs"]
    a = ax[1, 0]                                   # geometric LOS weight f(rhot)
    a.plot(xs, wp["f"], "g-", lw=1.4, label="f geometric (box)")
    a.fill_between(xs, 0.0, wp["f"], color="g", alpha=0.15)
    a.set_xlim(0, 1); a.set_ylim(0, 1.05); a.grid(True, ls=":")
    a.set_xlabel("rhot"); a.set_ylabel(r"$f(\rho_t)$")
    a.set_title("LOS weight function")
    if los_det is not None:                        # real-LOS detector coupling
        fd = np.asarray(los_det["f_det"])
        fdn = fd / np.nanmax(fd) if np.nanmax(fd) > 0 else fd
        a.plot(xs, fdn, "m--", lw=1.4, label="f_det (real LOS, norm.)")
        rc = los_det.get("rhot_crit")
        if rc is not None and np.isfinite(rc):
            a.axvline(rc, color="c", ls=":", lw=1.0,
                      label=f"rhot_crit={rc:.3f}")
        a.legend(fontsize=8)

    a = ax[1, 1]                                   # per-shell LOS rate eps*f*DVOL
    a.plot(xs, wp["shell_th"], "r-", lw=1.4, label=r"TH$\cdot f\cdot$DVOL")
    a.plot(xs, wp["shell_bt"], "b-", lw=1.4, label=r"BT$\cdot f\cdot$DVOL")
    a.set_xlim(0, 1); a.grid(True, ls=":")
    a.set_xlabel("rhot"); a.set_ylabel("per-shell LOS rate [n/s]")
    a.set_title("LOS contribution per flux shell")
    a.legend(fontsize=8)

    a = ax[1, 2]                                   # local TH/BT(rhot)
    a.plot(xs, wp["ratio_local"], "m-", lw=1.4)
    a.axhline(th_los / bt_los, color="k", ls="--", lw=0.9,
              label=f"(TH/BT)_LOS={th_los/bt_los:.3f}")
    a.set_xlim(0, 1); a.grid(True, ls=":")
    a.set_xlabel("rhot"); a.set_ylabel("TH/BT (local)")
    a.set_title(r"Local ratio  TH$(\rho_t)$/BT$(\rho_t)$")
    a.legend(fontsize=8)

    a = ax[1, 3]                                   # cumulative TH/BT(rhot)
    a.plot(xs, wp["ratio_cum"], "c-", lw=1.6,
           label=f"LOS geometric (rhot=1: {wp['ratio_cum'][-1]:.3f})")
    if los_det is not None:                        # real-LOS C(etendue)-weighted
        a.plot(xs, los_det["ratio_cum_det"], "m-", lw=1.8,
               label=f"real LOS C-weighted (rhot=1: {los_det['thbt_det']:.3f})")
    a.plot(xs, wp["ratio_cum_plasma"], color="0.35", ls="--", lw=1.4,
           label=f"TRANSP full plasma (rhot=1: {wp['ratio_cum_plasma'][-1]:.3f})")
    a.set_xlim(0, 1); a.grid(True, ls=":")
    a.set_xlabel("rhot")
    a.set_ylabel(r"$\sum$TH$\,$DVOL / $\sum$BT$\,$DVOL  ($\times f$: LOS)")
    a.set_title("Cumulative TH/BT  (LOS vs full plasma)")
    a.legend(fontsize=8)

    fig.suptitle(f"{run_id} idx{idx}  KM14 LOS TH/BT  {chan_label}")
    fig.tight_layout()

    # ---- second figure: effective per-shell weights vs rhot ----
    # Left panel: each weight (f*DVOL, f_det*DVOL, DVOL) normalized to a true
    # PDF in rhot (divided by its trapezoidal integral over [0, 1] so the area
    # under the curve is 1) -- shows where in rhot each weight concentrates.
    # Middle/right panels: the same three PDFs multiplied by TH and BT
    # emissivity -- the per-shell LOS contribution under each weighting,
    # apples-to-apples since all three weights now share the same normalisation.
    dvol = wp["dvol_m3"]
    fdvol = wp["f"] * dvol
    th_si = wp["th_si"]; bt_si = wp["bt_si"]
    def _pdf(w):
        a_int = float(np.trapezoid(w, xs))
        return w / a_int if a_int > 0 else w
    dvol_pdf = _pdf(dvol)
    fdvol_pdf = _pdf(fdvol)
    fd_arr = np.asarray(los_det["f_det"]) if los_det is not None else None
    fddvol_pdf = _pdf(fd_arr * dvol) if fd_arr is not None else None

    fig2, ax2 = plt.subplots(1, 3, figsize=(18, 5))

    a = ax2[0]
    a.plot(xs, dvol_pdf, color="0.4", ls="-", lw=1.4, label="DVOL")
    a.plot(xs, fdvol_pdf, "g-", lw=1.4, label=r"$f\cdot$DVOL (geom. LOS)")
    if fddvol_pdf is not None:
        a.plot(xs, fddvol_pdf, "m-", lw=1.4,
               label=r"$f_{\rm det}\cdot$DVOL (real LOS)")
        rc = los_det.get("rhot_crit")
        if rc is not None and np.isfinite(rc):
            a.axvline(rc, color="c", ls=":", lw=1.0,
                      label=f"rhot_crit={rc:.3f}")
    a.set_xlim(0, 1); a.grid(True, ls=":")
    a.set_xlabel("rhot"); a.set_ylabel(r"PDF in rhot  ($\int_0^1 w\,d\rho_t = 1$)")
    a.set_title("Effective per-shell weight (PDF in rhot)")
    a.legend(fontsize=9)

    for a, eps, name, color_em in [(ax2[1], th_si, "TH", "r"),
                                   (ax2[2], bt_si, "BT", "b")]:
        a.plot(xs, eps * dvol_pdf, color="0.4", ls="-", lw=1.4,
               label=name + r"$\cdot$DVOL$\,/\!\int\!$DVOL$\,d\rho_t$")
        a.plot(xs, eps * fdvol_pdf, color_em + "-", lw=1.4,
               label=name + r"$\cdot f\cdot$DVOL$\,/\!\int\!f\cdot$DVOL$\,d\rho_t$")
        if fddvol_pdf is not None:
            a.plot(xs, eps * fddvol_pdf, "m-", lw=1.4,
                   label=name + r"$\cdot f_{\rm det}\cdot$DVOL$\,/\!\int\!f_{\rm det}\cdot$DVOL$\,d\rho_t$")
            rc = los_det.get("rhot_crit")
            if rc is not None and np.isfinite(rc):
                a.axvline(rc, color="c", ls=":", lw=1.0)
        a.set_xlim(0, 1); a.grid(True, ls=":")
        a.set_xlabel("rhot")
        a.set_ylabel(name + r"$\,\cdot$ weight  [n/m$^3$/s]")
        a.set_title(name + " per-shell contribution (PDF-weighted)")
        a.legend(fontsize=8)

    fig2.suptitle(f"{run_id} idx{idx}  KM14 LOS  {chan_label}  -  effective weights")
    fig2.tight_layout()

    if save:
        fig.savefig(save, dpi=130); print(f"\nSaved plot to {save}")
        import os
        base, ext = os.path.splitext(save)
        save_w = f"{base}_weights{ext or '.png'}"
        fig2.savefig(save_w, dpi=130); print(f"Saved weights plot to {save_w}")
    else:
        plt.show()


if __name__ == "__main__":
    raise SystemExit(main())
