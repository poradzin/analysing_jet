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

CLI
---
    python los_th_bt_ratio.py 104614 M30 --idx 1 --plot          # total (default)
    python los_th_bt_ratio.py 104614 M30 --channel dd --bt-mode zone

Flags: --idx --data-dir --channel {total,dd,dt} --bt-mode {flux,zone}
       --Rmin --Rmax --wtor --nR --nZ --plot --save --no-plot
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

    def rhot_pinned(self, R, Z):
        """rhot at (R, Z) [m] via cubic griddata of psi_n with the magnetic axis
        pinned at psi_n = 0.

        Bilinear interpolation of the coarse PSIRZ grid (dR ~ 2 cm) cannot reach
        the true axis minimum (psi is paraboloidal with its vertex between
        nodes), flooring psi_n at ~3e-4, i.e. rhot ~ 0.02 -- so the innermost
        f(rhot) bins get no LOS volume. Injecting the axis as a scattered node
        removes the floor. Mirrors ``los_thermal_rate.EqCDF.rhot_on_grid``.
        """
        from scipy.interpolate import griddata
        Rg_n, Zg_n = np.meshgrid(self.RG, self.ZG)  # (nZ, nR), matches psin[iz,ir]
        pts_R = np.concatenate([Rg_n.ravel(), [self.Rmag]])
        pts_Z = np.concatenate([Zg_n.ravel(), [self.Zmag]])
        pts_psin = np.concatenate([self.psin.ravel(), [0.0]])
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

    with np.errstate(invalid="ignore", divide="ignore"):
        ratio_local = np.where(bt_si > 0, th_si / bt_si, np.nan)
        ratio_cum = np.where(cum_bt > 0, cum_th / cum_bt, np.nan)

    return dict(xs=xs, f=f, th_si=th_si, bt_si=bt_si,
                thkm14=thkm14, btkm14=btkm14, dvol_m3=dvol_m3,
                shell_th=shell_th, shell_bt=shell_bt,
                cum_th=cum_th, cum_bt=cum_bt,
                ratio_local=ratio_local, ratio_cum=ratio_cum,
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
    # rhot_pinned (griddata + pinned axis) instead of plain bilinear eq.rhot,
    # so f(rhot) reaches ~1 on axis (see rhot_pinned docstring).
    rhot = eq.rhot_pinned(Rg.ravel(), Zg.ravel()).reshape(Rg.shape)
    inside = MplPath(np.column_stack([Rb, Zb])).contains_points(
        np.column_stack([Rg.ravel(), Zg.ravel()])).reshape(Rg.shape)
    print(f"chord   : R[{R[0]:.2f},{R[-1]:.2f}] x Z[{Z[0]:.2f},{Z[-1]:.2f}] m, "
          f"{int(inside.sum())}/{inside.size} pts inside LCFS")

    # thermal emissivity on grid (N/CM3/SEC -> n/m3/s)
    th_grid = profile_on_grid(x_th, th_prof, rhot, inside) * 1.0e6

    # beam-target emissivity on grid. The flux-surface average is computed
    # unconditionally -- it feeds the rhot weight function below (a flux-surface
    # quantity) even when --bt-mode zone is used for the 2-D maps.
    x_bt, bt_prof = flux_avg_profile(bt_zone_sum, fi["x2d"], fi["bmvol"])
    if args.bt_mode == "flux":
        bt_grid = profile_on_grid(x_bt, bt_prof, rhot, inside) * 1.0e6
    else:
        bt_grid = zone_emis_on_grid(bt_zone_sum, fi["r2d"] / 100.0,
                                    fi["z2d"] / 100.0, Rg, Zg, inside) * 1.0e6

    # chord and full-torus integrals
    th_los, th_innerZ = chord_integral(th_grid, R, Z, args.wtor)
    bt_los, bt_innerZ = chord_integral(bt_grid, R, Z, args.wtor)
    th_tor = toroidal_integral(th_grid, R, Z)
    bt_tor = toroidal_integral(bt_grid, R, Z)

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

    if args.save:
        _save_weight_profile(wp, run_id, idx, chan_label, th_var, bt_keys,
                             fi["time"], th_los / bt_los)

    if (args.plot or args.save) and not args.no_plot:
        make_plot(R, Z, rhot, inside, th_grid, bt_grid, th_innerZ, bt_innerZ,
                  Rb, Zb, run_id, idx, chan_label, th_los, bt_los, wp,
                  save=args.save)
    return 0


def _save_weight_profile(wp, run_id, idx, chan_label, th_var, bt_keys,
                         time_s, th_bt_los):
    """Write the rhot weight-function / TH-BT profile to src/tmp/ as text."""
    import os
    src_dir = os.path.abspath(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
    outdir = os.path.join(src_dir, "tmp")
    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(
        outdir, f"{run_id}_KM14_LOS_THBT_weight_idx{idx}_t{time_s:.3f}s.txt")
    data = np.column_stack([
        wp["xs"], wp["f"], wp["th_si"], wp["bt_si"],
        wp["thkm14"], wp["btkm14"], wp["dvol_m3"],
        wp["shell_th"], wp["shell_bt"], wp["ratio_local"], wp["ratio_cum"],
    ])
    np.savetxt(
        out_path, data,
        header=(f"{run_id}  idx {idx}  channel {chan_label}  "
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
              Rb, Zb, run_id, idx, chan_label, th_los, bt_los, wp, save=None):
    import matplotlib
    if save is not None:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    Rg, Zg = np.meshgrid(R, Z, indexing="ij")
    ratio_grid = np.divide(th_grid, bt_grid, out=np.full_like(th_grid, np.nan),
                           where=bt_grid > 0)
    fig, ax = plt.subplots(2, 4, figsize=(20, 11))
    # ---- top row: (R, Z) emissivity maps + R-integrated vs Z ----
    for a, grid, ttl, cmap in [
            (ax[0, 0], np.where(inside, th_grid, np.nan), "eps_TH [n/m3/s]", "inferno"),
            (ax[0, 1], np.where(inside, bt_grid, np.nan), "eps_BT [n/m3/s]", "inferno"),
            (ax[0, 2], ratio_grid, "TH/BT (local)", "coolwarm")]:
        pc = a.pcolormesh(Rg, Zg, grid, shading="auto", cmap=cmap)
        a.plot(Rb, Zb, "c-", lw=1)
        fig.colorbar(pc, ax=a)
        a.set_aspect("equal"); a.set_xlabel("R [m]"); a.set_ylabel("Z [m]")
        a.set_title(ttl)
    ax[0, 3].plot(th_innerZ, Z, "r-", label="TH")
    ax[0, 3].plot(bt_innerZ, Z, "b-", label="BT")
    ax[0, 3].set_xlabel(r"$\int \epsilon\, dR$ [n/m2/s]"); ax[0, 3].set_ylabel("Z [m]")
    ax[0, 3].set_title(f"R-integrated vs Z\n(TH/BT)_LOS={th_los/bt_los:.3f}")
    ax[0, 3].legend(); ax[0, 3].grid(True, ls=":")

    # ---- bottom row: LOS weight function & TH/BT vs rhot ----
    xs = wp["xs"]
    a = ax[1, 0]                                   # geometric LOS weight f(rhot)
    a.plot(xs, wp["f"], "g-", lw=1.4)
    a.fill_between(xs, 0.0, wp["f"], color="g", alpha=0.15)
    a.set_xlim(0, 1); a.set_ylim(0, 1.05); a.grid(True, ls=":")
    a.set_xlabel("rhot"); a.set_ylabel(r"$f(\rho_t)$")
    a.set_title("LOS weight function (geometric)")

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
    a.plot(xs, wp["ratio_cum"], "c-", lw=1.6)
    a.axhline(th_los / bt_los, color="k", ls="--", lw=0.9,
              label=f"(TH/BT)_LOS={th_los/bt_los:.3f}")
    a.set_xlim(0, 1); a.grid(True, ls=":")
    a.set_xlabel("rhot")
    a.set_ylabel(r"$\sum$TH$f$DVOL / $\sum$BT$f$DVOL")
    a.set_title(f"Cumulative TH/BT  (rhot=1: {wp['ratio_cum'][-1]:.3f})")
    a.legend(fontsize=8)

    fig.suptitle(f"{run_id} idx{idx}  KM14 LOS TH/BT  {chan_label}")
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=130); print(f"\nSaved plot to {save}")
    else:
        plt.show()


if __name__ == "__main__":
    raise SystemExit(main())
