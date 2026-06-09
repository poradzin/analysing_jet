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
        self.psin = (psi - psi0) / (psibnd - psi0)           # normalized psi(R,Z)
        # rhot(psi_n): psi_n on the XB grid is monotonic 0->1; XB is rhot
        psin_b = (PLFLX - PLFLX[0]) / (PLFLX[-1] - PLFLX[0])
        self._psin_b = psin_b
        self._rhot_b = XB

    def rhot(self, R, Z):
        """Bilinear psi_n at (R, Z) [m] -> rhot (sqrt tor flux). Clipped [0,1]."""
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

    # chord grid
    R = np.linspace(args.Rmin, args.Rmax, args.nR)
    Z = np.linspace(float(eq.ZG.min()), float(eq.ZG.max()), args.nZ)
    Rg, Zg = np.meshgrid(R, Z, indexing="ij")
    rhot = eq.rhot(Rg.ravel(), Zg.ravel()).reshape(Rg.shape)
    inside = MplPath(np.column_stack([Rb, Zb])).contains_points(
        np.column_stack([Rg.ravel(), Zg.ravel()])).reshape(Rg.shape)
    print(f"chord   : R[{R[0]:.2f},{R[-1]:.2f}] x Z[{Z[0]:.2f},{Z[-1]:.2f}] m, "
          f"{int(inside.sum())}/{inside.size} pts inside LCFS")

    # thermal emissivity on grid (N/CM3/SEC -> n/m3/s)
    th_grid = profile_on_grid(x_th, th_prof, rhot, inside) * 1.0e6

    # beam-target emissivity on grid
    if args.bt_mode == "flux":
        x_bt, bt_prof = flux_avg_profile(bt_zone_sum, fi["x2d"], fi["bmvol"])
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

    if (args.plot or args.save) and not args.no_plot:
        make_plot(R, Z, rhot, inside, th_grid, bt_grid, th_innerZ, bt_innerZ,
                  Rb, Zb, run_id, idx, chan_label, th_los, bt_los, save=args.save)
    return 0


def make_plot(R, Z, rhot, inside, th_grid, bt_grid, th_innerZ, bt_innerZ,
              Rb, Zb, run_id, idx, chan_label, th_los, bt_los, save=None):
    import matplotlib
    if save is not None:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    Rg, Zg = np.meshgrid(R, Z, indexing="ij")
    ratio_grid = np.divide(th_grid, bt_grid, out=np.full_like(th_grid, np.nan),
                           where=bt_grid > 0)
    fig, ax = plt.subplots(1, 4, figsize=(20, 6))
    for a, grid, ttl, cmap in [
            (ax[0], np.where(inside, th_grid, np.nan), "eps_TH [n/m3/s]", "inferno"),
            (ax[1], np.where(inside, bt_grid, np.nan), "eps_BT [n/m3/s]", "inferno"),
            (ax[2], ratio_grid, "TH/BT (local)", "coolwarm")]:
        pc = a.pcolormesh(Rg, Zg, grid, shading="auto", cmap=cmap)
        a.plot(Rb, Zb, "c-", lw=1)
        fig.colorbar(pc, ax=a)
        a.set_aspect("equal"); a.set_xlabel("R [m]"); a.set_ylabel("Z [m]")
        a.set_title(ttl)
    ax[3].plot(th_innerZ, Z, "r-", label="TH")
    ax[3].plot(bt_innerZ, Z, "b-", label="BT")
    ax[3].set_xlabel(r"$\int \epsilon\, dR$ [n/m2/s]"); ax[3].set_ylabel("Z [m]")
    ax[3].set_title(f"R-integrated vs Z\n(TH/BT)_LOS={th_los/bt_los:.3f}")
    ax[3].legend(); ax[3].grid(True, ls=":")
    fig.suptitle(f"{run_id} idx{idx}  KM14 LOS TH/BT  {chan_label}")
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=130); print(f"\nSaved plot to {save}")
    else:
        plt.show()


if __name__ == "__main__":
    raise SystemExit(main())
