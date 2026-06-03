#!/usr/bin/env python
"""
bt_poloidal_distribution.py
---------------------------

Spatial poloidal-angle distribution of the beam-target (BT) neutron
emissivity from a TRANSP NUBEAM run, read directly from the
``<runid>_fi_<idx>.cdf`` and ``<runid>_neut_<idx>.cdf`` outputs.

Purpose
~~~~~~~
For the KM14 LOS analysis with thermal neutrons we assume the emissivity is
flux-function (THNTX(rhot) only). BT neutrons are not isotropic in
emission direction and their spatial distribution can also have poloidal
asymmetries (banana orbits, beam deposition shape, finite-orbit width).
This script answers the *spatial* part of the question:

    Is BT emissivity flux-function, or does it depend on poloidal angle?

It maps the per-zone BT emissivity (BTN4 + optional BTN1/5/7 for DT/TT/TD
campaigns) onto the NUBEAM zonal grid

    X2D  -- sqrt(toroidal flux)  ~ rhot
    TH2D -- poloidal angle  in [-pi, pi]
    BMVOL-- zone volume

contained in the ``_fi_<idx>.cdf`` companion file.

What it does NOT do
~~~~~~~~~~~~~~~~~~~
It does not produce the *emission-direction* angular distribution at each
cell (i.e. the LOS-direction-dependent emissivity needed for KM14). That
requires folding F_D_NBI(zone, pitch, E) with the local Maxwellian target,
DD/DT cross sections, and CM->lab kinematics -- which is a separate code
(NEMO / FIDASIM / GENESIS) and a separate script.

Data layout
~~~~~~~~~~~
TRANSP run directory is searched in this order:
    1. --data-dir if given,           <data-dir>/<pulse>/<run_suffix>
    2. ~/jet/data/<pulse>/<run_suffix>          (WSL local)
    3. /common/transp_shared/Data/result/JET/<pulse>/<run_suffix>  (heimdall)

Usage
-----
    python bt_poloidal_distribution.py 104614 M30
    python bt_poloidal_distribution.py 104614 M30 --idx 2
    python bt_poloidal_distribution.py 104614 M30 --data-dir /common/transp_shared/Data/result/JET --save bt_M30_idx1.png
"""

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset


DEFAULT_LOCAL_BASE = Path.home() / "jet" / "data"
HEIMDALL_BASE = Path("/common/transp_shared/Data/result/JET")


# ----------------------------------------------------------------------
# I/O
# ----------------------------------------------------------------------

def find_run_dir(pulse, run_suffix, data_dir=None):
    candidates = []
    if data_dir is not None:
        candidates.append(Path(data_dir) / str(pulse) / run_suffix)
    candidates.append(DEFAULT_LOCAL_BASE / str(pulse) / run_suffix)
    candidates.append(HEIMDALL_BASE / str(pulse) / run_suffix)
    for c in candidates:
        if c.is_dir():
            return c
    raise FileNotFoundError(
        f"No TRANSP run directory found for {pulse}{run_suffix}. "
        f"Tried: {[str(c) for c in candidates]}"
    )


def list_fbm_indices(run_dir, run_id):
    """Return sorted list of FBM time-window indices available in run_dir."""
    idxs = []
    for p in run_dir.glob(f"{run_id}_fi_*.cdf"):
        try:
            idxs.append(int(p.stem.split("_")[-1]))
        except ValueError:
            continue
    return sorted(idxs)


def read_fi_geometry(fi_path):
    """Read NUBEAM zone geometry from _fi_*.cdf (no F_D_NBI here)."""
    with Dataset(fi_path, "r") as d:
        out = {
            "x2d":   np.array(d["X2D"][:]),         # (n_zone,)
            "th2d":  np.array(d["TH2D"][:]),        # (n_zone,) rad
            "r2d":   np.array(d["R2D"][:]),         # (n_zone,) cm
            "z2d":   np.array(d["Z2D"][:]),         # (n_zone,) cm
            "bmvol": np.array(d["BMVOL"][:]),       # (n_zone,) cm^3
            "nthzsm": np.array(d["NTHZSM"][:]),     # (n_row+1,)
            "xsurf": np.array(d["XSURF"][:]),       # (n_row+1,) row boundaries
            "thsurf": np.array(d["THSURF"][:]),
            "rsurf": np.array(d["RSURF"][:]),       # (n_row+1, ntheta) cm
            "zsurf": np.array(d["ZSURF"][:]),
            "time":  float(d["TIME"][:]),
        }
    return out


def read_neut_emissivity(neut_path):
    """Read per-zone BT neutron emissivity components from _neut_*.cdf."""
    with Dataset(neut_path, "r") as d:
        keys = set(d.variables.keys())
        components = {}
        for var, label in [("BTN4", "DD"),
                           ("BTN1", "DT"),
                           ("BTN5", "TT"),
                           ("BTN7", "TD")]:
            if var in keys:
                components[label] = np.array(d[var][:])
        if not components:
            raise RuntimeError(f"No BT components found in {neut_path}")
        bt_total = sum(components.values())
        out = {
            "time":        float(d["TA"][:]),
            "components":  components,
            "bt_total":    bt_total,
            "R_neut":      np.array(d["R"][:]),
            "Z_neut":      np.array(d["Z"][:]),
        }
        if "THNTNT2d" in keys:
            out["th_total"] = np.array(d["THNTNT2d"][:])
        if "TOTNTNF2d" in keys:
            out["tot_total"] = np.array(d["TOTNTNF2d"][:])
    return out


# ----------------------------------------------------------------------
# Geometry helpers
# ----------------------------------------------------------------------

def row_slices(nthzsm):
    """Slices into the (n_zone,) zonal arrays, one per x-row."""
    return [slice(int(nthzsm[i]), int(nthzsm[i + 1]))
            for i in range(len(nthzsm) - 1)]


def x_row_centers(x2d, rows):
    """Centre of each x-row (NUBEAM x = sqrt(toroidal flux)).

    Each row has a constant X2D value; XSURF is a finer (typically 2x)
    flux-surface grid used for RSURF/ZSURF and does not match the zonal
    row partition. Take the row centre directly from X2D.
    """
    return np.array([float(np.mean(x2d[sl])) for sl in rows])


def magnetic_axis(rsurf, zsurf):
    """Return (Rmag, Zmag) in metres from RSURF/ZSURF (cm) at xsurf=0."""
    rmag = float(np.mean(rsurf[0, :])) / 100.0
    zmag = float(np.mean(zsurf[0, :])) / 100.0
    return rmag, zmag


# ----------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------

def print_summary(fi, neut, run_id, idx):
    bt = neut["bt_total"]
    bmvol = fi["bmvol"]
    rate_total = float(np.sum(bt * bmvol))
    rmag, zmag = magnetic_axis(fi["rsurf"], fi["zsurf"])
    rows = row_slices(fi["nthzsm"])
    xs = x_row_centers(fi["x2d"], rows)
    print(f"=== {run_id}  FBM idx={idx}  BT neutron summary ===")
    print(f"  TRANSP time           : {neut['time']:.4f} s "
          f"(_fi {fi['time']:.4f} s)")
    print(f"  Components present    : {list(neut['components'].keys())}")
    print(f"  Total BT rate         : {rate_total:.3e} 1/s   "
          f"(sum_zone BTN_total * BMVOL)")
    for label, comp in neut["components"].items():
        r = float(np.sum(comp * bmvol))
        print(f"      {label:3s}              : {r:.3e} 1/s "
              f"({100*r/rate_total:5.1f}%)")
    if "th_total" in neut:
        th_rate = float(np.sum(neut["th_total"] * bmvol))
        print(f"  Cross-check  TH rate  : {th_rate:.3e} 1/s")
    if "tot_total" in neut:
        tot_rate = float(np.sum(neut["tot_total"] * bmvol))
        print(f"  Cross-check  TOT rate : {tot_rate:.3e} 1/s")
    print(f"  Magnetic axis (R,Z)   : ({rmag:.3f}, {zmag:.3f}) m")
    print(f"  Zonal grid            : {len(rows)} x-rows, "
          f"{len(bt)} total zones; xrow centres {xs}")


# ----------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------

def plot_panels(fi, neut, run_id, idx, save=None):
    bt = neut["bt_total"]                         # 1/cm^3/s
    bmvol = fi["bmvol"]                           # cm^3
    rate = bt * bmvol                             # 1/s per zone
    rmag, zmag = magnetic_axis(fi["rsurf"], fi["zsurf"])
    R_m = fi["r2d"] / 100.0
    Z_m = fi["z2d"] / 100.0
    rows = row_slices(fi["nthzsm"])
    xs = x_row_centers(fi["x2d"], rows)
    n_rows = len(rows)
    cmap = plt.cm.viridis

    # LCFS from outermost flux surface in _fi
    r_lcfs = fi["rsurf"][-1, :] / 100.0
    z_lcfs = fi["zsurf"][-1, :] / 100.0

    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    fig.suptitle(
        f"{run_id}  FBM idx={idx}  t = {neut['time']:.3f} s   "
        f"BT neutron emissivity (poloidal distribution)"
    )

    # --- A: scatter ε_BT(R, Z) on the zone grid --------------------------
    ax = axes[0, 0]
    sc = ax.scatter(R_m, Z_m, c=bt, s=30, cmap="plasma", edgecolors="none")
    ax.plot(r_lcfs, z_lcfs, "k-", lw=1, alpha=0.6)
    ax.plot(rmag, zmag, "wx", ms=8, mew=1.5)
    ax.set_aspect("equal")
    ax.set_xlabel("R [m]")
    ax.set_ylabel("Z [m]")
    ax.set_title(r"$\varepsilon_{BT}(R, Z)$ on NUBEAM zones")
    fig.colorbar(sc, ax=ax, label="1/cm³/s")

    # --- B: ε_BT(θ_pol) line plot, one curve per x-row -------------------
    ax = axes[0, 1]
    for i, sl in enumerate(rows):
        th = fi["th2d"][sl]
        em = bt[sl]
        order = np.argsort(th)
        ax.plot(th[order], em[order],
                color=cmap(i / max(n_rows - 1, 1)),
                marker="o", ms=3, lw=1,
                label=f"x={xs[i]:.2f}")
    ax.set_xlabel(r"$\theta_{pol}$ [rad]   (0 = LFS midplane, $\pm\pi$ = HFS)")
    ax.set_ylabel(r"$\varepsilon_{BT}$ [1/cm³/s]")
    ax.set_title(r"$\varepsilon_{BT}(\theta_{pol})$  per x-row")
    ax.set_xlim(-np.pi, np.pi)
    ax.legend(fontsize=8, ncol=2, loc="best")
    ax.grid(alpha=0.3)

    # --- C: Up-down asymmetry vs θ at fixed x ---------------------------
    ax = axes[0, 2]
    for i, sl in enumerate(rows):
        th = fi["th2d"][sl]
        em = bt[sl]
        order = np.argsort(th)
        th_s, em_s = th[order], em[order]
        em_mirror = np.interp(-th_s, th_s, em_s)
        denom = em_s + em_mirror
        asym = np.where(denom > 0, (em_s - em_mirror) / denom, 0.0)
        mask = th_s >= 0
        ax.plot(th_s[mask], asym[mask],
                color=cmap(i / max(n_rows - 1, 1)),
                marker="o", ms=3, lw=1,
                label=f"x={xs[i]:.2f}")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel(r"$\theta_{pol}$ [rad, top half]")
    ax.set_ylabel(r"$[\varepsilon(\theta)-\varepsilon(-\theta)]\,/\,$sum")
    ax.set_title("Up-down asymmetry")
    ax.set_xlim(0, np.pi)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2)

    # --- D: 2D scatter ε_BT(x, θ_pol) ------------------------------------
    ax = axes[1, 0]
    sc = ax.scatter(fi["th2d"], fi["x2d"], c=bt, s=40, cmap="plasma",
                    edgecolors="none")
    ax.set_xlabel(r"$\theta_{pol}$ [rad]")
    ax.set_ylabel("x = sqrt(tor flux)")
    ax.set_xlim(-np.pi, np.pi)
    ax.set_ylim(0, 1)
    ax.set_title(r"$\varepsilon_{BT}(x, \theta_{pol})$")
    fig.colorbar(sc, ax=ax, label="1/cm³/s")

    # --- E: vol-avg emissivity and rate per x-row -----------------------
    ax = axes[1, 1]
    emis_vol_avg = np.array([np.sum(rate[sl]) / np.sum(bmvol[sl])
                             for sl in rows])
    rate_per_row = np.array([np.sum(rate[sl]) for sl in rows])
    ax.plot(xs, emis_vol_avg, "ko-", label=r"$\langle\varepsilon\rangle_\theta$")
    ax.set_xlabel("x = sqrt(tor flux)")
    ax.set_ylabel(r"$\langle\varepsilon\rangle$ [1/cm³/s]", color="k")
    ax.tick_params(axis="y", labelcolor="k")
    ax2 = ax.twinx()
    ax2.plot(xs, rate_per_row, "rs--", label="rate / x-row")
    ax2.set_ylabel("rate [1/s] per x-row", color="r")
    ax2.tick_params(axis="y", labelcolor="r")
    ax.set_title("Poloidally averaged profile")
    ax.grid(alpha=0.3)

    # --- F: BT component breakdown vs x-row -----------------------------
    ax = axes[1, 2]
    for label, comp in neut["components"].items():
        rate_comp = np.array([np.sum(comp[sl] * bmvol[sl]) for sl in rows])
        ax.plot(xs, rate_comp, "o-", label=label)
    ax.set_xlabel("x = sqrt(tor flux)")
    ax.set_ylabel("rate [1/s] per x-row")
    ax.set_yscale("log")
    ax.set_title("BT components")
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3, which="both")

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    if save:
        fig.savefig(save, dpi=150, bbox_inches="tight")
        print(f"Saved figure -> {save}")
    else:
        plt.show()


# ----------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("pulse", type=int, help="JET pulse number, e.g. 104614")
    p.add_argument("run_suffix", help="TRANSP run suffix, e.g. M30")
    p.add_argument("--idx", type=int, default=None,
                   help="FBM time-window index (1, 2, ...). "
                        "Default = first available.")
    p.add_argument("--data-dir", default=None,
                   help="Override base data directory.")
    p.add_argument("--no-plot", action="store_true",
                   help="Print summary only, no plot.")
    p.add_argument("--save", default=None,
                   help="Save figure to PNG/PDF instead of showing.")
    args = p.parse_args()

    run_dir = find_run_dir(args.pulse, args.run_suffix, args.data_dir)
    run_id = f"{args.pulse}{args.run_suffix}"
    available = list_fbm_indices(run_dir, run_id)
    if not available:
        print(f"No _fi_*.cdf in {run_dir}", file=sys.stderr)
        sys.exit(1)
    idx = args.idx if args.idx is not None else available[0]
    if idx not in available:
        print(f"FBM idx={idx} not in {available}", file=sys.stderr)
        sys.exit(1)

    fi_path   = run_dir / f"{run_id}_fi_{idx}.cdf"
    neut_path = run_dir / f"{run_id}_neut_{idx}.cdf"
    print(f"Reading {fi_path}")
    fi = read_fi_geometry(fi_path)
    print(f"Reading {neut_path}")
    neut = read_neut_emissivity(neut_path)

    if abs(fi["time"] - neut["time"]) > 1e-3:
        print(f"WARNING: _fi TIME = {fi['time']:.4f} s  "
              f"!= _neut TA = {neut['time']:.4f} s", file=sys.stderr)

    print_summary(fi, neut, run_id, idx)
    if not args.no_plot:
        plot_panels(fi, neut, run_id, idx, save=args.save)


if __name__ == "__main__":
    main()
