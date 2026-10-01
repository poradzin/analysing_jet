#!/usr/bin/env python
"""
Plot the NBI birth (deposition) footprint in R,Z from a TRANSP `_birth.cdf<n>`
file, and show the companion `BDEP_D` radial deposition profile from the main
run CDF for comparison.

Background: `_birth.cdf<n>` does NOT contain a pre-binned R,Z deposition
*profile*. It contains one row per Monte-Carlo birth marker:
  bs_r_D_MCBEAM, bs_z_D_MCBEAM   R,Z at deposition [cm]
  bs_wght_D_MCBEAM               statistical weight [real particles/s]
  bs_einj_D_MCBEAM               injection energy [eV]
  bs_xksid_D_MCBEAM              pitch v_||/v at deposition
  bs_time_D_MCBEAM               deposition time [s] (narrow window - one
                                  NUBEAM birth step, not the full discharge)
  bs_ib_D_MCBEAM                 beam/source index
To get an R,Z *rate density* map (what most people picture as "the
deposition profile in R,Z") you histogram bs_r/bs_z weighted by bs_wght.

`BDEP_D` ("D BEAM DEPOSITION (TOTAL)", N/CM3/SEC) lives in the *main* run
CDF (`<run>.CDF`), not in `_fi_<idx>.cdf`, and it is a 1-D flux-surface
profile vs (TIME3, X), not an R,Z map — it's the birth markers already
binned onto radius and flux-surface-averaged.

Self-contained on the TRANSP CDFs (netCDF4 + numpy + matplotlib only; the
LCFS overlay uses the local ``equilibrium.py`` helper, no other part of
analysing_jet needed).

Usage:
    python plot_birth_deposition.py 49365 M04 --idx 1
    python plot_birth_deposition.py 49365 M04 --idx 1 --data-dir ~/jet/data/MAST --save birth.png
"""
import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset

from equilibrium import load_lcfs


def find_run_dir(pulse, run, data_dir=None):
    """Search order: --data-dir, ~/jet/data/MAST/<pulse>/<run>, ~/jet/data/<pulse>/<run>."""
    home = Path(os.environ["HOME"])
    candidates = []
    if data_dir is not None:
        candidates.append(Path(data_dir).expanduser() / str(pulse) / run)
    candidates.append(home / "jet" / "data" / "MAST" / str(pulse) / run)
    candidates.append(home / "jet" / "data" / str(pulse) / run)
    for c in candidates:
        if (c / f"{pulse}{run}_birth.cdf1").exists() or (c / f"{pulse}{run}.CDF").exists():
            return c
    raise FileNotFoundError(
        f"No run directory with {pulse}{run}_birth.cdf1 / {pulse}{run}.CDF found in: "
        + ", ".join(str(c) for c in candidates)
    )


def load_birth(run_dir, pulse, run, idx):
    path = run_dir / f"{pulse}{run}_birth.cdf{idx}"
    ds = Dataset(path, "r")
    r_cm = np.array(ds.variables["bs_r_D_MCBEAM"][:])
    z_cm = np.array(ds.variables["bs_z_D_MCBEAM"][:])
    wght = np.array(ds.variables["bs_wght_D_MCBEAM"][:])  # real particles/s per marker
    ib = np.array(ds.variables["bs_ib_D_MCBEAM"][:])
    time = np.array(ds.variables["bs_time_D_MCBEAM"][:])
    ds.close()
    return path, r_cm / 100.0, z_cm / 100.0, wght, ib, time


def load_bdep_d(run_dir, pulse, run):
    path = run_dir / f"{pulse}{run}.CDF"
    ds = Dataset(path, "r")
    time3 = np.array(ds.variables["TIME3"][:])
    if "BDEP_D" not in ds.variables:
        ds.close()
        return path, time3, None, None
    bdep = np.array(ds.variables["BDEP_D"][:])   # (TIME3, X) N/CM3/SEC
    x = np.array(ds.variables["X"][:])            # (TIME3, X) normalized rho grid
    ds.close()
    return path, time3, x, bdep


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pulse", type=int)
    ap.add_argument("run", type=str)
    ap.add_argument("--idx", type=int, default=1, help="birth file index (_birth.cdf<idx>), default 1")
    ap.add_argument("--data-dir", default=None, help="override base data dir (default: search order below)")
    ap.add_argument("--bins", type=int, default=80, help="R,Z histogram bins per axis")
    ap.add_argument("--save", default=None, help="save PNG instead of showing interactively")
    args = ap.parse_args()

    run_dir = find_run_dir(args.pulse, args.run, args.data_dir)
    print(f"Run directory: {run_dir}")

    birth_path, r, z, wght, ib, time = load_birth(run_dir, args.pulse, args.run, args.idx)
    print(f"Birth file: {birth_path}")
    print(f"  markers: {len(r)}, beams present: {np.unique(ib)}")
    print(f"  deposition time window: [{time.min():.5f}, {time.max():.5f}] s")
    print(f"  R range: [{r.min():.3f}, {r.max():.3f}] m, Z range: [{z.min():.3f}, {z.max():.3f}] m")
    print(f"  total deposited rate (sum of marker weights): {wght.sum():.4e} particles/s")

    bdep_path, time3, xrho, bdep = load_bdep_d(run_dir, args.pulse, args.run)
    have_bdep = bdep is not None
    print(f"Main CDF: {bdep_path}")
    # nearest TIME3 slice to the birth-file deposition window -- used both for
    # BDEP_D (if present) and for the LCFS overlay
    t_mid = 0.5 * (time.min() + time.max())
    i_t = int(np.abs(time3 - t_mid).argmin())
    print(f"  nearest TIME3 slice to birth window midpoint {t_mid:.5f}s: t={time3[i_t]:.5f}s (index {i_t})")
    if have_bdep:
        print(f"  BDEP_D shape (TIME3, X): {bdep.shape} — a 1-D rho profile vs time, NOT an R,Z map")
    else:
        print("  no BDEP_D in this CDF — skipping radial-profile panel")

    Rb, Zb = load_lcfs(bdep_path, i_t)

    # --- R,Z weighted histogram: this is the actual "deposition rate in R,Z" map ---
    r_edges = np.linspace(r.min(), r.max(), args.bins + 1)
    z_edges = np.linspace(z.min(), z.max(), args.bins + 1)
    hist, _, _ = np.histogram2d(r, z, bins=[r_edges, z_edges], weights=wght)
    # per-bin footprint rate [particles/s per bin]; NOT divided by cell volume,
    # so this is a deposition *rate density in the R,Z plane*, not a true n/cm^3/s
    hist = hist.T  # imshow wants (Z, R)

    fig, axes = plt.subplots(1, 2 if have_bdep else 1, figsize=(12 if have_bdep else 6, 5))
    if not have_bdep:
        axes = [axes]

    ax0 = axes[0]
    im = ax0.imshow(
        hist,
        origin="lower",
        aspect="equal",
        extent=[r_edges[0], r_edges[-1], z_edges[0], z_edges[-1]],
        cmap="plasma",
    )
    ax0.scatter(r, z, s=1, c="white", alpha=0.15, linewidths=0)
    ax0.plot(np.append(Rb, Rb[0]), np.append(Zb, Zb[0]), color="cyan", lw=1.5,
              label=f"LCFS @ t={time3[i_t]:.4f}s")
    ax0.legend(loc="lower right", fontsize="small")
    fig.colorbar(im, ax=ax0, label="deposited rate per (R,Z) bin [particles/s]")
    ax0.set_xlabel("R [m]")
    ax0.set_ylabel("Z [m]")
    ax0.set_title(
        f"{args.pulse}{args.run} birth.cdf{args.idx}\n"
        f"NBI deposition footprint, t=[{time.min():.4f},{time.max():.4f}]s "
        f"({len(r)} markers)"
    )

    if have_bdep:
        ax1 = axes[1]
        ax1.plot(xrho[i_t], bdep[i_t], "k-", lw=2)
        ax1.set_xlabel(r"$\rho$ (X, normalized flux radius)")
        ax1.set_ylabel(r"BDEP_D [N/CM$^3$/SEC]")
        ax1.set_title(
            f"BDEP_D at TIME3={time3[i_t]:.4f}s\n"
            "main .CDF, flux-surface-averaged radial profile\n(NOT an R,Z map)"
        )

    fig.tight_layout()
    if args.save:
        fig.savefig(args.save, dpi=150)
        print(f"Saved: {args.save}")
    else:
        plt.show()


if __name__ == "__main__":
    main()
