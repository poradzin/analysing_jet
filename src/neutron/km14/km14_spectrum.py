#!/usr/bin/env python3
"""km14_spectrum.py -- forward-model the KM14 diamond neutron spectrum (TH + BT)
along the KM14 line of sight, for overlay on the measured #104614 spectrum.

KM14 is a single-crystal diamond detector. It records the *deposited* energy of
14 MeV DT neutrons via the 12C(n,alpha0)9Be reaction, so the neutron energy E_n
appears at E_dep = E_n + Q with Q ~ -5.7 MeV (the 14.03 MeV DT line -> ~8.4 MeV
deposited). This script builds the thermal-DT and beam-thermal-DT *neutron energy
spectra* weighted by the KM14 line of sight (the same geometric LOS weight
f(rhot) used by los_thermal_rate.py / los_th_bt_ratio.py), maps E_n -> E_dep,
convolves with the diamond energy resolution, and overlays the result on the
experimental figure.

Components
----------
* **Thermal DT** -- each TRANSP flux shell emits a Doppler-broadened Gaussian
  line: mean E0 (~14.03 MeV), FWHM = 177*sqrt(Ti[keV]) keV (DT Brysk width).
  Shells are weighted by the LOS-seen rate THNTX_DT(rhot) * f(rhot) * DVOL.
* **Beam-thermal DT** -- per NUBEAM zone, MC-sample the fast-D distribution F
  (speed, pitch, gyrophase via B_hat) + a thermal-T Maxwellian (+rotation),
  weight by sigma*v_rel, emit an isotropic-CM neutron -> lab velocity, keep the
  neutrons heading into the KM14 detector cone, and histogram their lab energy.
  Each zone is scaled by its LOS-seen rate BTN1(zone) * f(rhot_zone) * BMVOL.
  (Reuses the bt_los_emissivity / bt_kinematics machinery.)

The combined model (TH:BT fixed by the LOS-integrated rates) is scaled so its
total peak matches the measured peak (default 145 counts), i.e. only one global
amplitude is free -- the TH/BT split is the physics prediction under test.

CLI
---
    python km14_spectrum.py 104614 M29 --idx 1 --save
    python km14_spectrum.py 104614 M29 --idx 1 --det-fwhm 0.2 --nsamp 60000
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np
from netCDF4 import Dataset
from matplotlib.path import Path as MplPath
from scipy.ndimage import gaussian_filter1d

import bt_kinematics as btk
import bt_zone_integrator as bzi
import bt_los_emissivity as ble
import los_th_bt_ratio as ltb
from los_thermal_rate import los_shell_fraction

# ----------------------------------------------------------------------
# Physics / detector constants
# ----------------------------------------------------------------------
E0_DT_MEV = 14.03          # DT neutron birth energy (stationary reactants)
DT_WIDTH_C = 177.0         # FWHM[keV] = C*sqrt(Ti[keV]); DT thermal Doppler width
EDEP_PEAK_DEFAULT = 8.40   # measured TH peak position in E_dep [MeV]

# Experimental-figure axis calibration (px -> data), from the PNG axes/ticks.
PNG_DEFAULT = ("/home/mporad/jet/data/104614/figs/"
               "104614_KM14_spectral_analysis.png")
PNG_CAL = dict(x0=129, e0=7.2, x1=725, e1=9.2,    # x-ticks 7.2 .. 9.2 MeV
               y0=584, c0=0.0, y1=32, c1=200.0)   # y-spine 0 .. 200 counts
PNG_XLIM = (7.2, 9.5)
PNG_YLIM = (0.0, 200.0)


# ----------------------------------------------------------------------
# LOS weight f(rhot) on the TRANSP X grid (shared by TH and BT)
# ----------------------------------------------------------------------

def los_weight_on_x(eq, cdf_path, fi_path, args):
    """Geometric LOS shell fraction f(rhot) on the sorted TRANSP X grid.

    Builds the KM14 chord (R, Z) grid exactly as los_th_bt_ratio.main does and
    calls the shared los_shell_fraction. Returns (xs, f, dvol_cm3)."""
    Rb, Zb = ltb.read_lcfs(fi_path)
    R = np.linspace(args.Rmin, args.Rmax, args.nR)
    Z = np.linspace(float(eq.ZG.min()), float(eq.ZG.max()), args.nZ)
    if R[0] <= eq.Rmag <= R[-1]:
        R = np.unique(np.concatenate([R, [eq.Rmag]]))
    if Z[0] <= eq.Zmag <= Z[-1]:
        Z = np.unique(np.concatenate([Z, [eq.Zmag]]))
    Rg, Zg = np.meshgrid(R, Z, indexing="ij")
    rhot = eq.rhot_pinned(Rg.ravel(), Zg.ravel(), Rb=Rb, Zb=Zb).reshape(Rg.shape)
    inside = MplPath(np.column_stack([Rb, Zb])).contains_points(
        np.column_stack([Rg.ravel(), Zg.ravel()])).reshape(Rg.shape)

    with Dataset(str(cdf_path), "r") as d:
        X = np.array(d["X"][eq.tind])
        dvol_cm3 = np.array(d["DVOL"][eq.tind])             # CM**3
    order = np.argsort(X)
    xs = X[order]
    dvol_cm3 = dvol_cm3[order]
    dvol_m3 = dvol_cm3 * 1.0e-6
    f, _ = los_shell_fraction(rhot, inside, R, Z, xs, dvol_m3, subgrid=True,
                              rmag=eq.Rmag)
    return xs, f, dvol_cm3, int(inside.sum())


# ----------------------------------------------------------------------
# Thermal-DT neutron spectrum (Doppler-broadened, LOS-weighted)
# ----------------------------------------------------------------------

def thermal_spectrum(xs, f, dvol_cm3, cdf_path, tind, En_keV, e0_mev):
    """Sum of Ti-broadened Gaussian DT lines weighted by THNTX_DT*f*DVOL.

    Returns (spectrum [n/s/keV] on En_keV, total LOS rate [n/s])."""
    with Dataset(str(cdf_path), "r") as d:
        Xraw = np.array(d["X"][tind])
        thntx = np.array(d["THNTX_DT"][tind])              # 1/cm3/s
        ti_ev = np.array(d["TI"][tind])                    # eV
    order = np.argsort(Xraw)
    thntx = thntx[order]
    ti_keV = ti_ev[order] / 1000.0
    w_shell = thntx * dvol_cm3 * f                         # n/s per shell
    mean_keV = e0_mev * 1000.0
    spec = np.zeros_like(En_keV)
    rt2pi = np.sqrt(2.0 * np.pi)
    for wi, tk in zip(w_shell, ti_keV):
        if wi <= 0.0 or tk <= 0.0:
            continue
        sig = (DT_WIDTH_C * np.sqrt(tk)) / 2.35482
        spec += wi * np.exp(-0.5 * ((En_keV - mean_keV) / sig) ** 2) / (sig * rt2pi)
    return spec, float(np.sum(w_shell[w_shell > 0]))


# ----------------------------------------------------------------------
# Beam-thermal-DT neutron spectrum (kinematic, LOS- + detector-direction)
# ----------------------------------------------------------------------

def bt_spectrum(fi, neut_dt, cdf_path, xs, f, En_keV, det_R, det_Z, args, rng):
    """Per-zone MC lab-energy spectrum of BT-DT neutrons into the KM14 cone,
    each zone scaled by its LOS-seen rate BTN1*f*BMVOL.

    Returns (spectrum [n/s/keV] on En_keV, total LOS rate [n/s])."""
    therm = bzi.read_thermal_profiles(cdf_path, fi["time"], fi["x2d"])
    bfield = ble.BField(cdf_path, fi["time"])
    n_fast = bzi.fast_ion_density(fi)
    if args.fast_norm == "bdens":
        n_fast, _ = bzi.renormalize_to_bdens(n_fast, fi["x2d"], fi["bmvol"],
                                             therm["bdens_d"])
    Rz = fi["r2d"] / 100.0
    Zz = fi["z2d"] / 100.0
    bhat = bfield.bhat(Rz, Zz)
    if args.no_rotation:
        vrot = np.zeros_like(Rz)
    else:
        with Dataset(str(cdf_path), "r") as d:
            t3 = np.array(d["TIME3"][:])
            ti = int(np.argmin(np.abs(t3 - fi["time"])))
            X = np.array(d["X"][ti]); OM = np.array(d["OMEGA"][ti])
        vrot = np.interp(fi["x2d"], X, OM) * Rz

    # per-zone LOS direction to the detector and acceptance cone
    dR = det_R - Rz
    dZ = det_Z - Zz
    nrm = np.sqrt(dR * dR + dZ * dZ)
    n_los = np.column_stack([dR / nrm, np.zeros_like(dR), dZ / nrm])
    cos_alpha = np.cos(np.radians(args.cone_deg))

    # per-zone LOS-seen rate weight: BTN1 [1/cm3/s] * f(rhot_zone) * BMVOL [cm3]
    fz = np.clip(np.interp(fi["x2d"], xs, f), 0.0, 1.0)
    zone_los = neut_dt * fz * fi["bmvol"]                  # n/s per zone

    edges = np.concatenate([En_keV - 0.5 * (En_keV[1] - En_keV[0]),
                            [En_keV[-1] + 0.5 * (En_keV[1] - En_keV[0])]])
    dE = En_keV[1] - En_keV[0]
    spec = np.zeros_like(En_keV)
    n_th = therm["nt"]; ti_z = therm["ti"]
    nz = fi["F"].shape[0]
    for z in range(nz):
        if zone_los[z] <= 0 or n_th[z] <= 0 or ti_z[z] <= 0:
            continue
        v_speed, xi = bzi.sample_fast_speeds(fi["F"][z], fi["EB"], fi["AB"],
                                             fi["E"], fi["A"], args.nsamp, rng)
        if v_speed is None:
            continue
        v_f = ble.fast_velocity_vectors(v_speed, xi, bhat[z], rng)
        v_t = ble.thermal_velocity_vectors(ti_z[z], btk.M_T_AMU, vrot[z],
                                           args.nsamp, rng)
        E_cm = btk.relative_cm_energy_keV(v_f, v_t, btk.M_D_AMU, btk.M_T_AMU)
        w = btk.sigma_dt_mb(E_cm) * btk.MB_M2 * np.linalg.norm(v_f - v_t, axis=-1)
        cos_th, phi = btk.isotropic_cm(args.nsamp, rng)
        v_n = btk.neutron_lab_velocity(v_f, v_t, "DT", cos_th, phi)
        nh = v_n / np.linalg.norm(v_n, axis=-1, keepdims=True)
        sel = (nh @ n_los[z] >= cos_alpha) & (w > 0)
        if not np.any(sel):
            continue
        En = btk.neutron_lab_energy_keV(v_n[sel])
        h, _ = np.histogram(En, bins=edges, weights=w[sel])
        hsum = h.sum()
        if hsum <= 0:
            continue
        spec += zone_los[z] * h / (hsum * dE)              # n/s/keV
    return spec, float(np.sum(zone_los[zone_los > 0]))


# ----------------------------------------------------------------------
# Overlay
# ----------------------------------------------------------------------

def _fwhm_mev(edep_mev, y):
    """FWHM of a peak y(edep) by half-max crossings (linear interp)."""
    if y.max() <= 0:
        return float("nan")
    half = 0.5 * y.max()
    above = np.where(y >= half)[0]
    if above.size < 2:
        return float("nan")
    return float(edep_mev[above[-1]] - edep_mev[above[0]])


def make_overlay(edep_mev, th, bt, tot, png_path, save, peak_target,
                 run_id, idx, t_jet_label):
    import matplotlib
    if save is not None:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg

    c = PNG_CAL
    epx = lambda x: c["e0"] + (x - c["x0"]) * (c["e1"] - c["e0"]) / (c["x1"] - c["x0"])
    cpx = lambda y: c["c0"] + (y - c["y0"]) * (c["c1"] - c["c0"]) / (c["y1"] - c["y0"])
    img = mpimg.imread(png_path)
    H, W = img.shape[:2]
    extent = [epx(0), epx(W), cpx(H), cpx(0)]

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.imshow(img, extent=extent, aspect="auto", zorder=0)
    ax.plot(edep_mev, th, color="magenta", lw=2.0, ls="--",
            label="TH model (M29 LOS)", zorder=5)
    ax.plot(edep_mev, bt, color="lime", lw=2.0, ls=":",
            label="B-th model (M29 LOS)", zorder=5)
    ax.plot(edep_mev, tot, color="black", lw=2.2,
            label="Total model", zorder=6)
    ax.set_xlim(*PNG_XLIM); ax.set_ylim(*PNG_YLIM)
    ax.set_xlabel(r"$E_{dep}$ (MeV)"); ax.set_ylabel("Counts")
    ax.set_title(f"{run_id} idx{idx} {t_jet_label}: KM14 model TH/BT over data "
                 f"(peak->{peak_target:.0f})", fontsize=10)
    ax.legend(loc="upper left", fontsize=9, framealpha=0.85)
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=140); print(f"Saved overlay -> {save}")
    else:
        plt.show()


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
    # chord geometry (same defaults as los_th_bt_ratio)
    p.add_argument("--Rmin", type=float, default=ltb.R_MIN_DEFAULT)
    p.add_argument("--Rmax", type=float, default=ltb.R_MAX_DEFAULT)
    p.add_argument("--nR", type=int, default=ltb.NR_DEFAULT)
    p.add_argument("--nZ", type=int, default=ltb.NZ_DEFAULT)
    # BT detector / MC
    p.add_argument("--detector-R", type=float, default=ltb.DET_R_DEFAULT)
    p.add_argument("--detector-height", type=float, default=ltb.DET_HEIGHT_DEFAULT)
    p.add_argument("--cone-deg", type=float, default=20.0)
    p.add_argument("--nsamp", type=int, default=60000)
    p.add_argument("--fast-norm", choices=["bdens", "ntot"], default="bdens")
    p.add_argument("--no-rotation", action="store_true")
    p.add_argument("--seed", type=int, default=12345)
    # spectrum / detector response
    p.add_argument("--e0-dt", type=float, default=E0_DT_MEV,
                   help="DT neutron birth energy [MeV] (default 14.03)")
    p.add_argument("--edep-peak", type=float, default=EDEP_PEAK_DEFAULT,
                   help="E_dep [MeV] to align the thermal peak to (sets the "
                        "E_n->E_dep shift; default 8.40)")
    p.add_argument("--det-fwhm", type=float, default=0.20,
                   help="diamond detector energy resolution FWHM [MeV] (default 0.20)")
    p.add_argument("--de-kev", type=float, default=5.0, help="energy grid step [keV]")
    p.add_argument("--peak-counts", type=float, default=145.0,
                   help="scale the model total so its peak matches this [counts]")
    p.add_argument("--png", default=PNG_DEFAULT)
    p.add_argument("--save", nargs="?", const="__auto__", default=None,
                   help="save the overlay PNG (optional path; default into figs/)")
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

    fi = bzi.read_fi_distribution(fi_path)
    neut = bzi.read_neut_rates(neut_path)
    if "DT" not in neut:
        print(f"No DT beam-target (BTN1) in {neut_path} -- KM14 is a DT diagnostic",
              file=sys.stderr)
        return 1
    eq = ltb.CdfEquilibrium(cdf_path, fi["time"])

    print(f"Run dir : {run_dir}")
    print(f"FBM idx : {idx}   t_TRANSP = {fi['time']:.4f} s   (slice {eq.tind})")

    # LOS weight f(rhot)
    xs, f, dvol_cm3, n_in = los_weight_on_x(eq, cdf_path, fi_path, args)
    print(f"LOS     : {n_in} chord pts inside LCFS,  max f = {f.max():.3f}")

    # energy grid
    En = np.arange(11.5e3, 16.5e3 + args.de_kev, args.de_kev)   # keV
    shift_mev = args.edep_peak - args.e0_dt                     # E_dep = E_n + shift

    # spectra (n/s/keV on En)
    th_spec, th_rate = thermal_spectrum(xs, f, dvol_cm3, cdf_path, eq.tind, En,
                                        args.e0_dt)
    rng = np.random.default_rng(args.seed)
    det_Z = eq.Zmag + args.detector_height
    bt_spec, bt_rate = bt_spectrum(fi, neut["DT"], cdf_path, xs, f, En,
                                   args.detector_R, det_Z, args, rng)

    # E_n -> E_dep and detector-resolution convolution
    edep = En / 1000.0 + shift_mev                              # MeV
    sig_bins = (args.det_fwhm * 1000.0 / 2.35482) / args.de_kev
    th_c = gaussian_filter1d(th_spec, sig_bins, mode="constant")
    bt_c = gaussian_filter1d(bt_spec, sig_bins, mode="constant")
    tot_c = th_c + bt_c

    # normalize total peak to the measured peak (TH:BT kept = physics)
    scale = args.peak_counts / tot_c.max()
    th_n, bt_n, tot_n = th_c * scale, bt_c * scale, tot_c * scale

    print("\n==============  KM14 LOS neutron spectrum (DT)  ==============")
    print(f"  TH LOS rate  = {th_rate:.4e} n/s")
    print(f"  BT LOS rate  = {bt_rate:.4e} n/s")
    print(f"  (TH/BT)_LOS  = {th_rate / bt_rate:.4f}   "
          f"(TH fraction = {100 * th_rate / (th_rate + bt_rate):.1f}%)")
    print(f"  E_n->E_dep shift = {shift_mev:+.3f} MeV  "
          f"(diamond 12C(n,a0); thermal peak -> {args.edep_peak:.2f} MeV)")
    print(f"  detector FWHM = {args.det_fwhm:.2f} MeV;  model total peak "
          f"scaled to {args.peak_counts:.0f} counts (x{scale:.3e})")
    print(f"  E_dep peak   TH = {edep[np.argmax(th_c)]:.2f}  "
          f"BT = {edep[np.argmax(bt_c)]:.2f}  Total = {edep[np.argmax(tot_c)]:.2f} MeV")
    print(f"  E_dep FWHM   TH = {_fwhm_mev(edep, th_c):.2f}  "
          f"BT = {_fwhm_mev(edep, bt_c):.2f}  Total = {_fwhm_mev(edep, tot_c):.2f} MeV")
    print("=============================================================")

    if args.no_plot and args.save is None:
        return 0
    save = args.save
    if save == "__auto__":
        save = os.path.join(os.path.dirname(args.png),
                            f"{run_id}_KM14_spectrum_overlay_idx{idx}.png")
    t_lbl = "t=53.15-53.5s" if idx == 1 else f"window idx{idx}"
    make_overlay(edep, th_n, bt_n, tot_n, args.png, save, args.peak_counts,
                 run_id, idx, t_lbl)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
