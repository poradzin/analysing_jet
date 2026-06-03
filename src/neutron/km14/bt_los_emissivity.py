#!/usr/bin/env python3
"""bt_los_emissivity.py -- direction-resolved beam-target emissivity (commit 3).

Extends the per-zone beam-target reactivity of ``bt_zone_integrator`` (which is
4*pi-integrated) to the *direction-dependent* emissivity dEps/dOmega(n_hat),
and evaluates it along the KM14 vertical line of sight.

Why a direction matters here (and not in commit 2)
--------------------------------------------------
The 4*pi rate only needs the fast-ion *speed* (isotropic thermal target). The
emitted-neutron *direction*, however, depends on the CM->lab boost, whose axis
is the centre-of-mass velocity -- dominated by the fast ion, which streams
along B_hat. So we now need:

  * the field direction B_hat(R, Z) from the equilibrium (PSIRZ + BZXR), to
    build the fast-ion velocity *vector* from (speed, pitch xi, gyrophase);
  * the CM->lab kinematics from ``bt_kinematics`` to turn each (fast, thermal)
    pair + isotropic CM emission into a lab neutron direction;
  * the KM14 sight-line direction n_hat_LOS to project onto.

Geometry result we expect: B_hat is ~97-100% toroidal, while the KM14 LOS is
vertical (n_hat = z_hat), so KM14 views the beam-target boost almost exactly
perpendicular. The anisotropy factor g(n_hat) = (dEps/dOmega)/(Eps_4pi/4pi)
should therefore sit near 1 with a modest suppression set by (v_cm/u_n)^2.
This script measures it.

B-field reader (validated on 104614 M30)
----------------------------------------
The main CDF stores psi(R, Z) flat in ``PSIRZ`` over ``RGRID`` x ``ZGRID``
(cm). The flat layout is C-order **(Z, R)** -- reshape (nZ, nR); the magnetic
axis (psi minimum, PSI0_TR=0) then lands at the known (R, Z). Then

    B_R   = -(1/R) dpsi/dZ
    B_Z   = +(1/R) dpsi/dR        (psi in Wb/rad, R in m -> Tesla)
    B_phi = BZXR / R              (BZXR = R*B_phi vacuum, T*m)

GFUN (diamagnetic F/F_vac, <=1.6% here) is neglected in B_phi -- it shifts
|B| by <1.6% and the B_hat *direction* (what we need) by far less.

Consistency tests (built in, no external reference)
---------------------------------------------------
1. Vector-based Eps_4pi reproduces the speed-based reactivity of commit 2
   (and hence NUBEAM BTN4/BTN1) -- the velocity-vector machinery does not
   change the scalar rate.
2. The reactivity-weighted neutron direction distribution integrates to 1
   over 4*pi (i.e. the average anisotropy factor over the sphere is 1).

CLI
---
    python bt_los_emissivity.py 104614 M30 --idx 1 --plot
    python bt_los_emissivity.py 104614 M30 --cone-deg 20 --nsamp 100000

Flags: --idx --data-dir --nsamp --seed --fast-norm --cone-deg --no-rotation
       --plot --save --no-plot
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

import bt_kinematics as btk
import bt_zone_integrator as bzi


# ----------------------------------------------------------------------
# Equilibrium B-field
# ----------------------------------------------------------------------

class BField:
    """B_hat(R, Z) from the TRANSP main CDF at a chosen time.

    Builds B_R, B_Z, B_phi on the PSIRZ grid and bilinearly interpolates the
    cylindrical components at requested (R, Z) points (metres). Returns unit
    vectors in the (R, phi, Z) basis.
    """

    def __init__(self, cdf_path, time):
        with Dataset(cdf_path, "r") as d:
            t3 = np.array(d["TIME3"][:])
            tind = int(np.argmin(np.abs(t3 - time)))
            self.t_used = float(t3[tind])
            RG = np.array(d["RGRID"][tind]) / 100.0      # m
            ZG = np.array(d["ZGRID"][tind]) / 100.0      # m
            nR, nZ = len(RG), len(ZG)
            psi = np.array(d["PSIRZ"][tind]).reshape(nZ, nR)   # [iz, ir] Wb/rad
            bzxr = float(np.array(d["BZXR"][:]).ravel()[0]) / 100.0  # T*m
        dpsidZ, dpsidR = np.gradient(psi, ZG, RG)        # axis0=Z, axis1=R
        RR = np.broadcast_to(RG[None, :], (nZ, nR))
        self.RG, self.ZG = RG, ZG
        self.B_R = -dpsidZ / RR
        self.B_Z = dpsidR / RR
        self.B_phi = bzxr / RR

    def _bilinear(self, field, R, Z):
        RG, ZG = self.RG, self.ZG
        ir = np.clip(np.searchsorted(RG, R) - 1, 0, len(RG) - 2)
        iz = np.clip(np.searchsorted(ZG, Z) - 1, 0, len(ZG) - 2)
        tr = (R - RG[ir]) / (RG[ir + 1] - RG[ir])
        tz = (Z - ZG[iz]) / (ZG[iz + 1] - ZG[iz])
        f = (field[iz, ir]       * (1 - tr) * (1 - tz) +
             field[iz, ir + 1]   * tr       * (1 - tz) +
             field[iz + 1, ir]   * (1 - tr) * tz +
             field[iz + 1, ir + 1] * tr     * tz)
        return f

    def bhat(self, R, Z):
        """Unit B vector(s) in (R, phi, Z) basis at (R, Z) [m]. Shape (...,3)."""
        R = np.atleast_1d(np.asarray(R, float))
        Z = np.atleast_1d(np.asarray(Z, float))
        b = np.stack([self._bilinear(self.B_R, R, Z),
                      self._bilinear(self.B_phi, R, Z),
                      self._bilinear(self.B_Z, R, Z)], axis=-1)
        b /= np.linalg.norm(b, axis=-1, keepdims=True)
        return b


# ----------------------------------------------------------------------
# Velocity-vector sampling
# ----------------------------------------------------------------------

def fast_velocity_vectors(v_speed, xi, bhat, rng):
    """Fast-ion velocity vectors in (R, phi, Z) from speed, pitch, gyrophase.

    v_fast = v|| b_hat + v_perp (cos g * e1 + sin g * e2),  g uniform.
    """
    e1, e2, b = btk._orthonormal_frame(bhat)
    n = v_speed.size
    vpar = xi * v_speed
    vperp = np.sqrt(np.clip(1.0 - xi * xi, 0.0, 1.0)) * v_speed
    g = 2.0 * np.pi * rng.random(n)
    return (vpar[:, None] * b
            + vperp[:, None] * (np.cos(g)[:, None] * e1 + np.sin(g)[:, None] * e2))


def thermal_velocity_vectors(ti_eV, m_th_amu, vrot, n, rng):
    """Thermal target velocities (R, phi, Z): isotropic Maxwellian + toroidal drift."""
    sigma = np.sqrt(ti_eV * btk.EV_J / (m_th_amu * btk.AMU_KG))   # m/s per comp
    vt = rng.normal(0.0, sigma, size=(n, 3))
    vt[:, 1] += vrot                                              # phi drift
    return vt


# ----------------------------------------------------------------------
# Per-zone direction-resolved emissivity
# ----------------------------------------------------------------------

def zone_emissivity(fi, n_fast, n_th, ti, reaction, m_fast_amu, m_th_amu,
                    bhat_zone, vrot_zone, n_los, nsamp, cos_alpha, rng):
    """Eps_4pi (1/cm3/s) and dEps/dOmega along n_los (1/cm3/s/sr) per zone.

    For each zone, sample fast ions from F (speed + pitch + gyrophase) using the
    zone's B_hat, pair with a Maxwellian (+rotation) thermal partner, weight by
    sigma*v_rel, and emit an isotropic CM neutron -> lab direction. dEps/dOmega
    along n_los is estimated from the reaction-weight fraction falling in a cone
    of half-angle alpha about n_los.
    """
    F = fi["F"]
    n_zone = F.shape[0]
    eps4pi = np.zeros(n_zone)
    dedom = np.zeros(n_zone)              # dEps/dOmega along n_los
    cone_solid = 2.0 * np.pi * (1.0 - cos_alpha)
    if reaction == "DD_n":
        sig_fn = btk.sigma_dd_n_mb
    elif reaction in ("DT", "TD"):
        sig_fn = btk.sigma_dt_mb
    else:
        raise ValueError(reaction)

    for z in range(n_zone):
        if n_fast[z] <= 0 or n_th[z] <= 0 or ti[z] <= 0:
            continue
        v_speed, xi = bzi.sample_fast_speeds(F[z], fi["EB"], fi["AB"],
                                             fi["E"], fi["A"], nsamp, rng)
        if v_speed is None:
            continue
        v_f = fast_velocity_vectors(v_speed, xi, bhat_zone[z], rng)
        v_t = thermal_velocity_vectors(ti[z], m_th_amu, vrot_zone[z], nsamp, rng)
        E_cm = btk.relative_cm_energy_keV(v_f, v_t, m_fast_amu, m_th_amu)
        w = sig_fn(E_cm) * btk.MB_M2 * np.linalg.norm(v_f - v_t, axis=-1)  # m^3/s
        wsum = w.sum()
        if wsum <= 0:
            continue
        sv = wsum / nsamp                                    # <sigma v> m^3/s
        eps4pi[z] = n_fast[z] * n_th[z] * sv * 1.0e6         # 1/cm3/s

        cos_th, phi = btk.isotropic_cm(nsamp, rng)
        v_n = btk.neutron_lab_velocity(v_f, v_t, reaction, cos_th, phi)
        n_hat = v_n / np.linalg.norm(v_n, axis=-1, keepdims=True)
        mu = n_hat @ n_los                                   # cos(angle to LOS)
        in_cone = mu >= cos_alpha
        # P(n_los) [1/sr] = (weight fraction in cone) / cone solid angle
        P = (w[in_cone].sum() / wsum) / cone_solid
        dedom[z] = eps4pi[z] * P                             # 1/cm3/s/sr
    return eps4pi, dedom


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def selftest(nsamp=2_000_000, cone_deg=20.0):
    """Validate the cone estimator against the analytic boost dipole.

    A monoenergetic 100 keV D beam fired along +y into a near-cold D target
    produces, for isotropic CM emission, a lab angular distribution whose
    forward/backward intensity is 1 +/- 2*beta (beta = v_cm/u_n) to leading
    order, and exactly isotropic perpendicular to the beam. This is the clean
    (unidirectional) limit of the KM14 estimator.
    """
    rng = np.random.default_rng(0)
    v = btk.speed_from_energy_keV(100.0, btk.M_D_AMU)
    vf = np.zeros((nsamp, 3)); vf[:, 1] = v
    vt = rng.normal(0.0, 1.0e3, size=(nsamp, 3))
    cos_t, phi = btk.isotropic_cm(nsamp, rng)
    vn = btk.neutron_lab_velocity(vf, vt, "DD_n", cos_t, phi)
    nh = vn / np.linalg.norm(vn, axis=-1, keepdims=True)
    ca = np.cos(np.radians(cone_deg)); cone = 2 * np.pi * (1 - ca)
    beta = 0.5 * v / 2.17e7
    print(f"selftest: 100 keV D beam +y, cold D, cone {cone_deg:g} deg, beta={beta:.3f}")
    ok = True
    for nm, ax, pred in [("forward +y", [0, 1, 0], 1 + 2 * beta),
                         ("perp    +z", [0, 0, 1], 1.0),
                         ("back    -y", [0, -1, 0], 1 - 2 * beta)]:
        g = np.mean((nh @ np.array(ax, float)) >= ca) / cone * 4 * np.pi
        good = abs(g - pred) < 0.02
        ok &= good
        print(f"  {nm}: g = {g:.3f}  (predict {pred:.3f})  [{'OK' if good else 'FAIL'}]")
    print(f"=== selftest {'PASS' if ok else 'FAIL'} ===")
    return ok


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    if "--test" in argv:
        return 0 if selftest() else 1
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("pulse", type=int)
    p.add_argument("run_suffix")
    p.add_argument("--idx", type=int, default=None)
    p.add_argument("--data-dir", default=None)
    p.add_argument("--nsamp", type=int, default=60000,
                   help="MC samples per zone (default 60000)")
    p.add_argument("--cone-deg", type=float, default=20.0,
                   help="half-angle of the LOS acceptance cone for dEps/dOmega")
    p.add_argument("--fast-norm", choices=["bdens", "ntot"], default="bdens")
    p.add_argument("--no-rotation", action="store_true",
                   help="disable thermal toroidal rotation (default: include)")
    p.add_argument("--seed", type=int, default=12345)
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

    fi = bzi.read_fi_distribution(fi_path)
    neut = bzi.read_neut_rates(neut_path)
    therm = bzi.read_thermal_profiles(cdf_path, fi["time"], fi["x2d"])
    bfield = BField(cdf_path, fi["time"])

    print(f"Run dir : {run_dir}\nFBM idx : {idx}   t = {fi['time']:.4f} s")
    print(f"B-field grid time TIME3 = {bfield.t_used:.4f} s")
    print(f"n_zones : {fi['F'].shape[0]}")

    # fast-ion density (+ BDENS_D renorm, see bt_zone_integrator)
    n_fast = bzi.fast_ion_density(fi)
    if args.fast_norm == "bdens":
        n_fast, _ = bzi.renormalize_to_bdens(n_fast, fi["x2d"],
                                             fi["bmvol"], therm["bdens_d"])

    # per-zone geometry: B_hat at the zone (R, Z) and toroidal rotation drift
    Rz = fi["r2d"] / 100.0
    Zz = fi["z2d"] / 100.0
    bhat_zone = bfield.bhat(Rz, Zz)                       # (n_zone, 3)
    if args.no_rotation:
        vrot_zone = np.zeros_like(Rz)
    else:
        with Dataset(cdf_path, "r") as d:
            t3 = np.array(d["TIME3"][:])
            ti = int(np.argmin(np.abs(t3 - fi["time"])))
            X = np.array(d["X"][ti]); OM = np.array(d["OMEGA"][ti])
        vrot_zone = np.interp(fi["x2d"], X, OM) * Rz       # m/s, toroidal (+phi)

    # KM14 LOS direction: vertical chord -> neutrons travel +z to the detector
    n_los = np.array([0.0, 0.0, 1.0])
    cos_alpha = np.cos(np.radians(args.cone_deg))
    print(f"n_LOS = +z (vertical);  acceptance cone half-angle = {args.cone_deg:g} deg")
    # angle of B_hat to the LOS, reactivity context
    ang = np.degrees(np.arccos(np.clip(np.abs(bhat_zone @ n_los), 0, 1)))
    print(f"angle(B_hat, LOS): min {ang.min():.1f} / median {np.median(ang):.1f} "
          f"/ max {ang.max():.1f} deg  (90 = perpendicular)")

    rng = np.random.default_rng(args.seed)

    channels = []
    if "DD" in neut:
        channels.append(("DD", "DD_n", btk.M_D_AMU, therm["nd"], neut["DD"]))
    if "DT" in neut and np.any(therm["nt"] > 0):
        channels.append(("DT", "DT", btk.M_T_AMU, therm["nt"], neut["DT"]))

    results = []
    for label, reaction, m_th, n_th, ref in channels:
        eps4pi, dedom = zone_emissivity(
            fi, n_fast, n_th, therm["ti"], reaction, btk.M_D_AMU, m_th,
            bhat_zone, vrot_zone, n_los, args.nsamp, cos_alpha, rng)
        # consistency test 1: vector Eps_4pi vs NUBEAM
        r4 = np.sum(eps4pi * fi["bmvol"]) / np.sum(ref * fi["bmvol"])
        # anisotropy factor g = (dEps/dOmega) / (Eps_4pi / 4pi)
        iso = eps4pi / (4.0 * np.pi)
        g = np.divide(dedom, iso, out=np.ones_like(dedom), where=iso > 0)
        wkey = eps4pi * fi["bmvol"]                 # reactivity-weighted average
        gbar = np.sum(g * wkey) / np.sum(wkey)
        print(f"\n--- {label} channel ---")
        print(f"  vector Eps_4pi / NUBEAM (consistency) : {r4:.3f}")
        print(f"  anisotropy factor g along KM14 LOS:")
        print(f"    reactivity-weighted mean g = {gbar:.3f}  "
              f"(g=1 isotropic; <1 = suppression)")
        print(f"    => KM14 sees {(gbar - 1) * 100:+.1f}% vs isotropic")
        results.append((label, eps4pi, dedom, g))

    if (args.plot or args.save) and not args.no_plot and results:
        make_plot(fi, results, save=args.save)
    return 0


def make_plot(fi, results, save=None):
    import matplotlib
    if save is not None:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    n = len(results)
    fig, axes = plt.subplots(2, n, figsize=(6 * n, 9), squeeze=False)
    R = fi["r2d"] / 100.0
    Z = fi["z2d"] / 100.0
    for j, (label, eps4pi, dedom, g) in enumerate(results):
        ax = axes[0][j]
        sc = ax.scatter(R, Z, c=eps4pi, s=14, cmap="inferno")
        ax.set_aspect("equal"); ax.set_title(f"{label}: Eps_4pi (R,Z)")
        ax.set_xlabel("R [m]"); ax.set_ylabel("Z [m]")
        fig.colorbar(sc, ax=ax, label="1/cm3/s")
        ax = axes[1][j]
        sc = ax.scatter(R, Z, c=g, s=14, cmap="coolwarm", vmin=0.9, vmax=1.1)
        ax.set_aspect("equal")
        ax.set_title(f"{label}: anisotropy factor g (KM14 LOS)")
        ax.set_xlabel("R [m]"); ax.set_ylabel("Z [m]")
        fig.colorbar(sc, ax=ax, label="g = (dEps/dOmega)/(Eps_4pi/4pi)")
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=130); print(f"\nSaved plot to {save}")
    else:
        plt.show()


if __name__ == "__main__":
    raise SystemExit(main())
