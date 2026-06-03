#!/usr/bin/env python3
"""bt_zone_integrator.py -- per-zone beam-target reactivity (commit 2 of step 2).

Computes the 4*pi-integrated beam-target neutron reactivity for every NUBEAM
zone by folding the fast-ion distribution F_D_NBI against a thermal Maxwellian
target, using the Bosch-Hale cross sections from ``bt_kinematics``. The result
is validated against NUBEAM's own per-zone rates (BTN4 for D+D, BTN1 for D+T)
read from the ``_neut`` file -- this is the acceptance test for the whole
kinematics chain before the LOS-direction projection is added in commit 3.

Key physics simplification
--------------------------
For the *4*pi-integrated* reactivity the thermal target is an isotropic
Maxwellian, so <sigma * v_rel> depends only on the fast-ion *speed* -- not on
the pitch xi = v||/v, the gyrophase, or the magnetic-field direction. The
per-zone rate is therefore

    eps_BT(zone) = n_fast(zone) * n_th(zone) * <sigma(E_cm) * v_rel>      [1/cm3/s]

and no equilibrium / B-field reader is needed here. The PSIRZ -> B_hat
machinery only enters commit 3, where the *direction* of the emitted neutron
(and hence the LOS projection) finally depends on the fast-ion velocity vector.

Normalization (verified exact against NTOT_D_NBI)
-------------------------------------------------
F_D_NBI has units #/cm3/eV/d(Omega/4pi). Integrating over gyrophase
(d(Omega/4pi) -> dxi/2 for a gyrotropic distribution) gives the fast-ion
density

    n_fast(zone) = sum_{E,xi} F * dE * dxi * 0.5                         [1/cm3]

which reproduces sum(n_fast * BMVOL) = NTOT_D_NBI to machine precision.

The BDENS_D renormalization (default)
-------------------------------------
The binned distribution F_D_NBI integrates to NTOT_D_NBI = 8.18e19 ions
(104614 M30, idx 1), but the full beam-ion density profile BDENS_D integrates
to 1.01e20 -- about 24% more. NUBEAM's BTN4/BTN1 are computed from the full
beam density, so using the raw F normalization underestimates the per-zone
rate by a flat ~0.81 across radius (verified identical for the DD and DT
channels, with the kinematics validated independently). The *shape* of F
(poloidal asymmetry, energy distribution) is correct -- only its overall
normalization undercounts.

By default we therefore rescale n_fast(zone) per flux-surface row so that its
BMVOL-weighted flux average matches BDENS_D(x). This preserves F's validated
poloidal/energy structure while fixing the radial magnitude, and brings both
channels to ~0.98 of NUBEAM. Pass --fast-norm ntot to keep the raw F
normalization for comparison.

CLI
---
    python bt_zone_integrator.py 104614 M30 --idx 1 --plot
    python bt_zone_integrator.py 104614 M30 --nsamp 50000 --seed 0
    python bt_zone_integrator.py 104614 M30 --data-dir /common/transp_shared/Data/result/JET

Flags: --idx --data-dir --nsamp --seed --plot --save --no-plot
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset

import bt_kinematics as btk

DEFAULT_LOCAL_BASE = Path.home() / "jet" / "data"
HEIMDALL_BASE = Path("/common/transp_shared/Data/result/JET")


# ----------------------------------------------------------------------
# I/O  (mirrors bt_poloidal_distribution.py conventions)
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
    idxs = []
    for p in run_dir.glob(f"{run_id}_fi_*.cdf"):
        try:
            idxs.append(int(p.stem.split("_")[-1]))
        except ValueError:
            continue
    return sorted(idxs)


def read_fi_distribution(fi_path):
    """Read the fast-ion distribution and zone geometry from _fi_*.cdf."""
    with Dataset(fi_path, "r") as d:
        out = {
            "F":     np.array(d["F_D_NBI"][:]),   # (n_zone, n_xi, n_E) #/cm3/eV/d(Om/4pi)
            "E":     np.array(d["E_D_NBI"][:]),   # (n_E,)  eV  centres
            "EB":    np.array(d["EB_D_NBI"][:]),  # (n_E+1,) eV bdys
            "A":     np.array(d["A_D_NBI"][:]),   # (n_xi,) pitch centres v||/v
            "AB":    np.array(d["AB_D_NBI"][:]),  # (n_xi+1,) pitch bdys
            "x2d":   np.array(d["X2D"][:]),       # (n_zone,) zone sqrt(tor flux)
            "th2d":  np.array(d["TH2D"][:]),      # (n_zone,) rad
            "r2d":   np.array(d["R2D"][:]),       # (n_zone,) cm
            "z2d":   np.array(d["Z2D"][:]),       # (n_zone,) cm
            "bmvol": np.array(d["BMVOL"][:]),     # (n_zone,) cm^3
            "ntot":  float(d["NTOT_D_NBI"][:]),
            "time":  float(d["TIME"][:]),
        }
    return out


def read_neut_rates(neut_path):
    """Read per-zone beam-target neutron rates (1/cm3/s) from _neut_*.cdf."""
    with Dataset(neut_path, "r") as d:
        keys = set(d.variables.keys())
        out = {"time": float(d["TA"][:])}
        # BTN4 = D(d,n)3He, BTN1 = T(d,n)4He (fast D on thermal T)
        for var, label in [("BTN4", "DD"), ("BTN1", "DT"),
                           ("BTN5", "TT"), ("BTN7", "TD")]:
            if var in keys:
                out[label] = np.array(d[var][:])
    return out


def read_thermal_profiles(cdf_path, time, x2d):
    """Interpolate thermal n_D, n_T, T_i onto the zone x grid at ``time``.

    Profiles live on the TRANSP X grid (sqrt toroidal flux) in the main CDF:
        ND  [N/CM3], NT [N/CM3], TI [eV]
    The zone label X2D is the same flux coordinate, so a 1-D interpolation in
    x is appropriate (the poloidal variation of a thermal flux function is
    zero by construction).
    """
    with Dataset(cdf_path, "r") as d:
        t3 = np.array(d["TIME3"][:])
        tind = int(np.argmin(np.abs(t3 - time)))
        X = np.array(d["X"][tind])
        nd = np.array(d["ND"][tind])
        ti = np.array(d["TI"][tind])
        nt = np.array(d["NT"][tind]) if "NT" in d.variables else np.zeros_like(nd)
        bdens_d = np.array(d["BDENS_D"][tind])
        t_used = float(t3[tind])
    return {
        "tind":    tind,
        "t_used":  t_used,
        "nd":      np.interp(x2d, X, nd),       # (n_zone,) 1/cm3
        "nt":      np.interp(x2d, X, nt),       # (n_zone,) 1/cm3
        "ti":      np.interp(x2d, X, ti),       # (n_zone,) eV
        "bdens_d": np.interp(x2d, X, bdens_d),  # (n_zone,) 1/cm3 full beam-D density
    }


# ----------------------------------------------------------------------
# Physics
# ----------------------------------------------------------------------

def fast_ion_density(fi):
    """n_fast(zone) = sum_{E,xi} F * dE * dxi * 0.5   [1/cm3].

    Exact (non-MC); reproduces NTOT_D_NBI when volume-weighted.
    """
    dE = np.diff(fi["EB"])                       # (n_E,) eV
    dA = np.diff(fi["AB"])                        # (n_xi,)
    return np.einsum("zae,e,a->z", fi["F"], dE, dA) * 0.5


def renormalize_to_bdens(n_fast, x2d, bmvol, bdens_zone):
    """Rescale n_fast per flux-surface row so its flux average matches BDENS_D.

    F_D_NBI's normalization (NTOT_D_NBI) undercounts the full beam-ion density
    BDENS_D that NUBEAM uses for the beam-target rates. Scaling each x-row by
    BDENS_D(x) / <n_fast>_flux(x) preserves the poloidal asymmetry encoded in F
    while fixing the radial magnitude. Returns the scaled density and the
    per-row scale factors (for reporting).
    """
    scaled = n_fast.astype(float).copy()
    rows = {}
    for xv in np.unique(x2d):
        m = x2d == xv
        favg = np.sum(n_fast[m] * bmvol[m]) / np.sum(bmvol[m])
        bavg = np.sum(bdens_zone[m] * bmvol[m]) / np.sum(bmvol[m])
        if favg > 0:
            f = bavg / favg
            scaled[m] *= f
            rows[float(xv)] = f
    return scaled, rows


def sample_fast_speeds(F_zone, EB, AB, E, A, nsamp, rng):
    """Draw ``nsamp`` fast-ion speeds (m/s) from one zone's F(xi, E).

    Probability of bin (a, e) is proportional to F * dxi * dE. Within the
    chosen bin both E and xi are sampled uniformly (xi is returned for
    completeness / commit-3 reuse but is not needed for the scalar rate).
    """
    dE = np.diff(EB)
    dA = np.diff(AB)
    w = (F_zone * dA[:, None] * dE[None, :]).ravel()
    tot = w.sum()
    if tot <= 0:
        return None, None
    p = w / tot
    n_xi, n_E = F_zone.shape
    flat = rng.choice(p.size, size=nsamp, p=p)
    a_idx, e_idx = np.divmod(flat, n_E)
    # uniform within the selected bin
    E_s = EB[e_idx] + rng.random(nsamp) * (EB[e_idx + 1] - EB[e_idx])     # eV
    xi_s = AB[a_idx] + rng.random(nsamp) * (AB[a_idx + 1] - AB[a_idx])
    v_s = btk.speed_from_energy_keV(E_s / 1.0e3, btk.M_D_AMU)             # m/s
    return v_s, xi_s


def sigv_average(v_fast, ti_eV, reaction, m_fast_amu, m_th_amu, rng):
    """<sigma * v_rel> in m^3/s, averaged over an isotropic Maxwellian target.

    v_fast : (nsamp,) fast-ion speeds in m/s. One thermal partner is drawn per
    fast-ion sample from a 3-D Maxwellian at T_i. Because the thermal
    distribution is isotropic, the fast-ion direction is irrelevant to the
    scalar rate, so v_fast is placed along z.
    """
    n = v_fast.size
    sigma_v_th = np.sqrt(ti_eV * btk.EV_J / (m_th_amu * btk.AMU_KG))      # m/s
    v_f = np.zeros((n, 3))
    v_f[:, 2] = v_fast
    v_t = rng.normal(0.0, sigma_v_th, size=(n, 3))
    E_cm = btk.relative_cm_energy_keV(v_f, v_t, m_fast_amu, m_th_amu)     # keV
    if reaction == "DD_n":
        sig_mb = btk.sigma_dd_n_mb(E_cm)
    elif reaction in ("DT", "TD"):
        sig_mb = btk.sigma_dt_mb(E_cm)
    else:
        raise ValueError(f"unsupported reaction {reaction!r}")
    v_rel = np.linalg.norm(v_f - v_t, axis=-1)                           # m/s
    return float(np.mean(sig_mb * btk.MB_M2 * v_rel))                    # m^3/s


def zone_reactivity(fi, n_fast, n_th, ti, reaction, m_fast_amu, m_th_amu,
                    nsamp, rng):
    """Per-zone 4*pi reactivity eps_BT (1/cm3/s) for one reaction channel.

        eps = n_fast * n_th * <sigma v_rel>
    with n in 1/cm3 and <sigma v_rel> converted m^3/s -> cm^3/s.
    """
    F = fi["F"]
    n_zone = F.shape[0]
    eps = np.zeros(n_zone)
    for z in range(n_zone):
        if n_fast[z] <= 0 or n_th[z] <= 0 or ti[z] <= 0:
            continue
        v_fast, _ = sample_fast_speeds(F[z], fi["EB"], fi["AB"],
                                       fi["E"], fi["A"], nsamp, rng)
        if v_fast is None:
            continue
        sv_m3s = sigv_average(v_fast, ti[z], reaction,
                              m_fast_amu, m_th_amu, rng)
        eps[z] = n_fast[z] * n_th[z] * sv_m3s * 1.0e6   # m^3 -> cm^3
    return eps


# ----------------------------------------------------------------------
# Reporting / acceptance test
# ----------------------------------------------------------------------

def compare(label, eps, ref, bmvol):
    """Print zone-wise and volume-integrated comparison vs NUBEAM."""
    rate_mine = float(np.sum(eps * bmvol))
    rate_ref = float(np.sum(ref * bmvol))
    print(f"\n--- {label} channel ---")
    print(f"  volume-integrated rate (in-house) : {rate_mine:.4e} 1/s")
    print(f"  volume-integrated rate (NUBEAM)    : {rate_ref:.4e} 1/s")
    if rate_ref > 0:
        print(f"  ratio in-house / NUBEAM            : {rate_mine / rate_ref:.4f}")
    # per-zone where NUBEAM is non-negligible
    m = ref > ref.max() * 1e-3
    if m.any():
        r = eps[m] / ref[m]
        print(f"  per-zone ratio (ref>1e-3*max, n={m.sum()}): "
              f"min {r.min():.3f}  median {np.median(r):.3f}  max {r.max():.3f}")
    return rate_mine, rate_ref


def make_plot(fi, results, save=None):
    import matplotlib
    if save is not None:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_chan = len(results)
    fig, axes = plt.subplots(2, n_chan, figsize=(6 * n_chan, 9), squeeze=False)
    r = fi["r2d"] / 100.0
    z = fi["z2d"] / 100.0
    for j, (label, eps, ref) in enumerate(results):
        # row 0: scatter of in-house eps in (R,Z)
        ax = axes[0][j]
        sc = ax.scatter(r, z, c=eps, s=14, cmap="inferno")
        ax.set_aspect("equal")
        ax.set_title(f"{label}: in-house eps_BT (R,Z)")
        ax.set_xlabel("R [m]"); ax.set_ylabel("Z [m]")
        fig.colorbar(sc, ax=ax, label="1/cm3/s")
        # row 1: per-zone in-house vs NUBEAM scatter
        ax = axes[1][j]
        good = ref > ref.max() * 1e-4
        ax.loglog(ref[good], eps[good], ".", ms=4, alpha=0.6)
        lim = [min(ref[good].min(), eps[good][eps[good] > 0].min()),
               max(ref[good].max(), eps[good].max())]
        ax.plot(lim, lim, "k--", lw=1, label="y = x")
        ax.set_xlabel("NUBEAM BTN [1/cm3/s]")
        ax.set_ylabel("in-house eps_BT [1/cm3/s]")
        ax.set_title(f"{label}: per-zone agreement")
        ax.legend()
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=130)
        print(f"\nSaved plot to {save}")
    else:
        plt.show()


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("pulse", type=int, help="JET pulse number, e.g. 104614")
    p.add_argument("run_suffix", help="TRANSP run suffix, e.g. M30")
    p.add_argument("--idx", type=int, default=None,
                   help="FBM time-window index (default: first available)")
    p.add_argument("--data-dir", default=None,
                   help="explicit <base> for <base>/<pulse>/<run_suffix>")
    p.add_argument("--nsamp", type=int, default=20000,
                   help="MC samples per zone (default 20000)")
    p.add_argument("--fast-norm", choices=["bdens", "ntot"], default="bdens",
                   help="fast-ion density normalization: 'bdens' rescales F to "
                        "the full BDENS_D profile (matches NUBEAM, default); "
                        "'ntot' keeps the raw F/NTOT_D_NBI normalization")
    p.add_argument("--seed", type=int, default=12345, help="RNG seed")
    p.add_argument("--plot", action="store_true", help="show diagnostic plot")
    p.add_argument("--save", default=None, help="save plot to file")
    p.add_argument("--no-plot", action="store_true",
                   help="force no plotting even with --save")
    args = p.parse_args(argv)

    run_id = f"{args.pulse}{args.run_suffix}"
    run_dir = find_run_dir(args.pulse, args.run_suffix, args.data_dir)
    idxs = list_fbm_indices(run_dir, run_id)
    if not idxs:
        print(f"No {run_id}_fi_*.cdf found in {run_dir}", file=sys.stderr)
        return 1
    idx = args.idx if args.idx is not None else idxs[0]
    if idx not in idxs:
        print(f"idx {idx} not in available {idxs}", file=sys.stderr)
        return 1

    fi_path = run_dir / f"{run_id}_fi_{idx}.cdf"
    neut_path = run_dir / f"{run_id}_neut_{idx}.cdf"
    cdf_path = run_dir / f"{run_id}.CDF"

    print(f"Run dir : {run_dir}")
    print(f"FBM idx : {idx}")
    fi = read_fi_distribution(fi_path)
    neut = read_neut_rates(neut_path)
    therm = read_thermal_profiles(cdf_path, fi["time"], fi["x2d"])

    print(f"_fi  time : {fi['time']:.4f} s")
    print(f"_neut time: {neut['time']:.4f} s")
    print(f"thermal slice TIME3[{therm['tind']}] = {therm['t_used']:.4f} s")
    print(f"n_zones   : {fi['F'].shape[0]}  "
          f"(n_xi={fi['F'].shape[1]}, n_E={fi['F'].shape[2]})")

    # fast-ion density + normalization cross-check
    n_fast = fast_ion_density(fi)
    n_from_F = float(np.sum(n_fast * fi["bmvol"]))
    n_bdens = float(np.sum(therm["bdens_d"] * fi["bmvol"]))
    print(f"\nNormalization check:")
    print(f"  sum(n_fast*BMVOL)  = {n_from_F:.4e}")
    print(f"  NTOT_D_NBI         = {fi['ntot']:.4e}  "
          f"(ratio {n_from_F / fi['ntot']:.4f})")
    print(f"  sum(BDENS_D*BMVOL) = {n_bdens:.4e}  "
          f"(F undercounts by {n_from_F / n_bdens:.3f})")

    if args.fast_norm == "bdens":
        n_fast, rows = renormalize_to_bdens(n_fast, fi["x2d"],
                                            fi["bmvol"], therm["bdens_d"])
        f_arr = np.array(list(rows.values()))
        print(f"  fast-norm = bdens: per-row scale {f_arr.min():.3f}..{f_arr.max():.3f}")
    else:
        print(f"  fast-norm = ntot (raw F normalization)")

    rng = np.random.default_rng(args.seed)

    # Channels to validate. Beams are D-only for these runs, so the fast
    # reactant is always D; the thermal target sets the channel.
    channels = []
    if "DD" in neut:
        channels.append(("DD", "DD_n", btk.M_D_AMU, therm["nd"], neut["DD"]))
    if "DT" in neut and np.any(therm["nt"] > 0):
        channels.append(("DT", "DT", btk.M_T_AMU, therm["nt"], neut["DT"]))

    results = []
    for label, reaction, m_th, n_th, ref in channels:
        eps = zone_reactivity(fi, n_fast, n_th, therm["ti"],
                              reaction, btk.M_D_AMU, m_th,
                              args.nsamp, rng)
        compare(label, eps, ref, fi["bmvol"])
        results.append((label, eps, ref))

    if (args.plot or args.save) and not args.no_plot and results:
        make_plot(fi, results, save=args.save)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
