#!/usr/bin/env python
"""
bt_kinematics.py
----------------

Pure-physics module for beam-target fusion neutron production.

Provides
~~~~~~~~
* Bosch-Hale cross sections sigma(E_cm) for the BT neutron channels
  relevant to JET DT plasmas
      DD_n  :  D + D -> 3He + n + 3.27 MeV    (Bosch-Hale 1992)
      DT    :  D + T -> 4He + n + 17.59 MeV   (Bosch-Hale 1992)
      TD    :  same reaction as DT (cross section depends only on E_cm,
               not on which reactant is "fast"); kinematics differ in lab
      TT    :  T + T -> 4He + 2n + 11.33 MeV  (NOT implemented in this
               commit; 3-body final state, no Bosch-Hale parameterisation)
* CM-frame neutron emission distribution
      isotropic_cm()  -- DD and DT are approximately isotropic in CM at
                         beam energies up to ~200 keV (intrinsic
                         anisotropy ~5-10 %; lab anisotropy is dominated
                         by the CM->lab boost). Hooks left for an
                         anisotropic Legendre expansion later.
* Classical 2-body kinematics
      neutron_lab_velocity(...) -- given lab velocities of the two
      reactants and a CM-frame neutron emission direction, return the
      lab-frame neutron velocity vector. No relativistic correction
      (~0.1 % at 2.45 MeV).

No TRANSP / NUBEAM / equilibrium dependencies. Self-contained.

References
~~~~~~~~~~
* Bosch, H.-S. and Hale, G. M., "Improved formulas for fusion
  cross-sections and thermal reactivities", Nucl. Fusion 32, 611 (1992).
  Cross-section parameterisation in Table VII, valid E_cm ranges below.
* Pankin et al. (2025) Section 6.1.1 and Table 8 -- list of fusion
  product channels modelled by NUBEAM.
"""

from __future__ import annotations

import numpy as np


# ----------------------------------------------------------------------
# Physical constants
# ----------------------------------------------------------------------

AMU_KG    = 1.66053906660e-27    # kg
EV_J      = 1.602176634e-19      # J / eV
KEV_J     = 1.602176634e-16      # J / keV
MEV_J     = 1.602176634e-13      # J / MeV
BARN_M2   = 1.0e-28              # m^2 / barn
MB_M2     = 1.0e-31              # m^2 / mb

# Particle masses (amu, AME2020)
M_D_AMU   = 2.0141017781
M_T_AMU   = 3.0160492779
M_HE3_AMU = 3.0160293201
M_HE4_AMU = 4.0026032541
M_N_AMU   = 1.0086649157

# Q-values in MeV
Q_DD_N    = 3.26890
Q_DT      = 17.58928
Q_TT      = 11.33178


# ----------------------------------------------------------------------
# Bosch-Hale cross sections
# ----------------------------------------------------------------------
#
# Bosch & Hale (1992), Table VII. Define
#     S(E)     = (A1 + E (A2 + E (A3 + E (A4 + E A5))))
#                / (1 + E (B1 + E (B2 + E (B3 + E B4))))
#     sigma(E) = S(E) / (E exp(B_G / sqrt(E)))    [mb]
# with E in keV (CM kinetic energy of the reactants).

_BH = {
    # D(d,n)3He, valid E_cm = 0.5 .. 4900 keV
    "DD_n": dict(
        E_range=(0.5, 4900.0),
        B_G=31.3970,
        A=(5.3701e+4, 3.3027e+2, -1.2706e-1, 2.9327e-5, -2.5151e-9),
        B=(0.0, 0.0, 0.0, 0.0),
    ),
    # T(d,n)4He, low-energy regime 0.5 .. 550 keV
    "DT_low": dict(
        E_range=(0.5, 550.0),
        B_G=34.3827,
        A=(6.927e+4, 7.454e+8, 2.050e+6, 5.2002e+4, 0.0),
        B=(6.38e+1, -9.95e-1, 6.981e-5, 1.728e-4),
    ),
    # T(d,n)4He, high-energy regime 550 .. 4700 keV
    "DT_high": dict(
        E_range=(550.0, 4700.0),
        B_G=34.3827,
        A=(-1.4714e+6, 0.0, 0.0, 0.0, 0.0),
        B=(-8.4127e-3, 4.7983e-6, -1.0748e-9, 8.5184e-14),
    ),
}


def _bh_sigma_mb(E_keV, coeffs):
    """Bosch-Hale sigma in mb, scalar or array."""
    A = coeffs["A"]
    B = coeffs["B"]
    BG = coeffs["B_G"]
    E = np.asarray(E_keV, dtype=float)
    # Horner-style polynomial
    num = A[0] + E*(A[1] + E*(A[2] + E*(A[3] + E*A[4])))
    den = 1.0 + E*(B[0] + E*(B[1] + E*(B[2] + E*B[3])))
    S = num / den
    return S / (E * np.exp(BG / np.sqrt(E)))


def sigma_dd_n_mb(E_cm_keV):
    """D(d,n)3He cross section in mb. E_cm in keV.

    Bosch-Hale 1992 parameterisation, valid 0.5 .. 4900 keV.
    Outside this range the returned value is set to zero.
    """
    E = np.asarray(E_cm_keV, dtype=float)
    lo, hi = _BH["DD_n"]["E_range"]
    mask = (E >= lo) & (E <= hi)
    out = np.zeros_like(E)
    if mask.any():
        out[mask] = _bh_sigma_mb(E[mask], _BH["DD_n"])
    return out if out.ndim else float(out)


def sigma_dt_mb(E_cm_keV):
    """T(d,n)4He cross section in mb. E_cm in keV.

    Bosch-Hale 1992 parameterisation. Two regimes joined at 550 keV.
    Valid 0.5 .. 4700 keV; outside, zero.

    Same parameterisation is used for the TD channel (fast triton on
    thermal deuterium): the cross section depends only on E_cm.
    """
    E = np.asarray(E_cm_keV, dtype=float)
    out = np.zeros_like(E)
    low = (E >= _BH["DT_low"]["E_range"][0]) & (E <= _BH["DT_low"]["E_range"][1])
    high = (E > _BH["DT_low"]["E_range"][1]) & (E <= _BH["DT_high"]["E_range"][1])
    if low.any():
        out[low] = _bh_sigma_mb(E[low], _BH["DT_low"])
    if high.any():
        out[high] = _bh_sigma_mb(E[high], _BH["DT_high"])
    return out if out.ndim else float(out)


def sigma_tt_mb(E_cm_keV):
    """T(t,2n)4He cross section in mb.  *Not implemented in this commit.*

    T+T -> 4He + 2n is a 3-body final state with a continuous neutron
    spectrum (0 .. ~9 MeV) and is not part of Bosch-Hale 1992. To add it
    later, use either an ENDF/B-VIII tabulation or the Drosg
    parameterisation. For DTE3 (D-only beams, e.g. JET pulse 104614) the
    TT contribution is negligible; for DTE2 (D and T beams) it can reach
    a few percent.
    """
    raise NotImplementedError(
        "TT cross section not implemented yet (placeholder for future "
        "ENDF/Drosg-based addition)."
    )


# ----------------------------------------------------------------------
# Reaction metadata helpers
# ----------------------------------------------------------------------

# (m_fast, m_thermal, m_other_product, Q_MeV, sigma_function)
# 'm_other_product' is the partner of the neutron in the final state.
_REACTIONS = {
    "DD_n": dict(m_fast=M_D_AMU, m_th=M_D_AMU,  m_other=M_HE3_AMU,
                 Q=Q_DD_N, sigma=sigma_dd_n_mb),
    "DT":   dict(m_fast=M_D_AMU, m_th=M_T_AMU,  m_other=M_HE4_AMU,
                 Q=Q_DT,   sigma=sigma_dt_mb),
    "TD":   dict(m_fast=M_T_AMU, m_th=M_D_AMU,  m_other=M_HE4_AMU,
                 Q=Q_DT,   sigma=sigma_dt_mb),
    "TT":   dict(m_fast=M_T_AMU, m_th=M_T_AMU,  m_other=M_HE4_AMU,
                 Q=Q_TT,   sigma=sigma_tt_mb),
}


def reaction_info(reaction):
    """Return the metadata dict for one of 'DD_n', 'DT', 'TD', 'TT'."""
    try:
        return _REACTIONS[reaction]
    except KeyError as exc:
        raise KeyError(
            f"Unknown reaction {reaction!r}; choose from {list(_REACTIONS)}"
        ) from exc


# ----------------------------------------------------------------------
# Kinematics
# ----------------------------------------------------------------------

def relative_cm_energy_keV(v1, v2, m1_amu, m2_amu):
    """CM kinetic energy of two reactants given lab velocities (m/s).

    E_cm = 0.5 * mu * |v1 - v2|^2     with    mu = m1 m2 / (m1 + m2)

    Accepts vector or array inputs; the last axis is the 3 components.
    Returns E_cm in keV.
    """
    v1 = np.asarray(v1, dtype=float)
    v2 = np.asarray(v2, dtype=float)
    mu_kg = (m1_amu * m2_amu) / (m1_amu + m2_amu) * AMU_KG
    dv = v1 - v2
    v2sq = np.sum(dv * dv, axis=-1)
    return 0.5 * mu_kg * v2sq / KEV_J


def isotropic_cm(n, rng=None):
    """Sample ``n`` isotropic directions in the CM frame.

    Returns
    -------
    cos_theta : (n,) array in [-1, 1]
    phi       : (n,) array in [0, 2 pi)
    """
    rng = np.random.default_rng() if rng is None else rng
    cos_theta = 2.0 * rng.random(n) - 1.0
    phi = 2.0 * np.pi * rng.random(n)
    return cos_theta, phi


def _orthonormal_frame(n_axis):
    """Build an orthonormal triad (e_perp1, e_perp2, n_axis).

    n_axis : (..., 3) unit vectors. Picks a perpendicular by cross with the
    coordinate axis least aligned to n_axis, then completes the triad.
    """
    n_axis = np.asarray(n_axis, dtype=float)
    n_axis = n_axis / np.linalg.norm(n_axis, axis=-1, keepdims=True)
    # Pick reference axis: x if |n_z| > 0.9 else z (avoid degeneracy)
    z_ref = np.where(
        np.abs(n_axis[..., 2:3]) > 0.9,
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 0.0, 1.0]),
    )
    e1 = np.cross(n_axis, z_ref)
    e1 = e1 / np.linalg.norm(e1, axis=-1, keepdims=True)
    e2 = np.cross(n_axis, e1)
    return e1, e2, n_axis


def neutron_lab_velocity(v_fast, v_thermal, reaction,
                         cos_theta_cm, phi_cm):
    """Lab-frame neutron velocity vector for a single fusion event.

    Parameters
    ----------
    v_fast, v_thermal : array_like, shape (..., 3)
        Lab velocities of the fast and thermal reactants in m/s.
    reaction : str
        One of 'DD_n', 'DT', 'TD'. ('TT' not supported -- 3-body.)
    cos_theta_cm, phi_cm : array_like
        CM-frame emission direction of the neutron.  ``cos_theta_cm`` is
        measured relative to the *fast reactant's velocity in the CM
        frame*; for isotropic emission this choice is irrelevant.

    Returns
    -------
    v_n_lab : ndarray, shape (..., 3)
        Lab-frame neutron velocity in m/s.

    Notes
    -----
    Energy conservation in the CM frame:
        E_cm_final = E_cm_kin + Q
        E_n_cm     = E_cm_final * m_other / (m_n + m_other)
        |u_n|      = sqrt(2 * E_n_cm / m_n)
    Then v_n_lab = u_n_vec + v_CM, with v_CM = (m1 v1 + m2 v2)/(m1+m2).
    """
    info = reaction_info(reaction)
    if reaction == "TT":
        raise NotImplementedError("TT has a 3-body final state.")
    m1 = info["m_fast"] * AMU_KG
    m2 = info["m_th"]   * AMU_KG
    mn = M_N_AMU         * AMU_KG
    mo = info["m_other"] * AMU_KG
    Q_J = info["Q"] * MEV_J

    v1 = np.asarray(v_fast,    dtype=float)
    v2 = np.asarray(v_thermal, dtype=float)
    v_cm = (m1 * v1 + m2 * v2) / (m1 + m2)
    dv = v1 - v2
    mu = m1 * m2 / (m1 + m2)
    E_kin_cm = 0.5 * mu * np.sum(dv * dv, axis=-1)
    E_tot_cm = E_kin_cm + Q_J

    # Neutron CM speed from 2-body energy partition
    E_n_cm = E_tot_cm * mo / (mn + mo)
    u_n = np.sqrt(2.0 * E_n_cm / mn)              # scalar speed in m/s

    # Frame: polar axis along v_fast_in_CM (the "beam" direction in CM)
    v1_cm = v1 - v_cm
    e1, e2, n_axis = _orthonormal_frame(v1_cm)

    cos_t = np.asarray(cos_theta_cm, dtype=float)
    sin_t = np.sqrt(np.clip(1.0 - cos_t * cos_t, 0.0, None))
    phi   = np.asarray(phi_cm, dtype=float)
    dirn = (cos_t[..., None] * n_axis
            + sin_t[..., None] * (np.cos(phi)[..., None] * e1
                                  + np.sin(phi)[..., None] * e2))
    return v_cm + u_n[..., None] * dirn


def neutron_lab_energy_keV(v_n_lab):
    """Convert a lab neutron velocity vector to kinetic energy in keV."""
    v_n_lab = np.asarray(v_n_lab, dtype=float)
    vsq = np.sum(v_n_lab * v_n_lab, axis=-1)
    return 0.5 * (M_N_AMU * AMU_KG) * vsq / KEV_J


def speed_from_energy_keV(E_keV, m_amu):
    """Non-relativistic |v| in m/s from kinetic energy in keV."""
    return np.sqrt(2.0 * np.asarray(E_keV) * KEV_J / (m_amu * AMU_KG))


# ----------------------------------------------------------------------
# Unit tests / demo (run as a script)
# ----------------------------------------------------------------------

def _test_cross_sections():
    """Spot-check Bosch-Hale sigma at reference CM energies."""
    print("\n--- Bosch-Hale cross sections ---")
    # Reference values cross-checked by hand-evaluating the Bosch-Hale
    # parameterisation (Bosch & Hale 1992, Table VII). They are
    # regression anchors, not independent measurements.
    checks = [
        ("DD_n", sigma_dd_n_mb,  50.0,   16.5,  1.0),    # rising
        ("DD_n", sigma_dd_n_mb, 100.0,   37.0,  2.0),
        ("DT",   sigma_dt_mb,    50.0, 4220.0, 50.0),    # below peak
        ("DT",   sigma_dt_mb,    64.0, 5063.0, 50.0),    # peak in E_cm
        ("DT",   sigma_dt_mb,   200.0, 1138.0, 50.0),    # well past peak
    ]
    ok = True
    for name, fn, E, ref, tol in checks:
        s = float(fn(E))
        good = abs(s - ref) < tol
        status = "OK" if good else "FAIL"
        ok &= good
        print(f"  sigma_{name:4s}({E:6.1f} keV) = {s:9.2f} mb  "
              f"(expected {ref:7.1f} +/- {tol:5.1f})  [{status}]")
    return ok


def _test_kinematics():
    """Sanity-check kinematic limits."""
    print("\n--- Kinematics: 100 keV D on cold D, isotropic CM ---")
    rng = np.random.default_rng(0)
    n = 50000
    E_d_lab_keV = 100.0
    v_d = np.array([speed_from_energy_keV(E_d_lab_keV, M_D_AMU), 0.0, 0.0])
    v_t = np.zeros(3)                                 # cold target
    cos_t, phi = isotropic_cm(n, rng=rng)
    v1 = np.broadcast_to(v_d, (n, 3))
    v2 = np.broadcast_to(v_t, (n, 3))
    v_n = neutron_lab_velocity(v1, v2, "DD_n", cos_t, phi)
    E_n = neutron_lab_energy_keV(v_n)

    # Analytical reference (cold target, fast beam):
    #   E_cm_kin = E_beam * m_target / (m_beam + m_target)
    # For D on cold D this is just E_d_lab / 2.
    E_cm_kin = E_d_lab_keV * M_D_AMU / (M_D_AMU + M_D_AMU)
    E_tot_cm = E_cm_kin + Q_DD_N * 1e3                   # keV
    E_n_cm = E_tot_cm * M_HE3_AMU / (M_N_AMU + M_HE3_AMU)
    u_n = speed_from_energy_keV(E_n_cm, M_N_AMU)
    v_cm = (M_D_AMU * v_d) / (M_D_AMU + M_D_AMU)        # half v_d along x
    speed_cm = np.linalg.norm(v_cm)
    # Lab neutron energy in [(u_n - speed_cm)^2, (u_n + speed_cm)^2] * m_n / 2
    E_min = 0.5 * (M_N_AMU * AMU_KG) * (u_n - speed_cm)**2 / KEV_J
    E_max = 0.5 * (M_N_AMU * AMU_KG) * (u_n + speed_cm)**2 / KEV_J
    print(f"  E_d_lab = {E_d_lab_keV} keV  -> E_cm_kin = {E_cm_kin:.2f} keV")
    print(f"  Predicted neutron lab energy window: "
          f"{E_min:.1f} .. {E_max:.1f} keV")
    print(f"  Simulated mean E_n = {E_n.mean():.1f} keV  "
          f"(min {E_n.min():.1f},  max {E_n.max():.1f})")
    # Sampling won't reach the analytical extremes exactly; require
    # min/max within 2 keV (well under the 700 keV window).
    band_ok = (abs(E_n.min() - E_min) < 2.0 and abs(E_n.max() - E_max) < 2.0)
    mean_ok = abs(E_n.mean() - 0.5 * (E_min + E_max)) < 2.0
    status = "OK" if (band_ok and mean_ok) else "FAIL"
    print(f"  status: {status}")
    return band_ok and mean_ok


def _test_kinematics_dt():
    """Sanity-check DT: 100 keV D on cold T -> neutron at ~14.1 MeV."""
    print("\n--- Kinematics: 100 keV D on cold T, isotropic CM ---")
    rng = np.random.default_rng(0)
    n = 50000
    v_d = np.array([speed_from_energy_keV(100.0, M_D_AMU), 0.0, 0.0])
    v_t = np.zeros(3)
    cos_t, phi = isotropic_cm(n, rng=rng)
    v1 = np.broadcast_to(v_d, (n, 3))
    v2 = np.broadcast_to(v_t, (n, 3))
    v_n = neutron_lab_velocity(v1, v2, "DT", cos_t, phi)
    E_n = neutron_lab_energy_keV(v_n)
    # Expected lab spread around 14.1 MeV
    mu = M_D_AMU * M_T_AMU / (M_D_AMU + M_T_AMU)
    E_cm_kin = 100.0 * M_T_AMU / (M_D_AMU + M_T_AMU)
    E_tot_cm = E_cm_kin + Q_DT * 1e3                    # keV
    E_n_cm = E_tot_cm * M_HE4_AMU / (M_N_AMU + M_HE4_AMU)
    u_n = speed_from_energy_keV(E_n_cm, M_N_AMU)
    v_cm = M_D_AMU * v_d / (M_D_AMU + M_T_AMU)
    speed_cm = np.linalg.norm(v_cm)
    E_min = 0.5 * (M_N_AMU * AMU_KG) * (u_n - speed_cm)**2 / KEV_J
    E_max = 0.5 * (M_N_AMU * AMU_KG) * (u_n + speed_cm)**2 / KEV_J
    print(f"  Predicted lab neutron window: "
          f"{E_min/1000.0:.3f} .. {E_max/1000.0:.3f} MeV")
    print(f"  Simulated mean E_n = {E_n.mean()/1000.0:.3f} MeV  "
          f"(min {E_n.min()/1000.0:.3f},  max {E_n.max()/1000.0:.3f})")
    band_ok = (abs(E_n.min() - E_min) < 20.0
               and abs(E_n.max() - E_max) < 20.0)
    mean_ok = abs(E_n.mean() - 0.5 * (E_min + E_max)) < 20.0
    status = "OK" if (band_ok and mean_ok) else "FAIL"
    print(f"  status: {status}")
    return band_ok and mean_ok


def _plot_cross_sections(save=None):
    """Plot sigma(E_cm) for DD-n and DT."""
    import matplotlib.pyplot as plt
    E = np.logspace(0, np.log10(4000.0), 400)            # 1 .. 4000 keV
    s_dd = sigma_dd_n_mb(E)
    s_dt = sigma_dt_mb(E)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(E, s_dd, "C0-", label="D(d,n)$^3$He")
    ax.loglog(E, s_dt, "C3-", label="T(d,n)$^4$He  (also TD)")
    ax.set_xlabel(r"$E_{cm}$ [keV]")
    ax.set_ylabel(r"$\sigma$ [mb]")
    ax.set_title("Bosch-Hale cross sections")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(loc="lower right")
    # Mark peak of DT
    iE = np.argmax(s_dt)
    ax.axvline(E[iE], color="C3", ls=":", alpha=0.4)
    ax.text(E[iE]*1.1, s_dt[iE]*0.9,
            f"  peak: {s_dt[iE]:.0f} mb at {E[iE]:.0f} keV",
            color="C3", fontsize=9)
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=150, bbox_inches="tight")
        print(f"Saved {save}")
    else:
        plt.show()


def main():
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--plot", action="store_true",
                   help="Also plot sigma(E_cm) for the implemented reactions.")
    p.add_argument("--save", default=None,
                   help="Save sigma plot to this file instead of show().")
    args = p.parse_args()
    ok1 = _test_cross_sections()
    ok2 = _test_kinematics()
    ok3 = _test_kinematics_dt()
    overall = "PASS" if (ok1 and ok2 and ok3) else "FAIL"
    print(f"\n=== Overall: {overall} ===")
    if args.plot or args.save:
        _plot_cross_sections(save=args.save)


if __name__ == "__main__":
    main()
