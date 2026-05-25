#!/usr/bin/env python
"""
estimate_outer_th_fraction.py
-----------------------------

Estimate the fraction of the thermal-neutron line-of-sight (LOS) signal that
originates *outside* a central plasma volume around the magnetic axis.

The script is meant to quantify the error one makes when assuming that the
whole thermal-neutron signal recorded by a vertical LOS comes only from the
central core of the plasma.

Physical / geometrical assumptions (documented here so they are easy to find):
    * circular flux surfaces in the (R,Z) poloidal plane;
    * a single vertical LOS through the magnetic axis, R = RMAG, varying Z;
    * a constant toroidal LOS width, w = 0.4 m, used to define both the
      central cylinder volume V0 = w*pi*r0^2 and the geometrical weighting
      w_geom(Z) = w / (2*pi*|Z - ZMAG|) for off-axis surfaces;
    * up-down symmetry of the emissivity profile around Z = ZMAG (we
      integrate only one side and multiply by 2);
    * the finite detector solid angle is omitted - it cancels exactly in
      the relative fraction f_outer = Crest / (C0 + Crest).

Numerical safeguards:
    * the central exclusion |Z - ZMAG| <= r0 is removed from the outer
      integral, which also removes the 1/r geometric singularity at r=0.

Usage
-----
    python estimate_outer_th_fraction.py 104614 M30 dda=eftp uid=jetppf seq=405
    python estimate_outer_th_fraction.py 104614 M30 --dda eftp --uid jetppf --seq 405 \\
            --time 49.0 --plot

The convenience ``key=value`` form (e.g. ``dda='eftp'``) is also accepted.
"""

import sys
sys.append('../../')
import argparse

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

import profiles as ps
from change_rho import RZ_to_rhot


# -----------------------------------------------------------------------------
# Geometric constants (fixed by the problem statement)
# -----------------------------------------------------------------------------
R0_CENTRAL = 0.2   # central poloidal radius around the axis [m]
W_TOROIDAL = 0.4   # effective toroidal width of the LOS [m]


# -----------------------------------------------------------------------------
# CLI helpers
# -----------------------------------------------------------------------------
def _preprocess_argv(argv):
    """Convert ``key=value`` style tokens to ``--key value`` pairs.

    Accepts the user-friendly form shown in the docstring
    (``dda='eftp' uid='jetppf' seq=405``) and rewrites it to the canonical
    argparse form. Positional arguments (no '=') are left untouched.
    """
    out = []
    for tok in argv:
        if tok.startswith('-') or '=' not in tok:
            out.append(tok)
            continue
        k, v = tok.split('=', 1)
        # strip whitespace and matching surrounding quotes
        v = v.strip().strip("'").strip('"')
        out.extend([f'--{k}', v])
    return out


def parse_args(argv=None):
    """Build the CLI parser and return the parsed namespace."""
    if argv is None:
        argv = sys.argv[1:]
    argv = _preprocess_argv(argv)

    parser = argparse.ArgumentParser(
        description='Estimate the outer-region fraction of the thermal '
                    'neutron LOS signal.',
    )
    parser.add_argument('pulse', type=int, help='JET pulse number')
    parser.add_argument('runid', type=str, help='TRANSP runid (e.g. M30)')
    parser.add_argument('--dda', type=str, default='eftp',
                        help='Equilibrium PPF DDA (default: eftp)')
    parser.add_argument('--uid', type=str, default='jetppf',
                        help='Equilibrium PPF UID (default: jetppf)')
    parser.add_argument('--seq', type=int, default=0,
                        help='Equilibrium PPF sequence (default: 0)')
    parser.add_argument('-t', '--time', type=float, default=None,
                        help='Time in JET convention (s, t>40). '
                             'If omitted, the midpoint of the TRANSP run is used.')
    parser.add_argument('-n', '--nz', type=int, default=400,
                        help='Number of LOS Z samples (default 400)')
    parser.add_argument('--plot', action='store_true',
                        help='Show diagnostic plots.')
    return parser.parse_args(argv)


# -----------------------------------------------------------------------------
# Data extraction / LOS construction
# -----------------------------------------------------------------------------
def get_transp_slice(tr, time_jet):
    """Return THNTX, DVOL, rhot and TRANSP time index at the JET time *time_jet*.

    Parameters
    ----------
    tr : ps.Transp
        TRANSP wrapper.
    time_jet : float
        Requested time in JET convention (s).

    Returns
    -------
    thntx : 1D array, [n/s/m^3] on the rhot grid (`X`)
    dvol  : 1D array, [m^3] differential flux-surface volume on the rhot grid
    x_rhot : 1D array, rhot = sqrt(normalized toroidal flux)
    trind : int, TRANSP time index actually used
    """
    # TRANSP times are 0-based; JET convention adds 40 s.
    t_jet_grid = tr.t + 40.0
    trind = int(np.abs(t_jet_grid - time_jet).argmin())

    thntx_all = tr.profile('THNTX')
    dvol_all = tr.profile('DVOL')
    x_all = tr.x

    # All three may be 1D (static grid) or 2D (n_t, n_x).
    thntx = thntx_all[trind] if thntx_all.ndim == 2 else thntx_all
    dvol = dvol_all[trind] if dvol_all.ndim == 2 else dvol_all
    x_rhot = x_all[trind] if x_all.ndim == 2 else x_all

    return thntx, dvol, x_rhot, trind


def build_los(eq, tind_eq, n_z):
    """Build the vertical LOS sample points along R = RMAG.

    The Z extent is taken symmetric around ZMAG and reaches to the maximum
    vertical extent of the LCFS (whichever side is larger), so that the
    central point Z = ZMAG is always present and the upper/lower limbs are
    sampled identically.

    Returns
    -------
    R, Z : 1D arrays of LOS sample coordinates (length n_z)
    Rmag, Zmag : magnetic-axis coordinates [m] at this equilibrium time
    """
    Rmag = float(eq._Rmag[tind_eq])
    Zmag = float(eq._Zmag[tind_eq])

    # `_Zbnd` is shape (n_bnd_points, n_t); pull the LCFS for this slice.
    Zbnd = np.asarray(eq._Zbnd[:, tind_eq])
    z_top = float(np.nanmax(Zbnd))
    z_bot = float(np.nanmin(Zbnd))
    dz = max(z_top - Zmag, Zmag - z_bot)

    Z = np.linspace(Zmag - dz, Zmag + dz, int(n_z))
    R = np.full_like(Z, Rmag)
    return R, Z, Rmag, Zmag


def remap_profile_on_z(x_rhot, profile_on_x, R, Z, eq, time_jet):
    """Remap a TRANSP profile defined on rhot onto the LOS (R, Z).

    Uses ``change_rho.RZ_to_rhot`` to convert each LOS point to rhot, then
    PCHIP-interpolates the profile from the rhot grid. The rhot grid is
    padded with the axis value at rhot=0 and zero at rhot=1 so the LOS
    point at the magnetic axis (rhot~0) and any point outside the LCFS
    receive a sensible value instead of NaN.

    Returns
    -------
    profile_on_z : 1D array, same shape as Z
    rhot_on_z    : 1D array, rhot value at every LOS sample (clipped to [0, 1])
    """
    rhot_on_z = RZ_to_rhot(R, Z, eq, time_jet, method='cubic', clip_psin=True)
    rhot_on_z = np.asarray(rhot_on_z).ravel()

    # Sort x ascending for the interpolator and remove duplicate knots.
    order = np.argsort(x_rhot)
    xs = np.asarray(x_rhot)[order]
    ys = np.asarray(profile_on_x)[order]
    uniq = np.concatenate(([True], np.diff(xs) > 0))
    xs = xs[uniq]
    ys = ys[uniq]

    # Pad the rhot grid so we can evaluate everywhere in [0, 1]:
    #   * rhot = 0  -> assume the value of the innermost TRANSP zone
    #     (TRANSP X starts at the first half-cell, never at 0)
    #   * rhot = 1  -> profile vanishes outside the LCFS
    if xs[0] > 0.0:
        xs = np.concatenate(([0.0], xs))
        ys = np.concatenate(([ys[0]], ys))
    if xs[-1] < 1.0:
        xs = np.concatenate((xs, [1.0]))
        ys = np.concatenate((ys, [0.0]))

    f = PchipInterpolator(xs, ys, extrapolate=False)
    out = f(rhot_on_z)
    out = np.where(np.isnan(out), 0.0, out)
    return out, rhot_on_z


# -----------------------------------------------------------------------------
# Physics: central vs outer contributions
# -----------------------------------------------------------------------------
def central_contribution(x_rhot, thntx_x, dvol_x, Rmag, Zmag, eq, time_jet):
    """Cumulative volume-integrated thermal-neutron rate inside the central
    flux surface that the LOS intersects at Z = Zmag + r0.

    Algorithm (as suggested by the user):
        1.  rhot_edge = rhot at (R = Rmag, Z = Zmag + r0)
        2.  cum_emis(rhot) = cumsum(THNTX * DVOL) along the TRANSP rhot grid
        3.  C0 = cum_emis evaluated (interpolated) at rhot_edge

    The result is the *toroidally integrated* emission rate inside the central
    flux surface, i.e. the full volume integral of THNTX up to rhot_edge.
    The cylindrical reference volume V0_ref = w * pi * r0^2 and the actual
    flux-surface volume V_inside are returned for context.

    Returns
    -------
    C0          : scalar [n/s]
    rhot_edge   : scalar - rhot at the edge of the central region
    V_inside    : scalar [m^3] - flux-surface volume inside rhot_edge
    V0_ref      : scalar [m^3] - cylindrical reference volume (w * pi * r0^2)
    cum_emis    : 1D array - cumulative emission along the rhot grid
    """
    # 1) rhot at the top edge of the central region (R=Rmag, Z=Zmag+r0)
    rhot_edge = float(RZ_to_rhot(
        np.array([Rmag]),
        np.array([Zmag + R0_CENTRAL]),
        eq, time_jet,
        method='cubic', clip_psin=True,
    )[0])

    xs = np.asarray(x_rhot)
    th = np.asarray(thntx_x)
    dv = np.asarray(dvol_x)

    # 2) cumulative flux-surface integrals on the TRANSP rhot grid
    cum_emis = np.cumsum(th * dv)
    cum_vol = np.cumsum(dv)

    # 3) interpolate to rhot_edge
    def _interp_at_edge(cum):
        if rhot_edge <= xs[0]:
            return float(cum[0])
        if rhot_edge >= xs[-1]:
            return float(cum[-1])
        return float(PchipInterpolator(xs, cum, extrapolate=False)(rhot_edge))

    C0 = _interp_at_edge(cum_emis)
    V_inside = _interp_at_edge(cum_vol)
    V0_ref = W_TOROIDAL * np.pi * R0_CENTRAL ** 2

    return C0, rhot_edge, V_inside, V0_ref, cum_emis


def equivalent_rhot(target, x_rhot, cum_emis):
    """Find rhot such that ``cumsum(THNTX*DVOL)`` up to that rhot equals
    *target*.

    Used as a consistency check: if the LOS-based estimate of the central +
    outer signal is interpreted as a closed-volume integral, this returns
    the rhot boundary that would reproduce it. If *target* exceeds the
    full-plasma integral (``cum_emis[-1]``) no such rhot exists inside the
    LCFS, and NaN is returned together with the (target / total) ratio so
    the caller can report the over-count.

    Returns
    -------
    rhot_eq      : float, rhot in [x_rhot[0], x_rhot[-1]] or NaN
    total_plasma : float, full-plasma cumulative emission (= cum_emis[-1])
    """
    total_plasma = float(cum_emis[-1])
    xs = np.asarray(x_rhot)
    cum = np.asarray(cum_emis)

    if target <= cum[0]:
        return float(xs[0]), total_plasma
    if target >= total_plasma:
        return float('nan'), total_plasma

    # Make the (cum, x) pairs strictly monotonic before inverting with PCHIP.
    keep = np.concatenate(([True], np.diff(cum) > 0))
    cum_m = cum[keep]
    xs_m = xs[keep]
    f_inv = PchipInterpolator(cum_m, xs_m, extrapolate=False)
    return float(f_inv(target)), total_plasma


def outer_contribution(thntx_z, dvol_z, Z, Zmag):
    """Compute Crest, the off-axis contribution to the LOS signal.

    For every Z with |Z - Zmag| > r0 the geometrical weight is

        w_geom(Z) = w / (2*pi*|Z - Zmag|)

    representing the fraction of the poloidal flux-surface layer at radius
    r = |Z - Zmag| that is actually intercepted by the LOS of toroidal
    width w. Integration is performed on the upper half (Z > Zmag + r0) and
    multiplied by 2 to account for up-down symmetry.

    Returns
    -------
    crest        : scalar [n/s]
    w_geom_full  : 1D array, same length as Z, zero inside the central cut
    integrand    : 1D array, THNTX(Z) * DVOL(Z) * w_geom(Z) over the full Z grid
    """
    dz = np.abs(Z - Zmag)
    outside_central = dz > R0_CENTRAL

    w_geom_full = np.zeros_like(Z, dtype=float)
    # safe division: zero outside the kept band, no 1/0 inside the cut
    w_geom_full[outside_central] = W_TOROIDAL / (2.0 * np.pi * dz[outside_central])

    integrand = thntx_z * dvol_z * w_geom_full

    upper_half = (Z - Zmag) > R0_CENTRAL
    crest = 2.0 * float(np.sum(integrand[upper_half]))

    return crest, w_geom_full, integrand


# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
def diagnostic_plots(Z, Zmag, thntx_z, w_geom, integrand,
                     pulse, runid, time_jet, eq_time):
    """Three-panel diagnostic: THNTX(Z), w_geom(Z) and weighted integrand."""
    dz = Z - Zmag

    fig, axes = plt.subplots(3, 1, figsize=(7.5, 9.0), sharex=True)
    fig.suptitle(
        f'Outer thermal-neutron fraction estimate\n'
        f'pulse {pulse}  TRANSP {runid}  '
        f't_JET={time_jet:.3f}s  t_EQ={eq_time:.3f}s'
    )

    axes[0].plot(dz, thntx_z, lw=1.4)
    axes[0].set_ylabel(r'THNTX$(Z)$ [n s$^{-1}$ m$^{-3}$]')
    axes[0].set_yscale('log')

    axes[1].plot(dz, w_geom, lw=1.4)
    axes[1].set_ylabel(r'$w_{\rm geom}(Z) = w/(2\pi |Z-Z_{\rm mag}|)$')

    axes[2].plot(dz, integrand, lw=1.4)
    axes[2].set_ylabel(r'THNTX$\,\cdot$DVOL$\,\cdot w_{\rm geom}$')
    axes[2].set_xlabel(r'$Z - Z_{\rm mag}$ [m]')

    for ax in axes:
        ax.axvline(R0_CENTRAL, color='r', ls='--', lw=0.8)
        ax.axvline(-R0_CENTRAL, color='r', ls='--', lw=0.8)
        ax.grid(True, ls=':', lw=0.5)
    axes[0].legend(['THNTX(Z)', f'|Z-Zmag|={R0_CENTRAL} m'], loc='best')

    plt.tight_layout()
    plt.show()


# -----------------------------------------------------------------------------
# Main entry point
# -----------------------------------------------------------------------------
def main(argv=None):
    args = parse_args(argv)

    # ------- TRANSP -----------------------------------------------------
    print(f'Loading TRANSP pulse {args.pulse} run {args.runid} ...')
    tr = ps.Transp(args.pulse, args.runid)

    t_jet_grid = tr.t + 40.0
    if args.time is None:
        time_jet = float(t_jet_grid[len(t_jet_grid) // 2])
        print(f'No time supplied; using midpoint t_JET = {time_jet:.3f} s.')
    else:
        time_jet = float(args.time)

    thntx_x, dvol_x, x_rhot, trind = get_transp_slice(tr, time_jet)

    # Pick up the native TRANSP units so the printout is unambiguous
    # (TRANSP is CGS: DVOL in CM**3, THNTX in N/CM3/SEC -> product is n/s).
    tr.add_data('THNTX', 'DVOL')
    unit_th = tr.units('THNTX') or ''
    unit_dv = tr.units('DVOL') or ''
    # Convert any DVOL-derived volume to m^3 for display only.
    dvol_to_m3 = 1.0e-6 if 'CM' in unit_dv.upper() else 1.0

    print(f'TRANSP slice [{trind}]: t_JET = {t_jet_grid[trind]:.3f} s '
          f'(t_TRANSP = {tr.t[trind]:.3f} s)')
    print(f'  THNTX units = [{unit_th}],  DVOL units = [{unit_dv}]')
    print(f'  THNTX(rhot=0) = {thntx_x[0]:.3e}  [{unit_th}]')
    print(f'  sum(DVOL)     = {dvol_x.sum() * dvol_to_m3:.3e}  '
          f'[m^3]  (total plasma volume, converted from {unit_dv})')

    # ------- Equilibrium ------------------------------------------------
    print(f'Loading equilibrium PPF {args.dda}/{args.uid}/{args.seq} ...')
    eq = ps.Eq(args.pulse, dda=args.dda, uid=args.uid, seq=args.seq)
    tind_eq = int(np.abs(eq.t - time_jet).argmin())
    print(f'Equilibrium slice [{tind_eq}]: t = {eq.t[tind_eq]:.3f} s')

    # ------- LOS geometry ----------------------------------------------
    R, Z, Rmag, Zmag = build_los(eq, tind_eq, args.nz)
    print(f'  Rmag = {Rmag:.3f} m, Zmag = {Zmag:.3f} m')
    print(f'  LOS: R = {Rmag:.3f} m, Z in [{Z.min():.3f}, {Z.max():.3f}] m '
          f'({args.nz} samples)')

    # ------- Remap profiles onto the LOS Z grid -------------------------
    # (used for the outer integral and for the diagnostic plots; the
    # central contribution is computed directly on the TRANSP rhot grid.)
    thntx_z, _ = remap_profile_on_z(x_rhot, thntx_x, R, Z, eq, time_jet)
    dvol_z, _ = remap_profile_on_z(x_rhot, dvol_x, R, Z, eq, time_jet)

    iz0 = int(np.argmin(np.abs(Z - Zmag)))
    th0 = float(thntx_z[iz0])

    # ------- Central / outer contributions ------------------------------
    C0, rhot_edge, V_inside, V0_ref, cum_emis = central_contribution(
        x_rhot, thntx_x, dvol_x, Rmag, Zmag, eq, time_jet)
    Crest, w_geom, integrand = outer_contribution(thntx_z, dvol_z, Z, Zmag)

    total = C0 + Crest
    f_outer = (Crest / total) if total > 0 else float('nan')

    # ------- Report -----------------------------------------------------
    print()
    print('-------------------- Central region --------------------')
    print(f'  r0                  = {R0_CENTRAL:.3f} m')
    print(f'  w (toroidal width)  = {W_TOROIDAL:.3f} m')
    print(f'  (R, Z) at edge      = ({Rmag:.3f}, {Zmag + R0_CENTRAL:.3f}) m')
    print(f'  rhot at edge        = {rhot_edge:.4f}')
    print(f'  V_inside (flux vol) = {V_inside * dvol_to_m3:.4e} m^3')
    print(f'  V0_ref = w*pi*r0^2  = {V0_ref:.4e} m^3  (reference cylinder)')
    print(f'  TH0 (THNTX at axis) = {th0:.4e} [{unit_th}]')
    print(f'  C0 = cumsum(THNTX*DVOL) up to rhot_edge = {C0:.4e} n/s')
    print()
    print('--------------------- Outer region ---------------------')
    print(f'  |Z - Zmag| > {R0_CENTRAL:.3f} m,  '
          f'w_geom = w/(2*pi*|Z-Zmag|)')
    print(f'  Crest               = {Crest:.4e} n/s')
    print()
    print('---------------------- Summary -------------------------')
    print(f'  C0 + Crest          = {total:.4e} n/s')
    print(f'  f_outer             = {f_outer:.4f}')
    print(f'  outer fraction      = {100.0 * f_outer:.2f} %')

    # ------- Equivalent rhot consistency check --------------------------
    rhot_eq, total_plasma = equivalent_rhot(total, x_rhot, cum_emis)
    print()
    print('------------ Equivalent rhot_bnd check -----------------')
    print(f'  total plasma cumsum(THNTX*DVOL) = {total_plasma:.4e} n/s')
    if np.isnan(rhot_eq):
        ratio = total / total_plasma if total_plasma > 0 else float('nan')
        print(f'  C0 + Crest exceeds the full-plasma integral '
              f'(ratio = {ratio:.4f}),')
        print(f'  so no rhot_bnd inside the LCFS reproduces it.')
    else:
        print(f'  rhot_bnd such that cumsum up to it = C0 + Crest: '
              f'{rhot_eq:.4f}')

    if args.plot:
        diagnostic_plots(Z, Zmag, thntx_z, w_geom, integrand,
                         args.pulse, args.runid, time_jet, eq.t[tind_eq])

    return 0


if __name__ == '__main__':
    sys.exit(main())
