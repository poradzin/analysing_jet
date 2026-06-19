"""los_common.py -- shared LOS / equilibrium / binning library.

Hosts the symbols used by the KM14 and KM9 line-of-sight thermal-rate scripts,
plus the per-cell detector-rate calculation that works for any LOS geometry.

Migrated from src/neutron/km14/{los_thermal_rate,los_th_bt_ratio,
bt_zone_integrator}.py on 2026-06-18 so a KM9 sibling script can share the
implementation without code duplication. The km14 scripts now import from here
too; see src/neutron/km14/los_thermal_rate.py for the box-chord pipeline and
plots that remain KM14-specific.

Public surface
~~~~~~~~~~~~~~
I/O
    find_run_dir, read_time_grid, read_thermal_slice, read_los_file
    read_scalar_totals                            -- 0D BTNTS_* totals
    read_fbm_avg_window, resolve_time_selection   -- -t / --idx time handling

Beam-target / flux profiles (los_th_bt_ratio scripts)
    flux_avg_profile, interp_flux_to

Equilibrium (rhot(R, Z) from a TRANSP CDF or a JET PPF)
    CdfEquilibrium  -- low-level psi(R,Z) -> rhot from the main CDF
    EqCDF           -- high-level wrapper used by the box pipelines
    EqPPF          -- same interface backed by a PPF (lazy ppf import)
    _boundary_from_moments, _rhot_scatter, _lcfs_contour  -- helpers

Numerics
    thntx_on_grid, find_rho_bnd, _subgrid_bin, _running_median3

LOS rates
    los_shell_fraction          -- f(rhot) from a box-chord (R, Z) grid
    los_file_detector_rate      -- f_det(rhot) + Rate_det from real LOS cells
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
from netCDF4 import Dataset
from matplotlib.path import Path as MplPath
from scipy.interpolate import PchipInterpolator


# ----------------------------------------------------------------------
# Data tree
# ----------------------------------------------------------------------

DEFAULT_LOCAL_BASE = Path.home() / "jet" / "data"
HEIMDALL_BASE = Path("/common/transp_shared/Data/result/JET")


def find_run_dir(pulse, run_suffix, data_dir=None):
    """Locate a TRANSP run directory <base>/<pulse>/<run_suffix>.

    Search order: explicit ``data_dir`` (if given), the local WSL tree under
    ~/jet/data, then the heimdall results tree. Raises ``FileNotFoundError``
    if none of those exist.
    """
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


# ----------------------------------------------------------------------
# Thermal profiles from the TRANSP main CDF
# ----------------------------------------------------------------------

def read_time_grid(cdf_path):
    """Return TIME3 (TRANSP seconds) from the main CDF."""
    with Dataset(cdf_path, "r") as d:
        return np.array(d["TIME3"][:])


def read_thermal_slice(cdf_path, th_var, trind):
    """Return (thntx, x_rhot, dvol, thntx_to_si, dvol_to_m3, units) at *trind*.

    Reads the chosen thermal emissivity profile (``th_var``: THNTX, THNTX_DD,
    THNTX_DT, ...), the rhot grid X, and DVOL. Returns the CGS->SI conversion
    factors inferred from the CDF units attributes.

    ``trind`` may be a single integer (one time slice) **or a sequence of
    integers**, in which case the profiles ``THNTX``, ``X`` and ``DVOL`` are
    **averaged element-wise across the requested time slices** (a straight,
    unweighted mean over the TIME3 indices). The zone grids X and DVOL are
    nearly time-stationary so averaging them is well defined; this implements
    the time-window averaging requested by the ``-t t1 t2`` / ``--idx`` CLI
    modes (see :func:`resolve_time_selection`).

    NB: THNTX is a normal TIME3 output -- it is **not** in the run's ``SELAVG``
    list, so TRANSP never applied its own ACFILE/``MTHDAVG`` averaging to it
    (that governed only the fast-ion data ``FBM BMVOL ...`` written to
    ``_fi``/``_neut``). The windowed mean here is our own construct, intended to
    make the thermal profile time-consistent with the ``[OUTTIM-AVGTIM,
    OUTTIM]`` window the fast-ion / BT data was averaged over. A plain TIME3
    mean is the best achievable (THNTX is not stored on the heating-source
    timestep grid ``MTHDAVG=2`` sampled); on a uniform TIME3 grid it equals the
    trapezoidal time-average anyway.
    """
    trinds = np.atleast_1d(np.asarray(trind, dtype=int))
    with Dataset(cdf_path, "r") as d:
        if th_var not in d.variables:
            raise KeyError(f"{th_var} not present in {cdf_path}")
        thntx = np.asarray(d[th_var][trinds]).mean(axis=0)
        x_rhot = np.asarray(d["X"][trinds]).mean(axis=0)
        dvol = np.asarray(d["DVOL"][trinds]).mean(axis=0)
        unit_th = getattr(d[th_var], "units", "") or ""
        unit_dv = getattr(d["DVOL"], "units", "") or ""
    thntx_to_si = 1.0e6 if 'CM' in unit_th.upper() else 1.0   # N/CM3/SEC -> n/m3/s
    dvol_to_m3 = 1.0e-6 if 'CM' in unit_dv.upper() else 1.0   # CM**3 -> m^3
    return thntx, x_rhot, dvol, thntx_to_si, dvol_to_m3, unit_th


# ----------------------------------------------------------------------
# Beam-target / flux-profile helpers (shared by the los_th_bt_ratio scripts)
# ----------------------------------------------------------------------

def read_scalar_totals(cdf_path, tind):
    """0D beam-target neutron totals BTNTS_{DD,DT,TT,TD} [n/s] at time index tind."""
    with Dataset(cdf_path, "r") as d:
        out = {}
        for k in ("BTNTS_DD", "BTNTS_DT", "BTNTS_TT", "BTNTS_TD"):
            if k in d.variables:
                out[k] = float(d[k][tind])
    return out


def flux_avg_profile(values_zone, x2d, bmvol):
    """BMVOL-weighted flux-surface average of a per-zone NUBEAM quantity.

    Returns ``(x_rows, profile)`` -- the unique x-row labels and the
    volume-weighted average of ``values_zone`` over the poloidal zones in each
    x-row (collapsing the NUBEAM (x-row, theta) zone grid to a 1-D flux profile).
    """
    xu = np.unique(x2d)
    prof = np.array([np.sum(values_zone[x2d == xv] * bmvol[x2d == xv])
                     / np.sum(bmvol[x2d == xv]) for xv in xu])
    return xu, prof


def interp_flux_to(x_src, y_src, x_dst):
    """PCHIP a flux profile y_src(x_src) onto x_dst, padded to cover [0, 1].

    rhot=0 takes the on-axis value, rhot=1 is forced to zero, points outside the
    source support are set to 0 (never negative). Used to put a beam-target flux
    profile on the same rhot grid as THNTX/DVOL.
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


def read_fbm_avg_window(tr_dat_path, idx):
    """Fast-ion output averaging window for output index *idx* from the namelist.

    TRANSP writes ACFILE output (``<run>.DATA<idx>`` -> ``_fi_<idx>.cdf`` /
    ``_neut_<idx>.cdf``) at the times ``OUTTIM(idx)`` listed in the run namelist
    ``<pulse><runid>TR.DAT`` (group ``&ACFILE``). The variables named in
    ``SELAVG`` (the fast-ion set ``FBM BMVOL BDENS2 EBA2PL EBA2PP``) are
    averaged over the preceding ``AVGTIM`` seconds, i.e. the window is

        [OUTTIM(idx) - AVGTIM,  OUTTIM(idx)]

    (all in TRANSP seconds = JET time - 40). ``idx`` is 1-based to match the
    Fortran ``OUTTIM(idx)`` indexing and the ``.DATA<idx>``/``_fi_<idx>.cdf``
    file numbering. ``MTHDAVG`` selects how TRANSP sampled within that window
    (1=heating-source timestep [default], 2=on each heating-source timestep
    completion, 3=fixed ``AVGSAMP`` period); it governs the **fast-ion** average
    only -- THNTX (thermal) is not in ``SELAVG`` and is unaffected (see
    :func:`read_thermal_slice`). We reuse this window purely to average THNTX
    consistently with the fast-ion snapshot.

    Returns ``(t_lo, t_hi, t_out, avgtim, mthdavg)`` in TRANSP seconds, where
    ``t_out = OUTTIM(idx)`` is the nominal output time (window upper bound).
    """
    import f90nml
    nml = f90nml.read(str(tr_dat_path))
    if 'acfile' not in nml:
        raise KeyError(f"no &ACFILE group in {tr_dat_path}")
    ac = nml['acfile']
    outtim = ac.get('outtim')
    avgtim = ac.get('avgtim')
    if outtim is None or avgtim is None:
        raise KeyError(f"OUTTIM/AVGTIM missing from &ACFILE in {tr_dat_path}")
    outtim = np.atleast_1d(np.asarray(outtim, dtype=float))
    if idx < 1 or idx > outtim.size:
        raise IndexError(
            f"--idx {idx} out of range: OUTTIM has {outtim.size} entries "
            f"({outtim.tolist()})")
    t_out = float(outtim[idx - 1])
    avgtim = float(avgtim)
    mthdavg = ac.get('mthdavg')
    return t_out - avgtim, t_out, t_out, avgtim, mthdavg


def resolve_time_selection(cdf_path, run_dir, run_id, times=None, idx=None):
    """Resolve the ``-t`` / ``--idx`` CLI selection into slice indices + eq time.

    Selection modes (mutually exclusive ``times`` vs ``idx``):

    * ``idx`` (1-based) -- read ``[OUTTIM(idx)-AVGTIM, OUTTIM(idx)]`` from the
      run namelist; average all TIME3 slices in that window; map the
      equilibrium to the slice nearest ``OUTTIM(idx)`` (the window upper bound,
      the nominal fast-ion output time).
    * ``times = [t1, t2]`` (JET s) -- average all TIME3 slices in ``[t1, t2]``;
      map the equilibrium to the **window midpoint** ``(t1+t2)/2``.
    * ``times = [t1]`` -- snap to the single nearest TIME3 slice; equilibrium at
      that slice.
    * neither -- use the midpoint slice of the run.

    All times are handled in the JET convention (TRANSP seconds + 40). If a
    window encloses no TIME3 slice (narrower than the grid spacing) it falls
    back to the single nearest slice to the window centre, with a warning.

    Returns ``(trinds, eq_ref_jet, rep_jet, desc)``:
        trinds      1-D int array of TIME3 indices to average over
        eq_ref_jet  JET time the (single) equilibrium maps to [s]
        rep_jet     representative time of the selection (filenames/titles) [s]
        desc        human-readable description of the resolved selection
    """
    t_jet_grid = read_time_grid(cdf_path) + 40.0

    def _nearest(t):
        return int(np.abs(t_jet_grid - t).argmin())

    def _window(lo_jet, hi_jet, centre_jet):
        sel = np.where((t_jet_grid >= lo_jet) & (t_jet_grid <= hi_jet))[0]
        if sel.size == 0:
            j = _nearest(centre_jet)
            print(f'  WARNING: window [{lo_jet:.4f}, {hi_jet:.4f}] s encloses '
                  f'no TIME3 slice (dt~{np.median(np.diff(t_jet_grid)):.4f} s); '
                  f'falling back to nearest slice [{j}] @ {t_jet_grid[j]:.4f} s.')
            sel = np.array([j], dtype=int)
        return sel

    if idx is not None:
        if times:
            raise ValueError("--idx and -t are mutually exclusive")
        tr_dat = Path(run_dir) / f"{run_id}TR.DAT"
        if not tr_dat.exists():
            raise FileNotFoundError(f"namelist not found: {tr_dat}")
        t_lo, t_hi, t_out, avgtim, mthdavg = read_fbm_avg_window(tr_dat, idx)
        lo_jet, hi_jet, out_jet = t_lo + 40.0, t_hi + 40.0, t_out + 40.0
        trinds = _window(lo_jet, hi_jet, out_jet)
        eq_ref_jet = float(t_jet_grid[_nearest(out_jet)])
        desc = (
            f'--idx {idx}: OUTTIM={t_out:.4f}s, AVGTIM={avgtim:.4f}s '
            f'(fast-ion MTHDAVG={mthdavg}) -> TRANSP window '
            f'[{t_lo:.4f}, {t_hi:.4f}]s (JET [{lo_jet:.4f}, {hi_jet:.4f}]s); '
            f'mean of {trinds.size} TIME3 THNTX slice(s) over the same window; '
            f'eq @ OUTTIM (JET {out_jet:.4f}s)'
        )
        return trinds, eq_ref_jet, out_jet, desc

    if times and len(times) == 2:
        t1, t2 = sorted(float(t) for t in times)
        trinds = _window(t1, t2, 0.5 * (t1 + t2))
        eq_ref_jet = float(t_jet_grid[_nearest(0.5 * (t1 + t2))])
        rep = 0.5 * (t1 + t2)
        desc = (
            f'-t {t1:.4f} {t2:.4f} s: averaging {trinds.size} TIME3 slice(s) '
            f'in [{t1:.4f}, {t2:.4f}]s; eq @ window midpoint (JET {rep:.4f}s)'
        )
        return trinds, eq_ref_jet, rep, desc

    if times and len(times) == 1:
        t1 = float(times[0])
        j = _nearest(t1)
        desc = (f'-t {t1:.4f} s: nearest TIME3 slice [{j}] @ '
                f'{t_jet_grid[j]:.4f}s; eq @ same slice')
        return np.array([j], dtype=int), float(t_jet_grid[j]), float(t_jet_grid[j]), desc

    # No time supplied: midpoint of the run.
    j = len(t_jet_grid) // 2
    desc = (f'no time supplied: run-midpoint TIME3 slice [{j}] @ '
            f'{t_jet_grid[j]:.4f}s')
    return np.array([j], dtype=int), float(t_jet_grid[j]), float(t_jet_grid[j]), desc


# ----------------------------------------------------------------------
# Real line-of-sight cell file (_KM3.los, _KM9.los, ...)
# ----------------------------------------------------------------------

def read_los_file(path):
    """Read a Genesis line-of-sight cell file.

    Each row is one LOS cell with 8 whitespace-separated columns

        x, y, z, C, V, u, v, w

    (see ``KM3_LoS_readme.txt`` / ``KM9_LoS_readme.txt``). All in a fixed,
    discharge-/time-independent right-handed Cartesian frame: the poloidal
    plane is x-z, the toroidal plane is x-y, so

      * x       -- tangential / toroidal coordinate [m]
      * y       -- horizontal in-plane coordinate [m]; cell major radius is
                   R = sqrt(x^2 + y^2)
      * z       -- vertical coordinate Z [m]
      * C       -- etendue weight [m^3]. Detector solid angle seen from the
                   cell is Omega = 4*pi*C/V, so for an *isotropic* emitter of
                   volumetric emissivity eps the rate of neutrons that
                   actually reach the detector is eps * V * Omega/(4*pi) =
                   eps * C.
      * V       -- cell volume [m^3]
      * u, v, w -- unit vector of the emission direction a neutron born in
                   the cell must have to reach the detector.

    Returns a dict with 1-D arrays ``x, y, z, C, V, u, v, w`` plus derived
    ``R = sqrt(x^2 + y^2)`` and ``Z = z``. (For an almost-vertical chord like
    KM14, |x| << y so R ~= y; for a horizontal chord like KM9, x dominates
    and the chord wraps the torus.)
    """
    raw = np.loadtxt(path)
    if raw.ndim != 2 or raw.shape[1] != 8:
        raise ValueError(
            f"{path}: expected an (N, 8) table (x,y,z,C,V,u,v,w), "
            f"got shape {raw.shape}"
        )
    x, y, z, C, V, u, v, w = raw.T
    return {
        "x": x, "y": y, "z": z, "C": C, "V": V, "u": u, "v": v, "w": w,
        "R": np.hypot(x, y), "Z": z.copy(),
    }


# ----------------------------------------------------------------------
# Equilibrium: psi(R, Z) -> rhot from a TRANSP CDF (no ppf needed)
# ----------------------------------------------------------------------

def _boundary_from_moments(d, tind, ntheta=257):
    """Reconstruct the LCFS (R, Z) [m] from TRANSP asymmetric boundary moments.

    TRANSP stores the plasma boundary as a Fourier series in a poloidal angle:

        R(theta) = RMCB0 + sum_{n>=1} [RMCBn cos(n th) + RMSBn sin(n th)]
        Z(theta) = YMCB0 + sum_{n>=1} [YMCBn cos(n th) + YMSBn sin(n th)]

    (units cm in the CDF). Verified to reproduce the ``_fi`` RSURF/ZSURF
    outermost surface to ~1e-4 cm on 104614 M30.
    """
    th = np.linspace(0.0, 2.0 * np.pi, ntheta)
    R = np.zeros_like(th)
    Z = np.zeros_like(th)
    n = 0
    while f"RMCB{n}" in d.variables:
        R += float(d[f"RMCB{n}"][tind]) * np.cos(n * th)
        Z += float(d[f"YMCB{n}"][tind]) * np.cos(n * th)
        if n >= 1 and f"RMSB{n}" in d.variables:
            R += float(d[f"RMSB{n}"][tind]) * np.sin(n * th)
            Z += float(d[f"YMSB{n}"][tind]) * np.sin(n * th)
        n += 1
    return R / 100.0, Z / 100.0


class CdfEquilibrium:
    """rhot(R, Z) from the TRANSP main CDF at a chosen time.

    Low-level: holds the psi_n(R, Z) grid, the magnetic axis, and the
    psi_n -> rhot map (via PLFLX vs XB). Used by ``los_th_bt_ratio.py``'s
    own pipeline directly and by :class:`EqCDF` (which wraps it and adds the
    high-level ``rhot_on_grid`` / ``rhot_at_points`` interface).

    Migrated verbatim from ``los_th_bt_ratio.CdfEquilibrium``; the methods
    ``rhot`` (bilinear, kept for ad-hoc lookups) and ``rhot_pinned`` (cubic
    griddata with pinned axis + optional pinned LCFS) are unchanged.
    """

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
        psi_n = 1 and removes that crater.
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


def _rhot_scatter(Rq, Zq, pts_R, pts_Z, pts_psin, lcfs_RZ, psin_to_rhot):
    """Evaluate (rhot, inside, psin) at scattered query points (Rq, Zq).

    Shared core of :class:`EqCDF` / :class:`EqPPF` -- griddata-cubic of the
    psi_n node set (PSIRZ grid + pinned axis at psi_n=0 + pinned LCFS at
    psi_n=1), linear fallback for NaNs, point-in-LCFS-polygon mask, then
    psi_n -> rhot via the source-specific ``psin_to_rhot`` callable. Works
    for any-shape flat inputs, so it serves both the dense (R, Z) box grid
    (via ``rhot_on_grid``) and the irregular cloud of real LOS cells (via
    ``rhot_at_points``).
    """
    from scipy.interpolate import griddata
    Rq = np.asarray(Rq, dtype=float).ravel()
    Zq = np.asarray(Zq, dtype=float).ravel()
    query = np.column_stack([Rq, Zq])
    psin = griddata((pts_R, pts_Z), pts_psin, query, method='cubic')
    if np.any(np.isnan(psin)):
        psin_lin = griddata((pts_R, pts_Z), pts_psin, query, method='linear')
        psin = np.where(np.isnan(psin), psin_lin, psin)
    inside = MplPath(lcfs_RZ).contains_points(query)
    psin_clip = np.clip(np.where(np.isnan(psin), 1.0, psin), 0.0, 1.0)
    rhot = np.clip(np.asarray(psin_to_rhot(psin_clip)), 0.0, 1.0)
    return rhot, inside, psin


def _lcfs_contour(eq, tind_eq):
    """Return (Rbnd, Zbnd) at tind_eq using the correct time-first indexing.

    profiles.py reshapes RBND/ZBND with ``(-1, n_t)`` but the array actually
    ends up shaped ``(n_t, n_bnd)`` (confirmed empirically: on pulse 104614
    only ``_Zbnd[tind, :]`` gives the expected ~[-1.37, 1.68] range, while
    ``_Zbnd[:, tind]`` returns a single bnd point sampled across the
    discharge). The LCFS may be padded with zeros at unused slots, so we
    drop the trailing (0, 0) entries.
    """
    Rb = np.asarray(eq._Rbnd[tind_eq, :], dtype=float)
    Zb = np.asarray(eq._Zbnd[tind_eq, :], dtype=float)
    valid = ~((Rb == 0.0) & (Zb == 0.0)) & ~np.isnan(Rb) & ~np.isnan(Zb)
    return Rb[valid], Zb[valid]


class EqCDF:
    """Self-contained equilibrium from the TRANSP main CDF (no ppf).

    Common interface (used by both EqCDF and EqPPF):
        .tind        thermal/eq time index on the relevant time grid (int)
        .t_eq_jet    equilibrium time in JET convention [s]
        .label       human-readable source description
        .Rmag, .Zmag magnetic axis [m]
        .z_extent()  -> (z_bot, z_top) computational-box Z range [m]
        .lcfs()      -> (Rb, Zb) closed LCFS polygon [m]
        .rhot_on_grid(R, Z) -> (rhot, inside, psin_dense) on the (nR, nZ) grid
        .rhot_at_points(Rq, Zq) -> (rhot, inside, psin) at scattered points
        .native_rhot()      -> (Rg_n, Zg_n, rhot_native) masked to the LCFS
    """

    label = "TRANSP CDF (PSIRZ + boundary moments)"

    def __init__(self, cdf_path, time_jet):
        self.ceq = CdfEquilibrium(cdf_path, time_jet - 40.0)  # TRANSP time
        self.tind = self.ceq.tind
        self.t_eq_jet = self.ceq.t_used + 40.0
        self.RG = self.ceq.RG
        self.ZG = self.ceq.ZG
        with Dataset(cdf_path, "r") as d:
            self.Rmag = float(d["RAXIS"][self.tind]) / 100.0
            self.Zmag = float(d["YAXIS"][self.tind]) / 100.0
            self._Rb, self._Zb = _boundary_from_moments(d, self.tind)

    def z_extent(self):
        return float(self.ZG.min()), float(self.ZG.max())

    def lcfs(self):
        return self._Rb, self._Zb

    def _scatter_nodes(self):
        """psi_n node set for griddata: PSIRZ grid + pinned axis + pinned LCFS.

        The magnetic axis is pinned at psi_n = 0 because plain interpolation of
        the coarse PSIRZ grid (dR ~ 2 cm) cannot recover the true paraboloidal
        axis minimum (it floors psi_n ~ 3e-4, i.e. rhot ~ 0.02, zeroing the
        innermost f(rhot) bins). The LCFS polygon is pinned at psi_n = 1 because
        the cubic interpolation otherwise floors short of 1 inside the mask
        (max psi_n ~ 0.999), starving the outermost flux shell of LOS cells and
        cratering f(rhot) in the last 1-2 bins.
        """
        Rg_n, Zg_n = np.meshgrid(self.RG, self.ZG)  # (nZ, nR), matches psin[iz,ir]
        Rb, Zb = self.lcfs()
        pts_R = np.concatenate([Rg_n.ravel(), [self.Rmag], Rb])
        pts_Z = np.concatenate([Zg_n.ravel(), [self.Zmag], Zb])
        pts_psin = np.concatenate([self.ceq.psin.ravel(), [0.0],
                                   np.ones(Rb.size)])
        return pts_R, pts_Z, pts_psin

    def rhot_at_points(self, Rq, Zq):
        """(rhot, inside, psin) at scattered query points (flat arrays)."""
        pts_R, pts_Z, pts_psin = self._scatter_nodes()
        Rb, Zb = self.lcfs()
        return _rhot_scatter(
            Rq, Zq, pts_R, pts_Z, pts_psin, np.column_stack([Rb, Zb]),
            lambda p: np.interp(p, self.ceq._psin_b, self.ceq._rhot_b),
        )

    def rhot_on_grid(self, R, Z):
        Rg, Zg = np.meshgrid(R, Z, indexing='ij')
        rhot, inside, psin = self.rhot_at_points(Rg.ravel(), Zg.ravel())
        return (rhot.reshape(Rg.shape), inside.reshape(Rg.shape),
                psin.reshape(Rg.shape))

    def native_rhot(self):
        Rg_n, Zg_n = np.meshgrid(self.RG, self.ZG)  # (nZ, nR), matches psin[iz,ir]
        psin_clip = np.clip(self.ceq.psin, 0.0, 1.0)
        rhot = np.clip(np.interp(psin_clip.ravel(), self.ceq._psin_b,
                                 self.ceq._rhot_b), 0.0, 1.0).reshape(Rg_n.shape)
        mask = MplPath(np.column_stack(self.lcfs())).contains_points(
            np.column_stack([Rg_n.ravel(), Zg_n.ravel()])).reshape(Rg_n.shape)
        return Rg_n, Zg_n, np.where(mask, rhot, np.nan)


class EqPPF:
    """JET PPF equilibrium via profiles.Eq (imports ppf lazily).

    Same interface as :class:`EqCDF`; used by :func:`los_thermal_rate.main`
    when ``--eq-source ppf`` is selected. The ppf import is deferred inside
    ``__init__`` so the cdf-only default path doesn't need ppf at all.
    """

    label = "PPF equilibrium"

    def __init__(self, pulse, dda, uid, seq, time_jet):
        # profiles.py / change_rho.py live in src/ (two levels up from
        # src/neutron/common/), and the calling script's own dir is the only
        # one on sys.path. Add src/ before the lazy import.
        src_dir = os.path.abspath(
            os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
        if src_dir not in sys.path:
            sys.path.append(src_dir)
        import profiles as ps
        from change_rho import psin_to_sqrt_ftor_norm
        self._psin_to_rhot = psin_to_sqrt_ftor_norm
        self.eq = ps.Eq(pulse, dda=dda, uid=uid, seq=seq)
        self.time_jet = time_jet
        self.tind = int(np.abs(self.eq.t - time_jet).argmin())
        self.t_eq_jet = float(self.eq.t[self.tind])
        self.Rmag = float(self.eq._Rmag[self.tind])
        self.Zmag = float(self.eq._Zmag[self.tind])
        self.label = f"PPF equilibrium {dda}/{uid}/{seq}"

    def z_extent(self):
        PSIZ = np.asarray(self.eq._psiz)
        return float(np.nanmin(PSIZ)), float(np.nanmax(PSIZ))

    def lcfs(self):
        return _lcfs_contour(self.eq, self.tind)

    def _scatter_nodes(self):
        """psi_n node set: PSIRZ grid + pinned axis (psi_n=0) + pinned LCFS (1).

        Same pinning rationale as :meth:`EqCDF._scatter_nodes` (axis floor and
        edge crater), so both equilibrium sources share one scattered
        evaluator.
        """
        eq, tind_eq = self.eq, self.tind
        Rg_n = np.asarray(eq._psirzmg[0])
        Zg_n = np.asarray(eq._psirzmg[1])
        psin_flat = np.asarray(eq._psi_norm[tind_eq]).ravel()
        Rb, Zb = self.lcfs()
        pts_R = np.concatenate([Rg_n.ravel(), [self.Rmag], Rb])
        pts_Z = np.concatenate([Zg_n.ravel(), [self.Zmag], Zb])
        pts_psin = np.concatenate([psin_flat, [0.0], np.ones(Rb.size)])
        return pts_R, pts_Z, pts_psin

    def rhot_at_points(self, Rq, Zq):
        """(rhot, inside, psin) at scattered query points (flat arrays)."""
        pts_R, pts_Z, pts_psin = self._scatter_nodes()
        Rb, Zb = self.lcfs()
        return _rhot_scatter(
            Rq, Zq, pts_R, pts_Z, pts_psin, np.column_stack([Rb, Zb]),
            lambda p: self._psin_to_rhot(p, self.eq, self.time_jet),
        )

    def rhot_on_grid(self, R, Z):
        Rg, Zg = np.meshgrid(R, Z, indexing='ij')
        rhot, inside, psin = self.rhot_at_points(Rg.ravel(), Zg.ravel())
        return (rhot.reshape(Rg.shape), inside.reshape(Rg.shape),
                psin.reshape(Rg.shape))

    def native_rhot(self):
        eq, tind_eq = self.eq, self.tind
        Rg_n = np.asarray(eq._psirzmg[0])
        Zg_n = np.asarray(eq._psirzmg[1])
        psin_native = np.asarray(eq._psi_norm[tind_eq]).reshape(Rg_n.shape)
        psin_clip = np.clip(np.where(np.isnan(psin_native), 1.0, psin_native),
                            0.0, 1.0)
        rhot = self._psin_to_rhot(
            psin_clip.ravel(), eq, self.time_jet).reshape(Rg_n.shape)
        Rb, Zb = self.lcfs()
        mask = MplPath(np.column_stack([Rb, Zb])).contains_points(
            np.column_stack([Rg_n.ravel(), Zg_n.ravel()])).reshape(Rg_n.shape)
        return Rg_n, Zg_n, np.where(mask, rhot, np.nan)


# ----------------------------------------------------------------------
# Numerics: profile interpolation, anti-aliased binning, rho_bnd inversion
# ----------------------------------------------------------------------

def thntx_on_grid(x_rhot, thntx_x, rhot_grid, inside_mask):
    """PCHIP-interpolate THNTX(X) onto rhot(R, Z); zero outside the LCFS.

    The TRANSP rhot grid is padded with the axis value at rhot=0 and 0 at
    rhot=1 so the interpolator covers the full [0, 1] range. Works on any
    shape of ``rhot_grid`` / ``inside_mask`` (dense 2-D for the box, flat
    cloud for the cell file).
    """
    order = np.argsort(x_rhot)
    xs = np.asarray(x_rhot)[order]
    ys = np.asarray(thntx_x)[order]
    uniq = np.concatenate(([True], np.diff(xs) > 0))
    xs = xs[uniq]
    ys = ys[uniq]

    # Pad so the interpolator covers [0, 1]:
    if xs[0] > 0.0:
        xs = np.concatenate(([0.0], xs))
        ys = np.concatenate(([ys[0]], ys))
    if xs[-1] < 1.0:
        xs = np.concatenate((xs, [1.0]))
        ys = np.concatenate((ys, [0.0]))

    f = PchipInterpolator(xs, ys, extrapolate=False)
    out = f(rhot_grid)
    out = np.where(np.isnan(out), 0.0, out)
    out = np.where(inside_mask, out, 0.0)
    return out


def find_rho_bnd(target, xs, cum_emis):
    """Find rho_bnd such that cumsum(THNTX*DVOL)|_{rho_bnd} = target.

    Returns (rho_bnd, total_plasma). rho_bnd is NaN if target exceeds the
    full-plasma cumulative emission cum_emis[-1].
    """
    total = float(cum_emis[-1])
    if target <= cum_emis[0]:
        return float(xs[0]), total
    if target >= total:
        return float('nan'), total
    keep = np.concatenate(([True], np.diff(cum_emis) > 0))
    f_inv = PchipInterpolator(cum_emis[keep], xs[keep], extrapolate=False)
    return float(f_inv(target)), total


def normalized_coupling_weight(weight, dvol_m3, xs):
    """Normalize a per-shell LOS coupling weight ``w(rhot)`` for cross-diagnostic
    comparison: divide by ``integral(w * DVOL d rhot)`` so that
    ``integral(w_norm * DVOL d rhot) = 1``.

    The detector-coupling weight ``f_det = C_bin / DVOL`` (and the geometric box
    ``f``) carry an arbitrary *absolute* scale -- the detector etendue / solid
    angle, which differs between KM9 and KM14 -- so their raw magnitudes are not
    comparable across the two diagnostics. Every LOS-weighted rate *ratio*
    (notably ``TH/BT = sum(TH*w*DVOL)/sum(BT*w*DVOL)``) is invariant under a
    constant rescaling of ``w``, since the factor cancels between numerator and
    denominator, so normalizing is free: it pins the scale to a common
    convention without changing any reported ratio. Units of ``w_norm``: 1/m^3.
    """
    w = np.asarray(weight, dtype=float)
    dv = np.asarray(dvol_m3, dtype=float)
    norm = float(np.trapezoid(w * dv, np.asarray(xs, dtype=float)))
    return w / norm if norm > 0 else w


def _running_median3(a):
    """3-point running median; endpoints unchanged. Cosmetic de-speckle."""
    a = np.asarray(a, dtype=float)
    if a.size < 3:
        return a.copy()
    out = a.copy()
    out[1:-1] = np.median(np.stack([a[:-2], a[1:-1], a[2:]]), axis=0)
    return out


def _running_median(a, k=5):
    """Centered k-point running median (k forced odd); cosmetic de-speckle.

    Near the ends the window shrinks symmetrically so the first/last bins stay
    put (no edge bias). For ``k == 3`` this matches :func:`_running_median3`.
    Used to suppress the residual flux-shell aliasing wiggle in a horizontal
    chord's ``f_det`` (a regular cell lattice beating against the rhot bins);
    cosmetic only -- ``Rate_det`` is summed from the cells, not this profile.
    """
    a = np.asarray(a, dtype=float)
    n = a.size
    if k % 2 == 0:
        k += 1
    if n < 3 or k < 3:
        return a.copy()
    h = k // 2
    out = a.copy()
    for i in range(n):
        w = min(h, i, n - 1 - i)          # symmetric, shrinking at the ends
        if w > 0:
            out[i] = np.median(a[i - w:i + w + 1])
    return out


def _subgrid_bin(rhot_c, half, weight, edges, density=None):
    """Anti-aliased binning: spread each sample over [rhot-half, rhot+half].

    Shared by the box LOS path (:func:`los_shell_fraction`) and the real-cell
    detector path (:func:`los_file_detector_rate`). A sample at ``rhot_c``
    with half-range ``half`` distributes its ``weight`` across the rhot
    interval and adds to each bin in proportion to its overlap with that
    interval. Total weight is conserved (zero-half samples drop entirely
    into their containing bin). Vectorised with ``np.add.at``.

    Parameters
    ----------
    density : 1-D array (len = n_bins) or None
        Optional per-bin density for a **Jacobian-weighted** distribution.
        With ``density=None`` (default) the sample is spread *uniformly in
        rhot* (overlap-weighted) -- the original behaviour, used by the box
        path. When ``density`` is given the within-sample weight goes as
        ``overlap * density_bin`` (re-normalised per sample so the total is
        unchanged). Passing ``density = DVOL`` makes the split follow the
        flux-shell volume, which is the physically correct distribution of a
        Cartesian cell's volume across shells near the axis (where
        ``dV/drhot -> 0`` so a cell holds more volume in its outer-rhot
        part). This removes the near-axis over-fill/starve artifact in
        ``f_det`` without the geometric path's "pin enclosed shells to 1"
        trick (which is inapplicable because ``f_det`` is not unity there).

    Returns
    -------
    binned weight array of shape ``(len(edges) - 1,)``.
    """
    rhot_c = np.asarray(rhot_c, dtype=float)
    half = np.asarray(half, dtype=float)
    weight = np.asarray(weight, dtype=float)
    n_bins = len(edges) - 1
    out = np.zeros(n_bins)
    if rhot_c.size == 0:
        return out

    # Clip to [0, 1] so weight isn't spilled outside the physical range.
    rhot_lo = np.clip(rhot_c - half, 0.0, 1.0)
    rhot_hi = np.clip(rhot_c + half, 0.0, 1.0)
    span = rhot_hi - rhot_lo

    i_lo = np.clip(np.searchsorted(edges, rhot_lo, side='right') - 1,
                   0, n_bins - 1)
    i_hi = np.clip(np.searchsorted(edges, rhot_hi, side='right') - 1,
                   0, n_bins - 1)
    n_per_cell = i_hi - i_lo + 1

    total = int(n_per_cell.sum())
    flat_cell = np.repeat(np.arange(rhot_c.size), n_per_cell)
    starts = np.concatenate(([0], np.cumsum(n_per_cell[:-1])))
    within_offset = np.arange(total) - np.repeat(starts, n_per_cell)
    flat_bin = np.repeat(i_lo, n_per_cell) + within_offset

    cell_lo_b = rhot_lo[flat_cell]
    cell_hi_b = rhot_hi[flat_cell]
    cell_span_b = span[flat_cell]
    cell_w_b = weight[flat_cell]
    bin_lo = edges[flat_bin]
    bin_hi = edges[flat_bin + 1]
    ovl = np.maximum(
        0.0, np.minimum(cell_hi_b, bin_hi) - np.maximum(cell_lo_b, bin_lo))

    # Per-(cell, bin) contribution before normalisation.
    if density is None:
        contrib = ovl
    else:
        contrib = ovl * np.asarray(density, dtype=float)[flat_bin]

    # Per-cell normalisation so each sample's total weight is preserved. For
    # density=None this norm equals the (clipped) span, recovering the original
    # ovl/span formula exactly.
    norm = np.zeros(rhot_c.size)
    np.add.at(norm, flat_cell, contrib)
    norm_b = norm[flat_cell]

    zero_span = cell_span_b == 0
    good = (~zero_span) & (norm_b > 0.0)
    # Default: zero-span samples deposit their full weight in their one bin.
    w = np.where(zero_span, cell_w_b, 0.0)
    w = np.where(good, cell_w_b * contrib / np.where(norm_b > 0.0, norm_b, 1.0), w)
    # Degenerate span>0 but zero density over the touched bins: fall back to
    # the uniform (overlap/span) split so weight is never lost.
    bad = (~zero_span) & (norm_b <= 0.0)
    if np.any(bad):
        w = np.where(
            bad,
            cell_w_b * ovl / np.where(cell_span_b > 0.0, cell_span_b, 1.0),
            w,
        )
    np.add.at(out, flat_bin, w)
    return out


# ----------------------------------------------------------------------
# LOS shell fractions: box chord (R, Z) grid OR real LOS cell cloud
# ----------------------------------------------------------------------

def los_shell_fraction(rhot_los, inside, R, Z, x_rhot, dvol_x_m3,
                       subgrid=True, rmag=None):
    """LOS-weighted fraction f(rhot) of each TRANSP flux shell.

    For every (R, Z) cell of the LOS grid we accumulate the toroidal
    volume element ``dV = 2*pi*R*dR*dZ`` (zero outside the LCFS) into
    bins centred on the TRANSP rhot grid, then

        f(rhot_i) = LOS_vol_bin(i) / DVOL_TRANSP_bin(i)

    so ``f`` lies in [0, 1] and ``sum(THNTX * DVOL * f) == Rate_tor`` by
    construction (modulo discretisation).

    Binning modes
    -------------
    * ``subgrid=False`` -- point binning: each cell's volume goes to a
      single bin chosen by its centre rhot. Suffers from aliasing when
      the rhot change across one (R, Z) cell exceeds one TRANSP bin
      width (typical for the Z direction at modest nZ).
    * ``subgrid=True`` (default) -- anti-aliased: each cell spans a rhot
      range estimated from ``|grad rhot| * cell_size`` (L1 corner half-
      range). Its volume is distributed uniformly across that range and
      added to every TRANSP bin proportional to overlap, vectorised via
      ``np.add.at``. Eliminates the period-N oscillation in f(rhot) and
      correctly accounts for the mismatch between a rectangular (R, Z)
      grid and curved flux surfaces.

    Parameters
    ----------
    rhot_los  : 2D array, rhot on the LOS (R, Z) grid
    inside    : 2D bool array, inside-LCFS mask (same shape)
    R, Z      : 1D LOS grid axes [m]
    x_rhot    : 1D TRANSP rhot grid (sorted ascending)
    dvol_x_m3 : 1D TRANSP DVOL converted to m^3 (same shape as x_rhot)
    subgrid   : bool, see above
    rmag      : float or None. Magnetic-axis R [m]. When given and inside the
                LOS R-band, every flux shell whose full poloidal extent stays
                within [R[0], R[-1]] is swept in its entirety by the chord, so
                ``f == 1`` exactly there; those inner bins are set to 1
                analytically (see "Enclosed-shell correction" below).

    Enclosed-shell correction
    -------------------------
    Near the axis the flux-shell volume Jacobian dV/drhot -> 0, while a single
    Cartesian (R, Z) cell spans many rhot bins. Both binning modes then
    mis-share that volume among the innermost bins (over-filling the centre,
    starving the next ring), producing spurious f < 1 dips even though those
    shells are entirely inside the R-band and therefore fully seen. We compute
    ``rhot_crit`` = the rhot of the innermost flux surface that reaches either
    R boundary (min over Z, inside the LCFS, of rhot at the Rmin/Rmax columns)
    and pin ``f = 1`` for every bin lying fully inside it.

    Returns
    -------
    f            : 1D weight, shape (len(x_rhot),)
    los_vol_bin  : 1D LOS toroidal volume per TRANSP shell [m^3]
    """
    R_arr = np.asarray(R, dtype=float)
    Z_arr = np.asarray(Z, dtype=float)
    dR = np.gradient(R_arr)
    dZ = np.gradient(Z_arr)
    Rg, _ = np.meshgrid(R_arr, Z_arr, indexing='ij')
    cell_vol = np.where(
        inside,
        2.0 * np.pi * Rg * dR[:, None] * dZ[None, :],
        0.0,
    )

    xs = np.asarray(x_rhot, dtype=float)
    edges = np.concatenate(([0.0], 0.5 * (xs[:-1] + xs[1:]), [1.0]))
    n_bins = len(edges) - 1

    if not subgrid:
        los_vol_bin, _ = np.histogram(
            rhot_los.ravel(), bins=edges, weights=cell_vol.ravel()
        )
    else:
        # rhot half-extent across the cell (L1 corner half-range)
        grad_R, grad_Z = np.gradient(rhot_los, R_arr, Z_arr)
        half = (0.5 * np.abs(grad_R) * dR[:, None]
                + 0.5 * np.abs(grad_Z) * dZ[None, :])
        keep = cell_vol.ravel() > 0.0
        los_vol_bin = _subgrid_bin(
            rhot_los.ravel()[keep], half.ravel()[keep],
            cell_vol.ravel()[keep], edges,
        )

    dvol = np.asarray(dvol_x_m3, dtype=float)
    with np.errstate(invalid='ignore', divide='ignore'):
        f = np.where(dvol > 0.0, los_vol_bin / dvol, 0.0)
    f = np.clip(f, 0.0, 1.0)

    # Enclosed-shell correction: pin f = 1 for flux shells fully inside the
    # R-band (see docstring). Only meaningful when the axis sits in the band;
    # otherwise the innermost shells lie outside the chord (f -> 0, not 1).
    if rmag is not None and R_arr[0] <= rmag <= R_arr[-1]:
        col_lo = rhot_los[0, :][inside[0, :]]    # Rmin column, inside LCFS
        col_hi = rhot_los[-1, :][inside[-1, :]]  # Rmax column, inside LCFS
        cols = np.concatenate([col_lo, col_hi])
        if cols.size:
            rhot_crit = float(cols.min())
            enclosed = edges[1:] <= rhot_crit    # bin upper edge inside band
            f = np.where(enclosed, 1.0, f)
    return f, los_vol_bin


def los_file_detector_rate(cells, eqs, xs_sorted, thntx_sorted,
                           thntx_to_si, dvol_m3, subgrid=True,
                           chord_kind="vertical", median=True):
    """Thermal-neutron quantities from a real LOS cell cloud.

    Each cell carries its own volume ``V`` and etendue weight ``C`` (with
    ``C = V * Omega/4pi``, ``Omega`` = detector solid angle from the cell).
    For the isotropic thermal source the *detector-reaching* rate is

        Rate_det   = sum_cells  eps(rhot_cell) * C_cell      [n/s]

    and the bare real-chord emission (no solid-angle weighting, for comparison
    with the idealised ``Rate_LOS``) is

        Rate_chord = sum_cells  eps(rhot_cell) * V_cell      [n/s]

    The detector response is also binned onto the TRANSP flux shells to give a
    LOS detector-coupling weight ``f_det(rhot) = C_bin / DVOL`` (dimensionless,
    ~Omega/4pi times the swept volume fraction) and the LOS-weighted
    emissivity ``THLOS_det = THNTX * f_det`` with the same closure as the
    geometric path: ``sum(THNTX * DVOL * f_det) == Rate_det``.

    Parameters
    ----------
    cells       : dict from :func:`read_los_file`
    eqs         : equilibrium object exposing ``rhot_at_points``
    xs_sorted   : TRANSP rhot grid (ascending)
    thntx_sorted: THNTX on ``xs_sorted`` (native CGS units)
    thntx_to_si : CGS->SI factor for THNTX
    dvol_m3     : TRANSP DVOL on ``xs_sorted`` in m^3
    subgrid     : anti-aliased binning of C, V across shells (default).
    median      : apply the **cosmetic** running-median de-speckle of ``f_det``
                  (default True): the 3-bin median outside ``rhot_crit`` for
                  vertical chords and the 5-bin median on the crossed shells
                  for horizontal chords. Set False to return the raw binned
                  ``f_det`` (useful to estimate the median's effect). Purely
                  cosmetic -- ``Rate_det``/``rho_50`` are summed from the cells
                  and are identical either way.
    chord_kind  : ``"vertical"`` (KM14) / ``"horizontal"`` (KM9) / ``"auto"``.
                  Controls the **enclosed-shell flattening + outside median3**
                  of ``f_det``: only meaningful when the chord pierces the
                  magnetic axis along Z (so inner flux surfaces are fully
                  enclosed within an R-band), i.e. a vertical pencil chord.
                  For a horizontal chord that wraps tangentially the
                  enclosed-shell construction has no analogue and is skipped
                  (``rhot_crit`` returned as ``None``, raw binned ``f_det``).
                  ``"auto"`` infers from the cell versor: vertical if
                  ``median(|w|) > 0.9``, else horizontal. Default is
                  ``"vertical"`` to preserve historical KM14 behaviour
                  byte-for-byte; KM9 should pass ``"horizontal"`` or
                  ``"auto"``. For non-vertical chords ``f_det`` is additionally
                  zeroed on flux shells *entirely* below the chord's closest
                  approach to the axis (``rhot_min``), which a grazing chord
                  never crosses (see the inner-shell block below); vertical
                  chords handle their core via ``rhot_crit`` flattening instead.

    Returns
    -------
    dict with keys
        rate_det, rate_chord, R_cells, Z_cells, C_cells, rhot_cells,
        inside_cells, eps_si, c_bin, v_bin, f_det, thkm14_det_si,
        rate_from_profile, cum_rhot, cum_frac, rho_med, rhot_crit, rhot_min,
        vol_tot, c_tot, n_inside.

    ``rhot_min`` is the chord's closest-approach rhot (set for non-vertical
    chords; ``None`` for vertical, which use ``rhot_crit``).

    NB the ``thkm14_det_si`` key name is preserved for backward
    compatibility; conceptually it is now ``THLOS_det`` for any diagnostic.
    """
    R = cells["R"]
    Z = cells["Z"]
    C = np.asarray(cells["C"], dtype=float)
    V = np.asarray(cells["V"], dtype=float)

    if chord_kind == "auto":
        w = np.asarray(cells.get("w", np.array([])), dtype=float)
        chord_kind = "vertical" if (w.size and np.median(np.abs(w)) > 0.9) \
            else "horizontal"

    rhot_c, inside_c, _ = eqs.rhot_at_points(R, Z)

    # Emissivity at each cell (n/m3/s), zeroed outside the LCFS. Reuse the same
    # padded-PCHIP interpolation used for the box grid (works on any shape).
    eps_si = thntx_on_grid(xs_sorted, thntx_sorted, rhot_c, inside_c) * thntx_to_si

    sig = eps_si * C      # per-cell detector-reaching rate
    rate_det = float(np.sum(sig))
    rate_chord = float(np.sum(eps_si * V))

    # Bin the etendue / volume onto the TRANSP flux shells (same edges as
    # los_shell_fraction: midpoints of xs, plus 0 and 1).
    edges = np.concatenate(([0.0], 0.5 * (xs_sorted[:-1] + xs_sorted[1:]), [1.0]))
    w_in = inside_c.astype(float)
    rhot_crit = None
    rhot_min = None
    if not subgrid:
        c_bin, _ = np.histogram(rhot_c, bins=edges, weights=C * w_in)
        v_bin, _ = np.histogram(rhot_c, bins=edges, weights=V * w_in)
    else:
        # Anti-aliased binning. The LOS cells are Cartesian while the flux
        # surfaces are ~cylindrical, so one cell spans several TRANSP rhot bins
        # (mostly in Z) -> point binning aliases f_det. Give each cell a rhot
        # half-span 0.5*(|d rhot/dR|*dR + |d rhot/dZ|*dZ) and spread its weight
        # over that range (same scheme as the box subgrid path). Local
        # |grad rhot| comes from a dense structured rhot field; per-cell sizes
        # are dZ = z-slice spacing and dR = sqrt(V/dZ) (a cell is ~one slice
        # thick, so V ~= dR*dx*dZ with dR ~= dx in the poloidal plane).
        from scipy.interpolate import RegularGridInterpolator
        nrg, nzg = 200, 600
        Rg1 = np.linspace(float(R.min()), float(R.max()), nrg)
        Zg1 = np.linspace(float(Z.min()), float(Z.max()), nzg)
        rhot_g, inside_g, _ = eqs.rhot_on_grid(Rg1, Zg1)   # (nrg, nzg), 'ij'
        gR, gZ = np.gradient(rhot_g, Rg1, Zg1)
        igR = RegularGridInterpolator((Rg1, Zg1), np.abs(gR),
                                      bounds_error=False, fill_value=0.0)
        igZ = RegularGridInterpolator((Rg1, Zg1), np.abs(gZ),
                                      bounds_error=False, fill_value=0.0)
        keep = inside_c & (C > 0.0)
        pts = np.column_stack([R[keep], Z[keep]])
        dZc = float(np.median(np.diff(np.unique(Z))))
        dRc = np.sqrt(np.maximum(V[keep], 0.0) / dZc)
        half = 0.5 * (igR(pts) * dRc + igZ(pts) * dZc)
        # DVOL-weighted (Jacobian) split: a Cartesian cell's volume is spread
        # across the shells it spans in proportion to each shell's volume, not
        # uniformly in rhot. This is the physically correct distribution and
        # removes the near-axis over-fill/starve dip in f_det. Total C / V are
        # conserved.
        c_bin = _subgrid_bin(rhot_c[keep], half, C[keep], edges, density=dvol_m3)
        v_bin = _subgrid_bin(rhot_c[keep], half, V[keep], edges, density=dvol_m3)

        # rhot_crit = innermost flux surface that reaches an R-boundary of the
        # cell footprint (same construction as the geometric enclosed-shell
        # correction). Shells inside it are fully swept by the chord. Only
        # meaningful for a vertical pencil chord (KM14); skipped for chords
        # that don't pierce the axis along Z (KM9).
        if chord_kind == "vertical" and Rg1[0] <= eqs.Rmag <= Rg1[-1]:
            col_lo = rhot_g[0, :][inside_g[0, :]]
            col_hi = rhot_g[-1, :][inside_g[-1, :]]
            cols = np.concatenate([col_lo, col_hi])
            if cols.size:
                rhot_crit = float(cols.min())

    # Un-crossed inner shells (non-vertical / grazing chords only). A vertical
    # pencil chord (KM14) geometrically *encloses* the innermost flux shells
    # even when no discrete cell centre lands there, so its core is handled by
    # the enclosed-shell flattening below (rhot_crit). A horizontal / grazing
    # chord (KM9) instead has a finite closest approach to the magnetic axis:
    # its minimum sampled rhot is rhot_min = min rhot over in-LCFS cells, and
    # flux shells *entirely* below rhot_min are never crossed -> no thermal
    # emission from them reaches the detector, so f_det there must be exactly
    # zero. The DVOL-weighted subgrid splat otherwise leaks a little C into
    # those bins through the symmetric cell tails (a spurious near-axis shelf
    # at ~half the plateau). Zero those bins and fold their leaked C into the
    # first crossed bin so Rate_det closure (total C) is preserved. The first
    # *partially* crossed bin (which straddles rhot_min) keeps its reduced
    # value, giving an honest roll-up from the closest-approach radius.
    if chord_kind != "vertical":
        keep_in = inside_c & (C > 0.0)
        if keep_in.any():
            rhot_min = float(rhot_c[keep_in].min())
            below = edges[1:] <= rhot_min      # bins entirely below the chord
            if below.any() and not below.all():
                first = int(np.argmax(~below))  # first bin not fully below
                c_bin = c_bin.copy()
                v_bin = v_bin.copy()
                c_bin[first] += c_bin[below].sum()
                v_bin[first] += v_bin[below].sum()
                c_bin[below] = 0.0
                v_bin[below] = 0.0

    with np.errstate(invalid="ignore", divide="ignore"):
        f_det = np.where(dvol_m3 > 0.0, c_bin / dvol_m3, 0.0)

    # Enclosed-shell flattening. For shells fully inside the chord footprint
    # (rhot <= rhot_crit) the wide R-band captures the whole poloidal cross
    # section, so the chord samples a constant toroidal+solid-angle fraction
    # and f_det is physically *flat*. The few innermost bins are nonetheless
    # under-sampled (the ~13 mm cells barely resolve the tiny axis tube, where
    # DVOL -> 0), so f_det dips/rises there. Pin the enclosed region to its
    # DVOL-weighted mean -- the detector analogue of the geometric path
    # pinning f = 1 there, but to the real plateau value since f_det != 1.
    # The volume weighting naturally down-weights the corrupt tiny-DVOL inner
    # bins, and Sum(DVOL*f_det) over the region is unchanged so Rate_det
    # closure is preserved. Skipped for non-vertical chords (rhot_crit None).
    if rhot_crit is not None:
        enclosed = edges[1:] <= rhot_crit
        dsum = float(dvol_m3[enclosed].sum())
        if enclosed.any() and dsum > 0.0:
            f_det = np.where(enclosed, float(c_bin[enclosed].sum()) / dsum, f_det)
        # Light 3-bin running median *outside* rhot_crit only, to suppress the
        # ~1-2% finite-cell sampling wiggle (the physical bump/roll-off
        # survives; the flat enclosed plateau is left untouched). Cosmetic --
        # Rate_det (= sum eps*C) is computed directly from the cells and is
        # unchanged.
        if subgrid and median:
            outside = edges[1:] > rhot_crit
            f_det = np.where(outside, _running_median3(f_det), f_det)
    elif subgrid and median and chord_kind != "vertical":
        # Horizontal / grazing chord: no enclosed plateau, but the regular cell
        # lattice (uniform z-slices x regular along-chord steps) beats against
        # the curved flux shells, leaving a quasi-periodic ~15% wiggle that the
        # subgrid splat alone cannot remove (widening the splat is non-monotonic
        # -- a classic aliasing beat, not under-smearing). Apply a 5-bin running
        # median as a cosmetic de-speckle on the *crossed* shells only (above
        # rhot_min), so the zeroed un-crossed core is untouched. Rate_det
        # (= sum eps*C) is summed from the cells directly and is unchanged.
        idx = np.where(edges[1:] > (rhot_min if rhot_min is not None else 0.0))[0]
        if idx.size >= 3:
            f_det = f_det.copy()
            f_det[idx] = _running_median(f_det[idx], 5)

    thntx_si = thntx_sorted * thntx_to_si
    thkm14_det_si = thntx_si * f_det
    # Closure check: profile-reconstructed rate should match Rate_det.
    rate_from_profile = float(np.sum(thntx_si * dvol_m3 * f_det))

    # Cumulative detector-signal fraction vs rhot -> signal-median radius.
    order = np.argsort(rhot_c)
    cum_rhot = rhot_c[order]
    cum_sig = np.cumsum(sig[order])
    tot = cum_sig[-1]
    cum_frac = cum_sig / tot if tot > 0 else np.zeros_like(cum_sig)
    rho_med = float(np.interp(0.5, cum_frac, cum_rhot)) if tot > 0 else float("nan")

    return {
        "rate_det": rate_det,
        "rate_chord": rate_chord,
        "R_cells": R,
        "Z_cells": Z,
        "C_cells": C,
        "rhot_cells": rhot_c,
        "inside_cells": inside_c,
        "eps_si": eps_si,
        "c_bin": c_bin,
        "v_bin": v_bin,
        "f_det": f_det,
        "thkm14_det_si": thkm14_det_si,
        "rate_from_profile": rate_from_profile,
        "cum_rhot": cum_rhot,
        "cum_frac": cum_frac,
        "rho_med": rho_med,
        "rhot_crit": rhot_crit,
        "rhot_min": rhot_min,
        "vol_tot": float(np.sum(V * w_in)),
        "c_tot": float(np.sum(C * w_in)),
        "n_inside": int(np.count_nonzero(inside_c)),
        "chord_kind": chord_kind,
    }
