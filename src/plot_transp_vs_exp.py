#!/usr/bin/env python3
"""
plot_transp_vs_exp.py -- TRANSP output profiles (NE, NI, TE, TI) with measured
HRTS / KK3 (and optionally core CX) data overlaid, styled like the paper figure:

    (b) n_i (black) and n_e (red)   [1e20 m^-3]
    (c) T_i = T_e (black line)      [keV]          vs  r/a = sqrt(Psi_tor,norm)

Standalone: no OMFIT, no Ufiles, no input fits. Only needs numpy, matplotlib,
netCDF4 and (for the measured data) the JET `ppf` module.

Equilibrium
-----------
Measurement positions (R, Z) are mapped onto rho_tor_norm (TRANSP X / XB) using
the flux-surface geometry stored in the TRANSP CDF itself -- the Fourier moments
RMCmm/RMSmm/YMCmm/YMSmm on the XB grid:

    R(theta) = RMC00 + sum_m [RMCmm cos(m theta) + RMSmm sin(m theta)]
    Z(theta) = YMC00 + sum_m [YMCmm cos(m theta) + YMSmm sin(m theta)]

For each point the distance from the magnetic axis (RAXIS, YAXIS) is compared
with the distance of every flux surface along the same ray; rho is interpolated
between surfaces (linear extrapolation beyond the LCFS, so SOL points get
rho > 1). This is the same equilibrium TRANSP used, so measured data and TRANSP
profiles are consistent by construction. If the CDF has no moments, the script
falls back to an outboard-midplane mapping with RMJMP (Z ignored, with warning).

Use --plot-geometry to check the mapping (flux surfaces + measurement positions).

Examples
--------
    python plot_transp_vs_exp.py 104990 M01 --time 50.0
    python plot_transp_vs_exp.py 104990 M01 --time 50.0 --twin 0.2 \\
        --hrts-te-range 0.65 1.0 --kk3-range 0.0 0.7 --kk3-seq 2 --save-data
    python plot_transp_vs_exp.py 104990 M01 --time 50.0 --cdf ./104990M01.CDF --q --save
    python plot_transp_vs_exp.py 104990 M01 --time 50.0 --nlim 0.05 0.55 --tlim 0 8 --xlim 0 1.05
    python plot_transp_vs_exp.py 104990 M01 --time 50.0 --pub --save        # legend: HRTS, ECE
    python plot_transp_vs_exp.py 104990 M01 --time 50.0 --kk3-label 'ECE (KK3)'

The figure is shown on screen; it is written to file only with --save
(--save-data writes the plotted data as ASCII).

Time is JET time [s]; TRANSP time = JET time - 40 s.
"""

import argparse
import os
import re
import sys

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator, FormatStrFormatter, MaxNLocator

try:
    import ppf
except ImportError:  # only available on the JET analysis cluster
    ppf = None

JET_T_OFFSET = 40.0


# ----------------------------------------------------------------------------
# Locating the run (same search order as profiles.find_run_files)
# ----------------------------------------------------------------------------
def find_run_files(pulse, runid, verbose=True):
    """Return (cdf_path, namelist_path) for <pulse><runid>.

    Search order: $JET_DATA_DIR (default ~/jet/data)/<pulse>/<runid>/,
    Heimdall results tree (as-is, then lowercase), Heimdall warehouse.
    """
    tr_seq = runid.replace(str(pulse), '')
    local_base = os.environ.get('JET_DATA_DIR', os.path.join(os.path.expanduser('~'), 'jet', 'data'))
    local_dir = os.path.join(local_base, str(pulse), runid)
    result_dir = f'/common/transp_shared/Data/result/JET/{pulse}/{runid}'
    result_dir2 = f'/common/transp_shared/Data/result/JET/{pulse:d}/{tr_seq.upper()}'
    whshot = f'{pulse:07d}'
    whousedir = f'/common/transp_shared/Data/whouse/JET/{whshot[0:4]}___/{whshot}/{tr_seq.lower()}'
    candidates = [
        (os.path.join(local_dir, f'{pulse}{runid}.CDF'), os.path.join(local_dir, f'{pulse}{runid}TR.DAT'),
         f'local data tree ({local_base})'),
        (f'{result_dir}/{pulse}{runid}.CDF', f'{result_dir}/{pulse}{runid}TR.DAT', 'results directory'),
        (f'{result_dir2}/{pulse:d}{tr_seq.lower()}.cdf', f'{result_dir2}/{pulse:d}{tr_seq.lower()}TR.DAT',
         'results directory (lowercase)'),
        (f'{whousedir}/{whshot}{tr_seq.lower()}.cdf', f'{whousedir}/{whshot}{tr_seq.lower()}TR.DAT',
         'warehouse directory'),
    ]
    for cdf_c, nml_c, where in candidates:
        if os.path.isfile(cdf_c):
            if verbose:
                print(f'Using netCDF file found in {where}: {cdf_c}')
            return cdf_c, nml_c
    raise FileNotFoundError(f'CDF for {pulse}{runid} not found. Tried: ' + ', '.join(c[0] for c in candidates))


def warn(msg):
    print(f'WARNING: {msg}', file=sys.stderr)


def time_indices(t, t0, window):
    """Indices of t within [t0 - window/2, t0 + window/2]; nearest one if window<=0 or empty."""
    t = np.asarray(t, dtype=float)
    if window and window > 0:
        idx = np.where(np.abs(t - t0) <= window / 2.0)[0]
        if idx.size:
            return idx
    return np.array([int(np.argmin(np.abs(t - t0)))])


# ----------------------------------------------------------------------------
# TRANSP output
# ----------------------------------------------------------------------------
class TranspRun:
    def __init__(self, cdf_path):
        self.path = cdf_path
        self.ds = Dataset(cdf_path)
        self.keys = set(self.ds.variables.keys())
        self.t = self.get('TIME3')  # TRANSP time [s]

    def has(self, key):
        return key in self.keys

    def get(self, key):
        return np.asarray(np.ma.filled(self.ds.variables[key][:], np.nan), dtype=float)

    def units(self, key):
        return getattr(self.ds.variables[key], 'units', '').strip()

    def radial_grid_name(self, key):
        dims = self.ds.variables[key].dimensions
        return dims[-1] if dims[-1] in self.keys else 'X'

    def _weights(self, idx):
        if idx.size == 1:
            return np.ones(1)
        dt = np.gradient(self.t)
        return dt[idx]

    def profile(self, key, t_tr, window=0.0):
        """Profile averaged (dt-weighted) over the window. Returns (rho, y, t_used)."""
        data = self.get(key)
        rho = self.get(self.radial_grid_name(key))
        idx = time_indices(self.t, t_tr, window)
        w = self._weights(idx)
        y = np.average(data[idx], axis=0, weights=w)
        x = np.average(rho[idx], axis=0, weights=w)
        return x, y, self.t[idx]


class TranspEquilibrium:
    """(R, Z) [m] -> rho_tor_norm using the flux-surface moments in the TRANSP CDF."""

    MOM_RE = re.compile(r'^([RY])M([CS])(\d{2})$')

    def __init__(self, run, ntheta=720):
        self.run = run
        self.ntheta = ntheta
        self.t = run.t
        mom = {}
        for k in run.keys:
            m = self.MOM_RE.match(k)
            if m:
                mom[(m.group(1), m.group(2), int(m.group(3)))] = k
        self.has_moments = ('R', 'C', 0) in mom and ('Y', 'C', 0) in mom and \
            any(key[2] >= 1 for key in mom)

        if self.has_moments:
            ref = mom[('R', 'C', 0)]
            u = run.units(ref).lower()
            self.scale = 1.0 if u in ('m', 'meters', 'metres') else 0.01  # TRANSP default: cm
            self.mmax = max(key[2] for key in mom)
            self.mom = {key: run.get(name) * self.scale for key, name in mom.items()}
            self.rho_b = run.get(run.radial_grid_name(ref))  # XB
            print(f'Equilibrium: TRANSP flux-surface moments up to m={self.mmax} '
                  f'on {run.radial_grid_name(ref)} grid')
        else:
            warn('No RMCmm/YMSmm flux-surface moments in the CDF: falling back to outboard-midplane '
                 'mapping with RMJMP (Z of the measurement is ignored!).')
            if not run.has('RMJMP'):
                raise KeyError('Neither flux-surface moments nor RMJMP found in CDF: cannot map (R,Z).')
            self.scale = 0.01 if run.units('RMJMP').lower() != 'm' else 1.0
            self.rmjmp = run.get('RMJMP') * self.scale
            self.rho_b = run.get(run.radial_grid_name('RMJMP'))

        # magnetic axis
        if run.has('RAXIS') and run.has('YAXIS'):
            s = 0.01 if run.units('RAXIS').lower() != 'm' else 1.0
            self.Rax = run.get('RAXIS') * s
            self.Zax = run.get('YAXIS') * s
        elif self.has_moments:
            self.Rax = self.mom[('R', 'C', 0)][:, 0]
            self.Zax = self.mom[('Y', 'C', 0)][:, 0]
        else:
            self.Rax = self.rmjmp[:, 0]
            self.Zax = np.zeros_like(self.Rax)

    def it(self, t_tr):
        return int(np.argmin(np.abs(self.t - t_tr)))

    def axis(self, t_tr):
        i = self.it(t_tr)
        return self.Rax[i], self.Zax[i]

    def surfaces(self, t_tr, ntheta=None):
        """R, Z [m] of the XB flux surfaces, shape (nxb, ntheta)."""
        if not self.has_moments:
            raise RuntimeError('No flux-surface moments available.')
        i = self.it(t_tr)
        th = np.linspace(0.0, 2.0 * np.pi, ntheta or self.ntheta, endpoint=False)
        R = np.repeat(self.mom[('R', 'C', 0)][i][:, None], th.size, axis=1)
        Z = np.repeat(self.mom[('Y', 'C', 0)][i][:, None], th.size, axis=1)
        for m in range(1, self.mmax + 1):
            c, s = np.cos(m * th)[None, :], np.sin(m * th)[None, :]
            for key, trig, arr in ((('R', 'C', m), c, R), (('R', 'S', m), s, R),
                                   (('Y', 'C', m), c, Z), (('Y', 'S', m), s, Z)):
                if key in self.mom:
                    arr += self.mom[key][i][:, None] * trig
        return R, Z

    def rz_to_rhot(self, R, Z, t_tr):
        R = np.atleast_1d(np.asarray(R, dtype=float))
        Z = np.broadcast_to(np.atleast_1d(np.asarray(Z, dtype=float)), R.shape).copy()
        i = self.it(t_tr)
        rho_b = self.rho_b[i]
        out = np.full(R.shape, np.nan)
        good = np.isfinite(R) & np.isfinite(Z)

        if not self.has_moments:
            # outboard midplane only
            rm = self.rmjmp[i]
            out[good] = np.interp(R[good], rm, rho_b, left=np.nan, right=np.nan)
            return out

        Rax, Zax = self.Rax[i], self.Zax[i]
        Rs, Zs = self.surfaces(t_tr)
        phis = np.arctan2(Zs - Zax, Rs - Rax)
        dsurf = np.hypot(Rs - Rax, Zs - Zax)
        php = np.arctan2(Z[good] - Zax, R[good] - Rax)
        dp = np.hypot(R[good] - Rax, Z[good] - Zax)

        # distance of every surface along each point's ray from the axis
        D = np.empty((rho_b.size, php.size))
        for j in range(rho_b.size):
            o = np.argsort(phis[j])
            D[j] = np.interp(php, phis[j, o], dsurf[j, o], period=2.0 * np.pi)

        rr = np.concatenate(([0.0], rho_b))
        res = np.empty(php.size)
        for p in range(php.size):
            dd = np.concatenate(([0.0], D[:, p]))
            if dp[p] <= dd[-1]:
                res[p] = np.interp(dp[p], dd, rr)
            else:  # outside LCFS: linear extrapolation from the last two surfaces
                slope = (rr[-1] - rr[-2]) / (dd[-1] - dd[-2])
                res[p] = rr[-1] + slope * (dp[p] - dd[-1])
        out[good] = res
        return out


# ----------------------------------------------------------------------------
# PPF readers
# ----------------------------------------------------------------------------
def _require_ppf():
    if ppf is None:
        raise ImportError('`ppf` module not available: measured data cannot be read '
                          '(run on the JET cluster, or use --no-hrts --no-kk3).')


def ppfread(pulse, dda, dty, uid, seq):
    """ppf.ppfdata wrapper -> dict(data, x, t, seq, units) or None. 2D data reshaped to (t, x)."""
    _require_ppf()
    out = ppf.ppfdata(pulse, dda, dty, seq=seq, uid=uid, device='JET')
    data, x, t = (np.atleast_1d(np.asarray(a, dtype=float)) for a in out[:3])
    ier = out[-1]
    if ier != 0 or data.size == 0:
        return None
    if x.size > 1 and t.size > 1 and data.size == x.size * t.size:
        data = data.reshape(t.size, x.size)
    return {'data': data, 'x': x, 't': t, 'seq': out[11], 'units': str(out[6]).strip()}


def _window_average(values, errors):
    """Mean over the first axis; error of the mean propagated from the per-slice errors."""
    n = values.shape[0]
    if n == 1:
        return values[0], errors[0]
    with np.errstate(invalid='ignore'):
        y = np.nanmean(values, axis=0)
        e = np.sqrt(np.nansum(errors ** 2, axis=0)) / np.sum(np.isfinite(values), axis=0)
    return y, e


def read_hrts(pulse, dty, uid, seq, t_jet, window, eq, dda='HRTS'):
    """HRTS NE [m^-3] or TE [eV] mapped onto rho_tor_norm."""
    d = ppfread(pulse, dda, dty, uid, seq)
    if d is None:
        warn(f'No {dda}/{dty}/{uid}/{seq} for pulse {pulse}. Skipping.')
        return None
    de = ppfread(pulse, dda, f'D{dty}', uid, seq)
    z = ppfread(pulse, dda, 'Z', uid, seq)
    if z is None:
        warn(f'No {dda}/Z for pulse {pulse}: cannot map HRTS positions. Skipping.')
        return None

    idx = time_indices(d['t'], t_jet, window)
    vals = np.atleast_2d(d['data'])[idx]
    errs = np.atleast_2d(de['data'])[idx] if de is not None else 0.1 * vals
    y, e = _window_average(vals, errs)

    # positions: x vector of the Z signal is R, its data is Z (as in plot_validation)
    R = z['x']
    zdat = z['data']
    Zrow = zdat[time_indices(z['t'], t_jet, 0)[0]] if zdat.ndim == 2 else zdat
    if R.size != y.size:
        R = d['x']
        Zrow = np.interp(R, z['x'], Zrow)
    rho = eq.rz_to_rhot(R, Zrow, t_jet - JET_T_OFFSET)
    return {'rho': rho, 'y': y, 'err': e, 'R': R, 'Z': Zrow, 't': d['t'][idx],
            'label': f'{dda}/{uid}/{d["seq"]}'}


def read_kk3(pulse, uid, seq, t_jet, window, eq, side='lfs', dda='KK3'):
    """KK3 Te [eV] from TEnn/RCnn channels, mapped onto rho_tor_norm at Z = Z_axis."""
    gen = ppfread(pulse, dda, 'GEN', uid, seq)
    if gen is None:
        warn(f'No {dda}/GEN/{uid}/{seq} for pulse {pulse}. Skipping KK3.')
        return None
    nch = gen['x'].size
    TE, RC, tused = [], [], None
    for ich in range(1, nch + 1):
        te = ppfread(pulse, dda, f'TE{ich:02d}', uid, seq)
        rc = ppfread(pulse, dda, f'RC{ich:02d}', uid, seq)
        if te is None or rc is None:
            continue
        it_te = time_indices(te['t'], t_jet, window)
        it_rc = time_indices(rc['t'], t_jet, window)
        TE.append(np.nanmean(te['data'][it_te]))
        RC.append(np.nanmean(rc['data'][it_rc]))
        tused = te['t'][it_te]
    TE, RC = np.array(TE), np.array(RC)
    ok = np.isfinite(TE) & np.isfinite(RC) & (RC > 0) & (TE > 0)
    TE, RC = TE[ok], RC[ok]
    if TE.size == 0:
        warn(f'No valid {dda} channels for pulse {pulse}. Skipping KK3.')
        return None

    Rax, Zax = eq.axis(t_jet - JET_T_OFFSET)
    sel = {'lfs': RC > Rax, 'hfs': RC < Rax, 'both': np.ones(RC.size, bool)}[side]
    TE, RC = TE[sel], RC[sel]
    o = np.argsort(RC)
    TE, RC = TE[o], RC[o]
    Z = np.full(RC.size, Zax)
    rho = eq.rz_to_rhot(RC, Z, t_jet - JET_T_OFFSET)
    return {'rho': rho, 'y': TE, 'err': 0.1 * TE, 'R': RC, 'Z': Z, 't': tused,
            'label': f'{dda}/{uid}/{seq}'}


def read_cx_core(pulse, dda, uid, seq, t_jet, window, eq, dty='TI'):
    """Core CX Ti [eV]; position convention as in plot_validation.py (POS: x=R, data=Z)."""
    d = ppfread(pulse, dda, dty, uid, seq)
    if d is None and dty != 'TI':
        d = ppfread(pulse, dda, 'TI', uid, seq)
    if d is None:
        warn(f'No {dda}/{dty}/{uid}/{seq} for pulse {pulse}. Skipping.')
        return None
    hi = ppfread(pulse, dda, 'TIHI', uid, seq)
    pos = ppfread(pulse, dda, 'POS', uid, seq)
    if pos is None:
        warn(f'No {dda}/POS: cannot map CX positions. Skipping.')
        return None
    idx = time_indices(d['t'], t_jet, window)
    vals = np.atleast_2d(d['data'])[idx]
    errs = (np.atleast_2d(hi['data'])[idx] - vals) if hi is not None else 0.1 * vals
    y, e = _window_average(vals, errs)
    pdat = np.atleast_2d(pos['data'])
    Z = np.nanmean(pdat[idx], axis=0) if pdat.shape[0] == d['t'].size else pdat[0]
    R = pos['x']
    rho = eq.rz_to_rhot(R, Z, t_jet - JET_T_OFFSET)
    return {'rho': rho, 'y': y, 'err': e, 'R': R, 'Z': Z, 't': d['t'][idx],
            'label': f'{dda}/{uid}/{d["seq"]}'}


def read_efit_q(pulse, t_jet, window, run, dda='EFTP', uid='jetppf', seq=0):
    """EFIT q profile mapped onto rho_tor_norm.

    EFIT Q is given on normalised poloidal flux psi_N. It is mapped with EFIT's own
    toroidal flux FTOR (same equilibrium): rho_tor = sqrt(FTOR / FTOR(LCFS)).
    If FTOR is unavailable, TRANSP PLFLX/TRFLX at the nearest time is used instead.
    Returns (rho_tor, q, label) or None.
    """
    q = ppfread(pulse, dda, 'Q', uid, seq)
    if q is None:
        warn(f'No {dda}/Q/{uid}/{seq} for pulse {pulse}.')
        return None
    idx = time_indices(q['t'], t_jet, window)
    qv = np.nanmean(np.atleast_2d(q['data'])[idx], axis=0)
    psin = q['x']

    ftor = ppfread(pulse, dda, 'FTOR', uid, seq)
    if ftor is not None:
        it_f = time_indices(ftor['t'], t_jet, window)
        f = np.nanmean(np.atleast_2d(ftor['data'])[it_f], axis=0)
        f = np.interp(psin, ftor['x'], f) if ftor['x'].size == f.size else f
        rho = np.sqrt(np.clip(f / f[-1], 0.0, None))
        how = f'{dda}/FTOR'
    elif run.has('PLFLX') and run.has('TRFLX'):
        i = int(np.argmin(np.abs(run.t - (t_jet - JET_T_OFFSET))))
        pl, tr = run.get('PLFLX')[i], run.get('TRFLX')[i]
        rho = np.sqrt(np.interp(psin, (pl - pl[0]) / (pl[-1] - pl[0]), tr / tr[-1]))
        how = 'TRANSP PLFLX/TRFLX'
    else:
        warn(f'No {dda}/FTOR and no TRANSP PLFLX/TRFLX: cannot map EFIT q onto rho_tor.')
        return None
    # EFIT grid is coarse (~33 psi_N points) and rho_tor ~ sqrt(psi_N) stretches the first
    # interval over r/a ~ 0-0.25: resample on a fine psi_N grid (linear in psi_N) before mapping.
    o = np.argsort(psin)
    psin_f = np.linspace(psin[o][0], psin[o][-1], 1001)
    rho = np.interp(psin_f, psin[o], rho[o] ** 2) ** 0.5  # rho^2 ~ normalised toroidal flux, linear in psi_N near axis
    qv = np.interp(psin_f, psin[o], qv[o])
    o = np.argsort(rho)
    print(f'q profile: {dda}/Q/{uid}/{q["seq"]} at {q["t"][idx][0]:.4f} s'
          + (f' (avg of {idx.size})' if idx.size > 1 else '') + f', psi_N -> rho_tor via {how}')
    return rho[o], qv[o], f'{dda}/{uid}/{q["seq"]}'


def select_range(d, rng):
    """Keep finite points with rng[0] <= rho <= rng[1], sorted by rho."""
    if d is None:
        return None
    m = np.isfinite(d['rho']) & np.isfinite(d['y']) & (d['rho'] >= rng[0]) & (d['rho'] <= rng[1])
    o = np.argsort(d['rho'][m])
    out = dict(d)
    for k in ('rho', 'y', 'err', 'R', 'Z'):
        out[k] = np.asarray(d[k])[m][o]
    if out['rho'].size == 0:
        warn(f'{d["label"]}: no points in rho range {rng}.')
    return out


# ----------------------------------------------------------------------------
# Plotting
# ----------------------------------------------------------------------------
STYLE = {
    'HRTS': dict(marker='o', ms=5, mfc='none', mew=1.0, ls='none', capsize=2, elinewidth=0.8),
    'KK3': dict(marker='s', ms=5, mfc='none', mew=1.0, ls='none', capsize=2, elinewidth=0.8),
    'CX': dict(marker='^', ms=5, mfc='none', mew=1.0, ls='none', capsize=2, elinewidth=0.8),
}
C_E, C_I = 'tab:red', 'black'
C_DATA = {'HRTS': '#d62728', 'KK3': '#ff7f0e', 'CX': '#555555'}
FS_LABEL, FS_TICK, FS_PANEL, FS_CURVE = 19, 11, 17, 17

# Panel geometry measured from the reference figure (paper Fig. 1):
# each panel is 253 x 187 px -> height/width = 0.74, gap between panels = 0.20 x panel height.
PANEL_W = 3.0                   # inches
PANEL_H = 0.74 * PANEL_W        # 2.22 in
PANEL_GAP = 0.20 * PANEL_H      # 0.44 in
MARGIN = dict(left=0.85, right=0.20, bottom=0.62, top=0.38)  # inches


def _overlay(ax, d, kind, scale):
    if d is None or d['rho'].size == 0:
        return
    ax.errorbar(d['rho'], d['y'] * scale, yerr=d['err'] * scale, color=C_DATA[kind],
                label=d['label'], zorder=1, **STYLE[kind])


def _panel_label(ax, txt):
    ax.text(-0.02, 1.01, txt, transform=ax.transAxes, ha='right', va='bottom', fontsize=FS_PANEL)


def make_figure(prof, exp, title=None, with_q=False, limits=None, box=False):
    """limits: dict with optional 'x', 'q', 'n', 'T' -> (min, max); None entries keep autoscale.
    box=False reproduces the reference style (left/bottom axes, outward ticks);
    box=True draws a closed frame with inward ticks on all sides."""
    limits = limits or {}
    n = 3 if with_q else 2
    fw = MARGIN['left'] + PANEL_W + MARGIN['right']
    fh = MARGIN['bottom'] + n * PANEL_H + (n - 1) * PANEL_GAP + MARGIN['top']
    fig = plt.figure(figsize=(fw, fh))
    axes = []
    for i in range(n):  # top panel first, fixed size in inches
        y0 = MARGIN['bottom'] + (n - 1 - i) * (PANEL_H + PANEL_GAP)
        rect = [MARGIN['left'] / fw, y0 / fh, PANEL_W / fw, PANEL_H / fh]
        axes.append(fig.add_axes(rect, sharex=axes[0] if axes else None))
    axes = np.array(axes)
    for a in axes[:-1]:
        a.tick_params(labelbottom=False)
    k = 0
    if with_q:
        ax = axes[0]
        ax.plot(prof['q'][0], prof['q'][1], color='k', lw=1.5)
        ax.set_ylabel('q profile', fontsize=FS_LABEL)
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))  # 1, 2, 3, ... as in the reference
        _panel_label(ax, '(a)')
        k = 1

    # ---------------- density
    ax = axes[k]
    sc = 1e-20
    xe, ne = prof['NE']
    xi, ni = prof['NI']
    _overlay(ax, exp.get('HRTS_NE'), 'HRTS', sc)
    ax.plot(xe, ne * sc, color=C_E, lw=1.5, zorder=3)
    ax.plot(xi, ni * sc, color=C_I, lw=1.5, zorder=3)
    ymax = np.nanmax(ne[xe <= 1.0] * sc)
    ax.text(0.62, np.interp(0.62, xe, ne * sc) + 0.15 * ymax, r'n$_\mathrm{e}$', color=C_E,
            fontsize=FS_CURVE, ha='center', va='bottom')
    ax.text(0.33, np.interp(0.33, xi, ni * sc) - 0.06 * ymax, r'n$_\mathrm{i}$', color=C_I,
            fontsize=FS_CURVE, ha='center', va='top')
    ax.set_ylabel(r'n (10$^{20}$ m$^{-3}$)', fontsize=FS_LABEL)
    _panel_label(ax, '(b)')
    if exp.get('HRTS_NE') is not None and exp['HRTS_NE']['rho'].size:
        ax.legend(loc='upper right', fontsize=9, frameon=False, handletextpad=0.3)

    # ---------------- temperature
    ax = axes[k + 1]
    sc = 1e-3
    xte, te = prof['TE']
    xti, ti = prof['TI']
    for key, kind in (('HRTS_TE', 'HRTS'), ('KK3_TE', 'KK3'), ('CX_TI', 'CX')):
        for d in exp.get(key, []) if key == 'CX_TI' else [exp.get(key)]:
            _overlay(ax, d, kind, sc)
    inside = xte <= 1.0
    same = np.nanmax(np.abs(np.interp(xte, xti, ti) - te)[inside]) < 0.01 * np.nanmax(te[inside])
    tmax = np.nanmax(np.r_[te[inside], ti[xti <= 1.0]] * sc)
    if same:
        ax.plot(xti, ti * sc, color=C_I, lw=1.5, zorder=3)
        t1 = ax.text(0.45, 0.55 * tmax, r'T$_\mathrm{i}$ = ', color=C_I, fontsize=FS_CURVE,
                     ha='left', va='bottom')
        ax.annotate(r'T$_\mathrm{e}$', xy=(1, 0), xycoords=t1, color=C_E, fontsize=FS_CURVE,
                    ha='left', va='bottom')
    else:
        ax.plot(xte, te * sc, color=C_E, lw=1.5, zorder=3)
        ax.plot(xti, ti * sc, color=C_I, lw=1.5, zorder=3)
        ax.text(0.55, np.interp(0.55, xte, te * sc) + 0.06 * tmax, r'T$_\mathrm{e}$', color=C_E,
                fontsize=FS_CURVE, ha='left', va='bottom')
        ax.text(0.30, np.interp(0.30, xti, ti * sc) - 0.06 * tmax, r'T$_\mathrm{i}$', color=C_I,
                fontsize=FS_CURVE, ha='right', va='top')
    ax.set_ylabel('T (keV)', fontsize=FS_LABEL)
    ax.set_ylim(bottom=0.0)
    _panel_label(ax, '(c)')
    if ax.get_legend_handles_labels()[0]:
        ax.legend(loc='upper right', fontsize=9, frameon=False, handletextpad=0.3)

    # axis limits (None in a pair keeps that end automatic)
    panel_ax = {'n': axes[k], 'T': axes[k + 1]}
    if with_q:
        panel_ax['q'] = axes[0]
    for key, lim in limits.items():
        if lim is not None and key in panel_ax:
            panel_ax[key].set_ylim(*lim)
    axes[-1].set_xlim(*(limits.get('x') or (0.0, 1.0)))
    axes[-1].xaxis.set_major_locator(MultipleLocator(0.25))
    axes[-1].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    for a in axes:
        a.xaxis.set_minor_locator(AutoMinorLocator(2))   # one minor tick between majors, as in the reference
        a.yaxis.set_minor_locator(AutoMinorLocator(2))
        if box:
            a.tick_params(which='both', direction='in', top=True, right=True)
        else:
            a.spines['top'].set_visible(False)
            a.spines['right'].set_visible(False)
            a.tick_params(which='both', direction='out', top=False, right=False)
        a.tick_params(which='both', labelsize=FS_TICK)
        a.tick_params(which='major', length=5)
        a.tick_params(which='minor', length=2.5)
    axes[-1].set_xlabel('r/a', fontsize=FS_LABEL)
    if title:
        fig.text(0.5, 1 - 0.08 / fh, title, fontsize=9, ha='center', va='top')
    return fig


def plot_geometry(eq, t_tr, exp_all, title=''):
    fig, ax = plt.subplots(figsize=(5, 7))
    if eq.has_moments:
        R, Z = eq.surfaces(t_tr)
        nxb = R.shape[0]
        for j in list(range(nxb // 10 - 1, nxb, max(1, nxb // 10))) + [nxb - 1]:
            ax.plot(np.r_[R[j], R[j, 0]], np.r_[Z[j], Z[j, 0]], color='0.6', lw=0.7)
        ax.plot(np.r_[R[-1], R[-1, 0]], np.r_[Z[-1], Z[-1, 0]], color='k', lw=1.2, label='LCFS (TRANSP)')
    Rax, Zax = eq.axis(t_tr)
    ax.plot(Rax, Zax, 'k+', ms=10)
    for name, d in exp_all.items():
        if d is None:
            continue
        sc = ax.scatter(d['R'], d['Z'], c=d['rho'], s=14, cmap='viridis', vmin=0, vmax=1.1, label=name)
    if exp_all:
        fig.colorbar(sc, ax=ax, label=r'mapped $\rho_{tor}$')
    ax.set_aspect('equal')
    ax.set_xlabel('R (m)')
    ax.set_ylabel('Z (m)')
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=8, loc='lower left')
    fig.tight_layout()
    return fig


def save_data(prefix, prof, exp):
    """ASCII dump of everything plotted (for re-plotting in another tool)."""
    for key, (x, y) in prof.items():
        np.savetxt(f'{prefix}_TRANSP_{key}.txt', np.c_[x, y], header=f'rho_tor_norm  {key} [SI: m^-3 / eV]')
    for key, d in exp.items():
        for dd in (d if isinstance(d, list) else [d]):
            if dd is None or dd['rho'].size == 0:
                continue
            tag = dd.get('tag', dd['label']).replace('/', '_').replace(' ', '_')
            np.savetxt(f'{prefix}_{key}_{tag}.txt', np.c_[dd['rho'], dd['y'], dd['err'], dd['R'], dd['Z']],
                       header=f'{dd["label"]}  rho_tor_norm  value  error [SI: m^-3 / eV]  R[m]  Z[m]')
    print(f'Data written to {prefix}_*.txt')


# ----------------------------------------------------------------------------
def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('pulse', type=int)
    p.add_argument('runid', help='TRANSP run id, e.g. M01 (file <pulse><runid>.CDF)')
    p.add_argument('--time', type=float, required=True,
                   help='JET time [s] (a value < 40 is taken as TRANSP time)')
    p.add_argument('--twin', type=float, default=0.0,
                   help='averaging window [s] centred on --time for TRANSP and measured data (default 0: nearest slice)')
    p.add_argument('--cdf', default=None, help='explicit path to the TRANSP CDF (skips the search)')

    g = p.add_argument_group('HRTS (ne, Te)')
    g.add_argument('--no-hrts', action='store_true')
    g.add_argument('--hrts-dda', default='HRTS')
    g.add_argument('--hrts-uid', default='jetppf')
    g.add_argument('--hrts-seq', type=int, default=0)
    g.add_argument('--hrts-ne-range', type=float, nargs=2, default=[0.0, 1.0], metavar=('MIN', 'MAX'),
                   help='rho range of HRTS ne (default: whole profile 0-1)')
    g.add_argument('--hrts-te-range', '--hrts-range', type=float, nargs=2, default=[0.65, 1.0], metavar=('MIN', 'MAX'),
                   help='rho range of HRTS Te (default 0.65-1.0, complements KK3)')

    g = p.add_argument_group('KK3 ECE (Te)')
    g.add_argument('--no-kk3', action='store_true')
    g.add_argument('--kk3-dda', default='KK3')
    g.add_argument('--kk3-uid', default='jetppf')
    g.add_argument('--kk3-seq', type=int, default=0)
    g.add_argument('--kk3-range', type=float, nargs=2, default=[0.0, 0.7], metavar=('MIN', 'MAX'))
    g.add_argument('--kk3-side', choices=['lfs', 'hfs', 'both'], default='lfs')

    g = p.add_argument_group('core CX (Ti) -- off by default')
    g.add_argument('--cx', nargs='*', default=[], metavar='DDA', help='e.g. --cx CXG6 CXD6')
    g.add_argument('--cx-uid', default='jetppf')
    g.add_argument('--cx-seq', type=int, default=0)
    g.add_argument('--cx-dty', default='TI')
    g.add_argument('--cx-range', type=float, nargs=2, default=[0.0, 1.0], metavar=('MIN', 'MAX'))

    g = p.add_argument_group('output')
    g.add_argument('-q', '--q', action='store_true', help='add panel (a) with the q profile on top')
    g.add_argument('--q-source', choices=['efit', 'transp'], default='efit',
                   help='q profile from EFIT PPF (default) or TRANSP Q')
    g.add_argument('--q-dda', default='EFTP')
    g.add_argument('--q-uid', default='jetppf')
    g.add_argument('--q-seq', type=int, default=0)
    g.add_argument('--xlim', type=float, nargs=2, default=None, metavar=('MIN', 'MAX'), help='r/a range (default 0 1)')
    g.add_argument('--nlim', type=float, nargs=2, default=None, metavar=('MIN', 'MAX'), help='density axis [1e20 m^-3]')
    g.add_argument('--tlim', type=float, nargs=2, default=None, metavar=('MIN', 'MAX'), help='temperature axis [keV]')
    g.add_argument('--qlim', type=float, nargs=2, default=None, metavar=('MIN', 'MAX'), help='q axis (with --q)')
    g.add_argument('--out', default=None, help='output file prefix (default <pulse><runid>_<time>)')
    g.add_argument('--save', action='store_true', help='save the figure(s) to files (default: only show)')
    g.add_argument('--formats', nargs='+', default=['png', 'pdf'], help='file formats used with --save')
    g.add_argument('--no-title', action='store_true')
    g.add_argument('--box', action='store_true', help='closed frame with inward ticks (default: reference style)')
    g.add_argument('--pub', action='store_true',
                   help='publication labels: HRTS / ECE / CXRS instead of DDA/uid/seq, and no title')
    g.add_argument('--hrts-label', default=None, help='custom legend label for HRTS (overrides --pub)')
    g.add_argument('--kk3-label', default=None, help='custom legend label for KK3 (overrides --pub)')
    g.add_argument('--cx-label', default=None, help='custom legend label for CX (overrides --pub)')
    g.add_argument('--save-data', action='store_true', help='write plotted data as ASCII')
    g.add_argument('--plot-geometry', action='store_true', help='extra figure: flux surfaces + measurement positions')
    g.add_argument('--no-show', action='store_true', help='do not open the plot window (e.g. batch use with --save)')
    return p.parse_args(argv)


def main(argv=None):
    a = parse_args(argv)
    if a.time < JET_T_OFFSET:
        print(f'--time {a.time} < 40: interpreted as TRANSP time (JET {a.time + JET_T_OFFSET:.3f} s)')
        a.time += JET_T_OFFSET
    t_jet = a.time
    t_tr = t_jet - JET_T_OFFSET

    cdf = a.cdf or find_run_files(a.pulse, a.runid)[0]
    run = TranspRun(cdf)
    it = int(np.argmin(np.abs(run.t - t_tr)))
    print(f'Closest TRANSP time slice: {run.t[it] + JET_T_OFFSET:.4f} s (JET)')

    # TRANSP profiles in SI (NE, NI: cm^-3 -> m^-3; TE, TI: eV)
    prof = {}
    for key in ('NE', 'NI', 'TE', 'TI'):
        x, y, tt = run.profile(key, t_tr, a.twin)
        u = run.units(key).upper()
        if key in ('NE', 'NI') and ('CM' in u or u == ''):  # TRANSP: N/CM**3
            y = y * 1e6
        prof[key] = (x, y)
    if a.twin > 0:
        print(f'TRANSP profiles averaged over {tt.size} slices in [{tt[0] + 40:.3f}, {tt[-1] + 40:.3f}] s')
    if a.q:
        efit = None
        if a.q_source == 'efit':
            try:
                efit = read_efit_q(a.pulse, t_jet, a.twin, run, dda=a.q_dda, uid=a.q_uid, seq=a.q_seq)
            except ImportError as e:
                warn(str(e))
            if efit is None:
                warn('EFIT q not available: using TRANSP Q instead.')
        if efit is not None:
            prof['q'] = (efit[0], efit[1])
        else:
            x, y, _ = run.profile('Q', t_tr, a.twin)
            prof['q'] = (x, y)
            print('q profile: TRANSP Q')

    eq = TranspEquilibrium(run)

    exp = {}
    if not a.no_hrts:
        for dty in ('NE', 'TE'):
            d = read_hrts(a.pulse, dty, a.hrts_uid, a.hrts_seq, t_jet, a.twin, eq, dda=a.hrts_dda)
            exp[f'HRTS_{dty}'] = d
    if not a.no_kk3:
        exp['KK3_TE'] = read_kk3(a.pulse, a.kk3_uid, a.kk3_seq, t_jet, a.twin, eq, side=a.kk3_side, dda=a.kk3_dda)
    if a.cx:
        exp['CX_TI'] = [read_cx_core(a.pulse, dda, a.cx_uid, a.cx_seq, t_jet, a.twin, eq, dty=a.cx_dty) for dda in a.cx]

    exp_full = {k: v for k, v in exp.items() if not isinstance(v, list)}
    rng = {'HRTS_NE': a.hrts_ne_range, 'HRTS_TE': a.hrts_te_range, 'KK3_TE': a.kk3_range}
    for k in list(exp):
        if k == 'CX_TI':
            exp[k] = [select_range(d, a.cx_range) for d in exp[k]]
        else:
            exp[k] = select_range(exp[k], rng[k])

    # legend labels: DDA/uid/seq by default, short names with --pub, or user-defined
    pub = {'HRTS': 'HRTS', 'KK3': 'ECE', 'CX': 'CXRS'}
    custom = {'HRTS': a.hrts_label, 'KK3': a.kk3_label, 'CX': a.cx_label}
    for k, v in exp.items():
        kind = k.split('_')[0]
        items = v if isinstance(v, list) else [v]
        n_items = sum(d is not None for d in items)
        for d in items:
            if d is None:
                continue
            d['tag'] = d['label']  # keeps DDA/uid/seq for --save-data file names
            new = custom[kind] or (pub[kind] if a.pub else None)
            if new is not None:
                # several CX systems: keep them distinguishable
                d['label'] = f"{new} ({d['tag'].split('/')[0]})" if n_items > 1 else new

    title = None if (a.no_title or a.pub) else f'{a.pulse}{a.runid} @ {run.t[it] + 40:.3f} s' + \
        (f' (avg {a.twin:.3g} s)' if a.twin > 0 else '')
    prefix = a.out or f'{a.pulse}{a.runid}_{t_jet:.3f}s'
    fig = make_figure(prof, exp, title=title, with_q=a.q,
                      limits={'x': a.xlim, 'n': a.nlim, 'T': a.tlim, 'q': a.qlim}, box=a.box)
    if a.save:
        for fmt in a.formats:
            fig.savefig(f'{prefix}.{fmt}', dpi=300)
            print(f'Saved {prefix}.{fmt}')

    if a.plot_geometry:
        fg = plot_geometry(eq, t_tr, {k: v for k, v in exp_full.items() if v is not None},
                           title=f'{a.pulse}{a.runid} @ {t_jet:.3f} s')
        if a.save:
            fg.savefig(f'{prefix}_geometry.png', dpi=200)
            print(f'Saved {prefix}_geometry.png')
    if a.save_data:
        save_data(prefix, prof, exp)
    if not a.no_show:
        plt.show()


if __name__ == '__main__':
    main()
