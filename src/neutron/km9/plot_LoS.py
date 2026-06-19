#!/usr/bin/env python
"""
plot_LoS.py
-----------

Geometry views of a real line-of-sight (LOS) cell cloud, shared by the KM9
``los_thermal_rate.py`` and ``los_th_bt_ratio.py`` scripts (called optionally
via ``--plot-los``).

It draws three orthogonal projections of the LOS cells (coloured by the etendue
weight ``log10 C``, out-of-LCFS cells dimmed), which together show how a
semi-tangential horizontal chord like KM9 only *grazes* the magnetic axis:

  1. **Poloidal (R, Z)** -- the usual flux-coordinate cross-section, with the
     magnetic axis and the LCFS. The tangential chord's inbound and outbound
     legs both fall here, so the cloud brackets the axis.
  2. **Top (x, y)** -- the toroidal plane viewed from above, with the magnetic
     axis as the circle ``R = Rmag`` and the LCFS annulus. A tangential chord is
     a near-straight secant that cuts the axis circle at two toroidal points.
  3. **Side elevation (x, Z)** -- the torus viewed from the side (looking
     radially), the analogue of Fig. 9 (top) of Andersson Sunden et al., NIM A
     610 (2009) 682: the chord is a near-horizontal band; the magnetic-axis
     *height* ``Z = Zmag`` and the plasma vertical extent are drawn, so one sees
     the chord skims the axis height only over a short toroidal span.

The point of closest approach to the magnetic axis (min poloidal distance over
in-LCFS cells) is starred in the poloidal and side views.
"""

import numpy as np
import matplotlib.pyplot as plt


def _scatter_cells(ax, fig, hx, hy, Cc, inn, out, clabel=True):
    """Scatter in-LCFS cells (coloured by log10 C) and dimmed out-of-LCFS cells."""
    ax.scatter(hx[out], hy[out], c='0.8', s=2, alpha=0.4, linewidths=0,
               label='cells outside LCFS')
    sc = ax.scatter(hx[inn], hy[inn], c=np.log10(Cc[inn]), s=4, cmap='plasma',
                    alpha=0.7, linewidths=0, label='cells inside LCFS')
    if clabel:
        cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
        cb.set_label(r'$\log_{10}$ C  [m$^3$]  (etendue)')
    return sc


def plot_los_geometry(cells, inside_cells, Rb, Zb, Rmag, Zmag,
                      title="", rhot_min=None, save=None, show=True,
                      max_pts=40000):
    """Three-projection geometry figure of a real LOS cell cloud.

    Parameters
    ----------
    cells        : dict from ``los_common.read_los_file`` (x, y, z, R, Z, C, ...)
    inside_cells : bool array (in-LCFS mask), aligned with ``cells``
    Rb, Zb       : LCFS polygon [m]
    Rmag, Zmag   : magnetic axis [m]
    title        : figure suptitle
    rhot_min     : optional, annotated as the closest-approach rhot
    save         : path to save the figure (optional)
    show         : call plt.show() (default True)
    """
    x = np.asarray(cells['x']); y = np.asarray(cells['y'])
    Rc = np.asarray(cells['R']); Zc = np.asarray(cells['Z'])
    Cc = np.asarray(cells['C']); ins = np.asarray(inside_cells)
    Rb = np.asarray(Rb); Zb = np.asarray(Zb)

    idx = np.where(Cc > 0)[0]
    if idx.size > max_pts:
        idx = idx[np.linspace(0, idx.size - 1, max_pts).astype(int)]
    inn = idx[ins[idx]]
    out = idx[~ins[idx]]

    # Closest-approach cell to the axis (poloidal distance), over in-LCFS cells.
    star = None
    if ins.any():
        d = np.hypot(Rc - Rmag, Zc - Zmag)
        d = np.where(ins & (Cc > 0), d, np.inf)
        j = int(np.argmin(d))
        if np.isfinite(d[j]):
            star = (j, float(d[j]))

    fig, axes = plt.subplots(1, 3, figsize=(19.0, 5.6))
    if title:
        fig.suptitle(title)
    th = np.linspace(0.0, 2.0 * np.pi, 361)

    # ---- (0) poloidal (R, Z) ---------------------------------------------
    ax = axes[0]
    ax.plot(Rb, Zb, 'b-', lw=1.4, label='LCFS')
    _scatter_cells(ax, fig, Rc, Zc, Cc, inn, out)
    ax.plot([Rmag], [Zmag], 'r+', ms=11, label=f'axis ({Rmag:.2f}, {Zmag:.2f})')
    if star is not None:
        ax.plot([Rc[star[0]]], [Zc[star[0]]], marker='*', ms=13, mfc='lime',
                mec='k', mew=0.6, ls='none',
                label=f'closest approach ({star[1] * 100:.1f} cm)')
    ax.set_xlabel('R [m]'); ax.set_ylabel('Z [m]')
    ax.set_title('Poloidal (R, Z)')
    ax.legend(loc='best', fontsize=8); ax.grid(True, ls=':', lw=0.5)
    ax.set_aspect('equal', adjustable='box')

    # ---- (1) top (x, y) ---------------------------------------------------
    ax = axes[1]
    _scatter_cells(ax, fig, x, y, Cc, inn, out)
    ax.plot(Rmag * np.cos(th), Rmag * np.sin(th), 'r-', lw=1.4,
            label=f'magnetic axis  R={Rmag:.2f} m')
    ax.plot(Rb.min() * np.cos(th), Rb.min() * np.sin(th), 'b--', lw=0.9,
            label=f'LCFS R_in={Rb.min():.2f}')
    ax.plot(Rb.max() * np.cos(th), Rb.max() * np.sin(th), 'b--', lw=0.9,
            label=f'LCFS R_out={Rb.max():.2f}')
    ax.plot([0.0], [0.0], 'k+', ms=9, label='machine axis')
    ax.set_xlabel('x [m]  (toroidal)'); ax.set_ylabel('y [m]  (radial)')
    ax.set_title('Top view (x, y): LOS vs axis circle')
    ax.legend(loc='upper right', fontsize=7); ax.grid(True, ls=':', lw=0.5)
    ax.set_aspect('equal', adjustable='box')

    # ---- (2) side elevation (x, Z) ---------------------------------------
    # Torus viewed from the side (looking radially, along y). The magnetic axis
    # ring projects to the horizontal line Z = Zmag; the plasma vertical extent
    # is [Zb.min, Zb.max]; the LCFS toroidal half-width at the tangency is
    # R_out, so the plasma silhouette spans x in [-R_out, R_out].
    ax = axes[2]
    _scatter_cells(ax, fig, x, Zc, Cc, inn, out)
    ax.axhline(Zmag, color='r', lw=1.4, label=f'axis height Z={Zmag:.2f} m')
    ax.axhspan(float(Zb.min()), float(Zb.max()), color='b', alpha=0.07)
    ax.axhline(float(Zb.min()), color='b', ls='--', lw=0.8,
               label=f'plasma Z extent [{Zb.min():.2f}, {Zb.max():.2f}]')
    ax.axhline(float(Zb.max()), color='b', ls='--', lw=0.8)
    ax.axvline(-float(Rb.max()), color='0.6', ls=':', lw=0.8)
    ax.axvline(float(Rb.max()), color='0.6', ls=':', lw=0.8,
               label=f'|x| = R_out = {Rb.max():.2f} m')
    if star is not None:
        ax.plot([x[star[0]]], [Zc[star[0]]], marker='*', ms=13, mfc='lime',
                mec='k', mew=0.6, ls='none', label='closest approach')
    ax.set_xlabel('x [m]  (toroidal)'); ax.set_ylabel('Z [m]')
    ttl = 'Side view (x, Z): chord vs axis height'
    if rhot_min is not None and np.isfinite(rhot_min):
        ttl += f'   (min sampled rhot = {rhot_min:.4f})'
    ax.set_title(ttl)
    ax.legend(loc='best', fontsize=7); ax.grid(True, ls=':', lw=0.5)
    ax.set_aspect('equal', adjustable='box')

    plt.tight_layout()
    if save:
        fig.savefig(save, dpi=130)
        print(f'Saved LOS geometry figure to {save}')
    if show:
        plt.show()
    return fig
