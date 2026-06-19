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
                      max_pts=40000, axis_tol=0.05):
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

    # Where does the magnetic-axis *ring* (R=Rmag at every toroidal angle) pass
    # through the LOS? By axisymmetry the ring's 3-D distance to a cell equals
    # the poloidal distance hypot(R-Rmag, Z-Zmag). The cells within ``axis_tol``
    # of the ring are the ones the axis threads through -- recolour them so the
    # passage (the chord goes in one side of the axis and out the other) is
    # visible as a coloured cluster, rather than the chord looking like it runs
    # *along* the axis. The closest cell (the actual grazing point) is reported.
    d_axis = np.hypot(Rc - Rmag, Zc - Zmag)
    on_axis = ins & (Cc > 0) & (d_axis <= axis_tol)
    dmin_cm = float(d_axis[ins & (Cc > 0)].min()) * 100.0 if ins.any() else float('nan')
    on_lbl = (f'on axis ($\\leq${axis_tol * 100:.0f} cm; min {dmin_cm:.1f} cm)'
              if on_axis.any() else None)

    def _mark_axis(ax, hx, hy):
        """Recolour the cells the axis ring threads through (bright red)."""
        if on_axis.any():
            ax.scatter(hx[on_axis], hy[on_axis], c='red', s=16, alpha=0.9,
                       linewidths=0, label=on_lbl, zorder=5)

    fig, axes = plt.subplots(1, 3, figsize=(19.0, 5.6))
    if title:
        fig.suptitle(title)
    th = np.linspace(0.0, 2.0 * np.pi, 361)

    # ---- (0) poloidal (R, Z) ---------------------------------------------
    ax = axes[0]
    ax.plot(Rb, Zb, 'b-', lw=1.4, label='LCFS')
    _scatter_cells(ax, fig, Rc, Zc, Cc, inn, out)
    ax.plot([Rmag], [Zmag], 'k+', ms=11, label=f'axis ({Rmag:.2f}, {Zmag:.2f})')
    _mark_axis(ax, Rc, Zc)
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
    _mark_axis(ax, x, y)
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
    # Axis height is only a *reference* level (the ring sits at Z=Zmag for every
    # phi); the chord is NOT along the axis -- it merely runs near this height.
    # The actual axis-ring crossings (R=Rmag) are starred, with drop-lines at
    # their toroidal x, so the touch (lime) vs the far crossing (cyan) is clear.
    ax.axhline(Zmag, color='r', ls=':', lw=1.0,
               label=f'axis height Z={Zmag:.2f} m (ref., any $\\phi$)')
    ax.axhspan(float(Zb.min()), float(Zb.max()), color='b', alpha=0.07)
    ax.axhline(float(Zb.min()), color='b', ls='--', lw=0.8,
               label=f'plasma Z extent [{Zb.min():.2f}, {Zb.max():.2f}]')
    ax.axhline(float(Zb.max()), color='b', ls='--', lw=0.8)
    _mark_axis(ax, x, Zc)
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
