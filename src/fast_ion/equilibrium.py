"""Minimal equilibrium helper: LCFS reconstruction from a TRANSP main CDF.

``boundary_from_moments`` is copied verbatim (renamed from the private
``_boundary_from_moments``) from ``src/neutron/common/los_common.py`` -- it
has no other dependency on that library, so it's duplicated here rather than
imported, keeping this directory self-contained and shareable on its own
(no ppf, no rest of analysing_jet needed).
"""
from __future__ import annotations

import numpy as np
from netCDF4 import Dataset


def boundary_from_moments(d, tind, ntheta=257):
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


def load_lcfs(cdf_path, tind):
    """Last-closed-flux-surface (R, Z) [m] at TIME3 index ``tind`` of ``cdf_path``."""
    with Dataset(cdf_path, "r") as d:
        return boundary_from_moments(d, tind)
