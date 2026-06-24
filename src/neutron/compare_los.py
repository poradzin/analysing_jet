#!/usr/bin/env python3
"""compare_los.py -- overlay the KM9 and KM14 **real-LOS** TH/BT diagnostics.

Both ``km9/los_th_bt_ratio.py`` and ``km14/los_th_bt_ratio.py`` reduce, for the
*real* line of sight, to the identical machinery: the per-shell detector
coupling ``f_det = C_bin/DVOL`` from the diagnostic's own ``.los`` cell cloud
(``los_common.los_file_detector_rate``), applied to the thermal (THNTX) and the
flux-averaged beam-target (BTN) emissivity to form the detector-weighted
cumulative ``(TH/BT)_LOS = Sum(TH f_det DVOL)/Sum(BT f_det DVOL)``. The two
scripts differ only in the ``.los`` file and the (auto-detected) chord
orientation. This driver recomputes that **real-LOS** path for each diagnostic
via the shared ``los_common`` library -- so the numbers reproduce the standalone
scripts exactly -- and overlays them in a single 1 x (2 + N) figure (a 1x4 for
KM9 + KM14):

  (0) Detector-coupling weight  f_det(rhot), normalized so int(f_det DVOL)=1
                                (``normalized_coupling_weight`` -- the arbitrary
                                etendue scale differs between KM9 and KM14, so
                                the *normalized* weight is what's comparable),
                                plus (right axis, dotted) the cumulative fraction
                                of the detector-weighted BT signal banked by each
                                rhot -- the core/edge split that explains why the
                                panel-1 ratios diverge (see doc/compare_los.md)
  (1) TH/BT vs rhot             cumulative real-LOS C-weighted TH/BT per
                                diagnostic, over the shared full-plasma (0D)
                                cumulative and local TH(rhot)/BT(rhot)
  (2..) Per-shell detector signal  TH*f_det*DVOL / BT*f_det*DVOL -- **one panel
                                per diagnostic** (own y-axis): the absolute
                                detector-reaching rate scales with each
                                diagnostic's etendue, which differs by ~the
                                detector solid angle between KM9 and KM14, so a
                                shared axis would bury the smaller signal.

No idealised box chord is used for KM14 here -- only the real LOS, so KM9 and
KM14 are compared on the same footing.

The full-plasma cumulative and the local TH/BT depend only on the emissivity
profiles and DVOL (no LOS), so they are identical for both diagnostics at the
same pulse/run/idx/channel and are drawn once as a shared grey reference.

CLI
---
    python compare_los.py 104614 M29 --idx 2
    python compare_los.py 104614 M29 --idx 2 --channel dd --save cmp.png
    python compare_los.py 104614 M29 --idx 2 \
        --km9-los src/neutron/km9/KM9.los --km14-los src/neutron/km14/KM3.los
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

# Shared LOS / equilibrium / binning library (src/neutron/common/) and the
# NUBEAM _fi/_neut reader (currently in km14/ -- see CONTEXT Stage-B follow-ups).
_HERE = os.path.dirname(os.path.abspath(__file__))
for _sub in ("common", "km14"):
    _d = os.path.join(_HERE, _sub)
    if _d not in sys.path:
        sys.path.insert(0, _d)

from los_common import (
    find_run_dir,
    read_thermal_slice,
    read_los_file,
    flux_avg_profile,
    interp_flux_to,
    EqCDF,
    los_file_detector_rate,
    normalized_coupling_weight,
)
import bt_zone_integrator as bzi

# Default real LOS cell files (sit next to each diagnostic's script).
DEFAULT_KM9_LOS = os.path.join(_HERE, "km9", "KM9.los")
DEFAULT_KM14_LOS = os.path.join(_HERE, "km14", "KM3.los")

# (thermal profile var, list of beam-target _neut keys to sum, label).
# Empty bt-key list => sum every BT component present in _neut.
CHANNELS = {
    "total": ("THNTX",    [],      "total (DD+DT+...)"),
    "dd":    ("THNTX_DD", ["DD"],  "DD (2.45 MeV)"),
    "dt":    ("THNTX_DT", ["DT"],  "DT (14 MeV)"),
}


def real_los_profiles(pulse, runid, los_file, name, idx=None, channel="total",
                      data_dir=None, chord_kind="auto", subgrid=True,
                      median=True):
    """Real-LOS TH/BT detector profiles for one diagnostic.

    Generalises the ``km9/los_th_bt_ratio.py`` ``main()`` over an arbitrary
    ``.los`` file (KM9's path *is* the pure real-LOS path; KM14's ``--los-file``
    branch was built to mirror it), so the same routine serves both diagnostics
    and the results reproduce each standalone script.

    Returns a comparison dict on the sorted TRANSP rhot grid ``xs``:
      name, label, pulse, runid, idx, t_transp, chan_label
      xs, f_det, f_det_norm                  detector coupling weight
      th_si, bt_si                           emissivities [n/m3/s]
      dvol_m3                                shell volume [m3]
      ratio_local                            TH(rhot)/BT(rhot)            (shared)
      ratio_cum_det                          cum. real-LOS C-weighted TH/BT
      ratio_cum_plasma                       cum. full-plasma 0D TH/BT    (shared)
      shell_th_det, shell_bt_det             per-shell TH*f_det*DVOL etc. [n/s]
      rate_det_th, rate_det_bt, thbt_los     detector totals / endpoint
      rho_med_th, rho_med_bt                 signal-median rhot
      rhot_min, rhot_crit                    chord resolution markers
      n_inside, n_cells, chord_kind          geometry diagnostics
    """
    run_id = f"{pulse}{runid}"
    run_dir = find_run_dir(pulse, runid, data_dir)
    cdf_path = run_dir / f"{run_id}.CDF"
    th_var, bt_keys_req, chan_label = CHANNELS[channel]

    idxs = bzi.list_fbm_indices(run_dir, run_id)
    if not idxs:
        raise FileNotFoundError(f"No {run_id}_fi_*.cdf in {run_dir}")
    idx = idx if idx is not None else idxs[0]
    fi_path = run_dir / f"{run_id}_fi_{idx}.cdf"
    neut_path = run_dir / f"{run_id}_neut_{idx}.cdf"

    # --- BT zone data + flux-averaged BT profile -----------------------------
    fi = bzi.read_fi_distribution(fi_path)
    neut = bzi.read_neut_rates(neut_path)
    available = [k for k in ("DD", "DT", "TT", "TD") if k in neut]
    bt_keys = bt_keys_req if bt_keys_req else available
    missing = [k for k in bt_keys if k not in neut]
    if missing:
        raise KeyError(f"Beam-target key(s) {missing} not in {neut_path}")
    bt_zone_sum = sum(neut[k] for k in bt_keys)
    x_bt, bt_prof = flux_avg_profile(bt_zone_sum, fi["x2d"], fi["bmvol"])

    # --- Equilibrium (CDF) at the fast-ion output time -----------------------
    eqs = EqCDF(cdf_path, fi["time"] + 40.0)   # EqCDF wants JET time

    # --- Thermal profile + DVOL on the TRANSP X grid -------------------------
    (thntx_x, x_rhot, dvol_x, th_to_si,
     dvol_to_m3, _unit) = read_thermal_slice(cdf_path, th_var, eqs.tind)
    order = np.argsort(x_rhot)
    xs = np.asarray(x_rhot)[order]
    th_native = np.asarray(thntx_x)[order]                 # native CGS
    dvol_m3 = np.asarray(dvol_x)[order] * dvol_to_m3
    th_si = th_native * th_to_si                           # n/m3/s
    bt_native = interp_flux_to(x_bt, bt_prof, xs)          # 1/cm3/s on xs
    bt_si = bt_native * 1.0e6                              # n/m3/s

    # --- Real LOS cells: f_det (geometric) + per-shell detector signals ------
    cells = read_los_file(los_file)
    common = dict(subgrid=subgrid, chord_kind=chord_kind, median=median)
    fTH = los_file_detector_rate(cells, eqs, xs, th_native, th_to_si,
                                 dvol_m3, **common)
    fBT = los_file_detector_rate(cells, eqs, xs, bt_native, 1.0e6,
                                 dvol_m3, **common)
    f_det = fTH["f_det"]
    assert np.allclose(f_det, fBT["f_det"]), "f_det mismatch (should be geometric)"

    shell_th_det = th_si * f_det * dvol_m3                 # n/s per shell
    shell_bt_det = bt_si * f_det * dvol_m3
    cum_th_det = np.cumsum(shell_th_det)
    cum_bt_det = np.cumsum(shell_bt_det)
    cum_th_plasma = np.cumsum(th_si * dvol_m3)
    cum_bt_plasma = np.cumsum(bt_si * dvol_m3)
    with np.errstate(invalid="ignore", divide="ignore"):
        ratio_local = np.where(bt_si > 0, th_si / bt_si, np.nan)
        ratio_cum_det = np.where(cum_bt_det > 0, cum_th_det / cum_bt_det, np.nan)
        ratio_cum_plasma = np.where(cum_bt_plasma > 0,
                                    cum_th_plasma / cum_bt_plasma, np.nan)
    thbt_los = fTH["rate_det"] / fBT["rate_det"]           # = ratio_cum_det[-1]
    # cumulative fraction of the detector-weighted BT signal (the cumulative
    # ratio's denominator weight) -- the "ballast" curve: how much of the LOS
    # signal each diagnostic has banked by a given rhot. Its core vs edge split
    # is what drives the divergence of ratio_cum_det between KM9 and KM14 (see
    # doc/compare_los.md). Normalized to 1 at rhot=1 so the two are comparable.
    cum_bt_frac = cum_bt_det / cum_bt_det[-1] if cum_bt_det[-1] > 0 \
        else np.zeros_like(cum_bt_det)

    return dict(
        name=name, label=f"{name} ({os.path.basename(los_file)})",
        pulse=pulse, runid=runid, idx=idx, t_transp=fi["time"],
        chan_label=chan_label,
        xs=xs, f_det=f_det,
        f_det_norm=normalized_coupling_weight(f_det, dvol_m3, xs),
        th_si=th_si, bt_si=bt_si, dvol_m3=dvol_m3,
        ratio_local=ratio_local, ratio_cum_det=ratio_cum_det,
        ratio_cum_plasma=ratio_cum_plasma,
        shell_th_det=shell_th_det, shell_bt_det=shell_bt_det,
        cum_bt_frac=cum_bt_frac,
        rate_det_th=fTH["rate_det"], rate_det_bt=fBT["rate_det"],
        thbt_los=thbt_los,
        rho_med_th=fTH["rho_med"], rho_med_bt=fBT["rho_med"],
        rhot_min=fTH["rhot_min"], rhot_crit=fTH["rhot_crit"],
        n_inside=fTH["n_inside"], n_cells=int(cells["R"].size),
        chord_kind=fTH["chord_kind"],
    )


def compare_plot(results, save=None, show=True):
    """Overlay N real-LOS diagnostics: two shared panels (normalized weight,
    TH/BT vs rhot) followed by one **dedicated per-shell detector-signal panel
    per diagnostic**. The per-shell signal carries each diagnostic's absolute
    etendue scale, which differs by ~the detector solid angle between KM9 and
    KM14, so overlaying them buries the smaller one -- hence one panel each, with
    its own y-axis. Layout is therefore 1 x (2 + N) (a 1x4 for KM9 + KM14)."""
    import matplotlib.pyplot as plt

    colors = ["C3", "C0", "C2", "C4"]   # one per diagnostic (KM9, KM14, ...)
    n_panels = 2 + len(results)
    fig, ax = plt.subplots(1, n_panels, figsize=(5.5 * n_panels, 5.6))
    ref = results[0]
    fig.suptitle(
        f"KM9 vs KM14 real-LOS TH/BT   pulse {ref['pulse']}  "
        f"TRANSP {ref['runid']}   idx {ref['idx']}  "
        f"t={ref['t_transp']:.3f}s   {ref['chan_label']}")

    # ---- (0) detector-coupling weight (normalized) + cumulative signal ------
    # Solid = normalized f_det (left axis). Faint dotted = cumulative fraction of
    # the detector-weighted BT signal (right axis, 0->1): the "ballast" curve.
    # The two diagnostics' normalized weights nearly overlap for rhot > ~0.4, yet
    # they bank different fractions of the signal in the core -- and that core/
    # edge split, leveraged by the steep edge falloff of the local TH/BT, is what
    # drives the cumulative-ratio divergence in panel 1 (see doc/compare_los.md).
    a = ax[0]
    a2 = a.twinx()
    for r, c in zip(results, colors):
        a.plot(r["xs"], r["f_det_norm"], "-", color=c, lw=1.5, label=r["label"])
        a.fill_between(r["xs"], 0.0, r["f_det_norm"], color=c, alpha=0.10)
        a2.plot(r["xs"], r["cum_bt_frac"], ":", color=c, lw=1.3, alpha=0.9)
        rc = r["rhot_crit"] if r["rhot_crit"] is not None else r["rhot_min"]
        if rc is not None and np.isfinite(rc):
            a.axvline(rc, color=c, ls=":", lw=1.0)
    a.set_xlim(0, 1); a.set_ylim(bottom=0.0); a.grid(True, ls=":", lw=0.5)
    a2.set_ylim(0.0, 1.0)
    a2.set_ylabel(r"cum. BT signal fraction $\sum$BT$\,f_{\rm det}$DVOL "
                  "(dotted)", fontsize=9)
    a.set_xlabel("rhot")
    a.set_ylabel(r"$f_{\rm det}/\!\int\! f_{\rm det}\,$DVOL$\,d\rho_t$ (solid) "
                 r"[m$^{-3}$]")
    a.set_title("Detector-coupling weight + cumulative signal")
    a.legend(loc="upper right", fontsize=8)

    # ---- (1) TH/BT vs rhot --------------------------------------------------
    # full-plasma cumulative and the local ratio depend only on the emissivity
    # profiles + DVOL (no LOS), so they are shared between diagnostics -- drawn
    # once from the first result as a grey reference.
    a = ax[1]
    a.plot(ref["xs"], ref["ratio_local"], color="0.55", lw=1.0, alpha=0.9,
           label=r"local TH$(\rho_t)$/BT$(\rho_t)$")
    a.plot(ref["xs"], ref["ratio_cum_plasma"], color="0.35", ls="--", lw=1.3,
           label=f"cum. full plasma (rhot=1: {ref['ratio_cum_plasma'][-1]:.3f})")
    for r, c in zip(results, colors):
        a.plot(r["xs"], r["ratio_cum_det"], "-", color=c, lw=1.7,
               label=f"{r['name']} cum. real-LOS (rhot=1: {r['thbt_los']:.3f})")
    a.set_xlim(0, 1); a.set_ylim(bottom=0.0); a.grid(True, ls=":", lw=0.5)
    a.set_xlabel("rhot"); a.set_ylabel("TH / BT")
    a.set_title("TH/BT vs rhot  (cumulative real-LOS & local)")
    a.legend(loc="best", fontsize=8)

    # ---- (2..) per-shell detector signal -- one panel per diagnostic --------
    # Each diagnostic gets its own panel/y-axis: the absolute detector-reaching
    # rate scales with the (very different) KM9 vs KM14 etendue, so a shared axis
    # would render the smaller signal invisible. TH solid, BT dashed; colour =
    # diagnostic (matching panels 0-1).
    for i, (r, c) in enumerate(zip(results, colors)):
        a = ax[2 + i]
        a.plot(r["xs"], r["shell_th_det"], "-", color=c, lw=1.5,
               label=r"TH$\cdot f_{\rm det}\cdot$DVOL")
        a.plot(r["xs"], r["shell_bt_det"], "--", color=c, lw=1.5,
               label=r"BT$\cdot f_{\rm det}\cdot$DVOL")
        a.fill_between(r["xs"], 0.0, r["shell_th_det"], color=c, alpha=0.08)
        a.set_xlim(0, 1); a.set_ylim(bottom=0.0); a.grid(True, ls=":", lw=0.5)
        a.set_xlabel("rhot"); a.set_ylabel("per-shell LOS rate [n/s]")
        a.set_title(f"{r['name']} per-shell detector signal\n"
                    f"($\\int$: TH={r['rate_det_th']:.2e}, "
                    f"BT={r['rate_det_bt']:.2e} n/s)")
        a.legend(loc="best", fontsize=8)

    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=130)
        print(f"Saved comparison figure to {save}")
    if show:
        plt.show()
    return fig


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("pulse", type=int, help="JET pulse number")
    p.add_argument("runid", type=str, help="TRANSP run suffix (e.g. M29)")
    p.add_argument("--idx", type=int, default=None,
                   help="Fast-ion output index (1-based). Default: first available.")
    p.add_argument("--channel", choices=list(CHANNELS), default="total",
                   help="Channel: total (default), dd, dt.")
    p.add_argument("--data-dir", default=None,
                   help="Base data dir (<base>/<pulse>/<runid>).")
    p.add_argument("--km9-los", default=DEFAULT_KM9_LOS,
                   help=f"Real KM9 LOS cell file (default: {DEFAULT_KM9_LOS}).")
    p.add_argument("--km14-los", default=DEFAULT_KM14_LOS,
                   help=f"Real KM14 LOS cell file (default: {DEFAULT_KM14_LOS}).")
    p.add_argument("--save", default=None, help="Save the figure to this path.")
    p.add_argument("--no-plot", action="store_true",
                   help="Suppress the figure even with --save.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    common = dict(idx=args.idx, channel=args.channel, data_dir=args.data_dir)
    diags = [("KM9", args.km9_los), ("KM14", args.km14_los)]

    results = []
    for name, los in diags:
        r = real_los_profiles(args.pulse, args.runid, los, name, **common)
        results.append(r)
        print(f"{name:5s}  los={os.path.basename(los)}  "
              f"chord={r['chord_kind']}  cells_inside={r['n_inside']}/{r['n_cells']}  "
              f"(TH/BT)_LOS={r['thbt_los']:.4f}  "
              f"rho50(TH/BT)={r['rho_med_th']:.3f}/{r['rho_med_bt']:.3f}")
    print(f"  full-plasma 0D (TH/BT) = {results[0]['ratio_cum_plasma'][-1]:.4f}")

    if args.no_plot and args.save:
        import matplotlib
        matplotlib.use("Agg")               # headless: write the file, no window
        compare_plot(results, save=args.save, show=False)
    elif not args.no_plot:
        compare_plot(results, save=args.save)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
