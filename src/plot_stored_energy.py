#!/usr/bin/env python
"""
Compare TRANSP stored-energy time traces for several runs of one JET pulse.

Two panels:
  top    - plasma energy: EQ/WDIA, EQ/WP and per run
             WDIA (TRANSP) = 1.5*UFASTPP + UTHRM
             WMHD (TRANSP) = 0.75*UFASTPP + 1.5*(UFASTPA + UPHI) + UTHRM
  bottom - plasma thermal energy, per run:
             EQ/WDIA - 1.5*UFASTPP
             EQ/WP   - 0.75*UFASTPP - 1.5*(UFASTPA + UPHI)
             UTHRM (TRANSP)

Signal definitions follow OMFIT modules/TRANSP/PLOTS/JET/plot_summar_JET.py
("Diamagnetic Energy" and "Thermal Energy" plots). UPHI is used only when
present in the run.

Usage:
  python plot_stored_energy.py 99173 M01 M02 M03
  python plot_stored_energy.py 99173 M01 M02 --eq EFTP --eq-uid myuid --eq-seq 12
"""

import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np

import plotWindow as pw
import profiles as ps

try:
    import ppf
except ImportError:
    ppf = None


ENERGY_SIGNALS = ["UTHRM", "UFASTPP", "UFASTPA", "UPHI"]
ENERGY_SCALE = 1e-6  # J -> MJ
ENERGY_UNIT = "MJ"


def cornernote(pulse):
    plt.gcf().text(0.98, 0.01, str(pulse), fontsize=9, ha="right", va="bottom")


def volume_integral(transp, signal):
    """Volume integral of a [J/cm3] profile -> [J] time trace (DVOL in cm3)."""
    data = transp.transp(signal)
    if data is None:
        return None
    if data.ndim == 1:
        return data
    return np.nansum(data * np.asarray(transp.dvol), axis=1)


def load_run(pulse, runid):
    transp = ps.Transp(pulse, runid)
    if not hasattr(transp, "_dat"):
        print(f"Skipping run {runid}: CDF not found.")
        return None

    for sig in ENERGY_SIGNALS:
        if sig in transp.all_keys:
            transp.add_data(sig)

    missing = [s for s in ("UTHRM", "UFASTPP", "UFASTPA") if s not in transp.all_keys]
    if missing:
        print(f"Skipping run {runid}: missing {', '.join(missing)}.")
        return None

    w = {sig: volume_integral(transp, sig) for sig in ENERGY_SIGNALS if sig in transp.all_keys}
    uphi = w.get("UPHI", np.zeros_like(w["UTHRM"]))
    fast_mhd = 0.75 * w["UFASTPP"] + 1.5 * (w["UFASTPA"] + uphi)

    return {
        "runid": runid,
        "label": transp.transpcid,
        "t": np.asarray(transp.t),
        "uthrm": w["UTHRM"],
        "fast_mhd": fast_mhd,
        "fast_dia": 1.5 * w["UFASTPP"],
        "wmhd": fast_mhd + w["UTHRM"],
        "wdia": 1.5 * w["UFASTPP"] + w["UTHRM"],
        "has_uphi": "UPHI" in w,
    }


def read_eq(pulse, dda, dty, uid, seq):
    """Equilibrium energy time trace -> (t_jet, data) or None."""
    if ppf is None:
        return None
    out = ppf.ppfdata(pulse, dda, dty, seq=seq, uid=uid, device="JET")
    data, t, ier = np.asarray(out[0], dtype=float), np.asarray(out[2], dtype=float), out[-1]
    if ier != 0 or data.size == 0:
        print(f"No {dda}/{dty} (uid={uid}, seq={seq}) for pulse {pulse}.")
        return None
    return t, data


def on_grid(t_new, exp):
    """Interpolate an experimental trace onto TRANSP times (NaN outside its range)."""
    t, y = exp
    return np.interp(t_new, t, y, left=np.nan, right=np.nan)


def plot_stored_energy(pulse, runs, exp, eq):
    """Line style = quantity (dotted: WDIA-based, dashed: WMHD-based, solid: UTHRM);
    colour = source (black: equilibrium, one colour per TRANSP run)."""
    fig, (ax_dia, ax_th) = plt.subplots(2, 1, sharex=True, figsize=(9, 8))
    fig.suptitle(f"JET {pulse}: stored energy", fontsize=13)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    wp, wdia = exp.get("WP"), exp.get("WDIA")

    # Panel 1: plasma (diamagnetic / MHD) energy
    if wdia is not None:
        ax_dia.plot(wdia[0], wdia[1] * ENERGY_SCALE, color="black", ls=":", lw=1.8,
                    label=f"{eq}/WDIA")
    if wp is not None:
        ax_dia.plot(wp[0], wp[1] * ENERGY_SCALE, color="black", ls="--", lw=1.5,
                    label=f"{eq}/WP")

    for i, run in enumerate(runs):
        c = colors[i % len(colors)]
        t, lab = run["t"], run["label"]
        uphi = " + UPHI" if run["has_uphi"] else ""

        ax_dia.plot(t, run["wdia"] * ENERGY_SCALE, color=c, ls=":", lw=2,
                    label=f"{lab} WDIA = 1.5 UFASTPP + UTHRM")
        ax_dia.plot(t, run["wmhd"] * ENERGY_SCALE, color=c, ls="--", lw=2,
                    label=f"{lab} WMHD = 0.75 UFASTPP + 1.5 (UFASTPA{uphi}) + UTHRM")

        # Panel 2: thermal energy
        if wdia is not None:
            w_th = on_grid(t, wdia) - run["fast_dia"]
            ax_th.plot(t, w_th * ENERGY_SCALE, color=c, ls=":", lw=1.8,
                       label=f"{eq}/WDIA - 1.5 UFASTPP [{run['runid']}]")
        if wp is not None:
            w_th = on_grid(t, wp) - run["fast_mhd"]
            ax_th.plot(t, w_th * ENERGY_SCALE, color=c, ls="--", lw=1.8,
                       label=f"{eq}/WP - 0.75 UFASTPP - 1.5 (UFASTPA{uphi}) [{run['runid']}]")
        ax_th.plot(t, run["uthrm"] * ENERGY_SCALE, color=c, lw=2.2, label=f"{lab} UTHRM")

    t_min = min(r["t"][0] for r in runs)
    t_max = max(r["t"][-1] for r in runs)
    ax_th.set_xlim(t_min, t_max)

    ax_dia.set_ylabel(f"Energy [{ENERGY_UNIT}]")
    ax_dia.set_title("Plasma energy (diamagnetic / MHD)", fontsize=11)
    ax_th.set_ylabel(f"Energy [{ENERGY_UNIT}]")
    ax_th.set_title("Plasma thermal energy", fontsize=11)
    ax_th.set_xlabel("time [s] (TRANSP)")
    for ax in (ax_dia, ax_th):
        ax.set_ylim(bottom=0.0)
        ax.grid(alpha=0.3)
        ax.tick_params(axis="both", labelsize=12)
        leg = ax.legend(fontsize="x-small")
        leg.set_draggable(True)

    # y-limit from TRANSP traces so noisy EQ data at start/end doesn't squash the plot
    ax_dia.set_ylim(top=1.5 * max(np.nanmax(r["wdia"]) for r in runs) * ENERGY_SCALE)
    ax_th.set_ylim(top=1.5 * max(np.nanmax(r["uthrm"]) for r in runs) * ENERGY_SCALE)

    cornernote(pulse)
    fig.tight_layout()
    return fig


def parse_arguments():
    p = argparse.ArgumentParser(
        description="Compare thermal and MHD stored energy of several TRANSP runs."
    )
    p.add_argument("pulse", type=int, help="JET pulse number")
    p.add_argument("runids", nargs="+", help="TRANSP run ids, e.g. M01 M02 M03")
    p.add_argument("--eq", default="EFIT", help="equilibrium DDA for WP/WDIA (default: EFIT)")
    p.add_argument("--eq-uid", default="jetppf", help="equilibrium PPF uid (default: jetppf)")
    p.add_argument("--eq-seq", type=int, default=0, help="equilibrium PPF sequence (default: 0)")
    p.add_argument("--toffset", type=float, default=40.0,
                   help="JET time - TRANSP time offset [s] (default: 40)")
    p.add_argument("--no-exp", action="store_true",
                   help="TRANSP only: skip PPF equilibrium data")
    return p.parse_args()


def main():
    args = parse_arguments()

    runs = [r for r in (load_run(args.pulse, runid) for runid in args.runids) if r is not None]
    if not runs:
        print("No usable TRANSP runs.")
        return 1

    exp = {}
    if not args.no_exp:
        if ppf is None:
            print("ppf module not available: plotting TRANSP data only.")
        else:
            for dty in ("WP", "WDIA"):
                d = read_eq(args.pulse, args.eq, dty, args.eq_uid, args.eq_seq)
                if d is not None:
                    exp[dty] = (d[0] - args.toffset, d[1])

    fig = plot_stored_energy(args.pulse, runs, exp, args.eq)
    win = pw.plotWindow()
    win.addPlot("stored energy", fig)
    win.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
