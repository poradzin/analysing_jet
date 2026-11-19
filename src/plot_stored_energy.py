#!/usr/bin/env python
"""
Compare TRANSP stored-energy time traces for several runs of one JET pulse.

Two panels:
  top    - plasma energy:
             EQ/WDIA                                             (black)
             EQ/WP                                               (navy)
             WDIA (TRANSP) = 1.5*UFASTPP + UTHRM                 (red)
             WMHD (TRANSP) = 0.75*UFASTPP + 1.5*(UFASTPA + UPHI) + UTHRM  (dark orange)
  bottom - plasma thermal energy:
             EQ/WDIA - 1.5*UFASTPP                               (black)
             EQ/WP   - 0.75*UFASTPP - 1.5*(UFASTPA + UPHI)       (navy)
             UTHRM (TRANSP)                                      (red)

Colour = quantity, line style = TRANSP run (solid, dashed, dotted, ...).

The equilibrium (DDA/SEQ/UID) of each run is read from the MMX block of
<pulse><runid>_Provenance.DAT in the run directory. Runs sharing the same
equilibrium get its WDIA/WP plotted only once (in the style of the first such
run). --eq/--eq-uid/--eq-seq override the provenance for all runs.

Signal definitions follow OMFIT modules/TRANSP/PLOTS/JET/plot_summar_JET.py.
UPHI is used only when present in the run.

Usage:
  python plot_stored_energy.py 99173 M01 M02 M03
  python plot_stored_energy.py 99173 M01 M02 --eq EFTP --eq-uid myuid --eq-seq 12
"""

import argparse
import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

import plotWindow as pw
import profiles as ps

try:
    import ppf
except ImportError:
    ppf = None


ENERGY_SIGNALS = ["UTHRM", "UFASTPP", "UFASTPA", "UPHI"]
ENERGY_SCALE = 1e-6  # J -> MJ
ENERGY_UNIT = "MJ"
DEFAULT_EQ = ("EFIT", 0, "jetppf")

COLORS = {
    "eq_wdia": "black",
    "eq_wp": "navy",
    "tr_wdia": "red",
    "tr_wmhd": "darkorange",
    "th_wdia": "black",
    "th_wp": "navy",
    "th_uthrm": "red",
}

RUN_STYLES = ["-", "--", ":", "-.", (0, (5, 1, 1, 1, 1, 1)), (0, (8, 3)), (0, (1, 4))]


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


def read_provenance_eq(transp, pulse, runid):
    """(DDA, SEQ, UID) from the MMX block of <pulse><runid>_Provenance.DAT, or None."""
    rundir = os.path.dirname(transp._cdfloc)
    target = f"{pulse}{runid}_provenance.dat".lower()
    matches = [f for f in glob.glob(os.path.join(rundir, "*")) if os.path.basename(f).lower() == target]
    if not matches:
        print(f"{runid}: no Provenance file in {rundir}.")
        return None

    with open(matches[0]) as f:
        lines = f.read().splitlines()

    for i, line in enumerate(lines):
        if line.strip().upper() != "MMX:":
            continue
        info = {}
        for entry in lines[i + 1:i + 6]:
            if ":" not in entry:
                continue
            key, val = (s.strip() for s in entry.split(":", 1))
            info[key.upper()] = val
        try:
            return info["DDA"].upper(), int(info["SEQ"]), info["UID"].lower()
        except (KeyError, ValueError):
            break

    print(f"{runid}: could not parse MMX block in {matches[0]}.")
    return None


def load_run(pulse, runid, eq_override):
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

    if eq_override is not None:
        eq = eq_override
    else:
        eq = read_provenance_eq(transp, pulse, runid)
        if eq is None:
            eq = DEFAULT_EQ
            print(f"{runid}: falling back to {eq_label(eq)}.")
    print(f"{runid}: equilibrium {eq_label(eq)}")

    w = {sig: volume_integral(transp, sig) for sig in ENERGY_SIGNALS if sig in transp.all_keys}
    uphi = w.get("UPHI", np.zeros_like(w["UTHRM"]))
    fast_mhd = 0.75 * w["UFASTPP"] + 1.5 * (w["UFASTPA"] + uphi)

    return {
        "runid": runid,
        "label": transp.transpcid,
        "eq": eq,
        "t": np.asarray(transp.t),
        "uthrm": w["UTHRM"],
        "fast_mhd": fast_mhd,
        "fast_dia": 1.5 * w["UFASTPP"],
        "wmhd": fast_mhd + w["UTHRM"],
        "wdia": 1.5 * w["UFASTPP"] + w["UTHRM"],
        "has_uphi": "UPHI" in w,
    }


def eq_label(eq):
    dda, seq, uid = eq
    return f"{dda}/{uid}/{seq}"


def read_eq(pulse, eq, dty, toffset):
    """Equilibrium energy time trace -> (t_transp, data) or None."""
    if ppf is None:
        return None
    dda, seq, uid = eq
    out = ppf.ppfdata(pulse, dda, dty, seq=seq, uid=uid, device="JET")
    data, t, ier = np.asarray(out[0], dtype=float), np.asarray(out[2], dtype=float), out[-1]
    if ier != 0 or data.size == 0:
        print(f"No {dda}/{dty} (uid={uid}, seq={seq}) for pulse {pulse}.")
        return None
    return t - toffset, data


def on_grid(t_new, exp):
    """Interpolate an experimental trace onto TRANSP times (NaN outside its range)."""
    t, y = exp
    return np.interp(t_new, t, y, left=np.nan, right=np.nan)


def plot_stored_energy(pulse, runs, eq_data):
    fig, (ax_dia, ax_th) = plt.subplots(2, 1, sharex=True, figsize=(9, 8))
    fig.suptitle(f"JET {pulse}: stored energy", fontsize=13)
    lw = 2.0

    # Panel 1: equilibrium WDIA / WP, once per distinct equilibrium
    plotted_eqs = set()
    for i, run in enumerate(runs):
        ls = RUN_STYLES[i % len(RUN_STYLES)]
        eq = run["eq"]
        if eq in plotted_eqs:
            continue
        plotted_eqs.add(eq)
        wdia, wp = eq_data.get((eq, "WDIA")), eq_data.get((eq, "WP"))
        if wdia is not None:
            ax_dia.plot(wdia[0], wdia[1] * ENERGY_SCALE, color=COLORS["eq_wdia"], ls=ls, lw=lw)
        if wp is not None:
            ax_dia.plot(wp[0], wp[1] * ENERGY_SCALE, color=COLORS["eq_wp"], ls=ls, lw=lw)

    for i, run in enumerate(runs):
        ls = RUN_STYLES[i % len(RUN_STYLES)]
        t, eq = run["t"], run["eq"]
        wdia, wp = eq_data.get((eq, "WDIA")), eq_data.get((eq, "WP"))

        # Panel 1: TRANSP WDIA / WMHD
        ax_dia.plot(t, run["wdia"] * ENERGY_SCALE, color=COLORS["tr_wdia"], ls=ls, lw=lw)
        ax_dia.plot(t, run["wmhd"] * ENERGY_SCALE, color=COLORS["tr_wmhd"], ls=ls, lw=lw)

        # Panel 2: thermal energy
        if wdia is not None:
            ax_th.plot(t, (on_grid(t, wdia) - run["fast_dia"]) * ENERGY_SCALE,
                       color=COLORS["th_wdia"], ls=ls, lw=lw)
        if wp is not None:
            ax_th.plot(t, (on_grid(t, wp) - run["fast_mhd"]) * ENERGY_SCALE,
                       color=COLORS["th_wp"], ls=ls, lw=lw)
        ax_th.plot(t, run["uthrm"] * ENERGY_SCALE, color=COLORS["th_uthrm"], ls=ls, lw=lw)

    # Legends: colour = quantity (per panel), line style = run (top panel)
    ddas = {r["eq"][0] for r in runs}
    eq_name = ddas.pop() if len(ddas) == 1 else "EQ"  # per-run source is in the run legend
    any_uphi = any(r["has_uphi"] for r in runs)
    uphi = " + UPHI" if any_uphi else ""

    def colour_handles(items):
        return [Line2D([], [], color=c, lw=lw, label=lab) for c, lab in items]

    leg_q1 = ax_dia.legend(handles=colour_handles([
        (COLORS["eq_wdia"], f"{eq_name}/WDIA"),
        (COLORS["eq_wp"], f"{eq_name}/WP"),
        (COLORS["tr_wdia"], "WDIA (TRANSP) = 1.5 UFASTPP + UTHRM"),
        (COLORS["tr_wmhd"], f"WMHD (TRANSP) = 0.75 UFASTPP + 1.5 (UFASTPA{uphi}) + UTHRM"),
    ]), fontsize="x-small", loc="upper left")
    leg_q1.set_draggable(True)
    ax_dia.add_artist(leg_q1)

    run_handles = [
        Line2D([], [], color="0.3", lw=lw, ls=RUN_STYLES[i % len(RUN_STYLES)],
               label=f"{r['label']}  [{eq_label(r['eq'])}]")
        for i, r in enumerate(runs)
    ]
    leg_runs = ax_dia.legend(handles=run_handles, fontsize="x-small", loc="upper right",
                             title="TRANSP run [equilibrium]", title_fontsize="x-small")
    leg_runs.set_draggable(True)

    leg_q2 = ax_th.legend(handles=colour_handles([
        (COLORS["th_wdia"], f"WDIA ({eq_name}) - 1.5 UFASTPP (TRANSP)"),
        (COLORS["th_wp"], f"WP ({eq_name}) - 0.75 UFASTPP - 1.5 (UFASTPA{uphi}) (TRANSP)"),
        (COLORS["th_uthrm"], "UTHRM (TRANSP)"),
    ]), fontsize="x-small", loc="upper left")
    leg_q2.set_draggable(True)

    t_min = min(r["t"][0] for r in runs)
    t_max = max(r["t"][-1] for r in runs)
    ax_th.set_xlim(t_min, t_max)

    ax_dia.set_ylabel(f"Energy [{ENERGY_UNIT}]")
    ax_dia.set_title("Plasma energy (diamagnetic / MHD)", fontsize=11)
    ax_th.set_ylabel(f"Energy [{ENERGY_UNIT}]")
    ax_th.set_title("Plasma thermal energy", fontsize=11)
    ax_th.set_xlabel("time [s] (TRANSP)")
    for ax in (ax_dia, ax_th):
        ax.grid(alpha=0.3)
        ax.tick_params(axis="both", labelsize=12)

    # y-limit from TRANSP traces so noisy EQ data at start/end doesn't squash the plot
    ax_dia.set_ylim(0.0, 1.5 * max(np.nanmax(r["wdia"]) for r in runs) * ENERGY_SCALE)
    ax_th.set_ylim(0.0, 1.5 * max(np.nanmax(r["uthrm"]) for r in runs) * ENERGY_SCALE)

    cornernote(pulse)
    fig.tight_layout()
    return fig


def parse_arguments():
    p = argparse.ArgumentParser(
        description="Compare plasma and thermal stored energy of several TRANSP runs."
    )
    p.add_argument("pulse", type=int, help="JET pulse number")
    p.add_argument("runids", nargs="+", help="TRANSP run ids, e.g. M01 M02 M03")
    p.add_argument("--eq", default=None,
                   help="equilibrium DDA for all runs (default: from each run's Provenance MMX)")
    p.add_argument("--eq-uid", default="jetppf", help="equilibrium uid with --eq (default: jetppf)")
    p.add_argument("--eq-seq", type=int, default=0, help="equilibrium seq with --eq (default: 0)")
    p.add_argument("--toffset", type=float, default=40.0,
                   help="JET time - TRANSP time offset [s] (default: 40)")
    p.add_argument("--no-exp", action="store_true",
                   help="TRANSP only: skip PPF equilibrium data")
    return p.parse_args()


def main():
    args = parse_arguments()

    eq_override = (args.eq.upper(), args.eq_seq, args.eq_uid.lower()) if args.eq else None
    runs = [r for r in (load_run(args.pulse, runid, eq_override) for runid in args.runids)
            if r is not None]
    if not runs:
        print("No usable TRANSP runs.")
        return 1

    # Read each distinct equilibrium once
    eq_data = {}
    if not args.no_exp:
        if ppf is None:
            print("ppf module not available: plotting TRANSP data only.")
        else:
            for eq in dict.fromkeys(r["eq"] for r in runs):
                for dty in ("WDIA", "WP"):
                    d = read_eq(args.pulse, eq, dty, args.toffset)
                    if d is not None:
                        eq_data[(eq, dty)] = d

    fig = plot_stored_energy(args.pulse, runs, eq_data)
    win = pw.plotWindow()
    win.addPlot("stored energy", fig)
    win.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
