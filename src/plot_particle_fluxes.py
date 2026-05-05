#!/usr/bin/env python
"""
Plot TRANSP particle-flux divergence terms as profile snapshots and time traces.

Default particle-balance check terms:
  - direct divergence terms for electrons, thermal ions, impurity ions, and D+ ions
  - residual profile terms from TRANSP particle-balance outputs

Relevant source:
  - TRANSP particle balance overview:
    https://transp.pppl.gov/modules/ptcl_balance.html
  - TRANSP output variable table:
    https://transp.pppl.gov/nml/transp_output.html
"""

import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np

import plotWindow as pw
import profiles as ps


ALL_PARTICLE_FLUX_TERMS = [
    "GFLNC_E",
    "GFLNC_X",
    "GFLNC_I",
    "DIVFE",
    "EPTR_MOD",
    "EPTR_OBS",
    "DIVFI",
    "IPTR_MOD",
    "IPTR_OBS",
    "DFIMP",
    "XPTR_MOD",
    "XPTR_OBS",
    "DIVFD",
    "PTRD_MOD",
    "PTRD_OBS",
]

DEFAULT_PARTICLE_FLUX_TERMS = [
    "DIVFE",
    "DIVFI",
    "DFIMP",
    "DIVFD",
]

RESIDUAL_PROFILE_TERMS = [
    "RESPROFPE",
    "RESPROFPI",
    "RESPROFPX",
]

RESIDUAL_L2_TERMS = [
    "RESL2PE",
    "RESL2PI",
    "RESL2PX",
]

COLORS = list(plt.get_cmap("tab20").colors)
TERM_STYLES = {
    signal: {"color": COLORS[idx % len(COLORS)], "linestyle": "-", "linewidth": 1.9}
    for idx, signal in enumerate(ALL_PARTICLE_FLUX_TERMS)
}

RESIDUAL_STYLES = {
    "RESPROFPE": {"color": "tab:blue"},
    "RESPROFPI": {"color": "tab:orange"},
    "RESPROFPX": {"color": "tab:green"},
    "RESL2PE": {"color": "tab:blue"},
    "RESL2PI": {"color": "tab:orange"},
    "RESL2PX": {"color": "tab:green"},
}


def gi(val, data):
    """Return index of data closest to val."""
    return np.argmin(np.abs(data - val))


def convert_units(units):
    conversion = {
        "N/CM3/SEC": (r"$m^{-3}s^{-1}$", 1e6),
        "N/CM^3/SEC": (r"$m^{-3}s^{-1}$", 1e6),
        "WATTS/CM3": (r"$W/m^3$", 1e6),
        "W/CM3": (r"$W/m^3$", 1e6),
    }
    key = units.strip().upper()
    if key in conversion:
        return conversion[key]
    return (units, 1.0)


def cornernote(axis, pulse):
    axis.annotate(
        str(pulse),
        xy=(0.98, 0.02),
        xycoords="figure fraction",
        fontsize=9,
        ha="right",
        va="bottom",
        label="cornernote",
    )


def load_run(pulse, runid, signals):
    transp = ps.Transp(pulse, runid)
    loaded = []
    for signal in signals:
        if signal in transp.all_keys:
            transp.add_data(signal)
            loaded.append(signal)
        else:
            print(f"Skipping missing signal {signal} for run {runid}")
    return transp, loaded


def get_common_units(transp, signals):
    units = [transp.units(signal) for signal in signals if transp.units(signal) is not None]
    if not units:
        return (r"$units$", 1.0)
    if all(unit == units[0] for unit in units):
        return convert_units(units[0])
    return (r"$units$", 1.0)


def signed_symlog_threshold(values):
    flat = np.concatenate([np.ravel(np.asarray(v)) for v in values if v is not None]) if values else np.array([])
    flat = flat[np.isfinite(flat)]
    flat = flat[flat != 0.0]
    if flat.size == 0:
        return 1.0
    max_abs = np.max(np.abs(flat))
    return max(max_abs * 1e-3, 1e-12)


def plot_profile(transp, signals, time_value, pulse, win):
    fig = plt.figure()
    fig.suptitle(f"{transp.transpcid} particle-flux divergence at t = {time_value:.2f}", fontsize=13)
    ax = fig.add_subplot(111)

    unit, factor = get_common_units(transp, signals)
    time_index = gi(transp.t, time_value)
    plotted = []

    for signal in signals:
        data = transp.transp(signal)
        if data is None or data.ndim != 2:
            print(f"Skipping non-profile signal {signal} in profile plot")
            continue
        style = TERM_STYLES.get(signal, {})
        yvals = data[time_index, :] * factor
        plotted.append(yvals)
        ax.plot(
            transp.x[time_index, :],
            yvals,
            color=style.get("color", "tab:blue"),
            linestyle=style.get("linestyle", "-"),
            linewidth=style.get("linewidth", 1.9),
            label=signal,
        )

    ax.set_yscale("symlog", linthresh=signed_symlog_threshold(plotted))
    ax.tick_params(axis="both", labelsize=14)
    ax.grid(axis="y", linestyle=":", linewidth=0.8, color="0.8")
    ax.set_xlabel(r"$\rho_{tor}^{norm}$")
    ax.set_ylabel(unit)
    ax.set_xlim(0.0, 1.0)
    cornernote(ax, pulse)
    leg = ax.legend()
    leg.set_draggable(True)
    win.addPlot(f"particle flux profile @ {time_value:.2f}s", fig)


def plot_time_trace(transp, signals, x_pos, pulse, win):
    fig = plt.figure()
    fig.suptitle(f"{transp.transpcid} particle-flux divergence at x = {x_pos:.2f}", fontsize=13)
    ax = fig.add_subplot(111)

    unit, factor = get_common_units(transp, signals)
    x_index = gi(transp.x[0], x_pos)
    plotted = []

    for signal in signals:
        data = transp.transp(signal)
        if data is None or data.ndim != 2:
            print(f"Skipping non-profile signal {signal} in time trace plot")
            continue
        style = TERM_STYLES.get(signal, {})
        yvals = data[:, x_index] * factor
        plotted.append(yvals)
        ax.plot(
            transp.t,
            yvals,
            color=style.get("color", "tab:blue"),
            linestyle=style.get("linestyle", "-"),
            linewidth=style.get("linewidth", 1.9),
            label=signal,
        )

    ax.set_yscale("symlog", linthresh=signed_symlog_threshold(plotted))
    ax.tick_params(axis="both", labelsize=14)
    ax.grid(axis="y", linestyle=":", linewidth=0.8, color="0.8")
    ax.set_xlabel("time [s]")
    ax.set_ylabel(unit)
    cornernote(ax, pulse)
    leg = ax.legend()
    leg.set_draggable(True)
    win.addPlot(f"particle flux time trace @ {x_pos:.2f}", fig)


def plot_residual_profile(transp, signals, time_value, pulse, win):
    fig = plt.figure()
    fig.suptitle(f"{transp.transpcid} particle-balance residuals at t = {time_value:.2f}", fontsize=13)
    ax = fig.add_subplot(111)

    time_index = gi(transp.t, time_value)

    for signal in signals:
        data = transp.transp(signal)
        if data is None or data.ndim != 2:
            print(f"Skipping non-profile signal {signal} in residual profile plot")
            continue
        style = RESIDUAL_STYLES.get(signal, {})
        ax.plot(
            transp.x[time_index, :],
            data[time_index, :],
            color=style.get("color", "tab:blue"),
            linestyle="-",
            linewidth=1.9,
            label=signal,
        )

    ax.axhline(0.0, color="0.3", linewidth=1.0)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    ax.tick_params(axis="both", labelsize=14)
    ax.grid(axis="y", linestyle=":", linewidth=0.8, color="0.8")
    ax.set_xlabel(r"$\rho_{tor}^{norm}$")
    ax.set_ylabel("Residual")
    ax.set_xlim(0.0, 1.0)
    cornernote(ax, pulse)
    leg = ax.legend()
    leg.set_draggable(True)
    win.addPlot(f"particle residual profile @ {time_value:.2f}s", fig)


def plot_residual_l2(transp, signals, pulse, win):
    fig = plt.figure()
    fig.suptitle(f"{transp.transpcid} particle-balance L2 residuals", fontsize=13)
    ax = fig.add_subplot(111)

    for signal in signals:
        data = transp.transp(signal)
        if data is None:
            continue
        series = np.squeeze(np.asarray(data))
        if series.ndim != 1:
            print(f"Skipping non-1D residual signal {signal} in L2 plot")
            continue
        style = RESIDUAL_STYLES.get(signal, {})
        ax.plot(
            transp.t,
            series,
            color=style.get("color", "tab:blue"),
            linestyle="-",
            linewidth=1.9,
            label=signal,
        )

    ax.set_yscale("log")
    ax.tick_params(axis="both", labelsize=14)
    ax.grid(axis="y", linestyle=":", linewidth=0.8, color="0.8")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("L2 residual")
    cornernote(ax, pulse)
    leg = ax.legend()
    leg.set_draggable(True)
    win.addPlot("particle residual L2", fig)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plot TRANSP particle-flux divergence terms as profile and time-trace views."
    )
    parser.add_argument("pulse", type=int, help="JET discharge number")
    parser.add_argument("runid", type=str, help="TRANSP run id")
    parser.add_argument(
        "--time",
        type=float,
        nargs="+",
        default=[9.0],
        help="One or more TRANSP time slices for particle-flux profiles (default: 9.0)",
    )
    parser.add_argument(
        "--xpos",
        type=float,
        default=0.5,
        help="Normalized radius for the time-trace plot (default: 0.5)",
    )
    parser.add_argument(
        "--signals",
        nargs="*",
        default=DEFAULT_PARTICLE_FLUX_TERMS,
        help="Subset of particle-flux signals to plot",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()

    transp, loaded = load_run(args.pulse, args.runid, args.signals)
    _, residual_profiles = load_run(args.pulse, args.runid, RESIDUAL_PROFILE_TERMS)
    _, residual_l2 = load_run(args.pulse, args.runid, RESIDUAL_L2_TERMS)

    if not loaded and not residual_profiles and not residual_l2:
        print("No requested particle-flux signals were found in this run.")
        return 1

    win = pw.plotWindow()
    for time_value in args.time:
        if loaded:
            plot_profile(transp, loaded, time_value, args.pulse, win)
        if residual_profiles:
            plot_residual_profile(transp, residual_profiles, time_value, args.pulse, win)
    if loaded:
        plot_time_trace(transp, loaded, args.xpos, args.pulse, win)
    if residual_l2:
        plot_residual_l2(transp, residual_l2, args.pulse, win)
    win.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
