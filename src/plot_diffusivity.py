#!/usr/bin/env python
"""
Plot TRANSP diffusivity coefficients as profile snapshots and time traces.

The default signal list is the set of TRANSP output variables whose
descriptions contain "diffusivity", based on the user-supplied list.

A residual panel is included using TRANSP transport-equation residual profiles.
"""

import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np

import plotWindow as pw
import profiles as ps


ALL_DIFFUSIVITY_TERMS = [
    "BDIFBX_D",
    "CONDIWNC",
    "CONDWNCE",
    "CONDWNCX",
    "CONDNC_BE",
    "CONDNC_NE",
    "CONDNC_NI",
    "DFENC",
    "DFINC_BE",
    "PTDFI_BE",
    "DFINC_NE",
    "PTDFI_NE",
    "DFINC_NI",
    "PTDFI_NI",
    "DFINC_D",
    "DFINS_D",
    "CONDE",
    "DIFB",
    "DIFFE",
    "DIFFNE",
    "DIFWE",
    "DEINT",
    "CONDI",
    "DIFFIGLF",
    "ETPHIGLF",
    "DIFFINEO",
    "ETPHINEO",
    "XKHDRBM",
    "CHPHI",
    "DIFFX",
    "DIFFD",
    "CONDWNCD",
    "DIFFI",
    "DIFFI0",
    "DFI_D",
    "DIFBX",
    "KAPA",
    "KAPA6",
    "KAPAN",
]

DEFAULT_DIFFUSIVITY_TERMS = [
    "CONDE",
    "CONDI",
    "DIFFE",
    "DIFFI",
    "DIFWE",
]

RESIDUAL_TERMS = [
    "RESPROFPE",
    "RESPROFPI",
    "RESPROFTE",
    "RESPROFTI",
]

DIFFUSIVITY_SET = set(ALL_DIFFUSIVITY_TERMS)
COLORS = list(plt.get_cmap("tab20").colors)
TERM_STYLES = {
    signal: {"color": COLORS[idx % len(COLORS)], "linestyle": "-", "linewidth": 1.9}
    for idx, signal in enumerate(ALL_DIFFUSIVITY_TERMS)
}

RESIDUAL_STYLES = {
    "RESPROFPE": {"color": "tab:blue", "linestyle": "-", "linewidth": 1.9},
    "RESPROFPI": {"color": "tab:orange", "linestyle": "-", "linewidth": 1.9},
    "RESPROFTE": {"color": "tab:green", "linestyle": "-", "linewidth": 1.9},
    "RESPROFTI": {"color": "tab:red", "linestyle": "-", "linewidth": 1.9},
}


def gi(val, data):
    """Return index of data closest to val."""
    return np.argmin(np.abs(data - val))


def convert_units(units):
    conversion = {
        "CM**2/SEC": (r"$m^2/s$", 1e-4),
        "CM^2/SEC": (r"$m^2/s$", 1e-4),
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


def plot_profile(transp, signals, time_value, pulse, win):
    fig = plt.figure()
    fig.suptitle(f"{transp.transpcid} diffusivity at t = {time_value:.2f}", fontsize=13)
    ax = fig.add_subplot(111)

    unit, factor = get_common_units(transp, signals)
    time_index = gi(transp.t, time_value)

    for signal in signals:
        data = transp.transp(signal)
        if data is None or data.ndim != 2:
            print(f"Skipping non-profile signal {signal} in profile plot")
            continue
        style = TERM_STYLES.get(signal, {})
        ax.plot(
            transp.x[time_index, :],
            data[time_index, :] * factor,
            color=style.get("color", "tab:blue"),
            linestyle=style.get("linestyle", "-"),
            linewidth=style.get("linewidth", 1.9),
            label=signal,
        )

    ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    ax.tick_params(axis="both", labelsize=14)
    ax.grid(axis="y", linestyle=":", linewidth=0.8, color="0.8")
    ax.set_xlabel(r"$\rho_{tor}^{norm}$")
    ax.set_ylabel(unit)
    ax.set_xlim(0.0, 1.0)
    cornernote(ax, pulse)
    leg = ax.legend()
    leg.set_draggable(True)
    win.addPlot(f"diffusivity profile @ {time_value:.2f}s", fig)


def plot_time_trace(transp, signals, x_pos, pulse, win):
    fig = plt.figure()
    fig.suptitle(f"{transp.transpcid} diffusivity at x = {x_pos:.2f}", fontsize=13)
    ax = fig.add_subplot(111)

    unit, factor = get_common_units(transp, signals)
    x_index = gi(transp.x[0], x_pos)

    for signal in signals:
        data = transp.transp(signal)
        if data is None or data.ndim != 2:
            print(f"Skipping non-profile signal {signal} in time trace plot")
            continue
        style = TERM_STYLES.get(signal, {})
        ax.plot(
            transp.t,
            data[:, x_index] * factor,
            color=style.get("color", "tab:blue"),
            linestyle=style.get("linestyle", "-"),
            linewidth=style.get("linewidth", 1.9),
            label=signal,
        )

    ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    ax.tick_params(axis="both", labelsize=14)
    ax.grid(axis="y", linestyle=":", linewidth=0.8, color="0.8")
    ax.set_xlabel("time [s]")
    ax.set_ylabel(unit)
    cornernote(ax, pulse)
    leg = ax.legend()
    leg.set_draggable(True)
    win.addPlot(f"diffusivity time trace @ {x_pos:.2f}", fig)


def plot_residual_profile(transp, signals, time_value, pulse, win):
    fig = plt.figure()
    fig.suptitle(f"{transp.transpcid} residuals at t = {time_value:.2f}", fontsize=13)
    ax = fig.add_subplot(111)

    time_index = gi(transp.t, time_value)

    for signal in signals:
        data = transp.transp(signal)
        if data is None or data.ndim != 2:
            print(f"Skipping non-profile signal {signal} in residual plot")
            continue
        style = RESIDUAL_STYLES.get(signal, {})
        ax.plot(
            transp.x[time_index, :],
            data[time_index, :],
            color=style.get("color", "tab:blue"),
            linestyle=style.get("linestyle", "-"),
            linewidth=style.get("linewidth", 1.9),
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
    win.addPlot(f"residual profile @ {time_value:.2f}s", fig)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plot TRANSP diffusivity coefficients as profile and time-trace views."
    )
    parser.add_argument("pulse", type=int, help="JET discharge number")
    parser.add_argument("runid", type=str, help="TRANSP run id")
    parser.add_argument(
        "--time",
        type=float,
        nargs="+",
        default=[9.0],
        help="One or more TRANSP time slices for diffusivity profiles (default: 9.0)",
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
        default=DEFAULT_DIFFUSIVITY_TERMS,
        help="Subset of diffusivity signals to plot",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()

    transp, loaded = load_run(args.pulse, args.runid, args.signals)
    _, residuals = load_run(args.pulse, args.runid, RESIDUAL_TERMS)
    if not loaded and not residuals:
        print("No requested diffusivity signals were found in this run.")
        return 1

    win = pw.plotWindow()
    for time_value in args.time:
        plot_profile(transp, loaded, time_value, args.pulse, win)
        if residuals:
            plot_residual_profile(transp, residuals, time_value, args.pulse, win)
    plot_time_trace(transp, loaded, args.xpos, args.pulse, win)
    win.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
