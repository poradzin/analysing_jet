#!/usr/bin/env python
import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np

import plotWindow as pw
import profiles as ps


POWER_BALANCE_TERMS = [
    "PBTH",
    "GAINI",
    "PCOND",
    "QIE",
    "P0NET",
    "PCONV",
    "QROT",
    "IHEAT",
    "TIBAL",
]

TERM_STYLES = {
    "PBTH": {"color": "tab:blue", "linestyle": "-"},
    "GAINI": {"color": "tab:orange", "linestyle": "-"},
    "QROT": {"color": "tab:green", "linestyle": "-"},
    "PCOND": {"color": "tab:red", "linestyle": "--"},
    "QIE": {"color": "tab:purple", "linestyle": "--"},
    "P0NET": {"color": "tab:brown", "linestyle": "--"},
    "PCONV": {"color": "tab:pink", "linestyle": "--"},
    "IHEAT": {"color": "black", "linestyle": "-", "linewidth": 2.4},
    "TIBAL": {"color": "dimgray", "linestyle": "-.", "linewidth": 2.2},
}


def gi(val, data):
    """Return index of data closest to val."""
    return np.argmin(np.abs(data - val))


def convert_units(units):
    conversion = {
        "CM**2/SEC": (r"$m^2/s$", 1e-4),
        "WATTS/CM3": (r"$W/m^3$", 1e6),
        "W/CM3": (r"$W/m^3$", 1e6),
        "WATTS": (r"$W$", 1.0),
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
    fig.suptitle(f"{transp.transpcid} power balance at t = {time_value:.2f}", fontsize=13)
    ax = fig.add_subplot(111)

    unit, factor = get_common_units(transp, signals)
    time_index = gi(transp.t, time_value)

    for signal in signals:
        data = transp.transp(signal)
        if data is None:
            continue
        style = TERM_STYLES.get(signal, {})
        ax.plot(
            transp.x[time_index, :],
            data[time_index, :] * factor,
            color=style.get("color", "tab:blue"),
            linestyle=style.get("linestyle", "-"),
            linewidth=style.get("linewidth", 1.8),
            label=signal,
        )

    ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    ax.tick_params(axis="both", labelsize=14)
    ax.set_xlabel(r"$\rho_{tor}^{norm}$")
    ax.set_ylabel(unit)
    ax.set_xlim(0.0, 1.0)
    cornernote(ax, pulse)
    leg = ax.legend()
    leg.set_draggable(True)
    win.addPlot(f"power balance profile @ {time_value:.2f}s", fig)


def plot_time_trace(transp, signals, x_pos, pulse, win):
    fig = plt.figure()
    fig.suptitle(f"{transp.transpcid} power balance at x = {x_pos:.2f}", fontsize=13)
    ax = fig.add_subplot(111)

    unit, factor = get_common_units(transp, signals)
    x_index = gi(transp.x[0], x_pos)

    for signal in signals:
        data = transp.transp(signal)
        if data is None:
            continue
        style = TERM_STYLES.get(signal, {})
        ax.plot(
            transp.t,
            data[:, x_index] * factor,
            color=style.get("color", "tab:blue"),
            linestyle=style.get("linestyle", "-"),
            linewidth=style.get("linewidth", 1.8),
            label=signal,
        )

    ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    ax.tick_params(axis="both", labelsize=14)
    ax.set_xlabel("time [s]")
    ax.set_ylabel(unit)
    cornernote(ax, pulse)
    leg = ax.legend()
    leg.set_draggable(True)
    win.addPlot(f"power balance time trace @ {x_pos:.2f}", fig)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plot TRANSP ion power-balance terms as profiles and time traces."
    )
    parser.add_argument("pulse", type=int, help="JET discharge number")
    parser.add_argument("runid", type=str, help="TRANSP run id")
    parser.add_argument(
        "--time",
        type=float,
        nargs="+",
        default=[9.0],
        help="One or more TRANSP time slices for profile plots (default: 9.0)",
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
        default=POWER_BALANCE_TERMS,
        help="Subset of power-balance signals to plot",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()

    transp, loaded = load_run(args.pulse, args.runid, args.signals)
    if not loaded:
        print("No requested power-balance signals were found in this run.")
        return 1

    win = pw.plotWindow()
    for time_value in args.time:
        plot_profile(transp, loaded, time_value, args.pulse, win)
    plot_time_trace(transp, loaded, args.xpos, args.pulse, win)
    win.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
