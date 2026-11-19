#!/usr/bin/env python
"""
Plot TRANSP ion power-balance terms as volume-integrated balance checks.

Signed IEBAL convention used here:
  +PBTH -GAINI -PCOND +QIE -P0NET -PCONV +QROT +IHEAT +TIBAL

Relevant source:
  - TRANSP RPLOT IEBAL package definition with signs:
    https://w3.pppl.gov/~xshare/Rplot/nstx_multi.pdf
  - TRANSP output variable table:
    https://transp.pppl.gov/nml/transp_output.html
"""

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

BALANCE_SIGNS = {
    "PBTH": 1.0,
    "GAINI": -1.0,
    "PCOND": -1.0,
    "QIE": 1.0,
    "P0NET": -1.0,
    "PCONV": -1.0,
    "QROT": 1.0,
    "IHEAT": 1.0,
    "TIBAL": 1.0,
}

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

VOLUME_SCALE = 1e-6
VOLUME_UNIT = "MW"


def gi(val, data):
    """Return index of data closest to val."""
    return np.argmin(np.abs(data - val))


def window_indices(times, center, window):
    """Return indices of times within a window of total width `window` centred on `center`.

    Falls back to the single closest index when no window is given or no
    sample lies inside the window.
    """
    times = np.asarray(times)
    if window is None or window <= 0:
        return np.array([gi(times, center)])
    idx = np.where(np.abs(times - center) <= 0.5 * window)[0]
    if idx.size == 0:
        idx = np.array([gi(times, center)])
    return idx


def rolling_mean(times, values, window):
    """Centred rolling mean of `values` over a time window of total width `window`.

    Works on non-uniform time bases and ignores NaNs.
    """
    values = np.asarray(values, dtype=float)
    if window is None or window <= 0:
        return values
    times = np.asarray(times)
    valid = np.isfinite(values)
    csum = np.concatenate(([0.0], np.cumsum(np.where(valid, values, 0.0))))
    ccount = np.concatenate(([0], np.cumsum(valid)))
    lo = np.searchsorted(times, times - 0.5 * window, side="left")
    hi = np.searchsorted(times, times + 0.5 * window, side="right")
    counts = ccount[hi] - ccount[lo]
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(counts > 0, (csum[hi] - csum[lo]) / counts, np.nan)


def window_label(times, idx):
    """Describe the TRANSP time slices actually used for a snapshot."""
    used = np.asarray(times)[idx]
    if used.size == 1:
        return f"t = {used[0]:.3f} s"
    return f"<t> = {used.mean():.3f} s (mean over {used[0]:.3f}-{used[-1]:.3f} s, {used.size} slices)"


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


def integrate_volume(transp, signal):
    data = transp.transp(signal)
    if data is None:
        return None

    if data.ndim == 1:
        return data

    dvol = np.asarray(transp.dvol)
    if dvol.ndim == 1:
        dvol = np.broadcast_to(dvol, data.shape)
    elif dvol.shape != data.shape:
        dvol = np.broadcast_to(dvol, data.shape)

    return np.nansum(data * dvol, axis=1)


def zone_contribution(transp, signal):
    data = transp.transp(signal)
    if data is None:
        return None

    if data.ndim == 1:
        return data

    dvol = np.asarray(transp.dvol)
    if dvol.ndim == 1:
        dvol = np.broadcast_to(dvol, data.shape)
    elif dvol.shape != data.shape:
        dvol = np.broadcast_to(dvol, data.shape)

    return data * dvol


def get_volume_integrated_series(transp, signals):
    integrated = {}
    for signal in signals:
        series = integrate_volume(transp, signal)
        if series is not None:
            integrated[signal] = series
    return integrated


def plot_snapshot(transp, integrated, signals, time_value, pulse, win, window=None):
    time_idx = window_indices(transp.t, time_value, window)
    fig = plt.figure()
    prefix = "time-averaged " if window else ""
    fig.suptitle(
        f"{transp.transpcid} {prefix}volume-integrated power balance at {window_label(transp.t, time_idx)}",
        fontsize=13,
    )
    ax = fig.add_subplot(111)

    values = []
    labels = []
    colors = []
    for signal in signals:
        if signal not in integrated:
            continue
        signed_value = BALANCE_SIGNS.get(signal, 1.0) * np.nanmean(integrated[signal][time_idx]) * VOLUME_SCALE
        values.append(signed_value)
        labels.append(signal)
        colors.append(TERM_STYLES.get(signal, {}).get("color", "tab:blue"))

    ax.barh(labels, values, color=colors, edgecolor="black", linewidth=0.5)
    ax.axvline(0.0, color="0.3", linewidth=1.0)
    ax.tick_params(axis="both", labelsize=13)
    ax.set_xlabel(f"Volume-integrated power [{VOLUME_UNIT}]")
    closure = float(np.nansum(values))
    ax.set_title(f"Signed closure = {closure:.3f} {VOLUME_UNIT}", fontsize=11)

    if values:
        span = max(abs(min(values)), abs(max(values)))
        if span > 0:
            ax.set_xlim(-1.15 * span, 1.15 * span)

    cornernote(ax, pulse)
    win.addPlot(f"power balance @ {time_value:.2f}s", fig)


def plot_zone_snapshot(transp, signals, time_value, pulse, win, window=None):
    time_idx = window_indices(transp.t, time_value, window)
    fig = plt.figure()
    prefix = "time-averaged " if window else ""
    fig.suptitle(
        f"{transp.transpcid} {prefix}zone contributions at {window_label(transp.t, time_idx)}",
        fontsize=13,
    )
    ax = fig.add_subplot(111)

    x_mean = np.nanmean(np.asarray(transp.x)[time_idx, :], axis=0)
    unit = "W"

    for signal in signals:
        contrib = zone_contribution(transp, signal)
        if contrib is None:
            continue
        style = TERM_STYLES.get(signal, {})
        if contrib.ndim == 1:
            continue
        signed_contrib = BALANCE_SIGNS.get(signal, 1.0) * np.nanmean(contrib[time_idx, :], axis=0)
        ax.plot(
            x_mean,
            signed_contrib,
            color=style.get("color", "tab:blue"),
            linestyle=style.get("linestyle", "-"),
            linewidth=style.get("linewidth", 1.8),
            label=signal,
        )

    ax.axhline(0.0, color="0.3", linewidth=1.0)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    ax.tick_params(axis="both", labelsize=14)
    ax.set_xlabel(r"$\rho_{tor}^{norm}$")
    ax.set_ylabel(f"Zone contribution [{unit}]")
    ax.set_xlim(0.0, 1.0)
    cornernote(ax, pulse)
    leg = ax.legend()
    leg.set_draggable(True)
    win.addPlot(f"zone contributions @ {time_value:.2f}s", fig)


def plot_time_trace(transp, signals, pulse, win, window=None):
    fig = plt.figure()
    title = f"{transp.transpcid} volume-integrated power balance"
    if window:
        title += f" (rolling mean, {window:.3g} s window)"
    fig.suptitle(title, fontsize=13)
    ax = fig.add_subplot(111)

    integrated = get_volume_integrated_series(transp, signals)

    for signal in signals:
        if signal not in integrated:
            continue
        style = TERM_STYLES.get(signal, {})
        signed_series = rolling_mean(
            transp.t, BALANCE_SIGNS.get(signal, 1.0) * integrated[signal] * VOLUME_SCALE, window
        )
        ax.plot(
            transp.t,
            signed_series,
            color=style.get("color", "tab:blue"),
            linestyle=style.get("linestyle", "-"),
            linewidth=style.get("linewidth", 1.8),
            label=f"{'+' if BALANCE_SIGNS.get(signal, 1.0) > 0 else '-'}{signal}",
        )

    closure = np.zeros_like(transp.t, dtype=float)
    for signal in signals:
        if signal in integrated:
            closure += BALANCE_SIGNS.get(signal, 1.0) * integrated[signal] * VOLUME_SCALE
    closure = rolling_mean(transp.t, closure, window)
    ax.plot(
        transp.t,
        closure,
        color="black",
        linestyle="--",
        linewidth=2.2,
        label="signed sum",
    )

    ax.axhline(0.0, color="0.3", linewidth=1.0)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 3))
    ax.tick_params(axis="both", labelsize=14)
    ax.set_xlabel("time [s]")
    ax.set_ylabel(f"Power [{VOLUME_UNIT}]")
    cornernote(ax, pulse)
    leg = ax.legend()
    leg.set_draggable(True)
    win.addPlot("power balance time trace", fig)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Plot TRANSP ion power-balance terms as volume-integrated balance checks."
    )
    parser.add_argument("pulse", type=int, help="JET discharge number")
    parser.add_argument("runid", type=str, help="TRANSP run id")
    parser.add_argument(
        "--time",
        type=float,
        nargs="+",
        default=[9.0],
        help="One or more TRANSP time slices for volume-integrated snapshots (default: 9.0)",
    )
    parser.add_argument(
        "--window",
        type=float,
        default=None,
        help=(
            "Averaging time window width in seconds. Snapshots are averaged over "
            "[t - window/2, t + window/2]; the time trace is smoothed with a centred "
            "rolling mean of the same width (default: no averaging)"
        ),
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

    t_min, t_max = float(np.min(transp.t)), float(np.max(transp.t))
    print(f"TRANSP time range for {transp.transpcid}: {t_min:.3f} - {t_max:.3f} s")
    times = [tv for tv in args.time if t_min <= tv <= t_max]
    for tv in args.time:
        if tv not in times:
            print(f"WARNING: requested time {tv:.3f} s is outside the TRANSP time range; skipping.")
    if not times:
        print("None of the requested --time values lie within the run; only the time trace is plotted.")

    integrated = get_volume_integrated_series(transp, loaded)
    win = pw.plotWindow()
    for time_value in times:
        plot_snapshot(transp, integrated, loaded, time_value, args.pulse, win, args.window)
        plot_zone_snapshot(transp, loaded, time_value, args.pulse, win, args.window)
    plot_time_trace(transp, loaded, args.pulse, win, args.window)
    win.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
