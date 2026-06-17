#!/usr/bin/env python3
"""extract_km14_png.py -- pull the curves and data points out of the published
KM14 diamond-spectrum figure (Nocente analysis) and save each as a 2-column
text file (E_dep_MeV, counts) for downstream forward-model fitting.

The published PNG has the same axis layout for every pulse: a 7.2-9.2 MeV E_dep
x-axis and a 0-200 counts y-axis. The pixel-to-data calibration is the same as
in km14_spectrum.py (PNG_CAL). Each fit component is drawn in a distinct colour
(Th orange-dashed, B-th blue-dotted, Scatt dark-green dashdot, Total red), and
the measurement is black markers + error bars; we pick each one out by HSV-hue
filtering.

Output (default into the same directory as the input PNG, one file per curve):
    <pulse>_KM14_data.txt      # data points (E_dep, counts)
    <pulse>_KM14_data_err.txt  # data points + error bars (E_dep, counts, err)
    <pulse>_KM14_th.txt        # orange dashed
    <pulse>_KM14_bt.txt        # blue dotted
    <pulse>_KM14_scatt.txt     # green dash-dot
    <pulse>_KM14_total.txt     # red solid

Data points are found from the error-bar geometry (peak per column = stem),
so points whose marker disk is hidden under the red Total curve are recovered,
and the error bars are read off the caps (err = half the bar length).

Each text file: two columns, header line. Sorted by E_dep ascending; dashed/
dash-dot curves have gaps where no matching pixel was found in a column --
let the downstream code np.interp across them.

Usage
-----
    python extract_km14_png.py ~/jet/data/104614/figs/104614_KM14_spectral_analysis.png
    python extract_km14_png.py path/to/img.png --pulse 104614 --outdir /tmp
    python extract_km14_png.py img.png --legend-bbox 340,40,650,150
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
from matplotlib.colors import rgb_to_hsv
import matplotlib.image as mpimg


# Pixel <-> data calibration of the published PNG axes (same as km14_spectrum.py)
PNG_CAL = dict(x0=129, e0=7.2, x1=725, e1=9.2,
               y0=584, c0=0.0, y1=32, c1=200.0)

# Bounding box (xpix0, ypix0, xpix1, ypix1) of the legend block to suppress.
# The Nocente legend has 6 rows (JPN id, t window, Th, B-th, Scatt, Total); the
# "Total" red swatch extends down to ypix~175, just above the Total-curve peak
# (ypix~184 at the 8.4 MeV peak). Override with --legend-bbox if a different
# figure has the legend elsewhere.
DEFAULT_LEGEND_BBOX = (560, 30, 730, 240)

# Data-point detection. Each error bar leaves a tall black stem at the point's
# exact x, so peaks in the per-column black-pixel count locate every point --
# including those whose marker disk is hidden under a coloured fit curve (the
# red Total especially). Tuned for the published ~1300x705 px figures.
DATA_PEAK_MIN_HEIGHT = 6      # min black pixels in a column to be a candidate
DATA_PEAK_MIN_DISTANCE = 5    # min px between adjacent data points
DATA_PEAK_MIN_PROMINENCE = 3  # min peak prominence in the column-count profile
DATA_BAR_MIN_RUN = 5          # min black run (px) kept when reading an error bar

# The central value is read from the marker disk wherever it is visible. The disk
# is the only black feature both wider than the stem (~3 px) and taller than the
# caps (~3 rows), so a square binary erosion isolates it; its centroid is the
# plotted value. Where a fit curve (esp. the red Total) hides the disk, no blob
# survives within DISK_MATCH_PX of the point and we fall back to the cap midpoint.
DATA_MARKER_ERODE = 4   # square erosion size (px) isolating marker disks
DISK_MARKER_MERGE = 3   # merge disk blobs within this many px in x (split disks)
DISK_MATCH_PX = 4       # max x-distance (px) to bind a disk to a detected point

# Per-curve hue / saturation / value filter parameters.
# Hue is on [0, 1] as returned by matplotlib's rgb_to_hsv. Hue wraps around if
# hue_lo > hue_hi (e.g. red, 0.93 -> 0.07).
CURVES = {
    "th":    dict(hue=(0.05, 0.13), sat_min=0.40, val_min=0.40,  # orange
                  color_name="orange",   marker_like=False),
    "bt":    dict(hue=(0.55, 0.72), sat_min=0.40, val_min=0.30,  # blue
                  color_name="blue",     marker_like=False),
    "scatt": dict(hue=(0.30, 0.40), sat_min=0.35, val_min=0.20,  # green
                  color_name="green",    marker_like=False),
    "total": dict(hue=(0.93, 0.07), sat_min=0.55, val_min=0.55,  # red
                  color_name="red",      marker_like=False),
    "data":  dict(hue=(0.00, 1.00), sat_min=0.00, val_min=0.00,  # black markers
                  color_name="black",    marker_like=True),
}


def _load_image(path):
    img = mpimg.imread(str(path))
    if img.dtype.kind in "ui" and img.max() > 1.0:
        img = img.astype(float) / 255.0
    return img


def _hue_mask(rgb, hue_lo, hue_hi, sat_min, val_min):
    hsv = rgb_to_hsv(rgb)
    H, S, V = hsv[..., 0], hsv[..., 1], hsv[..., 2]
    if hue_lo <= hue_hi:
        h_mask = (H >= hue_lo) & (H <= hue_hi)
    else:
        h_mask = (H >= hue_lo) | (H <= hue_hi)
    return h_mask & (S >= sat_min) & (V >= val_min)


def _black_mask(rgb, bright_max=0.30, sat_max=0.15):
    """Low brightness AND low saturation -> data markers (and grid/tick text,
    suppressed via the plot-area window)."""
    bright = rgb.mean(axis=-1)
    sat = rgb.max(axis=-1) - rgb.min(axis=-1)
    return (bright < bright_max) & (sat < sat_max)


def _strip_axis_frame(win, ymin, ymax, xmin, xmax, frac=0.5):
    """Blank the axis spines: any row/column that is black across more than
    `frac` of the plot span is a frame line, not data. Done in place so the
    column-count peak finder and the cap scan don't latch onto the frame."""
    sub = win[ymin:ymax + 1, xmin:xmax + 1]
    win[np.where(sub.mean(axis=1) > frac)[0] + ymin, :] = False
    win[:, np.where(sub.mean(axis=0) > frac)[0] + xmin] = False
    return win


def _error_bar_extent(win, xc, min_run=DATA_BAR_MIN_RUN):
    """Top/bottom pixel rows of the error bar at the stem column `xc`, or None.

    Read from the single stem column, not a +-1 px band: near the peak the point
    spacing (~5-7 px) is smaller than the cap width (~13 px), so a band catches
    neighbouring points' caps and inflates the bar. Keep only black runs
    >= `min_run` px -- the point's own caps connect to its stem (one long run,
    or two long runs split by the red curve), whereas a neighbour cap that bleeds
    into this column is a short isolated run (no stem here) and is dropped."""
    ys = np.where(win[:, xc])[0]
    if ys.size == 0:
        return None
    gaps = np.where(np.diff(ys) > 1)[0]
    starts = np.concatenate([[0], gaps + 1])
    ends = np.concatenate([gaps + 1, [ys.size]])
    runs = [ys[a:b] for a, b in zip(starts, ends)]
    kept = [r for r in runs if r.size >= min_run] or runs
    return int(min(r.min() for r in kept)), int(max(r.max() for r in kept))


def _disk_centroids(win, erode=DATA_MARKER_ERODE, merge_px=DISK_MARKER_MERGE):
    """Centroid (row, col) of each visible marker disk in the windowed black mask.

    A square binary erosion removes the error-bar stems (too narrow) and caps
    (too short), leaving the disks; labelling gives one blob per visible marker.
    Disks split in two by a fit curve crossing their middle (two same-x arcs) are
    merged back by averaging blobs within `merge_px` in x, so the centroid lands
    on the true disk centre."""
    from scipy import ndimage  # available in repo env; local import keeps top light
    eroded = ndimage.binary_erosion(win, structure=np.ones((erode, erode)))
    lab, n = ndimage.label(eroded)
    if n == 0:
        return np.empty((0, 2))
    cen = np.array(ndimage.center_of_mass(eroded, lab, range(1, n + 1)))
    cen = cen[np.argsort(cen[:, 1])]
    groups = [[cen[0]]]
    for r in cen[1:]:
        if r[1] - groups[-1][-1][1] <= merge_px:
            groups[-1].append(r)
        else:
            groups.append([r])
    return np.array([np.mean(g, axis=0) for g in groups])


def _extract_data_points(rgb, cal, legend_bbox=None):
    """Return (E_dep_MeV, counts, err, from_disk) for every data point.

    Points are *located* by the peak in the per-column black-pixel count -- the
    error-bar stem is black at the point's exact x even when the marker disk is
    hidden under a fit curve -- so points buried under the red Total curve are
    recovered. The *central value* is read from the marker disk wherever it is
    visible (`from_disk=True`); where no disk survives the curve overlap we fall
    back to the cap midpoint (`from_disk=False`). The error bar is always read
    from the cap extent in a +-1 px band: 1-sigma = half the cap separation."""
    from scipy import signal  # available in repo env; local import keeps top light
    mask = _black_mask(rgb)
    mask = _apply_legend_bbox(mask, legend_bbox)
    ymin, ymax = sorted((int(cal["y0"]), int(cal["y1"])))
    xmin, xmax = sorted((int(cal["x0"]), int(cal["x1"])))
    win = np.zeros_like(mask)
    win[ymin:ymax + 1, xmin:xmax + 1] = mask[ymin:ymax + 1, xmin:xmax + 1]
    win = _strip_axis_frame(win, ymin, ymax, xmin, xmax)

    disks = _disk_centroids(win)            # (row, col) per visible marker
    disk_x = disks[:, 1] if disks.size else np.empty(0)

    colcount = win[ymin:ymax + 1, xmin:xmax + 1].sum(axis=0)
    peaks, _ = signal.find_peaks(colcount, height=DATA_PEAK_MIN_HEIGHT,
                                 distance=DATA_PEAK_MIN_DISTANCE,
                                 prominence=DATA_PEAK_MIN_PROMINENCE)
    e_list, c_list, err_list, src_list = [], [], [], []
    for xc in xmin + peaks:
        ext = _error_bar_extent(win, xc)
        if ext is None:
            continue
        y_top, y_bot = ext
        e_val, c_top = _pixel_to_data(float(xc), float(y_top), cal)
        _, c_bot = _pixel_to_data(float(xc), float(y_bot), cal)
        # Central value: marker disk if one is bound to this point, else cap mid.
        on_disk = disk_x.size and np.min(np.abs(disk_x - xc)) <= DISK_MATCH_PX
        if on_disk:
            row = disks[int(np.argmin(np.abs(disk_x - xc))), 0]
            _, c_val = _pixel_to_data(float(xc), float(row), cal)
        else:
            c_val = 0.5 * (c_top + c_bot)
        e_list.append(e_val)
        c_list.append(c_val)
        err_list.append(0.5 * (c_top - c_bot))
        src_list.append(bool(on_disk))
    e = np.array(e_list)
    c = np.clip(np.array(c_list), 0.0, None)
    err = np.array(err_list)
    from_disk = np.array(src_list, dtype=bool)
    order = np.argsort(e)
    return e[order], c[order], err[order], from_disk[order]


def _column_centroid(mask, ymin, ymax, xpix):
    """For each x-column, return the median row of matching pixels (or NaN)."""
    out = np.full(xpix.size, np.nan, dtype=float)
    for i, xp in enumerate(xpix):
        col = mask[ymin:ymax + 1, xp]
        idx = np.where(col)[0]
        if idx.size == 0:
            continue
        out[i] = ymin + float(np.median(idx))
    return out


def _apply_legend_bbox(mask, bbox):
    if bbox is None:
        return mask
    x0, y0, x1, y1 = bbox
    out = mask.copy()
    out[y0:y1 + 1, x0:x1 + 1] = False
    return out


def _pixel_to_data(xpix, ypix, cal):
    e = cal["e0"] + (xpix - cal["x0"]) * (cal["e1"] - cal["e0"]) / (cal["x1"] - cal["x0"])
    cv = cal["c0"] + (ypix - cal["y0"]) * (cal["c1"] - cal["c0"]) / (cal["y1"] - cal["y0"])
    return e, cv


def _hampel_filter(e, c, window=11, n_sigma=4.0):
    """Reject outliers from (e, c) by comparing each c to the rolling median over
    a sliding window of `window` neighbours (in index order, after sorting by E).
    Drops points where |c - median| > n_sigma * 1.4826 * MAD. Keeps both arrays
    in sync. No-op for n < window."""
    if e.size < window:
        return e, c
    half = window // 2
    keep = np.ones(e.size, dtype=bool)
    for i in range(e.size):
        a = max(0, i - half); b = min(e.size, i + half + 1)
        seg = c[a:b]
        med = np.median(seg)
        mad = np.median(np.abs(seg - med))
        if mad <= 0:
            continue
        if abs(c[i] - med) > n_sigma * 1.4826 * mad:
            keep[i] = False
    return e[keep], c[keep]


def extract_curve(img, cal, params, legend_bbox=None):
    """Run the extraction for one colour curve. Returns (E_dep_MeV, counts).

    Data markers are not handled here -- see _extract_data_points (the markers
    need error-bar logic, not a per-column hue centroid)."""
    rgb = img[..., :3]
    # Smooth colour curves: per-column centroid of the hue-matched pixels.
    hue_lo, hue_hi = params["hue"]
    mask = _hue_mask(rgb, hue_lo, hue_hi, params["sat_min"], params["val_min"])
    mask = _apply_legend_bbox(mask, legend_bbox)

    ymin, ymax = sorted((int(cal["y0"]), int(cal["y1"])))
    xpix = np.arange(int(cal["x0"]), int(cal["x1"]) + 1)
    ypix = _column_centroid(mask, ymin, ymax, xpix)

    keep = np.isfinite(ypix)
    e, cv = _pixel_to_data(xpix[keep].astype(float), ypix[keep], cal)
    order = np.argsort(e)
    e, cv = e[order], np.clip(cv[order], 0.0, None)
    e, cv = _hampel_filter(e, cv)
    return e, cv


def save_curve(path, e, c, header):
    arr = np.column_stack([e, c])
    np.savetxt(path, arr, fmt="%.4f %.3f", header=header, comments="# ")
    print(f"  saved {len(e):4d} pts -> {path}")


def save_data_with_err(path, e, c, err, header):
    arr = np.column_stack([e, c, err])
    np.savetxt(path, arr, fmt="%.4f %.3f %.3f", header=header, comments="# ")
    print(f"  saved {len(e):4d} pts -> {path}")


def _parse_bbox(s):
    parts = [int(p.strip()) for p in re.split(r"[ ,]+", s.strip()) if p.strip()]
    if len(parts) != 4:
        raise argparse.ArgumentTypeError("bbox needs 4 ints: x0,y0,x1,y1")
    return tuple(parts)


def _infer_pulse(png_path):
    m = re.search(r"(\d{4,6})", Path(png_path).name)
    return m.group(1) if m else "unknown"


def main(argv=None):
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("png", help="path to the KM14 spectrum PNG to extract from")
    p.add_argument("--pulse", default=None,
                   help="pulse number for filename prefix (default: inferred from PNG name)")
    p.add_argument("--outdir", default=None,
                   help="output directory (default: same as input PNG)")
    p.add_argument("--legend-bbox", type=_parse_bbox, default=DEFAULT_LEGEND_BBOX,
                   help=f"x0,y0,x1,y1 pixel bbox to suppress (the in-plot legend); "
                        f"default {DEFAULT_LEGEND_BBOX}. Pass --legend-bbox 0,0,0,0 "
                        f"to disable.")
    p.add_argument("--curves", nargs="+",
                   default=["data", "th", "bt", "scatt", "total"],
                   choices=list(CURVES.keys()),
                   help="which curves to extract (default: all)")
    p.add_argument("--preview", default=None,
                   help="optional output PNG visualising the extracted curves "
                        "on top of the input image (for sanity-checking)")
    args = p.parse_args(argv)

    png_path = Path(args.png).expanduser()
    if not png_path.exists():
        print(f"ERROR: {png_path} not found", file=sys.stderr)
        return 1

    pulse = args.pulse or _infer_pulse(png_path)
    outdir = Path(args.outdir).expanduser() if args.outdir else png_path.parent
    outdir.mkdir(parents=True, exist_ok=True)

    bbox = args.legend_bbox
    if bbox == (0, 0, 0, 0):
        bbox = None

    img = _load_image(png_path)
    print(f"PNG     : {png_path}  ({img.shape[1]}x{img.shape[0]} px)")
    print(f"Pulse   : {pulse}")
    print(f"Outdir  : {outdir}")
    print(f"Legend  : suppressed bbox = {bbox}")

    rgb = img[..., :3]
    results = {}
    data_err = None
    for name in args.curves:
        if name == "data":
            # Central value from the marker disk where visible; cap-midpoint
            # fallback recovers points hidden under a fit curve (see
            # _extract_data_points). Error bars from the cap extent.
            e, c, err, from_disk = _extract_data_points(rgb, PNG_CAL, legend_bbox=bbox)
            ndisk = int(from_disk.sum()); nrec = int((~from_disk).sum())
            hdr = (f"KM14 data points extracted from {png_path.name}\n"
                   f"# pulse {pulse}\n"
                   f"# columns: E_dep_MeV  counts")
            save_curve(outdir / f"{pulse}_KM14_data.txt", e, c, hdr)
            hdr_err = (f"KM14 data points + error bars extracted from {png_path.name}\n"
                       f"# pulse {pulse}; err = half the error-bar length (1-sigma)\n"
                       f"# value: marker disk where visible, else cap midpoint\n"
                       f"# columns: E_dep_MeV  counts  err")
            save_data_with_err(outdir / f"{pulse}_KM14_data_err.txt", e, c, err, hdr_err)
            print(f"    ({ndisk} from marker disk, {nrec} recovered from caps)")
            results[name] = (e, c)
            data_err = (e, c, err, from_disk)
        else:
            e, c = extract_curve(img, PNG_CAL, CURVES[name], legend_bbox=bbox)
            hdr = (f"KM14 {name} curve extracted from {png_path.name}\n"
                   f"# pulse {pulse}\n"
                   f"# columns: E_dep_MeV  counts")
            save_curve(outdir / f"{pulse}_KM14_{name}.txt", e, c, hdr)
            results[name] = (e, c)

    if args.preview:
        _save_preview(img, results, Path(args.preview).expanduser(),
                      legend_bbox=bbox, data_err=data_err)
    return 0


def _save_preview(img, results, out_path, legend_bbox=None, data_err=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    c = PNG_CAL
    epx = lambda x: c["e0"] + (x - c["x0"]) * (c["e1"] - c["e0"]) / (c["x1"] - c["x0"])
    cpx = lambda y: c["c0"] + (y - c["y0"]) * (c["c1"] - c["c0"]) / (c["y1"] - c["y0"])
    H, W = img.shape[:2]
    extent = [epx(0), epx(W), cpx(H), cpx(0)]
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.imshow(img, extent=extent, aspect="auto", zorder=0)
    styles = {"data": ("magenta", "o", "Data extracted"),
              "th": ("orange", "--", "Th extracted"),
              "bt": ("blue", ":", "B-th extracted"),
              "scatt": ("green", "-.", "Scatt extracted"),
              "total": ("red", "-", "Total extracted")}
    for name, (e, cv) in results.items():
        col, ls, lab = styles.get(name, ("magenta", "-", name))
        if name == "data":
            if data_err is not None:
                ed, cd, err, from_disk = data_err
                d = from_disk
                ax.errorbar(ed[d], cd[d], yerr=err[d], fmt="o", color=col, ms=4.0,
                            markeredgecolor="white", markeredgewidth=0.4,
                            ecolor=col, elinewidth=0.8, capsize=2,
                            label=lab, zorder=10)
                if (~d).any():
                    ax.errorbar(ed[~d], cd[~d], yerr=err[~d], fmt="s", color="cyan",
                                ms=5.0, markeredgecolor="black", markeredgewidth=0.5,
                                ecolor="cyan", elinewidth=0.9, capsize=2,
                                label="Recovered (under curve)", zorder=11)
            else:
                ax.plot(e, cv, "o", color=col, ms=4.0, alpha=1.0, label=lab,
                        markeredgecolor="white", markeredgewidth=0.4, zorder=10)
        else:
            ax.plot(e, cv, color=col, lw=1.6, ls=ls, label=lab, zorder=5)
    if legend_bbox is not None:
        bx0, by0, bx1, by1 = legend_bbox
        # bbox is in pixel coords; convert to data coords for the overlay.
        ex0, ex1 = epx(bx0), epx(bx1)
        cy0, cy1 = cpx(by0), cpx(by1)
        ax.plot([ex0, ex1, ex1, ex0, ex0], [cy0, cy0, cy1, cy1, cy0],
                color="grey", lw=0.8, ls="--", alpha=0.8, zorder=4,
                label=f"legend bbox (cuts counts > {min(cy0, cy1):.0f})")
    ax.set_xlim(7.2, 9.5); ax.set_ylim(0, 200)
    ax.set_xlabel(r"$E_{dep}$ (MeV)"); ax.set_ylabel("Counts")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    print(f"Preview : {out_path}")


if __name__ == "__main__":
    raise SystemExit(main())
