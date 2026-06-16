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
    <pulse>_KM14_data.txt    # black markers (data points)
    <pulse>_KM14_th.txt      # orange dashed
    <pulse>_KM14_bt.txt      # blue dotted
    <pulse>_KM14_scatt.txt   # green dash-dot
    <pulse>_KM14_total.txt   # red solid

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
DEFAULT_LEGEND_BBOX = (340, 30, 730, 240)

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


def _erode_horizontal(mask, width=2):
    """Keep only pixels with `width` black neighbours on each side at the same
    row. Suppresses vertical error-bar lines (1 col wide) and most axis text
    while leaving filled marker disks intact."""
    out = mask.copy()
    for dx in range(1, width + 1):
        shifted_left = np.zeros_like(mask); shifted_left[:, dx:] = mask[:, :-dx]
        shifted_right = np.zeros_like(mask); shifted_right[:, :-dx] = mask[:, dx:]
        out &= shifted_left & shifted_right
    return out


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


def _column_biggest_cluster_center(mask, ymin, ymax, xpix, min_size=3):
    """For each x-column, return the centre y of the largest contiguous run of
    True pixels (clusters separated by gaps > 1). Suitable for data markers
    after horizontal erosion: error-bar end caps shrink to small clusters,
    the marker survives as the largest cluster."""
    out = np.full(xpix.size, np.nan, dtype=float)
    for i, xp in enumerate(xpix):
        col = mask[ymin:ymax + 1, xp]
        idx = np.where(col)[0]
        if idx.size < min_size:
            continue
        gaps = np.where(np.diff(idx) > 1)[0]
        starts = np.concatenate([[0], gaps + 1])
        ends = np.concatenate([gaps + 1, [idx.size]])
        sizes = ends - starts
        b = int(np.argmax(sizes))
        if sizes[b] < min_size:
            continue
        cluster = idx[starts[b]:ends[b]]
        out[i] = ymin + float(np.mean(cluster))
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
    """Run the extraction for one curve. Returns (E_dep_MeV, counts) arrays."""
    rgb = img[..., :3]
    if params["color_name"] == "black":
        # Horizontally erode so vertical error-bar lines (1 col wide) and axis
        # tick text disappear; only the wider marker disks (and end caps,
        # rejected later by cluster-size selection) survive.
        mask = _erode_horizontal(_black_mask(rgb), width=2)
    else:
        hue_lo, hue_hi = params["hue"]
        mask = _hue_mask(rgb, hue_lo, hue_hi,
                         params["sat_min"], params["val_min"])
    mask = _apply_legend_bbox(mask, legend_bbox)

    ymin, ymax = sorted((int(cal["y0"]), int(cal["y1"])))
    xpix = np.arange(int(cal["x0"]), int(cal["x1"]) + 1)
    if params["marker_like"]:
        ypix = _column_biggest_cluster_center(mask, ymin, ymax, xpix)
    else:
        ypix = _column_centroid(mask, ymin, ymax, xpix)

    keep = np.isfinite(ypix)
    e, cv = _pixel_to_data(xpix[keep].astype(float), ypix[keep], cal)
    order = np.argsort(e)
    e, cv = e[order], np.clip(cv[order], 0.0, None)
    # Hampel rejection only on the smooth curves -- data markers are sparse
    # and a marker fluctuating off the trend is real signal, not an outlier.
    if not params["marker_like"]:
        e, cv = _hampel_filter(e, cv)
    return e, cv


def save_curve(path, e, c, header):
    arr = np.column_stack([e, c])
    np.savetxt(path, arr, fmt="%.4f %.3f", header=header, comments="# ")
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

    results = {}
    for name in args.curves:
        e, c = extract_curve(img, PNG_CAL, CURVES[name], legend_bbox=bbox)
        out_path = outdir / f"{pulse}_KM14_{name}.txt"
        hdr = (f"KM14 {name} curve extracted from {png_path.name}\n"
               f"# pulse {pulse}\n"
               f"# columns: E_dep_MeV  counts")
        save_curve(out_path, e, c, hdr)
        results[name] = (e, c)

    if args.preview:
        _save_preview(img, results, Path(args.preview).expanduser())
    return 0


def _save_preview(img, results, out_path):
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
    styles = {"data": ("k", "o", "Data extracted"),
              "th": ("orange", "--", "Th extracted"),
              "bt": ("blue", ":", "B-th extracted"),
              "scatt": ("green", "-.", "Scatt extracted"),
              "total": ("red", "-", "Total extracted")}
    for name, (e, cv) in results.items():
        col, ls, lab = styles.get(name, ("magenta", "-", name))
        if name == "data":
            ax.plot(e, cv, "o", color=col, ms=2.5, alpha=0.7, label=lab, zorder=5)
        else:
            ax.plot(e, cv, color=col, lw=1.6, ls=ls, label=lab, zorder=5)
    ax.set_xlim(7.2, 9.5); ax.set_ylim(0, 200)
    ax.set_xlabel(r"$E_{dep}$ (MeV)"); ax.set_ylabel("Counts")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    print(f"Preview : {out_path}")


if __name__ == "__main__":
    raise SystemExit(main())
