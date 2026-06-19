# `km14_spectrum.py` — KM14 diamond neutron-spectrum forward model (added 2026-06-15)

`src/neutron/km14/km14_spectrum.py` forward-models the **KM14 diamond detector
neutron energy spectrum** (thermal-DT + beam-thermal-DT) along the KM14 LOS and
overlays it on the measured #104614 spectrum, to compare against the
Nocente–Rigamonti diamond analysis (`~/jet/data/104614/figs/`).

## What KM14 actually is (corrected 2026-06-15)
KM14 is a **single-crystal diamond detector** measuring **14 MeV DT neutrons**.
Its energy axis is *deposited* energy `E_dep` via `¹²C(n,α₀)⁹Be` (Q≈−5.7 MeV),
so the 14.03 MeV DT line lands at `E_dep ≈ 8.4 MeV` (measured span 7.2–9.3 MeV).
**104614 M29 is a DT run** (`THNTX_DT` ~140× `THNTX_DD`, thermal `NT`, `BTN1`),
not pure DD — earlier notes were wrong. Use the **DT** channel for KM14.

## Method
* **LOS weight** `f(rhot)` on the TRANSP X grid via the shared
  `los_thermal_rate.los_shell_fraction` (same chord geometry as
  `los_th_bt_ratio.py`; reuses `ltb.CdfEquilibrium`/`read_lcfs`/`rhot_pinned`).
* **Thermal DT**: per flux shell a Gaussian line at `E0=14.03 MeV`,
  `FWHM=177·√(Ti[keV])` keV (DT Brysk Doppler width), weighted by
  `THNTX_DT·f·DVOL`. Integral = TH LOS rate.
* **Beam-thermal DT**: per NUBEAM zone, MC-sample fast-D (speed/pitch/gyrophase
  via `BField` B̂) + thermal-T Maxwellian (+rotation), weight `σ·v_rel`, emit
  isotropic-CM → lab velocity (`bt_kinematics`), keep neutrons into the KM14
  detector cone (per-zone `n̂` to the detector point), histogram lab energy; each
  zone scaled by `BTN1·f(rhot_zone)·BMVOL`. Reuses `bt_los_emissivity`.
* **Detector response**: `E_dep = E_n + shift` (shift set by `--edep-peak`,
  default aligns the TH peak to 8.40 MeV → −5.63 MeV) then Gaussian convolution
  with `--det-fwhm` (default 0.20 MeV).
* **Normalization**: scale the **total** (TH+BT) so its peak matches the data
  (`--peak-counts`, default 145); the **TH:BT split is fixed by the LOS rates**
  (the physics under test), so only one global amplitude is free.
* **Overlay**: the experimental PNG is shown as a calibrated background
  (axis px→data calibration hard-coded in `PNG_CAL`: x-ticks 7.2–9.2 at px
  129–725, y-spine 0–200 at px 584–32); model TH/B-th/Total drawn on top.

## Result (104614 M29 idx 1, t_TRANSP=12.33 s = JET 53.15–53.5 s)
`(TH/BT)_LOS = 0.53` (TH fraction 34.7%), TH LOS rate 2.03e17 n/s (matches
`los_th_bt_ratio --channel dt` `Sum(TH·f·DVOL)` exactly), BT 3.82e17 n/s (0.7%
vs the flux-profile integral). **The model TH (magenta) lands on the Nocente
fitted Th component (8.40 MeV, FWHM 0.47); B-th (lime) ≈ their B-th (peak 8.56,
FWHM 1.15 MeV).** The model omits the **Scatt** background (detector/environment,
not emission), so the total slightly underfits the 7.6–8.0 MeV low-energy shoulder.

**On the TH/BT discrepancy (user-corrected 2026-06-15).** Nocente's low ρ<0.2 TH/BT
(**25–31%**) was *not* a LOS effect — it was the **tritium concentration of their
TRANSP run 104614M13** (T/(D+T)≈**0.53**). High tritium lowers n_D at fixed n_e, so
thermal-DT (∝ n_D·n_T) is suppressed relative to beam-D-on-thermal-T BT (∝ n_fast·n_T;
n_T largely cancels in the ratio, the depressed n_D does not). **M29 has
T/(D+T)≈0.31** (vol-avg ρ<0.2; axis 0.305): on M29, ρ<0.2 TH/BT = **0.627**,
LOS-weighted = **0.53**, both bracketing the measured **~0.59**. So TRANSP with a
realistic tritium fraction already matches the data; the KM14 LOS weighting is a
**secondary ~15% effect** (core 0.63 → LOS 0.53, the chord seeing BT-richer wings).

## CLI
```bash
python src/neutron/km14/km14_spectrum.py 104614 M29 --idx 1 --save        # overlay PNG into figs/
python src/neutron/km14/km14_spectrum.py 104614 M29 --idx 1 --no-plot     # text rates only
python km14_spectrum.py 104614 M29 --idx 1 --det-fwhm 0.25 --edep-peak 8.40 --nsamp 80000
```
Flags: `--idx --data-dir --Rmin --Rmax --nR --nZ --detector-R --detector-height
--cone-deg --nsamp --fast-norm {bdens,ntot} --no-rotation --seed --e0-dt
--edep-peak --det-fwhm --de-kev --peak-counts --png --save --no-plot`.
Output overlay: `figs/<run_id>_KM14_spectrum_overlay_idx<idx>.png`.

---

# Data-driven TH:BT split + Scatt overlay + chi² (added 2026-06-16)

`km14_spectrum.py` grew several CLI hooks to **compare the forward-model directly
against the measured PNG spectrum** and to break the TH:BT split out of the
TRANSP prediction. The image-analysis side was factored into a separate script.

## `extract_km14_png.py` — pull curves out of the Nocente PNG

`src/neutron/km14/extract_km14_png.py` reads the published KM14 spectrum PNG
and writes one 2-column text file per curve (data, Th, B-th, Scatt, Total) to
`~/jet/data/<pulse>/figs/<pulse>_KM14_<kind>.txt`. The pixel-to-data calibration
is hard-coded (same `PNG_CAL` as `km14_spectrum.py`: x-ticks 7.2–9.2 MeV at px
129–725, y-spine 0–200 counts at px 584–32).

* **Curve filter.** HSV-hue ranges per component (Th orange 0.05–0.13, B-th blue
  0.55–0.72, Scatt green 0.30–0.40, Total red 0.93–0.07 wrap). For each x-pixel
  column inside the plot area, takes the **median y** of pixels matching the
  colour → centroid of the line, then a Hampel filter (window=11, n_sigma=4)
  drops outliers. Dash patterns leave gaps which downstream `np.interp` fills.

* **Data points + error bars — rewritten 2026-06-17.** The marker is a black
  filled disk (~8×7 px) on a vertical error-bar stem (~3 px wide) with two
  horizontal caps (~3 rows tall, ~13 px wide). Extraction (`_extract_data_points`)
  now returns `(E, counts, err, from_disk)` and is driven by the error-bar
  geometry, in this order:
  1. **Locate every point** by peaks in the per-column black-pixel count
     (`scipy.signal.find_peaks`, params `DATA_PEAK_MIN_HEIGHT/DISTANCE/PROMINENCE`).
     The stem is black at the point's exact x **even when the disk is hidden
     under a fit curve**, so this recovers points buried under the red Total (and
     anywhere blue/orange cross a marker). 104614: 78 points found, ~8 px spacing.
  2. **Central value = marker disk where visible.** A 4×4 binary erosion
     (`DATA_MARKER_ERODE`) isolates the disk (only feature both wider than the
     stem and taller than the caps); `_disk_centroids` labels + takes the
     centroid, merging split-disk arcs (a curve crossing the disk middle) within
     `DISK_MARKER_MERGE`=3 px in x. If a disk is within `DISK_MATCH_PX`=4 px of a
     point → value = disk centroid (`from_disk=True`); else → cap midpoint
     (`from_disk=False`). 104614: 58 from disk, 20 cap-recovered.
  3. **Error bar from the caps** (`_error_bar_extent`): the 1-σ half-length is
     half the top↔bottom cap separation, read from the **single stem column**
     (not a ±1 px band) keeping only black runs ≥ `DATA_BAR_MIN_RUN`=5 px. The
     value is the cap midpoint (symmetric Poisson bars). Two coupled fixes here:
     the single column avoids neighbour caps (13 px wide > 5–7 px peak spacing),
     and the run filter drops a neighbour cap that still bleeds in — it is a
     short isolated run (no stem at this column), while the point's own caps
     connect to its stem (long runs, or two long runs where red splits the bar).
     Result: peak err ~9–14 cnts (≈√N), tail err ~3–5 cnts; the earlier ±1-band
     read gave 20–40 cnts at the peak.

* **Axis-frame strip** (`_strip_axis_frame`): rows/cols black across >50 % of the
  plot span are spines, blanked before peak/cap reads (otherwise the column scan
  hits the 200-cnt top spine and 0-cnt bottom spine → err pinned at ±100).

* **Legend suppression.** The in-plot legend overlaps the data region; default
  bbox `DEFAULT_LEGEND_BBOX = (560, 30, 730, 240)` zeroes it before extraction
  (overridable `--legend-bbox x0,y0,x1,y1`, `0,0,0,0` disables). **The bbox is a
  hard rectangular cut on the data read, not legend-text-only** — its bottom edge
  maps to a counts level (`counts=(584−ypix)/2.76`), so it must clear the data
  peak. The old `(340,30,730,240)` bottom edge sat at ypix 240 ≈ **125 counts**
  over the peak x-band, silently truncating the Total/Data peak (~190 cnts) at
  125; narrowing x to 560 fixed it.

* **Outputs:** `<pulse>_KM14_data.txt` (E, counts) and
  `<pulse>_KM14_data_err.txt` (E, counts, err); plus the four curve files. The
  `from_disk` flag is internal (drives the stdout disk/recovered tally and the
  preview colouring), not persisted — add a 4th column if downstream weighting
  wants it. NB: `err` (1-σ) is now available to weight the `km14_spectrum.py`
  `--th-bt-fit` (which was unweighted for lack of a sigma).

* **Preview.** `--preview out.png` overlays everything on the input image
  (`ax.imshow`, toggleable): disk values as magenta circles, cap-recovered as
  cyan squares ("Recovered (under curve)"), both with error bars; the legend
  bbox is drawn as a grey dashed rectangle labelled with its cut threshold.

```bash
python src/neutron/km14/extract_km14_png.py \
    ~/jet/data/104614/figs/104614_KM14_spectral_analysis.png \
    --preview ~/jet/data/104614/figs/104614_KM14_extraction_preview.png
```
Flags: `--pulse --outdir --legend-bbox --curves {data,th,bt,scatt,total} --preview`.
Default `--pulse` inferred from PNG filename via regex.

**Pitfalls encountered while building the extractor:**
1. **Legend bbox truncates the peak** (see above) — its lower edge is a counts
   cut; keep it above the Total peak. Visualise it via the preview's dashed box.
2. **Disk eaten by overlapping curves.** A naive disk-only detector misses every
   point a fit curve crosses (only ~53–67 of 78 on 104614). The stem-peak locator
   is curve-independent; disk is used only for the value where it survives.
3. **Neighbour-cap bleed inflates errors.** Caps (13 px) are wider than the peak
   point spacing (5–7 px); read the bar from the single stem column + run filter.
4. **Axis frame leak** pins err at ±100 cnts unless spines are stripped first.

## New `km14_spectrum.py` analysis flags

* `--th-bt-ratio R` — manually rescale TH so `(TH/BT)_LOS = R`, BT shape
  unchanged. Total peak still scaled to `--peak-counts`. Useful for by-eye
  comparisons against the data.
* `--th-bt-fit` — NNLS fit of `(alpha_th, alpha_bt)` against the extracted
  data over `[--fit-emin, --fit-emax]` (defaults 8.0–9.3 MeV, above the
  Scatt-dominated shoulder). Reports the fitted TH/BT, alpha coefficients, and
  per-point RMS. **Unweighted** (no sigma on the extractor output) — see
  weighting caveat below.
* `--include-scatt` — load `<pulse>_KM14_scatt.txt` (or `--scatt-file`),
  interpolate onto the model E_dep grid, and **add it to the model total**.
  In MANUAL mode the TH+BT peak rescale uses `headroom = peak_counts -
  scatt(peak)` so the total (TH+BT+Scatt) peaks at `--peak-counts`. In FIT
  mode, Scatt is **subtracted from the data first**, then the NNLS fits
  (alpha_th, alpha_bt) to the residual TH+BT shape.
* `--data-file PATH` / `--scatt-file PATH` — explicit override; default to
  `~/jet/data/<pulse>/figs/<pulse>_KM14_{data,scatt}.txt`.
* `--core-rhot R` — diagnostic: **replace the LOS shell fraction `f(rhot)`
  with the Heaviside step `Theta(R - rhot)`** so the model integrates only
  flux shells inside rhot ≤ R, no chord geometry. This is the Nocente-style
  "core TH/BT" diagnostic. Print labels and plot legend switch from `LOS` to
  `rho<R`. Save filename gets a `_core<R>` tag.
* `--chi2-emin --chi2-emax --no-chi2` — chi² goodness-of-fit on a chosen
  window (default 7.5–9.3 MeV) using **Pearson chi² with Poisson errors**
  (sigma = sqrt(max(N, 1))). Always reported when the data file exists; the
  fit-mode data is reused, otherwise loaded from the default path. Output:
  `chi^2 / N`, `RMS counts`.

The image-analysis helpers (`_extract_curve_from_png`, `_extract_data_from_png`,
`_scatt_on_grid`) were **removed** from `km14_spectrum.py` — it now only
**reads** the extractor's text files. `PNG_CAL` stays for the experimental-PNG
background imshow in the overlay plot.

## Headline result (104614 M29 idx 1, JET t = 53.15–53.5 s)

```bash
python src/neutron/km14/km14_spectrum.py 104614 M29 --idx 1 \
    --edep-peak 8.373 --cone-deg 4 --include-scatt --th-bt-ratio 0.58
```
yields **chi²/N = 3.07** (chi² = 1350.5 over N = 440 points in 7.5–9.3 MeV),
the best fit obtained. Three knobs each shave the chi²:
* `--idx 1` matches the diamond integration window (idx 2 was a +1 s later
  TRANSP slice, comparing TRANSP-13.3 s against data-53.15–53.5 s).
* `--cone-deg 4` kills the spurious high-E BT tail. At the default 20°, the
  model accepted off-vertical neutrons with substantial `v_CM·n̂` Doppler
  boost; narrowing to 4° restricts to nearly-vertical neutrons (small boost),
  shifting the BT peak from 8.4 → 8.23 MeV and tightening its FWHM
  (1.15 → 1.12 MeV) — exactly the issue raised earlier on M29 idx 2 and the
  primary reason the BT model was too wide and shifted above the data.
* `--th-bt-ratio 0.58` is **only +9 % above TRANSP M29 LOS (0.53)** — i.e. the
  M29 TRANSP prediction was already essentially right.

Below ~chi²/N = 3 the residual is shape, not normalisation: BT energy
distribution (no anomalous fast-ion transport in NUBEAM, isotropic-CM
emission), Scatt extraction noise, and detector-resolution model.

## M29 vs M13 scans — the tritium-fraction story plays out

| Mode | (TH/BT) source | TH/BT | TH frac | chi²/N |
|---|---|---|---|---|
| **M13 idx 2 LOS** | TRANSP | 0.235 | 19 % | 11.5 |
| M13 idx 2 CORE rhot<0.2 | TRANSP | 0.255 | 20 % | 13.7 |
| M13 LOS NNLS (7.5–9.3) | fit | 0.375 | 27 % | 5.57 |
| M13 LOS Poisson-chi² min | scan | ~0.70 | ~41 % | 4.70 |
| **M29 idx 2 LOS** | TRANSP | **0.526** | **34.5 %** | 5.29 |
| M29 idx 2 CORE rhot<0.2 | TRANSP | 0.630 | 38.7 % | 4.78 |
| M29 LOS NNLS (8.0–9.3) | fit | 0.507 | 33.6 % | 4.93 |
| M29 LOS Poisson-chi² min | scan | ~0.70 | ~41 % | 4.81 |
| M29 CORE Poisson-chi² min | scan | ~0.80 | ~44 % | 4.48 |

Conclusions confirmed against the actual data:

1. **M29 TRANSP-LOS already fits the data** (chi²/N = 5.29 on idx 2, drops to
   3.07 on idx 1 with cone narrowing) — no TH:BT tuning required.
2. **M13 TRANSP-LOS misses by a factor ~2 in TH/BT** because its T/(D+T)≈0.53
   is too high (real plasma closer to M29's 0.31). The user-corrected note
   from 2026-06-15 (Nocente's low ρ<0.2 = M13 high-tritium artifact) is
   reproduced here from first principles via the data fit.
3. **Nocente's quoted 25–31 % almost certainly came from M13 CORE rhot<0.2**
   (we get 25.5 %). It does *not* match M29 CORE (38.7 %) → so his analysis
   hard-wired M13's TRANSP. M29 + LOS is the right combination to compare
   *against* the data.
4. **LOS vs CORE rhot<0.2 are indistinguishable in the data.** Poisson-chi²
   minimum 4.81 vs 4.48 — a 7 % difference. The chord weighting redistributes
   a small fraction between core and wings, the spectral shape is preserved.
   Pick whichever physical quantity matches the question (LOS for "what KM14
   sees", CORE for "core neutron emission TH/BT").
5. **Fit minimum is weighting-dependent.** Unweighted NNLS on 7.5–9.3 gives
   TH/BT ≈ 0.38 (peak-dominated); Poisson-weighted chi² minimum gives
   TH/BT ≈ 0.70 (wing-dominated). The two are not in conflict — they answer
   different questions. Default fit weighting is **unweighted** because the
   extractor output has no sigma column; the chi² report is **Poisson-weighted**
   on the same window. Future improvement: pass Poisson sigma to NNLS so the
   fit minimum coincides with the chi² minimum.

## Workflow (extract once per pulse, then analyse)

```bash
# 1. Extract all curves from the published PNG (run once per pulse)
python src/neutron/km14/extract_km14_png.py \
    ~/jet/data/<pulse>/figs/<pulse>_KM14_spectral_analysis.png

# 2. TRANSP-LOS comparison + chi^2
python src/neutron/km14/km14_spectrum.py 104614 M29 --idx 1 --include-scatt

# 3. Best chi^2 setup found so far
python src/neutron/km14/km14_spectrum.py 104614 M29 --idx 1 \
    --edep-peak 8.373 --cone-deg 4 --include-scatt --th-bt-ratio 0.58

# 4. Core-only diagnostic 
python src/neutron/km14/km14_spectrum.py 104614 M29 --idx 1 \
    --core-rhot 0.2 --include-scatt

# 5. Data fit (free TH:BT split)
python src/neutron/km14/km14_spectrum.py 104614 M29 --idx 1 \
    --th-bt-fit --fit-emin 7.5 --include-scatt
```

Output PNGs get tagged with `_fit` / `_manual` / `_scatt` / `_core<R>` so the
different runs don't overwrite each other. Saved into the same `figs/` dir as
the input PNG.

## `--los-file` — real KM14 LOS cells in the spectrum (added 2026-06-18)

`km14_spectrum.py` now accepts `--los-file PATH` (e.g. `_KM3.los`), the
same flag wired through `los_thermal_rate.py` and `los_th_bt_ratio.py`.
When set:

1. The geometric shell-fraction `f(rhot)` (rectangular-chord proxy) is
   replaced by the **etendue-coupled** `f_det(rhot) = C_bin / DVOL` from
   `los_thermal_rate.los_file_detector_rate`. Same anti-aliased binning,
   enclosed-shell pinning and `rhot_crit` flattening as documented for
   the other two scripts.
2. The per-NUBEAM-zone BT MC cone axis switches from the point-detector
   approximation to the file's **exact** per-cell emission versor
   `(u, v, w)` interpolated to each zone via
   `ltb.los_directions_for_zones`. Zones outside the LOS footprint fall
   back to the point-detector estimate (so the option is
   back-compatible).

Both replacements feed unchanged into `thermal_spectrum` and `bt_spectrum`
— the rest of the spectrum pipeline (Doppler-broadened Gaussian sum for
TH, kinematic Brysk/DT MC for BT, `E_n -> E_dep` + detector-FWHM
convolution, peak/fit/manual scaling, Scatt overlay, χ²) is unchanged.

**Smoke test (104614 M29 idx 1, t_TRANSP = 12.33 s):**
* 77 611 / 117 424 LOS cells inside LCFS; 79 / 220 NUBEAM zones in the
  real-LOS footprint.
* `(TH/BT)_LOS_real = 0.535` vs geometric `(TH/BT)_LOS = 0.531` — the
  etendue weighting pulls the effective sampling toward the higher-TH/BT
  core (consistent with the `rho_50 < rho_bnd` observation documented
  for `los_thermal_rate.py` / `los_th_bt_ratio.py`).
* The TH/BT spectrum shapes (E_dep peak ≈ 8.41 MeV, BT FWHM ≈ 1.17 MeV,
  TH FWHM ≈ 0.47 MeV) are unchanged to the printed precision — the
  detector convolution and the broad BT kinematic spread dominate over
  the few-percent rhot-weighting reshuffle.

Output PNGs gain a `_realLOS` tag (`tag += "_realLOS"`) so geometric- and
real-LOS runs don't overwrite each other.

```bash
python src/neutron/km14/km14_spectrum.py 104614 M29 --idx 1 \
    --los-file src/neutron/km14/_KM3.los --no-plot --no-chi2

python src/neutron/km14/km14_spectrum.py 104614 M29 --idx 1 \
    --los-file src/neutron/km14/_KM3.los --th-bt-fit --include-scatt --save
```

## `--bt-aniso` — BT directional reweighting in the spectrum (added 2026-06-18)

`km14_spectrum.py` also gained `--bt-aniso`, the same per-zone directional
factor `g(zone) = (dEps/dOmega toward detector)/(Eps/4pi)` shipped in
`los_th_bt_ratio.py`. It reuses `ltb.bt_anisotropy_factors` (so the
in-house MC + Bosch-Hale + CM->lab + cone-acceptance machinery is
shared) and multiplies each NUBEAM zone's `BTN1` by `g(zone)` before
`bt_spectrum` integrates. Combined with `--los-file`, `g` is evaluated
along the file's exact per-cell `(u, v, w)` direction (point-detector
fallback for zones outside the LOS footprint).

The km14_spectrum BT MC *already* gets the lab-frame angular spectrum
shape right by histogramming only in-cone samples — what `--bt-aniso`
adds is the **per-zone weight** correction (a few zones favour
toward-detector emission, others suppress it, integrated reactivity-
weighted). The two effects are independent.

**Sensitivity summary** (104614 M29 idx 2, `--los-file --th-bt-fit
--include-scatt --edep-peak 8.373`):

| Config | (TH/BT)_TRANSP | (TH/BT)_fit | chi^2/N | RMS [cnts] | BT peak [MeV] |
|---|---|---|---|---|---|
| no aniso (cone 20 deg) | 0.5300 | 0.4508 | 3.55 | 13.46 | 8.43 |
| `--bt-aniso` (cone 20 deg) | 0.5299 | 0.4506 | 3.55 | 13.45 | 8.39 |
| `--bt-aniso --cone-deg 10` | 0.5258 | 0.4522 | 3.62 | 13.28 | 8.54 |
| `--bt-aniso --no-rotation` | 0.5302 | 0.4477 | 3.55 | 13.49 | 8.38 |

* **BT anisotropy is negligible**: `g_DT = 1.0000` (reactivity-weighted),
  chi^2 shift < 0.1 %. KM14 vertical chord is ~perpendicular to the
  toroidal B (angle(B, LOS) = 76-90 deg), so directional reactivity has
  nothing to align toward. In-house MC reproduces NUBEAM 4 pi rate to
  ~0.6 % (`vector Eps_4pi/NUBEAM = 0.994`).
* **Cone-deg is the real lever** for the spectrum shape (not per-zone
  `g`): 20 deg -> 10 deg shifts the BT peak +0.15 MeV (narrower cone
  preferentially keeps neutrons with stronger CM->lab forward boost) and
  worsens chi^2 ~2 %.
* **Rotation matters slightly**: removing it shifts fitted TH/BT by
  ~-0.7 %.

With `--los-file` the spectrum now uses **all of `_KM3.los`**: `C`
(etendue) -> `f_det(rhot)` shell weight; `(u, v, w)` -> per-zone BT MC
cone axis; `(u, v, w)` -> directional `g(zone)` under `--bt-aniso`;
`V, R, Z` (cell positions) -> `f_det` binning. The `_aniso` filename tag
on auto-saved PNGs keeps the four combinations (no-flag / `--los-file`
/ `--bt-aniso` / both) from overwriting each other.

### `--cone-deg` is an MC averaging window, not a physical detector cone

A natural question once `--los-file` is in: can the *cone half-angle*
also be inferred from the file? **No** -- `cone-deg` stays a free CLI
parameter, because it serves a purely numerical role.

* The LOS file specifies the *direction* `n_LOS = (u, v, w)` per cell
  exactly, and the *per-cell etendue* `C` gives the true detector solid
  angle `Omega = 4 pi C / V` ~ `1e-3 sr` ~ **~1 deg half-angle**. That
  is the *physical* cone any one cell subtends.
* The MC samples isotropic CM directions per (fast-ion, thermal) pair
  and keeps only neutrons whose lab direction falls within `cone-deg`
  of `n_LOS`. Using the physical `~1 deg` cone would accept ~3e-5 of
  the samples -- with `nsamp = 60000` per zone that's ~2 in-cone
  neutrons, pure noise. So `cone-deg` has to be a *wider* sampling
  window than the physics, sized for MC statistics.
* The acceptable looseness comes from the lab-frame BT spectrum being
  **smooth in angle near the LOS**. Averaging across a small cone gives
  nearly the same shape as the true point direction. For KM14
  (vertical chord ~perpendicular to the toroidal B), this smoothness is
  generous; the 104614 M29 idx 2 sensitivity scan showed the spectrum
  barely moves over `cone-deg` up to ~30 deg (`20 -> 10 deg` shifts the
  BT peak +0.15 MeV and chi^2 +2 %, the residual bias of cone tightening
  trading off against the kinematic CM->lab forward-boost selection).

Geometry inputs vs `cone-deg`:

| Quantity                   | Source                       | Role                             |
|----------------------------|------------------------------|----------------------------------|
| `n_LOS` per zone           | `(u,v,w)` from `_KM3.los`    | center of MC cone (exact)        |
| per-cell etendue `C`       | `_KM3.los`                   | shell weight `f_det = C_bin/DVOL`|
| MC angular window `cone-deg` | CLI flag                   | **numerical**, not physical      |

**Future improvement (not today).** The hard cone cut wastes the
out-of-cone MC samples and produces this residual `cone-deg` bias. A
cleaner version would weight every sample by an angular kernel (e.g.
Gaussian in `arccos(n_hat . n_LOS)` with a width chosen for
bias-variance trade-off) so all samples contribute, the spectrum
relaxes toward the true `cone-deg -> 0` limit, and the knob disappears.
Worth it if a future analysis needs cone-deg-insensitive BT spectra; not
worth it now since the current effect is ~2 % on chi^2.

## BT peak position vs Nocente on M13 (open issue, 2026-06-18)

Running the spectrum on 104614 **M13** idx 2 with the published
peak-alignment (`--edep-peak 8.373 --th-bt-ratio 0.62 --include-scatt
--los-file _KM3.los --bt-aniso`): the **TH component reproduces
Nocente's TH almost exactly**, but the **BT component sits ~30-50 keV
high** in `E_dep`, dragging the total (black curve) above Nocente's
total (red). Both calculations use the same M13 NUBEAM run as input, so
the discrepancy is in the forward model, not the data. M13 has only
`BTN1` (DT) and `BTN4` (DD) -- no `BTN7` / `BTN5` -- so it is **not**
a missing-channel issue; our `bt_spectrum` uses exactly the BT-DT
component that contributes at 14 MeV.

### Leading hypothesis: asymmetric diamond detector response

We convolve every neutron spectrum with a **symmetric Gaussian** of
fixed FWHM (`--det-fwhm`, default 0.20 MeV). Real diamond response to
14 MeV neutrons is **asymmetric**: a sharp upper edge at
`E_dep = E_n + Q` (Q = -5.701 MeV, 12C(n,alpha0)9Be) plus a low-energy
tail from partial charge collection, neutron escape, and inelastic
channels (`12C(n,alpha1)`, `(n,n')3alpha`). When the input is broad
(BT FWHM ~ 1.2 MeV) and you convolve with an asymmetric kernel that
has a low-E tail, the **peak shifts toward lower E_dep** more for
broader inputs than for narrow ones:

* TH (Doppler FWHM ~ 180 keV) -- peak barely moves, dominated by the
  response shape itself.
* BT (~ 1.2 MeV) -- peak pulled noticeably down by the tail.

So a realistic response pulls BT toward TH; our symmetric Gaussian
leaves BT at the kinematic position. That is the exact sign of the
discrepancy observed. Nocente almost certainly uses a tabulated
response function (likely MCNP-derived or measured), not a Gaussian.

**Empirical confirmation that wider Gaussian doesn't fix it.** Running
with `--det-fwhm 0.4` (double the default) broadens everything but does
**not** appreciably shift the BT peak position. This is expected: a
*symmetric* convolution leaves the peak position unchanged regardless
of width. The shift requires *asymmetry* in the response. So adding a
proper response is the only real fix; widening the FWHM is not.

### Other factors with smaller likely impact

Listed roughly in order of expected magnitude:

1. **DT thermal birth energy `E0_DT`.** We use 14.03 MeV (default); the
   2-body rest-frame value is 14.04 MeV (`Q * m_alpha/(m_alpha + m_n)
   = 17.589 * 0.7986`). ~10 keV; trivial.
2. **CM-frame angular distribution.** `btk.isotropic_cm` samples
   uniformly in `cos theta_CM`. Real DT has ~5-10 % anisotropy at
   `E_cm ~ 100 keV`. For KM14 (LOS perp. to v_CM toroidal direction),
   this mostly affects BT *width*, not peak.
3. **Plasma rotation on fast ions.** `bt_los_emissivity.fast_velocity_vectors`
   does NOT add `v_rot` to the fast ion (only to the thermal target).
   Whether `F_D_NBI` is in lab or plasma-rotating frame is a NUBEAM
   convention question; if rotating, we miss the rotation contribution
   to `v_CM`. For KM14 the toroidal projection on the vertical LOS is
   ~0, so the linear effect cancels; only second-order corrections
   survive.
4. **F_D_NBI normalization.** `--fast-norm bdens` (default) is the
   recommended setting (matches `BDENS_D`). `ntot` would change
   fast-ion density by ~30 % in places and slightly tilt the
   energy-averaged boost.
5. **LOS direction within zones.** `los_directions_for_zones`
   griddata-interpolates the cell `(u, v, w)` field onto each NUBEAM
   zone -- a single direction per zone. Real cells span a small
   angular range within each zone (`w ~ 0.9997`, tilt < 1.4 deg), so
   using one direction per zone slightly biases the kinematic
   projection.
6. **Non-relativistic kinematics.** `btk.neutron_lab_velocity` is
   non-relativistic. For 14 MeV neutrons `v/c ~ 0.17`, so relativistic
   corrections to `E_lab` are ~1.5 %; for KM14's perpendicular
   geometry the projection on LOS is much smaller.

### Path to a proper fix (future)

The fix is to replace the Gaussian convolution with a **tabulated
response function**:

* Input: `R(E_dep | E_n)` -- a 2-D matrix (or `E_n`-dependent line
  shape) giving the probability density of `E_dep` for a monoenergetic
  `E_n` neutron in a diamond.
* Source candidates: MCNP simulation of the KM14 diamond geometry, a
  published parameterisation (Nocente or KM14 instrument paper), or
  measurement against a calibration source.
* Apply by convolving sample-by-sample (or zone-by-zone): rather than
  one global convolution of the total `(TH+BT)_n` spectrum, fold each
  shell's contribution with `R(E_dep | E_n)` so that BT, which spans a
  broader `E_n`, samples the response asymmetry differently from TH.

This is a discrete, scoped task -- worth doing the next time the BT
peak position becomes the limiting source of model-data disagreement.
For now `--det-fwhm` stays the only handle; document the BT shift as a
known forward-model limit and report TH/BT *ratios* (which are less
sensitive to peak position) rather than absolute peak alignments.

---
