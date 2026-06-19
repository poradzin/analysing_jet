# KM9 (MPRu) line-of-sight — context

KM9/MPRu is a *horizontal*, semi-tangential chord that crosses the vessel
twice (contrast KM14's near-vertical pencil chord). Its scripts live in
`src/neutron/km9/` and reuse the geometry-agnostic real-cell path from
`common/los_common.py`. See [shared_los_library.md](shared_los_library.md)
for the library that was factored out to support it.

## Stage B — tomorrow morning (planned 2026-06-19)

**Add `src/neutron/km9/los_thermal_rate.py`** as a thin sibling
script. Pure real-cell path; no box, no `rho_bnd`, no A/B/C
decomposition (the box rate has no clean analogue when the chord
wraps the torus tangentially and crosses the vessel twice).

Outline:

1. CLI: `pulse runid [--data-dir DIR] [--channel total|dd|dt]
   [-t TIME] [--los-file PATH] [--no-subgrid] [--plot] [--save]
   [--eq-source cdf|ppf …]`. Default `--los-file` points to
   `src/neutron/km9/KM9.los`.
2. Wiring: `sys.path.insert(0, ../common)`, then
   `from los_common import find_run_dir, read_time_grid,
   read_thermal_slice, read_los_file, EqCDF, EqPPF,
   los_file_detector_rate`. No box-pipeline imports
   (`los_shell_fraction`, `integrate_rate`, `toroidal_rate` are not
   needed).
3. Call `los_file_detector_rate(..., chord_kind="auto")` — the
   auto-classifier reads `median(|w|) ≈ 0.05` for KM9 cells and
   selects `"horizontal"`. Verify on day 1 by also running
   `chord_kind="horizontal"` explicitly and confirming identical
   results.
4. Report `Rate_chord = sum eps*V`, `Rate_det = sum eps*C`, closure
   `sum(THNTX*DVOL*f_det)`, `rho_50` (signal-median radius). No
   `rho_bnd` line. Cells-inside-LCFS count and `vol_tot/c_tot` for
   bookkeeping.
5. KM9-tailored plots (1x3): (a) poloidal scatter of cells coloured
   by `log10 C` over the LCFS — this is the diagnostic that visually
   demonstrates the twice-crossing horizontal geometry; (b)
   `f_det(rhot)` (single curve, no `rhot_crit` overlay since the
   construction does not apply; no enclosed plateau expected); (c)
   cumulative detector-signal fraction `cum_frac(rhot)` with `rho_50`
   marked.
6. `--save` writes `src/tmp/<pulse><runid>_KM9_LOS_profile_<chan>_<eq>_t<t>s.txt`,
   columns `rhot THNTX[n/m3/s] f_det THLOS_det[n/m3/s] DVOL[m3]`.
   Note the column header uses `THLOS_det` (the generic name); the
   common-lib dict key `thkm14_det_si` is kept for backward
   compatibility with km14 callers but interpret it as `THLOS_det` —
   refactor the key name to `thlos_det_si` in a follow-up that also
   updates `km14/los_thermal_rate.py` and `km14/km14_spectrum.py` (if
   it reads that key).

Expected first-cut results (no historical baseline — KM9 is new):

* ~329 696 cells, fraction inside LCFS will vary because the LOS
  extends well beyond the plasma both inboard and outboard
  (R∈[1.7, …] m at one end, dragged out by large |x| at the other).
* `Rate_det` should be of the same order as KM14's (1.8e+8 n/s) but
  not numerically comparable — different etendue normalisation
  altogether. The closure
  `sum(THNTX*DVOL*f_det) ≈ Rate_det` is the cross-check that the
  pipeline is consistent for the new chord; that's the acceptance
  test.
* `rho_50` will reflect where on the chord the LOS picks up most of
  its thermal signal. For a horizontal chord at `Z` near the
  midplane that should be midradius-ish; for an upper chord it will
  be larger. Worth comparing to the KM14 value as a sanity reading
  of the geometry.

Open question to flag for tomorrow: the KM9 LOS reaches well **outside
the LCFS** (and even outside the equilibrium computational box for
some cells). `los_file_detector_rate` zeros emissivity outside the
LCFS via the polygon mask, so those cells contribute zero to
`Rate_det` (correct), but `rhot_at_points` in `_rhot_scatter` uses
griddata-cubic which extrapolates over the `psi_n` node set
(PSIRZ grid + axis + LCFS) and may return NaN or unstable values for
the far-out cells before the inside mask zeros them. Likely benign
(the mask kills it) but worth a print of the fraction of NaN/clipped
rhot values on the first run.

## Stage B — DONE (2026-06-19)

`src/neutron/km9/los_thermal_rate.py` written as planned: a thin sibling to
the KM14 script, pure real-cell path (no box, no `rho_bnd`, no A/B/C). Imports
only `find_run_dir, read_time_grid, read_thermal_slice, read_los_file, EqCDF,
EqPPF, los_file_detector_rate` from `los_common`; default `--los-file` =
`src/neutron/km9/KM9.los`. KM9-tailored 1x3 plot ((a) cell scatter coloured by
`log10 C` over LCFS, (b) `f_det(rhot)`, (c) cumulative detector-signal fraction
with `rho_50`). `--save` writes `src/tmp/<run>_KM9_LOS_profile_<chan>_<eq>
_t<t>s.txt` with columns `rhot THNTX f_det THLOS_det DVOL` (reads the lib's
`thkm14_det_si` key, labels the column `THLOS_det`).

**Verified — 104614 M29, t=53.5268 s, CDF eq, total channel:**
```
329696 cells; median |w| = 0.083  -> chord_kind auto -> horizontal (rhot_crit None)
rhot NaN cells: 0 / 329696 (0.00 %)   <- the flagged open question: benign, none
Cells inside LCFS: 229687 / 329696 (69.7 %)
Rate_chord = 1.3465e+16 n/s,  Rate_det = 4.3415e+09 n/s
closure sum(THNTX*DVOL*f_det) = 4.3403e+09  (resid 2.75e-04)   <- acceptance test PASS
rho_50 = 0.2642   (cf. KM14 box rho_bnd 0.3585 / rho_50 0.2520)
```
`--chord-kind auto` and `--chord-kind horizontal` give byte-identical output
(day-1 verification, step 3 done). `--save` and `--channel dd` (rho_50=0.2717)
also exercised; `--plot` renders.

**Geometry note (from the cell scatter):** KM9's footprint here is `Z∈[-0.07,
1.337]`, `R∈[1.82, 4.30]` — an **upper, near-horizontal** oblique band, not a
midplane chord; it crosses the torus twice *toroidally* (|x|≤3.94 m), which is
why R reaches 4.30 m (outboard) and 1.82 m (inboard) while Z stays in the upper
half. `f_det` is peaked near-axis with a ~1-2% finite-cell wiggle on the
roll-off (no enclosed-shell plateau, no cosmetic median3 — both correctly
skipped for `horizontal`). The open-question NaN check came back clean (0
cells), so the far-out griddata extrapolation is benign as suspected.

## `-t` / `--idx` time-window averaging (added 2026-06-19)

The thermal-rate scripts (`km9/` and `km14/los_thermal_rate.py`) now support
**time-window averaging** of the TRANSP signals, driven by a shared helper
`los_common.resolve_time_selection(cdf_path, run_dir, run_id, times, idx)`:

* **`-t t1`** (one value) -- snap to the nearest TIME3 output slice (old
  single-value behaviour; KM14 baseline reproduced byte-for-byte).
* **`-t t1 t2`** (two values) -- average THNTX/X/DVOL over **all TIME3 slices in
  `[t1, t2]`** (straight unweighted mean); equilibrium = single slice nearest
  the **window midpoint** `(t1+t2)/2` (user decision 2026-06-19).
* **`--idx X`** (1-based, mutually exclusive with `-t`) -- read
  `OUTTIM`/`AVGTIM` from the `&ACFILE` group of the run namelist
  `<run_dir>/<run_id>TR.DAT` (via `f90nml`), average over the fast-ion output
  window `[OUTTIM(X)-AVGTIM, OUTTIM(X)]`, equilibrium nearest **OUTTIM(X)** (the
  window upper bound = nominal output time). **No `_fi_X.cdf` read needed.**
* none -- run-midpoint slice.

Mechanics: `read_thermal_slice` now accepts an int **or a sequence** of TIME3
indices and means THNTX/X/DVOL element-wise across them. `resolve_time_selection`
returns `(trinds, eq_ref_jet, rep_jet, desc)`; `main()` builds the equilibrium at
`eq_ref_jet` (a *single* slice -- the equilibrium is **not** averaged, per the
"map the averaged profile onto one equilibrium" decision) and tags
filenames/titles with `rep_jet`. A window narrower than the TIME3 spacing
(~0.02 s here) falls back to the nearest single slice with a warning. `-t`
arity (1 or 2) and `-t`/`--idx` mutual exclusion are validated in `parse_args`.

**Averaging semantics (clarified from the TRANSP ACFILE doc, 2026-06-19).** The
`[OUTTIM(j)-AVGTIM, OUTTIM(j)]` window and `MTHDAVG` govern only the variables
in the run's `SELAVG` list -- here `FBM BMVOL BDENS2 EBA2PL EBA2PP`, the
**fast-ion** data written to `_fi`/`_neut` (and hence the BT neutrons).
**THNTX (thermal) is NOT in `SELAVG`**, so TRANSP never applied its own
ACFILE/`MTHDAVG` averaging to it -- THNTX is a plain TIME3 output. Our windowed
mean of THNTX is therefore *our own construct*, done to make the thermal profile
time-consistent with the window the fast-ion/BT snapshot was averaged over (the
motivation for `--idx` in the first place: matching TH to the `_fi_X`/`_neut_X`
data). `MTHDAVG=2` ("sample on each heating-source timestep completion")
described the fast-ion sampling, not ours; it is read and printed (tagged
"fast-ion MTHDAVG") for provenance only. A straight TIME3 mean is the best
achievable for THNTX (not stored on the heating-source grid) and equals the
trapezoidal time-average on the uniform 0.02 s TIME3 grid. `.DATA<j>` /
`_fi_<j>` are written at `OUTTIM(j)`, confirming `--idx X <-> OUTTIM(X)`;
windows are non-overlapping (`OUTTIM` spacing >> `AVGTIM`).

**Verified 104614 M29:** `--idx 2` -> OUTTIM=13.5s, AVGTIM=0.35s, TRANSP window
[13.15,13.5]s = 17 TIME3 slices [48..64]; eq @ JET 53.500->53.509s. KM14 `--idx 2`:
rho_bnd 0.3584 / Rate_det 1.8233e8 / rho_50 0.2525 (vs single-slice 0.3585 /
1.8247e8 / 0.2520 -- negligible drift, plasma quasi-stationary over the window).
KM9 `--idx 2` == `-t 53.15 53.5` for the *rate* (same slices -> Rate_det
4.3361e9, closure resid 2.7e-4) but different eq time (53.509 vs midpoint
53.329); rho_50 0.2639 either way. Single `-t 53.5268` reproduces the
pre-change KM14 baseline exactly.

## KM9 near-axis `f_det` dip -> zero un-crossed inner shells (added 2026-06-19)

User saw a dip in the KM9 `f_det` detector-coupling weight at the very core
(104614 M29 `--idx 2`). **Diagnosed:** it is *not* a mapping-at-origin artifact.
The KM9 chord's closest approach to the magnetic axis is **0.71 cm** -> rhot
**0.0082**; no in-LCFS cell reaches below that. The psin->rhot table does floor
at rhot=0.005 at psin=0 (XB[0]=0.005), but that floor does **not** bind here
(geometry puts the chord at 0.0082, above it), so refining the origin mapping
would not move it. The visible dip was the two innermost TRANSP bins
(rhot < 0.01, below the chord's reach) being filled *only* by the DVOL-weighted
subgrid splat's symmetric tails from cells at rhot~0.008-0.02, depositing C
proportional to shell DVOL -> a flat shelf at exactly half the plateau
(`f_det[0]==f_det[1]`). Worth 0.055 % of `Rate_det` (bins 0-1).

**Physics:** KM9 is a *horizontal grazing* chord -- unlike vertical KM14 it does
**not** enclose whole core shells, so a roll-down toward the un-crossed core is
expected; only the flat half-shelf in the un-reached bins was an artifact.

**Fix (user chose the physically-honest option).** In `los_file_detector_rate`,
for **non-vertical** chords only (`chord_kind != "vertical"`), compute
`rhot_min` = min rhot over in-LCFS cells (the chord's closest approach) and zero
`f_det`/`c_bin`/`v_bin` on flux shells whose **upper edge <= rhot_min** (never
crossed -> no emission reaches the detector), folding their leaked C/V into the
**first crossed bin** so `Rate_det` closure (total C) is preserved. The bin that
*straddles* `rhot_min` keeps its reduced value, giving an honest roll-up from the
closest-approach radius. Returned/plotted as `rhot_min` (cyan marker on the KM9
`f_det` panel). **Vertical chords (KM14) are untouched** -- they geometrically
enclose the inner shells even when no discrete cell centre lands there, and keep
their `rhot_crit` enclosed-shell flattening (the plateau); `rhot_min` is `None`
for them. Gated on `chord_kind`, not `rhot_crit`, so it never collides with the
KM14 plateau.

**Verified 104614 M29 idx 2:** KM9 `f_det[:6]` = `[0, 1.53e-8, 1.75e-8, 2.12e-8,
2.26e-8, 2.16e-8]` (was `[1.15e-8, 1.15e-8, 1.75e-8, ...]`); `rhot_min=0.0082`;
`Rate_det=4.3361e9` and closure resid `2.72e-4` **unchanged**; `rho_50=0.2639`
unchanged. KM14 `--idx 2` and the single-slice baseline (`-t 53.5268`:
rho_bnd 0.3585 / Rate_det 1.8247e8 / rho_50 0.2520) reproduce **byte-for-byte**
(rhot_min None, flat plateau 6.87e-10, rhot_crit 0.1242 intact).

## KM9 outer-core `f_det` oscillation -> median5 de-speckle (added 2026-06-19)

User saw a ~15% peak-to-peak oscillation in KM9 `f_det` in the outer core.
**Source (diagnosed, not C/etendue structure and not shot noise):** geometric
**aliasing of the regular cell lattice against the curved flux shells**. The
KM9 cells sit on a regular Cartesian grid (71 uniform z-slices dZ=0.0201 m x
regular along-chord steps; only 34.5k unique rhot among 230k in-LCFS cells), so
the per-bin *cell count* itself oscillates ~22% (~770-1440 cells/bin, far above
the ~3% Poisson floor) -- a systematic moiré, matched by the C-weighted
oscillation. The DVOL-weighted subgrid splat already suppresses it ~5-8x (point
binning is 90-250% p2p) to ~2% std in the core, but ~15% p2p by rhot 0.5-0.7 and
worse at the edge. Confirmed aliasing (not under-smearing): scaling the splat
half-span is **non-monotonic** (x1.5 worse than x1 or x2) -- a beat between the
kernel width and the lattice period that a fixed-width linear splat cannot
cancel. KM14 hides the identical artifact with a 3-bin running median applied
*outside* `rhot_crit`, but that lives inside the vertical-only `rhot_crit` block,
so the horizontal KM9 path got **no** smoothing.

**Fix (user requested median5).** Added `_running_median(a, k=5)` (centered,
symmetric-shrinking ends; `_running_median3` left intact for the vertical path
so KM14 stays byte-for-byte). In `los_file_detector_rate`, a new
`elif subgrid and median and chord_kind != "vertical":` branch applies a 5-bin
running median to `f_det` on the **crossed** shells only (`edges > rhot_min`),
leaving the zeroed un-crossed core untouched. Cosmetic: `Rate_det` (= sum eps*C)
and `rho_50` are from the cells directly and don't change.

**`median` toggle / `--no-median` flag (added 2026-06-19).**
`los_file_detector_rate(..., median=True)` gates *both* cosmetic medians (the
horizontal 5-bin and the vertical outside-`rhot_crit` 3-bin). Both `km9/` and
`km14/los_thermal_rate.py` expose `--no-median` to return the raw binned
`f_det` (to estimate the median's effect). Default True preserves all baselines.
Effect on KM9 idx 2: mean |median-raw|/raw ~3.7% (peaks ~64% on a single noisy
edge bin where the raw spike is clamped); `Rate_det`/`rho_50` bit-identical
with or without (`--no-subgrid` is the *other* toggle, and unlike `--no-median`
it disables the anti-alias splat entirely -> the very noisy point-binned f_det).

**Verified 104614 M29 idx 2:** KM9 outer-core osc (std) 0.3-0.5 / 0.5-0.7 /
0.7-0.9 = 2.0/4.2/9.6% -> **1.1/1.3/1.6%**; `f_det[:6]` =
`[0, 1.53e-8, 1.75e-8, 2.12e-8, 2.16e-8, 2.16e-8]` (core zero + rhot_min 0.0082
preserved); `Rate_det=4.3361e9`, `rho_50=0.2639` unchanged, closure 2.72e-4 ->
8.3e-4 (still negligible). KM14 (`vertical`) fully unchanged: flat plateau
6.87e-10, single-slice baseline rho_bnd 0.3585 / Rate_det 1.8247e8 / rho_50
0.2520 reproduced byte-for-byte.

## KM9 plot: added top-of-machine (x-y) view (added 2026-06-19)

User wanted a top view to see whether the LOS crosses the magnetic axis. KM9
`diagnostic_plots` went **1x3 -> 2x2** (existing three panels kept): new panel
(0,1) scatters the cells in the toroidal **x-y plane** (coloured by `log10 C`,
out-of-LCFS dimmed) with the **R=Rmag magnetic-axis circle** (red) and the LCFS
annulus `R in [Rb.min, Rb.max]` (blue dashed) + machine axis at the origin.
`diagnostic_plots` now also takes `cells` (for `x, y`; `losf` only carries
`R=hypot(x,y)`). KM9 cells are a *thin* band `y in [1.70, 2.27]` m spanning
`x in [-3.94, 3.85]` m, so in x-y the chord is a near-horizontal secant; with
equal aspect it visibly crosses the R=Rmag circle at `x ~ +/-sqrt(Rmag^2-y^2)`
~ +/-2.29 m (two toroidal points), consistent with the poloidal `rhot_min=0.0082`
grazing. Panels (1,0) f_det and (1,1) cumulative unchanged.

**axis-in-cloud vs rhot_min, clarified (2026-06-19).** User noticed the axis
red-cross sits *inside* the poloidal cell scatter yet `rhot_min=0.0082` (>0).
Reconciliation: the poloidal (R,Z) panel is a **projection** of the 3-D
toroidally-extended LOS -- 394 cells project within 3 cm of the axis from 76
distinct toroidal x positions (`x in [-2.5, 2.4]`), surrounding it in all four
quadrants, so it *looks* enclosed; but `rhot_min` is the nearest discrete cell
**center** (0.71 cm, sub-grid vs the ~2 cm z-slice spacing). So `rhot_min` is a
**sampling-resolution floor**, not a hard geometric miss -- the continuous LOS
likely grazes/contains the axis, we just can't resolve below the cell spacing.
The f_det plot label was softened "closest approach" -> "min sampled rhot
(resolution floor)" accordingly (the code comment already said "minimum
*sampled* rhot").

**2x3: added per-shell emission panel (2026-06-19).** New panel (0,2) (like
KM14's THNTX/THKM14 profile) plots `THNTX*DVOL` (whole-plasma per-shell rate
[n/s], left axis) and `THNTX*f_det*DVOL` (detector-weighted per-shell, right
axis; integral = Rate_det) vs rhot, with `rho_50` marked. Twin y-axis because
`f_det ~ 1e-8` makes the two differ by ~1e8 (unlike KM14's box `f<=1`, single
axis). The detector-weighting pulls the per-shell peak inward, 0.292 -> 0.273
(104614 M29 idx 2); detector curve sums to 4.34e9 = Rate_det (closure). The
unused (1,2) slot is removed via `fig.delaxes`.

## KM9 `los_th_bt_ratio.py` (real-cell TH/BT, added 2026-06-19)

`src/neutron/km9/los_th_bt_ratio.py` -- the KM9 sibling of
`km14/los_th_bt_ratio.py`, built per the Stage-B follow-up. KM9 is horizontal so
the KM14 **box-chord** layer (Rmin/Rmax/wtor rates, `chord_integral`,
`los_shell_fraction`, and the **point-detector** `--bt-aniso` approximation) has
no analogue and is dropped; only the geometry-agnostic **real-cell** layer
remains (mirrors how `km9/los_thermal_rate.py` relates to its KM14 counterpart).

Pipeline: read TH (`THNTX` via `read_thermal_slice`) + BT (flux-averaged BTN via
`flux_avg_profile`/`interp_flux_to`) at `--idx` (the `_fi`/`_neut` index =
`OUTTIM(idx)`); `EqCDF(cdf, fi.time+40)`; `los_file_detector_rate` on the real
`KM9.los` cells (`chord_kind="auto"`->horizontal) **twice** (TH and BT) -- f_det
is purely geometric so both calls return identical f_det (asserted); the two
`rate_det` give `(TH/BT)_LOS = Rate_det(TH)/Rate_det(BT)` directly from the
cells. **BT treated isotropic in v1** (same f_det for TH and BT; directional
`--bt-aniso` deferred -- heavier for a horizontal chord). Reports local
TH/BT(rhot), C-weighted cumulative (TH/BT)_LOS, whole-plasma 0D TH/BT, core
enhancement, and TH/BT signal-median rho_50. 2x3 plot mirrors the KM9 thermal
script (poloidal cells, top view, f_det, TH/BT-vs-rhot, per-shell
TH*f_det*DVOL & BT*f_det*DVOL).

**Helper migration:** `read_scalar_totals`, `flux_avg_profile`, and
`interp_flux_to` (was `_interp_flux_to`) moved from `km14/los_th_bt_ratio.py`
into `los_common.py` (avoids a km9->km14 dependency for these; both scripts
import them). `bt_zone_integrator` (the `_fi`/`_neut` reader) still lives in
km14/, so km9 adds `../km14` to sys.path to import `bzi` -- a real shared-infra
dependency; moving `bzi` to common/ is a later refactor (many callers).

**Verified 104614 M29 idx 2 (total):** (TH/BT)_LOS = **0.5095** (BT/TH 1.963)
vs full-plasma 0D **0.4333** (core enhancement 1.18x; TH rho_50 0.264 < BT 0.321
-> TH more core-concentrated, as expected). Whole-plasma TH 3.8436e17 / BT
8.8700e17 / BTNTS 8.8287e17 match the KM14 TH/BT script exactly. `--channel dd`:
(TH/BT)_LOS 0.3601 vs plasma 0.3239. The cells-direct `thbt_los` is robust to
`--no-median` (endpoint check 0.5094 raw vs 0.5077 with median -- the ~0.4% gap
was just the cosmetic f_det smoothing; the headline ratio is from `rate_det`,
unaffected). **KM14 `los_th_bt_ratio.py` reproduces byte-for-byte** after the
migration (0.5294 / 0.5299).

## `km9/plot_LoS.py` -- LOS geometry views (added 2026-06-19)

Motivated by Fig. 9 (top) of Andersson Sunden et al., NIM A 610 (2009) 682
(`~/jet/documentation/2009_Anderrson_Sunden_MPRu_at_JET_...pdf`, PDF p.10) -- a
*side elevation* of the torus showing the semi-tangential MPRu chord skimming
the plasma centre. New shared module `src/neutron/km9/plot_LoS.py` with
`plot_los_geometry(cells, inside_cells, Rb, Zb, Rmag, Zmag, ...)` draws **three
orthogonal projections** of the real LOS cell cloud (coloured by `log10 C`,
out-of-LCFS dimmed):
1. **Poloidal (R, Z)** -- cells + LCFS + axis; closest-approach cell starred.
2. **Top (x, y)** -- cells + the `R=Rmag` axis circle + LCFS annulus (the
   tangential secant cuts the circle at two toroidal points).
3. **Side elevation (x, Z)** -- the Fig.9 analog: cells as a near-horizontal
   band, with the magnetic-axis *height* `Z=Zmag` (a thin dotted **reference**,
   not a line the chord follows), plasma Z-extent band, and `|x|=R_out`. All
   panels equal-aspect.

**Axis-passage highlight (added 2026-06-19, refined).** The magnetic axis is a
*ring* (`R=Rmag` at every phi), so by axisymmetry its 3-D distance to a cell is
the poloidal `hypot(R-Rmag, Z-Zmag)`. The cells within `axis_tol` (default
5 cm) of the ring are **recoloured red** in all three panels -- this is the
"change of colour" showing the axis *threads through* the LOS (in one toroidal
side, out the other), rather than the chord looking like it runs *along* the
axis. The `Z=Zmag` line in the side view is a thin dotted **reference** only.
On 104614 M29 idx 2 the red cells form **two symmetric clusters** at `x~-2.5`
and `x~+2.5` m (`y~2`, `Z~0.29`, both ~0.7 cm min distance) -- the axis enters
the chord at one toroidal lobe and leaves at the other, exactly the grazing
geometry of Fig. 9. (An earlier lime/cyan two-star version was replaced by this
colour highlight per user preference.)

Called optionally via **`--plot-los`** (lazy `import plot_LoS`) from both
`km9/los_thermal_rate.py` and `km9/los_th_bt_ratio.py`. The function is
geometry-agnostic (takes any cell cloud), so KM14 could call it too if wanted.

**Geometry panels removed from `--plot` (2026-06-19).** The `--plot` analysis
figures used to embed their own poloidal+top panels, so `--plot --plot-los`
showed the geometry twice. Both analysis figures are now **1x3** with no
geometry: `los_thermal_rate` -> [f_det, per-shell emission (twin-axis
THNTX*DVOL & THNTX*f_det*DVOL), cumulative signal fraction];
`los_th_bt_ratio` -> [f_det, TH/BT vs rhot, per-shell TH/BT*f_det*DVOL]. Geometry
(poloidal/top/side) lives only in `plot_LoS` / `--plot-los`. `diagnostic_plots`
dropped its geometry args (`cells`/`Rb`/`Zb`/`Rmag`/`Zmag`). `--plot` alone now
shows one figure; `--plot --plot-los` two distinct ones. `Rate_det`/
`(TH/BT)_LOS` unchanged (plot-only).

## Stage B follow-ups (not Stage B itself; later)

* When `km14_spectrum.py` is opened next, swap its
  `bzi.find_run_dir` to `from los_common import find_run_dir` for
  consistency; same for `bt_los_emissivity.py` and
  `bt_poloidal_distribution.py` (which carries its own duplicate).
* Consider renaming the dict key `thkm14_det_si ->  thlos_det_si` in
  `los_file_detector_rate` (currently kept for back-compat); requires
  touching `km14/los_thermal_rate.py` and any other reader.
* ~~Mirror the split for `los_th_bt_ratio.py`~~ **DONE 2026-06-19** --
  `km9/los_th_bt_ratio.py` built (real-cell TH/BT, isotropic BT v1); see
  the dedicated section above. Remaining: directional `--bt-aniso` for the
  KM9 horizontal chord (deferred), and moving `bt_zone_integrator` into
  common/ (km9 currently sys.paths ../km14 to import it).
