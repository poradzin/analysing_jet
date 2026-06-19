# Shared LOS library `src/neutron/common/los_common.py` (added 2026-06-18)

Motivation: a new neutron diagnostic — **KM9 (MPRu)** — is being added.
M. Nocente supplied `src/neutron/km9/KM9.los` (329 696 cells) using the
same Genesis 8-column cell format as `_KM3.los` (`KM9_LoS_readme.txt`).
**KM9 is *horizontal* and crosses the vessel twice** — the file's `x`
range is [-3.94, 3.85] m (toroidal dominant), `y` in [1.70, 2.27] m,
`z` in [-0.07, 1.34] m, versor `u≈0.99` — fundamentally different from
KM14's near-vertical pencil chord (`|x|≤0.098`, `w≈1.0`). The
"simplified box" (vertical chord, fixed `R∈[Rmin,Rmax]`, full Z extent,
`w_tor=0.4 m`) that frames the KM14 pipeline is physically meaningful
only for that chord; for KM9 it has no analogue. The *real-cell*
detector-rate path (`los_file_detector_rate`), however, is
geometry-agnostic and should be shared. Hence: factor the shared code,
keep two thin scripts (one per diagnostic).

## What moved into `los_common.py`

One module hosts every symbol used by more than one LOS script (the
directory `src/neutron/common/` is the unit of organisation; we will
add more modules there as scripts are built):

* **I/O.** `find_run_dir` (was `bt_zone_integrator.find_run_dir`),
  `read_time_grid`, `read_thermal_slice`, `read_los_file`.
* **Equilibrium.** `_boundary_from_moments`, `_rhot_scatter`,
  `_lcfs_contour`; the low-level `CdfEquilibrium` (was duplicated in
  `los_th_bt_ratio.py` — the canonical version is the one with the
  pinned-axis-+-LCFS `rhot_pinned`); the high-level `EqCDF` (wraps
  `CdfEquilibrium`) and `EqPPF`.
* **Numerics.** `thntx_on_grid`, `find_rho_bnd`, `_subgrid_bin`
  (overlap- and DVOL-density-weighted), `_running_median3`.
* **LOS rates.** `los_shell_fraction` (box-chord `f(rhot)` with the
  enclosed-shell `f=1` pin gated by `rmag`), `los_file_detector_rate`
  (real-cell `f_det(rhot)`, `Rate_det`, `rho_50`).

What stays in `km14/los_thermal_rate.py` (~840 lines, down from 1395):
the box-chord pipeline `build_rz_grid` / `integrate_rate` /
`toroidal_rate`, the A/B/C decomposition, the 2x3 diagnostic plot, the
CLI, channel constants, and KM14-specific defaults
(`R_MIN/MAX_DEFAULT=2.60/3.16`, `W_TOR_DEFAULT=0.4`).

## `chord_kind` kwarg on `los_file_detector_rate`

The only KM14-specific assumption in `los_file_detector_rate` is the
`rhot_crit` construction at the end of the subgrid block: it takes
the innermost flux surface that reaches the R-boundaries of the
structured rhot grid as the radius inside which the chord encloses
whole flux surfaces. That encodes "vertical pencil chord with axis
inside the R-band". For a horizontal chord like KM9 there is no
analogous "enclosed shell" — the chord doesn't fully sweep any flux
surface. So `los_file_detector_rate` now takes
`chord_kind: "vertical"|"horizontal"|"auto"`:

* `"vertical"` (**default**, preserves KM14 numbers byte-for-byte) —
  compute `rhot_crit`, apply enclosed-shell flattening, apply outside
  3-bin median.
* `"horizontal"` — skip all three; `rhot_crit` returned as `None`,
  `f_det` is the raw DVOL-weighted-subgrid binning, no cosmetic
  median3.
* `"auto"` — infer from `median(|cells["w"]|)` (vertical if > 0.9,
  else horizontal). KM14 cells (`w≈1.0`) → vertical; KM9 cells
  (`w≈0.05`) → horizontal. Same default produces the expected
  behaviour for each.

## How callers were threaded onto the lib

* `km14/los_thermal_rate.py`: `sys.path.insert(0, ../common)` at the
  top, then `from los_common import …`. All migrated function
  *bodies* removed; only box pipeline + plots + CLI remain. The script
  default `chord_kind` is unchanged (`"vertical"`).
* `km14/los_th_bt_ratio.py`: dropped its duplicated `CdfEquilibrium`,
  dropped the three lazy `from los_thermal_rate import …` clauses
  (cycle-break trick is no longer needed), pulls
  `CdfEquilibrium / EqCDF / read_los_file / los_shell_fraction /
  los_file_detector_rate / find_run_dir` from `los_common`. `bzi`
  import kept — it still uses `bzi.read_thermal_profiles /
  fast_ion_density / list_fbm_indices / read_fi_distribution /
  read_neut_rates / sample_fast_speeds / renormalize_to_bdens`.
* `km14/bt_zone_integrator.py`: its own `find_run_dir` body deleted,
  replaced with `from los_common import find_run_dir, DEFAULT_LOCAL_BASE,
  HEIMDALL_BASE`. The re-export keeps `bzi.find_run_dir(...)` working
  unchanged for **all** other callers — `bt_los_emissivity.py` and
  `km14_spectrum.py` still call `bzi.find_run_dir` and need no edit.

Files **intentionally not touched** (no breakage; clean up
opportunistically): `km14_spectrum.py` and `bt_los_emissivity.py`
still use `bzi.find_run_dir` (works via the re-export);
`bt_poloidal_distribution.py` still carries its own duplicate
`find_run_dir` — independent of the shared lib, also unchanged.

## Verification — M29 t=53.5268 s, idx 2 — every number bit-identical

`los_thermal_rate.py 104614 M29 -t 53.5268 --channel total --los-file
_KM3.los`:

```
Rate_LOS            = 5.5845e+15 n/s
Rate_tor            = 2.5608e+17 n/s
rho_bnd             = 0.3585
A / B               = 1.8821e+17 / 6.7862e+16  (73.50 / 26.50 %)
closure sum(THNTX*DVOL*f)  = 2.5808e+17 n/s
Rate_chord (real)   = 1.8573e+15 n/s
Rate_det            = 1.8247e+08 n/s
closure sum(THNTX*DVOL*f_det) = 1.8237e+08  (resid 5.26e-04)
rho_50              = 0.2520
```

`los_th_bt_ratio.py 104614 M29 --idx 2`:

```
(TH/BT)_LOS                       = 0.5294
cumulative (TH/BT)(rhot=1)        = 0.5299
TH chord = 4.4135e+15, BT chord = 8.3372e+15
TH whole-plasma = 3.8436e+17, BT whole-plasma = 8.8700e+17
```

All match the pre-refactor baseline captured immediately before the
migration (and reproduce the CONTEXT.md targets recorded for these
runs). `rho_bnd 0.3585` is the wide-chord (default `Rmin/Rmax
2.60/3.16`) value; the historical `0.3070` quoted earlier in CONTEXT
is the narrow-chord `2.70/3.10` value — both still reachable from CLI.

