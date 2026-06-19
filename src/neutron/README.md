# Neutron diagnostics — line-of-sight analysis

Tools for line-of-sight (LOS) analysis of JET neutron diagnostics from
TRANSP/NUBEAM output. The core deliverables are the **KM14** (vertical
diamond spectrometer) and **KM9 / MPRu** (horizontal magnetic proton
recoil) line-integrated thermal-neutron rates and thermal/beam-target
(TH/BT) ratios, plus a KM14 diamond neutron-spectrum forward model.

For the physics and numerics rationale behind each script — algorithms,
bug investigations, validation results — see [`doc/`](doc/) (index at the
bottom of this file). This README is the **how to run** reference.

## Directory layout

```
src/neutron/
  common/los_common.py     shared LOS library (equilibrium, I/O, rates)
  km14/                    KM14 vertical-chord scripts (+ BT physics chain)
  km9/                     KM9/MPRu horizontal-chord scripts
  doc/                     physics & numerics deliberations
```

`*.los` files (`km14/KM3.los`, `km9/KM9.los`) are M. Nocente's exact LOS
cell clouds (8-column Genesis format `x y z C V u v w`; see the adjacent
`*_LoS_readme.txt`). `C` is the per-cell etendue weight, `(u,v,w)` the
emission versor toward the detector.

## Environment

* **Self-contained on the TRANSP CDFs** (`netCDF4 + numpy + scipy +
  matplotlib` only) by default — no `ppf`, no `profiles.Eq` — so the
  scripts run in the WSL dev env, on freia and on heimdall alike.
  `los_thermal_rate.py` has an opt-in `--eq-source ppf` path that lazily
  imports `profiles`/`change_rho`/`ppf` (Heimdall only).
* **numpy ≥ 2.0** is required (`np.trapezoid`). On legacy numpy swap to
  `from scipy.integrate import trapezoid`.
* Run scripts **from inside their own directory** (they import each other
  as plain modules, e.g. `import bt_kinematics`):
  ```bash
  cd ~/jet/analysing_jet/src/neutron/km14   # or .../km9
  ```

## Data directory search order

Per `(pulse, run)` the scripts need three CDFs in the run directory:
`<run>.CDF` (main), `<run>_fi_<idx>.cdf` (fast-ion dist + zone geometry),
`<run>_neut_<idx>.cdf` (per-zone BT rates). Search order:

1. `--data-dir <base>/<pulse>/<run>`
2. `~/jet/data/<pulse>/<run>` (local WSL)
3. `/common/transp_shared/Data/result/JET/<pulse>/<run>` (heimdall)

Locally **104614 M29** and **M30** are present under `~/jet/data/`.
`--idx` selects the FBM time window (1-based); `-t` selects/averages
TIME3 slices (see "time-window averaging" below).

## Scripts

### KM14 (`km14/`)
| script | purpose |
|---|---|
| `bt_kinematics.py` | Pure-physics core: Bosch–Hale σ(E_cm), CM-frame sampler, classical CM→lab kinematics. Self-test: `python bt_kinematics.py`. |
| `bt_zone_integrator.py` | Per-zone 4π beam-target reactivity; acceptance-tested against NUBEAM `BTN4`/`BTN1` (≈1.00). |
| `bt_los_emissivity.py` | Direction-resolved BT emissivity `dε/dΩ(n̂)`; KM14 anisotropy factor g≈1.00. Self-test: `--test`. |
| `bt_poloidal_distribution.py` | Spatial poloidal-angle map of BT emissivity (step 1; standalone CDF reader). |
| `los_thermal_rate.py` | KM14 thermal LOS rate, equivalent `rho_bnd`, `THKM14(rhot)` weight, real-cell `f_det` (`--los-file`). |
| `los_th_bt_ratio.py` | KM14 LOS TH/BT ratio (`--channel total\|dd\|dt`, `--bt-mode flux\|zone`, `--bt-aniso`). |
| `km14_spectrum.py` | KM14 diamond neutron-spectrum forward model (TH + B-th DT), overlaid on the measured spectrum; χ²/fit. |
| `extract_km14_png.py` | Extract curves + data points/error bars from the published Nocente PNG. |
| `estimate_outer_th_fraction.py` | Early geometric outer-fraction estimate — **superseded** by `los_thermal_rate.py`; kept for reference. |

### KM9 / MPRu (`km9/`)
| script | purpose |
|---|---|
| `los_thermal_rate.py` | KM9 thermal LOS rate (real-cell path only; no box / `rho_bnd` / A/B/C). |
| `los_th_bt_ratio.py` | KM9 LOS TH/BT ratio (real-cell, isotropic BT v1). |
| `plot_LoS.py` | KM9 LOS geometry: three orthogonal projections (poloidal, top, side elevation). Called via `--plot-los`. |

## Quickstart (recommended order — each step validates the next)

```bash
cd ~/jet/analysing_jet/src/neutron/km14

# 0. Kinematics unit tests (no data needed)
python bt_kinematics.py                       # σ spot-checks + kinematic bands → PASS

# 1. Per-zone reactivity vs NUBEAM (acceptance test ≈ 1.00)
python bt_zone_integrator.py 104614 M30 --idx 1 --plot

# 2. LOS angular emissivity — KM14 sees BT as isotropic (g ≈ 1.00)
python bt_los_emissivity.py --test
python bt_los_emissivity.py 104614 M30 --idx 1 --plot

# 3. The deliverables
python los_thermal_rate.py 104614 M29 -t 53.5268 --channel total --los-file _KM3.los
python los_th_bt_ratio.py 104614 M29 --idx 2 --plot
python km14_spectrum.py   104614 M29 --idx 1 --include-scatt

# KM9
cd ../km9
python los_thermal_rate.py 104614 M29 --idx 2 --plot --plot-los
python los_th_bt_ratio.py  104614 M29 --idx 2 --plot
```

Every `--plot` opens an interactive window; `--save <file.png>` writes a
PNG headless (handy over SSH); `--no-plot` forces text-only. The older
scripts also accept `key=value` argument style (e.g. `dda=eftp`).

### Time selection (`los_thermal_rate.py`, both KM14 and KM9)
* `-t t1` — snap to the nearest TIME3 output slice.
* `-t t1 t2` — average THNTX/X/DVOL over all TIME3 slices in `[t1,t2]`;
  equilibrium = single slice nearest the window midpoint.
* `--idx X` — average over the fast-ion output window
  `[OUTTIM(X)-AVGTIM, OUTTIM(X)]` read from the run namelist; equilibrium
  nearest `OUTTIM(X)`. (Mutually exclusive with `-t`.)

See [`doc/km9.md`](doc/km9.md) for the averaging semantics (THNTX is not
in `SELAVG`, so the windowed THNTX mean is our own time-consistency
construct).

## Documentation index (`doc/`)

* [thermal_los_rate.md](doc/thermal_los_rate.md) — KM14 thermal LOS rate:
  algorithm, equilibrium handling, near-axis & LCFS-crater fixes, subgrid
  anti-aliased binning, real-cell detector path (`f_det`, `rho_50`).
* [beam_target.md](doc/beam_target.md) — BT physics chain: `bt_kinematics`,
  `bt_zone_integrator` (F_D_NBI-vs-BDENS_D normalization gap),
  `bt_los_emissivity` (KM14 isotropy), `bt_poloidal_distribution`, NUBEAM
  MC-grid reference.
* [th_bt_ratio.md](doc/th_bt_ratio.md) — KM14 LOS TH/BT ratio: method,
  cumulative-ratio plateau, `--bt-aniso`, effective-weight figure.
* [km14_spectrum.md](doc/km14_spectrum.md) — diamond spectrum forward
  model, tritium-fraction story, PNG extraction, BT-peak-position issue,
  `--cone-deg` discussion.
* [shared_los_library.md](doc/shared_los_library.md) — the `los_common.py`
  refactor: what moved, the `chord_kind` switch, caller wiring.
* [km9.md](doc/km9.md) — KM9/MPRu: horizontal-chord handling, near-axis
  `f_det` fixes, aliasing/median de-speckle, time-window averaging,
  `plot_LoS` geometry views.
