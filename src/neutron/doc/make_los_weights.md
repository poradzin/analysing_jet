# `make_los_weights.py` — time-resolved LoS weight `.npz` for KM9 / KM14 (added 2026-06-25)

`src/neutron/make_los_weights.py` is a thin driver that walks a TRANSP run's
TIME3 grid and dumps the per-flux-shell **line-of-sight detector coupling
weight** `f_det(rhot)` for one of the JET neutron diagnostics (KM9 / MPRu or
KM14 diamond) into a single `.npz` file. No new physics, no new plots — it
calls the validated `los_common.los_file_detector_rate` per slice and stacks
the result.

## Why this script exists

The consumer is `~/jet/jet_tritium`, which currently approximates the
LoS-integrated thermal rate with two crude proxies:

* **KM14**: volume integral from `rhot = 0` to `0.2` (the "core ~25 % proxy").
* **KM9**: full-plasma volume integral (no LoS reweighting at all).

Both are stand-ins for the same thing: a LoS-weighted volume integral

    Rate_det(t) = Σᵢ epsᵢ(rhot, t) · f_detᵢ(t) · DVOLᵢ(t)              [n/s]

where `eps` is any volumetric emissivity sampled on the TRANSP rhot grid
`X(t)` and `i` runs over TRANSP flux shells. `f_det = C_bin / DVOL` is the
dimensionless **per-shell detector coupling weight** the `los_common`
machinery already computes — see [`thermal_los_rate.md`](thermal_los_rate.md)
for the construction, [`km9.md`](km9.md) for the horizontal-chord variant.
Closure is exact by construction:

    Σ(THNTX · f_det · DVOL)  ≡  Σ_cells eps · C_cell  =  Rate_det.

`make_los_weights.py` precomputes `f_det` (+ DVOL + reference THNTX) on the
full TIME3 grid so `jet_tritium` can do **one `np.load` + one `np.sum`** per
time slice instead of re-running TRANSP-LOS analysis. The weight is purely
geometry × equilibrium (channel-independent) so the same file covers TH and
the ~isotropic BT (KM14 anisotropy g ≈ 1.00; KM9 BT treated isotropic in v1
per CONTEXT).

## Reuses the validated path verbatim

Per TIME3 index `ti` the loop does exactly what the standalone
`km9/los_thermal_rate.py` / `km14/los_thermal_rate.py` do for a single slice:

1. `read_thermal_slice(cdf, th_var, ti)` — `THNTX(X)`, `X`, `DVOL` at `ti`
   (no time averaging; one slice at a time, so the resulting `f_det` lives on
   the *native* TIME3 cadence).
2. `EqCDF(cdf, t_jet)` (default) or `EqPPF(...)` — `psi(R,Z) → rhot` map
   from the same slice; cubic griddata with the magnetic axis pinned at
   `psi_n = 0` and the LCFS at `psi_n = 1` (the standard `los_common`
   pinning, see [`thermal_los_rate.md`](thermal_los_rate.md)).
3. `los_file_detector_rate(cells, eqs, ...)` — DVOL-Jacobian subgrid binning
   of `C`/`V` onto the TRANSP shells, vertical/horizontal chord-kind logic
   intact (enclosed-shell flattening for KM14, un-crossed-shell zeroing for
   KM9), cosmetic median-3/5 on by default.

So the per-slice numbers are **byte-for-byte** what the single-slice scripts
report. Spot check on 104614 M29 around t = 53.5 s:

| diag | source | `Rate_det` [n/s] | `rho_50` |
|---|---|---|---|
| KM9  | `make_los_weights.py` | 4.34 × 10⁹ | 0.264 |
| KM9  | `km9/los_thermal_rate.py` | 4.34 × 10⁹ | 0.264 |
| KM14 | `make_los_weights.py` | 1.82 × 10⁸ | 0.252 |
| KM14 | CONTEXT headline (real-cell) | — | 0.252 |

Closure residual (max over slices) of the saved arrays:
`max |Σ THNTX·f_det·DVOL − Rate_det| / Rate_det ≈ 1e-3` — the same level the
single-slice scripts report (the small mismatch comes from the cosmetic
running-median on `f_det`, which `Rate_det` itself does not use; the median is
purely a `f_det`-display de-speckle, the cell-summed rate is unchanged).

## Output schema (`np.savez_compressed`)

2-D arrays indexed `[it, irho]` — `it` over selected TIME3 slices, `irho`
over TRANSP rhot bins (typically `n_rho = 200`):

| key | shape | meaning |
|---|---|---|
| `rhot`         | (n_t, n_rho) | TRANSP rhot grid `X` per slice (sorted ascending) |
| `dvol_m3`      | (n_t, n_rho) | Flux-shell volume `DVOL` [m³] |
| `f_det`        | (n_t, n_rho) | **The weight**: `C_bin / DVOL` [-] |
| `thntx_si`     | (n_t, n_rho) | Reference `THNTX` [n/m³/s] for the selected channel |
| `thlos_det_si` | (n_t, n_rho) | `THNTX · f_det` [n/m³/s] (Σ·DVOL = `rate_det`) |

1-D per slice:

| key | shape | meaning |
|---|---|---|
| `times_jet`    | (n_t,) | JET time [s] (= TIME3 + 40 s) |
| `trinds`       | (n_t,) | original TIME3 index in the main CDF |
| `rate_det`     | (n_t,) | LoS-weighted thermal rate [n/s] (canonical) |
| `rate_chord`   | (n_t,) | bare real-chord emission `Σ eps·V` [n/s] (no solid angle) |
| `rate_closure` | (n_t,) | `Σ THNTX·f_det·DVOL` [n/s] (≈ `rate_det` to ~1e-3) |
| `rho_50`       | (n_t,) | signal-median radius (50 % of `rate_det` enclosed) |
| `rhot_min`     | (n_t,) | chord closest-approach rhot (KM9; NaN for KM14) |
| `rhot_crit`    | (n_t,) | enclosed-shell rhot (KM14; NaN for KM9) |
| `ok`           | (n_t,) | per-slice success flag (False = slice skipped) |

0-D string / scalar metadata: `pulse`, `runid`, `diagnostic`, `channel`,
`th_var`, `eq_source`, `chord_kind`, `los_file`.

Default output path: `src/tmp/<run_id>_<DIAG>_los_weights.npz` (override with
`--output`).

## Note on per-slice rhot grid

`X(t)` in TRANSP is nominally time-dependent (rezoned). In practice for
104614 M29 the slice-to-slice shift is at the ~1e-3 level, but **`rhot` is
stored per slice** so the file is correct under rezone-heavy runs without
relying on a "close enough" common grid. Consumers index `z['rhot'][it]` per
slice; if you need a single 1-D `rhot`, average across `it` after checking
`rhot.std(axis=0).max()` is small enough for your tolerance (a `--common-rhot`
mode could be added if a use case appears — not done preemptively).

## Channel-independence and BT

`f_det` depends only on the LoS cells (geometry) and the equilibrium
(`rhot(R,Z)` map). The chosen `--channel` only controls which `THNTX_<chan>`
is **co-stored** as the reference profile / closure check — switching channel
**does not** change `f_det`. For BT in `jet_tritium`, apply the same `f_det`
to the BT emissivity profile; CONTEXT validates this for KM14 (vertical chord
sees BT as isotropic to <0.5 %) and KM9 (v1 BT isotropic).

If KM9 grows a directional BT model later, a separate `f_det_bt` would be
needed; the script's split-by-channel design already accommodates that
without breaking the existing key layout.

## CLI

```bash
cd ~/jet/analysing_jet/src/neutron

# Every TIME3 slice in the run, KM9, CDF equilibrium (default)
python make_los_weights.py 104614 M29 --diag km9

# KM14, restricted time window 53.0..54.0 s, every 5th slice
python make_los_weights.py 104614 M29 --diag km14 -t 53.0 54.0 --stride 5

# PPF equilibrium (Heimdall / freia)
python make_los_weights.py 104614 M29 --diag km9 \
    --eq-source ppf --dda eftp --uid jetppf --seq 0

# Custom output path
python make_los_weights.py 104614 M29 --diag km9 --output /tmp/km9_w.npz
```

Flags: `--diag {km9,km14}` (required); `--channel {total,dd,dt}` (default
`total`); `--eq-source {cdf,ppf}`; `--dda --uid --seq` (PPF mode);
`--data-dir`; `--los-file` (override the diagnostic-default `.los`);
`-t/--time` (1 value = nearest TIME3; 2 values = inclusive window in JET s);
`--stride N`; `--chord-kind {auto,vertical,horizontal}` (override the
diagnostic default); `--no-subgrid`; `--no-median`; `--output`.

A bad slice (equilibrium failure on a ramp slice, NaN griddata, ...) is
caught per-slice, marked `ok=False`, and the loop continues — masking with
`z['ok']` on the consumer side filters those out.

## Performance

The cubic-griddata equilibrium build is the per-slice bottleneck (~3–6 s on
the WSL dev env, dominated by the 200×600 grid `los_file_detector_rate`
constructs for its subgrid splat). A 100-slice TIME3 grid is therefore
~5–10 min; `--stride` thins the work proportionally if every slice is not
required (the underlying `f_det` shape varies slowly outside transients).

## How `jet_tritium` consumes it

Minimal recipe (no external deps beyond NumPy / SciPy):

```python
import numpy as np
from scipy.interpolate import PchipInterpolator

z = np.load('104614M29_KM9_los_weights.npz')
ok = z['ok']
t_w   = z['times_jet'][ok]    # weight time grid
rhot  = z['rhot'][ok]         # (n_t, n_rho)
fdet  = z['f_det'][ok]        # (n_t, n_rho)
dvol  = z['dvol_m3'][ok]      # (n_t, n_rho)

# For each consumer time slice, PCHIP eps(rhot_consumer) -> z['rhot'][it_w],
# then take the inner sum:
def los_weighted_rate(eps_rhot, rhot_eps, it_w):
    eps_i = PchipInterpolator(rhot_eps, eps_rhot, extrapolate=False)(rhot[it_w])
    eps_i = np.where(np.isnan(eps_i), 0.0, eps_i)
    return float(np.sum(eps_i * fdet[it_w] * dvol[it_w]))
```

If `jet_tritium`'s emissivity is itself read off the same TRANSP `X` grid no
PCHIP step is needed — just `np.sum(eps * fdet[it] * dvol[it])`.

For an apples-to-apples sanity check against the file's own
`rate_det[it]`, feed in `z['thntx_si'][it]` on `z['rhot'][it]` — the inner
sum returns `rate_closure[it]` (≡ `rate_det[it]` to ~1e-3).

## See also

* [`thermal_los_rate.md`](thermal_los_rate.md) — the `f_det = C_bin / DVOL`
  construction; algorithmic pinning and subgrid binning.
* [`km9.md`](km9.md) — horizontal-chord specifics (`rhot_min`, the 5-bin
  median, no `rhot_crit`).
* [`compare_los.md`](compare_los.md) — how the two `f_det` weights relate when
  the same TH/BT plasma is observed by both chords.
* [`shared_los_library.md`](shared_los_library.md) — the `los_common` surface
  this driver reuses.
