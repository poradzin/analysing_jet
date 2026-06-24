# CONTEXT — neutron LOS analysis (orientation)

Working notes for the **neutron line-of-sight (LOS)** analysis in
`src/neutron/`. This file is the lean session-startup orientation; the
detailed material that used to live here has moved:

* **How to run / code structure / per-script summaries** →
  [`src/neutron/README.md`](src/neutron/README.md).
* **Physics & numerics deliberations** (algorithms, bug investigations,
  validation results) → [`src/neutron/doc/`](src/neutron/doc/) — see the
  index at the bottom of the README.
* **Cross-session facts & pitfalls** → auto-memory (`MEMORY.md`): KM14
  diamond/Nocente analysis, ZBND axis-order bug, PSI reshape, CDF axis
  psin floor, F_D_NBI-vs-BDENS_D norm, NUBEAM zone-grid pitfalls, ACFILE
  averaging. Verify any file/flag a memory names still exists before
  acting on it.

`src/neutron/` is intended to become its own standalone repo.

## What this is

Line-integrated thermal-neutron rates and thermal/beam-target (TH/BT)
ratios for the JET neutron diagnostics **KM14** (vertical diamond
spectrometer) and **KM9 / MPRu** (horizontal MPR), from TRANSP/NUBEAM
output, plus a KM14 diamond neutron-spectrum forward model. Everything is
self-contained on the TRANSP CDFs (no `ppf` by default), so it runs in the
WSL dev env. Shared code lives in `common/los_common.py`; `km14/` and
`km9/` hold thin per-diagnostic scripts.

## Current state (validated)

* **BT physics chain validated end-to-end**: in-house 4π reactivity
  reproduces NUBEAM `BTN4`/`BTN1` to ≈1.00; BT emission is **isotropic to
  <0.5 % for KM14's vertical chord** (g≈1.00), so the 4π per-zone
  emissivity is the right thing to line-integrate.
* **Headline numbers (104614 M29, t≈53.53 s / idx 2, total channel):**
  KM14 `(TH/BT)_LOS = 0.529` vs full-plasma 0D 0.433; KM9
  `(TH/BT)_LOS = 0.510` vs 0.433. KM14 box `rho_bnd = 0.3585` (wide chord
  2.60/3.16) / `0.3070` (narrow 2.70/3.10); real-cell `rho_50 ≈ 0.252`.
  KM9 `Rate_det = 4.34e9 n/s`, `rho_50 = 0.264`.
* **KM14 spectrum vs data**: M29 TRANSP-LOS already fits the measured
  diamond spectrum (best χ²/N = 3.07) with no TH:BT tuning. The
  tritium-fraction story is settled — Nocente's low core TH/BT came from
  run M13's high T/(D+T)≈0.53; M29's realistic ≈0.31 matches the data
  (~0.59). LOS weighting is a secondary ~15 % effect. See
  [`doc/km14_spectrum.md`](src/neutron/doc/km14_spectrum.md).
* **Shared library + KM9** built (Stage B done 2026-06-19): `los_common.py`
  factored out; KM9 `los_thermal_rate.py` / `los_th_bt_ratio.py` /
  `plot_LoS.py` written; time-window averaging (`-t`/`--idx`) added. KM14
  numbers reproduce byte-for-byte after the refactor.
* **KM9-vs-KM14 real-LOS comparison** (`compare_los.py`, added 2026-06-24):
  one driver overlays both diagnostics' real-LOS weight, cumulative TH/BT and
  per-shell signal (no box); reuses `los_common`, reproduces both scripts
  exactly. Key insight: the two cumulative TH/BT coincide on the core and
  diverge only past rhot ≈ 0.4 — a weighted average of the *shared* local
  ratio is weight-insensitive where that ratio is flat (core), so the chords'
  different core/edge "ballast" (KM14 banks 69 % of its BT signal inside
  rhot 0.4 vs KM9's 63 %) only surfaces on the structured edge. See
  [`doc/compare_los.md`](src/neutron/doc/compare_los.md).

## Open items / next steps

* **Known forward-model limit**: BT spectral peak sits ~30–50 keV high on
  M13 — our symmetric-Gaussian detector convolution lacks the diamond's
  asymmetric low-energy tail. Fix = tabulated response `R(E_dep|E_n)`
  (MCNP/measured), folded per shell. Scoped, deferred. Report TH/BT
  *ratios* (peak-insensitive) until then. See `doc/km14_spectrum.md`.
* **KM9**: directional `--bt-aniso` for the horizontal chord (deferred;
  BT is isotropic-treated in v1).
* **Refactor follow-ups**: swap `bzi.find_run_dir` → `los_common` in
  `km14_spectrum.py` / `bt_los_emissivity.py` / `bt_poloidal_distribution.py`;
  rename the `los_file_detector_rate` dict key `thkm14_det_si` →
  `thlos_det_si`; move `bt_zone_integrator` into `common/` (km9 currently
  `sys.path`s `../km14` to import it). See `doc/km9.md` Stage-B follow-ups.
* **External validation**: FIDASIM benchmark of the absolute angular
  machinery (low priority for KM14 since g≈1, still wanted for tangential
  channels like KN3); compare KM9 results against measured signals.
