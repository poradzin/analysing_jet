# `los_th_bt_ratio.py` — context (added 2026-06-03, commit 3b of step 2)

`src/neutron/km14/los_th_bt_ratio.py` computes the **KM14 line-of-sight
TH/BT neutron ratio** for a TRANSP run — the quantity to compare against
the spectroscopically-measured KM14 TH/BT. It is the beam-target
extension of `los_thermal_rate.py`, but **fully self-contained on the
TRANSP CDFs** (no `profiles.Eq` / ppf), so it runs in the WSL dev env as
well as on freia. This is justified because commit 3 showed the BT
*emission direction* is isotropic to <0.5% for KM14, so the 4π per-zone
emissivity is the correct thing to line-integrate.

## Results

**104614 M30 idx 1, DD-only channel** (2.45 MeV, what KM14 spectroscopically separates):
```
(TH/BT)_LOS              = 0.353   (flux mode) / 0.367 (zone mode)
(TH/BT) whole plasma 0D  = 0.310
core enhancement LOS/0D  = 1.14x
```
KM14 sees a ~14 % higher thermal fraction than the volume-averaged ratio
— the chord weights the core differently from the total yield. The
poloidal-asymmetry correction (flux- vs zone-mode) is only ~4 %.

**104614 M30 idx 2 (t_TRANSP=13.33 s), `total` channel** (default,
added 2026-06-09 — see "Channels" below):
```
(TH/BT)_LOS              = 0.505
(TH/BT) whole plasma 0D  = 0.415
core enhancement LOS/0D  = 1.22x
TH whole-plasma          = 3.99e17 n/s  (matches los_thermal_rate.py exactly)
BT whole-plasma          = 9.61e17 n/s  (BTNTS_DD+BTNTS_DT = 9.63e17, 0.1% gap)
                                          BTNTS_DD = 8.11e15  (0.8%)
                                          BTNTS_DT = 9.54e17  (99.2%)
```
For a DT-campaign run THNTX is dominated by THNTX_DT (~159× THNTX_DD on
M30); the "total" channel therefore differs from the DD-only channel by
~160× on the chord and whole-plasma rates.

## Cross-validation with `los_thermal_rate.py` (2026-06-09)

At t_TRANSP=13.33 s on M30, `los_th_bt_ratio.py --channel total` and
`los_thermal_rate.py` agree on the **thermal whole-plasma rate to
4 significant figures** (both 3.988e17 n/s) and on the **TH chord rate
to ~1 %** (4.58e15 vs 4.63e15 n/s). The residual is the equilibrium
source: `los_th_bt_ratio.py` uses TRANSP `PSIRZ` self-contained;
`los_thermal_rate.py` uses PPF EFTP via `profiles.Eq`. This confirms
the chord-integration math is identical in the two scripts and that
the earlier 160× discrepancy was the channel-variable difference, not
a bug.

## Method (CDF-only)

* **Equilibrium** `CdfEquilibrium`: psi(R,Z) from `PSIRZ` (C-order (Z,R)),
  normalized by `PSI0_TR`/`PLFLXA`; rhot(psi_n) inverted from `PLFLX` vs
  `XB`. Magnetic axis from `RAXIS`/`YAXIS`. LCFS polygon from `_fi`
  `RSURF/ZSURF[-1]`. The chord grid uses **`rhot_pinned`** (cubic `griddata`
  of ψn with the axis pinned at ψn=0) + **axis insertion into the R,Z arrays**,
  not plain bilinear `rhot`, so `f(rhot)` reaches ~1 on axis (the 2026-06-12
  CDF axis-floor fix; identical recipe to `los_thermal_rate.EqCDF`). Plain
  bilinear `rhot` is kept for ad-hoc point lookups. Ratios change <0.1%.
* **Thermal** eps_TH(rhot): the flux-function `THNTX[_DD|_DT]` profile
  selected by `--channel` (see below).
* **Beam-target** eps_BT(rhot): per-zone `BTN4`/`BTN1`/`BTN5`/`BTN7`
  summed across the keys the active channel selects, then either
  `--bt-mode flux` (default, BMVOL-weighted flux-surface average) or
  `--bt-mode zone` (per-zone cubic griddata onto the chord, keeps
  poloidal asymmetry — ~4% effect).
* Both mapped onto the KM14 chord grid, chord-integrated. **For the ratio
  the chord geometry (w_tor, 2πR, solid angle) cancels**, so it's robust
  and equilibrium-choice-insensitive (TH and BT see the same surfaces).
* Cross-checks printed: full-torus ratio (== chord ratio for flux mode),
  0D whole-plasma ratio, BT whole-plasma vs `sum(BTN·BMVOL)` and the
  matching `BTNTS_*` scalars (sum + per-component breakdown when the
  channel spans more than one BT key).

## TH/BT weight function vs rhot (commit 3c, 2026-06-12)

`rhot_weight_profiles()` adds the rhot-resolved weight function / TH-BT
breakdown, mirroring `los_thermal_rate.py`. The LOS weight **`f(rhot)` is
purely geometric** — the fraction of each TRANSP flux shell's toroidal volume
inside the chord R-band — so it is **identical for TH and BT** (flux mode); it
reuses `los_thermal_rate.los_shell_fraction` (lazy import; the two modules
import each other, so the import is inside the function to avoid a cycle).
Everything is put on the TRANSP `X` grid; the BT flux profile (always the
BMVOL flux-average, even under `--bt-mode zone`, since the weight function is a
flux-surface quantity) is PCHIP-interpolated onto it via `_interp_flux_to`.

Quantities (also written to `src/tmp/<run_id>_KM14_LOS_THBT_weight_idx<idx>_t<t>s.txt`
on `--save`): `f(rhot)`; LOS-weighted `THKM14 = TH·f`, `BTKM14 = BT·f`;
per-shell LOS rate `TH·f·DVOL`, `BT·f·DVOL` [n/s]; **local ratio**
`TH(x)/BT(x)` (f and DVOL cancel); **cumulative ratio**
`cumΣ(TH·f·DVOL)/cumΣ(BT·f·DVOL)` whose **rhot=1 endpoint closes to
`(TH/BT)_LOS`** (verified ~0.1%: M30 idx1 total 0.5060 vs 0.5056; dd 0.3528 vs
0.3529 — per-shell sums match the (R,Z) toroidal integrals to ~0.7%, i.e. no
BT normalization gap). Kept in the **TH/BT** convention (BT/TH = inverse). In
`--bt-mode zone` the endpoint uses flux-averaged BT, so it differs from the
asymmetric zone `(TH/BT)_LOS` by the ~3% poloidal asymmetry (the report says so).

## Channels (updated 2026-06-09)

* **`--channel total` (default)** — unseparated `THNTX` (DD+DT+...) vs
  the sum of every BT component present in `_neut` (DD+DT+TT+TD,
  whichever exist). This is the unfolded total neutron rate KM14 sees
  if no spectroscopic separation is applied. **Use this to compare
  against `los_thermal_rate.py`** (which also reads the unseparated
  `THNTX`).
* **`--channel dd`** — `THNTX_DD` vs `BTN4` (2.45 MeV), the channel
  KM14 spectroscopically separates on DT-campaign pulses.
* **`--channel dt`** — `THNTX_DT` vs `BTN1` (14 MeV).

NB (corrected 2026-06-15): **104614 M29 is a DT run**, not pure DD — the local
`~/jet/data/104614/M29` files have `THNTX_DT` (~140× `THNTX_DD`), thermal `NT`,
`BTNTS_DT`, and `BTN1`+`BTN4` in `_neut`. The KM14 diamond detector measures the
14 MeV DT line (E_dep = E_n − 5.7 MeV via ¹²C(n,α₀), peak 8.4 MeV), so the
**`dt` channel is the physical one for KM14 on this pulse** (earlier "M29 pure DD"
notes were wrong).

## BT emission anisotropy `--bt-aniso` (added 2026-06-15)

Optional finite-height-detector treatment of the beam-target **emission
direction**, broadening the script beyond the 4π-isotropic default. KM14 is a
point detector ~10 m above `Zmag` on the chord axis, so BT neutrons reach it
within a **narrow upward cone** whose axis tilts a few degrees off vertical
across the chord (per emission point). Each NUBEAM zone's `BTN*` rate is
weighted by the directional factor `g(zone) = (dε/dΩ toward the detector) /
(ε/4π)`, computed by reusing `bt_los_emissivity` (`BField` + per-zone MC of the
fast-ion dist + CM→lab boost). **Only the dimensionless `g` is taken from the
MC**; the 4π magnitude stays NUBEAM's validated `BTN*`. In the TH/BT ratio the
common detector solid angle cancels, so the exact net effect is
`(TH/BT)_LOS → (TH/BT)_LOS / ⟨g⟩_BT`.

* **Geometry.** Per-zone LOS `n̂ = (R_det−R, 0, Z_det−Z)` (zero toroidal
  component, since detector and chord share a toroidal plane). Detector at
  `--detector-R` (default chord centre 2.90 m), `--detector-height` above Zmag
  (default 10 m). `zone_emissivity` in `bt_los_emissivity.py` was generalized to
  accept a per-zone `(n_zone,3)` `n_los` (back-compatible with the single `(3,)`
  vector its own `main` passes — selftest + M29 run unchanged).
* **Implementation.** `bt_anisotropy_factors()` builds `g_map[key]` per BT
  component (DD→`DD_n`, DT→`DT`; keys without a Bosch–Hale σ such as TT get
  `g=1`); `main` forms `bt_zone_eff = Σ_k BTN_k·g_k` and runs it through the
  *existing* flux/zone grid + chord machinery, so flux-/zone-mode, the `f(rhot)`
  weight function and the cumulative-ratio closure all carry through unchanged.
  The 0D whole-plasma cross-check keeps the isotropic `bt_zone_sum`.
* **Result (matches the commit-3 prediction g≈1.00).** The ~5° poloidal tilt
  stays ⟂ to the toroidal B̂, so the correction is negligible:
  - 104614 **M29 idx 1 dd**: `(TH/BT)_LOS` 0.3744→0.3742, `⟨g⟩_BT=1.0008`
    (−0.08%); reactivity-wtd `g_DD=1.0013` (≈ the fixed-+z 1.001 from
    `bt_los_emissivity`, confirming the finite-height tilt is immaterial).
  - 104614 **M30 idx 1 total**: 0.5056→0.5059, `g_DD=1.0016`, `g_DT=0.9981`,
    `⟨g⟩_BT=0.9995` (+0.05%, DT-weighted).
  Reporting prints both ratios, `⟨g⟩_BT`, the % shift, per-component `g` and the
  `vector Eps_4pi/NUBEAM` consistency (≈0.98–0.99). `--save` tags the weight
  `.txt` and plot label with `_aniso`/`[BT aniso]` so they don't overwrite the
  isotropic run. New flags: `--bt-aniso --detector-R --detector-height
  --cone-deg --nsamp --fast-norm {bdens,ntot} --no-rotation --seed`.

## CLI
```bash
python src/neutron/km14/los_th_bt_ratio.py 104614 M30 --idx 2 --plot       # total (default)
python src/neutron/km14/los_th_bt_ratio.py 104614 M30 --idx 1 --channel dd # DD-only
python src/neutron/km14/los_th_bt_ratio.py 104614 M30 --bt-mode zone
python src/neutron/km14/los_th_bt_ratio.py 104614 M29 --idx 1 --channel dd --bt-aniso  # finite-height BT view
```
Flags: `--idx --data-dir --channel {total,dd,dt} --bt-mode {flux,zone} --Rmin --Rmax --wtor --nR --nZ --plot --save --no-plot`
plus (anisotropy) `--bt-aniso --detector-R --detector-height --cone-deg --nsamp --fast-norm {bdens,ntot} --no-rotation --seed`.
Plot (2×4): top row eps_TH(R,Z), eps_BT(R,Z), local TH/BT, R-integrated vs Z;
bottom row (commit 3c) `f(rhot)`, per-shell LOS rate `eps·f·DVOL`, local
`TH/BT(rhot)`, cumulative `TH/BT(rhot)`. The cumulative panel overlays the
LOS-weighted curve (cyan, → `(TH/BT)_LOS` at rhot=1) and a dashed TRANSP
full-plasma reference `cumΣ(TH·DVOL)/cumΣ(BT·DVOL)` (no `f`, → whole-plasma 0D
ratio at rhot=1) so the LOS-narrowing reweighting is read off directly.
The local-TH/BT (R,Z) map masks near-LCFS cells where BT→0 and clips the colour
scale to the 2–98th percentile (else edge spikes wash it to one flat colour).
`--save <png>` also writes the weight-profile `.txt` to `src/tmp/`.

## Effective-weight figure & cumulative-ratio plateau (added 2026-06-18)

`make_plot` now also produces a **second figure** (1×3) with the three
weight PDFs and their action on TH/BT. Each weight (`DVOL`, `f·DVOL`,
`f_det·DVOL`) is normalized to a *true* PDF in rhot — divided by its
trapezoidal integral over [0, 1] so `∫₀¹ w(ρ) dρ = 1`. The three panels:

1. The three normalized weights (DVOL = grey, `f·DVOL` = green geometric
   LOS, `f_det·DVOL` = magenta real LOS) with `rhot_crit` marked.
2. `TH·PDF_w` for each of the three weights — per-shell TH contribution
   under each weighting, directly comparable since all PDFs share the
   same area-1 normalization.
3. Same for BT.

`--save` writes the figure to `<save>_weights.<ext>` next to the main PNG.

**Why cumulative TH/BT is identical up to rhot ~ 0.26** (104614 M29 dd,
observed across all three weightings). Inside `rhot_crit ≈ 0.12` the
geometric `f` is pinned to 1 (CONTEXT lines ~242-280) and `f_det` is
pinned to its flat plateau (lines ~520-540). For any *constant* weight
`w₀`, the cumulative ratio reduces algebraically:

    ∫₀^ρ TH·w₀·DVOL dρ' / ∫₀^ρ BT·w₀·DVOL dρ' = ∫₀^ρ TH·DVOL / ∫₀^ρ BT·DVOL

i.e. the constant `w₀` cancels and all three weightings give *exactly*
the same cumulative TH/BT inside any constant-weight region. The
agreement extends past `rhot_crit` to ~0.26 because `f` rolls off only
~1/r and `f_det` rises slowly with solid angle, so `w(ρ)` stays close
enough to constant that the cumulative ratio hasn't moved measurably yet.

**Why the cumulative drops 0.6 → 0.54 from ~0.26 to rhot=1 despite the
rate falling fast.** The drop is a *relative* 10 %, which by a simple
mass balance requires the outer contribution to BE non-negligible. Let
`A_TH/A_BT = 0.6` be the inner totals up to ρ_split and `δTH/δBT =
r_outer` the outer slab. Then

    (A_TH + δTH)/(A_BT + δBT) = 0.54
    => δBT / A_BT = 0.06 / (0.54 - r_outer)

so for `r_outer ≈ 0.1` the outer BT yield is ~14 % of the inner total,
for `r_outer ≈ 0.3` it's ~25 %. *Not* small. The outer BT is sizeable
because `BTN(ρ)·DVOL(ρ)` is broader than `THNTX·DVOL` — fast-ion
deposition is shoulder-like, not core-peaked, and the volume Jacobian
(DVOL grows ~ρ from perimeter × shell width) fights the emissivity
falloff. The right per-shell weight for the *ratio* is `ε·DVOL`, not
`ε`, and `BT·DVOL` peaks well off-axis on this kind of profile. Panel 3
of the new figure (BT × PDF) makes this visually obvious: the area
under those curves between ρ ~ 0.26 and 1 is what drives the drop.

Practical reading: the cumulative TH/BT *under any of these three
weightings* is insensitive to LOS-weight choices wherever the weight is
locally flat (i.e. on the enclosed core, ρ < ~0.26 here). The
LOS-vs-full-plasma differences in the ratio show up only in the outer
slab, weighted by `f`/`f_det` there — which is the regime in which
choice of weighting actually matters. For ratio-only analyses on
shells fully inside the chord, picking DVOL, `f·DVOL` or `f_det·DVOL`
gives the same answer; pick whichever is cheapest.

## Open / next

* **Run on the M29 pure-DD case** (the actual KM14 thermal-analysis run,
  pulse 104614 M29) and **compare against the measured KM14 TH/BT.** M29 is
  now **local** at `~/jet/data/104614/M29` (`.CDF` + `_fi`/`_neut` idx 1–3),
  so `python los_th_bt_ratio.py 104614 M29 --idx 1 --plot` runs in the dev env.
* Decide whether to fold the commit-3 angular factor g (≈1.00 for KM14) in
  explicitly — currently omitted as <0.5%.
* Optional: detector response (energy window/scattering) — same open
  question flagged for the thermal channel; currently pure emissivity.
* Time dependence: BT only exists at FBM idx times; for a TH/BT *trace*
  loop over available `_fi`/`_neut` indices.

---
