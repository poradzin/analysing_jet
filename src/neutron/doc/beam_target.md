# `bt_poloidal_distribution.py` — context (added 2026-05-27)

`src/neutron/km14/bt_poloidal_distribution.py` is a standalone script
(no `profiles.py` / `ppf` dependency, reads TRANSP CDFs directly via
`netCDF4`) that maps the **spatial** poloidal-angle distribution of
beam-target neutron emissivity. It is **step 1** of the plan to fold BT
neutrons into the KM14 LOS framework — it answers "is ε_BT flux-function
or does it vary with θ_pol on a flux surface?" and nothing else. It does
not produce the *emission-direction* angular distribution; that is
step 2 below.

## Inputs

* `<run_dir>/<runid>_fi_<idx>.cdf` — NUBEAM zonal geometry only:
  * `X2D`, `TH2D`, `R2D`, `Z2D`, `BMVOL` — per-zone scalars (220 zones
    on the standard NUBEAM 10-row × {4,8,…,40}-θ grid for 104614)
  * `NTHZSM` — cumulative zone counts per x-row
  * `RSURF`, `ZSURF` — flux surfaces for plotting (and `RSURF[0,:]`,
    `ZSURF[0,:]` gives the magnetic axis at xsurf=0)
  * `XSURF`, `THSURF` — the *finer* flux-surface grid (dx=0.05) used by
    `RSURF/ZSURF`, **not** the zonal row partition (which has dx=0.1)
  * F_D_NBI etc. are intentionally not read here
* `<run_dir>/<runid>_neut_<idx>.cdf` — per-zone neutron emissivity:
  * `BTN4` (DD-BT, always present), plus optionally `BTN1` (DT),
    `BTN5` (TT), `BTN7` (TD)
  * `THNTNT2d`, `TOTNTNF2d` — for sanity cross-checks
  * `TA` — time of the FBM averaging window

## Data-directory search order

1. `--data-dir <base>/<pulse>/<run_suffix>` if given
2. `~/jet/data/<pulse>/<run_suffix>` (local WSL)
3. `/common/transp_shared/Data/result/JET/<pulse>/<run_suffix>` (heimdall)

## Algorithm

1. List `<runid>_fi_*.cdf` files, pick the requested `--idx` (default =
   first available).
2. Read NUBEAM zone geometry from `_fi`, BT components from `_neut`. The
   per-zone (R, Z) values differ slightly (~10–17 cm) between the two
   files but the **zonal indexing is the same** — the 220-entry arrays
   in both files follow the same NUBEAM (x-row, θ-bin) ordering. We
   trust the indexing and pair `BTN4(zone)` directly with
   `TH2D(zone)` from `_fi`.
3. Sum BT components → `bt_total(zone)` in 1/cm³/s.
4. Compute per-zone rate `bt_total * BMVOL` (1/s), poloidally-averaged
   emissivity per x-row, and the up–down asymmetry
   `[ε(θ)−ε(−θ)]/(ε(θ)+ε(−θ))` per row (interp onto θ-sorted grid).

## CLI

```bash
python src/neutron/km14/bt_poloidal_distribution.py 104614 M30
python src/neutron/km14/bt_poloidal_distribution.py 104614 M30 --idx 2
python src/neutron/km14/bt_poloidal_distribution.py 104614 M30 \
       --data-dir /common/transp_shared/Data/result/JET --save bt.png
```

Flags: `--idx`, `--data-dir`, `--no-plot`, `--save`.

## Plot layout (2×3)

| | col 0 | col 1 | col 2 |
|---|---|---|---|
| row 0 | ε_BT(R,Z) scatter + LCFS + axis | ε_BT(θ_pol) curves, one per x-row | up–down asymmetry vs θ |
| row 1 | ε_BT(x, θ_pol) scatter | vol-averaged ε(x) + per-row rate (twin axis) | BT component breakdown vs x (log) |

## Pitfalls / lessons

1. **XSURF ≠ zone partition.** XSURF has 21 boundaries (dx=0.05) for the
   finer RSURF/ZSURF grid; zones live at 10 unique X2D values (dx=0.1),
   one per row defined by NTHZSM. Always take row centres from
   `X2D` itself, not from `XSURF`.
2. **(R,Z) mismatch between `_fi` and `_neut`.** Up to ~17 cm difference
   on the same zone index; the user's earlier `TRANSP_read_plot_GENERAL_DT.py`
   silently mixes them. Treat the per-zone *index* as the link, not the
   (R,Z) values; if you need (R,Z), use the `_fi` ones for plotting and
   accept the small mismatch, or interpolate.

## NUBEAM MC grid — authoritative reference

Pankin et al., *Comp. Phys. Comm.* **312** (2025) 109611, Section 6.1.1
("Monte Carlo grid"), Fig. 6
(`~/jet/documentation/2025_Pankin_TRANSP_overview_Comp_Phys_Comm_312_109611.pdf`).
Confirms the convention used above:

* Irregular 2-D grid of "zone rows" aligned with flux coordinates,
  equally spaced in ξ; each row subdivided into a different number of θ-zones,
  fewer near axis, more near edge — "designed so all zones have roughly
  equal volume and cross-sectional area" (explains the {4,8,…,40} pattern
  for the 10-row 104614 grid).
* Poloidal zones stored contiguously, index increasing with θ
  counter-clockwise in the plasma cross-section drawn to the right of
  the machine axis. θ-start can be `0`, `-π`, or default (`0` for
  up-down symmetric, `-π` for asymmetric); for 104614 (X-point divertor)
  TH2D ∈ [−π, π].
* Same grid hosts F_D_NBI, halo neutral sources, BT/BB fusion rates.

Table 8 of the same paper enumerates the fusion products and reactions
(`BTN4` = DD→³He + n, `BTN1` = DT, `BTN5` = TT→2n, `BTN7` = TD) — the
reference list we sum in `bt_total`. Useful as the cross-section ground
truth for step 2's Bosch–Hale folding.

## Smoke test (104614 M30, idx 1, t=12.33 s TRANSP time)

```
DD-BT  : 8.17e15 1/s  ( 0.8%)
DT-BT  : 9.57e17 1/s  (99.2%)
TH     : 3.99e17 1/s
Mag axis (R,Z): (3.037, 0.285) m   -- significant upward shift
```

DT dominance is expected for an M30 DT-campaign run; M29 (used by
`los_thermal_rate.py`) is pure DD.

---

# Step 2: in-house BT angular emissivity code (planned, 2026-05-27)

The next step beyond `bt_poloidal_distribution.py` is to compute
ε_BT,LOS(R, Z, n̂) — the BT emissivity along the KM14 LOS direction —
which `_neut` does *not* contain (BTN4 is 4π-integrated). The plan is to
**write a standalone in-house module** and benchmark it against FIDASIM
(git-cloned to `~/jet/FIDASIM*`). NEMO is not available (appears to live
in an ITER-restricted repo).

## Algorithm sketch (per NUBEAM zone)

1. Read `F_D_NBI(zone, ξ, E)` from `_fi_<idx>.cdf` (ξ = v∥/v in the
   local B frame, E in eV, units `#/cm³/eV/d(Ω/4π)`).
2. Read thermal target n_D(x), T_i(x) from the main `<runid>.CDF`
   (`ND`/`NDB` and `TI`); interpolate onto the per-zone x.
3. For each (ξ, E) bin: unfold uniform gyrophase → ring of fast-ion
   velocity vectors with axis along B̂(R, Z).
4. Convolve with thermal Maxwellian to get the relative-velocity
   distribution.
5. DD/DT cross section via Bosch–Hale parameterisation; CM-frame
   neutron angular distribution (DD ≈ isotropic in CM to ~5 %, DT has
   stronger angular dependence at the energies of interest).
6. CM→lab transform: lab neutron direction = v_CM + v_n,CM where
   |v_CM| ≈ v_beam/2 (v_CM/v_n,CM ≈ 0.07 for 100-keV D on a few-keV
   target — modest forward boost but 10–30 % lab anisotropy that KM14,
   roughly ⊥ to tangential beams, sees as a suppression vs isotropic).
7. Project onto n̂_LOS and accumulate.
8. Sanity-check: 4π integral per zone must reproduce `BTN4(zone)`.

## Benchmarks

* Per-zone 4π integral vs `BTN4` from `_neut`.
* Full LOS integral against FIDASIM `weights`-mode output for the same
  (run, time, channel). FIDASIM input prep needs the TRANSP CDF + an
  EFIT-derived equilibrium in FIDASIM's expected format — there are
  IDL/Python preprocessors in the FIDASIM tree (`lib/python/`).

## Open scoping questions

* Use the existing zone grid as-is, or resample F_D_NBI onto a regular
  (R, Z) grid before kinematics? (NUBEAM zone-centred is cheaper and
  natural for the BTN4 sanity check.)
* Whether to include finite-orbit corrections beyond NUBEAM's already
  orbit-averaged F_D_NBI (probably not; F_D_NBI already encodes
  banana-orbit topology via the (ξ, x, θ) sampling).
* KM14 detector response folding (energy window, scattering) — same
  open question as for the thermal channel.

## Build plan: three commits

The step-2 module is being built inside-out so each piece is
independently testable against an external reference (cross-section
tables → BTN4 → FIDASIM):

1. **Commit 1 — `bt_kinematics.py` (pure-physics module, done 2026-05-27).**
   Bosch–Hale σ(E_cm), CM-frame angular sampler, classical 2-body
   CM→lab kinematics. No TRANSP / equilibrium dependencies. See
   dedicated section below.
2. **Commit 2 — single-zone integrator (next).** Read F_D_NBI for *one*
   zone, equilibrium B̂(R,Z) from the main CDF, n_D and T_i for that
   zone; gyrophase-unfold; sample Maxwellian thermal partner; fold with
   `bt_kinematics`. Acceptance test: 4π integral reproduces
   `BTN4(zone)` to within MC statistics.
3. **Commit 3 — full LOS pipeline + FIDASIM benchmark.** Loop over
   zones, project onto n̂_KM14, compare against FIDASIM `weights`-mode.

---

# `bt_kinematics.py` — context (added 2026-05-27, commit 1 of step 2)

`src/neutron/km14/bt_kinematics.py` is the pure-physics core of the
in-house BT angular-emissivity module. Self-contained — only depends on
`numpy` (and `matplotlib` for the optional `--plot` flag). No
TRANSP/NUBEAM/equilibrium imports.

## What it provides

* **Cross sections** σ(E_cm) in mb, from Bosch & Hale 1992 Table VII
  parameterisation. Returns 0 outside the published validity range.
  * `sigma_dd_n_mb(E_cm_keV)` — D(d,n)³He, valid 0.5–4900 keV
  * `sigma_dt_mb(E_cm_keV)` — T(d,n)⁴He, valid 0.5–4700 keV. Same
    function covers the TD channel: σ depends only on E_cm, the
    "which-is-fast" label only enters the kinematics.
  * `sigma_tt_mb` — `NotImplementedError` placeholder. Negligible for
    DTE3 (D-only beams, e.g. 104614). For DTE2 (D and T beams),
    implement from ENDF/B-VIII or Drosg's parameterisation.
* **Reaction metadata** via `reaction_info("DD_n" | "DT" | "TD" | "TT")`
  — masses (amu), Q-value (MeV), σ function pointer.
* **CM-frame neutron emission sampler**
  * `isotropic_cm(n, rng=None)` — uniform samples of `(cos_θ, φ)`. DD
    and DT are isotropic in CM to ~5 % at our energies; the dominant
    lab anisotropy comes from the CM→lab boost, which is exact.
    Intrinsic-CM anisotropy hook left for later (Legendre-coefficient
    sampler swappable at the call site).
* **Kinematics** (classical, non-relativistic; ~0.1 % correction at
  14 MeV — ignored)
  * `neutron_lab_velocity(v_fast, v_th, reaction, cos_θ_cm, φ_cm)` —
    array-broadcast; polar axis for the CM emission is the fast
    reactant's velocity in the CM frame (so cos_θ_cm = 1 → forward
    along the beam-in-CM direction).
* **Helpers**: `relative_cm_energy_keV`, `neutron_lab_energy_keV`,
  `speed_from_energy_keV`.

## Unit tests (`python3 src/neutron/km14/bt_kinematics.py`)

Five σ spot-checks against hand-evaluated Bosch–Hale Table VII values
plus two kinematic-band tests:

```
sigma_DD_n( 50.0 keV) =   16.49 mb   (target  16.5 +/- 1.0)
sigma_DD_n(100.0 keV) =   37.01 mb   (target  37.0 +/- 2.0)
sigma_DT  ( 50.0 keV) = 4218.62 mb   (target 4220  +/-  50)
sigma_DT  ( 64.0 keV) = 5063.54 mb   (target 5063  +/-  50)   <- peak in E_cm
sigma_DT  (200.0 keV) = 1137.97 mb   (target 1138  +/-  50)

Kinematics: 100 keV D on cold D, 50k isotropic CM samples
  predicted lab E_n window: 2146.7 .. 2852.6 keV
  simulated min/mean/max  : 2146.7 / 2500.1 / 2852.6 keV   <- exact match

Kinematics: 100 keV D on cold T, 50k isotropic CM samples
  predicted lab E_n window: 13.432 .. 14.778 MeV
  simulated min/mean/max  : 13.432 / 14.106 / 14.778 MeV   <- exact match
```

`--plot` / `--save out.png` flag draws σ(E_cm) curves on log-log axes
with the DT peak marked.

## Caveats / things to revisit when integrating with commit 2

* DT cross section peaks at **E_cm ≈ 64 keV**, not E_lab — easy to
  mis-state. After the peak σ falls steeply (5.06 b at 64 keV → 1.14 b
  at 200 keV).
* `neutron_lab_velocity`'s polar axis convention (fast reactant in CM):
  irrelevant for isotropic sampling, but matters once we plug in an
  anisotropic Legendre sampler — most tabulations use this same
  convention (incident-particle direction in CM).
* The same `sigma_dt_mb` function and Q-value are used for both DT and
  TD; the kinematics correctly distinguish them via which reactant
  carries the fast velocity.

---

# Next steps (Thursday 2026-05-28)

User to start tomorrow morning by checking out commit 1 (read the
module, run the unit tests, optionally `--plot`). Then proceed to
**commit 2 — single-zone integrator**:

1. **Equilibrium B-field reader** from the main `<runid>.CDF`. Need
   B_R, B_Z, B_φ at any (R, Z) in the plasma. TRANSP's main CDF
   typically gives the poloidal flux ψ(R, Z) on its internal grid
   (variable name varies — likely `PLFLX` or similar; needs inspection
   like we did for `_fi`/`_neut`) plus the F = R·B_φ function vs ψ.
   From those:
     * B_R = -1/R · ∂ψ/∂Z
     * B_Z =  1/R · ∂ψ/∂R
     * B_φ = F(ψ)/R
   Return unit vector B̂(R, Z) for the gyrophase unfold.
2. **n_D(x), T_i(x) reader.** From main CDF: `ND` (or `NDB`) and `TI`,
   defined on the X-grid (sqrt-toroidal-flux). PCHIP-interpolate to the
   per-zone x from `X2D[zone]`. Watch for the "thermal D includes
   beam-slow-down fast ions" gotcha — TRANSP convention is usually
   `ND` = thermal, `NDB` = beam component; verify by checking the .CDF
   `long_name` attributes.
3. **Maxwellian thermal sampler.** 3-D Gaussian with σ_v = √(T_i/m_D),
   sampled per fast-ion sample. Vectorise.
4. **Per-zone integrator.** Read F_D_NBI for one zone. For each (ξ, E)
   bin with nonzero weight: sample N gyrophase angles; for each,
   compute v_fast vector. Pair with a Maxwellian thermal sample, get
   E_cm, σ, sample isotropic CM, get v_n_lab. Accumulate weighted
   contribution to ε_BT(zone, n̂).
5. **Acceptance test.** ∫(4π) of ε_BT(zone, n̂) dΩ ÷ n_D / σ_v
   should reproduce `BTN4(zone)` to within MC statistics. If it does,
   the kinematics chain is correct end-to-end.

Open questions to settle when starting commit 2:

* B̂ source: parse main `.CDF` directly, or via `profiles.Eq` (which
  needs ppf — but we have an alternative SAL path for local
  development)? Direct .CDF parse is more self-contained.
* MC sample budget per zone: 10⁴ enough for 4π sanity? Probably yes;
  the angular-distribution use case will need ≥10⁵.
* Where to land commit 2: same `src/neutron/km14/` directory, name
  `bt_zone_integrator.py` (or a clearer name suggesting "per-zone
  emissivity").

---

# `bt_zone_integrator.py` — context (added 2026-06-03, commit 2 of step 2)

`src/neutron/km14/bt_zone_integrator.py` is the per-zone beam-target
reactivity integrator. **Acceptance test PASSED**: the in-house 4π
reactivity reproduces NUBEAM's `BTN4` (DD) and `BTN1` (DT) on 104614 M30
idx 1 to **1.004 / 1.001** on the volume integral (per-zone medians
1.017 / 1.003). This validates the `bt_kinematics` chain end-to-end.

## Key simplification — no B-field needed in commit 2

For the *4π-integrated* per-zone rate the thermal target is an isotropic
Maxwellian, so `<σ·v_rel>` depends only on the fast-ion **speed**, not on
pitch ξ, gyrophase, or B̂. Hence:

```
ε_BT(zone) = n_fast(zone) · n_th(zone) · <σ(E_cm)·v_rel>     [1/cm³/s]
```

The whole PSIRZ→B̂ machinery (planned in the Thursday notes) is **not
needed here** — it only enters commit 3, where the emitted-neutron
*direction* (LOS projection) finally depends on the fast-ion velocity
vector. The B-field reader is therefore deferred to commit 3.

## The F_D_NBI vs BDENS_D normalization gap (the main finding)

`F_D_NBI` integrates (Σ F·dE·dξ·½) to `NTOT_D_NBI` *exactly*, but that is
**~19% below** the full beam-ion density `BDENS_D`:
`sum(n_fast·BMVOL)=8.18e19` vs `sum(BDENS_D·BMVOL)=1.01e20` (ratio 0.809).
NUBEAM computes `BTN4`/`BTN1` from the full beam density, so using the raw
F normalization gives a flat **0.81×** deficit across radius, *identical
for DD and DT* (different cross sections → rules out a kinematics error).

Diagnostic path that nailed it (all dead-ends documented so we don't
re-walk them):
* `<σv>` validated absolutely vs Bosch–Hale thermal DD reactivity
  (6.1e-19 / 2.6e-18 cm³/s at 10/20 keV) — kinematics correct.
* `BTNTOT4 = BTN4·BMVOL`; `BTN4` is the emissivity (1/cm³/s) we compare to.
* **Rotation is a red herring**: Ω≈5.3e4 rad/s → v_rot≈159 km/s, only ~6%
  of the fast-ion speed (2680 km/s). Toroidal-approx test gave only ±3%,
  and the co-injection (co-rotation) sign makes agreement *worse*
  (0.804→0.789). Not the cause.
* **Time-averaging is a red herring**: FBM window `DT_AVG=0.175 s` but
  thermal ND/TI/NE vary <1% across it.
* **`sum(BMVOL)==sum(DVOL)`** exactly (6.815e7 cm³) — not a volume-
  coverage issue.

Fix (default `--fast-norm bdens`): rescale `n_fast` per flux-surface row
by `BDENS_D(x)/<n_fast(F)>_flux(x)` (per-row factors 1.12–1.27). This
preserves F's validated poloidal/energy *shape* and fixes only the radial
*magnitude*. `--fast-norm ntot` keeps the raw F normalization (0.81) for
comparison.

## CLI / structure
```bash
python src/neutron/km14/bt_zone_integrator.py 104614 M30 --idx 1 --plot
python src/neutron/km14/bt_zone_integrator.py 104614 M30 --fast-norm ntot
```
Flags: `--idx --data-dir --nsamp --seed --fast-norm {bdens,ntot} --plot --save --no-plot`.
Inputs: `F_D_NBI`/`E_D_NBI`/`A_D_NBI`/`BMVOL`/`X2D` (`_fi`); `ND`/`NT`/`TI`/
`BDENS_D` (main CDF, interp to zone X2D); `BTN4`/`BTN1` (`_neut`, reference).
Reuses `bt_kinematics` for σ and `relative_cm_energy_keV`.

## Open items before commit 3 (LOS projection)

* Per-zone *scatter* (min 0.47, max 1.23 about the median ~1.0) is the
  residual poloidal mismatch: BDENS_D is a flux function so the per-row
  renorm matches the flux average, but individual θ-zones still differ
  from NUBEAM's poloidal F structure. Acceptable for the 4π test;
  irrelevant once we integrate along the LOS, but worth a look.
* Commit 3 still needs the **B-field reader** (`PSIRZ` 101×161 on
  `RGRID`/`ZGRID` in cm; `B_φ` from `BZXR`·`GFUN`/R, BZXR=854.8 T·cm) for
  the velocity-vector direction and the LOS projection.
* Whether to fold the same BDENS_D renorm into commit 3 (yes — it's a
  density normalization, independent of the angular treatment).

---

# `bt_los_emissivity.py` — context (added 2026-06-03, commit 3 of step 2)

`src/neutron/km14/bt_los_emissivity.py` is the direction-resolved
beam-target emissivity: it extends the 4π per-zone reactivity of
`bt_zone_integrator` to `dε/dΩ(n̂)` and evaluates the anisotropy along the
KM14 vertical sight line.

## Headline result

**For KM14's vertical sight line, BT neutron emission is isotropic to
<0.5%** — reactivity-weighted mean anisotropy factor g = 1.00 (DD +0.1%,
DT −0.3%) on 104614 M30 idx 1. So the commit-2 4π emissivity is an
excellent approximation for the KM14 BT contribution; **no angular
correction is needed**. BT can be folded into the KM14 LOS using the
commit-2 spatial emissivity directly (e.g. via the `los_thermal_rate`
chord-integral machinery).

Why: the CM→lab boost axis is the centre-of-mass velocity, which is
dominated by the fast ion streaming along B̂. B̂ is 97–100% toroidal, and
KM14 views vertically — `angle(B̂, LOS)` is 76–90° (median 84°), i.e.
nearly perpendicular to the boost. The kinematic anisotropy is a dipole
`1 ± 2β` (β = v_cm/u_n ≈ 0.07 for DD, 0.02 for DT) which **cancels at 90°**
to leading order; bidirectional fast ions (co+counter) cancel residuals
further. The old notes' "10–30% anisotropy" is the *forward/backward*
peak-to-peak (what a tangential viewer like KN3 sees), not KM14.

## Validation

* **Self-test** (`python bt_los_emissivity.py --test`): a clean 100 keV D
  beam fired +ŷ into cold D reproduces the analytic boost dipole — forward
  g=1.144 (predict 1.143), perp g=1.001, back g=0.869 (predict 0.857).
  Confirms the cone estimator is correct in absolute magnitude.
* **Consistency test 1** (built into the main run): the vector-based
  `Eps_4pi` reproduces NUBEAM BTN4/BTN1 (0.983/0.988) — the velocity-vector
  machinery doesn't change the scalar rate (matches commit 2).
* Estimator cross-check: projecting along toroidal ŷ (along boost) gives
  the largest g (DD +1.5%), vertical/radial ≈1 — directionally sensible.

## B-field reader (new, validated)

`BField(cdf_path, time)` reads `PSIRZ` (flat, **C-order (Z, R)** →
reshape (nZ, nR); axis = ψ-min = PSI0_TR=0 lands at the known magnetic
axis) on `RGRID`/`ZGRID` (cm). Then `B_R=-(1/R)∂ψ/∂Z`, `B_Z=+(1/R)∂ψ/∂R`
(ψ in Wb/rad → T), `B_phi=BZXR/R` (BZXR = R·Bφ vacuum = 854.8 T·cm; GFUN
diamagnetic correction ≤1.6% neglected — negligible for B̂ *direction*).
Bilinear interp returns B̂ in the (R, φ, Z) basis. Sanity-checked: B_pol→0
on axis, B_Z=+0.55 T outboard midplane, B_φ=2.82 T·(R_axis/R).

## Method

Per zone, MC-sample fast ions from F (speed + pitch ξ + uniform gyrophase,
built into a velocity *vector* using the zone B̂), pair with a Maxwellian
(+toroidal rotation Ω·R) thermal partner, weight by σ·v_rel, emit an
isotropic-CM neutron → lab direction via `bt_kinematics.neutron_lab_velocity`.
`dε/dΩ(n̂_LOS)` = `Eps_4pi` × (reaction-weight fraction in a cone of
half-angle `--cone-deg` about n̂_LOS) / cone solid angle. n̂_LOS = +ẑ.

## CLI
```bash
python src/neutron/km14/bt_los_emissivity.py 104614 M30 --idx 1 --plot
python src/neutron/km14/bt_los_emissivity.py --test
```
Flags: `--idx --data-dir --nsamp --cone-deg --fast-norm --no-rotation --seed --plot --save --no-plot`.
Plot (2×n): row 0 `Eps_4pi(R,Z)`, row 1 anisotropy factor g(R,Z).

## Still open (commit 3b / later)

* **FIDASIM benchmark** — the external cross-check from the original
  commit-3 scope. Needs FIDASIM input prep (TRANSP CDF + EFIT equilibrium
  in FIDASIM format; preprocessors in `~/jet/FIDASIM*/lib/python/`). Given
  the anisotropy is <0.5% for KM14, this is now lower priority for the
  *thermal-equivalent* goal but still wanted to validate the absolute
  angular machinery for tangential channels (KN3).
* **Full LOS chord integral** — fold the (now validated as ~isotropic) BT
  emissivity into the KM14 chord integral. Because g≈1, this reduces to
  running the commit-2 4π emissivity through the `los_thermal_rate`
  chord-integral / THKM14-weight pipeline. No new angular machinery needed.
* Per-zone g has visible MC scatter (~few % at nsamp=60k); the
  reactivity-weighted mean is robust. Raise `--nsamp` if a per-zone map is
  wanted.

---
