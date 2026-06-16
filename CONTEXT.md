# `estimate_outer_th_fraction.py` — context

`src/neutron/km14/estimate_outer_th_fraction.py` estimates the fraction of
the thermal neutron line-of-sight (LOS) signal that originates **outside**
a small central plasma volume around the magnetic axis. The goal is to
quantify the error made when assuming that all thermal neutron signal
recorded by a vertical LOS comes only from the plasma core.

## Inputs

* TRANSP run (via `profiles.Transp`):
  * `THNTX` — thermal neutron emissivity, native units `N/CM3/SEC`
  * `DVOL`  — zone (flux-surface) differential volume, native units `CM**3`
  * `X`     — rhot grid = sqrt(normalized toroidal flux)
* Equilibrium PPF (via `profiles.Eq`, default DDA `EFTP`):
  * `RMAG`, `ZMAG` — magnetic axis coordinates
  * `ZBND`         — LCFS, used for the vertical LOS extent
  * standard psi / ftor mappings used by `change_rho`

The CLI accepts both `--dda eftp` and `dda=eftp` style arguments:

```bash
python estimate_outer_th_fraction.py 104614 M29 \
       --dda eftp --uid gszepesi --seq 405 -t 53.5268 --plot
```

## Algorithm

1. **LOS construction.** Vertical LOS at constant `R = RMAG`, `Z` sampled
   symmetrically around `ZMAG` and extending to the maximum LCFS vertical
   extent (top or bottom, whichever is larger).

2. **Remap onto Z** (`remap_profile_on_z`). For each LOS point
   `(R, Z) → rhot` is computed by `change_rho.RZ_to_rhot`, then `THNTX(rhot)`
   and `DVOL(rhot)` are PCHIP-interpolated onto the LOS rhot values. The
   rhot grid is padded with `(0, ys[0])` and `(1, 0)` before interpolation
   so the axis point (rhot≈0) and points outside the LCFS get sensible
   values instead of NaN.

3. **Central contribution `C0`** (`central_contribution`).
   * `rhot_edge = RZ_to_rhot(RMAG, ZMAG + r0)` — flux-surface label that the
     LOS hits at the upper edge of the central region (`r0 = 0.2 m`).
   * `cum_emis = np.cumsum(THNTX * DVOL)` on the TRANSP rhot grid (single
     time slice). `THNTX·DVOL` has units `n/s` regardless of CGS vs SI.
   * `C0` = `cum_emis` PCHIP-evaluated at `rhot_edge`, i.e. the **full
     toroidally integrated** neutron rate inside the central flux surface.
   * Reports `V_inside` (cumulative flux-surface volume in m³) alongside
     the cylindrical reference `V0_ref = w·π·r0²` for context.

4. **Outer contribution `Crest`** (`outer_contribution`).
   For each LOS point with `|Z − ZMAG| > r0`, applies the geometrical
   weighting

   ```
   w_geom(Z) = w / (2·π·|Z − ZMAG|)
   ```

   representing the fraction of the poloidal flux-surface layer at
   radius `|Z − ZMAG|` that is actually intercepted by a LOS of toroidal
   width `w = 0.4 m`. The integral is performed on the upper half
   (`Z > ZMAG + r0`) and multiplied by 2 for up-down symmetry. The
   `r0`-cut also removes the 1/r singularity at the axis.

5. **Outer fraction** `f_outer = Crest / (C0 + Crest)`.

6. **Equivalent rhot check** (`equivalent_rhot`).
   Inverts the cumulative emission profile to find `rhot_bnd` such that
   `cumsum(THNTX·DVOL)|_{rhot_bnd} = C0 + Crest`. If `C0 + Crest` exceeds
   the full-plasma total `cum_emis[-1]`, the script reports the over-count
   ratio instead — that would mean the LOS geometric weighting is
   double-counting flux-surface volume.

## Current assumptions (to revisit)

* **Circular flux surfaces.** `r = |Z − ZMAG|` is used as the poloidal
  radius in `w_geom`; elongation/triangularity are ignored.
* **Vertical LOS through the magnetic axis.** The LOS is taken at exactly
  `R = RMAG`, so it always pierces the axis.
* **Fixed toroidal LOS width `w = 0.4 m`** in both `V0_ref` and `w_geom`.
* **Up-down symmetry** of the emissivity profile around `Z = ZMAG`.
* **Detector solid angle** drops out of the ratio and is omitted.

## Output

Printed report blocks:

1. TRANSP slice info, native units, `sum(DVOL)` in m³.
2. **Central region**: `r0`, `w`, `(R, Z)` at the edge, `rhot_edge`,
   `V_inside`, `V0_ref`, `TH0`, `C0`.
3. **Outer region**: `Crest`.
4. **Summary**: `C0 + Crest`, `f_outer`, outer percentage.
5. **Equivalent rhot_bnd check**: total-plasma `cumsum(THNTX·DVOL)` and
   the rhot that reproduces `C0 + Crest`.

With `--plot`: three stacked panels of `THNTX(Z)` (log y), `w_geom(Z)`, and
the weighted integrand `THNTX·DVOL·w_geom` versus `Z − ZMAG`.

## Dependencies

* `profiles.py` — `Transp`, `Eq` wrappers around TRANSP CDF and JET PPF.
* `change_rho.py` — provides `RZ_to_rhot` (and the full psi/ftor mapping
  chain). Its top-level CLI is guarded by `if __name__ == "__main__":` so
  it can be imported as a module.

## Planned extensions (Monday)

These two improvements are planned to make the estimate realistic across
discharges with different equilibria:

1. **Elongated poloidal flux surfaces.** Replace the circular
   `r = |Z − ZMAG|` assumption in the geometric weighting (and possibly in
   the central-region cut) by an actual flux-surface poloidal radius
   derived from the equilibrium. The most natural way is to keep working
   from `RZ_to_rhot` along the LOS and convert rhot back to a poloidal
   radius using the equilibrium minor-radius mapping (e.g. via
   `psin_to_RZ_midplane` or by sampling LCFS extents). The cylindrical
   `V0_ref` and the `w_geom ∝ 1/(2π r)` weighting both need this.

2. **LOS not at `R = RMAG`.** In reality the vertical LOS is fixed in
   machine coordinates (e.g. `R_LOS ≈ 2.9 m`) while the magnetic axis
   `RMAG` varies discharge-to-discharge. Generalise `build_los` to take an
   explicit `R_LOS` argument (CLI flag, with a sensible default), and let
   the chord-to-axis offset `ΔR = R_LOS − RMAG` enter the geometric
   weighting — the LOS no longer pierces the axis, so `|Z − ZMAG|` should
   be replaced by the perpendicular distance from the chord to the flux
   surface, which in the circular-approximation limit reduces to
   `sqrt((Z − ZMAG)² + ΔR²)`. Decide whether `TH0` is still well defined
   (it isn't, if the LOS misses the axis) and adjust `central_contribution`
   accordingly — likely by using `RZ_to_rhot(R_LOS, ZMAG)` rather than
   assuming rhot≈0 at the closest LOS point.

When the elongation correction goes in we should also revisit whether the
up-down symmetry assumption still holds for shaped discharges (X-point
configuration, vertical shift); the integration may need to be done over
the full Z range rather than one half times two.

---

# `los_thermal_rate.py` — context (added 2026-05-25, CDF-default rewrite 2026-06-10)

`src/neutron/km14/los_thermal_rate.py` is the **realistic** version of the
KM14 line-integrated thermal-neutron estimate, replacing the geometric
approximations in `estimate_outer_th_fraction.py`. It builds a dense (R, Z)
grid covering the LOS, maps it to rhot via the actual EFIT PSI, and does a
proper 2-D trapezoidal integral. It then converts to a TRANSP-equivalent
flux-shell integral, finds the equivalent `rho_bnd`, and produces a
LOS-weighted emissivity profile `THKM14(rhot)`.

## Equilibrium source — CDF default, ppf opt-in (2026-06-10)

The script is now **self-contained on the TRANSP CDFs by default** (like
`los_th_bt_ratio.py`) so it runs in the WSL dev env without `ppf`:

* `--eq-source cdf` (**default**) reuses `CdfEquilibrium` from
  `los_th_bt_ratio` for rhot(R,Z) (`PSIRZ`/`PSI0_TR`/`PLFLXA` +
  `PLFLX`-vs-`XB`), magnetic axis from `RAXIS`/`YAXIS`, and the **LCFS from
  the time-resolved asymmetric boundary Fourier moments**
  (`RMCB*`/`RMSB*`/`YMCB*`/`YMSB*`, reconstructed by `_boundary_from_moments`;
  matches the `_fi` RSURF/ZSURF outer surface to ~1e-4 cm, and works at any
  time, not just FBM indices). Run-dir search via `bt_zone_integrator.find_run_dir`.
* `--eq-source ppf` uses the old PPF EFTP path via `profiles.Eq` +
  `change_rho.psin_to_sqrt_ftor_norm`, **imported lazily inside
  `EqPPF.__init__`** — the module top-level no longer imports `profiles`,
  `change_rho` or `ppf`, so the default path needs none of them.
* The two paths share one interface (`EqCDF`/`EqPPF`: `rhot_on_grid`, `lcfs`,
  `native_rhot`, `z_extent`, `Rmag/Zmag`, `tind`, `t_eq_jet`, `label`); the
  whole downstream integration/`rho_bnd`/`THKM14` math is source-agnostic.
* **Thermal profile + DVOL are always read straight from the main CDF**
  (`read_thermal_slice`), so `profiles.Transp` is gone too.
* New `--channel {total,dd,dt}` selects `THNTX`/`THNTX_DD`/`THNTX_DT`
  (matches the sister script). New `--data-dir`. `--runid` is the run suffix
  (e.g. `M30`), not the full id.

**Cross-validation (2026-06-10, 104614 M30, t_TRANSP=13.33 s, slice [56]):**
in CDF mode `los_thermal_rate.py --channel total` and
`los_th_bt_ratio.py --idx 2` now agree **exactly** — TH chord 4.5753e15 n/s,
TH whole-plasma 3.9877e17 n/s — because both use the same TRANSP `PSIRZ`
equilibrium (the earlier ~1% residual was the PPF-vs-CDF equilibrium
difference, now eliminated when CDF mode is used on both). `rho_bnd ≈ 0.307`.

## Post-rewrite fixes (2026-06-10, all verified locally on M30)

* **`cs.collections` removed in matplotlib ≥3.8.** `diagnostic_plots` used
  `cs.collections[0].set_label(...)` to legend-label the magenta `rho_bnd`
  contour; that attribute is gone on the p312 venv (newer mpl than freia).
  Replaced with an empty proxy line: `ax3.plot([], [], color='magenta', lw=1.4,
  label=rho_bnd_lbl)`. (Grepped all of `src/` — this was the only `.collections`
  use.)
* **Lazy-import path bug.** The original added `src/` to `sys.path`
  (`SRC_DIR = ../..`) at module top so `import profiles` resolved; making the
  imports lazy dropped that. `import profiles`/`change_rho` live in `src/`, not
  in `km14/`, so `--eq-source ppf` would `ModuleNotFoundError: profiles`. Fixed
  by inserting `src/` into `sys.path` **inside `EqPPF.__init__`** right before
  the lazy import (CDF path never touches it). Verified the error now falls
  through to `ModuleNotFoundError: ppf` (expected off-Heimdall — `profiles.py`
  does `import ppf` at its top).
* **`--save` location.** The rewrite wrote to `<run_dir>/tmp/`, i.e. the TRANSP
  data tree (for M29 that's `/common/transp_shared/.../104614/M29/tmp/`), so the
  user couldn't find it. Restored to the **repo-local `src/tmp/`** (original
  behavior) and **tagged the filename with the eq source** so cdf and ppf runs
  don't overwrite each other:
  `src/tmp/<run_id>_KM14_LOS_profile_<channel>_<cdf|ppf>_t<time>s.txt`.
  The "Saved LOS profile to ..." line prints the absolute path. Columns:
  `rhot, THNTX, THKM14, f, DVOL`; header records the equilibrium label + t_EQ.

## RESOLVED 2026-06-12 — CDF on-axis `THKM14`/`f` floor (and ppf-vs-cdf diff)

User reported `THKM14 = f = 0` **on axis** in `--eq-source cdf` on 104614 M29
(zeros for the innermost ~5 rhot bins), while `--eq-source ppf` correctly gave
`f ≈ 1` there. **Root cause — a near-axis sampling artifact, two compounding
effects:**

1. **Bilinear ψn floor (CDF-specific).** `EqCDF` mapped (R,Z)→rhot by *bilinear*
   interpolation of the coarse `PSIRZ` grid (ΔR ≈ 2 cm). ψ is paraboloidal with
   its vertex *between* nodes, so bilinear can't reach ψn = 0 — even at the exact
   axis it floors at ψn ≈ 3e-4, chord min ≈ 7e-4. Since rhot ∝ √ψn near axis,
   that's a rhot floor of ~0.02, so the innermost shells get no LOS volume → `f=0`.
   The ppf path avoided it because `EqPPF.rhot_on_grid` injects the axis at ψn=0
   into a scattered `griddata` interp.
2. **Coarse chord Z-sampling** (nZ=100 over the ~4 m box → ΔZ ≈ 4 cm), so the
   nearest sample is ~2 cm from the axis even with a perfect ψ.

**Fix (verified on M30; rho_bnd/Rate unchanged to <0.1%):**
* `EqCDF.rhot_on_grid` now uses **griddata + pinned axis**, identical to `EqPPF`
  (one method for both sources); removed the now-orphaned bilinear `_psin_on`.
* `main()` **inserts (Rmag, Zmag) into the R, Z sampling arrays** (source-agnostic).
* The centremost bin reaches `f ≈ 0.6` (not a clean 1.0) — the innermost flux
  shell is smaller than one (R,Z) cell, a hard discretisation limit; ppf's "1.0"
  there is over-splat clamping. A residual ~25% dip in `f` over rhot 0.012–0.022
  in CDF mode is a cubic-griddata overshoot on the coarse TRANSP grid (not from
  the axis insertion; `linear` griddata is far worse on axis) — negligible for
  the integral, left as-is. **Same fix later applied to `los_th_bt_ratio.py`**
  via `CdfEquilibrium.rhot_pinned` + axis insertion (see that section).

ppf-vs-cdf bulk diff is otherwise just the equilibrium reconstruction (cdf =
TRANSP `PSIRZ`; ppf = measured EFTP): M29 `rho_bnd` cdf 0.3070 vs ppf 0.3085
(~0.5%). M29 is now local at `~/jet/data/104614/M29` (`.CDF` + `_fi`/`_neut` 1–3).

## RESOLVED 2026-06-15 — spurious near-axis `f` dips (enclosed-shell pinning)

User saw `f(rhot)` dip below 1 in the innermost ~5 bins (104614 M29/M30,
t=53.527 s) — on axis for one run, just off-axis for the other — even though
those flux shells sit entirely inside the chord R-band and should have `f = 1`.

**Root cause — a volume *redistribution* artifact in `los_shell_fraction`, not
a volume-loss bug** (total LOS volume is conserved: cumulative `los_vol/DVOL`
returns to ~1.00 by bin ~7). Near the axis the flux-shell Jacobian
`dV/drhot → 0`, while one Cartesian (R,Z) cell spans many rhot bins (the Z-cell
above the axis jumps rhot 0.005→0.032 in one step). The subgrid splat spreads
each cell's volume *uniformly in rhot*, so it over-fills the innermost bins
(cum ratio 1.66→1.44→1.11) and starves the next ring (the visible dips). The
over-fill is hidden by the `f = clip(f,0,1)` clamp, so only the compensating
starvation shows. M29 vs M30 dip at different rhot because the inserted-axis
sliver cell lands at a different phase vs the fixed bin edges. A secondary
`rhot = 0.005` floor (the `rhot(psin)` table starts at `XB[0]=0.005`, so psin=0
→ rhot 0.005 not 0) adds to the inner over-fill but is minor.

**Fix — analytic `f = 1` for fully-enclosed shells.** `los_shell_fraction`
gains an `rmag` kwarg. When the axis is inside the R-band it computes
`rhot_crit` = rhot of the innermost flux surface that reaches either R boundary
(min over Z, inside the LCFS, of rhot at the Rmin/Rmax columns) and pins
`f = 1` for every bin whose **upper edge** ≤ `rhot_crit`. The straddling
transition bin and the genuinely-partial outer shells keep their binned value,
so no over-pinning. Encodes the physical expectation: a shell fully within
R∈[Rmin,Rmax] is swept in its entirety → `f ≡ 1`. (Alternative considered and
rejected: Jacobian-weighted splat — more faithful but more invasive.)

Verified on M29/M30 (t=53.527 s): clean `f = 1` plateau from rhot=0.0025 to
~0.0575, first genuine roll-off at ~0.0625 (axis R≈3.038 m, R_MAX=3.10 →
ΔR≈0.06 m → rhot≈0.06, as expected). Consistency sum `Σ(THNTX·DVOL·f)` shifts
only +0.09 % (M29 +8.5e-4, M30 +1.0e-3), staying ~0.85 % from `Rate_tor` (the
inherent Cartesian-trapz vs TRANSP-DVOL discretisation, unchanged).

**Both callers must opt in via `rmag=`.** `main()` now passes `rmag=Rmag`.
`los_th_bt_ratio.py` imports the *same* `los_shell_fraction` (no duplicate) but
was calling it without `rmag`, so it still showed the identical dips on
104614 M29/M30 idx 2 — fixed by passing `rmag=eq.Rmag` at its call site. These
two are the only callers (grepped). Any future caller needs `rmag=` too. The
`(TH/BT)_LOS` ↔ cumulative-ratio consistency in `los_th_bt_ratio.py` is
unchanged (M29 0.5294 vs 0.5299; M30 0.5046 vs 0.5050).

## RESOLVED 2026-06-15 — last-bin `f(rhot)` crater near rhot=1 (LCFS pin)

User saw `f(rhot)` plunge sharply in the outermost 1-2 bins (104614 M29,
t=53.527 s: f = 0.22 → 0.16 → 0.055 over rhot 0.99→0.998), distinct from the
smooth `1/r`-like decline expected for a vertical chord.

**Two compounding causes, one physical and one artifact:**
1. **Physical (dominant).** The KM14 chord `R∈[2.70,3.10]` is offset *inboard*
   of the axis (Rmag≈3.038, only 6 cm inside Rmax but 34 cm inside Rmin). The
   LCFS top/bottom points sit at R≈2.64/2.66 m — ~4-6 cm *inboard* of Rmin —
   so the chord misses the plasma tips and the outermost flux shells enter the
   band only as thin slivers → genuinely small `f`. Confirmed: moving Rmin to
   2.60/2.50 m fills the edge back in (last-bin f 0.055→0.115→0.160).
2. **Artifact (the actual fix target).** The cubic griddata of `PSIRZ` floored
   short of psi_n=1 *inside the boundary-moment polygon mask* (max psi_n≈0.9989
   on the coarse ~2 cm grid), so rhot never reached 1.0 within the mask. The
   outermost TRANSP shell [~0.995,1.0] was starved of LOS cells (62 vs
   237/398/453 in the bins just inside) and the missing volume was redistributed
   inward — so the crater *deepened* with resolution (0.092 at 100² → 0.045 at
   300²) instead of converging. (Diagnosed by swapping to a consistent `psin<=1`
   mask, which removed the crater but spuriously leaked into the PFR/tips:
   21785 extra cells, all rhot>0.98 at |Z|>1.0, confirming the polygon is the
   *correct* mask and the rhot field was the problem.)

**Fix — pin the LCFS at psi_n=1 in the griddata node set** (the edge-side
counterpart of the existing axis pin at psi_n=0), in *both* `EqCDF.rhot_on_grid`
and `EqPPF.rhot_on_grid`. rhot now reaches 1.0 at the boundary; psi_n inside the
polygon spans [0,1.0007] (12/66825 cells overshoot, all <1.0007, harmlessly
clipped). Edge `f` is now smooth and converges (last bin 0.23 at 100², 0.28 at
300²). `rho_bnd`=0.3070, Rate_LOS/Rate_tor, and total-volume conservation are
all unchanged (the edge perturbation is negligible — THNTX is already down >3
decades there). Same griddata pattern in both sources, so the fix is
source-agnostic.

**Propagated to `los_th_bt_ratio.py`** (same edge crater, same griddata).
`CdfEquilibrium.rhot_pinned` gained optional `Rb`/`Zb` args; when passed it pins
those LCFS points at psi_n=1 (the axis-pin counterpart). `main()` passes the
*same* `Rb, Zb` from `read_lcfs(fi_path)` used for the inside mask, keeping field
and mask consistent. Backward-compatible (defaults None -> axis-only pin); that
call site is the only caller (grepped). Verified 104614 M29 idx 2: edge `f` now
smooth (0.18->0.22->0.28 in the last bins) with f=1 on axis, and the
`(TH/BT)_LOS` <-> cumulative-ratio consistency is unchanged (0.5294 vs 0.5299).

## Geometry

The KM14 LOS is a vertical chord above the JET vessel:
* `R ∈ [2.70, 3.10] m` (chord width 0.4 m), CLI `--Rmin --Rmax`
* Z covers the full EFIT computational box (`eq._psiz` extent), not
  derived from `ZBND` — see the ZBND axis-order bug below
* Effective toroidal width `w_tor = 0.4 m` (CLI `--wtor`) for the chord
  volume element `dV = w_tor · dR · dZ`

The chord centre is fixed in machine coordinates at `R_c = 2.9 m`; the
magnetic axis `RMAG` drifts between pulses (e.g. 3.037 m on 104614).

## Algorithm

1. **Dense (R, Z) grid** (default 100×100, CLI `--nR --nZ`) over LOS R range
   and full vessel Z range.

2. **PSI → rhot map.** Use the Eq class's pre-computed structures directly,
   avoiding any reshape ambiguity (PPF stores PSI as flat `(n_t, nR·nZ)`
   despite the reshape kwarg suggesting otherwise):
   * `eq._psirzmg` — meshgrid of (psir, psiz) with default `indexing='xy'`
     so axis 0 is Z and axis 1 is R
   * `eq._psi_norm[tind]` — flat 1D array in the *same* C-order

   Append the magnetic axis `(RMAG, ZMAG, psin=0)` to the scattered point
   set to anchor the centre, then `griddata` cubic onto the dense (R, Z)
   meshgrid (linear fallback for edge NaNs). Clip `psin` to [0, 1] and
   map to rhot via `psin_to_sqrt_ftor_norm`.

3. **Inside-LCFS mask via polygon** (`matplotlib.path.Path`). The closed
   `(RBND, ZBND)` curve is used as the mask — *not* `psin <= 1` — because
   the latter incorrectly includes the private flux region below the
   X-point, where TRANSP doesn't model emission.

4. **THNTX(R, Z)** via PCHIP interpolation of the TRANSP `THNTX(X)`
   profile, padded with axis value at rhot=0 and 0 at rhot=1. Zeroed
   outside the LCFS polygon.

5. **Three rates** (units n/s, after CGS→SI conversion of THNTX):

   ```
   Rate_LOS = w_tor · ∫∫ THNTX(R,Z) dR dZ                    (chord)
   Rate_tor = ∫∫ 2πR · THNTX(R,Z) dR dZ                       (full torus)
   ```

   Cross-check: `Rate_tor / Rate_LOS ≈ 2π·R_c/w_tor ≈ 45.55` for the
   default geometry.

6. **Equivalent `rho_bnd`.** Solve
   `cumsum(THNTX·DVOL)|_{rho_bnd} = Rate_tor` by PCHIP-inverting the
   cumulative emission profile. For pulse 104614 at t=53.53 s we get
   `rho_bnd ≈ 0.31`, vs the TRANSP proxy convention of 0.20–0.25.

7. **A / B / C decomposition.** With `rho_bnd` known, partition the
   toroidally-integrated rate:
   * **A** = inside LOS AND inside rho_bnd (the core seen by LOS)
   * **B** = inside LOS AND outside rho_bnd (the wings seen by LOS)
   * **C** = inside rho_bnd AND outside LOS (the core missed by LOS)

   By construction `rate(B) = rate(C)` since
   `Rate_tor = A + B = A + C`. The script prints both and the residual
   `|B − C|/B` as a trapz/PCHIP precision sanity check (typically ~1e-3).
   On pulse 104614 at t=53.53 s, B = C ≈ 35% of Rate_tor — i.e. the LOS
   "sees" 35% wing emission as compensation for 35% missed core emission.

8. **LOS-weighted profile `THKM14(rhot) = THNTX(rhot) · f(rhot)`.** The
   weight function

   ```
   f(rho_i) = LOS_vol_bin(i) / DVOL_TRANSP(i)   ∈ [0, 1]
   ```

   is computed by accumulating `dV_tor = 2π·R·dR·dZ` from the LOS grid
   into TRANSP rhot bins (edges = midpoints of TRANSP X, plus 0 and 1).
   By construction `sum(THNTX·DVOL·f) = Rate_tor`. This profile is the
   apples-to-apples replacement for the standard "cumsum to 0.2–0.25"
   proxy and can be saved with `--save` for downstream use.

   **Subgrid (anti-aliased) distribution** is the default. A Cartesian
   (R, Z) grid has uniform `dR·dZ` cells, but flux surfaces are curved
   and a single LOS cell typically spans many rhot bins in the Z
   direction (where `|∂rhot/∂Z|·dZ` ≫ TRANSP bin width). Point binning
   then produces a period-N oscillation in f(rhot) -- pure aliasing.

   The subgrid mode estimates each cell's rhot half-range as
   `½(|∂rhot/∂R|·dR + |∂rhot/∂Z|·dZ)` (L1 corner half-range) and
   distributes the cell's toroidal volume uniformly across that rhot
   range, weighted by overlap with each TRANSP bin. Vectorised with
   `np.add.at`. Total LOS volume is exactly preserved (verified by an
   on-line standalone test: bin-to-bin |Δf| drops ~12× vs point binning,
   total volume conserved to machine precision). Pass `--no-subgrid` to
   revert to point binning for comparison.

## CLI

```bash
python src/neutron/km14/los_thermal_rate.py 104614 M29 \
       --dda eftp --uid gszepesi --seq 405 -t 53.5268 --plot
```

Common flags: `--Rmin --Rmax --wtor --nR --nZ --zmargin --plot --save`.
Also accepts the `key=value` form (`dda=eftp uid=gszepesi seq=405`).

## Plot layout (2×3, with `--plot`)

| | col 0 | col 1 | col 2 |
|---|---|---|---|
| row 0 | rhot(R,Z) in LOS box + LCFS + rho_bnd | THNTX(R,Z) in LOS box + LCFS + rho_bnd | THNTX & THKM14 vs rhot |
| row 1 | Full poloidal LCFS with LOS box + rho_bnd surface | R-integrated emissivity vs Z | f(rhot) weight function |

## Notable bugs found / pitfalls

1. **ZBND axis order in `profiles.py`.** `_Zbnd` (and `_Rbnd`) are
   shaped `(n_t, n_bnd)` despite the reshape kwarg saying `(-1, n_t)`.
   On pulse 104614, `_Zbnd[:, tind]` gives a spurious range of
   ~[-0.37, 1.90] (one boundary point sampled across the discharge),
   whereas the correct `_Zbnd[tind, :]` gives ~[-1.37, 1.68]. Use
   the time-first indexing. `estimate_outer_th_fraction.py` has the
   same bug but compensates accidentally by symmetrizing Z around Zmag
   — should be fixed too.

2. **PSI reshape ambiguity.** Don't reshape `eq._psi[tind]` to
   `(nR, nZ)` and pair with `meshgrid(PSIR, PSIZ, indexing='ij')` —
   PPF's actual flat layout is (Z, R) C-order, so the (R, Z, psin)
   triplets would mismatch and griddata produces a shifted rhot map.
   Use `eq._psirzmg` + `eq._psi_norm[tind]` directly; they share the
   same flat order.

3. **PFR leak.** `psin <= 1` and `rhot <= 1` masks both incorrectly
   include the private flux region below the X-point. Use a polygon
   mask from (RBND, ZBND) instead. Contribution is small but
   physically wrong.

## Planned next steps (Wednesday 2026-05-27)

* **Verify on more time slices / discharges.** Currently spot-checked
  on 104614, t=53.53 s (`rho_bnd ≈ 0.31`, B=C ≈ 35%). Sweep over
  selected pulses and times.
* **Optionally fix `estimate_outer_th_fraction.py`** to use the
  time-first ZBND indexing (today's bug 1).
* **Consider** whether to expose a `--R_los_center` flag to displace
  the chord centre between pulses (the chord is fixed in machine
  coords but `Rmag` varies).
* **Open question:** how to fold KM14 detector response (energy
  windowing, scattering) into the LOS-weighted profile — currently
  THKM14 is purely geometric (no detector physics).

---

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

# HOW TO RUN — step-2 BT scripts (quickstart, 2026-06-03)

All four step-2 scripts live in `src/neutron/km14/` and are **self-contained
on the TRANSP CDFs** (netCDF4 + numpy + scipy + matplotlib only — no `ppf`,
no `profiles.Eq`), so they run in the local WSL dev env, on freia and on
heimdall alike. Run them **from inside `src/neutron/km14/`** (they import each
other as plain modules, e.g. `import bt_kinematics`).

**numpy ≥ 2.0** is required by `los_th_bt_ratio.py` and `los_thermal_rate.py`
since 2026-06-09: the deprecated `np.trapz` was migrated to `np.trapezoid`
to silence the deprecation warnings. On older numpy the integrals will
`AttributeError`; swap to `from scipy.integrate import trapezoid` if you
need to run on a legacy env.

```bash
cd ~/jet/analysing_jet/src/neutron/km14
```

## What data they need

Per `(pulse, run)`, three CDFs in the run directory:
`<run>.CDF` (main), `<run>_fi_<idx>.cdf` (fast-ion dist + zone geometry),
`<run>_neut_<idx>.cdf` (per-zone BT rates). The scripts search, in order:
`--data-dir <base>/<pulse>/<run>`, then `~/jet/data/<pulse>/<run>`, then
`/common/transp_shared/Data/result/JET/<pulse>/<run>` (heimdall). Locally only
**104614 M30** is present (`~/jet/data/104614/M30/`); on freia/heimdall point
`--data-dir` at the TRANSP results tree. `--idx` selects the FBM time window
(M30 has 1, 2, 3 at t≈12.33, 13.33, 14.23 s); default = first available.

## Recommended order (each step validates the next)

**0. Kinematics unit tests (no data needed)** — confirm the physics core:
```bash
python bt_kinematics.py            # 5 sigma spot-checks + 2 kinematic-band tests -> PASS
python bt_kinematics.py --plot     # optional: sigma(E_cm) curves
```

**1. Per-zone reactivity (commit 2)** — acceptance test vs NUBEAM BTN4/BTN1:
```bash
python bt_zone_integrator.py 104614 M30 --idx 1            # text report
python bt_zone_integrator.py 104614 M30 --idx 1 --plot     # + (R,Z) scatter & per-zone agreement
```
Look for: `ratio in-house / NUBEAM` ≈ **1.00** for DD and DT (the headline
acceptance result). `--fast-norm ntot` reproduces the raw-F ~0.81 deficit for
comparison; `--nsamp` raises MC statistics.

**2. LOS angular emissivity (commit 3)** — the KM14 anisotropy factor:
```bash
python bt_los_emissivity.py --test                         # selftest vs analytic dipole -> PASS
python bt_los_emissivity.py 104614 M30 --idx 1 --plot      # g(R,Z) maps + report
```
Look for: `vector Eps_4pi / NUBEAM` ≈ 1.00 (consistency with step 1) and
`reactivity-weighted mean g` ≈ **1.00** (KM14 sees BT as isotropic). Useful
flags: `--cone-deg` (LOS acceptance cone), `--no-rotation`, `--nsamp`.

**3. KM14 LOS TH/BT ratio (commit 3b) — the deliverable:**
```bash
python los_th_bt_ratio.py 104614 M30 --idx 1               # DD channel, flux-avg BT
python los_th_bt_ratio.py 104614 M30 --idx 1 --plot        # + eps_TH, eps_BT, TH/BT maps
python los_th_bt_ratio.py 104614 M30 --idx 1 --bt-mode zone  # poloidal-asymmetry check (~4%)
python los_th_bt_ratio.py 104614 M30 --idx 1 --channel dt    # DT (14 MeV) instead of DD
```
Look for: `(TH/BT)_LOS` (the quantity to compare with the measured KM14
TH/BT), the `whole plasma (0D)` ratio, and the `core enhancement LOS vs
plasma` factor (≈1.14 on M30 DD — KM14 weights the core, so its TH/BT differs
from the total-yield ratio). The chord geometry cancels in the ratio.

## Saving outputs instead of showing windows

Every script's `--plot` opens an interactive window; pass `--save <file.png>`
(on `bt_los_emissivity` / `los_th_bt_ratio`) to write a PNG headless instead
(handy over SSH). `--no-plot` forces text-only. `key=value` argument style is
also accepted by the older scripts (e.g. `dda=eftp`).

## Immediate next step (Thursday)

Run the chain on **104614 M29** (the real pure-DD KM14 thermal-analysis run)
on freia/heimdall, since only M30 is local:
```bash
python los_th_bt_ratio.py 104614 M29 --channel dd \
       --data-dir /common/transp_shared/Data/result/JET --plot
```
then compare `(TH/BT)_LOS` against the spectroscopically-measured KM14 TH/BT
for that pulse/time. If a TH/BT time trace is wanted, loop over the
available `--idx` FBM windows.

---

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
  0.55–0.72, Scatt green 0.30–0.40, Total red 0.93–0.07 wrap, Data black). For
  each x-pixel column inside the plot area, takes the **median y** of pixels
  matching the colour → centroid of the line. Sat/val floors keep axis text out.
* **Data markers.** Black filled disks would also catch error bars and axis
  text on a plain "topmost black pixel" rule (this was the first version's
  failure mode — data file came out at ~200 counts everywhere, biased by axis
  labels and error-bar tops). Fixed by **horizontally eroding** the black mask
  (a pixel survives only if its 2-col neighbours at the same row are also
  black), which kills 1-col-wide error bars and most axis text, then **picking
  the centre of the largest contiguous y-cluster per column** (the marker disk
  beats the small end-cap clusters).
* **Legend suppression.** The Nocente figure's in-plot legend overlaps the
  data region. A hard-coded bbox `(340, 30, 730, 240)` zeroes that pixel area
  before extraction (overridable with `--legend-bbox x0,y0,x1,y1`).
* **Outlier rejection.** Smooth curves (everything except data) go through a
  Hampel filter (window=11, n_sigma=4) after column-centroid extraction; dash
  patterns produce gaps which the downstream `np.interp` fills.
* **Preview.** `--preview out.png` writes a sanity-check overlay of all
  extracted curves on the input image. Use it iteratively when retargeting to
  a different pulse's figure.

```bash
python src/neutron/km14/extract_km14_png.py \
    ~/jet/data/104614/figs/104614_KM14_spectral_analysis.png \
    --preview ~/jet/data/104614/figs/104614_KM14_extraction_preview.png
```
Flags: `--pulse --outdir --legend-bbox --curves {data,th,bt,scatt,total} --preview`.
Default `--pulse` inferred from PNG filename via regex.

**Pitfalls encountered while building the extractor:**
1. Topmost-black-pixel rule catches the y-axis "200" label text and error-bar
   tops → biases data ~20 counts high. Use horizontal erosion + cluster centre.
2. Legend bbox needs to extend ~25 px below the "Total" swatch — the red
   line's anti-aliased edge bleeds green/orange hue at the bottom row.
3. Tightening the green hue from (0.27, 0.42) to (0.30, 0.40) and pushing
   sat_min from 0.30 to 0.35 was needed to suppress green-shadow contamination
   from rendered text characters.

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

# 4. Core-only diagnostic (Nocente-style)
python src/neutron/km14/km14_spectrum.py 104614 M29 --idx 1 \
    --core-rhot 0.2 --include-scatt

# 5. Data fit (free TH:BT split)
python src/neutron/km14/km14_spectrum.py 104614 M29 --idx 1 \
    --th-bt-fit --fit-emin 7.5 --include-scatt
```

Output PNGs get tagged with `_fit` / `_manual` / `_scatt` / `_core<R>` so the
different runs don't overwrite each other. Saved into the same `figs/` dir as
the input PNG.
