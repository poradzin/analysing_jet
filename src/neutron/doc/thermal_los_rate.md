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

## `--los-file` mode — real KM14 LOS cell file (added 2026-06-17)

M. Nocente supplied the **exact** KM14 line-of-sight geometry as a cell
cloud: `src/neutron/km14/_KM3.los` (117 424 rows), described in
`KM3_LoS_readme.txt`. Each row is one LOS cell with 8 columns
`x, y, z, C, V, u, v, w` in a fixed (discharge-/time-independent)
right-handed Cartesian frame — poloidal plane = x-z, toroidal plane = x-y:

* `x` tangential/toroidal offset [m] (|x|≤0.098); `y` horizontal in-plane
  [m]; `z` vertical Z [m]. Cell major radius `R = sqrt(x²+y²)` (x≪y so
  R≈y to <2 mm). Ranges: R∈[2.60, 3.16], Z∈[−2, 2], ΔZ≈0.020 m uniform,
  footprint tapers 654→514 cells/slice (real diverging collimated view).
* `V` cell volume [m³] (~2.4e-6); `C` **etendue weight [m³]** with
  detector solid angle `Ω = 4π·C/V`, so for the isotropic thermal source
  the rate of neutrons *reaching the detector* from a cell is
  `ε·V·Ω/4π = ε·C`. C spans 3+ decades.
* `u, v, w` unit vector of the emission direction a neutron must have to
  reach the detector (`w≈0.9997–1.0`, near-vertical). Read & returned but
  **unused for the isotropic thermal channel**; it is exactly the per-cell
  `n̂_LOS` the BT angular code (step-2 commit 3) needs.

**Box default updated** to the real footprint: `R_MIN/MAX = 2.60/3.16`
(was 2.70/3.10). The idealised-box pipeline (Rate_LOS/Rate_tor/`rho_bnd`/
geometric `f`) is unchanged in structure and still runs — it is the
comparison baseline. The old `rho_bnd=0.3070` reproduces exactly under
`--Rmin 2.70 --Rmax 3.10` (the scattered-point refactor below is faithful).

**New `--los-file PATH`** adds a detector-weighted block (`read_los_file`
+ `los_file_detector_rate`), computing on the *real* cells:
* `Rate_chord = Σ ε·V` (real narrow chord, no solid angle) and
  `Rate_det = Σ ε·C` — **the thermal rate actually reaching KM14**.
* C-binned shell response → detector-coupling weight
  `f_det(rhot) = C_bin/DVOL` (dimensionless ~Ω/4π × swept fraction) and
  `THKM14_det = THNTX·f_det`, with closure `Σ THNTX·DVOL·f_det = Rate_det`.
  **Anti-aliased binning (added 2026-06-17, default).** The cells are
  Cartesian, flux surfaces ~cylindrical, so one cell spans several rhot bins
  (mostly in Z) → point binning makes `f_det` noisy (same artifact as the box
  path, CONTEXT lines ~403-416). Each cell is given a rhot half-span
  `0.5(|∂rhot/∂R|·dR + |∂rhot/∂Z|·dZ)` and its C spread over that range via the
  shared `_subgrid_bin` helper (now also used by `los_shell_fraction`). Local
  `|∇rhot|` from a dense structured field (`RegularGridInterpolator`); per-cell
  `dZ` = z-slice spacing (≈0.0201 m, uniform), `dR = √(V/dZ)` (a cell is ~one
  slice thick, ≈11 mm). The cells are *not* a regular poloidal lattice — each
  carries its own volume `V≈2.4e-6 m³` (≈13 mm). Result: bin-to-bin
  `|Δf_det|/f_det` median drops ~14× (0.24→0.017), `Rate_det` unchanged,
  closure improves to 3.5e-5. `--no-subgrid` reverts to point binning.

  **Near-axis de-artifacting (added 2026-06-17).** Two further fixes so
  `f_det` is a clean flat plateau over the fully-enclosed core (matching the
  user's expectation that THKM14_det ∝ THNTX where whole surfaces are captured):
  1. *DVOL-weighted (Jacobian) splat.* `_subgrid_bin` gained an optional
     `density=` arg; passing `density=DVOL` distributes each Cartesian cell's
     weight across the shells it spans ∝ shell volume (the physically correct
     split — near axis `dV/dρ→0` so a cell holds more volume in its outer-rhot
     part) instead of uniform-in-rhot. `density=None` reproduces the old
     overlap/span formula *exactly* (geometric path untouched, 0.3070 preserved).
     This removes the over-fill(centre)/starve(next ring) oscillation.
  2. *Enclosed-shell flattening.* `rhot_crit` (innermost surface reaching an
     R-boundary of the footprint, same formula as the geometric f=1 pin) is
     computed from the structured rhot grid; bins with `edges<=rhot_crit` are
     pinned to the DVOL-weighted mean `ΣC_enc/ΣDVOL_enc`. `f_det` is genuinely
     flat there (verified 0.94–1.04× plateau over [0.018, rhot_crit]); the
     volume weighting down-weights the under-sampled tiny-DVOL inner bins (the
     ~13 mm cells barely resolve the axis tube), and `ΣDVOL·f_det` over the
     region is unchanged so `Rate_det` closure holds. The detector analogue of
     pinning f=1, but to the real plateau (`f_det≈6.9e-10`, not 1).
  Net: clean flat `f_det` from rhot=0 to `rhot_crit≈0.124`, then smooth roll-off
  (104614 M29); `Rate_det`/`rho_50` unchanged, closure 2.4e-4. `--plot` marks
  `rhot_crit` (cyan) and `rho_50` (black) on the weight panel.

  **Post-`rhot_crit` bump is physical (not an artifact).** `f_det = capture ×
  ⟨Ω/4π⟩` where capture = `V_bin/DVOL` (geometric shell fraction) and ⟨Ω/4π⟩ =
  `C_bin/V_bin` (mean detector solid angle from `C`). Just past `rhot_crit`
  capture *decreases* monotonically (0.83%→0.47% by rho 0.30, as expected — the
  chord stops enclosing the whole surface) but ⟨Ω/4π⟩ *rises* ~+45% to a peak
  near rho≈0.2 (real KM14 collimator etendue), so their product bumps ~+4.6%
  over rho 0.12–0.16 before the capture fall dominates and `f_det` decreases.
  The residual ±1–2% wiggle for rho>`rhot_crit` is finite-cell sampling noise
  (117k cells / ~200 shells), not aliasing. Both confirmed by decomposing
  `f_det` into the two factors. A cosmetic **3-bin running median** is applied
  to `f_det` *outside* `rhot_crit` (`_running_median3`; enclosed plateau and the
  physical bump/roll-off preserved, `Rate_det` unchanged as it is summed
  directly from cells) to suppress the wiggle.
* Signal-median radius `rho_50` (50 % of `Rate_det` enclosed) as the
  detector-weighted analogue of `rho_bnd`.

**Refactor enabling it:** `EqCDF`/`EqPPF` gained `_scatter_nodes()` +
`rhot_at_points(Rq,Zq)` (shared module helper `_rhot_scatter`), so rhot is
evaluated at the irregular cell cloud with the *same* pinned-axis(ψn=0)+
pinned-LCFS(ψn=1) griddata as the box; `rhot_on_grid` now just reshapes a
`rhot_at_points` call (behaviour identical — 0.3070 cross-check).

**First result (104614 M29, t=53.527 s, CDF eq):** 66.1 % of cells inside
LCFS; `Rate_chord=1.86e15`, `Rate_det=1.82e8 n/s` (closure 1.5e-4);
`rho_50=0.252` vs idealised-box `rho_bnd=0.3585` (wider chord) — the
solid-angle weighting pulls the effective radius **inward** (Ω largest
looking straight up through the core), as expected. Runs <1 s.
`--save` adds `f_det, THKM14_det` columns; `--plot` overlays the
normalised `f_det` and `rho_50` on the weight/profile panels, and panel
(1,0) scatters the real LOS cells (coloured by `log10 C` = etendue) over the
simplified box rectangle. That panel shows the real LOS is a **slanted,
narrow band** (tilts outboard with height; C peaks in a central strip and
tapers at the footprint edges), not the vertical box — which is why `rho_50`
is more core-concentrated than the box `rho_bnd`.

## `los_th_bt_ratio.py` `--los-file` — exact BT emission directions (added 2026-06-17)

`los_th_bt_ratio.py --bt-aniso` weights each NUBEAM zone's BT rate by the
directional factor `g = (dEps/dΩ toward detector)/(Eps/4π)`, MC-sampled along a
per-zone LOS direction. That direction was a **point-detector approximation**
(zone → a point at `--detector-R`, `--detector-height` above Zmag, zero toroidal
component). New **`--los-file _KM3.los`** replaces it with the file's **exact
per-cell versor `(u,v,w)`**: detector-Cartesian (x=toroidal, y=radial, z=vert)
→ `(R,φ,Z)` MC basis as **`(v, u, w)`** (restores the small toroidal tilt the
point model zeroed). `los_directions_for_zones` griddata-interpolates the versor
field onto each zone `(Rz,Zz)`; zones outside the LOS footprint keep the
point-detector estimate. `bt_anisotropy_factors` gained the override + a
diagnostic (zones-in-footprint, angular diff vs point model); `tilt`/`ang_bL`
now computed from the final direction.

**104614 M29 idx2 dd, --bt-aniso:** 79/220 zones in footprint; real dirs differ
from point-detector by med/max **0.6/1.4°**; `<g>_BT` 1.0044 → 1.0042,
`(TH/BT)_LOS` shift −0.43% → −0.42% — negligible (KM14 BT ≈ isotropic, LOS ≈
vertical, as commit 3 found), but the input is now exact. Default (no
`--bt-aniso`/`--los-file`) path unchanged: `(TH/BT)_LOS = 0.3708`.

**`--los-file` also adds real-LOS plots + C-weighted cumulative TH/BT (without
needing `--bt-aniso`).** `main()` computes the per-shell detector coupling
`f_det = C_bin/DVOL` by reusing `los_thermal_rate.los_file_detector_rate` (build
an `EqCDF(cdf_path, fi['time']+40)`; `f_det` is emissivity-independent so the TH
profile is passed just to satisfy the signature), then applies it to *both* TH
and BT (BT≈isotropic for KM14 → same etendue coupling) for a real-detector
cumulative `Σ(TH·f_det·DVOL)/Σ(BT·f_det·DVOL)`. `make_plot(..., los_det=)`:
(0,2) overlays the real LOS cells (lime) on the TH/BT-local map; (1,0) overlays
normalised `f_det` + `rhot_crit`; (1,3) overlays the C-weighted cumulative TH/BT
(magenta) beside the geometric-LOS (cyan) and full-plasma (grey). Local TH/BT
(1,2) is weight-independent (= TH/BT, cancels) so unchanged. 104614 M29 total:
C-weighted cum (TH/BT)_LOS=0.540 vs geometric-LOS 0.535 vs full-plasma 0.439 —
the etendue weights the core (higher TH/BT) slightly more. Default path
(no `--los-file`) unchanged.

Remaining idea (bigger, not done): a full per-cell BT directional LOS integral
`Σ eps_BT(rhot, n̂=(u,v,w))·C` mirroring the thermal `Σε·C`, instead of the
per-zone flux-averaged `g`. Only worth it if higher-fidelity BT anisotropy is
needed; current `g`-based treatment is adequate given the <0.5% effect.

---
