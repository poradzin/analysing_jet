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

# `los_thermal_rate.py` — context (added 2026-05-25)

`src/neutron/km14/los_thermal_rate.py` is the **realistic** version of the
KM14 line-integrated thermal-neutron estimate, replacing the geometric
approximations in `estimate_outer_th_fraction.py`. It builds a dense (R, Z)
grid covering the LOS, maps it to rhot via the actual EFIT PSI, and does a
proper 2-D trapezoidal integral. It then converts to a TRANSP-equivalent
flux-shell integral, finds the equivalent `rho_bnd`, and produces a
LOS-weighted emissivity profile `THKM14(rhot)`.

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
