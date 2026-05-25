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
