# `compare_los.py` — KM9 vs KM14 real-LOS comparison (added 2026-06-24)

`src/neutron/compare_los.py` overlays the **KM9 (MPRu)** and **KM14
(diamond)** real-line-of-sight TH/BT diagnostics on one figure so the two
chords can be read side by side. It is a thin driver: no new physics, no
idealised box chord — only the real `.los` cell clouds.

## Why one driver (not a flag in each script)

For the *real* LOS both `km9/los_th_bt_ratio.py` and the `--los-file`
branch of `km14/los_th_bt_ratio.py` reduce to the **same machinery**: the
per-shell detector coupling `f_det = C_bin/DVOL` from the diagnostic's own
`.los` cells (`los_common.los_file_detector_rate`), applied to the thermal
(`THNTX`) and flux-averaged beam-target (`BTN`) emissivity to form the
detector-weighted cumulative

    (TH/BT)_LOS = Σ(TH·f_det·DVOL) / Σ(BT·f_det·DVOL).

The two diagnostics differ only in the `.los` file and the (auto-detected)
chord orientation. So `compare_los.real_los_profiles()` is just
`km9/los_th_bt_ratio.main()` generalised over the `.los` file — it
reproduces both standalone scripts exactly (verified: KM9 `(TH/BT)_LOS =
0.5095`, KM14 `= 0.5353` on 104614 M29 idx 2, matching each script's own
output) without touching either of them.

## The figure (`--plot`, 1 × (2 + N); a 1×4 for KM9 + KM14)

* **(0) Detector-coupling weight + cumulative signal.** Solid (left axis):
  the normalized `f_det(rhot)` for each diagnostic
  (`normalized_coupling_weight`, `∫f_det·DVOL dρ = 1`, units 1/m³ — the raw
  etendue scale differs between KM9 and KM14, so only the *normalized* shape
  is comparable). Dotted (right axis, 0→1): the **cumulative fraction of the
  detector-weighted BT signal** `Σ BT·f_det·DVOL` banked by each rhot — the
  "ballast" curve that explains panel 1 (see below).
* **(1) TH/BT vs rhot.** Each diagnostic's cumulative real-LOS C-weighted
  TH/BT, over a *shared* grey full-plasma (0D) cumulative and the shared
  local `TH(ρ)/BT(ρ)`. The full-plasma cumulative and the local ratio depend
  only on the emissivity profiles + DVOL (no LOS), so they are identical for
  both diagnostics at the same pulse/run/idx/channel and are drawn once.
* **(2..) Per-shell detector signal — one panel per diagnostic.**
  `TH·f_det·DVOL` (solid) / `BT·f_det·DVOL` (dashed), titled with the
  integrated detector rates. Each diagnostic gets its **own y-axis**: the
  absolute detector-reaching rate scales with the detector etendue, which
  differs by ~the solid angle between KM9 (~8.5e9 n/s TH) and KM14
  (~4.3e8 n/s TH), so a shared axis buries the smaller signal.

## Why the cumulative TH/BT diverges past rhot ≈ 0.4 although the weights overlap

This is the subtle point the comparison surfaces. On 104614 M29 idx 2 the
two cumulative ratios sit on top of each other out to rhot ≈ 0.4 (gap
`R_KM14 − R_KM9 ≈ +0.004`) and then fan apart from 0.4 → 0.8 (gap → +0.025,
then plateau) — precisely the band where the *normalized* `f_det` weights
(panel 0, solid) nearly coincide. The mechanism:

**The cumulative ratio is a weighted average of the *local* ratio.**

    R_cum(ρ) = Σ₀^ρ TH·f·DVOL / Σ₀^ρ BT·f·DVOL = ⟨ r_local ⟩_w,
    r_local = TH(ρ)/BT(ρ),   weight w = BT · f_det · DVOL.

`r_local` is **identical** for the two diagnostics (same plasma); only the
weight `w` differs, through `f_det`. Two consequences combine:

1. **Where `r_local` is flat, the weighting cannot matter.** `r_local` stays
   in ~0.62–0.74 across the whole core (ρ ≲ 0.4), then collapses
   (0.53 → 0.28 → 0.15 → 0.08 by ρ 0.7). A weighted average of a
   near-constant function is that constant *regardless of the weights*, so
   even though the core weights differ a lot (KM14's normalized `f_det` runs
   +25–43 % above KM9's around ρ 0.15–0.35) both cumulative ratios sit at
   ~0.63–0.67 out to 0.4. The weight difference is present in the core but
   **invisible** — there is no radial contrast in `r_local` to expose it.

2. **The steep edge falloff of `r_local` is the lever, and the two LOS reach
   ρ 0.4 with different "ballast".** The cumulative average is dragged toward
   the low edge value in proportion to the edge weight *relative to the core
   weight already accumulated*. By ρ = 0.4:
   - **KM14** (vertical, core-peaked) has banked **69 %** of its BT signal →
     only 31 % sits in the low-ratio edge;
   - **KM9** (horizontal) only **63 %** → **37 %** in the edge.

   So KM9 carries more weight into the collapsing-ratio edge and its running
   average is pulled down further (→ 0.508); KM14, with more high-ratio core
   ballast, resists the pull (→ 0.535). In the far edge the normalized weights
   even cross (KM9 ~10–20 % *above* KM14 beyond ρ 0.8), reinforcing the same
   direction. Note this would happen **even if the edge weights were
   identical**: identical edge increments dilute a smaller accumulated
   denominator (KM9) more than a larger one (KM14).

**One-line summary.** The divergence is *seeded in the core* (ρ < 0.4),
where the two LOS genuinely weight the plasma differently, but stays *hidden*
there because `r_local` is flat; it only *becomes visible past 0.4*, where
the steep edge falloff of `r_local` levers the difference in how much signal
each LOS banked in the core. The normalized weights overlapping at ρ > 0.4 is
exactly why no *new* divergence is added there — the gap that appears is the
core difference finally being exposed, not created. The panel-0 dotted
cumulative-fraction curves (69 % vs 63 % at ρ 0.4) make this explicit.

This is the cross-diagnostic counterpart of the single-chord
"cumulative-ratio plateau / outer-slab drop" discussion in
[`th_bt_ratio.md`](th_bt_ratio.md): there one weighting is flat on the
enclosed core so DVOL / `f` / `f_det` give the same cumulative ratio until
the outer slab; here two *different* chords give the same cumulative ratio on
the flat-`r_local` core and separate only on the structured edge.

## CLI

```bash
python src/neutron/compare_los.py 104614 M29 --idx 2 --plot
python src/neutron/compare_los.py 104614 M29 --idx 2 --channel dd --save cmp.png
python src/neutron/compare_los.py 104614 M29 --idx 2 \
    --km9-los src/neutron/km9/KM9.los --km14-los src/neutron/km14/KM3.los
```

Flags: `--idx --channel {total,dd,dt} --data-dir --km9-los --km14-los --save
--no-plot`. Defaults: `km9/KM9.los`, `km14/KM3.los`. `--save` writes the PNG
headless; `--no-plot` text-only (prints the per-diagnostic `(TH/BT)_LOS`,
`rho_50` and the shared full-plasma 0D ratio).
