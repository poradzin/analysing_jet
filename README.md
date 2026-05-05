# analysing_jet
Scripts to analyse JET data

## Power balance plotting

`src/plot_power_balance.py` plots TRANSP ion power balance terms as volume-integrated balance checks and per-zone contributions using `DVOL`.

Signed IEBAL convention used in the script:

`+PBTH -GAINI -PCOND +QIE -P0NET -PCONV +QROT +IHEAT +TIBAL`

Relevant source:

- TRANSP RPLOT IEBAL package definition with signs: https://w3.pppl.gov/~xshare/Rplot/nstx_multi.pdf
- TRANSP output variable table: https://transp.pppl.gov/nml/transp_output.html

## Diffusivity plotting

`src/plot_diffusivity.py` plots TRANSP diffusivity coefficients as radial profiles and time traces.

The default signal list is the core set:

`CONDE`, `CONDI`, `DIFFE`, `DIFFI`, `DIFWE`

`DEINT` is excluded from the default list because it was zero in the pulses being checked.

The script still accepts `--signals` for the broader diffusivity list from the request, and it also adds a residual profile panel using `RESPROFPE`, `RESPROFPI`, `RESPROFTE`, and `RESPROFTI`.

## Particle flux plotting

`src/plot_particle_fluxes.py` plots TRANSP particle-flux divergence terms as radial profiles and time traces.

The default flux set is:

`DIVFE`, `DIVFI`, `DFIMP`, `DIVFD`

Residual checks use:

`RESPROFPE`, `RESPROFPI`, `RESPROFPX`

and the scalar L2 residuals:

`RESL2PE`, `RESL2PI`, `RESL2PX`
