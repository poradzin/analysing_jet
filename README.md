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

`CONDE`, `DIFB`, `DIFFE`, `DIFFNE`, `DIFWE`, `DEINT`, `CONDI`

The script still accepts `--signals` for the broader diffusivity list from the request.
