from netCDF4 import Dataset
import numpy as np

transp_path = '/common/transp_shared/Data/result/JET'
pulse = 104614
runid = 'M30'
cdf_file = str(pulse)+runid+'.CDF'
run_path = '/'.join([transp_path, str(pulse),runid,cdf_file])

f = Dataset(run_path)
v = lambda k: np.asarray(f.variables[k][:])
for h, ts in [("H98Y2","TAUE98Y2"), ("H98Y2E","TAUE98Y2E"),
              ("H89P","TAUE89P"), ("H97LG","TAUE97LG")]:
    num = v(h) * v(ts)
    for ref in ("TAUA1", "TAUEA", "TAUEE"):
        print(h, ref, np.nanmedian(num / v(ref)))

a, e = v("TAUA1"), v("TAUEA")
print(np.allclose(a, e), np.nanmax(np.abs(a - e)))
for k in ("TAUA1", "TAUEA", "TAUEE"):
    print(k, f.variables[k].long_name, np.nanmedian(v(k)))
for k in ("UFASTPA", "UFASTPP"):   # fast-ion energy density, if present
    print(k, k in f.variables)


t = v("TIME")
dv = v("DVOL")                      # zone volumes, (time, X)
for k in ("UTHRM", "UFASTPA", "UFASTPP", "DVOL"):
    print(k, f.variables[k].units)
W = {k: np.nansum(v(k) * dv, axis=1) for k in ("UTHRM", "UFASTPA", "UFASTPP")}
Wth  = W["UTHRM"]
Wfst = W["UFASTPA"] + W["UFASTPP"]
print("Wth  median", np.nanmedian(Wth))
print("Wfst median", np.nanmedian(Wfst))
print("Wfst/Wth median", np.nanmedian(Wfst / Wth))
Wtot = Wth + Wfst
num89 = v("H89P") * v("TAUE89P")
P89 = Wtot / num89
dWdt = np.gradient(Wtot, t)
X = np.column_stack([v("PINJ") + v("POHT"), v("PICHTOT"), dWdt, np.ones_like(dWdt)])
ok = np.all(np.isfinite(X), axis=1) & np.isfinite(P89)
coef, *_ = np.linalg.lstsq(X[ok], P89[ok], rcond=None)
print("NBI+OH %.3f  ICRH %.3f  dW/dt %.3f  const %.2e" % tuple(coef))
print("median dW/dt [MW]:", np.nanmedian(dWdt) * 1e-6)

dWth = np.gradient(Wth, t)

for name, P in [("TAUEA", Wth / v("TAUEA")), ("TAUA1", Wtot / v("TAUA1"))]:
    print(name, "corr(P_implied, dW/dt) = %.3f" % np.corrcoef(P[ok], (dWth if name=="TAUEA" else dWdt)[ok])[0,1])

dWth = np.gradient(Wth, t)
heat = v("BPTI") + v("BPTE") + v("PICHTOT") + v("POHT")   # time-varying heating proxy
for name, P, dW in [("TAUEA", Wth/v("TAUEA"), dWth), ("TAUA1", Wtot/v("TAUA1"), dWdt)]:
    X = np.column_stack([heat, dW, np.ones_like(dW)])
    m = np.all(np.isfinite(X), axis=1) & np.isfinite(P)
    c, *_ = np.linalg.lstsq(X[m], P[m], rcond=None)
    print(f"{name}: dW/dt coefficient {c[1]:+.2f}")      # ≈ -1 corrected, ≈ 0 not
