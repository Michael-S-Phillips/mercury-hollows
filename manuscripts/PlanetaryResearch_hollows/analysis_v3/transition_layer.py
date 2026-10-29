"""Thickness of the layer that crosses the orthorhombic-to-monoclinic sulfur transition (368.5 K)
every Mercury solar day, and whether it predicts hollow hosting better than dT.

Usage: python3 transition_layer.py <dir with thermal.py> <craters_dT.csv> <out_dir>
thermal.py physics is reused unchanged (constant conductivity, TI 47, Bond albedo 0.09); the column is
refined to DSKIN/10 so the layer thickness can be interpolated between nodes.
L368 = depth at which the daily maximum temperature falls to 368.5 K (the daily minimum is below 368.5 K
at every depth, because the deep temperature is <= 329 K), so the whole layer above L368 cycles daily.
"""
import sys, json, numpy as np, pandas as pd, statsmodels.formula.api as smf
TH, CR, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, TH); import thermal as th

def layer(lat, lon, thr=(368.5, 388.4), nz=140, nt=900, ndays=45):
    t, F = th._flux(lat, lon, nt)
    dt, dz = t[1]-t[0], th.DSKIN/10
    assert dt < 0.45*dz*dz/th.KAPPA
    c = th.KAPPA*dt/dz**2; ns = lat.size
    T = np.full((ns, nz), 300.0)
    for day in range(ndays):
        net = np.zeros(ns); Tx = np.full_like(T, -1.0); Tn = np.full_like(T, 1e4)
        for n in range(nt):
            T[:, 1:-1] += c*(T[:, 2:]-2*T[:, 1:-1]+T[:, :-2]); T[:, -1] = T[:, -2]
            Ts, T1 = T[:, 0].copy(), T[:, 1]
            for _ in range(8):
                f = th.EMIS*th.SIG*Ts**4 - F[:, n] + th.K*(Ts-T1)/dz
                Ts -= f/(4*th.EMIS*th.SIG*Ts**3 + th.K/dz)
            T[:, 0] = Ts; Tx = np.maximum(Tx, T); Tn = np.minimum(Tn, T)
            net += F[:, n] - th.EMIS*th.SIG*Ts**4
        net /= nt
        T += np.clip(net/(4*th.EMIS*th.SIG*np.maximum(T[:, -1], 50)**3), -40, 40)[:, None]
    z = np.arange(nz)*dz
    out = {}
    for thr_ in thr:
        L = np.zeros(ns)
        for k in range(ns):
            i = np.where(Tx[k] > thr_)[0]
            if i.size == 0: continue
            j = i.max()
            if j+1 < nz:   # linear interpolation of the crossing between nodes j and j+1
                L[k] = z[j] + dz*(Tx[k, j]-thr_)/(Tx[k, j]-Tx[k, j+1])
            else: L[k] = z[j]
        out[thr_] = L
    assert (Tn[:, 1:] < 368.5).all()
    return out, Tx[:, 0], T[:, -1], float(np.abs(net).max())

la = np.arange(0, 91, 5.0); lo = np.arange(0, 91, 5.0)
LA, LO = np.meshgrid(la, lo, indexing="ij")
Ls, Tsmax, Tdeep, res = layer(LA.ravel(), LO.ravel())
L368 = Ls[368.5].reshape(LA.shape); L388 = Ls[388.4].reshape(LA.shape)
np.savez(f"{OUT}/transition_layer_grid.npz", la=la, lo=lo, L368=L368, L388=L388,
         Tsmax=Tsmax.reshape(LA.shape), Tdeep=Tdeep.reshape(LA.shape))

s = pd.read_csv(CR)
s["L368"] = 100*th.lookup(la, lo, L368, s.latitude.values, s.longitude.values)   # cm
s["L388"] = 100*th.lookup(la, lo, L388, s.latitude.values, s.longitude.values)
R = dict(threshold_K=368.5, melt_K=388.4, skin_depth_m=float(th.DSKIN), dz_m=float(th.DSKIN/10),
         max_flux_residual=res, n_craters=int(len(s)), n_hosts=int(s.hosts.sum()),
         solar_days_per_Gyr=float(1e9*365.25/th.P_SOLAR))
R["grid_cm"] = {f"lat{a}": {f"lon{o}": round(float(100*L368[list(la).index(a), list(lo).index(o)]), 1) for o in (0, 30, 60, 90)}
                for a in (0, 30, 45, 60, 75, 85)}
R["corr_L368_dT"] = float(np.corrcoef(s.L368, s.dT)[0, 1])
fits = {}
for name, f in [("dT", "hosts ~ dT * logD + abslat"), ("L368", "hosts ~ L368 * logD + abslat"),
                ("L388", "hosts ~ L388 * logD + abslat"), ("both", "hosts ~ dT * logD + L368 * logD + abslat")]:
    m = smf.logit(f, data=s).fit(disp=0)
    fits[name] = dict(aic=float(m.aic), params={k: float(v) for k, v in m.params.items()},
                      p={k: float(v) for k, v in m.pvalues.items()})
R["fits_abslat_le60"] = fits
# odds ratio per cm of layer at 10 and 100 km, from the L368 model
p = fits["L368"]["params"]
R["L368_OR_per_cm"] = {"10km": float(np.exp(p["L368"]+p["L368:logD"])), "100km": float(np.exp(p["L368"]+2*p["L368:logD"]))}
# equator hot-vs-warm difference in L368 and the implied odds ratio at 10 km
dL = 100*(th.lookup(la, lo, L368, 0, 0)-th.lookup(la, lo, L368, 0, 90))
R["equator_dL_cm"] = float(dL)
json.dump(R, open(f"{OUT}/transition_layer_record.json", "w"), indent=1)
print(json.dumps({k: v for k, v in R.items() if k != "fits_abslat_le60"}, indent=1))
for k, v in fits.items(): print(k, round(v["aic"], 1), {a: round(b, 4) for a, b in v["params"].items()})
