#!/usr/bin/env python3
import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), "../../.."))
"""Overburden at the column thickness H that collapse (d = phi*H) needs, compared with published strengths.
phi: volume fraction of a volatile cement, here elemental sulfur at 2-4 wt% (density 2,070 kg/m3) in rock of grain density
2,900 kg/m3. d: hollow depth, raw median 20.8 m (Blewett et al. 2016) and corrected median 24-28 m equatorward of 60 deg.
Overburden sigma = rho*g*H with g = 3.70 m/s2 and bulk density 2,000-2,900 kg/m3.
Strengths: cemented Ries suevite (Heap et al. 2020, Icarus 349, 113873): UCS 7.3-48.6 MPa; hydrostatic pore-collapse pressure P*
22 +/- 5 MPa (porosity 0.39) and 90 +/- 5 MPa (porosity 0.24). Uncemented granular aggregates in one-dimensional compression:
yield stress 2 MPa (decomposed granite), 7 MPa (quartz sand), 10 MPa (glass beads) (Wu, Yamamoto & Izumi); yield stress scales
with grain tensile strength and rises as grains get smaller (McDowell 2002)."""
import json, numpy as np
g = 3.70; rho_S, rho_rock = 2070.0, 2900.0
def phi(w): return (w / rho_S) / (w / rho_S + (1 - w) / rho_rock)
R = {"phi_vol": {f"{w:.2f}": phi(w) for w in (0.02, 0.03, 0.04)}}
out = {}
for d in (20.8, 25.0, 28.0):
    for w in (0.02, 0.04):
        H = d / phi(w)
        out[f"d{d}_S{w}"] = dict(H_m=H, sigma_MPa=[2000 * g * H / 1e6, 2900 * g * H / 1e6])
R["required"] = out
# forward: H* = yield/(rho g), d = phi*H*
fw = {}
for name, sy in [("granular 2 MPa", 2e6), ("granular 7 MPa", 7e6), ("granular 10 MPa", 10e6),
                 ("suevite UCS 7.3 MPa", 7.3e6), ("suevite P* 22 MPa", 22e6), ("suevite P* 90 MPa", 90e6)]:
    Hs = [sy / (2900 * g), sy / (2000 * g)]
    fw[name] = dict(H_star_m=Hs, depth_m_for_phi_2_8_to_5_5pct=[phi(0.02) * Hs[0], phi(0.04) * Hs[1]])
R["forward"] = fw
json.dump(R, open(_REPO + "/manuscripts/GRL_insolation_hollows/depthctl/strength_calc_record.json", "w"), indent=1)
print("phi:", {k: round(v * 100, 2) for k, v in R["phi_vol"].items()}, "vol%")
for k, v in out.items(): print(k, "H = %.0f m" % v["H_m"], "sigma = %.1f-%.1f MPa" % tuple(v["sigma_MPa"]))
for k, v in fw.items(): print(k, "H* = %.0f-%.0f m" % tuple(v["H_star_m"]), "d = %.0f-%.0f m" % tuple(v["depth_m_for_phi_2_8_to_5_5pct"]))
