"""Empirical check of the sun-azimuth convention. For every Bickel box in the 130 validation frames, the image is sampled
along the label sun direction (u, from the sunward box edge, 0, to the anti-sunward edge, 1; light travels toward +u, taking
SUB_SOLAR_AZIMUTH as clockwise from +sample with lines down, pointing to the sun) and the across-box mean DN profile is
normalised by the surroundings median. Inside depressions the sunward part should be dark (cast shadow) and the anti-sunward
part bright (sun-facing wall); with the opposite convention the pattern would be mirrored. Writes sun_convention.json."""
import json, numpy as np, pandas as pd, dep_common as C
from scipy.ndimage import map_coordinates
fr = pd.read_csv("registration_frames.csv", dtype={"num": str}); k = C.load_bickel(with_corners=False)
prof = []; az_all = []
for n in fr[fr.accepted].num:
    img, info = C.read_edr(C.edr_path(n)); az_all.append(info["subsol_az"]); a = np.radians(info["subsol_az"])
    du = np.array([-np.cos(a), -np.sin(a)]); dw = np.array([-du[1], du[0]]); N, L = info["samples"], info["lines"]
    for _, r in k[k.num == n].iterrows():
        x0, x1, y0, y1 = r.x0 * N, r.x1 * N, r.y0 * L, r.y1 * L; cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        cs = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]]) - [cx, cy]; cu, cw = cs @ du, cs @ dw
        U = np.linspace(cu.min(), cu.max(), 21); W = np.linspace(cw.min(), cw.max(), 15); UU, WW = np.meshgrid(U, W)
        xs, ys = cx + UU * du[0] + WW * dw[0], cy + UU * du[1] + WW * dw[1]
        if xs.min() < 0 or ys.min() < 0 or xs.max() > N - 1 or ys.max() > L - 1: continue
        A = map_coordinates(img, [ys, xs], order=1); ring = img[max(0, int(y0) - 10):int(y1) + 10, max(0, int(x0) - 10):int(x1) + 10]
        prof.append(A.mean(0) / np.median(ring))
p = np.array(prof); m = p.mean(0)
R = dict(n_boxes=len(p), n_frames=int(fr.accepted.sum()), u=np.round(np.linspace(0, 1, 21), 2).tolist(), mean_profile=np.round(m, 4).tolist(),
         sunward_half_mean=float(m[:10].mean()), antisunward_half_mean=float(m[11:].mean()),
         frac_boxes_sunward_darker=float((p[:, :10].mean(1) < p[:, 11:].mean(1)).mean()),
         frac_boxes_min_in_sunward_half=float((p.argmin(1) < 10).mean()), frac_boxes_max_in_antisunward_half=float((p.argmax(1) > 10).mean()),
         sun_azimuth_deg_percentiles={q: float(np.percentile(az_all, q)) for q in (5, 25, 50, 75, 95)})
json.dump(R, open("sun_convention.json", "w"), indent=1); print(json.dumps(R))
