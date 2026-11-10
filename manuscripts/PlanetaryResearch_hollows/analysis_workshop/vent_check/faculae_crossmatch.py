#!/usr/bin/env python3
"""Cross-match catalog vents (swapped reading) with the 29 named faculae on Mercury (IAU names, center coordinates as
listed on Wikipedia 'List of geological features on Mercury', accessed 2 Nov 2026; longitudes there are 0-360 W, converted
to E). Also the as-labeled reading for comparison. Writes faculae_crossmatch.csv/.json."""
import os, json, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.expanduser("~/mnt/hollows-emissivity")
F = [("Abeeso",21.70,214.60),("Agwo",22.40,213.70),("Ahas",-21.37,89.40),("Amaru",-49.80,349.50),("Bibilava",16.40,202.80),
     ("Bitin",-51.55,28.45),("Coatl",-29.75,216.55),("Ejo",14.50,200.50),("Gata",-52.90,321.40),("Havu",-52.22,28.45),
     ("Ibab",14.50,199.20),("Inyoka",14.00,197.70),("Maciji",14.90,196.00),("Madu",-21.24,74.87),("Nakahi",-52.70,342.20),
     ("Nathair",36.00,295.50),("Neidr",35.90,302.70),("Ngu",-22.77,90.85),("Nzoka",15.40,194.70),("Orm",26.58,59.68),
     ("Pampu",-57.76,31.79),("Sarpa",-53.07,30.87),("Serp",3.90,309.06),("Slang",24.50,179.30),("Suge",26.10,300.40),
     ("Thueban",48.70,200.50),("Ular",-55.10,29.95),("Yinshe",-46.32,191.22),("Zmija",-37.35,267.75)]
f = pd.DataFrame(F, columns=["name", "lat", "lon_W"]); f["lon_E"] = (360 - f.lon_W) % 360
v = pd.read_csv(f"{ROOT}/data/data/MER_Vent_Locations_Classifications.csv"); v.columns = [c.strip() for c in v.columns]
R = 2439.4
def gc(la1, lo1, la2, lo2):
    la1, lo1, la2, lo2 = map(np.radians, (la1, lo1, la2, lo2))
    return R * np.arccos(np.clip(np.sin(la1) * np.sin(la2) + np.cos(la1) * np.cos(la2) * np.cos(lo1 - lo2), -1, 1))
out = []
for _, q in f.iterrows():
    ds = gc(q.lat, q.lon_E, v.Lon.values, v.Lat.values % 360)          # swapped: lat = Lon, lon = Lat
    da = gc(q.lat, q.lon_E, v.Lat.clip(-90, 90).values, v.Lon.values % 360)  # as labeled
    j = int(np.argmin(ds)); k = int(np.argmin(da))
    out.append(dict(name=q["name"], lat=q.lat, lon_E=q.lon_E, nearest_vent_swapped=int(v["Vent ID"].iloc[j]), dist_km_swapped=float(ds[j]),
                    nearest_vent_as_labeled=int(v["Vent ID"].iloc[k]), dist_km_as_labeled=float(da[k])))
o = pd.DataFrame(out); o.to_csv(f"{HERE}/faculae_crossmatch.csv", index=False)
S = {"n_faculae": int(len(o)),
     "swapped": dict(median_km=float(o.dist_km_swapped.median()), within_25km=int((o.dist_km_swapped <= 25).sum()), within_50km=int((o.dist_km_swapped <= 50).sum())),
     "as_labeled": dict(median_km=float(o.dist_km_as_labeled.median()), within_25km=int((o.dist_km_as_labeled <= 25).sum()), within_50km=int((o.dist_km_as_labeled <= 50).sum()))}
json.dump(S, open(f"{HERE}/faculae_crossmatch.json", "w"), indent=1); print(json.dumps(S, indent=1))
print(o[["name","lat","lon_E","dist_km_swapped","dist_km_as_labeled"]].round(1).to_string(index=False))
