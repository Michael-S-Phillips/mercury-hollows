import os as _os_repo  # added for the archive: repository root, or $MERCURY_HOLLOWS_ROOT
_REPO = _os_repo.environ.get("MERCURY_HOLLOWS_ROOT") or _os_repo.path.abspath(_os_repo.path.join(_os_repo.path.dirname(__file__), ".."))
"""Within-crater observed/expected hollow counts by regional slope bin, overall, by aspect, and by terrain.
Expected counts spread each crater's hollows over its completeness-weighted area in proportion to area, so differences
between craters cancel. 95% intervals from 2000 bootstrap resamples of craters. Reads occ_by_crater_bin.csv."""
import json, os, numpy as np, pandas as pd
HE = os.path.expanduser(_REPO + "/slope_test")
d = pd.read_csv(f"{HE}/occ_by_crater_bin.csv")
rng = np.random.default_rng(20261026)
def oe(sub, ref, col="n", groups=("sbin",)):
    """O/E per group; expected from each crater's own count and area over `ref` (rows defining the reference area)."""
    tot = ref.groupby("uid").agg(N=(col, "sum"), A=("A_km2", "sum"))
    s = sub.merge(tot, left_on="uid", right_index=True)
    s["E"] = s.N * s.A_km2 / s.A
    g = s.groupby(["uid", *groups]).agg(O=(col, "sum"), E=("E", "sum"), A=("A_km2", "sum")).reset_index()
    uids = g.uid.unique(); res = {}
    piv_O = g.pivot_table(index="uid", columns=list(groups), values="O", aggfunc="sum", fill_value=0)
    piv_E = g.pivot_table(index="uid", columns=list(groups), values="E", aggfunc="sum", fill_value=0)
    piv_A = g.pivot_table(index="uid", columns=list(groups), values="A", aggfunc="sum", fill_value=0)
    O, E = piv_O.values, piv_E.values
    bt = []
    for _ in range(2000):
        ix = rng.integers(0, len(uids), len(uids)); bt.append(O[ix].sum(0) / np.maximum(E[ix].sum(0), 1e-12))
    bt = np.array(bt)
    for k, col_ in enumerate(piv_O.columns):
        res[str(col_)] = dict(O=int(O[:, k].sum()), E=round(float(E[:, k].sum()), 1), ratio=round(float(O[:, k].sum() / max(E[:, k].sum(), 1e-12)), 3),
                             lo=round(float(np.percentile(bt[:, k], 2.5)), 3), hi=round(float(np.percentile(bt[:, k], 97.5)), 3),
                             A_km2=round(float(piv_A.values[:, k].sum())), n_craters_with_area=int((piv_A.values[:, k] > 0).sum()))
    return res
R = {}
for bl in ("s1", "s3"):
    x = d[d.baseline == bl]
    R[bl] = {"all_hollows": oe(x, x), "all_fields": oe(x, x, col="nf"),
             "by_aspect": oe(x[x.sbin >= 5], x, groups=("aspect", "sbin")),
             "wall_only": oe(x[x.terrain == "wall"], x[x.terrain == "wall"]),
             "peak_only": oe(x[x.terrain == "peak"], x[x.terrain == "peak"]),
             "floor_only": oe(x[x.terrain == "floor"], x[x.terrain == "floor"]),
             "by_terrain": oe(x, x, groups=("terrain",))}
    sl = x.groupby("sbin").A_km2.sum(); R[bl]["area_fraction"] = {str(k): round(float(v / sl.sum()), 4) for k, v in sl.items()}
json.dump(R, open(f"{HE}/occ_analysis_record.json", "w"), indent=1)
for bl in ("s1",):
    for k in ("all_hollows", "all_fields", "wall_only", "peak_only", "floor_only", "by_terrain"):
        print(bl, k); [print("   ", kk, v["O"], v["E"], v["ratio"], (v["lo"], v["hi"]), v["A_km2"]) for kk, v in R[bl][k].items()]
    print(bl, "by_aspect"); [print("   ", kk, v["O"], v["E"], v["ratio"], (v["lo"], v["hi"])) for kk, v in R[bl]["by_aspect"].items()]
print("s3 all_hollows"); [print("   ", kk, v["O"], v["E"], v["ratio"], (v["lo"], v["hi"])) for kk, v in R["s3"]["all_hollows"].items()]
