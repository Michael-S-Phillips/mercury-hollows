"""Validation of automatic shadow depths against the manual shadow-length depths of Blewett et al. (2016).

Usage: python3 validate.py register | fetch | measure | stats      (each stage reads the previous stage's files)

register: Blewett points and Bickel detections share 416 images (2,230 Blewett depths). Blewett coordinates sit a median
  2.2 km from the camera-model detection positions (only 3% fall inside a detection box; 30% fall outside the frame), so a
  per-image translation is estimated: each frame gets an affine lon/lat -> pixel map fitted to the camera-model box corners
  (dep_common.FrameMap; median rms 0.45 px), Blewett points are projected into the frame, and the pixel shift (within +-8 km)
  that puts the most Blewett points inside detection boxes (60 m margin) is found exhaustively. A frame is accepted when it has
  >= 3 Blewett points, the best shift captures >= 3 and >= 50% of them, the best shift is unique (the best count at shifts
  > 1 km away is lower by >= 2) and well defined (all optimal shifts within 0.5 km). Depths are NOT used in the registration.
  A Blewett point is matched to the detection box that contains it after the shift (nearest box centre if several).
  Writes registration_frames.csv, registration_points.csv.
fetch: downloads the EDRs of accepted frames (dep_common.fetch_edr).
measure: runs shadow_measure.measure on every detection in accepted frames; writes validation_measurements.csv.
stats: builds validation.csv (one row per matched detection: auto depth vs the Blewett depth(s) inside its box) and
  validation_record.json; draws qc/validation_sheet_*.png.
Comparison: per detection box, the auto depth (longest shadow in the box) is compared with the deepest Blewett depth in the
  box (primary; a box can hold several hollows and the longest shadow belongs to the deepest) and with the median Blewett depth;
  the subset with exactly one Blewett point in the box is reported separately. Statistics: n, Pearson r of ln depths, median
  ratio auto/manual with a 95% bootstrap interval (2,000 resamples of images), MAD of ln ratio; restricted to emission <= 30 deg
  and to measurements without the flags no_shadow / edge_sun / frame_edge (the 'clean' set), with all rows also reported.
"""
import sys, os, json
import numpy as np, pandas as pd
import dep_common as C
import shadow_measure as M

HERE = C.HERE
SHIFT_KM, MARGIN_M = 8.0, 60.0


def register():
    k = C.load_bickel(); b = C.load_blewett()
    common = sorted(set(b.num) & set(k.num)); gk = dict(list(k.groupby("num"))); gb = dict(list(b.groupby("num")))
    fr, pts = [], []
    for n in common:
        q = gk[n]; fm = C.FrameMap(q); N = q.frame_samples.values[0]; L = q.frame_lines.values[0]; bb = gb[n]
        p = fm.to_px(bb.longitude.values, bb.latitude.values); mpp = fm.km_per_px * 1000
        S = int(SHIFT_KM * 1000 / mpp); m = MARGIN_M / mpp
        X0, X1, Y0, Y1 = q.x0.values * N, q.x1.values * N, q.y0.values * L, q.y1.values * L
        tot = np.zeros((2 * S + 1, 2 * S + 1), np.int32)
        for u, v in p:
            G = np.zeros((2 * S + 2, 2 * S + 2), np.int32)
            for a0, a1, c0, c1 in zip(X0, X1, Y0, Y1):
                lu = max(int(np.ceil(a0 - m - u)) + S, 0); hu = min(int(np.floor(a1 + m - u)) + S, 2 * S)
                lv = max(int(np.ceil(c0 - m - v)) + S, 0); hv = min(int(np.floor(c1 + m - v)) + S, 2 * S)
                if lu > hu or lv > hv: continue
                G[lv, lu] += 1; G[lv, hu + 1] -= 1; G[hv + 1, lu] -= 1; G[hv + 1, hu + 1] += 1
            tot += (G.cumsum(0).cumsum(1)[:-1, :-1] > 0)
        best = int(tot.max()); ys, xs = np.where(tot == best); xm, ym = xs.mean(), ys.mean()
        spread = float(np.hypot(xs - xm, ys - ym).max() * mpp / 1000)
        far = np.hypot(np.arange(2 * S + 1)[None, :] - xm, np.arange(2 * S + 1)[:, None] - ym) * mpp / 1000 > 1.0
        second = int(tot[far].max()) if far.any() else 0
        su, sv = xm - S, ym - S
        ok = (len(bb) >= 3) and best >= 3 and best >= 0.5 * len(bb) and best - second >= 2 and spread <= 0.5
        fr.append(dict(num=n, n_blewett=len(bb), n_det=len(q), best=best, second=second, shift_u_px=su, shift_v_px=sv,
                       shift_km=float(np.hypot(su, sv) * mpp / 1000), spread_km=spread, map_rms_px=fm.rms, m_per_px=mpp,
                       accepted=bool(ok), datetime=q.datetime.values[0], emission=float(q.emission_rel.values[0])))
        pu, pv = p[:, 0] + su, p[:, 1] + sv
        for i in range(len(bb)):
            inside = (pu[i] >= X0 - m) & (pu[i] <= X1 + m) & (pv[i] >= Y0 - m) & (pv[i] <= Y1 + m)
            fid = -1
            if inside.any():
                cand = np.where(inside)[0]; dc = np.hypot((X0[cand] + X1[cand]) / 2 - pu[i], (Y0[cand] + Y1[cand]) / 2 - pv[i])
                fid = int(q.fid.values[cand[np.argmin(dc)]])
            pts.append(dict(bid=int(bb.bid.values[i]), num=n, u=pu[i], v=pv[i], fid=fid, frame_accepted=bool(ok)))
    fr = pd.DataFrame(fr); pts = pd.DataFrame(pts)
    fr.to_csv(f"{HERE}/registration_frames.csv", index=False); pts.to_csv(f"{HERE}/registration_points.csv", index=False)
    a = fr[fr.accepted]
    print(f"frames {len(fr)}, accepted {len(a)}; Blewett points in accepted frames {a.n_blewett.sum()}, matched "
          f"{((pts.fid >= 0) & pts.frame_accepted).sum()} to {pts[(pts.fid >= 0) & pts.frame_accepted].fid.nunique()} boxes; "
          f"median shift {a.shift_km.median():.2f} km")


def fetch_many(pairs, workers=8):
    """pairs: list of (num, datetime). Parallel download; returns the number of frames present afterwards."""
    from concurrent.futures import ThreadPoolExecutor
    todo = [(n, d) for n, d in pairs if not os.path.exists(C.edr_path(n))]
    with ThreadPoolExecutor(workers) as ex:
        res = list(ex.map(lambda nd: C.fetch_edr(*nd), todo))
    miss = [n for (n, _), r in zip(todo, res) if r is None]
    if miss: print("missing", miss)
    return sum(os.path.exists(C.edr_path(n)) for n, _ in pairs)


def fetch():
    fr = pd.read_csv(f"{HERE}/registration_frames.csv", dtype={"num": str}); todo = fr[fr.accepted]
    print(f"have {fetch_many(list(zip(todo.num, todo.datetime)))}/{len(todo)}")


def measure_frames(nums, k, out_csv, tmax=165):
    import time
    t = time.time(); done = set()
    if os.path.exists(out_csv):
        done = set(pd.read_csv(out_csv, dtype={"num": str}).num)
    rows = []
    for n in nums:
        if n in done: continue
        if time.time() - t > tmax: break
        p = C.edr_path(n)
        if not os.path.exists(p): continue
        img, info = C.read_edr(p); q = k[k.num == n]; fs = M.foreshortening(info)
        for _, r in q.iterrows():
            N, L = info["samples"], info["lines"]
            box = (r.x0 * N, r.x1 * N, r.y0 * L, r.y1 * L)
            o = M.measure(img, info, box)
            inc = float(C.local_incidence(r.lat, r.lon180, info["subsol_lat"], info["subsol_lon"]))
            Lm = o["L_px"] * info["px"] * fs
            az = np.radians(info["subsol_az"]); D = np.radians(info["subsol_az"] - info["subsc_az"])
            fs_across = float(np.sqrt((np.sin(D) / np.cos(np.radians(info["emi"]))) ** 2 + np.cos(D) ** 2))
            rows.append(dict(fid=int(r.fid), num=n, image=f"EN{n}M", confidence=r.confidence, lat=r.lat, lon=r.lon180,
                             px_m=info["px"], inc_frame=info["inc"], inc_local=inc, emission=info["emi"], sun_az=info["subsol_az"],
                             sc_az=info["subsc_az"], foreshort=fs, B=o["B"], S=o["S"], T=o["T"], shadow_px=o["L_px"], shadow_m=Lm,
                             depth_m=Lm / np.tan(np.radians(inc)) if np.isfinite(Lm) else np.nan,
                             width_along_m=o["along_px"] * info["px"] * fs, width_across_m=o["across_px"] * info["px"] * fs_across,
                             poly_eqd_m=r.poly_eqd_km * 1000, n_shadow_runs=o["n_runs"], flags=o["flags"],
                             run_w=o["run_w"], run_u0=o["run_u0"], run_u1=o["run_u1"], x0=box[0], x1=box[1], y0=box[2], y1=box[3]))
        done.add(n)
    if rows:
        df = pd.DataFrame(rows)
        df.to_csv(out_csv, mode="a", header=not os.path.exists(out_csv), index=False)
    return len(done), len(rows)


def measure():
    fr = pd.read_csv(f"{HERE}/registration_frames.csv", dtype={"num": str}); nums = list(fr[fr.accepted].num)
    k = C.load_bickel(with_corners=False); k = k[k.num.isin(nums)]
    nd, nr = measure_frames(nums, k, f"{HERE}/validation_measurements.csv")
    print(f"frames done {nd}/{len(nums)}, new rows {nr}")


BAD = ("no_shadow", "edge_sun", "frame_edge")


def boot_median_ratio(df, col_a, col_m, nb=2000, seed=20261102):
    rng = np.random.default_rng(seed); lr = np.log(df[col_a] / df[col_m]).values; g = df.image.values
    imgs = np.unique(g); idx = {i: np.where(g == i)[0] for i in imgs}
    bt = [np.median(lr[np.concatenate([idx[i] for i in rng.choice(imgs, len(imgs))])]) for _ in range(nb)]
    return float(np.exp(np.median(lr))), [float(np.exp(np.percentile(bt, 2.5))), float(np.exp(np.percentile(bt, 97.5)))]


def summarize(df, col_a="depth_m", col_m="blewett_max_m"):
    d = df[(df[col_a] > 0) & (df[col_m] > 0)]
    if len(d) < 3: return dict(n=int(len(d)))
    la, lm = np.log(d[col_a]), np.log(d[col_m]); med, ci = boot_median_ratio(d, col_a, col_m)
    lr = la - lm
    return dict(n=int(len(d)), n_images=int(d.image.nunique()), pearson_r_log=float(np.corrcoef(la, lm)[0, 1]),
                spearman_rho=float(pd.Series(la.values).corr(pd.Series(lm.values), method="spearman")),
                median_ratio=med, median_ratio_ci95=ci, mad_log_ratio=float(np.median(np.abs(lr - np.median(lr)))),
                frac_within_factor_1p5=float((np.abs(lr) <= np.log(1.5)).mean()),
                median_auto_m=float(d[col_a].median()), median_manual_m=float(d[col_m].median()))


def stats():
    b = C.load_blewett().set_index("bid"); pts = pd.read_csv(f"{HERE}/registration_points.csv", dtype={"num": str})
    me = pd.read_csv(f"{HERE}/validation_measurements.csv", dtype={"num": str})
    mp = pts[(pts.fid >= 0) & pts.frame_accepted].copy(); mp["bdepth"] = b.loc[mp.bid, "depth"].values
    mp["bpx"] = b.loc[mp.bid, "px"].values; mp["binc"] = b.loc[mp.bid, "inc"].values
    g = mp.groupby("fid").agg(n_blewett=("bid", "size"), blewett_max_m=("bdepth", "max"), blewett_median_m=("bdepth", "median"),
                              blewett_min_m=("bdepth", "min"), blewett_inc=("binc", "first"), blewett_px=("bpx", "first"),
                              blewett_ids=("bid", lambda s: ";".join(map(str, s)))).reset_index()
    v = g.merge(me, on="fid", how="left")
    v["clean"] = ~v["flags"].fillna("").str.contains("|".join(BAD)) & (v.depth_m > 0)
    v["emi_le30"] = v.emission <= 30
    v["ratio_max"] = v.depth_m / v.blewett_max_m
    v.to_csv(f"{HERE}/validation.csv", index=False)
    R = dict(stage="validation", method_version=os.environ.get("METHOD_TAG", "v1"), frac=M.FRAC,
             registration=dict(frames_common=int(len(pd.read_csv(f"{HERE}/registration_frames.csv"))),
                               frames_accepted=int(pd.read_csv(f"{HERE}/registration_frames.csv").accepted.sum()),
                               points_matched=int(len(mp)), boxes_matched=int(len(g))))
    sets = {"emi_le30_clean": v[v.emi_le30 & v.clean], "emi_le30_all_measured": v[v.emi_le30 & (v.depth_m > 0)],
            "all_emission_clean": v[v.clean], "emi_le30_clean_single": v[v.emi_le30 & v.clean & (v.n_blewett == 1)],
            "emi_le30_clean_nomulti": v[v.emi_le30 & v.clean & ~v["flags"].fillna("").str.contains("multi")]}
    for nm, s in sets.items():
        R[nm] = dict(vs_max=summarize(s, "depth_m", "blewett_max_m"), vs_median=summarize(s, "depth_m", "blewett_median_m"))
    e = v[v.emi_le30 & v.clean]
    R["by_latitude_band_emi_le30_clean"] = {str(kk): summarize(gg) for kk, gg in e.groupby(pd.cut(e.lat.abs(), [0, 30, 50, 60, 90]), observed=True)}
    R["by_incidence_emi_le30_clean"] = {str(kk): summarize(gg) for kk, gg in e.groupby(pd.cut(e.inc_local, [0, 60, 70, 75, 90]), observed=True)}
    R["by_pixel_scale_emi_le30_clean"] = {str(kk): summarize(gg) for kk, gg in e.groupby(pd.cut(e.px_m, [0, 8, 14, 17, 25]), observed=True)}
    w = v[v.emi_le30].copy(); w["Lb_px"] = w.blewett_max_m * np.tan(np.radians(w.blewett_inc)) / w.blewett_px
    R["yield_emi_le30"] = dict(boxes=int(len(w)), clean_measured=int(w.clean.sum()), frac=float(w.clean.mean()),
                               frac_by_blewett_shadow_px={str(kk): [int(len(gg)), float(gg.clean.mean())] for kk, gg in
                                                          w.groupby(pd.cut(w.Lb_px, [0, 3, 5, 8, 12, 100]), observed=True)},
                               note="Failures are mostly no_shadow (contrast below MIN_CONTRAST or no run passing the rules): "
                                    "short (<= 3-5 px) shadows and lower incidence. Measured sets are therefore truncated "
                                    "against shallow hollows.")
    R["pass_criterion"] = dict(rule="median ratio within ~0.8-1.25, Pearson r(ln) >= 0.5, n >= 30 (emi_le30_clean, vs_max)")
    s0 = R["emi_le30_clean"]["vs_max"]
    R["pass_criterion"]["result"] = bool(0.8 <= s0["median_ratio"] <= 1.25 and s0["pearson_r_log"] >= 0.5 and s0["n"] >= 30)
    R["pass_criterion"]["note"] = ("Point estimate inside the window but the 95% interval reaches above 1.25: a marginal pass with a "
                                   "systematic overestimate of about +24%, larger at |lat| >= 60 deg (see by_latitude_band).")
    R["flag_counts_emi_le30"] = v[v.emi_le30]["flags"].fillna("").str.split(";").explode().value_counts().to_dict()
    R["no_shadow_frac_emi_le30"] = float(v[v.emi_le30]["flags"].fillna("").str.contains("no_shadow").mean())
    json.dump(R, open(f"{HERE}/validation_record.json", "w"), indent=1)
    print(json.dumps({k_: R[k_] for k_ in ["registration", "emi_le30_clean", "emi_le30_clean_single", "emi_le30_all_measured"]}, indent=1))


if __name__ == "__main__":
    dict(register=register, fetch=fetch, measure=measure, stats=stats)[sys.argv[1]]()
