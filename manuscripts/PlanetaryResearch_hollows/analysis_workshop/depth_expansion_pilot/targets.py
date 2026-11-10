"""New shadow depths for under-sampled hollows (run only because validate.py passed, marginally; see validation_record.json).

Usage: python3 targets.py select | run | compare | qc     (run is resumable; repeat until it reports all frames done)

select: confident Bickel detections (release confidence >= 0.75) in NAC frames with pixel scale <= 20 m, incidence 60-80 deg
  and emission <= 30 deg (release per-image geometry) = the 'base' set (1,866 detections, 583 frames). Hosts are assigned with
  the camera-model positions (dep_common.assign_hosts: smallest Herrick et al. 2018 crater >= 5 km whose 1.1 radii contain the
  point, as in depth_nested.py). Targets: (a) host D < 20 km (75 detections, 32 hosts, 42 frames); (b) |lat| >= 60 deg (203
  detections, 134 frames). Writes targets_selection.csv.
run: measures every base detection with shadow_measure.measure (the validated v2 method), streaming EDRs from the PDS Imaging
  Node (frames outside the target and validation sets are deleted after use to save VM disk). Writes base_measurements.csv.
compare: writes new_depths.csv (targets) and targets_record.json. Depth summaries are medians with 95% intervals from 2,000
  bootstrap resamples of images, of depths corrected to the median pixel scale and incidence of the Blewett set with the
  paper's OLS ln h ~ ln px + ln tan(inc) fitted to the 2,518 Blewett depths (as depth_nested.py). Three comparisons:
  (1) method-internal: target groups vs the rest of the auto-measured base set (same method, same image selection, so the
      method's overall bias cancels);
  (2) against Blewett depths of comparable host size / latitude, with the auto depths divided by the validation median
      ratio (auto/manual = 1.236 overall; 1.62 at |lat| >= 60, n = 23, shown as a sensitivity);
  (3) raw auto depths against Blewett depths (uncalibrated, for reference).
qc: contact sheet of 30 random target measurements (qc/targets_sheet_random.png).
"""
import sys, os, json, time
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
import dep_common as C
import shadow_measure as M
import validate as V

HERE = C.HERE
F_SEL, F_BASE = f"{HERE}/targets_selection.csv", f"{HERE}/base_measurements.csv"


def select():
    k = C.load_bickel(with_corners=False)
    t = k[(k.confidence >= 0.75) & (k.px_rel <= 20) & (k.incidence_rel >= 60) & (k.incidence_rel <= 80) & (k.emission_rel <= 30)].copy()
    hid, hD, rho = C.assign_hosts(t.lat.values, t.lon180.values % 360)
    t["host_id"], t["host_D_km"], t["host_rho"] = hid, hD, rho
    t["small_host"] = t.host_D_km < 20; t["poleward"] = t.lat.abs() >= 60
    cols = ["fid", "mdis_id", "num", "datetime", "confidence", "lat", "lon180", "lat_rel", "lon_rel", "px_rel", "incidence_rel",
            "emission_rel", "subsol_az_rel", "x0", "x1", "y0", "y1", "poly_eqd_km", "host_id", "host_D_km", "host_rho", "small_host", "poleward"]
    t[cols].to_csv(F_SEL, index=False)
    a, b = t[t.small_host], t[t.poleward]
    print(f"base {len(t)} dets / {t.num.nunique()} frames; small hosts {len(a)} / {a.host_id.nunique()} hosts / {a.num.nunique()} frames; "
          f"poleward {len(b)} / {b.num.nunique()} frames")


def run(tmax=100, batch=16):
    t0 = time.time(); s = pd.read_csv(F_SEL, dtype={"num": str, "host_id": str})
    keep = set(s[s.small_host | s.poleward].num) | set(pd.read_csv(f"{HERE}/registration_frames.csv", dtype={"num": str}).query("accepted").num)
    done = set(pd.read_csv(F_BASE, dtype={"num": str}).num) if os.path.exists(F_BASE) else set()
    frames = [n for n in s.num.unique() if n not in done]
    k = C.load_bickel(with_corners=False); k = k[k.fid.isin(s.fid)]
    while frames and time.time() - t0 < tmax:
        chunk, frames = frames[:batch], frames[batch:]
        dts = s.drop_duplicates("num").set_index("num").datetime
        V.fetch_many([(n, dts[n]) for n in chunk])
        V.measure_frames(chunk, k, F_BASE, tmax=1e9)
        for n in chunk:
            if n not in keep and os.path.exists(C.edr_path(n)): os.remove(C.edr_path(n))
    done = set(pd.read_csv(F_BASE, dtype={"num": str}).num)
    print(f"frames measured {len(done)}/{s.num.nunique()}")


def correction():
    d = pd.read_csv(C.F_BLEWETT); d["ltan"] = np.log(np.tan(np.radians(d.inc)))
    c = smf.ols("lnh ~ lnpx + ltan", data=d).fit()
    return c.params, float(d.lnpx.median()), float(d.ltan.median())


def corrected(h, px, inc, P, mp, mt):
    return np.exp(np.log(h) - P["lnpx"] * (np.log(px) - mp) - P["ltan"] * (np.log(np.tan(np.radians(inc))) - mt))


def med_ci(vals, imgs, nb=2000, seed=20261102):
    vals, imgs = np.asarray(vals), np.asarray(imgs); ok = np.isfinite(vals) & (vals > 0); vals, imgs = vals[ok], imgs[ok]
    if len(vals) < 3: return dict(n=int(len(vals)))
    rng = np.random.default_rng(seed); u = np.unique(imgs); g = {i: vals[imgs == i] for i in u}
    bt = [np.median(np.concatenate([g[i] for i in rng.choice(u, len(u))])) for _ in range(nb)]
    return dict(n=int(len(vals)), n_images=int(len(u)), median=float(np.median(vals)),
                ci95=[float(np.percentile(bt, 2.5)), float(np.percentile(bt, 97.5))])


def ratio_ci(a_vals, a_imgs, b_vals, b_imgs, nb=2000, seed=7):
    """Ratio of medians (a/b) with an image-bootstrap interval (independent resampling of both groups)."""
    rng = np.random.default_rng(seed)
    def prep(v, i):
        v, i = np.asarray(v), np.asarray(i); ok = np.isfinite(v) & (v > 0); v, i = v[ok], i[ok]; u = np.unique(i)
        return v, {x: v[i == x] for x in u}, u
    va, ga, ua = prep(a_vals, a_imgs); vb, gb, ub = prep(b_vals, b_imgs)
    if len(va) < 3 or len(vb) < 3: return None
    bt = [np.median(np.concatenate([ga[x] for x in rng.choice(ua, len(ua))])) / np.median(np.concatenate([gb[x] for x in rng.choice(ub, len(ub))]))
          for _ in range(nb)]
    return dict(ratio=float(np.median(va) / np.median(vb)), ci95=[float(np.percentile(bt, 2.5)), float(np.percentile(bt, 97.5))])


def compare():
    s = pd.read_csv(F_SEL, dtype={"num": str, "host_id": str}); m = pd.read_csv(F_BASE, dtype={"num": str})
    vr = json.load(open(f"{HERE}/validation_record.json"))
    cal = vr["emi_le30_clean"]["vs_max"]["median_ratio"]; cal_pole = vr["by_latitude_band_emi_le30_clean"]["(60, 90]"]["median_ratio"]
    x = s.drop(columns=["x0", "x1", "y0", "y1"]).merge(m.drop(columns=["num", "confidence", "lat", "lon"]), on="fid", how="left")  # boxes in pixels from m
    x["qc_clean"] = ~x["flags"].fillna("").str.contains("|".join(V.BAD)) & (x.depth_m > 0)
    P, mp, mt = correction()
    x["depth_corr_m"] = np.where(x.qc_clean, corrected(x.depth_m.where(x.depth_m > 0), x.px_m, x.inc_local, P, mp, mt), np.nan)
    x["depth_cal_m"] = x.depth_m / cal; x["depth_corr_cal_m"] = x.depth_corr_m / cal
    x["group"] = np.where(x.small_host & x.poleward, "small_host;poleward", np.where(x.small_host, "small_host", np.where(x.poleward, "poleward", "")))
    tg = x[x.small_host | x.poleward].copy()
    out = tg.rename(columns={"lon180": "lon", "mdis_id": "bickel_cdr_id", "sun_az": "sun_az_deg", "shadow_px": "shadow_len_px",
                             "shadow_m": "shadow_len_m", "inc_local": "inc_deg", "inc_frame": "inc_frame_deg", "emission": "emission_deg",
                             "width_along_m": "width_along_sun_m", "width_across_m": "width_across_sun_m", "poly_eqd_m": "release_poly_eqdiam_m"})
    out.insert(1, "release_row", out.fid - 1)
    cols = ["fid", "release_row", "image", "bickel_cdr_id", "group", "confidence", "lat", "lon", "px_m", "inc_deg", "inc_frame_deg", "emission_deg",
            "sun_az_deg", "foreshort", "shadow_len_px", "shadow_len_m", "depth_m", "depth_cal_m", "depth_corr_m", "depth_corr_cal_m",
            "width_along_sun_m", "width_across_sun_m", "release_poly_eqdiam_m", "n_shadow_runs", "host_id", "host_D_km", "host_rho",
            "flags", "qc_clean", "B", "S", "T", "x0", "x1", "y0", "y1", "run_w", "run_u0", "run_u1"]
    out[cols].to_csv(f"{HERE}/new_depths.csv", index=False)
    # ---------------- Blewett reference groups (paper's host assignment on Blewett positions)
    b = C.load_blewett(); hid, hD, _ = C.assign_hosts(b.latitude.values, b.longitude.values)
    b["host_D_km"] = hD; b["hc"] = corrected(b.depth, b.px, b.inc, P, mp, mt); b["img"] = b.filename
    c = x[x.qc_clean]; R = dict(correction=dict(b_lnpx=float(P["lnpx"]), b_ltan=float(P["ltan"]), px_median_m=float(np.exp(mp)),
                                               inc_median_deg=float(np.degrees(np.arctan(np.exp(mt))))),
                               calibration=dict(overall=cal, poleward_validation=cal_pole))
    R["counts"] = dict(base_dets=int(len(x)), base_frames=int(x.num.nunique()), base_measured_clean=int(len(c)),
                       base_clean_frac=float(x.qc_clean.mean()),
                       small_host=dict(dets=int(x.small_host.sum()), hosts=int(x[x.small_host].host_id.nunique()), frames=int(x[x.small_host].num.nunique()),
                                       clean=int(c.small_host.sum()), clean_hosts=int(c[c.small_host].host_id.nunique())),
                       poleward=dict(dets=int(x.poleward.sum()), frames=int(x[x.poleward].num.nunique()), clean=int(c.poleward.sum()),
                                     clean_frames=int(c[c.poleward].num.nunique())))
    R["flag_counts_targets"] = tg["flags"].fillna("").str.split(";").explode().value_counts().to_dict()
    def grp(df, col, imgcol):
        return med_ci(df[col].values, df[imgcol].values)
    ref = c[~c.small_host & (c.host_D_km >= 20)]; lowlat = c[~c.poleward]
    R["auto_internal"] = dict(
        small_host=grp(c[c.small_host], "depth_corr_m", "num"), host_ge20=grp(ref, "depth_corr_m", "num"),
        no_host=grp(c[c.host_D_km.isna()], "depth_corr_m", "num"),
        ratio_small_vs_ge20=ratio_ci(c[c.small_host].depth_corr_m, c[c.small_host].num, ref.depth_corr_m, ref.num),
        poleward=grp(c[c.poleward], "depth_corr_m", "num"), lat_lt60=grp(lowlat, "depth_corr_m", "num"),
        lat_50_60=grp(c[(c.lat.abs() >= 50) & (c.lat.abs() < 60)], "depth_corr_m", "num"),
        ratio_poleward_vs_lt60=ratio_ci(c[c.poleward].depth_corr_m, c[c.poleward].num, lowlat.depth_corr_m, lowlat.num),
        ratio_poleward_vs_50_60=ratio_ci(c[c.poleward].depth_corr_m, c[c.poleward].num,
                                         c[(c.lat.abs() >= 50) & (c.lat.abs() < 60)].depth_corr_m, c[(c.lat.abs() >= 50) & (c.lat.abs() < 60)].num),
        small_host_raw=grp(c[c.small_host], "depth_m", "num"), poleward_raw=grp(c[c.poleward], "depth_m", "num"))
    bs, bl = b[b.host_D_km < 20], b[b.host_D_km >= 20]; bp = b[b.latitude.abs() >= 60]; b56 = b[(b.latitude.abs() >= 50) & (b.latitude.abs() < 60)]
    R["blewett"] = dict(small_host=grp(bs, "hc", "img"), host_ge20=grp(bl, "hc", "img"), poleward=grp(bp, "hc", "img"),
                        lat_50_60=grp(b56, "hc", "img"), lat_lt60=grp(b[b.latitude.abs() < 60], "hc", "img"), all=grp(b, "hc", "img"),
                        small_host_hosts=int(bs.drop_duplicates("host_D_km").shape[0]))
    sh, po = c[c.small_host], c[c.poleward]
    R["auto_calibrated_vs_blewett"] = dict(
        small_host=dict(auto_cal=grp(sh, "depth_corr_cal_m", "num"),
                        ratio_to_blewett_small=ratio_ci(sh.depth_corr_cal_m, sh.num, bs.hc, bs.img),
                        ratio_to_blewett_ge20=ratio_ci(sh.depth_corr_cal_m, sh.num, bl.hc, bl.img)),
        poleward=dict(auto_cal=grp(po, "depth_corr_cal_m", "num"),
                      auto_cal_poleward_factor=med_ci((po.depth_corr_m / cal_pole).values, po.num.values),
                      ratio_to_blewett_poleward=ratio_ci(po.depth_corr_cal_m, po.num, bp.hc, bp.img),
                      ratio_to_blewett_poleward_polecal=ratio_ci(po.depth_corr_m / cal_pole, po.num, bp.hc, bp.img),
                      ratio_to_blewett_50_60=ratio_ci(po.depth_corr_cal_m, po.num, b56.hc, b56.img)))
    R["auto_raw_vs_blewett"] = dict(small_host=ratio_ci(sh.depth_corr_m, sh.num, bs.hc, bs.img), poleward=ratio_ci(po.depth_corr_m, po.num, bp.hc, bp.img))
    json.dump(R, open(f"{HERE}/targets_record.json", "w"), indent=1)
    print(json.dumps(R, indent=1)[:6000])


def qc(n=30, seed=11):
    import qc_render as Q
    d = pd.read_csv(f"{HERE}/new_depths.csv", dtype={"image": str}); d = d[d.qc_clean].sample(n, random_state=seed)
    tiles = []
    for _, r in d.iterrows():
        num = r.image[2:12]; p = C.edr_path(num)
        if not os.path.exists(p): V.fetch_many([(num, s_dt(num))])
        img, info = C.read_edr(p)
        tiles.append(Q.tile(img, info, (r.x0, r.x1, r.y0, r.y1), dict(run_u0=r.run_u0, run_u1=r.run_u1, run_w=r.run_w), f"{num} f{r.fid} {r.group}\nh={r.depth_m:.0f}m ({r.shadow_len_px:.1f}px) px{r.px_m:.0f} i{r.inc_deg:.0f}\nlat{r.lat:.0f} hostD {r.host_D_km:.0f} {r['flags'] if isinstance(r['flags'], str) else ''}"))
    print(Q.sheet(tiles, ncol=6, path=f"{HERE}/qc/targets_sheet_random.png"))


def s_dt(num):
    s = pd.read_csv(F_SEL, dtype={"num": str}); return s[s.num == num].datetime.iloc[0]


if __name__ == "__main__":
    dict(select=select, run=run, compare=compare, qc=qc)[sys.argv[1]]()
