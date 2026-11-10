#!/usr/bin/env python3
"""g19 (7 Nov 2026), new script: knob hollow density / hollow density on the peaks of Kinczyk et al. (2020) Class 1-2 craters, with peak
outlines traced by hand (tracer groups A-D, traces_clean_r3k.json). Class 1-2 craters 30-250 km (2,160, as g13). Hosts (assigned hollows > 0)
are all traced: their hollows on traced peaks and their traced peak area enter exactly. Craters without hollows enter through the PPS sample
of tracer group D (g17 selection.json, inclusion probabilities pi): Horvitz-Thompson total of traced peak area, sum(A_i / pi_i). Craters
marked can't tell: hosts contribute no hollows and no area (as a crater with no peak); a sampled crater without hollows is dropped and the
non-certainty part of the sample is scaled by m / (m - dropped). Intervals: 4,000 bootstrap replicates resampling knobs, hosts, and the
non-certainty sample units independently (certainty units fixed), default_rng(1). Sensitivity: can't-tell craters with curvature v3 peaks;
curvature v3 for every crater (the g13 definition). Writes class12_peaks_traced_r3k.csv and knob_ratio_design.json."""
import os, glob, json, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); RR = os.path.dirname(os.path.dirname(HERE)); B = 4000
c = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{HERE}/_chunks_r3k/c12_*.csv"))]).drop_duplicates("uid").sort_values("uid").reset_index(drop=True)
c.to_csv(f"{HERE}/class12_peaks_traced_r3k.csv", index=False)
TRc = json.load(open(f"{os.path.dirname(HERE)}/traces/traces_clean_r3k.json"))["craters"]; status = {int(k): v["status"] for k, v in TRc.items()}
sel = json.load(open(f"{RR}/g17_tracer_round2/selection.json")); Dl = pd.DataFrame(sel["group_D"])
kh = pd.read_csv(f"{RR}/g9_draft_reruns/knobs/out_v3_gridv3/knob_hollows.csv"); KN, KA = kh.n_knob.values.astype(float), kh.Aw_knob.values.astype(float)
c["status"] = c.uid.map(status).fillna("untraced")
H = c[c.n_assigned > 0].copy(); NH = c[c.n_assigned == 0]
S = NH.merge(Dl[["uid", "pi"]], on="uid", how="inner"); assert len(S) == len(Dl), (len(S), len(Dl))
assert (H.status != "untraced").all(), H[H.status == "untraced"].uid.tolist()
def cols(df, mode):
    if mode == "curv": return df.n_peak_v3.values.astype(float), df.A_peak_v3.values.astype(float)
    n = np.where(df.status == "done", df.n_peak_tr, 0.0); a = np.where(df.status == "done", df.A_peak_tr, 0.0)
    if mode == "traced_skip_curv": sk = (df.status == "skip").values; n = np.where(sk, df.n_peak_v3, n); a = np.where(sk, df.A_peak_v3, a)
    return n.astype(float), a.astype(float)
def estimate(mode, rng=None):
    hn, ha = cols(H, mode); sn, sa = cols(S, mode); pi = S.pi.values; cert = pi >= 1
    keep = np.ones(len(S), bool) if mode != "traced" else (S.status == "done").values | cert   # drop sampled can't-tell units (none is a certainty unit)
    def total(idx_h, idx_s):
        ht_c = (sa[cert & keep] / pi[cert & keep]).sum()
        nc = np.where(~cert & keep)[0] if idx_s is None else idx_s
        m_all = (~cert).sum(); m = len(np.where(~cert & keep)[0])
        ht_n = (sa[nc] / pi[nc]).sum() * (m_all / m)
        ih = np.arange(len(hn)) if idx_h is None else idx_h
        return hn[ih].sum(), ha[ih].sum() + ht_c + ht_n, ha[ih].sum(), ht_c + ht_n
    N, A, Ah, Aht = total(None, None); kd = KN.sum() / KA.sum(); r = kd / (N / A)
    rng = np.random.default_rng(1); nc0 = np.where(~cert & keep)[0]; q = []
    for b in range(B):
        ik = rng.integers(0, len(KN), len(KN)); ih = rng.integers(0, len(hn), len(hn)); isd = nc0[rng.integers(0, len(nc0), len(nc0))]
        Nb, Ab, _, _ = total(ih, isd); q.append((KN[ik].sum() / KA[ik].sum()) / (Nb / Ab) if Nb > 0 else np.nan)
    q = np.array(q); q = q[np.isfinite(q)]
    return dict(peak_hollows=int(N), hosts=int(len(H)), hosts_with_peak_hollows=int((hn > 0).sum()), host_peak_area_km2=round(float(Ah)), nonhost_peak_area_HT_km2=round(float(Aht)),
                peak_density_per_1000km2=round(1000 * N / A, 3), knob_density_per_1000km2=round(1000 * kd, 3), knob_over_peak=round(float(r), 2),
                ci=[round(float(np.percentile(q, 2.5)), 2), round(float(np.percentile(q, 97.5)), 2)], n_boot=int(len(q)))
out = dict(created="2026-11-07", note=__doc__, n_class12=int(len(c)), n_hosts=int(len(H)), n_nonhosts=int(len(NH)), sample=int(len(S)), sample_certainty=int((S.pi >= 1).sum()),
           host_status={k: int(v) for k, v in H.status.value_counts().items()}, sample_status={k: int(v) for k, v in S.status.value_counts().items()},
           sample_with_traced_peak=int(((S.status == "done") & (S.A_peak_tr > 0)).sum()),
           design_check=dict(curv_nonhost_total=round(float(NH.A_peak_v3.sum())), curv_nonhost_HT_from_sample=round(float((S.A_peak_v3 / S.pi).sum()))),
           traced=estimate("traced"), traced_skip_curv=estimate("traced_skip_curv"), curvature_all=estimate("curv"))
json.dump(out, open(f"{HERE}/knob_ratio_design.json", "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "note"}, indent=1))
