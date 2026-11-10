# g19 (8 Nov 2026), new script: checks on the traced knob comparison that use no curvature classes. (1) the can't-tell craters: their hollows
# and the hollows and area within 0.15 R; a bound counting all of that as peak; (2) which bootstrap component dominates the interval
# (knobs only, hosts only, sample only). Reads the same inputs as knob_ratio_design.py. Writes knob_checks_nocurv.json.
import os, json, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); RR = os.path.dirname(os.path.dirname(HERE)); B = 4000
c = pd.read_csv(f"{HERE}/class12_peaks_traced_r3k.csv")
TRc = json.load(open(f"{os.path.dirname(HERE)}/traces/traces_clean_r3k.json"))["craters"]; c["status"] = c.uid.map({int(k): v["status"] for k, v in TRc.items()}).fillna("untraced")
c["note"] = c.uid.map({int(k): v.get("note", "") for k, v in TRc.items()}).fillna("")
Dl = pd.DataFrame(json.load(open(f"{RR}/g17_tracer_round2/selection.json"))["group_D"])
kh = pd.read_csv(f"{RR}/g9_draft_reruns/knobs/out_v3_gridv3/knob_hollows.csv"); KN, KA = kh.n_knob.values.astype(float), kh.Aw_knob.values.astype(float)
H = c[c.n_assigned > 0].copy(); S = c[c.n_assigned == 0].merge(Dl[["uid", "pi"]], on="uid")
sk_h = H[H.status == "skip"]; sk_s = S[S.status == "skip"]
def cols(df, bound):
    n = np.where(df.status == "done", df.n_peak_tr, 0.0); a = np.where(df.status == "done", df.A_peak_tr, 0.0)
    if bound: s = (df.status == "skip").values; n = np.where(s, df.n_p015, n); a = np.where(s, df.A_p015, a)
    return n.astype(float), a.astype(float)
def est(bound=False, comp=("k", "h", "s")):
    hn, ha = cols(H, bound); sn, sa = cols(S, bound); pi = S.pi.values; cert = pi >= 1
    keep = np.ones(len(S), bool) if bound else (S.status == "done").values | cert
    m_all = (~cert).sum(); nc0 = np.where(~cert & keep)[0]
    def ratio(ik, ih, isd):
        A = ha[ih].sum() + (sa[cert & keep] / pi[cert & keep]).sum() + (sa[isd] / pi[isd]).sum() * (m_all / len(nc0))
        N = hn[ih].sum(); return (KN[ik].sum() / KA[ik].sum()) / (N / A) if N > 0 else np.nan, N, A
    r0, N0, A0 = ratio(np.arange(len(KN)), np.arange(len(hn)), nc0)
    rng = np.random.default_rng(1); q = []
    for b in range(B):
        ik = rng.integers(0, len(KN), len(KN)) if "k" in comp else np.arange(len(KN))
        ih = rng.integers(0, len(hn), len(hn)) if "h" in comp else np.arange(len(hn))
        isd = nc0[rng.integers(0, len(nc0), len(nc0))] if "s" in comp else nc0
        q.append(ratio(ik, ih, isd)[0])
    q = np.array(q); q = q[np.isfinite(q)]
    return dict(ratio=round(float(r0), 2), ci=[round(float(np.percentile(q, 2.5)), 2), round(float(np.percentile(q, 97.5)), 2)], peak_hollows=int(N0), peak_area_km2=round(float(A0)),
                sd_log=round(float(np.std(np.log(q))), 3))
out = dict(created="2026-11-08",
           cant_tell_hosts=sk_h[["uid", "D", "kclass", "n_assigned", "n_p015", "A_p015", "note"]].to_dict("records"),
           cant_tell_sample=sk_s[["uid", "D", "kclass", "pi", "n_p015", "A_p015", "note"]].to_dict("records"),
           n_knob_hollows=int(KN.sum()), n_knobs=int(len(KN)), n_knobs_with_hollows=int((KN > 0).sum()),
           main=est(), bound_cant_tell_within_0p15R_as_peak=est(bound=True),
           components=dict(knobs_only=est(comp=("k",)), hosts_only=est(comp=("h",)), sample_only=est(comp=("s",))))
json.dump(out, open(f"{HERE}/knob_checks_nocurv.json", "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
