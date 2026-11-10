#!/usr/bin/env python3
"""Decision sheets for v3 vs v3b (6 Nov 2026): what the uncertain (U) detections look like.

U1: U detections that are the only hollows of a v3 host crater (the hosts that v3b loses), random 24, spread over size.
U2: other U detections, random 24, spread over size.
REF: 12 kept (H) and 12 removed (B/W/O) detections from screen set B, at the same sizes, for calibration.
Tiles are the screening tiles (qc/tiles/<hollow_id>.png: best EDR, yellow = detection box, cyan = toward the Sun).
Each tile gets the screener's note. Writes qc/decide_U1.png, decide_U2.png, decide_REF.png and decide_ids.csv."""
import os, textwrap, numpy as np, pandas as pd
from PIL import Image, ImageDraw, ImageFont
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.dirname(HERE); TD = f"{HERE}/tiles"
rng = np.random.default_rng(20261106)
s = pd.read_csv(f"{OUT}/screen_classes.csv"); ss = pd.read_csv(f"{OUT}/screen_set.csv")[["hollow_id", "lat", "lon", "host_uid", "host_n09"]]
s = s.merge(ss, on="hollow_id", how="left")
h = pd.read_csv(f"{OUT}/host_assignments_v3.csv", low_memory=False)
h = h[(h.rho <= 1.1)].merge(s[["hollow_id", "class"]], on="hollow_id", how="left"); h["class"] = h["class"].fillna("unscreened")
allU = h.groupby("host_uid")["class"].agg(lambda c: (c == "U").all())
lost_hosts = set(allU[allU].index); h_lost = h[h.host_uid.isin(lost_hosts)]
u = s[s["class"] == "U"].copy(); u["sole_host"] = u.hollow_id.isin(h_lost.hollow_id)
print("v3 hosts (rho<=1.1):", h.host_uid.nunique(), " hosts whose hollows are all U:", len(lost_hosts), " U detections in them:", int(u.sole_host.sum()), "of", len(u))

def spread(df, n):
    df = df[[os.path.exists(f"{TD}/{i}.png") for i in df.hollow_id]].sort_values("D_eq_km").reset_index(drop=True)
    if len(df) <= n: return df
    q = np.array_split(np.arange(len(df)), n); return df.iloc[[rng.choice(k) for k in q]].reset_index(drop=True)

try: F = ImageFont.truetype("DejaVuSans.ttf", 13); FB = ImageFont.truetype("DejaVuSans-Bold.ttf", 20)
except Exception: F = FB = ImageFont.load_default()
T, CAP, NC = 266, 46, 6
def sheet(df, title, path, label):
    nr = int(np.ceil(len(df) / NC)); W = NC * 244; S = Image.new("RGB", (W, 44 + nr * (T + CAP + 6)), (250, 250, 248)); d = ImageDraw.Draw(S)
    d.text((10, 10), title, fill=(10, 10, 10), font=FB)
    for i, r in df.iterrows():
        x, y = (i % NC) * 244 + 2, 44 + (i // NC) * (T + CAP + 6)
        S.paste(Image.open(f"{TD}/{int(r.hollow_id)}.png"), (x, y))
        d.rectangle([x + 2, y + 2, x + 30, y + 24], fill=(10, 10, 10)); d.text((x + 5, y + 3), f"{i + 1}", fill=(255, 255, 255), font=FB)
        note = str(r.note).replace(" [subagent]", "").replace("classified by Claude (main session)", "(no note)")
        lines = textwrap.wrap(f"{label(r)}: {note}", 34)[:3]
        d.multiline_text((x + 2, y + T + 2), "\n".join(lines), fill=(30, 30, 30), font=F, spacing=1)
    S.save(path)

u1 = spread(u[u.sole_host], 24); u2 = spread(u[~u.sole_host], 24)
sz = (u.D_eq_km.min(), u.D_eq_km.max())
b = s[(s.set == "B") & s.D_eq_km.between(*sz)]
ref = pd.concat([spread(b[b["class"] == "H"], 12).assign(kind="kept"), spread(b[b["class"].isin(["B", "W", "O"])], 12).assign(kind="removed")], ignore_index=True)
sheet(u1, f"U1  uncertain detections that are a host's only hollows (v3 keeps, v3b drops the host)  [{len(u1)} of {int(u.sole_host.sum())}]", f"{HERE}/decide_U1.png", lambda r: "U")
sheet(u2, f"U2  other uncertain detections (host keeps other hollows either way)  [{len(u2)} of {int((~u.sole_host).sum())}]", f"{HERE}/decide_U2.png", lambda r: "U")
sheet(ref, "REF  tiles 1-12 kept as hollows (H); tiles 13-24 removed as crater bowl, wall or other (B/W/O); same size range", f"{HERE}/decide_REF.png",
      lambda r: f"{'kept' if r['class'] == 'H' else 'removed'} {r['class']}")
pd.concat([u1.assign(sheet="U1"), u2.assign(sheet="U2"), ref.assign(sheet="REF")]).assign(tile=lambda x: x.groupby("sheet").cumcount() + 1)[
    ["sheet", "tile", "hollow_id", "class", "D_eq_km", "px_best", "lat", "lon", "host_uid", "note"]].to_csv(f"{HERE}/decide_ids.csv", index=False)
print("U px_best by sole_host:", u.groupby("sole_host").px_best.median().round(0).to_dict(), " D_eq:", u.groupby("sole_host").D_eq_km.median().round(2).to_dict())
print("H in set B px_best median:", round(b[b['class']=='H'].px_best.median()), " U median px:", round(u.px_best.median()))
