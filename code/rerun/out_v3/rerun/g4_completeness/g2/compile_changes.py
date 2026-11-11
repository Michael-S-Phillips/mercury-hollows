# g4 copy of code/rerun/out_v3/rerun/g2_insolation/compile_changes.py (6 Nov 2026). Changed: folder depth 5 -> 6; G4_MAP (env) maps the
#!/usr/bin/env python3
# columns rec / v3 / v3b to folders of out/ (e.g. rec=v3,v3=v3r,v3b=v3r compares v3r with old-grid v3); output file from G4_OUTCSV.
"""g2 (6 Nov 2026), new script: collect every catalog-dependent number of the g2 sections (main text 3.2 plains part and 3.5;
supplement S3, S5, S6, S8 and Tables S5-S6; Figures 1b and 3; Table 1 and the Discussion) for the stored v2 record and the
v3 / v3b reruns, and write changes_g2.csv.
v2_record comes from the stored record of the original script; where no record was stored (plains_hosts, wall fractions of
Table S6, LR p values of Table S5), it is the v2 rerun of the copied script, which reproduces the paper (see
out/v2_reproduction.json and README.md).
conclusion_change: yes when, relative to v2, a ratio or odds ratio crosses 1 or its 95% interval newly includes or
excludes 1 (0 for an amplitude), a p value crosses 0.05, or an AIC ranking changes (sign of an AIC difference flips);
the bracket names the catalog(s) where this happens."""
import os, json, csv, math
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, *[".."] * 6))
MAP = dict(x.split("=") for x in os.environ.get("G4_MAP", "").split(",") if x)   # g4: catalog column -> out/<folder>
PR = f"{ROOT}/manuscripts/PlanetaryResearch_hollows"; W6 = f"{PR}/analysis_workshop/insolation_v6"; GR = f"{ROOT}/manuscripts/insolation_thermal/records"
REC = {"latcorr_density_record.json": f"{GR}/latcorr_density_record.json", "latcorr_inc60_75_record.json": f"{GR}/latcorr_inc60_75_record.json",
       "dT_plains_record.json": f"{PR}/analysis_v3/dT_plains_record.json", "dT_record.json": f"{PR}/analysis_v3/dT_record.json",
       "dT_hosting_completeness_record.json": f"{PR}/analysis_v7/dT_hosting_completeness_record.json",
       "transition_layer_hilat_record.json": f"{PR}/analysis_v3/transition_layer_hilat_record.json",
       "hotwarm_record.json": f"{ROOT}/code/rerun/latitude/hotwarm/hotwarm_record.json",
       "aspect_record.json": f"{ROOT}/code/rerun/latitude/aspect/aspect_record.json"}
for n in ("cos2lon_fit_lat15", "size_quirk", "band_crater_hotwarm", "hollow_rich_hosts", "threshold", "stage_test", "stage_test_lowinc"):
    REC[f"{n}.json"] = f"{W6}/{n}.json"
_cache = {}
def J(cat, name):
    key = (cat, name)
    if key not in _cache:
        src = MAP.get(cat, cat)
        if src == "rec":
            if name == "wall_fractions.json": d = json.load(open(f"{GR}/figure_numbers_v2.json"))["fig2c_wall"]; d["n_wall_hollows_all"] = json.load(open(REC["aspect_record.json"]))["n_wall_hollows"]
            elif name in REC: d = json.load(open(REC[name]))
            else: d = json.load(open(f"{HERE}/out/v2/{name}"))      # no stored record: v2 rerun
        else: d = json.load(open(f"{HERE}/out/{src}/{name}"))
        _cache[key] = d
    return _cache[key]
CATS = ("rec", "v3", "v3b")
def f2(x, d=2): return f"{x:.{d}f}"
def ci(x, lo, hi, d=2): return f"{x:.{d}f} ({lo:.{d}f}-{hi:.{d}f})"
def pf(p): return f"{p:.2g}" if p >= 1e-3 else f"{p:.1e}"
rows = []
def flag(kind, vals):
    """vals: dict cat -> tuple. kind: 'ci1' (est, lo, hi vs 1), 'ci0' (vs 0), 'p' (p), 'sign' (difference), None."""
    def state(v):
        if kind == "ci1": return (v[1] > 1) - (v[2] < 1)          # +1 excludes 1 above, -1 below, 0 includes
        if kind == "ci0": return (v[1] > 0) - (v[2] < 0)
        if kind == "p": return v[0] < 0.05
        if kind == "sign": return v[0] > 0
        if kind == "pt1": return v[0] > 1
    if kind is None: return "no"
    s0 = state(vals["rec"]); ch = [c for c in ("v3", "v3b") if state(vals[c]) != s0]
    return f"yes ({', '.join(ch)})" if ch else "no"
def add(section, quantity, text, fmt, kind, script, get):
    vals = {c: get(c) for c in CATS}
    rows.append(dict(section=section, quantity=quantity, text_value=text, v2_record=fmt(vals["rec"]), v3=fmt(vals["v3"]),
                     v3b=fmt(vals["v3b"]), script=script, conclusion_change=flag(kind, vals)))
def manual(section, quantity, text, v2, v3, v3b, script, change):
    rows.append(dict(section=section, quantity=quantity, text_value=text, v2_record=v2, v3=v3, v3b=v3b, script=script, conclusion_change=change))

# ------------------------------------------------------------------ plains (3.2, S6, Figure 1b)
LD, LI = "latcorr_density_record.json", "latcorr_inc60_75_record.json"; SL = "latcorr.py + latcorr3.py"
for band, txt in (("40-60", "0.41"), ("60-80", "0.13")):
    add("3.2; S6", f"Field density inside / outside the plains, {band}N, cells <= 40 m/px", txt, lambda v: f2(v[0]), None, SL,
        lambda c, b=band: (J(c, LD)[f"N_le40_in_{b}"]["d"] / J(c, LD)[f"N_le40_out_{b}"]["d"],))
for zone, band, txt in (("in", "40-60", "36 (25-50)"), ("in", "60-80", "12 (7-19)"), ("out", "40-60", "88 (76-101)"), ("out", "60-80", "91 (72-112)"), ("out", "0-30", "61 (54-69)")):
    add("S6; Fig 1b", f"Field density {'inside' if zone == 'in' else 'outside'} the plains, {band}N, per 1e6 km2 (95% Poisson); fields", txt,
        lambda v: f"{v[0]:.0f} ({v[1]:.0f}-{v[2]:.0f}); k={v[3]}", None, SL,
        lambda c, z=zone, b=band: (J(c, LD)[f"N_le40_{z}_{b}"]["d"], J(c, LD)[f"N_le40_{z}_{b}"]["lo"], J(c, LD)[f"N_le40_{z}_{b}"]["hi"], J(c, LD)[f"N_le40_{z}_{b}"]["k"]))
add("S6", "Odds ratio of hosting (any hollow within 1.1 R) on vs off the plains, craters >= 10 km at 40-80N imaged <= 80 m/px (1,179 craters); hosts", "0.42 (0.24-0.74)",
    lambda v: f"{ci(*v[:3])}; p={pf(v[3])}; hosts={v[4]} ({v[5]} on plains)", "ci1", "dT_plains.py",
    lambda c: (J(c, "dT_plains_record.json")["dT"]["OR"], *J(c, "dT_plains_record.json")["dT"]["ci"], J(c, "dT_plains_record.json")["dT"]["p"], J(c, "dT_plains_record.json")["dT"]["hosts"], J(c, "dT_plains_record.json")["dT"]["hosts_nsp"]))
add("S6", "Craters on the plains with hollows within 0.9 R (of 1,849 craters >= 5 km)", "21", lambda v: str(v[0]), None, "plains_hosts.py", lambda c: (J(c, "plains_hosts.json")["plains_h09"],))
add("S6", "Of those, number >= 10 km across", "20 of 21", lambda v: f"{v[0]} of {v[1]}", None, "plains_hosts.py", lambda c: (J(c, "plains_hosts.json")["plains_h09_ge10km"], J(c, "plains_hosts.json")["plains_h09"]))
add("S6", "Median diameter of plains craters with / without hollows (km)", "32 vs 10", lambda v: f"{v[0]:.0f} vs {v[1]:.0f}", None, "plains_hosts.py",
    lambda c: (J(c, "plains_hosts.json")["median_D_plains_h09"], J(c, "plains_hosts.json")["median_D_plains_no_h09"]))
add("S6", "Outside the plains, density at 60-80N / 40-60N: all images; without the 90-150E province; 60-75 deg incidence only (range quoted)", "0.77-1.03",
    lambda v: f"{min(v[0], v[1], v[2]):.2f}-{max(v[0], v[1], v[2]):.2f} (all {v[0]:.2f} [{v[3][0]:.2f}-{v[3][1]:.2f}]; no province {v[1]:.2f}; incidence {v[2]:.2f})", None, "latcorr3.py; latcorr5.py",
    lambda c: (J(c, LD)["rr_60-80_vs_40-60_out"]["rr"], J(c, LD)["rr_60-80_vs_40-60_out_noprov"]["rr"], J(c, LI)["le40_out_rr"]["rr"], J(c, LD)["rr_60-80_vs_40-60_out"]["ci"]))
for key, band, lab in (("N_le40_out", "30-40", "outside, all images"), ("N_le40_out_40-60_noprov", None, "outside, 90-150E province removed, 40-60N"),
                       ("N_le40_out_60-80_noprov", None, "outside, 90-150E province removed, 60-80N")):
    add("Fig 1b", f"Field density {lab}{', ' + band + 'N' if band else ''} (plotted, not quoted)", "plotted",
        lambda v: f"{v[0]:.0f} ({v[1]:.0f}-{v[2]:.0f}); k={v[3]}", None, "latcorr3.py",
        lambda c, k=key, b=band: tuple((J(c, LD)[k][b] if b else J(c, LD)[k])[x] for x in ("d", "lo", "hi", "k")))
for band in ("0-30", "30-40", "40-60", "60-80"):
    add("Fig 1b", f"Field density outside the plains, 60-75 deg incidence images only, {band}N (plotted, not quoted)", "plotted",
        lambda v: f"{v[0]:.0f} ({v[1]:.0f}-{v[2]:.0f}); k={v[3]}", None, "latcorr5.py",
        lambda c, b=band: tuple(J(c, LI)["le40_out"][b][x] for x in ("d", "lo", "hi", "k")))

# ------------------------------------------------------------------ hot longitudes: fields (3.5, S8, Table S6, Figure 3a)
CL = "cos2lon_fit_lat15.json"
add("3.5; S8; Fig 3a", "Hollow fields within 15 deg of the equator, hot-to-warm ratio of observed/expected (Caloris < 1,200 km excluded); fields", "1.95 (1.59-2.47); 379 fields",
    lambda v: f"{ci(*v[:3])}; {v[3]} fields", "ci1", "cos2lon_lat15.py",
    lambda c: (J(c, CL)["lat15"]["sector_hot_over_warm"], *J(c, CL)["lat15"]["sector_hot_over_warm_ci"], J(c, CL)["lat15"]["n"]))
add("Fig 3a; Table S6", "Amplitude b of a + b cos 2 lon, fields within 15 deg; intercept a", "0.67 (0.55-0.78); a = 0.99",
    lambda v: f"{ci(*v[:3])}; a = {v[3]:.2f}", "ci0", "cos2lon_lat15.py",
    lambda c: (J(c, CL)["lat15"]["slope"], *J(c, CL)["lat15"]["slope_ci"], J(c, CL)["lat15"]["intercept"]))
HW = "hotwarm_record.json"
def hwb(c, var, hemi, band):
    r = next(x for x in J(c, HW)["variants"][var][hemi]["bands"] if x["band"] == band)
    return (r["hot_over_warm"], r["hw_ci"][0], r["hw_ci"][1], r["perm_p"], r["n"])
for sec, var, hemi, band, txt in (("Table S6", "caloris_excl", "north", "0-15", "1.85 (1.40-2.49); 197 fields"), ("Table S6", "caloris_excl", "south", "0-15", "2.10 (1.56-2.94); 182 fields"),
                                  ("Table S6", "caloris_excl", "north", "15-30", "1.57 (1.18-2.07); 199 fields"), ("Table S6", "caloris_excl", "north", "30-45", "0.75 (0.53-1.06); 140 fields"),
                                  ("3.5; S8; Table S6", "caloris_excl", "north", "45-60", "0.44 (0.29-0.61); p = 0.0005; 147 fields"),
                                  ("S8; Table S6", "caloris_excl", "north", "60-75", "0.65 (0.42-1.01); 91 fields"),
                                  ("S8; Table S6", "caloris_and_province_excl", "north", "45-60", "0.66 (0.42-1.00); p = 0.065; 86 fields"),
                                  ("S8; Table S6", "caloris_and_province_excl", "north", "60-75", "3.1 (1.6-9.8); 37 fields")):
    lab = f"Hot-to-warm ratio of hollow fields, {band}{'N' if hemi == 'north' else 'S'}{', without the 90-150E province' if 'province' in var else ''} (Caloris excluded); p; fields"
    add(sec, lab, txt, lambda v: f"{ci(*v[:3])}; p={v[3]:.2g}; {v[4]} fields", "ci1", "hotwarm.py", lambda c, a=var, h=hemi, b=band: hwb(c, a, h, b))
    _ci, _p = flag("ci1", {c: hwb(c, var, hemi, band) for c in CATS}), flag("p", {c: (hwb(c, var, hemi, band)[3],) for c in CATS})
    _ch = sorted({x.strip() for f_ in (_ci, _p) if f_.startswith("yes") for x in f_[5:-1].split(",")})
    rows[-1]["conclusion_change"] = f"yes ({', '.join(_ch)})" if _ch else "no"
add("3.5; S8", "Share of the 45-60N field reversal (log scale) that the 90-150E province accounts for", "about half",
    lambda v: f"{v[0]:.2f}", None, "hotwarm.py", lambda c: (1 - math.log(hwb(c, "caloris_and_province_excl", "north", "45-60")[0]) / math.log(hwb(c, "caloris_excl", "north", "45-60")[0]),))
def hh(c, hemi, band):
    r = next(x for x in J(c, HW)["individual_hollows"][hemi]["bands"] if x["band"] == band); return (r["hot_over_warm"], r["hw_ci"][0], r["hw_ci"][1], r["n"])
add("Table S6", "Hot-to-warm ratio of individual hollows, 0-15N", "1.08 (0.95-1.24)", lambda v: f"{ci(*v[:3])}; {v[3]} hollows", "ci1", "hotwarm.py", lambda c: hh(c, "north", "0-15"))
add("Table S6", "Hot-to-warm ratio of individual hollows, 0-15S", "5.4 (4.3-7.1)", lambda v: f"{ci(*v[:3], d=1)}; {v[3]} hollows", "ci1", "hotwarm.py", lambda c: hh(c, "south", "0-15"))
add("Table S6; S8", "Hot-to-warm ratio of individual hollows, 15-75N, four 15-deg bands (range)", "0.37-0.57",
    lambda v: f"{min(v):.2f}-{max(v):.2f}", "pt1", "hotwarm.py", lambda c: tuple(hh(c, "north", b)[0] for b in ("15-30", "30-45", "45-60", "60-75")))
add("S3", "Share of detected fields expected in the hot sectors under a uniform distribution, 15-deg northern bands; without fields near Caloris", "0.46-0.59 (0.43-0.52)",
    lambda v: f"{v[0]:.2f}-{v[1]:.2f} ({v[2]:.2f}-{v[3]:.2f})", None, "hotwarm.py",
    lambda c: (min(r["null_hot_frac"] for r in J(c, HW)["variants"]["all"]["north"]["bands"] if r["null_hot_frac"] is not None),
               max(r["null_hot_frac"] for r in J(c, HW)["variants"]["all"]["north"]["bands"] if r["null_hot_frac"] is not None),
               min(r["null_hot_frac"] for r in J(c, HW)["variants"]["caloris_excl"]["north"]["bands"] if r["null_hot_frac"] is not None),
               max(r["null_hot_frac"] for r in J(c, HW)["variants"]["caloris_excl"]["north"]["bands"] if r["null_hot_frac"] is not None)))
SQ = "size_quirk.json"
for k, lab, txt in (("N1_C1_weighted_crater_counts", "completeness-weighted crater counts", "1.71 (1.36-2.16); p < 0.001"),
                    ("N2_C1_weighted_crater_area", "completeness-weighted crater area", "1.84 (1.50-2.33); p < 0.001"),
                    ("N3_predicted_hosting_no_longitude", "hosting model without longitude", "2.08 (1.69-2.59); p < 0.001")):
    add("S8", f"Field hot-to-warm ratio within 15 deg with expectation from {lab}", txt, lambda v: f"{ci(*v[:3])}; p={v[3]:.2g}", "ci1", "size_quirk.py",
        lambda c, kk=k: (J(c, SQ)["ii_fields"]["nulls"][kk]["hot_over_warm"], *J(c, SQ)["ii_fields"]["nulls"][kk]["ci"], J(c, SQ)["ii_fields"]["nulls"][kk]["p"]))

# ------------------------------------------------------------------ hot longitudes: craters (3.5, S8, Table 1, Discussion)
add("3.5; S8", "Odds ratio of hosting (h09), hot vs warm sector, |lat| < 15, completeness and diameter controlled; h09 hosts", "2.2 (1.6-3.0) / 2.19 (1.63-2.96)",
    lambda v: f"{ci(*v[:3])}; {v[3]} hosts", "ci1", "size_quirk.py",
    lambda c: (J(c, SQ)["iii_crater_logit_lat15"]["all"]["all_sizes"]["OR"], J(c, SQ)["iii_crater_logit_lat15"]["all"]["all_sizes"]["lo"], J(c, SQ)["iii_crater_logit_lat15"]["all"]["all_sizes"]["hi"],
               J(c, "band_crater_hotwarm.json")["-15..15"]["h09"]))
add("3.5; S8", "Hot-sector x log D interaction within 15 deg (p)", "0.22", lambda v: f"{v[0]:.2f}", "p", "size_quirk.py",
    lambda c: (J(c, SQ)["iii_crater_logit_lat15"]["all"]["hot_x_logD"]["p"],))
BC = "band_crater_hotwarm.json"
for band, txt, sec in (("15..30", "3.2 (2.1-4.8) / 3.18 (2.10-4.80)", "3.5; S8"), ("30..45", "1.8 (1.1-2.8) / 1.76 (1.12-2.77)", "3.5; S8"),
                       ("45..60", "1.72 (0.95-3.13)", "S8"), ("60..75", "1.48 (0.66-3.33)", "S8")):
    add(sec, f"Odds ratio of hosting (h09), hot vs warm sector, {band.replace('..', '-')}N, completeness and diameter controlled; p; h09 hosts", txt,
        lambda v: f"{ci(*v[:3])}; p={v[3]:.2g}; {v[4]} hosts", "ci1", "band_crater_hotwarm.py", lambda c, b=band: (*J(c, BC)[b]["C1"], J(c, BC)[b]["h09"]))
add("S8; 3.5; Discussion 4.1", "Delta T_w x log D interaction (p) at 15-30N and 30-45N, the bands the size trend came from", "0.01-0.03",
    lambda v: f"{v[0]:.2g} and {v[1]:.2g}", None, "band_crater_hotwarm.py", lambda c: (J(c, BC)["15..30"]["dTxlogD"][3], J(c, BC)["30..45"]["dTxlogD"][3]))
rows[-1]["conclusion_change"] = "yes (v3, v3b)" if all(J(c, BC)[b]["dTxlogD"][3] >= 0.05 for c in ("v3", "v3b") for b in ("15..30", "30..45")) else "check"
HR = "hollow_rich_hosts.json"
add("3.5; S8", "Hollows within 0.9 R of h09 hosts at 15-75N near the warm / hot longitudes", "2,388 / 1,169", lambda v: f"{v[0]:,} / {v[1]:,}", None, "hollow_rich_hosts.py", lambda c: (J(c, HR)["hollows_warm"], J(c, HR)["hollows_hot"]))
add("S8", "h09 host craters at 15-75N near the hot / warm longitudes", "188 / 131", lambda v: f"{v[0]} / {v[1]}", None, "hollow_rich_hosts.py", lambda c: (J(c, HR)["hosts_hot"], J(c, HR)["hosts_warm"]))
add("3.5; S8", "Ten hosts with most hollows: share of the hollows at 15-75N; number near the warm longitudes", "49%; 8 of 10",
    lambda v: f"{100 * v[0]:.0f}%; {v[1]} of 10", None, "hollow_rich_hosts.py", lambda c: (J(c, HR)["top10_share_of_hollows"], sum(1 - x["hot"] for x in J(c, HR)["top10"])))
add("3.5; S8", "Without those ten hosts: hollows near the hot / warm longitudes", "1,026 / 792", lambda v: f"{v[0]:,} / {v[1]:,}", None, "hollow_rich_hosts.py",
    lambda c: (J(c, HR)["without_top10"]["hollows_hot"], J(c, HR)["without_top10"]["hollows_warm"]))

# ------------------------------------------------------------------ Delta T_w fits (3.5, S3, Table S6, Figure 3b)
DT, DC = "dT_record.json", "dT_hosting_completeness_record.json"
add("S3; 3.5; Fig 3b; Table S6", "Craters >= 5 km within 60 deg (any-hollow hosts within 1.1 R)", "27,271 (657 hosts)", lambda v: f"{v[0]:,} ({v[1]} hosts)", None, "dT_hosting.py",
    lambda c: (J(c, DT)["n_craters"], J(c, DT)["n_hosts"]))
add("S3", "Host craters among the 31,571 craters >= 5 km (any hollow within 1.1 R)", "707", lambda v: f"{v[0]}", None, "transition_layer_hilat.py",
    lambda c: (J(c, "transition_layer_hilat_record.json")["abslat_le90"]["hosts"],))
add("3.5; S3; Table S6", "Odds ratio per 10 K of Delta T_w, 10 km craters (no completeness control)", "1.60 (1.43-1.78)", lambda v: ci(*v), "ci1", "dT_hosting.py",
    lambda c: (J(c, DT)["dT_model"]["OR10K_10km"], *J(c, DT)["dT_model"]["OR10K_10km_ci"]))
add("3.5; S3; Table S6", "Odds ratio per 10 K of Delta T_w, 100 km craters (no completeness control)", "1.18 (1.08-1.29)", lambda v: ci(*v), "ci1", "dT_hosting.py",
    lambda c: (J(c, DT)["dT_model"]["OR10K_100km"], *J(c, DT)["dT_model"]["OR10K_100km_ci"]))
add("3.5; S3; Table S6", "Diameter at which the odds ratio reaches 1 (km)", "near 360", lambda v: f"{v[0]:,.0f}", None, "dT_hosting.py", lambda c: (J(c, DT)["dT_model"]["zero_km"],))
rows[-1]["conclusion_change"] = "yes (v3, v3b): beyond the largest crater; the interaction is not significant"
add("3.5; S3; Fig 3b", "Delta T_w x log D interaction p (no completeness control)", "< 0.0001 / 1.5e-5", lambda v: pf(v[0]), "p", "dT_hosting.py", lambda c: (J(c, DT)["dT_model"]["p"]["dT:logD"],))
add("S3", "AIC, Delta T_w model vs distance from the nearest warm longitude (difference = distance - Delta T_w)", "5,295 vs 5,333",
    lambda v: f"{v[0]:,.0f} vs {v[1]:,.0f} (diff {v[1] - v[0]:+.0f})", None, "dT_hosting.py", lambda c: (J(c, DT)["dT_model"]["aic"], J(c, DT)["hotness_model"]["aic"]))
rows[-1]["conclusion_change"] = flag("sign", {c: (J(c, DT)["hotness_model"]["aic"] - J(c, DT)["dT_model"]["aic"],) for c in CATS})
add("3.5; S3; Table S6", "Odds ratio per 10 K, 10 km craters, completeness controlled", "1.64 (1.46-1.83)", lambda v: ci(*v), "ci1", "dT_hosting_completeness.py",
    lambda c: tuple(J(c, DC)["with_C1"]["OR10K_10km"]))
add("3.5; S3; Table S6", "Odds ratio per 10 K, 100 km craters, completeness controlled", "1.26 (1.14-1.38)", lambda v: ci(*v), "ci1", "dT_hosting_completeness.py",
    lambda c: tuple(J(c, DC)["with_C1"]["OR10K_100km"]))
add("3.5; S3; Table S6", "Diameter at which the odds ratio reaches 1, completeness controlled (km)", "near 720", lambda v: f"{v[0]:,.0f}", None, "dT_hosting_completeness.py", lambda c: (J(c, DC)["with_C1"]["zero_km"],))
rows[-1]["conclusion_change"] = "yes (v3, v3b): beyond the largest crater; the interaction is not significant"
add("3.5; S3", "Delta T_w x log D interaction p, completeness controlled", "0.0003", lambda v: pf(v[0]), "p", "dT_hosting_completeness.py", lambda c: (J(c, DC)["with_C1"]["p_interaction"],))
add("S3", "AIC of the completeness-controlled fit", "4,658", lambda v: f"{v[0]:,.0f}", None, "dT_hosting_completeness.py", lambda c: (J(c, DC)["with_C1"]["aic"],))
add("Fig 3b", "Host counts in the six size bins 5-10, 10-20, 20-40, 40-80, 80-160, >= 160 km (plotted labels)", "130, 112, 116, 112, 122, 63",
    lambda v: ", ".join(str(x) for x in v), None, "dT_hosting.py", lambda c: tuple(b["hosts"] for b in J(c, DT)["bins"]))
add("Fig 3b", "Odds ratio per 10 K in the six size bins (plotted points)", "plotted",
    lambda v: ", ".join(f"{x:.2f}" for x in v), None, "dT_hosting.py", lambda c: tuple(b["OR10K"] for b in J(c, DT)["bins"]))
add("Fig 3b", "Odds ratio per 10 K in the 5-10 km bin (plotted point)", "plotted", lambda v: ci(*v), "ci1", "dT_hosting.py",
    lambda c: (J(c, DT)["bins"][0]["OR10K"], *J(c, DT)["bins"][0]["ci"]))

TH = "threshold.json"
def T(c, s="lat60_main"): return J(c, TH)[s]
add("S8", "Containment-test craters within 60 deg (h09 hosts)", "27,257 (573)", lambda v: f"{v[0]:,} ({v[1]})", None, "threshold.py", lambda c: (T(c)["n"], T(c)["hosts"]))
add("3.5; S8", "Odds ratio per 10 K of Delta T_w, h09, completeness controlled, within 60 deg", "1.43 (1.32-1.54)", lambda v: ci(*v), "ci1", "threshold.py", lambda c: tuple(T(c)["linear_OR_per10K"]))
add("S8", "AIC of restricted cubic splines (3-5 df) minus linear Delta T_w", "+1.8 to +5.5",
    lambda v: f"{min(v):+.1f} to {max(v):+.1f}", None, "threshold.py", lambda c: tuple(T(c)["AIC"][f"spline_dT_df{k}"] - T(c)["AIC"]["linear_dT"] for k in (3, 4, 5)))
rows[-1]["conclusion_change"] = flag("sign", {c: (min(T(c)["AIC"][f"spline_dT_df{k}"] - T(c)["AIC"]["linear_dT"] for k in (3, 4, 5)),) for c in CATS})
add("3.5; S8", "Best single step in Delta T_w: threshold (K) and AIC minus linear (incl. +2 for choosing t)", "41 K; +15.6",
    lambda v: f"{v[0]:.0f} K; {v[1]:+.1f}", None, "threshold.py", lambda c: (T(c)["step_dT_t"], T(c)["AIC"]["step_dT"] - T(c)["AIC"]["linear_dT"]))
rows[-1]["conclusion_change"] = flag("sign", {c: (T(c)["AIC"]["step_dT"] - T(c)["AIC"]["linear_dT"],) for c in CATS})
for s, lab, txt in (("lat60_main", "within 60 deg", "linear p = 0.0001; step p = 0.16"), ("global_main", "all latitudes", "step not significant (p = 0.13)")):
    add("S8", f"Step and linear term together, {lab}: p of the linear term; p of the step", txt, lambda v: f"linear p={v[0]:.2g}; step p={v[1]:.2g}", "p", "threshold.py",
        lambda c, ss=s: (T(c, ss)["step_plus_linear"]["step_p"], T(c, ss)["step_plus_linear"]["dT10_p"]))
    rows[-1]["v2_record"], rows[-1]["v3"], rows[-1]["v3b"] = [f"linear p={T(c, s)['step_plus_linear']['dT10_p']:.2g}; step p={T(c, s)['step_plus_linear']['step_p']:.2g}" for c in CATS]
    lin = {c: T(c, s)["step_plus_linear"]["dT10_p"] < 0.05 for c in CATS}; stp = {c: T(c, s)["step_plus_linear"]["step_p"] < 0.05 for c in CATS}
    ch = [c for c in ("v3", "v3b") if lin[c] != lin["rec"] or stp[c] != stp["rec"]]; rows[-1]["conclusion_change"] = f"yes ({', '.join(ch)})" if ch else "no"
add("S8", "AIC of a cos 2 lon term in place of Delta T_w, minus linear Delta T_w", "+22", lambda v: f"{v[0]:+.1f}", "sign", "threshold.py", lambda c: (T(c)["AIC"]["cos2lon"] - T(c)["AIC"]["linear_dT"],))

# ------------------------------------------------------------------ Stage 1 (3.5, S8, Discussion)
ST, SL7 = "stage_test.json", "stage_test_lowinc.json"
add("S8", "Hollows with stage labels (Stage 1 / Stage 2)", "2,468 (1,500 / 968)", lambda v: f"{v[0]:,} ({v[1]:,} / {v[2]})", None, "stage_test.py",
    lambda c: (J(c, ST)["n_labeled_stage1_2"], J(c, ST)["n_stage1"], J(c, ST)["n_stage2"]))
add("S8", "Odds of Stage 1 per 10 K of Delta T_w, all classifications", "1.10 (0.97-1.26)", lambda v: ci(*v), "ci1", "stage_test.py",
    lambda c: tuple(J(c, ST)["results"]["all"]["dT10"][k] for k in ("OR", "lo", "hi")))
L7 = lambda c: J(c, SL7)["inc_lt70"]
add("S8", "Hollows classified in images at incidence < 70 deg (|lat| <= 60)", "1,647", lambda v: f"{v[0]:,}", None, "stage_test_lowinc.py", lambda c: (L7(c)["n"],))
add("3.5; S8; Discussion 4.1", "Stage 1 fraction where Delta T_w < ~45 K (bins up to 50 K) / > 50 K", "56-71% / 96-100%",
    lambda v: f"{100 * v[0]:.0f}-{100 * v[1]:.0f}% / {100 * v[2]:.0f}-{100 * v[3]:.0f}%", None, "stage_test_lowinc.py",
    lambda c: (min(x["s1"] for x in L7(c)["s1_share_by_dT"][:4]), max(x["s1"] for x in L7(c)["s1_share_by_dT"][:4]),
               min(x["s1"] for x in L7(c)["s1_share_by_dT"][4:]), max(x["s1"] for x in L7(c)["s1_share_by_dT"][4:])))
add("3.5; S8", "Step in Stage 1 odds: best threshold; odds ratio (clustered by field)", "44 K; 3.0 (1.5-6.2)", lambda v: f"{v[0]:.0f} K; {ci(*v[1:])}", "ci1", "stage_test_lowinc.py",
    lambda c: (L7(c)["best_t"], L7(c)["step"]["OR"], L7(c)["step"]["lo"], L7(c)["step"]["hi"]))
rows[-1]["conclusion_change"] = flag("ci1", {c: (L7(c)["step"]["OR"], L7(c)["step"]["lo"], L7(c)["step"]["hi"]) for c in CATS})
add("3.5; S8", "AIC of linear minus step model for Stage 1", "27", lambda v: f"{v[0]:.1f}", "sign", "stage_test_lowinc.py", lambda c: (L7(c)["AIC"]["linear"] - L7(c)["AIC"]["step_best"],))
add("S8", "Step odds ratio without the Caloris region", "2.9", lambda v: ci(*v), "ci1", "stage_test_lowinc.py", lambda c: tuple(L7(c)["step_noCaloris"][k] for k in ("OR", "lo", "hi")))

# ------------------------------------------------------------------ walls (3.5, S8, Table S6, Table 1)
WF = "wall_fractions.json"
add("S8", "Northern+southern wall hollows (0.55-1.05 R of hosts >= 15 km)", "2,401", lambda v: f"{v[0]:,}", None, "aspect_stats.py; wall_fractions.py", lambda c: (J(c, WF)["n_wall_hollows_all"],))
for b, txt in (("0-15", "0.48 (0.42-0.54); 0.57 (0.46-0.67); 267 hollows, 63 craters"), ("15-30", "0.39 (0.36-0.42); 0.53 (0.42-0.64); 874, 60"),
               ("15-30_noCaloris", "0.45 (0.41-0.48); 0.54 (0.42-0.66); 747, 53"), ("30-45", "0.62 (0.58-0.65); 0.49 (0.39-0.60); 621, 60"),
               ("45-60", "0.71 (0.66-0.75); 0.51 (0.38-0.65); 330, 38"), ("60-75", "0.75 (0.64-0.86); 0.71 (0.54-0.88); 56, 26")):
    add("3.5; S8; Table S6; Table 1" if b == "60-75" else "S8; Table S6", f"Fraction of northern wall hollows on equator-facing walls, {b.replace('_noCaloris', ', Caloris excluded')}N: per hollow; per host crater; n", txt,
        lambda v: f"{ci(*v[:3])}; {ci(*v[3:6])}; {v[6]}, {v[7]}", None, "wall_fractions.py",
        lambda c, bb=b: tuple(J(c, WF)[bb][k] for k in ("f", "lo", "hi", "fc", "clo", "chi", "n", "nc")))
    v = {c: tuple(J(c, WF)[b][k] for k in ("fc", "clo", "chi")) for c in CATS}
    st = lambda x: (x[1] > .5) - (x[2] < .5); ch = [c for c in ("v3", "v3b") if st(v[c]) != st(v["rec"])]
    vh = {c: tuple(J(c, WF)[b][k] for k in ("f", "lo", "hi")) for c in CATS}; ch += [c for c in ("v3", "v3b") if st(vh[c]) != st(vh["rec"]) and c not in ch]
    rows[-1]["conclusion_change"] = f"yes ({', '.join(ch)})" if ch else "no"
manual("S8; Table S6", "Equator-facing / pole-facing wall hollow density, 352 host craters 30-250 km", "1.52 (0.95-2.47)", "1.52 (0.95-2.47)", "not rerun", "not rerun",
       "terrain analysis of Section S4 / Table S2 (another group)", "n/a")

# ------------------------------------------------------------------ transition layer (S5, Table S5, Discussion)
TL = "transition_layer_hilat_record.json"
from scipy.stats import chi2
def tl(c, k):
    r = J(c, TL)[k]; a = {m: r["fits"][m]["aic"] for m in ("dT", "L368", "both")}
    return a, chi2.sf(a["dT"] - a["both"] + 4, 2), chi2.sf(a["L368"] - a["both"] + 4, 2), r["hosts"], r["corr"]
for k, lab, txt, ptxt in (("abslat_le60", "|lat| <= 60", "5,295.0 / 5,297.7 / 5,291.8; 657 hosts", "0.03 / 0.007"),
                          ("abslat_le75", "|lat| <= 75", "5,715.2 / 5,728.7 / 5,715.2; 701 hosts", "0.14 / 0.0002"),
                          ("abslat_le90", "all latitudes", "5,789.6 / 5,803.9 / 5,792.5; 707 hosts", "0.58 / 0.0005")):
    add("S5; Table S5" + ("; Discussion 4.1" if k != "abslat_le75" else ""), f"AIC of hosting fits, {lab}: Delta T_w / L368 / both (L368 - Delta T_w); hosts", txt,
        lambda v: f"{v[0]['dT']:,.1f} / {v[0]['L368']:,.1f} / {v[0]['both']:,.1f} ({v[0]['L368'] - v[0]['dT']:+.1f}); {v[3]} hosts", None, "transition_layer_hilat.py", lambda c, kk=k: tl(c, kk))
    rows[-1]["conclusion_change"] = flag("sign", {c: (tl(c, k)[0]["L368"] - tl(c, k)[0]["dT"],) for c in CATS})
    add("S5; Table S5", f"LR test p, {lab}: layer added to Delta T_w / Delta T_w added to layer (2 df)", ptxt,
        lambda v: f"{v[1]:.2g} / {v[2]:.2g}", None, "transition_layer_hilat.py (p from the AICs, compile_changes.py)", lambda c, kk=k: tl(c, kk))
    ch = [c for c in ("v3", "v3b") if (tl(c, k)[1] < .05) != (tl("rec", k)[1] < .05) or (tl(c, k)[2] < .05) != (tl("rec", k)[2] < .05)]
    rows[-1]["conclusion_change"] = f"yes ({', '.join(ch)})" if ch else "no"

# ------------------------------------------------------------------ qualitative statements (Table 1, Discussion, Abstract)
bc = lambda c, b: J(c, BC)[b]["C1"]
def excl(v): return v[1] > 1
manual("Table 1; Discussion 4.1; 3.5", "Craters at the hot longitudes host hollows more often 'near the equator and to 45N' (crater-level odds ratio excludes 1 in each band)",
       "near the equator and to 45N",
       "; ".join(f"{b.replace('..', '-')}: {'yes' if excl(bc('rec', b)) else 'no'}" for b in ("-15..15", "15..30", "30..45")),
       "; ".join(f"{b.replace('..', '-')}: {'yes' if excl(bc('v3', b)) else 'no'}" for b in ("-15..15", "15..30", "30..45")),
       "; ".join(f"{b.replace('..', '-')}: {'yes' if excl(bc('v3b', b)) else 'no'}" for b in ("-15..15", "15..30", "30..45")),
       "band_crater_hotwarm.py; size_quirk.py", "yes (v3, v3b)" if not excl(bc("v3", "30..45")) and not excl(bc("v3b", "30..45")) else "check")
manual("3.5; Discussion 4.1", "The rise of hosting odds with Delta T_w weakens with crater size (Delta T_w x log D interaction p, without / with completeness)",
       "p < 0.0001 / 0.0003", *[f"{pf(J(c, DT)['dT_model']['p']['dT:logD'])} / {pf(J(c, DC)['with_C1']['p_interaction'])}" for c in CATS], "dT_hosting.py; dT_hosting_completeness.py",
       "yes (v3, v3b)" if all(J(c, DT)['dT_model']['p']['dT:logD'] > .05 and J(c, DC)['with_C1']['p_interaction'] > .05 for c in ("v3", "v3b")) else "check")
manual("Abstract; 3.5", "Hollows are more common at the hot longitudes near the equator (field ratio within 15 deg; crater odds ratio within 15 deg)",
       "qualitative", *[f"fields {J(c, CL)['lat15']['sector_hot_over_warm']:.2f} ({J(c, CL)['lat15']['sector_hot_over_warm_ci'][0]:.2f}-{J(c, CL)['lat15']['sector_hot_over_warm_ci'][1]:.2f}); craters {bc(c, '-15..15')[0]:.2f} ({bc(c, '-15..15')[1]:.2f}-{bc(c, '-15..15')[2]:.2f})" for c in CATS],
       "cos2lon_lat15.py; band_crater_hotwarm.py", "yes (v3b; field level only, crater level unchanged)")

with open(f"{HERE}/" + os.environ.get("G4_OUTCSV", "changes_g2_passA.csv"), "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["section", "quantity", "text_value", "v2_record", "v3", "v3b", "script", "conclusion_change"]); w.writeheader(); w.writerows(rows)
print(len(rows), "rows")
for r in rows: print(f"{r['conclusion_change']:15s} | {r['quantity'][:70]:70s} | {r['text_value'][:28]:28s} | {r['v2_record'][:40]:40s} | {r['v3'][:40]:40s} | {r['v3b'][:40]}")
