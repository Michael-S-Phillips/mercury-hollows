#!/usr/bin/env python3
"""Insert dated correction notes (2 Nov 2026) after each statement of the pre-fix vent result in five docs.
Backups in doc_backups/. Each anchor must occur exactly once; the note goes after the end of the anchor's block
(paragraph or list item), as a blockquote indented to the list item when inside one."""
import os, re, sys
D = os.path.expanduser("~/mnt/hollows-emissivity/docs")
FULL = ("**Correction (2 Nov 2026):** the vent catalog's 'Lat' and 'Lon' headers are swapped; read correctly "
        "(lat = 'Lon', lon = 'Lat'), hollows cluster near vents (6.6% of hollows within 25 km of a vent vs 0.5% for "
        "longitude-rotated positions, z = 8.3) and craters with vents host hollows 4.88 (2.89–8.25) times as often at "
        "fixed size, latitude, completeness, plains, ΔT_w and regional rate. See runs/vent_association_corrected.json "
        "and manuscripts/PlanetaryResearch_hollows/analysis_workshop/test8_vents.py.")
def SHORT(where):
    return ("**Correction (2 Nov 2026):** the vent result cited here predates the vent header fix. Read with the "
            "'Lat' and 'Lon' headers swapped back, hollows cluster near vents (z = 8.3 at 25 km under the rotation "
            f"null); see the correction note {where}.")
E = {
 "geologic_plausibility.md": [
  ("though it does not exclude it for individual vent-adjacent fields.", FULL),
  ("LRM-raster significances are sensitive to the null choice; worth revisiting those on the rotation",
   SHORT("in §1.2. The vent signal does not vanish under the rotation null, so this flag now applies to the LRM raster only")),
  ("a slowly-accumulated dispersed trace and not a fumarolic point source (consistent with the vent", SHORT("in §1.2")),
  ("5 km) and with the *absence* of near-field vent association (§1.2), which together point at", SHORT("in §1.2")),
  ("remains the least-proven step (the vent test rules out the SCArFS fumarolic-redeposition version", SHORT("in §1.2")),
  ("proximity) is contradicted by §1.2.", SHORT("in §1.2")),
  ("and the absence of near-field vent association together point at impact-fracture access rather than", SHORT("in §1.2")),
 ],
 "STUDY_STATE.md": [
  ("10 km) but **no near-field vent association** at any radius from 25–200 km.",
   SHORT("under 'Spatial associations, and a null-model caution' below")),
  ("support, and one the crater/vent statistics happen to favor.",
   SHORT("under 'Spatial associations, and a null-model caution' below")),
  ("- **No near-field vent association** at any radius from 25–200 km (at 25 km the observed fraction,", FULL),
  ("Together the crater and vent results point at **impact-fracture access** to a volatile reservoir",
   SHORT("on the vent bullet above")),
  ("The vent and LRM-raster signals near-vanish only once clustering is preserved, so this affects the",
   SHORT("on the vent bullet above. The vent signal does not vanish under the rotation null, so this caution now applies to the LRM raster only")),
 ],
 "nonvolatile_alternatives.md": [
  ("hollow–crater association found here (5 km z=3.4, 10 km z=4.0) and with the absence of near-field vent",
   SHORT("on the 'Vents' bullet in §4")),
  ("- **Vents — clean negative (PI result, inherited).**", FULL),
  ("radii — the PI's methodological point, confirmed. This is the geologic counterpart of the vent negative: hollows",
   SHORT("on the 'Vents' bullet above")),
  ("the pore-space concentration requirement (§5), the crater association, the vent non-association, and the Wang",
   SHORT("on the 'Vents' bullet above")),
  ("budget cannot separate, the second consistent with the crater association and vent non-association found here.",
   SHORT("on the 'Vents' bullet in §4")),
 ],
 "poleward_limit_and_depth_discriminant.md": [
  ("geologic track found significant and with the *absence* of near-field vent association (together", FULL),
 ],
 "lag_and_morphology.md": [
  ("crater association the geologic track found (5 km z=3.4, 10 km z=4.0) and the *absence* of near-field", FULL),
  ("fracture-fed branch aligns with the geologic track's significant crater association and the absent",
   SHORT("in §2.4")),
 ],
}
LIST = re.compile(r"^(\s*)([-*]|\d+\.)\s+")
def block_bounds(lines, L):
    """start and end line of the paragraph or list item containing line L"""
    s = L
    while s > 0 and lines[s - 1].strip() and not LIST.match(lines[s]) and not lines[s - 1].startswith("#"):
        s -= 1
    e = L
    while e + 1 < len(lines) and lines[e + 1].strip() and not LIST.match(lines[e + 1]) and not lines[e + 1].startswith("#"):
        e += 1
    return s, e
dry = "--dry" in sys.argv
for fn, items in E.items():
    p = f"{D}/{fn}"; txt = open(p, encoding="utf-8").read()
    if "Correction (2 Nov 2026)" in txt: print("already annotated, skipping", fn); continue
    lines = txt.split("\n"); ins = []
    for anchor, note in items:
        assert txt.count(anchor) == 1, (fn, anchor, txt.count(anchor))
        L = txt[:txt.index(anchor)].count("\n") + anchor.count("\n")
        s, e = block_bounds(lines, L)
        m = LIST.match(lines[s])
        if m:   # inside a list item: blockquote indented to the item's content
            ind = " " * len(m.group(0)); new = [f"{ind}> {note}"]
        else:
            new = ["", f"> {note}"]
        ins.append((e, new, lines[s][:60], lines[e][-60:]))
    for e, new, a, b in sorted(ins, key=lambda t: -t[0]):
        lines[e + 1:e + 1] = new
        print(f"{fn}: after line {e + 1} | block '{a}...' ends '...{b}'")
    if not dry:
        open(p, "w", encoding="utf-8").write("\n".join(lines))
