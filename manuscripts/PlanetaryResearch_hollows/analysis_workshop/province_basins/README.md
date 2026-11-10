# The hollow-rich province NW of Caloris (5 Nov 2026; workshop analysis, not in the paper)

**Catalog note.** `GRL_insolation_hollows/data/host_assignments.csv` holds only the 5,409 hollows that lie within 1.1 radii of a host crater >= 5 km. The full catalog of 6,667 non-trace hollows is `code/rerun/out_v2/hollows_canonical.csv`. Of the 1,258 hollows outside host craters, 765 are at 45-60N. This region has more of them than anywhere else.
- Crater-level hosting tests (h09) use host craters by design: `province_basins.py`, `b59_profile.py`, `caloris_sectors.py`, `basin_rank.py`.
- Maps and field-level tests now use the full catalog: `fig_fassett_region.py`, `hrp_edge_test.py` and `field_alignment.py` (v2), and `band_fit.py`.

**Results** (records in the matching *.json files):
- **b59 (Orgel et al. 2020).** Hosting is enhanced inside b59, by a factor of 6.6 at 0.5-1 radius, and at 1-1.5 radii. Among 25 testable old basins, b59 has the strongest signal. b63 shows none.
- **Caloris exterior** (rim to ~1,000 km beyond it). Hosting is enhanced only to the NW (8.0x), E (6.0x) and N (3.6x).
- **Band of hollows.** A NW-SE band runs from b59's NE rim toward Caloris: 257 hollows in 16 fields, about 600 km long and 40 km wide. 90% of these hollows are outside host craters. The band runs within ~21 deg of radial to Caloris, and its great circle passes 470 km from the Caloris center, inside the basin. It traces the SW arc of candidate basin A (possiblebasins.gpkg), so A's "rim" is probably this band.
- **Field orientation.** Individual fields (36 elongated) do not align with the band or with the radial direction from Caloris. They tend slightly concentric to Caloris.
- **Higher-albedo plains of Fassett et al. 2009 (Fig. 9).** They released no map files, so I built an automated bright/red mask from MDIS MDR color. It captures the NW arc inside b59 but misses the southern patches.
  - Within 1.2 b59 radii, 23 of 61 fields are on the mask or within 25 km of it, against 22% of the area (near-edge p = 0.03).
  - A hand-traced map from Michael would replace the mask.

## Tectonic map check (5 Nov 2026): tectonic_lines.py, fig_tectonic_lines.png
Klimczak et al. (2025) global tectonic map (GIS/vector/klimczak_tectonic_map_lines.gpkg), frame 38-62N, 98-145E.
- In this frame the map has 185 contractional landforms (8,939 km), 6 short extensional landforms (33 km) and no troughs.
- The contractional landforms run concentric to Caloris rather than radial: 44.5% of their length lies within 30 deg of concentric and 25.2% within 30 deg of radial (33% each if orientations were uniform).
- Hollow fields are no closer to them than the area is: 37 of 114 fields (32.5%) are within 25 km, against 39.7% of the area (p = 0.95). Within 10 km the figures are 14 of 114 (12.3%) against 15.7% of the area.
- The large Caloris-radial ridges (for example the bright ridge at ~45-49N, 135-145E, which runs parallel to the hollows of the NW-SE band) are not in the map. They are Van Eyck lineated terrain (impact sculpture), not tectonic landforms, so a test of hollows against them needs new tracing.
- Literature (searched 5 Nov): no published test of hollows against the Caloris radial lineations. Bickel 2025 and Deutsch 2025 do not address it. The closest are Thomas et al. 2014a (two groups of isolated hollows NW of Caloris, and pitted ground near the Timgad/Paestum Valles convergence at ~59N 122E, as summarized by Blewett et al. 2018) and Wright et al. 2020 (hollows on Odin knobs).

## Tracing page (5 Nov 2026): tracer/
build_tracer_assets.py writes the basemaps and overlays.json for the Caloris Lineation Tracer artifact (https://claude.ai/artifact/8o9BmuvhhKs9Vg87keMgkC). Traces save to its database (collection "traces": kind, note, lon[], lat[], created). Image mapping: x = (lon - 98)/47 * W, y = (62 - lat)/24 * H.

## Knobs (5 Nov 2026): knob_depths.py, wright2020_knobs/
Wright et al. (2020) knob shapefile, unzipped from Downloads into wright2020_knobs/ (read_dbf.py reads the attribute table).
- Of the 6,667 cataloged hollows, 27 lie on knobs, on 21 of the 2,011 knobs, and 18 of those 27 fall within host craters.
- Blewett depths on knobs: 2, both from one image at 62.2N on the catalog's northern edge. These are not usable.
- Automatic shadow depths: see ../depth_expansion_pilot/knob_measure.py and its README section.
