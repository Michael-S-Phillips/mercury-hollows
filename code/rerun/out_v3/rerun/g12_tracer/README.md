# g12: host-crater tracing page (6 Nov 2026)

Purpose: replace DEM-only terrain classes with outlines traced on the MDIS 166 m mosaic for the hosts that hold most
terrain-test hollows. Curvature v3 (g10_threshold/curv_classes_v3.py) supplies the starting outlines.

Hosts: the top 64 of 321 by n_peak + n_wall + n_floor + n_other in g11_curvature/terrain/out_v3/rows_curv_n24_h100_s0.03.csv.
They hold 3,011 of 3,760 hollows (80%); the top 30 hold 2,388 (64%).

Files
- build_tracer_assets.py (device): per crater, assets/img/<uid>.jpg (BDR 166 m resampled to a local azimuthal-equidistant grid,
  half-width 1.3 R), assets/dem/<uid>.png (665 m DEM on the same grid, z - zmin in R*256+G, B = 255 where valid),
  assets/meta/<uid>.json (v3 peak/ring mask polygons from DEM pixel edges, v3 wall-foot radius per 2 deg, hollows, masked craters).
- build_tracer_data.py (cloud): simplifies the v3 peak/ring polygons (0.5 km), drops parts and holes < 2 km2, writes tracer_data.json.
- tracer_page.html: the published page (https://claude.ai/artifact/Hnvf7kvMZBqPpkaKvRX6So). Traces save to the artifact's
  database, collection "hosts", one document per crater uid: status (todo/editing/done/skip), note, shapes[{kind peak|hole|foot|rim,
  origin v3|drawn, edited, lon[], lat[]}], updated. Longitudes are continuous about the crater center (may fall outside 0-360).

Units from traces (planned rerun): peak = inside peak outlines minus hole outlines; wall = wall foot to rim crest (1.05 R if no rim
traced), outside peak; floor = inside wall foot, outside peak. Traced hosts use these units; the other hosts keep curvature v3.
Then: v3 error against the tracings; terrain test, positions and slope tests, knob comparison, Figures S2 and 2d.

## Version 3 of the page (6 Nov 2026)
- Freehand drawing (version 2): in Draw mode, drag draws an outline and releasing closes it; in Reshape mode, a dragged line replaces the
  shorter stretch of the selected outline between its ends; Space or right-drag pans. "Clear peak outlines" removes peak and ring outlines.
- Slump (version 3): a fifth outline type, kind "slump", for wall material that has slid down the wall or out onto the floor. It may cross
  the wall foot; where it overlaps wall or floor it counts as slump. A second per-crater status, pass2 (todo/done, "Wall foot and slump:
  Not checked / Checked"), records which craters have been revisited for wall-foot cleanup and slump, so a crater checked with no slump
  outline means no slump rather than not looked at. The host list dot shows pass2.
- Planned analysis once the pass is done (revised after Michael's comments: gas could diffuse through a whole slump body, and charging
  from buried melt and freshness from slopes can both act): classes peak, slump, wall, floor (slump over wall and floor); slump hollow
  density against intact wall and floor at matched slope (Mantel-Haenszel by slope bin, as the peak/wall slope control), so slope
  and material are separated; slump bodies (one outline per body) that reach the floor level compared with bodies that stay on the
  wall; and one rate model with terrain class and slope together, so both effects are estimated rather than one chosen.
  The 665 m DEM slope is coarser than most slump scarps, which limits the slope term.

## Version 4 of the page (7 Nov 2026)
- Base image choice (Display > Base image): monochrome 166 m (img/), MD3 color 665 m (col/), enhanced color 665 m (ecol/), or none.
  Color images from build_color_assets.py (GIS/raster MD3Color and EnhancedColor global mosaics, lon_0 = 0, bilinear onto the tracer
  grid at 333 m, the three bands stretched together between their joint 1st and 99.5th percentiles). Added at Michael's request for
  craters with poor lighting and long shadows in the monochrome mosaic. Alignment checked on Sholem Aleichem and Eminescu.
