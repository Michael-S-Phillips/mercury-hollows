# g28_fig2a (10 Nov 2026): Figure 2a on the MDIS mosaic
Michael asked for panel a of Figure 2 (north polar map) with a monochrome basemap and the image footprints on a color ramp.
- `layers.py` writes `_cache/basemap_laea2km.npy` (MDIS 166 m BDR mosaic north of 27 N, reduced to 1.33 km and sampled on a 2 km polar
  Lambert equal-area grid) and `_cache/bestpx_laea2km.npy` (best pixel scale of the covering images, from the coverage grid of the analysis,
  code/rerun/out_v2/coverage_10ppd.npz). `layers_record.json`: 62% of the area north of 30 N is imaged at <= 150 m/pixel, 44% at <= 40.
  The MDIS index on disk has no NAC corner coordinates, so the footprints are those of the coverage grid (latitude-longitude boxes), the same
  coverage the analysis uses. The radial strips are narrow NAC image strips along MESSENGER's near-polar ground tracks.
- `make_figs_pr_v10.py` (copy of g18_traced_r3/terrain/make_figs_pr_v9_v3cat_traced_r3.py) draws Figure 2 with the new panel a: mosaic in gray,
  coverage tinted 60% with a one-hue aqua ramp on a log scale (darker = finer, 10-150 m/pixel), colorbar under the legend, white halos on the
  in-map labels; panels b-d unchanged; Caloris record v3r. Output fig1_exposure_v10.png/.pdf. FIG2A_BLEND=multiply gives the darker variant
  (trial renders in _cache/).
