# bickel2024_hollows_corrected.gpkg

Re-georeferenced Bickel et al. (2024) MDIS hollow detections (0.7-confidence release, 19,110 detections, 4,915 NAC images).

## Why
The released detection coordinates carry per-image georeferencing errors: median 4.4 km, 17% of images > 20 km, 1.3% > 100 km
(a few near the prime meridian are thousands of km off, e.g. CN1005052854M catalogued at lon -159 vs PDS frame centre -0.2).
Error grows with emission angle and pixel scale. See bickel_position_error_record.json and a5_*.png in data/bickel-2024/.

## How
The release's 'Image x/y coord (pix)' columns are exact full-frame normalized coordinates in raw EDR orientation (verified by
template-matching every thumbnail for four images: correlation 0.97-1.00, positions agree to 0.001). For each image:
EDR (PDS MSGRMDS_1001) -> mdis2isis -> spiceinit (ISIS 8.0, reconstructed SPK/CK, USGS MSGR DEM shape model) -> campt at the
bbox centre and 4 bbox corners. All 4,915 images processed; 0 failures after retries (kernel note: msgr_0810_v03 -> symlink to v04).

## Validation
- Recomputed positions lie inside their source frame for 97.7% (the rest are frame-edge points in oblique views where the nadir
  half-diagonal test is too tight); catalogue positions only 91.2%.
- Same hollow seen in two images (< 60 m/px): nearest cross-image detection within 0.5 km for 52% (catalogue 11%);
  median cross-image separation 0.41 km (catalogue 1.90 km).
- Visual: recomputed points sit on the hollow clusters in CN0254138067M / CN0259207642M (bickel_recompute_demo.png).
Residual accuracy is that of uncontrolled reconstructed SPICE: expect ~0.3-1 km against the USGS controlled basemap.

## Layers (CRS: Mercury sphere R = 2439400 m, geographic, east-positive, lon -180..180)
- hollows_points_corrected   : detection centres (camera model). Attributes: lat/lon180 (corrected), lat_orig/lon180_orig (release),
                               shift_east_km/shift_north_km/shift_km, img_x/img_y, frame_samples/lines, confidence, stage, incidence,
                               emission, phase, px_scale_m, thumbnail, in_08_release, georef ('camera_model' for all rows).
- hollows_polygons_corrected : release polygon shape translated by the centre shift (polygon_method attribute).
- hollows_bbox_cameramodel   : detection bounding box, all four corners through the camera model (true footprint).

## Caveats
- Confidence for 187 rows was parsed from a merged 'saturated pixels + confidence' field in the release CSV (last token).
- Coarse-resolution (> 100 m/px), low-incidence detections remain unreliable as hollows (bright ejecta look-alikes) regardless of position.
- The project's earlier spatial analyses used release positions; at 8 ppd (14 km) most are unaffected but 17% of images moved > 20 km,
  so hollows_canonical.csv and anything at 16 ppd should be rebuilt from this file.
Generated 2026-09-11.
