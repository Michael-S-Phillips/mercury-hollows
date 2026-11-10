# Why the released Bickel et al. (2025) coordinates are offset (4 Nov 2026)

The released latitudes and longitudes are reproduced to within 0–2 m (all 6 frames whose PDS labels are in
data/mdis_frames/IMG) by a linear mapping. The normalized image pixel coordinates are stretched onto the
latitude–longitude bounding box of the four corner (reticle) points in the PDS label, as if the image edges ran
north–south and east–west (bbox_test.py). The mapping is oriented one way in 412 of 933 images and turned 180° in 501; 19 are mirrored in one axis.

- The label corner points agree with the camera-model corners to 0.1–0.5 km, so the inputs were accurate. The
  error comes from ignoring the image's rotation relative to north (reticle_label.py).
- Consequence: released = recomputed rotated by the image's angle from north (90° − NORTH_AZIMUTH, modulo 180°)
  and stretched by |cos a| + |sin a| about the image center. Across 933 images with ≥ 6 detections, a rotation plus
  uniform stretch fits to a median of 0.14 km, against 3.5 km for a shift alone and 1.8 km for a rotation alone
  (released_error.py, rot_vs_north.py).
- Figure S1 rows: rotation 16.2°, 23.9°, 35.7° (90° − north azimuth: 16.3°, 22.9°, 35.8°); stretch 1.24, 1.31, 1.39
  (predicted 1.24, 1.31, 1.40).

## What the release left out (orientation_rule.py)
- **The rule.** The direction of each axis is set per axis by comparing corner coordinates: column 0 is west unless the
  left corners lie east of the right corners, and row 0 is north unless the top corners lie south of the bottom
  corners. This matches the best-fitting orientation in 927 of 933 images. Columns always map to longitude and rows to
  latitude, so images turned near 90° from north get the largest errors, and a few are mirrored.
- **What is missing.** The rotation of the image relative to north (NORTH_AZIMUTH and TWIST_ANGLE are in every label).
  This is equivalent to a north-up geotransform with no rotation terms.
- **What the same label corners would have given.** Used as the corners of the tilted frame, they give coordinates
  within 0.06–0.62 km (median per frame) of the camera-model coordinates in the 6 labeled frames. The release is off by
  2.5–22.7 km in the same frames.
- Not for the paper (Michael, 4 Nov 2026): documented here for reference only.
