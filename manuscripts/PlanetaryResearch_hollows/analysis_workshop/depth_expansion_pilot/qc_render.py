"""Contact-sheet rendering for QC: crop with the Bickel box (orange), sun direction (yellow arrow pointing toward the sun),
the measured shadow run (cyan line from the sunward end to the tip) and a caption (id, auto depth, reference depth, flags)."""
import numpy as np
from PIL import Image, ImageDraw


def tile(img, info, box, meas, caption, size=260, m_extra=0.6):
    H, W = img.shape; x0, x1, y0, y1 = box; bw = max(x1 - x0, y1 - y0); m = max(bw * m_extra, 15)
    X0 = int(max(0, x0 - m)); X1 = int(min(W, x1 + m)); Y0 = int(max(0, y0 - m)); Y1 = int(min(H, y1 + m))
    sub = img[Y0:Y1, X0:X1]; lo, hi = np.percentile(sub, [1, 99.5]); s = np.clip((sub - lo) / (hi - lo + 1e-6), 0, 1)
    z = (size - 40) / max(sub.shape)
    im = Image.fromarray((s * 255).astype(np.uint8)).resize((max(1, int(sub.shape[1] * z)), max(1, int(sub.shape[0] * z))), Image.NEAREST).convert("RGB")
    d = ImageDraw.Draw(im)
    d.rectangle([(x0 - X0) * z, (y0 - Y0) * z, (x1 - X0) * z, (y1 - Y0) * z], outline=(255, 140, 0))
    if meas is not None and np.isfinite(meas.get("run_u0", np.nan)):
        az = np.radians(info["subsol_az"]); du = np.array([-np.cos(az), -np.sin(az)]); dw = np.array([-du[1], du[0]])
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        p0 = np.array([cx, cy]) + meas["run_u0"] * du + meas["run_w"] * dw; p1 = np.array([cx, cy]) + meas["run_u1"] * du + meas["run_w"] * dw
        d.line([(p0[0] - X0) * z, (p0[1] - Y0) * z, (p1[0] - X0) * z, (p1[1] - Y0) * z], fill=(0, 255, 255), width=2)
    T = Image.new("RGB", (size, size)); T.paste(im, (0, 0)); d = ImageDraw.Draw(T)
    a = np.radians(info["subsol_az"]); cx, cy = size - 22, 18
    d.line([cx, cy, cx + 15 * np.cos(a), cy + 15 * np.sin(a)], fill=(255, 255, 0), width=3)
    d.ellipse([cx + 15 * np.cos(a) - 3, cy + 15 * np.sin(a) - 3, cx + 15 * np.cos(a) + 3, cy + 15 * np.sin(a) + 3], fill=(255, 255, 0))
    for i, line in enumerate(caption.split("\n")[:3]):
        d.text((2, size - 38 + 12 * i), line, fill=(0, 255, 0))
    return T


def sheet(tiles, ncol=6, size=260, path="sheet.png"):
    nrow = int(np.ceil(len(tiles) / ncol)); S = Image.new("RGB", (ncol * size, max(1, nrow) * size))
    for i, t in enumerate(tiles): S.paste(t, ((i % ncol) * size, (i // ncol) * size))
    S.save(path); return S.size
