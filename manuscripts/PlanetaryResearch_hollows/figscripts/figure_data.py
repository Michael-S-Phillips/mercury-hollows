"""Write the data drawn in a Matplotlib figure to CSV (one file per figure), for the 'tables behind each figure'.
Added for the archive (29 Oct 2026). Called just before a figure is saved; it reads the plotted artists and does not change
the figure. Columns: panel (axes in drawing order, a, b, c, ...), kind (line, points, bar, segment), label, x, y, and x2, y2
for segments (error bars and intervals). Images (maps) are not written; their source grids are in the repository."""
import csv, os, string
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.collections import PathCollection, LineCollection
from matplotlib.patches import Rectangle

def dump(fig, name, outdir):
    os.makedirs(outdir, exist_ok=True)
    rows = []
    axes = [a for a in fig.axes if a.get_visible() and a.get_label() != "<colorbar>"]
    for k, ax in enumerate(axes):
        panel = string.ascii_lowercase[k] if k < 26 else str(k)
        for ln in ax.get_lines():
            x, y = np.asarray(ln.get_xdata(), float), np.asarray(ln.get_ydata(), float)
            if x.size == 0 or not ln.get_visible(): continue
            lab = "" if str(ln.get_label()).startswith("_") else ln.get_label()
            rows += [(panel, "line", lab, xi, yi, "", "") for xi, yi in zip(x, y)]
        for c in ax.collections:
            lab = "" if str(c.get_label()).startswith("_") else c.get_label()
            if isinstance(c, PathCollection):
                off = c.get_offsets()
                rows += [(panel, "points", lab, float(p[0]), float(p[1]), "", "") for p in np.asarray(off)]
            elif isinstance(c, LineCollection):
                for seg in c.get_segments():
                    if len(seg) >= 2: rows.append((panel, "segment", lab, float(seg[0][0]), float(seg[0][1]), float(seg[-1][0]), float(seg[-1][1])))
        for p in ax.patches:
            if isinstance(p, Rectangle) and p.get_width() and p.get_height():
                lab = "" if str(p.get_label()).startswith("_") else p.get_label()
                rows.append((panel, "bar", lab, p.get_x() + p.get_width() / 2, p.get_y() + p.get_height(), p.get_y(), p.get_width()))
    path = os.path.join(outdir, f"{name}_data.csv")
    with open(path, "w", newline="") as f:
        w = csv.writer(f); w.writerow(["panel", "kind", "label", "x", "y", "x2_or_bar_bottom", "y2_or_bar_width"]); w.writerows(rows)
    return path
