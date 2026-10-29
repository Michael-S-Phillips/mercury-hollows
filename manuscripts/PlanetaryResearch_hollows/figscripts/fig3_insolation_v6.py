#!/usr/bin/env python3
"""Figure 3 of manuscript v6 (and Figure 3 of the formation-model draft): panels a-b of figures/fig2_insolation_v3.png,
the top 1,060 pixel rows of that 400 dpi figure. fig2_insolation_v3.png is make_figs_pr.py's 'fig2_insolation' output
(the v3 revision with the Delta-T odds ratio in panel b)."""
import os
from PIL import Image
F = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "figures")
im = Image.open(os.path.join(F, "fig2_insolation_v3.png"))
im.crop((0, 0, im.width, 1060)).save(os.path.join(F, "fig3_insolation_v6.png"), dpi=(400, 400))
print("wrote fig3_insolation_v6.png", im.width, "x 1060")
