#!/usr/bin/env python3
"""Plot a 2D wavefield snapshot (.dat) as an image."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fname = sys.argv[1] if len(sys.argv) > 1 else "snap_vx_it0800.dat"
data = np.loadtxt(fname)

# symmetric color scale, clipped to enhance weaker phases
vmax = np.percentile(np.abs(data), 99.5)

fig, ax = plt.subplots(figsize=(7, 6))
im = ax.imshow(data, cmap="seismic", vmin=-vmax, vmax=vmax,
               origin="upper", aspect="equal")
ax.set_title(fname)
ax.set_xlabel("x (grid)")
ax.set_ylabel("z (grid)")
fig.colorbar(im, ax=ax, label="amplitude")

out = fname.rsplit(".", 1)[0] + ".png"
fig.tight_layout()
fig.savefig(out, dpi=150)
print("saved", out, "shape", data.shape, "vmax", vmax)
