#!/usr/bin/env python3
"""Plot a receiver-gather trace (.dat): rows = positions, cols = time steps."""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fname = sys.argv[1] if len(sys.argv) > 1 else "trace_zaxis_txx.dat"
data = np.loadtxt(fname)          # shape (npos, ntime)
npos, ntime = data.shape

vmax = np.percentile(np.abs(data), 99.5)
if vmax == 0:
    vmax = np.abs(data).max() or 1.0

fig, ax = plt.subplots(figsize=(9, 5))
im = ax.imshow(data, cmap="seismic", vmin=-vmax, vmax=vmax,
               origin="upper", aspect="auto",
               extent=[0, ntime, npos, 0])
ax.set_title(f"{fname}   ({npos} pos x {ntime} steps)")
ax.set_xlabel("time step")
ax.set_ylabel("z position (grid)")
fig.colorbar(im, ax=ax, label="amplitude")

out = fname.rsplit(".", 1)[0] + ".png"
fig.tight_layout()
fig.savefig(out, dpi=150)
print("saved", out, "shape", data.shape, "vmax", vmax)
