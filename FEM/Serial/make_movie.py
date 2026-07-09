#!/usr/bin/env python3
"""Build an MP4 movie of wavefield snapshots from it0000 to it4050 (step 50).

Usage: python3 make_movie.py [field]   e.g. field = txx | txz | tzz | vx | vz
"""
import sys
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

FIELD = sys.argv[1] if len(sys.argv) > 1 else "txx"
PREFIX = f"snap_{FIELD}_it"
START, STOP, STEP = 0, 4050, 50

# collect existing files in iteration order
files = []
for it in range(START, STOP + 1, STEP):
    f = f"{PREFIX}{it:04d}.dat"
    if glob.glob(f):
        files.append((it, f))
print(f"found {len(files)} frames")

# load all frames
frames = [(it, np.loadtxt(f)) for it, f in files]

# fixed symmetric color scale. Pass an explicit vmax as the 2nd CLI arg
# (e.g. `make_movie.py txx 1e-4`); otherwise use a robust global percentile
# so relative amplitudes stay comparable across time.
if len(sys.argv) > 2:
    vmax = float(sys.argv[2])
else:
    allabs = np.concatenate([np.abs(d).ravel() for _, d in frames])
    vmax = np.percentile(allabs, 99.5)
print(f"global vmax = {vmax:.3e}")

fig, ax = plt.subplots(figsize=(7, 6))
it0, d0 = frames[0]
im = ax.imshow(d0, cmap="seismic", vmin=-vmax, vmax=vmax,
               origin="upper", aspect="equal", animated=True)
ax.set_xlabel("x (grid)")
ax.set_ylabel("z (grid)")
cb = fig.colorbar(im, ax=ax, label=f"{FIELD} amplitude")
title = ax.set_title(f"{FIELD}  it={it0:04d}")
fig.tight_layout()

def update(idx):
    it, d = frames[idx]
    im.set_array(d)
    title.set_text(f"{FIELD}  it={it:04d}")
    return im, title

ani = animation.FuncAnimation(fig, update, frames=len(frames),
                              interval=100, blit=False)

out = f"snap_{FIELD}_movie.mp4"
writer = animation.FFMpegWriter(fps=10, bitrate=2400)
ani.save(out, writer=writer, dpi=130)
print("saved", out)
