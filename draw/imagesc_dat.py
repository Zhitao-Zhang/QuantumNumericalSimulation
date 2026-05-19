#!/usr/bin/env python3
"""
Read a 2D numeric array from a whitespace-separated text file (e.g. snap_*.dat)
and render it like MATLAB imagesc: pseudocolor plot with auto color limits,
save as PNG.

Dependencies: numpy, matplotlib
  pip install numpy matplotlib
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_grid(path: Path) -> np.ndarray:
    data = np.loadtxt(path, dtype=np.float64)
    if data.ndim != 2:
        raise ValueError(f"expected 2D array after loadtxt, got shape {data.shape}")
    return data


def plot_imagesc(
    data: np.ndarray,
    out_path: Path,
    cmap: str = "viridis",
    vmin: float | None = None,
    vmax: float | None = None,
    dpi: float = 150,
    figsize: tuple[float, float] | None = None,
    aspect: str = "equal",
    colorbar: bool = True,
) -> None:
    # origin='upper': first row at top, like MATLAB imagesc default (YDir reverse).
    if figsize is None:
        figsize = (8, 6 + (0.35 if colorbar else 0))

    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    im = ax.imshow(
        data,
        cmap=cmap,
        aspect=aspect,
        origin="upper",
        interpolation="nearest",
        vmin=vmin,
        vmax=vmax,
    )
    ax.set_xlabel("column index")
    ax.set_ylabel("row index")
    if colorbar:
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Plot a 2D .dat file like MATLAB imagesc and save PNG."
    )
    p.add_argument(
        "input",
        type=Path,
        help="Input file: rows of whitespace-separated numbers",
    )
    p.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output PNG path (default: same basename as input with .png)",
    )
    p.add_argument(
        "--cmap",
        default="viridis",
        help="Matplotlib colormap name (default: viridis). Try jet, RdBu_r, seismic",
    )
    p.add_argument("--vmin", type=float, default=None, help="Fixed color scale minimum")
    p.add_argument("--vmax", type=float, default=None, help="Fixed color scale maximum")
    p.add_argument("--dpi", type=float, default=150, help="Figure DPI (default: 150)")
    p.add_argument(
        "--aspect",
        default="equal",
        choices=("equal", "auto"),
        help="Axes aspect: equal ~ square pixels like axis image (default: equal)",
    )
    p.add_argument("--no-colorbar", action="store_true", help="Omit colorbar")
    args = p.parse_args(argv)

    inp = args.input.expanduser().resolve()
    if not inp.is_file():
        print(f"error: not a file: {inp}", file=sys.stderr)
        return 1

    out = args.output
    if out is None:
        out = inp.with_suffix(".png")
    else:
        out = out.expanduser().resolve()

    try:
        data = load_grid(inp)
    except Exception as e:
        print(f"error reading {inp}: {e}", file=sys.stderr)
        return 1

    plot_imagesc(
        data,
        out,
        cmap=args.cmap,
        vmin=args.vmin,
        vmax=args.vmax,
        dpi=args.dpi,
        aspect=args.aspect,
        colorbar=not args.no_colorbar,
    )
    print(f"wrote {out}  shape={data.shape}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
