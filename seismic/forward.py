"""2D Acoustic Wave Forward Solver (finite-difference, explicit).

Solves the scalar wave equation:
    1/c(x,z)^2 * ∂²p/∂t² = ∂²p/∂x² + ∂²p/∂z²  + s(t)δ(x-xs,z-zs)

Using 2nd-order FD in space and time with optional simple absorbing boundaries.
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple


@dataclass
class AcousticModel:
    """2D acoustic velocity model."""
    vp: np.ndarray       # (nz, nx) P-wave velocity [m/s]
    dx: float            # grid spacing in x [m]
    dz: float            # grid spacing in z [m]
    nx: int
    nz: int


@dataclass
class Source:
    """Seismic source."""
    ix: int              # x index
    iz: int              # z index
    wavelet: np.ndarray  # source time function (nt,)
    dt: float            # time step [s]


@dataclass
class Receiver:
    """Seismic receiver."""
    ix: int
    iz: int


def ricker_wavelet(f0: float, dt: float, nt: int, t0: Optional[float] = None) -> np.ndarray:
    """Generate a Ricker wavelet.

    Args:
        f0: dominant frequency [Hz]
        dt: time step [s]
        nt: number of time steps
        t0: delay time (default: 1.2/f0)
    """
    if t0 is None:
        t0 = 1.2 / f0
    t = np.arange(nt) * dt
    a = (np.pi * f0 * (t - t0)) ** 2
    return (1.0 - 2.0 * a) * np.exp(-a)


def solve_acoustic_2d(
    model: AcousticModel,
    source: Source,
    receivers: list,
    nt: int,
    abc_width: int = 20,
) -> Tuple[np.ndarray, Optional[np.ndarray]]:
    """Solve 2D acoustic wave equation with finite differences.

    Args:
        model: velocity model
        source: source object
        receivers: list of Receiver objects
        nt: number of time steps
        abc_width: absorbing boundary width (simple exponential damping)

    Returns:
        seismogram: (n_receivers, nt) pressure at receivers
        snapshots: None (set to save memory)
    """
    nx, nz = model.nx, model.nz
    dx, dz = model.dx, model.dz
    dt = source.dt
    vp = model.vp

    # CFL check
    cfl = np.max(vp) * dt * np.sqrt(1.0 / dx**2 + 1.0 / dz**2)
    if cfl > 1.0:
        raise ValueError(f"CFL condition violated: {cfl:.3f} > 1.0")

    # Precompute velocity squared * dt^2 / h^2
    c2dt2_dx2 = (vp * dt / dx) ** 2
    c2dt2_dz2 = (vp * dt / dz) ** 2

    # Absorbing boundary damping
    damp = np.ones((nz, nx), dtype=np.float64)
    if abc_width > 0:
        for i in range(abc_width):
            val = np.exp(-((abc_width - i) * 0.015) ** 2)
            # left
            damp[:, i] = np.minimum(damp[:, i], val)
            # right
            damp[:, nx - 1 - i] = np.minimum(damp[:, nx - 1 - i], val)
            # top
            damp[i, :] = np.minimum(damp[i, :], val)
            # bottom
            damp[nz - 1 - i, :] = np.minimum(damp[nz - 1 - i, :], val)

    # Pressure fields: current, previous
    p_cur = np.zeros((nz, nx), dtype=np.float64)
    p_prev = np.zeros((nz, nx), dtype=np.float64)

    # Output seismogram
    n_rec = len(receivers)
    seismogram = np.zeros((n_rec, nt), dtype=np.float64)

    for it in range(nt):
        # Laplacian (2nd order FD)
        d2p_dx2 = np.zeros_like(p_cur)
        d2p_dz2 = np.zeros_like(p_cur)

        d2p_dx2[:, 1:-1] = p_cur[:, 2:] - 2.0 * p_cur[:, 1:-1] + p_cur[:, :-2]
        d2p_dz2[1:-1, :] = p_cur[2:, :] - 2.0 * p_cur[1:-1, :] + p_cur[:-2, :]

        # Time stepping
        p_next = (2.0 * p_cur - p_prev
                  + c2dt2_dx2 * d2p_dx2
                  + c2dt2_dz2 * d2p_dz2)

        # Add source
        if it < len(source.wavelet):
            p_next[source.iz, source.ix] += (vp[source.iz, source.ix] ** 2
                                              * dt ** 2 * source.wavelet[it])

        # Apply absorbing boundary
        p_next *= damp

        # Record at receivers
        for ir, rec in enumerate(receivers):
            seismogram[ir, it] = p_next[rec.iz, rec.ix]

        # Swap
        p_prev = p_cur
        p_cur = p_next

    return seismogram, None
