"""Ising problem definition and energy computation."""

import numpy as np
import torch
from dataclasses import dataclass
from typing import Optional


@dataclass
class IsingResult:
    """Result of a CIM solver run.

    Attributes:
        best_spins: spin configuration with the lowest energy found (N,)
        best_energy: the lowest Ising energy found
        spin_trajectories: spin amplitudes over time (batch, N, T)
        energy_evolution: Ising energy over time (batch, T)
        all_best_spins: best spins per run (batch, N)
        all_best_energies: best energy per run (batch,)
        wall_time: wall-clock time in seconds
    """
    best_spins: np.ndarray
    best_energy: float
    spin_trajectories: np.ndarray
    energy_evolution: np.ndarray
    all_best_spins: np.ndarray
    all_best_energies: np.ndarray
    wall_time: float


def ising_energy(J: np.ndarray, spins: np.ndarray,
                 h: Optional[np.ndarray] = None) -> float:
    """Compute Ising energy H = -1/2 * s^T J s - h^T s."""
    E = -0.5 * spins @ J @ spins
    if h is not None:
        E -= h @ spins
    return float(E)


def ising_energy_batch(J: torch.Tensor, sig: torch.Tensor,
                       h: Optional[torch.Tensor] = None) -> torch.Tensor:
    """Batch Ising energy: (batch,) from J (N,N), sig (batch, N)."""
    batch_size, N = sig.shape
    Js = sig @ J  # (batch, N)
    E = -0.5 * (sig * Js).sum(dim=1)  # -1/2 * s^T J s
    if h is not None:
        E -= sig @ h
    return E


def validate_ising(J: np.ndarray, h: Optional[np.ndarray] = None):
    """Validate Ising problem inputs."""
    assert J.ndim == 2 and J.shape[0] == J.shape[1], "J must be square"
    assert np.allclose(J, J.T), "J must be symmetric"
    if h is not None:
        assert h.ndim == 1 and h.shape[0] == J.shape[0], \
            "h must be 1-D with same size as J"
