"""CIM-based seismic velocity inversion via Travel-Time Tomography.

PENALTY-FREE formulation: Binary velocity choice per region.
  - Each region chooses between two velocity candidates (fast/slow)
  - One Ising spin per region: s_i = +1 → fast, s_i = -1 → slow
  - Travel time: T = L @ slowness(s), linear in spins!
  - Misfit ||T_obs - T_syn||² maps EXACTLY to Ising H = -1/2 s^T J s - h^T s
  - NO penalty terms, NO constraint violations, CIM works on pure physics!
"""

import time
import numpy as np
from itertools import product as iterproduct
from typing import Tuple, List

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from cimsim.solver import solve as cim_solve
from cimsim.ising import ising_energy


# ===================================================================
#  Ray geometry
# ===================================================================

def compute_ray_lengths(
    nz: int, nx: int, dz: float, dx: float,
    region_map: np.ndarray,
    source_pos: Tuple[float, float],
    receiver_positions: List[Tuple[float, float]],
    n_regions: int,
) -> np.ndarray:
    """Compute straight-ray path lengths through each region.

    Returns: L (n_receivers, n_regions)
    """
    n_rec = len(receiver_positions)
    L = np.zeros((n_rec, n_regions), dtype=np.float64)

    sx, sz = source_pos
    for ir, (rx, rz) in enumerate(receiver_positions):
        n_segments = 1000
        t_vals = np.linspace(0, 1, n_segments + 1)
        x_pts = sx + t_vals * (rx - sx)
        z_pts = sz + t_vals * (rz - sz)
        ray_length = np.sqrt((rx - sx) ** 2 + (rz - sz) ** 2)
        ds = ray_length / n_segments

        x_mid = 0.5 * (x_pts[:-1] + x_pts[1:])
        z_mid = 0.5 * (z_pts[:-1] + z_pts[1:])
        ix = np.clip((x_mid / dx).astype(int), 0, nx - 1)
        iz = np.clip((z_mid / dz).astype(int), 0, nz - 1)
        regions_hit = region_map[iz, ix]

        for r in range(n_regions):
            L[ir, r] = np.sum(regions_hit == r) * ds

    return L


def compute_ray_lengths_multi_source(
    nz, nx, dz, dx, region_map,
    source_positions, receiver_positions_list, n_regions,
):
    """Compute ray lengths for multiple source-receiver configurations."""
    all_L = []
    for src_pos, rec_list in zip(source_positions, receiver_positions_list):
        L = compute_ray_lengths(nz, nx, dz, dx, region_map,
                                src_pos, rec_list, n_regions)
        all_L.append(L)
    return np.vstack(all_L)


# ===================================================================
#  Binary velocity QUBO (penalty-free!)
# ===================================================================

def build_binary_qubo(
    L: np.ndarray,
    T_obs: np.ndarray,
    v_slow: np.ndarray,
    v_fast: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Build EXACT Ising Hamiltonian for binary velocity tomography.

    Each region i has two velocity choices:
      s_i = -1 → v_slow[i]  (slowness = 1/v_slow[i])
      s_i = +1 → v_fast[i]  (slowness = 1/v_fast[i])

    Slowness parameterization:
      slow_i(s) = (s_slow_i + s_fast_i)/2 + (s_fast_i - s_slow_i)/2 * s_i
                = mean_i + delta_i * s_i

    Travel time: T_syn = L @ slow(s) = L @ mean + L @ (delta * s)
                       = T_mean + L_delta @ s

    where L_delta[:, i] = L[:, i] * delta_i

    Misfit: f(s) = ||T_obs - T_syn||²
                 = ||T_obs - T_mean - L_delta @ s||²
                 = ||(T_obs - T_mean) - L_delta @ s||²
                 = ||r - L_delta @ s||²

    Let r = T_obs - T_mean, A = L_delta:
    f(s) = ||r - A@s||² = r^T r - 2 r^T A s + s^T A^T A s

    Since s_i² = 1: s^T (A^T A) s = Σ_{i≠j} (A^T A)_{ij} s_i s_j + Σ_i (A^T A)_{ii}

    f(s) = Σ_{i≠j} (A^T A)_{ij} s_i s_j - 2 (A^T r)^T s + const

    CIM: H = -(1/2) s^T J s - h^T s = -(1/2) Σ_{i≠j} J_{ij} s_i s_j - h^T s

    Match: -(1/2) J_{ij} = (A^T A)_{ij}  → J_{ij} = -2 (A^T A)_{ij}  for i≠j
           -h_i = -2 (A^T r)_i             → h_i = 2 (A^T r)_i

    Returns: (J, h) ready for CIM
    """
    n_regions = L.shape[1]

    s_slow = 1.0 / v_slow
    s_fast = 1.0 / v_fast
    mean_s = (s_slow + s_fast) / 2.0
    delta_s = (s_fast - s_slow) / 2.0

    T_mean = L @ mean_s
    r = T_obs - T_mean  # residual

    A = L * delta_s[np.newaxis, :]  # scale columns of L by delta

    AtA = A.T @ A
    Atr = A.T @ r

    J = -2.0 * AtA
    np.fill_diagonal(J, 0.0)
    h = 2.0 * Atr

    return J, h


def decode_spins(spins, v_slow, v_fast):
    """Decode Ising spins to velocities."""
    velocities = np.where(spins > 0, v_fast, v_slow)
    return velocities


def compute_travel_time(L, velocities):
    """Compute travel times given velocity model."""
    return L @ (1.0 / velocities)


# ===================================================================
#  Solvers
# ===================================================================

def cim_inversion(L, T_obs, v_slow, v_fast, **kwargs):
    """CIM-based inversion."""
    import torch
    torch.set_num_threads(4)

    t0 = time.time()
    J, h = build_binary_qubo(L, T_obs, v_slow, v_fast)

    num_timesteps = kwargs.get('num_timesteps', 3000)
    batch_size = kwargs.get('batch_size', 50)
    result = cim_solve(J, h=h, num_timesteps=num_timesteps, batch_size=batch_size)

    # Decode best solution
    best_vel = decode_spins(result.best_spins, v_slow, v_fast)
    T_syn = compute_travel_time(L, best_vel)
    misfit = float(np.sum((T_obs - T_syn) ** 2))

    # Also check other solutions from the batch
    for spins in result.all_best_spins:
        vel = decode_spins(spins, v_slow, v_fast)
        T_syn = compute_travel_time(L, vel)
        m = float(np.sum((T_obs - T_syn) ** 2))
        if m < misfit:
            misfit = m
            best_vel = vel.copy()

    total_time = time.time() - t0
    return best_vel, misfit, total_time


def exhaustive_search(L, T_obs, v_slow, v_fast):
    """Try all 2^N configurations."""
    n_regions = len(v_slow)
    t0 = time.time()
    best_misfit = float('inf')
    best_vel = None

    for bits in range(2 ** n_regions):
        spins = np.array([(bits >> i) & 1 for i in range(n_regions)]) * 2.0 - 1.0
        vel = decode_spins(spins, v_slow, v_fast)
        T_syn = compute_travel_time(L, vel)
        misfit = float(np.sum((T_obs - T_syn) ** 2))
        if misfit < best_misfit:
            best_misfit = misfit
            best_vel = vel.copy()

    return best_vel, best_misfit, time.time() - t0


def random_search(L, T_obs, v_slow, v_fast, n_trials=10000, seed=42):
    """Random search."""
    n_regions = len(v_slow)
    rng = np.random.RandomState(seed)
    t0 = time.time()
    best_misfit = float('inf')
    best_vel = None

    for _ in range(n_trials):
        spins = rng.choice([-1.0, 1.0], size=n_regions)
        vel = decode_spins(spins, v_slow, v_fast)
        T_syn = compute_travel_time(L, vel)
        misfit = float(np.sum((T_obs - T_syn) ** 2))
        if misfit < best_misfit:
            best_misfit = misfit
            best_vel = vel.copy()

    return best_vel, best_misfit, time.time() - t0


def greedy_search(L, T_obs, v_slow, v_fast, n_restarts=20, seed=99):
    """Multi-start greedy coordinate descent."""
    n_regions = len(v_slow)
    rng = np.random.RandomState(seed)
    t0 = time.time()
    best_misfit = float('inf')
    best_vel = None

    for restart in range(n_restarts):
        if restart == 0:
            spins = -np.ones(n_regions)
        else:
            spins = rng.choice([-1.0, 1.0], size=n_regions)

        for _ in range(n_regions * 3):
            improved = False
            for r in range(n_regions):
                cur_vel = decode_spins(spins, v_slow, v_fast)
                cur_misfit = np.sum((T_obs - compute_travel_time(L, cur_vel)) ** 2)
                spins[r] *= -1
                new_vel = decode_spins(spins, v_slow, v_fast)
                new_misfit = np.sum((T_obs - compute_travel_time(L, new_vel)) ** 2)
                if new_misfit < cur_misfit:
                    improved = True
                else:
                    spins[r] *= -1
            if not improved:
                break

        vel = decode_spins(spins, v_slow, v_fast)
        misfit = float(np.sum((T_obs - compute_travel_time(L, vel)) ** 2))
        if misfit < best_misfit:
            best_misfit = misfit
            best_vel = vel.copy()

    return best_vel, best_misfit, time.time() - t0


def simulated_annealing(L, T_obs, v_slow, v_fast, n_steps=50000, seed=123):
    """Simulated annealing baseline."""
    n_regions = len(v_slow)
    rng = np.random.RandomState(seed)
    t0 = time.time()

    spins = rng.choice([-1.0, 1.0], size=n_regions)
    vel = decode_spins(spins, v_slow, v_fast)
    cur_misfit = float(np.sum((T_obs - compute_travel_time(L, vel)) ** 2))
    best_misfit = cur_misfit
    best_vel = vel.copy()

    for step in range(n_steps):
        T = 1.0 * (1.0 - step / n_steps)  # linear cooling
        r = rng.randint(n_regions)
        spins[r] *= -1
        vel = decode_spins(spins, v_slow, v_fast)
        new_misfit = float(np.sum((T_obs - compute_travel_time(L, vel)) ** 2))

        delta = new_misfit - cur_misfit
        if delta < 0 or (T > 0 and rng.rand() < np.exp(-delta / max(T, 1e-10))):
            cur_misfit = new_misfit
            if cur_misfit < best_misfit:
                best_misfit = cur_misfit
                best_vel = vel.copy()
        else:
            spins[r] *= -1

    return best_vel, best_misfit, time.time() - t0
