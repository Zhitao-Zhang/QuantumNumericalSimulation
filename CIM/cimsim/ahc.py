"""Amplitude-Heterogeneity Correction (AHC) CIM solver.

Reference:
    T. Leleu, Y. Yamamoto, P. L. McMahon, and K. Aihara,
    Phys. Rev. Lett. 122, 040607 (2019).
"""

import time
import numpy as np
import torch
from typing import Optional, Callable

from .ising import IsingResult

# Small-matrix CIM loops are bottlenecked by thread-sync at high thread counts
if torch.get_num_threads() > 8:
    torch.set_num_threads(4)


@torch.no_grad()
def solve_ahc(
    J: np.ndarray,
    h: Optional[np.ndarray] = None,
    num_timesteps: int = 1000,
    dt: float = 0.05,
    batch_size: int = 1,
    r: float = 0.2,
    beta: float = 0.05,
    eps: float = 0.07,
    mu: float = 1.0,
    noise: float = 0.0,
    a0: float = 0.2,
    pump_schedule: Optional[Callable] = None,
    feedback_schedule: Optional[Callable] = None,
    device: Optional[torch.device] = None,
    save_trajectories: bool = True,
) -> IsingResult:
    if device is None:
        device = torch.device('cpu')

    J_t = torch.from_numpy(np.ascontiguousarray(J)).float().to(device)
    N = J_t.shape[0]
    has_h = h is not None
    h_t = torch.from_numpy(np.ascontiguousarray(h)).float().to(device) if has_h else None

    ticks = int(num_timesteps / dt)
    EMAX = 32.0

    # --- Build schedules ---
    tick_idx = torch.arange(ticks, device=device)
    r_sched = (pump_schedule(tick_idx * dt) if pump_schedule else
               torch.full((ticks,), r, device=device))
    eps_sched = (feedback_schedule(tick_idx * dt) if feedback_schedule else
                 torch.full((ticks,), eps, device=device))

    # --- Initialize state ---
    x = torch.empty(batch_size, N, device=device).uniform_(-5e-4, 5e-4)
    e = torch.ones(batch_size, N, device=device)

    # Pre-allocate scratch
    MVM = torch.empty_like(x)
    x_sq = torch.empty_like(x)

    best_energy = torch.full((batch_size,), 1e20, device=device)
    sig_opt = torch.zeros(batch_size, N, device=device)

    if save_trajectories:
        spin_traj = torch.zeros(batch_size, N, ticks, device=device)
        energy_traj = torch.zeros(batch_size, ticks, device=device)

    t0 = time.time()

    for t in range(ticks):
        sig = (2.0 * (x > 0) - 1.0).float()

        # Energy with cached sJ
        sJ = sig @ J_t
        H = -0.5 * (sig * sJ).sum(dim=1)
        if has_h:
            H -= sig @ h_t

        if save_trajectories:
            spin_traj[:, :, t] = x
            energy_traj[:, t] = H

        # Pre-pump MVM
        torch.mul(x, x, out=x_sq)
        torch.mm(x, J_t, out=MVM)

        # Pump + saturation
        r_t = r_sched[t]
        x += dt * (x * ((r_t - 1.0) - mu * x_sq))

        # Feedback (old MVM)
        eps_t = eps_sched[t]
        x += (dt * eps_t) * (MVM * e)

        # Noise
        if noise > 0:
            x += eps_t * noise * (torch.rand_like(x) - 0.5)

        # Adaptive target amplitude (from cached sJ)
        delta_a = eps_t * (sJ * sig).mean(dim=1)
        target_a = a0 + delta_a

        # Error variable (new x)
        e += dt * (-beta * (x.square() - target_a[:, None]) * e)
        e.clamp_(max=EMAX)

        # Track best
        improved = H < best_energy
        if improved.any():
            sig_opt[improved] = sig[improved]
            best_energy = torch.minimum(best_energy, H)

    wall_time = time.time() - t0

    # --- Results ---
    sig_opt_np = sig_opt.cpu().numpy()
    all_best_energies = np.empty(batch_size)
    for k in range(batch_size):
        s = sig_opt[k]
        ev = -0.5 * (s @ J_t @ s)
        if has_h:
            ev -= s @ h_t
        all_best_energies[k] = ev.item()

    best_idx = int(np.argmin(all_best_energies))
    return IsingResult(
        best_spins=sig_opt_np[best_idx],
        best_energy=float(all_best_energies[best_idx]),
        spin_trajectories=spin_traj.cpu().numpy() if save_trajectories else None,
        energy_evolution=energy_traj.cpu().numpy() if save_trajectories else None,
        all_best_spins=sig_opt_np,
        all_best_energies=all_best_energies,
        wall_time=wall_time,
    )
