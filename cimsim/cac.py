"""Chaotic Amplitude Control (CAC) CIM solver.

Reference:
    T. Leleu, F. Khoyratee, T. Levi, R. Hamerly, T. Kohno, K. Aihara,
    Commun. Phys. 4, 266 (2021).
"""

import time
import numpy as np
import torch
from typing import Optional, Callable

from .ising import IsingResult

if torch.get_num_threads() > 8:
    torch.set_num_threads(4)


@torch.no_grad()
def solve_cac(
    J: np.ndarray,
    h: Optional[np.ndarray] = None,
    num_timesteps: int = 1000,
    dt: float = 0.05,
    batch_size: int = 1,
    r: Optional[float] = None,
    alpha: float = 3.0,
    beta: float = 0.25,
    gamma: float = 1.1e-4,
    delta: float = 10.0,
    mu: float = 1.0,
    rho: float = 3.0,
    tau: float = 1000.0,
    noise: float = 0.0,
    pump_schedule: Optional[Callable] = None,
    feedback_schedule: Optional[Callable] = None,
    nonlinearity: Optional[Callable] = None,
    device: Optional[torch.device] = None,
    save_trajectories: bool = True,
) -> IsingResult:
    if device is None:
        device = torch.device('cpu')

    J_t = torch.from_numpy(np.ascontiguousarray(J)).float().to(device)
    N = J_t.shape[0]
    has_h = h is not None
    h_t = torch.from_numpy(np.ascontiguousarray(h)).float().to(device) if has_h else None

    if r is None:
        r = 0.8 - (N / 220.0) ** 2
    if nonlinearity is None:
        nonlinearity = torch.tanh

    EMAX = 32.0
    effective_tau = tau / dt
    gamma_dt = gamma * dt

    # --- Build schedules as contiguous tensors ---
    step_idx = torch.arange(num_timesteps, device=device)
    r_sched = (pump_schedule(step_idx) if pump_schedule else
               torch.full((num_timesteps,), r, device=device))
    beta_sched = (feedback_schedule(step_idx) if feedback_schedule else
                  torch.full((num_timesteps,), beta, device=device))

    # --- Initialize state ---
    x = torch.empty(batch_size, N, device=device).uniform_(-5e-4, 5e-4)
    e = torch.ones(batch_size, N, device=device)
    xi = torch.zeros(batch_size, device=device)
    a = torch.full((batch_size,), alpha, device=device)
    t_c = torch.zeros(batch_size, device=device)

    # Pre-allocate scratch tensors
    sig = torch.empty_like(x)
    MVM = torch.empty_like(x)
    x_sq = torch.empty_like(x)
    tmp = torch.empty_like(x)
    H = torch.empty(batch_size, device=device)

    # Initial energy & best tracking
    torch.sign(x, out=sig)  # approximate; correct sign below
    sig = (2.0 * (x > 0) - 1.0).float()
    H_opt = -0.5 * (sig * (sig @ J_t)).sum(dim=1)
    if has_h:
        H_opt -= sig @ h_t
    sig_opt = sig.clone()

    if save_trajectories:
        spin_traj = torch.zeros(batch_size, N, num_timesteps, device=device)
        energy_traj = torch.zeros(batch_size, num_timesteps, device=device)

    t0 = time.time()

    for t in range(num_timesteps):
        if save_trajectories:
            spin_traj[:, :, t] = x

        # Binarize & energy
        sig.copy_((2.0 * (x > 0) - 1.0).float())
        sJ = sig @ J_t
        H = -0.5 * (sig * sJ).sum(dim=1)
        if has_h:
            H -= sig @ h_t

        if save_trajectories:
            energy_traj[:, t] = H

        # Pre-pump MVM
        torch.mul(x, x, out=x_sq)
        torch.mm(x, J_t, out=MVM)

        # Pump + saturation: x += dt * x * ((r-1) - mu*x^2)
        r_t = r_sched[t]
        x += dt * (x * ((r_t - 1.0) - mu * x_sq))

        # Feedback: x += dt * beta * MVM * e
        beta_t = beta_sched[t]
        x += (dt * beta_t) * (MVM * e)

        # Noise
        if noise > 0:
            x += beta_t * noise * (torch.rand_like(x) - 0.5)

        # Error variable: e += dt * (-xi * (x^2 - a) * e)
        e += dt * ((-xi[:, None]) * (x.square() - a[:, None]) * e)
        e.clamp_(max=EMAX)

        # Chaotic amplitude control
        xi += gamma_dt
        dH = H - H_opt
        a = alpha + rho * nonlinearity(delta * dH)

        # Reset (baseline ordering: t_c update invalidates xi check)
        t_c[t_c < t - effective_tau] = t

        # Track best
        improved = H < H_opt
        if improved.any():
            sig_opt[improved] = sig[improved]
            t_c[improved] = t
            H_opt = torch.minimum(H_opt, H)

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
