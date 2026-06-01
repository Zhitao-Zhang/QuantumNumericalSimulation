"""Smart CIM solver — multi-configuration ensemble.

Runs AHC and CAC with several hyperparameter settings in parallel,
then returns the best solution found across all runs.  No J-scaling
or ramp schedule (empirically harmful); the advantage comes purely
from broader exploration of the parameter × random-seed landscape.
"""

import time
import numpy as np
from typing import Optional

from .ahc import solve_ahc
from .cac import solve_cac
from .ising import IsingResult


# Parameter grids chosen to cover the sweet spots for different problem families
_AHC_EPS_GRID = [0.04, 0.07, 0.12, 0.20]
_CAC_GAMMA_GRID = [5e-3, 1e-2, 5e-2]


def solve(
    J: np.ndarray,
    h: Optional[np.ndarray] = None,
    num_timesteps: int = 2000,
    batch_size: int = 10,
    algorithm: str = 'both',
    device=None,
) -> IsingResult:
    """Multi-configuration CIM solver.

    Runs AHC with 4 eps values and CAC with 3 gamma values, each using
    the given batch_size.  Returns the single best solution found.

    Total effective trials = batch_size * (4 + 3) = 7 * batch_size
    when algorithm='both'.
    """
    t0 = time.time()
    results = []

    if algorithm in ('ahc', 'both'):
        for eps in _AHC_EPS_GRID:
            r = solve_ahc(
                J, h=h, num_timesteps=num_timesteps, batch_size=batch_size,
                eps=eps, noise=0, save_trajectories=False, device=device,
            )
            results.append(r)

    if algorithm in ('cac', 'both'):
        for gamma in _CAC_GAMMA_GRID:
            r = solve_cac(
                J, h=h, num_timesteps=num_timesteps, batch_size=batch_size,
                gamma=gamma, noise=0, save_trajectories=False, device=device,
            )
            results.append(r)

    all_spins = np.concatenate([r.all_best_spins for r in results], axis=0)
    all_energies = np.concatenate([r.all_best_energies for r in results])
    total_time = time.time() - t0
    best_idx = int(np.argmin(all_energies))

    return IsingResult(
        best_spins=all_spins[best_idx],
        best_energy=float(all_energies[best_idx]),
        spin_trajectories=None,
        energy_evolution=None,
        all_best_spins=all_spins,
        all_best_energies=all_energies,
        wall_time=total_time,
    )
