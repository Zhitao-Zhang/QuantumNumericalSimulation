"""Scaling experiment: CIM vs Classical at Large N.

Fair comparison strategy:
  - CIM: FIXED parameters (timesteps=2000, batch=30) for all N
  - Greedy/SA: given the SAME wall-clock time as CIM, how good are they?
  - Also compare at their own natural time budgets

Key scaling:
  - CIM cost per run: O(N² × timesteps × batch) — dominated by matmul
  - Greedy (one restart): O(N² × sweeps) — dominated by local_field update
  - SA step: O(N) — one random flip + field update
"""

import sys
import time
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import torch
torch.set_num_threads(4)

from cimsim.solver import solve as cim_solve
from cimsim.ising import ising_energy


def make_frustrated_ising(N, density=0.3, field_scale=0.5, seed=42):
    """Generate a frustrated Ising instance mimicking seismic tomography.

    Mix of:
    - Dense antiferromagnetic couplings (from data term A^T A)
    - Sparse ferromagnetic couplings (from smoothness regularization)
    - External field (from data residuals)
    """
    rng = np.random.RandomState(seed)
    n_side = int(np.sqrt(N))

    # Dense part: random matrix simulating A^T A structure
    # A is (M × N) with M ~ N/2 observations
    M = max(N // 2, 10)
    A = rng.randn(M, N) * 0.1
    AtA = A.T @ A

    # Sparse ferromagnetic part: grid adjacency
    J_smooth = np.zeros((N, N))
    beta = 0.3 * np.mean(np.abs(AtA[np.triu_indices(N, 1)]))
    for i in range(N):
        row, col = i // n_side, i % n_side
        if col + 1 < n_side and i + 1 < N:
            j = i + 1
            J_smooth[i, j] = beta
            J_smooth[j, i] = beta
        if row + 1 < n_side:
            j = i + n_side
            if j < N:
                J_smooth[i, j] = beta
                J_smooth[j, i] = beta

    # Combined J
    J = -AtA + J_smooth
    np.fill_diagonal(J, 0.0)
    J = (J + J.T) / 2.0

    # External field
    h = rng.randn(N) * field_scale * np.sqrt(np.mean(np.diag(AtA)))

    return J, h


def solve_greedy_timed(J, h, time_budget):
    """Greedy with time budget."""
    n = len(h)
    rng = np.random.RandomState(42)
    t0 = time.time()
    best_E = float('inf')
    restarts = 0

    while time.time() - t0 < time_budget:
        if restarts == 0:
            spins = -np.ones(n)
        elif restarts == 1:
            spins = np.ones(n)
        else:
            spins = rng.choice([-1.0, 1.0], size=n)

        local_field = J @ spins + h
        for _ in range(n * 3):
            improved = False
            for i in range(n):
                delta_E = 2.0 * spins[i] * local_field[i]
                if delta_E < 0:
                    local_field += 2.0 * (-spins[i]) * J[:, i]
                    spins[i] *= -1
                    improved = True
            if not improved:
                break

        E = ising_energy(J, spins, h)
        if E < best_E:
            best_E = E
        restarts += 1

    return best_E, time.time() - t0, restarts


def solve_sa_timed(J, h, time_budget):
    """SA with time budget."""
    n = len(h)
    rng = np.random.RandomState(123)
    t0 = time.time()

    spins = rng.choice([-1.0, 1.0], size=n)
    local_field = J @ spins + h
    cur_E = float(ising_energy(J, spins, h))
    best_E = cur_E
    steps = 0
    total_steps_target = n * 10000  # will be cut by time

    while time.time() - t0 < time_budget and steps < total_steps_target:
        T = 2.0 * max(0, 1.0 - steps / total_steps_target) ** 2
        i = rng.randint(n)
        delta_E = 2.0 * spins[i] * local_field[i]

        if delta_E < 0 or (T > 1e-10 and rng.rand() < np.exp(-delta_E / T)):
            local_field += 2.0 * (-spins[i]) * J[:, i]
            spins[i] *= -1
            cur_E += delta_E
            if cur_E < best_E:
                best_E = cur_E
        steps += 1

    return best_E, time.time() - t0, steps


def solve_greedy_fixed(J, h, n_restarts=30):
    """Greedy with fixed restarts."""
    n = len(h)
    rng = np.random.RandomState(42)
    t0 = time.time()
    best_E = float('inf')

    for restart in range(n_restarts):
        if restart == 0:
            spins = -np.ones(n)
        elif restart == 1:
            spins = np.ones(n)
        else:
            spins = rng.choice([-1.0, 1.0], size=n)

        local_field = J @ spins + h
        for _ in range(n * 3):
            improved = False
            for i in range(n):
                delta_E = 2.0 * spins[i] * local_field[i]
                if delta_E < 0:
                    local_field += 2.0 * (-spins[i]) * J[:, i]
                    spins[i] *= -1
                    improved = True
            if not improved:
                break

        E = ising_energy(J, spins, h)
        if E < best_E:
            best_E = E

    return best_E, time.time() - t0


def solve_sa_fixed(J, h, n_steps):
    """SA with fixed steps."""
    n = len(h)
    rng = np.random.RandomState(123)
    t0 = time.time()

    spins = rng.choice([-1.0, 1.0], size=n)
    local_field = J @ spins + h
    cur_E = float(ising_energy(J, spins, h))
    best_E = cur_E

    for step in range(n_steps):
        T = 2.0 * (1.0 - step / n_steps) ** 2
        i = rng.randint(n)
        delta_E = 2.0 * spins[i] * local_field[i]
        if delta_E < 0 or (T > 1e-10 and rng.rand() < np.exp(-delta_E / T)):
            local_field += 2.0 * (-spins[i]) * J[:, i]
            spins[i] *= -1
            cur_E += delta_E
            if cur_E < best_E:
                best_E = cur_E

    return best_E, time.time() - t0


def main():
    print("=" * 70)
    print("  SCALING: CIM vs Greedy vs SA")
    print("  Fixed CIM params: timesteps=2000, batch=30")
    print("  Comparison: (1) fixed budget  (2) equal time budget")
    print("=" * 70)

    sizes = [50, 100, 200, 400, 800]
    results = []

    for N in sizes:
        print(f"\n{'─' * 60}")
        print(f"  N = {N}")
        print(f"{'─' * 60}")

        J, h = make_frustrated_ising(N, seed=42 + N)

        # CIM with fixed params
        t0 = time.time()
        cim_result = cim_solve(J, h=h, num_timesteps=2000, batch_size=30)
        cim_E = cim_result.best_energy
        cim_t = time.time() - t0
        print(f"  CIM:         E = {cim_E:>12.4f}   t = {cim_t:>7.2f}s")

        # Greedy with 30 restarts
        greedy_E, greedy_t = solve_greedy_fixed(J, h, 30)
        print(f"  Greedy(30):  E = {greedy_E:>12.4f}   t = {greedy_t:>7.2f}s")

        # SA with N*3000 steps
        sa_steps = N * 3000
        sa_E, sa_t = solve_sa_fixed(J, h, sa_steps)
        print(f"  SA({sa_steps:>7}): E = {sa_E:>12.4f}   t = {sa_t:>7.2f}s")

        # Equal-time comparison: give Greedy/SA same time as CIM
        print(f"\n  [Equal time = {cim_t:.1f}s]")
        greedy_eq_E, greedy_eq_t, greedy_eq_r = solve_greedy_timed(J, h, cim_t)
        sa_eq_E, sa_eq_t, sa_eq_steps = solve_sa_timed(J, h, cim_t)
        print(f"  Greedy(={cim_t:.0f}s): E = {greedy_eq_E:>12.4f}  "
              f"({greedy_eq_r} restarts)")
        print(f"  SA(={cim_t:.0f}s):     E = {sa_eq_E:>12.4f}  "
              f"({sa_eq_steps} steps)")

        best_E = min(cim_E, greedy_E, sa_E, greedy_eq_E, sa_eq_E)
        results.append({
            'N': N,
            'cim': (cim_E, cim_t),
            'greedy': (greedy_E, greedy_t),
            'sa': (sa_E, sa_t),
            'greedy_eq': (greedy_eq_E, greedy_eq_t, greedy_eq_r),
            'sa_eq': (sa_eq_E, sa_eq_t, sa_eq_steps),
            'best': best_E,
        })

    # Summary
    print("\n" + "=" * 70)
    print("  SUMMARY: Fixed Budget Comparison")
    print("=" * 70)
    print(f"\n  {'N':>5} | {'CIM E':>10} {'CIM t':>8} | "
          f"{'Greedy E':>10} {'Greedy t':>9} | "
          f"{'SA E':>10} {'SA t':>8}")
    print(f"  {'-'*5}-+-{'-'*10}-{'-'*8}-+-{'-'*10}-{'-'*9}-+-{'-'*10}-{'-'*8}")
    for res in results:
        print(f"  {res['N']:>5} | {res['cim'][0]:>10.2f} {res['cim'][1]:>7.1f}s | "
              f"{res['greedy'][0]:>10.2f} {res['greedy'][1]:>8.2f}s | "
              f"{res['sa'][0]:>10.2f} {res['sa'][1]:>7.1f}s")

    print(f"\n  EQUAL-TIME comparison (all methods get CIM's time budget):")
    print(f"  {'N':>5} | {'CIM E':>10} | {'Greedy(=t)':>12} {'restarts':>9} | "
          f"{'SA(=t)':>10} {'steps':>10}")
    print(f"  {'-'*5}-+-{'-'*10}-+-{'-'*12}-{'-'*9}-+-{'-'*10}-{'-'*10}")
    for res in results:
        print(f"  {res['N']:>5} | {res['cim'][0]:>10.2f} | "
              f"{res['greedy_eq'][0]:>12.2f} {res['greedy_eq'][2]:>9} | "
              f"{res['sa_eq'][0]:>10.2f} {res['sa_eq'][2]:>10}")

    # Quality gap at equal time
    print(f"\n  Energy gap vs best (at equal time, lower = better):")
    print(f"  {'N':>5} | {'CIM':>10} {'Greedy':>10} {'SA':>10} | {'Winner':>8}")
    print(f"  {'-'*5}-+-{'-'*10}-{'-'*10}-{'-'*10}-+-{'-'*8}")
    for res in results:
        best = min(res['cim'][0], res['greedy_eq'][0], res['sa_eq'][0])
        c = res['cim'][0] - best
        g = res['greedy_eq'][0] - best
        s = res['sa_eq'][0] - best
        winner = 'CIM' if c <= g and c <= s else ('Greedy' if g <= s else 'SA')
        print(f"  {res['N']:>5} | {c:>10.4f} {g:>10.4f} {s:>10.4f} | {winner:>8}")

    print()


if __name__ == "__main__":
    main()
