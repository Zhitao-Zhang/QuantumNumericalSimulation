"""Evaluation framework for comparing CIM solvers against the baseline."""

import sys
import os
import time
import numpy as np
import torch
from pathlib import Path
from dataclasses import dataclass, field

# Small-matrix matmul is bottlenecked by thread-sync overhead at high thread counts.
# 4 threads is the empirical sweet spot for batch_size ≤ 10 and N ≤ 1000.
torch.set_num_threads(4)

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from benchmark.instances import load_instance, INSTANCES_DIR
from baseline.cim_optimizer.CAC import CIM_CAC_GPU
from baseline.cim_optimizer.AHC import CIM_AHC_GPU
from cimsim.cac import solve_cac
from cimsim.ahc import solve_ahc
from cimsim.solver import solve


@dataclass
class SolverResult:
    instance_name: str
    solver_name: str
    best_energy: float
    best_spins: np.ndarray
    wall_time: float
    num_runs: int
    energies_per_run: list = field(default_factory=list)


def ising_energy(J, spins, h=None):
    E = -0.5 * spins @ J @ spins
    if h is not None:
        E -= h @ spins
    return float(E)


# ======================== Baseline Solvers ========================

def run_baseline_cac(J, num_timesteps=2000, batch_size=10, gamma=0.01, **kw):
    t0 = time.time()
    sig_opt, _, _, _, _ = CIM_CAC_GPU(
        T_time=num_timesteps, J=J, batch_size=batch_size, gamma=gamma, noise=0)
    wt = time.time() - t0
    energies = [ising_energy(J, sig_opt[i]) for i in range(batch_size)]
    bi = int(np.argmin(energies))
    return SolverResult("", "baseline_CAC", energies[bi], sig_opt[bi], wt, batch_size, energies)


def run_baseline_ahc(J, num_timesteps=2000, batch_size=10, eps=0.07, **kw):
    t0 = time.time()
    sig_opt, _, _, _, _ = CIM_AHC_GPU(
        T_time=num_timesteps, J=J, batch_size=batch_size, eps=eps, noise=0)
    wt = time.time() - t0
    energies = [ising_energy(J, sig_opt[i]) for i in range(batch_size)]
    bi = int(np.argmin(energies))
    return SolverResult("", "baseline_AHC", energies[bi], sig_opt[bi], wt, batch_size, energies)


# ======================== Our Solvers ========================

def run_our_cac(J, num_timesteps=2000, batch_size=10, gamma=0.01, **kw):
    r = solve_cac(J, num_timesteps=num_timesteps, batch_size=batch_size,
                  gamma=gamma, noise=0, save_trajectories=False)
    return SolverResult("", "ours_CAC", r.best_energy, r.best_spins,
                        r.wall_time, batch_size, r.all_best_energies.tolist())


def run_our_ahc(J, num_timesteps=2000, batch_size=10, eps=0.07, **kw):
    r = solve_ahc(J, num_timesteps=num_timesteps, batch_size=batch_size,
                  eps=eps, noise=0, save_trajectories=False)
    return SolverResult("", "ours_AHC", r.best_energy, r.best_spins,
                        r.wall_time, batch_size, r.all_best_energies.tolist())


def run_our_smart(J, num_timesteps=2000, batch_size=10, **kw):
    r = solve(J, num_timesteps=num_timesteps, batch_size=batch_size, algorithm='both')
    return SolverResult("", "ours_SMART", r.best_energy, r.best_spins,
                        r.wall_time, len(r.all_best_energies),
                        r.all_best_energies.tolist())


# ======================== Evaluation ========================

def print_comparison(all_results: dict):
    print("\n" + "=" * 110)
    print(f"{'Instance':<25} {'Solver':<20} {'Best Energy':>12} {'Mean Energy':>12} "
          f"{'Std Energy':>12} {'Time (s)':>10} {'Gap%':>8}")
    print("=" * 110)

    for instance_name, solvers in all_results.items():
        global_best = min(r.best_energy for r in solvers.values())
        for solver_name, result in solvers.items():
            energies = result.energies_per_run
            mean_e = np.mean(energies)
            std_e = np.std(energies)
            gap = 100.0 * (result.best_energy - global_best) / abs(global_best) if global_best != 0 else 0
            print(f"{instance_name:<25} {solver_name:<20} {result.best_energy:>12.2f} "
                  f"{mean_e:>12.2f} {std_e:>12.2f} {result.wall_time:>10.4f} {gap:>7.2f}%")
        print("-" * 110)

    # Summary: count wins for ours_SMART vs best baseline
    print("\n=== Summary: ours_SMART vs best_baseline ===")
    wins = ties = losses = 0
    for instance_name, solvers in all_results.items():
        bl_best = min(solvers[k].best_energy for k in solvers if k.startswith("baseline"))
        our_best = solvers.get("ours_SMART", solvers.get("ours_AHC"))
        if our_best is None:
            continue
        if our_best.best_energy < bl_best - 0.01:
            wins += 1
            tag = "WIN"
        elif our_best.best_energy > bl_best + 0.01:
            losses += 1
            tag = "LOSS"
        else:
            ties += 1
            tag = "TIE"
        print(f"  {instance_name:<25} baseline={bl_best:>10.2f}  ours={our_best.best_energy:>10.2f}  {tag}")
    print(f"\n  Wins={wins}  Ties={ties}  Losses={losses}")
    print()


def run_full_benchmark():
    instance_files = sorted(INSTANCES_DIR.glob("*.npz"))
    if not instance_files:
        from benchmark.instances import generate_benchmark_suite
        generate_benchmark_suite()
        instance_files = sorted(INSTANCES_DIR.glob("*.npz"))

    solvers = {
        "baseline_CAC": lambda J, **kw: run_baseline_cac(J, gamma=0.01, **kw),
        "baseline_AHC": lambda J, **kw: run_baseline_ahc(J, eps=0.07, **kw),
        "ours_CAC":     lambda J, **kw: run_our_cac(J, gamma=0.01, **kw),
        "ours_AHC":     lambda J, **kw: run_our_ahc(J, eps=0.07, **kw),
        "ours_SMART":   lambda J, **kw: run_our_smart(J, **kw),
    }

    all_results = {}
    for instance_file in instance_files:
        name = instance_file.stem
        print(f"\nEvaluating instance: {name} ...")
        J, h = load_instance(name)
        N = J.shape[0]
        ts = min(2000, max(500, N * 10))
        bs = min(10, max(1, 1000 // N))

        results = {}
        for sname, sfn in solvers.items():
            result = sfn(J, num_timesteps=ts, batch_size=bs)
            result.instance_name = name
            results[sname] = result
        all_results[name] = results

    print_comparison(all_results)
    return all_results


if __name__ == "__main__":
    run_full_benchmark()
