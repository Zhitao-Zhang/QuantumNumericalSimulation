"""CIM-Based Seismic Velocity Inversion — Frustrated Binary Tomography.

Key insight: Pure travel-time tomography is "too easy" for greedy methods.
Real seismic inversion uses REGULARIZATION (spatial smoothness) to handle
ill-posedness. This creates FRUSTRATION: data wants heterogeneity, regularization
wants smoothness. The resulting Ising model has BOTH positive and negative
couplings — exactly where CIM excels over classical heuristics.

Formulation:
  E(s) = α * ||T_obs - T_syn(s)||² + β * Σ_{<ij>} (1 - s_i * s_j)

  α: data fidelity, β: smoothness penalty (Potts model)
  This maps to: H = -(1/2) s^T J s - h^T s
    J_ij = -2α (A^T A)_{ij} - β * adjacent(i,j)   [frustration: opposite signs!]
    h_i = 2α (A^T r)_i
"""

import sys
import time
import numpy as np
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import torch
torch.set_num_threads(4)

from seismic.inversion import (
    compute_ray_lengths_multi_source, decode_spins, compute_travel_time,
)
from cimsim.solver import solve as cim_solve
from cimsim.ising import ising_energy


def make_grid_regions(nz, nx, n_vert, n_horiz):
    region_map = np.zeros((nz, nx), dtype=int)
    for iv in range(n_vert):
        z_start = iv * (nz // n_vert)
        z_end = (iv + 1) * (nz // n_vert) if iv < n_vert - 1 else nz
        for ih in range(n_horiz):
            x_start = ih * (nx // n_horiz)
            x_end = (ih + 1) * (nx // n_horiz) if ih < n_horiz - 1 else nx
            region_map[z_start:z_end, x_start:x_end] = iv * n_horiz + ih
    return region_map


def get_adjacency(n_vert, n_horiz):
    """Get list of adjacent region pairs (4-connected grid)."""
    pairs = []
    for iv in range(n_vert):
        for ih in range(n_horiz):
            idx = iv * n_horiz + ih
            if ih + 1 < n_horiz:
                pairs.append((idx, idx + 1))
            if iv + 1 < n_vert:
                pairs.append((idx, idx + n_horiz))
    return pairs


def setup_geometry(nz, nx, dx, dz, n_src=4, n_rec=20):
    sources = []
    receivers_list = []
    recs_bottom = [(i * nx * dx / (n_rec + 1), (nz - 2) * dz)
                   for i in range(1, n_rec + 1)]

    sources.append((nx * dx / 2, 2 * dz))
    receivers_list.append(recs_bottom)
    sources.append((3 * dx, 2 * dz))
    receivers_list.append(recs_bottom)
    if n_src >= 3:
        sources.append(((nx - 3) * dx, 2 * dz))
        receivers_list.append(recs_bottom)
    if n_src >= 4:
        sources.append((2 * dx, nz * dz / 2))
        recs_right = [((nx - 2) * dx, i * nz * dz / (n_rec + 1))
                      for i in range(1, n_rec + 1)]
        receivers_list.append(recs_right)
    return sources, receivers_list


def build_frustrated_ising(L, T_obs, v_slow, v_fast, adj_pairs, alpha, beta):
    """Build frustrated Ising Hamiltonian.

    H = -(1/2) s^T J s - h^T s

    Data term (antiferromagnetic tendencies from A^T A):
      J_data = -2α A^T A (off-diagonal)
      h_data = 2α A^T r

    Smoothness term (ferromagnetic: prefers aligned neighbors):
      J_smooth[i,j] += β for adjacent pairs
      (This is: -β Σ_{<ij>} s_i*s_j which penalizes anti-alignment)

    Combined: J = J_data + J_smooth → FRUSTRATION when signs conflict!
    """
    n_regions = L.shape[1]
    s_slow = 1.0 / v_slow
    s_fast = 1.0 / v_fast
    mean_s = (s_slow + s_fast) / 2.0
    delta_s = (s_fast - s_slow) / 2.0

    T_mean = L @ mean_s
    r = T_obs - T_mean
    A = L * delta_s[np.newaxis, :]

    AtA = A.T @ A
    Atr = A.T @ r

    # Data term
    J = -2.0 * alpha * AtA
    np.fill_diagonal(J, 0.0)
    h = 2.0 * alpha * Atr

    # Smoothness term (ferromagnetic coupling for adjacent regions)
    for i, j in adj_pairs:
        J[i, j] += beta
        J[j, i] += beta

    return J, h


def compute_objective(spins, L, T_obs, v_slow, v_fast, adj_pairs, alpha, beta):
    """Compute full objective: data + smoothness."""
    vel = decode_spins(spins, v_slow, v_fast)
    T_syn = compute_travel_time(L, vel)
    data_misfit = alpha * np.sum((T_obs - T_syn) ** 2)

    smooth_cost = 0.0
    for i, j in adj_pairs:
        if spins[i] != spins[j]:
            smooth_cost += beta

    return data_misfit + smooth_cost


def cim_inversion(J, h, v_slow, v_fast, num_timesteps=3000, batch_size=50):
    """CIM solver."""
    t0 = time.time()
    result = cim_solve(J, h=h, num_timesteps=num_timesteps, batch_size=batch_size)
    best_vel = decode_spins(result.best_spins, v_slow, v_fast)
    total_time = time.time() - t0
    return result.best_spins, best_vel, result.best_energy, total_time


def greedy_search(J, h, v_slow, v_fast, n_restarts=30, seed=42):
    """Multi-start greedy coordinate descent on Ising energy."""
    n = len(h)
    rng = np.random.RandomState(seed)
    t0 = time.time()
    best_energy = float('inf')
    best_spins = None

    for restart in range(n_restarts):
        if restart == 0:
            spins = -np.ones(n)
        elif restart == 1:
            spins = np.ones(n)
        else:
            spins = rng.choice([-1.0, 1.0], size=n)

        for _ in range(n * 5):
            improved = False
            for i in range(n):
                # Compute energy change from flipping spin i
                delta_E = 2.0 * spins[i] * (np.dot(J[i], spins) + h[i])
                if delta_E < 0:  # flip decreases energy
                    spins[i] *= -1
                    improved = True
            if not improved:
                break

        energy = ising_energy(J, spins, h)
        if energy < best_energy:
            best_energy = energy
            best_spins = spins.copy()

    best_vel = decode_spins(best_spins, v_slow, v_fast)
    return best_spins, best_vel, best_energy, time.time() - t0


def simulated_annealing(J, h, v_slow, v_fast, n_steps=100000, seed=123):
    """Simulated annealing on Ising energy."""
    n = len(h)
    rng = np.random.RandomState(seed)
    t0 = time.time()

    spins = rng.choice([-1.0, 1.0], size=n)
    cur_energy = ising_energy(J, spins, h)
    best_energy = cur_energy
    best_spins = spins.copy()

    for step in range(n_steps):
        T = 2.0 * (1.0 - step / n_steps) ** 2
        i = rng.randint(n)
        delta_E = 2.0 * spins[i] * (np.dot(J[i], spins) + h[i])

        if delta_E < 0 or (T > 0 and rng.rand() < np.exp(-delta_E / max(T, 1e-10))):
            spins[i] *= -1
            cur_energy += delta_E
            if cur_energy < best_energy:
                best_energy = cur_energy
                best_spins = spins.copy()

    best_vel = decode_spins(best_spins, v_slow, v_fast)
    return best_spins, best_vel, best_energy, time.time() - t0


def random_search(J, h, v_slow, v_fast, n_trials=50000, seed=99):
    """Random search."""
    n = len(h)
    rng = np.random.RandomState(seed)
    t0 = time.time()
    best_energy = float('inf')
    best_spins = None

    for _ in range(n_trials):
        spins = rng.choice([-1.0, 1.0], size=n)
        energy = ising_energy(J, spins, h)
        if energy < best_energy:
            best_energy = energy
            best_spins = spins.copy()

    best_vel = decode_spins(best_spins, v_slow, v_fast)
    return best_spins, best_vel, best_energy, time.time() - t0


def run_test(name, n_vert, n_horiz, nz, nx, dx, dz,
             v_slow_val, v_fast_val, true_choices,
             alpha, beta, n_src=4, n_rec=16, noise_std=0.0,
             cim_ts=3000, cim_bs=50):
    """Run frustrated inversion test."""
    n_regions = n_vert * n_horiz
    region_map = make_grid_regions(nz, nx, n_vert, n_horiz)
    adj_pairs = get_adjacency(n_vert, n_horiz)

    if isinstance(v_slow_val, (int, float)):
        v_slow = np.full(n_regions, float(v_slow_val))
        v_fast = np.full(n_regions, float(v_fast_val))
    else:
        v_slow = np.array(v_slow_val, dtype=np.float64)
        v_fast = np.array(v_fast_val, dtype=np.float64)

    true_spins = np.array(true_choices, dtype=np.float64) * 2.0 - 1.0
    true_vel = decode_spins(true_spins, v_slow, v_fast)

    print(f"\n{'=' * 70}")
    print(f"  {name}")
    print(f"  N = {n_regions} spins, α = {alpha}, β = {beta}")
    print(f"  Adjacent pairs: {len(adj_pairs)}")
    print(f"{'=' * 70}")

    # Geometry and observations
    sources, receivers_list = setup_geometry(nz, nx, dx, dz, n_src, n_rec)
    L = compute_ray_lengths_multi_source(nz, nx, dz, dx, region_map,
                                         sources, receivers_list, n_regions)
    T_obs = compute_travel_time(L, true_vel)
    if noise_std > 0:
        T_obs += np.random.RandomState(7).randn(len(T_obs)) * noise_std

    # Build frustrated Ising
    J, h = build_frustrated_ising(L, T_obs, v_slow, v_fast, adj_pairs, alpha, beta)

    # Check frustration level
    n_positive_J = np.sum(J > 0) // 2  # symmetric, count once
    n_negative_J = np.sum(J < 0) // 2
    print(f"  J couplings: {n_positive_J} ferro (+), {n_negative_J} antiferro (-)")
    print(f"  J range: [{J.min():.4e}, {J.max():.4e}]")
    print(f"  h range: [{h.min():.4e}, {h.max():.4e}]")

    E_true = ising_energy(J, true_spins, h)
    obj_true = compute_objective(true_spins, L, T_obs, v_slow, v_fast,
                                 adj_pairs, alpha, beta)
    print(f"  True Ising energy: {E_true:.6e}")

    results = {}

    # CIM
    print(f"\n  [CIM] num_timesteps={cim_ts}, batch_size={cim_bs}...")
    cim_spins, cim_vel, cim_E, cim_time = cim_inversion(
        J, h, v_slow, v_fast, cim_ts, cim_bs)
    cim_obj = compute_objective(cim_spins, L, T_obs, v_slow, v_fast,
                                adj_pairs, alpha, beta)
    cim_correct = int(np.sum(cim_spins == true_spins))
    results['CIM'] = (cim_E, cim_obj, cim_time, cim_correct)

    # Greedy
    print(f"  [Greedy] 30 restarts...")
    gr_spins, gr_vel, gr_E, gr_time = greedy_search(J, h, v_slow, v_fast, 30)
    gr_obj = compute_objective(gr_spins, L, T_obs, v_slow, v_fast,
                               adj_pairs, alpha, beta)
    gr_correct = int(np.sum(gr_spins == true_spins))
    results['Greedy'] = (gr_E, gr_obj, gr_time, gr_correct)

    # SA
    n_sa = n_regions * 5000
    print(f"  [SA] {n_sa} steps...")
    sa_spins, sa_vel, sa_E, sa_time = simulated_annealing(
        J, h, v_slow, v_fast, n_sa)
    sa_obj = compute_objective(sa_spins, L, T_obs, v_slow, v_fast,
                               adj_pairs, alpha, beta)
    sa_correct = int(np.sum(sa_spins == true_spins))
    results['SA'] = (sa_E, sa_obj, sa_time, sa_correct)

    # Random
    print(f"  [Random] 50000 trials...")
    rn_spins, rn_vel, rn_E, rn_time = random_search(J, h, v_slow, v_fast, 50000)
    rn_obj = compute_objective(rn_spins, L, T_obs, v_slow, v_fast,
                               adj_pairs, alpha, beta)
    rn_correct = int(np.sum(rn_spins == true_spins))
    results['Random'] = (rn_E, rn_obj, rn_time, rn_correct)

    # Summary
    print(f"\n  {'Method':<12} {'Ising E':>12} {'Objective':>12} {'Correct':>10} {'Time':>8} {'Gap%':>8}")
    print(f"  {'-'*12} {'-'*12} {'-'*12} {'-'*10} {'-'*8} {'-'*8}")
    print(f"  {'True':<12} {E_true:>12.4e} {obj_true:>12.4e} "
          f"{n_regions:>10}/{n_regions} {'-':>8} {'-':>8}")

    for method, (E, obj, t, correct) in results.items():
        gap = (E - E_true) / abs(E_true) * 100 if E_true != 0 else 0
        print(f"  {method:<12} {E:>12.4e} {obj:>12.4e} "
              f"{correct:>10}/{n_regions} {t:>7.2f}s {gap:>7.1f}%")

    return results, E_true


def main():
    print("=" * 70)
    print("  CIM-Based Seismic Velocity Inversion")
    print("  Frustrated Binary Tomography (Data + Smoothness Regularization)")
    print("=" * 70)
    print("\n  The smoothness penalty (Potts model) creates FRUSTRATION:")
    print("  - Data term → anti-ferromagnetic (wants contrast)")
    print("  - Smoothness → ferromagnetic (wants homogeneity)")
    print("  - Competition = NP-hard optimization, ideal for CIM!")

    rng = np.random.RandomState(42)

    # Test 1: N=20
    n1 = 20
    true1 = rng.randint(0, 2, size=n1).tolist()
    res1, E1 = run_test(
        "TEST 1: 4×5 grid (N=20)", n_vert=4, n_horiz=5,
        nz=80, nx=100, dx=10.0, dz=10.0,
        v_slow_val=2500.0, v_fast_val=4500.0,
        true_choices=true1,
        alpha=100.0, beta=0.05,
        n_src=3, n_rec=12,
        cim_ts=2000, cim_bs=40,
    )

    # Test 2: N=36
    n2 = 36
    true2 = rng.randint(0, 2, size=n2).tolist()
    res2, E2 = run_test(
        "TEST 2: 6×6 grid (N=36)", n_vert=6, n_horiz=6,
        nz=120, nx=120, dx=10.0, dz=10.0,
        v_slow_val=2000.0, v_fast_val=5000.0,
        true_choices=true2,
        alpha=50.0, beta=0.08,
        n_src=3, n_rec=15,
        cim_ts=3000, cim_bs=50,
    )

    # Test 3: N=64
    n3 = 64
    true3 = rng.randint(0, 2, size=n3).tolist()
    res3, E3 = run_test(
        "TEST 3: 8×8 grid (N=64)", n_vert=8, n_horiz=8,
        nz=160, nx=160, dx=10.0, dz=10.0,
        v_slow_val=2500.0, v_fast_val=4500.0,
        true_choices=true3,
        alpha=30.0, beta=0.1,
        n_src=4, n_rec=20,
        cim_ts=4000, cim_bs=80,
    )

    # Test 4: N=100
    n4 = 100
    true4 = rng.randint(0, 2, size=n4).tolist()
    res4, E4 = run_test(
        "TEST 4: 10×10 grid (N=100)", n_vert=10, n_horiz=10,
        nz=200, nx=200, dx=10.0, dz=10.0,
        v_slow_val=2000.0, v_fast_val=5000.0,
        true_choices=true4,
        alpha=20.0, beta=0.15,
        n_src=4, n_rec=25,
        cim_ts=5000, cim_bs=100,
    )

    # Final summary
    print("\n" + "=" * 70)
    print("  FINAL SUMMARY — Ising Energy Gap vs True Ground State")
    print("=" * 70)
    all_res = [("N=20", 20, res1, E1), ("N=36", 36, res2, E2),
               ("N=64", 64, res3, E3), ("N=100", 100, res4, E4)]

    print(f"\n  {'Problem':<10} | {'CIM gap%':>10} {'Greedy gap%':>12} {'SA gap%':>10} {'Random gap%':>12}")
    print(f"  {'-'*10}-+-{'-'*10}-{'-'*12}-{'-'*10}-{'-'*12}")
    for name, N, res, E_true in all_res:
        def gap(E):
            return (E - E_true) / abs(E_true) * 100 if E_true != 0 else 0
        c = gap(res['CIM'][0])
        g = gap(res['Greedy'][0])
        s = gap(res['SA'][0])
        r = gap(res['Random'][0])
        print(f"  {name:<10} | {c:>9.2f}% {g:>11.2f}% {s:>9.2f}% {r:>11.2f}%")

    print(f"\n  {'Problem':<10} | {'CIM time':>10} {'Greedy time':>12} {'SA time':>10} {'Random time':>12}")
    print(f"  {'-'*10}-+-{'-'*10}-{'-'*12}-{'-'*10}-{'-'*12}")
    for name, N, res, _ in all_res:
        print(f"  {name:<10} | {res['CIM'][2]:>9.1f}s {res['Greedy'][2]:>11.2f}s "
              f"{res['SA'][2]:>9.1f}s {res['Random'][2]:>11.2f}s")

    print()


if __name__ == "__main__":
    main()
