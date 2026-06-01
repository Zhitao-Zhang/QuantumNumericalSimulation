"""Scientific computing benchmarks for the CIM Ising solver.

Covers five domains:
  1. Condensed-matter physics  — 2D Edwards-Anderson spin glass ground states
  2. Graph theory              — Max-Cut on classic named graphs
  3. Combinatorial optimization — Number partitioning problem
  4. Statistical mechanics      — Frustrated antiferromagnet (triangular lattice)
  5. Operations research        — Weighted MAX-2-SAT
"""

import sys
import time
import numpy as np
import torch
from pathlib import Path
from dataclasses import dataclass
from typing import Optional, Tuple, List

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
torch.set_num_threads(4)

from cimsim.solver import solve
from cimsim.ising import ising_energy


# ===================================================================
#  Exact solver  (brute-force enumeration, feasible for N ≤ 24)
# ===================================================================

def exact_ground_state(J: np.ndarray, h: Optional[np.ndarray] = None,
                       max_n: int = 25) -> Tuple[np.ndarray, float]:
    """Find the exact Ising ground state by exhaustive enumeration.

    Processes configurations in batches for speed. Feasible up to N~25.
    """
    N = J.shape[0]
    assert N <= max_n, f"Brute force infeasible for N={N} (limit {max_n})"
    Jf = J.astype(np.float64)
    hf = h.astype(np.float64) if h is not None else None
    total = 1 << N
    BATCH = min(total, 1 << 16)

    best_E = float('inf')
    best_bits = 0

    bit_masks = (1 << np.arange(N)).astype(np.int64)

    for start in range(0, total, BATCH):
        end = min(start + BATCH, total)
        indices = np.arange(start, end, dtype=np.int64)
        spins = 2.0 * ((indices[:, None] & bit_masks[None, :]) > 0).astype(np.float64) - 1.0
        Js = spins @ Jf                     # (batch, N)
        E = -0.5 * np.sum(spins * Js, axis=1)  # (batch,)
        if hf is not None:
            E -= spins @ hf
        idx = np.argmin(E)
        if E[idx] < best_E:
            best_E = E[idx]
            best_bits = int(indices[idx])

    best_s = 2.0 * ((np.int64(best_bits) & bit_masks) > 0).astype(np.float64) - 1.0
    return best_s, float(best_E)


# ===================================================================
#  Problem 1 — 2D Edwards-Anderson Spin Glass
# ===================================================================

def make_2d_spin_glass(L: int, seed: int = 0) -> Tuple[np.ndarray, str]:
    """Square lattice spin glass with random ±1 nearest-neighbor couplings.

    Returns (J, description).  N = L*L.
    """
    rng = np.random.RandomState(seed)
    N = L * L
    J = np.zeros((N, N))
    for y in range(L):
        for x in range(L):
            i = y * L + x
            # right neighbor
            if x + 1 < L:
                j = y * L + (x + 1)
                c = rng.choice([-1.0, 1.0])
                J[i, j] = c; J[j, i] = c
            # down neighbor
            if y + 1 < L:
                j = (y + 1) * L + x
                c = rng.choice([-1.0, 1.0])
                J[i, j] = c; J[j, i] = c
    desc = f"2D Spin Glass {L}×{L} (N={N})"
    return J, desc


# ===================================================================
#  Problem 2 — Max-Cut on Classic Graphs
# ===================================================================

def _graph_to_maxcut_J(edges: List[Tuple[int, int]], N: int,
                       weights=None) -> np.ndarray:
    """Convert an edge list to the Ising J for Max-Cut.

    Max-Cut:  maximize Σ_{(i,j)∈E} (1 - s_i s_j)/2
    Equivalent to minimizing H = -Σ_{(i,j)} w_{ij}/2 * (1 - s_i s_j)
     which maps to Ising with J_ij = -w_ij  (we want to minimize s^T J s
     with a sign flip: cut size = const - 1/2 s^T J s).
    Simplest: set J_ij = +1 for edges. Then H = -1/2 s^T J s,
    and the minimum H corresponds to maximum cut.
    """
    J = np.zeros((N, N))
    for idx, (i, j) in enumerate(edges):
        w = weights[idx] if weights is not None else 1.0
        J[i, j] = w
        J[j, i] = w
    return J


def make_petersen_graph() -> Tuple[np.ndarray, str, int]:
    """Petersen graph: 10 vertices, 15 edges, max-cut = 12."""
    edges = [
        (0,1),(1,2),(2,3),(3,4),(4,0),       # outer pentagon
        (5,7),(7,9),(9,6),(6,8),(8,5),       # inner pentagram
        (0,5),(1,6),(2,7),(3,8),(4,9),       # spokes
    ]
    J = _graph_to_maxcut_J(edges, 10)
    return J, "Max-Cut: Petersen Graph (N=10)", 12


def make_cube_graph() -> Tuple[np.ndarray, str, int]:
    """3D Hypercube (cube) graph: 8 vertices, 12 edges, max-cut = 12."""
    edges = [
        (0,1),(1,3),(3,2),(2,0),  # bottom face
        (4,5),(5,7),(7,6),(6,4),  # top face
        (0,4),(1,5),(2,6),(3,7),  # vertical edges
    ]
    J = _graph_to_maxcut_J(edges, 8)
    return J, "Max-Cut: Cube Graph (N=8)", 12


def make_dodecahedron_graph() -> Tuple[np.ndarray, str, None]:
    """Dodecahedral graph: 20 vertices, 30 edges. Max-cut found by exact solver."""
    edges = [
        (0,1),(1,2),(2,3),(3,4),(4,0),
        (0,5),(1,6),(2,7),(3,8),(4,9),
        (5,10),(6,11),(7,12),(8,13),(9,14),
        (10,11),(11,12),(12,13),(13,14),(14,10),
        (5,15),(6,16),(7,17),(8,18),(9,19),
        (15,16),(16,17),(17,18),(18,19),(19,15),
    ]
    J = _graph_to_maxcut_J(edges, 20)
    return J, "Max-Cut: Dodecahedron (N=20)", None


# ===================================================================
#  Problem 3 — Number Partitioning
# ===================================================================

def make_number_partition(numbers: np.ndarray) -> Tuple[np.ndarray, np.ndarray, str]:
    """Map the number-partitioning problem to Ising.

    Goal: find s_i ∈ {-1,+1} minimizing (Σ a_i s_i)^2.
    Ising: J_ij = -2 a_i a_j (i≠j), H = (Σ a_i s_i)^2 - Σ a_i^2 + const.
    Returns (J, numbers, description).
    """
    a = numbers.astype(np.float64)
    N = len(a)
    J = -2.0 * np.outer(a, a)
    np.fill_diagonal(J, 0.0)
    return J, a, f"Number Partitioning (N={N}, sum={int(a.sum())})"


# ===================================================================
#  Problem 4 — Frustrated Triangular Antiferromagnet
# ===================================================================

def make_triangular_antiferromagnet(L: int) -> Tuple[np.ndarray, str]:
    """Triangular lattice antiferromagnet — all nearest-neighbor couplings = -1.

    This is a classic frustrated system from statistical mechanics.
    The ground state has extensive degeneracy. N = L*L.
    """
    N = L * L
    J = np.zeros((N, N))
    for y in range(L):
        for x in range(L):
            i = y * L + x
            # right
            if x + 1 < L:
                j = y * L + (x + 1)
                J[i, j] = -1.0; J[j, i] = -1.0
            # down
            if y + 1 < L:
                j = (y + 1) * L + x
                J[i, j] = -1.0; J[j, i] = -1.0
            # diagonal (down-right): triangular lattice
            if x + 1 < L and y + 1 < L:
                j = (y + 1) * L + (x + 1)
                J[i, j] = -1.0; J[j, i] = -1.0
    desc = f"Triangular AF {L}×{L} (N={N})"
    return J, desc


# ===================================================================
#  Problem 5 — Weighted MAX-2-SAT
# ===================================================================

def make_weighted_max2sat(num_vars: int, num_clauses: int,
                          seed: int = 0) -> Tuple[np.ndarray, np.ndarray, list, str]:
    """Random weighted MAX-2-SAT instance mapped to Ising.

    Each clause (l_a OR l_b) with weight w is mapped to:
      w/4 * (1 + J_a s_a + J_b s_b - J_a J_b s_a s_b)
    where J_a = -1 if l_a is negated, else +1 (similarly for J_b).
    Returns (J, h, clauses, description).
    """
    rng = np.random.RandomState(seed)
    N = num_vars
    J = np.zeros((N, N))
    h = np.zeros(N)
    clauses = []
    for _ in range(num_clauses):
        i, j = rng.choice(N, 2, replace=False)
        neg_i = rng.rand() < 0.5
        neg_j = rng.rand() < 0.5
        w = rng.uniform(1.0, 10.0)
        clauses.append((i, neg_i, j, neg_j, w))

        si = -1.0 if neg_i else 1.0
        sj = -1.0 if neg_j else 1.0

        # Ising contribution: minimize -(satisfied clauses)
        # (l_i OR l_j) is unsatisfied only when l_i=False AND l_j=False
        # Penalty for both false: w * 1/4 * (1 + si*s_i + sj*s_j + si*sj*s_i*s_j)
        # → h_i -= w/4 * si,  h_j -= w/4 * sj,  J_ij -= w/2 * si * sj
        h[i] -= (w / 4.0) * si
        h[j] -= (w / 4.0) * sj
        J[i, j] -= (w / 2.0) * si * sj
        J[j, i] -= (w / 2.0) * si * sj

    desc = f"Weighted MAX-2-SAT ({num_vars} vars, {num_clauses} clauses)"
    return J, h, clauses, desc


def satisfied_clauses(spins: np.ndarray, clauses: list) -> Tuple[int, float]:
    """Count satisfied clauses and total satisfied weight."""
    count = 0
    total_w = 0.0
    for (i, ni, j, nj, w) in clauses:
        li = (spins[i] == -1) if ni else (spins[i] == 1)
        lj = (spins[j] == -1) if nj else (spins[j] == 1)
        if li or lj:
            count += 1
            total_w += w
    return count, total_w


# ===================================================================
#  Runner
# ===================================================================

@dataclass
class ScientificResult:
    name: str
    N: int
    exact_energy: float
    cim_energy: float
    gap_pct: float
    exact_time: float
    cim_time: float
    extra: str = ""


def run_problem(name: str, J: np.ndarray, h: Optional[np.ndarray] = None,
                known_optimal: Optional[float] = None,
                extra_fn=None) -> ScientificResult:
    """Solve a problem with exact solver and CIM, return comparison."""
    N = J.shape[0]

    # --- Exact ---
    t0 = time.time()
    if N <= 25:
        exact_spins, exact_E = exact_ground_state(J, h)
        exact_time = time.time() - t0
    elif known_optimal is not None:
        exact_E = known_optimal
        exact_spins = None
        exact_time = 0.0
    else:
        exact_E = None
        exact_spins = None
        exact_time = 0.0

    # --- CIM ---
    ts = min(3000, max(500, N * 20))
    bs = min(20, max(2, 2000 // N))
    cim_result = solve(J, h=h, num_timesteps=ts, batch_size=bs)
    cim_E = cim_result.best_energy
    cim_spins = cim_result.best_spins
    cim_time = cim_result.wall_time

    # --- Gap ---
    if exact_E is not None and exact_E != 0:
        gap = 100.0 * (cim_E - exact_E) / abs(exact_E)
    elif exact_E is not None:
        gap = abs(cim_E - exact_E)
    else:
        gap = float('nan')

    extra = ""
    if extra_fn is not None:
        extra = extra_fn(cim_spins, exact_spins)

    return ScientificResult(name, N, exact_E if exact_E is not None else float('nan'),
                            cim_E, gap, exact_time, cim_time, extra)


def maxcut_size(J: np.ndarray, spins: np.ndarray) -> int:
    """Compute the cut size from Ising J (adjacency matrix) and spins."""
    cut = 0.0
    N = J.shape[0]
    for i in range(N):
        for j in range(i+1, N):
            if J[i, j] != 0 and spins[i] != spins[j]:
                cut += abs(J[i, j])
    return int(cut)


def partition_residue(a: np.ndarray, spins: np.ndarray) -> int:
    return int(abs(np.dot(a, spins)))


def run_all():
    results: List[ScientificResult] = []

    # ---- 1. 2D Spin Glass ----
    print("=" * 70)
    print("  Domain 1: Condensed Matter Physics — 2D Spin Glass")
    print("=" * 70)
    for L in [3, 4, 5]:
        J, desc = make_2d_spin_glass(L, seed=42)
        r = run_problem(desc, J)
        results.append(r)
        print(f"  [{r.name}]  exact={r.exact_energy:.4f}  "
              f"CIM={r.cim_energy:.4f}  gap={r.gap_pct:.2f}%  "
              f"(exact {r.exact_time:.2f}s, CIM {r.cim_time:.2f}s)")

    # ---- 2. Max-Cut on Classic Graphs ----
    print()
    print("=" * 70)
    print("  Domain 2: Graph Theory — Max-Cut on Famous Graphs")
    print("=" * 70)

    for make_fn in [make_petersen_graph, make_cube_graph, make_dodecahedron_graph]:
        ret = make_fn()
        J, desc = ret[0], ret[1]
        known_cut = ret[2] if len(ret) > 2 else None

        def cut_info(cim_s, exact_s, _J=J, _kc=known_cut):
            c = maxcut_size(_J, cim_s)
            s = f"CIM cut={c}"
            if exact_s is not None:
                s += f", exact cut={maxcut_size(_J, exact_s)}"
            if _kc is not None:
                s += f", known optimal={_kc}"
            return s

        # For max-cut, optimal energy = -(total_edges - max_cut)/2... 
        # Actually with J_ij=1 for edges: H = -1/2 s^T J s.
        # If all edges cut: s_i s_j = -1 for all, H = +|E|/2.
        # If no edges cut: H = -|E|/2.
        # Max-cut maximizes the number of disagreeing edges.
        # Our Ising minimizes H. Fewer agreeing edges = lower H... actually:
        # H = -1/2 Σ J_{ij} s_i s_j. With J=+1 on edges,
        # H = -1/2 Σ_{edges} s_i s_j = -1/2 (agree - disagree)
        # = -1/2 ((|E| - cut) - cut) = -1/2 (|E| - 2*cut) = cut - |E|/2.
        # So H = cut - |E|/2 → minimize H = maximize cut? No!
        # H = -|E|/2 + cut → more cut = higher H.  
        # To maximize cut we need to MAXIMIZE H... but we minimize.
        # With J=+1: we'd want J=-1 on edges for max-cut.
        # Let me fix: set J_ij = -1 for edges. Then:
        # H = -1/2 Σ (-1) s_i s_j = 1/2 Σ s_i s_j
        # = 1/2 (agree - disagree) = 1/2 (|E| - 2*cut) = |E|/2 - cut.
        # Minimize H → maximize cut. ✓

        # Actually I need to fix _graph_to_maxcut_J...
        # The function currently sets J[i,j] = +w. But for max-cut minimization:
        # We want J[i,j] = -w so that disagreeing spins (max-cut) gives lower energy.
        # Let me not change the function but negate J here.
        J_mc = -J  # negate: J=-1 on edges → min H = max cut
        known_E = -known_cut + 0.5 * abs(J).sum() / 2 if known_cut is not None else None
        # H_min = |E|/2 - max_cut. With weights=1: |E| = number of edges.
        num_edges = int(np.sum(np.abs(J)) / 2)
        if known_cut is not None:
            known_E = num_edges / 2.0 - known_cut
        else:
            known_E = None

        r = run_problem(desc, J_mc, known_optimal=known_E, extra_fn=cut_info)
        # Recover cut size from energy: cut = |E|/2 - H
        cim_cut = num_edges / 2.0 - r.cim_energy
        exact_cut = num_edges / 2.0 - r.exact_energy if not np.isnan(r.exact_energy) else "?"
        results.append(r)
        print(f"  [{r.name}]  |E|={num_edges}")
        print(f"    exact cut={exact_cut}  CIM cut={cim_cut:.0f}  gap={r.gap_pct:.2f}%")
        print(f"    (exact {r.exact_time:.2f}s, CIM {r.cim_time:.2f}s)")
        if r.extra:
            print(f"    {r.extra}")

    # ---- 3. Number Partitioning ----
    print()
    print("=" * 70)
    print("  Domain 3: Combinatorial Optimization — Number Partitioning")
    print("=" * 70)

    for N, seed in [(15, 100), (20, 200), (24, 300)]:
        rng = np.random.RandomState(seed)
        numbers = rng.randint(1, 100, size=N)
        J, a, desc = make_number_partition(numbers)

        def part_info(cim_s, exact_s, _a=a):
            s = f"CIM residue=|{int(np.dot(_a, cim_s))}|"
            if exact_s is not None:
                s += f", exact residue=|{int(np.dot(_a, exact_s))}|"
            return s

        r = run_problem(desc, J, extra_fn=part_info)
        results.append(r)
        nums_str = str([int(x) for x in numbers])
        print(f"  [{r.name}]  numbers={nums_str}")
        print(f"    exact_E={r.exact_energy:.2f}  CIM_E={r.cim_energy:.2f}  "
              f"gap={r.gap_pct:.2f}%")
        print(f"    {r.extra}")
        print(f"    (exact {r.exact_time:.2f}s, CIM {r.cim_time:.2f}s)")

    # ---- 4. Frustrated Triangular Antiferromagnet ----
    print()
    print("=" * 70)
    print("  Domain 4: Statistical Mechanics — Frustrated Antiferromagnet")
    print("=" * 70)

    for L in [3, 4, 5]:
        J, desc = make_triangular_antiferromagnet(L)
        r = run_problem(desc, J)
        results.append(r)
        print(f"  [{r.name}]  exact={r.exact_energy:.4f}  "
              f"CIM={r.cim_energy:.4f}  gap={r.gap_pct:.2f}%  "
              f"(exact {r.exact_time:.2f}s, CIM {r.cim_time:.2f}s)")

    # ---- 5. Weighted MAX-2-SAT ----
    print()
    print("=" * 70)
    print("  Domain 5: Operations Research — Weighted MAX-2-SAT")
    print("=" * 70)

    for nv, nc, seed in [(12, 30, 10), (16, 50, 20), (20, 80, 30)]:
        J, h, clauses, desc = make_weighted_max2sat(nv, nc, seed)
        total_weight = sum(w for *_, w in clauses)

        def sat_info(cim_s, exact_s, _cls=clauses, _tw=total_weight):
            cc, cw = satisfied_clauses(cim_s, _cls)
            s = f"CIM: {cc}/{len(_cls)} clauses (weight {cw:.1f}/{_tw:.1f})"
            if exact_s is not None:
                ec, ew = satisfied_clauses(exact_s, _cls)
                s += f" | exact: {ec}/{len(_cls)} (weight {ew:.1f}/{_tw:.1f})"
            return s

        r = run_problem(desc, J, h=h, extra_fn=sat_info)
        results.append(r)
        print(f"  [{r.name}]")
        print(f"    exact_E={r.exact_energy:.4f}  CIM_E={r.cim_energy:.4f}  "
              f"gap={r.gap_pct:.2f}%")
        print(f"    {r.extra}")
        print(f"    (exact {r.exact_time:.2f}s, CIM {r.cim_time:.2f}s)")

    # ---- Summary Table ----
    print()
    print("=" * 100)
    print(f"{'Problem':<45} {'N':>4} {'Exact E':>12} {'CIM E':>12} "
          f"{'Gap%':>8} {'Exact(s)':>9} {'CIM(s)':>8}")
    print("=" * 100)
    for r in results:
        exact_str = f"{r.exact_energy:.4f}" if not np.isnan(r.exact_energy) else "N/A"
        gap_str = f"{r.gap_pct:.2f}%" if not np.isnan(r.gap_pct) else "N/A"
        print(f"{r.name:<45} {r.N:>4} {exact_str:>12} {r.cim_energy:>12.4f} "
              f"{gap_str:>8} {r.exact_time:>9.2f} {r.cim_time:>8.2f}")
    print("=" * 100)

    optimal_count = sum(1 for r in results if abs(r.gap_pct) < 0.01)
    total = len(results)
    print(f"\nCIM found exact optimal: {optimal_count}/{total} problems")
    near_optimal = sum(1 for r in results if not np.isnan(r.gap_pct) and abs(r.gap_pct) < 1.0)
    print(f"CIM within 1% of optimal: {near_optimal}/{total} problems")


if __name__ == "__main__":
    run_all()
