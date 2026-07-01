"""Benchmark Ising problem instance generators."""

import numpy as np
from pathlib import Path

INSTANCES_DIR = Path(__file__).resolve().parent.parent / "instances"


def random_ising(N: int, density: float = 1.0, seed: int = 42) -> np.ndarray:
    """Generate a random symmetric Ising coupling matrix J with entries in [-1, 1].

    Args:
        N: number of spins
        density: fraction of non-zero off-diagonal couplings (1.0 = fully connected)
        seed: random seed for reproducibility
    Returns:
        J: N x N symmetric matrix with zero diagonal
    """
    rng = np.random.RandomState(seed)
    J = rng.uniform(-1, 1, size=(N, N))
    J = (J + J.T) / 2
    np.fill_diagonal(J, 0)
    if density < 1.0:
        mask = rng.random((N, N)) < density
        mask = mask | mask.T
        np.fill_diagonal(mask, False)
        J *= mask
    return J


def random_maxcut(N: int, density: float = 0.5, seed: int = 42) -> np.ndarray:
    """Generate a random MAX-CUT instance as an Ising coupling matrix.

    For MAX-CUT, J_ij = w_ij/2 where w_ij are edge weights.
    Minimizing H = -sum J_ij sigma_i sigma_j maximizes the cut.

    Args:
        N: number of nodes
        density: edge probability
        seed: random seed
    Returns:
        J: N x N symmetric coupling matrix
    """
    rng = np.random.RandomState(seed)
    adj = (rng.random((N, N)) < density).astype(float)
    adj = np.triu(adj, k=1)
    adj = adj + adj.T
    J = -adj / 2
    return J


def sk_model(N: int, seed: int = 42) -> np.ndarray:
    """Generate a Sherrington-Kirkpatrick model instance.

    J_ij ~ N(0, 1/N) for i < j, symmetric, zero diagonal.
    """
    rng = np.random.RandomState(seed)
    J = rng.normal(0, 1.0 / np.sqrt(N), size=(N, N))
    J = np.triu(J, k=1)
    J = J + J.T
    return J


def save_instance(name: str, J: np.ndarray, h: np.ndarray = None):
    """Save an Ising instance to the instances directory."""
    INSTANCES_DIR.mkdir(parents=True, exist_ok=True)
    data = {"J": J}
    if h is not None:
        data["h"] = h
    np.savez(INSTANCES_DIR / f"{name}.npz", **data)


def load_instance(name: str):
    """Load an Ising instance from the instances directory."""
    data = np.load(INSTANCES_DIR / f"{name}.npz")
    J = data["J"]
    h = data.get("h", None)
    return J, h


def generate_benchmark_suite():
    """Generate a standard suite of benchmark instances."""
    configs = [
        ("random_N20_d100", 20, 1.0, 100),
        ("random_N50_d100", 50, 1.0, 200),
        ("random_N100_d100", 100, 1.0, 300),
        ("random_N200_d050", 200, 0.5, 400),
        ("random_N500_d050", 500, 0.5, 500),
        ("random_N1000_d010", 1000, 0.1, 600),
    ]
    for name, N, density, seed in configs:
        J = random_ising(N, density=density, seed=seed)
        save_instance(name, J)

    maxcut_configs = [
        ("maxcut_N20_d050", 20, 0.5, 700),
        ("maxcut_N50_d050", 50, 0.5, 800),
        ("maxcut_N100_d050", 100, 0.5, 900),
        ("maxcut_N200_d050", 200, 0.5, 1000),
        ("maxcut_N500_d020", 500, 0.2, 1100),
    ]
    for name, N, density, seed in maxcut_configs:
        J = random_maxcut(N, density=density, seed=seed)
        save_instance(name, J)

    sk_configs = [
        ("sk_N50", 50, 42),
        ("sk_N100", 100, 43),
        ("sk_N200", 200, 44),
        ("sk_N500", 500, 45),
    ]
    for name, N, seed in sk_configs:
        J = sk_model(N, seed=seed)
        save_instance(name, J)

    print(f"Generated {len(configs) + len(maxcut_configs) + len(sk_configs)} benchmark instances in {INSTANCES_DIR}")


if __name__ == "__main__":
    generate_benchmark_suite()
