"""Deterministic comparison with a direct reduced least-squares reference.

Run from this project directory with ``python benchmarks/benchmark_rrgmres.py``.
This is a local microbenchmark, not a comparison against MATLAB or a claim of
hardware-independent speedup.
"""

import sys
import time
from pathlib import Path

import numpy as np

# Also allow running this script directly from a clean source checkout without
# first installing the local package.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from rrgmres import rrgmres_iter


def direct_reference(A, b, maxiter):
    """Rebuild each range basis and solve each full projected LS problem."""
    n = b.size
    v = b / np.linalg.norm(b)
    V = np.zeros((n, maxiter + 1), dtype=float)
    V[:, 0] = v
    Z = np.zeros((n, maxiter), dtype=float)
    residuals = []
    for k in range(maxiter):
        z = A @ V[:, k]
        Z[:, k] = z
        q, r = np.linalg.qr(Z[:, : k + 1], mode="reduced")
        y = np.linalg.lstsq(A @ q, b, rcond=None)[0]
        x = q @ y
        residuals.append(np.linalg.norm(b - A @ x))
        w = z.copy()
        for _ in range(2):
            w -= V[:, : k + 1] @ (V[:, : k + 1].T @ w)
        norm_w = np.linalg.norm(w)
        if norm_w <= 1e-14 * np.linalg.norm(z):
            break
        V[:, k + 1] = w / norm_w
    return np.asarray(residuals)


def main():
    rng = np.random.default_rng(20261030)
    n, maxiter = 240, 18
    left, _ = np.linalg.qr(rng.normal(size=(n, n)))
    right, _ = np.linalg.qr(rng.normal(size=(n, n)))
    singular_values = np.geomspace(1.0, 1e-5, n)
    A = (left * singular_values) @ right.T
    truth = rng.normal(size=n)
    b = A @ truth + 1e-4 * rng.normal(size=n)

    # Warm the BLAS and NumPy paths before timing.
    rrgmres_iter(A, b, 2)
    direct_reference(A, b, 2)

    fast_times = []
    reference_times = []
    result = None
    reference_residuals = None
    for _ in range(5):
        start = time.perf_counter()
        result = rrgmres_iter(A, b, maxiter)
        fast_times.append(time.perf_counter() - start)

        start = time.perf_counter()
        reference_residuals = direct_reference(A, b, maxiter)
        reference_times.append(time.perf_counter() - start)
    fast_time = float(np.median(fast_times))
    reference_time = float(np.median(reference_times))

    print("problem: dense {}x{}, requested iterates {}".format(n, n, maxiter))
    print("incremental implementation: median {:.6f}s, {} iterates, {} matvecs, residual {:.6e}".format(
        fast_time, result.iterations, result.matvec_count, result.residual_norm))
    print("recomputed QR/lstsq reference: median {:.6f}s, {} iterates, residual {:.6e}".format(
        reference_time, reference_residuals.size, reference_residuals[-1]))
    print("timing ratio (reference / incremental): {:.2f}x".format(reference_time / fast_time))


if __name__ == "__main__":
    main()
