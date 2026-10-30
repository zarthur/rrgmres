"""Range-restricted GMRES for square linear systems.

The iterates minimize ``||A x - b||_2`` over
``span{A b, A**2 b, ..., A**k b}``.  The implementation uses NumPy only and
accepts dense arrays, sparse/matrix-like objects implementing ``@``,
``scipy.sparse.linalg.LinearOperator``-style objects, or a callable matvec.

The historic MATLAB implementation is by A. Neuman, L. Reichel, and
H. Sadok (2011).  See the repository README for its retained license notice
and the paper reference.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Optional, Union

import numpy as np


MatvecOperator = Union[Any, Callable[[np.ndarray], np.ndarray]]


@dataclass(frozen=True)
class RRGMRESResult:
    """Result returned by :func:`rrgmres` and its stopping-rule wrappers.

    ``residual_norms`` contains the initial residual norm at index zero,
    followed by one value per completed outer iteration. Discrepancy runs
    explicitly verify changed iterates with ``A @ x`` before stopping. Fixed
    iteration runs report residuals from the reduced QR representation. If
    history was requested, ``x_history`` has the same convention and starts
    with x_0=0.

    ``converged`` means that the discrepancy threshold was reached.  A
    fixed-iteration run has no discrepancy target, so this flag is false even
    if its final residual happens to be small.  ``status`` distinguishes
    discrepancy termination, the iteration limit, and invariant-subspace
    (Arnoldi breakdown) termination.
    """

    x: np.ndarray
    residual_norms: np.ndarray
    iterations: int
    converged: bool
    status: str
    matvec_count: int
    x_history: Optional[np.ndarray] = None

    @property
    def residual_norm(self) -> float:
        """Return the final residual norm."""

        return float(self.residual_norms[-1])


def _vector(value: Any, name: str) -> np.ndarray:
    arr = np.asarray(value)
    if arr.ndim == 2 and 1 in arr.shape:
        arr = arr.reshape(-1)
    if arr.ndim != 1:
        raise ValueError("{} must be a one-dimensional vector".format(name))
    if arr.dtype.kind not in "fciu":
        raise TypeError("{} must contain real or complex numeric values".format(name))
    return arr


def _norm(vector: np.ndarray) -> float:
    """A scaled Euclidean norm that avoids overflow in sum(abs(x)**2)."""

    return float(np.hypot.reduce(np.abs(vector))) if vector.size else 0.0


def _matvec(operator: MatvecOperator, vector: np.ndarray, n: int) -> np.ndarray:
    if callable(operator) and not hasattr(operator, "__matmul__"):
        result = operator(vector)
    else:
        result = operator @ vector
    result = np.asarray(result)
    if result.ndim == 2 and 1 in result.shape:
        result = result.reshape(-1)
    if result.ndim != 1 or result.size != n:
        raise ValueError("A must map vectors of length {} to length {}".format(n, n))
    if result.dtype.kind not in "fciu":
        raise TypeError("A @ v must return real or complex numeric values")
    if not np.all(np.isfinite(result)):
        raise FloatingPointError("A @ v returned NaN or infinity")
    return result


def _orthogonalize(basis: np.ndarray, count: int, vector: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Two-pass modified Gram-Schmidt against a column-oriented basis."""

    work = np.array(vector, copy=True)
    coefficients = np.zeros(count, dtype=np.result_type(basis.dtype, vector.dtype))
    if count == 0:
        return work, coefficients
    q = basis[:, :count]
    # Two block projections are a vectorized, reorthogonalized MGS step.  The
    # second pass is important for clustered/ill-conditioned Krylov spaces.
    for _ in range(2):
        correction = q.conj().T @ work
        work -= q @ correction
        coefficients += correction
    return work, coefficients


def _back_substitute(upper: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """Solve a nonsingular upper-triangular system in O(k**2) work."""

    n = rhs.size
    solution = np.zeros(n, dtype=np.result_type(upper.dtype, rhs.dtype))
    for i in range(n - 1, -1, -1):
        solution[i] = (rhs[i] - upper[i, i + 1 : n] @ solution[i + 1 : n]) / upper[i, i]
    return solution


def _iteration_limit(value: Optional[int], n: int) -> int:
    if value is None:
        return n
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise TypeError("maxiter must be a non-negative integer")
    if value < 0:
        raise ValueError("maxiter must be non-negative")
    return min(int(value), n)


def rrgmres(
    A: MatvecOperator,
    b: Any,
    *,
    maxiter: Optional[int] = None,
    discrepancy: Optional[float] = None,
    eta: float = 1.0,
    store_history: bool = False,
    rank_tol: Optional[float] = None,
) -> RRGMRESResult:
    """Solve a square system with range-restricted GMRES.

    At outer iteration ``k``, this computes

    ``argmin ||A x - b||_2  subject to  x in span{A b, ..., A**k b}``.

    Parameters
    ----------
    A
        Square dense array, sparse/matrix-like object with ``A @ v``,
        matrix-free operator with ``A @ v``, or a callable ``A(v)``.
        Only matrix-vector products are required.
    b
        Right-hand side. A row or column vector is accepted and flattened.
    maxiter
        Maximum number of completed RRGMRES iterates. Defaults to the matrix
        dimension; values above the dimension are capped. A value of zero
        returns x_0=0 without applying A.
    discrepancy
        Optional bound ``delta`` for the data error. With this argument, stop
        at the first iterate satisfying ``||A x - b|| <= eta * delta``.
        The initial iterate is checked too. If omitted, exactly ``maxiter``
        steps are attempted unless the Krylov subspace breaks down.
    eta
        Discrepancy safety factor, positive and finite. Set to 1 to use
        ``discrepancy`` as the direct residual threshold.
    store_history
        If true, retain solution vectors. ``x_history[:, 0]`` is x_0 and
        later columns are completed iterates. Keeping this disabled avoids
        the O(n * maxiter) history allocation for large problems.
    rank_tol
        Optional relative threshold for declaring a new Krylov/QR direction
        numerically dependent. The default is ``10 * eps * n``.

    Returns
    -------
    RRGMRESResult
        Final solution, residual history, status, and number of operator
        applications. ``residual_norms`` starts with ``||b||``. Iterations
        count completed RRGMRES iterates, not Arnoldi/operator steps. In
        discrepancy mode, candidate residuals used for stopping are verified
        with an explicit operator application.

    Notes
    -----
    The method uses Arnoldi vectors for ``K_k(A,b)`` and incrementally
    orthogonalizes ``A v_k`` to build a basis for the required range-restricted
    search space. A second incrementally updated QR factorization represents
    ``A`` applied to that basis. Full-rank least-squares solves use triangular
    back-substitution; rank-deficient reduced systems fall back to
    ``numpy.linalg.lstsq``. No normal equations or explicit matrix inverse are
    formed.

    The dense working bases use O(n * maxiter) storage. For a dense matrix,
    each ordinary iteration performs an Arnoldi matvec and, when the search
    space grows, one additional matvec for its new range basis vector. A
    discrepancy run also verifies each changed solution with one additional
    matvec, so numerical rank truncation cannot falsely satisfy the threshold.
    """

    b_arr = _vector(b, "b")
    n = b_arr.size
    if n == 0:
        raise ValueError("b must not be empty")
    if not np.all(np.isfinite(b_arr)):
        raise ValueError("b must contain only finite values")
    shape = getattr(A, "shape", None)
    if shape is not None and tuple(shape) != (n, n):
        raise ValueError("A must be square with shape ({}, {})".format(n, n))

    steps = _iteration_limit(maxiter, n)
    if not np.isfinite(eta) or eta <= 0:
        raise ValueError("eta must be positive and finite")
    if discrepancy is not None:
        if not np.isfinite(discrepancy) or discrepancy < 0:
            raise ValueError("discrepancy must be non-negative and finite")
        target = float(eta * discrepancy)
    else:
        target = None
    if rank_tol is not None and (not np.isfinite(rank_tol) or rank_tol < 0):
        raise ValueError("rank_tol must be non-negative and finite")

    norm_b = _norm(b_arr)
    initial_dtype = np.result_type(b_arr.dtype, np.float64)
    x = np.zeros(n, dtype=initial_dtype)
    residual_history = [norm_b]
    solution_history = [x.copy()] if store_history else None

    if target is not None and norm_b <= target:
        return RRGMRESResult(x, np.asarray(residual_history), 0, True, "discrepancy", 0,
                             np.column_stack(solution_history) if store_history else None)
    if norm_b == 0:
        return RRGMRESResult(x, np.asarray(residual_history), 0, False, "zero_rhs", 0,
                             np.column_stack(solution_history) if store_history else None)
    if steps == 0:
        return RRGMRESResult(x, np.asarray(residual_history), 0, False, "iteration_limit", 0,
                             np.column_stack(solution_history) if store_history else None)

    v0 = b_arr / norm_b
    first_product = _matvec(A, v0, n)
    dtype = np.result_type(b_arr.dtype, first_product.dtype, np.float64)
    real_dtype = np.empty((), dtype=dtype).real.dtype
    eps = np.finfo(real_dtype).eps
    relative_rank_tol = float(rank_tol) if rank_tol is not None else 10.0 * eps * max(1, n)

    # Fortran-order storage makes each Krylov vector a contiguous column and
    # keeps the tall-thin BLAS projections efficient.
    V = np.zeros((n, steps + 1), dtype=dtype, order="F")
    U = np.zeros((n, steps), dtype=dtype, order="F")
    Q = np.zeros((n, steps), dtype=dtype, order="F")
    T = np.zeros((steps, steps), dtype=dtype, order="F")
    q_rhs = np.zeros(steps, dtype=dtype)
    V[:, 0] = v0

    v_count = 1
    u_count = 0
    q_count = 0
    matvec_count = 1
    rank_deficient = False
    current_residual = norm_b
    status = "iteration_limit"
    iterations = 0
    z = np.asarray(first_product, dtype=dtype)

    for k in range(steps):
        z_norm = _norm(z)

        # Stage 1: range(A V_k) = span{A b, ..., A**k b}. Two-pass
        # orthogonalization preserves the range-restricted nested spaces.
        u_candidate, _ = _orthogonalize(U, u_count, z)
        u_norm = _norm(u_candidate)
        added_u = u_norm > relative_rank_tol * z_norm
        if added_u:
            u = u_candidate / u_norm
            U[:, u_count] = u
            new_u = u_count
            u_count += 1

            # Stage 2: incremental QR for A U. The coefficient column is
            # formed in the existing Q basis, then a stable residual is added.
            au = _matvec(A, u, n).astype(dtype, copy=False)
            matvec_count += 1
            au_norm = _norm(au)
            q_candidate, coefficients = _orthogonalize(Q, q_count, au)
            T[:q_count, new_u] = coefficients
            q_norm = _norm(q_candidate)
            if q_norm > relative_rank_tol * au_norm:
                Q[:, q_count] = q_candidate / q_norm
                T[q_count, new_u] = q_norm
                q_rhs[q_count] = np.vdot(Q[:, q_count], b_arr)
                q_count += 1
            else:
                rank_deficient = True

            if q_count == u_count and not rank_deficient:
                y = _back_substitute(T[:u_count, :u_count], q_rhs[:u_count])
            elif q_count and u_count:
                y = np.linalg.lstsq(T[:q_count, :u_count], q_rhs[:q_count], rcond=None)[0]
            else:
                y = np.zeros(u_count, dtype=dtype)

            x = U[:, :u_count] @ y
            projected_ax = Q[:, :q_count] @ (T[:q_count, :u_count] @ y)
            current_residual = _norm(b_arr - projected_ax)
            if target is not None:
                # The reduced QR may intentionally discard a direction below
                # rank_tol. Verify the actual candidate before accepting a
                # discrepancy crossing, and include this operator call in the
                # public count. If no range direction was added, x is unchanged
                # and its last verified residual remains valid.
                true_ax = _matvec(A, x, n)
                matvec_count += 1
                current_residual = _norm(b_arr - true_ax)

        iterations = k + 1
        residual_history.append(current_residual)
        if solution_history is not None:
            solution_history.append(x.copy())

        if target is not None and current_residual <= target:
            status = "discrepancy"
            break

        # Arnoldi expansion of K_k(A,b), required even when this iteration's
        # range basis did not grow; a later Krylov vector can add a direction.
        v_candidate, _ = _orthogonalize(V, v_count, z)
        v_norm = _norm(v_candidate)
        arnoldi_breakdown = v_norm <= relative_rank_tol * z_norm
        if arnoldi_breakdown:
            status = "invariant_subspace"
            break

        V[:, v_count] = v_candidate / v_norm
        v_count += 1
        if k + 1 < steps:
            z = _matvec(A, V[:, k + 1], n).astype(dtype, copy=False)
            matvec_count += 1

    return RRGMRESResult(
        x=np.asarray(x),
        residual_norms=np.asarray(residual_history, dtype=float),
        iterations=iterations,
        converged=(status == "discrepancy"),
        status=status,
        matvec_count=matvec_count,
        x_history=np.column_stack(solution_history) if solution_history is not None else None,
    )


def rrgmres_dp(
    A: MatvecOperator,
    b: Any,
    discrepancy: float,
    *,
    eta: float = 1.0,
    maxiter: Optional[int] = None,
    store_history: bool = False,
) -> RRGMRESResult:
    """RRGMRES stopped by ``||A x - b|| <= eta * discrepancy``."""

    return rrgmres(A, b, maxiter=maxiter, discrepancy=discrepancy, eta=eta,
                   store_history=store_history)


def rrgmres_iter(
    A: MatvecOperator,
    b: Any,
    iterations: int,
    *,
    store_history: bool = False,
) -> RRGMRESResult:
    """Run RRGMRES for a fixed number of outer iterations (unless it breaks down)."""

    return rrgmres(A, b, maxiter=iterations, store_history=store_history)


def sym_rrgmres_dp(
    A: MatvecOperator,
    b: Any,
    discrepancy: float,
    *,
    eta: float = 1.0,
    maxiter: Optional[int] = None,
    store_history: bool = False,
) -> RRGMRESResult:
    """Symmetric-problem entry point using the stable shared Arnoldi engine.

    Symmetry is not required for correctness and is not inferred or checked.
    A short-recurrence Lanczos implementation could reduce storage, but this
    NumPy port shares the reorthogonalized engine to retain robust behavior.
    """

    return rrgmres_dp(A, b, discrepancy, eta=eta, maxiter=maxiter,
                      store_history=store_history)


def sym_rrgmres_iter(
    A: MatvecOperator,
    b: Any,
    iterations: int,
    *,
    store_history: bool = False,
) -> RRGMRESResult:
    """Symmetric-problem fixed-iteration entry point; see :func:`sym_rrgmres_dp`."""

    return rrgmres_iter(A, b, iterations, store_history=store_history)


__all__ = [
    "RRGMRESResult",
    "rrgmres",
    "rrgmres_dp",
    "rrgmres_iter",
    "sym_rrgmres_dp",
    "sym_rrgmres_iter",
]
