# Range-restricted GMRES in NumPy

This Python port implements the range-restricted GMRES (RRGMRES) method from
the original MATLAB library in this repository. It supports discrepancy
principle stopping and a fixed iteration limit, and also accepts matrix-free
matvec operators.

## Mathematical problem

At iteration `k`, RRGMRES computes

```text
minimize ||A x - b||_2  over  x in span{A b, A^2 b, ..., A^k b}.
```

That range restriction is the defining property of RRGMRES. The implementation
starts an Arnoldi process with `b`, incrementally builds an orthonormal basis
for `range(A V_k)`, and solves the reduced least-squares problem for `A` applied
to that basis. It does not substitute ordinary GMRES, form normal equations, or
compute an inverse.

## Quick start

```python
import numpy as np
from rrgmres import rrgmres_dp, rrgmres_iter

A = np.diag(np.geomspace(1.0, 1e-6, 100))
x_true = np.ones(100)
b = A @ x_true

# Stop at the first iterate with residual <= eta * data_error_bound.
dp = rrgmres_dp(A, b, discrepancy=1e-4, eta=1.01, maxiter=40)

# Or request a fixed number of RRGMRES iterates.
fixed = rrgmres_iter(A, b, iterations=12, store_history=True)

print(dp.iterations, dp.residual_norm, dp.status)
print(fixed.x_history.shape)  # (n, iterations + 1), including x_0
```

The main `rrgmres(A, b, ...)` entry point accepts either stopping option.
`RRGMRESResult.residual_norms` starts with the initial residual `||b||`, then
contains one residual per completed RRGMRES iteration. `iterations` counts
completed outer iterates, not matrix-vector products. `matvec_count` reports
operator applications. A fixed-iteration run returns `status="iteration_limit"`
unless the Krylov sequence reaches an invariant subspace first. For discrepancy
runs, `converged` is true exactly when the specified threshold is reached.

## Operators and dtypes

`A` may be a NumPy array, a sparse or other matrix-like object implementing
`A @ v`, a `LinearOperator`-style object, or a callable `A(v)`. Only square
systems are supported. Vector inputs may be one-dimensional or row/column
vectors. Integer and single-precision inputs are promoted to at least float64;
complex inputs remain complex. `store_history=False` avoids retaining every
solution vector.

## Numerical behavior and limits

The Arnoldi, range-basis, and reduced least-squares bases use two-pass
reorthogonalization. Full-rank reduced systems are solved by triangular
back-substitution. If a reduced column is dependent, the implementation uses
`numpy.linalg.lstsq` on the small projected system to handle singular and
rank-deficient cases without normal equations. Fixed-iteration residuals come
from the reduced QR representation. Discrepancy runs verify each changed
candidate with a fresh `A @ x` before accepting a threshold crossing, so a
direction dropped at the numerical rank threshold cannot produce a false
convergence result. These verification products are included in
`matvec_count`.

The implementation stores O(n k) basis vectors and normally uses two operator
applications per fixed iteration: one Arnoldi product and one product for a new
range-basis vector. Discrepancy mode adds a verification product whenever the
candidate changes. It accepts sparse/matrix-free operators, but each vector
basis is dense. The historical symmetric entry points share this reorthogonalized
Arnoldi engine for robustness; they do not use the original short-recurrence
Lanczos storage reduction.

Compared with the MATLAB routines, the Python API checks the initial residual,
handles zero right-hand sides, detects Arnoldi breakdown, and returns structured
status/history instead of relying on preallocated output shapes. The original
discrepancy code's residual variable was scalar despite its help text describing
a history; the Python result returns a consistent history including `x_0`.

## Differences from the MATLAB API

The MATLAB routines return multiple positional outputs (`X`, residual data,
and sometimes an iteration count). Python returns one `RRGMRESResult`: `.x` is
the final solution, `.residual_norms` includes the initial residual, and
`.x_history` is optional and includes `x_0` in its first column. `maxiter`
counts completed RRGMRES iterates; the MATLAB implementation internally took
one additional Arnoldi step to form each iterate. Discrepancy stopping accepts
equality (`<= eta * discrepancy`) and checks the initial residual. If an
invariant subspace is reached, Python returns the best computed iterate and
reports `status="invariant_subspace"` instead of dividing by a zero Arnoldi
norm.

`examples/basic_usage.py` is a runnable example. Its captured output is in
`examples/output.txt`.

## Validation and benchmark

Create an environment and run the tests and benchmark from this directory:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install .
.venv/bin/python -m unittest discover -s tests -v
.venv/bin/python benchmarks/benchmark_rrgmres.py
```

The tests compare iterates against an independent QR-plus-least-squares
construction of the RRGMRES subspace, check true residuals, and cover
discrepancy thresholds, near-rank-threshold residual verification, a multi-column
rank-deficient reduced system, breakdown, singular systems, complex/integer
inputs, scaling, and matrix-free operator calls. A captured usage example is in
`examples/output.txt`.

## Reference

A. Neuman, L. Reichel, and H. Sadok, “Algorithms for range restricted iterative
methods for linear discrete ill-posed problems,” *Numerical Algorithms*,
59 (2012), 325–331, doi:10.1007/s11075-011-9491-4. The implementation follows
the range-restricted minimization problem in that paper and the MATLAB code in
this repository.

## Original software notice

This library retains the permission and warranty notice carried in the original
README. See [LICENSE.txt](LICENSE.txt). Original MATLAB software copyright
(c) 2011, A. Neuman, L. Reichel, and H. Sadok.
