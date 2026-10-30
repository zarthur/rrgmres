"""Small reproducible examples for fixed and discrepancy stopping."""

import numpy as np

from rrgmres import rrgmres_dp, rrgmres_iter


def main():
    A = np.diag([1.0, 2.0])
    b = np.ones(2)

    fixed = rrgmres_iter(A, b, iterations=1, store_history=True)
    print("fixed x:", np.array2string(fixed.x, precision=8))
    print("fixed residuals:", np.array2string(fixed.residual_norms, precision=8))
    print("fixed iterations:", fixed.iterations)
    print("fixed history shape:", fixed.x_history.shape)

    discrepancy = rrgmres_dp(A, b, discrepancy=0.75, maxiter=2)
    print("discrepancy:", discrepancy.iterations,
          "residual:", "{:.8f}".format(discrepancy.residual_norm),
          "status:", discrepancy.status,
          "matvecs:", discrepancy.matvec_count)


if __name__ == "__main__":
    main()
