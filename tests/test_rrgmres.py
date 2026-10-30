import unittest
from unittest import mock

import numpy as np

from rrgmres import (
    rrgmres,
    rrgmres_dp,
    rrgmres_iter,
    sym_rrgmres_dp,
    sym_rrgmres_iter,
)


def _orthogonalize(columns, vector):
    work = vector.copy()
    for _ in range(2):
        for q in columns:
            work -= q * np.vdot(q, work)
    return work


def direct_reference(A, b, maxiter):
    """Independent RRGMRES oracle: explicit search QR and fresh least squares."""
    b = np.asarray(b, dtype=np.result_type(np.asarray(b).dtype, np.float64)).reshape(-1)
    beta = np.linalg.norm(b)
    if beta == 0:
        return [np.zeros_like(b)]
    V = [b / beta]
    U = []
    iterates = []
    for k in range(maxiter):
        z = A @ V[k]
        u = _orthogonalize(U, z)
        unorm = np.linalg.norm(u)
        if unorm > 1e-13 * np.linalg.norm(z):
            U.append(u / unorm)
            Q = np.column_stack(U)
            y = np.linalg.lstsq(A @ Q, b, rcond=None)[0]
            iterates.append(Q @ y)
        elif iterates:
            iterates.append(iterates[-1].copy())
        else:
            iterates.append(np.zeros_like(b))

        v = _orthogonalize(V, z)
        vnorm = np.linalg.norm(v)
        if vnorm <= 1e-13 * np.linalg.norm(z):
            break
        V.append(v / vnorm)
    return iterates


class RRGMRESTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(1205)
        self.A = rng.normal(size=(9, 9)) + 0.2 * np.triu(rng.normal(size=(9, 9)), 1)
        self.b = rng.normal(size=9)

    def test_iterates_match_independent_range_restricted_lstsq(self):
        expected = direct_reference(self.A, self.b, 5)
        result = rrgmres_iter(self.A, self.b, 5, store_history=True)
        self.assertEqual(result.iterations, 5)
        self.assertEqual(result.x_history.shape, (9, 6))
        for k, expected_x in enumerate(expected, start=1):
            np.testing.assert_allclose(result.x_history[:, k], expected_x, rtol=2e-10, atol=2e-11)
            true_residual = np.linalg.norm(self.b - self.A @ expected_x)
            self.assertAlmostEqual(result.residual_norms[k], true_residual, places=10)
        np.testing.assert_allclose(result.x, expected[-1], rtol=2e-10, atol=2e-11)

    def test_search_space_is_range_restricted_not_ordinary_gmres(self):
        A = np.diag([1.0, 2.0])
        b = np.ones(2)
        result = rrgmres_iter(A, b, 1)
        np.testing.assert_allclose(result.x, [5.0 / 17.0, 10.0 / 17.0], atol=1e-12)
        self.assertNotAlmostEqual(result.x[0], result.x[1])

    def test_discrepancy_stops_at_first_crossing(self):
        reference = direct_reference(self.A, self.b, 5)
        residuals = [np.linalg.norm(self.b - self.A @ x) for x in reference]
        threshold = 0.5 * (residuals[0] + residuals[1])
        result = rrgmres_dp(self.A, self.b, threshold, maxiter=5)
        self.assertTrue(result.converged)
        self.assertEqual(result.iterations, 2)
        self.assertLessEqual(result.residual_norm, threshold)
        self.assertGreater(result.residual_norms[1], threshold)

    def test_initial_discrepancy_including_equality_and_zero(self):
        A = np.eye(3)
        b = np.array([1.0, -2.0, 3.0])
        result = rrgmres_dp(A, b, np.linalg.norm(b), maxiter=2)
        self.assertEqual(result.iterations, 0)
        self.assertTrue(result.converged)
        self.assertEqual(result.status, "discrepancy")
        result_zero = rrgmres_dp(np.eye(2), np.zeros(2), 0.0, maxiter=2)
        self.assertEqual(result_zero.iterations, 0)
        self.assertTrue(result_zero.converged)
        exact = rrgmres_dp(A, b, 0.0, maxiter=2)
        self.assertEqual(exact.iterations, 1)
        self.assertTrue(exact.converged)
        self.assertEqual(exact.residual_norm, 0.0)

    def test_zero_budget_and_invalid_iteration_limits(self):
        result = rrgmres(self.A, self.b, maxiter=0, store_history=True)
        self.assertEqual(result.iterations, 0)
        self.assertEqual(result.matvec_count, 0)
        np.testing.assert_array_equal(result.x, np.zeros_like(self.b))
        self.assertEqual(result.x_history.shape, (9, 1))
        with self.assertRaises(ValueError):
            rrgmres(self.A, self.b, maxiter=-1)
        with self.assertRaises(TypeError):
            rrgmres(self.A, self.b, maxiter=1.5)

    def test_zero_rhs_and_zero_matrix(self):
        zero_b = rrgmres(np.eye(3), np.zeros(3), maxiter=2)
        self.assertEqual(zero_b.status, "zero_rhs")
        self.assertEqual(zero_b.iterations, 0)
        zero_A = rrgmres(np.zeros((3, 3)), np.ones(3), maxiter=3)
        self.assertEqual(zero_A.status, "invariant_subspace")
        self.assertEqual(zero_A.iterations, 1)
        self.assertAlmostEqual(zero_A.residual_norm, np.sqrt(3))

    def test_singular_inconsistent_and_nilpotent_breakdown(self):
        diagonal = np.diag([1.0, 0.0])
        result = rrgmres_iter(diagonal, np.array([1.0, 1.0]), 5, store_history=True)
        self.assertEqual(result.status, "invariant_subspace")
        self.assertAlmostEqual(result.residual_norm, 1.0, places=12)
        np.testing.assert_allclose(result.x, np.array([1.0, 0.0]), atol=1e-12)
        self.assertAlmostEqual(np.linalg.norm(diagonal @ result.x - [1.0, 1.0]), 1.0, places=12)

        nilpotent = np.array([[0.0, 1.0], [0.0, 0.0]])
        nil_result = rrgmres_iter(nilpotent, np.array([0.0, 1.0]), 4)
        self.assertEqual(nil_result.status, "invariant_subspace")
        self.assertAlmostEqual(nil_result.residual_norm, 1.0, places=12)

    def test_multi_column_rank_deficient_reduced_operator_fallback(self):
        # The search basis has two independent directions, while A applied to
        # those directions has rank one. This exercises the lstsq fallback for
        # a genuinely multi-column rank-deficient projected operator.
        A = np.array([[0.0, 1.0, 0.0],
                      [0.0, 0.0, 1.0],
                      [0.0, 0.0, 0.0]])
        b = np.array([1.0, 2.0, 3.0])
        expected = direct_reference(A, b, 3)
        with mock.patch("numpy.linalg.lstsq", wraps=np.linalg.lstsq) as lstsq:
            result = rrgmres_iter(A, b, 3, store_history=True)
        self.assertGreaterEqual(lstsq.call_count, 1)
        self.assertEqual(result.status, "invariant_subspace")
        for k, expected_x in enumerate(expected, start=1):
            np.testing.assert_allclose(result.x_history[:, k], expected_x,
                                       rtol=2e-10, atol=2e-11)
            self.assertAlmostEqual(
                result.residual_norms[k], np.linalg.norm(b - A @ expected_x), places=11
            )

    def test_discrepancy_verifies_true_residual_after_rank_truncation(self):
        # With this rank_tol, the reduced QR estimate at iteration 3 is about
        # 1.79086 while the actual residual is about 1.79959. Earlier true
        # residuals exceed the midpoint threshold, so only an explicit check
        # prevents a false discrepancy stop at iteration 3.
        A = np.array([
            [1.0094856917558223e-1, 2.0765360379822578e-1,
             -3.7677936827534281e-1, -2.2578624544235731e-1],
            [-6.3185644745677519e-2, 1.1512022147562230e-2,
             -1.9200698124635651e-4, -1.5871600408605357e-1],
            [3.5096052026639923e-1, 5.8548569050767815e-2,
             -2.3645375310496769e-1, 4.8143997756804807e-1],
            [-2.9514837319659581e-1, -1.6171008199996698e-1,
             3.0952556345745502e-1, -6.0102744540902175e-1],
        ])
        b = np.array([0.7981906868272103, 2.0365058013148760,
                      -0.14057751100983934, 0.8269254800023204])
        rank_tol = 0.2
        fixed = rrgmres(A, b, maxiter=4, rank_tol=rank_tol, store_history=True)
        actual = np.array([
            np.linalg.norm(b - A @ fixed.x_history[:, k])
            for k in range(fixed.iterations + 1)
        ])
        mismatch_k = 3
        projected = fixed.residual_norms[mismatch_k]
        self.assertGreater(actual[mismatch_k] - projected, 0.005)
        threshold = 0.5 * (projected + actual[mismatch_k])
        self.assertTrue(np.all(actual[:mismatch_k] > threshold))

        calls = []

        def apply(v):
            calls.append(1)
            return A @ v

        result = rrgmres(apply, b, maxiter=4, discrepancy=threshold,
                         rank_tol=rank_tol, store_history=True)
        self.assertFalse(result.converged)
        self.assertGreater(result.residual_norm, threshold)
        self.assertEqual(result.iterations, fixed.iterations)
        self.assertEqual(result.matvec_count, len(calls))
        for k in range(result.iterations + 1):
            true_residual = np.linalg.norm(b - A @ result.x_history[:, k])
            self.assertAlmostEqual(result.residual_norms[k], true_residual, places=11)

    def test_integer_and_complex_inputs(self):
        integer_result = rrgmres_iter(np.eye(3, dtype=int), np.array([1, 2, 3]), 2)
        self.assertEqual(integer_result.x.dtype, np.dtype(float))
        np.testing.assert_allclose(integer_result.x, [1, 2, 3], atol=1e-12)

        A = np.array([[2.0, 1.0j], [-1.0j, 3.0]], dtype=complex)
        b = np.array([1.0 + 2.0j, -0.5j])
        expected = direct_reference(A, b, 2)
        result = rrgmres_iter(A, b, 2)
        self.assertEqual(result.x.dtype, np.dtype(complex))
        np.testing.assert_allclose(result.x, expected[-1], rtol=2e-10, atol=2e-11)

    def test_matrix_free_operator_and_matvec_count(self):
        A = np.array([[2.0, 1.0], [0.0, 3.0]])
        calls = []

        def apply(v):
            calls.append(1)
            return A @ v

        result = rrgmres(apply, np.array([1.0, 2.0]), maxiter=1)
        self.assertEqual(result.iterations, 1)
        self.assertEqual(result.matvec_count, len(calls))
        self.assertEqual(result.matvec_count, 2)

    def test_scaled_matrix(self):
        A = 1.0e-100 * np.eye(3)
        b = np.array([1.0, -2.0, 0.5])
        result = rrgmres_dp(A, b, discrepancy=1e-12, maxiter=2)
        self.assertTrue(result.converged)
        np.testing.assert_allclose(A @ result.x, b, rtol=1e-12, atol=1e-12)
        large_result = rrgmres_iter(1.0e100 * np.eye(3), b, 2)
        np.testing.assert_allclose((1.0e100 * np.eye(3)) @ large_result.x, b,
                                   rtol=1e-12, atol=1e-12)

    def test_square_shape_and_threshold_validation(self):
        with self.assertRaises(ValueError):
            rrgmres(np.ones((2, 3)), np.ones(2), maxiter=1)
        with self.assertRaises(ValueError):
            rrgmres(self.A, self.b, discrepancy=-1.0)
        with self.assertRaises(ValueError):
            rrgmres(self.A, self.b, eta=0.0)

    def test_symmetric_entry_points_share_same_minimization(self):
        A = np.diag([1.0, 0.2, 0.01])
        b = np.array([1.0, 2.0, 3.0])
        general = rrgmres_iter(A, b, 2, store_history=True)
        symmetric = sym_rrgmres_iter(A, b, 2, store_history=True)
        np.testing.assert_allclose(general.x, symmetric.x)
        dp = sym_rrgmres_dp(A, b, 1e-6, maxiter=3)
        self.assertTrue(dp.converged)


if __name__ == "__main__":
    unittest.main()
