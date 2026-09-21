import itertools
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from scipy.optimize import nnls as scipy_nnls
from threadpoolctl import threadpool_info, threadpool_limits

from depth_nnls import (
    _adelie_nnls, _solve_working_set, aggregate_identical_depth_rows, certify_depth_nnls,
    predict_depth_from_weights, solve_depth_nnls,
)


def enumerated_nnls(matrix, target):
    """Independent small-problem reference: enumerate every active face."""
    best = np.zeros(matrix.shape[1])
    best_sse = target @ target
    for size in range(1, matrix.shape[1] + 1):
        for indices in itertools.combinations(range(matrix.shape[1]), size):
            values = np.linalg.lstsq(matrix[:, indices], target, rcond=None)[0]
            if np.any(values < 0):
                continue
            candidate = np.zeros(matrix.shape[1])
            candidate[list(indices)] = values
            sse = np.linalg.norm(matrix @ candidate - target)**2
            if sse < best_sse:
                best, best_sse = candidate, sse
    return best


class DepthNnlsTests(unittest.TestCase):
    def test_noisy_repeated_rows_preserve_full_objective(self):
        matrix = np.array([[1, 0], [1, 0], [0, -1], [0, -1], [0, -1], [1, 1]], dtype=float)
        target = np.array([2, 4, -2, -1, 0, 4], dtype=float)
        reduced, reduced_target = aggregate_identical_depth_rows(matrix, target)
        constant = target @ target - reduced_target @ reduced_target
        for weights in (np.zeros(2), np.array([3, 1]), np.array([0.2, 7])):
            self.assertAlmostEqual(np.linalg.norm(matrix @ weights-target)**2,
                                   np.linalg.norm(reduced @ weights-reduced_target)**2 + constant, places=11)
            original_kkt = certify_depth_nnls(matrix, target, weights)["kkt_scaled_max"]
            aggregated_kkt = certify_depth_nnls(
                reduced, reduced_target, weights, target_norm=np.linalg.norm(target)
            )["kkt_scaled_max"]
            self.assertAlmostEqual(aggregated_kkt, original_kkt, places=14)
        fit, info = solve_depth_nnls(matrix, target, return_diagnostics=True)
        reference = enumerated_nnls(matrix, target)
        np.testing.assert_allclose(matrix @ fit, matrix @ reference, atol=1e-12)
        self.assertEqual(info["aggregated_rows"], 3)
        self.assertIsNone(info["kkt_matrix"])
        self.assertIsNone(info["kkt_scaled_max"])
        residual = matrix @ fit - target
        self.assertAlmostEqual(info["sse"], residual @ residual, places=12)
        self.assertAlmostEqual(info["relative_error"], np.linalg.norm(residual) / np.linalg.norm(target), places=12)

    def test_nested_correlated_features_cannot_worsen_raw_fit(self):
        rng = np.random.default_rng(57)
        base = rng.uniform(0.2, 1.5, (25, 3))
        perturbation = rng.normal(0, 0.001, 25)
        larger = np.column_stack([base, base[:, 0] + perturbation, np.zeros(25), base[:, 1]])
        target = 2 * larger[:, 3] + 0.4 * base[:, 1] + rng.normal(0, 0.0001, 25)
        low, low_info = solve_depth_nnls(base, target, return_diagnostics=True)
        high, high_info = solve_depth_nnls(larger, target, return_diagnostics=True)
        reference = enumerated_nnls(larger, target)
        self.assertLessEqual(high_info["sse"], low_info["sse"] + 1e-12)
        np.testing.assert_allclose(larger @ high, larger @ reference, atol=1e-9)
        self.assertEqual(high[4], 0)
        self.assertTrue(np.all(high >= 0))

    def test_full_scipy_result_is_returned_without_a_kkt_recheck(self):
        matrix = np.repeat(np.eye(2), 3, axis=0)
        target = np.array([1., 2., 3., 2., 3., 4.])
        # An inaccurate, finite nonnegative solver result is accepted as returned.
        returned = np.zeros(2)
        with patch("depth_nnls.nnls", return_value=(returned, np.linalg.norm(target))) as solver, \
             patch("depth_nnls.certify_depth_nnls", side_effect=AssertionError("redundant KKT scan")):
            weights, info = solve_depth_nnls(matrix, target, return_diagnostics=True)
        np.testing.assert_array_equal(weights, returned)
        self.assertEqual(solver.call_count, 1)
        self.assertEqual(info["strategy"], "full")
        self.assertEqual(info["convergence_check"], "solver_default")
        self.assertIsNone(info["kkt_scaled_max"])
        self.assertIsNone(info["kkt_tolerance"])
        self.assertIsNone(info["kkt_matrix"])
        self.assertIsNone(info["converged"])

    def test_working_set_reuses_convergence_without_another_certificate(self):
        matrix = np.repeat(np.array([[1., -1.], [0., 1.]]), 3, axis=0)
        target = np.array([0.9, 1., 1.1, 0.8, 1., 1.2])
        with patch("depth_nnls.FULL_SOLVE_COLUMN_THRESHOLD", 0), \
             patch("depth_nnls.certify_depth_nnls", side_effect=AssertionError("redundant KKT scan")):
            weights, info = solve_depth_nnls(matrix, target, return_diagnostics=True)
        reference = enumerated_nnls(matrix, target)
        np.testing.assert_allclose(matrix @ weights, matrix @ reference, atol=1e-12)
        self.assertEqual(info["strategy"], "working_set")
        self.assertEqual(info["kkt_matrix"], "aggregated")
        self.assertEqual(info["kkt_tolerance"], 1e-12)
        self.assertEqual(info["convergence_check"], "aggregated_kkt")
        self.assertIs(info["converged"], True)
        self.assertLessEqual(info["kkt_scaled_max"], 1e-12)
        self.assertAlmostEqual(info["sse"], np.linalg.norm(matrix @ weights - target)**2, places=12)

    def test_adelie_fallback_uses_aggregated_rows_without_kkt_recheck(self):
        matrix = np.repeat(np.array([[1., -1.], [0., 1.]]), 3, axis=0)
        target = np.array([0.9, 1., 1.1, 0.8, 1., 1.2])

        def force_fallback(*args, **kwargs):
            return _solve_working_set(*args, **kwargs, max_rounds=0)

        with patch("depth_nnls.FULL_SOLVE_COLUMN_THRESHOLD", 0), \
             patch("depth_nnls._solve_working_set", side_effect=force_fallback), \
             patch("depth_nnls.nnls", side_effect=AssertionError("full SciPy fallback")), \
             patch("depth_nnls.certify_depth_nnls", side_effect=AssertionError("fallback KKT recheck")), \
             patch("depth_nnls._adelie_nnls", wraps=_adelie_nnls) as fallback:
            weights, info = solve_depth_nnls(matrix, target, return_diagnostics=True)
        reference = enumerated_nnls(matrix, target)
        np.testing.assert_allclose(matrix @ weights, matrix @ reference, atol=1e-3)
        self.assertEqual(info["strategy"], "adelie_fallback")
        self.assertEqual(info["solver"], "adelie.solver.bvls")
        self.assertEqual(info["solver_tolerance"], 1e-7)
        self.assertEqual(info["convergence_check"], "solver_default")
        self.assertIsNone(info["converged"])
        self.assertIsNone(info["kkt_matrix"])
        self.assertIsNone(info["kkt_scaled_max"])
        self.assertIsNone(info["kkt_tolerance"])
        self.assertEqual(fallback.call_count, 1)
        self.assertEqual(fallback.call_args.args[0].shape, (2, 2))
        self.assertAlmostEqual(info["sse"], np.linalg.norm(matrix @ weights - target)**2, places=12)

    def test_adelie_result_above_strict_kkt_tolerance_is_accepted(self):
        matrix = np.eye(2)
        target = np.array([2., 3.])
        approximate = np.array([1.8, 2.7])
        self.assertGreater(certify_depth_nnls(matrix, target, approximate)["kkt_scaled_max"], 1e-8)
        with patch("depth_nnls.FULL_SOLVE_COLUMN_THRESHOLD", 0), \
             patch("depth_nnls.nnls", side_effect=RuntimeError("iteration limit")) as scipy_solve, \
             patch("depth_nnls._adelie_nnls", return_value=approximate), \
             patch("depth_nnls.certify_depth_nnls", side_effect=AssertionError("fallback KKT recheck")):
            weights, info = solve_depth_nnls(matrix, target, return_diagnostics=True)
        np.testing.assert_array_equal(weights, approximate)
        self.assertEqual(scipy_solve.call_count, 1)
        self.assertEqual(info["strategy"], "adelie_fallback")
        self.assertEqual(info["working_set_fallback"], "subproblem exceeded its iteration limit")
        self.assertIsNone(info["converged"])

    def test_invalid_adelie_results_are_rejected(self):
        matrix = np.eye(2)
        target = np.array([2., 3.])
        for values in ([np.nan, 3.], [np.inf, 3.], [-1., 3.]):
            with self.subTest(values=values), \
                 patch("depth_nnls.FULL_SOLVE_COLUMN_THRESHOLD", 0), \
                 patch("depth_nnls.nnls", side_effect=RuntimeError("iteration limit")), \
                 patch("depth_nnls._adelie_nnls", return_value=np.array(values)):
                with self.assertRaisesRegex(RuntimeError, "non-finite or negative coefficients"):
                    solve_depth_nnls(matrix, target)

    def test_adelie_reported_error_is_not_silently_accepted(self):
        state = SimpleNamespace(beta=np.array([2., 3.]), error="iteration limit reached")
        with patch("adelie.solver.bvls", return_value=state):
            with self.assertRaisesRegex(RuntimeError, "Adelie NNLS fallback failed: iteration limit reached"):
                _adelie_nnls(np.eye(2), np.array([2., 3.]))

    def test_zero_target_and_zero_columns(self):
        fit, info = solve_depth_nnls(np.zeros((4, 3)), np.zeros(4), return_diagnostics=True)
        np.testing.assert_array_equal(fit, np.zeros(3))
        self.assertEqual(info["sse"], 0)
        self.assertIsNone(info["kkt_scaled_max"])

    def test_invalid_inputs_do_not_reach_solver(self):
        for matrix, target in [(np.ones((2, 2)), [1]), (np.empty((0, 2)), []),
                               (np.array([[np.nan]]), [1]), (np.ones((1, 1)), [np.inf])]:
            with self.subTest(matrix=matrix, target=target):
                with self.assertRaises(ValueError):
                    solve_depth_nnls(matrix, target)

    def test_column_with_no_initial_violation_can_enter_later(self):
        matrix = np.array([[1., -1.], [0., 1.]])
        target = np.array([1., 1.])
        # The second column has zero initial correlation with the target,
        # but its coefficient must be positive at the full optimum.
        self.assertEqual(matrix[:, 1] @ target, 0)
        weights, info = _solve_working_set(matrix, target, np.linalg.norm(target), batch_size=1)
        np.testing.assert_allclose(weights, [2., 1.], atol=1e-12)
        self.assertEqual(info["strategy"], "working_set")
        self.assertEqual(info["working_set_rounds"], 2)

    def test_working_set_matches_independent_active_face_enumeration(self):
        rng = np.random.default_rng(74)
        matrix = rng.normal(size=(9, 6))
        matrix = np.column_stack([matrix, matrix[:, 0], np.zeros(9)])
        target = rng.normal(size=9)
        weights, info = _solve_working_set(matrix, target, np.linalg.norm(target), batch_size=1)
        expected = enumerated_nnls(matrix, target)
        np.testing.assert_allclose(matrix @ weights, matrix @ expected, atol=1e-10)
        self.assertLess(certify_depth_nnls(matrix, target, weights)["kkt_scaled_max"], 1e-10)
        self.assertEqual(info["strategy"], "working_set")

    def test_working_round_limit_falls_back_to_adelie(self):
        matrix = np.array([[1., -1., 0.], [0., 1., 1.], [1., 0., -1.]])
        target = np.array([1., 2., 0.])
        with patch("depth_nnls.nnls", side_effect=AssertionError("full SciPy fallback")):
            weights, info = _solve_working_set(matrix, target, np.linalg.norm(target), max_rounds=0)
        expected = enumerated_nnls(matrix, target)
        np.testing.assert_allclose(matrix @ weights, matrix @ expected, atol=1e-3)
        self.assertEqual(info["strategy"], "adelie_fallback")
        self.assertEqual(info["working_set_fallback"], "working-set round limit")

    def test_inaccurate_subproblem_is_replaced_by_adelie(self):
        matrix = np.eye(3)
        target = np.array([2., 3., 4.])
        with patch("depth_nnls.FULL_SOLVE_COLUMN_THRESHOLD", 0), \
             patch("depth_nnls.nnls", return_value=(np.zeros(3), np.linalg.norm(target))) as scipy_solve:
            weights, info = solve_depth_nnls(matrix, target, return_diagnostics=True)
        np.testing.assert_allclose(weights, target, atol=1e-12)
        self.assertEqual(scipy_solve.call_count, 1)
        self.assertEqual(info["strategy"], "adelie_fallback")
        self.assertEqual(info["working_set_fallback"], "subproblem left an admitted-column optimality violation")

    def test_sparse_weight_prediction_preserves_original_feature_positions(self):
        matrix = np.asarray(np.arange(5*520).reshape(5, 520), dtype=np.float32, order="F")
        weights = np.zeros(520)
        weights[[0, 257, 519]] = [1., 2., 3.]
        expected = matrix[:, 0].astype(float) + 2*matrix[:, 257] + 3*matrix[:, 519]
        np.testing.assert_array_equal(predict_depth_from_weights(matrix, weights), expected)
        np.testing.assert_array_equal(predict_depth_from_weights(matrix, np.zeros(520)), np.zeros(5))
        self.assertEqual(predict_depth_from_weights(matrix[:0], weights).shape, (0,))

    def test_solver_respects_callers_blas_limit(self):
        observed = []
        def capture(*args, **kwargs):
            observed.extend(r["num_threads"] for r in threadpool_info() if r["user_api"] == "blas")
            return scipy_nnls(*args, **kwargs)
        with threadpool_limits(limits=4, user_api="blas"):
            before = [(r["filepath"], r["num_threads"]) for r in threadpool_info() if r["user_api"] == "blas"]
            with patch("depth_nnls.nnls", side_effect=capture):
                weights = solve_depth_nnls(np.eye(2), np.array([2., 3.]))
            after = [(r["filepath"], r["num_threads"]) for r in threadpool_info() if r["user_api"] == "blas"]
            self.assertEqual(before, after)
        self.assertTrue(observed)
        self.assertTrue(all(n == 4 for n in observed), observed)
        np.testing.assert_array_equal(weights, [2., 3.])


if __name__ == "__main__":
    unittest.main()
