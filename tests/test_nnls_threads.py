import contextlib
import importlib.util
import io
import json
import pickle
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import h5py
import numpy as np
from scipy.optimize import nnls
from threadpoolctl import threadpool_info, threadpool_limits

from depth_nnls import predict_depth_from_weights


SPEC = importlib.util.spec_from_file_location(
    "run_nnls_threads_test", Path(__file__).resolve().parents[1] / "23_run_nnls.py"
)
stage = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stage)


def blas_limits():
    return {pool["filepath"]: pool["num_threads"] for pool in threadpool_info()
            if pool["user_api"] == "blas"}


class NnlsThreadTests(unittest.TestCase):
    def test_thread_cli_defaults_aliases_and_invalid_counts(self):
        parser = stage.build_parser()
        self.assertEqual(parser.parse_args(["result"]).thread, 1)
        for option in ("-t", "--thread", "--threads"):
            self.assertEqual(parser.parse_args(["result", option, "7"]).thread, 7)
        for value in ("0", "-1", "1.5", "abc"):
            with self.subTest(value=value), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit):
                    parser.parse_args(["result", "-t", value])

    def test_stage_applies_requested_threads_to_fit_and_predictions_then_restores(self):
        self._check_stage_threads(adelie_fallback=False)

    def test_stage_records_adelie_fallback_and_writes_predictions(self):
        self._check_stage_threads(adelie_fallback=True)

    def _check_stage_threads(self, *, adelie_fallback):
        seen_fit, seen_predictions, seen_adelie = [], [], []
        if adelie_fallback:
            from adelie.solver import bvls

        def capture_adelie(*args, **kwargs):
            seen_adelie.append(blas_limits())
            self.assertEqual(kwargs["n_threads"], 1)
            return bvls(*args, **kwargs)

        def capture_fit(*args, **kwargs):
            seen_fit.append(blas_limits())
            if adelie_fallback:
                raise RuntimeError("forced working-set failure")
            return nnls(*args, **kwargs)

        def capture_prediction(*args, **kwargs):
            seen_predictions.append(blas_limits())
            return predict_depth_from_weights(*args, **kwargs)

        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary)
            matrix = np.array([[1., 0.], [0., 1.], [1., 1.]], dtype=np.float32)
            target = np.array([2., 3., 5.], dtype=np.float32)
            with h5py.File(prefix / "matrix.h5", "w") as handle:
                handle.attrs["matrix_contract"] = "depth_only_v1"
                handle["A"] = matrix.T
                handle["B"] = target
                handle["A_fail"] = np.array([[1.], [2.]], dtype=np.float32)
                handle["B_fail"] = np.array([8.], dtype=np.float32)
            (prefix / "23_input.pkl").write_bytes(pickle.dumps(dict(
                matrix_contract="depth_only_v1", B_depth_start=0, B_depth_end=3,
                chr_filt_st_list=[("chr1", i * 100000) for i in range(3)],
            )))
            status = {"stages": {"initial": {"active_columns": [0, 1]}}}
            with contextlib.ExitStack() as stack, threadpool_limits(limits=3, user_api="blas"):
                if adelie_fallback:
                    stack.enter_context(patch("depth_nnls.FULL_SOLVE_COLUMN_THRESHOLD", 0))
                    stack.enter_context(patch("adelie.solver.bvls", side_effect=capture_adelie))
                before = blas_limits()
                with patch("depth_nnls.nnls", side_effect=capture_fit), \
                     patch.object(stage, "predict_depth_from_weights", side_effect=capture_prediction), \
                     patch.object(stage, "load_path_usage", return_value=[{}, {}]), \
                     patch.object(stage, "load_filter_status", return_value=status), \
                     patch.object(stage, "record_filter_stage"), \
                     patch.object(stage, "save_filter_status"):
                    stage.main([str(prefix), "--threads", "4"])
                self.assertEqual(blas_limits(), before)
            self.assertEqual(len(seen_fit), 1)
            self.assertEqual(len(seen_predictions), 2)
            self.assertEqual(len(seen_adelie), int(adelie_fallback))
            for limits in seen_fit + seen_predictions + seen_adelie:
                self.assertTrue(limits)
                self.assertEqual(set(limits.values()), {4})
            diagnostics = json.loads((prefix / "nnls_diagnostics.json").read_text())
            self.assertEqual(diagnostics["blas_threads_requested"], 4)
            self.assertEqual({p["num_threads"] for p in diagnostics["blas_threadpools"]}, {4})
            self.assertEqual(diagnostics["convergence_check"], "solver_default")
            self.assertIsNone(diagnostics["kkt_scaled_max"])
            self.assertIsNone(diagnostics["kkt_tolerance"])
            self.assertIsNone(diagnostics["kkt_matrix"])
            self.assertIsNone(diagnostics["converged"])
            if adelie_fallback:
                self.assertEqual(diagnostics["solver"], "adelie.solver.bvls")
                self.assertEqual(diagnostics["strategy"], "adelie_fallback")
                self.assertEqual(diagnostics["working_set_fallback"], "subproblem exceeded its iteration limit")
                self.assertEqual(diagnostics["solver_tolerance"], 1e-7)
            else:
                self.assertEqual(diagnostics["solver"], "scipy.optimize.nnls")
                self.assertEqual(diagnostics["strategy"], "full")
                self.assertIsNone(diagnostics["solver_tolerance"])
            fit = np.load(prefix / "weight.npy")
            prediction = np.load(prefix / "predict_B.npy")
            expected = np.r_[matrix @ fit, fit[0] + 2 * fit[1]]
            np.testing.assert_allclose(prediction, expected, atol=1e-12)
            np.testing.assert_allclose(prediction, [2., 3., 5., 8.], atol=5e-3 if adelie_fallback else 1e-12)


if __name__ == "__main__":
    unittest.main()
