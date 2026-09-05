import ast
import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np
from scipy.optimize import nnls

from skype_output_files import build_matrix_column_locations


SKYPE_ROOT = Path(__file__).resolve().parents[1]
RUN_NNLS_PATH = SKYPE_ROOT / "23_run_nnls.py"
SAVE_MATRIX_PATH = SKYPE_ROOT / "22_save_matrix.py"
BUILD_GRAPH_PATH = SKYPE_ROOT / "10_Graph_Find_Paths.py"


def load_functions(*names):
    tree = ast.parse(RUN_NNLS_PATH.read_text(encoding="utf-8"))
    functions = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in names
    ]
    if {node.name for node in functions} != set(names):
        raise AssertionError(f"Missing raw-NNLS helpers: {names}")
    calls = []

    def fake_bvls(matrix, target, lower, upper, n_threads):
        calls.append((matrix.copy(), target.copy(), n_threads))
        return SimpleNamespace(beta=nnls(matrix, target)[0].astype(matrix.dtype))

    namespace = {
        "MATRIX_CONTRACT": "depth_only_v1",
        "h5py": h5py,
        "np": np,
        "os": os,
        "bvls": fake_bvls,
    }
    module = ast.Module(body=functions, type_ignores=[])
    exec(compile(module, str(RUN_NNLS_PATH), "exec"), namespace)
    return namespace, calls


class RawNnlsTests(unittest.TestCase):
    def test_synthetic_matrix_matches_direct_full_column_nnls_once(self):
        namespace, calls = load_functions("fit_raw_nnls")
        matrix = np.asarray(
            [[1, 0, 1], [0, 1, 1], [1, 1, 0], [2, 0, 0]],
            dtype=np.float32,
        )
        target = np.asarray([2, 3, 4, 2], dtype=np.float32)

        observed = namespace["fit_raw_nnls"](matrix, target)
        expected = nnls(matrix, target)[0]

        np.testing.assert_allclose(observed, expected, rtol=1e-5, atol=1e-5)
        self.assertEqual(len(observed), matrix.shape[1])
        self.assertEqual(len(calls), 1)
        np.testing.assert_array_equal(calls[0][0], matrix)

    def test_explicitly_excluded_zero_column_stays_in_the_solve(self):
        namespace, calls = load_functions("fit_raw_nnls")
        matrix = np.asarray(
            [[1, 0], [1, 0], [0, 0]],
            dtype=np.float32,
        )
        target = np.asarray([2, 2, 0], dtype=np.float32)

        weights = namespace["fit_raw_nnls"](matrix, target)

        np.testing.assert_allclose(weights, [2, 0], atol=1e-6)
        self.assertEqual(calls[0][0].shape[1], 2)

    def test_main_contains_exactly_one_raw_solver_call(self):
        tree = ast.parse(RUN_NNLS_PATH.read_text(encoding="utf-8"))
        main = next(
            node
            for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name == "main"
        )
        calls = [
            node
            for node in ast.walk(main)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "fit_raw_nnls"
        ]
        self.assertEqual(len(calls), 1)

    def test_depth_only_contract_accepts_all_feature_rows(self):
        namespace, _ = load_functions("_contract_text", "load_depth_only_matrix")
        with tempfile.TemporaryDirectory() as prefix:
            with h5py.File(Path(prefix) / "matrix.h5", "w") as handle:
                handle.attrs["matrix_contract"] = "depth_only_v1"
                handle.create_dataset("A", data=np.ones((3, 4), dtype=np.float32))
                handle.create_dataset("A_fail", data=np.ones((3, 2), dtype=np.float32))
                handle.create_dataset("B", data=np.ones(4, dtype=np.float32))
                handle.create_dataset("B_fail", data=np.ones(2, dtype=np.float32))

            feature_depth, feature_fail, target_depth, target_fail = namespace[
                "load_depth_only_matrix"
            ](prefix)

        self.assertEqual(feature_depth.shape, (3, 4))
        self.assertEqual(feature_fail.shape, (3, 2))
        self.assertEqual(target_depth.shape, (4,))
        self.assertEqual(target_fail.shape, (2,))

    def test_old_prior_matrix_is_rejected(self):
        namespace, _ = load_functions("_contract_text", "load_depth_only_matrix")
        with tempfile.TemporaryDirectory() as prefix:
            with h5py.File(Path(prefix) / "matrix.h5", "w") as handle:
                handle.create_dataset("A", data=np.ones((2, 3), dtype=np.float32))
                handle.create_dataset("A_fail", data=np.ones((2, 1), dtype=np.float32))
                handle.create_dataset("B", data=np.ones(3, dtype=np.float32))
                handle.create_dataset("B_fail", data=np.ones(1, dtype=np.float32))
            with self.assertRaisesRegex(ValueError, "depth-only"):
                namespace["load_depth_only_matrix"](prefix)

    def test_stage22_has_no_prior_solve_or_prior_artifact(self):
        source = SAVE_MATRIX_PATH.read_text(encoding="utf-8")
        self.assertNotIn("normal_prior", source)
        self.assertNotIn("scipy.optimize import nnls", source)
        self.assertNotIn("normal_prior_data.pkl", source)
        self.assertIn("matrix_contract", source)
        self.assertIn(
            "set(range(n)) - explicit_excluded_columns",
            source,
        )
        self.assertIn('f"{PREFIX}/tot_loc_list.pkl"', source)

    def test_matrix_column_locations_follow_stage22_feature_order(self):
        observed = build_matrix_column_locations(
            "/result",
            [("/paths/1.index.txt", [1]), ("/paths/2.index.txt", [2])],
            [
                (1, "/out/front/1.win.stat.gz", "/out/front/1_base.win.stat.gz", -1, None),
            ],
            [
                (1, "/out/back/1.win.stat.gz", "/out/back/1_base.win.stat.gz", -1, None),
            ],
            [(1, "/out/ecdna/1.win.stat.gz")],
            [("chr1", {"dir": False}), ("chr2", {"dir": True})],
        )

        self.assertEqual(
            observed,
            [
                "/paths/1.index.txt",
                "/paths/2.index.txt",
                "/out/front/1_base.paf",
                "/out/back/1_base.paf",
                "/out/ecdna/1.paf",
                "/result/12_cent_fragment/chr1/left.fragment",
                "/result/12_cent_fragment/chr2/right.fragment",
            ],
        )

    def test_mode_cli_and_postprocessing_stages_are_removed(self):
        source = BUILD_GRAPH_PATH.read_text(encoding="utf-8")
        self.assertNotIn('"--karyotype_mode"', source)
        self.assertNotIn('"--variant_mode"', source)
        self.assertFalse((SKYPE_ROOT / "24_cluster_weight.py").exists())
        self.assertFalse((SKYPE_ROOT / "30_virtual_sky.py").exists())


if __name__ == "__main__":
    unittest.main()
