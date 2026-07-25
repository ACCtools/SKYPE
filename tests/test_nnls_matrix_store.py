import ast
import unittest
from pathlib import Path

import numpy as np


RUN_NNLS_PATH = Path(__file__).resolve().parents[1] / "23_run_nnls.py"


def load_solver_matrix_from_store():
    tree = ast.parse(RUN_NNLS_PATH.read_text(encoding="utf-8"))
    function = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name == "solver_matrix_from_store"
    )
    module = ast.Module(body=[function], type_ignores=[])
    namespace = {}
    exec(compile(module, str(RUN_NNLS_PATH), "exec"), namespace)
    return namespace["solver_matrix_from_store"]


solver_matrix_from_store = load_solver_matrix_from_store()


class SolverMatrixStoreTests(unittest.TestCase):
    def test_compacted_workspace_prefix_is_an_exact_zero_copy_view(self):
        feature_major = np.arange(8 * 4, dtype=np.float32).reshape(8, 4)

        solver_matrix = solver_matrix_from_store(
            feature_major,
            n_features=5,
        )

        np.testing.assert_array_equal(
            solver_matrix,
            feature_major[:5, :].T,
        )
        self.assertTrue(np.shares_memory(solver_matrix, feature_major))

    def test_initial_fit_uses_direct_rescue_prefix_slice(self):
        source = RUN_NNLS_PATH.read_text(encoding="utf-8")
        self.assertIn(
            "A_initial = A_store[num_rescue_paths:, :].T",
            source,
        )
        self.assertNotIn("initial_A_idx_list", source)

    def test_only_added_pairs_use_info_for_censat_pair_logs(self):
        tree = ast.parse(RUN_NNLS_PATH.read_text(encoding="utf-8"))
        rescue_logs = []
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Call)
                and isinstance(node.func, ast.Attribute)
                and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "logging"
                and node.args
                and isinstance(node.args[0], ast.Constant)
                and isinstance(node.args[0].value, str)
                and "censat-pair" in node.args[0].value.lower()
            ):
                rescue_logs.append((node.func.attr, node.args[0].value))

        info_logs = [
            message
            for level, message in rescue_logs
            if level == "info"
        ]
        self.assertEqual(
            info_logs,
            ["Censat-pair rescue added: %s"],
        )
        self.assertTrue(
            all(
                level == "debug"
                or message == "Censat-pair rescue added: %s"
                for level, message in rescue_logs
            )
        )


if __name__ == "__main__":
    unittest.main()
