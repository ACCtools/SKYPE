import ast
import gzip
import logging
import pickle
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from nclose_tracking import _build_bnd_event
from skype_utils import require_fixed_depth_matrix


ROOT = Path(__file__).resolve().parents[1]


def write_depth_fixture(prefix):
    """One path, a cutoff-boundary bin, a high bin, a CenSat bin and chrM."""
    coordinates = [("chr1", i * 100000 + 1) for i in range(6)]
    coordinates += [("chrY", 1), ("chrM", 1)]
    observed = [10., 20., 20., 60., 300., 40., 20., 10000.]
    feature = [1., 1., 1., 1., 1., 1., 0., 1.]
    depth_dir = prefix / "21_pat_depth"
    depth_dir.mkdir()
    for path, values in ((prefix / "observed.win.stat.gz", observed),
                         (depth_dir / "1.win.stat.gz", feature)):
        with gzip.open(path, "wt") as handle:
            for (chrom, start), value in zip(coordinates, values):
                fields = (chrom, start, start + 99999, 100000, 100000,
                          value * 100000, 100, value)
                handle.write("\t".join(map(str, fields)) + "\n")
    (prefix / "censat.bed").write_text("chr1\t500000\t600000\nchrY\t500000\t600000\n")
    path = [("chr1f", 0, 0), (2, 0, 0, 0), (3, 1, 0, 0), ("chr1b", 0, 0)]
    path_file = str(prefix / "20_depth" / "chr1f_chr1b" / "1.paf")
    nodes = [["u1", 200000, 0, 100000, "+", "chr1", 600000, 0, 100000],
             ["u1", 200000, 100000, 200000, "+", "chr1", 600000, 500000, 600000]]
    artifacts = {
        "nclose_nodes.pkl": {"u1": [(0, 1)]},
        "nclose_event_catalog.pkl": [_build_bnd_event((0, 1), nodes, "compressed_nclose_nodes_list.txt")],
        "path_data.pkl": {"chr1f_chr1b": [(path, [("chr1", 600000)])]},
        "contig_pat_vec_data.pkl": ([(path_file, [1])], [1], {}, [0.]),
        "ecdna_circuit_data.pkl": ([], {}),
        "cen_fragment_data.pkl": {},
        "indel_exclude_idx_set.pkl": set(),
        "conjoined_type4_ins_del.pkl": ([], []),
        "01_nclose_data.pkl": {"contig_data": nodes},
    }
    for name, value in artifacts.items():
        (prefix / name).write_bytes(pickle.dumps(value))


class FixedDepthMatrixTests(unittest.TestCase):
    def run_stage(self, name, *args):
        result = subprocess.run([sys.executable, str(ROOT / name), *map(str, args)],
                                capture_output=True, text=True, timeout=60)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_matrix_and_fit_use_original_depth_exclusions(self):
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary)
            write_depth_fixture(prefix)
            stale = [prefix / "high_depth_gate.tsv", prefix / "unexplained_high_depth.npy"]
            for path in stale:
                path.write_bytes(b"old diagnostic")
            self.run_stage("22_save_matrix.py", prefix / "censat.bed",
                           prefix / "observed.win.stat.gz", prefix, "-t", 1)
            metadata = pickle.loads((prefix / "23_input.pkl").read_bytes())
            expected_coordinates = [("chr1", i * 100000 + 1) for i in range(4)] + [("chrY", 1)]
            self.assertEqual(metadata, dict(matrix_contract="depth_only_v1",
                                           chr_filt_st_list=expected_coordinates,
                                           B_depth_start=0, B_depth_end=5))
            with h5py.File(prefix / "matrix.h5") as matrix:
                np.testing.assert_array_equal(matrix["A"][:], [[1., 1., 1., 1., 0.]])
                np.testing.assert_array_equal(matrix["B"][:], [10., 20., 20., 60., 20.])
                np.testing.assert_array_equal(matrix["A_fail"][:], [[1., 1.]])
                np.testing.assert_array_equal(matrix["B_fail"][:], [300., 40.])
            self.assertTrue(all(not path.exists() for path in stale))
            # A cached coefficient must not affect the ordinary NNLS solution.
            np.save(prefix / "weight.npy", [999.])
            self.run_stage("23_run_nnls.py", prefix, "-t", 1)
            np.testing.assert_allclose(np.load(prefix / "weight.npy"), [27.5], atol=1e-12)
            np.testing.assert_allclose(np.load(prefix / "predict_B.npy"),
                                       [27.5, 27.5, 27.5, 27.5, 0., 27.5, 27.5], atol=1e-12)
            self.assertTrue(all(not path.exists() for path in stale))

    def test_old_fixed_depth_metadata_is_still_accepted(self):
        for metadata in ({}, {"depth_policy": "legacy"}):
            with self.subTest(metadata=metadata):
                require_fixed_depth_matrix(metadata)

    def test_refit_rejects_retired_matrix_before_touching_saved_weights(self):
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary)
            np.save(prefix / "weight.npy", [999.])
            for policy in ("all", "nclose_l2", "nclose_huber"):
                with self.subTest(policy=policy):
                    (prefix / "23_input.pkl").write_bytes(pickle.dumps(dict(
                        matrix_contract="depth_only_v1", depth_policy=policy)))
                    result = subprocess.run([sys.executable, str(ROOT / "23_run_nnls.py"), str(prefix)],
                                            capture_output=True, text=True, timeout=60)
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn("rerun SKYPE from stage 22", result.stderr)
                    np.testing.assert_array_equal(np.load(prefix / "weight.npy"), [999.])


class Stage31CoordinateTests(unittest.TestCase):
    def test_output_uses_saved_fitted_and_excluded_depth_row_order(self):
        source = ROOT / "31_depth_analysis.py"
        text = source.read_text()
        tree = ast.parse(text)
        begin = next(i for i, node in enumerate(tree.body) if isinstance(node, ast.With)
                     and "23_input.pkl" in ast.get_source_segment(text, node))
        end = next(i for i in range(begin, len(tree.body)) if isinstance(tree.body[i], ast.FunctionDef))
        block = ast.Module(body=tree.body[begin:end], type_ignores=[])
        frame = pd.DataFrame(dict(chr=["chr1", "chr1", "chr2", "chr2"], st=[1, 100001, 1, 100001]))
        clean = [("chr1", 1), ("chr2", 1), ("chr2", 100001)]
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary)
            (prefix / "23_input.pkl").write_bytes(pickle.dumps(dict(
                matrix_contract="depth_only_v1", B_depth_start=0, B_depth_end=3,
                chr_filt_st_list=clean)))
            np.save(prefix / "predict_B.npy", np.arange(4.))
            namespace = dict(PREFIX=str(prefix), pkl=pickle, np=np, pd=pd, logging=logging,
                             df=frame, B=np.arange(4.), DEPTH_ONLY_MATRIX_CONTRACT="depth_only_v1")
            exec(compile(block, str(source), "exec"), namespace)
            self.assertEqual(namespace["chr_filt_st_list"], clean)
            self.assertEqual(namespace["chr_no_filt_st_list"], [("chr1", 100001)])
            self.assertEqual(namespace["chr_filt_idx_list"], [0, 2, 3])
            saved = pd.read_csv(prefix / "stage31_depth_coordinates.tsv", sep="\t")
            self.assertEqual(list(zip(saved.chrom, saved.start)), clean + [("chr1", 100001)])


if __name__ == "__main__":
    unittest.main()
