import ast
import importlib.util
import json
import logging
import pickle
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import h5py
import numpy as np
import pandas as pd
from scipy.optimize import nnls
from threadpoolctl import threadpool_limits

import high_depth


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('high_depth_stage23', ROOT / '23_run_nnls.py')
stage23 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stage23)


def depth_frame(values):
    return pd.DataFrame(dict(chr=['chr1'] * len(values),
                             st=[i*100000+1 for i in range(len(values))],
                             nd=[(i+1)*100000 for i in range(len(values))],
                             meandepth=values))


class HighDepthGateTests(unittest.TestCase):
    def check_empty_gate(self, frame, censat):
        clean = list(zip(frame.chr, frame.st))
        with patch('high_depth.nclose_intervals', side_effect=AssertionError('unneeded topology read')), \
             patch('high_depth.estimate_noise_sigma_by_chromosome', side_effect=AssertionError('unneeded noise estimate')):
            result = high_depth.high_depth_gate('missing-prefix', frame, clean, censat)
        self.assertTrue(result.empty)
        self.assertEqual(list(result.columns), list(high_depth.HIGH_DEPTH_COLUMNS))
        self.assertEqual(result.nclose_gate.dtype, np.dtype(bool))
        self.assertEqual(result.direct_nclose_gate.dtype, np.dtype(bool))
        self.assertEqual(list(result.loc[result.nclose_gate].columns), list(result.columns))
        self.assertEqual(list(result[result.nclose_gate].columns), list(result.columns))

    def test_no_high_depth_bins_has_typed_empty_schema(self):
        self.check_empty_gate(depth_frame([60., 61., 59.]), {})

    def test_censat_can_exclude_every_high_depth_candidate(self):
        self.check_empty_gate(depth_frame([60., 60., 60., 300.]), {'chr1': [(300000, 400000)]})

    def test_high_depth_without_nclose_selects_no_rows(self):
        frame = depth_frame([60., 60., 60., 300.])
        with patch('high_depth.nclose_intervals', return_value={}):
            result = high_depth.high_depth_gate('unused', frame, list(zip(frame.chr[:3], frame.st[:3])), {})
        self.assertEqual(len(result), 1)
        self.assertFalse(result.nclose_gate.any())
        self.assertTrue(result.loc[result.nclose_gate].empty)

    def test_nclose_presence_extends_only_through_contiguous_high_bins(self):
        frame = depth_frame([60., 60., 60., 60., 300., 310., 60., 300.])
        clean = list(zip(frame.loc[frame.meandepth <= 180, 'chr'], frame.loc[frame.meandepth <= 180, 'st']))
        with patch('high_depth.nclose_intervals', return_value={'chr1': [(400000, 420000, 'NClose0')]}):
            result = high_depth.high_depth_gate('unused', frame, clean, {})
        self.assertEqual(result.direct_nclose_gate.tolist(), [True, False, False])
        self.assertEqual(result.nclose_gate.tolist(), [True, True, False])
        self.assertEqual(result.loc[result.nclose_gate, 'start'].tolist(), [400001, 500001])


class RobustEmptyRowsTests(unittest.TestCase):
    def test_no_high_rows_calls_original_nnls_and_preserves_predictions(self):
        matrix = np.array([[1., -1., 0.], [0., 1., 0.], [1., 0., 0.], [2., 1., 0.]])
        target = np.array([1., 1., 2., 5.])
        expected, expected_info = high_depth.solve_depth_nnls(matrix, target, return_diagnostics=True)
        with patch('high_depth.minimize', side_effect=AssertionError('robust solver must be skipped')), \
             patch('high_depth.solve_depth_nnls', wraps=high_depth.solve_depth_nnls) as fallback:
            weights, info = high_depth.solve_robust(matrix, target, np.zeros(4, dtype=bool), np.empty(0))
        self.assertEqual(fallback.call_count, 1)
        np.testing.assert_array_equal(weights, expected)
        self.assertEqual(info['strategy'], expected_info['strategy'])
        self.assertEqual(info['high_depth_fallback'], 'no_high_depth_rows')
        self.assertEqual(info['unexplained_high_depth_sum'], 0)

    def test_nonempty_robust_loss_still_limits_an_unexplainable_spike(self):
        matrix = np.ones((101, 1)); target = np.r_[np.ones(100), 100.]
        high = np.arange(101) == 100
        with threadpool_limits(limits=1):
            weights, info = high_depth.solve_robust(matrix, target, high, np.array([3.]))
        self.assertAlmostEqual(weights[0], 1.03, places=7)
        self.assertLessEqual(info['kkt_scaled_max'], 1e-10)

    def test_stage23_empty_policy_does_not_need_cached_weights(self):
        matrix = np.array([[1., 0.], [0., 1.], [1., 1.]], dtype=np.float32)
        target = np.array([2., 3., 5.], dtype=np.float32)
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary)
            with h5py.File(prefix / 'matrix.h5', 'w') as f:
                f.attrs['matrix_contract'] = 'depth_only_v1'
                f['A'], f['B'] = matrix.T, target
                f['A_fail'], f['B_fail'] = np.array([[1.], [2.]], dtype=np.float32), np.array([8.], dtype=np.float32)
            (prefix / '23_input.pkl').write_bytes(pickle.dumps(dict(
                matrix_contract='depth_only_v1', B_depth_start=0, B_depth_end=3,
                chr_filt_st_list=[('chr1', i*100000+1) for i in range(3)],
                depth_policy='nclose_huber', high_depth_rows=[], high_depth_tau=[])))
            with patch.object(stage23, 'fit_raw_nnls', wraps=stage23.fit_raw_nnls) as raw, \
                 patch('high_depth.solve_robust', side_effect=AssertionError('no high rows')), \
                 patch.object(stage23, 'load_path_usage', return_value=[{}, {}]), \
                 patch.object(stage23, 'load_filter_status', return_value={'stages': {'initial': {'active_columns': [0, 1]}}}), \
                 patch.object(stage23, 'record_filter_stage'), patch.object(stage23, 'save_filter_status'):
                stage23.main([str(prefix)])
            self.assertEqual(raw.call_count, 1)
            np.testing.assert_array_equal(np.load(prefix / 'weight.npy'), nnls(matrix, target)[0])
            np.testing.assert_allclose(np.load(prefix / 'predict_B.npy'), [2., 3., 5., 8.])
            self.assertEqual(np.load(prefix / 'unexplained_high_depth.npy').shape, (0,))
            info = json.loads((prefix / 'nnls_diagnostics.json').read_text())
            self.assertEqual(info['high_depth_fallback'], 'no_high_depth_rows')


class Stage31CoordinateTests(unittest.TestCase):
    def test_output_uses_saved_depth_row_order_after_policy_changes(self):
        source = ROOT / '31_depth_analysis.py'
        tree = ast.parse(source.read_text())
        begin = next(i for i, node in enumerate(tree.body) if isinstance(node, ast.With)
                     and '23_input.pkl' in ast.get_source_segment(source.read_text(), node))
        end = next(i for i in range(begin, len(tree.body)) if isinstance(tree.body[i], ast.FunctionDef))
        block = ast.Module(body=tree.body[begin:end], type_ignores=[])
        frame = pd.DataFrame(dict(chr=['chr1','chr1','chr2','chr2'], st=[1,100001,1,100001]))
        clean = [('chr1',1), ('chr2',1), ('chr2',100001)]
        with tempfile.TemporaryDirectory() as tmp:
            prefix = Path(tmp)
            (prefix/'23_input.pkl').write_bytes(pickle.dumps(dict(matrix_contract='depth_only_v1',
                B_depth_start=0, B_depth_end=3, chr_filt_st_list=clean)))
            np.save(prefix/'predict_B.npy', np.arange(4.))
            namespace = dict(PREFIX=tmp, pkl=pickle, np=np, pd=pd, logging=logging,
                             df=frame, B=np.arange(4.), DEPTH_ONLY_MATRIX_CONTRACT='depth_only_v1')
            exec(compile(block, str(source), 'exec'), namespace)
            self.assertEqual(namespace['chr_filt_st_list'], clean)
            self.assertEqual(namespace['chr_no_filt_st_list'], [('chr1',100001)])
            self.assertEqual(namespace['chr_filt_idx_list'], [0,2,3])
            saved = pd.read_csv(prefix/'stage31_depth_coordinates.tsv', sep='\t')
            self.assertEqual(list(zip(saved.chrom, saved.start)), clean+[('chr1',100001)])


if __name__ == '__main__':
    unittest.main()
