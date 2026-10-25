import copy
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import circos_plotting  # noqa: E402
import depth_scale  # noqa: E402
from structure_nclose import StructureWeights  # noqa: E402
from test_stage_nclose_output import alias_fixture, export, read_tsv  # noqa: E402


DEPTH = np.array([10.0, 20.0, 30.0])


def write_table(path, depth, chrom='chr1'):
    with gzip.open(path, 'wt') as fh:
        for i, d in enumerate(depth):
            fh.write(f'{chrom}\t{1 + 100_000 * i}\t{100_000 * (i + 1)}\t100000\t100000\t0\t100\t{d!r}\n')
        fh.write('chrM\t1\t16569\t16569\t16569\t0\t100\t999.0\n')


class CopyUnitTests(unittest.TestCase):
    def write(self, prefix, n, depth=DEPTH, depth_file='depth.win.stat.gz'):
        return depth_scale.write(prefix, depth_file, np.array(depth),
                                 {'status': 'corrected'}, {'N': n, 'status': 'supported'})

    def test_baf_unit_is_used_for_the_same_depth_table(self):
        with tempfile.TemporaryDirectory() as prefix:
            scale = self.write(prefix, 7.0)
            self.assertEqual((scale['copy_unit'], scale['copy_unit_source']), (7.0, 'baf'))
            self.assertEqual(depth_scale.copy_unit(prefix, DEPTH[::-1].copy()), (7.0, 'baf'))

    def test_falls_back_to_median_half(self):
        with tempfile.TemporaryDirectory() as prefix:
            self.assertEqual(depth_scale.copy_unit(prefix, DEPTH), (10.0, 'median'))  # no json
            self.write(prefix, None)
            self.assertEqual(depth_scale.copy_unit(prefix, DEPTH), (10.0, 'median'))  # no BAF N
            self.write(prefix, 7.0)
            for other in ([12.0, 24.0, 36.0], [10.0, 20.0, 31.0], [10.0, 20.0, 30.0, 20.0]):
                with self.assertLogs(level='WARNING'):  # another table, even with the same median
                    self.assertEqual(depth_scale.copy_unit(prefix, np.array(other)),
                                     (np.median(other) / 2, 'median'))
            (Path(prefix) / depth_scale.JSON_NAME).write_text('{')
            with self.assertLogs(level='WARNING'):
                self.assertEqual(depth_scale.copy_unit(prefix, DEPTH), (10.0, 'median'))  # unreadable

    def test_viewer_unit_checks_the_recorded_table_and_the_shown_depths(self):
        coords = [('chr1', 1), ('chr1', 100_001), ('chr1', 200_001)]
        shown = DEPTH[[2, 0, 1]]
        order = [coords[i] for i in (2, 0, 1)]   # stage 22 orders bins its own way
        with tempfile.TemporaryDirectory() as prefix:
            table = Path(prefix) / 'sample_normalized.win.stat.gz'
            self.assertIsNone(depth_scale.recorded_copy_unit(prefix, order, shown))  # no json
            write_table(table, DEPTH)
            self.assertEqual(list(depth_scale.read_table(table).items()), list(zip(coords, DEPTH)))
            self.write(prefix, 7.0, depth_file=str(table))
            self.assertEqual(depth_scale.recorded_copy_unit(prefix, order, shown), 7.0)
            self.assertEqual(depth_scale.recorded_copy_unit(prefix, order[:2], shown[:2]), 7.0)
            # A float32 midpoint: Python reads 20.000000953674316 as 20.0 in float32, pandas
            # (stage 22) one float32 step above.
            mid = np.array([10.0, 20.000000953674316, 30.0])
            write_table(table, mid)
            self.write(prefix, 7.0, mid, depth_file=str(table))
            stage22 = np.float32([30.0, 10.0, 20.000001907348633])
            self.assertEqual(depth_scale.recorded_copy_unit(prefix, order, stage22), 7.0)
            write_table(table, DEPTH)
            self.write(prefix, 7.0, depth_file=str(table))
            for other in ([30.0, 10.0, 21.0], [20.0, 10.0, 30.0]):   # not the table's depths
                with self.assertLogs(level='WARNING'):
                    self.assertIsNone(depth_scale.recorded_copy_unit(prefix, order, np.array(other)))
            with self.assertLogs(level='WARNING'):
                self.assertIsNone(depth_scale.recorded_copy_unit(prefix, [('chr2', 1)], [10.0]))
            write_table(table, [10.0, 20.0, 31.0])
            with self.assertLogs(level='WARNING'):
                self.assertIsNone(depth_scale.recorded_copy_unit(prefix, order, shown))
            table.unlink()
            with self.assertLogs(level='WARNING'):
                self.assertIsNone(depth_scale.recorded_copy_unit(prefix, order, shown))
            self.write(prefix, None, depth_file=str(table))
            write_table(table, DEPTH)
            self.assertIsNone(depth_scale.recorded_copy_unit(prefix, order, shown))

    def test_remove(self):
        with tempfile.TemporaryDirectory() as prefix:
            depth_scale.remove(prefix)
            self.write(prefix, 7.0)
            depth_scale.remove(prefix)
            self.assertFalse((Path(prefix) / depth_scale.JSON_NAME).exists())

    def test_json_has_no_non_finite_values(self):
        with tempfile.TemporaryDirectory() as prefix:
            depth_scale.write(prefix, 'd', np.array([1.0, 3.0]), {'x': np.nan},
                              {'N': None, 'status': 'insufficient_data', 'v': np.float64(np.inf)})
            scale = json.loads((Path(prefix) / depth_scale.JSON_NAME).read_text())
        self.assertIsNone(scale['gc_correction']['x'])
        self.assertIsNone(scale['baf']['v'])
        self.assertEqual(scale['copy_unit'], 1.0)
        self.assertEqual((scale['depth_rows'], scale['depth_sum']), (2, 4.0))


class ThresholdUnitTests(unittest.TestCase):
    """Reported copy numbers follow n_unit; the min_cn cuts follow threshold_unit."""

    def contexts(self, weight, unit=8.):
        _, model = alias_fixture()
        return (StructureWeights(copy.deepcopy(model), [weight, weight], 20., threshold_unit=20.),
                StructureWeights(copy.deepcopy(model), [weight, weight], unit, threshold_unit=20.))

    def test_cut_does_not_move_with_the_reporting_unit(self):
        # Two structures of weight 3 give the NClose 6 = 0.3 threshold-unit copies: shown.
        same, scaled = self.contexts(3.)
        self.assertEqual(len(same.visible_ncloses()), len(scaled.visible_ncloses()))
        self.assertGreater(len(scaled.visible_ncloses()), 0)
        totals = [t for _, t in scaled.visible_ncloses()]
        inv, trans, indel = scaled.cn_lists()
        self.assertEqual(sorted(inv + trans + indel), sorted(t / 8. for t in totals))
        # Total 1.5 is 0.075 threshold-unit copies (0.19 in the reporting unit 8): hidden.
        same, scaled = self.contexts(0.75)
        self.assertEqual(same.visible_ncloses(), [])
        self.assertEqual(scaled.visible_ncloses(), [])
        self.assertEqual(scaled.cn_lists(), ([], [], []))

    def test_reports_use_the_reporting_unit(self):
        nodes = alias_fixture()[0]
        _, scaled = self.contexts(3.)
        with tempfile.TemporaryDirectory() as prefix:
            records = export(Path(prefix), scaled, nodes)
            bed = read_tsv(Path(prefix) / 'SKYPE_result.bed')
        self.assertTrue(records)
        self.assertTrue(bed)
        totals = {round(t / 8., 9) for _, t in scaled.visible_ncloses()}
        self.assertTrue(all(round(float(r['weight (N)']), 9) in totals for r in bed))
        _, hidden = self.contexts(0.75)
        with tempfile.TemporaryDirectory() as prefix:
            self.assertEqual(export(Path(prefix), hidden, nodes), [])
            self.assertEqual(read_tsv(Path(prefix) / 'SKYPE_result.bed'), [])

    def test_weight_at_the_cut_is_hidden_in_every_report(self):
        # Total 2 is exactly 0.1 threshold-unit copies. In the reporting unit 3.9 the rounded
        # copy number re-expressed in the threshold unit would land just above 0.1.
        nodes = alias_fixture()[0]
        _, scaled = self.contexts(1., unit=3.9)
        self.assertEqual(scaled.visible_ncloses(), [])
        with tempfile.TemporaryDirectory() as prefix:
            self.assertEqual(export(Path(prefix), scaled, nodes), [])
            self.assertEqual(read_tsv(Path(prefix) / 'SKYPE_result.bed'), [])


class CircosUnitTests(unittest.TestCase):
    def test_levels_above_the_depth_axis_are_skipped(self):
        # 4N = 64 lies above the 3 x median = 60 axis; pycirclize rejects such a line.
        st = 1 + 100_000 * np.arange(10)
        depth = pd.DataFrame({'chr': 'chr1', 'st': st, 'meandepth': 20.0})
        drawn = []
        line = circos_plotting.Track.line

        def record(track, x, y, **kwargs):
            if kwargs.get('color') in ('blue', 'red'):
                drawn.append((y[0], kwargs['color']))
            return line(track, x, y, **kwargs)

        with tempfile.TemporaryDirectory() as tmp, \
                mock.patch.object(circos_plotting.Track, 'line', record), \
                mock.patch.object(circos_plotting.Circos, 'plotfig', side_effect=RuntimeError('drawn')):
            (Path(tmp) / 'cyto.bed').write_text('chr1\t0\t1000000\tp1\tgneg\n')
            with self.assertRaisesRegex(RuntimeError, 'drawn'):
                circos_plotting.render_total_coverage_circos(
                    output_prefix=tmp, fig_prefix='', chromosome_lengths={'chr1': 1_000_000},
                    cytobands_path=str(Path(tmp) / 'cyto.bed'), observed_depth=depth,
                    prediction_coords=[('chr1', int(s)) for s in st], predicted_depth=np.full(10, 20.0),
                    prediction_labels=np.zeros(10), mean_depth=20.0, copy_unit=16.0)
        self.assertEqual(drawn, [(16.0, 'blue')])


if __name__ == '__main__':
    unittest.main()
