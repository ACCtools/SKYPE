import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import depth_scale  # noqa: E402


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


if __name__ == '__main__':
    unittest.main()
