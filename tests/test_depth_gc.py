import gzip
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import depth_gc  # noqa: E402

CHROMS = [f'chr{i}' for i in range(1, 23)] + ['chrX']


def genome(rows_per_chrom=120, slope=-1.5, seed=0):
    """Rows of 100 kb on 23 chromosomes; two copy levels per chromosome; log2 depth carries a
    linear GC response of the given slope (per unit GC)."""
    rng = np.random.default_rng(seed)
    frames = []
    for k, chrom in enumerate(CHROMS):
        st = 1 + 100_000 * np.arange(rows_per_chrom)
        gc = np.clip(0.42 + 0.05 * np.sin(np.arange(rows_per_chrom) / 3 + k) + rng.normal(0, 0.03, rows_per_chrom),
                     0.30, 0.60)
        level = np.where(np.arange(rows_per_chrom) < rows_per_chrom // 2, 2.0, 3.0 if k % 2 else 2.0)
        depth = 10 * level * 2 ** (slope * (gc - 0.42)) * np.exp(rng.normal(0, 0.03, rows_per_chrom))
        frames.append(pd.DataFrame({'chr': chrom, 'st': st, 'nd': st + 99_999, 'meandepth': depth,
                                    'gc': gc, 'level': level}))
    return pd.concat(frames, ignore_index=True)


def write_gc_table(path, df):
    with gzip.open(path, 'wt') as fh:
        fh.write('#chr\tst\tnd\tgc\tnfrac\n')
        for r in df.itertuples(index=False):
            fh.write(f'{r.chr}\t{r.st}\t{r.nd}\t{r.gc!r}\t0.0\n')


def gc_slope(depth, df):
    resid = np.log2(depth / (10 * df.level.to_numpy()))
    return np.polyfit(df.gc.to_numpy(), resid, 1)[0]


class GcCorrectionTests(unittest.TestCase):
    def test_removes_the_gc_response_and_keeps_the_median(self):
        df = genome()
        out, train, info = depth_gc.correct(df, df.gc.to_numpy(), {})
        self.assertEqual(info['status'], 'corrected', info.get('reason'))
        self.assertGreater(info['lam'], 0.8)
        self.assertLess(abs(gc_slope(out, df)), 0.15 * abs(gc_slope(df.meandepth.to_numpy(), df)))
        self.assertAlmostEqual(np.median(out), np.median(df.meandepth), places=6)
        self.assertEqual(int(train.sum()), info['n_train'])

    def test_exempt_rows_and_unknown_gc_stay_unchanged(self):
        df = genome()
        gc = df.gc.to_numpy().copy()
        gc[:5] = np.nan
        chry = df.chr == 'chrX'
        df.loc[chry, 'chr'] = 'chrY'
        out, _, info = depth_gc.correct(df, gc, {})
        self.assertEqual(info['status'], 'corrected', info.get('reason'))
        np.testing.assert_array_equal(out[:5], df.meandepth.to_numpy()[:5])
        np.testing.assert_array_equal(out[chry.to_numpy()], df.meandepth.to_numpy()[chry.to_numpy()])

    def test_few_training_rows_fall_back_to_identity(self):
        df = genome(rows_per_chrom=40)
        out, _, info = depth_gc.correct(df, df.gc.to_numpy(), {})
        self.assertEqual(info['status'], 'identity')
        self.assertIn('training rows', info['reason'])
        np.testing.assert_array_equal(out, df.meandepth.to_numpy())

    def test_gc_table_of_another_reference_is_rejected(self):
        df = genome(rows_per_chrom=30)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, 'gc.tsv.gz')
            write_gc_table(path, df)
            gc, reason = depth_gc.row_gc(path, df)
            self.assertEqual(reason, '')
            np.testing.assert_allclose(gc, df.gc.to_numpy())
            other = df.copy()
            other.loc[other.chr == 'chr1', 'nd'] += 1234
            gc, reason = depth_gc.row_gc(path, other)
            self.assertIsNone(gc)
            self.assertIn('another reference', reason)


class Stage00CliTests(unittest.TestCase):
    def test_writes_gc_corrected_depth_and_can_skip_it(self):
        df = genome()
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            cols = ['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth']

            def stat(path, depth):
                t = df.assign(length=100_000, covsite=100_000, cov=100.0, meandepth=depth,
                              totaldepth=np.round(depth * 100_000).astype(int))
                with gzip.open(path, 'wt') as fh:
                    fh.write('#Chr\tStart\tEnd\tLength\tCoveredSite\tTotalDepth\tCoverage(%)\tMeanDepth\n')
                    t[cols].to_csv(fh, sep='\t', header=False, index=False)

            def run(*extra):
                subprocess.run([sys.executable, str(ROOT / '00_depth_norm.py'), str(tmp / 'sample.win.stat.gz'),
                                str(tmp / 'ref.win.stat.gz'), str(tmp / 'censat.bed'),
                                '--gc-table', str(tmp / 'gc.tsv.gz'), *extra], check=True, capture_output=True)
                return pd.read_csv(tmp / 'sample_normalized.win.stat.gz', sep='\t', header=None, names=cols)

            stat(tmp / 'sample.win.stat.gz', df.meandepth.to_numpy())
            stat(tmp / 'ref.win.stat.gz', np.full(len(df), 30.0))
            (tmp / 'censat.bed').write_text('chr1\t0\t1000\tcensat\n')
            write_gc_table(str(tmp / 'gc.tsv.gz'), df)
            out = run()
            off = run('--no-gc-correction')
        self.assertAlmostEqual(np.median(out.meandepth), np.median(df.meandepth), places=6)
        self.assertLess(abs(gc_slope(out.meandepth.to_numpy(), df)), 0.15 * abs(gc_slope(df.meandepth.to_numpy(), df)))
        np.testing.assert_array_equal(out.totaldepth, np.round(out.meandepth * 100_000).astype(int))
        np.testing.assert_allclose(off.meandepth, df.meandepth, rtol=1e-12)


if __name__ == '__main__':
    unittest.main()
