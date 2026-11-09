import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import baf_ploidy as PB  # noqa: E402

N_TRUE = 14.0
HET_RATE, CALLABLE = 1e-3, 15000


def simulate(pieces, seed):
    """pieces: (rows, major, minor), four to a chromosome of 100 kb rows. Depth = N x copies
    with 1.5% per-piece and 2% per-row noise; germline het sites 1 per kb of 15 kb callable
    sequence, read depth ~ Poisson(row depth), minor reads binomial at minor / copies, kept by
    the het-site calling rule. LOH pieces keep non-germline sites at 3e-5 per bp."""
    rng = np.random.default_rng(seed)
    rows, plan, used = [], [], {}
    for k, (n_rows, a, b) in enumerate(pieces):
        chrom = f'chr{k // 4 + 1}'
        first = used.get(chrom, 0)
        used[chrom] = first + n_rows
        st = 1 + PB.BIN * (first + np.arange(n_rows))
        depth = N_TRUE * (a + b) * (1 + rng.normal(0, 0.015)) * (1 + rng.normal(0, 0.02, n_rows))
        rows.append(pd.DataFrame({'chr': chrom, 'st': st, 'nd': st + PB.BIN - 1, 'depth': depth}))
        plan.append((chrom, st, depth, min(a, b) / (a + b)))
    df = pd.concat(rows, ignore_index=True)
    cap = 3 * np.median(df.depth)
    sites, bins = [], []
    for chrom, st, depth, q in plan:
        het = np.zeros(len(st), dtype=int)
        for i in range(len(st)):
            m = rng.poisson(CALLABLE * (HET_RATE if q > 0 else 3e-5))
            frac = np.full(m, q) if q > 0 else rng.uniform(0.10, 0.50, m)
            n = rng.poisson(depth[i], m)
            k = rng.binomial(n, np.where(rng.random(m) < 0.5, frac, 1 - frac))
            mn = np.minimum(k, n - k)
            keep = (n >= 8) & (n <= cap) & (mn >= PB.MIN_MINOR) & (mn >= PB.MIN_FRAC * n)
            het[i] = keep.sum()
            sites.append(pd.DataFrame({'chrom': chrom, 'pos': st[i] + 40_000 + np.arange(keep.sum()),
                                       'depth': n[keep], 'minor': mn[keep]}))
        bins.append(pd.DataFrame({'chrom': chrom, 'start': st, 'callable': CALLABLE, 'het': het}))
    return df, pd.concat(sites, ignore_index=True), pd.concat(bins, ignore_index=True)


def run(pieces, seed=7):
    df, sites, bins = simulate(pieces, seed)
    d = df.depth.to_numpy()
    train = d <= 3 * np.median(d)
    return PB.estimate(df, d, train, np.ones(len(df), dtype=bool), sites, bins), np.median(d)


class BafPloidyTests(unittest.TestCase):
    def test_near_triploid_scale_is_recovered_where_the_median_is_three_copies(self):
        pieces = ([(400, 2, 1)] * 8 + [(300, 1, 1)] * 4 + [(250, 2, 2)] * 3 + [(200, 2, 0)] * 3
                  + [(150, 3, 1)] * 2 + [(150, 1, 0)])
        (N, segments, diag), median = run(pieces)
        self.assertEqual(diag['status'], 'supported', diag)
        self.assertAlmostEqual(N / N_TRUE, 1.0, delta=0.03)
        self.assertGreater(median / 2 / N_TRUE, 1.3)   # median / 2 would be 1.5 N
        self.assertAlmostEqual(diag['ploidy'], 2.9, delta=0.15)
        self.assertTrue({'balanced', 'third', 'loh'} <= set(segments.state))
        self.assertEqual(diag['fit_quality'], 'adequate')
        self.assertTrue(diag['scale_identified'], diag['scale_rivals'])
        self.assertIsNone(diag['scale_unidentified_by'])
        self.assertGreaterEqual(diag['scale_margin_fewer_copies'], PB.FIT_MARGIN)
        relations = {r['relation']: r for r in diag['scale_rivals']}
        self.assertTrue({'N/2', '2N', 'median/2'} <= set(relations))
        self.assertAlmostEqual(relations['N/2']['N'], N / 2)
        self.assertEqual([x['sigma_c'] for x in diag['sensitivity_wide']],
                         [PB.SIGMA_C * f for f in PB.SENS_WIDE])

    def test_near_diploid_matches_the_median(self):
        pieces = [(400, 1, 1)] * 10 + [(200, 2, 1)] * 3 + [(150, 1, 0)] * 3 + [(150, 2, 0)] * 2
        (N, _, diag), median = run(pieces)
        self.assertEqual(diag['status'], 'supported', diag)
        self.assertAlmostEqual(N / N_TRUE, 1.0, delta=0.03)
        self.assertAlmostEqual(median / 2 / N_TRUE, 1.0, delta=0.05)

    def test_identified_uses_the_selection_gaps(self):
        self.assertTrue(PB.identified(None, None))
        self.assertTrue(PB.identified(PB.FIT_MARGIN, -PB.AMBIG_GAP))
        self.assertTrue(PB.identified(0.5, 0.90 - 0.93))                # exactly the gap, after rounding
        self.assertFalse(PB.identified(PB.FIT_MARGIN - 0.01, 0.2))    # fewer copies nearly as good
        self.assertFalse(PB.identified(0.5, -PB.AMBIG_GAP - 0.01))     # more copies fit better

    def test_scale_rivals_sides_and_window(self):
        scored = [(10.0, 'N', 0.90, 0.9, 0.9), (5.0, 'N/2', 0.95, 0.9, 0.9), (20.0, '2N', 0.30, 0.3, 0.3),
                  (9.5, 'median/2', 0.99, 0.9, 0.9), (10.8, 'candidate', 0.99, 0.9, 0.9),
                  (14.0, 'candidate', 0.85, 0.8, 0.8)]
        rivals, fewer, more = PB.scale_rivals(10.0, scored, 30.0)
        self.assertEqual([(r['relation'], r['N']) for r in rivals],   # 9.5 and 10.8: same family as 10
                         [('N/2', 5.0), ('candidate', 14.0), ('2N', 20.0)])
        self.assertAlmostEqual(fewer, 0.90 - 0.85)
        self.assertAlmostEqual(more, 0.90 - 0.95)
        self.assertAlmostEqual(rivals[1]['median_copies'], 30.0 / 14.0)
        self.assertEqual(PB.scale_rivals(10.0, scored[:1], 30.0), ([], None, None))

    def test_refine_keeps_the_grid_point_when_the_median_raises_the_cost(self):
        # LOH at 2, 3, 4 and 5 copies of 10; the 5-copy segment sits 1% high. The weighted median
        # of D / C is exactly 10, a higher mean cost than the grid minimum just above 10.
        s = pd.DataFrame({'depth': [20.0, 30.0, 40.0, 50.5], 'bins': [400.0, 300.0, 200.0, 300.0],
                          'state': ['loh'] * 4})
        bc, cred = np.zeros((4, PB.CMAX)), np.ones((4, PB.CMAX), dtype=bool)
        N, _ = PB.fit_scale(s, bc, cred, 30.0)
        grid = np.exp(np.linspace(np.log(30.0 / PB.CMAX), np.log(30.0 * 1.2), 1500))

        def mean_cost(n):
            return float((PB.costs(s.depth.to_numpy(), bc, np.array([n]))[0][:, 0] * s.bins).sum() / s.bins.sum())

        self.assertIn(N, set(grid.tolist()))
        self.assertNotEqual(N, 10.0)
        self.assertLess(mean_cost(N), mean_cost(10.0))
        self.assertAlmostEqual(N, 10.0, delta=0.1)

    def test_failing_rival_is_on_the_failing_side(self):
        rivals = [{'N': 5.0, 'anchor_fit': 1.0}, {'N': 20.0, 'anchor_fit': 0.8}]
        self.assertIsNone(PB.failing_rival(rivals, 10.0, 0.5, 0.0))
        # Fewer-copy margin exactly at FIT_MARGIN after rounding passes; the more-copy side fails.
        self.assertEqual(PB.failing_rival(rivals, 10.0, 0.90 - 0.80, -0.10)['N'], 5.0)
        self.assertEqual(PB.failing_rival(rivals, 10.0, 0.05, 0.5)['N'], 20.0)

    def test_segments_that_fit_no_state_have_no_call(self):
        # 2.5 copies with minor fraction 0.2 fits no whole state at the near-diploid scale.
        pieces = ([(400, 1, 1)] * 10 + [(200, 2, 1)] * 3 + [(150, 1, 0)] * 3
                  + [(150, 2, 0.5)] * 2)
        (N, segments, diag), _ = run(pieces)
        self.assertAlmostEqual(N / N_TRUE, 1.0, delta=0.03)
        off = segments[~segments.fits]
        self.assertGreater(len(off), 0)
        self.assertTrue(off[['C', 'minor']].isna().all().all())
        self.assertTrue(segments.loc[segments.fits & (segments.state != 'loh'), ['C', 'minor']].notna().all().all())

    def test_no_reference_state_gives_no_scale(self):
        pieces = [(300, 4, 1)] * 6 + [(300, 5, 1)] * 6
        (N, _, diag), _ = run(pieces)
        self.assertIsNone(N)
        self.assertEqual(diag['status'], 'insufficient_anchors')

    def test_repeat_free_marks_rows_touching_the_bed(self):
        df = pd.DataFrame({'chr': ['chr1'] * 3 + ['chr2'], 'st': [1, 100_001, 200_001, 1],
                           'nd': [100_000, 200_000, 300_000, 100_000]})
        with tempfile.NamedTemporaryFile('w', suffix='.bed') as bed:
            bed.write('chr1\t150000\t150010\tsd\nchr1\t299999\t400000\tsd\n')
            bed.flush()
            self.assertEqual(PB.repeat_free(df, bed.name).tolist(), [True, False, False, True])

    def test_reads_without_base_qualities_count_nothing(self):
        import pysam
        with tempfile.TemporaryDirectory() as tmp:
            path = str(Path(tmp) / 'sample.bam')
            header = {'HD': {'VN': '1.6', 'SO': 'coordinate'}, 'SQ': [{'SN': 'chr1', 'LN': 1000}]}
            with pysam.AlignmentFile(path, 'wb', header=header) as out:
                for name, qual in (('with_q', [30] * 8), ('no_q', None)):
                    r = pysam.AlignedSegment(out.header)
                    r.query_name, r.reference_id, r.reference_start = name, 0, 10
                    r.mapping_quality, r.cigartuples = 60, [(0, 8)]
                    r.query_sequence = 'ACGTACGT'
                    if qual is not None:
                        r.query_qualities = pysam.qualitystring_to_array(''.join(chr(q + 33) for q in qual))
                    out.write(r)
            pysam.index(path)
            with pysam.AlignmentFile(path) as bam:
                counts = PB.allele_counts(bam, 'chr1', 0, 30)
                expected = np.array(bam.count_coverage('chr1', 0, 30, quality_threshold=PB.MIN_BQ,
                                                       read_callback='all'))
        np.testing.assert_array_equal(counts, expected)
        self.assertEqual(int(counts.sum()), 8)

    def test_missing_or_unindexed_bam_is_reported(self):
        self.assertIn('no sample BAM', PB.bam_problem('/nonexistent/sample.bam'))
        with tempfile.TemporaryDirectory() as tmp:
            bam = Path(tmp) / 'sample.bam'
            bam.write_bytes(b'')
            self.assertIn('not indexed', PB.bam_problem(str(bam)))
            (Path(tmp) / 'sample.bam.bai').write_bytes(b'')
            self.assertEqual(PB.bam_problem(str(bam)), '')


if __name__ == '__main__':
    unittest.main()
