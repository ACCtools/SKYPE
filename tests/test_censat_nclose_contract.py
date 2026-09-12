"""Endpoint workflow boundaries replacing the removed legacy BOTH filter."""
import tempfile
import unittest
from pathlib import Path

from censat_endpoints import consistency, partition_unitigs, terminal_position, trace_source


class CensatEndpointTests(unittest.TestCase):
    def test_strand_aware_terminal_bases_are_half_open(self):
        row = dict(rs=100, re=200, strand='+')
        self.assertEqual(terminal_position(row, 'left'), 100)
        self.assertEqual(terminal_position(row, 'right'), 199)
        row['strand'] = '-'
        self.assertEqual(terminal_position(row, 'left'), 199)
        self.assertEqual(terminal_position(row, 'right'), 100)

    def test_partition_uses_all_bed_intervals_and_excludes_same_state_owners(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            bed = root / 'censat.bed'
            bed.write_text('chr1\t100\t200\nchr2\t300\t400\n')
            paf = root / 'utg.aln.paf'
            # Reverse query row order; low MAPQ, short intervals, and query
            # overhangs do not invalidate the new route. Internal non-censat
            # chunks do not affect partition membership.
            paf.write_text(
                'both\t1000\t800\t900\t-\tchr2\t10000\t300\t400\t100\t100\t0\n'
                'both\t1000\t100\t200\t+\tchr1\t10000\t100\t200\t100\t100\t0\n'
                'both\t1000\t300\t400\t+\tchr3\t10000\t500\t600\t100\t100\t0\n'
                'same\t100\t0\t100\t+\tchr1\t10000\t100\t200\t100\t100\t60\n'
                'legacy\t100\t0\t100\t+\tchr1\t10000\t150\t250\t100\t100\t60\n'
            )
            partition, candidates = partition_unitigs(paf, bed)
        self.assertEqual({r['unitig']: r['route'] for r in partition},
                         {'both': 'censat', 'same': 'censat', 'legacy': 'legacy'})
        self.assertEqual([r['unitig'] for r in candidates], ['both'])
        self.assertEqual(candidates[0]['left']['index'], 1)
        self.assertEqual(candidates[0]['right']['index'], 0)

    def test_half_chunk_per_alignment_and_source_offset(self):
        chunk = dict(qs=1000, qe=1200, chrom='chr1', strand='+')
        hit = dict(qs=100, qe=200, chrom='chr1', strand='+', tp='P')
        self.assertEqual(consistency(chunk, [hit], 900)['status'], 'consistent')
        conflict = dict(hit, strand='-', tp='S')
        self.assertEqual(consistency(chunk, [hit, conflict], 900)['status'], 'conflict')
        # Two short pieces are not added together to manufacture 50% coverage.
        short = [dict(hit, qs=100, qe=150), dict(hit, qs=150, qe=200)]
        self.assertEqual(consistency(chunk, short, 900)['status'], 'no_alignment')
        self.assertEqual(consistency(chunk, [], 900)['status'], 'no_alignment')

    def test_xi_A_missing_and_mismatched_sources_are_unassessed(self):
        chunk = dict(name='u', chrom='chr1', strand='+', qs=100, qe=200, xi='P_7')
        source = dict(chunk, qs=50, qe=250)
        self.assertEqual(trace_source(chunk, {7: source}), (source, None))
        for xi in ('A_7', None, 'P_8'):
            self.assertIsNotNone(trace_source(dict(chunk, xi=xi), {7: source})[1])
        self.assertEqual(trace_source(chunk, {7: dict(source, strand='-')})[1],
                         'source_mismatch')


if __name__ == '__main__':
    unittest.main()
