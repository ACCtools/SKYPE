import unittest

from nclose_read_geometry import (
    aligned_depth_interval,
    interval_weighted_mean,
    reference_support_points,
    reference_support_span_rows,
    split_nclose_flanks,
)


class NCloseReadGeometryTests(unittest.TestCase):
    def test_user_minus_to_plus_example(self):
        self.assertEqual(
            reference_support_points(3000, "-", 0, 1000, 1_000_000),
            (3000, 2000),
        )
        self.assertEqual(
            reference_support_points(2000, "+", 1, 1000, 1_000_000),
            (1000, 2000),
        )
        self.assertEqual(
            aligned_depth_interval(3000, "-", 0, 500_000, 1_000_000),
            (3001, 503000),
        )
        self.assertEqual(
            aligned_depth_interval(2000, "+", 1, 500_000, 1_000_000),
            (2001, 502000),
        )

    def test_fixed_window_flank_split(self):
        self.assertEqual(split_nclose_flanks(0, 5000), (2500, 2500))
        self.assertEqual(split_nclose_flanks(533, 5000), (2233, 2234))
        self.assertEqual(split_nclose_flanks(1000, 5000), (2000, 2000))
        self.assertEqual(split_nclose_flanks(5000, 5000), (0, 0))
        self.assertIsNone(split_nclose_flanks(5001, 5000))

    def test_production_window_reuses_shifted_d2_anchors(self):
        # With a 1 kb gap in the 5 kb count window, the unchanged d2 logic
        # places both anchors 2 kb into their aligned reference chunks.
        first_d2_anchor = 1_002_000  # first endpoint is '-'
        second_d2_anchor = 2_002_000  # second endpoint is '+'

        self.assertEqual(
            reference_support_points(
                first_d2_anchor, "-", 0, 5000, 3_000_000
            ),
            (1_002_000, 997000),
        )
        self.assertEqual(
            reference_support_points(
                second_d2_anchor, "+", 1, 5000, 3_000_000
            ),
            (1_997_000, 2_002_000),
        )

    def test_all_strand_combinations_follow_query_direction(self):
        shared_a = 1_000_000
        shared_b = 2_000_000
        for strand in ("+", "-"):
            with self.subTest(side="a", strand=strand):
                start, end = reference_support_points(
                    shared_a, strand, 0, 5000, 3_000_000
                )
                self.assertEqual(start, shared_a)
                self.assertEqual(end > start, strand == "+")
            with self.subTest(side="b", strand=strand):
                start, end = reference_support_points(
                    shared_b, strand, 1, 5000, 3_000_000
                )
                self.assertEqual(end, shared_b)
                self.assertEqual(end > start, strand == "+")

    def test_span_rows_keep_query_order_and_expected_strand(self):
        shared_a = [7, 102, 1, "chr2", 999_000, 1_002_000, 1_002_000, False]
        d1_rows = reference_support_span_rows(
            7, 101, shared_a, "-", 0, 5000, 3_000_000
        )
        self.assertEqual([row[6] for row in d1_rows], [1_002_000, 997000])
        self.assertEqual([row[7] for row in d1_rows], [False, False])
        self.assertEqual(d1_rows[0][6], shared_a[6])

        shared_b = [7, 102, 2, "chr3", 2_000_000, 2_002_000, 2_002_000, True]
        d3_rows = reference_support_span_rows(
            7, 103, shared_b, "+", 1, 5000, 3_000_000
        )
        self.assertEqual([row[6] for row in d3_rows], [1_997_000, 2_002_000])
        self.assertEqual([row[7] for row in d3_rows], [True, True])
        self.assertEqual(d3_rows[1][6], shared_b[6])

    def test_chromosome_boundary_clips_only_nonshared_point(self):
        self.assertEqual(
            reference_support_points(3000, "-", 0, 5000, 1_000_000),
            (3000, 1),
        )
        self.assertEqual(
            reference_support_points(2000, "+", 1, 5000, 1_000_000),
            (1, 2000),
        )
        self.assertEqual(
            reference_support_points(999_000, "+", 0, 5000, 1_000_000),
            (999_000, 1_000_000),
        )

    def test_local_depth_uses_only_the_aligned_reference_side(self):
        depth_rows = [
            (1, 500_000, 10.0),
            (500_001, 1_000_000, 30.0),
        ]
        cases = [
            ("+", 0, (1, 500_000), 10.0),
            ("-", 0, (500_001, 1_000_000), 30.0),
            ("+", 1, (500_001, 1_000_000), 30.0),
            ("-", 1, (1, 500_000), 10.0),
        ]
        for strand, side_index, expected_interval, expected_depth in cases:
            with self.subTest(strand=strand, side_index=side_index):
                interval = aligned_depth_interval(
                    500_000, strand, side_index, 500_000, 1_000_000
                )
                self.assertEqual(interval, expected_interval)
                self.assertEqual(
                    interval_weighted_mean(depth_rows, interval), expected_depth
                )


if __name__ == "__main__":
    unittest.main()
