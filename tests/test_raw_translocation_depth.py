import unittest

from raw_translocation_depth import (
    breakpoint_flank_depths,
    build_raw_point_depths,
    depth_pair_is_balanced,
    depth_pair_mean,
    weighted_mean_depth,
)


class RawTranslocationDepthTests(unittest.TestCase):
    def test_breakpoint_flanks_are_disjoint_one_sided_windows(self):
        depth_by_chrom = {
            "chr1": [
                (1, 500_000, 48.0),
                (500_001, 1_000_000, 96.0),
            ],
        }

        depths = breakpoint_flank_depths(
            depth_by_chrom,
            "chr1",
            500_001,
            500_001,
            1_000_000,
            flank=500_000,
        )

        self.assertEqual(depths, {"front": 48.0, "back": 96.0})
        self.assertFalse(depth_pair_is_balanced(depths, 40.0))
        self.assertTrue(depth_pair_is_balanced(depths, 50.0))
        self.assertEqual(depth_pair_mean(depths), 72.0)

    def test_breakpoint_interval_is_excluded_from_both_flanks(self):
        depth_by_chrom = {
            "chr1": [
                (1, 100, 10.0),
                (101, 199, 1_000.0),
                (200, 300, 20.0),
            ],
        }

        depths = breakpoint_flank_depths(
            depth_by_chrom,
            "chr1",
            101,
            199,
            300,
            flank=100,
        )

        self.assertEqual(depths, {"front": 10.0, "back": 20.0})

    def test_candidate_depths_use_inner_breakpoint_intervals(self):
        candidates = [{
            "pair_id": 3,
            "side_records": [
                {
                    "chrom": "chr1",
                    "inner_st": 500_001,
                    "inner_nd": 500_015,
                    "ref_count_st": 497_500,
                    "ref_count_nd": 502_500,
                },
                {
                    "chrom": "chr8",
                    "inner_st": 200_001,
                    "inner_nd": 200_100,
                    "ref_count_st": 197_500,
                    "ref_count_nd": 202_500,
                },
            ],
        }]
        depth_by_chrom = {
            "chr1": [(1, 500_000, 48.0), (500_001, 1_000_000, 96.0)],
            "chr8": [(1, 200_000, 160.0), (200_001, 400_000, 200.0)],
        }

        observed = build_raw_point_depths(
            candidates,
            depth_by_chrom,
            {"chr1": 1_000_000, "chr8": 400_000},
            flank=100_000,
        )

        self.assertEqual(observed[3]["point_a"], {
            "front": 48.0,
            "back": 96.0,
        })
        self.assertEqual(observed[3]["point_b"], {
            "front": 160.0,
            "back": 200.0,
        })

    def test_weighted_mean_and_balance_handle_partial_or_missing_flanks(self):
        depth_by_chrom = {
            "chr1": [(1, 50, 10.0), (51, 100, 20.0)],
        }
        self.assertEqual(
            weighted_mean_depth(depth_by_chrom, "chr1", 26, 75),
            15.0,
        )
        self.assertIsNone(
            weighted_mean_depth(depth_by_chrom, "chr2", 1, 100)
        )
        self.assertFalse(depth_pair_is_balanced({"front": None, "back": 10}, 15.0))
        self.assertTrue(depth_pair_is_balanced({"front": 100, "back": 109}, 15.0))


if __name__ == "__main__":
    unittest.main()
