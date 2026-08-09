import importlib.util
import os
import unittest
from collections import Counter


MODULE_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "26_nclose_read_depth_plot.py",
)
SPEC = importlib.util.spec_from_file_location("nclose_read_depth_plot", MODULE_PATH)
depth_plot = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(depth_plot)


class NCloseReadDepthPlotTests(unittest.TestCase):
    def test_normalized_read_depth_uses_raw_coverage_and_supports_old_pickles(self):
        self.assertAlmostEqual(
            depth_plot.normalized_read_depth(
                {"weighted_read_depth": 69.0}, 34.5
            ),
            2.0,
        )
        self.assertAlmostEqual(
            depth_plot.normalized_read_depth(
                {"weighted_read_depth_N": 1.75}, 34.5
            ),
            1.75,
        )

    def test_scatter_category_marks_problematic_rows(self):
        self.assertEqual(
            depth_plot.scatter_category({"status": "AGREE_WITHIN_2X"}),
            "ASSESSABLE",
        )
        self.assertEqual(
            depth_plot.scatter_category({
                "status": "AGREE_WITHIN_2X",
                "current_depth_N": 1.0,
                "read_depth_N": 1.3,
            }),
            "ASSESSABLE_DISCORDANT",
        )
        self.assertEqual(
            depth_plot.scatter_category({
                "status": "UNASSESSABLE_REPEAT_OR_LOWMAPQ"
            }),
            "REPEAT_OR_LOWMAPQ",
        )
        self.assertEqual(
            depth_plot.scatter_category({"status": "HIGH_QUALITY_ZERO"}),
            "ZERO_SUPPORT",
        )

    def test_scatter_excludes_any_censat_overlap(self):
        self.assertFalse(depth_plot.scatter_is_included({
            "overlaps_censat": True,
        }))
        self.assertTrue(depth_plot.scatter_is_included({
            "overlaps_censat": False,
        }))

    def test_connected_components_follow_path_co_carriage(self):
        event_a = (1, 2)
        event_b = (3, 4)
        event_c = (5, 6)
        path_usage = [
            Counter({event_a: 1, event_b: 1}),
            Counter({event_c: 1}),
        ]

        components, carriers = depth_plot.connected_components(
            {event_a, event_b, event_c}, [0, 1], path_usage
        )

        self.assertIn({event_a, event_b}, components)
        self.assertIn({event_c}, components)
        self.assertEqual(carriers[event_a], {0})
        self.assertEqual(carriers[event_c], {1})

    def test_proportionality_metrics_use_high_quality_rows(self):
        rows = [
            {
                "status": "AGREE_WITHIN_2X",
                "current_depth_N": 1.0,
                "read_depth_N": 1.1,
                "junction_reads": 11,
            },
            {
                "status": "AGREE_WITHIN_2X",
                "current_depth_N": 2.0,
                "read_depth_N": 2.2,
                "junction_reads": 22,
            },
            {
                "status": "UNASSESSABLE_REPEAT_OR_LOWMAPQ",
                "current_depth_N": 1.0,
                "read_depth_N": 0.0,
                "junction_reads": 100,
            },
        ]

        metrics = depth_plot.proportionality_metrics(rows)

        self.assertEqual(metrics["high_quality_positive_event_count"], 2)
        self.assertAlmostEqual(metrics["slope_through_origin"], 1.1)
        self.assertAlmostEqual(metrics["junction_reads_per_fitted_N_slope"], 11.0)
        self.assertAlmostEqual(metrics["median_read_to_current_ratio"], 1.1)
        self.assertEqual(metrics["rough_agreement_fraction"], 1.0)

if __name__ == "__main__":
    unittest.main()
