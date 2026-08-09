import importlib.util
import os
import unittest

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt

from virtual_sky_plotting import plot_virtual_chromosome, weighted_vaf_read_depth


MODULE_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "25_cluster_nclose_read_count.py",
)
SPEC = importlib.util.spec_from_file_location("cluster_nclose_read_count", MODULE_PATH)
cluster_count = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(cluster_count)


def ppc_row(name, q_start, q_end, strand, chrom, ref_start, ref_end, censat="0"):
    return [
        name, 1000, q_start, q_end, strand, chrom, 250_000_000,
        ref_start, ref_end, 60, 1, 0, 1, "0", "0", "0", "0", "0",
        censat, strand, chrom, "0.0",
    ]


class ClusterNCloseCandidateTests(unittest.TestCase):
    def test_candidates_follow_catalog_order_and_preserve_cluster_namespace(self):
        contig_data = [
            ppc_row("junction", 0, 100, "+", "chr1", 100, 200),
            ppc_row("junction", 100, 1000, "-", "chr5", 300, 400, "rin"),
        ]
        catalog = [
            {"event_key": ("front_jump", 1, -1), "kind": "indel"},
            {
                "event_key": (0, 1),
                "kind": "bnd",
                "contig_name": "junction",
            },
        ]

        candidates = cluster_count.build_cluster_nclose_candidates(
            contig_data, catalog, {(0, 1)}
        )

        self.assertEqual(len(candidates), 1)
        self.assertEqual(candidates[0]["pair_id"], 0)
        self.assertEqual(candidates[0]["nclose_id"], "SKYPE.nclose.2")
        self.assertEqual(candidates[0]["nclose_key"], (0, 1))
        self.assertTrue(candidates[0]["overlaps_censat"])
        self.assertEqual(candidates[0]["layout"]["chroms"], ("chr1", "chr5"))
        self.assertEqual(candidates[0]["layout"]["coords"], (200, 400))

    def test_graph_only_indel_uses_type4_node_pair_for_raw_counting(self):
        contig_data = [
            ppc_row("type4", 0, 100, "-", "chrX", 100, 200),
            ppc_row("type4", 100, 1000, "-", "chrX", 300, 400),
        ]
        event_key = ("type4_graph_indel", "type4_tuple", 0, 1)
        catalog = [{
            "event_key": event_key,
            "kind": "indel",
            "graph_only": True,
            "type4_tuple": (0, 1),
            "contig_name": "type4",
        }]

        candidates = cluster_count.build_cluster_nclose_candidates(
            contig_data, catalog, {event_key}
        )

        self.assertEqual(len(candidates), 1)
        self.assertEqual(candidates[0]["event_key"], event_key)
        self.assertEqual(candidates[0]["event_kind"], "indel")
        self.assertEqual(candidates[0]["nclose_key"], (0, 1))


class ClusterNClosePlotTests(unittest.TestCase):
    def test_weighted_vaf_read_depth_uses_local_coverage(self):
        record = {
            "read_counts": {"d1": 0, "d2": 30, "d3": 35},
            "vaf": {"chr_a": 1.0, "chr_b": 30 / 65},
            "read_count_depth_estimate": {
                "chr_a_local_depth": 72.07306802799827,
                "chr_b_local_depth": 102.4727703820461,
            },
        }

        depth = weighted_vaf_read_depth(record)

        self.assertAlmostEqual(depth, 55.119738445277164)

    def test_junction_label_prefers_vaf_read_depth(self):
        figure, axis = plt.subplots()
        try:
            plot_virtual_chromosome(
                axis,
                ("path", [(('chr18', '+'), 25.0), (('chrX', '-'), 75.0)]),
                maxh=100,
                cell_col=3,
                default_cell_col=1.8,
                junction_read_counts=[(1, 30)],
                junction_read_depths=[(1, 55.1197)],
            )
            self.assertIn("t(18;X) : 55.1x", [text.get_text() for text in axis.texts])
        finally:
            plt.close(figure)

    def test_junction_label_includes_raw_read_count(self):
        figure, axis = plt.subplots()
        try:
            plot_virtual_chromosome(
                axis,
                ("path", [(('chr1', '+'), 25.0), (('chr5', '-'), 75.0)]),
                maxh=100,
                cell_col=3,
                default_cell_col=1.8,
                junction_read_counts=[(1, 38)],
            )
            self.assertIn("t(1;5) : 38", [text.get_text() for text in axis.texts])
        finally:
            plt.close(figure)

    def test_uncountable_junction_is_not_reported_as_zero(self):
        figure, axis = plt.subplots()
        try:
            plot_virtual_chromosome(
                axis,
                ("path", [(('chr6', '+'), 25.0), (('chr8', '-'), 75.0)]),
                maxh=100,
                cell_col=3,
                default_cell_col=1.8,
                junction_read_counts=[(1, None)],
            )
            self.assertIn("t(6;8) : NA", [text.get_text() for text in axis.texts])
        finally:
            plt.close(figure)

    def test_uncountable_junction_depth_is_reported_as_na(self):
        figure, axis = plt.subplots()
        try:
            plot_virtual_chromosome(
                axis,
                ("path", [(('chr6', '+'), 25.0), (('chr8', '-'), 75.0)]),
                maxh=100,
                cell_col=3,
                default_cell_col=1.8,
                junction_read_depths=[(1, None)],
            )
            self.assertIn("t(6;8) : NA", [text.get_text() for text in axis.texts])
        finally:
            plt.close(figure)


if __name__ == "__main__":
    unittest.main()
