from __future__ import annotations

import csv
import sys
import tempfile
import unittest
from pathlib import Path

import matplotlib


matplotlib.use("Agg")

SKYPE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKYPE_ROOT))

from reference_path_clustering import (  # noqa: E402
    ReferenceSpan,
    build_path_overlap_edges,
    build_shared_reference_blocks,
    cluster_reference_paths,
    extract_ordered_reference_spans,
    is_full_reference_path,
    merge_intervals,
)
from reference_path_plotting import render_reference_path_clusters  # noqa: E402


def ppc_row(strand, chrom, start, end):
    row = [None] * 9
    row[4] = strand
    row[5] = chrom
    row[7] = start
    row[8] = end
    return tuple(row)


class CompactReferenceSpanTests(unittest.TestCase):
    def test_extracts_ordered_forward_and_reverse_runs_without_mutation(self) -> None:
        compact_path = [
            ("chr1f", 0, 0),
            (0, 0, 0, 0),
            (0, 1, 0, 0),
            (0, 2, 1, 0),
            ("chr2b", 1, 0),
        ]
        original = list(compact_path)
        ppc_data = [
            ppc_row("+", "chr1", 100, 200),
            ppc_row("+", "chr1", 300, 400),
            ppc_row("-", "chr2", 800, 900),
        ]

        spans = extract_ordered_reference_spans(compact_path, ppc_data)

        self.assertEqual(compact_path, original)
        self.assertEqual(spans, (
            ReferenceSpan("chr1", 100, 400, "+"),
            ReferenceSpan("chr2", 800, 900, "-"),
        ))

    def test_forced_break_edge_splits_a_counter_stable_path(self) -> None:
        compact_path = [
            ("chr1f", 0, 0),
            (0, 0, 0, 0),
            (1, 1, 0, 0),
            ("chr1b", 0, 0),
        ]
        ppc_data = [
            ppc_row("+", "chr1", 0, 100),
            ppc_row("+", "chr1", 300, 400),
        ]

        spans = extract_ordered_reference_spans(
            compact_path,
            ppc_data,
            forced_break_edges={((0, 0), (1, 1))},
        )

        self.assertEqual(spans, (
            ReferenceSpan("chr1", 0, 100, "+"),
            ReferenceSpan("chr1", 300, 400, "+"),
        ))

    def test_interrupt_spans_are_inserted_at_their_graph_edge(self) -> None:
        compact_path = [
            ("chr1f", 0, 0),
            (0, 0, 0, 0),
            (1, 1, 0, 0),
            ("chr1b", 0, 0),
        ]
        ppc_data = [
            ppc_row("+", "chr1", 0, 100),
            ppc_row("+", "chr1", 300, 400),
        ]
        interrupt = ReferenceSpan("chr6", 50_000, 250_000, "+", "interrupt")

        spans = extract_ordered_reference_spans(
            compact_path,
            ppc_data,
            interrupt_spans_by_edge={((0, 0), (1, 1)): (interrupt,)},
        )

        self.assertEqual(spans, (
            ReferenceSpan("chr1", 0, 100, "+"),
            interrupt,
            ReferenceSpan("chr1", 300, 400, "+"),
        ))


class ReferenceOverlapTests(unittest.TestCase):
    def test_merge_prevents_repeated_traversal_from_inflating_coverage(self) -> None:
        self.assertEqual(
            merge_intervals([(0, 200_000), (50_000, 150_000), (200_000, 220_000)]),
            ((0, 220_000),),
        )

    def test_exact_100kb_qualifies_but_99999_does_not(self) -> None:
        exact = {
            1: (ReferenceSpan("chr1", 0, 100_000, "+"),),
            2: (ReferenceSpan("chr1", 0, 100_000, "-"),),
        }
        short = {
            1: (ReferenceSpan("chr1", 0, 99_999, "+"),),
            2: (ReferenceSpan("chr1", 0, 99_999, "-"),),
        }

        self.assertEqual(len(build_path_overlap_edges(exact, minimum_overlap=100_000)), 1)
        self.assertEqual(build_path_overlap_edges(short, minimum_overlap=100_000), ())

    def test_disjoint_small_intersections_are_not_summed(self) -> None:
        path_spans = {
            1: (
                ReferenceSpan("chr1", 0, 60_000, "+"),
                ReferenceSpan("chr1", 200_000, 260_000, "+"),
            ),
            2: (
                ReferenceSpan("chr1", 0, 60_000, "-"),
                ReferenceSpan("chr1", 200_000, 260_000, "-"),
            ),
        }

        self.assertEqual(
            build_path_overlap_edges(path_spans, minimum_overlap=100_000),
            (),
        )

    def test_adjacent_reference_runs_separated_in_path_are_not_combined(self) -> None:
        path_spans = {
            1: (
                ReferenceSpan("chr1", 0, 60_000, "+"),
                ReferenceSpan("chr2", 0, 200_000, "+"),
                ReferenceSpan("chr1", 60_000, 120_000, "+"),
            ),
            2: (ReferenceSpan("chr1", 0, 120_000, "+"),),
        }

        self.assertEqual(
            build_path_overlap_edges(path_spans, minimum_overlap=100_000),
            (),
        )

    def test_overlap_requires_the_same_chromosome(self) -> None:
        path_spans = {
            1: (ReferenceSpan("chr1", 0, 1_000_000, "+"),),
            2: (ReferenceSpan("chr2", 0, 1_000_000, "+"),),
        }
        self.assertEqual(
            build_path_overlap_edges(path_spans, minimum_overlap=100_000),
            (),
        )

    def test_transitive_bridge_forms_one_component(self) -> None:
        path_spans = {
            1: (ReferenceSpan("chr1", 0, 200_000, "+"),),
            2: (
                ReferenceSpan("chr1", 0, 200_000, "+"),
                ReferenceSpan("chr2", 0, 200_000, "+"),
            ),
            3: (ReferenceSpan("chr2", 0, 200_000, "+"),),
            4: (ReferenceSpan("chr4", 0, 200_000, "+"),),
        }

        clusters, singletons = cluster_reference_paths(
            path_spans, minimum_overlap=100_000
        )

        self.assertEqual([(item.cluster_id, item.members) for item in clusters], [
            ("C01", (1, 2, 3)),
        ])
        self.assertEqual(singletons, (4,))

    def test_shared_block_does_not_add_an_incidental_subthreshold_carrier(self) -> None:
        path_spans = {
            1: (ReferenceSpan("chr1", 0, 200_000, "+"),),
            2: (
                ReferenceSpan("chr1", 0, 200_000, "+"),
                ReferenceSpan("chr2", 0, 200_000, "+"),
            ),
            3: (
                ReferenceSpan("chr1", 0, 50_000, "+"),
                ReferenceSpan("chr2", 0, 200_000, "+"),
            ),
        }
        clusters, _singletons = cluster_reference_paths(
            path_spans, minimum_overlap=100_000
        )

        blocks = build_shared_reference_blocks(clusters[0], path_spans)

        chr1_blocks = [block for block in blocks if block.chrom == "chr1"]
        self.assertEqual(len(chr1_blocks), 1)
        self.assertEqual(chr1_blocks[0].carriers, (1, 2))

    def test_full_reference_requires_empty_nclose_and_ratio_cutoff(self) -> None:
        spans = (ReferenceSpan("chr1", 10, 910, "+"),)

        self.assertEqual(
            is_full_reference_path(
                spans,
                has_nclose=False,
                chromosome_lengths={"chr1": 1000},
                minimum_ratio=0.9,
            ),
            (True, "chr1", 0.9),
        )
        self.assertFalse(is_full_reference_path(
            spans,
            has_nclose=True,
            chromosome_lengths={"chr1": 1000},
            minimum_ratio=0.9,
        )[0])
        self.assertFalse(is_full_reference_path(
            (ReferenceSpan("chr1", 10, 909, "+"),),
            has_nclose=False,
            chromosome_lengths={"chr1": 1000},
            minimum_ratio=0.9,
        )[0])
        self.assertFalse(is_full_reference_path(
            (
                ReferenceSpan("chr1", 0, 1_000, "+"),
                ReferenceSpan("chr2", 0, 200, "+"),
            ),
            has_nclose=False,
            chromosome_lengths={"chr1": 1000, "chr2": 1000},
            minimum_ratio=0.9,
        )[0])


class HG008ReferenceSpanRegressionTests(unittest.TestCase):
    def test_expected_four_non_singleton_components(self) -> None:
        s = ReferenceSpan
        path_spans = {
            80: (s("chr12", 1_356, 35_119_729, "+"), s("chr16", 36_118_033, 96_330_238, "+")),
            481: (
                s("chr18", 28_677_604, 80_324_602, "-"),
                s("chrX", 8_062_973, 94_256_226, "-"),
                s("chr2", 80_401_865, 89_333_754, "-"),
                s("chr2", 80_401_586, 82_711_740, "+"),
                s("chr20", 40_155_829, 66_209_441, "+"),
            ),
            802: (
                s("chr3", 196_600_634, 201_103_480, "-"),
                s("chr3", 142_746_299, 191_918_353, "+"),
                s("chr3", 93_035_054, 194_884_273, "+"),
                s("chr11", 68_236_657, 105_906_877, "+"),
                s("chr6", 159_219_001, 172_125_873, "+"),
            ),
            1700: (
                s("chrX", 362, 93_652_941, "+"),
                s("chrX", 89_080_411, 94_256_226, "+"),
                s("chr18", 28_677_604, 80_324_602, "+"),
            ),
            1893: (
                s("chr2", 80_401_865, 242_685_071, "-"),
                s("chr2", 80_401_586, 82_711_740, "+"),
                s("chr20", 40_155_829, 66_209_441, "+"),
            ),
            1934: (
                s("chr18", 176, 19_757_839, "+"),
                s("chr18", 21_678_664, 28_644_975, "-"),
                s("chr20", 25_564_158, 66_209_441, "+"),
            ),
            1963: (s("chr12", 1_356, 121_434_211, "+"), s("chr15", 38_258_149, 99_752_809, "+")),
            1969: (s("chr18", 176, 28_676_458, "+"), s("chrX", 94_258_293, 154_258_652, "+")),
            1974: (
                s("chr7", 1_000, 61_559_000, "+"),
                s("chr6", 58_286_685, 58_535_570, "+"),
                s("chr11", 125_690_000, 135_127_000, "+"),
            ),
            1997: (
                s("chr6", 9_000, 73_200_000, "+"),
                s("chr8", 45_875_000, 146_203_000, "+"),
            ),
            2004: (s("chr11", 59, 68_235_811, "+"), s("chr3", 194_884_277, 201_103_480, "+")),
            2008: (s("chr3", 907, 142_724_004, "+"), s("chr13", 43_549, 113_562_282, "-")),
        }

        clusters, singletons = cluster_reference_paths(
            path_spans, minimum_overlap=100_000
        )

        self.assertEqual([cluster.members for cluster in clusters], [
            (481, 1700, 1893, 1934, 1969),
            (802, 2004, 2008),
            (80, 1963),
            (1974, 1997),
        ])
        self.assertEqual(singletons, ())
        self.assertEqual([len(cluster.edges) for cluster in clusters], [5, 2, 1, 1])
        shared_blocks = build_shared_reference_blocks(clusters[0], path_spans)
        self.assertTrue(any(
            block.chrom == "chr20"
            and block.carriers == (481, 1893, 1934)
            for block in shared_blocks
        ))

class ReferencePathPlottingTests(unittest.TestCase):
    def test_renderer_writes_one_figure_per_cluster_and_auditable_tables(self) -> None:
        path_spans = {
            1: (ReferenceSpan("chr1", 0, 300_000, "+"),),
            2: (ReferenceSpan("chr1", 100_000, 400_000, "-"),),
        }
        clusters, singletons = cluster_reference_paths(
            path_spans, minimum_overlap=100_000
        )
        metadata = {
            1: {"location": "/run/a/1.paf", "raw_depth": 30.0, "depth_n": 1.0, "nclose_count": 1},
            2: {"location": "/run/b/2.paf", "raw_depth": 15.0, "depth_n": 0.5, "nclose_count": 1},
        }

        with tempfile.TemporaryDirectory() as tmpdir:
            result = render_reference_path_clusters(
                output_prefix=tmpdir,
                cell_line="TEST",
                clusters=clusters,
                singletons=singletons,
                path_spans=path_spans,
                path_metadata=metadata,
                excluded_full_reference=[{
                    "path_column": 9,
                    "location": "/run/full/9.paf",
                    "raw_depth": 60.0,
                    "depth_n": 2.0,
                    "nclose_count": 0,
                    "full_reference_chrom": "chr9",
                    "full_reference_ratio": 0.99,
                }],
                minimum_overlap=100_000,
                full_reference_ratio=0.9,
            )

            self.assertEqual(result["cluster_count"], 1)
            self.assertTrue(Path(result["figures"][0]["png_path"]).is_file())
            self.assertTrue(Path(result["figures"][0]["pdf_path"]).is_file())
            with open(result["membership_path"], newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(
                [row["status"] for row in rows],
                ["clustered", "clustered", "excluded_full_reference"],
            )
            with open(result["overlap_path"], newline="", encoding="utf-8") as handle:
                overlap_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(overlap_rows[0]["overlap_bp"], "200000")


if __name__ == "__main__":
    unittest.main()
