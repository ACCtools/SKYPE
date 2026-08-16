from __future__ import annotations

import ast
import unittest
from pathlib import Path


BUILD_GRAPH_PATH = (
    Path(__file__).resolve().parents[1] / "02_Build_Breakend_Graph_Limited.py"
)


def load_nclose_compression_helpers():
    """Load the pure helpers without running stage 02's argparse path."""
    function_names = {
        "_flip_ctg_dir",
        "nclose_canonical_directions",
        "nclose_compression_candidate_matches",
        "nclose_compression_layout",
        "nclose_merge_direction_matches",
        "nclose_merge_alignment",
    }
    tree = ast.parse(BUILD_GRAPH_PATH.read_text(encoding="utf-8"))
    functions = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in function_names
    ]
    module = ast.Module(body=functions, type_ignores=[])
    namespace = {
        "CTG_DIR": 4,
        "chr2int": lambda chrom: int(chrom[3:]),
        "distance_checker": lambda _first, _second: 0,
    }
    exec(compile(module, str(BUILD_GRAPH_PATH), "exec"), namespace)
    return namespace


HELPERS = load_nclose_compression_helpers()
nclose_compression_layout = HELPERS["nclose_compression_layout"]
nclose_merge_alignment = HELPERS["nclose_merge_alignment"]
nclose_canonical_directions = HELPERS["nclose_canonical_directions"]
nclose_compression_candidate_matches = HELPERS[
    "nclose_compression_candidate_matches"
]
nclose_merge_direction_matches = HELPERS["nclose_merge_direction_matches"]


def contig_rows(path_nodes, path_directions):
    rows = [[None] * 5 for _ in range(max(path_nodes) + 1)]
    for node_idx, direction in zip(path_nodes, path_directions):
        rows[node_idx][4] = direction
    return rows


class NCloseCompressionHelperTests(unittest.TestCase):
    def test_cross_chromosome_layout_stores_lower_chromosome_first(self):
        chr2 = ("=", "chr2")
        chr10 = ("=", "chr10")

        self.assertEqual(
            nclose_compression_layout(
                chr2, chr10, 4, 9, (900, 1000), (100, 200)
            ),
            ((chr2, chr10), (4, 9), False, True),
        )
        self.assertEqual(
            nclose_compression_layout(
                chr10, chr2, 4, 9, (100, 200), (900, 1000)
            ),
            ((chr2, chr10), (9, 4), False, False),
        )

    def test_same_chromosome_layout_uses_interval_order(self):
        chrom = ("=", "chr7")

        self.assertEqual(
            nclose_compression_layout(
                chrom, chrom, 4, 9, (100, 200), (900, 1000)
            ),
            ((chrom, chrom), (4, 9), True, True),
        )
        self.assertEqual(
            nclose_compression_layout(
                chrom, chrom, 4, 9, (100, 200), (100, 200)
            ),
            ((chrom, chrom), (4, 9), True, True),
        )
        self.assertEqual(
            nclose_compression_layout(
                chrom, chrom, 4, 9, (900, 1000), (100, 200)
            ),
            ((chrom, chrom), (9, 4), True, False),
        )

    def test_cross_chromosome_merge_uses_stored_endpoint_order(self):
        representative = [
            "rep", [10, 20], [30, 40], 9, 4, ("+", "+"), False
        ]

        self.assertEqual(
            nclose_merge_alignment((4, 9), (9, 4), False, representative),
            ((9, 4), (1, 2)),
        )

    def test_same_chromosome_merge_follows_representative_path_order(self):
        forward_rep = [
            "forward", [10, 20], [30, 40], 4, 9, ("+", "+"), True
        ]
        reverse_rep = [
            "reverse", [30, 40], [10, 20], 9, 4, ("+", "+"), False
        ]

        self.assertEqual(
            nclose_merge_alignment((11, 15), (11, 15), True, forward_rep),
            ((11, 15), (1, 2)),
        )
        self.assertEqual(
            nclose_merge_alignment((11, 15), (15, 11), True, reverse_rep),
            ((11, 15), (2, 1)),
        )

    def test_merge_alignment_uses_explicit_representative_path_flag(self):
        descending_ids_but_forward = [
            "forward", [10, 20], [30, 40], 9, 4, ("+", "+"), True
        ]
        ascending_ids_but_reverse = [
            "reverse", [10, 20], [30, 40], 4, 9, ("+", "+"), False
        ]

        self.assertEqual(
            nclose_merge_alignment(
                (11, 15), (11, 15), True, descending_ids_but_forward
            ),
            ((11, 15), (1, 2)),
        )
        self.assertEqual(
            nclose_merge_alignment(
                (11, 15), (11, 15), True, ascending_ids_but_reverse
            ),
            ((11, 15), (2, 1)),
        )

    def test_canonical_direction_matrix_uses_explicit_path_order(self):
        cases = [
            ((4, 9), (4, 9), ("+", "+"), ("+", "+")),
            ((4, 9), (4, 9), ("+", "-"), ("+", "-")),
            ((4, 9), (4, 9), ("-", "+"), ("-", "+")),
            ((4, 9), (4, 9), ("-", "-"), ("-", "-")),
            ((9, 4), (4, 9), ("+", "+"), ("-", "-")),
            ((9, 4), (4, 9), ("+", "-"), ("+", "-")),
            ((9, 4), (4, 9), ("-", "+"), ("-", "+")),
            ((9, 4), (4, 9), ("-", "-"), ("+", "+")),
        ]

        for path_nodes, stored_nodes, path_dirs, expected in cases:
            with self.subTest(
                path_nodes=path_nodes,
                path_dirs=path_dirs,
            ):
                self.assertEqual(
                    nclose_canonical_directions(
                        contig_rows(path_nodes, path_dirs),
                        path_nodes,
                        stored_nodes,
                    ),
                    expected,
                )

    def test_equal_interval_direction_signature_is_rc_canonical(self):
        forward = nclose_canonical_directions(
            contig_rows((4, 9), ("+", "+")),
            (4, 9),
            (4, 9),
            canonical_order_tied=True,
        )
        reverse = nclose_canonical_directions(
            contig_rows((9, 4), ("-", "-")),
            (9, 4),
            (9, 4),
            canonical_order_tied=True,
        )

        self.assertEqual(forward, ("+", "+"))
        self.assertEqual(reverse, forward)

        tied_forward_inversion = nclose_canonical_directions(
            contig_rows((4, 9), ("+", "-")),
            (4, 9),
            (4, 9),
            canonical_order_tied=True,
        )
        tied_reverse_inversion = nclose_canonical_directions(
            contig_rows((4, 9), ("-", "+")),
            (4, 9),
            (4, 9),
            canonical_order_tied=True,
        )
        self.assertEqual(tied_forward_inversion, ("+", "-"))
        self.assertEqual(tied_reverse_inversion, ("-", "+"))
        self.assertNotEqual(tied_forward_inversion, tied_reverse_inversion)

    def test_merge_direction_is_checked_per_endpoint(self):
        rows = [[None] * 5 for _ in range(12)]
        rows[4][4] = "+"
        rows[9][4] = "-"
        rows[10][4] = "+"
        rows[11][4] = "-"
        representative = [
            "representative",
            [100, 200],
            [300, 400],
            4,
            9,
            ("+", "-"),
            True,
        ]

        self.assertTrue(
            nclose_merge_direction_matches(
                rows,
                10,
                representative,
                1,
            )
        )
        self.assertFalse(
            nclose_merge_direction_matches(
                rows,
                11,
                representative,
                1,
            )
        )

        # A reverse representative aligns its path start to payload slot 2.
        representative[6] = False
        self.assertTrue(
            nclose_merge_direction_matches(
                rows,
                11,
                representative,
                representative_payload_slot=2,
            )
        )

    def test_coordinate_match_does_not_compress_opposite_signatures(self):
        rows = contig_rows((4, 9), ("+", "-"))
        representative = [
            "representative",
            [100, 200],
            [300, 400],
            4,
            9,
            ("-", "+"),
            True,
        ]

        self.assertFalse(
            nclose_compression_candidate_matches(
                rows,
                (4, 9),
                ("+", "-"),
                representative,
                100_000,
            )
        )
        representative[5] = ("+", "-")
        self.assertTrue(
            nclose_compression_candidate_matches(
                rows,
                (4, 9),
                ("+", "-"),
                representative,
                100_000,
            )
        )


if __name__ == "__main__":
    unittest.main()
