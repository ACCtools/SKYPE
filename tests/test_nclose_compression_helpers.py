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
        "nclose_compression_layout",
        "nclose_merge_alignment",
    }
    tree = ast.parse(BUILD_GRAPH_PATH.read_text(encoding="utf-8"))
    functions = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in function_names
    ]
    module = ast.Module(body=functions, type_ignores=[])
    namespace = {"chr2int": lambda chrom: int(chrom[3:])}
    exec(compile(module, str(BUILD_GRAPH_PATH), "exec"), namespace)
    return namespace


HELPERS = load_nclose_compression_helpers()
nclose_compression_layout = HELPERS["nclose_compression_layout"]
nclose_merge_alignment = HELPERS["nclose_merge_alignment"]


class NCloseCompressionHelperTests(unittest.TestCase):
    def test_cross_chromosome_layout_stores_lower_chromosome_first(self):
        chr2 = ("=", "chr2")
        chr10 = ("=", "chr10")

        self.assertEqual(
            nclose_compression_layout(
                chr2, chr10, 4, 9, (900, 1000), (100, 200)
            ),
            ((chr2, chr10), (4, 9), False),
        )
        self.assertEqual(
            nclose_compression_layout(
                chr10, chr2, 4, 9, (100, 200), (900, 1000)
            ),
            ((chr2, chr10), (9, 4), False),
        )

    def test_same_chromosome_layout_uses_interval_order(self):
        chrom = ("=", "chr7")

        self.assertEqual(
            nclose_compression_layout(
                chrom, chrom, 4, 9, (100, 200), (900, 1000)
            ),
            ((chrom, chrom), (4, 9), True),
        )
        self.assertEqual(
            nclose_compression_layout(
                chrom, chrom, 4, 9, (100, 200), (100, 200)
            ),
            ((chrom, chrom), (4, 9), True),
        )
        self.assertEqual(
            nclose_compression_layout(
                chrom, chrom, 4, 9, (900, 1000), (100, 200)
            ),
            ((chrom, chrom), (9, 4), True),
        )

    def test_cross_chromosome_merge_uses_stored_endpoint_order(self):
        representative = ["rep", [10, 20], [30, 40], 9, 4]

        self.assertEqual(
            nclose_merge_alignment((4, 9), (9, 4), False, representative),
            ((9, 4), (1, 2)),
        )

    def test_same_chromosome_merge_follows_representative_path_order(self):
        forward_rep = ["forward", [10, 20], [30, 40], 4, 9]
        reverse_rep = ["reverse", [30, 40], [10, 20], 9, 4]

        self.assertEqual(
            nclose_merge_alignment((11, 15), (11, 15), True, forward_rep),
            ((11, 15), (1, 2)),
        )
        self.assertEqual(
            nclose_merge_alignment((11, 15), (15, 11), True, reverse_rep),
            ((11, 15), (2, 1)),
        )


if __name__ == "__main__":
    unittest.main()
