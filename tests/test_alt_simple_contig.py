from __future__ import annotations

import ast
import unittest
from collections import Counter
from pathlib import Path


BUILD_GRAPH_PATH = (
    Path(__file__).resolve().parents[1] / "02_Build_Breakend_Graph_Limited.py"
)


def load_alt_simple_helpers():
    """Load the pure simple_ctg_alt helpers without running stage 02."""
    function_names = {
        "alt_simple_ref_len",
        "alt_simple_chr_sort_key",
        "select_alt_simple_major_chroms",
        "alt_simple_alignment_state",
        "find_alt_simple_inward_bounds",
    }
    tree = ast.parse(BUILD_GRAPH_PATH.read_text(encoding="utf-8"))
    functions = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in function_names
    ]
    module = ast.Module(body=functions, type_ignores=[])
    namespace = {
        "Counter": Counter,
        "INF": 1_000_000_000,
        "CTG_DIR": 4,
        "CHR_NAM": 5,
        "CHR_STR": 7,
        "CHR_END": 8,
        "ALT_SIMPLE_MIN_SEGMENT_LEN": 10_000,
        "ALT_SIMPLE_MAJOR_CHR_RATIO": 0.90,
        "chr2int": lambda chrom: int(chrom[3:]),
    }
    exec(compile(module, str(BUILD_GRAPH_PATH), "exec"), namespace)
    return namespace


HELPERS = load_alt_simple_helpers()
select_alt_simple_major_chroms = HELPERS["select_alt_simple_major_chroms"]
find_alt_simple_inward_bounds = HELPERS["find_alt_simple_inward_bounds"]


def chunk(
    query_start: int,
    query_end: int,
    strand: str,
    chrom: str,
    ref_span: int,
) -> list:
    return [
        "ctg",
        1_000_000,
        query_start,
        query_end,
        strand,
        chrom,
        250_000_000,
        1_000_000,
        1_000_000 + ref_span,
        60,
        0,
    ]


class AltSimpleContigTests(unittest.TestCase):
    def test_terminal_chromosomes_are_union_with_90_percent_set(self) -> None:
        indexed = [
            (0, chunk(0, 950_000, "+", "chr6", 950_000)),
            (1, chunk(950_000, 990_000, "+", "chr7", 40_000)),
            (2, chunk(990_000, 1_010_000, "+", "chr11", 20_000)),
        ]

        selected = select_alt_simple_major_chroms(
            indexed, required_chroms=("chr7", "chr11")
        )

        self.assertEqual(selected, {"chr6", "chr7", "chr11"})

    def test_inward_search_ignores_noise_but_retains_it_inside_bounds(self) -> None:
        chunks = [
            chunk(0, 20_000, "+", "chr7", 20_000),
            chunk(20_000, 70_000, "-", "chr2", 50_000),   # non-major: ignore
            chunk(70_000, 100_000, "+", "chr7", 30_000),
            chunk(100_000, 109_000, "-", "chr7", 9_000),  # short: ignore
            chunk(109_000, 149_000, "+", "chr7", 40_000),
            chunk(149_000, 169_000, "-", "chr7", 20_000), # first real change
            chunk(169_000, 249_000, "+", "chr6", 80_000),
            chunk(249_000, 258_000, "-", "chr11", 9_000), # short: ignore
            chunk(258_000, 288_000, "+", "chr11", 30_000),
        ]
        indexed = list(enumerate(chunks))

        bounds = find_alt_simple_inward_bounds(
            indexed, {"chr6", "chr7", "chr11"}
        )

        self.assertIsNotNone(bounds)
        left, right = bounds
        self.assertEqual((left[0], right[0]), (4, 8))
        # All raw chunks between the inclusive bounds go to PPC, including
        # the direction-changing and short chunks at indices 5 and 7.
        self.assertEqual(
            list(range(left[0], right[0] + 1)),
            [4, 5, 6, 7, 8],
        )
        self.assertEqual((left[1][5], left[1][4]), ("chr7", "+"))
        self.assertEqual((right[1][5], right[1][4]), ("chr11", "+"))

    def test_same_post_telomere_terminal_chromosome_is_rejected(self) -> None:
        indexed = [
            (0, chunk(0, 20_000, "+", "chr7", 20_000)),
            (1, chunk(20_000, 40_000, "+", "chr6", 20_000)),
            (2, chunk(40_000, 60_000, "-", "chr7", 20_000)),
        ]

        self.assertIsNone(
            find_alt_simple_inward_bounds(indexed, {"chr6", "chr7"})
        )


if __name__ == "__main__":
    unittest.main()
