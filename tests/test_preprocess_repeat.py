from __future__ import annotations

import ast
import unittest
from pathlib import Path


BUILD_GRAPH_PATH = (
    Path(__file__).resolve().parents[1] / "02_Build_Breakend_Graph_Limited.py"
)


def load_preprocess_repeat():
    """Load preprocess_repeat without executing the stage-02 entry point."""
    tree = ast.parse(BUILD_GRAPH_PATH.read_text(encoding="utf-8"))
    function = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name == "preprocess_repeat"
    )
    namespace = {
        "CTG_NAM": 0,
        "CTG_MAPQ": 9,
        "CTG_TYP": 10,
        "CTG_STRND": 11,
        "CTG_ENDND": 12,
        "CTG_TELCON": 15,
        "CTG_RPTCHR": 16,
        "MAPQ_BOUND": 60,
    }
    module = ast.Module(body=[function], type_ignores=[])
    exec(compile(module, str(BUILD_GRAPH_PATH), "exec"), namespace)
    return namespace["preprocess_repeat"]


preprocess_repeat = load_preprocess_repeat()


def contig(*specs):
    """Build one type-1 contig from (chromosome, MAPQ, is_repeat) rows."""
    rows = []
    end_index = len(specs) - 1
    for index, (chromosome, mapq, is_repeat) in enumerate(specs):
        rows.append([
            "query",
            1_000,
            index * 100,
            (index + 1) * 100,
            "+",
            chromosome,
            1_000_000,
            10_000 + index * 100,
            10_100 + index * 100,
            mapq,
            1,
            0,
            end_index,
            "0",
            "0",
            "0",
            chromosome if is_repeat else "0",
            "r" if is_repeat else "0",
            "0",
            "+",
            chromosome,
            f"1.{index}",
        ])
    return rows


class PreprocessRepeatTests(unittest.TestCase):
    def test_back_high_mapq_repeat_restores_adjacent_low_mapq_repeat(self):
        rows = contig(
            ("chrX", 60, True),
            ("chr2", 43, True),
        )

        kept = preprocess_repeat(rows)

        self.assertEqual(
            [(row[5], row[9]) for row in kept],
            [("chrX", 60), ("chr2", 43)],
        )

    def test_front_high_mapq_repeat_restores_adjacent_low_mapq_repeat(self):
        rows = contig(
            ("chr2", 43, True),
            ("chrX", 60, True),
        )

        kept = preprocess_repeat(rows)

        self.assertEqual(
            [(row[5], row[9]) for row in kept],
            [("chr2", 43), ("chrX", 60)],
        )

    def test_nonrepeat_anchor_still_restores_only_one_repeat(self):
        rows = contig(
            ("chr1", 20, True),
            ("chr1", 30, True),
            ("chr1", 60, False),
        )

        kept = preprocess_repeat(rows)

        self.assertEqual(
            [(row[9], row[16]) for row in kept],
            [(30, "chr1"), (60, "0")],
        )


if __name__ == "__main__":
    unittest.main()
