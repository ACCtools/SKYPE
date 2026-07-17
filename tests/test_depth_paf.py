from __future__ import annotations

import ast
import unittest
from pathlib import Path


RUN_DEPTH_PATH = Path(__file__).resolve().parents[1] / "21_run_depth.py"


def load_depth_paf_formatter():
    """Load the pure helper without executing the stage-21 pipeline script."""
    tree = ast.parse(RUN_DEPTH_PATH.read_text(encoding="utf-8"))
    function = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name == "format_nonzero_depth_paf_row"
    )
    module = ast.Module(body=[function], type_ignores=[])
    namespace = {}
    exec(compile(module, str(RUN_DEPTH_PATH), "exec"), namespace)
    return namespace["format_nonzero_depth_paf_row"]


format_nonzero_depth_paf_row = load_depth_paf_formatter()


class DepthPafOutputTests(unittest.TestCase):
    def make_row(self) -> list:
        return [
            "contig",
            1000,
            10,
            20,
            "+",
            "chr1",
            100_000,
            100,
            110,
            10,
            10,
            60,
            "cs:Z::10",
        ]

    def test_formats_positive_length_alignment(self) -> None:
        output = format_nonzero_depth_paf_row(self.make_row(), "10M")

        self.assertIsNotNone(output)
        self.assertTrue(output.endswith("\tcg:Z:10M"))

    def test_omits_h1437_zero_length_overlap_result(self) -> None:
        row = self.make_row()
        row[2:4] = [1416, 1416]
        row[7:9] = [72_179_591, 72_179_591]
        row[9:11] = [0, 0]
        row[-1] = "cs:Z:"

        self.assertIsNone(format_nonzero_depth_paf_row(row, ""))

    def test_omits_zero_length_in_any_paf_length_field(self) -> None:
        for field_pair in ((2, 3), (7, 8)):
            with self.subTest(field_pair=field_pair):
                row = self.make_row()
                row[field_pair[1]] = row[field_pair[0]]
                self.assertIsNone(format_nonzero_depth_paf_row(row, "10M"))

        row = self.make_row()
        row[10] = 0
        self.assertIsNone(format_nonzero_depth_paf_row(row, "10M"))

    def test_omits_empty_cigar_even_with_positive_coordinates(self) -> None:
        self.assertIsNone(format_nonzero_depth_paf_row(self.make_row(), ""))


if __name__ == "__main__":
    unittest.main()
