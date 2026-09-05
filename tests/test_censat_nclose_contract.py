from __future__ import annotations

import ast
import logging
import unittest
from collections import Counter
from enum import Enum
from pathlib import Path

from nclose_candidate import NCloseCandidate, apply_nclose_filter


BUILD_GRAPH_PATH = (
    Path(__file__).resolve().parents[1] / "nclose_preprocess.py"
)


CTG_NAM = 0
CTG_DIR = 4
CHR_NAM = 5
CHR_LEN = 6
CHR_STR = 7
CHR_END = 8
CTG_MAPQ = 9
CTG_CENSAT = 18


def load_censat_helpers():
    """Load CenSat contracts without executing stage 01's CLI entry point."""

    definition_names = {
        "CensatPairClass",
        "node_is_censat",
        "classify_censat_pair",
        "_censat_at_chromosome_end",
        "apply_censat_censat_filter",
    }
    tree = ast.parse(BUILD_GRAPH_PATH.read_text(encoding="utf-8"))
    definitions = [
        node
        for node in tree.body
        if isinstance(node, (ast.ClassDef, ast.FunctionDef))
        and node.name in definition_names
    ]
    found_names = {node.name for node in definitions}
    missing_names = definition_names - found_names
    if missing_names:
        raise AssertionError(
            f"Missing CenSat contract definitions: {sorted(missing_names)}"
        )

    namespace = {
        "Counter": Counter,
        "Enum": Enum,
        "logging": logging,
        "NCloseCandidate": NCloseCandidate,
        "apply_nclose_filter": apply_nclose_filter,
        "CTG_NAM": CTG_NAM,
        "CTG_DIR": CTG_DIR,
        "CHR_NAM": CHR_NAM,
        "CHR_LEN": CHR_LEN,
        "CHR_STR": CHR_STR,
        "CHR_END": CHR_END,
        "CTG_MAPQ": CTG_MAPQ,
        "CTG_CENSAT": CTG_CENSAT,
    }
    module = ast.Module(body=definitions, type_ignores=[])
    exec(compile(module, str(BUILD_GRAPH_PATH), "exec"), namespace)
    return namespace


HELPERS = load_censat_helpers()
CensatPairClass = HELPERS["CensatPairClass"]
classify_censat_pair = HELPERS["classify_censat_pair"]
apply_censat_censat_filter = HELPERS["apply_censat_censat_filter"]


def contig_row(
    name: str,
    chrom: str,
    start: int,
    end: int,
    direction: str = "+",
    censat: str = "0",
    mapq: int = 60,
):
    row = ["0"] * 22
    row[CTG_NAM] = name
    row[1] = max(end - start, 1)
    row[2] = 0
    row[3] = max(end - start, 1)
    row[CTG_DIR] = direction
    row[CHR_NAM] = chrom
    row[CHR_LEN] = 10_000_000
    row[CHR_STR] = start
    row[CHR_END] = end
    row[CTG_MAPQ] = mapq
    row[CTG_CENSAT] = censat
    return tuple(row)


class CensatPairClassificationTests(unittest.TestCase):
    def test_none_one_and_both_depend_only_on_endpoint_labels(self):
        contig_data = [
            contig_row("none-a", "chr1", 10, 20, censat="0"),
            contig_row("none-b", "chr2", 10, 20, censat="0"),
            contig_row("one-a", "chr1", 10, 20, censat="r"),
            contig_row("one-b", "chr2", 10, 20, censat="0"),
            contig_row("both-a", "chr1", 10, 20, censat="r"),
            contig_row("both-b", "chr2", 10, 20, censat="rin"),
        ]

        self.assertIs(
            classify_censat_pair(contig_data, (0, 1)),
            CensatPairClass.NONE,
        )
        self.assertIs(
            classify_censat_pair(contig_data, (2, 3)),
            CensatPairClass.ONE,
        )
        self.assertIs(
            classify_censat_pair(contig_data, (4, 5)),
            CensatPairClass.BOTH,
        )

class BothCensatFilterTests(unittest.TestCase):
    def test_both_only_filter_keeps_none_and_one_even_with_low_mapq(self):
        contig_data = [
            contig_row("none-a", "chr1", 1_000, 1_100, censat="0", mapq=1),
            contig_row("none-b", "chr2", 1_000, 1_100, censat="0", mapq=1),
            contig_row("one-a", "chr3", 1_000, 1_100, censat="r", mapq=1),
            contig_row("one-b", "chr4", 1_000, 1_100, censat="0", mapq=1),
            contig_row("both-a", "chr5", 1_000, 1_100, censat="r", mapq=59),
            contig_row("both-b", "chr6", 1_000, 1_100, censat="rin", mapq=60),
        ]
        candidates = [
            NCloseCandidate("none", (0, 1)),
            NCloseCandidate("one", (2, 3)),
            NCloseCandidate("both", (4, 5)),
        ]

        kept, rejected = apply_censat_censat_filter(
            candidates,
            contig_data,
            {f"chr{i}": 10_000_000 for i in range(1, 7)},
            {},
        )

        self.assertEqual([candidate.contig_name for candidate in kept], [
            "none",
            "one",
        ])
        self.assertEqual(
            [(item.candidate.contig_name, item.reason) for item in rejected],
            [("both", "mapq")],
        )

    def test_stage01_keeps_cross_chromosome_opposite_direction_pairs(self):
        contig_data = [
            contig_row("opposite-a", "chr1", 1_000, 1_100, "+", "r"),
            contig_row("opposite-a", "chr2", 2_000, 2_100, "+", "r"),
            contig_row("opposite-b", "chr1", 1_200, 1_300, "-", "r"),
            contig_row("opposite-b", "chr2", 2_200, 2_300, "-", "r"),
            contig_row("same-chrom", "chr3", 1_000, 1_100, "+", "r", 1),
            contig_row("same-chrom", "chr3", 2_000, 2_100, "-", "r", 1),
            contig_row("mapq", "chr4", 1_000, 1_100, "+", "r", 59),
            contig_row("mapq", "chr5", 2_000, 2_100, "+", "r", 60),
            contig_row("terminal", "chr6", 100, 200, "+", "r", 1),
            contig_row("terminal", "chr7", 2_000, 2_100, "+", "r", 1),
        ]
        candidates = [
            NCloseCandidate("opposite-a", (0, 1)),
            NCloseCandidate("opposite-b", (2, 3)),
            NCloseCandidate("same-chrom", (4, 5)),
            NCloseCandidate("mapq", (6, 7)),
            NCloseCandidate("terminal", (8, 9)),
        ]

        kept, rejected = apply_censat_censat_filter(
            candidates,
            contig_data,
            {f"chr{i}": 10_000_000 for i in range(1, 8)},
            {"chr6": [(0, 1_000)]},
        )

        self.assertEqual(kept, candidates[:2])
        self.assertEqual(
            {item.candidate.contig_name: item.reason for item in rejected},
            {
                "same-chrom": "same_chrom_opposite_dir",
                "mapq": "mapq",
                "terminal": "terminal",
            },
        )

if __name__ == "__main__":
    unittest.main()
