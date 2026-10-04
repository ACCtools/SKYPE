"""Depth conservation at reference connections that meet at one boundary."""
from __future__ import annotations

import ast
from collections import Counter
import copy
from itertools import product
from pathlib import Path
import re
import unittest

from path_geometry import walked_strand


SOURCE = Path(__file__).resolve().parents[1] / "21_run_depth.py"


def load_renderer(rows):
    """Load the real renderer and its pure helpers without executing its CLI."""
    names = {
        "node_dir_transform", "form_virtual_contig", "form_normal_contig",
        "form_adjusted_contig", "cs_to_cigar", "compute_aln_stats", "parse_cs",
        "construct_cs", "trim_cs_right_extended", "trim_cs_left_extended",
        "adjust_paf_overlap", "format_nonzero_depth_paf_row",
        "virtual_bridge_follows_strand", "process_raw_contig_list", "distance_checker",
    }
    tree = ast.parse(SOURCE.read_text())
    functions = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    constants = {"CTG_NAM": 0, "CTG_LEN": 1, "CTG_STR": 2, "CTG_END": 3,
                 "CTG_DIR": 4, "CHR_NAM": 5, "CHR_LEN": 6, "CHR_STR": 7,
                 "CHR_END": 8, "CTG_GLOBALIDX": 21, "DIR_FOR": 1, "DIR_BAK": 0,
                 "VIRTUAL_CONTIG_PREFIX": "virtual_contig"}
    namespace = dict(constants, re=re, copy=copy, walked_strand=walked_strand,
                     paf_file=[rows])
    nodes = []
    for index, row in enumerate(rows):
        # Stage-01 rows retain a pointer to the original PAF carrying its CS.
        nodes.append(row[:9] + [row[11], 1, index, index, "0", "0", "0",
                                "0", "0", "0", row[4], row[5], f"0.{index}"])
    namespace["contig_data"] = nodes
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(SOURCE), "exec"), namespace)
    return namespace["process_raw_contig_list"]


def row(name, start, end, strand="+", cs=None):
    cs = cs if cs is not None else f":{end-start}"
    query_span = sum(int(n) for n in re.findall(r":(\d+)", cs))
    query_span += sum(len(s) for s in re.findall(r"\+([a-z]+)", cs))
    query_span += len(re.findall(r"\*[a-z]{2}", cs))
    # The fixtures use exact matches plus insertions; reference consumption is
    # prescribed independently by start/end, while query length includes I.
    return [name, query_span+1000, 500, 500+query_span, strand, "chr1", 100000,
            start, end, end-start, query_span, 60, "tp:A:P", "cs:Z:"+cs]


def traversal(strand, walked):
    return int(strand == walked)


def reference_coverage(lines):
    coverage = Counter()
    for line in lines:
        fields = line.split("\t")
        coverage.update(range(int(fields[7]), int(fields[8])))
    return coverage


class TouchingReferenceDepthTests(unittest.TestCase):
    def render(self, rows, path):
        original = copy.deepcopy(rows)
        result = load_renderer(rows)(path)
        self.assertEqual(rows, original, "Depth rendering must not edit source alignments")
        return result

    def test_zero_contact_and_ordinary_contact_for_all_strands(self):
        # In the inward-facing order, entry and exit are both 11000: the
        # reference walk has zero length. The opposite order spans 2000 bp.
        for first_strand, last_strand, walked, reverse_order in product(
                ["+", "-"], ["+", "-"], ["+", "-"], [False, True]):
            with self.subTest(first=first_strand, last=last_strand,
                              walked=walked, reverse_order=reverse_order):
                intervals = [(10000, 11000), (11000, 12000)]
                if reverse_order:
                    intervals.reverse()
                rows = [row("first", *intervals[0], first_strand),
                        row("last", *intervals[1], last_strand)]
                path = [(traversal(first_strand, walked), 0),
                        (traversal(last_strand, walked), 1)]
                lines, skipped = self.render(rows, path)
                low, high = ((intervals[0][0], intervals[1][1]) if walked == "+"
                             else (intervals[1][0], intervals[0][1]))
                expected = Counter(range(low, high))
                self.assertEqual(reference_coverage(lines), expected)
                self.assertEqual(skipped, 2 if low == high else 0)
                reverse = [(1-d, i) for d, i in reversed(path)]
                reverse_lines, _ = self.render(rows, reverse)
                self.assertEqual(reference_coverage(reverse_lines), expected)

    def test_contact_as_limit_of_one_base_overlap(self):
        # A 1-bp overlapping inward walk contains one base, and touching
        # contains none. Do not jump discontinuously to 2000 bp at equality.
        for overlap in [0, 1, 2]:
            with self.subTest(overlap=overlap):
                rows = [row("first", 10000, 11000),
                        row("last", 11000-overlap, 12000-overlap)]
                lines, _ = self.render(rows, [(0, 0), (0, 1)])
                self.assertEqual(reference_coverage(lines), Counter(range(11000-overlap, 11000)))

    def test_one_base_forward_gap_is_filled_once(self):
        rows = [row("first", 10000, 11000), row("last", 11001, 12001)]
        lines, _ = self.render(rows, [(1, 0), (1, 1)])
        self.assertEqual(reference_coverage(lines), Counter(range(10000, 12001)))

    def test_zero_contact_embedded_in_a_longer_walk(self):
        rows = [row("before", 12000, 13000), row("first", 10000, 11000),
                row("last", 11000, 12000), row("after", 9000, 10000)]
        lines, skipped = self.render(rows, [(0, i) for i in range(4)])
        self.assertEqual(reference_coverage(lines), Counter(range(9000, 13000)))
        self.assertEqual(skipped, 2)

    def test_terminal_insertions_do_not_give_zero_reference_walk_depth(self):
        for first_strand, last_strand in product(["+", "-"], repeat=2):
            with self.subTest(first=first_strand, last=last_strand):
                rows = [row("first", 10000, 11000, first_strand, "+aa:1000+tt"),
                        row("last", 11000, 12000, last_strand, "+gg:1000+cc")]
                lines, skipped = self.render(rows, [(traversal(first_strand, "-"), 0),
                                                    (traversal(last_strand, "-"), 1)])
                self.assertEqual(lines, [])
                self.assertEqual(skipped, 2)

    def test_ordinary_touch_preserves_original_cs_and_query_intervals(self):
        rows = [row("first", 10000, 11000, cs=":1000+aaa"),
                row("last", 11000, 12000, cs="+ttt:1000")]
        lines, skipped = self.render(rows, [(1, 0), (1, 1)])
        self.assertEqual(skipped, 0)
        for original, line in zip(rows, lines, strict=True):
            self.assertEqual(line.split("\t")[:len(original)], list(map(str, original)))

    def test_same_owner_assembly_connection_is_not_a_reference_bridge(self):
        # Source continuity can encode a variant even when its reference
        # pieces appear in the opposite order. It must not be erased here.
        rows = [row("one_observed_source", 10000, 11000),
                row("one_observed_source", 11000, 12000)]
        lines, skipped = self.render(rows, [(0, 0), (0, 1)])
        self.assertEqual(reference_coverage(lines), Counter(range(10000, 12000)))
        self.assertEqual(skipped, 0)


if __name__ == "__main__":
    unittest.main()
