"""Reject zero-span NClose bridges before graph enumeration."""
import itertools
import unittest

import breakend_graph as bg
import nclose_preprocess as pre
from path_geometry import assert_nclose_reference_connections


def node(owner, start, end, strand, chrom="chr7"):
    return [owner, 1000, 0, end-start, strand, chrom, 200_000_000, start, end]


def stage10_graph(rows, pairs):
    _, reverse = bg.chr_correlation_maker(rows)
    return bg.initialize_bnd_graph(rows, pairs, {}, len(rows), reverse)


class ReferenceConnectionAssertionTests(unittest.TestCase):
    BUILDERS = (pre.initialize_inversion_only_graph, stage10_graph)

    def test_both_builders_reject_real_zero_span_circuits(self):
        cases = {
            "HCC1937.AMP.39": [
                node("utg031127l", 145446834, 145467425, "+"),
                node("utg031127l", 145522733, 145544272, "-"),
                node("utg074142l", 145544272, 145551108, "-"),
                node("utg074142l", 145467425, 145482807, "+"),
            ],
            "HCC1954.AMP.570": [
                node("utg044948l", 78251482, 78290651, "+"),
                node("utg044948l", 76505222, 76531313, "-"),
                node("utg072358l", 76531313, 76546415, "-"),
                node("utg072358l", 78290651, 78315807, "+"),
            ],
        }
        for label, rows in cases.items():
            pairs = {rows[0][0]: [(0, 1)], rows[2][0]: [(2, 3)]}
            for builder in self.BUILDERS:
                with self.subTest(case=label, builder=builder.__name__):
                    graph = builder(rows, pairs)
                    self.assertNotIn([1, 2], graph.get((1, 1), []))
                    self.assertNotIn([1, 0], graph.get((1, 3), []))
                    self.assertIn([1, 1], graph[1, 0])
                    self.assertIn([1, 3], graph[1, 2])
                    # The independent assertion still catches a bad edge if
                    # another producer introduces one in the future.
                    with self.assertRaisesRegex(AssertionError, "reference_span_bp=0") as caught:
                        assert_nclose_reference_connections(
                            rows, {(1, 1): [(1, 2)]}, context=builder.__name__,
                        )
                    for detail in (rows[0][0], rows[2][0], "chr7:", "PAF=", "walked="):
                        self.assertIn(detail, str(caught.exception))
            self.assertEqual(pre.find_ecdna_circuits(rows, pairs), [])

    def test_normal_touch_gap_and_positive_overlap_survive_both_builders(self):
        for first, last in [((1000, 1100), (900, 1000)),
                            ((1000, 1100), (800, 900)),
                            ((950, 1050), (1000, 1100))]:
            rows = [node("a", 700, 750, "+"), node("a", *first, "-"),
                    node("b", *last, "-"), node("b", 1300, 1350, "+")]
            for builder in self.BUILDERS:
                with self.subTest(first=first, last=last, builder=builder.__name__):
                    graph = builder(rows, {"a": [(0, 1)], "b": [(2, 3)]})
                    self.assertIn([1, 2], graph[1, 1])
                    self.assertIn([0, 1], graph[0, 2])

    def test_ecDNA_rejects_zero_bridge_even_if_other_components_have_positive_length(self):
        # One invalid bridge and one ordinary bridge: checking only total > 0
        # would let this circuit pass.
        rows = [node("a", 100, 200, "+"), node("a", 900, 1000, "-"),
                node("b", 1000, 1100, "-"), node("b", 50, 100, "+")]
        self.assertEqual(pre.circuit_length_calculator((0, 1, 2, 3), rows), float("inf"))
        self.assertEqual(pre.find_ecdna_circuits(rows, {"a": [(0, 1)], "b": [(2, 3)]}), [])

    def test_boundary_rules_for_all_NClose_ports_and_both_walked_strands(self):
        # Direction and endpoint role are independent of the input PAF strand.
        # Include both owner orders to exercise both reciprocal case branches.
        cases = [
            ("normal_touch", (900, 1000), (1000, 1100), True),
            ("inward_touch", (1000, 1100), (900, 1000), False),
            ("normal_gap", (800, 900), (1000, 1100), True),
            ("reverse_gap", (1100, 1200), (900, 1000), False),
            ("positive_overlap", (1000, 1100), (950, 1050), True),
        ]
        for a_side, b_side, walked, swap in itertools.product((0, 1), (0, 1), ("+", "-"), (False, True)):
            a_type, b_type = a_side, 1-b_side
            other = "-" if walked == "+" else "+"
            for label, first, last, expected in cases:
                if walked == "-":
                    first, last = last, first
                rows = [node("a", 10000, 10100, "+"), node("a", 11000, 11100, "-"),
                        node("b", 20000, 20100, "-"), node("b", 21000, 21100, "+")]
                rows[a_side] = node("a", *first, walked if a_type else other)
                rows[2+b_side] = node("b", *last, walked if b_type else other)
                pairs = {"a": [(0, 1)], "b": [(2, 3)]}
                if swap:
                    pairs = dict(reversed(list(pairs.items())))
                for builder in self.BUILDERS:
                    with self.subTest(case=label, a=a_side, b=b_side, walked=walked, swap=swap,
                                      builder=builder.__name__):
                        graph = builder(rows, pairs)
                        self.assertEqual([b_type, 2+b_side] in graph.get((a_type, a_side), []), expected)
                        self.assertEqual([1-a_type, a_side] in graph.get((1-b_type, 2+b_side), []), expected)

    def test_assertion_uses_walked_strand_and_is_reverse_symmetric(self):
        for a_strand, b_strand, walked in itertools.product(("+", "-"), repeat=3):
            intervals = [(900, 1000), (1000, 1100)]
            if walked == "-":
                intervals.reverse()
            rows = [node("a", *intervals[0], a_strand), node("b", *intervals[1], b_strand)]
            first, last = (int(a_strand == walked), 0), (int(b_strand == walked), 1)
            with self.subTest(a=a_strand, b=b_strand, walked=walked):
                assert_nclose_reference_connections(rows, {first: [last]}, context="test")
                assert_nclose_reference_connections(
                    rows, {(1-last[0], 1): [(1-first[0], 0)]}, context="reverse",
                )

    def test_rejects_negative_span_and_incompatible_chromosome_or_strand(self):
        for last, expected in [
            (node("b", 1100, 1200, "-"), "reference_span_bp=-100"),
            (node("b", 800, 900, "-", chrom="chr8"), "same chromosome"),
            (node("b", 800, 900, "+"), "same chromosome and walked strand"),
        ]:
            with self.subTest(last=last):
                with self.assertRaisesRegex(AssertionError, expected):
                    assert_nclose_reference_connections(
                        [node("a", 900, 1000, "-"), last], {(1, 0): [(1, 1)]}, context="test",
                    )

    def test_internal_nclose_jumps_and_named_telomeres_are_exempt(self):
        rows = [node("a", 1000, 1100, "+"), node("a", 900, 1000, "-", chrom="chr8")]
        assert_nclose_reference_connections(
            rows, {"chr7f": [(1, 0)], (1, 0): [(1, 1)], (1, 1): ["chr8b"]}, context="test",
        )


if __name__ == "__main__":
    unittest.main()
