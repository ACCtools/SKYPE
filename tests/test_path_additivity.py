from __future__ import annotations

import sys
import unittest
from collections import Counter
from pathlib import Path

import numpy as np


SKYPE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKYPE_ROOT))

from path_additivity import (  # noqa: E402
    choose_removal,
    find_additive_path_triples,
    find_chained_path_triples,
    find_path_conflicts,
)


class AdditivePathTripleTests(unittest.TestCase):
    def test_repeated_and_shared_nclose_counts_are_added(self) -> None:
        x = (1, 2)
        y = (3, 4)
        usage = [
            Counter({x: 2}),
            Counter({x: 1, y: 1}),
            Counter({x: 3, y: 1}),
        ]

        triples = find_additive_path_triples(
            {0, 1, 2}, usage, np.asarray([3.0, 2.0, 4.0])
        )

        self.assertEqual(len(triples), 1)
        self.assertEqual(
            (triples[0].path_a, triples[0].path_b, triples[0].path_ab),
            (0, 1, 2),
        )
        self.assertEqual(triples[0].total_weight, 9.0)

    def test_plain_set_match_with_wrong_multiplicity_is_rejected(self) -> None:
        x = (1, 2)
        y = (3, 4)
        usage = [
            Counter({x: 2}),
            Counter({y: 1}),
            Counter({x: 1, y: 1}),
        ]

        self.assertEqual(
            find_additive_path_triples({0, 1, 2}, usage, [1.0, 1.0, 1.0]),
            [],
        )

    def test_a_and_b_must_have_different_usage_multisets(self) -> None:
        x = (1, 2)
        usage = [
            Counter({x: 1}),
            Counter({x: 1}),
            Counter({x: 2}),
        ]

        self.assertEqual(
            find_additive_path_triples({0, 1, 2}, usage, [1.0, 1.0, 1.0]),
            [],
        )

    def test_empty_usage_is_not_an_additive_component(self) -> None:
        x = (1, 2)
        usage = [Counter(), Counter({x: 1}), Counter({x: 1})]

        self.assertEqual(
            find_additive_path_triples({0, 1, 2}, usage, [1.0, 1.0, 1.0]),
            [],
        )

    def test_multiple_ab_paths_are_distinct_and_sorted_by_weight_sum(self) -> None:
        x = (1, 2)
        y = (3, 4)
        z = (5, 6)
        usage = [
            Counter({x: 1}),
            Counter({y: 1}),
            Counter({x: 1, y: 1}),
            Counter({x: 1, y: 1}),
            Counter({z: 1}),
            Counter({x: 1, z: 1}),
        ]
        weights = [5.0, 4.0, 3.0, 1.0, 2.0, 0.5]

        triples = find_additive_path_triples(range(6), usage, weights)

        self.assertEqual(
            [(item.path_a, item.path_b, item.path_ab) for item in triples],
            [(0, 4, 5), (0, 1, 3), (0, 1, 2)],
        )


class ChainedPathTripleTests(unittest.TestCase):
    def test_exact_repeated_counts_form_a_bc_ab_chain(self) -> None:
        a = (1, 2)
        b = (3, 4)
        c = (5, 6)
        usage = [
            Counter({a: 2}),
            Counter({b: 2, c: 1}),
            Counter({a: 2, b: 2}),
        ]

        triples = find_chained_path_triples(
            {0, 1, 2}, usage, np.asarray([3.0, 2.0, 4.0])
        )

        self.assertEqual(len(triples), 1)
        self.assertEqual(
            (triples[0].path_a, triples[0].path_bc, triples[0].path_ab),
            (0, 1, 2),
        )
        self.assertEqual(triples[0].total_weight, 9.0)

    def test_wrong_shared_b_multiplicity_is_rejected(self) -> None:
        a = (1, 2)
        b = (3, 4)
        c = (5, 6)
        usage = [
            Counter({a: 1}),
            Counter({b: 1, c: 1}),
            Counter({a: 1, b: 2}),
        ]

        self.assertEqual(
            find_chained_path_triples({0, 1, 2}, usage, [1.0, 1.0, 1.0]),
            [],
        )

    def test_a_overlap_with_bc_is_rejected(self) -> None:
        a = (1, 2)
        b = (3, 4)
        c = (5, 6)
        usage = [
            Counter({a: 1}),
            Counter({a: 1, b: 1, c: 1}),
            Counter({a: 1, b: 1}),
        ]

        self.assertEqual(
            find_chained_path_triples({0, 1, 2}, usage, [1.0, 1.0, 1.0]),
            [],
        )

    def test_b_or_c_cannot_be_empty(self) -> None:
        a = (1, 2)
        b = (3, 4)
        usage = [
            Counter({a: 1}),
            Counter({b: 1}),
            Counter({a: 1, b: 1}),
        ]

        self.assertEqual(
            find_chained_path_triples({0, 1, 2}, usage, [1.0, 1.0, 1.0]),
            [],
        )

    def test_both_conflict_types_share_global_weight_order(self) -> None:
        p = (1, 2)
        q = (3, 4)
        a = (5, 6)
        b = (7, 8)
        c = (9, 10)
        usage = [
            Counter({p: 1}),
            Counter({q: 1}),
            Counter({p: 1, q: 1}),
            Counter({a: 1}),
            Counter({b: 1, c: 1}),
            Counter({a: 1, b: 1}),
        ]
        weights = [4.0, 4.0, 4.0, 1.0, 1.0, 1.0]

        conflicts = find_path_conflicts(range(6), usage, weights)

        self.assertEqual(
            [conflict.conflict_type for conflict in conflicts],
            ['A_BC_AB', 'A_B_AB'],
        )
        self.assertEqual(
            (conflicts[0].path_a, conflicts[0].path_b, conflicts[0].path_ab),
            (3, 4, 5),
        )


class RemovalChoiceTests(unittest.TestCase):
    def test_only_feasible_single_is_removed(self) -> None:
        self.assertEqual(
            choose_removal(
                10, 20,
                a_feasible=True,
                b_feasible=False,
                a_error=5.0,
                b_error=1.0,
            ).removed_columns,
            (10,),
        )

    def test_joint_removal_wins_when_feasible(self) -> None:
        choice = choose_removal(
            10, 20,
            a_feasible=True,
            b_feasible=True,
            both_feasible=True,
            a_error=5.0,
            b_error=1.0,
        )

        self.assertEqual(choice.removed_columns, (10, 20))
        self.assertEqual(choice.reason, "joint_removal_feasible")

    def test_lower_error_single_wins_after_joint_failure(self) -> None:
        choice = choose_removal(
            10, 20,
            a_feasible=True,
            b_feasible=True,
            both_feasible=False,
            a_error=5.0,
            b_error=1.0,
        )

        self.assertEqual(choice.removed_columns, (20,))
        self.assertEqual(choice.reason, "lower_error_after_joint_failure")

    def test_lower_error_single_is_forced_when_neither_is_feasible(self) -> None:
        choice = choose_removal(
            10, 20,
            a_feasible=False,
            b_feasible=False,
            a_error=2.0,
            b_error=3.0,
        )

        self.assertEqual(choice.removed_columns, (10,))
        self.assertEqual(choice.reason, "lower_error_forced")

    def test_equal_error_is_deterministic(self) -> None:
        choice = choose_removal(
            20, 10,
            a_feasible=False,
            b_feasible=False,
            a_error=2.0,
            b_error=2.0,
        )

        self.assertEqual(choice.removed_columns, (10,))


if __name__ == "__main__":
    unittest.main()
