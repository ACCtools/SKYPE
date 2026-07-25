import os
import tempfile
import unittest
from collections import Counter

import numpy as np
from scipy.optimize import nnls

from censat_pair_rescue import (
    Alignment,
    allowed_rescue_reactivation_columns,
    canonicalize_pair,
    classify_rescue_cf_swap_failures,
    empty_rescue_artifact,
    eligible_rescue_event_keys,
    load_rescue_artifact,
    plan_rescue_cf_swap,
    prepare_reactivation_warm_start,
    rescue_record,
    save_rescue_artifact,
    select_rescue_candidates_from_groups,
    select_rescue_reactivation_columns,
    stable_rescue_prefix_permutation,
    validate_rescue_path_prefix,
)


def alignment(
    query_name,
    query_length,
    query_start,
    query_end,
    strand,
    chromosome,
    source_order=0,
):
    return Alignment(
        query_name=query_name,
        query_length=query_length,
        query_start=query_start,
        query_end=query_end,
        strand=strand,
        target_name=chromosome,
        target_length=10_000_000,
        target_start=100,
        target_end=200,
        matches=100,
        block_length=100,
        mapq=20,
        source_order=source_order,
    )


def candidate_rows(
    query_name,
    query_length=1_000,
    left_state=("chr1", "+"),
    right_state=("chr2", "-"),
    left_span=(0, 600),
    right_span=(400, 1_000),
):
    left = alignment(
        query_name,
        query_length,
        left_span[0],
        min(left_span[1], query_length),
        left_state[1],
        left_state[0],
        0,
    )
    right = alignment(
        query_name,
        query_length,
        right_span[0],
        min(right_span[1], query_length),
        right_state[1],
        right_state[0],
        1,
    )
    return [left, right], [left, right]


def artifact_record(event_key, chroms):
    return {
        "event_key": event_key,
        "query_name": f"q{event_key[0]}",
        "query_length": 100,
        "canonical_states": (
            (chroms[0], "+"),
            (chroms[1], "+"),
        ),
        "chroms": chroms,
    }


class CensatPairCandidateTests(unittest.TestCase):
    def setUp(self):
        self.censat = {
            f"chr{i}": [(0, 2_000_001)]
            for i in range(1, 20)
        }
        # Raw chr1+ -> chr2- normalizes to chr1+, chr2+.
        self.meta = {
            f"chr{i}": {"dir": False}
            for i in range(1, 20)
        }

    def select(self, aligned, original, meta=None):
        return select_rescue_candidates_from_groups(
            aligned,
            original,
            self.censat,
            self.meta if meta is None else meta,
        )

    def test_original_state_set_must_be_exact(self):
        good_aln, good_original = candidate_rows("good")
        missing_aln, missing_original = candidate_rows("missing")
        extra_aln, extra_original = candidate_rows("extra")
        missing_original = missing_original[:1]
        extra_original = extra_original + [
            alignment("extra", 1_000, 200, 300, "+", "chr3", 2)
        ]

        selected = self.select(
            {
                "good": good_aln,
                "missing": missing_aln,
                "extra": extra_aln,
            },
            {
                "good": good_original,
                "missing": missing_original,
                "extra": extra_original,
            },
        )

        self.assertEqual([row.query_name for row in selected], ["good"])

    def test_bounding_start_and_end_order_are_both_strict(self):
        good_aln, good_original = candidate_rows("good")
        equal_start_aln, equal_start_original = candidate_rows(
            "equal_start",
            left_span=(0, 600),
            right_span=(0, 900),
        )
        equal_end_aln, equal_end_original = candidate_rows(
            "equal_end",
            left_span=(0, 600),
            right_span=(200, 600),
        )

        selected = self.select(
            {
                "good": good_aln,
                "equal_start": equal_start_aln,
                "equal_end": equal_end_aln,
            },
            {
                "good": good_original,
                "equal_start": equal_start_original,
                "equal_end": equal_end_original,
            },
        )

        self.assertEqual([row.query_name for row in selected], ["good"])

    def test_reverse_complements_share_group_and_longest_then_name_wins(self):
        short_aln, short_original = candidate_rows(
            "short",
            query_length=1_000,
            left_state=("chr1", "+"),
            right_state=("chr2", "-"),
        )
        long_z_aln, long_z_original = candidate_rows(
            "z_long",
            query_length=2_000,
            left_state=("chr2", "+"),
            right_state=("chr1", "-"),
            left_span=(0, 1_200),
            right_span=(800, 2_000),
        )
        long_a_aln, long_a_original = candidate_rows(
            "a_long",
            query_length=2_000,
            left_state=("chr1", "+"),
            right_state=("chr2", "-"),
            left_span=(0, 1_200),
            right_span=(800, 2_000),
        )

        selected = self.select(
            {
                "short": short_aln,
                "z_long": long_z_aln,
                "a_long": long_a_aln,
            },
            {
                "short": short_original,
                "z_long": long_z_original,
                "a_long": long_a_original,
            },
        )

        self.assertEqual(len(selected), 1)
        self.assertEqual(selected[0].query_name, "a_long")
        self.assertEqual(
            selected[0].canonical_states,
            canonicalize_pair((("chr1", "+"), ("chr2", "-"))),
        )

    def test_both_normalized_endpoint_directions_must_match(self):
        aligned, original = candidate_rows("query")
        mismatch_meta = dict(self.meta)
        mismatch_meta["chr2"] = {"dir": True}

        matched = self.select({"query": aligned}, {"query": original})
        self.assertEqual([row.query_name for row in matched], ["query"])
        self.assertEqual(
            self.select(
                {"query": aligned},
                {"query": original},
                mismatch_meta,
            ),
            [],
        )


class CensatPairArtifactAndPathTests(unittest.TestCase):
    def test_missing_and_round_trip_artifact(self):
        with tempfile.TemporaryDirectory() as prefix:
            self.assertEqual(load_rescue_artifact(prefix), empty_rescue_artifact())
            record = artifact_record((4, 5), ("chr12", "chr16"))
            save_rescue_artifact(prefix, [record])
            loaded = load_rescue_artifact(prefix)

        self.assertEqual(loaded["version"], 1)
        self.assertEqual(loaded["records"], [record])

    def test_path_and_depth_use_the_same_stable_prefix_permutation(self):
        paths = ["depth10", "depth9", "depth8", "depth7"]
        depths = [10, 9, 8, 7]
        permutation = stable_rescue_prefix_permutation(
            [False, True, True, False]
        )

        self.assertEqual(permutation, [1, 2, 0, 3])
        self.assertEqual([paths[i] for i in permutation], [
            "depth9", "depth8", "depth10", "depth7"
        ])
        self.assertEqual([depths[i] for i in permutation], [9, 8, 10, 7])

    def test_rescue_carriers_form_exact_contiguous_path_prefix(self):
        records = [artifact_record((10, 11), ("chr12", "chr16"))]
        good_usage = [
            Counter({(10, 11): 1}),
            Counter({(10, 11): 1}),
            Counter(),
            Counter(),
        ]
        bad_usage = [
            Counter({(10, 11): 1}),
            Counter(),
            Counter({(10, 11): 1}),
            Counter(),
        ]

        self.assertEqual(
            validate_rescue_path_prefix(good_usage, records, 4),
            2,
        )
        with self.assertRaisesRegex(ValueError, "contiguous prefix"):
            validate_rescue_path_prefix(bad_usage, records, 4)

    def test_only_pair_with_both_retained_chromosomes_reactivates(self):
        records = [
            artifact_record((10, 11), ("chr12", "chr16")),
            artifact_record((20, 21), ("chr6", "chr16")),
            artifact_record((30, 31), ("chr1", "chr5")),
        ]
        usage = [
            Counter({(10, 11): 1}),
            Counter({(20, 21): 1}),
            Counter({(30, 31): 1}),
            Counter({(10, 11): 1, (20, 21): 1}),
        ]

        eligible = eligible_rescue_event_keys(
            records, {"chr12", "chr16", "chr8"}
        )
        columns, selected_events = select_rescue_reactivation_columns(
            records,
            usage,
            {"chr12", "chr16", "chr8"},
            allowed_columns={0, 1, 2, 3},
        )

        self.assertEqual(eligible, {(10, 11)})
        self.assertEqual(selected_events, {(10, 11)})
        # A multi-rescue path is restored when any one event is eligible.
        self.assertEqual(columns, [0, 3])

    def test_unrelated_filters_still_block_a_shared_rescue_path(self):
        rescue_key = (10, 11)
        unrelated_key = (50, 51)
        usage = [
            Counter({rescue_key: 1}),
            Counter({rescue_key: 1, unrelated_key: 1}),
            Counter({rescue_key: 1}),
        ]

        allowed = allowed_rescue_reactivation_columns(
            deferred_columns={0, 1, 2},
            path_usage=usage,
            rescue_keys={rescue_key},
            blocked_columns={2},
            blocked_events={rescue_key, unrelated_key},
        )

        # The rescue event's own pre-CF result is ignored. An unrelated event
        # or reciprocal-path rejection remains authoritative.
        self.assertEqual(allowed, {0})

    def test_noncontiguous_global_indices_keep_rows_and_weights_aligned(self):
        current_indices = [2, 7]
        current_weights = np.array([3.0, 4.0])
        new_indices, warm_start, full_warm = prepare_reactivation_warm_start(
            current_weights,
            current_indices,
            reactivated_global_indices={1, 9},
            feature_count=10,
        )

        self.assertEqual(new_indices, [1, 2, 7, 9])
        np.testing.assert_array_equal(warm_start, [0.0, 3.0, 4.0, 0.0])
        np.testing.assert_array_equal(
            full_warm,
            [0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0],
        )

        # matrix.h5 is feature-major. Re-reading these global rows in sorted
        # order and transposing must preserve the global column identity.
        feature_major = np.eye(10, dtype=float)
        target = np.zeros(10)
        target[[1, 2, 7, 9]] = [5.0, 3.0, 4.0, 6.0]
        fitted, _ = nnls(feature_major[new_indices].T, target)
        full_fitted = np.zeros(10)
        full_fitted[new_indices] = fitted
        np.testing.assert_allclose(
            full_fitted[[1, 2, 7, 9]],
            [5.0, 3.0, 4.0, 6.0],
        )
        np.testing.assert_array_equal(
            np.flatnonzero(full_fitted),
            [1, 2, 7, 9],
        )

    def test_multiple_rescue_events_are_planned_as_one_cf_swap(self):
        records = [
            artifact_record((10, 11), ("chr12", "chr16")),
            artifact_record((20, 21), ("chr6", "chr7")),
        ]
        usage = [
            Counter({(10, 11): 1}),
            Counter({(20, 21): 1}),
            Counter({(10, 11): 1, (20, 21): 1}),
        ]

        plan = plan_rescue_cf_swap(
            records,
            usage,
            eligible_event_keys={(10, 11), (20, 21)},
            reactivated_columns={0, 1, 2},
            cent_fragment_col2chrom={
                20: "chr6",
                21: "chr7",
                22: "chr8",
                23: "chr12",
                24: "chr16",
            },
            active_columns={0, 1, 2, 20, 21, 22, 23, 24},
        )

        self.assertEqual(
            plan.event_keys,
            frozenset({(10, 11), (20, 21)}),
        )
        self.assertEqual(
            plan.active_event_chroms,
            frozenset({"chr6", "chr7", "chr12", "chr16"}),
        )
        self.assertEqual(plan.disabled_cf_columns, (20, 21, 23, 24))

    def test_swap_does_not_disable_cf_for_event_without_a_carrier(self):
        records = [
            artifact_record((10, 11), ("chr12", "chr16")),
            artifact_record((20, 21), ("chr6", "chr7")),
        ]
        usage = [Counter({(10, 11): 1})]

        plan = plan_rescue_cf_swap(
            records,
            usage,
            eligible_event_keys={(10, 11), (20, 21)},
            reactivated_columns={0},
            cent_fragment_col2chrom={
                20: "chr6",
                21: "chr7",
                23: "chr12",
                24: "chr16",
            },
            active_columns={0, 20, 21, 23, 24},
        )

        self.assertEqual(plan.event_keys, frozenset({(10, 11)}))
        self.assertEqual(
            plan.active_event_chroms,
            frozenset({"chr12", "chr16"}),
        )
        self.assertEqual(plan.disabled_cf_columns, (23, 24))

    def test_failed_event_drops_shared_paths_and_restores_all_its_cf(self):
        event_6_9 = (10, 11)
        event_6_11 = (20, 21)
        event_2_3 = (30, 31)
        event_9_12 = (40, 41)
        records = [
            artifact_record(event_6_9, ("chr6", "chr9")),
            artifact_record(event_6_11, ("chr6", "chr11")),
            artifact_record(event_2_3, ("chr2", "chr3")),
            artifact_record(event_9_12, ("chr9", "chr12")),
        ]
        usage = [
            Counter({event_6_9: 1}),
            Counter({event_6_11: 1}),
            Counter({event_2_3: 1}),
            Counter({event_9_12: 1}),
            Counter({event_6_9: 1, event_2_3: 1}),
        ]
        cf_col2chrom = {
            20: "chr2",
            21: "chr3",
            22: "chr6",
            23: "chr9",
            24: "chr11",
            25: "chr12",
        }
        active_columns = set(range(5)) | set(cf_col2chrom)
        initial_plan = plan_rescue_cf_swap(
            records,
            usage,
            eligible_event_keys={
                event_6_9,
                event_6_11,
                event_2_3,
                event_9_12,
            },
            reactivated_columns=range(5),
            cent_fragment_col2chrom=cf_col2chrom,
            active_columns=active_columns,
        )

        rejected, unrelated = classify_rescue_cf_swap_failures(
            records,
            initial_plan.event_keys,
            failed_chromosomes={"chr6"},
        )
        retry_plan = plan_rescue_cf_swap(
            records,
            usage,
            eligible_event_keys=initial_plan.event_keys,
            reactivated_columns=range(5),
            cent_fragment_col2chrom=cf_col2chrom,
            active_columns=active_columns,
            rejected_event_keys=rejected,
        )

        self.assertFalse(unrelated)
        self.assertEqual(rejected, {event_6_9, event_6_11})
        # Both chr6 paths, including the one shared with (chr2, chr3), drop.
        self.assertEqual(retry_plan.rescue_columns, (2, 3))
        self.assertEqual(
            retry_plan.event_keys,
            frozenset({event_2_3, event_9_12}),
        )
        # chr9 is restored even though the retained (chr9, chr12) uses it.
        self.assertEqual(
            retry_plan.restored_cf_chroms,
            frozenset({"chr6", "chr9", "chr11"}),
        )
        self.assertEqual(
            retry_plan.disabled_cf_chroms,
            frozenset({"chr2", "chr3", "chr12"}),
        )
        self.assertEqual(retry_plan.disabled_cf_columns, (20, 21, 25))

    def test_failure_on_unrelated_chromosome_requests_full_cancel(self):
        records = [
            artifact_record((10, 11), ("chr6", "chr9")),
            artifact_record((20, 21), ("chr2", "chr3")),
        ]

        rejected, unrelated = classify_rescue_cf_swap_failures(
            records,
            active_event_keys={(10, 11), (20, 21)},
            failed_chromosomes={"chr6", "chr20"},
        )

        self.assertEqual(rejected, {(10, 11)})
        self.assertEqual(unrelated, {"chr20"})

    def test_related_failures_can_be_pruned_over_multiple_rounds(self):
        event_6_9 = (10, 11)
        event_2_3 = (20, 21)
        event_4_5 = (30, 31)
        records = [
            artifact_record(event_6_9, ("chr6", "chr9")),
            artifact_record(event_2_3, ("chr2", "chr3")),
            artifact_record(event_4_5, ("chr4", "chr5")),
        ]
        usage = [
            Counter({event_6_9: 1}),
            Counter({event_2_3: 1}),
            Counter({event_4_5: 1}),
        ]
        cf_col2chrom = {
            20: "chr2",
            21: "chr3",
            22: "chr4",
            23: "chr5",
            24: "chr6",
            25: "chr9",
        }
        candidates = {event_6_9, event_2_3, event_4_5}
        rejected = set()

        first_rejected, first_unrelated = classify_rescue_cf_swap_failures(
            records,
            candidates,
            failed_chromosomes={"chr6"},
        )
        self.assertFalse(first_unrelated)
        rejected.update(first_rejected)
        second_plan = plan_rescue_cf_swap(
            records,
            usage,
            eligible_event_keys=candidates,
            reactivated_columns=range(3),
            cent_fragment_col2chrom=cf_col2chrom,
            active_columns=set(cf_col2chrom),
            rejected_event_keys=rejected,
        )

        second_rejected, second_unrelated = classify_rescue_cf_swap_failures(
            records,
            second_plan.event_keys,
            failed_chromosomes={"chr2"},
        )
        self.assertFalse(second_unrelated)
        rejected.update(second_rejected)
        third_plan = plan_rescue_cf_swap(
            records,
            usage,
            eligible_event_keys=candidates,
            reactivated_columns=range(3),
            cent_fragment_col2chrom=cf_col2chrom,
            active_columns=set(cf_col2chrom),
            rejected_event_keys=rejected,
        )

        self.assertEqual(rejected, {event_6_9, event_2_3})
        self.assertEqual(third_plan.event_keys, frozenset({event_4_5}))
        self.assertEqual(third_plan.rescue_columns, (2,))
        self.assertEqual(
            third_plan.restored_cf_chroms,
            frozenset({"chr2", "chr3", "chr6", "chr9"}),
        )
        self.assertEqual(
            third_plan.disabled_cf_chroms,
            frozenset({"chr4", "chr5"}),
        )

    def test_reactivation_warm_start_can_disable_cf_columns_atomically(self):
        new_indices, warm_start, full_warm = prepare_reactivation_warm_start(
            current_weights=np.array([2.0, 12.0, 16.0]),
            current_global_indices=[2, 12, 16],
            reactivated_global_indices={0, 1},
            feature_count=20,
            deactivated_global_indices={12, 16},
        )

        self.assertEqual(new_indices, [0, 1, 2])
        np.testing.assert_array_equal(warm_start, [0.0, 0.0, 2.0])
        np.testing.assert_array_equal(
            np.flatnonzero(full_warm),
            [2, 12, 16],
        )


if __name__ == "__main__":
    unittest.main()
