from __future__ import annotations

import ast
import logging
import os
import pickle
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock

from nclose_candidate import (
    NCloseCandidate,
    apply_nclose_filter,
    candidates_to_legacy,
    endpoint_metadata,
    iter_legacy_candidates,
)


BUILD_GRAPH_PATH = (
    Path(__file__).resolve().parents[1] / "nclose_preprocess.py"
)


def ordered_direct_calls(function_name):
    tree = ast.parse(BUILD_GRAPH_PATH.read_text(encoding="utf-8"))
    function = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == function_name
    )
    calls = sorted(
        (node for node in ast.walk(function) if isinstance(node, ast.Call)),
        key=lambda node: (node.lineno, node.col_offset),
    )
    return [
        (call.func.id, call.lineno)
        for call in calls
        if isinstance(call.func, ast.Name)
    ]


def load_stage01_function(function_name, **extra_globals):
    tree = ast.parse(BUILD_GRAPH_PATH.read_text(encoding="utf-8"))
    function = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == function_name
    )
    namespace = {
        "apply_nclose_filter": apply_nclose_filter,
        "logging": logging,
        **extra_globals,
    }
    exec(
        compile(
            ast.Module(body=[function], type_ignores=[]),
            str(BUILD_GRAPH_PATH),
            "exec",
        ),
        namespace,
    )
    return namespace[function_name]


def node(repeat_chrom="0", repeat_case="0", censat="0"):
    row = ["0"] * 19
    row[16] = repeat_chrom
    row[17] = repeat_case
    row[18] = censat
    return row


class NCloseCandidateTests(unittest.TestCase):
    def test_path_pair_is_distinct_from_sorted_event_key(self):
        signature = (("chr4", 100, 200), ("chr9", 300, 400))
        reverse = NCloseCandidate(
            "ctg",
            (9, 4),
            signatures={"compression_layout": signature},
        )
        forward = NCloseCandidate(
            "ctg",
            (4, 9),
            signatures={"compression_layout": signature},
        )

        self.assertEqual(reverse.path_pair, (9, 4))
        self.assertEqual(reverse.identity, ("ctg", (9, 4)))
        self.assertEqual(reverse.event_key, (4, 9))
        self.assertEqual(forward.event_key, reverse.event_key)
        self.assertNotEqual(forward.path_pair, reverse.path_pair)
        self.assertEqual(
            forward.signatures["compression_layout"],
            reverse.signatures["compression_layout"],
        )

    def test_invalid_path_pair_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "Invalid NClose path pair"):
            NCloseCandidate("ctg", (1,))

    def test_legacy_round_trip_preserves_all_orders(self):
        legacy = {
            "ctg_b": [(9, 4), (2, 8)],
            "combined_name_is_not_implicit_provenance": [(7, 1)],
        }
        candidates = list(
            iter_legacy_candidates(
                legacy,
                origin_by_identity={("ctg_b", (2, 8)): "paf"},
            )
        )

        self.assertEqual(
            [candidate.identity for candidate in candidates],
            [
                ("ctg_b", (9, 4)),
                ("ctg_b", (2, 8)),
                ("combined_name_is_not_implicit_provenance", (7, 1)),
            ],
        )
        self.assertEqual(candidates[0].origin, "legacy_unknown")
        self.assertEqual(candidates[1].origin, "paf")
        self.assertEqual(candidates[2].origin, "legacy_unknown")

        round_trip = candidates_to_legacy(candidates)
        self.assertIsInstance(round_trip, defaultdict)
        self.assertEqual(list(round_trip.items()), list(legacy.items()))

    def test_legacy_adapter_preserves_duplicate_occurrences_and_pair_order(self):
        candidates = [
            NCloseCandidate(
                "owner_b",
                (9, 4),
                provenance={"source": "primary"},
            ),
            NCloseCandidate(
                "owner_b",
                (9, 4),
                provenance={"source": "secondary"},
            ),
            NCloseCandidate("owner_b", (4, 9)),
            NCloseCandidate("owner_a", (7, 1)),
        ]

        legacy = candidates_to_legacy(candidates)

        self.assertEqual(
            list(legacy.items()),
            [
                ("owner_b", [(9, 4), (9, 4), (4, 9)]),
                ("owner_a", [(7, 1)]),
            ],
        )
        self.assertEqual(
            [candidate.identity for candidate in iter_legacy_candidates(legacy)],
            [candidate.identity for candidate in candidates],
        )

    def test_endpoint_metadata_is_live_and_in_path_order(self):
        contig_data = [
            node("chr1", "r", "0"),
            node("chr2", "rin", "rin"),
        ]
        candidate = NCloseCandidate("ctg", (1, 0))

        self.assertEqual(
            endpoint_metadata(candidate, contig_data),
            (
                {
                    "node_idx": 1,
                    "repeat_chrom": "chr2",
                    "repeat_case": "rin",
                    "censat": "rin",
                },
                {
                    "node_idx": 0,
                    "repeat_chrom": "chr1",
                    "repeat_case": "r",
                    "censat": "0",
                },
            ),
        )

        contig_data[1][18] = "r"
        self.assertEqual(endpoint_metadata(candidate, contig_data)[0]["censat"], "r")

    def test_synthetic_origin_and_parent_provenance_are_explicit(self):
        candidate = NCloseCandidate(
            "combined_ctg",
            (20, 21),
            origin="synthetic_combined_test",
            synthetic=True,
            provenance={"parent_event_keys": ((3, 8),), "type2_query": "utg7"},
        )

        self.assertTrue(candidate.synthetic)
        self.assertEqual(candidate.origin, "synthetic_combined_test")
        self.assertEqual(candidate.provenance["parent_event_keys"], ((3, 8),))

    def test_append_and_remove_are_stable_for_duplicates_and_synthetic_provenance(self):
        primary = NCloseCandidate(
            "shared-owner",
            (8, 3),
            origin="paf",
            provenance={"source": "primary"},
        )
        secondary = NCloseCandidate(
            "shared-owner",
            (8, 3),
            origin="paf",
            provenance={"source": "secondary"},
        )
        combined = NCloseCandidate(
            "combined-owner",
            (20, 21),
            origin="synthetic_combined_test",
            synthetic=True,
            provenance={"parent_event_keys": ((3, 8),)},
        )
        synthetic = NCloseCandidate(
            "synthetic-owner",
            (30, 31),
            origin="synthetic_test",
            synthetic=True,
            provenance={"source_query": "unitig-7"},
        )

        candidates = [primary, secondary]
        candidates.extend([combined, synthetic])
        without_secondary, first_rejections = apply_nclose_filter(
            candidates,
            "source_filter",
            lambda candidate: (
                "secondary"
                if candidate.provenance.get("source") == "secondary"
                else None
            ),
        )
        final, second_rejections = apply_nclose_filter(
            without_secondary,
            "synthetic_filter",
            lambda candidate: (
                "combined"
                if candidate.origin == "synthetic_combined_test"
                else None
            ),
        )

        self.assertEqual(
            [candidate.identity for candidate in candidates],
            [
                ("shared-owner", (8, 3)),
                ("shared-owner", (8, 3)),
                ("combined-owner", (20, 21)),
                ("synthetic-owner", (30, 31)),
            ],
        )
        self.assertIs(first_rejections[0].candidate, secondary)
        self.assertIs(second_rejections[0].candidate, combined)
        self.assertEqual(final, [primary, synthetic])
        self.assertIs(final[0], primary)
        self.assertIs(final[1], synthetic)
        self.assertEqual(final[1].origin, "synthetic_test")
        self.assertTrue(final[1].synthetic)
        self.assertEqual(final[1].provenance, {"source_query": "unitig-7"})

    def test_event_key_filter_removes_every_duplicate_and_reverse_occurrence(self):
        first = NCloseCandidate("owner-a", (8, 3), origin="paf")
        duplicate = NCloseCandidate(
            "owner-a",
            (8, 3),
            origin="paf",
            provenance={"occurrence": 2},
        )
        reverse = NCloseCandidate("owner-b", (3, 8), origin="paf")
        survivor = NCloseCandidate("owner-c", (9, 7), origin="paf")

        kept, rejected = apply_nclose_filter(
            [first, duplicate, reverse, survivor],
            "bam_event_key",
            lambda candidate: (
                "bam_rejected" if candidate.event_key == (3, 8) else None
            ),
        )

        self.assertEqual(kept, [survivor])
        self.assertIs(kept[0], survivor)
        self.assertEqual(
            [rejection.candidate for rejection in rejected],
            [first, duplicate, reverse],
        )
        self.assertEqual(
            [rejection.candidate.event_key for rejection in rejected],
            [(3, 8), (3, 8), (3, 8)],
        )

    def test_raw_virtual_inversion_filter_uses_event_key_without_deduping(self):
        subprocess = SimpleNamespace(
            run=Mock(return_value=SimpleNamespace(returncode=0))
        )
        collect_pairs = Mock(return_value={(3, 8)})
        apply_filter = load_stage01_function(
            "apply_raw_virtual_inversion_filter",
            pd=SimpleNamespace(read_csv=Mock(return_value=Mock())),
            TYPE2_NOISE_SIGMA_MULTIPLIER=1.5,
            estimate_type2_global_noise_sigma=Mock(return_value=1.0),
            import_censat_repeat_data=Mock(return_value={}),
            CENSAT_PATH="censat.bed",
            SKIP_BAM_ANAL=False,
            JULIA_BAM_THREAD_LIM=2,
            THREAD=4,
            args=SimpleNamespace(progress=False),
            subprocess=subprocess,
            os=os,
            PREFIX="/output",
            read_bam_loc="reads.bam",
            CHROMOSOME_INFO_FILE_PATH="reference.fai",
            main_stat_loc="depth.txt",
            collect_raw_virtual_inv_nclose_pairs=collect_pairs,
            __file__=str(BUILD_GRAPH_PATH),
        )
        first = NCloseCandidate("owner-a", (8, 3), origin="paf")
        duplicate = NCloseCandidate(
            "owner-a",
            (8, 3),
            origin="paf",
            provenance={"occurrence": 2},
        )
        reverse = NCloseCandidate("owner-b", (3, 8), origin="paf")
        survivor = NCloseCandidate(
            "owner-c",
            (9, 7),
            origin="synthetic_test",
            synthetic=True,
        )

        kept = apply_filter([first, duplicate, reverse, survivor])

        self.assertEqual(kept, [survivor])
        self.assertIs(kept[0], survivor)
        collect_pairs.assert_called_once_with("/output")

    def test_user_exclusion_removes_all_owner_occurrences_in_file_order(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            prefix = Path(tmp_dir)
            exclusion_path = prefix / "exclude.txt"
            exclusion_path.write_text(
                "owner-b\nowner-b\nINDEL_INDEX_7\nmissing\nowner-a\n",
                encoding="utf-8",
            )
            logger = Mock()
            apply_exclusions = load_stage01_function(
                "apply_user_nclose_exclusions",
                args=SimpleNamespace(
                    exclude_nclose_list_loc=str(exclusion_path)
                ),
                PREFIX=str(prefix),
                pkl=pickle,
                logging=logger,
            )
            owner_a_first = NCloseCandidate("owner-a", (1, 2))
            owner_b = NCloseCandidate("owner-b", (3, 4))
            owner_a_second = NCloseCandidate(
                "owner-a",
                (5, 6),
                provenance={"occurrence": 2},
            )
            survivor = NCloseCandidate(
                "owner-c",
                (7, 8),
                origin="synthetic_combined_test",
                synthetic=True,
            )

            kept, indel_indices = apply_exclusions(
                [owner_a_first, owner_b, owner_a_second, survivor]
            )

            self.assertEqual(kept, [survivor])
            self.assertIs(kept[0], survivor)
            self.assertEqual(indel_indices, {7})
            self.assertEqual(
                logger.warning.call_args_list[0].args[0],
                "Skipped contig : owner-b, owner-a",
            )
            with (prefix / "indel_exclude_idx_set.pkl").open("rb") as handle:
                self.assertEqual(pickle.load(handle), {7})

    def test_filter_order_and_first_rejection_are_stable(self):
        candidates = [
            NCloseCandidate("a", (0, 1)),
            NCloseCandidate("b", (2, 3)),
            NCloseCandidate("c", (4, 5)),
        ]

        kept, first_rejections = apply_nclose_filter(
            candidates,
            "censat",
            lambda candidate: (
                "FILTERED_02_CENSAT_MAPQ"
                if candidate.contig_name == "b"
                else None
            ),
        )
        seen_by_second_filter = []

        def second_reason(candidate):
            seen_by_second_filter.append(candidate.contig_name)
            if candidate.contig_name in {"a", "b"}:
                return "FILTERED_02_SUBTELO"
            return None

        final, second_rejections = apply_nclose_filter(
            kept,
            "subtelo",
            second_reason,
        )

        self.assertEqual([candidate.contig_name for candidate in kept], ["a", "c"])
        self.assertEqual(seen_by_second_filter, ["a", "c"])
        self.assertEqual([candidate.contig_name for candidate in final], ["c"])
        self.assertEqual(first_rejections[0].candidate.contig_name, "b")
        self.assertEqual(first_rejections[0].stage, "censat")
        self.assertEqual(
            first_rejections[0].reason,
            "FILTERED_02_CENSAT_MAPQ",
        )
        self.assertEqual(second_rejections[0].candidate.contig_name, "a")

    def test_nclose_calc_materializes_legacy_only_at_final_output_boundary(self):
        calls = ordered_direct_calls("nclose_calc")
        call_names = [name for name, _ in calls]

        self.assertEqual(
            call_names.count("candidates_to_legacy"),
            1,
            "Candidate-native processing must have one final legacy conversion",
        )
        self.assertNotIn(
            "refresh_nclose_candidates",
            call_names,
            "nclose_calc must not reconcile repeated Candidate/legacy round trips",
        )
        self.assertNotIn("iter_legacy_candidates", call_names)

        conversion_pos = call_names.index("candidates_to_legacy")
        final_output_pos = call_names.index("finalize_nclose_outputs")
        self.assertLess(conversion_pos, final_output_pos)
        self.assertLess(call_names.index("run_nclose_pipeline"), conversion_pos)


if __name__ == "__main__":
    unittest.main()
