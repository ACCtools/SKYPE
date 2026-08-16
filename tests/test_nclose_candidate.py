from __future__ import annotations

import unittest
from collections import defaultdict

from nclose_candidate import (
    NCloseCandidate,
    apply_nclose_filter,
    candidates_to_legacy,
    endpoint_metadata,
    iter_legacy_candidates,
)


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
            origin="combined_censat_noncensat",
            synthetic=True,
            provenance={"parent_event_keys": ((3, 8),), "type2_query": "utg7"},
        )

        self.assertTrue(candidate.synthetic)
        self.assertEqual(candidate.origin, "combined_censat_noncensat")
        self.assertEqual(candidate.provenance["parent_event_keys"], ((3, 8),))

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


if __name__ == "__main__":
    unittest.main()
