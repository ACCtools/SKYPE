from __future__ import annotations

import ast
import logging
import unittest
from collections import Counter, defaultdict
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from nclose_candidate import NCloseCandidate, apply_nclose_filter


BUILD_GRAPH_PATH = (
    Path(__file__).resolve().parents[1] / "02_Build_Breakend_Graph_Limited.py"
)


CTG_NAM = 0
CTG_DIR = 4
CHR_NAM = 5
CHR_LEN = 6
CHR_STR = 7
CHR_END = 8
CTG_MAPQ = 9
CTG_CENSAT = 18
NCLOSE_COMPRESS_LIMIT = 100_000
OFFSET_DIR_GROUP_LIMIT = 100_000


def load_censat_helpers():
    """Load CenSat contracts without executing stage 02's CLI entry point."""

    definition_names = {
        "PregraphSourceMode",
        "NCloseSourceConfig",
        "CensatPairClass",
        "CensatNoncensatCandidate",
        "chr2int",
        "distance_checker",
        "get_corr_dir",
        "_flip_ctg_dir",
        "node_is_censat",
        "classify_censat_pair",
        "nclose_has_simple_alt",
        "build_censat_noncensat_candidate",
        "group_censat_noncensat_candidates",
        "iter_nclose_owner_pairs",
        "iter_censat_noncensat_candidates",
        "_cen_fragment_target_dir_from_meta",
        "collect_missing_cen_fragment_dir_censat_noncensat",
        "apply_offset_direction_mismatched_censat_noncensat_filter",
        "filter_offset_direction_mismatched_censat_noncensat",
        "_censat_at_chromosome_end",
        "_canonical_censat_pair_info",
        "_normalized_censat_endpoint_dirs",
        "apply_censat_censat_filter",
        "apply_censat_fragment_direction_filter",
        "should_add_censat_pair_rescue",
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
        "defaultdict": defaultdict,
        "dataclass": dataclass,
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
        "NCLOSE_COMPRESS_LIMIT": NCLOSE_COMPRESS_LIMIT,
        "OFFSET_DIR_GROUP_LIMIT": OFFSET_DIR_GROUP_LIMIT,
        "PIPELINE_MODE_KARYOTYPE": "karyotype",
    }
    module = ast.Module(body=definitions, type_ignores=[])
    exec(compile(module, str(BUILD_GRAPH_PATH), "exec"), namespace)
    return namespace


HELPERS = load_censat_helpers()
PregraphSourceMode = HELPERS["PregraphSourceMode"]
NCloseSourceConfig = HELPERS["NCloseSourceConfig"]
CensatPairClass = HELPERS["CensatPairClass"]
classify_censat_pair = HELPERS["classify_censat_pair"]
build_censat_noncensat_candidate = HELPERS[
    "build_censat_noncensat_candidate"
]
group_censat_noncensat_candidates = HELPERS[
    "group_censat_noncensat_candidates"
]
collect_missing_cen_fragment_dir_censat_noncensat = HELPERS[
    "collect_missing_cen_fragment_dir_censat_noncensat"
]
filter_offset_direction_mismatched_censat_noncensat = HELPERS[
    "filter_offset_direction_mismatched_censat_noncensat"
]
apply_censat_censat_filter = HELPERS["apply_censat_censat_filter"]
apply_censat_fragment_direction_filter = HELPERS[
    "apply_censat_fragment_direction_filter"
]
should_add_censat_pair_rescue = HELPERS["should_add_censat_pair_rescue"]


def ordered_call_names(function_name: str):
    tree = ast.parse(BUILD_GRAPH_PATH.read_text(encoding="utf-8"))
    function = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == function_name
    )
    calls = [node for node in ast.walk(function) if isinstance(node, ast.Call)]
    calls.sort(key=lambda node: (node.lineno, node.col_offset))
    return [
        call.func.id
        for call in calls
        if isinstance(call.func, ast.Name)
    ]


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


def append_one_censat_event(
    contig_data,
    name,
    noncensat_pos,
    *,
    censat_side=0,
    censat_dir="-",
    censat_chrom="chr1",
    noncensat_chrom="chr2",
    simple_alt=False,
):
    endpoint_name = f"simple_ctg_alt_{name}" if simple_alt else name
    censat = contig_row(
        endpoint_name,
        censat_chrom,
        2_000_000,
        2_000_100,
        direction=censat_dir,
        censat="rin",
    )
    noncensat = contig_row(
        endpoint_name,
        noncensat_chrom,
        noncensat_pos - 50,
        noncensat_pos + 50,
        direction="+",
        censat="0",
    )
    pair_rows = (censat, noncensat) if censat_side == 0 else (noncensat, censat)
    pair = (len(contig_data), len(contig_data) + 1)
    contig_data.extend(pair_rows)
    return pair


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

    def test_one_censat_view_preserves_path_order_and_normalizes_direction(self):
        contig_data = []
        first_pair = append_one_censat_event(
            contig_data,
            "first",
            100_000,
            censat_side=0,
            censat_dir="-",
        )
        second_pair = append_one_censat_event(
            contig_data,
            "second",
            200_001,
            censat_side=1,
            censat_dir="-",
            simple_alt=True,
        )

        first = build_censat_noncensat_candidate(
            contig_data, "first", first_pair
        )
        second = build_censat_noncensat_candidate(
            contig_data, "second", second_pair
        )

        self.assertEqual(first.pair, first_pair)
        self.assertEqual(first.censat_side, 0)
        self.assertEqual(first.censat_idx, first_pair[0])
        self.assertEqual(first.noncensat_idx, first_pair[1])
        self.assertEqual(first.censat_norm_dir, "-")
        self.assertEqual(first.noncensat_pos, 100_000)
        self.assertFalse(first.is_simple_alt)

        self.assertEqual(second.pair, second_pair)
        self.assertEqual(second.censat_side, 1)
        self.assertEqual(second.censat_idx, second_pair[1])
        self.assertEqual(second.noncensat_idx, second_pair[0])
        self.assertEqual(second.censat_norm_dir, "+")
        self.assertEqual(second.noncensat_pos, 200_001)
        self.assertTrue(second.is_simple_alt)

    def test_one_censat_builder_rejects_none_and_both(self):
        contig_data = [
            contig_row("none-a", "chr1", 10, 20, censat="0"),
            contig_row("none-b", "chr2", 10, 20, censat="0"),
            contig_row("both-a", "chr1", 10, 20, censat="r"),
            contig_row("both-b", "chr2", 10, 20, censat="rin"),
        ]

        self.assertIsNone(
            build_censat_noncensat_candidate(contig_data, "none", (0, 1))
        )
        self.assertIsNone(
            build_censat_noncensat_candidate(contig_data, "both", (2, 3))
        )


class CensatNoncensatGroupingTests(unittest.TestCase):
    def test_grouping_is_transitive_but_the_100kb_boundary_is_exclusive(self):
        contig_data = []
        pairs = {
            "a": append_one_censat_event(contig_data, "a", 100_000),
            "b": append_one_censat_event(contig_data, "b", 199_999),
            "c": append_one_censat_event(contig_data, "c", 299_998),
            "boundary": append_one_censat_event(
                contig_data, "boundary", 399_998
            ),
            "other-key": append_one_censat_event(
                contig_data,
                "other-key",
                150_000,
                noncensat_chrom="chr3",
            ),
        }
        candidates = [
            build_censat_noncensat_candidate(contig_data, name, pairs[name])
            for name in ("c", "boundary", "b", "other-key", "a")
        ]

        clusters = group_censat_noncensat_candidates(candidates)
        observed = sorted(
            tuple(candidate.contig_name for candidate in cluster)
            for cluster in clusters
        )

        self.assertEqual(
            observed,
            sorted([("a", "b", "c"), ("boundary",), ("other-key",)]),
        )

    def test_missing_direction_and_final_arbitration_share_group_contract(self):
        contig_data = []
        wrong_a = append_one_censat_event(
            contig_data, "wrong-a", 100_000, censat_dir="-"
        )
        wrong_b = append_one_censat_event(
            contig_data,
            "wrong-b",
            199_999,
            censat_side=1,
            censat_dir="+",  # reverse-path normalization makes this '-'
        )
        correct_at_boundary = append_one_censat_event(
            contig_data,
            "correct",
            299_999,
            censat_side=1,
            censat_dir="-",  # reverse-path normalization makes this '+'
        )
        simple_alt_wrong = append_one_censat_event(
            contig_data,
            "ignored-simple",
            120_000,
            censat_dir="-",
            simple_alt=True,
        )
        nclose_nodes = defaultdict(
            list,
            {
                "wrong-a": [wrong_a],
                "wrong-b": [wrong_b],
                "correct": [correct_at_boundary],
                "ordinary-owner": [simple_alt_wrong],
            },
        )
        meta = {"chr1": {"dir": False}}

        missing = collect_missing_cen_fragment_dir_censat_noncensat(
            contig_data, nclose_nodes, meta
        )

        self.assertEqual(len(missing), 1)
        self.assertEqual(
            [candidate.contig_name for candidate in missing[0]],
            ["wrong-a", "wrong-b"],
        )
        self.assertEqual(
            {candidate.censat_norm_dir for candidate in missing[0]},
            {"-"},
        )

        filtered, removed = filter_offset_direction_mismatched_censat_noncensat(
            contig_data, nclose_nodes, meta
        )
        self.assertEqual(removed, 0)
        self.assertEqual(dict(filtered), dict(nclose_nodes))

        # Move the correct candidate one base inside the neighboring group.
        correct_idx = correct_at_boundary[0]
        correct_row = list(contig_data[correct_idx])
        correct_row[CHR_STR] -= 1
        correct_row[CHR_END] -= 1
        contig_data[correct_idx] = tuple(correct_row)

        filtered, removed = filter_offset_direction_mismatched_censat_noncensat(
            contig_data, nclose_nodes, meta
        )
        self.assertEqual(removed, 3)
        self.assertEqual(filtered["wrong-a"], [])
        self.assertEqual(filtered["wrong-b"], [])
        self.assertEqual(filtered["correct"], [correct_at_boundary])
        self.assertEqual(filtered["ordinary-owner"], [])


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

    def test_both_filter_rejection_precedence_and_competitor_rule(self):
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

        self.assertEqual(kept, [])
        self.assertEqual(
            {item.candidate.contig_name: item.reason for item in rejected},
            {
                "opposite-a": "opposite_dir",
                "opposite-b": "opposite_dir",
                "same-chrom": "same_chrom_opposite_dir",
                "mapq": "mapq",
                "terminal": "terminal",
            },
        )

    def test_fragment_direction_uses_path_normalized_both_endpoints(self):
        contig_data = [
            contig_row("matched", "chr1", 1_000, 1_100, "+", "r"),
            contig_row("matched", "chr2", 2_000, 2_100, "-", "r"),
            contig_row("mismatch", "chr1", 3_000, 3_100, "+", "r"),
            contig_row("mismatch", "chr2", 4_000, 4_100, "+", "r"),
            contig_row("one", "chr1", 5_000, 5_100, "-", "r"),
            contig_row("one", "chr2", 6_000, 6_100, "+", "0"),
        ]
        candidates = [
            NCloseCandidate("matched", (0, 1)),
            NCloseCandidate("mismatch", (2, 3)),
            NCloseCandidate("one", (4, 5)),
        ]

        kept, rejected = apply_censat_fragment_direction_filter(
            candidates,
            contig_data,
            {"chr1": {"dir": False}, "chr2": {"dir": False}},
        )

        self.assertEqual(
            [candidate.contig_name for candidate in kept],
            ["matched", "one"],
        )
        self.assertEqual(
            [(item.candidate.contig_name, item.reason) for item in rejected],
            [("mismatch", "direction_mismatch")],
        )


class CensatPipelineOrderingTests(unittest.TestCase):
    def test_pair_rescue_source_gate_truth_table(self):
        eligible = NCloseSourceConfig(
            mode=PregraphSourceMode.CONFIGURED_PAF,
            paf_file_paths=("primary.paf", "unitig.paf"),
            original_paf_paths=("primary.raw.paf", "unitig.raw.paf"),
            is_unitig_reduced=False,
            secondary_candidate_paf="unitig.paf",
        )
        self.assertTrue(should_add_censat_pair_rescue("karyotype", eligible))

        ineligible_sources = {
            "primary retry": NCloseSourceConfig(
                mode=PregraphSourceMode.PRIMARY_ONLY_RETRY,
                paf_file_paths=eligible.paf_file_paths,
                original_paf_paths=eligible.original_paf_paths,
                is_unitig_reduced=False,
                secondary_candidate_paf="unitig.paf",
            ),
            "no secondary": NCloseSourceConfig(
                mode=PregraphSourceMode.CONFIGURED_PAF,
                paf_file_paths=eligible.paf_file_paths,
                original_paf_paths=eligible.original_paf_paths,
                is_unitig_reduced=False,
                secondary_candidate_paf=None,
            ),
            "unitig reduced": NCloseSourceConfig(
                mode=PregraphSourceMode.CONFIGURED_PAF,
                paf_file_paths=eligible.paf_file_paths,
                original_paf_paths=eligible.original_paf_paths,
                is_unitig_reduced=True,
                secondary_candidate_paf="unitig.paf",
            ),
            "single aligned PAF": NCloseSourceConfig(
                mode=PregraphSourceMode.CONFIGURED_PAF,
                paf_file_paths=("primary.paf",),
                original_paf_paths=eligible.original_paf_paths,
                is_unitig_reduced=False,
                secondary_candidate_paf="unitig.paf",
            ),
            "single original PAF": NCloseSourceConfig(
                mode=PregraphSourceMode.CONFIGURED_PAF,
                paf_file_paths=eligible.paf_file_paths,
                original_paf_paths=("primary.raw.paf",),
                is_unitig_reduced=False,
                secondary_candidate_paf="unitig.paf",
            ),
        }
        for case_name, source in ineligible_sources.items():
            with self.subTest(case_name=case_name):
                self.assertFalse(
                    should_add_censat_pair_rescue("karyotype", source)
                )

        self.assertFalse(should_add_censat_pair_rescue("other", eligible))

    def test_pair_rescue_stays_after_ordinary_censat_arbitration(self):
        calls = ordered_call_names("nclose_calc")
        ordered_contract = [
            "apply_censat_censat_filter",
            "apply_censat_fragment_direction_filter",
            "apply_simple_alt_preference_filter",
            "collect_missing_cen_fragment_dir_censat_noncensat",
            "add_nearest_combined_censat_noncensat_ncloses",
            "apply_offset_direction_mismatched_censat_noncensat_filter",
            "add_censat_pair_rescue_ncloses",
            "apply_nclose_count_vaf_filter",
        ]

        positions = [calls.index(call_name) for call_name in ordered_contract]
        self.assertEqual(
            positions,
            sorted(positions),
            "CenSat repair/filter/rescue lifecycle changed",
        )


if __name__ == "__main__":
    unittest.main()
