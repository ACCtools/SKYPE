from __future__ import annotations

import ast
import os
import unittest
from dataclasses import dataclass, fields
from enum import Enum
from pathlib import Path
from unittest.mock import Mock, patch


BUILD_GRAPH_PATH = (
    Path(__file__).resolve().parents[1] / "02_Build_Breakend_Graph_Limited.py"
)
VCF_SYNTHETIC_PAF_NAME = "vcf_synthetic.paf"


def load_pregraph_helpers():
    """Load the pregraph data contracts without running stage 02."""

    definition_names = {
        "PregraphSourceMode",
        "NCloseSourceConfig",
        "NCloseBuildResult",
        "extract_source_censat_type2_candidates",
        "resolve_pregraph_source",
    }
    tree = ast.parse(BUILD_GRAPH_PATH.read_text(encoding="utf-8"))
    definitions = [
        node
        for node in tree.body
        if isinstance(node, (ast.ClassDef, ast.FunctionDef))
        and node.name in definition_names
    ]
    module = ast.Module(body=definitions, type_ignores=[])
    namespace = {
        "dataclass": dataclass,
        "Enum": Enum,
        "os": os,
        "VCF_SYNTHETIC_PAF_NAME": VCF_SYNTHETIC_PAF_NAME,
    }
    exec(compile(module, str(BUILD_GRAPH_PATH), "exec"), namespace)
    return namespace


HELPERS = load_pregraph_helpers()
PregraphSourceMode = HELPERS["PregraphSourceMode"]
NCloseSourceConfig = HELPERS["NCloseSourceConfig"]
NCloseBuildResult = HELPERS["NCloseBuildResult"]
resolve_pregraph_source = HELPERS["resolve_pregraph_source"]
extract_source_censat_type2_candidates = HELPERS[
    "extract_source_censat_type2_candidates"
]


def args_alt_accesses(function_name: str) -> list[ast.Attribute]:
    """Return direct ``args.alt`` accesses in one top-level function."""

    tree = ast.parse(BUILD_GRAPH_PATH.read_text(encoding="utf-8"))
    function = next(
        node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == function_name
    )
    return [
        node
        for node in ast.walk(function)
        if isinstance(node, ast.Attribute)
        and node.attr == "alt"
        and isinstance(node.value, ast.Name)
        and node.value.id == "args"
    ]


class PregraphSourceTests(unittest.TestCase):
    def test_configured_singleton_uses_only_primary_paf(self) -> None:
        source = resolve_pregraph_source(
            PregraphSourceMode.CONFIGURED_PAF,
            "/inputs/primary.paf",
            None,
            ["/originals/primary.paf"],
            "/output",
        )

        self.assertIsInstance(source, NCloseSourceConfig)
        self.assertIs(source.mode, PregraphSourceMode.CONFIGURED_PAF)
        self.assertEqual(source.paf_file_paths, ("/inputs/primary.paf",))
        self.assertEqual(
            source.original_paf_paths, ("/originals/primary.paf",)
        )
        self.assertFalse(source.is_unitig_reduced)
        self.assertIsNone(source.secondary_candidate_paf)

    def test_configured_dual_preserves_primary_and_alternative_order(self) -> None:
        source = resolve_pregraph_source(
            PregraphSourceMode.CONFIGURED_PAF,
            "/inputs/primary.paf",
            "/inputs/alternative.paf",
            ["/originals/primary.paf", "/originals/alternative.paf"],
            "/output",
        )

        self.assertEqual(
            source.paf_file_paths,
            ("/inputs/primary.paf", "/inputs/alternative.paf"),
        )
        self.assertEqual(
            source.original_paf_paths,
            ("/originals/primary.paf", "/originals/alternative.paf"),
        )
        self.assertFalse(source.is_unitig_reduced)
        self.assertEqual(
            source.secondary_candidate_paf, "/inputs/alternative.paf"
        )

    def test_vcf_uses_synthetic_paf_then_assembly_anchor(self) -> None:
        source = resolve_pregraph_source(
            PregraphSourceMode.VCF,
            "/inputs/assembly.paf",
            "/inputs/ins-alt.paf",
            ["/originals/ignored.paf"],
            "/output",
        )

        self.assertIs(source.mode, PregraphSourceMode.VCF)
        self.assertEqual(
            source.paf_file_paths,
            ("/output/vcf_synthetic.paf", "/inputs/assembly.paf"),
        )
        self.assertEqual(source.original_paf_paths, ())
        self.assertFalse(source.is_unitig_reduced)
        self.assertIsNone(source.secondary_candidate_paf)

    def test_primary_only_retry_exactly_duplicates_primary_sources(self) -> None:
        source = resolve_pregraph_source(
            PregraphSourceMode.PRIMARY_ONLY_RETRY,
            "/inputs/primary.paf",
            "/inputs/alternative.paf",
            ["/originals/primary.paf", "/originals/alternative.paf"],
            "/output",
        )

        self.assertIs(source.mode, PregraphSourceMode.PRIMARY_ONLY_RETRY)
        self.assertEqual(
            source.paf_file_paths,
            ("/inputs/primary.paf", "/inputs/primary.paf"),
        )
        self.assertEqual(
            source.original_paf_paths,
            ("/originals/primary.paf", "/originals/primary.paf"),
        )
        self.assertTrue(source.is_unitig_reduced)
        self.assertIsNone(source.secondary_candidate_paf)

    def test_configured_source_delegates_censat_type2_to_explicit_alt(self) -> None:
        source = resolve_pregraph_source(
            PregraphSourceMode.CONFIGURED_PAF,
            "/inputs/primary.paf",
            "/inputs/alternative.paf",
            ["/originals/primary.paf", "/originals/alternative.paf"],
            "/output",
        )
        repeat_censat_data = {"chr1": [(100, 200)]}
        expected = [{"ctg_name": "unitig-candidate"}]
        delegate = Mock(return_value=expected)

        with patch.dict(
            extract_source_censat_type2_candidates.__globals__,
            {"extract_raw_censat_type2_candidates": delegate},
        ):
            observed = extract_source_censat_type2_candidates(
                source, repeat_censat_data
            )

        self.assertIs(observed, expected)
        delegate.assert_called_once_with(
            "/inputs/alternative.paf", repeat_censat_data
        )

    def test_primary_only_retry_skips_censat_type2_delegate(self) -> None:
        source = resolve_pregraph_source(
            PregraphSourceMode.PRIMARY_ONLY_RETRY,
            "/inputs/primary.paf",
            "/inputs/alternative.paf",
            ["/originals/primary.paf", "/originals/alternative.paf"],
            "/output",
        )
        delegate = Mock(side_effect=AssertionError("unitig evidence was read"))

        with patch.dict(
            extract_source_censat_type2_candidates.__globals__,
            {"extract_raw_censat_type2_candidates": delegate},
        ):
            observed = extract_source_censat_type2_candidates(source, {})

        self.assertEqual(observed, [])
        delegate.assert_not_called()

    def test_paf_pregraph_builders_do_not_access_raw_cli_alt(self) -> None:
        for function_name in ("contig_preprocessing_00", "nclose_calc"):
            with self.subTest(function_name=function_name):
                accesses = args_alt_accesses(function_name)
                self.assertEqual(
                    accesses,
                    [],
                    f"{function_name} must use NCloseSourceConfig, not args.alt",
                )

    def test_nclose_build_result_has_explicit_stable_contract(self) -> None:
        self.assertEqual(
            [field.name for field in fields(NCloseBuildResult)],
            [
                "df",
                "no_chrY",
                "repeat_censat_data",
                "chr_len",
                "contig_data",
                "contig_data_size",
                "chr_corr",
                "chr_rev_corr",
                "telo_contig",
                "telo_node_count",
                "telo_set",
                "rpt_con",
                "rpt_censat_con",
                "bnd_contig",
                "raw_nclose_nodes",
                "nclose_nodes",
                "nclose_start_compress",
                "nclose_end_compress",
                "vctg_dict",
                "all_nclose_comp",
                "nclose_coverage",
                "nclose_compress_track",
                "st_compress",
                "ed_compress",
                "uncomp_node_count",
                "nclose_node_count",
                "transloc_nclose_pair_count",
                "indel_exclude_idx_set",
                "telo_coverage",
            ],
        )


if __name__ == "__main__":
    unittest.main()
