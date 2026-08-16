from __future__ import annotations

import ast
import os
import unittest
from dataclasses import dataclass, fields
from enum import Enum
from pathlib import Path


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
