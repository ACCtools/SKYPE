from __future__ import annotations

import ast
import collections
import csv
import os
from pathlib import Path
import pickle
import re
import tempfile
import unittest

import vcfpy

from bp_step_depth_ratio import (
    BP_STEP_DEPTH_RATIO_B, BP_STEP_DEPTH_RATIO_PREDICT_B, add_ratio_info,
)
from native_bnd import (
    constituent_pairs, merged_type4_pairs, layout_endpoints,
)


def node(name, qs, qe, strand, chrom, rs, re):
    return [name, 1000, qs, qe, strand, chrom, 1_000_000, rs, re]


def export_namespace():
    # Stage 31 is an executable script. Load only the native export functions,
    # exercising the real vcfpy serialization without running plots or NNLS.
    source = Path(__file__).resolve().parents[1] / "31_depth_analysis.py"
    names = {
        "bnd_alt", "choose_alt_forms", "make_strands", "invert_strand",
        "build_vcf_header", "vcf_filter_values", "write_bnd_vcf_pair",
        "write_symbolic_vcf_record",
    }
    tree = ast.parse(source.read_text())
    functions = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    import logging
    namespace = dict(
        collections=collections, vcfpy=vcfpy, re=re, pkl=pickle, logging=logging,
        os=os, csv=csv,
        add_ratio_info=add_ratio_info,
        BP_STEP_DEPTH_RATIO_B=BP_STEP_DEPTH_RATIO_B,
        BP_STEP_DEPTH_RATIO_PREDICT_B=BP_STEP_DEPTH_RATIO_PREDICT_B,
        SKYPE_VCF_SOURCE="SKYPE", N=10.0,
        RAW_TRANSLOCATION_RESULT_PKL="raw_translocation_result.pkl",
    )
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(source), "exec"), namespace)
    return namespace


class NativeBndTests(unittest.TestCase):
    def setUp(self):
        self.nodes = [
            node("u1", 0, 100, "+", "chr1", 100, 200),
            node("u1", 100, 200, "-", "chr1", 400, 500),
            node("u2", 0, 100, "-", "chr1", 210, 310),
            node("u2", 100, 200, "+", "chr1", 410, 510),
        ]

    def test_compound_pairs_remain_original_ncloses(self):
        pairs = merged_type4_pairs({"type2_merge_idx": 2}, ([(2, 3, 0, 1)], [(0, 1, 2, 3)]))
        self.assertEqual(pairs, ((0, 1), (2, 3)))
        self.assertEqual(constituent_pairs((0, 1, 0, 1)), ((0, 1), (0, 1)))

    def test_virtual_uses_saved_junction_facing_coordinates(self):
        record = virtual_record()
        endpoints, _ = layout_endpoints(record["layout_a"])
        self.assertEqual(endpoints, (("chr2", 105573362, "L"), ("chr6", 35361412, "L")))
        endpoints, _ = layout_endpoints(record["layout_b"])
        self.assertEqual(endpoints, (("chr2", 105573372, "R"), ("chr6", 35361267, "R")))

    def test_all_four_strand_combinations_use_retained_anchor_bases(self):
        ns = export_namespace()
        for da in "+-":
            for db in "+-":
                with self.subTest(da=da, db=db), tempfile.TemporaryDirectory() as tmp:
                    path = Path(tmp) / "calls.vcf"
                    header = ns["build_vcf_header"]({"chr1": 1000, "chr2": 1000})
                    with vcfpy.Writer.from_path(path, header) as writer:
                        ns["write_bnd_vcf_pair"](
                            writer, "check", "chr1", 200, da, "chr2", 300, db, 1.0, "u")
                    with vcfpy.Reader.from_path(path) as reader:
                        a, b = list(reader)
                    self.assertEqual(a.POS, 200 + (da == "-"))
                    self.assertEqual(b.POS, 300 + (db == "+"))
                    self.assertEqual(a.ALT[0].mate_pos, b.POS)
                    self.assertEqual(b.ALT[0].mate_pos, a.POS)
                    self.assertEqual(a.INFO["STRANDS"], b.INFO["STRANDS"][::-1])


def virtual_record():
    def endpoint(chrom, coord, direction, name):
        return dict(chrom=chrom, coord=coord, dir=direction, ctg_name=name,
                    ref_st=coord-1000, ref_nd=coord+1000)
    return dict(
        pair_id=7, nclose_key_a=(1944, 1945), nclose_key_b=(976, 977),
        layout_a=dict(ordered_endpoints=(
            endpoint("chr2", 105573362, "+", "utg014247l"),
            endpoint("chr6", 35361412, "-", "utg014247l"))),
        layout_b=dict(ordered_endpoints=(
            endpoint("chr2", 105573372, "-", "utg006773l"),
            endpoint("chr6", 35361267, "+", "utg006773l"))),
    )


if __name__ == "__main__":
    unittest.main()
