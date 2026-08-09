from __future__ import annotations

import ast
import gzip
import json
import pickle
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import networkx as nx


SKYPE_ROOT = Path(__file__).resolve().parents[1]
BUILD_GRAPH_PATH = SKYPE_ROOT / "02_Build_Breakend_Graph_Limited.py"


def load_indel_graph_helpers():
    tree = ast.parse(BUILD_GRAPH_PATH.read_text(encoding="utf-8"))
    function_names = {
        "should_add_indel_graph",
        "get_graph_dimension_increment",
        "point_overlaps_censat",
        "type4_indel_expected_high_sides",
        "type4_indel_depth_supported",
        "make_graph",
    }
    functions = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name in function_names
    ]
    namespace = {
        "PIPELINE_MODE_KARYOTYPE": "karyotype",
        "TYPE4_INDEL_GRAPH_DIMENSION_INCREMENT": 1,
        "TYPE4_INDEL_GRAPH_DEPTH_DIFF_RATIO": 0.2,
        "CTG_NAM": 0,
        "CTG_DIR": 4,
        "CHR_NAM": 5,
        "CHR_STR": 7,
        "CHR_END": 8,
        "nx": nx,
    }
    module = ast.Module(body=functions, type_ignores=[])
    exec(compile(module, str(BUILD_GRAPH_PATH), "exec"), namespace)
    return namespace


INDEL_GRAPH_HELPERS = load_indel_graph_helpers()
should_add_indel_graph = INDEL_GRAPH_HELPERS["should_add_indel_graph"]
get_graph_dimension_increment = INDEL_GRAPH_HELPERS[
    "get_graph_dimension_increment"
]
type4_indel_depth_supported = INDEL_GRAPH_HELPERS[
    "type4_indel_depth_supported"
]


def graph_node(name: str, direction: str, chromosome: str) -> list:
    node = [None] * 9
    node[0] = name
    node[4] = direction
    node[5] = chromosome
    node[7] = 0
    node[8] = 100
    return node


def build_indel_dimension_test_graph(graph_depth: int):
    namespace = INDEL_GRAPH_HELPERS
    namespace["contig_data"] = [
        graph_node("ctg", "+", "chr1"),
        graph_node("ctg", "+", "chr1"),
    ]
    namespace["bnd_graph_adjacency"] = {
        (1, 0): [[1, 1]],
        (1, 1): [],
    }
    namespace["type4_indel_dimension_edge_set"] = {
        ((1, 0), (1, 1)),
    }
    namespace["distance_checker"] = lambda _left, _right: 0
    namespace["overlap_calculator"] = lambda _left, _right: 0
    return namespace["make_graph"](graph_depth, 0)


class IndelGraphModeTests(unittest.TestCase):
    def test_karyotype_mode_enables_indel_graph_automatically(self) -> None:
        self.assertTrue(should_add_indel_graph("karyotype"))

    def test_variant_mode_keeps_indel_graph_opt_in(self) -> None:
        self.assertFalse(should_add_indel_graph("variant"))
        self.assertTrue(
            should_add_indel_graph("variant", explicitly_enabled=True)
        )

    def test_indel_edge_consumes_one_graph_depth_dimension(self) -> None:
        left = graph_node("ctg", "+", "chr1")
        right = graph_node("ctg", "+", "chr1")

        self.assertEqual(
            get_graph_dimension_increment(
                left,
                right,
                is_type4_indel_edge=True,
            ),
            (1, 0),
        )
        self.assertEqual(get_graph_dimension_increment(left, right), (0, 0))

    def test_existing_graph_dimension_rules_are_preserved(self) -> None:
        left = graph_node("ctg", "+", "chr1")
        translocation = graph_node("other", "+", "chr2")
        inversion = graph_node("ctg", "-", "chr1")

        self.assertEqual(
            get_graph_dimension_increment(left, translocation),
            (1, 0),
        )
        self.assertEqual(
            get_graph_dimension_increment(left, inversion),
            (0, 1),
        )

    def test_make_graph_increments_indel_state_and_honors_depth_limit(self) -> None:
        graph = build_indel_dimension_test_graph(graph_depth=1)

        self.assertTrue(graph.has_edge((1, 0, 0, 0), (1, 1, 1, 0)))
        self.assertFalse(graph.has_edge((1, 0, 0, 0), (1, 1, 0, 0)))
        self.assertFalse(graph.has_edge((1, 0, 1, 0), (1, 1, 2, 0)))

        zero_depth_graph = build_indel_dimension_test_graph(graph_depth=0)
        self.assertEqual(zero_depth_graph.number_of_edges(), 0)

    def test_single_censat_endpoint_remains_eligible(self) -> None:
        INDEL_GRAPH_HELPERS["expected_depth_step"] = (
            lambda *_args, **_kwargs: 10.0
        )
        for censat_intervals in (
            [(90, 110)],
            [(190, 210)],
        ):
            with self.subTest(censat_intervals=censat_intervals):
                supported, steps, reason = type4_indel_depth_supported(
                    {},
                    {"chr3": censat_intervals},
                    "chr3",
                    100,
                    200,
                    "insertion",
                    10.0,
                )
                self.assertTrue(supported)
                self.assertEqual(steps, [10.0, 10.0])
                self.assertEqual(reason, "pass")

    def test_both_censat_endpoints_are_rejected(self) -> None:
        supported, steps, reason = type4_indel_depth_supported(
            {},
            {"chr3": [(90, 110), (190, 210)]},
            "chr3",
            100,
            200,
            "insertion",
            10.0,
        )

        self.assertFalse(supported)
        self.assertIsNone(steps)
        self.assertEqual(reason, "both_endpoints_censat")


class NativeType4GraphStage11Tests(unittest.TestCase):
    def test_graph_selected_native_type4_is_not_a_step11_column(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "output"
            output.mkdir()
            fai = root / "reference.fai"
            original_paf = root / "original.paf"
            ppc = root / "assembly.ppc.paf"
            fai.write_text("chr1\t5000000\n", encoding="utf-8")
            original_paf.write_text("", encoding="utf-8")

            def ppc_row(query_start: int, query_end: int,
                        ref_start: int, ref_end: int) -> str:
                return "\t".join(map(str, [
                    "native_type4", 200000, query_start, query_end, "+",
                    "chr1", 5000000, ref_start, ref_end, 60, 4, 0, 1,
                    "0", "0", "0", "0", "0", "0", "+", "chr1",
                    "99.0",
                ])) + "\n"

            ppc.write_text(
                ppc_row(0, 100000, 0, 100000)
                + ppc_row(100000, 200000, 1100000, 1200000),
                encoding="utf-8",
            )
            with (output / "paf_file_path.pkl").open("wb") as handle:
                pickle.dump([str(original_paf)], handle)
            with (output / "conjoined_type4_ins_del.pkl").open("wb") as handle:
                pickle.dump(([], []), handle)
            with (output / "nclose_event_catalog.pkl").open("wb") as handle:
                pickle.dump([], handle)
            graph_edges = [
                {
                    "src": (1, 0), "dst": (1, 1),
                    "type4_tuple": (0, 1),
                    "indel_kind": "deletion",
                    "base_chrom": "chr1",
                    "base_st": 100000,
                    "base_nd": 1100000,
                    "contig_name": "native_type4",
                    "graph_dimension_increment": 1,
                },
                {
                    "src": (0, 1), "dst": (0, 0),
                    "type4_tuple": (0, 1),
                    "indel_kind": "deletion",
                    "base_chrom": "chr1",
                    "base_st": 100000,
                    "base_nd": 1100000,
                    "contig_name": "native_type4",
                    "graph_dimension_increment": 1,
                },
            ]
            with (output / "type4_indel_graph_edges.pkl").open("wb") as handle:
                pickle.dump(graph_edges, handle)

            result = subprocess.run(
                [
                    sys.executable,
                    str(SKYPE_ROOT / "11_Ref_Outlier_Contig_Modify.py"),
                    str(fai),
                    str(ppc),
                    str(output),
                ],
                cwd=SKYPE_ROOT,
                capture_output=True,
                text=True,
            )

            self.assertEqual(
                result.returncode,
                0,
                msg=result.stdout + "\n" + result.stderr,
            )
            outlier_root = output / "11_ref_ratio_outliers"
            self.assertFalse(any(outlier_root.rglob("*_base.paf")))
            with (output / "nclose_event_catalog.pkl").open("rb") as handle:
                catalog = pickle.load(handle)
            self.assertEqual(len(catalog), 1)
            self.assertTrue(catalog[0]["graph_only"])
            self.assertEqual(
                catalog[0]["source"],
                "type4_indel_graph_edges.pkl",
            )


class VcfAddIndelGraphIntegrationTests(unittest.TestCase):
    def test_vcf_mode_adds_del_graph_edges_but_keeps_ins_step11_only(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "output"
            fai = root / "reference.fai"
            paf = root / "assembly.paf"
            ins_paf = root / "ins.paf"
            telomere_bed = root / "telomere.bed"
            repeat_bed = root / "repeat.bed"
            censat_bed = root / "censat.bed"
            depth_path = root / "depth.win.stat.gz"
            vcf = root / "input.vcf"
            native_limits = root / "limit_combinations.json"

            chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
            fai.write_text(
                "".join(f"{chrom}\t5000000\n" for chrom in chromosomes),
                encoding="utf-8",
            )
            paf.write_text("", encoding="utf-8")
            telomere_bed.write_text("", encoding="utf-8")
            repeat_bed.write_text("", encoding="utf-8")
            censat_bed.write_text(
                "chr1\t1100000\t2900000\n", encoding="utf-8"
            )
            vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "##source=Sniffles2_2.2\n"
                "##contig=<ID=chr1,length=5000000>\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "chr1\t200000\tinv\tN\t<INV>\t60\tPASS\t"
                "SVTYPE=INV;END=800000\n"
                "chr1\t1000000\tdel\tN\t<DEL>\t60\tPASS\t"
                "SVTYPE=DEL;END=3000000;SVLEN=-2000000\n"
                "chr1\t500000\tins\tN\t<INS>\t60\tPASS\t"
                "SVTYPE=INS;SVLEN=1200000\n",
                encoding="utf-8",
            )
            native_limits.write_text(
                json.dumps({"limit_combinations": [1, 0]}),
                encoding="utf-8",
            )
            ins_paf.write_text(
                "vcf_ins_7_ins\t1200000\t0\t1200000\t+\tchr1\t5000000\t"
                "500000\t501000\t1000\t1200000\t60\ttp:A:P\tcs:Z::1000\n",
                encoding="utf-8",
            )
            with gzip.open(depth_path, "wt", encoding="utf-8") as handle:
                for st in range(0, 5_000_000, 100_000):
                    nd = st + 100_000
                    mean_depth = 10 if 1_000_000 <= st < 3_000_000 else 30
                    handle.write(
                        f"chr1\t{st}\t{nd}\t100000\t100000\t"
                        f"{mean_depth * 100000}\t100\t{mean_depth}\n"
                    )

            build_command = [
                sys.executable,
                str(SKYPE_ROOT / "02_Build_Breakend_Graph_Limited.py"),
                str(paf),
                str(fai),
                str(telomere_bed),
                str(repeat_bed),
                str(censat_bed),
                str(depth_path),
                str(output),
                str(root / "unused.bam"),
                "--alt",
                str(ins_paf),
                "--vcf_input",
                str(vcf),
                "--limit_combinations",
                str(native_limits),
                "--add_indel_graph",
                "--skip_bam_analysis",
                "-t",
                "1",
                "-d",
                "0",
            ]
            build_result = subprocess.run(
                build_command,
                cwd=SKYPE_ROOT,
                capture_output=True,
                text=True,
            )
            self.assertEqual(
                build_result.returncode,
                0,
                msg=build_result.stdout + "\n" + build_result.stderr,
            )
            self.assertIn(
                "Using fixed limit_combinations",
                build_result.stdout + build_result.stderr,
            )
            self.assertEqual(
                json.loads(
                    (output / "limit_combinations.json").read_text(
                        encoding="utf-8"
                    )
                )["limit_combinations"],
                [1, 0],
            )

            summary = json.loads(
                (output / "vcf_mode_summary.json").read_text(encoding="utf-8")
            )
            self.assertEqual(summary["vcf_type4_graph_events"], 1)
            self.assertEqual(summary["vcf_type4_graph_nodes"], 2)
            self.assertEqual(summary["vcf_type4_graph_ins_excluded"], 1)
            self.assertEqual(summary["nclose_pairs"], 2)

            with (output / "vcf_type4_events.pkl").open("rb") as handle:
                step11_events = pickle.load(handle)
            self.assertEqual(
                [event["svtype"] for event in step11_events], ["DEL", "INS"]
            )

            with (output / "type4_indel_graph_edges.pkl").open("rb") as handle:
                graph_edges = pickle.load(handle)
            self.assertEqual(len(graph_edges), 2)
            self.assertEqual(
                {edge["graph_dimension_increment"] for edge in graph_edges},
                {1},
            )
            self.assertEqual(
                {edge["indel_kind"] for edge in graph_edges}, {"deletion"}
            )
            self.assertTrue(
                all(
                    str(edge["contig_name"]).startswith("vcf_type4_graph_")
                    for edge in graph_edges
                )
            )

            with (output / "nclose_chunk_data.pkl").open("rb") as handle:
                canonical_nclose, _, _ = pickle.load(handle)
            ppc_rows = [
                line.rstrip("\n").split("\t")
                for line in (output / "assembly.paf.ppc.paf").read_text(
                    encoding="utf-8"
                ).splitlines()
            ]
            self.assertEqual(sum(map(len, canonical_nclose.values())), 2)
            for pair_list in canonical_nclose.values():
                for node_a_idx, node_b_idx in pair_list:
                    node_a = ppc_rows[node_a_idx]
                    node_b = ppc_rows[node_b_idx]
                    self.assertFalse(
                        node_a[5] == node_b[5] and node_a[4] == node_b[4],
                        msg=f"Indel-like pair leaked into canonical nclose: "
                        f"{node_a_idx}, {node_b_idx}",
                    )

            ppc_path = output / "assembly.paf.ppc.paf"
            stale_outlier_root = output / "11_ref_ratio_outliers"
            for folder_name in (
                "front_jump",
                "back_jump",
                "ecdna",
                "type2_ins",
            ):
                stale_folder = stale_outlier_root / folder_name
                stale_folder.mkdir(parents=True, exist_ok=True)
                (stale_folder / "999.paf").write_text(
                    "stale\n", encoding="utf-8"
                )
                (stale_folder / "999.win.stat.gz").write_bytes(b"stale")

            step11_result = subprocess.run(
                [
                    sys.executable,
                    str(SKYPE_ROOT / "11_Ref_Outlier_Contig_Modify.py"),
                    str(fai),
                    str(ppc_path),
                    str(output),
                ],
                cwd=SKYPE_ROOT,
                capture_output=True,
                text=True,
            )
            self.assertEqual(
                step11_result.returncode,
                0,
                msg=step11_result.stdout + "\n" + step11_result.stderr,
            )
            with (output / "vcf_type4_outlier_index.pkl").open("rb") as handle:
                outlier_index = pickle.load(handle)
            self.assertEqual(
                {event["svtype"] for event in outlier_index.values()},
                {"INS"},
            )
            self.assertIn(
                "Graph-selected type4 Indels excluded from step 11",
                step11_result.stdout + step11_result.stderr,
            )
            self.assertFalse(
                any((stale_outlier_root / "front_jump").glob("*_base.paf")),
                msg="graph-selected DEL leaked into a standalone step-11 column",
            )
            with (output / "nclose_event_catalog.pkl").open("rb") as handle:
                event_catalog = pickle.load(handle)
            graph_catalog_events = [
                event for event in event_catalog
                if event.get("source") == "type4_indel_graph_edges.pkl"
            ]
            step11_catalog_events = [
                event for event in event_catalog
                if event.get("kind") == "indel" and not event.get("graph_only")
            ]
            self.assertEqual(
                [event.get("svtype") for event in graph_catalog_events],
                ["DEL"],
            )
            self.assertEqual(
                [event.get("svtype") for event in step11_catalog_events],
                ["INS"],
            )
            self.assertFalse(
                any(stale_outlier_root.rglob("999.paf")),
                msg="stage 11 retained an outlier PAF from a previous run",
            )
            self.assertFalse(
                any(stale_outlier_root.rglob("999.win.stat.gz")),
                msg="stage 11 retained a depth file from a previous run",
            )
            self.assertFalse((stale_outlier_root / "ecdna").exists())
            self.assertFalse((stale_outlier_root / "type2_ins").exists())


if __name__ == "__main__":
    unittest.main()
