from __future__ import annotations

import ast
import unittest
from pathlib import Path


RUN_DEPTH_PATH = Path(__file__).resolve().parents[1] / "21_run_depth.py"


def load_depth_paf_formatter():
    """Load the pure helper without executing the stage-21 pipeline script."""
    tree = ast.parse(RUN_DEPTH_PATH.read_text(encoding="utf-8"))
    function = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name == "format_nonzero_depth_paf_row"
    )
    module = ast.Module(body=[function], type_ignores=[])
    namespace = {}
    exec(compile(module, str(RUN_DEPTH_PATH), "exec"), namespace)
    return namespace["format_nonzero_depth_paf_row"]


format_nonzero_depth_paf_row = load_depth_paf_formatter()


def load_virtual_bridge_check():
    """Load the stage-21 virtual bridge strand check with its constants."""
    tree = ast.parse(RUN_DEPTH_PATH.read_text(encoding="utf-8"))
    function = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name == "virtual_bridge_follows_strand"
    )
    module = ast.Module(body=[function], type_ignores=[])
    namespace = {"DIR_FOR": 1, "DIR_BAK": 0, "CTG_DIR": 4, "CHR_STR": 7, "CHR_END": 8}
    exec(compile(module, str(RUN_DEPTH_PATH), "exec"), namespace)
    return namespace["virtual_bridge_follows_strand"]


virtual_bridge_follows_strand = load_virtual_bridge_check()


def load_ecdna_circuit_segments():
    """Load the stage-21 ecDNA circuit layout with its constants."""
    tree = ast.parse(RUN_DEPTH_PATH.read_text(encoding="utf-8"))
    function = next(
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name == "ecdna_circuit_segments"
    )
    module = ast.Module(body=[function], type_ignores=[])
    namespace = {"DIR_FOR": 1, "DIR_BAK": 0, "BND_TYPE": 0, "CTG_IN_TYPE": 1}
    exec(compile(module, str(RUN_DEPTH_PATH), "exec"), namespace)
    return namespace["ecdna_circuit_segments"]


ecdna_circuit_segments = load_ecdna_circuit_segments()


def load_bnd_raw_contig_list(contig_data, graph):
    """Load the stage-21 BND path builder bound to the given rows and graph."""
    import networkx as nx

    tree = ast.parse(RUN_DEPTH_PATH.read_text(encoding="utf-8"))
    functions = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef)
        and node.name in ("bnd_raw_contig_list", "distance_checker")
    ]
    module = ast.Module(body=functions, type_ignores=[])
    namespace = {
        "nx": nx, "G": graph, "contig_data": contig_data,
        "CTG_NAM": 0, "CHR_NAM": 5, "CHR_STR": 7, "CHR_END": 8,
        "DIR_FOR": 1, "DIR_BAK": 0, "DIR_IN": 3, "DIR_OUT": 2,
    }
    exec(compile(module, str(RUN_DEPTH_PATH), "exec"), namespace)
    return namespace["bnd_raw_contig_list"]
DIR_FOR, DIR_BAK, DIR_OUT, DIR_IN = 1, 0, 2, 3


def aligned_row(name, strand, start, end):
    return [name, 100_000, 0, end - start, strand, "chr5", 182_045_439, start, end]


class DepthPafOutputTests(unittest.TestCase):
    def make_row(self) -> list:
        return [
            "contig",
            1000,
            10,
            20,
            "+",
            "chr1",
            100_000,
            100,
            110,
            10,
            10,
            60,
            "cs:Z::10",
        ]

    def test_formats_positive_length_alignment(self) -> None:
        output = format_nonzero_depth_paf_row(self.make_row(), "10M")

        self.assertIsNotNone(output)
        self.assertTrue(output.endswith("\tcg:Z:10M"))

    def test_omits_h1437_zero_length_overlap_result(self) -> None:
        row = self.make_row()
        row[2:4] = [1416, 1416]
        row[7:9] = [72_179_591, 72_179_591]
        row[9:11] = [0, 0]
        row[-1] = "cs:Z:"

        self.assertIsNone(format_nonzero_depth_paf_row(row, ""))

    def test_omits_zero_length_in_any_paf_length_field(self) -> None:
        for field_pair in ((2, 3), (7, 8)):
            with self.subTest(field_pair=field_pair):
                row = self.make_row()
                row[field_pair[1]] = row[field_pair[0]]
                self.assertIsNone(format_nonzero_depth_paf_row(row, "10M"))

        row = self.make_row()
        row[10] = 0
        self.assertIsNone(format_nonzero_depth_paf_row(row, "10M"))

    def test_omits_empty_cigar_even_with_positive_coordinates(self) -> None:
        self.assertIsNone(format_nonzero_depth_paf_row(self.make_row(), ""))


class VirtualBridgeStrandTests(unittest.TestCase):
    def test_hcc1954_type2_446_bridge_runs_backwards(self) -> None:
        # utg007115l tail (-) exits at 113,335,109; utg115842l (-) enters at 115,142,823.
        curr = aligned_row("utg007115l", "-", 113_335_109, 113_336_742)
        nxt = aligned_row("utg115842l", "-", 115_133_474, 115_142_823)
        self.assertFalse(virtual_bridge_follows_strand(curr, curr, DIR_FOR, nxt, nxt, DIR_FOR))

    def test_hcc1954_type2_447_bridge_follows_minus_strand(self) -> None:
        curr = aligned_row("utg007115l", "-", 113_335_109, 113_336_742)
        nxt = aligned_row("utg115842l", "+", 110_820_022, 110_823_544)
        self.assertTrue(virtual_bridge_follows_strand(curr, curr, DIR_FOR, nxt, nxt, DIR_BAK))

    def test_plus_strand_gap_and_reverse_traversal(self) -> None:
        left = aligned_row("a", "+", 1_000, 2_000)
        right = aligned_row("b", "+", 5_000, 6_000)
        self.assertTrue(virtual_bridge_follows_strand(left, left, DIR_FOR, right, right, DIR_FOR))
        self.assertFalse(virtual_bridge_follows_strand(right, right, DIR_FOR, left, left, DIR_FOR))
        # Walking both rows backwards turns them into a valid minus-strand continuation.
        self.assertTrue(virtual_bridge_follows_strand(right, right, DIR_BAK, left, left, DIR_BAK))

    def test_mismatched_walked_strands_are_rejected(self) -> None:
        left = aligned_row("a", "+", 1_000, 2_000)
        right = aligned_row("b", "-", 5_000, 6_000)
        self.assertFalse(virtual_bridge_follows_strand(left, left, DIR_FOR, right, right, DIR_FOR))

    def test_graph_traversal_labels_are_not_checked(self) -> None:
        left = aligned_row("a", "+", 5_000, 6_000)
        right = aligned_row("b", "+", 1_000, 2_000)
        self.assertTrue(virtual_bridge_follows_strand(left, left, DIR_OUT, right, right, DIR_FOR))


class EcdnaCircuitSegmentTests(unittest.TestCase):
    BND_TYPE, CTG_IN_TYPE = 0, 1

    def test_crossed_circuit_keeps_the_legacy_layout(self) -> None:
        # Stage 01 stores a crossed circuit as (s1, e1, e2, s2).
        s1, e1, s2, e2 = 8, 9, 1472, 1473
        self.assertEqual(
            ecdna_circuit_segments((s1, e1, e2, s2)),
            [
                (self.CTG_IN_TYPE, ((DIR_FOR, s1), (DIR_FOR, e1))),
                (self.BND_TYPE, ((DIR_FOR, e1), (DIR_BAK, e2))),
                (self.CTG_IN_TYPE, ((DIR_BAK, e2), (DIR_BAK, s2))),
                (self.BND_TYPE, ((DIR_BAK, s2), (DIR_FOR, s1))),
            ],
        )

    def test_tandem_circuit_walks_the_second_unitig_forward(self) -> None:
        # H1437 circuit 28: utg019208l (476, 477) then utg041993l (890, 892).
        self.assertEqual(
            ecdna_circuit_segments((476, 477, 890, 892)),
            [
                (self.CTG_IN_TYPE, ((DIR_FOR, 476), (DIR_FOR, 477))),
                (self.BND_TYPE, ((DIR_FOR, 477), (DIR_FOR, 890))),
                (self.CTG_IN_TYPE, ((DIR_FOR, 890), (DIR_FOR, 892))),
                (self.BND_TYPE, ((DIR_FOR, 892), (DIR_FOR, 476))),
            ],
        )

    def test_tandem_bridges_follow_the_walked_strand(self) -> None:
        rows = {
            476: aligned_row("utg019208l", "+", 99_318_297, 99_382_654),
            477: aligned_row("utg019208l", "-", 99_314_518, 99_315_331),
            890: aligned_row("utg041993l", "-", 99_296_135, 99_313_461),
            892: aligned_row("utg041993l", "+", 99_318_896, 99_339_418),
        }
        segments = ecdna_circuit_segments((476, 477, 890, 892))
        (curr_type, curr), (next_type, nxt) = segments[1][1]
        self.assertTrue(virtual_bridge_follows_strand(
            rows[curr], rows[curr], curr_type, rows[nxt], rows[nxt], next_type,
        ))
        # The old fixed crossed layout bridged 477 into 892 walked backwards.
        self.assertFalse(virtual_bridge_follows_strand(
            rows[477], rows[477], DIR_FOR, rows[892], rows[892], DIR_BAK,
        ))


class BndRawContigListTests(unittest.TestCase):
    def make_graph(self, *edges):
        import networkx as nx

        graph = nx.DiGraph()
        graph.add_weighted_edges_from(edges)
        return graph

    def test_overlapping_endpoints_skip_the_graph_detour(self) -> None:
        # HCC1954 chr5: utg007115l (254) and utg023986l (744) overlap at
        # 262,461-277,788; G offers a strandless detour through the contained
        # utg110820l (4668), which stage 21 used to bridge with a 227 kb virtual.
        rows = {
            254: aligned_row("utg007115l", "+", 19_331, 277_788),
            744: aligned_row("utg023986l", "+", 262_461, 299_276),
            4668: aligned_row("utg110820l", "+", 20_563, 49_779),
        }
        graph = self.make_graph(
            ((DIR_OUT, 254), (DIR_FOR, 4668), 0),
            ((DIR_FOR, 4668), (DIR_IN, 744), 212_682),
        )
        build = load_bnd_raw_contig_list(rows, graph)
        self.assertEqual(
            build(((DIR_BAK, 254), (DIR_BAK, 744))),
            [(DIR_BAK, 254), (DIR_BAK, 744)],
        )
        self.assertEqual(
            build(((DIR_FOR, 744), (DIR_FOR, 254))),
            [(DIR_FOR, 744), (DIR_FOR, 254)],
        )

    def test_disjoint_endpoints_are_bridged_through_the_graph(self) -> None:
        rows = {
            1: aligned_row("a", "+", 1_000, 2_000),
            2: aligned_row("b", "+", 3_000, 4_000),
            3: aligned_row("c", "+", 5_000, 6_000),
        }
        graph = self.make_graph(
            ((DIR_OUT, 1), (DIR_FOR, 2), 1_000),
            ((DIR_FOR, 2), (DIR_IN, 3), 1_000),
        )
        build = load_bnd_raw_contig_list(rows, graph)
        self.assertEqual(
            build(((DIR_FOR, 1), (DIR_FOR, 3))),
            [(DIR_FOR, 1), (DIR_FOR, 2), (DIR_FOR, 3)],
        )

    def test_endpoints_without_graph_path_are_joined_directly(self) -> None:
        rows = {
            1: aligned_row("a", "+", 1_000, 2_000),
            3: aligned_row("c", "+", 5_000, 6_000),
        }
        # Stage 21 adds IN/OUT nodes for every endpoint even without edges.
        graph = self.make_graph()
        graph.add_nodes_from([(DIR_OUT, 1), (DIR_IN, 3)])
        build = load_bnd_raw_contig_list(rows, graph)
        self.assertEqual(
            build(((DIR_FOR, 1), (DIR_FOR, 3))),
            [(DIR_FOR, 1), (DIR_FOR, 3)],
        )


if __name__ == "__main__":
    unittest.main()
