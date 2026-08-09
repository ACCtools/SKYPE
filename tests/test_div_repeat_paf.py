from __future__ import annotations

import ast
import tempfile
import unittest
from pathlib import Path


BUILD_GRAPH_PATH = (
    Path(__file__).resolve().parents[1] / "02_Build_Breakend_Graph_Limited.py"
)


def load_div_repeat_helpers():
    """Load div-repeat helpers without executing the stage-02 entry point."""
    function_names = {
        "is_stepwise_nonoverlapping",
        "is_alt_contained_in_segments",
        "max_overlap",
        "get_qry_cord_data",
        "div_repeat_paf",
    }
    tree = ast.parse(BUILD_GRAPH_PATH.read_text(encoding="utf-8"))
    functions = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name in function_names
    ]
    namespace = {
        "CTG_NAM": 0,
        "CTG_ENDND": 12,
        "CTG_GLOBALIDX": 21,
        "MAX_OVERLAP_SCORE": 3,
    }
    module = ast.Module(body=functions, type_ignores=[])
    exec(compile(module, str(BUILD_GRAPH_PATH), "exec"), namespace)
    return namespace


HELPERS = load_div_repeat_helpers()
div_repeat_paf = HELPERS["div_repeat_paf"]


def paf_row(
    query_start: int,
    query_end: int,
    chromosome: str,
    mapq: int,
    *tags: str,
) -> str:
    fields = [
        "utg013626l",
        "43480",
        str(query_start),
        str(query_end),
        "-",
        chromosome,
        "250000000",
        "100",
        str(100 + query_end - query_start),
        str(query_end - query_start),
        str(query_end - query_start),
        str(mapq),
        *tags,
    ]
    return "\t".join(fields) + "\n"


def ppc_node(global_index: str, end_index: int) -> list:
    row = ["0"] * 22
    row[0] = "utg013626l"
    row[12] = end_index
    row[21] = global_index
    return row


class DivRepeatPafTests(unittest.TestCase):
    def test_strong_chain_is_independent_of_original_paf_row_order(self):
        # Reproduce utg013626l: the two selected primary records form a valid
        # query-ordered chain, but the original PAF stores the right half first.
        # A secondary alignment duplicates that right half, yielding overlap 3.
        original_text = "".join([
            paf_row(21654, 43480, "chr2", 43, "tp:A:P"),
            paf_row(0, 21658, "chrX", 60, "tp:A:P"),
            paf_row(21654, 43480, "chr2", 0, "tp:A:S"),
        ])
        aligned_text = "".join([
            paf_row(0, 21654, "chrX", 60, "tp:A:P", "xi:Z:P_1"),
            paf_row(21654, 43480, "chr2", 43, "tp:A:P", "xi:Z:P_0"),
        ])

        with tempfile.TemporaryDirectory() as temporary:
            temporary_path = Path(temporary)
            original_path = temporary_path / "original.paf"
            aligned_path = temporary_path / "aligned.paf"
            original_path.write_text(original_text, encoding="utf-8")
            aligned_path.write_text(aligned_text, encoding="utf-8")

            HELPERS["ori_ctg_name_data"] = [
                ["utg013626l", "utg013626l"]
            ]
            excluded = div_repeat_paf(
                [str(original_path)],
                [str(aligned_path)],
                [ppc_node("0.0", 1), ppc_node("0.1", 1)],
            )

        self.assertNotIn("utg013626l", excluded)


if __name__ == "__main__":
    unittest.main()
