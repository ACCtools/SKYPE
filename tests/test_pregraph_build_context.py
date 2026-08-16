from __future__ import annotations

import ast
import unittest
from dataclasses import dataclass, fields
from enum import Enum
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, Mock


BUILD_GRAPH_PATH = (
    Path(__file__).resolve().parents[1] / "02_Build_Breakend_Graph_Limited.py"
)


def parse_stage02() -> ast.Module:
    return ast.parse(BUILD_GRAPH_PATH.read_text(encoding="utf-8"))


def top_level_definition(name: str):
    return next(
        node
        for node in parse_stage02().body
        if isinstance(node, (ast.ClassDef, ast.FunctionDef))
        and node.name == name
    )


def load_definitions(*definition_names: str, **extra_globals):
    """Load selected contracts without executing stage 02's CLI body."""

    requested = set(definition_names)
    definitions = [
        node
        for node in parse_stage02().body
        if isinstance(node, (ast.ClassDef, ast.FunctionDef))
        and node.name in requested
    ]
    found = {node.name for node in definitions}
    missing = requested - found
    if missing:
        raise AssertionError(
            f"Missing pregraph build contracts: {sorted(missing)}"
        )

    namespace = {
        "dataclass": dataclass,
        "Enum": Enum,
        **extra_globals,
    }
    module = ast.Module(body=definitions, type_ignores=[])
    exec(compile(module, str(BUILD_GRAPH_PATH), "exec"), namespace)
    return namespace


def direct_name_calls(function_name: str) -> list[tuple[str, int]]:
    function = top_level_definition(function_name)
    calls = sorted(
        (node for node in ast.walk(function) if isinstance(node, ast.Call)),
        key=lambda node: (node.lineno, node.col_offset),
    )
    return [
        (call.func.id, call.lineno)
        for call in calls
        if isinstance(call.func, ast.Name)
    ]


def loaded_names(function_name: str) -> set[str]:
    function = top_level_definition(function_name)
    return {
        node.id
        for node in ast.walk(function)
        if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load)
    }


def function_arg_names(function_name: str) -> list[str]:
    function = top_level_definition(function_name)
    return [argument.arg for argument in function.args.args]


class PregraphContextContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.contracts = load_definitions(
            "PregraphSourceMode",
            "NCloseSourceConfig",
            "PregraphBuildContext",
            "ContigPreprocessResult",
            "PafSourceKind",
            "PafPreprocessPolicy",
        )

    def test_context_has_one_explicit_immutable_build_contract(self) -> None:
        context_type = self.contracts["PregraphBuildContext"]

        self.assertEqual(
            [field.name for field in fields(context_type)],
            [
                "source",
                "ori_ctg_name_data",
                "prefix",
                "preprocessed_paf_path",
                "reference_fai_path",
                "telomere_bed_path",
                "repeat_bed_path",
                "censat_bed_path",
                "main_stat_path",
                "asm2cov",
                "disable_alt_ctg_simple",
            ],
        )
        self.assertTrue(
            context_type.__dataclass_params__.frozen,
            "A build attempt must not mutate its source/configuration context",
        )

    def test_source_paths_exist_only_inside_source_config(self) -> None:
        context_fields = {
            field.name
            for field in fields(self.contracts["PregraphBuildContext"])
        }

        self.assertTrue(
            {
                "paf_file_paths",
                "original_paf_paths",
                "is_unitig_reduced",
                "secondary_candidate_paf",
            }.isdisjoint(context_fields)
        )

    def test_contig_result_carries_only_data_discovered_during_preprocessing(
        self,
    ) -> None:
        result_type = self.contracts["ContigPreprocessResult"]

        self.assertEqual(
            [field.name for field in fields(result_type)],
            ["telo_coverage", "depth_df", "no_chrY"],
        )

    def test_paf_policy_contract_distinguishes_primary_and_secondary(self) -> None:
        policy_type = self.contracts["PafPreprocessPolicy"]
        source_kind = self.contracts["PafSourceKind"]

        self.assertEqual(
            [field.name for field in fields(policy_type)],
            ["kind", "source_index", "global_index_prefix", "log_label"],
        )
        self.assertEqual(
            list(source_kind.__members__), ["PRIMARY", "SECONDARY"]
        )

        primary = policy_type(source_kind.PRIMARY, 0, "0", "primary PAF")
        secondary = policy_type(
            source_kind.SECONDARY,
            1,
            "1",
            "alternative PAF",
        )
        self.assertEqual(
            (
                primary.source_index,
                primary.global_index_prefix,
                primary.log_label,
            ),
            (0, "0", "primary PAF"),
        )
        self.assertEqual(
            (
                secondary.source_index,
                secondary.global_index_prefix,
                secondary.log_label,
            ),
            (1, "1", "alternative PAF"),
        )


class PregraphSourceIsolationTests(unittest.TestCase):
    SOURCE_GLOBALS = {
        "PAF_FILE_PATH",
        "ORIGINAL_PAF_LOC_LIST",
        "ori_ctg_name_data",
        "is_unitig_reduced",
    }
    CONTEXT_CONFIG_GLOBALS = {
        "PREFIX",
        "PREPROCESSED_PAF_FILE_PATH",
        "CHROMOSOME_INFO_FILE_PATH",
        "TELOMERE_INFO_FILE_PATH",
        "REPEAT_INFO_FILE_PATH",
        "CENSAT_PATH",
        "main_stat_loc",
        "asm2cov",
    }

    def test_pregraph_build_does_not_activate_source_before_success(self) -> None:
        calls = direct_name_calls("build_pregraph_state")
        called_names = [name for name, _ in calls]

        self.assertNotIn("activate_nclose_source", called_names)
        self.assertNotIn("activate_pregraph_state", called_names)
        self.assertIn("PregraphBuildContext", called_names)
        self.assertIn("contig_preprocessing_00", called_names)
        self.assertIn("nclose_calc", called_names)
        self.assertIn("PregraphState", called_names)

        positions = {
            name: next(line for called, line in calls if called == name)
            for name in (
                "PregraphBuildContext",
                "contig_preprocessing_00",
                "nclose_calc",
                "PregraphState",
            )
        }
        self.assertLess(
            positions["PregraphBuildContext"],
            positions["contig_preprocessing_00"],
        )
        self.assertLess(
            positions["contig_preprocessing_00"], positions["nclose_calc"]
        )
        self.assertLess(positions["nclose_calc"], positions["PregraphState"])

    def test_build_workers_use_context_instead_of_source_globals(self) -> None:
        self.assertEqual(function_arg_names("contig_preprocessing_00"), ["context"])
        self.assertEqual(
            function_arg_names("nclose_calc"),
            ["context", "preprocess_result"],
        )

        for function_name in ("contig_preprocessing_00", "nclose_calc"):
            with self.subTest(function_name=function_name):
                leaked = loaded_names(function_name) & self.SOURCE_GLOBALS
                self.assertEqual(
                    leaked,
                    set(),
                    f"{function_name} leaked selected source globals: {leaked}",
                )

        contig_leaks = loaded_names("contig_preprocessing_00") & (
            self.CONTEXT_CONFIG_GLOBALS | {"args"}
        )
        self.assertEqual(
            contig_leaks,
            set(),
            "contig preprocessing must read paths/options from its context",
        )

        nclose_leaks = loaded_names("nclose_calc") & {
            "PAF_FILE_PATH",
            "ORIGINAL_PAF_LOC_LIST",
            "PREPROCESSED_PAF_FILE_PATH",
            "PREFIX",
            "asm2cov",
        }
        self.assertEqual(
            nclose_leaks,
            set(),
            "NClose discovery must read source/build paths from its context",
        )

    def test_failed_build_cannot_publish_partial_source_state(self) -> None:
        namespace = load_definitions("build_pregraph_state")
        build_pregraph_state = namespace["build_pregraph_state"]
        source = SimpleNamespace(
            mode=SimpleNamespace(name="CONFIGURED_PAF"),
            paf_file_paths=("primary.paf",),
            original_paf_paths=("primary.original.paf",),
            secondary_candidate_paf=None,
            is_unitig_reduced=False,
        )
        activation = Mock(
            side_effect=AssertionError("partial build source was activated")
        )
        context_builder = Mock(return_value=MagicMock())
        original_name_loader = Mock(return_value=[["ctg"]])

        namespace.update(
            {
                "PregraphSourceMode": SimpleNamespace(
                    VCF=object(), CONFIGURED_PAF=source.mode
                ),
                "PregraphBuildContext": context_builder,
                "get_ori_ctg_name_data": original_name_loader,
                "contig_preprocessing_00": Mock(
                    side_effect=RuntimeError("preprocessing failed")
                ),
                "nclose_calc": Mock(),
                "build_vcf_mode_inputs": Mock(),
                "PregraphState": Mock(),
                "activate_nclose_source": activation,
                "activate_pregraph_state": activation,
                "PREFIX": "/output",
                "PREPROCESSED_PAF_FILE_PATH": "/output/preprocessed.paf",
                "CHROMOSOME_INFO_FILE_PATH": "/ref/reference.fai",
                "TELOMERE_INFO_FILE_PATH": "/ref/telomere.bed",
                "REPEAT_INFO_FILE_PATH": "/ref/repeat.bed",
                "CENSAT_PATH": "/ref/censat.bed",
                "main_stat_loc": "/output/depth.tsv.gz",
                "asm2cov": {"primary": 30},
                "args": SimpleNamespace(
                    censat_bed_path="/ref/censat.bed",
                    disable_alt_ctg_simple=False,
                ),
            }
        )

        with self.assertRaisesRegex(RuntimeError, "preprocessing failed"):
            build_pregraph_state(source)

        activation.assert_not_called()
        namespace["nclose_calc"].assert_not_called()
        namespace["PregraphState"].assert_not_called()
        original_name_loader.assert_called_once_with(["primary.paf"])
        context_builder.assert_called_once_with(
            source=source,
            ori_ctg_name_data=[["ctg"]],
            prefix="/output",
            preprocessed_paf_path="/output/preprocessed.paf",
            reference_fai_path="/ref/reference.fai",
            telomere_bed_path="/ref/telomere.bed",
            repeat_bed_path="/ref/repeat.bed",
            censat_bed_path="/ref/censat.bed",
            main_stat_path="/output/depth.tsv.gz",
            asm2cov={"primary": 30},
            disable_alt_ctg_simple=False,
        )


class ContigPreprocessingOrchestrationTests(unittest.TestCase):
    REQUIRED_STAGES = [
        "preprocess_paf_source",
        "build_split_contigs",
        "build_simple_alt_contigs",
        "finalize_preprocessed_contigs",
    ]

    @staticmethod
    def run_mocked_orchestrator(source):
        namespace = load_definitions(
            "PafSourceKind",
            "PafPreprocessPolicy",
            "ContigPreprocessResult",
            "contig_preprocessing_00",
        )
        function = namespace["contig_preprocessing_00"]
        source_kind = namespace["PafSourceKind"]
        events = []
        resources = SimpleNamespace(depth_df="depth", no_chrY=True)
        primary = SimpleNamespace(
            contigs=["primary-0", "primary-1"],
            kept_primary_names={"primary-kept"},
            excluded_telomere_origins={(0, 7)},
            original_node_count=11,
        )
        secondary = SimpleNamespace(
            contigs=["secondary-0"],
            kept_primary_names=set(),
            excluded_telomere_origins={(1, 3)},
            original_node_count=5,
        )

        def load_resources(context):
            events.append(("resources", context))
            return resources

        def preprocess_source(
            context,
            observed_resources,
            paf_path,
            policy,
            index_offset=0,
        ):
            events.append(
                (
                    "source",
                    policy.kind,
                    paf_path,
                    index_offset,
                    policy,
                )
            )
            if policy.kind == source_kind.PRIMARY:
                return primary
            return secondary

        def split_contigs(
            context,
            observed_resources,
            contigs,
            excluded_origins,
            original_node_count,
        ):
            events.append(
                (
                    "split",
                    tuple(contigs),
                    frozenset(excluded_origins),
                    original_node_count,
                )
            )
            return list(contigs) + ["split"], "telomere-bounds"

        def simple_alt(
            context,
            observed_resources,
            contigs,
            primary_kept_names,
        ):
            events.append(
                (
                    "simple-alt",
                    tuple(contigs),
                    frozenset(primary_kept_names),
                )
            )
            return list(contigs) + ["simple-alt"]

        def finalize(
            context,
            observed_resources,
            contigs,
            telo_bound_dict,
            excluded_origins,
        ):
            events.append(
                (
                    "finalize",
                    tuple(contigs),
                    telo_bound_dict,
                    frozenset(excluded_origins),
                )
            )
            return "telomere-coverage"

        namespace.update(
            {
                "load_contig_preprocess_resources": Mock(
                    side_effect=load_resources
                ),
                "preprocess_paf_source": Mock(side_effect=preprocess_source),
                "build_split_contigs": Mock(side_effect=split_contigs),
                "build_simple_alt_contigs": Mock(side_effect=simple_alt),
                "finalize_preprocessed_contigs": Mock(side_effect=finalize),
            }
        )
        context = SimpleNamespace(source=source)
        result = function(context)
        return namespace, result, events

    def test_dual_source_policies_and_stage_order(self) -> None:
        source = SimpleNamespace(
            paf_file_paths=("primary.paf", "unitig.paf"),
            secondary_candidate_paf="unitig.paf",
            is_unitig_reduced=False,
        )

        namespace, result, events = self.run_mocked_orchestrator(source)
        source_kind = namespace["PafSourceKind"]

        self.assertEqual(
            [event[0] for event in events],
            [
                "resources",
                "source",
                "source",
                "split",
                "simple-alt",
                "finalize",
            ],
        )
        primary_call, secondary_call = events[1:3]
        self.assertEqual(
            primary_call[1:4],
            (source_kind.PRIMARY, "primary.paf", 0),
        )
        self.assertEqual(
            secondary_call[1:4],
            (source_kind.SECONDARY, "unitig.paf", 2),
        )
        self.assertEqual(
            (
                primary_call[4].source_index,
                primary_call[4].global_index_prefix,
                primary_call[4].log_label,
            ),
            (0, "0", "primary PAF"),
        )
        self.assertEqual(
            (
                secondary_call[4].source_index,
                secondary_call[4].global_index_prefix,
                secondary_call[4].log_label,
            ),
            (1, "1", "alternative PAF"),
        )
        self.assertEqual(
            events[3],
            (
                "split",
                ("primary-0", "primary-1", "secondary-0"),
                frozenset({(0, 7), (1, 3)}),
                16,
            ),
        )
        self.assertEqual(
            events[4][-1],
            frozenset({"primary-kept"}),
            "simple-alt exclusion must come only from kept primary names",
        )
        self.assertEqual(result.telo_coverage, "telomere-coverage")
        self.assertEqual(result.depth_df, "depth")
        self.assertTrue(result.no_chrY)

    def test_singleton_and_primary_retry_never_run_secondary_policy(self) -> None:
        sources = {
            "singleton": SimpleNamespace(
                paf_file_paths=("primary.paf",),
                secondary_candidate_paf=None,
                is_unitig_reduced=False,
            ),
            "primary-only-retry": SimpleNamespace(
                paf_file_paths=("primary.paf", "primary.paf"),
                secondary_candidate_paf=None,
                is_unitig_reduced=True,
            ),
        }

        for label, source in sources.items():
            with self.subTest(label=label):
                namespace, _, events = self.run_mocked_orchestrator(source)
                source_events = [
                    event for event in events if event[0] == "source"
                ]
                self.assertEqual(len(source_events), 1)
                self.assertIs(
                    source_events[0][1], namespace["PafSourceKind"].PRIMARY
                )
                self.assertEqual(source_events[0][2], "primary.paf")
                self.assertEqual(
                    [event[0] for event in events],
                    [
                        "resources",
                        "source",
                        "split",
                        "simple-alt",
                        "finalize",
                    ],
                )

    def test_stage_calls_are_present_in_dependency_order(self) -> None:
        calls = direct_name_calls("contig_preprocessing_00")

        stage_lines = {}
        for stage_name in self.REQUIRED_STAGES:
            matching_lines = [line for name, line in calls if name == stage_name]
            self.assertTrue(matching_lines, f"missing stage {stage_name}")
            stage_lines[stage_name] = matching_lines

        preprocess_lines = stage_lines["preprocess_paf_source"]
        self.assertEqual(
            len(preprocess_lines),
            2,
            "primary and optional secondary must share one source helper",
        )
        self.assertLess(preprocess_lines[0], preprocess_lines[1])
        self.assertLess(
            preprocess_lines[1], stage_lines["build_split_contigs"][0]
        )
        self.assertLess(
            stage_lines["build_split_contigs"][0],
            stage_lines["build_simple_alt_contigs"][0],
        )
        self.assertLess(
            stage_lines["build_simple_alt_contigs"][0],
            stage_lines["finalize_preprocessed_contigs"][0],
        )

    def test_secondary_source_is_gated_by_explicit_source_config(self) -> None:
        function = top_level_definition("contig_preprocessing_00")
        preprocess_calls = [
            node
            for node in ast.walk(function)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "preprocess_paf_source"
        ]
        preprocess_calls.sort(key=lambda node: (node.lineno, node.col_offset))
        self.assertEqual(
            len(preprocess_calls),
            2,
            "primary and optional secondary require two explicit call sites",
        )
        secondary_call = preprocess_calls[1]

        guarding_ifs = [
            node
            for node in ast.walk(function)
            if isinstance(node, ast.If)
            and node.lineno < secondary_call.lineno <= node.end_lineno
        ]
        self.assertTrue(
            any(
                "secondary_candidate_paf" in {
                    child.attr
                    for child in ast.walk(guard.test)
                    if isinstance(child, ast.Attribute)
                }
                for guard in guarding_ifs
            ),
            "secondary preprocessing must be gated by context.source.secondary_candidate_paf",
        )


class PrimaryRetryCompatibilityTests(unittest.TestCase):
    def test_primary_only_retry_keeps_legacy_duplicated_path_layout(self) -> None:
        helpers = load_definitions(
            "PregraphSourceMode",
            "NCloseSourceConfig",
            "resolve_pregraph_source",
            VCF_SYNTHETIC_PAF_NAME="vcf_synthetic.paf",
        )
        source_mode = helpers["PregraphSourceMode"]
        resolve_source = helpers["resolve_pregraph_source"]

        source = resolve_source(
            source_mode.PRIMARY_ONLY_RETRY,
            "/inputs/primary.paf",
            "/inputs/unitig.paf",
            ["/originals/primary.paf", "/originals/unitig.paf"],
            "/output",
        )

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


if __name__ == "__main__":
    unittest.main()
