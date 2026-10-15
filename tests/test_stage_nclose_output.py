import ast
import copy
import csv
import pickle
import tempfile
import unittest
from pathlib import Path

import vcfpy

from native_bnd import node_pair_endpoints
from native_structure_output import display_events, write_native_bed, write_native_vcf
from nclose_tracking import _build_bnd_event, carrier_columns, save_filter_status
from structure_nclose import (
    MODEL_FILE, StructureWeights, add_virtual_structure,
    build_stage_weight_contexts, load_or_build_structure_model, save_structure_model,
    stage_report_values, write_structure_reports,
)
from test_native_bnd import export_namespace, node
from test_structure_nclose import PAIR, OTHER, fixture, write_internal_chromosome_fixture


def stage_status(filter_reasons=None, cluster_reasons=None, base_reasons=None):
    """Stage-25 source reasons as stage 31 reads them."""
    filter_reasons = filter_reasons or {}
    return {"stages": {
        "base": dict(reasons=base_reasons or {}),
        "filter": dict(reasons=filter_reasons),
        "cluster": dict(reasons=filter_reasons if cluster_reasons is None else cluster_reasons),
    }}


def alias_fixture():
    nodes, model = fixture([dict(kind="PATH", nclose_counts={PAIR: 1}),
                            dict(kind="PATH", nclose_counts={OTHER: 1})])
    nodes[1] = node("u1", 100, 200, "-", "chr2", 400, 500)
    nodes[2] = node("u2", 0, 100, "+", "chr1", 100, 200)
    nodes[3] = node("u2", 100, 200, "-", "chr2", 400, 500)
    for key in (PAIR, OTHER):
        event = _build_bnd_event(key, nodes, "compressed_nclose_nodes_list.txt")
        event["nclose_id"] = model["ncloses"][key]["nclose_id"]
        event["endpoints"], event["contig_names"] = node_pair_endpoints(nodes, key)
        model["ncloses"][key] = event
    return nodes, model


def read_tsv(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def export(prefix, context, nodes, suffix=""):
    namespace = export_namespace()

    class Ratios:
        def pair(self, *args):
            return {namespace["BP_STEP_DEPTH_RATIO_B"]: [None, None],
                    namespace["BP_STEP_DEPTH_RATIO_PREDICT_B"]: [None, None]}

    vcf = prefix / f"SV_call_result{suffix}.vcf"
    write_native_vcf(context, {row[5]: 1_000_000 for row in nodes}, vcf, Ratios(),
                     lambda lengths: namespace["build_vcf_header"](lengths, suffix),
                     namespace["write_bnd_vcf_pair"],
                     namespace["write_symbolic_vcf_record"], nodes=nodes)
    write_native_bed(prefix, display_events(context, nodes, {}), suffix)
    with vcfpy.Reader.from_path(str(vcf)) as reader:
        records = list(reader)
    return records


class CanonicalStageStatusTests(unittest.TestCase):
    def test_live_alias_overrides_removed_representative_in_report_vcf_and_bed(self):
        nodes, model = alias_fixture()
        raw = StructureWeights(copy.deepcopy(model), [10., 10.], 10.)
        filtered = StructureWeights(copy.deepcopy(model), [0., 20.], 10.)
        status = stage_status({PAIR: "FILTERED_24_GREEDY_DEPTH"})
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary)
            write_structure_reports(prefix, raw, status, {"filter": filtered})
            rows = read_tsv(prefix / "nclose_report.tsv")
            row = next(r for r in rows if r["nclose_key"] == "0:1")
            self.assertEqual((row["nclose_filter"], row["nclose_filter_reason"]), ("2.0", "PASS"))
            records = export(prefix, filtered, nodes, "_filter")
            self.assertTrue(any(row["nclose_id"] in r.INFO["NCLOSE_IDS"] for r in records))
            bed = read_tsv(prefix / "SKYPE_result_filter.bed")
            self.assertTrue(any(r["nclose_id"] == row["nclose_id"] and float(r["weight (N)"]) == 2. for r in bed))
        self.assertEqual(status["stages"]["filter"]["reasons"][PAIR], "FILTERED_24_GREEDY_DEPTH")

    def test_natural_zero_alias_is_still_selected(self):
        _, model = alias_fixture()
        context = StructureWeights(model, [0., 0.], 10.)
        status = stage_status({PAIR: "FILTERED_24_GREEDY_DEPTH"})
        self.assertEqual(stage_report_values(status, "filter", PAIR, context), (0., "PASS"))
        self.assertEqual(stage_report_values(status, "cluster", PAIR, context), (0., "PASS"))

    def test_preprocessing_history_does_not_hide_a_retained_compound(self):
        nodes, model = fixture([dict(kind="MERGE_TYPE4", nclose_counts={PAIR: 1},
                                    event=dict(chrom="chr1", st=100, nd=500, event_type="front_jump"))])
        # Stage 25 passes PAIR through the retained compound; its own paths stay history.
        status = stage_status(base_reasons={PAIR: "FILTERED_02_NO_ELIGIBLE_PATH"})
        for weight in (0., 20.):
            with self.subTest(weight=weight):
                context = StructureWeights(copy.deepcopy(model), [weight], 10.)
                self.assertEqual(stage_report_values(status, "filter", PAIR, context), (weight / 10., "PASS"))
                self.assertEqual(stage_report_values(status, "cluster", PAIR, context), (weight / 10., "PASS"))
        with tempfile.TemporaryDirectory() as temporary:
            context = StructureWeights(model, [20.], 10.)
            write_structure_reports(temporary, context, status, {"filter": context, "cluster": context})
            row = read_tsv(Path(temporary) / "nclose_report.tsv")[0]
            self.assertEqual(row["preprocessing_reason"], "FILTERED_02_NO_ELIGIBLE_PATH")
            self.assertEqual(row["nclose_cluster_reason"], "PASS")
            export(Path(temporary), context, nodes, "_cluster")

    def test_all_alias_carriers_removed_reports_the_filter_reason(self):
        _, model = alias_fixture()
        context = StructureWeights(model, [0., 0.], 10.)
        status = stage_status({PAIR: "FILTERED_24_GREEDY_DEPTH", OTHER: "FILTERED_24_COFILTERED_PATH"})
        for stage in ("filter", "cluster"):
            self.assertEqual(stage_report_values(status, stage, PAIR, context), (-1, "FILTERED_24_GREEDY_DEPTH"))

    def test_floor_cut_is_distinct_from_an_nnls_zero_and_live_alias_support(self):
        _, model = alias_fixture()
        status = stage_status(cluster_reasons={PAIR: "FILTERED_24_CLUSTER_MIN_WEIGHT"})
        for survivor, expected in ((0., (-1, "FILTERED_24_CLUSTER_MIN_WEIGHT")),
                                   (20., (2., "PASS"))):
            with self.subTest(survivor=survivor):
                context = StructureWeights(copy.deepcopy(model), [0., survivor], 10.)
                self.assertEqual(stage_report_values(status, "cluster", PAIR, context), expected)

    def test_qualified_virtual_support_survives_removed_matrix_carriers(self):
        nodes, model = fixture([dict(kind="PATH", nclose_counts={PAIR: 1})])
        record = dict(pair_id=1, nclose_key_a=PAIR, nclose_key_b=OTHER, layout_a={}, layout_b={})
        add_virtual_structure(model, record, 5., [("chr1", 100, 500)], nodes)
        context = StructureWeights(model, [0.], 10.)
        status = stage_status({PAIR: "FILTERED_24_GREEDY_DEPTH"})
        self.assertEqual(stage_report_values(status, "filter", PAIR, context), (.5, "PASS"))

    def test_alias_sources_removed_at_different_stages_report_the_latest_cause(self):
        _, model = alias_fixture()
        context = StructureWeights(model, [0., 0.], 10.)
        status = stage_status({PAIR: "FILTERED_02_NO_ELIGIBLE_PATH", OTHER: "FILTERED_24_COFILTERED_PATH"},
                              {PAIR: "FILTERED_02_NO_ELIGIBLE_PATH", OTHER: "FILTERED_24_CLUSTER_SELECTION"})
        self.assertEqual(stage_report_values(status, "filter", PAIR, context), (-1, "FILTERED_24_COFILTERED_PATH"))
        self.assertEqual(stage_report_values(status, "cluster", PAIR, context),
                         (-1, "FILTERED_24_CLUSTER_SELECTION"))

    def test_stage25_records_structure_carriers_including_restored_sources(self):
        source = Path(__file__).resolve().parents[1] / "25_cluster_weight.py"
        names = {"record_structure_stage", "filter_reason", "cluster_reason"}
        functions = [node for node in ast.parse(source.read_text()).body
                     if isinstance(node, ast.FunctionDef) and node.name in names]
        restored, greedy, unselected, early = (10, 11), (12, 13), (14, 15), (16, 17)
        # Column 0 is a retained compound carrying PAIR, whose own paths were excluded.
        base_reasons = {PAIR: "FILTERED_02_NO_ELIGIBLE_PATH"}
        status = {"version": 1, "event_keys": [PAIR, greedy, unselected, early],
                  "stages": {"base": dict(active_columns=[0, 1, 2, 3], reasons=base_reasons)}}
        with tempfile.TemporaryDirectory() as temporary:
            namespace = dict(
                PREFIX=temporary, nclose_filter_status=status, carrier_columns=carrier_columns,
                save_filter_status=save_filter_status,
                provenance_usage=[{PAIR: 1}, {restored: 1}, {greedy: 1}, {unselected: 1}, {early: 1}],
                tracked_keys=[PAIR, greedy, unselected, early, restored],
                base_active_columns={0, 1, 2, 3}, base_reasons=base_reasons,
                filter_active_column_set={0, 1, 3}, not_essential_nclose={greedy},
                selected_columns={0, 1}, cluster_thresholded_columns={1},
                cluster_positive_columns=set(), cluster_active_columns={0})
            exec(compile(ast.Module(body=functions, type_ignores=[]), str(source), "exec"), namespace)
            namespace["record_structure_stage"]("filter", {0, 1, 3}, namespace["filter_reason"])
            namespace["filter_reasons"] = status["stages"]["filter"]["reasons"]
            namespace["record_structure_stage"]("cluster", {0}, namespace["cluster_reason"])
            saved = pickle.loads((Path(temporary) / "nclose_filter_status.pkl").read_bytes())
        removed_by_filter = {greedy: "FILTERED_24_GREEDY_DEPTH", early: "FILTERED_02_COFILTERED_PATH"}
        self.assertEqual(saved["stages"]["filter"]["reasons"], removed_by_filter)
        self.assertEqual(saved["stages"]["cluster"]["reasons"], {
            **removed_by_filter, unselected: "FILTERED_24_CLUSTER_SELECTION",
            restored: "FILTERED_24_CLUSTER_MIN_WEIGHT"})


def virtual_fixture():
    nodes, model = fixture([dict(kind="PATH", nclose_counts={PAIR: 1})])
    records = {}
    for pair_id in (1, 2):
        record = {"pair_id": pair_id}
        for side in ("a", "b"):
            index = len(nodes)
            pos = 1000 * (index + 1)
            name = f"v{pair_id}{side}"
            nodes.extend([node(name, 0, 100, "+", "chr1", pos, pos + 100),
                          node(name, 100, 200, "-", "chr1", pos + 500, pos + 600)])
            record[f"nclose_key_{side}"] = (index, index + 1)
            record[f"layout_{side}"] = dict(ordered_endpoints=(
                dict(chrom="chr1", coord=pos + 100, dir="+", ctg_name=name),
                dict(chrom="chr1", coord=pos + 600, dir="-", ctg_name=name)))
        records[pair_id] = record
    return nodes, model, records


class StageVirtualIdentityTests(unittest.TestCase):
    def test_stage_only_virtuals_preserve_raw_files_and_share_stage_ids(self):
        nodes, model, records = virtual_fixture()
        pristine = copy.deepcopy(model)
        views = [("chr1", 100, 500, 5.)]
        raw_virtuals = {2: views}
        baseline = copy.deepcopy(model)
        add_virtual_structure(baseline, records[2], 5., [("chr1", 100, 500)], nodes)
        baseline = StructureWeights(baseline, [10.], 10.)
        contexts = build_stage_weight_contexts(
            model, {"": [10.], "_filter": [20.], "_cluster": [30.]},
            {"": raw_virtuals, "_filter": {1: views}, "_cluster": {1: views, 2: views}},
            records, nodes, 10.)
        self.assertEqual(model, pristine)
        raw_key = records[2]["nclose_key_a"]
        extra_key = records[1]["nclose_key_a"]
        self.assertEqual(contexts[""].ncloses[raw_key]["nclose_id"], "SKYPE.nclose.4")
        self.assertNotIn(extra_key, contexts[""].ncloses)
        self.assertEqual(contexts["_filter"].ncloses[extra_key]["nclose_id"], "SKYPE.nclose.6")
        for key in (raw_key, extra_key):
            self.assertEqual(contexts["_filter"].ncloses[key]["nclose_id"],
                             contexts["_cluster"].ncloses[key]["nclose_id"])
        self.assertEqual(contexts["_filter"].totals[raw_key], 0.)
        self.assertEqual(contexts["_cluster"].totals[raw_key], 5.)
        self.assertEqual(contexts["_filter"].totals[extra_key], 5.)

        with tempfile.TemporaryDirectory() as temporary:
            old, new = Path(temporary) / "old", Path(temporary) / "new"
            old.mkdir()
            new.mkdir()
            write_structure_reports(old, baseline, {})
            export(old, baseline, nodes)
            reports = {stage[1:]: context for stage, context in contexts.items() if stage}
            status = stage_status()
            write_structure_reports(new, contexts[""], status, reports)
            export(new, contexts[""], nodes)
            for path in old.iterdir():
                if path.name == "nclose_report.tsv":
                    omit = {"nclose_filter", "nclose_filter_reason", "nclose_cluster", "nclose_cluster_reason"}
                    self.assertEqual([{k: v for k, v in r.items() if k not in omit} for r in read_tsv(path)],
                                     [{k: v for k, v in r.items() if k not in omit} for r in read_tsv(new / path.name)])
                else:
                    self.assertEqual(path.read_bytes(), (new / path.name).read_bytes(), path.name)
            for suffix in ("_filter", "_cluster"):
                write_structure_reports(new, contexts[suffix], status, reports, suffix=suffix)
                export(new, contexts[suffix], nodes, suffix)
                rows = read_tsv(new / f"nclose_report{suffix}.tsv")
                extra = next(r for r in rows if r["nclose_id"] == "SKYPE.nclose.6")
                self.assertEqual(float(extra["nclose_cn"]), .5)
                self.assertEqual(extra[f"nclose{suffix}_reason"], "PASS")
                self.assertTrue((new / f"nclose_sources{suffix}.tsv").exists())
                self.assertTrue((new / f"structure_report{suffix}.tsv").exists())
                self.assertTrue((new / f"structure_nclose_usage{suffix}.tsv").exists())
                vcf_text = (new / f"SV_call_result{suffix}.vcf").read_text()
                self.assertIn(f"linked to structure_report{suffix}.tsv", vcf_text)
                self.assertIn(f"linked to nclose_report{suffix}.tsv", vcf_text)

    def test_old_cache_with_a_virtual_registry_is_rebuilt(self):
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary)
            nodes, _, _, _ = write_internal_chromosome_fixture(prefix)
            pristine = load_or_build_structure_model(prefix)
            record = dict(pair_id=1, nclose_key_a=(0, 1), nclose_key_b=(1, 2))
            for side, pair in (("a", (0, 1)), ("b", (1, 2))):
                a, b = (nodes[i] for i in pair)
                record[f"layout_{side}"] = dict(ordered_endpoints=(
                    dict(chrom=a[5], coord=a[8], dir=a[4], ctg_name=a[0]),
                    dict(chrom=b[5], coord=b[7], dir=b[4], ctg_name=b[0])))
            for strip_structures in (False, True):
                with self.subTest(stripped=strip_structures):
                    dirty = copy.deepcopy(pristine)
                    add_virtual_structure(dirty, record, 5., [], nodes)
                    dirty["version"] = 3
                    if strip_structures:
                        dirty["structures"] = [s for s in dirty["structures"] if s["kind"] != "VIRTUAL_INV"]
                    (prefix / MODEL_FILE).write_bytes(pickle.dumps(dirty))
                    rebuilt = load_or_build_structure_model(prefix)
                    self.assertEqual(rebuilt, pristine)

    def test_cached_topology_and_ids_survive_repeated_output_mode_changes(self):
        nodes, model, records = virtual_fixture()
        views = [("chr1", 100, 500, 5.)]
        with tempfile.TemporaryDirectory() as temporary:
            prefix = Path(temporary)
            save_structure_model(prefix, model)
            saved = (prefix / MODEL_FILE).read_bytes()
            ids = []
            for include_stages in (True, False, True):
                base = pickle.loads(saved)
                weights, virtuals = {"": [10.]}, {"": {2: views}}
                if include_stages:
                    weights.update(_filter=[20.], _cluster=[30.])
                    virtuals.update(_filter={1: views}, _cluster={1: views, 2: views})
                contexts = build_stage_weight_contexts(base, weights, virtuals, records, nodes, 10.)
                ids.append(contexts[""].ncloses[records[2]["nclose_key_a"]]["nclose_id"])
                self.assertEqual((prefix / MODEL_FILE).read_bytes(), saved)
                with self.assertRaises(AssertionError):
                    save_structure_model(prefix, contexts[""].model)
            self.assertEqual(ids, ["SKYPE.nclose.4"] * 3)


if __name__ == "__main__":
    unittest.main()
