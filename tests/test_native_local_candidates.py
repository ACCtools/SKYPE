"""Candidate completeness, multiplicity, saved-state checks and matrix limits."""
from collections import Counter
import copy
import hashlib
import json
from pathlib import Path
import pickle
import tempfile
import unittest
from unittest import mock

import h5py
import numpy as np

import native_local_candidates as candidates_module
from native_local_candidates import extract_candidates, read_completed_context
from structure_nclose import StructureWeights, register_nclose, set_structure_membership, source_signature


def node(name, query_start, chrom, start, strand="+"):
    return [name, 10000, query_start, query_start + 100, strand, chrom, 1000000,
            start, start + 100, 60, 100, 100, 100, 100, 100, "tp:A:P", "cm:i:1",
            "s1:i:100", "s2:i:0", "dv:f:0", "cs:Z::100", 0]


def save_pickle(path, value):
    with path.open("wb") as handle:
        pickle.dump(value, handle)


class Fixture:
    def __init__(self, prefix):
        self.prefix = prefix
        self.nodes = [node("chain", 0, "chr1", 0), node("chain", 100, "chr2", 200),
                      node("chain", 200, "chr1", 0), node("chain", 300, "chr2", 200),
                      node("alias", 0, "chr1", 0), node("alias", 100, "chr2", 200),
                      node("orphan", 0, "chr3", 50), node("orphan", 100, "chr4", 500),
                      node("virtual", 0, "chrV", 800, "-"), node("virtual", 100, "chrV", 100, "-"),
                      node("vcf_synthetic", 0, "chrS", 50), node("vcf_synthetic", 100, "chrS", 300)]
        self.keys = [(0, 3), (4, 5), (6, 7), (8, 9), (10, 11), ("front_jump", 1),
                     ("front_jump", 2), ("front_jump", 3)]
        self.model = dict(version=3, ncloses={}, structures=[], next_nclose_id=1, column_count=4)
        geometries = [(('chr1',100,'L'),('chr2',200,'R')),
                      (('chr1',100,'L'),('chr2',200,'R')),
                      (('chr3',-5,'L'),('chr4',1000001,'R')),
                      (('chrV',0,'R'),('chrV',999999,'L')),
                      (('chrS',150,'L'),('chrS',300,'R'))]
        for key, endpoints in zip(self.keys, geometries):
            register_nclose(self.model, dict(event_key=key, kind="bnd", endpoints=endpoints,
                                            contig_names=("cen" if key == (0,3) else "other",)))
        for key, kind, chrom, start, end in [(self.keys[5], "deletion", "chr5", 0, 100),
                                            (self.keys[6], "insertion", "chr6", 10, 20),
                                            (self.keys[7], "duplication", "chr6", 10, 20)]:
            register_nclose(self.model, dict(event_key=key, kind="indel", indel_kind=kind,
                                            chrom=chrom, st=start, nd=end))
        self.locations = [f"component_{i}.paf" for i in range(4)]
        incidences = [{self.keys[0]:2, self.keys[1]:1}, {self.keys[5]:1},
                      {self.keys[3]:1}, {self.keys[6]:1, self.keys[7]:1}]
        for index, incidence in enumerate(incidences):
            item = dict(structure_id=f"S{index}", feature_index=index, kind="PATH" if index == 0 else "TYPE4",
                        source=self.locations[index], telomere_counts=Counter(), legacy_ids=[])
            set_structure_membership(self.model, item, incidence)
            self.model["structures"].append(item)
        item = dict(structure_id="V", feature_index=None, kind="VIRTUAL_INV", source="virtual:1",
                    raw_weight=.4, telomere_counts=Counter(), legacy_ids=[])
        set_structure_membership(self.model, item, {self.keys[3]:2})
        self.model["structures"].append(item)
        self.weights = np.asarray([0., .1, 1., 0.])
        self.context = StructureWeights(self.model, self.weights, 2.)
        save_pickle(prefix / "01_nclose_data.pkl", {"contig_data": self.nodes})
        save_pickle(prefix / "tot_loc_list.pkl", self.locations)
        (prefix / "censat_endpoint_candidates.tsv").write_text("unitig\tstatus\ncen\trealigned_consistent\n")
        self.model["source_signature"] = source_signature(prefix)
        self.save_model()
        np.save(prefix / "weight.npy", self.weights)
        self.A = np.asarray([[1,-1,0], [0,0,0], [2,-2,1], [0,0,0]], dtype=np.float32)
        self.A_fail = np.asarray([[0,0],[1,0],[0,0],[0,0]], dtype=np.float32)
        self.B = np.asarray([5,10,15,20,25], dtype=np.float32)
        with h5py.File(prefix / "matrix.h5", "w") as handle:
            handle.attrs["matrix_contract"] = "depth_only_v1"
            for name, value in [("A",self.A),("A_fail",self.A_fail),("B",self.B[:3]),("B_fail",self.B[3:])]:
                handle[name] = value
        np.save(prefix / "B.npy", self.B)
        np.save(prefix / "predict_B.npy", np.concatenate([self.A.T @ self.weights, self.A_fail.T @ self.weights]))
        save_pickle(prefix / "23_input.pkl", dict(matrix_contract="depth_only_v1", B_depth_start=0,
                    B_depth_end=3, chr_filt_st_list=[("chr1",i) for i in range(3)]))

    def save_model(self):
        save_pickle(self.prefix / "structure_nclose_model.pkl", self.model)

    def extract(self, matrix_path=None):
        return extract_candidates(self.prefix, self.context, self.nodes, matrix_path)


class CandidateTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.prefix = Path(self.temp.name)
        self.fx = Fixture(self.prefix)

    def test_all_sources_repeated_primitives_aliases_and_no_carriers(self):
        result = self.fx.extract()
        rows = {tuple(tuple(ep) for ep in r["endpoints"]):r for r in result["candidates"]}
        self.assertEqual(len(rows), 7)
        repeated = rows[(('chr1',100,'L'),('chr2',200,'R'))]
        self.assertEqual(repeated["source_primitive_counts"], {"0:3":2,"4:5":1})
        self.assertEqual(repeated["carrier_primitive_counts"], {"S0":5})
        self.assertEqual(len(repeated["source_alignment_occurrences"]), 3)
        self.assertEqual(len(repeated["nclose_ids"]), 1)
        self.assertEqual(repeated["conditional_model_contribution_N"], 0)
        self.assertEqual(repeated["original_candidate_state"], "candidate_NNLS_zero")
        self.assertEqual([r["CEN_source"] for r in repeated["source_alignment_occurrences"]], [True,True,False])
        orphan = rows[(('chr3',-5,'L'),('chr4',1000001,'R'))]
        self.assertEqual(orphan["carrier_count"], 0)
        self.assertEqual(orphan["original_candidate_state"], "candidate_no_modeled_structure")
        symbolic = rows[(('chr6',10,'R'),('chr6',20,'L'))]
        self.assertEqual(symbolic["source_kinds"], ["symbolic_duplication", "symbolic_insertion"])
        self.assertEqual(symbolic["carrier_primitive_counts"], {"S3":2})
        self.assertEqual(symbolic["source_alignment_occurrences"], [])
        self.assertFalse(symbolic["interpret_reference_endpoint_distance_as_variant_size"])
        for row in rows.values():
            expected = sum(result["features"][sid]["weight_N"] * count
                           for sid,count in row["carrier_primitive_counts"].items())
            self.assertAlmostEqual(row["conditional_model_contribution_N"], expected)

    def test_signed_masked_zero_and_virtual_depth_states(self):
        result = self.fx.extract()
        features = result["features"]
        self.assertEqual(features["S0"]["depth_column_state"], "nonzero_fitted_depth")
        self.assertEqual(features["S0"]["negative_fitted_bins"], 1)
        self.assertTrue(features["S0"]["has_signed_depth_entries"])
        self.assertEqual(features["S1"]["depth_column_state"], "entirely_masked_depth_design")
        self.assertEqual(features["S3"]["depth_column_state"], "zero_depth_design")
        self.assertEqual(features["V"]["depth_column_state"], "postfit_qualified_virtual_structure")
        self.assertIsNone(features["V"]["nonzero_fitted_bins"])
        self.assertIsNone(features["V"]["has_signed_depth_entries"])
        virtual = next(r for r in result["candidates"] if "V" in r["carrier_primitive_counts"])
        self.assertAlmostEqual(virtual["conditional_model_contribution_N"], .9)
        self.assertIsNone(virtual["source_alignment_occurrences"][0]["query_gap_bp"])
        low = next(r for r in result["candidates"] if "symbolic_deletion" in r["source_kinds"])
        self.assertEqual(low["original_candidate_state"], "candidate_positive_below_export_threshold")
        synthetic = next(r for r in result["candidates"] if r["endpoints"][0][0] == "chrS")
        self.assertEqual(synthetic["source_alignment_occurrences"][0]["geometry_mode"], "SYNTHETIC_OUTER")
        self.assertIsNone(synthetic["source_alignment_occurrences"][0]["query_gap_bp"])
        json.dumps(result, allow_nan=False)

    def test_missing_matrix_does_not_become_zero_design(self):
        result = self.fx.extract(self.prefix / "absent.h5")
        self.assertEqual(result["provenance"]["matrix"]["state"], "unavailable")
        for sid in ("S0", "S1", "S2", "S3"):
            self.assertEqual(result["features"][sid]["depth_column_state"], "unavailable")
            self.assertIsNone(result["features"][sid]["nonzero_fitted_bins"])
        self.assertEqual(len(result["candidates"]), 7)
        (self.prefix / "matrix.h5").unlink()
        self.assertEqual(self.fx.extract()["provenance"]["matrix"]["state"], "unavailable")

    def test_input_files_and_in_memory_context_are_unchanged(self):
        before = {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in self.prefix.iterdir()}
        model_before = pickle.dumps(self.fx.context.model)
        self.fx.extract()
        self.assertEqual(model_before, pickle.dumps(self.fx.context.model))
        self.assertEqual(before, {p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in self.prefix.iterdir()})
        context, nodes = read_completed_context(self.prefix)
        self.assertEqual(model_before, pickle.dumps(context.model))
        self.assertEqual(nodes, self.fx.nodes)

    def test_stale_saved_model_and_context_nodes_weights_fail(self):
        original = copy.deepcopy(self.fx.context)
        self.fx.context.model["structures"][0]["raw_weight"] = 1
        with self.assertRaisesRegex(ValueError, "context differs"):
            self.fx.extract()
        self.fx.context = original
        self.fx.nodes = copy.deepcopy(self.fx.nodes)
        self.fx.nodes[0][7] += 1
        with self.assertRaisesRegex(ValueError, "nodes differ"):
            self.fx.extract()
        self.fx.nodes[0][7] -= 1
        np.save(self.prefix / "weight.npy", self.fx.weights + 1)
        with self.assertRaisesRegex(ValueError, "coefficient differs"):
            read_completed_context(self.prefix)
        np.save(self.prefix / "weight.npy", self.fx.weights)
        (self.prefix / "telomere_connected_list.txt").write_text("changed\n")
        with self.assertRaisesRegex(ValueError, "source_signature is stale"):
            read_completed_context(self.prefix)

    def test_occurrence_registry_and_source_order_fail_without_rebuild(self):
        self.fx.model["structures"][0]["nclose_counts"][(0,3)] += 1
        self.fx.save_model()
        with self.assertRaisesRegex(ValueError, "occurrence counts disagree"):
            read_completed_context(self.prefix)
        self.fx.model["structures"][0]["nclose_counts"][(0,3)] -= 1
        self.fx.model["structures"][0]["source"] = "wrong_source.paf"
        self.fx.save_model()
        with self.assertRaisesRegex(ValueError, "column-source order"):
            read_completed_context(self.prefix)

    def test_matrix_target_prediction_shape_and_nonfinite_fail(self):
        path = self.prefix / "matrix.h5"
        with h5py.File(path, "r+") as handle:
            handle["B"][0] = 100
        with self.assertRaisesRegex(ValueError, "differs from completed B"):
            self.fx.extract()
        with h5py.File(path, "r+") as handle:
            handle["B"][0] = self.fx.B[0]
            handle["A"][2,0] = 3
        with self.assertRaisesRegex(ValueError, "do not reconstruct"):
            self.fx.extract()
        with h5py.File(path, "r+") as handle:
            handle["A"][2,0] = self.fx.A[2,0]
            handle["A"][0,0] = np.nan
        with self.assertRaisesRegex(ValueError, "nonfinite"):
            self.fx.extract()
        with h5py.File(path, "r+") as handle:
            del handle["A"]
            handle["A"] = self.fx.A[:3]
        with self.assertRaisesRegex(ValueError, "shape differs"):
            self.fx.extract()

    def test_matrix_association_cannot_certify_zero_weight_row_identity(self):
        # A row swap among zero coefficients is invisible in B and prediction.
        # The report must expose this historical binding limit even if it passes.
        with h5py.File(self.prefix / "matrix.h5", "r+") as handle:
            data = handle["A"][:]
            data[[0,3]] = data[[3,0]]
            handle["A"][:] = data
        result = self.fx.extract()
        info = result["provenance"]["matrix"]
        self.assertEqual(info["state"], "current_numerical_association_verified")
        self.assertEqual(info["historical_generation_binding"], "unavailable")
        self.assertIn("no embedded", info["row_source_association_assumption"])

    def test_changed_input_during_extraction_and_changed_implementation_fail(self):
        original = candidates_module._matrix_counts

        def mutate(*args, **kwargs):
            result = original(*args, **kwargs)
            np.save(self.prefix / "weight.npy", self.fx.weights + 1)
            return result

        with mock.patch.object(candidates_module, "_matrix_counts", side_effect=mutate):
            with self.assertRaisesRegex(ValueError, "Input changed during candidate extraction"):
                self.fx.extract()
        with mock.patch.object(candidates_module, "_IMPLEMENTATION_SHA256", "outdated"):
            with self.assertRaisesRegex(ValueError, "implementation changed since import"):
                self.fx.extract()

    def test_finite_saved_values_cannot_overflow_reported_arithmetic(self):
        self.fx.weights[2] = 1e300
        self.fx.model["structures"][2].update(raw_weight=1e300, weight_N=5e299)
        self.fx.save_model()
        np.save(self.prefix / "weight.npy", self.fx.weights)
        with h5py.File(self.prefix / "matrix.h5", "r+") as handle:
            handle["A"][2,0] = 3e38
        with np.errstate(over="ignore", invalid="ignore"):
            with self.assertRaisesRegex(ValueError, "prediction arithmetic is nonfinite"):
                self.fx.extract()
        self.fx.model["n_unit"] = 1e-309
        for structure in self.fx.model["structures"]:
            structure["weight_N"] = structure["raw_weight"] / 1e-309
        self.fx.save_model()
        with self.assertRaisesRegex(ValueError, "weights/n_unit disagree"):
            read_completed_context(self.prefix)


if __name__ == "__main__":
    unittest.main()
