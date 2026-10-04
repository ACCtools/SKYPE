import argparse
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import pysam

import local_hifi_evidence as evidence


class LocalHiFiEvidenceTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.path = Path(self.temp.name)
        self.header = pysam.AlignmentHeader.from_dict({"HD": {"SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 20000}, {"SN": "chr2", "LN": 20000}]})

    def tearDown(self):
        self.temp.cleanup()

    def record(self, name, cigar="1500M1000D1500M", start=1000, flag=0, mapq=60, chrom=0):
        row = pysam.AlignedSegment(self.header)
        row.query_name, row.flag = name, flag
        row.reference_id, row.reference_start, row.mapping_quality = chrom, start, mapq
        row.cigarstring = cigar
        row.query_sequence = "A" * row.infer_query_length()
        return row

    def bam(self, rows, filename="reads.bam", index=False):
        path = self.path / filename
        with pysam.AlignmentFile(str(path), "wb", header=self.header) as handle:
            for row in sorted(rows, key=lambda item: (item.reference_id, item.reference_start)):
                handle.write(row)
        if index:
            pysam.index(str(path))
        return path

    def target(self, tid="target", endpoints=None, contribution=0):
        return dict(id=tid, endpoints=endpoints or [["chr1", 2500, "L"], ["chr1", 3500, "R"]],
            conditional_model_contribution_N=contribution, source_nclose_keys=["1:2"],
            nclose_ids=["N0001"], source_kinds=["bnd"], original_candidate_state="NNLS_zero",
            carrier_count=1, carrier_depth_state_counts={"nonzero_fitted_depth": 1})

    def test_actual_bam_support_rc_duplicates_short_barrier_and_exposure(self):
        rows = [self.record("f1"), self.record("f1"), self.record("f2"), self.record("reverse", flag=16),
                self.record("low", mapq=19), self.record("secondary", flag=256),
                self.record("short_barrier", "1500M500D200M500D1500M")]
        bam = self.bam(rows)
        result, names, stats = evidence.assess_cache([self.target()], bam)
        for row in result:
            self.assertEqual(row["support_read_names"], ["f1", "f2", "reverse"])
            self.assertEqual(row["local_evidence_state"], "LOCAL_HIFI_SUPPORTED")
            self.assertEqual(row["conditional_model_contribution_N"], 0)
            self.assertIsNone(row["normal_evidence"])
            self.assertIsNone(row["ONT_evidence"])
        self.assertEqual(stats["support_exposure_subset_violations"], 0)
        for row in names:
            self.assertTrue(all({"f1", "f2", "reverse"} <= set(group)
                                for group in row["endpoint_exposure_read_names"]))

    def test_observed_cigar_overrides_redundant_SA(self):
        first = self.record("read", "1000M738D1000M2000S", start=10000)
        second = self.record("read", "2000H2000M", start=10000, flag=2048, chrom=1)
        second.set_tag("SA", "chr1,10001,+,2000M2000S,60,0;" * 2)
        target = self.target(endpoints=[["chr1", 11000, "L"], ["chr1", 11738, "R"]])
        rows, _, _ = evidence.assess_cache([target], self.bam([second, first]))
        self.assertTrue(all(row["support_read_names"] == ["read"] for row in rows))
        self.assertTrue(all(row["local_evidence_state"] == "LOCAL_HIFI_BELOW_GATE" for row in rows))

    def test_missing_and_impossible_boundary_are_not_zero(self):
        targets = [self.target("missing", [["chr1", 2500, "L"], ["absent", 500, "R"]]),
                   self.target("invalid", [["chr1", 0, "L"], ["chr1", 3500, "R"]])]
        rows, _, _ = evidence.assess_cache(targets, self.bam([]))
        for row in rows:
            self.assertIsNone(row["distinct_support"])
            self.assertIsNone(row["passes_local_molecule_gate"])
            self.assertEqual(row["local_evidence_state"], "UNASSESSABLE_WITH_REASON")
        lengths = {"chr1": 20000}
        self.assertEqual(evidence.endpoint_status(["chr1", 0, "R"], lengths), "assessable")
        self.assertEqual(evidence.endpoint_status(["chr1", 20000, "L"], lengths), "assessable")
        self.assertEqual(evidence.endpoint_status(["chr1", 20000, "R"], lengths), "unencodable_reference_boundary")

    def test_fixed_tolerance_and_nearby_shared_read_relations(self):
        targets = [self.target("a"), self.target("b", [["chr1", 2650, "L"], ["chr1", 3650, "R"]], 1)]
        rows, _, _ = evidence.assess_cache(targets, self.bam([self.record("r" + str(i)) for i in range(3)]))
        match = {(row["id"], row["tolerance_bp"]): row for row in rows}
        self.assertEqual(match[("b", 100)]["distinct_support"], 0)
        self.assertEqual(match[("b", 500)]["distinct_support"], 3)
        relation = next(row for row in evidence.relations(rows) if row["tolerance_bp"] == 500)
        self.assertFalse(relation["oriented_neighbor_100bp"])
        self.assertTrue(relation["oriented_neighbor_500bp"])
        self.assertEqual(relation["shared_supporting_read_names"], 3)
        self.assertEqual(relation["conditional_model_contributions_N"], [0, 1])

    def test_region_query_deduplicates_overlapping_windows(self):
        bam = self.bam([self.record("r"+str(i)) for i in range(3)], index=True)
        cache, stats = evidence.query_bam(bam, evidence.find_bam_index(bam), [self.target()], self.path)
        self.assertEqual(stats, {"query_region_count": 1, "queried_records": 3})
        rows, _, _ = evidence.assess_cache([self.target()], cache)
        self.assertTrue(all(row["distinct_support"] == 3 for row in rows))

    def test_reference_identity_and_M5_mismatch(self):
        fasta = self.path / "reference.fa"
        fasta.write_text(">chr1\n" + "A"*20000 + "\n>chr2 description\n" + "C"*20000 + "\n")
        reference = evidence.reference_identity(fasta)
        self.assertEqual(reference["sequences"]["chr2"]["md5"], hashlib.md5(b"C"*20000).hexdigest())
        bam = self.bam([], index=True)
        with pysam.AlignmentFile(str(bam), "rb") as handle:
            checked = evidence.check_bam_reference(handle, reference)
        self.assertEqual(checked["sequence_MD5_unavailable_contigs"], ["chr1", "chr2"])
        self.header = pysam.AlignmentHeader.from_dict({"HD": {"SO": "coordinate"},
            "SQ": [{"SN": "chr1", "LN": 20000, "M5": "0"*32}]})
        bad = self.bam([], "bad.bam", index=True)
        with pysam.AlignmentFile(str(bad), "rb") as handle:
            with self.assertRaisesRegex(ValueError, "MD5 mismatch"):
                evidence.check_bam_reference(handle, reference)
        before = evidence.signature(fasta)
        fasta.write_text(">different\nAAA\n")
        with self.assertRaisesRegex(ValueError, "changed"):
            evidence.unchanged(before)

    def test_VCF_all_orientations_and_reference_end_roundtrip(self):
        targets = [self.target("J"+a+b, [["chr1", 0 if a == "R" else 20000, a],
                                        ["chr2", 0 if b == "R" else 20000, b]])
                   for a in "LR" for b in "LR"]
        rows = [dict(row, tolerance_bp=tol, distinct_support=3,
                     conditional_endpoint_exposure=[4, 5], passes_local_molecule_gate=True)
                for row in targets for tol in evidence.TOLERANCES]
        report = evidence.export_vcf(self.path, rows, {"chr1": 20000, "chr2": 20000}, "test")
        for tol in evidence.TOLERANCES:
            self.assertEqual(report[str(tol)]["record_count"], 8)
            with pysam.VariantFile(str(self.path / f"ExperimentalEvidence.{tol}bp.vcf")) as handle:
                self.assertFalse(list(handle.header.samples))
                for row in handle:
                    self.assertIsNone(row.qual)
                    self.assertEqual(list(row.filter), ["ExperimentalEvidence"])
                    self.assertNotIn("GT", row.format)
                    self.assertNotIn("NORMAL_SR", row.info)

    def test_options_are_explicit_and_native_only(self):
        parser = argparse.ArgumentParser()
        evidence.add_arguments(parser)
        args = parser.parse_args([])
        self.assertFalse(evidence.validate_arguments(args))
        self.assertIsNone(evidence.run_from_arguments(args, self.path, None, []))
        self.assertFalse((self.path / "local_hifi_evidence").exists())
        with self.assertRaises(ValueError):
            evidence.validate_arguments(parser.parse_args(["--local-hifi-export"]))
        enabled = parser.parse_args(["--local-hifi-bam", "a", "--local-hifi-reference", "b"])
        self.assertTrue(evidence.validate_arguments(enabled))
        with self.assertRaises(ValueError):
            evidence.validate_arguments(enabled, native=False)

    def test_unassessable_shared_support_remains_missing(self):
        targets = [self.target("missing_a", [["absent", 1000, "L"], ["absent", 2000, "R"]]),
                   self.target("missing_b", [["absent", 1020, "L"], ["absent", 2020, "R"]])]
        rows, _, _ = evidence.assess_cache(targets, self.bam([]))
        for row in evidence.relations(rows):
            self.assertIsNone(row["shared_supporting_read_names"])
            self.assertEqual(row["shared_support_assessment"], "unassessable")
            self.assertTrue(row["oriented_neighbor_100bp"])

    def test_model_reference_mismatch_is_local_and_export_explained(self):
        targets = [self.target("valid"), self.target("affected", [["chr1", 2500, "L"], ["chr2", 3500, "R"]])]
        reference = {"sequences": {"chr1": {"length": 20000}, "chr2": {"length": 20000}}}
        nodes = [[None]*5+["chr1", 20000], [None]*5+["chr2", 19999]]
        association = evidence.associate_model_reference(targets, nodes, reference)
        self.assertEqual(association["affected_candidate_ids"], ["affected"])
        self.assertEqual(association["mismatched_contigs"],
                         {"chr2": {"native_node_lengths": [19999], "supplied_reference_length": 20000}})
        rows, witnesses, _ = evidence.assess_cache(targets, self.bam([self.record("r"+str(i)) for i in range(3)]))
        for row in rows:
            if row["id"] == "valid":
                self.assertEqual(row["distinct_support"], 3)
                self.assertTrue(row["passes_local_molecule_gate"])
            else:
                self.assertIsNone(row["distinct_support"])
                self.assertIsNone(row["passes_local_molecule_gate"])
                self.assertEqual(row["conditional_endpoint_exposure"], [3, None])
                self.assertIsNone(row["exposure_union"])
                self.assertIsNone(row["exposure_intersection"])
        output = evidence.export_vcf(self.path, rows, {"chr1": 20000, "chr2": 20000}, "test")
        for row in output.values():
            self.assertEqual(row["candidate_count"], 1)
            self.assertEqual(row["unencodable"], [{"id": "affected", "reason": ["assessable", "model_reference_length_mismatch"]}])

    def test_pipeline_routes_only_explicit_options_to_stage31(self):
        import pipeline
        common = dict(CELL_LINE="test", PREFIX=str(self.path / "output"),
            ctg_paf="ctg.paf", ctg_aln_paf="ctg.aln.paf", utg_paf="utg.paf",
            utg_aln_paf="utg.aln.paf", depth_loc="test.win.stat.gz", thread=1,
            dep_folder=str(self.path), is_progress=False, skype_force=False, graph_depth=2,
            skype_start_at=31, print_args=True, alignasm_ref="ref.fa", chr_fai="ref.fa.fai",
            tel_bed="tel.bed", rpt_bed="repeat.bed", rcs_bed="censat.bed", cyt_bed="cyt.bed",
            raw_rescue_method="off")
        with patch.object(pipeline, "subprocess_print") as command:
            pipeline.run_skype(**common)
        report = next(call.args[0] for call in command.call_args_list if "31_depth_analysis.py" in call.args[0][1])
        self.assertFalse(any(part.startswith("--local-hifi") for part in report))
        options = "--local-hifi-bam 'read input.bam' --local-hifi-reference ref.fa --local-hifi-export --local-hifi-same-input-as-assembly yes"
        with patch.object(pipeline, "subprocess_print") as command:
            pipeline.run_skype(**common, option_skype=options)
        report = next(call.args[0] for call in command.call_args_list if "31_depth_analysis.py" in call.args[0][1])
        self.assertEqual(report[report.index("--local-hifi-bam")+1], "read input.bam")
        self.assertEqual(report[report.index("--local-hifi-same-input-as-assembly")+1], "yes")
        self.assertIn("--local-hifi-export", report)
        with self.assertRaisesRegex(pipeline.SkypeArgumentError, "native assembly"):
            pipeline.run_skype(**common, option_skype=options, benchmark_vcf_loc="input.vcf")
        output = Path(common["PREFIX"])
        output.mkdir()
        for filename in ("total_cov.png", "nclose_report.tsv", "SV_call_result.vcf", "SKYPE_result.bed"):
            (output / filename).write_text("existing")
        common["skype_start_at"] = 0
        with patch.object(pipeline, "subprocess_print") as command:
            pipeline.run_skype(**common, option_skype=options)
        self.assertEqual(command.call_count, 1)
        self.assertEqual(Path(command.call_args.args[0][1]).name, "local_hifi_evidence.py")


if __name__ == "__main__":
    unittest.main()
