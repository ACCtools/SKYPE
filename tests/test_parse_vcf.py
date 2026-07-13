from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


SKYPE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKYPE_ROOT))

from parse_vcf import (  # noqa: E402
    VcfParseError,
    parse_vcf_events,
    select_vcf_type4_graph_events,
)


REFERENCE_LENGTHS = {"chr1": 1_000_000, "chr2": 1_000_000,
                     "chr3": 1_000_000, "chr4": 1_000_000}


def vcf_text(
    source: str | None,
    records: list[str],
    *,
    extra_headers: tuple[str, ...] = (),
) -> str:
    lines = ["##fileformat=VCFv4.2"]
    if source is not None:
        lines.append(f"##source={source}")
    lines.extend(extra_headers)
    lines.extend(
        f"##contig=<ID={chrom},length={length}>"
        for chrom, length in REFERENCE_LENGTHS.items()
    )
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    lines.extend(records)
    return "\n".join(lines) + "\n"


def bnd_record(
    chrom: str,
    pos: int,
    record_id: str,
    alt: str,
    *,
    mate_id: str | None = None,
    mate_key: str = "MATEID",
    strands: str | None = None,
    filt: str = "PASS",
    extra_info: tuple[str, ...] = (),
) -> str:
    info = ["SVTYPE=BND"]
    if mate_id is not None:
        info.append(f"{mate_key}={mate_id}")
    if strands is not None:
        info.append(f"STRANDS={strands}")
    info.extend(extra_info)
    return (
        f"{chrom}\t{pos}\t{record_id}\tN\t{alt}\t60\t{filt}\t"
        + ";".join(info)
    )


class ParseVcfTests(unittest.TestCase):
    def parse(
        self,
        text: str,
        *,
        pass_filters=("PASS", "."),
        ins_alt_alignments=None,
    ):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        path = Path(temporary.name) / "input.vcf"
        path.write_text(text, encoding="utf-8")
        return parse_vcf_events(
            path,
            REFERENCE_LENGTHS,
            pass_filters=pass_filters,
            ins_alt_alignments=ins_alt_alignments,
        )

    def test_unknown_source_warns_and_generic_singleton_is_used(self) -> None:
        text = vcf_text(
            "MindaV0.0.2",
            [bnd_record("chr1", 100, "bnd1", "N[chr2:200[")],
        )
        with self.assertLogs(level="WARNING") as captured:
            result = self.parse(text)

        self.assertIsNone(result.source.caller)
        self.assertEqual(len(result.nclose_specs), 1)
        self.assertEqual(
            (result.nclose_specs[0].dir_a, result.nclose_specs[0].dir_b),
            ("+", "+"),
        )
        self.assertIn("compatibility is not guaranteed", "\n".join(captured.output))

    def test_custom_filter_replaces_defaults(self) -> None:
        text = vcf_text(
            "SAVANAv1.0.3",
            [
                "chr1\t100\tinv1\tN\t<INV>\t60\tCandidate\t"
                "SVTYPE=INV;END=300"
            ],
        )
        default_result = self.parse(text)
        custom_result = self.parse(text, pass_filters=("Candidate",))

        self.assertEqual(len(default_result.nclose_specs), 0)
        self.assertEqual(default_result.summary["filtered_records"], 1)
        self.assertEqual(len(custom_result.nclose_specs), 2)
        self.assertEqual(custom_result.summary["used_inv_events"], 1)

    def test_half_filtered_bnd_pair_is_dropped_together(self) -> None:
        text = vcf_text(
            "SAVANAv1.0.3",
            [
                bnd_record(
                    "chr1", 100, "pair_1", "N[chr2:200[", mate_id="pair_2"
                ),
                bnd_record(
                    "chr2",
                    200,
                    "pair_2",
                    "]chr1:100]N",
                    mate_id="pair_1",
                    filt="LowQual",
                ),
            ],
        )
        result = self.parse(text)

        self.assertEqual(result.nclose_specs, [])
        self.assertEqual(result.summary["filtered_bnd_pairs"], 1)
        self.assertEqual(result.summary["filtered_records"], 2)
        reasons = {row[3] for row in result.skipped_records}
        self.assertEqual(reasons, {"mate_filtered", "FILTER=LowQual"})

    def test_sniffles_singleton_bnd_needs_no_mate(self) -> None:
        text = vcf_text(
            "Sniffles2_2.2",
            [bnd_record("chr1", 100, "Sniffles2.BND.1", "N[chr2:200[")],
        )
        result = self.parse(text)

        self.assertEqual(result.source.caller, "sniffles2")
        self.assertEqual(result.source.version, "2.2")
        self.assertEqual(result.summary["singleton_bnd_events"], 1)
        self.assertEqual(result.summary["missing_mates"], 1)
        self.assertEqual(len(result.nclose_specs), 1)

    def test_duplicate_mate_id_uses_unique_coordinate_fallback(self) -> None:
        text = vcf_text(
            "svaba(v1.2.0)",
            [
                bnd_record("chr1", 100, "a", "N[chr2:200[", mate_id="dup"),
                bnd_record("chr2", 200, "dup", "]chr1:100]N", mate_id="a"),
                bnd_record(
                    "chr3", 300, "dup", "N[chr4:400[", filt="LowQual"
                ),
            ],
        )
        result = self.parse(text)

        self.assertEqual(result.source.caller, "svaba")
        self.assertEqual(result.source.version, "1.2.0")
        self.assertEqual(result.summary["coordinate_fallback_pairs"], 1)
        self.assertEqual(len(result.nclose_specs), 1)
        self.assertEqual(result.nclose_specs[0].record_ids, ("a", "dup"))

    def test_gripss_header_and_single_breakend_without_remote_locus(self) -> None:
        text = vcf_text(
            None,
            [
                bnd_record("chr1", 100, "gridss_1", "N]chr2:200]",
                           mate_id="gridss_2"),
                bnd_record("chr2", 200, "gridss_2", "N]chr1:100]",
                           mate_id="gridss_1"),
                bnd_record(
                    "chr3",
                    300,
                    "gridss_sgl",
                    ".ACGT",
                    extra_info=("EVENTTYPE=SGL",),
                ),
            ],
            extra_headers=(
                "##gridssVersion=2.13.2-gridss",
                "##gripssVersion=2.2",
            ),
        )
        with self.assertLogs(level="WARNING"):
            result = self.parse(text)

        self.assertEqual(result.source.caller, "gripss")
        self.assertEqual(result.source.version, "2.2")
        self.assertEqual(len(result.nclose_specs), 1)
        self.assertEqual(
            result.summary["single_breakends_without_remote_locus"], 1
        )
        self.assertEqual(result.summary["malformed_records"], 0)
        self.assertIn(
            "single_breakend_without_remote_locus",
            {row[3] for row in result.skipped_records},
        )

    def test_manta_mateid_pairs_microhomology_shifted_coordinates(self) -> None:
        text = vcf_text(
            "GenerateSVCandidates UNKNOWN",
            [
                bnd_record(
                    "chr1",
                    100,
                    "MantaBND:1:0",
                    "N]chr2:202]",
                    mate_id="MantaBND:1:1",
                    extra_info=("HOMLEN=2",),
                ),
                bnd_record(
                    "chr2",
                    200,
                    "MantaBND:1:1",
                    "N]chr1:102]",
                    mate_id="MantaBND:1:0",
                    extra_info=("HOMLEN=2",),
                ),
                "chr1\t500\tMantaDEL:1\tAAAA\tA\t60\tPASS\t"
                "SVTYPE=DEL;SVLEN=-3",
            ],
        )
        result = self.parse(text)

        self.assertEqual(result.source.caller, "manta")
        self.assertIsNone(result.source.version)
        self.assertEqual(len(result.nclose_specs), 1)
        self.assertEqual(result.nclose_specs[0].record_ids,
                         ("MantaBND:1:0", "MantaBND:1:1"))
        self.assertEqual(result.summary["missing_mates"], 0)
        self.assertEqual(result.summary["mate_id_coordinate_disagreements"], 1)
        self.assertEqual(result.summary["orientation_mismatches"], 0)
        self.assertEqual(result.summary["skipped_small_indel_size"], 1)
        self.assertEqual(result.summary["malformed_records"], 0)

    def test_manta_inv3_and_inv5_each_emit_one_junction(self) -> None:
        text = vcf_text(
            "GenerateSVCandidates UNKNOWN",
            [
                "chr1\t100\tMantaINV:1:INV3\tN\t<INV>\t60\tPASS\t"
                "END=300;SVTYPE=INV;SVLEN=200;INV3",
                "chr1\t500\tMantaINV:1:INV5\tN\t<INV>\t60\tPASS\t"
                "END=800;SVTYPE=INV;SVLEN=300;INV5",
            ],
        )
        result = self.parse(text)

        self.assertEqual(result.source.caller, "manta")
        self.assertEqual(len(result.nclose_specs), 2)
        observed = {
            spec.record_ids[0]: (
                spec.event_name.rsplit("_", 1)[-1],
                spec.dir_a,
                spec.dir_b,
            )
            for spec in result.nclose_specs
        }
        self.assertEqual(
            observed,
            {
                "MantaINV:1:INV3": ("left", "+", "-"),
                "MantaINV:1:INV5": ("right", "-", "+"),
            },
        )
        self.assertEqual(result.summary["used_inv_events"], 2)
        self.assertEqual(result.summary["used_inv_junctions"], 2)
        self.assertEqual(result.summary["used_manta_inv3_records"], 1)
        self.assertEqual(result.summary["used_manta_inv5_records"], 1)

    def test_ambiguous_manta_inv_flags_are_skipped(self) -> None:
        text = vcf_text(
            "GenerateSVCandidates UNKNOWN",
            [
                "chr1\t100\tmissing_flags\tN\t<INV>\t60\tPASS\t"
                "END=300;SVTYPE=INV;SVLEN=200",
                "chr1\t500\tconflicting_flags\tN\t<INV>\t60\tPASS\t"
                "END=800;SVTYPE=INV;SVLEN=300;INV3;INV5",
            ],
        )
        result = self.parse(text)

        self.assertEqual(result.nclose_specs, [])
        self.assertEqual(result.summary["used_inv_events"], 0)
        self.assertEqual(result.summary["ambiguous_manta_inv_records"], 2)
        self.assertEqual(result.summary["malformed_records"], 2)
        self.assertEqual(
            {row[3] for row in result.skipped_records},
            {"missing_Manta_INV3_INV5", "conflicting_Manta_INV3_INV5"},
        )

    def test_absent_reference_contig_drops_only_affected_pair(self) -> None:
        text = vcf_text(
            "SAVANAv1.0.3",
            [
                bnd_record("chr1", 100, "keep_1", "N[chr2:200[",
                           mate_id="keep_2"),
                bnd_record("chr2", 200, "keep_2", "]chr1:100]N",
                           mate_id="keep_1"),
                bnd_record("chr1", 300, "drop_1", "N[chrAlt:400[",
                           mate_id="drop_2"),
                bnd_record("chrAlt", 400, "drop_2", "]chr1:300]N",
                           mate_id="drop_1"),
            ],
            extra_headers=("##contig=<ID=chrAlt,length=500000>",),
        )
        with self.assertLogs(level="WARNING") as captured:
            result = self.parse(text)

        self.assertEqual(len(result.nclose_specs), 1)
        self.assertEqual(result.nclose_specs[0].record_ids, ("keep_1", "keep_2"))
        self.assertEqual(result.summary["reference_header_contigs_absent"], 1)
        self.assertEqual(result.summary["reference_used_contigs_absent"], 1)
        self.assertEqual(result.summary["unknown_reference_contig_events"], 1)
        self.assertEqual(result.summary["unknown_reference_contig_records"], 2)
        self.assertEqual(result.summary["orientation_mismatches"], 0)
        self.assertIn("affected events will be excluded", "\n".join(captured.output))

    def test_shared_reference_contig_length_mismatch_remains_fatal(self) -> None:
        text = vcf_text(
            "SAVANAv1.0.3",
            [bnd_record("chr1", 100, "bnd", "N[chr2:200[")],
        ).replace(
            "##contig=<ID=chr1,length=1000000>",
            "##contig=<ID=chr1,length=999999>",
        )

        with self.assertRaisesRegex(VcfParseError, "VCF length 999999"):
            self.parse(text)

    def test_existing_mate_with_orientation_mismatch_is_excluded(self) -> None:
        text = vcf_text(
            "SAVANAv1.0.3",
            [
                bnd_record("chr1", 100, "a", "N[chr2:200[", mate_id="b"),
                bnd_record("chr2", 200, "b", "[chr1:100[N", mate_id="a"),
            ],
        )
        result = self.parse(text)

        self.assertEqual(result.nclose_specs, [])
        self.assertEqual(result.summary["orientation_mismatches"], 1)
        self.assertEqual(len(result.orientation_mismatches), 1)
        self.assertEqual(
            {row[3] for row in result.skipped_records},
            {"BND_ALT_orientation_mismatch"},
        )

    def test_raw_severus_pre17_normalizes_all_strand_forms(self) -> None:
        records = [
            bnd_record(
                "chr1", 100, "pp_1", "N]chr2:200]", mate_id="pp_2",
                mate_key="MATE_ID", strands="++"
            ),
            bnd_record(
                "chr2", 200, "pp_2", "N]chr1:100]", mate_id="pp_1",
                mate_key="MATE_ID", strands="++"
            ),
            bnd_record(
                "chr1", 300, "pm_1", "]chr2:400]N", mate_id="pm_2",
                mate_key="MATE_ID", strands="+-"
            ),
            bnd_record(
                "chr2", 400, "pm_2", "N[chr1:300[", mate_id="pm_1",
                mate_key="MATE_ID", strands="+-"
            ),
            bnd_record(
                "chr1", 500, "mp_1", "N[chr2:600[", mate_id="mp_2",
                mate_key="MATE_ID", strands="-+"
            ),
            bnd_record(
                "chr2", 600, "mp_2", "]chr1:500]N", mate_id="mp_1",
                mate_key="MATE_ID", strands="-+"
            ),
            bnd_record(
                "chr1", 700, "mm_1", "[chr2:800[N", mate_id="mm_2",
                mate_key="MATE_ID", strands="--"
            ),
            bnd_record(
                "chr2", 800, "mm_2", "[chr1:700[N", mate_id="mm_1",
                mate_key="MATE_ID", strands="--"
            ),
        ]
        result = self.parse(vcf_text("Severus_v1.6", records))

        observed = {
            spec.record_ids[0]: (spec.dir_a, spec.dir_b)
            for spec in result.nclose_specs
        }
        self.assertEqual(
            observed,
            {
                "pp_1": ("+", "-"),
                "pm_1": ("+", "+"),
                "mp_1": ("-", "-"),
                "mm_1": ("-", "+"),
            },
        )
        self.assertEqual(result.summary["severus_corrected_records"], 4)
        self.assertEqual(result.summary["orientation_mismatches"], 0)

    def test_old_severus_missing_strands_drops_pair_and_singleton(self) -> None:
        text = vcf_text(
            "Severus_v1.2",
            [
                bnd_record(
                    "chr1", 100, "pair_1", "]chr2:200]N", mate_id="pair_2",
                    mate_key="MATE_ID"
                ),
                bnd_record(
                    "chr2", 200, "pair_2", "N[chr1:100[", mate_id="pair_1",
                    mate_key="MATE_ID", strands="+-"
                ),
                bnd_record("chr3", 300, "single", "N[chr4:400["),
            ],
        )
        result = self.parse(text)

        self.assertEqual(result.nclose_specs, [])
        self.assertEqual(result.summary["missing_or_invalid_severus_strands"], 2)
        self.assertEqual(result.summary["severus_strands_skipped_events"], 2)
        self.assertEqual(
            {row[3] for row in result.skipped_records},
            {"missing_or_invalid_Severus_STRANDS"},
        )

    def test_raw_old_severus_singleton_is_corrected_with_strands(self) -> None:
        text = vcf_text(
            "Severus_v1.2",
            [
                bnd_record(
                    "chr1", 100, "single", "]chr2:200]N", strands="+-"
                )
            ],
        )
        result = self.parse(text)

        self.assertEqual(len(result.nclose_specs), 1)
        self.assertEqual(
            (result.nclose_specs[0].dir_a, result.nclose_specs[0].dir_b),
            ("+", "+"),
        )
        self.assertEqual(result.summary["severus_corrected_records"], 1)

    def test_severus_17_uses_alt_without_strands(self) -> None:
        text = vcf_text(
            "Severus_v1.7",
            [
                bnd_record(
                    "chr1", 100, "pair_1", "N[chr2:200[", mate_id="pair_2",
                    mate_key="MATE_ID"
                ),
                bnd_record(
                    "chr2", 200, "pair_2", "]chr1:100]N", mate_id="pair_1",
                    mate_key="MATE_ID"
                ),
            ],
        )
        result = self.parse(text)

        self.assertEqual(len(result.nclose_specs), 1)
        self.assertEqual(result.summary["severus_corrected_records"], 0)

    def test_severus_17_junction_del_pairs_by_alt_without_mateid(self) -> None:
        text = vcf_text(
            "Severus_v1.7",
            [
                "chr1\t100\tseverus_BND1_1\tN\tN[chr1:100100[\t60\tPASS\t"
                "SVTYPE=DEL;SVLEN=100000;END=100100;STRANDS=+-",
                "chr1\t100100\tseverus_BND1_2\tN\t]chr1:100]N\t60\tPASS\t"
                "SVTYPE=DEL;SVLEN=100000;END=100;STRANDS=-+",
            ],
        )
        result = self.parse(text)

        self.assertEqual(result.nclose_specs, [])
        self.assertEqual(len(result.type4_events), 1)
        self.assertEqual(result.type4_events[0]["event_type"], "front_jump")
        self.assertEqual(result.type4_events[0]["svlen"], 100_000)
        self.assertEqual(result.summary["used_bnd_type4_events"], 1)
        self.assertEqual(result.summary["coordinate_fallback_pairs"], 1)
        self.assertEqual(result.summary["severus_bracket_alt_overrides"], 2)
        self.assertEqual(result.summary["missing_mates"], 0)

    def test_severus_17_junction_inv_pair_is_one_nclose(self) -> None:
        text = vcf_text(
            "Severus_v1.7",
            [
                "chr1\t200000\tseverus_BND2_1\tN\tN]chr1:300000]\t60\tPASS\t"
                "SVTYPE=INV;SVLEN=100000;END=300000;STRANDS=++",
                "chr1\t300000\tseverus_BND2_2\tN\tN]chr1:200000]\t60\tPASS\t"
                "SVTYPE=INV;SVLEN=100000;END=200000;STRANDS=++",
            ],
        )
        result = self.parse(text)

        self.assertEqual(len(result.nclose_specs), 1)
        spec = result.nclose_specs[0]
        self.assertEqual(
            spec.record_ids,
            ("severus_BND2_1", "severus_BND2_2"),
        )
        self.assertEqual((spec.dir_a, spec.dir_b), ("+", "-"))
        self.assertEqual(result.type4_events, [])
        self.assertEqual(result.summary["used_inv_events"], 0)
        self.assertEqual(result.summary["used_bnd_events"], 1)
        self.assertEqual(result.summary["coordinate_fallback_pairs"], 1)
        self.assertEqual(result.summary["severus_bracket_alt_overrides"], 2)

    def test_bracket_alt_override_is_severus_only(self) -> None:
        text = vcf_text(
            "Sniffles2_2.2",
            [
                "chr1\t100\tinv\tN\tN]chr1:300]\t60\tPASS\t"
                "SVTYPE=INV;END=300"
            ],
        )
        result = self.parse(text)

        self.assertEqual(result.summary["severus_bracket_alt_overrides"], 0)
        self.assertEqual(result.summary["used_inv_events"], 1)
        self.assertEqual(len(result.nclose_specs), 2)

    def test_unknown_severus_version_drops_only_bnd(self) -> None:
        text = vcf_text(
            "Severus_unknown",
            [
                bnd_record(
                    "chr1", 100, "bnd", "N[chr2:200[", strands="+-"
                ),
                "chr1\t300\tinv\tN\t<INV>\t60\tPASS\tSVTYPE=INV;END=500",
            ],
        )
        with self.assertLogs(level="WARNING") as captured:
            result = self.parse(text)

        self.assertEqual(result.summary["unknown_severus_version_bnd_records"], 1)
        self.assertEqual(result.summary["used_inv_events"], 1)
        self.assertEqual(len(result.nclose_specs), 2)
        self.assertIn("version could not be determined", "\n".join(captured.output))

    def test_nanomonsv_anchor_and_svinslen_are_parsed_without_rewrite(self) -> None:
        text = vcf_text(
            "nanomonsv-0.7.0",
            [
                bnd_record("chr1", 100, "r_0", "[chr2:200[TCGG", mate_id="r_1"),
                bnd_record("chr2", 200, "r_1", "[chr1:100[G", mate_id="r_0"),
                "chr1\t300\tins1\tA\t<INS>\t60\tPASS\t"
                "SVTYPE=INS;SVINSLEN=150000",
            ],
        )
        # fileformat + source + 4 contigs + #CHROM => first record is line 8.
        alignments = {
            "vcf_ins_10_ins1": {
                "chrom": "chr1",
                "st": 400_000,
                "nd": 550_000,
                "paf_rows": [],
            }
        }
        result = self.parse(text, ins_alt_alignments=alignments)

        self.assertEqual(len(result.nclose_specs), 1)
        self.assertEqual(len(result.type4_events), 1)
        self.assertEqual(result.type4_events[0]["svlen"], 150_000)

    def test_symbolic_del_dup_threshold_uses_end_span(self) -> None:
        text = vcf_text(
            "Sniffles2_2.2",
            [
                "chr1\t100\tlarge_del\tN\t<DEL>\t60\tPASS\t"
                "SVTYPE=DEL;END=100100;SVLEN=-10",
                "chr1\t200000\tsmall_dup\tN\t<DUP>\t60\tPASS\t"
                "SVTYPE=DUP;END=299999;SVLEN=200000",
            ],
        )
        result = self.parse(text)

        self.assertEqual(result.nclose_specs, [])
        self.assertEqual(len(result.type4_events), 1)
        event = result.type4_events[0]
        self.assertEqual(event["vcf_id"], "large_del")
        self.assertEqual(event["svlen"], 100_000)
        self.assertEqual(event["reported_svlen"], -10)
        self.assertIn(
            ("small_dup", "small_indel_size"),
            {(row[1], row[3]) for row in result.skipped_records},
        )

    def test_bnd_del_dup_follow_the_same_100kb_threshold(self) -> None:
        records = [
            bnd_record(
                "chr1", 100, "small_del_1", "N[chr1:100099[",
                mate_id="small_del_2",
            ),
            bnd_record(
                "chr1", 100099, "small_del_2", "]chr1:100]N",
                mate_id="small_del_1",
            ),
            bnd_record(
                "chr1", 200000, "large_del_1", "N[chr1:300000[",
                mate_id="large_del_2",
            ),
            bnd_record(
                "chr1", 300000, "large_del_2", "]chr1:200000]N",
                mate_id="large_del_1",
            ),
            bnd_record(
                "chr1", 400000, "small_dup_1", "]chr1:499999]N",
                mate_id="small_dup_2",
            ),
            bnd_record(
                "chr1", 499999, "small_dup_2", "N[chr1:400000[",
                mate_id="small_dup_1",
            ),
            bnd_record(
                "chr1", 500000, "large_dup_1", "]chr1:600000]N",
                mate_id="large_dup_2",
            ),
            bnd_record(
                "chr1", 600000, "large_dup_2", "N[chr1:500000[",
                mate_id="large_dup_1",
            ),
            bnd_record(
                "chr1", 700000, "inv_1", "N]chr1:800000]",
                mate_id="inv_2",
            ),
            bnd_record(
                "chr1", 800000, "inv_2", "N]chr1:700000]",
                mate_id="inv_1",
            ),
            bnd_record(
                "chr1", 900000, "bnd_1", "N[chr2:100000[",
                mate_id="bnd_2",
            ),
            bnd_record(
                "chr2", 100000, "bnd_2", "]chr1:900000]N",
                mate_id="bnd_1",
            ),
        ]
        result = self.parse(vcf_text("SAVANAv1.0.3", records))

        self.assertEqual(
            [(event["vcf_id"], event["event_type"]) for event in result.type4_events],
            [("large_del_1", "front_jump"), ("large_dup_1", "back_jump")],
        )
        self.assertEqual(
            [spec.record_ids for spec in result.nclose_specs],
            [("inv_1", "inv_2"), ("bnd_1", "bnd_2")],
        )
        self.assertEqual(result.summary["used_bnd_type4_events"], 2)
        self.assertEqual(result.summary["skipped_bnd_small_indel_events"], 2)
        self.assertEqual(
            {
                row[1]
                for row in result.skipped_records
                if row[3] == "small_indel_size"
            },
            {"small_del_1", "small_del_2", "small_dup_1", "small_dup_2"},
        )

    def test_ins_threshold_uses_only_svlen_not_aligned_span(self) -> None:
        text = vcf_text(
            "Sniffles2_2.2",
            [
                "chr1\t100\tins_small\tN\t<INS>\t60\tPASS\t"
                "SVTYPE=INS;SVLEN=99999",
                "chr1\t200\tins_large\tN\t<INS>\t60\tPASS\t"
                "SVTYPE=INS;SVLEN=100000",
            ],
        )
        # fileformat + source + 4 contigs + #CHROM => records are lines 8/9.
        alignments = {
            "vcf_ins_8_ins_small": {
                "chrom": "chr1",
                "st": 100_000,
                "nd": 400_000,
                "paf_rows": [["small_with_large_reference_span"]],
            },
            "vcf_ins_9_ins_large": {
                "chrom": "chr1",
                "st": 500_000,
                "nd": 501_000,
                "paf_rows": [["large_with_small_reference_span"]],
            },
        }

        result = self.parse(text, ins_alt_alignments=alignments)

        self.assertEqual(
            [event["vcf_id"] for event in result.type4_events],
            ["ins_large"],
        )
        self.assertEqual(result.type4_events[0]["svlen"], 100_000)
        self.assertEqual(
            (result.type4_events[0]["st"], result.type4_events[0]["nd"]),
            (500_000, 501_000),
        )
        self.assertIn(
            ("ins_small", "small_indel_size"),
            {(row[1], row[3]) for row in result.skipped_records},
        )

    def test_all_vcf_event_classes_are_partitioned_and_ins_is_not_graph_eligible(
        self,
    ) -> None:
        records = [
            "chr1\t10000\tsymbolic_inv\tN\t<INV>\t60\tPASS\t"
            "SVTYPE=INV;END=210000",
            "chr1\t220000\tsymbolic_del\tN\t<DEL>\t60\tPASS\t"
            "SVTYPE=DEL;END=320000;SVLEN=-100000",
            "chr1\t330000\tsymbolic_dup\tN\t<DUP>\t60\tPASS\t"
            "SVTYPE=DUP;END=430000;SVLEN=100000",
            bnd_record(
                "chr1", 440000, "bnd_del_1", "N[chr1:540000[",
                mate_id="bnd_del_2",
            ),
            bnd_record(
                "chr1", 540000, "bnd_del_2", "]chr1:440000]N",
                mate_id="bnd_del_1",
            ),
            bnd_record(
                "chr1", 550000, "bnd_dup_1", "]chr1:650000]N",
                mate_id="bnd_dup_2",
            ),
            bnd_record(
                "chr1", 650000, "bnd_dup_2", "N[chr1:550000[",
                mate_id="bnd_dup_1",
            ),
            bnd_record(
                "chr1", 660000, "bnd_inv_1", "N]chr1:760000]",
                mate_id="bnd_inv_2",
            ),
            bnd_record(
                "chr1", 760000, "bnd_inv_2", "N]chr1:660000]",
                mate_id="bnd_inv_1",
            ),
            bnd_record(
                "chr1", 770000, "translocation_1", "N[chr2:100000[",
                mate_id="translocation_2",
            ),
            bnd_record(
                "chr2", 100000, "translocation_2", "]chr1:770000]N",
                mate_id="translocation_1",
            ),
            "chr1\t800000\tins\tN\t<INS>\t60\tPASS\t"
            "SVTYPE=INS;SVLEN=100000",
        ]
        ins_line_no = 7 + len(records)
        alignments = {
            f"vcf_ins_{ins_line_no}_ins": {
                "chrom": "chr1",
                "st": 800000,
                "nd": 801000,
                "paf_rows": [["ins_alignment"]],
            }
        }

        result = self.parse(
            vcf_text("SAVANAv1.0.3", records),
            ins_alt_alignments=alignments,
        )

        self.assertEqual(
            [spec.record_ids for spec in result.nclose_specs],
            [
                ("symbolic_inv",),
                ("symbolic_inv",),
                ("bnd_inv_1", "bnd_inv_2"),
                ("translocation_1", "translocation_2"),
            ],
        )
        self.assertEqual(
            [event["svtype"] for event in result.type4_events],
            ["DEL", "DUP", "BND", "BND", "INS"],
        )
        self.assertEqual(
            [
                event["svtype"]
                for event in select_vcf_type4_graph_events(result.type4_events)
            ],
            ["DEL", "DUP", "BND", "BND"],
        )


if __name__ == "__main__":
    unittest.main()
