import gzip
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import h5py
import numpy as np

SKYPE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKYPE_DIR))

import full_assembly_pipeline as full  # noqa: E402
from skype_utils import (  # noqa: E402
    load_nclose_nodes,
    load_pipeline_input,
    pipeline_input_is_full_assembly,
    save_pipeline_input,
)


def write_depth(path: Path, values):
    with gzip.open(path, "wt") as handle:
        for chrom, st, nd, depth in values:
            length = nd - st
            handle.write(
                f"{chrom}\t{st}\t{nd}\t{length}\t{length}\t"
                f"{depth * length}\t1\t{depth}\n"
            )


class FullAssemblyModeTests(unittest.TestCase):
    def test_pipeline_input_round_trip(self):
        with tempfile.TemporaryDirectory() as temporary:
            save_pipeline_input(
                temporary,
                full_assembly_input=True,
                full_assembly_path="/data/truth.fa",
                full_assembly_paf_path="/data/truth.hs1.aln.paf",
                reference="hs1",
            )
            config = load_pipeline_input(temporary)
            self.assertTrue(pipeline_input_is_full_assembly(config))
            self.assertEqual(config["full_assembly_path"], "/data/truth.fa")
            self.assertEqual(config["reference"], "hs1")

    def test_contig_paths_matrix_solver_and_plots(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prefix = root / "result"
            prefix.mkdir()
            assembly = root / "truth.fasta"
            assembly.write_text(
                ">ctgA\n" + "A" * 1000 + "\n"
                ">ctgB\n" + "C" * 1000 + "\n"
                ">ctgC\n" + "G" * 1000 + "\n",
                encoding="ascii",
            )
            paf = root / "truth.hs1.aln.paf"
            # PAF order is deliberately different from FASTA order. Secondary
            # alignments must not create paths or contribute to depth.
            paf.write_text(
                "ctgB\t1000\t0\t500\t+\tchr1\t200000\t100000\t100500\t500\t500\t60\ttp:A:P\tcs:Z::500\n"
                "ctgB\t1000\t500\t1000\t-\tchr2\t200000\t0\t500\t500\t500\t60\ttp:A:P\tcs:Z::500\n"
                "ctgB\t1000\t0\t1000\t+\tchr2\t200000\t100000\t101000\t1000\t1000\t0\ttp:A:S\n"
                "ctgA\t1000\t0\t1000\t+\tchr1\t200000\t0\t1000\t1000\t1000\t60\ttp:A:P\tcs:Z::1000\n",
                encoding="ascii",
            )

            per_contig_depth = {
                "0": [
                    ("chr1", 0, 100000, 1),
                    ("chr1", 100000, 200000, 1),
                ],
                "1": [
                    ("chr1", 100000, 200000, 1),
                    ("chr2", 0, 100000, 1),
                ],
            }

            def fake_pandepth(_binary, paf_path):
                index = Path(paf_path).stem
                write_depth(Path(paf_path).with_suffix(".win.stat.gz"), per_contig_depth[index])

            with mock.patch.object(full, "_run_pandepth", side_effect=fake_pandepth):
                records = full.run_full_assembly_stage21(
                    str(assembly), str(paf), str(prefix), "pandepth", 2
                )

            self.assertEqual([record["query_name"] for record in records], ["ctgA", "ctgB"])
            self.assertEqual([Path(record["paf_path"]).name for record in records], ["0.paf", "1.paf"])
            self.assertNotIn("tp:A:S", (prefix / "21_pat_depth" / "1.paf").read_text())
            self.assertIn("cg:Z:1000M", (prefix / "21_pat_depth" / "0.paf").read_text())
            paths_tsv = (prefix / full.FULL_ASSEMBLY_PATHS_TSV).read_text()
            self.assertIn("ctgC\t1000\tunmapped", paths_tsv)

            main_depth = root / "sample.win.stat.gz"
            write_depth(
                main_depth,
                [
                    ("chr1", 0, 100000, 10),
                    ("chr1", 100000, 200000, 20),
                    ("chr2", 0, 100000, 10),
                    ("chr2", 100000, 200000, 0),
                ],
            )
            censat = root / "censat.bed"
            censat.write_text("", encoding="ascii")

            matrix_info = full.build_full_assembly_matrix(
                str(prefix), str(main_depth), str(censat)
            )
            self.assertEqual(matrix_info["features"], 2)
            self.assertEqual(load_nclose_nodes(prefix), {})
            self.assertFalse((prefix / "nclose_chunk_data.pkl").exists())
            with h5py.File(prefix / "matrix.h5", "r") as handle:
                np.testing.assert_array_equal(
                    handle["A"][:],
                    np.asarray([[1, 1, 0, 0], [0, 1, 1, 0]], dtype=np.float32),
                )

            solve_info = full.solve_full_assembly_nnls(str(prefix), str(main_depth))
            np.testing.assert_allclose(solve_info["weights"], [10, 10], atol=1e-5)
            reference_fai = root / "reference.fa.fai"
            reference_fai.write_text(
                "chr1\t200000\t0\t0\t0\n"
                "chr2\t200000\t0\t0\t0\n"
                "chrM\t1000\t0\t0\t0\n",
                encoding="ascii",
            )
            cytobands = root / "cytobands.bed"
            cytobands.write_text(
                "chr1\t0\t200000\tp1\tgneg\n"
                "chr2\t0\t200000\tp1\tgneg\n",
                encoding="ascii",
            )
            telomeres = root / "telomeres.bed"
            telomeres.write_text("", encoding="ascii")
            save_pipeline_input(
                str(prefix),
                full_assembly_input=True,
                full_assembly_path=str(assembly),
                full_assembly_paf_path=str(paf),
                reference="hs1",
            )
            virtual_sky_info = full.draw_full_assembly_virtual_sky(
                str(prefix),
                str(main_depth),
                str(reference_fai),
                "toy-cell",
            )
            self.assertEqual(virtual_sky_info["displayed_contigs"], 2)
            coverage_info = full.draw_full_assembly_coverage(
                str(prefix),
                str(main_depth),
                str(reference_fai),
                str(cytobands),
            )
            self.assertLess(coverage_info["relative_error"], 1e-6)
            for filename in (
                "virtual_sky.png",
                "virtual_sky.pdf",
                "karyotype.txt",
                "total_cov.png",
                "total_cov.pdf",
                "report.txt",
                "matrix.h5",
            ):
                self.assertTrue((prefix / filename).is_file(), filename)
            self.assertFalse((prefix / "SV_call_result.vcf").exists())
            self.assertFalse((prefix / "SKYPE_result.bed").exists())

    def test_piece_merging_uses_query_order_and_strand(self):
        rows = [
            full.parse_paf_line(
                "q\t300\t100\t200\t-\tchr2\t1000\t100\t200\t100\t100\t60\tcg:Z:100M"
            ),
            full.parse_paf_line(
                "q\t300\t0\t100\t+\tchr1\t1000\t0\t100\t100\t100\t60\tcg:Z:100M"
            ),
            full.parse_paf_line(
                "q\t300\t200\t300\t-\tchr2\t1000\t0\t100\t100\t100\t60\tcg:Z:100M"
            ),
        ]
        pieces = full.paf_records_to_pieces(sorted(rows, key=lambda row: row.query_start))
        self.assertEqual(pieces, [("chr1", "+", 100), ("chr2", "-", 200)])

    def test_display_piece_filter_removes_short_noise_then_remerges(self):
        pieces = [
            ("chr1", "+", 2_000_000),
            ("chr13", "+", 50_000),
            ("chr1", "+", 3_000_000),
            ("chr2", "-", 1_500_000),
        ]
        self.assertEqual(
            full.filter_full_assembly_display_pieces(pieces),
            [("chr1", "+", 5_000_000), ("chr2", "-", 1_500_000)],
        )
        self.assertEqual(
            full.filter_full_assembly_display_pieces([("chr1", "+", 1000)]),
            [("chr1", "+", 1000)],
        )

    def test_alignasm_cs_is_converted_for_pandepth(self):
        self.assertEqual(
            full.cs_to_cigar(":10*ac:5+aa:3-tt"),
            "10M1X5M2I3M2D",
        )
        with self.assertRaisesRegex(ValueError, "Unsupported cs:Z operation"):
            full.cs_to_cigar(":10~gt5ag:3")

    def test_single_pipeline_entrypoint_owns_all_full_stages(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            prefix = root / "result"
            aligned_paf = root / "truth.hs1.aln.paf"
            assembly = root / "truth.fa"
            depth = root / "depth.win.stat.gz"
            censat = root / "censat.bed"
            telomere = root / "telomere.bed"
            reference_fai = root / "reference.fa.fai"
            cytobands = root / "cytobands.bed"
            for path in (
                aligned_paf,
                assembly,
                depth,
                censat,
                telomere,
                reference_fai,
                cytobands,
            ):
                path.write_text("placeholder", encoding="ascii")

            with mock.patch.object(
                full, "run_full_assembly_stage21", return_value=[]
            ) as stage21, mock.patch.object(
                full, "build_full_assembly_matrix", return_value={}
            ) as stage22, mock.patch.object(
                full, "solve_full_assembly_nnls", return_value={}
            ) as stage23, mock.patch.object(
                full, "draw_full_assembly_virtual_sky", return_value={}
            ) as stage30, mock.patch.object(
                full, "draw_full_assembly_coverage", return_value={}
            ) as stage31:
                full.run_full_assembly_pipeline(
                    aligned_paf_path=str(aligned_paf),
                    assembly_path=str(assembly),
                    main_depth_path=str(depth),
                    censat_bed_path=str(censat),
                    telomere_bed_path=str(telomere),
                    reference_fai_path=str(reference_fai),
                    reference_cytobands_path=str(cytobands),
                    prefix=str(prefix),
                    cell_line="toy",
                    pandepth_path="pandepth",
                    thread=2,
                    reference="hs1",
                )
            stage21.assert_called_once()
            stage22.assert_called_once()
            stage23.assert_called_once()
            stage30.assert_called_once()
            stage31.assert_called_once()

            raw_paf = root / "truth.hs1.paf"
            raw_paf.write_text("placeholder", encoding="ascii")
            with self.assertRaisesRegex(ValueError, "alignasm output"):
                full.run_full_assembly_pipeline(
                    aligned_paf_path=str(raw_paf),
                    assembly_path=str(assembly),
                    main_depth_path=str(depth),
                    censat_bed_path=str(censat),
                    telomere_bed_path=str(telomere),
                    reference_fai_path=str(reference_fai),
                    reference_cytobands_path=str(cytobands),
                    prefix=str(prefix),
                    cell_line="toy",
                    pandepth_path="pandepth",
                    thread=2,
                    reference="hs1",
                )


if __name__ == "__main__":
    unittest.main()
