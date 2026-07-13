from __future__ import annotations

import gzip
import json
import pickle
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SKYPE_ROOT = Path(__file__).resolve().parents[1]


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
                {"DEL", "INS"},
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
