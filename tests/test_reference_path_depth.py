from __future__ import annotations

import csv
import gzip
import pickle as pkl
import sys
import tempfile
import unittest
from pathlib import Path

import matplotlib


matplotlib.use("Agg")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

SKYPE_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKYPE_ROOT))

from reference_path_clustering import (  # noqa: E402
    ReferenceSpan,
    cluster_reference_paths,
)
from reference_path_depth import (  # noqa: E402
    PathDepthProfiles,
    build_cluster_depth_tracks,
    build_depth_bin_index,
    load_model_copy_number,
    observed_copy_number,
    summarize_path_span,
    summarize_region,
)
from reference_path_plotting import render_reference_path_clusters  # noqa: E402


WINDOW = 100_000


def depth_frame(rows):
    """Build a win.stat frame from (chrom, window index, mean depth) rows."""

    return pd.DataFrame([
        {
            "chr": chrom,
            "st": index * WINDOW + 1,
            "nd": (index + 1) * WINDOW,
            "length": WINDOW,
            "covsite": WINDOW,
            "totaldepth": int(mean_depth * WINDOW),
            "cov": 100.0,
            "meandepth": mean_depth,
        }
        for chrom, index, mean_depth in rows
    ])


def write_piece(folder: Path, key_int: int, rows) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    frame = depth_frame(rows)
    with gzip.open(folder / f"{key_int}.win.stat.gz", "wt") as handle:
        frame.to_csv(handle, sep="\t", header=False, index=False)


class DepthBinIndexTests(unittest.TestCase):
    def setUp(self) -> None:
        self.frame = depth_frame([
            ("chr1", 0, 20.0), ("chr1", 1, 40.0), ("chr1", 2, 10.0),
            ("chr2", 0, 30.0),
        ])
        self.index = build_depth_bin_index(self.frame)

    def test_one_based_windows_become_zero_based_half_open(self) -> None:
        self.assertEqual(self.index.starts.tolist(), [0, 100_000, 200_000, 0])
        self.assertEqual(self.index.ends.tolist(), [100_000, 200_000, 300_000, 100_000])
        self.assertEqual(self.index.row_by_key[("chr1", 100_001)], 1)

    def test_span_selects_only_intersecting_windows_of_that_chromosome(self) -> None:
        rows = self.index.rows_for_span("chr1", 100_000, 200_000)
        self.assertEqual(rows.tolist(), [1])
        # A span touching two windows keeps both, and never leaks to chr2.
        rows = self.index.rows_for_span("chr1", 150_000, 250_000)
        self.assertEqual(rows.tolist(), [1, 2])
        self.assertEqual(self.index.rows_for_span("chr3", 0, WINDOW).tolist(), [])

    def test_empty_and_inverted_spans_return_no_rows(self) -> None:
        self.assertEqual(self.index.rows_for_span("chr1", 500, 500).tolist(), [])
        self.assertEqual(self.index.rows_for_span("chr1", 900, 100).tolist(), [])

    def test_partial_window_overlap_is_measured_in_bases(self) -> None:
        rows = self.index.rows_for_span("chr1", 150_000, 250_000)
        lengths = self.index.overlap_lengths(rows, 150_000, 250_000)
        self.assertEqual(lengths.tolist(), [50_000, 50_000])

    def test_observed_copy_number_normalises_by_mean_depth(self) -> None:
        observed = observed_copy_number(self.frame, 20.0)
        self.assertEqual(observed.tolist(), [2.0, 4.0, 1.0, 3.0])


class PathDepthProfileTests(unittest.TestCase):
    def test_pieces_sum_and_unknown_windows_are_dropped(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir)
            frame = depth_frame([("chr1", 0, 20.0), ("chr1", 1, 20.0)])
            index = build_depth_bin_index(frame)
            folder = prefix / "21_pat_depth"
            write_piece(folder, 3, [("chr1", 0, 1.0), ("chr1", 1, 1.0)])
            # The second piece doubles back over window 0 and also reports a
            # window outside the main stat grid, which must be ignored.
            write_piece(folder, 4, [("chr1", 0, 1.0), ("chr9", 0, 1.0)])

            profiles = PathDepthProfiles(str(prefix), {"p": [3, 4]}, index)
            self.assertEqual(profiles.unit_profile("p").tolist(), [2.0, 1.0])
            self.assertEqual(profiles.missing_pieces, ())

    def test_missing_piece_is_reported_instead_of_raising(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir)
            index = build_depth_bin_index(depth_frame([("chr1", 0, 20.0)]))
            profiles = PathDepthProfiles(str(prefix), {"p": [7]}, index)
            self.assertEqual(profiles.unit_profile("p").tolist(), [0.0])
            self.assertEqual(profiles.missing_pieces, (7,))

    def test_unknown_path_raises(self) -> None:
        index = build_depth_bin_index(depth_frame([("chr1", 0, 20.0)]))
        profiles = PathDepthProfiles("/nonexistent", {}, index)
        with self.assertRaises(KeyError):
            profiles.unit_profile("missing")


class ModelCopyNumberTests(unittest.TestCase):
    def test_filtered_then_excluded_order_is_restored_to_window_order(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir)
            frame = depth_frame([
                ("chr1", 0, 20.0), ("chr1", 1, 20.0), ("chr2", 0, 20.0),
            ])
            index = build_depth_bin_index(frame)
            # Stage 22 puts the CENSAT/over-coverage window last; here that is
            # ("chr1", 1), so predict_B order is chr1:0, chr2:0, chr1:1.
            with open(prefix / "23_input.pkl", "wb") as handle:
                pkl.dump(([("chr1", 1), ("chr2", 1)], None, None, None, None), handle)
            np.save(prefix / "predict_B.npy", np.asarray([20.0, 60.0, 40.0]))

            model = load_model_copy_number(str(prefix), "", index, 20.0)

            self.assertEqual(model.tolist(), [2.0, 4.0, 6.0])

    def test_missing_inputs_degrade_to_none(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            index = build_depth_bin_index(depth_frame([("chr1", 0, 20.0)]))
            self.assertIsNone(load_model_copy_number(tmpdir, "", index, 20.0))

    def test_length_mismatch_degrades_to_none(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir)
            index = build_depth_bin_index(
                depth_frame([("chr1", 0, 20.0), ("chr1", 1, 20.0)])
            )
            with open(prefix / "23_input.pkl", "wb") as handle:
                pkl.dump(([("chr1", 1)], None, None, None, None), handle)
            np.save(prefix / "predict_B.npy", np.asarray([20.0]))
            self.assertIsNone(load_model_copy_number(str(prefix), "", index, 20.0))


class ClusterTrackTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tmpdir = tempfile.TemporaryDirectory()
        prefix = Path(self.tmpdir.name)
        self.prefix = prefix
        # chr1 has four windows; the observed depth is 3N over the shared half.
        self.frame = depth_frame([
            ("chr1", 0, 40.0), ("chr1", 1, 60.0),
            ("chr1", 2, 60.0), ("chr1", 3, 20.0),
        ])
        self.index = build_depth_bin_index(self.frame)
        folder = prefix / "21_pat_depth"
        write_piece(folder, 1, [("chr1", 1, 1.0), ("chr1", 2, 1.0)])
        write_piece(folder, 2, [("chr1", 0, 1.0), ("chr1", 1, 1.0), ("chr1", 2, 1.0)])
        self.profiles = PathDepthProfiles(
            str(prefix), {"/run/a.paf": [1], "/run/b.paf": [2]}, self.index
        )
        self.metadata = {
            11: {"location": "/run/a.paf", "depth_n": 2.0, "raw_depth": 40.0,
                 "path_column": 11, "nclose_count": 1},
            12: {"location": "/run/b.paf", "depth_n": 1.0, "raw_depth": 20.0,
                 "path_column": 12, "nclose_count": 1},
        }
        self.spans = {
            11: (ReferenceSpan("chr1", 100_000, 300_000, "+"),),
            12: (ReferenceSpan("chr1", 0, 300_000, "+"),),
        }
        self.tracks = build_cluster_depth_tracks(
            members=(11, 12),
            path_metadata=self.metadata,
            path_spans=self.spans,
            bin_index=self.index,
            observed_cn=observed_copy_number(self.frame, 20.0),
            model_cn=np.asarray([1.5, 3.5, 3.5, 0.0]),
            profiles=self.profiles,
        )

    def tearDown(self) -> None:
        self.tmpdir.cleanup()

    def test_member_track_is_the_matrix_column_times_fitted_copy_number(self) -> None:
        self.assertEqual(self.tracks.member_cn[11].tolist(), [0.0, 2.0, 2.0, 0.0])
        self.assertEqual(self.tracks.member_cn[12].tolist(), [1.0, 1.0, 1.0, 0.0])
        self.assertEqual(self.tracks.cluster_cn.tolist(), [1.0, 3.0, 3.0, 0.0])

    def test_other_paths_never_go_negative(self) -> None:
        self.assertEqual(self.tracks.other_cn.tolist(), [0.5, 0.5, 0.5, 0.0])

    def test_cluster_rows_cover_the_union_of_member_spans(self) -> None:
        self.assertEqual(self.tracks.cluster_rows.tolist(), [0, 1, 2])

    def test_region_summary_is_length_weighted_over_partial_windows(self) -> None:
        summary = summarize_region(
            label="R01", chrom="chr1", start=100_000, end=300_000,
            tracks=self.tracks, carriers=(11, 12),
        )
        self.assertAlmostEqual(summary.observed_cn, 6.0)
        self.assertAlmostEqual(summary.cluster_cn, 3.0)
        self.assertAlmostEqual(summary.member_cn[11], 2.0)
        self.assertAlmostEqual(summary.member_cn[12], 1.0)
        self.assertAlmostEqual(summary.other_cn, 0.5)
        self.assertAlmostEqual(summary.observed_sd, 0.0)

    def test_half_window_region_still_weights_by_overlapping_bases(self) -> None:
        summary = summarize_region(
            label="R02", chrom="chr1", start=150_000, end=250_000,
            tracks=self.tracks,
        )
        self.assertAlmostEqual(summary.observed_cn, 6.0)

    def test_path_footprint_summary_reports_the_merged_span_length(self) -> None:
        summary = summarize_path_span(
            path=11, spans=self.spans[11], tracks=self.tracks
        )
        self.assertEqual(summary.end, 200_000)
        self.assertAlmostEqual(summary.member_cn[11], 2.0)
        self.assertAlmostEqual(summary.observed_cn, 6.0)


class DepthFigureRenderTests(unittest.TestCase):
    def test_every_depth_mode_writes_figures_and_the_audit_table(self) -> None:
        frame = depth_frame([("chr1", index, 60.0) for index in range(6)])
        bin_index = build_depth_bin_index(frame)
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir)
            folder = prefix / "21_pat_depth"
            write_piece(folder, 1, [("chr1", index, 1.0) for index in range(4)])
            write_piece(folder, 2, [("chr1", index, 1.0) for index in range(2, 6)])
            profiles = PathDepthProfiles(
                str(prefix), {"/run/a.paf": [1], "/run/b.paf": [2]}, bin_index
            )
            metadata = {
                11: {"location": "/run/a.paf", "depth_n": 1.5, "raw_depth": 45.0,
                     "path_column": 11, "nclose_count": 2},
                12: {"location": "/run/b.paf", "depth_n": 1.5, "raw_depth": 45.0,
                     "path_column": 12, "nclose_count": 2},
            }
            path_spans = {
                11: (ReferenceSpan("chr1", 0, 400_000, "+"),),
                12: (ReferenceSpan("chr1", 200_000, 600_000, "+"),),
            }
            clusters, singletons = cluster_reference_paths(
                path_spans, minimum_overlap=100_000
            )
            tracks = {
                clusters[0].cluster_id: build_cluster_depth_tracks(
                    members=clusters[0].members,
                    path_metadata=metadata,
                    path_spans=path_spans,
                    bin_index=bin_index,
                    observed_cn=observed_copy_number(frame, 20.0),
                    model_cn=np.full(6, 3.0),
                    profiles=profiles,
                )
            }

            result = render_reference_path_clusters(
                output_prefix=str(prefix),
                cell_line="TEST",
                clusters=clusters,
                singletons=singletons,
                path_spans=path_spans,
                path_metadata=metadata,
                excluded_full_reference=[],
                minimum_overlap=100_000,
                full_reference_ratio=0.9,
                modes=("plain", "lane", "profile", "attribution"),
                cluster_tracks=tracks,
            )

            self.assertEqual(
                [figure["mode"] for figure in result["figures"]],
                ["plain", "lane", "profile", "attribution"],
            )
            for figure in result["figures"]:
                self.assertTrue(Path(figure["png_path"]).is_file())
                self.assertTrue(Path(figure["pdf_path"]).is_file())
            with open(result["depth_path"], newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            block_rows = [row for row in rows if row["region_type"] == "shared_block"]
            self.assertEqual({row["region"] for row in block_rows}, {"R01"})
            self.assertAlmostEqual(float(block_rows[0]["observed_CN"]), 6.0)
            self.assertAlmostEqual(float(block_rows[0]["cluster_CN"]), 3.0)
            self.assertEqual(
                {row["path_column"] for row in block_rows}, {"11", "12"}
            )

    def test_depth_modes_are_skipped_when_no_tracks_are_available(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            path_spans = {
                11: (ReferenceSpan("chr1", 0, 400_000, "+"),),
                12: (ReferenceSpan("chr1", 200_000, 600_000, "+"),),
            }
            metadata = {
                column: {"location": f"/run/{column}.paf", "depth_n": 1.0,
                         "raw_depth": 30.0, "path_column": column,
                         "nclose_count": 1}
                for column in path_spans
            }
            clusters, singletons = cluster_reference_paths(
                path_spans, minimum_overlap=100_000
            )

            result = render_reference_path_clusters(
                output_prefix=tmpdir,
                cell_line="TEST",
                clusters=clusters,
                singletons=singletons,
                path_spans=path_spans,
                path_metadata=metadata,
                excluded_full_reference=[],
                minimum_overlap=100_000,
                full_reference_ratio=0.9,
                modes=("plain", "lane", "profile", "attribution"),
                cluster_tracks={},
            )

            self.assertEqual([figure["mode"] for figure in result["figures"]], ["plain"])
            self.assertNotIn("depth_path", result)


if __name__ == "__main__":
    unittest.main()
