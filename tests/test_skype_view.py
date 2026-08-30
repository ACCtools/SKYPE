import json
import pickle
import socket
import tempfile
import threading
import unittest
import urllib.error
import urllib.request
from pathlib import Path

import numpy as np

import skype_view


class SkypeViewDataTest(unittest.TestCase):
    def setUp(self):
        self.tempdir = tempfile.TemporaryDirectory()
        self.run_dir = Path(self.tempdir.name) / "example_run"
        self.run_dir.mkdir()

    def tearDown(self):
        self.tempdir.cleanup()

    def write_run(
        self,
        total_bins=31_185,
        coords=None,
        observed=None,
        predicted=None,
    ):
        if coords is None:
            coords = [("chr1", 1), ("chr1", 100_001), ("chrY", 62_400_001)]
        if observed is None:
            observed = np.full(total_bins, 10, dtype=np.float32)
        if predicted is None:
            predicted = np.full(total_bins, 12, dtype=np.float32)
        with (self.run_dir / "23_input.pkl").open("wb") as handle:
            pickle.dump({"chr_filt_st_list": coords}, handle)
        np.save(self.run_dir / "B.npy", observed)
        np.save(self.run_dir / "predict_B.npy", predicted)
        return self.run_dir

    def test_reference_layout_selection(self):
        self.assertEqual(skype_view._select_reference(31_185).name, "hs1")
        self.assertEqual(skype_view._select_reference(30_894).name, "hg38")

    def test_clean_values_are_mapped_and_prediction_is_clamped(self):
        observed = np.full(31_185, 10, dtype=np.float32)
        predicted = np.full(31_185, 12, dtype=np.float32)
        observed[:3] = [4.25, 5.5, 6.75]
        predicted[:3] = [3.5, -2.0, 9.25]
        run_dir = self.write_run(observed=observed, predicted=predicted)

        data = skype_view.load_viewer_data(run_dir)
        payload = json.loads(data.json_bytes)

        self.assertEqual(data.reference, "hs1")
        self.assertEqual(data.clean_bins, 3)
        self.assertEqual(data.total_bins, 31_185)
        chr1 = payload["chromosomes"][0]
        self.assertEqual(chr1["starts"], [1, 100_001])
        self.assertEqual(chr1["reference_depth"], [4.25, 5.5])
        self.assertEqual(chr1["predicted_depth"], [3.5, 0.0])
        self.assertEqual(chr1["predicted_depth_raw"], [3.5, -2.0])

    def test_missing_required_file_is_rejected(self):
        self.run_dir.mkdir(exist_ok=True)
        with self.assertRaisesRegex(skype_view.ViewerDataError, "Missing required"):
            skype_view.load_viewer_data(self.run_dir)

    def test_array_length_mismatch_is_rejected(self):
        observed = np.ones(31_185, dtype=np.float32)
        predicted = np.ones(31_184, dtype=np.float32)
        run_dir = self.write_run(observed=observed, predicted=predicted)
        with self.assertRaisesRegex(skype_view.ViewerDataError, "length mismatch"):
            skype_view.load_viewer_data(run_dir)

    def test_unsupported_layout_is_rejected(self):
        run_dir = self.write_run(
            total_bins=20,
            coords=[("chr1", 1)],
            observed=np.ones(20, dtype=np.float32),
            predicted=np.ones(20, dtype=np.float32),
        )
        with self.assertRaisesRegex(skype_view.ViewerDataError, "Unsupported"):
            skype_view.load_viewer_data(run_dir)

    def test_server_routes_and_port_collision(self):
        data = skype_view.load_viewer_data(self.write_run())
        server = skype_view.create_server(data, 0)
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        base_url = f"http://127.0.0.1:{server.server_address[1]}"
        try:
            with urllib.request.urlopen(base_url + "/", timeout=5) as response:
                self.assertEqual(response.status, 200)
                self.assertIn(b"SKYPE Interactive Depth Viewer", response.read())
            with urllib.request.urlopen(base_url + "/api/depth", timeout=5) as response:
                payload = json.load(response)
                self.assertEqual(payload["sample"], "example_run")
                self.assertEqual(payload["clean_bins"], 3)
            with self.assertRaises(urllib.error.HTTPError) as caught:
                urllib.request.urlopen(base_url + "/not-a-route", timeout=5)
            self.assertEqual(caught.exception.code, 404)
        finally:
            server.shutdown()
            server.server_close()
            thread.join(timeout=5)

        with socket.socket() as occupied:
            occupied.bind((skype_view.HOST, 0))
            port = occupied.getsockname()[1]
            occupied.listen(1)
            with self.assertRaises(OSError):
                skype_view.create_server(data, port)


class Hcc1937WorkspaceTest(unittest.TestCase):
    def test_workspace_sample_mapping(self):
        workspace = Path(__file__).resolve().parents[3]
        run_dir = workspace / "skype" / "HCC1937_21_53_34"
        if not run_dir.is_dir():
            self.skipTest("workspace HCC1937 sample is unavailable")

        data = skype_view.load_viewer_data(run_dir)
        self.assertEqual(data.reference, "hs1")
        self.assertEqual(data.clean_bins, 26_652)
        self.assertEqual(data.total_bins, 31_185)

        with (run_dir / "23_input.pkl").open("rb") as handle:
            coords = pickle.load(handle)["chr_filt_st_list"]
        observed = np.load(run_dir / "B.npy")
        predicted = np.load(run_dir / "predict_B.npy")
        lookup = {
            (chromosome["name"], start): (reference, prediction, raw)
            for chromosome in data.payload["chromosomes"]
            for start, reference, prediction, raw in zip(
                chromosome["starts"],
                chromosome["reference_depth"],
                chromosome["predicted_depth"],
                chromosome["predicted_depth_raw"],
            )
        }
        for index in (0, len(coords) // 2, len(coords) - 1):
            reference, prediction, raw = lookup[coords[index]]
            self.assertEqual(reference, float(observed[index]))
            self.assertEqual(prediction, max(float(predicted[index]), 0.0))
            self.assertEqual(raw, float(predicted[index]))


if __name__ == "__main__":
    unittest.main()
