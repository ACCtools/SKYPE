#!/usr/bin/env python3
"""Serve an interactive circular viewer for SKYPE depth results."""

from __future__ import annotations

import argparse
import json
import os
import pickle
from dataclasses import dataclass
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any

import numpy as np


BIN_SIZE = 100_000
HOST = "0.0.0.0"
CANONICAL_CHROMS = tuple(
    [f"chr{index}" for index in range(1, 23)] + ["chrX", "chrY"]
)
REFERENCE_FILES = {
    "hs1": "chm13v2.0.fa.fai",
    "hg38": "hg38.fa.fai",
}
REQUIRED_FILES = ("23_input.pkl", "B.npy", "predict_B.npy")
STATIC_DIR = Path(__file__).resolve().with_name("skype_view_static")
PUBLIC_DATA_DIR = Path(__file__).resolve().with_name("public_data")


class ViewerDataError(ValueError):
    """Raised when a SKYPE result cannot satisfy the viewer data contract."""


@dataclass(frozen=True)
class ReferenceLayout:
    name: str
    chromosomes: tuple[tuple[str, int], ...]
    total_bins: int


@dataclass(frozen=True)
class ViewerData:
    run_dir: Path
    sample: str
    reference: str
    clean_bins: int
    total_bins: int
    median_depth: float
    scale_max: float
    payload: dict[str, Any]
    json_bytes: bytes


def _load_reference_layouts() -> tuple[ReferenceLayout, ...]:
    layouts = []
    for reference, filename in REFERENCE_FILES.items():
        path = PUBLIC_DATA_DIR / filename
        if not path.is_file():
            raise ViewerDataError(f"Reference FAI is missing: {path}")

        lengths: dict[str, int] = {}
        with path.open("rt", encoding="utf-8") as handle:
            for line in handle:
                fields = line.split()
                if len(fields) >= 2 and fields[0] in CANONICAL_CHROMS:
                    lengths[fields[0]] = int(fields[1])

        missing = [chrom for chrom in CANONICAL_CHROMS if chrom not in lengths]
        if missing:
            raise ViewerDataError(
                f"Reference FAI {path} is missing canonical chromosomes: "
                + ", ".join(missing)
            )
        chromosomes = tuple((chrom, lengths[chrom]) for chrom in CANONICAL_CHROMS)
        total_bins = sum((length + BIN_SIZE - 1) // BIN_SIZE for _, length in chromosomes)
        layouts.append(ReferenceLayout(reference, chromosomes, total_bins))
    return tuple(layouts)


def _select_reference(total_bins: int) -> ReferenceLayout:
    matches = [
        layout for layout in _load_reference_layouts() if layout.total_bins == total_bins
    ]
    if len(matches) == 1:
        return matches[0]
    supported = ", ".join(
        f"{layout.name}={layout.total_bins}"
        for layout in _load_reference_layouts()
    )
    if matches:
        raise ViewerDataError(
            f"Depth-vector length {total_bins} matches multiple references"
        )
    raise ViewerDataError(
        f"Unsupported depth-vector length {total_bins}; expected {supported} bins"
    )


def _load_clean_coords(metadata_path: Path) -> list[tuple[str, int]]:
    try:
        with metadata_path.open("rb") as handle:
            metadata = pickle.load(handle)
    except (OSError, pickle.PickleError, EOFError) as exc:
        raise ViewerDataError(f"Could not read {metadata_path}: {exc}") from exc

    if not isinstance(metadata, dict):
        raise ViewerDataError(f"{metadata_path} must contain a metadata dictionary")
    raw_coords = metadata.get("chr_filt_st_list")
    if not isinstance(raw_coords, list):
        raise ViewerDataError(
            f"{metadata_path} does not contain chr_filt_st_list"
        )

    coords = []
    for index, coord in enumerate(raw_coords):
        if not isinstance(coord, (tuple, list)) or len(coord) != 2:
            raise ViewerDataError(
                f"Malformed clean coordinate at index {index}: {coord!r}"
            )
        chrom, start = str(coord[0]), int(coord[1])
        coords.append((chrom, start))
    return coords


def _load_depth_vector(path: Path, label: str) -> np.ndarray:
    try:
        values = np.load(path, allow_pickle=False)
    except (OSError, ValueError) as exc:
        raise ViewerDataError(f"Could not read {label} array {path}: {exc}") from exc
    if values.ndim != 1:
        raise ViewerDataError(f"{label} must be one-dimensional, got {values.shape}")
    if not np.issubdtype(values.dtype, np.number):
        raise ViewerDataError(f"{label} must contain numeric values")
    if not np.isfinite(values).all():
        raise ViewerDataError(f"{label} contains NaN or infinite values")
    return values


def load_viewer_data(run_dir: str | os.PathLike[str]) -> ViewerData:
    """Load, validate, and serialize one SKYPE depth result directory."""

    resolved = Path(run_dir).expanduser().resolve()
    if not resolved.is_dir():
        raise ViewerDataError(f"SKYPE result directory not found: {resolved}")

    missing = [name for name in REQUIRED_FILES if not (resolved / name).is_file()]
    if missing:
        raise ViewerDataError(
            f"Missing required SKYPE result file(s) in {resolved}: "
            + ", ".join(missing)
        )

    coords = _load_clean_coords(resolved / "23_input.pkl")
    observed_all = _load_depth_vector(resolved / "B.npy", "B.npy")
    predicted_all = _load_depth_vector(
        resolved / "predict_B.npy", "predict_B.npy"
    )
    if len(observed_all) != len(predicted_all):
        raise ViewerDataError(
            "B.npy and predict_B.npy length mismatch: "
            f"{len(observed_all)} != {len(predicted_all)}"
        )
    if len(coords) > len(observed_all):
        raise ViewerDataError(
            "Clean coordinate count exceeds depth-vector length: "
            f"{len(coords)} > {len(observed_all)}"
        )

    layout = _select_reference(len(observed_all))
    chromosome_lengths = dict(layout.chromosomes)
    seen: set[tuple[str, int]] = set()
    grouped: dict[str, list[tuple[int, float, float, float]]] = {
        chrom: [] for chrom, _ in layout.chromosomes
    }
    for index, (chrom, start) in enumerate(coords):
        if chrom not in chromosome_lengths:
            raise ViewerDataError(
                f"Clean coordinate uses unsupported chromosome at index {index}: {chrom}"
            )
        if start < 1 or start > chromosome_lengths[chrom]:
            raise ViewerDataError(
                f"Clean coordinate is outside {chrom} at index {index}: {start}"
            )
        key = (chrom, start)
        if key in seen:
            raise ViewerDataError(f"Duplicate clean coordinate: {chrom}:{start}")
        seen.add(key)
        observed = float(observed_all[index])
        predicted_raw = float(predicted_all[index])
        grouped[chrom].append(
            (start, observed, max(predicted_raw, 0.0), predicted_raw)
        )

    chromosomes = []
    for chrom, length in layout.chromosomes:
        rows = sorted(grouped[chrom], key=lambda row: row[0])
        chromosomes.append(
            {
                "name": chrom,
                "length": length,
                "starts": [row[0] for row in rows],
                "reference_depth": [row[1] for row in rows],
                "predicted_depth": [row[2] for row in rows],
                "predicted_depth_raw": [row[3] for row in rows],
            }
        )

    median_depth = float(np.median(observed_all))
    if median_depth <= 0:
        raise ViewerDataError(
            f"B.npy median depth must be positive, got {median_depth:g}"
        )
    payload = {
        "sample": resolved.name,
        "reference": layout.name,
        "bin_size": BIN_SIZE,
        "clean_bins": len(coords),
        "total_bins": len(observed_all),
        "ignored_bins": len(observed_all) - len(coords),
        "median_depth": median_depth,
        "scale_max": 3.0 * median_depth,
        "chromosomes": chromosomes,
    }
    json_bytes = json.dumps(
        payload,
        ensure_ascii=False,
        allow_nan=False,
        separators=(",", ":"),
    ).encode("utf-8")
    return ViewerData(
        run_dir=resolved,
        sample=resolved.name,
        reference=layout.name,
        clean_bins=len(coords),
        total_bins=len(observed_all),
        median_depth=median_depth,
        scale_max=3.0 * median_depth,
        payload=payload,
        json_bytes=json_bytes,
    )


class ViewerRequestHandler(BaseHTTPRequestHandler):
    """Serve only the fixed viewer assets and the preloaded depth payload."""

    api_bytes = b"{}"
    static_assets: dict[str, tuple[str, bytes]] = {}

    def do_GET(self) -> None:  # noqa: N802 - BaseHTTPRequestHandler API
        path = self.path.partition("?")[0]
        if path == "/api/depth":
            self._send("application/json; charset=utf-8", self.api_bytes)
            return
        asset = self.static_assets.get(path)
        if asset is None:
            self.send_error(404, "Not found")
            return
        self._send(asset[0], asset[1])

    def _send(self, content_type: str, body: bytes) -> None:
        self.send_response(200)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, fmt: str, *args: object) -> None:
        print(f"[skype_view] {self.address_string()} - {fmt % args}")


def _read_static_assets() -> dict[str, tuple[str, bytes]]:
    routes = {
        "/": ("text/html; charset=utf-8", "index.html"),
        "/app.js": ("text/javascript; charset=utf-8", "app.js"),
        "/style.css": ("text/css; charset=utf-8", "style.css"),
    }
    assets = {}
    for route, (content_type, filename) in routes.items():
        path = STATIC_DIR / filename
        if not path.is_file():
            raise ViewerDataError(f"Viewer asset is missing: {path}")
        assets[route] = (content_type, path.read_bytes())
    return assets


def create_server(
    data: ViewerData,
    port: int,
    host: str = HOST,
) -> ThreadingHTTPServer:
    """Create a viewer server without entering its event loop."""

    if not 0 <= port <= 65_535:
        raise ViewerDataError(f"Port must be between 0 and 65535, got {port}")
    handler = type(
        "BoundViewerRequestHandler",
        (ViewerRequestHandler,),
        {
            "api_bytes": data.json_bytes,
            "static_assets": _read_static_assets(),
        },
    )
    return ThreadingHTTPServer((host, port), handler)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Serve an interactive circular viewer for a SKYPE result"
    )
    parser.add_argument("run_dir", help="SKYPE result directory")
    parser.add_argument(
        "--host",
        default=HOST,
        help=f"bind address (default: {HOST})",
    )
    parser.add_argument(
        "--port",
        type=int,
        default=3000,
        help="server port (default: 3000)",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        data = load_viewer_data(args.run_dir)
        server = create_server(data, args.port, args.host)
    except ViewerDataError as exc:
        parser.error(str(exc))
    except OSError as exc:
        parser.error(f"Could not bind {args.host}:{args.port}: {exc}")

    actual_port = server.server_address[1]
    local_host = "127.0.0.1" if args.host == "0.0.0.0" else args.host
    local_url = f"http://{local_host}:{actual_port}"
    print(
        f"SKYPE depth viewer: {data.sample} "
        f"({data.clean_bins:,}/{data.total_bins:,} clean bins, {data.reference})"
    )
    print(f"Listening on {args.host}:{actual_port}")
    print(f"Open {local_url}")
    if args.host == "0.0.0.0":
        print(f"Remote URL: http://<server-ip>:{actual_port}")
    print("Press Ctrl-C to stop.")

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopping SKYPE depth viewer.")
    finally:
        server.server_close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
