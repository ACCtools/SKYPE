"""Standalone full-assembly SKYPE pipeline and reusable rendering adapters.

The full assembly is aligned once by the outer ACCtools orchestrator.  This
module treats every mapped FASTA record as one matrix feature and deliberately
does not fabricate native stage-01/10/11 graph paths.
"""

from __future__ import annotations

import csv
import argparse
import json
import logging
import math
import os
import pickle
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Iterable, Sequence

import h5py
import numpy as np
import pandas as pd
from scipy.optimize import nnls

from nclose_tracking import (
    initialise_filter_status,
    save_event_catalog,
    save_filter_status,
    save_path_usage,
)
from skype_utils import (
    load_pipeline_input,
    pipeline_input_is_full_assembly,
    save_nclose_nodes,
    save_pipeline_input,
)
from virtual_sky_plotting import render_karyotype_diagram


FULL_ASSEMBLY_PATHS_PKL = "full_assembly_paths.pkl"
FULL_ASSEMBLY_PATHS_TSV = "full_assembly_paths.tsv"
FULL_ASSEMBLY_SUMMARY_JSON = "full_assembly_summary.json"
FULL_ASSEMBLY_WEIGHTS_TSV = "full_assembly_weights.tsv"
FULL_ASSEMBLY_DEPTH_LAYOUT_PKL = "full_assembly_depth_layout.pkl"
CONTIG_PATH_DATA_PKL = "contig_pat_vec_data.pkl"
DEPTH_WINDOW = 100_000
PANDEPTH_RETRIES = 3
ABS_MAX_COVERAGE_RATIO = 3.0
FULL_ASSEMBLY_PLOT_MIN_PIECE = 1_000_000

DEPTH_COLUMNS = [
    "chr",
    "st",
    "nd",
    "length",
    "covsite",
    "totaldepth",
    "cov",
    "meandepth",
]

@dataclass(frozen=True)
class PafRecord:
    query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    matches: int
    block_length: int
    mapq: int
    tags: tuple[str, ...]
    raw_line: str

    @property
    def is_secondary(self) -> bool:
        return "tp:A:S" in self.tags

    @property
    def has_alignment_trace(self) -> bool:
        return any(
            tag.startswith(("cs:Z:", "cg:Z:"))
            for tag in self.tags
        )


def parse_paf_line(line: str, source: str = "PAF", line_number: int = 0) -> PafRecord:
    values = line.rstrip("\r\n").split("\t")
    location = f"{source}:{line_number}" if line_number else source
    if len(values) < 12:
        raise ValueError(f"{location}: expected at least 12 PAF fields, found {len(values)}")
    try:
        numbers = [int(values[index]) for index in (1, 2, 3, 6, 7, 8, 9, 10, 11)]
    except ValueError as exc:
        raise ValueError(f"{location}: invalid numeric PAF field") from exc
    qlen, qst, qend, tlen, tst, tend, matches, block_length, mapq = numbers
    if values[4] not in {"+", "-"}:
        raise ValueError(f"{location}: invalid PAF strand {values[4]!r}")
    if not 0 <= qst <= qend <= qlen:
        raise ValueError(f"{location}: invalid query interval {qst}-{qend}/{qlen}")
    if not 0 <= tst <= tend <= tlen:
        raise ValueError(f"{location}: invalid target interval {tst}-{tend}/{tlen}")
    return PafRecord(
        query_name=values[0],
        query_length=qlen,
        query_start=qst,
        query_end=qend,
        strand=values[4],
        target_name=values[5],
        target_length=tlen,
        target_start=tst,
        target_end=tend,
        matches=matches,
        block_length=block_length,
        mapq=mapq,
        tags=tuple(values[12:]),
        raw_line=line.rstrip("\r\n"),
    )


def cs_to_cigar(cs_tag: str) -> str:
    """Convert minimap2's short ``cs`` payload to a PanDepth ``cg`` CIGAR."""

    pattern = re.compile(r"(:\d+|\*[A-Za-z]{2}|[+\-][A-Za-z]+)")
    operations = []
    pending_op = "M"
    pending_length = 0

    def flush() -> None:
        nonlocal pending_length
        if pending_length > 0:
            operations.append(f"{pending_length}{pending_op}")
            pending_length = 0

    consumed = 0
    for match in pattern.finditer(cs_tag):
        if match.start() != consumed:
            raise ValueError(f"Unsupported cs:Z operation near {cs_tag[consumed:]!r}")
        token = match.group(0)
        consumed = match.end()
        op = token[0]
        if op == ":":
            length = int(token[1:])
            if pending_op == "M":
                pending_length += length
            else:
                flush()
                pending_op = "M"
                pending_length = length
        elif op in {"+", "-"}:
            flush()
            operations.append(f"{len(token) - 1}{'I' if op == '+' else 'D'}")
            pending_op = "M"
        elif op == "*":
            if pending_op == "X":
                pending_length += 1
            else:
                flush()
                pending_op = "X"
                pending_length = 1
    if consumed != len(cs_tag):
        raise ValueError(f"Unsupported cs:Z operation near {cs_tag[consumed:]!r}")
    flush()
    if not operations:
        raise ValueError("Cannot derive cg:Z CIGAR from an empty/unsupported cs:Z tag")
    return "".join(operations)


def pandepth_paf_line(record: PafRecord) -> str:
    """Return an alignment row with the cg tag required by PanDepth."""

    if any(tag.startswith("cg:Z:") for tag in record.tags):
        return record.raw_line
    cs_tag = next(
        (tag[len("cs:Z:"):] for tag in record.tags if tag.startswith("cs:Z:")),
        None,
    )
    if cs_tag is None:
        raise ValueError(f"{record.query_name}: alignment has neither cs:Z nor cg:Z")
    return f"{record.raw_line}\tcg:Z:{cs_to_cigar(cs_tag)}"


def read_fasta_contigs(fasta_path: str | os.PathLike) -> list[tuple[str, int]]:
    """Return unique FASTA records in source order, using an FAI when present."""

    fasta_path = os.fspath(fasta_path)
    fai_path = f"{fasta_path}.fai"
    records: list[tuple[str, int]] = []
    if os.path.isfile(fai_path):
        with open(fai_path, "rt", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                values = line.rstrip("\r\n").split("\t")
                if len(values) < 2:
                    raise ValueError(f"{fai_path}:{line_number}: malformed FAI row")
                records.append((values[0], int(values[1])))
    else:
        name = None
        length = 0
        with open(fasta_path, "rt", encoding="ascii") as handle:
            for line_number, line in enumerate(handle, start=1):
                if line.startswith(">"):
                    if name is not None:
                        records.append((name, length))
                    name = line[1:].strip().split(maxsplit=1)[0]
                    if not name:
                        raise ValueError(f"{fasta_path}:{line_number}: empty FASTA record name")
                    length = 0
                else:
                    if name is None and line.strip():
                        raise ValueError(f"{fasta_path}:{line_number}: sequence before FASTA header")
                    length += len(line.strip())
        if name is not None:
            records.append((name, length))

    if not records:
        raise ValueError(f"Full assembly FASTA has no records: {fasta_path}")
    names = [name for name, _ in records]
    duplicates = sorted(name for name, count in Counter(names).items() if count > 1)
    if duplicates:
        raise ValueError(f"Full assembly FASTA has duplicate record names: {duplicates[:10]}")
    return records


def _union_length(intervals: Iterable[tuple[int, int]]) -> int:
    ordered = sorted((int(st), int(nd)) for st, nd in intervals if nd > st)
    if not ordered:
        return 0
    total = 0
    current_st, current_nd = ordered[0]
    for st, nd in ordered[1:]:
        if st <= current_nd:
            current_nd = max(current_nd, nd)
        else:
            total += current_nd - current_st
            current_st, current_nd = st, nd
    return total + current_nd - current_st


def _run_pandepth(pandepth_path: str, paf_path: str) -> None:
    output_prefix = os.path.splitext(paf_path)[0]
    command = [
        pandepth_path,
        "--fast-paf-window",
        "-w",
        str(DEPTH_WINDOW),
        "-t",
        "1",
        "-i",
        paf_path,
        "-o",
        output_prefix,
    ]
    last_returncode = None
    for attempt in range(1, PANDEPTH_RETRIES + 1):
        result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE, text=True)
        last_returncode = result.returncode
        if result.returncode == 0:
            return
        logging.warning(
            "PanDepth failed for %s (attempt %d/%d): %s",
            paf_path,
            attempt,
            PANDEPTH_RETRIES,
            result.stderr.strip(),
        )
    raise RuntimeError(f"PanDepth failed for {paf_path} with exit code {last_returncode}")


def run_full_assembly_stage21(
    assembly_path: str,
    paf_path: str,
    prefix: str,
    pandepth_path: str,
    thread: int,
) -> list[dict]:
    """Split a whole-assembly PAF by FASTA record and run PanDepth per record."""

    assembly_path = os.path.abspath(assembly_path)
    paf_path = os.path.abspath(paf_path)
    prefix = os.path.abspath(prefix)
    contigs = read_fasta_contigs(assembly_path)
    contig_index = {name: (index, length) for index, (name, length) in enumerate(contigs)}

    output_folder = os.path.join(prefix, "21_pat_depth")
    if os.path.isdir(output_folder):
        shutil.rmtree(output_folder)
    os.makedirs(output_folder, exist_ok=True)

    alignment_count = Counter()
    query_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    temporary_output_by_query = {
        name: os.path.join(output_folder, f".fasta_{index}.paf.tmp")
        for index, (name, _) in enumerate(contigs)
    }

    with open(paf_path, "rt", encoding="utf-8") as paf_handle:
        for line_number, line in enumerate(paf_handle, start=1):
            if not line.strip():
                continue
            record = parse_paf_line(line, paf_path, line_number)
            if record.query_name not in contig_index:
                raise ValueError(
                    f"{paf_path}:{line_number}: query {record.query_name!r} is absent from {assembly_path}"
                )
            _, expected_length = contig_index[record.query_name]
            if record.query_length != expected_length:
                raise ValueError(
                    f"{paf_path}:{line_number}: query length mismatch for {record.query_name}: "
                    f"PAF={record.query_length}, FASTA={expected_length}"
                )
            if record.is_secondary:
                continue
            if not record.has_alignment_trace:
                raise ValueError(
                    f"{paf_path}:{line_number}: primary alignment is missing a cs:Z/cg:Z trace tag"
                )
            with open(temporary_output_by_query[record.query_name], "at", encoding="utf-8") as out_handle:
                out_handle.write(pandepth_paf_line(record) + "\n")
            alignment_count[record.query_name] += 1
            query_intervals[record.query_name].append((record.query_start, record.query_end))

    mapped_records = []
    all_rows = []
    for fasta_index, (query_name, query_length) in enumerate(contigs):
        count = int(alignment_count[query_name])
        aligned_bases = _union_length(query_intervals[query_name])
        matrix_index = len(mapped_records) if count else None
        final_paf_path = (
            os.path.join(output_folder, f"{matrix_index}.paf") if count else None
        )
        if count:
            os.replace(temporary_output_by_query[query_name], final_paf_path)
        status = "mapped" if count else "unmapped"
        row = {
            "fasta_index": fasta_index,
            "matrix_index": matrix_index,
            "query_name": query_name,
            "query_length": query_length,
            "status": status,
            "alignment_count": count,
            "aligned_query_bases": aligned_bases,
            "aligned_fraction": aligned_bases / query_length if query_length else 0.0,
            "paf_path": final_paf_path,
            "stat_path": (
                os.path.splitext(final_paf_path)[0] + ".win.stat.gz"
                if count
                else None
            ),
        }
        all_rows.append(row)
        if count:
            mapped_records.append(row)

    if not mapped_records:
        raise RuntimeError("Full assembly PAF contains no primary mapped FASTA contigs")

    max_workers = max(1, min(int(thread or 1), len(mapped_records)))
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(_run_pandepth, pandepth_path, row["paf_path"]): row
            for row in mapped_records
        }
        for future in as_completed(futures):
            future.result()

    with open(os.path.join(prefix, FULL_ASSEMBLY_PATHS_PKL), "wb") as handle:
        pickle.dump(mapped_records, handle)
    with open(os.path.join(prefix, FULL_ASSEMBLY_PATHS_TSV), "wt", newline="") as handle:
        fieldnames = [
            "fasta_index",
            "matrix_index",
            "query_name",
            "query_length",
            "status",
            "alignment_count",
            "aligned_query_bases",
            "aligned_fraction",
            "paf_path",
            "stat_path",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)

    paf_ans_list = [(row["paf_path"], [row["matrix_index"]]) for row in mapped_records]
    key_list = [row["matrix_index"] for row in mapped_records]
    int2key = {
        row["matrix_index"]: ("full_assembly", row["query_name"])
        for row in mapped_records
    }
    dep_list = [0] * len(mapped_records)
    with open(os.path.join(prefix, CONTIG_PATH_DATA_PKL), "wb") as handle:
        pickle.dump((paf_ans_list, key_list, int2key, dep_list), handle)

    summary = {
        "assembly_path": assembly_path,
        "paf_path": paf_path,
        "fasta_contigs": len(contigs),
        "mapped_contigs": len(mapped_records),
        "unmapped_contigs": len(contigs) - len(mapped_records),
        "depth_window": DEPTH_WINDOW,
    }
    with open(os.path.join(prefix, FULL_ASSEMBLY_SUMMARY_JSON), "wt") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
    logging.info(
        "Full-assembly stage 21: %d/%d contigs mapped; %d excluded as unmapped",
        len(mapped_records),
        len(contigs),
        len(contigs) - len(mapped_records),
    )
    return mapped_records


def _read_depth_frame(path: str) -> pd.DataFrame:
    frame = pd.read_csv(
        path,
        compression="gzip" if str(path).endswith(".gz") else None,
        comment="#",
        sep="\t",
        names=DEPTH_COLUMNS,
    )
    frame = frame.query('chr != "chrM"').copy()
    if frame.empty:
        raise ValueError(f"Depth file has no non-chrM rows: {path}")
    return frame


def _read_bed(path: str) -> dict[str, list[tuple[int, int]]]:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    if not path or not os.path.isfile(path):
        return intervals
    with open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            values = line.rstrip("\r\n").split("\t")
            if len(values) < 3:
                raise ValueError(f"{path}:{line_number}: malformed BED row")
            intervals[values[0]].append((int(values[1]), int(values[2])))
    return intervals


def _depth_layout(main_depth_path: str, censat_bed_path: str) -> tuple[pd.DataFrame, dict]:
    frame = _read_depth_frame(main_depth_path)
    bed = _read_bed(censat_bed_path)
    median_depth = float(np.median(frame["meandepth"].to_numpy(dtype=float)))
    clean_indices = []
    failed_indices = []
    for index, row in enumerate(frame.itertuples(index=False)):
        failed = float(row.meandepth) > ABS_MAX_COVERAGE_RATIO * median_depth
        if not failed:
            failed = any(
                bed_start <= int(row.st) and int(row.nd) <= bed_end
                for bed_start, bed_end in bed.get(row.chr, ())
            )
        (failed_indices if failed else clean_indices).append(index)
    if not clean_indices:
        raise RuntimeError("Full-assembly matrix has no clean depth windows")
    ordered_indices = clean_indices + failed_indices
    coords = [
        (str(frame.iloc[index]["chr"]), int(frame.iloc[index]["st"]), int(frame.iloc[index]["nd"]))
        for index in ordered_indices
    ]
    layout = {
        "clean_indices": clean_indices,
        "failed_indices": failed_indices,
        "ordered_indices": ordered_indices,
        "coords": coords,
        "clean_count": len(clean_indices),
        "failed_count": len(failed_indices),
        "median_depth": median_depth,
    }
    return frame, layout


def _depth_vector(stat_path: str, coords: Sequence[tuple[str, int, int]]) -> np.ndarray:
    frame = _read_depth_frame(stat_path)
    by_coord = {
        (str(row.chr), int(row.st)): float(row.meandepth)
        for row in frame.itertuples(index=False)
    }
    return np.asarray([by_coord.get((chrom, st), 0.0) for chrom, st, _ in coords], dtype=np.float32)


def _noise_amplitude(frame: pd.DataFrame, clean_indices: Sequence[int]) -> float:
    clean = frame.iloc[list(clean_indices)]
    differences = []
    for _, chrom_frame in clean.groupby("chr", sort=False):
        values = chrom_frame.sort_values("st")["meandepth"].to_numpy(dtype=float)
        if len(values) > 1:
            differences.extend(np.diff(values).tolist())
    if not differences:
        return 0.0
    return float(np.std(np.asarray(differences)) / math.sqrt(2.0))


def build_full_assembly_matrix(
    prefix: str,
    main_depth_path: str,
    censat_bed_path: str,
) -> dict:
    """Build a feature-major matrix with exactly one feature per mapped contig."""

    prefix = os.path.abspath(prefix)
    paths_path = os.path.join(prefix, FULL_ASSEMBLY_PATHS_PKL)
    if not os.path.isfile(paths_path):
        raise FileNotFoundError(
            f"Missing full-assembly stage-21 metadata: {paths_path}. Rerun from stage 21."
        )
    with open(paths_path, "rb") as handle:
        path_records = pickle.load(handle)
    if not path_records:
        raise RuntimeError("Full-assembly stage 21 produced no mapped contigs")

    frame, layout = _depth_layout(main_depth_path, censat_bed_path)
    coords = layout["coords"]
    vectors = []
    for record in path_records:
        stat_path = record["stat_path"]
        if not os.path.isfile(stat_path):
            raise FileNotFoundError(
                f"Missing PanDepth output for {record['query_name']}: {stat_path}. Rerun from stage 21."
            )
        vectors.append(_depth_vector(stat_path, coords))
    A = np.vstack(vectors).astype(np.float32, copy=False)

    ordered_frame = frame.iloc[layout["ordered_indices"]]
    B = ordered_frame["meandepth"].to_numpy(dtype=np.float32)
    clean_count = int(layout["clean_count"])
    A_clean = A[:, :clean_count]
    A_fail = A[:, clean_count:]
    B_clean = B[:clean_count]
    B_fail = B[clean_count:]

    matrix_path = os.path.join(prefix, "matrix.h5")
    with h5py.File(matrix_path, "w") as handle:
        handle.attrs["matrix_contract"] = "full_assembly_depth_v1"
        handle.create_dataset("A", data=A_clean)
        handle.create_dataset("A_fail", data=A_fail)
        handle.create_dataset("B", data=B_clean)
        handle.create_dataset("B_fail", data=B_fail)
        handle.create_dataset("B_depth_start", data=0)
        handle.create_dataset("B_depth_end", data=clean_count)

    np.save(os.path.join(prefix, "B.npy"), B)
    with open(os.path.join(prefix, FULL_ASSEMBLY_DEPTH_LAYOUT_PKL), "wb") as handle:
        pickle.dump(layout, handle)

    path_nclose_usage = [Counter() for _ in path_records]
    save_event_catalog(prefix, [])
    save_path_usage(prefix, path_nclose_usage)
    filter_status = initialise_filter_status([], path_nclose_usage, range(len(path_records)))
    save_filter_status(prefix, filter_status)

    empty_path_sets = {index: set() for index in range(len(path_records))}
    amplitude = _noise_amplitude(frame, layout["clean_indices"])
    bed_data = _read_bed(censat_bed_path)
    matrix_meta = {
        "B_depth_start": 0,
        "B_depth_end": clean_count,
        "full_assembly": True,
    }
    with open(os.path.join(prefix, "23_input.pkl"), "wb") as handle:
        pickle.dump(
            (
                [(chrom, st) for chrom, st, _ in coords[:clean_count]],
                empty_path_sets,
                amplitude,
                dict(bed_data),
                matrix_meta,
            ),
            handle,
        )
    for filename, value in (
        ("path_data.pkl", {}),
        ("cen_fragment_data.pkl", {}),
        ("ecdna_circuit_data.pkl", ([], {})),
        ("conjoined_type4_ins_del.pkl", ([], [])),
        ("indel_exclude_idx_set.pkl", set()),
    ):
        with open(os.path.join(prefix, filename), "wb") as handle:
            pickle.dump(value, handle)
    save_nclose_nodes(prefix, {})
    with open(os.path.join(prefix, "tot_loc_list.pkl"), "wb") as handle:
        pickle.dump([record["paf_path"] for record in path_records], handle)

    logging.info(
        "Full-assembly stage 22: matrix=%s, contigs=%d, clean_windows=%d, failed_windows=%d",
        A.shape,
        len(path_records),
        clean_count,
        len(B_fail),
    )
    return {
        "features": len(path_records),
        "clean_windows": clean_count,
        "failed_windows": len(B_fail),
        "matrix_path": matrix_path,
    }


def solve_full_assembly_nnls(prefix: str, main_depth_path: str) -> dict:
    prefix = os.path.abspath(prefix)
    matrix_path = os.path.join(prefix, "matrix.h5")
    if not os.path.isfile(matrix_path):
        raise FileNotFoundError(f"Missing full-assembly matrix: {matrix_path}. Rerun from stage 22.")
    with h5py.File(matrix_path, "r") as handle:
        A = handle["A"][:]
        A_fail = handle["A_fail"][:]
        B = handle["B"][:]
        B_fail = handle["B_fail"][:]
    if A.shape[0] == 0:
        raise RuntimeError("Full-assembly matrix has no contig features")
    weights, _ = nnls(A.T, B)
    weights = weights.astype(np.float32, copy=False)
    predict_clean = A.T.dot(weights)
    predict_fail = A_fail.T.dot(weights) if A_fail.shape[1] else np.asarray([], dtype=np.float32)
    predict = np.concatenate((predict_clean, predict_fail)).astype(np.float32, copy=False)
    np.save(os.path.join(prefix, "weight.npy"), weights)
    np.save(os.path.join(prefix, "predict_B.npy"), predict)
    with open(os.path.join(prefix, "A_idx_list.pkl"), "wb") as handle:
        pickle.dump(list(range(len(weights))), handle)

    with open(os.path.join(prefix, FULL_ASSEMBLY_PATHS_PKL), "rb") as handle:
        path_records = pickle.load(handle)
    depth_frame = _read_depth_frame(main_depth_path)
    median_depth = float(np.median(depth_frame["meandepth"].to_numpy(dtype=float)))
    haploid_depth = median_depth / 2.0
    with open(os.path.join(prefix, FULL_ASSEMBLY_WEIGHTS_TSV), "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["matrix_index", "query_name", "query_length", "raw_weight", "copy_number_N", "paf_path"]
        )
        for record, weight in zip(path_records, weights, strict=True):
            writer.writerow(
                [
                    record["matrix_index"],
                    record["query_name"],
                    record["query_length"],
                    f"{float(weight):.8g}",
                    f"{float(weight / haploid_depth) if haploid_depth else 0.0:.8g}",
                    record["paf_path"],
                ]
            )

    norm = float(np.linalg.norm(B))
    error = float(np.linalg.norm(predict_clean - B))
    relative_error = error / norm if norm else 0.0
    logging.info("Full-assembly stage 23 error: %.4f (relative %.6f)", error, relative_error)
    return {
        "weights": weights,
        "error": error,
        "relative_error": relative_error,
        "failed_target_windows": len(B_fail),
    }


def _load_primary_paf(path: str) -> list[PafRecord]:
    records = []
    with open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            record = parse_paf_line(line, path, line_number)
            if not record.is_secondary:
                records.append(record)
    return sorted(
        records,
        key=lambda record: (
            record.query_start,
            record.query_end,
            record.target_name,
            record.target_start,
        ),
    )


def paf_records_to_pieces(records: Sequence[PafRecord]) -> list[tuple[str, str, int]]:
    """Convert query-ordered alignments into merged (chrom, strand, query-bp) pieces."""

    pieces: list[list] = []
    consumed_query_end = 0
    for record in records:
        usable_start = max(record.query_start, consumed_query_end)
        usable_length = record.query_end - usable_start
        if usable_length <= 0:
            continue
        consumed_query_end = max(consumed_query_end, record.query_end)
        state = (record.target_name, record.strand)
        if pieces and tuple(pieces[-1][:2]) == state:
            pieces[-1][2] += usable_length
        else:
            pieces.append([record.target_name, record.strand, usable_length])
    return [(str(chrom), str(strand), int(length)) for chrom, strand, length in pieces]


def filter_full_assembly_display_pieces(
    pieces: Sequence[tuple[str, str, int]],
    minimum_length: int = FULL_ASSEMBLY_PLOT_MIN_PIECE,
) -> list[tuple[str, str, int]]:
    """Remove sub-ISCN noise pieces and merge the newly adjacent states.

    Minimap2 can emit hundreds of short, overlapping primary alignments around
    repetitive sequence even with ``--secondary=no``.  The existing Virtual
    SKY/ISCN code treats sections below 1 Mb as non-karyotypic, so applying the
    same cutoff before rendering prevents those alignments from becoming a
    wall of false junction labels.  Tiny toy contigs retain their longest piece
    so the renderer remains useful in tests and small references.
    """

    retained = [piece for piece in pieces if int(piece[2]) >= int(minimum_length)]
    if not retained and pieces:
        retained = [max(pieces, key=lambda piece: int(piece[2]))]

    merged: list[list] = []
    for chrom, strand, length in retained:
        state = (str(chrom), str(strand))
        if merged and tuple(merged[-1][:2]) == state:
            merged[-1][2] += int(length)
        else:
            merged.append([state[0], state[1], int(length)])
    return [tuple(piece) for piece in merged]


def pieces_to_karyotype(pieces: Sequence[tuple[str, str, int]]) -> str:
    visible = [piece for piece in pieces if piece[2] >= 100_000]
    if not visible:
        return "unaligned"
    if len(visible) == 1:
        return visible[0][0].removeprefix("chr")
    tokens = []
    for (left_chrom, left_strand, _), (right_chrom, right_strand, _) in zip(
        visible, visible[1:]
    ):
        if left_chrom != right_chrom:
            tokens.append(
                f"t({left_chrom.removeprefix('chr')};{right_chrom.removeprefix('chr')})"
            )
        elif left_strand != right_strand:
            tokens.append(f"inv({right_chrom.removeprefix('chr')})")
    return "".join(tokens) or visible[0][0].removeprefix("chr")


def prepare_full_assembly_virtual_sky(
    prefix: str,
    main_depth_path: str,
    minimum_depth_n: float = 0.2,
) -> dict:
    """Prepare full-assembly paths for the canonical stage-30 renderer."""

    prefix = os.path.abspath(prefix)
    with open(os.path.join(prefix, FULL_ASSEMBLY_PATHS_PKL), "rb") as handle:
        path_records = pickle.load(handle)
    weights = np.load(os.path.join(prefix, "weight.npy"))
    if len(weights) != len(path_records):
        raise ValueError(
            f"Full-assembly weight/path mismatch: {len(weights)} != {len(path_records)}"
        )
    frame = _read_depth_frame(main_depth_path)
    median_depth = float(np.median(frame["meandepth"].to_numpy(dtype=float)))
    haploid_depth = median_depth / 2.0
    threshold = float(minimum_depth_n) * haploid_depth

    display_rows = []
    for record, weight in zip(path_records, weights, strict=True):
        if float(weight) <= threshold:
            continue
        raw_pieces = paf_records_to_pieces(_load_primary_paf(record["paf_path"]))
        pieces = filter_full_assembly_display_pieces(raw_pieces)
        if not pieces:
            continue
        display_rows.append(
            {
                "record": record,
                "weight": float(weight),
                "copy_number": float(weight / haploid_depth) if haploid_depth else 0.0,
                "pieces": pieces,
                "raw_piece_count": len(raw_pieces),
                "karyotype": pieces_to_karyotype(pieces),
            }
        )

    logging.info(
        "Full-assembly stage 30: prepared %d contigs for the shared renderer",
        len(display_rows),
    )
    return {
        "display_rows": display_rows,
        "displayed_contigs": len(display_rows),
        "median_depth": median_depth,
    }


def _read_chromosome_lengths(path: str) -> dict[str, int]:
    lengths = {}
    with open(path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            values = line.rstrip("\r\n").split("\t")
            if len(values) < 2:
                raise ValueError(f"{path}:{line_number}: malformed FASTA index row")
            lengths[values[0]] = int(values[1])
    # Stage 31 removes chrM from the depth frame before plotting; keep the
    # canonical renderer's sector list in sync with that behavior.
    lengths.pop("chrM", None)
    if not lengths:
        raise ValueError(f"Reference FASTA index is empty: {path}")
    return lengths


def _full_assembly_karyotype(
    pieces: Sequence[tuple[str, str, int]],
    chromosome_lengths: dict[str, int],
) -> str | None:
    """Apply the normal stage-30 ISCN rules to assembly-derived pieces."""

    visible = [piece for piece in pieces if int(piece[2]) >= FULL_ASSEMBLY_PLOT_MIN_PIECE]
    if not visible:
        return None
    junctions = []
    for (left_chrom, left_strand, _), (right_chrom, right_strand, _) in zip(
        visible, visible[1:]
    ):
        if left_chrom != right_chrom:
            junctions.append(
                f"t({left_chrom.removeprefix('chr')};{right_chrom.removeprefix('chr')})"
            )
        elif left_strand != right_strand:
            junctions.append(f"inv({right_chrom.removeprefix('chr')})")
    if junctions:
        return "".join(junctions)

    chrom = visible[0][0]
    total_length = sum(length for piece_chrom, _strand, length in visible if piece_chrom == chrom)
    reference_length = chromosome_lengths.get(chrom)
    if reference_length and total_length / reference_length < 0.9:
        return f"del({chrom.removeprefix('chr')})"
    return chrom.removeprefix("chr")


def draw_full_assembly_virtual_sky(
    prefix: str,
    main_depth_path: str,
    reference_fai_path: str,
    cell_line: str,
    minimum_depth_n: float = 0.2,
) -> dict:
    """Render assembly paths directly with the canonical stage-30 function."""

    prefix = os.path.abspath(prefix)
    chromosome_lengths = _read_chromosome_lengths(reference_fai_path)
    prepared = prepare_full_assembly_virtual_sky(
        prefix,
        main_depth_path,
        minimum_depth_n=minimum_depth_n,
    )

    path_pieces = {}
    path_depth_n = {}
    path_karyotype = {}
    for row in prepared["display_rows"]:
        path = row["record"]["query_name"]
        pieces = [
            ((chrom, strand), int(length))
            for chrom, strand, length in row["pieces"]
        ]
        if not pieces:
            continue
        path_pieces[path] = pieces
        path_depth_n[path] = float(row["copy_number"])
        path_karyotype[path] = _full_assembly_karyotype(
            row["pieces"], chromosome_lengths
        )

    maxh = max(chromosome_lengths.values())
    for pieces in path_pieces.values():
        maxh = max(maxh, sum(length for _state, length in pieces))

    grouped_norm_data = defaultdict(list)
    for path, pieces in path_pieces.items():
        normalized = [
            (state, length / maxh * 100)
            for state, length in pieces
        ]
        chromosome_weight = Counter()
        for state, length in normalized:
            chromosome_weight[state] += length
        primary_chromosome = sorted(
            chromosome_weight.items(), key=lambda value: -value[1]
        )[0][0][0]
        grouped_norm_data[primary_chromosome].append((path, normalized))

    result = render_karyotype_diagram(
        output_prefix=prefix,
        cell_line=cell_line,
        chromosome_lengths=chromosome_lengths,
        grouped_norm_data=grouped_norm_data,
        display_indel=defaultdict(list),
        virtual_inv_display=[],
        maxh=maxh,
        path_depth_n=path_depth_n,
        path_karyotype=path_karyotype,
        path_event_labels={},
    )
    result["displayed_contigs"] = len(path_pieces)
    logging.info(
        "Full-assembly stage 30: rendered %d contigs with the canonical Virtual SKY figure",
        len(path_pieces),
    )
    return result


def draw_full_assembly_coverage(
    prefix: str,
    main_depth_path: str,
    reference_fai_path: str,
    reference_cytobands_path: str,
) -> dict:
    """Render full-assembly coverage with the canonical stage-31 circos code."""

    # Import lazily: seaborn configures matplotlib's global savefig bbox on
    # import.  Keeping it out of stage 30 preserves the legacy Virtual SKY
    # canvas dimensions exactly.
    from circos_plotting import (
        rebin_depth_dataframe,
        render_total_coverage_circos,
        smooth_depth_for_coords,
    )

    prefix = os.path.abspath(prefix)
    with open(os.path.join(prefix, FULL_ASSEMBLY_DEPTH_LAYOUT_PKL), "rb") as handle:
        layout = pickle.load(handle)
    observed = np.load(os.path.join(prefix, "B.npy"))
    predicted = np.maximum(np.load(os.path.join(prefix, "predict_B.npy")), 0)
    if len(observed) != len(layout["coords"]) or len(predicted) != len(observed):
        raise ValueError("Full-assembly coverage vector/layout length mismatch")

    clean_count = int(layout["clean_count"])
    clean_observed = observed[:clean_count]
    clean_predicted = predicted[:clean_count]
    error = float(np.linalg.norm(clean_predicted - clean_observed))
    norm = float(np.linalg.norm(clean_observed))
    relative_error = error / norm if norm else 0.0
    frame = _read_depth_frame(main_depth_path)
    coords = [(str(chrom), int(start)) for chrom, start, _end in layout["coords"]]
    smoothed_observed = smooth_depth_for_coords(frame, coords, 10)
    amplitude = _noise_amplitude(frame, layout["clean_indices"])
    miss = np.abs(predicted - smoothed_observed) > amplitude * 1.959963984540054
    labels = np.ones(len(predicted), dtype=float)
    labels[miss] = 3
    labels[clean_count:] = 5

    figure_paths = render_total_coverage_circos(
        output_prefix=prefix,
        fig_prefix="",
        chromosome_lengths=_read_chromosome_lengths(reference_fai_path),
        cytobands_path=reference_cytobands_path,
        observed_depth=rebin_depth_dataframe(frame, 2),
        prediction_coords=coords,
        predicted_depth=predicted,
        prediction_labels=labels,
        mean_depth=float(layout["median_depth"]),
        fragment_depth_per_chrom={},
        centromere_fragment_meta={},
        breakend_links=(),
        telomere_arrows=(),
    )
    with open(os.path.join(prefix, "cn_data.pkl"), "wb") as handle:
        pickle.dump(([], [], [], []), handle)

    weights = np.load(os.path.join(prefix, "weight.npy"))
    with open(os.path.join(prefix, "report.txt"), "wt") as handle:
        handle.write("mode\tfull_assembly\n")
        handle.write(f"matrix_columns\t{len(weights)}\n")
        handle.write(f"clean_windows\t{clean_count}\n")
        handle.write(f"failed_windows\t{layout['failed_count']}\n")
        handle.write(f"median_depth\t{layout['median_depth']:.8g}\n")
        handle.write(f"error\t{error:.8g}\n")
        handle.write(f"relative_error\t{relative_error:.8g}\n")
    logging.info("Full-assembly stage 31 relative error: %.6f", relative_error)
    return {
        "error": error,
        "relative_error": relative_error,
        **figure_paths,
    }


def _validate_restart_metadata(
    prefix: str,
    assembly_path: str,
    aligned_paf_path: str,
    reference: str,
) -> None:
    config = load_pipeline_input(prefix)
    if not pipeline_input_is_full_assembly(config):
        raise FileNotFoundError(
            f"Full-assembly restart metadata is missing from {prefix}. "
            "Restart from the beginning of full_assembly_pipeline.py."
        )
    expected = {
        "full_assembly_path": os.path.abspath(assembly_path),
        "full_assembly_paf_path": os.path.abspath(aligned_paf_path),
        "reference": reference,
    }
    mismatches = {
        key: (config.get(key), value)
        for key, value in expected.items()
        if config.get(key) != value
    }
    if mismatches:
        raise ValueError(
            "Full-assembly restart metadata does not match the requested inputs: "
            f"{mismatches}"
        )


def run_full_assembly_pipeline(
    *,
    aligned_paf_path: str,
    assembly_path: str,
    main_depth_path: str,
    censat_bed_path: str,
    telomere_bed_path: str,
    reference_fai_path: str,
    reference_cytobands_path: str,
    prefix: str,
    cell_line: str,
    pandepth_path: str,
    thread: int,
    reference: str,
    start_at: int = 21,
) -> dict:
    """Run the complete full-assembly workflow behind one dedicated entrypoint."""

    aligned_paf_path = os.path.abspath(aligned_paf_path)
    assembly_path = os.path.abspath(assembly_path)
    main_depth_path = os.path.abspath(main_depth_path)
    prefix = os.path.abspath(prefix)
    if not aligned_paf_path.endswith(".aln.paf"):
        raise ValueError(
            "full_assembly_pipeline.py requires an alignasm output ending in .aln.paf"
        )
    for path, label in (
        (aligned_paf_path, "alignasm PAF"),
        (assembly_path, "full assembly FASTA"),
        (main_depth_path, "sample depth statistics"),
        (reference_fai_path, "reference FAI"),
        (reference_cytobands_path, "reference cytobands"),
    ):
        if not os.path.isfile(path):
            raise FileNotFoundError(f"Missing {label}: {path}")
    if start_at not in {21, 22, 23, 30, 31}:
        raise ValueError("--start-at must be one of 21, 22, 23, 30, or 31")

    os.makedirs(prefix, exist_ok=True)
    if start_at == 21:
        save_pipeline_input(
            prefix,
            full_assembly_input=True,
            full_assembly_path=assembly_path,
            full_assembly_paf_path=aligned_paf_path,
            reference=reference,
        )
    else:
        _validate_restart_metadata(
            prefix, assembly_path, aligned_paf_path, reference
        )

    result = {}
    if start_at <= 21:
        result["depth"] = run_full_assembly_stage21(
            assembly_path=assembly_path,
            paf_path=aligned_paf_path,
            prefix=prefix,
            pandepth_path=pandepth_path,
            thread=thread,
        )
    if start_at <= 22:
        result["matrix"] = build_full_assembly_matrix(
            prefix=prefix,
            main_depth_path=main_depth_path,
            censat_bed_path=censat_bed_path,
        )
    if start_at <= 23:
        result["nnls"] = solve_full_assembly_nnls(prefix, main_depth_path)
    if start_at <= 30:
        result["virtual_sky"] = draw_full_assembly_virtual_sky(
            prefix=prefix,
            main_depth_path=main_depth_path,
            reference_fai_path=reference_fai_path,
            cell_line=cell_line,
        )
    if start_at <= 31:
        result["coverage"] = draw_full_assembly_coverage(
            prefix,
            main_depth_path,
            reference_fai_path,
            reference_cytobands_path,
        )
    return result


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run SKYPE using only contigs from a complete assembly"
    )
    parser.add_argument("aligned_paf", help="alignasm-produced *.aln.paf")
    parser.add_argument("assembly_fasta", help="Complete assembly FASTA")
    parser.add_argument("main_depth", help="Sample PanDepth *.win.stat.gz")
    parser.add_argument("censat_bed", help="Reference censat BED")
    parser.add_argument("telomere_bed", help="Reference telomere BED")
    parser.add_argument("reference_fai", help="Reference FASTA index")
    parser.add_argument("reference_cytobands", help="Reference cytoband BED")
    parser.add_argument("prefix", help="Full-assembly output directory")
    parser.add_argument("cell_line", help="Sample/cell-line label")
    parser.add_argument("--pandepth-loc", default="pandepth")
    parser.add_argument("-t", "--thread", type=int, default=1)
    parser.add_argument("--reference", choices=("hs1", "hg38"), default="hs1")
    parser.add_argument(
        "--start-at",
        type=int,
        choices=(21, 22, 23, 30, 31),
        default=21,
    )
    return parser


def main() -> None:
    args = get_parser().parse_args()
    run_full_assembly_pipeline(
        aligned_paf_path=args.aligned_paf,
        assembly_path=args.assembly_fasta,
        main_depth_path=args.main_depth,
        censat_bed_path=args.censat_bed,
        telomere_bed_path=args.telomere_bed,
        reference_fai_path=args.reference_fai,
        reference_cytobands_path=args.reference_cytobands,
        prefix=args.prefix,
        cell_line=args.cell_line,
        pandepth_path=args.pandepth_loc,
        thread=args.thread,
        reference=args.reference,
        start_at=args.start_at,
    )


if __name__ == "__main__":
    main()
