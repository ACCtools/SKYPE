"""Karyotype-only censat--censat NClose rescue helpers.

The stage-02 selector intentionally uses only the aligned and original unitig
PAFs.  Selected records are persisted separately from the existing NClose
pickle schemas so older runs remain a no-op in stages 21 and 23.
"""

from __future__ import annotations

import csv
import os
import pickle
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np


CENSAT_PAIR_RESCUE_PKL = "censat_pair_rescue.pkl"
CENSAT_PAIR_RESCUE_VERSION = 1
DEFAULT_MIN_CENSAT_LENGTH = 1_000_000


@dataclass(frozen=True)
class Alignment:
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
    source_order: int

    @property
    def state(self) -> tuple[str, str]:
        return self.target_name, self.strand


@dataclass(frozen=True)
class RescueCandidate:
    query_name: str
    query_length: int
    left: Alignment
    right: Alignment
    canonical_states: tuple[tuple[str, str], tuple[str, str]]

    @property
    def chroms(self) -> tuple[str, str]:
        return tuple(state[0] for state in self.canonical_states)


@dataclass(frozen=True)
class RescueCFSwapPlan:
    event_keys: frozenset[tuple[int, int]]
    rescue_columns: tuple[int, ...]
    active_event_chroms: frozenset[str]
    restored_cf_chroms: frozenset[str]
    disabled_cf_chroms: frozenset[str]
    disabled_cf_columns: tuple[int, ...]


def _as_int(value: str, field: str, path: Path, line_number: int) -> int:
    try:
        return int(value)
    except ValueError as exc:
        raise ValueError(
            f"{path}:{line_number}: invalid integer in {field}: {value!r}"
        ) from exc


def _make_alignment(
    values: Sequence[str],
    path: Path,
    line_number: int,
    source_order: int,
) -> Alignment:
    if len(values) < 12:
        raise ValueError(
            f"{path}:{line_number}: expected at least 12 alignment columns, "
            f"found {len(values)}"
        )

    alignment = Alignment(
        query_name=values[0],
        query_length=_as_int(values[1], "query length", path, line_number),
        query_start=_as_int(values[2], "query start", path, line_number),
        query_end=_as_int(values[3], "query end", path, line_number),
        strand=values[4],
        target_name=values[5],
        target_length=_as_int(values[6], "target length", path, line_number),
        target_start=_as_int(values[7], "target start", path, line_number),
        target_end=_as_int(values[8], "target end", path, line_number),
        matches=_as_int(values[9], "matching bases", path, line_number),
        block_length=_as_int(values[10], "alignment block length", path, line_number),
        mapq=_as_int(values[11], "MAPQ", path, line_number),
        source_order=int(source_order),
    )
    if alignment.strand not in {"+", "-"}:
        raise ValueError(
            f"{path}:{line_number}: invalid strand {alignment.strand!r}"
        )
    if not (
        0
        <= alignment.query_start
        <= alignment.query_end
        <= alignment.query_length
    ):
        raise ValueError(
            f"{path}:{line_number}: invalid query interval "
            f"[{alignment.query_start}, {alignment.query_end})"
        )
    if not (
        0
        <= alignment.target_start
        <= alignment.target_end
        <= alignment.target_length
    ):
        raise ValueError(
            f"{path}:{line_number}: invalid target interval "
            f"[{alignment.target_start}, {alignment.target_end})"
        )
    return alignment


def read_alignments(path_value: str | os.PathLike) -> list[Alignment]:
    """Read standard PAF or the project's twelve-column PAF CSV.

    For PAF input, ``source_order`` is the physical zero-based line index.  It
    can therefore be used directly in ``CTG_GLOBALIDX`` by stage 21.
    """

    path = Path(path_value)
    rows = []
    with path.open("rt", newline="") as handle:
        first_line = handle.readline()
        handle.seek(0)
        is_csv = first_line.startswith("qry_name,")

        if is_csv:
            fields = [
                "qry_name",
                "qry_len",
                "qry_str",
                "qry_end",
                "aln_dir",
                "ref_name",
                "ref_len",
                "ref_str",
                "ref_end",
                "mat_len",
                "aln_len",
                "mapq",
            ]
            reader = csv.DictReader(handle)
            if reader.fieldnames != fields:
                raise ValueError(
                    f"{path}: unexpected CSV header {reader.fieldnames!r}; "
                    f"expected {fields!r}"
                )
            for source_order, row in enumerate(reader):
                rows.append(
                    _make_alignment(
                        [row[field] for field in fields],
                        path,
                        reader.line_num,
                        source_order,
                    )
                )
            return rows

        for line_index, line in enumerate(handle):
            if not line.strip() or line.startswith("#"):
                continue
            rows.append(
                _make_alignment(
                    line.rstrip("\n").split("\t"),
                    path,
                    line_index + 1,
                    line_index,
                )
            )
    return rows


def group_alignments(
    alignments_or_path: Iterable[Alignment] | str | os.PathLike,
) -> dict[str, list[Alignment]]:
    if isinstance(alignments_or_path, (str, os.PathLike)):
        alignments = read_alignments(alignments_or_path)
    else:
        alignments = list(alignments_or_path)

    grouped = defaultdict(list)
    query_lengths = {}
    for alignment in alignments:
        previous = query_lengths.setdefault(
            alignment.query_name, alignment.query_length
        )
        if previous != alignment.query_length:
            raise ValueError(
                f"Inconsistent query length for {alignment.query_name}: "
                f"{previous} versus {alignment.query_length}"
            )
        grouped[alignment.query_name].append(alignment)
    for rows in grouped.values():
        rows.sort(
            key=lambda row: (
                row.query_start,
                row.query_end,
                row.source_order,
            )
        )
    return dict(grouped)


def read_censat_bed(
    path_value: str | os.PathLike,
    min_censat_length: int = DEFAULT_MIN_CENSAT_LENGTH,
) -> dict[str, list[tuple[int, int]]]:
    path = Path(path_value)
    intervals = defaultdict(list)
    with path.open("rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(
                    f"{path}:{line_number}: expected at least three BED columns"
                )
            start = _as_int(fields[1], "BED start", path, line_number)
            end = _as_int(fields[2], "BED end", path, line_number)
            if not 0 <= start <= end:
                raise ValueError(
                    f"{path}:{line_number}: invalid BED interval [{start}, {end})"
                )
            if end - start > int(min_censat_length):
                intervals[fields[0]].append((start, end))
    for chrom_intervals in intervals.values():
        chrom_intervals.sort()
    return dict(intervals)


def alignment_overlaps_censat(
    alignment: Alignment,
    censat_intervals: Mapping[str, Sequence[tuple[int, int]]],
) -> bool:
    return any(
        max(alignment.target_start, start)
        < min(alignment.target_end, end)
        for start, end in censat_intervals.get(alignment.target_name, ())
    )


def merge_intervals(
    intervals: Iterable[tuple[int, int]],
) -> list[tuple[int, int]]:
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = merged[-1][0], max(merged[-1][1], end)
        else:
            merged.append((start, end))
    return merged


def chromosome_sort_key(chromosome: str) -> tuple[int, int | str]:
    match = re.fullmatch(r"chr(\d+)", chromosome)
    if match:
        return 0, int(match.group(1))
    special = {"chrX": 23, "chrY": 24, "chrM": 25}
    if chromosome in special:
        return 0, special[chromosome]
    return 1, chromosome


def state_sort_key(
    state: tuple[str, str],
) -> tuple[tuple[int, int | str], str]:
    return chromosome_sort_key(state[0]), state[1]


def flip_state(state: tuple[str, str]) -> tuple[str, str]:
    chrom, strand = state
    if strand not in {"+", "-"}:
        raise ValueError(f"Invalid state strand: {state!r}")
    return chrom, "-" if strand == "+" else "+"


def reverse_complement_pair(
    pair: tuple[tuple[str, str], tuple[str, str]],
) -> tuple[tuple[str, str], tuple[str, str]]:
    return flip_state(pair[1]), flip_state(pair[0])


def _pair_sort_key(pair):
    return state_sort_key(pair[0]), state_sort_key(pair[1])


def canonicalize_pair(
    pair: tuple[tuple[str, str], tuple[str, str]],
) -> tuple[tuple[str, str], tuple[str, str]]:
    reverse_complement = reverse_complement_pair(pair)
    return min(pair, reverse_complement, key=_pair_sort_key)


def cen_fragment_target_direction(
    cen_fragment_meta: Mapping,
    chromosome: str,
) -> str | None:
    info = cen_fragment_meta.get(chromosome)
    if not isinstance(info, Mapping) or "dir" not in info:
        return None
    return "-" if bool(info["dir"]) else "+"


def terminal_directions_match_cen_fragment(
    left: Alignment,
    right: Alignment,
    cen_fragment_meta: Mapping,
) -> bool:
    left_target = cen_fragment_target_direction(
        cen_fragment_meta, left.target_name
    )
    right_target = cen_fragment_target_direction(
        cen_fragment_meta, right.target_name
    )
    normalized_right = flip_state(right.state)[1]
    return (
        left_target is not None
        and right_target is not None
        and left.strand == left_target
        and normalized_right == right_target
    )


def _strict_original_state_order(
    original_rows: Sequence[Alignment],
    left_state: tuple[str, str],
    right_state: tuple[str, str],
) -> bool:
    expected_states = {left_state, right_state}
    if {row.state for row in original_rows} != expected_states:
        return False

    rows_by_state = defaultdict(list)
    for row in original_rows:
        rows_by_state[row.state].append(row)
    left_union = merge_intervals(
        (row.query_start, row.query_end)
        for row in rows_by_state[left_state]
    )
    right_union = merge_intervals(
        (row.query_start, row.query_end)
        for row in rows_by_state[right_state]
    )
    left_span = left_union[0][0], left_union[-1][1]
    right_span = right_union[0][0], right_union[-1][1]
    return (
        left_span[0] < right_span[0]
        and left_span[1] < right_span[1]
    )


def select_rescue_candidates_from_groups(
    aligned_by_query: Mapping[str, Sequence[Alignment]],
    original_by_query: Mapping[str, Sequence[Alignment]],
    censat_intervals: Mapping[str, Sequence[tuple[int, int]]],
    cen_fragment_meta: Mapping,
) -> list[RescueCandidate]:
    """Select one longest query per reverse-complement-equivalent state pair."""

    grouped_candidates = defaultdict(list)
    for query_name, aligned_rows_value in aligned_by_query.items():
        aligned_rows = sorted(
            aligned_rows_value,
            key=lambda row: (
                row.query_start,
                row.query_end,
                row.source_order,
            ),
        )
        if not aligned_rows:
            continue
        left = aligned_rows[0]
        right = aligned_rows[-1]
        if left.target_name == right.target_name:
            continue
        if not (
            alignment_overlaps_censat(left, censat_intervals)
            and alignment_overlaps_censat(right, censat_intervals)
        ):
            continue
        original_rows = original_by_query.get(query_name, ())
        if not original_rows:
            continue
        if not _strict_original_state_order(
            original_rows, left.state, right.state
        ):
            continue
        if not terminal_directions_match_cen_fragment(
            left, right, cen_fragment_meta
        ):
            continue

        canonical_states = canonicalize_pair((left.state, right.state))
        grouped_candidates[canonical_states].append(
            RescueCandidate(
                query_name=query_name,
                query_length=left.query_length,
                left=left,
                right=right,
                canonical_states=canonical_states,
            )
        )

    selected = []
    for candidates in grouped_candidates.values():
        candidates.sort(
            key=lambda candidate: (
                -candidate.query_length,
                candidate.query_name,
            )
        )
        selected.append(candidates[0])
    selected.sort(key=lambda candidate: _pair_sort_key(candidate.canonical_states))
    return selected


def select_rescue_candidates(
    aligned_paf: str | os.PathLike,
    original_paf: str | os.PathLike,
    censat_bed: str | os.PathLike,
    cen_fragment_meta: Mapping,
    min_censat_length: int = DEFAULT_MIN_CENSAT_LENGTH,
) -> list[RescueCandidate]:
    return select_rescue_candidates_from_groups(
        group_alignments(aligned_paf),
        group_alignments(original_paf),
        read_censat_bed(censat_bed, min_censat_length),
        cen_fragment_meta,
    )


def rescue_record(
    candidate: RescueCandidate,
    event_key: Sequence[int],
) -> dict:
    event_key = tuple(sorted(map(int, event_key)))
    if len(event_key) != 2:
        raise ValueError(f"Invalid rescue event key: {event_key!r}")
    return {
        "event_key": event_key,
        "query_name": str(candidate.query_name),
        "query_length": int(candidate.query_length),
        "canonical_states": tuple(
            (str(chrom), str(strand))
            for chrom, strand in candidate.canonical_states
        ),
        "chroms": tuple(map(str, candidate.chroms)),
    }


def _normalize_record(record: Mapping) -> dict:
    event_key = tuple(sorted(map(int, record["event_key"])))
    canonical_states = tuple(
        (str(state[0]), str(state[1]))
        for state in record["canonical_states"]
    )
    chroms = tuple(map(str, record["chroms"]))
    if len(event_key) != 2:
        raise ValueError(f"Invalid rescue event key: {event_key!r}")
    if (
        len(canonical_states) != 2
        or any(strand not in {"+", "-"} for _, strand in canonical_states)
    ):
        raise ValueError(
            f"Invalid rescue canonical states: {canonical_states!r}"
        )
    if len(chroms) != 2:
        raise ValueError(f"Invalid rescue chromosome pair: {chroms!r}")
    if chroms != tuple(state[0] for state in canonical_states):
        raise ValueError(
            "Rescue chromosome pair does not match canonical states: "
            f"{chroms!r} versus {canonical_states!r}"
        )
    return {
        "event_key": event_key,
        "query_name": str(record["query_name"]),
        "query_length": int(record["query_length"]),
        "canonical_states": canonical_states,
        "chroms": chroms,
    }


def empty_rescue_artifact() -> dict:
    return {
        "version": CENSAT_PAIR_RESCUE_VERSION,
        "records": [],
    }


def save_rescue_artifact(prefix: str | os.PathLike, records: Sequence[Mapping]):
    normalized = [_normalize_record(record) for record in records]
    event_keys = [record["event_key"] for record in normalized]
    if len(event_keys) != len(set(event_keys)):
        raise ValueError("Duplicate event_key in censat-pair rescue artifact")
    artifact = {
        "version": CENSAT_PAIR_RESCUE_VERSION,
        "records": normalized,
    }
    path = os.path.join(os.fspath(prefix), CENSAT_PAIR_RESCUE_PKL)
    with open(path, "wb") as handle:
        pickle.dump(artifact, handle)
    return artifact


def load_rescue_artifact(prefix: str | os.PathLike) -> dict:
    path = os.path.join(os.fspath(prefix), CENSAT_PAIR_RESCUE_PKL)
    if not os.path.isfile(path):
        return empty_rescue_artifact()
    with open(path, "rb") as handle:
        artifact = pickle.load(handle)
    if (
        not isinstance(artifact, Mapping)
        or artifact.get("version") != CENSAT_PAIR_RESCUE_VERSION
        or not isinstance(artifact.get("records"), (list, tuple))
    ):
        raise ValueError(
            f"Unsupported censat-pair rescue artifact schema: {path}"
        )
    return {
        "version": CENSAT_PAIR_RESCUE_VERSION,
        "records": [
            _normalize_record(record)
            for record in artifact["records"]
        ],
    }


def rescue_event_keys(records: Sequence[Mapping]) -> set[tuple[int, int]]:
    return {
        tuple(sorted(map(int, record["event_key"])))
        for record in records
    }


def stable_rescue_prefix_permutation(
    contains_rescue: Sequence[bool],
) -> list[int]:
    """Return a stable true-prefix/false-suffix permutation."""

    flags = [bool(value) for value in contains_rescue]
    return (
        [index for index, flag in enumerate(flags) if flag]
        + [index for index, flag in enumerate(flags) if not flag]
    )


def rescue_carrier_columns(
    path_usage: Sequence[Mapping],
    event_keys: Iterable[tuple[int, int]],
    column_limit: int | None = None,
) -> set[int]:
    event_keys = set(event_keys)
    limit = len(path_usage) if column_limit is None else int(column_limit)
    return {
        column_index
        for column_index, usage in enumerate(path_usage[:limit])
        if any(int(usage.get(event_key, 0)) > 0 for event_key in event_keys)
    }


def validate_rescue_path_prefix(
    path_usage: Sequence[Mapping],
    records: Sequence[Mapping],
    path_column_count: int,
) -> int:
    carriers = sorted(
        rescue_carrier_columns(
            path_usage,
            rescue_event_keys(records),
            column_limit=path_column_count,
        )
    )
    expected = list(range(len(carriers)))
    if carriers != expected:
        raise ValueError(
            "Censat-pair rescue path columns are not a contiguous prefix: "
            f"found {carriers[:20]}, expected {expected[:20]}. "
            "Rerun SKYPE from stage 21."
        )
    return len(carriers)


def eligible_rescue_event_keys(
    records: Sequence[Mapping],
    retained_cent_fragment_chroms: Iterable[str],
) -> set[tuple[int, int]]:
    retained = set(map(str, retained_cent_fragment_chroms))
    return {
        tuple(record["event_key"])
        for record in records
        if set(record["chroms"]) <= retained
    }


def allowed_rescue_reactivation_columns(
    deferred_columns: Iterable[int],
    path_usage: Sequence[Mapping],
    rescue_keys: Iterable[tuple[int, int]],
    blocked_columns: Iterable[int] = (),
    blocked_events: Iterable = (),
) -> set[int]:
    """Keep deferred carriers unless an unrelated earlier filter rejected them."""

    rescue_keys = set(rescue_keys)
    blocked_columns = set(map(int, blocked_columns))
    blocked_events = set(blocked_events) - rescue_keys
    return {
        col_idx
        for col_idx in map(int, deferred_columns)
        if col_idx not in blocked_columns
        and not (set(path_usage[col_idx]) & blocked_events)
    }


def select_rescue_reactivation_columns(
    records: Sequence[Mapping],
    path_usage: Sequence[Mapping],
    retained_cent_fragment_chroms: Iterable[str],
    allowed_columns: Iterable[int],
) -> tuple[list[int], set[tuple[int, int]]]:
    eligible_events = eligible_rescue_event_keys(
        records, retained_cent_fragment_chroms
    )
    allowed = set(map(int, allowed_columns))
    columns = sorted(
        column_index
        for column_index in allowed
        if any(
            int(path_usage[column_index].get(event_key, 0)) > 0
            for event_key in eligible_events
        )
    )
    return columns, eligible_events


def plan_rescue_cf_swap(
    records: Sequence[Mapping],
    path_usage: Sequence[Mapping],
    eligible_event_keys: Iterable[tuple[int, int]],
    reactivated_columns: Iterable[int],
    cent_fragment_col2chrom: Mapping[int, str],
    active_columns: Iterable[int],
    rejected_event_keys: Iterable[tuple[int, int]] = (),
) -> RescueCFSwapPlan:
    """Plan an atomic rescue-on/cent-fragment-off retry.

    A rejected event removes every candidate path that carries it.  Every CF
    chromosome belonging to a rejected event is restored even when a retained
    rescue event also uses that chromosome.
    """

    eligible_events = {
        tuple(sorted(map(int, event_key)))
        for event_key in eligible_event_keys
    }
    reactivated = sorted(set(map(int, reactivated_columns)))
    carried_events = {
        event_key
        for column_index in reactivated
        for event_key in eligible_events
        if int(path_usage[column_index].get(event_key, 0)) > 0
    }
    rejected_events = {
        tuple(sorted(map(int, event_key)))
        for event_key in rejected_event_keys
    }
    unknown_rejected_events = rejected_events - carried_events
    if unknown_rejected_events:
        raise ValueError(
            "Cannot reject events outside the rescue swap candidates: "
            f"{sorted(unknown_rejected_events)}"
        )

    event_chroms = {
        tuple(sorted(map(int, record["event_key"]))): {
            str(chrom) for chrom in record["chroms"]
        }
        for record in records
    }
    missing_events = carried_events - set(event_chroms)
    if missing_events:
        raise ValueError(
            "Censat-pair rescue events are missing from the artifact: "
            f"{sorted(missing_events)}"
        )

    retained_events = carried_events - rejected_events
    retained_columns = tuple(
        column_index
        for column_index in reactivated
        if not any(
            int(path_usage[column_index].get(event_key, 0)) > 0
            for event_key in rejected_events
        )
    )
    events_without_standalone_carriers = {
        event_key
        for event_key in retained_events
        if not any(
            int(path_usage[column_index].get(event_key, 0)) > 0
            for column_index in retained_columns
        )
    }
    if events_without_standalone_carriers:
        raise ValueError(
            "Retained rescue events have no carrier after rejecting shared "
            f"paths: {sorted(events_without_standalone_carriers)}"
        )

    active_event_chroms = {
        chrom
        for event_key in retained_events
        for chrom in event_chroms[event_key]
    }
    restored_cf_chroms = {
        chrom
        for event_key in rejected_events
        for chrom in event_chroms[event_key]
    }
    disabled_cf_chroms = active_event_chroms - restored_cf_chroms

    active = set(map(int, active_columns))
    active_cf_col2chrom = {
        int(column_index): str(chrom)
        for column_index, chrom in cent_fragment_col2chrom.items()
        if int(column_index) in active
    }
    cf_columns = sorted(
        column_index
        for column_index, chrom in active_cf_col2chrom.items()
        if chrom in disabled_cf_chroms
    )

    covered_chroms = {
        active_cf_col2chrom[column_index]
        for column_index in cf_columns
    }
    missing_chroms = disabled_cf_chroms - covered_chroms
    if missing_chroms:
        raise ValueError(
            "No active cent_fragment column exists for rescue chromosomes: "
            f"{sorted(missing_chroms, key=chromosome_sort_key)}"
        )

    return RescueCFSwapPlan(
        event_keys=frozenset(retained_events),
        rescue_columns=retained_columns,
        active_event_chroms=frozenset(active_event_chroms),
        restored_cf_chroms=frozenset(restored_cf_chroms),
        disabled_cf_chroms=frozenset(disabled_cf_chroms),
        disabled_cf_columns=tuple(cf_columns),
    )


def classify_rescue_cf_swap_failures(
    records: Sequence[Mapping],
    active_event_keys: Iterable[tuple[int, int]],
    failed_chromosomes: Iterable[str],
) -> tuple[set[tuple[int, int]], set[str]]:
    """Return events to reject and failures unrelated to the active rescue."""

    active_events = {
        tuple(sorted(map(int, event_key)))
        for event_key in active_event_keys
    }
    event_chroms = {
        tuple(sorted(map(int, record["event_key"]))): {
            str(chrom) for chrom in record["chroms"]
        }
        for record in records
    }
    missing_events = active_events - set(event_chroms)
    if missing_events:
        raise ValueError(
            "Active rescue events are missing from the artifact: "
            f"{sorted(missing_events)}"
        )

    failed = set(map(str, failed_chromosomes))
    active_chroms = {
        chrom
        for event_key in active_events
        for chrom in event_chroms[event_key]
    }
    unrelated = failed - active_chroms
    rejected = {
        event_key
        for event_key in active_events
        if event_chroms[event_key] & failed
    }
    return rejected, unrelated


def prepare_reactivation_warm_start(
    current_weights: Sequence[float],
    current_global_indices: Sequence[int],
    reactivated_global_indices: Iterable[int],
    feature_count: int,
    deactivated_global_indices: Iterable[int] = (),
) -> tuple[list[int], np.ndarray, np.ndarray]:
    """Scatter local weights before atomically adding/removing global rows."""

    current_indices = list(map(int, current_global_indices))
    weights = np.asarray(current_weights)
    if len(weights) != len(current_indices):
        raise ValueError(
            "Current weight/index length mismatch: "
            f"{len(weights)} != {len(current_indices)}"
        )
    full_weights = np.zeros(int(feature_count), dtype=weights.dtype)
    full_weights[current_indices] = weights
    new_indices = sorted(
        (
            set(current_indices)
            | set(map(int, reactivated_global_indices))
        )
        - set(map(int, deactivated_global_indices))
    )
    return new_indices, full_weights[new_indices], full_weights
