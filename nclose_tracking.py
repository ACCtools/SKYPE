"""Shared NClose event, path-usage, and filter-provenance helpers.

The pipeline intentionally keeps the existing internal event keys:

* ordinary BND: ``(left_contig_idx, right_contig_idx)``
* step-11 INDEL: ``(event_type, event_idx, type2_merge_idx)``

This lets stages 22, 23, 24, and 31 consume one path-usage representation
without translating the keys used by the existing filtering code.
"""

from __future__ import annotations

import csv
import glob
import os
import pickle
from collections import Counter
from typing import Iterable, Mapping, Sequence


EVENT_CATALOG_PKL = "nclose_event_catalog.pkl"
PATH_USAGE_PKL = "nclose_path_usage.pkl"
FILTER_STATUS_PKL = "nclose_filter_status.pkl"
TYPE4_INDEL_GRAPH_EDGE_PKL = "type4_indel_graph_edges.pkl"
VCF_TYPE4_OUTLIER_INDEX_PKL = "vcf_type4_outlier_index.pkl"
INDEL_MERGE_TOLERANCE = 10_000
NCLOSE_ID_PREFIX = "SKYPE.nclose."

REPORT_COLUMNS = (
    "nclose_id",
    "start_chr", "start_pos", "start_dir",
    "end_chr", "end_pos", "end_dir",
    "nclose_cn", "nclose_cn_reason",
    "nclose_filter", "nclose_filter_reason",
    "nclose_cluster", "nclose_cluster_reason",
)

CTG_NAM = 0
CTG_DIR = 4
CHR_NAM = 5
CHR_STR = 7
CHR_END = 8


def _artifact_path(prefix: str, filename: str) -> str:
    return os.path.join(prefix, filename)


def _require_artifact(prefix: str, filename: str) -> str:
    path = _artifact_path(prefix, filename)
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"Missing NClose tracking artifact: {path}. "
            "Rerun the SKYPE pipeline from stage 02 so NClose filter "
            "provenance is not inferred from zero weights."
        )
    return path


def save_event_catalog(prefix: str, catalog: Sequence[dict]) -> None:
    with open(_artifact_path(prefix, EVENT_CATALOG_PKL), "wb") as handle:
        pickle.dump(list(catalog), handle)


def load_event_catalog(prefix: str) -> list[dict]:
    with open(_require_artifact(prefix, EVENT_CATALOG_PKL), "rb") as handle:
        catalog = pickle.load(handle)
    validate_event_catalog(catalog)
    return catalog


def save_path_usage(prefix: str, path_usage: Sequence[Counter]) -> None:
    with open(_artifact_path(prefix, PATH_USAGE_PKL), "wb") as handle:
        pickle.dump([Counter(usage) for usage in path_usage], handle)


def load_path_usage(prefix: str, expected_len: int | None = None) -> list[Counter]:
    with open(_require_artifact(prefix, PATH_USAGE_PKL), "rb") as handle:
        path_usage = [Counter(usage) for usage in pickle.load(handle)]
    if expected_len is not None and len(path_usage) != int(expected_len):
        raise ValueError(
            "NClose path-usage column count mismatch: "
            f"expected {expected_len}, found {len(path_usage)}. "
            "Rerun SKYPE from stage 22."
        )
    return path_usage


def save_filter_status(prefix: str, status: dict) -> None:
    with open(_artifact_path(prefix, FILTER_STATUS_PKL), "wb") as handle:
        pickle.dump(status, handle)


def load_filter_status(prefix: str) -> dict:
    with open(_require_artifact(prefix, FILTER_STATUS_PKL), "rb") as handle:
        status = pickle.load(handle)
    if status.get("version") != 1 or "stages" not in status:
        raise ValueError(
            "Unsupported NClose filter-status schema. Rerun SKYPE from stage 22."
        )
    return status


def corrected_direction(is_forward: bool, direction: str) -> str:
    if direction not in {"+", "-"}:
        raise ValueError(f"Unknown path direction: {direction}")
    if is_forward:
        return direction
    return "+" if direction == "-" else "-"


def _build_bnd_event(pair: Sequence[int], contig_data: Sequence, source: str) -> dict:
    pair = tuple(map(int, pair))
    if len(pair) != 2:
        raise ValueError(f"Invalid BND nclose pair: {pair}")
    contig_a = contig_data[pair[0]]
    contig_b = contig_data[pair[1]]
    is_forward = pair[0] < pair[1]
    return {
        "event_key": tuple(sorted(pair)),
        "kind": "bnd",
        "start_chr": contig_a[CHR_NAM],
        # Deliberately mirror the third field in each compressed-txt list.
        "start_pos": int(contig_a[CHR_STR]),
        "start_dir": corrected_direction(is_forward, contig_a[CTG_DIR]),
        "end_chr": contig_b[CHR_NAM],
        "end_pos": int(contig_b[CHR_STR]),
        "end_dir": corrected_direction(is_forward, contig_b[CTG_DIR]),
        "contig_name": contig_a[CTG_NAM],
        "source": source,
    }


def build_bnd_event_catalog(nclose_type: Mapping, contig_data: Sequence) -> list[dict]:
    """Build BND metadata in exactly the order used by the compressed txt."""
    catalog = []
    seen = set()
    for _chrom_pair, pair_list in nclose_type.items():
        for pair_value in pair_list:
            pair = tuple(pair_value)
            event = _build_bnd_event(
                pair, contig_data, "compressed_nclose_nodes_list.txt"
            )
            event_key = event["event_key"]
            if event_key in seen:
                raise ValueError(f"Duplicate compressed BND event: {event_key}")
            seen.add(event_key)
            catalog.append(event)
    validate_event_catalog(catalog)
    return catalog


def ecdna_circuit_event_keys(circuit: Sequence[int]) -> tuple[tuple, tuple]:
    if len(circuit) != 4:
        raise ValueError(f"Unexpected ecDNA circuit shape: {circuit}")
    s1, e1, s2, e2 = map(int, circuit)
    return tuple(sorted((s1, e1))), tuple(sorted((s2, e2)))


def append_ecdna_events_to_catalog(
    prefix: str,
    circuits: Sequence[Sequence[int]],
    contig_data: Sequence,
    circuit_indices: Sequence[int] | None = None,
) -> list[dict]:
    """Append missing ecDNA inversion NCloses after compressed BND rows.

    The two inversion pairs making each Amplicon need stable report IDs even
    when one was removed from the final compressed graph. Existing compressed
    events retain their original order and IDs.
    """
    catalog = load_event_catalog(prefix)
    event_by_key = event_catalog_by_key(catalog)
    if circuit_indices is None:
        circuit_indices = range(len(circuits))
    circuit_indices = list(map(int, circuit_indices))
    if len(circuit_indices) != len(circuits):
        raise ValueError(
            "ecDNA circuit/index count mismatch: "
            f"{len(circuits)} != {len(circuit_indices)}"
        )
    for circuit_idx, circuit in zip(circuit_indices, circuits):
        for pair_idx, event_key in enumerate(
            ecdna_circuit_event_keys(circuit), start=1
        ):
            event = event_by_key.get(event_key)
            ref = (int(circuit_idx) + 1, int(pair_idx))
            if event is not None:
                refs = event.setdefault("ecdna_circuit_refs", [])
                if ref not in refs:
                    refs.append(ref)
                continue
            event = _build_bnd_event(
                event_key,
                contig_data,
                "ecdna_circuit_data.pkl",
            )
            event["ecdna_circuit_refs"] = [ref]
            catalog.append(event)
            event_by_key[event_key] = event
    validate_event_catalog(catalog)
    save_event_catalog(prefix, catalog)
    return catalog


def replace_catalog_ecdna_events(
    prefix: str,
    circuits: Sequence[Sequence[int]],
    contig_data: Sequence,
    circuit_indices: Sequence[int] | None = None,
) -> list[dict]:
    """Replace ecDNA-only catalog rows while preserving base event IDs."""
    catalog = []
    for original_event in load_event_catalog(prefix):
        if original_event.get("source") == "ecdna_circuit_data.pkl":
            continue
        event = dict(original_event)
        event.pop("ecdna_circuit_refs", None)
        catalog.append(event)
    save_event_catalog(prefix, catalog)
    return append_ecdna_events_to_catalog(
        prefix,
        circuits,
        contig_data,
        circuit_indices=circuit_indices,
    )


def _read_paf_rows(path: str) -> list[list]:
    rows = []
    with open(path, "rt") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) <= CHR_END:
                raise ValueError(f"Malformed PAF row in {path}: {line.rstrip()}")
            row[CHR_STR] = int(row[CHR_STR])
            row[CHR_END] = int(row[CHR_END])
            rows.append(row)
    return rows


def _primary_indel_paf(folder: str, event_idx: int) -> tuple[str, int]:
    ordinary = os.path.join(folder, f"{event_idx}.paf")
    type2_paths = sorted(glob.glob(os.path.join(
        folder, f"{event_idx}_type2_merge_*.paf"
    )))
    candidates = ([ordinary] if os.path.isfile(ordinary) else []) + type2_paths
    if len(candidates) != 1:
        raise ValueError(
            f"Expected one primary step-11 PAF for event {event_idx} in "
            f"{folder}, found {candidates}. Rerun SKYPE from stage 11."
        )
    if not type2_paths:
        return candidates[0], -1
    suffix = os.path.basename(type2_paths[0]).removesuffix(".paf").rsplit("_", 1)[-1]
    if not suffix.isdigit():
        raise ValueError(f"Malformed type2-merge PAF name: {type2_paths[0]}")
    return type2_paths[0], int(suffix)


def indel_event_key(event_type: str, event_idx: int, type2_merge_idx: int = -1):
    if event_type not in {"front_jump", "back_jump"}:
        raise ValueError(f"Unknown step-11 INDEL type: {event_type}")
    return event_type, int(event_idx), int(type2_merge_idx)


def _vcf_outlier_index(prefix: str) -> dict:
    path = _artifact_path(prefix, VCF_TYPE4_OUTLIER_INDEX_PKL)
    if not os.path.isfile(path):
        return {}
    with open(path, "rb") as handle:
        return pickle.load(handle)


def discover_step11_indel_events(prefix: str) -> list[dict]:
    """Return one ordered event per emitted step-11 base PAF."""
    outlier_root = os.path.join(prefix, "11_ref_ratio_outliers")
    vcf_index = _vcf_outlier_index(prefix)
    events = []
    report_index = 0

    for event_type in ("front_jump", "back_jump"):
        folder = os.path.join(outlier_root, event_type)
        base_paths = glob.glob(os.path.join(folder, "*_base.paf"))
        base_by_idx = {}
        for base_path in base_paths:
            stem = os.path.basename(base_path).removesuffix("_base.paf")
            if not stem.isdigit():
                raise ValueError(f"Unexpected step-11 base PAF: {base_path}")
            base_by_idx[int(stem)] = base_path

        indices = sorted(base_by_idx)
        expected = list(range(1, len(indices) + 1))
        if indices != expected:
            raise ValueError(
                f"Non-contiguous step-11 {event_type} indices: "
                f"expected {expected}, found {indices}."
            )

        for event_idx in indices:
            primary_path, type2_merge_idx = _primary_indel_paf(folder, event_idx)
            primary_rows = _read_paf_rows(primary_path)
            base_rows = _read_paf_rows(base_by_idx[event_idx])
            metadata = (
                vcf_index.get((event_type, event_idx))
                or vcf_index.get(f"{event_type}:{event_idx}")
                or {}
            )

            if primary_rows:
                first = primary_rows[0]
                last = primary_rows[-1]
                start_dir = first[CTG_DIR]
                end_dir = last[CTG_DIR]
                start_pos = first[CHR_END] if start_dir == "+" else first[CHR_STR]
                end_pos = last[CHR_STR] if end_dir == "+" else last[CHR_END]
                start_chr = first[CHR_NAM]
                end_chr = last[CHR_NAM]
            else:
                if not metadata:
                    raise ValueError(
                        f"Empty primary step-11 PAF without VCF metadata: {primary_path}"
                    )
                lo, hi = sorted((int(metadata["st"]), int(metadata["nd"])))
                start_chr = end_chr = metadata["chrom"]
                start_dir = end_dir = "+"
                if event_type == "front_jump":
                    start_pos, end_pos = lo, hi
                else:
                    start_pos, end_pos = hi, lo

            if metadata:
                span_chrom = metadata.get("chrom", start_chr)
                span_st, span_nd = sorted((int(metadata["st"]), int(metadata["nd"])))
            elif base_rows:
                span_chrom = base_rows[0][CHR_NAM]
                same_chrom_rows = [row for row in base_rows if row[CHR_NAM] == span_chrom]
                span_st = min(
                    min(row[CHR_STR], row[CHR_END]) for row in same_chrom_rows
                )
                span_nd = max(
                    max(row[CHR_STR], row[CHR_END]) for row in same_chrom_rows
                )
            else:
                span_chrom = start_chr
                span_st, span_nd = sorted((int(start_pos), int(end_pos)))

            event = {
                "event_key": indel_event_key(event_type, event_idx, type2_merge_idx),
                "kind": "indel",
                "event_type": event_type,
                "indel_kind": "deletion" if event_type == "front_jump" else "insertion",
                "start_chr": start_chr,
                "start_pos": int(start_pos),
                "start_dir": start_dir,
                "end_chr": end_chr,
                "end_pos": int(end_pos),
                "end_dir": end_dir,
                "chrom": span_chrom,
                "st": int(span_st),
                "nd": int(span_nd),
                "event_idx": int(event_idx),
                "type2_merge_idx": int(type2_merge_idx),
                "primary_paf": primary_path,
                "base_paf": base_by_idx[event_idx],
                "report_index": report_index,
                "source": metadata.get("source", f"INDEL_INDEX_{report_index}"),
            }
            for field in (
                "event_id", "vcf_id", "mate_id", "merged_vcf_ids", "svtype", "nodes"
            ):
                if field in metadata:
                    event[field] = metadata[field]
            events.append(event)
            report_index += 1

    return events


def replace_catalog_indels(prefix: str) -> list[dict]:
    catalog = [
        event for event in load_event_catalog(prefix)
        if event["kind"] == "bnd"
        and event.get("source") == "compressed_nclose_nodes_list.txt"
    ]
    catalog.extend(discover_step11_indel_events(prefix))
    validate_event_catalog(catalog)
    save_event_catalog(prefix, catalog)
    return catalog


def validate_event_catalog(catalog: Sequence[dict]) -> None:
    seen = set()
    required = {
        "event_key", "kind", "start_chr", "start_pos", "start_dir",
        "end_chr", "end_pos", "end_dir",
    }
    for event in catalog:
        missing = required - set(event)
        if missing:
            raise ValueError(f"NClose event is missing fields {sorted(missing)}: {event}")
        key = event["event_key"]
        if key in seen:
            raise ValueError(f"Duplicate NClose catalog key: {key}")
        seen.add(key)
        if event["kind"] not in {"bnd", "indel"}:
            raise ValueError(f"Unknown NClose event kind: {event['kind']}")
        if event["start_dir"] not in {"+", "-"} or event["end_dir"] not in {"+", "-"}:
            raise ValueError(f"Invalid NClose event direction: {event}")


def event_catalog_by_key(catalog: Sequence[dict]) -> dict:
    return {event["event_key"]: event for event in catalog}


def bnd_event_keys(catalog: Sequence[dict]) -> set:
    return {event["event_key"] for event in catalog if event["kind"] == "bnd"}


def compressed_bnd_event_keys(catalog: Sequence[dict]) -> set:
    return {
        event["event_key"]
        for event in catalog
        if event["kind"] == "bnd"
        and event.get("source") == "compressed_nclose_nodes_list.txt"
    }


def indel_event_keys(catalog: Sequence[dict]) -> set:
    return {event["event_key"] for event in catalog if event["kind"] == "indel"}


def load_type4_edge_event_map(
    prefix: str,
    catalog: Sequence[dict],
    tolerance: int = INDEL_MERGE_TOLERANCE,
) -> dict:
    edge_path = _artifact_path(prefix, TYPE4_INDEL_GRAPH_EDGE_PKL)
    if not os.path.isfile(edge_path):
        return {}
    with open(edge_path, "rb") as handle:
        edges = pickle.load(handle)
    if not edges:
        return {}

    indels = [event for event in catalog if event["kind"] == "indel"]
    edge_to_event = {}
    for edge in edges:
        indel_kind = edge.get("indel_kind")
        if indel_kind not in {"deletion", "insertion"}:
            raise ValueError(f"Unknown type4 graph INDEL kind: {indel_kind}")
        event_type = "front_jump" if indel_kind == "deletion" else "back_jump"
        chrom = edge.get("base_chrom") or edge.get("chrom")
        st, nd = sorted((
            int(edge.get("base_st", edge.get("span_st", 0))),
            int(edge.get("base_nd", edge.get("span_nd", 0))),
        ))
        candidates = [
            event for event in indels
            if event.get("event_type") == event_type
            and event.get("chrom") == chrom
            and abs(int(event.get("st", 0)) - st) <= tolerance
            and abs(int(event.get("nd", 0)) - nd) <= tolerance
        ]
        if not candidates:
            raise ValueError(
                "Type4 graph edge has no matching emitted step-11 INDEL: "
                f"{event_type} {chrom}:{st}-{nd}. Rerun stages 02 and 11."
            )
        candidates.sort(key=lambda event: (
            abs(int(event["st"]) - st) + abs(int(event["nd"]) - nd),
            event["event_key"],
        ))
        edge_key = (tuple(edge["src"]), tuple(edge["dst"]))
        event_key = candidates[0]["event_key"]
        previous = edge_to_event.setdefault(edge_key, event_key)
        if previous != event_key:
            raise ValueError(f"Type4 graph edge maps to multiple INDELs: {edge_key}")
    return edge_to_event


def count_index_path_events(
    index_path: Sequence,
    bnd_keys: set,
    type4_edge_to_event: Mapping | None = None,
) -> Counter:
    usage = Counter()
    s = 1
    while s < len(index_path) - 2:
        left = index_path[s][1]
        right = index_path[s + 1][1]
        if isinstance(left, int) and isinstance(right, int):
            candidate = tuple(sorted((left, right)))
            if candidate in bnd_keys:
                usage[candidate] += 1
                s += 2
                continue
        s += 1

    if type4_edge_to_event:
        compact_path = [tuple(list(node)[:2]) for node in index_path]
        for left, right in zip(compact_path, compact_path[1:]):
            event_key = type4_edge_to_event.get((left, right))
            if event_key is not None:
                usage[event_key] += 1
    return usage


def count_ecdna_circuit_events(circuit: Sequence[int], bnd_keys: set) -> Counter:
    usage = Counter()
    for event_key in ecdna_circuit_event_keys(circuit):
        if event_key in bnd_keys:
            usage[event_key] += 1
    return usage


def calculate_event_weights(
    path_usage: Sequence[Mapping],
    weights: Sequence,
    active_columns: Iterable[int] | None = None,
) -> dict:
    if len(path_usage) != len(weights):
        raise ValueError(
            f"NClose usage/weight length mismatch: {len(path_usage)} != {len(weights)}"
        )
    active = None if active_columns is None else set(map(int, active_columns))
    event_weights = Counter()
    for col_idx, usage in enumerate(path_usage):
        if active is not None and col_idx not in active:
            continue
        weight = float(weights[col_idx])
        for event_key, count in usage.items():
            event_weights[event_key] += int(count) * weight
    return dict(event_weights)


def carrier_columns(path_usage: Sequence[Mapping], event_key) -> set[int]:
    return {
        col_idx for col_idx, usage in enumerate(path_usage)
        if int(usage.get(event_key, 0)) > 0
    }


def _has_active_carrier(path_usage: Sequence[Mapping], event_key, active: set[int]) -> bool:
    return any(
        col_idx in active and int(usage.get(event_key, 0)) > 0
        for col_idx, usage in enumerate(path_usage)
    )


def initialise_filter_status(
    catalog: Sequence[dict],
    path_usage: Sequence[Mapping],
    active_columns: Iterable[int],
    explicit_reasons: Mapping | None = None,
) -> dict:
    active = set(map(int, active_columns))
    explicit_reasons = dict(explicit_reasons or {})
    reasons = {}
    for event in catalog:
        event_key = event["event_key"]
        if event_key in explicit_reasons:
            reasons[event_key] = explicit_reasons[event_key]
            continue
        all_carriers = carrier_columns(path_usage, event_key)
        if not all_carriers:
            reasons[event_key] = "FILTERED_02_NO_ELIGIBLE_PATH"
        elif not (all_carriers & active):
            reasons[event_key] = "FILTERED_02_COFILTERED_PATH"
    return {
        "version": 1,
        "event_keys": [event["event_key"] for event in catalog],
        "stages": {
            "initial": {
                "active_columns": sorted(active),
                "reasons": reasons,
            }
        },
    }


def record_filter_stage(
    status: dict,
    stage: str,
    previous_stage: str,
    path_usage: Sequence[Mapping],
    active_columns: Iterable[int],
    direct_reasons: Mapping | None,
    cofiltered_reason: str,
    forced_reasons: Mapping | None = None,
) -> dict:
    if previous_stage not in status["stages"]:
        raise ValueError(f"Missing previous NClose filter stage: {previous_stage}")
    active = set(map(int, active_columns))
    previous_active = set(
        map(int, status["stages"][previous_stage].get("active_columns", []))
    )
    if not active <= previous_active:
        reintroduced = sorted(active - previous_active)
        raise ValueError(
            f"NClose filter stage {stage} reintroduced columns removed by "
            f"{previous_stage}: {reintroduced[:20]}"
        )
    previous_reasons = status["stages"][previous_stage].get("reasons", {})
    direct_reasons = dict(direct_reasons or {})
    forced_reasons = dict(forced_reasons or {})
    reasons = {}
    for event_key in status["event_keys"]:
        if event_key in previous_reasons:
            reasons[event_key] = previous_reasons[event_key]
        # Cluster thresholding is event-level provenance: an exact-zero
        # carrier can remain selected while another positive carrier is cut.
        elif event_key in forced_reasons:
            reasons[event_key] = forced_reasons[event_key]
        elif _has_active_carrier(path_usage, event_key, active):
            continue
        else:
            reasons[event_key] = direct_reasons.get(event_key, cofiltered_reason)
    status["stages"][stage] = {
        "active_columns": sorted(active),
        "reasons": reasons,
    }
    return status


def stage_reason(status: Mapping, stage: str, event_key) -> str:
    return status["stages"][stage].get("reasons", {}).get(event_key, "PASS")


def reconcile_filter_status_catalog(
    status: Mapping,
    catalog: Sequence[dict],
    path_usage: Sequence[Mapping],
) -> dict:
    """Reconcile late report-only events with already completed filter stages."""
    event_keys = [event["event_key"] for event in catalog]
    reconciled = {
        "version": 1,
        "event_keys": event_keys,
        "stages": {},
    }
    stage_defaults = {
        "initial": "FILTERED_02_COFILTERED_PATH",
        "base": "FILTERED_23_COFILTERED_PATH",
        "filter": "FILTERED_24_COFILTERED_PATH",
        "cluster": "FILTERED_24_CLUSTER_COFILTERED_PATH",
    }
    previous_reasons = {}
    for stage in ("initial", "base", "filter", "cluster"):
        if stage not in status.get("stages", {}):
            continue
        original_stage = status["stages"][stage]
        original_reasons = original_stage.get("reasons", {})
        active = set(map(int, original_stage.get("active_columns", [])))
        reasons = {}
        for event_key in event_keys:
            if event_key in previous_reasons:
                reasons[event_key] = previous_reasons[event_key]
            elif event_key in original_reasons:
                reasons[event_key] = original_reasons[event_key]
            elif _has_active_carrier(path_usage, event_key, active):
                continue
            elif stage == "initial" and not carrier_columns(path_usage, event_key):
                reasons[event_key] = "FILTERED_02_NO_ELIGIBLE_PATH"
            else:
                reasons[event_key] = stage_defaults[stage]
        reconciled["stages"][stage] = {
            "active_columns": sorted(active),
            "reasons": reasons,
        }
        previous_reasons = reasons
    return reconciled


def nclose_event_id(index: int) -> str:
    index = int(index)
    if index < 1:
        raise ValueError(f"NClose ID index must be 1-based: {index}")
    return f"{NCLOSE_ID_PREFIX}{index}"


def nclose_event_id_by_key(catalog: Sequence[dict]) -> dict:
    return {
        event["event_key"]: nclose_event_id(index)
        for index, event in enumerate(catalog, start=1)
    }


def format_nclose_ids(event_keys: Sequence, event_id_by_key: Mapping) -> str:
    missing = [event_key for event_key in event_keys if event_key not in event_id_by_key]
    if missing:
        raise ValueError(f"NClose IDs are missing for event keys: {missing}")
    return ";".join(event_id_by_key[event_key] for event_key in event_keys)


def bed_visible_ecdna_indices(
    weights: Sequence,
    ecdna_column_by_index: Mapping[int, int],
    raw_weight_cutoff: float,
) -> list[int]:
    """Select one BED stage's Amplicons using the strict raw-weight gate."""
    return sorted(
        int(ecdna_idx)
        for ecdna_idx, column_idx in ecdna_column_by_index.items()
        if float(weights[int(column_idx)]) > float(raw_weight_cutoff)
    )


def bed_visible_ecdna_indices_across_stages(
    stage_weights: Mapping[str, Sequence],
    ecdna_column_by_index: Mapping[int, int],
    raw_weight_cutoff: float,
) -> list[int]:
    """Return the Amplicon union displayed by any generated BED stage."""
    visible_indices = set()
    for weights in stage_weights.values():
        visible_indices.update(bed_visible_ecdna_indices(
            weights,
            ecdna_column_by_index,
            raw_weight_cutoff,
        ))
    return sorted(visible_indices)


def _format_cn(value: float) -> str:
    value = float(value)
    if abs(value) < 1e-12:
        return "0"
    return f"{value:.10g}"


def build_nclose_report_rows(
    catalog: Sequence[dict],
    path_usage: Sequence[Mapping],
    status: Mapping,
    n_unit: float,
    stage_weights: Mapping[str, Sequence],
    run_filter_cluster: bool,
) -> list[dict]:
    """Build the fixed 12-column report without inferring provenance."""
    validate_event_catalog(catalog)
    catalog_keys = [event["event_key"] for event in catalog]
    if list(status.get("event_keys", [])) != catalog_keys:
        raise ValueError(
            "NClose catalog/filter-status event order mismatch. "
            "Rerun the pipeline from stage 22."
        )
    n_unit = float(n_unit)
    if not n_unit > 0:
        raise ValueError(f"NClose CN normalisation unit must be positive: {n_unit}")
    if "base" not in stage_weights:
        raise ValueError("Missing base weights for NClose report")

    required_stages = ["base"]
    if run_filter_cluster:
        required_stages.extend(["filter", "cluster"])
    aggregates = {}
    for stage in required_stages:
        if stage not in status.get("stages", {}):
            raise ValueError(
                f"Missing NClose filter-status stage {stage}. "
                "Rerun the pipeline from the corresponding NNLS stage."
            )
        if stage not in stage_weights:
            raise ValueError(f"Missing {stage} weights for NClose report")
        aggregates[stage] = calculate_event_weights(
            path_usage, stage_weights[stage]
        )

    rows = []
    for event_index, event in enumerate(catalog, start=1):
        event_key = event["event_key"]
        row = {
            "nclose_id": nclose_event_id(event_index),
            "start_chr": event["start_chr"],
            "start_pos": str(int(event["start_pos"])),
            "start_dir": event["start_dir"],
            "end_chr": event["end_chr"],
            "end_pos": str(int(event["end_pos"])),
            "end_dir": event["end_dir"],
        }
        for stage, value_col, reason_col in (
            ("base", "nclose_cn", "nclose_cn_reason"),
            ("filter", "nclose_filter", "nclose_filter_reason"),
            ("cluster", "nclose_cluster", "nclose_cluster_reason"),
        ):
            if stage != "base" and not run_filter_cluster:
                row[value_col] = "NA"
                row[reason_col] = "NOT_RUN_VARIANT_MODE"
                continue
            reason = stage_reason(status, stage, event_key)
            if reason != "PASS":
                row[value_col] = "-1"
                row[reason_col] = reason
            else:
                raw_weight = aggregates[stage].get(event_key, 0.0)
                row[value_col] = _format_cn(raw_weight / n_unit)
                row[reason_col] = "PASS"
        rows.append(row)
    return rows


def write_nclose_report(
    prefix: str,
    catalog: Sequence[dict],
    path_usage: Sequence[Mapping],
    status: Mapping,
    n_unit: float,
    stage_weights: Mapping[str, Sequence],
    run_filter_cluster: bool,
) -> str:
    rows = build_nclose_report_rows(
        catalog,
        path_usage,
        status,
        n_unit,
        stage_weights,
        run_filter_cluster,
    )
    output_path = _artifact_path(prefix, "nclose_report.tsv")
    with open(output_path, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=REPORT_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    return output_path
