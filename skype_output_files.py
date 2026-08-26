from __future__ import annotations

import glob
import os


DEPTH_STAT_SUFFIX = ".win.stat.gz"


def depth_stat_to_paf_path(depth_path):
    """Return the PAF paired with one stage-21 depth-stat artifact."""

    depth_path = os.fspath(depth_path)
    if not depth_path.endswith(DEPTH_STAT_SUFFIX):
        raise ValueError(
            f"Expected a {DEPTH_STAT_SUFFIX} depth path, found: {depth_path}"
        )
    return depth_path[:-len(DEPTH_STAT_SUFFIX)] + ".paf"


def build_matrix_column_locations(
    prefix,
    paf_ans_list,
    front_depth_inputs,
    back_depth_inputs,
    ecdna_depth_inputs,
    cen_fragment_list,
):
    """Build stage-31 location metadata in exact stage-22 column order."""

    locations = [path for path, _ in paf_ans_list]
    locations.extend(
        depth_stat_to_paf_path(base_depth_path)
        for _, _, base_depth_path, _, _ in front_depth_inputs
    )
    locations.extend(
        depth_stat_to_paf_path(base_depth_path)
        for _, _, base_depth_path, _, _ in back_depth_inputs
    )
    locations.extend(
        depth_stat_to_paf_path(depth_path)
        for _, depth_path in ecdna_depth_inputs
    )
    for chrom, info in cen_fragment_list:
        side = "right" if info["dir"] else "left"
        locations.append(
            os.path.join(
                os.fspath(prefix),
                "12_cent_fragment",
                chrom,
                f"{side}.fragment",
            )
        )
    return locations


def _validate_contiguous_depth_indices(index_to_path, folder_path, file_kind):
    indices = sorted(index_to_path)
    expected = list(range(1, len(indices) + 1))
    if indices != expected:
        raise ValueError(
            f"Non-contiguous {file_kind} indices in {folder_path}: "
            f"expected {expected}, found {indices}. Rerun SKYPE from stage 11 "
            "to replace stale intermediate files."
        )
    return indices


def discover_jump_depth_inputs(folder_path, type2_folder_path):
    """Return one validated depth-input record per jump event."""
    base_suffix = "_base.win.stat.gz"
    base_paths = glob.glob(os.path.join(folder_path, f"*{base_suffix}"))
    base_by_index = {}
    for base_path in base_paths:
        index_text = os.path.basename(base_path)[:-len(base_suffix)]
        if not index_text.isdigit():
            raise ValueError(f"Unexpected jump base-depth file: {base_path}")
        base_by_index[int(index_text)] = base_path

    indices = _validate_contiguous_depth_indices(
        base_by_index, folder_path, "jump depth"
    )
    records = []
    for event_idx in indices:
        ordinary_path = os.path.join(
            folder_path, f"{event_idx}.win.stat.gz"
        )
        type2_paths = sorted(glob.glob(os.path.join(
            folder_path, f"{event_idx}_type2_merge_*.win.stat.gz"
        )))
        primary_paths = (
            ([ordinary_path] if os.path.isfile(ordinary_path) else [])
            + type2_paths
        )
        if len(primary_paths) != 1:
            raise ValueError(
                f"Expected exactly one primary depth file for jump event "
                f"{event_idx} in {folder_path}, found {primary_paths}. "
                "This usually means outputs from different runs were mixed; "
                "rerun SKYPE from stage 11."
            )

        primary_path = primary_paths[0]
        type2_idx = -1
        type2_path = None
        if type2_paths:
            type2_stem = os.path.basename(primary_path).removesuffix(
                ".win.stat.gz"
            )
            type2_index_text = type2_stem.rsplit("_", 1)[-1]
            if not type2_index_text.isdigit():
                raise ValueError(f"Malformed type2 merge depth file: {primary_path}")
            type2_idx = int(type2_index_text)
            type2_path = os.path.join(
                type2_folder_path, f"{type2_idx}.win.stat.gz"
            )
            if not os.path.isfile(type2_path):
                raise FileNotFoundError(
                    f"Missing type2 insertion depth file paired with "
                    f"{primary_path}: {type2_path}"
                )

        records.append((
            event_idx,
            primary_path,
            base_by_index[event_idx],
            type2_idx,
            type2_path,
        ))
    return records


def discover_ecdna_depth_inputs(folder_path):
    suffix = ".win.stat.gz"
    depth_by_index = {}
    for depth_path in glob.glob(os.path.join(folder_path, f"*{suffix}")):
        index_text = os.path.basename(depth_path)[:-len(suffix)]
        if not index_text.isdigit():
            raise ValueError(f"Unexpected ecDNA depth file: {depth_path}")
        depth_by_index[int(index_text)] = depth_path

    indices = _validate_contiguous_depth_indices(
        depth_by_index, folder_path, "ecDNA depth"
    )
    return [(event_idx, depth_by_index[event_idx]) for event_idx in indices]
