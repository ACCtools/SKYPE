"""Read every retained primitive and its unchanged, conditional model dosage.

This module never builds a graph, refits weights, or selects on local evidence.
An optional depth matrix is checked against the completed model numerically;
legacy files do not supply a historical matrix-row/source generation binding.
"""
from __future__ import annotations

from collections import Counter, defaultdict
import csv
import hashlib
from math import fsum, isfinite
from pathlib import Path
import pickle

import h5py
import numpy as np

from native_vcf_bnd import nclose_bnd_memberships
from structure_nclose import MODEL_VERSION, StructureWeights, source_signature


SOURCE_FILES = (
    "tot_loc_list.pkl", "path_data.pkl", "contig_pat_vec_data.pkl",
    "nclose_event_catalog.pkl", "ecdna_circuit_data.pkl",
    "conjoined_type4_ins_del.pkl", "type4_indel_graph_edges.pkl",
    "01_nclose_data.pkl", "telomere_connected_list.txt",
)
MATRIX_CONTRACT = "depth_only_v1"
BLOCK_ROWS = 128
_IMPLEMENTATION_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def _stat(path):
    if not path.exists():
        return None
    info = path.stat()
    return (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns)


class _Inputs:
    """Hash inputs before reading, then reject changes during this extraction."""
    def __init__(self):
        self.records, self.stats = {}, {}

    def add(self, path, required=True):
        path = Path(path).resolve()
        if path in self.stats:
            if required and self.stats[path] is None:
                raise ValueError(f"Completed native artifact is missing: {path}")
            return path
        before = _stat(path)
        if before is None and required:
            raise ValueError(f"Completed native artifact is missing: {path}")
        digest = None
        if before is not None:
            value = hashlib.sha256()
            with path.open("rb") as handle:
                for block in iter(lambda: handle.read(8 << 20), b""):
                    value.update(block)
            digest = value.hexdigest()
        if _stat(path) != before:
            raise ValueError(f"Input changed while hashing: {path}")
        self.stats[path] = before
        self.records[str(path)] = dict(path=str(path), sha256=digest,
                                       exists=before is not None)
        return path

    def pickle(self, path):
        with self.add(path).open("rb") as handle:
            return pickle.load(handle)

    def array(self, path):
        return np.load(self.add(path), allow_pickle=False)

    def finish(self):
        for path, before in self.stats.items():
            if _stat(path) != before:
                raise ValueError(f"Input changed during candidate extraction: {path}")
        return self.records


def _same(a, b, strip_newlines=False):
    """Compare saved and in-memory views without coercing geometry or weights."""
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_same(a[k], b[k], strip_newlines) for k in a)
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(_same(x, y, strip_newlines) for x, y in zip(a, b))
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        return np.array_equal(a, b)
    if strip_newlines and isinstance(a, str) and isinstance(b, str):
        return a.rstrip("\r\n") == b.rstrip("\r\n")
    return a == b


def _integer(value, label, minimum=None):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{label} must be an integer: {value!r}")
    result = int(value)
    if minimum is not None and result < minimum:
        raise ValueError(f"{label} must be at least {minimum}: {result}")
    return result


def _counts(value, label):
    return Counter({key: _integer(count, label, 1) for key, count in value.items()})


def _validate_model(model, weights, locations):
    if model.get("version") != MODEL_VERSION:
        raise ValueError("Unsupported completed native model version; extraction cannot rebuild it")
    column_count = _integer(model["column_count"], "Model column count", 0)
    if weights.shape != (column_count,) or len(locations) != column_count:
        raise ValueError("Saved model, weights and column-source order have different sizes")
    if not np.isfinite(weights).all() or np.any(weights < 0):
        raise ValueError("Saved NNLS weights must be finite and nonnegative")
    n_unit = float(model["n_unit"])
    if not isfinite(n_unit) or n_unit <= 0:
        raise ValueError("Completed model n_unit must be finite and positive")
    sources, aliases, ncloses = (model[name] for name in
                                 ("nclose_sources", "nclose_aliases", "ncloses"))
    if set(sources) != set(aliases) or any(key not in ncloses for key in aliases.values()):
        raise ValueError("Saved NClose source/alias registry is inconsistent")
    if len({_key(key) for key in sources}) != len(sources):
        raise ValueError("Saved NClose sources have colliding text keys")
    for key, event in sources.items():
        if tuple(event["event_key"]) != key or event["kind"] not in {"bnd", "indel"}:
            raise ValueError(f"Invalid saved NClose source: {key}")
    seen_ids, seen_columns = set(), set()
    for structure in model["structures"]:
        sid, column = structure["structure_id"], structure["feature_index"]
        if sid in seen_ids:
            raise ValueError(f"Duplicate saved structure ID: {sid}")
        seen_ids.add(sid)
        raw, norm = float(structure["raw_weight"]), float(structure["weight_N"])
        if not isfinite(raw) or not isfinite(norm) or raw < 0 or norm != raw / n_unit:
            raise ValueError(f"Saved structure weights/n_unit disagree: {sid}")
        if column is None:
            if structure["kind"] != "VIRTUAL_INV":
                raise ValueError(f"Nonvirtual structure has no fitted column: {sid}")
        else:
            column = _integer(column, "Feature index", 0)
            if column >= column_count or column in seen_columns:
                raise ValueError(f"Invalid or duplicate fitted feature index: {sid}")
            seen_columns.add(column)
            if raw != float(weights[column]):
                raise ValueError(f"Saved model coefficient differs from weight.npy: {sid}")
            if str(structure["source"]) != str(locations[column]):
                raise ValueError(f"Saved model column-source order differs from tot_loc_list: {sid}")
        incidence = _counts(structure["source_nclose_counts"], "Source occurrence count")
        canonical = Counter()
        for key, count in incidence.items():
            if key not in sources:
                raise ValueError(f"Structure refers to an absent NClose source: {sid}: {key}")
            canonical[aliases[key]] += count
        if canonical != _counts(structure["nclose_counts"], "NClose occurrence count"):
            raise ValueError(f"Saved source and canonical occurrence counts disagree: {sid}")
    if seen_columns != set(range(column_count)):
        raise ValueError("Completed model is missing fitted columns")


def _load_completed(prefix, inputs, context=None, nodes=None):
    prefix = Path(prefix)
    for name in SOURCE_FILES:
        inputs.add(prefix / name, required=False)
    model = inputs.pickle(prefix / "structure_nclose_model.pkl")
    if model.get("source_signature") != source_signature(prefix):
        raise ValueError("Saved model source_signature is stale; candidate extraction cannot rebuild it")
    weights = inputs.array(prefix / "weight.npy")
    locations = inputs.pickle(prefix / "tot_loc_list.pkl")
    saved_nodes = inputs.pickle(prefix / "01_nclose_data.pkl")["contig_data"]
    _validate_model(model, weights, locations)
    if context is not None:
        if (not _same(context.model, model) or float(context.n_unit) != model["n_unit"]
                or not _same(context.ncloses, model["ncloses"])
                or not _same(context.structures, {s["structure_id"]: s for s in model["structures"]})):
            raise ValueError("Provided native context differs from the completed saved model")
    if nodes is not None and not _same(nodes, saved_nodes, strip_newlines=True):
        raise ValueError("Provided alignment nodes differ from the completed saved nodes")
    return model, weights, saved_nodes


def read_completed_context(prefix):
    """Load only an already completed, internally consistent native model.

    StructureWeights joins private in-memory copies only after every saved
    coefficient and occurrence count has been checked. It cannot rebuild/save.
    """
    inputs = _Inputs()
    model, weights, nodes = _load_completed(prefix, inputs)
    context = StructureWeights(model, weights, model["n_unit"])
    inputs.finish()
    return context, nodes


def _matrix_counts(prefix, matrix_path, model, weights, inputs):
    path = inputs.add(matrix_path, required=False)
    if not path.exists():
        return None, dict(state="unavailable", reason="matrix_file_missing", path=str(path),
                          historical_generation_binding="unavailable")
    B = inputs.array(Path(prefix) / "B.npy")
    prediction = inputs.array(Path(prefix) / "predict_B.npy")
    meta = inputs.pickle(Path(prefix) / "23_input.pkl")
    if B.ndim != 1 or prediction.shape != B.shape:
        raise ValueError("Completed target and prediction shapes disagree")
    if not np.isfinite(B).all() or not np.isfinite(prediction).all():
        raise ValueError("Completed target or prediction is not finite")
    counts, sizes, errors, bounds = {}, {}, [], []
    active_count = int(np.count_nonzero(weights))
    # Two summations of the same products, with arbitrary block grouping.
    eps = np.finfo(np.float64).eps
    operations = 2 * active_count + 2
    gamma = operations * eps / (1 - operations * eps)
    with h5py.File(path, "r") as handle:
        contract = handle.attrs.get("matrix_contract", "")
        if isinstance(contract, bytes):
            contract = contract.decode()
        if contract != MATRIX_CONTRACT or meta.get("matrix_contract") != MATRIX_CONTRACT:
            raise ValueError("Candidate depth reporting requires the saved depth_only_v1 contract")
        if not {"A", "A_fail", "B", "B_fail"}.issubset(handle):
            raise ValueError("Depth matrix is missing A/A_fail/B/B_fail")
        offset = 0
        for name, target_name in (("A", "B"), ("A_fail", "B_fail")):
            data, target = handle[name], np.asarray(handle[target_name])
            if (data.ndim != 2 or data.shape[0] != model["column_count"]
                    or target.shape != (data.shape[1],)):
                raise ValueError(f"Depth matrix {name} shape differs from the completed model/target")
            end = offset + data.shape[1]
            if not np.array_equal(target, B[offset:end]):
                raise ValueError(f"Depth matrix {target_name} differs from completed B.npy")
            sizes[name] = data.shape[1]
            fields = {field: np.zeros(data.shape[0], dtype=np.int64)
                      for field in ("nonzero", "negative", "positive")}
            reconstructed = np.zeros(data.shape[1], dtype=np.float64)
            absolute = np.zeros(data.shape[1], dtype=np.float64)
            for start in range(0, data.shape[0], BLOCK_ROWS):
                block = data[start:start + BLOCK_ROWS]
                if not np.isfinite(block).all():
                    raise ValueError(f"Depth matrix {name} contains nonfinite entries")
                sl = slice(start, start + len(block))
                fields["nonzero"][sl] = np.count_nonzero(block, axis=1)
                fields["negative"][sl] = np.count_nonzero(block < 0, axis=1)
                fields["positive"][sl] = np.count_nonzero(block > 0, axis=1)
                local = weights[sl]
                active = np.flatnonzero(local)
                if len(active):
                    products = np.asarray(block[active], dtype=np.float64) * local[active, None]
                    reconstructed += products.sum(axis=0, dtype=np.float64)
                    absolute += np.abs(products).sum(axis=0, dtype=np.float64)
            if not np.isfinite(reconstructed).all() or not np.isfinite(absolute).all():
                raise ValueError(f"Depth matrix {name} prediction arithmetic is nonfinite")
            allowance = 2 * gamma * absolute + 8 * eps * np.maximum(1, np.abs(reconstructed))
            error = np.abs(reconstructed - prediction[offset:end])
            if np.any(error > allowance):
                raise ValueError(f"Depth matrix {name} and saved weights do not reconstruct predict_B.npy")
            errors.append(float(error.max(initial=0)))
            bounds.append(float(allowance.max(initial=0)))
            counts[name] = fields
            offset = end
        if offset != len(B):
            raise ValueError("Depth matrix bin count differs from completed B.npy")
        if (meta.get("B_depth_start") != 0 or meta.get("B_depth_end") != sizes["A"]
                or len(meta.get("chr_filt_st_list", ())) != sizes["A"]):
            raise ValueError("Saved fitted-bin metadata differs from the depth matrix")
    return counts, dict(state="current_numerical_association_verified", path=str(path),
        historical_generation_binding="unavailable", fitted_bins=sizes["A"],
        excluded_bins=sizes["A_fail"], feature_count=model["column_count"],
        target_equality=True, stored_source_order_equality=True,
        prediction_max_absolute_error=max(errors), prediction_max_roundoff_bound=max(bounds),
        row_source_association_assumption="Feature-major row order follows tot_loc_list.pkl; the matrix has no embedded source-row generation binding.")


def _feature_catalog(model, observed):
    result = {}
    for s in model["structures"]:
        col = s["feature_index"]
        fields = {f"{sign}_{group}_bins": None for sign in ("nonzero", "negative", "positive")
                  for group in ("fitted", "excluded")}
        if col is None:
            state = "postfit_qualified_virtual_structure"
        elif observed is None:
            state = "unavailable"
        else:
            fields = {f"{sign}_{group}_bins": int(observed[name][sign][col])
                      for name, group in (("A", "fitted"), ("A_fail", "excluded"))
                      for sign in ("nonzero", "negative", "positive")}
            state = ("nonzero_fitted_depth" if fields["nonzero_fitted_bins"] else
                     "entirely_masked_depth_design" if fields["nonzero_excluded_bins"] else
                     "zero_depth_design")
        signed = (None if fields["negative_fitted_bins"] is None else
                  bool(fields["negative_fitted_bins"] + fields["negative_excluded_bins"]))
        result[s["structure_id"]] = dict(structure_id=s["structure_id"], feature_index=col,
            kind=s["kind"], source=str(s["source"]), raw_weight=float(s["raw_weight"]),
            weight_N=float(s["weight_N"]), depth_column_state=state,
            source_nclose_counts={_key(k): int(n) for k, n in s["source_nclose_counts"].items()},
            has_signed_depth_entries=signed, **fields)
    return result


def _key(key):
    return ":".join(map(str, key))


def _geometry(endpoints):
    if len(endpoints) != 2:
        raise ValueError("A retained primitive must have two endpoints")
    result = []
    for chrom, boundary, side in endpoints:
        if not isinstance(chrom, str) or not chrom or side not in {"L", "R"}:
            raise ValueError(f"Malformed retained primitive endpoint: {(chrom, boundary, side)!r}")
        # Unencodable/out-of-reference boundaries remain in the candidate set.
        result.append((chrom, _integer(boundary, "Reference boundary"), side))
    return tuple(sorted(result))


def _source_occurrence(key, event, junction, nodes, censat_names):
    a, b = [nodes[i] for i in junction["node_pair"]]
    direct = _geometry(((a[5], int(a[8]) if a[4] == "+" else int(a[7]), "L" if a[4] == "+" else "R"),
                        (b[5], int(b[7]) if b[4] == "+" else int(b[8]), "R" if b[4] == "+" else "L")))
    exact = direct == _geometry(junction["endpoints"]) and junction["mode"] == "ADJACENT_ALIGNMENT"
    query_gap = int(b[2]) - int(a[3]) if exact else None
    same_flow = a[5] == b[5] and a[4] == b[4]
    reference_gap = ((int(b[7]) - int(a[8]) if a[4] == "+" else int(a[7]) - int(b[8]))
                     if same_flow and exact else None)
    names = set(event.get("contig_names", ())) | {event.get("contig_name", "")}
    return dict(source_nclose_key=_key(key), contig_names=sorted(n for n in names if n),
        CEN_source=None if censat_names is None else bool(names & censat_names),
        junction_index=junction["junction_index"], node_pair=list(junction["node_pair"]),
        geometry_mode=junction["mode"], endpoints_are_direct_alignment_boundaries=exact,
        query_gap_bp=query_gap, signed_reference_gap_bp=reference_gap,
        reference_minus_query_gap_bp=None if reference_gap is None else reference_gap-query_gap,
        query_intervals=[[int(a[2]), int(a[3])], [int(b[2]), int(b[3])]],
        reference_intervals=[[a[5], int(a[7]), int(a[8]), a[4]],
                             [b[5], int(b[7]), int(b[8]), b[4]]])


def extract_candidates(prefix, context, nodes, matrix_path=None):
    """Return all exact primitives, their carriers, and current-file provenance."""
    prefix, inputs = Path(prefix), _Inputs()
    implementation = inputs.add(Path(__file__))
    if inputs.records[str(implementation)]["sha256"] != _IMPLEMENTATION_SHA256:
        raise ValueError("Candidate extractor implementation changed since import")
    model, weights, saved_nodes = _load_completed(prefix, inputs, context, nodes)
    observed, matrix_info = _matrix_counts(prefix, prefix / "matrix.h5" if matrix_path is None else matrix_path,
                                            model, weights, inputs)
    features = _feature_catalog(model, observed)
    cen_path = inputs.add(prefix / "censat_endpoint_candidates.tsv", required=False)
    censat_names = None
    if cen_path.exists():
        with cen_path.open() as handle:
            censat_names = {r["unitig"] for r in csv.DictReader(handle, delimiter="\t")
                            if r["status"] == "realigned_consistent"}
    groups = {}

    def add(endpoints, key, kind, nclose_id, occurrence=None):
        geometry = _geometry(endpoints)
        group = groups.setdefault(geometry, dict(source_counts=Counter(), ids=set(), kinds=set(), occurrences=[]))
        group["source_counts"][_key(key)] += 1
        group["ids"].add(nclose_id)
        group["kinds"].add(kind)
        if occurrence is not None:
            group["occurrences"].append(occurrence)

    for key, junctions in nclose_bnd_memberships(context, saved_nodes).items():
        event = model["nclose_sources"][key]
        for junction in junctions:
            add(junction["endpoints"], key, "primitive_BND", junction["nclose_id"],
                _source_occurrence(key, event, junction, saved_nodes, censat_names))
    for key, event in model["nclose_sources"].items():
        if event["kind"] != "indel":
            continue
        kind = event["indel_kind"]
        if kind not in {"deletion", "insertion", "duplication"}:
            raise ValueError(f"Unsupported retained symbolic event kind: {kind}")
        sides = ("L", "R") if kind == "deletion" else ("R", "L")
        canonical = model["nclose_aliases"][key]
        add(((event["chrom"], event["st"], sides[0]), (event["chrom"], event["nd"], sides[1])),
            key, "symbolic_" + kind, model["ncloses"][canonical]["nclose_id"])
    source_carriers = defaultdict(list)
    for sid, feature in features.items():
        for key, count in feature["source_nclose_counts"].items():
            source_carriers[key].append((sid, count))
    rows = []
    for geometry, group in sorted(groups.items()):
        carriers = Counter()
        for key, primitive_count in group["source_counts"].items():
            for sid, source_count in source_carriers[key]:
                carriers[sid] += primitive_count * source_count
        contribution = fsum(features[sid]["weight_N"] * count for sid, count in carriers.items())
        if not isfinite(contribution):
            raise ValueError("Primitive multiplicity-weighted contribution is nonfinite")
        state = ("candidate_no_modeled_structure" if not carriers else
                 "candidate_NNLS_zero" if contribution <= 0 else
                 "candidate_positive_below_export_threshold" if contribution <= .1 else
                 "candidate_positive_above_export_threshold")
        states = Counter(features[sid]["depth_column_state"] for sid in carriers)
        signed_values = [features[sid]["has_signed_depth_entries"] for sid in carriers]
        any_signed = True if any(v is True for v in signed_values) else None if None in signed_values else False
        rows.append(dict(geometry_key=geometry, endpoints=[list(ep) for ep in geometry],
            source_nclose_keys=sorted(group["source_counts"]), nclose_ids=sorted(group["ids"]),
            source_kinds=sorted(group["kinds"]), source_primitive_counts=dict(sorted(group["source_counts"].items())),
            source_alignment_occurrences=sorted(group["occurrences"], key=lambda r: (r["source_nclose_key"], r["junction_index"])),
            carrier_primitive_counts=dict(sorted(carriers.items())), carrier_count=len(carriers),
            conditional_model_contribution_N=contribution, original_candidate_state=state,
            carrier_depth_state_counts=dict(states), any_carrier_signed_depth_entries=any_signed,
            reference_endpoint_distance_bp=abs(geometry[0][1]-geometry[1][1]) if geometry[0][0] == geometry[1][0] else None,
            interpret_reference_endpoint_distance_as_variant_size=False,
            dosage_identifiability="not_assessed"))
    return dict(candidates=rows, features=features, provenance=dict(
        schema="SKYPE.native_local_candidates.v1", coordinate_convention="zero-based half-open boundaries with retained sides L/R",
        selection="All exact retained BND primitives and symbolic adjacencies; no weight, truth, support or normal selection.",
        source_signature_verified=True, source_signature=model["source_signature"], n_unit=float(model["n_unit"]),
        matrix=matrix_info, input_files=inputs.finish(), implementation_sha256=_IMPLEMENTATION_SHA256,
        input_stability_verified=True, historical_generation_binding="unavailable",
        assumptions=["N is the original source/primitive-multiplicity-weighted conditional contribution, not absolute CN or VAF.",
                     "A nonzero whole-carrier column does not establish local junction observability.",
                     "Signed design entries are retained; only fitted NNLS coefficients are constrained nonnegative.",
                     "Symbolic insertion/duplication adjacency does not establish inserted sequence identity."]))
