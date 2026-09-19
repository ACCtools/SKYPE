"""Native structure weights and original NClose multiplicities.

Stage 22 writes topology, stage 31 joins the fitted coefficients and qualified
virtual inversions. Renderers never feed an aggregated NClose back into a
structure coefficient. The legacy event catalog remains preprocessing history.
"""
from __future__ import annotations

import csv
import hashlib
import pickle
from collections import Counter, defaultdict
from math import fsum
from pathlib import Path

from native_bnd import constituent_pairs, layout_endpoints, node_pair_endpoints
from nclose_tracking import (
    REPORT_COLUMNS, _build_bnd_event, compressed_bnd_event_keys,
    load_type4_edge_event_map, nclose_event_id_by_key,
)

MODEL_FILE = "structure_nclose_model.pkl"
MODEL_VERSION = 2


def read_pickle(prefix, name):
    with (Path(prefix) / name).open("rb") as handle:
        return pickle.load(handle)


def source_signature(prefix):
    digest = hashlib.sha256()
    for name in ("tot_loc_list.pkl", "path_data.pkl", "contig_pat_vec_data.pkl",
                 "nclose_event_catalog.pkl", "ecdna_circuit_data.pkl",
                 "conjoined_type4_ins_del.pkl", "type4_indel_graph_edges.pkl",
                 "01_nclose_data.pkl", "telomere_connected_list.txt"):
        path = Path(prefix) / name
        digest.update(name.encode())
        if path.exists():
            digest.update(path.read_bytes())
    return digest.hexdigest()


def key_text(key):
    return ":".join(map(str, key))


def initialize_nclose_registry(model):
    """Keep source pairs, but reuse one identity per exact native BND geometry."""
    if "nclose_sources" in model:
        return
    events = list(model["ncloses"].values())
    model.update(ncloses={}, nclose_sources={}, nclose_aliases={}, bnd_by_endpoints={})
    # Prefer the existing compressed-graph identity over restored AMP sources.
    for event in sorted(events, key=lambda e: e.get("source") != "compressed_nclose_nodes_list.txt"):
        register_nclose(model, event)


def register_nclose(model, event):
    initialize_nclose_registry(model)
    key = event["event_key"]
    if key in model["nclose_sources"]:
        return model["nclose_aliases"][key]
    geometry = tuple(sorted(event["endpoints"])) if event["kind"] == "bnd" else None
    canonical = model["bnd_by_endpoints"].get(geometry, key) if geometry is not None else key
    event = dict(event)
    if canonical == key:
        if "nclose_id" not in event:
            event["nclose_id"] = f"SKYPE.nclose.{model['next_nclose_id']}"
            model["next_nclose_id"] += 1
        model["ncloses"][key] = event
        if geometry is not None:
            model["bnd_by_endpoints"][geometry] = key
    model["nclose_sources"][key] = event
    model["nclose_aliases"][key] = canonical
    return canonical


def set_structure_membership(model, structure, source_counts):
    """Combine identities after source-position counting; preserve real repeats."""
    structure["source_nclose_counts"] = Counter(source_counts)
    counts = Counter()
    for key, count in source_counts.items():
        counts[model["nclose_aliases"][key]] += count
    structure["nclose_counts"] = counts


def count_path_members(path, bnd_keys, edge_events, constituents, aliases=None):
    """Count traversals, resolving aliases at the same traversal position.

    Distinct positions add, even when they visit the same original NClose.
    A direct tag and a Type4 tag at one position describe one traversal.
    """
    at_position = defaultdict(Counter)
    s = 1
    while s < len(path) - 2:
        a, b = path[s][1], path[s + 1][1]
        if isinstance(a, int) and isinstance(b, int):
            key = tuple(sorted((a, b)))
            if key in bnd_keys:
                at_position[s][key] = 1
                s += 2
                continue
        s += 1
    for pos, (a, b) in enumerate(zip(path, path[1:])):
        event = edge_events.get((tuple(a[:2]), tuple(b[:2])))
        if event is not None:
            if aliases is None:
                at_position[pos] |= constituents[event]
            else:
                # Direct and compound descriptions of one traversal can use
                # different source pairs for the same exact NClose. Take their
                # maximum multiplicity at this position, not their sum.
                existing, incoming = Counter(), Counter()
                for key, count in at_position[pos].items():
                    existing[aliases[key]] += count
                for key, count in constituents[event].items():
                    incoming[aliases[key]] += count
                for key, count in constituents[event].items():
                    canonical = aliases[key]
                    extra = min(max(0, count - at_position[pos][key]),
                                max(0, incoming[canonical] - existing[canonical]))
                    if extra:
                        at_position[pos][key] += extra
                        existing[canonical] += extra
    total = Counter()
    for counts in at_position.values():
        total.update(counts)
    return total


def add_bnd_metadata(model, key, nodes, source, layout=None):
    initialize_nclose_registry(model)
    key = tuple(sorted(key))
    if key in model["nclose_sources"]:
        return model["nclose_aliases"][key]
    event = _build_bnd_event(key, nodes, source)
    endpoints, names = (layout_endpoints(layout) if layout is not None
                        else node_pair_endpoints(nodes, key))
    event.update(endpoints=endpoints, contig_names=names)
    return register_nclose(model, event)


def build_structure_model(prefix, nodes=None):
    """Same builder for a new stage-22 model and a legacy stage-31 replay."""
    catalog = read_pickle(prefix, "nclose_event_catalog.pkl")
    locations = read_pickle(prefix, "tot_loc_list.pkl")
    paths = read_pickle(prefix, "contig_pat_vec_data.pkl")[0]
    path_dict = read_pickle(prefix, "path_data.pkl")
    circuits = read_pickle(prefix, "ecdna_circuit_data.pkl")[0]
    ins, dels = read_pickle(prefix, "conjoined_type4_ins_del.pkl")
    merged_circuits = list(ins) + list(dels)
    if nodes is None:
        nodes = read_pickle(prefix, "01_nclose_data.pkl")["contig_data"]
    ids = nclose_event_id_by_key(catalog)
    model = dict(version=MODEL_VERSION, source_signature=source_signature(prefix),
                 ncloses={}, structures=[], path_splits=[],
                 next_nclose_id=len(catalog) + 1, column_count=len(locations))
    constituents = {}
    event_by_location = {}
    for original in sorted(catalog, key=lambda e: e.get("source") != "compressed_nclose_nodes_list.txt"):
        event = dict(original)
        key = event["event_key"]
        if event["kind"] == "indel":
            event_by_location[(event["event_type"], event["event_idx"])] = event
        if event["kind"] == "indel" and event["type2_merge_idx"] >= 0:
            constituents[key] = Counter(constituent_pairs(
                merged_circuits[event["type2_merge_idx"] - 1]))
            continue
        event["nclose_id"] = ids[key]
        if event["kind"] == "bnd":
            event["endpoints"], event["contig_names"] = node_pair_endpoints(nodes, key)
        register_nclose(model, event)
        constituents[key] = Counter({key: 1})
    # IDs for all original catalog events are reserved, including compound
    # aliases that are now structure metadata rather than actual NClose rows.
    for counts in constituents.values():
        for key in counts:
            if len(key) == 2:
                add_bnd_metadata(model, key, nodes, "MERGE_TYPE4")
    for circuit in circuits:
        for key in constituent_pairs(circuit):
            add_bnd_metadata(model, key, nodes, "AMP")
    edge_events = load_type4_edge_event_map(prefix, catalog)
    bnd_keys = compressed_bnd_event_keys(catalog)
    for column, location in enumerate(locations):
        loc = Path(location)
        structure = dict(structure_id=f"SKYPE.STRUCTURE.{column + 1}",
                         feature_index=column, source=str(location),
                         nclose_counts=Counter(), telomere_counts=Counter(),
                         legacy_ids=[])
        if column < len(paths):
            path_loc = Path(paths[column][0])
            path = path_dict[path_loc.parent.name][int(path_loc.stem) - 1][0]
            structure["kind"] = "PATH"
            structure["nclose_counts"] = count_path_members(
                path, bnd_keys, edge_events, constituents, model["nclose_aliases"])
            # Actual path ends identify telomere-bearing graph nodes.
            structure["telomere_counts"].update([path[1][1], path[-2][1]])
        elif loc.parent.name in {"front_jump", "back_jump"}:
            event = event_by_location[(loc.parent.name, int(loc.stem.split('_')[0]))]
            key = event["event_key"]
            structure.update(kind="MERGE_TYPE4" if event["type2_merge_idx"] >= 0 else "TYPE4",
                             nclose_counts=constituents[key].copy(),
                             event=dict(event), legacy_ids=[ids[key]])
        elif loc.parent.name == "ecdna":
            index = int(loc.stem) - 1
            circuit = circuits[index]
            structure.update(kind="AMP", circuit_index=index,
                             legacy_ids=[f"SKYPE.AMP.{index + 1}"],
                             nclose_counts=Counter(constituent_pairs(circuit)),
                             span=(nodes[circuit[0]][5], min(nodes[i][7] for i in circuit),
                                   max(nodes[i][8] for i in circuit)))
        elif loc.parent.parent.name == "12_cent_fragment":
            structure.update(kind="CENTROMERE", chrom=loc.parent.name,
                             side=loc.stem)
        else:
            raise ValueError(f"Unknown native matrix feature: {location}")
        set_structure_membership(model, structure, structure["nclose_counts"])
        model["structures"].append(structure)
    return model


def save_structure_model(prefix, model):
    with (Path(prefix) / MODEL_FILE).open("wb") as handle:
        pickle.dump(model, handle)


def load_or_build_structure_model(prefix, nodes=None):
    path = Path(prefix) / MODEL_FILE
    if path.exists():
        model = read_pickle(prefix, MODEL_FILE)
        if (model.get("version") == MODEL_VERSION
                and model.get("source_signature") == source_signature(prefix)):
            # Virtual qualification is re-evaluated against the final weights.
            model["structures"] = [s for s in model["structures"] if s["kind"] != "VIRTUAL_INV"]
            return model
    return build_structure_model(prefix, nodes)


def add_virtual_structure(model, record, raw_weight, views, nodes):
    pair_id = int(record["pair_id"])
    counts = Counter()
    for side in ("a", "b"):
        key = tuple(sorted(record[f"nclose_key_{side}"]))
        add_bnd_metadata(model, key, nodes, "VIRTUAL_INV", record[f"layout_{side}"])
        counts[key] += 1
    structure = dict(
        structure_id=f"SKYPE.STRUCTURE.VIRTUAL_INV.{pair_id}", kind="VIRTUAL_INV",
        feature_index=None, source=f"raw_translocation_result.pkl:{pair_id}",
        legacy_ids=[f"RAW_TRANSLOCATION_PAIR_{pair_id}"], nclose_counts=counts,
        telomere_counts=Counter(), raw_weight=float(raw_weight), views=views)
    set_structure_membership(model, structure, counts)
    model["structures"].append(structure)


def set_path_splits(model, feature_usage, parent_usage, idx_to_key):
    initialize_nclose_registry(model)
    splits = []
    for number, (key, counts) in enumerate(sorted(feature_usage.items()), 1):
        parent, split_index, ca, pa, da, cb, pb, db, name = key
        parent_key = idx_to_key[parent]
        canonical = model["nclose_aliases"][parent_key]
        splits.append(dict(
            projection_id=f"{model['ncloses'][canonical]['nclose_id']}.PATH_SPLIT.{number}",
            parent_key=parent_key, split_index=split_index, contig_names=(name,),
            endpoints=((ca, int(pa), "L" if da == "+" else "R"),
                       (cb, int(pb), "R" if db == "+" else "L")),
            column_counts=dict(counts)))
    model["path_splits"] = splits
    model["split_parent_usage"] = {idx_to_key[k]: dict(v) for k, v in parent_usage.items()}


class StructureWeights:
    def __init__(self, model, weights, n_unit):
        initialize_nclose_registry(model)
        if len(weights) != model["column_count"]:
            raise ValueError("Structure model and fitted column count differ; rebuild stage 22")
        self.model, self.n_unit = model, float(n_unit)
        self.ncloses = model["ncloses"]
        self.structures = {s["structure_id"]: s for s in model["structures"]}
        self.contributions = []
        self.by_nclose = defaultdict(list)
        self.telomere_weights = defaultdict(float)
        for structure in self.structures.values():
            set_structure_membership(model, structure,
                                     structure.get("source_nclose_counts", structure["nclose_counts"]))
            col = structure["feature_index"]
            if col is not None:
                structure["raw_weight"] = float(weights[col])
            structure["weight_N"] = structure["raw_weight"] / self.n_unit
            for key, count in structure["nclose_counts"].items():
                sources = {source: n for source, n in structure["source_nclose_counts"].items()
                           if model["nclose_aliases"][source] == key}
                row = self.contribution(structure, key, count, source_counts=sources)
                self.contributions.append(row)
                self.by_nclose[key].append(row)
            for node, count in structure["telomere_counts"].items():
                self.telomere_weights[node] += count * structure["raw_weight"]
        self.totals = {key: fsum(r["contribution"] for r in self.by_nclose[key])
                       for key in self.ncloses}
        self.projections, self.projected_rows = self._project()
        self.projected_totals = {
            key: fsum(r["contribution"] for r in rows)
            for key, rows in self.projected_rows.items()}
        model["n_unit"] = self.n_unit

    def contribution(self, structure, key, count, projection_id=None, source_counts=None):
        source_counts = {key: count} if source_counts is None else source_counts
        key = self.model["nclose_aliases"][key]
        sources = sorted(source_counts)
        return dict(structure_id=structure["structure_id"], kind=structure["kind"],
                    feature_index=structure["feature_index"], nclose_key=key,
                    nclose_id=self.ncloses[key]["nclose_id"], occurrence_count=int(count),
                    structure_weight=structure["raw_weight"],
                    structure_weight_N=structure["weight_N"],
                    contribution=int(count) * structure["raw_weight"],
                    contribution_N=int(count) * structure["weight_N"],
                    projection_id=projection_id or self.ncloses[key]["nclose_id"],
                    source_nclose_keys=tuple(sources),
                    source_occurrence_counts=tuple(int(source_counts[k]) for k in sources))

    def _project(self):
        metadata = {e["nclose_id"]: dict(e, projection_id=e["nclose_id"],
                                       parent_key=key, is_split=False)
                    for key, e in self.ncloses.items()}
        rows = defaultdict(list)
        removed = self.model.get("split_parent_usage", {})
        for structure in self.structures.values():
            remaining = defaultdict(Counter)
            for source, original_count in structure["source_nclose_counts"].items():
                count = original_count - removed.get(source, {}).get(structure["feature_index"], 0)
                if count < 0:
                    raise ValueError(f"PATH_SPLIT exceeds source NClose count: {source}")
                if count:
                    remaining[self.model["nclose_aliases"][source]][source] = count
            for key, sources in remaining.items():
                pid = self.ncloses[key]["nclose_id"]
                rows[pid].append(self.contribution(structure, key, sum(sources.values()),
                                                   source_counts=sources))
        by_column = {s["feature_index"]: s for s in self.structures.values()
                     if s["feature_index"] is not None}
        for split in self.model.get("path_splits", []):
            pid = split["projection_id"]
            canonical = self.model["nclose_aliases"][split["parent_key"]]
            metadata[pid] = dict(split, kind="bnd", is_split=True,
                                 nclose_id=self.ncloses[canonical]["nclose_id"])
            for col, count in split["column_counts"].items():
                rows[pid].append(self.contribution(by_column[col], split["parent_key"], count, pid))
        return metadata, rows

    def visible_projections(self, min_cn=0.1):
        return [(self.projections[pid], total) for pid, total in self.projected_totals.items()
                if total / self.n_unit > min_cn]

    def structure_summaries(self, min_cn=0.1):
        return [s for s in self.structures.values()
                if s["kind"] in {"MERGE_TYPE4", "AMP", "VIRTUAL_INV"}
                and s["weight_N"] > min_cn]

    def cn_lists(self, min_cn=0.1):
        inv, trans, indel = [], [], []
        for key, event in self.ncloses.items():
            cn = self.totals[key] / self.n_unit
            if cn <= min_cn:
                continue
            if event["kind"] == "indel":
                indel.append(cn)
            else:
                a, b = event["endpoints"]
                if a[0] != b[0]:
                    trans.append(cn)
                elif a[2] == b[2]:
                    inv.append(cn)
        return inv, trans, indel


def write_tsv(path, columns, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_structure_reports(prefix, context, status):
    structures = []
    for s in context.structures.values():
        structures.append(dict(
            structure_id=s["structure_id"], kind=s["kind"],
            feature_index=s["feature_index"] if s["feature_index"] is not None else ".",
            raw_weight=s["raw_weight"], weight_N=s["weight_N"], source=s["source"],
            legacy_ids=";".join(s["legacy_ids"]),
            nclose_occurrence_count=sum(s["nclose_counts"].values())))
    write_tsv(Path(prefix) / "structure_report.tsv",
              ["structure_id", "kind", "feature_index", "raw_weight", "weight_N",
               "source", "legacy_ids", "nclose_occurrence_count"], structures)
    usage = []
    for row in context.contributions:
        usage.append({k: (key_text(v) if k == "nclose_key" else
                         ";".join(key_text(key) for key in v) if k == "source_nclose_keys" else
                         ";".join(map(str, v)) if k == "source_occurrence_counts" else
                         "." if v is None else v)
                      for k, v in row.items() if k != "projection_id"})
    write_tsv(Path(prefix) / "structure_nclose_usage.tsv",
              ["structure_id", "kind", "feature_index", "nclose_key", "nclose_id",
               "occurrence_count", "structure_weight", "structure_weight_N",
               "contribution", "contribution_N", "source_nclose_keys", "source_occurrence_counts"], usage)
    sources = []
    for key, event in context.model["nclose_sources"].items():
        canonical = context.model["nclose_aliases"][key]
        sources.append(dict(source_nclose_key=key_text(key),
                            source_nclose_id=event.get("nclose_id", "."),
                            nclose_id=context.ncloses[canonical]["nclose_id"],
                            canonical_nclose_key=key_text(canonical), is_alias=int(key != canonical),
                            contig_names=";".join(event.get("contig_names", ())),
                            source=event.get("source", ".")))
    write_tsv(Path(prefix) / "nclose_sources.tsv",
              ["source_nclose_key", "source_nclose_id", "nclose_id", "canonical_nclose_key",
               "is_alias", "contig_names", "source"], sources)
    reports = []
    history = status.get("stages", {}).get("base", status.get("stages", {}).get("initial", {})).get("reasons", {})
    for key, event in context.ncloses.items():
        rows = context.by_nclose[key]
        cn = context.totals[key] / context.n_unit
        row = {k: event[k] for k in REPORT_COLUMNS[:7]}
        row.update(nclose_cn=cn, nclose_cn_reason="SUPPORTED" if cn > 0 else "ZERO_SUPPORT",
                   nclose_filter="NA", nclose_filter_reason="NOT_RUN_RAW_NNLS",
                   nclose_cluster="NA", nclose_cluster_reason="NOT_RUN_RAW_NNLS",
                   preprocessing_reason=history.get(key, "PASS" if event.get("source") not in {"AMP", "MERGE_TYPE4", "VIRTUAL_INV"} else "RESTORED_CONSTITUENT"),
                   kind=event["kind"], nclose_key=key_text(key),
                   model_cn=fsum(r["contribution_N"] for r in rows if r["feature_index"] is not None),
                   virtual_cn=fsum(r["contribution_N"] for r in rows if r["feature_index"] is None))
        reports.append(row)
    write_tsv(Path(prefix) / "nclose_report.tsv",
              list(REPORT_COLUMNS) + ["preprocessing_reason", "kind", "nclose_key", "model_cn", "virtual_cn"], reports)
