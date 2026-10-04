"""Native PATH terminal evidence, separate from fitted structure dosage.

Stage 21 records its actual component traversal before depth-row omission.
Stage 31 only joins that record to source sequence and final structure weights.
No function here changes a graph, PAF, matrix, coefficient, VCF or BED.
"""
from __future__ import annotations

import argparse
import ast
from collections import Counter, defaultdict
import csv
import hashlib
import json
from math import fsum
from pathlib import Path
import pickle
import re

from terminal_repeat import (RULES, RULES_PATH, RULES_SHA256, IMPLEMENTATION_SHA256 as REPEAT_IMPLEMENTATION_SHA256,
                             describe_source, repeat_tracts, reverse_complement)

IMPLEMENTATION_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

COMPONENT_FILE = "terminal_component_context.json"
INPUT_FILE = "terminal_source_inputs.json"
REPORT_FILE = "terminal_evidence.json"
USAGE_FILE = "structure_terminal_usage.tsv"
INVENTORY_FILE = "terminal_host_inventory.tsv"
COMPONENT_SCHEMA = "SKYPE.native_terminal_components.v1"
BINDING_SCHEMA = "SKYPE.assembly_alignment_source.v1"
REPORT_SCHEMA = "SKYPE.native_terminal_evidence.v1"


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _node_digest(node):
    return _json_digest([v.rstrip("\r\n") if isinstance(v, str) else v for v in node])


def _write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def _read_pickle(prefix, name):
    with (Path(prefix) / name).open("rb") as handle:
        return pickle.load(handle)


def _write_tsv(path, columns, rows):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "." if row.get(key) is None else row[key] for key in columns})


def depth_row_record(row):
    if row is None:
        return None
    fields = row.split("\t")
    return dict(query=fields[0], query_start=int(fields[2]), query_end=int(fields[3]),
                strand=fields[4], chrom=fields[5], start=int(fields[7]), end=int(fields[8]),
                row_sha256=hashlib.sha256(row.encode()).hexdigest())


def capture_component(component_id, key, raw_path, nodes, rendered_terminals):
    """Called on the exact resolver output, including un-emitted raw endpoints."""
    record = dict(component_id=int(component_id), component_key=key,
                  status="resolved" if raw_path else "empty_raw_component", endpoints={})
    for end, raw in (("start", raw_path[:1]), ("end", raw_path[-1:])):
        if not raw:
            record["endpoints"][end] = None
            continue
        traversal, index = raw[0]
        emitted = rendered_terminals.get(end)
        record["endpoints"][end] = dict(
            node_index=int(index), raw_query_traversal=int(traversal),
            node_sha256=_node_digest(nodes[index]),
            host_depth_row_emitted=emitted is not None,
            emitted_depth_row=depth_row_record(emitted))
    return record


def component_inputs_snapshot(prefix, ppc_paf):
    """Input bytes that must stay fixed while stage 21 resolves components."""
    prefix = Path(prefix)
    dependencies = {name: sha256(prefix / name) for name in
                    ("path_data.pkl", "paf_file_path.pkl")}
    selected = _read_pickle(prefix, "paf_file_path.pkl")
    return dict(ppc_paf=dict(name=Path(ppc_paf).name, sha256=sha256(ppc_paf)),
                dependencies=dependencies,
                selected_pafs=[dict(path=str(Path(path).resolve()), sha256=sha256(path)) for path in selected])


def begin_component_capture(prefix, ppc_paf):
    """A failed rebuild must not leave an apparently complete older capture."""
    (Path(prefix) / COMPONENT_FILE).unlink(missing_ok=True)
    return component_inputs_snapshot(prefix, ppc_paf)


def write_component_context(prefix, records, ppc_paf, *, expected_inputs):
    """Bind a newly captured sidecar to the inputs actually used by its resolver."""
    prefix = Path(prefix)
    current = component_inputs_snapshot(prefix, ppc_paf)
    if current != expected_inputs:
        raise ValueError("Terminal component inputs changed during stage 21; capture was not published")
    current["dependencies"]["contig_pat_vec_data.pkl"] = sha256(prefix / "contig_pat_vec_data.pkl")
    value = dict(schema=COMPONENT_SCHEMA, capture="actual_stage21_component_resolver", **current,
                 input_stability_verified=True,
                 components={str(row["component_id"]): row for row in sorted(records, key=lambda row: row["component_id"])})
    _write_json(prefix / COMPONENT_FILE, value)
    return value


def load_component_context(prefix):
    prefix = Path(prefix)
    path = prefix / COMPONENT_FILE
    if not path.exists():
        return {}, "missing_stage21_context"
    try:
        value = json.loads(path.read_text())
        if value.get("schema") != COMPONENT_SCHEMA:
            return {}, "unsupported_stage21_context"
        if value.get("input_stability_verified") is not True:
            return {}, "unverified_stage21_input_stability"
        if any(sha256(prefix / name) != digest for name, digest in value["dependencies"].items()):
            return {}, "stale_stage21_context"
        ppc = value["ppc_paf"]
        if sha256(prefix / ppc["name"]) != ppc["sha256"]:
            return {}, "stale_stage21_nodes"
        if any(sha256(row["path"]) != row["sha256"] for row in value["selected_pafs"]):
            return {}, "stale_stage21_selected_alignment"
        return value["components"], "captured_stage21_context"
    except (OSError, ValueError, TypeError, KeyError):
        return {}, "unreadable_stage21_context"


def _physical_endpoint(node, into):
    if into not in (0, 1) or node[4] not in ("+", "-"):
        return None
    plus = (node[4] == "+") == (into == 1)
    return [node[5], int(node[7] if plus else node[8]), "R" if plus else "L"]


def project_source_boundary(piece, boundary):
    """All endpoint projections; an insertion boundary is explicitly ambiguous."""
    if not piece["start"] <= boundary <= piece["end"]:
        return [], "outside_source_alignment"
    cs = piece.get("cs")
    if not cs:
        if boundary in (piece["start"], piece["end"]):
            q = (piece["qstart"] if (boundary == piece["start"]) == (piece["strand"] == "+") else piece["qend"])
            return [q], "source_alignment_endpoint"
        return [], "missing_cs_projection"
    tokens = re.findall(r":\d+|=[A-Za-z]+|\*[A-Za-z]{2}|[+-][A-Za-z]+|~[A-Za-z]{2}\d+[A-Za-z]{2}", cs)
    if "".join(tokens) != cs:
        return [], "unsupported_cs_projection"
    reference, query = piece["start"], 0
    candidates, inside_gap = set(), False
    for token in tokens:
        op = token[0]
        length = int(token[1:]) if op == ":" else (1 if op == "*" else
                 int(token[3:-2]) if op == "~" else len(token) - 1)
        dr, dq = (0, length) if op == "+" else ((length, 0) if op in ("-", "~") else (length, length))
        if reference <= boundary <= reference + dr:
            if dr == 0:
                candidates.update((query, query + dq))
            elif dq == 0:
                candidates.add(query)
                inside_gap |= reference < boundary < reference + dr
            else:
                candidates.add(query + boundary - reference)
        reference += dr
        query += dq
    if reference != piece["end"] or query != piece["qend"] - piece["qstart"]:
        return [], "cs_span_mismatch"
    points = sorted(piece["qstart"] + q if piece["strand"] == "+" else piece["qend"] - q for q in candidates)
    status = "inside_reference_gap" if inside_gap else "ambiguous_source_boundary" if len(points) != 1 else "resolved"
    return points, status


def _alias_frame_matches(node, piece):
    """The explicit source-row pointer must also prove the original query frame.

    No alias name is parsed or stripped. Local cut coordinates and guessed
    offsets cannot establish this contract; they remain unverified/unavailable.
    """
    if int(node[1]) != piece["qlen"]:
        return False
    for reference, query in ((int(node[7]), int(node[2] if node[4] == "+" else node[3])),
                             (int(node[8]), int(node[3] if node[4] == "+" else node[2]))):
        candidates, status = project_source_boundary(piece, reference)
        if status not in ("resolved", "source_alignment_endpoint") or candidates != [query]:
            return False
    return True


def _read_selected(path, names, source_indices=()):
    # A preprocessing query alias still points at an original selected PAF row.
    # Discover its original name first, then retain that name's complete chain.
    names = set(names)
    if source_indices:
        with Path(path).open("rb") as handle:
            for index, line in enumerate(handle):
                if index in source_indices:
                    names.add(line.split(b"\t", 1)[0].decode())
    digest, by_index, chains = hashlib.sha256(), {}, defaultdict(list)
    before = Path(path).stat()
    with Path(path).open("rb") as handle:
        for index, line in enumerate(handle):
            digest.update(line)
            if line.split(b"\t", 1)[0].decode() not in names:
                continue
            fields = line.decode().rstrip().split("\t")
            tags = {part[:2]: part[5:] for part in fields[12:] if len(part) >= 5 and part[2] == ":"}
            row = dict(query=fields[0], qlen=int(fields[1]), qstart=int(fields[2]), qend=int(fields[3]),
                       strand=fields[4], chrom=fields[5], target_length=int(fields[6]),
                       start=int(fields[7]), end=int(fields[8]), mapq=int(fields[11]), row_index=index,
                       cs=tags.get("cs"), xi=tags.get("xi"), row_sha256=hashlib.sha256(line).hexdigest())
            by_index[index] = row
            chains[row["query"]].append(row)
    for chain in chains.values():
        chain.sort(key=lambda row: (row["qstart"], row["qend"], row["row_index"]))
    after = Path(path).stat()
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (
            after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns):
        raise ValueError("Selected PAF changed during evidence reading")
    return dict(path=str(Path(path).resolve()), sha256=digest.hexdigest(), rows=by_index, chains=chains)


def read_fasta(path, names=None, regions=None):
    """Read selected records/regions and hash bytes without using/writing an index."""
    path = Path(path)
    before = path.stat()
    digest, lengths, duplicates = hashlib.sha256(), {}, set()
    chunks = defaultdict(list)
    by_chrom = defaultdict(list)
    for key, (chrom, start, end) in (regions or {}).items():
        by_chrom[chrom].append((key, start, end))
    for intervals in by_chrom.values():
        intervals.sort(key=lambda item: (item[1], item[2], item[0]))
    name, position, intervals, first_interval, selected_name = None, 0, [], 0, False
    with path.open("rb") as handle:
        for line in handle:
            digest.update(line)
            if line.startswith(b">"):
                if name is not None:
                    lengths[name] = position
                name = line[1:].split()[0].decode()
                if name in lengths:
                    duplicates.add(name)
                position = 0
                intervals = by_chrom.get(name, [])
                first_interval = 0
                selected_name = names is not None and name in names
                continue
            if name is None:
                if line.strip():
                    raise ValueError("Sequence before FASTA header")
                continue
            bases = line.strip()
            line_end = position + len(bases)
            while first_interval < len(intervals) and intervals[first_interval][2] <= position:
                first_interval += 1
            sequence = bases.decode().upper() if (selected_name or first_interval < len(intervals) and intervals[first_interval][1] < line_end) else None
            if selected_name:
                chunks[name].append(sequence)
            if sequence is not None:
                for key, start, end in intervals[first_interval:]:
                    if start >= line_end:
                        break
                    lo, hi = max(start, position), min(end, line_end)
                    if lo < hi:
                        chunks[key].append(sequence[lo - position:hi - position])
            position = line_end
    if name is not None:
        lengths[name] = position
    after = path.stat()
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (
            after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns):
        raise ValueError("FASTA changed during evidence reading")
    return dict(path=str(path.resolve()), sha256=digest.hexdigest(), lengths=lengths,
                duplicate_names=sorted(duplicates), sequences={key: "".join(parts) for key, parts in chunks.items()})


def validate_upstream_binding(binding_path, source, reference, raw_paf, selected):
    """Validate a producer attestation; never create or upgrade one downstream.

    The supported producer records before/after FASTA/reference hashes during
    fresh raw AND alternate alignment generation. The existing alignasm cache
    then connects those outputs to the selected PAF. Its cache alone is not a
    source attestation. Hashes computed here only check the attested content.
    """
    if not binding_path:
        return dict(status="unverified_no_upstream_source_binding")
    try:
        binding_path = Path(binding_path)
        data = json.loads(binding_path.read_text())
        result = dict(binding_path=str(binding_path.resolve()), binding_sha256=sha256(binding_path))
        if (data.get("schema") != BINDING_SCHEMA or data.get("status") != "complete" or
                data.get("producer_stage") != "raw_and_alternate_alignment_generation" or
                not data.get("generation_commands")):
            return dict(result, status="unsupported_upstream_source_binding")
        expected = dict(fasta={key: source[key] for key in ("path", "sha256")},
                        reference={key: reference[key] for key in ("path", "sha256")})
        if data.get("inputs_before") != expected or data.get("inputs_after") != expected:
            return dict(result, status="source_or_reference_content_mismatch")
        index_binding = data.get("reference_index_binding", {})
        if (index_binding.get("schema") != "SKYPE.reference_index_source.v1" or
                index_binding.get("status") != "complete" or
                index_binding.get("reference") != expected["reference"] or
                index_binding.get("preset") != "asm20"):
            return dict(result, status="unverified_reference_index_source_binding")
        index = index_binding["index"]
        mapper = index_binding["minimap2"]
        if (str(Path(index["path"]).resolve()) != index["path"] or sha256(index["path"]) != index["sha256"] or
                not mapper.get("version") or not mapper.get("sha256") or
                str(Path(mapper["path"]).resolve()) != mapper["path"]):
            return dict(result, status="reference_index_content_mismatch")
        generation = index_binding.get("generation_command", [])
        commands = data["generation_commands"]
        generated = index_binding.get("generation_output", {})
        if (not generation or generation[0] != mapper["path"] or generation.count("-d") != 1 or
                generation[generation.index("-d") + 1] != generated.get("path") or
                generated.get("sha256") != index["sha256"] or
                expected["reference"]["path"] not in generation or
                any(not isinstance(command, list) or not command for command in commands)):
            return dict(result, status="reference_index_command_binding_mismatch")
        outputs = data["outputs"]
        if data.get("outputs_at_generation") != outputs:
            return dict(result, status="alignment_output_generation_binding_mismatch")
        if not raw_paf or outputs["primary_paf"]["path"] != str(Path(raw_paf).resolve()):
            return dict(result, status="raw_paf_identity_mismatch")
        for key in ("primary_paf", "alternate_paf"):
            signature = outputs[key]
            if str(Path(signature["path"]).resolve()) != signature["path"] or sha256(signature["path"]) != signature["sha256"]:
                return dict(result, status="raw_or_alternate_content_mismatch")
        mapping_commands = [command for command in commands if command[0] == mapper["path"]]
        auxiliary_commands = [command for command in commands if command[0] != mapper["path"]]
        if (not mapping_commands or any(index["path"] not in command or command.count("-o") != 1 for command in mapping_commands)):
            return dict(result, status="mapping_command_source_binding_mismatch")
        mapped_outputs = {command[command.index("-o") + 1]: command for command in mapping_commands}
        primary_command = mapped_outputs.get(outputs["primary_paf"]["path"], [])
        if expected["fasta"]["path"] not in primary_command:
            return dict(result, status="mapping_command_source_binding_mismatch")
        if auxiliary_commands:
            extractor = data.get("gap_extractor", {})
            if (len(auxiliary_commands) != 1 or len(auxiliary_commands[0]) != 5 or
                    auxiliary_commands[0][1:4] != [extractor.get("path"), expected["fasta"]["path"], outputs["primary_paf"]["path"]] or
                    sha256(extractor["path"]) != extractor["sha256"]):
                return dict(result, status="gap_extraction_source_binding_mismatch")
            alternate_command = mapped_outputs.get(outputs["alternate_paf"]["path"], [])
            if auxiliary_commands[0][4] not in alternate_command:
                return dict(result, status="alternate_mapping_source_binding_mismatch")
        elif (outputs["alternate_paf"]["sha256"] != hashlib.sha256(b"").hexdigest() or len(mapping_commands) != 1):
            return dict(result, status="missing_alternate_generation_binding")
        link_path = Path(selected["path"] + ".alignasm.json")
        link = json.loads(link_path.read_text())
        if (link.get("output") != {key: selected[key] for key in ("path", "sha256")} or
                any(link["inputs"].get(key) != outputs[key] for key in ("primary_paf", "alternate_paf"))):
            return dict(result, status="selected_paf_content_binding_mismatch")
        return dict(result, status="verified_upstream_content_binding", selected_binding_path=str(link_path),
                    selected_binding_sha256=sha256(link_path), producer_stage=data["producer_stage"])
    except (OSError, ValueError, TypeError, KeyError, IndexError):
        return dict(status="unavailable_or_invalid_upstream_source_binding", binding_path=str(binding_path))


def display_class(port_status, evidence_state, source_context, reference_status, *, anchor_qualified=True):
    """Deterministic unknown precedence; repeat topology remains a separate axis."""
    if (port_status != "resolved" or not anchor_qualified or evidence_state != "assembly_evidence_bound" or
            source_context in ("no_outward_sequence_observed", "repeat_context_ambiguous", "repeat_context_outward_unqualified", "unavailable")):
        return "unavailable_sequence_or_ambiguous_physical_port"
    if source_context == "repeat_reached_through_ordered_aligned_source_pieces":
        return "donor_route_repeat"
    if source_context == "direct_repeat_extension" and reference_status == "no_qualifying_outward_reference_tract":
        return "observed_nonreference_telomeric_extension"
    if source_context == "direct_repeat_extension" and reference_status in (
            "reference_unavailable", "reference_ambiguous", "reference_comparison_limit", "reference_coordinate_mismatch"):
        return "unavailable_sequence_or_ambiguous_physical_port"
    return "ordinary_graph_anchor"


def _source_request(prefix, source_fasta, reference_fasta, raw_paf, source_binding):
    prefix = Path(prefix)
    supplied = dict(source_fasta=source_fasta, reference_fasta=reference_fasta,
                    raw_paf=raw_paf, source_binding=source_binding)
    if any(value is not None for value in supplied.values()):
        request = {key: str(Path(value).resolve()) if value else None for key, value in supplied.items()}
        if not request["source_binding"] and raw_paf:
            candidate = Path(str(raw_paf) + ".source_binding.json")
            if candidate.exists():
                request["source_binding"] = str(candidate.resolve())
        request["meaning"] = "Requested files only; this downstream record does not attest historical alignment generation."
        _write_json(prefix / INPUT_FILE, request)
        return request
    try:
        return json.loads((prefix / INPUT_FILE).read_text())
    except (OSError, ValueError):
        return supplied


def write_terminal_evidence(prefix, context, nodes, *, source_fasta=None,
                            reference_fasta=None, raw_paf=None, source_binding=None):
    """Write all native PATH ends and a separate complete graph-host inventory."""
    prefix = Path(prefix)
    components, component_status = load_component_context(prefix)
    paths = _read_pickle(prefix, "path_data.pkl")
    vectors = _read_pickle(prefix, "contig_pat_vec_data.pkl")
    by_source = {source: ids for source, ids in vectors[0]}
    graph_hosts = {}
    for line in (prefix / "telomere_connected_list.txt").read_text().splitlines():
        label, encoded = line.split("\t", 1)
        port = ast.literal_eval(encoded)
        graph_hosts[(label, int(port[1]))] = dict(terminal_label=label, node_index=int(port[1]),
                                                direct_graph_port=list(port),
                                                direct_graph_endpoint=_physical_endpoint(nodes[int(port[1])], port[0]))
    usages, contexts = [], {}
    for structure in context.structures.values():
        if structure["kind"] != "PATH":
            continue
        source = structure["source"]
        loc = Path(source)
        path = paths[loc.parent.name][int(loc.stem) - 1][0]
        ids = by_source[source]
        counted = Counter()
        for end, path_node, label, component_id in (("start", path[1], path[0][0], ids[0]),
                                                   ("end", path[-2], path[-1][0], ids[-1])):
            index = int(path_node[1])
            counted[index] += 1
            node = nodes[index]
            captured = components.get(str(component_id), {}).get("endpoints", {}).get(end)
            port_status, into, endpoint = component_status, None, None
            if captured is not None:
                if captured["node_index"] != index or captured["node_sha256"] != _node_digest(node):
                    port_status = "terminal_component_node_mismatch"
                elif captured["raw_query_traversal"] not in (0, 1):
                    port_status = "ambiguous_raw_query_traversal"
                else:
                    into = captured["raw_query_traversal"] if end == "start" else 1 - captured["raw_query_traversal"]
                    endpoint = _physical_endpoint(node, into)
                    port_status = "resolved" if endpoint is not None else "ambiguous_alignment_strand"
            key = (index, into, tuple(endpoint or ()), port_status)
            if key not in contexts:
                context_id = f"SKYPE.TERMINAL_CONTEXT.{len(contexts) + 1}"
                contexts[key] = dict(context_id=context_id, node_index=index,
                    physical_port_status=port_status, physical_endpoint=endpoint,
                    entry_query_traversal=into, query_toward_terminal=None if into is None else (-1 if into == 1 else 1),
                    node_alignment_strand=node[4], node_reference_interval=[node[5], int(node[7]), int(node[8])],
                    source_query=node[0], source_row_identity=str(node[21]).strip(),
                    assembly_evidence_state="unavailable", source_repeat_context="unavailable",
                    reference_comparison_status="reference_unavailable",
                    internal_repeat_with_farther_aligned_flank=None)
            row = contexts[key]
            usages.append(dict(structure_id=structure["structure_id"], feature_index=structure["feature_index"],
                terminal_end=end, terminal_label=label, node_index=index, context_id=row["context_id"],
                path_port=":".join(map(str, path_node[:2])), terminal_component_id=int(component_id),
                conditional_structure_weight=float(structure["raw_weight"]),
                conditional_structure_weight_N=float(structure["weight_N"]),
                host_depth_row_emitted=None if captured is None else captured["host_depth_row_emitted"],
                emitted_depth_row=None if captured is None else captured["emitted_depth_row"]))
        if counted != structure["telomere_counts"]:
            raise ValueError(f"Terminal occurrence mismatch for {structure['structure_id']}")
    request = _source_request(prefix, source_fasta, reference_fasta, raw_paf, source_binding)
    wanted = {nodes[index][0] for _, index in graph_hosts} | {row["source_query"] for row in contexts.values()}
    selected, source_sequences, source_failure = {}, None, None
    requested_indices = defaultdict(set)
    for index in {index for _, index in graph_hosts} | {row["node_index"] for row in contexts.values()}:
        try:
            file_index, source_index = map(int, str(nodes[index][21]).strip().split("."))
            if file_index < 2:
                requested_indices[file_index].add(source_index)
        except ValueError:
            pass
    for index, path in enumerate(_read_pickle(prefix, "paf_file_path.pkl")):
        try:
            selected[index] = _read_selected(path, wanted, requested_indices[index])
            wanted.update(selected[index]["chains"])
        except (OSError, ValueError, IndexError):
            selected[index] = dict(path=str(path), status="selected_alignment_unavailable")
    try:
        if request.get("source_fasta"):
            source_sequences = read_fasta(request["source_fasta"], names=wanted)
        else:
            source_failure = "source_fasta_not_supplied"
    except (OSError, ValueError, UnicodeError) as error:
        source_failure = type(error).__name__ + ": " + str(error)
    reference_regions, tract_cache = {}, {}
    for row in contexts.values():
        node = nodes[row["node_index"]]
        try:
            file_index, source_index = map(int, row["source_row_identity"].split("."))
        except ValueError:
            row["assembly_evidence_state"] = "invalid_source_row_identity"
            continue
        if file_index < 0 or source_index < 0:
            row["assembly_evidence_state"] = "invalid_source_row_identity"
            continue
        if file_index >= 2:
            row["assembly_evidence_state"] = "inapplicable_virtual_reference_node"
            continue
        paf = selected.get(file_index, {})
        piece = paf.get("rows", {}).get(source_index)
        if piece is None or piece["chrom"] != node[5] or piece["strand"] != node[4]:
            row["assembly_evidence_state"] = "source_alignment_identity_mismatch"
            row["source_pointer_candidate"] = piece
            continue
        if piece["query"] != node[0] and not _alias_frame_matches(node, piece):
            row["assembly_evidence_state"] = "source_query_frame_mismatch"
            row["source_pointer_candidate"] = piece
            continue
        query_name = piece["query"]
        row.update(source_paf_index=file_index, selected_paf_path=paf["path"], selected_paf_sha256=paf["sha256"],
                   source_alignment=piece, selected_source_chain=paf["chains"][query_name],
                   graph_query_label=node[0], source_query=query_name, source_query_alias_resolved=query_name != node[0],
                   alias_resolution="explicit_selected_PAF_row_identity_and_CIGAR_query_frame" if query_name != node[0] else None)
        if source_sequences is None:
            row["source_unavailable_reason"] = source_failure
            continue
        sequence = source_sequences["sequences"].get(query_name)
        if sequence is None or query_name in source_sequences["duplicate_names"] or len(sequence) != piece["qlen"]:
            row["assembly_evidence_state"] = "source_query_missing_duplicate_or_length_mismatch"
            continue
        row.update(assembly_evidence_state="assembly_evidence_unverified", source_sequence_length=len(sequence),
                   current_source_sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest())
        if row["physical_port_status"] != "resolved":
            continue
        projections, projection_status = project_source_boundary(piece, row["physical_endpoint"][1])
        row.update(source_boundary_candidates=projections, source_boundary_status=projection_status)
        if len(projections) != 1 or projection_status not in ("resolved", "source_alignment_endpoint"):
            row["physical_port_status"] = "ambiguous_source_projection"
            continue
        boundary, toward = projections[0], row["query_toward_terminal"]
        row["source_query_boundary"] = boundary
        if query_name not in tract_cache:
            tract_cache[query_name] = repeat_tracts(sequence)
        description = describe_source(sequence, boundary, toward, paf["chains"][query_name], piece, tract_cache[query_name])
        row.update(source_evidence=description, source_repeat_context=description["context"],
                   internal_repeat_with_farther_aligned_flank=description["internal_repeat_with_farther_aligned_flank"])
        retained_query = piece["qend"] - boundary if toward == -1 else boundary - piece["qstart"]
        chrom, position, side = row["physical_endpoint"]
        retained_reference = piece["end"] - position if side == "R" else position - piece["start"]
        row["retained_source_anchor"] = dict(query_bp=retained_query, reference_bp=retained_reference, mapq=piece["mapq"])
        row["source_anchor_qualified"] = (piece["mapq"] >= RULES["minimum_host_mapq"] and
            min(retained_query, retained_reference) >= RULES["minimum_host_query_and_reference_span_bp"])
        distances = [abs((t["tract"]["end"] if toward == 1 else t["tract"]["start"]) - boundary)
                     for t in description["adjacent_tracts"] if t["tract"]["primary"]]
        length = max(distances + [0]) + RULES["reference_padding_bp"]
        if length > RULES["maximum_reference_comparison_bp"]:
            row["reference_comparison_status"] = "reference_comparison_limit"
            continue
        start, end = (max(0, position - length), position) if side == "R" else (position, min(int(node[6]), position + length))
        row["reference_requested_interval"] = [chrom, start, end]
        row["reference_requested_outward_length"] = length
        row["reference_available_outward_length"] = end - start
        reference_regions[row["context_id"]] = (chrom, start, end)
    reference_sequences, reference_failure = None, None
    try:
        if request.get("reference_fasta"):
            reference_sequences = read_fasta(request["reference_fasta"], regions=reference_regions)
        else:
            reference_failure = "reference_fasta_not_supplied"
    except (OSError, ValueError, UnicodeError) as error:
        reference_failure = type(error).__name__ + ": " + str(error)
    bindings = {}
    if source_sequences is not None and reference_sequences is not None:
        for index, paf in selected.items():
            bindings[index] = (validate_upstream_binding(request.get("source_binding"), source_sequences,
                reference_sequences, request.get("raw_paf"), paf) if "sha256" in paf else
                dict(status="selected_alignment_unavailable"))
    for row in contexts.values():
        binding = bindings.get(row.get("source_paf_index"), dict(status="source_or_reference_unavailable"))
        row["source_binding_status"] = binding["status"]
        if row["assembly_evidence_state"] == "assembly_evidence_unverified" and binding["status"] == "verified_upstream_content_binding":
            row["assembly_evidence_state"] = "assembly_evidence_bound"
        if row["context_id"] in reference_regions and reference_sequences is not None:
            chrom, start, end = reference_regions[row["context_id"]]
            length = reference_sequences["lengths"].get(chrom)
            sequence = reference_sequences["sequences"].get(row["context_id"], "")
            if length is None or chrom in reference_sequences["duplicate_names"] or not 0 <= start <= end <= length or len(sequence) != end - start:
                status = "reference_coordinate_mismatch"
            elif start == end:
                status = "reference_end_no_outward_sequence"
            elif end - start < row["reference_requested_outward_length"]:
                status = "reference_end_truncated_comparison"
            elif any(base not in "ACGT" for base in sequence):
                status = "reference_ambiguous"
            else:
                if row["physical_endpoint"][2] == "R":
                    sequence = reverse_complement(sequence)
                tracts = repeat_tracts(sequence)
                row["reference_repeat_tracts"] = tracts
                status = "qualifying_outward_reference_repeat" if any(t["primary"] for t in tracts) else "no_qualifying_outward_reference_tract"
            row.update(reference_comparison_status=status, reference_interval_sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest())
        row["display_class"] = display_class(row["physical_port_status"], row["assembly_evidence_state"],
                                               row["source_repeat_context"], row["reference_comparison_status"],
                                               anchor_qualified=row.get("source_anchor_qualified") is True)
        row["evidence_origin"] = "assembly_sequence_and_selected_alignment"
        row["functional_chromosome_cap_inferred"] = False
    usage_by_host = defaultdict(list)
    for row in usages:
        usage_by_host[(row["terminal_label"], row["node_index"])].append(row)
        emitted = row.pop("emitted_depth_row")
        row["emitted_depth_row"] = json.dumps(emitted, separators=(",", ":")) if emitted is not None else None
    inventory = []
    for key, host in sorted(graph_hosts.items()):
        rows = usage_by_host.get(key, [])
        node = nodes[host["node_index"]]
        inventory.append(dict(terminal_label=key[0], node_index=key[1], source_query=node[0],
            graph_only=not rows, modeled_occurrences=len(rows), positive_weight_occurrences=sum(r["conditional_structure_weight"] > 0 for r in rows),
            conditional_occurrence_weight_N=fsum(r["conditional_structure_weight_N"] for r in rows),
            direct_graph_port=json.dumps(host["direct_graph_port"], separators=(",", ":")),
            direct_graph_endpoint=json.dumps(host["direct_graph_endpoint"], separators=(",", ":")),
            modeled_context_ids=";".join(sorted({r["context_id"] for r in rows})),
            source_sequence_applicability="inapplicable_virtual_reference_node" if str(node[21]).strip().split(".")[0] not in ("0", "1") else "assembly_source",
            inventory_status="graph_only_no_modeled_physical_use" if not rows else "modeled_path_terminal"))
    if set(usage_by_host) - set(graph_hosts):
        raise ValueError("Modeled terminal use is missing from graph host inventory")
    result = dict(schema=REPORT_SCHEMA, scope="all_native_PATH_terminal_occurrences_and_graph_host_inventory",
        source_evidence_origin="assembly-derived; independent molecule corroboration is not supplied",
        component_context_status=component_status, source_request=request,
        rules_sha256=RULES_SHA256, rules=RULES,
        implementation_sha256={"terminal_evidence.py": IMPLEMENTATION_SHA256, "terminal_repeat.py": REPEAT_IMPLEMENTATION_SHA256},
        source_files={"fasta": {k: source_sequences[k] for k in ("path", "sha256")} if source_sequences else None,
                      "reference": {k: reference_sequences[k] for k in ("path", "sha256")} if reference_sequences else None},
        source_unavailable_reason=source_failure, reference_unavailable_reason=reference_failure,
        upstream_bindings={str(k): v for k, v in bindings.items()},
        counts=dict(path_structures=len(usages) // 2, terminal_occurrences=len(usages), graph_hosts=len(inventory),
                    graph_only_hosts=sum(r["graph_only"] for r in inventory), physical_source_contexts=len(contexts),
                    display_classes=dict(Counter(r["display_class"] for r in contexts.values()))),
        contexts=list(contexts.values()),
        limitations=[
            "All repeat descriptions are assembly-derived and may be unverified against historical alignment inputs.",
            "An assembly end cannot distinguish a chromosome end from an internal or expanded repeat with an unobserved distal flank.",
            "An internal-repeat/farther-flank flag remains part of every repeat context, including a nonreference extension display.",
            "A donor-route label preserves an ordered selected-alignment path; it does not establish donor identity or independent-molecule support.",
            "No cap function, chromosome healing, somatic origin, clonality, absolute copy number or complete chromosome path is inferred.",
            "Occurrence weights are conditional structure-fit contributions, not an identified telomere dose or a count of chromosome ends.",
            "Graph-only direct edges are inventory metadata and are not assigned a modeled physical use.",
            "Full-assembly-only and VCF-only workflows are outside the native PATH annotation contract."])
    _write_json(prefix / REPORT_FILE, result)
    _write_tsv(prefix / USAGE_FILE, ["structure_id", "feature_index", "terminal_end", "terminal_label", "node_index", "context_id", "path_port", "terminal_component_id", "conditional_structure_weight", "conditional_structure_weight_N", "host_depth_row_emitted", "emitted_depth_row"], usages)
    _write_tsv(prefix / INVENTORY_FILE, ["terminal_label", "node_index", "source_query", "graph_only", "modeled_occurrences", "positive_weight_occurrences", "conditional_occurrence_weight_N", "direct_graph_port", "direct_graph_endpoint", "modeled_context_ids", "source_sequence_applicability", "inventory_status"], inventory)
    return result


def main():
    parser = argparse.ArgumentParser(description="Annotate a completed native model without changing its fitted artifacts.")
    parser.add_argument("prefix")
    for option in ("source-fasta", "reference-fasta", "raw-paf", "source-binding"):
        parser.add_argument("--" + option)
    args = parser.parse_args()
    from types import SimpleNamespace
    model = _read_pickle(args.prefix, "structure_nclose_model.pkl")
    nodes = _read_pickle(args.prefix, "01_nclose_data.pkl")["contig_data"]
    context = SimpleNamespace(structures={s["structure_id"]: s for s in model["structures"]})
    result = write_terminal_evidence(args.prefix, context, nodes, source_fasta=args.source_fasta,
        reference_fasta=args.reference_fasta, raw_paf=args.raw_paf, source_binding=args.source_binding)
    print(json.dumps(result["counts"]))


if __name__ == "__main__":
    main()
