"""Opt-in local HiFi evidence, separate from fitted dosage and default calls.

The versioned gate measures local alignment support. It does not establish
somatic origin, a complete chromosome path, or independent validation.
"""
import argparse
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import time
import uuid

import pysam

import local_hifi_decoder as decoder

HERE = Path(__file__).resolve().parent
POLICY_PATH = HERE / "local_hifi_policy.json"
TOLERANCES = (100, 500)


def digest(path):
    value = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 << 20), b""):
            value.update(block)
    return value.hexdigest()


def identity(path):
    path = Path(path).resolve(strict=True)
    stat = path.stat()
    return dict(path=str(path), device=stat.st_dev, inode=stat.st_ino,
                size=stat.st_size, mtime_ns=stat.st_mtime_ns, ctime_ns=stat.st_ctime_ns)


def signature(path):
    before = identity(path)
    value = dict(before, sha256=digest(path))
    if before != identity(path):
        raise ValueError(f"Input changed during checksum: {path}")
    return value


def unchanged(value):
    if identity(value["path"]) != {k: v for k, v in value.items() if k != "sha256"}:
        raise ValueError(f"Input changed during evidence assessment: {value['path']}")


def save(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + ".partial")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def save_gzip(path, value):
    temporary = Path(str(path) + ".partial")
    with gzip.open(temporary, "wt") as handle:
        json.dump(value, handle, allow_nan=False)
        handle.write("\n")
    temporary.replace(path)


def policy_and_sources():
    policy = json.loads(POLICY_PATH.read_text())
    sources = {name: signature(HERE / name) for name in (
        "local_hifi_policy.json", "local_hifi_decoder.py", "local_hifi_evidence.py",
        "native_local_candidates.py", "native_vcf_bnd.py", "structure_nclose.py",
        "pipeline.py", "31_depth_analysis.py")}
    if sources["local_hifi_decoder.py"]["sha256"] != policy["decoder_sha256"]:
        raise ValueError("The local HiFi decoder differs from the frozen evidence policy")
    return policy, sources


def reference_identity(path):
    """Read an uncompressed FASTA once; bind file and named sequence identities."""
    before, full = identity(path), hashlib.sha256()
    sequences, name, length, md5, sha = {}, None, 0, None, None

    def finish():
        if name is not None:
            if name in sequences or not length:
                raise ValueError("Reference has duplicate names or an empty sequence")
            sequences[name] = dict(length=length, md5=md5.hexdigest(), sha256=sha.hexdigest())

    with Path(path).open("rb") as handle:
        for line in handle:
            full.update(line)
            if line.startswith(b">"):
                finish()
                name = line[1:].split()[0].decode("ascii")
                length, md5, sha = 0, hashlib.md5(), hashlib.sha256()
            else:
                seq = b"".join(line.split()).upper()
                if seq and name is None:
                    raise ValueError("Expected an uncompressed reference FASTA")
                if seq:
                    length += len(seq)
                    md5.update(seq)
                    sha.update(seq)
    finish()
    if not sequences or identity(path) != before:
        raise ValueError("Reference is empty or changed while reading")
    sequence_identity = hashlib.sha256(json.dumps(sequences, sort_keys=True,
                                                   separators=(",", ":")).encode()).hexdigest()
    return dict(file=dict(before, sha256=full.hexdigest()), sequences=sequences,
                sequence_identity_sha256=sequence_identity)


def find_bam_index(path):
    path = Path(path)
    for candidate in (Path(str(path) + ".bai"), path.with_suffix(".bai"),
                      Path(str(path) + ".csi"), path.with_suffix(".csi")):
        if candidate.is_file():
            return candidate
    raise ValueError("An existing BAM BAI/CSI index is required")


def check_bam_reference(bam, reference):
    header = bam.header.to_dict()
    if header.get("HD", {}).get("SO") != "coordinate" or not bam.has_index():
        raise ValueError("Evidence BAM must be coordinate sorted and indexed")
    checked, missing_md5 = [], []
    for row in header.get("SQ", []):
        chrom = row["SN"]
        expected = reference["sequences"].get(chrom)
        if expected is None or int(row["LN"]) != expected["length"]:
            raise ValueError(f"BAM/reference name or length mismatch: {chrom}")
        if "M5" in row:
            if row["M5"].lower() != expected["md5"]:
                raise ValueError(f"BAM/reference sequence MD5 mismatch: {chrom}")
            checked.append(chrom)
        else:
            missing_md5.append(chrom)
    if not header.get("SQ"):
        raise ValueError("Evidence BAM has no reference dictionary")
    return dict(sequence_MD5_verified_contigs=checked,
                sequence_MD5_unavailable_contigs=missing_md5,
                missing_BAM_contigs=sorted(set(reference["sequences"])-set(bam.references)),
                historical_alignment_reference_generation="not_attested",
                reference_check="SQ names/lengths and every available M5; missing M5 remains unverified")


def endpoint_status(endpoint, lengths):
    chrom, position, side = endpoint
    if chrom not in lengths:
        return "missing_reference"
    if side not in {"L", "R"} or not 0 <= position-(side == "L") < lengths[chrom]:
        return "unencodable_reference_boundary"
    return "assessable"


def assign_ids(candidates, reference):
    """Identity includes exact oriented geometry and named reference sequence identity."""
    result, seen = [], set()
    for candidate in candidates:
        endpoints = sorted([list(ep) for ep in candidate["endpoints"]])
        payload = [reference["sequence_identity_sha256"], endpoints]
        key = hashlib.sha256(json.dumps(payload, separators=(",", ":")).encode()).hexdigest()
        if key in seen:
            raise ValueError("Duplicate exact candidate geometry")
        seen.add(key)
        result.append(dict(candidate, id="J" + key, endpoints=endpoints))
    return sorted(result, key=lambda row: row["endpoints"])


def associate_model_reference(targets, nodes, reference):
    """Localize uncertain model/reference coordinate joins without renumbering."""
    observed = defaultdict(set)
    for node in nodes:
        observed[str(node[5])].add(int(node[6]))
    mismatches = {}
    for chrom, lengths in observed.items():
        expected = reference["sequences"].get(chrom, {}).get("length")
        if lengths != {expected}:
            mismatches[chrom] = dict(native_node_lengths=sorted(lengths), supplied_reference_length=expected)
    affected = []
    for target in targets:
        states = ["model_reference_length_mismatch" if ep[0] in mismatches else "assessable"
                  for ep in target["endpoints"]]
        target["model_reference_endpoint_assessment"] = states
        if any(state != "assessable" for state in states):
            affected.append(target["id"])
    return dict(check="Current native node chromosome names/lengths against supplied FASTA",
        mismatched_contigs=mismatches, affected_candidate_ids=affected,
        action="Affected endpoint measurements remain null; other candidate endpoints use the unchanged local gate",
        historical_model_reference_generation="not_attested",
        sequence_identity_basis="user_supplied_reference; historical generation is not inferred from matching lengths")


def merged_regions(targets, lengths):
    intervals = defaultdict(list)
    for target in targets:
        for ep in target["endpoints"]:
            if endpoint_status(ep, lengths) == "assessable":
                chrom, position, _ = ep
                intervals[chrom].append((max(0, position-2000), min(lengths[chrom], position+2000)))
    regions = []
    for chrom, values in sorted(intervals.items()):
        merged = []
        for start, end in sorted(values):
            if merged and start <= merged[-1][1]:
                merged[-1][1] = max(end, merged[-1][1])
            else:
                merged.append([start, end])
        regions.extend((chrom, start, end) for start, end in merged)
    return regions


def query_bam(bam_path, index_path, targets, destination):
    seen = set()
    with pysam.AlignmentFile(str(bam_path), "rb", index_filename=str(index_path)) as bam:
        lengths = dict(zip(bam.references, bam.lengths))
        regions = merged_regions(targets, lengths)
        save(destination / "query_regions.json", regions)
        path = destination / "queried_alignments.bam"
        header = bam.header.to_dict()
        header.setdefault("HD", {})["SO"] = "unsorted"
        with pysam.AlignmentFile(str(path), "wb", header=header) as cache:
            for number, (chrom, start, end) in enumerate(regions, 1):
                for record in bam.fetch(chrom, start, end):
                    if record.is_unmapped or record.is_secondary:
                        continue
                    key = (record.query_name, record.reference_id, record.reference_start,
                           record.flag, record.cigarstring)
                    if key not in seen:
                        seen.add(key)
                        cache.write(record)
                if number % 100 == 0:
                    print("LOCAL_HIFI_QUERY", number, "/", len(regions), flush=True)
    return path, dict(query_region_count=len(regions), queried_records=len(seen))


def assess_cache(targets, bam_path):
    """Assess a regional BAM cache using the frozen decoder and distinct names."""
    reads, record_count = defaultdict(list), 0
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        lengths = dict(zip(bam.references, bam.lengths))
        for record in bam:
            if record.is_unmapped or record.is_secondary:
                continue
            record_count += 1
            reads[record.query_name].extend(decoder.read_pieces(record))
    groups, boundaries = defaultdict(list), defaultdict(list)
    states = [[endpoint_status(ep, lengths) for ep in row["endpoints"]] for row in targets]
    for index, target in enumerate(targets):
        model_states = target.get("model_reference_endpoint_assessment", ["assessable", "assessable"])
        if len(model_states) != 2:
            raise ValueError("A candidate must have two model-reference endpoint assessments")
        for endpoint_index, state in enumerate(model_states):
            if state != "assessable":
                states[index][endpoint_index] = state
    for index, target in enumerate(targets):
        a, b = target["endpoints"]
        if all(state == "assessable" for state in states[index]):
            groups[(a[0], a[2], b[0], b[2])].append((a[1], b[1], index))
        for side_index, (chrom, position, side) in enumerate(target["endpoints"]):
            if states[index][side_index] == "assessable":
                boundaries[(chrom, side)].append((position, index, side_index))
    for collection in (groups, boundaries):
        for values in collection.values():
            values.sort()
    positions = {key: [row[0] for row in values] for key, values in groups.items()}
    boundary_positions = {key: [row[0] for row in values] for key, values in boundaries.items()}
    support = {tol: [set() for _ in targets] for tol in TOLERANCES}
    exposure = {tol: [[set(), set()] for _ in targets] for tol in TOLERANCES}
    for name, pieces in reads.items():
        for piece in decoder.normalized_pieces(pieces):
            for side, (start, end) in decoder.flank_boundary_ranges(piece).items():
                group = (piece["chrom"], side)
                coords = boundary_positions.get(group, [])
                for tol in TOLERANCES:
                    lo, hi = bisect_left(coords, start-tol), bisect_right(coords, end+tol)
                    for _, index, side_index in boundaries.get(group, [])[lo:hi]:
                        exposure[tol][index][side_index].add(name)
        for a, b in decoder.connections(pieces):
            group = (a[0], a[2], b[0], b[2])
            coords = positions.get(group, [])
            for tol in TOLERANCES:
                lo, hi = bisect_left(coords, a[1]-tol), bisect_right(coords, a[1]+tol)
                for _, second, index in groups.get(group, [])[lo:hi]:
                    if abs(second-b[1]) <= tol:
                        support[tol][index].add(name)
    evidence, witnesses = [], []
    for index, target in enumerate(targets):
        valid = all(state == "assessable" for state in states[index])
        for tol in TOLERANCES:
            names, exposed = support[tol][index], exposure[tol][index]
            if valid and not all(names <= values for values in exposed):
                raise ValueError("Supporting molecule lacks conditional endpoint exposure")
            counts = [len(values) if state == "assessable" else None
                      for values, state in zip(exposed, states[index])]
            passes = len(names) >= 3 if valid else None
            state = ("UNASSESSABLE_WITH_REASON" if not valid else
                     "LOCAL_HIFI_SUPPORTED" if passes else "LOCAL_HIFI_BELOW_GATE")
            evidence.append(dict(target, tolerance_bp=tol, endpoint_assessment=states[index],
                assay_status="assessable" if valid else "unassessable_reference",
                local_evidence_state=state, distinct_support=len(names) if valid else None,
                support_read_names=sorted(names) if valid else None, passes_local_molecule_gate=passes,
                conditional_endpoint_exposure=counts,
                conditional_min_endpoint_exposure=min(counts) if valid else None,
                exposure_union=len(exposed[0] | exposed[1]) if valid else None,
                exposure_intersection=len(exposed[0] & exposed[1]) if valid else None,
                independent_validation_state="not_performed", normal_evidence=None,
                ONT_evidence=None, whole_chain_evidence=None, somatic_state="not_assessed",
                dosage_identifiability="not_assessed", unique_locus_origin="not_assessed",
                insertion_sequence_identity="not_assessed"))
            witnesses.append(dict(id=target["id"], tolerance_bp=tol,
                endpoint_exposure_read_names=[sorted(values) if state == "assessable" else None
                    for values, state in zip(exposed, states[index])]))
    return evidence, witnesses, dict(cached_records=record_count, queried_molecules=len(reads),
                                    reference_lengths=lengths, support_exposure_subset_violations=0)


def relations(rows):
    """Describe resolution neighbors and shared molecules; never collapse alleles."""
    result = []
    for tol in TOLERANCES:
        selected = [row for row in rows if row["tolerance_bp"] == tol]
        by_id = {row["id"]: row for row in selected}
        pairs, by_name = set(), defaultdict(list)
        for row in selected:
            for name in row["support_read_names"] or []:
                by_name[name].append(row["id"])
        shared = Counter()
        for ids in by_name.values():
            ids = sorted(set(ids))
            for i, left in enumerate(ids):
                for right in ids[i+1:]:
                    shared[(left, right)] += 1
        pairs.update(shared)
        for i, left in enumerate(selected):
            for right in selected[i+1:]:
                if decoder.matches(left["endpoints"], right["endpoints"], 500):
                    pairs.add(tuple(sorted((left["id"], right["id"]))))
        for left_id, right_id in sorted(pairs):
            left, right = by_id[left_id], by_id[right_id]
            shared_assessable = (left["support_read_names"] is not None
                                 and right["support_read_names"] is not None)
            result.append(dict(tolerance_bp=tol, left=left_id, right=right_id,
                oriented_neighbor_100bp=decoder.matches(left["endpoints"], right["endpoints"], 100),
                oriented_neighbor_500bp=decoder.matches(left["endpoints"], right["endpoints"], 500),
                shared_supporting_read_names=shared[(left_id, right_id)] if shared_assessable else None,
                shared_support_assessment="assessable" if shared_assessable else "unassessable",
                shared_support_missing_reason=None if shared_assessable else "At least one candidate has unassessable local support",
                shared_source_nclose_keys=sorted(set(left["source_nclose_keys"]) & set(right["source_nclose_keys"])),
                conditional_model_contributions_N=[left["conditional_model_contribution_N"], right["conditional_model_contribution_N"]],
                interpretation="Resolution/shared-molecule relation, not independent event counting or proof of the same allele"))
    return result


def summary(rows):
    result = {}
    for tol in TOLERANCES:
        selected = [row for row in rows if row["tolerance_bp"] == tol]
        strata = Counter()
        for row in selected:
            hifi = row["passes_local_molecule_gate"]
            model = row["conditional_model_contribution_N"] > .1
            label = ("unassessable" if hifi is None else "both" if hifi and model
                     else "HiFi_only" if hifi else "model_only" if model else "neither")
            strata[label] += 1
        result[str(tol)] = dict(candidate_count=len(selected), strata_counts=dict(strata),
            supported=sum(row["passes_local_molecule_gate"] is True for row in selected),
            supported_zero_model_contribution=sum(row["passes_local_molecule_gate"] is True and
                row["conditional_model_contribution_N"] == 0 for row in selected),
            model_comparison_gate_N=0.1, state_counts=dict(Counter(row["local_evidence_state"] for row in selected)))
    return result


def write_table(path, rows):
    columns = ["id", "endpoints", "tolerance_bp", "local_evidence_state", "endpoint_assessment",
               "distinct_support", "conditional_endpoint_exposure", "conditional_model_contribution_N",
               "original_candidate_state", "carrier_count", "carrier_depth_state_counts",
               "source_nclose_keys", "nclose_ids", "source_kinds", "independent_validation_state",
               "same_input_as_assembly", "normal_evidence", "ONT_evidence", "whole_chain_evidence",
               "dosage_identifiability", "any_carrier_signed_depth_entries"]
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(row.get(key), separators=(",", ":"), allow_nan=False)
                             for key in columns})


def vcf_records(row):
    common = (f"SVTYPE=BND;CANDIDATE_ID={row['id']};HIFI_SR={row['distinct_support']}"
              f";HIFI_TOL={row['tolerance_bp']}"
              f";HIFI_FLANK_EXPOSURE={','.join(map(str, row['conditional_endpoint_exposure']))}"
              f";MODEL_CONTRIB_N={row['conditional_model_contribution_N']:.17g}")
    for i, (local, remote) in enumerate((row["endpoints"], row["endpoints"][::-1]), 1):
        chrom, boundary, side = local
        mate_chrom, mate_boundary, mate_side = remote
        pos, mate_pos = boundary+(side == "R"), mate_boundary+(mate_side == "R")
        bracket = "[" if mate_side == "R" else "]"
        mate = f"{bracket}{mate_chrom}:{mate_pos}{bracket}"
        alt = "N"+mate if side == "L" else mate+"N"
        yield (chrom, pos, f"{row['id']}.{i}", "N", alt, ".", "ExperimentalEvidence",
               common+f";MATEID={row['id']}.{3-i}")


def export_vcf(destination, rows, lengths, reference_id):
    from parse_vcf import parse_vcf_bnd_alt
    rank = {chrom: i for i, chrom in enumerate(lengths)}
    export = {}
    definitions = [
        ("SVTYPE", "1", "String", "Primitive adjacency represented as paired BND"),
        ("MATEID", "1", "String", "Reciprocal BND record ID"),
        ("CANDIDATE_ID", "1", "String", "Exact geometry and reference sequence identity; not an allele ID"),
        ("HIFI_SR", "1", "Integer", "Distinct supporting local HiFi read names"),
        ("HIFI_TOL", "1", "Integer", "Boundary matching tolerance in bp"),
        ("HIFI_FLANK_EXPOSURE", "2", "Integer", "Conditional exposure in canonical endpoint order; not VAF denominator"),
        ("MODEL_CONTRIB_N", "1", "Float", "Multiplicity-weighted conditional model contribution; not absolute CN")]
    for tol in TOLERANCES:
        selected, excluded = [], []
        for row in rows:
            if row["tolerance_bp"] != tol:
                continue
            reason = ([endpoint_status(ep, lengths) for ep in row["endpoints"]]
                      if any(endpoint_status(ep, lengths) != "assessable" for ep in row["endpoints"])
                      else None)
            if row["passes_local_molecule_gate"] is None:
                reason = row.get("endpoint_assessment", ["local_evidence_unassessable"])
            if any(re.search(r"[\s,:;<>\[\]{}]", ep[0]) for ep in row["endpoints"]):
                reason = ["unencodable_VCF_contig_name"]
            if reason:
                excluded.append(dict(id=row["id"], reason=reason))
            elif row["passes_local_molecule_gate"]:
                selected.append(row)
        path = destination / f"ExperimentalEvidence.{tol}bp.vcf"
        header = ["##fileformat=VCFv4.3", "##source=SKYPE_local_HiFi_evidence_v1",
            f"##skype_reference_sequence_identity_sha256={reference_id}",
            "##skype_scope=Local alignment evidence; independent validation and somatic origin not assessed",
            "##skype_resolution=Nearby geometries and shared molecules need not be independent alleles",
            '##FILTER=<ID=ExperimentalEvidence,Description="Experimental local evidence; not a validated callset">']
        header += [f'##contig=<ID={chrom},length={length}>' for chrom, length in lengths.items()
                   if not re.search(r"[\s,:;<>\[\]{}]", chrom)]
        header += [f'##INFO=<ID={name},Number={number},Type={kind},Description="{description}">'
                   for name, number, kind, description in definitions]
        emitted = [record for row in selected for record in vcf_records(row)]
        with path.open("x") as handle:
            handle.write("\n".join(header)+"\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for record in sorted(emitted, key=lambda row: (rank[row[0]], row[1], row[2])):
                handle.write("\t".join(map(str, record))+"\n")
        expected = {row["id"]: row for row in selected}
        mates = {}
        with pysam.VariantFile(str(path)) as handle:
            if list(handle.header.samples):
                raise ValueError("Experimental VCF must not contain sample/genotype columns")
            for record in handle:
                target = expected[record.info["CANDIDATE_ID"]]
                alt = parse_vcf_bnd_alt(record.alts[0])
                local_side = "L" if alt.dir_a == "+" else "R"
                remote_side = "R" if alt.dir_b == "+" else "L"
                decoded = sorted([[record.chrom, record.pos-(local_side == "R"), local_side],
                    [alt.mate_chrom, alt.mate_pos-(remote_side == "R"), remote_side]])
                if decoded != target["endpoints"] or record.qual is not None or list(record.filter) != ["ExperimentalEvidence"]:
                    raise ValueError("Experimental VCF boundary/filter round-trip failed")
                if record.id in mates:
                    raise ValueError("Duplicate experimental BND record ID")
                mates[record.id] = record.info["MATEID"]
        if len(mates) != 2*len(selected) or any(mates.get(mate) != rid for rid, mate in mates.items()):
            raise ValueError("Experimental VCF reciprocal mate validation failed")
        export[str(tol)] = dict(candidate_count=len(selected), record_count=len(mates),
                               sha256=digest(path), unencodable=excluded,
                               reciprocal_boundary_roundtrip=True, PASS_records=0, genotypes=0)
    save(destination / "export.json", export)
    return export


def write_local_hifi_evidence(prefix, context, nodes, *, bam_path, reference_path,
                             matrix_path=None, export=False, same_input_as_assembly=None):
    from native_local_candidates import extract_candidates
    prefix = Path(prefix)
    root = prefix / "local_hifi_evidence"
    root.mkdir(exist_ok=True)
    destination = root / (datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "-" + uuid.uuid4().hex)
    destination.mkdir()
    manifest = dict(schema="SKYPE.local_hifi_evidence.v1", status="running",
                    started_utc=datetime.now(timezone.utc).isoformat(),
                    same_input_as_assembly=same_input_as_assembly,
                    assembly_input_relationship_basis="user_declaration" if same_input_as_assembly is not None else "unknown",
                    independent_validation_state="not_performed")
    save(destination / "manifest.json", manifest)
    started = time.monotonic()
    try:
        policy, sources = policy_and_sources()
        model_inputs = {name: signature(prefix / name) for name in (
            "structure_nclose_model.pkl", "weight.npy", "B.npy", "predict_B.npy",
            "tot_loc_list.pkl", "01_nclose_data.pkl", "path_data.pkl", "contig_pat_vec_data.pkl",
            "nclose_event_catalog.pkl", "ecdna_circuit_data.pkl", "conjoined_type4_ins_del.pkl",
            "type4_indel_graph_edges.pkl", "telomere_connected_list.txt", "stage31_depth_coordinates.tsv")
            if (prefix / name).is_file()}
        selected_matrix = Path(matrix_path) if matrix_path is not None else prefix / "matrix.h5"
        if selected_matrix.is_file():
            model_inputs["matrix"] = signature(selected_matrix)
        reference = reference_identity(reference_path)
        bam_signature, index_path = signature(bam_path), find_bam_index(bam_path)
        index_signature = signature(index_path)
        with pysam.AlignmentFile(str(bam_path), "rb", index_filename=str(index_path)) as bam:
            reference_check = check_bam_reference(bam, reference)
        extracted = extract_candidates(prefix, context, nodes, matrix_path=matrix_path)
        targets = assign_ids(extracted["candidates"], reference)
        model_reference = associate_model_reference(targets, nodes, reference)
        manifest.update(policy=policy, implementation=sources, model_inputs=model_inputs,
                        candidate_extraction=extracted["provenance"], reference=reference,
                        model_reference_association=model_reference,
                        BAM=bam_signature, BAM_index=index_signature, BAM_reference_check=reference_check,
                        candidate_count=len(targets))
        save(destination / "manifest.json", manifest)
        save(destination / "candidates.json", targets)
        save(destination / "features.json", extracted["features"])
        cache, query_stats = query_bam(bam_path, index_path, targets, destination)
        rows, witnesses, stats = assess_cache(targets, cache)
        for row in rows:
            row["same_input_as_assembly"] = same_input_as_assembly
        save(destination / "evidence.json", rows)
        save_gzip(destination / "endpoint_exposure_names.json.gz", witnesses)
        save(destination / "resolution_and_shared_molecules.json", relations(rows))
        write_table(destination / "evidence.tsv", rows)
        if export:
            manifest["export"] = export_vcf(destination, rows, stats["reference_lengths"],
                                             reference["sequence_identity_sha256"])
        for value in [*sources.values(), *model_inputs.values(), reference["file"], bam_signature, index_signature]:
            unchanged(value)
        manifest.update(status="complete", elapsed_seconds=time.monotonic()-started,
                        summary=summary(rows), **query_stats, **stats)
        manifest["outputs"] = {path.name: digest(path) for path in destination.iterdir()
                               if path.is_file() and path.name != "manifest.json"}
        save(destination / "manifest.json", manifest)
        save(root / "latest.json", dict(run=destination.name, manifest_sha256=digest(destination / "manifest.json")))
        return manifest
    except Exception as exc:
        manifest.update(status="failed", error=f"{type(exc).__name__}: {exc}",
                        elapsed_seconds=time.monotonic()-started)
        save(destination / "manifest.json", manifest)
        raise


def add_arguments(parser):
    parser.add_argument("--local-hifi-bam", help="Opt-in local HiFi alignment evidence; coordinate-sorted indexed BAM")
    parser.add_argument("--local-hifi-reference", help="Matching uncompressed reference FASTA")
    parser.add_argument("--local-hifi-matrix", help="Optional retained matrix.h5 for design-state annotation")
    parser.add_argument("--local-hifi-export", action="store_true", help="Also write separate ExperimentalEvidence BND VCFs")
    parser.add_argument("--local-hifi-same-input-as-assembly", choices=("yes", "no", "unknown"), default="unknown")


def validate_arguments(args, native=True):
    enabled = bool(args.local_hifi_bam)
    if bool(args.local_hifi_reference) != enabled:
        raise ValueError("--local-hifi-bam and --local-hifi-reference must be supplied together")
    if not enabled and (args.local_hifi_export or args.local_hifi_matrix or args.local_hifi_same_input_as_assembly != "unknown"):
        raise ValueError("Local HiFi options require --local-hifi-bam and --local-hifi-reference")
    if enabled and not native:
        raise ValueError("Local HiFi evidence is supported for native assembly models only")
    return enabled


def run_from_arguments(args, prefix, context, nodes):
    if not validate_arguments(args, native=context is not None):
        return None
    return write_local_hifi_evidence(prefix, context, nodes, bam_path=args.local_hifi_bam,
        reference_path=args.local_hifi_reference, matrix_path=args.local_hifi_matrix,
        export=args.local_hifi_export,
        same_input_as_assembly={"yes": True, "no": False, "unknown": None}[args.local_hifi_same_input_as_assembly])


def main():
    from native_local_candidates import read_completed_context
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prefix", help="Completed native SKYPE artifact directory")
    add_arguments(parser)
    args = parser.parse_args()
    if not validate_arguments(args):
        parser.error("--local-hifi-bam and --local-hifi-reference are required")
    context, nodes = read_completed_context(args.prefix)
    run_from_arguments(args, args.prefix, context, nodes)


if __name__ == "__main__":
    main()
