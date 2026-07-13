"""Parse caller VCF records into graph-independent SKYPE events.

The module deliberately stops at ``NCloseSpec`` and type4 event dictionaries.
Synthetic PAF nodes and graph construction remain in
``02_Build_Breakend_Graph_Limited.py``.

For Severus releases older than 1.7, paired BND ALT fields emitted for mixed
``STRANDS`` values were assigned to the wrong endpoints.  Such inputs are
assumed to be raw, uncorrected Severus output.  ``STRANDS`` is therefore
required for their BND records and is used to normalize directions in memory;
the input VCF is never rewritten.

GRIPSS single breakends without a remote locus are reported and skipped, while
paired records use standard BND ALT semantics.  Unique MATEID links are
authoritative so Manta microhomology shifts do not duplicate one junction as
two singleton events.  VCF contigs absent from a primary-only reference FAI
are handled event-by-event instead of rejecting the whole file.
"""

from __future__ import annotations

import logging
import re
from collections import defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Collection, Mapping, Sequence

from skype_utils import VCF_TYPE4_MIN_SPAN, invert_vcf_strand


DEFAULT_PASS_FILTERS = ("PASS", ".")
SUPPORTED_CALLERS = (
    "sniffles2",
    "severus",
    "nanomonsv",
    "savana",
    "gripss",
    "manta",
    "svaba",
)
SOURCE_PREFIX_CALLERS = (
    ("generatesvcandidates", "manta"),
    *((caller, caller) for caller in SUPPORTED_CALLERS),
)
VERSION_HEADER_CALLERS = {
    "gripssversion": "gripss",
}
VERSION_HEADER_NAMES = {
    "gripssversion": "gripssVersion",
}
VALID_SEVERUS_STRANDS = {"++", "+-", "-+", "--"}
SEVERUS_FIXED_VERSION = (1, 7)


class VcfParseError(ValueError):
    """Raised when the VCF as a whole cannot be parsed safely."""


@dataclass(frozen=True)
class SourceMetadata:
    source_lines: tuple[str, ...]
    caller: str | None
    version: str | None
    version_parts: tuple[int, ...] | None


@dataclass(frozen=True)
class BreakendAlt:
    mate_chrom: str
    mate_pos: int
    dir_a: str
    dir_b: str
    shape: str


@dataclass
class VcfRecord:
    row_index: int
    line_no: int
    chrom: str = ""
    pos: int = 0
    record_id: str = "."
    explicit_id: bool = False
    ref: str = ""
    alt: str = ""
    qual: str = "."
    filt: str = "."
    info: dict[str, str | bool] | None = None
    raw: str = ""
    malformed_reason: str | None = None

    @property
    def svtype(self) -> str:
        if self.info is None:
            return ""
        value = self.info.get("SVTYPE", "")
        return str(value).upper() if value is not True else ""


@dataclass(frozen=True)
class NCloseSpec:
    event_name: str
    record_ids: tuple[str, ...]
    line_nos: tuple[int, ...]
    svtype: str
    chrom_a: str
    pos_a: int
    dir_a: str
    chrom_b: str
    pos_b: int
    dir_b: str


@dataclass
class VcfParseResult:
    source: SourceMetadata
    nclose_specs: list[NCloseSpec]
    type4_events: list[dict]
    summary: dict[str, object]
    skipped_records: list[tuple[object, ...]]
    orientation_mismatches: list[tuple[object, ...]]


@dataclass(frozen=True)
class _NormalizedBreakend:
    alt: BreakendAlt | None
    error: str | None = None
    severus_corrected: bool = False


def sanitize_vcf_id(value: object) -> str:
    value = str(value) if value not in (None, "") else "unknown"
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def parse_vcf_info(info_text: str) -> dict[str, str | bool]:
    info: dict[str, str | bool] = {}
    if info_text in {"", "."}:
        return info
    for item in info_text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True
    return info


def parse_optional_int(value: str | bool | None) -> int | None:
    if value in {None, True, ""}:
        return None
    try:
        return int(float(str(value).split(",", 1)[0]))
    except ValueError:
        return None


def parse_vcf_bnd_alt(alt: str | None) -> BreakendAlt | None:
    """Parse the first explicit BND ALT without validating its local anchor."""

    if not alt:
        return None
    alt = str(alt).split(",", 1)[0]
    if alt in {".", "<BND>"}:
        return None
    brackets = [(idx, char) for idx, char in enumerate(alt) if char in "[]"]
    if len(brackets) != 2:
        return None
    (first_idx, first_bracket), (second_idx, second_bracket) = brackets
    if first_bracket != second_bracket:
        return None
    mate = alt[first_idx + 1 : second_idx]
    if ":" not in mate:
        return None
    chrom, pos_text = mate.rsplit(":", 1)
    if not chrom:
        return None
    try:
        pos = int(pos_text)
    except ValueError:
        return None
    if pos < 1:
        return None
    return BreakendAlt(
        mate_chrom=chrom,
        mate_pos=pos,
        dir_a="+" if first_idx > 0 else "-",
        dir_b="+" if first_bracket == "[" else "-",
        shape=("prefix" if first_idx > 0 else "suffix") + first_bracket,
    )


def is_vcf_single_breakend_alt(alt: str | None) -> bool:
    """Return whether ALT is the VCF single-breakend ``A.``/``.A`` form."""

    if not alt:
        return False
    allele = str(alt).split(",", 1)[0]
    if any(bracket in allele for bracket in "[]") or allele.count(".") != 1:
        return False
    return len(allele) > 1 and (allele.startswith(".") ^ allele.endswith("."))


def vcf_ins_query_name(record: VcfRecord) -> str:
    return f"vcf_ins_{record.line_no}_{sanitize_vcf_id(record.record_id)}"


def select_vcf_type4_graph_events(
    type4_events: Sequence[Mapping[str, object]],
) -> list[Mapping[str, object]]:
    """Return step-11 VCF events eligible for ``--add_indel_graph``.

    Sequence-resolved INS records need their ALT-derived PAF rows in step 11,
    but must not create graph rescue edges. DEL, DUP, and indel-like BND
    events retain their normalized parser order and remain eligible.
    """

    return [
        event
        for event in type4_events
        if str(event.get("svtype", "")).upper() in {"DEL", "DUP", "BND"}
    ]


def _parse_contig_header(line: str) -> tuple[str, int] | None:
    if not line.startswith("##contig=<") or not line.endswith(">"):
        return None
    fields: dict[str, str] = {}
    for item in line[len("##contig=<") : -1].split(","):
        if "=" in item:
            key, value = item.split("=", 1)
            fields[key] = value
    if "ID" not in fields or "length" not in fields:
        return None
    try:
        return fields["ID"], int(fields["length"])
    except ValueError:
        return None


def _read_vcf(
    vcf_path: str | Path,
) -> tuple[dict[str, int], list[VcfRecord], list[str], dict[str, list[str]]]:
    path = Path(vcf_path)
    if not path.is_file():
        raise VcfParseError(f"VCF input does not exist or is not a file: {path}")
    if path.suffix == ".gz":
        raise VcfParseError("SKYPE VCF mode requires a plain-text, uncompressed VCF")

    header_contigs: dict[str, int] = {}
    records: list[VcfRecord] = []
    source_lines: list[str] = []
    version_headers: dict[str, list[str]] = defaultdict(list)
    column_header_count = 0
    saw_record = False
    try:
        handle = path.open("rt", encoding="utf-8")
    except OSError as exc:
        raise VcfParseError(f"cannot open VCF {path}: {exc}") from exc

    try:
        with handle:
            for line_no, raw_line in enumerate(handle, start=1):
                line = raw_line.rstrip("\r\n")
                if not line:
                    continue
                if line.startswith("#"):
                    if saw_record:
                        raise VcfParseError(
                            f"line {line_no}: VCF header occurs after records"
                        )
                    if line.startswith("##") and "=" in line:
                        key, value = line[2:].split("=", 1)
                        key_lower = key.lower()
                        if key_lower == "source":
                            source_lines.append(value)
                        if key_lower in VERSION_HEADER_CALLERS:
                            version_headers[key_lower].append(value)
                    if line.startswith("##contig=<"):
                        contig = _parse_contig_header(line)
                        if contig is not None:
                            header_contigs[contig[0]] = contig[1]
                    if line.startswith("#CHROM\t"):
                        column_header_count += 1
                    continue

                saw_record = True
                row_index = len(records)
                columns = line.split("\t")
                if len(columns) < 8:
                    records.append(
                        VcfRecord(
                            row_index=row_index,
                            line_no=line_no,
                            raw=line,
                            malformed_reason="fewer_than_8_columns",
                        )
                    )
                    continue

                chrom, pos_text, raw_id, ref, alt, qual, filt, info_text = columns[:8]
                try:
                    pos = int(pos_text)
                except ValueError:
                    records.append(
                        VcfRecord(
                            row_index=row_index,
                            line_no=line_no,
                            raw=line,
                            malformed_reason="bad_pos",
                        )
                    )
                    continue
                if pos < 1:
                    records.append(
                        VcfRecord(
                            row_index=row_index,
                            line_no=line_no,
                            raw=line,
                            malformed_reason="nonpositive_pos",
                        )
                    )
                    continue

                explicit_id = raw_id not in {"", "."}
                record_id = raw_id if explicit_id else f"VCF_RECORD_{line_no}"
                records.append(
                    VcfRecord(
                        row_index=row_index,
                        line_no=line_no,
                        chrom=chrom,
                        pos=pos,
                        record_id=record_id,
                        explicit_id=explicit_id,
                        ref=ref,
                        alt=alt.split(",", 1)[0],
                        qual=qual,
                        filt=filt,
                        info=parse_vcf_info(info_text),
                        raw=line,
                    )
                )
    except (OSError, UnicodeError) as exc:
        raise VcfParseError(f"cannot read VCF {path}: {exc}") from exc

    if column_header_count != 1:
        raise VcfParseError(
            f"expected exactly one #CHROM header, found {column_header_count}"
        )
    return header_contigs, records, source_lines, dict(version_headers)


def _version_parts(value: str) -> tuple[int, ...] | None:
    match = re.search(r"\d+(?:\.\d+)*", value)
    if match is None:
        return None
    try:
        return tuple(int(part) for part in match.group(0).split("."))
    except ValueError:
        return None


def _version_before(parts: tuple[int, ...], target: tuple[int, ...]) -> bool:
    width = max(len(parts), len(target))
    padded_parts = parts + (0,) * (width - len(parts))
    padded_target = target + (0,) * (width - len(target))
    return padded_parts < padded_target


def _detect_source(
    source_lines: Sequence[str],
    version_headers: Mapping[str, Sequence[str]],
) -> SourceMetadata:
    matches: list[tuple[str, str | None, tuple[int, ...] | None]] = []
    for source in source_lines:
        source_value = source.strip()
        lower = source_value.lower()
        matched_prefix_and_caller = next(
            (
                (prefix, caller)
                for prefix, caller in SOURCE_PREFIX_CALLERS
                if lower.startswith(prefix)
            ),
            None,
        )
        if matched_prefix_and_caller is None:
            continue
        prefix, caller = matched_prefix_and_caller
        # Sniffles2 contains a digit in the caller name itself, so version
        # extraction must start after the matched caller prefix.
        parts = _version_parts(source_value[len(prefix) :])
        version = ".".join(map(str, parts)) if parts is not None else None
        matches.append((caller, version, parts))

    callers = {item[0] for item in matches}
    if len(callers) > 1:
        found = ", ".join(source_lines) if source_lines else "none"
        logging.warning(
            "VCF source is not a uniquely supported caller "
            f"({', '.join(SUPPORTED_CALLERS)}); compatibility is not "
            f"guaranteed. Continuing with generic ALT parsing. Found: {found}"
        )
        return SourceMetadata(tuple(source_lines), None, None, None)

    if len(callers) == 1:
        caller = next(iter(callers))
        versions = {(item[1], item[2]) for item in matches if item[0] == caller}
        if len(versions) == 1:
            version, parts = next(iter(versions))
        else:
            version, parts = None, None
            logging.warning(
                f"VCF has ambiguous {caller} source versions; caller-specific "
                "version handling is disabled"
            )
        logging.info(
            f"VCF caller detected: {caller}"
            + (f" {version}" if version is not None else " (version unknown)")
        )
        if caller == "severus" and parts is None:
            logging.warning(
                "Severus source version could not be determined; Severus BND "
                "records will be excluded because pre-1.7 direction correction "
                "cannot be selected safely"
            )
        return SourceMetadata(tuple(source_lines), caller, version, parts)

    # GRIPSS commonly preserves an upstream GRIDSS header without adding a
    # ##source line, but does add its own unambiguous version key.
    version_key = next(
        (
            key
            for key in ("gripssversion",)
            if version_headers.get(key)
        ),
        None,
    )
    if version_key is not None:
        caller = VERSION_HEADER_CALLERS[version_key]
        version_header_name = VERSION_HEADER_NAMES[version_key]
        parsed_versions = {
            (
                ".".join(map(str, parsed)) if parsed is not None else None,
                parsed,
            )
            for raw_version in version_headers[version_key]
            for parsed in (_version_parts(raw_version),)
        }
        if len(parsed_versions) == 1:
            version, parts = next(iter(parsed_versions))
        else:
            version, parts = None, None
            logging.warning(
                f"VCF has ambiguous ##{version_header_name} values; caller "
                "version is "
                "reported as unknown"
            )
        found = ", ".join(source_lines) if source_lines else "none"
        logging.warning(
            f"VCF ##source is missing or unsupported (found: {found}); detected "
            f"{caller} from ##{version_header_name}. Standard BND ALT parsing "
            "will be used."
        )
        logging.info(
            f"VCF caller detected: {caller}"
            + (f" {version}" if version is not None else " (version unknown)")
        )
        return SourceMetadata(tuple(source_lines), caller, version, parts)

    found = ", ".join(source_lines) if source_lines else "none"
    logging.warning(
        "VCF source is not a uniquely supported caller "
        f"({', '.join(SUPPORTED_CALLERS)}); compatibility is not guaranteed. "
        f"Continuing with generic ALT parsing. Found: {found}"
    )
    return SourceMetadata(tuple(source_lines), None, None, None)


def _validate_reference(
    header_contigs: Mapping[str, int],
    records: Sequence[VcfRecord],
    reference_lengths: Mapping[str, int],
) -> tuple[set[str], set[str]]:
    errors: list[str] = []
    for chrom, length in sorted(header_contigs.items()):
        if chrom in reference_lengths and int(reference_lengths[chrom]) != int(length):
            errors.append(
                f"{chrom}: VCF length {length} != FAI length {reference_lengths[chrom]}"
            )

    used_chroms = {
        record.chrom
        for record in records
        if record.malformed_reason is None and record.chrom
    }
    for record in records:
        if record.malformed_reason is not None:
            continue
        parsed = parse_vcf_bnd_alt(record.alt)
        if parsed is not None:
            used_chroms.add(parsed.mate_chrom)
    if errors:
        raise VcfParseError("VCF/reference mismatch:\n" + "\n".join(errors))

    missing_header_contigs = set(header_contigs) - set(reference_lengths)
    missing_used_contigs = used_chroms - set(reference_lengths)
    if missing_header_contigs:
        names = ", ".join(sorted(missing_header_contigs)[:8])
        if len(missing_header_contigs) > 8:
            names += f", ... (+{len(missing_header_contigs) - 8})"
        logging.warning(
            f"VCF header contains {len(missing_header_contigs)} contig(s) absent "
            f"from reference_fai ({names}). Records using absent contigs will be "
            "excluded; shared-contig length mismatches remain fatal."
        )
    if missing_used_contigs:
        names = ", ".join(sorted(missing_used_contigs)[:8])
        if len(missing_used_contigs) > 8:
            names += f", ... (+{len(missing_used_contigs) - 8})"
        logging.warning(
            f"VCF records use {len(missing_used_contigs)} contig(s) absent from "
            f"reference_fai; affected events will be excluded: {names}"
        )
    return missing_header_contigs, missing_used_contigs


def _info_string(record: VcfRecord, key: str) -> str | None:
    if record.info is None:
        return None
    value = record.info.get(key)
    if value in {None, True, ""}:
        return None
    return str(value)


def _record_mate_id(record: VcfRecord) -> tuple[str | None, bool]:
    values = {
        value.split(",", 1)[0]
        for key in ("MATE_ID", "MATEID")
        if (value := _info_string(record, key)) is not None
    }
    if len(values) == 1:
        return next(iter(values)), False
    return None, len(values) > 1


def _normalize_breakend(
    record: VcfRecord,
    parsed: BreakendAlt | None,
    source: SourceMetadata,
) -> _NormalizedBreakend:
    if parsed is None:
        return _NormalizedBreakend(None, "malformed_BND_ALT")
    if source.caller != "severus":
        return _NormalizedBreakend(parsed)
    if source.version_parts is None:
        return _NormalizedBreakend(None, "unknown_Severus_version")
    if not _version_before(source.version_parts, SEVERUS_FIXED_VERSION):
        return _NormalizedBreakend(parsed)

    strands = _info_string(record, "STRANDS")
    if strands not in VALID_SEVERUS_STRANDS:
        return _NormalizedBreakend(
            None,
            "missing_or_invalid_Severus_STRANDS",
        )
    if strands in {"+-", "-+"}:
        # Before 1.7, Severus selected the mixed-strand ALT branch in bp2/bp1
        # order.  The two ALT records remain reciprocal, but each endpoint's
        # semantic direction is reversed.  Inverting both parsed directions
        # corrects either endpoint without depending on record order.
        return _NormalizedBreakend(
            replace(
                parsed,
                dir_a=invert_vcf_strand(parsed.dir_a),
                dir_b=invert_vcf_strand(parsed.dir_b),
            ),
            severus_corrected=True,
        )
    return _NormalizedBreakend(parsed)


def _reciprocal_coordinates(
    first: VcfRecord,
    second: VcfRecord,
    raw_alts: Mapping[int, BreakendAlt | None],
) -> bool:
    first_alt = raw_alts.get(first.row_index)
    second_alt = raw_alts.get(second.row_index)
    if first_alt is None or second_alt is None:
        return False
    return (
        (first_alt.mate_chrom, first_alt.mate_pos) == (second.chrom, second.pos)
        and (second_alt.mate_chrom, second_alt.mate_pos) == (first.chrom, first.pos)
    )


def _resolve_mate(
    record: VcfRecord,
    by_id: Mapping[str, list[VcfRecord]],
    by_coordinate: Mapping[tuple[str, int], list[VcfRecord]],
    raw_alts: Mapping[int, BreakendAlt | None],
) -> tuple[VcfRecord | None, str | None, bool]:
    mate_id, id_conflict = _record_mate_id(record)
    if mate_id is not None:
        id_candidates = [
            candidate
            for candidate in by_id.get(mate_id, [])
            if candidate.row_index != record.row_index
        ]
        if len(id_candidates) == 1:
            # MATEID is authoritative when it resolves uniquely.  Manta can
            # shift each mate record across microhomology, so its ALT remote
            # coordinate need not equal the mate record POS even though the
            # two IDs and orientations describe one event.
            method = (
                "id"
                if _reciprocal_coordinates(record, id_candidates[0], raw_alts)
                else "id_coordinate_disagreement"
            )
            return id_candidates[0], method, id_conflict

    parsed = raw_alts.get(record.row_index)
    if parsed is None:
        return None, None, id_conflict
    coordinate_candidates = [
        candidate
        for candidate in by_coordinate.get((parsed.mate_chrom, parsed.mate_pos), [])
        if candidate.row_index != record.row_index
        and _reciprocal_coordinates(record, candidate, raw_alts)
    ]
    if len(coordinate_candidates) == 1:
        return coordinate_candidates[0], "coordinate", id_conflict
    return None, None, id_conflict


def _bnd_directions_are_reciprocal(
    first: BreakendAlt, second: BreakendAlt
) -> bool:
    return (second.dir_a, second.dir_b) == (
        invert_vcf_strand(first.dir_b),
        invert_vcf_strand(first.dir_a),
    )


def _unique_event_name(base: str, line_no: int, used: set[str]) -> str:
    if base not in used:
        used.add(base)
        return base
    candidate = f"{base}_L{line_no}"
    serial = 2
    while candidate in used:
        candidate = f"{base}_L{line_no}_{serial}"
        serial += 1
    used.add(candidate)
    return candidate


def _inversion_junction_roles(
    record: VcfRecord,
    source: SourceMetadata,
) -> tuple[tuple[tuple[str, str, str], ...], str | None]:
    """Return the junction(s) represented by one symbolic INV record.

    A conventional symbolic INV record represents both reciprocal inversion
    junctions.  Manta's ``convertInversion.py`` instead converts each BND pair
    into one ``<INV>`` record and marks that single junction with ``INV3`` or
    ``INV5``.  Expanding each of those records into both junctions would count
    the same inversion boundaries twice.
    """

    full_inversion = (
        ("left", "+", "-"),
        ("right", "-", "+"),
    )
    if source.caller != "manta":
        return full_inversion, None

    info = record.info or {}
    has_inv3 = "INV3" in info
    has_inv5 = "INV5" in info
    if has_inv3 and has_inv5:
        return (), "conflicting_Manta_INV3_INV5"
    if not has_inv3 and not has_inv5:
        return (), "missing_Manta_INV3_INV5"
    if has_inv3:
        return (("left", "+", "-"),), None
    return (("right", "-", "+"),), None


def _bnd_event_name(
    record: VcfRecord,
    mate: VcfRecord | None,
    used_names: set[str],
) -> str:
    base = f"vcf_bnd_{sanitize_vcf_id(record.record_id)}"
    if mate is not None:
        base += f"_{sanitize_vcf_id(mate.record_id)}"
    return _unique_event_name(base, record.line_no, used_names)


def _bnd_type4_event(
    record: VcfRecord,
    mate: VcfRecord | None,
    parsed: BreakendAlt,
    event_name: str,
) -> tuple[dict | None, str | None]:
    if record.chrom != parsed.mate_chrom:
        return None, None
    st, nd = sorted((record.pos, parsed.mate_pos))
    span = nd - st
    if parsed.dir_a == "+" and parsed.dir_b == "+":
        event_type = "front_jump" if record.pos < parsed.mate_pos else "back_jump"
    elif parsed.dir_a == "-" and parsed.dir_b == "-":
        event_type = "back_jump" if record.pos < parsed.mate_pos else "front_jump"
    else:
        return None, None
    if span < VCF_TYPE4_MIN_SPAN:
        return None, "small_indel_size"
    mate_id = mate.record_id if mate is not None else None
    return {
        "event_id": event_name,
        "vcf_id": record.record_id,
        "mate_id": mate_id,
        "svtype": "BND",
        "event_type": event_type,
        "indel_kind": "deletion" if event_type == "front_jump" else "insertion",
        "chrom": record.chrom,
        "st": st,
        "nd": nd,
        "svlen": span,
        "line_no": record.line_no,
        "source": f"VCF_BND_{event_name}",
    }, None


def _type4_event_from_record(
    record: VcfRecord,
    reference_lengths: Mapping[str, int],
    caller: str | None,
    ins_alt_alignments: Mapping[str, Mapping[str, object]] | None,
) -> tuple[dict | None, str | None]:
    info = record.info or {}
    svtype = record.svtype
    parsed_svlen = parse_optional_int(info.get("SVLEN"))
    if parsed_svlen is None and caller == "nanomonsv" and svtype == "INS":
        parsed_svlen = parse_optional_int(info.get("SVINSLEN"))
    svlen = abs(parsed_svlen or 0)

    if svtype in {"DEL", "DUP"}:
        end = parse_optional_int(info.get("END"))
        if end is None:
            # Some callers omit END for short sequence-resolved alleles.  They
            # can still be rejected safely by SVLEN, but a large event cannot
            # be handed to step 11 without both reference breakpoints.
            if 0 < svlen < VCF_TYPE4_MIN_SPAN:
                return None, "small_indel_size"
            return None, "missing_END"
        if record.chrom not in reference_lengths:
            return None, "unknown_chrom"
        st, nd = sorted((record.pos, end))
        span = nd - st
        # END/POS describe the two reference breakpoints consumed by step 11,
        # so their span is authoritative when SVLEN disagrees with them.
        if span < VCF_TYPE4_MIN_SPAN:
            return None, "small_indel_size"
        event_type = "front_jump" if svtype == "DEL" else "back_jump"
        return {
            "event_id": sanitize_vcf_id(record.record_id),
            "vcf_id": record.record_id,
            "svtype": svtype,
            "event_type": event_type,
            "indel_kind": "deletion" if svtype == "DEL" else "insertion",
            "chrom": record.chrom,
            "st": st,
            "nd": nd,
            "svlen": span,
            "reported_svlen": parsed_svlen,
            "line_no": record.line_no,
            "source": f"VCF_{sanitize_vcf_id(record.record_id)}",
        }, None

    if svtype == "INS":
        # INS eligibility is determined only by the caller-reported insertion
        # length.  The ALT alignment supplies the reference path used by the
        # type4 event, but its reference span must not promote a smaller (or
        # unknown-size) insertion across the 100 kb threshold.
        if svlen < VCF_TYPE4_MIN_SPAN:
            return None, "small_indel_size"
        query_name = vcf_ins_query_name(record)
        aligned = (ins_alt_alignments or {}).get(query_name)
        if aligned is None:
            return None, "INS_no_alt_alignment"
        chrom = str(aligned["chrom"])
        st = int(aligned["st"])
        nd = int(aligned["nd"])
        if chrom not in reference_lengths:
            return None, "unknown_alt_alignment_chrom"
        return {
            "event_id": sanitize_vcf_id(record.record_id),
            "vcf_id": record.record_id,
            "svtype": svtype,
            "event_type": "back_jump",
            "indel_kind": "insertion",
            "chrom": chrom,
            "st": st,
            "nd": nd,
            "svlen": svlen,
            "line_no": record.line_no,
            "query_name": query_name,
            "paf_rows": aligned.get("paf_rows", []),
            "source": f"VCF_{sanitize_vcf_id(record.record_id)}",
        }, None
    return None, "not_type4_svtype"


def _initial_summary(
    records: Sequence[VcfRecord],
    source: SourceMetadata,
    pass_filters: Sequence[str],
) -> dict[str, object]:
    counters = {
        metric: 0
        for metric in (
            "used_bnd_events",
            "used_bnd_type4_events",
            "used_inv_events",
            "used_inv_junctions",
            "used_manta_inv3_records",
            "used_manta_inv5_records",
            "ambiguous_manta_inv_records",
            "used_type4_events",
            "skipped_ins",
            "skipped_pure_ins",
            "skipped_ins_small_size",
            "skipped_ins_no_alt_alignment",
            "skipped_small_indel_size",
            "skipped_bnd_small_indel_events",
            "skipped_not_type4_svtype",
            "missing_mates",
            "malformed_records",
            "orientation_mismatches",
            "filtered_records",
            "filtered_bnd_pairs",
            "singleton_bnd_events",
            "coordinate_fallback_pairs",
            "mate_id_coordinate_disagreements",
            "mate_id_conflicts",
            "duplicate_mate_records",
            "single_breakends_without_remote_locus",
            "reference_header_contigs_absent",
            "reference_used_contigs_absent",
            "unknown_reference_contig_records",
            "unknown_reference_contig_events",
            "severus_corrected_records",
            "severus_bracket_alt_overrides",
            "missing_or_invalid_severus_strands",
            "severus_strands_skipped_events",
            "unknown_severus_version_bnd_records",
        )
    }
    counters.update(
        {
            "vcf_records": sum(
                record.malformed_reason is None for record in records
            ),
            "vcf_source": " | ".join(source.source_lines) or "none",
            "vcf_caller": source.caller or "unknown",
            "vcf_caller_version": source.version or "unknown",
            "vcf_filter_pass": ",".join(pass_filters),
        }
    )
    return counters


def parse_vcf_events(
    vcf_path: str | Path,
    reference_lengths: Mapping[str, int],
    pass_filters: Collection[str] = DEFAULT_PASS_FILTERS,
    ins_alt_alignments: Mapping[str, Mapping[str, object]] | None = None,
) -> VcfParseResult:
    """Parse a plain-text VCF into normalized SKYPE event specifications."""

    pass_filter_values = tuple(dict.fromkeys(str(value) for value in pass_filters))
    allowed_filters = set(pass_filter_values)
    header_contigs, records, source_lines, version_headers = _read_vcf(vcf_path)
    source = _detect_source(source_lines, version_headers)
    missing_header_contigs, missing_used_contigs = _validate_reference(
        header_contigs, records, reference_lengths
    )

    summary = _initial_summary(records, source, pass_filter_values)
    summary["reference_header_contigs_absent"] = len(missing_header_contigs)
    summary["reference_used_contigs_absent"] = len(missing_used_contigs)
    skipped_records: list[tuple[object, ...]] = []
    orientation_mismatches: list[tuple[object, ...]] = []
    nclose_specs: list[NCloseSpec] = []
    type4_events: list[dict] = []
    used_event_names: set[str] = set()

    def skip(record: VcfRecord, reason: str) -> None:
        skipped_records.append(
            (record.line_no, record.record_id, record.svtype or ".", reason, record.raw)
        )

    for record in records:
        if record.malformed_reason is not None:
            summary["malformed_records"] += 1
            skip(record, record.malformed_reason)

    valid_records = [record for record in records if record.malformed_reason is None]
    severus_bracket_alts = {
        record.row_index: parsed
        for record in valid_records
        if source.caller == "severus" and record.svtype != "BND"
        if (parsed := parse_vcf_bnd_alt(record.alt)) is not None
    }
    bnd_row_indices = {
        record.row_index
        for record in valid_records
        if record.svtype == "BND" or record.row_index in severus_bracket_alts
    }
    bnd_records = [
        record for record in valid_records if record.row_index in bnd_row_indices
    ]
    summary["severus_bracket_alt_overrides"] = len(severus_bracket_alts)
    summary["mate_id_conflicts"] = sum(
        _record_mate_id(record)[1] for record in bnd_records
    )
    raw_alts = {
        record.row_index: (
            severus_bracket_alts[record.row_index]
            if record.row_index in severus_bracket_alts
            else parse_vcf_bnd_alt(record.alt)
        )
        for record in bnd_records
    }
    normalized_alts = {
        record.row_index: _normalize_breakend(
            record, raw_alts[record.row_index], source
        )
        for record in bnd_records
    }
    by_id: dict[str, list[VcfRecord]] = defaultdict(list)
    by_coordinate: dict[tuple[str, int], list[VcfRecord]] = defaultdict(list)
    for record in bnd_records:
        if record.explicit_id:
            by_id[record.record_id].append(record)
        by_coordinate[(record.chrom, record.pos)].append(record)

    consumed_bnd_rows: set[int] = set()
    for record in valid_records:
        if record.row_index not in bnd_row_indices:
            if record.filt not in allowed_filters:
                summary["filtered_records"] += 1
                skip(record, f"FILTER={record.filt}")
                continue

            if record.chrom not in reference_lengths:
                summary["unknown_reference_contig_records"] += 1
                summary["unknown_reference_contig_events"] += 1
                skip(record, f"unknown_reference_contig={record.chrom}")
                continue

            if record.svtype == "INV":
                end = parse_optional_int((record.info or {}).get("END"))
                if end is None:
                    summary["malformed_records"] += 1
                    skip(record, "missing_END")
                    continue
                pos_a, pos_b = sorted((record.pos, end))
                junction_roles, inv_error = _inversion_junction_roles(
                    record, source
                )
                if inv_error is not None:
                    summary["malformed_records"] += 1
                    summary["ambiguous_manta_inv_records"] += 1
                    skip(record, inv_error)
                    continue
                common = {
                    "record_ids": (record.record_id,),
                    "line_nos": (record.line_no,),
                    "svtype": "INV",
                    "chrom_a": record.chrom,
                    "pos_a": pos_a,
                    "chrom_b": record.chrom,
                    "pos_b": pos_b,
                }
                for role, dir_a, dir_b in junction_roles:
                    event_name = _unique_event_name(
                        f"vcf_inv_{sanitize_vcf_id(record.record_id)}_{role}",
                        record.line_no,
                        used_event_names,
                    )
                    nclose_specs.append(
                        NCloseSpec(
                            event_name=event_name,
                            dir_a=dir_a,
                            dir_b=dir_b,
                            **common,
                        )
                    )
                summary["used_inv_events"] += 1
                summary["used_inv_junctions"] += len(junction_roles)
                if source.caller == "manta":
                    metric = (
                        "used_manta_inv3_records"
                        if "INV3" in (record.info or {})
                        else "used_manta_inv5_records"
                    )
                    summary[metric] += 1
                continue

            event, reason = _type4_event_from_record(
                record,
                reference_lengths,
                source.caller,
                ins_alt_alignments,
            )
            if event is not None:
                type4_events.append(event)
                summary["used_type4_events"] += 1
            elif record.svtype == "INS" and reason == "INS_no_alt_alignment":
                summary["skipped_ins"] += 1
                summary["skipped_ins_no_alt_alignment"] += 1
                skip(record, reason)
            elif record.svtype in {"DEL", "DUP", "INS"} and reason == "small_indel_size":
                if record.svtype == "INS":
                    summary["skipped_ins"] += 1
                    summary["skipped_ins_small_size"] += 1
                summary["skipped_small_indel_size"] += 1
                skip(record, reason)
            elif reason == "not_type4_svtype":
                summary["skipped_not_type4_svtype"] += 1
                skip(record, reason)
            elif reason in {"unknown_chrom", "unknown_alt_alignment_chrom"}:
                summary["unknown_reference_contig_records"] += 1
                summary["unknown_reference_contig_events"] += 1
                skip(record, reason)
            elif record.svtype in {"DEL", "DUP", "INS"}:
                summary["malformed_records"] += 1
                skip(record, reason or "malformed_type4_record")
            continue

        if record.row_index in consumed_bnd_rows:
            continue
        raw_alt = raw_alts.get(record.row_index)
        if raw_alt is None:
            consumed_bnd_rows.add(record.row_index)
            if is_vcf_single_breakend_alt(record.alt):
                summary["single_breakends_without_remote_locus"] += 1
                skip(record, "single_breakend_without_remote_locus")
            else:
                summary["malformed_records"] += 1
                skip(record, "malformed_BND_ALT")
            continue

        mate, mate_method, _ = _resolve_mate(
            record, by_id, by_coordinate, raw_alts
        )
        if mate is not None and mate.row_index in consumed_bnd_rows:
            consumed_bnd_rows.add(record.row_index)
            summary["duplicate_mate_records"] += 1
            skip(record, "duplicate_mate_event")
            continue

        members = [record] if mate is None else [record, mate]
        consumed_bnd_rows.update(member.row_index for member in members)
        if mate is None:
            summary["missing_mates"] += 1
        elif mate_method == "coordinate":
            summary["coordinate_fallback_pairs"] += 1
        elif mate_method == "id_coordinate_disagreement":
            summary["mate_id_coordinate_disagreements"] += 1

        if any(member.filt not in allowed_filters for member in members):
            summary["filtered_records"] += len(members)
            if mate is not None:
                summary["filtered_bnd_pairs"] += 1
            for member in members:
                reason = (
                    f"FILTER={member.filt}"
                    if member.filt not in allowed_filters
                    else "mate_filtered"
                )
                skip(member, reason)
            continue

        event_chroms = {member.chrom for member in members}
        for member in members:
            member_alt = raw_alts.get(member.row_index)
            if member_alt is not None:
                event_chroms.add(member_alt.mate_chrom)
        unknown_event_chroms = sorted(
            chrom for chrom in event_chroms if chrom not in reference_lengths
        )
        if unknown_event_chroms:
            summary["unknown_reference_contig_records"] += len(members)
            summary["unknown_reference_contig_events"] += 1
            reason = "unknown_reference_contig=" + ",".join(unknown_event_chroms)
            for member in members:
                skip(member, reason)
            continue

        normalized_members = [normalized_alts[member.row_index] for member in members]
        errors = [item.error for item in normalized_members if item.error is not None]
        if errors:
            reason = errors[0]
            if reason == "missing_or_invalid_Severus_STRANDS":
                summary["missing_or_invalid_severus_strands"] += sum(
                    item.error == reason for item in normalized_members
                )
                summary["severus_strands_skipped_events"] += 1
            elif reason == "unknown_Severus_version":
                summary["unknown_severus_version_bnd_records"] += len(members)
            else:
                summary["malformed_records"] += len(members)
            for member in members:
                skip(member, reason)
            continue

        first_alt = normalized_members[0].alt
        assert first_alt is not None
        if mate is not None:
            second_alt = normalized_members[1].alt
            assert second_alt is not None
            if not _bnd_directions_are_reciprocal(first_alt, second_alt):
                summary["orientation_mismatches"] += 1
                orientation_mismatches.append(
                    (
                        record.line_no,
                        record.record_id,
                        mate.line_no,
                        mate.record_id,
                        first_alt.dir_a + first_alt.dir_b,
                        second_alt.dir_a + second_alt.dir_b,
                        record.alt,
                        mate.alt,
                    )
                )
                skip(record, "BND_ALT_orientation_mismatch")
                skip(mate, "BND_ALT_orientation_mismatch")
                continue

        summary["severus_corrected_records"] += sum(
            item.severus_corrected for item in normalized_members
        )
        event_name = _bnd_event_name(record, mate, used_event_names)
        type4_event, bnd_indel_reason = _bnd_type4_event(
            record, mate, first_alt, event_name
        )
        if type4_event is not None:
            type4_events.append(type4_event)
            summary["used_type4_events"] += 1
            summary["used_bnd_type4_events"] += 1
            if mate is None:
                summary["singleton_bnd_events"] += 1
            continue
        if bnd_indel_reason == "small_indel_size":
            summary["skipped_small_indel_size"] += 1
            summary["skipped_bnd_small_indel_events"] += 1
            for member in members:
                skip(member, bnd_indel_reason)
            continue

        nclose_specs.append(
            NCloseSpec(
                event_name=event_name,
                record_ids=tuple(member.record_id for member in members),
                line_nos=tuple(member.line_no for member in members),
                svtype="BND",
                chrom_a=record.chrom,
                pos_a=record.pos,
                dir_a=first_alt.dir_a,
                chrom_b=first_alt.mate_chrom,
                pos_b=first_alt.mate_pos,
                dir_b=first_alt.dir_b,
            )
        )
        summary["used_bnd_events"] += 1
        if mate is None:
            summary["singleton_bnd_events"] += 1

    return VcfParseResult(
        source=source,
        nclose_specs=nclose_specs,
        type4_events=type4_events,
        summary=summary,
        skipped_records=skipped_records,
        orientation_mismatches=orientation_mismatches,
    )
