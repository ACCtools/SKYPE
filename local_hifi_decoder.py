"""Isolated molecule decoder v2; the shared v1 decoder remains immutable.

Changes from frozen v1 dbf16221a148b18235991f34b02efd2727fb8327faf8d7d2c6588a007e8f8bac:
1. Explicit CIGAR I/D gaps >=50bp are represented as separate pieces.
2. MAPQ>=20 short pieces remain in query order before anchor eligibility is
   checked, so they cannot be bypassed to invent an outer adjacency.

Observed records precede redundant SA summaries. Query coordinates include
hard clips and are expressed on the full forward molecule. Endpoint positions
are 0-based reference boundaries with retained sides L/R. This is alignment
evidence, not a calibrated variant probability or independent validation.
"""
import re

POLICY = dict(version="candidate_molecule_decoder_v2",
              cigar_gap_min_bp=50, mapq_min=20, anchor_min_bp=500,
              query_gap_range_bp=[-500, 1000], minimum_gap_discrepancy_bp=50,
              intervening_short_high_mapq_pieces="preserve_before_anchor_filter",
              SA_precedence="observed_record_with_same_source_alignment_identity",
              anchor_definition="min(aligned_query_span,reference_span)")


def endpoint(piece, leaving):
    left = (piece["strand"] == "+") == leaving
    return (piece["chrom"], piece["end0"] if left else piece["start0"], "L" if left else "R")


def cigar_alignment(chrom, start0, strand, cigar, mapq):
    ops = [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]
    if not ops:
        return None
    qlen = sum(n for n, op in ops if op in "MISH=X")
    ref_len = sum(n for n, op in ops if op in "MDN=X")
    aligned = sum(n for n, op in ops if op in "MI=X")
    leading = 0
    for n, op in ops:
        if op not in "SH":
            break
        leading += n
    qstart = leading if strand == "+" else qlen-leading-aligned
    return dict(chrom=chrom, start0=int(start0), end0=int(start0)+ref_len,
                strand=strand, qstart=qstart, qend=qstart+aligned, mapq=int(mapq),
                anchor=min(aligned, ref_len), cigar=cigar)


def cigar_pieces(chrom, start0, strand, cigar, mapq, gap_min=50):
    row = cigar_alignment(chrom, start0, strand, cigar, mapq)
    if row is None:
        return []
    ops = [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar)]
    if not any(op == "N" or (op in "ID" and n >= gap_min) for n, op in ops):
        return [row]
    qlen = sum(n for n, op in ops if op in "MISH=X")
    qpos, rpos = 0, int(start0)
    qbegin, rbegin = qpos, rpos
    result, block = [], []

    def flush():
        if qpos <= qbegin or rpos <= rbegin:
            return
        qa, qb = (qbegin, qpos) if strand == "+" else (qlen-qpos, qlen-qbegin)
        result.append(dict(chrom=chrom, start0=rbegin, end0=rpos, strand=strand,
                           qstart=qa, qend=qb, mapq=int(mapq), anchor=min(qb-qa, rpos-rbegin),
                           cigar="".join(f"{n}{op}" for n, op in block)))

    for n, op in ops:
        boundary = op in "SHN" or (op in "ID" and n >= gap_min)
        if boundary:
            flush()
        if op in "MISH=X":
            qpos += n
        if op in "MDN=X":
            rpos += n
        if boundary:
            qbegin, rbegin, block = qpos, rpos, []
        else:
            block.append((n, op))
    flush()
    return result


def read_pieces(record):
    pieces = []

    def add(chrom, start, strand, cigar, mapq, origin):
        whole = cigar_alignment(chrom, start, strand, cigar, mapq)
        if whole is None:
            return
        identity = (chrom, start, strand, whole["qstart"], whole["qend"])
        pieces.extend(dict(p, origin=origin, source_alignment=identity)
                      for p in cigar_pieces(chrom, start, strand, cigar, mapq))

    add(record.reference_name, record.reference_start,
        "-" if record.is_reverse else "+", record.cigarstring, record.mapping_quality, "record")
    if record.has_tag("SA"):
        for entry in record.get_tag("SA").rstrip(";").split(";"):
            fields = entry.split(",")
            if len(fields) >= 6:
                add(fields[0], int(fields[1])-1, fields[2], fields[3], int(fields[4]), "SA")
    return pieces


def normalized_pieces(pieces, min_mapq=20):
    """Return ordered unique high-MAPQ pieces, including short separators."""
    observed = {p["source_alignment"] for p in pieces if p.get("origin") == "record"}
    authoritative = [p for p in pieces if p.get("origin") != "SA" or p.get("source_alignment") not in observed]
    unique = {(p["chrom"], p["start0"], p["end0"], p["strand"], p["qstart"], p["qend"]): p
              for p in authoritative if p["mapq"] >= min_mapq}
    return sorted(unique.values(), key=lambda p: (p["qstart"], p["qend"], p["chrom"], p["start0"]))


def connections(pieces, min_mapq=20, min_anchor=500):
    ordered = normalized_pieces(pieces, min_mapq)
    links = []
    for a, b in zip(ordered, ordered[1:]):
        gap = b["qstart"]-a["qend"]
        if a["anchor"] < min_anchor or b["anchor"] < min_anchor or not -500 <= gap <= 1000:
            continue
        if a["chrom"] == b["chrom"] and a["strand"] == b["strand"]:
            rgap = b["start0"]-a["end0"] if a["strand"] == "+" else a["start0"]-b["end0"]
            if abs(rgap-gap) < 50:
                continue
        links.append(sorted([endpoint(a, True), endpoint(b, False)]))
    return links


def matches(left, right, tolerance):
    return all(a[0] == b[0] and a[2] == b[2] and abs(a[1]-b[1]) <= tolerance
               for a, b in zip(left, right))


def flank_boundary_ranges(piece, anchor=500):
    """Boundaries with >=anchor query AND reference span on the retained side.

    M/I/=/X consume aligned query, M/D/N/=/X consume reference, and clips
    contribute to neither flank. This follows the decoder's span-based anchor
    definition; it does not require anchor exact nucleotide matches.
    """
    ops = [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", piece["cigar"])]
    qspan = sum(n for n, op in ops if op in "MI=X")
    rspan = sum(n for n, op in ops if op in "MDN=X")
    assert qspan == piece["qend"]-piece["qstart"]
    assert rspan == piece["end0"]-piece["start0"]
    if min(qspan, rspan) < anchor:
        return {}
    thresholds = []
    for reverse in [False, True]:
        cursor = piece["end0"] if reverse else piece["start0"]
        direction, query = (-1 if reverse else 1), 0
        threshold = None
        for length, op in (ops[::-1] if reverse else ops):
            if op in "MI=X":
                if query < anchor <= query+length:
                    threshold = cursor + (direction*(anchor-query) if op in "M=X" else 0)
                    break
                query += length
            if op in "MDN=X":
                cursor += direction*length
        assert threshold is not None
        thresholds.append(threshold)
    return dict(L=(max(piece["start0"]+anchor, thresholds[0]), piece["end0"]),
                R=(piece["start0"], min(piece["end0"]-anchor, thresholds[1])))
