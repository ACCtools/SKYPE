"""Fixed, descriptive assembly-repeat annotations; no calling or fit decisions."""
from __future__ import annotations

import json
import hashlib
from pathlib import Path
import re

import numpy as np

RULES_PATH = Path(__file__).with_name("terminal_evidence_rules.json")
_RULES_BYTES = RULES_PATH.read_bytes()
RULES = json.loads(_RULES_BYTES)
RULES_SHA256 = hashlib.sha256(_RULES_BYTES).hexdigest()
IMPLEMENTATION_SHA256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(sequence):
    return sequence.translate(_COMPLEMENT)[::-1]


def motif_pattern(motifs):
    rotations = {word[i:] + word[:i] for word in motifs for i in range(len(word))}
    return re.compile("(?=(" + "|".join(sorted(rotations, key=lambda s: (-len(s), s))) + "))")


PATTERNS = {"C": motif_pattern(RULES["motifs"]),
            "G": motif_pattern([reverse_complement(s) for s in RULES["motifs"]])}
CANONICAL = {"C": motif_pattern([RULES["canonical_motif"]]),
             "G": motif_pattern([reverse_complement(RULES["canonical_motif"])])}


def coverage(sequence, pattern):
    changes = np.zeros(len(sequence) + 1, dtype=np.int32)
    for hit in pattern.finditer(sequence):
        changes[hit.start()] += 1
        changes[hit.start() + len(hit.group(1))] -= 1
    return np.cumsum(changes[:-1]) > 0


def repeat_tracts(sequence):
    """Keep tract span and actual covered bases distinct, including overlap."""
    sequence = sequence.upper()
    rule = RULES["tract_rule"]
    window = rule["window_bp"]
    if len(sequence) < window:
        return []
    result = []
    for orientation in ("C", "G"):
        covered = coverage(sequence, PATTERNS[orientation])
        canonical = coverage(sequence, CANONICAL[orientation])
        integral = np.concatenate(([0], np.cumsum(covered)))
        starts = np.flatnonzero(integral[window:] - integral[:-window] >=
                                window * rule["minimum_motif_covered_fraction"])
        intervals = []
        for start in starts.tolist():
            if intervals and start <= intervals[-1][1] + rule["join_qualifying_window_gaps_bp"]:
                intervals[-1][1] = start + window
            else:
                intervals.append([start, start + window])
        for lo, hi in intervals:
            observed = np.flatnonzero(covered[lo:hi])
            if not len(observed):
                continue
            start, end = lo + int(observed[0]), lo + int(observed[-1]) + 1
            c, total = int(canonical[start:end].sum()), int(covered[start:end].sum())
            result.append(dict(start=start, end=end, length=end - start,
                               orientation=orientation, canonical_bases=c,
                               motif_covered_bases=total,
                               primary=(end - start >= rule["primary_min_tract_bp"] and
                                        c >= rule["primary_min_canonical_bases"]),
                               sensitivity=(end - start >= rule["sensitivity_min_tract_bp"] and
                                            c >= rule["sensitivity_min_canonical_bases"])))
    return sorted(result, key=lambda row: (row["start"], row["end"], row["orientation"]))


def interval_coverage(sequence, start, end, orientation):
    """Count coverage in the full sequence so boundary-crossing hits survive."""
    sequence = sequence.upper()
    end = max(start, end)
    return dict(span_bp=end - start,
                canonical_bases=int(coverage(sequence, CANONICAL[orientation])[start:end].sum()),
                motif_covered_bases=int(coverage(sequence, PATTERNS[orientation])[start:end].sum()))


def _endpoint(piece, leaving):
    left = (piece["strand"] == "+") == leaving
    return [piece["chrom"], piece["end"] if left else piece["start"], "L" if left else "R"]


def describe_source(sequence, boundary, toward, pieces, host, tract_rows=None):
    """Describe every outward repeat/route, without assuming a chromosome end."""
    sequence = sequence.upper()
    rows = repeat_tracts(sequence) if tract_rows is None else tract_rows
    observed = len(sequence) - boundary if toward == 1 else boundary
    result = dict(observed_outward_bases=observed, adjacent_tracts=[], ordered_repeat_routes=[],
                  internal_repeat_with_farther_aligned_flank=False)
    if observed == 0:
        result["context"] = "no_outward_sequence_observed"
        return result
    lo_gap, hi_gap = RULES["adjacent_query_gap_bp"]
    for tract in rows:
        inner = tract["start"] if toward == 1 else tract["end"]
        gap = (inner - boundary) * toward
        if toward == 1 and tract["end"] <= boundary or toward == -1 and tract["start"] >= boundary:
            continue
        out_lo, out_hi = (max(boundary, tract["start"]), tract["end"]) if toward == 1 else (
            tract["start"], min(boundary, tract["end"]))
        overlap_lo, overlap_hi = (tract["start"], min(boundary, tract["end"])) if toward == 1 else (
            max(boundary, tract["start"]), tract["end"])
        low, high = sorted((boundary, inner))
        intervening, low_confidence, farther = [], [], []
        for piece in pieces:
            if piece is not host and min(high, piece["qend"]) > max(low, piece["qstart"]):
                # A repetitive alignment wholly inside the repeat is not donor identity.
                wholly_repeat = tract["start"] <= piece["qstart"] <= piece["qend"] <= tract["end"]
                if not wholly_repeat:
                    (intervening if piece["mapq"] >= RULES["minimum_host_mapq"] else low_confidence).append(piece)
            span = (piece["qend"] - max(piece["qstart"], tract["end"]) if toward == 1 else
                    min(piece["qend"], tract["start"]) - piece["qstart"])
            if (piece["mapq"] >= RULES["minimum_host_mapq"] and
                    min(span, piece["end"] - piece["start"]) >= RULES["minimum_farther_aligned_flank_bp"]):
                farther.append(piece)
        detail = dict(tract=tract, query_gap_bp=gap,
                      outward_coverage=interval_coverage(sequence, out_lo, out_hi, tract["orientation"]),
                      aligned_side_overlap_coverage=interval_coverage(sequence, overlap_lo, overlap_hi, tract["orientation"]),
                      intervening_source_pieces=intervening, low_confidence_intervening_pieces=low_confidence,
                      farther_aligned_flanks=farther,
                      internal_repeat_with_farther_aligned_flank=bool(farther),
                      sequence_tail_bp=len(sequence) - tract["end"] if toward == 1 else tract["start"])
        detail["outward_primary"] = (
            detail["outward_coverage"]["span_bp"] >= RULES["tract_rule"]["primary_min_tract_bp"] and
            detail["outward_coverage"]["canonical_bases"] >= RULES["tract_rule"]["primary_min_canonical_bases"])
        detail["sequence_terminal_like"] = detail["sequence_tail_bp"] <= RULES["source_terminal_tail_max_bp"]
        if lo_gap <= gap <= hi_gap:
            result["adjacent_tracts"].append(detail)
        if not tract["primary"] or gap < lo_gap:
            continue
        chain = [host] + sorted(intervening, key=lambda p: (p["qstart"], p["qend"], p["row_index"]), reverse=toward == -1)
        transitions = []
        for first, last in zip(chain, chain[1:]):
            qgap = last["qstart"] - first["qend"] if toward == 1 else first["qstart"] - last["qend"]
            a, b = _endpoint(first, toward == 1), _endpoint(last, toward != 1)
            reference_gap = (b[1] - a[1]) * (1 if a[2] == "L" else -1)
            transitions.append(dict(directed_endpoints=[a, b], query_gap_bp=qgap,
                                    source_row_indices=[first["row_index"], last["row_index"]],
                                    reference_continuation_compatible=(first["chrom"] == last["chrom"] and
                                        first["strand"] == last["strand"] and abs(reference_gap - qgap) <
                                        RULES["reference_continuation_discrepancy_exclusive_bp"])))
        if intervening:
            result["ordered_repeat_routes"].append(dict(**detail, ordered_source_pieces=chain,
                                                        ordered_transitions=transitions,
                                                        evidence_origin="assembly_sequence_and_selected_alignment"))
    primary = [row for row in result["adjacent_tracts"] if row["tract"]["primary"]]
    result["internal_repeat_with_farther_aligned_flank"] = any(
        row["internal_repeat_with_farther_aligned_flank"] for row in primary + result["ordered_repeat_routes"])
    direct = [row for row in primary if not row["intervening_source_pieces"] and row["outward_primary"]]
    if any(row["low_confidence_intervening_pieces"] for row in direct):
        result["context"] = "repeat_context_ambiguous"
    elif direct:
        result["context"] = "direct_repeat_extension"
    elif any(row["outward_primary"] and any(not t["reference_continuation_compatible"] for t in row["ordered_transitions"])
             for row in result["ordered_repeat_routes"]):
        result["context"] = "repeat_reached_through_ordered_aligned_source_pieces"
    elif any(row["outward_primary"] for row in result["ordered_repeat_routes"]):
        result["context"] = "repeat_reached_through_reference_compatible_source_pieces"
    elif primary or result["ordered_repeat_routes"]:
        result["context"] = "repeat_context_outward_unqualified"
    else:
        result["context"] = "no_qualifying_repeat_observed"
    return result
