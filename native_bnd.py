"""Native BND representation of selected NClose and compound events.

This layer does not discover or split NCloses. A compound event contributes
its original two NCloses, never an inferred bridge between different unitigs.
Weight accounting lives in structure_nclose; these helpers only describe
original junction geometry.
"""

from __future__ import annotations


Endpoint = tuple[str, int, str]  # chromosome, zero-based boundary, retained L/R


def node_pair_endpoints(contig_data, pair):
    a, b = sorted((contig_data[i] for i in pair), key=lambda n: (n[2], n[3]))
    return (
        (a[5], int(a[8] if a[4] == "+" else a[7]), "L" if a[4] == "+" else "R"),
        (b[5], int(b[7] if b[4] == "+" else b[8]), "R" if b[4] == "+" else "L"),
    ), tuple(dict.fromkeys((str(a[0]), str(b[0]))))


def layout_endpoints(layout):
    # Use the original junction-facing coordinates, not display_points or the
    # opposite ends of their alignment anchors.
    a, b = layout["ordered_endpoints"]
    return (
        (a["chrom"], int(a["coord"]), "L" if a["dir"] == "+" else "R"),
        (b["chrom"], int(b["coord"]), "R" if b["dir"] == "+" else "L"),
    ), tuple(dict.fromkeys((a["ctg_name"], b["ctg_name"])))


def constituent_pairs(circuit):
    s1, e1, s2, e2 = circuit
    return (tuple(sorted((s1, e1))), tuple(sorted((s2, e2))))


def merged_type4_pairs(event, type2_ins_del):
    # Stage 11 and 21 number insertion circuits first, then deletion circuits,
    # starting at one, including circuits suppressed by step-11 deduplication.
    insertion, deletion = type2_ins_del
    circuit = (list(insertion) + list(deletion))[int(event["type2_merge_idx"]) - 1]
    return constituent_pairs(circuit)
