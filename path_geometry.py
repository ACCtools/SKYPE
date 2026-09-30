"""Reference geometry shared by path construction and compound INDEL output."""


def walked_strand(strand, traversal):
    """FOR=1 / BAK=0; graph-only IN/OUT ports do not encode a strand."""
    if traversal not in (0, 1):
        return None
    return strand if traversal == 1 else ("-" if strand == "+" else "+")


def reference_connection_allowed(first, first_type, last, last_type):
    """Can these pieces join by trimming overlap or filling a reference gap?"""
    if first[5] != last[5]:
        return False
    a = walked_strand(first[4], first_type)
    b = walked_strand(last[4], last_type)
    if a is not None and b is not None and a != b:
        return False
    strand = a or b
    # The overlap itself can be split at a common boundary, including when
    # one anchor contains the other. The walked strands must still agree.
    if strand == "+":
        return int(first[7]) <= int(last[8])
    if strand == "-":
        return int(last[7]) <= int(first[8])
    return True


def conjoined_outer_geometry(first, last, circuit):
    """Outer breakends and depth baseline for the ordered A,B,C,D path.

    PAF strands describe the original alignments. Row order inside each
    unitig determines whether that alignment is actually walked in reverse.
    The depth baseline sign is independent of the DEL/DUP reporting label:
    overlapping duplication anchors need their union subtracted once.
    """
    a, b, c, d = circuit
    start_dir = walked_strand(first[4], int(a < b))
    end_dir = walked_strand(last[4], int(c < d))
    start_pos = int(first[8] if start_dir == "+" else first[7])
    end_pos = int(last[7] if end_dir == "+" else last[8])
    signed_gap = (end_pos - start_pos) * (1 if start_dir == "+" else -1)
    event_type = "front_jump" if signed_gap > 0 else "back_jump"
    left, right = sorted((first, last), key=lambda row: (int(row[7]), int(row[8])))
    overlapping = int(right[7]) <= int(left[8])
    subtract_base = event_type == "front_jump" or overlapping
    base_st, base_nd = (
        (min(int(first[7]), int(last[7])), max(int(first[8]), int(last[8])))
        if subtract_base else (int(left[8]), int(right[7]))
    )
    return dict(
        start_chr=first[5], start_pos=start_pos, start_dir=start_dir,
        end_chr=last[5], end_pos=end_pos, end_dir=end_dir,
        chrom=first[5], st=min(start_pos, end_pos), nd=max(start_pos, end_pos),
        event_type=event_type, signed_gap=signed_gap,
        depth_base_st=base_st, depth_base_nd=base_nd,
        depth_base_sign=-1 if subtract_base else 1,
    )


def indel_depth_vector(observed, baseline, event):
    """Preserve compound depth arithmetic independently of its output label."""
    sign = event.get("depth_base_sign", -1 if event["event_type"] == "front_jump" else 1)
    if sign not in (-1, 1):
        raise ValueError(f"Invalid INDEL baseline sign: {sign}")
    return observed + sign * baseline
