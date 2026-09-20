"""VCF-only projection from original NClose sources to exact BND adjacencies.

The depth model and its NClose identities are read-only inputs. Projection uses
each source's own alignment interval and occurrence count, before merging exact
endpoint/retained-side pairs across sources or structures.
"""

from collections import defaultdict
from math import fsum

from native_bnd import node_pair_endpoints


def is_reference_continuation(a, b):
    if a[5] != b[5] or a[4] != b[4]:
        return False
    query_gap = int(b[2]) - int(a[3])
    reference_gap = (int(b[7]) - int(a[8]) if a[4] == "+"
                     else int(a[7]) - int(b[8]))
    return reference_gap == query_gap


def source_junctions(event, nodes):
    """Describe adjacent retained alignments inside one source NClose only."""
    start, end = sorted(event["event_key"])
    if not 0 <= start < end < len(nodes):
        raise ValueError(f"Invalid NClose alignment interval: {event['event_key']}")
    owner = str(nodes[start][0])
    names = (owner, str(nodes[end][0]))
    if (event.get("rescue_method") in {"read", "olc"}
            or any(name.startswith(("raw_rescue_read_", "raw_rescue_olc_")) for name in names)):
        mode = "RAW_RESCUE_OUTER"
    elif any(name.startswith(("debug_forced_nclose_", "vcf_")) for name in names):
        mode = "SYNTHETIC_OUTER"
    else:
        mode = "ADJACENT_ALIGNMENT"
    if mode != "ADJACENT_ALIGNMENT":
        return [dict(junction_index=1, node_pair=(start, end), mode=mode,
                     endpoints=tuple(sorted(event["endpoints"])))]

    chain = list(range(start, end + 1))
    if any(str(nodes[i][0]) != owner for i in chain):
        raise ValueError(f"NClose crosses different alignment owners: {event['event_key']}")
    chain.sort(key=lambda i: (int(nodes[i][2]), int(nodes[i][3]), i))
    junctions = []
    for index, (a, b) in enumerate(zip(chain, chain[1:]), 1):
        if is_reference_continuation(nodes[a], nodes[b]):
            continue
        # A two-anchor source already has authoritative junction geometry,
        # including saved virtual-inversion endpoints.
        endpoints = (event["endpoints"] if len(chain) == 2
                     else node_pair_endpoints(nodes, (a, b))[0])
        junctions.append(dict(junction_index=index, node_pair=(a, b), mode=mode,
                              endpoints=tuple(sorted(endpoints))))
    return junctions


def nclose_bnd_memberships(context, nodes):
    """Keep all source-to-junction links, including zero/below-threshold sources."""
    memberships = {}
    for key, event in sorted(context.model["nclose_sources"].items(), key=lambda item: str(item[0])):
        if event["kind"] != "bnd":
            continue
        canonical = context.model["nclose_aliases"][key]
        memberships[key] = [dict(junction, source_nclose_key=key,
                                 nclose_id=context.ncloses[canonical]["nclose_id"])
                            for junction in source_junctions(event, nodes)]
    return memberships


def aggregate_bnd_calls(context, memberships, min_cn=0.1):
    by_source = {}
    for key, junctions in memberships.items():
        by_source[key] = defaultdict(list)
        for junction in junctions:
            by_source[key][junction["endpoints"]].append(junction["junction_index"])

    grouped = {}
    for row in context.contributions:
        if row["contribution"] <= 0:
            continue
        for key, count in zip(row["source_nclose_keys"], row["source_occurrence_counts"], strict=True):
            if count <= 0:
                continue
            for endpoints, indices in by_source.get(key, {}).items():
                call = grouped.setdefault(endpoints, dict(endpoints=endpoints, rows=[], names=set(),
                                                         nclose_ids=set(), keys=set(), classes=set()))
                occurrences = int(count) * len(indices)
                contribution = dict(row, source_nclose_keys=(key,), source_occurrence_counts=(int(count),),
                                    source_nclose_key=key, nclose_occurrence_count=int(count),
                                    bnd_occurrences_per_nclose=len(indices), junction_indices=tuple(indices),
                                    occurrence_count=occurrences,
                                    contribution=occurrences * row["structure_weight"],
                                    contribution_N=occurrences * row["structure_weight_N"])
                call["rows"].append(contribution)
                call["nclose_ids"].add(row["nclose_id"])
                call["keys"].add(key)
                call["names"].update(context.model["nclose_sources"][key].get("contig_names", ()))
                call["classes"].add({
                    "PATH": "NCLOSE", "AMP": "AMPLICON", "MERGE_TYPE4": "MERGED_TYPE4",
                    "VIRTUAL_INV": "VIRTUAL_INV", "TYPE4": "NCLOSE",
                }[row["kind"]])
    result = []
    for endpoints, call in sorted(grouped.items()):
        rows = call["rows"]
        weight = fsum(r["contribution_N"] for r in rows)
        if weight <= min_cn:
            continue
        parents = sorted({r["structure_id"] for r in rows})
        own = {r["structure_id"]: r["structure_weight_N"] for r in rows}
        model = [r for r in rows if r["feature_index"] is not None]
        virtual = [r for r in rows if r["feature_index"] is None]
        call.update(weight=weight, info=dict(
            SVCLASS=sorted(call["classes"]), PARENT_IDS=parents,
            PARENT_WEIGHTS=[own[p] for p in parents],
            PARENT_MULTIPLICITY=[sum(r["occurrence_count"] for r in rows if r["structure_id"] == p) for p in parents],
            NCLOSE_IDS=sorted(call["nclose_ids"]),
            NCLOSE_KEYS=[f"{a}:{b}" for a, b in sorted(call["keys"])],
            WEIGHT_METHOD="STRUCTURE_SUM",
            MODEL_WEIGHT=fsum(r["contribution_N"] for r in model),
            VIRTUAL_WEIGHT=fsum(r["contribution_N"] for r in virtual),
            MODEL_FEATURE_COUNT=len({r["feature_index"] for r in model}),
            MODEL_OCCURRENCE_COUNT=sum(r["occurrence_count"] for r in model)))
        result.append(call)
    return result


def bnd_calls(context, nodes, min_cn=0.1):
    return aggregate_bnd_calls(context, nclose_bnd_memberships(context, nodes), min_cn)
