"""Native output views of the shared structure/NClose accounting table."""
from collections import defaultdict
from math import fsum
from pathlib import Path
import csv

import vcfpy


def bnd_calls(context, min_cn=0.1):
    grouped = {}
    for pid, rows in context.projected_rows.items():
        event = context.projections[pid]
        if event["kind"] != "bnd":
            continue
        endpoints = tuple(sorted(event["endpoints"]))
        call = grouped.setdefault(endpoints, dict(endpoints=endpoints, rows={}, names=set(),
                                                 nclose_ids=set(), keys=set(), classes=set()))
        for row in rows:
            if row["contribution"] <= 0:
                continue
            call["rows"][(row["structure_id"], pid)] = row
            call["nclose_ids"].add(row["nclose_id"])
            call["keys"].add(row["nclose_key"])
            call["names"].update(event.get("contig_names", ()))
            call["classes"].add("PATH_SPLIT" if event["is_split"] else {
                "PATH": "NCLOSE", "AMP": "AMPLICON", "MERGE_TYPE4": "MERGED_TYPE4",
                "VIRTUAL_INV": "VIRTUAL_INV", "TYPE4": "NCLOSE",
            }[row["kind"]])
    result = []
    for endpoints, call in sorted(grouped.items()):
        rows = list(call["rows"].values())
        weight = fsum(r["contribution_N"] for r in rows)
        if weight <= min_cn:
            continue
        parents = sorted({r["structure_id"] for r in rows})
        own = {r["structure_id"]: r["structure_weight_N"] for r in rows}
        model = [r for r in rows if r["feature_index"] is not None]
        virtual = [r for r in rows if r["feature_index"] is None]
        call.update(weight=weight, rows=rows, info=dict(
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


def write_native_vcf(context, contig_lengths, path, ratios, build_header, write_bnd, write_symbolic):
    records, contributions = [], []
    class Buffer:
        def write_record(self, record):
            records.append(record)
    writer = Buffer()
    calls = bnd_calls(context)
    for index, call in enumerate(calls, 1):
        name = f"SKYPE.BND.{index}"
        a, b = call["endpoints"]
        ratio_info = ratios.pair([(c, p, "left" if side == "L" else "right")
                                  for c, p, side in call["endpoints"]], call["weight"] * context.n_unit)
        write_bnd(writer, name, a[0], a[1], "+" if a[2] == "L" else "-",
                  b[0], b[1], "+" if b[2] == "R" else "-", call["weight"],
                  "|".join(sorted(call["names"])), bp_ratio_info=ratio_info, extra_info=call["info"])
        for row in call["rows"]:
            contributions.append((name, row["feature_index"] if row["feature_index"] is not None else ".",
                                  row["projection_id"], row["occurrence_count"],
                                  row["structure_weight_N"], row["contribution_N"], "STRUCTURE_SUM",
                                  row["structure_id"], row["nclose_id"], row["kind"]))
    counts = defaultdict(int)
    for event, raw_weight in context.visible_projections():
        if event["kind"] != "indel":
            continue
        chrom, st, nd = event["chrom"], event["st"], event["nd"]
        svtype = "DEL" if event["event_type"] == "front_jump" else "DUP"
        counts[svtype] += 1
        sides = ("left", "right") if svtype == "DEL" else ("right", "left")
        ratio_info = ratios.pair([(chrom, st, sides[0]), (chrom, nd, sides[1])], raw_weight)
        write_symbolic(writer, chrom, st, f"SKYPE.{svtype}.{counts[svtype]}", svtype,
                       nd, (nd - st) * (-1 if svtype == "DEL" else 1),
                       raw_weight / context.n_unit, event.get("source", event["nclose_id"]),
                       bp_ratio_info=ratio_info,
                       extra_info={"NCLOSE_IDS": [event["nclose_id"]], "WEIGHT_METHOD": "STRUCTURE_SUM"})
    rank = {chrom: i for i, chrom in enumerate(contig_lengths)}
    records.sort(key=lambda r: (rank[r.CHROM], r.POS, r.ID[0]))
    with vcfpy.Writer.from_path(str(path), build_header(contig_lengths)) as output:
        for record in records:
            output.write_record(record)
    with Path(path).with_suffix(".bnd_weights.tsv").open("w") as handle:
        output = csv.writer(handle, delimiter="\t", lineterminator="\n")
        output.writerow(("bnd_id", "feature_index", "source_occurrence", "occurrence_count",
                         "feature_weight_N", "contribution_N", "weight_method", "structure_id", "nclose_id", "structure_kind"))
        output.writerows(contributions)
    return len(calls), len(records)


def display_events(context, nodes, centromeres, min_cn=0.1):
    events = []
    for projection, weight in context.visible_projections(min_cn):
        event = dict(weight=weight, weight_N=weight / context.n_unit,
                     source_id=projection["projection_id"], nclose_ids=projection["nclose_id"],
                     weight_scope="NCLOSE", kind="NCLOSE")
        if projection["kind"] == "indel":
            chrom, st, nd = projection["chrom"], projection["st"], projection["nd"]
            event.update(type="Deletion" if projection["event_type"] == "front_jump" else "Duplication",
                         spans=[(chrom, st, nd)], link=[(chrom, st), (chrom, nd)], link_type="indel")
        else:
            a, b = projection["endpoints"]
            event.update(type="Breakend", link_type="inversion" if a[0] == b[0] and a[2] == b[2] else "breakend")
            if projection["is_split"]:
                event.update(spans=[(c, max(0, p - 1), p) for c, p, _ in (a, b)],
                             link=[a[:2], b[:2]])
            else:
                left, right = (nodes[i] for i in projection["parent_key"])
                event.update(spans=[(n[5], n[7], n[8]) for n in (left, right)],
                             link=[(n[5], n[7]) for n in (left, right)])
        events.append(event)
    for structure in context.structure_summaries(min_cn):
        kind = structure["kind"]
        if kind == "MERGE_TYPE4":
            e = structure["event"]
            spans = [(e["chrom"], e["st"], e["nd"])]
            label = "Deletion" if e["event_type"] == "front_jump" else "Duplication"
            link_type = "indel"
        elif kind == "AMP":
            spans, label, link_type = [structure["span"]], "Amplicon", "amplicon"
        else:
            spans, label, link_type = structure["views"], "Virtual_inversion", "virtual_inv"
        ids = ";".join(context.ncloses[k]["nclose_id"] for k in structure["nclose_counts"])
        for chrom, st, nd in spans:
            events.append(dict(kind=kind, type=label, link_type=link_type,
                               spans=[(chrom, st, nd)], link=[(chrom, st), (chrom, nd)],
                               weight=structure["raw_weight"], weight_N=structure["weight_N"],
                               source_id=structure["structure_id"], nclose_ids=ids, weight_scope="STRUCTURE"))
    for structure in context.structures.values():
        if structure["kind"] != "CENTROMERE" or structure["weight_N"] <= min_cn:
            continue
        chrom = structure["chrom"]
        info = centromeres[chrom]
        st, nd = (info["mid"], info["chr_len"]) if info["dir"] else (0, info["mid"])
        events.append(dict(kind="CENTROMERE", type="Centromere", spans=[(chrom, st, nd)],
                           weight=structure["raw_weight"], weight_N=structure["weight_N"],
                           source_id=structure["structure_id"], nclose_ids=".", weight_scope="STRUCTURE"))
    return events


def write_native_bed(prefix, events):
    with (Path(prefix) / "SKYPE_result.bed").open("w") as handle:
        output = csv.writer(handle, delimiter="\t", lineterminator="\n")
        output.writerow(("#chrom", "cordst", "cordnd", "type", "weight (N)", "nclose_id", "weight_scope", "source_id"))
        for event in events:
            for chrom, st, nd in event["spans"]:
                output.writerow((chrom, st, nd, event["type"], event["weight_N"],
                                 event["nclose_ids"], event["weight_scope"], event["source_id"]))
