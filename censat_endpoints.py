"""Independent CEN-SAT endpoint workflow and its reusable alignment inputs.

Partition on raw alignasm terminal bases, then require chromosome/strand
consistency on each original end chunk in both raw and cropped P+S alignments.
The CLI prepares sample/reference artifacts; stage 01 evaluates and imports
the accepted unitigs without passing them through legacy NClose filters.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import json
import logging
from pathlib import Path
import subprocess


ALIGN_OPTIONS = ["--cs", "-x", "asm20", "--no-long-join", "-r2k", "-K10G",
                 "-N", "5000", "-p", "0.5"]


def signature(path):
    path = Path(path).resolve()
    stat = path.stat()
    return dict(path=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns)


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + "\n")


def write_tsv(path, records, columns):
    with Path(path).open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def read_bed(path):
    """All BED intervals, without the legacy 1 Mb size cutoff."""
    intervals = defaultdict(list)
    for line in Path(path).open():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        chrom, start, end = line.split()[:3]
        intervals[chrom].append((int(start), int(end)))
    return intervals


def read_paf(path, names=None):
    """Keep original zero-based row indices, including when selecting queries."""
    with Path(path).open() as handle:
        for index, line in enumerate(handle):
            if names is not None and line.split("\t", 1)[0] not in names:
                continue
            fields = line.rstrip("\n").split("\t")
            tags = dict(tag.split(":", 2)[::2] for tag in fields[12:])
            yield dict(name=fields[0], qlen=int(fields[1]), qs=int(fields[2]),
                       qe=int(fields[3]), strand=fields[4], chrom=fields[5],
                       rlen=int(fields[6]), rs=int(fields[7]), re=int(fields[8]),
                       mapq=int(fields[11]), index=index, xi=tags.get("xi"),
                       tp=tags.get("tp"))


def terminal_position(row, side):
    use_start = (side == "left") == (row["strand"] == "+")
    return row["rs"] if use_start else row["re"] - 1


def partition_unitigs(paf, bed):
    intervals = read_bed(bed)
    ends = {}
    for row in read_paf(paf):
        name = row["name"]
        if name not in ends:
            ends[name] = dict(left=row, right=row)
        pair = ends[name]
        left_key = lambda r: (r["qs"], r["qe"], r["index"])
        right_key = lambda r: (r["qe"], r["qs"], -r["index"])
        pair["left"] = min(pair["left"], row, key=left_key)
        pair["right"] = max(pair["right"], row, key=right_key)
    partition, candidates = [], []
    for name, pair in ends.items():
        both = all(any(start <= terminal_position(row, side) < end
                       for start, end in intervals[row["chrom"]])
                   for side, row in pair.items())
        left, right = pair["left"], pair["right"]
        same = (left["chrom"], left["strand"]) == (right["chrom"], right["strand"])
        reason = ("same_chrom_same_strand" if same else "candidate") if both else "legacy"
        partition.append(dict(unitig=name, route="censat" if both else "legacy",
                              reason=reason, left_index=left["index"],
                              right_index=right["index"]))
        if both and not same:
            candidates.append(dict(unitig=name, **pair))
    return partition, candidates


def consistency(chunk, alignments, shift=0):
    hits = [row for row in alignments if row["tp"] in ("P", "S")
            and 2 * max(0, min(chunk["qe"], row["qe"] + shift)
                        - max(chunk["qs"], row["qs"] + shift))
            >= chunk["qe"] - chunk["qs"]]
    conflicts = sum((row["chrom"], row["strand"]) !=
                    (chunk["chrom"], chunk["strand"]) for row in hits)
    return dict(status="conflict" if conflicts else "consistent" if hits else "no_alignment",
                hits=len(hits), conflicts=conflicts)


def pair_status(results):
    states = {result["status"] for result in results}
    return "conflict" if "conflict" in states else "no_alignment" if "no_alignment" in states else "consistent"


def trace_source(chunk, sources):
    xi = chunk.get("xi") or ""
    if not xi.startswith("P_"):
        return None, "unsupported_xi" if xi else "missing_xi"
    try:
        source = sources.get(int(xi[2:]))
    except ValueError:
        source = None
    if source is None:
        return None, "missing_source"
    if (source["name"], source["chrom"], source["strand"]) != (
            chunk["name"], chunk["chrom"], chunk["strand"]) or not (
            source["qs"] <= chunk["qs"] < chunk["qe"] <= source["qe"]):
        return None, "source_mismatch"
    return source, None


def prepare_inputs(aln_paf, raw_paf, fasta, reference, bed, outdir, thread=1, force=False):
    """Cache expensive preparation independently of a SKYPE result directory."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    metadata = dict(version=1, inputs={key: signature(path) for key, path in
                    dict(aln_paf=aln_paf, raw_paf=raw_paf, fasta=fasta,
                         reference=reference, bed=bed).items()},
                    implementation=signature(__file__), options=ALIGN_OPTIONS,
                    minimap2=subprocess.check_output(["minimap2", "--version"], text=True).strip())
    meta_path = outdir / "cache_metadata.json"
    artifacts = ("partition.tsv", "candidates.json", "manifest.tsv",
                 "source_end_regions.fa", "realigned.paf")
    if not force and meta_path.exists() and all((outdir / name).exists() for name in artifacts):
        if json.loads(meta_path.read_text()) == metadata:
            logging.info("Reusing CEN-SAT endpoint alignments")
            return outdir
    # An interrupted rebuild must never leave an apparently current cache.
    meta_path.unlink(missing_ok=True)
    partition, candidates = partition_unitigs(aln_paf, bed)
    write_tsv(outdir / "partition.tsv", partition,
              ["unitig", "route", "reason", "left_index", "right_index"])
    names = {candidate["unitig"] for candidate in candidates}
    raw = defaultdict(list)
    sources = {}
    for row in read_paf(raw_paf, names):
        raw[row["name"]].append(row)
        sources[row["index"]] = row
    manifest = []
    for candidate in candidates:
        results = [consistency(candidate[side], raw[candidate["unitig"]])
                   for side in ("left", "right")]
        candidate["raw_results"] = results
        candidate["raw_status"] = pair_status(results)
        candidate["trace_reason"] = ""
        if candidate["raw_status"] != "consistent":
            continue
        traced = [trace_source(candidate[side], sources) for side in ("left", "right")]
        errors = [reason for _, reason in traced if reason]
        if errors:
            candidate["trace_reason"] = ";".join(errors)
            continue
        for side, (source, _) in zip(("left", "right"), traced):
            chunk = candidate[side]
            query_id = f'{candidate["unitig"]}__{side}__{chunk["xi"]}'
            manifest.append(dict(query_id=query_id, unitig=candidate["unitig"],
                                 side=side, xi=chunk["xi"], source_qstart=source["qs"],
                                 source_qend=source["qe"], chunk_qstart=chunk["qs"],
                                 chunk_qend=chunk["qe"], chrom=chunk["chrom"],
                                 strand=chunk["strand"]))
    write_json(outdir / "candidates.json", candidates)
    write_tsv(outdir / "manifest.tsv", manifest,
              ["query_id", "unitig", "side", "xi", "source_qstart", "source_qend",
               "chunk_qstart", "chunk_qend", "chrom", "strand"])
    fa_path = outdir / "source_end_regions.fa"
    paf_path = outdir / "realigned.paf"
    with fa_path.open("w") as output:
        if manifest:
            fai = Path(str(fasta) + ".fai")
            if not fai.exists() or fai.stat().st_mtime_ns < Path(fasta).stat().st_mtime_ns:
                subprocess.run(["samtools", "faidx", str(fasta)], check=True)
            regions = [f'{row["unitig"]}:{row["source_qstart"] + 1}-{row["source_qend"]}'
                       for row in manifest]
            region_path = outdir / "regions.txt"
            region_path.write_text("\n".join(regions) + "\n")
            extracted = subprocess.check_output(
                ["samtools", "faidx", str(fasta), "-r", str(region_path)], text=True)
            records = extracted.split(">")[1:]
            if len(records) != len(manifest):
                raise ValueError("CEN-SAT source extraction returned an unexpected number of records")
            for row, record in zip(manifest, records):
                _, sequence = record.split("\n", 1)
                output.write(f'>{row["query_id"]}\n{sequence}')
    logging.info(
        "CEN-SAT preparation: %d base candidates, %d unitigs for realignment",
        len(candidates), len(manifest) // 2,
    )
    if manifest:
        subprocess.run(["minimap2", *ALIGN_OPTIONS, "-t", str(thread),
                        str(reference), str(fa_path), "-o", str(paf_path)], check=True)
    else:
        paf_path.write_text("")
    write_json(meta_path, metadata)
    return outdir


def evaluate_inputs(outdir, aln_paf, raw_paf, bed):
    """Read prepared inputs in stage 01; no realignment or legacy filtering."""
    outdir = Path(outdir)
    metadata = json.loads((outdir / "cache_metadata.json").read_text())
    for key, path in dict(aln_paf=aln_paf, raw_paf=raw_paf, bed=bed).items():
        if metadata["inputs"][key] != signature(path):
            raise ValueError(f"Stale CEN-SAT endpoint cache ({key}); rebuild through SKYPE.py")
    with (outdir / "partition.tsv").open() as handle:
        partition = list(csv.DictReader(handle, delimiter="\t"))
    candidates = json.loads((outdir / "candidates.json").read_text())
    with (outdir / "manifest.tsv").open() as handle:
        manifest = {(row["unitig"], row["side"]): row
                    for row in csv.DictReader(handle, delimiter="\t")}
    realigned = defaultdict(list)
    for row in read_paf(outdir / "realigned.paf"):
        realigned[row["name"]].append(row)
    accepted, diagnostics = [], []
    for candidate in candidates:
        results = []
        if candidate["raw_status"] != "consistent":
            status = "raw_" + candidate["raw_status"]
        elif candidate["trace_reason"]:
            status = candidate["trace_reason"]
        else:
            for side in ("left", "right"):
                row = manifest[candidate["unitig"], side]
                results.append(consistency(candidate[side], realigned[row["query_id"]],
                                           int(row["source_qstart"])))
            status = "realigned_" + pair_status(results)
            if status == "realigned_consistent":
                accepted.append(candidate)
        diagnostics.append(dict(unitig=candidate["unitig"], status=status,
                                left_hits=results[0]["hits"] if results else "",
                                right_hits=results[1]["hits"] if results else "",
                                left_conflicts=results[0]["conflicts"] if results else "",
                                right_conflicts=results[1]["conflicts"] if results else ""))
    excluded = {row["unitig"] for row in partition if row["route"] == "censat"}
    summary = dict(cache_dir=str(outdir.resolve()), total_unitigs=len(partition),
                   censat_unitigs=len(excluded), legacy_unitigs=len(partition) - len(excluded),
                   base_candidates=len(candidates),
                   raw_statuses=dict(Counter(c["raw_status"] for c in candidates)),
                   realigned_unitigs=len(manifest) // 2, accepted_unitigs=len(accepted),
                   statuses=dict(Counter(row["status"] for row in diagnostics)))
    return excluded, accepted, diagnostics, summary


def main():
    import skype_utils  # Apply the shared pipeline logging configuration.

    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("aln-paf", "raw-paf", "fasta", "reference", "bed", "outdir"):
        parser.add_argument("--" + name, required=True)
    parser.add_argument("-t", "--thread", type=int, default=1)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    prepare_inputs(**vars(args))


if __name__ == "__main__":
    main()
