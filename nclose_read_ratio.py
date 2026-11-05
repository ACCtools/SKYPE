"""Sequence-resolved REF/ALT read ratios for assembly NCloses.

An NClose stands for every unitig junction that stage 01 compressed into it
(both endpoints within the compression limit, same directions), so each member
unitig is measured on its own.  At each NClose end the unitig leaves the
reference trajectory of its anchor.  The unitig sequence around that point
(retained flank, then the unitig continuation) is the ALT allele; the
reference across the same breakpoint is the REF allele.  ALT sequences that
are nearly identical across the breakpoint (other members, repeat copies at
other coordinates, both ends of a direct junction) form one allele class.
Reads fetched at the breakpoint are aligned to every REF/ALT sequence at once
and assigned to the class whose sequence explains the read around the
breakpoint best by a clear score margin.

ALT sequences come from the sample's own assembly and REF from the reference
genome, so the sample's divergence from the reference favours ALT wherever the
two alleles are the same sequence.  Where the ALT junction itself occurs in the
reference near the breakpoint at near-identity -- the continuation stays
collinear with the reference past the breakpoint (inverted repeats, segmental
duplications), or another copy within the reach of the reads fetched there
carries the same sequence across the junction (tandem and nearby repeats) --
reads are assigned only once they cover sequence beyond that homology; if the
whole ALT continuation is homologous, the end is not testable.

Reads are fetched by their position in the read alignment (BAM).  When the
unitig's anchor on the retained side is shorter than FLANK, the junction
molecule holds only a short piece of that locus (typically a templated
insertion); aligners place such reads at the host locus with the piece as an
insertion, so the reads carrying the junction are mostly absent at this
position.  Such an end would undercount ALT, so it takes part in read
assignment but not in the estimate.

The retained flank of a breakpoint is present in CN_local copies (N units):
the junction copies J, the reference-adjacency copies and the copies of any
other junction at the same point.  CN_local is the normalized depth of the
nearest 100 kb windows inside the retained segment; when another breakpoint
bounds that segment before a whole window fits (insertions, complex events),
the windows would measure other sequence, and the read depth at the breakpoint
itself, calibrated against the windows at the sample's unbounded sites, is used.  Reads crossing the breakpoint sample those
copies in proportion, so J = ALT / (ALT + REF + OTHER) * CN_local.  Capture
efficiency and local coverage cancel in the ratio.  An end is not testable
when its retained flank is not uniquely placed, when ALT and REF are not
distinguishable by sequence, or when too few reads cross it; such NCloses are
kept and reported as NOT_TESTABLE.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import logging
import math
import os
import pickle
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pysam
from scipy.integrate import quad
from scipy.optimize import linprog
from scipy.special import pdtr

import depth_scale

CTG_NAM = 0
CTG_LEN = 1
CTG_STR = 2
CTG_END = 3
CTG_DIR = 4
CHR_NAM = 5
CHR_STR = 7
CHR_END = 8
CTG_MAPQ = 9
CTG_STRND = 11
CTG_ENDND = 12
CTG_GLOBALIDX = 21

FLANK = 5_000                  # ALT/REF sequence on each side of a breakpoint
MIN_SIDE = 200                 # shortest retained flank / continuation tested
MAX_REQUIRED_SIDE = 1_000      # crossing reads must cover this much per side
REQUIRED_SIDE_FRACTION = 0.8   # ... or this fraction of a shorter flank
FETCH_PAD = 1_000
MIN_ANCHOR_MAPQ = 20
MIN_READ_MAPQ = 20
MIN_UNIQUE_FRACTION = 0.8
MAX_SITE_READS = 400           # deterministic subsample of deep loci
SCORE_MARGIN_BASES = 50        # winner must beat every other class by this
OTHER_STOP_TOLERANCE = 300
SITE_TOLERANCE = 1_000         # ALT classes compete at breakpoints this close
CLUSTER_IDENTITY = 0.98
CLUSTER_BP_TOLERANCE = 200
HOMOLOGY_IDENTITY = 0.97       # ALT continuation this similar to REF is not discriminating
MAX_QUERY_GAP = 50             # reference-built ALT: largest unaligned query gap between pieces
SV_INDEL = 50                  # an indel this long on a read is a structural variant, not error
CN_WINDOWS = 3
CN_SEARCH = 300_000
MIN_LOCAL_CN = 0.1
MIN_SITE_READS = 3
MAX_AMBIGUOUS_FRACTION = 0.5
DISPLAY_N = 0.2                # SKYPE display floor and filter effect size
MIN_EXPECTED_ALT_AT_DISPLAY = 3.0
POSTHOC_TAIL = 0.01
MIN_SCALE_SD = 0.1             # bounds on the per-sample reads-per-copy scatter (log scale)
MAX_SCALE_SD = 1.0

MATCH_SCORE = {"map-hifi": 1, "lr:hq": 1, "map-ont": 2, "map-pb": 2}
OUTPUT_DIR = "nclose_read_ratio"
MEMBERS_JSON = "members.json"
NCLOSE_TSV = "nclose_read_ratio.tsv"
POSTHOC_TSV = "nclose_read_ratio_posthoc.tsv"
REMOVE_TSV = "posthoc_remove.tsv"
PREFIT_REMOVE_TSV = "prefit_remove.tsv"

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq):
    return seq.translate(_COMPLEMENT)[::-1]


def identity_key(identity):
    contig, pair = identity
    return f"{contig}:{int(pair[0])}:{int(pair[1])}"


def parse_identity_key(text):
    contig, a, b = text.rsplit(":", 2)
    return contig, (int(a), int(b))


@dataclass
class MemberEnd:
    """One NClose end of one member unitig and its ALT sequence."""

    end_id: int
    nclose: tuple
    member: tuple
    end: str
    row: int
    chrom: str
    bp: int
    retained: str
    anchor_mapq: int
    seq_source: str = "none"
    alt_seq: str = ""
    alt_bp: int = 0
    cont_len: int = 0
    reason: str = ""
    cls: str = ""
    site: int = -1
    discriminating: bool = True
    homology: int = 0
    req_cont: int = 0


@dataclass
class Site:
    """A reference breakpoint at which crossing reads are counted."""

    site_id: int
    chrom: str
    bp: int
    retained: str
    req_ret: int
    req_cont: int
    ref_target: str = ""
    ref_offset: int = 0
    ctx_target: str = ""
    reach_lo: int = -1
    reach_hi: int = -1
    ctx_lo: int = 0
    ctx_hi: int = 0
    local_targets: set = field(default_factory=set)
    reads: list = field(default_factory=list)
    sampled_fraction: float = 1.0
    unique_fraction: float = float("nan")
    mapq_records: int = 0
    cn_local: float = float("nan")
    cn_source: str = ""
    counts: dict = field(default_factory=dict)
    alt_counts: dict = field(default_factory=dict)
    alt_names: dict = field(default_factory=dict)


def anchor_breakpoint(row, end):
    """Reference breakpoint and retained side where the unitig leaves ``row``.

    End A leaves the anchor at its query end, end B at its query start.
    """
    forward = row[CTG_DIR] == "+"
    if end == "A":
        return (int(row[CHR_END]), "left") if forward else (int(row[CHR_STR]), "right")
    return (int(row[CHR_STR]), "right") if forward else (int(row[CHR_END]), "left")


def required_side(length):
    return int(min(MAX_REQUIRED_SIDE, REQUIRED_SIDE_FRACTION * length))


class _Sequences:
    """Unitig FASTA plus named in-memory sequences (raw-rescue template reads)."""

    def __init__(self, fasta, extra):
        self.fasta = pysam.FastaFile(fasta) if fasta else None
        self.extra = dict(extra or {})
        self.references = set(self.fasta.references if self.fasta else ()) | set(self.extra)

    def get_reference_length(self, name):
        if name in self.extra:
            return len(self.extra[name])
        return self.fasta.get_reference_length(name)

    def fetch(self, name, start, end):
        if name in self.extra:
            return self.extra[name][start:end]
        return self.fasta.fetch(name, start, end)


def load_rescue_templates(prefix):
    """Template sequence of every stage-24 rescue owner (its chain's query).

    Read-method chains are raw reads (raw_realign.fa); OLC-method chains are
    unitigs assembled from them (olc_realign.fa).
    """
    rescue = Path(prefix) / "24_raw_rescue"
    candidates_path = rescue / "raw_candidates.json"
    if not candidates_path.exists():
        return {}
    with open(candidates_path) as handle:
        candidates = json.load(handle)
    owner_read = {e["owner"]: e["name"] for e in candidates if e.get("owner") and e.get("name")}
    wanted = set(owner_read.values())
    reads = {}
    for fasta in ("raw_realign.fa", "olc_realign.fa"):
        if not (rescue / fasta).exists():
            continue
        with pysam.FastxFile(str(rescue / fasta)) as handle:
            for entry in handle:
                if entry.name in wanted and entry.name not in reads:
                    reads[entry.name] = entry.sequence.upper()
    missing = sorted(name for name in wanted if name not in reads)
    if missing:
        logging.warning("No template sequence for %d rescue chain(s), e.g. %s", len(missing), missing[0])
    return {owner: reads[name] for owner, name in owner_read.items() if name in reads}


class SequenceSource:
    """Unitig sequence by name, with PAF global-index fallback for renames."""

    def __init__(self, unitig_fasta, reference_fasta, paf_paths, extra_sequences=None):
        self.unitigs = _Sequences(unitig_fasta, extra_sequences) if (unitig_fasta or extra_sequences) else None
        self.reference = pysam.FastaFile(reference_fasta)
        self.paf_paths = list(paf_paths or ())
        self._paf_names = {}

    def _paf_name(self, global_idx):
        try:
            source, row = str(global_idx).split(".", 1)
            source, row = int(source), int(row)
        except ValueError:
            return None
        if source >= len(self.paf_paths):
            return None
        if source not in self._paf_names:
            names = []
            with open(self.paf_paths[source]) as handle:
                for line in handle:
                    names.append(line.split("\t", 1)[0])
            self._paf_names[source] = names
        names = self._paf_names[source]
        return names[row] if row < len(names) else None

    def unitig_name(self, row):
        if self.unitigs is None:
            return None
        names = self.unitigs.references
        name = str(row[CTG_NAM])
        if name in names:
            return name
        original = self._paf_name(row[CTG_GLOBALIDX])
        return original if original in names else None

    def oriented_reference(self, row, start=None, end=None):
        """Reference sequence of ``row`` in query orientation, sliced by query offset."""
        seq = self.reference.fetch(str(row[CHR_NAM]), int(row[CHR_STR]), int(row[CHR_END])).upper()
        if row[CTG_DIR] == "-":
            seq = revcomp(seq)
        return seq[start:end]


def _unitig_alt(source, name, row, end):
    length = source.unitigs.get_reference_length(name)
    if end == "A":
        q = int(row[CTG_END])
        ret = min(FLANK, q - int(row[CTG_STR]))
        cont = min(FLANK, length - q)
        if ret <= 0 or cont <= 0:
            return "", ret, cont
        return source.unitigs.fetch(name, q - ret, q + cont).upper(), ret, cont
    q = int(row[CTG_STR])
    ret = min(FLANK, int(row[CTG_END]) - q)
    cont = min(FLANK, q)
    if ret <= 0 or cont <= 0:
        return "", ret, cont
    return revcomp(source.unitigs.fetch(name, q - cont, q + ret).upper()), ret, cont


def _reference_alt(source, contig_data, pair, end):
    """ALT from the reference pieces between the anchors (no unitig sequence).

    Only where consecutive pieces are adjacent on the query: an unaligned stretch
    between them (inserted sequence) is unknown, and a join without it would be
    a junction the molecule does not have.
    """
    first, last = pair
    order = range(first, last + 1) if end == "A" else range(last, first - 1, -1)
    used, cont_len = [order[0]], 0
    for idx in order[1:]:
        used.append(idx)
        cont_len += int(contig_data[idx][CTG_END]) - int(contig_data[idx][CTG_STR])
        if cont_len >= FLANK:
            break
    for a, b in zip(sorted(used), sorted(used)[1:]):
        if abs(int(contig_data[b][CTG_STR]) - int(contig_data[a][CTG_END])) > MAX_QUERY_GAP:
            return "", 0, 0
    if end == "A":
        retained = source.oriented_reference(contig_data[first])[-FLANK:]
        cont = "".join(source.oriented_reference(contig_data[idx]) for idx in used[1:])[:FLANK]
        return retained + cont, len(retained), len(cont)
    retained = source.oriented_reference(contig_data[last])[:FLANK]
    cont = "".join(source.oriented_reference(contig_data[idx]) for idx in reversed(used[1:]))[-FLANK:]
    return revcomp(cont + retained), len(retained), len(cont)


def build_member_ends(contig_data, members_by_nclose, source):
    ends = []
    for nclose, members in members_by_nclose.items():
        for member in members:
            contig, pair = member
            pair = tuple(int(x) for x in pair)
            for end, row_idx in (("A", pair[0]), ("B", pair[1])):
                row = contig_data[row_idx]
                bp, retained = anchor_breakpoint(row, end)
                item = MemberEnd(
                    end_id=len(ends), nclose=nclose, member=(contig, pair), end=end,
                    row=row_idx, chrom=str(row[CHR_NAM]), bp=bp, retained=retained,
                    anchor_mapq=int(row[CTG_MAPQ]),
                )
                name = source.unitig_name(row)
                try:
                    if name is not None:
                        seq, ret, cont = _unitig_alt(source, name, row, end)
                        item.seq_source = "unitig"
                    else:
                        seq, ret, cont = _reference_alt(source, contig_data, pair, end)
                        item.seq_source = "reference"
                except (KeyError, ValueError, IndexError) as exc:
                    logging.warning("No ALT sequence for %s end %s: %s", contig, end, exc)
                    seq, ret, cont = "", 0, 0
                item.alt_seq, item.alt_bp, item.cont_len = seq, ret, cont
                if not seq:
                    item.reason = "no_sequence"
                elif ret < MIN_SIDE:
                    item.reason = "short_retained_flank"
                elif cont < MIN_SIDE:
                    item.reason = "short_continuation"
                ends.append(item)
    return ends


def build_sites(ends):
    """One site per distinct (chrom, bp, retained side) of a usable member end."""
    site_by_key = {}
    sites = []
    for item in ends:
        if item.reason:
            continue
        key = (item.chrom, item.bp, item.retained)
        site = site_by_key.get(key)
        if site is None:
            site = Site(len(sites), item.chrom, item.bp, item.retained,
                        required_side(item.alt_bp), required_side(item.cont_len))
            site_by_key[key] = site
            sites.append(site)
        site.req_ret = min(site.req_ret, required_side(item.alt_bp))
        site.req_cont = min(site.req_cont, required_side(item.cont_len))
        item.site = site.site_id
    return sites


def _merged_windows(sites, reference, window, label):
    """Reference windows ``window(site)`` merged per chromosome; (targets, site -> (name, start))."""
    by_chrom = defaultdict(list)
    for site in sites:
        length = reference.get_reference_length(site.chrom)
        lo, hi = window(site)
        by_chrom[site.chrom].append((max(0, lo), min(length, hi), site))
    targets, placement = {}, {}
    for chrom, items in by_chrom.items():
        items.sort(key=lambda x: x[0])
        groups = []
        for st, nd, site in items:
            if groups and st <= groups[-1][1]:
                groups[-1][1] = max(groups[-1][1], nd)
                groups[-1][2].append(site)
            else:
                groups.append([st, nd, [site]])
        for st, nd, members in groups:
            name = f"{label}|{chrom}|{st}|{nd}"
            targets[name] = reference.fetch(chrom, st, nd).upper()
            for site in members:
                placement[site.site_id] = (name, st)
    return targets, placement


def build_ref_targets(sites, reference):
    """Merge breakpoint windows per chromosome so REF targets never overlap."""
    targets, placement = _merged_windows(sites, reference, lambda s: (s.bp - FLANK, s.bp + FLANK), "REF")
    for site in sites:
        site.ref_target, start = placement[site.site_id]
        site.ref_offset = site.bp - start
    return targets


def build_context_targets(sites, reference):
    """Reference around each site as far as the reads fetched there reach.

    A reference read fetched at a breakpoint comes from this stretch, so a copy
    of the ALT junction inside it would let such reads pass for junction reads.
    """
    def window(site):
        lo = min(site.bp, site.reach_lo) if site.reach_lo >= 0 else site.bp
        hi = max(site.bp, site.reach_hi)
        return lo - FLANK, hi + FLANK

    targets, placement = _merged_windows(sites, reference, window, "CTX")
    for site in sites:
        site.ctx_target = placement[site.site_id][0]
        lo, hi = window(site)
        site.ctx_lo, site.ctx_hi = max(0, lo), min(reference.get_reference_length(site.chrom), hi)
    return targets


def _within_reach(aln, site):
    """The alignment lies in this site's own read reach (merged targets hold other sites' reach too)."""
    start = int(aln.target.rsplit("|", 2)[1])
    return site.ctx_lo <= start + aln.ts and start + aln.te <= site.ctx_hi


def _hash_rank(name):
    return hashlib.md5(name.encode()).hexdigest()


def fetch_site_reads(bam_path, sites):
    """Collect read names/sequences crossing each site and local mappability."""
    sequences = {}
    need_primary = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for site in sites:
            names = set()
            probe = site.bp - 100 if site.retained == "left" else site.bp + 100
            mapq_ok = mapq_all = 0
            for rec in bam.fetch(site.chrom, max(0, site.bp - FETCH_PAD), site.bp + FETCH_PAD):
                if rec.is_secondary or rec.is_qcfail or rec.is_duplicate or rec.is_unmapped:
                    continue
                site.reach_lo = rec.reference_start if site.reach_lo < 0 else min(site.reach_lo, rec.reference_start)
                site.reach_hi = max(site.reach_hi, rec.reference_end)
                if rec.reference_start <= probe < rec.reference_end:
                    mapq_all += 1
                    mapq_ok += rec.mapping_quality >= MIN_READ_MAPQ
                name = rec.query_name
                names.add(name)
                if not rec.is_supplementary:
                    if name not in sequences and rec.query_sequence:
                        sequences[name] = rec.query_sequence
                elif name not in sequences and rec.has_tag("SA"):
                    need_primary.setdefault(name, rec.get_tag("SA"))
            site.mapq_records = mapq_all
            site.unique_fraction = mapq_ok / mapq_all if mapq_all else float("nan")
            ranked = sorted(names, key=_hash_rank)
            site.sampled_fraction = (min(len(ranked), MAX_SITE_READS) / len(ranked)) if ranked else 1.0
            site.reads = ranked[:MAX_SITE_READS]
        wanted = {name for site in sites for name in site.reads}
        for name, sa in need_primary.items():
            if name in sequences or name not in wanted:
                continue
            for entry in sa.rstrip(";").split(";"):
                chrom, pos = entry.split(",")[:2]
                pos = int(pos) - 1
                for rec in bam.fetch(chrom, pos, pos + 1):
                    if rec.query_name == name and not rec.is_supplementary and not rec.is_secondary:
                        sequences[name] = rec.query_sequence
                        break
                if name in sequences:
                    break
    for site in sites:
        site.reads = [name for name in site.reads if name in sequences]
    return sequences


def detect_preset(bam_path):
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for pg in bam.header.to_dict().get("PG", []):
            # -x map-ont, -ax map-ont, -cx lr:hq ...
            match = re.search(r"(?:^|\s)-[A-Za-z]*x\s*(\S+)", pg.get("CL", ""))
            if match and match.group(1) in MATCH_SCORE:
                return match.group(1)
    return "map-hifi"


@dataclass
class Alignment:
    target: str
    strand: str
    qlen: int
    qs: int
    qe: int
    tlen: int
    ts: int
    te: int
    matches: int
    alnlen: int
    score: float
    cigar: str = ""
    # Read alignments: (target pos, query offset in target order, "I"/"D", length) of indels >= SV_INDEL.
    indels: tuple = ()


def run_minimap2(targets_fa, queries_fa, out_paf, preset, threads):
    cmd = ["minimap2", "-c", "--eqx", "-x", preset, "--secondary=yes", "-N", "100", "-p", "0.1",
           "-t", str(threads), targets_fa, queries_fa]
    with open(out_paf, "w") as out:
        subprocess.run(cmd, stdout=out, stderr=subprocess.DEVNULL, check=True)


def read_indels(aln, cigar):
    """Structural indels (SV_INDEL or longer) of a read alignment, for ``Alignment.indels``."""
    out, t, q = [], aln.ts, 0
    for n, op in re.findall(r"(\d+)([MIDNX=])", cigar):
        n = int(n)
        if op in "IDN" and n >= SV_INDEL:
            out.append((t, q, "I" if op == "I" else "D", n))
        if op in "MX=":
            t, q = t + n, q + n
        elif op == "I":
            q += n
        else:
            t += n
    return tuple(out)


def _indel_at(indel, pos):
    t, _, op, n = indel
    if op == "I":
        return abs(t - pos) <= OTHER_STOP_TOLERANCE
    return t - OTHER_STOP_TOLERANCE <= pos <= t + n + OTHER_STOP_TOLERANCE


def cut_at(aln, pos):
    """``aln`` cut at its structural indels at target position ``pos``.

    A read with an insertion or deletion of SV_INDEL or more at a breakpoint
    (or a deletion over it) does not carry that sequence's adjacency there: it
    leaves it, as a split-read aligner reports it, so it never passes for a read
    crossing the breakpoint contiguously.  Indels elsewhere on the read (germline
    variants near the breakpoint) stay inside a piece.  Each piece keeps its
    share of the score.
    """
    cuts = [x for x in aln.indels if _indel_at(x, pos)]
    if not cuts:
        return [aln]
    bounds = []
    t0, q0 = aln.ts, 0
    for t, q, op, n in cuts:
        bounds.append((t0, t, q0, q))
        t0, q0 = (t, q + n) if op == "I" else (t + n, q)
    bounds.append((t0, aln.te, q0, aln.qe - aln.qs))
    bounds = [b for b in bounds if b[1] > b[0] and b[3] > b[2]]
    aligned = sum(q1 - q0 for _, _, q0, q1 in bounds)
    pieces = []
    for t0, t1, q0, q1 in bounds:
        share = (q1 - q0) / aligned
        qs, qe = (aln.qs + q0, aln.qs + q1) if aln.strand == "+" else (aln.qe - q1, aln.qe - q0)
        rest = tuple((t, q - q0, op, n) for t, q, op, n in aln.indels if t0 <= t < t1 and not _indel_at((t, q, op, n), pos))
        pieces.append(Alignment(aln.target, aln.strand, aln.qlen, qs, qe, aln.tlen, t0, t1,
                                round(aln.matches * share), round(aln.alnlen * share), aln.score * share,
                                indels=rest))
    return pieces


def parse_paf(path):
    by_query = defaultdict(list)
    with open(path) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            score = next((int(x[5:]) for x in f[12:] if x.startswith("AS:i:")), 0)
            cigar = next((x[5:] for x in f[12:] if x.startswith("cg:Z:")), "")
            aln = Alignment(f[5], f[4], int(f[1]), int(f[2]), int(f[3]), int(f[6]), int(f[7]),
                            int(f[8]), int(f[9]), int(f[10]), score)
            if f[0].startswith("ALT|"):
                # ALT-vs-reference comparisons read the CIGAR itself.
                aln.cigar = cigar
                by_query[f[0]].append(aln)
            else:
                # Read CIGARs are not kept, only their structural indels.
                aln.indels = read_indels(aln, cigar)
                by_query[f[0]].append(aln)
    return by_query


def alt_target_name(item):
    return f"ALT|{item.end_id}"


def cluster_alt_targets(ends, alignments):
    """Union ALT sequences that are nearly identical across both sides of the breakpoint."""
    by_name = {alt_target_name(e): e for e in ends if e.alt_seq and not e.reason}
    parent = {name: name for name in by_name}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for qname, q in by_name.items():
        for aln in alignments.get(qname, ()):
            if aln.target == qname or aln.target not in by_name:
                continue
            t = by_name[aln.target]
            if aln.alnlen == 0 or aln.matches / aln.alnlen < CLUSTER_IDENTITY:
                continue
            ret = min(required_side(q.alt_bp), required_side(t.alt_bp))
            cont = min(required_side(q.cont_len), required_side(t.cont_len))
            if not (aln.qs <= q.alt_bp - ret and aln.qe >= q.alt_bp + cont):
                continue
            if not (aln.ts <= t.alt_bp - ret and aln.te >= t.alt_bp + cont):
                continue
            offset = q.alt_bp - aln.qs
            projected = aln.ts + offset if aln.strand == "+" else aln.te - offset
            if abs(projected - t.alt_bp) > CLUSTER_BP_TOLERANCE:
                continue
            parent[find(qname)] = find(aln.target)
    roots = {}
    for name, item in by_name.items():
        root = find(name)
        roots.setdefault(root, f"C{len(roots)}")
        item.cls = roots[root]


def mark_discriminating(ends, sites, alignments, match_score):
    """ALT is distinguishable only if its own sequence beats REF by the margin."""
    margin = SCORE_MARGIN_BASES * match_score
    for item in ends:
        if item.reason or item.site < 0:
            continue
        site = sites[item.site]
        alns = alignments.get(alt_target_name(item), ())
        self_score = max((a.score for a in alns if a.target == alt_target_name(item)),
                         default=len(item.alt_seq) * match_score)
        ref_score = max((a.score for a in alns if a.target in (site.ref_target, site.ctx_target)
                         and _within_reach(a, site)), default=0)
        item.discriminating = self_score - ref_score >= margin
        if not item.discriminating:
            item.reason = "no_discriminating_sequence"


def junction_carriers(ends, alignments):
    """ALT sequences of other classes that contain each class's junction.

    A sequence that aligns across an end's breakpoint, with that end's required
    flanks, at cluster identity carries the same adjacency (for example the
    other end of a junction with inserted sequence, or a neighbouring junction
    of the same complex haplotype on one read).  It explains a junction read
    consistently, so it is not a competing allele for that class.
    """
    by_name = {alt_target_name(e): e for e in ends if e.alt_seq and e.cls and not e.reason}
    carriers = defaultdict(set)
    for qname, q in by_name.items():
        ret = required_side(q.alt_bp)
        cont = q.req_cont or required_side(q.cont_len)
        for aln in alignments.get(qname, ()):
            t = by_name.get(aln.target)
            if t is None or t.cls == q.cls:
                continue
            if aln.alnlen == 0 or aln.matches / aln.alnlen < CLUSTER_IDENTITY:
                continue
            if aln.qs <= q.alt_bp - ret and aln.qe >= q.alt_bp + cont:
                carriers[q.cls].add(aln.target)
    return carriers


def _identity_beyond(aln, qpos):
    """(identical bases, alignment columns) of ``aln`` at query positions >= qpos."""
    same = cols = 0
    step = 1 if aln.strand == "+" else -1
    q = aln.qs if step == 1 else aln.qe
    for n, op in re.findall(r"(\d+)([MIDNX=])", aln.cigar):
        n = int(n)
        if op in "MX=I":
            lo, hi = (q, q + n) if step == 1 else (q - n, q)
            inside = max(0, hi - max(lo, qpos))
            cols += inside
            if op in "M=":  # without --eqx an M is counted as identical
                same += inside
            q += step * n
        elif op in "DN" and q >= qpos:  # a deletion at the breakpoint itself breaks the junction
            cols += n
    return same, cols


def junction_homology(item, site, alignments):
    """ALT bases past the breakpoint that the reference near the site also has.

    An ALT-to-reference alignment that runs through the end of the retained
    flank and on past the breakpoint means the reference already holds the
    junction sequence there: as the reference allele itself (the continuation
    stays collinear past the reference breakpoint) or as another copy within
    the reach of the reads fetched at this site (a tandem or nearby repeat;
    reads from that copy cross the junction without carrying it).  When that
    copy is near-identical, reads inside it cannot tell the alleles apart.
    A copy counts only if it also holds the retained flank that an ALT read
    must cover (less MIN_SIDE of alignment extension into diverged sequence);
    reads from a copy that shares less of it never qualify as ALT.
    """
    ret = min(site.req_ret, required_side(item.alt_bp))
    need = max(MIN_SIDE, ret - MIN_SIDE)
    best = 0
    for aln in alignments.get(alt_target_name(item), ()):
        if aln.target not in (site.ref_target, site.ctx_target) or not aln.cigar:
            continue
        if not _within_reach(aln, site):
            continue
        if not (aln.qs <= item.alt_bp - need and aln.qe >= item.alt_bp + MIN_SIDE):
            continue
        same, cols = _identity_beyond(aln, item.alt_bp)
        if cols and same / cols >= HOMOLOGY_IDENTITY:
            best = max(best, aln.qe - item.alt_bp)
    return best


def apply_junction_homology(ends, sites, alignments):
    """Move each site's discriminating point past ALT/REF homology.

    Every class at a site, and REF, must be covered beyond the longest
    homologous continuation there, so all alleles are counted from reads of the
    same reach.  An end whose continuation does not extend MIN_SIDE past that
    point has no discriminating sequence.
    """
    by_site = defaultdict(list)
    for item in ends:
        if not item.reason and item.site >= 0:
            item.homology = junction_homology(item, sites[item.site], alignments)
            by_site[item.site].append(item)
    for site_id, items in by_site.items():
        site = sites[site_id]
        while items:
            reach = max(e.homology for e in items)
            short = [e for e in items if e.cont_len - reach < MIN_SIDE]
            if not short:
                break
            for e in short:
                e.reason = "homologous_junction"
                e.discriminating = False
            items = [e for e in items if not e.reason]
        if not items:
            continue
        reach = max(e.homology for e in items)
        site.req_cont = reach + required_side(min(e.cont_len - reach for e in items))
        for e in items:
            e.req_cont = site.req_cont


def _covers(aln, bp, left_need, right_need):
    return aln.ts <= bp - left_need and aln.te >= bp + right_need


def _ref_needs(site):
    if site.retained == "left":
        return site.req_ret, site.req_cont
    return site.req_cont, site.req_ret


def _anchored_crossing(site, alns):
    """Read leaves the reference at this breakpoint into unexplained sequence."""
    off = site.ref_offset
    for aln in alns:
        if aln.target != site.ref_target:
            continue
        if site.retained == "left":
            anchored = aln.ts <= off - site.req_ret and off - OTHER_STOP_TOLERANCE <= aln.te < off + site.req_cont
            beyond = aln.qlen - aln.qe if aln.strand == "+" else aln.qs
        else:
            anchored = aln.te >= off + site.req_ret and off - site.req_cont < aln.ts <= off + OTHER_STOP_TOLERANCE
            beyond = aln.qs if aln.strand == "+" else aln.qlen - aln.qe
        if anchored and beyond >= site.req_cont:
            return True
    return False


def _query_position(aln, target_pos):
    """Read coordinate of a target position (through structural indels; small ones ignored)."""
    target_pos = min(max(target_pos, aln.ts), aln.te)
    offset = target_pos - aln.ts
    for t, _, op, n in aln.indels:
        if t <= target_pos:
            offset += n if op == "I" else -min(n, target_pos - t)
    offset = min(max(offset, 0), aln.qe - aln.qs)
    return aln.qs + offset if aln.strand == "+" else aln.qe - offset


def _window_score(aln, lo, hi):
    """Share of an alignment's score that falls in read interval [lo, hi)."""
    overlap = min(aln.qe, hi) - max(aln.qs, lo)
    if overlap <= 0:
        return 0.0
    return aln.score * overlap / max(1, aln.qe - aln.qs)


def call_read(site, alns, alt_info, target_class, margin, carriers=None):
    """Return REF, ('ALT', class), OTHER, AMB or NA for one read at one site.

    A class qualifies when one of its sequences from an end at this breakpoint
    (``site.local_targets``, same orientation as the site) covers the junction
    with the required flanks: the site's, and the sequence's own reach past its
    homology, whichever is longer.  Explanations are compared on the same stretch of the
    read, FLANK on each side of the breakpoint, so other junctions elsewhere on
    a long read and the differing lengths of the target sequences do not
    compete; the winner must beat every other target by ``margin``, except
    sequences that carry the winning junction themselves (``carriers``).
    A structural indel at the breakpoint cuts the read's alignment there
    (``cut_at``).
    """
    pieces = []
    for aln in alns:
        if aln.target == site.ref_target:
            pieces += cut_at(aln, site.ref_offset)
        elif aln.target in site.local_targets:
            pieces += cut_at(aln, alt_info[aln.target].alt_bp)
        else:
            pieces.append(aln)
    alns = pieces
    covering = []
    ref_class = site.ref_target
    left_need, right_need = _ref_needs(site)
    for aln in alns:
        cls = target_class.get(aln.target, aln.target)
        if aln.target == ref_class:
            if _covers(aln, site.ref_offset, left_need, right_need):
                covering.append((cls, aln, site.ref_offset))
        elif aln.target in site.local_targets:
            item = alt_info[aln.target]
            ret = min(site.req_ret, required_side(item.alt_bp))
            cont = max(site.req_cont, item.req_cont or required_side(item.cont_len))
            if _covers(aln, item.alt_bp, ret, cont):
                covering.append((cls, aln, item.alt_bp))
    if not covering:
        return ("OTHER", None) if _anchored_crossing(site, alns) else ("NA", None)
    _, top, top_bp = max(covering, key=lambda c: c[1].score)
    q_bp = _query_position(top, top_bp)
    lo, hi = q_bp - FLANK, q_bp + FLANK
    class_score = defaultdict(lambda: -math.inf)
    for cls, aln, _ in covering:
        class_score[cls] = max(class_score[cls], _window_score(aln, lo, hi))
    best = max(sorted(class_score), key=class_score.get)
    allies = (carriers or {}).get(best, ())
    runner = max((_window_score(aln, lo, hi) for aln in alns
                  if target_class.get(aln.target, aln.target) != best and aln.target not in allies),
                 default=-math.inf)
    if class_score[best] - runner < margin:
        return "AMB", best
    return ("REF", None) if best == ref_class else ("ALT", best)


def load_depth(main_stat_path):
    rows = defaultdict(list)
    depths = []
    opener = gzip.open if str(main_stat_path).endswith(".gz") else open
    with opener(main_stat_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] == "chrM":
                continue
            st, nd, depth = int(f[1]), int(f[2]), float(f[7])
            rows[f[0]].append((st, nd, depth))
            depths.append(depth)
    for chrom in rows:
        rows[chrom].sort()
    return rows, float(np.median(depths)) / 2


def local_cn(depth_rows, n_unit, chrom, bp, retained, bound=None):
    """Median depth (N) of the nearest whole windows inside the retained segment.

    ``bound`` is the nearest other breakpoint on the retained side; windows
    beyond it may hold another copy number.  NaN when no whole window fits.
    """
    windows = depth_rows.get(chrom, ())
    if retained == "left":
        lo = bp - CN_SEARCH if bound is None else max(bp - CN_SEARCH, bound)
        side = [w for w in windows if w[1] < bp and w[0] > lo]
        side = side[-CN_WINDOWS:]
    else:
        hi = bp + CN_SEARCH if bound is None else min(bp + CN_SEARCH, bound)
        side = [w for w in windows if w[0] > bp and w[1] < hi]
        side = side[:CN_WINDOWS]
    if not side:
        return float("nan")
    return float(np.median([w[2] for w in side])) / n_unit


def segment_bounds(sites):
    """Nearest other breakpoint on each site's retained side (None if none).

    Breakpoints closer than MIN_SIDE are the other side of the same event
    (target-site duplication, reciprocal pair), not the end of a segment.
    """
    by_chrom = defaultdict(list)
    for site in sites:
        by_chrom[site.chrom].append(site.bp)
    bounds = {}
    for site in sites:
        others = by_chrom[site.chrom]
        if site.retained == "left":
            before = [x for x in others if x <= site.bp - MIN_SIDE]
            bounds[site.site_id] = max(before) if before else None
        else:
            after = [x for x in others if x >= site.bp + MIN_SIDE]
            bounds[site.site_id] = min(after) if after else None
    return bounds


def assign_local_cn(sites, depth_rows, n_unit):
    """CN_local per site: segment windows, else calibrated depth at the breakpoint."""
    bounds = segment_bounds(sites)
    for site in sites:
        site.cn_local = local_cn(depth_rows, n_unit, site.chrom, site.bp, site.retained,
                                 bounds[site.site_id])
        site.cn_source = "windows" if site.cn_local == site.cn_local else ""
    rates = [s.mapq_records / s.cn_local for s in sites
             if s.cn_source == "windows" and bounds[s.site_id] is None and s.cn_local >= 0.5
             and s.unique_fraction >= MIN_UNIQUE_FRACTION and s.mapq_records > 0]
    per_copy = float(np.median(rates)) if len(rates) >= 5 else n_unit
    for site in sites:
        if not site.cn_source and site.mapq_records > 0:
            site.cn_local = site.mapq_records / per_copy
            site.cn_source = "point"
    return per_copy


def _site_testable(site, ends_here):
    if not any(e.anchor_mapq >= MIN_ANCHOR_MAPQ for e in ends_here):
        return "anchor_not_unique"
    if not (site.unique_fraction >= MIN_UNIQUE_FRACTION):
        return "reads_not_unique"
    if not (site.cn_local >= MIN_LOCAL_CN):
        return "no_local_copy"
    informative = site.counts.get("informative", 0)
    if informative < MIN_SITE_READS:
        return "too_few_crossing_reads"
    amb = site.counts.get("AMB", 0)
    if amb / (amb + informative) > MAX_AMBIGUOUS_FRACTION:
        return "ambiguous_reads"
    return ""


def measure(contig_data, members_by_nclose, *, unitig_fasta, reference_fasta, bam_path,
            main_stat_path, paf_paths, threads, workdir, preset=None, extra_sequences=None):
    """Measure every member end and return (ends, sites, n_unit, read_calls)."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    source = SequenceSource(unitig_fasta, reference_fasta, paf_paths, extra_sequences)
    ends = build_member_ends(contig_data, members_by_nclose, source)
    sites = build_sites(ends)
    ref_targets = build_ref_targets(sites, source.reference)
    preset = preset or detect_preset(bam_path)
    match_score = MATCH_SCORE.get(preset, 1)
    margin = SCORE_MARGIN_BASES * match_score
    logging.info("Read ratio: %d member end(s), %d site(s), preset %s", len(ends), len(sites), preset)

    sequences = fetch_site_reads(bam_path, sites)
    context_targets = build_context_targets(sites, source.reference)
    scratch = Path(tempfile.mkdtemp(prefix="reads_", dir=workdir))
    try:
        targets_fa = scratch / "targets.fa"
        queries_fa = scratch / "queries.fa"
        alt_info = {}
        with open(targets_fa, "w") as out:
            for name, seq in ref_targets.items():
                out.write(f">{name}\n{seq}\n")
            for item in ends:
                if item.alt_seq and not item.reason:
                    alt_info[alt_target_name(item)] = item
                    out.write(f">{alt_target_name(item)}\n{item.alt_seq}\n")
        with open(queries_fa, "w") as out:
            for name in alt_info:
                out.write(f">{name}\n{alt_info[name].alt_seq}\n")
            for name, seq in sequences.items():
                out.write(f">{name}\n{seq}\n")
        paf = scratch / "alignments.paf"
        run_minimap2(str(targets_fa), str(queries_fa), str(paf), preset, threads)
        alignments = parse_paf(paf)
        # Reference copies of the ALT junctions anywhere within read reach of their sites.
        context_fa = scratch / "context.fa"
        alt_fa = scratch / "alt.fa"
        with open(context_fa, "w") as out:
            for name, seq in context_targets.items():
                out.write(f">{name}\n{seq}\n")
        with open(alt_fa, "w") as out:
            for name in alt_info:
                out.write(f">{name}\n{alt_info[name].alt_seq}\n")
        context_paf = scratch / "context.paf"
        if alt_info:
            run_minimap2(str(context_fa), str(alt_fa), str(context_paf), preset, threads)
            for name, alns in parse_paf(context_paf).items():
                alignments[name].extend(alns)
        for src_path, dst_name in ((paf, "alignments.paf.gz"), (context_paf, "context_alignments.paf.gz")):
            if not src_path.exists():
                continue
            with open(src_path) as src, gzip.open(workdir / dst_name, "wt") as dst:
                for line in src:
                    fields = line.rstrip("\n").split("\t")
                    score = [x for x in fields[12:] if x.startswith("AS:i:")]
                    dst.write("\t".join(fields[:12] + score) + "\n")
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    cluster_alt_targets(ends, alignments)
    mark_discriminating(ends, sites, alignments, match_score)
    apply_junction_homology(ends, sites, alignments)
    carriers = junction_carriers(ends, alignments)
    target_class = {name: item.cls for name, item in alt_info.items()}
    usable = [e for e in ends if e.cls and e.site >= 0 and not e.reason]
    for site in sites:
        for item in usable:
            other = sites[item.site]
            if (other.chrom == site.chrom and other.retained == site.retained
                    and abs(other.bp - site.bp) <= SITE_TOLERANCE):
                site.local_targets.add(alt_target_name(item))
    for item in usable:
        if item.alt_bp < FLANK:
            item.reason = "short_retained_segment"

    depth_rows, n_unit = load_depth(main_stat_path)
    assign_local_cn(sites, depth_rows, n_unit)
    read_calls = []
    for site in sites:
        counts = defaultdict(int)
        alt_names = defaultdict(set)
        for name in site.reads:
            call, cls = call_read(site, alignments.get(name, ()), alt_info, target_class, margin, carriers)
            counts[call] += 1
            if call == "ALT":
                alt_names[cls].add(name)
            read_calls.append((site.site_id, name, call, cls or "."))
        counts["informative"] = counts["REF"] + counts["ALT"] + counts["OTHER"]
        site.counts = dict(counts)
        site.alt_names = dict(alt_names)
        site.alt_counts = {cls: len(names) for cls, names in alt_names.items()}
    return ends, sites, n_unit, read_calls, preset


def _component_estimate(component, members, class_est):
    """Most junction copies the reads allow for members linked by shared classes.

    Each member is a junction molecule, and every molecule carrying an allele
    class contributes reads to it, so the copies of the members carrying a
    tested class add up to at most that class's estimate (once per member,
    whichever ends carry it); an untested class bounds nothing.  The bound is
    the largest total these constraints allow (a small linear program).
    Returns (copies, reads per N): each member is limited by a tight class at
    the optimum, the best-sampled one counts, and the least-sampled member sets
    the component's reads per N.  Copies are math.inf when a member carries no
    tested class.
    """
    carried = [sorted({e.cls for e in members[m] if e.cls in class_est}) for m in component]
    if any(not classes for classes in carried):
        return math.inf, float("nan")
    classes = sorted({c for cs in carried for c in cs})
    a_ub = np.array([[1.0 if c in cs else 0.0 for cs in carried] for c in classes])
    b_ub = np.array([class_est[c]["j"] for c in classes])
    res = linprog(-np.ones(len(component)), A_ub=a_ub, b_ub=b_ub, bounds=(0, None), method="highs")
    if res.status != 0:
        raise RuntimeError(f"member copy bound failed: {res.message}")
    slack = dict(zip(classes, b_ub - a_ub @ res.x))
    reads = math.inf
    for cs in carried:
        tight = [c for c in cs if slack[c] <= 1e-9 + 1e-7 * class_est[c]["j"]] or cs
        reads = min(reads, max(class_est[c]["r_eff"] for c in tight))
    return float(-res.fun), reads


def summarize(ends, sites):
    """Per-NClose junction copy estimates from the member-end measurements.

    An allele class is estimated only at this NClose's own testable sites, so a
    class shared with another NClose never lends it reads from elsewhere.  A
    junction molecule has both ends, so a member is bounded by each of its
    tested classes; members linked by shared classes split those bounds
    (``_component_estimate``) and unlinked groups add up.

    For the post-fit test a read is one molecule: a junction read is fetched at
    both of its breakpoints, so the ALT reads are the distinct reads called for
    the NClose's tested classes at its tested sites.  The exposure (reads per N)
    is the least-sampled tested end of any member.  Both ends of a junction see
    the same molecules, so ends whose reads per copy differ disagree about the
    calibration (CN_local or read collection at one of them), not about the
    junction, and the end predicting fewer reads is the one that cannot
    overstate the weight's expectation.  A member without a tested end of its
    own leaves the exposure undefined.
    """
    ends_by_site = defaultdict(list)
    for item in ends:
        if item.site >= 0:
            ends_by_site[item.site].append(item)
    site_reason = {s.site_id: _site_testable(s, ends_by_site[s.site_id]) for s in sites}
    owners = defaultdict(set)
    for item in ends:
        if item.cls:
            owners[item.cls].add(item.nclose)

    by_nclose = defaultdict(list)
    for item in ends:
        by_nclose[item.nclose].append(item)
    def reads_per_n(sid):
        return sites[sid].counts["informative"] / sites[sid].cn_local

    results = {}
    for nclose, items in by_nclose.items():
        class_sites = defaultdict(set)
        member_sites = defaultdict(set)
        for item in items:
            # The end's own anchor must be uniquely placed: another NClose's anchor at the same
            # breakpoint does not place this one.
            if (item.cls and item.site >= 0 and not item.reason and not site_reason[item.site]
                    and item.anchor_mapq >= MIN_ANCHOR_MAPQ):
                class_sites[item.cls].add(item.site)
                member_sites[item.member].add(item.site)
        class_est = {}
        for cls, site_ids in class_sites.items():
            alt = sum(sites[sid].alt_counts.get(cls, 0) for sid in site_ids)
            rates = [reads_per_n(sid) for sid in site_ids]
            class_est[cls] = dict(j=alt / sum(rates), r_eff=max(rates))

        members = defaultdict(list)
        for item in items:
            members[item.member].append(item)
        member_classes, discordant = {}, False
        for member, member_ends in members.items():
            classes = {e.cls for e in member_ends if e.cls}
            member_classes[member] = classes
            tested = [class_est[c]["j"] for c in classes if c in class_est]
            if tested:
                discordant |= min(tested) <= DISPLAY_N and max(tested) >= 2.5 * DISPLAY_N
        groups = []
        for member, classes in member_classes.items():
            merged = [g for g in groups if g["classes"] & classes]
            group = {"members": [member], "classes": set(classes)}
            for g in merged:
                group["members"] = g["members"] + group["members"]
                group["classes"] |= g["classes"]
                groups.remove(g)
            groups.append(group)
        j_total = 0.0
        r_eff = math.inf
        untested = 0
        for group in groups:
            copies, reads = _component_estimate(group["members"], members, class_est)
            if copies == math.inf:
                untested += 1
                continue
            j_total += copies
            r_eff = min(r_eff, reads)
        alt_names = set()
        for sid in set().union(*class_sites.values()):
            for cls in class_est:
                alt_names |= sites[sid].alt_names.get(cls, set())
        exposure = (min(reads_per_n(sid) for m in members for sid in member_sites[m])
                    if all(member_sites[m] for m in members) else float("nan"))
        shared = sorted({c for classes in member_classes.values() for c in classes
                         if len(owners[c]) > 1})
        if untested == len(groups):
            status = "NOT_TESTABLE"
            reasons = sorted({e.reason or site_reason.get(e.site, "")
                              or ("anchor_not_unique" if e.anchor_mapq < MIN_ANCHOR_MAPQ else "")
                              for e in items} - {""})
            reason = ";".join(reasons) or "no_testable_end"
            j_total, r_eff = float("nan"), float("nan")
        elif untested:
            status, reason = "PARTIAL", f"{untested}_untested_member_group(s)"
        else:
            status, reason = "TESTABLE", ""
        results[nclose] = dict(
            status=status, reason=reason, j=j_total, r_eff=r_eff,
            alt_reads=len(alt_names), exposure=exposure if status == "TESTABLE" else float("nan"),
            members=len(members), groups=len(groups), shared=";".join(shared),
            classes=";".join(sorted({e.cls for e in items if e.cls})),
            end_discordant=int(discordant),
        )
    return results, site_reason


def prefit_decision(est):
    """(A): remove a powered NClose whose read-implied copy number is <= 0.2N."""
    if est["status"] != "TESTABLE":
        return ""
    if DISPLAY_N * est["r_eff"] < MIN_EXPECTED_ALT_AT_DISPLAY:
        return ""
    return "read_ratio_low" if est["j"] <= DISPLAY_N else ""


def poisson_cdf(k, mu):
    if mu <= 0:
        return 1.0
    k = int(math.floor(k + 1e-9))
    if k < 0:
        return 0.0
    # Regularized incomplete gamma: no underflow of exp(-mu) at large means.
    return min(1.0, float(pdtr(k, mu)))


def marginal_poisson_cdf(k, mu, scale_sd):
    """P(X <= k) for X ~ Poisson(mu * exp(scale_sd * Z)), Z standard normal.

    Converting a copy number into an expected read count uses the reads per
    copy at the site (informative reads / CN_local).  That rate is itself
    uncertain (CN_local from depth windows, local coverage, required reach), so
    the read-sampling tail is averaged over the sample's scatter of that rate.
    """
    if scale_sd <= 0 or mu <= 0:
        return poisson_cdf(k, mu)

    def integrand(z):
        return poisson_cdf(k, mu * math.exp(scale_sd * z)) * math.exp(-0.5 * z * z)

    # The tail steps from 1 to 0 where mu * exp(scale_sd * z) passes k; a fixed
    # quadrature misses a sharp step, so integrate adaptively around it.
    bound = 12.0
    step = min(max(math.log(max(k, 0.5) / mu) / scale_sd, -bound + 1e-6), bound - 1e-6)
    total, _ = quad(integrand, -bound, bound, points=[step], limit=200, epsabs=1e-13, epsrel=1e-9)
    return min(1.0, total / math.sqrt(2 * math.pi))


def scale_scatter(sites, site_reason):
    """Robust SD of log(reads per copy) over the testable sites of a sample."""
    logs = [math.log(s.counts["informative"] / s.cn_local) for s in sites
            if not site_reason.get(s.site_id) and s.cn_local > 0 and s.counts.get("informative", 0) > 0]
    if len(logs) < 5:
        return MAX_SCALE_SD
    centre = float(np.median(logs))
    mad = float(np.median([abs(x - centre) for x in logs])) * 1.4826
    return min(MAX_SCALE_SD, max(MIN_SCALE_SD, mad))


def posthoc_decision(est, weight, scale_sd=0.0):
    """(B): remove an NClose whose fitted weight exceeds what its reads allow.

    The observed ALT reads are tested against the fitted weight placed where it
    would predict the fewest of them (the least-sampled member), so any split
    of the weight among the members is at least as unlikely.
    """
    if est["status"] != "TESTABLE" or not weight > 0 or not est["exposure"] > 0:
        return "", float("nan")
    p_value = marginal_poisson_cdf(est["alt_reads"], weight * est["exposure"], scale_sd)
    if weight - est["j"] >= DISPLAY_N and p_value < POSTHOC_TAIL:
        return "read_ratio_overweight", p_value
    return "", p_value


END_COLUMNS = ["nclose", "member", "end", "chrom", "bp", "retained", "anchor_mapq", "seq_source",
               "retained_len", "continuation_len", "homology", "class", "site", "end_reason"]
SITE_COLUMNS = ["site", "site_chrom", "site_bp", "site_retained", "ref_target", "ctx_target", "req_ret", "req_cont",
                "reads", "sampled_fraction", "unique_fraction", "cn_local", "cn_source", "REF", "ALT_class",
                "ALT_other_class", "OTHER", "AMB", "NA", "informative", "alt_ratio", "j_site",
                "r_site", "site_reason"]


def write_tables(outdir, ends, sites, site_reason, read_calls):
    outdir = Path(outdir)
    with open(outdir / "member_ends.tsv", "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(END_COLUMNS + SITE_COLUMNS)
        for item in ends:
            row = [identity_key(item.nclose), identity_key(item.member), item.end, item.chrom,
                   item.bp, item.retained, item.anchor_mapq, item.seq_source, item.alt_bp,
                   item.cont_len, item.homology, item.cls or ".", item.site, item.reason or "."]
            if item.site >= 0:
                s = sites[item.site]
                inf = s.counts.get("informative", 0)
                alt_cls = s.alt_counts.get(item.cls, 0)
                alt_other = sum(v for c, v in s.alt_counts.items() if c != item.cls)
                ratio = alt_cls / inf if inf else float("nan")
                row += [s.site_id, s.chrom, s.bp, s.retained, s.ref_target, s.ctx_target, s.req_ret, s.req_cont,
                        len(s.reads), f"{s.sampled_fraction:.3f}", f"{s.unique_fraction:.3f}",
                        f"{s.cn_local:.3f}", s.cn_source or ".", s.counts.get("REF", 0), alt_cls, alt_other,
                        s.counts.get("OTHER", 0), s.counts.get("AMB", 0), s.counts.get("NA", 0), inf,
                        f"{ratio:.4f}", f"{ratio * s.cn_local:.4f}",
                        f"{(inf / s.cn_local) if s.cn_local > 0 else float('nan'):.3f}",
                        site_reason.get(s.site_id) or "."]
            else:
                row += ["."] * len(SITE_COLUMNS)
            writer.writerow(row)
    with gzip.open(outdir / "read_calls.tsv.gz", "wt") as handle:
        handle.write("site\tread\tcall\tclass\n")
        for row in read_calls:
            handle.write("\t".join(map(str, row)) + "\n")


NCLOSE_COLUMNS = ["identity", "contig", "start_chr", "start_pos", "end_chr", "end_pos", "members",
                  "groups", "classes", "shared_classes", "status", "status_reason", "j_hat_N",
                  "reads_per_N", "power_ok", "end_discordant", "prefit_decision", "junctions"]


def member_junctions(contig_data, members):
    """Breakpoints of every member junction, 'chrom:bp:side>chrom:bp:side' joined by ';'.

    The side is the retained reference side (L/R), as stage 24 writes endpoints.
    """
    junctions = []
    for _, pair in members:
        parts = []
        for end, idx in (("A", pair[0]), ("B", pair[1])):
            row = contig_data[int(idx)]
            bp, retained = anchor_breakpoint(row, end)
            parts.append(f"{row[CHR_NAM]}:{bp}:{'L' if retained == 'left' else 'R'}")
        junctions.append(">".join(parts))
    return ";".join(junctions)


def nclose_rows(contig_data, results, members_by_nclose):
    rows = []
    for nclose, est in results.items():
        contig, pair = nclose
        a, b = contig_data[pair[0]], contig_data[pair[1]]
        power = (est["status"] == "TESTABLE"
                 and DISPLAY_N * est["r_eff"] >= MIN_EXPECTED_ALT_AT_DISPLAY)
        rows.append(dict(
            identity=identity_key(nclose), contig=contig, start_chr=a[CHR_NAM],
            start_pos=a[CHR_STR], end_chr=b[CHR_NAM], end_pos=b[CHR_STR],
            members=est["members"], groups=est["groups"], classes=est["classes"] or ".",
            shared_classes=est["shared"] or ".", status=est["status"],
            status_reason=est["reason"] or ".", j_hat_N=f"{est['j']:.4f}",
            reads_per_N=f"{est['r_eff']:.3f}", power_ok=int(power),
            end_discordant=est["end_discordant"],
            prefit_decision=prefit_decision(est) or "KEEP",
            junctions=member_junctions(contig_data, members_by_nclose.get(nclose, [nclose])),
        ))
    return rows


def write_rows(path, columns, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def save_members(prefix, members_by_nclose):
    outdir = Path(prefix) / OUTPUT_DIR
    outdir.mkdir(parents=True, exist_ok=True)
    payload = {identity_key(n): [identity_key(m) for m in members]
               for n, members in members_by_nclose.items()}
    with open(outdir / MEMBERS_JSON, "w") as handle:
        json.dump(payload, handle, indent=1)
        handle.write("\n")


def load_members(prefix):
    path = Path(prefix) / OUTPUT_DIR / MEMBERS_JSON
    with open(path) as handle:
        payload = json.load(handle)
    return {parse_identity_key(k): [parse_identity_key(m) for m in v] for k, v in payload.items()}


def read_remove_list(path):
    """Stage-01 NCloses of a post-fit removal list (stage-24 rows go to stage 24)."""
    if not path:
        return {}
    removed = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("origin", "stage01") == "stage01":
                removed[parse_identity_key(row["identity"])] = row
    return removed


def evaluate_stage01(prefix, contig_data, members_by_nclose, *, unitig_fasta, reference_fasta,
                     bam_path, main_stat_path, paf_paths, threads):
    """Stage-01 entry: measure, write tables, and return {identity: (decision, est)}."""
    outdir = Path(prefix) / OUTPUT_DIR
    outdir.mkdir(parents=True, exist_ok=True)
    # The post-fit round reruns stage 01 on identical inputs; reuse that measurement.
    def stamp(path):
        return None if not path else (str(path), os.path.getsize(path), os.path.getmtime(path))

    # Whole rows: the global index (PAF row) names the unitig a renamed child is read from.
    key = hashlib.sha256(json.dumps([
        sorted((identity_key(n), [identity_key(m) for m in ms]) for n, ms in members_by_nclose.items()),
        [list(map(str, contig_data[idx])) for ms in members_by_nclose.values()
         for _, pair in ms for idx in range(min(map(int, pair)), max(map(int, pair)) + 1)],
        [stamp(x) for x in (unitig_fasta, reference_fasta, bam_path, main_stat_path)],
        [stamp(x) for x in paf_paths or ()],
        hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    ]).encode()).hexdigest()
    cache = outdir / "stage01" / "measure_cache.pkl"
    cached = None
    if cache.exists():
        with open(cache, "rb") as handle:
            cached = pickle.load(handle)
    if cached is not None and cached.get("key") == key:
        logging.info("NClose read ratio: reusing the stage-01 measurement")
        ends, sites, n_unit, read_calls, preset = cached["value"]
    else:
        measured = measure(
            contig_data, members_by_nclose, unitig_fasta=unitig_fasta,
            reference_fasta=reference_fasta, bam_path=bam_path, main_stat_path=main_stat_path,
            paf_paths=paf_paths, threads=threads, workdir=outdir / "stage01")
        with open(cache, "wb") as handle:
            pickle.dump(dict(key=key, value=measured), handle)
        ends, sites, n_unit, read_calls, preset = measured
    results, site_reason = summarize(ends, sites)
    write_tables(outdir / "stage01", ends, sites, site_reason, read_calls)
    rows = nclose_rows(contig_data, results, members_by_nclose)
    write_rows(outdir / NCLOSE_TSV, NCLOSE_COLUMNS, rows)
    write_rows(Path(prefix) / NCLOSE_TSV, NCLOSE_COLUMNS, rows)
    # Stage 24 does not rescue these junctions from reads again.
    write_rows(outdir / PREFIT_REMOVE_TSV, NCLOSE_COLUMNS, [r for r in rows if r["prefit_decision"] != "KEEP"])
    summary = dict(preset=preset, n_unit=n_unit, member_ends=len(ends), sites=len(sites),
                   unitig_fasta=os.path.abspath(unitig_fasta) if unitig_fasta else None,
                   ncloses=len(results),
                   status={s: sum(r["status"] == s for r in results.values())
                           for s in ("TESTABLE", "PARTIAL", "NOT_TESTABLE")},
                   prefit_removed=sum(bool(prefit_decision(r)) for r in results.values()))
    with open(outdir / "stage01_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    return {n: (prefit_decision(est), est) for n, est in results.items()}, summary


def _node_key(text):
    try:
        return tuple(sorted(int(x) for x in text.split(":")))
    except (AttributeError, ValueError):
        return None


def report_unit_ratio(prefix, main_stat_path):
    """(report 1N / (median depth / 2), report 1N, its source) of stage 31's NClose CN.

    Stage 31 writes ``nclose_cn`` in the copy unit of depth_scale.json: the
    B-allele-fraction 1N when stage 00 estimated one, otherwise median depth / 2.
    This module measures local copy numbers, reads per copy and its 0.2N cut in
    median depth / 2, the threshold unit of depth_scale.  depth_scale.copy_unit
    on the same non-chrM depths returns the unit stage 31 used, so multiplying a
    reported CN by the ratio puts it in this module's unit.  The ratio is
    exactly 1 when the report used median depth / 2.
    """
    rows, half = load_depth(main_stat_path)
    depth = [window[2] for windows in rows.values() for window in windows]
    unit, source = depth_scale.copy_unit(str(prefix), depth)
    return unit / half, unit, source


def _report_weights(prefix, scale=1.0):
    """Fitted NClose weight per sorted node key, the max over fit stages.

    Weights are attributed to NCloses by stage 31 (nclose_report*.tsv).  Stage
    31 reports NCloses of identical BND geometry under one canonical key and
    lists the others as its sources (nclose_sources*.tsv); each source NClose
    carries that weight.  A report older than the NClose set it would describe
    is refused.  ``scale`` multiplies every reported CN (report_unit_ratio
    converts stage 31's copy unit into median depth / 2).
    """
    weights = defaultdict(float)
    per_stage = defaultdict(dict)
    prefix = Path(prefix)
    if not (prefix / "nclose_report.tsv").exists():
        raise FileNotFoundError(f"{prefix / 'nclose_report.tsv'}: run stage 31 before the post-fit check")
    built = (prefix / "01_nclose_data.pkl").stat().st_mtime
    for stage, suffix in (("raw", ""), ("filter", "_filter"), ("cluster", "_cluster")):
        path = prefix / f"nclose_report{suffix}.tsv"
        if not path.exists():
            continue
        sources = prefix / f"nclose_sources{suffix}.tsv"
        for item in (path, sources):
            if item.exists() and item.stat().st_mtime < built:
                raise RuntimeError(f"{item} predates 01_nclose_data.pkl; rerun stage 31")
        stage_weight = {}
        with open(path) as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row.get("kind", "bnd") != "bnd":
                    continue
                key = _node_key(row.get("nclose_key", ""))
                if key is None or len(key) != 2:
                    continue
                try:
                    value = max(0.0, float(row["nclose_cn"])) * scale
                except ValueError:
                    value = 0.0
                stage_weight[key] = max(stage_weight.get(key, 0.0), value)
        if sources.exists():
            with open(sources) as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    source = _node_key(row.get("source_nclose_key", ""))
                    canonical = _node_key(row.get("canonical_nclose_key", ""))
                    if source and canonical in stage_weight and source not in stage_weight:
                        stage_weight[source] = stage_weight[canonical]
        for key, value in stage_weight.items():
            per_stage[key][stage] = value
            weights[key] = max(weights[key], value)
    return weights, per_stage


POSTHOC_COLUMNS = ["identity", "contig", "start_chr", "start_pos", "end_chr", "end_pos", "origin",
                   "members", "status", "status_reason", "end_discordant", "j_hat_N", "reads_per_N", "w_raw",
                   "w_filter", "w_cluster", "w_fit", "expected_alt", "observed_alt", "p_lower",
                   "posthoc_decision", "removable_at_stage01", "junctions"]


def posthoc(prefix, *, unitig_fasta, reference_fasta, bam_path, main_stat_path, threads, tag=""):
    """(B) after the fits: measure the final NCloses and list over-weighted ones.

    A non-empty ``tag`` writes a report-only copy (no removal list) under that suffix.
    """
    suffix = f"_{tag}" if tag else ""
    from breakend_graph import load_stage10_input

    prefix = Path(prefix)
    outdir = prefix / OUTPUT_DIR
    stage = load_stage10_input(prefix / "01_nclose_data.pkl")
    contig_data = stage.contig_data
    stage01_members = load_members(prefix)
    members_by_nclose = {}
    origin = {}
    for owner, pairs in stage.nclose_nodes.items():
        for pair in pairs:
            identity = (str(owner), tuple(int(x) for x in pair))
            members_by_nclose[identity] = stage01_members.get(identity, [identity])
            origin[identity] = "stage01" if identity in stage01_members else "stage24"
    with open(prefix / "paf_file_path.pkl", "rb") as handle:
        paf_paths = pickle.load(handle)
    ends, sites, n_unit, read_calls, preset = measure(
        contig_data, members_by_nclose, unitig_fasta=unitig_fasta,
        reference_fasta=reference_fasta, bam_path=bam_path, main_stat_path=main_stat_path,
        paf_paths=paf_paths, threads=threads, workdir=outdir / f"posthoc{suffix}",
        extra_sequences=load_rescue_templates(prefix))
    results, site_reason = summarize(ends, sites)
    write_tables(outdir / f"posthoc{suffix}", ends, sites, site_reason, read_calls)
    # Fitted weights in the unit of j, the exposure and the 0.2N cut (median depth / 2).
    to_threshold, report_unit, report_unit_source = report_unit_ratio(prefix, main_stat_path)
    weights, per_stage = _report_weights(prefix, to_threshold)
    scale_sd = scale_scatter(sites, site_reason)
    rows, removals = [], []
    for nclose, est in results.items():
        contig, pair = nclose
        key = tuple(sorted(pair))
        weight = weights.get(key, 0.0)
        decision, p_value = posthoc_decision(est, weight, scale_sd)
        a, b = contig_data[pair[0]], contig_data[pair[1]]
        row = dict(
            identity=identity_key(nclose), contig=contig, start_chr=a[CHR_NAM],
            start_pos=a[CHR_STR], end_chr=b[CHR_NAM], end_pos=b[CHR_STR], origin=origin[nclose],
            members=est["members"], status=est["status"], status_reason=est["reason"] or ".",
            end_discordant=est["end_discordant"],
            j_hat_N=f"{est['j']:.4f}", reads_per_N=f"{est['r_eff']:.3f}",
            w_raw=per_stage[key].get("raw", 0.0), w_filter=per_stage[key].get("filter", 0.0),
            w_cluster=per_stage[key].get("cluster", 0.0), w_fit=f"{weight:.4f}",
            expected_alt=f"{weight * est['exposure']:.2f}" if est["exposure"] > 0 else ".",
            observed_alt=est["alt_reads"] if est["exposure"] > 0 else ".",
            p_lower=f"{p_value:.3g}", posthoc_decision=decision or "KEEP",
            removable_at_stage01=int(origin[nclose] == "stage01"),
            junctions=member_junctions(contig_data, members_by_nclose[nclose]),
        )
        rows.append(row)
        if decision:
            removals.append(row)
    table = POSTHOC_TSV.replace(".tsv", f"{suffix}.tsv")
    write_rows(prefix / table, POSTHOC_COLUMNS, rows)
    write_rows(outdir / table, POSTHOC_COLUMNS, rows)
    summary = dict(preset=preset, n_unit=n_unit, report_unit=report_unit,
                   report_unit_source=report_unit_source, report_to_n_unit=to_threshold,
                   scale_sd=scale_sd, ncloses=len(results), flagged=len(removals),
                   removable_at_stage01=sum(r["removable_at_stage01"] for r in removals))
    if not tag:
        write_rows(outdir / REMOVE_TSV, POSTHOC_COLUMNS, removals)
        summary["remove_list"] = str(outdir / REMOVE_TSV)
    with open(outdir / f"posthoc{suffix}_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")
    return summary


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    sub = parser.add_subparsers(dest="command", required=True)
    post = sub.add_parser("posthoc", help="(B) compare fitted NClose weights with read ratios")
    post.add_argument("prefix")
    post.add_argument("--unitig-fasta", required=True)
    post.add_argument("--reference", required=True)
    post.add_argument("--bam", required=True)
    post.add_argument("--main-stat", required=True)
    post.add_argument("-t", "--threads", type=int, default=1)
    post.add_argument("--tag", default="", help="Report-only output suffix")
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.command == "posthoc":
        summary = posthoc(args.prefix, unitig_fasta=args.unitig_fasta, reference_fasta=args.reference,
                          bam_path=args.bam, main_stat_path=args.main_stat, threads=args.threads,
                          tag=args.tag)
        logging.info("Read-ratio posthoc: %s", json.dumps(summary))
    return 0


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s:%(message)s")
    raise SystemExit(main())
