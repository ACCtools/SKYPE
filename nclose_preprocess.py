"""Stage-01 NClose/telomere preprocessing shared by the native SKYPE CLI.

Assembly NCloses deliberately follow a small pipeline: legacy telomere/source
preprocessing, repeat trimming, one terminal pair per type-1/2 unitig, and one
global coordinate/direction clustering pass.  The graph-facing handoff remains
limited to ``contig_data``, ``nclose_nodes``, and ``telo_contig``.
"""

import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from skype_utils import *
from parse_vcf import parse_vcf_events, select_vcf_type4_graph_events
from nclose_tracking import (
    build_bnd_event_catalog,
    save_event_catalog,
)
from nclose_candidate import (
    NCloseCandidate,
    NCloseRejection,
    apply_nclose_filter,
    candidates_to_legacy,
)
from breakend_graph import save_stage10_input

import subprocess
import csv
import json
from dataclasses import dataclass, field, replace
from enum import Enum
from types import SimpleNamespace

import pickle as pkl
import pandas as pd
import numpy as np

import ast

import copy
import bisect

import graph_tool.all as gt
import networkx as nx

from collections import defaultdict, Counter

import logging

# logging 설정(레벨/포맷)은 skype_utils 에서 중앙 관리한다 (LOG_LEVEL).

CTG_NAM = 0
CTG_LEN = 1
CTG_STR = 2
CTG_END = 3
CTG_DIR = 4
CHR_NAM = 5
CHR_LEN = 6
CHR_STR = 7
CHR_END = 8
CTG_MAPQ = 9
CTG_TYP = 10
CTG_STRND = 11
CTG_ENDND = 12
CTG_TELCHR = 13
CTG_TELDIR = 14
CTG_TELCON = 15
CTG_RPTCHR = 16
CTG_RPTCASE = 17
CTG_CENSAT = 18
CTG_MAINFLOWDIR = 19
CTG_MAINFLOWCHR = 20
CTG_GLOBALIDX = 21

DIR_FOR = 1
DIR_BAK = 0
DIR_IN = 3
DIR_OUT = 2
INF = 1000000000
BUFFER = 10000000

CHUKJI_LIMIT = -1
BND_CONTIG_BOUND = 0.1
TYPE2_CONTIG_MINIMUM_LENGTH = 200*K
SUBTELOMERE_REPEAT_LENGTH = 0 # 100*K
TYPE2_FLANKING_LENGTH = 500 * K
TYPE2_SIM_COMPARE_RAITO = 1.3
TYPE34_BREAK_CHUKJI_LIMIT = 1*M
CHUKJI_FAIL_TYPE2_RESCUE_THRESHOLD = 2*K
CIRCUIT_ECDNA_LENGTH_LIMIT = 80*M

NCLOSE_SIM_DIFF_THRESHOLD = 5

NCLOSE_COMPRESS_LIMIT = 100*K
ALL_REPEAT_NCLOSE_COMPRESS_LIMIT = 500*K
SUBTELO_TIP_LIMIT = 500*K
OFFSET_DIR_GROUP_LIMIT = 100*K

NON_REPEAT_NOISE_RATIO=0.1

CONTIG_MINIMUM_SIZE = 100*K
BND_CONTIG_BOUND = 0.1
RPT_BND_CONTIG_BOUND = 0.2

MAPQ_BOUND = 60

TELOMERE_EXPANSION = 5 * K
TELOMERE_COMPRESS_RANGE = 100*K
CENSAT_COMPRESSABLE_THRESHOLD = 1000*K
VIRTUAL_TELOMERE_FLANK = 1 * K
VIRTUAL_TELOMERE_NODE_PREFIX = "virtual_telo_"

FORCE_TELOMERE_THRESHOLD = 10*K
TELOMERE_CLUSTER_THRESHOLD = 500*K
SUBTELOMERE_LENGTH = 500*K
MULTI_END_ALIGNMENT_WINDOW = 500*K

CENSAT_OUT_DIFF_RATIO = 0.30

MIN_FLANK_SIZE_BP = 1*M

REPEAT_MERGE_GAP = 0

MAX_OVERLAP_SCORE = 3

SPLIT_CTG_LEN_LIMIT = 100 * K
TRUST_SPLIT_CTG_LEN_LIMIT = 20 * K

TYPE2_CHUKJI_AS_TYPE4 = 5 * M
TYPE2_CONJOIN_COMPRESS_LIMIT = 1 * M
TYPE2_DIST_FLIP_THRESHOLD = 100 * K

RAW_TRANSLOCATION_CANDIDATE_PKL = 'raw_translocation_candidates.pkl'
RAW_TRANSLOCATION_RESULT_PKL = 'raw_translocation_result.pkl'
RAW_TRANSLOCATION_WINDOW = 5 * K
RAW_TRANSLOCATION_MIN_SAME_CHROM_SPAN = 10 * M

NCLOSE_COUNT_CANDIDATE_PKL = 'nclose_count_candidates.pkl'
NCLOSE_COUNT_RESULT_PKL = 'nclose_count_result.pkl'
NCLOSE_COUNT_DEFAULT_VAF_THRESHOLD = 0.1

STAGE01_SUMMARY_JSON = "stage01_nclose_summary.json"
STAGE01_REJECTIONS_TSV = "stage01_nclose_rejections.tsv"

ALT_SIMPLE_MIN_SEGMENT_LEN = 10 * K
ALT_SIMPLE_MAJOR_CHR_RATIO = 0.90
ALT_SIMPLE_EXISTING_NCLOSE_DIST = 10 * K

JULIA_BAM_THREAD_LIM = 32


class PregraphSourceMode(Enum):
    """Input layouts supported by the stage-01 pregraph builder."""

    CONFIGURED_PAF = "configured_paf"
    PRIMARY_ONLY_RETRY = "primary_only_retry"
    VCF = "vcf"


class PafSourceKind(Enum):
    """PAF source roles with intentionally different legacy preprocessing."""

    PRIMARY = "primary"
    SECONDARY = "secondary"


class CensatPairClass(Enum):
    """Number of CEN-SAT-labelled endpoints in one path-ordered NClose."""

    NONE = 0
    ONE = 1
    BOTH = 2


@dataclass(frozen=True)
class NCloseSourceConfig:
    """Resolved PAF inputs for one pregraph build attempt."""

    mode: PregraphSourceMode
    paf_file_paths: tuple
    original_paf_paths: tuple
    is_unitig_reduced: bool
    secondary_candidate_paf: str | None


@dataclass(frozen=True)
class PregraphBuildContext:
    """Immutable inputs for one pregraph build attempt."""

    source: NCloseSourceConfig
    ori_ctg_name_data: list
    prefix: str
    preprocessed_paf_path: str
    reference_fai_path: str
    telomere_bed_path: str
    repeat_bed_path: str
    censat_bed_path: str
    main_stat_path: str
    asm2cov: object
    disable_alt_ctg_simple: bool


@dataclass(frozen=True)
class PafPreprocessPolicy:
    """Source-specific parts of the otherwise shared PAF preprocessing."""

    kind: PafSourceKind
    source_index: int
    global_index_prefix: str
    log_label: str


@dataclass
class PafPreprocessResult:
    """One source after the two legacy PAF preprocessing passes."""

    contigs: list
    kept_primary_names: set
    original_node_count: int
    excluded_telomere_origins: set


@dataclass
class ContigPreprocessResult:
    """Explicit outputs consumed by NClose discovery."""

    telo_coverage: object
    depth_df: object
    no_chrY: bool
    stage_records: tuple = ()


@dataclass(frozen=True)
class ContigPreprocessResources:
    """Reference annotations shared by all contig preprocessing stages."""

    chr_len: object
    telo_data: object
    telo_dict: object
    repeat_data: object
    repeat_censat_data: object
    cen_fragment_meta: object
    depth_df: object
    no_chrY: bool


@dataclass
class NCloseBuildResult:
    """Complete stage-01 state used to persist graph and downstream inputs."""

    df: object
    no_chrY: bool
    repeat_censat_data: object
    chr_len: object
    contig_data: list
    contig_data_size: int
    chr_corr: object
    chr_rev_corr: object
    telo_contig: object
    telo_node_count: int
    telo_set: set
    rpt_con: set
    bnd_contig: set
    raw_nclose_nodes: object
    nclose_nodes: object
    vctg_dict: dict
    all_nclose_comp: object
    uncomp_node_count: int
    nclose_node_count: int
    transloc_nclose_pair_count: int
    indel_exclude_idx_set: set
    telo_coverage: object
    contig_stage_records: tuple = ()
    nclose_stage_records: tuple = ()
    nclose_rejections: tuple = ()


@dataclass(frozen=True)
class NCloseCandidateBuildContext:
    """Inputs shared by primary-contig and unitig candidate discovery."""

    contig_data: list
    bnd_contig: set
    repeat_contig_names: set
    repeat_censat_data: object
    paf_file_paths: tuple
    original_paf_paths: tuple
    telo_set: set
    telo_contig: object
    chr_len: object
    original_contig_names: list


@dataclass(frozen=True)
class CensatNoncensatCandidate:
    """Path-aware view of an NClose with exactly one CEN-SAT endpoint."""

    contig_name: object
    pair: tuple[int, int]
    censat_idx: int
    noncensat_idx: int
    censat_side: int
    censat_chrom: str
    noncensat_chrom: str
    noncensat_pos: int
    censat_norm_dir: str
    is_simple_alt: bool


@dataclass
class NCloseCandidateBuildResult:
    """Unitig builder output, excluding removed st/ed compression maps."""

    candidates: tuple
    raw_nodes: object
    virtual_contigs: dict
    all_nodes: object


@dataclass(frozen=True)
class ContigPipelineContext:
    """Read-only inputs shared by ordered contig augmentation stages.

    A future splitter only needs to implement ``run(context, state)`` and can
    be inserted before or after a named stage.  Keeping the reference inputs
    here avoids adding another positional argument to the stage-01 driver for
    every new unitig rescue policy.
    """

    build: PregraphBuildContext
    resources: ContigPreprocessResources
    excluded_telomere_origins: frozenset
    original_node_count: int
    primary_kept_contig_names: frozenset


@dataclass
class ContigPipelineState:
    """Mutable contig collection passed between augmentation stages."""

    contigs: list
    telomere_bounds: dict
    stage_metadata: dict = field(default_factory=dict)
    stage_records: list = field(default_factory=list)


@dataclass(frozen=True)
class ContigPipelineStage:
    """One named, ordered unitig split/augmentation operation."""

    name: str
    run: object

    def __post_init__(self):
        if not self.name:
            raise ValueError("Contig pipeline stage name cannot be empty")
        if not callable(self.run):
            raise TypeError(f"Contig pipeline stage {self.name!r} is not callable")


@dataclass(frozen=True)
class NClosePipelineContext:
    """Read-only data available to every NClose candidate stage."""

    build: PregraphBuildContext
    source: NCloseSourceConfig
    contig_data: list
    repeat_censat_data: object
    chr_len: object
    cen_fragment_meta: object
    telo_contig: object
    raw_nclose_nodes: object
    all_nclose_nodes: object
    rejected_pairs: set
    rescued_pairs: set


@dataclass
class NClosePipelineState:
    """Candidate population and derived state passed between NClose stages."""

    candidates: list
    contig_data_size: int
    chr_corr: object
    chr_rev_corr: object
    indel_exclude_idx_set: set = field(default_factory=set)
    rejections: list = field(default_factory=list)
    stage_metadata: dict = field(default_factory=dict)
    stage_records: list = field(default_factory=list)


@dataclass(frozen=True)
class NClosePipelineStage:
    """One named candidate filter, augmentation, or artifact operation."""

    name: str
    run: object

    def __post_init__(self):
        if not self.name:
            raise ValueError("NClose pipeline stage name cannot be empty")
        if not callable(self.run):
            raise TypeError(f"NClose pipeline stage {self.name!r} is not callable")


def insert_pipeline_stage(stages, stage, *, before=None, after=None):
    """Return a stage tuple with ``stage`` inserted at one named boundary.

    This is the public composition primitive for future unitig splitters and
    NClose filters.  The caller chooses the biological ordering explicitly;
    no module-level registry or import side effect is involved.
    """

    if (before is None) == (after is None):
        raise ValueError("Specify exactly one of before= or after=")
    stages = tuple(stages)
    names = [item.name for item in stages]
    if len(names) != len(set(names)):
        raise ValueError(f"Pipeline contains duplicate stage names: {names}")
    if stage.name in names:
        raise ValueError(f"Pipeline stage already exists: {stage.name}")
    anchor = before if before is not None else after
    if anchor not in names:
        raise KeyError(f"Pipeline stage anchor does not exist: {anchor}")
    index = names.index(anchor) + (after is not None)
    return stages[:index] + (stage,) + stages[index:]


@dataclass(frozen=True)
class Stage01Pipeline:
    """Composable implementation policy for one stage-01 run."""

    candidate_builder: object
    contig_stages: tuple
    nclose_stages: tuple

    def with_contig_stage(self, stage, *, before=None, after=None):
        return replace(
            self,
            contig_stages=insert_pipeline_stage(
                self.contig_stages,
                stage,
                before=before,
                after=after,
            ),
        )

    def with_nclose_stage(self, stage, *, before=None, after=None):
        return replace(
            self,
            nclose_stages=insert_pipeline_stage(
                self.nclose_stages,
                stage,
                before=before,
                after=after,
            ),
        )


@dataclass(frozen=True)
class ExtractedNClose:
    """One path-ordered terminal pair before global NClose clustering."""

    contig_name: str
    path_pair: tuple[int, int]
    outer_pair: tuple[int, int]


@dataclass(frozen=True)
class NCloseClusterRepresentative:
    """The first extracted pair retained for one spatial cluster."""

    contig_name: str
    path_pair: tuple[int, int]
    stored_pair: tuple[int, int]
    canonical_directions: tuple[str, str]


@dataclass(frozen=True)
class Stage01Config:
    """CLI-independent inputs for one stage-01 preprocessing run."""

    paf_file_path: str
    reference_fai_path: str
    telomere_bed_path: str
    repeat_bed_path: str
    censat_bed_path: str
    main_stat_path: str
    prefix: str
    read_bam_path: str
    alt_path: str | None = None
    original_paf_paths: tuple[str, ...] = ()
    thread: int = 1
    progress: bool = False
    exclude_nclose_list_path: str | None = None
    skip_bam_analysis: bool = False
    check_nclose_count: bool = False
    nclose_count_vaf_threshold: float = NCLOSE_COUNT_DEFAULT_VAF_THRESHOLD
    disable_alt_ctg_simple: bool = False
    vcf_input_path: str | None = None
    vcf_filter_pass: tuple[str, ...] = ("PASS", ".")


def iter_contig_ranges(contig_data):
    """Yield the contiguous node-index range owned by each preprocessed contig."""

    start = 0
    while start < len(contig_data):
        end = int(contig_data[start][CTG_ENDND])
        yield start, end
        start = end + 1


def trim_nclose_terminal_pair(contig_data, start, end, max_gap=0):
    """Move each terminal inward across redundant terminal-state alignments.

    ``max_gap`` counts consecutive nonmatching states.  Stage 01 intentionally
    uses zero, so the first chromosome/direction change stops each scan.
    """

    start_state = (
        contig_data[start][CTG_DIR],
        contig_data[start][CHR_NAM],
    )
    end_state = (
        contig_data[end][CTG_DIR],
        contig_data[end][CHR_NAM],
    )

    gap = 0
    trimmed_start = start
    cursor = start + 1
    while cursor < end:
        state = (
            contig_data[cursor][CTG_DIR],
            contig_data[cursor][CHR_NAM],
        )
        if state == start_state:
            trimmed_start = cursor
            gap = 0
        else:
            gap += 1
            if gap > max_gap:
                break
        cut_ratio, _ = calculate_single_contig_ref_ratio(
            contig_data[start:trimmed_start + 1]
        )
        if abs(cut_ratio - 1) > BND_CONTIG_BOUND:
            break
        cursor += 1

    gap = 0
    trimmed_end = end
    cursor = end - 1
    while cursor > trimmed_start:
        state = (
            contig_data[cursor][CTG_DIR],
            contig_data[cursor][CHR_NAM],
        )
        if state == end_state:
            trimmed_end = cursor
            gap = 0
        else:
            gap += 1
            if gap > max_gap:
                break
        cut_ratio, _ = calculate_single_contig_ref_ratio(
            contig_data[trimmed_end:end + 1]
        )
        if abs(cut_ratio - 1) > BND_CONTIG_BOUND:
            break
        cursor -= 1

    return trimmed_start, trimmed_end


def extract_unitig_nclose_pairs(contig_data, max_gap=0):
    """Extract exactly one terminal pair from every post-trim type-1/2 unitig."""

    extracted = []
    all_nodes = defaultdict(list)
    for start, end in iter_contig_ranges(contig_data):
        if int(contig_data[start][CTG_TYP]) not in (1, 2):
            continue

        pair = trim_nclose_terminal_pair(
            contig_data,
            start,
            end,
            max_gap=max_gap,
        )
        if pair[0] >= pair[1]:
            logging.warning(
                "Skipping collapsed type-%s NClose owner %s: %s",
                contig_data[start][CTG_TYP],
                contig_data[start][CTG_NAM],
                pair,
            )
            continue

        contig_name = str(contig_data[start][CTG_NAM])
        extracted.append(
            ExtractedNClose(
                contig_name=contig_name,
                path_pair=pair,
                outer_pair=(start, end),
            )
        )
        chrom_pair = (
            contig_data[pair[0]][CHR_NAM],
            contig_data[pair[1]][CHR_NAM],
        )
        all_nodes[chrom_pair].append(pair)

    return tuple(extracted), all_nodes


def nclose_compression_layout(
    st_chr,
    ed_chr,
    st,
    ed,
    st_interval,
    ed_interval,
):
    """Return the chromosome bucket and canonical endpoint order."""

    st_chr_order = chr2int(st_chr[1])
    ed_chr_order = chr2int(ed_chr[1])
    same_chrom = st_chr_order == ed_chr_order
    store_path_order = (
        st_chr_order < ed_chr_order
        or (same_chrom and st_interval <= ed_interval)
    )
    if store_path_order:
        return (st_chr, ed_chr), (st, ed), same_chrom
    return (ed_chr, st_chr), (ed, st), same_chrom


def nclose_canonical_directions(
    contig_data,
    path_pair,
    stored_pair,
    canonical_order_tied=False,
):
    """Return endpoint directions in the canonical comparison frame."""

    path_directions = tuple(
        contig_data[node_idx][CTG_DIR]
        for node_idx in path_pair
    )
    if stored_pair == path_pair:
        canonical_directions = path_directions
    elif stored_pair == tuple(reversed(path_pair)):
        canonical_directions = tuple(
            _flip_ctg_dir(direction)
            for direction in reversed(path_directions)
        )
    else:
        raise ValueError(
            f"Stored NClose nodes {stored_pair} do not match path {path_pair}"
        )

    if not canonical_order_tied:
        return canonical_directions
    reverse_complement = tuple(
        _flip_ctg_dir(direction)
        for direction in reversed(canonical_directions)
    )
    return min(canonical_directions, reverse_complement)


def nclose_cluster_candidate_matches(
    contig_data,
    stored_pair,
    canonical_directions,
    representative,
    compress_limit,
):
    """Return whether both endpoints and the full direction signature match."""

    return (
        distance_checker(
            contig_data[stored_pair[0]],
            contig_data[representative.stored_pair[0]],
        ) < compress_limit
        and distance_checker(
            contig_data[stored_pair[1]],
            contig_data[representative.stored_pair[1]],
        ) < compress_limit
        and canonical_directions == representative.canonical_directions
    )


def _contained_censat_locus(node, repeat_censat_data):
    for interval_idx, interval in enumerate(
        repeat_censat_data[node[CHR_NAM]]
    ):
        if inclusive_checker_tuple(
            interval,
            (node[CHR_STR], node[CHR_END]),
        ):
            return f"{node[CHR_NAM]}.{interval_idx}"
    return None


def cluster_unitig_nclose_pairs(
    contig_data,
    extracted,
    repeat_contig_names,
    repeat_censat_data,
):
    """Apply global spatial clustering and the legacy CEN-SAT-locus dedup."""

    representatives_by_bucket = defaultdict(list)
    seen_censat_pairs = set()
    kept = []

    for candidate in extracted:
        path_pair = candidate.path_pair
        start_idx, end_idx = path_pair
        start_node = contig_data[start_idx]
        end_node = contig_data[end_idx]
        start_chr = ("=", start_node[CHR_NAM])
        end_chr = ("=", end_node[CHR_NAM])
        start_interval = (start_node[CHR_STR], start_node[CHR_END])
        end_interval = (end_node[CHR_STR], end_node[CHR_END])
        bucket, stored_pair, same_chrom = nclose_compression_layout(
            start_chr,
            end_chr,
            start_idx,
            end_idx,
            start_interval,
            end_interval,
        )
        canonical_directions = nclose_canonical_directions(
            contig_data,
            path_pair,
            stored_pair,
            canonical_order_tied=(
                same_chrom and start_interval == end_interval
            ),
        )

        representatives = representatives_by_bucket[bucket]
        duplicate = False
        for representative in representatives:
            if (
                candidate.contig_name in repeat_contig_names
                and representative.contig_name in repeat_contig_names
            ):
                compress_limit = ALL_REPEAT_NCLOSE_COMPRESS_LIMIT
            else:
                compress_limit = NCLOSE_COMPRESS_LIMIT
            if nclose_cluster_candidate_matches(
                contig_data,
                stored_pair,
                canonical_directions,
                representative,
                compress_limit,
            ):
                duplicate = True
                break
        if duplicate:
            continue

        outer_start = contig_data[candidate.outer_pair[0]]
        outer_end = contig_data[candidate.outer_pair[1]]
        if node_is_censat(outer_start) and node_is_censat(outer_end):
            start_locus = _contained_censat_locus(
                start_node,
                repeat_censat_data,
            )
            end_locus = _contained_censat_locus(
                end_node,
                repeat_censat_data,
            )
            if start_locus is not None and end_locus is not None:
                # Preserve the legacy path-ordered, direction-agnostic key.
                censat_key = (("=", start_locus), ("=", end_locus))
                if censat_key in seen_censat_pairs:
                    continue
                seen_censat_pairs.add(censat_key)

        representatives.append(
            NCloseClusterRepresentative(
                contig_name=candidate.contig_name,
                path_pair=path_pair,
                stored_pair=stored_pair,
                canonical_directions=canonical_directions,
            )
        )
        kept.append(candidate)

    return tuple(kept)


def build_unitig_nclose_candidates(context):
    """Extract and globally cluster the minimal assembly NClose population."""

    extracted, all_nodes = extract_unitig_nclose_pairs(
        context.contig_data,
        max_gap=0,
    )
    clustered = cluster_unitig_nclose_pairs(
        context.contig_data,
        extracted,
        context.repeat_contig_names,
        context.repeat_censat_data,
    )
    candidates = tuple(
        NCloseCandidate(
            contig_name=candidate.contig_name,
            path_pair=candidate.path_pair,
            origin="paf",
        )
        for candidate in clustered
    )
    raw_nodes = candidates_to_legacy(candidates)
    logging.info(
        "Minimal unitig NClose build: %d extracted pair(s), %d clustered pair(s)",
        len(extracted),
        len(clustered),
    )
    return NCloseCandidateBuildResult(
        candidates=candidates,
        raw_nodes=raw_nodes,
        virtual_contigs={},
        all_nodes=all_nodes,
    )

def import_data(file_path : str) -> list :
    contig_data = []
    int_induce_idx = [1, 2, 3, 6, 7, 8, 9]
    idx = 0
    with open(file_path, "r") as paf_file:
        for curr_contig in paf_file:
            a = curr_contig.split("\t")
            temp_list = a[:9]
            temp_list.append(a[11])
            for i in int_induce_idx:
                temp_list[i] = int(temp_list[i])
            temp_list.append(idx)
            contig_data.append(temp_list)
            idx+=1
    return contig_data


def terminal_alignment_anchors(row, window=MULTI_END_ALIGNMENT_WINDOW) -> set:
    """Return chromosome ends that fully contain this PAF alignment."""
    chrom = row[CHR_NAM]
    chrom_len = int(row[CHR_LEN])
    ref_st = int(row[CHR_STR])
    ref_nd = int(row[CHR_END])
    anchors = set()

    if ref_nd <= min(int(window), chrom_len):
        anchors.add((chrom, 'f'))
    if ref_st >= max(0, chrom_len - int(window)):
        anchors.add((chrom, 'b'))
    return anchors


def find_multi_end_aligned_contigs(contig_data, window=MULTI_END_ALIGNMENT_WINDOW):
    """Find contigs whose every PAF row lies within at least two chromosome ends."""
    row_indices_by_contig = defaultdict(list)
    anchors_by_contig = defaultdict(set)
    nonterminal_contigs = set()

    for row_idx, row in enumerate(contig_data):
        contig_name = row[CTG_NAM]
        row_indices_by_contig[contig_name].append(row_idx)
        anchors = terminal_alignment_anchors(row, window)
        if anchors:
            anchors_by_contig[contig_name].update(anchors)
        else:
            nonterminal_contigs.add(contig_name)

    excluded_contigs = {
        contig_name
        for contig_name, anchors in anchors_by_contig.items()
        if contig_name not in nonterminal_contigs and len(anchors) >= 2
    }
    excluded_row_indices = {
        row_idx
        for contig_name in excluded_contigs
        for row_idx in row_indices_by_contig[contig_name]
    }
    return excluded_contigs, excluded_row_indices


def import_data2(file_path : str) -> list :
    paf_file = open(file_path, "r")
    contig_data = []
    for curr_contig in paf_file:
        curr_contig = curr_contig.rstrip()
        temp_list = curr_contig.split("\t")
        int_induce_idx = [CTG_LEN, CTG_STR, CTG_END, \
                          CHR_LEN, CHR_STR, CHR_END, \
                          CTG_MAPQ, CTG_TYP, CTG_STRND, CTG_ENDND,]
        for i in int_induce_idx:
            temp_list[i] = int(temp_list[i])
        contig_data.append(tuple(temp_list))
    paf_file.close()
    return contig_data

def import_telo_data(file_path : str, chr_len : dict) -> dict :
    fai_file = open(file_path, "r")
    telo_data = []
    for curr_data in fai_file:
        temp_list = curr_data.split("\t")
        int_induce_idx = [1, 2]
        for i in int_induce_idx:
            temp_list[i] = int(temp_list[i])
        if temp_list[1]>chr_len[temp_list[0]]/2:
            temp_list[1]-=TELOMERE_EXPANSION
            temp_list.append('b')
        else:
            temp_list.append('f')
            temp_list[2]+=TELOMERE_EXPANSION
        telo_data.append(tuple(temp_list))
    fai_file.close()
    return telo_data


def import_repeat_data_00(file_path : str) -> dict :
    fai_file = open(file_path, "r")
    repeat_data = defaultdict(list)
    for curr_data in fai_file:
        temp_list = curr_data.split("\t")
        curr_patch = (int(temp_list[1]), int(temp_list[2]))
        if len(repeat_data[temp_list[0]])==0:
            repeat_data[temp_list[0]].append(curr_patch)
        else:
            if abs(curr_patch[0] - repeat_data[temp_list[0]][-1][1]) < REPEAT_MERGE_GAP:
                repeat_data[temp_list[0]][-1] = (repeat_data[temp_list[0]][-1][0], curr_patch[1])
            else:
                repeat_data[temp_list[0]].append(curr_patch)
    # for k, v in repeat_data.items():
    #     print(k, len(v))
    fai_file.close()
    return repeat_data

def import_censat_repeat_data(file_path : str) -> dict :
    fai_file = open(file_path, "r")
    repeat_data = defaultdict(list)
    for curr_data in fai_file:
        temp_list = curr_data.split("\t")
        ref_data = (int(temp_list[1]), int(temp_list[2]))
        if abs(ref_data[1] - ref_data[0]) > CENSAT_COMPRESSABLE_THRESHOLD:
            repeat_data[temp_list[0]].append(ref_data)
    fai_file.close()
    return repeat_data

def is_stepwise_nonoverlapping(intervals: list) -> bool:
    n = len(intervals)
    if n < 2:
        return True

    for i in range(1, n):
        st_i, nd_i       = intervals[i]
        st_prev, nd_prev = intervals[i-1]

        if st_i <= st_prev:
            return False

        if nd_i <= nd_prev:
            return False

        if i >= 2:
            _, nd_prev2 = intervals[i-2]
            if st_i < nd_prev2:
                return False

    return True

def is_alt_contained_in_segments(segments: list, alt_segments: list) -> bool:
    n = len(segments)
    for st_b, nd_b in alt_segments:
        ok = False
        for i, (st_a, nd_a) in enumerate(segments):
            cond_start = True if i == 0 else (st_a <= st_b)
            cond_end   = True if i == n-1 else (nd_b <= nd_a)

            if cond_start and cond_end:
                ok = True
                break
        if not ok:
            return False
    return True

def max_overlap(intervals, target_intervals):
    events = []
    temp_inf = 2 * len(intervals) if target_intervals else 0

    for l in target_intervals:
        st, nd = l[1], l[2]
        events.append((st, +temp_inf))
        events.append((nd, -temp_inf))

    for l in intervals:
        st, nd = l[1], l[2]
        events.append((st, +1))
        events.append((nd, -1))

    events.sort(key=lambda x: (x[0], x[1]))

    curr = max_cnt = 0
    for _, delta in events:
        curr += delta
        max_cnt = max(max_cnt, curr)

    return max_cnt - temp_inf

def get_qry_cord_data(paf_path: str, get_ori_cord: bool = False) -> tuple:
    end_idx_dict = dict()
    paf_data_list = []

    with open(paf_path, "r") as paf_file:
        for line in paf_file:
            cols = line.rstrip("\n").split("\t")
            # [contig, start, end, mapq]
            tar_data = [
                cols[0],
                int(cols[2]),
                int(cols[3]),
                int(cols[11])
            ]

            if get_ori_cord:
                tar_col = None
                for col in cols:
                    if col.startswith('xi:Z:'):
                        tar_col = col
                        break

                ori_type, ori_cord = tar_col.split(':')[-1].split('_')
                tar_data.append(int(ori_cord) if ori_type == 'P' else -1)

            paf_data_list.append(tar_data)

    if paf_data_list:
        current_contig = paf_data_list[0][0]
        start_idx = 0
        for i, entry in enumerate(paf_data_list):
            contig = entry[0]
            if contig != current_contig:
                end_idx_dict[current_contig] = (start_idx, i - 1)
                current_contig = contig
                start_idx = i

        end_idx_dict[current_contig] = (start_idx, len(paf_data_list) - 1)

    return paf_data_list, end_idx_dict

def div_repeat_paf(
    original_paf_path_list: list,
    aln_paf_path_list: list,
    contig_data: list,
    original_contig_names: list,
) -> set:
    not_using_contig = set()

    ori_paf_data_list_data = []
    ori_end_idx_dict_data = []

    for original_paf_path in original_paf_path_list:
        ori_paf_data_list, ori_end_idx_dict = get_qry_cord_data(original_paf_path)

        ori_paf_data_list_data.append(ori_paf_data_list)
        ori_end_idx_dict_data.append(ori_end_idx_dict)

    aln_paf_data_list_data = []
    aln_end_idx_dict_data = []

    for aln_paf_path in aln_paf_path_list:
        aln_paf_data_list, aln_end_idx_dict = get_qry_cord_data(aln_paf_path, get_ori_cord=True)

        aln_paf_data_list_data.append(aln_paf_data_list)
        aln_end_idx_dict_data.append(aln_end_idx_dict)

    ppc_aln_end_idx_dict = dict()

    s = 0
    contig_data_size = len(contig_data)
    while s < contig_data_size:
        e = contig_data[s][CTG_ENDND]
        ori_ctg_name = contig_data[s][CTG_NAM]

        s_cl_ind, s_c_ind = map(int, contig_data[s][CTG_GLOBALIDX].split('.'))
        e_cl_ind, e_c_ind = map(int, contig_data[e][CTG_GLOBALIDX].split('.'))

        assert(s_cl_ind == e_cl_ind)

        if s_cl_ind < 2:
            assert(original_contig_names[s_cl_ind][s_c_ind] == original_contig_names[e_cl_ind][e_c_ind])
            ppc_aln_end_idx_dict[ori_ctg_name] = (s_cl_ind, s_c_ind, e_c_ind)

        s = e+1

    for ctg_name, (cl_ind, ppc_st, ppc_nd) in ppc_aln_end_idx_dict.items():
        ori_paf_data_list = ori_paf_data_list_data[cl_ind]
        ori_end_idx_dict = ori_end_idx_dict_data[cl_ind]

        aln_paf_data_list = aln_paf_data_list_data[cl_ind]
        aln_end_idx_dict = aln_end_idx_dict_data[cl_ind]

        ori_ctg_name = original_contig_names[cl_ind][ppc_st]

        ori_st, ori_nd = ori_end_idx_dict[ori_ctg_name]
        aln_st, aln_nd = aln_end_idx_dict[ori_ctg_name]
        assert(aln_st <= ppc_st and ppc_nd <= aln_nd)

        overlap_score = max_overlap(ori_paf_data_list[ori_st:ori_nd+1], [aln_paf_data_list[ppc_st], aln_paf_data_list[ppc_nd]])

        is_strong_paf = False

        ori_aln_ind = []
        for aln_ind in range(ppc_st, ppc_nd + 1):
            ori_ind = aln_paf_data_list[aln_ind][-1]

            if ori_ind != -1:
                ori_aln_ind.append(ori_ind)

        ori_aln_ind_set = set(ori_aln_ind)

        ori_aln_intervals = []
        ori_alt_intervals = []
        for ori_ind in range(ori_st, ori_nd + 1):
            _, st, nd, _ = ori_paf_data_list[ori_ind]

            if ori_ind in ori_aln_ind_set:
                ori_aln_intervals.append((st, nd))
            else:
                ori_alt_intervals.append((st, nd))

        if is_stepwise_nonoverlapping(ori_aln_intervals):
            is_strong_paf = is_alt_contained_in_segments(ori_aln_intervals, ori_alt_intervals)

        if overlap_score >= MAX_OVERLAP_SCORE and not is_strong_paf:
            not_using_contig.add(ctg_name)

    return not_using_contig

def get_overlap_bothend_score_dict(
    original_paf_path_list: list,
    aln_paf_path_list: list,
    contig_data: list,
    original_contig_names: list,
) -> dict:
    ori_paf_data_list_data = []
    ori_end_idx_dict_data = []

    for original_paf_path in original_paf_path_list:
        ori_paf_data_list, ori_end_idx_dict = get_qry_cord_data(original_paf_path)

        ori_paf_data_list_data.append(ori_paf_data_list)
        ori_end_idx_dict_data.append(ori_end_idx_dict)

    aln_paf_data_list_data = []
    aln_end_idx_dict_data = []

    for aln_paf_path in aln_paf_path_list:
        aln_paf_data_list, aln_end_idx_dict = get_qry_cord_data(aln_paf_path, get_ori_cord=True)

        aln_paf_data_list_data.append(aln_paf_data_list)
        aln_end_idx_dict_data.append(aln_end_idx_dict)

    ppc_aln_end_idx_dict = dict()

    s = 0
    contig_data_size = len(contig_data)
    while s < contig_data_size:
        e = contig_data[s][CTG_ENDND]
        ori_ctg_name = contig_data[s][CTG_NAM]

        s_cl_ind, s_c_ind = map(int, contig_data[s][CTG_GLOBALIDX].split('.'))
        e_cl_ind, e_c_ind = map(int, contig_data[e][CTG_GLOBALIDX].split('.'))

        assert(s_cl_ind == e_cl_ind)

        if s_cl_ind < 2:
            assert(original_contig_names[s_cl_ind][s_c_ind] == original_contig_names[e_cl_ind][e_c_ind])
            ppc_aln_end_idx_dict[ori_ctg_name] = (s_cl_ind, s_c_ind, e_c_ind)
        s = e+1

    ctgname2overlap = dict()
    for ctg_name, (cl_ind, ppc_st, ppc_nd) in ppc_aln_end_idx_dict.items():
        ori_paf_data_list = ori_paf_data_list_data[cl_ind]
        ori_end_idx_dict = ori_end_idx_dict_data[cl_ind]

        aln_paf_data_list = aln_paf_data_list_data[cl_ind]
        aln_end_idx_dict = aln_end_idx_dict_data[cl_ind]

        ori_ctg_name = original_contig_names[cl_ind][ppc_st]

        ori_st, ori_nd = ori_end_idx_dict[ori_ctg_name]
        aln_st, aln_nd = aln_end_idx_dict[ori_ctg_name]
        assert(aln_st <= ppc_st and ppc_nd <= aln_nd)

        overlap_score = max_overlap(ori_paf_data_list[ori_st:ori_nd+1], [aln_paf_data_list[ppc_st], aln_paf_data_list[ppc_nd]])
        ctgname2overlap[ctg_name] = overlap_score

    return ctgname2overlap

def get_overlap_total_score_dict(original_paf_path_list: list) -> dict:
    ori_paf_data_list_data = []
    ori_end_idx_dict_data = []

    for original_paf_path in original_paf_path_list:
        ori_paf_data_list, ori_end_idx_dict = get_qry_cord_data(original_paf_path)

        ori_paf_data_list_data.append(ori_paf_data_list)
        ori_end_idx_dict_data.append(ori_end_idx_dict)

    ctgname2overlap = dict()
    for cl_ind, ori_end_idx_dict in enumerate(ori_end_idx_dict_data):
        for ori_ctg_name in ori_end_idx_dict.keys():
            ori_paf_data_list = ori_paf_data_list_data[cl_ind]
            ori_end_idx_dict = ori_end_idx_dict_data[cl_ind]

            ori_st, ori_nd = ori_end_idx_dict[ori_ctg_name]

            overlap_score = max_overlap(ori_paf_data_list[ori_st:ori_nd+1], [])
            ctgname2overlap[ori_ctg_name] = overlap_score

    return ctgname2overlap

def get_not_trust_contig_name(original_paf_path_list: list) -> set:
    not_using_contig = set()
    for file_path in original_paf_path_list:
        with open(file_path, "r") as paf_file:
            for curr_contig in paf_file:
                a = curr_contig.split("\t")
                flag = None
                for i in a:
                    if i.startswith("tp:A:"):
                        if i[5] == 'P':
                            flag = False
                        else:
                            flag = True
                        break
                assert(flag is not None)

                if flag:
                    not_using_contig.add(a[0])
                if int(a[11]) < 60:
                    not_using_contig.add(a[0])

    return not_using_contig

def find_chr_len(file_path : str) -> dict:
    chr_data_file = open(file_path, "r")
    chr_len = {}
    for curr_data in chr_data_file:
        curr_data = curr_data.split("\t")
        chr_len[curr_data[0]] = int(curr_data[1])
    chr_data_file.close()
    return chr_len

def extract(contig : list) -> list:
    return [contig[CTG_NAM], contig[CHR_STR], contig[CHR_END]]

def chr2int(x):
    chrXY2int = {'chrX' : 24, 'chrY' : 25}
    if x in chrXY2int:
        return chrXY2int[x]
    else:
        return int(x[3:])

def distance_checker(node_a : tuple, node_b : tuple) -> int :
    if max(int(node_a[CHR_STR]), int(node_b[CHR_STR])) < min(int(node_a[CHR_END]), int(node_b[CHR_END])):
        return 0
    else:
        return min(abs(int(node_b[CHR_STR]) - int(node_a[CHR_END])), abs(int(node_b[CHR_END]) - int(node_a[CHR_STR])))

def distance_checker_tuple(node_a : tuple, node_b : tuple) -> int :
    if max(int(node_a[0]), int(node_b[0])) < min(int(node_a[1]), int(node_b[1])):
        return 0
    else:
        return min(abs(int(node_b[0]) - int(node_a[1])), abs(int(node_b[1]) - int(node_a[0])))

def alt_simple_ref_len(node) -> int:
    return abs(int(node[CHR_END]) - int(node[CHR_STR]))

def alt_simple_chr_sort_key(chrom: str):
    try:
        return chr2int(chrom)
    except (ValueError, TypeError):
        return INF

def alt_simple_terminal_repeat(node, repeat_label: tuple, chr_len: dict) -> bool:
    if repeat_label[1] not in ("r", "rin"):
        return False
    curr_chr_len = chr_len.get(node[CHR_NAM])
    if curr_chr_len is None:
        return False
    return node[CHR_STR] <= TELOMERE_EXPANSION or node[CHR_END] >= curr_chr_len - TELOMERE_EXPANSION

def trim_alt_simple_terminal_indices(chunks: list, telo_labels: list, repeat_labels: list, chr_len: dict) -> list:
    def is_terminal_chunk(idx: int) -> bool:
        return telo_labels[idx][0] != '0' or alt_simple_terminal_repeat(chunks[idx], repeat_labels[idx], chr_len)

    left = 0
    right = len(chunks) - 1
    while left <= right and is_terminal_chunk(left):
        left += 1
    while right >= left and is_terminal_chunk(right):
        right -= 1
    return list(range(left, right + 1))

def alt_simple_terminal_censat_included(trimmed_indices: list, censat_labels: list) -> bool:
    if not trimmed_indices:
        return False
    terminal_indices = (trimmed_indices[0], trimmed_indices[-1])
    return any(censat_labels[idx][1] == "rin" for idx in terminal_indices)

def select_alt_simple_major_chroms(indexed_chunks: list, required_chroms=()) -> set:
    chrom_len = Counter()
    total_len = 0
    for _, node in indexed_chunks:
        node_len = alt_simple_ref_len(node)
        if node_len <= ALT_SIMPLE_MIN_SEGMENT_LEN:
            continue
        chrom_len[node[CHR_NAM]] += node_len
        total_len += node_len

    if total_len == 0:
        return set(required_chroms)

    selected = set()
    acc_len = 0
    for chrom, chrom_span in sorted(chrom_len.items(), key=lambda x: (-x[1], alt_simple_chr_sort_key(x[0]), x[0])):
        selected.add(chrom)
        acc_len += chrom_span
        if acc_len / total_len >= ALT_SIMPLE_MAJOR_CHR_RATIO:
            break
    selected.update(required_chroms)
    return selected

def alt_simple_alignment_state(node) -> tuple:
    return (node[CHR_NAM], node[CTG_DIR])

def find_alt_simple_inward_bounds(indexed_chunks: list, selected_chroms: set):
    """
    Narrow both ends while preserving the chromosome/strand states observed
    immediately after terminal telomere/repeat trimming.

    Short (<=10 kb) and non-major-chromosome chunks are invisible to boundary
    discovery.  They are not deleted: the caller keeps every original chunk
    between the returned inclusive indices.
    """
    if len(indexed_chunks) < 2:
        return None

    left_terminal = indexed_chunks[0]
    right_terminal = indexed_chunks[-1]
    left_state = alt_simple_alignment_state(left_terminal[1])
    right_state = alt_simple_alignment_state(right_terminal[1])
    if left_state[0] == right_state[0]:
        return None

    left_bound = left_terminal
    for curr in indexed_chunks[1:]:
        node = curr[1]
        if alt_simple_ref_len(node) <= ALT_SIMPLE_MIN_SEGMENT_LEN:
            continue
        if node[CHR_NAM] not in selected_chroms:
            continue
        if alt_simple_alignment_state(node) != left_state:
            break
        left_bound = curr

    right_bound = right_terminal
    for curr in reversed(indexed_chunks[:-1]):
        node = curr[1]
        if alt_simple_ref_len(node) <= ALT_SIMPLE_MIN_SEGMENT_LEN:
            continue
        if node[CHR_NAM] not in selected_chroms:
            continue
        if alt_simple_alignment_state(node) != right_state:
            break
        right_bound = curr

    if left_bound[0] >= right_bound[0]:
        return None
    if alt_simple_alignment_state(left_bound[1]) != left_state:
        return None
    if alt_simple_alignment_state(right_bound[1]) != right_state:
        return None
    return left_bound, right_bound

def find_alt_simple_transition_candidates(indexed_chunks: list, selected_chroms: set) -> tuple:
    same_dir_candidates = []
    diff_dir_transitions = []
    prev = None
    for chunk_idx, node in indexed_chunks:
        if node[CHR_NAM] not in selected_chroms:
            continue
        if alt_simple_ref_len(node) <= ALT_SIMPLE_MIN_SEGMENT_LEN:
            continue

        curr = (chunk_idx, node)
        if prev is None:
            prev = curr
            continue

        prev_node = prev[1]
        if prev_node[CHR_NAM] == node[CHR_NAM]:
            prev = curr
            continue

        if prev_node[CTG_DIR] == node[CTG_DIR]:
            same_dir_candidates.append((prev, curr))
        else:
            diff_dir_transitions.append((prev, curr))
        prev = curr
    return same_dir_candidates, diff_dir_transitions

def alt_simple_node_locus(node) -> tuple:
    return (node[CHR_NAM], int(node[CHR_STR]), int(node[CHR_END]))

def alt_simple_pair_loci(candidate: tuple) -> tuple:
    return (alt_simple_node_locus(candidate[0][1]), alt_simple_node_locus(candidate[1][1]))

def alt_simple_locus_close(locus_a: tuple, locus_b: tuple, max_dist: int) -> bool:
    if locus_a[0] != locus_b[0]:
        return False
    return distance_checker_tuple((locus_a[1], locus_a[2]), (locus_b[1], locus_b[2])) <= max_dist

def alt_simple_pair_close(pair_a: tuple, pair_b: tuple, max_dist: int) -> bool:
    return (
        alt_simple_locus_close(pair_a[0], pair_b[0], max_dist)
        and alt_simple_locus_close(pair_a[1], pair_b[1], max_dist)
    )

def alt_simple_candidate_near_existing(candidate: tuple, existing_pairs: list, max_dist: int) -> bool:
    cand_pair = alt_simple_pair_loci(candidate)
    cand_pair_rev = (cand_pair[1], cand_pair[0])
    for existing_pair in existing_pairs:
        if alt_simple_pair_close(cand_pair, existing_pair, max_dist):
            return True
        if alt_simple_pair_close(cand_pair_rev, existing_pair, max_dist):
            return True
    return False

def collect_existing_alt_simple_nclose_loci(contig_data: list) -> list:
    existing_pairs = []
    st = 0
    while st < len(contig_data):
        ed = int(contig_data[st][CTG_ENDND])
        curr_contig = contig_data[st:ed + 1]
        st = ed + 1
        if not curr_contig:
            continue
        if curr_contig[0][CTG_NAM].startswith("simple_ctg_alt_"):
            continue
        try:
            curr_type = int(curr_contig[0][CTG_TYP])
        except (TypeError, ValueError):
            continue
        if curr_type not in (1, 2):
            continue

        indexed_chunks = list(enumerate(curr_contig))
        long_chunks = [
            (idx, node)
            for idx, node in indexed_chunks
            if alt_simple_ref_len(node) > ALT_SIMPLE_MIN_SEGMENT_LEN
        ]
        if len(long_chunks) >= 2 and long_chunks[0][1][CHR_NAM] != long_chunks[-1][1][CHR_NAM]:
            existing_pairs.append((
                alt_simple_node_locus(long_chunks[0][1]),
                alt_simple_node_locus(long_chunks[-1][1]),
            ))

        selected_chroms = select_alt_simple_major_chroms(indexed_chunks)
        if len(selected_chroms) < 2:
            continue
        same_dir_candidates, _ = find_alt_simple_transition_candidates(indexed_chunks, selected_chroms)
        for candidate in same_dir_candidates:
            existing_pairs.append(alt_simple_pair_loci(candidate))

    return existing_pairs

def overlap_calculator(node_a : tuple, node_b : tuple) -> int :
    return min(abs(node_a[CHR_END] - node_b[CHR_STR]), abs(node_b[CHR_END] - node_a[CHR_STR]))

def inclusive_checker_tuple(tuple_a : tuple, tuple_b : tuple) -> bool :
    if int(tuple_a[0]) <= int(tuple_b[0]) and int(tuple_b[1]) <= int(tuple_a[1]):
        return True
    else:
        return False

def chr_correlation_maker(contig_data):
    chr_corr = {}
    chr_rev_corr = {}
    contig_data_size = len(contig_data)
    for i in range(1, CHROMOSOME_COUNT):
        chr_corr['chr'+str(i)+'f'] = contig_data_size + i - 1
        chr_rev_corr[contig_data_size + i - 1] = 'chr'+str(i)+'f'
    chr_corr['chrXf'] = contig_data_size + CHROMOSOME_COUNT - 1
    chr_corr['chrYf'] = contig_data_size + CHROMOSOME_COUNT - 1
    chr_rev_corr[contig_data_size + CHROMOSOME_COUNT - 1] = 'chrXf'
    for i in range(1, CHROMOSOME_COUNT):
        chr_corr['chr'+str(i)+'b'] = contig_data_size + CHROMOSOME_COUNT + i - 1
        chr_rev_corr[contig_data_size + CHROMOSOME_COUNT + i - 1] = 'chr'+str(i)+'b'
    chr_corr['chrXb'] = contig_data_size + 2*CHROMOSOME_COUNT - 1
    chr_corr['chrYb'] = contig_data_size + 2*CHROMOSOME_COUNT - 1
    chr_rev_corr[contig_data_size + 2*CHROMOSOME_COUNT - 1] = 'chrXb'

    return chr_corr, chr_rev_corr

def extract_all_repeat_contig(contig_data : list, repeat_data : dict, ctg_index : int, baseline : float = 0) -> set:
    contig_data_size = len(contig_data)
    rpt_con = set()
    ends_map = {
        chrom: [iv[1] for iv in intervals]
        for chrom, intervals in repeat_data.items()
    }
    s = 0
    while s<contig_data_size:
        e = contig_data[s][CTG_ENDND]
        flag = True
        terminal_repeat = True
        total_ref_len = 0
        non_rpt_ref_len = 0
        for i in range(s, e+1):
            total_ref_len += contig_data[i][CHR_END] - contig_data[i][CHR_STR]
            if contig_data[i][ctg_index] == '0':
                non_rpt_ref_len += contig_data[i][CHR_END] - contig_data[i][CHR_STR]

        str_dir = contig_data[s][CTG_DIR]
        str_ref = contig_data[s][CHR_STR if str_dir == '+' else CHR_END]
        str_chr = contig_data[s][CHR_NAM]
        end_dir = contig_data[e][CTG_DIR]
        end_ref = contig_data[e][CHR_END if end_dir == '+' else CHR_STR]
        end_chr = contig_data[e][CHR_NAM]

        for ref, chrom in ((str_ref, str_chr), (end_ref, end_chr)):
            intervals = repeat_data.get(chrom, [])
            ends = ends_map.get(chrom, [])
            if not intervals:
                continue
            c_start, c_end = ref, ref
            idx = bisect.bisect_left(ends, c_start)
            if idx < len(intervals):
                iv_start, iv_end = intervals[idx]
                if iv_start > ref or ref > iv_end:
                    terminal_repeat = False

        if ctg_index == CTG_CENSAT:
            terminal_repeat = True
        if non_rpt_ref_len == 0 or baseline > non_rpt_ref_len/total_ref_len:
            if terminal_repeat:
                rpt_con.add(contig_data[s][CTG_NAM])
            # else:
            #     print(contig_data[s][CTG_NAM])
        s = e+1
    return rpt_con

def extract_bnd_contig(contig_data : list) -> set:
    s = 0
    contig_data_size = len(contig_data)
    bnd_contig = set()
    while s<contig_data_size:
        e = contig_data[s][CTG_ENDND]
        if contig_data[s][CTG_TYP] in {1, 2}: # 2 넣기
            bnd_contig.add(contig_data[s][CTG_NAM])
        s = e+1
    return bnd_contig

def calculate_single_contig_ref_ratio(contig_data : list) -> tuple:
    total_ref_len = 0
    if contig_data[0][CTG_DIR] == '+':
        estimated_ref_len = contig_data[-1][CHR_END] - contig_data[0][CHR_STR]
    else:
        estimated_ref_len = contig_data[0][CHR_END] - contig_data[-1][CHR_STR]
    for node in contig_data:
        total_ref_len += node[CHR_END] - node[CHR_STR]
    return estimated_ref_len/total_ref_len, total_ref_len

# 00_contig_preprocessing
def inclusive_checker_node(contig_node : tuple, telomere_node : tuple) -> bool :
    if int(telomere_node[CHR_STR]) <= int(contig_node[CHR_STR]) and int(contig_node[CHR_END]) <= int(telomere_node[CHR_END]):
        return True
    else:
        return False

def telo_distance_checker(node: tuple, telo: tuple) -> int :
    return min(abs(telo[CHR_STR] - node[CHR_END]), abs(telo[CHR_END] - node[CHR_STR]))

def telo_distance_checker_cord(node_st, node_nd, telo_st, telo_nd) -> int :
    return min(abs(telo_st - node_nd), abs(telo_nd - node_st))

def label_node(contig_data : list, telo_data) -> list :
    label = []
    contig_data_size = len(contig_data)
    for i in range(contig_data_size):
        checker = 0
        for j in telo_data[contig_data[i][CHR_NAM]]:
            if distance_checker(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0], j[1])) == 0:
                inclusive_label = ""
                if inclusive_checker_node(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0], j[1])):
                    inclusive_label = "in"
                label.append((contig_data[i][CHR_NAM], j[2]+inclusive_label))
                checker = 1
                break
        if checker==0:
            label.append(('0', '0'))
    return label

def label_subtelo_node(contig_data : list, telo_data) -> list :
    label = []
    contig_data_size = len(contig_data)
    for i in range(contig_data_size):
        checker = 0
        for j in telo_data[contig_data[i][CHR_NAM]]:
            if j[2]=='f':
                if distance_checker(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0], j[1]+SUBTELOMERE_LENGTH)) == 0:
                    inclusive_label = ""
                    if inclusive_checker_node(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0], j[1]+SUBTELOMERE_LENGTH)):
                        inclusive_label = "in"
                    label.append((contig_data[i][CHR_NAM], j[2]+inclusive_label))
                    checker = 1
                    break
            else:
                if distance_checker(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0] - SUBTELOMERE_LENGTH, j[1])) == 0:
                    inclusive_label = ""
                    if inclusive_checker_node(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0] - SUBTELOMERE_LENGTH, j[1])):
                        inclusive_label = "in"
                    label.append((contig_data[i][CHR_NAM], j[2]+inclusive_label))
                    checker = 1
                    break
        if checker==0:
            label.append(('0', '0'))
    return label

def label_repeat_node(contig_data: list, repeat_data : dict, chr_len : dict) -> list:
    labels = []
    ends_map = {
        chrom: [iv[1] for iv in intervals]
        for chrom, intervals in repeat_data.items()
    }

    for contig in contig_data:
        chrom = contig[CHR_NAM]
        intervals = repeat_data.get(chrom, [])
        curr_chr_len = chr_len.get(chrom, 0)
        ends = ends_map.get(chrom, [])
        if not intervals:
            labels.append(('0','0'))
            continue
        c_start, c_end = contig[7], contig[8]
        idx = bisect.bisect_left(ends, c_start)
        interval_list = []
        if SUBTELOMERE_REPEAT_LENGTH > 0:
            interval_list.append((0, SUBTELOMERE_REPEAT_LENGTH))
            interval_list.append((curr_chr_len - SUBTELOMERE_REPEAT_LENGTH, curr_chr_len))
        if idx < len(intervals):
            interval_list.append(intervals[idx])
        r_chk = False
        rin_chk = False
        for iv in interval_list:
            iv_start, iv_end = iv[0], iv[1]
            if distance_checker(contig, (0,0,0,0,0,0,0, iv_start, iv_end)) == 0:
                r_chk = True
                rin_chk = True if inclusive_checker_node(contig, (0,0,0,0,0,0,0, iv_start, iv_end)) else False
        if rin_chk:
            labels.append((chrom, "rin"))
        elif r_chk:
            labels.append((chrom, "r"))
        else:
            labels.append(('0','0'))

    return labels

def preprocess_telo(contig_data : list, node_label : list) -> tuple :
    telo_preprocessed_contig = []
    telo_connect_info = {}
    semi_telomere_contig = []
    report_case = {'A':[], 'B':[], 'C':[], 'ESC':[], 'ALL_TELO_NON_ESC':[]}
    contig_data_size = len(contig_data)
    curr_contig_st = 0
    contig_data.append((0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    for i in range(1, contig_data_size+1):
        if contig_data[i-1][CTG_NAM] != contig_data[i][CTG_NAM]:
            curr_contig_ed = i-1
            front_telo_bound = curr_contig_st
            end_telo_bound = curr_contig_ed
            while front_telo_bound<=curr_contig_ed and node_label[front_telo_bound][0] != '0' :
                front_telo_bound+=1
            while end_telo_bound>=curr_contig_st and node_label[end_telo_bound][0] != '0':
                end_telo_bound-=1
            st = curr_contig_st
            ed = curr_contig_ed
            all_telo = False
            escape = False
            if front_telo_bound == curr_contig_ed+1:
                all_telo = True
                if len(node_label[curr_contig_st][1])==1 \
                and (node_label[curr_contig_st][1] + contig_data[curr_contig_st][CTG_DIR] in ("b+", "f-")):
                    escape = True
                elif len(node_label[curr_contig_ed][1])==1 \
                and (node_label[curr_contig_ed][1] + contig_data[curr_contig_ed][CTG_DIR] in ("b-", "f+")):
                    escape = True
            else:
                escape = True
            if escape:
                if front_telo_bound > curr_contig_st:
                    front_telo_bound-=1
                    # Check if first telomere node is 'Nin' -> Else: escape
                    if len(node_label[curr_contig_st][1])==1 \
                    and (node_label[curr_contig_st][1] + contig_data[curr_contig_st][CTG_DIR] in ("b+", "f-")):
                        report_case['ESC'].append(curr_contig_st)
                        st = curr_contig_st
                    elif len(node_label[curr_contig_st][1])>1:
                        # Next node is connected with telomere node.
                        front_telo_bound+=1
                        if front_telo_bound <= curr_contig_ed:
                            dest = contig_data[front_telo_bound][CHR_NAM]
                            if contig_data[front_telo_bound][CTG_DIR] == '+':
                                dest += 'f'
                                telo_connect_info[front_telo_bound] = dest
                                report_case['A'].append([dest, front_telo_bound])
                            else:
                                dest += 'b'
                                telo_connect_info[front_telo_bound] = dest
                                report_case['A'].append([dest, front_telo_bound])
                        st = front_telo_bound
                    # If boundary node is not "Nin"
                    else:
                        if node_label[curr_contig_st][1] == 'f':
                            if contig_data[front_telo_bound][CTG_DIR]=='+':
                                # boundary node is connected with telomere node.
                                dest = contig_data[front_telo_bound][CHR_NAM]+'f'
                                telo_connect_info[front_telo_bound] = dest
                                report_case['B'].append([dest, front_telo_bound])
                            else:
                                # Treat as "Nin" Case
                                front_telo_bound+=1
                                if front_telo_bound <= curr_contig_ed:
                                    dest = contig_data[front_telo_bound][CHR_NAM]
                                    if contig_data[front_telo_bound][CTG_DIR] == '+':
                                        dest += 'f'
                                        telo_connect_info[front_telo_bound] = dest
                                        report_case['C'].append([dest, front_telo_bound])
                                    else:
                                        dest += 'b'
                                        telo_connect_info[front_telo_bound] = dest
                                        report_case['C'].append([dest, front_telo_bound])
                        else:
                            if contig_data[front_telo_bound][CTG_DIR]=='-':
                                dest = contig_data[front_telo_bound][CHR_NAM]+'b'
                                telo_connect_info[front_telo_bound] = dest
                                report_case['B'].append([dest, front_telo_bound])
                            else:
                                front_telo_bound+=1
                                if front_telo_bound <= curr_contig_ed:
                                    dest = contig_data[front_telo_bound][CHR_NAM]
                                    if contig_data[front_telo_bound][CTG_DIR] == '+':
                                        dest += 'f'
                                        telo_connect_info[front_telo_bound] = dest
                                        report_case['C'].append([dest, front_telo_bound])
                                    else:
                                        dest += 'b'
                                        telo_connect_info[front_telo_bound] = dest
                                        report_case['C'].append([dest, front_telo_bound])
                        st = front_telo_bound

                if end_telo_bound < curr_contig_ed:
                    end_telo_bound+=1
                    # Check if first telomere node is 'Nin' -> Else: escape
                    if len(node_label[curr_contig_ed][1])==1 \
                    and (node_label[curr_contig_ed][1] + contig_data[curr_contig_ed][CTG_DIR] in ("b-", "f+")):
                        report_case['ESC'].append(curr_contig_ed)
                        ed = curr_contig_ed
                    # If boundary node is 'Nin'
                    elif len(node_label[curr_contig_ed][1])>1:
                        # Next node is connected with telomere node.
                        end_telo_bound-=1
                        if end_telo_bound>=curr_contig_st:
                            dest = contig_data[end_telo_bound][CHR_NAM]
                            if contig_data[end_telo_bound][CTG_DIR] == '+':
                                dest += 'b'
                                telo_connect_info[end_telo_bound] = dest
                                report_case['A'].append([dest, end_telo_bound])
                            else:
                                dest += 'f'
                                telo_connect_info[end_telo_bound] = dest
                                report_case['A'].append([dest, end_telo_bound])
                        ed = end_telo_bound
                    # If boundary node is not "Nin"
                    else:
                        if node_label[curr_contig_ed][1] == 'b':
                            if contig_data[end_telo_bound][CTG_DIR]=='+':
                                # boundary node is connected with telomere node.
                                dest = contig_data[end_telo_bound][CHR_NAM]+'b'
                                telo_connect_info[end_telo_bound] = dest
                                report_case['B'].append([dest, end_telo_bound])
                            else:
                                # Treat as "Nin" Case
                                end_telo_bound+=1
                                if end_telo_bound>=curr_contig_st:
                                    dest = contig_data[end_telo_bound][CHR_NAM]
                                    if contig_data[end_telo_bound][CTG_DIR] == '+':
                                        dest += 'b'
                                        telo_connect_info[end_telo_bound] = dest
                                        report_case['C'].append([dest, end_telo_bound])
                                    else:
                                        dest += 'f'
                                        telo_connect_info[end_telo_bound] = dest
                                        report_case['C'].append([dest, end_telo_bound])
                        else:
                            if contig_data[end_telo_bound][CTG_DIR]=='-':
                                dest = contig_data[end_telo_bound][CHR_NAM]+'f'
                                telo_connect_info[end_telo_bound] = dest
                                report_case['B'].append([dest, end_telo_bound])
                            else:
                                end_telo_bound+=1
                                if end_telo_bound>=curr_contig_st:
                                    dest = contig_data[end_telo_bound][CHR_NAM]
                                    if contig_data[end_telo_bound][CTG_DIR] == '+':
                                        dest += 'b'
                                        telo_connect_info[end_telo_bound] = dest
                                        report_case['C'].append([dest, end_telo_bound])
                                    else:
                                        dest += 'f'
                                        telo_connect_info[end_telo_bound] = dest
                                        report_case['C'].append([dest, end_telo_bound])
                        ed = end_telo_bound
                if all_telo and st>ed:
                    if curr_contig_st <= st and st <= curr_contig_ed:
                        for j in range(curr_contig_st, st+1):
                            telo_preprocessed_contig.append(j)
                    elif curr_contig_st <= ed and ed <= curr_contig_ed:
                        for j in range(ed, curr_contig_ed+1):
                            telo_preprocessed_contig.append(j)
                else:
                    for j in range(st, ed+1):
                        telo_preprocessed_contig.append(j)

            else:
                report_case['ALL_TELO_NON_ESC'].append(curr_contig_st)
            # initialize
            curr_contig_st = curr_contig_ed+1
    contig_data = contig_data[0:-1]
    return telo_preprocessed_contig, report_case, telo_connect_info

def subtelo_cut(contig_data : list, node_label : list, subnode_label : list) -> list :
    telo_preprocessed_contig = []
    contig_data_size = len(contig_data)
    curr_contig_st = 0
    tcnt = 0
    idx = 0
    while curr_contig_st < contig_data_size:
        curr_contig_ed = contig_data[curr_contig_st][CTG_ENDND]
        front_telo_bound = curr_contig_st
        end_telo_bound = curr_contig_ed
        front_telcon = False
        end_telcon = False
        front_block_contain_telo = False
        end_block_contain_telo = False
        st = curr_contig_st
        ed = curr_contig_ed
        if contig_data[curr_contig_st][CTG_TELCON] != '0':
            front_telcon = True
        if contig_data[curr_contig_ed][CTG_TELCON] != '0':
            end_telcon = True
        front_dest = '0'
        end_dest = '0'
        if not front_telcon:
            while front_telo_bound<=curr_contig_ed and subnode_label[front_telo_bound][0] != '0' :
                if node_label[front_telo_bound][0] != '0':
                    front_block_contain_telo = True
                front_telo_bound+=1
        if not end_telcon:
            while end_telo_bound>=curr_contig_st and subnode_label[end_telo_bound][0] != '0':
                if node_label[end_telo_bound][0] != '0':
                    end_block_contain_telo = True
                end_telo_bound-=1
        for_cut = False
        bak_cut = False
        if front_telo_bound > curr_contig_st and front_block_contain_telo:
            st = front_telo_bound
            if contig_data[st][CTG_DIR]=='+':
                front_dest = contig_data[st][CHR_NAM] + 'f'
                for_cut = True
            else:
                front_dest = contig_data[st][CHR_NAM] + 'b'
                for_cut = True
        if end_telo_bound < curr_contig_ed and end_block_contain_telo:
            ed = end_telo_bound
            if contig_data[ed][CTG_DIR]=='+':
                end_dest = contig_data[ed][CHR_NAM] + 'b'
                bak_cut = True
            else:
                end_dest = contig_data[ed][CHR_NAM] + 'f'
                bak_cut = True
        if (for_cut or bak_cut) and st < ed:
            l = ed - st
            tcnt+=1
            if for_cut:
                temp_list = copy.deepcopy(contig_data[st])
                temp_list[CTG_NAM] = f'subtelomere_cut_contig_{tcnt}'
                temp_list[CTG_STRND] = idx
                temp_list[CTG_ENDND] = idx + l
                temp_list[CTG_TELCON] = front_dest
                telo_preprocessed_contig.append(temp_list)
            else:
                temp_list = copy.deepcopy(contig_data[st])
                temp_list[CTG_NAM] = f'subtelomere_cut_contig_{tcnt}'
                temp_list[CTG_STRND] = idx
                temp_list[CTG_ENDND] = idx + l
                telo_preprocessed_contig.append(temp_list)
            for i in range(st+1, ed):
                temp_list = copy.deepcopy(contig_data[i])
                temp_list[CTG_NAM] = f'subtelomere_cut_contig_{tcnt}'
                temp_list[CTG_STRND] = idx
                temp_list[CTG_ENDND] = idx + l
                telo_preprocessed_contig.append(temp_list)
            if bak_cut:
                temp_list = copy.deepcopy(contig_data[ed])
                temp_list[CTG_NAM] = f'subtelomere_cut_contig_{tcnt}'
                temp_list[CTG_STRND] = idx
                temp_list[CTG_ENDND] = idx + l
                temp_list[CTG_TELCON] = end_dest
                telo_preprocessed_contig.append(temp_list)
            else:
                temp_list = copy.deepcopy(contig_data[ed])
                temp_list[CTG_NAM] = f'subtelomere_cut_contig_{tcnt}'
                temp_list[CTG_STRND] = idx
                temp_list[CTG_ENDND] = idx + l
                telo_preprocessed_contig.append(temp_list)
            idx += l + 1

        curr_contig_st = curr_contig_ed+1
    return telo_preprocessed_contig


def initial_graph_build(contig_data : list, telo_data : dict, no_chrY : bool) -> list :
    '''
    Initialize
    '''
    contig_data_size = len(contig_data)
    chr_corr, chr_rev_corr = chr_correlation_maker(contig_data)

    adjacency = [[[] for _ in range(contig_data_size+CHROMOSOME_COUNT*2)], [[] for _ in range(contig_data_size+CHROMOSOME_COUNT*2)]]
    '''
    Algorithm
    '''
    curr_contig_st = 0
    while curr_contig_st<contig_data_size:
        curr_contig_ed = contig_data[curr_contig_st][CTG_ENDND]
        if curr_contig_st != curr_contig_ed:
            if contig_data[curr_contig_st][CTG_TELCON] != '0':
                dest = chr_corr[contig_data[curr_contig_st][CTG_TELCON]]
                adjacency[DIR_FOR][dest].append([DIR_FOR, curr_contig_st, 0])
                adjacency[DIR_FOR][curr_contig_st].append([DIR_FOR, dest, 0])
                adjacency[DIR_BAK][curr_contig_st].append([DIR_FOR, dest, 0])
            if contig_data[curr_contig_ed][CTG_TELCON] != '0':
                dest = chr_corr[contig_data[curr_contig_ed][CTG_TELCON]]
                adjacency[DIR_FOR][curr_contig_ed].append([DIR_FOR, dest, 0])
                adjacency[DIR_BAK][curr_contig_ed].append([DIR_FOR, dest, 0])
                adjacency[DIR_FOR][dest].append([DIR_BAK, curr_contig_ed, 0])
        else:
            if contig_data[curr_contig_st][CTG_TELCON] != '0':
                if contig_data[curr_contig_st][CTG_TELCON][-1]=='f':
                    dest = chr_corr[contig_data[curr_contig_st][CTG_TELCON]]
                    adjacency[DIR_FOR][dest].append([DIR_FOR, curr_contig_st, 0])
                    adjacency[DIR_BAK][curr_contig_st].append([DIR_FOR, dest, 0])
                    adjacency[DIR_FOR][curr_contig_st].append([DIR_FOR, dest, 0])
                else:
                    dest = chr_corr[contig_data[curr_contig_ed][CTG_TELCON]]
                    adjacency[DIR_FOR][curr_contig_ed].append([DIR_FOR, dest, 0])
                    adjacency[DIR_BAK][curr_contig_ed].append([DIR_FOR, dest, 0])
                    adjacency[DIR_FOR][dest].append([DIR_BAK, curr_contig_ed, 0])

        curr_contig_st = curr_contig_ed + 1
    '''
    If telomere node has 0 connection, connect with closest node with same chromosome type.
    '''
    end_node_list = set()
    for i in contig_data:
        end_node_list.add((i[CTG_STRND], i[CTG_ENDND]))

    now_telo_list = list(chr_corr.keys())

    if no_chrY:
        now_telo_list = [i for i in now_telo_list if 'chrY' not in i]

    for now_telo in now_telo_list:
        now_telo_chr = now_telo[:-1]
        curr_telo_set = set()
        i = chr_corr[now_telo]
        flag = False

        mini_telo_dist = INF
        for j in range(2):
            for connected_node in adjacency[j][i]:
                if contig_data[connected_node[1]][CHR_NAM] == now_telo_chr:
                    # 10K 이내면 없애기
                    if i < contig_data_size + CHROMOSOME_COUNT:
                        mini_telo_dist = min(mini_telo_dist, 0 if contig_data[connected_node[1]][CHR_STR] < telo_data[now_telo][0] else contig_data[connected_node[1]][CHR_STR] - telo_data[now_telo][0])
                        if contig_data[connected_node[1]][CHR_STR] < telo_data[now_telo][0]+FORCE_TELOMERE_THRESHOLD:
                            flag = True
                    else:
                        mini_telo_dist = min(mini_telo_dist, 0 if contig_data[connected_node[1]][CHR_END] > telo_data[now_telo][1] else telo_data[now_telo][1] - contig_data[connected_node[1]][CHR_END])
                        if contig_data[connected_node[1]][CHR_END] > telo_data[now_telo][1]-FORCE_TELOMERE_THRESHOLD:
                            flag = True
                    curr_telo_set.add(connected_node[1])

        if flag:
            continue

        telo_dist = mini_telo_dist
        telo_connect_node = INF
        telo_dir = 0
        temp_contig = (0, 0, 0, 0, 0, 0, 0, telo_data[now_telo][0], telo_data[now_telo][1])

        #우리가 연결하는 텔로미어 노드
        if now_telo[-1]=='f':
            for st, ed in end_node_list:
                # 텔로미어와 겹치지 않으며, 배제된 노드가 아니고, 중복이 아니며, + 방향이고, 거리가 갱신가능한 경우
                if contig_data[st][CHR_NAM] == now_telo_chr \
                and st not in curr_telo_set \
                and telo_dist > telo_distance_checker(contig_data[st], temp_contig) \
                and contig_data[st][CTG_DIR]=='+':
                    telo_connect_node = st
                    telo_dist = telo_distance_checker(contig_data[st], temp_contig)
                    telo_dir = '+'
                if contig_data[ed][CHR_NAM] == now_telo_chr \
                and ed not in curr_telo_set \
                and telo_dist > telo_distance_checker(contig_data[ed], temp_contig) \
                and contig_data[ed][CTG_DIR]=='-':
                    telo_connect_node = ed
                    telo_dist = telo_distance_checker(contig_data[ed], temp_contig)
                    telo_dir = '-'
            if telo_connect_node != INF:
                if telo_dir == '+':
                    adjacency[DIR_FOR][i].append([DIR_FOR, telo_connect_node, telo_dist])
                    adjacency[DIR_BAK][telo_connect_node].append([DIR_FOR, i, telo_dist])
                else:
                    adjacency[DIR_FOR][i].append([DIR_BAK, telo_connect_node, telo_dist])
                    adjacency[DIR_FOR][telo_connect_node].append([DIR_FOR, i, telo_dist])
        else:
            for st, ed in end_node_list:
                if contig_data[st][CHR_NAM] == now_telo_chr \
                and st not in curr_telo_set \
                and telo_dist > telo_distance_checker(contig_data[st], temp_contig) \
                and contig_data[st][CTG_DIR]=='-':
                    telo_connect_node = st
                    telo_dist = telo_distance_checker(contig_data[st], temp_contig)
                    telo_dir = '-'
                if contig_data[ed][CHR_NAM] == now_telo_chr \
                and ed not in curr_telo_set \
                and telo_dist > telo_distance_checker(contig_data[ed], temp_contig) \
                and contig_data[ed][CTG_DIR]=='+':
                    telo_connect_node = ed
                    telo_dist = telo_distance_checker(contig_data[ed], temp_contig)
                    telo_dir = '+'
            if telo_connect_node != INF:
                if telo_dir == '+':
                    adjacency[DIR_FOR][i].append([DIR_BAK, telo_connect_node, telo_dist])
                    adjacency[DIR_FOR][telo_connect_node].append([DIR_FOR, i, telo_dist])
                    adjacency[DIR_BAK][telo_connect_node].append([DIR_FOR, i, telo_dist])
                else:
                    adjacency[DIR_FOR][i].append([DIR_FOR, telo_connect_node, telo_dist])
                    adjacency[DIR_BAK][telo_connect_node].append([DIR_FOR, i, telo_dist])
                    adjacency[DIR_FOR][telo_connect_node].append([DIR_FOR, i, telo_dist])

    return adjacency

def edge_optimization(contig_data : list, contig_adjacency : list, telo_dict : dict,
                      asm2cov : dict, original_contig_names : list,
                      excluded_telomere_origins=None) -> tuple :

    contig_data_size = len(contig_data)
    excluded_telomere_origins = excluded_telomere_origins or set()
    excluded_candidate_nodes = set()
    chr_corr, chr_rev_corr = chr_correlation_maker(contig_data)
    contig_pair_nodes = defaultdict(list)
    telo_coverage = Counter()
    optimized_adjacency = [[[] for _ in range(contig_data_size + CHROMOSOME_COUNT*2)], [[] for _ in range(contig_data_size + CHROMOSOME_COUNT*2)]]
    for _ in range(2):
        for i in range(len(contig_data)+CHROMOSOME_COUNT*2):
            for edge in contig_adjacency[_][i]:
                if edge[1] >= contig_data_size or i >= contig_data_size:
                    optimized_adjacency[_][i].append(edge)
                    continue
                if contig_data[i][CTG_STRND]>contig_data[edge[1]][CTG_STRND]:
                    contig_pair_nodes[(contig_data[edge[1]][CTG_NAM], contig_data[i][CTG_NAM])] \
                    .append([edge[1], i])
                else:
                    contig_pair_nodes[(contig_data[i][CTG_NAM], contig_data[edge[1]][CTG_NAM])]\
                    .append([i, edge[1]])
    minmaxnode_contigpair = {}
    for pair in contig_pair_nodes:
        min_first_contig = INF
        max_first_contig = 0
        min_second_contig = INF
        max_second_contig = 0
        for edge in contig_pair_nodes[pair]:
            min_first_contig = min(min_first_contig, edge[0])
            max_first_contig = max(max_first_contig, edge[0])
            min_second_contig = min(min_second_contig, edge[1])
            max_second_contig = max(max_second_contig, edge[1])
            minmaxnode_contigpair[pair] = [min_first_contig, max_first_contig, min_second_contig, max_second_contig]

    contig_data_size = len(contig_data)
    for i in range(2):
        for j in range(contig_data_size):
            first_contig_name = contig_data[j][CTG_NAM]
            fcn_int = contig_data[j][CTG_STRND]
            for edge in contig_adjacency[i][j]:
                if edge[1] >= contig_data_size:
                    continue
                second_contig_name = contig_data[edge[1]][CTG_NAM]
                scn_int = contig_data[edge[1]][CTG_STRND]
                if fcn_int == scn_int:
                    optimized_adjacency[i][j].append(edge)
                elif fcn_int < scn_int:
                    if j in minmaxnode_contigpair[(first_contig_name, second_contig_name)][0:2]\
                    or edge[1] in minmaxnode_contigpair[(first_contig_name, second_contig_name)][2:4]:
                        optimized_adjacency[i][j].append(edge)
                else:
                    if edge[1] in minmaxnode_contigpair[(second_contig_name, first_contig_name)][0:2] \
                    or j in minmaxnode_contigpair[(second_contig_name, first_contig_name)][2:4]:
                        optimized_adjacency[i][j].append(edge)


    for i in range(contig_data_size, contig_data_size+2*CHROMOSOME_COUNT):
        telo_name = chr_rev_corr[i]
        telo_range = telo_dict[telo_name]
        for j in range(2):
            using_edge = []
            now_edge = [-1, 0, 0]
            curr_coverage = 0
            for edge in optimized_adjacency[j][i]:
                cl_ind, c_ind = map(int, contig_data[edge[1]][CTG_GLOBALIDX].split('.'))
                if (cl_ind, c_ind) in excluded_telomere_origins:
                    excluded_candidate_nodes.add(edge[1])
                    continue
                if cl_ind < 2:
                    name = original_contig_names[cl_ind][c_ind]
                else:
                    name = contig_data[edge[1]][CTG_NAM]
                    asm2cov[name] = -1
                if telo_name[-1]=='f':
                    if contig_data[edge[1]][CHR_STR]<=TELOMERE_CLUSTER_THRESHOLD:
                        if now_edge[0]<0:
                            curr_coverage += asm2cov[name]
                            now_edge = edge
                        else:
                            if contig_data[now_edge[1]][CHR_STR] > contig_data[edge[1]][CHR_STR]:
                                curr_coverage += asm2cov[name]
                                now_edge = edge
                            elif contig_data[now_edge[1]][CHR_STR] == contig_data[edge[1]][CHR_STR]:
                                if contig_data[now_edge[1]][CTG_LEN] > contig_data[edge[1]][CTG_LEN]:
                                    curr_coverage += asm2cov[name]
                                    now_edge = edge
                    else:
                        telo_compress_flag = False
                        for existing_edge in using_edge:
                            if distance_checker(contig_data[existing_edge[1]], contig_data[edge[1]]) < TELOMERE_COMPRESS_RANGE:
                                telo_compress_flag = True
                                telo_coverage[existing_edge[1]] += asm2cov[name]
                                break
                        if not telo_compress_flag:
                            telo_coverage[edge[1]] += asm2cov[name]
                            using_edge.append(edge)
                else:
                    if contig_data[edge[1]][CHR_END]>=telo_range[1]-TELOMERE_CLUSTER_THRESHOLD:
                        if now_edge[0]<0:
                            curr_coverage += asm2cov[name]
                            now_edge = edge
                        else:
                            if contig_data[now_edge[1]][CHR_END] < contig_data[edge[1]][CHR_END]:
                                curr_coverage += asm2cov[name]
                                now_edge = edge
                            elif contig_data[now_edge[1]][CHR_END] == contig_data[edge[1]][CHR_END]:
                                if contig_data[now_edge[1]][CTG_LEN] > contig_data[edge[1]][CTG_LEN]:
                                    curr_coverage += asm2cov[name]
                                    now_edge = edge
                    else:
                        telo_compress_flag = False
                        for existing_edge in using_edge:
                            if distance_checker(contig_data[existing_edge[1]], contig_data[edge[1]]) < TELOMERE_COMPRESS_RANGE:
                                telo_compress_flag = True
                                telo_coverage[existing_edge[1]] += asm2cov[name]
                                break
                        if not telo_compress_flag:
                            telo_coverage[edge[1]] += asm2cov[name]
                            using_edge.append(edge)
            if now_edge == [-1, 0, 0]:
                pass
            else:
                # telo_coverage[now_edge[1]] += curr_coverage # norm telomere coverage
                using_edge.append(now_edge)
            optimized_adjacency[j][i] = using_edge

    telo_connected_set = set()
    telo_connected_dict = dict()
    telo_connected_graph_dict = defaultdict(list)
    for j in range(2):
        for i in range(contig_data_size, contig_data_size + 2*(CHROMOSOME_COUNT)):
            for k in optimized_adjacency[j][i]:
                telo_connected_set.add(k[1])
                telo_connected_dict[k[1]] = chr_rev_corr[i]
                telo_connected_graph_dict[chr_rev_corr[i]].append(k)

    if excluded_candidate_nodes:
        logging.info(
            f"Skipped {len(excluded_candidate_nodes)} multi-end-aligned PAF nodes "
            "during telomere edge optimization"
        )

    return telo_connected_set, telo_connected_dict, telo_connected_graph_dict, telo_coverage

def _flip_ctg_dir(ctg_dir):
    return '-' if ctg_dir == '+' else '+'


def node_is_censat(node):
    """Return whether preprocessing labelled a node as overlapping CEN-SAT."""

    return node[CTG_CENSAT] != '0'


def classify_censat_pair(contig_data, pair):
    """Classify an NClose without losing its biologically meaningful path order."""

    censat_count = sum(node_is_censat(contig_data[node_idx]) for node_idx in pair)
    return CensatPairClass(censat_count)


def nclose_has_simple_alt(contig_name, pair, contig_data):
    """Recognize simple-alt evidence at either the owner or endpoint boundary."""

    return (
        str(contig_name).startswith('simple_ctg_alt_')
        or contig_data[pair[0]][CTG_NAM].startswith('simple_ctg_alt_')
        or contig_data[pair[1]][CTG_NAM].startswith('simple_ctg_alt_')
    )


def build_censat_noncensat_candidate(contig_data, contig_name, pair):
    """Build the shared one-CEN-SAT view used by repair and arbitration."""

    pair = tuple(pair)
    if classify_censat_pair(contig_data, pair) != CensatPairClass.ONE:
        return None

    start_is_censat = node_is_censat(contig_data[pair[0]])
    censat_side = 0 if start_is_censat else 1
    censat_idx = pair[censat_side]
    noncensat_idx = pair[1 - censat_side]
    censat_dir = contig_data[censat_idx][CTG_DIR]
    if censat_side == 1:
        censat_dir = _flip_ctg_dir(censat_dir)

    noncensat_node = contig_data[noncensat_idx]
    return CensatNoncensatCandidate(
        contig_name=contig_name,
        pair=pair,
        censat_idx=censat_idx,
        noncensat_idx=noncensat_idx,
        censat_side=censat_side,
        censat_chrom=contig_data[censat_idx][CHR_NAM],
        noncensat_chrom=noncensat_node[CHR_NAM],
        noncensat_pos=(noncensat_node[CHR_STR] + noncensat_node[CHR_END]) // 2,
        censat_norm_dir=censat_dir,
        is_simple_alt=nclose_has_simple_alt(contig_name, pair, contig_data),
    )


def group_censat_noncensat_candidates(candidates):
    """Cluster one-CEN-SAT candidates by chromosome pair and nearby offset."""

    grouped = defaultdict(list)
    for candidate in candidates:
        grouped[(candidate.censat_chrom, candidate.noncensat_chrom)].append(
            candidate
        )

    clusters = []
    for items in grouped.values():
        items.sort(key=lambda candidate: candidate.noncensat_pos)
        start = 0
        while start < len(items):
            end = start + 1
            while (
                end < len(items)
                and items[end].noncensat_pos - items[end - 1].noncensat_pos
                < OFFSET_DIR_GROUP_LIMIT
            ):
                end += 1
            clusters.append(items[start:end])
            start = end
    return clusters


def iter_nclose_owner_pairs(nclose_source):
    """Yield owner/pair occurrences from candidates or the legacy mapping."""

    if hasattr(nclose_source, "items"):
        for contig_name, pair_list in nclose_source.items():
            for pair in pair_list:
                yield contig_name, tuple(pair)
        return

    for candidate in nclose_source:
        yield candidate.contig_name, candidate.path_pair


def iter_censat_noncensat_candidates(contig_data, nclose_source):
    """Yield live one-CEN-SAT views without changing candidate order."""

    for contig_name, pair in iter_nclose_owner_pairs(nclose_source):
        candidate = build_censat_noncensat_candidate(
            contig_data,
            contig_name,
            pair,
        )
        if candidate is not None:
            yield candidate


def _cen_fragment_target_dir_from_meta(cen_fragment_meta, chrom):
    return '-' if cen_fragment_meta[chrom]['dir'] else '+'

def _normalized_telo_censat_dir(telo_name, ctg_dir):
    if telo_name[-1] == 'f':
        return ctg_dir
    return _flip_ctg_dir(ctg_dir)

def filter_telomere_connected_cen_fragment_mismatch(
    contig_data,
    telo_connected_graph_dict,
    cen_fragment_meta,
    stage_name,
):
    filtered_graph_dict = defaultdict(list)
    filtered_dict = {}
    filtered_set = set()
    removed = 0

    for telo_name, edge_list in telo_connected_graph_dict.items():
        for edge in edge_list:
            node_idx = edge[1]
            contig = contig_data[node_idx]
            chrom = contig[CHR_NAM]
            mismatch = False
            if node_is_censat(contig) and chrom in cen_fragment_meta and telo_name[-1] in ('f', 'b'):
                norm_dir = _normalized_telo_censat_dir(telo_name, contig[CTG_DIR])
                mismatch = norm_dir != _cen_fragment_target_dir_from_meta(cen_fragment_meta, chrom)

            if mismatch:
                removed += 1
                continue

            filtered_graph_dict[telo_name].append(edge)
            filtered_dict[node_idx] = telo_name
            filtered_set.add(node_idx)

    logging.info(
        f"Removed {removed} {stage_name} telomere-connected censat nodes "
        f"with cen_fragment direction mismatch"
    )
    return filtered_set, filtered_dict, filtered_graph_dict

def calc_ratio(contig_data : list) -> dict:
    contig_data_size = len(contig_data)
    contig_data.append((0, 0, 0, 0, 0, 0, 0, 0, 0))
    total_ref_len = 0
    curr_contig_first_fragment = contig_data[0]
    ref_qry_ratio = {}
    for i in range(1, contig_data_size+1):
        total_ref_len += contig_data[i-1][CHR_END] - contig_data[i-1][CHR_STR]
        if contig_data[i][CTG_NAM] != contig_data[i-1][CTG_NAM]:
            curr_contig_name = contig_data[i-1][CTG_NAM]
            curr_contig_end_fragment = contig_data[i-1]
            if curr_contig_first_fragment[CTG_DIR] == '+':
                estimated_ref_len = curr_contig_end_fragment[CHR_END] - curr_contig_first_fragment[CHR_STR]
            else:
                estimated_ref_len = curr_contig_first_fragment[CHR_END] - curr_contig_end_fragment[CHR_STR]
            try:
                ref_qry_ratio[curr_contig_name] = estimated_ref_len / total_ref_len
            except:
                if total_ref_len >=0:
                    ref_qry_ratio[curr_contig_name] = INF
                else:
                    ref_qry_ratio[curr_contig_name] = -INF
            total_ref_len = 0
            curr_contig_first_fragment = contig_data[i]
    contig_data = contig_data[:-1]
    return ref_qry_ratio

def calc_chukji(contig_data : list) -> tuple:
    contig_data_size = len(contig_data)
    contig_data.append((0, 0, 0, 0, 0, 0, 0, 0, 0))
    total_ref_len = 0
    curr_contig_first_fragment = contig_data[0]
    ref_qry_ratio = {}
    chukji = {}
    for i in range(1, contig_data_size+1):
        total_ref_len += contig_data[i-1][CHR_END] - contig_data[i-1][CHR_STR]
        if contig_data[i][CTG_NAM] != contig_data[i-1][CTG_NAM]:
            curr_contig_name = contig_data[i-1][CTG_NAM]
            curr_contig_end_fragment = contig_data[i-1]
            if curr_contig_first_fragment[CTG_DIR] == '+':
                estimated_ref_len = curr_contig_end_fragment[CHR_END] - curr_contig_first_fragment[CHR_STR]
            else:
                estimated_ref_len = curr_contig_first_fragment[CHR_END] - curr_contig_end_fragment[CHR_STR]
            try:
                ref_qry_ratio[curr_contig_name] = estimated_ref_len / total_ref_len
                chukji[curr_contig_name] = abs(estimated_ref_len)
            except:
                if total_ref_len >=0:
                    ref_qry_ratio[curr_contig_name] = INF
                else:
                    ref_qry_ratio[curr_contig_name] = -INF
            total_ref_len = 0
            curr_contig_first_fragment = contig_data[i]
    contig_data = contig_data[:-1]
    return ref_qry_ratio, chukji

def find_mainflow(contig_data : list) -> dict:
    contig_data_size = len(contig_data)
    mainflow_dict = {}
    st = 0
    while st<contig_data_size:
        ed = contig_data[st][CTG_ENDND]
        ref_length_counter = Counter()
        for i in range(st, ed+1):
            ref_length_counter[(contig_data[i][CTG_DIR], contig_data[i][CHR_NAM])]\
            +=contig_data[i][CTG_END]-contig_data[i][CTG_STR]
        max_count = 0
        max_chr = (0, 0)
        total_len = 0
        for i in ref_length_counter:
            total_len += ref_length_counter[i]
            if ref_length_counter[i]>max_count:
                max_count = ref_length_counter[i]
                max_chr = i
        mainflow_dict[contig_data[st][CTG_NAM]] = max_chr
        st = ed+1
    return mainflow_dict

def pipeline_preprocess_contig(contig_data : list, telo_label : list, ref_qry_ratio : dict, repeat_label : list, telo_connect_info : set) -> list :
    checker = 0
    contig_data_size = len(contig_data)-1
    curr_contig_first_fragment = contig_data[0]
    using_contig_list = []
    contig_terminal_node = {}
    contig_type = {}
    is_telo = False
    is_front_back_repeat = False
    chrM_flag = False
    idx = 0
    cnt = 0
    len_count = Counter()
    for i in range(1, contig_data_size+1):
        len_count[contig_data[i-1][CTG_NAM]]+=contig_data[i-1][CHR_END]-contig_data[i-1][CHR_STR]
        cnt+=1
        if (i-1) in telo_connect_info:
            is_telo = True
        if contig_data[i-1][CHR_NAM]=='chrM':
            chrM_flag = True
        # contig 넘어갈 때:
        if contig_data[i][CTG_NAM] != contig_data[i-1][CTG_NAM]:
            curr_contig_name = contig_data[i-1][CTG_NAM]
            if repeat_label[i-1][0]!='0':
                is_front_back_repeat = True
            curr_contig_end_fragment = contig_data[i-1]
            if curr_contig_first_fragment[CHR_NAM] != curr_contig_end_fragment[CHR_NAM]:
                checker = 1
            elif curr_contig_first_fragment[CTG_DIR] != curr_contig_end_fragment[CTG_DIR]:
                checker = 2
            else:
                if is_front_back_repeat:
                    bound = RPT_BND_CONTIG_BOUND
                else:
                    bound = BND_CONTIG_BOUND
                if abs(ref_qry_ratio[curr_contig_name]-1)< bound:
                    if contig_data[i-1][CTG_LEN] > CONTIG_MINIMUM_SIZE or is_telo:
                        checker = 3
                    elif cnt==1:
                        checker = 3
                else:
                    checker = 4
            if chrM_flag:
                checker = 0
            if checker>0:
                using_contig_list.append(curr_contig_name)
                contig_type[curr_contig_name] = checker
                contig_terminal_node[curr_contig_name] = (idx, idx+cnt-1)
                idx+=cnt
            # initialize
            curr_contig_first_fragment = contig_data[i]
            cnt = 0
            checker = 0
            is_telo = False
            is_front_back_repeat = False
            chrM_flag = False
            if i < contig_data_size and repeat_label[i][0]!='0':
                is_front_back_repeat = True
    contig_data = contig_data[:-1]
    return [using_contig_list, contig_type, contig_terminal_node, len_count]


def preprocess_contig(contig_data : list, telo_label : list, ref_qry_ratio : dict, repeat_label : list, telo_connect_info : set) -> list :
    checker = 0
    contig_data_size = len(contig_data)-1
    curr_contig_first_fragment = contig_data[0]
    using_contig_list = []
    contig_terminal_node = {}
    contig_type = {}
    is_telo = False
    is_front_back_repeat = False
    chrM_flag = False
    first_idx = 0
    idx = 0
    cnt = 0
    len_count = Counter()
    for i in range(1, contig_data_size+1):
        len_count[contig_data[i-1][CTG_NAM]]+=contig_data[i-1][CHR_END]-contig_data[i-1][CHR_STR]
        cnt+=1
        if (i-1) in telo_connect_info:
            is_telo = True
        if contig_data[i-1][CHR_NAM]=='chrM':
            chrM_flag = True
        # contig 넘어갈 때:
        if contig_data[i][CTG_NAM] != contig_data[i-1][CTG_NAM]:
            splited_type4 = False
            curr_contig_name = contig_data[i-1][CTG_NAM]
            if repeat_label[i-1][0]!='0':
                is_front_back_repeat = True
            curr_contig_end_fragment = contig_data[i-1]
            if curr_contig_first_fragment[CHR_NAM] != curr_contig_end_fragment[CHR_NAM]:
                checker = 1
            elif curr_contig_first_fragment[CTG_DIR] != curr_contig_end_fragment[CTG_DIR]:
                checker = 2
            else:
                if is_front_back_repeat:
                    bound = RPT_BND_CONTIG_BOUND
                else:
                    bound = BND_CONTIG_BOUND
                if abs(ref_qry_ratio[curr_contig_name]-1)< bound:
                    if contig_data[i-1][CTG_LEN] > CONTIG_MINIMUM_SIZE or is_telo:
                        checker = 3
                    elif cnt==1:
                        checker = 3
                else:
                    checker = 4
            if chrM_flag:
                checker = 0
            if checker==3:
                using_contig_list.append(curr_contig_name)
                contig_type[curr_contig_name] = checker
                contig_terminal_node[curr_contig_name] = (idx, idx+cnt-1)
                idx+=cnt
            # initialize
            curr_contig_first_fragment = contig_data[i]
            first_idx = i
            cnt = 0
            checker = 0
            is_telo = False
            is_front_back_repeat = False
            chrM_flag = False
            if i < contig_data_size and repeat_label[i][0]!='0':
                is_front_back_repeat = True
    contig_data = contig_data[:-1]
    return [using_contig_list, contig_type, contig_terminal_node, len_count]

def similar_check(v1, v2, ratio=TYPE2_SIM_COMPARE_RAITO):
    assert(v1 >= 0 and v2 >= 0)
    mi, ma = sorted([v1, v2])
    return False if mi == 0 else (ma / mi <= ratio) or ma-mi < NCLOSE_SIM_DIFF_THRESHOLD

def exist_near_bnd_point(depth_df, chrom, inside_st):
    # subset of depth_df for the given chromosome
    df_chr = depth_df[depth_df['chr'] == chrom]

    def mean_depth(start, end):
        """Return mean meandepth over windows overlapping [start, end)."""
        mask = (df_chr['nd'] > start) & (df_chr['st'] < end)
        return df_chr.loc[mask, 'meandepth'].mean()

    # for inside_st
    st_depth = mean_depth(inside_st - TYPE2_FLANKING_LENGTH, inside_st)
    nd_depth = mean_depth(inside_st, inside_st + TYPE2_FLANKING_LENGTH)

    # print(chrom, inside_st, inside_nd, not similar_check(st_depth, nd_depth))
    return not similar_check(st_depth, nd_depth)

def censat_overlap_check(censat_dict, chrom, inside_st, inside_nd):
    if chrom not in censat_dict.keys():
        return False
    for st, nd in censat_dict[chrom]:
        if not (nd < inside_st or inside_nd < st):
            return True
    return False

def alt_preprocess_contig(contig_data : list, telo_label : list, ref_qry_ratio : dict, repeat_label : list, telo_connect_info : set, telo_dict : dict, censat_dict : dict) -> list :
    checker = 0
    contig_data_size = len(contig_data)-1
    curr_contig_first_fragment = contig_data[0]
    using_contig_list = []
    contig_terminal_node = {}
    contig_type = {}
    is_telo = False
    is_front_back_repeat = False
    chrM_flag = False
    len_count = Counter()
    idx = 0
    cnt = 0
    telo_node_count = 0
    all_mapq_0 = True

    telo_inside_dict = dict()
    for chr_name, data in telo_dict.items():
        for s, e, d in data:
            key = chr_name + d
            func = max if d == 'f' else min
            if key not in telo_inside_dict:
                telo_inside_dict[key] = (s, e)

            telo_inside_dict[key] = func(telo_inside_dict[key], (s, e), key=lambda t: t[0])

    telo_safe_dict = defaultdict(list)
    for i in range(1, contig_data_size+1):
        len_count[contig_data[i-1][CTG_NAM]]+=contig_data[i-1][CHR_END]-contig_data[i-1][CHR_STR]
        cnt+=1
        if (i-1) in telo_connect_info:
            telo_node_count += 1
            is_telo = True
        if contig_data[i-1][CHR_NAM] == 'chrM':
            chrM_flag = True
        if contig_data[i-1][CTG_MAPQ] != 0:
            all_mapq_0 = False
        # contig 넘어갈 때:
        if contig_data[i][CTG_NAM] != contig_data[i-1][CTG_NAM]:
            if repeat_label[i-1][0]!='0':
                is_front_back_repeat = True
            curr_contig_name = contig_data[i-1][CTG_NAM]
            curr_contig_end_fragment = contig_data[i-1]
            if curr_contig_first_fragment[CHR_NAM] != curr_contig_end_fragment[CHR_NAM]:
                checker = 1
            elif curr_contig_first_fragment[CTG_DIR] != curr_contig_end_fragment[CTG_DIR]:
                checker = 2
            else:
                if is_telo:
                    checker = 5
                else:
                    if is_front_back_repeat:
                        bound = RPT_BND_CONTIG_BOUND
                    else:
                        bound = BND_CONTIG_BOUND

                    if abs(ref_qry_ratio[curr_contig_name]-1) >= bound:
                        checker = 4
                    else:
                        checker = 3

            if chrM_flag:
                checker = 0
            elif checker == 2:
                using_contig_list.append(curr_contig_name)
                contig_type[curr_contig_name] = checker
                contig_terminal_node[curr_contig_name] = (idx, idx+cnt-1)
                idx+=cnt

            elif (checker>0 and checker != 3 and checker != 5):
                using_contig_list.append(curr_contig_name)
                contig_type[curr_contig_name] = checker
                contig_terminal_node[curr_contig_name] = (idx, idx+cnt-1)
                idx+=cnt
            elif is_telo:
                using_contig_list.append(curr_contig_name)
                contig_type[curr_contig_name] = checker
                contig_terminal_node[curr_contig_name] = (idx, idx+telo_node_count-1)
                idx+=telo_node_count

            if checker == 3 and cnt == 1:
                now_chr = contig_data[i-1][CHR_NAM]
                st_nd = telo_inside_dict[now_chr + 'f'][1]
                nd_st = telo_inside_dict[now_chr + 'b'][0]

                st, nd = contig_data[i-1][CHR_STR], contig_data[i-1][CHR_END]
                if st_nd < st and nd < nd_st:
                    s, n = telo_inside_dict[now_chr + 'f']
                    telo_safe_dict[now_chr + 'f'].append((telo_distance_checker_cord(st, nd, s, n), curr_contig_name, i-1))

                    s, n = telo_inside_dict[now_chr + 'b']
                    telo_safe_dict[now_chr + 'b'].append((telo_distance_checker_cord(st, nd, s, n), curr_contig_name, i-1))


            curr_contig_first_fragment = contig_data[i]
            all_mapq_0 = True
            cnt = 0
            checker = 0
            telo_node_count = 0
            is_telo = False
            is_front_back_repeat = False
            chrM_flag = False
            if i < contig_data_size and repeat_label[i][0]!='0':
                is_front_back_repeat = True
    contig_data = contig_data[:-1]

    using_type3_contig_data = set()
    for k, node_list in telo_safe_dict.items():
        dis, curr_contig_name, ind = min(node_list, key=lambda t: t[0])
        using_type3_contig_data.add((curr_contig_name, ind))
    using_type3_contig_data = list(using_type3_contig_data)
    using_type3_contig_list = []
    using_type3_contig_data.sort(key=lambda t: t[1])
    for curr_contig_name, _ in using_type3_contig_data:
        if curr_contig_name not in using_contig_list:
            using_type3_contig_list.append(curr_contig_name)

            cnt = 1
            contig_type[curr_contig_name] = 3
            contig_terminal_node[curr_contig_name] = (idx, idx+cnt-1)
            idx+=cnt

    return [using_contig_list, using_type3_contig_list, contig_type, contig_terminal_node, len_count]

def preprocess_repeat(contig_data : list) -> list:
    repeat_preprocessed_contig = []
    contig_data_size = len(contig_data)
    curr_contig_st = 0
    idx = 0
    while curr_contig_st<contig_data_size:
        curr_contig_ed = contig_data[curr_contig_st][CTG_ENDND]
        front_telo_bound = curr_contig_st
        end_telo_bound = curr_contig_ed
        st = front_telo_bound
        ed = end_telo_bound
        if contig_data[curr_contig_st][CTG_TYP] != 3:
            while front_telo_bound <= curr_contig_ed and contig_data[front_telo_bound][CTG_TELCON] != '0':
                front_telo_bound+=1
            while end_telo_bound >= curr_contig_st and contig_data[end_telo_bound][CTG_TELCON] != '0':
                end_telo_bound-=1
            front_repeat_bound = front_telo_bound
            end_repeat_bound = end_telo_bound
            while front_repeat_bound<=curr_contig_ed \
                  and contig_data[front_repeat_bound][CTG_RPTCHR] != '0' \
                  and contig_data[front_repeat_bound][CTG_MAPQ] < MAPQ_BOUND:
                front_repeat_bound+=1
            if front_repeat_bound <= curr_contig_ed and front_repeat_bound > front_telo_bound:
                if contig_data[front_repeat_bound][CTG_RPTCHR] == '0':
                    front_repeat_bound-=1
            while end_repeat_bound>=curr_contig_st \
                  and contig_data[end_repeat_bound][CTG_RPTCHR] != '0' \
                  and contig_data[end_repeat_bound][CTG_MAPQ] < MAPQ_BOUND:
                end_repeat_bound-=1
            if end_repeat_bound >= curr_contig_st and end_repeat_bound < end_telo_bound:
                if contig_data[end_repeat_bound][CTG_RPTCHR] == '0':
                    end_repeat_bound+=1
            if front_repeat_bound >= front_telo_bound:
                st = front_repeat_bound
            if end_repeat_bound <= end_telo_bound:
                ed = end_repeat_bound
        repeat_temp_set = set()
        for i in range(curr_contig_st, front_telo_bound):
            repeat_temp_set.add(i)
        for i in range(st, ed+1):
            repeat_temp_set.add(i)
        for i in range(end_telo_bound+1, curr_contig_ed+1):
            repeat_temp_set.add(i)
        for i in sorted(list(repeat_temp_set)):
            repeat_preprocessed_contig.append(contig_data[i])
        for i in range(idx, len(repeat_preprocessed_contig)):
            repeat_preprocessed_contig[i][CTG_STRND] = idx
            repeat_preprocessed_contig[i][CTG_ENDND] = len(repeat_preprocessed_contig)-1

        idx = len(repeat_preprocessed_contig)
        curr_contig_st = curr_contig_ed+1
    return repeat_preprocessed_contig

def weighted_avg_meandepth(chrom_df, region_start, region_end):
    overlapping = chrom_df[(chrom_df['nd'] >= region_start) & (chrom_df['st'] <= region_end)]
    if overlapping.empty:
        return None
    total_weight = 0
    weighted_sum = 0
    for _, row in overlapping.iterrows():
        overlap_start = max(row['st'], region_start)
        overlap_end = min(row['nd'], region_end)
        if overlap_start <= overlap_end:
            length = overlap_end - overlap_start + 1
            weighted_sum += row['meandepth'] * length
            total_weight += length
    return weighted_sum / total_weight if total_weight > 0 else None


def find_breakend_centromere(
    repeat_censat_data : dict,
    chr_len : dict,
    df : pd.DataFrame,
    raw_nclose_nodes : dict = None,
    contig_data : list = None,
    log_context : str = "",
):
    # Remove unused Y
    ydf = df.query('chr == "chrY"')
    ydepth = np.mean(ydf['meandepth'].to_numpy())

    depth = np.mean(df['meandepth'].to_numpy())

    sd, bd = sorted([depth, ydepth])
    yratio = bd / sd
    if yratio > 2:
        df = df.query('chr != "chrY"')

    meandepth = np.median(df['meandepth'])

    def _node_overlaps_region(node, region_start, region_end):
        return max(int(node[CHR_STR]), region_start) <= min(int(node[CHR_END]), region_end)

    def _normalized_pair_endpoint_dir(pair, endpoint_order):
        ctg_dir = contig_data[pair[endpoint_order]][CTG_DIR]
        if endpoint_order == 0:
            return ctg_dir
        return _flip_ctg_dir(ctg_dir)

    def _expected_intervention_dir(side, right_high):
        # For a right-high jump, chr1:31M+ as pair[1] and chr1:33M- as pair[0]
        # normalize to '-' on the right flank, while chr1:13M+ normalizes to '+'
        # on the left flank and can explain the left-side drop. Reverse the roles
        # for a left-high jump.
        if right_high:
            return '+' if side == 'left' else '-'
        return '-' if side == 'left' else '+'

    def _has_depth_compatible_raw_nclose(chrom, rep_start_0, rep_end_0, flank_bp, right_high):
        if raw_nclose_nodes is None or contig_data is None:
            return False

        left_start = max(1, (rep_start_0 + 1) - flank_bp)
        left_end = rep_start_0
        right_start = rep_end_0 + 1
        right_end = min(chr_len[chrom], rep_end_0 + flank_bp)

        for pair_list in raw_nclose_nodes.values():
            for pair in pair_list:
                for endpoint_order, node_idx in enumerate(pair):
                    node = contig_data[node_idx]
                    if node[CHR_NAM] != chrom:
                        continue

                    side = None
                    if left_end >= left_start and _node_overlaps_region(node, left_start, left_end):
                        side = 'left'
                    elif right_end >= right_start and _node_overlaps_region(node, right_start, right_end):
                        side = 'right'
                    if side is None:
                        continue

                    norm_dir = _normalized_pair_endpoint_dir(pair, endpoint_order)
                    if norm_dir == _expected_intervention_dir(side, right_high):
                        return True
        return False

    depth_diff_data = dict()
    depth_dir_data = dict()
    relaxed_depth_chroms = set()
    for chrom, intervals in repeat_censat_data.items():
        chrom_length = chr_len.get(chrom)
        if chrom_length is None:
            continue
        chrom_df = df[df['chr'] == chrom]
        if chrom_df.empty:
            continue

        rep_start_0, rep_end_0 = intervals[0]  # 0-indexed 좌표

        if rep_start_0 == 0 or rep_end_0 == chrom_length - 1:
            continue

        tmp_depth_data = []
        for FLANK_SIZE_BP in [1*M, 5*M, 10*M]:
            # 좌측 flanking: repeat의 1-indexed 시작은 rep_start_0 + 1
            if rep_start_0 > 0:
                left_flank_end = rep_start_0  # repeat 시작 전 마지막 base (1-indexed)
                left_flank_start = max(1, (rep_start_0 + 1) - FLANK_SIZE_BP)

                if left_flank_end < left_flank_start or (left_flank_end - left_flank_start + 1) < FLANK_SIZE_BP:
                    left_flank_start = None
                    left_flank_end = None
                    left_weighted = None
                else:
                    left_weighted = weighted_avg_meandepth(chrom_df, left_flank_start, left_flank_start + FLANK_SIZE_BP)
            else:
                left_flank_start = None
                left_flank_end = None
                left_weighted = None

            # 우측 flanking: repeat의 끝 이후 첫 base부터
            if rep_end_0 < chrom_length:
                right_flank_start = rep_end_0 + 1
                right_flank_end = min(chrom_length, rep_end_0 + FLANK_SIZE_BP)

                if right_flank_end < right_flank_start or (right_flank_end - right_flank_start + 1) < FLANK_SIZE_BP:
                    right_flank_start = None
                    right_flank_end = None
                    right_weighted = None
                else:
                    right_weighted = weighted_avg_meandepth(chrom_df, right_flank_end - FLANK_SIZE_BP, right_flank_end)
            else:
                right_flank_start = None
                right_flank_end = None
                right_weighted = None

            if not (right_weighted is None or left_weighted is None):
                tmp_depth_data.append(right_weighted - left_weighted)
            else:
                break

        if not tmp_depth_data or tmp_depth_data[0] == 0:
            continue

        first_sign = tmp_depth_data[0] > 0

        rules = [(1, 0.20), (2, 0.15), (3, 0.10)]

        for k, thr in rules:
            if len(tmp_depth_data) < k:
                continue

            head = tmp_depth_data[:k]
            if any(x == 0 for x in head):
                continue

            same_sign = all((x > 0) == first_sign for x in head)
            strong_enough = all(abs(x) >= abs(thr * meandepth) for x in head)
            any_strong = any(abs(x) >= abs(thr * meandepth) for x in head)
            localized = all(
                abs(x - tmp_depth_data[0]) <= abs(CENSAT_OUT_DIFF_RATIO * meandepth)
                for x in head[1:]
            )

            if same_sign and strong_enough and localized:
                depth_diff_data[chrom] = abs(tmp_depth_data[0])
                depth_dir_data[chrom] = first_sign
                break
            elif raw_nclose_nodes is not None and contig_data is not None and k > 1 and same_sign and any_strong and localized:
                # A compatible breakpoint outside the inner 5 Mb can still change
                # the 5 Mb average by altering the segment toward the censat edge.
                nclose_flank_bp = 10*M
                has_nclose_intervention = _has_depth_compatible_raw_nclose(
                    chrom,
                    rep_start_0,
                    rep_end_0,
                    nclose_flank_bp,
                    first_sign,
                )
                if not has_nclose_intervention:
                    depth_diff_data[chrom] = max(
                        abs(x) for x in head if abs(x) >= abs(thr * meandepth)
                    )
                    depth_dir_data[chrom] = first_sign
                    relaxed_depth_chroms.add(chrom)
                    break


    cen_fragment_meta = {}
    for chrom, diff in depth_diff_data.items():
        intervals = repeat_censat_data[chrom]
        rep_start_0, rep_end_0 = intervals[0]
        mid_censat = (rep_start_0 + rep_end_0) // 2

        cen_fragment_meta[chrom] = {
            "dir": depth_dir_data[chrom],
            "mid": mid_censat,
            "depth_diff": diff,
            "chr_len": chr_len[chrom],
        }

    # censat(centromere) breakend 로 인식된 염색체 목록을 로그로 출력
    detected_summary = ', '.join(sorted(cen_fragment_meta, key=chr2int))
    log_prefix = f'{log_context} ' if log_context else ''
    logging.info(
        f'{log_prefix}Censat breakend chromosome: '
        f'{detected_summary if cen_fragment_meta else "(none)"}'
    )
    if relaxed_depth_chroms:
        relaxed_summary = ', '.join(sorted(relaxed_depth_chroms, key=chr2int))
        logging.info(
            f'{log_prefix}Censat breakend chromosome relaxed by raw nclose direction: '
            f'{relaxed_summary}'
        )

    return cen_fragment_meta

def break_double_telomere_contig(contig_data : list, telo_connected_set : set):
    s = 0
    vtg_list = []
    idx = 0
    count = 0
    tot_cur_ctg_cnt = 0
    contig_data_size = len(contig_data)
    while s<contig_data_size:
        e = contig_data[s][CTG_ENDND]
        cnt = 0
        for i in range(s+1, e):
            if contig_data[i][CTG_TELDIR][0] in ('f', 'b')\
            and s not in telo_connected_set \
            and e not in telo_connected_set:
                cnt+=1
        if cnt >= 2:
            st = s
            ed = e
            while st <= e and contig_data[st][CTG_TELDIR] == '0':
                st+=1
            while ed >= s and contig_data[ed][CTG_TELDIR] == '0':
                ed-=1
            if s <= st:
                tot_cur_ctg_cnt+=1
            for i in range(s, st+1):
                temp_list = copy.deepcopy(contig_data[i])
                temp_list[CTG_NAM] = f'telomere_middle_cut_contig_{tot_cur_ctg_cnt}'
                temp_list[CTG_STRND] = idx
                temp_list[CTG_ENDND] = idx + st - s
                vtg_list.append(temp_list)
                count+=1
            idx += st-s+1
            if ed<=e:
                tot_cur_ctg_cnt+=1
            for i in range(ed, e+1):
                temp_list = copy.deepcopy(contig_data[i])
                temp_list[CTG_NAM] = f'telomere_middle_cut_contig_{tot_cur_ctg_cnt}'
                temp_list[CTG_STRND] = idx
                temp_list[CTG_ENDND] = idx + e - ed
                vtg_list.append(temp_list)
                count+=1
        s = e+1
    return vtg_list

def delete_contig(contig_data : list, to_delete_contig_set : set) -> list:
    new_contig_data = []
    delete_count = 0
    for chunk in contig_data:
        if chunk[CTG_NAM] in to_delete_contig_set:
            delete_count += 1
        else:
            chunk[CTG_STRND] -= delete_count
            chunk[CTG_ENDND] -= delete_count
            new_contig_data.append(chunk)

    return new_contig_data

def break_type34_contig(contig_data : list):
    s = 0
    vtg_list = []
    broken_contig_set = set()

    while s < len(contig_data):
        e = contig_data[s][CTG_ENDND]
        curr_contig_name = contig_data[s][CTG_NAM]
        broken_chk = False  # reset per contig: avoid UnboundLocalError on non-type3 first contig and stale carry-over
        if contig_data[s][CTG_TYP] == 3:
            cnt = 0
            for i in range(s, e):
                front_chr = (contig_data[i][CHR_NAM], contig_data[i][CTG_DIR])
                back_chr = (contig_data[i+1][CHR_NAM], contig_data[i+1][CTG_DIR])
                ref_ratio, chukji = calculate_single_contig_ref_ratio(contig_data[i:i+2])
                if contig_data[i][CTG_RPTCHR] != '0' or contig_data[i+1][CTG_RPTCHR] != '0':
                    bnd_bound = RPT_BND_CONTIG_BOUND
                else:
                    bnd_bound = BND_CONTIG_BOUND
                if front_chr == back_chr and abs(1-ref_ratio) > bnd_bound and abs(chukji) > TYPE34_BREAK_CHUKJI_LIMIT:
                    broken_chk = True
                    cnt += 1
                    vtg = copy.deepcopy(contig_data[i:i+2])
                    for j in vtg:
                        j[CTG_NAM] = f'{curr_contig_name}_split34_{cnt}'
                        vtg_list.append(j)
        if broken_chk:
            broken_contig_set.add(curr_contig_name)
        s = e+1
    return vtg_list, broken_contig_set

def pass_pipeline(pre_contig_data, telo_dict, telo_bound_dict, repeat_data, repeat_censat_data, telo_ppc_passed, chr_len):
    if not telo_ppc_passed:
        if len(pre_contig_data)==0:
            return []
        contig_data = []

        for i in pre_contig_data:
            contig_data.append(i[:9] + [i[CTG_MAPQ], i[CTG_GLOBALIDX]])

        node_label = label_node(contig_data, telo_dict)

        repeat_label = label_repeat_node(contig_data, repeat_data, chr_len)

        telo_preprocessed_contig, report_case, telo_connect_info = preprocess_telo(contig_data, node_label)

        new_contig_data = []
        telcon_set = set()
        idxcnt = 0
        for i in telo_preprocessed_contig:
            temp_list = contig_data[i]
            if i in telo_connect_info:
                temp_list.append(telo_connect_info[i])
                telcon_set.add(idxcnt)
            else:
                temp_list.append("0")
            new_contig_data.append(temp_list)
            idxcnt+=1


        new_node_repeat_label = label_repeat_node(new_contig_data, repeat_data, chr_len)
        new_node_repeat_censat_label = label_repeat_node(new_contig_data, repeat_censat_data, chr_len)
        new_node_telo_label = label_node(new_contig_data, telo_dict)

        ref_qry_ratio = calc_ratio(new_contig_data)
        preprocess_result, \
        preprocess_contig_type, \
        preprocess_terminal_nodes, \
        len_counter = pipeline_preprocess_contig(new_contig_data, new_node_telo_label, ref_qry_ratio, new_node_repeat_label, telcon_set)
        preprocess_result = set(preprocess_result)
        new_contig_data = new_contig_data[0:-1]
        contig_data_size = len(new_contig_data)

        telo_ppc_contig = []
        cnt = 0
        for i in range(contig_data_size):
            if new_contig_data[i][CTG_NAM] in preprocess_result:
                temp_list = new_contig_data[i][:10]
                temp_list.append(preprocess_contig_type[new_contig_data[i][CTG_NAM]])
                temp_list.append(preprocess_terminal_nodes[new_contig_data[i][CTG_NAM]][0])
                temp_list.append(preprocess_terminal_nodes[new_contig_data[i][CTG_NAM]][1])
                temp_list.append(new_node_telo_label[i][0])
                temp_list.append(new_node_telo_label[i][1])
                temp_list.append(new_contig_data[i][11])
                temp_list.append(new_node_repeat_label[i][0])
                temp_list.append(new_node_repeat_label[i][1])
                temp_list.append(new_node_repeat_censat_label[i][1])
                temp_list.append(new_contig_data[i][10])
                telo_ppc_contig.append(temp_list)
                cnt+=1
        mainflow_dict = find_mainflow(telo_ppc_contig)
        telo_ppc_size = len(telo_ppc_contig)

        for i in range(telo_ppc_size):
            max_chr = mainflow_dict[telo_ppc_contig[i][CTG_NAM]]
            temp = [max_chr[0], max_chr[1], telo_ppc_contig[i][-1]]
            telo_ppc_contig[i] = telo_ppc_contig[i][:-1]
            telo_ppc_contig[i] += temp
    else:
        if len(pre_contig_data)==0:
            return []
        telo_ppc_contig = pre_contig_data
        mainflow_dict = find_mainflow(telo_ppc_contig)
        telo_ppc_size = len(telo_ppc_contig)

        for i in range(telo_ppc_size):
            max_chr = mainflow_dict[telo_ppc_contig[i][CTG_NAM]]
            telo_ppc_contig[i][CTG_MAINFLOWDIR] = max_chr[0]
            telo_ppc_contig[i][CTG_MAINFLOWCHR] = max_chr[1]



    final_contig = preprocess_repeat(telo_ppc_contig)

    final_contig_repeat_label = label_repeat_node(final_contig, repeat_data, chr_len)

    final_telo_node_label = label_node(final_contig, telo_dict)

    final_ref_qry_ratio = calc_ratio(final_contig)

    final_telo_connect = set()

    for v, i in enumerate(final_contig):
        try:
            if i[CTG_TELCON] != '0':
                final_telo_connect.add(v)
        except:
            pass

    final_using_contig, final_ctg_typ, final_preprocess_terminal_nodes, _ = pipeline_preprocess_contig(final_contig, final_telo_node_label, final_ref_qry_ratio, final_contig_repeat_label, final_telo_connect)

    final_contig = final_contig[:-1]

    real_final_contig = []

    final_using_contig = set(final_using_contig)
    for i in range(0, len(final_contig)):
        if final_contig[i][CTG_NAM] in final_using_contig:
            final_contig[i][CTG_TYP] = final_ctg_typ[final_contig[i][CTG_NAM]]
            final_contig[i][CTG_STRND] = final_preprocess_terminal_nodes[final_contig[i][CTG_NAM]][0]
            final_contig[i][CTG_ENDND] = final_preprocess_terminal_nodes[final_contig[i][CTG_NAM]][1]
            real_final_contig.append(final_contig[i])

    return real_final_contig
















def make_virtual_ord_ctg(contig_data, fake_bnd):
    contig_data_size = len(contig_data)
    s = 0
    virtual_ctg = []
    idx = 0
    while s < contig_data_size:
        e = contig_data[s][CTG_ENDND]
        curr_ctg = contig_data[s][CTG_NAM]
        if curr_ctg in fake_bnd.keys():
            vctg_rng = fake_bnd[curr_ctg]
            for i in range(vctg_rng[0], vctg_rng[1] +1):
                temp_list = list(contig_data[i])
                temp_list[CTG_NAM] = temp_list[CTG_NAM][:-1] + "l"
                temp_list[CTG_TYP] = 3
                temp_list[CTG_STRND] = idx + contig_data_size
                temp_list[CTG_ENDND] = idx + vctg_rng[1] - vctg_rng[0] + contig_data_size
                virtual_ctg.append(temp_list)
            idx+=vctg_rng[1] - vctg_rng[0]+1
        s = e+1
    return virtual_ctg

def extract_telomere_connect_contig(telo_info_path : str) -> dict:
    telomere_connect_contig = defaultdict(list)
    with open(telo_info_path) as f:
        for curr_data in f:
            curr_data = curr_data.rstrip()
            temp_list = curr_data.split("\t")
            chr_info = temp_list[0]
            contig_id = ast.literal_eval(temp_list[1])
            telomere_connect_contig[chr_info].append(contig_id)

    return telomere_connect_contig

def extract_telomere_connect_contig_bytuple(telo_info_path : str) -> list:
    telomere_connect_contig = []
    with open(telo_info_path) as f:
        for curr_data in f:
            curr_data = curr_data.rstrip()
            temp_list = curr_data.split("\t")
            chr_info = temp_list[0]
            contig_id = ast.literal_eval(temp_list[1])
            telomere_connect_contig.append((chr_info, contig_id[1]))

    return telomere_connect_contig


def telomere_name_sort_key(telo_name):
    telo_name = str(telo_name)
    side = telo_name[-1] if telo_name.endswith(('f', 'b')) else ''
    chrom = telo_name[:-1] if side else telo_name
    try:
        chrom_order = chr2int(chrom)
    except (TypeError, ValueError):
        chrom_order = INF
    return ({'f': 0, 'b': 1}.get(side, 2), chrom_order, chrom, telo_name)


def write_telomere_connected_outputs(prefix: str, telo_edges: list,
                                      contig_data: list) -> None:
    telo_edges_by_name = defaultdict(list)
    with open(f"{prefix}/telomere_connected_list.txt", "wt") as f:
        for telo_name, edge in telo_edges:
            edge = tuple(edge)
            telo_edges_by_name[telo_name].append(edge)
            print(telo_name, edge, sep="\t", file=f)

    with open(f"{prefix}/telomere_connected_list_readable.txt", "wt") as f:
        for telo_name in sorted(telo_edges_by_name, key=telomere_name_sort_key):
            edges = telo_edges_by_name[telo_name]
            print(telo_name, file=f)
            for edge in edges:
                print(tuple(contig_data[edge[1]]), file=f)
            print("", file=f)


def group_nclose_nodes_by_chrom(contig_data: list, nclose_nodes: dict) -> dict:
    nclose_type = defaultdict(list)
    for pair_list in nclose_nodes.values():
        for pair in pair_list:
            contig_a = contig_data[pair[0]]
            contig_b = contig_data[pair[1]]
            if chr2int(contig_a[CHR_NAM]) <= chr2int(contig_b[CHR_NAM]):
                chrom_pair = (contig_a[CHR_NAM], contig_b[CHR_NAM])
                normalized_pair = tuple(pair)
            else:
                chrom_pair = (contig_b[CHR_NAM], contig_a[CHR_NAM])
                normalized_pair = (pair[1], pair[0])
            nclose_type[chrom_pair].append(normalized_pair)
    return nclose_type


def write_nclose_nodes_list(path: str, nclose_type: dict, contig_data: list,
                            repeat_contig_names=()) -> int:
    node_count = 0
    with open(path, "wt") as f:
        for chrom_pair, pair_list in nclose_type.items():
            print(f"{chrom_pair[0]}, {chrom_pair[1]}, {len(pair_list)}", file=f)
            for pair in pair_list:
                node_count += 2
                contig_a = contig_data[pair[0]]
                contig_b = contig_data[pair[1]]
                is_for = pair[0] < pair[1]
                list_a = [
                    contig_a[CTG_NAM],
                    get_corr_dir(is_for, contig_a[CTG_DIR]),
                    contig_a[CHR_STR],
                    contig_a[CHR_END],
                ]
                list_b = [
                    contig_b[CTG_NAM],
                    get_corr_dir(is_for, contig_b[CTG_DIR]),
                    contig_b[CHR_STR],
                    contig_b[CHR_END],
                ]
                if contig_a[CTG_NAM] in repeat_contig_names:
                    print(list_a, list_b, "all_repeat", file=f)
                else:
                    print(list_a, list_b, file=f)
            print("", file=f)
    return node_count


def write_nclose_nodes_index(path: str, nclose_nodes: dict,
                             contig_data: list) -> int:
    node_count = 0
    with open(path, "wt") as f:
        for contig_name, pair_list in nclose_nodes.items():
            for pair in pair_list:
                node_count += 2
                print(
                    contig_name,
                    pair[0],
                    pair[1],
                    contig_data[pair[0]][CTG_TYP],
                    file=f,
                )
    return node_count


def write_virtual_ordinary_contig(path: str, contig_data: list) -> None:
    with open(path, "wt") as f:
        for contig in contig_data:
            for field in contig:
                print(field, end="\t", file=f)
            print("", file=f)


def initialize_inversion_only_graph(contig_data : list, nclose_nodes_original : dict) -> dict:
    bnd_adjacency = defaultdict(list)
    nclose_nodes = defaultdict(list)

    for j in nclose_nodes_original:
        for i in nclose_nodes_original[j]:
            if contig_data[i[0]][CHR_NAM] == contig_data[i[1]][CHR_NAM]:
                nclose_nodes[j].append(i)

    for j in nclose_nodes:
        for i in nclose_nodes[j]:
            bnd_adjacency[(DIR_FOR, i[0])].append([DIR_FOR, i[1]])
            bnd_adjacency[(DIR_BAK, i[1])].append([DIR_BAK, i[0]])
            if len(i)==4:
                bnd_adjacency[(DIR_FOR, i[2])].append([DIR_FOR, i[3]])
                bnd_adjacency[(DIR_BAK, i[3])].append([DIR_BAK, i[2]])
    nclose_nodes_key = list(nclose_nodes.keys())
    key_len = len(nclose_nodes)
    # nclose nodes connection.
    for i1 in range(key_len):
        for i2 in range(i1+1, key_len):
            for n1 in nclose_nodes[nclose_nodes_key[i1]]:
                for n2 in nclose_nodes[nclose_nodes_key[i2]]:
                    n1_len = len(n1)
                    assert(n1_len==2)
                    n2_len = len(n2)
                    assert(n2_len==2)
                    for i in range(n1_len):
                        for j in range(n2_len):
                            if contig_data[n1[i]][CHR_NAM]==contig_data[n2[j]][CHR_NAM]:
                                contig_i = contig_data[n1[i]]
                                contig_j = contig_data[n2[j]]
                                i_ind = 'b' if i%2 else 'f'
                                j_ind = 'b' if j%2 else 'f'
                                i_ind += contig_i[CTG_DIR]
                                j_ind += contig_j[CTG_DIR]
                                i_str = contig_i[CHR_STR]
                                j_str = contig_j[CHR_STR]
                                i_end = contig_i[CHR_END]
                                j_end = contig_j[CHR_END]
                                is_overlap = distance_checker(contig_i, contig_j)==0
                                if i_ind == 'f+':
                                    if j_ind == 'f-':
                                        if i_str > j_str or is_overlap:
                                            bnd_adjacency[(DIR_BAK, n1[i])].append([DIR_FOR, n2[j]])
                                            bnd_adjacency[(DIR_BAK, n2[j])].append([DIR_FOR, n1[i]])
                                    elif j_ind == 'b+':
                                        if i_str > j_end or is_overlap:
                                            bnd_adjacency[(DIR_FOR, n2[j])].append([DIR_FOR, n1[i]])
                                            bnd_adjacency[(DIR_BAK, n1[i])].append([DIR_BAK, n2[j]])
                                elif i_ind == 'b+':
                                    if j_ind == 'b-':
                                        if j_end > i_end or is_overlap:
                                            bnd_adjacency[(DIR_FOR, n1[i])].append([DIR_BAK, n2[j]])
                                            bnd_adjacency[(DIR_FOR, n2[j])].append([DIR_BAK, n1[i]])
                                    elif j_ind == 'f+':
                                        if j_str > i_end or is_overlap:
                                            bnd_adjacency[(DIR_BAK, n2[j])].append([DIR_BAK, n1[i]])
                                            bnd_adjacency[(DIR_FOR, n1[i])].append([DIR_FOR, n2[j]])
                                elif i_ind == 'f-':
                                    if j_ind == 'f+':
                                        if j_str > i_str or is_overlap:
                                            bnd_adjacency[(DIR_BAK, n1[i])].append([DIR_FOR, n2[j]])
                                            bnd_adjacency[(DIR_BAK, n2[j])].append([DIR_FOR, n1[i]])
                                    elif j_ind == 'b-':
                                        if j_end > i_str or is_overlap:
                                            bnd_adjacency[(DIR_FOR, n2[j])].append([DIR_FOR, n1[i]])
                                            bnd_adjacency[(DIR_BAK, n1[i])].append([DIR_BAK, n2[j]])
                                elif i_ind == 'b-':
                                    if j_ind == 'b+':
                                        if i_end > j_end or is_overlap:
                                            bnd_adjacency[(DIR_FOR, n1[i])].append([DIR_BAK, n2[j]])
                                            bnd_adjacency[(DIR_FOR, n2[j])].append([DIR_BAK, n1[i]])
                                    elif j_ind == 'f-':
                                        if i_end > j_str or is_overlap:
                                            bnd_adjacency[(DIR_BAK, n2[j])].append([DIR_BAK, n1[i]])
                                            bnd_adjacency[(DIR_FOR, n1[i])].append([DIR_FOR, n2[j]])

    return bnd_adjacency


def convert_all_nclose_comp_to_nclose_nodes(contig_data : list, all_nclose_comp : dict) -> dict:
    all_nclose_nodes = defaultdict(list)
    seen_nclose = set()

    for pairs in all_nclose_comp.values():
        for pair in pairs:
            pair = tuple(pair)
            if len(pair) != 2 or pair in seen_nclose:
                continue
            all_nclose_nodes[contig_data[pair[0]][CTG_NAM]].append(pair)
            seen_nclose.add(pair)

    return all_nclose_nodes


def build_ecdna_nclose_nodes(raw_nclose_nodes : dict, all_nclose_nodes : dict) -> dict:
    ecdna_nclose_nodes = defaultdict(list)
    seen_nclose = set()

    for nclose_source in (raw_nclose_nodes, all_nclose_nodes):
        for contig_name, pairs in nclose_source.items():
            for pair in pairs:
                pair = tuple(pair)
                if len(pair) != 2 or pair in seen_nclose:
                    continue
                ecdna_nclose_nodes[contig_name].append(pair)
                seen_nclose.add(pair)

    return ecdna_nclose_nodes


def make_inversion_nx_graph(bnd_graph_adjacency, contig_rows):
    G = nx.DiGraph()
    for i in bnd_graph_adjacency:
        G.add_node(i)
    for node in bnd_graph_adjacency:
        for edge in bnd_graph_adjacency[node]:
            if type(edge)==str:
                G.add_weighted_edges_from([(tuple(node), edge, 0)])
            elif type(node)==str:
                G.add_weighted_edges_from([(node, tuple(edge), 0)])
            else:
                d = 0
                contig_s = contig_rows[node[1]]
                contig_e = contig_rows[edge[1]]
                if contig_s[CTG_NAM] == contig_e[CTG_NAM]:
                    if node[1] < edge[1]:
                        d = 0
                        for i in range(node[1]+1, edge[1]):
                            d += contig_rows[i][CHR_END] - contig_rows[i][CHR_STR]
                    else:
                        d = 0
                        for i in range(edge[1]+1, node[1]):
                            d += contig_rows[i][CHR_END] - contig_rows[i][CHR_STR]

                else:
                    ds = contig_s[CHR_END] - contig_s[CHR_STR]
                    de = contig_e[CHR_END] - contig_e[CHR_STR]
                    if distance_checker(contig_s, contig_e)==0:
                        d = ds + de - overlap_calculator(contig_s, contig_e)
                    else:
                        d = distance_checker(contig_s, contig_e) + ds + de
                G.add_weighted_edges_from([(tuple(node), tuple(edge), d)])
    return G

def circuit_length_calculator(circuit, contig_rows):
    circuit_len = 0
    for i in range(len(circuit)):
        curr_node = circuit[i]
        next_node = circuit[i+1] if i+1 < len(circuit) else circuit[0]
        circuit_len += abs(contig_rows[curr_node][CHR_END] - contig_rows[curr_node][CHR_STR])
        circuit_len += abs(contig_rows[curr_node][CHR_END] - contig_rows[next_node][CHR_STR])
    return circuit_len


def build_ecdna_circuits(contig_rows, raw_nclose_nodes, all_nclose_comp):
    """Reproduce the legacy four-node inversion-circuit ecDNA candidates."""

    all_nclose_nodes = convert_all_nclose_comp_to_nclose_nodes(
        contig_rows,
        all_nclose_comp,
    )
    ecdna_nclose_nodes = build_ecdna_nclose_nodes(
        raw_nclose_nodes,
        all_nclose_nodes,
    )
    logging.debug(
        "ecDNA nclose input: %d raw nclose, %d all-file nclose, %d merged nclose",
        sum(len(value) for value in raw_nclose_nodes.values()),
        sum(len(value) for value in all_nclose_nodes.values()),
        sum(len(value) for value in ecdna_nclose_nodes.values()),
    )

    adjacency = initialize_inversion_only_graph(
        contig_rows,
        ecdna_nclose_nodes,
    )
    nx_graph = make_inversion_nx_graph(adjacency, contig_rows)
    tool_graph = gt.Graph(directed=True)
    node_to_vertex = {}
    vertex_to_node = {}
    for node in nx_graph.nodes:
        vertex = tool_graph.add_vertex()
        node_to_vertex[node] = vertex
        vertex_to_node[int(vertex)] = node
    for source, target in nx_graph.edges:
        tool_graph.add_edge(node_to_vertex[source], node_to_vertex[target])

    circuit_candidates = set()
    for circuit in gt.all_circuits(tool_graph, max_length=4):
        if len(circuit) != 4:
            continue
        circuit_candidates.add(
            tuple(
                sorted(
                    int(vertex_to_node[int(vertex)][1])
                    for vertex in circuit
                )
            )
        )

    circuits = sorted(
        circuit
        for circuit in circuit_candidates
        if circuit_length_calculator(circuit, contig_rows)
        < CIRCUIT_ECDNA_LENGTH_LIMIT
    )
    return circuits, ecdna_nclose_nodes


def save_ecdna_outputs(prefix, contig_rows, raw_nclose_nodes, all_nclose_comp):
    circuits, ecdna_nclose_nodes = build_ecdna_circuits(
        contig_rows,
        raw_nclose_nodes,
        all_nclose_comp,
    )
    with open(os.path.join(prefix, "ecdna_circuit_data.pkl"), "wb") as handle:
        pkl.dump((circuits, ecdna_nclose_nodes), handle)
    return circuits, ecdna_nclose_nodes


def load_contig_preprocess_resources(context):
    """Load reference/depth inputs once for one preprocessing attempt."""

    chr_len = find_chr_len(context.reference_fai_path)
    telo_data = import_telo_data(context.telomere_bed_path, chr_len)
    repeat_data = import_repeat_data_00(context.repeat_bed_path)
    repeat_censat_data = import_censat_repeat_data(context.censat_bed_path)
    depth_df = pd.read_csv(
        context.main_stat_path,
        compression='gzip',
        comment='#',
        sep='\t',
        names=[
            'chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov',
            'meandepth',
        ],
    ).query('chr != "chrM"')
    cen_fragment_meta = find_breakend_centromere(
        repeat_censat_data,
        chr_len,
        depth_df,
        log_context="Strict",
    )
    with open(f'{context.prefix}/cen_fragment_data.pkl', 'wb') as f:
        pkl.dump(cen_fragment_meta, f)

    telo_dict = defaultdict(list)
    for row in telo_data:
        telo_dict[row[0]].append(row[1:])

    bed_data = import_bed(context.censat_bed_path)
    ydf = depth_df.query('chr == "chrY"')
    bed_intervals = pd.IntervalIndex.from_tuples(
        bed_data['chrY'],
        closed='left',
    )
    y_mask = ydf.apply(
        lambda row: bed_intervals.overlaps(
            pd.Interval(row['st'], row['nd'], closed='left')
        ),
        axis=1,
    )
    correct_mask = y_mask.apply(any)
    ydf_not_censat = ydf[~correct_mask]
    chry_nz_len = len(ydf_not_censat.query('meandepth != 0'))
    no_chrY = (chry_nz_len / len(ydf_not_censat)) < chrY_MINIMUM_RATIO

    return ContigPreprocessResources(
        chr_len=chr_len,
        telo_data=telo_data,
        telo_dict=telo_dict,
        repeat_data=repeat_data,
        repeat_censat_data=repeat_censat_data,
        cen_fragment_meta=cen_fragment_meta,
        depth_df=depth_df,
        no_chrY=no_chrY,
    )


def _prepare_paf_source_rows(resources, paf_path, policy):
    """Run the source-independent PAF/telomere labelling prefix."""

    contig_data = import_data(paf_path)
    original_node_count = len(contig_data)
    excluded_contigs, excluded_rows = find_multi_end_aligned_contigs(contig_data)
    if excluded_contigs:
        logging.info(
            f"Detected {len(excluded_contigs)} multi-end-aligned contigs "
            f"({len(excluded_rows)} PAF rows)"
        )

    node_label = label_node(contig_data, resources.telo_dict)
    telo_preprocessed_contig, _, telo_connect_info = preprocess_telo(
        contig_data,
        node_label,
    )
    excluded_telo_candidates = sum(
        row_idx in excluded_rows for row_idx in telo_connect_info
    )
    telo_connect_info = {
        row_idx: telo_name
        for row_idx, telo_name in telo_connect_info.items()
        if row_idx not in excluded_rows
    }
    if excluded_telo_candidates:
        logging.info(
            f"Excluded {excluded_telo_candidates} {policy.log_label} "
            "telomere boundary candidates from multi-end-aligned contigs"
        )

    new_contig_data = []
    telcon_set = set()
    for row_idx in telo_preprocessed_contig:
        row = contig_data[row_idx]
        if row_idx in telo_connect_info:
            row.append(telo_connect_info[row_idx])
            telcon_set.add(len(new_contig_data))
        else:
            row.append("0")
        new_contig_data.append(row)

    return new_contig_data, telcon_set, excluded_rows, original_node_count


def _append_source_row(
    output,
    row,
    contig_type,
    terminal_nodes,
    telo_label,
    repeat_label,
    censat_label,
    global_index_prefix,
):
    """Append one legacy 22-column PPC row without changing its schema."""

    temp_list = row[:10]
    temp_list.append(contig_type)
    temp_list.append(terminal_nodes[0])
    temp_list.append(terminal_nodes[1])
    temp_list.append(telo_label[0])
    temp_list.append(telo_label[1])
    temp_list.append(row[11])
    temp_list.append(repeat_label[0])
    temp_list.append(repeat_label[1])
    temp_list.append(censat_label[1])
    temp_list.append(f'{global_index_prefix}.{row[10]}')
    output.append(temp_list)


def _attach_mainflow(contigs):
    """Add the legacy main-flow fields while retaining row order."""

    mainflow_dict = find_mainflow(contigs)
    for row_idx in range(len(contigs)):
        max_chr = mainflow_dict[contigs[row_idx][CTG_NAM]]
        global_index = contigs[row_idx][-1]
        contigs[row_idx] = contigs[row_idx][:-1] + [
            max_chr[0],
            max_chr[1],
            global_index,
        ]
    return contigs


def preprocess_paf_source(
    context,
    resources,
    paf_path,
    policy,
    index_offset=0,
):
    """Preprocess one PAF source under an explicit primary/secondary policy."""

    (
        new_contig_data,
        telcon_set,
        excluded_rows,
        original_node_count,
    ) = _prepare_paf_source_rows(resources, paf_path, policy)
    new_node_telo_label = label_node(new_contig_data, resources.telo_dict)
    new_node_repeat_label = label_repeat_node(
        new_contig_data,
        resources.repeat_data,
        resources.chr_len,
    )
    new_node_repeat_censat_label = label_repeat_node(
        new_contig_data,
        resources.repeat_censat_data,
        resources.chr_len,
    )
    ref_qry_ratio = calc_ratio(new_contig_data)

    first_pass_rows = []
    if policy.kind == PafSourceKind.PRIMARY:
        (
            preprocess_result,
            preprocess_contig_type,
            preprocess_terminal_nodes,
            _,
        ) = preprocess_contig(
            new_contig_data,
            new_node_telo_label,
            ref_qry_ratio,
            new_node_repeat_label,
            telcon_set,
        )
        selected_groups = [set(preprocess_result)]
    else:
        (
            preprocess_result,
            preprocess_type3_result,
            preprocess_contig_type,
            preprocess_terminal_nodes,
            _,
        ) = alt_preprocess_contig(
            new_contig_data,
            new_node_telo_label,
            ref_qry_ratio,
            new_node_repeat_label,
            telcon_set,
            resources.telo_dict,
            resources.repeat_censat_data,
        )
        selected_groups = [set(preprocess_result), set(preprocess_type3_result)]

    new_contig_data = new_contig_data[:-1]
    for selected_names in selected_groups:
        for row_idx, row in enumerate(new_contig_data):
            if row[CTG_NAM] not in selected_names:
                continue
            contig_type = preprocess_contig_type[row[CTG_NAM]]
            if policy.kind == PafSourceKind.SECONDARY and contig_type == 5:
                if row_idx not in telcon_set:
                    continue
                contig_type = 3
            _append_source_row(
                first_pass_rows,
                row,
                contig_type,
                preprocess_terminal_nodes[row[CTG_NAM]],
                new_node_telo_label[row_idx],
                new_node_repeat_label[row_idx],
                new_node_repeat_censat_label[row_idx],
                policy.global_index_prefix,
            )

    first_pass_rows = _attach_mainflow(first_pass_rows)
    final_contig = preprocess_repeat(first_pass_rows)
    final_repeat_label = label_repeat_node(
        final_contig,
        resources.repeat_data,
        resources.chr_len,
    )
    final_telo_label = label_node(final_contig, resources.telo_dict)
    final_telo_connect = {
        row_idx
        for row_idx, row in enumerate(final_contig)
        if row[CTG_TELCON] != '0'
    }
    final_ref_qry_ratio = calc_ratio(final_contig)

    if policy.kind == PafSourceKind.PRIMARY:
        (
            final_using_contig,
            final_ctg_typ,
            final_terminal_nodes,
            _,
        ) = preprocess_contig(
            final_contig,
            final_telo_label,
            final_ref_qry_ratio,
            final_repeat_label,
            final_telo_connect,
        )
        final_selected_groups = [set(final_using_contig)]
        kept_primary_names = set(final_using_contig)
    else:
        (
            final_using_contig,
            final_using_type3_contig,
            final_ctg_typ,
            final_terminal_nodes,
            _,
        ) = alt_preprocess_contig(
            final_contig,
            final_telo_label,
            final_ref_qry_ratio,
            final_repeat_label,
            final_telo_connect,
            resources.telo_dict,
            resources.repeat_censat_data,
        )
        final_selected_groups = [
            set(final_using_contig),
            set(final_using_type3_contig),
        ]
        kept_primary_names = set()

    final_contig = final_contig[:-1]
    selected_contigs = []
    for selected_names in final_selected_groups:
        for row in final_contig:
            if row[CTG_NAM] not in selected_names:
                continue
            contig_type = final_ctg_typ[row[CTG_NAM]]
            if policy.kind == PafSourceKind.SECONDARY and contig_type == 5:
                if row[CTG_TELCON] == '0':
                    continue
                contig_type = 3
            row[CTG_TYP] = contig_type
            row[CTG_STRND] = final_terminal_nodes[row[CTG_NAM]][0] + index_offset
            row[CTG_ENDND] = final_terminal_nodes[row[CTG_NAM]][1] + index_offset
            selected_contigs.append(row)

    return PafPreprocessResult(
        contigs=selected_contigs,
        kept_primary_names=kept_primary_names,
        original_node_count=original_node_count,
        excluded_telomere_origins={
            (policy.source_index, row_idx) for row_idx in excluded_rows
        },
    )


def _extend_contigs_with_index_bias(contigs, added_contigs):
    """Append one locally indexed synthetic population at the current tail."""

    index_bias = len(contigs)
    for row in added_contigs:
        row[CTG_STRND] += index_bias
        row[CTG_ENDND] += index_bias
        contigs.append(row)


def _build_telomere_split_contigs(
    context,
    resources,
    contigs,
    excluded_telomere_origins,
    original_node_count,
):
    """Create one generation of legacy telomere-derived child contigs.

    Type-3/4, low-overlap, and simple-alt candidate generators intentionally do
    not participate in the minimal stage-01 assembly route.  Both telomere
    child populations are collected from the post-repeat-trim parent snapshot,
    processed once, and appended without recursively splitting the children.
    """

    telo_fb_dict = defaultdict(list)
    for chrom, intervals in resources.telo_dict.items():
        for interval in intervals:
            telo_fb_dict[chrom + interval[-1]].append([interval[0], interval[1]])

    telo_bound_dict = {}
    for telo_name, intervals in telo_fb_dict.items():
        if telo_name[-1] == 'f':
            telo_bound_dict[telo_name] = min(intervals, key=lambda value: value[0])
        else:
            telo_bound_dict[telo_name] = max(intervals, key=lambda value: value[1])

    telo_label = label_node(contigs, resources.telo_dict)
    subtelo_label = label_subtelo_node(contigs, resources.telo_dict)
    subtelo_ppc_node = subtelo_cut(contigs, telo_label, subtelo_label)

    adjacency = initial_graph_build(
        contigs,
        telo_bound_dict,
        resources.no_chrY,
    )
    (
        telo_connected_node,
        _,
        telo_connected_graph_dict,
        _,
    ) = edge_optimization(
        contigs,
        adjacency,
        telo_bound_dict,
        context.asm2cov,
        context.ori_ctg_name_data,
        excluded_telomere_origins,
    )
    (
        telo_connected_node,
        _,
        _,
    ) = filter_telomere_connected_cen_fragment_mismatch(
        contigs,
        telo_connected_graph_dict,
        resources.cen_fragment_meta,
        "pre-break",
    )

    break_contig = break_double_telomere_contig(contigs, telo_connected_node)

    final_break_contig = (
        pass_pipeline(
            break_contig,
            resources.telo_dict,
            telo_bound_dict,
            resources.repeat_data,
            resources.repeat_censat_data,
            False,
            resources.chr_len,
        )
        if break_contig
        else []
    )
    final_subtelo_ppc_node = (
        pass_pipeline(
            subtelo_ppc_node,
            resources.telo_dict,
            telo_bound_dict,
            resources.repeat_data,
            resources.repeat_censat_data,
            True,
            resources.chr_len,
        )
        if subtelo_ppc_node
        else []
    )
    _extend_contigs_with_index_bias(contigs, final_break_contig)
    _extend_contigs_with_index_bias(contigs, final_subtelo_ppc_node)

    add_node_count = len(final_break_contig) + len(final_subtelo_ppc_node)
    logging.info(f"Original PAF file length : {original_node_count}")
    logging.info(f"Final preprocessed PAF file length: {len(contigs)}")
    logging.info(
        "Number of telomere-derived contig nodes added on preprocessing : "
        f"{add_node_count}"
    )
    return contigs, telo_bound_dict


def _run_telomere_children_stage(context, state):
    """Default stage: add the two non-recursive telomere child populations."""

    state.contigs, state.telomere_bounds = _build_telomere_split_contigs(
        context.build,
        context.resources,
        state.contigs,
        set(context.excluded_telomere_origins),
        context.original_node_count,
    )
    state.stage_metadata["telomere_children"] = {
        "node_count": len(state.contigs),
    }
    return state


def default_contig_pipeline_stages():
    """Return the minimal stage-01 contig policy in biological order."""

    return (
        ContigPipelineStage(
            name="telomere_children",
            run=_run_telomere_children_stage,
        ),
    )


def run_contig_pipeline(context, state, stages):
    """Run an explicit sequence of unitig split/augmentation stages."""

    names = [stage.name for stage in stages]
    if len(names) != len(set(names)):
        raise ValueError(f"Duplicate contig pipeline stages: {names}")
    for stage in stages:
        before_count = len(state.contigs)
        updated = stage.run(context, state)
        if updated is not None:
            if not isinstance(updated, ContigPipelineState):
                raise TypeError(
                    f"Contig stage {stage.name!r} must return "
                    "ContigPipelineState or None"
                )
            state = updated
        after_count = len(state.contigs)
        state.stage_records.append({
            "name": stage.name,
            "before_nodes": before_count,
            "after_nodes": after_count,
            "added_nodes": max(after_count - before_count, 0),
            "removed_nodes": max(before_count - after_count, 0),
            "metadata": dict(state.stage_metadata.get(stage.name, {})),
        })
        logging.info(
            "Contig stage %s: %d -> %d node(s)",
            stage.name,
            before_count,
            after_count,
        )
    return state


def build_split_contigs(
    context,
    resources,
    contigs,
    excluded_telomere_origins,
    original_node_count,
    *,
    primary_kept_contig_names=(),
    stages=None,
    return_stage_records=False,
):
    """Run the configurable contig augmentation portion of stage 01.

    New non-telomere unitig splitters should be expressed as a
    :class:`ContigPipelineStage` and inserted relative to
    ``telomere_children``.  Existing callers retain the legacy two-value
    return contract.
    """

    pipeline_context = ContigPipelineContext(
        build=context,
        resources=resources,
        excluded_telomere_origins=frozenset(excluded_telomere_origins),
        original_node_count=int(original_node_count),
        primary_kept_contig_names=frozenset(primary_kept_contig_names),
    )
    state = ContigPipelineState(contigs=contigs, telomere_bounds=None)
    state = run_contig_pipeline(
        pipeline_context,
        state,
        default_contig_pipeline_stages() if stages is None else tuple(stages),
    )
    if state.telomere_bounds is None:
        raise ValueError(
            "Contig pipeline did not initialize telomere bounds; retain the "
            "telomere_children stage or provide equivalent bounds"
        )
    result = (state.contigs, state.telomere_bounds)
    if return_stage_records:
        return result + (tuple(state.stage_records),)
    return result


def build_simple_alt_contigs(
    context,
    resources,
    contigs,
    primary_kept_contig_names,
):
    """Append primary-PAF simple-alt contigs under the legacy policy."""

    if context.disable_alt_ctg_simple:
        return contigs

    existing_names = {row[CTG_NAM] for row in contigs}
    existing_nclose_loci = collect_existing_alt_simple_nclose_loci(contigs)
    chunks_per_contig = defaultdict(list)
    primary_paf_path = context.source.paf_file_paths[0]
    with open(primary_paf_path, "rt") as f:
        for line_idx, line in enumerate(f):
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            contig_name = cols[0]
            if contig_name in primary_kept_contig_names:
                continue
            try:
                chunk = [
                    contig_name,
                    int(cols[1]),
                    int(cols[2]),
                    int(cols[3]),
                    cols[4],
                    cols[5],
                    int(cols[6]),
                    int(cols[7]),
                    int(cols[8]),
                    int(cols[11]),
                    line_idx,
                ]
            except (IndexError, ValueError):
                continue
            chunks_per_contig[contig_name].append(chunk)

    simple_alt_count = 0
    simple_alt_skip_existing_count = 0
    simple_alt_skip_terminal_censat_count = 0
    simple_alt_keep_diff_dir_count = 0
    for contig_name, chunks in sorted(chunks_per_contig.items()):
        if len(chunks) < 2:
            continue
        chunks.sort(key=lambda row: (row[CTG_STR], row[CTG_END]))
        telo_labels = label_node(chunks, resources.telo_dict)
        repeat_labels = label_repeat_node(
            chunks,
            resources.repeat_data,
            resources.chr_len,
        )
        censat_labels = label_repeat_node(
            chunks,
            resources.repeat_censat_data,
            resources.chr_len,
        )
        trimmed_indices = trim_alt_simple_terminal_indices(
            chunks,
            telo_labels,
            repeat_labels,
            resources.chr_len,
        )
        if len(trimmed_indices) < 2:
            continue

        indexed_chunks = [
            (chunk_idx, chunks[chunk_idx]) for chunk_idx in trimmed_indices
        ]
        terminal_left = indexed_chunks[0][1]
        terminal_right = indexed_chunks[-1][1]
        if terminal_left[CHR_NAM] == terminal_right[CHR_NAM]:
            continue
        if alt_simple_terminal_censat_included(trimmed_indices, censat_labels):
            simple_alt_skip_terminal_censat_count += 1
            continue

        selected_chroms = select_alt_simple_major_chroms(
            indexed_chunks,
            required_chroms=(
                terminal_left[CHR_NAM],
                terminal_right[CHR_NAM],
            ),
        )
        candidate = find_alt_simple_inward_bounds(indexed_chunks, selected_chroms)
        if candidate is None:
            continue

        same_terminal_dir = (
            candidate[0][1][CTG_DIR] == candidate[1][1][CTG_DIR]
        )
        if same_terminal_dir:
            if alt_simple_candidate_near_existing(
                candidate,
                existing_nclose_loci,
                ALT_SIMPLE_EXISTING_NCLOSE_DIST,
            ):
                simple_alt_skip_existing_count += 1
                continue
        else:
            simple_alt_keep_diff_dir_count += 1

        synthetic_name = f"simple_ctg_alt_{contig_name}"
        if synthetic_name in existing_names:
            continue
        existing_names.add(synthetic_name)

        left_chunk_idx = candidate[0][0]
        right_chunk_idx = candidate[1][0]
        start_idx = len(contigs)
        end_idx = start_idx + right_chunk_idx - left_chunk_idx
        for chunk_idx in range(left_chunk_idx, right_chunk_idx + 1):
            raw = chunks[chunk_idx]
            line_idx = raw[10]
            synthetic_row = list(raw[:10]) + [
                1,
                start_idx,
                end_idx,
                telo_labels[chunk_idx][0],
                telo_labels[chunk_idx][1],
                '0',
                repeat_labels[chunk_idx][0],
                repeat_labels[chunk_idx][1],
                censat_labels[chunk_idx][1],
                '0',
                '0',
                f'0.{line_idx}',
            ]
            synthetic_row[CTG_NAM] = synthetic_name
            contigs.append(synthetic_row)

        if same_terminal_dir:
            existing_nclose_loci.append(alt_simple_pair_loci(candidate))
        simple_alt_count += 1

    if simple_alt_count > 0:
        logging.info(
            f"Number of simple ctg alt contigs added : {simple_alt_count}"
        )
    if simple_alt_skip_existing_count > 0:
        logging.info(
            "Number of simple ctg alt candidates skipped as existing nclose : "
            f"{simple_alt_skip_existing_count}"
        )
    if simple_alt_skip_terminal_censat_count > 0:
        logging.info(
            "Number of simple ctg alt contigs skipped because a "
            "post-terminal-trim endpoint is included in censat : "
            f"{simple_alt_skip_terminal_censat_count}"
        )
    if simple_alt_keep_diff_dir_count > 0:
        logging.info(
            "Number of simple ctg alt contigs with different terminal "
            f"directions kept unfiltered : {simple_alt_keep_diff_dir_count}"
        )
    return contigs


def finalize_preprocessed_contigs(
    context,
    resources,
    contigs,
    telo_bound_dict,
    excluded_telomere_origins,
):
    """Build final telomere state and persist the preprocessed PAF contract."""

    adjacency = initial_graph_build(
        contigs,
        telo_bound_dict,
        resources.no_chrY,
    )
    (
        _,
        _,
        telo_connected_graph_dict,
        telo_coverage,
    ) = edge_optimization(
        contigs,
        adjacency,
        telo_bound_dict,
        context.asm2cov,
        context.ori_ctg_name_data,
        excluded_telomere_origins,
    )
    (
        _,
        _,
        telo_connected_graph_dict,
    ) = filter_telomere_connected_cen_fragment_mismatch(
        contigs,
        telo_connected_graph_dict,
        resources.cen_fragment_meta,
        "final",
    )
    telo_edges = [
        (telo_name, tuple(edge))
        for telo_name, edges in telo_connected_graph_dict.items()
        for edge in edges
    ]
    virtual_telomere_nodes = add_missing_virtual_telomeres(
        contigs,
        telo_edges,
        resources.chr_len,
        resources.telo_data,
        resources.repeat_censat_data,
    )
    if virtual_telomere_nodes:
        logging.info(
            f"Added {virtual_telomere_nodes} missing chromosome-end "
            "virtual telomere nodes"
        )
    write_telomere_connected_outputs(context.prefix, telo_edges, contigs)

    real_final_mainflow = find_mainflow(contigs)
    for row in contigs:
        try:
            row[CTG_MAINFLOWDIR] = real_final_mainflow[row[CTG_NAM]][0]
            row[CTG_MAINFLOWCHR] = real_final_mainflow[row[CTG_NAM]][1]
        except Exception:
            pass

    with open(context.preprocessed_paf_path, "wt") as f:
        for row in contigs:
            for value in row:
                print(value, end="\t", file=f)
            print("", file=f)
    return telo_coverage


def contig_preprocessing_00(context, *, contig_stages=None):
    """Preprocess only the configured unitig source and its derived children.

    The stage-01 CLI retains the legacy primary PAF positional argument for
    command compatibility and PPC naming, but ``context.source`` contains
    only the effective unitig PAF in assembly mode.  Treating that single
    source with the legacy secondary policy preserves the unitig-specific
    selection rules without ever opening the primary PAF.
    """

    resources = load_contig_preprocess_resources(context)
    unitig_policy = PafPreprocessPolicy(
        kind=PafSourceKind.SECONDARY,
        source_index=0,
        global_index_prefix="0",
        log_label="unitig PAF",
    )
    unitig_result = preprocess_paf_source(
        context,
        resources,
        context.source.paf_file_paths[0],
        unitig_policy,
    )
    contigs = unitig_result.contigs
    excluded_telomere_origins = set(
        unitig_result.excluded_telomere_origins
    )
    original_node_count = unitig_result.original_node_count

    contigs, telo_bound_dict, stage_records = build_split_contigs(
        context,
        resources,
        contigs,
        excluded_telomere_origins,
        original_node_count,
        primary_kept_contig_names=(),
        stages=contig_stages,
        return_stage_records=True,
    )
    telo_coverage = finalize_preprocessed_contigs(
        context,
        resources,
        contigs,
        telo_bound_dict,
        excluded_telomere_origins,
    )
    return ContigPreprocessResult(
        telo_coverage=telo_coverage,
        depth_df=resources.depth_df,
        no_chrY=resources.no_chrY,
        stage_records=stage_records,
    )


def get_corr_dir(is_for : bool, dir_str : str) -> str:
    assert(dir_str == '+' or dir_str == '-')

    if is_for:
        return dir_str
    else:
        if dir_str == '+':
            return '-'
        else:
            return '+'

def get_left_right_centromere(repeat_censat_data : dict, chr_len : dict):
    left_start_cent = {}
    right_end_cent = {}
    for chrom, intervals in repeat_censat_data.items():
        chrom_length = chr_len[chrom]
        rep_start_0, rep_end_0 = intervals[0]  # 0-indexed 좌표
        if rep_start_0 == 0:
            left_start_cent[chrom] = rep_end_0
        elif rep_end_0 == chrom_length - 1:
            right_end_cent[chrom] = rep_start_0

    return left_start_cent, right_end_cent

def similar_centromere_nclose_cluster(nclose_dict : dict, contig_data : list, repeat_censat_data : dict, chr_len : dict):
    nclose_nodes = set()
    for j in nclose_dict:
        for k in nclose_dict[j]:
            nclose_nodes.add(k)

    left_cent, right_cent = get_left_right_centromere(repeat_censat_data, chr_len)
    centromere_nclose_master = set()
    slave_dict = dict()
    for (s, e) in nclose_nodes:
        contig_s = contig_data[s]
        contig_e = contig_data[e]
        if contig_s[CHR_NAM] in left_cent and contig_s[CHR_STR] <= left_cent[contig_s[CHR_NAM]] and contig_s[CTG_DIR] == '+':
            chk = True
            for (sm, em) in centromere_nclose_master:
                contig_sm = contig_data[sm]
                contig_em = contig_data[em]
                if contig_sm[CHR_NAM] in left_cent and contig_e[CHR_NAM] == contig_em[CHR_NAM] \
                        and contig_s[CTG_DIR] == contig_sm[CTG_DIR] and contig_e[CTG_DIR] == contig_em[CTG_DIR] \
                        and distance_checker(contig_e, contig_em) <= NCLOSE_COMPRESS_LIMIT:
                    slave_dict[(s, e)] = (sm, em)
                    chk = False
                    break
                elif contig_em[CHR_NAM] in left_cent and contig_e[CHR_NAM] == contig_sm[CHR_NAM] \
                        and contig_s[CTG_DIR] != contig_em[CTG_DIR] and contig_e[CTG_DIR] != contig_sm[CTG_DIR] \
                        and distance_checker(contig_e, contig_sm) <= NCLOSE_COMPRESS_LIMIT:
                    slave_dict[(s, e)] = (sm, em)
                    chk = False
                    break
            if chk:
                centromere_nclose_master.add((s, e))

        elif contig_e[CHR_NAM] in left_cent and contig_e[CHR_STR] <= left_cent[contig_e[CHR_NAM]] and contig_e[CTG_DIR] == '-':
            chk = True
            for (sm, em) in centromere_nclose_master:
                contig_sm = contig_data[sm]
                contig_em = contig_data[em]
                if contig_em[CHR_NAM] in left_cent and contig_s[CHR_NAM] == contig_sm[CHR_NAM] \
                        and contig_s[CTG_DIR] == contig_sm[CTG_DIR] and contig_e[CTG_DIR] == contig_em[CTG_DIR] \
                        and distance_checker(contig_s, contig_sm) <= NCLOSE_COMPRESS_LIMIT:
                    slave_dict[(s, e)] = (sm, em)
                    chk = False
                    break
                elif contig_sm[CHR_NAM] in left_cent and contig_s[CHR_NAM] == contig_em[CHR_NAM] \
                        and contig_s[CTG_DIR] != contig_em[CTG_DIR] and contig_e[CTG_DIR] != contig_sm[CTG_DIR] \
                        and distance_checker(contig_s, contig_em) <= NCLOSE_COMPRESS_LIMIT:
                    slave_dict[(s, e)] = (sm, em)
                    chk = False
                    break
            if chk:
                centromere_nclose_master.add((s, e))

        elif contig_s[CHR_NAM] in right_cent and contig_s[CHR_END] >= right_cent[contig_s[CHR_NAM]] and contig_s[CTG_DIR] == '-':
            chk = True
            for (sm, em) in centromere_nclose_master:
                contig_sm = contig_data[sm]
                contig_em = contig_data[em]
                if contig_sm[CHR_NAM] in right_cent and contig_e[CHR_NAM] == contig_em[CHR_NAM] \
                        and contig_s[CTG_DIR] == contig_sm[CTG_DIR] and contig_e[CTG_DIR] == contig_em[CTG_DIR] \
                        and distance_checker(contig_e, contig_em) <= NCLOSE_COMPRESS_LIMIT:
                    slave_dict[(s, e)] = (sm, em)
                    chk = False
                    break
                elif contig_em[CHR_NAM] in right_cent and contig_e[CHR_NAM] == contig_sm[CHR_NAM] \
                        and contig_s[CTG_DIR] != contig_em[CTG_DIR] and contig_e[CTG_DIR] != contig_sm[CTG_DIR] \
                        and distance_checker(contig_e, contig_sm) <= NCLOSE_COMPRESS_LIMIT:
                    slave_dict[(s, e)] = (sm, em)
                    chk = False
                    break
            if chk:
                centromere_nclose_master.add((s, e))
        elif contig_e[CHR_NAM] in right_cent and contig_e[CHR_END] >= right_cent[contig_e[CHR_NAM]] and contig_e[CTG_DIR] == '+':
            chk = True
            for (sm, em) in centromere_nclose_master:
                contig_sm = contig_data[sm]
                contig_em = contig_data[em]
                if contig_em[CHR_NAM] in right_cent and contig_s[CHR_NAM] == contig_sm[CHR_NAM] \
                        and contig_s[CTG_DIR] == contig_sm[CTG_DIR] and contig_e[CTG_DIR] == contig_em[CTG_DIR] \
                        and distance_checker(contig_s, contig_sm) <= NCLOSE_COMPRESS_LIMIT:
                    slave_dict[(s, e)] = (sm, em)
                    chk = False
                    break
                elif contig_sm[CHR_NAM] in right_cent and contig_s[CHR_NAM] == contig_em[CHR_NAM] \
                        and contig_s[CTG_DIR] != contig_em[CTG_DIR] and contig_e[CTG_DIR] != contig_sm[CTG_DIR] \
                        and distance_checker(contig_s, contig_em) <= NCLOSE_COMPRESS_LIMIT:
                    slave_dict[(s, e)] = (sm, em)
                    chk = False
                    break

            if chk:
                centromere_nclose_master.add((s, e))

    return centromere_nclose_master, slave_dict

def conjoined_type4(contig_data, type2_nclose_node):
    """
    Detect type4 insertions/deletions by pairing type2 nodes within each chromosome.
    If a contig's internal span >= 100 kb, also test its flipped orientation in a nested loop.
    """
    conjoined_type4_ins = set()
    conjoined_type4_del = set()

    def flip_dir(ctg):
        flipped = list(ctg)
        flipped[CTG_DIR] = '+' if ctg[CTG_DIR] == '-' else '-'
        return flipped

    for chrom, type2_list in type2_nclose_node.items():
        L = len(type2_list)
        for i in range(L):
            for j in range(i + 1, L):
                t2n1 = type2_list[i]
                t2n2 = type2_list[j]

                c1fi, c1bi = t2n1
                c2fi, c2bi = t2n2

                c1f = contig_data[c1fi]
                c1b = contig_data[c1bi]
                c2f = contig_data[c2fi]
                c2b = contig_data[c2bi]

                # Internal spans; switch to CHR_END if that's your convention
                len_c1 = abs(c1f[CHR_STR] - c1b[CHR_STR])
                len_c2 = abs(c2f[CHR_STR] - c2b[CHR_STR])

                # Candidate orientations per contig (original, and flipped if span >= threshold)
                c1_flips = [(c1f, c1fi, c1b, c1bi)]
                c2_flips = [(c2f, c2fi, c2b, c2bi)]

                if len_c1 >= TYPE2_DIST_FLIP_THRESHOLD:
                    c1_flips.append((flip_dir(c1b), c1bi, flip_dir(c1f), c1fi))
                if len_c2 >= TYPE2_DIST_FLIP_THRESHOLD:
                    c2_flips.append((flip_dir(c2b), c2bi, flip_dir(c2f), c2fi))

                # Test all combinations (original/flip × original/flip)
                for c1f_use, c1fi_, c1b_use, c1bi_ in c1_flips:
                    for c2f_use, c2fi_, c2b_use, c2bi_ in c2_flips:
                        # Forward direction
                        if c1b_use[CTG_DIR] == c2f_use[CTG_DIR]:
                            dist_for = distance_checker(c1b_use, c2f_use)
                            if dist_for is not None and dist_for < TYPE2_CONJOIN_COMPRESS_LIMIT:
                                ratio, total_ref_len = calculate_single_contig_ref_ratio([c1f_use, c2b_use])
                                abs_estimated = total_ref_len * abs(ratio)
                                if abs(ratio - 1) > BND_CONTIG_BOUND and abs_estimated > TYPE2_CHUKJI_AS_TYPE4:
                                    # Insertion
                                    if ratio < 0:
                                        if contig_data[c1fi_][CHR_STR] > contig_data[c2bi_][CHR_STR]:
                                            conjoined_type4_ins.add((c1fi_, c1bi_, c2fi_, c2bi_))
                                        else:
                                            conjoined_type4_ins.add((c2bi_, c2fi_, c1bi_, c1fi_))
                                    # Deletion
                                    else:
                                        if contig_data[c1fi_][CHR_STR] < contig_data[c2bi_][CHR_STR]:
                                            conjoined_type4_del.add((c1fi_, c1bi_, c2fi_, c2bi_))
                                        else:
                                            conjoined_type4_del.add((c2bi_, c2fi_, c1bi_, c1fi_))


    conjoined_type4_ins = sorted(conjoined_type4_ins)
    conjoined_type4_del = sorted(conjoined_type4_del)

    return conjoined_type4_ins, conjoined_type4_del

def get_breakend_coord(contig, side_idx):
    nclose_loc = side_idx == 0
    nclose_dir = contig[CTG_DIR] == '+'
    return contig[CHR_STR if nclose_dir ^ nclose_loc else CHR_END]

def get_expected_high_side_from_contig(contig, side_idx):
    nclose_loc = side_idx == 0
    nclose_dir = contig[CTG_DIR] == '+'
    return 'right' if nclose_dir ^ nclose_loc else 'left'

def get_raw_endpoint(contig, side_idx, contig_idx):
    nclose_loc = side_idx == 0
    is_front_dir = contig[CTG_DIR] == '+'
    return {
        'chrom': contig[CHR_NAM],
        'coord': int(get_breakend_coord(contig, side_idx)),
        'dir': contig[CTG_DIR],
        'match_forward': bool(is_front_dir == nclose_loc),
        'contig_idx': int(contig_idx),
        'ctg_name': contig[CTG_NAM],
        'ctg_st': int(contig[CTG_STR]),
        'ctg_nd': int(contig[CTG_END]),
        'ref_st': int(contig[CHR_STR]),
        'ref_nd': int(contig[CHR_END]),
        'expected_high_side': get_expected_high_side_from_contig(contig, side_idx),
    }

def canonical_raw_nclose_layout(contig_data, nclose_pair):
    ordered_endpoints = [
        get_raw_endpoint(contig_data[contig_idx], side_idx, contig_idx)
        for side_idx, contig_idx in enumerate(nclose_pair)
    ]

    order = [0, 1]
    keys = [(chr2int(ordered_endpoints[i]['chrom']), ordered_endpoints[i]['coord']) for i in order]
    if keys[1] < keys[0]:
        order = [1, 0]

    endpoints = [ordered_endpoints[i] for i in order]
    return {
        'nclose_key': tuple(sorted(int(x) for x in nclose_pair)),
        'chroms': tuple(ep['chrom'] for ep in endpoints),
        'coords': tuple(ep['coord'] for ep in endpoints),
        'sides': tuple(ep['expected_high_side'] for ep in endpoints),
        'endpoints': tuple(endpoints),
        'ordered_endpoints': tuple(ordered_endpoints),
    }

def get_inner_boundary_interval(region_a, region_b):
    st_a, nd_a = int(region_a['ref_st']), int(region_a['ref_nd'])
    st_b, nd_b = int(region_b['ref_st']), int(region_b['ref_nd'])

    if (st_b, nd_b) < (st_a, nd_a):
        st_a, nd_a, st_b, nd_b = st_b, nd_b, st_a, nd_a

    if nd_a <= st_b:
        return nd_a, st_b

    return st_b, min(nd_a, nd_b)

def get_outer_reference_interval(region_a, region_b):
    return (
        min(int(region_a['ref_st']), int(region_b['ref_st'])),
        max(int(region_a['ref_nd']), int(region_b['ref_nd'])),
    )

def centered_reference_interval(chrom, inner_st, inner_nd, chr_len, size=RAW_TRANSLOCATION_WINDOW):
    center = int(round((int(inner_st) + int(inner_nd)) / 2))
    chrom_len = int(chr_len.get(chrom, center + size))
    st = center - size // 2
    nd = st + size

    if st < 0:
        st = 0
        nd = min(size, chrom_len)
    elif nd > chrom_len:
        nd = chrom_len
        st = max(0, nd - size)

    return int(st), int(nd)

def layout_has_split_contig_endpoint(layout):
    return any(
        str(endpoint['ctg_name']).startswith('split_contig_')
        for endpoint in layout['ordered_endpoints']
    )

def layout_endpoint_dirs_are_opposite(layout):
    endpoints = layout['ordered_endpoints']
    return len(endpoints) == 2 and endpoints[0]['dir'] != endpoints[1]['dir']

def is_large_same_chrom_raw_candidate(chrom_pair, layout_a, layout_b):
    if chrom_pair[0] != chrom_pair[1]:
        return False
    if not layout_endpoint_dirs_are_opposite(layout_a):
        return False
    if not layout_endpoint_dirs_are_opposite(layout_b):
        return False
    if layout_has_split_contig_endpoint(layout_a) or layout_has_split_contig_endpoint(layout_b):
        return False

    span_a = abs(int(layout_a['coords'][1]) - int(layout_a['coords'][0]))
    span_b = abs(int(layout_b['coords'][1]) - int(layout_b['coords'][0]))
    return max(span_a, span_b) >= RAW_TRANSLOCATION_MIN_SAME_CHROM_SPAN

def node_fully_inside_censat(node, repeat_censat_data):
    for censat_st, censat_nd in repeat_censat_data.get(node[CHR_NAM], []):
        if censat_st <= node[CHR_STR] and node[CHR_END] <= censat_nd:
            return True
    return False

def extract_raw_censat_type2_candidates(paf_path, repeat_censat_data):
    candidates = []
    if paf_path is None or not os.path.exists(paf_path):
        return candidates

    chunks_per_contig = defaultdict(list)
    with open(paf_path, "rt") as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            try:
                chunk = [
                    cols[0],
                    int(cols[1]),
                    int(cols[2]),
                    int(cols[3]),
                    cols[4],
                    cols[5],
                    int(cols[6]),
                    int(cols[7]),
                    int(cols[8]),
                    int(cols[11]),
                ]
            except (IndexError, ValueError):
                continue
            chunks_per_contig[chunk[CTG_NAM]].append(chunk)

    for ctg_name, chunks in chunks_per_contig.items():
        if len(chunks) < 2:
            continue
        chunks.sort(key=lambda c: (c[CTG_STR], c[CTG_END]))
        if any(chunk[CTG_MAPQ] <= 0 for chunk in chunks):
            continue
        if any(not node_fully_inside_censat(chunk, repeat_censat_data) for chunk in chunks):
            continue

        start = chunks[0]
        end = chunks[-1]
        if start[CHR_NAM] != end[CHR_NAM]:
            continue
        if start[CTG_DIR] == end[CTG_DIR]:
            continue

        candidates.append({
            "ctg_name": ctg_name,
            "chrom": start[CHR_NAM],
            "pair": (tuple(start), tuple(end)),
        })

    return candidates


def extract_source_censat_type2_candidates(source, repeat_censat_data):
    """Read optional type2 evidence only from the active source contract."""

    if source.secondary_candidate_paf is None:
        return []
    return extract_raw_censat_type2_candidates(
        source.secondary_candidate_paf,
        repeat_censat_data,
    )

def connected_type2_endpoint_for_template(template_contig, template_side, type2_pair):
    type2_front, type2_back = type2_pair
    if template_side == 0:
        if template_contig[CTG_DIR] == '+':
            dist = distance_checker(type2_back, template_contig)
            if type2_back[CHR_END] < template_contig[CHR_STR] or dist == 0:
                return type2_back, dist
        else:
            dist = distance_checker(type2_front, template_contig)
            if type2_front[CHR_END] > template_contig[CHR_STR] or dist == 0:
                return type2_front, dist
    else:
        if template_contig[CTG_DIR] == '+':
            dist = distance_checker(type2_front, template_contig)
            if type2_front[CHR_STR] > template_contig[CHR_END] or dist == 0:
                return type2_front, dist
        else:
            dist = distance_checker(type2_back, template_contig)
            if type2_back[CHR_STR] < template_contig[CHR_END] or dist == 0:
                return type2_back, dist
    return None, None

def append_ppc_rows(ppc_path, rows):
    if not rows:
        return
    with open(ppc_path, "at") as f:
        for row in rows:
            for value in row:
                print(value, end="\t", file=f)
            print("", file=f)


def collect_missing_cen_fragment_dir_censat_noncensat(
    contig_data,
    nclose_nodes,
    cen_fragment_meta,
):
    """Return one-CEN-SAT offset clusters lacking the depth-supported direction."""

    candidates = [
        candidate
        for candidate in iter_censat_noncensat_candidates(
            contig_data,
            nclose_nodes,
        )
        if not candidate.is_simple_alt
        and candidate.censat_chrom in cen_fragment_meta
    ]

    missing = []
    for group in group_censat_noncensat_candidates(candidates):
        target_norm_dir = _cen_fragment_target_dir_from_meta(
            cen_fragment_meta,
            group[0].censat_chrom,
        )
        if all(
            candidate.censat_norm_dir != target_norm_dir
            for candidate in group
        ):
            missing.append(group)

    return missing

def apply_offset_direction_mismatched_censat_noncensat_filter(
    candidates,
    contig_data,
    cen_fragment_meta,
):
    """Prefer the depth-supported direction in stable candidate order."""

    one_censat_candidates = [
        one_censat
        for one_censat in iter_censat_noncensat_candidates(
            contig_data,
            candidates,
        )
        if one_censat.censat_chrom in cen_fragment_meta
    ]

    offset_to_remove = set()
    for group in group_censat_noncensat_candidates(one_censat_candidates):
        directions = {candidate.censat_norm_dir for candidate in group}
        if len(directions) < 2:
            continue
        target_norm_dir = _cen_fragment_target_dir_from_meta(
            cen_fragment_meta,
            group[0].censat_chrom,
        )
        offset_to_remove.update(
            (candidate.contig_name, candidate.pair)
            for candidate in group
            if candidate.censat_norm_dir != target_norm_dir
        )

    return apply_nclose_filter(
        candidates,
        "censat_noncensat_offset_direction",
        lambda candidate: (
            "direction_mismatch"
            if candidate.identity in offset_to_remove
            else None
        ),
    )


def filter_offset_direction_mismatched_censat_noncensat(
    contig_data,
    nclose_nodes,
    cen_fragment_meta,
):
    """Compatibility adapter for callers that still own the legacy mapping."""

    candidates = [
        NCloseCandidate(contig_name=contig_name, path_pair=pair)
        for contig_name, pair in iter_nclose_owner_pairs(nclose_nodes)
    ]
    kept, rejections = apply_offset_direction_mismatched_censat_noncensat_filter(
        candidates,
        contig_data,
        cen_fragment_meta,
    )
    filtered = defaultdict(list)
    for candidate in kept:
        filtered[candidate.contig_name].append(candidate.path_pair)
    return filtered, len(rejections)


def add_nearest_combined_censat_noncensat_ncloses(
    contig_data,
    missing_candidate_groups,
    raw_censat_type2_candidates,
    repeat_censat_data,
    chr_len,
    cen_fragment_meta,
    ppc_path,
):
    if not missing_candidate_groups or not raw_censat_type2_candidates:
        return []

    type2_by_chrom = defaultdict(list)
    for candidate in raw_censat_type2_candidates:
        type2_by_chrom[candidate["chrom"]].append(candidate)

    rows_to_append = []
    seen_signature = set()
    added_candidates = []

    for group in missing_candidate_groups:
        best = None
        for cand in group:
            target_dir = _cen_fragment_target_dir_from_meta(
                cen_fragment_meta,
                cand.censat_chrom,
            )
            censat_node = contig_data[cand.censat_idx]
            for type2 in type2_by_chrom.get(cand.censat_chrom, []):
                type2_endpoint, dist = connected_type2_endpoint_for_template(
                    censat_node,
                    cand.censat_side,
                    type2["pair"],
                )
                if type2_endpoint is None:
                    continue
                if type2_endpoint[CTG_DIR] != target_dir:
                    continue
                best_key = (
                    dist,
                    abs(censat_node[CHR_STR] - type2_endpoint[CHR_STR]),
                    cand.noncensat_pos,
                )
                if best is None or best_key < best[0]:
                    best = (best_key, cand, type2, type2_endpoint)

        if best is None:
            continue

        _, cand, type2, remaining_endpoint = best
        noncensat_node = contig_data[cand.noncensat_idx]
        # The synthetic contig is written as censat -> noncensat. If the
        # original nclose had censat on the second side, preserve the
        # noncensat endpoint in that reverse-complemented frame.
        noncensat_dir = noncensat_node[CTG_DIR]
        if cand.censat_side == 1:
            noncensat_dir = _flip_ctg_dir(noncensat_dir)
        signature = (
            remaining_endpoint[CTG_NAM],
            remaining_endpoint[CHR_NAM],
            remaining_endpoint[CHR_STR],
            remaining_endpoint[CHR_END],
            remaining_endpoint[CTG_DIR],
            cand.noncensat_idx,
            noncensat_dir,
        )
        if signature in seen_signature:
            continue
        seen_signature.add(signature)

        added = len(added_candidates)
        new_name = f"combined_{type2['ctg_name']}_{noncensat_node[CTG_NAM]}_{added}"
        new_idx0 = len(contig_data)
        new_idx1 = new_idx0 + 1
        len0 = remaining_endpoint[CHR_END] - remaining_endpoint[CHR_STR]
        len1 = noncensat_node[CHR_END] - noncensat_node[CHR_STR]
        total_len = len0 + len1
        sv_type = 1 if remaining_endpoint[CHR_NAM] != noncensat_node[CHR_NAM] else 2
        censat_label = label_repeat_node([remaining_endpoint], repeat_censat_data, chr_len)[0][1]

        node0 = [
            new_name,
            total_len,
            0,
            len0,
            remaining_endpoint[CTG_DIR],
            remaining_endpoint[CHR_NAM],
            remaining_endpoint[CHR_LEN],
            remaining_endpoint[CHR_STR],
            remaining_endpoint[CHR_END],
            remaining_endpoint[CTG_MAPQ],
            sv_type,
            new_idx0,
            new_idx1,
            '0',
            '0',
            '0',
            '0',
            '0',
            censat_label,
            remaining_endpoint[CTG_DIR],
            remaining_endpoint[CHR_NAM],
            f'4.{2 * added}',
        ]
        node1 = [
            new_name,
            total_len,
            len0,
            total_len,
            noncensat_dir,
            noncensat_node[CHR_NAM],
            noncensat_node[CHR_LEN],
            noncensat_node[CHR_STR],
            noncensat_node[CHR_END],
            noncensat_node[CTG_MAPQ],
            sv_type,
            new_idx0,
            new_idx1,
            '0',
            '0',
            '0',
            noncensat_node[CTG_RPTCHR],
            noncensat_node[CTG_RPTCASE],
            noncensat_node[CTG_CENSAT],
            noncensat_dir,
            noncensat_node[CHR_NAM],
            f'4.{2 * added + 1}',
        ]

        contig_data.append(tuple(node0))
        contig_data.append(tuple(node1))
        rows_to_append.extend([node0, node1])
        added_candidates.append(NCloseCandidate(
            contig_name=new_name,
            path_pair=(new_idx0, new_idx1),
            origin="combined_censat_noncensat",
            synthetic=True,
            signatures={"synthetic_dedup": signature},
            provenance={
                "parent_identity": (cand.contig_name, cand.pair),
                "type2_contig": type2["ctg_name"],
            },
        ))

    append_ppc_rows(ppc_path, rows_to_append)
    return added_candidates


def build_raw_translocation_candidates(contig_data, nclose_source, chr_len,
                                       distance_threshold=RAW_TRANSLOCATION_WINDOW,
                                       candidate_filter=None):
    layout_by_key = {}
    for _, pair in iter_nclose_owner_pairs(nclose_source):
        nclose_key = tuple(sorted(int(x) for x in pair))
        if len(nclose_key) != 2:
            continue
        layout_by_key[nclose_key] = canonical_raw_nclose_layout(
            contig_data,
            tuple(int(x) for x in pair),
        )

    groups = defaultdict(list)
    for nclose_key, layout in layout_by_key.items():
        groups[layout['chroms']].append((nclose_key, layout))

    candidates = []
    seen_candidate_signatures = set()
    pair_id = 0
    for chrom_pair, items in groups.items():
        for i in range(len(items)):
            key_i, layout_i = items[i]
            for j in range(i + 1, len(items)):
                key_j, layout_j = items[j]
                if layout_i['sides'][0] == layout_j['sides'][0] or layout_i['sides'][1] == layout_j['sides'][1]:
                    continue
                if abs(layout_i['coords'][0] - layout_j['coords'][0]) > distance_threshold:
                    continue
                if abs(layout_i['coords'][1] - layout_j['coords'][1]) > distance_threshold:
                    continue

                if layout_i['sides'][0] == 'left':
                    key_a, layout_a = key_i, layout_i
                    key_b, layout_b = key_j, layout_j
                else:
                    key_a, layout_a = key_j, layout_j
                    key_b, layout_b = key_i, layout_i

                if candidate_filter is not None and not candidate_filter(chrom_pair, layout_a, layout_b):
                    continue

                candidate_signature = (
                    chrom_pair,
                    layout_a['coords'],
                    layout_a['sides'],
                    layout_b['coords'],
                    layout_b['sides'],
                )
                if candidate_signature in seen_candidate_signatures:
                    continue
                seen_candidate_signatures.add(candidate_signature)

                side_records = []
                for side_idx, chrom in enumerate(chrom_pair):
                    inner_st, inner_nd = get_inner_boundary_interval(
                        layout_a['endpoints'][side_idx],
                        layout_b['endpoints'][side_idx],
                    )
                    ref_st, ref_nd = centered_reference_interval(chrom, inner_st, inner_nd, chr_len)
                    path_drop_st, path_drop_nd = get_outer_reference_interval(
                        layout_a['endpoints'][side_idx],
                        layout_b['endpoints'][side_idx],
                    )
                    side_records.append({
                        'side_id': f'R{pair_id}.{side_idx}',
                        'chrom': chrom,
                        'inner_st': int(inner_st),
                        'inner_nd': int(inner_nd),
                        'ref_count_st': int(ref_st),
                        'ref_count_nd': int(ref_nd),
                        'path_drop_st': int(path_drop_st),
                        'path_drop_nd': int(path_drop_nd),
                        'no_spanning_rawread': False,
                        'no_spanning_utg': False,
                        'crossing_cols': [],
                        'accepted_forbid': False,
                    })

                candidates.append({
                    'pair_id': pair_id,
                    'nclose_key_a': key_a,
                    'nclose_key_b': key_b,
                    'chrom_pair': chrom_pair,
                    'layout_a': layout_a,
                    'layout_b': layout_b,
                    'split_junctions': [
                        {'count_name': 'd1', 'count_idx': 1, 'nclose_key': key_a, 'endpoints': layout_a['ordered_endpoints']},
                        {'count_name': 'd4', 'count_idx': 4, 'nclose_key': key_b, 'endpoints': layout_b['ordered_endpoints']},
                    ],
                    'span_junctions': [
                        {'count_name': 'd2', 'count_idx': 2, 'side_idx': 0, 'side_id': side_records[0]['side_id']},
                        {'count_name': 'd3', 'count_idx': 3, 'side_idx': 1, 'side_id': side_records[1]['side_id']},
                    ],
                    'side_records': side_records,
                    'pair_protected_from_1M_drop': False,
                })
                pair_id += 1

    return candidates

def raw_translocation_candidate_signature(candidate):
    return (
        candidate['chrom_pair'],
        candidate['layout_a']['coords'],
        candidate['layout_a']['sides'],
        candidate['layout_b']['coords'],
        candidate['layout_b']['sides'],
    )

def renumber_raw_translocation_candidates(candidates):
    for pair_id, candidate in enumerate(candidates):
        candidate['pair_id'] = pair_id
        for side_idx, side_record in enumerate(candidate['side_records']):
            side_record['side_id'] = f'R{pair_id}.{side_idx}'
        for span_junction in candidate['span_junctions']:
            side_idx = span_junction.get('side_idx')
            if side_idx is not None and 0 <= side_idx < len(candidate['side_records']):
                span_junction['side_id'] = candidate['side_records'][side_idx]['side_id']
    return candidates

def merge_raw_translocation_candidates(*candidate_lists):
    merged = []
    seen = set()
    for candidate_list in candidate_lists:
        for candidate in candidate_list:
            signature = raw_translocation_candidate_signature(candidate)
            if signature in seen:
                continue
            seen.add(signature)
            merged.append(candidate)
    return renumber_raw_translocation_candidates(merged)

def build_nclose_count_candidates(contig_data, nclose_source):
    candidates = []
    seen = set()
    pair_id = 0
    for ctg_name, pair in iter_nclose_owner_pairs(nclose_source):
        nclose_key = tuple(sorted(int(x) for x in pair))
        if len(nclose_key) != 2:
            continue
        if nclose_key in seen:
            continue
        seen.add(nclose_key)
        layout = canonical_raw_nclose_layout(
            contig_data,
            tuple(int(x) for x in pair),
        )
        overlaps_censat = any(
            node_is_censat(contig_data[int(idx)])
            for idx in nclose_key
        )
        candidates.append({
            'pair_id': pair_id,
            'nclose_key': nclose_key,
            'ctg_name': ctg_name,
            'layout': layout,
            'overlaps_censat': overlaps_censat,
        })
        pair_id += 1
    return candidates

def collect_nclose_count_keep_pairs(prefix, vaf_threshold):
    result_path = f'{prefix}/{NCLOSE_COUNT_RESULT_PKL}'
    if not os.path.isfile(result_path):
        logging.warning(f'NClose raw-count result not found: {result_path}')
        return None, 0, 0

    with open(result_path, 'rb') as f:
        records = pkl.load(f)

    keep_pairs = set()
    for record in records:
        if 'keep_nclose' in record:
            if record['keep_nclose']:
                keep_pairs.add(tuple(sorted(int(x) for x in record['nclose_key'])))
            continue
        if record.get('filter_eligible') is False:
            keep_pairs.add(tuple(sorted(int(x) for x in record['nclose_key'])))
            continue
        vaf = record.get('vaf', {})
        chr_a_vaf = vaf.get('chr_a')
        chr_b_vaf = vaf.get('chr_b')
        if (
            (chr_a_vaf is not None and chr_a_vaf >= vaf_threshold) or
            (chr_b_vaf is not None and chr_b_vaf >= vaf_threshold)
        ):
            keep_pairs.add(tuple(sorted(int(x) for x in record['nclose_key'])))

    return keep_pairs, len(records), len(keep_pairs)

def raw_false_value(value):
    return value is False or value == 0 or value == 'False' or value == 'false'

def raw_record_has_both_point_spans(record):
    no_span = record.get('raw_point_no_spanning')
    if isinstance(no_span, dict):
        return raw_false_value(no_span.get('point_a')) and raw_false_value(no_span.get('point_b'))

    side_records = record.get('side_records', [])
    if len(side_records) < 2:
        return False
    side_flags = [
        side.get('raw_point_no_spanning', side.get('no_spanning_rawread'))
        for side in side_records[:2]
    ]
    return raw_false_value(side_flags[0]) and raw_false_value(side_flags[1])

def raw_record_is_depth_balanced(record):
    if 'depth_balanced_translocation' not in record:
        return True
    return bool(record.get('depth_balanced_translocation'))

def collect_raw_virtual_inv_nclose_pairs(prefix):
    result_path = f'{prefix}/{RAW_TRANSLOCATION_RESULT_PKL}'
    if not os.path.isfile(result_path):
        return set()

    with open(result_path, 'rb') as f:
        records = pkl.load(f)

    pairs_to_remove = set()
    for record in records:
        if not raw_record_has_both_point_spans(record):
            continue
        if not raw_record_is_depth_balanced(record):
            continue
        for key_name in ('nclose_key_a', 'nclose_key_b'):
            key = record.get(key_name)
            if key is None:
                continue
            pairs_to_remove.add(tuple(sorted(int(x) for x in key)))

    return pairs_to_remove

















VCF_SYNTHETIC_FLANK = 1 * K
VCF_SYNTHETIC_PAF_NAME = "vcf_synthetic.paf"
VCF_SKIPPED_RECORDS_TSV = "vcf_mode_skipped_records.tsv"
VCF_ORIENTATION_MISMATCH_TSV = "vcf_mode_orientation_mismatches.tsv"
VCF_TELOMERE_PAF_NODES_TSV = "vcf_telomere_paf_nodes.tsv"






def load_vcf_ins_alt_alignment_spans(paf_path):
    best_row_by_query = {}
    best_key_by_query = {}
    rows_by_query = defaultdict(list)
    if not paf_path:
        return {}
    if not os.path.isfile(paf_path):
        logging.warning(f"VCF INS --alt PAF does not exist: {paf_path}")
        return {}

    with open(paf_path, "rt") as paf:
        for line in paf:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 12:
                continue
            try:
                qname = cols[0]
                chrom = cols[5]
                st, nd = sorted((int(cols[7]), int(cols[8])))
                matches = int(cols[9])
                block_len = int(cols[10])
                mapq = int(cols[11])
            except (IndexError, ValueError):
                continue
            if nd <= st:
                continue
            rows_by_query[qname].append(cols)
            # Keep every row for the INS type4 base PAF, but derive the
            # required single-chromosome representative coordinates from the
            # largest individual alignment block only.
            row_key = (block_len, mapq, matches, nd - st)
            if qname not in best_key_by_query or row_key > best_key_by_query[qname]:
                best_key_by_query[qname] = row_key
                best_row_by_query[qname] = {
                    "chrom": chrom,
                    "st": st,
                    "nd": nd,
                    "score": block_len,
                    "rows": 1,
                    "mapq": mapq,
                }

    best_by_query = {}
    for qname, best in best_row_by_query.items():
        if best["chrom"] is not None and best["nd"] > best["st"]:
            best = dict(best)
            best["paf_rows"] = rows_by_query.get(qname, [])
            best_by_query[qname] = best
    return best_by_query


def endpoint_interval(pos, chrom_len, path_dir, endpoint_role):
    pos = max(1, min(int(pos), int(chrom_len)))
    flank = min(VCF_SYNTHETIC_FLANK, max(1, int(chrom_len)))
    if endpoint_role == "exit":
        if path_dir == "+":
            ref_end = pos
            ref_st = max(0, ref_end - flank)
        else:
            ref_st = pos
            ref_end = min(int(chrom_len), ref_st + flank)
    else:
        if path_dir == "+":
            ref_st = pos
            ref_end = min(int(chrom_len), ref_st + flank)
        else:
            ref_end = pos
            ref_st = max(0, ref_end - flank)

    if ref_end <= ref_st:
        if ref_st <= 0:
            ref_st, ref_end = 0, 1
        else:
            ref_st, ref_end = ref_st - 1, ref_st
    return int(ref_st), int(ref_end)


def make_synthetic_vcf_node(contig_name, chrom, pos, path_dir, endpoint_role,
                            chr_len, ctg_typ, idx, global_idx):
    ref_st, ref_end = endpoint_interval(pos, chr_len[chrom], path_dir, endpoint_role)
    ref_len = max(1, ref_end - ref_st)
    return [
        contig_name,
        ref_len,
        0,
        ref_len,
        path_dir,
        chrom,
        int(chr_len[chrom]),
        ref_st,
        ref_end,
        60,
        ctg_typ,
        idx,
        idx,
        "0",
        "0",
        "0",
        "0",
        "0",
        "0",
        path_dir,
        chrom,
        f"2.{global_idx}",
    ]


def make_synthetic_span_node(contig_name, chrom, st, nd, path_dir,
                             chr_len, ctg_typ, idx, global_idx):
    st, nd = sorted((int(st), int(nd)))
    if nd <= st:
        nd = st + 1
    st = max(0, min(st, int(chr_len[chrom]) - 1))
    nd = max(st + 1, min(nd, int(chr_len[chrom])))
    ref_len = nd - st
    return [
        contig_name,
        ref_len,
        0,
        ref_len,
        path_dir,
        chrom,
        int(chr_len[chrom]),
        st,
        nd,
        60,
        ctg_typ,
        idx,
        idx,
        "0",
        "0",
        "0",
        "0",
        "0",
        "0",
        path_dir,
        chrom,
        f"2.{global_idx}",
    ]


def synthetic_node_to_paf_row(node):
    ref_len = max(1, int(node[CHR_END]) - int(node[CHR_STR]))
    return list(node[:9]) + [ref_len, ref_len, int(node[CTG_MAPQ]), "tp:A:P", f"cs:Z::{ref_len}"]


def telo_data_to_dict(telo_data):
    telo_dict = defaultdict(list)
    for entry in telo_data:
        telo_dict[entry[0]].append(entry[1:])
    return telo_dict


def telomere_label_to_name(label):
    if label[0] == "0" or label[1] == "0":
        return None
    side = str(label[1])[0]
    if side not in {"f", "b"}:
        return None
    return f"{label[0]}{side}"


def no_repeat_labels(contig_data):
    return [("0", "0") for _ in contig_data]


def telomere_bounds_by_name(telo_dict):
    telo_bounds = {}
    for chrom, intervals in telo_dict.items():
        intervals_by_side = defaultdict(list)
        for st, nd, side in intervals:
            intervals_by_side[side].append((int(st), int(nd)))
        for side, side_intervals in intervals_by_side.items():
            if side == "f":
                telo_bounds[f"{chrom}{side}"] = min(
                    side_intervals,
                    key=lambda interval: interval[0],
                )
            else:
                telo_bounds[f"{chrom}{side}"] = max(
                    side_intervals,
                    key=lambda interval: interval[1],
                )
    return telo_bounds


def is_terminal_telomere_candidate(row, telo_name, telo_bounds, chr_len):
    if telo_name.endswith("f"):
        return int(row[CHR_STR]) <= TELOMERE_CLUSTER_THRESHOLD

    telo_bound = telo_bounds.get(telo_name)
    if telo_bound is not None:
        telo_end = int(telo_bound[1])
    else:
        telo_end = int(chr_len[telo_name[:-1]])
    return int(row[CHR_END]) >= telo_end - TELOMERE_CLUSTER_THRESHOLD


def compress_paf_telomere_candidates(raw_rows, telo_connect_info, telo_dict,
                                     chr_len):
    """Cluster telomere anchors derived only from the assembly PAF."""
    telo_bounds = telomere_bounds_by_name(telo_dict)
    candidates_by_telo = defaultdict(list)
    for raw_idx, telo_name in telo_connect_info.items():
        candidates_by_telo[telo_name].append(raw_idx)

    compressed_candidates = []
    for telo_name, raw_indices in candidates_by_telo.items():
        nonterminal_indices = []
        terminal_idx = None
        for raw_idx in raw_indices:
            row = raw_rows[raw_idx]
            if is_terminal_telomere_candidate(
                row,
                telo_name,
                telo_bounds,
                chr_len,
            ):
                if terminal_idx is None:
                    terminal_idx = raw_idx
                    continue

                terminal_row = raw_rows[terminal_idx]
                if telo_name.endswith("f"):
                    is_better = (
                        int(row[CHR_STR]) < int(terminal_row[CHR_STR])
                        or (
                            int(row[CHR_STR]) == int(terminal_row[CHR_STR])
                            and int(row[CTG_LEN]) < int(terminal_row[CTG_LEN])
                        )
                    )
                else:
                    is_better = (
                        int(row[CHR_END]) > int(terminal_row[CHR_END])
                        or (
                            int(row[CHR_END]) == int(terminal_row[CHR_END])
                            and int(row[CTG_LEN]) < int(terminal_row[CTG_LEN])
                        )
                    )
                if is_better:
                    terminal_idx = raw_idx
                continue

            if any(
                distance_checker(raw_rows[kept_idx], row)
                < TELOMERE_COMPRESS_RANGE
                for kept_idx in nonterminal_indices
            ):
                continue
            nonterminal_indices.append(raw_idx)

        compressed_candidates.extend(
            (telo_name, raw_idx)
            for raw_idx in nonterminal_indices
        )
        if terminal_idx is not None:
            compressed_candidates.append((telo_name, terminal_idx))

    # break_double_telomere_contig prevents an original contig name from being
    # emitted for more than one final telomere edge.  This PAF-only path keeps
    # the strongest boundary anchor instead of manufacturing split contig paths.
    candidates_by_contig = defaultdict(list)
    for telo_name, raw_idx in compressed_candidates:
        candidates_by_contig[raw_rows[raw_idx][CTG_NAM]].append(
            (telo_name, raw_idx)
        )

    selected_candidates = set()
    duplicate_contigs_removed = 0
    for contig_candidates in candidates_by_contig.values():
        def candidate_rank(candidate):
            telo_name, raw_idx = candidate
            row = raw_rows[raw_idx]
            is_terminal = is_terminal_telomere_candidate(
                row,
                telo_name,
                telo_bounds,
                chr_len,
            )
            return (
                int(is_terminal),
                -int(row[CTG_MAPQ]),
                -(int(row[CHR_END]) - int(row[CHR_STR])),
                raw_idx,
                telo_name,
            )

        selected_candidates.add(min(contig_candidates, key=candidate_rank))
        duplicate_contigs_removed += len(contig_candidates) - 1

    return [
        candidate
        for candidate in compressed_candidates
        if candidate in selected_candidates
    ], duplicate_contigs_removed


def build_paf_telomere_nodes(paf_path, telo_data, repeat_censat_data,
                             chr_len, base_idx):
    """Build telomere nodes from PAF alignments without consulting VCF calls."""
    raw_rows = import_data(paf_path)
    if not raw_rows:
        return [], [], [], Counter()
    excluded_contigs, excluded_rows = find_multi_end_aligned_contigs(raw_rows)

    telo_dict = telo_data_to_dict(telo_data)
    raw_telo_labels = label_node(raw_rows, telo_dict)
    _, preprocess_report, telo_connect_info = preprocess_telo(
        raw_rows,
        raw_telo_labels,
    )
    # preprocess_telo appends a sentinel in-place while scanning contig groups.
    raw_rows.pop()

    boundary_candidate_count = len(telo_connect_info)
    excluded_boundary_candidate_count = sum(
        raw_idx in excluded_rows for raw_idx in telo_connect_info
    )
    telo_connect_info = {
        raw_idx: telo_name
        for raw_idx, telo_name in telo_connect_info.items()
        if raw_idx not in excluded_rows
    }
    if excluded_contigs:
        logging.info(
            f"Excluded {len(excluded_contigs)} multi-end-aligned contigs "
            f"({len(excluded_rows)} PAF rows, "
            f"{excluded_boundary_candidate_count} telomere boundary candidates) "
            "from VCF-mode telomere PAF"
        )
    selected_candidates, duplicate_contigs_removed = \
        compress_paf_telomere_candidates(
            raw_rows,
            telo_connect_info,
            telo_dict,
            chr_len,
        )
    if not selected_candidates:
        return [], [], [], Counter({
            "telomere_paf_rows": len(raw_rows),
            "telomere_paf_contigs": 0,
            "telomere_paf_nodes": 0,
            "telomere_paf_edges": 0,
            "telomere_paf_boundary_candidates": boundary_candidate_count,
            "telomere_paf_duplicate_contigs_removed": duplicate_contigs_removed,
            "telomere_paf_multi_end_contigs_excluded": len(excluded_contigs),
            "telomere_paf_multi_end_rows_excluded": len(excluded_rows),
            "telomere_paf_multi_end_boundary_candidates_excluded": excluded_boundary_candidate_count,
        })

    selected_rows = [raw_rows[raw_idx] for _, raw_idx in selected_candidates]
    selected_censat_labels = label_repeat_node(
        selected_rows,
        repeat_censat_data,
        chr_len,
    )
    nodes = []
    edges = []
    report_rows = []
    for local_idx, (telo_name, raw_idx) in enumerate(selected_candidates):
        row = raw_rows[raw_idx]
        global_node_idx = base_idx + local_idx
        node = [
            row[CTG_NAM],
            int(row[CTG_LEN]),
            int(row[CTG_STR]),
            int(row[CTG_END]),
            row[CTG_DIR],
            row[CHR_NAM],
            int(row[CHR_LEN]),
            int(row[CHR_STR]),
            int(row[CHR_END]),
            int(row[CTG_MAPQ]),
            3,
            global_node_idx,
            global_node_idx,
            raw_telo_labels[raw_idx][0],
            raw_telo_labels[raw_idx][1],
            telo_name,
            "0",
            "0",
            selected_censat_labels[local_idx][1],
            row[CTG_DIR],
            row[CHR_NAM],
            f"1.{raw_idx}",
        ]
        nodes.append(node)
        edges.append((telo_name, (DIR_OUT, global_node_idx, 0)))
        report_rows.append((
            "preprocess_telo_boundary",
            telo_name,
            global_node_idx,
            row[CTG_NAM],
            row[CTG_STR],
            row[CTG_END],
            row[CHR_NAM],
            row[CHR_STR],
            row[CHR_END],
            row[CTG_DIR],
            raw_idx,
        ))

    metrics = Counter({
        "telomere_paf_rows": len(raw_rows),
        "telomere_paf_contigs": len(nodes),
        "telomere_paf_nodes": len(nodes),
        "telomere_paf_edges": len(edges),
        "telomere_paf_boundary_candidates": boundary_candidate_count,
        "telomere_paf_duplicate_contigs_removed": duplicate_contigs_removed,
        "telomere_paf_multi_end_contigs_excluded": len(excluded_contigs),
        "telomere_paf_multi_end_rows_excluded": len(excluded_rows),
        "telomere_paf_multi_end_boundary_candidates_excluded": excluded_boundary_candidate_count,
    })
    metrics.update({
        f"telomere_paf_preprocess_{case.lower()}": len(rows)
        for case, rows in preprocess_report.items()
    })
    return nodes, edges, report_rows, metrics


def annotate_synthetic_nodes(contig_data, telo_data, repeat_censat_data,
                             chr_len, annotate_telomeres=True):
    if not hasattr(telo_data, "items"):
        telo_data = telo_data_to_dict(telo_data)
    if annotate_telomeres:
        telo_labels = label_node(contig_data, telo_data)
    else:
        telo_labels = [("0", "0") for _ in contig_data]
    repeat_labels = no_repeat_labels(contig_data)
    censat_labels = label_repeat_node(contig_data, repeat_censat_data, chr_len)
    for idx, node in enumerate(contig_data):
        node[CTG_TELCHR] = telo_labels[idx][0]
        node[CTG_TELDIR] = telo_labels[idx][1]
        node[CTG_RPTCHR] = repeat_labels[idx][0]
        node[CTG_RPTCASE] = repeat_labels[idx][1]
        node[CTG_CENSAT] = censat_labels[idx][1]
        if (
            annotate_telomeres
            and node[CTG_TELCON] == "0"
            and str(node[CTG_NAM]).startswith(VIRTUAL_TELOMERE_NODE_PREFIX)
        ):
            node[CTG_TELCON] = str(node[CTG_NAM])[len(VIRTUAL_TELOMERE_NODE_PREFIX):]


def has_subtelomeric_telomere_node(contig_data, telo_edges, telo_name,
                                   chrom, chrom_len, side):
    for edge_telo_name, edge in telo_edges:
        if edge_telo_name != telo_name:
            continue
        node_idx = int(edge[1])
        if node_idx < 0 or node_idx >= len(contig_data):
            continue
        node = contig_data[node_idx]
        if node[CHR_NAM] != chrom:
            continue
        if side == "f":
            if int(node[CHR_STR]) <= SUBTELOMERE_LENGTH:
                return True
        elif int(chrom_len) - int(node[CHR_END]) <= SUBTELOMERE_LENGTH:
            return True
    return False


def add_missing_virtual_telomeres(contig_data, telo_edges, chr_len, telo_data,
                                  repeat_censat_data):
    virtual_nodes = []

    for chrom in [f"chr{i}" for i in range(1, 23)] + ["chrX"]:
        if chrom not in chr_len:
            continue
        chrom_len = int(chr_len[chrom])
        for side, path_dir, st, nd in (
            ("f", "+", 0, min(VIRTUAL_TELOMERE_FLANK, chrom_len)),
            ("b", "-", max(0, chrom_len - VIRTUAL_TELOMERE_FLANK), chrom_len),
        ):
            telo_name = f"{chrom}{side}"
            if has_subtelomeric_telomere_node(
                contig_data,
                telo_edges,
                telo_name,
                chrom,
                chrom_len,
                side,
            ):
                continue

            node_idx = len(contig_data) + len(virtual_nodes)
            node = make_synthetic_span_node(
                f"{VIRTUAL_TELOMERE_NODE_PREFIX}{telo_name}",
                chrom,
                st,
                nd,
                path_dir,
                chr_len,
                3,
                node_idx,
                node_idx,
            )
            node[CTG_TELCON] = telo_name
            virtual_nodes.append(node)
            telo_edges.append((telo_name, (DIR_OUT, node_idx, 0)))

    if virtual_nodes:
        annotate_synthetic_nodes(
            virtual_nodes,
            telo_data,
            repeat_censat_data,
            chr_len,
        )
        contig_data.extend(virtual_nodes)
    return len(virtual_nodes)


def add_vcf_nclose_pair(contig_data, nclose_nodes, event_name, chr_a, pos_a, dir_a,
                        chr_b, pos_b, dir_b, chr_len, global_idx_start):
    ctg_typ = 1 if chr_a != chr_b else 2
    s_idx = len(contig_data)
    e_idx = s_idx + 1
    node_a = make_synthetic_vcf_node(
        event_name, chr_a, pos_a, dir_a, "exit", chr_len, ctg_typ, s_idx, global_idx_start
    )
    node_b = make_synthetic_vcf_node(
        event_name, chr_b, pos_b, dir_b, "entry", chr_len, ctg_typ, e_idx, global_idx_start + 1
    )
    node_a[CTG_STRND] = node_b[CTG_STRND] = s_idx
    node_a[CTG_ENDND] = node_b[CTG_ENDND] = e_idx
    contig_data.extend([node_a, node_b])
    nclose_nodes[event_name].append((s_idx, e_idx))
    return global_idx_start + 2


def canonical_nclose_snapshot(nclose_nodes):
    return {
        key: tuple(tuple(int(node_idx) for node_idx in pair) for pair in pair_list)
        for key, pair_list in nclose_nodes.items()
    }


def assert_vcf_nclose_has_no_indel_like(contig_data, nclose_nodes, context):
    """Fail if canonical VCF nclose contains a same-chrom/same-dir pair."""

    offenders = []
    for event_name, pair_list in nclose_nodes.items():
        for pair in pair_list:
            if len(pair) != 2:
                raise AssertionError(
                    f"{context}: invalid canonical VCF nclose pair for "
                    f"{event_name}: {pair}"
                )
            node_a = contig_data[int(pair[0])]
            node_b = contig_data[int(pair[1])]
            if (
                node_a[CHR_NAM] == node_b[CHR_NAM]
                and node_a[CTG_DIR] == node_b[CTG_DIR]
            ):
                offenders.append(
                    (
                        event_name,
                        tuple(int(node_idx) for node_idx in pair),
                        node_a[CHR_NAM],
                        node_a[CTG_DIR],
                    )
                )

    if offenders:
        preview = ", ".join(
            f"{name}:{pair}:{chrom}{direction}/{chrom}{direction}"
            for name, pair, chrom, direction in offenders[:8]
        )
        raise AssertionError(
            f"{context}: canonical VCF nclose contains {len(offenders)} "
            f"Indel-like same-chrom/same-dir pair(s): {preview}"
        )


def add_vcf_type4_graph_pair(contig_data, event, event_index, chr_len,
                             global_idx_start):
    """Add graph-only breakpoint flanks for one non-INS VCF type4 event."""

    chrom = str(event["chrom"])
    st, nd = sorted((int(event["st"]), int(event["nd"])))
    event_type = str(event["event_type"])
    if event_type == "front_jump":
        exit_pos, entry_pos = st, nd
    elif event_type == "back_jump":
        exit_pos, entry_pos = nd, st
    else:
        raise ValueError(f"Invalid VCF type4 graph event_type: {event_type}")

    event_id = str(event.get("event_id", event.get("vcf_id", event_index)))
    contig_name = f"{VCF_TYPE4_GRAPH_NODE_PREFIX}{event_index}_{event_id}"
    s_idx = len(contig_data)
    e_idx = s_idx + 1
    node_a = make_synthetic_vcf_node(
        contig_name, chrom, exit_pos, "+", "exit", chr_len, 4,
        s_idx, global_idx_start,
    )
    node_b = make_synthetic_vcf_node(
        contig_name, chrom, entry_pos, "+", "entry", chr_len, 4,
        e_idx, global_idx_start + 1,
    )
    node_a[CTG_STRND] = node_b[CTG_STRND] = s_idx
    node_a[CTG_ENDND] = node_b[CTG_ENDND] = e_idx

    len_a = int(node_a[CTG_END]) - int(node_a[CTG_STR])
    len_b = int(node_b[CTG_END]) - int(node_b[CTG_STR])
    total_len = len_a + len_b
    node_a[CTG_LEN] = node_b[CTG_LEN] = total_len
    node_a[CTG_STR], node_a[CTG_END] = 0, len_a
    node_b[CTG_STR], node_b[CTG_END] = len_a, total_len
    contig_data.extend([node_a, node_b])
    return global_idx_start + 2


def build_vcf_mode_inputs(context):
    chr_len = find_chr_len(context.reference_fai_path)
    vcf_ins_alt_alignments = load_vcf_ins_alt_alignment_spans(args.alt)
    parsed_vcf = parse_vcf_events(
        args.vcf_input,
        chr_len,
        pass_filters=args.vcf_filter_pass,
        ins_alt_alignments=vcf_ins_alt_alignments,
    )

    depth_df = pd.read_csv(
        context.main_stat_path,
        compression="gzip",
        comment="#",
        sep="\t",
        names=["chr", "st", "nd", "length", "covsite", "totaldepth", "cov", "meandepth"],
    ).query('chr != "chrM"')

    repeat_censat_data = import_censat_repeat_data(context.censat_bed_path)
    telo_data = import_telo_data(context.telomere_bed_path, chr_len)

    summary = dict(parsed_vcf.summary)
    skipped_records = list(parsed_vcf.skipped_records)
    orientation_mismatches = list(parsed_vcf.orientation_mismatches)
    contig_data = []
    nclose_nodes = defaultdict(list)
    # VCF indels are never added to the canonical stage-01 NClose mapping.  The parser
    # places only >=100 kb indels here; step 11 consumes this handoff file.
    step11_vcf_type4_events = list(parsed_vcf.type4_events)
    global_idx = 0

    for spec in parsed_vcf.nclose_specs:
        global_idx = add_vcf_nclose_pair(
            contig_data,
            nclose_nodes,
            spec.event_name,
            spec.chrom_a,
            spec.pos_a,
            spec.dir_a,
            spec.chrom_b,
            spec.pos_b,
            spec.dir_b,
            chr_len,
            global_idx,
        )

    assert_vcf_nclose_has_no_indel_like(
        contig_data,
        nclose_nodes,
        "VCF parser handoff",
    )

    # Prepare every eligible graph-only node in stage 01.  Stage 10 decides
    # whether to activate its edge via --add-indel-graph.
    vcf_type4_graph_events = select_vcf_type4_graph_events(
        step11_vcf_type4_events
    )
    for event_index, event in enumerate(vcf_type4_graph_events):
        global_idx = add_vcf_type4_graph_pair(
            contig_data,
            event,
            event_index,
            chr_len,
            global_idx,
        )
    invalid_graph_svtypes = sorted({
        str(event.get("svtype", "")).upper()
        for event in vcf_type4_graph_events
        if str(event.get("svtype", "")).upper() not in {"DEL", "DUP", "BND"}
    })
    if invalid_graph_svtypes:
        raise AssertionError(
            "VCF graph preparation received ineligible SVTYPE(s): "
            + ", ".join(invalid_graph_svtypes)
        )
    graph_excluded_ins = sum(
        str(event.get("svtype", "")).upper() == "INS"
        for event in step11_vcf_type4_events
    )
    logging.info(
        f"VCF graph preparation: prepared {len(vcf_type4_graph_events)} "
        f"DEL/DUP/BND Indel-like event(s); excluded {graph_excluded_ins} INS"
    )

    # Everything accumulated so far was created from VCF records.  Keep its
    # telomere labels at zero; only PAF/virtual nodes may carry telomere data.
    vcf_node_count = len(contig_data)
    telo_edges = []
    telomere_paf_report_rows = []
    telomere_paf_metrics = Counter()
    if len(context.source.paf_file_paths) > 1:
        # Telomere nodes come exclusively from the assembly PAF.  VCF records
        # above contribute nclose/type4 events but never telomere candidates.
        telomere_paf_nodes, telomere_paf_edges, telomere_paf_report_rows, telomere_paf_metrics = \
            build_paf_telomere_nodes(
                context.source.paf_file_paths[1],
                telo_data,
                repeat_censat_data,
                chr_len,
                len(contig_data),
            )
        contig_data.extend(telomere_paf_nodes)
        telo_edges.extend(telomere_paf_edges)

    virtual_telomere_nodes = add_missing_virtual_telomeres(
        contig_data,
        telo_edges,
        chr_len,
        telo_data,
        repeat_censat_data,
    )
    if virtual_telomere_nodes:
        logging.info(
            f"Added {virtual_telomere_nodes} missing chromosome-end virtual telomere nodes"
        )

    annotate_synthetic_nodes(
        contig_data[:vcf_node_count],
        telo_data,
        repeat_censat_data,
        chr_len,
        annotate_telomeres=False,
    )
    contig_data = [tuple(row) for row in contig_data]

    with open(context.preprocessed_paf_path, "wt") as f:
        for row in contig_data:
            print("\t".join(map(str, row)), file=f)
    with open(context.source.paf_file_paths[0], "wt") as f:
        for row in contig_data:
            print("\t".join(map(str, synthetic_node_to_paf_row(row))), file=f)

    write_telomere_connected_outputs(context.prefix, telo_edges, contig_data)
    with open(f"{context.prefix}/{VCF_TELOMERE_PAF_NODES_TSV}", "wt") as f:
        print("kind\ttelomere\tnode_idx\tcontig\tquery_start\tquery_end\tchrom\tref_start\tref_end\tdir\tpaf_row_idx", file=f)
        for row in telomere_paf_report_rows:
            print("\t".join(map(str, row)), file=f)

    raw_nclose_nodes = nclose_nodes
    all_nclose_comp = defaultdict(list)
    for key, pair_list in nclose_nodes.items():
        for pair in pair_list:
            all_nclose_comp[key].append(tuple(pair))

    with open(f"{context.prefix}/conjoined_type4_ins_del.pkl", "wb") as f:
        pkl.dump(([], []), f)
    with open(f"{context.prefix}/{VCF_TYPE4_EVENTS_PKL}", "wb") as f:
        pkl.dump(step11_vcf_type4_events, f)
    with open(f"{context.prefix}/indel_exclude_idx_set.pkl", "wb") as f:
        pkl.dump(set(), f)

    cen_fragment_meta = find_breakend_centromere(
        repeat_censat_data,
        chr_len,
        depth_df,
        raw_nclose_nodes=raw_nclose_nodes,
        contig_data=contig_data,
        log_context="VCF mode",
    )
    with open(f"{context.prefix}/cen_fragment_data.pkl", "wb") as f:
        pkl.dump(cen_fragment_meta, f)

    chr_corr, chr_rev_corr = chr_correlation_maker(contig_data)
    telo_contig = extract_telomere_connect_contig(
        f"{context.prefix}/telomere_connected_list.txt"
    )
    telo_set = {edge[1] for edge_list in telo_contig.values() for edge in edge_list}
    telo_node_count = len(telo_set)

    write_virtual_ordinary_contig(
        f"{context.prefix}/virtual_ordinary_contig.txt",
        [],
    )
    all_nclose_type = group_nclose_nodes_by_chrom(contig_data, all_nclose_comp)
    uncomp_node_count = write_nclose_nodes_list(
        f"{context.prefix}/all_nclose_nodes_list.txt",
        all_nclose_type,
        contig_data,
    )
    compressed_nclose_type = group_nclose_nodes_by_chrom(contig_data, nclose_nodes)
    write_nclose_nodes_list(
        f"{context.prefix}/compressed_nclose_nodes_list.txt",
        compressed_nclose_type,
        contig_data,
    )
    save_event_catalog(
        context.prefix,
        build_bnd_event_catalog(compressed_nclose_type, contig_data),
    )
    nclose_node_count = write_nclose_nodes_index(
        f"{context.prefix}/nclose_nodes_index.txt",
        nclose_nodes,
        contig_data,
    )

    summary.update({
        "nclose_pairs": sum(len(v) for v in nclose_nodes.values()),
        "synthetic_nodes": len(contig_data),
        "virtual_telomere_nodes": virtual_telomere_nodes,
        "type4_min_span": VCF_TYPE4_MIN_SPAN,
        "vcf_ins_alt_alignment_queries": len(vcf_ins_alt_alignments),
        "vcf_type4_graph_events": len(vcf_type4_graph_events),
        "vcf_type4_graph_nodes": 2 * len(vcf_type4_graph_events),
        "vcf_type4_graph_ins_excluded": graph_excluded_ins,
    })
    summary.update(telomere_paf_metrics)
    with open(f"{context.prefix}/{VCF_MODE_SUMMARY_JSON}", "wt") as f:
        json.dump(dict(summary), f, indent=2, sort_keys=True)
    with open(f"{context.prefix}/{VCF_MODE_SUMMARY_TSV}", "wt") as f:
        print("metric\tvalue", file=f)
        for key, value in sorted(summary.items()):
            print(f"{key}\t{value}", file=f)
    with open(f"{context.prefix}/{VCF_SKIPPED_RECORDS_TSV}", "wt") as f:
        print("line_no\tid\tsvtype\treason\traw", file=f)
        for row in skipped_records:
            print("\t".join(map(str, row)), file=f)
    with open(f"{context.prefix}/{VCF_ORIENTATION_MISMATCH_TSV}", "wt") as f:
        print("line_no\tid\tmate_line_no\tmate_id\tdirs\tmate_dirs\talt\tmate_alt", file=f)
        for row in orientation_mismatches:
            print("\t".join(map(str, row)), file=f)

    no_chrY = False
    if "chrY" in set(depth_df["chr"]):
        ydf = depth_df.query('chr == "chrY"')
        if len(ydf) > 0:
            no_chrY = (len(ydf.query("meandepth != 0")) / len(ydf)) < chrY_MINIMUM_RATIO

    if VCF_TYPE4_MIN_SPAN % M == 0:
        vcf_indel_min_size_label = f"{VCF_TYPE4_MIN_SPAN // M}Mbp"
    elif VCF_TYPE4_MIN_SPAN % K == 0:
        vcf_indel_min_size_label = f"{VCF_TYPE4_MIN_SPAN // K}Kbp"
    else:
        vcf_indel_min_size_label = f"{VCF_TYPE4_MIN_SPAN}bp"

    ins_skip_reasons = []
    if summary["skipped_ins_small_size"] > 0:
        reason = f"SIZE <{vcf_indel_min_size_label}"
        if summary["skipped_ins_small_size"] != summary["skipped_ins"]:
            reason = f"{reason}: {summary['skipped_ins_small_size']}"
        ins_skip_reasons.append(reason)
    if summary["skipped_ins_no_alt_alignment"] > 0:
        reason = "no ALT alignment"
        if summary["skipped_ins_no_alt_alignment"] != summary["skipped_ins"]:
            reason = f"{reason}: {summary['skipped_ins_no_alt_alignment']}"
        ins_skip_reasons.append(reason)
    ins_skip_detail = f" ({', '.join(ins_skip_reasons)})" if ins_skip_reasons else ""

    logging.info(
        f"VCF mode: {summary['used_bnd_events']} BND, {summary['used_inv_events']} INV, "
        f"{summary['used_type4_events']} Indel events "
        f"({summary['used_bnd_type4_events']} BND Indel), "
        f"{summary['skipped_ins']} INS skipped{ins_skip_detail}"
    )

    return NCloseBuildResult(
        df=depth_df,
        no_chrY=no_chrY,
        repeat_censat_data=repeat_censat_data,
        chr_len=chr_len,
        contig_data=contig_data,
        contig_data_size=len(contig_data),
        chr_corr=chr_corr,
        chr_rev_corr=chr_rev_corr,
        telo_contig=telo_contig,
        telo_node_count=telo_node_count,
        telo_set=telo_set,
        rpt_con=set(),
        bnd_contig={key for key in nclose_nodes},
        raw_nclose_nodes=raw_nclose_nodes,
        nclose_nodes=nclose_nodes,
        vctg_dict={},
        all_nclose_comp=all_nclose_comp,
        uncomp_node_count=uncomp_node_count,
        nclose_node_count=nclose_node_count,
        transloc_nclose_pair_count=sum(
            1
            for pair_list in nclose_nodes.values()
            for a, b in pair_list
            if contig_data[a][CHR_NAM] != contig_data[b][CHR_NAM]
        ),
        indel_exclude_idx_set=set(),
        telo_coverage=Counter(),
    )


def plan_initial_nclose_rejections(
    contig_data,
    nclose_nodes,
    repeat_censat_data,
    chr_len,
    depth_df,
):
    """Mark depth-like type2 and centromere-slave candidates for removal."""

    not_using_nclose_node = set()
    type1_nclose_node = []
    type2_nclose_node = defaultdict(list)

    for pair_list in nclose_nodes.values():
        for s, e in pair_list:
            curr_contig_first_fragment = contig_data[s]
            curr_contig_end_fragment = contig_data[e]

            if curr_contig_first_fragment[CHR_NAM] == curr_contig_end_fragment[CHR_NAM]:
                type2_nclose_node[curr_contig_first_fragment[CHR_NAM]].append((s, e))
                if curr_contig_first_fragment[CTG_DIR] == '+':
                    inside_st, inside_nd = sorted([
                        curr_contig_end_fragment[CHR_END],
                        curr_contig_first_fragment[CHR_STR],
                    ])
                else:
                    inside_st, inside_nd = sorted([
                        curr_contig_first_fragment[CHR_END],
                        curr_contig_end_fragment[CHR_STR],
                    ])
                chukji_chrom = curr_contig_first_fragment[CHR_NAM]

                if (
                    not exist_near_bnd_point(depth_df, chukji_chrom, inside_st)
                    and not exist_near_bnd_point(depth_df, chukji_chrom, inside_nd)
                    and not censat_overlap_check(
                        repeat_censat_data,
                        chukji_chrom,
                        inside_st,
                        inside_nd,
                    )
                ):
                    not_using_nclose_node.add((s, e))
            else:
                type1_nclose_node.append((s, e))

    _, centromere_slave = similar_centromere_nclose_cluster(
        nclose_nodes,
        contig_data,
        repeat_censat_data,
        chr_len,
    )

    saved_not_using_nclose_node = set()
    for candidate_pair in not_using_nclose_node:
        chrom = contig_data[candidate_pair[0]][CHR_NAM]
        for type1_pair in type1_nclose_node:
            if (
                contig_data[type1_pair[0]][CHR_NAM] != chrom
                and contig_data[type1_pair[1]][CHR_NAM] != chrom
            ):
                continue
            type2_contig_front = contig_data[candidate_pair[0]]
            type2_contig_back = contig_data[candidate_pair[1]]
            if contig_data[type1_pair[0]][CHR_NAM] == chrom:
                template_contig = contig_data[type1_pair[0]]
                if template_contig[CTG_DIR] == '+':
                    dist = distance_checker(type2_contig_back, template_contig)
                    if (
                        type2_contig_back[CHR_END] < template_contig[CHR_STR]
                        and dist < TYPE2_CONTIG_MINIMUM_LENGTH
                    ) or dist == 0:
                        saved_not_using_nclose_node.add(candidate_pair)
                else:
                    dist = distance_checker(type2_contig_front, template_contig)
                    if (
                        type2_contig_front[CHR_END] > template_contig[CHR_STR]
                        and dist < TYPE2_CONTIG_MINIMUM_LENGTH
                    ) or dist == 0:
                        saved_not_using_nclose_node.add(candidate_pair)
            elif contig_data[type1_pair[1]][CHR_NAM] == chrom:
                template_contig = contig_data[type1_pair[1]]
                if template_contig[CTG_DIR] == '+':
                    dist = distance_checker(type2_contig_front, template_contig)
                    if (
                        type2_contig_front[CHR_STR] > template_contig[CHR_END]
                        and dist < TYPE2_CONTIG_MINIMUM_LENGTH
                    ) or dist == 0:
                        saved_not_using_nclose_node.add(candidate_pair)
                else:
                    dist = distance_checker(type2_contig_back, template_contig)
                    if (
                        type2_contig_back[CHR_STR] < template_contig[CHR_END]
                        and dist < TYPE2_CONTIG_MINIMUM_LENGTH
                    ) or dist == 0:
                        saved_not_using_nclose_node.add(candidate_pair)

    not_using_nclose_node |= set(centromere_slave.keys())
    return (
        not_using_nclose_node,
        saved_not_using_nclose_node,
        type2_nclose_node,
    )




def apply_initial_nclose_rejections(candidates, rejected_pairs, rescued_pairs):
    """Apply the planned initial rejection with its historical rescue override."""

    return apply_nclose_filter(
        candidates,
        "initial_type2_and_centromere",
        lambda candidate: (
            "FILTERED_02_INITIAL"
            if candidate.path_pair in rejected_pairs
            and candidate.path_pair not in rescued_pairs
            else None
        ),
    )


def _censat_at_chromosome_end(
    contig_data,
    chr_len,
    repeat_censat_data,
    node_idx,
):
    chrom = contig_data[node_idx][CHR_NAM]
    chromosome_length = chr_len.get(chrom, 0)
    chunk_start = contig_data[node_idx][CHR_STR]
    chunk_end = contig_data[node_idx][CHR_END]
    for interval_start, interval_end in repeat_censat_data.get(chrom, []):
        if max(interval_start, chunk_start) < min(interval_end, chunk_end):
            if interval_start == 0 or interval_end == chromosome_length:
                return True
    return False


def _canonical_censat_pair_info(contig_data, pair):
    node_a_idx, node_b_idx = pair
    node_a = contig_data[node_a_idx]
    node_b = contig_data[node_b_idx]
    key_a = (
        chr2int(node_a[CHR_NAM]),
        node_a[CHR_STR],
        node_a[CHR_END],
    )
    key_b = (
        chr2int(node_b[CHR_NAM]),
        node_b[CHR_STR],
        node_b[CHR_END],
    )
    if key_b < key_a:
        node_a_idx, node_b_idx = node_b_idx, node_a_idx
        node_a, node_b = node_b, node_a

    is_for = node_a_idx < node_b_idx
    return {
        "chroms": (node_a[CHR_NAM], node_b[CHR_NAM]),
        "idxs": (node_a_idx, node_b_idx),
        "dirs": (
            get_corr_dir(is_for, node_a[CTG_DIR]),
            get_corr_dir(is_for, node_b[CTG_DIR]),
        ),
    }


def _normalized_censat_endpoint_dirs(contig_data, pair):
    endpoint_dirs = []
    for order_idx, node_idx in enumerate(pair):
        if not node_is_censat(contig_data[node_idx]):
            continue
        ctg_dir = contig_data[node_idx][CTG_DIR]
        normalized_dir = ctg_dir if order_idx == 0 else _flip_ctg_dir(ctg_dir)
        endpoint_dirs.append((node_idx, normalized_dir))
    return endpoint_dirs


def apply_censat_censat_filter(
    candidates,
    contig_data,
    chr_len,
    repeat_censat_data,
):
    """Apply terminal, direction, competitor, then MAPQ CENSAT rules."""

    mapq60_groups = defaultdict(list)
    for candidate in candidates:
        pair = candidate.path_pair
        if classify_censat_pair(contig_data, pair) != CensatPairClass.BOTH:
            continue
        if (
            contig_data[pair[0]][CTG_MAPQ] < 60
            or contig_data[pair[1]][CTG_MAPQ] < 60
        ):
            continue
        pair_info = _canonical_censat_pair_info(contig_data, pair)
        mapq60_groups[pair_info["chroms"]].append((candidate, pair_info))

    opposite_pairs = set()
    for items in mapq60_groups.values():
        for i in range(len(items)):
            candidate_i, info_i = items[i]
            node_a_i, node_b_i = info_i["idxs"]
            for j in range(i + 1, len(items)):
                candidate_j, info_j = items[j]
                if not (
                    info_i["dirs"][0] != info_j["dirs"][0]
                    and info_i["dirs"][1] != info_j["dirs"][1]
                ):
                    continue
                node_a_j, node_b_j = info_j["idxs"]
                if (
                    distance_checker(
                        contig_data[node_a_i],
                        contig_data[node_a_j],
                    ) <= NCLOSE_COMPRESS_LIMIT
                    and distance_checker(
                        contig_data[node_b_i],
                        contig_data[node_b_j],
                    ) <= NCLOSE_COMPRESS_LIMIT
                ):
                    opposite_pairs.add(candidate_i.identity)
                    opposite_pairs.add(candidate_j.identity)

    def reject_reason(candidate):
        pair = candidate.path_pair
        if classify_censat_pair(contig_data, pair) != CensatPairClass.BOTH:
            return None
        if (
            _censat_at_chromosome_end(
                contig_data,
                chr_len,
                repeat_censat_data,
                pair[0],
            )
            or _censat_at_chromosome_end(
                contig_data,
                chr_len,
                repeat_censat_data,
                pair[1],
            )
        ):
            return "terminal"
        if (
            contig_data[pair[0]][CHR_NAM]
            == contig_data[pair[1]][CHR_NAM]
            and contig_data[pair[0]][CTG_DIR]
            != contig_data[pair[1]][CTG_DIR]
        ):
            return "same_chrom_opposite_dir"
        if candidate.identity in opposite_pairs:
            return "opposite_dir"
        if (
            contig_data[pair[0]][CTG_MAPQ] < 60
            or contig_data[pair[1]][CTG_MAPQ] < 60
        ):
            return "mapq"
        return None

    kept, rejections = apply_nclose_filter(
        candidates,
        "censat_censat",
        reject_reason,
    )
    removal_counts = Counter(rejection.reason for rejection in rejections)
    logging.info(
        f'Removed {removal_counts["mapq"]} censat-censat nclose where either endpoint MAPQ < 60, '
        f'{removal_counts["terminal"]} censat-censat nclose with a terminal-censat endpoint, '
        f'{removal_counts["same_chrom_opposite_dir"]} same-chromosome opposite-direction censat-censat nclose, '
        f'{removal_counts["opposite_dir"]} opposite-direction censat-censat nclose'
    )
    return kept, rejections


def apply_censat_fragment_direction_filter(
    candidates,
    contig_data,
    cen_fragment_meta,
):
    """Keep two-CENSAT pairs only when every inferred CEN direction matches."""

    cent_fragment_chroms = set(cen_fragment_meta.keys())

    def reject_reason(candidate):
        pair = candidate.path_pair
        if classify_censat_pair(contig_data, pair) != CensatPairClass.BOTH:
            return None
        for node_idx, normalized_dir in _normalized_censat_endpoint_dirs(
            contig_data,
            pair,
        ):
            chrom = contig_data[node_idx][CHR_NAM]
            if chrom not in cent_fragment_chroms:
                continue
            if normalized_dir != _cen_fragment_target_dir_from_meta(
                cen_fragment_meta,
                chrom,
            ):
                return "direction_mismatch"
        return None

    kept, rejections = apply_nclose_filter(
        candidates,
        "censat_fragment_direction",
        reject_reason,
    )
    logging.info(
        f'Removed {len(rejections)} censat-censat '
        f'nclose pairs with cen_fragment direction mismatch'
    )
    return kept, rejections


def apply_simple_alt_preference_filter(
    candidates,
    contig_data,
    cen_fragment_meta,
):
    """Prefer simple-alt evidence within strictly overlapping target loci."""

    cent_fragment_chroms = set(cen_fragment_meta.keys())
    candidate_groups = defaultdict(list)
    for candidate in candidates:
        pair = candidate.path_pair
        one_censat = build_censat_noncensat_candidate(
            contig_data,
            candidate.contig_name,
            pair,
        )
        if one_censat is None:
            continue
        if one_censat.censat_chrom in cent_fragment_chroms:
            continue
        noncensat_node = contig_data[one_censat.noncensat_idx]
        candidate_groups[
            (one_censat.censat_chrom, one_censat.noncensat_chrom)
        ].append((
            candidate,
            noncensat_node[CHR_STR],
            noncensat_node[CHR_END],
            one_censat.is_simple_alt,
        ))

    identities_to_remove = set()
    for items in candidate_groups.values():
        items.sort(key=lambda item: (item[1], item[2]))
        i = 0
        while i < len(items):
            j = i + 1
            group_end = items[i][2]
            while j < len(items) and items[j][1] < group_end:
                group_end = max(group_end, items[j][2])
                j += 1
            overlapping_items = items[i:j]
            if any(item[3] for item in overlapping_items):
                identities_to_remove.update(
                    item[0].identity
                    for item in overlapping_items
                    if not item[3]
                )
            i = j

    kept, rejections = apply_nclose_filter(
        candidates,
        "simple_alt_preference",
        lambda candidate: (
            "overlapping_simple_alt"
            if candidate.identity in identities_to_remove
            else None
        ),
    )
    logging.info(
        f'Removed {len(rejections)} non-simple censat-noncensat nclose pairs '
        f'overlapping a simple_ctg_alt non-censat locus'
    )
    return kept, rejections


def apply_subtelomeric_orientation_filter(
    candidates,
    contig_data,
    telo_contig,
    chr_len,
):
    """Remove pairs whose two endpoints have the same telomere orientation."""

    telo_anchor_by_chrom = defaultdict(list)
    telo_anchor_seen = set()

    def add_telo_anchor(telo_name, node_idx):
        if (
            node_idx < 0
            or node_idx >= len(contig_data)
            or telo_name == '0'
        ):
            return
        key = (telo_name, node_idx)
        if key in telo_anchor_seen:
            return
        telo_anchor_seen.add(key)
        telo_anchor_by_chrom[contig_data[node_idx][CHR_NAM]].append(
            (telo_name, node_idx)
        )

    for telo_name, edge_list in telo_contig.items():
        for edge in edge_list:
            add_telo_anchor(telo_name, edge[1])
    for node_idx, chunk in enumerate(contig_data):
        if chunk[CTG_TELCON] != '0':
            add_telo_anchor(chunk[CTG_TELCON], node_idx)

    def orientation_from_telo_side(telo_side, ctg_dir):
        if telo_side == 'f':
            return 'inward' if ctg_dir == '+' else 'outward'
        return 'outward' if ctg_dir == '+' else 'inward'

    def endpoint_orientations(node_idx):
        chrom = contig_data[node_idx][CHR_NAM]
        chromosome_length = chr_len.get(chrom, 0)
        ctg_dir = contig_data[node_idx][CTG_DIR]
        orientations = set()
        if contig_data[node_idx][CHR_STR] < SUBTELO_TIP_LIMIT:
            orientations.add(orientation_from_telo_side('f', ctg_dir))
        if contig_data[node_idx][CHR_END] > chromosome_length - SUBTELO_TIP_LIMIT:
            orientations.add(orientation_from_telo_side('b', ctg_dir))
        for telo_name, anchor_idx in telo_anchor_by_chrom.get(chrom, []):
            if (
                distance_checker(
                    contig_data[node_idx],
                    contig_data[anchor_idx],
                ) < SUBTELO_TIP_LIMIT
            ):
                orientations.add(
                    orientation_from_telo_side(telo_name[-1], ctg_dir)
                )
        return orientations

    def reject_reason(candidate):
        start_orientations = endpoint_orientations(candidate.path_pair[0])
        end_orientations = endpoint_orientations(candidate.path_pair[1])
        if start_orientations and not start_orientations.isdisjoint(
            end_orientations
        ):
            return "same_orientation"
        return None

    kept, rejections = apply_nclose_filter(
        candidates,
        "subtelomeric_orientation",
        reject_reason,
    )
    logging.info(
        f"Removed {len(rejections)} same-orientation subtelomeric tip nclose "
        f"(both junctions within {SUBTELO_TIP_LIMIT//1000}kb of chromosome end or telomere anchor, "
        f"same telo-in/out direction)"
    )
    return kept, rejections


def apply_nclose_count_vaf_filter(contig_data, candidates):
    """Optionally remove candidates without sufficient raw-read VAF support."""

    if not args.check_nclose_count:
        return candidates
    if SKIP_BAM_ANAL:
        logging.warning("NClose raw-count VAF filter requested but BAM analysis is skipped")
        return candidates

    nclose_count_candidates = build_nclose_count_candidates(
        contig_data,
        candidates,
    )
    with open(f"{PREFIX}/{NCLOSE_COUNT_CANDIDATE_PKL}", "wb") as f:
        pkl.dump(nclose_count_candidates, f)
    logging.info(
        f"NClose raw-count VAF filter candidates : {len(nclose_count_candidates)} "
        f"(threshold={args.nclose_count_vaf_threshold})"
    )
    if not nclose_count_candidates:
        return candidates

    thread_lim = min(JULIA_BAM_THREAD_LIM, THREAD)
    progress_args = ['--progress'] if args.progress else []
    nclose_count_result = subprocess.run(
        [
            'python',
            "-X",
            f"juliacall-threads={thread_lim}",
            "-X",
            "juliacall-handle-signals=yes",
            os.path.join(
                os.path.dirname(os.path.abspath(__file__)),
                '03_Anal_bam.py',
            ),
            PREFIX,
            read_bam_loc,
            CHROMOSOME_INFO_FILE_PATH,
            main_stat_loc,
            '--nclose_count_only',
            '--nclose_count_vaf_threshold',
            str(args.nclose_count_vaf_threshold),
        ] + progress_args
    )
    if nclose_count_result.returncode != 0:
        logging.warning(
            f"NClose raw-count BAM analysis failed with exit code {nclose_count_result.returncode}; "
            "skip nclose raw-count VAF filter"
        )
        return candidates

    keep_pairs, record_count, keep_count = collect_nclose_count_keep_pairs(
        PREFIX,
        args.nclose_count_vaf_threshold,
    )
    if keep_pairs is None:
        return candidates
    kept, rejections = apply_nclose_filter(
        candidates,
        "raw_count_vaf",
        lambda candidate: (
            None if candidate.event_key in keep_pairs else "below_threshold"
        ),
    )
    logging.info(
        f"NClose raw-count VAF filter : kept {keep_count}/{record_count} "
        f"candidate nclose pairs, removed {len(rejections)}"
    )
    return kept


def write_raw_translocation_candidate_artifact(
    contig_data,
    candidates,
    raw_nclose_nodes,
    all_nclose_comp,
    chr_len,
):
    """Serialize BAM-analysis inputs before virtual-inversion/user filtering."""

    base_candidates = build_raw_translocation_candidates(
        contig_data,
        candidates,
        chr_len,
    )
    all_raw_nodes = convert_all_nclose_comp_to_nclose_nodes(
        contig_data,
        all_nclose_comp,
    )
    raw_candidate_nodes = build_ecdna_nclose_nodes(
        raw_nclose_nodes,
        all_raw_nodes,
    )
    logging.info(
        f"Raw-read translocation candidate input : "
        f"{sum(len(v) for v in raw_nclose_nodes.values())} raw nclose, "
        f"{sum(len(v) for v in all_raw_nodes.values())} all nclose, "
        f"{sum(len(v) for v in raw_candidate_nodes.values())} merged nclose"
    )
    large_same_chrom_candidates = build_raw_translocation_candidates(
        contig_data,
        raw_candidate_nodes,
        chr_len,
        candidate_filter=is_large_same_chrom_raw_candidate,
    )
    raw_translocation_candidates = merge_raw_translocation_candidates(
        base_candidates,
        large_same_chrom_candidates,
    )
    with open(f"{PREFIX}/{RAW_TRANSLOCATION_CANDIDATE_PKL}", "wb") as f:
        pkl.dump(raw_translocation_candidates, f)
    logging.info(
        f"Raw-read translocation candidates : {len(raw_translocation_candidates)} "
        f"({len(base_candidates)} final-nclose, "
        f"{len(large_same_chrom_candidates)} large same-chrom)"
    )


def apply_raw_virtual_inversion_filter(candidates):
    """Run the raw-read analysis and remove its virtual-inversion pairs."""

    if SKIP_BAM_ANAL:
        logging.info("Raw-read translocation BAM analysis skipped")
        raw_virtual_inv_nclose_pairs = set()
    else:
        thread_lim = min(JULIA_BAM_THREAD_LIM, THREAD)
        progress_args = ['--progress'] if args.progress else []
        raw_bam_result = subprocess.run(
            [
                'python',
                "-X",
                f"juliacall-threads={thread_lim}",
                "-X",
                "juliacall-handle-signals=yes",
                os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    '03_Anal_bam.py',
                ),
                PREFIX,
                read_bam_loc,
                CHROMOSOME_INFO_FILE_PATH,
                main_stat_loc,
                '--skip_nclose_count',
            ] + progress_args
        )
        if raw_bam_result.returncode == 0:
            raw_virtual_inv_nclose_pairs = \
                collect_raw_virtual_inv_nclose_pairs(PREFIX)
        else:
            logging.warning(
                f"Raw-read translocation BAM analysis failed with exit code {raw_bam_result.returncode}; "
                "skip raw virtual-inversion nclose removal"
            )
            raw_virtual_inv_nclose_pairs = set()

    if not raw_virtual_inv_nclose_pairs:
        return candidates
    kept, rejections = apply_nclose_filter(
        candidates,
        "raw_virtual_inversion",
        lambda candidate: (
            "balanced_reference_spans"
            if candidate.event_key in raw_virtual_inv_nclose_pairs
            else None
        ),
    )
    logging.info(
        f"Removed {len(rejections)} raw-read virtual-inversion nclose pairs "
        f"({len(raw_virtual_inv_nclose_pairs)} unique pairs) from final compressed nclose"
    )
    return kept


def apply_user_nclose_exclusions(candidates):
    """Apply whole-event and step-11 INDEL exclusions after BAM analysis."""

    indel_exclude_idx_set = set()
    excluded_owner_names = set()
    if args.exclude_nclose_list_loc is not None:
        active_owner_names = {
            candidate.contig_name for candidate in candidates
        }
        with open(args.exclude_nclose_list_loc, 'r') as f:
            excluded_names = []
            for line in f:
                name = line.strip()
                if name.startswith('INDEL_INDEX_'):
                    indel_exclude_idx_set.add(
                        int(name.removeprefix('INDEL_INDEX_'))
                    )
                elif name in active_owner_names:
                    active_owner_names.remove(name)
                    excluded_owner_names.add(name)
                    excluded_names.append(name)

        if excluded_names:
            logging.warning(f"Skipped contig : {', '.join(excluded_names)}")
        if indel_exclude_idx_set:
            logging.warning(
                "Skipped indel index : "
                + ', '.join(map(str, indel_exclude_idx_set))
            )

    with open(f'{PREFIX}/indel_exclude_idx_set.pkl', 'wb') as f:
        pkl.dump(indel_exclude_idx_set, f)
    kept, _ = apply_nclose_filter(
        candidates,
        "user_exclusion",
        lambda candidate: (
            "excluded_owner"
            if candidate.contig_name in excluded_owner_names
            else None
        ),
    )
    return kept, indel_exclude_idx_set


def write_compressed_nclose_outputs(
    nclose_nodes,
    contig_data,
    not_using_nclose_node,
    saved_not_using_nclose_node,
    rpt_con,
):
    """Write the compressed list/catalog and return translocation pair count."""

    nclose_type = group_nclose_nodes_by_chrom(contig_data, nclose_nodes)
    translocation_pair_count = 0
    for chrom_pair, pair_list in nclose_type.items():
        if chrom_pair[0] != chrom_pair[1]:
            translocation_pair_count += len(pair_list)

        for pair in pair_list:
            original_nclose = tuple(sorted(pair))
            assert (
                original_nclose not in not_using_nclose_node
                or original_nclose in saved_not_using_nclose_node
            )

    write_nclose_nodes_list(
        f"{PREFIX}/compressed_nclose_nodes_list.txt",
        nclose_type,
        contig_data,
        rpt_con,
    )
    save_event_catalog(
        PREFIX,
        build_bnd_event_catalog(nclose_type, contig_data),
    )
    return translocation_pair_count


def finalize_nclose_outputs(
    nclose_nodes,
    contig_data,
    not_using_nclose_node,
    saved_not_using_nclose_node,
    rpt_con,
    uncomp_node_count,
    telo_node_count,
):
    """Finalize graph-facing NClose artifacts."""

    translocation_pair_count = write_compressed_nclose_outputs(
        nclose_nodes,
        contig_data,
        not_using_nclose_node,
        saved_not_using_nclose_node,
        rpt_con,
    )
    node_count = write_nclose_nodes_index(
        f"{PREFIX}/nclose_nodes_index.txt",
        nclose_nodes,
        contig_data,
    )
    logging.info(f"Uncompressed NClose node count : {uncomp_node_count}")
    logging.info(f"NClose node count : {node_count}")
    logging.info(f"Telomere connected node count : {telo_node_count}")
    return translocation_pair_count, node_count


def _store_nclose_filter_result(state, result, stage_name):
    try:
        kept, rejections = result
    except (TypeError, ValueError) as exc:
        raise TypeError(
            f"NClose filter stage {stage_name!r} must return "
            "(kept_candidates, rejections)"
        ) from exc
    kept = list(kept)
    rejections = list(rejections)
    if not all(isinstance(candidate, NCloseCandidate) for candidate in kept):
        raise TypeError(
            f"NClose filter stage {stage_name!r} returned a non-NCloseCandidate"
        )
    if not all(
        isinstance(rejection, NCloseRejection)
        for rejection in rejections
    ):
        raise TypeError(
            f"NClose filter stage {stage_name!r} returned a non-NCloseRejection"
        )
    rejections = [
        rejection
        if rejection.stage == stage_name
        else NCloseRejection(
            rejection.candidate,
            stage_name,
            rejection.reason,
        )
        for rejection in rejections
    ]
    state.candidates = kept
    state.rejections.extend(rejections)
    return state


def make_nclose_filter_stage(name, filter_candidates):
    """Adapt a small candidate filter to the named NClose stage contract.

    ``filter_candidates`` receives ``(NClosePipelineContext, candidates)`` and
    returns ``(kept_candidates, rejections)``.  This is the intended API for
    adding another filter without editing :func:`nclose_calc`.
    """

    def run(context, state):
        return _store_nclose_filter_result(
            state,
            filter_candidates(context, tuple(state.candidates)),
            name,
        )

    return NClosePipelineStage(name=name, run=run)


def _run_initial_rejection_stage(context, state):
    return _store_nclose_filter_result(
        state,
        apply_initial_nclose_rejections(
            state.candidates,
            context.rejected_pairs,
            context.rescued_pairs,
        ),
        "initial_rejections",
    )


def _run_censat_pair_stage(context, state):
    return _store_nclose_filter_result(
        state,
        apply_censat_censat_filter(
            state.candidates,
            context.contig_data,
            context.chr_len,
            context.repeat_censat_data,
        ),
        "censat_pair",
    )


def _run_censat_fragment_direction_stage(context, state):
    return _store_nclose_filter_result(
        state,
        apply_censat_fragment_direction_filter(
            state.candidates,
            context.contig_data,
            context.cen_fragment_meta,
        ),
        "censat_fragment_direction",
    )


def _run_simple_alt_preference_stage(context, state):
    return _store_nclose_filter_result(
        state,
        apply_simple_alt_preference_filter(
            state.candidates,
            context.contig_data,
            context.cen_fragment_meta,
        ),
        "simple_alt_preference",
    )


def _run_combined_censat_noncensat_stage(context, state):
    missing_groups = collect_missing_cen_fragment_dir_censat_noncensat(
        context.contig_data,
        state.candidates,
        context.cen_fragment_meta,
    )
    source_candidates = extract_source_censat_type2_candidates(
        context.source,
        context.repeat_censat_data,
    )
    combined = add_nearest_combined_censat_noncensat_ncloses(
        context.contig_data,
        missing_groups,
        source_candidates,
        context.repeat_censat_data,
        context.chr_len,
        context.cen_fragment_meta,
        context.build.preprocessed_paf_path,
    )
    state.candidates.extend(combined)
    state.stage_metadata["combined_censat_noncensat"] = {
        "added": len(combined),
        "source_candidates": len(source_candidates),
        "missing_groups": len(missing_groups),
    }
    logging.info(
        "Added %d nearest combined censat-noncensat nclose pairs from %d "
        "raw censat-internal type2 candidates across %d target-missing groups",
        len(combined),
        len(source_candidates),
        len(missing_groups),
    )
    if combined:
        state.contig_data_size = len(context.contig_data)
        state.chr_corr, state.chr_rev_corr = chr_correlation_maker(
            context.contig_data
        )
    return state


def _run_subtelomeric_orientation_stage(context, state):
    return _store_nclose_filter_result(
        state,
        apply_subtelomeric_orientation_filter(
            state.candidates,
            context.contig_data,
            context.telo_contig,
            context.chr_len,
        ),
        "subtelomeric_orientation",
    )


def _run_offset_direction_stage(context, state):
    result = apply_offset_direction_mismatched_censat_noncensat_filter(
        state.candidates,
        context.contig_data,
        context.cen_fragment_meta,
    )
    state = _store_nclose_filter_result(state, result, "offset_direction")
    logging.info(
        "Removed %d offset-direction-mismatched censat-noncensat nclose pairs",
        len(result[1]),
    )
    return state


def _run_raw_count_vaf_stage(context, state):
    before = list(state.candidates)
    kept = list(
        apply_nclose_count_vaf_filter(context.contig_data, before)
    )
    kept_identities = {candidate.identity for candidate in kept}
    state.rejections.extend(
        NCloseRejection(candidate, "raw_count_vaf", "below_threshold")
        for candidate in before
        if candidate.identity not in kept_identities
    )
    state.candidates = kept
    return state


def _run_raw_translocation_artifact_stage(context, state):
    write_raw_translocation_candidate_artifact(
        context.contig_data,
        state.candidates,
        context.raw_nclose_nodes,
        context.all_nclose_nodes,
        context.chr_len,
    )
    return state


def _run_raw_virtual_inversion_stage(context, state):
    before = list(state.candidates)
    kept = list(apply_raw_virtual_inversion_filter(before))
    kept_identities = {candidate.identity for candidate in kept}
    state.rejections.extend(
        NCloseRejection(
            candidate,
            "raw_virtual_inversion",
            "balanced_reference_spans",
        )
        for candidate in before
        if candidate.identity not in kept_identities
    )
    state.candidates = kept
    return state


def _run_user_exclusion_stage(context, state):
    before = list(state.candidates)
    kept, state.indel_exclude_idx_set = apply_user_nclose_exclusions(before)
    state.candidates = list(kept)
    kept_identities = {
        candidate.identity for candidate in state.candidates
    }
    state.rejections.extend(
        NCloseRejection(candidate, "user_exclusion", "excluded_owner")
        for candidate in before
        if candidate.identity not in kept_identities
    )
    return state


def default_nclose_pipeline_stages():
    """Return the 02-compatible NClose postprocessing sequence."""

    return (
        NClosePipelineStage("initial_rejections", _run_initial_rejection_stage),
        NClosePipelineStage("censat_pair", _run_censat_pair_stage),
        NClosePipelineStage(
            "censat_fragment_direction",
            _run_censat_fragment_direction_stage,
        ),
        NClosePipelineStage(
            "simple_alt_preference",
            _run_simple_alt_preference_stage,
        ),
        NClosePipelineStage(
            "combined_censat_noncensat",
            _run_combined_censat_noncensat_stage,
        ),
        NClosePipelineStage(
            "subtelomeric_orientation",
            _run_subtelomeric_orientation_stage,
        ),
        NClosePipelineStage("offset_direction", _run_offset_direction_stage),
        NClosePipelineStage("raw_count_vaf", _run_raw_count_vaf_stage),
        NClosePipelineStage(
            "raw_translocation_artifact",
            _run_raw_translocation_artifact_stage,
        ),
        NClosePipelineStage(
            "raw_virtual_inversion",
            _run_raw_virtual_inversion_stage,
        ),
        NClosePipelineStage("user_exclusion", _run_user_exclusion_stage),
    )


def run_nclose_pipeline(context, state, stages):
    """Run ordered NClose stages and retain per-stage rejection provenance."""

    names = [stage.name for stage in stages]
    if len(names) != len(set(names)):
        raise ValueError(f"Duplicate NClose pipeline stages: {names}")
    for stage in stages:
        before_count = len(state.candidates)
        before_rejections = len(state.rejections)
        updated = stage.run(context, state)
        if updated is not None:
            if not isinstance(updated, NClosePipelineState):
                raise TypeError(
                    f"NClose stage {stage.name!r} must return "
                    "NClosePipelineState or None"
                )
            state = updated
        if not all(
            isinstance(candidate, NCloseCandidate)
            for candidate in state.candidates
        ):
            raise TypeError(
                f"NClose stage {stage.name!r} returned a non-NCloseCandidate"
            )
        if not all(
            isinstance(rejection, NCloseRejection)
            for rejection in state.rejections
        ):
            raise TypeError(
                f"NClose stage {stage.name!r} recorded a non-NCloseRejection"
            )
        new_rejections = len(state.rejections) - before_rejections
        after_count = len(state.candidates)
        state.stage_records.append({
            "name": stage.name,
            "before_candidates": before_count,
            "after_candidates": after_count,
            "added_candidates": max(after_count - before_count, 0),
            "removed_candidates": max(before_count - after_count, 0),
            "recorded_rejections": new_rejections,
            "metadata": dict(state.stage_metadata.get(stage.name, {})),
        })
        logging.info(
            "NClose stage %s: %d -> %d candidate(s)",
            stage.name,
            before_count,
            after_count,
        )
    return state


def default_stage01_pipeline():
    """Build the default, side-effect-free stage composition object."""

    return Stage01Pipeline(
        candidate_builder=build_unitig_nclose_candidates,
        contig_stages=default_contig_pipeline_stages(),
        nclose_stages=default_nclose_pipeline_stages(),
    )


def nclose_calc(
    context,
    preprocess_result,
    candidate_builder=build_unitig_nclose_candidates,
    *,
    nclose_stages=None,
):
    source = context.source
    repeat_censat_data = import_censat_repeat_data(context.censat_bed_path)
    chr_len = find_chr_len(context.reference_fai_path)
    contig_data = import_data2(context.preprocessed_paf_path)

    TELO_CONNECT_NODES_INFO_PATH = context.prefix+"/telomere_connected_list.txt"

    os.makedirs(context.prefix, exist_ok=True)

    contig_data_size = len(contig_data)

    chr_corr, chr_rev_corr = chr_correlation_maker(contig_data)

    telo_contig = extract_telomere_connect_contig(TELO_CONNECT_NODES_INFO_PATH)

    telo_node_count = 0
    telo_set = set()

    for i in telo_contig:
        for j in telo_contig[i]:
            telo_node_count+=1
            telo_set.add(j[1])

    repeat_data = import_repeat_data_00(context.repeat_bed_path)
    rpt_con = extract_all_repeat_contig(contig_data, repeat_data, CTG_RPTCASE, NON_REPEAT_NOISE_RATIO)

    bnd_contig = extract_bnd_contig(contig_data)


    # Type-1/2 candidates from the effective unitig source and its children.
    candidate_build = candidate_builder(NCloseCandidateBuildContext(
        contig_data=contig_data,
        bnd_contig=bnd_contig,
        repeat_contig_names=rpt_con,
        repeat_censat_data=repeat_censat_data,
        paf_file_paths=source.paf_file_paths,
        original_paf_paths=source.original_paf_paths,
        telo_set=telo_set,
        telo_contig=telo_contig,
        chr_len=chr_len,
        original_contig_names=context.ori_ctg_name_data,
    ))
    if not isinstance(candidate_build, NCloseCandidateBuildResult):
        raise TypeError(
            "Stage-01 candidate builder must return "
            "NCloseCandidateBuildResult"
        )
    if not all(
        isinstance(candidate, NCloseCandidate)
        for candidate in candidate_build.candidates
    ):
        raise TypeError(
            "Stage-01 candidate builder returned a non-NCloseCandidate"
        )
    raw_nclose_nodes = candidate_build.raw_nodes
    vctg_dict = candidate_build.virtual_contigs
    all_nclose_comp = candidate_build.all_nodes

    depth_df = preprocess_result.depth_df
    cen_fragment_meta = find_breakend_centromere(
        repeat_censat_data,
        chr_len,
        depth_df,
        raw_nclose_nodes=raw_nclose_nodes,
        contig_data=contig_data,
        log_context="Raw-nclose adjusted",
    )
    with open(f'{context.prefix}/cen_fragment_data.pkl', 'wb') as f:
        pkl.dump(cen_fragment_meta, f)

    (
        not_using_nclose_node,
        saved_not_using_nclose_node,
        type2_nclose_node,
    ) = plan_initial_nclose_rejections(
        contig_data,
        raw_nclose_nodes,
        repeat_censat_data,
        chr_len,
        depth_df,
    )

    with open(f"{context.prefix}/conjoined_type4_ins_del.pkl", "wb") as f:
        pkl.dump(conjoined_type4(contig_data, type2_nclose_node), f)

    virtual_ordinary_contig = make_virtual_ord_ctg(contig_data, vctg_dict)
    write_virtual_ordinary_contig(
        f"{context.prefix}/virtual_ordinary_contig.txt",
        virtual_ordinary_contig,
    )

    all_nclose_type = group_nclose_nodes_by_chrom(contig_data, all_nclose_comp)
    uncomp_node_count = write_nclose_nodes_list(
        f"{context.prefix}/all_nclose_nodes_list.txt",
        all_nclose_type,
        contig_data,
        rpt_con,
    )

    nclose_pipeline_context = NClosePipelineContext(
        build=context,
        source=source,
        contig_data=contig_data,
        repeat_censat_data=repeat_censat_data,
        chr_len=chr_len,
        cen_fragment_meta=cen_fragment_meta,
        telo_contig=telo_contig,
        raw_nclose_nodes=raw_nclose_nodes,
        all_nclose_nodes=all_nclose_comp,
        rejected_pairs=not_using_nclose_node,
        rescued_pairs=saved_not_using_nclose_node,
    )
    nclose_pipeline_state = NClosePipelineState(
        candidates=list(candidate_build.candidates),
        contig_data_size=contig_data_size,
        chr_corr=chr_corr,
        chr_rev_corr=chr_rev_corr,
    )
    nclose_pipeline_state = run_nclose_pipeline(
        nclose_pipeline_context,
        nclose_pipeline_state,
        (
            default_nclose_pipeline_stages()
            if nclose_stages is None
            else tuple(nclose_stages)
        ),
    )
    nclose_candidates = nclose_pipeline_state.candidates
    contig_data_size = nclose_pipeline_state.contig_data_size
    chr_corr = nclose_pipeline_state.chr_corr
    chr_rev_corr = nclose_pipeline_state.chr_rev_corr
    indel_exclude_idx_set = nclose_pipeline_state.indel_exclude_idx_set

    nclose_nodes = candidates_to_legacy(nclose_candidates)
    transloc_nclose_pair_count, nclose_node_count = finalize_nclose_outputs(
        nclose_nodes,
        contig_data,
        not_using_nclose_node,
        saved_not_using_nclose_node,
        rpt_con,
        uncomp_node_count,
        telo_node_count,
    )
    return NCloseBuildResult(
        df=depth_df,
        no_chrY=preprocess_result.no_chrY,
        repeat_censat_data=repeat_censat_data,
        chr_len=chr_len,
        contig_data=contig_data,
        contig_data_size=contig_data_size,
        chr_corr=chr_corr,
        chr_rev_corr=chr_rev_corr,
        telo_contig=telo_contig,
        telo_node_count=telo_node_count,
        telo_set=telo_set,
        rpt_con=rpt_con,
        bnd_contig=bnd_contig,
        raw_nclose_nodes=raw_nclose_nodes,
        nclose_nodes=nclose_nodes,
        vctg_dict=vctg_dict,
        all_nclose_comp=all_nclose_comp,
        uncomp_node_count=uncomp_node_count,
        nclose_node_count=nclose_node_count,
        transloc_nclose_pair_count=transloc_nclose_pair_count,
        indel_exclude_idx_set=indel_exclude_idx_set,
        telo_coverage=preprocess_result.telo_coverage,
        contig_stage_records=tuple(preprocess_result.stage_records),
        nclose_stage_records=tuple(nclose_pipeline_state.stage_records),
        nclose_rejections=tuple(nclose_pipeline_state.rejections),
    )

def get_ori_ctg_name_data(PAF_FILE_PATH : list) -> list:
    ori_ctg_name_data = []

    for ori_paf_loc in PAF_FILE_PATH:
        ori_ctg_name_list = []
        with open(ori_paf_loc, 'r') as f:
            for paf_line in f:
                paf_line = paf_line.split('\t')
                ori_ctg_name_list.append(paf_line[CTG_NAM])

        ori_ctg_name_data.append(ori_ctg_name_list)

    return ori_ctg_name_data


def resolve_pregraph_source(
    mode,
    primary_paf,
    alt_paf,
    original_paf_paths,
    prefix,
):
    """Resolve the effective source without reading compatibility-only primary data.

    Assembly mode keeps the legacy primary/``--alt`` command shape, but only
    the aligned unitig supplied by ``--alt`` and the last original PAF are
    active inputs.  Earlier original PAF entries are accepted solely for CLI
    compatibility with the legacy 02 command.
    """

    original_paf_paths = tuple(original_paf_paths or ())
    if mode == PregraphSourceMode.VCF:
        return NCloseSourceConfig(
            mode=mode,
            paf_file_paths=(
                f"{prefix}/{VCF_SYNTHETIC_PAF_NAME}",
                primary_paf,
            ),
            original_paf_paths=(),
            is_unitig_reduced=False,
            secondary_candidate_paf=None,
        )
    if mode == PregraphSourceMode.PRIMARY_ONLY_RETRY:
        primary_original_paf = original_paf_paths[0]
        return NCloseSourceConfig(
            mode=mode,
            paf_file_paths=(primary_paf, primary_paf),
            original_paf_paths=(primary_original_paf, primary_original_paf),
            is_unitig_reduced=True,
            secondary_candidate_paf=None,
        )

    if alt_paf is None:
        raise ValueError(
            "Assembly stage 01 requires --alt with the aligned unitig PAF; "
            "the positional primary PAF is compatibility-only"
        )
    if not original_paf_paths:
        raise ValueError(
            "Assembly stage 01 requires --original-paf-loc with the "
            "original unitig PAF as its last value"
        )
    return NCloseSourceConfig(
        mode=mode,
        paf_file_paths=(alt_paf,),
        original_paf_paths=(original_paf_paths[-1],),
        is_unitig_reduced=False,
        secondary_candidate_paf=alt_paf,
    )


def import_bed(bed_path: str) -> dict:
    intervals = defaultdict(list)
    with open(bed_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            fields = line.split("\t")
            intervals[fields[0]].append((int(fields[1]), int(fields[2])))
    return intervals


_INVALIDATED_STAGE01_OUTPUTS = (
    STAGE01_HANDOFF_PKL,
    NCLOSE_NODES_PKL,
    LEGACY_NCLOSE_CHUNK_PKL,
    PIPELINE_INPUT_PKL,
    "ecdna_circuit_data.pkl",
    "paf_file_path.pkl",
    "compressed_nclose_nodes_list.txt",
    "all_nclose_nodes_list.txt",
    "nclose_nodes_index.txt",
    "nclose_event_catalog.pkl",
    "virtual_ordinary_contig.txt",
    "conjoined_type4_ins_del.pkl",
    "cen_fragment_data.pkl",
    "telomere_connected_list.txt",
    "telomere_connected_list_readable.txt",
    "indel_exclude_idx_set.pkl",
    RAW_TRANSLOCATION_CANDIDATE_PKL,
    RAW_TRANSLOCATION_RESULT_PKL,
    NCLOSE_COUNT_CANDIDATE_PKL,
    NCLOSE_COUNT_RESULT_PKL,
    VCF_TYPE4_EVENTS_PKL,
    VCF_MODE_SUMMARY_JSON,
    VCF_MODE_SUMMARY_TSV,
    VCF_SKIPPED_RECORDS_TSV,
    VCF_ORIENTATION_MISMATCH_TSV,
    VCF_TELOMERE_PAF_NODES_TSV,
    STAGE01_SUMMARY_JSON,
    STAGE01_REJECTIONS_TSV,
)


def remove_stale_stage01_outputs(prefix):
    """Invalidate final NClose/graph artifacts before rebuilding preprocessing."""

    for filename in _INVALIDATED_STAGE01_OUTPUTS:
        path = os.path.join(prefix, filename)
        if os.path.isfile(path):
            os.remove(path)


def _stage01_preprocessed_paf_path(config: Stage01Config) -> str:
    """Return the compatibility PPC name chosen from the positional PAF."""

    return os.path.join(
        config.prefix,
        f"{os.path.basename(config.paf_file_path)}.ppc.paf",
    )


def _activate_stage01_runtime(config: Stage01Config) -> None:
    """Populate the small legacy-global surface still used by post-filters."""

    global args, PREFIX, THREAD, SKIP_BAM_ANAL, read_bam_loc
    global CHROMOSOME_INFO_FILE_PATH, TELOMERE_INFO_FILE_PATH
    global REPEAT_INFO_FILE_PATH, CENSAT_PATH, main_stat_loc
    global PREPROCESSED_PAF_FILE_PATH, PAF_FILE_PATH
    global ORIGINAL_PAF_LOC_LIST, ORIGINAL_PAF_LOC_LIST_, asm2cov

    PREFIX = config.prefix
    THREAD = int(config.thread)
    SKIP_BAM_ANAL = bool(config.skip_bam_analysis)
    read_bam_loc = config.read_bam_path
    CHROMOSOME_INFO_FILE_PATH = config.reference_fai_path
    TELOMERE_INFO_FILE_PATH = config.telomere_bed_path
    REPEAT_INFO_FILE_PATH = config.repeat_bed_path
    CENSAT_PATH = config.censat_bed_path
    main_stat_loc = config.main_stat_path
    PREPROCESSED_PAF_FILE_PATH = _stage01_preprocessed_paf_path(config)
    ORIGINAL_PAF_LOC_LIST_ = list(config.original_paf_paths)
    ORIGINAL_PAF_LOC_LIST = list(config.original_paf_paths)
    PAF_FILE_PATH = []
    asm2cov = Counter()
    args = SimpleNamespace(
        paf_file_path=config.paf_file_path,
        alt=config.alt_path,
        vcf_input=config.vcf_input_path,
        vcf_filter_pass=list(config.vcf_filter_pass),
        check_nclose_count=config.check_nclose_count,
        nclose_count_vaf_threshold=config.nclose_count_vaf_threshold,
        progress=config.progress,
        exclude_nclose_list_loc=config.exclude_nclose_list_path,
        disable_alt_ctg_simple=config.disable_alt_ctg_simple,
    )


def _make_stage01_context(config: Stage01Config, source: NCloseSourceConfig):
    ori_names = (
        []
        if source.mode == PregraphSourceMode.VCF
        else get_ori_ctg_name_data(list(source.paf_file_paths))
    )
    return PregraphBuildContext(
        source=source,
        ori_ctg_name_data=ori_names,
        prefix=config.prefix,
        preprocessed_paf_path=PREPROCESSED_PAF_FILE_PATH,
        reference_fai_path=config.reference_fai_path,
        telomere_bed_path=config.telomere_bed_path,
        repeat_bed_path=config.repeat_bed_path,
        censat_bed_path=config.censat_bed_path,
        main_stat_path=config.main_stat_path,
        asm2cov=asm2cov,
        disable_alt_ctg_simple=config.disable_alt_ctg_simple,
    )


def _json_safe(value):
    """Convert extension metadata into deterministic JSON-compatible data."""

    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, dict):
        return {
            str(key): _json_safe(item)
            for key, item in sorted(value.items(), key=lambda item: str(item[0]))
        }
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (set, frozenset)):
        return sorted((_json_safe(item) for item in value), key=str)
    return repr(value)


def write_stage01_provenance(
    config: Stage01Config,
    source: NCloseSourceConfig,
    result: NCloseBuildResult,
) -> None:
    """Persist compact stage counts and first-rejection provenance."""

    if source.mode == PregraphSourceMode.VCF:
        effective_unitig_paf = (
            source.paf_file_paths[1]
            if len(source.paf_file_paths) > 1
            else None
        )
        ignored_primary_paf = None
    else:
        effective_unitig_paf = source.paf_file_paths[0]
        ignored_primary_paf = config.paf_file_path

    nclose_records = list(result.nclose_stage_records)
    initial_candidate_count = (
        int(nclose_records[0]["before_candidates"])
        if nclose_records
        else sum(len(pairs) for pairs in result.raw_nclose_nodes.values())
    )
    summary = {
        "version": 1,
        "source_mode": source.mode.value,
        "ignored_primary_paf": ignored_primary_paf,
        "effective_unitig_paf": effective_unitig_paf,
        "effective_original_unitig_paf": (
            source.original_paf_paths[-1]
            if source.original_paf_paths
            else None
        ),
        "preprocessed_paf": _stage01_preprocessed_paf_path(config),
        "contig_stages": _json_safe(result.contig_stage_records),
        "nclose_stages": _json_safe(nclose_records),
        "initial_candidate_count": initial_candidate_count,
        "final_candidate_count": sum(
            len(pairs) for pairs in result.nclose_nodes.values()
        ),
        "recorded_rejection_count": len(result.nclose_rejections),
    }
    with open(
        os.path.join(config.prefix, STAGE01_SUMMARY_JSON),
        "wt",
        encoding="utf-8",
    ) as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")

    header = (
        "contig_name",
        "node_a",
        "node_b",
        "chrom_a",
        "breakend_a",
        "direction_a",
        "chrom_b",
        "breakend_b",
        "direction_b",
        "origin",
        "stage",
        "reason",
    )
    seen = set()
    with open(
        os.path.join(config.prefix, STAGE01_REJECTIONS_TSV),
        "wt",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for rejection in result.nclose_rejections:
            if not isinstance(rejection, NCloseRejection):
                raise TypeError(
                    "Stage-01 provenance contains a non-NCloseRejection"
                )
            candidate = rejection.candidate
            if candidate.identity in seen:
                continue
            seen.add(candidate.identity)
            node_a_idx, node_b_idx = candidate.path_pair
            node_a = result.contig_data[node_a_idx]
            node_b = result.contig_data[node_b_idx]
            writer.writerow((
                candidate.contig_name,
                node_a_idx,
                node_b_idx,
                node_a[CHR_NAM],
                get_breakend_coord(node_a, 0),
                node_a[CTG_DIR],
                node_b[CHR_NAM],
                get_breakend_coord(node_b, 1),
                node_b[CTG_DIR],
                candidate.origin,
                rejection.stage,
                rejection.reason,
            ))


def persist_stage01_outputs(
    config: Stage01Config,
    source: NCloseSourceConfig,
    result: NCloseBuildResult,
):
    """Write ecDNA, downstream NClose state, and the exact stage-10 handoff."""

    if config.vcf_input_path is not None:
        assert_vcf_nclose_has_no_indel_like(
            result.contig_data,
            result.nclose_nodes,
            "final nclose_nodes.pkl",
        )

    save_ecdna_outputs(
        config.prefix,
        result.contig_data,
        result.raw_nclose_nodes,
        result.all_nclose_comp,
    )
    save_pipeline_input(
        config.prefix,
        vcf_input=config.vcf_input_path is not None,
        vcf_input_path=(
            os.path.abspath(config.vcf_input_path)
            if config.vcf_input_path is not None
            else None
        ),
    )
    save_nclose_nodes(config.prefix, result.nclose_nodes)

    handoff_path = os.path.join(config.prefix, STAGE01_HANDOFF_PKL)
    save_stage10_input(
        handoff_path,
        result.contig_data,
        result.nclose_nodes,
        result.telo_contig,
    )
    with open(os.path.join(config.prefix, "paf_file_path.pkl"), "wb") as handle:
        pkl.dump(list(source.paf_file_paths), handle)
    write_stage01_provenance(config, source, result)
    return handoff_path


def run_stage01(
    config: Stage01Config,
    candidate_builder=None,
    *,
    pipeline=None,
):
    """Run stage 01 with an explicit, composable implementation policy."""

    if config.thread < 1:
        raise ValueError("thread must be positive")
    if pipeline is None:
        pipeline = default_stage01_pipeline()
    if not isinstance(pipeline, Stage01Pipeline):
        raise TypeError("pipeline must be a Stage01Pipeline")
    if candidate_builder is not None:
        pipeline = replace(pipeline, candidate_builder=candidate_builder)
    if not callable(pipeline.candidate_builder):
        raise TypeError("Stage-01 candidate builder must be callable")
    if config.vcf_input_path is not None:
        source = resolve_pregraph_source(
            PregraphSourceMode.VCF,
            config.paf_file_path,
            config.alt_path,
            (),
            config.prefix,
        )
    else:
        source = resolve_pregraph_source(
            PregraphSourceMode.CONFIGURED_PAF,
            config.paf_file_path,
            config.alt_path,
            config.original_paf_paths,
            config.prefix,
        )
        if len(source.paf_file_paths) != len(source.original_paf_paths):
            raise ValueError(
                "Assembly input requires one --original-paf-loc per aligned PAF"
            )

    os.makedirs(config.prefix, exist_ok=True)
    remove_stale_stage01_outputs(config.prefix)
    preprocessed_paf_path = _stage01_preprocessed_paf_path(config)
    if os.path.isfile(preprocessed_paf_path):
        os.remove(preprocessed_paf_path)
    _activate_stage01_runtime(config)

    global PAF_FILE_PATH, ORIGINAL_PAF_LOC_LIST, is_unitig_reduced
    PAF_FILE_PATH = list(source.paf_file_paths)
    ORIGINAL_PAF_LOC_LIST = list(source.original_paf_paths)
    is_unitig_reduced = source.is_unitig_reduced
    context = _make_stage01_context(config, source)

    if source.mode == PregraphSourceMode.VCF:
        result = build_vcf_mode_inputs(context)
    else:
        preprocess_result = contig_preprocessing_00(
            context,
            contig_stages=pipeline.contig_stages,
        )
        result = nclose_calc(
            context,
            preprocess_result,
            candidate_builder=pipeline.candidate_builder,
            nclose_stages=pipeline.nclose_stages,
        )

    persist_stage01_outputs(config, source, result)
    return result
