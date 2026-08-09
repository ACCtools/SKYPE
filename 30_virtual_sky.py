import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from skype_utils import *
from nclose_tracking import bnd_event_keys, load_event_catalog, load_path_usage
from reference_path_clustering import (
    ReferenceSpan,
    cluster_reference_paths,
    extract_ordered_reference_spans,
    is_full_reference_path,
)
from reference_path_depth import (
    PathDepthProfiles,
    build_cluster_depth_tracks,
    build_depth_bin_index,
    load_model_copy_number,
    observed_copy_number,
)
from reference_path_plotting import DEPTH_MODES, render_reference_path_clusters

import numpy as np
import pandas as pd
from virtual_sky_plotting import (
    render_karyotype_diagram as render_virtual_sky,
    weighted_vaf_read_depth,
)
import pickle as pkl

import csv
import re
import ast
import glob
import logging
import argparse

from collections import defaultdict
from collections import Counter

# logging 설정(레벨/포맷)은 skype_utils 에서 중앙 관리한다 (LOG_LEVEL).
logging.info("30_virtual_sky start")

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
CTG_MAINFLOWDIR = 18
CTG_MAINFLOWCHR = 19

DIR_FOR = 1
DIR_BAK = 1

ABS_MAX_COVERAGE_RATIO = 3
MAX_PATH_CNT = 100
INF = 1000000000

DIR_FOR = 1
TELOMERE_EXPANSION = 5 * K

BND_TYPE = 0
CTG_IN_TYPE = 1
TEL_TYPE = 2

MAJOR_BASELINE = 0.6

TARGET_DEPTH = 0.2

MEANDEPTH_FLANKING_LENGTH = 5*M
TYPE34_BREAK_CHUKJI_LIMIT = 1*M

TYPE4_CLUSTER_SIZE = 10 * M
TYPE4_MEANDEPTH_FLANKING_LENGTH = 500 * K
NCLOSE_SIM_COMPARE_RATIO = 1.2
NCLOSE_SIM_DIFF_THRESHOLD = 5
RAW_TRANSLOCATION_RESULT_PKL = 'raw_translocation_result.pkl'
RAW_TRANSLOCATION_REPORT_TSV = 'raw_translocation_read_counts.tsv'
CLUSTER_NCLOSE_COUNT_RESULT_PKL = 'cluster_nclose_count_result.pkl'
NCLOSE_COUNT_MAX_QUERY_GAP = 5 * K

JOIN_BASELINE = 0.8
KARYOTYPE_SECTION_MINIMUM_LENGTH = 100 * K
KARYOTYPE_MIN_SEGMENT_LENGTH = 1 * M  # karyotype 텍스트 표기 시 이보다 짧은 segment/indel 은 무시 (1Mb)
KARYOTYPE_NORMAL_RATIO = 0.9  # 단일 염색체 path 가 reference 길이의 이 비율 미만이면 normal('1') 대신 del 로 취급
REFERENCE_PATH_CLUSTER_MIN_OVERLAP = 100 * K
# reference-sharing cluster 그림 종류. SKYPE_REFSPAN_MODES 로 하나만 골라 그릴 수 있다.
# plain: 기존 branch graph / lane: lane 안 depth ribbon /
# profile: chromosome 별 CN 패널 / attribution: region 별 기여도 막대
REFERENCE_PATH_CLUSTER_FIGURE_MODES = tuple(
    mode
    for mode in os.environ.get(
        'SKYPE_REFSPAN_MODES', 'plain,lane,profile,attribution'
    ).replace(' ', '').split(',')
    if mode
)

CTG_INTYPE_CHECK_MIN_LENGTH = 100 * K
CTG_INTYPE_INSERT_MIN_SEGMENT_LENGTH = 10 * K
CTG_INTYPE_INSERT_MIN_RATIO = 0.2



NODE_NAME = 1
CHR_CHANGE_IDX = 2
DIR_CHANGE_IDX = 3

def similar_check(v1, v2, ratio):
    try:
        assert(v1 >= 0 and v2 >= 0)
    except:
        logging.error(f"Invalid values for similarity check: v1={v1}, v2={v2}")
        assert(False)
    mi, ma = sorted([v1, v2])
    return False if mi == 0 else (ma / mi <= ratio) or ma-mi < NCLOSE_SIM_DIFF_THRESHOLD

def exist_near_bnd(chrom, inside_st, inside_nd, ratio=NCLOSE_SIM_COMPARE_RATIO):
    # subset of df for the given chromosome
    df_chr = df[df['chr'] == chrom]

    def mean_depth(start, end):
        """Return mean meandepth over windows overlapping [start, end)."""
        mask = (df_chr['nd'] > start) & (df_chr['st'] < end)
        return df_chr.loc[mask, 'meandepth'].mean()

    # for inside_st
    st_depth = mean_depth(inside_st - MEANDEPTH_FLANKING_LENGTH, inside_st)
    nd_depth = mean_depth(inside_nd, inside_nd + MEANDEPTH_FLANKING_LENGTH)
    
    if np.isnan(st_depth) or np.isnan(nd_depth):
        return True

    # print(chrom, inside_st, inside_nd, not similar_check(st_depth, nd_depth))
    return not similar_check(st_depth, nd_depth, ratio)

def check_near_type4(chrom, inside_st, inside_nd):
    # subset of df for the given chromosome
    df_chr = df[df['chr'] == chrom]

    def mean_depth(start, end):
        """Return mean meandepth over windows overlapping [start, end)."""
        mask = (df_chr['nd'] > start) & (df_chr['st'] < end)
        return df_chr.loc[mask, 'meandepth'].mean()

    # for inside_st
    st_depth = mean_depth(inside_st - TYPE4_MEANDEPTH_FLANKING_LENGTH, inside_st)
    nd_depth = mean_depth(inside_nd, inside_nd + TYPE4_MEANDEPTH_FLANKING_LENGTH)
    if np.isnan(st_depth) or np.isnan(nd_depth):
        return True

    # print(chrom, inside_st, inside_nd, not similar_check(st_depth, nd_depth))
    return not similar_check(st_depth, nd_depth, NCLOSE_SIM_COMPARE_RATIO)

def chr2int(x):
    if x.startswith('chr'):
        chrXY2int = {'chrX' : 24, 'chrY' : 25}
        if x in chrXY2int:
            return chrXY2int[x]
        else:
            return int(x[3:])
    else:
        return INF

def find_chr_len(file_path : str) -> dict:
    chr_data_file = open(file_path, "r")
    chr_len = {}
    for curr_data in chr_data_file:
        curr_data = curr_data.split("\t")
        chr_len[curr_data[0]] = int(curr_data[1])
    chr_data_file.close()
    return chr_len

def import_ppc_data(file_path : str) -> list :
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

def import_paf_data(file_path : str) -> list :
    contig_data = []
    int_induce_idx = [1, 2, 3, 6, 7, 8, 9]
    idx = 0
    with open(file_path, "r") as paf_file:
        for curr_contig in paf_file:
            curr_contig = curr_contig.rstrip()
            a = curr_contig.split("\t")
            temp_list = a[:9]
            temp_list.append(a[11])
            for i in int_induce_idx:
                temp_list[i] = int(temp_list[i])
            temp_list.append(idx)
            contig_data.append(temp_list)
            idx+=1
    return contig_data

def import_index_path(file_path : str) -> list:
    file_path_list = file_path.split('/')
    key = file_path_list[-2]
    cnt = int(file_path_list[-1].split('.')[0]) - 1

    return path_list_dict[key][cnt][0]

def import_path_paf_rows(file_path : str) -> list:
    rows = []
    with open(file_path, "r") as paf_file:
        for curr_contig in paf_file:
            curr_contig = curr_contig.rstrip()
            if not curr_contig:
                continue
            temp_list = curr_contig.split("\t")
            int_induce_idx = [
                CTG_LEN, CTG_STR, CTG_END,
                CHR_LEN, CHR_STR, CHR_END,
                CTG_MAPQ,
            ]
            for i in int_induce_idx:
                temp_list[i] = int(temp_list[i])
            rows.append(tuple(temp_list))
    return rows

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

def extract_telomere_connect_contig(telo_info_path : str) -> list:
    telomere_connect_contig = []
    with open(telo_info_path) as f:
        for curr_data in f:
            curr_data = curr_data.rstrip()
            temp_list = curr_data.split("\t")
            chr_info = temp_list[0]
            contig_id = ast.literal_eval(temp_list[1])
            telomere_connect_contig.append((chr_info, contig_id[1]))
    
    return telomere_connect_contig

def distance_checker(node_a : tuple, node_b : tuple) -> int :
    if max(int(node_a[0]), int(node_b[0])) < min(int(node_a[1]), int(node_b[1])):
        return 0   
    else:
        return min(abs(int(node_b[0]) - int(node_a[1])), abs(int(node_b[1]) - int(node_a[0])))

def telo_condition(node : list, need_label_index : dict) -> bool:
    return node in need_label_index

def virtual_event_label(event_type : str) -> str:
    return {'d': 'del', 'i': 'ins', 'v': 'inv'}[event_type]

def parse_chromosome_labels(s):
    """
    Parse '...<f|b>_...<f|b>' into a canonical tuple:
      (left_label, left_is_f, right_label, right_is_f)

    Rules:
    - The two ends must end with 'f' or 'b' (assert if not).
    - Keep each '...' label string intact (e.g., 'chr12', 'scaf_007', etc.).
    - Canonicalize by sorting so the lexicographically smaller label comes first.
      If labels are equal, put 'f' (True) before 'b' (False).
    - When swapping due to sorting, directions stay attached to their original labels.
      (So 'chr12f_chr1b' becomes ('chr1', True, 'chr12', False).)
    """
    m = re.fullmatch(r'(.+?)([fb])_(.+?)([fb])', s)
    assert m is not None, "Input must match ...(f|b)_...(f|b) pattern"

    a_label, a_dir_ch, b_label, b_dir_ch = m.groups()
    a_is_f = (a_dir_ch == 'f')
    b_is_f = (b_dir_ch == 'f')

    # Canonical order by label; if same label, 'f' (True) first.
    if (a_label > b_label) or (a_label == b_label and not a_is_f and b_is_f):
        # Swap ends to enforce canonical order; keep directions with their labels.
        a_label, b_label = b_label, a_label
        a_is_f, b_is_f = b_is_f, a_is_f

    return (a_label, a_is_f, b_label, b_is_f)


def max_aligned_match_length(
    seq_a: list[tuple[tuple[str, str], int]],
    seq_b: list[tuple[tuple[str, str], int]],
) -> int:
    """
    Return the maximum total matched length after sliding two piecewise-constant
    label sequences along one axis. A and B are lists of ((chrom, strand), length).
    Only regions with exactly the same (chrom, strand) contribute to the score.

    Algorithm:
      1) Convert each sequence into absolute intervals [(start, end, label)].
      2) Consider candidate shifts = {a_ep - b_ep | a_ep in endpoints(A), b_ep in endpoints(B)}.
         (The overlap configuration only changes when an endpoint meets another.)
      3) For each shift, line-sweep over the two interval lists and accumulate
         overlap length where labels are equal.
      4) Return the maximum accumulated length across all shifts.

    Time complexity:
      Let n, m be #segments. Endpoints ~ (n+1), (m+1).
      Candidates O((n+1)*(m+1)); each evaluation O(n+m). Works well for tens~hundreds of segments.
    """
    # --- build absolute intervals: [(start, end, label)] and endpoint lists ---
    def build_intervals(seq):
        intervals = []
        endpoints = []
        pos = 0
        endpoints.append(pos)
        for (label, length) in seq:
            start = pos
            end = pos + length
            intervals.append((start, end, label))
            pos = end
            endpoints.append(pos)
        return intervals, endpoints

    A, A_ep = build_intervals(seq_a)
    B, B_ep = build_intervals(seq_b)

    if not A or not B:
        return 0

    # --- generate candidate shifts (all endpoint differences) ---
    # shift d means: compare A intervals with B intervals shifted by +d
    candidates = set()
    for a_e in A_ep:
        for b_e in B_ep:
            candidates.add(a_e - b_e)

    # --- overlap length for a given shift ---
    def match_length_for_shift(d: int) -> int:
        i, j = 0, 0
        total = 0
        # Two-pointer sweep over A and shifted-B
        while i < len(A) and j < len(B):
            a_s, a_e, a_lab = A[i]
            b_s, b_e, b_lab = B[j]
            b_s += d
            b_e += d

            # If no overlap, advance the one that ends earlier / starts later
            if a_e <= b_s:
                i += 1
                continue
            if b_e <= a_s:
                j += 1
                continue

            # Overlapping segment
            ov_s = a_s if a_s > b_s else b_s
            ov_e = a_e if a_e < b_e else b_e
            if ov_e > ov_s and a_lab == b_lab:
                total += (ov_e - ov_s)

            # Advance the interval that ends first
            if a_e <= b_e:
                i += 1
            else:
                j += 1
        return total

    best = 0
    # (Optional) small heuristic: iterate over sorted candidates for deterministic behavior
    for d in sorted(candidates):
        val = match_length_for_shift(d)
        if val > best:
            best = val

    return best

def should_join_by_baseline(
    seq_a: list[tuple[tuple[str, str], int]],
    seq_b: list[tuple[tuple[str, str], int]]
) -> bool:
    """
    Decide if two sequences should be joined based on:
      max_aligned_match_length(seq_a, seq_b) / max(total_len_a, total_len_b) >= JOIN_BASELINE

    Notes:
      - Returns False if both sequences have total length 0 (to avoid 0-division).
      - Assumes non-negative lengths.
      - Threshold is inclusive (>=).

    """

    total_a = sum(length for (_, length) in seq_a)
    total_b = sum(length for (_, length) in seq_b)
    denom = total_a if total_a >= total_b else total_b
    if denom == 0:
        return False

    score = max_aligned_match_length(seq_a, seq_b)
    return (score / denom) >= JOIN_BASELINE

def append_karyotype_piece(pieces : list, chrom : str, strand : str, length : int, merge : bool = True):
    if length <= 0:
        return
    key = (chrom, strand)
    if merge and pieces and pieces[-1][0] == key:
        pieces[-1] = (key, pieces[-1][1] + length)
    else:
        pieces.append((key, length))

def append_ctg_intype_interrupt_piece(pieces : list, chrom : str, strand : str, length : int):
    if length <= 0:
        return
    if pieces and pieces[-1][0][0] == chrom:
        prev_chrom, prev_strand = pieces[-1][0]
        pieces[-1] = ((prev_chrom, prev_strand), pieces[-1][1] + length)
    else:
        pieces.append(((chrom, strand), length))

def infer_path_strand_from_paf_rows(rows : list, idx : int) -> str:
    chrom = rows[idx][CHR_NAM]
    curr_mid = (rows[idx][CHR_STR] + rows[idx][CHR_END]) // 2

    if idx + 1 < len(rows) and rows[idx + 1][CHR_NAM] == chrom:
        next_mid = (rows[idx + 1][CHR_STR] + rows[idx + 1][CHR_END]) // 2
        if next_mid != curr_mid:
            return '+' if next_mid > curr_mid else '-'

    if idx > 0 and rows[idx - 1][CHR_NAM] == chrom:
        prev_mid = (rows[idx - 1][CHR_STR] + rows[idx - 1][CHR_END]) // 2
        if curr_mid != prev_mid:
            return '+' if curr_mid > prev_mid else '-'

    return rows[idx][CTG_DIR] if rows[idx][CTG_DIR] in {'+', '-'} else '+'

def get_ctg_intype_interrupt_rows(key_int : int, endpoint_chroms : set) -> list:
    paf_loc = f"{output_folder}/{key_int}.paf"
    if not os.path.isfile(paf_loc):
        return []

    rows = import_path_paf_rows(paf_loc)
    if not rows:
        return []

    row_lengths = [abs(row[CHR_END] - row[CHR_STR]) for row in rows]
    total_length = sum(row_lengths)
    if total_length <= CTG_INTYPE_CHECK_MIN_LENGTH:
        return []

    foreign_chrom_lengths = Counter()
    for row, length in zip(rows, row_lengths):
        chrom = row[CHR_NAM]
        if chrom not in endpoint_chroms:
            foreign_chrom_lengths[chrom] += length

    interrupt_chroms = {
        chrom for chrom, length in foreign_chrom_lengths.items()
        if length / total_length >= CTG_INTYPE_INSERT_MIN_RATIO
    }
    if not interrupt_chroms:
        return []

    selected_rows = []
    for i, (row, length) in enumerate(zip(rows, row_lengths)):
        chrom = row[CHR_NAM]
        if chrom not in interrupt_chroms:
            continue
        if length < CTG_INTYPE_INSERT_MIN_SEGMENT_LENGTH:
            continue
        strand = infer_path_strand_from_paf_rows(rows, i)
        selected_rows.append((row, strand, length))
    return selected_rows

def get_ctg_intype_interrupt_pieces(key_int : int, endpoint_chroms : set) -> list:
    pieces = []
    for row, strand, length in get_ctg_intype_interrupt_rows(
        key_int, endpoint_chroms
    ):
        chrom = row[CHR_NAM]
        append_ctg_intype_interrupt_piece(pieces, chrom, strand, length)
    return pieces

def get_ctg_intype_interrupt_reference_spans(
    key_int : int, endpoint_chroms : set
) -> tuple:
    return tuple(
        ReferenceSpan(
            str(row[CHR_NAM]),
            min(int(row[CHR_STR]), int(row[CHR_END])),
            max(int(row[CHR_STR]), int(row[CHR_END])),
            strand,
            "interrupt",
        )
        for row, strand, _length in get_ctg_intype_interrupt_rows(
            key_int, endpoint_chroms
        )
    )

def get_ctg_intype_key_by_edge(path_path : str) -> dict:
    key_by_edge = {}
    for key_int in path2key_int_list.get(path_path, []):
        key_type, key_value = int2key[key_int]
        if key_type == CTG_IN_TYPE:
            key_by_edge[key_value] = key_int
    return key_by_edge

def get_path_ctg_intype_interrupt_spans(path_path : str) -> dict:
    path = import_index_path(path_path)
    key_by_edge = get_ctg_intype_key_by_edge(path_path)
    spans_by_edge = {}
    for previous_step, current_step in zip(path, path[1:]):
        previous_id = previous_step[NODE_NAME] if len(previous_step) >= 4 else None
        current_id = current_step[NODE_NAME] if len(current_step) >= 4 else None
        if not isinstance(previous_id, int) or not isinstance(current_id, int):
            continue
        edge_key = (
            (previous_step[0], previous_id),
            (current_step[0], current_id),
        )
        key_int = key_by_edge.get(edge_key)
        if key_int is None:
            continue
        endpoint_chroms = {
            ppc_data[previous_id][CHR_NAM],
            ppc_data[current_id][CHR_NAM],
        }
        interrupt_spans = get_ctg_intype_interrupt_reference_spans(
            key_int, endpoint_chroms
        )
        if interrupt_spans:
            spans_by_edge[edge_key] = interrupt_spans
    return spans_by_edge

def get_karyotype_summary_from_index(path_path : str, type4_edge_to_event_key=None,
                                     type4_event_by_key=None, junction_event_keys=None,
                                     nclose_boundaries=None) -> list:
    """
    Fallback summary from the compact index path. This keeps the old behavior
    for prefixes that do not have the 21_pat_depth PAF fragments available.
    """
    pieces = []
    junction_event_keys = set(junction_event_keys or ())
    path = import_index_path(path_path)
    ctg_intype_key_by_edge = get_ctg_intype_key_by_edge(path_path)
    
    # Padding for easier calculation
    if len(path[0]) < 4:
        path[0] = tuple([0] + list(path[0]))
    if len(path[-1]) < 4:
        path[-1] = tuple([0] + list(path[-1]))
        
    # Initialize direction and chromosome from the first dummy node
    curr_incr = '+' if path[0][NODE_NAME][-1] == 'f' else '-'
    
    # Set the starting reference using the first real node (index 1)
    # instead of assuming the absolute ends of the chromosome (0 or chr_len)
    first_real_node = ppc_data[path[1][NODE_NAME]]
    # Take chr from the real node's CHR_NAM, not the telomere endpoint name:
    # the shared chrX/chrY telomere node is always labeled chrXf/chrXb, so a
    # pure-chrY path with no chr/dir-change transition would be mislabeled chrX.
    curr_chr = [first_real_node[CHR_NAM], curr_incr]
    curr_ref = first_real_node[CHR_STR] if curr_incr == '+' else first_real_node[CHR_END]
    
    for i in range(1, len(path)-1):
        prev_node_name = path[i-1][NODE_NAME]
        curr_node_name = path[i][NODE_NAME]
        if not isinstance(prev_node_name, int) or not isinstance(curr_node_name, int):
            continue

        last_node = ppc_data[prev_node_name]
        curr_node = ppc_data[curr_node_name]
        edge_key = (
            (path[i-1][0], prev_node_name),
            (path[i][0], curr_node_name),
        )
        type4_event_key = None
        type4_event = None
        if type4_edge_to_event_key is not None:
            type4_event_key = type4_edge_to_event_key.get(edge_key)
        if type4_event_key is not None and type4_event_by_key is not None:
            type4_event = type4_event_by_key.get(type4_event_key)
        type4_deletion_edge = (
            type4_event is not None and type4_event.get("event_type") == "d"
        )
        ctg_intype_key_int = ctg_intype_key_by_edge.get(edge_key)
        interrupt_pieces = []
        if ctg_intype_key_int is not None:
            endpoint_chroms = {last_node[CHR_NAM], curr_node[CHR_NAM]}
            interrupt_pieces = get_ctg_intype_interrupt_pieces(ctg_intype_key_int, endpoint_chroms)

        if path[i][CHR_CHANGE_IDX] > path[i-1][CHR_CHANGE_IDX] \
        or path[i][DIR_CHANGE_IDX] > path[i-1][DIR_CHANGE_IDX] \
        or interrupt_pieces \
        or type4_deletion_edge:
            
            # Add last piece
            if curr_incr == '+':
                append_karyotype_piece(pieces, curr_chr[0], curr_chr[1], abs(last_node[CHR_END] - curr_ref), merge=False)
            else:
                append_karyotype_piece(pieces, curr_chr[0], curr_chr[1], abs(curr_ref - last_node[CHR_STR]), merge=False)

            ordinary_event_key = tuple(sorted((prev_node_name, curr_node_name)))
            junction_event_key = None
            if ordinary_event_key in junction_event_keys:
                junction_event_key = ordinary_event_key
            elif type4_event_key in junction_event_keys:
                junction_event_key = type4_event_key
            if nclose_boundaries is not None and junction_event_key is not None:
                nclose_boundaries.append((
                    len(pieces),
                    junction_event_key,
                ))

            for piece_chr, piece_length in interrupt_pieces:
                append_karyotype_piece(pieces, piece_chr[0], piece_chr[1], piece_length, merge=True)
                        
            # Update info of new piece (starting ref, chromosome type, increment ..)
            if path[i][NODE_NAME] > path[i-1][NODE_NAME]:
                curr_incr = curr_node[CTG_DIR]
                curr_chr = [curr_node[CHR_NAM], curr_incr]
                curr_ref = curr_node[CHR_STR] if curr_incr == '+' else curr_node[CHR_END]
            else:
                curr_incr = '-' if curr_node[CTG_DIR] == '+' else '+'
                curr_chr = [curr_node[CHR_NAM], curr_incr]
                curr_ref = curr_node[CHR_STR] if curr_incr == '+' else curr_node[CHR_END]
    
    # Process the final piece using the last real node (index -2)
    # instead of extending it to absolute end
    last_real_node = ppc_data[path[-2][NODE_NAME]]
    if curr_incr == '+':
        final_length = last_real_node[CHR_END] - curr_ref
    else:
        final_length = curr_ref - last_real_node[CHR_STR]
        
    append_karyotype_piece(pieces, curr_chr[0], curr_chr[1], abs(final_length), merge=False)
    return pieces

def get_karyotype_summary(non_type4_path_list: list, type4_edge_to_event_key=None,
                          type4_event_by_key=None, junction_event_keys=None,
                          path_nclose_boundaries=None):
    """
    Summarizes karyotype data from the compact index path. Only CTG_IN_TYPE
    edges are inspected in the expanded PAF fragments to reveal long inserted
    sequence from chromosomes outside the edge endpoints.
    """
    karyotypes_data_direction_include = {}
    
    for path_path in non_type4_path_list:
        nclose_boundaries = []
        pieces = get_karyotype_summary_from_index(
            path_path, type4_edge_to_event_key, type4_event_by_key,
            junction_event_keys=junction_event_keys,
            nclose_boundaries=nclose_boundaries,
        )
        karyotypes_data_direction_include[path_path] = pieces
        if path_nclose_boundaries is not None:
            path_nclose_boundaries[path_path] = nclose_boundaries
        
    return karyotypes_data_direction_include


def ecdna_format(x:int) -> str:
    if x >= 1_000_000_000:
        return f"{x / 1_000_000_000:.2f}G"
    elif x >= 1_000_000:
        return f"{x / 1_000_000:.2f}M"
    elif x >= 1_000:
        return f"{x / 1_000:.2f}K"
    else:
        return str(x)

def chrom_to_iscn(chrom : str) -> str:
    """'chr1' -> '1', 'chrX' -> 'X' (plot_virtual_chromosome 라벨과 동일 규칙)."""
    return chrom[3:] if chrom.startswith('chr') else chrom

def parse_optional_float(value):
    if value is None or value == '*':
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if np.isfinite(parsed) else None

def load_cluster_nclose_read_evidence(prefix):
    result_path = f'{prefix}/{CLUSTER_NCLOSE_COUNT_RESULT_PKL}'
    if not os.path.isfile(result_path):
        return None

    with open(result_path, 'rb') as f:
        records = pkl.load(f)

    evidence = {}
    for record in records:
        raw_event_key = record.get('event_key')
        if raw_event_key is not None:
            event_key = tuple(raw_event_key)
        else:
            try:
                event_key = tuple(sorted(int(index) for index in record['nclose_key']))
            except (KeyError, TypeError, ValueError):
                continue
        query_gap = record.get('nclose_query_gap')
        count_available = (
            query_gap is not None and
            int(query_gap) <= NCLOSE_COUNT_MAX_QUERY_GAP
        )
        if count_available:
            try:
                read_count = int(record.get('read_counts', {}).get('d2', 0))
            except (TypeError, ValueError):
                read_count = None
            read_depth = weighted_vaf_read_depth(record)
        else:
            read_count = None
            read_depth = None
        value = {
            'junction_reads': read_count,
            'vaf_read_depth': read_depth,
        }
        previous = evidence.setdefault(event_key, value)
        if previous != value:
            raise ValueError(
                f'Conflicting cluster raw-read evidence for NClose {event_key}: '
                f'{previous} != {value}'
            )
        if raw_event_key is not None:
            node_pair = tuple(sorted(int(index) for index in record.get('nclose_key', ())))
            if len(node_pair) == 2:
                previous = evidence.setdefault(node_pair, value)
                if previous != value:
                    raise ValueError(
                        f'Conflicting cluster raw-read evidence for node pair {node_pair}: '
                        f'{previous} != {value}'
                    )
    return evidence

def load_cluster_nclose_read_counts(prefix):
    evidence = load_cluster_nclose_read_evidence(prefix)
    if evidence is None:
        return None
    return {
        event_key: value['junction_reads']
        for event_key, value in evidence.items()
    }

def load_raw_translocation_report(prefix):
    report_path = f'{prefix}/{RAW_TRANSLOCATION_REPORT_TSV}'
    if not os.path.isfile(report_path):
        return {}

    rows_by_pair = {}
    with open(report_path, 'r') as f:
        for row in csv.DictReader(f, delimiter='\t'):
            try:
                pair_id = int(row['pair_id'])
            except (KeyError, TypeError, ValueError):
                continue
            rows_by_pair[pair_id] = row
    return rows_by_pair

def finite_mean(values):
    finite_values = [
        float(value) for value in values
        if value is not None and np.isfinite(value)
    ]
    if not finite_values:
        return None
    return sum(finite_values) / len(finite_values)

def estimate_raw_virtual_inv_depth(record, report_row=None):
    estimate = record.get('depth_weighted_nclose_estimate', {})
    expected_depth = parse_optional_float(estimate.get('weighted_expected_nclose_depth'))
    if expected_depth is not None:
        return expected_depth

    if report_row is not None:
        expected_depth = parse_optional_float(report_row.get('weighted_expected_nclose_depth'))
        if expected_depth is not None:
            return expected_depth

    counts = record.get('read_counts', {})
    point_depth = record.get('point_500k_depth', {})
    side_inputs = [
        (point_depth.get('point_a', {}), counts.get('d1', 0), counts.get('d2', 0)),
        (point_depth.get('point_b', {}), counts.get('d4', 0), counts.get('d3', 0)),
    ]

    weighted_sum = 0.0
    weight_sum = 0
    for depth_pair, nclose_count, point_count in side_inputs:
        point_mean = finite_mean([depth_pair.get('front'), depth_pair.get('back')])
        if point_mean is None:
            continue
        try:
            weight = int(nclose_count) + int(point_count)
            nclose_count = int(nclose_count)
        except (TypeError, ValueError):
            continue
        if weight <= 0:
            continue
        weighted_sum += point_mean * (nclose_count / weight) * weight
        weight_sum += weight

    if weight_sum == 0:
        return None
    return weighted_sum / weight_sum

def record_depth_is_balanced(record, report_row=None):
    if 'depth_balanced_translocation' in record:
        return bool(record.get('depth_balanced_translocation'))
    # The TSV is emitted after the depth-balance filter, so an old-schema pkl can
    # be treated as balanced only when its pair_id still exists in the TSV.
    return report_row is not None

def raw_false_value(value):
    return value is False or value == 0 or value == 'False' or value == 'false'

def record_has_both_point_spans(record):
    raw_no_span = record.get('raw_point_no_spanning')
    if isinstance(raw_no_span, dict):
        return raw_false_value(raw_no_span.get('point_a')) and raw_false_value(raw_no_span.get('point_b'))

    side_records = record.get('side_records', [])
    if len(side_records) < 2:
        return False
    side_flags = [
        side.get('raw_point_no_spanning', side.get('no_spanning_rawread'))
        for side in side_records[:2]
    ]
    return raw_false_value(side_flags[0]) and raw_false_value(side_flags[1])

def raw_virtual_inv_vafs(record, report_row=None):
    estimate = record.get('depth_weighted_nclose_estimate', {})
    point_a_vaf = parse_optional_float(estimate.get('point_a_nclose_vaf'))
    point_b_vaf = parse_optional_float(estimate.get('point_b_nclose_vaf'))
    if point_a_vaf is not None and point_b_vaf is not None:
        return point_a_vaf, point_b_vaf

    if report_row is not None:
        if point_a_vaf is None:
            point_a_vaf = parse_optional_float(report_row.get('point_a_nclose_vaf'))
        if point_b_vaf is None:
            point_b_vaf = parse_optional_float(report_row.get('point_b_nclose_vaf'))
        if point_a_vaf is not None and point_b_vaf is not None:
            return point_a_vaf, point_b_vaf

    counts = record.get('read_counts', {})
    try:
        d1 = int(counts.get('d1', 0))
        d2 = int(counts.get('d2', 0))
        d3 = int(counts.get('d3', 0))
        d4 = int(counts.get('d4', 0))
    except (TypeError, ValueError):
        return point_a_vaf, point_b_vaf

    if point_a_vaf is None and d1 + d2 > 0:
        point_a_vaf = d1 / (d1 + d2)
    if point_b_vaf is None and d4 + d3 > 0:
        point_b_vaf = d4 / (d4 + d3)
    return point_a_vaf, point_b_vaf

def record_passes_virtual_inv_vaf(record, report_row=None, min_vaf=RAW_VIRTUAL_INV_MIN_VAF):
    point_a_vaf, point_b_vaf = raw_virtual_inv_vafs(record, report_row)
    if point_a_vaf is None or point_b_vaf is None:
        return False
    return point_a_vaf > min_vaf and point_b_vaf > min_vaf

def directed_entry_coord(endpoint):
    return int(endpoint['ref_st']) if endpoint['dir'] == '+' else int(endpoint['ref_nd'])

def record_display_points(record, report_row=None):
    if report_row is not None:
        chrom_a = report_row.get('chrom_a')
        chrom_b = report_row.get('chrom_b')
        try:
            point_a = int(report_row.get('point_a'))
            point_b = int(report_row.get('point_b'))
        except (TypeError, ValueError):
            chrom_a = chrom_b = None
        else:
            if chrom_a and chrom_b:
                return chrom_a, point_a, chrom_b, point_b

    chrom_pair = record.get('chrom_pair', ('*', '*'))
    layout_a = record.get('layout_a', {})
    layout_b = record.get('layout_b', {})
    endpoints_a = list(layout_a.get('endpoints', ()))
    endpoints_b = list(layout_b.get('endpoints', ()))
    side_records = record.get('side_records', [])

    chrom_a = chrom_pair[0] if len(chrom_pair) > 0 else '*'
    chrom_b = chrom_pair[1] if len(chrom_pair) > 1 else '*'
    coord_a = None
    coord_b = None

    if len(endpoints_b) > 0:
        coord_a = directed_entry_coord(endpoints_b[0])
    elif len(side_records) > 0:
        coord_a = int(side_records[0]['inner_st'])

    if len(endpoints_a) > 1:
        coord_b = directed_entry_coord(endpoints_a[1])
    elif len(side_records) > 1:
        coord_b = int(side_records[1]['inner_nd'])

    return chrom_a, coord_a, chrom_b, coord_b

def point_to_chrom_end_interval(chrom, point, side, chrom_lengths):
    chrom_len = int(chrom_lengths[chrom])
    point = max(0, min(int(point), chrom_len))
    if side == 'left':
        st, nd = 0, point
    else:
        st, nd = point, chrom_len
    if nd <= st:
        return None
    return st, nd

def layout_side(record, layout_name, side_idx, default='right'):
    sides = record.get(layout_name, {}).get('sides', ())
    if side_idx < len(sides):
        return sides[side_idx]
    return default

def read_component_ref_intervals(prefix, key_int, cache):
    if key_int in cache:
        return cache[key_int]

    by_chrom = defaultdict(list)
    paf_path = f'{prefix}/21_pat_depth/{key_int}.paf'
    if not os.path.isfile(paf_path):
        cache[key_int] = by_chrom
        return by_chrom

    with open(paf_path, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) <= CHR_END:
                continue
            try:
                st = int(fields[CHR_STR])
                nd = int(fields[CHR_END])
            except ValueError:
                continue
            if nd < st:
                continue
            by_chrom[fields[CHR_NAM]].append((st, nd))

    cache[key_int] = by_chrom
    return by_chrom

def merge_ref_intervals(intervals):
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for st, nd in intervals[1:]:
        if st <= merged[-1][1]:
            if nd > merged[-1][1]:
                merged[-1][1] = nd
        else:
            merged.append([st, nd])
    return [tuple(x) for x in merged]

def interval_strictly_contains_any(merged_intervals, st, nd):
    return any(intv_st < st and nd < intv_nd for intv_st, intv_nd in merged_intervals)

def raw_true_value(value):
    return value is True or value == 1 or value == 'True' or value == 'true'

def record_has_any_point_no_span(record):
    raw_no_span = record.get('raw_point_no_spanning')
    if isinstance(raw_no_span, dict):
        return raw_true_value(raw_no_span.get('point_a')) or raw_true_value(raw_no_span.get('point_b'))
    side_records = record.get('side_records', [])
    return any(
        raw_true_value(side.get('raw_point_no_spanning', side.get('no_spanning_rawread')))
        for side in side_records
    )

def side_inner_interval(side):
    return (
        int(side.get('inner_st', side.get('path_drop_st', 0))),
        int(side.get('inner_nd', side.get('path_drop_nd', 0))),
    )

def record_has_same_chrom_contiguous_span_path(record, prefix, weights, meandepth, min_depth_N):
    if weights is None:
        return False
    chrom_pair = record.get('chrom_pair', ())
    side_records = record.get('side_records', [])
    if len(side_records) < 2:
        return False
    chrom_a = chrom_pair[0] if len(chrom_pair) > 0 else side_records[0].get('chrom')
    chrom_b = chrom_pair[1] if len(chrom_pair) > 1 else side_records[1].get('chrom')
    if chrom_a != chrom_b:
        return False
    if not record_has_any_point_no_span(record):
        return False

    a_st, a_nd = side_inner_interval(side_records[0])
    b_st, b_nd = side_inner_interval(side_records[1])
    span_st = min(a_st, b_st)
    span_nd = max(a_nd, b_nd)
    if span_nd <= span_st:
        return False

    path_records = globals().get('paf_ans_list', [])
    if not path_records:
        return False
    min_weight = float(min_depth_N) * float(meandepth) / 2.0
    component_cache = {}
    for col_idx, (_, key_int_list) in enumerate(path_records):
        if col_idx >= len(weights) or float(weights[col_idx]) <= min_weight:
            continue
        intervals = []
        for key_int in key_int_list:
            intervals.extend(read_component_ref_intervals(prefix, key_int, component_cache).get(chrom_a, []))
        if interval_strictly_contains_any(merge_ref_intervals(intervals), span_st, span_nd):
            return True
    return False

def build_virtual_inv_events(prefix, meandepth, chrom_lengths, min_depth_N=0.0, weights=None):
    result_path = f'{prefix}/{RAW_TRANSLOCATION_RESULT_PKL}'
    if not os.path.isfile(result_path) or meandepth <= 0:
        return []

    report_by_pair = load_raw_translocation_report(prefix)
    with open(result_path, 'rb') as f:
        records = pkl.load(f)

    display_inv = []
    for record in records:
        try:
            pair_id = int(record.get('pair_id'))
        except (TypeError, ValueError):
            continue

        report_row = report_by_pair.get(pair_id)
        if not record_depth_is_balanced(record, report_row):
            continue
        if not record_passes_virtual_inv_vaf(record, report_row):
            continue
        if not (
            record_has_both_point_spans(record) or
            record_has_same_chrom_contiguous_span_path(record, prefix, weights, meandepth, min_depth_N)
        ):
            continue

        expected_depth = estimate_raw_virtual_inv_depth(record, report_row)
        if expected_depth is None:
            continue
        depth_N = expected_depth / meandepth * 2
        if depth_N < min_depth_N:
            continue

        chrom_a, point_a, chrom_b, point_b = record_display_points(record, report_row)
        if chrom_a not in chrom_lengths or chrom_b not in chrom_lengths:
            continue
        if point_a is None or point_b is None:
            continue

        if chrom_a == chrom_b:
            st, nd = sorted([int(point_a), int(point_b)])
            st = max(0, min(st, int(chrom_lengths[chrom_a])))
            nd = max(0, min(nd, int(chrom_lengths[chrom_a])))
            if nd > st:
                display_inv.append(('v', st, nd, depth_N, chrom_a, f'RAW_TRANSLOCATION_PAIR_{pair_id}'))
            continue

        side_a = layout_side(record, 'layout_b', 0)
        side_b = layout_side(record, 'layout_a', 1)
        interval_a = point_to_chrom_end_interval(chrom_a, point_a, side_a, chrom_lengths)
        interval_b = point_to_chrom_end_interval(chrom_b, point_b, side_b, chrom_lengths)
        if interval_a is not None:
            display_inv.append(('v', interval_a[0], interval_a[1], depth_N, chrom_a, f'RAW_TRANSLOCATION_PAIR_{pair_id}_A'))
        if interval_b is not None:
            display_inv.append(('v', interval_b[0], interval_b[1], depth_N, chrom_b, f'RAW_TRANSLOCATION_PAIR_{pair_id}_B'))

    return display_inv

def type4_indel_graph_source_label(event, event_key):
    contig_name = event.get("contig_name")
    if contig_name:
        return f"TYPE4_INDEL_GRAPH_{contig_name}"
    type4_tuple = event.get("type4_tuple")
    if type4_tuple:
        return "TYPE4_INDEL_GRAPH_" + "_".join(map(str, type4_tuple))
    if isinstance(event_key, tuple):
        return "TYPE4_INDEL_GRAPH_" + "_".join(map(str, event_key))
    return f"TYPE4_INDEL_GRAPH_{event_key}"


def get_path_type4_indel_events(path, type4_event_by_key, type4_path_event_usage):
    return [
        type4_event_by_key[event_key]
        for event_key in type4_path_event_usage.get(path, {})
        if event_key in type4_event_by_key
    ]


def type4_indel_karyotype_labels(type4_indel_events):
    labels = []
    seen = set()
    for event in type4_indel_events or []:
        if event.get("span_len", 0) < KARYOTYPE_MIN_SEGMENT_LENGTH:
            continue
        key = (
            event.get("event_type"),
            event.get("chrom"),
            event.get("st"),
            event.get("nd"),
        )
        if key in seen:
            continue
        seen.add(key)
        labels.append(
            f"{virtual_event_label(event['event_type'])}({chrom_to_iscn(event['chrom'])})"
        )
    return labels


def append_extra_karyotype_labels(base_iscn, extra_labels):
    remaining = list(extra_labels)
    if base_iscn in remaining:
        remaining.remove(base_iscn)
    return base_iscn + ''.join(remaining)


def get_type4_indel_boundary_labels(path_path, type4_edge_to_event_key,
                                    type4_event_by_key, maxh):
    if not type4_edge_to_event_key or not type4_event_by_key:
        return []

    path = import_index_path(path_path)
    ctg_intype_key_by_edge = {}
    for key_int in path2key_int_list.get(path_path, []):
        key_type, key_value = int2key[key_int]
        if key_type == CTG_IN_TYPE:
            ctg_intype_key_by_edge[key_value] = key_int

    if len(path[0]) < 4:
        path[0] = tuple([0] + list(path[0]))
    if len(path[-1]) < 4:
        path[-1] = tuple([0] + list(path[-1]))

    curr_incr = '+' if path[0][NODE_NAME][-1] == 'f' else '-'
    first_real_node = ppc_data[path[1][NODE_NAME]]
    curr_chr = [first_real_node[CHR_NAM], curr_incr]
    curr_ref = first_real_node[CHR_STR] if curr_incr == '+' else first_real_node[CHR_END]
    current_y_bp = 0
    labels = []
    seen = set()

    for i in range(1, len(path) - 1):
        prev_node_name = path[i - 1][NODE_NAME]
        curr_node_name = path[i][NODE_NAME]
        if not isinstance(prev_node_name, int) or not isinstance(curr_node_name, int):
            continue

        last_node = ppc_data[prev_node_name]
        curr_node = ppc_data[curr_node_name]
        edge_key = (
            (path[i - 1][0], prev_node_name),
            (path[i][0], curr_node_name),
        )
        type4_event_key = type4_edge_to_event_key.get(edge_key)
        type4_event = type4_event_by_key.get(type4_event_key) if type4_event_key is not None else None
        type4_deletion_edge = (
            type4_event is not None and type4_event.get("event_type") == "d"
        )

        ctg_intype_key_int = ctg_intype_key_by_edge.get(edge_key)
        interrupt_pieces = []
        if ctg_intype_key_int is not None:
            endpoint_chroms = {last_node[CHR_NAM], curr_node[CHR_NAM]}
            interrupt_pieces = get_ctg_intype_interrupt_pieces(ctg_intype_key_int, endpoint_chroms)

        if not (
            path[i][CHR_CHANGE_IDX] > path[i - 1][CHR_CHANGE_IDX]
            or path[i][DIR_CHANGE_IDX] > path[i - 1][DIR_CHANGE_IDX]
            or interrupt_pieces
            or type4_deletion_edge
        ):
            continue

        if curr_incr == '+':
            current_y_bp += abs(last_node[CHR_END] - curr_ref)
        else:
            current_y_bp += abs(curr_ref - last_node[CHR_STR])

        if type4_deletion_edge:
            label_key = (
                type4_event.get("event_type"),
                type4_event.get("chrom"),
                type4_event.get("st"),
                type4_event.get("nd"),
            )
            if label_key not in seen:
                seen.add(label_key)
                labels.append((
                    current_y_bp / maxh * 100,
                    f"{virtual_event_label(type4_event['event_type'])}({chrom_to_iscn(type4_event['chrom'])})",
                ))

        for _piece_chr, piece_length in interrupt_pieces:
            current_y_bp += piece_length

        if path[i][NODE_NAME] > path[i - 1][NODE_NAME]:
            curr_incr = curr_node[CTG_DIR]
            curr_chr = [curr_node[CHR_NAM], curr_incr]
            curr_ref = curr_node[CHR_STR] if curr_incr == '+' else curr_node[CHR_END]
        else:
            curr_incr = '-' if curr_node[CTG_DIR] == '+' else '+'
            curr_chr = [curr_node[CHR_NAM], curr_incr]
            curr_ref = curr_node[CHR_STR] if curr_incr == '+' else curr_node[CHR_END]

    return labels


def type4_indel_graph_weights_for_paths(paths, loc2weight, type4_path_event_usage):
    weights_by_event = Counter()
    for path in paths:
        path_weight = float(loc2weight.get(path, 0.0))
        if path_weight <= 0:
            continue
        for event_key, count in type4_path_event_usage.get(path, {}).items():
            weights_by_event[event_key] += path_weight * count
    return weights_by_event


def build_aggregated_indel_events(weights, min_weight=0.0, min_span=0,
                                  type4_usage_data=None,
                                  type4_excluded_weights=None):
    indel_events = {}
    rpll = len(paf_ans_list)

    for i in range(rpll, min(len(weights), len(tot_loc_list))):
        paf_loc = tot_loc_list[i]
        if paf_loc.split('/')[-3] == '12_cent_fragment':
            continue

        indel_ind = paf_loc.split('/')[-2]
        if indel_ind not in {'front_jump', 'back_jump'}:
            continue

        with open(paf_loc, "r") as f:
            l = f.readline().rstrip().split("\t")
            chrom = l[CHR_NAM]
            pos1 = int(l[CHR_STR])
            pos2 = int(l[CHR_END])

        event_type = 'd' if indel_ind == 'front_jump' else 'i'
        add_weighted_indel_event(
            indel_events, event_type, chrom, pos1, pos2, float(weights[i]),
            source=f"INDEL_INDEX_{i - rpll}",
        )

    if type4_usage_data is None:
        type4_usage_data = summarize_type4_indel_graph_usage(
            PREFIX, paf_ans_list, weights, import_index_path
        )
    type4_events, type4_weights, _ = type4_usage_data
    if type4_excluded_weights is None:
        type4_excluded_weights = {}
    for event_key, raw_weight in type4_weights.items():
        raw_weight = float(raw_weight) - float(type4_excluded_weights.get(event_key, 0.0))
        if raw_weight <= 0:
            continue
        event = type4_events[event_key]
        add_weighted_indel_event(
            indel_events,
            event["event_type"],
            event["chrom"],
            event["st"],
            event["nd"],
            float(raw_weight),
            source=type4_indel_graph_source_label(event, event_key),
        )

    return [
        event for event in indel_events.values()
        if event["weight"] > min_weight and (event["nd"] - event["st"]) > min_span
    ]


def karyotype_path_to_iscn(pieces : list, type4_indel_events=None):
    """get_karyotype_summary 가 만든 pieces([((chrom, strand), length_bp), ...]) 를
    ISCN 문자열로 변환한다. KARYOTYPE_MIN_SEGMENT_LENGTH(1Mb) 미만 segment 는 무시.
    단, CTG_IN_TYPE 내부에서 끼어든 것으로 남긴 중간 chr 조각은 10Kb 이상이면 유지.
      - segment 1개(정상 단일 염색체): '1'
      - junction 마다: 다른 염색체면 t(a;b), 같은 염색체 방향(strand)전환이면 inv(a)
    표기할 segment 가 하나도 없으면 None 반환.
    """
    filtered = []
    for i, (seg_chr, length) in enumerate(pieces):
        keep = length >= KARYOTYPE_MIN_SEGMENT_LENGTH
        if not keep and length >= CTG_INTYPE_INSERT_MIN_SEGMENT_LENGTH and 0 < i < len(pieces) - 1:
            prev_chrom = pieces[i - 1][0][0]
            curr_chrom = seg_chr[0]
            next_chrom = pieces[i + 1][0][0]
            keep = curr_chrom != prev_chrom and curr_chrom != next_chrom
        if keep:
            filtered.append((seg_chr, length))
    extra_labels = type4_indel_karyotype_labels(type4_indel_events)
    if not filtered:
        return ''.join(extra_labels) if extra_labels else None
    tokens = []
    for i in range(1, len(filtered)):
        prev_chrom, prev_strand = filtered[i - 1][0]
        curr_chrom, curr_strand = filtered[i][0]
        if prev_chrom != curr_chrom:
            tokens.append(f"t({chrom_to_iscn(prev_chrom)};{chrom_to_iscn(curr_chrom)})")
        elif prev_strand != curr_strand:
            tokens.append(f"inv({chrom_to_iscn(curr_chrom)})")
        # 같은 염색체·같은 strand 인접(작은 segment 제거로 병합)은 junction 아님 -> 생략
    if tokens:
        base_iscn = ''.join(tokens)
        return append_extra_karyotype_labels(base_iscn, extra_labels)
    # junction 이 없는 단일 염색체: reference 길이의 90%(KARYOTYPE_NORMAL_RATIO) 미만이면
    # 일부 결실로 보아 del 로 취급, 그렇지 않으면 정상('1')
    chrom = filtered[0][0][0]
    total_len = sum(length for _, length in filtered)
    ratio = total_len / chr_len[chrom] if chrom in chr_len else None
    if ratio is not None:
        logging.debug(f"single-chrom path {chrom}: len={total_len} ref_ratio={ratio:.3f}")
    if ratio is not None and ratio < KARYOTYPE_NORMAL_RATIO:
        base_iscn = f"del({chrom_to_iscn(chrom)})"
        return append_extra_karyotype_labels(base_iscn, extra_labels)
    base_iscn = chrom_to_iscn(chrom)
    if extra_labels:
        return ''.join(extra_labels)
    return base_iscn


def build_reference_cluster_figure_modes() -> tuple:
    """Validate the requested reference-cluster figure layouts."""

    modes = tuple(
        mode for mode in REFERENCE_PATH_CLUSTER_FIGURE_MODES if mode in DEPTH_MODES
    )
    unknown = set(REFERENCE_PATH_CLUSTER_FIGURE_MODES) - set(DEPTH_MODES)
    if unknown:
        logging.warning(
            'Ignoring unknown SKYPE_REFSPAN_MODES values: %s (known: %s)',
            ','.join(sorted(unknown)), ','.join(DEPTH_MODES),
        )
    return modes or ('plain',)


def build_reference_cluster_depth_tracks(
    *, clusters, path_spans, path_metadata, fig_prefix: str
) -> dict:
    """Attach observed, fitted, and per-path depth to each reference cluster.

    Member paths reuse the same separated depth pieces the NNLS matrix was
    built from, so a route's contribution is read at exactly the windows the
    solver fitted rather than being re-estimated from the karyotype.
    """

    bin_index = build_depth_bin_index(df)
    observed_cn = observed_copy_number(df, meandepth)
    model_cn = load_model_copy_number(PREFIX, fig_prefix, bin_index, meandepth)
    profiles = PathDepthProfiles(PREFIX, path2key_int_list, bin_index)
    tracks = {}
    for cluster in clusters:
        try:
            tracks[cluster.cluster_id] = build_cluster_depth_tracks(
                members=cluster.members,
                path_metadata=path_metadata,
                path_spans=path_spans,
                bin_index=bin_index,
                observed_cn=observed_cn,
                model_cn=model_cn,
                profiles=profiles,
            )
        except KeyError as error:
            logging.warning(
                'Reference cluster %s has no depth track: %s',
                cluster.cluster_id, error,
            )
    if profiles.missing_pieces:
        logging.warning(
            'Reference cluster depth: %d separated depth files missing (%s...)',
            len(profiles.missing_pieces),
            ','.join(str(piece) for piece in profiles.missing_pieces[:5]),
        )
    return tracks


def build_karyotype_diagram(fig_prefix : str = '', filter_depth_N : float = TARGET_DEPTH):
    weights = np.load(f'{PREFIX}/weight{fig_prefix}.npy')
    loc2weight = dict(zip(tot_loc_list, weights))
    weights_sorted_data = sorted(enumerate(weights), key=lambda t:t[1], reverse=True)
    type4_usage_data = summarize_type4_indel_graph_usage(
        PREFIX, paf_ans_list, weights, import_index_path
    )
    type4_event_by_key, _, type4_path_event_usage = type4_usage_data
    _type4_loaded_events, type4_edge_to_event_key = load_type4_indel_graph_event_index(PREFIX)
    
    non_type4_top_path = []
    for ind, w in weights_sorted_data:
        paf_loc = tot_loc_list[ind]
        key = paf_loc.split('/')[-3]
        if key not in {'11_ref_ratio_outliers', '12_cent_fragment'}:
            if w > filter_depth_N * meandepth / 2:
                non_type4_top_path.append(paf_loc)

    nclose_event_catalog = load_event_catalog(PREFIX)
    junction_event_keys = bnd_event_keys(nclose_event_catalog) | {
        event['event_key']
        for event in nclose_event_catalog
        if event.get('kind') == 'indel'
        and event.get('graph_only')
        and len(event.get('type4_tuple', ())) == 2
    }
    junction_event_keys.update(
        tuple(event['type4_tuple'])
        for event in nclose_event_catalog
        if event.get('kind') == 'indel'
        and event.get('graph_only')
        and len(event.get('type4_tuple', ())) == 2
    )
    path_nclose_boundaries = {}
    karyotypes_data = get_karyotype_summary(
        non_type4_top_path, type4_edge_to_event_key, type4_event_by_key,
        junction_event_keys=junction_event_keys,
        path_nclose_boundaries=path_nclose_boundaries,
    )
    shown_type4_indel_weights = type4_indel_graph_weights_for_paths(
        karyotypes_data.keys(), loc2weight, type4_path_event_usage
    )

    all_ecdna_path = []
    ecdna_path_dict = {}

    for ind, w in weights_sorted_data:
        paf_loc = tot_loc_list[ind]
        key = paf_loc.split('/')[-3]
        indel_ind = paf_loc.split('/')[-2]
        if key == '11_ref_ratio_outliers':
            if indel_ind == 'ecdna' and w > NCLOSE_SIM_DIFF_THRESHOLD:
                ecdna_path_dict[paf_loc] = w
                all_ecdna_path.append(paf_loc)

    # long_ecdna_path = {}
    # for i in all_ecdna_path:
    #     min_pos = 1e9
    #     max_pos = -1
    #     chr_nam = ''
    #     depth = ecdna_path_dict[i]/meandepth * 2
    #     with open(i, "r") as f:
    #         for line in f:
    #             line = line.rstrip()
    #             line = line.split("\t")
    #             if chr_nam == '':
    #                 chr_nam = line[CHR_NAM]
    #             pos1 = int(line[CHR_STR])
    #             pos2 = int(line[CHR_END])
    #             min_pos = min(min_pos, pos1, pos2)
    #             max_pos = max(max_pos, pos1, pos2)
    #     long_ecdna_path[i] = (chr_nam, min_pos, max_pos, depth)

    display_indel = defaultdict(list)

    for event in build_aggregated_indel_events(
        weights, NCLOSE_SIM_DIFF_THRESHOLD, TYPE4_CLUSTER_SIZE, type4_usage_data,
        shown_type4_indel_weights
    ):
        chrom = event["chrom"]
        display_indel[chrom].append((
            event["event_type"],
            event["st"],
            event["nd"],
            event["weight"] / meandepth * 2,
            chrom,
            indel_event_source_label(event),
        ))

    virtual_inv_display = build_virtual_inv_events(
        PREFIX, meandepth, chr_len, min_depth_N=filter_depth_N, weights=weights
    )
        
    loc2weight = dict(zip(tot_loc_list, weights))
    
    maxh = max(chr_len.values())
    for i in karyotypes_data.values():
        h = 0
        for j in i:
            h += j[1]
        maxh = max(maxh, h)

    karyotypes_norm_data = dict()
    for path, i in karyotypes_data.items():
        temp_list = []
        for j in i:
            temp_list.append((j[0], j[1] / maxh * 100))
        karyotypes_norm_data[path] = temp_list

    grouped_norm_data = defaultdict(list)
    
    for path, data in karyotypes_norm_data.items():
        cnt = Counter()
        for c, w in data:
            cnt[c] += w
        sorted_cnt_data = sorted(cnt.items(), key=lambda t: -t[1])
        grouped_norm_data[sorted_cnt_data[0][0][0]].append((path, data))

    fragment_display = []
    for ind, w in weights_sorted_data:
        paf_loc = tot_loc_list[ind]
        if paf_loc.split('/')[-3] == '12_cent_fragment':
            if w > filter_depth_N * meandepth / 2:
                chrom = paf_loc.split('/')[-2]
                side = paf_loc.split('/')[-1].split('.')[0]
                info = cen_fragment_meta[chrom]
                fragment_display.append((chrom, side, info["mid"], info["chr_len"], w / meandepth * 2))

    path_depth_n = {
        path: float(loc2weight[path] / meandepth * 2)
        for path in karyotypes_data
    }
    path_karyotype = {}
    path_event_labels = {}
    for path, pieces in karyotypes_data.items():
        path_type4_events = get_path_type4_indel_events(
            path, type4_event_by_key, type4_path_event_usage
        )
        path_karyotype[path] = karyotype_path_to_iscn(pieces, path_type4_events)
        path_event_labels[path] = get_type4_indel_boundary_labels(
            path, type4_edge_to_event_key, type4_event_by_key, maxh
        )

    path_junction_read_counts = {}
    path_junction_read_depths = {}
    path_depth_raw = None
    if fig_prefix == '_cluster':
        cluster_nclose_read_evidence = load_cluster_nclose_read_evidence(PREFIX)
        if cluster_nclose_read_evidence is not None:
            path_junction_read_depths = {
                path: [
                    (
                        boundary_index,
                        (
                            cluster_nclose_read_evidence.get(event_key, {})
                            .get('vaf_read_depth')
                        ),
                    )
                    for boundary_index, event_key in path_nclose_boundaries.get(path, [])
                ]
                for path in karyotypes_data
            }
        path_depth_raw = {
            path: float(loc2weight[path])
            for path in karyotypes_data
        }

    render_result = render_virtual_sky(
        output_prefix=PREFIX,
        cell_line=CELL_LINE,
        chromosome_lengths=chr_len,
        grouped_norm_data=grouped_norm_data,
        display_indel=display_indel,
        virtual_inv_display=virtual_inv_display,
        fragment_display=fragment_display,
        maxh=maxh,
        path_depth_n=path_depth_n,
        path_karyotype=path_karyotype,
        path_event_labels=path_event_labels,
        path_junction_read_counts=path_junction_read_counts,
        path_junction_read_depths=path_junction_read_depths,
        path_depth_raw=path_depth_raw,
        fig_prefix=fig_prefix,
    )

    if fig_prefix == '_cluster':
        ordinary_column_by_path = {
            path: column
            for column, (path, _key_int_list) in enumerate(paf_ans_list)
        }
        path_usage = load_path_usage(PREFIX, expected_len=len(weights))
        type4_deletion_edges = {
            edge
            for edge, event_key in type4_edge_to_event_key.items()
            if type4_event_by_key.get(event_key, {}).get('event_type') == 'd'
        }
        candidate_spans = {}
        path_cluster_metadata = {}
        excluded_full_reference = []
        for path in karyotypes_data:
            column = ordinary_column_by_path.get(path)
            if column is None:
                raise ValueError(
                    f'Virtual SKY ordinary path is missing its matrix column: {path}'
                )
            spans = extract_ordered_reference_spans(
                import_index_path(path),
                ppc_data,
                forced_break_edges=type4_deletion_edges,
                interrupt_spans_by_edge=get_path_ctg_intype_interrupt_spans(path),
            )
            nclose_count = sum(int(count) for count in path_usage[column].values())
            full_reference, full_reference_chrom, full_reference_ratio = (
                is_full_reference_path(
                    spans,
                    has_nclose=bool(path_usage[column]),
                    chromosome_lengths=chr_len,
                    minimum_ratio=KARYOTYPE_NORMAL_RATIO,
                )
            )
            metadata = {
                'path_column': column,
                'location': path,
                'raw_depth': float(loc2weight[path]),
                'depth_n': float(path_depth_n[path]),
                'nclose_count': nclose_count,
            }
            if full_reference:
                excluded_full_reference.append({
                    **metadata,
                    'full_reference_chrom': full_reference_chrom,
                    'full_reference_ratio': full_reference_ratio,
                })
                continue
            candidate_spans[column] = spans
            path_cluster_metadata[column] = metadata

        reference_clusters, reference_singletons = cluster_reference_paths(
            candidate_spans,
            minimum_overlap=REFERENCE_PATH_CLUSTER_MIN_OVERLAP,
        )
        figure_modes = build_reference_cluster_figure_modes()
        cluster_tracks = {}
        if reference_clusters and set(figure_modes) - {'plain'}:
            cluster_tracks = build_reference_cluster_depth_tracks(
                clusters=reference_clusters,
                path_spans=candidate_spans,
                path_metadata=path_cluster_metadata,
                fig_prefix=fig_prefix,
            )
            if not cluster_tracks:
                figure_modes = ('plain',)
        reference_cluster_result = render_reference_path_clusters(
            output_prefix=PREFIX,
            cell_line=CELL_LINE,
            clusters=reference_clusters,
            singletons=reference_singletons,
            path_spans=candidate_spans,
            path_metadata=path_cluster_metadata,
            excluded_full_reference=excluded_full_reference,
            minimum_overlap=REFERENCE_PATH_CLUSTER_MIN_OVERLAP,
            full_reference_ratio=KARYOTYPE_NORMAL_RATIO,
            modes=figure_modes,
            cluster_tracks=cluster_tracks,
        )
        logging.info(
            'Reference-overlap path clusters: %d clusters, %d singleton paths, '
            '%d full-reference paths excluded',
            reference_cluster_result['cluster_count'],
            reference_cluster_result['singleton_count'],
            reference_cluster_result['excluded_full_reference_count'],
        )
        render_result['reference_path_clusters'] = reference_cluster_result

    return render_result


parser = argparse.ArgumentParser(description="SKYPE depth analysis")

parser.add_argument("ppc_paf_file_path", 
                    help="Path to the preprocessed PAF file.")

parser.add_argument("main_stat_loc", 
                    help="Cancer coverage location file")

parser.add_argument("telomere_bed_path", 
                    help="Path to the telomere information file.")

parser.add_argument("reference_fai_path", 
                    help="Path to the chromosome information file.")

parser.add_argument("prefix", 
                    help="Pefix for pipeline")

parser.add_argument("cell_line_name", 
                    help="Path to the cytoband information file.")

args = parser.parse_args()

# t = """
# 30_virtual_sky.py /home/hyunwoo/ACCtools-pipeline/90_skype_run/U2OS/20_alignasm/U2OS.ctg.aln.paf.ppc.paf /home/hyunwoo/ACCtools-pipeline/90_skype_run/U2OS/01_depth/U2OS_normalized.win.stat.gz public_data/chm13v2.0_telomere.bed public_data/chm13v2.0.fa.fai 30_skype_pipe/U2OS_20_58_22 U2OS
# """
# args = parser.parse_args(t.strip().split()[1:])

PREFIX = args.prefix
CHROMOSOME_INFO_FILE_PATH = args.reference_fai_path
main_stat_loc = args.main_stat_loc
TELOMERE_INFO_FILE_PATH = args.telomere_bed_path
PREPROCESSED_PAF_FILE_PATH = args.ppc_paf_file_path
CELL_LINE = args.cell_line_name
pipeline_mode_config = load_pipeline_mode(PREFIX)

df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
df = df.query('chr != "chrM"')
meandepth = np.median(df['meandepth'])
chr_len = find_chr_len(CHROMOSOME_INFO_FILE_PATH)

RATIO_OUTLIER_FOLDER = f"{PREFIX}/11_ref_ratio_outliers/"
front_contig_path = RATIO_OUTLIER_FOLDER+"front_jump/"
back_contig_path = RATIO_OUTLIER_FOLDER+"back_jump/"
ecdna_contig_path = RATIO_OUTLIER_FOLDER+"ecdna/"
TELO_CONNECT_NODES_INFO_PATH = PREFIX+"/telomere_connected_list.txt"
output_folder = f'{PREFIX}/21_pat_depth'
NCLOSE_FILE_PATH = f"{args.prefix}/nclose_nodes_index.txt"

ppc_data = import_ppc_data(PREPROCESSED_PAF_FILE_PATH)

with open(f'{PREFIX}/path_data.pkl', 'rb') as f:
    path_list_dict = pkl.load(f)

fclen = len(glob.glob(front_contig_path+"*"))
bclen = len(glob.glob(back_contig_path+"*"))
eclen = len(glob.glob(ecdna_contig_path+"*"))

tot_loc_list = []

with open(f'{PREFIX}/contig_pat_vec_data.pkl', 'rb') as f:
    paf_ans_list, key_list, int2key, _ = pkl.load(f)
path2key_int_list = dict(paf_ans_list)

for loc, ll in paf_ans_list:
    tot_loc_list.append(loc)

for i in range(1, fclen//4 + 1):
    bv_paf_loc = front_contig_path+f"{i}_base.paf"
    tot_loc_list.append(bv_paf_loc)

for i in range(1, bclen//4 + 1):
    bv_paf_loc = back_contig_path+f"{i}_base.paf"
    tot_loc_list.append(bv_paf_loc)

for i in range(1, eclen//2 + 1):
    ec_paf_loc = ecdna_contig_path+f"{i}.paf"
    tot_loc_list.append(ec_paf_loc)

with open(f'{PREFIX}/cen_fragment_data.pkl', 'rb') as f:
    cen_fragment_meta = pkl.load(f)
cen_fragment_list = sorted(cen_fragment_meta.items(), key=lambda kv: chr2int(kv[0]))
for chrom, info in cen_fragment_list:
    side = 'right' if info["dir"] else 'left'
    tot_loc_list.append(f'{PREFIX}/12_cent_fragment/{chrom}/{side}.fragment')

with open(f'{PREFIX}/tot_loc_list.pkl', 'wb') as f:
    pkl.dump(tot_loc_list, f)

tot_loc_list2nclosecnt = dict()
for paf_loc in tot_loc_list:
    if paf_loc.split('/')[-3] not in {'11_ref_ratio_outliers', '12_cent_fragment'}:
        path = import_index_path(paf_loc)

        if len(path[0]) < 4:
            path[0] = tuple([0] + list(path[0])) # padding for easier calculation
        if len(path[-1]) < 4:
            path[-1] = tuple([0] + list(path[-1]))

        nclose_use_cnt = 0
        for i in range(1, len(path)-1):
            if path[i][CHR_CHANGE_IDX] > path[i-1][CHR_CHANGE_IDX] \
            or path[i][DIR_CHANGE_IDX] > path[i-1][DIR_CHANGE_IDX]:
                nclose_use_cnt += 1

        tot_loc_list2nclosecnt[paf_loc] = nclose_use_cnt

telo_data = import_telo_data(TELOMERE_INFO_FILE_PATH, chr_len)
telo_connected_node = extract_telomere_connect_contig(TELO_CONNECT_NODES_INFO_PATH)
telo_dict = defaultdict(list)

for _ in telo_data:
    telo_dict[_[0]].append(_[1:])

telo_fb_dict = defaultdict(list)
for k, v in telo_dict.items():
    for i in v:
        telo_fb_dict[k+i[-1]].append([i[0], i[1]])

chr_inf = max(chr_len.values())
chr_fb_len_dict = defaultdict(list)

for chr_dir, node_id in telo_connected_node:
    telo_len = chr_inf
    for telo_bed in telo_fb_dict[chr_dir]:
        telo_len = min(telo_len, distance_checker(tuple(telo_bed), (ppc_data[node_id][CHR_STR], ppc_data[node_id][CHR_END])))
    chr_fb_len_dict[chr_dir].append((node_id, telo_len, chr_dir))

telo_len_data = []
for chr_dir, telo_len_list in chr_fb_len_dict.items():
    s_telo_len_list = sorted(telo_len_list, key=lambda t: t[1])
    telo_len_data.extend(filter(lambda t: t[1] > 0, s_telo_len_list[1:]))

need_label = defaultdict(list)
need_label_index = dict()
for node_id, telo_len, chr_dir in telo_len_data:
    need_label[chr_dir[:-1]].append((node_id, chr_dir[-1]))
    need_label_index[node_id] = (chr_dir, telo_len)

weights = np.load(f'{PREFIX}/weight.npy')

build_karyotype_diagram()

use_julia_solver = pipeline_mode_is_karyotype(pipeline_mode_config)

if use_julia_solver:
    build_karyotype_diagram(fig_prefix='_filter')
    build_karyotype_diagram(fig_prefix='_cluster')
