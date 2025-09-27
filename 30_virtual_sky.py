import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle as pkl
import matplotlib.patches as patches

import os
import re
import ast
import glob
import logging
import argparse

from matplotlib.ticker import MultipleLocator
from collections import defaultdict
from collections import Counter

logging.basicConfig(
    format='%(asctime)s %(levelname)s:%(message)s',
    level=logging.INFO,
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
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
K = 1000
M = K * 1000
MAX_PATH_CNT = 100
INF = 1000000000

CHROMOSOME_COUNT = 23
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

JOIN_BASELINE = 0.8
KARYOTYPE_SECTION_MINIMUM_LENGTH = 100 * K

HARD_PATH_COUNT_BASELINE = 100 * K

NODE_NAME = 1
CHR_CHANGE_IDX = 2
DIR_CHANGE_IDX = 3

CHR_COLORS = {
    "chr1":  "#FFC000",
    "chr2":  "#FF0000",
    "chr3":  "#00B050",
    "chr4":  "#0070C0",
    "chr5":  "#7030A0",
    "chr6":  "#92D050",
    "chr7":  "#C00000",
    "chr8":  "#00B0F0",
    "chr9":  "#FFFF00",
    "chr10": "#F08080",
    "chr11": "#48D1CC",
    "chr12": "#D2B48C",
    "chr13": "#B0C4DE",
    "chr14": "#FFA07A",
    "chr15": "#CD5C5C",
    "chr16": "#E9967A",
    "chr17": "#BDB76B",
    "chr18": "#D8BFD8",
    "chr19": "#BC8F8F",
    "chr20": "#FFD700",
    "chr21": "#7FFFD4",
    "chr22": "#ADFF2F",
    "chrX":  "#6495ED",
    "chrY":  "#EE82EE",
}

def similar_check(v1, v2, ratio):
    try:
        assert(v1 >= 0 and v2 >= 0)
    except:
        logging.error(f"Invalid values for similarity check: v1={v1}, v2={v2}")
        assert(False)
    mi, ma = sorted([v1, v2])
    return False if mi == 0 else (ma / mi <= ratio) or ma-mi < NCLOSE_SIM_DIFF_THRESHOLD

def check_near_bnd(chrom, inside_st, inside_nd, ratio=NCLOSE_SIM_COMPARE_RATIO):
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

def import_telo_data(file_path : str, chr_len : dict) -> dict :
    fai_file = open(file_path, "r")
    telo_data = [(0,1)]
    for curr_data in fai_file:
        temp_list = curr_data.split("\t")
        int_induce_idx = [1, 2]
        for i in int_induce_idx:
            temp_list[i] = int(temp_list[i])
        if temp_list[0]!=telo_data[-1][0]:
            temp_list[2]+=TELOMERE_EXPANSION
            temp_list.append('f')
        else:
            if temp_list[1]>chr_len[temp_list[0]]/2:
                temp_list[1]-=TELOMERE_EXPANSION
                temp_list.append('b')
            else:
                temp_list.append('f')
                temp_list[2]+=TELOMERE_EXPANSION
        telo_data.append(tuple(temp_list))
    fai_file.close()
    return telo_data[1:]

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

# SKY Figure
def plot_virtual_chromosome(ax, segments_data, maxh, cell_col, def_cell_col, label=None):
    """
    주어진 segments 리스트로 가상 염색체를 그립니다.
    """
    path, segments = segments_data

    x_center = 0  # 좌우 중앙값 (필요에 따라 조정)
    width = 10     # 염색체 폭

    # 전체 가상 염색체 길이 계산
    total_length = sum(seg_length for (_, seg_length) in segments)
    radius = width / 2.0
    round_radius = min(radius, total_length / 4)

    # 전체 영역을 둥근 모서리로 클리핑할 패치 생성
    clip_patch = patches.FancyBboxPatch(
        (x_center - radius, 0),
        width,
        total_length,
        boxstyle=f"round,pad=0,rounding_size={round_radius}",
        facecolor='none',
        edgecolor='none'
    )
    ax.add_patch(clip_patch)

    # 각 segment를 누적해서 그리며, 경계가 바뀔 때마다 텍스트 추가
    current_y = 0
    text_obj_list = []
    text_x = x_center + radius * 2
    decoration_args = {
        'ha': 'left',
        'va': 'center',
        'fontsize': 10,
        'color': 'black'
    }
    for i, (real_chr, seg_length) in enumerate(segments):
        color = CHR_COLORS.get(real_chr[0], "gray")
        rect = patches.Rectangle(
            (x_center - radius, current_y),
            width,
            seg_length,
            facecolor=color,
            edgecolor='none'
        )
        rect.set_clip_path(clip_patch)
        ax.add_patch(rect)

        if 0 < i and i < len(segments):
            last_chr = segments[i-1][0]
            if real_chr[0] != last_chr[0]:
                text_label = f"t({last_chr[0][3:]};{real_chr[0][3:]})"
            else:
                text_label = f"t({last_chr[0][3:]+str(last_chr[1])};{real_chr[0][3:] + str(real_chr[1])})"
            text_obj = ax.text(text_x, current_y, text_label, **decoration_args)
            text_obj_list.append(text_obj)

        current_y += seg_length
        
    mark_overlapping_texts_with_arrows(ax, text_obj_list, min_gap=5)
    outline_patch = patches.FancyBboxPatch(
        (x_center - radius, 0),
        width,
        total_length,
        boxstyle=f"round,pad=0,rounding_size={round_radius}",
        facecolor='none',
        edgecolor='black',
        linewidth=1.2,
        clip_on=False
    )
    ax.add_patch(outline_patch)


    # 가상 염색체 라벨 (옵션)
    if label:
        ax.text(x_center, -5, label, ha='center', va='top', fontsize=10)

    ax.set_xlim(0, 60 * cell_col / def_cell_col)
    ax.set_ylim(0, 100)
    ax.axis('off')

def plot_indel(ax, indel, maxh, cell_col, def_cell_col, chr_len, label=None):

    x_center = 0  # 좌우 중앙값 (필요에 따라 조정)
    width = 10     # 염색체 
    real_chr = indel[4]
    total_length = chr_len
    radius = width / 2.0
    round_radius = min(radius, total_length / 4)
    chr_color = CHR_COLORS.get(real_chr, "gray")
    background_color = 'white' if indel[0]=='i' else chr_color
    face_color = chr_color if indel[0]=='i' else 'white'

    clip_patch = patches.FancyBboxPatch(
        (x_center - radius, 0),
        width,
        total_length,
        boxstyle=f"round,pad=0,rounding_size={round_radius}",
        facecolor='none',
        edgecolor='none'
    )
    ax.add_patch(clip_patch)

    outline_patch = patches.FancyBboxPatch(
        (x_center - radius, 0),
        width,
        total_length,
        boxstyle=f"round,pad=0,rounding_size={round_radius}",
        facecolor=background_color,
        edgecolor='black',
        linestyle='dashed',
        linewidth=1.2,
        clip_on=False
    )
    ax.add_patch(outline_patch)


    current_y = indel[1] / maxh * 100
    y1 = indel[1] / maxh * 100
    y2 = indel[2] / maxh * 100
    midy = (y1+y2)/2
    seg_length = (indel[2] - indel[1]) / maxh * 100

    indel_rect = patches.Rectangle(
                (x_center - radius, current_y),
                width,
                seg_length,
                facecolor = face_color,
                edgecolor='none'
            )
    
    indel_rect.set_clip_path(clip_patch)
    ax.add_patch(indel_rect)
        
    x_start = x_center + radius
    x_end   = x_center + 2 * radius

    # 위쪽 선
    ax.plot([x_start, x_end],
            [y1,     midy],
            linewidth=1,
            color='black')

    # 아래쪽 선
    ax.plot([x_start, x_end],
            [y2,     midy],
            linewidth=1,
            color='black')
    
    text_obj = ax.text(x_end+radius/5, midy, f"{'del' if indel[0] == 'd' else 'ins'}({indel[-2][3:]})", ha='left', va='center', fontsize=10, color='black')
    
    mark_overlapping_texts_with_arrows(ax, [text_obj], min_gap=5)

    if label:
        ax.text(x_center, -5, label, ha='center', va='top', fontsize=10)

    ax.set_xlim(0, 60 * cell_col / def_cell_col)
    ax.set_ylim(0, 100)
    ax.axis('off')

def plot_scale_bar(ax, chr_name, maxh):
    ax.set_xlim(-5, 5)
    ax.axis(True)

    ax.set_ylim(0, maxh / M)
    ax.set_ylabel("Scale bar (Mb)", rotation=90)

    ax.yaxis.set_major_locator(MultipleLocator(25))  # Major ticks every 25 Mb
    ax.yaxis.set_minor_locator(MultipleLocator(5))   # Minor ticks every 5 Mb
    ax.tick_params(axis='y', which='major', length=5)
    ax.tick_params(axis='y', which='minor', length=3)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.spines['left'].set_position(('data', 0))
    ax.xaxis.set_visible(False)

    if chr_name in chr_len:
        chrom_length = chr_len[chr_name] / M
        line = ax.vlines(x=0.5, ymin=0, ymax=chrom_length, colors=CHR_COLORS[chr_name], linewidth=4, alpha=0.6)

def plot_chr_name(ax, chr_name):
    # ax.axis('on')
    ax.text(0, 0.5, chr_name,
            transform=ax.transAxes,  # 축 비율 기준
            ha='center', va='center',  # 수평/수직 정렬
            fontsize=15, rotation=90)
    
def mark_overlapping_texts_with_arrows(ax, texts, min_gap=5):
    """
    Adjust positions of overlapping text labels on the given axis and draw arrows from the original
    position to the new position if an adjustment is made.
    
    Parameters:
    - ax: matplotlib Axes object where texts are drawn.
    - texts: list of matplotlib.text.Text objects.
    - min_gap: float, minimum vertical gap (in data coordinates) required between texts.
    """
    # texts를 y 좌표 기준으로 정렬 (대부분 동일한 x 좌표라 가정)
    texts_sorted = sorted(texts, key=lambda t: t.get_position()[1])
    
    # 첫 번째 텍스트는 그대로 두고, 이후 text들과 비교하여 겹치는 경우 위치 조정
    for i in range(1, len(texts_sorted)):
        prev_text = texts_sorted[i-1]
        curr_text = texts_sorted[i]
        x, prev_y = prev_text.get_position()
        x, curr_y = curr_text.get_position()
        
        # 만약 이전 text와의 y 간격이 min_gap보다 작으면 조정
        if curr_y - prev_y < min_gap:
            new_y = prev_y + min_gap  # 충분한 간격 확보
            original_y = curr_y       # 원래 위치 저장
            curr_text.set_position((x, new_y))
            # 원래 위치에서 새 위치로 화살표 표시
            ax.annotate("",
                        xy=(x, new_y), xycoords='data',
                        xytext=(5, original_y), textcoords='data',
                        arrowprops=dict(arrowstyle="-|>", color="black", lw=0.5))

def extract_nclose_node(nclose_path: str) -> list:
    nclose_list = []
    with open(nclose_path, "r") as f:
        for line in f:
            line = line.split()
            nclose_list.append((int(line[1]), int(line[2])))
    return nclose_list


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

def get_karyotype_summary(non_type4_path_list : list):
    karyotypes_data_direction_include = {}
    
    for path_path in non_type4_path_list:
        pieces = []
        path = import_index_path(path_path)

        # padding for easier calculation
        if len(path[0]) < 4:
            path[0] = tuple([0] + list(path[0]))
        if len(path[-1]) < 4:
            path[-1] = tuple([0] + list(path[-1]))
        
        curr_ref = 0 if path[0][NODE_NAME][-1] =='f' else chr_len[path[0][NODE_NAME][:-1]]
        curr_incr = '+' if path[0][NODE_NAME][-1] =='f' else '-'
        curr_chr = [path[0][NODE_NAME][:-1], curr_incr]
        
        nclose_use_cnt = 0
        for i in range(1, len(path)-1):
            if path[i][CHR_CHANGE_IDX] > path[i-1][CHR_CHANGE_IDX] \
            or path[i][DIR_CHANGE_IDX] > path[i-1][DIR_CHANGE_IDX]:
                nclose_use_cnt += 1
                last_node = ppc_data[path[i-1][NODE_NAME]]
                curr_node = ppc_data[path[i][NODE_NAME]]
                # add last piece
                if curr_incr == '+':
                    pieces.append((tuple(curr_chr), last_node[CHR_END] - curr_ref))
                else:
                    pieces.append((tuple(curr_chr), curr_ref - last_node[CHR_STR]))
                
                # update info of new piece (starting ref, chromosome type, increment ..)
                if path[i][NODE_NAME] > path[i-1][NODE_NAME]:
                    curr_incr = curr_node[CTG_DIR]
                    curr_chr = [curr_node[CHR_NAM], curr_incr]
                    curr_ref = curr_node[CHR_STR] if curr_incr == '+' else curr_node[CHR_END]
                else:
                    curr_incr = '-' if curr_node[CTG_DIR] == '+' else '+'
                    curr_chr = [curr_node[CHR_NAM], curr_incr]
                    curr_ref = curr_node[CHR_STR] if curr_incr == '+' else curr_node[CHR_END]

        pieces.append((tuple(curr_chr), chr_len[curr_chr[0]] - curr_ref if curr_incr == '+' else curr_ref))
        karyotypes_data_direction_include[path_path] = pieces

    return karyotypes_data_direction_include

def build_karyotype_diagram(weights, fig_prefix : str = '', filter_depth_N : float = 0):
    loc2weight = dict(zip(tot_loc_list, weights))
    weights_sorted_data = sorted(enumerate(weights), key=lambda t:t[1], reverse=True)
    
    non_type4_top_path = []
    for ind, w in weights_sorted_data:
        paf_loc = tot_loc_list[ind]
        key = paf_loc.split('/')[-3]
        if paf_loc.split('/')[-3] != '11_ref_ratio_outliers':
            if w > filter_depth_N * meandepth / 2:
                non_type4_top_path.append(paf_loc)

    karyotypes_data = get_karyotype_summary(non_type4_top_path)

    all_insertion_path = []
    insertion_path_dict = {}

    all_deletion_path = []
    deletion_path_dict = {}

    for ind, w in weights_sorted_data:
        paf_loc = tot_loc_list[ind]
        key = paf_loc.split('/')[-3]
        indel_ind = paf_loc.split('/')[-2]
        if key == '11_ref_ratio_outliers':
            if indel_ind == 'front_jump' and w > NCLOSE_SIM_DIFF_THRESHOLD:
                deletion_path_dict[paf_loc] = w
                all_deletion_path.append(paf_loc)
            if indel_ind == 'back_jump' and w > NCLOSE_SIM_DIFF_THRESHOLD:
                insertion_path_dict[paf_loc] = w
                all_insertion_path.append(paf_loc)

    long_deletion_path = {}
    long_insertion_path = {}

    for i in all_deletion_path:
        with open(i, "r") as f:
            l = f.readline()
            l = l.rstrip()
            l = l.split("\t")
            chr_nam1 = l[CHR_NAM]
            chr_nam2 = l[CHR_NAM]
            pos1 = int(l[CHR_STR])
            pos2 = int(l[CHR_END])
            if abs(pos1-pos2) > TYPE4_CLUSTER_SIZE:
                long_deletion_path[i] = (chr_nam1, pos1, pos2)

    for i in all_insertion_path:
        with open(i, "r") as f:
            l = f.readline()
            l = l.rstrip()
            l = l.split("\t")
            chr_nam1 = l[CHR_NAM]
            chr_nam2 = l[CHR_NAM]
            pos1 = int(l[CHR_STR])
            pos2 = int(l[CHR_END])
            if abs(pos1-pos2) > TYPE4_CLUSTER_SIZE:
                long_insertion_path[i] = (chr_nam1, pos1, pos2)

    display_indel = defaultdict(list)

    for path, check_arg in long_deletion_path.items():
        chrom, pos1, pos2 = check_arg
        # if check_near_type4(chrom, pos1, pos1) or \
        #    check_near_type4(chrom, pos2, pos2):
        display_indel[chrom].append(("d", pos1, pos2, deletion_path_dict[path]/meandepth * 2, chrom, path))

    for path, check_arg in long_insertion_path.items():
        chrom, pos1, pos2 = check_arg
        # if check_near_type4(chrom, pos1, pos1) or \
        #    check_near_type4(chrom, pos2, pos2):
        display_indel[chrom].append(("i", pos1, pos2, insertion_path_dict[path]/meandepth * 2, chrom, path))

        
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

    cols = 10
    rows = 0
    for chr_name, data_list in grouped_norm_data.items():
        chr_indel = display_indel.get(chr_name, [])
        indel_len = len(chr_indel)
        
        rows += ((len(data_list) + indel_len - 1) // cols + 1) if len(data_list) + indel_len > 0 else 0

    cell_col = 3
    cell_row = 4.5
    def_cell_col = 1.8

    prefix_ratios = [0.5, 1]
    ratio_ratios = prefix_ratios + [cell_col] * cols
    height_ratios = [1.5] + [cell_row] * rows

    fig, ax_array = plt.subplots(len(height_ratios), len(ratio_ratios), figsize=(sum(ratio_ratios), sum(height_ratios)), 
                                 gridspec_kw={"width_ratios" : ratio_ratios,
                                            "height_ratios" : height_ratios,
                                            "wspace": 0,
                                            "hspace": 0.2})
                                            
    for ax_list in ax_array:
        for ax in ax_list:
            ax.axis('off')

    sorted_grouped_norm_data_items = sorted(grouped_norm_data.items(), key=lambda t: chr2int(t[0]))

    now_col = 1
    for chr_name, data_list in sorted_grouped_norm_data_items:
        chr_indel = display_indel.get(chr_name, [])
        indel_len = len(chr_indel)

        bef_now_col = now_col
        plot_chr_name(ax_array[now_col][0], chr_name)
        for i, data in enumerate(data_list):
            plot_virtual_chromosome(ax_array[i // cols + now_col][i % cols + len(prefix_ratios)], data, maxh,
                                    cell_col, def_cell_col,
                                    label=f"{round(float(loc2weight[data[0]] / meandepth * 2), 2)}N")
        for j, data in enumerate(chr_indel):
            i = j + len(data_list)
            plot_indel(ax_array[i // cols + now_col][i % cols + len(prefix_ratios)], data, maxh,
                                    cell_col, def_cell_col, chr_len[chr_name] / maxh * 100,
                                    label=f"{round(data[3], 2)}N")

        now_col += ((len(data_list) + indel_len - 1) // cols + 1) if len(data_list) + indel_len > 0 else 0
        for col in range(bef_now_col, now_col):
            plot_scale_bar(ax_array[col][len(prefix_ratios) - 1], chr_name, maxh)



    legend_handles = []
    sorted_chr_items = sorted(CHR_COLORS.items(), key=lambda item: chr2int(item[0]))

    for chr_name, color in sorted_chr_items:
        patch = patches.Patch(facecolor=color, edgecolor='black', linewidth=0.5, label=chr_name)
        legend_handles.append(patch)

    # 범례 생성 및 배치
    gs = ax_array[0, 0].get_gridspec()
    for ax in ax_array[0]:
        ax.remove()

    label_ncol = 8
    reorder = lambda l, nc: sum((l[i::nc] for i in range(nc)), [])
    legend_ax = fig.add_subplot(gs[0, :])
    legend_ax.legend(handles=reorder(legend_handles, label_ncol),
                    labels=reorder([h.get_label() for h in legend_handles], label_ncol),
                    ncol=label_ncol,
                    loc='center',
                    fontsize=9,
                    frameon=True,       # 범례 테두리
                    handlelength=1.0,    # 핸들(색상 패치) 길이
                    handleheight=1.0,    # 핸들(색상 패치) 높이
                    labelspacing=0.8)    # 라벨 간 수직 간격

    legend_ax.set_title(f'{CELL_LINE} Virtual SKY result', fontsize=15)
    legend_ax.axis('off')

    fig.savefig(f'{PREFIX}/virtual_sky{fig_prefix}.pdf')
    fig.savefig(f'{PREFIX}/virtual_sky{fig_prefix}.png')

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
# 30_virtual_sky.py /home/hyunwoo/ACCtools-pipeline/90_skype_run/Caki-1/20_alignasm/Caki-1.ctg.aln.paf.ppc.paf /home/hyunwoo/ACCtools-pipeline/90_skype_run/Caki-1/01_depth/Caki-1_normalized.win.stat.gz public_data/chm13v2.0_telomere.bed public_data/chm13v2.0.fa.fai 30_skype_pipe/Caki-1_14_35_45 Caki-1
# """
# args = parser.parse_args(t.strip().split()[1:])

PREFIX = args.prefix
CHROMOSOME_INFO_FILE_PATH = args.reference_fai_path
main_stat_loc = args.main_stat_loc
TELOMERE_INFO_FILE_PATH = args.telomere_bed_path
PREPROCESSED_PAF_FILE_PATH = args.ppc_paf_file_path
CELL_LINE = args.cell_line_name

RATIO_OUTLIER_FOLDER = f"{PREFIX}/11_ref_ratio_outliers/"
front_contig_path = RATIO_OUTLIER_FOLDER+"front_jump/"
back_contig_path = RATIO_OUTLIER_FOLDER+"back_jump/"
TELO_CONNECT_NODES_INFO_PATH = PREFIX+"/telomere_connected_list.txt"
output_folder = f'{PREFIX}/21_pat_depth'
NCLOSE_FILE_PATH = f"{args.prefix}/nclose_nodes_index.txt"

nclose_nodes = set(extract_nclose_node(NCLOSE_FILE_PATH))
ppc_data = import_ppc_data(PREPROCESSED_PAF_FILE_PATH)

with open(f'{PREFIX}/path_data.pkl', 'rb') as f:
    path_list_dict = pkl.load(f)

df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
df = df.query('chr != "chrM"')

meandepth = np.median(df['meandepth'])
chr_len = find_chr_len(CHROMOSOME_INFO_FILE_PATH)

fclen = len(glob.glob(front_contig_path+"*"))
bclen = len(glob.glob(back_contig_path+"*"))

tot_loc_list = []

with open(f'{PREFIX}/contig_pat_vec_data.pkl', 'rb') as f:
    paf_ans_list, key_list, int2key, _ = pkl.load(f)

for loc, ll in paf_ans_list:
    tot_loc_list.append(loc)

for i in range(1, fclen//4 + 1):
    bv_paf_loc = front_contig_path+f"{i}_base.paf"
    tot_loc_list.append(bv_paf_loc)

for i in range(1, bclen//4 + 1):
    bv_paf_loc = back_contig_path+f"{i}_base.paf"
    tot_loc_list.append(bv_paf_loc)

with open(f'{PREFIX}/tot_loc_list.pkl', 'wb') as f:
    pkl.dump(tot_loc_list, f)

tot_loc_list2nclosecnt = dict()
for paf_loc in tot_loc_list:
    if paf_loc.split('/')[-3] != '11_ref_ratio_outliers':
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

build_karyotype_diagram(weights, "_raw")

build_karyotype_diagram(weights, filter_depth_N=TARGET_DEPTH)

with open(f"{PREFIX}/report.txt", 'r') as f:
    f.readline()
    path_cnt = int(f.readline().strip())

use_julia_solver = path_cnt <= HARD_PATH_COUNT_BASELINE

if use_julia_solver:
    weights_cluster = np.load(f'{PREFIX}/weight_cluster.npy')
    build_karyotype_diagram(weights_cluster, '_cluster')
