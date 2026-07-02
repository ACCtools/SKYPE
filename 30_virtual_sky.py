import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from skype_utils import *

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle as pkl
import matplotlib.patches as patches

import csv
import re
import ast
import glob
import logging
import argparse

from matplotlib.ticker import MultipleLocator
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

JOIN_BASELINE = 0.8
KARYOTYPE_SECTION_MINIMUM_LENGTH = 100 * K
KARYOTYPE_MIN_SEGMENT_LENGTH = 1 * M  # karyotype 텍스트 표기 시 이보다 짧은 segment/indel 은 무시 (1Mb)
KARYOTYPE_NORMAL_RATIO = 0.9  # 단일 염색체 path 가 reference 길이의 이 비율 미만이면 normal('1') 대신 del 로 취급

CTG_INTYPE_CHECK_MIN_LENGTH = 100 * K
CTG_INTYPE_INSERT_MIN_SEGMENT_LENGTH = 10 * K
CTG_INTYPE_INSERT_MIN_RATIO = 0.2



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

def virtual_event_label(event_type : str) -> str:
    return {'d': 'del', 'i': 'ins', 'v': 'inv'}[event_type]

# SKY Figure
def plot_virtual_chromosome(ax, segments_data, maxh, cell_col, def_cell_col,
                            label=None, karyotype_str=None, event_labels=None):
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
    event_labels = event_labels or []
    event_label_y_positions = [event_y for event_y, _event_label in event_labels]

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

        is_annotated_event_boundary = any(
            abs(current_y - event_y) < 1e-6 for event_y in event_label_y_positions
        )
        if 0 < i and i < len(segments) and not is_annotated_event_boundary:
            last_chr = segments[i-1][0]
            if real_chr[0] != last_chr[0]:
                text_label = f"t({last_chr[0][3:]};{real_chr[0][3:]})"
            else:
                text_label = f"t({last_chr[0][3:]+str(last_chr[1])};{real_chr[0][3:] + str(real_chr[1])})"
            text_obj = ax.text(text_x, current_y, text_label, **decoration_args)
            text_obj_list.append(text_obj)

        current_y += seg_length
        
    for event_y, event_label in event_labels:
        event_y = max(0, min(100, event_y))
        text_obj = ax.text(text_x, event_y, event_label, **decoration_args)
        text_obj_list.append(text_obj)

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
    # 경로(path) 대신 karyotype ISCN 표기를 표시
    if karyotype_str:
        ax.text(x_center, -10, karyotype_str, ha='center', va='top', fontsize=10)

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
    background_color = 'white' if indel[0] in {'i', 'v'} else chr_color
    face_color = chr_color if indel[0] in {'i', 'v'} else 'white'

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
    
    event_label = virtual_event_label(indel[0])
    text_obj = ax.text(x_end+radius/5, midy, f"{event_label}({indel[-2][3:]})", ha='left', va='center', fontsize=10, color='black')

    mark_overlapping_texts_with_arrows(ax, [text_obj], min_gap=5)

    if label:
        ax.text(x_center, -5, label, ha='center', va='top', fontsize=10)
    # 아래: 다른 염색체와 동일하게 ISCN indel 표기
    ax.text(x_center, -10, f"{event_label}({chrom_to_iscn(indel[4])})", ha='center', va='top', fontsize=10)

    ax.set_xlim(0, 60 * cell_col / def_cell_col)
    ax.set_ylim(0, 100)
    ax.axis('off')

def plot_centromere_fragment(ax, frag, maxh, cell_col, def_cell_col, chr_len_norm, label=None):
    """frag = (chrom, side, mid_bp, chr_len_bp).
    점선 outline 으로 reference 염색체 윤곽을 그리고, centromere mid 를 기준으로
    한쪽 arm fragment 만 색칠한다.
    """
    chrom, side, mid_bp, chr_len_bp = frag
    x_center = 0
    width = 10
    radius = width / 2.0
    total_length = chr_len_norm
    round_radius = min(radius, total_length / 4) if total_length > 0 else radius
    chr_color = CHR_COLORS.get(chrom, "gray")

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
        facecolor='white',
        edgecolor='black',
        linestyle='dashed',
        linewidth=1.2,
        clip_on=False
    )
    ax.add_patch(outline_patch)

    mid_norm = mid_bp / chr_len_bp * total_length
    if side == 'right':
        y1, y2 = mid_norm, total_length
    else:
        y1, y2 = 0, mid_norm
    seg_length = y2 - y1

    frag_rect = patches.Rectangle(
        (x_center - radius, y1),
        width,
        seg_length,
        facecolor=chr_color,
        edgecolor='none'
    )
    frag_rect.set_clip_path(clip_patch)
    ax.add_patch(frag_rect)

    midy = (y1 + y2) / 2
    x_start = x_center + radius
    x_end = x_center + 2 * radius
    ax.plot([x_start, x_end], [y1, midy], linewidth=1, color='black')
    ax.plot([x_start, x_end], [y2, midy], linewidth=1, color='black')

    arm = 'q' if side == 'right' else 'p'  # right=q arm, left=p arm
    text_obj = ax.text(
        x_end + radius / 5, midy,
        f"{chrom}{arm}",   # 옆 라벨: chr13q / chr13p (한쪽 arm 이름)
        ha='left', va='center', fontsize=10, color='black'
    )
    mark_overlapping_texts_with_arrows(ax, [text_obj], min_gap=5)

    if label:
        ax.text(x_center, -5, label, ha='center', va='top', fontsize=10)
    # 아래: 다른 염색체와 동일하게 ISCN isochromosome 표기
    ax.text(x_center, -10, f"i({chrom_to_iscn(chrom)})({arm}10)", ha='center', va='top', fontsize=10)

    ax.set_xlim(0, 60 * cell_col / def_cell_col)
    ax.set_ylim(0, 100)
    ax.axis('off')

def plot_ecdna(ax, chr_name, cell_col, def_cell_col, label = None):
    axis_xlim = 60
    axis_ylim = 100
    x_center = 30 * cell_col / def_cell_col  # 좌우 중앙값 (필요에 따라 조정)
    width = 10     # 염색체 
    radius = width / 1.3
    chr_color = CHR_COLORS.get(chr_name[0])

    ax.set_xlim(0, axis_xlim * cell_col / def_cell_col)
    ax.set_ylim(0, axis_ylim)
    ax.axis('off')

    outline_circle = patches.Ellipse(
        (x_center, 50),
        radius*2*1.2, # Todo : 이 상수가 뭔지 이해하기 (일단 ylim/xlim은 아님)
        radius*2,
        facecolor=chr_color,
        edgecolor='black',
        linewidth=1.2
    )
    ax.add_patch(outline_circle)
    hole_circle = patches.Ellipse(
        (x_center, 50),
        radius*1.2,
        radius,
        facecolor='white',
        edgecolor='black',
        linewidth=1.2
    )
    ax.add_patch(hole_circle)
    if label:
        ax.text(x_center, -5, label, ha='center', va='top', fontsize=10)

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

def get_ctg_intype_interrupt_pieces(key_int : int, endpoint_chroms : set) -> list:
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

    pieces = []
    for i, (row, length) in enumerate(zip(rows, row_lengths)):
        chrom = row[CHR_NAM]
        if chrom not in interrupt_chroms:
            continue
        if length < CTG_INTYPE_INSERT_MIN_SEGMENT_LENGTH:
            continue
        strand = infer_path_strand_from_paf_rows(rows, i)
        append_ctg_intype_interrupt_piece(pieces, chrom, strand, length)

    return pieces

def get_karyotype_summary_from_index(path_path : str, type4_edge_to_event_key=None,
                                     type4_event_by_key=None) -> list:
    """
    Fallback summary from the compact index path. This keeps the old behavior
    for prefixes that do not have the 21_pat_depth PAF fragments available.
    """
    pieces = []
    path = import_index_path(path_path)
    ctg_intype_key_by_edge = {}
    for key_int in path2key_int_list.get(path_path, []):
        key_type, key_value = int2key[key_int]
        if key_type == CTG_IN_TYPE:
            ctg_intype_key_by_edge[key_value] = key_int
    
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
                          type4_event_by_key=None):
    """
    Summarizes karyotype data from the compact index path. Only CTG_IN_TYPE
    edges are inspected in the expanded PAF fragments to reveal long inserted
    sequence from chromosomes outside the edge endpoints.
    """
    karyotypes_data_direction_include = {}
    
    for path_path in non_type4_path_list:
        pieces = get_karyotype_summary_from_index(
            path_path, type4_edge_to_event_key, type4_event_by_key
        )
        karyotypes_data_direction_include[path_path] = pieces
        
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

    karyotypes_data = get_karyotype_summary(
        non_type4_top_path, type4_edge_to_event_key, type4_event_by_key
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

    cols = 10
    display_chroms = sorted(
        set(grouped_norm_data.keys()) | set(display_indel.keys()),
        key=chr2int,
    )

    rows = 0
    for chr_name in display_chroms:
        data_list = grouped_norm_data.get(chr_name, [])
        chr_indel = display_indel.get(chr_name, [])
        indel_len = len(chr_indel)

        rows += ((len(data_list) + indel_len - 1) // cols + 1) if len(data_list) + indel_len > 0 else 0

    if len(fragment_display) > 0:
        rows += (len(fragment_display) - 1) // cols + 1

    if len(virtual_inv_display) > 0:
        rows += (len(virtual_inv_display) - 1) // cols + 1

    # if len(long_ecdna_path) > 0:
    #     rows += ((len(long_ecdna_path)-1) // cols + 1)

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

    sorted_grouped_norm_data_items = [
        (chr_name, grouped_norm_data.get(chr_name, []))
        for chr_name in display_chroms
    ]

    now_col = 1
    for chr_name, data_list in sorted_grouped_norm_data_items:
        chr_indel = display_indel.get(chr_name, [])
        indel_len = len(chr_indel)

        bef_now_col = now_col
        plot_chr_name(ax_array[now_col][0], chr_name)
        for i, data in enumerate(data_list):
            path_type4_events = get_path_type4_indel_events(
                data[0], type4_event_by_key, type4_path_event_usage
            )
            iscn = karyotype_path_to_iscn(karyotypes_data[data[0]], path_type4_events)
            event_labels = get_type4_indel_boundary_labels(
                data[0], type4_edge_to_event_key, type4_event_by_key, maxh
            )
            plot_virtual_chromosome(ax_array[i // cols + now_col][i % cols + len(prefix_ratios)], data, maxh,
                                    cell_col, def_cell_col,
                                    label=f"{round(float(loc2weight[data[0]] / meandepth * 2), 2)}N",
                                    karyotype_str=iscn if iscn is not None else '',
                                    event_labels=event_labels)
        for j, data in enumerate(chr_indel):
            i = j + len(data_list)
            plot_indel(ax_array[i // cols + now_col][i % cols + len(prefix_ratios)], data, maxh,
                                    cell_col, def_cell_col, chr_len[chr_name] / maxh * 100,
                                    label=f"{round(data[3], 2)}N")

        now_col += ((len(data_list) + indel_len - 1) // cols + 1) if len(data_list) + indel_len > 0 else 0
        for col in range(bef_now_col, now_col):
            plot_scale_bar(ax_array[col][len(prefix_ratios) - 1], chr_name, maxh)

    if len(virtual_inv_display) > 0:
        bef_now_col = now_col
        plot_chr_name(ax_array[now_col][0], 'virtual inv')
        for i, data in enumerate(virtual_inv_display):
            chrom = data[4]
            plot_indel(
                ax_array[i // cols + now_col][i % cols + len(prefix_ratios)],
                data,
                maxh,
                cell_col,
                def_cell_col,
                chr_len[chrom] / maxh * 100,
                label=f"{round(data[3], 2)}N",
            )
        now_col += (len(virtual_inv_display) - 1) // cols + 1
        for col in range(bef_now_col, now_col):
            plot_scale_bar(ax_array[col][len(prefix_ratios) - 1], 'virtual inv', maxh)

    if len(fragment_display) > 0:
        bef_now_col = now_col
        plot_chr_name(ax_array[now_col][0], 'cen frag')
        for i, frag in enumerate(fragment_display):
            chrom, side, mid_bp, chr_len_bp, depth_N = frag
            chr_len_norm = chr_len_bp / maxh * 100
            plot_centromere_fragment(
                ax_array[i // cols + now_col][i % cols + len(prefix_ratios)],
                (chrom, side, mid_bp, chr_len_bp),
                maxh, cell_col, def_cell_col, chr_len_norm,
                label=f"{round(depth_N, 2)}N"
            )
        now_col += (len(fragment_display) - 1) // cols + 1
        for col in range(bef_now_col, now_col):
            plot_scale_bar(ax_array[col][len(prefix_ratios) - 1], 'cen', maxh)

    # if len(long_ecdna_path) > 0:
    #     plot_chr_name(ax_array[now_col][0], 'ecDNA')
    #     for i, data in enumerate(long_ecdna_path.values()):
    #         plot_ecdna(ax_array[i // cols + now_col][i % cols + len(prefix_ratios)], data[0], cell_col, def_cell_col, label=f"{data[0]} : {ecdna_format(data[1])}-{ecdna_format(data[2])}\n({round(data[3], 2)}N)")



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

    # ----- ISCN karyotype 텍스트 출력: 그림에 그린 breakend 를 표기법으로 정리 -----
    # 그림과 동일한 순서(염색체 그룹별 path -> 같은 염색체 indel -> centromere fragment).
    karyotype_rows = []  # (iscn_str, depth_N)
    for chr_name, data_list in sorted_grouped_norm_data_items:
        for path, _norm in data_list:
            path_type4_events = get_path_type4_indel_events(
                path, type4_event_by_key, type4_path_event_usage
            )
            iscn = karyotype_path_to_iscn(karyotypes_data[path], path_type4_events)
            if iscn is not None:
                karyotype_rows.append((iscn, loc2weight[path] / meandepth * 2))
        for indel in display_indel.get(chr_name, []):
            typ, pos1, pos2, depth_N, chrom, _path = indel
            if abs(pos2 - pos1) < KARYOTYPE_MIN_SEGMENT_LENGTH:
                continue
            name = chrom_to_iscn(chrom)
            karyotype_rows.append((f"{virtual_event_label(typ)}({name})", depth_N))
    for inv_event in virtual_inv_display:
        typ, pos1, pos2, depth_N, chrom, _name = inv_event
        if abs(pos2 - pos1) < KARYOTYPE_MIN_SEGMENT_LENGTH:
            continue
        karyotype_rows.append((f"{virtual_event_label(typ)}({chrom_to_iscn(chrom)})", depth_N))
    for chrom, side, mid_bp, chr_len_bp, depth_N in fragment_display:
        arm = 'q' if side == 'right' else 'p'  # right=q arm, left=p arm
        karyotype_rows.append((f"i({chrom_to_iscn(chrom)})({arm}10)", depth_N))

    with open(f'{PREFIX}/karyotype{fig_prefix}.txt', 'w') as kf:
        kf.write("karyotype\tdepth\n")
        for iscn, depth_N in karyotype_rows:
            kf.write(f"{iscn}\t{round(float(depth_N), 2)}N\n")

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
logging.info(describe_pipeline_mode(pipeline_mode_config))

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

df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
df = df.query('chr != "chrM"')

meandepth = np.median(df['meandepth'])
chr_len = find_chr_len(CHROMOSOME_INFO_FILE_PATH)

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
