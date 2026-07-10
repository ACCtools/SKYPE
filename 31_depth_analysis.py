import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from skype_utils import *

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import pickle as pkl

import csv
import ast
import logging
import argparse
import collections
import scipy.stats
import glob
import re
import vcfpy


from scipy.signal import butter, filtfilt
from pycirclize import Circos
from matplotlib.lines import Line2D
from matplotlib.projections.polar import PolarAxes
from pycirclize.track import Track
from collections import defaultdict
from collections import Counter
from collections import OrderedDict

# logging 설정(레벨/포맷)은 skype_utils 에서 중앙 관리한다 (LOG_LEVEL).
logging.info("31_depth_analysis start")

BREAKEND_REMARKABLE_CN_RATIO = 0.05
TELOMERE_REMARKABLE_CN_RATIO = 0.05
VCF_FILTER_DEPTH_N = 0.1
RAW_TRANSLOCATION_RESULT_PKL = 'raw_translocation_result.pkl'
RAW_TRANSLOCATION_REPORT_TSV = 'raw_translocation_read_counts.tsv'

BND_TYPE = 0
CTG_IN_TYPE = 1
TEL_TYPE = 2

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

NODE_NAME = 1
CHR_CHANGE_IDX = 2
DIR_CHANGE_IDX = 3

ABS_MAX_COVERAGE_RATIO = 3
MAX_PATH_CNT = 100

DIR_FOR = 1
DIR_BAK = 0
TELOMERE_EXPANSION = 5 * K

CONJOINED_CONTIG_MINIMUM_LENGTH = 200*K
SIM_COMPARE_RAITO = 1.2

TYPE2_FLANKING_LENGTH = 5*M
TYPE2_SIM_COMPARE_RAITO = 1.5
TYPE34_BREAK_CHUKJI_LIMIT = 1*M

NCLOSE_SIM_COMPARE_RAITO = 1.2

CTG_INTYPE_CHECK_MIN_LENGTH = 100 * K
CTG_INTYPE_INSERT_MIN_SEGMENT_LENGTH = 10 * K
CTG_INTYPE_INSERT_MIN_RATIO = 0.2


def similar_check(v1, v2, ratio=TYPE2_SIM_COMPARE_RAITO):
    try:
        assert(v1 >= 0 and v2 >= 0)
    except:
        logging.error(f"Invalid values for similarity check: v1={v1}, v2={v2}")
        assert(False)
    mi, ma = sorted([v1, v2])
    return False if mi == 0 else (ma / mi <= ratio) or ma - mi < NCLOSE_SIM_DIFF_THRESHOLD

def exist_near_bnd(chrom, inside_st, inside_nd):
    # subset of df for the given chromosome
    df_chr = df[df['chr'] == chrom]

    def mean_depth(start, end):
        """Return mean meandepth over windows overlapping [start, end)."""
        mask = (df_chr['nd'] > start) & (df_chr['st'] < end)
        return df_chr.loc[mask, 'meandepth'].mean()

    # for inside_st
    st_depth = mean_depth(inside_st - TYPE2_FLANKING_LENGTH, inside_st)
    nd_depth = mean_depth(inside_nd, inside_nd + TYPE2_FLANKING_LENGTH)
    
    if np.isnan(st_depth) or np.isnan(nd_depth):
        return True

    # print(chrom, inside_st, inside_nd, not similar_check(st_depth, nd_depth))
    return not similar_check(st_depth, nd_depth, NCLOSE_SIM_COMPARE_RAITO)


def import_index_path(file_path : str) -> list:
    file_path_list = file_path.split('/')
    key = file_path_list[-2]
    cnt = int(file_path_list[-1].split('.')[0]) - 1

    return path_list_dict[key][cnt][0]

def extract_groups(lst):
    if not lst:
        return []
    
    result = []
    seen = set()
    current = lst[0]
    result.append(current)
    seen.add(current)
    
    for num in lst[1:]:
        if num != current:
            # 새로운 숫자가 등장했는데 이미 이전에 등장한 적이 있다면 에러 처리
            if num in seen:
                raise ValueError(f"Error: {num}가 연속된 구간 이후에 다시 등장합니다.")
            result.append(num)
            seen.add(num)
            current = num
    return result

def line_track_plot(
    line_track: Track,
    x: list[float] | np.ndarray,
    y: list[float] | np.ndarray,
    *,
    vmin: float = 0,
    vmax: float | None = None,
    arc: bool = True,
    **kwargs,) -> None:
    # Check x, y list length
    if len(x) != len(y):
        err_msg = f"List length is not match ({len(x)=}, {len(y)=})"
        raise ValueError(err_msg)

    # Convert (x, y) to (rad, r)
    rad = list(map(line_track.x_to_rad, x))
    vmax = max(y) if vmax is None else vmax
    # line_track._check_value_min_max(y, vmin, vmax)
    r = [line_track._y_to_r(v, vmin, vmax) for v in y]

    if arc:
        # Convert linear line to arc line (rad, r) points
        plot_rad, plot_r = line_track._to_arc_radr(rad, r)
    else:
        plot_rad, plot_r = rad, r

    def plot_line(ax: PolarAxes) -> None:
        line, = ax.plot(plot_rad, plot_r, **kwargs)
        clip_patch = plt.Circle((0, 0), max(line_track.r_plot_lim), transform=ax.transProjectionAffine + ax.transAxes)

        clip_patch.set_visible(False)
        ax.add_patch(clip_patch)
        line.set_clip_path(clip_patch)

    line_track._plot_funcs.append(plot_line)

def line_track_circos(
    line_track: Track,
    x: list[float] | np.ndarray,
    y1: list[float] | np.ndarray,
    *,
    vmin: float = 0,
    vmax: float | None = None,
    arc: bool = True,
    **kwargs,
) -> None:
    y2 = 0
    rad = list(map(line_track.x_to_rad, x))
    if isinstance(y2, (list, tuple, np.ndarray)):
        y_all = list(y1) + list(y2)
    else:
        y_all = list(y1) + [y2]
        y2 = [float(y2)] * len(x)
    vmin = min(y_all) if vmin is None else vmin
    vmax = max(y_all) if vmax is None else vmax
    # line_track._check_value_min_max(y_all, vmin, vmax)

    r2 = [line_track._y_to_r(v, vmin, vmax) for v in y2]
    r = [line_track._y_to_r(v, vmin, vmax) for v in y1]
    if arc:
        plot_rad, plot_r2 = line_track._to_arc_radr(rad, r2)
        _, plot_r = line_track._to_arc_radr(rad, r)
    else:
        plot_rad, plot_r, plot_r2 = rad, r, r2

    def plot_fill_between(ax: PolarAxes) -> None:
        line = ax.fill_between(plot_rad, plot_r, plot_r2, **kwargs)  # type: ignore

        clip_patch = plt.Circle((0, 0), max(line_track.r_plot_lim), transform=ax.transProjectionAffine + ax.transAxes)

        clip_patch.set_visible(False)
        ax.add_patch(clip_patch)
        line.set_clip_path(clip_patch)

    line_track._plot_funcs.append(plot_fill_between)

def line_track_circos_color(
    line_track: Track,
    x: list[float] | np.ndarray,
    y1: list[float] | np.ndarray,
    *,
    vmin: float = 0,
    vmax: float | None = None,
    target_value: list[float],
    target_color: list[str],
    **kwargs,
) -> None:
    y2 = 0
    arc = False

    # Verify that target_value is sorted in ascending order
    if any(np.diff(target_value) < 0):
        raise ValueError("target_value must be sorted in ascending order")

    # Verify that the number of target_color is exactly one more than the number of target_value boundaries
    if len(target_color) != len(target_value) + 1:
        raise ValueError("The number of target_color must be equal to len(target_value) + 1")

    rad = list(map(line_track.x_to_rad, x))
    if isinstance(y2, (list, tuple, np.ndarray)):
        y_all = list(y1) + list(y2)
    else:
        y_all = list(y1) + [y2]
        y2 = [float(y2)] * len(x)
    vmin = min(y_all) if vmin is None else vmin
    vmax = max(y_all) if vmax is None else vmax
    # line_track._check_value_min_max(y_all, vmin, vmax)

    r2 = [line_track._y_to_r(v, vmin, vmax) for v in y2]
    r = [line_track._y_to_r(v, vmin, vmax) for v in y1]
    if arc:
        plot_rad, plot_r2 = line_track._to_arc_radr(rad, r2)
        _, plot_r = line_track._to_arc_radr(rad, r)
    else:
        plot_rad, plot_r, plot_r2 = rad, r, r2

    plot_r = np.asarray(plot_r)
    plot_r2 = np.asarray(plot_r2)
    def plot_fill_between(ax: PolarAxes) -> None:
        # Convert y1 to a numpy array for element-wise operations
        y1_arr = np.array(y1)

        # Define segments based on target_value boundaries:
        # Segment 0: y1 < target_value[0]
        # Segment i: target_value[i-1] <= y1 < target_value[i] for i = 1,...,len(target_value)-1
        # Segment last: y1 >= target_value[-1]
        segments = []
        segments.append((None, target_value[0]))  # lower segment (no lower bound)
        for i in range(1, len(target_value)):
            segments.append((target_value[i-1], target_value[i]))
        segments.append((target_value[-1], None))  # upper segment (no upper bound)

        for (seg_lower, seg_upper), color in zip(segments, target_color):
            if seg_lower is None:
                mask = y1_arr < seg_upper
            elif seg_upper is None:
                mask = y1_arr >= seg_lower
            else:
                mask = (y1_arr >= seg_lower) & (y1_arr < seg_upper)
            
            # mask 경계에 한 칸씩 확장하여 gap을 없앰 (dilation)
            mask[1:] = mask[1:] | mask[:-1]
            # 채우기 영역 그리기
            line = ax.fill_between(plot_rad, plot_r, plot_r2, where=mask, color=color, **kwargs)  # type: ignore

            # 채우기 영역 제한을 위한 보이지 않는 클리핑 패치 생성
            clip_patch = plt.Circle(
                (0, 0),
                max(line_track.r_plot_lim),
                transform=ax.transProjectionAffine + ax.transAxes
            )
            clip_patch.set_visible(False)
            ax.add_patch(clip_patch)
            line.set_clip_path(clip_patch)
    line_track._plot_funcs.append(plot_fill_between)

def find_chr_len(file_path : str) -> dict:
    chr_data_file = open(file_path, "r")
    chr_len = {}
    for curr_data in chr_data_file:
        curr_data = curr_data.split("\t")
        chr_len[curr_data[0]] = int(curr_data[1])
    chr_data_file.close()
    return chr_len

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
    
def chr2int(x):
    chrXY2int = {'chrX' : 24, 'chrY' : 25}
    if x in chrXY2int:
        return chrXY2int[x]
    else:
        return int(x[3:])

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

def highpass_filter(data, cutoff, fs, order=3):
    """
    Butterworth high-pass filter를 이용하여 저주파 성분을 제거합니다.
    
    Parameters:
    - data: 1D array, 입력 신호 (CN 값)
    - cutoff: float, 컷오프 주파수
    - fs: float, 샘플링 주파수 (데이터 포인트 간 간격)
    - order: int, 필터 차수
    
    Returns:
    - 필터링된 신호 (노이즈 성분 추출)
    """
    nyq = 0.5 * fs  # Nyquist 주파수
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='high', analog=False)
    return filtfilt(b, a, data)

def rebin_dataframe(df: pd.DataFrame, n: int) -> pd.DataFrame:
    """
    Group rows in the DataFrame into bins spanning n consecutive units.
    The unit is determined from the 'length' of the first row in each chromosome group.
    
    For each group:
      - The new start (st) is taken from the first row.
      - The new end (nd) is taken from the last row in the group.
      - For a complete group (n rows), values are summed.
      - For an incomplete group (fewer than n rows), the new bin's length is the sum of the available lengths,
        and new covsite and totaldepth are computed as the weighted average:
            new_value = sum( value_i * length_i ) / sum(length_i)
      - New coverage (cov) is computed as (new_covsite / new_length) * 100.
      - New mean depth (meandepth) is computed as new_totaldepth / new_length.
    
    Parameters:
        df (pd.DataFrame): Input DataFrame with columns ['chr', 'st', 'nd', 'length',
                          'covsite', 'totaldepth', 'cov', 'meandepth'].
        n (int): Number of consecutive units (rows) to combine into each bin.
    
    Returns:
        pd.DataFrame: New DataFrame with binned rows.
    """
    new_rows = []
    
    # Process each chromosome separately.
    for chrom, sub_df in df.groupby('chr'):
        # Sort rows by starting position.
        sub_df = sub_df.sort_values('st').reset_index(drop=True)
        
        # Group rows in chunks of size n.
        for i in range(0, len(sub_df), n):
            chunk = sub_df.iloc[i:i+n]
            new_st = chunk['st'].iloc[0]
            new_nd = chunk['nd'].iloc[-1]
            
            # Use the actual sum of lengths in the chunk.
            sum_length = chunk['length'].sum()
            new_meandepth = np.sum(chunk['totaldepth']) / sum_length
            
            new_rows.append({
                'chr': chrom,
                'st': new_st,
                'nd': new_nd,
                'meandepth': new_meandepth
            })
    
    return pd.DataFrame(new_rows)

def rebin_dataframe_B(df: pd.DataFrame, n: int) -> np.array:
    new_rows = dict()
    
    # Process each chromosome separately.
    for chrom, sub_df in df.groupby('chr'):
        chr_mean_list = []
        sub_df = sub_df.sort_values('st').reset_index(drop=True)
        
        # Group rows in chunks of size n.
        for i in range(0, len(sub_df), n):
            chunk = sub_df.iloc[i:i+n]
            new_st = chunk['st'].iloc[0]
            new_nd = chunk['nd'].iloc[-1]
            
            # Use the actual sum of lengths in the chunk.
            sum_length = chunk['length'].sum()
            new_meandepth = np.sum(chunk['totaldepth']) / sum_length
            
            chr_mean_list.extend([new_meandepth for _ in range(len(chunk))])
        
        new_rows[chrom] = np.asarray(chr_mean_list)
    
    ans_B = []
    for c in chr_order_list:
        ans_B.append(new_rows[c])

    return np.hstack(ans_B)

def import_data2(file_path : str) -> list :
    paf_file = open(file_path, "r")
    contig_data = []
    for curr_contig in paf_file:
        temp_list = curr_contig.split("\t")
        int_induce_idx = [
            CTG_LEN, CTG_STR, CTG_END, \
                          CHR_LEN, CHR_STR, CHR_END, \
                          CTG_MAPQ, CTG_TYP, CTG_STRND, CTG_ENDND,]
        for i in int_induce_idx:
            temp_list[i] = int(temp_list[i])
        contig_data.append(tuple(temp_list))
    paf_file.close()
    return contig_data

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

def import_bed(bed_path : str) -> dict:
    bed_data_file = open(bed_path, "r")
    chr_len = collections.defaultdict(list)
    for curr_data in bed_data_file:
        curr_data = curr_data.split("\t")
        chr_len[curr_data[0]].append((int(curr_data[1]), int(curr_data[2])))
    bed_data_file.close()
    return chr_len

def inclusive_checker(tuple_a : tuple, tuple_b : tuple) -> bool :
    if int(tuple_a[0]) <= int(tuple_b[0]) and int(tuple_b[1]) <= int(tuple_a[1]):
        return True
    else:
        return False

def ordered_mapping(*items):
    return OrderedDict(items)


def add_alt_header_line(header, alt_id, description):
    header.add_line(vcfpy.AltAlleleHeaderLine.from_mapping(ordered_mapping(
        ("ID", alt_id),
        ("Description", description),
    )))


def add_info_header_line(header, info_id, number, type_, description):
    header.add_info_line(ordered_mapping(
        ("ID", info_id),
        ("Number", number),
        ("Type", type_),
        ("Description", description),
    ))


def build_skype_vcf_header(contig_lengths):
    header = vcfpy.Header(lines=[vcfpy.HeaderLine("fileformat", "VCFv4.3")])
    add_alt_header_line(header, "BND", "Breakend")
    add_alt_header_line(header, "INV", "Inversion")
    add_alt_header_line(header, "DEL", "Deletion")
    add_alt_header_line(header, "DUP", "Duplication")
    add_info_header_line(header, "SVTYPE", 1, "String", "Type of structural variant")
    add_info_header_line(header, "END", 1, "Integer", "End position of SV")
    add_info_header_line(header, "SVLEN", 1, "Integer", "Length of the SV")
    add_info_header_line(header, "WEIGHT", 1, "Float", "Depth for breakend")
    add_info_header_line(header, "CTG_NAME", 1, "String", "Name of contig for supporting variant")
    add_info_header_line(header, "SVCLASS", 1, "String", "SKYPE event class")
    add_info_header_line(header, "STRANDS", 1, "String", "Breakpoint strandedness")
    add_info_header_line(header, "MATEID", ".", "String", "ID of mate breakends")
    add_info_header_line(header, "MERGE_MATEID", ".", "String", "IDs of merged breakends")
    for chrom, length in contig_lengths.items():
        header.add_contig_line(ordered_mapping(
            ("ID", chrom),
            ("length", int(length)),
        ))
    return header


def vcf_filter_list(filter_str):
    if filter_str in (None, "", "."):
        return []
    if isinstance(filter_str, str):
        return filter_str.split(";")
    return list(filter_str)


def make_vcf_record(chrom, pos, sv_id, ref, alt, quality, filter_str, info):
    return vcfpy.Record(
        CHROM=str(chrom),
        POS=int(pos),
        ID=[] if sv_id in (None, "", ".") else [str(sv_id)],
        REF=str(ref),
        ALT=[alt],
        QUAL=quality,
        FILTER=vcf_filter_list(filter_str),
        INFO=info,
    )


BND_ALT_FORM_TO_ORIENTATION = {
    "t[p[": ("-", "+"),
    "t]p]": ("-", "-"),
    "]p]t": ("+", "-"),
    "[p[t": ("+", "+"),
}


def breakend_alt(mate_chr, mate_pos, form, sequence="N"):
    orientation, mate_orientation = BND_ALT_FORM_TO_ORIENTATION[form]
    return vcfpy.BreakEnd(
        str(mate_chr),
        int(mate_pos),
        orientation,
        mate_orientation,
        str(sequence),
        True,
    )

def choose_alt_forms(dir_a, dir_b):
    """
    Map SKYPE path directions (dir_a, dir_b) to a VCF 4.3 breakend ALT pair.
    Convention: a = arm1 (junction at its contig-3' end), b = arm2 (junction at its contig-5' end).
    - dir_a decides whether a's retained sequence is left ('+': t-prefix) or right ('-': t-suffix).
    - dir_b decides whether the joined arm2 piece extends right of p forward ('+': '[') or
      is the reverse-complement extending left of p ('-': ']').
    Verified against VCFv4.3 §5.4 (Fig.1 all-orientations, Fig.7 RR0 transloc, Fig.8 INV0 inversion);
    valid reciprocal mate pairs are  t[p[ <-> ]p]t ,  t]p] <-> t]p] ,  [p[t <-> [p[t .
    """
    if dir_a == '+' and dir_b == '+':
        return ("t[p[", "]p]t")
    if dir_a == '+' and dir_b == '-':
        return ("t]p]", "t]p]")
    if dir_a == '-' and dir_b == '+':
        return ("[p[t", "[p[t")
    if dir_a == '-' and dir_b == '-':
        return ("]p]t", "t[p[")
    assert(False)

def make_strands(dir_a, dir_b):
    """
    Compose VCF breakend STRANDS from SKYPE path directions.

    SKYPE nclose notation records traversal as A_dir => B_dir.  VCF STRANDS
    records the two breakpoint sides that are newly adjacent.  The A exit side
    is A_dir, while the B entry side is the opposite of B_dir.
    """
    a = dir_a if dir_a in ('+', '-') else '.'
    b = invert_strand(dir_b) if dir_b in ('+', '-') else '.'
    return f"{a}{b}"

def write_bnd_vcf_pair(
    writer,
    sv_id_base,
    chr_a,
    pos_a,
    dir_a,
    chr_b,
    pos_b,
    dir_b,
    weight_N,
    ctg_name,
    quality=60,
    filter_str='.',
    merge_mate_ids=None,
):
    pos_a = max(1, int(pos_a))
    pos_b = max(1, int(pos_b))
    strands = make_strands(dir_a, dir_b)
    form_a, form_b = choose_alt_forms(dir_a, dir_b)
    sv_id_a = f"{sv_id_base}_1"
    sv_id_b = f"{sv_id_base}_2"
    ref = "N"
    merge_mate_id_value = None
    if merge_mate_ids:
        merge_mate_id_value = [str(mate_id) for mate_id in merge_mate_ids]

    info_a = ordered_mapping(
        ("SVTYPE", "BND"),
        ("WEIGHT", round(weight_N, 2)),
        ("CTG_NAME", str(ctg_name)),
        ("STRANDS", strands),
        ("MATEID", [sv_id_b]),
    )
    info_b = ordered_mapping(
        ("SVTYPE", "BND"),
        ("WEIGHT", round(weight_N, 2)),
        ("CTG_NAME", str(ctg_name)),
        ("STRANDS", strands),
        ("MATEID", [sv_id_a]),
    )
    if merge_mate_id_value is not None:
        info_a["MERGE_MATEID"] = merge_mate_id_value
        info_b["MERGE_MATEID"] = merge_mate_id_value

    writer.write_record(make_vcf_record(
        chr_a, pos_a, sv_id_a, ref, breakend_alt(chr_b, pos_b, form_a, ref),
        quality, filter_str, info_a
    ))
    writer.write_record(make_vcf_record(
        chr_b, pos_b, sv_id_b, ref, breakend_alt(chr_a, pos_a, form_b, ref),
        quality, filter_str, info_b
    ))

def invert_strand(strand):
    return '-' if strand == '+' else '+'

def make_breakend_endpoint(chrom, pos, strand):
    return {
        'chrom': chrom,
        'pos': int(pos),
        'strand': strand,
    }

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

def paf_row_entry_endpoint(row, strand):
    pos = row[CHR_STR] if strand == '+' else row[CHR_END]
    return make_breakend_endpoint(row[CHR_NAM], pos, strand)

def paf_row_exit_endpoint(row, strand):
    pos = row[CHR_END] if strand == '+' else row[CHR_STR]
    return make_breakend_endpoint(row[CHR_NAM], pos, strand)

def append_ctg_intype_interrupt_segment(
    segments : list,
    row,
    strand : str,
    length : int,
    row_idx : int,
):
    if length <= 0:
        return
    chrom = row[CHR_NAM]
    if segments and segments[-1]['chrom'] == chrom:
        segments[-1]['exit'] = paf_row_exit_endpoint(row, strand)
        segments[-1]['length'] += length
        segments[-1]['last_idx'] = row_idx
    else:
        segments.append({
            'chrom': chrom,
            'strand': strand,
            'entry': paf_row_entry_endpoint(row, strand),
            'exit': paf_row_exit_endpoint(row, strand),
            'length': length,
            'first_idx': row_idx,
            'last_idx': row_idx,
        })

def nearest_different_chrom_row(rows : list, row_idx : int, step : int, chrom : str):
    idx = row_idx + step
    while 0 <= idx < len(rows):
        if rows[idx][CHR_NAM] != chrom:
            return idx
        idx += step
    return None

def attach_ctg_intype_interrupt_flanks(segments : list, rows : list):
    for segment in segments:
        left_idx = nearest_different_chrom_row(
            rows, segment['first_idx'], -1, segment['chrom']
        )
        if left_idx is not None:
            left_strand = infer_path_strand_from_paf_rows(rows, left_idx)
            segment['left_flank_exit'] = paf_row_exit_endpoint(rows[left_idx], left_strand)

        right_idx = nearest_different_chrom_row(
            rows, segment['last_idx'], 1, segment['chrom']
        )
        if right_idx is not None:
            right_strand = infer_path_strand_from_paf_rows(rows, right_idx)
            segment['right_flank_entry'] = paf_row_entry_endpoint(rows[right_idx], right_strand)

def get_ctg_intype_interrupt_segments(key_int : int, endpoint_chroms : set) -> list:
    """
    Coordinate-bearing mirror of 30_virtual_sky.get_ctg_intype_interrupt_pieces().
    Keep the same chromosome selection/filtering; only attach entry/exit coords.
    """
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

    segments = []
    for i, (row, length) in enumerate(zip(rows, row_lengths)):
        chrom = row[CHR_NAM]
        if chrom not in interrupt_chroms:
            continue
        if length < CTG_INTYPE_INSERT_MIN_SEGMENT_LENGTH:
            continue

        strand = infer_path_strand_from_paf_rows(rows, i)
        append_ctg_intype_interrupt_segment(segments, row, strand, length, i)

    attach_ctg_intype_interrupt_flanks(segments, rows)
    return segments

def append_chrom_change_break(ordered_breaks : list, left : dict, right : dict):
    if left['chrom'] != right['chrom']:
        ordered_breaks.append((left, right))

def get_karyotype_next_state(curr_node, prev_node_name : int, curr_node_name : int):
    if curr_node_name > prev_node_name:
        next_incr = curr_node[CTG_DIR]
    else:
        next_incr = '-' if curr_node[CTG_DIR] == '+' else '+'
    next_ref = curr_node[CHR_STR] if next_incr == '+' else curr_node[CHR_END]
    return [curr_node[CHR_NAM], next_incr], next_ref

def build_ctg_intype_split_bnds(weights, min_weight):
    split_weight_by_key = defaultdict(float)
    split_parent_weight = defaultdict(float)

    for path_idx, (paf_loc, key_int_list) in enumerate(paf_ans_list):
        if path_idx >= len(weights):
            break
        path_weight = float(weights[path_idx])
        if path_weight <= min_weight:
            continue

        path = import_index_path(paf_loc)
        if len(path[0]) < 4:
            path[0] = tuple([0] + list(path[0]))
        if len(path[-1]) < 4:
            path[-1] = tuple([0] + list(path[-1]))

        ctg_intype_key_by_edge = {}
        for key_int in key_int_list:
            key_type, key_value = int2key[key_int]
            if key_type == CTG_IN_TYPE:
                ctg_intype_key_by_edge[key_value] = key_int

        curr_incr = '+' if path[0][NODE_NAME][-1] == 'f' else '-'
        first_real_node = contig_data[path[1][NODE_NAME]]
        curr_chr = [first_real_node[CHR_NAM], curr_incr]
        curr_ref = first_real_node[CHR_STR] if curr_incr == '+' else first_real_node[CHR_END]

        for i in range(1, len(path) - 1):
            prev_node_name = path[i - 1][1]
            curr_node_name = path[i][1]
            if not isinstance(prev_node_name, int) or not isinstance(curr_node_name, int):
                continue

            last_node = contig_data[prev_node_name]
            curr_node = contig_data[curr_node_name]
            edge_key = (
                (path[i - 1][0], prev_node_name),
                (path[i][0], curr_node_name),
            )
            key_int = ctg_intype_key_by_edge.get(edge_key)

            interrupt_segments = []
            if key_int is not None:
                endpoint_chroms = {last_node[CHR_NAM], curr_node[CHR_NAM]}
                interrupt_segments = get_ctg_intype_interrupt_segments(key_int, endpoint_chroms)

            if not (
                path[i][CHR_CHANGE_IDX] > path[i - 1][CHR_CHANGE_IDX]
                or path[i][DIR_CHANGE_IDX] > path[i - 1][DIR_CHANGE_IDX]
                or interrupt_segments
            ):
                continue

            if interrupt_segments:
                parent_nclose = tuple(sorted([prev_node_name, curr_node_name]))
                parent_nclose_idx = nclose2idx.get(parent_nclose)
                if parent_nclose_idx is not None:
                    left_pos = last_node[CHR_END] if curr_incr == '+' else last_node[CHR_STR]
                    left_endpoint = interrupt_segments[0].get('left_flank_exit')
                    if left_endpoint is None or left_endpoint['chrom'] != curr_chr[0]:
                        left_endpoint = make_breakend_endpoint(curr_chr[0], left_pos, curr_chr[1])

                    next_chr, _next_ref = get_karyotype_next_state(
                        curr_node, prev_node_name, curr_node_name
                    )
                    right_pos = curr_node[CHR_STR] if next_chr[1] == '+' else curr_node[CHR_END]
                    right_endpoint = interrupt_segments[-1].get('right_flank_entry')
                    if right_endpoint is None or right_endpoint['chrom'] != next_chr[0]:
                        right_endpoint = make_breakend_endpoint(next_chr[0], right_pos, next_chr[1])

                    ordered_breaks = []
                    append_chrom_change_break(
                        ordered_breaks, left_endpoint, interrupt_segments[0]['entry']
                    )
                    for segment_idx in range(len(interrupt_segments) - 1):
                        append_chrom_change_break(
                            ordered_breaks,
                            interrupt_segments[segment_idx]['exit'],
                            interrupt_segments[segment_idx + 1]['entry'],
                        )
                    append_chrom_change_break(
                        ordered_breaks, interrupt_segments[-1]['exit'], right_endpoint
                    )

                    if ordered_breaks:
                        split_parent_weight[parent_nclose_idx] += path_weight
                        ctg_name = last_node[CTG_NAM]
                        for split_idx, (left, right) in enumerate(ordered_breaks, start=1):
                            split_key = (
                                parent_nclose_idx,
                                split_idx,
                                left['chrom'],
                                int(left['pos']),
                                left['strand'],
                                right['chrom'],
                                int(right['pos']),
                                right['strand'],
                                ctg_name,
                            )
                            split_weight_by_key[split_key] += path_weight

            curr_chr, curr_ref = get_karyotype_next_state(
                curr_node, prev_node_name, curr_node_name
            )

    return split_weight_by_key, split_parent_weight

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

def point_to_chrom_end_interval(chrom, point, side, contig_lengths):
    chrom_len = int(contig_lengths[chrom])
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

def build_virtual_inv_events(prefix, meandepth, contig_lengths, min_depth_N=0.0, weights=None):
    result_path = f'{prefix}/{RAW_TRANSLOCATION_RESULT_PKL}'
    if not os.path.isfile(result_path) or meandepth <= 0:
        return []

    report_by_pair = load_raw_translocation_report(prefix)
    with open(result_path, 'rb') as f:
        records = pkl.load(f)

    events = []
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
        if chrom_a not in contig_lengths or chrom_b not in contig_lengths:
            continue
        if point_a is None or point_b is None:
            continue

        if chrom_a == chrom_b:
            st, nd = sorted([int(point_a), int(point_b)])
            st = max(0, min(st, int(contig_lengths[chrom_a])))
            nd = max(0, min(nd, int(contig_lengths[chrom_a])))
            if nd > st:
                events.append((chrom_a, st, nd, expected_depth, depth_N, f'RAW_TRANSLOCATION_PAIR_{pair_id}'))
            continue

        side_a = layout_side(record, 'layout_b', 0)
        side_b = layout_side(record, 'layout_a', 1)
        interval_a = point_to_chrom_end_interval(chrom_a, point_a, side_a, contig_lengths)
        interval_b = point_to_chrom_end_interval(chrom_b, point_b, side_b, contig_lengths)
        if interval_a is not None:
            events.append((chrom_a, interval_a[0], interval_a[1], expected_depth, depth_N, f'RAW_TRANSLOCATION_PAIR_{pair_id}_A'))
        if interval_b is not None:
            events.append((chrom_b, interval_b[0], interval_b[1], expected_depth, depth_N, f'RAW_TRANSLOCATION_PAIR_{pair_id}_B'))

    return events

def _pairs_to_vcf(nclose_pairs, contig_data, contig_lengths, display_indel, amplicon_events, virtual_inv_events,
                  split_bnd_weights, out_vcf_path, nclose_cn_std):
    with vcfpy.Writer.from_path(out_vcf_path, build_skype_vcf_header(contig_lengths)) as vcf_writer:
        # 2) Translocation, Inversion 처리
        for nclose in nclose_pairs:
            # 깊이 유의성(significant_nclose) 기반 PASS/FAIL 구분 제거: balanced(copy-number-neutral)
            # translocation은 flanking depth 단차가 없어 FAIL로 강등됐지만 실제로는 진짜 junction이므로,
            # 모든 breakend에 필터를 적용하지 않는다(FILTER='.'). 깊이 값 자체는 INFO의 WEIGHT로 남는다.
            quality = 60
            filter_str = '.'

            a_idx, b_idx = nclose
            bnd_nclose_ind = nclose2idx[nclose]
            nclose_weight = nclose_cn_std[bnd_nclose_ind]

            a = contig_data[a_idx]
            b = contig_data[b_idx]

            # Order anchors by contig coordinate so a = arm1 (junction at its contig-3' end)
            # and b = arm2 (junction at its contig-5' end). The breakend adjacency is symmetric,
            # so this canonical order yields a valid reciprocal mate pair regardless of how the
            # derivative path actually traverses the contig (FOR/BAK).
            if a[CTG_STR] > b[CTG_STR]:
                a, b = b, a

            dir_a = a[CTG_DIR]
            dir_b = b[CTG_DIR]

            chr_a, chr_b = a[CHR_NAM], b[CHR_NAM]
            # Junction base depends on each anchor's SKYPE path direction (CTG_DIR):
            #   arm1 exits at its contig-3' end  -> CHR_END if '+', else CHR_STR
            #   arm2 enters at its contig-5' end -> CHR_STR if '+', else CHR_END
            pos_a = a[CHR_END] if dir_a == '+' else a[CHR_STR]
            pos_b = b[CHR_STR] if dir_b == '+' else b[CHR_END]

            # if nclose not in nclose_set:
            #     i1, i2 = conjoined_track_data[nclose]
            #     merge_mate_id_str = f';MERGE_MATEID=SKYPE.BND.{i1},SKYPE.BND.{i2}'
            #     ctg_name = f"{a[CTG_NAM]},{b[CTG_NAM]}"
            # else:
            ctg_name = a[CTG_NAM]

            write_bnd_vcf_pair(
                vcf_writer,
                f"SKYPE.BND.{bnd_nclose_ind}",
                chr_a,
                pos_a,
                dir_a,
                chr_b,
                pos_b,
                dir_b,
                nclose_weight / N,
                ctg_name,
                quality=quality,
                filter_str=filter_str,
            )

        for split_key, split_weight in sorted(split_bnd_weights.items()):
            (
                parent_nclose_idx,
                split_idx,
                chr_a,
                pos_a,
                dir_a,
                chr_b,
                pos_b,
                dir_b,
                ctg_name,
            ) = split_key
            write_bnd_vcf_pair(
                vcf_writer,
                f"SKYPE.BND.{parent_nclose_idx}.{split_idx}",
                chr_a,
                pos_a,
                dir_a,
                chr_b,
                pos_b,
                dir_b,
                split_weight / N,
                ctg_name,
                merge_mate_ids=[
                    f"SKYPE.BND.{parent_nclose_idx}_1",
                    f"SKYPE.BND.{parent_nclose_idx}_2",
                ],
            )

        # 3) indel 처리
        #   front_jump -> DEL (reference 구간이 어셈블리에서 사라짐)
        #   back_jump  -> DUP (reference 구간 [lo,hi]가 추가 카피로 재방문됨 = tandem duplication)
        # depth 기반 검출은 novel insertion이 아니라 reference 구간의 copy 변화를 보므로,
        # END/SVLEN은 reference span을 그대로 사용하고 POS<=END가 되도록 정규화한다.
        dup_counter = 1
        del_dounter = 1

        for chrom, indel_list in display_indel.items():
            for indel in indel_list:
                indel_type, start, end, w, _, indel_idx = indel
                lo, hi = min(start, end), max(start, end)
                if indel_type == 'd':
                    sv_id = f"SKYPE.DEL.{del_dounter}"
                    svlen = -(hi - lo)
                    alt = "<DEL>"

                    del_dounter += 1
                elif indel_type == 'i':
                    sv_id = f"SKYPE.DUP.{dup_counter}"
                    svlen = hi - lo
                    alt = "<DUP>"

                    dup_counter += 1
                else:
                    assert(False)

                ctg_name = (
                    f"INDEL_INDEX_{indel_idx}"
                    if isinstance(indel_idx, int)
                    else str(indel_idx)
                )
                vcf_writer.write_record(make_vcf_record(
                    chrom,
                    lo,
                    sv_id,
                    "N",
                    vcfpy.SymbolicAllele(alt[1:-1]),
                    60,
                    ".",
                    ordered_mapping(
                        ("SVTYPE", alt[1:-1]),
                        ("END", int(hi)),
                        ("SVLEN", int(svlen)),
                        ("WEIGHT", round(w, 2)),
                        ("CTG_NAME", str(ctg_name)),
                    ),
                ))

        for amplicon_counter, (chrom, st, nd, _, depth_N, amplicon_idx) in enumerate(amplicon_events, start=1):
            pos = max(1, int(min(st, nd)))
            end = max(pos, int(max(st, nd)))
            vcf_writer.write_record(make_vcf_record(
                chrom,
                pos,
                f"SKYPE.AMP.{amplicon_counter}",
                "N",
                vcfpy.SymbolicAllele("DUP"),
                60,
                ".",
                ordered_mapping(
                    ("SVTYPE", "DUP"),
                    ("END", end),
                    ("SVLEN", end - pos),
                    ("WEIGHT", round(depth_N, 2)),
                    ("CTG_NAME", f"AMPLICON_INDEX_{amplicon_idx}"),
                    ("SVCLASS", "AMPLICON"),
                ),
            ))

        for inv_counter, (chrom, st, nd, _, depth_N, name) in enumerate(virtual_inv_events, start=1):
            pos = max(1, int(st))
            end = max(pos, int(nd))
            vcf_writer.write_record(make_vcf_record(
                chrom,
                pos,
                f"SKYPE.VINV.{inv_counter}",
                "N",
                vcfpy.SymbolicAllele("INV"),
                60,
                ".",
                ordered_mapping(
                    ("SVTYPE", "INV"),
                    ("END", end),
                    ("SVLEN", end - pos),
                    ("WEIGHT", round(depth_N, 2)),
                    ("CTG_NAME", str(name)),
                ),
            ))

parser = argparse.ArgumentParser(description="SKYPE depth analysis")

parser.add_argument("censat_bed_path", 
                    help="Path to the censat repeat information file.")

parser.add_argument("ppc_paf_file_path", 
                    help="Path to the preprocessed PAF file.")

parser.add_argument("main_stat_loc", 
                    help="Cancer coverage location file")

parser.add_argument("telomere_bed_path", 
                    help="Path to the telomere information file.")

parser.add_argument("reference_fai_path", 
                    help="Path to the chromosome information file.")

parser.add_argument("reference_cytobands_path", 
                    help="Path to the cytoband information file.")

parser.add_argument("prefix", 
                    help="Pefix for pipeline")

parser.add_argument("-t", "--thread", 
                    help="Number of thread", type=int)

parser.add_argument("--progress", 
                    help="Show progress bar", action='store_true')

args = parser.parse_args()

# t = "31_depth_analysis.py public_data/chm13v2.0_censat_v2.1.m.bed /home/hyunwoo/ACCtools-pipeline/90_skype_run/Caki-1/20_alignasm/Caki-1.ctg.aln.paf.ppc.paf /home/hyunwoo/ACCtools-pipeline/90_skype_run/Caki-1/01_depth/Caki-1_normalized.win.stat.gz public_data/chm13v2.0_telomere.bed public_data/chm13v2.0.fa.fai public_data/chm13v2.0_cytobands_allchrs.bed 30_skype_pipe/Caki-1_18_22_40 -t 128"
# args = parser.parse_args(t.split()[1:])

bed_data = import_bed(args.censat_bed_path)

PREFIX = args.prefix
THREAD = args.thread
CHROMOSOME_INFO_FILE_PATH = args.reference_fai_path
main_stat_loc = args.main_stat_loc
TELOMERE_INFO_FILE_PATH = args.telomere_bed_path
PREPROCESSED_PAF_FILE_PATH = args.ppc_paf_file_path
pipeline_mode_config = load_pipeline_mode(PREFIX)

RATIO_OUTLIER_FOLDER = f"{PREFIX}/11_ref_ratio_outliers/"
front_contig_path = RATIO_OUTLIER_FOLDER+"front_jump/"
back_contig_path = RATIO_OUTLIER_FOLDER+"back_jump/"
ecdna_contig_path = RATIO_OUTLIER_FOLDER + "ecdna/"
output_folder = f'{PREFIX}/21_pat_depth'

TELO_CONNECT_NODES_INFO_PATH = PREFIX+"/telomere_connected_list.txt"

contig_data = import_data2(PREPROCESSED_PAF_FILE_PATH)
telo_connected_node = extract_telomere_connect_contig(TELO_CONNECT_NODES_INFO_PATH)

with open(f'{PREFIX}/path_data.pkl', 'rb') as f:
    path_list_dict = pkl.load(f)

B = np.load(f'{PREFIX}/B.npy')

df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
df = df.query('chr != "chrM"')

meandepth = np.median(df['meandepth'])

BREAKEND_REMARKABLE_CN = meandepth * BREAKEND_REMARKABLE_CN_RATIO
TELOMERE_REMARKABLE_CN = meandepth * TELOMERE_REMARKABLE_CN_RATIO
N = meandepth / 2
NCLOSE_SIM_DIFF_THRESHOLD = 0.1 * N

chr_order_list = extract_groups(list(df['chr']))

chr_filt_st_list = []
chr_no_filt_st_list = []

chr_filt_idx_list = []
chr_no_filt_idx_list = []

for ind, l in enumerate(df.itertuples(index=False)):
    flag = True
    if l.meandepth > ABS_MAX_COVERAGE_RATIO * meandepth:
        flag = False
    for i in bed_data:
        if l.chr == i:
            for j in bed_data[i]:
                if inclusive_checker(j, (l.st, l.nd)):
                    flag = False
                    break
        if not flag:
            break
    
    if flag:
        chr_filt_st_list.append((l.chr, l.st))
        chr_filt_idx_list.append(ind)
    else:
        chr_no_filt_st_list.append((l.chr, l.st))
        chr_no_filt_idx_list.append(ind)

filter_len = len(chr_filt_st_list)

def get_vec_from_stat_loc(stat_loc_):
    df = pd.read_csv(stat_loc_, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
    df = df.query('chr != "chrM"')

    chr_st_data = dict()
    for l in df.itertuples(index=False):
         chr_st_data[(l.chr, l.st)] = l.meandepth
    
    v = []
    for cs in chr_filt_st_list:
        if cs not in chr_st_data:
            v.append(0)
        else:
            v.append(chr_st_data[cs])
    
    for cs in chr_no_filt_st_list:
        if cs not in chr_st_data:
            v.append(0)
        else:
            v.append(chr_st_data[cs])
    
    return np.asarray(v, dtype=np.float32)

def get_vec_from_ki(ki):
    stat_loc_ = f'{output_folder}/{ki}.win.stat.gz'
    df = pd.read_csv(stat_loc_, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
    df = df.query('chr != "chrM"')

    chr_st_data = dict()
    for l in df.itertuples(index=False):
         chr_st_data[(l.chr, l.st)] = l.meandepth
    
    v = []
    for cs in chr_filt_st_list:
        if cs not in chr_st_data:
            v.append(0)
        else:
            v.append(chr_st_data[cs])
    
    for cs in chr_no_filt_st_list:
        if cs not in chr_st_data:
            v.append(0)
        else:
            v.append(chr_st_data[cs])
    
    return ki, np.asarray(v, dtype=np.float32)

with open(f'{PREFIX}/contig_pat_vec_data.pkl', 'rb') as f:
    paf_ans_list, key_list, int2key, _ = pkl.load(f)

raw_path_list = []
for i, ll in paf_ans_list:
    il = i.split("/")
    cnt = il[-1].split('.')[0]
    chr2chr = il[-2]
    raw_path_list.append(f'{PREFIX}/00_raw/{chr2chr}/{cnt}.index.txt')
rpll = len(raw_path_list)

with open(f'{PREFIX}/tot_loc_list.pkl', 'rb') as f:
    tot_loc_list = pkl.load(f)

with open(f'{PREFIX}/cen_fragment_data.pkl', 'rb') as f:
    cen_fragment_meta = pkl.load(f)

grouped_data = defaultdict(lambda: {"positions": [], "values": []})
for i, (chrom, pos) in enumerate(chr_filt_st_list):
    grouped_data[chrom]["positions"].append(pos)
    grouped_data[chrom]["values"].append(B[i])

# 필터링 파라미터 (데이터 특성에 따라 조정)
cutoff_frequency = 0.1  # 컷오프 주파수 (예시 값)
fs = 1.0                 # 샘플링 주파수 (데이터 포인트 간 간격이 1이라고 가정)
order = 3                # 필터 차수

# 각 염색체별로 필터링 수행
filtered_values_list = []
for chrom, data in grouped_data.items():
    positions = np.array(data["positions"])
    values = np.array(data["values"])
    
    # 만약 위치가 정렬되어 있지 않다면 정렬 (예시에서는 이미 정렬된 상태)
    sort_idx = np.argsort(positions)
    positions = positions[sort_idx]
    values = values[sort_idx]
    
    # 고주파 필터 적용: 기본 CN 트렌드를 제거하여 노이즈 성분 추출
    filtered_values = highpass_filter(values, cutoff=cutoff_frequency, fs=fs, order=order)

    filtered_values_list.append(filtered_values)

noise_array = np.hstack(filtered_values_list)
amplitude = np.std(noise_array)

chr_len = find_chr_len(CHROMOSOME_INFO_FILE_PATH)
del chr_len['chrM']

with open(f'{PREFIX}/nclose_chunk_data.pkl', 'rb') as f:
    nclose_nodes_pkl, _, _ = pkl.load(f)

nclose_list = []
for vl in nclose_nodes_pkl.values():
    for v in vl:
        nclose_list.append(v)

nclose_set = set(nclose_list)

nclose_idx = len(nclose_list)
nclose2idx = dict(zip(nclose_list, range(1, len(nclose_list) + 1)))
idx2nclose = dict(zip(range(1, len(nclose_list) + 1), nclose_list))

nclose_str_pos = dict()
for k, v in idx2nclose.items():
    nclose_str_pos[k] = (contig_data[v[0]][CHR_STR], contig_data[v[1]][CHR_STR])

telo_node_set = set()

for chr_dir, node_id in telo_connected_node:
    telo_node_set.add(node_id)

path_nclose_usage = []
path_telo_usage = []
for i, raw_path in enumerate(raw_path_list):
    using_nclose = Counter()
    using_telo = Counter()
    path = import_index_path(raw_path)
    s = 1
    while s < len(path)-2:
        nclose_cand = tuple(sorted([path[s][1], path[s+1][1]]))
        if nclose_cand in nclose_set:
            using_nclose[nclose2idx[nclose_cand]]+=1
            s+=2
        else:
            s+=1
    if path[1][1] in telo_node_set:
        using_telo[path[1][1]]+=1
    if path[len(path)-2][1] in telo_node_set:
        using_telo[path[len(path)-2][1]]+=1

    path_nclose_usage.append(using_nclose)
    path_telo_usage.append(using_telo)          

fclen = len(glob.glob(front_contig_path + "*"))
bclen = len(glob.glob(back_contig_path + "*"))
eclen = len(glob.glob(ecdna_contig_path + "*"))

n = len(paf_ans_list) + fclen // 4 + bclen // 4 + eclen // 2 + len(cen_fragment_meta)

with open(f'{PREFIX}/ecdna_circuit_data.pkl', 'rb') as f:
    ecdna_circuit, _ = pkl.load(file=f)

def build_amplicon_events(weights, min_weight=BREAKEND_REMARKABLE_CN):
    events = []
    ecdna_offset = rpll + fclen // 4 + bclen // 4
    for i in range(rpll, len(weights)):
        paf_loc = tot_loc_list[i]
        if paf_loc.split('/')[-3] == '12_cent_fragment':
            continue
        if paf_loc.split('/')[-2] != 'ecdna':
            continue

        raw_weight = weights[i]
        if raw_weight <= min_weight:
            continue

        ecdna_ind = i - ecdna_offset
        if ecdna_ind < 0 or ecdna_ind >= len(ecdna_circuit):
            logging.warning(f"Skip amplicon with invalid ecdna index: {ecdna_ind}")
            continue

        circuit = ecdna_circuit[ecdna_ind]
        if len(circuit) == 0:
            continue

        ctg_st_list = []
        ctg_nd_list = []
        for ctg_ind in circuit:
            ctg_st_list.append(contig_data[ctg_ind][CHR_STR])
            ctg_nd_list.append(contig_data[ctg_ind][CHR_END])

        ctg_ind = circuit[0]
        events.append((
            contig_data[ctg_ind][CHR_NAM],
            min(ctg_st_list),
            max(ctg_nd_list),
            raw_weight,
            raw_weight / N * 2,
            ecdna_ind,
        ))

    return events


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


def build_aggregated_indel_events(weights, min_weight=0.0):
    indel_events = {}

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
            source=f"INDEL_INDEX_{i - rpll}"
        )

    type4_events, type4_weights, _ = summarize_type4_indel_graph_usage(
        PREFIX, paf_ans_list, weights, import_index_path
    )
    for event_key, raw_weight in type4_weights.items():
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
        if event["weight"] > min_weight
    ]


def input_vcf_record_key(line_no, rec_id):
    if rec_id in ("", "."):
        return f"VCF_RECORD_{line_no}"
    return rec_id


def sanitize_vcf_id(value):
    value = str(value) if value not in (None, "") else "unknown"
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value)


def sanitize_vcf_info_token(value):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "unknown"))


def format_skype_cn(value):
    value = float(value)
    if abs(value) < 1e-12:
        value = 0.0
    return f"{value:.6g}"


SKYPE_ANNOTATION_INFO_FIELDS = (
    "SKYPE_CN",
    "SKYPE_STATUS",
    "SKYPE_CN_DETAIL",
    "SKYPE_STATUS_DETAIL",
)


def add_info_header_line_if_missing(header, info_id, number, type_, description):
    if not header.has_header_line("INFO", info_id):
        add_info_header_line(header, info_id, number, type_, description)


def add_skype_annotation_headers(header, include_detail_records):
    add_info_header_line_if_missing(
        header,
        "SKYPE_CN",
        1,
        "Float",
        "SKYPE normalized copy-number support",
    )
    add_info_header_line_if_missing(
        header,
        "SKYPE_STATUS",
        1,
        "String",
        "SKYPE VCF input mode evaluation status",
    )
    if include_detail_records:
        add_info_header_line_if_missing(
            header,
            "SKYPE_CN_DETAIL",
            1,
            "String",
            "Pipe-delimited SKYPE normalized copy-number support detail. "
            "Only emitted for records with side-specific detail; for INV records the order is left|right",
        )
        add_info_header_line_if_missing(
            header,
            "SKYPE_STATUS_DETAIL",
            1,
            "String",
            "Pipe-delimited SKYPE evaluation status detail. "
            "Only emitted for records with side-specific detail; for INV records the order is left|right",
        )


def clear_skype_annotation_info(record):
    for key in SKYPE_ANNOTATION_INFO_FIELDS:
        record.INFO.pop(key, None)


def read_input_vcf_records(vcf_path):
    records = []
    sanitized_to_keys = defaultdict(list)
    with open(vcf_path, "rt") as f:
        for line_no, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 8:
                records.append({
                    "line_no": line_no,
                    "cols": cols,
                    "key": f"VCF_RECORD_{line_no}",
                    "malformed": True,
                })
                continue
            key = input_vcf_record_key(line_no, cols[2])
            records.append({
                "line_no": line_no,
                "cols": cols,
                "key": key,
                "malformed": False,
            })
            sanitized_to_keys[sanitize_vcf_id(key)].append(key)
    return records, sanitized_to_keys


def load_vcf_skip_reasons():
    skip_by_key = {}
    skip_by_line = {}
    skipped_path = f"{PREFIX}/vcf_mode_skipped_records.tsv"
    if not os.path.isfile(skipped_path):
        return skip_by_key, skip_by_line
    with open(skipped_path, "rt") as f:
        next(f, None)
        for line in f:
            cols = line.rstrip("\n").split("\t", 4)
            if len(cols) < 4:
                continue
            line_no_text, rec_id, _svtype, reason = cols[:4]
            try:
                line_no = int(line_no_text)
            except ValueError:
                continue
            skip_by_line[line_no] = reason
            key = input_vcf_record_key(line_no, rec_id)
            skip_by_key[key] = reason
    return skip_by_key, skip_by_line


def resolve_sanitized_vcf_id(sanitized_id, sanitized_to_keys):
    return list(sanitized_to_keys.get(sanitized_id, [sanitized_id]))


VCF_DETAIL_SIDE_ORDER = ("left", "right")


def vcf_record_keys_from_node_name(node_name, sanitized_to_keys):
    node_name = str(node_name)
    if node_name.startswith("vcf_inv_"):
        body = node_name[len("vcf_inv_"):]
        for suffix in ("_left", "_right"):
            if body.endswith(suffix):
                body = body[:-len(suffix)]
                break
        return resolve_sanitized_vcf_id(body, sanitized_to_keys)

    if node_name.startswith("vcf_bnd_"):
        body = node_name[len("vcf_bnd_"):]
        sanitized_ids = sorted(sanitized_to_keys, key=len, reverse=True)
        for sanitized_id in sanitized_ids:
            if body == sanitized_id:
                return resolve_sanitized_vcf_id(sanitized_id, sanitized_to_keys)
            prefix = sanitized_id + "_"
            if not body.startswith(prefix):
                continue
            mate_id = body[len(prefix):]
            if mate_id in sanitized_to_keys:
                return (
                    resolve_sanitized_vcf_id(sanitized_id, sanitized_to_keys)
                    + resolve_sanitized_vcf_id(mate_id, sanitized_to_keys)
                )
        return resolve_sanitized_vcf_id(body, sanitized_to_keys)

    return []


def vcf_record_key_details_from_node_name(node_name, sanitized_to_keys):
    node_name = str(node_name)
    if not node_name.startswith("vcf_inv_"):
        return []

    body = node_name[len("vcf_inv_"):]
    side = None
    for suffix in ("_left", "_right"):
        if body.endswith(suffix):
            body = body[:-len(suffix)]
            side = suffix[1:]
            break

    if side not in VCF_DETAIL_SIDE_ORDER:
        return []
    return [
        (key, side)
        for key in resolve_sanitized_vcf_id(body, sanitized_to_keys)
    ]


def event_vcf_record_keys(event):
    keys = []
    for value in event.get("merged_vcf_ids", []):
        if value not in (None, "", True):
            keys.append(str(value))
    for key in ("vcf_id", "mate_id"):
        value = event.get(key)
        if value not in (None, "", True):
            keys.append(str(value))
    seen = set()
    unique_keys = []
    for key in keys:
        if key in seen:
            continue
        seen.add(key)
        unique_keys.append(key)
    return unique_keys


def add_vcf_cn(cn_by_key, evaluated_keys, record_keys, cn):
    for key in record_keys:
        if key in (None, "", True):
            continue
        key = str(key)
        cn_by_key[key] = max(float(cn), cn_by_key.get(key, 0.0))
        evaluated_keys.add(key)


def get_vcf_type4_outlier_event(event_type, outlier_idx, vcf_type4_outlier_index):
    return (
        vcf_type4_outlier_index.get((event_type, outlier_idx))
        or vcf_type4_outlier_index.get(f"{event_type}:{outlier_idx}")
    )


def build_vcf_annotation_cn(weights, sanitized_to_keys):
    cn_by_key = {}
    evaluated_keys = set()
    detail_cn_by_key = defaultdict(dict)
    detail_evaluated_sides_by_key = defaultdict(set)
    detail_known_sides_by_key = defaultdict(set)
    n_unit = float(N) if N else 0.0

    def weight_to_cn(weight):
        if n_unit == 0:
            return 0.0
        return max(0.0, float(weight) / n_unit)

    vcf_type4_outlier_index_path = f"{PREFIX}/{VCF_TYPE4_OUTLIER_INDEX_PKL}"
    if os.path.isfile(vcf_type4_outlier_index_path):
        with open(vcf_type4_outlier_index_path, "rb") as f:
            vcf_type4_outlier_index = pkl.load(f)
    else:
        vcf_type4_outlier_index = {}

    for i, paf_loc in enumerate(tot_loc_list[:len(weights)]):
        event_type = os.path.basename(os.path.dirname(paf_loc))
        if event_type not in {"front_jump", "back_jump"}:
            continue
        basename = os.path.basename(paf_loc)
        if not basename.endswith("_base.paf"):
            continue
        try:
            outlier_idx = int(basename.split("_", 1)[0])
        except ValueError:
            continue
        event = get_vcf_type4_outlier_event(event_type, outlier_idx, vcf_type4_outlier_index)
        if event is None:
            continue
        add_vcf_cn(cn_by_key, evaluated_keys, event_vcf_record_keys(event), weight_to_cn(weights[i]))

    for nclose in idx2nclose.values():
        for node_idx in nclose:
            if 0 <= int(node_idx) < len(contig_data):
                for key, side in vcf_record_key_details_from_node_name(
                    contig_data[int(node_idx)][CTG_NAM], sanitized_to_keys
                ):
                    detail_known_sides_by_key[key].add(side)

    nclose_cn = defaultdict(float)
    for i, ctr in enumerate(path_nclose_usage):
        if i >= len(weights):
            break
        for nclose_idx, count in ctr.items():
            nclose_cn[nclose_idx] += count * float(weights[i])

    for nclose_idx, raw_weight in nclose_cn.items():
        nclose = idx2nclose.get(nclose_idx)
        if nclose is None:
            continue
        record_keys = []
        detail_entries = []
        for node_idx in nclose:
            if 0 <= int(node_idx) < len(contig_data):
                record_keys.extend(
                    vcf_record_keys_from_node_name(contig_data[int(node_idx)][CTG_NAM], sanitized_to_keys)
                )
                detail_entries.extend(
                    vcf_record_key_details_from_node_name(
                        contig_data[int(node_idx)][CTG_NAM], sanitized_to_keys
                    )
                )
        if record_keys:
            cn = weight_to_cn(raw_weight)
            add_vcf_cn(cn_by_key, evaluated_keys, record_keys, cn)
            for key, side in detail_entries:
                detail_known_sides_by_key[key].add(side)
                detail_cn_by_key[key][side] = max(
                    float(cn),
                    float(detail_cn_by_key[key].get(side, 0.0)),
                )
                detail_evaluated_sides_by_key[key].add(side)

    return (
        cn_by_key,
        evaluated_keys,
        detail_cn_by_key,
        detail_evaluated_sides_by_key,
        detail_known_sides_by_key,
    )


def build_vcf_detail_fields(
    key,
    cn,
    status,
    record,
    detail_cn_by_key,
    detail_evaluated_sides_by_key,
    detail_known_sides_by_key,
    skip_by_key,
    skip_by_line,
):
    known_sides = detail_known_sides_by_key.get(key, set())
    if not known_sides:
        return None, None

    cn_detail = []
    status_detail = []
    skipped_reason = (
        skip_by_key.get(key)
        or skip_by_line.get(record["line_no"])
        or "NOT_EVALUATED"
    )
    skipped_status = "SKIPPED_" + sanitize_vcf_info_token(skipped_reason)
    evaluated_sides = detail_evaluated_sides_by_key.get(key, set())
    side_cn = detail_cn_by_key.get(key, {})
    for side in VCF_DETAIL_SIDE_ORDER:
        cn_detail.append(format_skype_cn(side_cn.get(side, 0.0)))
        status_detail.append("EVALUATED" if side in evaluated_sides else skipped_status)

    return "|".join(cn_detail), "|".join(status_detail)


def write_annotated_input_vcf(weights):
    vcf_input_path = pipeline_mode_config.get("vcf_input_path")
    if not vcf_input_path:
        raise FileNotFoundError("VCF input mode is missing vcf_input_path in pipeline_mode.pkl")
    if not os.path.isfile(vcf_input_path):
        raise FileNotFoundError(f"VCF input file does not exist: {vcf_input_path}")

    records, sanitized_to_keys = read_input_vcf_records(vcf_input_path)
    (
        cn_by_key,
        evaluated_keys,
        detail_cn_by_key,
        detail_evaluated_sides_by_key,
        detail_known_sides_by_key,
    ) = build_vcf_annotation_cn(weights, sanitized_to_keys)
    skip_by_key, skip_by_line = load_vcf_skip_reasons()

    has_detail_records = any(detail_known_sides_by_key.values())

    with vcfpy.Reader.from_path(vcf_input_path) as reader:
        header = reader.header.copy()
        add_skype_annotation_headers(header, has_detail_records)
        metadata_iter = iter(records)
        with vcfpy.Writer.from_path(f"{PREFIX}/SV_benchmark_result.vcf", header) as writer:
            for vcf_record in reader:
                try:
                    record = next(metadata_iter)
                except StopIteration:
                    raise ValueError(
                        "VCF metadata scan produced fewer records than vcfpy.Reader"
                    )

                key = record["key"]
                cn = cn_by_key.get(key, 0.0)
                if key in evaluated_keys:
                    status = "EVALUATED"
                else:
                    reason = skip_by_key.get(key) or skip_by_line.get(record["line_no"]) or "NOT_EVALUATED"
                    status = "SKIPPED_" + sanitize_vcf_info_token(reason)

                cn_detail, status_detail = build_vcf_detail_fields(
                    key,
                    cn,
                    status,
                    record,
                    detail_cn_by_key,
                    detail_evaluated_sides_by_key,
                    detail_known_sides_by_key,
                    skip_by_key,
                    skip_by_line,
                )
                clear_skype_annotation_info(vcf_record)
                vcf_record.INFO["SKYPE_CN"] = format_skype_cn(cn)
                vcf_record.INFO["SKYPE_STATUS"] = status
                if cn_detail is not None and status_detail is not None:
                    vcf_record.INFO["SKYPE_CN_DETAIL"] = cn_detail
                    vcf_record.INFO["SKYPE_STATUS_DETAIL"] = status_detail
                writer.write_record(vcf_record)

        try:
            next(metadata_iter)
        except StopIteration:
            pass
        else:
            raise ValueError(
                "VCF metadata scan produced more records than vcfpy.Reader"
            )


def draw_circos_plot(fig_prefix=''):
    weights = np.load(f'{PREFIX}/weight{fig_prefix}.npy')
    smooth_B = rebin_dataframe_B(df, 10)[chr_filt_idx_list + chr_no_filt_idx_list]
    predicted_B = np.maximum(np.load(f'{PREFIX}/predict_B{fig_prefix}.npy'), 0)
    diff = predicted_B - smooth_B

    threshold = amplitude * scipy.stats.norm.ppf(0.975)
    miss_B = np.abs(diff) > threshold
    over_B = np.asarray([False] * len(chr_filt_st_list) + [True] * len(chr_no_filt_st_list))

    color_label = np.ones_like(B)
    color_label[miss_B] = 3
    color_label[over_B] = 5

    if fig_prefix == '':
        msg_prefix = 'Fail'
    else:
        msg_prefix = f'{fig_prefix[1:].capitalize()} fail'

    logging.info(f'{msg_prefix} ratio : {round(sum(color_label[miss_B] == 3) / len(B) * 100, 3)}%')

    nclose_cn = defaultdict(float)
    for i, ctr in enumerate(path_nclose_usage):
        for j, v in ctr.items():
            nclose_cn[j] += v*weights[i]

    split_bnd_weights, split_parent_weight = build_ctg_intype_split_bnds(
        weights, BREAKEND_REMARKABLE_CN
    )
    adjusted_nclose_cn = defaultdict(float, nclose_cn)
    for parent_nclose_idx, split_weight in split_parent_weight.items():
        adjusted_nclose_cn[parent_nclose_idx] = max(
            0.0, adjusted_nclose_cn[parent_nclose_idx] - split_weight
        )
    
    # with open(f"{PREFIX}/nclose_cn.txt", "wt") as f:
    #     for k, v in nclose_cn.items():
    #         if v > BREAKEND_REMARKABLE_CN:
    #             curr_nclose = reverse_nclose_dict[k]
    #             print(contig_data[curr_nclose[0]], "\n", contig_data[curr_nclose[1]], "\n", v, "\n", file=f)


    rdf = rebin_dataframe(df, 2)

    fragment_depth_per_chrom = {}
    for i, paf_loc in enumerate(tot_loc_list):
        if paf_loc.split('/')[-3] == '12_cent_fragment':
            chrom = paf_loc.split('/')[-2]
            side = paf_loc.split('/')[-1].split('.')[0]
            fragment_depth_per_chrom[chrom] = (side, float(weights[i]))

    circos = Circos(chr_len, space=3)
    circos.add_cytoband_tracks((95, 100), args.reference_cytobands_path)

    target_color = ['blue', 'red', 'gray']

    for sector in circos.sectors:
        sector.text(sector.name, r=107, size=8)
        sector.get_track("cytoband").xticks_by_interval(
            25 * M,
            label_size=6,
            tick_length=1,
            label_orientation="vertical",
            label_formatter=lambda v: f"{int(v / M)}",
        )

        sector.get_track("cytoband").xticks_by_interval(
            5 * M,
            tick_length=0.5,
            label_formatter=lambda v: "",
        )

        tdf = rdf.query(f'chr == "{sector.name}"')

        tchr_st_list = []
        tchr_idx_list = []
        for i, v in enumerate(chr_filt_st_list + chr_no_filt_st_list):
            if v[0] == sector.name:
                tchr_st_list.append(v[1])
                tchr_idx_list.append(i)
        
        tchr_st_list = np.array(tchr_st_list)
        tchr_idx_list = np.array(tchr_idx_list)

        sort_ind = np.argsort(tchr_st_list)
        tchr_st_list = tchr_st_list[sort_ind]
        tchr_idx_list = tchr_idx_list[sort_ind]

        cn_track = sector.add_track((60, 93))
        cn_track.axis()
        cn_track.grid(y_grid_num=7, zorder=-2)
        cn_track.line([0, max(tdf['st'])], [0.5 * meandepth] * 2, color='blue', vmax=ABS_MAX_COVERAGE_RATIO * meandepth, zorder=-1, alpha=0.5)
        cn_track.line([0, max(tdf['st'])], [2 * meandepth] * 2, color='red', vmax=ABS_MAX_COVERAGE_RATIO * meandepth, zorder=-1, alpha=0.5)
        line_track_plot(cn_track, x=tchr_st_list, y=predicted_B[tchr_idx_list], vmax=ABS_MAX_COVERAGE_RATIO * meandepth, color='green', linewidth=0.5, zorder=10)
        line_track_plot(cn_track, x=tdf['st'].to_list(), y=tdf['meandepth'].to_list(), vmax=ABS_MAX_COVERAGE_RATIO * meandepth, color='black', linewidth=0.2)
        line_track_circos_color(cn_track, x=tdf['st'].to_list(), y1=tdf['meandepth'].to_list(), vmax=ABS_MAX_COVERAGE_RATIO * meandepth,
                                target_color=['#6baed6', '#969696', '#fb6a4a'],
                                target_value=[0.5 * meandepth, 2 * meandepth])

        if sector.name in fragment_depth_per_chrom:
            side, frag_w = fragment_depth_per_chrom[sector.name]
            info = cen_fragment_meta[sector.name]
            if side == 'right':
                st_lo, st_hi = info['mid'], info['chr_len']
            else:
                st_lo, st_hi = 0, info['mid']
            cn_track.line(
                [st_lo, st_lo, st_hi, st_hi],
                [0, frag_w, frag_w, 0],
                color='purple',
                linestyle='--',
                linewidth=0.7,
                vmax=ABS_MAX_COVERAGE_RATIO * meandepth,
                zorder=11,
            )

        color_track = sector.add_track((53, 58))
        color_track.axis()
        line_track_circos_color(color_track, x=tchr_st_list, y1=color_label[tchr_idx_list], vmax=0.5,
                                target_color=target_color,
                                target_value=[2, 4])  

    telo_cn = defaultdict(float)
    for i, ctr in enumerate(path_telo_usage):
        for j, v in ctr.items():
            telo_cn[j] += v*weights[i]

    inv_val_list = []
    transloc_val_list = []
    bnd_cn_data = []

    for k in nclose_cn:
        v = adjusted_nclose_cn[k]
        pos1, pos2 = nclose_str_pos[k]
        idx1, idx2 = idx2nclose[k]
        chr_nam1 = contig_data[idx1][CHR_NAM]
        chr_nam2 = contig_data[idx2][CHR_NAM]

        event_type = 'breakend'
        if chr_nam1 != chr_nam2 and v > BREAKEND_REMARKABLE_CN:
            transloc_val_list.append(v / meandepth * 2)
        elif chr_nam1 == chr_nam2 and contig_data[idx1][CTG_DIR] != contig_data[idx2][CTG_DIR] and v > BREAKEND_REMARKABLE_CN:
            inv_val_list.append(v / meandepth * 2)
            event_type = 'inversion'
        
        if v > BREAKEND_REMARKABLE_CN:
            bnd_cn_data.append([(chr_nam1, pos1), (chr_nam2, pos2), v, event_type])

    for split_key, v in sorted(split_bnd_weights.items()):
        (
            _parent_nclose_idx,
            _split_idx,
            chr_a,
            pos_a,
            dir_a,
            chr_b,
            pos_b,
            dir_b,
            _ctg_name,
        ) = split_key
        if v <= BREAKEND_REMARKABLE_CN:
            continue

        event_type = 'breakend'
        if chr_a != chr_b:
            transloc_val_list.append(v / meandepth * 2)
        elif dir_a != dir_b:
            inv_val_list.append(v / meandepth * 2)
            event_type = 'inversion'
        bnd_cn_data.append([(chr_a, pos_a), (chr_b, pos_b), v, event_type])

    indel_val_list = []
    for event in build_aggregated_indel_events(weights, BREAKEND_REMARKABLE_CN):
        bnd_cn_data.append([
            (event["chrom"], event["st"]),
            (event["chrom"], event["nd"]),
            event["weight"],
            'indel',
        ])
        indel_val_list.append(event["weight"] / meandepth * 2)

    for chrom, st, nd, raw_weight, _, _ in build_amplicon_events(weights):
        bnd_cn_data.append([(chrom, st), (chrom, nd), raw_weight, 'amplicon'])

    virtual_inv_events = build_virtual_inv_events(
        PREFIX, meandepth, chr_len, min_depth_N=BREAKEND_REMARKABLE_CN / N, weights=weights
    )
    for chrom, st, nd, expected_depth, depth_N, _ in virtual_inv_events:
        bnd_cn_data.append([(chrom, st), (chrom, nd), expected_depth, 'virtual_inv'])
        inv_val_list.append(depth_N)

    a = sorted(list(telo_cn.items()), key = lambda t:t[1])
    telo_zorder_dict = {}
    for i, v in enumerate(a):
        telo_zorder_dict[v[0]] = i

    off = 0.03
    lower = 0
    upper = meandepth * 2

    cmap = sns.color_palette("rocket_r", as_cmap=True)

    bnd_cn_data = sorted(bnd_cn_data, key=lambda t: t[2])

    line_style_by_event = {
        'breakend': '-',
        'inversion': '-.',
        'indel': '--',
        'amplicon': (0, (5, 2)),
        'virtual_inv': (0, (2, 2)),
    }

    for i, (bnd_loc1, bnd_loc2, cn, event_type) in enumerate(bnd_cn_data):
        norm_cn = np.clip((cn - lower) / (upper - lower), 0, 1)

        color = cmap(norm_cn)
        linestyle = line_style_by_event[event_type]
        
        if cn >= 0:
            circos.link_line(bnd_loc1, bnd_loc2, color=color, linestyle=linestyle, lw=1, zorder=i)

    circos.colorbar(vmin=0, vmax=2, bounds=(1.01 + off, 0.825, 0.02, 0.1), cmap=cmap,
                    label='Breakend CN', orientation='vertical',
                    colorbar_kws=dict(format=mpl.ticker.FixedFormatter(['0', '2N', '4N']),
                                    ticks=mpl.ticker.FixedLocator([0, 1, 2])),
                    label_kws=dict(fontsize=9,),
                    tick_kws=dict(labelsize=8,))

    telo_data = import_telo_data(TELOMERE_INFO_FILE_PATH, chr_len)
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
            telo_len = min(telo_len, distance_checker_tuple(tuple(telo_bed), (contig_data[node_id][CHR_STR], contig_data[node_id][CHR_END])))
        chr_fb_len_dict[chr_dir].append((node_id, telo_len, chr_dir))

    telo_len_data = []
    for chr_dir, telo_len_list in chr_fb_len_dict.items():
        s_telo_len_list = sorted(telo_len_list, key=lambda t: t[1])
        telo_len_data.extend(filter(lambda t: t[1] > 0, s_telo_len_list[1:]))

    need_label = defaultdict(list)
    for node_id, telo_len, chr_dir in telo_len_data:
        need_label[chr_dir[:-1]].append((node_id, chr_dir[-1]))

    new_telo_val_list = []
    for sector in circos.sectors:
        track1 = sector.add_track((58.2, 59.8))
        for j in need_label[sector.name]:
            cn = telo_cn[j[0]]
            norm_cn = np.clip((cn - lower) / (upper - lower), 0, 1)

            new_telo_val_list.append(cn / meandepth * 2)
            color = cmap(norm_cn)
            if j[1]=='f':
                # Prevent out of range
                try:
                    track1.arrow(contig_data[j[0]][CHR_STR], contig_data[j[0]][CHR_STR] + 15*M, fc=color, zorder=telo_zorder_dict[j[0]])
                except:
                    pass
            else:
                try:
                    track1.arrow(contig_data[j[0]][CHR_END], contig_data[j[0]][CHR_END] - 15*M, fc=color, zorder=telo_zorder_dict[j[0]])
                except:
                    pass

    with open(f'{PREFIX}/cn_data.pkl', 'wb') as f:
        pkl.dump((inv_val_list, transloc_val_list, indel_val_list, new_telo_val_list), f)

    fig = circos.plotfig(figsize=(10, 10), dpi=500)

    loc = 'upper left'

    breakend_line_handles = [
        Line2D([], [], color="black", label="Real CN", linestyle = '-'),
        Line2D([], [], color="green", label="SKYPE CN", linestyle = '-'),
    ]

    breakend_line_legend = circos.ax.legend(
        handles=breakend_line_handles,
        loc=loc,
        bbox_to_anchor=(0.85 + off, 1.03),
        fontsize=8,
        title="CN lines",
        handlelength=2,
    )
    circos.ax.add_artist(breakend_line_legend)

    scatter_handles = []
    for color, label in zip(target_color, ['Sucess', 'Failed', 'Ignored']):
        scatter_handles.append(mpl.patches.Patch(color=color, label=label))

    scatter_legend = circos.ax.legend(
        loc=loc,
        handles=scatter_handles,
        bbox_to_anchor=(0.97 + off, 1.03),
        fontsize=8,
        title="Predict label",
        handlelength=2,
    )
    circos.ax.add_artist(scatter_legend)

    cn_line_handles = [
        Line2D([], [], color="black", label="Breakend", linestyle = '-'),
        Line2D([], [], color="black", label="Inversion", linestyle = '-.'),
        Line2D([], [], color="black", label="Indel", linestyle = '--'),
        Line2D([], [], color="black", label="Amplicon", linestyle = (0, (5, 2))),
        Line2D([], [], color="black", label="Virtual inversion", linestyle = (0, (2, 2))),
    ]

    cn_line_legend = circos.ax.legend(
        loc=loc,
        handles=cn_line_handles,
        bbox_to_anchor=(0.85 + off, 0.94),
        fontsize=8,
        title="Breakend lines",
        handlelength=2,
    )
    circos.ax.add_artist(cn_line_legend)

    fig.savefig(f'{PREFIX}/total_cov{fig_prefix}.pdf')
    fig.savefig(f"{PREFIX}/total_cov{fig_prefix}.png")

weights = np.load(f'{PREFIX}/weight.npy')

use_julia_solver = pipeline_mode_is_karyotype(pipeline_mode_config)

draw_circos_plot()
if use_julia_solver:
    draw_circos_plot('_filter')
    draw_circos_plot('_cluster')

# Parse as vcf
def pairs_to_vcf(out_prefix=''):
    weights = np.load(f'{PREFIX}/weight{out_prefix}.npy')
    nclose_cn_std = defaultdict(float)
    for i, ctr in enumerate(path_nclose_usage):
        for j, v in ctr.items():
            nclose_cn_std[j] += v*weights[i]
    split_bnd_weights, split_parent_weight = build_ctg_intype_split_bnds(
        weights, VCF_FILTER_DEPTH_N * N
    )
    adjusted_nclose_cn_std = defaultdict(float, nclose_cn_std)
    for parent_nclose_idx, split_weight in split_parent_weight.items():
        adjusted_nclose_cn_std[parent_nclose_idx] = max(
            0.0, adjusted_nclose_cn_std[parent_nclose_idx] - split_weight
        )

    type1_nclose_node = set()
    type2_nclose_node = set()

    for nclose in nclose_set:
        s, e = nclose
        if contig_data[s][CHR_NAM] != contig_data[e][CHR_NAM]:
            type1_nclose_node.add(nclose)
        else:
            type2_nclose_node.add(nclose)

    display_indel = defaultdict(list)
    amplicon_events = build_amplicon_events(weights)
    virtual_inv_events = build_virtual_inv_events(
        PREFIX, meandepth, chr_len, min_depth_N=VCF_FILTER_DEPTH_N, weights=weights
    )

    # Gate must match make_bed_output() so SKYPE_result.bed and SV_call_result.vcf
    # report the SAME jump (DEL/DUP) call set.  This includes type4 indel graph
    # edges embedded inside ordinary paths by --add_indel_graph.
    for event in build_aggregated_indel_events(weights, BREAKEND_REMARKABLE_CN):
        chrom = event["chrom"]
        display_indel[chrom].append((
            event["event_type"],
            event["st"],
            event["nd"],
            event["weight"] / N,
            chrom,
            indel_event_source_label(event),
        ))
    
    all_nclose = []
    for nclose in nclose_set:
        if adjusted_nclose_cn_std[nclose2idx[nclose]] > VCF_FILTER_DEPTH_N * N:
            all_nclose.append(nclose)
                
    indel_event_count = sum(len(indel_list) for indel_list in display_indel.values())
    amplicon_event_count = len(amplicon_events)
    logging.info(
        f"{out_prefix[1:].capitalize()}{' ' if out_prefix else ''}"
        f"Total called breakends (DUP, DEL, BND, INV, AMP) : "
        f"{len(all_nclose) + len(split_bnd_weights) + indel_event_count + amplicon_event_count + len(virtual_inv_events)}"
    )

    vcf_path = f"{PREFIX}/SV_call_result{out_prefix}.vcf"
    _pairs_to_vcf(
        all_nclose, contig_data, chr_len, display_indel, amplicon_events, virtual_inv_events,
        split_bnd_weights, vcf_path, adjusted_nclose_cn_std
    )

if pipeline_mode_is_vcf_input(pipeline_mode_config):
    write_annotated_input_vcf(weights)
else:
    pairs_to_vcf()
    if use_julia_solver:
        pairs_to_vcf('_filter')
        pairs_to_vcf('_cluster')

# Bed output for further analysis

def make_bed_output(output_prefix=''):
    weights = np.load(f'{PREFIX}/weight{output_prefix}.npy')
    nclose_cn = defaultdict(float)

    for i, ctr in enumerate(path_nclose_usage):
        for j, v in ctr.items():
            nclose_cn[j] += v * weights[i]

    split_bnd_weights, split_parent_weight = build_ctg_intype_split_bnds(
        weights, BREAKEND_REMARKABLE_CN
    )
    adjusted_nclose_cn = defaultdict(float, nclose_cn)
    for parent_nclose_idx, split_weight in split_parent_weight.items():
        adjusted_nclose_cn[parent_nclose_idx] = max(
            0.0, adjusted_nclose_cn[parent_nclose_idx] - split_weight
        )

    def point_interval(pos):
        pos = int(pos)
        return max(0, pos - 1), pos
    
    with open(f'{PREFIX}/SKYPE_result{output_prefix}.bed', 'w') as f:
        cf = csv.writer(f, delimiter='\t')
        cf.writerow(['#chrom', 'cordst', 'cordnd', 'type', 'weight (N)'])

    
        for nclose_ind in nclose_cn:
            v = adjusted_nclose_cn[nclose_ind]
            if v > BREAKEND_REMARKABLE_CN:
                for ind in idx2nclose[nclose_ind]:
                    cf.writerow([contig_data[ind][CHR_NAM], contig_data[ind][CHR_STR], contig_data[ind][CHR_END],
                                'Breakend', round(v / N, 2)])

        for split_key, v in sorted(split_bnd_weights.items()):
            (
                _parent_nclose_idx,
                _split_idx,
                chr_a,
                pos_a,
                _dir_a,
                chr_b,
                pos_b,
                _dir_b,
                _ctg_name,
            ) = split_key
            if v > BREAKEND_REMARKABLE_CN:
                st, nd = point_interval(pos_a)
                cf.writerow([chr_a, st, nd, 'Breakend', round(v / N, 2)])
                st, nd = point_interval(pos_b)
                cf.writerow([chr_b, st, nd, 'Breakend', round(v / N, 2)])

        for event in build_aggregated_indel_events(weights, BREAKEND_REMARKABLE_CN):
            cf.writerow([
                event["chrom"], event["st"], event["nd"],
                'Deletion' if event["event_type"] == 'd' else 'Duplication',
                round(event["weight"] / N, 2),
            ])

        for i in range(rpll, min(len(weights), len(tot_loc_list))):
            v = weights[i]
            if v <= BREAKEND_REMARKABLE_CN:
                continue

            paf_loc = tot_loc_list[i]
            if paf_loc.split('/')[-3] == '12_cent_fragment':
                chrom = paf_loc.split('/')[-2]
                info = cen_fragment_meta[chrom]
                if info["dir"]:
                    st, nd = info["mid"], info["chr_len"]
                else:
                    st, nd = 0, info["mid"]
                cf.writerow([chrom, st, nd, 'Centromere', round(v / N, 2)])

        for chrom, st, nd, _, depth_N, _ in build_amplicon_events(weights):
            cf.writerow([chrom, st, nd, 'Amplicon', round(depth_N, 2)])

        for chrom, st, nd, _, depth_N, _name in build_virtual_inv_events(
            PREFIX, meandepth, chr_len, min_depth_N=BREAKEND_REMARKABLE_CN / N, weights=weights
        ):
            cf.writerow([chrom, st, nd, 'Virtual_inversion', round(depth_N, 2)])


if not pipeline_mode_is_vcf_input(pipeline_mode_config):
    make_bed_output()
    if use_julia_solver:
        make_bed_output('_cluster')
        make_bed_output('_filter')

# os.remove(f'{PREFIX}/matrix.h5')
logging.info("SKYPE pipeline end")
