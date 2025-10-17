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
import h5py
import logging
import argparse
import collections
import scipy.stats
import glob


from scipy.signal import butter, filtfilt
from pycirclize import Circos
from matplotlib.lines import Line2D
from matplotlib.projections.polar import PolarAxes
from pycirclize.track import Track
from collections import defaultdict
from collections import Counter

logging.basicConfig(
    format='%(asctime)s %(levelname)s:%(message)s',
    level=logging.INFO,
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logging.info("31_depth_analysis start")

BREAKEND_REMARKABLE_CN_RATIO = 0.05
TELOMERE_REMARKABLE_CN_RATIO = 0.05

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

ABS_MAX_COVERAGE_RATIO = 3
MAX_PATH_CNT = 100

DIR_FOR = 1
TELOMERE_EXPANSION = 5 * K

CONJOINED_CONTIG_MINIMUM_LENGTH = 200*K
SIMILAR_DEPTH_NCLOSE_RATIO = 0.2

TYPE2_FLANKING_LENGTH = 5*M
TYPE2_SIM_COMPARE_RAITO = 1.5
TYPE34_BREAK_CHUKJI_LIMIT = 1*M

NCLOSE_SIM_COMPARE_RAITO = 1.2
NCLOSE_SIM_DIFF_THRESHOLD = 5


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

def bnd_alt(t, mate_chr, mate_pos, form):
    """
    Build ALT per the 4 canonical breakend forms:
        "t[p[", "t]p]", "]p]t", "[p[t"
    """
    p = f"{mate_chr}:{mate_pos}"
    if form == "t[p[":
        return f"{t}[{p}["
    if form == "t]p]":
        return f"]{p}]{t}"
    if form == "]p]t":
        return f"]{p}]{t}"
    if form == "[p[t":
        return f"[{p}[{t}"
    assert(False)

def choose_alt_forms(dir_a, dir_b):
    """
    Map (dir_a, dir_b) to a complementary ALT pair.
    You can tweak this mapping to your exact strand semantics.
    """
    if dir_a == '+' and dir_b == '+':
        return ("t[p[", "]p]t")
    if dir_a == '-' and dir_b == '-':
        return ("t]p]", "[p[t")
    if dir_a == '+' and dir_b == '-':
        return ("t[p[", "[p[t")
    if dir_a == '-' and dir_b == '+':
        return ("t]p]", "]p]t")
    assert(False)

def make_strands(dir_a, dir_b):
    """
    Compose STRANDS string as '<A><B>' where each is '+' or '-'.
    Fallback to '.' if a direction is unknown.
    """
    a = dir_a if dir_a in ('+', '-') else '.'
    b = dir_b if dir_b in ('+', '-') else '.'
    return f"{a}{b}"

def write_vcf_header(fh, contig_lengths):
    fh.write("##fileformat=VCFv4.3\n")
    fh.write('##ALT=<ID=BND,Description="Breakend">\n')
    fh.write('##ALT=<ID=INV,Description="Inversion">\n')
    fh.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
    fh.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of SV">\n')
    fh.write('##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of mate breakend">\n')
    for chrom, length in contig_lengths.items():
        fh.write(f"##contig=<ID={chrom},length={length}>\n")
    fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

def pairs_to_vcf(pairs, contig_data, contig_lengths, display_indel, out_vcf_path):
    with open(out_vcf_path, 'w') as fo:
        # 1) 헤더 쓰기
        write_vcf_header(fo, contig_lengths)

        # 2) Translocation, Inversion 처리
        tr_counter = 1
        for a_idx, b_idx in pairs:
            a = contig_data[a_idx]
            b = contig_data[b_idx]

            chr_a, chr_b = a[CHR_NAM], b[CHR_NAM]
            pos_a = a[CHR_END]
            pos_b = b[CHR_STR]

            dir_a = a[CTG_DIR]
            dir_b = b[CTG_DIR]

            strands = make_strands(dir_a, dir_b)

            sv_id_a = f"SKYPE.BND.{tr_counter}_1"
            sv_id_b = f"SKYPE.BND.{tr_counter}_2"
            tr_counter += 1

            form_a, form_b = choose_alt_forms(dir_a, dir_b)

            ref = "N"
            alt_a = bnd_alt("N", chr_b, pos_b, form_a)
            alt_b = bnd_alt("N", chr_a, pos_a, form_b)

            fo.write(
                f"{chr_a}\t{pos_a}\t{sv_id_a}\t{ref}\t{alt_a}\t60\tPASS\t"
                f"SVTYPE=BND;MATEID={sv_id_b};STRANDS={strands}\n"
            )
            fo.write(
                f"{chr_b}\t{pos_b}\t{sv_id_b}\t{ref}\t{alt_b}\t60\tPASS\t"
                f"SVTYPE=BND;MATEID={sv_id_a};STRANDS={strands}\n"
            )

        # 3) indel 처리
        ins_counter = 1
        del_dounter = 1
        
        for chrom, indel_list in display_indel.items():
            for indel in indel_list:
                indel_type, start, end, _, _ = indel
                if indel_type == 'd':
                    sv_id = f"SKYPE.DEL.{del_dounter}"
                    svlen = -(end - start)
                    alt = "<DEL>"

                    del_dounter += 1
                elif indel_type == 'i':
                    sv_id = f"SKYPE.INS.{ins_counter}"
                    svlen = end - start
                    alt = "<INS>"

                    ins_counter += 1
                else:
                    assert(False)

                fo.write(
                    f"{chrom}\t{start}\t{sv_id}\tN\t{alt}\t60\tPASS\t"
                    f"SVTYPE={alt[1:-1]};END={end};SVLEN={svlen}\n"
                )

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

RATIO_OUTLIER_FOLDER = f"{PREFIX}/11_ref_ratio_outliers/"
front_contig_path = RATIO_OUTLIER_FOLDER+"front_jump/"
back_contig_path = RATIO_OUTLIER_FOLDER+"back_jump/"
ecdna_contig_path = RATIO_OUTLIER_FOLDER + "ecdna/"

TELO_CONNECT_NODES_INFO_PATH = PREFIX+"/telomere_connected_list.txt"
vcf_path = f"{PREFIX}/SV_call_result.vcf"

contig_data = import_data2(PREPROCESSED_PAF_FILE_PATH)
telo_connected_node = extract_telomere_connect_contig(TELO_CONNECT_NODES_INFO_PATH)

with open(f'{PREFIX}/path_data.pkl', 'rb') as f:
    path_list_dict = pkl.load(f)

with h5py.File(f"{PREFIX}/matrix.h5", "r") as f:
    dB = f["B"]
    B = np.empty(dB.shape, dtype=dB.dtype)
    dB.read_direct(B)

    dBf = f["B_fail"]
    B_fail = np.empty(dBf.shape, dtype=dBf.dtype)
    dBf.read_direct(B_fail)

    b_start_ind = int(f["B_depth_start"][()])

B = np.hstack([B, B_fail])[b_start_ind:]

df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
df = df.query('chr != "chrM"')

meandepth = np.median(df['meandepth'])

BREAKEND_REMARKABLE_CN = meandepth * BREAKEND_REMARKABLE_CN_RATIO
TELOMERE_REMARKABLE_CN = meandepth * TELOMERE_REMARKABLE_CN_RATIO
N = meandepth / 2

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


with open(f'{PREFIX}/tot_loc_list.pkl', 'rb') as f:
    tot_loc_list = pkl.load(f)

grouped_data = defaultdict(lambda: {"positions": [], "values": []})
for i, (chrom, pos) in enumerate(chr_filt_st_list + chr_no_filt_st_list):
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

def extract_nclose_node(nclose_path: str) -> list:
        nclose_list = []
        with open(nclose_path, "r") as f:
            for line in f:
                line = line.split()
                nclose_list.append((int(line[1]), int(line[2])))
        return nclose_list

NCLOSE_FILE_PATH = f"{args.prefix}/nclose_nodes_index.txt"

nclose_list = extract_nclose_node(NCLOSE_FILE_PATH)
nclose_set = set(nclose_list)
nclose_dict = dict(zip(nclose_list, range(len(nclose_list))))
reverse_nclose_dict = dict(zip(range(len(nclose_list)), nclose_list))

nclose_str_pos = dict()
for k, v in reverse_nclose_dict.items():
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
            using_nclose[nclose_dict[nclose_cand]]+=1
            s+=2
        else:
            s+=1
    if path[1][1] in telo_node_set:
        using_telo[path[1][1]]+=1
    if path[len(path)-2][1] in telo_node_set:
        using_telo[path[len(path)-2][1]]+=1

    path_nclose_usage.append(using_nclose)
    path_telo_usage.append(using_telo)          


def draw_circos_plot(weights, fig_prefix=''):
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
        msg_prefix = 'Cluster fail'

    logging.info(f'{msg_prefix} ratio : {round(sum(color_label[miss_B] == 3) / len(B) * 100, 3)}%')

    nclose_cn = defaultdict(float)
    for i, ctr in enumerate(path_nclose_usage):
        for j, v in ctr.items():
            nclose_cn[j] += v*weights[i]
    
    # with open(f"{PREFIX}/nclose_cn.txt", "wt") as f:
    #     for k, v in nclose_cn.items():
    #         if v > BREAKEND_REMARKABLE_CN:
    #             curr_nclose = reverse_nclose_dict[k]
    #             print(contig_data[curr_nclose[0]], "\n", contig_data[curr_nclose[1]], "\n", v, "\n", file=f)


    rdf = rebin_dataframe(df, 2)

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

    for k, v in nclose_cn.items():
        pos1, pos2 = nclose_str_pos[k]
        idx1, idx2 = reverse_nclose_dict[k]
        chr_nam1 = contig_data[idx1][CHR_NAM]
        chr_nam2 = contig_data[idx2][CHR_NAM]

        if chr_nam1 != chr_nam2 and v > BREAKEND_REMARKABLE_CN:
            transloc_val_list.append(v / meandepth * 2)
        elif chr_nam1 == chr_nam2 and contig_data[idx1][CTG_DIR] != contig_data[idx2][CTG_DIR] and v > BREAKEND_REMARKABLE_CN:
            inv_val_list.append(v / meandepth * 2)

        virtual_flag = False
        if contig_data[idx1][CTG_NAM].startswith('v'):
            virtual_flag = True
        
        if v > BREAKEND_REMARKABLE_CN:
            bnd_cn_data.append([(chr_nam1, pos1), (chr_nam2, pos2), v, virtual_flag])

    display_indel = defaultdict(list)

    indel_val_list = []
    non_type4_cnt = len(bnd_cn_data)
    rpll = len(raw_path_list)
    for i in range(rpll, len(weights)):
        paf_loc = tot_loc_list[i]
        indel_ind = paf_loc.split('/')[-2]
        with open(tot_loc_list[i], "r") as f:
            l = f.readline()
            l = l.rstrip()
            l = l.split("\t")
            chr_nam1 = l[CHR_NAM]
            chr_nam2 = l[CHR_NAM]
            pos1 = int(l[CHR_STR])
            pos2 = int(l[CHR_END])
        v = weights[i]
        chrom = chr_nam1
        if abs(pos1-pos2) > TYPE2_FLANKING_LENGTH and v > NCLOSE_SIM_DIFF_THRESHOLD:
            if exist_near_bnd(chrom, pos1, pos1) or \
            exist_near_bnd(chrom, pos2, pos2):
                if indel_ind == 'front_jump':
                    display_indel[chrom].append(("d", pos1, pos2, v/meandepth * 2, chrom))
                if indel_ind == 'back_jump':
                    display_indel[chrom].append(("i", pos1, pos2, v/meandepth * 2, chrom))
        indel_val_list.append(v / meandepth * 2)
        bnd_cn_data.append([(chr_nam1, pos1), (chr_nam2, pos2), v, False])

    a = sorted(list(telo_cn.items()), key = lambda t:t[1])
    telo_zorder_dict = {}
    for i, v in enumerate(a):
        telo_zorder_dict[v[0]] = i

    off = 0.03
    lower = 0
    upper = meandepth * 2

    cmap = sns.color_palette("rocket_r", as_cmap=True)

    # bnd_cn_data = [[('chr1', 60 * M), ('chr2', 30 * M), 4.6, False],
    #                [('chr13', 10 * M), ('chr21', 25 * M), 21.6, False],
    #                [('chr1', 125 * M), ('chr9', 60 * M), 13.6, True]]
    bnd_cn_data = sorted(bnd_cn_data, key=lambda t: t[2])

    for i, (bnd_loc1, bnd_loc2, cn, is_vir) in enumerate(bnd_cn_data):
        norm_cn = np.clip((cn - lower) / (upper - lower), 0, 1)

        color = cmap(norm_cn)
        linestyle = '--' if is_vir else '-'
        
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
        Line2D([], [], color="black", label="Virtual breakend", linestyle = '--'),
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
draw_circos_plot(weights)

with open(f"{PREFIX}/report.txt", 'r') as f:
    f.readline()
    path_cnt = int(f.readline().strip())

use_julia_solver = path_cnt <= HARD_PATH_COUNT_BASELINE

if use_julia_solver:
    weights_cluster = np.load(f'{PREFIX}/weight_cluster.npy')
    draw_circos_plot(weights_cluster, '_cluster')

# Parse as vcf

nclose_cn_std = defaultdict(float)
for i, ctr in enumerate(path_nclose_usage):
    for j, v in ctr.items():
        nclose_cn_std[j] += v*weights[i]

conjoined_nclose_node_set = set()

type1_nclose_node = set()
type2_nclose_node = set()

for nclose in nclose_set:
    s, e = nclose
    if contig_data[s][CHR_NAM] != contig_data[e][CHR_NAM]:
        type1_nclose_node.add(nclose)
    else:
        type2_nclose_node.add(nclose)

for i, nclose_a in enumerate(nclose_list):
    for nclose_b in nclose_list[i+1:]: 
        if abs(1 - nclose_cn_std[nclose_dict[nclose_a]]/(nclose_cn_std[nclose_dict[nclose_b]]+1e-9)) < SIMILAR_DEPTH_NCLOSE_RATIO: # Conjoin only if depth is similar # Todo 03 Similar로 바꾸기
            flag = True # To check if nclose_a and nclose_b can be conjoined
            connecting_nclosea_node = -1
            chrom = None
            for idx in (0, 1):
                curr_chrom = contig_data[nclose_a[idx]][CHR_NAM]
                if contig_data[nclose_b[0]][CHR_NAM] == curr_chrom or contig_data[nclose_b[1]][CHR_NAM] == curr_chrom:
                    flag = False
                    if connecting_nclosea_node < 0:
                        connecting_nclosea_node = idx
                        chrom = curr_chrom  # chrom : sharing chromosome type
                    else:
                        flag = True # Pass if same chromosome component
            if flag:
                continue
            # Only two cases:
            # 1. nclose_a is type2, nclose_b is type1
            # 2. nclose_a is type1, nclose_b is type1
            nclosea_contig_front = contig_data[nclose_a[0]]
            nclosea_contig_back = contig_data[nclose_a[1]]
            front_is_connectable = nclosea_contig_front[CHR_NAM] == chrom
            back_is_connectable = nclosea_contig_back[CHR_NAM] == chrom
            if contig_data[nclose_b[0]][CHR_NAM] == chrom:
                template_contig = contig_data[nclose_b[0]]
                if template_contig[CTG_DIR] == '+':
                    if nclosea_contig_back[CTG_DIR] == '+' and back_is_connectable:
                        dist = distance_checker(nclosea_contig_back, template_contig)
                        if (nclosea_contig_back[CHR_END] < template_contig[CHR_STR] and dist < CONJOINED_CONTIG_MINIMUM_LENGTH) or dist == 0:
                            conjoined_nclose_node_set.add((nclose_b[1], nclose_a[0]))
                    elif nclosea_contig_front[CTG_DIR] == '-' and front_is_connectable:
                        dist = distance_checker(nclosea_contig_front, template_contig)
                        if (nclosea_contig_front[CHR_END] < template_contig[CHR_STR] and dist < CONJOINED_CONTIG_MINIMUM_LENGTH) or dist == 0:
                            conjoined_nclose_node_set.add((nclose_b[1], nclose_a[1]))
                else:
                    if nclosea_contig_back[CTG_DIR] == '-' and back_is_connectable:
                        dist = distance_checker(nclosea_contig_back, template_contig)
                        if (nclosea_contig_back[CHR_STR] > template_contig[CHR_END] and dist < CONJOINED_CONTIG_MINIMUM_LENGTH) or dist == 0:
                            conjoined_nclose_node_set.add((nclose_b[1], nclose_a[0]))
                    elif nclosea_contig_front[CTG_DIR] == '+' and front_is_connectable:
                        dist = distance_checker(nclosea_contig_front, template_contig)
                        if (nclosea_contig_front[CHR_STR] > template_contig[CHR_END] and dist < CONJOINED_CONTIG_MINIMUM_LENGTH) or dist == 0:
                            conjoined_nclose_node_set.add((nclose_b[1], nclose_a[1]))
            elif contig_data[nclose_b[1]][CHR_NAM] == chrom:
                template_contig = contig_data[nclose_b[1]]  
                if template_contig[CTG_DIR] == '+':
                    if nclosea_contig_back[CTG_DIR] == '-' and back_is_connectable:
                        dist = distance_checker(nclosea_contig_back, template_contig)
                        if (nclosea_contig_back[CHR_STR] > template_contig[CHR_END] and dist < CONJOINED_CONTIG_MINIMUM_LENGTH) or dist == 0:
                            conjoined_nclose_node_set.add((nclose_b[0], nclose_a[0]))
                    elif nclosea_contig_front[CTG_DIR] == '+' and front_is_connectable:
                        dist = distance_checker(nclosea_contig_front, template_contig)
                        if (nclosea_contig_front[CHR_STR] > template_contig[CHR_END] and dist < CONJOINED_CONTIG_MINIMUM_LENGTH) or dist == 0:
                            conjoined_nclose_node_set.add((nclose_b[0], nclose_a[1]))
                else:
                    if nclosea_contig_back[CTG_DIR] == '+' and back_is_connectable:
                        dist = distance_checker(nclosea_contig_back, template_contig)
                        if (nclosea_contig_back[CHR_END] < template_contig[CHR_STR] and dist < CONJOINED_CONTIG_MINIMUM_LENGTH) or dist == 0:
                            conjoined_nclose_node_set.add((nclose_b[0], nclose_a[0]))
                    elif nclosea_contig_front[CTG_DIR] == '-' and front_is_connectable:
                        dist = distance_checker(nclosea_contig_front, template_contig)
                        if (nclosea_contig_front[CHR_END] < template_contig[CHR_STR] and dist < CONJOINED_CONTIG_MINIMUM_LENGTH) or dist == 0:
                            conjoined_nclose_node_set.add((nclose_b[0], nclose_a[1]))

significant_nclose = []
for nclose in (nclose_set | conjoined_nclose_node_set):
    st, ed = nclose
    if exist_near_bnd(contig_data[st][CHR_NAM], contig_data[st][CHR_STR], contig_data[st][CHR_END]) or \
       exist_near_bnd(contig_data[ed][CHR_NAM], contig_data[ed][CHR_STR], contig_data[ed][CHR_END]):
        significant_nclose.append(nclose)

significant_nclose = set(significant_nclose)

display_indel = defaultdict(list)

inv_val_list = []
transloc_val_list = []
bnd_cn_data = []

for k, v in nclose_cn_std.items():
    pos1, pos2 = nclose_str_pos[k]
    idx1, idx2 = reverse_nclose_dict[k]
    chr_nam1 = contig_data[idx1][CHR_NAM]
    chr_nam2 = contig_data[idx2][CHR_NAM]

    if chr_nam1 != chr_nam2 and v > BREAKEND_REMARKABLE_CN:
        transloc_val_list.append(v / meandepth * 2)
    elif chr_nam1 == chr_nam2 and contig_data[idx1][CTG_DIR] != contig_data[idx2][CTG_DIR] and v > BREAKEND_REMARKABLE_CN:
        inv_val_list.append(v / meandepth * 2)

    virtual_flag = False
    if contig_data[idx1][CTG_NAM].startswith('v'):
        virtual_flag = True
    
    if v > BREAKEND_REMARKABLE_CN:
        bnd_cn_data.append([(chr_nam1, pos1), (chr_nam2, pos2), v, virtual_flag])

indel_val_list = []
non_type4_cnt = len(bnd_cn_data)

rpll = len(raw_path_list)
for i in range(rpll, len(weights)):
    paf_loc = tot_loc_list[i]
    indel_ind = paf_loc.split('/')[-2]
    with open(tot_loc_list[i], "r") as f:
        l = f.readline()
        l = l.rstrip()
        l = l.split("\t")
        chr_nam1 = l[CHR_NAM]
        chr_nam2 = l[CHR_NAM]
        pos1 = int(l[CHR_STR])
        pos2 = int(l[CHR_END])
    v = weights[i]
    chrom = chr_nam1
    if abs(pos1-pos2) > TYPE2_FLANKING_LENGTH and v > NCLOSE_SIM_DIFF_THRESHOLD:
        if exist_near_bnd(chrom, pos1, pos1) or \
        exist_near_bnd(chrom, pos2, pos2):
            if indel_ind == 'front_jump':
                display_indel[chrom].append(("d", pos1, pos2, v/meandepth * 2, chrom))
            if indel_ind == 'back_jump':
                display_indel[chrom].append(("i", pos1, pos2, v/meandepth * 2, chrom))
    indel_val_list.append(v / meandepth * 2)
    bnd_cn_data.append([(chr_nam1, pos1), (chr_nam2, pos2), v, False])

logging.info(f"Total called breakends (INS, DEL, BND, INV) : {len(significant_nclose) + len(display_indel)}")

pairs_to_vcf(significant_nclose, contig_data, chr_len, display_indel, vcf_path)

# Bed output for further analysis

fclen = len(glob.glob(front_contig_path + "*"))
bclen = len(glob.glob(back_contig_path + "*"))
eclen = len(glob.glob(ecdna_contig_path + "*"))

n = len(paf_ans_list) + fclen // 4 + bclen // 4 + eclen // 2

with open(f'{PREFIX}/ecdna_circuit_data.pkl', 'rb') as f:
    ecdna_circuit, _ = pkl.load(file=f)

def make_bed_output(weights, output_prefix=''):
    nclose_cn = defaultdict(float)

    for i, ctr in enumerate(path_nclose_usage):
        for j, v in ctr.items():
            nclose_cn[j] += v * weights[i]
    
    with open(f'{PREFIX}/SKYPE_result{output_prefix}.bed', 'w') as f:
        cf = csv.writer(f, delimiter='\t')
        cf.writerow(['#chrom', 'cordst', 'cordnd', 'type', 'weight (N)'])

    
        for nclose_ind, v in nclose_cn.items():
            if v > BREAKEND_REMARKABLE_CN:
                for ind in reverse_nclose_dict[nclose_ind]:
                    cf.writerow([contig_data[ind][CHR_NAM], contig_data[ind][CHR_STR], contig_data[ind][CHR_END],
                                'Breakend', round(v / N, 2)])

        for i in range(rpll, len(weights)):
            v = weights[i]

            if v > BREAKEND_REMARKABLE_CN:
                paf_loc = tot_loc_list[i]
                event_type = paf_loc.split('/')[-2]

                if event_type in {'front_jump', 'back_jump'}:
                    with open(tot_loc_list[i], "r") as f:
                        l = f.readline()
                        l = l.rstrip()
                        l = l.split("\t")
                        

                        cf.writerow([l[CHR_NAM], l[CHR_STR], l[CHR_END],
                                        'Insertion' if event_type == 'front_jump' else 'Deletion',
                                        round(v / N, 2)])

                elif event_type == 'ecdna':
                    ecdna_ind = i - (rpll + fclen // 4 + bclen // 4)

                    ctg_st_list = []
                    ctg_nd_list = []

                    for ctg_ind in ecdna_circuit[ecdna_ind]:
                        ctg_st_list.append(contig_data[ctg_ind][CHR_STR])
                        ctg_nd_list.append(contig_data[ctg_ind][CHR_END])

                    ctg_ind = ecdna_circuit[ecdna_ind][0]
                    
                    ctg_st = min(ctg_st_list)
                    ctg_nd = max(ctg_nd_list)

                    cf.writerow([contig_data[ctg_ind][CHR_NAM], ctg_st, ctg_nd,
                                     'Amplicon', round(v / N * 2, 2)])

                else:
                    assert(False)


make_bed_output(weights)

if use_julia_solver:
    make_bed_output(weights_cluster, '_cluster')

os.remove(f'{PREFIX}/matrix.h5')
logging.info("SKYPE pipeline end")
