import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sys
import logging
import argparse
import collections
import scipy.stats
import glob

from tqdm import tqdm
from scipy.optimize import nnls
from scipy.signal import butter, filtfilt
from concurrent.futures import ProcessPoolExecutor, as_completed

from pycirclize import Circos
from matplotlib.projections.polar import PolarAxes
from pycirclize.track import Track
from collections import defaultdict

logging.basicConfig(
    format='%(asctime)s %(levelname)s:%(message)s',
    level=logging.INFO,
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logging.info("22_depth_analysis start")

np.seterr(invalid='ignore')

ABS_MAX_COVERAGE_RATIO = 3
K = 1000
M = K * 1000
MAX_PATH_CNT = 100

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

def rebin_dataframe(df: pd.DataFrame, n: int = 2) -> pd.DataFrame:
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

def chr2int(x):
    chrXY2int = {'chrX' : 24, 'chrY' : 25}
    if x in chrXY2int:
        return chrXY2int[x]
    else:
        return int(x[3:])

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

def find_chr_len(file_path : str) -> dict:
    chr_data_file = open(file_path, "r")
    chr_len = {}
    for curr_data in chr_data_file:
        curr_data = curr_data.split("\t")
        chr_len[curr_data[0]] = int(curr_data[1])
    return chr_len

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


bpath = 'public_data/chm13v2.0_censat_v2.1.m.bed'


def import_bed(bed_path : str) -> dict:
    bed_data_file = open(bed_path, "r")
    chr_len = collections.defaultdict(list)
    for curr_data in bed_data_file:
        curr_data = curr_data.split("\t")
        chr_len[curr_data[0]].append((int(curr_data[1]), int(curr_data[2])))
    return chr_len

def inclusive_checker(tuple_a : tuple, tuple_b : tuple) -> bool :
    if int(tuple_a[0]) <= int(tuple_b[0]) and int(tuple_b[1]) <= int(tuple_a[1]):
        return True
    else:
        return False

bed_data = import_bed(bpath)

parser = argparse.ArgumentParser(description="SKYPE depth analysis")

parser.add_argument("main_stat_loc", 
                    help="Cancer coverage location file")

parser.add_argument("reference_fai_path", 
                    help="Path to the chromosome information file.")

parser.add_argument("reference_cytobands_path", 
                    help="Path to the cytoband information file.")


parser.add_argument("prefix", 
                    help="Pefix for pipeline")

parser.add_argument("-t", "--thread", 
                    help="Number of thread", type=int)

args = parser.parse_args()

# t = "python 22_depth_analysis.py /home/hyunwoo/51g_cancer_denovo/51_depth_data/KKU-100.win.stat.gz public_data/chm13v2.0.fa.fai public_data/chm13v2.0_cytobands_allchrs.bed 30_skype_pipe/KKU-100_05_36_56 -t 25"
# args = parser.parse_args(t.split()[2:])

PREFIX = args.prefix
THREAD = args.thread
CHROMOSOME_INFO_FILE_PATH = args.reference_fai_path
main_stat_loc = args.main_stat_loc

# HuH-28의 /home/hyunwoo/51g_cancer_denovo/51_depth_data/HuH-28.win.stat.gz를 읽어서 모든 경우 chr, st set을 가지고온다 (나머지 없는 경우를 0으로 채우게)
# chrM은 뺀다
RATIO_OUTLIER_FOLDER = f"{PREFIX}/11_ref_ratio_outliers/"
front_contig_path = RATIO_OUTLIER_FOLDER+"front_jump/"
back_contig_path = RATIO_OUTLIER_FOLDER+"back_jump/"


df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
df = df.query('chr != "chrM"')

meandepth = np.median(df['meandepth'])
chr_order_list = extract_groups(list(df['chr']))

chr_st_list = []
chr_flag_list = []
chr_filt_st_list = []

for l in df.itertuples(index=False):
    chr_st_list.append((l.chr, l.st))
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
    
    chr_flag_list.append(flag)

def get_vec_from_stat_loc(stat_loc_):
    df = pd.read_csv(stat_loc_, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
    df = df.query('chr != "chrM"')

    chr_st_data = dict()
    for l in df.itertuples(index=False):
         chr_st_data[(l.chr, l.st)] = l.meandepth
    
    fv = []
    for cs in chr_filt_st_list:
        if cs not in chr_st_data:
            fv.append(0)
        else:
            fv.append(chr_st_data[cs])
    
    v = []
    for cs in chr_st_list:
        if cs not in chr_st_data:
            v.append(0)
        else:
            v.append(chr_st_data[cs])
    
    return np.asarray(fv, dtype=np.float64), np.asarray(v, dtype=np.float64)

PATH_FILE_FOLDER = f"{PREFIX}/20_depth"
chr_chr_folder_path = sorted(glob.glob(PATH_FILE_FOLDER+"/*"))

def get_two_vec_list_folder(folder_path):
    filter_vec_list = []
    vec_list = []
    loc_list = []

    paf_paths = sorted(glob.glob(folder_path + "/*.win.stat.gz"))
    for stat_loc in paf_paths:
        fv, v = get_vec_from_stat_loc(stat_loc)

        filter_vec_list.append(fv)
        vec_list.append(v)
        loc_list.append(stat_loc)

    return filter_vec_list, vec_list, loc_list

main_filter_vec, main_vec = get_vec_from_stat_loc(main_stat_loc)

filter_vec_list = []
vec_list = []
tot_loc_list = []

with ProcessPoolExecutor(max_workers=THREAD) as executor:
    futures = [executor.submit(get_two_vec_list_folder, folder_path) for folder_path in chr_chr_folder_path]
    for future in tqdm(as_completed(futures), total=len(futures), desc='Parse coverage from gz files', disable=not sys.stdout.isatty()):
        fvl, vl, ll = future.result()

        filter_vec_list.extend(fvl)
        vec_list.extend(vl)
        tot_loc_list.extend(ll)

fclen = len(glob.glob(front_contig_path+"*"))
for i in tqdm(range(1, fclen//4 + 1), desc='Parse coverage from forward-directed outlier contig gz files', disable=not sys.stdout.isatty()):
    ov_loc = front_contig_path+f"{i}.win.stat.gz"
    bv_loc = front_contig_path+f"{i}_base.win.stat.gz"
    ofv, ov = get_vec_from_stat_loc(ov_loc)
    bfv, bv = get_vec_from_stat_loc(bv_loc)
    filter_vec_list.append(ofv-bfv)
    vec_list.append(ov-bv)
    tot_loc_list.append(ov_loc)

bclen = len(glob.glob(back_contig_path+"*"))
for i in tqdm(range(1, bclen//4 + 1), desc='Parse coverage from backward-directed outlier contig gz files', disable=not sys.stdout.isatty()):
    ov_loc = back_contig_path+f"{i}.win.stat.gz"
    bv_loc = back_contig_path+f"{i}_base.win.stat.gz"
    ofv, ov = get_vec_from_stat_loc(ov_loc)
    bfv, bv = get_vec_from_stat_loc(bv_loc)
    filter_vec_list.append(ofv+bfv)
    vec_list.append(ov+bv)
    tot_loc_list.append(ov_loc)

logging.info("Regression analysis is ongoing...")

A = np.vstack(filter_vec_list).T
B = main_filter_vec

weights, loss = nnls(A, B)

error = np.linalg.norm(A @ weights - B)
b_norm = np.linalg.norm(B)

logging.info(f'Error : {round(error, 4)}')
logging.info(f'Norm error : {round(error / b_norm, 4)}')

logging.info("Forming result images...")

# Calculate noise

grouped_data = defaultdict(lambda: {"positions": [], "values": []})
for i, (chrom, pos) in enumerate(chr_st_list):
    grouped_data[chrom]["positions"].append(pos)
    grouped_data[chrom]["values"].append(main_vec[i])

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

A = np.vstack(vec_list)
B = main_vec

smooth_B = rebin_dataframe_B(df, 10)
predicted_B = np.maximum(A.T.dot(weights), 0)
diff = predicted_B - smooth_B

threshold = amplitude * scipy.stats.norm.ppf(0.975)
miss_B = np.abs(diff) > threshold
over_B = ~np.asarray(chr_flag_list)

color_label = np.ones_like(B)
color_label[miss_B] = 3
color_label[over_B] = 5

logging.info(f'Fail ratio : {round(sum(color_label[miss_B] == 3) / len(B) * 100, 3)}%')

chr_len = find_chr_len(CHROMOSOME_INFO_FILE_PATH)
del chr_len['chrM']

rdf = rebin_dataframe(df, 2)

circos = Circos(chr_len, space=3)
circos.add_cytoband_tracks((95, 100), '../00_build_graph/public_data/chm13v2.0_cytobands_allchrs.bed')

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
    for i, v in enumerate(chr_st_list):
        if v[0] == sector.name:
            tchr_st_list.append(v[1])
            tchr_idx_list.append(i)
    
    cn_track = sector.add_track((50, 90))

    cn_track.axis()
    cn_track.grid(y_grid_num=7, zorder=-2)
    cn_track.line([0, max(tdf['st'])], [0.5 * meandepth] * 2, color='blue', vmax=ABS_MAX_COVERAGE_RATIO * meandepth, zorder=-1, alpha=0.5)
    cn_track.line([0, max(tdf['st'])], [2 * meandepth] * 2, color='red', vmax=ABS_MAX_COVERAGE_RATIO * meandepth, zorder=-1, alpha=0.5)
    line_track_plot(cn_track, x=tchr_st_list, y=predicted_B[tchr_idx_list], vmax=ABS_MAX_COVERAGE_RATIO * meandepth, color='green', linewidth=0.5, zorder=10)
    line_track_plot(cn_track, x=tdf['st'].to_list(), y=tdf['meandepth'].to_list(), vmax=ABS_MAX_COVERAGE_RATIO * meandepth, color='black', linewidth=0.2)
    line_track_circos_color(cn_track, x=tdf['st'].to_list(), y1=tdf['meandepth'].to_list(), vmax=ABS_MAX_COVERAGE_RATIO * meandepth,
                            target_color=['#6baed6', '#969696', '#fb6a4a'],
                            target_value=[0.5 * meandepth, 2 * meandepth])
    

    color_track = sector.add_track((40, 45))
    color_track.axis()
    line_track_circos_color(color_track, x=tchr_st_list, y1=color_label[tchr_idx_list], vmax=0.5,
                            target_color=['blue', 'red', 'gray'],
                            target_value=[2, 4])

fig = circos.plotfig(figsize=(10, 10), dpi=500)
fig.savefig(f"{PREFIX}/total_cov.png")
