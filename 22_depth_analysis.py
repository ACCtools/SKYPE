import numpy as np
import pandas as pd
import sys
import collections
import matplotlib.pyplot as plt
import glob

from tqdm import tqdm
from scipy.optimize import nnls

ABS_MAX_COVERAGE_RATIO = int(sys.argv[3])

def chr2int(x):
    chrXY2int = {'chrX' : 24, 'chrY' : 25}
    if x in chrXY2int:
        return chrXY2int[x]
    else:
        return int(x[3:])

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


# HuH-28의 /home/hyunwoo/51g_cancer_denovo/51_depth_data/HuH-28.win.stat.gz를 읽어서 모든 경우 chr, st set을 가지고온다 (나머지 없는 경우를 0으로 채우게)
# chrM은 뺀다
RATIO_OUTLIER_FOLDER = f"{sys.argv[2]}/11_ref_ratio_outliers/"
main_stat_loc = sys.argv[1]
front_contig_path = RATIO_OUTLIER_FOLDER+"front_jump/"
back_contig_path = RATIO_OUTLIER_FOLDER+"back_jump/"


df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
df = df.query('chr != "chrM"')
med_meandepth = np.median(df['meandepth'])

chr_st_list = []
chr_filt_st_list = []

for l in df.itertuples(index=False):
    chr_st_list.append((l.chr, l.st))
    flag = True
    if l.meandepth > ABS_MAX_COVERAGE_RATIO * med_meandepth:
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

PATH_FILE_FOLDER = f"{sys.argv[2]}/20_depth"
chr_chr_folder_path = glob.glob(PATH_FILE_FOLDER+"/*")

main_filter_vec, main_vec = get_vec_from_stat_loc(main_stat_loc)

filter_vec_list = []
vec_list = []
for folder_path in tqdm(chr_chr_folder_path, desc='Parse coverage from gz files', disable=not sys.stdout.isatty()):
    paf_paths = glob.glob(folder_path + "/*.win.stat.gz")
    for stat_loc in paf_paths:
        fv, v = get_vec_from_stat_loc(stat_loc)

        filter_vec_list.append(fv)
        vec_list.append(v)

fclen = len(glob.glob(front_contig_path+"*"))
for i in tqdm(range(1, fclen//4 + 1), desc='Parse coverage from forward-directed outlier contig gz files', disable=not sys.stdout.isatty()):
    ov_loc = front_contig_path+f"{i}.win.stat.gz"
    bv_loc = front_contig_path+f"{i}_base.win.stat.gz"
    ofv, ov = get_vec_from_stat_loc(ov_loc)
    bfv, bv = get_vec_from_stat_loc(bv_loc)
    filter_vec_list.append(ofv-bfv)
    vec_list.append(ov-bv)

bclen = len(glob.glob(back_contig_path+"*"))
for i in tqdm(range(1, bclen//4 + 1), desc='Parse coverage from backward-directed outlier contig gz files', disable=not sys.stdout.isatty()):
    ov_loc = back_contig_path+f"{i}.win.stat.gz"
    bv_loc = back_contig_path+f"{i}_base.win.stat.gz"
    ofv, ov = get_vec_from_stat_loc(ov_loc)
    bfv, bv = get_vec_from_stat_loc(bv_loc)
    filter_vec_list.append(ofv+bfv)
    vec_list.append(ov+bv)


print("Regression analysis is ongoing...")

A = np.vstack(filter_vec_list).T
B = main_filter_vec

weights, loss = nnls(A, B)

error = np.linalg.norm(A @ weights - B)
b_norm = np.linalg.norm(B)

print(f'Error : {round(error, 4)}')
print(f'Norm error : {round(error / b_norm, 4)}')

print("Forming result images...")

A = np.vstack(vec_list)
B = main_vec

# NNLS 예측 coverage: 각 구성요소의 기여도를 합산하여 하나의 벡터로 만듦
predicted_B = A.T.dot(weights)

# 4. 각 window의 각도 계산 (circos-like 배치를 위해)
total_windows = len(chr_st_list)

# 각 window의 (chr, st)를 그룹화 (염색체별 인덱스 목록)
chr_to_indices = collections.defaultdict(list)
for i, (chr_name, st) in enumerate(chr_st_list):
    chr_to_indices[chr_name].append(i)

# 염색체 이름들을 정렬 (원하는 순서로 조정 가능)
chromosomes = sorted(chr_to_indices.keys(), key=lambda x:chr2int(x))

# 전체 원(2π)을 각 염색체의 window 비율에 따라 나누기
chr_angles = {}
current_angle = 0
for chr_name in chromosomes:
    n_wins = len(chr_to_indices[chr_name])
    angle_span = 2 * np.pi * n_wins / total_windows
    chr_angles[chr_name] = (current_angle, current_angle + angle_span)
    current_angle += angle_span

# 각 window에 대해 각도를 할당 (해당 염색체 내에서 균등 분포)
window_angles = np.zeros(total_windows)
for chr_name in chromosomes:
    start_angle, end_angle = chr_angles[chr_name]
    indices = chr_to_indices[chr_name]
    n = len(indices)
    angles = np.linspace(start_angle, end_angle, n, endpoint=False)
    for idx, angle in zip(indices, angles):
        window_angles[idx] = angle

# 5. Polar plot (Circos-like plot) 생성: 원래의 coverage와 NNLS 예측 coverage 비교
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(15, 15), dpi=500)
ax.set_theta_direction(-1)       # 시계 반대 방향
ax.set_theta_zero_location('N')  # 0도가 북쪽에 위치
ax.set_facecolor('white')

# 원래의 coverage(B)와 예측 coverage를 선으로 플롯합니다.
ax.scatter(window_angles, B, color='gray', linewidth=1.5, label='Observed Coverage (B)', s=2, alpha=0.5,)
ax.scatter(window_angles, predicted_B, color='red', linewidth=1.5, label='Predicted Coverage (NNLS)', s=2, alpha=0.5,)

# 염색체 경계 표시 및 라벨 추가
for chr_name in chromosomes:
    start_angle, end_angle = chr_angles[chr_name]
    # 염색체 경계선 (시작 각도에 선을 그림)
    ax.axvline(x=start_angle, color='black', linewidth=0.8, alpha=0.7)
    # 염색체 이름 라벨을 중앙 각도에 배치 (최대 반지름보다 약간 바깥쪽)
    mid_angle = (start_angle + end_angle) / 2
    max_r = max(B.max(), predicted_B.max())
    ax.text(mid_angle, np.median(B) * 2.5, chr_name, ha='center', va='center', fontsize=10, fontweight='bold')

ax.set_yticklabels([])
ax.set_xticklabels([])

ax.set_ylim(ax.get_ylim()[0], np.median(B) * 3)

plt.title("Circos-like Plot: Observed vs. NNLS Predicted Coverage", fontsize=14)
plt.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize='small')
plt.tight_layout()
plt.savefig(f"{sys.argv[2]}/total_cov.png")

# print("NNLS weights:", weights)

# ----------------------------
# 상위 50개 구성요소 선택 및 각 구성요소의 coverage 계산
# ----------------------------

# weight가 높은 순으로 상위 50개 구성요소의 인덱스 (가상 염색체)
top_indices = np.argsort(weights)[-50:]
top_indices = top_indices[::-1]  # 내림차순 정렬

# 각 구성요소의 contribution profile: weight * 해당 A 행
top_components_coverage = {}
for idx in top_indices:
    # A[idx, :]는 각 window에 대한 값, 여기에 weight를 곱함
    top_components_coverage[idx] = weights[idx] * A[idx, :]

# ----------------------------
# chr_st_list를 이용해 각 윈도우의 위치(각도)를 계산
# ----------------------------

# 각 윈도우에 해당하는 (chr, st) 정보를 그룹화
chr_to_indices = collections.defaultdict(list)
for i, (chr_name, st) in enumerate(chr_st_list):
    chr_to_indices[chr_name].append(i)

# 염색체 목록 (알파벳 순 혹은 원하는 순서로 정렬)
chromosomes = sorted(chr_to_indices.keys(), key=lambda x:chr2int(x))

# 전체 윈도우 수
total_windows = len(chr_st_list)
# 각 염색체에 할당할 원형 아크의 각도는 그 염색체의 윈도우 수 비율에 따라 결정
chr_angles = {}
current_angle = 0
for chr_name in chromosomes:
    n_wins = len(chr_to_indices[chr_name])
    angle_span = 2 * np.pi * n_wins / total_windows
    chr_angles[chr_name] = (current_angle, current_angle + angle_span)
    current_angle += angle_span

# 각 윈도우별 각도 계산 (각 염색체 내에서는 균등 분포)
window_angles = np.zeros(total_windows)
for chr_name in chromosomes:
    start_angle, end_angle = chr_angles[chr_name]
    indices = chr_to_indices[chr_name]
    n = len(indices)
    angles = np.linspace(start_angle, end_angle, n, endpoint=False)
    for idx, angle in zip(indices, angles):
        window_angles[idx] = angle

# ----------------------------
# Circos-like plot 생성 (matplotlib의 polar plot 사용)
# ----------------------------

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'}, figsize=(15, 15), dpi=500)
ax.set_theta_direction(-1)
ax.set_theta_zero_location('N')
ax.set_facecolor('white')

# 기본 coverage(B)를 배경으로 scatter plot (회색, 투명도 적용)
ax.scatter(window_angles, B, s=5, color='gray', alpha=0.5, label='Coverage B')

# 상위 50개 구성요소의 coverage profile을 선으로 그리기 (각 구성요소별 다른 색상)
cmap = plt.get_cmap('tab20')
for j, idx in enumerate(top_indices):
    comp_cov = top_components_coverage[idx]
    # comp_cov: 각 window에 대한 값; window_angles: 각 window의 각도
    ax.plot(window_angles, comp_cov, linewidth=1.5, 
            alpha=0.8, color=cmap(j % 20), label=f'Comp {idx}')

# 염색체 경계 그리기 및 라벨 추가
for chr_name in chromosomes:
    start_angle, end_angle = chr_angles[chr_name]
    # 경계선
    ax.axvline(x=start_angle, color='black', linewidth=0.8, alpha=0.7)
    # 라벨: 각 염색체 아크의 중간 위치에 표시
    mid_angle = (start_angle + end_angle) / 2
    # 반지름은 최대 B 값의 약간 위로 배치
    ax.text(mid_angle, np.median(B) * 2.5, chr_name, ha='center', va='center', fontsize=10, fontweight='bold')

# 옵션: 불필요한 그리드/눈금 제거
ax.set_yticklabels([])
ax.set_xticklabels([])

ax.set_ylim(ax.get_ylim()[0], np.median(B) * 3) 

plt.title("Circos-like Plot: Coverage and Top 50 Components")
plt.legend(loc='upper right', bbox_to_anchor=(1.3, 1.1), fontsize='small')
plt.tight_layout()
plt.savefig(f'{sys.argv[2]}/virtual_chromosome_circosplot.png')