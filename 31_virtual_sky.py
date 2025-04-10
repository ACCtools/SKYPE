import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle as pkl
import matplotlib.patches as patches

import ast
import logging
import argparse
from datetime import datetime

from matplotlib.ticker import MultipleLocator
from collections import defaultdict
from collections import Counter

logging.basicConfig(
    format='%(asctime)s %(levelname)s:%(message)s',
    level=logging.INFO,
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logging.info("31_virtual_sky start")

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

def telo_condition(node : list, need_label_index:dict) -> bool:
    return node in need_label_index

# SKY Figure
def plot_virtual_chromosome(ax, segments_data, cut_dict, label=None):
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
    for i, (real_chr, seg_length) in enumerate(segments):
        color = CHR_COLORS.get(real_chr, "gray")
        rect = patches.Rectangle(
            (x_center - radius, current_y),
            width,
            seg_length,
            facecolor=color,
            edgecolor='none'
        )
        rect.set_clip_path(clip_patch)
        ax.add_patch(rect)

        current_y += seg_length
    
    text_obj_list = []
    for j in cut_dict[path]:
        for i in j:
            if len(i)==0:
                continue
            if len(i) > 2:
                real_chr = i[0]
                next_chr = i[1]
                if next_chr != real_chr:
                    # "chr" 접두어 제거 후 숫자나 문자만 추출 (예: chr1 -> 1, chrX -> X)
                    from_val = real_chr[3:] if real_chr.lower().startswith("chr") else real_chr
                    to_val = next_chr[3:] if next_chr.lower().startswith("chr") else next_chr
                    text_label = f"t({to_val};{from_val})"
                    # 텍스트는 가상 염색체 오른쪽 (x_center + radius + offset) 경계에 수직 중앙 정렬
                    text_x = x_center + radius * 2
                    text_obj = ax.text(text_x, i[2]/maxh*100, text_label, ha='left', va='center', fontsize=10, color='black')
                    text_obj_list.append(text_obj)
                else:
                    from_val = real_chr[3:] if real_chr.lower().startswith("chr") else real_chr
                    to_val = next_chr[3:] if next_chr.lower().startswith("chr") else next_chr
                    text_label = f"t({to_val+i[2]};{from_val+i[3]})"
                    text_x = x_center + radius * 2
                    text_obj = ax.text(text_x, i[4]/maxh*100, text_label, ha='left', va='center', fontsize=10, color='black')
                    text_obj_list.append(text_obj)
            else:
                ctg_idx = i[0]
                telo_dist = need_label_index[ctg_idx]
                text_label = f"{telo_dist[0]}, {round(ppc_data[ctg_idx][CHR_END]/M, 1)}Mb"
                text_x = x_center + radius * 2
                text_obj = ax.text(text_x, i[1]/maxh*100, text_label, ha='left', va='center', fontsize=10, color='black')
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

    ax.set_xlim(0, 60 * cell_col / def_cell_col)
    ax.set_ylim(0, 100)
    ax.axis('off')

def plot_scale_bar(ax, chr_name):
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

ppc_data = import_ppc_data(PREPROCESSED_PAF_FILE_PATH)

with open(f'{PREFIX}/path_data.pkl', 'rb') as f:
    path_list_dict = pkl.load(f)

df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
df = df.query('chr != "chrM"')

meandepth = np.median(df['meandepth'])
chr_len = find_chr_len(CHROMOSOME_INFO_FILE_PATH)

with open(f'{PREFIX}/contig_pat_vec_data.pkl', 'rb') as f:
    paf_ans_list, key_list, int2key = pkl.load(f)
paf_ans_dict = dict(paf_ans_list)

with open(f'{PREFIX}/depth_weight.pkl', 'rb') as f:
    tot_loc_list, weights = pkl.load(f)

with open(f'{PREFIX}/for_dir_data.pkl', 'rb') as f:
    for_dir_data = pkl.load(f)

# Normalize type 4
for def_ind, tind, w in for_dir_data:
    weights[def_ind] += weights[tind] * w

weights_sorted_data = sorted(enumerate(weights), key=lambda t:t[1], reverse=True)

path_dict = {}
non_type4_top_path = []
for ind, w in weights_sorted_data:
    paf_loc = tot_loc_list[ind]
    key = paf_loc.split('/')[-3]
    if paf_loc.split('/')[-3] != '11_ref_ratio_outliers' and w > meandepth / 10: # weight > 1/2N
        non_type4_top_path.append(paf_loc)
        path_dict[paf_loc] = w

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

cut_dict={}
karyotypes_data = {}
for path_path in non_type4_top_path:
    paf_combo = []
    path_contig_data = []
    curr_path_cut_list = []
    curr_path_inverse_list = []
    curr_path_telo_list = []

    tot_len = 0
    for idx, ki in enumerate(paf_ans_dict[path_path]):
        temp_list = import_paf_data(f'{output_folder}/{ki}.paf')
        curr_len = 0
        for i in temp_list:
            curr_len += i[CHR_END] - i[CHR_STR]
        if int2key[ki][0]!=TEL_TYPE:
            s = int2key[ki][1][0][1]
            e = int2key[ki][1][1][1]
            if idx == 0 and telo_condition(s, need_label_index):
                curr_path_telo_list.append((s, 0))
            if idx == len(paf_ans_dict[path_path])-1 and telo_condition(e, need_label_index):
                curr_path_telo_list.append((e, tot_len + curr_len))
            if ppc_data[s][CHR_NAM]!=ppc_data[e][CHR_NAM]:
                curr_path_cut_list.append((ppc_data[s][CHR_NAM], ppc_data[e][CHR_NAM], tot_len + curr_len//2))
            elif ppc_data[s][CTG_DIR] != ppc_data[e][CTG_DIR] and int2key[ki][0] == CTG_IN_TYPE:
                curr_path_inverse_list.append((ppc_data[s][CHR_NAM], ppc_data[e][CHR_NAM], ppc_data[s][CTG_DIR], ppc_data[e][CTG_DIR], tot_len + curr_len//2))
        else:
            s = int2key[ki][1]
            if telo_condition(s, need_label_index):
                curr_path_telo_list.append((s, tot_len + curr_len//2))
        tot_len += curr_len
        path_contig_data += temp_list
    cut_dict[path_path] = [curr_path_cut_list, curr_path_inverse_list, curr_path_telo_list]
            
    
    curr_chr = path_contig_data[0][CHR_NAM]
    curr_combo_len = path_contig_data[0][CHR_END] - path_contig_data[0][CHR_STR]
    curr_chr_range = [path_contig_data[0][CHR_STR], path_contig_data[0][CHR_END]]
    for i in range(1, len(path_contig_data)):
        if path_contig_data[i][CHR_NAM] != curr_chr:
            paf_combo.append((curr_chr, curr_combo_len))
            curr_chr = path_contig_data[i][CHR_NAM]
            curr_combo_len = path_contig_data[i][CHR_END] - path_contig_data[i][CHR_STR]
        else:
            curr_combo_len += path_contig_data[i][CHR_END] - path_contig_data[i][CHR_STR]
    paf_combo.append((curr_chr, curr_combo_len))
    karyotypes_data[path_path] = paf_combo

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

loc2weight = dict(zip(tot_loc_list, weights))

grouped_norm_data = defaultdict(list)
for path, data in karyotypes_norm_data.items():
    cnt = Counter()
    for c, w in data:
        cnt[c] += w
    sorted_cnt_data = sorted(cnt.items(), key=lambda t: -t[1])

    if sorted_cnt_data[0][1] / sum(cnt.values()) > MAJOR_BASELINE:
        grouped_norm_data[sorted_cnt_data[0][0]].append((path, data))
    else:
        grouped_norm_data['Mixed'].append((path, data))

cols = 10
rows = 0
for chr_name, data_list in grouped_norm_data.items():
    rows += (len(data_list) // cols + 1)

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
    bef_now_col = now_col
    plot_chr_name(ax_array[now_col][0], chr_name)
    for i, data in enumerate(data_list):
        plot_virtual_chromosome(ax_array[i // cols + now_col][i % cols + len(prefix_ratios)], data, cut_dict,
                                label=f"{round(float(loc2weight[data[0]] / meandepth * 2), 2)}N")

    now_col += (len(data_list) // cols + 1)
    for col in range(bef_now_col, now_col):
        plot_scale_bar(ax_array[col][len(prefix_ratios) - 1], chr_name)



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

fig.savefig(f'{PREFIX}/virtual_sky.pdf')
fig.savefig(f'{PREFIX}/virtual_sky.png')

logging.info('SKYPE pipeline end')
