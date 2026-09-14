import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from skype_utils import *

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from collections import defaultdict 
import argparse

CENSAT_COMPRESSABLE_THRESHOLD = 1000*K
INF = 1e9

def chr2int(x):
    if x.startswith('chr'):
        chrXY2int = {'chrX' : 24, 'chrY' : 25}
        if x in chrXY2int:
            return chrXY2int[x]
        else:
            return int(x[3:])
    else:
        return INF

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



def is_in_repeat(chr_name: str, st: int, nd: int, repeat_dict: dict) -> bool:
    """
    주어진 (chr, st, nd) 구간이 repeat_dict[chr] 안의
    어떤 (r_start, r_end) 구간과라도 겹치면 True 반환.
    """
    if chr_name == 'chrY':
        return True
    for r_start, r_end in repeat_dict.get(chr_name, []):
        # 겹친다는 건: not (bin 완전 좌측 or bin 완전 우측)
        if not (nd < r_start or st > r_end):
            return True
    return False

parser = argparse.ArgumentParser(description="Find breakend contigs with contig data and map data")
parser.add_argument("main_stat_loc", type=str, help="Path to the main depth statistics file")
parser.add_argument("ref_stat_loc", type=str, help="Path to the reference depth statistics file")
parser.add_argument("censat_path", type=str, help="Path to the Censat repeat data file")
args = parser.parse_args()

main_stat_loc = args.main_stat_loc
ref_stat_loc = args.ref_stat_loc
censat_path = args.censat_path
censat_repeat_data = import_censat_repeat_data(censat_path)

main_df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
main_df2 = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
ref_df = pd.read_csv(ref_stat_loc, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
ref_meandepth = np.median(ref_df['meandepth'])
main_df = main_df[main_df['chr'] != 'chrM']
ref_df = ref_df[ref_df['chr'] != 'chrM']

main_df = main_df[main_df['meandepth'] != 'meandepth']
ref_df  = ref_df[ ref_df['meandepth']  != 'meandepth']
main_df = main_df[main_df['st'] != 'st']
ref_df  = ref_df[ ref_df['st']  != 'st']
main_df = main_df[main_df['nd'] != 'nd']
ref_df  = ref_df[ ref_df['nd']  != 'nd']
main_df = main_df[main_df['totaldepth'] != 'totaldepth']
ref_df  = ref_df[ ref_df['totaldepth']  != 'totaldepth']

for df in (main_df, ref_df):
    df['st']         = pd.to_numeric(df['st'], errors='coerce').astype('Int64')
    df['nd']         = pd.to_numeric(df['nd'], errors='coerce').astype('Int64')
    df['meandepth']  = pd.to_numeric(df['meandepth'], errors='coerce').astype(float)
dfs = [main_df, ref_df]

ref_df = ref_df.assign(
    offset = ref_df['meandepth'] / ref_meandepth
)

# Smooth only the log control factors: sigma=1 bin (100 kb), kernel +/-4 bins.
# Keep chromosome boundaries, coordinate gaps, and normalization exemptions
# separate so they cannot contribute to a neighboring factor.
offsets = ref_df['offset'].to_numpy(copy=True)
excluded = np.array([
    is_in_repeat(row.chr, row.st, row.nd, censat_repeat_data)
    for row in ref_df.itertuples(index=False)
], dtype=bool)
eligible = (~excluded) & np.isfinite(offsets) & (offsets > 0)
indices = np.flatnonzero(eligible)
chroms = ref_df['chr'].to_numpy()
starts = ref_df['st'].to_numpy(dtype=np.int64)
ends = ref_df['nd'].to_numpy(dtype=np.int64)
breaks = np.flatnonzero(
    (np.diff(indices) != 1)
    | (chroms[indices[1:]] != chroms[indices[:-1]])
    | (starts[indices[1:]] != ends[indices[:-1]] + 1)
) + 1
for block in np.split(indices, breaks):
    if len(block):
        offsets[block] = np.exp(gaussian_filter1d(
            np.log(offsets[block]), sigma=1.0, mode='reflect', truncate=4.0
        ))
ref_df['offset'] = offsets

# --- main_df에 offset 적용 ---
# 1) ref의 offset을 (chr,st,nd) 기준으로 매핑할 dict 생성
offset_dict = {
    (c, s, e): off
    for c, s, e, off in zip(
        ref_df['chr'], ref_df['st'], ref_df['nd'], ref_df['offset']
    )
}

# 2) main_df를 보정 계수로 나누고, 대응 control이 없으면 계수 1 사용
def apply_offset(row):
    key = (row['chr'], row['st'], row['nd'])
    # 1) 이 bin 이 repeat 구간에 속하면 보정하지 않음
    if is_in_repeat(row['chr'], row['st'], row['nd'], censat_repeat_data):
        return row['meandepth']
    # 2) 아닐 경우, smoothed ref offset으로 나누고, 없으면 계수 1 사용
    return row['meandepth'] / offset_dict.get(key, 1.0)

main_df['meandepth'] = main_df.apply(apply_offset, axis=1).astype(float)

main_df['totaldepth'] = (main_df['meandepth'] * 100000).round().astype(int)


cols = ['chr','st','nd','length','covsite','totaldepth','cov','meandepth']
output_path = "/".join(main_stat_loc.split("/")[:-1]) + f"/{main_stat_loc.split("/")[-1].split(".")[0]}_normalized.win.stat.gz"
main_df.to_csv(
    output_path,
    sep='\t',
    columns=cols,
    header=False,
    index=False,
    compression='gzip'
)
