import numpy as np
import pandas as pd
from collections import defaultdict 
import argparse

K = 1000
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

# --- main_df에 offset 적용 ---
# 1) ref의 offset을 (chr,st,nd) 기준으로 매핑할 dict 생성
offset_dict = {
    (c, s, e): off
    for c, s, e, off in zip(
        ref_df['chr'], ref_df['st'], ref_df['nd'], ref_df['offset']
    )
}

# 2) main_df에 보정값이 있으면 빼고, 없으면 0으로 처리
def apply_offset(row):
    key = (row['chr'], row['st'], row['nd'])
    # 1) 이 bin 이 repeat 구간에 속하면 보정하지 않음
    if is_in_repeat(row['chr'], row['st'], row['nd'], censat_repeat_data):
        return row['meandepth']
    # 2) 아닐 경우, ref offset 이 있으면 빼고, 없으면 0
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

