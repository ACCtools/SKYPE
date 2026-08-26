import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from skype_utils import *
from skype_output_files import (
    build_matrix_column_locations,
    discover_ecdna_depth_inputs,
    discover_jump_depth_inputs,
)
from nclose_tracking import (
    bnd_event_keys,
    compressed_bnd_event_keys,
    count_ecdna_circuit_events,
    count_index_path_events,
    event_catalog_by_key,
    indel_event_key,
    initialise_filter_status,
    load_event_catalog,
    load_type4_edge_event_map,
    save_filter_status,
    save_path_usage,
)

import numpy as np
import pandas as pd
import pickle as pkl

import h5py
import logging
import argparse
import collections
import glob

from collections import Counter, defaultdict

from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

# logging 설정(레벨/포맷)은 skype_utils 에서 중앙 관리한다 (LOG_LEVEL).
logging.info("22_save_matrix start")


ABS_MAX_COVERAGE_RATIO = 3



def import_index_path(file_path : str) -> list:
    file_path_list = file_path.split('/')
    key = file_path_list[-2]
    cnt = int(file_path_list[-1].split('.')[0]) - 1

    return path_list_dict[key][cnt][0]

def chr2int(x):
    chrXY2int = {'chrX': 24, 'chrY': 25}
    if x in chrXY2int:
        return chrXY2int[x]
    return int(x[3:])


def import_bed(bed_path: str) -> dict:
    bed_data_file = open(bed_path, "r")
    chr_len = collections.defaultdict(list)
    for curr_data in bed_data_file:
        curr_data = curr_data.split("\t")
        chr_len[curr_data[0]].append((int(curr_data[1]), int(curr_data[2])))
    bed_data_file.close()
    return chr_len


def inclusive_checker(tuple_a: tuple, tuple_b: tuple) -> bool:
    if int(tuple_a[0]) <= int(tuple_b[0]) and int(tuple_b[1]) <= int(tuple_a[1]):
        return True
    else:
        return False


parser = argparse.ArgumentParser(description="SKYPE depth analysis")

parser.add_argument("censat_bed_path",
                    help="Path to the censat repeat information file.")

parser.add_argument("main_stat_loc",
                    help="Cancer coverage location file")

parser.add_argument("prefix",
                    help="Pefix for pipeline")

parser.add_argument("-t", "--thread",
                    help="Number of thread", type=int)

parser.add_argument("--progress",
                    help="Show progress bar", action='store_true')

args = parser.parse_args()

bed_data = import_bed(args.censat_bed_path)

PREFIX = args.prefix
THREAD = args.thread
main_stat_loc = args.main_stat_loc

RATIO_OUTLIER_FOLDER = f"{PREFIX}/11_ref_ratio_outliers/"
front_contig_path = RATIO_OUTLIER_FOLDER + "front_jump/"
back_contig_path = RATIO_OUTLIER_FOLDER + "back_jump/"
ecdna_contig_path = RATIO_OUTLIER_FOLDER + "ecdna/"
type2_ins_contig_path = RATIO_OUTLIER_FOLDER + "type2_ins/"

df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t',
                 names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
df = df.query('chr != "chrM"')

nclose_nodes_pkl = load_nclose_nodes(PREFIX)

nclose_set = set()
for vl in nclose_nodes_pkl.values():
    for v in vl:
        nclose_set.add(tuple(sorted(v)))

nclose_event_catalog = load_event_catalog(PREFIX)
nclose_event_by_key = event_catalog_by_key(nclose_event_catalog)
nclose_bnd_keys = bnd_event_keys(nclose_event_catalog)
compressed_nclose_bnd_keys = compressed_bnd_event_keys(nclose_event_catalog)
if nclose_set != compressed_nclose_bnd_keys:
    raise ValueError(
        "NClose BND catalog does not match nclose_nodes.pkl. "
        "Rerun SKYPE from its NClose preprocessing stage."
    )
type4_edge_to_event = load_type4_edge_event_map(PREFIX, nclose_event_catalog)

with open(f'{PREFIX}/path_data.pkl', 'rb') as f:
    path_list_dict = pkl.load(f)

ydf = df.query('chr == "chrY"')

bed_intervals = pd.IntervalIndex.from_tuples(bed_data['chrY'], closed='left')
y_mask = ydf.apply(
    lambda row: bed_intervals.overlaps(pd.Interval(row['st'], row['nd'], closed='left')),
    axis=1
)

correct_mask = y_mask.apply(any)
ydf_not_censat = ydf[~correct_mask]

chry_nz_len = len(ydf_not_censat.query('meandepth != 0'))
no_chrY = (chry_nz_len / len(ydf_not_censat)) < chrY_MINIMUM_RATIO

meandepth = np.median(df['meandepth'])
chr_filt_st_list = []
chr_no_filt_st_list = []

for l in df.itertuples(index=False):
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
    else:
        chr_no_filt_st_list.append((l.chr, l.st))

filter_len = len(chr_filt_st_list)

def get_vec_from_stat_loc(stat_loc_):
    df = pd.read_csv(stat_loc_, compression='gzip', comment='#', sep='\t',
                     names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
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
    df = pd.read_csv(stat_loc_, compression='gzip', comment='#', sep='\t',
                     names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
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


PATH_FILE_FOLDER = f"{PREFIX}/20_depth"
chr_chr_folder_path = sorted(glob.glob(PATH_FILE_FOLDER + "/*"))

chr_st_data = dict()
for l in df.itertuples(index=False):
    if l.chr == 'chrY' and no_chrY:
        chr_st_data[(l.chr, l.st)] = 0
    else:
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

B = np.asarray(v, dtype=np.float32)

with open(f'{PREFIX}/contig_pat_vec_data.pkl', 'rb') as f:
    paf_ans_list, key_list, _, dep_list = pkl.load(f)

vec_dict = [None] * (len(key_list) + 1)
output_folder = f'{PREFIX}/21_pat_depth'
with ProcessPoolExecutor(max_workers=THREAD) as executor:
    futures = []
    for ki in key_list:
        futures.append(executor.submit(get_vec_from_ki, ki))
    for future in tqdm(as_completed(futures), total=len(futures), desc='Parse depth for each seperated paths\' gz file',
                       disable=not sys.stdout.isatty() and not args.progress):
        i, v = future.result()
        vec_dict[i] = v

front_depth_inputs = discover_jump_depth_inputs(
    front_contig_path, type2_ins_contig_path
)
back_depth_inputs = discover_jump_depth_inputs(
    back_contig_path, type2_ins_contig_path
)
ecdna_depth_inputs = discover_ecdna_depth_inputs(ecdna_contig_path)

with open(f'{PREFIX}/ecdna_circuit_data.pkl', 'rb') as f:
    ecdna_circuits, _ = pkl.load(f)
if len(ecdna_circuits) != len(ecdna_depth_inputs):
    raise ValueError(
        "ecDNA circuit/depth-input count mismatch: "
        f"{len(ecdna_circuits)} != {len(ecdna_depth_inputs)}. "
        "Rerun SKYPE from stage 21."
    )

with open(f'{PREFIX}/cen_fragment_data.pkl', 'rb') as f:
    cen_fragment_meta = pkl.load(f)
cen_fragment_list = sorted(cen_fragment_meta.items(), key=lambda kv: chr2int(kv[0]))

m = np.shape(B)[0]
n = (
    len(paf_ans_list)
    + len(front_depth_inputs)
    + len(back_depth_inputs)
    + len(ecdna_depth_inputs)
    + len(cen_fragment_list)
)
ncnt = 0

shape = (n, m)
A_arr = np.empty(shape, dtype=np.float32, order='C')
fm = filter_len

filter_vec_list = []
tot_loc_list = build_matrix_column_locations(
    PREFIX,
    paf_ans_list,
    front_depth_inputs,
    back_depth_inputs,
    ecdna_depth_inputs,
    cen_fragment_list,
)
if len(tot_loc_list) != n:
    raise ValueError(
        "Stage-22 matrix-column location mismatch: "
        f"{len(tot_loc_list)} locations for {n} matrix columns"
    )

tmp_v = np.zeros(m, dtype=np.float32)

ncnt = 0
# Matrix column index -> canonical event tags.
# ordinary nclose: (left_contig_idx, right_contig_idx), sorted tuple of ints
# step-11 INDEL: (event_type, event_idx, type2_merge_idx)
# ecDNA columns carry their two compressed BND keys.
path_nclose_usage = [Counter() for _ in range(n)]
explicit_filter_reasons = {}
for path, key_int_list in tqdm(paf_ans_list, desc='Recover depth from separated paths',
                                    disable=not sys.stdout.isatty() and not args.progress):
    ki = key_int_list[0]
    np.copyto(tmp_v, vec_dict[ki])
    for ki in key_int_list[1:]:
        tmp_v += vec_dict[ki]

    idx_path = import_index_path(path)

    usage = count_index_path_events(
        idx_path, compressed_nclose_bnd_keys, type4_edge_to_event
    )
    path_nclose_usage[ncnt].update(usage)
    A_arr[ncnt, :m] = tmp_v
    ncnt += 1

with open(f'{PREFIX}/indel_exclude_idx_set.pkl', 'rb') as f:
    indel_exclude_idx_set = pkl.load( f)
tv_empty = np.zeros(len(chr_filt_st_list) + len(chr_no_filt_st_list), dtype=np.float32)

indel_idx = 0
# Process forward-directed outlier contigs
for i, ov_loc, bv_loc, type2_ins_idx, type2_ins_loc in tqdm(
              front_depth_inputs,
              desc='Parse coverage from forward-directed outlier contig gz files',
              disable=not sys.stdout.isatty() and not args.progress):
    ov = get_vec_from_stat_loc(ov_loc)
    if type2_ins_loc is not None:
        ov += get_vec_from_stat_loc(type2_ins_loc)
    bv = get_vec_from_stat_loc(bv_loc)
    
    event_key = indel_event_key('front_jump', i, type2_ins_idx)
    if event_key not in nclose_event_by_key:
        raise ValueError(f"Missing step-11 INDEL catalog event: {event_key}")
    if indel_idx in indel_exclude_idx_set:
        tv = tv_empty
        explicit_filter_reasons[event_key] = 'FILTERED_02_EXCLUDE_LIST'
    else:
        tv = ov - bv
    
    A_arr[ncnt, :m] = tv
        
    path_nclose_usage[ncnt][event_key] += 1
    ncnt += 1
    indel_idx += 1
# Process backward-directed outlier contigs
for i, ov_loc, bv_loc, type2_ins_idx, type2_ins_loc in tqdm(
              back_depth_inputs,
              desc='Parse coverage from backward-directed outlier contig gz files',
              disable=not sys.stdout.isatty() and not args.progress):
    ov = get_vec_from_stat_loc(ov_loc)
    if type2_ins_loc is not None:
        ov += get_vec_from_stat_loc(type2_ins_loc)
    bv = get_vec_from_stat_loc(bv_loc)

    event_key = indel_event_key('back_jump', i, type2_ins_idx)
    if event_key not in nclose_event_by_key:
        raise ValueError(f"Missing step-11 INDEL catalog event: {event_key}")
    if indel_idx in indel_exclude_idx_set:
        tv = tv_empty
        explicit_filter_reasons[event_key] = 'FILTERED_02_EXCLUDE_LIST'
    else:
        tv = ov + bv
    
    A_arr[ncnt, :m] = tv
        
    path_nclose_usage[ncnt][event_key] += 1
    ncnt += 1
    indel_idx += 1

for i, ov_loc in ecdna_depth_inputs:
    ov = get_vec_from_stat_loc(ov_loc)
    type2_ins_idx = -1

    A_arr[ncnt, :m] = ov
        
    ecdna_usage = count_ecdna_circuit_events(
        ecdna_circuits[i - 1], nclose_bnd_keys
    )
    path_nclose_usage[ncnt].update(ecdna_usage)
    ncnt += 1

# Centromere depth-difference fragment chromosomes
# 02번에서 detect 된 (chrom, dir, mid_censat) 기반으로 fragment 1개 = column 1개.
# 단위 indicator vector(0/1)로 채우면 NNLS weight ≈ 추가 copy 수.
all_st_list = chr_filt_st_list + chr_no_filt_st_list
for chrom, info in cen_fragment_list:
    mid = info["mid"]
    chr_length = info["chr_len"]
    if info["dir"]:
        st_lo, st_hi = mid, chr_length
    else:
        st_lo, st_hi = 0, mid

    tv = np.zeros(m, dtype=np.float32)
    for idx, (c, st) in enumerate(all_st_list):
        if c == chrom and st_lo <= st < st_hi:
            tv[idx] = 1.0

    A_arr[ncnt, :m] = tv

    ncnt += 1

dep_list.extend([0] * (
    len(front_depth_inputs)
    + len(back_depth_inputs)
    + len(ecdna_depth_inputs)
    + len(cen_fragment_list)
))

assert(len(dep_list) == n)
assert ncnt == n

# INDEL_INDEX_* exclusions are resolved completely while constructing the
# matrix. Any matrix column carrying an explicitly excluded event becomes a
# zero column. Stage 23 still passes every matrix column to the single raw
# solve, so the explicit exclusion remains visible in provenance.
explicit_excluded_columns = {
    col_idx
    for event_key in explicit_filter_reasons
    for col_idx, usage in enumerate(path_nclose_usage)
    if usage.get(event_key, 0) > 0
}
if explicit_excluded_columns:
    A_arr[np.asarray(sorted(explicit_excluded_columns)), :] = 0.0
    logging.info(
        "Zeroed %d matrix columns carrying explicitly excluded INDEL events",
        len(explicit_excluded_columns),
    )

initial_active_columns = set(range(n)) - explicit_excluded_columns

nclose_filter_status = initialise_filter_status(
    nclose_event_catalog,
    path_nclose_usage,
    initial_active_columns,
    explicit_filter_reasons,
)
save_path_usage(PREFIX, path_nclose_usage)
save_filter_status(PREFIX, nclose_filter_status)
logging.info(
    "NClose path usage : %d events across %d/%d matrix columns",
    len({key for usage in path_nclose_usage for key in usage}),
    sum(bool(usage) for usage in path_nclose_usage),
    len(path_nclose_usage),
)

MATRIX_CONTRACT = DEPTH_ONLY_MATRIX_CONTRACT

with h5py.File(f'{PREFIX}/matrix.h5', 'w') as hf:
    hf.attrs['matrix_contract'] = MATRIX_CONTRACT
    hf.create_dataset('A', data=A_arr[:, :fm])
    hf.create_dataset('A_fail', data=A_arr[:, fm:m])
    hf.create_dataset('B', data=B[:fm])
    hf.create_dataset('B_fail', data=B[fm:])
    hf.create_dataset('B_depth_start', data=0)
    hf.create_dataset('B_depth_end', data=fm)

with open(f"{PREFIX}/23_input.pkl", "wb") as f:
    pkl.dump({
        "matrix_contract": MATRIX_CONTRACT,
        "chr_filt_st_list": chr_filt_st_list,
        "B_depth_start": 0,
        "B_depth_end": fm,
    }, f)

with open(f"{PREFIX}/tot_loc_list.pkl", "wb") as f:
    pkl.dump(tot_loc_list, f)

np.save(f'{PREFIX}/B.npy', B)
