import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import h5py
import glob
import psutil
import logging
import warnings
import argparse

import pickle as pkl
import numpy as np
import pandas as pd

from types import SimpleNamespace
from collections import defaultdict, Counter
from scipy.stats import mannwhitneyu, ttest_ind

from adelie.solver import bvls

from h5py._hl import selections as sel

from skype_utils import *

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
CTG_CENSAT = 18
CTG_MAINFLOWDIR = 19
CTG_MAINFLOWCHR = 20
CTG_GLOBALIDX = 21

FILTER_P_VALUE = 0.01

CHROM_ERROR_FAIL_RATE = 2.0

RECIPROCAL_PAIR_DISTANCE = 10 * K
RECIPROCAL_FORBID_MAX_INTERVAL = 10 * K

RAW_TRANSLOCATION_RESULT_PKL = 'raw_translocation_result.pkl'

def chr2int(x):
    chrXY2int = {'chrX': 24, 'chrY': 25}
    if x in chrXY2int:
        return chrXY2int[x]
    return int(x[3:])

def import_ppc_contig_data(file_path: str) -> list:
    paf_file = open(file_path, "r")
    contig_data = []
    for curr_contig in paf_file:
        temp_list = curr_contig.split("\t")
        int_induce_idx = [CTG_LEN, CTG_STR, CTG_END, \
                          CHR_LEN, CHR_STR, CHR_END, \
                          CTG_MAPQ, CTG_TYP, CTG_STRND, CTG_ENDND, ]
        for i in int_induce_idx:
            temp_list[i] = int(temp_list[i])
        contig_data.append(tuple(temp_list))
    paf_file.close()
    return contig_data

# Hypothesis test nclose filtering

def get_contig_pair(i, cd):
    nclose_loc = i == 0
    nclose_dir = '+' == cd[CTG_DIR]

    chr_name = cd[CHR_NAM]
    chr_cord = cd[CHR_STR if nclose_dir ^ nclose_loc else CHR_END]

    return chr_name, chr_cord

def get_expected_high_side(i, cd):
    nclose_loc = i == 0
    nclose_dir = '+' == cd[CTG_DIR]

    # If the breakend is at CHR_STR, the included/reference side is the
    # higher-coordinate side; if it is at CHR_END, it is the lower side.
    return 'right' if nclose_dir ^ nclose_loc else 'left'

def get_chr_cord_data(contig_data, nclose_nodes, type4_nclose_nodes, type4_expected_high_side):
    chr_cord_dict = defaultdict(list)
    path_tag2cord = {}
    path_tag2nclose_key = {}
    path_tag2expected_high_side = {}

    def add_nclose_cords(path_tag, nclose_key, cords, expected_high_side):
        path_tag2cord[path_tag] = cords
        path_tag2nclose_key[path_tag] = nclose_key
        path_tag2expected_high_side[path_tag] = expected_high_side

        for side_idx, (chr_name, chr_cord) in enumerate(cords):
            chr_cord_dict[chr_name].append((chr_cord, (path_tag, side_idx)))
    
    for ctg_name, nclose_idx_pair_list in nclose_nodes.items():
        for j, nclose_idx_pair in enumerate(nclose_idx_pair_list):
            path_tag = (ctg_name, j)
            nclose_key = tuple(nclose_idx_pair)
            assert len(nclose_key) == 2 and nclose_key == tuple(sorted(nclose_key))
            cords = [get_contig_pair(i, contig_data[idx]) for i, idx in enumerate(nclose_idx_pair)]
            expected_high_side = [
                get_expected_high_side(i, contig_data[idx])
                for i, idx in enumerate(nclose_idx_pair)
            ]
            add_nclose_cords(path_tag, nclose_key, cords, expected_high_side)
                
    
    for type4_nclose_key, (c1, i1, c2, i2) in type4_nclose_nodes.items():
        add_nclose_cords(
            type4_nclose_key, type4_nclose_key, [(c1, i1), (c2, i2)],
            type4_expected_high_side.get(type4_nclose_key, [None, None])
        )


    for l in chr_cord_dict.values():
        l.sort(key=lambda t: t[0])

    return chr_cord_dict, path_tag2cord, path_tag2nclose_key, path_tag2expected_high_side

def build_chr_cord_dict(path_tag2cord, path_tag2nclose_key, skip_nclose_set=None):
    skip_nclose_set = skip_nclose_set or set()
    chr_cord_dict = defaultdict(list)
    for path_tag, cords in path_tag2cord.items():
        if path_tag2nclose_key[path_tag] in skip_nclose_set:
            continue
        for side_idx, (chr_name, chr_cord) in enumerate(cords):
            chr_cord_dict[chr_name].append((chr_cord, (path_tag, side_idx)))

    for l in chr_cord_dict.values():
        l.sort(key=lambda t: t[0])

    return chr_cord_dict

def collect_nclose_test_data(contig_data, nclose_nodes):
    type4_tot_loc_list = []
    for bv_paf_loc in glob.glob(f'{RATIO_OUTLIER_FOLDER}/front_jump/*_base.paf'):
        type4_tot_loc_list.append(bv_paf_loc)

    for bv_paf_loc in glob.glob(f'{RATIO_OUTLIER_FOLDER}/back_jump/*_base.paf'):
        type4_tot_loc_list.append(bv_paf_loc)

    type4_nclose_nodes = dict()
    type4_expected_high_side = dict()

    def parse_type4_event_breakend_cords(paf_loc):
        paf_rows = []
        with open(paf_loc, 'r') as f:
            for line in f:
                if line.strip():
                    paf_rows.append(line.rstrip().split('\t'))
        if not paf_rows:
            return None

        def event_breakend_coord(row, side_idx):
            nclose_loc = side_idx == 0
            nclose_dir = row[CTG_DIR] == '+'
            coord_idx = CHR_STR if nclose_dir ^ nclose_loc else CHR_END
            return row[CHR_NAM], int(row[coord_idx])

        cords = [
            event_breakend_coord(paf_rows[0], 0),
            event_breakend_coord(paf_rows[-1], 1),
        ]
        cords.sort(key=lambda x: (chr2int(x[0]), x[1]))
        return cords

    for paf_loc in type4_tot_loc_list:
        event_type = paf_loc.split('/')[-2]
        event_idx = int(paf_loc.split('/')[-1].split('_')[0])

        type2_merge_paf_idx = -1
        event_paf_loc = f'{RATIO_OUTLIER_FOLDER}/{event_type}/{event_idx}.paf'
        if not os.path.isfile(event_paf_loc):
            type2_merge_paf_loc = glob.glob(f'{RATIO_OUTLIER_FOLDER}/{event_type}/{event_idx}_type2_merge_*.paf')[0]
            type2_merge_paf_idx = int(type2_merge_paf_loc.split('/')[-1].split('.')[0].split('_')[-1])
            event_paf_loc = type2_merge_paf_loc

        cords = parse_type4_event_breakend_cords(event_paf_loc)
        if cords is None:
            continue

        type4_tag = (event_type, event_idx, type2_merge_paf_idx)
        type4_nclose_nodes[type4_tag] = (cords[0][0], cords[0][1], cords[1][0], cords[1][1])
        if event_type == 'front_jump':
            # deletion: the outside of the breakend span should be higher.
            type4_expected_high_side[type4_tag] = ['left', 'right']
        elif event_type == 'back_jump':
            # insertion/duplication: the inside of the breakend span should be higher.
            type4_expected_high_side[type4_tag] = ['right', 'left']

    return get_chr_cord_data(contig_data, nclose_nodes, type4_nclose_nodes, type4_expected_high_side)

def build_path_nclose_count(weights_len, ordinary_candidate_nclose_keys, type4_candidate_nclose_keys,
                            path_nclose_set_dict, path_list_dict, paf_sort_ans_list):
    path_nclose_count = [Counter() for _ in range(weights_len)]
    rpll = len(paf_sort_ans_list)

    for col_idx in range(rpll):
        paf_loc = paf_sort_ans_list[col_idx][0]
        key = paf_loc.split('/')[-2]
        cnt = int(paf_loc.split('/')[-1].split('.')[0]) - 1
        idx_path = path_list_dict[key][cnt][0]
        ctr = path_nclose_count[col_idx]
        s = 1
        while s < len(idx_path) - 2:
            cand = tuple(sorted([idx_path[s][1], idx_path[s+1][1]]))
            if cand in ordinary_candidate_nclose_keys:
                ctr[cand] += 1
                s += 2
            else:
                s += 1

    for col_idx in range(weights_len):
        ctr = Counter()
        for nclose_key in path_nclose_set_dict.get(col_idx, set()):
            if nclose_key in type4_candidate_nclose_keys:
                ctr[nclose_key] += 1
            elif col_idx >= rpll and nclose_key in ordinary_candidate_nclose_keys:
                ctr[nclose_key] += 1
        path_nclose_count[col_idx].update(ctr)

    return path_nclose_count

def calculate_nclose_depth(path_nclose_count, weight_full):
    nclose_depth = defaultdict(float)
    for col_idx, ctr in enumerate(path_nclose_count):
        for nc, v in ctr.items():
            nclose_depth[nc] += v * weight_full[col_idx]
    return nclose_depth

def canonical_nclose_layout(cords, expected_high_side):
    order = [0, 1]
    keys = [(chr2int(cords[i][0]), cords[i][1]) for i in order]
    if keys[1] < keys[0]:
        order = [1, 0]

    return {
        'chroms': tuple(cords[i][0] for i in order),
        'coords': tuple(cords[i][1] for i in order),
        'sides': tuple(expected_high_side[i] if i < len(expected_high_side) else None for i in order),
    }

def find_opposite_direction_depth_drop(path_tag2cord, path_tag2nclose_key,
                                       path_tag2expected_high_side, nclose_depth,
                                       max_depth_diff,
                                       distance_threshold=1 * M):
    layout_by_key = {}
    for path_tag, cords in path_tag2cord.items():
        nclose_key = path_tag2nclose_key[path_tag]
        if len(nclose_key) != 2:
            continue
        layout = canonical_nclose_layout(
            cords,
            path_tag2expected_high_side.get(path_tag, [None, None])
        )
        if None in layout['sides']:
            continue
        layout_by_key[nclose_key] = layout

    groups = defaultdict(list)
    for nclose_key, layout in layout_by_key.items():
        groups[layout['chroms']].append((nclose_key, layout))

    drop_set = set()
    dropped_pair_count = 0
    for items in groups.values():
        for i in range(len(items)):
            key_i, layout_i = items[i]
            depth_i = nclose_depth.get(key_i, 0.0)
            for j in range(i + 1, len(items)):
                key_j, layout_j = items[j]
                if layout_i['sides'][0] == layout_j['sides'][0] or layout_i['sides'][1] == layout_j['sides'][1]:
                    continue
                if abs(layout_i['coords'][0] - layout_j['coords'][0]) > distance_threshold:
                    continue
                if abs(layout_i['coords'][1] - layout_j['coords'][1]) > distance_threshold:
                    continue

                depth_j = nclose_depth.get(key_j, 0.0)
                if max(depth_i, depth_j) < 1e-6:
                    continue
                if abs(depth_i - depth_j) < max_depth_diff:
                    drop_set.add(key_i)
                    drop_set.add(key_j)
                    dropped_pair_count += 1

    return drop_set, dropped_pair_count

def canonical_nclose_region_layout(contig_data, nclose_key, expected_high_side):
    regions = []
    cords = []
    for side_idx, contig_idx in enumerate(nclose_key):
        cd = contig_data[contig_idx]
        cords.append(get_contig_pair(side_idx, cd))
        regions.append((
            cd[CHR_NAM],
            int(cd[CHR_STR]),
            int(cd[CHR_END]),
            cd[CTG_NAM],
            cd[CTG_DIR],
        ))

    order = [0, 1]
    keys = [(chr2int(cords[i][0]), cords[i][1]) for i in order]
    if keys[1] < keys[0]:
        order = [1, 0]

    return {
        'chroms': tuple(cords[i][0] for i in order),
        'coords': tuple(cords[i][1] for i in order),
        'regions': tuple(regions[i] for i in order),
        'sides': tuple(expected_high_side[i] if i < len(expected_high_side) else None for i in order),
    }

def get_inner_boundary_interval(region_a, region_b):
    chrom_a, st_a, nd_a, _, _ = region_a
    chrom_b, st_b, nd_b, _, _ = region_b
    assert chrom_a == chrom_b

    if (st_b, nd_b) < (st_a, nd_a):
        st_a, nd_a, st_b, nd_b = st_b, nd_b, st_a, nd_a

    if nd_a <= st_b:
        return nd_a, st_b, 'gap'

    return st_b, min(nd_a, nd_b), 'overlap'

def get_outer_reference_interval(region_a, region_b):
    chrom_a, st_a, nd_a, _, _ = region_a
    chrom_b, st_b, nd_b, _, _ = region_b
    assert chrom_a == chrom_b
    return min(st_a, st_b), max(nd_a, nd_b)

def load_utg_spans(paf_utg_path):
    utg_spans = defaultdict(list)
    if not paf_utg_path or not os.path.isfile(paf_utg_path):
        logging.warning(f'Reciprocal translocation check skipped raw utg PAF load: {paf_utg_path}')
        return utg_spans

    with open(paf_utg_path, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) <= CHR_END:
                continue
            qname = fields[CTG_NAM]
            if not qname.startswith('utg'):
                continue
            try:
                st = int(fields[CHR_STR])
                nd = int(fields[CHR_END])
            except ValueError:
                continue
            if nd < st:
                continue
            utg_spans[fields[CHR_NAM]].append({
                'qname': qname,
                'st': st,
                'nd': nd,
            })

    for spans in utg_spans.values():
        spans.sort(key=lambda x: (x['st'], x['nd'], x['qname']))

    return utg_spans

def find_strict_spanning_utgs(utg_spans, chrom, st, nd):
    hits = []
    for span in utg_spans.get(chrom, []):
        if span['st'] >= st:
            break
        if span['nd'] > nd:
            hits.append(span)
    return hits

def build_censat_interval_dict(cdf):
    censat_intervals = defaultdict(list)
    for chrom, st, nd in cdf[['chr', 'st', 'nd']].itertuples(index=False, name=None):
        censat_intervals[chrom].append((int(st), int(nd)))

    for chrom in censat_intervals:
        censat_intervals[chrom].sort()

    return censat_intervals

def overlaps_censat(censat_intervals, chrom, st, nd):
    for censat_st, censat_nd in censat_intervals.get(chrom, []):
        if nd <= st:
            if censat_st <= st < censat_nd:
                return True
            if censat_st > st:
                break
            continue

        if censat_nd <= st:
            continue
        if censat_st >= nd:
            break
        return True

    return False

def build_reciprocal_no_span_records(contig_data, path_tag2cord, path_tag2nclose_key,
                                     path_tag2expected_high_side, utg_spans, censat_intervals,
                                     distance_threshold=RECIPROCAL_PAIR_DISTANCE):
    layout_by_key = {}
    for path_tag, cords in path_tag2cord.items():
        nclose_key = path_tag2nclose_key[path_tag]
        if len(nclose_key) != 2:
            continue
        expected_high_side = path_tag2expected_high_side.get(path_tag, [None, None])
        layout = canonical_nclose_layout(cords, expected_high_side)
        if None in layout['sides']:
            continue
        layout.update(canonical_nclose_region_layout(contig_data, nclose_key, expected_high_side))
        layout_by_key[nclose_key] = layout

    groups = defaultdict(list)
    for nclose_key, layout in layout_by_key.items():
        groups[layout['chroms']].append((nclose_key, layout))

    records = []
    pair_id = 0
    for chrom_pair, items in groups.items():
        for i in range(len(items)):
            key_i, layout_i = items[i]
            for j in range(i + 1, len(items)):
                key_j, layout_j = items[j]
                if layout_i['sides'][0] == layout_j['sides'][0] or layout_i['sides'][1] == layout_j['sides'][1]:
                    continue
                if abs(layout_i['coords'][0] - layout_j['coords'][0]) > distance_threshold:
                    continue
                if abs(layout_i['coords'][1] - layout_j['coords'][1]) > distance_threshold:
                    continue

                side_records = []
                for side_idx, chrom in enumerate(chrom_pair):
                    inner_st, inner_nd, _ = get_inner_boundary_interval(
                        layout_i['regions'][side_idx],
                        layout_j['regions'][side_idx],
                    )
                    path_drop_st, path_drop_nd = get_outer_reference_interval(
                        layout_i['regions'][side_idx],
                        layout_j['regions'][side_idx],
                    )
                    inner_len = max(0, inner_nd - inner_st)
                    no_spanning_utg = False
                    if (
                        inner_len <= RECIPROCAL_FORBID_MAX_INTERVAL and
                        not overlaps_censat(censat_intervals, chrom, inner_st, inner_nd)
                    ):
                        no_spanning_utg = len(find_strict_spanning_utgs(utg_spans, chrom, inner_st, inner_nd)) == 0
                    side_records.append({
                        'side_id': f'R{pair_id}.{side_idx}',
                        'chrom': chrom,
                        'inner_st': inner_st,
                        'inner_nd': inner_nd,
                        'path_drop_st': path_drop_st,
                        'path_drop_nd': path_drop_nd,
                        'no_spanning_utg': no_spanning_utg,
                        'crossing_cols': [],
                        'accepted_forbid': False,
                    })

                records.append({
                    'pair_id': pair_id,
                    'nclose_key_a': key_i,
                    'nclose_key_b': key_j,
                    'chrom_pair': chrom_pair,
                    'layout_a': layout_i,
                    'layout_b': layout_j,
                    'side_records': side_records,
                    'pair_protected_from_1M_drop': False,
                })
                pair_id += 1

    return records

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

def side_inner_interval(side):
    return (
        int(side.get('inner_st', side['path_drop_st'])),
        int(side.get('inner_nd', side['path_drop_nd'])),
    )

def same_chrom_partner_by_side_id(records):
    partner_by_id = {}
    for record in records:
        side_records = record.get('side_records', [])
        if len(side_records) != 2:
            continue
        chrom_pair = record.get('chrom_pair')
        if chrom_pair and len(chrom_pair) >= 2:
            same_chrom = chrom_pair[0] == chrom_pair[1]
        else:
            same_chrom = side_records[0].get('chrom') == side_records[1].get('chrom')
        if not same_chrom:
            continue
        partner_by_id[side_records[0]['side_id']] = side_records[1]
        partner_by_id[side_records[1]['side_id']] = side_records[0]
    return partner_by_id

def same_chrom_contiguous_span(merged_intervals, side, partner):
    side_st, side_nd = side_inner_interval(side)
    partner_st, partner_nd = side_inner_interval(partner)
    span_st = min(side_st, partner_st)
    span_nd = max(side_nd, partner_nd)
    return interval_strictly_contains_any(merged_intervals, span_st, span_nd)

def assign_reciprocal_crossing_path_cols(records, prefix, paf_sort_ans_list):
    forbidden_by_chrom = defaultdict(list)
    side_by_id = {}
    for record in records:
        for side in record['side_records']:
            side_by_id[side['side_id']] = side
            side.setdefault('same_chrom_spanning_cols', [])
            if side_is_forbidden_by_raw_read(side):
                forbidden_by_chrom[side['chrom']].append(side)

    if not forbidden_by_chrom:
        return side_by_id

    partner_by_id = same_chrom_partner_by_side_id(records)
    component_cache = {}
    target_chroms = set(forbidden_by_chrom)
    for col_idx, (_, key_int_list) in enumerate(paf_sort_ans_list):
        col_intervals_by_chrom = defaultdict(list)
        for key_int in key_int_list:
            key_intervals = read_component_ref_intervals(prefix, key_int, component_cache)
            for chrom in target_chroms:
                col_intervals_by_chrom[chrom].extend(key_intervals.get(chrom, []))

        for chrom, side_list in forbidden_by_chrom.items():
            merged = merge_ref_intervals(col_intervals_by_chrom.get(chrom, []))
            if not merged:
                continue
            for side in side_list:
                path_drop_st = side['path_drop_st']
                path_drop_nd = side['path_drop_nd']
                if not interval_strictly_contains_any(merged, path_drop_st, path_drop_nd):
                    continue
                partner = partner_by_id.get(side['side_id'])
                if partner is not None and same_chrom_contiguous_span(merged, side, partner):
                    side['same_chrom_spanning_cols'].append(col_idx)
                    continue
                side['crossing_cols'].append(col_idx)

    same_chrom_forbidden_sides = [
        side for side in side_by_id.values()
        if side['side_id'] in partner_by_id and side_is_forbidden_by_raw_read(side)
    ]
    same_chrom_drop_cols = {
        col for side in same_chrom_forbidden_sides
        for col in side.get('crossing_cols', [])
    }
    same_chrom_spanning_cols = {
        col for side in same_chrom_forbidden_sides
        for col in side.get('same_chrom_spanning_cols', [])
    }
    if same_chrom_drop_cols or same_chrom_spanning_cols:
        logging.info(
            'Reciprocal same-chrom no-span path classification : '
            f'drop_cols={len(same_chrom_drop_cols)}, '
            f'contiguous_span_cols={len(same_chrom_spanning_cols)}'
        )

    return side_by_id

def side_is_forbidden_by_raw_read(side):
    return bool(side.get('no_spanning_rawread', side.get('no_spanning_utg', False)))

def record_has_raw_point_no_span(record):
    no_span = record.get('raw_point_no_spanning')
    if isinstance(no_span, dict):
        return bool(no_span.get('point_a', False) or no_span.get('point_b', False))

    side_flags = [
        side.get('raw_point_no_spanning')
        for side in record.get('side_records', [])
        if 'raw_point_no_spanning' in side
    ]
    return any(bool(flag) for flag in side_flags)

def load_raw_translocation_records(prefix):
    result_path = f'{prefix}/{RAW_TRANSLOCATION_RESULT_PKL}'
    if not os.path.isfile(result_path):
        logging.warning(f'Raw-read translocation result not found: {result_path}')
        return []

    with open(result_path, 'rb') as f:
        records = pkl.load(f)

    for record in records:
        for side in record.get('side_records', []):
            side.setdefault('crossing_cols', [])
            side.setdefault('accepted_forbid', False)
            if 'no_spanning_rawread' not in side:
                side['no_spanning_rawread'] = bool(side.get('no_spanning_utg', False))

    return records

def get_chrom_acc_sum_max(chr_filt_st_list, predict_a, predict_b):
    chrom_acc_sum_dict = defaultdict(float)
    chrom_acc_sum_dict_max = defaultdict(float)
    for i, (chrom, _) in enumerate(chr_filt_st_list):
        chrom_acc_sum_dict[chrom] += predict_a[i] - predict_b[i]
        if abs(chrom_acc_sum_dict[chrom]) > chrom_acc_sum_dict_max[chrom]:
            chrom_acc_sum_dict_max[chrom] = abs(chrom_acc_sum_dict[chrom])
    return chrom_acc_sum_dict_max

def get_fail_chrom_list(chrom_acc_sum_dict_max, chrom_acc_sum_dict_max_base, no_chrY):
    fail_chrom_list = []
    for chrom, acc_sum_max in chrom_acc_sum_dict_max.items():
        if chrom == 'chrY' and no_chrY:
            continue
        acc_sum_max_base = chrom_acc_sum_dict_max_base[chrom]
        if acc_sum_max_base == 0:
            chrom_error_rate = float('inf') if acc_sum_max > 0 else 0
        else:
            chrom_error_rate = acc_sum_max / acc_sum_max_base
        if chrom_error_rate > CHROM_ERROR_FAIL_RATE:
            fail_chrom_list.append(chrom)
    return fail_chrom_list

def group_close_coordinates(data_list, distance_threshold=50 * K):
    """
    Groups a list of tuples by their coordinate values.
    
    This function assumes the input data_list is already sorted by the 
    coordinate (the first element of the tuple). Tuples are grouped together 
    if the distance between consecutive coordinates is less than or equal 
    to the distance_threshold (default is 100k).
    
    Args:
        data_list: A list of tuples formatted as (cord, ctg_name, nclose_idx).
        distance_threshold: Maximum distance between consecutive coordinates to be grouped.
        
    Returns:
        A list of grouped tuples formatted as (group_start_cord, [(ctg_name, nclose_idx), ...]).
    """
    if not data_list:
        return []

    result = []
    
    # Initialize the first group using the first element in the list
    first_cord = data_list[0][0]
    last_cord = data_list[0][0]
    current_group_items = [data_list[0][1]]

    for i in range(1, len(data_list)):
        cord, nclose_data = data_list[i]

        # Check if the current coordinate is within the threshold from the previous one
        if cord - last_cord <= distance_threshold:
            current_group_items.append(nclose_data)
        else:
            # Save the completed group and start a new one
            result.append((first_cord, current_group_items))
            first_cord = cord
            current_group_items = [nclose_data]
            
        # Update the last seen coordinate for the next iteration
        last_cord = cord

    # Append the final group after the loop finishes
    result.append((first_cord, current_group_items))

    return result

def merge_regions_iteratively(df_cov, censat_intervals,
                              target_chr, chr_cord_dict,
                              p_value):
    """
    Iteratively merges adjacent regions based on t-test p-values.
    It prevents merging only if the specific split point falls within any excluded areas.
    """
    split_points = group_close_coordinates(chr_cord_dict[target_chr])
    out_split_points = []
    
    # Filter dataframes for the specific chromosome to speed up operations
    df_cov_chr = df_cov[df_cov['chr'] == target_chr].copy()
    
    # Filter out the excluded regions from the coverage data upfront
    # Bins that overlap with any excluded region are completely removed
    keep_mask = [
        not overlaps_censat(censat_intervals, target_chr, int(st), int(nd))
        for st, nd in df_cov_chr[['st', 'nd']].itertuples(index=False, name=None)
    ]
    df_cov_chr = df_cov_chr[keep_mask]
        
    # Copy the split points to safely modify the active list during iteration
    active_splits = split_points.copy()
    
    # Determine the absolute maximum boundary from the remaining valid coverage data
    max_nd = df_cov_chr['nd'].max() if not df_cov_chr.empty else 0
    
    while True:
        merged_in_this_round = False
        
        for i in range(len(active_splits)):
            current_split = active_splits[i][0]
            
            # Check if the specific split point itself falls within any excluded region
            # Skip the merge process if the exact point is located inside an excluded area
            is_split_point_excluded = overlaps_censat(
                censat_intervals, target_chr, current_split, current_split + 1
            )
            
            if is_split_point_excluded:
                continue
            
            # Define boundaries for the left and right regions
            left_start = active_splits[i-1][0] if i > 0 else 0
            left_end = current_split
            
            right_start = current_split
            right_end = active_splits[i+1][0] if i < len(active_splits) - 1 else max_nd
            
            # Extract coverage values for left and right regions based on bin coordinates
            # Since df_cov_chr has already been filtered, excluded regions naturally won't be included
            left_cov = list(df_cov_chr[
                (df_cov_chr['nd'] > left_start) & 
                (df_cov_chr['st'] < left_end)
            ]['meandepth'])
            
            right_cov = list(df_cov_chr[
                (df_cov_chr['nd'] > right_start) & 
                (df_cov_chr['st'] < right_end)
            ]['meandepth'])
            
            # Ensure both regions have at least 2 data points to perform the statistical test
            # This naturally avoids statistical calculation errors
            if len(left_cov) >= 5 and len(right_cov) >= 5:
                # stat, p_val = ttest_ind(left_cov, right_cov)
                stat, p_val = mannwhitneyu(left_cov, right_cov)

                # If distributions are not significantly different, merge them
                if p_val >= p_value and (not np.isnan(p_val)):
                    out_split_points.append(active_splits[i])
                    active_splits.pop(i)
                    merged_in_this_round = True
                    break # Break the inner loop to restart the while loop with updated regions
                    
        # If no merges occurred in the entire pass, the regions are stable
        if not merged_in_this_round:
            break
            
    return active_splits, out_split_points

def get_censet_count(censat_intervals, c1, i1, c2, i2):
    cnt = 0
    for chr_name, chr_cord in [(c1, i1), (c2, i2)]:
        if overlaps_censat(censat_intervals, chr_name, chr_cord, chr_cord + 1):
            cnt += 1

    return cnt

def filter_nclose_by_test(contig_data, nclose_nodes, df, censat_intervals,
                          pre_drop_nclose_set=None, nclose_test_data=None):
    if nclose_test_data is None:
        nclose_test_data = collect_nclose_test_data(contig_data, nclose_nodes)

    _, path_tag2cord, path_tag2nclose_key, path_tag2expected_high_side = nclose_test_data
    chr_cord_dict = build_chr_cord_dict(path_tag2cord, path_tag2nclose_key, pre_drop_nclose_set)

    in_cnt = defaultdict(lambda : [0, 0])
    nclose_key_side_out = defaultdict(lambda: [False, False])
    for chr_key in chr_cord_dict:
        in_split, out_split = merge_regions_iteratively(
            df, censat_intervals, chr_key, chr_cord_dict, p_value=FILTER_P_VALUE
        )

        for ctg_idx, nclose_side_list in in_split:
            for path_tag, side_idx in nclose_side_list:
                in_cnt[path_tag][0] += 1

        for ctg_idx, nclose_side_list in out_split:
            for path_tag, side_idx in nclose_side_list:
                in_cnt[path_tag][1] += 1
                nclose_key = path_tag2nclose_key[path_tag]
                nclose_key_side_out[nclose_key][side_idx] = True

    return in_cnt, path_tag2cord, path_tag2nclose_key, nclose_key_side_out, path_tag2expected_high_side

def _group_consecutive(sorted_idx_list):
    """
    Group a sorted list of indices into (start, end) ranges
    where end is inclusive.
    e.g. [0,1,2,5,6,9] -> [(0,2), (5,6), (9,9)]
    """
    if not sorted_idx_list:
        return []
    ranges = []
    start = sorted_idx_list[0]
    prev = start
    for idx in sorted_idx_list[1:]:
        if idx == prev + 1:
            prev = idx
        else:
            ranges.append((start, prev))
            start = idx
            prev = idx
    ranges.append((start, prev))
    return ranges


def read_selected_feature_rows_direct(dA, A_store, A_idx_list):
    """
    Read specific feature rows from a transposed HDF5 matrix into the top rows
    of a pre-allocated NumPy array using the low-level HDF5 API.

    A_idx_list must be sorted. Consecutive indices are batched into
    contiguous slice reads for much better I/O performance.
    """
    ranges = _group_consecutive(A_idx_list)

    dest_row = 0
    for start, end in ranges:
        nrows = end - start + 1

        # Source: contiguous feature-row slice in the file
        source_sel = sel.select(dA.shape, np.s_[start:end + 1, :], dA)
        fspace = source_sel.id

        # Destination: corresponding slice in the output array
        dest_sel = sel.select(A_store.shape, np.s_[dest_row:dest_row + nrows, :])

        for mspace in dest_sel.broadcast(source_sel.array_shape):
            dA.id.read(mspace, fspace, A_store, dxpl=dA._dxpl)

        dest_row += nrows


def solver_matrix_from_store(A_store, n_features=None):
    if n_features is not None:
        A_store = A_store[:n_features, :]
    return A_store.T


def make_bvls_warm_start(beta, lower, upper):
    beta = np.asarray(beta, dtype=lower.dtype)
    if beta.shape != lower.shape:
        raise ValueError(
            f"warm_start has shape {beta.shape}, expected {lower.shape}"
        )
    beta = np.maximum(beta, lower)
    beta = np.minimum(beta, upper).copy()
    is_active = (lower < beta) & (beta < upper)
    active_set_size = int(is_active.sum())
    active_set = np.empty(beta.shape[0], dtype=int)
    active_set[:active_set_size] = np.flatnonzero(is_active)
    return SimpleNamespace(
        beta=beta,
        active_set=active_set,
        active_set_size=active_set_size,
        is_active=is_active,
    )


def fit_nnls_warm(X, y, warm_start=None):
    p = X.shape[1]
    lower = np.zeros(p, dtype=X.dtype)
    upper = np.full(p, np.finfo(X.dtype).max, dtype=X.dtype)
    warm_state = None
    if warm_start is not None:
        warm_state = make_bvls_warm_start(warm_start, lower, upper)

    state = bvls(
        X,
        y,
        lower,
        upper,
        n_threads=1,
        warm_start=warm_state,
    )
    return state.beta


def scatter_subset_weights(weights, A_idx_list, n_features, dtype=None):
    out = np.zeros(n_features, dtype=dtype or weights.dtype)
    out[A_idx_list] = weights
    return out

# logging 설정(레벨/포맷)은 skype_utils 에서 중앙 관리한다 (LOG_LEVEL).
logging.info("23_run_nnls start")

parser = argparse.ArgumentParser(description="SKYPE depth analysis")
parser.add_argument("ppc_paf_file_path", 
                    help="Path to the preprocessed PAF file.")

parser.add_argument("prefix", help="Prefix for pipeline")

parser.add_argument("main_stat_path", help="Path to the main depth statistics file")

parser.add_argument("censat_bed_path", 
                    help="Path to the censat repeat information file.")

parser.add_argument("paf_utg_path",
                    help="Path to the raw unitig PAF file.")

parser.add_argument("-t", "--thread", help="Number of threads", type=int)

args = parser.parse_args()

# t = "23_run_nnls.py /Data/hyunwoo/00_skype_run_data/HG008-T/20_alignasm/HG008-T.ctg.aln.paf.ppc.paf /hyunwoo/63_skype_test/skype/HG008-T_18_08_06 /Data/hyunwoo/00_skype_run_data/HG008-T/01_depth/HG008-T_normalized.win.stat.gz /hyunwoo/63_skype_test/deps/SKYPE/public_data/chm13v2.0_censat_v2.1.m.bed -t 72"
# t = "23_run_nnls.py /Data/hyunwoo/00_skype_run_data/COLO829/20_alignasm/COLO829.ctg.aln.paf.ppc.paf /hyunwoo/63_skype_test/skype/COLO829_15_06_54 /Data/hyunwoo/00_skype_run_data/COLO829/01_depth/COLO829_normalized.win.stat.gz /hyunwoo/63_skype_test/deps/SKYPE/public_data/chm13v2.0_censat_v2.1.m.bed -t 16"
# args = parser.parse_args(t.split()[1:])

PREPROCESSED_PAF_FILE_PATH = args.ppc_paf_file_path
PREFIX = args.prefix
main_stat_loc = args.main_stat_path
THREAD = args.thread
CENSAT_PATH = args.censat_bed_path
PAF_UTG_LOC = args.paf_utg_path

RATIO_OUTLIER_FOLDER = f"{PREFIX}/11_ref_ratio_outliers/"

with open(f"{PREFIX}/report.txt", 'r') as f:
    f.readline()
    path_cnt = int(f.readline().strip())

use_julia_solver = path_cnt <= HARD_PATH_COUNT_BASELINE

core_num = psutil.cpu_count(logical=False)
if core_num is None:
    THREAD = THREAD
else:
    THREAD = min(int(THREAD), core_num)

contig_data = import_ppc_contig_data(PREPROCESSED_PAF_FILE_PATH)

with open(f'{PREFIX}/nclose_chunk_data.pkl', 'rb') as f:
    nclose_nodes, _, _ = pkl.load(f)

with open(f'{PREFIX}/cen_fragment_data.pkl', 'rb') as f:
    cen_fragment_meta = pkl.load(f)

df = pd.read_csv(main_stat_loc, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])
df = df.query('chr != "chrM"')
meandepth = np.median(df['meandepth'])
N = meandepth / 2

cdf = pd.read_csv(CENSAT_PATH, sep='\t', names=['chr', 'st', 'nd'])

with open(f'{PREFIX}/23_input.pkl', 'rb') as f:
    chr_filt_st_list, path_nclose_set_dict, amplitude, bed_data = pkl.load(f)
# path_nclose_set_dict is keyed by matrix column index. Values are canonical event tags:
# ordinary nclose: (left_contig_idx, right_contig_idx), sorted tuple of ints
# type4/ecdna: (event_type, event_idx, type2_merge_idx)
# cent_fragment: ('cent_fragment', chrom, direction_bool)

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

with h5py.File(f"{PREFIX}/matrix.h5", "r") as f:
    dA = f["A"]
    A_store = np.empty(dA.shape, dtype=dA.dtype)
    dA.read_direct(A_store)

    dB = f["B"]
    B = np.empty(dB.shape, dtype=dB.dtype)
    dB.read_direct(B)

    dAf = f["A_fail"]
    A_fail_store = np.empty(dAf.shape, dtype=dAf.dtype)
    dAf.read_direct(A_fail_store)

    b_start_ind = int(f["B_depth_start"][()])

A = solver_matrix_from_store(A_store)

weights_fullsize = fit_nnls_warm(A, B)
predict_suc_B_base = A.dot(weights_fullsize)

with open(f'{PREFIX}/path_data.pkl', 'rb') as f:
    path_list_dict = pkl.load(f)
with open(f'{PREFIX}/contig_pat_vec_data.pkl', 'rb') as f:
    paf_sort_ans_list, _, _, _ = pkl.load(f)

rpll = len(paf_sort_ans_list)
assert len(weights_fullsize) >= rpll

nclose_test_data = collect_nclose_test_data(contig_data, nclose_nodes)
_, base_path_tag2cord, base_path_tag2nclose_key, base_path_tag2expected_high_side = nclose_test_data
base_candidate_nclose_keys = set(base_path_tag2nclose_key.values())
base_ordinary_candidate_nclose_keys = {
    nclose_key for nclose_key in base_candidate_nclose_keys
    if len(nclose_key) == 2
}
base_type4_candidate_nclose_keys = {
    nclose_key for nclose_key in base_candidate_nclose_keys
    if len(nclose_key) == 3
}
path_nclose_count = build_path_nclose_count(
    len(weights_fullsize), base_ordinary_candidate_nclose_keys, base_type4_candidate_nclose_keys,
    path_nclose_set_dict, path_list_dict, paf_sort_ans_list
)

base_nclose_depth = calculate_nclose_depth(path_nclose_count, weights_fullsize)
chrom_acc_sum_dict_max_base = get_chrom_acc_sum_max(chr_filt_st_list, predict_suc_B_base, B)

censat_intervals = build_censat_interval_dict(cdf)
reciprocal_records = load_raw_translocation_records(PREFIX)
reciprocal_side_by_id = assign_reciprocal_crossing_path_cols(
    reciprocal_records, PREFIX, paf_sort_ans_list
)

active_reciprocal_side_ids = {
    side['side_id']
    for record in reciprocal_records
    if record_has_raw_point_no_span(record)
    for side in record.get('side_records', [])
    if side_is_forbidden_by_raw_read(side)
}
reciprocal_path_drop_cols = set()
reciprocal_retry_iter = 0
while active_reciprocal_side_ids:
    reciprocal_retry_iter += 1
    next_drop_cols = set()
    for side_id in active_reciprocal_side_ids:
        next_drop_cols.update(reciprocal_side_by_id[side_id]['crossing_cols'])

    A_idx_test = [
        col_idx for col_idx in range(len(weights_fullsize))
        if col_idx not in next_drop_cols
    ]
    if not A_idx_test:
        active_reciprocal_side_ids = set()
        reciprocal_path_drop_cols = set()
        break

    with h5py.File(f"{PREFIX}/matrix.h5", "r") as f:
        dA = f["A"]
        read_selected_feature_rows_direct(dA, A_store, A_idx_test)

    A_test = solver_matrix_from_store(A_store, len(A_idx_test))
    reciprocal_weights = fit_nnls_warm(
        A_test, B, warm_start=weights_fullsize[A_idx_test]
    )
    reciprocal_predict_suc_B = A_test.dot(reciprocal_weights)

    chrom_acc_sum_dict_max = get_chrom_acc_sum_max(
        chr_filt_st_list, reciprocal_predict_suc_B, predict_suc_B_base
    )
    reciprocal_fail_chrom_list = get_fail_chrom_list(
        chrom_acc_sum_dict_max, chrom_acc_sum_dict_max_base, no_chrY
    )
    if not reciprocal_fail_chrom_list:
        reciprocal_path_drop_cols = next_drop_cols
        break

    fail_chrom_set = set(reciprocal_fail_chrom_list)
    remove_side_ids = {
        side_id for side_id in active_reciprocal_side_ids
        if reciprocal_side_by_id[side_id]['chrom'] in fail_chrom_set
    }
    if not remove_side_ids:
        logging.warning(
            'Reciprocal translocation retry failed on chromosomes without matching '
            f'candidate intervals: {",".join(sorted(fail_chrom_set))}'
        )
        active_reciprocal_side_ids = set()
        reciprocal_path_drop_cols = set()
        break

    active_reciprocal_side_ids -= remove_side_ids
    reciprocal_path_drop_cols = set()

for side_id in active_reciprocal_side_ids:
    reciprocal_side_by_id[side_id]['accepted_forbid'] = True

reciprocal_transloc_nclose_set = set()
for record in reciprocal_records:
    if (
        record_has_raw_point_no_span(record) and
        any(side['accepted_forbid'] for side in record['side_records'])
    ):
        record['pair_protected_from_1M_drop'] = True
        reciprocal_transloc_nclose_set.add(record['nclose_key_a'])
        reciprocal_transloc_nclose_set.add(record['nclose_key_b'])

logging.info(
    f'Reciprocal translocation no-span candidates : {len(active_reciprocal_side_ids)}'
)

opposite_dir_depth_drop_set, opposite_dir_depth_pair_count = find_opposite_direction_depth_drop(
    base_path_tag2cord, base_path_tag2nclose_key, base_path_tag2expected_high_side,
    base_nclose_depth, max_depth_diff=0.1 * N
)
opposite_dir_depth_drop_set -= reciprocal_transloc_nclose_set
logging.info(
    f'Opposite-direction depth prefilter nclose count : {len(opposite_dir_depth_drop_set)} '
    f'(protected reciprocal nclose count: {len(reciprocal_transloc_nclose_set)})'
)

# Loop start
prev_fail_chrom_list = None
fail_chrom_list = []
nclose_filter_warm_fullsize = weights_fullsize.copy()

in_cnt, path_tag2cord, path_tag2nclose_key, nclose_key_side_out, path_tag2expected_high_side = \
    filter_nclose_by_test(
        contig_data, nclose_nodes, df, censat_intervals,
        pre_drop_nclose_set=opposite_dir_depth_drop_set,
        nclose_test_data=nclose_test_data
    )

while True:
    not_using_nclose_set = set(opposite_dir_depth_drop_set)
    hyp_test_nclose_set = set()

    for path_tag, nclose_cnt in in_cnt.items():
        (c1, i1), (c2, i2) = path_tag2cord[path_tag]
        nclose_key = path_tag2nclose_key[path_tag]
        
        if nclose_key in reciprocal_transloc_nclose_set:
            continue

        censat_count = get_censet_count(censat_intervals, c1, i1, c2, i2)
        if (not (censat_count == 2 or nclose_cnt[0] - censat_count >= 1) and
            c1 not in fail_chrom_list and c2 not in fail_chrom_list):
            hyp_test_nclose_set.add(nclose_key)

    not_using_nclose_set |= hyp_test_nclose_set

    fail_chrom_text = ','.join(fail_chrom_list) if fail_chrom_list else 'none'
    logging.info(f'Hyp. test filtering nclose count : {len(hyp_test_nclose_set)} (fail chroms: {fail_chrom_text})')
    logging.debug(
        f'Hyp. test filtering nclose breakdown : '
        f'{dict(Counter("ordinary" if len(tag) == 2 else tag[0] for tag in hyp_test_nclose_set))}'
    )

    A_idx_list = []
    for k, path_nclose_list in path_nclose_set_dict.items():
        if k in reciprocal_path_drop_cols:
            continue
        if not path_nclose_list or len(not_using_nclose_set & set(path_nclose_list)) == 0:
            A_idx_list.append(k)

    A_idx_list.sort()

    with h5py.File(f"{PREFIX}/matrix.h5", "r") as f:
        dA = f["A"]
        dAf = f["A_fail"]

        read_selected_feature_rows_direct(dA, A_store, A_idx_list)
        read_selected_feature_rows_direct(dAf, A_fail_store, A_idx_list)

    A_new = solver_matrix_from_store(A_store, len(A_idx_list))
    A_fail_new = solver_matrix_from_store(A_fail_store, len(A_idx_list))

    filter_weights = fit_nnls_warm(
        A_new, B, warm_start=nclose_filter_warm_fullsize[A_idx_list]
    )
    nclose_filter_warm_fullsize = scatter_subset_weights(
        filter_weights, A_idx_list, len(weights_fullsize), dtype=weights_fullsize.dtype
    )

    predict_suc_B = A_new.dot(filter_weights)
    predict_fal_B = A_fail_new.dot(filter_weights)

    chrom_acc_sum_dict_max = get_chrom_acc_sum_max(
        chr_filt_st_list, predict_suc_B, predict_suc_B_base
    )
    new_fail_chrom_list = get_fail_chrom_list(
        chrom_acc_sum_dict_max, chrom_acc_sum_dict_max_base, no_chrY
    )

    # 새로 발견된 fail chrom을 누적
    fail_chrom_list = sorted(set(fail_chrom_list + new_fail_chrom_list))

    # Break: 새로 추가된 fail chrom이 없으면 종료
    if not new_fail_chrom_list or (prev_fail_chrom_list is not None and fail_chrom_list == prev_fail_chrom_list):
        break

    prev_fail_chrom_list = fail_chrom_list[:]

# === Cent_fragment greedy off-test ===
# nclose 필터링 직후를 fail 기준점(고정 baseline)으로 잡고, weight 작은 cent_fragment
# 컬럼부터 누적 제거 시도. PASS 시 컬럼 누적 제거(=offset 누적)되며, fail 비교 대상
# baseline 은 항상 nclose 필터링 직후로 고정 — 누적 drift 도 base 대비 한도 안이면 OK.
CF_GREEDY_MIN_WEIGHT_N = 0.1
CF_GREEDY_MIN_WEIGHT = CF_GREEDY_MIN_WEIGHT_N * N

cent_fragment_col2chrom = {}
for k, tags in path_nclose_set_dict.items():
    for tag in tags:
        if isinstance(tag, tuple) and len(tag) > 0 and tag[0] == 'cent_fragment':
            cent_fragment_col2chrom[k] = tag[1]
            break

def _compute_acc_max(predict_a, predict_b):
    acc = defaultdict(float)
    acc_max = defaultdict(float)
    for i, (chrom, _) in enumerate(chr_filt_st_list):
        acc[chrom] += predict_a[i] - predict_b[i]
        if abs(acc[chrom]) > acc_max[chrom]:
            acc_max[chrom] = abs(acc[chrom])
    return acc_max

# 고정 baseline (nclose-filter 직후) — 누적 제거 후 비교 대상
predict_suc_B_post = predict_suc_B.copy()
cf_acc_max_base_fixed = _compute_acc_max(predict_suc_B_post, B)

# 누적 상태 (PASS마다 갱신)
A_idx_list_curr = sorted(A_idx_list)
predict_suc_B_curr = predict_suc_B.copy()
predict_fal_B_curr = predict_fal_B.copy()
filter_weights_curr = filter_weights.copy()

logging.info(
    f'CF greedy start: pool={len(cent_fragment_col2chrom)} cent_fragment cols '
    f'(candidate weight >= {CF_GREEDY_MIN_WEIGHT_N}N)'
)

removed_cf = []
while True:
    pos_of = {col: i for i, col in enumerate(A_idx_list_curr)}
    pool = []
    for col, chrom in cent_fragment_col2chrom.items():
        if col not in pos_of:
            continue
        w = filter_weights_curr[pos_of[col]]
        if w >= CF_GREEDY_MIN_WEIGHT:
            pool.append((w, col, chrom))
    if not pool:
        break
    pool.sort(key=lambda x: x[0])

    progressed = False
    for w, col, chrom in pool:
        A_idx_test = sorted(set(A_idx_list_curr) - {col})
        with h5py.File(f"{PREFIX}/matrix.h5", "r") as f:
            read_selected_feature_rows_direct(f["A"], A_store, A_idx_test)
            read_selected_feature_rows_direct(f["A_fail"], A_fail_store, A_idx_test)
        
        A_test = solver_matrix_from_store(A_store, len(A_idx_test))
        A_fail_test = solver_matrix_from_store(A_fail_store, len(A_idx_test))

        cf_warm_fullsize = scatter_subset_weights(
            filter_weights_curr, A_idx_list_curr, len(weights_fullsize),
            dtype=weights_fullsize.dtype
        )
        weights_test = fit_nnls_warm(
            A_test, B, warm_start=cf_warm_fullsize[A_idx_test]
        )
        predict_suc_B_test = A_test.dot(weights_test)
        predict_fal_B_test = A_fail_test.dot(weights_test)

        # numerator: 누적 제거된 test prediction vs FIXED baseline (nclose-filter 직후)
        acc_max_diff = _compute_acc_max(predict_suc_B_test, predict_suc_B_post)

        ok = True
        worst_chrom = None
        worst_ratio = 0.0
        for chk, max_diff in acc_max_diff.items():
            if chk == 'chrY' and no_chrY:
                continue
            base_max = cf_acc_max_base_fixed.get(chk, 0.0)
            if base_max == 0:
                continue
            ratio = max_diff / base_max
            if ratio > CHROM_ERROR_FAIL_RATE:
                ok = False
                if ratio > worst_ratio:
                    worst_ratio = ratio
                    worst_chrom = chk

        if ok:
            logging.debug(f'CF PASS: drop col {col} chrom={chrom} w={w:.4f}')
            A_idx_list_curr = A_idx_test
            predict_suc_B_curr = predict_suc_B_test
            predict_fal_B_curr = predict_fal_B_test
            filter_weights_curr = weights_test
            removed_cf.append((col, chrom))
            progressed = True
            break
        else:
            logging.debug(f'CF FAIL: keep col {col} chrom={chrom} w={w:.4f} '
                         f'(worst chrom={worst_chrom} ratio={worst_ratio:.3f})')

    if not progressed:
        break

if removed_cf:
    logging.info(f'CF greedy result chromosomes: {", ".join(chrom for _, chrom in removed_cf)}')

A_idx_list = A_idx_list_curr
filter_weights = filter_weights_curr
predict_suc_B = predict_suc_B_curr
predict_fal_B = predict_fal_B_curr

# === BP step / depth ratio 기반 nclose 추가 필터링 ===
# predict_suc_B와 B 각각에서, ordinary nclose와 type4 indel의 두 BP 점 양옆 ±500kb를 확인.
# ordinary nclose 는 side/방향에서, type4 는 front_jump(del)/back_jump(ins)에서
# 예상되는 high-depth 쪽까지 확인한다.
# 각 side가 hyp-test out_split이거나 BP step/depth ratio 기준을 통과하고, 양쪽 side가 모두 통과하면 drop set에 추가.
# 최종 제거 대상은 predict_suc_B 기반 drop set과 B 기반 drop set의 union.
# BP 가 censat 영역이면 cen_fragment_meta 방향으로 판정하고, 방향 정보가 없거나
# ±500kb window 못 만들면 BP 기준은 false 처리 (=hyp-test out_split만 반영).
# cent_fragment greedy off-test 가 끝나 weight 가 정리된 상태에서 적용 (cent_fragment 의 잔여 영향 배제).
STEP_WIN = 500 * K
STEP_RATIO_THRESHOLD = 0.25

candidate_nclose_keys = set(path_tag2nclose_key.values())
ordinary_candidate_nclose_keys = {nclose_key for nclose_key in candidate_nclose_keys if len(nclose_key) == 2}
type4_candidate_nclose_keys = {nclose_key for nclose_key in candidate_nclose_keys if len(nclose_key) == 3}

# 각 nclose 의 두 BP (chr, cord)
nclose_cords = {}
nclose_expected_high_side = {}
for path_tag, cords in path_tag2cord.items():
    nclose_key = path_tag2nclose_key[path_tag]
    if nclose_key in candidate_nclose_keys:
        nclose_cords[nclose_key] = cords
        nclose_expected_high_side[nclose_key] = path_tag2expected_high_side.get(path_tag, [None, None])

# 현재 filter_weights 를 전체 column 길이로 unfold
cf_weight_fullsize = np.zeros_like(weights_fullsize)
cf_weight_fullsize[A_idx_list] = filter_weights

# nclose 별 depth = Σ (count × weight)
cf_nclose_depth = calculate_nclose_depth(path_nclose_count, cf_weight_fullsize)

# chr -> sorted [(bin_idx, st), ...]
chr_to_bins = defaultdict(list)
for i, (chrom, st) in enumerate(chr_filt_st_list):
    chr_to_bins[chrom].append((i, st))
for chrom in chr_to_bins:
    chr_to_bins[chrom].sort(key=lambda x: x[1])

def in_censat(chrom, cord):
    return overlaps_censat(censat_intervals, chrom, cord, cord + 1)

def censat_expected_high_side(chrom, cord):
    if not in_censat(chrom, cord):
        return None
    info = cen_fragment_meta.get(chrom)
    if info is None:
        return None
    return 'right' if info['dir'] else 'left'

def censat_is_droppable_by_ratio(chrom, cord, depth_val, expected_high_side):
    if expected_high_side is None:
        return False
    actual_high_side = censat_expected_high_side(chrom, cord)
    if actual_high_side is None:
        return False
    depth_diff = cen_fragment_meta[chrom]['depth_diff']
    step = depth_diff if expected_high_side == actual_high_side else -depth_diff
    return (step / depth_val) < STEP_RATIO_THRESHOLD

def signal_is_droppable_by_ratio(signal_B, left_bins, right_bins, depth_val, expected_high_side=None):
    left_rows = [b_start_ind + i for i in left_bins]
    right_rows = [b_start_ind + i for i in right_bins]
    left_mean = signal_B[left_rows].mean()
    right_mean = signal_B[right_rows].mean()
    if expected_high_side == 'left':
        step = left_mean - right_mean
    elif expected_high_side == 'right':
        step = right_mean - left_mean
    else:
        step = abs(right_mean - left_mean)
    return (step / depth_val) < STEP_RATIO_THRESHOLD

def side_passes_ratio(chrom, cord, depth_val, signal_B, expected_high_side=None):
    if in_censat(chrom, cord):
        return censat_is_droppable_by_ratio(chrom, cord, depth_val, expected_high_side)

    bins = chr_to_bins.get(chrom, [])
    left = [i for i, st in bins if cord - STEP_WIN <= st < cord]
    right = [i for i, st in bins if cord <= st < cord + STEP_WIN]
    if not left or not right:
        return False
    return signal_is_droppable_by_ratio(signal_B, left, right, depth_val, expected_high_side)

def side_is_droppable(nclose_key, side_idx, chrom, cord, depth_val, signal_B, expected_high_side=None):
    ttest_out = nclose_key_side_out[nclose_key][side_idx]
    return ttest_out or side_passes_ratio(chrom, cord, depth_val, signal_B, expected_high_side)

def find_nclose_to_drop(signal_B, nclose_depth):
    nclose_to_drop = set()
    for nclose_key, cords in nclose_cords.items():
        if nclose_key in reciprocal_transloc_nclose_set:
            continue
        
        depth = nclose_depth.get(nclose_key, 0.0)
        if depth < 1e-6:
            continue
        expected_high_side = nclose_expected_high_side.get(nclose_key, [None] * len(cords))
        if all(
            side_is_droppable(
                nclose_key, side_idx, chrom, pos, depth, signal_B,
                expected_high_side[side_idx] if side_idx < len(expected_high_side) else None
            )
            for side_idx, (chrom, pos) in enumerate(cords)
        ):
            nclose_to_drop.add(nclose_key)
    return nclose_to_drop

nclose_to_drop_predict = find_nclose_to_drop(predict_suc_B, cf_nclose_depth)
nclose_to_drop_observed = find_nclose_to_drop(B, cf_nclose_depth)
nclose_to_drop = nclose_to_drop_predict | nclose_to_drop_observed

logging.info(
    f'BP step/depth filter : nclose to remove = {len(nclose_to_drop)} '
    f'(Source: predict={len(nclose_to_drop_predict)}, depth={len(nclose_to_drop_observed)})'
)

drop_cols = set()
for col_idx, ctr in enumerate(path_nclose_count):
    if any(nc in nclose_to_drop for nc in ctr):
        drop_cols.add(col_idx)

if drop_cols:
    new_A_idx_list = sorted(c for c in A_idx_list if c not in drop_cols)
    with h5py.File(f"{PREFIX}/matrix.h5", "r") as f:
        dA = f["A"]
        dAf = f["A_fail"]
        read_selected_feature_rows_direct(dA, A_store, new_A_idx_list)
        read_selected_feature_rows_direct(dAf, A_fail_store, new_A_idx_list)
    
    A_new = solver_matrix_from_store(A_store, len(new_A_idx_list))
    A_fail_new = solver_matrix_from_store(A_fail_store, len(new_A_idx_list))
    bp_warm_fullsize = scatter_subset_weights(
        filter_weights, A_idx_list, len(weights_fullsize), dtype=weights_fullsize.dtype
    )
    filter_weights = fit_nnls_warm(
        A_new, B, warm_start=bp_warm_fullsize[new_A_idx_list]
    )
    predict_suc_B = A_new.dot(filter_weights)
    predict_fal_B = A_fail_new.dot(filter_weights)
    A_idx_list = new_A_idx_list
else:
    logging.info('BP step/depth filter : no nclose removed, skip refit')

final_weights_fullsize = np.zeros_like(weights_fullsize)
final_weights_fullsize[A_idx_list] = filter_weights

b_norm = np.linalg.norm(B)

error = np.linalg.norm(predict_suc_B - B)
predict_B = np.concatenate((predict_suc_B, predict_fal_B))[b_start_ind:]

logging.info(f'Error : {error:.4f}')
logging.info(f'Relative error : {error/b_norm:.4f}')

np.save(f'{PREFIX}/weight.npy', final_weights_fullsize)
np.save(f'{PREFIX}/predict_B.npy', predict_B)

with open(f'{PREFIX}/A_idx_list.pkl', 'wb') as f:
    pkl.dump(A_idx_list, f)
