import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from skype_utils import *

import argparse
import copy
import gzip
import logging
import pickle as pkl

import numpy as np

from juliacall import Main as jl

# logging 설정(레벨/포맷)은 skype_utils 에서 중앙 관리한다 (LOG_LEVEL).
logging.info("03_anal_bam start")

RAW_TRANSLOCATION_CANDIDATE_PKL = 'raw_translocation_candidates.pkl'
RAW_TRANSLOCATION_RESULT_PKL = 'raw_translocation_result.pkl'
RAW_TRANSLOCATION_REPORT_TSV = 'raw_translocation_read_counts.tsv'
NCLOSE_COUNT_CANDIDATE_PKL = 'nclose_count_candidates.pkl'
NCLOSE_COUNT_RESULT_PKL = 'nclose_count_result.pkl'
NCLOSE_COUNT_REPORT_TSV = 'nclose_read_counts.tsv'
RAW_TRANSLOCATION_WINDOW = 5000
RAW_TRANSLOCATION_DEPTH_FLANK = 500000
RAW_TRANSLOCATION_DEPTH_BALANCED_RATIO = 1.1
NCLOSE_COUNT_DEFAULT_VAF_THRESHOLD = 0.1
RAW_COUNT_NAMES = {
    'd1': 'nclose_a_read_count',
    'd2': 'point_a_norm_read_count',
    'd3': 'point_b_norm_read_count',
    'd4': 'nclose_b_read_count',
}
RAW_BASIS_NAMES = {
    'd1': 'nclose_a_basis',
    'd2': 'point_a_norm_basis',
    'd3': 'point_b_norm_basis',
    'd4': 'nclose_b_basis',
}
NCLOSE_COUNT_INDEX = {
    'd1': 101,
    'd2': 102,
    'd3': 103,
}

parser = argparse.ArgumentParser(description="SKYPE raw-read translocation analysis")

parser.add_argument("prefix",
                    help="Prefix for pipeline")

parser.add_argument("read_bam_loc",
                    help="Raw read alignment bam location")

parser.add_argument("reference_fai_path",
                    help="Path to the chromosome information file.")

parser.add_argument("depth_stat_path", nargs="?", default=None,
                    help="Optional normalized depth win.stat.gz for 500kb coordinate depth reporting.")

parser.add_argument("--progress",
                    help="Show progress bar", action='store_true')
parser.add_argument("--nclose_count_only",
                    help="Only run single-nclose raw-read VAF counting.",
                    action='store_true')
parser.add_argument("--skip_nclose_count",
                    help="Skip single-nclose raw-read VAF counting even if candidates exist.",
                    action='store_true')
parser.add_argument("--nclose_count_vaf_threshold",
                    help="VAF threshold used for the pass/fail column in nclose_read_counts.tsv.",
                    type=float, default=NCLOSE_COUNT_DEFAULT_VAF_THRESHOLD)

args = parser.parse_args()

PREFIX = args.prefix
read_bam_loc = args.read_bam_loc
reference_fai_path = args.reference_fai_path
depth_stat_path = args.depth_stat_path

raw_translocation_candidates = []
if not args.nclose_count_only:
    raw_candidate_path = f"{PREFIX}/{RAW_TRANSLOCATION_CANDIDATE_PKL}"
    if os.path.isfile(raw_candidate_path):
        with open(raw_candidate_path, "rb") as f:
            raw_translocation_candidates = pkl.load(f)

nclose_count_candidates = []
if not args.skip_nclose_count:
    nclose_count_candidate_path = f"{PREFIX}/{NCLOSE_COUNT_CANDIDATE_PKL}"
    if os.path.isfile(nclose_count_candidate_path):
        with open(nclose_count_candidate_path, "rb") as f:
            nclose_count_candidates = pkl.load(f)

def read_chr_len(fai_path):
    chr_len = {}
    with open(fai_path, "r") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                chr_len[parts[0]] = int(parts[1])
    return chr_len

def read_depth_by_chrom(depth_path):
    if depth_path is None or not os.path.isfile(depth_path):
        return {}

    depth_by_chrom = {}
    opener = gzip.open if str(depth_path).endswith(".gz") else open
    with opener(depth_path, "rt") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            try:
                chrom = parts[0]
                st = int(parts[1])
                nd = int(parts[2])
                meandepth = float(parts[7])
            except ValueError:
                continue
            depth_by_chrom.setdefault(chrom, []).append((st, nd, meandepth))

    return depth_by_chrom

def mean_depth_around_coord(depth_by_chrom, chrom, coord, chr_len,
                            flank=RAW_TRANSLOCATION_DEPTH_FLANK):
    if coord == '*' or chrom not in depth_by_chrom:
        return None

    chrom_len = int(chr_len.get(chrom, int(coord) + flank))
    st = max(1, int(coord) - flank)
    nd = min(chrom_len, int(coord) + flank)
    if nd < st:
        return None

    weighted_sum = 0.0
    weight_sum = 0
    for win_st, win_nd, meandepth in depth_by_chrom[chrom]:
        overlap_st = max(st, win_st)
        overlap_nd = min(nd, win_nd)
        if overlap_nd < overlap_st:
            continue
        weight = overlap_nd - overlap_st + 1
        weighted_sum += meandepth * weight
        weight_sum += weight

    if weight_sum == 0:
        return None
    return weighted_sum / weight_sum

def format_depth(value):
    return '*' if value is None or not np.isfinite(value) else round(float(value), 6)

def depth_pair_is_balanced(depth_pair, ratio=RAW_TRANSLOCATION_DEPTH_BALANCED_RATIO):
    if not depth_pair:
        return False
    front = depth_pair.get('front')
    back = depth_pair.get('back')
    if front is None or back is None:
        return False
    if not (np.isfinite(front) and np.isfinite(back)):
        return False
    if front < 0 or back < 0:
        return False
    low, high = sorted([float(front), float(back)])
    if low == 0:
        return high == 0
    return (high / low) <= ratio

def depth_pair_mean(depth_pair):
    if not depth_pair:
        return None
    front = depth_pair.get('front')
    back = depth_pair.get('back')
    if front is None or back is None:
        return None
    if not (np.isfinite(front) and np.isfinite(back)):
        return None
    return (float(front) + float(back)) / 2

def estimate_nclose_depth(point_depth, nclose_count, point_count):
    weight = int(nclose_count) + int(point_count)
    if point_depth is None or weight <= 0:
        return None, None, 0
    vaf = int(nclose_count) / weight
    return vaf, float(point_depth) * vaf, weight

def weighted_expected_depth(*estimate_data):
    weighted_sum = 0.0
    weight_sum = 0
    for _, expected_depth, weight in estimate_data:
        if expected_depth is None or not np.isfinite(expected_depth) or weight <= 0:
            continue
        weighted_sum += expected_depth * weight
        weight_sum += weight
    if weight_sum == 0:
        return None
    return weighted_sum / weight_sum

def clamp_break_coord(coord, chrom_len):
    return max(1, min(int(coord), int(chrom_len)))

def half_open_to_1based_interval(chrom, st0, nd0, chr_len):
    chrom_len = int(chr_len.get(chrom, max(int(st0), int(nd0), 1)))
    st0 = max(0, min(int(st0), chrom_len))
    nd0 = max(0, min(int(nd0), chrom_len))
    if nd0 < st0:
        st0, nd0 = nd0, st0
    if nd0 == st0:
        if nd0 < chrom_len:
            nd0 += 1
        elif st0 > 0:
            st0 -= 1
        else:
            return None
    st1 = st0 + 1
    nd1 = nd0
    if nd1 < st1:
        return None
    return st1, nd1

def endpoint_span_row(pair_id, count_idx, anchor_idx, chrom, coord0,
                      expected_positive, role, flank, chr_len):
    flank = max(0, int(flank))
    coord0 = int(coord0)
    lower_reference_side = (
        (role == 1 and expected_positive) or
        (role == 2 and not expected_positive)
    )
    if lower_reference_side:
        st0 = coord0 - flank if flank > 0 else coord0 - 1
        nd0 = coord0
    else:
        st0 = coord0
        nd0 = coord0 + flank if flank > 0 else coord0 + 1

    interval = half_open_to_1based_interval(chrom, st0, nd0, chr_len)
    if interval is None:
        return None
    point_coord = interval[0] if lower_reference_side else interval[1]
    chrom_len = int(chr_len.get(chrom, max(interval[1], point_coord, 1)))
    return [
        int(pair_id), int(count_idx), int(anchor_idx), chrom,
        int(interval[0]), int(interval[1]),
        clamp_break_coord(point_coord, chrom_len),
        bool(expected_positive),
    ]

def breakend_query_gap(endpoints):
    if len(endpoints) != 2:
        return None
    e1, e2 = endpoints
    if e1.get('ctg_name') != e2.get('ctg_name'):
        return None
    try:
        gap = int(e2['ctg_st']) - int(e1['ctg_nd'])
    except KeyError:
        return None
    return max(0, gap)

def add_breakend_span_rows(span_rows, pair_id, split, chr_len):
    endpoints = split.get('endpoints', [])
    if len(endpoints) != 2:
        return
    query_gap = breakend_query_gap(endpoints)
    if query_gap is None or query_gap > RAW_TRANSLOCATION_WINDOW:
        return

    remaining = RAW_TRANSLOCATION_WINDOW - query_gap
    flank1 = remaining // 2
    flank2 = remaining - flank1
    e1, e2 = endpoints
    row1 = endpoint_span_row(
        pair_id, int(split['count_idx']), 1,
        e1['chrom'], int(e1['coord']), e1['dir'] == '+',
        1, flank1, chr_len,
    )
    row2 = endpoint_span_row(
        pair_id, int(split['count_idx']), 2,
        e2['chrom'], int(e2['coord']), e2['dir'] == '+',
        2, flank2, chr_len,
    )
    if row1 is not None and row2 is not None:
        span_rows.extend([row1, row2])

def add_reference_span_rows(span_rows, pair_id, count_idx, side, chr_len):
    chrom = side['chrom']
    coord1 = int(side['inner_st'])
    coord2 = int(side['inner_nd'])
    if coord2 < coord1:
        coord1, coord2 = coord2, coord1
    ref_gap = coord2 - coord1
    if ref_gap > RAW_TRANSLOCATION_WINDOW:
        return

    remaining = RAW_TRANSLOCATION_WINDOW - ref_gap
    flank1 = remaining // 2
    flank2 = remaining - flank1
    row1 = endpoint_span_row(
        pair_id, count_idx, 1, chrom, coord1, True, 1, flank1, chr_len,
    )
    row2 = endpoint_span_row(
        pair_id, count_idx, 2, chrom, coord2, True, 2, flank2, chr_len,
    )
    if row1 is not None and row2 is not None:
        span_rows.extend([row1, row2])

def endpoint_segment_text(endpoint):
    return repr([
        endpoint['ctg_name'],
        endpoint['dir'],
        int(endpoint['ref_st']),
        int(endpoint['ref_nd']),
    ])

def layout_segment_text(layout):
    return " ".join(endpoint_segment_text(endpoint) for endpoint in layout['ordered_endpoints'])

def directed_entry_coord(endpoint):
    return int(endpoint['ref_st']) if endpoint['dir'] == '+' else int(endpoint['ref_nd'])

def candidate_display_points(record):
    chrom_pair = record.get('chrom_pair', ('*', '*'))
    layout_a = record.get('layout_a', {})
    layout_b = record.get('layout_b', {})
    endpoints_a = list(layout_a.get('endpoints', ()))
    endpoints_b = list(layout_b.get('endpoints', ()))
    side_records = record.get('side_records', [])

    chrom_a = chrom_pair[0] if len(chrom_pair) > 0 else '*'
    chrom_b = chrom_pair[1] if len(chrom_pair) > 1 else '*'
    coord_a = '*'
    coord_b = '*'

    if len(endpoints_b) > 0:
        coord_a = directed_entry_coord(endpoints_b[0])
    elif len(side_records) > 0:
        coord_a = int(side_records[0]['inner_st'])

    if len(endpoints_a) > 1:
        coord_b = directed_entry_coord(endpoints_a[1])
    elif len(side_records) > 1:
        coord_b = int(side_records[1]['inner_nd'])

    return chrom_a, coord_a, chrom_b, coord_b

def span_row_basis_text(row):
    chrom = row[3]
    point_coord = int(row[6])
    direction = '+' if bool(row[7]) else '-'
    return f"{chrom}:{point_coord}({direction})"

def build_count_basis_texts(span_rows):
    count_names = {1: 'd1', 2: 'd2', 3: 'd3', 4: 'd4'}
    rows_by_count = {}
    for row in span_rows:
        pair_id = int(row[0])
        count_idx = int(row[1])
        if count_idx not in count_names:
            continue
        rows_by_count.setdefault((pair_id, count_idx), []).append(row)

    basis_by_pair = {}
    for (pair_id, count_idx), rows in rows_by_count.items():
        rows.sort(key=lambda row: int(row[2]))
        basis_by_pair.setdefault(pair_id, {})[count_names[count_idx]] = "->".join(
            span_row_basis_text(row) for row in rows
        )
    return basis_by_pair

def build_nclose_count_basis_texts(span_rows):
    count_names = {v: k for k, v in NCLOSE_COUNT_INDEX.items()}
    rows_by_count = {}
    for row in span_rows:
        pair_id = int(row[0])
        count_idx = int(row[1])
        if count_idx not in count_names:
            continue
        rows_by_count.setdefault((pair_id, count_idx), []).append(row)

    basis_by_pair = {}
    for (pair_id, count_idx), rows in rows_by_count.items():
        rows.sort(key=lambda row: int(row[2]))
        basis_by_pair.setdefault(pair_id, {})[count_names[count_idx]] = "->".join(
            span_row_basis_text(row) for row in rows
        )
    return basis_by_pair

def build_point_depths(span_rows, depth_by_chrom, chr_len):
    count_names = {2: 'point_a', 3: 'point_b'}
    rows_by_count = {}
    for row in span_rows:
        pair_id = int(row[0])
        count_idx = int(row[1])
        if count_idx not in count_names:
            continue
        rows_by_count.setdefault((pair_id, count_idx), []).append(row)

    point_depths = {}
    for (pair_id, count_idx), rows in rows_by_count.items():
        rows.sort(key=lambda row: int(row[2]))
        depths = []
        for row in rows:
            chrom = row[3]
            point_coord = int(row[6])
            depths.append(mean_depth_around_coord(depth_by_chrom, chrom, point_coord, chr_len))
        while len(depths) < 2:
            depths.append(None)
        point_depths.setdefault(pair_id, {})[count_names[count_idx]] = {
            'front': depths[0],
            'back': depths[1],
        }
    return point_depths

def nclose_query_ordered_endpoints(record):
    endpoints = list(record.get('layout', {}).get('ordered_endpoints', ()))
    if len(endpoints) != 2:
        return endpoints
    if endpoints[0].get('ctg_name') != endpoints[1].get('ctg_name'):
        return endpoints
    try:
        return sorted(endpoints, key=lambda endpoint: (
            min(int(endpoint['ctg_st']), int(endpoint['ctg_nd'])),
            max(int(endpoint['ctg_st']), int(endpoint['ctg_nd'])),
        ))
    except KeyError:
        return endpoints

def nclose_candidate_query_gap(record):
    endpoints = nclose_query_ordered_endpoints(record)
    if len(endpoints) != 2:
        return None
    if endpoints[0].get('ctg_name') != endpoints[1].get('ctg_name'):
        return None
    try:
        left_nd = max(int(endpoints[0]['ctg_st']), int(endpoints[0]['ctg_nd']))
        right_st = min(int(endpoints[1]['ctg_st']), int(endpoints[1]['ctg_nd']))
    except KeyError:
        return None
    return max(0, right_st - left_nd)

def add_nclose_count_rows(span_rows, candidate, chr_len):
    query_gap = nclose_candidate_query_gap(candidate)
    if query_gap is None or query_gap > RAW_TRANSLOCATION_WINDOW:
        return

    endpoints = nclose_query_ordered_endpoints(candidate)
    if len(endpoints) != 2:
        return

    pair_id = int(candidate['pair_id'])
    add_breakend_span_rows(
        span_rows, pair_id,
        {'count_idx': NCLOSE_COUNT_INDEX['d2'], 'endpoints': endpoints},
        chr_len,
    )
    add_reference_span_rows(
        span_rows, pair_id, NCLOSE_COUNT_INDEX['d1'],
        {
            'chrom': endpoints[0]['chrom'],
            'inner_st': int(endpoints[0]['coord']),
            'inner_nd': int(endpoints[0]['coord']),
        },
        chr_len,
    )
    add_reference_span_rows(
        span_rows, pair_id, NCLOSE_COUNT_INDEX['d3'],
        {
            'chrom': endpoints[1]['chrom'],
            'inner_st': int(endpoints[1]['coord']),
            'inner_nd': int(endpoints[1]['coord']),
        },
        chr_len,
    )

def nclose_candidate_display_points(record):
    endpoints = nclose_query_ordered_endpoints(record)
    chrom_a = endpoints[0]['chrom'] if len(endpoints) > 0 else '*'
    chrom_b = endpoints[1]['chrom'] if len(endpoints) > 1 else '*'
    coord_a = int(endpoints[0]['coord']) if len(endpoints) > 0 else '*'
    coord_b = int(endpoints[1]['coord']) if len(endpoints) > 1 else '*'
    return chrom_a, coord_a, chrom_b, coord_b

def nclose_count_vaf(nclose_count, normal_count):
    total = int(nclose_count) + int(normal_count)
    if total <= 0:
        return None
    return int(nclose_count) / total

chr_len = read_chr_len(reference_fai_path)
depth_by_chrom = read_depth_by_chrom(depth_stat_path)
raw_junction_span_rows = []
for candidate in raw_translocation_candidates:
    pair_id = int(candidate['pair_id'])
    for split in candidate.get('split_junctions', []):
        add_breakend_span_rows(raw_junction_span_rows, pair_id, split, chr_len)
    side_by_id = {side['side_id']: side for side in candidate.get('side_records', [])}
    for span in candidate.get('span_junctions', []):
        side = side_by_id.get(span['side_id'])
        if side is None:
            continue
        add_reference_span_rows(raw_junction_span_rows, pair_id, int(span['count_idx']), side, chr_len)
raw_count_basis_texts = build_count_basis_texts(raw_junction_span_rows)
raw_point_depths = build_point_depths(raw_junction_span_rows, depth_by_chrom, chr_len)

nclose_count_span_rows = []
for candidate in nclose_count_candidates:
    add_nclose_count_rows(nclose_count_span_rows, candidate, chr_len)
nclose_count_basis_texts = build_nclose_count_basis_texts(nclose_count_span_rows)

all_junction_span_rows = raw_junction_span_rows + nclose_count_span_rows

jl.seval("global raw_junction_span_vec = Vector{Vector{Any}}()")
for row in all_junction_span_rows:
    jl.push_b(jl.raw_junction_span_vec, jl.Vector[jl.Any](row))

is_progress_bar = sys.stdout.isatty() or args.progress
jl.include(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'anal_bam.jl'))

raw_translocation_count_list = jl.anal_bam(
    read_bam_loc, reference_fai_path, is_progress_bar,
    jl.raw_junction_span_vec
)

raw_translocation_counts = {
    (int(pair_id), int(count_idx)): int(count)
    for pair_id, count_idx, count in raw_translocation_count_list
}

raw_translocation_records = []
if not args.nclose_count_only:
    for candidate in raw_translocation_candidates:
        pair_id = int(candidate['pair_id'])
        d1 = raw_translocation_counts.get((pair_id, 1), 0)
        d2 = raw_translocation_counts.get((pair_id, 2), 0)
        d3 = raw_translocation_counts.get((pair_id, 3), 0)
        d4 = raw_translocation_counts.get((pair_id, 4), 0)
        point_a_no_spanning = d2 == 0
        point_b_no_spanning = d3 == 0
        point_depth = raw_point_depths.get(pair_id, {})
        point_a_depth_balanced = depth_pair_is_balanced(point_depth.get('point_a', {}))
        point_b_depth_balanced = depth_pair_is_balanced(point_depth.get('point_b', {}))
        depth_balanced_translocation = point_a_depth_balanced and point_b_depth_balanced
        point_a_depth_mean = depth_pair_mean(point_depth.get('point_a', {}))
        point_b_depth_mean = depth_pair_mean(point_depth.get('point_b', {}))
        point_a_estimate = estimate_nclose_depth(point_a_depth_mean, d1, d2)
        point_b_estimate = estimate_nclose_depth(point_b_depth_mean, d4, d3)
        weighted_nclose_depth = weighted_expected_depth(point_a_estimate, point_b_estimate)

        record = copy.deepcopy(candidate)
        record['read_counts'] = {'d1': d1, 'd2': d2, 'd3': d3, 'd4': d4}
        record['raw_point_no_spanning'] = {'point_a': point_a_no_spanning, 'point_b': point_b_no_spanning}
        record['point_500k_depth'] = point_depth
        record['raw_point_depth_balanced'] = {
            'point_a': point_a_depth_balanced,
            'point_b': point_b_depth_balanced,
        }
        record['depth_balanced_translocation'] = depth_balanced_translocation
        record['depth_weighted_nclose_estimate'] = {
            'point_a_mean_500k_depth': point_a_depth_mean,
            'point_b_mean_500k_depth': point_b_depth_mean,
            'point_a_nclose_vaf': point_a_estimate[0],
            'point_b_nclose_vaf': point_b_estimate[0],
            'nclose_a_expected_depth': point_a_estimate[1],
            'nclose_b_expected_depth': point_b_estimate[1],
            'weighted_expected_nclose_depth': weighted_nclose_depth,
        }
        for side_idx, side in enumerate(record.get('side_records', [])):
            point_no_spanning = point_a_no_spanning if side_idx == 0 else point_b_no_spanning
            point_depth_balanced = point_a_depth_balanced if side_idx == 0 else point_b_depth_balanced
            side['no_spanning_rawread'] = point_no_spanning
            side['no_spanning_utg'] = point_no_spanning
            side['raw_point_no_spanning'] = point_no_spanning
            side['raw_depth_balanced_500k'] = point_depth_balanced
            side['depth_balanced_translocation'] = depth_balanced_translocation
            side['accepted_forbid'] = False
            side.setdefault('crossing_cols', [])
        raw_translocation_records.append(record)

    with open(f"{PREFIX}/{RAW_TRANSLOCATION_RESULT_PKL}", "wb") as f:
        pkl.dump(raw_translocation_records, f)

    with open(f"{PREFIX}/{RAW_TRANSLOCATION_REPORT_TSV}", "wt") as f:
        print(*[
            'chrom_a', 'point_a', 'chrom_b', 'point_b',
            'point_a_front_500k_mean_depth', 'point_a_back_500k_mean_depth',
            'point_b_front_500k_mean_depth', 'point_b_back_500k_mean_depth',
            'nclose_a', 'nclose_b',
            RAW_COUNT_NAMES['d1'], RAW_COUNT_NAMES['d2'],
            RAW_COUNT_NAMES['d4'], RAW_COUNT_NAMES['d3'],
            'point_a_no_spanning_rawread', 'point_b_no_spanning_rawread',
            'point_a_nclose_vaf', 'point_b_nclose_vaf',
            'nclose_a_expected_depth', 'nclose_b_expected_depth',
            'weighted_expected_nclose_depth',
            RAW_BASIS_NAMES['d1'], RAW_BASIS_NAMES['d2'],
            RAW_BASIS_NAMES['d4'], RAW_BASIS_NAMES['d3'],
            'pair_id',
        ], sep='\t', file=f)
        for record in raw_translocation_records:
            if not record.get('depth_balanced_translocation', False):
                continue
            counts = record['read_counts']
            chrom_a, coord_a, chrom_b, coord_b = candidate_display_points(record)
            basis = raw_count_basis_texts.get(int(record['pair_id']), {})
            point_depth = record.get('point_500k_depth', {})
            point_a_depth = point_depth.get('point_a', {})
            point_b_depth = point_depth.get('point_b', {})
            depth_estimate = record.get('depth_weighted_nclose_estimate', {})
            print(
                chrom_a, coord_a, chrom_b, coord_b,
                format_depth(point_a_depth.get('front')),
                format_depth(point_a_depth.get('back')),
                format_depth(point_b_depth.get('front')),
                format_depth(point_b_depth.get('back')),
                layout_segment_text(record['layout_a']),
                layout_segment_text(record['layout_b']),
                counts['d1'], counts['d2'],
                counts['d4'], counts['d3'],
                record.get('raw_point_no_spanning', {}).get('point_a', False),
                record.get('raw_point_no_spanning', {}).get('point_b', False),
                format_depth(depth_estimate.get('point_a_nclose_vaf')),
                format_depth(depth_estimate.get('point_b_nclose_vaf')),
                format_depth(depth_estimate.get('nclose_a_expected_depth')),
                format_depth(depth_estimate.get('nclose_b_expected_depth')),
                format_depth(depth_estimate.get('weighted_expected_nclose_depth')),
                basis.get('d1', '*'), basis.get('d2', '*'),
                basis.get('d4', '*'), basis.get('d3', '*'),
                record['pair_id'],
                sep='\t', file=f
            )

    logging.info(
        f"Raw-read translocation point no-span sides : "
        f"{sum(sum(1 for s in r.get('side_records', []) if s.get('no_spanning_rawread', False)) for r in raw_translocation_records)}"
    )

nclose_count_records = []
if not args.skip_nclose_count:
    for candidate in nclose_count_candidates:
        pair_id = int(candidate['pair_id'])
        d1 = raw_translocation_counts.get((pair_id, NCLOSE_COUNT_INDEX['d1']), 0)
        d2 = raw_translocation_counts.get((pair_id, NCLOSE_COUNT_INDEX['d2']), 0)
        d3 = raw_translocation_counts.get((pair_id, NCLOSE_COUNT_INDEX['d3']), 0)
        chr_a_vaf = nclose_count_vaf(d2, d1)
        chr_b_vaf = nclose_count_vaf(d2, d3)
        query_gap = nclose_candidate_query_gap(candidate)
        filter_eligible = query_gap is not None and query_gap <= RAW_TRANSLOCATION_WINDOW
        keep_by_vaf = filter_eligible and (
            (chr_a_vaf is not None and chr_a_vaf >= args.nclose_count_vaf_threshold) or
            (chr_b_vaf is not None and chr_b_vaf >= args.nclose_count_vaf_threshold)
        )
        keep_nclose = (not filter_eligible) or keep_by_vaf

        record = copy.deepcopy(candidate)
        record['read_counts'] = {'d1': d1, 'd2': d2, 'd3': d3}
        record['nclose_query_gap'] = query_gap
        record['filter_eligible'] = filter_eligible
        record['vaf'] = {'chr_a': chr_a_vaf, 'chr_b': chr_b_vaf}
        record['keep_by_vaf'] = keep_by_vaf
        record['keep_nclose'] = keep_nclose
        nclose_count_records.append(record)

    with open(f"{PREFIX}/{NCLOSE_COUNT_RESULT_PKL}", "wb") as f:
        pkl.dump(nclose_count_records, f)

    with open(f"{PREFIX}/{NCLOSE_COUNT_REPORT_TSV}", "wt") as f:
        print(*[
            'chrom_a', 'point_a', 'chrom_b', 'point_b',
            'nclose_query_gap', 'filter_eligible',
            'nclose',
            'chr_a_norm_read_count', 'nclose_read_count', 'chr_b_norm_read_count',
            'chr_a_nclose_vaf', 'chr_b_nclose_vaf', 'keep_by_vaf', 'keep_nclose',
            'chr_a_norm_basis', 'nclose_basis', 'chr_b_norm_basis',
            'pair_id', 'nclose_key', 'contig_name',
        ], sep='\t', file=f)
        for record in nclose_count_records:
            counts = record['read_counts']
            basis = nclose_count_basis_texts.get(int(record['pair_id']), {})
            chrom_a, coord_a, chrom_b, coord_b = nclose_candidate_display_points(record)
            print(
                chrom_a, coord_a, chrom_b, coord_b,
                record.get('nclose_query_gap', '*'),
                record.get('filter_eligible', False),
                layout_segment_text(record['layout']),
                counts['d1'], counts['d2'], counts['d3'],
                format_depth(record.get('vaf', {}).get('chr_a')),
                format_depth(record.get('vaf', {}).get('chr_b')),
                record.get('keep_by_vaf', False),
                record.get('keep_nclose', False),
                basis.get('d1', '*'), basis.get('d2', '*'), basis.get('d3', '*'),
                record['pair_id'],
                ",".join(str(x) for x in record.get('nclose_key', ())),
                record.get('ctg_name', '*'),
                sep='\t', file=f
            )

    logging.info(
        f"NClose raw-count VAF kept : "
        f"{sum(1 for r in nclose_count_records if r.get('keep_nclose', False))}/"
        f"{len(nclose_count_records)}"
    )
