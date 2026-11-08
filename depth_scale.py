"""Depth per copy (1N) written by stage 00 and read by the stages that show copy numbers.

Two units are kept apart:
  threshold unit  median depth / 2 over all non-chrM rows of the normalized depth. Every
                  cut that only asks whether a weight or depth difference is small (0.1N,
                  0.2N, 5% of the median, ...) keeps this unit.
  copy unit       the 1N of copy numbers shown to the user: Virtual SKY labels,
                  karyotype.txt, Circos levels and colours, and the CN values of the VCF,
                  BED and TSV reports. Stage 00 stores it in depth_scale.json: the
                  B-allele-fraction 1N of baf_ploidy.py whenever that returned one (whatever
                  its status), otherwise median depth / 2.
A stage falls back to median depth / 2 when the json is missing or unreadable, or was written
for another depth table (row count, median or sum of the depths differ from the table the stage
read). Stage 00 removes the json before it writes a new depth table.
"""
import gzip
import json
import logging
import math
import os

import numpy as np

JSON_NAME = 'depth_scale.json'
DETAIL_DIR = 'depth_scale'


def finite_or_none(x):
    """JSON-safe copy: non-finite floats become null, numpy scalars and arrays plain values."""
    if isinstance(x, dict):
        return {str(k): finite_or_none(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [finite_or_none(v) for v in x]
    if isinstance(x, np.ndarray):
        return finite_or_none(x.tolist())
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (float, np.floating)):
        return float(x) if math.isfinite(x) else None
    if isinstance(x, np.integer):
        return int(x)
    return x


def fingerprint(depth):
    """Row count, median and exact sum of a depth table's non-chrM mean depths."""
    depth = np.asarray(depth, dtype=float)
    return {'depth_rows': int(depth.size), 'median_depth': float(np.median(depth)),
            'depth_sum': math.fsum(depth)}


def read_table(path):
    """{(chrom, start): mean depth} of the non-chrM rows of a win.stat.gz table, in file order."""
    table = {}
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if line.startswith('#') or len(fields) < 8 or fields[0] == 'chrM':
                continue
            table[(fields[0], int(fields[1]))] = float(fields[7])
    return table


def remove(prefix):
    path = os.path.join(prefix, JSON_NAME)
    if os.path.lexists(path):
        os.remove(path)


def write(prefix, depth_file, depth, gc_info, baf_info):
    """Write depth_scale.json for the normalized depth (all non-chrM rows) and return it."""
    table = fingerprint(depth)
    median = table['median_depth']
    n = baf_info.get('N')
    use_baf = n is not None and math.isfinite(n) and n > 0
    scale = {
        'depth_file': os.path.abspath(depth_file),
        **table,
        'median_half': median / 2,
        'copy_unit': float(n) if use_baf else median / 2,
        'copy_unit_source': 'baf' if use_baf else 'median',
        'gc_correction': gc_info,
        'baf': baf_info,
    }
    path = os.path.join(prefix, JSON_NAME)
    tmp = f'{path}.tmp'
    with open(tmp, 'w') as fh:
        json.dump(finite_or_none(scale), fh, indent=1, allow_nan=False)
    os.replace(tmp, path)
    return scale


def _matches(scale, depth):
    now = fingerprint(depth)
    return (int(scale['depth_rows']) == now['depth_rows']
            and math.isclose(float(scale['median_depth']), now['median_depth'], rel_tol=1e-9)
            and math.isclose(float(scale['depth_sum']), now['depth_sum'], rel_tol=1e-9))


def _load(path):
    with open(path) as fh:
        scale = json.load(fh)
    float(scale['copy_unit'])
    return scale


def _baf_unit(path, scale):
    if scale.get('copy_unit_source') != 'baf':
        return None
    unit = float(scale['copy_unit'])
    if not (math.isfinite(unit) and unit > 0):
        logging.warning(f'{path} holds an invalid 1N ({unit}); copy numbers use median depth / 2')
        return None
    return unit


def copy_unit(prefix, depth):
    """(1N for shown copy numbers, 'baf' or 'median') for a stage that read these normalized
    depths (all non-chrM rows)."""
    half = float(np.median(depth)) / 2
    path = os.path.join(prefix, JSON_NAME)
    if not os.path.isfile(path):
        return half, 'median'
    try:
        scale = _load(path)
        same = _matches(scale, depth)
    except (OSError, ValueError, KeyError, TypeError) as exc:
        logging.warning(f'{path} is unreadable ({exc}); copy numbers use median depth / 2')
        return half, 'median'
    if not same:
        logging.warning(f'{path} was written for another depth table; copy numbers use median depth / 2')
        return half, 'median'
    unit = _baf_unit(path, scale)
    return (half, 'median') if unit is None else (unit, 'baf')


def _shows_table(table, coords, depth):
    """True when `depth` at `coords` holds the table's depths as stage 22 stores them: float32,
    within a few float32 steps (pandas and Python parse a few decimals one bit apart), and on
    chrY either that depth or 0 (stage 22 zeroes chrY when the sample has none)."""
    keys = [(str(c), int(s)) for c, s in coords]
    if len(keys) != len(depth) or any(k not in table for k in keys):
        return False
    want = np.asarray([table[k] for k in keys], dtype=np.float32)
    got = np.asarray(depth, dtype=np.float32)
    same = np.isclose(got, want, rtol=1e-6, atol=0)
    zeroed = np.array([c == 'chrY' for c, _ in keys], dtype=bool) & (got == 0)
    return bool((same | zeroed).all())


def recorded_copy_unit(prefix, coords, depth):
    """BAF 1N of depth_scale.json for a reader without the depth table (the viewer), or None.
    The table named in the json must still match it, and `depth` at `coords` must be the
    depths of that table."""
    path = os.path.join(prefix, JSON_NAME)
    if not os.path.isfile(path):
        return None
    try:
        scale = _load(path)
        table = read_table(scale['depth_file'])
        same = _matches(scale, np.fromiter(table.values(), dtype=float, count=len(table)))
    except (OSError, ValueError, KeyError, TypeError, EOFError) as exc:
        logging.warning(f'{path} cannot be checked ({exc}); copy numbers use median depth / 2')
        return None
    if not same:
        logging.warning(f'{path} no longer matches {scale["depth_file"]}; copy numbers use median depth / 2')
        return None
    if not _shows_table(table, coords, depth):
        logging.warning(f'the depths shown are not those of {scale["depth_file"]}; '
                        'copy numbers use median depth / 2')
        return None
    return _baf_unit(path, scale)
