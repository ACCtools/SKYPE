"""Per-sample GC correction of the stage-00 normalized window depth.

Stage 00 divides the sample depth by a CHM13 reference profile, which removes the GC
response the two libraries share but not the part that differs between libraries. This
module fits the remaining response on the sample's own depth and the GC of each window and
divides it out. Nothing from a later stage (no fitted karyotype) is used.

Rows
  exempt     chrY and rows overlapping a CEN-SAT interval > 1 Mb (00_depth_norm.is_in_repeat).
             Stage 00 does not normalize them, so they are left unchanged here too.
  training   the stage-22 row rule on the depth before GC correction (not fully inside a
             CEN-SAT interval, depth <= 3 x median of all rows), not exempt, depth > 0 and GC
             known (N fraction of the window <= 0.5).
  corrected  every non-exempt row with known GC.
Curve: log2 depth = segment level + f(GC). Segments come from an exact optimal partition of
the GC-corrected log2 depth per chromosome; after each fit the corrected depth is segmented
again and new breakpoints are added (at most 10 rounds). Within fixed segments f is fitted by
damped backfitting: f = median of (log2 depth - segment level) per 0.01 GC bin.
Safeguards
  support    a GC bin informs f only with >= 30 training rows from >= 10 segments; sparser GC
             bins are pooled with their neighbours until the pool has that support. Beyond
             the outermost supported GC f is flat (no extrapolation).
  shrinkage  f is refitted on each half of the chromosomes. lambda = 1 - Var(fA - fB) /
             (4 Var(f)) over the training rows, clipped to [0, 1]; the applied curve is
             lambda * f. Var(fA - fB) / 4 estimates the error variance of the full-data
             curve, so lambda is its reliability.
  fallback   no correction (identity) with a reason when the depth is not finite, there are
             < 2000 training rows, < 8 supported GC bins, a fit fails or does not converge,
             or lambda = 0.
Level: the corrected rows are scaled by one common factor so the median over all rows, which
stage 22 (3 x median cut) and later stages (1N = median / 2 for thresholds) take, is kept.
"""
import gzip

import numpy as np
import pandas as pd

GC_EDGES = np.arange(0.28, 0.66, 0.01)
MIN_GC_BIN_ROWS = 30
MIN_GC_BIN_SEGMENTS = 10
MIN_GROUP = 5
ITERS = 500
TOL = 1e-5
DAMP = 0.5
MAX_ROUNDS = 10
MIN_SEG = 3

CENSAT_EXEMPT_BP = 1_000_000  # 00_depth_norm CENSAT_COMPRESSABLE_THRESHOLD
MAX_RATIO = 3                 # 22_save_matrix ABS_MAX_COVERAGE_RATIO
MAX_NFRAC = 0.5
MIN_TRAIN = 2000
MIN_SUPPORT = 8
HALF_A = frozenset([f'chr{i}' for i in range(1, 23, 2)] + ['chrX'])


def binned_curve(x, y, grp):
    """Median of y per 0.01 GC bin with >= 30 points from >= 10 groups; a GC bin short of
    that support is pooled with the next ones (low to high GC) until the pool has it, and a
    pool left short at the high end is dropped. Constant beyond the outermost supported
    points. f.cx holds the supported GC positions."""
    cx, cy = [], []
    pool = np.zeros(len(x), dtype=bool)
    for lo, hi in zip(GC_EDGES[:-1], GC_EDGES[1:]):
        k = ((x >= lo) & (x < hi)) | pool
        if k.sum() >= MIN_GC_BIN_ROWS and len(np.unique(grp[k])) >= MIN_GC_BIN_SEGMENTS:
            cx.append(np.median(x[k]))
            cy.append(np.median(y[k]))
            pool[:] = False
        else:
            pool = k
    if len(cx) < 3:
        raise ValueError(f'{len(cx)} supported GC bins')
    cx, cy = np.array(cx), np.array(cy)

    def f(g):
        return np.interp(np.asarray(g, dtype=float), cx, cy)
    f.cx = cx
    return f


def opt_partition(y, pen):
    """Exact optimal partition of y into constant-mean pieces of >= MIN_SEG values (L2 cost
    + pen per piece); returns a bool array marking piece starts."""
    n = len(y)
    out = np.zeros(n, dtype=bool)
    if n == 0:
        return out
    out[0] = True
    if n < 2 * MIN_SEG:
        return out
    cs, cs2 = np.r_[0, np.cumsum(y)], np.r_[0, np.cumsum(y * y)]
    F = np.full(n + 1, np.inf)
    F[0] = -pen
    last = np.zeros(n + 1, dtype=int)
    for t in range(MIN_SEG, n + 1):
        s = np.arange(0, t - MIN_SEG + 1)
        cost = np.maximum(cs2[t] - cs2[s] - (cs[t] - cs[s]) ** 2 / (t - s), 0)  # >= 0 despite rounding
        vals = F[s] + cost + pen
        j = int(np.argmin(vals))
        F[t], last[t] = vals[j], s[j]
    t = n
    while t > 0:
        out[last[t]] = True
        t = last[t]
    return out


def seg_breaks(code, y):
    """Segment starts for values already sorted by (chromosome code, position). The noise
    scale of the penalty uses adjacent differences within chromosomes; when more than half of
    them are exactly 0 (a feature that is 0 over long stretches), the non-zero differences
    set it. A feature without within-chromosome variation gives no breaks inside
    chromosomes."""
    same = code[1:] == code[:-1]
    d = np.abs(np.diff(y)[same])
    if len(d) and np.median(d) == 0 and (d > 0).any():
        d = d[d > 0]
    sig = (np.median(d) if len(d) else 0.0) / 0.6745 / np.sqrt(2)
    pen = 3 * sig ** 2 * np.log(len(y))
    out = np.zeros(len(y), dtype=bool)
    if sig == 0:
        out[np.r_[0, np.flatnonzero(~same) + 1][:len(y)]] = len(y) > 0
        return out
    for k in np.unique(code):
        idx = np.flatnonzero(code == k)
        out[idx] = opt_partition(y[idx], pen)
    return out


def backfit(y, x, grp, fx):
    for it in range(ITERS):
        s = pd.Series(y - fx)
        a = s.groupby(grp).transform('median').to_numpy()
        size = s.groupby(grp).transform('size').to_numpy()
        use = size >= MIN_GROUP
        f = binned_curve(x[use], (y - a)[use], grp[use])
        new = f(x)
        new = fx + DAMP * (new - np.median(new) - fx)
        delta = np.max(np.abs(new - fx))
        fx = new
        if delta < TOL:
            break
    return f, fx, it + 1, delta


def fit_seg(y, x, code, start):
    """Segment-based fit of log2 depth y on GC x (code/start: chromosome code and row start,
    any order). Returns (curve, segmentation rounds, last backfit change, converged)."""
    order = np.lexsort((start, code))
    brk = np.zeros(len(y), dtype=bool)
    fx = np.zeros(len(y))
    stable = False
    log = []
    for n in range(1, MAX_ROUNDS + 1):
        new = brk | seg_breaks(code[order], (y - fx)[order])
        if n > 1 and (new == brk).all():
            stable = True
            break
        brk = new
        grp = np.empty(len(y), dtype=int)
        grp[order] = np.cumsum(brk)
        f, fx, iters, delta = backfit(y, x, grp, fx)
        log.append((int(brk.sum()), iters, float(delta)))
    f.log = log  # per round: segments, backfit iterations, last curve change (log2)
    return f, n, delta, stable and delta < TOL


def read_censat(path):
    """CEN-SAT intervals per chromosome as an (n, 2) array of the bed's start and end."""
    b = pd.read_csv(path, sep='\t', header=None, usecols=[0, 1, 2], names=['chr', 's', 'e'], comment='#')
    return {ch: g[['s', 'e']].to_numpy() for ch, g in b.groupby('chr')}


def censat_masks(chrom, st, nd, censat):
    """(exempt, inside): exempt = chrY or overlap with a CEN-SAT interval > 1 Mb (stage 00);
    inside = fully inside any CEN-SAT interval (stage 22)."""
    exempt = chrom == 'chrY'
    inside = np.zeros(len(chrom), dtype=bool)
    for ch, iv in censat.items():
        m = np.flatnonzero(chrom == ch)
        if not len(m):
            continue
        s, e = st[m][:, None], nd[m][:, None]
        big = iv[np.abs(iv[:, 1] - iv[:, 0]) > CENSAT_EXEMPT_BP]
        exempt[m] |= (~((e < big[:, 0]) | (s > big[:, 1]))).any(axis=1)
        inside[m] = ((iv[:, 0] <= s) & (e <= iv[:, 1])).any(axis=1)
    return exempt, inside


def row_gc(path, df):
    """GC per row of df (NaN where unknown or the window is more than half N), or (None,
    reason) when the table was built on another reference: the last window end of a
    chromosome present in both differs."""
    with gzip.open(path, 'rt') if path.endswith('.gz') else open(path) as fh:
        g = pd.read_csv(fh, sep='\t', comment='#', names=['chr', 'st', 'nd', 'gc', 'nfrac'])
    ends = g.groupby('chr').nd.max()
    mine = df.groupby('chr').nd.max()
    shared = ends.index.intersection(mine.index)
    if not len(shared):
        return None, 'GC table has no chromosome of the depth table'
    differ = [c for c in shared if int(ends[c]) != int(mine[c])]
    if differ:
        return None, f'GC table is for another reference (window ends differ on {", ".join(differ[:3])})'
    g = g.assign(gc=g.gc.where(g.nfrac <= MAX_NFRAC))
    k = pd.Series(g.gc.to_numpy(dtype=float),
                  index=pd.MultiIndex.from_arrays([g.chr.to_numpy(), g.st.to_numpy(dtype=np.int64)]))
    key = pd.MultiIndex.from_arrays([df.chr.to_numpy(), df.st.to_numpy(dtype=np.int64)])
    return k.reindex(key).to_numpy(dtype=float), ''


def level_shift(d, base, apply, target):
    """log2 factor s for the corrected rows so that median(all rows) == target (bisection;
    the median is continuous and non-decreasing in s)."""
    lo, hi = -3.0, 3.0
    for _ in range(100):
        s = (lo + hi) / 2
        if np.median(np.where(apply, base * 2 ** s, d)) < target:
            lo = s
        else:
            hi = s
    return (lo + hi) / 2


def training_mask(depth, gc, exempt, inside):
    med = np.median(depth)
    return ~inside & (depth <= MAX_RATIO * med) & ~exempt & (depth > 0) & np.isfinite(gc)


def correct(df, gc, censat):
    """df: chr, st, nd, meandepth (chrM dropped); gc: per-row GC (NaN = unknown).
    Returns (corrected meandepth, training mask, diagnostics)."""
    d = df.meandepth.to_numpy(dtype=float)
    chrom, st, nd = df.chr.to_numpy(), df.st.to_numpy(dtype=np.int64), df.nd.to_numpy(dtype=np.int64)
    info = {'status': 'identity', 'reason': '', 'n_rows': len(d)}
    if not len(d) or not np.isfinite(d).all():
        return d.copy(), np.zeros(len(d), dtype=bool), dict(info, reason='non-finite or empty depth')
    med0 = float(np.median(d))
    exempt, inside = censat_masks(chrom, st, nd, censat)
    known = np.isfinite(gc)
    train = training_mask(d, gc, exempt, inside)
    apply = ~exempt & known
    info.update(n_train=int(train.sum()), n_exempt=int(exempt.sum()), n_apply=int(apply.sum()),
                median_before=med0)
    if train.sum() < MIN_TRAIN:
        return d.copy(), train, dict(info, reason=f'{int(train.sum())} training rows < {MIN_TRAIN}')
    y, x = np.log2(d[train]), gc[train]
    code, start = pd.factorize(chrom[train])[0], st[train]
    half = np.isin(chrom[train], list(HALF_A))
    try:
        f, rounds, delta, ok = fit_seg(y, x, code, start)
        fa = fit_seg(y[half], x[half], code[half], start[half])
        fb = fit_seg(y[~half], x[~half], code[~half], start[~half])
    except ValueError as e:
        return d.copy(), train, dict(info, reason=f'fit failed: {e}')
    info.update(rounds=rounds, delta=delta, n_support=len(f.cx), round_log=f.log,
                support_lo=float(f.cx[0]), support_hi=float(f.cx[-1]))
    if not (ok and fa[3] and fb[3]):
        return d.copy(), train, dict(info, reason='fit did not converge')
    if len(f.cx) < MIN_SUPPORT:
        return d.copy(), train, dict(info, reason=f'{len(f.cx)} supported GC bins < {MIN_SUPPORT}')
    signal = float(np.var(f(x)))
    noise = float(np.var(fa[0](x) - fb[0](x)) / 4)
    lam = float(np.clip(1 - noise / signal, 0, 1)) if signal > 0 else 0.0
    info.update(lam=lam, sd_curve=signal ** 0.5, sd_noise=noise ** 0.5)
    if lam == 0:
        return d.copy(), train, dict(info, reason='curve not reproducible between chromosome halves')
    base = d.copy()
    base[apply] = d[apply] / 2 ** (lam * f(gc[apply]))
    s = level_shift(d, base, apply, med0)
    out = d.copy()
    out[apply] = base[apply] * 2 ** s
    g = gc[apply]
    keep_before = ~inside & (d <= MAX_RATIO * med0)
    keep_after = ~inside & (out <= MAX_RATIO * np.median(out))
    info.update(status='corrected', shift=s, median_after=float(np.median(out)),
                response=float(2 ** (lam * (f(0.50) - f(0.36)))), sd_applied=lam * signal ** 0.5,
                curve_gc=[float(v) for v in f.cx], curve_log2=[float(lam * v - s) for v in f(f.cx)],
                n_outside_support=int(((g < f.cx[0]) | (g > f.cx[-1])).sum()),
                n_mask_flip=int((keep_before != keep_after).sum()))
    return out, train, info
