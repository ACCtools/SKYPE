"""Depth per copy (1N) from depth + B-allele fraction, instead of 1N = median depth / 2.

The median-depth 1N assumes the median window holds two copies, which is wrong for
near-triploid or near-tetraploid samples. Here 1N is the depth scale at which segments with
a recognizable allelic state (balanced, one third, one quarter, loss of heterozygosity) sit
on whole copy numbers.

Heterozygous sites (collect_sites): one window of 20 kb in the middle of each stage-22 row
(not fully inside a CEN-SAT interval, depth <= 3 x median), excluding stage-00 exempt rows
and positions inside the repeat track (segmental duplications and CEN-SAT), where paralogous
sequence makes false heterozygous sites. Reads: MAPQ >= 20, no secondary, supplementary,
duplicate or QC-fail. Bases: BQ >= 20. A position is callable when its depth is 8 .. 3 x
median row depth; it is heterozygous when the second allele has >= 3 reads and >= 10% of
the depth.

1. Segments: the training rows of each chromosome are split by depth (log2, clean rows
   only), by het-site rate and by mean minor-allele fraction, each with the exact optimal
   partition of depth_gc.seg_breaks; the breakpoints of the three are joined. Clean rows are
   training rows that do not overlap the repeat track.
2. Per segment of >= 10 rows: depth D = median over its clean rows, used only with >= 3
   clean rows and >= 20% of the segment (a row overlapping a duplication can carry
   paralogous reads, so such rows never set D). Het rate = het sites / callable bp, with a
   one-sided 95% Poisson upper bound. With >= 50 kb callable bp a segment is
     LOH      upper bound < LOH_RATE (2e-4 per bp). Absolute, so a genome that is nearly
              all LOH is still called.
     BAF      observed rate >= HET_RATE_MIN (3e-4) and >= 50 het sites: likelihood curve of
              the minor-allele fraction p (3.).
     unknown  otherwise; not used, like segments with less callable sequence.
3. BAF likelihood of one het site (depth n, minor count k): one allele has K reads,
   K ~ beta-binomial(n, p, rho), k = min(K, n - K), conditioned on the calling rule, mixed
   with weight eps with a uniform distribution over the callable k for non-germline sites.
   Per segment the curve LL(p) on a 0.005 grid profiles eps over EPS_GRID (at most 20000
   sites per segment). rho (allelic overdispersion, one per sample): profile likelihood
   over the BAF segments, then refitted on the pooled sites of the balanced segments with
   p fixed at 0.5. Cost of an allelic fraction q: cost(q) = min over p of
   2 (LLmax - LL(p)) + ((p - q) / SYS_P)^2 with a systematic offset of sd SYS_P = 0.03.
4. Allelic states (total C = 1..16, minor m):
     balanced     cost(1/2) <= 1 and p_hat >= 0.40: m = C/2 (C even)
     one third    otherwise cost(1/3) <= 1 and p_hat < 0.40: m = C/3
     one quarter  otherwise p_hat 0.20-0.28, cost(1/4) <= 1/4 and cost(1/3) > 4: m = C/4
                  (a subclonal 2+0 / 2+1 mixture also lands near 0.2, hence the tolerance)
     imbalanced   any other BAF segment: the best m = 1..C/2
     LOH          m = 0, credible at C only when one remaining minor copy of C would have
                  been seen (detection power of a site with fraction 1/C at the segment's
                  read depth).
   Anchors: balanced, one-third, one-quarter and credible LOH segments.
5. For N on a grid (median/16 .. median x 1.2): segment cost = min over allowed (C, m) of
   ((D/N - C) / SIGMA_C)^2 + cost(m/C) (SIGMA_C = 0.15 copies), capped at 16; a segment fits
   when its cost < 4. Candidates = local minima of the length-weighted mean cost. Halving N
   doubles every state, so N = the largest candidate whose anchor fit is within 0.10 of the
   best and for which >= 80% of the anchor and het-anchor length is within 3 sd: the fewest
   copies consistent with the anchors. N is refined as the length-weighted median of D/C
   over the segments that fit.
6. Status:
     insufficient_data     N = None: too few training rows or non-finite depth
     insufficient_anchors  N = None: too few usable or anchor segments
     poor_fit              no candidate meets the thresholds; N = the largest candidate
                           within 0.10 of the best anchor fit
     ambiguous_scale       a smaller candidate fits > 3 points more anchor length
                           (N_alternative)
     scale_sensitive       N moves by > 25% with SIGMA_C x 0.75 or x 1.25
     supported
   Warnings (status unchanged): BAF conflict > 20%, fit < 0.6, germline het rate (90th
   percentile) < 2.5 x LOH_RATE.
Limits: whole-genome doubling (every allele count shares a factor) gives an N that is that
factor too large; subclonal anchors over ~10% of the anchor length are fitted as clonal at
the smaller N.
"""
import os
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy import optimize, sparse, special, stats

import depth_gc

BIN = 100_000
MIN_SEG_BINS = 10
MIN_CLEAN, MIN_CLEAN_FRAC = 3, 0.20
MIN_SEG_CALLABLE = 50_000
MIN_SITES = 50
LOH_RATE = 2e-4
HET_RATE_MIN = 3e-4
CMAX = 16
SIGMA_C, SYS_P, CAP = 0.15, 0.03, 16.0
FIT_COST = 4.0
ADEQ_COST = 9.0
ANCHOR_COST = 1.0
QUARTER_TOL = 0.25
SEED_COST = 4.0
FIT_MARGIN = 0.10
ANCHOR_MIN = 0.80
AMBIG_GAP = 0.03
MIN_ANCHOR = 0.05
MIN_ANCHOR_SEGS = 3
SENS_FACTORS, SENS_MAX = (0.75, 1.25), 0.25
CONFLICT_WARN, FIT_WARN, REF_WARN = 0.20, 0.60, 2.5
EPS_GRID = (0.02, 0.05, 0.10, 0.15, 0.20, 0.30)
MAX_SITES = 20000
PGRID = np.round(np.arange(0.02, 0.5001, 0.005), 4)
PCOARSE = np.round(np.arange(0.02, 0.5001, 0.01), 4)
RHO_GRID = (0.0, 0.001, 0.0025, 0.005, 0.01, 0.02, 0.035, 0.05, 0.075, 0.10, 0.15)
RHO_MAX = 0.2
REF_QUANTILE = 0.90
HET_ANCHORS = ('balanced', 'third', 'quarter')
HET_KINDS = HET_ANCHORS + ('imbalanced',)
FAMILY = {'balanced': 2, 'third': 3, 'quarter': 4}
QS = sorted({round(m / C, 6) for C in range(2, CMAX + 1) for m in range(1, C // 2 + 1)})

# Het-site calling
SITE_WINDOW = 20_000
MIN_MAPQ, MIN_BQ, MIN_DEPTH, MIN_MINOR, MIN_FRAC = 20, 20, 8, 3, 0.10
SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x400 | 0x800
CHUNK = 100
BASE_CODE = np.full(256, 4, dtype=np.int64)
for _i, _b in enumerate(b'ACGT'):
    BASE_CODE[_b] = _i
MATCH_OPS = (0, 7, 8)       # M, =, X
QUERY_OPS = (1, 4)          # I, S
REF_OPS = (2, 3)            # D, N


# ---------------------------------------------------------------- het sites from the BAM

def allele_counts(bam, ch, s0, e0):
    """A/C/G/T counts per position of [s0, e0) (0-based) from reads with MAPQ >= MIN_MAPQ and
    none of SKIP_FLAGS, bases with BQ >= MIN_BQ (pysam count_coverage, vectorised per read;
    like it, reads without base qualities count nothing)."""
    width = e0 - s0
    keys = []
    for r in bam.fetch(ch, s0, e0):
        if (r.mapping_quality < MIN_MAPQ or r.flag & SKIP_FLAGS or r.query_sequence is None
                or r.query_qualities is None):
            continue
        rs, qs, ln = [], [], []
        rp, qp = r.reference_start, 0
        for op, n in r.cigartuples:
            if op in MATCH_OPS:
                a, b = max(rp, s0), min(rp + n, e0)
                if a < b:
                    rs.append(a - s0)
                    qs.append(qp + a - rp)
                    ln.append(b - a)
                rp += n
                qp += n
            elif op in QUERY_OPS:
                qp += n
            elif op in REF_OPS:
                rp += n
            if rp >= e0:
                break
        if not ln:
            continue
        ln = np.array(ln)
        cum = np.cumsum(ln)
        k = np.arange(cum[-1]) - np.repeat(cum - ln, ln)
        ref = np.repeat(rs, ln) + k
        qry = np.repeat(qs, ln) + k
        base = BASE_CODE[np.frombuffer(r.query_sequence.encode(), dtype=np.uint8)[qry]]
        qual = np.frombuffer(r.query_qualities, dtype=np.uint8)[qry]
        ok = (base < 4) & (qual >= MIN_BQ)
        keys.append(ref[ok] * 4 + base[ok])
    if not keys:
        return np.zeros((4, width), dtype=np.int64)
    return np.bincount(np.concatenate(keys), minlength=4 * width).reshape(width, 4).T


def merged_intervals(path):
    """Repeat-track intervals per chromosome, merged, as an (n, 2) array (0-based, half-open)."""
    b = pd.read_csv(path, sep='\t', header=None, usecols=[0, 1, 2], names=['chr', 's', 'e'], comment='#')
    out = {}
    for ch, g in b.groupby('chr'):
        iv = g.sort_values('s')[['s', 'e']].to_numpy()
        merged = [list(iv[0])]
        for s, e in iv[1:]:
            if s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[ch] = np.array(merged)
    return out


def masked(iv, s0, e0):
    """Boolean array over [s0, e0): True inside a repeat interval."""
    m = np.zeros(e0 - s0, dtype=bool)
    if iv is None or not len(iv):
        return m
    i = np.searchsorted(iv[:, 1], s0, side='right')
    while i < len(iv) and iv[i, 0] < e0:
        a, b = max(iv[i, 0], s0) - s0, min(iv[i, 1], e0) - s0
        m[a:b] = True
        i += 1
    return m


def _scan(args):
    import pysam
    bam_path, ch, rows, iv, cap = args
    sites, bins = [], []
    with pysam.AlignmentFile(bam_path) as bam:
        for st, s0, e0 in rows:
            cov = allele_counts(bam, ch, int(s0), int(e0))
            n = cov.sum(axis=0)
            minor = np.sort(cov, axis=0)[2]
            good = (n >= MIN_DEPTH) & (n <= cap) & ~masked(iv, int(s0), int(e0))
            het = good & (minor >= MIN_MINOR) & (minor >= MIN_FRAC * n)
            idx = np.flatnonzero(het)
            sites.append(np.c_[s0 + idx + 1, n[idx], minor[idx]])
            bins.append((st, int(good.sum()), len(idx)))
    s = np.vstack(sites) if sites else np.zeros((0, 3), dtype=np.int64)
    sdf = pd.DataFrame(s, columns=['pos', 'depth', 'minor']).assign(chrom=ch)
    bdf = pd.DataFrame(bins, columns=['start', 'callable', 'het']).assign(chrom=ch)
    return sdf[['chrom', 'pos', 'depth', 'minor']], bdf[['chrom', 'start', 'callable', 'het']]


def site_rows(df, depth, exempt, inside):
    """Rows sampled for het sites (stage-22 rows that stage 00 normalizes, depth > 0) and the
    per-position depth cap."""
    med = float(np.median(depth))
    ok = ~inside & ~exempt & (depth > 0) & (depth <= depth_gc.MAX_RATIO * med)
    w = pd.DataFrame({'chr': df.chr.to_numpy()[ok], 'st': df.st.to_numpy(dtype=np.int64)[ok],
                      'nd': df.nd.to_numpy(dtype=np.int64)[ok]})
    mid = (w.st - 1 + w.nd) // 2
    w['s0'] = (mid - SITE_WINDOW // 2).clip(lower=w.st - 1)
    w['e0'] = (w.s0 + SITE_WINDOW).clip(upper=w.nd)
    return w, int(depth_gc.MAX_RATIO * med)


def collect_sites(bam_path, df, depth, exempt, inside, repeat_bed, threads=1):
    """Het sites (chrom, pos, depth, minor) and per-row callable bp and het count (chrom,
    start, callable, het) of the sample BAM."""
    w, cap = site_rows(df, depth, exempt, inside)
    rpt = merged_intervals(repeat_bed)
    jobs = []
    for ch, wc in w.groupby('chr', sort=False):
        rows = wc[['st', 's0', 'e0']].to_numpy()
        for i in range(0, len(rows), CHUNK):
            jobs.append((bam_path, ch, rows[i:i + CHUNK], rpt.get(ch), cap))
    if not jobs:
        return (pd.DataFrame(columns=['chrom', 'pos', 'depth', 'minor']),
                pd.DataFrame(columns=['chrom', 'start', 'callable', 'het']))
    workers = max(1, min(int(threads), len(jobs)))
    if workers == 1:
        res = [_scan(j) for j in jobs]
    else:
        with Pool(workers) as pool:
            res = pool.map(_scan, jobs)
    return (pd.concat([r[0] for r in res], ignore_index=True),
            pd.concat([r[1] for r in res], ignore_index=True))


def bam_problem(bam_path):
    """Why the BAM cannot be read for het sites, or '' when it can."""
    if not bam_path or not os.path.isfile(bam_path):
        return f'no sample BAM ({bam_path})'
    if not any(os.path.isfile(bam_path + ext) for ext in ('.bai', '.csi')) and \
            not os.path.isfile(os.path.splitext(bam_path)[0] + '.bai'):
        return f'sample BAM is not indexed ({bam_path})'
    return ''


# ---------------------------------------------------------------- estimate of N

def bin_key(chrom, pos):
    return pd.MultiIndex.from_arrays([chrom, (np.asarray(pos) - 1) // BIN * BIN + 1])


def chrom_rank(c):
    c = c[3:] if c.startswith('chr') else c
    return int(c) if c.isdigit() else {'X': 23, 'Y': 24}.get(c, 25)


def features(df, depth, train, clean, sites, bins):
    """Per training row: chrom, start, depth, clean, callable, het, het rate, mean minor
    fraction."""
    t = pd.DataFrame({'chrom': df.chr.to_numpy()[train], 'start': df.st.to_numpy(dtype=np.int64)[train],
                      'depth': depth[train], 'clean': np.asarray(clean)[train]})
    key = pd.MultiIndex.from_arrays([t.chrom, t.start])
    b = bins.set_index(['chrom', 'start']).reindex(key)
    t['callable'] = b.callable.fillna(0).to_numpy(dtype=float)
    t['het'] = b.het.fillna(0).to_numpy(dtype=float)
    s = sites.assign(start=bin_key(sites.chrom, sites.pos).get_level_values(1), mf=sites.minor / sites.depth)
    mf = s.groupby(['chrom', 'start']).mf.agg(['mean', 'size']).reindex(key)
    t['mf'] = np.where(mf['size'].fillna(0).to_numpy() >= 3, mf['mean'].to_numpy(), np.nan)
    t['rate'] = np.where(t.callable >= 2000, t.het / np.maximum(t.callable, 1), np.nan)
    return t.sort_values(['chrom', 'start'], key=lambda c: c.map(chrom_rank) if c.name == 'chrom' else c,
                         ignore_index=True)


def breaks_of(t, values):
    """Segment starts (bool per row of t) from one feature; rows with NaN are skipped."""
    ok = np.isfinite(values)
    out = np.zeros(len(t), dtype=bool)
    code = pd.factorize(t.chrom)[0]
    out[np.flatnonzero(ok)] = depth_gc.seg_breaks(code[ok], values[ok])
    return out


def segment(t):
    d = np.where(t.clean.to_numpy(bool), np.log2(np.maximum(t.depth.to_numpy(), 1e-3)), np.nan)
    brk = breaks_of(t, d)
    brk |= breaks_of(t, np.sqrt(1000 * t.rate.to_numpy()))
    brk |= breaks_of(t, t.mf.to_numpy())
    brk |= np.r_[True, t.chrom.to_numpy()[1:] != t.chrom.to_numpy()[:-1]]
    return np.cumsum(brk)


def log_bb(k, n, p, rho):
    """log P(K = k), K ~ beta-binomial(n, p, rho); binomial at rho = 0."""
    if rho <= 1e-7:
        return stats.binom.logpmf(k, n, p)
    a, b = p * (1 / rho - 1), (1 - p) * (1 / rho - 1)
    return (special.gammaln(n + 1) - special.gammaln(k + 1) - special.gammaln(n - k + 1)
            + special.betaln(k + a, n - k + b) - special.betaln(a, b))


def call_prob(n, p, rho):
    """P(a het site at depth n with allele fraction p passes the calling rule)."""
    n = np.asarray(n, dtype=float)
    t = np.maximum(MIN_MINOR, np.ceil(MIN_FRAC * n))
    d = stats.binom(n, p) if rho <= 1e-7 else stats.betabinom(n, p * (1 / rho - 1), (1 - p) * (1 / rho - 1))
    return np.where(t <= n - t, d.cdf(n - t) - d.cdf(t - 1), 0.0)


def site_loglik(pairs, p, rho):
    """Log-likelihood of called het sites, unique (n, k) rows of pairs, at minor fractions p:
    array (len(EPS_GRID), len(pairs), len(p))."""
    n, k = pairs[:, :1].astype(float), pairs[:, 1:].astype(float)
    p = np.asarray(p, dtype=float)[None, :]
    l1, l2 = log_bb(k, n, p, rho), log_bb(n - k, n, p, rho)
    f = np.where(2 * k == n, l1, np.logaddexp(l1, l2))
    un, inv = np.unique(n[:, 0], return_inverse=True)
    f = f - np.log(np.maximum(call_prob(un[:, None], p, rho), 1e-300))[inv.ravel()]
    t = np.maximum(MIN_MINOR, np.ceil(MIN_FRAC * n))
    lu = -np.log(np.maximum(np.floor(n / 2) - t + 1, 1))
    return np.stack([np.logaddexp(np.log1p(-e) + f, np.log(e) + lu) for e in EPS_GRID])


def curves_at(counts, pairs, rho, p=PGRID):
    """LL(p) per segment (rows of the segment x pair count matrix), eps profiled."""
    return np.max([np.asarray(counts @ ll) for ll in site_loglik(pairs, p, rho)], axis=0)


def state_cost(curve, q):
    """cost(q) per segment (rows of curve) and allelic fraction q."""
    L = 2 * (curve.max(axis=1, keepdims=True) - curve)
    q = np.atleast_1d(np.asarray(q, dtype=float))
    return (L[:, :, None] + ((PGRID[None, :, None] - q[None, None, :]) / SYS_P) ** 2).min(axis=1)


def fit_rho(counts, pairs, n_sites):
    """Sample overdispersion; returns (rho, record)."""
    w = np.minimum(1.0, 200 / np.maximum(n_sites, 1))

    def profile(r):
        return -float((curves_at(counts, pairs, r, PCOARSE).max(axis=1) * w).sum())

    vals = [profile(r) for r in RHO_GRID]
    i = int(np.argmin(vals))
    opt = optimize.minimize_scalar(profile, bounds=(RHO_GRID[max(i - 1, 0)], RHO_GRID[min(i + 1, len(RHO_GRID) - 1)]),
                                   method='bounded', options={'xatol': 2e-4})
    rho = float(opt.x) if opt.fun < vals[i] else float(RHO_GRID[i])
    rec = {'rho_profile': rho, 'rounds': []}
    for _ in range(5):
        cv = curves_at(counts, pairs, rho)
        seed = (PGRID[cv.argmax(axis=1)] >= 0.40) & (state_cost(cv, 0.5)[:, 0] <= SEED_COST)
        if n_sites[seed].sum() < 100:
            break
        pooled = np.asarray(counts[seed].sum(axis=0)).ravel()

        def balanced(r):
            return -float(max(pooled @ ll[:, 0] for ll in site_loglik(pairs, [0.5], r)))

        f = optimize.minimize_scalar(balanced, bounds=(0.0, RHO_MAX), method='bounded', options={'xatol': 2e-4})
        new = float(f.x) if f.fun < balanced(0.0) else 0.0
        rec['rounds'].append({'segments': int(seed.sum()), 'sites': int(n_sites[seed].sum()), 'rho': new})
        done = abs(new - rho) < 2e-4
        rho = new
        if done:
            break
    return rho, rec


def segment_table(t, sites):
    """Segments with depth, het rate and its upper bound; returns (table, BAF curves (one row
    per segment, NaN where none), calibration)."""
    t = t.assign(seg=segment(t))
    st = t.groupby('seg').agg(chrom=('chrom', 'first'), start=('start', 'min'), end=('start', 'max'),
                              bins=('depth', 'size'), n_clean=('clean', 'sum'), callable=('callable', 'sum'),
                              het=('het', 'sum')).reset_index()
    st['depth'] = st.seg.map(t[t.clean].groupby('seg').depth.median())
    st.loc[st.n_clean < np.maximum(MIN_CLEAN, MIN_CLEAN_FRAC * st.bins), 'depth'] = np.nan
    st['end'] = st.end + BIN - 1
    st['rate'] = st.het / np.maximum(st.callable, 1)
    st['rate_upper'] = stats.chi2.ppf(0.95, 2 * (st.het + 1)) / 2 / np.maximum(st.callable, 1)
    s = sites.assign(start=bin_key(sites.chrom, sites.pos).get_level_values(1))
    s = s.merge(t[['chrom', 'start', 'seg', 'depth', 'clean']].rename(columns={'depth': 'bin_depth'}),
                on=['chrom', 'start'])
    pos = (s.clean & (s.bin_depth > 0)).to_numpy()   # read depth per unit of row depth, clean rows
    ratio = float(np.median(s.depth[pos] / s.bin_depth[pos])) if pos.any() else 1.0
    st['read_depth'] = st.depth * ratio
    st['n_sites'] = st.seg.map(s.groupby('seg').size()).fillna(0).astype(int)
    st['ok'] = (st.bins >= MIN_SEG_BINS) & (st.callable >= MIN_SEG_CALLABLE) & st.depth.notna()
    big = st.ok.to_numpy()
    w, r = st.bins.to_numpy(float)[big], st.rate.to_numpy()[big]
    o = np.argsort(r)
    ref = float(r[o][min(np.searchsorted(np.cumsum(w[o]) / w.sum(), REF_QUANTILE), len(w) - 1)]) if big.any() else np.nan
    cal = {'rho': 0.0, 'status': 'too few BAF segments; rho = 0', 'read_depth_ratio': ratio, 'ref_het_rate': ref}
    curves = np.full((len(st), len(PGRID)), np.nan)
    has = np.flatnonzero(big & (st.n_sites >= MIN_SITES).to_numpy())
    if not len(has):
        return st, curves, cal
    x = s[s.seg.isin(st.seg.to_numpy()[has])]
    x = x[x.assign(r=np.random.default_rng(0).random(len(x))).groupby('seg').r.rank(method='first') <= MAX_SITES]
    pairs, inv = np.unique(np.c_[x.depth.to_numpy(int), x.minor.to_numpy(int)], axis=0, return_inverse=True)
    row = pd.Series(np.arange(len(has)), index=st.seg.to_numpy()[has]).loc[x.seg.to_numpy()].to_numpy()
    counts = sparse.csr_matrix((np.ones(len(x)), (row, inv.ravel())), shape=(len(has), len(pairs)))
    n_sites = np.asarray(counts.sum(axis=1)).ravel()
    fit = ((st.rate_upper.to_numpy()[has] >= LOH_RATE) & (st.rate.to_numpy()[has] >= HET_RATE_MIN))
    rho = 0.0
    if fit.sum() >= 2:
        rho, rec = fit_rho(counts[fit], pairs, n_sites[fit])
        cal.update(rho=rho, status='fitted', n_sites=int(n_sites[fit].sum()), n_segments=int(fit.sum()), **rec)
    curves[has] = curves_at(counts, pairs, rho)
    st['p_rho0'] = np.nan
    st.loc[has, 'p_rho0'] = PGRID[curves_at(counts, pairs, 0.0).argmax(axis=1)]
    return st, curves, cal


def classify(st, curves):
    """Allelic state per segment; adds state and p_hat."""
    st = st.copy()
    ok = st.ok.to_numpy(bool)
    loh = ok & (st.rate_upper.to_numpy() < LOH_RATE)
    baf = (ok & ~loh & (st.rate.to_numpy() >= HET_RATE_MIN) & (st.n_sites.to_numpy() >= MIN_SITES)
           & np.isfinite(curves).all(axis=1))
    state = np.where(loh, 'loh', 'unknown').astype(object)
    st['p_hat'] = np.nan
    if baf.any():
        cv = curves[baf]
        ph = PGRID[cv.argmax(axis=1)]
        c = state_cost(cv, [0.5, 1 / 3, 0.25])
        bal = (c[:, 0] <= ANCHOR_COST) & (ph >= 0.40)
        third = ~bal & (c[:, 1] <= ANCHOR_COST) & (ph < 0.40)
        quarter = (~bal & ~third & (ph >= 0.20) & (ph <= 0.28) & (c[:, 2] <= ANCHOR_COST * QUARTER_TOL)
                   & (c[:, 1] > FIT_COST))
        state[baf] = np.select([bal, third, quarter], ['balanced', 'third', 'quarter'], 'imbalanced')
        st.loc[baf, 'p_hat'] = ph
    st['state'] = state
    return st


def allele_costs(s, cv, free=False):
    """(BAF cost, minor copies) per segment and total C = 1..CMAX. free: every BAF segment
    may take any minor m (used for the conflict check)."""
    kind = s.state.to_numpy()
    out = np.full((len(s), CMAX), CAP)
    minor = np.zeros((len(s), CMAX), dtype=int)
    out[kind == 'loh'] = 0.0
    het = np.flatnonzero(np.isin(kind, HET_KINDS))
    if not len(het):
        return out, minor
    Q = state_cost(cv[het], QS)
    col = {q: j for j, q in enumerate(QS)}
    for C in range(2, CMAX + 1):
        ms = np.arange(1, C // 2 + 1)
        qc = Q[:, [col[round(m / C, 6)] for m in ms]]
        b = qc.argmin(axis=1)
        out[het, C - 1] = qc[np.arange(len(het)), b]
        minor[het, C - 1] = ms[b]
        if free:
            continue
        for k, d in FAMILY.items():
            h = kind[het] == k
            out[het[h], C - 1] = CAP if C % d else Q[h, col[round(1 / d, 6)]]
            minor[het[h], C - 1] = C // d
    return out, minor


def loh_credible(s, rho):
    """Per segment and C = 1..CMAX: False where an LOH segment could still hold one minor
    copy of C that the calling rule would mostly miss."""
    cred = np.ones((len(s), CMAX), dtype=bool)
    loh = np.flatnonzero(s.state.to_numpy() == 'loh')
    if len(loh):
        n = np.maximum(np.round(s.read_depth.to_numpy()[loh]), 1)[:, None]
        power = call_prob(n, 1.0 / np.arange(2, CMAX + 1)[None, :], rho)
        cred[loh, 1:] = s.rate_upper.to_numpy()[loh, None] / np.maximum(power, 1e-9) < LOH_RATE
    return cred


def costs(D, bc, N, sigma_c=SIGMA_C):
    """Segment x N matrices of the capped best-state cost and its total C."""
    c = D[:, None, None] / N[None, :, None]
    C = np.arange(1, CMAX + 1)[None, None, :]
    total = ((c - C) / sigma_c) ** 2 + bc[:, None, :]
    return np.minimum(total.min(axis=2), CAP), total.argmin(axis=2) + 1


def scores(cost, Cg, kind, cred, w, limit=FIT_COST):
    """Anchor fit A, het-anchor fit H, fit F (segments with cost < limit), mean cost E per
    column, and the anchor mask."""
    fits = cost < limit
    hanc = np.isin(kind, HET_ANCHORS)[:, None]
    anc = hanc | ((kind == 'loh')[:, None] & np.take_along_axis(cred, Cg - 1, axis=1))
    wa, wh = w[:, None] * anc, w[:, None] * hanc
    A = (fits * wa).sum(axis=0) / np.maximum(wa.sum(axis=0), 1e-12)
    H = (fits * wh).sum(axis=0) / max(float(wh.sum()), 1e-12)
    F = (fits * w[:, None]).sum(axis=0) / w.sum()
    E = (cost * w[:, None]).sum(axis=0) / w.sum()
    return A, H, F, E, anc


def fit_scale(s, bc, cred, median_depth, sigma_c=SIGMA_C):
    """N and the candidate scores."""
    D, w, kind = s.depth.to_numpy(), s.bins.to_numpy(float), s.state.to_numpy()
    grid = np.exp(np.linspace(np.log(median_depth / CMAX), np.log(median_depth * 1.2), 1500))
    cost, Cg = costs(D, bc, grid, sigma_c)
    A, H, F, E, _ = scores(cost, Cg, kind, cred, w)
    A3, H3 = scores(cost, Cg, kind, cred, w, ADEQ_COST)[:2]
    cand = sorted({int(np.argmin(E))} | {i for i in range(1, len(E) - 1) if E[i] <= E[i - 1] and E[i] <= E[i + 1]})
    best = max(A[i] for i in cand)
    ok = [i for i in cand if A[i] >= best - FIT_MARGIN and A3[i] >= ANCHOR_MIN and H3[i] >= ANCHOR_MIN]
    pick = max(ok) if ok else max(i for i in cand if A[i] >= best - FIT_MARGIN)
    better = [i for i in cand if i < pick and A[i] > A[pick] + AMBIG_GAP]
    alt = max(better, key=lambda i: (A[i], i)) if better else None

    def refine(i):
        c0, C0 = costs(D, bc, grid[[i]], sigma_c)
        f0 = c0[:, 0] < FIT_COST
        if not f0.any():
            return float(grid[i])
        r = D[f0] / C0[f0, 0]
        o = np.argsort(r)
        return float(r[o][np.searchsorted(np.cumsum(w[f0][o]) / w[f0].sum(), 0.5)])

    top = sorted(cand, key=lambda i: (-A[i], -F[i]))[:5]
    return refine(pick), {
        'adequate': bool(ok), 'N_alternative': refine(alt) if alt is not None else None,
        'anchor_fit_alternative': float(A[alt]) if alt is not None else None,
        'candidates': [(round(float(grid[i]), 3), round(float(A[i]), 3), round(float(H[i]), 3), round(float(F[i]), 3),
                        round(float((D / grid[i] * w).sum() / w.sum()), 2)) for i in top]}


def solve(st, curves, cal, median_depth):
    """States, N, status and diagnostics from segment_table output."""
    st = classify(st, curves)
    use = st.state.isin(('loh',) + HET_KINDS).to_numpy()
    s = st[use].reset_index(drop=True)
    cv = curves[use]
    kind, w = s.state.to_numpy(), s.bins.to_numpy(float)
    diag = {'median_depth': median_depth, 'N_median_half': median_depth / 2, 'calibration': cal,
            'n_segments': int(len(s)), 'states': {k: int((kind == k).sum()) for k in ('loh',) + HET_KINDS}}
    n_het = int(np.isin(kind, HET_ANCHORS).sum())
    anchor_w = float(w[np.isin(kind, ('loh',) + HET_ANCHORS)].sum())
    reason = ('< 5 usable segments' if len(s) < 5 or w.sum() < 300 else
              'no balanced / one-third / one-quarter segment' if n_het < 1 else
              f'< {MIN_ANCHOR_SEGS} anchor segments' if n_het + (kind == 'loh').sum() < MIN_ANCHOR_SEGS else
              f'anchors < {MIN_ANCHOR:.0%} of the length' if anchor_w < MIN_ANCHOR * w.sum() else '')
    if reason:
        return None, s, dict(diag, N=None, status='insufficient_anchors', reason=reason)
    bc, minor = allele_costs(s, cv)
    cred = loh_credible(s, cal['rho'])
    N, fit = fit_scale(s, bc, cred, median_depth)
    D = s.depth.to_numpy()
    cN, CN = costs(D, bc, np.array([N, N / 2, 2 * N]))
    A, H, F, E, anc = scores(cN, CN, kind, cred, w)
    A3, H3 = scores(cN, CN, kind, cred, w, ADEQ_COST)[:2]
    C = CN[:, 0]
    m = minor[np.arange(len(s)), C - 1].astype(float)
    m[kind == 'loh'] = np.where(cred[np.arange(len(s)), C - 1], 0, np.nan)[kind == 'loh']
    s = s.assign(C=C, minor=m, cost=cN[:, 0], fits=cN[:, 0] < FIT_COST, anchor=anc[:, 0])
    free, _ = allele_costs(s, cv, free=True)
    depth_ok = ((D[:, None] / N - np.arange(1, CMAX + 1)[None, :]) / SIGMA_C) ** 2 <= FIT_COST
    het = np.isin(kind, HET_KINDS)
    conflict = het & ~(depth_ok & (free < FIT_COST)).any(axis=1)
    sens = []
    for f in SENS_FACTORS:
        Nf, _ = fit_scale(s, bc, cred, median_depth, SIGMA_C * f)
        sens.append({'sigma_c': SIGMA_C * f, 'N': Nf, 'change': Nf / N - 1})
    flags = (['poor_fit'] if not fit['adequate'] else []) \
        + (['ambiguous_scale'] if fit['N_alternative'] is not None else []) \
        + (['scale_sensitive'] if any(abs(x['change']) > SENS_MAX for x in sens) else [])
    cf = float(w[conflict].sum() / w[het].sum()) if het.any() else 0.0
    warn = ([f'BAF conflict {cf:.2f} > {CONFLICT_WARN}'] if cf > CONFLICT_WARN else []) \
        + ([f'fit {F[0]:.2f} < {FIT_WARN}'] if F[0] < FIT_WARN else []) \
        + ([f'germline het rate {cal["ref_het_rate"]:.2e} < {REF_WARN} x LOH_RATE']
           if not cal['ref_het_rate'] >= REF_WARN * LOH_RATE else [])
    diag.update(fit, N=N, status=flags[0] if flags else 'supported', flags=flags, warnings=warn,
                ploidy=float((D / N * w).sum() / w.sum()), median_copies=median_depth / N,
                anchor_fit=float(A[0]), het_anchor_fit=float(H[0]), fit_fraction=float(F[0]),
                anchor_fit_3sd=float(A3[0]), het_anchor_fit_3sd=float(H3[0]),
                anchor_fit_half=float(A[1]), anchor_fit_double=float(A[2]),
                anchor_fraction=float(w[anc[:, 0]].sum() / w.sum()), loh_fraction=float(w[kind == 'loh'].sum() / w.sum()),
                baf_conflict_fraction=cf, sensitivity=sens)
    return N, s, diag


def estimate(df, depth, train, clean, sites, bins):
    """Returns (N or None, per-segment table, diagnostics). depth: per-row depth of df;
    train / clean: row masks (clean = training rows outside the repeat track)."""
    if train.sum() < 5 * MIN_SEG_BINS or not np.isfinite(depth[train]).all():
        return None, pd.DataFrame(), {'N': None, 'status': 'insufficient_data',
                                      'reason': 'too few training rows or non-finite depth'}
    t = features(df, depth, train, clean, sites, bins)
    st, curves, cal = segment_table(t, sites)
    return solve(st, curves, cal, float(np.median(depth)))


def repeat_free(df, path):
    """Per row of df (1-based inclusive st..nd): True when it does not overlap the bed."""
    b = pd.read_csv(path, sep='\t', header=None, usecols=[0, 1, 2], names=['chr', 's', 'e'], comment='#')
    out = np.ones(len(df), dtype=bool)
    chrom = df.chr.to_numpy()
    s0, e0 = df.st.to_numpy(dtype=np.int64) - 1, df.nd.to_numpy(dtype=np.int64)
    for ch, g in b.groupby('chr'):
        m = np.flatnonzero(chrom == ch)
        if not len(m):
            continue
        g = g.sort_values('s')
        s, e = g.s.to_numpy(), np.maximum.accumulate(g.e.to_numpy())
        i = np.searchsorted(e, s0[m], side='right')   # first interval ending after the row start
        hit = i < len(s)
        out[m[hit]] = s[i[hit]] >= e0[m[hit]]
    return out
