#!/usr/bin/env python
"""One depth-triggered NClose rescue pass, using raw reads or local OLC.

D detection follows analysis/HCC1937_raw_cusum_20260910/D_method.ko.md.
This stage prepares an augmented stage-01 handoff. pipeline.py installs it and
runs stages 10--23 once; this script never recursively invokes the pipeline.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from functools import lru_cache
import json
import logging
import math
import multiprocessing
from pathlib import Path
import pickle
import re
import shutil
import subprocess
import time

import numpy as np
import pandas as pd
import pysam

from local_olc import OLCConfig, assemble
from skype_utils import VCF_TYPE4_MIN_SPAN

BIN = 100000


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2)+'\n')


def write_tsv(path, records, columns=None):
    frame = pd.DataFrame(records)
    if frame.empty and columns:
        frame = pd.DataFrame(columns=columns)
    frame.to_csv(path, sep='\t', index=False)


def mad_sigma(values):
    return float(np.median(np.abs(values-np.median(values)))/0.6745)


def optimal_partition(y, penalty=96., min_size=8):
    n = len(y)
    if n < 2*min_size:
        return []
    y = y-np.mean(y)
    s, q = np.r_[0., np.cumsum(y)], np.r_[0., np.cumsum(y*y)]
    cost, previous = np.full(n+1, np.inf), np.full(n+1, -1, dtype=int)
    cost[0] = -penalty
    for end in range(min_size, n+1):
        starts = np.arange(end-min_size+1)
        values = cost[starts]+np.maximum(0., q[end]-q[starts]-(s[end]-s[starts])**2/(end-starts))+penalty
        best = int(np.argmin(values))
        cost[end], previous[end] = values[best], best
    cuts, end = [], n
    while previous[end] > 0:
        end = previous[end]
        cuts.append(int(end))
    return sorted(cuts)


def depth_runs(coords, observed, predicted):
    starts = [0]+[i for i in range(1, len(coords))
                  if coords[i][0] != coords[i-1][0] or int(coords[i][1])-int(coords[i-1][1]) != BIN]+[len(coords)]
    return [dict(chrom=str(coords[a][0]), run_id=i, start0=int(coords[a][1])-1,
                 b=np.asarray(observed[a:b], dtype=float), p=np.asarray(predicted[a:b], dtype=float))
            for i, (a, b) in enumerate(zip(starts[:-1], starts[1:])) if b > a]


def boundary_stats(run, cut, left, right, sigma):
    b, p = run['b'], run['p']
    yl, yr = 2*np.sqrt(np.maximum(b[left:cut], 0)), 2*np.sqrt(np.maximum(b[cut:right], 0))
    nl, nr = len(yl), len(yr)
    bl, br = float(np.median(b[left:cut])), float(np.median(b[cut:right]))
    y = np.r_[yl, yr]/sigma
    x = np.arange(len(y), dtype=float)-(len(y)-1)/2
    step_gain = (np.mean(y[nl:])-np.mean(y[:nl]))**2*nl*nr/(nl+nr)
    line_gain = np.dot(x, y)**2/np.dot(x, x)
    return dict(observed_left=bl, observed_right=br, observed_step=br-bl,
        pred_range=float(np.ptp(p[left:right])),
        z2=float((np.median(yr)-np.median(yl))**2/sigma**2*nl*nr/(nl+nr)),
        hom=max(mad_sigma(yl), mad_sigma(yr))/sigma,
        nl=nl, nr=nr, step_over_line=float(step_gain-line_gain))


@lru_cache(maxsize=64)
def query_model_templates(n, flank):
    """Centered step/ramp templates, with ramps averaged over each depth bin."""
    x=np.arange(n,dtype=float)
    templates=[];starts=[];ends=[];parameters=[]
    for cut in range(flank,n-flank+1):
        templates.append((x>=cut).astype(float));starts.append(cut);ends.append(cut);parameters.append(3)
    for left in range(flank,n-flank):
        for right in range(left+1,n-flank+1):
            def primitive(z):
                return np.where(z<=left,0.,np.where(z>=right,z-(left+right)/2,(z-left)**2/(2*(right-left))))
            templates.append(primitive(x+1)-primitive(x))
            starts.append(left);ends.append(right);parameters.append(4)
    matrix=np.asarray(templates)
    matrix-=matrix.mean(axis=1,keepdims=True)
    return matrix,np.sum(matrix*matrix,axis=1),np.asarray(starts),np.asarray(ends),np.asarray(parameters)


def adaptive_query_interval(run,cut,left,right,sigma,direction,args):
    """Union competitive transition locations across adjacent-segment contexts.

    This is a descriptive location envelope, not a calibrated confidence interval.
    It changes the read query only; the D boundary/evidence gates are unchanged.
    """
    lower=upper=cut;windows=[];seen=set()
    for context in args.query_context_bins:
        a,b=max(left,cut-context),min(right,cut+context)
        if (a,b) in seen or b-a<2*args.query_flank_bins:
            continue
        seen.add((a,b))
        y=2*np.sqrt(np.maximum(run['b'][a:b],0))/sigma;y-=y.mean()
        matrix,norm,starts,ends,parameters=query_model_templates(len(y),args.query_flank_bins)
        cross=matrix@y
        cost=np.maximum(0.,np.dot(y,y)-cross*cross/norm)+parameters*np.log(len(y))
        cost[cross*direction<=0]=np.inf
        if not np.isfinite(cost).any():
            continue
        best=int(np.argmin(cost));keep=cost<=cost[best]+args.query_model_delta
        lo=a+int(starts[keep].min());hi=a+int(ends[keep].max())
        lower=min(lower,lo);upper=max(upper,hi)
        windows.append(dict(context_start0=run['start0']+a*BIN,context_end0=run['start0']+b*BIN,
            best_model='step' if parameters[best]==3 else 'ramp',best_start0=run['start0']+(a+int(starts[best]))*BIN,
            best_end0=run['start0']+(a+int(ends[best]))*BIN,best_cost=float(cost[best]),
            envelope_start0=run['start0']+lo*BIN,envelope_end0=run['start0']+hi*BIN,
            competitive_models=int(keep.sum())))
    return dict(transition_start0=run['start0']+lower*BIN,transition_end0=run['start0']+upper*BIN,
        query_start0=max(run['start0'],run['start0']+lower*BIN-args.query_padding),
        query_end0=min(run['start0']+len(run['b'])*BIN,run['start0']+upper*BIN+args.query_padding),
        query_model_windows=windows)


def detect_depth(coords, observed, predicted, args):
    runs = depth_runs(coords, observed, predicted)
    differences = defaultdict(list)
    for run in runs:
        b = run['b']
        differences[run['chrom']].append(np.diff(2*np.sqrt(np.maximum(b, 0)))[(b[:-1] > 0) & (b[1:] > 0)])
    all_diff = np.concatenate([d for values in differences.values() for d in values]) if runs else np.array([])
    global_sigma = mad_sigma(all_diff)/np.sqrt(2) if len(all_diff) else 0.
    sigmas = {}
    for chrom, values in differences.items():
        d = np.concatenate(values)
        sigma = mad_sigma(d)/np.sqrt(2) if len(d) >= 100 else float('nan')
        sigmas[chrom] = sigma if np.isfinite(sigma) and sigma > 0 else global_sigma
    candidates, audit = [], []
    for run in runs:
        sigma = sigmas[run['chrom']]
        if not np.isfinite(sigma) or sigma <= 0:
            continue
        cuts = optimal_partition(2*np.sqrt(np.maximum(run['b'], 0))/sigma, args.dp_penalty, args.min_segment_bins)
        edges = [0]+cuts+[len(run['b'])]
        for i, cut in enumerate(cuts, 1):
            local = boundary_stats(run, cut, max(edges[i-1], cut-args.local_bins), min(edges[i+1], cut+args.local_bins), sigma)
            full = boundary_stats(run, cut, edges[i-1], edges[i+1], sigma)
            shape = lambda r: r['pred_range'] <= args.prediction_ratio*abs(r['observed_step']) and r['hom'] <= args.max_hom
            route, reason = '', 'local_shape'
            if shape(local):
                reason = 'insufficient_evidence'
                if local['z2'] >= args.dp_penalty:
                    route = 'local'
                elif shape(full) and np.sign(full['observed_step']) == np.sign(local['observed_step']) and full['z2'] >= args.dp_penalty and full['step_over_line'] > np.log(full['nl']+full['nr']):
                    route = 'long'
            pos = run['start0']+cut*BIN
            row = dict(chrom=run['chrom'], run_id=run['run_id'], representative0=pos,
                query_start0=max(run['start0'], pos-args.query_padding),
                query_end0=min(run['start0']+len(run['b'])*BIN, pos+args.query_padding),
                direction='rise' if local['observed_step'] > 0 else 'drop', sigma=sigma,
                route=route, status='candidate' if route else reason, **local,
                **{'long_'+k: v for k, v in full.items()})
            audit.append(row)
            if route:
                row['candidate_id'] = f'D{len(candidates)+1:03d}'
                row['fixed_query_start0'],row['fixed_query_end0']=row['query_start0'],row['query_end0']
                row['query_mode']=args.query_mode
                if args.query_mode=='adaptive':
                    row.update(adaptive_query_interval(run,cut,edges[i-1],edges[i+1],sigma,np.sign(local['observed_step']),args))
                else:
                    row.update(transition_start0=pos,transition_end0=pos,query_model_windows=[])
                candidates.append(row.copy())
    queries = []
    for candidate in sorted(candidates,key=lambda c:(c['run_id'],c['query_start0'],c['representative0'])):
        if queries and queries[-1]['run_id'] == candidate['run_id'] and candidate['query_start0'] <= queries[-1]['end0']:
            queries[-1]['end0'] = max(queries[-1]['end0'], candidate['query_end0'])
            queries[-1]['candidate_ids'].append(candidate['candidate_id'])
        else:
            queries.append(dict(query_id=f'Q{len(queries)+1:03d}', chrom=candidate['chrom'],
                start0=candidate['query_start0'], end0=candidate['query_end0'], run_id=candidate['run_id'],
                candidate_ids=[candidate['candidate_id']]))
    return candidates, queries, audit, runs, sigmas


def save_depth_detection(prefix, outdir, args):
    with (prefix/'23_input.pkl').open('rb') as f:
        metadata = pickle.load(f)
    sl = slice(metadata['B_depth_start'], metadata['B_depth_end'])
    b = np.load(prefix/'B.npy').reshape(-1)[sl]
    p = np.load(prefix/'predict_B.npy').reshape(-1)[sl]
    coords = metadata['chr_filt_st_list']
    candidates, queries, audit, runs, sigmas = detect_depth(coords, b, p, args)
    write_tsv(outdir/'depth_candidates.tsv', candidates, ['candidate_id','chrom','representative0','route'])
    write_tsv(outdir/'depth_boundary_audit.tsv', audit)
    write_tsv(outdir/'query_regions.tsv', queries, ['query_id','chrom','start0','end0'])
    with (outdir/'query_regions.bed').open('w') as f:
        for q in queries:
            f.write(f"{q['chrom']}\t{q['start0']}\t{q['end0']}\t{q['query_id']}\n")
    write_json(outdir/'queries.json', queries)
    write_json(outdir/'query_models.json',[{key:c[key] for key in ('candidate_id','chrom','representative0',
        'fixed_query_start0','fixed_query_end0','transition_start0','transition_end0','query_start0','query_end0','query_mode','query_model_windows')}
        for c in candidates])
    write_json(outdir/'depth_summary.json', dict(valid_bins=len(b), valid_runs=len(runs),
        candidate_boundaries=len(candidates), query_regions=len(queries), noise=sigmas,
        query_mode=args.query_mode,query_bases=sum(q['end0']-q['start0'] for q in queries),
        local=sum(c['route'] == 'local' for c in candidates), long=sum(c['route'] == 'long' for c in candidates)))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(outdir/'depth_candidates.pdf') as pdf:
        if not candidates:
            fig, ax = plt.subplots(); ax.text(.1,.5,'D method: no candidate boundaries'); ax.axis('off'); pdf.savefig(fig); plt.close(fig)
        for c in candidates:
            run = runs[c['run_id']]
            edges = run['start0']+np.arange(len(run['b'])+1)*BIN
            fig, axes = plt.subplots(2, 1, figsize=(12, 7), layout='constrained')
            for ax, pad in zip(axes, (2_000_000, 600_000)):
                ax.stairs(run['b'], edges/1e6, baseline=None, label='Normalized observed', color='#287aa5')
                ax.stairs(run['p'], edges/1e6, baseline=None, label='Refitted prediction', color='#d78026')
                ax.axvspan(c['query_start0']/1e6,c['query_end0']/1e6,color='#39a298',alpha=.15,label='Read query')
                ax.axvspan(c['fixed_query_start0']/1e6,c['fixed_query_end0']/1e6,color='#d4b44c',alpha=.25,label='Original +/-200 kb')
                ax.set(xlim=(min(c['representative0']-pad,c['query_start0']-BIN)/1e6,
                    max(c['representative0']+pad,c['query_end0']+BIN)/1e6), ylabel='Depth', xlabel=c['chrom']+' (Mb)')
                ax.legend()
            fig.suptitle(f"{c['candidate_id']} | {c['chrom']}:{c['representative0']} | {c['route']} | z²={c['z2']:.1f}, hom={c['hom']:.2f}")
            pdf.savefig(fig); plt.close(fig)
    return queries


def cigar_alignment(chrom, start, strand, cigar, mapq, nm, lengths):
    ops = [(int(n), op) for n, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]
    qlen = sum(n for n, op in ops if op in 'MIS H=X'.replace(' ', ''))
    left = right = 0
    for n, op in ops:
        if op not in 'SH':
            break
        left += n
    for n, op in reversed(ops):
        if op not in 'SH':
            break
        right += n
    qst, qend = left, qlen-right
    if strand == '-':
        qst, qend = qlen-qend, qlen-qst
    refspan = sum(n for n, op in ops if op in 'MDN=X')
    aligned = sum(n for n, op in ops if op in 'MI=X')
    return dict(chrom=chrom, start0=int(start), end0=int(start)+refspan, strand=strand,
        qstart=qst, qend=qend, qlen=qlen, mapq=int(mapq), nm=int(nm),
        identity=1-int(nm)/max(refspan, aligned, 1), cigar=cigar,
        chrom_length=int(lengths.get(chrom, start+refspan)))


def alignment_key(row):
    return row['chrom'], row['start0'], row['end0'], row['strand'], row['qstart'], row['qend']


def extract_region_reads(bam, query, outdir):
    """Collect full molecules and primary/SA chains, including remote primaries."""
    lengths = dict(zip(bam.references, bam.lengths))
    chains, sequences, partners = defaultdict(dict), {}, defaultdict(list)
    def collect(record):
        if record.is_unmapped or record.is_secondary:
            return
        name = record.query_name
        row = cigar_alignment(record.reference_name, record.reference_start,
            '-' if record.is_reverse else '+', record.cigarstring, record.mapping_quality,
            record.get_tag('NM') if record.has_tag('NM') else 0, lengths)
        chains[name][alignment_key(row)] = row
        if record.query_sequence and 'H' not in record.cigarstring:
            sequences[name] = record.get_forward_sequence().upper()
        if record.has_tag('SA'):
            for entry in record.get_tag('SA').strip(';').split(';'):
                chrom, pos, strand, cigar, mapq, nm = entry.split(',')
                row = cigar_alignment(chrom, int(pos)-1, strand, cigar, int(mapq), int(nm), lengths)
                chains[name][alignment_key(row)] = row
                partners[chrom].append((int(pos)-1, int(pos), name))
    for record in bam.fetch(query['chrom'], query['start0'], query['end0']):
        collect(record)
    missing = set(chains)-set(sequences)
    for chrom, entries in list(partners.items()):
        intervals = sorted((a, b) for a, b, name in entries if name in missing)
        merged = []
        for a, b in intervals:
            if merged and a <= merged[-1][1]+1000:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        for a, b in merged:
            for record in bam.fetch(chrom, max(0, a), min(lengths[chrom], b)):
                if record.query_name in missing:
                    collect(record)
    chains = {name: sorted(rows.values(), key=lambda r: (r['qstart'],r['qend'])) for name, rows in chains.items()}
    with (outdir/'reads.fa').open('w') as f:
        for name, sequence in sorted(sequences.items()):
            f.write(f'>{name}\n{sequence}\n')
    write_json(outdir/'bam_chains.json', chains)
    write_json(outdir/'extraction.json', dict(reads=len(chains), full_sequences=len(sequences),
        sequence_missing=sorted(set(chains)-set(sequences)), query=query))
    return sequences, chains


def split_large_indels(row, minimum):
    """Keep large biological gaps out of anchor identity and terminal geometry.

PAF cs is in reference-forward order even for a reverse-strand query. Split
pieces retain the original molecule's query coordinates and exact cs strings.
"""
    cigar_ops=[(int(n),op) for n,op in re.findall(r'(\d+)([MIDNSHP=X])',row.get('cigar',''))]
    if not any(op=='N' or (op in 'ID' and n>=minimum) for n,op in cigar_ops):
        return [row]
    fields=row['paf'].split('\t') if 'paf' in row else None
    tokens=[]
    if fields:
        cs=next(x[5:] for x in fields[12:] if x.startswith('cs:Z:'))
        for token in re.findall(r':\d+|\*[A-Za-z]{2}|[=+\-][A-Za-z]+|~[A-Za-z]{2}\d+[A-Za-z]{2}',cs):
            symbol=token[0]
            if symbol in ':=':
                n=int(token[1:]) if symbol==':' else len(token)-1
                tokens.append((n,'M',':'+str(n),n))
            elif symbol=='*': tokens.append((1,'X',token,0))
            elif symbol in '+-': tokens.append((len(token)-1,'I' if symbol=='+' else 'D',token,0))
            else: tokens.append((int(re.search(r'\d+',token).group()),'N',token,0))
    else:
        tokens=[(n,op,'',n if op in 'M=' else 0) for n,op in cigar_ops if op not in 'SH P'.replace(' ','')]
    large_error=sum(n for n,op,_,_ in tokens if op in 'ID' and n>=minimum)
    remaining_span=sum(n for n,op,_,_ in tokens if op!='N' and not (op in 'ID' and n>=minimum))
    bam_identity=1-max(0,row.get('nm',0)-large_error)/max(remaining_span,1)
    qpos=0;rpos=row['start0'];qbegin=0;rbegin=rpos;parts=[];result=[]
    def flush():
        if qpos<=qbegin or rpos<=rbegin:
            return
        qst,qend=(row['qstart']+qbegin,row['qstart']+qpos) if row['strand']=='+' else (row['qend']-qpos,row['qend']-qbegin)
        matches=sum(x[3] for x in parts)
        span=sum(x[0] for x in parts)
        cg=''.join(f'{n}{op}' for n,op,_,_ in parts)
        piece=dict(row,start0=rbegin,end0=rpos,qstart=qst,qend=qend,cigar=cg,
            identity=matches/max(span,1) if fields else bam_identity,
            split_large_indel=True,original_alignment=[row['qstart'],row['qend'],row['start0'],row['end0']])
        if fields:
            segment=[fields[0],row['qlen'],qst,qend,row['strand'],row['chrom'],row['chrom_length'],rbegin,rpos,matches,span,row['mapq']]
            piece['paf']='\t'.join(map(str,segment))+f'\ttp:A:P\tcg:Z:{cg}\tcs:Z:'+''.join(x[2] for x in parts)
        result.append(piece)
    for n,op,cs,match in tokens:
        split=op=='N' or (op in 'ID' and n>=minimum)
        if split:
            flush()
        if op in 'MI=X': qpos+=n
        if op in 'MDN=X': rpos+=n
        if split:
            qbegin,rbegin=qpos,rpos;parts=[]
        else:
            parts.append((n,op,cs,match))
    flush()
    return result


def reliable_chain(rows, args):
    rows=[piece for row in rows for piece in split_large_indels(row,args.min_sv_size)]
    eligible = [r for r in rows if r['mapq'] >= args.min_mapq
                and min(r['qend']-r['qstart'], r['end0']-r['start0']) >= args.min_anchor
                and r['identity'] >= args.min_alignment_identity]
    eligible.sort(key=lambda r: (r['qstart'], r['qend']))
    # Maximum covered-query chain. Large query overlaps are competing mappings,
    # not additional sequence. Short repeat hits cannot become terminal nodes.
    scores, previous = [], []
    for i, row in enumerate(eligible):
        score, prev = (row['qend']-row['qstart'])*row['identity'], -1
        for j, other in enumerate(eligible[:i]):
            overlap = max(0, other['qend']-row['qstart'])
            if row['qend'] <= other['qend'] or overlap > args.max_query_overlap:
                continue
            value = scores[j]+(row['qend']-row['qstart']-overlap)*row['identity']
            if value > score:
                score, prev = value, j
        scores.append(score); previous.append(prev)
    if not scores:
        return []
    i, selected = int(np.argmax(scores)), []
    while i >= 0:
        selected.append(eligible[i]); i = previous[i]
    return selected[::-1]


def endpoint(row, first):
    coordinate = row['end0'] if (row['strand'] == '+') == first else row['start0']
    # Side is the retained reference side, independent of query orientation.
    side = 'L' if (row['strand'] == '+') == first else 'R'
    return dict(chrom=row['chrom'], pos0=int(coordinate), side=side)


def find_outer_event(name, rows, query, args):
    chain = reliable_chain(rows, args)
    if len(chain) < 2:
        return None, 'fewer_than_two_reliable_alignments'
    first, last = chain[0], chain[-1]
    a, b = endpoint(first, True), endpoint(last, False)
    inside = lambda ep: ep['chrom'] == query['chrom'] and query['start0'] <= ep['pos0'] < query['end0']
    if not (inside(a) or inside(b)):
        return None, 'neither_outer_breakend_in_query'
    qgap = last['qstart']-first['qend']
    if first['chrom'] == last['chrom'] and first['strand'] == last['strand']:
        rgap = last['start0']-first['end0'] if first['strand'] == '+' else first['start0']-last['end0']
        if rgap==0 and qgap>=args.min_sv_size:
            return None,'pure_insertion_without_new_reference_edge'
        if abs(rgap-qgap) < args.min_sv_size:
            return None, 'collinear_outer_alignments'
    canonical = sorted([a, b], key=lambda e: (e['chrom'],e['pos0'],e['side']))
    return dict(name=name, query_id=query['query_id'], endpoint_a=canonical[0], endpoint_b=canonical[1],
        ordered_endpoints=[a,b], chain=chain, query_gap=qgap,
        minimum_anchor=min(first['qend']-first['qstart'],last['qend']-last['qstart']),
        minimum_mapq=min(first['mapq'],last['mapq']), reliable_alignment_count=len(chain)), 'outer_candidate'


def same_event(a, b, distance):
    for key in ('endpoint_a','endpoint_b'):
        x, y = a[key], b[key]
        if x['chrom'] != y['chrom'] or x['side'] != y['side'] or abs(x['pos0']-y['pos0']) > distance:
            return False
    return True


def pair_event(first, last):
    ends=sorted([endpoint(first,True),endpoint(last,False)],key=lambda e:(e['chrom'],e['pos0'],e['side']))
    return dict(endpoint_a=ends[0],endpoint_b=ends[1])


def junction_support(chain, raw_chains, names, args):
    """Support both outer junctions of a possibly longer assembled path.

A long unitig need not be spanned by one molecule. Each terminal junction must
have independent read support, while overlap layout links its interior.
"""
    targets=[pair_event(chain[0],chain[1]),pair_event(chain[-2],chain[-1])]
    support=[set(),set()]
    for name in names:
        rows=reliable_chain(raw_chains.get(name,[]),args)
        for a,b in zip(rows[:-1],rows[1:]):
            event=pair_event(a,b)
            for i,target in enumerate(targets):
                if same_event(event,target,args.breakpoint_cluster):
                    support[i].add(name)
    return [sorted(s) for s in support]


def cluster_events(events, distance):
    groups = []
    for event in sorted(events, key=lambda e: (-e['minimum_anchor'], -e['minimum_mapq'], e['name'])):
        group = next((g for g in groups if same_event(event, g[0], distance)), None)
        if group is None:
            groups.append([event])
        else:
            group.append(event)
    return groups


def paf_chains(path):
    chains = defaultdict(list)
    with Path(path).open() as f:
        for line in f:
            fields = line.rstrip().split('\t')
            tags = {x[:2]: x[5:] for x in fields[12:]}
            if tags.get('tp') == 'S':
                continue
            qlen, qst, qend = map(int, fields[1:4])
            rlen, st, end, matches, span, mapq = map(int, fields[6:12])
            row = dict(chrom=fields[5], start0=st, end0=end, strand=fields[4],
                qstart=qst, qend=qend, qlen=qlen, mapq=mapq, identity=matches/max(span,1),
                chrom_length=rlen, cigar=tags.get('cg',''), paf=line.rstrip())
            chains[fields[0]].append(row)
    return dict(chains)


def align_sequences(fasta, paf, args, outdir):
    if not Path(fasta).stat().st_size:
        Path(paf).write_text('')
        return {}
    index = Path(args.reference_index) if args.reference_index else outdir/'reference.mmi'
    if not index.exists():
        with (outdir/'reference_index.log').open('w') as log:
            subprocess.run(['minimap2','-x',args.minimap_preset,'-t',str(args.thread),
                '-d',str(index),str(args.reference)], stdout=log, stderr=subprocess.STDOUT, check=True)
    with Path(paf).open('w') as output, Path(str(paf)+'.log').open('w') as log:
        subprocess.run(['minimap2','-x',args.minimap_preset,'-t',str(args.thread),
            '-c','--cs','--secondary=yes','-N','10',str(index),str(fasta)], stdout=output,stderr=log,check=True)
    return paf_chains(paf)


def event_table_row(event):
    a,b = event['endpoint_a'],event['endpoint_b']
    return dict(candidate_id=event.get('candidate_id',''),method=event.get('method',''),
        query_id=event['query_id'],chrom_a=a['chrom'],pos_a0=a['pos0'],side_a=a['side'],
        chrom_b=b['chrom'],pos_b0=b['pos0'],side_b=b['side'],
        representative=event['name'],minimum_anchor=event['minimum_anchor'],
        minimum_mapq=event['minimum_mapq'],alignment_count=event['reliable_alignment_count'],
        raw_support=len(event.get('support_reads',[])),
        outer_pair_raw_support=event.get('outer_pair_raw_support',len(event.get('support_reads',[]))),
        left_junction_support=event.get('left_junction_support',''),
        right_junction_support=event.get('right_junction_support',''),
        bam_support=len(event.get('bam_support_reads',[])),
        olc_input_reads=event.get('olc_input_reads',0),extension_bp=event.get('extension_bp',0),
        status=event.get('status','candidate'),duplicate_of=event.get('duplicate_of',''),
        handoff_type=event.get('handoff_type',''),indel_event_type=event.get('indel_event_type',''),
        indel_span_bp=event.get('indel_span_bp',''),
        support_reads=','.join(event.get('support_reads',[])))


def discover_events(queries, outdir, args):
    sequences_by_query, bam_events, rejections = {}, {}, []
    raw_sequences = {}
    extraction_changed=False
    with pysam.AlignmentFile(args.bam, 'rb') as bam:
        for query in queries:
            qid = query['query_id']
            folder = outdir/qid
            folder.mkdir(exist_ok=True)
            previous=json.loads((folder/'extraction.json').read_text()).get('query',{}) if (folder/'extraction.json').exists() else {}
            same_bounds=all(previous.get(key)==query[key] for key in ('chrom','start0','end0'))
            if same_bounds:
                with pysam.FastxFile(str(folder/'reads.fa')) as fasta:
                    sequences = {r.name:r.sequence for r in fasta}
                chains = json.loads((folder/'bam_chains.json').read_text())
            else:
                for name in ('reads.fa','bam_chains.json','extraction.json'):
                    if (folder/name).is_symlink():
                        (folder/name).unlink()
                sequences, chains = extract_region_reads(bam, query, folder)
                extraction_changed=True
                if (folder/'olc').is_symlink():
                    (folder/'olc').unlink()
                elif (folder/'olc').exists():
                    shutil.rmtree(folder/'olc')
            sequences_by_query[qid] = sequences
            events = []
            for name, rows in chains.items():
                event, reason = find_outer_event(name, rows, query, args)
                if event:
                    events.append(event)
                else:
                    rejections.append(dict(query_id=qid,name=name,status=reason))
                # Realign all split/clip molecules, before applying terminal
                # reliability gates. This can recover misleading BAM chains.
                clipped = any('S' in r.get('cigar','') or 'H' in r.get('cigar','') or
                    any(int(n) >= args.min_sv_size for n in re.findall(r'(\d+)[ID]',r.get('cigar',''))) for r in rows)
                if name in sequences and (len(rows) > 1 or clipped):
                    raw_sequences[name] = sequences[name]
            bam_events[qid] = events
            write_json(folder/'bam_outer_events.json', events)
            logging.info('%s: %d full reads, %d BAM outer candidates', qid,len(sequences),len(events))
    if extraction_changed:
        # A reused query ID can refer to a different adaptive interval. Rebuild
        # the read union and alignments instead of silently using the old scope.
        for name in ('raw_realign.fa','raw_realign.paf','olc_realign.paf'):
            (outdir/name).unlink(missing_ok=True)
    raw_fasta = outdir/'raw_realign.fa'
    if not raw_fasta.exists():
        with raw_fasta.open('w') as f:
            for name, sequence in sorted(raw_sequences.items()):
                f.write(f'>{name}\n{sequence}\n')
    raw_paf = outdir/'raw_realign.paf'
    raw_chains = paf_chains(raw_paf) if raw_paf.exists() else align_sequences(raw_fasta,raw_paf,args,outdir)
    raw_events, raw_groups = {}, {}
    for query in queries:
        qid = query['query_id']
        events = []
        for name in sequences_by_query[qid]:
            if name in raw_chains:
                event, reason = find_outer_event(name,raw_chains[name],query,args)
                if event:
                    events.append(event)
        raw_events[qid] = events
        raw_groups[qid] = cluster_events(events,args.breakpoint_cluster)
        write_json(outdir/qid/'raw_outer_events.json',events)
    candidates = []
    if args.method == 'read':
        for query in queries:
            qid=query['query_id']
            for group in raw_groups[qid]:
                event=dict(group[0],method='read',support_reads=sorted({e['name'] for e in group}))
                event['bam_support_reads'] = sorted({e['name'] for e in bam_events[qid] if same_event(e,event,args.breakpoint_cluster)})
                event['status'] = 'supported' if len(event['support_reads']) >= args.min_support else 'insufficient_raw_support'
                candidates.append(event)
    else:
        config=OLCConfig(min_overlap=args.olc_min_overlap,min_identity=args.olc_min_identity,
            max_overhang=args.olc_max_overhang,band=args.olc_band,max_overlap_indel=args.olc_max_indel,
            k=args.olc_k,window=args.olc_window,max_seed_copies_per_read=args.olc_max_seed_copies)
        assembly_records={}
        records_by_query={}
        pending=[]
        for query in queries:
            qid=query['query_id'];folder=outdir/qid/'olc'
            cached=json.loads((folder/'assembly.json').read_text()) if (folder/'assembly.json').exists() else None
            if cached and cached['config'] == asdict(config):
                records_by_query[qid]=cached['records']
            else:
                pending.append((qid,folder))
        assembly_changed=bool(pending)
        workers=max(1,min(args.olc_jobs,args.thread))
        if workers==1 or len(pending)<=1:
            for qid,folder in pending:
                logging.info('%s: custom OLC on %d reads',qid,len(sequences_by_query[qid]))
                records_by_query[qid]=assemble(sequences_by_query[qid],folder,config)
        else:
            logging.info('Custom OLC: %d query regions, %d workers',len(pending),workers)
            with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as pool:
                futures={pool.submit(assemble,sequences_by_query[qid],folder,config):qid for qid,folder in pending}
                for future in as_completed(futures):
                    qid=futures[future]
                    records_by_query[qid]=future.result()
                    logging.info('%s: custom OLC completed (%d unitigs)',qid,len(records_by_query[qid]))
        with (outdir/'olc_realign.fa').open('w') as output:
            for query in queries:
                qid=query['query_id']; folder=outdir/qid/'olc'
                records=records_by_query[qid]
                records={r['name']:r for r in records}
                with pysam.FastxFile(str(folder/'unitigs.fa')) as fasta:
                    for unitig in fasta:
                        record=records[unitig.name]
                        if record['assembled']:
                            name=qid+'_'+unitig.name
                            assembly_records[name]=dict(record,query_id=qid)
                            output.write(f'>{name}\n{unitig.sequence}\n')
        olc_paf=outdir/'olc_realign.paf'
        olc_chains=paf_chains(olc_paf) if olc_paf.exists() and not assembly_changed else align_sequences(outdir/'olc_realign.fa',olc_paf,args,outdir)
        for query in queries:
            qid=query['query_id']; events=[]
            for name,record in assembly_records.items():
                if record['query_id'] != qid:
                    continue
                event,reason=find_outer_event(name,olc_chains.get(name,[]),query,args)
                if not event:
                    rejections.append(dict(query_id=qid,name=name,status=reason))
                    continue
                event.update(method='olc',olc_input_reads=len(record['input_reads']),
                    extension_bp=record['extension_bp'],assembly_record=record)
                events.append(event)
            for group in cluster_events(events,args.breakpoint_cluster):
                event=group[0]
                exact_support={e['name'] for e in raw_events[qid] if same_event(e,event,args.breakpoint_cluster)}
                terminal_support=junction_support(event['chain'],raw_chains,sequences_by_query[qid],args)
                event['outer_pair_raw_support']=len(exact_support)
                event['left_junction_support'],event['right_junction_support']=map(len,terminal_support)
                event['support_reads']=sorted(set(terminal_support[0])|set(terminal_support[1]))
                event['bam_support_reads']=sorted({e['name'] for e in bam_events[qid] if same_event(e,event,args.breakpoint_cluster)})
                event['status']='supported' if min(map(len,terminal_support)) >= args.min_support else 'insufficient_raw_support'
                candidates.append(event)
    for i,event in enumerate(candidates,1):
        event['candidate_id']=f'R{i:04d}'
    write_json(outdir/'raw_candidates.json',candidates)
    write_tsv(outdir/'read_rejections.tsv',rejections,['query_id','name','status'])
    write_tsv(outdir/'raw_candidates.tsv',[event_table_row(e) for e in candidates],['candidate_id','query_id','status'])
    return candidates


def filter_small_indel_candidates(candidates):
    """Apply the existing 100 kb minimum to standalone same-strand DUP/DELs.

    Keep small alignment pieces during discovery: they can be internal evidence
    for a longer compound NClose. Only the standalone augmentation is filtered.
    """
    for event in candidates:
        if event['status'] != 'supported':
            continue
        chain=event['chain']
        if len(chain) < 2 or len({(row['chrom'],row['strand']) for row in chain}) != 1:
            continue
        a,b=event['endpoint_a'],event['endpoint_b']
        if a['side'] == b['side']:
            continue
        span=abs(a['pos0']-b['pos0'])
        event['indel_span_bp']=span
        if span < VCF_TYPE4_MIN_SPAN:
            event['status']='below_min_indel_span'


def prepare_augmented_handoff(prefix, outdir, candidates, args):
    import nclose_preprocess as pre
    from breakend_graph import load_stage10_input, save_stage10_input
    from nclose_tracking import load_event_catalog, save_event_catalog, build_bnd_event_catalog
    from skype_utils import save_nclose_nodes
    from nclose_tracking import make_indel_candidate, same_indel_candidate, nclose_event_id_by_key

    filter_small_indel_candidates(candidates)
    stage=load_stage10_input(prefix/'01_nclose_data.pkl')
    nodes=list(stage.contig_data)
    ncloses={owner:list(pairs) for owner,pairs in stage.nclose_nodes.items()}
    catalog=load_event_catalog(str(prefix))
    event_ids=nclose_event_id_by_key(catalog)
    indel_representatives=[dict(e,source=event_ids[tuple(e['event_key'])])
                          for e in catalog if e['kind']=='indel']
    def indel_layout(first, last, source):
        # Same breakpoint geometry and DEL/DUP sign as stage 11, including RC.
        a=first[8] if first[4]=='+' else first[7]
        b=last[7] if last[4]=='+' else last[8]
        signed_span=(b-a) if first[4]=='+' else (a-b)
        return make_indel_candidate('front_jump' if signed_span>0 else 'back_jump',
                                    first[5],a,b,source)
    for i,node in enumerate(nodes):
        if node[10]==4 and node[11]==i:
            indel_representatives.append(indel_layout(node,nodes[node[12]],node[0]))
    repeats=pre.import_repeat_data_00(args.repeat_bed) if args.repeat_bed else {}
    censat=pre.import_censat_repeat_data(args.censat_bed)
    chr_lengths={row[5]:int(row[6]) for row in nodes}
    repeat_names=pre.extract_all_repeat_contig(nodes,repeats,pre.CTG_RPTCASE,pre.NON_REPEAT_NOISE_RATIO)
    representatives=defaultdict(list)
    def layout(pair):
        a,b=(nodes[i] for i in pair)
        ai,bi=(a[7],a[8]),(b[7],b[8])
        bucket,stored,same=pre.nclose_compression_layout(('=',a[5]),('=',b[5]),pair[0],pair[1],ai,bi)
        directions=pre.nclose_canonical_directions(nodes,pair,stored,same and ai == bi)
        return bucket,stored,directions
    for owner,pairs in ncloses.items():
        for pair in pairs:
            bucket,stored,directions=layout(tuple(pair))
            representatives[bucket].append(pre.NCloseClusterRepresentative(owner,tuple(pair),stored,directions))
    with (prefix/'paf_file_path.pkl').open('rb') as f:
        original_pafs=pickle.load(f)
    augmented=outdir/'augmented'
    augmented.mkdir(exist_ok=True)
    rescue_paf=augmented/'raw_rescue.paf'
    paf_rows=[]
    added=[]
    for event in candidates:
        if event['status'] != 'supported':
            continue
        start=len(nodes); owner=f"raw_rescue_{args.method}_{event['candidate_id']}"
        chain=event['chain']; end=start+len(chain)-1
        is_type4=len({(row['chrom'],row['strand']) for row in chain})==1
        handoff_type=4 if is_type4 else 1 if len({row['chrom'] for row in chain})>1 else 2
        event['handoff_type']=handoff_type
        chrom_span=defaultdict(int)
        for row in chain:
            chrom_span[(row['chrom'],row['strand'])] += row['end0']-row['start0']
        main_chrom,main_dir=max(chrom_span,key=chrom_span.get)
        current=[]
        for i,row in enumerate(chain):
            node=[owner,row['qlen'],row['qstart'],row['qend'],row['strand'],row['chrom'],row['chrom_length'],
                  row['start0'],row['end0'],row['mapq'],handoff_type,
                  start,end,'0','0','0','0','0','0',main_dir,main_chrom,f'{len(original_pafs)}.{len(paf_rows)+i}']
            label=pre.label_repeat_node([node],repeats,chr_lengths)[0]
            node[16],node[17]=label
            node[18]=pre.label_repeat_node([node],censat,chr_lengths)[0][1]
            current.append(tuple(node))
        if is_type4:
            indel=indel_layout(current[0],current[-1],owner)
            event['indel_event_type']=indel['event_type']
            duplicate=next((r for r in indel_representatives if same_indel_candidate(indel,r)),None)
            if duplicate is not None:
                event.update(status='compressed_duplicate',duplicate_of=duplicate['source'])
                continue
            indel_representatives.append(indel)
        nodes.extend(current)
        pair=(start,end)
        if not is_type4:
            # BNDs retain stage 01's repeat-aware anchor-interval compression.
            local_nodes=[list(n) for n in current]
            for n in local_nodes:
                n[11],n[12]=0,len(local_nodes)-1
            is_repeat=owner in pre.extract_all_repeat_contig(local_nodes,repeats,pre.CTG_RPTCASE,pre.NON_REPEAT_NOISE_RATIO)
            bucket,stored,directions=layout(pair)
            duplicate=None
            for representative in representatives[bucket]:
                limit=pre.ALL_REPEAT_NCLOSE_COMPRESS_LIMIT if is_repeat and representative.contig_name in repeat_names else pre.NCLOSE_COMPRESS_LIMIT
                if pre.nclose_cluster_candidate_matches(nodes,stored,directions,representative,limit):
                    duplicate=representative.contig_name
                    break
            if duplicate:
                del nodes[start:]
                event.update(status='compressed_duplicate',duplicate_of=duplicate)
                continue
            if is_repeat:
                repeat_names.add(owner)
            representatives[bucket].append(pre.NCloseClusterRepresentative(owner,pair,stored,directions))
            ncloses[owner]=[pair]
        for row in chain:
            fields=row['paf'].split('\t'); fields[0]=owner
            # Stage 21 expects cs to be the final PAF tag.
            cs=[x for x in fields[12:] if x.startswith('cs:Z:')]
            tags=['tp:A:P' if x=='tp:A:S' else x for x in fields[12:]
                  if not x.startswith(('cs:Z:','cg:Z:'))]
            paf_rows.append('\t'.join(fields[:12]+tags+cs))
        event.update(status='added',owner=owner,node_pair=list(pair))
        added.append(event)
    write_tsv(outdir/'raw_candidates.tsv',[event_table_row(e) for e in candidates],['candidate_id','query_id','status'])
    write_json(outdir/'raw_candidates.json',candidates)
    if not added:
        return [],[]
    rescue_paf.write_text('\n'.join(paf_rows)+'\n')
    original_pafs.append(str(rescue_paf.resolve()))
    with (augmented/'paf_file_path.pkl').open('wb') as f:
        pickle.dump(original_pafs,f)
    save_stage10_input(augmented/'01_nclose_data.pkl',nodes,ncloses,stage.telo_contig)
    save_nclose_nodes(str(augmented),ncloses)
    ppc=Path(args.ppc_paf) if args.ppc_paf else next(prefix.glob('*.ppc.paf'))
    with (augmented/ppc.name).open('w') as f:
        for node in nodes:
            f.write('\t'.join(map(str,node)).rstrip()+'\n')
    groups=pre.group_nclose_nodes_by_chrom(nodes,ncloses)
    pre.write_nclose_nodes_list(str(augmented/'compressed_nclose_nodes_list.txt'),groups,nodes,repeat_names)
    pre.write_nclose_nodes_index(str(augmented/'nclose_nodes_index.txt'),ncloses,nodes)
    old_keys={tuple(e['event_key']) for e in catalog}
    for event in build_bnd_event_catalog(groups,nodes):
        if tuple(event['event_key']) not in old_keys:
            event['rescue_method']=args.method
            catalog.append(event)
    save_event_catalog(str(augmented),catalog)
    files=[p.name for p in augmented.iterdir() if p.is_file() and p.name != rescue_paf.name]
    return added,files


def build_parser():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prefix',type=Path)
    parser.add_argument('--method',choices=('read','olc'),default='read')
    parser.add_argument('--bam')
    parser.add_argument('--reference')
    parser.add_argument('--reference-index')
    parser.add_argument('--censat-bed')
    parser.add_argument('--repeat-bed',default='')
    parser.add_argument('--ppc-paf')
    parser.add_argument('--outdir',type=Path)
    parser.add_argument('--detect-only',action='store_true')
    parser.add_argument('-t','--thread',type=int,default=8)
    parser.add_argument('--dp-penalty',type=float,default=96.)
    parser.add_argument('--min-segment-bins',type=int,default=8)
    parser.add_argument('--local-bins',type=int,default=20)
    parser.add_argument('--prediction-ratio',type=float,default=.2)
    parser.add_argument('--max-hom',type=float,default=2.5)
    parser.add_argument('--query-padding',type=int,default=200000)
    parser.add_argument('--query-mode',choices=('adaptive','fixed'),default='adaptive',
        help='Query the transition/location envelope, or the original fixed padding')
    parser.add_argument('--query-context-bins',type=int,nargs='+',default=[15,20,30])
    parser.add_argument('--query-flank-bins',type=int,default=5)
    parser.add_argument('--query-model-delta',type=float,default=6.,
        help='Allowed increase in standardized SSE + parameter-count*log(n) for query-location models')
    parser.add_argument('--min-mapq',type=int,default=20)
    parser.add_argument('--min-anchor',type=int,default=500)
    parser.add_argument('--min-alignment-identity',type=float,default=.9)
    parser.add_argument('--max-query-overlap',type=int,default=300)
    parser.add_argument('--min-sv-size',type=int,default=1000)
    parser.add_argument('--min-support',type=int,default=3)
    parser.add_argument('--breakpoint-cluster',type=int,default=200)
    parser.add_argument('--minimap-preset',default='map-hifi')
    parser.add_argument('--olc-min-overlap',type=int,default=1000)
    parser.add_argument('--olc-min-identity',type=float,default=.98)
    parser.add_argument('--olc-max-overhang',type=int,default=250)
    parser.add_argument('--olc-band',type=int,default=128)
    parser.add_argument('--olc-max-indel',type=int,default=50)
    parser.add_argument('--olc-k',type=int,default=19)
    parser.add_argument('--olc-window',type=int,default=40)
    parser.add_argument('--olc-max-seed-copies',type=int,default=2)
    parser.add_argument('--olc-jobs',type=int,default=4)
    return parser


def main(argv=None):
    args=build_parser().parse_args(argv)
    prefix=args.prefix.resolve()
    outdir=(args.outdir or prefix/'24_raw_rescue').resolve()
    outdir.mkdir(parents=True,exist_ok=True)
    summary_path=outdir/'summary.json'
    if summary_path.exists() and not args.detect_only:
        summary=json.loads(summary_path.read_text())
        if summary['method'] != args.method:
            raise ValueError('Use a separate baseline copy to compare rescue methods')
        logging.info('Stage 24 already completed once: %s',summary_path)
        return summary
    started=time.monotonic()
    queries=save_depth_detection(prefix,outdir,args)
    write_json(outdir/'parameters.json',{k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()})
    if args.detect_only:
        return
    if not args.bam or not args.reference or not args.censat_bed:
        raise ValueError('--bam, --reference and --censat-bed are required for raw rescue')
    before=outdir/'before'
    before.mkdir(exist_ok=True)
    for name in ('B.npy','predict_B.npy','weight.npy','23_input.pkl','nclose_report.tsv'):
        if (prefix/name).exists() and not (before/name).exists():
            shutil.copyfile(prefix/name,before/name)
    candidates=discover_events(queries,outdir,args) if queries else []
    added,files=prepare_augmented_handoff(prefix,outdir,candidates,args) if candidates else ([],[])
    summary=dict(method=args.method,query_regions=len(queries),raw_candidates=len(candidates),
        added_count=len(added),added_ids=[e['candidate_id'] for e in added],
        added_type4_count=sum(e['handoff_type']==4 for e in added),
        added_bnd_count=sum(e['handoff_type']!=4 for e in added),
        compressed_duplicates=sum(e['status'] == 'compressed_duplicate' for e in candidates),
        insufficient_support=sum(e['status'] == 'insufficient_raw_support' for e in candidates),
        minimum_indel_span=VCF_TYPE4_MIN_SPAN,
        small_indels_excluded=sum(e['status'] == 'below_min_indel_span' for e in candidates),
        augmented_dir=str(outdir/'augmented'),install_files=files,
        seconds=time.monotonic()-started,rounds=1)
    write_json(summary_path,summary)
    logging.info('Stage 24 result: %s',json.dumps(summary))
    return summary


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO,format='%(asctime)s %(levelname)s:%(message)s')
    main()
