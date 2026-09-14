"""Small local overlap-layout-consensus assembler used only by stage 24.

No external assembler or overlapper is invoked. Canonical minimizers propose
read pairs; a banded edit alignment verifies overlaps. The oriented graph is
reduced only through consistent two-edge alternatives, never by tip removal or
bubble collapse. Paths stop at branches. Consensus votes retain the path base
when the evidence is tied. All coordinates are zero-based, half-open.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
import json
import logging
from pathlib import Path
import time

import numpy as np
from numba import njit


@dataclass
class OLCConfig:
    k: int = 19
    window: int = 40
    min_seeds: int = 4
    max_seed_occurrences: int = 256
    max_seed_copies_per_read: int = 2
    min_overlap: int = 1000
    min_identity: float = 0.98
    max_overhang: int = 250
    band: int = 128
    max_overlap_indel: int = 50
    diagonal_bin: int = 100
    transitive_slack: int = 200
    consensus_fraction: float = 0.6
    collapse_contained_layout: bool = True


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGTNacgtn', 'TGCANtgcan'))[::-1]


@njit(cache=True)
def minimizers(sequence, k, window):
    mask = (np.uint64(1) << np.uint64(2*k)) - np.uint64(1)
    forward, reverse, valid = np.uint64(0), np.uint64(0), 0
    hashes = np.full(len(sequence), np.uint64(0xffffffffffffffff), dtype=np.uint64)
    directions = np.zeros(len(sequence), dtype=np.int8)
    for i in range(len(sequence)):
        base = sequence[i]
        code = 0 if base == 65 else 1 if base == 67 else 2 if base == 71 else 3 if base == 84 else -1
        if code < 0:
            forward, reverse, valid = np.uint64(0), np.uint64(0), 0
            continue
        forward = ((forward << np.uint64(2)) | np.uint64(code)) & mask
        reverse = (reverse >> np.uint64(2)) | (np.uint64(3-code) << np.uint64(2*(k-1)))
        valid += 1
        if valid >= k:
            value = min(forward, reverse)
            directions[i-k+1] = 1 if reverse < forward else 0
            value ^= value >> np.uint64(30)
            value *= np.uint64(0xbf58476d1ce4e5b9)
            value ^= value >> np.uint64(27)
            value *= np.uint64(0x94d049bb133111eb)
            hashes[i-k+1] = value ^ (value >> np.uint64(31))
    positions = []
    previous = -1
    for start in range(max(0, len(sequence)-k-window+2)):
        best = start
        for j in range(start+1, start+window):
            if hashes[j] < hashes[best]:
                best = j
        if best != previous and hashes[best] != np.uint64(0xffffffffffffffff):
            positions.append(best)
            previous = best
    return hashes[np.array(positions, dtype=np.int64)], np.array(positions, dtype=np.int64), directions[np.array(positions, dtype=np.int64)]


@njit(cache=True, nogil=True)
def banded_alignment(a, b, band):
    """Global unit-cost alignment. Ops: 0 diagonal, 1 delete a, 2 insert b."""
    n, m = len(a), len(b)
    band = max(band, abs(n-m)+2)
    width, inf = 2*band+1, n+m+1
    trace = np.full((n+1, width), 255, dtype=np.uint8)
    previous = np.full(width, inf, dtype=np.int32)
    for j in range(min(m, band)+1):
        previous[j+band] = j
        trace[0, j+band] = 2
    for i in range(1, n+1):
        current = np.full(width, inf, dtype=np.int32)
        for j in range(max(0, i-band), min(m, i+band)+1):
            col = j-i+band
            value, op = inf, 255
            if j > 0:
                value, op = previous[col] + (a[i-1] != b[j-1]), 0
            if col+1 < width and previous[col+1]+1 < value:
                value, op = previous[col+1]+1, 1
            if j > 0 and col > 0 and current[col-1]+1 < value:
                value, op = current[col-1]+1, 2
            current[col] = value
            trace[i, col] = op
        previous = current
    distance = previous[m-n+band]
    ops = np.empty(n+m, dtype=np.uint8)
    i, j, count = n, m, 0
    while i > 0 or j > 0:
        op = trace[i, j-i+band]
        if op == 255:
            return inf, np.empty(0, dtype=np.uint8)
        ops[count] = op
        count += 1
        if op != 2:
            i -= 1
        if op != 1:
            j -= 1
    return distance, ops[:count][::-1]


@njit(cache=True)
def longest_indel(ops):
    longest = count = 0
    previous = -1
    for value in ops:
        if value and value == previous:
            count += 1
        else:
            count = 1 if value else 0
        longest = max(longest, count)
        previous = value
    return longest


@njit(cache=True, nogil=True)
def verify_seed_chain(a, b, chain, k, band, max_distance, max_allowed_indel):
    """Align between exact anchors instead of a wide band over the whole read."""
    distance, max_indel = 0, 0
    for i in range(1, len(chain)):
        ast, bst = chain[i-1,0]+k, chain[i-1,1]+k
        aend, bend = chain[i,0], chain[i,1]
        equal = aend-ast == bend-bst
        if equal:
            for j in range(aend-ast):
                if a[ast+j] != b[bst+j]:
                    equal = False
                    break
        if equal:
            continue
        if distance+abs((aend-ast)-(bend-bst)) > max_distance:
            return max_distance+1,max_indel
        gap_band = min(band, max(8, abs((aend-ast)-(bend-bst))+8))
        errors, ops = banded_alignment(a[ast:aend], b[bst:bend], gap_band)
        distance += errors
        max_indel = max(max_indel,longest_indel(ops))
        if distance > max_distance or max_indel > max_allowed_indel:
            return distance,max_indel
    return distance,max_indel


def seed_pairs(seeds, lengths, config):
    index = defaultdict(list)
    for i, (hashes, positions, directions) in enumerate(seeds):
        counts = Counter(map(int, hashes))
        for h, p, d in zip(hashes, positions, directions):
            if counts[int(h)] <= config.max_seed_copies_per_read:
                index[int(h)].append((i, int(p), int(d)))
    index = {h: rows for h, rows in index.items()
             if 1 < len(rows) <= config.max_seed_occurrences}
    for i, (hashes, positions, directions) in enumerate(seeds):
        votes = defaultdict(list)
        for h, p, d in zip(hashes, positions, directions):
            for j, q, e in index.get(int(h), ()):
                if j <= i:
                    continue
                flip = int(d) != e
                q = lengths[j]-q-config.k if flip else q
                votes[(j, flip)].append((int(p), q))
        for (j, flip), anchors in votes.items():
            if len(anchors) < config.min_seeds:
                continue
            diagonals = Counter(round((p-q)/config.diagonal_bin) for p, q in anchors)
            best = max(diagonals, key=lambda d: diagonals[d-1]+diagonals[d]+diagonals[d+1])
            anchors = sorted((p, q) for p, q in anchors
                             if abs(p-q-best*config.diagonal_bin) <= config.diagonal_bin*1.5)
            # Keep a collinear chain; repeated/reversed seeds cannot support it.
            chain = []
            for p, q in anchors:
                if not chain or (p > chain[-1][0] and q > chain[-1][1]):
                    chain.append((p, q))
            if len(chain) >= config.min_seeds:
                yield i, j, flip, chain


def verified_overlaps(sequences, config):
    encoded = [np.frombuffer(s.encode('ascii'), dtype=np.uint8) for s in sequences]
    seeds = [minimizers(s, config.k, config.window) for s in encoded]
    oriented = {(i, 0): s for i, s in enumerate(encoded)}
    for i, s in enumerate(sequences):
        oriented[(i, 1)] = np.frombuffer(reverse_complement(s).encode('ascii'), dtype=np.uint8)
    edges, contained, considered = {}, [], 0
    for i, j, flip, chain in seed_pairs(seeds, list(map(len, sequences)), config):
        spaced = []
        for p,q in chain:
            if not spaced or (p >= spaced[-1][0]+config.k and q >= spaced[-1][1]+config.k):
                spaced.append((p,q))
        chain = spaced
        if len(chain) < config.min_seeds:
            continue
        a, b = (i, 0), (j, int(flip))
        ast, bst = chain[0]
        aend, bend = (chain[-1][0]+config.k, chain[-1][1]+config.k)
        an, bn = len(oriented[a]), len(oriented[b])
        span = min(aend-ast, bend-bst)
        if span < config.min_overlap:
            continue
        if abs((aend-ast)-(bend-bst)) > (1-config.min_identity)*max(aend-ast,bend-bst):
            continue
        left, right = ast-bst, (an-aend)-(bn-bend)
        # Dovetails or full contained alignments only, not internal repeats.
        if min(ast, bst) > config.max_overhang or min(an-aend, bn-bend) > config.max_overhang:
            continue
        considered += 1
        distance,max_indel=verify_seed_chain(oriented[a],oriented[b],np.asarray(chain,dtype=np.int64),config.k,config.band,
            int(np.ceil((1-config.min_identity)*max(aend-ast,bend-bst))),config.max_overlap_indel)
        identity = 1-distance/max(aend-ast, bend-bst)
        if identity < config.min_identity or max_indel > config.max_overlap_indel:
            continue
        if left >= 0 and right <= 0 and bn-bend > an-aend:
            pass
        elif left <= 0 and right >= 0 and an-aend > bn-bend:
            a, b, ast, bst, aend, bend = b, a, bst, ast, bend, aend
        else:
            # Preserve contained reads as independent sequences and consensus
            # evidence. They are never deleted from the input population.
            contained.append((a, b, ast-bst, identity))
            continue
        offset = ast-bst
        if offset <= 0 or len(oriented[b])+offset <= len(oriented[a]):
            contained.append((a,b,offset,identity))
            continue
        edge = dict(source=a, target=b, offset=offset, source_end=aend,
                    target_end=bend, identity=identity, overlap=span)
        edges[(a, b)] = edge
        ra, rb = (b[0], 1-b[1]), (a[0], 1-a[1])
        edges[(ra, rb)] = dict(source=ra, target=rb,
            offset=len(oriented[b])+offset-len(oriented[a]),
            source_end=len(oriented[b])-bst, target_end=len(oriented[a])-ast,
            identity=identity, overlap=span)
    return oriented, edges, contained, considered


def reduce_and_layout(edges, nreads, config, contained=(), lengths=()):
    contained_ids=set()
    if config.collapse_contained_layout:
        for a,b,offset,identity in contained:
            if (lengths[a[0]],-a[0]) >= (lengths[b[0]],-b[0]):
                contained_ids.add(b[0])
            else:
                contained_ids.add(a[0])
        edges={key:value for key,value in edges.items()
               if key[0][0] not in contained_ids and key[1][0] not in contained_ids}
    outgoing = defaultdict(dict)
    for (a, b), edge in edges.items():
        outgoing[a][b] = edge
    redundant = set()
    for (a, c), direct in edges.items():
        for b, first in outgoing[a].items():
            second = outgoing.get(b, {}).get(c)
            if second and abs(first['offset']+second['offset']-direct['offset']) <= config.transitive_slack:
                redundant.add((a, c))
                redundant.add(((c[0], 1-c[1]), (a[0], 1-a[1])))
                break
    kept = {k: v for k, v in edges.items() if k not in redundant}
    out, incoming = defaultdict(list), defaultdict(list)
    for a, b in kept:
        out[a].append(b)
        incoming[b].append(a)
    used, paths = set(), []
    def walk(a, b):
        path = [a]
        seen_reads = {a[0]}
        while (a, b) not in used and b[0] not in seen_reads:
            used.add((a, b))
            used.add(((b[0], 1-b[1]), (a[0], 1-a[1])))
            path.append(b)
            seen_reads.add(b[0])
            if len(incoming[b]) != 1 or len(out[b]) != 1:
                break
            a, b = b, out[b][0]
        if len(path) > 1:
            paths.append(path)
    for a, b in sorted(kept):
        if len(incoming[a]) != 1 or len(out[a]) != 1:
            walk(a, b)
    for a, b in sorted(kept):
        walk(a, b)
    used_reads = {i for path in paths for i, _ in path}
    paths.extend([[(i, 0)] for i in range(nreads) if i not in used_reads and i not in contained_ids])
    return kept, paths, len(redundant)


def make_scaffold(path, sequences, edges):
    scaffold = sequences[path[0][0]] if path[0][1] == 0 else reverse_complement(sequences[path[0][0]])
    placements, source_start = [(path[0], 0)], 0
    for a, b in zip(path[:-1], path[1:]):
        edge = edges[(a, b)]
        seq = sequences[b[0]] if b[1] == 0 else reverse_complement(sequences[b[0]])
        splice = source_start+edge['source_end']
        scaffold = scaffold[:splice]+seq[edge['target_end']:]
        source_start = splice-edge['target_end']
        placements.append((b, source_start))
    return scaffold, placements


@njit(cache=True, nogil=True)
def add_consensus_votes(votes, sequence, ops, offset):
    positions,starts,ends=[],[],[]
    p=q=cursor=0
    while cursor < len(ops):
        op=ops[cursor]
        if op == 0:
            base=sequence[q]
            code=0 if base==65 else 1 if base==67 else 2 if base==71 else 3 if base==84 else -1
            if code >= 0:
                votes[offset+p,code]+=1
            p+=1;q+=1
        elif op == 1:
            votes[offset+p,4]+=1;p+=1
        else:
            start=q
            while cursor < len(ops) and ops[cursor] == 2:
                q+=1;cursor+=1
            positions.append(offset+p);starts.append(start);ends.append(q)
            continue
        cursor+=1
    return positions,starts,ends


def consensus(scaffold, placements, oriented, config):
    template = np.frombuffer(scaffold.encode('ascii'), dtype=np.uint8)
    votes = np.zeros((len(template), 5), dtype=np.int32)
    insertions = defaultdict(Counter)
    supported = []
    for read, offset in placements:
        seq = oriented[read]
        lo, qlo = max(0, offset), max(0, -offset)
        n = min(len(template)-lo, len(seq)-qlo)
        if n < config.min_overlap:
            continue
        a, b = template[lo:lo+n], seq[qlo:qlo+n]
        distance, ops = banded_alignment(a, b, config.band)
        if not len(ops) or 1-distance/max(len(a), len(b)) < config.min_identity:
            continue
        supported.append(read[0])
        positions,starts,ends=add_consensus_votes(votes,b,ops,lo)
        for position,start,end in zip(positions,starts,ends):
            insertions[position][bytes(b[start:end])]+=1
    depths=votes.sum(axis=1)
    winners=votes.argmax(axis=1)
    confident=(depths>=2)&(votes[np.arange(len(votes)),winners]/np.maximum(depths,1)>=config.consensus_fraction)
    keep=~(confident&(winners==4))
    bases=template.copy()
    changed=confident&(winners<4)
    bases[changed]=np.frombuffer(b'ACGT',dtype=np.uint8)[winners[changed]]
    result=[];start=0
    for position,counts in sorted(insertions.items()):
        if position >= len(depths) or depths[position] == 0:
            continue
        insertion,count=counts.most_common(1)[0]
        if count>=2 and count/depths[position]>=config.consensus_fraction:
            result.extend([bytes(bases[start:position][keep[start:position]]),insertion])
            start=position
    result.append(bytes(bases[start:][keep[start:]]))
    return b''.join(result).decode('ascii'),sorted(set(supported))


def assemble(reads, outdir, config=None):
    config = config or OLCConfig()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    names = sorted(reads)
    sequences = [reads[name].upper() for name in names]
    oriented, edges, contained, considered = verified_overlaps(sequences, config)
    overlap_seconds=time.monotonic()-started
    logging.info('Local OLC: %d reads, %d oriented overlaps, %.1fs for overlaps',len(names),len(edges),overlap_seconds)
    kept, paths, reduced = reduce_and_layout(edges,len(names),config,contained,list(map(len,sequences)))
    layout_seconds=time.monotonic()-started-overlap_seconds
    logging.info('Local OLC: %d paths after reduction and containment layout',len(paths))
    layout_ids={i for path in paths for i,_ in path}
    extra = defaultdict(list)
    for a, b, offset, _ in contained:
        if (len(oriented[a]),-a[0]) >= (len(oriented[b]),-b[0]):
            extra[a].append((b, offset))
            extra[(a[0], 1-a[1])].append(((b[0], 1-b[1]), len(oriented[a])-offset-len(oriented[b])))
        else:
            extra[b].append((a, -offset))
            extra[(b[0], 1-b[1])].append(((a[0], 1-a[1]), len(oriented[b])+offset-len(oriented[a])))
    records = []
    with (outdir/'unitigs.fa').open('w') as fasta:
        for idx, path in enumerate(paths, 1):
            scaffold, placements = make_scaffold(path, sequences, kept)
            seen = {r[0] for r, _ in placements}
            cursor=0
            while cursor < len(placements):
                read,offset=placements[cursor]
                cursor+=1
                for additional, delta in extra[read]:
                    if additional[0] not in seen:
                        placements.append((additional, offset+delta))
                        seen.add(additional[0])
            sequence, supporters = consensus(scaffold, placements, oriented, config)
            name = f'olc_{idx:06d}'
            record = dict(name=name, length=len(sequence), layout_reads=[names[i] for i, _ in path],
                input_reads=[names[i] for i in sorted(seen)], consensus_reads=[names[i] for i in supporters],
                extension_bp=len(sequence)-max(len(sequences[i]) for i in seen),
                assembled=len(set(supporters)) > 1, path_orientations=[d for _, d in path])
            records.append(record)
            fasta.write(f'>{name}\n{sequence}\n')
    with (outdir/'overlaps.tsv').open('w') as handle:
        handle.write('source\tsource_strand\ttarget\ttarget_strand\toffset\toverlap\tidentity\tremoval_reason\n')
        for (a, b), edge in sorted(edges.items()):
            reason='contained_layout' if a[0] not in layout_ids or b[0] not in layout_ids else 'transitive' if (a,b) not in kept else ''
            handle.write(f"{names[a[0]]}\t{a[1]}\t{names[b[0]]}\t{b[1]}\t{edge['offset']}\t{edge['overlap']}\t{edge['identity']:.6f}\t{reason}\n")
    summary = dict(config=asdict(config), reads=len(reads), verified_pairs=considered,
        directed_overlaps=len(edges), transitive_removed=reduced, contained_alignments=len(contained),
        contained_layout_reads=len(names)-len(layout_ids),
        unitigs=len(records), assembled_unitigs=sum(r['assembled'] for r in records),
        overlap_seconds=overlap_seconds,layout_seconds=layout_seconds,
        seconds=time.monotonic()-started, records=records)
    (outdir/'assembly.json').write_text(json.dumps(summary, indent=2)+'\n')
    return records
