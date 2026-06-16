using XAM
using BioAlignments
using ProgressMeter
using DataStructures

SplitHit = Tuple{Int, Bool}
Tv = Vector{SplitHit}
AlnInfo = Tuple{Bool, Int, Int, Int}
RawCountKey = Tuple{Int, Int}
RawCounter = DefaultDict{RawCountKey, Int, Int}
RawAnchor = Tuple{Int, Int, Int, Int, Int, Int, Bool}

const QRY_OPS = (
    OP_MATCH,        # M
    OP_SEQ_MATCH,    # =
    OP_SEQ_MISMATCH, # X
    OP_INSERT,       # I
    OP_SOFT_CLIP,    # S
    OP_HARD_CLIP     # H
)

const REF_OPS = (
    OP_MATCH,        # M
    OP_SEQ_MATCH,    # =
    OP_SEQ_MISMATCH, # X
    OP_DELETE,       # D
    OP_SKIP          # N
)

const RAW_TRANSLOCATION_WINDOW::Int = 5 * 1e3
const PANDEPTH_EXCLUDE_FLAGS::UInt16 = UInt16(1796)
const PANDEPTH_MIN_MAPQ::UInt8 = UInt8(0)

function get_chr2int(fai_loc::String)
    chr2int = Dict{String, Int}()
    chr2intnum = 1

    open(fai_loc, "r") do f
        for line in eachline(f)
            chr = split(line, '\t')[1]
            chr2int[chr] = chr2intnum
            chr2intnum += 1
        end
    end

    return chr2int
end

function get_alignment_info(rec::XAM.BAM.Record)
    cigar_ops, cigar_lens = XAM.BAM.cigar_rle(rec)

    # compute reference-aligned length (M, =, X, D, N)
    rlen = sum(len for (op, len) in zip(cigar_ops, cigar_lens) if op in REF_OPS)

    # XAM.BAM.position is already 1-based; keep an inclusive reference span.
    rstart = XAM.BAM.position(rec)
    rend   = rstart + rlen - 1

    # strand symbol and reference name
    aln_dir = XAM.BAM.ispositivestrand(rec)
    rname   = chr2int[XAM.BAM.refname(rec)]

    return aln_dir, rname, rstart, rend
end

function pass_raw_record(rec::XAM.BAM.Record)
    return XAM.ismapped(rec) &&
           (XAM.flags(rec) & PANDEPTH_EXCLUDE_FLAGS) == 0 &&
           XAM.BAM.mappingquality(rec) >= PANDEPTH_MIN_MAPQ
end

function query_pos_at_ref(rec::XAM.BAM.Record, coord::Int, aln::Bool)
    cigar_ops, cigar_lens = XAM.BAM.cigar_rle(rec)
    qlen = sum(len for (op, len) in zip(cigar_ops, cigar_lens) if op in QRY_OPS)
    ref_pos = XAM.BAM.position(rec)
    query_pos = 1

    for (op, len) in zip(cigar_ops, cigar_lens)
        consumes_ref = op in REF_OPS
        consumes_query = op in QRY_OPS

        if consumes_ref
            ref_end = ref_pos + len - 1
            if ref_pos <= coord <= ref_end
                if !consumes_query
                    return nothing
                end
                qpos = query_pos + (coord - ref_pos)
                return aln ? qpos : (qlen - qpos + 1)
            end
            ref_pos = ref_end + 1
        end

        if consumes_query
            query_pos += len
        end
    end

    return nothing
end

function prepare_raw_translocation(raw_junction_spans::Vector{Vector{Any}})
    raw_anchor_data = Vector{Vector{RawAnchor}}([RawAnchor[] for _ in 1:length(chr2int)])
    raw_fetch_intervals = Vector{Tuple{String, Int, Int}}()

    for row in raw_junction_spans
        pair_id = Int(row[1])
        count_idx = Int(row[2])
        anchor_idx = Int(row[3])
        chr = String(row[4])
        span_st = Int(row[5])
        span_nd = Int(row[6])
        point_coord = Int(row[7])
        expected_positive = Bool(row[8])

        if !haskey(chr2int, chr) || anchor_idx < 1 || anchor_idx > 2
            continue
        end
        if span_nd < span_st
            span_st, span_nd = span_nd, span_st
        end
        if span_nd < span_st
            continue
        end

        push!(
            raw_anchor_data[chr2int[chr]],
            (span_st, span_nd, point_coord, pair_id, count_idx, anchor_idx, expected_positive)
        )
        push!(raw_fetch_intervals, (chr, span_st, span_nd))
    end

    return raw_anchor_data, unique(raw_fetch_intervals)
end

function add_raw_anchor_hits!(read_hits, rec::XAM.BAM.Record, info::AlnInfo)
    aln, rid, rst, rnd = info
    if rid > length(raw_anchor_data)
        return
    end

    local_hits = Vector{Tuple{Int, Int, Int, Int, Bool}}()
    for (_span_st, _span_nd, point_coord, pair_id, count_idx, anchor_idx, expected_positive) in raw_anchor_data[rid]
        if rst <= point_coord <= rnd
            qpos = query_pos_at_ref(rec, point_coord, aln)
            if qpos === nothing
                continue
            end
            push!(local_hits, (pair_id, count_idx, anchor_idx, qpos, aln == expected_positive))
        end
    end

    if isempty(local_hits)
        return
    end

    read_name = XAM.BAM.tempname(rec)
    hits_by_key = get!(read_hits, read_name) do
        DefaultDict{RawCountKey, Tuple{Tv, Tv}}(() -> (Tv(), Tv()))
    end
    for (pair_id, count_idx, anchor_idx, qpos, matches_expected) in local_hits
        push!(hits_by_key[(pair_id, count_idx)][anchor_idx], (qpos, matches_expected))
    end
end

function point_pair_is_supported!(hits1::Tv, hits2::Tv)
    if isempty(hits1) || isempty(hits2)
        return false
    end

    for (q1, m1) in hits1
        for (q2, m2) in hits2
            if m1 && m2 && q1 <= q2
                return true
            end
            if !m1 && !m2 && q2 <= q1
                return true
            end
        end
    end

    return false
end

function count_raw_anchor_hits(read_hits)
    raw_count = RawCounter(0)
    for (_, hits_by_key) in read_hits
        for (key, (hits1, hits2)) in pairs(hits_by_key)
            if point_pair_is_supported!(hits1, hits2)
                raw_count[key] += 1
            end
        end
    end
    return raw_count
end

function find_bam_index(bam_loc::String)
    candidates = String[bam_loc * ".bai"]
    if endswith(bam_loc, ".bam")
        push!(candidates, bam_loc[1:end-4] * ".bai")
    end
    for candidate in candidates
        if isfile(candidate)
            return candidate
        end
    end
    return nothing
end

function scan_bam_all!(read_hits, bam_loc::String, is_progress_bar::Bool)
    reader = open(BAM.Reader, bam_loc)
    record = BAM.Record()

    if is_progress_bar
        p = ProgressUnknown(desc="Analysis bam alignments:", dt=0.1, showspeed=true)
    end

    while !eof(reader)
        empty!(record)
        read!(reader, record)

        if is_progress_bar
            next!(p)
        end

        if pass_raw_record(record)
            add_raw_anchor_hits!(read_hits, record, get_alignment_info(record))
        end
    end

    if is_progress_bar
        println()
    end
    close(reader)
end

function scan_bam_indexed!(read_hits, bam_loc::String, index_loc::String,
                           is_progress_bar::Bool)
    reader = open(BAM.Reader, bam_loc; index=index_loc)

    if is_progress_bar
        p = Progress(length(raw_fetch_intervals), desc="Analysis bam candidate intervals:", dt=0.1)
    end

    for (chr, span_st, span_nd) in raw_fetch_intervals
        for record in BAM.eachoverlap(reader, chr, span_st:span_nd)
            if pass_raw_record(record)
                add_raw_anchor_hits!(read_hits, record, get_alignment_info(record))
            end
        end
        if is_progress_bar
            next!(p)
        end
    end

    close(reader)
end

function anal_bam(bam_loc::String, fai_loc::String, is_progress_bar::Bool,
                  raw_junction_spans::Vector{Vector{Any}}=Vector{Vector{Any}}())
    global chr2int
    global raw_anchor_data, raw_fetch_intervals
    chr2int = get_chr2int(fai_loc)

    raw_anchor_data, raw_fetch_intervals = prepare_raw_translocation(raw_junction_spans)
    if all(isempty, raw_anchor_data)
        return Tuple{Int, Int, Int}[]
    end

    read_hits = Dict{String, DefaultDict{RawCountKey, Tuple{Tv, Tv}}}()
    index_loc = find_bam_index(bam_loc)
    if index_loc === nothing
        scan_bam_all!(read_hits, bam_loc, is_progress_bar)
    else
        scan_bam_indexed!(read_hits, bam_loc, index_loc, is_progress_bar)
    end

    raw_translocation_count = count_raw_anchor_hits(read_hits)
    raw_count_list = [(k[1], k[2], v) for (k, v) in raw_translocation_count]
    return raw_count_list
end
