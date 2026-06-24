using XAM
using BioAlignments
using ProgressBars
using DataStructures

SplitHit = Tuple{Int, Bool}
Tv = Vector{SplitHit}
AlnInfo = Tuple{Bool, Int, Int, Int}
RawCountKey = Tuple{Int, Int}
RawCounter = DefaultDict{RawCountKey, Int, Int}
RawAnchor = Tuple{Int, Int, Int, Int, Int, Int, Bool}
ReadHits = Dict{String, DefaultDict{RawCountKey, Tuple{Tv, Tv}}}

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
const RAW_FETCH_MERGE_GAP::Int = 5_000
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

function merge_fetch_intervals(intervals::Vector{Tuple{String, Int, Int}},
                               merge_gap::Int=RAW_FETCH_MERGE_GAP)
    if isempty(intervals)
        return intervals
    end

    sorted_intervals = sort(unique(intervals), by=x -> (x[1], x[2], x[3]))
    merged = Vector{Tuple{String, Int, Int}}()

    curr_chr, curr_st, curr_nd = sorted_intervals[1]
    for (chr, span_st, span_nd) in sorted_intervals[2:end]
        if chr == curr_chr && span_st <= curr_nd + merge_gap
            curr_nd = max(curr_nd, span_nd)
        else
            push!(merged, (curr_chr, curr_st, curr_nd))
            curr_chr, curr_st, curr_nd = chr, span_st, span_nd
        end
    end
    push!(merged, (curr_chr, curr_st, curr_nd))

    return merged
end

function make_progress(desc::String; total=nothing)
    progress = total === nothing ?
               ProgressBar(printing_delay=0.1) :
               ProgressBar(total=Int64(total), printing_delay=0.1)
    set_description(progress, desc)
    return progress
end

function update_progress!(progress, progress_lock=nothing)
    if progress === nothing
        return
    end

    if progress_lock === nothing
        update(progress)
        return
    end

    lock(progress_lock)
    try
        update(progress)
    finally
        unlock(progress_lock)
    end
end

function finish_progress!(progress, progress_lock=nothing)
    if progress === nothing
        return
    end

    if progress_lock === nothing
        update(progress, 0; force_print=true)
        return
    end

    lock(progress_lock)
    try
        update(progress, 0; force_print=true)
    finally
        unlock(progress_lock)
    end
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

    return raw_anchor_data, merge_fetch_intervals(raw_fetch_intervals)
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

function empty_read_hits()
    return ReadHits()
end

function merge_read_hits!(dst::ReadHits, src::ReadHits)
    for (read_name, src_hits_by_key) in src
        dst_hits_by_key = get!(dst, read_name) do
            DefaultDict{RawCountKey, Tuple{Tv, Tv}}(() -> (Tv(), Tv()))
        end

        for (key, (src_hits1, src_hits2)) in pairs(src_hits_by_key)
            dst_hits1, dst_hits2 = dst_hits_by_key[key]
            append!(dst_hits1, src_hits1)
            append!(dst_hits2, src_hits2)
        end
    end

    return dst
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

    progress = nothing
    if is_progress_bar
        progress = make_progress("Analysis bam alignments:")
    end

    while !eof(reader)
        empty!(record)
        read!(reader, record)

        if is_progress_bar
            update_progress!(progress)
        end

        if pass_raw_record(record)
            add_raw_anchor_hits!(read_hits, record, get_alignment_info(record))
        end
    end

    if is_progress_bar
        finish_progress!(progress)
    end
    close(reader)
end

function split_interval_ranges(n_intervals::Int, n_tasks::Int)
    n_chunks = min(n_intervals, n_tasks)
    ranges = Vector{UnitRange{Int}}()

    for chunk_id in 1:n_chunks
        i0 = fld((chunk_id - 1) * n_intervals, n_chunks) + 1
        i1 = fld(chunk_id * n_intervals, n_chunks)
        if i0 <= i1
            push!(ranges, i0:i1)
        end
    end

    return ranges
end

function scan_bam_indexed_range(bam_loc::String, index_loc::String,
                                intervals::Vector{Tuple{String, Int, Int}},
                                interval_range::UnitRange{Int},
                                progress, progress_lock)
    local_hits = empty_read_hits()
    reader = open(BAM.Reader, bam_loc; index=index_loc)

    try
        for i in interval_range
            chr, span_st, span_nd = intervals[i]
            for record in BAM.eachoverlap(reader, chr, span_st:span_nd)
                if pass_raw_record(record)
                    add_raw_anchor_hits!(local_hits, record, get_alignment_info(record))
                end
            end
            update_progress!(progress, progress_lock)
        end
    finally
        close(reader)
    end

    return local_hits
end

function scan_bam_indexed!(read_hits, bam_loc::String, index_loc::String,
                           is_progress_bar::Bool)
    n_intervals = length(raw_fetch_intervals)
    if n_intervals == 0
        return
    end

    progress = nothing
    progress_lock = nothing
    if is_progress_bar
        progress = make_progress(
            "Analysis bam candidate intervals:";
            total=n_intervals
        )
        progress_lock = ReentrantLock()
    end

    n_tasks = min(Threads.nthreads(), n_intervals)
    if n_tasks <= 1
        merge_read_hits!(
            read_hits,
            scan_bam_indexed_range(
                bam_loc, index_loc, raw_fetch_intervals, 1:n_intervals,
                progress, progress_lock
            )
        )
    else
        tasks = [
            Threads.@spawn scan_bam_indexed_range(
                bam_loc, index_loc, raw_fetch_intervals, interval_range,
                progress, progress_lock
            )
            for interval_range in split_interval_ranges(n_intervals, n_tasks)
        ]

        for task in tasks
            merge_read_hits!(read_hits, fetch(task))
        end
    end

    if is_progress_bar
        finish_progress!(progress, progress_lock)
    end
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

    read_hits = empty_read_hits()
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
