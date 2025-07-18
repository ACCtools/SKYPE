using XAM
using BioAlignments
using ProgressMeter
using DataStructures

Tv = Vector{Tuple{Int, Bool}}
AlnInfo = Tuple{Int, Int, Bool, Int, Int, Int}
Counter = DefaultDict{Int, Int, Int}

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

const DIR2INT = Dict(
    (true,  false) => 1,
    (false, true)  => 2,
    (true,  true)  => 3,
    (false, false) => 4
)

const NCLOSE_COMPRESS_LIMIT::Int = 100 * 1e3

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
    # destructure into op-codes and lengths
    cigar_ops, cigar_lens = XAM.BAM.cigar_rle(rec)

    # count front and tail clips (hard/soft)
    front_clip = 0
    if !isempty(cigar_ops) && (cigar_ops[1] == OP_HARD_CLIP || cigar_ops[1] == OP_SOFT_CLIP)
        front_clip = cigar_lens[1]
        if length(cigar_ops) > 1 && cigar_ops[2] == OP_SOFT_CLIP
            front_clip += cigar_lens[2]
        end
    end
    tail_clip = 0
    n = length(cigar_ops)
    if n > 0 && (cigar_ops[n] == OP_HARD_CLIP || cigar_ops[n] == OP_SOFT_CLIP)
        tail_clip = cigar_lens[n]
        if n > 1 && cigar_ops[n-1] == OP_SOFT_CLIP
            tail_clip += cigar_lens[n-1]
        end
    end

    # compute full query length (including soft/hard clips and insertions)
    qlen = sum(len for (op, len) in zip(cigar_ops, cigar_lens) if op in QRY_OPS)

    # determine 1-based aligned span on query, accounting for strand
    aln_dir_flag = XAM.BAM.ispositivestrand(rec)
    if aln_dir_flag
        qstart = front_clip + 1
        qend   = qlen - tail_clip
    else
        qstart = tail_clip + 1
        qend   = qlen - front_clip
    end

    # compute reference-aligned length (M, =, X, D, N)
    rlen = sum(len for (op, len) in zip(cigar_ops, cigar_lens) if op in REF_OPS)

    # 1-based reference start and end (inclusive)
    rstart = XAM.BAM.position(rec) + 1        # POS is 0-based in SAM
    rend   = rstart + rlen - 1

    # strand symbol and reference name
    aln_dir = aln_dir_flag
    rname   = chr2int[XAM.BAM.refname(rec)]

    return qstart, qend, aln_dir, rname, rstart, rend
end

function analyze_alignments!(read_name::String, align_infos::Vector{AlnInfo},
                             nclose_cnt_vec::Vector{Counter}, wa_nclose_cnt::Counter)

    tar_data = DefaultDict{Int, Tuple{Tv, Tv}}(() -> (Tv(), Tv()))

    for (qst, qnd, aln, rid, rst, rnd) in align_infos
        qlen = qnd - qst + 1
        search_range = NCLOSE_COMPRESS_LIMIT + qlen

        cord_data = chr_cord_data[rid]
        cord_info_data = chr_cord_info_data[rid]

        i = searchsortedfirst(cord_data, rst - search_range)
        j = searchsortedlast(cord_data, rnd + search_range)

        for idx in i:j
            nidx, pairid, isfront = cord_info_data[idx] 

            push!(tar_data[nidx][pairid], (qst, isfront == aln))
        end
    end

    nclose_total_set = DefaultDict{Int, Set{Int}}(() -> Set{Int}())
    for (nidx, (qry_data_vec1, qry_data_vec2)) in pairs(tar_data)
        result = Set{Tuple{Bool,Bool}}()

        if !isempty(qry_data_vec1) && !isempty(qry_data_vec2)
            sort!(qry_data_vec1, by = x -> x[1])
            sort!(qry_data_vec2, by = x -> x[1])

            i = 1
            j = 1
            n1 = length(qry_data_vec1)
            n2 = length(qry_data_vec2)

            while i <= n1 && j <= n2
                q1cord, q1dir = qry_data_vec1[i]
                q2cord, q2dir = qry_data_vec2[j]

                if q1cord < q2cord
                    push!(result, (q1dir, q2dir))
                    i += 1
                elseif q2cord < q1cord
                    push!(result, (q2dir, q1dir))
                    j += 1
                else
                    i += 1
                    j += 1
                end
            end
        end
        
        if length(result) > 0
            ncl = idx2nclose[nidx]
            if length(result) == 1
                push!(nclose_total_set[ncl], DIR2INT[first(result)])
            else
                push!(nclose_total_set[ncl], 0)
            end
        end
    end

    for (ncl, type_set) in nclose_total_set
        if length(type_set) == 1
            d = first(type_set)
            if d == 0
                wa_nclose_cnt[ncl] += 1
            else
                nclose_cnt_vec[d][ncl] += 1
            end
        else
            wa_nclose_cnt[ncl] += 1
        end
    end
end

function anal_bam(bam_loc::String, fai_loc::String, nclose_cord_list::Vector{Vector{Any}}, is_progress_bar::Bool)
    global chr2int, chr_cord_data, chr_cord_info_data, idx2nclose
    chr2int = get_chr2int(fai_loc)

    chr_cord_data = Vector{Vector{Int}}([Vector{Int}() for _ in 1:length(chr2int)])
    chr_cord_info_data = Vector{Vector{Tuple{Int, Int, Bool}}}([Vector{Tuple{Int, Bool}}() for _ in 1:length(chr2int)])

    idx2nclose = Dict{Int, Int}()
    for (i, (chr1, cord1, dir1, chr2, cord2, dir2, num)) in enumerate(nclose_cord_list)
        for (chr, cord, dir, is_front) in [(chr1, cord1, dir1, true), (chr2, cord2, dir2, false)]
            dir = string(dir)
            is_front_dir = dir == "+"
            ci = chr2int[chr]

            push!(chr_cord_data[ci], cord)
            push!(chr_cord_info_data[ci], (i, is_front ? 1 : 2, is_front_dir == is_front))
        end
        
        idx2nclose[i] = num
    end

    for i in eachindex(chr_cord_data)
        si = sortperm(chr_cord_data[i])

        chr_cord_data[i] = chr_cord_data[i][si]
        chr_cord_info_data[i] = chr_cord_info_data[i][si]
    end

    current_read_name = nothing
    align_infos = Vector{AlnInfo}()

    reader = open(BAM.Reader, bam_loc)
    record = BAM.Record()

    nclose_cnt_vec = [Counter(0) for _ in 1:4]
    wa_nclose_cnt = Counter(0)

    if is_progress_bar
        p = ProgressUnknown(desc="Analysis bam alignments:", dt=0.1, showspeed=true)
    end

    while !eof(reader)
        empty!(record)
        read!(reader, record)
        
        if is_progress_bar
            next!(p)
        end
        read_name = XAM.BAM.tempname(record)

        if read_name != current_read_name
            if current_read_name !== nothing && length(align_infos) >= 2
                analyze_alignments!(current_read_name, align_infos, nclose_cnt_vec, wa_nclose_cnt)
            end

            current_read_name = read_name
            empty!(align_infos)
        end

        if XAM.ismapped(record)
            push!(align_infos, get_alignment_info(record))
        end
    end

    if is_progress_bar
        println()
    end
    
    if current_read_name !== nothing && length(align_infos) >= 2
        analyze_alignments!(current_read_name, align_infos, nclose_cnt_vec, wa_nclose_cnt)
    end

    return [[(k, v) for (k, v) in nclose_cnt_vec[i]] for i in 1:4], [(k, v) for (k, v) in wa_nclose_cnt]
end