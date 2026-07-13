import os
import sys
import pickle as pkl
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from skype_utils import *

import re
import logging
import argparse
import shutil

# logging 설정(레벨/포맷)은 skype_utils 에서 중앙 관리한다 (LOG_LEVEL).
logging.info("11_Ref_Outlier_Contig_Modify start")

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

CHUKJI_LIMIT = 100*K
INDEL_MERGE_TOLERANCE = 10*K

def import_origin_data(file_path : list) -> list :
    contig_data = []
    int_induce_idx = [1, 2, 3, 6, 7, 8, 9, 10, 11]
    idx = 0
    with open(file_path, "r") as paf_file:
        for curr_contig in paf_file:
            curr_contig = curr_contig.rstrip()
            temp_list = curr_contig.split("\t")
            for i in int_induce_idx:
                temp_list[i] = int(temp_list[i])
            contig_data.append(temp_list)
    return contig_data


def import_data(file_path : str) -> list :
    paf_file = open(file_path, "r")
    contig_data = []
    for curr_contig in paf_file:
        temp_list = curr_contig.rstrip().split("\t")
        int_induce_idx = [CTG_LEN, CTG_STR, CTG_END, \
                          CHR_LEN, CHR_STR, CHR_END, \
                          CTG_MAPQ, CTG_TYP, CTG_STRND, CTG_ENDND,]
        for i in int_induce_idx:
            temp_list[i] = int(temp_list[i])
        contig_data.append(tuple(temp_list))
    paf_file.close()
    return contig_data

def find_chr_len(file_path : str) -> dict:
    chr_data_file = open(file_path, "r")
    chr_len = {}
    for curr_data in chr_data_file:
        curr_data = curr_data.split("\t")
        chr_len[curr_data[0]] = int(curr_data[1])
    chr_data_file.close()
    return chr_len


def get_breakend_coord(contig, side_idx):
    nclose_loc = side_idx == 0
    nclose_dir = contig[CTG_DIR] == '+'
    return contig[CHR_STR if nclose_dir ^ nclose_loc else CHR_END]


def calculate_single_contig_ref_ratio(contig_data : list):
    total_ref_len = 0
    for node in contig_data:
        total_ref_len += node[CHR_END] - node[CHR_STR]

    if len(contig_data) == 1:
        estimated_ref_len = contig_data[0][CHR_END] - contig_data[0][CHR_STR]
        ref_st_ed = (contig_data[0][CHR_STR], contig_data[0][CHR_END])
        return estimated_ref_len/total_ref_len, ref_st_ed

    breakend_st = get_breakend_coord(contig_data[0], 0)
    breakend_ed = get_breakend_coord(contig_data[-1], 1)

    if contig_data[0][CTG_DIR] == '+':
        estimated_ref_len = breakend_ed - breakend_st
        ref_st_ed = (breakend_st, breakend_ed)
    else:
        estimated_ref_len = breakend_st - breakend_ed
        ref_st_ed = (breakend_ed, breakend_st)

    return estimated_ref_len/total_ref_len, ref_st_ed

def distance_checker(node_a : tuple, node_b : tuple) -> int :
    if max(int(node_a[0]), int(node_b[0])) < min(int(node_a[1]), int(node_b[1])):
        return 0   
    else:
        return min(abs(int(node_b[0]) - int(node_a[1])), abs(int(node_b[1]) - int(node_a[0])))
    

def inclusive_checker(tuple_a : tuple, tuple_b : tuple) -> bool :
    if int(tuple_a[0]) <= int(tuple_b[0]) and int(tuple_b[1]) <= int(tuple_a[1]):
        return True
    else:
        return False

def make_indel_candidate(event_type, chrom, ref_a, ref_b, source):
    st, nd = sorted((int(ref_a), int(ref_b)))
    return {
        'event_type': event_type,
        'chrom': chrom,
        'st': st,
        'nd': nd,
        'source': source,
    }

def same_indel_candidate(a, b, tolerance=INDEL_MERGE_TOLERANCE):
    return (
        a['event_type'] == b['event_type'] and
        a['chrom'] == b['chrom'] and
        abs(a['st'] - b['st']) <= tolerance and
        abs(a['nd'] - b['nd']) <= tolerance
    )
    
def cs_to_cigar(cs_tag: str) -> str:
   """
   cs:Z: 태그의 값(접두어 "cs:Z:" 제외)을 받아서 표준 CIGAR 문자열로 변환합니다.

   변환 규칙:
     - ":<number>" : 매칭(=, M) 연산. 이전 연산이 'M'이면 길이를 누적하고,
                    그렇지 않으면 이전 연산을 플러시한 후 'M'으로 시작합니다.
     - "+<seq>"    : 삽입(Insertion). 이전에 누적된 연산이 있으면 먼저 플러시하고,
                    삽입 길이는 해당 문자열의 길이로 계산합니다.
     - "-<seq>"    : 결실(Deletion). 위와 유사하게 처리합니다.
     - "*<base1><base2>": 치환(Mismatch). 연속되는 치환은 길이를 누적합니다.
   """
   pattern = re.compile(r"(:\d+|\*[a-z]{2}|[+\-][A-Za-z]+)")
   cigar = ""
   last_op = 'M'
   last_len = 0

   for match in pattern.finditer(cs_tag):
       part = match.group(0)
       op = part[0]
       if op == ":":
           # 예: ":10" → 매칭 10개
           length = int(part[1:])
           if last_op == 'M':
               last_len += length
           else:
               if last_len > 0:
                   cigar += f"{last_len}{last_op}"
               last_op = 'M'
               last_len = length
       elif op == "-":
           # 예: "-acgt" → deletion, 길이는 4
           length = len(part[1:])
           if last_len > 0:
               cigar += f"{last_len}{last_op}"
           cigar += f"{length}D"
           last_len = 0
           last_op = 'M'
       elif op == "+":
           # 예: "+ac" → insertion, 길이는 2
           length = len(part[1:])
           if last_len > 0:
               cigar += f"{last_len}{last_op}"
           cigar += f"{length}I"
           last_len = 0
           last_op = 'M'
       elif op == "*":
           # 예: "*at" : mismatch (치환)
           if last_op == 'X':
               last_len += 1
           else:
               if last_len > 0:
                   cigar += f"{last_len}{last_op}"
               last_op = 'X'
               last_len = 1
       else:
           # 정의되지 않은 연산은 무시
           continue

   if last_len > 0:
       cigar += f"{last_len}{last_op}"

   return cigar


def node_to_paf_row(node):
    ref_st = min(int(node[CHR_STR]), int(node[CHR_END]))
    ref_nd = max(int(node[CHR_STR]), int(node[CHR_END]))
    N = max(1, ref_nd - ref_st)
    return [
        node[CTG_NAM], int(node[CTG_LEN]), int(node[CTG_STR]), int(node[CTG_END]),
        node[CTG_DIR], node[CHR_NAM], int(node[CHR_LEN]), ref_st, ref_nd,
        N, N, int(node[CTG_MAPQ]), "tp:A:P", "cs:Z:" + f":{N}"
    ]


def node_original_or_synthetic_paf_row(node):
    global_idx = str(node[CTG_GLOBALIDX])
    try:
        paf_idx_str, row_idx_str = global_idx.split(".", 1)
        paf_idx = int(paf_idx_str)
        row_idx = int(row_idx_str)
        if 0 <= paf_idx < len(paf_file) and 0 <= row_idx < len(paf_file[paf_idx]):
            return paf_file[paf_idx][row_idx]
    except (ValueError, IndexError):
        pass
    return node_to_paf_row(node)


def row_with_cigar(row):
    row = list(row)
    cs_tag = next((str(item) for item in reversed(row) if str(item).startswith("cs:Z:")), None)
    if cs_tag is None:
        match_len = max(1, int(row[8]) - int(row[7]))
        cs_tag = "cs:Z:" + f":{match_len}"
        row.append(cs_tag)
    return row + ["cg:Z:" + cs_to_cigar(cs_tag[5:])]


def write_rows_as_paf(path, rows):
    with open(path, "wt") as f:
        for row in rows:
            for value in row_with_cigar(row):
                print(value, end="\t", file=f)
            print("", file=f)


def write_empty_paf(path):
    open(path, "wt").close()


def make_base_paf_row(chrom, ref_st, ref_nd):
    ref_st, ref_nd = sorted((int(ref_st), int(ref_nd)))
    N = max(1, ref_nd - ref_st)
    return [
        "base_contig_1", N, 0, N,
        "+", chrom, int(chr_data[chrom]), ref_st, ref_nd,
        N, N, 0, "tp:A:P", "cs:Z:" + f":{N}"
    ]


def write_base_paf(path, chrom, ref_st, ref_nd):
    write_rows_as_paf(path, [make_base_paf_row(chrom, ref_st, ref_nd)])


def classify_vcf_bnd_type4_pair(nclose_key, node_a_idx, node_b_idx):
    node_a = contig_data[node_a_idx]
    node_b = contig_data[node_b_idx]
    if not str(node_a[CTG_NAM]).startswith("vcf_bnd_"):
        return None
    if node_a[CHR_NAM] != node_b[CHR_NAM]:
        return None

    pos_a = int(get_breakend_coord(node_a, 0))
    pos_b = int(get_breakend_coord(node_b, 1))
    st, nd = sorted((pos_a, pos_b))
    span = nd - st
    if span < CHUKJI_LIMIT:
        return None

    dir_a = node_a[CTG_DIR]
    dir_b = node_b[CTG_DIR]
    if dir_a == "+" and dir_b == "+":
        event_type = "front_jump" if pos_a < pos_b else "back_jump"
    elif dir_a == "-" and dir_b == "-":
        event_type = "back_jump" if pos_a < pos_b else "front_jump"
    else:
        return None

    indel_kind = "deletion" if event_type == "front_jump" else "insertion"
    return {
        "event_id": str(nclose_key),
        "vcf_id": str(nclose_key).replace("vcf_bnd_", "", 1),
        "svtype": "BND",
        "event_type": event_type,
        "indel_kind": indel_kind,
        "chrom": node_a[CHR_NAM],
        "st": st,
        "nd": nd,
        "svlen": span,
        "source": f"VCF_BND_{nclose_key}",
        "nodes": (int(node_a_idx), int(node_b_idx)),
    }


def collect_vcf_bnd_type4_events():
    nclose_chunk_path = f"{args.prefix}/nclose_chunk_data.pkl"
    if not os.path.isfile(nclose_chunk_path):
        return []
    with open(nclose_chunk_path, "rb") as f:
        nclose_nodes, _, _ = pkl.load(f)

    events = []
    seen_spans = set()
    for nclose_key, pair_list in nclose_nodes.items():
        for node_a_idx, node_b_idx in pair_list:
            event = classify_vcf_bnd_type4_pair(nclose_key, node_a_idx, node_b_idx)
            if event is None:
                continue
            span_key = (
                event["event_type"],
                event["chrom"],
                event["st"],
                event["nd"],
            )
            if span_key in seen_spans:
                continue
            seen_spans.add(span_key)
            events.append(event)
    return events



parser = argparse.ArgumentParser(description="Find reference depth of reverse contig")

parser.add_argument("reference_fai_path", 
                        help="Path to the chromosome information file.")

parser.add_argument("ppc_paf_file_path", 
                        help="Path to the preprocessed PAF file.")

parser.add_argument("prefix", 
                    help="Pefix for pipeline")

parser.add_argument("--alt", 
                        help="Path to an alternative PAF file (optional).")

args = parser.parse_args()

PAF_FILE_PATH_PKL = f'{args.prefix}/paf_file_path.pkl'
with open(PAF_FILE_PATH_PKL, 'rb') as f:
    PAF_FILE_PATH = pkl.load(f)

paf_file = [import_origin_data(paf_loc) for paf_loc in PAF_FILE_PATH]

PREPROCESSED_PAF_FILE_PATH = args.ppc_paf_file_path
TYPE_4_VECTOR_PATH = f"{args.prefix}/11_ref_ratio_outliers"
CHROMOSOME_INFO_FILE_PATH = args.reference_fai_path

# This stage rebuilds the outlier set from scratch.  Keeping files from an
# earlier run is unsafe because event indices are reused and later stages infer
# matrix columns from the directory contents.
if os.path.isdir(TYPE_4_VECTOR_PATH):
    shutil.rmtree(TYPE_4_VECTOR_PATH)
os.makedirs(f"{TYPE_4_VECTOR_PATH}/front_jump", exist_ok=True)
os.makedirs(f"{TYPE_4_VECTOR_PATH}/back_jump", exist_ok=True)

chr_data = find_chr_len(CHROMOSOME_INFO_FILE_PATH)
contig_data = import_data(PREPROCESSED_PAF_FILE_PATH)
contig_data_size = len(contig_data)

with open(f'{args.prefix}/conjoined_type4_ins_del.pkl', 'rb') as f:
    type4_ins, type4_del = pkl.load(file=f)

vcf_type4_events_path = f"{args.prefix}/{VCF_TYPE4_EVENTS_PKL}"
if os.path.isfile(vcf_type4_events_path):
    with open(vcf_type4_events_path, "rb") as f:
        vcf_type4_events = pkl.load(f)
else:
    vcf_type4_events = []
vcf_bnd_type4_events = collect_vcf_bnd_type4_events()

unique_vcf_type4_events = []
seen_vcf_type4_spans = set()
vcf_span_to_event = {}
for event in vcf_type4_events + vcf_bnd_type4_events:
    span_key = (
        event.get("event_type"),
        event.get("chrom"),
        int(event.get("st", 0)),
        int(event.get("nd", 0)),
    )
    if span_key in seen_vcf_type4_spans:
        existing = vcf_span_to_event[span_key]
        for key in ("vcf_id", "mate_id"):
            if event.get(key):
                existing.setdefault("merged_vcf_ids", []).append(event[key])
        continue
    seen_vcf_type4_spans.add(span_key)
    event = dict(event)
    merged_vcf_ids = []
    for key in ("vcf_id", "mate_id"):
        if event.get(key):
            merged_vcf_ids.append(event[key])
    if merged_vcf_ids:
        event["merged_vcf_ids"] = merged_vcf_ids
    vcf_span_to_event[span_key] = event
    unique_vcf_type4_events.append(event)

emitted_indel_candidates = []
merged_indel_candidates = []
vcf_type4_outlier_index = {}

def should_emit_indel_candidate(candidate):
    for prev in emitted_indel_candidates:
        if same_indel_candidate(candidate, prev):
            merged_indel_candidates.append((candidate, prev))
            return False
    emitted_indel_candidates.append(candidate)
    return True

s = 0
cntfj = 0
cntbj = 0
while s<contig_data_size:
    e = contig_data[s][CTG_ENDND]

    # These flanks make the same VCF event traversable by step 02; the
    # canonical event below still owns the single step-11 outlier entry.
    is_vcf_graph_only_type4 = str(contig_data[s][CTG_NAM]).startswith(
        VCF_TYPE4_GRAPH_NODE_PREFIX
    )
    if contig_data[s][CTG_TYP] == 4 and not is_vcf_graph_only_type4:
        chr_name = contig_data[s][CHR_NAM]
        chr_len = chr_data[chr_name]
        rat, ref_st_ed = calculate_single_contig_ref_ratio(contig_data[s:e+1])
        if abs(ref_st_ed[1]-ref_st_ed[0]) >= CHUKJI_LIMIT:
            if rat > 0:
                candidate = make_indel_candidate('front_jump', chr_name, ref_st_ed[0], ref_st_ed[1], f'type4:{s}-{e}')
                if not should_emit_indel_candidate(candidate):
                    s = e+1
                    continue
                cntfj += 1
                write_rows_as_paf(
                    f"{TYPE_4_VECTOR_PATH}/front_jump/{cntfj}.paf",
                    [node_original_or_synthetic_paf_row(contig_data[i]) for i in range(s, e+1)]
                )
                with open(f"{TYPE_4_VECTOR_PATH}/front_jump/{cntfj}_base.paf", "wt") as f:
                    N = ref_st_ed[1] - ref_st_ed[0]
                    virtual_contig = ["base_contig_1", N, 0, N]
                    virtual_contig += ['+', chr_name, chr_len, ref_st_ed[0], ref_st_ed[1]]
                    virtual_contig += [N, N, 0]
                    virtual_contig += ['tp:A:P', 'cs:Z:'+f":{N}"]
                    cigar_str = 'cg:Z:' + cs_to_cigar(virtual_contig[-1][5:])
                    virtual_contig += [cigar_str]
                    for j in virtual_contig:
                        print(j, end="\t", file=f)
                    print("", file=f)
            else:
                candidate = make_indel_candidate('back_jump', chr_name, ref_st_ed[0], ref_st_ed[1], f'type4:{s}-{e}')
                if not should_emit_indel_candidate(candidate):
                    s = e+1
                    continue
                cntbj += 1
                write_rows_as_paf(
                    f"{TYPE_4_VECTOR_PATH}/back_jump/{cntbj}.paf",
                    [node_original_or_synthetic_paf_row(contig_data[i]) for i in range(s, e+1)]
                )
                with open(f"{TYPE_4_VECTOR_PATH}/back_jump/{cntbj}_base.paf", "wt") as f:
                    N = ref_st_ed[0] - ref_st_ed[1]
                    virtual_contig = ["base_contig_1", N, 0, N]
                    virtual_contig += ['+', chr_name, chr_len, ref_st_ed[1], ref_st_ed[0]]
                    virtual_contig += [N, N, 0]
                    virtual_contig += ['tp:A:P', 'cs:Z:'+f":{N}"]
                    cigar_str = 'cg:Z:' + cs_to_cigar(virtual_contig[-1][5:])
                    virtual_contig += [cigar_str]
                    for j in virtual_contig:
                        print(j, end="\t", file=f)
                    print("", file=f)
            
    s = e+1


for event in unique_vcf_type4_events:
    event = dict(event)
    event_type = event.get("event_type")
    chrom = event.get("chrom")
    ref_st = int(event.get("st", 0))
    ref_nd = int(event.get("nd", 0))
    if event_type not in {"front_jump", "back_jump"} or chrom not in chr_data or ref_st == ref_nd:
        logging.warning(f"Skipping malformed VCF Indel event in 11: {event}")
        continue

    svtype = str(event.get("svtype", "")).upper()
    event_size = (
        abs(int(event.get("svlen", 0)))
        if svtype == "INS"
        else abs(ref_nd - ref_st)
    )
    if event_size < VCF_TYPE4_MIN_SPAN:
        logging.warning(
            "Skipping VCF Indel event below the step-11 threshold "
            f"({event_size} < {VCF_TYPE4_MIN_SPAN}): {event}"
        )
        continue

    base_row = make_base_paf_row(chrom, ref_st, ref_nd)
    if svtype == "INS":
        base_rows = event.get("paf_rows") or []
        if not base_rows:
            raise ValueError(
                "VCF INS type4 event is missing query PAF rows: "
                f"{event.get('vcf_id', event)}"
            )
    else:
        base_rows = [base_row]

    if event_type == "front_jump":
        cntfj += 1
        outlier_idx = cntfj
        write_empty_paf(f"{TYPE_4_VECTOR_PATH}/front_jump/{outlier_idx}.paf")
        write_rows_as_paf(f"{TYPE_4_VECTOR_PATH}/front_jump/{outlier_idx}_base.paf", base_rows)
    else:
        cntbj += 1
        outlier_idx = cntbj
        write_empty_paf(f"{TYPE_4_VECTOR_PATH}/back_jump/{outlier_idx}.paf")
        write_rows_as_paf(f"{TYPE_4_VECTOR_PATH}/back_jump/{outlier_idx}_base.paf", base_rows)

    event["outlier_index"] = outlier_idx
    vcf_type4_outlier_index[(event_type, outlier_idx)] = event


type2_indel_cnt = 0

for s1, e1, s2, e2 in type4_ins:
    type2_indel_cnt += 1
    chr_name = contig_data[s1][CHR_NAM]
    chr_len = chr_data[chr_name]

    ref_st = min(contig_data[s1][CHR_STR], contig_data[s1][CHR_END])
    ref_nd = max(contig_data[e2][CHR_STR], contig_data[e2][CHR_END])

    candidate = make_indel_candidate('back_jump', chr_name, ref_st, ref_nd, f'type2_merge:{type2_indel_cnt}')
    if not should_emit_indel_candidate(candidate):
        continue

    cntbj+=1
    write_rows_as_paf(
        f"{TYPE_4_VECTOR_PATH}/back_jump/{cntbj}_type2_merge_{type2_indel_cnt}.paf",
        [node_original_or_synthetic_paf_row(contig_data[i]) for i in (s1, e2)]
    )
        
    with open(f"{TYPE_4_VECTOR_PATH}/back_jump/{cntbj}_base.paf", "wt") as f:
            N = ref_st - ref_nd
            virtual_contig = ["base_contig_1", N, 0, N]
            virtual_contig += ['+', chr_name, chr_len, ref_nd, ref_st]
            virtual_contig += [N, N, 0]
            virtual_contig += ['tp:A:P', 'cs:Z:'+f":{N}"]
            cigar_str = 'cg:Z:' + cs_to_cigar(virtual_contig[-1][5:])
            virtual_contig += [cigar_str]
            for j in virtual_contig:
                print(j, end="\t", file=f)
            print("", file=f)

for s1, e1, s2, e2 in type4_del:
    type2_indel_cnt += 1
    chr_name = contig_data[s1][CHR_NAM]
    chr_len = chr_data[chr_name]

    ref_st = min(contig_data[s1][CHR_STR], contig_data[s1][CHR_END])
    ref_nd = max(contig_data[e2][CHR_STR], contig_data[e2][CHR_END])

    candidate = make_indel_candidate('front_jump', chr_name, ref_st, ref_nd, f'type2_merge:{type2_indel_cnt}')
    if not should_emit_indel_candidate(candidate):
        continue

    cntfj+=1
    write_rows_as_paf(
        f"{TYPE_4_VECTOR_PATH}/front_jump/{cntfj}_type2_merge_{type2_indel_cnt}.paf",
        [node_original_or_synthetic_paf_row(contig_data[i]) for i in (s1, e2)]
    )
        
    with open(f"{TYPE_4_VECTOR_PATH}/front_jump/{cntfj}_base.paf", "wt") as f:
            N = ref_nd - ref_st
            virtual_contig = ["base_contig_1", N, 0, N]
            virtual_contig += ['+', chr_name, chr_len, ref_st, ref_nd]
            virtual_contig += [N, N, 0]
            virtual_contig += ['tp:A:P', 'cs:Z:'+f":{N}"]
            cigar_str = 'cg:Z:' + cs_to_cigar(virtual_contig[-1][5:])
            virtual_contig += [cigar_str]
            for j in virtual_contig:
                print(j, end="\t", file=f)
            print("", file=f)

with open(f"{args.prefix}/{VCF_TYPE4_OUTLIER_INDEX_PKL}", "wb") as f:
    pkl.dump(vcf_type4_outlier_index, f)

logging.info(f"Forward-directed outlier contig count : {cntfj}")
logging.info(f"Backward-directed outlier contig count : {cntbj}")
logging.info(f"VCF-derived Indel outlier count : {len(vcf_type4_outlier_index)}")
logging.info(f"VCF BND-derived Indel outlier count : {len(vcf_bnd_type4_events)}")
logging.info(f"Merged duplicate indel candidate count : {len(merged_indel_candidates)}")
for candidate, prev in merged_indel_candidates[:20]:
    logging.debug(
        "Merged duplicate indel candidate "
        f"{candidate['source']} {candidate['event_type']}:{candidate['chrom']}:{candidate['st']}-{candidate['nd']} "
        "into "
        f"{prev['source']} {prev['event_type']}:{prev['chrom']}:{prev['st']}-{prev['nd']} "
        f"(tolerance={INDEL_MERGE_TOLERANCE})"
    )
logging.info(f"Total count : {cntfj + cntbj}")
