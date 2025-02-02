import os
import argparse

from collections import defaultdict

CTG_NAM = 0
CTG_LEN = 1
CTG_STR = 2
CTG_END = 3
CTG_DIR = 4
CHR_NAM = 5
CHR_LEN = 6
CHR_STR = 7
CHR_END = 8
CTG_NUM = 9
K = 1

CONTIG_MINIMUM_SIZE = 500 * K

def dbg() :
    print("hi")

def import_data(file_path : list) -> list :
    contig_data = []
    with open(file_path, "r") as paf_file:
        for curr_contig in paf_file:
            contig_data.append(curr_contig.split("\t")[:9])
    return contig_data

def find_chr_len(file_path : str) -> dict:
    chr_data_file = open(file_path, "r")
    chr_len = {}
    for curr_data in chr_data_file:
        curr_data = curr_data.split("\t")
        chr_len[curr_data[0]] = int(curr_data[1])
    return chr_len

def import_telo_data(file_path : str, chr_len : dict) -> dict :
    fai_file = open(file_path, "r")
    telo_data = [(0,1)]
    for curr_data in fai_file:
        temp_list = curr_data.split("\t")
        int_induce_idx = [1, 2]
        for i in int_induce_idx:
            temp_list[i] = int(temp_list[i])
        if temp_list[0]!=telo_data[-1][0]:
            temp_list.append('f')
        else:
            if temp_list[1]>chr_len[temp_list[0]]/2:
                temp_list.append('b')
            else:
                temp_list.append('f')
        telo_data.append(tuple(temp_list))
    return telo_data[1:]


def distance_checker(node_a : tuple, node_b : tuple) -> int :
    # check if overlap
    if max(int(node_a[CHR_STR]), int(node_b[CHR_STR])) < min(int(node_a[CHR_END]), int(node_b[CHR_END])):
        return 0   
    else:
        # return distance of two ranges
        return min(abs(int(node_b[CHR_STR]) - int(node_a[CHR_END])), abs(int(node_b[CHR_END]) - int(node_a[CHR_STR])))

def inclusive_checker(contig_node : tuple, telomere_node : tuple) -> bool :
    if int(telomere_node[CHR_STR]) <= int(contig_node[CHR_STR]) and int(contig_node[CHR_END]) <= int(telomere_node[CHR_END]):
        return True
    else:
        return False

def label_node(contig_data : list, telo_data) -> list :
    label = []
    contig_data_size = len(contig_data)
    for i in range(contig_data_size):
        checker = 0
        for j in telo_data[contig_data[i][CHR_NAM]]:
            if distance_checker(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0], j[1])) == 0:
                inclusive_label = ""
                if inclusive_checker(contig_data[i], (0, 0, 0, 0, 0, 0, 0, j[0], j[1])):
                    inclusive_label = "in"
                label.append((contig_data[i][CHR_NAM], j[2]+inclusive_label))
                checker = 1
                break
        if checker==0:
            label.append(('0', '0'))
    return label

def preprocess_contig(contig_data : list, telo_label : list) -> list :
    curr_contig_name = contig_data[0][CTG_NAM]
    curr_contig_first_fragment = contig_data[0]
    curr_contig_first_fragment_index = 0
    prev_fragment = contig_data[0]
    using_contig_list = []
    contig_terminal_node = {}
    contig_type = {}
    cnt = 0
    idx = 0
    k = 0
    telo_exist = 0
    for i in contig_data:
        if i[CTG_NAM] != curr_contig_name: # 이제 다음 contig로 넘어갈 예정
            checker = 0
            # 길이가 짧을 경우,
            # 1) 텔로미어가 존재하는지
            # 2) 양 말단 방향이 다른지
            # 3) 양 말단 CHR 종류가 다른지 확인.
            if int(prev_fragment[CTG_LEN]) < CONTIG_MINIMUM_SIZE and telo_exist == 0 \
            and prev_fragment[CHR_NAM] == curr_contig_first_fragment[CHR_NAM] \
            and prev_fragment[CTG_DIR] == curr_contig_first_fragment[CTG_DIR]:
                checker = -1
            # 1. Check if contig's 1st & last chromosome is same
            elif prev_fragment[CHR_NAM] != curr_contig_first_fragment[CHR_NAM]:
                checker = 1
            # 2. Check if contig's 1st & last fragment have same direction
            elif prev_fragment[CTG_DIR] != curr_contig_first_fragment[CTG_DIR]:
                checker = 2
            # 3. Check if contig's direction has consistency
            else:
                if prev_fragment[CTG_DIR] == '-':
                    if int(curr_contig_first_fragment[CHR_STR]) >= int(prev_fragment[CHR_END]):
                        checker = 3
                if prev_fragment[CTG_DIR] == '+':
                    if int(curr_contig_first_fragment[CHR_END]) <= int(prev_fragment[CHR_STR]):
                        checker = 3
            if checker>0:
                using_contig_list.append(curr_contig_name)
                contig_type[curr_contig_name] = checker
                contig_terminal_node[curr_contig_name] = (idx, idx+cnt-1)
                idx+=cnt
            if cnt==1 and checker>=0:
                using_contig_list.append(curr_contig_name)
                contig_type[curr_contig_name] = 3
                contig_terminal_node[curr_contig_name] = (idx, idx+cnt-1)
                idx+=cnt
            curr_contig_name = i[0]
            curr_contig_first_fragment = i
            cnt=0
            telo_exist = 0
        if telo_label[k][0] != '0':
            telo_exist = 1
        prev_fragment = i
        cnt+=1
        k+=1
    checker = 0
    # Last element check
    if int(prev_fragment[CTG_LEN]) < CONTIG_MINIMUM_SIZE and telo_exist == 0 \
    and prev_fragment[CHR_NAM] == curr_contig_first_fragment[CHR_NAM] \
    and prev_fragment[CTG_DIR] == curr_contig_first_fragment[CTG_DIR]:
        checker = -1
    elif prev_fragment[CHR_NAM] != curr_contig_first_fragment[CHR_NAM]:
        checker = 1
    elif prev_fragment[CTG_DIR] != curr_contig_first_fragment[CTG_DIR]:
        checker = 2
    else:
        if prev_fragment[CTG_DIR] == '-':
            if int(curr_contig_first_fragment[CHR_STR]) >= int(prev_fragment[CHR_END]):
                checker = 3
        if prev_fragment[CTG_DIR] == '+':
            if int(curr_contig_first_fragment[CHR_END]) <= int(prev_fragment[CHR_STR]):
                checker = 3
    if checker>0:
        using_contig_list.append(curr_contig_name)
        contig_type[curr_contig_name] = checker
        contig_terminal_node[curr_contig_name] = (idx, idx+cnt-1)
        idx+=cnt
    if cnt==1 and checker>=0:
        print("hi")
        using_contig_list.append(curr_contig_name)
        contig_type[curr_contig_name] = 3
        contig_terminal_node[curr_contig_name] = (idx, idx+cnt-1)
        idx+=cnt
    return [using_contig_list, contig_type, contig_terminal_node]

def alt_preprocess_contig(contig_data : list, telo_label : list) -> list :
    curr_contig_name = contig_data[0][CTG_NAM]
    curr_contig_first_fragment = contig_data[0]
    curr_contig_first_fragment_index = 0
    prev_fragment = contig_data[0]
    using_contig_list = []
    contig_terminal_node = {}
    contig_type = {}
    cnt = 0
    idx = 0
    k = 0
    telo_exist = 0
    for i in contig_data:
        if i[CTG_NAM] != curr_contig_name: # 이제 다음 contig로 넘어갈 예정
            checker = 0
            # 길이가 짧을 경우,
            # 1) 텔로미어가 존재하는지
            # 2) 양 말단 방향이 다른지
            # 3) 양 말단 CHR 종류가 다른지 확인.
            if int(prev_fragment[CTG_LEN]) < CONTIG_MINIMUM_SIZE:
                checker = 0
            # 1. Check if contig's 1st & last chromosome is same
            elif prev_fragment[CHR_NAM] != curr_contig_first_fragment[CHR_NAM]:
                checker = 1
            # 2. Check if contig's 1st & last fragment have same direction
            elif prev_fragment[CTG_DIR] != curr_contig_first_fragment[CTG_DIR]:
                checker = 2
            # 3. Check if contig's direction has consistency
            else:
                if prev_fragment[CTG_DIR] == '-':
                    if int(curr_contig_first_fragment[CHR_STR]) >= int(prev_fragment[CHR_END]):
                        checker = 3
                if prev_fragment[CTG_DIR] == '+':
                    if int(curr_contig_first_fragment[CHR_END]) <= int(prev_fragment[CHR_STR]):
                        checker = 3
            if checker>0:
                if telo_exist \
                or curr_contig_first_fragment[CHR_NAM] != prev_fragment[CHR_NAM] \
                or curr_contig_first_fragment[CTG_DIR] != prev_fragment[CTG_DIR]:
                    using_contig_list.append(curr_contig_name)
                    contig_type[curr_contig_name] = checker
                    contig_terminal_node[curr_contig_name] = (idx, idx+cnt-1)
                    idx+=cnt
            curr_contig_name = i[0]
            curr_contig_first_fragment = i
            cnt=0
            telo_exist = 0
        if telo_label[k][0] != '0':
            telo_exist = 1
        prev_fragment = i
        cnt+=1
        k+=1
    checker = 0
    # Last element check
    if int(prev_fragment[CTG_LEN]) < CONTIG_MINIMUM_SIZE:
        checker = 0
    # 1. Check if contig's 1st & last chromosome is same
    elif prev_fragment[CHR_NAM] != curr_contig_first_fragment[CHR_NAM]:
        checker = 1
    # 2. Check if contig's 1st & last fragment have same direction
    elif prev_fragment[CTG_DIR] != curr_contig_first_fragment[CTG_DIR]:
        checker = 2
    # 3. Check if contig's direction has consistency
    else:
        if prev_fragment[CTG_DIR] == '-':
            if int(curr_contig_first_fragment[CHR_STR]) >= int(prev_fragment[CHR_END]):
                checker = 3
        if prev_fragment[CTG_DIR] == '+':
            if int(curr_contig_first_fragment[CHR_END]) <= int(prev_fragment[CHR_STR]):
                checker = 3
    if checker>0:
        if telo_exist \
        or curr_contig_first_fragment[CHR_NAM] != prev_fragment[CHR_NAM] \
        or curr_contig_first_fragment[CTG_DIR] != prev_fragment[CTG_DIR]:
            using_contig_list.append(curr_contig_name)
            contig_type[curr_contig_name] = checker
            contig_terminal_node[curr_contig_name] = (idx, idx+cnt-1)
            idx+=cnt
    return [using_contig_list, contig_type, contig_terminal_node]

def main():
    parser = argparse.ArgumentParser(description="Process file paths for telomere analysis.")

    # 위치 인자 정의
    parser.add_argument("paf_file_path", 
                        help="Path to the original PAF file.")
    parser.add_argument("telomere_bed_path", 
                        help="Path to the telomere information file.")
    parser.add_argument("reference_fai_path", 
                        help="Path to the chromosome information file.")
    
    parser.add_argument("--alt", 
                        help="Path to an alternative PAF file (optional).")

    # 인자 파싱
    args = parser.parse_args()

    PAF_FILE_PATH = []
    if args.alt is None or args.alt == "a":
        PAF_FILE_PATH = [args.paf_file_path]
    else:
        PAF_FILE_PATH = [args.paf_file_path, args.alt]
    
    TELOMERE_INFO_FILE_PATH = args.telomere_bed_path
    CHROMOSOME_INFO_FILE_PATH = args.reference_fai_path
    PREPROCESSED_PAF_FILE_PATH = (PAF_FILE_PATH[0] +'.ppc.paf')

    contig_data = import_data(PAF_FILE_PATH[0])

    chr_len = find_chr_len(CHROMOSOME_INFO_FILE_PATH)

    telo_data = import_telo_data(TELOMERE_INFO_FILE_PATH, chr_len)

    telo_dict = defaultdict(list)

    for _ in telo_data:
        telo_dict[_[0]].append(_[1:])

    node_label = label_node(contig_data, telo_dict)
    preprocess_result, preprocess_contig_type, preprocess_terminal_nodes = preprocess_contig(contig_data, node_label)
    contig_data_size = len(contig_data)
    if args.alt is not None and args.alt != "a":
        alt_contig_data = import_data(PAF_FILE_PATH[1])
        alt_node_label = label_node(alt_contig_data, telo_dict)
        alt_preprocess_result, alt_preprocess_contig_type, alt_preprocess_terminal_nodes = alt_preprocess_contig(alt_contig_data, alt_node_label)
        alt_contig_data_size = len(alt_contig_data)
        print(f"Before Preprocessing: {contig_data_size+alt_contig_data_size} nodes exist")
    else:
        print(f"Before Preprocessing: {contig_data_size} nodes exist")

    cnt = 0
    idx = 0
    with open(PREPROCESSED_PAF_FILE_PATH, "wt") as processed_paf:
        for i in range(contig_data_size):
            if contig_data[i][CTG_NAM] in preprocess_result:
                for __ in contig_data[i]: 
                    processed_paf.write(__+"\t")
                processed_paf.write(str(preprocess_contig_type[contig_data[i][CTG_NAM]])+"\t")
                processed_paf.write(str(preprocess_terminal_nodes[contig_data[i][CTG_NAM]][0])+"\t")
                processed_paf.write(str(preprocess_terminal_nodes[contig_data[i][CTG_NAM]][1])+"\t")
                processed_paf.write(node_label[i][0]+"\t"+node_label[i][1]+"\t")
                processed_paf.write("0." + str(idx) + "\n")
                cnt+=1
            idx+=1
        idx = 0
        bias = cnt
        if args.alt is not None and args.alt != "a":
            for i in range(alt_contig_data_size):
                if alt_contig_data[i][CTG_NAM] in alt_preprocess_result:
                    for __ in alt_contig_data[i]: 
                        processed_paf.write(__+"\t")
                    processed_paf.write(str(alt_preprocess_contig_type[alt_contig_data[i][CTG_NAM]])+"\t")
                    processed_paf.write(str(bias+alt_preprocess_terminal_nodes[alt_contig_data[i][CTG_NAM]][0])+"\t")
                    processed_paf.write(str(bias+alt_preprocess_terminal_nodes[alt_contig_data[i][CTG_NAM]][1])+"\t")
                    processed_paf.write(alt_node_label[i][0]+"\t"+alt_node_label[i][1]+"\t")
                    processed_paf.write("1." + str(idx) + "\n")
                    cnt+=1
                idx+=1

                

    print(f"After Preprocessing: {len(preprocess_result)} contigs, {cnt} nodes exist")

if __name__ == "__main__":
    main()