import os
import argparse
import pandas as pd

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
K = 1000
INF = 10000000000

CONTIG_MINIMUM_SIZE = 100*K

BND_CONTIG_BOUND = 0.1

def dbg() :
    print("hi")

def import_data(file_path : list) -> list :
    contig_data = []
    int_induce_idx = [1, 2, 3, 6, 7, 8]
    with open(file_path, "r") as paf_file:
        for curr_contig in paf_file:
            temp_list = curr_contig.split("\t")[:9]
            for i in int_induce_idx:
                temp_list[i] = int(temp_list[i])
            contig_data.append(temp_list)
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

def calc_ratio(contig_data : list) -> dict:
    contig_data_size = len(contig_data)
    contig_data.append((0, 0, 0, 0, 0, 0, 0, 0, 0))
    total_ref_len = 0
    curr_contig_first_fragment = contig_data[0]
    ref_qry_ratio = {}
    for i in range(1, contig_data_size+1):
        if contig_data[i][CTG_NAM] != contig_data[i-1][CTG_NAM]:
            curr_contig_name = contig_data[i-1][CTG_NAM]
            curr_contig_end_fragment = contig_data[i-1]
            total_ref_len = contig_data[i-1][CTG_LEN]
            if curr_contig_first_fragment[CTG_DIR] == '+':
                estimated_ref_len = curr_contig_end_fragment[CHR_END] - curr_contig_first_fragment[CHR_STR]
            else:
                estimated_ref_len = curr_contig_first_fragment[CHR_END] - curr_contig_end_fragment[CHR_STR]
            try:
                ref_qry_ratio[curr_contig_name] = total_ref_len / estimated_ref_len
            except:
                if total_ref_len >=0:
                    ref_qry_ratio[curr_contig_name] = INF
                else:
                    ref_qry_ratio[curr_contig_name] = -INF
            total_ref_len = 0
            curr_contig_first_fragment = contig_data[i]
    contig_data = contig_data[:-1]
    return ref_qry_ratio        
    

def preprocess_contig(contig_data : list, telo_label : list, ref_qry_ratio : dict) -> list :
    checker = 0
    contig_data_size = len(contig_data)
    contig_data.append((0, 0, 0, 0, 0, 0, 0, 0, 0))
    curr_contig_first_fragment = contig_data[0]
    using_contig_list = []
    contig_terminal_node = {}
    contig_type = {}
    is_telo = False
    idx = 0
    cnt = 0
    f = open("test.txt", "wt")
    for i in range(1, contig_data_size+1):
        cnt+=1
        if telo_label[i-1]:
            is_telo = True
        # contig 넘어갈 때:
        if contig_data[i][CTG_NAM] != contig_data[i-1][CTG_NAM]:
            curr_contig_name = contig_data[i-1][CTG_NAM]
            curr_contig_end_fragment = contig_data[i-1]
            if curr_contig_first_fragment[CHR_NAM] != curr_contig_end_fragment[CHR_NAM]:
                checker = 1
            elif curr_contig_first_fragment[CTG_DIR] != curr_contig_end_fragment[CTG_DIR]:
                checker = 2
            else:
                print(curr_contig_name, ref_qry_ratio[curr_contig_name], file = f)
                if abs(ref_qry_ratio[curr_contig_name]-1)< BND_CONTIG_BOUND:
                    if contig_data[i-1][CTG_LEN] > CONTIG_MINIMUM_SIZE or is_telo:
                        checker = 3
                    elif cnt==1:
                        checker = 3
                else:
                    checker = 4
            if checker>0:
                using_contig_list.append(curr_contig_name)
                contig_type[curr_contig_name] = checker
                contig_terminal_node[curr_contig_name] = (idx, idx+cnt-1)
                idx+=cnt
            # initialize
            curr_contig_first_fragment = contig_data[i]
            cnt = 0
            checker = 0
            is_telo = False
    contig_data = contig_data[:-1]
    f.close()
    return [using_contig_list, contig_type, contig_terminal_node]

def alt_preprocess_contig(contig_data : list, telo_label : list, ref_qry_ratio : dict) -> list :
    checker = 0
    contig_data = contig_data[:-1]
    contig_data_size = len(contig_data)
    contig_data.append((0, 0, 0, 0, 0, 0, 0, 0, 0))
    curr_contig_first_fragment = contig_data[0]
    using_contig_list = []
    contig_terminal_node = {}
    contig_type = {}
    is_telo = False
    idx = 0
    cnt = 0
    for i in range(1, contig_data_size+1):
        cnt+=1
        try:
            if telo_label[i-1]:
                is_telo = True
        except:
            print(contig_data_size, i-1)
            print(len(telo_label))
            exit(1)
        # contig 넘어갈 때:
        if contig_data[i][CTG_NAM] != contig_data[i-1][CTG_NAM]:
            curr_contig_name = contig_data[i-1][CTG_NAM]
            curr_contig_end_fragment = contig_data[i-1]
            if curr_contig_first_fragment[CHR_NAM] != curr_contig_end_fragment[CHR_NAM]:
                checker = 1
            elif curr_contig_first_fragment[CTG_DIR] != curr_contig_end_fragment[CTG_DIR]:
                checker = 2
            else:
                if is_telo:
                    checker = 3
                elif abs(ref_qry_ratio[curr_contig_name]-1) >= BND_CONTIG_BOUND:
                    checker = 4
            if checker>0:
                using_contig_list.append(curr_contig_name)
                contig_type[curr_contig_name] = checker
                contig_terminal_node[curr_contig_name] = (idx, idx+cnt-1)
                idx+=cnt
            # initialize
            curr_contig_first_fragment = contig_data[i]
            cnt = 0
            checker = 0
    contig_data = contig_data[:-1]
    return [using_contig_list, contig_type, contig_terminal_node]

def find_mid_breakend(contig_data : list, using_contig_nam : list, ref_qry_ratio : dict, contig_type : dict) -> list:
    contig_info = dict()
    contig_nam = []
    curr_contig_first_fragment = contig_data[0]
    contig_data_size = len(contig_data)
    contig_data.append((0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    for i in range(1, contig_data_size+1):
        if contig_data[i][CTG_NAM] != contig_data[i-1][CTG_NAM]:
            contig_nam.append(contig_data[i][CTG_NAM])
            curr_contig_end_fragment = contig_data[i-1]
            if contig_type[contig_data[i-1][CTG_NAM]] == 4:
                if contig_data[i-1][CTG_DIR]=='+':
                    contig_info[contig_data[i-1][CTG_NAM]] = (contig_data[i-1][CHR_NAM], curr_contig_first_fragment[CTG_DIR], curr_contig_first_fragment[CHR_STR], curr_contig_end_fragment[CHR_END])
                else:
                    contig_info[contig_data[i-1][CTG_NAM]] = (contig_data[i-1][CHR_NAM], curr_contig_first_fragment[CTG_DIR], curr_contig_end_fragment[CHR_STR], curr_contig_first_fragment[CHR_END])
            curr_contig_first_fragment = contig_data[i]
    contig_data = contig_data[0:-1]
    delete_contig = []
    # 뒤로 돌아가는 contig 찾기
    for i in contig_nam:
        if ref_qry_ratio[i]<0 and contig_type[i]==4:
            # 이 중에서 안에 들어가는 걸 찾아야 함
            for j in contig_info:
                if contig_info[i][0] == contig_info[j][0] and i!=j:
                    if contig_info[i][2] < contig_info[j][1] and contig_info[j][2] < contig_info[i][1]:
                        delete_contig.append(i)
    return delete_contig

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
    if args.alt is None:
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

    ref_qry_ratio = calc_ratio(contig_data)

    for _ in telo_data:
        telo_dict[_[0]].append(_[1:])

    node_label = label_node(contig_data, telo_dict)
    preprocess_result, preprocess_contig_type, preprocess_terminal_nodes = preprocess_contig(contig_data, node_label, ref_qry_ratio)
    contig_data_size = len(contig_data)
    if args.alt is not None:
        alt_contig_data = import_data(PAF_FILE_PATH[1])
        alt_node_label = label_node(alt_contig_data, telo_dict)
        alt_ref_qry_ratio = calc_ratio(alt_contig_data)
        alt_preprocess_result, alt_preprocess_contig_type, alt_preprocess_terminal_nodes = alt_preprocess_contig(alt_contig_data, alt_node_label, alt_ref_qry_ratio)
        alt_contig_data_size = len(alt_contig_data)
        test = []
        for _ in preprocess_contig_type:
            test.append(preprocess_contig_type[_])
        for _ in alt_preprocess_contig_type:
            test.append(alt_preprocess_contig_type[_])
        test = pd.Series(test)
        print(test.value_counts())
        
        print(f"Before Preprocessing: {contig_data_size+alt_contig_data_size} nodes exist")
    else:
        print(f"Before Preprocessing: {contig_data_size} nodes exist")

    cnt = 0
    idx = 0
    with open(PREPROCESSED_PAF_FILE_PATH, "wt") as processed_paf:
        for i in range(contig_data_size):
            if contig_data[i][CTG_NAM] in preprocess_result:
                for __ in contig_data[i]: 
                    processed_paf.write(str(__)+"\t")
                processed_paf.write(str(preprocess_contig_type[contig_data[i][CTG_NAM]])+"\t")
                processed_paf.write(str(preprocess_terminal_nodes[contig_data[i][CTG_NAM]][0])+"\t")
                processed_paf.write(str(preprocess_terminal_nodes[contig_data[i][CTG_NAM]][1])+"\t")
                processed_paf.write(node_label[i][0]+"\t"+node_label[i][1]+"\t")
                processed_paf.write("0." + str(idx) + "\n")
                cnt+=1
            idx+=1
        idx = 0
        bias = cnt
        if args.alt is not None:
            for i in range(alt_contig_data_size):
                if alt_contig_data[i][CTG_NAM] in alt_preprocess_result:
                    for __ in alt_contig_data[i]: 
                        processed_paf.write(str(__)+"\t")
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