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
CTG_CNT = 9
CTG_CONTEL = 10

K = 1000
INF = 10000000000
CONTIG_MINIMUM_SIZE = 100*K
BND_CONTIG_BOUND = 0.1
CHROMOSOME_COUNT = 23


def dbg() :
    print("hi")

def import_data(file_path : list) -> list :
    contig_data = []
    int_induce_idx = [1, 2, 3, 6, 7, 8]
    idx = 0
    with open(file_path, "r") as paf_file:
        for curr_contig in paf_file:
            temp_list = curr_contig.split("\t")[:9]
            for i in int_induce_idx:
                temp_list[i] = int(temp_list[i])
            temp_list.append(idx)
            contig_data.append(temp_list)
            idx+=1
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
    if max(int(node_a[CHR_STR]), int(node_b[CHR_STR])) < min(int(node_a[CHR_END]), int(node_b[CHR_END])):
        return 0   
    else:
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

def preprocess_telo(contig_data : list, node_label : list) -> list :
    telo_preprocessed_contig = []
    telo_connect_info = {}
    contig_endpoint = {}
    report_case = {'A':[], 'B':[], 'C':[], 'ESC':[], 'ALL_TELO_NON_ESC':[]}
    chr_corr = {}
    chr_rev_corr = {}
    total_contig_count = len(contig_data)
    for i in range(1, CHROMOSOME_COUNT):
        chr_corr['chr'+str(i)+'f'] = total_contig_count + i - 1
        chr_rev_corr[total_contig_count + i - 1] = 'chr'+str(i)+'f'
    chr_corr['chrXf'] = total_contig_count + CHROMOSOME_COUNT - 1
    chr_corr['chrYf'] = total_contig_count + CHROMOSOME_COUNT - 1
    chr_rev_corr[total_contig_count + CHROMOSOME_COUNT - 1] = 'chrXf'
    for i in range(1, CHROMOSOME_COUNT):
        chr_corr['chr'+str(i)+'b'] = total_contig_count + CHROMOSOME_COUNT + i - 1
        chr_rev_corr[total_contig_count + CHROMOSOME_COUNT + i - 1] = 'chr'+str(i)+'b'
    chr_corr['chrXb'] = total_contig_count + 2*CHROMOSOME_COUNT - 1
    chr_corr['chrYb'] = total_contig_count + 2*CHROMOSOME_COUNT - 1
    chr_rev_corr[total_contig_count + 2*CHROMOSOME_COUNT - 1] = 'chrXb'
    contig_data_size = len(contig_data)
    curr_contig_st = 0
    contig_data.append((0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    for i in range(1, contig_data_size+1):
        if contig_data[i-1][CTG_NAM] != contig_data[i][CTG_NAM]:
            curr_contig_ed = i-1
            front_telo_bound = curr_contig_st
            end_telo_bound = curr_contig_ed
            while node_label[front_telo_bound][0] != '0' and front_telo_bound<=curr_contig_ed:
                front_telo_bound+=1
            while node_label[end_telo_bound][0] != '0' and end_telo_bound>=curr_contig_st:
                end_telo_bound-=1
            st = curr_contig_st
            ed = curr_contig_ed
            escape = False
            if front_telo_bound == curr_contig_ed+1:
                if len(node_label[curr_contig_st][1])==1 \
                and (node_label[curr_contig_st][1] + contig_data[curr_contig_st][CTG_DIR] in ("b+", "f-")):
                    escape = True
                elif len(node_label[curr_contig_ed][1])==1 \
                and (node_label[curr_contig_ed][1] + contig_data[curr_contig_ed][CTG_DIR] in ("b-", "f+")):
                    escape = True
            else:
                escape = True
            if escape:
                if front_telo_bound > curr_contig_st:
                    front_telo_bound-=1
                    # Check if first telomere node is 'Nin' -> Else: escape
                    if len(node_label[curr_contig_st][1])==1 \
                    and (node_label[curr_contig_st][1] + contig_data[curr_contig_st][CTG_DIR] in ("b+", "f-")):
                        report_case['ESC'].append(curr_contig_st)
                        st = curr_contig_st
                    elif len(node_label[curr_contig_st][1])>1:
                        # Next node is connected with telomere node.
                        front_telo_bound+=1
                        if front_telo_bound <= curr_contig_ed:
                            dest = contig_data[front_telo_bound][CHR_NAM]
                            if contig_data[front_telo_bound][CTG_DIR] == '+':
                                dest += 'f'
                                telo_connect_info[front_telo_bound] = dest
                                report_case['A'].append([dest, front_telo_bound])
                            else:
                                dest += 'b'
                                telo_connect_info[front_telo_bound] = dest
                                report_case['A'].append([dest, front_telo_bound])
                        st = front_telo_bound
                    # If boundary node is not "Nin"
                    else:
                        if node_label[curr_contig_st][1] == 'f':
                            if contig_data[front_telo_bound][CTG_DIR]=='+':
                                # boundary node is connected with telomere node.
                                dest = contig_data[front_telo_bound][CHR_NAM]+'f'
                                telo_connect_info[front_telo_bound] = dest
                                report_case['B'].append([dest, front_telo_bound])
                            else:
                                # Treat as "Nin" Case
                                front_telo_bound+=1
                                if front_telo_bound <= curr_contig_ed:
                                    dest = contig_data[front_telo_bound][CHR_NAM]
                                    if contig_data[front_telo_bound][CTG_DIR] == '+':
                                        dest += 'f'
                                        telo_connect_info[front_telo_bound] = dest
                                        report_case['C'].append([dest, front_telo_bound])
                                    else:
                                        dest += 'b'
                                        telo_connect_info[front_telo_bound] = dest
                                        report_case['C'].append([dest, front_telo_bound])
                        else:
                            if contig_data[front_telo_bound][CTG_DIR]=='-':
                                dest = contig_data[front_telo_bound][CHR_NAM]+'b'
                                telo_connect_info[front_telo_bound] = dest
                                report_case['B'].append([dest, front_telo_bound])
                            else:
                                front_telo_bound+=1
                                if front_telo_bound <= curr_contig_ed:
                                    dest = contig_data[front_telo_bound][CHR_NAM]
                                    if contig_data[front_telo_bound][CTG_DIR] == '+':
                                        dest += 'f'
                                        telo_connect_info[front_telo_bound] = dest
                                        report_case['C'].append([dest, front_telo_bound])
                                    else:
                                        dest += 'b'
                                        telo_connect_info[front_telo_bound] = dest
                                        report_case['C'].append([dest, front_telo_bound])
                        st = front_telo_bound
                
                if end_telo_bound < curr_contig_ed:
                    end_telo_bound+=1
                    # Check if first telomere node is 'Nin' -> Else: escape
                    if len(node_label[curr_contig_ed][1])==1 \
                    and (node_label[curr_contig_ed][1] + contig_data[curr_contig_ed][CTG_DIR] in ("b-", "f+")):
                        report_case['ESC'].append(curr_contig_ed)
                        ed = curr_contig_ed
                    # If boundary node is 'Nin'
                    elif len(node_label[curr_contig_ed][1])>1:
                        # Next node is connected with telomere node.
                        end_telo_bound-=1
                        if end_telo_bound>=curr_contig_st:
                            dest = contig_data[end_telo_bound][CHR_NAM]
                            if contig_data[end_telo_bound][CTG_DIR] == '+':
                                dest += 'b'
                                telo_connect_info[end_telo_bound] = dest
                                report_case['A'].append([dest, end_telo_bound])
                            else:
                                dest += 'f'
                                telo_connect_info[end_telo_bound] = dest
                                report_case['A'].append([dest, end_telo_bound])
                        ed = end_telo_bound
                    # If boundary node is not "Nin"
                    else:
                        if node_label[curr_contig_ed][1] == 'b':
                            if contig_data[end_telo_bound][CTG_DIR]=='+':
                                # boundary node is connected with telomere node.
                                dest = contig_data[end_telo_bound][CHR_NAM]+'b'
                                telo_connect_info[end_telo_bound] = dest
                                report_case['B'].append([dest, end_telo_bound])
                            else:
                                # Treat as "Nin" Case
                                end_telo_bound+=1
                                if end_telo_bound>=curr_contig_st:
                                    dest = contig_data[end_telo_bound][CHR_NAM]
                                    if contig_data[end_telo_bound][CTG_DIR] == '+':
                                        dest += 'b'
                                        telo_connect_info[end_telo_bound] = dest
                                        report_case['C'].append([dest, end_telo_bound])
                                    else:
                                        dest += 'f'
                                        telo_connect_info[end_telo_bound] = dest
                                        report_case['C'].append([dest, end_telo_bound])
                        else:
                            if contig_data[end_telo_bound][CTG_DIR]=='-':
                                dest = contig_data[end_telo_bound][CHR_NAM]+'f'
                                telo_connect_info[end_telo_bound] = dest
                                report_case['B'].append([dest, end_telo_bound])
                            else:
                                end_telo_bound+=1
                                if end_telo_bound>=curr_contig_st:
                                    dest = contig_data[end_telo_bound][CHR_NAM]
                                    if contig_data[end_telo_bound][CTG_DIR] == '+':
                                        dest += 'b'
                                        telo_connect_info[end_telo_bound] = dest
                                        report_case['C'].append([dest, end_telo_bound])
                                    else:
                                        dest += 'f'
                                        telo_connect_info[end_telo_bound] = dest
                                        report_case['C'].append([dest, end_telo_bound])
                        ed = end_telo_bound
                for j in range(st, ed+1):
                    telo_preprocessed_contig.append(j)
            else:
                report_case['ALL_TELO_NON_ESC'].append(curr_contig_st)
            # initialize
            curr_contig_st = curr_contig_ed+1
    contig_data = contig_data[0:-1]
    return telo_preprocessed_contig, report_case, telo_connect_info

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

def calc_ratio(contig_data : list) -> dict:
    contig_data_size = len(contig_data)
    contig_data.append((0, 0, 0, 0, 0, 0, 0, 0, 0))
    total_ref_len = 0
    curr_contig_first_fragment = contig_data[0]
    ref_qry_ratio = {}
    for i in range(1, contig_data_size+1):
        total_ref_len += contig_data[i-1][CHR_END] - contig_data[i-1][CHR_STR]
        if contig_data[i][CTG_NAM] != contig_data[i-1][CTG_NAM]:
            curr_contig_name = contig_data[i-1][CTG_NAM]
            curr_contig_end_fragment = contig_data[i-1]
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
    contig_data_size = len(contig_data)-1
    curr_contig_first_fragment = contig_data[0]
    using_contig_list = []
    contig_terminal_node = {}
    contig_type = {}
    is_telo = False
    idx = 0
    cnt = 0
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
    return [using_contig_list, contig_type, contig_terminal_node]

def alt_preprocess_contig(contig_data : list, telo_label : list, ref_qry_ratio : dict) -> list :
    checker = 0
    contig_data_size = len(contig_data)-1
    curr_contig_first_fragment = contig_data[0]
    using_contig_list = []
    contig_terminal_node = {}
    contig_type = {}
    is_telo = False
    idx = 0
    cnt = 0
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
            is_telo = False
    contig_data = contig_data[:-1]
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
    if args.alt is None:
        PAF_FILE_PATH = [args.paf_file_path]
    else:
        PAF_FILE_PATH = [args.paf_file_path, args.alt]
    
    TELOMERE_INFO_FILE_PATH = args.telomere_bed_path
    CHROMOSOME_INFO_FILE_PATH = args.reference_fai_path
    PREPROCESSED_PAF_FILE_PATH = (PAF_FILE_PATH[0] +'.ppc.paf')
    chr_len = find_chr_len(CHROMOSOME_INFO_FILE_PATH)
    telo_data = import_telo_data(TELOMERE_INFO_FILE_PATH, chr_len)
    telo_dict = defaultdict(list)
    for _ in telo_data:
        telo_dict[_[0]].append(_[1:])


    contig_data = import_data(PAF_FILE_PATH[0])

    node_label = label_node(contig_data, telo_dict)

    telo_preprocessed_contig, report_case, telo_connect_info = preprocess_telo(contig_data, node_label)

    new_contig_data = []

    for i in telo_preprocessed_contig:
        temp_list = contig_data[i]
        if i in telo_connect_info:
            temp_list.append(telo_connect_info[i])
        else:
            temp_list.append("0")
        new_contig_data.append(temp_list)

    with open("telo_preprocess_contig.txt", "wt") as f:
        for i in new_contig_data:
            for j in i:
                print(j, end="\t", file=f)
            print("", file=f)
    
    new_node_label = label_node(new_contig_data, telo_dict)
    ref_qry_ratio = calc_ratio(new_contig_data)
    preprocess_result, preprocess_contig_type, preprocess_terminal_nodes = preprocess_contig(new_contig_data, new_node_label, ref_qry_ratio)
    new_contig_data = new_contig_data[0:-1]
    contig_data_size = len(new_contig_data)
    cnt = 0
    with open(PREPROCESSED_PAF_FILE_PATH, "wt") as f:
        for i in range(contig_data_size):
            if new_contig_data[i][CTG_NAM] in preprocess_result:
                for j in new_contig_data[i][:9]: 
                    print(j, end="\t", file=f)
                print(preprocess_contig_type[new_contig_data[i][CTG_NAM]], end="\t", file=f)
                print(preprocess_terminal_nodes[new_contig_data[i][CTG_NAM]][0], end="\t", file=f)
                print(preprocess_terminal_nodes[new_contig_data[i][CTG_NAM]][1], end="\t", file=f)
                print(new_node_label[i][0]+"\t"+new_node_label[i][1], end="\t", file=f)
                print(new_contig_data[i][CTG_CONTEL], end = "\t", file=f)
                print('0.', new_contig_data[i][CTG_CNT], sep = "", file=f)
                cnt+=1


    # alt process
    if args.alt != None:
        contig_data = import_data(PAF_FILE_PATH[1])

        node_label = label_node(contig_data, telo_dict)

        telo_preprocessed_contig, report_case, telo_connect_info = preprocess_telo(contig_data, node_label)

        new_contig_data = []

        for i in telo_preprocessed_contig:
            temp_list = contig_data[i]
            if i in telo_connect_info:
                temp_list.append(telo_connect_info[i])
            else:
                temp_list.append("0")
            new_contig_data.append(temp_list)

        with open("alt_telo_preprocess_contig.txt", "wt") as f:
            for i in new_contig_data:
                for j in i:
                    print(j, end="\t", file=f)
                print("", file=f)
        
        new_node_label = label_node(new_contig_data, telo_dict)
        ref_qry_ratio = calc_ratio(new_contig_data)
        preprocess_result, preprocess_contig_type, preprocess_terminal_nodes = alt_preprocess_contig(new_contig_data, new_node_label, ref_qry_ratio)
        new_contig_data = new_contig_data[0:-1]
        contig_data_size = len(new_contig_data)
        bias = cnt
        with open(PREPROCESSED_PAF_FILE_PATH, "a") as f:
            for i in range(contig_data_size):
                if new_contig_data[i][CTG_NAM] in preprocess_result:
                    for j in new_contig_data[i][:9]: 
                        print(j, end="\t", file=f)
                    print(preprocess_contig_type[new_contig_data[i][CTG_NAM]], end="\t", file=f)
                    print(bias+preprocess_terminal_nodes[new_contig_data[i][CTG_NAM]][0], end="\t", file=f)
                    print(bias+preprocess_terminal_nodes[new_contig_data[i][CTG_NAM]][1], end="\t", file=f)
                    print(new_node_label[i][0]+"\t"+new_node_label[i][1], end="\t", file=f)
                    print(new_contig_data[i][CTG_CONTEL], end = "\t", file=f)
                    print('1.', new_contig_data[i][CTG_CNT], sep="", file=f)
                    

if __name__ == "__main__":
    main()