
import ast
import pickle
import os
import argparse

from collections import defaultdict
from collections import Counter
from datetime import datetime

import networkx as nx
from concurrent.futures import ProcessPoolExecutor, as_completed

from tqdm import tqdm

CTG_NAM = 0
CTG_LEN = 1
CTG_STR = 2
CTG_END = 3
CTG_DIR = 4
CHR_NAM = 5
CHR_LEN = 6
CHR_STR = 7
CHR_END = 8
CTG_TYP = 9
CTG_STRND = 10
CTG_ENDND = 11
CTG_TELCHR = 12
CTG_TELDIR = 13

CHROMOSOME_COUNT = 23
K = 1e3
graph_num = 1000

def import_data(file_path : str) -> dict :
    graph_file = open(file_path, "r")
    graph_adjacency = {}
    cnt = 0
    for curr_edge in graph_file:
        l, r = curr_edge.split(":")
        r.lstrip()
        r.rstrip(',')
        r = ast.literal_eval('['+r+']')
        if cnt==0:
            cnt+=1
        l = ast.literal_eval(l)
        graph_adjacency[l] = r

    return graph_adjacency

def import_data2(file_path : str) -> list :
    paf_file = open(file_path, "r")
    contig_data = []
    for curr_contig in paf_file:
        temp_list = curr_contig.split("\t")
        int_induce_idx = [1, 2, 3, 6, 7, 8, 9, 10, 11]
        for i in int_induce_idx:
            temp_list[i] = int(temp_list[i])
        contig_data.append(tuple(temp_list))
    return contig_data

# Define a function to process a single pair
# def process_pair(args):
#     st, nd, G = args
#     local_path_list = []
#     for st_b, nd_b in [(1, 0), (1, 1)]:
#         # G_new = G.deepcopy()
#         # 그래프 수정 (st_b, st), (nd_b, nd) => 나머지 전부 없애기
#         if nx.has_path(G, (st_b, st), (nd_b, nd)):
#             i = 0
#             for path in nx.shortest_simple_paths(G, source=(st_b, st), target=(nd_b, nd), weight='weight'):
#                 i+=1
#                 if i > graph_num:
#                     break
#                 local_path_list.append(path)
#     return local_path_list

def process_pair(args):
    st, nd, G = args
    local_path_list = []
    for st_b, nd_b in [(1, 0), (1, 1)]:
        # G_new = G.deepcopy()
        # 그래프 수정 (st_b, st), (nd_b, nd) => 나머지 전부 없애기
        if nx.has_path(G, (st_b, st), (nd_b, nd)):
            local_path_list.append(nx.shortest_path(G, source=(st_b, st), target=(nd_b, nd), weight='weight'))
    return local_path_list

def main():
    parser = argparse.ArgumentParser(description="Find path for chrmosome to chrmosome.")
    
    # 위치 인자 정의
    parser.add_argument("ppc_paf_file_path", 
                        help="Path to the preprocessed PAF file.")
    parser.add_argument("graph_file_txt", 
                        help="Path to the graph text file.")

    args = parser.parse_args()

    graph_data = args.graph_file_txt
    PREPROCESSED_PAF_FILE_PATH = args.ppc_paf_file_path
    PREFIX = os.path.basename(PREPROCESSED_PAF_FILE_PATH).split('.')[0] + '.path.' + datetime.now().strftime("%H:%M:%S")

    graph_adjacency = import_data(graph_data)
    contig_data = import_data2(PREPROCESSED_PAF_FILE_PATH)
    print(len(contig_data))

    G = nx.DiGraph()
    for i in range(len(contig_data) + CHROMOSOME_COUNT*2):
        G.add_node((0, i)) # Backward 방향으로 가는 i번째 노드
        G.add_node((1, i)) # Forward 방향으로 가는 i번째 노드 

    cnt = 0
    for node in graph_adjacency:
        for edge in graph_adjacency[node]:
            if cnt==0:
                cnt+=1
                
            G.add_weighted_edges_from([(node, tuple([edge[0], edge[1]]), edge[2])])

    chr_corr = {}
    chr_rev_corr = {}
    total_contig_count = len(contig_data)
    for i in range(1, CHROMOSOME_COUNT):
        chr_corr['chr'+str(i)+'f'] = total_contig_count + i - 1
        chr_rev_corr[total_contig_count + i - 1] = 'chr'+str(i)+'f'
    chr_corr['chrXf'] = total_contig_count + CHROMOSOME_COUNT - 1
    chr_rev_corr[total_contig_count + CHROMOSOME_COUNT - 1] = 'chrXf'
    for i in range(1, CHROMOSOME_COUNT):
        chr_corr['chr'+str(i)+'b'] = total_contig_count + CHROMOSOME_COUNT + i - 1
        chr_rev_corr[total_contig_count + CHROMOSOME_COUNT + i - 1] = 'chr'+str(i)+'b'
    chr_corr['chrXb'] = total_contig_count + 2*CHROMOSOME_COUNT - 1
    chr_rev_corr[total_contig_count + 2*CHROMOSOME_COUNT - 1] = 'chrXb'


    print(G)

    nodes_in_each_contig = {}
    for i in contig_data:
        nodes_in_each_contig[i[0]] = [i[CTG_STRND], i[CTG_ENDND]]
    cnt = 0
    result = defaultdict(list)
    for i in nodes_in_each_contig:
        for j in range(nodes_in_each_contig[i][0], nodes_in_each_contig[i][1]+1):
            result[i].append(j)

    total_contig_count = len(contig_data)

    tar_node_list = list(range(total_contig_count, total_contig_count + CHROMOSOME_COUNT*2))
    # with open('data/HCC1954.telo.csv', 'r') as f:
    #     for l in f:
    #         ctg, tel_dir = l[:-1].split(',')
    #         tar_ind = 0
    #         if result[ctg]:
    #             tar_node_list.append(result[ctg][tar_ind])
    #         else:
    #             print(f'Contig : {ctg} is gone!')

    # test = [('15', '1'), ('16', '6'), ('6', '4'), ('1', '8'), ('11', '5'), ('21', '12'), ('4', '6'), ('19', '21'), ('2', '20')]
    # args_list = []
    # for s, t in test:
    #     for j in ('f', 'b'):
    #         for k in ('f', 'b'):
    #             args_list.append((chr_corr['chr'+s+j], chr_corr['chr'+t+k], G))
    args_list = []
    for i in tar_node_list:
        for j in tar_node_list:
            if i < j:
                args_list.append((i, j, G))

    # Use ProcessPoolExecutor for multiprocessing
    path_list = []
    path_counter = Counter()
    with ProcessPoolExecutor(max_workers=100) as executor:
        futures = [executor.submit(process_pair, args) for args in args_list]
        # Collect results using tqdm for progress bar
        for future in tqdm(as_completed(futures), total=len(futures)):
            want_to_print = future.result()
            for path in want_to_print:
                path_compress= (path[0][1], path[-1][1])
                path_counter[(path[0][1], path[-1][1])]+=1
                """
                first_chr = contig_data[path[0][1]][CHR_NAM]
                second_chr = contig_data[path[-1][1]][CHR_NAM]
                """
                folder_name = f"{PREFIX}/{chr_rev_corr[path_compress[0]]}_{chr_rev_corr[path_compress[1]]}"
                file_name = folder_name + f"/{path_counter[path_compress]}.paf"
                file_name2 = folder_name + f"/{path_counter[path_compress]}.index.txt"
                """
                if chr_counter['chr1'] >= 5*1e6 and chr_counter['chrX']>=5*1e6:
                    print(file_name)
                """
                os.makedirs(folder_name, exist_ok=True)
                f = open(file_name, "wt")
                g = open(file_name2, "wt")
                for nodes in path:
                    try:
                        f.write("\t".join(map(str, contig_data[nodes[1]])))
                    except:
                        f.write(chr_rev_corr[nodes[1]]+"\n")
                    g.write(" ".join(map(str, nodes))+"\n")

                f.close()
                g.close()

    # Save the path_list to a file
    # with open('path_list.pkl', 'wb') as f:
    #     pickle.dump(path_list, f)

    print(f"Completed processing. Total paths found: {len(path_list)}")

if __name__ == "__main__":
    main()