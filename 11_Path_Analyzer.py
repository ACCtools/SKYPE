import glob
import pandas as pd
from collections import defaultdict
from collections import Counter

START_NODE = "chr1f"

END_NODE = "chr1b"

QUERY = ""

CHROMOSOME_COUNT = 23

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

INF = int(1e9)

path_path_list = glob.glob(f"./HCC1954.path.15:50:07/{START_NODE}_{END_NODE}/*.paf")

chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

path_df = pd.DataFrame(columns=["index"] + chromosomes)

def overlap_checker(node_a : list, node_b : list) -> int :
    # check if overlap
    if max(node_a[0], node_b[0]) < min(node_a[1], node_b[1]):
        return 1   
    else:
        return 0
cnt=1
for _ in path_path_list:
    with open(_, "r") as f:
        chr_range_dict = defaultdict(list)
        path_components_dict = Counter()
        for contig_data in f:
            contig_data = contig_data.split('\t')
            if len(contig_data)==1:
                continue
            chr_range_dict[contig_data[CHR_NAM]].append([int(contig_data[CHR_STR]), int(contig_data[CHR_END])])
        for i in chr_range_dict:
            sorted_range = sorted(chr_range_dict[i])
            merge_range = []
            for _ in sorted_range:
                if merge_range==[]:
                    merge_range.append(_)
                else:
                    if overlap_checker(merge_range[-1], _):
                        merge_range[-1] = [min(merge_range[-1][0], _[0]), max(merge_range[-1][1], _[1])]
                    else:
                        merge_range.append(_)
            net_length = 0
            for j in merge_range:
                net_length+=j[1]-j[0]
            real_length = merge_range[-1][1]-merge_range[0][0]
            path_components_dict[i] = net_length # 이 부분 질문할것;
        record = [cnt]
        for i in chromosomes:
            record.append(path_components_dict[i])
        path_df.loc[len(path_df)] = record
    cnt+=1
    f.close()


user_query = input() # or "QUERY" 변수 활용
    
output_file = "query_result.csv"

filtered_df = path_df.query(user_query)
filtered_df.to_csv(output_file, index=False)
print(f"Query saved at {output_file}")