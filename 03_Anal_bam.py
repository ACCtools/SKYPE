import os
import sys
import logging
import argparse
import pickle as pkl

from juliacall import Main as jl

logging.basicConfig(
    format='%(asctime)s %(levelname)s:%(message)s',
    level=logging.INFO,
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logging.info("03_anal_bam start")


parser = argparse.ArgumentParser(description="SKYPE depth analysis")

parser.add_argument("prefix", 
                    help="Pefix for pipeline")

parser.add_argument("read_bam_loc", 
                    help="read_aln_loc")

parser.add_argument("reference_fai_path", 
                    help="Path to the chromosome information file.")

parser.add_argument("--progress", 
                    help="Show progress bar", action='store_true')

args = parser.parse_args()

PREFIX = args.prefix
read_bam_loc = args.read_bam_loc
reference_fai_path = args.reference_fai_path

with open(f"{PREFIX}/nclose_cord_list.pkl", "rb") as f:
    nclose_cord_list, nclose_idx_corr, total_nclose_cord_list_contig_name = pkl.load(f)

jl.nclose_cord_vec = jl.Vector[jl.Vector[jl.Any]]()
for l in nclose_cord_list:
    jl.push_b(jl.nclose_cord_vec, jl.Vector[jl.Any](l))

is_progress_bar = sys.stdout.isatty() or args.progress
jl.include(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'anal_bam.jl'))
ac_nclose_cnt_list, wa_nclose_cnt_list = jl.anal_bam(read_bam_loc, reference_fai_path, jl.nclose_cord_vec, is_progress_bar)

ac_nclose_cnt_dict = dict(list(ac_nclose_cnt_list))
wa_nclose_cnt_dict = dict(list(wa_nclose_cnt_list))

nclose2cov = dict()
for k, v in ac_nclose_cnt_dict.items():
    if k not in wa_nclose_cnt_dict:
        nclose2cov[nclose_idx_corr[k]] = v

with open(f"{PREFIX}/nclose2cov.pkl", "wb") as f:
    pkl.dump(nclose2cov, f)

with open(f"{PREFIX}/nclose_cord_list.txt", "wt") as f2:
    for l in total_nclose_cord_list_contig_name:
        print(*(l + [ac_nclose_cnt_dict.get(l[-2], -1), wa_nclose_cnt_dict.get(l[-2], -1)]),
              sep="\t", file=f2)

with open(f"{PREFIX}/bam_nclose_cnt.pkl", "wb") as f:
    pkl.dump((ac_nclose_cnt_dict, wa_nclose_cnt_dict), f)
