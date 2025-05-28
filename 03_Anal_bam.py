import os
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



args = parser.parse_args()

PREFIX = args.prefix
read_bam_loc = args.read_bam_loc
reference_fai_path = args.reference_fai_path

with open(f"{PREFIX}/nclose_cord_list.pkl", "rb") as f:
    nclose_cord_list, nclose_idx_corr = pkl.load(f)

jl.nclose_cord_vec = jl.Vector[jl.Vector[jl.Any]]()
for l in nclose_cord_list:
    jl.push_b(jl.nclose_cord_vec, jl.Vector[jl.Any](l))

jl.include(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'anal_bam.jl'))
ac_nclose_cnt_list, wa_nclose_cnt_list = jl.anal_bam(read_bam_loc, reference_fai_path, jl.nclose_cord_vec)

ac_nclose_cnt_dict = dict(list(ac_nclose_cnt_list))
wa_nclose_cnt_dict = dict(list(wa_nclose_cnt_list))

nclose2cov = dict()
for k, v in ac_nclose_cnt_dict.items():
    if k not in wa_nclose_cnt_dict:
        nclose2cov[nclose_idx_corr[k]] = v

with open(f"{PREFIX}/nclose2cov.pkl", "wb") as f:
    pkl.dump(nclose2cov, f)
