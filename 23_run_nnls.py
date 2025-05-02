import psutil
import logging
import argparse

import numpy as np
import pickle as pkl

from juliacall import Main as jl

logging.basicConfig(
    format='%(asctime)s %(levelname)s:%(message)s',
    level=logging.INFO,
    datefmt='%m/%d/%Y %I:%M:%S %p',
)
logging.info("23_run_nnls start")


parser = argparse.ArgumentParser(description="SKYPE depth analysis")

parser.add_argument("prefix", 
                    help="Pefix for pipeline")

parser.add_argument("-t", "--thread", 
                    help="Number of thread", type=int)

args = parser.parse_args()

PREFIX = args.prefix
THREAD = args.thread

# Run julia for NNLS
# jl.ENV["MKL_NUM_THREADS"] = str(min(THREAD, psutil.cpu_count(logical=False)))
# jl.seval('using MKLSparse')
# jl.seval('using MKL')

# Julia openblas setting
jl.seval("using HDF5, LinearAlgebra, SINNLS")
jl.BLAS.set_num_threads(min(THREAD, psutil.cpu_count(logical=False)))
jl.seval("""
function load_nnls_array(filename::String)
    h5open(filename, "r") do file
        return read(file["A"]), read(file["B"])
    end
end
         
function load_fail_array(filename::String)
    h5open(filename, "r") do file
        return read(file["A_fail"])
    end
end
""")

with open(f'{PREFIX}/pri_weight_data.pkl', 'rb') as f:
    init_cols, w_pri = pkl.load(f)

A_jl, B_jl = jl.load_nnls_array(f'{PREFIX}/matrix.h5')

logging.info('Regression analysis is ongoing...')

n = jl.size(A_jl)[1]
A_T = jl.eltype(A_jl)
jl.x0_ = jl.zeros(A_T, n)

jl.init_cols = jl.Vector[jl.Int64](init_cols)
jl.w_pri = jl.Vector[jl.Float64](w_pri)

jl.seval("x0_[init_cols] = w_pri")

H = 3600.0
weights_jl = jl.vec(jl.SI_NNLS(A_jl, B_jl,
                               x0_= jl.x0_,
                               total_time=24 * H,
                               restart_ratio=0.8,
                               epi=2))

predict_suc_B_jl = A_jl * weights_jl
error = jl.norm(predict_suc_B_jl - B_jl)
b_norm = jl.norm(B_jl)

A_jl = jl.nothing
jl.GC.gc()

logging.info(f'Error : {round(error, 4)}')
logging.info(f'Norm error : {round(error / b_norm, 4)}')

A_fail_jl = jl.load_fail_array(f'{PREFIX}/matrix.h5')
predict_B_jl = jl.vcat(predict_suc_B_jl, A_fail_jl * weights_jl)

weights = np.asarray(weights_jl)
np.save(f'{PREFIX}/weight.npy', weights)

predict_B = np.asarray(predict_B_jl)
np.save(f'{PREFIX}/predict_B.npy', predict_B)
