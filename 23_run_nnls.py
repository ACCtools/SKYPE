import psutil
import logging
import argparse

import numpy as np

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
jl.seval("using JLD2, LinearAlgebra, SINNLS")
jl.BLAS.set_num_threads(min(THREAD, psutil.cpu_count(logical=False)))
jl.seval("""
function load_nnls_array(filename::String)
    file = jldopen(filename, "r")
    return Matrix{Float32}(file["A"]), Vector{Float32}(file["B"])
end
""")

A_jl, B_jl = jl.load_nnls_array(f'{PREFIX}/matrix.h5')

logging.info('Regression analysis is ongoing...')
H = 3600.0
weights_jl = jl.vec(jl.SI_NNLS(A_jl, B_jl,
                               total_time=12 * H,
                               epi=5))

error = jl.norm(A_jl * weights_jl - B_jl)
b_norm = jl.norm(B_jl)

logging.info(f'Error : {round(error, 4)}')
logging.info(f'Norm error : {round(error / b_norm, 4)}')

weights = np.asarray(weights_jl)
np.save(f'{PREFIX}/weight.npy', weights)
