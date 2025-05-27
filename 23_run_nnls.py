import os
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

jl.include(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'run_nnls.jl'))
error, b_norm, weights_jl, predict_B_jl, b_start_ind = jl.run_nnls(PREFIX, min(THREAD, psutil.cpu_count(logical=False)))

logging.info(f'Error : {round(error, 4)}')
logging.info(f'Norm error : {round(error / b_norm, 4)}')

weights = np.asarray(weights_jl)
np.save(f'{PREFIX}/weight.npy', weights)

predict_B = np.asarray(predict_B_jl)[b_start_ind:]
np.save(f'{PREFIX}/predict_B.npy', predict_B)
