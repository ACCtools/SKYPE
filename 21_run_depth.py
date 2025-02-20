import os
import sys
import glob
import pandas as pd
import subprocess

from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor, as_completed

PATH_FILE_FOLDER = f"{sys.argv[1]}/20_depth/"
chr_chr_folder_path = glob.glob(PATH_FILE_FOLDER+"*")
DEPTH_WINDOW=100 * 1e3

DEPTH_THREAD=2
TOTAL_THREAD=128


def get_paf_run(paf_loc):
    paf_base = os.path.splitext(paf_loc)[0]
    result = subprocess.run(['./PanDepth/bin/pandepth', '-w', str(DEPTH_WINDOW),'-t', str(DEPTH_THREAD), '-i', paf_loc, '-o', paf_base],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, capture_output=False)
    
    if result.returncode != 0:
        print(paf_loc, 'Error!')
    
    # cov_file_loc = paf_base + '.win.stat.gz' 
    # df = pd.read_csv(cov_file_loc, compression='gzip', comment='#', sep='\t', names=['chr', 'st', 'nd', 'length', 'covsite', 'totaldepth', 'cov', 'meandepth'])

    # return paf_loc


with ProcessPoolExecutor(max_workers=int(TOTAL_THREAD/DEPTH_THREAD)) as executor:
    futures = []
    for folder_path in chr_chr_folder_path:
        paf_paths = glob.glob(folder_path + "/*.paf")
        for paf_loc in paf_paths:
            futures.append(executor.submit(get_paf_run, paf_loc))
    
    # 제출된 작업들이 완료될 때까지 진행 상황을 tqdm으로 표시합니다.
    for future in tqdm(as_completed(futures), total=len(futures), desc='Run PanDepth for each path file', disable=not sys.stdout.isatty()):
        pass
