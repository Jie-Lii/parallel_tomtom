# -*- coding: utf-8 -*-
'''
@ Author Li Jie
@ Version 1.0
@ FileName: benchmark.py
@ Date 2025/4/30
@ Description: 
'''

import gc
import time
import random
import numpy as np
import pandas as pd

import parallel_tomtom


def generate_acgt_motifs(n, length=10):
    bases = ['A', 'C', 'G', 'T']
    motifs = []

    for _ in range(n):
        motif = ''.join(random.choices(bases, k=length))
        motifs.append(motif)
    return motifs


if __name__ == '__main__':
    out_list = []

    for num in [100,500,1000,2500,5000]:
        tmp_out_list=[]
        motif_list = generate_acgt_motifs(num)
        for i in [1, 8, 32, 64, 128]:
            start_time = time.time()

            parallel_tomtom.run_mutiprocess_tomtom_suite(motif_list=motif_list,
                                                         open_database_path="../Data/db/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme",
                                                         save_file_path=f"../Result/benchmark/{num}_{i}.csv",
                                                         command="tomtom -no-ssc -oc . -verbosity 1 -min-overlap 5 -text -dist pearson -evalue -thresh 0.5",
                                                         num_processes=i,
                                                         header_list=None,
                                                         file_type="csv")

            end_time = time.time()
            tmp_out_list.append(end_time - start_time)
            print(end_time - start_time)
            gc.collect()
        out_list.append(tmp_out_list)

    df = pd.DataFrame(out_list, columns=['1', '8', '32', '64', '128'])
    df.to_csv("../Result/benchmark.csv",index=False)
