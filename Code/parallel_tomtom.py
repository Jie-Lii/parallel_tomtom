# -*- coding: utf-8 -*-
'''
@ author LiJ
@ version 1.0
@ date 2024 / 11 / 06 20: 35
@ Description: process pdi data
'''

import os
import shutil
import argparse
from multiprocessing import Pool

import numpy as np
import pandas as pd

dict_DnaRna = {'A': "1.0 0.0 0.0 0.0\n",
               'C': "0.0 1.0 0.0 0.0\n",
               'G': "0.0 0.0 1.0 0.0\n",
               'T': "0.0 0.0 0.0 1.0\n",
               'U': "0.0 0.0 0.0 1.0\n",
               'N': "0.25 0.25 0.25 0.25\n",
               '.': "0.25 0.25 0.25 0.25\n",
               'X': "0.25 0.25 0.25 0.25\n",
               'V': "0.333333 0.333333 0.333333 0\n",
               'H': "0.333333 0.333333 0 0.333333\n",
               'D': "0.333333 0 0.333333 0.333333\n",
               'B': "0 0.333333 0.333333 0.333333\n",
               'M': "0.5 0.5 0 0\n",
               'R': "0.5 0 0.5 0\n",
               'W': "0.5 0 0 0.5\n",
               'S': "0 0.5 0.5 0\n",
               'Y': "0 0.5 0 0.5\n",
               'K': "0 0 0.5 0.5\n"}


def split_works(works, num_processor=12):
    if len(works) < num_processor:
        return [[i] for i in works]
    else:
        avg_size = len(works) // num_processor  # 每个块的基本大小
        remainder = len(works) % num_processor  # 剩余元素数目

        result = []
        start = 0
        for i in range(num_processor):
            end = start + avg_size
            result.append(works[start:end])
            start = end
        for i in range(remainder):
            result[i].append(works[start + i])
        return result


def load_moitf_list_by_csv(open_file_path):
    file = pd.read_csv(open_file_path, encoding='utf-8').values[:, 0]
    return file


def load_motif_list_by_txt_meme(open_file_path):
    header_list = []
    motif_list = []
    with open(open_file_path, "r") as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith("MOTIF"):
            break
        header_list.append(lines[i])
        i += 1

    while i < len(lines):
        if lines[i].strip().startswith("MOTIF"):
            motif_block = []
            while i < len(lines):
                line = lines[i].strip()
                if line.startswith("MOTIF") and motif_block:
                    break
                motif_block.append(lines[i])
                i += 1
            motif_list.append(motif_block)
        else:
            i += 1
    return header_list, motif_list


def build_write_list(motif_list, header_list=None, file_type="csv"):
    write_list = []
    if file_type == "csv":
        write_list = ['MEME version 4\n', '\n',
                      'ALPHABET= ACGT\n', '\n',
                      'strands: + -\n', '\n',
                      'Background letter frequencies (from uniform):\n',
                      'A 0.25 C 0.25 G 0.25 T 0.25\n', ]
        for i in range(len(motif_list)):
            write_list.append('\n')
            write_list.append(f'MOTIF {motif_list[i]} {motif_list[i]}\n')
            write_list.append('\n')
            write_list.append(f'letter-probability matrix: alength= 4 w= {len(motif_list[i])} nsites= 1 E= 0e+0\n')
            for j in motif_list[i]:
                write_list.append(dict_DnaRna[j])
        write_list.append('\n')
    elif file_type == "txt" or file_type == "meme":
        write_list += header_list
        for i in range(len(motif_list)):
            for j in range(len(motif_list[i])):
                write_list.append(motif_list[i][j])
    return write_list


def write_tomeme(motif_list, save_file_path, header_list=None, file_type="csv"):
    # name = 'py_normall_8'
    write_list = build_write_list(motif_list, header_list, file_type)
    with open(save_file_path, 'w') as f:
        f.writelines(write_list)
        f.close()


def apart_of_motiflist(motif_list, num_apart, save_file_path='./tmp_data/meme_file/', header_list=None, file_type="csv"):
    meme_file_path_list = []
    motif_list = split_works(motif_list, num_processor=num_apart)
    for i in range(num_apart):
        write_tomeme(motif_list=motif_list[i], save_file_path=f'{save_file_path}/{i}.meme', header_list=header_list, file_type=file_type)
        meme_file_path_list.append(f'{save_file_path}/{i}.meme')
    return meme_file_path_list


def build_meme_result_path(num_file, build_file_path):
    save_file_path_list = []
    for i in range(num_file):
        try:
            os.makedirs('%s/%s' % (build_file_path, i))
        except:
            print("ERROR: file exists!")
        save_file_path_list.append("%s/%s" % (build_file_path, i))
    return save_file_path_list


def run_tomtom_command(command, open_meme_path, open_database_path, save_file_path):
    command = f"{command} {open_meme_path} {open_database_path} -oc {save_file_path} > {save_file_path}/tomtom.tsv"
    os.system(command=command)


def muti_run_tomtom_command(open_meme_path_list, open_database_path, save_file_path_list, command, show_log=False):
    if show_log == True:
        print(' TOMTOM compare begin '.center(100, '-') + '\n')

    mulpros = []

    pool = Pool(processes=len(open_meme_path_list))
    for i in range(len(open_meme_path_list)):
        res = pool.apply_async(run_tomtom_command, args=(command, open_meme_path_list[i], open_database_path, save_file_path_list[i],))
        mulpros.append(res)
    pool.close()
    pool.join()
    if show_log == True:
        print(' TOMTOM compare finish '.center(100, '-') + '\n')


def concatenate_tomtom_result(open_file_path, save_file_path, show_log=False):
    if show_log == True:
        print(' concatenating tomtom result '.center(100, '-') + '\n')
    compare_list = []
    file_path_list = os.listdir(open_file_path)
    for i in range(len(file_path_list)):
        try:
            file = pd.read_csv("%s/%s/tomtom.tsv" % (open_file_path, i), sep='\t', encoding='utf-8').values[:-3]
            compare_list.append(file)
        except:
            pass
    if len(compare_list) > 0:
        dir_path = os.path.dirname(save_file_path)
        os.makedirs(dir_path, exist_ok=True)

        compare_list = np.vstack(compare_list)
        df = pd.DataFrame(compare_list, columns=['Query_ID', 'Target_ID', 'Optimal_offset', 'p-value', 'E-value', 'q-value', 'Overlap', 'Query_consensus', 'Target_consensus', 'Orientation'])
        df.to_csv(save_file_path, index=False)
    else:
        print("There're no matching motifs!")
    if show_log == True:
        print(' concatenation completed'.center(100, '-') + '\n')


def remove_files(del_file_path):
    shutil.rmtree(path=del_file_path)


def mkfile(make_file_path):
    try:
        shutil.rmtree(path=make_file_path)
    except:
        pass
    os.mkdir(make_file_path)
    os.mkdir("%s/%s" % (make_file_path, 'meme_file'))
    os.mkdir("%s/%s" % (make_file_path, 'result_file'))


def run_mutiprocess_tomtom_suite(motif_list, open_database_path, save_file_path, command, num_processes=16, header_list=None, file_type="csv", show_log=False):
    if len(motif_list) < num_processes:
        num_processes = len(motif_list)

    mkfile(make_file_path='./tmp_data/')
    meme_file_path_list = apart_of_motiflist(motif_list=motif_list, num_apart=num_processes, header_list=header_list, file_type=file_type)
    save_file_path_list = build_meme_result_path(num_file=len(meme_file_path_list), build_file_path='./tmp_data/result_file/')
    muti_run_tomtom_command(command=command, open_meme_path_list=meme_file_path_list, open_database_path=open_database_path, save_file_path_list=save_file_path_list, show_log=False)
    concatenate_tomtom_result(open_file_path='./tmp_data/result_file/', save_file_path=save_file_path)
    remove_files(del_file_path='./tmp_data/')


def tomtom_motif(command, open_motif_file_path, open_database_path, save_file_path, num_processes=16):
    file_type = open_motif_file_path.split('.')[-1]
    header_list, motif_list = [], []
    if file_type == 'csv':
        header_list = None
        motif_list = load_moitf_list_by_csv(open_motif_file_path)
    elif file_type == 'meme' or file_type == 'txt':
        header_list, motif_list = load_motif_list_by_txt_meme(open_motif_file_path)
    run_mutiprocess_tomtom_suite(motif_list=motif_list,
                                 open_database_path=open_database_path,
                                 save_file_path=save_file_path,
                                 command=command,
                                 num_processes=num_processes,
                                 header_list=header_list,
                                 file_type=file_type)


def build_argparse():
    parser = argparse.ArgumentParser()
    parser.add_argument('-TOMTOM', type=str, help='Usage: tomtom [options] <query motif file> <target motif file>+  see https://meme-suite.org/meme/doc/tomtom.html', required=True)
    parser.add_argument('-motif_path', type=str, help='motif file path', default=None, required=True)
    parser.add_argument('-database_path', type=str, help='database file path', default=None, required=True)
    parser.add_argument('-save_path', type=str, help='save file path', default=None, required=True)
    parser.add_argument('-num_processes', type=int, help='num_processes', default=8, required=True)
    return parser


if __name__ == '__main__':
    """
    example running sample

    CSV file input
    python parallel_tomtom.py -TOMTOM "tomtom -no-ssc -oc . -verbosity 1 -min-overlap 5 -text -dist pearson -evalue -thresh 0.5" -motif_path ../Data/motif/Motifs.csv -database_path ../Data/db/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme -save_path ../Result/csv_result.csv -num_processes 2

    TXT file input
    python parallel_tomtom.py -TOMTOM "tomtom -no-ssc -oc . -verbosity 1 -min-overlap 5 -text -dist pearson -evalue -thresh 0.5" -motif_path ../Data/motif/Motifs.txt -database_path ../Data/db/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme -save_path ../Result/txt_result.csv -num_processes 2

    MEME file input
    python parallel_tomtom.py -TOMTOM "tomtom -no-ssc -oc . -verbosity 1 -min-overlap 5 -text -dist pearson -evalue -thresh 0.5" -motif_path ../Data/motif/Motifs.meme -database_path ../Data/db/JASPAR2022_CORE_vertebrates_non-redundant_v2.meme -save_path ../Result/meme_result.csv -num_processes 2
    """
    parser = build_argparse()
    args = parser.parse_args()
    tomtom_motif(args.TOMTOM, args.motif_path, args.database_path, args.save_path, args.num_processes)
