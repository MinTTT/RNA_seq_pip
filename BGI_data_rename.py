# -*- coding: utf-8 -*-

"""

 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""

# Built-in/Generic Imports
import os
import sys
# […]

# Libs
import pandas as pd
import numpy as np  # Or any other
# […]
import subprocess as sbps
import json
# Own modules




data_fold = r'/media/fulab/fulab-nas/chupan/fulab_zc_1/seq_data/2022_20_RNA_seq/PE'
target_dir = r'/media/fulab/fulab-nas/chupan/fulab_zc_1/seq_data/2022_20_RNA_seq/PE_data_rename'

if not os.path.exists(target_dir):
    os.mkdir(target_dir)

"""
------------------
|____________ PE
    |________ Sample_fold
        |____ *.fq.gz
        
"""
sample_folders = [os.path.join(data_fold, fold) for fold in os.listdir(data_fold)
               if os.path.isdir(os.path.join(data_fold, fold))]

origin_name_new_name = {}
for folder in sample_folders:
    fasta_files = [file for file in os.listdir(folder)
                   if os.path.isfile(os.path.join(folder, file)) and file[-5:] == 'fq.gz']
    for file in fasta_files:
        file_name_elem = file.split('.')[0].split('_')
        if file_name_elem[-1] == '1':
            read_index = 'R1'
        else:
            read_index = 'R2'
        file_new_name = "_".join([file_name_elem[0], file_name_elem[4], read_index]) + '.fastq.gz'
        origin_name_new_name[file] = file_new_name
        cmd = f"mv {os.path.join(folder, file)} {os.path.join(target_dir, file_new_name)}"
        print(cmd)
        status0 = sbps.Popen(cmd, shell=True, stdout=sbps.PIPE, cwd=os.getcwd())

with open(os.path.join(target_dir, 'origin_name_new_name.json'), 'w') as file:
    json.dump(origin_name_new_name, file)

