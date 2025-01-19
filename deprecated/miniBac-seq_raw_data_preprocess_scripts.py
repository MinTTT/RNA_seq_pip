# -*- coding: utf-8 -*-

"""

 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""

# Built-in/Generic Imports
import os

import subprocess as sbps
# […]

# […]
from threading import Thread, active_count
from time import sleep
def run_cmd(cmd):
    cwd = os.getcwd()
    stat = sbps.Popen(cmd, shell=True, cwd=cwd)
    return stat

# Own modules
#%%

"""
folder structure
--------- work folder -----------------
   |_____Raw data folder
      |_____ folders containing *.fastq.gz
   |_____Trimmed data folder
      |_____ *.fastq.gz
"""

cpu_num = 6

# raw data folder
raw_data_folder = r'./20240407_RNA-seq/PE_raw_data'
folders = os.listdir(raw_data_folder)
samples = [fold.split('-')[-1] for fold in folders]
# crate trimmed data folder
savdir = './20240407_RNA-seq/Trimmed_seq'
if not os.path.exists(savdir):
    os.makedirs(savdir)

workers = []
for filer_i, folder in enumerate(folders):

    # folder = folders[0]
    folder_dir = os.path.join(raw_data_folder, folder)
    print(f"FASTQ folder scan: {folder_dir}")
    files = [ file.name for file in os.scandir(folder_dir) if file.is_file()]
    suffiex = '.fastq.gz'
    readsfile_name = [file for file in files if file[-len(suffiex):] == suffiex]
    for file in readsfile_name:
        if file[:-len(suffiex)].split('_')[-1] == 'R1':
            read1 = os.path.join(folder_dir,file)
        elif file[:-len(suffiex)].split('_')[-1] == 'R2':
            read2 = os.path.join(folder_dir,file)
    sample_name = samples[filer_i]
    output_read1 = os.path.join(savdir, sample_name + "_R1.fastq.gz")
    output_read2 = os.path.join(savdir, sample_name + "_R2.fastq.gz")

    # 1. cut the sequence of the flank sequence (insertions the smaller than 150 bp).
    # 2. cut the fixed pooling barcode in read 1.

    # =================== cut adapter specified GTGTGAA ================================
    # cutadapt_command = f'cutadapt -j {cpu_num} ' \
    #                     + '-a CTGTCTCTTATACACATCTCCGAGCCCACGAGACTCGCCTTAATCTCGTATGCCGTCTTCTGCTTGX ' \
    #                     + '-A TTCACACAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGGCTATAGTGTAGATCTCGGTGGTCGCCGTATCATTX ' \
    #                     + f'-q 20,20 --minimum-length 100:100 --max-n 3 ' \
    #                     + '--interleaved ' \
    #                     + f'--json trim_r1_report_{sample_name}.json ' \
    #                     + f'{read1} {read2} | ' \
    #                     + f'cutadapt -j {cpu_num} ' \
    #                     + '-g "^GTGTGAA" --interleaved ' \
    #                     + f'--json trim_r2_report_{sample_name}.json ' \
    #                     + f'-o {output_read1} -p {output_read2}  -' \
    # =================== cut fixed length ================================
    log_r1 = os.path.join(savdir, 'trimmed_log', f'trim_r1_report_{sample_name}.json')
    log_r2 = os.path.join(savdir, 'trimmed_log', f'trim_r2_report_{sample_name}.json')
    log_folder = os.path.join(savdir, 'trimmed_log')
    if not os.path.exists(log_folder):
        os.makedirs(log_folder)
    cutadapt_command = f'cutadapt -j {cpu_num} ' \
                        + '-a CTGTCTCTTATACACATCTCCGAGCCCACGAGACTCGCCTTAATCTCGTATGCCGTCTTCTGCTTGX ' \
                        + '-A TTCACACAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGGCTATAGTGTAGATCTCGGTGGTCGCCGTATCATTX ' \
                        + f'-q 20,20 --minimum-length 100:100 --max-n 3 --pair-filter=any ' \
                        + '--interleaved ' \
                        + f'--json {log_r1} ' \
                        + f'{read1} {read2} | ' \
                        + f'cutadapt -j {cpu_num} ' \
                        + '-u 7 --interleaved ' \
                        + f'--json {log_r2} ' \
                        + f'-o {output_read1} -p {output_read2}  -' \

    print(cutadapt_command)

    workers.append(Thread(target=run_cmd, args=(cutadapt_command, )))

max_active_num = int(64/cpu_num)
for worker in workers:
    while active_count() >= max_active_num:
        sleep(5)
    worker.start()

for worker in workers:
    worker.join()

print('All sequences are trimmed.')