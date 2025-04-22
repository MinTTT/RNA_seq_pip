# -*- coding: utf-8 -*-

"""

 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""
# %%
# Built-in/Generic Imports
import os
import re
from time import sleep
import pandas as pd
import subprocess as sbps
from seq_utility import find_fq
from RNA_seq_analyzer import RNASeqAnalyzer
# import _thread as thread
from threading import Thread


def stat_thread(obj: RNASeqAnalyzer, thread_index: int, thread_exit):
    obj.counts_statistic()
    thread_exit[thread_index] = True
    return None

# ANSI color codes
RED = "\033[31m"
GREEN = "\033[32m"
YELLOW = "\033[33m"
BLUE = "\033[34m"
RESET = "\033[0m"  # Resets color to default

# %%
# gff_ps = r'./annotation_file/L3_strain/L3_strain.gff'
# fasta_ps = r'./annotation_file/L3_strain/L3_strain.fa'
gff_ps = r'./annotation_file/xcd001_reference/xcd001.1.gff'
fasta_ps = r'./annotation_file/xcd001_reference/xcd001.1.fa'
fastq_dir = r'/media/fulab/fulab-nas/chupan/fulab_zc_1/seq_data/20241201_RNA-seq/cleanedData'
output_dir = r'/media/fulab/fulab-nas/chupan/fulab_zc_1/seq_data/20241201_RNA-seq/outputData'
# fastq_dir = r'/media/fulab/fulab-nas/chupan/fulab_zc_1/seq_data/20250409_RNA-seq/BGI_seq/cleaned_data'
# output_dir = r'/media/fulab/fulab-nas/chupan/fulab_zc_1/seq_data/20250409_RNA-seq/BGI_seq/'
prefix = ''
threading_max = 3



seq_file_suffix = ['.fastq.gz', '.fq.gz', '.fastq', '.fq']
# pause the program and check the folder, press Enter for continue
while True:
    input_folder = input(f'{GREEN}Please check the folder that contains seq data (press Enter to continue or specify dir): \n{RESET}{fastq_dir}\n')
    if input_folder == '':
        break
    else:
        # check input is a folder
        if os.path.isdir(input_folder):
            fastq_dir = input_folder
            print(f'{GREEN}Folder is set to: {RESET}{fastq_dir}')
            break
        else:
            print(f'{GREEN}Please press Enter to continue.{RESET}')

# check the save folder
# print(f'Please check the save folder (press Enter to continue): \n{savdir}')
while True:
    input_savdir = input(f'{GREEN}Please check the folder for saving results (press Enter to continue or specify dir):{RESET}\n {output_dir}')

    if input_savdir == '':
        break
    else:
        # check input is a folder
        if os.path.isdir(input_savdir):
            savdir = input_savdir
            print(f'{GREEN}Folder is set to:{RESET} \n{savdir}')
            break
        else:
            print(f'{GREEN}The folder specified is not exist, will create it or not (y/n)?{RESET}')
            while True:
                create_folder = input(f'{GREEN}Please press Enter to continue or specify (y: create folder; n: exit): {RESET}')
                if create_folder == 'y':
                    os.makedirs(input_savdir)
                    print(f'Folder is set to: {output_dir}')
                    break
                elif create_folder == 'n':
                    print('The folder is not created, exit.')
                    exit(0)
                else:
                    pass

thread_init = 0
sample_names_rep = []
fastq_files = [file for file in os.listdir(fastq_dir)
               if os.path.isfile(os.path.join(fastq_dir, file))]

samples_dict = find_fq(fastq_dir)
samples = list(samples_dict.keys())
samples.sort()

for i, sample in enumerate(samples):
    print(f'{GREEN}[{i+1}] Sample name: {RESET} {sample}')
# let user check the sample names, run all samples by default, else user can choose the samples need to run, samples may
# be seperated by space or comma
while True:
    input_samples = input(f'{GREEN}Please check the sample names (press Enter to continue or specify samples): \n{RESET}')
    if input_samples == '':
        break
    else:

        # check input is a folder
        if len(input_samples.split(',')) > 0:
            samples = [sample.strip() for sample in input_samples.split(',')]
            print(f'{GREEN}Samples are set to: {RESET}{samples}')
            break
        else:
            print(f'{GREEN}Please press Enter to continue.{RESET}')

# ================== RNA-seq analysis Start ==================
thread_exit = []
statistic_list = []
for sample in samples:
    read1 = samples_dict[sample]['R1']
    read2 = samples_dict[sample]['R2']
    process_pip = RNASeqAnalyzer(prefix + sample, seq_ps1=read1, seq_ps2=read2, bowtie_pars={"-p": 60},
                                 output_dir=output_dir,
                                 ref_ps=fasta_ps, gff_ps=gff_ps)
    # mapping
    process_pip.seq_data_align()
    thread_exit.append(False)
    # statistic
    # thread.start_new_thread(stat_thread, (process_pip, thread_init, thread_exit))
    statistics_work = Thread(target=stat_thread, args=(process_pip, thread_init, thread_exit))
    statistics_work.start()
    statistic_list.append(statistics_work)
    thread_init += 1
    while sum([True if (sg is False) else False for sg in thread_exit]) >= threading_max:
        sleep(10)

# while False in thread_exit:
#     sleep(5)
for statistic in statistic_list:
    statistic.join()
# ================= RNA-seq analysis End ==================

# ================== Collect statistics Start ==================
# copy statistics files to the output directory, and summarize the alignment quality of all samples
"""
Read Log file
Line: "Input files DNA, FASTA:, "
extract one line after it  -> reference genome
Line: "27008777 reads; of these:"
extract 27008777 as reads number
Line: "87.55% overall alignment rate"
extract 87.55% as alignment rate

| sample name | fastq file 1 | fastq file 2 | fasta file name | gff file name | reads number | alignment rate|
"""
sample_name = []
fastq_1 = []
fastq_2 = []
fasta_file_name = []
gff_file_name = []
reads_num_list = []
alignment_rate_list = []
reads_num_map = 'reads; of these:'
alignment_rate_map = 'overall alignment rate'

for sample in samples:
    cmd_cp = f"cp {os.path.join(output_dir, prefix + sample + '_output', prefix + sample + '.expression_statistic.csv')} " + \
             f"{os.path.join(output_dir, prefix + sample + '.expression_statistic.csv')}"
    print(f'Collect {prefix + sample} statistics.')
    sbps.run(cmd_cp, shell=True)
    # read log file
    read_log_file = os.path.join(output_dir, prefix + sample + '_output', prefix + sample + '.log')
    with open(read_log_file, 'r') as f:
        lines = f.readlines()
        reads_num, alignment_rate = None, None
        for line in lines:
            if reads_num_map in line:
                # 27008777 reads; of these:
                reads_num = re.search(r'\d+', line)
                if reads_num:
                    reads_num = reads_num.group()
                    reads_num = int(reads_num)
            elif alignment_rate_map in line:
                # alignment rate
                alignment_rate = re.search(r'\d+.+%', line)
                if alignment_rate:
                    alignment_rate = alignment_rate.group()
        reads_num_list.append(reads_num)
        alignment_rate_list.append(alignment_rate)
        sample_name.append(prefix+sample)
        fastq_1.append(os.path.basename(samples_dict[sample]['R1']))
        fastq_2.append(os.path.basename(samples_dict[sample]['R2']))
        fasta_file_name.append(os.path.basename(fasta_ps))
        gff_file_name.append(os.path.basename(gff_ps))
# create a dataframe
# | sample name | fastq file 1 | fastq file 2 | fasta file name | gff file name | reads number | alignment rate|
align_table = pd.DataFrame({'sample name': sample_name,
                            'fastq file 1': fastq_1,
                            'fastq file 2': fastq_2,
                            'fasta file name': fasta_file_name,
                            'gff file name': gff_file_name,
                            'reads number': reads_num_list,
                            'alignment rate': alignment_rate_list})
# write to csv file
today = pd.Timestamp.now().strftime('%Y%m%d-%H%M')
align_table.to_csv(os.path.join(output_dir, f'{today}_alignment_statistic.csv'), index=False)
# =================== Collect statistics End ==================