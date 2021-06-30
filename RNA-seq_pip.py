# -*- coding: utf-8 -*-

"""

 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""

# Built-in/Generic Imports
import os
import sys
from time import sleep
import subprocess as sbps
# […]

# Libs
import pandas as pd
import numpy as np  # Or any other
from RNA_seq_analyzer import RNASeqAnalyzer
import _thread as thread

# […]

# Own modules

save_dir = r'/home/fulab/data/andong_rna'

sample_dir = r'/home/fulab/data/andong_rna/Clean'

sample_names = os.listdir(sample_dir)

sample_pars_list = [dict(name=sample, fasta_file=os.listdir(os.path.join(sample_dir, sample))) for sample in
                    sample_names]

genome_fa = r'./example_data/annotation_file/GCA_000005845.2_ASM584v2_genomic.fasta'
gff = r'./example_data/annotation_file/GCA_000005845.2_ASM584v2_genomic.gff'


def stat_thread(obj: RNASeqAnalyzer, thread_index: int):
    obj.counts_statistic()
    thread_exit[thread_index] = True
    return None


thread_exit = []
thread_init = 0
for sample in sample_pars_list:
    sample_fa = sample['fasta_file']
    sample_name = sample['name']
    sample_replicates = list(set([fafile.strip(fafile.split('_')[-1])[0:-1] for fafile in sample_fa]))
    for replicate in sample_replicates:

        print(f'[{replicate}] -> Processing Now')
        fa1 = os.path.join(sample_dir, sample_name, replicate + '_1.fq.gz')
        fa2 = os.path.join(sample_dir, sample_name, replicate + '_2.fq.gz')
        process_pip = RNASeqAnalyzer(replicate, genome_fa, gff_ps=gff, seq_ps1=fa1, seq_ps2=fa2,
                                     bowtie_pars={"-p": 64}, output_dir=save_dir)
        process_pip.seq_data_align()
        thread_exit.append(False)
        thread.start_new_thread(stat_thread, (process_pip, thread_init))
        thread_init += 1
        while sum([True if sg == False else False for sg in thread_exit]) >= 3:
            sleep(10)


while False in thread_exit:
    sleep(5)

for sample in sample_pars_list:
    sample_fa = sample['fasta_file']
    sample_name = sample['name']
    sample_replicates = list(set([fafile.strip(fafile.split('_')[-1])[0:-1] for fafile in sample_fa]))
    for replicate in sample_replicates:
        cmd_cp = f"cp {os.path.join(save_dir, replicate + '_output', replicate + '.counts_statistic.csv')} " + \
                 f"{os.path.join(save_dir, replicate + '.counts_statistic.csv')}"
        sbps.run(cmd_cp, shell=True)
