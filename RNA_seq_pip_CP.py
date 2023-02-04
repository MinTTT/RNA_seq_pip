# -*- coding: utf-8 -*-

"""

 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""
#%%
# Built-in/Generic Imports
import os
import sys
from time import sleep
import subprocess as sbps
import json
# […]

# Libs
import pandas as pd
import numpy as np  # Or any other
from RNA_seq_analyzer import RNASeqAnalyzer
import _thread as thread


def stat_thread(obj: RNASeqAnalyzer, thread_index: int):
    obj.counts_statistic()
    thread_exit[thread_index] = True
    return None

#%%
# […]
# sample info:
# samples_pars = {"1321_1": dict(ref_ps=r"./example_data/annotation_file/NH3.23.fasta",
#                                gff_ps=r"./example_data/annotation_file/NH3.23.gff"),
#                 "1321_2": dict(ref_ps=r"./example_data/annotation_file/NH3.23.fasta",
#                                gff_ps=r"./example_data/annotation_file/NH3.23.gff"),
#                 "1321_3": dict(ref_ps=r"./example_data/annotation_file/NH3.23.fasta",
#                                gff_ps=r"./example_data/annotation_file/NH3.23.gff"),
#                 "1321_4": dict(ref_ps=r"./example_data/annotation_file/NH3.23.fasta",
#                                gff_ps=r"./example_data/annotation_file/NH3.23.gff"),
#                 "1321_5": dict(ref_ps=r"./example_data/annotation_file/NH3.23.fasta",
#                                gff_ps=r"./example_data/annotation_file/NH3.23.gff"),
#                 "1321_6": dict(ref_ps=r"./example_data/annotation_file/NH3.24.fasta",
#                                gff_ps=r"./example_data/annotation_file/NH3.24.gff"),
#                 "1321_7": dict(ref_ps=r"./example_data/annotation_file/NH3.24.fasta",
#                                gff_ps=r"./example_data/annotation_file/NH3.24.gff"),
#                 "1321_8": dict(ref_ps=r"./example_data/annotation_file/NH3.24.fasta",
#                                gff_ps=r"./example_data/annotation_file/NH3.24.gff"),
#                 "1321_9": dict(ref_ps=r"./example_data/annotation_file/NH3.24.fasta",
#                                gff_ps=r"./example_data/annotation_file/NH3.24.gff"),
#                 "1321_10": dict(ref_ps=r"./example_data/annotation_file/NH3.24.fasta",
#                                 gff_ps=r"./example_data/annotation_file/NH3.24.gff")
#                 }
# save_dir = r"/home/fulab/data2/20210321_RNA_seq
# sample_dir = r"/home/fulab/data2/20210321_RNA_seq/CleanData"

sampleinfo = json.load(open(r'./RNASeqSampleInfo.json'))


sample_names = list(sampleinfo['samples'].keys())

sample_fa_files = [dict(name=sample, fasta_file=os.listdir(os.path.join(sampleinfo['samples'][sample][1], sample))) for sample in
                   sample_names]

# genome_fa = r'./example_data/annotation_file/CLB_strain.fa'
# gff = r'./example_data/annotation_file/CLB_strain.gff'
# adapters = ['AGATCGGAAGAGC', 'AGATCGGAAGAGC']
seq_data_suffix = '.fq.gz'

threading_max = 1
#%%
thread_exit = []
thread_init = 0
sample_names_rep = []
for sample in sample_fa_files:
    sample_fa = sample['fasta_file']  # type: list
    sample_name = sample['name']
    sample_type = sampleinfo['samples'][sample_name][0]
    sample_annotation = sampleinfo['genome_annotations'][sample_type]
    sample_replicates = list(set([fa_fl.strip(fa_fl.split('_')[-1])[0:-1] for fa_fl in sample_fa]))

    for rep_index, replicate in enumerate(sample_replicates):
        sample_name_rep = f"{sample_name}_{rep_index}"
        print(f'[{sample_name_rep}] -> Processing Now')
        fa1 = os.path.join(sampleinfo['samples'][sample_name][1], sample_name, replicate + '_1' + seq_data_suffix)  # single end seq


        fa2 = os.path.join(sampleinfo['samples'][sample_name][1], sample_name, replicate + '_2' + seq_data_suffix)
        if not os.path.exists(fa2):
           fa2 = None     

        process_pip = RNASeqAnalyzer(sample_name_rep, seq_ps1=fa1, seq_ps2=fa2, bowtie_pars={"-p": 32}, output_dir=sampleinfo['save_dir'],
                                     ref_ps=sample_annotation['fasta_ps'], gff_ps=sample_annotation['gff_ps'])
        sample_names_rep.append(sample_name_rep)
        process_pip.seq_data_align()
        thread_exit.append(False)
        thread.start_new_thread(stat_thread, (process_pip, thread_init))
        thread_init += 1
        while sum([True if sg is False else False for sg in thread_exit]) >= threading_max:
            sleep(10)


while False in thread_exit:
    sleep(5)

for sample_rep in sample_names_rep:
    # sample_fa = sample['fasta_file']
    # sample_name = sample['name']
    # sample_replicates = list(set([fafile.strip(fafile.split('_')[-1])[0:-1] for fafile in sample_fa]))
    # for rep_index, replicate in enumerate(sample_replicates):
    cmd_cp = f"cp {os.path.join(save_dir, sample_rep + '_output', sample_rep + '.expression_statistic.csv')} " + \
             f"{os.path.join(save_dir, sample_rep + '.expression_statistic.csv')}"
    sbps.run(cmd_cp, shell=True)
