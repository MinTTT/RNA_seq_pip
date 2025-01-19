# -*- coding: utf-8 -*-

"""

 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""
# %%
# Built-in/Generic Imports
import os
from time import sleep
import subprocess as sbps

from RNA_seq_analyzer import RNASeqAnalyzer
# import _thread as thread
from threading import Thread

def stat_thread(obj: RNASeqAnalyzer, thread_index: int, thread_exit):
    obj.counts_statistic()
    thread_exit[thread_index] = True
    return None




# %%


gff_ps = r'./annotation_file/L3_strain/L3_strain.gff'
fasta_ps = r'./annotation_file/L3_strain/L3_strain.fa'
fastq_dir = r'./20240407_RNA-seq/Trimmed_seq/'
output_dir = r'./20240407_RNA-seq/'
threading_max = 3


thread_init = 0
sample_names_rep = []
fastq_files = [file for file in os.listdir(fastq_dir)
               if os.path.isfile(os.path.join(fastq_dir, file))]
samples = list(set([ '_'.join((file.replace('.fastq.gz', '').split('_')[:-1]))  for file in fastq_files]))
samples.sort()
# force run this sample
# samples = ['1126_3s']

print('''all samples' name: ''', samples)

thread_exit = []
statistic_list = []
for sample in samples:
    read1 = os.path.join(fastq_dir, f'{sample}_R1.fastq.gz')
    read2 = os.path.join(fastq_dir, f'{sample}_R2.fastq.gz')
    process_pip = RNASeqAnalyzer(sample, seq_ps1=read1, seq_ps2=read2, bowtie_pars={"-p": 64}, output_dir=output_dir,
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

for sample in samples:

    cmd_cp = f"cp {os.path.join(output_dir, sample + '_output', sample + '.expression_statistic.csv')} " + \
             f"{os.path.join(output_dir, sample + '.expression_statistic.csv')}"
    print('Collect all statistics.')
    sbps.run(cmd_cp, shell=True)
#
# for sample in sample_fa_files:
#     sample_fa = sample['fasta_file']  # type: list
#     sample_name = sample['name']
#     sample_type = sampleinfo['samples'][sample_name][0]
#     sample_annotation = sampleinfo['genome_annotations'][sample_type]
#     sample_replicates = list(set([fa_fl.strip(fa_fl.split('_')[-1])[0:-1] for fa_fl in sample_fa]))
#
#     for rep_index, replicate in enumerate(sample_replicates):
#         sample_name_rep = f"{sample_name}_{rep_index}"
#         print(f'[{sample_name_rep}] -> Processing Now')
#         fa1 = os.path.join(sampleinfo['samples'][sample_name][1], sample_name,
#                            replicate + '_1' + seq_data_suffix)  # single end seq
#
#         fa2 = os.path.join(sampleinfo['samples'][sample_name][1], sample_name, replicate + '_2' + seq_data_suffix)
#         if not os.path.exists(fa2):
#             fa2 = None
#         sample_names_rep.append(sample_name_rep)
#         # process_pip = RNASeqAnalyzer(sample_name_rep, seq_ps1=fa1, seq_ps2=fa2, bowtie_pars={"-p": 32}, output_dir=sampleinfo['save_dir'],
#         #                              ref_ps=sample_annotation['fasta_ps'], gff_ps=sample_annotation['gff_ps'])
#
#         # process_pip.seq_data_align()
#         # thread_exit.append(False)
#         # thread.start_new_thread(stat_thread, (process_pip, thread_init))
#         # thread_init += 1
#         # while sum([True if sg is False else False for sg in thread_exit]) >= threading_max:
#         #     sleep(10)
#
# while False in thread_exit:
#     sleep(5)
#
# for sample_rep in sample_names_rep:
#     # sample_fa = sample['fasta_file']
#     # sample_name = sample['name']
#     # sample_replicates = list(set([fafile.strip(fafile.split('_')[-1])[0:-1] for fafile in sample_fa]))
#     # for rep_index, replicate in enumerate(sample_replicates):
#     cmd_cp = f"cp {os.path.join(sampleinfo['save_dir'], sample_rep + '_output', sample_rep + '.expression_statistic.csv')} " + \
#              f"{os.path.join(sampleinfo['save_dir'], sample_rep + '.expression_statistic.csv')}"
#     sbps.run(cmd_cp, shell=True)

# %%
