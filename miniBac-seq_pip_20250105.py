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
from seq_utility import find_fq
from RNA_seq_analyzer import RNASeqAnalyzer
# import _thread as thread
from threading import Thread

def stat_thread(obj: RNASeqAnalyzer, thread_index: int, thread_exit):
    obj.counts_statistic()
    thread_exit[thread_index] = True
    return None




# %%
# gff_ps = r'./annotation_file/L3_strain/L3_strain.gff'
# fasta_ps = r'./annotation_file/L3_strain/L3_strain.fa'
gff_ps = r'./annotation_file/xcd001_reference/xcd001.1.gff'
fasta_ps = r'./annotation_file/xcd001_reference/xcd001.1.fa'
fastq_dir = r'/media/fulab/fulab-nas/chupan/fulab_zc_1/seq_data/20241201_RNA-seq/E2_data'
output_dir = r'/media/fulab/fulab-nas/chupan/fulab_zc_1/seq_data/20241201_RNA-seq/E2_data'
threading_max = 3


thread_init = 0
sample_names_rep = []
fastq_files = [file for file in os.listdir(fastq_dir)
               if os.path.isfile(os.path.join(fastq_dir, file))]


samples_dict = find_fq(fastq_dir)
samples = list(samples_dict.keys())
samples.sort()

print('''all samples' name: ''', samples)

thread_exit = []
statistic_list = []
for sample in samples:
    read1 = samples_dict[sample]['R1']
    read2 = samples_dict[sample]['R2']
    process_pip = RNASeqAnalyzer(sample, seq_ps1=read1, seq_ps2=read2, bowtie_pars={"-p": 60}, output_dir=output_dir,
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

