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
    # stat = sbps.Popen('source ~/.bashrc && conda activate bioinfo && '+cmd,
    #                   shell=True, cwd=cwd)
    stat = sbps.run(cmd, shell=True, cwd=cwd)

    return stat

def find_fq(dir_name, suffix=None):
    """
    Find the fastq files in the folder.
    :param dir_name: str
        The folder name.
    :return: list
        The fastq files in the folder.
    """
    if suffix is None:
        suffix = ['.fastq.gz', '.fq.gz', '.fastq', '.fq']
    files = [file.name for file in os.scandir(dir_name) if file.is_file()]
    redas_files = []
    for sufix in suffix:
        redas_files += [file for file in files if file[-len(sufix):] == sufix]
    samples = list(set([file.split('.')[0] for file in redas_files]))
    samples_dict = {}
    for sample in samples:
        files_of_sample = [seqfile for seqfile in redas_files
                           if seqfile.split('.')[0] == sample]
        sample_dic = {}
        if len(files_of_sample) == 1:
            sample_dic['R1'] = os.path.join(dir_name, files_of_sample[0])
            sample_dic['R2'] = None
            sample_dic['paired'] = False
        elif len(files_of_sample) == 2:
            if '1' in files_of_sample[0].strip('sample'):
                sample_dic['R1'] = os.path.join(dir_name, files_of_sample[0])
                sample_dic['R2'] = os.path.join(dir_name, files_of_sample[1])
            else:
                sample_dic['R1'] = os.path.join(dir_name, files_of_sample[1])
                sample_dic['R2'] = os.path.join(dir_name, files_of_sample[0])

            sample_dic['paired'] = True
        samples_dict[sample] = sample_dic

    return samples_dict

# Own modules
#%%

"""
folder structure
--------- work folder -----------------
   |_____Raw data folder
      |_____ folders containing *.fastq.gz
   |_____Cleaned data folder
      |_____ *.fastq.gz
"""

cpu_num = 6

# raw data folder
raw_data_folder = r'/media/fulab/fulab-nas/chupan/fulab_zc_1/seq_data/20241201_RNA-seq/rawdata_mixed'
savdir = '/media/fulab/fulab-nas/chupan/fulab_zc_1/seq_data/20241201_RNA-seq/rawdata_mixed'
seq_file_suffix = '.fastq.gz'
# search all files in the folder
folders = os.listdir(raw_data_folder)
folders = [folder for folder in folders
           if os.path.isdir(os.path.join(raw_data_folder, folder))]
# sample in main folder
samples = find_fq(raw_data_folder, [seq_file_suffix])
# if have subfolder, find all samples
if len(folders) > 0:
    for folder in folders:
        samples.update(find_fq(os.path.join(raw_data_folder, folder), [seq_file_suffix]))


# crate trimmed data folder
if not os.path.exists(savdir):
    os.makedirs(savdir)

# clean the data
workers = []
for sample, seq_files in samples.items():
    read1 = seq_files['R1']
    read2 = seq_files['R2']
    if seq_files['paired']:
        output_read1 = os.path.join(savdir, sample + ".R1.fastq.gz")
        output_read2 = os.path.join(savdir, sample + ".R2.fastq.gz")
    else:
        output_read1 = os.path.join(savdir, sample + ".fastq.gz")
        output_read2 = None

    # =================== fastp for QC and cut adapter ================================
    # this command will remove the adapter directly with args: -U --umi_loc=read1 --umi_len=7
    # Attention! Not for demultiplexing the samples.
    fastp_command = (f'fastp -i {read1} -I {read2} -o {output_read1} -O {output_read2} ' +
                     f'-U --umi_loc=read1 --umi_len=7 ' +  # remove the adapter directly
                     f'-h {os.path.join(savdir, sample + "_fastp_report.html")} ' +
                     f'-j {os.path.join(savdir, sample + "_fastp_report.json")} ' + f'-w {cpu_num}')
    print(fastp_command)
    workers.append(Thread(target=run_cmd, args=(fastp_command, )))


# workers = []
# for filer_i, folder in enumerate(folders):
#
#     # folder = folders[0]
#     folder_dir = os.path.join(raw_data_folder, folder)
#     print(f"FASTQ folder scan: {folder_dir}")
#     files = [ file.name for file in os.scandir(folder_dir) if file.is_file()]
#     suffiex = '.fastq.gz'
#     readsfile_name = [file for file in files if file[-len(suffiex):] == suffiex]
#     for file in readsfile_name:
#         if file[:-len(suffiex)].split('_')[-1] == 'R1':
#             read1 = os.path.join(folder_dir,file)
#         elif file[:-len(suffiex)].split('_')[-1] == 'R2':
#             read2 = os.path.join(folder_dir,file)
#     sample_name = samples[filer_i]
#     output_read1 = os.path.join(savdir, sample_name + "_R1.fastq.gz")
#     output_read2 = os.path.join(savdir, sample_name + "_R2.fastq.gz")
#
#     # =================== fastp for QC and cut adapter ================================
#     fastp_command = f'fastp -i {read1} -I {read2} -o {output_read1} -O {output_read2} '
#
#     # 1. cut the sequence of the flank sequence (insertions the smaller than 150 bp).
#     # 2. cut the fixed pooling barcode in read 1.
#
#     # =================== cut adapter specified GTGTGAA ================================
#     # cutadapt_command = f'cutadapt -j {cpu_num} ' \
#     #                     + '-a CTGTCTCTTATACACATCTCCGAGCCCACGAGACTCGCCTTAATCTCGTATGCCGTCTTCTGCTTGX ' \
#     #                     + '-A TTCACACAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGGCTATAGTGTAGATCTCGGTGGTCGCCGTATCATTX ' \
#     #                     + f'-q 20,20 --minimum-length 100:100 --max-n 3 ' \
#     #                     + '--interleaved ' \
#     #                     + f'--json trim_r1_report_{sample_name}.json ' \
#     #                     + f'{read1} {read2} | ' \
#     #                     + f'cutadapt -j {cpu_num} ' \
#     #                     + '-g "^GTGTGAA" --interleaved ' \
#     #                     + f'--json trim_r2_report_{sample_name}.json ' \
#     #                     + f'-o {output_read1} -p {output_read2}  -' \
#     # =================== cut fixed length ================================
#     log_r1 = os.path.join(savdir, 'trimmed_log', f'trim_r1_report_{sample_name}.json')
#     log_r2 = os.path.join(savdir, 'trimmed_log', f'trim_r2_report_{sample_name}.json')
#     log_folder = os.path.join(savdir, 'trimmed_log')
#     if not os.path.exists(log_folder):
#         os.makedirs(log_folder)
#     cutadapt_command = f'cutadapt -j {cpu_num} ' \
#                         + '-a CTGTCTCTTATACACATCTCCGAGCCCACGAGACTCGCCTTAATCTCGTATGCCGTCTTCTGCTTGX ' \
#                         + '-A TTCACACAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGGCTATAGTGTAGATCTCGGTGGTCGCCGTATCATTX ' \
#                         + f'-q 20,20 --minimum-length 100:100 --max-n 3 --pair-filter=any ' \
#                         + '--interleaved ' \
#                         + f'--json {log_r1} ' \
#                         + f'{read1} {read2} | ' \
#                         + f'cutadapt -j {cpu_num} ' \
#                         + '-u 7 --interleaved ' \
#                         + f'--json {log_r2} ' \
#                         + f'-o {output_read1} -p {output_read2}  -' \
#
#
#     print(cutadapt_command)
#
#     workers.append(Thread(target=run_cmd, args=(cutadapt_command, )))

max_active_num = int(64/cpu_num)
for worker in workers:
    while active_count() >= max_active_num:
        sleep(5)
    worker.start()
    worker.join()




print('All sequences are trimmed.')

