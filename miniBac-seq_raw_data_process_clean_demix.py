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
    print(f'[Running] -> {cmd}')
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
        # determine the sample name by the 1st part of the file name
        files_of_sample = [seqfile for seqfile in redas_files
                           if seqfile.split('.')[0] == sample]
        # determine the sequence files and types: paired or single.
        sample_dic = {}
        if len(files_of_sample) == 1:
            sample_dic['R1'] = os.path.join(dir_name, files_of_sample[0])
            sample_dic['R2'] = None
            sample_dic['paired'] = False
        elif len(files_of_sample) == 2:
            # whatever the samples file how to name their sequence files, I identify the file type by the number
            # if 1 in the file name, it is R1, otherwise it is R2.
            if '1' in files_of_sample[0].strip(sample):
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
raw_data_folder = r'/media/fulab/fulab-nas/chupan/fulab_zc_1/seq_data/20250303_RNA-seq/rawdata'
savdir = '/media/fulab/fulab-nas/chupan/fulab_zc_1/seq_data/20250303_RNA-seq/cleaned_data'
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
                     f'-U --umi_loc=read1 --umi_len=7 ' +  # remove the adapter directly # ATTENTION!
                     f'-h {os.path.join(savdir, sample + "_fastp_report.html")} ' +
                     f'-j {os.path.join(savdir, sample + "_fastp_report.json")} ' + f'-w {cpu_num}')
    # print(fastp_command)
    workers.append(Thread(target=run_cmd, args=(fastp_command, )))



max_active_num = int(64/cpu_num)
for worker in workers:
    while active_count() >= max_active_num:
        sleep(5)
    worker.start()
for worker in workers:
    worker.join()

print('All sequences are trimmed.')

