# -*- coding: utf-8 -*-

"""

 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""

# Built-in/Generic Imports
import os
import subprocess as sbps
import re
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


# ANSI color codes
RED = "\033[31m"
GREEN = "\033[32m"
YELLOW = "\033[33m"
BLUE = "\033[34m"
RESET = "\033[0m"  # Resets color to default


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
    # files form different sequencing company may have different naming rules, so I need to check the file name
    # some use _ to identify the paired files, some use ., use regular expression to find the paired files
    # file names compose sample_name, paired number, and suffix
    # paired number can be: .1/.2; R1/R2; _1/_2; _R1/_R2; _1_/_2_; _R1_/_R2_; 1/2; R1/R2;
    for sufix in suffix:
        redas_files += [file for file in files if file[-len(sufix):] == sufix]
    # check the file name
    samples_dict = {}
    for file in redas_files:
        # check the file name
        match = re.match(r'(.+?)([._][Rr]?[12])(\..+)?$', file)
        if match:
            sample_name, paired_num, suffix = match.groups()
            if sample_name not in samples_dict:
                samples_dict[sample_name] = {}
                samples_dict[sample_name][paired_num] = (sample_name, paired_num, suffix)
            else:
                samples_dict[sample_name][paired_num] = (sample_name, paired_num, suffix)
    # check the paired number
    for sample_name, paired_files in samples_dict.items():
        if len(paired_files) == 1:
            file_name = list(paired_files.values())[0]
            samples_dict[sample_name]['R1'] = os.path.join(dir_name, ''.join(file_name))
            samples_dict[sample_name]['R2'] = None
            samples_dict[sample_name]['paired'] = False
        elif len(paired_files) == 2:
            # whatever the samples file how to name their sequence files, I identify the file type by the number
            file_keys = list(paired_files.keys())
            if '1' in file_keys[0]:
                samples_dict[sample_name]['R1'] = os.path.join(dir_name, ''.join(paired_files[file_keys[0]]))
                samples_dict[sample_name]['R2'] = os.path.join(dir_name, ''.join(paired_files[file_keys[1]]))
            else:
                samples_dict[sample_name]['R1'] = os.path.join(dir_name, ''.join(paired_files[file_keys[1]]))
                samples_dict[sample_name]['R2'] = os.path.join(dir_name, ''.join(paired_files[file_keys[0]]))

            samples_dict[sample_name]['paired'] = True


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

cpu_num = 8

# raw data folder
raw_data_folder = r'/media/fulab/fulab-nas/chupan/fulab_zc_1/seq_data/20251222_RNA-seq/RawData/N2601923_80-2084156168_20260125134435/260123-A00199B'
savdir = '/media/fulab/fulab-nas/chupan/fulab_zc_1/seq_data/20251222_RNA-seq/cleaned_data'
deduplication = False
seq_file_suffix = ['.fastq.gz', '.fq.gz', '.fastq', '.fq']
# pause the program and check the folder, press Enter for continue
while True:
    input_folder = input(f'{GREEN}Please check the folder (press Enter to continue or specify dir): \n{RESET}{raw_data_folder}\n')
    if input_folder == '':
        break
    else:
        # check input is a folder
        if os.path.isdir(input_folder):
            raw_data_folder = input_folder
            print(f'{GREEN}Folder is set to: {RESET}{raw_data_folder}')
            break
        else:
            print(f'{GREEN}Please press Enter to continue.{RESET}')

# check the save folder
# print(f'Please check the save folder (press Enter to continue): \n{savdir}')
while True:
    input_savdir = input(f'{GREEN}Please check the folder (press Enter to continue or specify dir):{RESET}\n {savdir}')

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
                    print(f'Folder is set to: {savdir}')
                    break
                elif create_folder == 'n':
                    print('The folder is not created, exit.')
                    exit(0)
                else:
                    pass


# search all files in the folder
folders = os.listdir(raw_data_folder)
folders = [folder for folder in folders
           if os.path.isdir(os.path.join(raw_data_folder, folder))]
# sample in main folder
samples = find_fq(raw_data_folder, seq_file_suffix)
# if have subfolder, find all samples
if len(folders) > 0:
    for folder in folders:
        samples.update(find_fq(os.path.join(raw_data_folder, folder), seq_file_suffix))

# crate trimmed data folder
if not os.path.exists(savdir):
    os.makedirs(savdir)

# check all files found
if len(samples) == 0:
    print(f'{RED}No fastq files found in {raw_data_folder}, please check the folder.{RESET}')
    exit(0)
else:
    print('''all samples' name: ''', samples)
    while True:
        for sample, seq_files in samples.items():
            print(f"{YELLOW}Sample: {sample}{RESET}")
            print(f"-> R1: {seq_files['R1']}")
            if seq_files['paired']:
                print(f"-> R2: {seq_files['R2']}")
            else:
                print("-> R2: None")
            print(f"-> Paired: {seq_files['paired']}")
        check_samples = input('Please check the samples (press Enter to continue or specify sample name, seperate by comma): ')
        if check_samples == '':
            break
        else:
            check_samples = check_samples.split(',')
            for sample in check_samples:
                if sample not in samples.keys():
                    print(f'Sample {sample} not found, please check the folder.')
                    break
            else:
                break

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
    if deduplication:
        fastp_command += ' --dedup'
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

