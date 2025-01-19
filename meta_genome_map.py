#%%
from RNA_seq_analyzer import DNASeqAnalyzer
import os
from seq_utility import BAMile
import os
import sys
import pysam
from tqdm import tqdm
import numpy as np
# […]
import numpy as np
from scipy.stats import binned_statistic


# %%  # set the path
fastq_file_dir = r'/media/fulab/fulab-nas/chupan/fulab_zc_1/seq_data/ZhenHai/mouse_gut_metagenomes/fastq_files'
reference_file_path = r'/media/fulab/fulab-nas/chupan/fulab_zc_1/seq_data/ZhenHai/mouse_gut_metagenomes/Alistipes_shahii_WAL_8301_complete_sequence.fasta'
# fastq_file_dir = r'Y:\chupan\fulab_zc_1\seq_data\ZhenHai\mouse gut metagenomes\fastq_files'
# reference_file_path = r'Y:\chupan\fulab_zc_1\seq_data\ZhenHai\mouse gut metagenomes\Alistipes shahii WAL 8301 complete sequence.gb'
fastq_files_list = [file for file in os.listdir(fastq_file_dir) if file.endswith('.fastq.gz') or file.endswith('.fq.gz')]

for fa_file in fastq_files_list:
    name_space = fa_file.split('.')
    if 'pair' in name_space:
        if name_space[-3] == '2':
            fastq_files_list.remove(fa_file)

fastq_files_path = [os.path.join(fastq_file_dir, file) for file in fastq_files_list]

# %%  # Perform alignment
pars = {'meta_low': {"-p": 60, "score-min": "G,1,0.1", 'match_bonus':'1', 'mismatch_penalty_max-min': '1,1'},
        'meta_high': {"-p": 60, "score-min": "G,9,10", }}
# try files one by on
for file in fastq_files_path:
    sample = DNASeqAnalyzer(sample_name=os.path.basename(file), seq_file_path=file,
                            ref_ps=reference_file_path, output_dir='./', bowtie_pars={"-p": 60, "score-min": "G,9,10", },)
    sample.seq_data_align()

    bam_file = BAMile(sample.bam_sorted_ps, gff_ps=None, reference_ps=sample.reference_file_path,
                      paired_flag=False)
    bam_file.calculate_depth()


#%% Statistics
# -*- coding: utf-8 -*-

"""
This code is used to clean the rew sequencing data.
 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""
#%%
# Built-in/Generic Imports

# Own modules
align_dir = r'/media/fulab/fulab-nas/chupan/fulab_zc_1/seq_data/ZhenHai/precise_align'

all_dir = os.listdir(align_dir)
#%%
# find all bam files
bam_files = []
for _dir in all_dir:
    if os.path.isdir(os.path.join(align_dir, _dir)):
        bam_files_ = [file for file in os.listdir(os.path.join(align_dir, _dir)) if file.endswith('.bam')]
        for bam_file in bam_files_:
            suffix = 'sorted.bam'
            if suffix in bam_file[-len(suffix):]:
                bam_files.append(os.path.join(align_dir, _dir, bam_file))

# merge them to one bam file
bam_files_str = ' '.join(bam_files)
cmd = f'samtools merge -f {align_dir}/all_samples.bam {bam_files_str}'
os.system(cmd)

#%%
# statistics of the bam file

# .a calculate the average quality of the alignment along the genome
sam_file = pysam.AlignmentFile(f'{align_dir}/all_samples.bam', 'rb')
index_name = [index[0] for index in sam_file.get_index_statistics()][0]


# calculate the coverage of the genome
coverage = sam_file.count_coverage(index_name)
total_coverage = np.sum(coverage, axis=0)

_itr = sam_file.pileup(index_name, 1)
average_mapping_quality_list = []
mapping_number_list = []
mapping_quality_list = []
for _i in _itr:
    # number_align = _i.get_num_aligned()
    # align_pos = _i.get_query_positions()
    mapping_quality = _i.get_mapping_qualities()
    mapping_number_list.append(len(mapping_quality))
    mapping_quality_list.append(mapping_quality)
    average_align_quality = sum(mapping_quality) / len(mapping_quality)
    average_mapping_quality_list.append(average_align_quality)

#%%


# .a sum of reads in each bin
bin_length = 15200  # 29024 for 125 bp; 15200 for 250 bp; 15052 for shift 125 bp
binned_mapping_number = binned_statistic(np.arange(1, len(total_coverage)+1), total_coverage,
                                         bins=bin_length, statistic='sum')
location_axis = (binned_mapping_number.bin_edges[1:] + binned_mapping_number.bin_edges[:-1]) / 2
binned_mapping_number_array = np.vstack((location_axis, binned_mapping_number.statistic)).T

# .b calculate the average quality of the alignment along the genome
average_mapping_quality_list = []
for region_i in tqdm(range(len(binned_mapping_number.bin_edges)-1)):
    # fetch reads mapped to the region
    reads_iter = sam_file.fetch(index_name,
                                int(binned_mapping_number.bin_edges[region_i]),
                                int(binned_mapping_number.bin_edges[region_i+1]))
    reads_qual_mean = np.mean([read.mapping_quality for read in reads_iter])
    average_mapping_quality_list.append(reads_qual_mean)
average_mapping_quality_array = np.vstack((location_axis, average_mapping_quality_list)).T



# .b calculate the number of high quality alignment along the genome
high_quality_mapping_number_list = []
high_quality_mapping_threshold = 30
for mapping_quality in tqdm(mapping_quality_list):
    high_quality_mapping_number = len([quality for quality in mapping_quality if quality >= high_quality_mapping_threshold])
    high_quality_mapping_number_list.append(high_quality_mapping_number)

bin_length = 4535
binned_high_quality_mapping = binned_statistic(np.arange(len(high_quality_mapping_number_list)), high_quality_mapping_number_list, bins=bin_length)
location_axis = (binned_high_quality_mapping.bin_edges[1:] + binned_high_quality_mapping.bin_edges[:-1]) / 2
binned_high_quality_mapping_array = np.vstack((location_axis, binned_high_quality_mapping.statistic)).T

# .c calculate the number of low quality alignment along the genome
low_quality_mapping_number_list = []
low_quality_mapping_threshold = 20
for mapping_quality in tqdm(mapping_quality_list):
    low_quality_mapping_number = len([quality for quality in mapping_quality if quality <= low_quality_mapping_threshold])
    low_quality_mapping_number_list.append(low_quality_mapping_number)

bin_length = 800
binned_low_quality_mapping = binned_statistic(np.arange(len(low_quality_mapping_number_list)), low_quality_mapping_number_list, bins=bin_length)
location_axis = (binned_low_quality_mapping.bin_edges[1:] + binned_low_quality_mapping.bin_edges[:-1]) / 2
binned_low_quality_mapping_array = np.vstack((location_axis, binned_low_quality_mapping.statistic)).T

# .c binned average of the alignment quality

bin_length = 4535
average_mapping_quality_list = np.array(average_mapping_quality_list)
binned_quality = binned_statistic(np.arange(len(average_mapping_quality_list)), average_mapping_quality_list, bins=bin_length)
location_axis = (binned_quality.bin_edges[1:] + binned_quality.bin_edges[:-1]) / 2
binned_quality_array = np.vstack((location_axis, binned_quality.statistic)).T

# .d high quality reads ratio of the alignment quality
bin_length = 4535
high_quality_ratio = np.array(high_quality_mapping_number_list) / np.array(mapping_number_list)
binned_high_quality_ratio = binned_statistic(np.arange(len(high_quality_ratio)), high_quality_ratio, bins=bin_length)
location_axis = (binned_high_quality_ratio.bin_edges[1:] + binned_high_quality_ratio.bin_edges[:-1]) / 2
binned_high_quality_ratio_array = np.vstack((location_axis, binned_high_quality_ratio.statistic)).T

high_quality_region_index = np.where(binned_high_quality_ratio_array[:, 1] > 0.9)[0]
high_quality_region_loaction = binned_high_quality_ratio_array[high_quality_region_index, 0]
for loc in high_quality_region_loaction:
    print(f'{int(loc)}')

# .e binned mapping number
bin_length = 4535
binned_mapping_number = binned_statistic(np.arange(len(mapping_number_list)), mapping_number_list, bins=bin_length)
location_axis = (binned_mapping_number.bin_edges[1:] + binned_mapping_number.bin_edges[:-1]) / 2
binned_mapping_number_array = np.vstack((location_axis, binned_mapping_number.statistic)).T
