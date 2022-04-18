# -*- coding: utf-8 -*-

"""

 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""
# %%
# Built-in/Generic Imports
import os
import sys
# […]

# Libs
import time

import pandas as pd
import numpy as np  # Or any other
# […]

# Own modules


from scipy.stats import binned_statistic, linregress
import matplotlib.pyplot as plt
import sciplot as splt
from RNA_seq_analyzer import RNASeqAnalyzer
from seq_utility import BAMile
from tqdm import tqdm
import _thread as thread
import subprocess as sbp

splt.whitegrid()

global_lock = thread.allocate_lock()


def deep_seq_pip(sample_name, ref_ps, reads, ori_site, bin_length, export_dir, index) -> int:
    """

    Parameters
    ----------
    sample_name : string
        sample name
    ref_ps : string
        reference genome path
    reads : str or list

    ori_site : int
        location of oriC
    bin_length : int
        binned size
    export_dir : str
        save data dir

    Returns
    -----------
    None

    """
    if reads is str:
        read1, read2 = reads, None
    else:
        read1, read2 = reads
    sample = RNASeqAnalyzer(sample_name=sample_name, ref_ps=ref_ps, gff_ps=None,
                            seq_ps1=read1, seq_ps2=read2, bowtie_pars={"-p": 64}, output_dir=export_dir)
    sample.seq_data_align()

    def coverage_process(sample: RNASeqAnalyzer, ori_site, bin_length, index):
        thread_state[index] = True
        bam_file = BAMile(sample.bam_sorted_ps, sample.gff_ps, sample.reference_file_path,
                          paired_flag=sample.paired_flag)
        bam_file.separate_bam_by_strand(clean_rtRNA=False)
        bam_file.count_coverage()
        coverage = bam_file.fetch_coverage(bam_file.genome_set[0], ori_site, ori_site - 1, move_average=150)

        genome_length = len(bam_file.genomes[bam_file.genome_set[0]])
        coverage_binned = binned_statistic(np.arange(len(coverage)), coverage, 'mean',
                                           bins=int(genome_length / bin_length))

        coverage_binned_mean = coverage_binned.statistic
        zerio_index = round(len(coverage_binned_mean) / 2)
        coverage_binned_mean = np.roll(coverage_binned_mean, round(zerio_index))

        left_pos = np.linspace(-1, 0, num=zerio_index, endpoint=False)
        right_pos = np.linspace(0, 1, num=(len(coverage_binned_mean) - zerio_index), endpoint=True)
        relative_pos = np.concatenate([left_pos, right_pos])

        genome_index = np.arange(1, genome_length)
        genome_index = np.roll(genome_index, genome_length - ori_site)[::bin_length][:-1]
        inf_filter = coverage_binned_mean > 0
        log2_coverage = np.zeros(len(coverage_binned_mean))
        log2_coverage[inf_filter] = np.log2(coverage_binned_mean[inf_filter])

        data_exp = pd.DataFrame(data=dict(Relative_position=relative_pos,
                                          genome_position=genome_index,
                                          Count=coverage_binned_mean,
                                          Log2_count=log2_coverage))
        data_exp.to_csv(os.path.join(sample.output_dir, f'{sample_name}_depth_statistic.csv'))

        x_fliter = relative_pos >= 0
        filter = np.logical_and(x_fliter, inf_filter)

        x_fliter2 = relative_pos <= 0
        filter2 = np.logical_and(x_fliter2, inf_filter)

        filters = [filter, filter2]

        fig1, ax2 = plt.subplots(1, 1, figsize=(12, 12))
        results = []

        for flt in filters:
            ret = linregress(relative_pos[flt], np.log2(coverage_binned_mean[flt]))
            results.append(ret)
            ax2.scatter(relative_pos[flt], np.log2(coverage_binned_mean[flt]), c='#85C1E9')
            ax2.plot(relative_pos[flt], ret.intercept + ret.slope * relative_pos[flt],
                     '--r', label='Slope: %.3f' % ret.slope, c='#F1948A')

        ax2.set_title('%s Average Slope: %.3f' %
                      (sample.sample_name, np.mean([np.abs(ret.slope) for ret in results])),
                      pad=12)
        ax2.set_ylabel('$\mathrm{log}_{2}X_c$', labelpad=7)
        ax2.set_xlabel('$m^{\prime}$', labelpad=7)
        ax2.legend()
        global_lock.acquire()
        fig1.savefig(os.path.join(sample.output_dir, f'{sample_name}_depth_statistic.svg'), transparent=True)
        global_lock.release()
        print(f'[{sample.sample_name}] -> Successful!')
        thread_state[index] = False
        return None

    t_id = thread.start_new_thread(coverage_process, (sample, ori_site, bin_length, index))
    return t_id


# %%
if __name__ == '__main__':
    # %%
    parent_dir = r'/media/fulab/fulab_zc_1/seq_data/LLW_data/20220321_data/soapnuke/clean'
    ref_ps = '/media/fulab/fulab_zc_1/seq_data/Genome_ref/1655_genome_Liu_lab_20220322.fa'
    exp_dir = r'/media/fulab/fulab_zc_1/seq_data/LLW_data/20220321_data/20220321_data_deep_seq_results'
    ori_site = 1  # 3925859
    bin_length = 5000

    sample_dir = [fold.name for fold in os.scandir(parent_dir) if fold.is_dir()]
    sample_msg = {}

    for dir in tqdm(sample_dir, desc=f'[Dir Scanning]'):
        reads = [os.path.join(parent_dir, dir, fa_file.name)
                 for fa_file in os.scandir(os.path.join(parent_dir, dir))
                 if fa_file.name.split('.')[-1] == 'gz']
        sample_msg[dir] = reads

    thread_state = [False] * len(list(sample_msg.keys()))

    for index, (sample, reads) in enumerate(tqdm(sample_msg.items())):
        th_id = deep_seq_pip(sample, ref_ps, reads, ori_site, bin_length, exp_dir, index)

    while True in thread_state:
        time.sleep(5)

    output_dirs = [os.path.join(exp_dir, dir.name)
                   for dir in os.scandir(exp_dir)
                   if dir.name.split('_')[-1] == 'output' and dir.is_dir()]
    rets_file = []
    for dir in output_dirs:
        rets_flies_list = [os.path.join(dir, file.name)
                           for file in os.scandir(dir)
                           if file.name.split('.')[-1] in ['csv', 'svg']]
        rets_file += rets_flies_list

    all_rets_ps = os.path.join(exp_dir, 'all_rests')
    try:
        os.mkdir(all_rets_ps)
    except FileExistsError:
        pass

    cmd = f"cp {' '.join(rets_file)} {all_rets_ps}"

    sbp.run(cmd, shell=True)



#%%
    # # deep_seq_pip('WT', ref_ps, sample_msg['WT'], ori_site, bin_length, exp_dir, 0)
    # sample_name = 'WT'
    # ref_ps = '/media/fulab/Fu_lab_data1/seq_data/20211101_dnaA_datA/1655_genome_Liu_lab.fa'
    # exp_dir = '/media/fulab/Fu_lab_data1/seq_data/20211101_dnaA_datA/dnaAdatA/'
    # sample = RNASeqAnalyzer(sample_name=sample_name, ref_ps=ref_ps, gff_ps=None,
    #                         seq_ps1=sample_msg['WT'][0], seq_ps2=sample_msg['WT'][1], bowtie_pars={"-p": 32},
    #                         output_dir=exp_dir)
    # sample.seq_data_align()
    #
    # bam_file = BAMile(sample.bam_sorted_ps, sample.gff_ps, sample.reference_file_path,
    #                   paired_flag=sample.paired_flag)
    # bam_file.separate_bam_by_strand(clean_rtRNA=False)
    # bam_file.count_coverage()
    # coverage = bam_file.fetch_coverage(bam_file.genome_set[0], ori_site, ori_site - 1, move_average=150)
    #
    # genome_length = len(bam_file.genomes[bam_file.genome_set[0]])
    # coverage_binned = binned_statistic(np.arange(len(coverage)), coverage, 'mean',
    #                                    bins=int(genome_length / bin_length))
    #
    # coverage_binned_mean = coverage_binned.statistic
    # zerio_index = round(len(coverage_binned_mean) / 2)
    # coverage_binned_mean = np.roll(coverage_binned_mean, round(zerio_index))
    #
    # left_pos = np.linspace(-1, 0, num=zerio_index, endpoint=False)
    # right_pos = np.linspace(0, 1, num=(len(coverage_binned_mean) - zerio_index), endpoint=True)
    # relative_pos = np.concatenate([left_pos, right_pos])
    #
    # genome_index = np.arange(1, genome_length)
    # genome_index = np.roll(genome_index, genome_length - ori_site)[::bin_length][:-1]
    # inf_filter = coverage_binned_mean > 0
    # log2_coverage = np.zeros(len(coverage_binned_mean))
    # log2_coverage[inf_filter] = np.log2(coverage_binned_mean[inf_filter])