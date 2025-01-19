# -*- coding: utf-8 -*-

"""

 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""

# Built-in/Generic Imports
import os
import sys
# […]

# Libs
import pandas as pd
import numpy as np  # Or any other
import subprocess as sbps
import datetime
from typing import Union, Optional
from seq_utility import count_feature_reads, BAMile


# […]

# Own modules


def file_prefix(prefix, ps):
    return prefix + os.path.basename(ps)


class DNASeqAnalyzer:
    def __init__(self, sample_name: str, ref_ps: str, gff_ps: str = None,
                 seq_file_path: Union[str, list] = None,
                 bowtie_pars: dict = None, output_dir: Optional[str] = None):
        self.sample_name = sample_name
        
        # high throughput sequencing data related.
        self.seq_file_path = seq_file_path
        # reference file related.
        self.reference_file_path = ref_ps
        self.reference_file_name = os.path.basename(self.reference_file_path)  # type: str
        self.reference_file_dir = os.path.dirname(self.reference_file_path)
        self.indexed_base_name = '.'.join(self.reference_file_name.split('.'))
        self.gff_ps = gff_ps
        # bowtie2 parameters
        self.bowtie_pars = {'-N': 1, '-q': '--phred64', '-p': 8, 'score-min': 'G,9,8'}
        if bowtie_pars is not None:
            for key in list(bowtie_pars.keys()):
                self.bowtie_pars[key] = bowtie_pars[key]
        # init output_dir
        if output_dir is None:
            self.output_dir = os.path.join(os.getcwd(), f'{self.sample_name}_output')
        else:
            self.output_dir = os.path.join(output_dir, f'{self.sample_name}_output')
        try:
            os.makedirs(self.output_dir)
        except FileExistsError:
            print(f'[{self.sample_name}] -> Attention! Dir {self.output_dir} already existed!')

        self.file_in_dir = os.listdir(self.output_dir)
        self.log_file_ps = os.path.join(self.output_dir, self.sample_name + '.log')
        # out put files
        self.sam_file_ps = os.path.join(self.output_dir, self.sample_name + '.sam')
        self.bam_ps = os.path.join(self.output_dir, self.sample_name + '.bam')
        self.bam_sorted_ps = os.path.join(self.output_dir, self.sample_name + '.sorted.bam')
        self.bam_index_ps = self.bam_sorted_ps + '.bai'
    def seq_data_align(self):
        """ Mapping the sequence data to the reference genome. 
        """

        # Step 1. make index files

        if os.path.basename(self.reference_file_name) not in self.file_in_dir:
            cmd_copy_ref = f'cp {self.reference_file_path} ' \
                           f'{os.path.join(self.output_dir, self.reference_file_name)}'
            # update the reference file ps
            self.reference_file_path = os.path.join(self.output_dir, self.reference_file_name)
            self.append_to_log(f'[{self.sample_name}] -> Copy Reference: {cmd_copy_ref}')
            status1 = self.cmd_shell(cmd_copy_ref)
            cmd_index = f'bowtie2-build -f {self.reference_file_name} {self.indexed_base_name}'
            print(f'[{self.sample_name}] -> Generate indexed reference: {cmd_index}')
            status2 = self.cmd_shell(cmd_index, cwd=self.output_dir)
        # Step 2. mapping reads, Currently only support single end reads.
        if os.path.exists(self.sam_file_ps):
            print(f'[{self.sample_name}] -> Skip mapping reads.')
            self.append_to_log(f'[{self.sample_name}] -> Skip mapping reads.')
        else:
            os.environ['BOWTIE2_INDEXES'] = self.output_dir
            seq_file_path_string = self.seq_file_path if isinstance(self.seq_file_path, str) else ','.join(self.seq_file_path)
            cmd_align = f'bowtie2 -p {self.bowtie_pars["-p"]} ' + \
                f'--local --no-unal ' \
                f'-N {self.bowtie_pars["-N"]} -x {self.indexed_base_name} ' \
                f' -U {seq_file_path_string} ' \
                f'-S {self.sam_file_ps}'
            # match bonus, int default 2
            if self.bowtie_pars.get('match_bonus') is not None:
                cmd_align = cmd_align + f' --ma {self.bowtie_pars["match_bonus"]}'
            # Sets the maximum (MX) and minimum (MN) mismatch penalties.  MX = 6, MN = 2
            if self.bowtie_pars.get('mismatch_penalty_max-min') is not None:
                cmd_align = cmd_align + f' --mp {self.bowtie_pars["mismatch_penalty_max-min"]}'
            # score-min, default G,20,8
            if self.bowtie_pars.get('score-min') is not None:
                cmd_align = cmd_align + f' --score-min {self.bowtie_pars["score-min"]}'

            print(f"[{self.sample_name}] -> Mapping reads: " + cmd_align)
            self.append_to_log(f"[{self.sample_name}] -> Mapping reads: " + cmd_align)
            status3 = self.cmd_shell(cmd_align)
            del os.environ['BOWTIE2_INDEXES']
        # Step 3. generate bam file
        # sort alignments
        if os.path.exists(self.bam_sorted_ps):
            print(f'[{self.sample_name}] -> Skip generating BAM.')
            self.append_to_log(f'[{self.sample_name}] -> Skip generating BAM.')
        else:
            cmd_gen_bam = f'samtools view -bS -@ {self.bowtie_pars["-p"]} {self.sam_file_ps} | samtools sort -o {self.bam_sorted_ps} -@ {self.bowtie_pars["-p"]}'
            print(f'[{self.sample_name}] -> Generating BAM: {cmd_gen_bam}')
            self.append_to_log(f'[{self.sample_name}] -> Generating BAM: {cmd_gen_bam}')
            status4 = self.cmd_shell(cmd_gen_bam)
            cmd_index_bam = f'samtools index -@ {self.bowtie_pars["-p"]} {self.bam_sorted_ps} {self.bam_index_ps}'
            print(f'[{self.sample_name}] -> Generating .bai: {cmd_index_bam}')
            status5 = self.cmd_shell(cmd_index_bam)
    def cmd_shell(self, cmd: str, cwd=None):
        """"
        execute command in shell and append output into log file.
        """
        if cwd is None:
            cwd = os.getcwd()
        stat = sbps.Popen(cmd, shell=True, stdout=sbps.PIPE, stderr=sbps.PIPE, cwd=cwd)
        std_out_log = stat.stdout.read()
        std_err_log = stat.stderr.read()

        self.append_to_log(std_out_log)  # write log
        if len(std_out_log) > 0:
            print(f'[{self.sample_name}] -> ' )
            print(std_out_log.decode('ascii'))

        self.append_to_log(std_err_log)  # write log
        if len(std_err_log) > 0:
            print(f'[{self.sample_name}] -> {std_err_log}' )
            print(std_err_log.decode('ascii'))

        return stat
    def append_to_log(self, stdout: Union[str, bytes]):
        with open(self.log_file_ps, 'a') as self.log_file:
            try:
                self.log_file.write(stdout.decode("utf-8") + '\n')
            except AttributeError:
                self.log_file.write(stdout + '\n')

class RNASeqAnalyzer:
    def __init__(self, sample_name: str, ref_ps: str, gff_ps: str = None,
                 seq_ps1: str = None, seq_ps2: str = None, adapter: list = None,
                 bowtie_pars: dict = None,
                 output_dir: str = None):
        self.sample_name = sample_name  # type: str # sample name
        self.reference_file_path = ref_ps  # type: str # genome fasta file
        self.gff_ps = gff_ps
        self.raw_seq_data_ps1 = seq_ps1
        self.raw_seq_data_ps2 = seq_ps2
        if self.raw_seq_data_ps2 is None:
            self.paired_flag = False
        else:
            self.paired_flag = True
        self.reference_file_name = os.path.basename(self.reference_file_path)  # type: str
        self.reference_file_dir = os.path.dirname(self.reference_file_path)
        self.indexed_base_name = '.'.join(self.reference_file_name.split('.')[:-1])
        self.adapter = adapter  # type: Union[str, None]  # adapter sequence for cutting.
        self.bowtie_pars = {'-N': 1, '-q': '--phred64', '-p': 8}
        self.time_now = datetime.datetime.now()
        if bowtie_pars is not None:
            for key in list(bowtie_pars.keys()):
                self.bowtie_pars[key] = bowtie_pars[key]
        # init output_dir
        if output_dir is None:
            self.output_dir = os.path.join(os.getcwd(), f'{self.sample_name}_output')
        else:
            self.output_dir = os.path.join(output_dir, f'{self.sample_name}_output')

        try:
            os.makedirs(self.output_dir)
        except FileExistsError:
            print(f'[{self.sample_name}] -> Attention! Dir {self.output_dir} already existed!')

        if self.adapter is None:
            self.clean_adapter = False
            self.seq_data_ps1 = self.raw_seq_data_ps1
            self.seq_data_ps2 = self.raw_seq_data_ps2
        else:
            self.clean_adapter = True
            self.seq_data_ps1 = os.path.join(self.output_dir, file_prefix('trimmed_adapter_', self.raw_seq_data_ps1))
            self.seq_data_ps2 = os.path.join(self.output_dir, file_prefix('trimmed_adapter_', self.raw_seq_data_ps2))
            self.trim_log_ps = os.path.join(self.output_dir, 'trim_log.json')

        self.sam_file_ps = os.path.join(self.output_dir, self.sample_name + '.sam')
        self.bam_ps = os.path.join(self.output_dir, self.sample_name + '.bam')
        self.bam_sorted_ps = os.path.join(self.output_dir, self.sample_name + '.sorted.bam')
        self.bam_index_ps = self.bam_sorted_ps + '.bai'
        self.counts_statistic_ps = os.path.join(self.output_dir, self.sample_name + '.expression_statistic.csv')
        self.file_in_dir = os.listdir(self.output_dir)
        self.log_file_ps = os.path.join(self.output_dir, self.sample_name + '.log')
        if not os.path.exists(self.log_file_ps):
            self.log_file = open(self.log_file_ps, 'w')
        else:
            self.log_file = open(self.log_file_ps, 'a')

        self.log_file.write(f'============= {self.time_now.year}-{self.time_now.month}-{self.time_now.hour}' + \
                            f'-{self.time_now.minute} ==============\n')
        self.log_file.close()

    def append_to_log(self, stdout: Union[str, bytes]):
        with open(self.log_file_ps, 'a') as self.log_file:
            try:
                self.log_file.write(stdout.decode("utf-8") + '\n')
            except AttributeError:
                self.log_file.write(stdout + '\n')

    def seq_data_align(self):

        # clean raw data
        if os.path.basename(self.seq_data_ps1) not in self.file_in_dir:
            if self.clean_adapter:
                # cmd_trime = f'cutadapt -a {self.adapter[0]} -A {self.adapter[1]}' + ' --discard-untrimmed -m 10 -o ' + \
                #             self.seq_data_ps1 + ' -p ' + self.seq_data_ps2 + ' ' + self.raw_seq_data_ps1 + ' ' + self.raw_seq_data_ps2
                cmd_trime = f"cutadapt -g X{self.adapter[0]} -m 150 -o " \
                            f"{self.seq_data_ps1} -p {self.seq_data_ps2}" \
                            f"--json {self.trim_log_ps}" \
                            f" {self.raw_seq_data_ps1} {self.raw_seq_data_ps2}"
                print(f"[{self.sample_name}] -> Removing linker: " + cmd_trime)
                status0 = sbps.Popen(cmd_trime, shell=True, stdout=sbps.PIPE, cwd=os.getcwd())
                self.append_to_log(status0.stdout.read())
            else:
                print(f'[{self.sample_name}] -> Pass linker removing.')
                self.append_to_log(f'[{self.sample_name}] -> Pass linker removing.')

        # make index files
        if os.path.basename(self.reference_file_name) not in self.file_in_dir:
            cmd_copy_ref = f'cp {self.reference_file_path} ' \
                           f'{os.path.join(self.output_dir, self.reference_file_name)}'
            # update the reference file ps
            self.reference_file_path = os.path.join(self.output_dir, self.reference_file_name)
            self.append_to_log(f'[{self.sample_name}] -> Copy Reference: {cmd_copy_ref}')
            status1 = self.cmd_shell(cmd_copy_ref)
            cmd_index = f'bowtie2-build -f {self.reference_file_name} {self.indexed_base_name}'
            print(f'[{self.sample_name}] -> Generate indexed reference: {cmd_index}')
            status2 = self.cmd_shell(cmd_index, cwd=self.output_dir)

        # mapping
        if os.path.basename(self.bam_index_ps) not in self.file_in_dir:
            os.environ['BOWTIE2_INDEXES'] = self.output_dir
            if self.seq_data_ps2 is None:  # unpaired reads
                cmd_align = f'bowtie2 -p {self.bowtie_pars["-p"]} --un-gz {self.output_dir} ' + \
                            f'-N {self.bowtie_pars["-N"]} -x {self.indexed_base_name}' \
                            f' -U {self.seq_data_ps1} ' \
                            f'-S {self.sam_file_ps}'
            else:  # paired reads
                cmd_align = f'bowtie2 -p {self.bowtie_pars["-p"]} --un-gz {self.output_dir} ' + \
                            f'--very-sensitive-local -X 1000 -I 18 --no-1mm-upfront --score-min G,9,8 --no-mixed --no-discordant ' \
                            f'-N {self.bowtie_pars["-N"]} -x {self.indexed_base_name} ' \
                            f' -1 {self.seq_data_ps1} -2 {self.seq_data_ps2} ' \
                            f'-S {self.sam_file_ps}'
                """
                change log
                20240113 add new parameters: 
                    --very-sensitive-local -X 1000 -I 18 --no-1mm-upfront 
                    --score-min G,9,8 --no-mixed --no-discordant 
                """

            print(f"[{self.sample_name}] -> Mapping reads: " + cmd_align)
            self.append_to_log(f"[{self.sample_name}] -> Mapping reads: " + cmd_align)
            status3 = self.cmd_shell(cmd_align)
            del os.environ['BOWTIE2_INDEXES']

            # sort alignments
            cmd_gen_bam = f'samtools view -bS -@ {self.bowtie_pars["-p"]} {self.sam_file_ps} | samtools sort -o {self.bam_sorted_ps} -@ {self.bowtie_pars["-p"]}'
            print(f'[{self.sample_name}] -> Generating BAM: {cmd_gen_bam}')
            self.append_to_log(f'[{self.sample_name}] -> Generating BAM: {cmd_gen_bam}')
            status4 = self.cmd_shell(cmd_gen_bam)
            cmd_index_bam = f'samtools index -@ {self.bowtie_pars["-p"]} {self.bam_sorted_ps} {self.bam_index_ps}'
            print(f'[{self.sample_name}] -> Generating .bai: {cmd_index_bam}')
            status5 = self.cmd_shell(cmd_index_bam)

    def counts_statistic(self):
        print(f'[{self.sample_name}] -> Calculate gene expression level.')
        counts_stat, bam = count_feature_reads(self.bam_sorted_ps, self.gff_ps, self.reference_file_path,
                                               paired_flag=self.paired_flag)
        self.__dict__['counts_stat'] = counts_stat
        self.__dict__['gene_dict'] = bam.gene_features
        self.__dict__['bam'] = bam
        counts_stat.to_csv(self.counts_statistic_ps)

    def cmd_shell(self, cmd: str, cwd=None):
        """"
        execute command in shell and append output into log file.
        """
        if cwd is None:
            cwd = os.getcwd()
        stat = sbps.Popen(cmd, shell=True, stdout=sbps.PIPE, stderr=sbps.PIPE, cwd=cwd)
        self.append_to_log(stat.stdout.read())  # write log
        self.append_to_log(stat.stderr.read())  # write log
        return stat


# %%
if __name__ == '__main__':
    #%%
    read1 = r'./example_data/seq_data/A1/A1.raw_1.fastq.gz'
    read2 = r'./example_data/seq_data/A1/A1.raw_2.fastq.gz'
    gff_file = r'./example_data/annotation_file/NC_000913.3_liulab.gff'
    sample_name = 'TestAD'
    ref_ps = r'./example_data/annotation_file/NC_000913.3_liulab.fa'
    # adapters = ['AGATCGGAAGAGC', 'AGATCGGAAGAGC']
    sample = RNASeqAnalyzer(sample_name=sample_name, ref_ps=ref_ps, gff_ps=gff_file,
                            seq_ps1=read1, seq_ps2=read2, bowtie_pars={"-p": 32})
    sample.seq_data_align()
    sample.counts_statistic()

    # %% Deep Seq
    #
    # from scipy.stats import binned_statistic, linregress
    # import matplotlib.pyplot as plt
    # import scipy.stats as stats
    # import sciplot as splt
    #
    # splt.whitegrid()
    # read1 = './example_data/seq_data/A1/A1.raw_1.fastq.gz'
    # read2 = './example_data/seq_data/A1/A1.raw_2.fastq.gz'
    # ref_ps = '/home/fulab/home/fulab/tmp/pycharm_project_757/example_data/annotation_file/GCA_000005845.2_ASM584v2_genomic.fasta'
    # sample_name = 'MG1655'
    # ori_site = 3925859
    # bin_length = 5000
    # sample = RNASeqAnalyzer(sample_name=sample_name, ref_ps=ref_ps, gff_ps=None,
    #                         seq_ps1=read1, seq_ps2=read2, bowtie_pars={"-p": 32})
    # sample.seq_data_align()
    # bam_file = BAMile(sample.bam_sorted_ps, sample.gff_ps, sample.reference_file_path,
    #                   paired_flag=sample.paired_flag)
    # bam_file.separate_bam_by_strand(clean_rtRNA=False)
    # bam_file.count_coverage()
    # coverage = bam_file.fetch_coverage(bam_file.genome_set[0], ori_site, ori_site - 1)
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
    #
    # data_exp = pd.DataFrame(data=dict(Relative_position=relative_pos,
    #                                   genome_position=genome_index,
    #                                   Count=coverage_binned_mean))
    # data_exp.to_csv(os.path.join(sample.output_dir, f'{sample_name}_depth_statistic.csv'))
    #
    # x_fliter = relative_pos > 0
    # inf_filter = ~np.isinf(np.log(coverage_binned_mean))
    # filter = np.logical_and(x_fliter, inf_filter)
    #
    # x_fliter = relative_pos <= 0
    # inf_filter = ~np.isinf(np.log(coverage_binned_mean))
    # filter2 = np.logical_and(x_fliter, inf_filter)
    #
    # filters = [filter, filter2]
    #
    # fig1, ax2 = plt.subplots(1, 1, figsize=(10, 10))
    # ax2.scatter(relative_pos, np.log2(coverage_binned_mean), c='#85C1E9')
    # ax2.plot()
    # results = []
    # for flt in filters:
    #     ret = linregress(relative_pos[flt], np.log2(coverage_binned_mean)[flt])
    #
    #     results.append(ret)
    #
    #     ax2.plot(relative_pos[flt], ret.intercept + ret.slope * relative_pos[flt],
    #              '--r', label='Slope: %.3f' % ret.slope, c='#F1948A')
    #
    # ax2.set_title('Average Slope %.3f' % np.mean([np.abs(ret.slope) for ret in results]))
    # ax2.legend()
    # fig1.show()
    #
    # fig1.savefig(os.path.join(sample.output_dir, f'{sample_name}_depth_statistic.svg'), transparent=True)
