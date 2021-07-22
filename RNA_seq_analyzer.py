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
from seq_utility import count_feature_reads

# […]

# Own modules


def file_prefix(prefix, ps):
    return prefix + os.path.basename(ps)


class RNASeqAnalyzer:
    def __init__(self, sample_name: str, ref_ps: str, gff_ps: str = None,
                 seq_ps1: str=None, seq_ps2: str=None, adapter: list = None,
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
        self.indexed_base_name = self.reference_file_name.split('.')[0]
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

        self.sam_file_ps = os.path.join(self.output_dir, self.sample_name + '.sam')
        self.bam_ps = os.path.join(self.output_dir, self.sample_name + '.bam')
        self.bam_sorted_ps = os.path.join(self.output_dir, self.sample_name + '.sorted.bam')
        self.bam_index_ps = self.bam_sorted_ps + '.bai'
        self.counts_statistic_ps = os.path.join(self.output_dir, self.sample_name + '.expression_statistic.csv')
        self.file_in_dir = os.listdir(self.output_dir)
        self.log_file_ps = os.path.join(self.output_dir, self.sample_name + '.log')
        self.log_file = open(self.log_file_ps, 'w')
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
        if os.path.basename(self.seq_data_ps1) not in self.file_in_dir:
            if self.clean_adapter:
                # cmd_trime = f'cutadapt -a {self.adapter[0]} -A {self.adapter[1]}' + ' --discard-untrimmed -m 10 -o ' + \
                #             self.seq_data_ps1 + ' -p ' + self.seq_data_ps2 + ' ' + self.raw_seq_data_ps1 + ' ' + self.raw_seq_data_ps2
                cmd_trime = f"cutadapt -a {self.adapter[0]} -A {self.adapter[1]} --discard-untrimmed -m 10 -o " \
                            f"{self.seq_data_ps1} -p {self.seq_data_ps2}" \
                            f" {self.raw_seq_data_ps1} {self.raw_seq_data_ps2}"
                print(f"[{self.sample_name}] -> Removing linker: " + cmd_trime)
                status0 = sbps.Popen(cmd_trime, shell=True, stdout=sbps.PIPE)
                self.append_to_log(status0.stdout.read())
            else:
                print(f'[{self.sample_name}] -> Pass linker removing.')
                self.append_to_log(f'[{self.sample_name}] -> Pass linker removing.')

        # make index files
        if os.path.basename(self.reference_file_name) not in self.file_in_dir:
            cmd_copy_ref = f'cp {self.reference_file_path} {os.path.join(self.output_dir, self.reference_file_name)}'
            # update the reference file ps
            self.reference_file_path = os.path.join(self.output_dir, self.reference_file_name)
            self.append_to_log(f'[{self.sample_name}] -> Copy Reference: {cmd_copy_ref}')
            status1 = self.cmd_shell(cmd_copy_ref)
            cmd_index = f'cd {self.output_dir}&&bowtie2-build -f {self.reference_file_path} {self.indexed_base_name}'
            print(f'[{self.sample_name}] -> Generate indexed reference: {cmd_index}')
            status2 = self.cmd_shell(cmd_index)

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
                            f'-N {self.bowtie_pars["-N"]} -x {self.indexed_base_name}' \
                            f' -1 {self.seq_data_ps1} -2 {self.seq_data_ps2} '\
                            f'-S {self.sam_file_ps}'
            print(f"[{self.sample_name}] -> Mapping reads: " + cmd_align)
            self.append_to_log(f"[{self.sample_name}] -> Mapping reads: " + cmd_align)
            status3 = self.cmd_shell(cmd_align)
            del os.environ['BOWTIE2_INDEXES']

            # sort alignments
            cmd_gen_bam = f'samtools view -bS -@ 8 {self.sam_file_ps} | samtools sort -o {self.bam_sorted_ps} -@ 8'
            print(f'[{self.sample_name}] -> Generating BAM: {cmd_gen_bam}')
            self.append_to_log(f'[{self.sample_name}] -> Generating BAM: {cmd_gen_bam}')
            status4 = self.cmd_shell(cmd_gen_bam)
            cmd_index_bam = f'samtools index -@ 8 {self.bam_sorted_ps} {self.bam_index_ps}'
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

    def cmd_shell(self, cmd: str):
        """"
        execute command in shell and append output into log file.
        """
        stat = sbps.Popen(cmd, shell=True, stdout=sbps.PIPE)
        self.append_to_log(stat.stdout.read())
        return stat


# %%
if __name__ == '__main__':
    read1 = r'./example_data/seq_data/A1.raw_1.fastq.gz'
    read2 = r'./example_data/seq_data/A1.raw_2.fastq.gz'
    gff_file = r'./example_data/annotation_file/CLB_strain.gff'
    sample_name = 'CaoLB'
    ref_ps = r'./example_data/annotation_file/CLB_strain.fa'
    adapters = ['AGATCGGAAGAGC', 'AGATCGGAAGAGC']
    sample = RNASeqAnalyzer(sample_name=sample_name, ref_ps=ref_ps, gff_ps=gff_file,
                            seq_ps1=read1, seq_ps2=read2, bowtie_pars={"-p": 32}, adapter=adapters)
    sample.seq_data_align()
    sample.counts_statistic()
