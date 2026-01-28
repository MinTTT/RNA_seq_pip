# -*- coding: utf-8 -*-

"""

 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""

# Built-in/Generic Imports
import os
# Libs
import subprocess as sbps
import datetime
from typing import Union, Optional
from seq_utility import count_feature_reads, BAMile, gff_parser


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
        """
        Append the output of the command to the log file.
        Parameters
        ----------
        stdout : str or bytes
            The output of the command.

        Returns
        -------

        """
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
        self.indexed_reference_dir = None

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

        # create index files
        # index dir name: Index_{reference file name without suffix}
        index_dir_name = f"Index_{os.path.splitext(os.path.basename(self.reference_file_name))[0]}"
        index_dir_path = os.path.join(self.output_dir, index_dir_name)
        index_flag = True
        if os.path.exists(index_dir_path):
            # check index files
            index_files = os.listdir(index_dir_path)
            for file_name in index_files:
                if file_name.startswith(self.indexed_base_name) and \
                        file_name.split('.')[-1] in ['1.bt2', '2.bt2', '3.bt2', '4.bt2', 'rev.1.bt2', 'rev.2.bt2']:
                    index_flag = False
                    break
        else:
            os.makedirs(os.path.join(self.output_dir, index_dir_name))
        if index_flag:
            self.indexed_reference_dir = os.path.join(self.output_dir, index_dir_name)
            self.append_to_log(f'[{self.sample_name}] -> Create index dir: {self.indexed_reference_dir}')
            lnk_path = os.path.join(self.indexed_reference_dir, 'ref.fa')
            if 'ref.fa' in os.listdir(self.indexed_reference_dir):
                os.remove(os.path.join(self.indexed_reference_dir, 'ref.fa'))

            cmd_lnk_ref = f'ln -s {self.reference_file_path} ' + \
                          f'{lnk_path}'
            print(f'[{self.sample_name}] -> Link Reference: {cmd_lnk_ref}')

            self.append_to_log(f'[{self.sample_name}] -> Link Reference: {cmd_lnk_ref}')
            status1 = self.cmd_shell(cmd_lnk_ref)
            if status1.returncode is not None:
                print(status1.returncode)
                print(f'[{self.sample_name}] -> Error! Link reference failed!')
                self.append_to_log(f'[{self.sample_name}] -> Error! Link reference failed!')
                exit(1)
            self.reference_file_path = lnk_path

            cmd_index = f'bowtie2-build -f {self.reference_file_path} {self.indexed_base_name}'
            print(f'[{self.sample_name}] -> Generate indexed reference: {cmd_index}')
            status2 = self.cmd_shell(cmd_index, cwd=self.indexed_reference_dir)
            # remove link file
            if os.path.exists(lnk_path):
                os.remove(lnk_path)

        else:
            self.indexed_reference_dir = os.path.join(self.output_dir, index_dir_name)
            self.append_to_log(f'[{self.sample_name}] -> Index dir existed: {self.indexed_reference_dir}')

        # cmd_copy_ref = f'cp {self.reference_file_path} ' \
        #                f'{os.path.join(self.output_dir, self.reference_file_name)}'

        # # update the reference file ps
        # self.reference_file_path = os.path.join(self.output_dir, self.reference_file_name)
        # self.append_to_log(f'[{self.sample_name}] -> Copy Reference: {cmd_copy_ref}')
        # status1 = self.cmd_shell(cmd_copy_ref)
        # if status1.returncode != 0:
        #     print(f'[{self.sample_name}] -> Error! Copy reference failed!')
        #     self.append_to_log(f'[{self.sample_name}] -> Error! Copy reference failed!')
        #     exit(1)

        # mapping
        # find bam file first, if not found, do mapping
        if os.path.basename(self.bam_index_ps) not in self.file_in_dir:
            os.environ['BOWTIE2_INDEXES'] = self.indexed_reference_dir
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
            if status3.returncode is not None:
                print(f'[{self.sample_name}] -> Error! Mapping reads failed!')
                self.append_to_log(f'[{self.sample_name}] -> Error! Mapping reads failed!')
                exit(1)
            del os.environ['BOWTIE2_INDEXES']

            # sort alignments
            # Sorts reads by chromosome and then position, crucial for downstream tools.
            cmd_gen_bam = f'samtools view -bS -@ {self.bowtie_pars["-p"]} {self.sam_file_ps} | samtools sort -o {self.bam_sorted_ps} -@ {self.bowtie_pars["-p"]}'
            print(f'[{self.sample_name}] -> Generating BAM: {cmd_gen_bam}')
            self.append_to_log(f'[{self.sample_name}] -> Generating BAM: {cmd_gen_bam}')
            status4 = self.cmd_shell(cmd_gen_bam)
            cmd_index_bam = f'samtools index -@ {self.bowtie_pars["-p"]} {self.bam_sorted_ps} {self.bam_index_ps}'
            print(f'[{self.sample_name}] -> Generating .bai: {cmd_index_bam}')
            status5 = self.cmd_shell(cmd_index_bam)
            if status5.returncode is not None:
                print(f'[{self.sample_name}] -> Error! Generating .bai failed!')
                self.append_to_log(f'[{self.sample_name}] -> Error! Generating .bai failed!')
                exit(1)
            # remove sam file to save space
            if os.path.exists(self.sam_file_ps):
                os.remove(self.sam_file_ps)


    def counts_statistic(self, count_mode='full'):
        """
        Calculate gene expression level from the alignment bam file.
        Parameters
        ----------
        count_mode: str
            'full', 'htseq', 'custom'

        Returns
        -------

        """
        print(f'[{self.sample_name}] -> Calculate gene expression level.')
        counts_stat, bam = count_feature_reads(self.bam_sorted_ps, self.gff_ps, self.reference_file_path,
                                               paired_flag=self.paired_flag, stranded='reverse',
                                               log_writer=self.append_to_log, mode=count_mode)
        # statistics tpm, fpkm, etc.
        # 1. calculate rRNA ratio
        gene_dcit = gff_parser(self.gff_ps)
        rRNA_list = gene_dcit['rRNA']
        tRNA_list = gene_dcit['tRNA']
        # rRNA counts
        rRNA_counts = 0
        for gene_feature in rRNA_list:
            gene_id = gene_feature.__dict__.get('locus_tag', None)
            rRNA_counts += counts_stat[counts_stat['locus_tag'] == gene_id]['htseq_counts'].values[0]
        tRNA_counts = 0
        for gene_feature in tRNA_list:
            gene_id = gene_feature.__dict__.get('locus_tag', None)
            tRNA_counts += counts_stat[counts_stat['locus_tag'] == gene_id]['htseq_counts'].values[0]

        total_counts = counts_stat['htseq_counts'].sum()
        rRNA_ratio = rRNA_counts / total_counts * 100
        tRNA_ratio = tRNA_counts / total_counts * 100
        self.append_to_log(f'[{self.sample_name}] -> rRNA counts: {rRNA_counts}, total counts: {total_counts}, rRNA ratio: {rRNA_ratio:.2f}%')
        self.append_to_log(f'[{self.sample_name}] -> tRNA counts: {tRNA_counts}, total counts: {total_counts}, tRNA ratio: {tRNA_ratio:.2f}%')
        with open(os.path.join(self.output_dir, f'{self.sample_name}_rRNA_tRNA_ratio.txt'), 'w') as ratio_file:
            ratio_file.write(f'rRNA counts\t{rRNA_counts}\ttotal counts\t{total_counts}\trRNA ratio\t{rRNA_ratio:.2f}\n')
            ratio_file.write(f'tRNA counts\t{tRNA_counts}\ttotal counts\t{total_counts}\ttRNA ratio\t{tRNA_ratio:.2f}\n')

        # # statistics CDS counts, and calculate FPKM, TPM, RPM
        # get CDS features
        cds_list = gene_dcit['CDS']
        cds_ids = [cds.__dict__.get('locus_tag', None) for cds in cds_list]

        cds_df = counts_stat[counts_stat['locus_tag'].isin(cds_ids)].copy()
        cds_reads_counts = cds_df['htseq_counts'].sum()
        cds_df['htseq_FPKM'] = cds_df['htseq_counts'] * (1000 / cds_df['length']) * (1e6 / cds_reads_counts)
        cds_df['htseq_TPM'] = 1e6 * cds_df['htseq_counts'] / cds_df['length'] / \
                               (cds_df['htseq_counts'] / cds_df['length']).sum()
        cds_df['htseq_RPM'] = cds_df['htseq_counts'] * (1e6 / cds_reads_counts)
        # write cds info to table
        for cds in cds_list:
            gene_id = cds.__dict__.get('locus_tag', None)
            cds_product = cds.__dict__.get('product', None)
            cds_df.loc[cds_df['locus_tag'] == gene_id, 'product'] = cds_product
        cds_df.to_csv(os.path.join(self.output_dir, f'{self.sample_name}_CDS_expression_statistic.csv'))
        # all_counts = htseq_counts['htseq_counts'].sum()
        # htseq_counts['htseq_FPKM'] = htseq_counts['htseq_counts'] * (1000 / htseq_counts['length']) * (
        #         1e6 / all_counts)
        # # htseq_counts['htseq_TPM'] = htseq_counts['htseq_FPKM'] * 1e6 / htseq_counts['htseq_FPKM'].sum()
        # htseq_counts['htseq_TPM'] = 1e6 * htseq_counts['htseq_counts'] / htseq_counts['length'] / \
        #                             np.sum(htseq_counts['htseq_counts'] / htseq_counts['length'])
        # htseq_counts['htseq_RPM'] = htseq_counts['htseq_counts'] * (1e6 / all_counts)


        self.__dict__['counts_stat'] = counts_stat
        if bam:
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
