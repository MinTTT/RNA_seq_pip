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
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np  # Or any other
import pysam
# import pysamstats
from tqdm import tqdm
from joblib import Parallel, delayed
from typing import Union, Tuple, List, Dict
import subprocess as sbps


# […]


class GeneFeature:
    def __init__(self, feature, genome, start, end, strand, annotation):
        self.genome = genome
        self.feature = feature
        self.start = start
        self.end = end
        self.strand = strand
        self.parse_annotation(annotation)
        self.length = self.end - self.start + 1

    def parse_annotation(self, annotation):
        for anno in annotation.split(';'):
            self.__dict__[anno.split('=')[0]] = anno.split('=')[-1]


def bed_writer(ps, features: List[GeneFeature]) -> None:
    with open(ps, 'w') as bed_fil:
        lines = [f'{feature.genome}\t{feature.start}\t{feature.end}\n' for feature in features]
        bed_fil.writelines(lines)
    return None


def gff_parser(gff_ps: str) -> dict:
    with open(gff_ps, 'r') as file:
        content = file.readlines()

    lines = [line.strip('\n') for line in content]
    gene_features = []
    for line in lines:
        if len(line) > 0:
            if line[0] == '#':
                pass
            else:
                tokens = line.split('\t')
                genome, feature, start, end, strand, annotate = tokens[0], tokens[2], int(tokens[3]), int(tokens[4]), \
                                                                tokens[6], tokens[-1]
                gene_features.append(GeneFeature(feature, genome, start, end, strand, annotate))

    feature_set = list(set([gene.feature for gene in gene_features]))
    gene_dict = {feature: [] for feature in feature_set}

    for gene in gene_features:
        gene_dict[gene.feature].append(gene)

    return gene_dict


def check_reverse(reads):
    if reads.is_proper_pair:
        read1 = reads.is_read1 and reads.is_reverse
        read2 = reads.is_read2 and (not reads.is_reverse)
        return read1 or read2
    else:
        return False


def check_forward(reads):
    if reads.is_proper_pair:
        read1 = reads.is_read1 and (not reads.is_reverse)
        read2 = reads.is_read2 and reads.is_reverse
        return read1 or read2
    else:
        return False


class BAMile:
    def __init__(self, bam_ps, gff_ps=None, threads=16):
        self.bam_ps = bam_ps
        self.bam_name = os.path.basename(self.bam_ps)
        self.dir = os.path.split(self.bam_ps)[0]
        self.threads = threads
        self.bam = pysam.AlignmentFile(self.bam_ps, 'rb', threads=threads, check_sq=False)  # type: pysam.AlignmentFile
        self.reads_forward_strand_ps = None
        self.reads_reverse_strand_ps = None
        self.forward_coverage_ps = None
        self.reverse_coverage_ps = None
        self.forward_coverage_data = None
        self.reverse_coverage_data = None
        self.genome_set = None
        self.coverage_all = None
        self.mapped_reads = self.bam.mapped
        self.cleaned_reads = None
        self.remove_fwd_bed_ps = None
        self.remove_rev_bed_ps = None
        self.cleaned_reads_forward_strand_ps = None
        self.cleaned_reads_reverse_strand_ps = None
        self.removed_fwd_bam_ps = None
        self.removed_rev_bam_ps = None
        self.rtRNA_clean_flag = False
        self.files = None

        if gff_ps is not None:
            self.gene_features = gff_parser(gff_ps)  # type: dict
        else:
            self.gene_features = None  # type: dict or None
        self.check_path()

    def check_path(self):
        """
        check that whether the files for processing mapped reads is in file.
        :return:
        """
        self.fmt_print('Checking process files.')
        self.files = os.listdir(self.dir)
        if self.bam_name + '.fwd_strand.bam' in self.files:
            self.fmt_print('Detected forward strand mapped sam file')
            self.reads_forward_strand_ps = self.bam_ps + '.fwd_strand.bam'
        if self.bam_name + '.rvs_strand.bam' in self.files:
            self.fmt_print('Detected reverse strand mapped sam file')
            self.reads_reverse_strand_ps = self.bam_ps + '.rvs_strand.bam'
        if self.bam_name + '.fwd_depth.tsv' in self.files:
            self.fmt_print('Loading forward strand coverage.')
            self.forward_coverage_ps = self.bam_ps + '.fwd_depth.tsv'
            fwd_tsv = pd.read_csv(self.forward_coverage_ps, sep='\t', names=['genome', 'location', 'coverage'])
            genome_set = list(set(fwd_tsv['genome'].tolist()))
            self.forward_coverage_data = {genome: fwd_tsv[fwd_tsv['genome'] == genome] for genome in genome_set}
        if self.bam_name + '.rev_depth.tsv' in self.files:
            self.fmt_print('Loading reverse strand coverage.')
            self.reverse_coverage_ps = self.bam_ps + '.rev_depth.tsv'
            rev_tsv = pd.read_csv(self.reverse_coverage_ps, sep='\t', names=['genome', 'location', 'coverage'])
            genome_set = list(set(rev_tsv['genome'].tolist()))
            self.reverse_coverage_data = {genome: rev_tsv[rev_tsv['genome'] == genome] for genome in genome_set}
            self.genome_set = genome_set

    def clean_rtRNA(self):
        # all coverage
        self.rtRNA_clean_flag = True
        self.fmt_print("Cleaning rRNA and tRNA reads.")
        rtRNA_list = []
        try:
            rRNAs = self.gene_features['rRNA']  # type List[GeneFeature]
            rtRNA_list += rRNAs
            rRNA_reads = np.sum(
                [self.count_cds_reads(rRNA.genome, rRNA.start, rRNA.end, rRNA.strand).sum() for rRNA in rRNAs])
        except AttributeError:
            rRNA_reads = 0
            self.fmt_print("Attention: Don't find rRNA annotations.")
        try:
            tRNAs = self.gene_features['tRNA']  # type List[GeneFeature]
            rtRNA_list += tRNAs
            tRNA_reads = np.sum(
                [self.count_cds_reads(tRNAs.genome, tRNAs.start, tRNAs.end, tRNAs.strand) for tRNAs in tRNAs])
        except AttributeError:
            tRNA_reads = 0
            self.fmt_print("Attention: Don't find tRNA annotations.")

        rtRNA_in_forward = [feature for feature in rtRNA_list if feature.strand == '+']
        rtRNA_in_reverse = [feature for feature in rtRNA_list if feature.strand == '-']

        self.remove_fwd_bed_ps = os.path.join(self.dir, 'out_fwd.bed')
        self.remove_rev_bed_ps = os.path.join(self.dir, 'out_rev.bed')
        self.cleaned_reads_forward_strand_ps = self.reads_forward_strand_ps + '.cleaned.bam'
        self.cleaned_reads_reverse_strand_ps = self.reads_reverse_strand_ps + '.cleaned.bam'
        self.removed_fwd_bam_ps = self.reads_forward_strand_ps + '.removed.bam'
        self.removed_rev_bam_ps = self.reads_reverse_strand_ps + '.removed.bam'
        if os.path.basename(self.remove_fwd_bed_ps) not in self.files:
            bed_writer(self.remove_fwd_bed_ps, rtRNA_in_forward)
            bed_writer(self.remove_rev_bed_ps, rtRNA_in_reverse)
        # cmd for clean the reads
        if os.path.basename(self.cleaned_reads_forward_strand_ps) not in self.files:
            if self.reads_forward_strand_ps:
                cmd = f'samtools view {self.reads_forward_strand_ps} -b -h -o {self.removed_fwd_bam_ps} ' \
                      f'-U {self.cleaned_reads_forward_strand_ps} -L {self.remove_rev_bed_ps} -@ {self.threads}'
                status = sbps.run(cmd, shell=True)
        if os.path.basename(self.cleaned_reads_reverse_strand_ps) not in self.files:
            if self.reads_reverse_strand_ps:
                cmd = f'samtools view {self.reads_reverse_strand_ps} -b -h -o {self.removed_rev_bam_ps} ' \
                      f'-U {self.cleaned_reads_reverse_strand_ps} -L {self.remove_fwd_bed_ps} -@ {self.threads}'
                status = sbps.run(cmd, shell=True)
        self.cleaned_reads = self.mapped_reads - rRNA_reads - tRNA_reads

    def fmt_print(self, msg):
        print(f'[{self.bam_name}] -> {msg}')

    def count_cds_reads(self, genome, start, end, strand=None):
        """
        the method return the reads number mapped to the gene feature.
        :param genome: chromosome name
        :param start: gene start site, 1 based
        :param end: gene end site, 1 based
        :param strand: '+' ,'-'. default is None, counting counts along both strand.
        :return: int
        """
        if strand == '+':
            return self.bam.count(contig=genome, start=start - 1, end=end, read_callback=check_reverse)
        elif strand == '-':
            return self.bam.count(contig=genome, start=start - 1, end=end, read_callback=check_forward)
        else:
            return self.bam.count(contig=genome, start=start - 1, end=end)

    def separate_bam_by_strand(self, clean_rtRNA=True):
        if (self.reads_forward_strand_ps is None) or (self.reads_reverse_strand_ps is None):
            self.reads_forward_strand_ps = self.bam_ps + '.fwd_strand.bam'
            cmd_sep = f"samtools view -f 99 -@ {self.threads} {self.bam_ps} -o {self.bam_ps + '.fwd_1.bam'} ;" \
                      f" samtools view -f 147 -@ {self.threads} {self.bam_ps} -o {self.bam_ps + '.fwd_2.bam'} ; " \
                      f"samtools merge -@ {self.threads} {self.reads_forward_strand_ps} " \
                      f"{self.bam_ps + '.fwd_1.bam'} {self.bam_ps + '.fwd_2.bam'}; " \
                      f"rm {self.bam_ps + '.fwd_1.bam'}; rm {self.bam_ps + '.fwd_2.bam'}"
            self.fmt_print(f'Separate bam file: {cmd_sep}.')
            status = sbps.run(cmd_sep, shell=True)

            self.reads_reverse_strand_ps = self.bam_ps + '.rvs_strand.bam'
            cmd_sep = f"samtools view -f 83 -@ {self.threads} {self.bam_ps} -o {self.bam_ps + '.rev_1.bam'} ;" \
                      f" samtools view -f 163 -@ {self.threads} {self.bam_ps} -o {self.bam_ps + '.rev_2.bam'} ; " \
                      f"samtools merge -@ {self.threads} {self.reads_reverse_strand_ps} " \
                      f"{self.bam_ps + '.rev_1.bam'} {self.bam_ps + '.rev_2.bam'}; " \
                      f"rm {self.bam_ps + '.rev_1.bam'}; rm {self.bam_ps + '.rev_2.bam'}"
            self.fmt_print(f'Separate bam file: {cmd_sep}.')
            status = sbps.run(cmd_sep, shell=True)
        if clean_rtRNA:
            self.rtRNA_clean_flag = True
            self.clean_rtRNA()

        return None

    def count_coverage(self):
        if self.forward_coverage_ps is None:
            self.forward_coverage_ps = self.bam_ps + '.fwd_depth.tsv'
            if self.rtRNA_clean_flag:
                cmd_depth = f"samtools mpileup -a {self.cleaned_reads_forward_strand_ps} " \
                            f"| cut -f 1,2,4 > {self.forward_coverage_ps}"
            else:
                cmd_depth = f"samtools mpileup -a {self.reads_forward_strand_ps} " \
                            f"| cut -f 1,2,4 > {self.forward_coverage_ps}"
            print(f"[{os.path.basename(self.bam_ps)}] -> Forward strand reads coverage: {cmd_depth}")
            status = sbps.run(cmd_depth, shell=True)
        fwd_tsv = pd.read_csv(self.forward_coverage_ps, sep='\t', names=['genome', 'location', 'coverage'])
        genome_set = list(set(fwd_tsv['genome'].tolist()))
        self.forward_coverage_data = {genome: fwd_tsv[fwd_tsv['genome'] == genome] for genome in genome_set}
        if self.reverse_coverage_ps is None:
            self.reverse_coverage_ps = self.bam_ps + '.rev_depth.tsv'
            if self.rtRNA_clean_flag:
                cmd_depth = f"samtools mpileup -a {self.cleaned_reads_reverse_strand_ps} " \
                            f"| cut -f 1,2,4 > {self.reverse_coverage_ps}"
            else:
                cmd_depth = f"samtools mpileup -a {self.reads_reverse_strand_ps} " \
                            f"| cut -f 1,2,4 > {self.reverse_coverage_ps}"
            print(f"[{os.path.basename(self.bam_ps)}] -> Reverse strand reads coverage: {cmd_depth}")
            status = sbps.run(cmd_depth, shell=True)
        rev_tsv = pd.read_csv(self.reverse_coverage_ps, sep='\t', names=['genome', 'location', 'coverage'])
        genome_set = list(set(rev_tsv['genome'].tolist()))
        self.reverse_coverage_data = {genome: rev_tsv[rev_tsv['genome'] == genome] for genome in genome_set}
        self.genome_set = genome_set
        self.coverage_all = np.sum(
            [self.reverse_coverage_data[genome]['coverage'].sum() for genome in self.genome_set]) + \
                            np.sum([self.forward_coverage_data[genome]['coverage'].sum() for genome in self.genome_set])
        for genome in self.genome_set:
            self.reverse_coverage_data[genome]['coverage'] /= (self.coverage_all / 1e9)
            self.forward_coverage_data[genome]['coverage'] /= (self.coverage_all / 1e9)
        return None

    def fetch_coverage(self, genome, start, end, strand=None, move_average: int = None):
        """
        fetch coverage from coverage data, execute count_coverage before this method.
        :param genome:
        :param start:
        :param end:
        :param strand:
        :param move_average:
        :return:
        """
        length = end - start + 1
        if strand == '+':
            genome_coverage = self.reverse_coverage_data[genome]['coverage']  # type: pd.Series
        elif strand == '-':
            genome_coverage = self.forward_coverage_data[genome]['coverage']  # type: pd.Series
        else:
            genome_coverage = self.reverse_coverage_data[genome]['coverage'] + \
                              self.forward_coverage_data[genome]['coverage']
        if (move_average is None) or (move_average == 0):
            return np.roll(genome_coverage.values, -(start - 1))[:length]  # genome_coverage.values[start - 1:end]
        else:
            raw_coverage = np.roll(genome_coverage.values, -(start - 1 - move_average))[
                           :length + move_average]  # genome_coverage.values[start - 1 - move_average:end + move_average]
            avg_coverage = np.array(
                [np.mean(raw_coverage[i:i + move_average]) for i in range(move_average + end - start + 1)])
            return avg_coverage[move_average:move_average + end - start + 1]

    def plot_coverage(self, genome, start, end, strand=None, window=10, ax: plt.axes = None, **kwargs):
        if ax is None:
            ax = plt.gca()
        feature_x = np.arange(start, end + 1, step=1)
        if strand is not None:
            coverage = self.fetch_coverage(genome, start, end, strand, move_average=window)
            if strand == '+':
                ax.bar(feature_x, coverage, width=1, **kwargs)
            if strand == '-':
                ax.bar(feature_x, -coverage, width=1, **kwargs)
            return np.hstack([feature_x.reshape(-1, 1), coverage.reshape(-1, 1)])
        else:
            coverage_sense = self.fetch_coverage(genome, start, end, strand='+', move_average=window)
            coverage_antisense = self.fetch_coverage(genome, start, end, strand='-', move_average=window)
            ax.bar(feature_x, coverage_sense, width=1, **kwargs)
            ax.bar(feature_x, -coverage_antisense, width=1, **kwargs)
            return np.hstack(
                [feature_x.reshape(-1, 1), coverage_sense.reshape(-1, 1), coverage_antisense.reshape(-1, 1)])


def sum_of_coverage(cds: GeneFeature, bam: BAMile):
    return np.sum(bam.fetch_coverage(cds.genome, cds.start, cds.end, cds.strand))


def count_reads_custom(bam_ps: str, gff_ps: str, feature: str = 'CDS') -> Tuple[pd.DataFrame, BAMile]:
    bamflie = BAMile(bam_ps, gff_ps)
    cds_list = bamflie.gene_features[feature]
    bamflie.separate_bam_by_strand()  # separate the sam file according to the strands
    bamflie.count_coverage()  # counting coverage along the genome
    reads_stat = []
    print(f"[{os.path.basename(bam_ps)}] -> Counting reads")
    for cds in tqdm(cds_list):
        info = [cds.gene, cds.locus_tag, cds.product, cds.length, cds.strand,
                bamflie.count_cds_reads(cds.genome, cds.start, cds.end, cds.strand)]
        reads_stat.append(info)
    print(f"[{os.path.basename(bam_ps)}] -> Counting coverage")
    coverage = Parallel(n_jobs=-1, require='sharedmem')(
        delayed(sum_of_coverage)(cds, bamflie) for cds in tqdm(cds_list))

    data_frame = pd.DataFrame(data=reads_stat,
                              columns=['gene', 'locus_tag', 'product', 'length', 'strand', 'counts'])
    data_frame['coverage'] = np.array(coverage)
    all_counts = bamflie.cleaned_reads
    data_frame['RPKM'] = data_frame['counts'] * (1000 / data_frame['length']) * (1e6 / all_counts)
    data_frame['TPM'] = data_frame['RPKM'] * 1e6 / data_frame['RPKM'].sum()
    data_frame['coverage_FPKM'] = data_frame['coverage'] / data_frame['length']
    data_frame['coverage_TPM'] = data_frame['coverage_FPKM'] * 1e6 / data_frame['coverage_FPKM'].sum()
    # data_frame.to_csv(bam_ps + '.custom.cds_counts.csv')
    return data_frame, bamflie


def count_reads_htseq(bam_ps: Union[str, list], gff_ps: str, feature: str = 'CDS', reverse: str = 'reverse'):
    if isinstance(bam_ps, str):
        bam_ps = [bam_ps]

    id = 'locus_tag'
    export_name = os.path.join(os.path.split(bam_ps[0])[0], 'htseq.cds_counts.tsv')
    if 'htseq.cds_counts.tsv' not in os.listdir(os.path.split(bam_ps[0])[0]):
        cmd_htseq = f"htseq-count -f bam -s {reverse} -r pos -t {feature} -i {id} -m union " \
                    f"--nonunique fraction -a 10 -n 32 -q {' '.join(bam_ps)} {gff_ps}" \
                    + f" > {export_name}"  # the prefix was used for init env.
        print(cmd_htseq)
        status = sbps.run(cmd_htseq, shell=True)
    columns_title = [id]
    for i in range(len(bam_ps)):
        columns_title.append(f'htseq_counts_{i}')
    count_tsv = pd.read_csv(export_name, sep='\t', names=columns_title)
    count_tsv = pd.concat([count_tsv.iloc[:, 0], count_tsv.iloc[:, 1:].apply(np.sum, axis=1)],
                          axis=1)
    count_tsv.columns = [id, "htseq_counts"]
    reads_info = ['__no_feature', '__ambiguous', '__too_low_aQual', '__not_aligned', '__alignment_not_unique']
    cds_filter = lambda row: False if row in reads_info else True
    cds_mask = count_tsv[id].apply(cds_filter)
    reverse_cds_mask = [i for i in map(lambda x: not x, cds_mask)]
    cds_tsv = count_tsv.loc[cds_mask]
    reads_info_tsv = count_tsv.loc[reverse_cds_mask]
    # print(cds_tsv)
    # print(reads_info_tsv)
    all_reads = np.sum(cds_tsv['htseq_counts']) + reads_info_tsv[reads_info_tsv[id] == '__no_feature']['htseq_counts'].values \
                + reads_info_tsv[reads_info_tsv[id] == '__ambiguous']['htseq_counts'].values
    return cds_tsv, all_reads


def count_feature_reads(bam_ps: str, gff_ps: str, feature: str = 'CDS') -> Tuple[pd.DataFrame, BAMile]:
    custom_counts, bamflie = count_reads_custom(bam_ps, gff_ps, feature)
    bamflie.fmt_print('HTseq counting.')
    htseq_counts, all_reads = count_reads_htseq([bamflie.cleaned_reads_reverse_strand_ps,
                                                 bamflie.cleaned_reads_forward_strand_ps],
                                                gff_ps, feature)
    data_frame = pd.merge(custom_counts, htseq_counts, on='locus_tag', how='left')
    data_frame['htseq_FPKM'] = data_frame['htseq_counts'] * (1000 / data_frame['length']) * (
            1e6 / all_reads)
    data_frame['htseq_TPM'] = data_frame['htseq_FPKM'] * 1e6 / data_frame['htseq_FPKM'].sum()
    # data_frame.to_csv(bam_ps + '.combined.cds_counts.csv')
    return data_frame, bamflie


# %%
if __name__ == '__main__':
    # %%
    bam_ps = r'./A1_andong_output/A1_andong.sorted.bam'
    gff_ps = r'./example_data/annotation_file/GCA_000005845.2_ASM584v2_genomic.gff'
    bams_ps = [r"./A1_andong_output/A1_andong.sorted.bam.fwd_strand.bam.cleaned.bam",
               r"./A1_andong_output/A1_andong.sorted.bam.rvs_strand.bam.cleaned.bam"]
    data, all_counts = count_feature_reads(bam_ps, gff_ps)
    # print(all_counts)
    data.to_csv(bam_ps + 'htseq.stats.csv')

    # # %% test the density plot function
    # bamfile = BAMile(bam_ps, gff_ps=gff_ps)
    # genes_dict = gff_parser(gff_ps)
    #
    # fig, ax = plt.subplots(1, 1)
    # gene = genes_dict['CDS'][50]  # type: GeneFeature
    # a = bamfile.plot_coverage(gene.genome, gene.start - 500, gene.end + 500, window=0)
    # # ax.hlines(y=100, xmax=gene.start, xmin=gene.end, linestyles='dotted', color='r', label=gene.gene)
    # try:
    #     name = gene.gene
    # except AttributeError:
    #     name = gene.gbkey
    # if gene.strand == '+':
    #     ax.arrow(gene.start, 0, gene.length, 0, width=10, color='r', label=name)
    # else:
    #     ax.arrow(gene.end, 0, -gene.length, 0, width=10, color='r', label=name)
    # ax.legend()
    # # ax.set_ylim(-1500, 1500)
    # fig.show()

    # # %% test custom-made file
    # genes_dict = gff_parser(gff_ps)
    # bamflie = BAMile(bam_ps)
    # bamflie.separate_bam_by_strand()
    # cds_list = genes_dict['CDS']
    # id = 'locus_tag'
    # reads_stat = []
    # for cds in tqdm(cds_list):
    #     info = [cds.gene, cds.locus_tag, cds.product, cds.length, cds.strand,
    #             bamflie.count_cds_reads(cds.genome, cds.start, cds.end, cds.strand)]
    #     # # print(cds.gene, (cds.genome, cds.start, cds.end, cds.strand))
    #     # count_reads = bamflie.count_cds_reads(cds.genome, cds.start, cds.end, cds.strand)
    #     reads_stat.append(info)
    #
    # data_frame = pd.DataFrame(data=reads_stat,
    #                           columns=['gene', 'locus_tag', 'product', 'length', 'strand', 'counts'])
    # all_counts = np.sum(data_frame['counts'])
    # data_frame['FPKM'] = data_frame['counts'] * (1000 / data_frame['length']) * (1e6 / all_counts)
    # data_frame['TPM'] = data_frame['FPKM'] * 1e6 / data_frame['FPKM'].sum()
    # data_frame.to_csv(bam_ps + '.custom.cds_counts.csv')
    # # %% test htseq
    # reverse = 'reverse'
    # feature = 'CDS'
    # export_name = bam_ps + '.htseq.cds_counts.tsv'
    # id = 'locus_tag'
    # cmd_htseq = f'htseq-count -f bam -s {reverse} -r pos -t {feature} -i {id} -m union -a 10 -n 32 {bam_ps} {gff_ps}' + \
    #             f' > {export_name}'  # the prefix was used for init env.
    # print(cmd_htseq)
    # pipe_cmd = sbps.run(cmd_htseq, shell=True, executable='/bin/bash')
    # count_tsv = pd.read_csv(export_name, sep='\t', names=[id, 'htseq_counts'], comment='_')
    # data_frame = pd.merge(data_frame, count_tsv, on='locus_tag', how='left')
    # data_frame.to_csv(bam_ps + '.combined.cds_counts.csv')
