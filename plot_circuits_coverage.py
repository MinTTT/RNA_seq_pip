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
import matplotlib.pyplot as plt
# […]

# Own modules
from seq_utility import BAMile, gff_parser, GeneFeature

#%%
ori_coverage_dict = dict(sample_name=[], RNAII_covg=[], RNAI_covg=[], RNAI_CPKM=[], RNAII_CPKM=[])
for sample in np.arange(1, 6):
    bam_ps = f"/home/fulab/data2/20210321_RNA_seq/1321_{sample}_0_output/1321_{sample}_0.sorted.bam.bin"
    gff_ps = r'./example_data/annotation_file/NH3.23.gff'
    ref_ps = r"./example_data/annotation_file/NH3.23.fasta"
    gene_names_list = ['tetR-sfgfp', 'kanR', 'Ptrc', 'gfp', 'lacI', 'ColE1']

    # bamfile = BAMile(bam_ps, gff_ps, ref_ps, paired_flag=False)
    bamfile = BAMile()
    bamfile.load_data(bam_ps)
    bamfile.set_gene_features(gff_ps)
    # bamfile.count_coverage()
    # bamfile.dump_data()
    genes_dict = bamfile.gene_features
    gene_features_dict = {}
    for typ_key, genes in genes_dict.items():
        for ge in genes:
            try:
                if ge.gene in gene_names_list:
                    gene_features_dict[ge.gene] = ge
            except AttributeError:
                pass
    # fetch ori coverage
    ori = gene_features_dict['ColE1']  # type: GeneFeature

    RNAII = bamfile.fetch_coverage(ori.genome, ori.start, ori.end, strand=ori.strand, move_average=100)
    if ori.strand == '-':
        reverse_strand = '+'
    else:
        reverse_strand = '-'
    RNAI = bamfile.fetch_coverage(ori.genome, ori.start, ori.end, strand=reverse_strand, move_average=100)
    ori_coverage_dict['sample_name'].append(sample)
    ori_coverage_dict['RNAII_CPKM'].append(RNAII.mean())
    ori_coverage_dict['RNAI_CPKM'].append(RNAI.mean())
    ori_coverage_dict['RNAII_covg'].append(RNAII.sum())
    ori_coverage_dict['RNAI_covg'].append(RNAI.sum())


    # fig, ax = plt.subplots(2, 1, figsize=(10, 4), gridspec_kw={'height_ratios': [3, 1]})
    # # gene = genes_dict['CDS'][50]  # type: GeneFeature
    # a = bamfile.plot_coverage('NH3.23_plasmid.1', 0, 5602, window=200, ax=ax[0])
    # a[:, -1] *= -1
    # ax[0].plot(a[:, 0], a[:, 1:], color='k')
    # # ax.hlines(y=100, xmax=gene.start, xmin=gene.end, linestyles='dotted', color='r', label=gene.gene)
    # # min, max, length = a[:, 0].min(), a[:, 0].max(), np.ptp(a[:, 0])
    # min, max, length = ori.start - 500, ori.end + 500, ori.length + 1000
    # ax[0].set_xlim(min, max)
    # ax[0].set_ylim(-20000, 20000)
    # ax[0].set_yscale('symlog', linthresh=3000)
    # ax[0].set_yticks([10000, 1000, 0, -1000, -10000])
    # ax[0].set_title(f'medium:{sample}')
    #
    # ax[1].hlines(0, min, max, color='k', zorder=0)
    # for gene_name, gene in gene_features_dict.items():
    #     try:
    #         name = gene.gene
    #     except AttributeError:
    #         name = gene.gbkey
    #     if gene.strand == '+':
    #         start, end = gene.start, gene.length
    #     else:
    #         start, end = gene.end, -gene.length
    #     if gene.feature == 'CDS':
    #         ax[1].arrow(start, 0, end, 0, width=3, head_length=length*0.05, label=name, head_width=6,
    #                     length_includes_head=True)
    #     else:
    #         ax[1].arrow(start, 0, end, 0, width=3, label=name, head_width=3, head_length=length*0.05,
    #                     length_includes_head=True)
    #
    # ax[1].set_xlim(ax[0].get_xlim())
    # ax[1].set_ylim(-10, 10)
    # ax[1].set_yticks([])
    # # ax[1].legend()
    # for skey, spine in ax[1].spines.items():
    #     spine.set_linewidth(0)
    # # ax.set_ylim(-1500, 1500)
    # fig.show()
    # fig.savefig(f'./example_data/annotation_file/sample_circuts_{sample}.png', transparent=True)

ori_coverage = pd.DataFrame(data=ori_coverage_dict)
ori_coverage.to_csv('./example_data/annotation_file/ori_covg_NH3.23.csv')