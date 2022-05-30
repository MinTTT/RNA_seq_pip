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
import sciplot as splt
splt.whitegrid()
# Own modules
from seq_utility import BAMile, gff_parser, GeneFeature

#%%
bam_objs = []
for sample in np.arange(1, 6):
    bam_ps = f"/home/fulab/data2/20210321_RNA_seq/1321_{sample}_0_output/1321_{sample}_0.sorted.bam.bin"
    gff_ps = r'./example_data/annotation_file/NH3.23.gff'
    ref_ps = r"./example_data/annotation_file/NH3.23.fasta"

    # bamfile = BAMile(bam_ps, gff_ps, ref_ps, paired_flag=False)
    bamfile = BAMile()
    bamfile.load_data(bam_ps)
    bamfile.set_gene_features(gff_ps)
    # bamfile.count_coverage()
    # bamfile.dump_data()
    genes_dict = bamfile.gene_features
    bam_objs.append(bamfile)
    # fetch ori coverage


#%%

# genome_range = [561894, 565429]  # oriC aslA
# genome_range = [561735, 564227]  # aslA



# genome_range = [1593243, 1600695]  # oriC yaiT
# genome_range = [1595019, 1598010]  # oriC yaiT

# genome_range = [2684545, 2693896]  # ter ynaE
# genome_range = [2685221, 2687760]  # ter ynaE
# genome_range = [2685812, 2687322]  # ter ynaE

# genome_range = [573827, 582244]  # oriC yigE
# genome_range = [575507, 579920]  # oriC yigE
# genome_range = [577687, 578968]  # oriC yigE

# genome_range = [653675, 660369]  # oriC yiiD (fabY)
# genome_range = [655979, 657049]  # oriC yiiD (fabY)
# genome_range = [655979, 656830]  # oriC yiiD
genome_range = [2410109, 2412697]  # ptsG
# genome_range = [2057376, 2061654]  # attB lambda

for bamfile in bam_objs:
    # post_cov = bamfile.fetch_coverage('NCM3722', 402000, 602000, '+', 500)
    # neg_cov = bamfile.fetch_coverage('NCM3722', 402000, 602000, '-', 500)
    fig, ax = plt.subplots(2, 1, figsize=(40, 20))
    bamfile.plot_coverage('NCM3722', *genome_range, '+', 200, ax=ax[0])
    bamfile.plot_coverage('NCM3722', *genome_range, '-', 200, ax=ax[0])
    genes_dict = bamfile.gene_features
    print(bamfile.bam_name)
    gene_features_dict = {}
    for typ_key, genes in genes_dict.items():
        for ge in genes:
            try:
                if genome_range[0] < ge.start < genome_range[-1]:
                    gene_features_dict[ge.gene] = ge
            except AttributeError:
                pass

    # ax[1].hlines(0, min, max, color='k', zorder=0)
    for gene_name, gene in gene_features_dict.items():
        try:
            name = gene.gene
        except AttributeError:
            name = gene.gbkey
        if gene.strand == '+':
            start, end = gene.start, gene.length
            color = '#F89388'
        else:
            start, end = gene.end, -gene.length
            color = '#88F8CB'

        if gene.feature == 'CDS':
            ax[1].arrow(start, 0, end, 0, width=3, head_length=gene.length*0.05, label=name, head_width=3,
                        length_includes_head=True, color=color)
        else:
            ax[1].arrow(start, 0, end, 0, width=3, label=name, head_width=3, head_length=gene.length*0.05,
                        length_includes_head=True, color=color)
        text = ax[1].text(np.mean([gene.start, gene.end]), 0, name, alpha=None, fontsize=45,
                          horizontalalignment='center', verticalalignment='center')
    # ax[0].set_ylim(-100, 100)
    ax[0].set_xlim(*genome_range)

    ax[1].set_xlim(ax[0].get_xlim())
    ax[1].set_ylim(-2, 2)
    ax[1].set_yticks([])
    splt.aspect_ratio(1/20, ax[1])
    # ax[1].legend()
    # ax[1].axis()
    ax[1].set_axis_off()

    fig.show()


#%%
from seq_utility import base_position2relative_pos
relat_oric, theta_oric = base_position2relative_pos(502169, 4678046, 502169, 2857171)
relat_ter, theta_ter = base_position2relative_pos(2857171, 4678046, 502169, 2857171)
relat_ynaE, theta_ynaE = base_position2relative_pos(2686363, 4678046, 502169, 2857171)  # ter
# relat_yaiT, theta_yaiT = base_position2relative_pos(1596523, 4678046, 502169, 2857171)
# relat_yigE, theta_yigE = base_position2relative_pos(568598, 4678046, 502169, 2857171)
# relat_yigE, theta_yigE = base_position2relative_pos(656388, 4678046, 502169, 2857171)  # oriC
relat_yiiD, theta_yiiD = base_position2relative_pos(656388, 4678046, 502169, 2857171)
relat_attB, theta_attB = base_position2relative_pos(2059716, 4678046, 502169, 2857171)  # attB

fig2, ax2 = plt.subplots(1, 1, subplot_kw={'projection': 'polar'})
ax2.plot(theta_oric, 1.1, 'o', color='k', zorder=1)
ax2.plot(theta_ter, 1.1, '*', color='k', zorder=1)
ax2.plot(theta_ynaE, 1.1, 's', color='#1f77b4', zorder=1)

ax2.plot(theta_yiiD, 1.1, 's', color='#2ba02b', zorder=1)
ax2.plot(theta_attB, 1.1, 'v', color='r', zorder=1)

ax2.grid(False)
ax2.set_rticks([])
ax2.set_xticks([])

ax2.set_rlim(0, 1.2)
# ax2.set_thetaticks([])
fig2.show()


