# -*- coding: utf-8 -*-

"""

 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""

# Built-in/Generic Imports
import os
import sys
import getopt

# Own modules

from BCBio import GFF
from Bio import SeqIO

def gb2gff(in_file, fasta=True):
    dir, basename = os.path.split(in_file)
    basename = basename[:-len(basename.split('.')[-1])]
    gff_file = os.path.join(dir, basename+'gff')
    fasta_file = os.path.join(dir, basename+'fasta')
    with open(in_file) as gbfile:
        gb = SeqIO.parse(gbfile, "genbank")
        with open(gff_file, 'w') as gfffile:
            GFF.write(gb, gfffile)
    if fasta:
        with open(in_file) as gbfile:
            gb = SeqIO.parse(gbfile, "genbank")
            features = [feature for feature in gb]
            print(f"features number: {len(features)}.")
            with open(fasta_file, 'w') as fafile:
                SeqIO.write(features[0], fafile, 'fasta')




#%%
if __name__ == '__main__':
    args = sys.argv[1:]
    opt, args = getopt.getopt(args, 'f:')
    opt_dict = {i[0]: i[1] for i in opt}
    try:
        iffasta = opt_dict['-f']  # type: str
        if iffasta.lower() == 'false':
            fasta_flag = False
        else:
            fasta_flag = True
    except KeyError:
        fasta_flag = True

    gb2gff(args[0], fasta=fasta_flag)

