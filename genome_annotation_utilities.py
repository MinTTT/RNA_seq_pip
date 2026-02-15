# -*- coding: utf-8 -*-

"""

 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""
#%%          
# Built-in/Generic Imports
import os
import sys
import getopt
import xml.etree.ElementTree as ET
# Own modules
from BCBio import GFF
from Bio import SeqIO, Entrez
import re
# read email and api key from a json file
try:
    entrez_config = os.path.join(os.path.dirname(__file__), 'entrez.json')
except NameError:
    entrez_config = os.path.join(os.getcwd(), 'entrez.json')
import json
if os.path.exists(entrez_config):
    with open(entrez_config) as config_file:
        config = json.load(config_file)
        # load to env
        os.environ['ENTREZ_CONFIG'] = config['ENTREZ_EMAIL']
        os.environ['ENTREZ_API_KEY'] = config['ENTREZ_API_KEY']


# try to get email from environment variable
if 'ENTREZ_EMAIL' in os.environ:
    Entrez.email = os.environ['ENTREZ_EMAIL']
else:
    Entrez.email = ''
# try to get api key from environment variable
if 'ENTREZ_API_KEY' in os.environ:
    Entrez.api_key = os.environ['ENTREZ_API_KEY']

def read_gb(in_file):
    """
    read .gb file and return the SeqRecord list.

    Parameters
    ----------
    in_file : str
        input genbank file path.

    Returns
    -------
    features : list
        list of SeqRecord objects.

    """
    with open(in_file) as gbfile:
        gb = SeqIO.parse(gbfile, "genbank")
        features = [feature for feature in gb]
    return features



def gb2gff(in_file, fasta=True):
    """
    convert .gb file to .gff file and .fa file.

    Parameters
    ----------
    in_file : str
        input genbank file path.
    fasta : bool, optional
        whether to export fasta file. The default is True.
    Returns
    -------
    Notes
    """
    dir, basename = os.path.split(in_file)
    basename = basename[:-len(basename.split('.')[-1])]
    gff_file = os.path.join(dir, basename + 'gff')
    fasta_file = os.path.join(dir, basename + 'fasta')
    with open(in_file) as gbfile:
        gb = SeqIO.parse(gbfile, "genbank")
        # # check number of records
        # gb_list = [record for record in gb]
        # if len(gb_list) == 1:
        #     print(f"Source length: {len(gb_list[0])}.")
        # else:
        #     print(f"Source length: {len(gb_list[0])}, but there are {len(gb_list)} records in the genbank file.")
        # # check features
        # for record in gb_list:
        #     print(f"Record {record.id} has {len(record.features)} features.")

        with open(gff_file, 'w') as gfffile:
            GFF.write(gb, gfffile)
    if fasta:
        with open(in_file) as gbfile:
            gb = SeqIO.parse(gbfile, "genbank")
            features = [feature for feature in gb]
            print(f"Source length: {len(features[0])}.")
            with open(fasta_file, 'w') as fafile:
                SeqIO.write(features, fafile, 'fasta')
    return None


def get_uniprot_from_ncbi_gene(gene_id:str):
    """
    Find the UniProtKB ID from NCBI Gene database using the GeneID.

    Parameters
    ----------
    gene_id: str
        gene_id from NCBI Gene database.

    Returns
    -------
    uniprot_id: str
        UniProtKB ID corresponding to the GeneID, or None if not found.

    """
    try:
        # Step 1: Fetch the Gene database record in XML format
        handle = Entrez.efetch(db="gene", id=gene_id, retmode="xml")
        record_xml = handle.read()
        handle.close()

        # Step 2: Parse the XML to find the db_xref for UniProtKB
        root = ET.fromstring(record_xml)

        # We look for Dbtag_db == 'UniProtKB/Swiss-Prot'
        for dbtag in root.findall(".//Dbtag"):
            db_name = dbtag.find("Dbtag_db")
            if db_name is not None and db_name.text == "UniProtKB/Swiss-Prot":
                # The actual ID is in the Object-id_str field
                tag_id = dbtag.find(".//Object-id_str")
                if tag_id is not None:
                    return tag_id.text

        return None

    except Exception as e:
        print(f"Error: {e}")
        return None




def gff3togff2(gff3_path, gff2_path):
    """
    # ref https://nbisweden.github.io/AGAT/gxf/#main-points-and-differences-between-gff-formats for
    # the differences between GFF2 and GFF3

    Convert GFF3 file to GFF2/GTF file.

    Parameters
    ----------
    gff3_path
    gff2_path

    Returns
    -------

    """
    print(f'Reading GFF3 file from {gff3_path}...')
    gff3 = GFF.parse(gff3_path)


    with open(gff2_path, 'w') as gff2_file:
        key_excludes = ['ID', 'Name', 'Parent', 'source', 'phase']
        for rec in gff3:
            print(f'Processing record: {rec.id} with {len(rec.features)} features.')
            for feature in rec.features:
                seqname = rec.id
                source = feature.qualifiers.get('source', ['.'])[0]
                ftype = feature.__dict__.get('type', '.')
                start = int(feature.location.start) + 1  # GFF is 1-based
                end = int(feature.location.end)
                score = feature.__dict__.get('score', ['.'])[0]
                strand = feature.location.strand
                if strand == 1:
                    strand = '+'
                elif strand == -1:
                    strand = '-'
                else:
                    strand = '.'
                phase = feature.qualifiers.get('phase', ['.'])[0]
                # attributes
                attributes_keys = list(feature.qualifiers.keys())
                attributes_keys = [k for k in attributes_keys if k not in key_excludes]
                attribute_string = ''
                for k in attributes_keys:
                    values = feature.qualifiers[k]
                    values_str = ','.join(values)
                    # find ECOCYC ID as gene_id
                    if k == 'db_xref':
                        eco_id = re.findall(r'ECOCYC:([A-Za-z0-9_.-]+)', values_str)
                        if eco_id:
                            gene_id = eco_id[0]
                            attribute_string += f'gene_id "{gene_id}"; '
                            continue

                    # rename gene to gene_name
                    if k == 'gene':
                        gene_name = values_str
                        attribute_string += f'gene_name "{gene_name}"; '
                        continue
                    # others only change the format
                    attribute_string += f'{k} "{values_str}"; '
                # raise warning, if have no gene_id, it may cause downstream tools errors
                if ('gene_id' not in attribute_string) and (ftype == 'gene' or ftype == 'CDS'):
                    # if no ECOCYC ID, use the locus_tag as gene_id
                    if 'locus_tag' in attributes_keys:
                        locus_tag = feature.qualifiers['locus_tag'][0]
                        attribute_string += f'gene_id "{locus_tag}"; '
                    else:
                        print(f"Warning: No gene_id found in attributes for feature {attribute_string}")
                # remove trailing space and semicolon
                attribute_string = attribute_string.strip().rstrip(';')
                # replace potential \n in attributes
                attribute_string = attribute_string.replace('\n', ' ')

                gff2_line = f'{seqname}\t{source}\t{ftype}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attribute_string}'
                # print(gff2_line)
                gff2_file.write(gff2_line + '\n')
    return None

# %%
if __name__ == '__main__':
    # args = sys.argv[1:]
    # opt, args = getopt.getopt(args, 'f:')
    # opt_dict = {i[0]: i[1] for i in opt}

    import argparse
    parser = argparse.ArgumentParser(description='Convert genbank file to gff file and fasta file.')
    parser.add_argument('infile', type=str, help='Input genbank file path.')
    parser.add_argument('-f', '--fasta', action='store_true', help='export fasta file.')
    parser.add_argument('-g2', '--gff2', action='store_true', help='export gtf file.')
    parser.add_argument('-2g', '--gff3togff2', action='store_true', help='convert gff3 to gff2/gtf file.')

    parsed_args = parser.parse_args()

    # convert gff3 to gff2/gtf file
    if parsed_args.gff3togff2:
        dir, basename = os.path.split(parsed_args.infile)
        basename = basename[:-len(basename.split('.')[-1])]
        gff3_path = os.path.join(dir, basename + 'gff')
        gff2_path = os.path.join(dir, basename + 'gtf')
        gff3togff2(gff3_path, gff2_path)
        sys.exit(0)
    # convert genbank to gff and fasta
    gb_path = parsed_args.infile
    fasta_flag = parsed_args.fasta
    gb2gff(gb_path, fasta=fasta_flag)
    if parsed_args.gff2:
        dir, basename = os.path.split(gb_path)
        basename = basename[:-len(basename.split('.')[-1])]
        gff3_path = os.path.join(dir, basename + 'gff')
        gff2_path = os.path.join(dir, basename + 'gtf')
        gff3togff2(gff3_path, gff2_path)
