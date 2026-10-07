# -*- coding: utf-8 -*-

"""

 @author: Pan M. CHU
 @Email: pan_chu@outlook.com
"""
#%%          
import os
import argparse
import xml.etree.ElementTree as ET
from pathlib import Path

try:
    from BCBio import GFF
except ImportError:
    GFF = None
try:
    from Bio import SeqIO, Entrez
except ImportError:
    SeqIO = None
    Entrez = None
import re
import json

# Read email and API key from a json file.
try:
    entrez_config = os.path.join(os.path.dirname(__file__), 'entrez.json')
except NameError:
    entrez_config = os.path.join(os.getcwd(), 'entrez.json')
if os.path.exists(entrez_config):
    with open(entrez_config) as config_file:
        config = json.load(config_file)
        if config.get('ENTREZ_EMAIL'):
            os.environ['ENTREZ_EMAIL'] = config['ENTREZ_EMAIL']
        if config.get('ENTREZ_API_KEY'):
            os.environ['ENTREZ_API_KEY'] = config['ENTREZ_API_KEY']


# Try to get Entrez credentials from environment variables.
if Entrez is not None:
    if 'ENTREZ_EMAIL' in os.environ:
        Entrez.email = os.environ['ENTREZ_EMAIL']
    else:
        Entrez.email = ''
    if 'ENTREZ_API_KEY' in os.environ:
        Entrez.api_key = os.environ['ENTREZ_API_KEY']


def require_biopython():
    if SeqIO is None or Entrez is None:
        raise ImportError('Missing dependency: install biopython to use this function.')


def require_gff():
    if GFF is None:
        raise ImportError('Missing dependency: install bcbio-gff to convert GFF/GenBank annotations.')

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
    require_biopython()
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
    require_biopython()
    require_gff()
    in_path = Path(in_file)
    gff_file = in_path.with_suffix('.gff')
    fasta_file = in_path.with_suffix('.fasta')
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


def convert_genbank(in_file, output_dir=None, prefix=None, fasta=False, gff2=False):
    """
    Convert a GenBank file to GFF3, optionally FASTA and GTF/GFF2.

    Parameters
    ----------
    in_file : str or Path
        Input GenBank file path.
    output_dir : str or Path, optional
        Output directory. Defaults to the input file directory.
    prefix : str, optional
        Output filename prefix. Defaults to the input filename stem.
    fasta : bool, optional
        Export FASTA sequence file.
    gff2 : bool, optional
        Export GTF/GFF2 annotation derived from the generated GFF3.

    Returns
    -------
    dict
        Output paths keyed by format name.
    """
    require_biopython()
    require_gff()
    in_path = Path(in_file)
    out_dir = Path(output_dir) if output_dir else in_path.parent
    out_prefix = prefix if prefix else in_path.stem
    out_dir.mkdir(parents=True, exist_ok=True)

    gff3_path = out_dir / f'{out_prefix}.gff'
    fasta_path = out_dir / f'{out_prefix}.fasta'
    gtf_path = out_dir / f'{out_prefix}.gtf'

    with open(in_path) as gbfile:
        gb = SeqIO.parse(gbfile, "genbank")
        with open(gff3_path, 'w') as gfffile:
            GFF.write(gb, gfffile)

    outputs = {'gff3': gff3_path}

    if fasta:
        with open(in_path) as gbfile:
            records = list(SeqIO.parse(gbfile, "genbank"))
        if not records:
            raise ValueError(f'No GenBank records found in {in_path}')
        print(f"Source length: {len(records[0])}.")
        with open(fasta_path, 'w') as fafile:
            SeqIO.write(records, fafile, 'fasta')
        outputs['fasta'] = fasta_path

    if gff2:
        gff3togff2(gff3_path, gtf_path)
        outputs['gtf'] = gtf_path

    return outputs


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
    require_biopython()
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
    require_gff()
    print(f'Reading GFF3 file from {gff3_path}...')

    with open(gff3_path) as gff3_file, open(gff2_path, 'w') as gff2_file:
        gff3 = GFF.parse(gff3_file)
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


def convert_gff3(gff3_path, output_dir=None, prefix=None):
    """
    Convert a GFF3 file to a GTF/GFF2-style file.

    Parameters
    ----------
    gff3_path : str or Path
        Input GFF3 file.
    output_dir : str or Path, optional
        Output directory. Defaults to the input file directory.
    prefix : str, optional
        Output filename prefix. Defaults to the input filename stem.

    Returns
    -------
    pathlib.Path
        Output GTF/GFF2 path.
    """
    in_path = Path(gff3_path)
    out_dir = Path(output_dir) if output_dir else in_path.parent
    out_prefix = prefix if prefix else in_path.stem
    out_dir.mkdir(parents=True, exist_ok=True)
    gtf_path = out_dir / f'{out_prefix}.gtf'
    gff3togff2(in_path, gtf_path)
    return gtf_path


def infer_conversion_mode(infile, explicit_mode=None):
    """
    Infer conversion mode from CLI input when the user did not choose a mode.
    """
    if explicit_mode:
        return explicit_mode

    suffix = Path(infile).suffix.lower()
    if suffix in {'.gb', '.gbk', '.gbff', '.genbank'}:
        return 'genbank'
    if suffix in {'.gff', '.gff3'}:
        return 'gff3'
    raise ValueError(
        f"Cannot infer input type from '{suffix}'. Use --mode genbank or --mode gff3."
    )


def build_parser():
    parser = argparse.ArgumentParser(
        description='Convert genome annotation files between GenBank, GFF3, GTF/GFF2, and FASTA.'
    )
    parser.add_argument('infile', type=Path, help='Input GenBank or GFF3 file path.')
    parser.add_argument(
        '-m', '--mode',
        choices=['auto', 'genbank', 'gff3'],
        default='auto',
        help='Input format. Default: auto-detect from file extension.'
    )
    parser.add_argument(
        '-o', '--output-dir',
        type=Path,
        default=None,
        help='Output directory. Default: same directory as input.'
    )
    parser.add_argument(
        '-p', '--prefix',
        default=None,
        help='Output filename prefix. Default: input filename without final extension.'
    )
    parser.add_argument(
        '-f', '--fasta',
        action='store_true',
        help='For GenBank input, also export FASTA.'
    )
    parser.add_argument(
        '-g2', '--gff2',
        action='store_true',
        help='For GenBank input, also export GTF/GFF2 after generating GFF3.'
    )
    parser.add_argument(
        '-2g', '--gff3togff2',
        action='store_true',
        help='Backward-compatible shortcut for --mode gff3.'
    )
    return parser


def main(argv=None):
    parser = build_parser()
    parsed_args = parser.parse_args(argv)

    if not parsed_args.infile.exists():
        parser.error(f'Input file does not exist: {parsed_args.infile}')
    if not parsed_args.infile.is_file():
        parser.error(f'Input path is not a file: {parsed_args.infile}')

    if parsed_args.gff3togff2:
        if parsed_args.fasta or parsed_args.gff2:
            parser.error('--gff3togff2 cannot be combined with --fasta or --gff2.')
        gff3_input = parsed_args.infile
        if parsed_args.infile.suffix.lower() in {'.gb', '.gbk', '.gbff', '.genbank'}:
            gff3_input = parsed_args.infile.with_suffix('.gff')
        if not gff3_input.exists():
            parser.error(f'GFF3 input file does not exist: {gff3_input}')
        try:
            output_path = convert_gff3(
                gff3_input,
                output_dir=parsed_args.output_dir,
                prefix=parsed_args.prefix,
            )
        except Exception as exc:
            parser.error(str(exc))
        print(f'Wrote gtf: {output_path}')
        return

    requested_mode = parsed_args.mode
    if requested_mode == 'auto':
        requested_mode = None

    try:
        mode = infer_conversion_mode(parsed_args.infile, requested_mode)
        if mode == 'genbank':
            outputs = convert_genbank(
                parsed_args.infile,
                output_dir=parsed_args.output_dir,
                prefix=parsed_args.prefix,
                fasta=parsed_args.fasta,
                gff2=parsed_args.gff2,
            )
            for output_type, output_path in outputs.items():
                print(f'Wrote {output_type}: {output_path}')
        elif mode == 'gff3':
            if parsed_args.fasta or parsed_args.gff2:
                parser.error('--fasta and --gff2 are only valid for GenBank input.')
            output_path = convert_gff3(
                parsed_args.infile,
                output_dir=parsed_args.output_dir,
                prefix=parsed_args.prefix,
            )
            print(f'Wrote gtf: {output_path}')
        else:
            parser.error(f'Unsupported conversion mode: {mode}')
    except Exception as exc:
        parser.error(str(exc))

# %%
if __name__ == '__main__':
    main()
