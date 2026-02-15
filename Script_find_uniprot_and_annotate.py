

#%%
from genome_annotation_utilities import read_gb, get_uniprot_from_ncbi_gene, GFF

import os


#%%
lambda_gb_file = r'G:\OnDrive-CP\OneDrive\Research\Data\gene\Genome\E.coli_K12\NC_001416_phage_lambda_sequence.gb'

gb_features = read_gb(lambda_gb_file)
gb_number = len(gb_features)
geneid_uniprot_dict = {}
if gb_number == 1:
    print(f"Source length: {len(gb_features[0])}.")
else:
    print(f"Source length: {len(gb_features[0])}, but there are {gb_number} records in the genbank file.")
for feature in gb_features[0].features:
    if 'qualifiers' not in feature.__dict__:
        continue
    feature_db_xref = feature.qualifiers.get('db_xref', "")
    if feature_db_xref:
        if 'GeneID' in feature_db_xref[0]:
            print(feature_db_xref[0])
            gene_id = feature_db_xref[0].split(':')[1]
            if gene_id not in geneid_uniprot_dict:
                uniprot_id = get_uniprot_from_ncbi_gene(gene_id)
                geneid_uniprot_dict[gene_id] = uniprot_id
                print(f"GeneID: {gene_id} -> UniProtKB: {uniprot_id}")

    feature_gene_name = feature.qualifiers.get('gene', "")
    if feature_gene_name:
        print(feature_gene_name[0])

# write geneid_uniprot_dict info to features

for feature in gb_features[0].features:
    if 'qualifiers' not in feature.__dict__:
        continue
    feature_db_xref = feature.qualifiers.get('db_xref', "")
    if feature_db_xref:
        if 'GeneID' in feature_db_xref[0]:
            gene_id = feature_db_xref[0].split(':')[1]
            uniprot_id = geneid_uniprot_dict.get(gene_id, None)
            if uniprot_id and uniprot_id != 'No Swiss-Prot link found in Gene record':
                if 'db_xref' in feature.qualifiers:
                    feature.qualifiers['db_xref'].append(f"UniProtKB/Swiss-Prot:{uniprot_id}")
                else:
                    feature.qualifiers['db_xref'] = [f"UniProtKB/Swiss-Prot:{uniprot_id}"]

# reate gff file with updated features
updated_gff_file = os.path.join(os.path.dirname(lambda_gb_file), 'NC_001416_phage_lambda_sequence.gff')
with open(updated_gff_file, 'w') as gfffile:
    GFF.write([gb_features[0]], gfffile)