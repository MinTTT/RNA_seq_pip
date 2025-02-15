"""
Custom-made code `reannotte_bg.py` was used to compare and modify the annotations in the gb file generated by prokka.
@author: CHU Pan
@mail: pan_chu@outlook.com

"""

#%%
from Bio import SeqIO
from Bio.Seq import Seq
import typing
import types
# from Bio import Align
from Bio.Blast.Applications import NcbiblastnCommandline
from io import StringIO 
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from tqdm import tqdm
from typing import Union


def blastseq(sequence: Union[SeqRecord, list], subject, saveBlast=None, queryFileps="./query.fasta"):
    SeqIO.write(sequence, queryFileps, 'fasta')
    output = NcbiblastnCommandline(query=queryFileps, subject=subject, outfmt=5)()[0]
    if isinstance(sequence, list):
        blast_result_record = NCBIXML.parse(StringIO(output))
    else:
        blast_result_record = NCBIXML.read(StringIO(output))
    
    if saveBlast:
        with open(saveBlast, 'w') as file:
            file.write(output)

    return blast_result_record


def findMinindex(expectList, maxexpect=10)->list:
    EXCEPTS = maxexpect
    minim_indxes = []
    for i, excepts in enumerate(expectList):
        if excepts < EXCEPTS:
            minim_indxes = [i]
            EXCEPTS = excepts
        elif excepts == EXCEPTS:
            minim_indxes.append(i)
    return minim_indxes


def findIdenticals(hsplist, geneLength):
    identicalPctg = [hsp.identities/geneLength for hsp in hsplist]
    identical_index = []
    for i, pctg in enumerate(identicalPctg):
        if pctg == 1.:
            identical_index.append(i)
    return identical_index


def getRank(loc):
    '''Return the loc rank
    '''
    rank = loc.copy()
    rank.sort()
    rank = [rank.index(i) for i in loc]
    return rank
#%%
gbRef_ps = r'/root/RNA_seq_pip/annotation_file/GCF_000005845.2_ASM584v2_genomic.gbff'

gbTarget_ps = r'/root/RNA_seq_pip/annotation_file/NCM3722_fulab_01/NCM3722.001.gbk'

genomeRef = SeqIO.parse(gbRef_ps, 'genbank')
genomeTarget = SeqIO.parse(gbTarget_ps, 'genbank')
genomeTarget = [fet for fet in genomeTarget][0]
genomeRef = [fet for fet in genomeRef][0]
ANNOTATION_OFFSET = 8  # the offset of bases number between target and reference when identify whether two genes are idential features.

#%%
'''
Ref:
{'gene': ['thrL'],
 'locus_tag': ['b0001'],
 'gene_synonym': ['ECK0001'],
 'codon_start': ['1'],
 'transl_table': ['11'],
 'product': ['thr operon leader peptide'],
 'protein_id': ['NP_414542.1'],
 'db_xref': ['UniProtKB/Swiss-Prot:P0AD86',
  'ASAP:ABE-0000006',
  'ECOCYC:EG11277',
  'GeneID:944742'],
 'translation': ['MKRISTTITTTITITTGNGAG']}
 
 Target:
 {'gene': ['aas'],
 'locus_tag': ['b_02884'],
 'EC_number': ['2.3.1.40'],
 'inference': ['ab initio prediction:Prodigal:002006',
  'similar to AA sequence:GCF_000005845.2_ASM584v2_genomic.gbff:NP_417313.1 '],
 'codon_start': ['1'],
 'transl_table': ['11'],
 'product': ['fused 2-acylglycerophospho-ethanolamine acyltransferase/acyl-acyl carrier protein synthetase'],
 'translation': ['MLFSFFRNLCRVLYRVRVTGDTQALKGERVLITPNHVSFIDGILLGLFLPVRPVFAVYTSISQQWYMRWLKSFIDFVPLDPTQPMAIKHLVRLVEQGRPVVIFPEGRITTTGSLMKIYDGAGFVAAKSGATVIPVRIEGAELTHFSRLKGLVKRRLFPQITLHILPPTQVAMPDAPRARDRRKIAGEMLHQIMMEARMAVRPRETLYESLLSAMYRFGAGKKCVEDVNFTPDSYRKLLTKTLFVGRILEKYSVEGERIGLMLPNAGISAAVIFGAIARRRMPAMMNYTAGVKGLTSAITAAEIKTIFTSRQFLDKGKLWHLPEQLTQVRWVYLEDLKADVTTADKVWIFAHLLMPRLAQVKQQPEEEALILFTSGSEGHPKGVVHSHKSILANVEQIKTIADFTTNDRFMSALPLFHSFGLTVGLFTPLLTGAEVFLYPSPLHYRIVPELVYDRSCTVLFGTSTFLGHYARFANPYDFYRLRYVVAGAEKLQESTKQLWQDKFGLRILEGYGVTECAPVVSINVPMAAKPGTVGRILPGMDARLLSVPGIEEGGRLQLKGPNIMNGYLRVEKPGVLEVPTAENVRGEMERGWYDTGDIVRFDEQGFVQIQGRAKRFAKIAGEMVSLEMVEQLALGVSPDKVHATAIKSDASKGEALVLFTTDNELTRDKLQQYAREHGVPELAVPRDIRYLKQMPLLGSGKPDFVTLKSWVDEAEQHDE']}

add db_xref, protein_id, gene_synonym to target
'''

geneTarget = {}
otherTarget = []
for fet in genomeTarget.features:
    try:
        name = fet.qualifiers['gene'][0]
        geneTarget[name] = fet
    except KeyError:
        try:
            name = fet.qualifiers['product'][0]
            geneTarget[name] = fet
        except KeyError:
            try:
                if fet.type == 'rep_origin':
                    name = fet.qualifiers['note'][0]
                    geneTarget[name] = fet
            except KeyError:
                otherTarget.append(fet)
geneRef = {}
otherRef = []
for fet in genomeRef.features:   
    try:
        name = fet.qualifiers['gene'][0]
        geneRef[name] = fet
    except KeyError:
        try:
            name = fet.qualifiers['product'][0]
            geneRef[name] = fet
        except KeyError:
            try:
                if fet.type == 'rep_origin':
                    name = fet.qualifiers['note'][0]
                    geneRef[name] = fet
            except KeyError:
                otherRef.append(fet)

        
# %%
uniqRef = []
for name in list(geneRef.keys()):
    if name not in geneTarget.keys():
        uniqRef.append(name)
print(f'Number of unique genes in Ref: {len(uniqRef)}')


uniqTarget = []
for name in list(geneTarget.keys()):
    if name not in geneRef.keys():
        uniqTarget.append(name)
print(f'Number of unique genes in Target: {len(uniqTarget)}')

blastret = {}
fa_target = r"./annotation_file/NCM3722_fulab_01/NCM3722.001.fna"
fa_ref = r'./annotation_file/GCF_000005845.2_ASM584v2_genomic.fa'
identicalGenesSeq = {}
identicalLocsinTarget = {}
identicalLocRankTarget = {}
identicalLocHSPTarget = {}


alluniqSeq = {}
deletedSeq = {}
for geneName in tqdm(uniqRef):
    # geneName = 'insH10'
    geneSequence = genomeRef[geneRef[geneName].location.start:geneRef[geneName].location.end]
    if geneRef[geneName].strand == -1:
        geneSequence = geneSequence.reverse_complement()
    geneSequence.id = geneName
    geneSequence.name = geneName
    geneSequence.description = ''
    alluniqSeq[geneName]=geneSequence


blast_result_records = blastseq([seq for _, seq in alluniqSeq.items()], fa_target, saveBlast='./blast_target.xml', queryFileps='./query_reference.fa')
blast_result_records = list(blast_result_records)

#%%
allnuiqSeqBlast = {}
for blast_result_record in tqdm(blast_result_records):
    # geneName = 'insH10'
    geneName = blast_result_record.query
    allnuiqSeqBlast[geneName] = blast_result_record
    geneLoc = (geneRef[geneName].location.start + geneRef[geneName].location.end) / 2
    try:
        align = blast_result_record.alignments[0]

        minim_indxes = findIdenticals(align.hsps, blast_result_record.query_length)
        idthspsTarget = [align.hsps[i] for i in minim_indxes]
        targetLoc = [(hsp.sbjct_start + hsp.sbjct_end)/2 for hsp in idthspsTarget]
        targetRank = targetLoc.copy()
        targetRank.sort()
        targetRank = [targetRank.index(i) for i in targetRank]
    

        if len(minim_indxes) == 1:
            hsp = align.hsps[0]
            strand = hsp.strand
            start, end = hsp.sbjct_start, hsp.sbjct_end
            blastret[geneName] = dict(strand=strand, loaction=(start, end))
        elif len(minim_indxes) == 0:
            print(f'{geneName} was deleted in Target.')
            deletedSeq[geneName] = dict(info='mutations or no hits', blast=blast_result_record) # case 1
        else:
            print(f'{geneName} find {len(minim_indxes)} duplication in Target')
            identicalLocHSPTarget[geneName] = idthspsTarget
            identicalGenesSeq[geneName] = alluniqSeq[geneName]
            identicalLocsinTarget[geneName] = targetLoc
            identicalLocRankTarget[geneName] = targetRank
    except IndexError:
        print(f'{geneName} doesn\'t match in Reference.')
        deletedSeq[geneName] = dict(info='no hits.') # case 2

#%%
idetGeneName = list(identicalGenesSeq.keys())
idtcGeneTag = [None] * len(idetGeneName)

taginit = 0
for i, tag in enumerate(idtcGeneTag):
    if tag is None:
        compareName = idetGeneName[i]
        compareSeq = identicalGenesSeq[compareName]
        for j, name in enumerate(idetGeneName):
            seq = identicalGenesSeq[name]
            if seq.seq == compareSeq.seq:
                idtcGeneTag[j] = taginit
        taginit += 1
idtcGeneLoc = {}
classfiedGeneDict = {}
idtcGeneLocRank = {}     
for i in range(max(idtcGeneTag)+1):
    locus = []
    nameList = []
    for jindex, j in enumerate(idtcGeneTag):
        if j == i:
            nameList.append(idetGeneName[jindex])
            geneRecord = geneRef[idetGeneName[jindex]]
            locus.append((geneRecord.location.start + geneRecord.location.end)/2)
    idtcGeneLoc[i] = locus
    idtcGeneLocRank[i] = getRank(locus)
    classfiedGeneDict[i] = nameList
    
for i, geneName in enumerate(idetGeneName):
    geneTag = idtcGeneTag[i]
    genesList = classfiedGeneDict[geneTag]
    geneRank = idtcGeneLocRank[geneTag][genesList.index(geneName)]
    if geneRank <= max(identicalLocRankTarget[geneName]):
        geneindexTarget = identicalLocRankTarget[geneName].index(geneRank)
        hsp = identicalLocHSPTarget[geneName][geneindexTarget]
        strand = hsp.strand
        start, end = hsp.sbjct_start, hsp.sbjct_end
        blastret[geneName] = dict(strand=strand, loaction=(start, end))
    else:
        print(f'{geneName} was deleted in Target.')
        deletedSeq[geneName] = dict(info='identical genes more than locus', blast=allnuiqSeqBlast[geneName]) # case 3         
    
#%% update features in target

addkeys = ['db_xref', 'protein_id', 'gene_synonym']

unaddfeatures = {}
mapped_features = dict(target=[], reference=[])
for name, fet in geneRef.items():
    
    try:
        # if gene have identical neame in target, add annotation
        fetinTarg = geneTarget[name]  # type: SeqFeature
        if fetinTarg.type == 'CDS':
            for key in addkeys:
                fetinTarg.qualifiers[key] = fet.qualifiers[key]
        
    except KeyError:
        try:
            gene_lac = blastret[name]
            # print(f'add feature: {name}')
            fetNew = SeqFeature(**fet.__dict__)  # copy feature
            if gene_lac['strand'][0] != gene_lac['strand'][1]:
                if gene_lac['strand'][0] == 'Plus':  # Plus/ Minus
                    fetNew.location = FeatureLocation(gene_lac['loaction'][1]-1, gene_lac['loaction'][0], strand=-1)
                else:  # Minus / Plus
                    fetNew.location = FeatureLocation(gene_lac['loaction'][1]-1, gene_lac['loaction'][0], strand=1)
            else:  # + / +
                    fetNew.location = FeatureLocation(gene_lac['loaction'][0]-1, gene_lac['loaction'][1], strand=1)
                
            # check if the gene have been identified in Target
            inTarget = False
            
            for fetindex, fetTarget in enumerate(genomeTarget.features):
                if abs(fetTarget.location.start - fetNew.location.start) < ANNOTATION_OFFSET and \
                abs(fetTarget.location.end - fetNew.location.end) < ANNOTATION_OFFSET:
                    fetNew.qualifiers['locus_tag'] = fetTarget.qualifiers['locus_tag']
                    for key in fetNew.qualifiers.keys():
                        fetTarget.qualifiers[key] = fetNew.qualifiers[key]
                    inTarget = True
                    mapped_features['reference'].append(fet)   
                    mapped_features['target'].append(fetTarget)
                    try:
                        if fetTarget.qualifiers['gene'] in uniqTarget:
                            uniqTarget.remove(fetTarget.qualifiers['gene'])
                    except KeyError:
                        if fetTarget.qualifiers['product'] in uniqTarget:
                            uniqTarget.remove(fetTarget.qualifiers['product'])
                    
            if inTarget == False:
                # print(fetNew.qualifiers)
                genomeTarget.features.append(fetNew)
            # print(len(genomeTarget.features))
        except KeyError:
            unaddfeatures[name] = fet


#%% Statistics deleted / mapped features in Ref. and unique features in Target

deletedSeqNum = len(deletedSeq)
deletedSeqName = list(deletedSeq.keys())
deletedGeneData = dict( FeatureName=[None]*deletedSeqNum, 
                        Description=[None]*deletedSeqNum, 
                        Info=[None]*deletedSeqNum, 
                        BLASTInFo=[None]*deletedSeqNum)
for i in range(deletedSeqNum):
    name = deletedSeqName[i]
    deletedGeneData['FeatureName'][i] = deletedSeqName[i]
    try:
        deletedGeneData['Description'][i] = geneRef[name].qualifiers['product'][0]
    except KeyError:
        deletedGeneData['Description'][i] = geneRef[name].qualifiers['locus_tag'][0]
    deletedGeneData['Info'][i] = deletedSeq[name]['info']
    if 'blast' in deletedSeq[name].keys():
        blasthsps = deletedSeq[name]['blast'].alignments[0].hsps
        blastMsg = ''
        for hsp in blasthsps:
            blastMsg += f'''Alignment Loc: {hsp.sbjct_start}:{hsp.sbjct_end}, ident: { '%.2f' % ((hsp.identities / (hsp.query_end - hsp.query_start + 1))*100)} %.\n'''
        deletedGeneData['BLASTInFo'][i] = blastMsg

# mapped_features
numMappedFets = len(mapped_features['reference'])
mappedFeatsData = dict(Target=[None]*numMappedFets, Description= [None]*numMappedFets, LocationTarget=[None]*numMappedFets,
                        Reference=[None]*numMappedFets, LocationReference=[None]*numMappedFets)
for i in range(numMappedFets):
    try:
        mappedFeatsData['Target'][i] = mapped_features['target'][i].qualifiers['gene'][0]
        mappedFeatsData['Description'][i] = mapped_features['reference'][i].qualifiers['product'][0]
    except KeyError:
        mappedFeatsData['Target'][i] = mapped_features['reference'][i].qualifiers['locus_tag'][0]

    mappedFeatsData['LocationTarget'][i] = f'''{mapped_features['target'][i].location} '''
    try:
        mappedFeatsData['Reference'][i] = mapped_features['reference'][i].qualifiers['gene'][0]
    except KeyError:
        mappedFeatsData['Reference'][i] = mapped_features['reference'][i].qualifiers['product'][0]
    mappedFeatsData['LocationReference'][i] = f'''{mapped_features['reference'][i].location} '''

alluniqSeqTarget = {}
uniqTargetNum = len(uniqTarget)
for uniqName in uniqTarget:
    uniqfet = geneTarget[uniqName]
    uiqSeq = genomeTarget[uniqfet.location.start:uniqfet.location.end]
    if uniqfet.strand == -1:
        uiqSeq = uiqSeq.reverse_complement()
    uiqSeq.name = uniqName
    uiqSeq.id = uniqName
    uiqSeq.description = ''

    alluniqSeqTarget[uniqName] = uiqSeq

blast_result_records = blastseq([seq for _, seq in alluniqSeqTarget.items()], fa_ref, saveBlast='./blast_reference.xml', queryFileps='./query_uiqtarget.fa')
blast_result_records = list(blast_result_records)

uniqGeneData = dict( FeatureName=[None]*uniqTargetNum, 
                        BLASTInFo=[None]*uniqTargetNum)
for i in range(uniqTargetNum):
    name = uniqTarget[i]
    uniqGeneData['FeatureName'][i] = name
    if len(blast_result_records[i].alignments) != 0:
        blasthsps = blast_result_records[i].alignments[0].hsps
        blastMsg = ''
        for hsp in blasthsps:
            fetName = 'Unknown Seq'
            for fetindex, fetRef in enumerate(genomeRef.features):
                if abs(hsp.sbjct_start - fetRef.location.start) < ANNOTATION_OFFSET and abs(hsp.sbjct_end - fetRef.location.end) < ANNOTATION_OFFSET:
                    try:
                        fetName = fetRef.qualifiers['gene'][0]
                    except KeyError:
                        try:
                            fetName = fetRef.qualifiers['product'][0]
                        except KeyError:
                            try:
                                fetName = fetRef.qualifiers['locus_tag'][0]
                            except KeyError:
                                fetName = f'{[it for _,it in  fetRef.qualifiers.items()]}'

                    
        blastMsg += f'''{fetName}, Alignment Loc: {hsp.sbjct_start}:{hsp.sbjct_end}, ident: { '%.2f' % ((hsp.identities / (hsp.query_end - hsp.query_start + 1))*100)} %.\n'''
        uniqGeneData['BLASTInFo'][i] = blastMsg
    else:
        uniqGeneData['BLASTInFo'][i] = 'No hits.'



    
# %%
SeqIO.write([genomeTarget], r'/root/RNA_seq_pip/annotation_file/NCM3722_fulab_01/NCM3722.002.gbk', 'genbank')

import pandas as pd
with pd.ExcelWriter("/root/RNA_seq_pip/annotation_file/NCM3722_fulab_01/NCM3722.002.compareSummary.xlsx", mode='w') as writer:
    pd.DataFrame(deletedGeneData).to_excel(writer, 'Deleted Features')
    pd.DataFrame(mappedFeatsData).to_excel(writer, 'Mapped Features')
    pd.DataFrame(uniqGeneData).to_excel(writer, 'Unique Features in NCM3722')


# %%
