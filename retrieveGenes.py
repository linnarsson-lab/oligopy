from pyensembl import EnsemblRelease
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
import os

def generate_fasta(inputexcel, species, ensembl_release, transcript_reference_file, log, folder=''):
    Ensembl = EnsemblRelease(ensembl_release, species) #Change to desired ensembl release version in variables.ini
    log.info(f'Using Ensembl realease: {ensembl_release}, for species: {species.latin_name}')
    genesbarcodes = pd.read_excel(inputexcel)
    genesbarcodes = genesbarcodes[genesbarcodes['Gene'] != 'nan']
    genesbarcodes = genesbarcodes[pd.isna(genesbarcodes['Gene']) == False]

    #Check if gene names are the main gene names
    error = False
    #log.info(genesbarcodes)
    for x in genesbarcodes['Gene']:
        try:
            Ensembl.transcript_ids_of_gene_name(x)
        except:
            error = True
            log.info(f'ERROR Failed to retrieve gene: "(x)". Probably this is not the main gene name, try finding the main gene at ensembl and replace it.')
            Ensembl.gene_names
    if error == True:
        log.info("""\nCould not find all the genes in Ensembl. Fix the name for the above mentioned genes.
If all genes return an error, you probably did not install the species in pyensemble.
Please run: `pyensembl install --release <release_number> --species <species_latin_name>`
Replace <release_number> with the required Ensembl release.
Replace <species_latin_name> with species name, like: mus_musculus
For more information see the `variables.ini` file in the oligopy folder.""")

    dic_seq = {}

    if folder != '':
        if not folder.endswith('/'):
            folder = f'{folder}/'
            
    fileout = open(f'{folder}{inputexcel.split(".")[0]}_Markers.fasta', 'w')

    #First try to fetch the sequence from the transcript_reference_file
    count = 0
    target_gene_list = list(genesbarcodes['Gene'])
    log.info(sorted(target_gene_list))
    for record in SeqIO.parse(transcript_reference_file, "fasta"):
        gene = record.id.split('|')[0].split('_')[0]
        if gene in target_gene_list:
            dic_seq[gene] = record.seq
            fileout.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
            count += 1
    log.info(f'Found sequences for {len(dic_seq)} out of {len(target_gene_list)} target genes in: {transcript_reference_file}')
    log.info(f'Could not find the following genes in the reference file: {[g for g in target_gene_list if g not in dic_seq.keys()]}')
    log.info('The sequences for these genes will be fetched from Ensemble.')
    
    #If gene was not in the transcript_reference_file fetch the longest sequence from Ensembl
    tsl = {1:100,2:80,3:70,4:60,5:50,None:0}
    isproteincoding = {True:100,False:0}
    length = 0
    for x in target_gene_list: #genesbarcodes['Gene']:
        if x not in dic_seq.keys():
            try:
                tids = Ensembl.transcript_ids_of_gene_name(x)
            except Exception as e:
                raise Exception(f'Gene "{x}" is not a valid gene name. Error message: {e}')
            longest = Ensembl.transcript_sequence(tids[0])
            tmp = Ensembl.transcript_by_id(tids[0])
            supp = tmp.support_level
            cds =tmp.is_protein_coding
            points1 = tsl[supp]+isproteincoding[cds]+length

            for y in tids[1:]:
                newseq = Ensembl.transcript_sequence(y)
                tmp = Ensembl.transcript_by_id(y)
                
                supp = tmp.support_level
                cds =tmp.is_protein_coding

                if len(newseq) > len(longest):
                    length = 5
                    
                points2 = tsl[supp]+isproteincoding[cds]+length
                    
                if points2 > points1:
                    longest = newseq
                    points1 = points2

            dic_seq[x] = Seq(longest)
            fileout.write('>'+str(x)+'\n'+str(longest)+'\n')
