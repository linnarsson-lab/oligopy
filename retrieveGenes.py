from pyensembl import EnsemblRelease
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO
import os

def generate_fasta(inputexcel, species, folder=''):
    Ensembl = EnsemblRelease(100,species)
    genesbarcodes = pd.read_excel(inputexcel)
    genesbarcodes = genesbarcodes[genesbarcodes['Gene'] != 'nan']
    genesbarcodes = genesbarcodes[pd.isna(genesbarcodes['Gene']) == False]

    print(genesbarcodes)
    for x in genesbarcodes['Gene']:
        try:
            Ensembl.transcript_ids_of_gene_name(x)
        except:
            print('{}: Failed to retrieve gene. Probably this is not the main gene name, try finding the main gene at ensembl and replace it.'.format(x))
            Ensembl.gene_names

    dic_seq = {}

    if folder != '':
        if not folder.endswith('/'):
            folder = f'{folder}/'
            
    fileout = open(f'{folder}{inputexcel.split(".")[0]}_Markers.fasta', 'w')

    if species == 'human':
        p = os.path.dirname(os.path.abspath(__file__))
        p = os.path.join(p,"HUMAN_NCBI_GENES_retrieved.fasta.masked")
        for record in SeqIO.parse(p, "fasta"):
            gene = record.id.split('|')[0].split('_')[0]
            if gene in list(genesbarcodes['Gene']):
                #print(gene)
                dic_seq[gene] = record.seq
                fileout.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
    elif species == 'mouse':
        p = os.path.dirname(os.path.abspath(__file__))
        p = os.path.join(p,"MOUSE_NCBI_GENES_retrieved.fasta.masked")
        for record in SeqIO.parse(p, "fasta"):
            gene = record.id.split('|')[0].split('_')[0]
            if gene in list(genesbarcodes['Gene']):
                #print(gene)
                dic_seq[gene] = record.seq
                fileout.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
    elif species == 'drosophila':
        print('\n\nNEW SPECIES: DROSOPHILA!!!\n\n')
        p = os.path.dirname(os.path.abspath(__file__))
        p = os.path.join(p,"DROSOPHILA_NCBI_GENES_retrieved.fasta.masked") #TODO: NOT SURE THIS IS CORRECT
        for record in SeqIO.parse(p, "fasta"):
            gene = record.id.split('|')[0].split('_')[0]
            if gene in list(genesbarcodes['Gene']):
                #print(gene)
                dic_seq[gene] = record.seq
                fileout.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
    
    tsl = {1:100,2:80,3:70,4:60,5:50,None:0}
    isproteincoding = {True:100,False:0}
    length = 0
    for x in genesbarcodes['Gene']:
        if x not in dic_seq:
            #print(x)
            '''tids = Ensembl.transcript_ids_of_gene_name(x)
            longest = Ensembl.transcript_sequence(tids[0])
            for x in tids[1:]:
                newseq = Ensembl.transcript_sequence(x)
                if len(newseq) > len(longest):
                    longest = newseq

            '''
            tids = Ensembl.transcript_ids_of_gene_name(x)
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




        