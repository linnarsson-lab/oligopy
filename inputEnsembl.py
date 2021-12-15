import sys
input_file = sys.argv[1]
import os
from pkg_resources import Requirement, resource_filename
filename_MOUSE = resource_filename(Requirement.parse("oligopy"),"docs/Mus_musculus_annotation_ensembl89.txt")
filename_HUMAN = resource_filename(Requirement.parse("oligopy"),"docs/Homo_sapiens_annotation_ensembl89.txt")

###Genes to obtain sequences: input fasta
genes_to_retrieve = []
genes_with_seq = {}
from Bio import SeqIO
for seq in SeqIO.parse(input_file, "fasta"):
    if len(seq) == 0:
        genes_to_retrieve.append(seq.name)
    else:
        genes_with_seq[seq.name] = seq
#Create class Transcript
if sys.argv[2] == "human":
    input_ensembl = open(filename_HUMAN, "r")
if sys.argv[2] == "mouse":
    input_ensembl = open(filename_MOUSE, "r")

class Transcript:
    def __init__(self,line):
        self.geneID = line[0]
        self.transcriptID = line[1]
        self.gencode = line[2]
        self.appris = line[3]
        self.geneName = line[4]
        self.len = line[5]
        self.level = line[6]
        self.source = line[7]
        self.version = line[8]
        self.count = line[9]
        self.type = line[10]
        self.NCBI = line[11]
transcripts = {}

for elem in input_ensembl:
    line = elem.strip("\n").split("\t")
    gene = line[4]
    if gene not in transcripts:
        transcripts[gene] = [Transcript(line)]
    else:
        transcripts[gene].append(Transcript(line))

def pick_transcript(list):
    type_levels = {"":0,"tslNA": 1000, "tsl5":2000, "tsl4": 3000, "tsl3": 4000, "tsl2": 5000, "tsl1": 6000}
    APPRIS_levels = {"": 0, "alternative2": 25000, "alternative1": 50000, "principal5":100000, "principal4": 200000, "principal3":300000, "principal2": 400000, "principal1": 500000}
    picked_transcript_NCBI = list[0].NCBI
    length = int(list[0].len)
    if list[0].type == "protein_coding":
        prot_score = 1000000
    else:
        prot_score = 0
    if list[0].source == "ensembl_havana":
        source_score = 10000
    else:
        source_score = 0
    try:
        appris_score = APPRIS_levels[list[0].appris]
    except:
        appris_score = 0
    try:
        tsl_score = type_levels[list[0].level]
    except:
        tsl_score = 0
    score = prot_score + source_score + appris_score + tsl_score
    for t in list:
        new_length = int(t.len)
        if t.type == "protein_coding":
            prot_score = 1000000
        else:
            prot_score = 0
        if t.source == "ensembl_havana":
            source_score = 10000
        else:
            source_score = 0
        try:
            appris_score = APPRIS_levels[t.appris]
        except:
            appris_score = 0
        try:
            tsl_score = type_levels[t.level]
        except:
            tsl_score = 0
        if new_length > length:
            length_score = 3000
        else:
            length_score = 0
        new_score = prot_score + source_score + appris_score + tsl_score + length_score
        if new_score >= score and t.NCBI != "":
            picked_transcript_NCBI = t.NCBI
            score = new_score
        if picked_transcript_NCBI == "":
            picked_transcript_NCBI = t.NCBI
            score = new_score
    return picked_transcript_NCBI

NCBI_IDs = {}
for gene in genes_to_retrieve:
    if len(transcripts[gene]) > 1:
        NCBI = pick_transcript(transcripts[gene])
        NCBI_IDs[gene] = NCBI
    else:
        NCBI_IDs[gene] = transcripts[gene][0].NCBI
###Pick the desired sequence according to the criteria in mart_export_ensembl_info.txt


###Retrieve sequences from ensembl database
from Bio import Entrez
from Bio import SeqIO
obtained_seqs = {}
Entrez.email = "A.N.Other@example.com"
for gene in NCBI_IDs:
    print "Gene is: " + gene
    handle =  Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=NCBI_IDs[gene])
    seq_record = SeqIO.read(handle, "fasta")
    print seq_record
    obtained_seqs[gene] = seq_record

out = open(input_file.split(".")[0]+"_retrieved.fasta", "w")
for seq in obtained_seqs:
    out.write(">"+ seq+"_"+ str(obtained_seqs[seq].id) +"\n"+ str(obtained_seqs[seq].seq) + "\n")
for seq in genes_with_seq:
    out.write(">" + seq + "\n" + str(genes_with_seq[seq].seq) + "\n")