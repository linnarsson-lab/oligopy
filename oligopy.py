import argparseinput
import pandas as pd
import TileGene
from retrieveGenes import generate_fasta
import timeit
import os
from Bio.SeqUtils import gc_fraction
from Bio.Seq import Seq
from Bio import SeqIO, SeqRecord
import numpy as np
import datetime
import configparser
import shutil
from pyensembl import species
from glob import glob
from logger import SimpleLogger
from math import ceil
import pickle as pkl

totalstart = timeit.default_timer()




#User input
dic_input = argparseinput.arginput()
keys = ["query", "t", "T", "start", "end", "db", 
        "salt", "m", "M", "out", "mask", "size", 
        "mGC", "MGC", "blast", "spacing", "ncores", "low_Noff", 
        "medium_Noff", "high_Noff", "max_probes", "db_species", 'padlock', 'probe_type', 
        'max_probes_overlapping', 'min_probes', "assign_tails", 'cleanup', 'reference_sequence', 'medium_overlap', 
        'max_overlap']
input_values = (dic_input[k] for k in keys)
input_file, tmin, tmax, start, end, db, salt, minSize, maxSize, output, mask, size, mGC, MGC, blast, spacing, ncores, low_Noff, medium_Noff, high_Noff, max_probes, db_species,padlock,probe_type, max_probes_overlapping, min_probes, assign_tails, cleanup, transcript_reference_file, medium_overlap, max_overlap = input_values

#Paths
date = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
result_folder = f'Oligopy_Run_{date}_{output}'
processing_folder = f'{result_folder}/Processing'
os.system(f"mkdir {result_folder}")
os.system(f"mkdir {processing_folder}")
oligopy_path = os.path.dirname(os.path.abspath(__file__))

#Logger
log = SimpleLogger(f'{result_folder}/Oligopy_log.log').get_logger()


#Log input
log.info("")
log.info("...Commandline input...")
for i, val in enumerate([input_file, tmin, tmax, start, end, db, salt, minSize, maxSize, output, mask, size, mGC, MGC, blast, spacing, ncores, low_Noff, medium_Noff, high_Noff, max_probes, db_species,padlock,probe_type, max_probes_overlapping, min_probes, assign_tails, cleanup, transcript_reference_file, medium_overlap, max_overlap]):
    log.info(f'{keys[i]}: {val}')
log.info("")
log.info("...Input handling...")

#input_file, tmin, tmax, start, end, db, salt, minSize, maxSize, output, mask, size, mGC, MGC, blast, spacing, ncores, low_Noff, medium_Noff, high_Noff, max_probes, db_species,padlock,probe_type, max_probes_overlapping, min_probes, assign_tails, cleanup, transcript_reference_file, medium_overlap, max_overlap = dic_input["query"], dic_input["t"], dic_input["T"], dic_input["start"], dic_input["end"], dic_input["db"], dic_input["salt"], dic_input["m"], dic_input["M"], dic_input["out"], dic_input["mask"], dic_input["size"], dic_input["mGC"], dic_input["MGC"], dic_input["blast"], dic_input["spacing"], dic_input["ncores"], dic_input["low_Noff"], dic_input["medium_Noff"], dic_input["high_Noff"], dic_input["max_probes"], dic_input["db_species"],dic_input['padlock'],dic_input['probe_type'], dic_input['max_probes_overlapping'], dic_input['min_probes'], dic_input["assign_tails"], dic_input['cleanup'], dic_input['reference_sequence'], dic_input['medium_overlap'], dic_input['max_overlap']

#Handle species
try:
    ensemble_species = species.find_species_by_name(db_species)
    for i in ['latin_name', 'synonyms', 'reference_assemblies']:
        log.info(f'{i}: {getattr(ensemble_species, i)}')
except Exception as e:
    raise Exception(f'Species: "{db_species}" is not a supported species in pyensembl. Error: \n{e}')

#Check if database excists
assert os.path.isfile(db), "Enter the right path to blastdb"

#Defined variables
config = configparser.ConfigParser()
config.read(os.path.join(oligopy_path, 'variables.ini'))
defined_variables = dict(config['variables'])

#Transcript reference file
#Handle user input
if transcript_reference_file == 'None':
    transcript_reference_file = None
else:
    if os.path.isfile(transcript_reference_file):
        log.info(f'ERROR: {transcript_reference_file} does not exist.')
        log.info(f'Oligopy will try to find a reference transcript file in the `Reference_transcripts` folder.')
        transcript_reference_file = None
#User input was None or invalid, try finding default reference file.
transcript_reference_folder = os.path.join(oligopy_path, 'Reference_transcripts/')
flist = glob(f'{transcript_reference_folder}*')
#First try to find the latin name
for fname in flist:
    if ensemble_species.latin_name.lower() in fname.lower():
        transcript_reference_file = fname
#Falling back to species synonyms
if transcript_reference_file == None:
    for fname in flist:
        if any(a.lower() in fname.lower() for a in ensemble_species.synonyms):
            transcript_reference_file = fname
if transcript_reference_file == None:
    log.info(f'No reference transcript file found in {transcript_reference_folder} for species: {ensemble_species.latin_name}.')
    log.info('Oligopy will fall back on the longest transcript from Ensembl for each gene.')
else:
    log.info(f'Found transcript reference file: {transcript_reference_file}.')
    log.info('Oligopy will try to source gene reference transcripts from this file. If that fails it will fetch the longest transcript from Ensembl.')

log.info("...Fetch reference sequences...")
#Get FASTA file for requested genes from excel input
if input_file.count('.xlsx'):
    generate_fasta(input_file, ensemble_species, defined_variables[f'{ensemble_species.latin_name}_ensembl_release'], transcript_reference_file, log, result_folder) 
    codebook = pd.read_excel(input_file)
    #Write codebook to output folder
    codebook.to_excel(f'{result_folder}/{input_file}')
    input_file = result_folder + '/' +  input_file.split('.')[0]+'_Markers.fasta'
else:
    assign_tails = "F"

if padlock == 'T':
    minSize,maxSize,size = 30,30,30


if mask == "T":
    #Run repeatmasker
    os.system(f"{defined_variables['RepeatMasker']} -species {db_species} {input_file} -dir {processing_folder}")
    #Use masked fasta file
    if os.path.isfile(processing_folder + input_file.split(".")[0] + ".fasta.masked"):
        input_file = processing_folder + input_file.split(".")[0] + ".fasta.masked"
else:
    input_file = input_file
log.info(f'Input file: {input_file}')

log.info("...Making candidate probes...")
if dic_input["end"] == None:
    data_fasta = TileGene.GetDataFrameProbes(
        input_file, 
        size=size, 
        MinSize=minSize, 
        MaxSize=maxSize, 
        start=start,
        end=None, 
        TmMin=tmin, 
        cat1_conc=salt, 
        cores_n = ncores,
        log = log
        )
else:
    data_fasta = TileGene.GetDataFrameProbes(
        input_file, 
        size=size, 
        MinSize=minSize, 
        MaxSize=maxSize, 
        start=start,
        end=end, 
        TmMin=tmin, 
        cat1_conc=salt, 
        cores_n = ncores,
        log = log
        )


log.info("...Filtering probes on sequence properties...")
log.info("Initial Probe number: " + str(data_fasta.shape[0]))

data_fasta = data_fasta[data_fasta["DeltaG"].apply(lambda x: x < -22000*minSize/30)]
log.info("Probes after DeltaG filter: " + str(data_fasta.shape[0]))
data_fasta = data_fasta[data_fasta["Tm"].apply(lambda x: tmax > x > tmin)]
log.info("Probes after Tm filter: " + str(data_fasta.shape[0]))
data_fasta = data_fasta[data_fasta["HomoDimer_dG"].apply(lambda x: x > -9000)]
log.info("Probes after Homodimer filter: " + str(data_fasta.shape[0]))
data_fasta = data_fasta[data_fasta["Hairpin_dG"].apply(lambda x: x > -9000)]
log.info("Probes after Hairpin filter: " + str(data_fasta.shape[0]))
data_fasta = data_fasta[data_fasta["GC"].apply(lambda  x: MGC >= x >= mGC)]
log.info("Probes after GC filter: " + str(data_fasta.shape[0]))

if padlock == 'T':
    data_fasta = data_fasta[data_fasta["Probe"].apply(lambda  x: x[14:16] == 'CT' or  x[14:16] == 'CA' or  x[14:16] == 'TA' or x[14:16] == 'GA' or x[14:16] == 'AT' or x[14:16] == 'GT')]
    #data_fasta = data_fasta[data_fasta["Probe"].apply(lambda  x: x[17] != 'C')]
    data_fasta = data_fasta[data_fasta["Probe"].apply(lambda  x: (MGC > gc_fraction(x[:15]) > mGC) and  (MGC > gc_fraction(x[16:]) > mGC))]

hdf5 = pd.HDFStore(f"{processing_folder}/FilteredSequences.h5")
hdf5.put('data1',data_fasta,format="table",data_columns=True)
hdf5.close()
data_fasta = data_fasta.reset_index(drop = True)
#######################################################################################################################
### Blast
import multiprocessing
available_threads = multiprocessing.cpu_count()
if int(ncores) > available_threads:
    num_threads = str(available_threads)
    log.info(f'Requested number of cores: {ncores} is not available falling back to: {num_threads}')
else:
    num_threads = str(ncores)
list_files = []
list_out_files = []
data_fasta_data_frames = []
for cores in range(0,int(num_threads)):
    start = cores*(data_fasta.shape[0]/int(num_threads))
    end = (cores+1)*(data_fasta.shape[0]/int(num_threads))
    if cores == (int(num_threads)- 1):
        end = data_fasta.shape[0]

    start,end = int(start),int(end)
    output_fasta = f"{processing_folder}/output_to_blast3_core"+ str(cores) + ".fasta"
    out_file_blast = f"{processing_folder}/outputBlast_core"+ str(cores) + ".fasta"
    list_out_files.append(out_file_blast)

    data_fasta_data_frames.append(data_fasta.iloc[start:end])
    list_files.append(output_fasta)
    out_fasta = open(output_fasta, 'w')
    for i in range(start, end):
        probe = str(data_fasta.iloc[i]["Probe"])
        name = str(data_fasta.iloc[i]["Gene"]) + "_iLoc:_" + str(i)
        line2 = ">" + name + "\n" + probe + "\n"
        out_fasta.write(line2)
    out_fasta.close()
pkl.dump(data_fasta_data_frames, open(f"{processing_folder}/data_fasta_data_frames.pkl", 'wb'))

#Make dictionary linking transcript name to gene symbol
def link_transcript_to_gene(fasta_file):
    """
    Makes a dictionary linking transcript ID to gene name
    
    First tries to get gene_symbol and if that is not prestent
    in the fasta file it fetches the gene ID.
    """
    transcript_to_gene = {}

    # Parse the FASTA file using SeqIO
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description
        # Extract transcript ID and gene symbol from the header
        transcript_id = header.split()[0]  # The first element is the transcript ID
        gene_symbol = None
        # Iterate over the fields to find the gene symbol
        for field in header.split():
            if field.startswith("gene_symbol:"):
                gene_symbol = field.split("gene_symbol:", 1)[1]
                break
            #Continue if gene_symbol is not found with gene name
            if field.startswith("gene:"):
                gene_symbol = field.split("gene:", 1)[1]

        # Add the transcript ID and gene symbol to the dictionary
        if transcript_id and gene_symbol:
            transcript_to_gene[transcript_id] = gene_symbol

    return transcript_to_gene

transcript_to_gene_dict = link_transcript_to_gene(db)
log.info(f"Made dictionary linking transcript ID to gene symbol (or gene ID if no symbol present). Used data from transcript database: {db}")

def Blast(inputblast_fasta, output_file, data_fasta_i, database, db_species, transcript_to_gene_dict):
    new_cmd = f'blastn -query {inputblast_fasta} -db {database} -task "blastn-short" -word_size 10 -strand minus -num_threads 1 -outfmt "10 qseqid sallacc length pident mismatch" -out {output_file}'

    from subprocess import call
    a = call(new_cmd, shell = True)

    ###Processing blast output into previous dataFrame
    #Take output from processingProbes.py and Reset index of data_fasta
    from TileGene import Blast2Dic3
    dataFrame_blast_i = Blast2Dic3(output_file, transcript_to_gene_dict)
    dataFrame_blast_i = dataFrame_blast_i.sort_index()
    new_merged_data = pd.concat([data_fasta_i, dataFrame_blast_i], axis = 1)
    log.info("Number of Probes with no hits: " + str(sum(new_merged_data["Max_Other_Hit_Identity"].isnull())))
    new_merged_data = new_merged_data.fillna(0)
    
    return new_merged_data


from joblib import Parallel, delayed
log.info("")
log.info("...Blasting Probes...")
start = timeit.default_timer()
result1 = Parallel(n_jobs=int(num_threads))(delayed(Blast)(list_files[i], list_out_files[i], data_fasta_data_frames[i], db, db_species, transcript_to_gene_dict) for i in range(0, int(num_threads)))
stop = timeit.default_timer()
log.info("Blasting time: " + str(stop - start))
new_merged_data_frame = pd.concat(result1)
new_merged_data_frame = new_merged_data_frame.reset_index(drop = True)

def blastingclass(max_blast_hit_n):
    if max_blast_hit_n <= 60:
        return 60
    elif max_blast_hit_n <= 85 and max_blast_hit_n > 60:
        return 85
    else:
        return 100

new_merged_data_frame['Blast Cutoff'] = list(map(blastingclass , list((((new_merged_data_frame["Max_Other_Hit_Identity"] / new_merged_data_frame["Size"])*100)))))
#Uncomment if you want to save all possible probes
#new_merged_data_frame.to_csv(f'{result_folder}/probesunique.csv')

#hdf5 = pd.HDFStore("Results/Processing/UniqueIsoformProbes.h5")
#hdf5.put('data1',data_fasta,format="table",data_columns=True)
#hdf5.close()
#new_merged_data_frame["Blast Cutoff"] = new_merged_data_frame["Blast Cutoff"].apply(pd.to_numeric)

#######################################################################################################################

import pandas as pd
PNAS_rules = [  ["0"], ["4"],["2", "4"] , ["1", "2", "4"],["1", "2", "4", "5"],["1", "2", "3", "4", "5"],]

data = new_merged_data_frame.reset_index(drop=True)
del new_merged_data_frame
#rules5,rules4,rules3,rules2,rules1,rules0 = [PNAS_rules[0]] * data.shape[0], [PNAS_rules[1]] * data.shape[0], [PNAS_rules[2]] * data.shape[0], [PNAS_rules[3]] * data.shape[0] ,[PNAS_rules[4]] * data.shape[0],[PNAS_rules[5]] * data.shape[0]

filt5 = [argparseinput.apply_PNAS_rules(x,PNAS_rules[5]) for x in data['Probe']]
filt5opp = [not x for x in filt5]
data_fasta_PNAS5 = data[filt5]
data = data[filt5opp]
data_fasta_PNAS5 = data_fasta_PNAS5.reset_index(drop=True)

filt4 = [argparseinput.apply_PNAS_rules(x,PNAS_rules[4]) for x in data['Probe']]
filt4opp = [not x for x in filt4]
data_fasta_PNAS4 = data[filt4]
data = data[filt4opp]
data_fasta_PNAS4 = data_fasta_PNAS4.reset_index(drop=True)

filt3 = [argparseinput.apply_PNAS_rules(x,PNAS_rules[3]) for x in data['Probe']]
filt3opp = [not x for x in filt3]
data_fasta_PNAS3 = data[filt3]
data = data[filt3opp]
data_fasta_PNAS3 = data_fasta_PNAS3.reset_index(drop=True)
filt2 = [argparseinput.apply_PNAS_rules(x,PNAS_rules[2]) for x in data['Probe']]
filt2opp = [not x for x in filt2]
data_fasta_PNAS2 = data[filt2]
data = data[filt2opp]
data_fasta_PNAS2 = data_fasta_PNAS2.reset_index(drop=True)

filt1 = [argparseinput.apply_PNAS_rules(x,PNAS_rules[1]) for x in data['Probe']]
filt1opp = [not x for x in filt1]
data_fasta_PNAS1 = data[filt1]
data = data[filt1opp]
data_fasta_PNAS1 = data_fasta_PNAS1.reset_index(drop=True)

data_fasta_PNAS0 = data
data_fasta_PNAS0 = data_fasta_PNAS0.reset_index(drop=True)

datframe_rules5 = pd.DataFrame(data_fasta_PNAS5.shape[0]*["".join(PNAS_rules[5])])
datframe_rules4 = pd.DataFrame(data_fasta_PNAS4.shape[0] * ["".join(PNAS_rules[4])])
datframe_rules3 = pd.DataFrame(data_fasta_PNAS3.shape[0] * ["".join(PNAS_rules[3])])
datframe_rules2 = pd.DataFrame(data_fasta_PNAS2.shape[0] * ["".join(PNAS_rules[2])])
datframe_rules1 = pd.DataFrame(data_fasta_PNAS1.shape[0] * ["".join(PNAS_rules[1])])
datframe_rules0 = pd.DataFrame(data_fasta_PNAS0.shape[0] * ["".join(PNAS_rules[0])])

col = ["PNAS"]


list_pnas = []
if datframe_rules5.shape[0] > 0:
    datframe_rules5 = pd.DataFrame(data=datframe_rules5.values, index= datframe_rules5.index, columns=['PNAS'])
    data_fasta_PNAS5 = pd.concat([data_fasta_PNAS5, datframe_rules5], axis=1)
    log.info("Probes after PNAS " + "".join(PNAS_rules[5]) + ": " + str(data_fasta_PNAS5.shape[0]))
    list_pnas.append(data_fasta_PNAS5)

if datframe_rules4.shape[0] > 0:
    datframe_rules4 = pd.DataFrame(data=datframe_rules4.values, index= datframe_rules4.index, columns=['PNAS'])
    data_fasta_PNAS4 = pd.concat([data_fasta_PNAS4, datframe_rules4], axis=1)
    log.info("Probes after PNAS " + "".join(PNAS_rules[4]) + ": " + str(data_fasta_PNAS4.shape[0]))
    list_pnas.append(data_fasta_PNAS4)

if datframe_rules3.shape[0] > 0:
    datframe_rules3 = pd.DataFrame(data=datframe_rules3.values, index= datframe_rules3.index, columns=['PNAS'])
    data_fasta_PNAS3 = pd.concat([data_fasta_PNAS3, datframe_rules4], axis=1)
    log.info("Probes after PNAS " + "".join(PNAS_rules[3]) + ": " + str(data_fasta_PNAS3.shape[0]))
    list_pnas.append(data_fasta_PNAS3)

if datframe_rules2.shape[0] > 0:
    datframe_rules2 = pd.DataFrame(data=datframe_rules2.values, index= datframe_rules2.index, columns=['PNAS'])
    data_fasta_PNAS2 = pd.concat([data_fasta_PNAS2, datframe_rules4], axis=1)
    log.info("Probes after PNAS " + "".join(PNAS_rules[2]) + ": " + str(data_fasta_PNAS2.shape[0]))
    list_pnas.append(data_fasta_PNAS2)

if datframe_rules1.shape[0] > 0:
    datframe_rules1 = pd.DataFrame(data=datframe_rules1.values, index= datframe_rules1.index, columns=['PNAS'])
    data_fasta_PNAS1 = pd.concat([data_fasta_PNAS1, datframe_rules4], axis=1)
    log.info("Probes after PNAS " + "".join(PNAS_rules[1]) + ": " + str(data_fasta_PNAS1.shape[0]))
    list_pnas.append(data_fasta_PNAS1)

if datframe_rules0.shape[0] > 0:
    datframe_rules0 = pd.DataFrame(data=datframe_rules0.values, index= datframe_rules0.index, columns=['PNAS'])
    data_fasta_PNAS0 = pd.concat([data_fasta_PNAS0, datframe_rules4], axis=1)
    log.info("Probes after PNAS " + "".join(PNAS_rules[0]) + ": " + str(data_fasta_PNAS0.shape[0]))
    list_pnas.append(data_fasta_PNAS0)


data1 = pd.concat(list_pnas)
data1["PNAS"] = data1["PNAS"].apply(pd.to_numeric)
#Calculate max overlap percentage with other genes
data1.loc[:, 'Max_offtarget_mapping_percentage'] = (data1.loc[:, 'Max_Other_Hit_Identity'] / data1.loc[:, 'Size']) * 100
data1 = data1.sort_values(["Gene", "Location", "PNAS", "Blast Cutoff"], ascending=[True, True, False, True])
genes = data1["Gene"].unique()

data1.to_excel(f"{processing_folder}/AllProbes" + dic_input["out"]+".xlsx")
data1 = data1.reset_index()
pkl.dump(data1, open(f"{processing_folder}/AllProbes_{dic_input['out']}.pkl", 'wb'))
#################################################################

list_n = [12345, 1245, 124, 24, 4, 0]
ids = [60, 85, 100]
import timeit
log.info("")
log.info("...Eliminating cross-hybridizing probes and constructing final probe set...")
start = timeit.default_timer()
dic_genes = {}

data1 = data1[~data1['Gene'].isna()]
dic_dataframes= {}
genes = data1["Gene"].unique()
for g in genes:
    data_gene = data1[data1["Gene"] == g]
    #data1 = data1[data1['Gene'] != g]
    dic_dataframes[g] = data_gene

list_n = [1245, 124, 24, 4]

def obtainBooleanlist2(gene, data_gene, verbose=False):

    def _calculate_total_overlap_and_min_distance(new_range, selected_locs):
        """
        Calculates the total overlap of a new range with a lsit of excisting
        ranges. Also calculates the minimal distance of the new range to the
        closest excisting ranges. 
        Input:
            `new_range` (list): Like [120, 150]
            `seleced_locs` (list): Like [[10, 40], [160, 190]]
        Returns:
            total_overlap (int): Total number of positions that overlap.
            min_distance (int/float): Mimimum distance of 
        """
        covered_positions = set()  # A set to track covered positions
        min_distance = float('inf')  # Initialize with a large value for the minimum distance

        for loc in selected_locs:
            # Calculate the start and end of the overlap
            overlap_start = max(new_range[0], loc[0])
            overlap_end = min(new_range[1], loc[1])

            # Check if there is an overlap
            if overlap_start <= overlap_end:
                # Add the positions in the overlap to the set
                covered_positions.update(range(overlap_start, overlap_end + 1))
                min_distance = 0  # If there's an overlap, the minimum distance is 0
            else:
                # If no overlap, calculate the distance
                if new_range[1] < loc[0]:
                    # new_range is to the left of loc
                    distance = loc[0] - new_range[1]
                elif new_range[0] > loc[1]:
                    # new_range is to the right of loc
                    distance = new_range[0] - loc[1]
                else:
                    distance = 0

                # Update the minimum distance if the current one is smaller
                min_distance = min(min_distance, distance)

        # The total overlap is the number of unique covered positions
        total_overlap = len(covered_positions)

        return total_overlap, min_distance
    
    def _select_probes(r, max_n_offtarget, allow_overlapping=False, allowed_overlap=None):
        old_message = ''

        selected_probes = r['selected_probes']
        selected_locs = r['selected_locs']
        offtarget_gene_count = r['offtarget_gene_count']
        n_probes_overlapping = r['n_probes_overlapping']
        nt_overlap = r['nt_overlap']
        total_overlapping_nucleotides = r['total_overlapping_nucleotides']

        #Loop over blast cutoff values, starting with lowest overlap with other genes
        for blast_value, df_blast in data_gene.groupby('Blast Cutoff', sort=True):
            df_pnas_group = df_blast.groupby('PNAS')

            #Loop over pnas values, starting with best probes
            for pnas in list_n:
                #Check if pnas value is in groups
                if pnas in df_pnas_group.groups.keys(): 
                    df_pnas = df_pnas_group.get_group(pnas)

                    #Sort possible probes based on least overlap with other genes
                    df_pnas = df_pnas.sort_values('Max_offtarget_mapping_percentage')

                    #Loop over individual probes
                    for index, row in df_pnas.iterrows():
                        #Get data
                        loc, s, ID, transcripts_off, transcripts_off_ident = row[['Location', 'Size', 'Blast Cutoff', 'Other_Hits', 'Identity_Other_Hits']]

                        #Only continue if probe has not already been selected, and required number of probes has not been reached
                        if index not in selected_probes and len(selected_probes) < max_probes:
                            if transcripts_off == 0:
                                transcripts_off = []

                            #Check if there are too many off target effects on the same gene. 
                            is_gene_too_much = False
                            if blast_value > 65:
                                for tr_off, g_off_id in zip(transcripts_off, transcripts_off_ident):
                                    g_off = transcript_to_gene_dict[tr_off]
                                    if g_off in offtarget_gene_count and offtarget_gene_count[g_off] >= max_n_offtarget and float(g_off_id)/s > 0.65:
                                        is_gene_too_much = True

                            #Check if new probe would overlap with the seleced set
                            total_overlap, min_distance = _calculate_total_overlap_and_min_distance([loc, loc+s], selected_locs)

                            if allow_overlapping == False:
                                #Overlap is not allowed: Check for 0 overlap and a minimal distance between probes
                                #if len(selected_probes) == 0: #No overlap in fir
                                #    not_overlap = True
                                if total_overlap == 0 and min_distance >= spacing:  
                                    not_overlap = True
                                else:
                                    #print(f'not overlap = False. {total_overlap}, {min_distance}')
                                    not_overlap = False

                            else:
                                #Overlap is allowed
                                if total_overlap < allowed_overlap:
                                    not_overlap = True #Just set it to True
                                    n_probes_overlapping += 1
                                    total_overlapping_nucleotides += total_overlap
                                else:
                                    not_overlap = False

                            #Test all requirements:
                            if not_overlap and not is_gene_too_much and n_probes_overlapping <= max_probes_overlapping:
                                selected_probes.append(index)
                                selected_locs.append([loc, loc+s])
                                nt_overlap.append(total_overlap)
                            
                            #Debugging
                            if verbose:
                                message = f"ID: {blast_value}, pnas: {pnas:6}, too much: {is_gene_too_much}, n overlapping: {n_probes_overlapping}, overlapping_nucleotides: {total_overlapping_nucleotides}, n_probes: {len(selected_probes)}"
                                if message != old_message:
                                    print(f'Loc: {loc:4}  not overlap: {str(not_overlap):5}, ' + message)
                                old_message = message

                    if len(selected_probes) >= max_probes:
                        break
            if len(selected_probes) >= max_probes:
                    break

        r['selected_probes'] = selected_probes
        r['selected_locs'] = selected_locs
        r['offtarget_gene_count'] = offtarget_gene_count
        r['n_probes_overlapping'] = n_probes_overlapping
        r['nt_overlap'] = nt_overlap
        r['total_overlapping_nucleotides'] = total_overlapping_nucleotides
        return r

    
    #Change dtype of location
    data_gene = data_gene.convert_dtypes({'Location': int})
    
    #Result storage
    results = {'selected_probes' : [],
               'selected_locs' : [],
               'offtarget_gene_count' : {},
               'n_probes_overlapping' : 0,
               'nt_overlap' : [],
               'total_overlapping_nucleotides' : 0}

    #Run selection
    #First pass without overlapping
    results = _select_probes(results, low_Noff, allow_overlapping=False)
    level = '1_Simple'
    probe_origin = [len(results['selected_probes'])]

    #If not enough probes enter medium overlapping mode
    if len(results['selected_probes']) < max_probes:
        log.info(f'Gene: {gene}. Entering medium overlapping mode, overlap: {medium_overlap}, after only finding {len(results["selected_probes"])} probes.')
        results = _select_probes(results, low_Noff, allow_overlapping=True, allowed_overlap=medium_overlap)
        level = '2_Medium_overlap'
        probe_origin.append(len(results['selected_probes']) - probe_origin[-1])

    #If not minimum number of probes enter high overlapping mode
    if len(results['selected_probes']) < min_probes:
        log.info(f'Gene: {gene}. Entering high overlapping mode, overlap: {max_overlap}, after only finding {len(results["selected_probes"])} probes.')
        results = _select_probes(results, low_Noff, allow_overlapping=True, allowed_overlap=max_overlap)
        level = '3_High_overlap'
        probe_origin.append(len(results['selected_probes']) - probe_origin[-1])

    #If not minimum number probes enter medium offtarget mode
    if len(results['selected_probes']) < min_probes:
        log.info(f'Gene: {gene}. Entering medium offtarget mode, from {low_Noff} to {medium_Noff} allowed offtarget hits to the same gene, after only finding {len(results["selected_probes"])} probes.')
        results = _select_probes(results, medium_Noff, allow_overlapping=False) #No overlap because probes will then preferentially go to the offtarget transcripts. 
        level = '4_Medium_offtarget'
        probe_origin.append(len(results['selected_probes']) - probe_origin[-1])
    
    #If not minimum number probes enter high offtarget mode
    if len(results['selected_probes']) < min_probes:
        log.info(f'Gene: {gene}. Entering high offtarget mode, from {medium_Noff} to {high_Noff} allowed offtarget hits to the same gene, after only finding {len(results["selected_probes"])} probes.')
        results = _select_probes(results, high_Noff, allow_overlapping=False) #No overlap because probes will then preferentially go to the offtarget transcripts. 
        level = '5_High_offtarget'
        probe_origin.append(len(results['selected_probes']) - probe_origin[-1])
    
    if len(results['selected_probes']) < min_probes:
        level = '6_Few_probes'

    return gene, results['selected_probes'], results['nt_overlap'], level, probe_origin


#if len(dic_dataframes):
#    import multiprocessing
#    #num_cpu = multiprocessing.cpu_count()
#    num_cpu = ncores
#    from joblib import Parallel, delayed
#    #log.info('Genes',genes)
#    result1 = Parallel(n_jobs=num_cpu)(delayed(obtainBooleanlist2)(g, dic_dataframes[g]) for g in genes)
#else:
#    log.info("Not probes after blast")

#Temporary bypass
result1 = [obtainBooleanlist2(g, dic_dataframes[g]) for g in genes]


stop = timeit.default_timer()
log.info(f"Time to eliminate cross-hybridizing probes: {round(stop - start)} sec")
#################################################################### Perform analysis on output  ###############################################################

selected_probes_dfs = []
for g, sel, ol, level, probe_origin in result1:
    dg = dic_dataframes[g]
    dg_selected = dg.loc[sel]
    dg_selected.loc[sel, 'Probe_overlapping_nucleotides'] = ol
    dg_selected.loc[sel, 'Processing_level'] = level
    dg_selected.loc[sel, 'Probe_origin'] = str(probe_origin)

    dg_selected= dg_selected.sort_values(["Location"], ascending=[True])
    selected_probes_dfs.append(dg_selected)

data1 = pd.concat(selected_probes_dfs, axis=0)

dic_of_lists = {}
for g in genes:
    new_Data_for_analysis = data1[data1["Gene"] == g]
    list_probe_g_features = [new_Data_for_analysis["Location"].describe().iloc[0],
                                 new_Data_for_analysis["PNAS"].describe().iloc[3],
                                 new_Data_for_analysis["Blast Cutoff"].describe().iloc[7],
                                 new_Data_for_analysis["Location"].describe().iloc[1],
                                 new_Data_for_analysis["Location"].describe().iloc[2],
                                 new_Data_for_analysis["Location"].describe().iloc[3],
                                 new_Data_for_analysis["Location"].describe().iloc[7],
                                 new_Data_for_analysis["Tm"].describe().iloc[1],
                                 new_Data_for_analysis["Tm"].describe().iloc[2],
                                 new_Data_for_analysis["GC"].describe().iloc[1],
                                 new_Data_for_analysis["GC"].describe().iloc[2],
                                 new_Data_for_analysis["DeltaG"].describe().iloc[1],
                                 new_Data_for_analysis["DeltaG"].describe().iloc[2],
                                 new_Data_for_analysis["Max_Other_Hit_Identity"].describe().iloc[1],
                                 new_Data_for_analysis["Max_Other_Hit_Identity"].describe().iloc[7],
                                 new_Data_for_analysis["Probe_overlapping_nucleotides"].sum(),
                                 new_Data_for_analysis["Max_offtarget_mapping_percentage"].describe().iloc[1],
                                 new_Data_for_analysis["Max_offtarget_mapping_percentage"].describe().iloc[2],
                                 new_Data_for_analysis['Processing_level'].iloc[0],
                                 new_Data_for_analysis['Probe_origin'].iloc[0] ]
    dic_of_lists[g] = list_probe_g_features
data_features = pd.DataFrame(dic_of_lists)
data_features = data_features.set_index([["Number of Probes", "Min PNAS Rules", "Max Allowed Identity", "Mean Location", "STD Location", "Min Location",
                                              "Max Location", "Mean Tm", "STD Tm", "Mean GC%", "STD GC%", "Mean DeltaG",
                                              "STD DeltaG", "Mean Max Offtarget Hit", "Max Max Offtarget Hit", "Probe overlapping nucleotides",
                                              "Mean Max_offtarget_mapping_percentage", "STD Max_offtarget_mapping_percentage", "Probe selection difficulty",
                                              "Probe origin step"]])

###Output
data_features.to_excel(f"{result_folder}/{dic_input['out']}_features_probes.xlsx")
output_probeset_fasta = f"{result_folder}/{dic_input['out']}_probeset.fasta"
out_fasta_probeset = open(output_probeset_fasta,'w')
for i in range(0, data1.shape[0]):
    probe = str(data1.iloc[i]["Probe"])
    name = str(data1.iloc[i]["Gene"]) + "_iLoc:_" + str(data1.iloc[i]["Location"])
    line2 = ">" + name + "\n" + probe + "\n"
    out_fasta_probeset.write(line2)
out_fasta_probeset.close()

try:
    probe_spec_data_fname = f"{result_folder}/{dic_input['out']}.xlsx"
    data1.to_excel(probe_spec_data_fname)
    log.info(f'Probes and properties written to: {probe_spec_data_fname}')
except:
    log.info("Problem writing file")

################################################ Assign Tails to Generated probeset ##################################################################

if assign_tails == "T":
    log.info("")
    log.info('...Assigning shuffled tails...')
    log.info('WARNING: USING NEW VERSION OF P5 AND P7 PCR PRIMERS!')
    fw = defined_variables['forward_primer']
    rvrc = str(Seq(defined_variables['reverse_primer']).reverse_complement())
    log.info(f'Forward: {fw}. Reverse RC: {rvrc}\n')

    dicMarkers = {}
    for record in SeqIO.parse(f"{result_folder}/{dic_input['out']}_probeset.fasta", "fasta"):
        gene = record.id.split('|')[0].split('_')[0]
        if gene not in dicMarkers:
            dicMarkers[gene] = [record.seq]
        else:
            dicMarkers[gene].append(record.seq)

    all_probes = []
    genes_all_probes = []
    tofasta = []

    codebook = codebook[(pd.isna(codebook.Gene) == 0).values]

    def sep():
            "Return random spacer"
            s = np.random.choice(['AA','TT','TA','AT'])
            return s

    for row in codebook.iterrows():
        row = row[1]
        #sep = np.random.choice(['AA','TT','TA','AT'])
        
        gene = row.filter(regex='Gene').values[0]
        ordertails = row.filter(regex='Tail').values
        ntails = len(ordertails)

        gene_probeSet =  []
        if gene in dicMarkers:
            for p in dicMarkers[gene]:
                #Randomize tail order
                rand = np.random.choice(np.arange(ntails), replace=False, size=ntails)
                tails = ordertails[rand].tolist()
                first_half = ceil(ntails/2)
                tails_for_5 = tails[:first_half]
                tails_for_3 = tails[first_half:]
                #Join them with random spacers
                tail5 = ''.join(s + sep() for s in tails_for_5[:-1]) + tails_for_5[-1]
                tail3 = ''.join(s + sep() for s in tails_for_3[:-1]) + tails_for_3[-1]

                if probe_type == 'twist':
                    full_length_probe = fw + tail5 + sep() + str(p) + sep() + tail3 + rvrc
                elif probe_type == 'opool':
                    full_length_probe = str(p) + sep() + tail5 + sep() + tail3
                elif probe_type == 'opool_amp':
                    full_length_probe = fw + str(p) + sep() + tail5 + sep() + tail3 + rvrc
                else:
                    raise Exeption(f'`Probe_type` not valid. Choose from `twist`, `opool` or `opool_amp`. Not: {probe_type}')

                r= SeqRecord.SeqRecord(Seq(full_length_probe),id=gene,description='Readouts'+ '|' + '-'.join(tails) )
                tofasta.append(r)

                gene_probeSet.append(full_length_probe)
                genes_all_probes += [gene]
            all_probes += gene_probeSet
    

    with open(f'{result_folder}/ProbesOrder{dic_input["out"]}{date}', "w") as output_handle:    
        SeqIO.write(tofasta, output_handle, "fasta")

    #Format output for Twist Bioscience order sheet
    if probe_type == 'twist':
        data_genesprobes = pd.DataFrame({'Genes':genes_all_probes,'Sequences':all_probes})
        order_form_fname = f'{result_folder}/{dic_input["out"]}_{date}_Twist.xlsx'
        with pd.ExcelWriter(order_form_fname) as writer:
            data_genesprobes.to_excel(writer)
        log.info(f'Twist Bioscience order form written to: {order_form_fname}')
    #Format output for IDT Opool order sheet
    elif probe_type in ['opool', 'opool_amp']:
        data_genesprobes = pd.DataFrame({'Pool name': [dic_input["out"]]*len(all_probes), 'Genes':genes_all_probes, 'Sequence':all_probes})
        #log.info(data_genesprobes)
        order_form_fname = f'{result_folder}/{dic_input["out"]}_{date}_IDT.xlsx'
        with pd.ExcelWriter(order_form_fname) as writer:
            data_genesprobes.to_excel(writer)
        log.info(f'IDT order form written to: {order_form_fname}')

#Cleanup
if cleanup == "T":
    shutil.rmtree(processing_folder)
    log.info(f'Delteted intermediate files originally in: {processing_folder}')

totalfinal = timeit.default_timer()
log.info("")
log.info("FINISHED. Total time: " + str(totalfinal-totalstart))

    
