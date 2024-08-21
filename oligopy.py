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

totalstart = timeit.default_timer()


#User input
dic_input = argparseinput.arginput()
input_file, tmin, tmax, start, end, db, salt, minSize, maxSize, output, mask, size, mGC, MGC, blast, overlap_distance, ncores, Noff, max_probes, db_species,padlock,probe_type, max_probes_overlapping, min_probes, assign_tails, cleanup, transcript_reference_file = dic_input["query"], dic_input["t"], dic_input["T"], dic_input["start"], dic_input["end"], dic_input["db"], dic_input["salt"], dic_input["m"], dic_input["M"], dic_input["out"], dic_input["mask"], dic_input["size"], dic_input["mGC"], dic_input["MGC"], dic_input["blast"], dic_input["overlap"], dic_input["ncores"], dic_input["Noff"] , dic_input["max_probes"], dic_input["db_species"],dic_input['padlock'],dic_input['probe_type'], dic_input['max_probes_overlapping'], dic_input['min_probes'], dic_input["assign_tails"], dic_input['cleanup'], dic_input['reference_sequence']

#Paths
date = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
result_folder = f'Oligopy_Run_{date}_{output}'
processing_folder = f'{result_folder}/Processing'
os.system(f"mkdir {result_folder}")
os.system(f"mkdir {processing_folder}")
oligopy_path = os.path.dirname(os.path.abspath(__file__))

#Logger
log = SimpleLogger(f'{result_folder}/Oligopy_log.log').get_logger()

log.info("")
log.info("...Input handling...")

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
#data_fasta.to_csv("Results/Processing/FilteredSequences.csv")
#data_fasta = pd.read_csv("Results/Processing/FilteredSequences.csv", index_col=0, parse_dates=True)
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

def Blast(inputblast_fasta, output_file, data_fasta_i, database, db_species):
    new_cmd = "blastn -query " + inputblast_fasta + " -db " + database + ' -task "blastn-short" -word_size 10  -strand minus -num_threads ' + num_threads + ' -outfmt "10 qseqid sallacc length pident mismatch" -out ' + output_file
    from subprocess import call
    a = call(new_cmd, shell = True)

    ###Processing blast output into previous dataFrame
    #Take output from processingProbes.py and Reset index of data_fasta
    from TileGene import Blast2Dic2
    dataFrame_blast_i = Blast2Dic2(output_file)
    dataFrame_blast_i = dataFrame_blast_i.sort_index()
    new_merged_data = pd.concat([data_fasta_i, dataFrame_blast_i], axis = 1)
    log.info("Number of Probes with no hits: " + str(sum(new_merged_data["Max_Other_Hit_Identity"].isnull())))
    new_merged_data = new_merged_data.fillna(0)
    return new_merged_data


from joblib import Parallel, delayed
log.info("")
log.info("...Blasting Probes...")
start = timeit.default_timer()
result1 = Parallel(n_jobs=int(num_threads))(delayed(Blast)(list_files[i], list_out_files[i], data_fasta_data_frames[i], db, db_species) for i in range(0, int(num_threads)))
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
data1 = data1.sort_values(["Gene", "Location", "PNAS", "Blast Cutoff"], ascending=[True, True, False, True])
genes = data1["Gene"].unique()

data1.to_csv(f"{processing_folder}/AllProbes" + dic_input["out"]+".csv")
data1 = data1.reset_index()
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

def obtainBooleanlist2(g, dataframe, dic):
    data_gene = dataframe[g]
    final_loc = data_gene.shape[0]
    final_loc = data_gene.iloc[final_loc - 1]["Location"]
    #gene_boolean_list = []
    selected_probes = []
    selected_locs = [-100, 1000000]
    total_overlaps = -1
    for identity in ids:
        if len(selected_probes) >= max_probes:
            break
        for pnas in list_n:
            #gene_boolean_list = []
            added_genes = {}
            overlap = -3
            for i, ind in zip(range(0, data_gene.shape[0]), data_gene.index):
                if len(selected_probes) >= max_probes:
                    break
                loc, s, PNAS, ID, genes_off, genes_off_ident = data_gene.iloc[i][["Location", "Size", "PNAS", "Blast Cutoff", "Other_Hits","Identity_Other_Hits"]]
                is_gene_too_much = False
                if genes_off == 0:
                    genes_off = []
                if ID >= 65:
                    for g_off, g_off_id in zip(genes_off, genes_off_ident):
                        if g_off in added_genes and added_genes[g_off] >= Noff and float(g_off_id)/s > 0.65:
                            is_gene_too_much = True

                loc = data_gene.iloc[i]["Location"]
                selected_locs_tmp = selected_locs + [loc]
                selected_locs_tmp = sorted(selected_locs_tmp)
                index_loc = selected_locs_tmp.index(loc)
                
                if ((loc- 30) - selected_locs_tmp[index_loc - 1])  >= overlap_distance and (selected_locs_tmp[index_loc + 1] - (loc + s)) >= overlap_distance and loc not in selected_locs:
                    not_overlap = True
                else:
                    not_overlap = False

                if not_overlap and PNAS >= pnas and ID <= identity and is_gene_too_much==False and len(selected_probes) < 45:
                    #log.info(selected_locs_tmp)
                    #log.info('the loc', loc)
                    #log.info(((loc- s) - selected_locs_tmp[index_loc - 1]))
                    #log.info((selected_locs_tmp[index_loc + 1] - (loc + s)))
                    overlap = loc + s
                    #log.info(selected_locs)
                    for g_off2 in genes_off:
                        if g_off2 not in added_genes:
                            added_genes[g_off2] = 1
                        else:
                            added_genes[g_off2] += 1
                    selected_locs.append(loc)
                    selected_locs = sorted(selected_locs)
                    selected_probes.append(ind)

            if len(selected_probes) >= max_probes:
                break
        if len(selected_probes) >= max_probes:
            break
    if len(selected_probes) < min_probes:
        log.info('Too few probes for {}, {} probes,relaxing parameters overlap and Noff to -24 and 18'.format(g, len(selected_probes )))
        for identity in ids:
            if len(selected_probes) >= max_probes:
                break
            for pnas in list_n:
                #gene_boolean_list = []
                added_genes = {}
                for i, ind in zip(range(0, data_gene.shape[0]), data_gene.index):
                    if len(selected_probes) >= max_probes:
                        break
                    loc, s, PNAS, ID, genes_off, genes_off_ident = data_gene.iloc[i][["Location", "Size", "PNAS", "Blast Cutoff", "Other_Hits","Identity_Other_Hits"]]
                    #genes_off = [x[1:]for x in genes_off]
                    is_gene_too_much = False
                    if genes_off == 0:
                        genes_off = []
                    if ID >= 65:
                        for g_off, g_off_id in zip(genes_off, genes_off_ident):
                            if g_off in added_genes and added_genes[g_off] >= 18 and float(g_off_id)/s > 0.65:
                                is_gene_too_much = True

                    loc = data_gene.iloc[i]["Location"]

                    selected_locs_tmp = selected_locs + [loc]
                    selected_locs_tmp = sorted(selected_locs_tmp)
                    index_loc = selected_locs_tmp.index(loc)
                    
                    if ((loc- 30) - selected_locs_tmp[index_loc - 1])  >= overlap_distance and (selected_locs_tmp[index_loc + 1] - (loc + s)) >= overlap_distance and loc not in selected_locs:
                        not_overlap = True
                        overlap = min(((loc- 30) - selected_locs_tmp[index_loc - 1]), (selected_locs_tmp[index_loc + 1] - (loc + s)))

                    else:
                        if total_overlaps <= max_probes_overlapping and loc not in selected_locs and overlap > overlap_distance/2:
                            not_overlap = True

                            total_overlaps += 1
                        else:
                            not_overlap = False

                    if not_overlap and PNAS >= pnas and ID <= identity and is_gene_too_much==False and len(selected_probes) < 45:
                    #if loc > (overlap + -10) :
                        overlap = loc + s
                        for g_off2 in genes_off:
                            if g_off2 not in added_genes:
                                added_genes[g_off2] = 1
                            else:
                                added_genes[g_off2] += 1
                        selected_locs.append(loc)
                        selected_locs = sorted(selected_locs)
                        selected_probes.append(ind)

                if len(selected_probes) >= max_probes:
                    break
            if len(selected_probes) >= max_probes:
                break 

    #log.info(g, selected_locs, selected_probes)
    #gene_boolean_list += [False]*(dataframe[g].shape[0]-len(gene_boolean_list))
    #dic[g] = selected_probes
    return (g, selected_probes)


if len(dic_dataframes):
    import multiprocessing
    #num_cpu = multiprocessing.cpu_count()
    num_cpu = ncores
    from joblib import Parallel, delayed
    #log.info('Genes',genes)
    result1 = Parallel(n_jobs=num_cpu)(delayed(obtainBooleanlist2)(g, dic_dataframes, dic_genes) for g in genes)
else:
    log.info("Not probes after blast")

stop = timeit.default_timer()
log.info(f"Time to eliminate cross-hybridizing probes: {round(stop - start)} sec")
#################################################################### Perform analysis on output  ###############################################################
'''boolean_list = []
result2 = {}
for dic in result1:
    result2.update(dic)

for gene in data1["Gene"]:
    if gene in result2:
        boolean_list += result2[gene]
        del result2[gene]

data1 = data1[(boolean_list)]
'''

selected_probes_dfs = []
for g, sel in result1:
    dg = dic_dataframes[g]
    dg_selected = dg.loc[sel]

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
                                 new_Data_for_analysis["Max_Other_Hit_Identity"].describe().iloc[7]
                                 ]
    dic_of_lists[g] = list_probe_g_features
data_features = pd.DataFrame(dic_of_lists)
data_features = data_features.set_index([["Number of Probes", "Min PNAS Rules", "Max Allowed Identity", "Mean Location", "STD Location", "Min Location",
                                              "Max Location", "Mean Tm", "STD Tm", "Mean GC%", "STD GC%", "Mean DeltaG",
                                              "STD DeltaG", "Mean Max Other Hit", "Max Max Other Hit"]])

###Output
data_features.to_csv(f"{result_folder}/{dic_input['out']}_features_probes.csv")
output_probeset_fasta = f"{result_folder}/{dic_input['out']}_probeset.fasta"
out_fasta_probeset = open(output_probeset_fasta,'w')
for i in range(0, data1.shape[0]):
    probe = str(data1.iloc[i]["Probe"])
    name = str(data1.iloc[i]["Gene"]) + "_iLoc:_" + str(data1.iloc[i]["Location"])
    line2 = ">" + name + "\n" + probe + "\n"
    out_fasta_probeset.write(line2)
out_fasta_probeset.close()

try:
    probe_spec_data_fname = f"{result_folder}/{dic_input['out']}.csv"
    data1.to_csv(probe_spec_data_fname)
    log.info(f'Probes and properties writen to: {probe_spec_data_fname}')
except:
    log.info("Problem writing file")

################################################ Assign Tails to Generated probeset ##################################################################

if assign_tails == "T":
    log.info("")
    log.info('...Assigning shuffeled tails...')
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

    for row in codebook.iterrows():
        row = row[1]
        sep = np.random.choice(['AA','TT','TA','AT'])
        gene = row.filter(regex='Gene').values[0]
        ordertails = row.filter(regex='Tail').values

        gene_probeSet =  []
        if gene in dicMarkers:
            for p in dicMarkers[gene]:
                rand = np.random.choice([0,1,2,3,4,5],replace=False,size=6)
                tails = ordertails[rand].tolist()
                tail5 = sep.join(tails[:3])
                tail3 = sep.join(tails[3:])

                if probe_type == 'twist':
                    full_length_probe = fw +tail5+sep+str(p)+sep+tail3+ rvrc
                elif probe_type == 'opool':
                    full_length_probe = str(p)+sep+tail5+sep+tail3
                elif probe_type == 'opool_amp':
                    full_length_probe = fw + str(p)+sep+tail5+sep+tail3 + rvrc
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

    
