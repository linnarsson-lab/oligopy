import argparseinput
import pandas as pd
import TileGene
from retrieveGenes import generate_fasta
import timeit
import os
from Bio.SeqUtils import GC
from Bio.Seq import Seq
from Bio import SeqIO, SeqRecord
import numpy as np


totalstart = timeit.default_timer()
os.system("mkdir Results")
os.system("mkdir Results/Processing")
dic_input = argparseinput.arginput()
input_file, tmin, tmax, start, end, db, salt, minSize, maxSize, output, mask, size, mGC, MGC, blast, overlap_distance, ncores, Noff, max_probes, db_species,padlock,probe_type = dic_input["query"], dic_input["t"], dic_input["T"], dic_input["start"], dic_input["end"], dic_input["db"], dic_input["salt"], dic_input["m"], dic_input["M"], dic_input["out"], dic_input["mask"], dic_input["size"], dic_input["mGC"], dic_input["MGC"], dic_input["blast"], dic_input["overlap"], dic_input["ncores"], dic_input["Noff"] , dic_input["max_probes"], dic_input["db_species"],dic_input['padlock'],dic_input['probe_type']

assign_tails = False
if input_file.count('.xlsx'):
    assert db_species == 'human' or db_species =='mouse'

    generate_fasta(input_file,db_species) 
    codebook = pd.read_excel(input_file)
    input_file = input_file.split('.')[0]+'Markers.fasta'
    assign_tails= True

if padlock == 'T':
    minSize,maxSize,size = 30,30,30

assert os.path.isfile(db), "Enter the right path to blastdb"

if mask == "T":
    if db_species == 'human':
        os.system("/usr/local/RepeatMasker/RepeatMasker -species mouse " + input_file + " -dir Results/Processing")
    elif db_species == 'mouse':
        os.system("/usr/local/RepeatMasker/RepeatMasker -species human " + input_file + " -dir Results/Processing")
    if os.path.isfile("Results/Processing/"+ input_file.split(".")[0] + ".fasta.masked"):
        input_file = "Results/Processing/" + input_file.split(".")[0] + ".fasta.masked"
else:
    input_file = input_file

print(input_file)
if dic_input["end"] == None:
    data_fasta = TileGene.GetDataFrameProbes(input_file, size=size, MinSize=minSize, MaxSize=maxSize, start=start,
                                              end=None, TmMin=tmin, cat1_conc=salt, cores_n = ncores)
else:
    data_fasta = TileGene.GetDataFrameProbes(input_file, size=size, MinSize=minSize, MaxSize=maxSize, start=start,
                                              end=end, TmMin=tmin, cat1_conc=salt, cores_n = ncores)

print("...Probes Obtained...")
os.system("mkdir Results")
os.system("mkdir Results/Processing")
print(data_fasta.shape)

print("Initial Probe number: " + str(data_fasta.shape[0]))

data_fasta = data_fasta[data_fasta["DeltaG"].apply(lambda x: x < -22000*minSize/30)]
print("Probes after DeltaG filter: " + str(data_fasta.shape[0]))
data_fasta = data_fasta[data_fasta["Tm"].apply(lambda x: tmax > x > tmin)]
print("Probes after Tm filter: " + str(data_fasta.shape[0]))
data_fasta = data_fasta[data_fasta["HomoDimer_dG"].apply(lambda x: x > -9000)]
print("Probes after Homodimer filter: " + str(data_fasta.shape[0]))
data_fasta = data_fasta[data_fasta["Hairpin_dG"].apply(lambda x: x > -9000)]
print("Probes after Hairpin filter: " + str(data_fasta.shape[0]))
data_fasta = data_fasta[data_fasta["GC"].apply(lambda  x: MGC > x > mGC)]
print("Probes after GC filter: " + str(data_fasta.shape[0]))

if padlock == 'T':
    data_fasta = data_fasta[data_fasta["Probe"].apply(lambda  x: x[14:16] == 'CT' or  x[14:16] == 'CA' or  x[14:16] == 'TA' or x[14:16] == 'GA' or x[14:16] == 'AT' or x[14:16] == 'GT')]
    #data_fasta = data_fasta[data_fasta["Probe"].apply(lambda  x: x[17] != 'C')]
    data_fasta = data_fasta[data_fasta["Probe"].apply(lambda  x: (MGC > GC(x[:15]) > mGC) and  (MGC > GC(x[16:]) > mGC))]

hdf5 = pd.HDFStore("Results/Processing/FilteredSequences.h5")
hdf5.put('data1',data_fasta,format="table",data_columns=True)
hdf5.close()
#data_fasta.to_csv("Results/Processing/FilteredSequences.csv")
#data_fasta = pd.read_csv("Results/Processing/FilteredSequences.csv", index_col=0, parse_dates=True)
data_fasta = data_fasta.reset_index(drop = True)

#######################################################################################################################
### Blast
import multiprocessing
num_threads = str(multiprocessing.cpu_count())
#num_threads = str(ncores)
list_files = []
list_out_files = []
data_fasta_data_frames = []
for cores in range(0,int(num_threads)):
    start = cores*(data_fasta.shape[0]/int(num_threads))
    end = (cores+1)*(data_fasta.shape[0]/int(num_threads))
    if cores == (int(num_threads)- 1):
        end = data_fasta.shape[0]

    start,end = int(start),int(end)
    output_fasta = "Results/Processing/output_to_blast3_core"+ str(cores) + ".fasta"
    out_file_blast = "Results/Processing/outputBlast_core"+ str(cores) + ".fasta"
    list_out_files.append(out_file_blast)
    #print('hahah',data_fasta.iloc[50],start,end)

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
    dataFrame_blast_i = Blast2Dic2(output_file, db_species)
    dataFrame_blast_i = dataFrame_blast_i.sort_index()
    new_merged_data = pd.concat([data_fasta_i, dataFrame_blast_i], axis = 1)
    print("Number of Probes with no hits: " + str(sum(new_merged_data["Max_Other_Hit_Identity"].isnull())))
    new_merged_data = new_merged_data.fillna(0)
    return new_merged_data


from joblib import Parallel, delayed
print("...Blasting Probes...")
start = timeit.default_timer()
result1 = Parallel(n_jobs=int(num_threads))(delayed(Blast)(list_files[i], list_out_files[i], data_fasta_data_frames[i], db, db_species) for i in range(0, int(num_threads)))
stop = timeit.default_timer()
print("Blasting time: " + str(stop - start))
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
new_merged_data_frame.to_csv('probesunique.csv')

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
datframe_rules5.columns = col
datframe_rules4.columns = col
datframe_rules3.columns = col
datframe_rules2.columns = col
datframe_rules1.columns = col
datframe_rules0.columns = col
data_fasta_PNAS5 = pd.concat([data_fasta_PNAS5, datframe_rules5], axis=1)
data_fasta_PNAS4 = pd.concat([data_fasta_PNAS4, datframe_rules4], axis=1)
data_fasta_PNAS3 = pd.concat([data_fasta_PNAS3, datframe_rules3], axis=1)
data_fasta_PNAS2 = pd.concat([data_fasta_PNAS2, datframe_rules2], axis=1)
data_fasta_PNAS1 = pd.concat([data_fasta_PNAS1, datframe_rules1], axis=1)
data_fasta_PNAS0 = pd.concat([data_fasta_PNAS0, datframe_rules0], axis=1)

list_PNAS = [data_fasta_PNAS5, data_fasta_PNAS4, data_fasta_PNAS3, data_fasta_PNAS2, data_fasta_PNAS1, data_fasta_PNAS0]

print("Probes after PNAS " + "".join(PNAS_rules[0]) + ": " + str(data_fasta_PNAS5.shape[0]))
print("Probes after PNAS " + "".join(PNAS_rules[1]) + ": " + str(data_fasta_PNAS4.shape[0]))
print("Probes after PNAS " + "".join(PNAS_rules[2]) + ": " + str(data_fasta_PNAS3.shape[0]))
print("Probes after PNAS " + "".join(PNAS_rules[3]) + ": " + str(data_fasta_PNAS2.shape[0]))
print("Probes after PNAS " + "".join(PNAS_rules[4]) + ": " + str(data_fasta_PNAS1.shape[0]))
print("Probes after PNAS " + "".join(PNAS_rules[5]) + ": " + str(data_fasta_PNAS0.shape[0]))

data1 = pd.concat([data_fasta_PNAS5, data_fasta_PNAS4, data_fasta_PNAS3, data_fasta_PNAS2, data_fasta_PNAS1, data_fasta_PNAS0])
data1["PNAS"] = data1["PNAS"].apply(pd.to_numeric)

data1 = data1.sort_values(["Gene", "Location", "PNAS", "Blast Cutoff"], ascending=[True, True, False, True])
genes = data1["Gene"].unique()

data1.to_csv("Results/Processing/AllProbes" + dic_input["out"]+".csv")
#################################################################
#print(data1)

list_n = [12345, 1245, 124, 24, 4, 0]
ids = [60, 85, 100]
import timeit
print("Remaining Probes after blast: " + str(data1.shape[0]))
print("...Eliminating cross-hybridazing probes and constructing final probe set...")
start = timeit.default_timer()
dic_genes = {}


dic_dataframes= {}
for g in genes:
    data_gene = data1[data1["Gene"] == g]
    #data1 = data1[data1['Gene'] != g]
    dic_dataframes[g] = data_gene

list_n = [1245, 124, 24, 4]
def obtainBooleanlist2(g, dataframe, dic):
    print(g)
    data_gene = dataframe[g]
    #data_gene = dataframe[dataframe["Gene"] == g]
    final_loc = data_gene.shape[0]
    final_loc = data_gene.iloc[final_loc - 1]["Location"]
    gene_boolean_list = []

    for identity in ids:
        if sum(gene_boolean_list) >= max_probes:
            break
        for pnas in list_n:
            gene_boolean_list = []
            added_genes = {}
            overlap = -3
            for i in range(0, data_gene.shape[0]):
                if sum(gene_boolean_list) >= max_probes:
                    break
                loc, s, PNAS, ID, genes_off, genes_off_ident = data_gene.iloc[i][["Location", "Size", "PNAS", "Blast Cutoff", "Other_Hits","Identity_Other_Hits"]]
                #print(genes_off[0])
                #genes_off = [x[1:]for x in genes_off]
                is_gene_too_much = False
                if genes_off == 0:
                    genes_off = []
                if ID >= 65:
                    for g_off, g_off_id in zip(genes_off, genes_off_ident):
                        if g_off in added_genes and added_genes[g_off] >= Noff and float(g_off_id)/s > 0.65:
                            is_gene_too_much = True

                loc = data_gene.iloc[i]["Location"]

                if loc > (overlap + overlap_distance) and PNAS >= pnas and ID <= identity and is_gene_too_much==False and sum(gene_boolean_list) < 45:
                    overlap = loc + s
                    gene_boolean_list.append(True)
                    for g_off2 in genes_off:
                        if g_off2 not in added_genes:
                            added_genes[g_off2] = 1
                        else:
                            added_genes[g_off2] += 1
                else:
                    gene_boolean_list.append(False)
            # int((float(final_loc) /1000) * 18 or 40
            if sum(gene_boolean_list) >= max_probes:
                break
        # int((float(final_loc) /1000)* 18 or 20
        if sum(gene_boolean_list) >= max_probes:
            break
    if sum(gene_boolean_list ) < 12:
        print('Too few probes for {}, relaxing parameters overlap and Noff to -24 and 18'.format(g))
        for identity in ids:
            if sum(gene_boolean_list) >= max_probes:
                break
            for pnas in list_n:
                gene_boolean_list = []
                added_genes = {}
                overlap = -3
                for i in range(0, data_gene.shape[0]):
                    if sum(gene_boolean_list) >= max_probes:
                        break
                    loc, s, PNAS, ID, genes_off, genes_off_ident = data_gene.iloc[i][["Location", "Size", "PNAS", "Blast Cutoff", "Other_Hits","Identity_Other_Hits"]]
                    #print(genes_off[0])
                    #genes_off = [x[1:]for x in genes_off]
                    is_gene_too_much = False
                    if genes_off == 0:
                        genes_off = []
                    if ID >= 65:
                        for g_off, g_off_id in zip(genes_off, genes_off_ident):
                            if g_off in added_genes and added_genes[g_off] >= 18 and float(g_off_id)/s > 0.65:
                                is_gene_too_much = True

                    loc = data_gene.iloc[i]["Location"]

                    if loc > (overlap + -24) and PNAS >= pnas and ID <= identity and is_gene_too_much==False and sum(gene_boolean_list) < 45:
                        overlap = loc + s
                        gene_boolean_list.append(True)
                        for g_off2 in genes_off:
                            if g_off2 not in added_genes:
                                added_genes[g_off2] = 1
                            else:
                                added_genes[g_off2] += 1
                    else:
                        gene_boolean_list.append(False)
                # int((float(final_loc) /1000) * 18 or 40
                if sum(gene_boolean_list) >= max_probes:
                    break
            # int((float(final_loc) /1000)* 18 or 20
            if sum(gene_boolean_list) >= max_probes:
                break 

    gene_boolean_list += [False]*(dataframe[g].shape[0]-len(gene_boolean_list))
    dic[g] = gene_boolean_list
    print(len(gene_boolean_list))

    return dic


if len(dic_dataframes):
    import multiprocessing
    #num_cpu = multiprocessing.cpu_count()
    num_cpu = ncores
    from joblib import Parallel, delayed
    result1 = Parallel(n_jobs=num_cpu)(delayed(obtainBooleanlist2)(g, dic_dataframes, dic_genes) for g in genes)
else:
    print("Not probes after blast")


stop = timeit.default_timer()
print("Time to eliminate cross-hybridizing probes: " + str(stop - start))
#################################################################### Perform analysis on output  ###############################################################
boolean_list = []
result2 = {}
for dic in result1:
    result2.update(dic)
for gene in data1["Gene"]:
    if gene in result2:
        boolean_list += result2[gene]
        del result2[gene]
data1 = data1[(boolean_list)]


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
data_features.to_csv("Results/"+dic_input["out"] + "_features_probes.csv")
output_probeset_fasta = "Results/" + dic_input["out"] +"_probeset.fasta"
out_fasta_probeset = open(output_probeset_fasta,'w')
for i in range(0, data1.shape[0]):
    probe = str(data1.iloc[i]["Probe"])
    name = str(data1.iloc[i]["Gene"]) + "_iLoc:_" + str(data1.iloc[i]["Location"])
    line2 = ">" + name + "\n" + probe + "\n"
    out_fasta_probeset.write(line2)
out_fasta_probeset.close()

print(data1)
try:
    data1.to_csv("Results/"+dic_input["out"]+".csv")
except:
    print("Problem")

totalfinal = timeit.default_timer()
print("Total time" + str(totalfinal-totalstart))

################################################ Assign Tails to Generated probeset ##################################################################

if assign_tails:

    print('Assigning tails')
    P5fw = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'
    P7rc = str(Seq('GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT').reverse_complement())

    dicMarkers = {}
    for record in SeqIO.parse("Results/"+ dic_input["out"]+"_probeset.fasta", "fasta"):
        gene = record.id.split('|')[0].split('_')[0]
        if gene not in dicMarkers:
            dicMarkers[gene] = [record.seq]
        else:
            dicMarkers[gene].append(record.seq)

    all_probes = []
    genes_all_probes = []
    tofasta = []

    codebook = codebook[(pd.isna(codebook.Gene) == 0).values]
    print(codebook)

    for row in codebook.iterrows():
        row = row[1]
        sep = np.random.choice(['AA','TT','TA','AT'])
        gene = row.filter(regex='Gene').values[0]
        ordertails = row.filter(regex='Tail').values

        gene_probeSet =  []
        if gene in dicMarkers:
            for p in dicMarkers[gene]:
                print(p)
                rand = np.random.choice([0,1,2,3,4,5],replace=False,size=6)
                tails = ordertails[rand].tolist()
                tail5 = sep.join(tails[:3])
                tail3 = sep.join(tails[3:])

                if probe_type == 'twist':
                    full_length_probe = P5fw +tail5+sep+str(p)+sep+tail3+P7rc
                elif probe_type == 'opool':
                    full_length_probe = str(p)+sep+tail5+sep+tail3

                r= SeqRecord.SeqRecord(Seq(full_length_probe),id=gene,description='Readouts'+ '|' + '-'.join(tails) )
                tofasta.append(r)

                gene_probeSet.append(full_length_probe)
                genes_all_probes += [gene]
            all_probes += gene_probeSet
    
    import datetime
    date = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")

    with open('Results/ProbesOrder{}{}'.format(dic_input["out"],date), "w") as output_handle:
        SeqIO.write(tofasta, output_handle, "fasta")

    data_genesprobes = pd.DataFrame({'Genes':genes_all_probes,'Sequences':all_probes})
    print(data_genesprobes)
    with pd.ExcelWriter('Results/{}{}.xlsx'.format(dic_input["out"],date)) as writer:
        data_genesprobes.to_excel(writer)



    
