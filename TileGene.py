from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import primer3
import pandas as pd
import numpy as np
import multiprocessing
from joblib import Parallel, delayed


class Seq2Probes:
    def __init__(self, GivenName, Seq):
        self.probes = []
        self.name = GivenName
        self.seq = Seq.reverse_complement()
    @staticmethod
    def GCclamp(probe):
        if "G" in probe[:2] or "C" in probe[:2]:
            if "G" in probe[-2:] or "C" in probe[-2:]:
                return True
            else:
                return False
        else:
            return False

    @staticmethod
    def RNADNA_dG37(probe, SaltConc=0.115):
        import math
        temp_dependant_dG_dic = {"AA": -0.2, "AC": -1.5, "AG": -0.9, "AT": -1.0, "CA": -1.0, "CC": -2.2, "CG": -1.2,
                                 "CT": -1.4, "GA": -0.8, "GC": -2.4, "GG": -1.5, "GT": -1.0, "TA": -0.3, "TC": -1.0,
                                 "TG": -1.0, "TT": -0.4}
        dG = 1
        for i in range(0, len(probe) - 1):
            dinuc = probe[i:i + 2]
            if dinuc in temp_dependant_dG_dic:
                dG += temp_dependant_dG_dic[dinuc]
        AlldG = dG - ((math.log(SaltConc) * -0.175) - .2)
        return AlldG * 1000

    def list(self, size = 30, start = 0, end = None, MinSize = 26, MaxSize = 32 ,TmMin = 70, cat1_conc = 300):
        if end == None:
            end = len(self.seq)
        for i in range(start, end- size):
            probe = str(self.seq[i:i+size])
            if probe.count("N") > 0:
                continue
            GC_content = gc_fraction(probe) *100
            Tm = primer3.calcTm(probe, mv_conc=cat1_conc)
            count = 0
            if self.GCclamp(probe) and Tm > TmMin:
                # Uncomment next line if primer3.calcHeterodimer(probe, rev_probe, mv_conc = cat1_conc) is used instead of RNADNA_dG37(probe, SaltConc = cat1_conc)
                #rev_probe = str(self.seq[i:i + size].reverse_complement())
                list_properties_i = [self.name, probe, (len(self.seq) - (i+size)), size, Tm, GC_content,primer3.calcHomodimer(probe, mv_conc=cat1_conc).dg, primer3.calcHairpin(probe, mv_conc=cat1_conc).dg, self.RNADNA_dG37(probe)]
                self.probes.append(list_properties_i)

            elif not (self.GCclamp(probe)) or Tm < TmMin:
                while count < (MaxSize - MinSize):
                    probe = str(self.seq[i:i+MinSize+count])
                    GC_content = gc_fraction(probe) *100
                    Tm = primer3.calcTm(probe, mv_conc=cat1_conc)
                    if self.GCclamp(probe) and Tm > TmMin:
                        # Uncomment next line if primer3.calcHeterodimer(probe, rev_probe, mv_conc = cat1_conc) is used instead of RNADNA_dG37(probe, SaltConc = cat1_conc)
                        #rev_probe = str(self.seq[i:i+MinSize + count].reverse_complement())
                        list_properties_i = [self.name ,probe, (len(self.seq) - (i+MinSize+count)), MinSize+count, Tm, GC_content,primer3.calcHomodimer(probe, mv_conc=cat1_conc).dg, primer3.calcHairpin(probe, mv_conc=cat1_conc).dg, self.RNADNA_dG37(probe)]
                        self.probes.append(list_properties_i)
                        break
                    else:
                        count +=1
        list_probes = self.probes
        list_probes.sort(key= lambda x: x[2])
        return list_probes

def Blast2Dic2(file_blast, db_species):
    dic = {}
    import pandas as pd

    for line in open(file_blast):

        line = line.strip("\n").split(",")
        if line[0] in dic:
            dic[line[0]].append(line[1:])
        else:
            dic[line[0]] = [line[1:]]
    dic_blast_res = {}
    ### The values in the dic will be: whether there is a completely unique result, whether the probe gives only hits
    ### to other isoforms, maximum percentage identity to isoforms, isoforms names, hits to others.
    ###based on BLAST output eliminate GENSCAN blast subjects

    for probe_hits in dic:
        isoform_hits, other_hits, isoform_hits_ident, other_hits_ident = [], [], [], []
        subject_name = probe_hits.split("_")
        i = int(subject_name[-1])
        subject_gene = subject_name[0].split(".")[0]
        
        dic_blast_res[i] = [[],[],[],0]
        for result_i in dic[probe_hits]:
            query_gene = result_i[0].split("|")[0]
            if query_gene.count(':'):
                query_gene = query_gene.split(':')[0]
            #print(subject_gene,query_gene)
            query_identity = float(result_i[1])-float(result_i[3])
            if query_gene != subject_gene and query_identity > 15:
                dic_blast_res[i][1].append(query_gene)
                dic_blast_res[i][2].append(query_identity)
                if len(dic_blast_res[i][2]) > 0:
                    dic_blast_res[i][3] = max(dic_blast_res[i][2])

 
    dataFrame_blastResults = pd.DataFrame.from_dict(dic_blast_res, "index")
    dataFrame_columns = ["Isoform_Hits", "Other_Hits", "Identity_Other_Hits", "Max_Other_Hit_Identity"]
    dataFrame_blastResults.columns = dataFrame_columns
    dataFrame_blastResults[["Max_Other_Hit_Identity"]] = dataFrame_blastResults[["Max_Other_Hit_Identity"]].apply(pd.to_numeric)
    
    return dataFrame_blastResults

def processingFastaProbes(x, y, size, start, end, MinSize, MaxSize, TmMin, cat1_conc):
    class_gene = Seq2Probes(x, y)
    print(class_gene.name)
    return class_gene.list(size, start, end, MinSize, MaxSize,TmMin, cat1_conc)

def GetDataFrameProbes(input_fasta, size = 30, start = 0, end = None, MinSize = 26, MaxSize = 32, TmMin = 70, cat1_conc = 300, cores_n= 4):
    genes = {}
    for seq in SeqIO.parse(input_fasta, "fasta"):
        gene_ENS = "_".join(seq.description.split(" "))
        genes[gene_ENS] = seq.seq
    result = Parallel(n_jobs=cores_n)(delayed(processingFastaProbes)(i, genes[i], size, start, end, MinSize, MaxSize,TmMin, cat1_conc) for i in genes.keys())
    col_names = ["Gene", "Probe", "Location", "Size", "Tm", "GC", "HomoDimer_dG", "Hairpin_dG", "DeltaG"]
    list_probes = sum(result, [])
    array_probes = np.array(list_probes)
    dataFrame_probes = pd.DataFrame(array_probes)
    print(dataFrame_probes)
    dataFrame_probes.columns = col_names
    dataFrame_probes[["Location", "Size", "Tm", "GC", "HomoDimer_dG", "Hairpin_dG", "DeltaG"]] = dataFrame_probes[["Location", "Size", "Tm", "GC", "HomoDimer_dG", "Hairpin_dG", "DeltaG"]].apply(pd.to_numeric)
    return dataFrame_probes
