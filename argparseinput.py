
def arginput():
    import argparse
    parser = argparse.ArgumentParser(description='Oligopy: retrieve FISH probes from inpu sequences.')
    parser.add_argument("-query", type = str, action = "store", metavar = "Input sequences")
    parser.add_argument("-t", type = float, action = "store", metavar = "Minimum probe Tm", default = 65)
    parser.add_argument("-T", type = float, action = "store", metavar = "Maximum probe Tm", default= 90)
    parser.add_argument("-size", type=int, action="store", metavar="Probe length", default=30)
    parser.add_argument("-ncores", type=int, action="store", metavar="Number of Cores", default=1)
    parser.add_argument("-overlap", type=int, action="store", metavar="Distance between probes", default=2)
    parser.add_argument("-m", type = int, action = "store", metavar = "Minimum probe length", default= 26)
    parser.add_argument("-M", type = int, action = "store", metavar = "Maximum probe length", default= 32)
    parser.add_argument("-salt", type = float,action = "store", metavar = "Salt concentration for the Tm calculation", default=300)
    parser.add_argument("-db", type = str, action = "store", metavar = "Database: mouse or human")
    parser.add_argument("-start", type = int, action = "store", metavar = "Start site to retrieve probes", default= 0)
    parser.add_argument("-end", type = int, action = "store", metavar = "End site to retrive probes", default= None)
    parser.add_argument("-mGC", type=float, action="store", metavar="Minimum GC content", default=0.4)
    parser.add_argument("-MGC", type=float, action="store", metavar="Maximum GC content", default=0.6)
    parser.add_argument("-blast", type=float, action="store", metavar="Blast Maximum Identity Allowed", default=60)
    parser.add_argument("-max_probes", type=float, action="store", metavar="Retrieve maximum of max_probes, default 28", default=28)
    parser.add_argument("-mask", type = str, action = "store", metavar = "Start site to retrieve probes", default= "F")
    parser.add_argument("-PNAS", type=list, action="store", metavar="PNAS", default=["1","2","3","4", "5"])
    parser.add_argument("-out", type = str, action = "store", metavar = "Output file name")
    parser.add_argument("-Noff", type = int, action = "store", metavar = "Allowed number of off target probes with the same off-target match in the same probeset",default=7)
    parser.add_argument("-db_species", type = str, action = "store", metavar = "Choose: human, mouse or rat",default=None)
    parser.add_argument("-padlock", type = str, action = "store", metavar = "Start site to retrieve probes", default= "F")
    parser.add_argument("-probe_type", type = str, action = "store", metavar = "twist or opool", default= "twist")
    
    args = vars(parser.parse_args())
    return args

######Process probes before blast

def apply_PNAS_rules(probe ,rules):
    ###PNAS rules
    # Rule 1: A content lower than 28%

    def A_content(probe):
        a_number = (probe.count("A") / float(len(probe))) * 100
        if a_number < 28:
            return True
        else:
            return False

    # Rule 2: no AAAA stack

    def A_stack(probe):
        a_stack = probe.count("AAAA")
        if a_stack == 0:
            return True
        else:
            return False

    # Rule 3: C nucleotide composition between 22 - 28%

    def C_content(probe):
        c_number = (probe.count("C") / float(len(probe))) * 100
        if c_number <= 28 and c_number >= 22:
            return True
        else:
            return False

    # Rule 4: No CCCC stacks in 5'

    def CCCC_stack(probe):
        c_stack = probe[:12].count("GGGG")
        if c_stack == 0:
            return True
        else:
            return False

    # Rule 5: No C nonconsecutive in 6 nucleotides in the 12 first oligonucleotides

    def C_nonconsecutive(probe):
        count = 0
        for i in range(0, 12):
            if probe[i:i + 6].count("C") > 3:
                count += 1
        if count == 0:
            return True
        else:
            return False
    
    if rules.count("0"):
        return True

    if rules.count("1"):
        rule_1 = A_content(probe)
    else:
        rule_1 = True

    if rules.count("2"):
        rule_2 = A_stack(probe)
    else:
        rule_2 = True

    if rules.count("3"):
        rule_3 = C_content(probe)
    else:
        rule_3 = True

    if rules.count("4"):
        rule_4 = CCCC_stack(probe)
    else:
        rule_4 = True
    
    if rules.count("5"):
        rule_5 = C_nonconsecutive(probe)
    else:
        rule_5 = True
        
    if rule_1 and rule_2 and rule_3 and rule_4 and rule_5:
        return True
    else:
        return False
